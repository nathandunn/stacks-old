// -*-mode:c++; c-style:k&r; c-basic-offset:4;-*-
//
// Copyright 2011, Julian Catchen <jcatchen@uoregon.edu>
//
// This file is part of Stacks.
//
// Stacks is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// Stacks is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with Stacks.  If not, see <http://www.gnu.org/licenses/>.
//

//
// process_shortreads -- clean raw reads using a sliding window approach;
//   split reads by barcode if barcodes provided, correct barcodes
//   within one basepair, truncate reads on request.
//
// Julian Catchen
// jcatchen@uoregon.edu
// University of Oregon
//
// $Id: process_shortreads.cc 2099 2011-04-30 22:04:37Z catchen $
//

#include "process_shortreads.h"

//
// Global variables to hold command-line options.
//
file_type in_file_type  = unknown;
file_type out_file_type = fastq;
string in_file;
string in_file_p1;
string in_file_p2;
string in_path_1;
string in_path_2;
string out_path;
string barcode_file;
bool   paired          = false;
bool   clean           = false;
bool   quality         = false;
bool   recover         = false;
bool   interleave      = false;
bool   discards        = false;
bool   overhang        = true;
bool   matepair        = false;
bool   filter_illumina = false;
bool   ill_barcode     = false;
int    truncate_seq = 0;
int    barcode_size = 0;
double win_size     = 0.15;
int    score_limit  = 10;
int    len_limit    = 31;
int    num_threads  = 1;

//
// How to shift FASTQ-encoded quality scores from ASCII down to raw scores
//     score = encoded letter - 64; Illumina version 1.3 - 1.5
//     score = encoded letter - 33; Sanger / Illumina version 1.6+
int qual_offset  = 64;


int main (int argc, char* argv[]) {

    parse_command_line(argc, argv);

    cerr << "Using Phred+" << qual_offset << " encoding for quality scores.\n";
    if (filter_illumina)
	cerr << "Discarding reads marked as 'failed' by Illumina's chastity/purity filters.\n";

    vector<pair<string, string> > files;
    vector<string> barcodes;
    map<string, ofstream *> pair_1_fhs, pair_2_fhs, rem_fhs;
    map<string, map<string, long> > counters, barcode_log;

    build_file_list(files);
    barcode_size = load_barcodes(barcode_file, barcodes);

    open_files(files, barcodes, pair_1_fhs, pair_2_fhs, rem_fhs, counters);

    for (uint i = 0; i < files.size(); i++) {
	cerr << "Processing file " << i+1 << " of " << files.size() << " [" << files[i].first.c_str() << "]\n";

	counters[files[i].first]["total"]        = 0;
	counters[files[i].first]["ill_filtered"] = 0;
	counters[files[i].first]["low_quality"]  = 0;
	counters[files[i].first]["trimmed"]      = 0;
	counters[files[i].first]["ambiguous"]    = 0;
	counters[files[i].first]["retained"]     = 0;
	counters[files[i].first]["orphaned"]     = 0;
	counters[files[i].first]["recovered"]    = 0;

	if (paired)
	    process_paired_reads(files[i].first, files[i].second, 
				 pair_1_fhs, pair_2_fhs, rem_fhs,
				 counters[files[i].first], barcode_log);
	else
	    process_reads(files[i].first, 
	    		  pair_1_fhs,
	    		  counters[files[i].first], barcode_log);

	cerr <<	"  " 
	     << counters[files[i].first]["total"] << " total reads; ";
	if (filter_illumina)
	    cerr << "-" << counters[files[i].first]["ill_filtered"] << " failed Illumina reads; ";
	cerr 
	     << "-" << counters[files[i].first]["ambiguous"]   << " ambiguous barcodes; "
	     << "+" << counters[files[i].first]["recovered"]   << " recovered; "
	     << "+" << counters[files[i].first]["trimmed"]     << " trimmed reads; "
	     << "-" << counters[files[i].first]["low_quality"] << " low quality reads; "
	     << counters[files[i].first]["orphaned"]    << " orphaned paired-ends; "
	     << counters[files[i].first]["retained"] << " retained reads.\n";
    }

    cerr << "Closing files, flushing buffers...\n";
    close_file_handles(pair_1_fhs);
    if (paired) {
	close_file_handles(pair_2_fhs);
	close_file_handles(rem_fhs);
    }

    print_results(argc, argv, barcodes, counters, barcode_log);

    return 0;
}

int process_paired_reads(string prefix_1,
			 string prefix_2,
			 map<string, ofstream *> &pair_1_fhs,
			 map<string, ofstream *> &pair_2_fhs,
			 map<string, ofstream *> &rem_fhs,
			 map<string, long> &counter,
			 map<string, map<string, long> > &barcode_log) {
    Input *fh_1, *fh_2;
    Read  *r_1, *r_2;
    ofstream *discard_fh_1, *discard_fh_2;

    string path_1 = in_path_1 + prefix_1;
    string path_2 = in_path_2 + prefix_2;

    cerr << "  Reading data from\n    " << path_1 << " and\n    " << path_2 << "\n";

    if (in_file_type == fastq) {
        fh_1 = new Fastq(path_1.c_str());
	fh_2 = new Fastq(path_2.c_str());
    } else if (in_file_type == bustard) {
        fh_1 = new Bustard(path_1.c_str());
        fh_2 = new Bustard(path_2.c_str());
    }

    //
    // Open a file for recording discarded reads
    //
    if (discards) {
	path_1 = out_path + prefix_1 + ".discards";
	discard_fh_1 = new ofstream(path_1.c_str(), ifstream::out);

	if (discard_fh_1->fail()) {
	    cerr << "Error opening discard output file '" << path_1 << "'\n";
	    exit(1);
	}

	path_2 = out_path + prefix_2 + ".discards";
	discard_fh_2 = new ofstream(path_2.c_str(), ifstream::out);

	if (discard_fh_1->fail()) {
	    cerr << "Error opening discard output file '" << path_2 << "'\n";
	    exit(1);
	}
    }

    //
    // Read in the first record, initializing the Seq object s. Then 
    // initialize the Read object r, then loop, using the same objects.
    //
    Seq *s_1 = fh_1->next_seq();
    Seq *s_2 = fh_2->next_seq();
    if (s_1 == NULL || s_2 == NULL) {
	cerr << "Unable to allocate Seq object.\n";
	exit(1);
    }

    int buf_len = truncate_seq > 0 ? barcode_size + truncate_seq : strlen(s_1->seq);

    r_1 = new Read;
    r_1->barcode    = new char[id_len  + 1];
    r_1->machine    = new char[id_len  + 1];
    r_1->seq        = new char[buf_len + 1];
    r_1->phred      = new char[buf_len + 1];
    r_1->int_scores = new  int[buf_len];
    r_1->read       = 1;

    //
    // If no barcodes were specified, set r->barcode to be the input file name so
    // that reads are written to an output file of the same name as the input file.
    //
    if (barcode_size == 0)
	strncpy(r_1->barcode, prefix_1.c_str(), id_len);

    //
    // Set the parameters for checking read quality later in processing.
    // Window length is 15% (rounded) of the sequence length.
    //
    r_1->len      = buf_len - barcode_size;
    r_1->win_len  = round(r_1->len * win_size);

    if (r_1->win_len < 1) r_1->win_len = 1;

    r_1->len     += barcode_size;
    r_1->stop_pos = r_1->len - r_1->win_len;

    r_2 = new Read;
    r_2->barcode    = new char[id_len  + 1];
    r_2->machine    = new char[id_len  + 1];
    r_2->seq        = new char[buf_len + 1];
    r_2->phred      = new char[buf_len + 1];
    r_2->int_scores = new  int[buf_len];
    r_2->read       = 2;
    r_2->len        = r_1->len;
    r_2->win_len    = r_1->win_len;
    r_2->stop_pos   = r_1->stop_pos;

    if (barcode_size == 0)
	strncpy(r_2->barcode, prefix_2.c_str(), id_len);

    long i = 1;

    do {
        if (i % 10000 == 0) cerr << "  Processing short read " << i << "       \r";

	parse_input_record(s_1, r_1);
	parse_input_record(s_2, r_2);
	counter["total"] += 2;

	process_singlet(pair_1_fhs, r_1, barcode_log, counter, false);
	process_singlet(pair_2_fhs, r_2, barcode_log, counter, true);

 	if (matepair) {
	    rev_complement(r_1->seq, true, overhang);
	    reverse_qual(r_1->phred, true, overhang);
	}

  	if (r_1->retain && r_2->retain) {
	    //
	    // Check to make sure the barcodes from the pair match. If not,
	    // we assume the wrong molecules have been matched up on the flow cell.
	    //
	    if (barcode_size > 0 && strcmp(r_1->barcode, r_2->barcode) != 0) {
		out_file_type == fastq ? 
		    write_fastq(rem_fhs, r_1, true, overhang) : 
		    write_fasta(rem_fhs, r_1, true, overhang);
		out_file_type == fastq ? 
		    write_fastq(rem_fhs, r_2, true, overhang) : 
		    write_fasta(rem_fhs, r_2, true, overhang);
	    } else {
		out_file_type == fastq ? 
		    write_fastq(pair_1_fhs, r_1, true, overhang) : 
		    write_fasta(pair_1_fhs, r_1, true, overhang);
		out_file_type == fastq ?
		    write_fastq(pair_2_fhs, r_2, true, overhang) :
		    write_fasta(pair_2_fhs, r_2, true, overhang);
	    }

	} else if (r_1->retain && !r_2->retain) {
	    // Write to a remainder file.
	    out_file_type == fastq ? 
		write_fastq(rem_fhs, r_1, true, overhang) : 
		write_fasta(rem_fhs, r_1, true, overhang);

	} else if (!r_1->retain && r_2->retain) {
	    // Write to a remainder file.
	    out_file_type == fastq ? 
		write_fastq(rem_fhs, r_2, true, overhang) : 
		write_fasta(rem_fhs, r_2, true, overhang);
	}

	if (discards && !r_1->retain)
	    out_file_type == fastq ? 
		write_fastq(discard_fh_1, s_1) : 
		write_fasta(discard_fh_1, s_1);
	if (discards && !r_2->retain)
	    out_file_type == fastq ? 
		write_fastq(discard_fh_2, s_2) : 
		write_fasta(discard_fh_2, s_2);

	delete s_1;
	delete s_2;

	i++;
    } while ((s_1 = fh_1->next_seq()) != NULL && 
	     (s_2 = fh_2->next_seq()) != NULL);


    if (discards) {
	delete discard_fh_1;
	delete discard_fh_2;
    }

    delete fh_1;
    delete fh_2;

    return 0;
}

int process_reads(string prefix, 
		  map<string, ofstream *> &pair_1_fhs, 
		  map<string, long> &counter, 
		  map<string, map<string, long> > &barcode_log) {
    Input *fh;
    Read  *r;
    ofstream *discard_fh;

    string path = in_path_1 + prefix;

    if (in_file_type == fastq)
        fh = new Fastq(path.c_str());
    else if (in_file_type == bustard)
        fh = new Bustard(path.c_str());

    //
    // Open a file for recording discarded reads
    //
    if (discards) {
	path = path + ".discards";
	discard_fh = new ofstream(path.c_str(), ifstream::out);

	if (discard_fh->fail()) {
	    cerr << "Error opening discard output file '" << path << "'\n";
	    exit(1);
	}
    }

    //
    // Read in the first record, initializing the Seq object s. Then 
    // initialize the Read object r, then loop, using the same objects.
    //
    Seq *s = fh->next_seq();
    if (s == NULL) {
    	cerr << "Unable to allocate Seq object.\n";
    	exit(1);
    }

    int buf_len = truncate_seq > 0 ? barcode_size + truncate_seq : strlen(s->seq);

    r = new Read;
    r->barcode    = new char [id_len  + 1];
    r->machine    = new char [id_len  + 1];
    r->seq        = new char [buf_len + 1];
    r->phred      = new char [buf_len + 1];
    r->int_scores = new  int [buf_len];
    r->read       = 1;

    //
    // Set the parameters for checking read quality later in processing.
    // Window length is 15% (rounded) of the sequence length.
    //
    r->len      = buf_len - barcode_size;
    r->win_len  = round(r->len * win_size);

    if (r->win_len < 1) r->win_len = 1;

    r->len     += barcode_size;
    r->stop_pos = r->len - r->win_len;

    //
    // If no barcodes were specified, set r->barcode to be the input file name so
    // that reads are written to an output file of the same name as the input file.
    //
    if (barcode_size == 0)
	strncpy(r->barcode, prefix.c_str(), id_len);

    //cerr << "Length: " << r->len << "; Window length: " << r->win_len << "; Stop position: " << r->stop_pos << "\n";

    long i = 1;
    do {
	if (i % 10000 == 0) cerr << "  Processing short read " << i << "       \r";
	counter["total"]++;

	parse_input_record(s, r);

	process_singlet(pair_1_fhs, r, barcode_log, counter, false);

	 if (r->retain)
	     out_file_type == fastq ? 
		 write_fastq(pair_1_fhs, r, true, overhang) : 
		 write_fasta(pair_1_fhs, r, true, overhang);

	 if (discards && !r->retain)
	     out_file_type == fastq ? 
		 write_fastq(discard_fh, s) : 
		 write_fasta(discard_fh, s);

	 delete s;

	i++;
    } while ((s = fh->next_seq()) != NULL);

    if (discards) delete discard_fh;

    //
    // Close the file and delete the Input object.
    //
    delete fh;

    return 0;
}

inline
int process_singlet(map<string, ofstream *> &fhs, Read *href, 
		    map<string, map<string, long> > &barcode_log, map<string, long> &counter, 
		    bool paired_end) {

    if (filter_illumina && href->filter) {
	counter["ill_filtered"]++;
	href->retain = 0;
	return 0;
    }

    if (barcode_size > 0) {
    	//
    	// Log the barcodes we receive.
    	//
    	if (barcode_log.count(href->barcode) == 0) {
    	    barcode_log[href->barcode]["total"]    = 0;
	    barcode_log[href->barcode]["retained"] = 0;
    	}
    	barcode_log[href->barcode]["total"]++;
    
    	//
    	// Is this a legitimate barcode?
    	//
    	if (fhs.count(href->barcode) == 0) {
    	    //
    	    // Try to correct the barcode.
    	    //
    	    if (!correct_barcode(fhs, href, counter, barcode_log)) {
    		counter["ambiguous"]++;
    		href->retain = 0;
    		return 0;
    	    }
    	}
    }

    // Drop this sequence if it has any uncalled nucleotides
    if (clean) {
	for (char *p = href->seq; *p != '\0'; p++)
	    if (*p == '.' || *p == 'N') {
		counter["low_quality"]++;
		href->retain = 0;
		return 0;
	    }
    }

    // Drop this sequence if it has low quality scores
    if(quality) {
	int res = check_quality_scores(href, paired_end);
	if (res == 0) {
	    counter["low_quality"]++;
	    href->retain = 0;
	    return 0;

	} else if (res < 0) {
	    counter["trimmed"]++;
	}
    }

    if (barcode_size > 0) barcode_log[href->barcode]["retained"]++;
    counter["retained"]++;

    return 0;
}

int correct_barcode(map<string, ofstream *> &fhs, Read *href, map<string, long> &counter, map<string, map<string, long> > &barcode_log) {
    if (recover == false)
	return 0;

    //
    // If the barcode sequence is off by no more than a single nucleotide, correct it.
    //
    int d, close;
    string barcode, b, old_barcode;
    map<string, ofstream *>::iterator it;

    close = 0;

    for (it = fhs.begin(); it != fhs.end(); it++) {
	d = dist(it->first.c_str(), href->barcode); 

	if (d <= 1) {
	    close++;
	    b = it->first;
	    break;
	}
    }

    if (close == 1) {
	//
	// Correct the barcode.
	//
	old_barcode = string(href->barcode);
	strcpy(href->barcode, b.c_str());
	for (int i = 0; i < barcode_size; i++)
	    href->seq[i] = href->barcode[i];

	counter["recovered"]++;
	barcode_log[old_barcode]["total"]--;

        if (barcode_log.count(href->barcode) == 0) {
            barcode_log[href->barcode]["total"]    = 0;
            barcode_log[href->barcode]["retained"] = 0;
        }
	barcode_log[href->barcode]["total"]++;
	return 1;
    }

    return 0;
}

int dist(const char *res_enz, char *seq) {
    const char *p; char *q;
    
    int dist = 0;

    for (p = res_enz, q = seq; *p != '\0'; p++, q++)
	if (*p != *q) dist++;

    return dist;
}

int check_quality_scores(Read *href, bool paired_end) {
    //
    // Phred quality scores are discussed here:
    //  http://en.wikipedia.org/wiki/FASTQ_format
    //
    // Illumina 1.3+ encodes phred scores between ASCII values 64 (0 quality) and 104 (40 quality)
    //
    //   @ABCDEFGHIJKLMNOPQRSTUVWXYZ[\]^_`abcdefgh
    //   |         |         |         |         |
    //  64        74        84        94       104
    //   0        10(90%)   20(99%)   30(99.9%) 40(99.99%)
    //
    //
    // Drop sequence if the average phred quality score drops below a threshold within a sliding window.
    //

    //
    // Convert the encoded quality scores to their integer values
    //
    for (uint j = 0; j < href->len; j++)
        href->int_scores[j] = href->phred[j] - qual_offset;

    // for (int j = barcode_size; j <= href->stop_pos; j++) {
    // 	double mean = 0;
    // 	int    stop = j + href->win_len;

    // 	for (int k = j; k < stop; k++)
    // 	    mean += href->int_scores[k];

    // 	mean = mean / href->win_len;

    // 	if (mean < score_limit) {
    // 	    href->retain = 0;
    // 	    return 0;
    // 	}
    // }

    int offset;
    offset = barcode_size > 0 ? barcode_size - 1 : 0;

    double mean        = 0.0;
    double working_sum = 0.0;
    int *p, *q, j;
    //
    // Populate the sliding window.
    //
    for (j = offset; j < href->win_len + offset; j++)
    	working_sum += href->int_scores[j];

    //
    // Set pointers to one position before the first element in the window, and to the last element in the window.
    //
    p = href->int_scores + offset;
    q = p + (int) href->win_len;
    j = offset + 1;
    do {
    	//
    	// Add the score from the front edge of the window, subtract the score
    	// from the back edge of the window.
    	//
    	working_sum -= (double) *p;
    	working_sum += (double) *q;

    	mean = working_sum / href->win_len;

    	if (mean < score_limit) {

	    if (j < len_limit) {
		href->retain = 0;
		return 0;
	    } else {
		href->seq[j]   = '\0';
		href->phred[j] = '\0';
		return -1;
	    }
    	}

    	p++;
    	q++;
    	j++;
    } while (j <= href->stop_pos);

    return 1;
}

int print_results(int argc, char **argv, vector<string> &barcodes, map<string, map<string, long> > &counters, map<string, map<string, long> > &barcode_log) {
    map<string, map<string, long> >::iterator it;

    string log_path = out_path + "process_shortreads.log";
    ofstream log(log_path.c_str());

    if (log.fail()) {
	cerr << "Unable to open log file '" << log_path << "'\n";
	return 0;
    }

    cerr << "Outputing details to log: '" << log_path << "'\n\n";

    //
    // Obtain the current date.
    //
    time_t     rawtime;
    struct tm *timeinfo;
    char       date[32];
    time(&rawtime);
    timeinfo = localtime(&rawtime);
    strftime(date, 32, "%F %T", timeinfo);

    for (int i = 0; i < argc; i++) {
	log << argv[i]; 
	if (i < argc - 1) log << " ";
    }
    log << "\n" << "process_shortreads executed " << date << "\n\n";

    log << "File\t"
	<< "Retained Reads\t";
    if (filter_illumina)
	log << "Illumina Filtered\t";
    log << "Low Quality\t"
	<< "Ambiguous Barcodes\t"
	<< "Trimmed Reads\t"
	<< "Orphaned paired-end reads\t"
	<< "Total\n";

    for (it = counters.begin(); it != counters.end(); it++) {
	log << it->first                 << "\t"
	    << it->second["retained"]    << "\t";
	if (filter_illumina)
	    log << it->second["ill_filtered"] << "\t";
	log << it->second["low_quality"] << "\t"
	    << it->second["ambiguous"]   << "\t"
	    << it->second["trimmed"]    << "\t"
	    << it->second["orphaned"]    << "\t"
	    << it->second["total"]       << "\n";
    }

    map<string, long> c;
    c["total"]       = 0;
    c["low_quality"] = 0;
    c["ill_filtered"] = 0;
    c["ambiguous"]   = 0;
    c["trimmed"]     = 0;
    c["orphaned"]    = 0;

    //
    // Total up the individual counters
    //
    for (it = counters.begin(); it != counters.end(); it++) {
	c["total"]       += it->second["total"];
	c["ill_filtered"] += it->second["ill_filtered"];
	c["low_quality"] += it->second["low_quality"];
	c["ambiguous"]   += it->second["ambiguous"];
	c["trimmed"]     += it->second["trimmed"];
	c["orphaned"]    += it->second["orphaned"];
	c["retained"]    += it->second["retained"];
    }

    cerr << c["total"] << " total sequences;\n";
    if (filter_illumina)
	cerr << "  " << c["ill_filtered"] << " failed Illumina filtered reads;\n";
    cerr << "  " << c["ambiguous"]   << " ambiguous barcode drops;\n"
	 << "  " << c["low_quality"] << " low quality read drops;\n"
	 << "  " << c["trimmed"]     << " trimmed reads;\n"
	 << "  " << c["orphaned"]    << " orphaned paired-end reads;\n"
	 << c["retained"] << " retained reads.\n";

    log	<< "\n" 
	<< "Total Sequences\t"      << c["total"]       << "\n";
    if (filter_illumina)
	log << "Failed Illumina filtered reads\t" << c["ill_filtered"] << "\n";
    log 
	<< "Ambiguous Barcodes\t"   << c["ambiguous"]   << "\n"
	<< "Low Quality\t"          << c["low_quality"] << "\n"
	<< "Trimmed Reads\t"        << c["trimmed"]     << "\n"
	<< "Orphaned Paired-ends\t" << c["orphaned"]    << "\n"
	<< "Retained Reads\t"       << c["retained"]      << "\n";

    if (barcode_size > 0) {
	//
	// Print out barcode information.
	//
	log << "\n"
	    << "Barcode\t" 
	    << "Total\t"
	    << "Retained\n";

	set<string> barcode_list;

	for (uint i = 0; i < barcodes.size(); i++) {
	    barcode_list.insert(barcodes[i]);

	    if (barcode_log.count(barcodes[i]) == 0)
		log << barcodes[i] << "\t" << "0\t" << "0\t" << "0\n";
	    else
		log << barcodes[i] << "\t"
		    << barcode_log[barcodes[i]]["total"]    << "\t"
		    << barcode_log[barcodes[i]]["retained"] << "\n";
	}

	log << "\n"
	    << "Sequences not recorded\n"
	    << "Barcode\t"
	    << "Total\n";

	//
	// Sort unused barcodes by number of occurances.
	//
	vector<pair<string, int> > bcs;
	for (it = barcode_log.begin(); it != barcode_log.end(); it++)
	    bcs.push_back(make_pair(it->first, it->second["total"]));
	sort(bcs.begin(), bcs.end(), compare_barcodes);

	for (uint i = 0; i < bcs.size(); i++) {
	    if (barcode_list.count(bcs[i].first)) continue;
	    if (bcs[i].second == 0) continue;

	    log << bcs[i].first << "\t"
		<< bcs[i].second << "\n";
	}
    }

    log.close();

    return 0;
}

int  compare_barcodes(pair<string, int> a, pair<string, int> b) {
    return a.second > b.second;
}

int parse_command_line(int argc, char* argv[]) {
    file_type ftype;
    int c;
     
    while (1) {
	static struct option long_options[] = {
	    {"help",               no_argument, NULL, 'h'},
            {"version",            no_argument, NULL, 'v'},
            {"quality",            no_argument, NULL, 'q'},
            {"clean",              no_argument, NULL, 'c'},
            {"recover",            no_argument, NULL, 'r'},
	    {"discards",           no_argument, NULL, 'D'},
	    {"paired",             no_argument, NULL, 'P'},
	    {"mate-pair",          no_argument, NULL, 'M'},
	    {"no_overhang",        no_argument, NULL, 'O'},
	    {"filter_illumina",    no_argument, NULL, 'F'},
	    {"illumina_barcodes",  no_argument, NULL, 'I'},
	    {"infile_type",  required_argument, NULL, 'i'},
	    {"outfile_type", required_argument, NULL, 'y'},
	    {"file",         required_argument, NULL, 'f'},
	    {"file_p1",      required_argument, NULL, '1'},
	    {"file_p2",      required_argument, NULL, '2'},
	    {"path",         required_argument, NULL, 'p'},
	    {"outpath",      required_argument, NULL, 'o'},
	    {"truncate",     required_argument, NULL, 't'},
	    {"barcodes",     required_argument, NULL, 'b'},
	    {"window_size",  required_argument, NULL, 'w'},
	    {"score_limit",  required_argument, NULL, 's'},
	    {"encoding",     required_argument, NULL, 'E'},
	    {0, 0, 0, 0}
	};
	
	// getopt_long stores the option index here.
	int option_index = 0;

	c = getopt_long(argc, argv, "hvcqrFIOPDi:y:f:o:t:b:1:2:p:s:w:E:", long_options, &option_index);

	// Detect the end of the options.
	if (c == -1)
	    break;
     
	switch (c) {
	case 'h':
	    help();
	    break;
     	case 'i':
            if (strcasecmp(optarg, "bustard") == 0)
                in_file_type = bustard;
            else
                in_file_type = fastq;
	    break;
     	case 'y':
            if (strcasecmp(optarg, "fasta") == 0)
                out_file_type = fasta;
	    else 
		out_file_type = fastq;
	    break;
     	case 'E':
            if (strcasecmp(optarg, "phred64") == 0)
                qual_offset = 64;
	    else if (strcasecmp(optarg, "phred33") == 0)
		qual_offset = 33;
	    break;
     	case 'f':
	    in_file = optarg;
	    ftype   = fastq;
	    break;
	case 'p':
	    in_path_1 = optarg;
	    in_path_2 = in_path_1;
	    ftype     = fastq;
	    break;
	case '1':
	    paired     = true;
	    in_file_p1 = optarg;
	    ftype      = fastq;
	    break;
	case '2':
	    paired     = true;
	    in_file_p2 = optarg;
	    ftype      = fastq;
	    break;
	case 'P':
	    paired = true;
	    break;
	case 'o':
	    out_path = optarg;
	    break;
	case 'M':
	    matepair = true;
	    break;
	case 'D':
	    discards = true;
	    break;
	case 'q':
	    quality = true;
	    break;
	case 'c':
	    clean = true;
	    break;
	case 'r':
	    recover = true;
	    break;
	case 'O':
	    overhang = false;
	    break;
	case 'F':
	    filter_illumina = true;
	    break;
	case 'I':
	    ill_barcode = true;
	    break;
	case 't':
	    truncate_seq = atoi(optarg);
	    break;
	case 'b':
	    barcode_file = optarg;
	    break;
 	case 'w':
	    win_size = atof(optarg);
	    break;
	case 's':
	    score_limit = atoi(optarg);
	    break;
        case 'v':
            version();
            break;
	case '?':
	    // getopt_long already printed an error message.
	    help();
	    break;
     
	default:
	    cerr << "Unknown command line option '" << (char) c << "'\n";
	    help();
	    abort();
	}
    }

    if (in_file.length() == 0 && in_path_1.length() == 0 && in_file_p1.length() == 0) {
	cerr << "You must specify an input file of a directory path to a set of input files.\n";
	help();
    }

    if (in_file.length() > 0 && in_path_1.length() > 0) {
	cerr << "You must specify either a single input file (-f) or a directory path (-p), not both.\n";
	help();
    }

    if (in_file.length() > 0 && (in_file_p1.length() > 0 || in_file_p2.length() > 0)) {
	cerr << "You must specify either a single input file (-f) or a set of paired files (-1, -2), not both.\n";
	help();
    }

    if (in_path_1.length() > 0 && (in_file_p1.length() > 0 || in_file_p2.length() > 0)) {
	cerr << "You must specify either a file path (-p) or a set of paired files (-1, -2), not both.\n";
	help();
    }

    if (in_path_1.length() > 0 && in_path_1.at(in_path_1.length() - 1) != '/') 
	in_path_1 += "/";

    if (in_path_2.length() > 0 && in_path_2.at(in_path_2.length() - 1) != '/') 
	in_path_2 += "/";

    if (out_path.length() == 0) 
	out_path = ".";

    if (out_path.at(out_path.length() - 1) != '/') 
	out_path += "/";

    if (barcode_file.length() == 0)
	cerr << "No barcodes specified, files will not be demultiplexed.\n";

    if (in_file_type == unknown)
	in_file_type = ftype;

    if (score_limit < 0 || score_limit > 40) {
	cerr << "Score limit must be between 0 and 40.\n";
	help();
    }

    if (win_size < 0 || win_size >= 1) {
	cerr << "Window size is a fraction between 0 and 1.\n";
	help();
    }

    return 0;
}

void version() {
    std::cerr << "process_shortreads " << VERSION << "\n\n";

    exit(0);
}

void help() {
    std::cerr << "process_shortreads " << VERSION << "\n"
              << "process_shortreads [-f in_file | -p in_dir [-P] | -1 pair_1 -2 pair_2] -b barcode_file -o out_dir [-i type] [-y type] [-c] [-q] [-r] [-E encoding] [-t len] [-D] [-w size] [-s lim] [-h]\n"
	      << "  f: path to the input file if processing single-end seqeunces.\n"
	      << "  i: input file type, either 'bustard' for the Illumina BUSTARD output files, or 'fastq' (default 'fastq').\n"
	      << "  p: path to a directory of single-end Illumina files.\n"
	      << "  1: first input file in a set of paired-end sequences.\n"
	      << "  2: second input file in a set of paired-end sequences.\n"
	      << "  P: specify that input is paired (for use with '-p').\n"
	      << "  o: path to output the processed files.\n"
	      << "  y: output type, either 'fastq' or 'fasta' (default fastq).\n"
	      << "  b: a list of barcodes for this run.\n"
	      << "  c: clean data, remove any read with an uncalled base.\n"
	      << "  q: discard reads with low quality scores.\n"
	      << "  r: rescue barcodes.\n"
	      << "  t: truncate final read length to this value.\n"
	      << "  E: specify how quality scores are encoded, 'phred33' (Illumina 1.8+, Sanger) or 'phred64' (Illumina 1.3 - 1.5, default).\n"
	      << "  D: capture discarded reads to a file.\n"
	      << "  w: set the size of the sliding window as a fraction of the read length, between 0 and 1 (default 0.15).\n"
	      << "  s: set the score limit. If the average score within the sliding window drops below this value, the read is discarded (default 10).\n"
	      << "  h: display this help messsage.\n\n"
	      << "  --filter_illumina: discard reads that have been marked by Illumina's chastity/purity filter as failing.\n"
	      << "  --illumina_barcodes: barcodes are not inline, but instead are part of Illumina's FASTQ header.\n"
	      << "  --mate-pair: raw reads are circularized mate-pair data, first read will be reverse complemented.\n"
	      << "  --no_overhang: data does not contain an overhang nucleotide between barcode and seqeunce.\n";

    exit(0);
}
