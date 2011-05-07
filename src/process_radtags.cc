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
// process_radtags -- clean raw reads using a sliding window approach;
//   split reads by barcode, check RAD cutsite is intact, correct barcodes/cutsites 
///  within one basepair, truncate reads on request.
//
// Julian Catchen
// jcatchen@uoregon.edu
// University of Oregon
//
// $Id: process_radtags.cc 2099 2011-04-30 22:04:37Z catchen $
//

#include "process_radtags.h"

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
string enz;
bool   paired       = false;
bool   clean        = false;
bool   quality      = false;
bool   recover      = false;
bool   interleave   = false;
int    truncate_seq = 0;
int    barcode_size = 0;
double win_size     = 0.15;
int    score_limit  = 10;
int    num_threads  = 1;

map<string, const char **> renz;
map<string, int>           renz_cnt;
map<string, int>           renz_len;

int main (int argc, char* argv[]) {

    parse_command_line(argc, argv);

    renz["sbfI"]  = sbfI;  // CCTGCA/GG, SbfI
    renz["pstI"]  = pstI;  // CTGCA/G, PstI
    renz["notI"]  = notI;  // GC/GGCCGC, NotI
    renz["ecoRI"] = ecoRI; // G/AATTC, EcoRI
    renz["sgrAI"] = sgrAI; // CR/CCGGYG, SgrAI; R=A or G; Y=C or T
    renz_cnt["sbfI"]  = 1;
    renz_cnt["pstI"]  = 1;
    renz_cnt["notI"]  = 1;
    renz_cnt["ecoRI"] = 1;
    renz_cnt["sgrAI"] = 2;
    renz_len["sbfI"]  = 6;
    renz_len["pstI"]  = 5;
    renz_len["notI"]  = 6;
    renz_len["ecoRI"] = 5;
    renz_len["sgrAI"] = 6;

    vector<pair<string, string> > files;
    vector<string> barcodes;
    map<string, ofstream *> pair_1_fhs, pair_2_fhs;
    map<string, map<string, long> > counters, barcode_log;

    build_file_list(files);
    barcode_size = load_barcodes(barcodes);
    open_files(barcodes, pair_1_fhs, pair_2_fhs, counters);

    for (uint i = 0; i < files.size(); i++) {
	cerr << "Processing file " << i+1 << " of " << files.size() << " [" << files[i].first.c_str() << "]\n";

	counters[files[i].first]["total"]       = 0;
	counters[files[i].first]["low_quality"] = 0;
	counters[files[i].first]["noradtag"]    = 0;
	counters[files[i].first]["ambiguous"]   = 0;
	counters[files[i].first]["retained"]    = 0;
	counters[files[i].first]["orphaned"]    = 0;
	counters[files[i].first]["recovered"]   = 0;

	if (paired)
	    process_paired_reads(files[i].first, files[i].second, 
				 pair_1_fhs, pair_2_fhs, 
				 counters[files[i].first], barcode_log);
	else
	    process_reads(files[i].first, 
	    		  pair_1_fhs, 
	    		  counters[files[i].first], barcode_log);

	cerr <<	"  " 
	     << counters[files[i].first]["total"] << " total reads; "
	     << "-" << counters[files[i].first]["ambiguous"]   << " ambiguous barcodes; "
	     << "-" << counters[files[i].first]["noradtag"]    << " ambiguous RAD-Tags; "
	     << "+" << counters[files[i].first]["recovered"]   << " recovered; "
	     << "-" << counters[files[i].first]["low_quality"] << " low quality reads; "
	     << "-" << counters[files[i].first]["orphaned"]    << " orphaned paired-end reads; "
	     << counters[files[i].first]["retained"] << " retained reads.\n";
    }

    cerr << "Closing files, flushing buffers...\n";
    close_file_handles(pair_1_fhs);
    if (paired)
	close_file_handles(pair_2_fhs);

    print_results(barcodes, counters, barcode_log);

    return 0;
}

int process_paired_reads(string prefix_1, 
			 string prefix_2,
			 map<string, ofstream *> &pair_1_fhs, 
			 map<string, ofstream *> &pair_2_fhs, 
			 map<string, long> &counter, 
			 map<string, map<string, long> > &barcode_log) {
    Input *fh_1, *fh_2;
    Read  *r_1, *r_2;

    string path_1 = in_path_1 + "/" + prefix_1;
    string path_2 = in_path_2 + "/" + prefix_2;

    cerr << "Reading data from " << path_1 << " and " << path_2 << "\n";

    if (in_file_type == fastq) {
        fh_1 = new Fastq(path_1.c_str());
	fh_2 = new Fastq(path_2.c_str());
    } else if (in_file_type == bustard) {
        fh_1 = new Bustard(path_1.c_str());
        fh_2 = new Bustard(path_2.c_str());
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

    r_1 = new Read;
    r_1->barcode    = new char[id_len + 1];
    r_1->machine    = new char[id_len + 1];
    r_1->seq        = new char[strlen(s_1->seq) + 1];
    r_1->phred      = new char[strlen(s_1->seq) + 1];
    r_1->int_scores = new  int[strlen(s_1->seq)];
    //
    // Set the parameters for checking read quality later in processing.
    // Window length is 15% (rounded) of the sequence length.
    //
    r_1->len      = strlen(s_1->seq) - barcode_size;
    r_1->win_len  = round(r_1->len * win_size);
    r_1->len     += barcode_size;
    r_1->stop_pos = r_1->len - r_1->win_len;

    r_2 = new Read;
    r_2->barcode    = new char[id_len + 1];
    r_2->machine    = new char[id_len + 1];
    r_2->seq        = new char[strlen(s_2->seq) + 1];
    r_2->phred      = new char[strlen(s_2->seq) + 1];
    r_2->int_scores = new  int[strlen(s_2->seq)];
    r_2->len        = r_1->len;
    r_2->win_len    = r_1->win_len;
    r_2->stop_pos   = r_1->stop_pos;

    long i = 1;

    do {
        if (i % 10000 == 0) cerr << "  Processing RAD-Tag " << i << "       \r";

	parse_input_record(s_1, r_1);
	parse_input_record(s_2, r_2);
	counter["total"] += 2;

	process_singlet(pair_1_fhs, r_1, barcode_log, counter, false);

	if (!r_1->retain) {
	    r_2->retain = 0;
	    counter["orphaned"]++;
	} else {
	    strcpy(r_2->barcode, r_1->barcode);
	    process_singlet(pair_2_fhs, r_2, barcode_log, counter, true);
	}

	if (r_1->retain && r_2->retain) {
	    out_file_type == fastq ? 
		write_fastq(pair_1_fhs, r_1, false) : write_fasta(pair_1_fhs, r_1, false);
            out_file_type == fasta ?
                write_fasta(pair_1_fhs, r_2, true) :
                write_fastq(pair_2_fhs, r_2, true);

	} else if (r_1->retain && !r_2->retain) {
	    // write to a remainder file.
	    out_file_type == fastq ? 
		write_fastq(pair_1_fhs, r_1, false) : write_fasta(pair_1_fhs, r_2, false);
	}

	delete s_1;
	delete s_2;

	i++;
    } while ((s_1 = fh_1->next_seq()) != NULL && 
	     (s_2 = fh_2->next_seq()) != NULL);

    return 0;
}

int process_reads(string prefix, 
		  map<string, ofstream *> &pair_1_fhs, 
		  map<string, long> &counter, 
		  map<string, map<string, long> > &barcode_log) {
    Input *fh;
    Read  *r;

    string path = in_path_1 + "/" + prefix;

    if (in_file_type == fastq)
        fh = new Fastq(path.c_str());
    else if (in_file_type == bustard)
        fh = new Bustard(path.c_str());

    //
    // Read in the first record, initializing the Seq object s. Then 
    // initialize the Read object r, then loop, using the same objects.
    //
    Seq *s = fh->next_seq();
    if (s == NULL) {
    	cerr << "Unable to allocate Seq object.\n";
    	exit(1);
    }

    r = new Read;
    r->barcode    = new char [id_len + 1];
    r->machine    = new char [id_len + 1];
    r->seq        = new char [strlen(s->seq) + 1];
    r->phred      = new char [strlen(s->seq) + 1];
    r->int_scores = new  int [strlen(s->seq)];
    //
    // Set the parameters for checking read quality later in processing.
    // Window length is 15% (rounded) of the sequence length.
    //
    r->len      = strlen(s->seq) - barcode_size;
    r->win_len  = round(r->len * win_size);
    r->len     += barcode_size;
    r->stop_pos = r->len - r->win_len;

    long i = 1;
    do {
	if (i % 10000 == 0) cerr << "  Processing RAD-Tag " << i << "       \r";
	counter["total"]++;

	parse_input_record(s, r);

	process_singlet(pair_1_fhs, r, barcode_log, counter, false);

	 if (r->retain)
	     out_file_type == fastq ? 
	 	write_fastq(pair_1_fhs, r, false) : write_fasta(pair_1_fhs, r, false);

	 delete s;

	i++;
    } while ((s = fh->next_seq()) != NULL);

    //
    // Close the file and delete the Input object.
    //
    delete fh;

    return 0;
}

int process_singlet(map<string, ofstream *> &fhs, Read *href, map<string, map<string, long> > &barcode_log, map<string, long> &counter, bool paired_end) {
    char *p;
    //
    // If requested, truncate this read to $truncate nucleotides
    //
    if (truncate_seq > 0) {
    	href->seq[barcode_size + truncate_seq]   = '\0';
    	href->phred[barcode_size + truncate_seq] = '\0';
    }

    if (paired_end == false) {
    	//
    	// Log the barcodes we receive.
    	//
    	if (barcode_log.count(href->barcode) == 0) {
    	    barcode_log[href->barcode]["noradtag"] = 0;
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

    	//
    	// Is the RADTAG intact?
    	//
	bool rad_cor = false;
        for (int i = 0; i < renz_cnt[enz]; i++) {
	    p = href->seq + barcode_size;

    	    if (strncmp(p, renz[enz][i], renz_len[enz]) == 0)
    		rad_cor = true;
        }
        if (rad_cor == false) {
    	    //
    	    // Try to correct the RAD-Tag.
    	    //
    	    if (!correct_radtag(href, counter)) {
    		barcode_log[href->barcode]["noradtag"]++;
    		counter["noradtag"]++;
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
    if(quality && !check_quality_scores(href, paired_end)) {
    	counter["low_quality"]++;
    	href->retain = 0;
    	return 0;
    } 

    barcode_log[href->barcode]["retained"]++;
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
            barcode_log[href->barcode]["noradtag"] = 0;
        }
	barcode_log[href->barcode]["total"]++;
	return 1;
    }

    return 0;
}

int correct_radtag(Read *href, map<string, long> &counter) {
    if (recover == false)
	return 0;
    //
    // If the RAD-Tag sequence is off by no more than a single nucleotide, correct it.
    //
    int d = 0;

    for (int i = 0; i < renz_cnt[enz]; i++) {
	
        d = dist(renz[enz][i], href->seq+barcode_size);

        if (d <= 1) {
            //
            // Correct the read.
            //
	    strncpy(href->seq + barcode_size, renz[enz][i], renz_len[enz]);
            counter["recovered"]++;

            return 1;
        }
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
    for (int j = 0; j < href->len; j++)
        href->int_scores[j] = href->phred[j] - 64;

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

    int    offset      = paired_end ? 0 : barcode_size - 1;
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
    	    href->retain = 0;
    	    return 0;
    	}

    	p++;
    	q++;
    	j++;
    } while (j <= href->stop_pos);

    return 1;
}

int parse_input_record(Seq *s, Read *r) {
    char *p, *q;
    //
    // Parse FASTQ header that looks like this:
    //  @HWI-ST0747_0141:4:1101:1240:2199#0/1
    //
    char *stop = s->id + strlen(s->seq);

    for (p = s->id, q = p; *q != ':' && q < stop; q++);
    *q = '\0';
    strcpy(r->machine, p);

    for (p = q+1, q = p; *q != ':' && q < stop; q++);
    *q = '\0';
    r->lane = atoi(p);

    for (p = q+1, q = p; *q != ':' && q < stop; q++);
    *q = '\0';
    r->tile = atoi(p);

    for (p = q+1, q = p; *q != ':' && q < stop; q++);
    *q = '\0';
    r->x = atoi(p);

    for (p = q+1, q = p; *q != '#' && q < stop; q++);
    *q = '\0';
    r->y = atoi(p);

    for (p = q+1, q = p; *q != '/' && q < stop; q++);
    *q = '\0';
    r->index = atoi(p);

    for (p = q+1, q = p; *q != '\0' && q < stop; q++);
    r->read = atoi(p);

    strcpy(r->seq, s->seq);
    strcpy(r->phred, s->qual);
    strncpy(r->barcode, r->seq, barcode_size);
    r->retain = 1;
    r->filter = 0;

    return 0;
}

int write_fasta(map<string, ofstream *> &fhs, Read *href, bool paired_end) {
    char tile[id_len];
    sprintf(tile, "%04d", href->tile);

    int offset = paired_end ? 0 : barcode_size;

    *(fhs[href->barcode]) <<
	">" << href->barcode <<
	"_" << href->lane <<
	"_" << tile << 
	"_" << href->x <<
	"_" << href->y <<
	"_" << href->read << "\n" <<
	href->seq + offset << "\n";

    return 0;
}

int write_fastq(map<string, ofstream *> &fhs, Read *href, bool paired_end) {
    //
    // Write the sequence and quality scores in FASTQ format. 
    //
    char tile[id_len];
    sprintf(tile, "%04d", href->tile);

    int offset = paired_end ? 0 : barcode_size;

    if (fhs.count(href->barcode) == 0)
	cerr << "Writing to unknown barcode: '" << href->barcode << "'\n";

    *(fhs[href->barcode]) <<
	"@" << href->barcode <<
	"_" << href->lane << 
	"_" << tile << 
	"_" << href->x << 
	"_" << href->y << 
	"_" << href->read << "\n" <<
	href->seq + offset << "\n" <<
	"+\n" <<
	href->phred + offset << "\n";

    return 0;
}

int print_results(vector<string> &barcodes, map<string, map<string, long> > &counters, map<string, map<string, long> > &barcode_log) {
    map<string, map<string, long> >::iterator it;

    string log_path = out_path + "process_radtags.log";
    ofstream log(log_path.c_str());

    if (log.fail()) {
	cerr << "Unable to open log file '" << log_path << "'\n";
	return 0;
    }

    cerr << "Outputing details to log: '" << log_path << "'\n\n";

    log << "File\t"
	<< "Retained Reads\t"
	<< "Low Quality\t"
	<< "Ambiguous Barcodes\t"
	<< "Ambiguous RAD-Tag\t"
	<< "Orphaned paired-end reads\t"
	<< "Total\n";

    for (it = counters.begin(); it != counters.end(); it++) {
	log << it->first                 << "\t"
	    << it->second["retained"]    << "\t"
	    << it->second["low_quality"] << "\t"
	    << it->second["ambiguous"]   << "\t"
	    << it->second["noradtag"]    << "\t"
	    << it->second["orphaned"]    << "\t"
	    << it->second["total"]       << "\n";
    }

    map<string, long> c;
    c["total"]       = 0;
    c["low_quality"] = 0;
    c["ambiguous"]   = 0;
    c["noradtag"]    = 0;
    c["orphaned"]    = 0;

    //
    // Total up the individual counters
    //
    for (it = counters.begin(); it != counters.end(); it++) {
	c["total"]       += it->second["total"];
	c["low_quality"] += it->second["low_quality"];
	c["ambiguous"]   += it->second["ambiguous"];
	c["orphaned"]    += it->second["orphaned"];
	c["noradtag"]    += it->second["noradtag"];
	c["retained"]    += it->second["retained"]; 
    }

    cerr << 
	c["total"] << " total sequences;\n"
	 << "  " << c["ambiguous"]   << " ambiguous barcode drops;\n"
	 << "  " << c["low_quality"] << " low quality read drops;\n"
	 << "  " << c["noradtag"]    << " ambiguous RAD-Tag drops;\n"
	 << "  " << c["orphaned"]    << " orphaned paired-end reads;\n"
	 << c["retained"] << " retained reads.\n";

    log	<< "Total Sequences\t"      << c["total"]       << "\n"
	<< "Ambiguous Barcodes\t"   << c["ambiguous"]   << "\n"
	<< "Low Quality\t"          << c["low_quality"] << "\n"
	<< "Ambiguous RAD-Tag\t"    << c["noradtag"]    << "\n"
	<< "Orphaned Paired-ends\t" << c["orphaned"]    << "\n"
	<< "Retained Reads\t"       << c["retained"]      << "\n";

    //
    // Print out barcode information.
    //
    log << "\n"
	<< "Barcode\t" 
	<< "Total\t"
	<< "No RadTag\t"
	<< "Retained\n";

    for (uint i = 0; i < barcodes.size(); i++) {
        if (barcode_log.count(barcodes[i]) == 0)
            log << barcodes[i] << "\t" << "0\t" << "0\t" << "0\n";
        else
	    log << barcodes[i] << "\t"
                << barcode_log[barcodes[i]]["total"]    << "\t"
		<< barcode_log[barcodes[i]]["noradtag"] << "\t"
                << barcode_log[barcodes[i]]["retained"] << "\n";
    }

    log << "\n"
	<< "Sequences not recorded\n"
	<< "Barcode\t"
	<< "Total\n";

    //
    // We need to sort by barcode hits.
    //
    for (it = barcode_log.begin(); it != barcode_log.end(); it++) {
	if (counters.count(it->first)) continue;
	if (barcode_log[it->first]["total"] == 0) continue;

	log << it->first << "\t"
	    << it->second["total"] << "\n";
    }

    log.close();

    return 0;
}

int open_files(vector<string> &barcodes, 
	       map<string, ofstream *> &pair_1_fhs, 
	       map<string, ofstream *> &pair_2_fhs, 
	       map<string, map<string, long> > &counters) {
    string path, suffix_1, suffix_2;

    if (out_file_type == fastq) {
	suffix_1 = ".fq";
	suffix_2 = ".fq";
    } else {
	suffix_1 = ".fa";
	suffix_2 = ".fa";
    }
    if (paired && interleave == false) {
	suffix_1 += "_1";
	suffix_2 += "_2";
    }

    for (uint i = 0; i < barcodes.size(); i++) {
	ofstream *fh;

        if (interleave == true) {
	    path = out_path + "sample_" + barcodes[i] + suffix_1;
	    fh = new ofstream(path.c_str());
            pair_1_fhs[barcodes[i]] = fh;
 
	    if (pair_1_fhs[barcodes[i]]->fail()) {
		cerr << "Error opening output file '" << path << "'\n";
		exit(1);
	    }

        } else {
	    path = out_path + "sample_" + barcodes[i] + suffix_1;
	    fh = new ofstream(path.c_str(), ifstream::out);
            pair_1_fhs[barcodes[i]] = fh;

	    if (pair_1_fhs[barcodes[i]]->fail()) {
		cerr << "Error opening output file '" << path << "'\n";
		exit(1);
	    }

            if (paired) {
		path = out_path + "sample_" + barcodes[i] + suffix_2;
		fh = new ofstream(path.c_str(), ifstream::out);
                pair_2_fhs[barcodes[i]] = fh;

		if (pair_2_fhs[barcodes[i]]->fail()) {
		    cerr << "Error opening output file '" << path << "'\n";
		    exit(1);
		}
            }
        }
    }

    return 0;
}

int close_file_handles(map<string, ofstream *> &fhs) {
    map<string, ofstream*>::iterator i;

    for (i = fhs.begin(); i != fhs.end(); i++) {
	i->second->close();
	delete i->second;
    }

    return 0;
}

int load_barcodes(vector<string> &barcodes) {
    char     line[id_len];
    ifstream fh(barcode_file.c_str(), ifstream::in);

    if (fh.fail()) {
        cerr << "Error opening barcode file '" << barcode_file << "'\n";
	exit(1);
    }

    while (fh.good()) {
	fh.getline(line, id_len);

	if (strlen(line) == 0) continue;

	//
	// Check that barcode is legitimate
	//
	for (char *p = line; *p != '\0'; p++)
	    switch (*p) {
	    case 'A':
	    case 'C':
	    case 'G':
	    case 'T':
		break;
	    case 'a':
		*p = 'A';
		break;
	    case 'c':
		*p = 'C';
		break;
	    case 'g':
		*p = 'G';
		break;
	    case 't':
		*p = 'T';
		break;
	    default:
		cerr << "Invalid barcode: " << line << "\n";
		exit(1);
	    }

	barcodes.push_back(string(line));
    }

    fh.close();

    if (barcodes.size() == 0) {
	cerr << "Unable to load any barcodes from '" << barcode_file << "'\n";
	help();
    }

    //
    // Determine the barcode length
    //
    int prev, blen;
    prev = barcodes[0].length();
    for (uint i = 0; i < barcodes.size(); i++) {
        blen = barcodes[i].length();

        if (prev != blen) {
            cerr << "Barcodes must all be the same length. Place different barcode lengths in separate runs.\n";
            help();
        }
    }

    cerr << "Loaded " << barcodes.size() << ", " << blen << "bp barcodes.\n";

    return blen;
}

int build_file_list(vector<pair<string, string> > &files) {

    if (in_path_1.length() > 0) {
	int    pos;
	string file;
	struct dirent *direntry;

	DIR *dir = opendir(in_path_1.c_str());

	if (dir == NULL) {
	    cerr << "Unable to open directory '" << in_path_1 << "' for reading.\n";
	    exit(1);
	}

	while ((direntry = readdir(dir)) != NULL) {
	    file = direntry->d_name;

	    if (file == "." || file == "..")
		continue;

	    pos  = file.find_last_of("/");
	    
	    // Do we want to check the filename to make sure it is legal?

	    files.push_back(make_pair(file.substr(pos+1), ""));
	}

	if (files.size() == 0) {
	    cerr << "Unable to locate any input files to process within '" << in_path_1 << "'\n";
	}
    } else {
	//
	// Break off file path and store path and file name.
	//
	if (paired) {
	    int pos_1   = in_file_p1.find_last_of("/");
	    in_path_1 = in_file_p1.substr(0, pos_1);
	    int pos_2   = in_file_p2.find_last_of("/");
	    in_path_2 = in_file_p2.substr(0, pos_2);
	    files.push_back(make_pair(in_file_p1.substr(pos_1+1), in_file_p2.substr(pos_2+1)));
	} else {
	    int pos   = in_file.find_last_of("/");
	    in_path_1 = in_file.substr(0, pos);
	    files.push_back(make_pair(in_file.substr(pos+1), ""));
	}
    }

    cerr << "Found " << files.size() << " input file(s).\n";

    return 0;
}

int parse_command_line(int argc, char* argv[]) {
    file_type ftype;
    int c;
     
    while (1) {
	static struct option long_options[] = {
	    {"help",         no_argument,       NULL, 'h'},
            {"version",      no_argument,       NULL, 'v'},
            {"quality",      no_argument,       NULL, 'q'},
            {"clean",        no_argument,       NULL, 'c'},
            {"recover",      no_argument,       NULL, 'r'},
	    {"infile_type",  required_argument, NULL, 'i'},
	    {"outfile_type", required_argument, NULL, 'y'},
	    {"file",         required_argument, NULL, 'f'},
	    {"file_p1",      required_argument, NULL, '1'},
	    {"file_p2",      required_argument, NULL, '2'},
	    {"path",         required_argument, NULL, 'p'},
	    {"outpath",      required_argument, NULL, 'o'},
	    {"truncate",     required_argument, NULL, 't'},
	    {"enzyme",       required_argument, NULL, 'e'},
	    {"barcodes",     required_argument, NULL, 'b'},
	    {"window_size",  required_argument, NULL, 'w'},
	    {"score_limit",  required_argument, NULL, 's'},
	    {0, 0, 0, 0}
	};
	
	// getopt_long stores the option index here.
	int option_index = 0;

	c = getopt_long(argc, argv, "hvcqri:y:f:o:t:e:b:1:2:p:s:w:", long_options, &option_index);
     
	// Detect the end of the options.
	if (c == -1)
	    break;
     
	switch (c) {
	case 'h':
	    help();
	    break;
     	case 'i':
            if (strcmp(optarg, "bustard") == 0)
                in_file_type = bustard;
            else
                in_file_type = fastq;
	    break;
     	case 'y':
            if (strcmp(optarg, "fasta") == 0)
                out_file_type = fasta;
	    else 
		out_file_type = fastq;
	    break;
     	case 'f':
	    in_file = optarg;
	    ftype   = fastq;
	    break;
	case 'p':
	    in_path_1 = optarg;
	    ftype     = bustard;
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
	case 'o':
	    out_path = optarg;
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
	case 't':
	    truncate_seq = atoi(optarg);
	    break;
	case 'e':
	    enz = optarg;
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
	cerr << "You must specify either a single input file (-f) or a directory path (-P), not both.\n";
	help();
    }

    if (in_file.length() > 0 && (in_file_p1.length() > 0 || in_file_p2.length() > 0)) {
	cerr << "You must specify either a single input file (-f) or a set of paired files (-1, -2), not both.\n";
	help();
    }

    if (in_path_1.length() > 0 && (in_file_p1.length() > 0 || in_file_p2.length() > 0)) {
	cerr << "You must specify either a file path (-P) or a set of paired files (-1, -2), not both.\n";
	help();
    }

    if (in_path_1.length() > 0 && in_path_1.at(in_path_1.length() - 1) != '/') 
	in_path_1 += "/";

    if (out_path.length() == 0) 
	out_path = ".";

    if (out_path.at(out_path.length() - 1) != '/') 
	out_path += "/";

    if (barcode_file.length() == 0) {
	cerr << "You must specify a file containing barcodes.\n";
	help();
    }

    if (in_file_type == unknown)
	in_file_type = ftype;

    if (enz.length() == 0) {
	cerr << "You must specify the restriction enzyme used.\n";
	help();
    }

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
    std::cerr << "process_radtags " << VERSION << "\n\n";

    exit(0);
}

void help() {
    std::cerr << "process_radtags " << VERSION << "\n"
              << "process_radtags [-f in_file | -p in_dir | -1 pair_1 -2 pair_2] -b barcode_file -o out_dir -e enz [-i type] [-y type] [-c] [-q] [-r] [-w] [-s] [-h]\n"
	      << "  f: path to the input file if processing single-end seqeunces.\n"
	      << "  i: input file type, either 'bustard' for the Illumina BUSTARD output files, or 'fastq' (default 'fastq').\n"
	      << "  p: path to a directory of single-end Bustard files.\n"
	      << "  1: first input file in a set of paired-end sequences.\n"
	      << "  2: second input file in a set of paired-end sequences.\n"
	      << "  o: path to output the processed files.\n"
	      << "  y: output type, either 'fastq' or 'fasta' (default fastq).\n"
	      << "  b: a list of barcodes for this run.\n"
	      << "  e: specify the restriction enzyme to look for (either 'sbfI', 'pstI', 'ecoRI', or 'sgrAI').\n"
	      << "  c: clean data, remove any read with an uncalled base.\n"
	      << "  q: discard reads with low quality scores.\n"
	      << "  r: rescue barcodes and RAD-Tags.\n"
	      << "  t: truncate final read length to this value.\n"
	      << "  w: set the size of the sliding window as a fraction of the read length, between 0 and 1 (default 0.15).\n"
	      << "  s: set the score limit. If the average score within the sliding window drops below this value, the read is discarded (default 10).\n"
	      << "  h: display this help messsage." << "\n\n";

    exit(0);
}
