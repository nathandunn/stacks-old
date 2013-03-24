// -*-mode:c++; c-style:k&r; c-basic-offset:4;-*-
//
// Copyright 2011-2013, Julian Catchen <jcatchen@uoregon.edu>
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
//   within one basepair, truncate reads on request.
//
// Julian Catchen
// jcatchen@uoregon.edu
// University of Oregon
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
string renz_1;
string renz_2;
char  *adapter_1;
char  *adapter_2;
barcodet barcode_type    = null_null;
bool     filter_adapter  = false;
bool     paired          = false;
bool     clean           = false;
bool     quality         = false;
bool     recover         = false;
bool     interleave      = false;
bool     merge           = false;
bool     discards        = false;
bool     overhang        = false;
bool     filter_illumina = false;
bool     check_radtag    = true;
int      truncate_seq = 0;
int      bc_size_1    = 0;
int      bc_size_2    = 0;
int      barcode_dist = 2;
double   win_size     = 0.15;
int      score_limit  = 10;
int      len_limit    = 100000;
int      num_threads  = 1;

//
// How to shift FASTQ-encoded quality scores from ASCII down to raw scores
//     score = encoded letter - 64; Illumina version 1.3 - 1.5
//     score = encoded letter - 33; Sanger / Illumina version 1.6+
int qual_offset  = 33;     

//
// Kmer data for adapter filtering.
//
int kmer_size = 5;
int distance  = 1;
int adp_1_len = 0;
int adp_2_len = 0;
vector<char *> adp_1_keys, adp_2_keys;
AdapterHash adp_1_kmers, adp_2_kmers;

map<string, const char **> renz;
map<string, int>           renz_cnt;
map<string, int>           renz_len;

int main (int argc, char* argv[]) {

    initialize_renz(renz, renz_cnt, renz_len);

    parse_command_line(argc, argv);

    cerr << "Using Phred+" << qual_offset << " encoding for quality scores.\n";
    if (truncate_seq > 0)
	cerr << "Reads will be truncated to " << truncate_seq << "bp\n";
    if (filter_illumina)
	cerr << "Discarding reads marked as 'failed' by Illumina's chastity/purity filters.\n";
    if (filter_adapter) {
	cerr << "Filtering reads for adapter sequence:\n";
	if (paired)
	    cerr << "  " << adapter_1 << "\n"
		 << "  " << adapter_2 << "\n";
	else
	    cerr << "  " << adapter_1 << "\n";
	cerr << "    " << distance << " mismatches allowed to adapter sequence.\n";
	init_adapter_seq(kmer_size, adapter_1, adp_1_len, adp_1_kmers, adp_1_keys);
	if (paired)
	    init_adapter_seq(kmer_size, adapter_2, adp_2_len, adp_2_kmers, adp_2_keys);
    }

    vector<pair<string, string> >        files;
    vector<BarcodePair>                  barcodes;
    set<string>                          se_bc, pe_bc;
    map<BarcodePair, ofstream *>         pair_1_fhs, pair_2_fhs, rem_1_fhs, rem_2_fhs;
    map<string, map<string, long> >      counters;
    map<BarcodePair, map<string, long> > barcode_log;

    build_file_list(files);
    load_barcodes(barcode_file, barcodes, se_bc, pe_bc, bc_size_1, bc_size_2);

    open_files(files, barcodes, pair_1_fhs, pair_2_fhs, rem_1_fhs, rem_2_fhs, counters);

    for (uint i = 0; i < files.size(); i++) {
	cerr << "Processing file " << i+1 << " of " << files.size() << " [" << files[i].first.c_str() << "]\n";

	counters[files[i].first]["total"]        = 0;
	counters[files[i].first]["ill_filtered"] = 0;
	counters[files[i].first]["adapter"]      = 0;
	counters[files[i].first]["low_quality"]  = 0;
	counters[files[i].first]["noradtag"]     = 0;
	counters[files[i].first]["ambiguous"]    = 0;
	counters[files[i].first]["retained"]     = 0;
	counters[files[i].first]["recovered"]    = 0;

	if (paired)
	    process_paired_reads(files[i].first, files[i].second, 
				 se_bc, pe_bc,
				 pair_1_fhs, pair_2_fhs, rem_1_fhs, rem_2_fhs,
				 counters[files[i].first], barcode_log);
	else
	    process_reads(files[i].first, 
			  se_bc, pe_bc,
	    		  pair_1_fhs, 
	    		  counters[files[i].first], barcode_log);

	cerr <<	"  " 
	     << counters[files[i].first]["total"] << " total reads; ";
	if (filter_illumina)
	    cerr << "-" << counters[files[i].first]["ill_filtered"] << " failed Illumina reads; ";
	cerr << "-" << counters[files[i].first]["ambiguous"]   << " ambiguous barcodes; "
	     << "-" << counters[files[i].first]["noradtag"]    << " ambiguous RAD-Tags; "
	     << "+" << counters[files[i].first]["recovered"]   << " recovered; "
	     << "-" << counters[files[i].first]["low_quality"] << " low quality reads; "
	     << counters[files[i].first]["retained"] << " retained reads.\n"
	     << "    ";
	if (filter_adapter)
	    cerr << counters[files[i].first]["adapter"] << " reads with adapter sequence";
	cerr << "\n";
    }

    cerr << "Closing files, flushing buffers...\n";
    close_file_handles(pair_1_fhs);
    if (paired) {
	close_file_handles(rem_1_fhs);
	close_file_handles(rem_2_fhs);
	close_file_handles(pair_2_fhs);
    }

    print_results(argc, argv, barcodes, counters, barcode_log);

    if (filter_adapter) {
	free_adapter_seq(adp_1_keys);
	if (paired)
	    free_adapter_seq(adp_2_keys);
    }

    return 0;
}

int process_paired_reads(string prefix_1, 
			 string prefix_2,
			 set<string> &se_bc, set<string> &pe_bc,
			 map<BarcodePair, ofstream *> &pair_1_fhs, 
			 map<BarcodePair, ofstream *> &pair_2_fhs, 
			 map<BarcodePair, ofstream *> &rem_1_fhs,
			 map<BarcodePair, ofstream *> &rem_2_fhs,
			 map<string, long> &counter, 
			 map<BarcodePair, map<string, long> > &barcode_log) {
    Input *fh_1, *fh_2;
    Read  *r_1, *r_2;
    ofstream *discard_fh_1, *discard_fh_2;

    string path_1 = in_path_1 + prefix_1;
    string path_2 = in_path_2 + prefix_2;

    cerr << "  Reading data from:\n  " << path_1 << " and\n  " << path_2 << "\n";

    if (in_file_type == fastq) {
        fh_1 = new Fastq(path_1.c_str());
	fh_2 = new Fastq(path_2.c_str());
    } else if (in_file_type == gzfastq) {
        fh_1 = new GzFastq(path_1.c_str());
	fh_2 = new GzFastq(path_2.c_str());
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
    	cerr << "Attempting to read first pair of input records, unable to allocate " 
	     << "Seq object (Was the correct input type specified?).\n";
	exit(1);
    }

    // 
    // If there is an inline barcode on either of the reads, we will need to set an 
    // offset to know the start of true sequence for checking quality.
    //
    int se_offset = 0;
    int pe_offset = 0;

    if (barcode_type == inline_null ||
	barcode_type == inline_inline ||
	barcode_type == inline_index)
	se_offset = bc_size_1;
    if (barcode_type == inline_inline || 
	barcode_type == index_inline)
	pe_offset = bc_size_2;

    int buf_len = truncate_seq > 0 ? se_offset + truncate_seq : strlen(s_1->seq);
    r_1 = new Read(buf_len, 1, se_offset, win_size);

    //
    // If no barcodes were specified, set r->barcode to be the input file name so
    // that reads are written to an output file of the same name as the input file.
    //
    if (bc_size_1 == 0)
    	strncpy(r_1->barcode, prefix_1.c_str(), id_len);

    //
    // Compute the parameters for the second read.
    //
    buf_len = truncate_seq > 0 ? pe_offset + truncate_seq : strlen(s_2->seq);
    r_2     = new Read(buf_len, 2, pe_offset, win_size);

    if (bc_size_1 == 0)
	strncpy(r_2->barcode, prefix_2.c_str(), id_len);
	
    BarcodePair bc;
    long i = 1;

    do {
        if (i % 10000 == 0) cerr << "  Processing RAD-Tag " << i << "       \r";

	parse_input_record(s_1, r_1);
	parse_input_record(s_2, r_2);
	counter["total"] += 2;

	bc.set(r_1->barcode, r_2->barcode);

	process_barcode(r_1, r_2, bc, pair_1_fhs, se_bc, pe_bc, barcode_log, counter);

	if (r_1->retain) 
	    process_singlet(r_1, renz_1, se_offset, false, barcode_log[bc], counter);
	if (r_2->retain) 
	    process_singlet(r_2, renz_2, pe_offset, true,  barcode_log[bc], counter);

	if (r_1->retain && r_2->retain) {
	    out_file_type == fastq ?
		write_fastq(pair_1_fhs[bc], r_1, overhang) :
 		write_fasta(pair_1_fhs[bc], r_1, overhang);
            out_file_type == fastq ?
                write_fastq(pair_2_fhs[bc], r_2, overhang) :
                write_fasta(pair_2_fhs[bc], r_2, overhang);

	} else if (r_1->retain && !r_2->retain) {
	    //
	    // Write to the remainder file.
	    //
	    out_file_type == fastq ? 
		write_fastq(rem_1_fhs[bc], r_1, overhang) : 
		write_fasta(rem_1_fhs[bc], r_1, overhang);
	} else if (!r_1->retain && r_2->retain) {
	    //
	    // Write to the remainder file.
	    //
	    out_file_type == fastq ? 
		write_fastq(rem_2_fhs[bc], r_2, overhang) : 
		write_fasta(rem_2_fhs[bc], r_2, overhang);
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
		  set<string> &se_bc, set<string> &pe_bc,
		  map<BarcodePair, ofstream *> &pair_1_fhs, 
		  map<string, long> &counter, 
		  map<BarcodePair, map<string, long> > &barcode_log) {
    Input *fh;
    Read  *r;
    ofstream *discard_fh;

    string path = in_path_1 + prefix;

    if (in_file_type == fastq)
        fh = new Fastq(path.c_str());
    else if (in_file_type == gzfastq)
        fh = new GzFastq(path.c_str());
    else if (in_file_type == bustard)
        fh = new Bustard(path.c_str());

    //
    // Open a file for recording discarded reads
    //
    if (discards) {
	path = out_path + prefix + ".discards";
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
    	cerr << "Attempting to read first input record, unable to allocate " 
	     << "Seq object (Was the correct input type specified?).\n";
    	exit(1);
    }

    int buf_len = truncate_seq > 0 ? bc_size_1 + truncate_seq : strlen(s->seq);

    r = new Read(buf_len, 1, bc_size_1, win_size);

    int se_offset = 0;

    if (barcode_type == inline_null)
	se_offset = bc_size_1;

    //
    // If no barcodes were specified, set r->barcode to be the input file name so
    // that reads are written to an output file of the same name as the input file.
    //
    if (bc_size_1 == 0)
	strncpy(r->barcode, prefix.c_str(), id_len);

    //cerr << "Length: " << r->len << "; Window length: " << r->win_len << "; Stop position: " << r->stop_pos << "\n";

    BarcodePair bc;
    long i = 1;
    do {
	if (i % 10000 == 0) cerr << "  Processing RAD-Tag " << i << "       \r";
	counter["total"]++;

	parse_input_record(s, r);

	bc.set(r->barcode);

	process_barcode(r, NULL, bc, pair_1_fhs, se_bc, pe_bc, barcode_log, counter);

	if (r->retain) 
	    process_singlet(r, renz_1, se_offset, false, barcode_log[bc], counter);

	 if (r->retain)
	     out_file_type == fastq ? 
		 write_fastq(pair_1_fhs[bc], r, overhang) : 
		 write_fasta(pair_1_fhs[bc], r, overhang);

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
int 
process_singlet(Read *href, 
		string res_enz, int offset, bool paired_end,
		map<string, long> &bc_log, map<string, long> &counter) 
{
    char *p;

    if (filter_illumina && href->filter) {
	counter["ill_filtered"]++;
	href->retain = 0;
	return 0;
    }

    //
    // Is the RADTAG intact?
    //
    if (check_radtag && res_enz.length() > 0) {
	bool rad_cor = false;

	for (int i = 0; i < renz_cnt[res_enz]; i++) {
	    p = href->seq + offset;

	    if (strncmp(p, renz[res_enz][i], renz_len[res_enz]) == 0)
		rad_cor = true;
	}
	if (rad_cor == false) {
	    //
	    // Try to correct the RAD-Tag.
	    //
	    if (!correct_radtag(href, offset, res_enz, counter)) {
	        if (barcode_type != null_null) bc_log["noradtag"]++;
	        counter["noradtag"]++;
	        href->retain = 0;
	        return 0;
	    }
	}
    }

    //
    // Drop this sequence if it has any uncalled nucleotides.
    //
    if (clean) {
	for (char *p = href->seq + offset; *p != '\0'; p++)
	    if (*p == '.' || *p == 'N') {
		counter["low_quality"]++;
		href->retain = 0;
		return 0;
	    }
    }

    //
    // Drop this sequence if it has low quality scores.
    //
    if (quality && 
	check_quality_scores(href, qual_offset, score_limit, len_limit, offset) <= 0) {
    	counter["low_quality"]++;
    	href->retain = 0;
    	return 0;
    }

    //
    // Drop this sequence if it contains adapter sequence.
    //
    if (filter_adapter) {
	int res;
	if (paired_end)
	    res = filter_adapter_seq(href, adapter_2, adp_2_len, adp_2_kmers,
				     kmer_size, distance, len_limit);
	else
	    res = filter_adapter_seq(href, adapter_1, adp_1_len, adp_1_kmers,
				     kmer_size, distance, len_limit);
	if (res <= 0) {
	    counter["adapter"]++;
	    href->retain = 0;
	    return 0;
	}
    }

    if (barcode_type != null_null) 
	bc_log["retained"]++;
    counter["retained"]++;

    return 0;
}

int 
correct_radtag(Read *href, int offset, string res_enz, map<string, long> &counter) 
{
    if (recover == false)
	return 0;
    //
    // If the RAD-Tag sequence is off by no more than a single nucleotide, correct it.
    //
    int d = 0;

    for (int i = 0; i < renz_cnt[res_enz]; i++) {
	
        d = dist(renz[res_enz][i], href->seq + offset);

        if (d <= 1) {
            //
            // Correct the read.
            //
	    strncpy(href->seq + offset, renz[res_enz][i], renz_len[res_enz]);
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

int 
print_results(int argc, char **argv, 
	      vector<BarcodePair> &barcodes, 
	      map<string, map<string, long> > &counters, 
	      map<BarcodePair, map<string, long> > &barcode_log) 
{
    map<string, map<string, long> >::iterator it;

    string log_path = out_path + "process_radtags.log";
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
    log << "\n" << "process_radtags executed " << date << "\n\n";

    log << "File\t"
	<< "Retained Reads\t";
    if (filter_illumina)
	log << "Illumina Filtered\t";
    if (filter_adapter)
	log << "Adapter Seq" << "\t";
    log << "Low Quality\t"
	<< "Ambiguous Barcodes\t"
	<< "Ambiguous RAD-Tag\t"
	<< "Total\n";

    for (it = counters.begin(); it != counters.end(); it++) {
	log << it->first                 << "\t"
	    << it->second["retained"]    << "\t";
	if (filter_illumina)
	    log << it->second["ill_filtered"] << "\t";
	if (filter_adapter)
	    log << it->second["adapter"] << "\t";
	log << it->second["low_quality"] << "\t"
	    << it->second["ambiguous"]   << "\t"
	    << it->second["noradtag"]    << "\t"
	    << it->second["total"]       << "\n";
    }

    map<string, long> c;
    c["total"]        = 0;
    c["low_quality"]  = 0;
    c["adapter"]      = 0;
    c["ill_filtered"] = 0;
    c["ambiguous"]    = 0;
    c["noradtag"]     = 0;

    //
    // Total up the individual counters
    //
    for (it = counters.begin(); it != counters.end(); it++) {
	c["total"]        += it->second["total"];
	c["ill_filtered"] += it->second["ill_filtered"];
	c["low_quality"]  += it->second["low_quality"];
	c["adapter"]      += it->second["adapter"];
	c["ambiguous"]    += it->second["ambiguous"];
	c["noradtag"]     += it->second["noradtag"];
	c["retained"]     += it->second["retained"]; 
    }

    cerr << c["total"] << " total sequences;\n";
    if (filter_illumina)
	cerr << "  " << c["ill_filtered"] << " failed Illumina filtered reads;\n";
    if (filter_adapter)
	cerr << "  " << c["adapter"] << " reads contained adapter sequence;\n";
    cerr << "  " << c["ambiguous"]   << " ambiguous barcode drops;\n"
	 << "  " << c["low_quality"] << " low quality read drops;\n"
	 << "  " << c["noradtag"]    << " ambiguous RAD-Tag drops;\n"
	 << c["retained"] << " retained reads.\n";

    log	<< "\n" 
	 << "Total Sequences\t"      << c["total"]       << "\n";
    if (filter_illumina)
	log << "Failed Illumina filtered reads\t" << c["ill_filtered"] << "\n";
    if (filter_adapter)
	log << "Reads containing adapter sequence\t" << c["adapter"] << "\n";
    log << "Ambiguous Barcodes\t"   << c["ambiguous"]   << "\n"
	<< "Low Quality\t"          << c["low_quality"] << "\n"
	<< "Ambiguous RAD-Tag\t"    << c["noradtag"]    << "\n"
	<< "Retained Reads\t"       << c["retained"]    << "\n";

    if (bc_size_1 == 0) return 0;

    //
    // Print out barcode information.
    //
    log << "\n"
	<< "Barcode\t" 
	<< "Total\t"
	<< "No RadTag\t"
	<< "Retained\n";

    set<BarcodePair> barcode_list;

    for (uint i = 0; i < barcodes.size(); i++) {
	barcode_list.insert(barcodes[i]);

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
    // Sort unused barcodes by number of occurances.
    //
    map<BarcodePair, map<string, long> >::iterator bit;
    vector<pair<BarcodePair, int> > bcs;
    for (bit = barcode_log.begin(); bit != barcode_log.end(); bit++)
    	bcs.push_back(make_pair(bit->first, bit->second["total"]));
    sort(bcs.begin(), bcs.end(), compare_barcodes);

    for (uint i = 0; i < bcs.size(); i++) {
    	if (barcode_list.count(bcs[i].first)) continue;
    	if (bcs[i].second == 0) continue;

    	log << bcs[i].first << "\t"
    	    << bcs[i].second << "\n";
    }

    log.close();

    return 0;
}

int  compare_barcodes(pair<BarcodePair, int> a, pair<BarcodePair, int> b) {
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
	    {"merge",              no_argument, NULL, 'm'},
	    {"disable_rad_check",  no_argument, NULL, 'R'},
	    {"filter_illumina",    no_argument, NULL, 'F'},
	    {"index_null",         no_argument, NULL, 'u'},
	    {"inline_null",        no_argument, NULL, 'V'},
	    {"index_index",        no_argument, NULL, 'W'},
	    {"inline_inline",      no_argument, NULL, 'x'},
	    {"index_inline",       no_argument, NULL, 'Y'},
	    {"inline_index",       no_argument, NULL, 'Z'},
	    {"barcode_dist", required_argument, NULL, 'B'},
	    {"infile_type",  required_argument, NULL, 'i'},
	    {"outfile_type", required_argument, NULL, 'y'},
	    {"file",         required_argument, NULL, 'f'},
	    {"file_p1",      required_argument, NULL, '1'},
	    {"file_p2",      required_argument, NULL, '2'},
	    {"path",         required_argument, NULL, 'p'},
	    {"outpath",      required_argument, NULL, 'o'},
	    {"truncate",     required_argument, NULL, 't'},
	    {"renz_1",       required_argument, NULL, 'e'},
	    {"renz_2",       required_argument, NULL, 'z'},
	    {"barcodes",     required_argument, NULL, 'b'},
	    {"window_size",  required_argument, NULL, 'w'},
	    {"score_limit",  required_argument, NULL, 's'},
	    {"encoding",     required_argument, NULL, 'E'},
	    {"adapter_1",    required_argument, NULL, 'A'},
	    {"adapter_2",    required_argument, NULL, 'G'},
	    {"adapter_mm",   required_argument, NULL, 'T'},
	    {0, 0, 0, 0}
	};
	
	// getopt_long stores the option index here.
	int option_index = 0;

	c = getopt_long(argc, argv, "uVWxYZhvRFIcqrDPmB:i:y:f:o:t:e:z:b:1:2:p:s:w:E:A:G:T:", long_options, &option_index);

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
	    else if (strcasecmp(optarg, "gzfastq") == 0)
                in_file_type = gzfastq;
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
	    else {
		cerr << "Unknown quality score encoding, '" << optarg << "'\n";
		help();
	    }
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
	case 'B':
	    barcode_dist = is_integer(optarg);
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
	    truncate_seq = is_integer(optarg);
	    break;
	case 'e':
	    renz_1 = optarg;
	    break;
	case 'z':
	    renz_2 = optarg;
	    break;
	case 'b':
	    barcode_file = optarg;
	    if (barcode_type == null_null)
		barcode_type = inline_null;
	    break;
	case 'm':
	    merge = true;
	    break;
	case 'D':
	    discards = true;
	    break;
	case 'R':
	    check_radtag = false;
	    break;
	case 'F':
	    filter_illumina = true;
	    break;
	case 'u':
	    barcode_type = index_null;
	    break;
	case 'V':
	    barcode_type = inline_null;
	    break;
	case 'W':
	    barcode_type = index_index;
	    break;
	case 'x':
	    barcode_type = inline_inline;
	    break;
	case 'Y':
	    barcode_type = index_inline;
	    break;
	case 'Z':
	    barcode_type = inline_index;
	    break;
     	case 'A':
	    adapter_1 = new char[strlen(optarg) + 1];
	    strcpy(adapter_1, optarg);
	    filter_adapter = true;
	    break;
     	case 'G':
	    adapter_2 = new char[strlen(optarg) + 1];
	    strcpy(adapter_2, optarg);
	    filter_adapter = true;
	    break;
     	case 'T':
	    distance = is_integer(optarg);
	    break;
 	case 'w':
	    win_size = is_double(optarg);
	    break;
	case 's':
	    score_limit = is_integer(optarg);
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

    if (in_file_type == unknown)
	in_file_type = ftype;

    if (barcode_file.length() == 0)
	cerr << "No barcodes specified, files will not be demultiplexed.\n";

    if (barcode_file.length() > 0 && merge) {
	cerr << "You may specify a set of barcodes, or that all files should be merged, not both.\n";
	help();
    }

    if (check_radtag && renz_1.length() == 0) {
	cerr << "You must specify the restriction enzyme used.\n";
	help();
    }

    if (check_radtag && renz.count(renz_1) == 0) {
	cerr << "Unrecognized restriction enzyme specified: '" << renz_1.c_str() << "'.\n";
	help();
    }

    if (check_radtag && renz_2.length() > 0 && renz.count(renz_2) == 0) {
	cerr << "Unrecognized restriction enzyme specified: '" << renz_2.c_str() << "'.\n";
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
              << "process_radtags [-f in_file | -p in_dir [-P] | -1 pair_1 -2 pair_2] -b barcode_file -o out_dir -e enz [-c] [-q] [-r] [-t len] [-D] [-w size] [-s lim] [-h]\n"
	      << "  f: path to the input file if processing single-end seqeunces.\n"
	      << "  i: input file type, either 'bustard' for the Illumina BUSTARD output files, 'fastq', or 'gzfastq' for gzipped Fastq (default 'fastq').\n"
	      << "  p: path to a directory of files.\n"
	      << "  P: files contained within directory specified by '-p' are paired.\n"
	      << "  1: first input file in a set of paired-end sequences.\n"
	      << "  2: second input file in a set of paired-end sequences.\n"
	      << "  o: path to output the processed files.\n"
	      << "  y: output type, either 'fastq' or 'fasta' (default fastq).\n"
	      << "  b: path to a file containing barcodes for this run.\n"
	      << "  c: clean data, remove any read with an uncalled base.\n"
	      << "  q: discard reads with low quality scores.\n"
	      << "  r: rescue barcodes and RAD-Tags.\n"
	      << "  t: truncate final read length to this value.\n"
	      << "  E: specify how quality scores are encoded, 'phred33' (Illumina 1.8+, Sanger, default) or 'phred64' (Illumina 1.3 - 1.5).\n"
	      << "  D: capture discarded reads to a file.\n"
	      << "  w: set the size of the sliding window as a fraction of the read length, between 0 and 1 (default 0.15).\n"
	      << "  s: set the score limit. If the average score within the sliding window drops below this value, the read is discarded (default 10).\n"
	      << "  h: display this help messsage." << "\n\n"
	      << "  Barcode options:\n"
	      << "    --inline_null:   barcode is inline with sequence, occurs only on single-end read (default).\n"
	      << "    --index_null:    barcode is provded in FASTQ header, occurs only on single-end read.\n"
	      << "    --inline_inline: barcode is inline with sequence, occurs on single and paired-end read.\n"
	      << "    --index_index:   barcode is provded in FASTQ header, occurs on single and paired-end read.\n"
	      << "    --inline_index:  barcode is inline with sequence on single-end read, occurs in FASTQ header for paired-end read.\n"
	      << "    --index_inline:  barcode occurs in FASTQ header for single-end read, is inline with sequence on paired-end read.\n\n"
	      << "  Restriction enzyme options:\n"
	      << "    -e <enz>, --renz_1 <enz>: provide the restriction enzyme used (cut site occurs on single-end read)\n"
	      << "    --renz_2 <enz>: if a double digest was used, provide the second restriction enzyme used (cut site occurs on the paired-end read).\n"
	      << "    Currently supported enzymes include:\n"
	      << "      ";

    map<string, int>::iterator it;
    uint cnt = renz_cnt.size();
    it = renz_cnt.begin();
    for (uint i = 1; i <= cnt; i++) {
	std::cerr << "'" << it->first << "'";
	if (i < cnt - 1) 
	    std::cerr << ", ";
	else if (i == cnt - 1)
	    std::cerr << ", or ";

	if (i % 6 == 0)
	    std::cerr << "\n      ";

	it++;
    }

    std::cerr << ".\n\n" 
	      << "  Adapter options:\n"
	      << "    --adapter_1 <sequence>: provide adaptor sequence that may occur on the single-end read for filtering.\n"
	      << "    --adapter_2 <sequence>: provide adaptor sequence that may occur on the paired-read for filtering.\n"
	      << "      --adapter_mm <mismatches>: number of mismatches allowed in the adapter sequence.\n\n"
	      << "  Output options:\n"
	      << "    --merge: if no barcodes are specified, merge all input files into a single output file.\n\n"
	      << "  Advanced options:\n"
	      << "    --filter_illumina: discard reads that have been marked by Illumina's chastity/purity filter as failing.\n"
	      << "    --disable_rad_check: disable checking if the RAD site is intact.\n"
	      << "    --barcode_dist: provide the distace between barcodes to allow for barcode rescue (default 2)\n";

    exit(0);
}
