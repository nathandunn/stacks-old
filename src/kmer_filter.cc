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
// kmer_filter -- 
//
// Julian Catchen
// jcatchen@uoregon.edu
// University of Oregon
//
// $Id: kmer_fileter.cc 2099 2011-04-30 22:04:37Z catchen $
//

#include "kmer_filter.h"

//
// Global variables to hold command-line options.
//
file_type in_file_type  = unknown;
file_type out_file_type = fastq;
vector<string> in_files;
vector<string> in_pair_files;
string in_path;
string out_path;
bool   discards     = false;
bool   record_kmers = false;
bool   kmer_distr   = false;
bool   ill_barcode  = false;
int    kmer_len     = 19;
int    min_k_freq   = 0;
int    max_k_freq   = 0;
int    min_lim      = 19;
int    max_lim      = 19;
int    num_threads  = 1;
int    barcode_size = 0;

int main (int argc, char* argv[]) {

    parse_command_line(argc, argv);

    vector<pair<string, string> > files, pair_files;
    map<string, map<string, long> > counters;
    SeqKmerHash kmers;

    build_file_list(in_files, files);
    cerr << "Found " << files.size() << " input file(s).\n";

    build_file_list(in_pair_files, pair_files);
    cerr << "Found " << pair_files.size() << " paired input file(s).\n";

    cerr << "Using a kmer size of " << kmer_len << "\n"
	 << "  Kmer is considered rare when it occurs " << min_k_freq << " or less times.\n"
	 << "  Kmer is considered abundant when it occurs " << max_k_freq << " or more times.\n"
	 << "  A read is dropped when it contains " << min_lim << " or more rare kmers.\n"
	 << "  A read is dropped when it contains " << max_lim << " or more abundant kmers.\n";

    populate_kmers(pair_files, files, kmers);

    if (record_kmers) write_rare_abundant_kmers(kmers);

    for (uint i = 0; i < pair_files.size(); i += 2) {
    	cerr << "Processing paired file " << i+1 << " of " << pair_files.size() << " [" << pair_files[i].second << "]\n";

    	counters[pair_files[i].second]["total"]      = 0;
    	counters[pair_files[i].second]["retained"]   = 0;
    	counters[pair_files[i].second]["rare_k"]     = 0;
    	counters[pair_files[i].second]["abundant_k"] = 0;

	process_paired_reads(pair_files[i].first,
			     pair_files[i].second, 
			     pair_files[i+1].first,
			     pair_files[i+1].second, 
			     kmers,
			     counters[pair_files[i].second]);

    	cerr <<	"  " 
    	     << counters[pair_files[i].second]["total"] << " total reads; "
    	     << "-" << counters[pair_files[i].second]["rare_k"] << " rare k-mer reads; "
    	     << "-" << counters[pair_files[i].second]["abundant_k"] << " abundant k-mer reads; "
    	     << counters[pair_files[i].second]["retained"] << " retained reads.\n";
    }

    for (uint i = 0; i < files.size(); i++) {
    	cerr << "Processing file " << i+1 << " of " << files.size() << " [" << files[i].second << "]\n";

    	counters[files[i].second]["total"]      = 0;
    	counters[files[i].second]["retained"]   = 0;
    	counters[files[i].second]["rare_k"]     = 0;
    	counters[files[i].second]["abundant_k"] = 0;

	process_reads(files[i].first,
		      files[i].second, 
		      kmers,
		      counters[files[i].second]);

    	cerr <<	"  " 
    	     << counters[files[i].second]["total"] << " total reads; "
    	     << "-" << counters[files[i].second]["rare_k"] << " rare k-mer reads; "
    	     << "-" << counters[files[i].second]["abundant_k"] << " abundant k-mer reads; "
    	     << counters[files[i].second]["retained"] << " retained reads.\n";
    }

    print_results(counters);

    return 0;
}

int process_paired_reads(string in_path_1, 
			 string in_file_1,
			 string in_path_2, 
			 string in_file_2,
			 SeqKmerHash &kmers,
			 map<string, long> &counter) {
    Input *fh_1, *fh_2;
    ofstream *discard_fh_1, *discard_fh_2;

    string path_1 = in_path_1 + in_file_1;
    string path_2 = in_path_2 + in_file_2;

    if (in_file_type == fastq) {
        fh_1 = new Fastq(path_1.c_str());
        fh_2 = new Fastq(path_2.c_str());
    } else if (in_file_type == fasta) {
        fh_1 = new Fasta(path_1.c_str());
        fh_2 = new Fasta(path_2.c_str());
    } else if (in_file_type == bustard) {
        fh_1 = new Bustard(path_1.c_str());
        fh_2 = new Bustard(path_2.c_str());
    }

    //
    // Open the output files.
    //
    int       pos_1  = in_file_1.find_last_of(".");
    string    path   = out_path + in_file_1.substr(0, pos_1) + ".fil" + in_file_1.substr(pos_1);
    ofstream *ofh_1  = new ofstream(path.c_str(), ifstream::out);
    pos_1 = in_file_2.find_last_of(".");
    path  = out_path + in_file_2.substr(0, pos_1) + ".fil" + in_file_2.substr(pos_1);
    ofstream *ofh_2  = new ofstream(path.c_str(), ifstream::out);
    path  = out_path + in_file_2.substr(0, pos_1) + ".fil.rem";
    path += out_file_type == fastq ? ".fq" : ".fa";
    ofstream *rem_fh = new ofstream(path.c_str(), ifstream::out);

    //
    // Open a file for recording discarded reads
    //
    if (discards) {
	int pos = in_file_1.find_last_of(".");
	path = out_path + in_file_1.substr(0, pos) + ".discards" + in_file_1.substr(pos);
	discard_fh_1 = new ofstream(path.c_str(), ifstream::out);

	if (discard_fh_1->fail()) {
	    cerr << "Error opening discard output file '" << path << "'\n";
	    exit(1);
	}
	pos  = in_file_2.find_last_of(".");
	path = out_path + in_file_2.substr(0, pos) + ".discards" + in_file_2.substr(pos);
	discard_fh_2 = new ofstream(path.c_str(), ifstream::out);

	if (discard_fh_2->fail()) {
	    cerr << "Error opening discard output file '" << path << "'\n";
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

    long i = 1;
    do {
	if (i % 10000 == 0) cerr << "  Processing short read " << i << "       \r";
	counter["total"] += 2;
	bool retain_1 = true;
	bool retain_2 = true;
	stringstream msg_1, msg_2;

	int rare_k     = 0; 
	int abundant_k = 0;
	char *kmer     = new char[kmer_len + 1];
	kmer[kmer_len] = '\0';
	int num_kmers  = strlen(s_1->seq) - kmer_len + 1;

	//
	// Drop the first sequence if it has too many rare or abundant kmers.
	//
	kmer_lookup(kmers, s_1->seq, kmer, num_kmers, &rare_k, &abundant_k);

	if (min_k_freq > 0 && min_lim && rare_k > min_lim) {
	    counter["rare_k"]++;
	    retain_1 = false;
	    msg_1 << "rare_k_" << rare_k;
	}

	if (retain_1 && max_k_freq > 0 && max_lim && abundant_k > max_lim) {
	    counter["abundant_k"]++;
	    retain_1 = false;
	    msg_1 << "abundant_k_" << abundant_k;
	}

	//
	// Drop the second sequence if it has too many rare or abundant kmers.
	//
	rare_k     = 0;
	abundant_k = 0;
	num_kmers  = strlen(s_2->seq) - kmer_len + 1;
	kmer_lookup(kmers, s_2->seq, kmer, num_kmers, &rare_k, &abundant_k);

	if (min_k_freq > 0 && min_lim && rare_k > min_lim) {
	    counter["rare_k"]++;
	    retain_2 = false;
	    msg_2 << "rare_k_" << rare_k;
	}

	if (retain_2 && max_k_freq > 0 && max_lim && abundant_k > max_lim) {
	    counter["abundant_k"]++;
	    retain_2 = false;
	    msg_2 << "abundant_k_" << abundant_k;
	}

	if (retain_1 && retain_2) {
	    counter["retained"] += 2;
	    out_file_type == fastq ? 
	 	write_fastq(ofh_1, s_1) : write_fasta(ofh_1, s_1);
	    out_file_type == fastq ? 
	 	write_fastq(ofh_2, s_2) : write_fasta(ofh_2, s_2);
	}

	if (retain_1 && !retain_2) {
	    counter["retained"]++;
	    out_file_type == fastq ? 
	 	write_fastq(rem_fh, s_1) : write_fasta(rem_fh, s_1);
	}

	if (!retain_1 && retain_2) {
	    counter["retained"]++;
	    out_file_type == fastq ? 
	 	write_fastq(rem_fh, s_2) : write_fasta(rem_fh, s_2);
	}

	if (discards && !retain_1)
	    out_file_type == fastq ? 
		write_fastq(discard_fh_1, s_1, msg_1.str()) : write_fasta(discard_fh_1, s_1, msg_1.str());
	if (discards && !retain_2)
	    out_file_type == fastq ? 
		write_fastq(discard_fh_2, s_2, msg_2.str()) : write_fasta(discard_fh_2, s_2, msg_2.str());

	delete s_1;
	delete s_2;

	i++;
    } while ((s_1 = fh_1->next_seq()) != NULL &&
	     (s_2 = fh_2->next_seq()) != NULL);

    if (discards) {
	delete discard_fh_1;
	delete discard_fh_2;
    }

    //
    // Close the file and delete the Input object.
    //
    delete fh_1;
    delete fh_2;
    delete ofh_1;
    delete ofh_2;
    delete rem_fh;

    return 0;
}

int process_reads(string in_path, 
		  string in_file,
		  SeqKmerHash &kmers,
		  map<string, long> &counter) {
    Input *fh;
    ofstream *discard_fh;

    string path = in_path + in_file;

    if (in_file_type == fastq)
        fh = new Fastq(path.c_str());
    else if (in_file_type == fasta)
        fh = new Fasta(path.c_str());
    else if (in_file_type == bustard)
        fh = new Bustard(path.c_str());

    //
    // Open the output file.
    //
    int pos = in_file.find_last_of(".");
    path    = out_path + in_file.substr(0, pos) + ".fil" + in_file.substr(pos);
    ofstream *out_fh = new ofstream(path.c_str(), ifstream::out);    

    //
    // Open a file for recording discarded reads
    //
    if (discards) {
	int pos = in_file.find_last_of(".");
	path = out_path + in_file.substr(0, pos) + ".discards" + in_file.substr(pos);
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

    long i = 1;
    do {
	if (i % 10000 == 0) cerr << "  Processing short read " << i << "       \r";
	counter["total"]++;
	bool retain = true;
	stringstream msg;

	//
	// Drop this sequence if it has too many rare or abundant kmers.
	//
	int rare_k, abundant_k, num_kmers;
	char *kmer     = new char[kmer_len + 1];
	kmer[kmer_len] = '\0';
	num_kmers      = strlen(s->seq) - kmer_len + 1;

	kmer_lookup(kmers, s->seq, kmer, num_kmers, &rare_k, &abundant_k);

	if (min_k_freq > 0 && min_lim && rare_k > min_lim) {
	    counter["rare_k"]++;
	    retain = false;
	    msg << "rare_k_" << rare_k;
	}

	if (retain && max_k_freq > 0 && max_lim && abundant_k > max_lim) {
	    counter["abundant_k"]++;
	    retain = false;
	    msg << "abundant_k_" << abundant_k;
	}

	if (retain) {
	    counter["retained"]++;
	    out_file_type == fastq ? 
	 	write_fastq(out_fh, s) : write_fasta(out_fh, s);
	}

	if (discards && !retain)
	    out_file_type == fastq ? 
		write_fastq(discard_fh, s, msg.str()) : write_fasta(discard_fh, s, msg.str());

	delete s;

	i++;
    } while ((s = fh->next_seq()) != NULL);

    if (discards) delete discard_fh;

    //
    // Close the file and delete the Input object.
    //
    delete fh;
    delete out_fh;

    return 0;
}

int populate_kmers(vector<pair<string, string> > &pair_files, 
		   vector<pair<string, string> > &files, 
		   SeqKmerHash &kmers) {
    //
    // Break each read down into k-mers and create a hash map of those k-mers
    // recording in which sequences they occur.
    //
    uint j   = 1;
    uint cnt = files.size() + pair_files.size();
    for (uint i = 0; i < files.size(); i++) {
	cerr << "Generating kmers from file " << j << " of " << cnt << " [" << files[i].second << "]\n";
	process_file_kmers(files[i].first + files[i].second, kmers);
	j++;
    }

    for (uint i = 0; i < pair_files.size(); i++) {
	cerr << "Generating kmers from file " << j << " of " << cnt << " [" << pair_files[i].second << "]\n";
	process_file_kmers(pair_files[i].first + pair_files[i].second, kmers);
	j++;
    }

    cerr << kmers.size() << " unique k-mers recorded.\n";

    if (kmer_distr) {
	generate_kmer_dist(kmers);
	exit(0);
    }

    return 0;
}

int write_rare_abundant_kmers(SeqKmerHash &kmer_map) {

    cerr << "Writing rare and abundant kmers...";

    string   path     = out_path + "rare_kmers.tsv";
    ofstream *rare_fh = new ofstream(path.c_str(), ifstream::out);
    if (rare_fh->fail()) {
	cerr << "Error opening rare kmer output file '" << path << "'\n";
	exit(1);
    }

    path = out_path + "abundant_kmers.tsv";
    ofstream *abun_fh = new ofstream(path.c_str(), ifstream::out);
    if (abun_fh->fail()) {
	cerr << "Error opening abundant kmer output file '" << path << "'\n";
	exit(1);
    }

    SeqKmerHash::iterator i;

    *rare_fh << "Kmer\tCount\n";
    *abun_fh << "Kmer\tCount\n";

    for (i = kmer_map.begin(); i != kmer_map.end(); i++) {

	if (i->second <= min_k_freq)
	    *rare_fh << i->first << "\t" << i->second << "\n";

	if (i->second >= max_k_freq)
	    *abun_fh << i->first << "\t" << i->second << "\n";
    }

    delete rare_fh;
    delete abun_fh;

    cerr << "done.\n";

    return 0;
}

int process_file_kmers(string path, SeqKmerHash &kmer_map) {
    vector<char *> kmers;
    char          *hash_key;
    bool           exists;
    int            j;
    Input         *fh;

    if (in_file_type == fastq)
        fh = new Fastq(path.c_str());
    else if (in_file_type == fasta)
        fh = new Fasta(path.c_str());
    else if (in_file_type == bustard)
        fh = new Bustard(path.c_str());

    //
    // Read in the first record, initializing the Seq object s.
    //
    Seq *s = fh->next_seq();
    if (s == NULL) {
    	cerr << "Unable to allocate Seq object.\n";
    	exit(1);
    }

    int   num_kmers;
    char *kmer = new char [kmer_len + 1];

    long i = 1;
    do {
	if (i % 10000 == 0) cerr << "  Processing short read " << i << "       \r";

	num_kmers = strlen(s->seq) - kmer_len + 1;

	//
	// Generate and hash the kmers for this raw read
	//
	kmer[kmer_len] = '\0';

	for (j = 0; j < num_kmers; j++) {
	    strncpy(kmer, s->seq + j, kmer_len);

	    exists = kmer_map.count(kmer) == 0 ? false : true;

	    if (exists) {
	    	hash_key = kmer;
	    } else {
	    	hash_key = new char [kmer_len + 1];
	    	strcpy(hash_key, kmer);
	    }

	    kmer_map[hash_key]++;
	}

	delete s;

	i++;
    } while ((s = fh->next_seq()) != NULL);

    //
    // Close the file and delete the Input object.
    //
    delete fh;

    return 0;
}

int generate_kmer_dist(SeqKmerHash &kmer_map) {
    SeqKmerHash::iterator i;
    map<uint, uint> bins;

    cerr << "Generating kmer distribution...\n";

    for (i = kmer_map.begin(); i != kmer_map.end(); i++)
	bins[i->second]++;

    map<uint, uint>::iterator j;
    vector<pair<uint, uint> > sorted_kmers;

    for (j = bins.begin(); j != bins.end(); j++)
	sorted_kmers.push_back(make_pair(j->first, j->second));

    cout << "KmerFrequency\tCount\n";

    for (unsigned long k = 0; k < sorted_kmers.size(); k++)
	cout << sorted_kmers[k].first << "\t" << sorted_kmers[k].second << "\n";

    return 0;
}

int kmer_map_cmp(pair<char *, long> a, pair<char *, long> b) {
    return (a.second < b.second);
}

inline
int kmer_lookup(SeqKmerHash &kmer_map, 
		char *read, char *kmer, 
		int num_kmers, 
		int *rare_k, int *abundant_k) {
    //
    // Generate and hash the kmers for this raw read
    //
    *rare_k        = 0;
    *abundant_k    = 0;
    kmer[kmer_len] = '\0';
    int cnt        = 0;

    for (int j = 0; j < num_kmers; j++) {
	strncpy(kmer, read + j, kmer_len);

	cnt = kmer_map[kmer];

	if (cnt <= min_k_freq) (*rare_k)++;
	if (cnt >= max_k_freq) (*abundant_k)++;
    }

    return 0;
}

int print_results(map<string, map<string, long> > &counters) {
    map<string, map<string, long> >::iterator it;

    string log_path = out_path + "kmer_filter.log";
    ofstream log(log_path.c_str());

    if (log.fail()) {
	cerr << "Unable to open log file '" << log_path << "'\n";
	return 0;
    }

    cerr << "Outputing details to log: '" << log_path << "'\n\n";

    log << "File\t"
	<< "Retained Reads\t"
	<< "Rare K\t"
	<< "Abundant K\t"
	<< "Total\n";

    for (it = counters.begin(); it != counters.end(); it++) {
	log << it->first                 << "\t"
	    << it->second["retained"]    << "\t"
	    << it->second["rare_k"]      << "\t"
	    << it->second["abundant_k"]  << "\t"
	    << it->second["total"]       << "\n";
    }

    map<string, long> c;
    c["total"]       = 0;

    //
    // Total up the individual counters
    //
    for (it = counters.begin(); it != counters.end(); it++) {
	c["total"]       += it->second["total"];
	c["retained"]    += it->second["retained"];
 	c["rare_k"]      += it->second["rare_k"];
 	c["abundant_k"]  += it->second["abundant_k"];
    }

    cerr << 
	c["total"] << " total sequences;\n"
	 << "  " << c["rare_k"]      << " rare k-mer reads;\n"
	 << "  " << c["abundant_k"]  << " abundant k-mer reads;\n"
	 << c["retained"] << " retained reads.\n";

    log	<< "Total Sequences\t"      << c["total"]       << "\n"
	<< "Retained Reads\t"       << c["retained"]      << "\n";

    log.close();

    return 0;
}

int build_file_list(vector<string> &in_files, vector<pair<string, string> > &files) {
    string file, suffix;
    int    pos;

    //
    // Scan a directory for a list of files.
    //
    if (in_path.length() > 0) {
	struct dirent *direntry;

	DIR *dir = opendir(in_path.c_str());

	if (dir == NULL) {
	    cerr << "Unable to open directory '" << in_path << "' for reading.\n";
	    exit(1);
	}

	while ((direntry = readdir(dir)) != NULL) {
	    file = direntry->d_name;

	    if (file == "." || file == "..")
		continue;

	    //
	    // Check that the file has the right suffix.
	    //
	    pos     = file.find_last_of(".");
	    suffix = file.substr(pos + 1);
	    if (in_file_type == fastq && (suffix.substr(0, 2) == "fq" || suffix.substr(0, 5) == "fastq"))
		files.push_back(make_pair(in_path, file));
	    else if (in_file_type == fasta && (suffix.substr(0, 2) == "fa" || suffix.substr(0, 5) == "fasta"))
		files.push_back(make_pair(in_path, file));
	}

	if (files.size() == 0)
	    cerr << "Unable to locate any input files to process within '" << in_path << "'\n";

    } else {
	string path;

	for (uint i = 0; i < in_files.size(); i++) {
	    //
	    // Files specified directly:
	    //    Break off file path and store path and file name.
	    //
 	    file = in_files[i];
	    pos  = file.find_last_of("/");
	    path = file.substr(0, pos + 1);
	    files.push_back(make_pair(path, file.substr(pos+1)));
	}
    }

    return 0;
}

int parse_command_line(int argc, char* argv[]) {
    string pair_1, pair_2;
    int c;
     
    while (1) {
	static struct option long_options[] = {
	    {"help",         no_argument,       NULL, 'h'},
            {"version",      no_argument,       NULL, 'v'},
	    {"discards",     no_argument,       NULL, 'D'},
	    {"record_kmers", no_argument,       NULL, 'R'},
	    {"pair_1",       required_argument, NULL, '1'},
	    {"pair_2",       required_argument, NULL, '2'},
	    {"infile_type",  required_argument, NULL, 'i'},
	    {"outfile_type", required_argument, NULL, 'y'},
	    {"file",         required_argument, NULL, 'f'},
	    {"path",         required_argument, NULL, 'p'},
	    {"outpath",      required_argument, NULL, 'o'},
	    {"k_dist",       no_argument,       NULL, 'I'},
	    {"k_len",        required_argument, NULL, 'K'},
	    {"min_k_freq",   required_argument, NULL, 'm'},
	    {"max_k_freq",   required_argument, NULL, 'M'},
	    {"min_lim",      required_argument, NULL, 'F'},
	    {"max_lim",      required_argument, NULL, 'G'},
	    {0, 0, 0, 0}
	};
	
	// getopt_long stores the option index here.
	int option_index = 0;

	c = getopt_long(argc, argv, "hvRDkI::K:F:G:M:m:i:y:f:o:t:p:1:2:", long_options, &option_index);

	// Detect the end of the options.
	if (c == -1)
	    break;
     
	switch (c) {
	case 'h':
	    help();
	    break;
     	case 'i':
             if (strcasecmp(optarg, "fasta") == 0)
                in_file_type = fasta;
            else
                in_file_type = fastq;
	    break;
     	case 'y':
            if (strcasecmp(optarg, "fasta") == 0)
                out_file_type = fasta;
	    else 
		out_file_type = fastq;
	    break;
     	case 'f':
	    in_files.push_back(optarg);
	    break;
	case '1':
	    pair_1 = optarg;
	    break;
	case '2':
	    pair_2 = optarg;
	    if (pair_1.length() == 0) help();
	    in_pair_files.push_back(pair_1);
	    in_pair_files.push_back(pair_2);
	    pair_1 = "";
	    pair_2 = "";
	    break;
	case 'p':
	    in_path = optarg;
	    break;
	case 'o':
	    out_path = optarg;
	    break;
	case 'D':
	    discards = true;
	    break;
	case 'I':
	    kmer_distr = true;
	    break;
	case 'K':
	    kmer_len = atoi(optarg);
	    break;
	case 'm':
	    min_k_freq = atoi(optarg);
	    break;
	case 'M':
	    max_k_freq = atoi(optarg);
	    break;
	case 'F':
	    min_lim = atoi(optarg);
	    break;
	case 'G':
	    max_lim = atoi(optarg);
	    break;
	case 'R':
	    record_kmers = true;
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

    if (in_files.size() == 0 && in_pair_files.size() == 0 && in_path.length() == 0) {
	cerr << "You must specify an input file of a directory path to a set of input files.\n";
	help();
    }

    if (in_files.size() > 0 && in_path.length() > 0) {
	cerr << "You must specify either a single input file (-f) or a directory path (-p), not both.\n";
	help();
    }

    if (in_path.length() > 0 && in_path.at(in_path.length() - 1) != '/') 
	in_path += "/";

    if (out_path.length() == 0) 
	out_path = ".";

    if (out_path.at(out_path.length() - 1) != '/') 
	out_path += "/";

    if (in_file_type == unknown)
	in_file_type = fastq;

    return 0;
}

void version() {
    std::cerr << "kmer_filter " << VERSION << "\n\n";

    exit(0);
}

void help() {
    std::cerr << "kmer_filter " << VERSION << "\n"
              << "kmer_filter [-f in_file_1 [-f in_file_2...] | -p in_dir] [-1 pair_1 -2 pair_2 [-1 pair_1...]] -o out_dir [-i type] [-y type] [-D] [-h]\n"
	      << "  f: path to the input file if processing single-end seqeunces.\n"
	      << "  i: input file type, either 'bustard' for the Illumina BUSTARD output files, 'fasta', or 'fastq' (default 'fastq').\n"
	      << "  p: path to a directory of files (for single-end files only).\n"
	      << "  1: specify the first in a pair of files to be processed together.\n"
	      << "  2: specify the second in a pair of files to be processed together.\n"
	      << "  o: path to output the processed files.\n"
	      << "  y: output type, either 'fastq' or 'fasta' (default fastq).\n"
	      << "  D: capture discarded reads to a file.\n"
	      << "  h: display this help messsage.\n\n"
	      << "  --k_dist: print k-mer frequency distribution and exit.\n"
	      << "  --k_len: specify k-mer size (default 19).\n"
	      << "  --min_k_freq: specify minimum limit for a kmer to be considered rare.\n"
	      << "  --max_k_freq: specify maximum limit for a kmer to be considered abundant.\n"
	      << "  --min_lim: specify number of rare kmers required to discard a read.\n"
	      << "  --max_lim: specify number of abundant kmers required to discard a read.\n"
	      << "  --record_kmers: record kmers above/below rare and abundant limits.\n"
	      << "\n";

    exit(0);
}
