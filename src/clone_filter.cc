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
// clone_filter -- 
//
// Julian Catchen
// jcatchen@uoregon.edu
// University of Oregon
//
// $Id: clone_filter.cc 2099 2011-04-30 22:04:37Z catchen $
//

#include "clone_filter.h"

//
// Global variables to hold command-line options.
//
file_type in_file_type  = unknown;
file_type out_file_type = fastq;
string in_path_1;
string in_path_2;
string out_path;
int    barcode_size = 0;


int main (int argc, char* argv[]) {

    parse_command_line(argc, argv);

    Input *fh_1, *fh_2;

    cerr << "Reading data from:\n  " 
	 << in_path_1 << " and\n  " 
	 << in_path_2 << "\n";

    if (in_file_type == fastq) {
        fh_1 = new Fastq(in_path_1.c_str());
	fh_2 = new Fastq(in_path_2.c_str());
    } else if (in_file_type == fasta) {
        fh_1 = new Fasta(in_path_1.c_str());
        fh_2 = new Fasta(in_path_2.c_str());
    } else if (in_file_type == bustard) {
        fh_1 = new Bustard(in_path_1.c_str());
        fh_2 = new Bustard(in_path_2.c_str());
    }

    //
    // Open output files.
    //
    ofstream *ofh_1, *ofh_2;

    int    pos_1  = in_path_1.find_last_of("/");
    int    pos_2  = in_path_1.find_last_of(".");
    string path   = out_path + in_path_1.substr(pos_1 + 1, pos_2 - pos_1 - 1) + ".fil" + in_path_1.substr(pos_2);
    ofh_1 = new ofstream(path.c_str(), ifstream::out);
    if (ofh_1->fail()) {
	cerr << "Error opening output file '" << path << "'\n";
	exit(1);
    }
    cerr << "Writing data to:\n  " 
	 << path << " and\n  ";

    pos_1  = in_path_2.find_last_of("/");
    pos_2  = in_path_2.find_last_of(".");
    path   = out_path + in_path_2.substr(pos_1 + 1, pos_2 - pos_1 - 1) + ".fil" + in_path_2.substr(pos_2);
    ofh_2 = new ofstream(path.c_str(), ifstream::out);
    if (ofh_2->fail()) {
	cerr << "Error opening output file '" << path << "'\n";
	exit(1);
    }
    cerr << path << "\n";

    CloneHash clone_map;

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

    long  i = 1;
    bool  exists;
    char *hash_key;
    uint  seq_len   = strlen(s_1->seq);
    uint  tot_reads = 0;
    uint  red_reads = 0;

    do {
        if (i % 10000 == 0) cerr << "  Processing short read " << i << "       \r";

	tot_reads++;

	exists = clone_map.count(s_1->seq) == 0 ? false : true;

	if (exists) {
	    hash_key = s_1->seq;
	    clone_map[hash_key].pairs[s_2->seq]++;
	} else {
	    hash_key = new char [seq_len + 1];
	    strcpy(hash_key, s_1->seq);
	    clone_map[hash_key].id = s_1->id;
	    clone_map[hash_key].pairs[s_2->seq] = 1;
	}

	delete s_1;
	delete s_2;

	i++;
    } while ((s_1 = fh_1->next_seq()) != NULL && 
	     (s_2 = fh_2->next_seq()) != NULL);

    cerr << "\n";

    delete fh_1;
    delete fh_2;

    CloneHash::iterator it;
    map<string, int>::iterator j;
    uint k;

    cerr << "Writing filtered data...";

    for (it = clone_map.begin(); it != clone_map.end(); it++) {
	k = 1;

	for (j = it->second.pairs.begin(); j != it->second.pairs.end(); j++) {

	    *ofh_1 << ">" << it->second.id << "|1\n"
		   << it->first << "\n";
	    *ofh_2 << ">" << it->second.id << "|2-" << k << "|" << j->second << "\n"
		   << j->first << "\n";
	    k++;
	    red_reads++;
	}
    }

    cerr << "done.\n";

    delete ofh_1;
    delete ofh_2;

    char buf[32];
    sprintf(buf, "%0.2f%%", ((float) (tot_reads - red_reads) / (float) tot_reads) * 100);
    cerr << tot_reads << " pairs of reads input. " << red_reads << " pairs of reads output, " << buf << " clone reads.\n";

    return 0;
}

int parse_command_line(int argc, char* argv[]) {
    int c;
     
    while (1) {
	static struct option long_options[] = {
	    {"help",         no_argument,       NULL, 'h'},
            {"version",      no_argument,       NULL, 'v'},
	    {"infile_type",  required_argument, NULL, 'i'},
	    {"outfile_type", required_argument, NULL, 'y'},
	    {"file_p1",      required_argument, NULL, '1'},
	    {"file_p2",      required_argument, NULL, '2'},
	    {"outpath",      required_argument, NULL, 'o'},
	    {0, 0, 0, 0}
	};
	
	// getopt_long stores the option index here.
	int option_index = 0;

	c = getopt_long(argc, argv, "hvi:y:o:1:2:", long_options, &option_index);

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
	    else if (strcasecmp(optarg, "fasta") == 0)
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
	case '1':
	    in_path_1 = optarg;
	    break;
	case '2':
	    in_path_2 = optarg;
	    break;
	case 'o':
	    out_path = optarg;
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

    if (in_path_1.length() == 0 || in_path_2.length() == 0) {
	cerr << "You must specify a pair of input files.\n";
	help();
    }

    if (out_path.length() == 0) 
	out_path = ".";

    if (out_path.at(out_path.length() - 1) != '/') 
	out_path += "/";

    if (in_file_type == unknown)
	in_file_type = fastq;

    return 0;
}

void version() {
    std::cerr << "clone_filter " << VERSION << "\n\n";

    exit(0);
}

void help() {
    std::cerr << "clone_filter " << VERSION << "\n"
              << "clone_filter -1 pair_1 -2 pair_2 -o out_dir [-i type] [-y type] [-h]\n"
	      << "  1: first input file in a set of paired-end sequences.\n"
	      << "  2: second input file in a set of paired-end sequences.\n"
	      << "  i: input file type, either 'bustard' for the Illumina BUSTARD output files, or 'fastq' (default 'fastq').\n"
	      << "  o: path to output the processed files.\n"
	      << "  y: output type, either 'fastq' or 'fasta' (default fastq).\n"
	      << "  h: display this help messsage.\n";

    exit(0);
}
