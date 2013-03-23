// -*-mode:c++; c-style:k&r; c-basic-offset:4;-*-
//
// Copyright 2012, Julian Catchen <jcatchen@uoregon.edu>
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
// file_io.cc -- common routines for opening groups of files and processing barcode lists.
//
// Julian Catchen
// jcatchen@uoregon.edu
// University of Oregon
//
// $Id: file_io.cc 2146 2011-08-02 22:11:50Z catchen $
//

#include "file_io.h"

int open_files(vector<pair<string, string> > &files,
	       vector<BarcodePair> &barcodes, 
	       map<BarcodePair, ofstream *> &pair_1_fhs, 
	       map<BarcodePair, ofstream *> &pair_2_fhs, 
	       map<BarcodePair, ofstream *> &rem_1_fhs, 
	       map<BarcodePair, ofstream *> &rem_2_fhs, 
	       map<string, map<string, long> > &counters) {
    string path, suffix_1, suffix_2;

    if (paired && interleave == false) {
	suffix_1 = ".1";
	suffix_2 = ".2";
    }
    if (out_file_type == fastq) {
	suffix_1 += ".fq";
	suffix_2 += ".fq";
    } else {
	suffix_1 += ".fa";
	suffix_2 += ".fa";
    }

    ofstream   *fh;
    BarcodePair bc;
    //
    // If the size of the barcodes vector is 0, then no barcodes
    // were submitted. In this case, we want to open output files
    // of the same name as input files, but in out_path.
    //
    if (barcodes.size() == 0 && merge == false) {

	struct stat sb_1, sb_2;

	for (uint i = 0; i < files.size(); i++) {

	    bc.se = files[i].first;
	    if (paired)
		bc.pe = files[i].second;

	    path = out_path + files[i].first;

	    if (stat((in_path_1 + files[i].first).c_str(), &sb_1) == -1) {
		cerr << "Unable to stat input file " << in_path_1 + files[i].first << "\n";
		exit(1);
	    }
	    if (stat(path.c_str(), &sb_2) == 0 &&
		sb_2.st_dev == sb_1.st_dev     &&
		sb_2.st_ino == sb_1.st_ino) {
		cerr << "Input and output files ('" << path << "') are the same and will cause the input " 
		     << "file to be overwritten. Please specify a separate output directory using '-o'.\n";
		help();
	    }

	    fh = new ofstream(path.c_str(), ifstream::out);
            pair_1_fhs[bc] = fh;

	    if (pair_1_fhs[bc]->fail()) {
		cerr << "Error opening output file '" << path << "'\n";
		exit(1);
	    }

            if (paired) {
		path = out_path + files[i].second;

		if (stat((in_path_2 + files[i].second).c_str(), &sb_1) == -1) {
		    cerr << "Unable to stat input file " << in_path_2 + files[i].second << "\n";
		    exit(1);
		}
		if (stat(path.c_str(), &sb_2) == 0 &&
		    sb_2.st_dev == sb_1.st_dev     &&
		    sb_2.st_ino == sb_1.st_ino) {
		    cerr << "Input and output file names ('" << path << "') are the same and will cause the input " 
			 << "file to be overwritten. Please specify a separate output directory using '-o'.\n";
		    help();
		}

		fh = new ofstream(path.c_str(), ifstream::out);
                pair_2_fhs[bc] = fh;

		if (pair_2_fhs[bc]->fail()) {
		    cerr << "Error opening output file '" << path << "'\n";
		    exit(1);
		}

		int pos = files[i].first.find_last_of(".");
		path = out_path + files[i].first.substr(0, pos) + ".rem" + files[i].first.substr(pos);

		fh = new ofstream(path.c_str(), ifstream::out);
                rem_1_fhs[bc] = fh;

 		if (rem_1_fhs[bc]->fail()) {
		    cerr << "Error opening remainder output file '" << path << "'\n";
		    exit(1);
		}

		pos  = files[i].second.find_last_of(".");
		path = out_path + files[i].second.substr(0, pos) + ".rem" + files[i].second.substr(pos);

		fh = new ofstream(path.c_str(), ifstream::out);
		rem_2_fhs[bc] = fh;

		if (rem_2_fhs[bc]->fail()) {
		    cerr << "Error opening remainder output file '" << path << "'\n";
		    exit(1);
		}
            }
	}

	return 0;

    } else if (barcodes.size() == 0 && merge == true) {

	path = out_path + "sample_unbarcoded" + suffix_1;
	fh   = new ofstream(path.c_str(), ifstream::out);

	if (fh->fail()) {
	    cerr << "Error opening output file '" << path << "'\n";
	    exit(1);
	}

	for (uint i = 0; i < files.size(); i++) {
	    bc.se = files[i].first;
	    if (paired)
		bc.pe = files[i].second;

            pair_1_fhs[bc] = fh;
	}

	if (paired) {
	    path = out_path + "sample_unbarcoded" + suffix_2;
	    fh   = new ofstream(path.c_str(), ifstream::out);

	    if (fh->fail()) {
		cerr << "Error opening output file '" << path << "'\n";
		exit(1);
	    }

	    for (uint i = 0; i < files.size(); i++) {
		bc.se = files[i].first;
		if (paired) bc.pe = files[i].second;

		pair_2_fhs[bc] = fh;
	    }

	    path = out_path + "sample_unbarcoded.rem" + suffix_1;
	    fh   = new ofstream(path.c_str(), ifstream::out);

	    if (fh->fail()) {
		cerr << "Error opening remainder output file '" << path << "'\n";
		exit(1);
	    }

	    for (uint i = 0; i < files.size(); i++) {
		rem_1_fhs[bc] = fh;
	    }

	    path = out_path + "sample_unbarcoded.rem" + suffix_2;
	    fh   = new ofstream(path.c_str(), ifstream::out);

	    if (fh->fail()) {
		cerr << "Error opening remainder output file '" << path << "'\n";
		exit(1);
	    }

	    for (uint i = 0; i < files.size(); i++) {
		rem_2_fhs[bc] = fh;
	    }
	}

	return 0;
    }

    for (uint i = 0; i < barcodes.size(); i++) {

        if (interleave == true) {
	    path = out_path + "sample_" + barcodes[i].str() + suffix_1;
	    fh = new ofstream(path.c_str());
            pair_1_fhs[barcodes[i]] = fh;
 
	    if (pair_1_fhs[barcodes[i]]->fail()) {
		cerr << "Error opening output file '" << path << "'\n";
		exit(1);
	    }

        } else {
	    path = out_path + "sample_" + barcodes[i].str() + suffix_1;
	    fh = new ofstream(path.c_str(), ifstream::out);
            pair_1_fhs[barcodes[i]] = fh;

	    if (pair_1_fhs[barcodes[i]]->fail()) {
		cerr << "Error opening output file '" << path << "'\n";
		exit(1);
	    }

            if (paired) {
		path = out_path + "sample_" + barcodes[i].str() + suffix_2;
		fh = new ofstream(path.c_str(), ifstream::out);
                pair_2_fhs[barcodes[i]] = fh;

		if (pair_2_fhs[barcodes[i]]->fail()) {
		    cerr << "Error opening output file '" << path << "'\n";
		    exit(1);
		}

		path = out_path + "sample_" + barcodes[i].str() + ".rem" + suffix_1;
		fh = new ofstream(path.c_str(), ifstream::out);
                rem_1_fhs[barcodes[i]] = fh;

		if (rem_1_fhs[barcodes[i]]->fail()) {
		    cerr << "Error opening remainder output file '" << path << "'\n";
		    exit(1);
		}

		path = out_path + "sample_" + barcodes[i].str() + ".rem" + suffix_2;
		fh = new ofstream(path.c_str(), ifstream::out);
                rem_2_fhs[barcodes[i]] = fh;

		if (rem_2_fhs[barcodes[i]]->fail()) {
		    cerr << "Error opening remainder output file '" << path << "'\n";
		    exit(1);
		}
            }
        }
    }

    return 0;
}

int 
close_file_handles(map<BarcodePair, ofstream *> &fhs) 
{
    map<BarcodePair, ofstream*>::iterator i;
    set<ofstream*> ptrs;
    set<ofstream*>::iterator j;

    for (i = fhs.begin(); i != fhs.end(); i++) {
	i->second->close();
	ptrs.insert(i->second);
    }

    for (j = ptrs.begin(); j != ptrs.end(); j++) {
	delete *j;
    }

    return 0;
}

int 
load_barcodes(string barcode_file, vector<BarcodePair> &barcodes, 
	      set<string> &se_bc, set<string> &pe_bc,
	      int &se_len, int &pe_len) 
{
    switch(barcode_type) {
    case null_null:
	cerr << "Barcode type unspecified, assuming unbarcoded data.\n";
	break;
    case index_null:
	cerr << "Searching for single-end, indexed barcodes.\n";
	break;
    case inline_null:
	cerr << "Searching for single-end, inlined barcodes.\n";
	break;
    case index_index:
	cerr << "Searching for single and paired-end, indexed barcodes.\n";
	break;
    case inline_inline:
	cerr << "Searching for single and paired-end, inlined barcodes.\n";
	break;
    case inline_index:
	cerr << "Searching for single-end, inlined and paired-end, indexed barcodes.\n";
	break;
    case index_inline:
	cerr << "Searching for single-end, indexed and paired-end, inlined barcodes.\n";
	break;
    }

    if (barcode_file.length() == 0)
	return 0;

    char     line[id_len];
    ifstream fh(barcode_file.c_str(), ifstream::in);

    if (fh.fail()) {
        cerr << "Error opening barcode file '" << barcode_file << "'\n";
	exit(1);
    }

    char *p, *q, *r;

    while (fh.good()) {
	memset(line, 0, id_len);
	fh.getline(line, id_len);

	if (strlen(line) == 0) continue;

	//
	// Identify the first barcode and check that it's legitimate.
	//
	p = line;
	q = p;
	while (*q != '\0') {
	    switch (*q) {
	    case 'A':
	    case 'C':
	    case 'G':
	    case 'T':
		break;
	    case 'a':
		*q = 'A';
		break;
	    case 'c':
		*q = 'C';
		break;
	    case 'g':
		*q = 'G';
		break;
	    case 't':
		*q = 'T';
		break;
	    case '\r':
	    case '\t':
		*q = '\0';
		break;
	    default:
		cerr << "Invalid barcode: '" << line << "'\n";
		exit(1);
	    }
	    if (*q != '\0') q++;
	}

	//
	// Identify the second barcode and check that it's legitimate.
	//
	if (q - p < id_len)
	    q++;
	r = q;
	while (*q != '\0') {
	    switch (*q) {
	    case 'A':
	    case 'C':
	    case 'G':
	    case 'T':
		break;
	    case 'a':
		*q = 'A';
		break;
	    case 'c':
		*q = 'C';
		break;
	    case 'g':
		*q = 'G';
		break;
	    case 't':
		*q = 'T';
		break;
	    case '\r':
	    case '\t':
		*q = '\0';
		break;
	    default:
		cerr << "Invalid barcode: '" << line << "'\n";
		exit(1);
	    }
	    if (*q != '\0') q++;
	}

	barcodes.push_back(BarcodePair(p, r));
	if (strlen(p) > 0) se_bc.insert(string(p));
	if (strlen(r) > 0) pe_bc.insert(string(r));
    }

    fh.close();

    if (barcodes.size() == 0) {
	cerr << "Unable to load any barcodes from '" << barcode_file << "'\n";
	help();
    }

    //
    // Make sure barcodes are properly paired up.
    //
    int pe_cnt = 0;
    int se_cnt = 0;
    for (uint i = 0; i < barcodes.size(); i++) {
	se_cnt += (barcodes[i].se.length() > 0) ? 1 : 0;
	pe_cnt += (barcodes[i].pe.length() > 0) ? 1 : 0;
    }

    if (pe_cnt > 0 && se_cnt != pe_cnt) {
	cerr << "Single and paired-end barcodes must be properly paired.\n";
	help();
    }

    //
    // Determine the barcode length for the single-end barcodes
    //
    int prev;
    prev = barcodes[0].se.length();
    for (uint i = 0; i < barcodes.size(); i++) {
        se_len = barcodes[i].se.length();

        if (prev != se_len) {
            cerr << "Single-end barcodes must all be the same length. Place different barcode lengths in separate runs.\n";
            help();
        }
    }

    //
    // Determine the barcode length for the paired-end barcodes
    //
    prev = barcodes[0].pe.length();
    for (uint i = 0; i < barcodes.size(); i++) {
        pe_len = barcodes[i].pe.length();

        if (prev != pe_len) {
            cerr << "Paired-end barcodes must all be the same length. Place different barcode lengths in separate runs.\n";
            help();
        }
    }

    //
    // If paired barcodes were supplied chech that a paired barcode type was 
    // specified and vice versa.
    //
    if (se_bc.size() > 0 && pe_bc.size() > 0) {
	if (barcode_type != inline_inline &&
	    barcode_type != index_index &&
	    barcode_type != inline_index &&
	    barcode_type != index_inline) {
	    cerr << "You provided paried barcodes but did not specify a paired barcode type.\n";
	    help();
	}
    } else {
	if (barcode_type != inline_null &&
	    barcode_type != index_null) {
	    cerr << "You provided single-end barcodes but did not specify a single-end barcode type.\n";
	    help();
	}
    }

    cerr << "Loaded " << barcodes.size() << " barcodes ";

    if (pe_bc.size() > 0) 
	cerr << "(" << se_len << "bp / " << pe_len << "bp).\n";
    else
	cerr << "(" << se_len << "bp).\n";

    return 0;
}

int 
build_file_list(vector<pair<string, string> > &files) 
{

    //
    // Scan a directory for a list of files.
    //
    if (in_path_1.length() > 0) {
	string file, paired_file;
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
	    // 
	    // If paired-end specified, parse file names to sort out which is which.
	    //
	    if (paired) {
		int res;

		if ((res = parse_illumina_v1(file.c_str())) > 0 ||
		    (res = parse_illumina_v2(file.c_str())) > 0) {
		    paired_file = file;
		    paired_file.replace(res, 1, "2");
		    files.push_back(make_pair(file, paired_file));
		}
	    } else {
		files.push_back(make_pair(file, ""));
	    }
	}

	if (files.size() == 0) {
	    cerr << "Unable to locate any input files to process within '" << in_path_1 << "'\n";
	}
    } else {
	//
	// Files specified directly:
	//   Break off file path and store path and file name.
	//
	if (paired) {
	    int pos_1   = in_file_p1.find_last_of("/");
	    in_path_1 = in_file_p1.substr(0, pos_1 + 1);
	    int pos_2   = in_file_p2.find_last_of("/");
	    in_path_2 = in_file_p2.substr(0, pos_2 + 1);
	    files.push_back(make_pair(in_file_p1.substr(pos_1+1), in_file_p2.substr(pos_2+1)));
	} else {
	    int pos   = in_file.find_last_of("/");
	    in_path_1 = in_file.substr(0, pos + 1);
	    files.push_back(make_pair(in_file.substr(pos+1), ""));
	}
    }

    cerr << "Found " << files.size(); 
    if (paired)
	cerr << " paired input file(s).\n";
    else 
	cerr << " input file(s).\n";

    return 0;
}

