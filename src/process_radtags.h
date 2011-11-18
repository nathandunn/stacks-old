// -*-mode:c++; c-style:k&r; c-basic-offset:4;-*-
//
// Copyright 2010, Julian Catchen <jcatchen@uoregon.edu>
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

#ifndef __PROCESS_RADTAGS_H__
#define __PROCESS_RADTAGS_H__

#include <stdlib.h>
#include <getopt.h> // Process command-line options
#include <dirent.h> // Open/Read contents of a directory
#include <string.h>
#include <math.h>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <sstream>
using std::stringstream;
using std::istream;
using std::ofstream;
using std::cin;
using std::cout;
using std::cerr;
using std::endl;

#include <vector>
using std::vector;
#include <map>
using std::map;
#include <set>
using std::set;
#include <utility>
using std::pair;

#include "constants.h" 
#include "renz.h"
#include "clean.h"
#include "Bustard.h"   // Reading input files in Tab-separated Bustard format
#include "Fastq.h"     // Reading input files in FASTQ format

struct eqstr {
    bool operator()(const char* s1, const char* s2) const {
	return strcmp(s1, s2) == 0;
    }
};

void help( void );
void version( void );
int  parse_command_line(int, char**);
int  build_file_list(vector<pair<string, string> > &);
int  load_barcodes(vector<string> &);
int  open_files(vector<string> &, 
		map<string, ofstream *> &, 
		map<string, ofstream *> &, 
		map<string, ofstream *> &,
		map<string, map<string, long> > &);
int  close_file_handles(map<string, ofstream *> &);
int  process_reads(string, 
		   map<string, ofstream *> &, 
		   map<string, long> &, map<string, map<string, long> > &);
int  process_paired_reads(string, string, 
			  map<string, ofstream *> &, 
			  map<string, ofstream *> &, 
			  map<string, ofstream *> &,
			  map<string, long> &, map<string, map<string, long> > &);
int  process_singlet(map<string, ofstream *> &, Read *, map<string, map<string, long> > &, map<string, long> &, bool);
int  correct_barcode(map<string, ofstream *> &, Read *, map<string, long> &, map<string, map<string, long> > &);
int  correct_radtag(Read *, map<string, long> &);
int  check_quality_scores(Read *, bool);
int  dist(const char *, char *);
int  print_results(vector<string> &, map<string, map<string, long> > &, map<string, map<string, long> > &);

int  compare_barcodes(pair<string, int>, pair<string, int>);

#endif // __PROCESS_RADTAGS_H__
