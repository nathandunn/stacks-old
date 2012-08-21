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

#ifndef __PROCESS_SHORTREADS_H__
#define __PROCESS_SHORTREADS_H__

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

#ifdef __GNUC__
#include <ext/hash_map>
using __gnu_cxx::hash_map;
using __gnu_cxx::hash;
#else
#include <hash_map>
#endif

#include "constants.h" 
#include "clean.h"
#include "kmers.h"
#include "file_io.h"
#include "Bustard.h"   // Reading input files in Tab-separated Bustard format
#include "Fastq.h"     // Reading input files in FASTQ format

void help( void );
void version( void );
int  parse_command_line(int, char **);
int  process_reads(string, map<string, ofstream *> &, 
		   map<string, long> &, map<string, map<string, long> > &);
int  process_paired_reads(string, string, 
			  map<string, ofstream *> &, map<string, ofstream *> &, map<string, ofstream *> &, 
			  map<string, long> &, map<string, map<string, long> > &);
int  process_singlet(map<string, ofstream *> &, Read *, map<string, map<string, long> > &, map<string, long> &, bool);
int  correct_barcode(map<string, ofstream *> &, Read *, map<string, long> &, map<string, map<string, long> > &);
int  check_quality_scores(Read *, bool);
int  dist(const char *, char *);
int  print_results(int, char **, vector<string> &, map<string, map<string, long> > &, map<string, map<string, long> > &);
int  compare_barcodes(pair<string, int>, pair<string, int>);

#endif // __PROCESS_SHORTREADS_H__
