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

#ifndef __PSTACKS_H__
#define __PSTACKS_H__

#ifdef _OPENMP
#include <omp.h>    // OpenMP library
#endif
#include <getopt.h> // Process command-line options
#include <string.h>
#include <iostream>
#include <fstream>
#include <sstream>
using std::stringstream;
using std::cin;
using std::cout;
using std::cerr;
using std::endl;

#ifdef __GNUC__
#include <ext/hash_map>
using __gnu_cxx::hash_map;
using __gnu_cxx::hash;
#else
#include <hash_map>
#endif

#include <vector>
using std::vector;
#include <map>
using std::map;
#include <set>
using std::set;
#include <utility>
using std::pair;

#include "constants.h" 
#include "stacks.h"    // Major data structures for holding stacks
#include "mstack.h"
#include "utils.h"
#include "models.h"    // Contains maximum likelihood statistical models.
#include "Tsv.h"       // Reading input files in Tab-separated values format
#include "Bowtie.h"    // Reading input files in Bowtie format
#include "Sam.h"       // Reading input files in SAM format
#include "Bam.h"       // Reading input files in BAM format
#include "DNANSeq.h"

const int barcode_size = 5;

struct eqstr {
    bool operator()(const char* s1, const char* s2) const {
	return strcmp(s1, s2) == 0;
    }
};

typedef hash_map<const char *, vector<Seq *>, hash<const char *>, eqstr> HashMap;

void help( void );
void version( void );
int  parse_command_line(int, char**);
int  load_radtags(string, HashMap &);
int  reduce_radtags(HashMap &, map<int, PStack *> &);
int  populate_merged_tags(map<int, PStack *> &, map<int, MergedStack *> &);
int  call_consensus(map<int, MergedStack *> &, map<int, PStack *> &, bool);
int  call_alleles(MergedStack *, vector<DNANSeq *> &);
int  count_raw_reads(map<int, PStack *> &, map<int, MergedStack *> &);
int  write_sql(map<int, MergedStack *> &, map<int, PStack *> &);
int  write_sam(map<int, MergedStack *> &, map<int, PStack *> &);

//
// Debugging
//
int  dump_stacks(map<int, PStack *> &);
int  dump_merged_stacks(map<int, MergedStack *> &);


#endif // __PSTACKS_H__
