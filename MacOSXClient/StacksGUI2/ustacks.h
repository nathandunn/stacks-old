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

#ifndef __USTACKS_H__
#define __USTACKS_H__

#ifdef _OPENMP
#include <omp.h>    // OpenMP library
#endif

#include <getopt.h> // Process command-line options
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <algorithm>
#include <utility>
using std::pair;
using std::make_pair;

#include <string>
using std::string;

#include <iostream>
#include <fstream>
#include <sstream>
using std::ofstream;
using std::stringstream;
using std::cin;
using std::cout;
using std::cerr;
using std::endl;

#include <vector>
using std::vector;

#include <map>
using std::map;

#include <queue>
using std::queue;

#include <set>
using std::set;

#ifdef __GNUC__
#include <ext/hash_map>
using __gnu_cxx::hash_map;
using __gnu_cxx::hash;
#else
#include <hash_map>
#endif

#include <unistd.h>


#include "constants.h" 
#include "kmers.h"
#include "utils.h"
#include "DNASeq.h"    // Class for storing two-bit compressed DNA sequences
#include "stacks.h"    // Major data structures for holding stacks
#include "mstack.h"
#include "mst.h"       // Minimum spanning tree implementation
#include "models.h"    // Contains maximum likelihood statistical models.
#include "Fasta.h"     // Reading input files in FASTA format
#include "Fastq.h"     // Reading input files in FASTQ format
#include "gzFasta.h"   // Reading gzipped input files in FASTA format
#include "gzFastq.h"   // Reading gzipped input files in FASTQ format

typedef unsigned int uint;

const int barcode_size = 5;

class HVal {
 public:
    vector<char *> ids;
    int  count;
    HVal() { this->count = 0; }
    ~HVal() { 
	for (uint i = 0; i < this->ids.size(); i++) 
	    delete [] this->ids[i];
	this->ids.clear(); 
    }

    int add_id(const char *id) {
	char *f = new char[strlen(id) + 1];
    	strcpy(f, id);
    	this->ids.push_back(f);
	return 0;
    }
};

typedef hash_map<DNASeq *, HVal, hash<DNASeq *>, dnaseq_eqstr> DNASeqHashMap;
typedef hash_map<const char *, HVal, hash<const char *>, eqstr> HashMap;

void help( void );
void version( void );
//int  parse_command_line(int, char**);
int  load_radtags(string, DNASeqHashMap &, vector<DNASeq *> &);
int  reduce_radtags(DNASeqHashMap &, map<int, Stack *> &, map<int, Rem *> &);
int  free_radtags_hash(DNASeqHashMap &, vector<DNASeq *> &);
int  populate_merged_tags(map<int, Stack *> &, map<int, MergedStack *> &);
int  merge_stacks(map<int, Stack *> &, map<int, Rem *> &, map<int, MergedStack *> &, set<int> &, int);
int  call_consensus(map<int, MergedStack *> &, map<int, Stack *> &, map<int, Rem *> &, bool);
int  call_alleles(MergedStack *, vector<DNASeq *> &, vector<read_type> &);
int  merge_remainders(map<int, MergedStack *> &, map<int, Rem *> &);
int  write_results(map<int, MergedStack *> &, map<int, Stack *> &, map<int, Rem *> &);

//
// Match MergedStacks using a k-mer hashing algorithm
//
int  calc_kmer_distance(map<int, MergedStack *> &, int);

//
// Calculate depth of coverage statistics for stacks
//
int    calc_coverage_distribution(map<int, Stack *> &, double &, double &);
double calc_merged_coverage_distribution(map<int, Stack *> &, map<int, MergedStack *> &);
int    count_raw_reads(map<int, Stack *> &, map<int, MergedStack *> &);

//
// Dealing with lumberjack (huge) stacks
//
int  calc_triggers(double, double, int &, int &);
int  remove_repetitive_stacks(map<int, Stack *> &, map<int, MergedStack *> &);
int  deleverage(map<int, Stack *> &, map<int, Rem *> &, map<int, MergedStack *> &, set<int> &, int, vector<MergedStack *> &);

//
// Debugging
//
int  dump_unique_tags(map<int, Stack *> &);
int  dump_merged_tags(map<int, MergedStack *> &);
int  dump_stack_graph(string, map<int, Stack *> &, map<int, MergedStack *> &, vector<int> &, map<int, map<int, double> > &, map<int, set<int> > &);

//
// Utilities
//
MergedStack *merge_tags(map<int, MergedStack *> &, set<int> &, int);
MergedStack *merge_tags(map<int, MergedStack *> &, int *, int, int);
long double factorial(int);

//
// Deprecated
//
int  calc_distance(map<int, MergedStack *> &, int);

#endif // __USTACKS_H__
