// -*-mode:c++; c-style:k&r; c-basic-offset:4;-*-
//
// Copyright 2010-2016, Julian Catchen <jcatchen@illinois.edu>
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

#include "constants.h" 

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
#include <iomanip> // std::setprecision

#include <vector>
using std::vector;
#include <map>
using std::map;
#include <unordered_map>
using std::unordered_map;
#include <queue>
using std::queue;
#include <set>
using std::set;
#include <unistd.h>

#ifdef HAVE_SPARSEHASH
#include <sparsehash/sparse_hash_map>
using google::sparse_hash_map;
#endif

#include "kmers.h"
#include "utils.h"
#include "DNASeq.h"     // Class for storing two-bit compressed DNA sequences
#include "stacks.h"     // Major data structures for holding stacks
#include "mstack.h"
#include "mst.h"        // Minimum spanning tree implementation
#include "models.h"     // Contains maximum likelihood statistical models.
#include "FastaI.h"     // Reading input files in FASTA format
#include "FastqI.h"     // Reading input files in FASTQ format
#include "gzFasta.h"    // Reading gzipped input files in FASTA format
#include "gzFastq.h"    // Reading gzipped input files in FASTQ format

class HVal {
 public:
    vector<int> ids;

    int count() {
	return this->ids.size();
    }
    int add_id(int id) {
    	this->ids.push_back(id);
	return 0;
    }
};

enum dynprog {dynp_down, dynp_right, dynp_diag};

const int    barcode_size   = 5;
const double gapopen_score  = -10;
const double gapext_score   = -0.5;
const double mismatch_score = -4;
const double match_score    =  5;

class AlignPath {
public:
    bool diag;
    bool left;
    bool up;

    AlignPath() {
	diag = false;
	left = false;
	up   = false;
    }
    int count() {
	int cnt;
	cnt  = this->up   ? 1 : 0;
	cnt += this->left ? 1 : 0;
	cnt += this->diag ? 1 : 0;

	return cnt;
    }
};

#ifdef HAVE_SPARSEHASH
typedef sparse_hash_map<DNANSeq *, HVal, hash_dnanseq, dnanseq_eqstr> DNASeqHashMap;
#else
typedef unordered_map<DNANSeq *, HVal, hash_dnanseq, dnanseq_eqstr> DNASeqHashMap;
#endif

void help( void );
void version( void );
int  parse_command_line(int, char**);
int  load_radtags(string, DNASeqHashMap &, vector<DNANSeq *> &);
int  load_seq_ids(vector<char *> &);
int  reduce_radtags(DNASeqHashMap &, map<int, Stack *> &, map<int, Rem *> &);
int  free_radtags_hash(DNASeqHashMap &, vector<DNANSeq *> &);
int  populate_merged_tags(map<int, Stack *> &, map<int, MergedStack *> &);
int  merge_stacks(map<int, Stack *> &, map<int, Rem *> &, map<int, MergedStack *> &, set<int> &, int);
int  call_consensus(map<int, MergedStack *> &, map<int, Stack *> &, map<int, Rem *> &, bool);
int  call_alleles(MergedStack *, vector<DNANSeq *> &, vector<read_type> &);
int  merge_remainders(map<int, MergedStack *> &, map<int, Rem *> &);
int  write_results(map<int, MergedStack *> &, map<int, Stack *> &, map<int, Rem *> &);

//
// Match MergedStacks using a k-mer hashing algorithm
//
int  calc_kmer_distance(map<int, MergedStack *> &, int);
int  search_for_gaps(map<int, MergedStack *> &, double);
int  merge_gapped_alns(map<int, Stack *> &, map<int, Rem *> &, map<int, MergedStack *> &);
int  parse_cigar(const char *, vector<pair<char, uint> > &);
int  edit_gapped_seqs(map<int, Stack *> &, map<int, Rem *> &, MergedStack *, vector<pair<char, uint> > &);
int  edit_gaps(vector<pair<char, uint> > &, char *);

//
// Needleman-Wunsch Alignment
//
int init_alignment(int, double ***, AlignPath ***);
int free_alignment(int, double **, AlignPath **);
int align(MergedStack *, MergedStack *, double **, AlignPath **);
inline int swap(double *, dynprog *, int, int);
int trace_alignment(MergedStack *, MergedStack *, AlignPath **);
int dump_alignment(MergedStack *, MergedStack *, double **, AlignPath **);

//
// Calculate depth of coverage statistics for stacks
//
int    calc_coverage_distribution(map<int, Stack *> &, double &, double &);
double calc_merged_coverage_distribution(map<int, Stack *> &, map<int, MergedStack *> &);
int    count_raw_reads(map<int, Stack *> &, map<int, Rem *> &, map<int, MergedStack *> &);

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
MergedStack *merge_tags(MergedStack *, MergedStack *, int);
MergedStack *merge_tags(map<int, MergedStack *> &, set<int> &, int);
MergedStack *merge_tags(map<int, MergedStack *> &, int *, int, int);
long double factorial(int);

//
// Deprecated
//
int  calc_distance(map<int, MergedStack *> &, int);

#endif // __USTACKS_H__
