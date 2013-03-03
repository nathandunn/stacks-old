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

#ifndef __HSTACKS_H__
#define __HSTACKS_H__

#ifdef _OPENMP
#include <omp.h>    // OpenMP library
#endif
#include <getopt.h> // Process command-line options
#include <string.h>
#include <math.h>
#include <stdlib.h>
#include <dirent.h>
#include <errno.h>
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
#include <algorithm>

#include "constants.h"
#include "stacks.h"
#include "locus.h"
#include "kmers.h"
#include "models.h"
#include "sql_utilities.h"

//
// A  map holding k-mer permutation strings. For use when generating fuzzy k-mers.
//
map<int, char **> pstrings;

//
// Homologous Locus Class
//
class HLocus : public Locus {
 public:
    int             uniq_id;   // An ID that is unique among all samples in this analysis.
    vector<Match *> matches;   // Matching tags found for the catalog.

    HLocus(): Locus() {}
    ~HLocus();

    int populate_alleles();
    int add_match(int, int);
};

void help( void );
void version( void );
int  parse_command_line(int, char**);
int  build_file_list(string, vector<string> &);
int  calc_kmer_distance(map<int, HLocus *> &, int);
int  calc_distance(map<int, HLocus *> &, int);
int  dist(HLocus *, HLocus *);
int  write_homologous_loci(map<int, HLocus *> &);
int  trace_stack_graph(HLocus *, map<int, HLocus *> &, set<int> &);
int  call_consensus(map<int, HLocus *> &, set<int> &, string &, vector<SNP *> &, vector<string> &);
int  call_alleles(vector<char *> &, vector<SNP *> &, vector<string> &);

int  populate_kmer_hash(map<int, HLocus *> &, CatKmerHashMap &, int);

bool compare_mdist(Match *, Match *);
bool compare_pair(pair<char, int>, pair<char, int>);

#endif // __HSTACKS_H__
