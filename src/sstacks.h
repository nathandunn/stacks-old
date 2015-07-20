// -*-mode:c++; c-style:k&r; c-basic-offset:4;-*-
//
// Copyright 2010-2013, Julian Catchen <jcatchen@uoregon.edu>
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

#ifndef __SSTACKS_H__
#define __SSTACKS_H__

#include "constants.h"

#ifdef _OPENMP
#include <omp.h>    // OpenMP library
#endif
#include <getopt.h> // Process command-line options
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <utility>
using std::pair;
using std::make_pair;

#include <string>
using std::string;

#include <iostream>
#include <fstream>
using std::ifstream;
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
#include <queue>
using std::queue;

#include <unordered_map>
using std::unordered_map;

#ifdef HAVE_SPARSEHASH
#include <sparsehash/sparse_hash_map>
using google::sparse_hash_map;
#endif

#include "kmers.h"
#include "stacks.h"
#include "locus.h"
#include "sql_utilities.h"
#include "utils.h"

#ifdef HAVE_SPARSEHASH
typedef sparse_hash_map<const char *, vector<pair<int, allele_type> >, hash_charptr, eqstr> HashMap;
#else
typedef unordered_map<const char *, vector<pair<int, allele_type> >, hash_charptr, eqstr> HashMap;
#endif

void help( void );
void version( void );
int  parse_command_line(int, char**);
int  populate_hash(map<int, Locus *> &, HashMap &, int);
int  find_matches_by_sequence(map<int, Locus *> &, map<int, QLocus *> &);
int  find_matches_by_genomic_loc(map<int, Locus *> &, map<int, QLocus *> &);
int  verify_sequence_match(map<int, Locus *> &, QLocus *, set<int> &, map<string, vector<string> > &, uint, unsigned long &, unsigned long &);
int  verify_genomic_loc_match(Locus *, QLocus *, set<string> &, unsigned long &);
int  generate_query_haplotypes(Locus *, QLocus *, set<string> &);
int  impute_haplotype(string, vector<pair<allele_type, string> > &, string &);
bool compare_dist(pair<int, int>, pair<int, int>);
int  write_matches(string, map<int, QLocus *> &);

#endif // __SSTACKS_H__
