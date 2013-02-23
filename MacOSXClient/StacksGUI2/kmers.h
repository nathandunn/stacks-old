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

#ifndef __KMERS_H__
#define __KMERS_H__

#include <math.h>
#include <string.h>
#include <string>
using std::string;
#include <vector>
using std::vector;
#include <map>
using std::map;
#include <utility>
using std::pair;
using std::make_pair;
#include <iostream>
using std::ifstream;
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

#include "stacks.h"
#include "mstack.h"
#include "input.h"

struct eqstr {
    bool operator()(const char* s1, const char* s2) const {
	return strcmp(s1, s2) == 0;
    }
};

typedef hash_map<const char *, vector<int>, hash<const char *>, eqstr> KmerHashMap;
typedef hash_map<const char *, vector<pair<string, int> >, hash<const char *>, eqstr> CatKmerHashMap;

int  determine_kmer_length(int, int);
int  calc_min_kmer_matches(int, int, int, bool);
int  generate_kmers(const char *, int, int, vector<char *> &);
int  populate_kmer_hash(map<int, MergedStack *> &, KmerHashMap &, vector<char *> &, int);
int  populate_kmer_hash(map<int, Locus *> &, CatKmerHashMap &, vector<char *> &, int);
int  free_kmer_hash(KmerHashMap &, vector<char *> &);
int  free_kmer_hash(CatKmerHashMap &, vector<char *> &);

int  generate_permutations(map<int, char **> &, int);

//
// Utilities
//
int dist(const char *, Locus *, allele_type);
int dist(Locus *, Locus *);
int dist(MergedStack *, MergedStack *);
int dist(MergedStack *, char *);

//
// For sorting functions.
//
bool compare_dist(pair<int, int>, pair<int, int>);

//
// Debugging
//
int  dump_kmer_map(KmerHashMap &);

#endif // __KMERS_H__
