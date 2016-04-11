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

#include "constants.h"

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
#include <unordered_map>
using std::unordered_map;

#ifdef HAVE_SPARSEHASH
#include <sparsehash/sparse_hash_map>
using google::sparse_hash_map;
#endif

#include "stacks.h"
#include "locus.h"
#include "mstack.h"
#include "input.h"

struct hash_charptr {
    size_t operator()(const char *__s) const
    {
        size_t __result = static_cast<size_t>(14695981039346656037ULL);
        unsigned int __len = strlen(__s);
        for (unsigned int i = 0; i < __len; i++) {
            __result ^= static_cast<size_t>(__s[i]);
            __result *= static_cast<size_t>(1099511628211ULL);
        }

        return __result;
    }
};
struct eqstr {
    bool operator()(const char* s1, const char* s2) const {
        return strcmp(s1, s2) == 0;
    }
};

#ifdef HAVE_SPARSEHASH
typedef sparse_hash_map<const char *, vector<int>, hash_charptr, eqstr> KmerHashMap;
typedef sparse_hash_map<const char *, vector<pair<string, int> >, hash_charptr, eqstr> CatKmerHashMap;
#else
typedef unordered_map<const char *, vector<int>, hash_charptr, eqstr> KmerHashMap;
typedef unordered_map<const char *, vector<pair<string, int> >, hash_charptr, eqstr> CatKmerHashMap;
#endif

int  determine_kmer_length(int, int);
int  calc_min_kmer_matches(int, int, int, bool);
int  initialize_kmers(int, int, vector<char *> &);
int  generate_kmers(const char *, int, int, vector<char *> &);
int  generate_kmers_lazily(const char *, uint, uint, vector<char *> &);

int  populate_kmer_hash(map<int, MergedStack *> &, KmerHashMap &, vector<char *> &, int);
int  populate_kmer_hash(map<int, Locus *> &, CatKmerHashMap &, vector<char *> &, int);
int  populate_kmer_hash(map<int, Locus *> &, KmerHashMap &, vector<char *> &, map<int, pair<allele_type, int> > &, int);
int  populate_kmer_hash(map<int, CLocus *> &, KmerHashMap &, vector<char *> &, map<int, pair<allele_type, int> > &, int);

int  free_kmer_hash(KmerHashMap &, vector<char *> &);
int  free_kmer_hash(CatKmerHashMap &, vector<char *> &);

int  generate_permutations(map<int, char **> &, int);

//
// Utilities
//
int dist(const char *, const char *, vector<pair<char, uint> > &);
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
