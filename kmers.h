// -*-mode:c++; c-style:k&r; c-basic-offset:4;-*-

#ifndef __KMERS_H__
#define __KMERS_H__

#include <math.h>
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
#include <ext/hash_fun.h>
using __gnu_cxx::hash_map;
using __gnu_cxx::hash;
#else
#include <hash_map>
#endif

#include "stacks.h"
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
int  populate_kmer_hash(map<int, MergedStack *> &, KmerHashMap &, int);
int  populate_kmer_hash(map<int, Locus *> &, CatKmerHashMap &, int);

int  generate_permutations(map<int, char **> &, int);

//
// Utilities
//
int dist(const char *, Locus *, allele_type);
int dist(Locus *, Locus *);
int dist(MergedStack *, MergedStack *);
int dist(MergedStack *, Seq *);

//
// For sorting functions.
//
bool compare_dist(pair<int, int>, pair<int, int>);

//
// Debugging
//
int  dump_kmer_map(KmerHashMap &);

#endif // __KMERS_H__
