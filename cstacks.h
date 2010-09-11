// -*-mode:c++; c-style:k&r; c-basic-offset:4;-*-

#ifndef __CSTACKS_H__
#define __CSTACKS_H__

#define VERSION 0.90

#include <omp.h>    // OpenMP library
#include <getopt.h> // Process command-line options
#include <math.h>
#include <stdlib.h>
#include <utility>
using std::pair;
using std::make_pair;

#include <string>
using std::string;

#include <iostream>
#include <fstream>
#include <sstream>
using std::ifstream;
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

#include <set>
using std::set;

#include <queue>
using std::queue;

#include "constants.h"
#include "stacks.h"
#include "kmers.h"
#include "sql_utilities.h"

enum searcht {sequence, genomic_loc};

typedef struct match {
    uint        cat_id;
    allele_type cat_type;
    allele_type query_type;
    uint        dist;
} Match;

//
// Query Locus Class
//
class QLocus : public Locus {
 public:
    vector<Match *> matches;   // Matching tags found for the catalog. Stored as catalog ID/allele_type pair.

    QLocus(): Locus() {}
    ~QLocus();

    int merge_snp(SNP *snp);
    int add_match(int, allele_type, allele_type, int);
};

QLocus::~QLocus() {
    vector<Match *>::iterator it;

    for (it = this->matches.begin(); it != this->matches.end(); it++)
        delete *it;
}

int QLocus::add_match(int catalog_id, allele_type cat_type, allele_type query_type, int distance) {
    Match *m = new Match;

    m->cat_id     = catalog_id;
    m->cat_type   = cat_type;
    m->query_type = query_type;
    m->dist       = distance;

    this->matches.push_back(m);

    return 0;
}

//
// Catalog Locus Class
//
class CLocus : public Locus {
 public:
    vector<pair<int, int> > sources;   // Sample/ID pairs for the sources contributing to this catalog entry

    int merge_snps(QLocus *);
};

void help( void );
void version( void );
int  parse_command_line(int, char**);
int  parse_tsv(const char *, vector<string> &);
int  initialize_catalog(pair<int, string> &, map<int, CLocus *> &);
int  find_kmer_matches_by_sequence(map<int, CLocus *> &, map<int, QLocus *> &, int);
int  find_matches_by_sequence(map<int, CLocus *> &, map<int, QLocus *> &);
int  find_matches_by_genomic_loc(map<int, CLocus *> &, map<int, QLocus *> &);
int  characterize_mismatch_snps(CLocus *, QLocus *);
int  merge_matches(map<int, CLocus *> &, map<int, QLocus *> &, pair<int, string> &, int);
int  add_unique_tag(pair<int, string> &, map<int, CLocus *> &, QLocus *);
bool compare_dist(pair<int, int>, pair<int, int>);
int  write_catalog(map<int, CLocus *> &);
int  write_simple_output(CLocus *, ofstream &, ofstream &, ofstream &);
bool compare_pair(pair<string, SNP *>, pair<string, SNP *>);
bool compare_matches(Match *, Match *);

int  populate_kmer_hash(map<int, CLocus *> &, CatKmerHashMap &, int);

#endif // __CSTACKS_H__
