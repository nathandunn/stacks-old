// -*-mode:c++; c-style:k&r; c-basic-offset:4;-*-

#ifndef __HSTACKS_H__
#define __HSTACKS_H__

#include <omp.h>    // OpenMP library
#include <getopt.h> // Process command-line options
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
#include <ext/hash_fun.h>
using __gnu_cxx::hash_map;
using __gnu_cxx::hash;
#else
#include <hash_map>
#endif

#include "constants.h"
#include "stacks.h"
#include "sql_utilities.h"

typedef struct match {
    uint id;
    uint dist;
} Match;

//
// Homologous Locus Class
//
class HLocus : public Locus {
 public:
    int             uniq_id;   // An ID that is unique among all samples in this analysis.
    vector<Match *> matches;   // Matching tags found for the catalog.

    HLocus(): Locus() {}
    ~HLocus();

    int add_match(int, int);
};

HLocus::~HLocus() {
    vector<Match *>::iterator it;

    for (it = this->matches.begin(); it != this->matches.end(); it++)
        delete *it;
}

int HLocus::add_match(int id, int distance) {
    Match *m = new Match;

    m->id   = id;
    m->dist = distance;

    this->matches.push_back(m);

    return 0;
}

class MergedTag {
 public:
    int   id;     // Identifier for the merged tag. 
    char *con;    // Consensus sequence
    int   count;  // Number of merged unique tags (UTags)
    vector<int>      utags;       // Other UTags that have been merged into this unique tag
    vector<pair<int, int> > dist; // Vector describing the distance between this radtag and other radtags.
    vector<int>      remtags;     // Remainder tags that have been merged into this unique tag
    vector<SNP *>    snps;        // Single Nucleotide Polymorphisms found in this Radtag
    vector<string>   alleles;
    bool             deleveraged;
    string           sample_id;

    MergedTag()  { id = 0; count = 0; con = NULL; deleveraged = false; }
    ~MergedTag() { delete [] con; }
    int add_consensus(const char *);
    int add_dist(const int id, const int dist);
};

int MergedTag::add_consensus(const char *seq) {
    if (this->con != NULL)
	delete [] this->con;

    this->con = new char[strlen(seq) + 1];
    strcpy(this->con, seq);

    return 0;
}

int MergedTag::add_dist(const int id, const int dist) {
    //
    // Store the ID and distance as a pair, ID in the first position,
    // dist in the second.
    //
    pair<int, int> p(id, dist);
    this->dist.push_back(p);

    return 0;
}


void help( void );
void version( void );
int  parse_command_line(int, char**);
int  build_file_list(string, vector<string> &);
int  calc_distance(map<int, HLocus *> &, int);
int  dist(HLocus *, HLocus *);
bool compare_dist(Match *, Match *);
int  write_homologous_loci(map<int, HLocus *> &);
int  trace_stack_graph(HLocus *, map<int, HLocus *> &, set<int> &);

#endif // __HSTACKS_H__
