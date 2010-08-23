// -*-mode:c++; c-style:k&r; c-basic-offset:4;-*-

#ifndef __STACKS_H__
#define __STACKS_H__

#include <string>
using std::string;
#include <vector>
using std::vector;
#include <map>
using std::map;
#include <utility>
using std::pair;
using std::make_pair;

#include "constants.h"

typedef string allele_type;

typedef struct seqid {
    char id[id_len];
} SeqId;

typedef struct locus {
    char chr[id_len];
    uint bp;
} PhyLoc;

class SNP {
 public:
    uint  col;
    float lratio;
    char  rank_1;
    char  rank_2;
};

class Stack {
 public:
    uint   id;
    uint   count;         // Number of identical reads forming this stack
    char  *seq;           // Sequence read
    vector<SeqId *> map; // List of sequence read IDs merged into this stack
    PhyLoc loc;           // Physical genome location of this stack.

    Stack()  { id = 0; count = 0; seq = NULL; }
    ~Stack() { delete [] seq; }
    int add_id(const char *);
    int add_seq(const char *);
};

class MergedStack {
 public:
    int   id;     // Identifier for the merged stack.
    char *con;    // Consensus sequence
    //
    // Stack component parts
    //
    int                     count;   // Number of merged stacks
    vector<int>             utags;   // Other stacks that have been merged into this MergedStack
    vector<pair<int, int> > dist;    // Vector describing the distance between this stack and other stacks.
    vector<int>             remtags; // Remainder tags that have been merged into this Stack
    //
    // Mapping components
    //
    PhyLoc           loc;     // Physical genome location of this Stack.
    vector<SNP *>    snps;    // Single Nucleotide Polymorphisms found in this Stack
    map<string, int> alleles; // Set of alleles defined by the SNPs found in this Stack
    //
    // K-mers generated from the consensus sequence
    //
    vector<char *> kmers;
    //
    // Flags
    //
    bool deleveraged;
    bool masked;
    bool blacklisted;
    bool lumberjackstack;

    MergedStack()  { 
        id          = 0;
        count       = 0;
        con         = NULL;
        deleveraged     = false;
        masked          = false;
        blacklisted     = false;
        lumberjackstack = false;
    }
    ~MergedStack() { 
        delete [] con;

        for (uint i = 0; i < kmers.size(); i++)
            delete [] kmers[i];
        for (uint i = 0; i < snps.size(); i++)
            delete snps[i];
    }
    int add_consensus(const char *);
    int add_dist(const int id, const int dist);
};

class Locus {
 public:
    int   id;   // Locus ID
    char *con;  // Consensus sequence

    PhyLoc               loc;   // Physical genome location of this Stack.
    vector<SNP *>       snps;   // Single Nucleotide Polymorphisms in this stack
    map<string, int> alleles;   // Map of the allelic configuration of SNPs in this stack along with the count of each
    vector<pair<allele_type, string> > strings; // Strings for matching (representing the various allele combinations)

    Locus()  { id = 0; con = NULL; }
    ~Locus() { delete [] con; }
    int add_consensus(const char *);
    int populate_alleles();
};

#endif // __STACKS_H__
