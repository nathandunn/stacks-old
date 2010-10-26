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

#ifndef __STACKS_H__
#define __STACKS_H__

#include <string>
using std::string;
#include <vector>
using std::vector;
#include <map>
using std::map;
#include <set>
using std::set;
#include <utility>
using std::pair;
using std::make_pair;
#include<iostream>
using std::cerr;
#include<sstream>
using std::stringstream;

#include "constants.h"
#include "input.h"

typedef unsigned int uint;
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
    uint   count;        // Number of identical reads forming this stack
    char  *seq;          // Sequence read
    vector<SeqId *> map; // List of sequence read IDs merged into this stack
    PhyLoc loc;          // Physical genome location of this stack.

    Stack()  { id = 0; count = 0; seq = NULL; loc.bp = 0; loc.chr[0] = '\0'; }
    ~Stack() { delete [] seq; }
    int add_id(const char *);
    int add_seq(const char *);
};

class Locus {
 public:
    int         id; // Locus ID
    int  sample_id; // Sample ID
    int      depth; // Stack depth
    char      *con; // Consensus sequence

    PhyLoc               loc;   // Physical genome location of this Stack.
    vector<SNP *>       snps;   // Single Nucleotide Polymorphisms in this stack
    map<string, int> alleles;   // Map of the allelic configuration of SNPs in this stack along with the count of each
    vector<pair<allele_type, string> > strings; // Strings for matching (representing the various allele combinations)

    Locus()  { id = 0; sample_id = 0; depth = 0; con = NULL; loc.bp = 0; loc.chr[0] = '\0'; }
    virtual ~Locus() { 
        delete [] con; 
        for (uint i = 0; i < snps.size(); i++)
            delete snps[i];
    }
    int add_consensus(const char *);
    virtual int populate_alleles();
};

#endif // __STACKS_H__
