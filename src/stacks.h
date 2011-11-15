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

#include <string.h>
#include <string>
using std::string;
#include <vector>
using std::vector;
#include <map>
using std::map;
#include <set>
using std::set;
#include <algorithm>
#include <utility>
using std::pair;
using std::make_pair;
#include<iostream>
using std::cerr;
#include<sstream>
using std::stringstream;

#include "constants.h"
#include "DNASeq.h"

typedef unsigned int uint;
typedef string allele_type;

enum snp_type    {het, hom};
enum read_type   {primary, secondary};
enum strand_type {plus, minus};

class PhyLoc {
public:
    char       *chr;
    uint        bp;
    strand_type strand;

    void set(const char *chr, uint bp, strand_type strand) {
	this->chr    = new char[strlen(chr)  + 1];
	this->bp     = bp;
	this->strand = strand;
	strcpy(this->chr,  chr);
    }
    PhyLoc() {
	chr    = NULL;
	bp     = 0;
	strand = plus;
    }
    PhyLoc(const char *chr, uint bp) {
	this->chr    = new char[strlen(chr)  + 1];
	this->bp     = bp;
	this->strand = plus;
	strcpy(this->chr,  chr);
    }
    PhyLoc(const char *chr, uint bp, strand_type strnd) {
	this->chr    = new char[strlen(chr)  + 1];
	this->bp     = bp;
	this->strand = strnd;
	strcpy(this->chr,  chr);
    }
    ~PhyLoc() {
	delete [] chr;
    }
};

class SNP {
 public:
    snp_type type;   // Heterozygous or homozygous
    uint     col;
    float    lratio;
    char     rank_1;
    char     rank_2;
    char     rank_3;
    char     rank_4;

    SNP() {
	col    = 0;
	lratio = 0.0;
	rank_1 = 0;
	rank_2 = 0;
	rank_3 = 0;
	rank_4 = 0;
    }
};

class Stack {
 public:
    uint    id;
    uint    count;       // Number of identical reads forming this stack
    DNASeq *seq;         // Sequence read
    uint    len;         // Read length
    vector<char *> map;  // List of sequence read IDs merged into this stack
    PhyLoc  loc;         // Physical genome location of this stack.

    Stack()  { 
	id     = 0; 
	count  = 0; 
	seq    = NULL; 
	len    = 0; 
    }
    ~Stack() { 
	delete [] seq; 
	for (unsigned int i = 0; i < this->map.size(); i++) 
	    delete [] map[i]; 
    }
    int add_id(const char *);
    int add_seq(const char *);
    int add_seq(const DNASeq *);
};

class Rem {
 public:
    uint    id;
    char   *seq_id;
    DNASeq *seq;     // Sequence read
    PhyLoc  loc;     // Physical genome location of this stack.
    bool    utilized;

    Rem();
    Rem(int, char *, DNASeq *);
    ~Rem() { delete [] seq_id; delete seq; }
    int add_id(const char *);
    int add_seq(const char *);
};

class MergedStack {
 public:
    int   id;     // Identifier for the merged stack.
    char *con;    // Consensus sequence
    uint  len;    // Sequence length
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
        id    = 0;
        count = 0;
        con   = NULL;
        len   = 0;

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
    int add_consensus(DNASeq *);
    int add_dist(const int id, const int dist);
};

class Locus {
 public:
    int         id; // Locus ID
    int  sample_id; // Sample ID
    int      depth; // Stack depth
    char      *con; // Consensus sequence
    uint       len; // Sequence length

    vector<char *>      comp;   // Raw components in this stack.
    vector<char *>     reads;   // Sequence reads contributing to this stack.
    PhyLoc               loc;   // Physical genome location of this stack.
    vector<SNP *>       snps;   // Single Nucleotide Polymorphisms in this stack.
    map<string, int> alleles;   // Map of the allelic configuration of SNPs in this stack along with the count of each
    vector<pair<allele_type, string> > strings; // Strings for matching (representing the various allele combinations)

    Locus()  { 
	id        = 0; 
	sample_id = 0; 
	depth     = 0; 
	con       = NULL; 
	len       = 0;
    }
    virtual ~Locus() { 
        delete [] con; 
        for (uint i = 0; i < snps.size(); i++)
            delete snps[i];
        for (uint i = 0; i < comp.size(); i++)
            delete [] comp[i];
        for (uint i = 0; i < reads.size(); i++)
            delete [] reads[i];
    }
    int add_consensus(const char *);
    virtual int populate_alleles();
};

class CatMatch {
public:
    int   batch_id;
    int   cat_id;
    int   sample_id;
    int   tag_id;
    int   depth;
    char *haplotype;

    CatMatch() { 
	batch_id  = 0; 
	cat_id    = 0; 
	sample_id = 0; 
	tag_id    = 0; 
	depth     = 0; 
	haplotype = NULL; 
    }
    ~CatMatch() { 
	delete [] haplotype; 
    }
};

bool bp_compare(Locus *, Locus *);

#endif // __STACKS_H__
