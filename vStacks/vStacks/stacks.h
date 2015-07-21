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
using std::ofstream;
using std::cerr;
#include<sstream>
using std::stringstream;

#include "constants.h"
#include "DNASeq.h"
#include "DNANSeq.h"

typedef unsigned int uint;
typedef string allele_type;

enum snp_type    {snp_type_het, snp_type_hom, snp_type_unk};
enum read_type   {primary, secondary};
enum strand_type {plus, minus};
enum searcht     {sequence, genomic_loc};

class PhyLoc {
public:
    char       *chr;
    uint        bp;
    strand_type strand;

    void set(const char *chr, uint bp, strand_type strand) {
	if (this->chr != NULL) 
	    delete [] this->chr;
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

class PStack {
 public:
    uint     id;
    uint     count;      // Number of identical reads forming this stack
    DNANSeq *seq;        // Sequence read
    uint     len;        // Read length
    vector<char *> map;  // List of sequence read IDs merged into this stack
    PhyLoc   loc;        // Physical genome location of this stack.

    PStack()  { 
	id     = 0; 
	count  = 0; 
	seq    = NULL; 
	len    = 0; 
    }
    ~PStack() { 
	delete this->seq; 
	for (unsigned int i = 0; i < this->map.size(); i++) 
	    delete [] this->map[i]; 
    }
    int  add_id(const char *);
    int  add_seq(const char *);
    int  add_seq(DNANSeq *);
};

class Stack {
 public:
    uint         id;
    DNASeq     *seq;  // Sequence read
    vector<uint> map; // List of sequence read IDs merged into this stack

    Stack()  {
	id  = 0;
	seq = NULL;
    }
    ~Stack() { 
	delete this->seq; 
    }
    uint count() { return this->map.size(); }
    int  add_id(uint);
    int  add_seq(const char *);
    int  add_seq(const DNASeq *);
};

class Rem {
 public:
    uint          id;
    vector<uint> map; // List of sequence read IDs merged into this stack
    DNASeq      *seq; // Sequence read
    bool    utilized;

    Rem();
    Rem(int, uint, DNASeq *);
    ~Rem() { 
	delete this->seq;
    }
    uint count() { return this->map.size(); }
    int  add_id(uint);
    int  add_seq(const char *);
    int  add_seq(const DNASeq *);
};

class CatMatch {
public:
    int    batch_id;
    int    cat_id;
    int    sample_id;
    int    tag_id;
    int    depth;
    double lnl;
    char  *haplotype;

    CatMatch() { 
	batch_id  = 0; 
	cat_id    = 0; 
	sample_id = 0; 
	tag_id    = 0; 
	depth     = 0; 
	lnl       = 0.0;
	haplotype = NULL; 
    }
    ~CatMatch() { 
	delete [] haplotype; 
    }
};

class ModRes {
public:
    int   sample_id;
    int   tag_id;
    char *model;

    ModRes(int samp_id, int tag_id, const char *model) { 
	this->sample_id = samp_id; 
	this->tag_id    = tag_id;
	this->model     = new char[strlen(model) + 1];
	strcpy(this->model, model);
    }
    ~ModRes() { 
	delete [] this->model; 
    }
};

class SNPRes {
public:
    int   sample_id;
    int   tag_id;
    vector<SNP *> snps;

    SNPRes(int samp_id, int tag_id) { 
	this->sample_id = samp_id; 
	this->tag_id    = tag_id;
    }
    ~SNPRes() { 
	for (uint i = 0; i < this->snps.size(); i++)
	    delete this->snps[i];
	this->snps.clear();
    }
};

#endif // __STACKS_H__
