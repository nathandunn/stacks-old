// -*-mode:c++; c-style:k&r; c-basic-offset:4;-*-
//
// Copyright 2012, Julian Catchen <jcatchen@uoregon.edu>
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

#ifndef __MSTACK_H__
#define __MSTACK_H__

#include <string>
using std::string;
#include <vector>
using std::vector;
#include <map>
using std::map;
#include <utility>
using std::pair;
using std::make_pair;
#include<iostream>
using std::cerr;

#include "stacks.h"

class MergedStack {
 public:
    int   id;     // Identifier for the merged stack.
    char *con;    // Consensus sequence
    uint  len;    // Sequence length

    //
    // Stack component parts
    //
    int                    count; // Number of merged stacks
    vector<int>            utags; // Stack IDs that have been merged into this MergedStack
    vector<pair<int, int> > dist; // Vector describing the distance between this stack and other stacks.
    vector<int>          remtags; // Remainder tag IDs that have been merged into this Stack
    DNASeq              **matrix; // Two-dimensional array for iterating over the combined stack (stacks and remainders).
    DNANSeq            **pmatrix; // Two-dimensional array for iterating over the combined stack aligned to a reference..

    int cohort_id; // Group ID of all stacks that were originally part of the same subgraph
    double    lnl; // Log likelihood of this stack

    //
    // Mapping components
    //
    PhyLoc               loc; // Physical genome location of this Stack.
    vector<SNP *>       snps; // Single Nucleotide Polymorphisms found in this Stack
    map<string, int> alleles; // Set of alleles defined by the SNPs found in this Stack
    //
    // Flags
    //
    bool deleveraged;
    bool masked;
    bool blacklisted;
    bool lumberjackstack;

    MergedStack();
    ~MergedStack();
    int       add_consensus(const char *);
    int       add_consensus(DNASeq *);
    int       add_consensus(DNANSeq *);
    int       add_dist(const int id, const int dist);
    DNASeq  **gen_matrix(map<int, Stack *> &, map<int, Rem *> &);
    DNANSeq **gen_matrix(map<int, PStack *> &);
    double    calc_likelihood();
    double    calc_likelihood_pstacks();
    string    write_cmb();
};

#endif // __MSTACK_H__
