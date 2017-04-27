// -*-mode:c++; c-style:k&r; c-basic-offset:4;-*-
//
// Copyright 2017, Julian Catchen <jcatchen@illinois.edu>
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

#ifndef __SUFFIXTREE_H__
#define __SUFFIXTREE_H__

#include <queue>
using std::queue;

#include "constants.h"
#include "utils.h"
#include "DNASeq4.h"

const size_t NT4cnt = 6;

class STEdge;

class STNode {
    size_t  id_;
    STEdge *edges_[NT4cnt]; // Indexed according to Nt4
    STNode *suffix_link_;

public:
    STNode(size_t id): id_(id), suffix_link_(NULL) { for (uint i = 0; i < NT4cnt; i++) this->edges_[i] = NULL; }
    ~STNode();

    size_t  id()          { return this->id_; }
    STEdge *edge(Nt4 n)   { return this->edges_[n.index()]; }
    STEdge *edge(uint i)  { return this->edges_[i]; }
    STNode *suffix_link() { return this->suffix_link_; }
    STEdge *add_edge(Nt4 n, int pos);
    STNode *add_suffix_link(STNode *);
};


class STEdge {
    int     s_;
    int     e_;
    STNode *succ_;

public:
    STEdge() : s_(-1), e_(-1), succ_(NULL) {}
    STEdge(int s) : s_(s), e_(-1), succ_(NULL) {}
    STEdge(int s, int e, STNode *succ) : s_(s), e_(e), succ_(succ) {}
    ~STEdge() { delete this->succ_; }

    int     start()         { return this->s_; }
    int     start(int s)    { this->s_ = s; return this->s_; }
    int     end()           { return this->e_; }
    int     end(int e)      { this->e_ = e; return this->e_; }
    STNode *succ()          { return this->succ_; }
    STNode *succ(STNode *s) { this->succ_ = s; return this->succ_; }
};


class SuffixTree {
    DNASeq4 seq_;
    STNode *root;

public:
    SuffixTree(DNASeq4 s): seq_(s)  { this->root = new STNode(1); }
    SuffixTree()  { this->root = new STNode(1); }
    ~SuffixTree() { delete this->root; }

    size_t build_tree();
    size_t write_dot(ofstream &);
    size_t write_suffixes(ostream &);
private:
    size_t write_suffix(vector<string> &, string, STNode *);
    int    forward_nodes(STNode **, Nt4 &, int &, int, int &, int);
};

#endif
