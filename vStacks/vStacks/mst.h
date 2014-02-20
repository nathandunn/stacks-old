// -*-mode:c++; c-style:k&r; c-basic-offset:4;-*-
//
// Copyright (c) 2014 University of Oregon
// Created by Julian Catchen <jcatchen@uoregon.edu>
//

#ifndef __MST_H__
#define __MST_H__

#include <stdio.h>
#include <string>
using std::string;
#include <map>
using std::map;
#include <vector>
using std::vector;
#include <set>
using std::set;
#include <queue>
using std::queue;
#include <sstream>
using std::stringstream;
#include<iostream>
using std::cerr;

#include "constants.h"

typedef unsigned int uint;

class Node;

class Edge {
 public:
    uint  dist;   // Distance or weight
    Node *child;
};

class Node {
public:
    uint id;
    vector<Edge *> edges;

    Node *parent;
    bool  update;
    uint  min_dist;

    //
    // List of adjacent nodes that are connected by minimal distance
    //
    vector<Node *> min_adj_list;

    Node(uint id) {
        this->id       = id; 
        this->parent   = NULL; 
        this->update   = true;
        this->min_dist = 1000000;
    }

    ~Node()  {
        for (uint i = 0; i < this->edges.size(); i++)
            delete this->edges[i];
    }

    Edge *add_edge(Node *, int);
};

class MinSpanTree {
    map<int, Node *> nodes;

 public:
    MinSpanTree()  { }
    ~MinSpanTree() {
        for (uint i = 0; i < this->nodes.size(); i++)
            delete this->nodes[i];
    }

    Node  *add_node(int id);
    int    build_tree();
    int    node_count();
    Node  *node(int id);
    Node  *head();
    bool   connected(int *, int);
    string vis(bool);
};

bool min_span_tree_cmp(const Node *, const Node *);

#endif // __MST_H__
