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

#ifndef __MST_H__
#define __MST_H__

#include <vector>
using std::vector;
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
        for (uint i = 0; i < edges.size(); i++)
            delete this->edges[i];
    }
};

bool min_span_tree_cmp(const Node *lhs, const Node *rhs) {
    return (lhs->min_dist > rhs->min_dist);
}

class MinSpanTree {
    vector<Node *> nodes;
    Node *head;

 public:
    MinSpanTree()  { head = NULL; }
    ~MinSpanTree() {
        for (uint i = 0; i < this->nodes.size(); i++)
            delete this->nodes[i];
    }

    int   add_node(Node *);
    int   build_tree();
    Node *head_node();
};

int MinSpanTree::add_node(Node *n) {
    this->nodes.push_back(n);

    return 0;
}

Node *MinSpanTree::head_node() {
    return this->head;
}

//
// Build a minimum spanning tree using Prim's alogorithm. Assume all necessary
// nodes have been added using the add_node function.
//
int MinSpanTree::build_tree() {
    //
    // Vector, which we treat as a binary heap to access nodes that are of minimal distance
    //
    vector<Node *> q;

    //
    // Select an initial node to process and initialize its minimum distance.
    //
    Node *n = this->nodes.front();
    n->min_dist = 0;

    this->head = n;

    //
    // Add all of the nodes to the binary heap; process them in order of min_dist
    //
    for (uint i = 0; i < this->nodes.size(); i++)
        q.push_back(this->nodes[i]);
    make_heap(q.begin(), q.end(), min_span_tree_cmp);

    while (q.size() > 0) {
        n = q.front();
        pop_heap(q.begin(), q.end());
        q.pop_back();
        n->update = false;

        //cerr << "Examining node: " << n->id << " (" << n->min_dist << ")\n";

        //
        // Record the minimum connection between parent and n.
        //
        if (n->parent != NULL) {
            n->parent->min_adj_list.push_back(n);
            //n->min_adj_list.push_back(n->parent);
        }

        //
        // Iterate through all of the edges of n and update the 
        // minimum distance to the proper nodes.
        //
        Edge *e;
        for (uint i = 0; i < n->edges.size(); i++) {
            e = n->edges[i];

            if (e->child->update == true && e->dist < e->child->min_dist) {
                e->child->parent   = n;
                e->child->min_dist = e->dist;

                //cerr << "  Updating node: " << e->child->id << " to have distance: " << e->child->min_dist << "\n";
            }
        }

        //
        // Resort the heap after possibly changing many min_dist values
        //
        make_heap(q.begin(), q.end(), min_span_tree_cmp);
    }

    return 0;
}

#endif // __MST_H__
