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

//
// stacks.cc -- routines for the stack-holding containers
//
// Julian Catchen
// jcatchen@uoregon.edu
// University of Oregon
//
// $Id$
//
#include "stacks.h"

int Stack::add_id(const char *id) {
    SeqId *f = new SeqId;
    strncpy(f->id, id, id_len - 1);
    f->id[id_len - 1] = '\0';
    this->map.push_back(f);

    return 0;
}

int Stack::add_seq(const char *seq) {
    if (this->seq != NULL)
	delete [] this->seq;

    size_t len = strlen(seq);
    this->seq = new char[len + 1];
    strncpy(this->seq, seq, len);
    this->seq[len] = '\0';

    return 0;
}

int MergedStack::add_consensus(const char *seq) {
    if (this->con != NULL)
	delete [] this->con;

    size_t len = strlen(seq);
    this->con = new char[len + 1];
    strncpy(this->con, seq, len);
    this->con[len] = '\0';

    return 0;
}

int MergedStack::add_dist(const int id, const int dist) {
    //
    // Store the ID and distance as a pair, ID in the first position,
    // dist in the second.
    //
    pair<int, int> p(id, dist);
    this->dist.push_back(p);

    return 0;
}

char **MergedStack::gen_matrix(map<int, Stack *> &unique, map<int, Seq *> &rem) {
    Stack *tag;

    //
    // Create a two-dimensional array, each row containing one read. For
    // each unique tag that has been merged together, add the sequence for
    // that tag into our array as many times as it originally occurred. 
    //
    // We do not allocate memory for the second dimension of the array, we simply
    // reuse the existing char arrays in the unique and rem maps
    //
    uint cnt = tag->utags.size() + tag->remtags.size();
    if (this->matrix != NULL)
        delete [] this->matrix;
    this->matrix = int * [cnt];

    vector<int>::iterator j;
    int i = 0;
    for (j = tag->utags.begin(); j != tag->utags.end(); j++) {
        tag = unique[*j];

        for (uint k = 0; k < tag->count; k++) {
            this->matrix[i] = utag->seq;
        }

        i++;
    }

    // For each remainder tag that has been merged into this Stack, add the sequence. 
    for (j = mtag->remtags.begin(); j != mtag->remtags.end(); j++) {
        this->matrix[i] = rem[*j]->seq;
        i++;
    }

    return this->matrix;
}

int MergedStack::calc_likelihood() {

    if (this->matrix == NULL)
        return 0;

    //
    // Iterate over each column of the array and call the consensus base.
    //
    int row, col;
    int length = strlen(this->matrix[0]);
    int height = this->utags.size() + this->remtags.size();
    string con;
    map<char, int> nuc;
    map<char, int>::iterator max, n;
    char *base;

    for (col = 0; col < length; col++) {
        nuc['A'] = 0; 
        nuc['C'] = 0;
        nuc['G'] = 0;
        nuc['T'] = 0;

        for (row = 0; row < height; row++) {
            base = this->matrix[row][col];
            //cerr << "    Row: " << row << " Col: " << col << " Base: " << *base << "\n";
            nuc[*base]++;
        }

        //
        // Find the base with a plurality of occurances and call it.
        //
        max = nuc.end();

        for (n = nuc.begin(); n != nuc.end(); n++) {

            if (max == nuc.end() || n->second > max->second)
                max = n;
        }
        con += max->first;

        // Search this column for the presence of a SNP
        call_multinomial_snp(this, col, nuc);
    }

    return 0;
}

int Locus::add_consensus(const char *seq) {
    if (this->con != NULL)
	delete [] this->con;

    this->con = new char[strlen(seq) + 1];
    strcpy(this->con, seq);

    return 0;
}

int Locus::populate_alleles() {
    vector<SNP *>::iterator  i;
    map<string, int>::iterator j;
    string s;
    int    k;

    this->strings.clear();

    if (this->snps.size() == 0) {
	this->strings.push_back(make_pair("consensus", this->con));
	return 0;
    }

    for (j = this->alleles.begin(); j != this->alleles.end(); j++) {
	s = this->con;
	k = 0;

	for (i = this->snps.begin(); i != this->snps.end(); i++) {
	    s.replace((*i)->col, 1, 1, j->first[k]);
	    k++;
	}

	this->strings.push_back(make_pair(j->first, s));
    }

    return 0;
}
