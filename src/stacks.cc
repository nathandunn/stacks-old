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

Rem::Rem() { 
    this->id         = 0;
    this->seq_id     = NULL; 
    this->seq        = NULL; 
    this->utilized   = false;
    this->loc.bp     = 0; 
    this->loc.chr[0] = '\0'; 
}

Rem::Rem(int id, char *seq_id, DNASeq *seq) { 
    this->id = id;

    uint len = strlen(seq_id);
    this->seq_id = new char[len + 1];
    strcpy(this->seq_id, seq_id);

    this->seq = new DNASeq(seq->size, seq->s);

    this->utilized   = false;
    this->loc.bp     = 0; 
    this->loc.chr[0] = '\0'; 
}

int Rem::add_id(const char *id) {
    if (this->seq_id != NULL)
	delete [] this->seq_id;

    uint len = strlen(id);
    this->seq_id = new char[len + 1];
    strcpy(this->seq_id, id);

    return 0;
}

int Rem::add_seq(const char *seq) {
    if (this->seq != NULL)
	delete [] this->seq;

    this->seq = new DNASeq(strlen(seq), seq);

    return 0;
}

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

    this->len = strlen(seq);
    this->seq = new DNASeq(this->len, seq);

    return 0;
}

int Stack::add_seq(const DNASeq *seq) {
    if (this->seq != NULL)
	delete [] this->seq;

    this->seq = new DNASeq(seq->size, seq->s);

    return 0;
}

int MergedStack::add_consensus(const char *seq) {
    if (this->con != NULL)
	delete [] this->con;

    this->len = strlen(seq);
    this->con = new char[this->len + 1];
    strncpy(this->con, seq, this->len);
    this->con[this->len] = '\0';

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

int Locus::add_consensus(const char *seq) {
    if (this->con != NULL)
	delete [] this->con;

    this->len = strlen(seq);
    this->con = new char[this->len + 1];
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

bool bp_compare(Locus *a, Locus *b) {
    return (a->loc.bp < b->loc.bp);
}
