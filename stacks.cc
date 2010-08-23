// -*-mode:c++; c-style:k&r; c-basic-offset:4;-*-
//
// stacks.cc -- routines for the stack-holding containers
//
// Julian Catchen
// jcatchen@uoregon.edu
// University of Oregon
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
