// -*-mode:c++; c-style:k&r; c-basic-offset:4;-*-
//
// Copyright (c) 2014 University of Oregon
// Created by Julian Catchen <jcatchen@uoregon.edu>
//

//
// locus.cc -- routines for the Locus class and its derivatives.
//
// Julian Catchen
// jcatchen@uoregon.edu
// University of Oregon
//
#include "locus.h"

uint Locus::sort_bp(uint k) {
    if (this->loc.strand == plus)
	return this->loc.bp + k;
    else
	return k == 0 ? this->loc.bp - this->len : this->loc.bp - k;
}

int Locus::add_consensus(const char *seq) {
    if (this->con != NULL)
	delete [] this->con;

    this->len = strlen(seq);
    this->con = new char[this->len + 1];
    strcpy(this->con, seq);

    return 0;
}

int 
Locus::populate_alleles() 
{
    vector<SNP *>::iterator  i;
    map<string, int>::iterator j;
    string s;
    int    k;

    //
    // Is this effective?
    //
    for (uint n = 0; n < this->strings.size(); n++) {
	this->strings[n].first.clear();
	this->strings[n].second.clear();
    }
    this->strings.clear();

    if (this->snps.size() == 0) {
	this->strings.push_back(make_pair("consensus", this->con));
	return 0;
    }

    for (j = this->alleles.begin(); j != this->alleles.end(); j++) {
	s = this->con;
	k = 0;

	for (i = this->snps.begin(); i != this->snps.end(); i++) {
	    if ((*i)->col < this->len)
		s.replace((*i)->col, 1, 1, j->first[k]);
	    k++;
	}

	this->strings.push_back(make_pair(j->first, s));
    }

    return 0;
}

bool 
bp_compare(Locus *a, Locus *b) 
{
    return (a->sort_bp() < b->sort_bp());
}

QLocus::~QLocus() 
{
    vector<Match *>::iterator it;

    for (it = this->matches.begin(); it != this->matches.end(); it++)
        delete *it;
}

int 
QLocus::add_match(int catalog_id, allele_type cat_type, allele_type query_type, int distance) 
{
    Match *m = new Match;

    m->cat_id     = catalog_id;
    m->cat_type   = cat_type;
    m->query_type = query_type;
    m->dist       = distance;

    this->matches.push_back(m);

    return 0;
}

int 
QLocus::add_match(int catalog_id, allele_type cat_type) 
{
    Match *m = new Match;

    m->cat_id     = catalog_id;
    m->cat_type   = cat_type;
    m->query_type = "";
    m->dist       = 0;

    this->matches.push_back(m);

    return 0;
}
