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
// mstack.cc -- implementation of the MergedStack Class
//
// Julian Catchen
// jcatchen@uoregon.edu
// University of Oregon
//
// $Id: mstack.cc 1987 2010-11-01 05:43:50Z catchen $
//

#include "mstack.h"
#include "models.h"

MergedStack::MergedStack()  { 
    this->id         = 0;
    this->count      = 0;
    this->len        = 0;
    this->con        = NULL;
    this->matrix     = NULL;
    this->lnl        = 0.0;
    this->cohort_id  = -1;

    this->deleveraged     = false;
    this->masked          = false;
    this->blacklisted     = false;
    this->lumberjackstack = false;
}

MergedStack::~MergedStack() { 
    delete [] this->con;

    for (uint i = 0; i < snps.size(); i++)
        delete this->snps[i];

    delete [] this->matrix;
}

int MergedStack::add_consensus(const char *seq) {
    if (this->con != NULL)
	delete [] this->con;

    this->len = strlen(seq);
    this->con = new char[len + 1];
    strncpy(this->con, seq, len);
    this->con[len] = '\0';

    return 0;
}

int MergedStack::add_consensus(DNASeq *seq) {
    if (this->con != NULL)
	delete [] this->con;

    this->len = seq->size;
    this->con = new char[this->len + 1];
    this->con = seq->seq(this->con);

    return 0;
}

int MergedStack::add_consensus(DNANSeq *seq) {
    if (this->con != NULL)
	delete [] this->con;

    this->len = seq->size();
    this->con = new char[this->len + 1];
    this->con = seq->seq(this->con);

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

DNASeq **MergedStack::gen_matrix(map<int, Stack *> &unique, map<int, Rem *> &rem) {
    Stack *tag;

    //
    // Create a two-dimensional array, each row containing one read. For
    // each unique tag that has been merged together, add the sequence for
    // that tag into our array as many times as it originally occurred. 
    //
    // We do not allocate memory for the second dimension of the array, we simply
    // reuse the existing char arrays in the unique and rem maps
    //
    uint cnt = this->count + this->remtags.size();
    if (this->matrix != NULL)
        delete [] this->matrix;
    this->matrix = new DNASeq * [cnt];

    vector<int>::iterator j;
    int i = 0;
    for (j = this->utags.begin(); j != this->utags.end(); j++) {
        tag = unique[*j];

        for (uint k = 0; k < tag->count(); k++) {
            this->matrix[i] = tag->seq;
            i++;
        }
    }

    // For each remainder tag that has been merged into this Stack, add the sequence. 
    for (j = this->remtags.begin(); j != this->remtags.end(); j++) {
        this->matrix[i] = rem[*j]->seq;
        i++;
    }

    return this->matrix;
}

DNANSeq **
MergedStack::gen_matrix(map<int, PStack *> &unique) 
{
    PStack *tag;
    //
    // Create a two-dimensional array, each row containing one read. For
    // each unique tag that has been merged together, add the sequence for
    // that tag into our array as many times as it originally occurred. 
    //
    // We do not allocate memory for the second dimension of the array, we simply
    // reuse the existing char arrays in the unique and rem maps
    //
    uint cnt = this->count;
    if (this->pmatrix != NULL)
        delete [] this->matrix;
    this->pmatrix = new DNANSeq * [cnt];

    vector<int>::iterator j;
    int i = 0;
    for (j = this->utags.begin(); j != this->utags.end(); j++) {
        tag = unique[*j];

        for (uint k = 0; k < tag->count; k++) {
            this->pmatrix[i] = tag->seq;
            i++;
        }
    }

    return this->pmatrix;
}

double 
MergedStack::calc_likelihood() 
{
    if (this->matrix == NULL || this->snps.size() == 0)
        return 0;

    //
    // Iterate over each column of the array and call the consensus base.
    //
    int row, col, tot;
    int length = this->matrix[0]->size;
    int height = this->count + this->remtags.size();
    map<char, int> nuc;
    map<char, int>::iterator max, n;
    DNASeq *d;

    this->lnl = 0;

    for (col = 0; col < length; col++) {
        nuc['A'] = 0; 
        nuc['G'] = 0;
        nuc['C'] = 0;
        nuc['T'] = 0;

        //
        // Count the nucleotide type at each position in the column.
        //
        for (row = 0; row < height; row++) {
	    d = this->matrix[row];
            nuc[(*d)[col]]++;
	}
        //
        // Find the base with a plurality of occurances and call it.
        //
        max = nuc.end();
        tot = 0;
        for (n = nuc.begin(); n != nuc.end(); n++) {
            tot += n->second;
            if (max == nuc.end() || n->second > max->second)
                max = n;
        }

	//
	// For nucleotide positions with potential polymorphism (i.e. two or more alleles at 
	// the locus that differ at that position), first find the ML genotype (call_multinomial_snp). 
	// If it returns 'het' calculate the heterozygous_likelihood(), otherwise calculate homozygous
	// likelihood.
	//
	snp_type res = this->snps[col]->type;

	if (res == snp_type_het) 
	    this->lnl += heterozygous_likelihood(col, nuc);
	else if (res == snp_type_hom)
	    this->lnl += homozygous_likelihood(col, nuc);
	else {
	    double homlnl = homozygous_likelihood(col, nuc);
	    double hetlnl = heterozygous_likelihood(col, nuc);
	    this->lnl += hetlnl > homlnl ? hetlnl : homlnl;
	}
    }

    return this->lnl;
}

double 
MergedStack::calc_likelihood_pstacks() 
{
    if (this->pmatrix == NULL || this->snps.size() == 0)
        return 0;

    //
    // Iterate over each column of the array and call the consensus base.
    //
    int row, col, tot;
    int length = this->pmatrix[0]->size();
    int height = this->count;
    map<char, int> nuc;
    map<char, int>::iterator max, n;
    DNANSeq *d;

    this->lnl = 0;

    for (col = 0; col < length; col++) {
        nuc['A'] = 0; 
        nuc['G'] = 0;
        nuc['C'] = 0;
        nuc['T'] = 0;

        //
        // Count the nucleotide type at each position in the column.
        //
        for (row = 0; row < height; row++) {
	    d = this->pmatrix[row];
            nuc[(*d)[col]]++;
	}
        //
        // Find the base with a plurality of occurances and call it.
        //
        max = nuc.end();
        tot = 0;
        for (n = nuc.begin(); n != nuc.end(); n++) {
            tot += n->second;
            if (max == nuc.end() || n->second > max->second)
                max = n;
        }

	//
	// For nucleotide positions with potential polymorphism (i.e. two or more alleles at 
	// the locus that differ at that position), first find the ML genotype (call_multinomial_snp). 
	// If it returns 'het' calculate the heterozygous_likelihood(), otherwise calculate homozygous
	// likelihood.
	//
	snp_type res = this->snps[col]->type;

	if (res == snp_type_het) 
	    this->lnl += heterozygous_likelihood(col, nuc);
	else if (res == snp_type_hom)
	    this->lnl += homozygous_likelihood(col, nuc);
	else {
	    double homlnl = homozygous_likelihood(col, nuc);
	    double hetlnl = heterozygous_likelihood(col, nuc);
	    this->lnl += hetlnl > homlnl ? hetlnl : homlnl;
	}
    }

    return this->lnl;
}

string MergedStack::write_cmb() {

    stringstream s;
    uint size = this->utags.size();

    s << "{";
    for (uint i = 0; i < size; i++) {
        s << this->utags[i];
        if (i < size - 1)
            s << ", ";
    }
    s << "}";

    return s.str();
}
