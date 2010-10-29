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
// $Id$
//

#include "mstack.h"
#include "models.h"

MergedStack::MergedStack()  { 
    id         = 0;
    count      = 0;
    con        = NULL;
    matrix     = NULL;
    likelihood = 0.0;
    loc.bp     = 0; 
    loc.chr[0] = '\0';
    deleveraged     = false;
    masked          = false;
    blacklisted     = false;
    lumberjackstack = false;
}

MergedStack::~MergedStack() { 
    delete [] con;

    for (uint i = 0; i < kmers.size(); i++)
        delete [] kmers[i];
    for (uint i = 0; i < snps.size(); i++)
        delete snps[i];

    delete [] matrix;
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
    uint cnt = this->count + this->remtags.size();
    if (this->matrix != NULL)
        delete [] this->matrix;
    this->matrix = new char * [cnt];

    vector<int>::iterator j;
    int i = 0;
    for (j = this->utags.begin(); j != this->utags.end(); j++) {
        tag = unique[*j];

        for (uint k = 0; k < tag->count; k++) {
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

double MergedStack::calc_likelihood() {

    if (this->matrix == NULL)
        return 0;

    //
    // Iterate over each column of the array and call the consensus base.
    //
    int row, col, tot;
    int length = strlen(this->matrix[0]);
    int height = this->count + this->remtags.size();
    string con;
    map<char, int> nuc;
    map<char, int>::iterator max, n;

    this->likelihood = 0;

    for (col = 0; col < length; col++) {
        nuc['A'] = 0; 
        nuc['C'] = 0;
        nuc['G'] = 0;
        nuc['T'] = 0;

        //
        // Count the nucleotide type at each position in the column.
        //
        for (row = 0; row < height; row++)
            nuc[this->matrix[row][col]]++;

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
        con += max->first;

        if (max->second == tot) {
            //
            // For nucleotide positions with no polymorphism (i.e. only one allele at the locus in 
            // this permutation, or alleles are identical at the nucleotide position), the lnL for 
            // an individual is just the output of homozygous_likelihood()
            //
            this->likelihood += homozygous_likelihood(col, nuc);

        } else {
            //
            // For nucleotide positions with potential polymorphism (i.e. two or more alleles at 
            // the locus that differ at that position), first find the ML genotype (call_multinomial_snp). 
            // If it returns 'het' calculate the heterozygous_likelihood(), otherwise calculate homozygous
            // likelihood.
            //
            allelet res = call_multinomial_snp(this, col, nuc);

            if (res == het) 
                this->likelihood += heterozygous_likelihood(col, nuc);
            else if (res == hom)
                this->likelihood += homozygous_likelihood(col, nuc);
            else {
                double homlln = homozygous_likelihood(col, nuc);
                double hetlln = heterozygous_likelihood(col, nuc);
                this->likelihood += hetlln > homlln ? hetlln : homlln;
            }
        }
    }

    return this->likelihood;
}

