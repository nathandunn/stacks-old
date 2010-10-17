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

// -*-mode:c++; c-style:k&r; c-basic-offset:4;-*-
//
// models.cc -- routines to detect polymorphism (snp) and detect a lack of polymorphism (fixed).
//
// Julian Catchen
// jcatchen@uoregon.edu
// University of Oregon
//
// $Id$
//
#include "models.h"

int call_multinomial_snp (MergedStack *tag, int col, map<char, int> &n) {
    const float heterozygote_limit = -3.84;
    const float homozygote_limit   =  3.84;

    vector<pair<char, int> > nuc;
    map<char, int>::iterator i;

    int total = 0;
    for (i = n.begin(); i != n.end(); i++) {
	total += i->second;
	nuc.push_back(make_pair(i->first, i->second));
    }

    if (total < snp_min)
	return 0;

    sort(nuc.begin(), nuc.end(), compare_pair);

    //
    // Method of Paul Hohenlohe <hohenlo@uoregon.edu>, personal communication.
    //
    // For a diploid individual, there are ten possible genotypes
    // (four homozygous and six heterozygous genotypes).  We calculate
    // the likelihood of each possible genotype by using a multinomial
    // sampling distribution, which gives the probability of observing
    // a set of read counts (n1,n2,n3,n4) given a particular genotype.
    //
    float nuc_1   = nuc[0].second;
    float nuc_2   = nuc[1].second;
    float nuc_3   = nuc[2].second;
    float nuc_4   = nuc[3].second;
    float l_ratio = 0;

    l_ratio = (nuc_1 * log(nuc_1 / total));

    if (total - nuc_1 > 0)
    	l_ratio += ((total - nuc_1) * log((total - nuc_1) / (3 * total)));

    if (nuc_1 + nuc_2 > 0)
	l_ratio -= ((nuc_1 + nuc_2) * log((nuc_1 + nuc_2) / (2 * total)));

    if (nuc_3 + nuc_4 > 0)
	l_ratio -= ((nuc_3 + nuc_4) * log((nuc_3 + nuc_4) / (2 * total)));

    l_ratio *= 2;

    //cerr << "Nuc_1: " << nuc_1 << " Nuc_2: " << nuc_2 << " Nuc_3: " << nuc_3 << " Nuc_4: " << nuc_4 << " Likelihood ratio: " << l_ratio << "\n";

    if (l_ratio <= heterozygote_limit) {
        // This locus is a heterozygote.
	// Record this SNP
	SNP *snp = new SNP;
	snp->col    = col;
	snp->lratio = l_ratio;
	snp->rank_1 = nuc[0].first;
	snp->rank_2 = nuc[1].first;

	tag->snps.push_back(snp);

    } else if (l_ratio >= homozygote_limit) {
        // This locus is a homozygote.

    } else {
        // Unknown whether this is a heterozygote or homozygote.

    }

    return 0;
}

int call_multinomial_fixed (MergedStack *tag, int col, map<char, int> &n) {
    const double nucleotide_fixed_limit = 1.92;

    vector<pair<char, int> > nuc;
    map<char, int>::iterator i;

    int total = 0;
    for (i = n.begin(); i != n.end(); i++) {
	total += i->second;
	nuc.push_back(make_pair(i->first, i->second));
    }

    if (total < snp_min)
	return 0;

    sort(nuc.begin(), nuc.end(), compare_pair);

    //
    // Method of Paul Hohenlohe <hohenlo@uoregon.edu>, personal communication.
    //
    // Each population sample contains DNA from 6 individuals, so a
    // sample of 12 alleles from the population. We want to assign a
    // nucleotide (A,C,G,T) to each position where the population is
    // fixed or nearly so, and N to each position that is either
    // polymorphic within the population or has insufficient coverage
    // depth to make a call. We can do this with a likelihood ratio
    // test of the read counts, testing whether the allele frequency
    // of the dominant allele is significantly larger than some
    // threshold p) , stepping through each nucleotide position across
    // RAD tags.
    //
    double nuc_1   = nuc[0].second;
    double nuc_2   = nuc[1].second;
    double n_ratio = 0.0;
    double l_ratio = 0.0;
    double epsilon = -1 * (log(1 - barcode_err_freq) / barcode_size);

    n_ratio = nuc_1 / (nuc_1 + nuc_2);

    l_ratio  = 
	nuc_1 * log( ((4 * nuc_1 * (1 - epsilon)) + ((nuc_1 + nuc_2) * epsilon)) / 
		     ((4 * p_freq * (nuc_1 + nuc_2) * (1 - epsilon)) + ((nuc_1 + nuc_2) * epsilon)) );

    l_ratio += 
	nuc_2 * log( ((4 * nuc_2 * (1 - epsilon)) + ((nuc_1 + nuc_2) * epsilon)) / 
		     ((4 * (1 - p_freq) * (nuc_1 + nuc_2) * (1 - epsilon)) + ((nuc_1 + nuc_2) * epsilon)) );

    //cerr << "Nuc_1: " << nuc_1 << " Nuc_2: " << nuc_2 << " Likelihood ratio: " << l_ratio << "\n";

    if (n_ratio < p_freq || l_ratio < nucleotide_fixed_limit) {
	// Record this likely-SNP
	SNP *snp = new SNP;
	snp->col    = col;
	snp->lratio = l_ratio;
	snp->rank_1 = nuc[0].first;
	snp->rank_2 = nuc[1].first;

	tag->snps.push_back(snp);
    }

    return 0;
}

bool compare_pair(pair<char, int> a, pair<char, int> b) {
    return (a.second > b.second);
}
