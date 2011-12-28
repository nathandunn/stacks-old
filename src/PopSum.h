// -*-mode:c++; c-style:k&r; c-basic-offset:4;-*-
//
// Copyright 2011, Julian Catchen <jcatchen@uoregon.edu>
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

#ifndef __POPSUM_H__
#define __POPSUM_H__

#include <string.h>
#include <string>
using std::string;
#include <vector>
using std::vector;
#include <map>
using std::map;
#include <set>
using std::set;
#include <algorithm>
#include <utility>
using std::pair;
using std::make_pair;

#include "stacks.h"

extern int progeny_limit;

class SumStat {
public:
    double num_indv;
    double p;
    double obs_het;
    double obs_hom;
    double exp_het;
    double exp_hom;

    SumStat() { 
	num_indv = 0.0;
	p        = 0.0;
	obs_het  = 0.0;
	obs_hom  = 0.0;
	exp_het  = 0.0;
	exp_hom  = 0.0;
    }
};

class LocSum {
public:
    SumStat *nucs; // Array containing summary statistics for 
                   // each nucleotide position at this locus.

    LocSum(int len)  { 
	this->nucs = new SumStat[len]; 
    }
    ~LocSum() {
	delete [] this->nucs;
    }
};

//
// Population Summary class contains a two dimensional array storing the 
// summary statistics for each locus in each of a set of populations:
//
//          Pop1    Pop2    Pop3  
// Locus1  SumStat SumStat SumStat
// Locus2  SumStat SumStat SumStat
// Locus3  SumStat SumStat SumStat
// Locus4  SumStat SumStat SumStat
// ...
//
template<class LocusT=Locus>
class PopSum {
    int       num_loci;
    int       num_pops;
    LocSum ***data;
    map<int, int> locus_order;  // LocusID => ArrayIndex; map catalog IDs to their first dimension 
                                // position in the LocSum array.
    map<int, int> rev_locus_order;
    map<int, int> pop_order;    // PopulationID => ArrayIndex; map defining at what position in 
                                // the second dimension of the LocSum array each population is stored.
    map<int, int> rev_pop_order;

public:
    PopSum(int, int);
    ~PopSum();

    int initialize(PopMap<LocusT> *);
    int add_population(map<int, LocusT *> &, PopMap<LocusT> *, uint, uint, uint);

    int loci_cnt() { return this->num_loci; }
    int rev_locus_index(int index) { return this->rev_locus_order[index]; }
    int pop_cnt()  { return this->num_pops; }
    int rev_pop_index(int index) { return this->rev_pop_order[index]; }

    LocSum **locus(int);
    LocSum  *pop(int, int);

private:
    int tally_heterozygous_pos(LocusT *, Datum **, LocSum *, int, int, uint, uint);
    int tally_fixed_pos(LocusT *, Datum **, LocSum *, int, uint, uint);
};

template<class LocusT>
PopSum<LocusT>::PopSum(int num_loci, int num_populations) {
    this->data = new LocSum **[num_loci];

    for (int i = 0; i < num_loci; i++) {
	this->data[i] = new LocSum *[num_populations];

	for (int j = 0; j < num_populations; j++)
	    this->data[i][j] = NULL;
    }

    this->num_pops = num_populations;
    this->num_loci = num_loci;
}

template<class LocusT>
PopSum<LocusT>::~PopSum() {
    for (int i = 0; i < this->num_loci; i++) {
	for (int j = 0; j < this->num_pops; j++)
	    delete this->data[i][j];
	delete [] this->data[i];
    }
    delete [] this->data;
}

template<class LocusT>
int PopSum<LocusT>::initialize(PopMap<LocusT> *pmap) {
    int locus_id;

    for (int i = 0; i < this->num_loci; i++) {
	locus_id = pmap->rev_locus_index(i);
	this->locus_order[locus_id] = i;
	this->rev_locus_order[i]    = locus_id;
    }

    return 0;
}

template<class LocusT>
int PopSum<LocusT>::add_population(map<int, LocusT *> &catalog,
				   PopMap<LocusT> *pmap, 
				   uint population_id,
				   uint start_index, uint end_index) {
    LocusT  *loc;
    Datum  **d;
    LocSum **s;
    uint locus_id, len;
    set<int> snp_cols;

    //
    // Determine the index for this population
    //
    uint pop_index = this->pop_order.size() == 0 ? 0 : this->pop_order.size();
    this->pop_order[population_id] = pop_index;
    this->rev_pop_order[pop_index] = population_id;

    for (int i = 0; i < this->num_loci; i++) {
	locus_id = pmap->rev_locus_index(i);
	d   = pmap->locus(locus_id);
	s   = this->locus(locus_id);
	loc = catalog[locus_id];
	//
	// Create an array of SumStat objects
	//
	len = strlen(loc->con);
	s[pop_index] = new LocSum(len);

	//
	// The catalog records which nucleotides are heterzygous. For these nucleotides we much
	// calculate observed genotype frequencies, allele frequencies, and expected genotype frequencies.
	//
	for (uint k = 0; k < loc->snps.size(); k++) {
	    this->tally_heterozygous_pos(loc, d, s[pop_index], 
					 loc->snps[k]->col, k, start_index, end_index);
	    snp_cols.insert(loc->snps[k]->col);
	}
	//
	// For all other fixed sites, we just need to record them.
	//
	for (uint k = 0; k < len; k++) {
	    if (snp_cols.count(k)) continue;
	    this->tally_fixed_pos(loc, d, s[pop_index], 
				  k, start_index, end_index);
	    
	}

	snp_cols.clear();
    }

    return 0;
}

template<class LocusT>
int PopSum<LocusT>::tally_fixed_pos(LocusT *locus, Datum **d, LocSum *s, int pos, uint start, uint end) {
    double num_indv = 0.0;

    for (uint i = start; i <= end; i++) {
	if (d[i] == NULL) continue;
	//
	// Before counting this individual, make sure the model definitively called this 
	// position as hEterozygous or hOmozygous.
	//
	if (d[i]->model[pos] == 'E' || d[i]->model[pos] == 'O')
	    num_indv++;
    }
    //
    // Record the results in the PopSum object.
    //
    s->nucs[pos].num_indv = num_indv;
    s->nucs[pos].p        = 1.0;
    s->nucs[pos].obs_hom  = 1.0;
    s->nucs[pos].obs_het  = 0.0;
    s->nucs[pos].exp_hom  = 1.0;
    s->nucs[pos].exp_het  = 0.0;

    return 0;
}

template<class LocusT>
int PopSum<LocusT>::tally_heterozygous_pos(LocusT *locus, Datum **d, LocSum *s, 
					   int pos, int snp_index, uint start, uint end) {
    int  nucs[4] = {0};
    uint i;
    char nuc;

    //cerr << "  Calculating summary stats at het locus " << locus->id << " position " << pos << "; snp_index: " << snp_index << "\n";

    //
    // Iterate over each individual in this population.
    //
    for (i = start; i <= end; i++) {
	if (d[i] == NULL) continue;

	//
	// Pull each allele for this SNP from the observed haplotype.
	//
	for (uint j = 0; j < d[i]->obshap.size(); j++) {
	    nuc = d[i]->obshap[j][snp_index];

	    switch(nuc) {
	    case 'A':
	    case 'a':
		nucs[0]++;
		break;
	    case 'C':
	    case 'c':
		nucs[1]++;
		break;
	    case 'G':
	    case 'g':
		nucs[2]++;
		break;
	    case 'T':
	    case 't':
		nucs[3]++;
		break;
	    }
	}
    }

    //
    // Determine how many alleles are present at this position in this population.
    // We cannot deal with more than two alternative alleles, if there are more than two
    // in a single population, print a warning and exclude this nucleotide position.
    //
    int allele_cnt = 0;
    for (i = 0; i < 4; i++)
	if (nucs[i] > 0) allele_cnt++;

    if (allele_cnt > 2) {
	cerr << "Locus " << locus->id << ", bp: " << pos << "; more than two alleles present.\n";
	return 0;
    }

    //
    // Record which nucleotide is the P allele and which is the Q allele.
    //
    char p_allele = 0;
    char q_allele = 0;

    i = 0;
    while (p_allele == 0 && i < 4) {
	if (nucs[i] > 0) {
	    switch(i) {
	    case 0:
		p_allele = 'A';
		break;
	    case 1:
		p_allele = 'C';
		break;
	    case 2:
		p_allele = 'G';
		break;
	    case 3:
		p_allele = 'T';
		break;
	    }
	}
	i++;
    }
    while (q_allele == 0 && i < 4) {
	if (nucs[i] > 0) {
	    switch(i) {
	    case 1:
		q_allele = 'C';
		break;
	    case 2:
		q_allele = 'G';
		break;
	    case 3:
		q_allele = 'T';
		break;
	    }
	}
	i++;
    }
    cerr << "  P Allele: " << p_allele << "; Q Allele: " << q_allele << "\n";

    //
    // Calculate obaserved genotype frequencies.
    //
    double num_indv = 0.0;
    double obs_het  = 0.0;
    double obs_p    = 0.0;
    double obs_q    = 0.0;

    for (i = start; i <= end; i++) {
	if (d[i] == NULL) continue;
	//
	// Before counting this individual, make sure the model definitively called this 
	// position as hEterozygous or hOmozygous.
	//
	if (d[i]->model[pos] == 'E' || d[i]->model[pos] == 'O')
	    num_indv++;
	else
	    continue;

	if (d[i]->obshap.size() > 1 &&
	    d[i]->obshap[0][snp_index] != d[i]->obshap[1][snp_index])
	    obs_het++;
	else if (d[i]->obshap[0][snp_index] == p_allele)
	    obs_p++;
	else if (d[i]->obshap[0][snp_index] == q_allele)
	    obs_q++;
    }
    cerr << "  Num Individuals: " << num_indv << "; Obs Hets: " << obs_het << "; Obs P: " << obs_p << "; Obs Q: " << obs_q << "\n";

    if (num_indv == 0 || num_indv < progeny_limit) return 0;

    //
    // Calculate allele frequencies: (1/2 Het + Hom) / N
    //
    double allele_p = ((0.5 * obs_het) + obs_p) / num_indv;
    double allele_q = ((0.5 * obs_het) + obs_q) / num_indv;

    cerr << "  P allele frequency: " << allele_p << "; Q allele frequency: " << allele_q << "\n";

    //
    // Calculate expected genotype frequencies.
    //
    double exp_het = 2 * allele_p * allele_q; // 2pq
    double exp_p   = allele_p * allele_p;     // p^2
    double exp_q   = allele_q * allele_q;     // q^2; 

    cerr << "  Expected Het: " << exp_het << "; Expected P: " << exp_p << "; Expected Q: " << exp_q << "\n";

    obs_het = obs_het / num_indv;
    obs_p   = obs_p   / num_indv;
    obs_q   = obs_q   / num_indv;

    cerr << "  Obs Hets Freq: " << obs_het << "; Obs P Freq: " << obs_p << "; Obs Q Freq: " << obs_q << "\n";

    //
    // Record the results in the PopSum object.
    //
    s->nucs[pos].num_indv = num_indv;
    s->nucs[pos].p        = allele_p > allele_q ? allele_p : allele_q;
    s->nucs[pos].obs_hom  = obs_p    > obs_q    ? obs_p    : obs_q;
    s->nucs[pos].obs_het  = obs_het;
    s->nucs[pos].exp_hom  = obs_p    > obs_q    ? exp_p    : exp_q;
    s->nucs[pos].exp_het  = exp_het;

    return 0;
}

template<class LocusT>
LocSum **PopSum<LocusT>::locus(int locus) {
    return this->data[this->locus_order[locus]];
}

template<class LocusT>
LocSum  *PopSum<LocusT>::pop(int locus, int pop_id) {
    return this->data[this->locus_order[locus]][this->pop_order[pop_id]];
}

#endif // __POPSUM_H__
