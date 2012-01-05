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

class PopPair {
public:
    int    bp;
    double pi;
    double fst;
    double wfst; // Weigted Fst (kernel-smoothed)

    PopPair() { 
	bp   = 0;
	pi   = 0.0;
	fst  = 0.0;
	wfst = 0.0;
    }
};

class SumStat {
public:
    double num_indv;
    char   p_nuc;
    char   q_nuc;
    double p;
    double obs_het;
    double obs_hom;
    double exp_het;
    double exp_hom;
    double H;
    double pi;
    double Fis;

    SumStat() { 
	num_indv = 0.0;
	p        = 0.0;
	p_nuc    = 0;
	q_nuc    = 0;
	obs_het  = 0.0;
	obs_hom  = 0.0;
	exp_het  = 0.0;
	exp_hom  = 0.0;
	pi       = 0.0;
	Fis      = 0.0;
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
    PopPair *Fst(int, int, int, int);

private:
    int    tally_heterozygous_pos(LocusT *, Datum **, LocSum *, int, int, uint, uint);
    int    tally_fixed_pos(LocusT *, Datum **, LocSum *, int, uint, uint);
    double pi(double, double, double);
    double binomial_coeff(double, double);
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
PopPair *PopSum<LocusT>::Fst(int locus, int pop_1, int pop_2, int pos) 
{
    LocSum *s_1 = this->pop(locus, pop_1);
    LocSum *s_2 = this->pop(locus, pop_2);

    //
    // If this locus only appears on one population, do not calculate Fst.
    //
    if (s_1->nucs[pos].num_indv == 0 || s_2->nucs[pos].num_indv == 0) return NULL;

    //
    // Calculate Fst at a locus, sub-population relative to that found in the entire population
    //   Fst = 1 - (Sum_j( (n_j choose 2) * pi_j)) / (pi_all * Sum_j( (n_j choose 2) ))
    //
    double n_1, n_2, pi_1, pi_2;

    n_1  = s_1->nucs[pos].num_indv * 2;
    n_2  = s_2->nucs[pos].num_indv * 2;
    pi_1 = s_1->nucs[pos].pi;
    pi_2 = s_2->nucs[pos].pi;

    if (pi_1 == 0 && pi_2 == 0 && s_1->nucs[pos].p_nuc == s_2->nucs[pos].p_nuc)
	return NULL;

    //
    // Calculate Pi over the entire pooled population.
    //
    // First, make sure this site is compatible between the two populations (no more than two alleles).
    //
    char nucs[4];
    int  ncnt[4] = {0};

    nucs[0] = s_1->nucs[pos].p_nuc;
    nucs[1] = s_1->nucs[pos].q_nuc;
    nucs[2] = s_2->nucs[pos].p_nuc;
    nucs[3] = s_2->nucs[pos].q_nuc;

    for (int i = 0; i < 4; i++) 
	switch(nucs[i]) {
	case 'A':
	    ncnt[0]++;
	    break;
	case 'C':
	    ncnt[1]++;
	    break;
	case 'G':
	    ncnt[2]++;
	    break;
	case 'T':
	    ncnt[3]++;
	    break;
	}

    int allele_cnt = 0;
    for (int i = 0; i < 4; i++)
	if (ncnt[i] > 0) allele_cnt++;

    if (allele_cnt > 2) {
	cerr << "Incompatible site for Fst calculation at locus " << locus << " between population " 
	     << pop_1 << " and population " << pop_2 << " position " << pos << "\n";
	return NULL;
    }

    double tot_alleles = n_1 + n_2;
    double p_1 = n_1 * s_1->nucs[pos].p;
    double q_1 = n_1 - p_1;
    double p_2 = 
	s_1->nucs[pos].p_nuc == s_2->nucs[pos].p_nuc ? 
	s_2->nucs[pos].p : (1 - s_2->nucs[pos].p);
    p_2 = n_2 * p_2;
    double q_2 = n_2 - p_2;

    double pi_all = this->pi(tot_alleles, p_1 + p_2, q_1 + q_2);

    double bcoeff_1 = this->binomial_coeff(n_1, 2);
    double bcoeff_2 = this->binomial_coeff(n_2, 2);

    double num = (bcoeff_1 * pi_1) + (bcoeff_2 * pi_2);
    double den = pi_all * (bcoeff_1 + bcoeff_2);

    double Fst = 1 - (num / den);

    PopPair *pair = new PopPair();
    pair->fst = Fst;
    pair->pi  = pi_all;

    // cerr << "Locus: " << locus << ", pos: " << pos << "\n"
    // 	 << "    p_1.nuc: " << s_1->nucs[pos].p_nuc << "; q_1.nuc: " << s_1->nucs[pos].q_nuc 
    // 	 << "; p_2.nuc: " << s_2->nucs[pos].p_nuc << "; q_2.nuc: " << s_2->nucs[pos].q_nuc << "\n"
    // 	 << "    Total alleles: " << tot_alleles << "; " << " s_1.p: " << s_1->nucs[pos].p 
    // 	 << "; s_2.p: " << s_2->nucs[pos].p << "\n"
    // 	 << "    p_1: " << p_1 << "; q_1: " << q_1 << " p_2: " << p_2 << "; q_2: " << q_2 << "\n"
    // 	 << "    Pi1: " << pi_1 << "; Pi2: " << pi_2 << "; PiAll: " << pi_all << "\n"
    // 	 << "    N1: " << n_1 << "; N1 choose 2: " << bcoeff_1 << "\n"
    // 	 << "    N2: " << n_2 << "; N2 choose 2: " << bcoeff_2 << "\n"
    // 	 << "  Fst: " << Fst << "\n";

    return pair;
}

template<class LocusT>
int PopSum<LocusT>::tally_fixed_pos(LocusT *locus, Datum **d, LocSum *s, int pos, uint start, uint end) 
{
    double num_indv = 0.0;
    char   p_nuc;

    for (uint i = start; i <= end; i++) {
	if (d[i] == NULL) continue;
	//
	// Before counting this individual, make sure the model definitively called this 
	// position as hEterozygous or hOmozygous.
	//
	if (d[i]->model[pos] == 'E') {
	    cerr << "Warning: heterozygous model call at fixed nucleotide position: " 
		 << "locus " << locus->id << " individual " << d[i]->id << "; position: " << pos << "\n";
	} else if (d[i]->model[pos] == 'O') {
	    num_indv++;
	    p_nuc = locus->con[pos];
	}
    }
    //
    // Record the results in the PopSum object.
    //
    s->nucs[pos].num_indv = num_indv;
    s->nucs[pos].p        = 1.0;
    s->nucs[pos].p_nuc    = p_nuc;
    s->nucs[pos].obs_hom  = 1.0;
    s->nucs[pos].obs_het  = 0.0;
    s->nucs[pos].exp_hom  = 1.0;
    s->nucs[pos].exp_het  = 0.0;
    s->nucs[pos].pi       = 0.0;

    return 0;
}

template<class LocusT>
int PopSum<LocusT>::tally_heterozygous_pos(LocusT *locus, Datum **d, LocSum *s, 
					   int pos, int snp_index, uint start, uint end) 
{
    //
    // Tally up the genotype frequencies.
    //
    int  nucs[4] = {0};
    uint i;
    char nuc;

    //cerr << "  Calculating summary stats at het locus " << locus->id << " position " << pos << "; snp_index: " << snp_index << "\n";

    //
    // Iterate over each individual in this sub-population.
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
    //cerr << "  P Allele: " << p_allele << "; Q Allele: " << q_allele << "\n";

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
    //cerr << "  Num Individuals: " << num_indv << "; Obs Hets: " << obs_het << "; Obs P: " << obs_p << "; Obs Q: " << obs_q << "\n";

    if (num_indv == 0 || num_indv < progeny_limit) return 0;

    //
    // Calculate total number of alleles
    //
    double tot_alleles = num_indv * 2;
    double allele_p    = obs_het + (2 * obs_p);
    double allele_q    = obs_het + (2 * obs_q);

    //
    // Calculate Pi, equivalent to expected heterozygosity (exp_het)
    //
    s->nucs[pos].pi = this->pi(tot_alleles, allele_p, allele_q);

    //
    // Convert to allele frequencies
    //
    allele_p = allele_p / tot_alleles;
    allele_q = allele_q / tot_alleles;

    //cerr << "  P allele frequency: " << allele_p << "; Q allele frequency: " << allele_q << "\n";

    //
    // Calculate expected genotype frequencies.
    //
    double exp_het = 2 * allele_p * allele_q; // 2pq
    double exp_p   = allele_p * allele_p;     // p^2
    double exp_q   = allele_q * allele_q;     // q^2; 

    //cerr << "  Expected Het: " << exp_het << "; Expected P: " << exp_p << "; Expected Q: " << exp_q << "\n";

    obs_het = obs_het / num_indv;
    obs_p   = obs_p   / num_indv;
    obs_q   = obs_q   / num_indv;

    //cerr << "  Obs Hets Freq: " << obs_het << "; Obs P Freq: " << obs_p << "; Obs Q Freq: " << obs_q << "\n";

    //
    // Record the results in the PopSum object.
    //
    s->nucs[pos].num_indv = num_indv;
    s->nucs[pos].p        = allele_p > allele_q ? allele_p : allele_q;
    s->nucs[pos].p_nuc    = allele_p > allele_q ? p_allele : q_allele;
    s->nucs[pos].q_nuc    = allele_p > allele_q ? q_allele : p_allele;
    s->nucs[pos].obs_hom  = obs_p    > obs_q    ? obs_p    : obs_q;
    s->nucs[pos].obs_het  = obs_het;
    s->nucs[pos].exp_hom  = obs_p    > obs_q    ? exp_p    : exp_q;
    s->nucs[pos].exp_het  = exp_het;

    //
    // Calculate F_is, the inbreeding coefficient of an individual (I) relative to the subpopulation:
    //   Fis = (exp_het - obs_het) / exp_het
    //
    s->nucs[pos].Fis = s->nucs[pos].pi == 0 ? 0 : (s->nucs[pos].pi - obs_het) / s->nucs[pos].pi;

    return 0;
}

template<class LocusT>
double PopSum<LocusT>::pi(double tot_alleles, double p, double q)
{
    //
    // Calculate Pi, equivalent to expected heterozygosity:
    //  pi = 1 - Sum_i( (n_i choose 2) ) / (n choose 2)
    //
    double pi = 
	this->binomial_coeff(p, 2) + 
	this->binomial_coeff(q, 2);
    pi = pi / binomial_coeff(tot_alleles, 2);
    pi = 1 - pi;

    return pi;
}

template<class LocusT>
double PopSum<LocusT>::binomial_coeff(double n, double k)
{
    if (n < k) return 0.0;
    //
    // Compute the binomial coefficient using the method of:
    // Y. Manolopoulos, "Binomial coefficient computation: recursion or iteration?", 
    // ACM SIGCSE Bulletin, 34(4):65-67, 2002.
    //
    double r = 1.0;
    double s = (k < n - k) ? n - k + 1 : k + 1;

    for (double i = n; i >= s; i--)
        r = r * i / (n - i + 1);

    return r;
}

template<class LocusT>
LocSum **PopSum<LocusT>::locus(int locus) 
{
    return this->data[this->locus_order[locus]];
}

template<class LocusT>
LocSum  *PopSum<LocusT>::pop(int locus, int pop_id) 
{
    return this->data[this->locus_order[locus]][this->pop_order[pop_id]];
}

#endif // __POPSUM_H__
