// -*-mode:c++; c-style:k&r; c-basic-offset:4;-*-
//
// Copyright 2011-2012, Julian Catchen <jcatchen@uoregon.edu>
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
#include <cstdint>
#include <math.h>

#include "stacks.h"

extern bool   log_fst_comp;
extern double minor_allele_freq;
extern map<int, string> pop_key;
const  uint   PopStatSize = 5;

class PopStat {
public:
    int    loc_id;
    int    bp;
    bool   fixed;
    double alleles;    // Number of alleles sampled at this location.
    uint   snp_cnt;    // Number of SNPs in kernel-smoothed window centered on this SNP.
    double stat[PopStatSize];
    double smoothed[PopStatSize];
    double bs[PopStatSize];

    PopStat() {
	this->loc_id  = 0;
	this->bp      = 0;
	this->fixed   = false;
	this->alleles = 0.0;
	this->snp_cnt = 0;

	for (uint i = 0; i < PopStatSize; i++) {
	    this->stat[i]     = 0.0;
	    this->smoothed[i] = 0.0;
	    this->bs[i]       = 0.0;
	}
     }
    virtual ~PopStat() {
    }
};

class HapStat: public PopStat {
    // PopStat[0]: Phi_st
    // PopStat[1]: Phi_ct
    // PopStat[2]: Phi_sc
    // PopStat[3]: Fst'
    // PopStat[4]: D_est
public:
    double *comp;
    uint    popcnt;

    HapStat(): PopStat() {
	comp = NULL;
    }
    ~HapStat() {
	if (this->comp != NULL)
	    delete [] comp;
    }
};

class LocStat: public PopStat {
    // PopStat[0]: gene diversity
    // PopStat[1]: haplotype diversity (Pi)
public:
    uint     hap_cnt; // Number of unique haplotypes at this locus. 
    string   hap_str; // Human-readable string of haplotype counts.

    LocStat(): PopStat() { 
	this->hap_cnt = 0;
    }
    ~LocStat() {};
};

class PopPair: public PopStat {
    // PopStat[0]: corrected Fst, (by p-value or Bonferroni p-value).
    // PopStat[1]: corrected AMOVA Fst
public:
    int     col;
    double  pi;
    double  fst;
    double  fet_p;      // Fisher's Exact Test p-value.
    double  fet_or;     // Fisher's exact test odds ratio.
    double  or_se;      // Fisher's exact test odds ratio standard error.
    double  lod;        // base 10 logarithm of odds score.
    double  ci_low;     // Fisher's exact test lower confidence interval.
    double  ci_high;    // Fisher's exact test higher confidence interval.
    double  amova_fst;  // AMOVA Fst method, from Weir, Genetic Data Analysis II .
    double *comp;

    PopPair() { 
	col       = 0;
	pi        = 0.0;
	fst       = 0.0;
	fet_p     = 0.0;
	fet_or    = 0.0;
	or_se     = 0.0;
	lod       = 0.0;
	ci_low    = 0.0;
	ci_high   = 0.0;
	amova_fst = 0.0;
	comp      = NULL;
    }
    ~PopPair() {
	if (this->comp != NULL)
	    delete [] comp;
    }
};

class SumStat: public PopStat {
    // PopStat[0]: pi
    // PopStat[1]: fis
public:
    bool    incompatible_site;
    bool    filtered_site;
    double  num_indv;
    char    p_nuc;
    char    q_nuc;
    double  p;
    double  obs_het;
    double  obs_hom;
    double  exp_het;
    double  exp_hom;
    double &pi;

    SumStat(): PopStat(), pi(this->stat[0]) {
	num_indv  = 0.0;
	p         = 0.0;
	p_nuc     = 0;
	q_nuc     = 0;
	obs_het   = 0.0;
	obs_hom   = 0.0;
	exp_het   = 0.0;
	exp_hom   = 0.0;
	snp_cnt   = 0;
	incompatible_site = false;
	filtered_site     = false;
    }
};

class LocSum {
public:
    SumStat *nucs;    // Array containing summary statistics for 
                      // each nucleotide position at this locus.
    LocSum(int len) {
	this->nucs = new SumStat[len]; 
    }
    ~LocSum() {
	delete [] this->nucs;
    }
};

class NucTally {
public:
    int      loc_id;
    int      bp;
    uint16_t col;
    uint16_t num_indv;
    uint16_t pop_cnt;
    uint16_t allele_cnt;
    char     p_allele;
    char     q_allele;
    double   p_freq;
    double   obs_het;
    bool     fixed;
    int      priv_allele;

    NucTally() { 
	loc_id      = 0;
	bp          = 0;
	col         = 0;
	num_indv    = 0;
	pop_cnt     = 0;
	allele_cnt  = 0;
	p_allele    = 0;
	q_allele    = 0;
	p_freq      = 0.0;
	obs_het     = 0.0;
	priv_allele = -1;
	fixed       = true;
    }
};

class LocTally {
public:
    NucTally *nucs;

    LocTally(int len)  { 
	this->nucs = new NucTally[len]; 
    }
    ~LocTally() {
	delete [] this->nucs;
    }
};

//
// Population Summary class contains a two dimensional array storing the 
// summary statistics for each locus in each of a set of populations:
//
//           Pop1        Pop2        Pop3
// Locus1  +-LocSum----+-LocSum----+-LocSum
//         |  |           |           |
//         |  +-SumStat   +-SumStat   +-SumStat (Nuc0)
//         |  +-SumStat   +-SumStat   +-SumStat (Nuc1)
//         |  ...
// Locus2  +-LocSum----+-LocSum----+-LocSum
//         |  |           |           |
//         |  +-SumStat   +-SumStat   +-SumStat (Nuc0)
//         |  +-SumStat   +-SumStat   +-SumStat (Nuc1)
//         |  ...
// Locus3  +-LocSum----+-LocSum----+-LocSum
//         |  |           |           |
//         |  +-SumStat   +-SumStat   +-SumStat (Nuc0)
//         |  +-SumStat   +-SumStat   +-SumStat (Nuc1)
//         |  ...
//         ...
//
template<class LocusT=Locus>
class PopSum {
    int           num_loci;
    int           num_pops;
    LocSum     ***data;
    LocTally    **loc_tally;
    map<int, int> locus_order;  // LocusID => ArrayIndex; map catalog IDs to their first dimension 
                                // position in the LocSum array.
    map<int, int> rev_locus_order;
    map<int, int> pop_order;    // PopulationID => ArrayIndex; map defining at what position in 
                                // the second dimension of the LocSum array each population is stored.
    map<int, int> rev_pop_order;
    map<int, int> pop_sizes;    // The maximum size of each separate population.

public:
    PopSum(int, int);
    ~PopSum();

    int initialize(PopMap<LocusT> *);
    int add_population(map<int, LocusT *> &, PopMap<LocusT> *, uint, uint, uint, bool, ofstream &);
    int tally(map<int, LocusT *> &);

    int loci_cnt() { return this->num_loci; }
    int rev_locus_index(int index) { return this->rev_locus_order[index]; }
    int pop_cnt()  { return this->num_pops; }
    int pop_index(int index)     { return this->pop_order[index]; }
    int rev_pop_index(int index) { return this->rev_pop_order[index]; }
    int pop_size(int pop_id)     { return this->pop_sizes[pop_id]; }

    LocSum  **locus(int);
    LocSum   *pop(int, int);
    LocTally *locus_tally(int);
    PopPair  *Fst(int, int, int, int);
    int       fishers_exact_test(PopPair *, double, double, double, double);

private:
    int    tally_heterozygous_pos(LocusT *, Datum **, LocSum *, int, int, uint, uint);
    int    tally_fixed_pos(LocusT *, Datum **, LocSum *, int, uint, uint);
    int    tally_ref_alleles(LocSum **, int, short unsigned int &, char &, char &, short unsigned int &, short unsigned int &); 
    int    tally_observed_haplotypes(vector<char *> &, int);
    double pi(double, double, double);
    double binomial_coeff(double, double);
};

template<class LocusT>
PopSum<LocusT>::PopSum(int num_loci, int num_populations) {
    this->loc_tally = new LocTally *[num_loci];
    this->data      = new LocSum  **[num_loci];

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
	delete this->loc_tally[i];
    }
    delete [] this->data;
    delete [] this->loc_tally;
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
			       uint start_index, uint end_index, 
			       bool verbose, ofstream &log_fh) {
    LocusT  *loc;
    Datum  **d;
    LocSum **s;
    uint locus_id, len;
    int res;
    set<int> snp_cols;

    int incompatible_loci = 0;

    if (verbose)
	log_fh << "\n#\n# Recording sites that have incompatible loci -- loci with too many alleles present.\n"
	       << "#\n"
	       << "# Level\tAction\tLocus ID\tChr\tBP\tColumn\tPopID\n#\n";

    //
    // Determine the index for this population
    //
    uint pop_index = this->pop_order.size() == 0 ? 0 : this->pop_order.size();
    this->pop_order[population_id] = pop_index;
    this->rev_pop_order[pop_index] = population_id;

    //
    // Record the maximal size of this population.
    //
    this->pop_sizes[population_id] = end_index - start_index + 1;

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
	// Check if this locus has already been filtered and is NULL in all individuals.
	//
	bool filtered = true;
	for (uint k = start_index; k <= end_index; k++) {
	    if (d[k] != NULL) filtered = false;
	}
	if (filtered == true) {
	    for (uint k = 0; k < len; k++) {
		s[pop_index]->nucs[k].filtered_site = true;
	    }
	    continue;
	}

	//
	// The catalog records which nucleotides are heterozygous. For these nucleotides we will
	// calculate observed genotype frequencies, allele frequencies, and expected genotype frequencies.
	//
	for (uint k = 0; k < loc->snps.size(); k++) {
	    res = this->tally_heterozygous_pos(loc, d, s[pop_index], 
					       loc->snps[k]->col, k, start_index, end_index);
	    //
	    // If site is incompatible (too many alleles present), log it.
	    //
	    if (res < 0) {
		s[pop_index]->nucs[loc->snps[k]->col].incompatible_site = true;

		incompatible_loci++;
		if (verbose)
		    log_fh << "within_population\t"
			   << "incompatible_locus\t"
			   << loc->id << "\t"
			   << loc->loc.chr << "\t"
			   << loc->sort_bp(loc->snps[k]->col) << "\t"
			   << loc->snps[k]->col << "\t" 
			   << pop_key[population_id] << "\n";
	    }

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

    cerr << "Population '" << pop_key[population_id] << "' contained " << incompatible_loci << " incompatible loci -- more than two alleles present.\n";
    log_fh <<  "Population " << population_id << " contained " << incompatible_loci << " incompatible loci -- more than two alleles present.\n";

    return 0;
}

template<class LocusT>
int PopSum<LocusT>::tally(map<int, LocusT *> &catalog) 
{
    LocusT   *loc;
    LocSum  **s;
    LocTally *ltally;
    int       locus_id, variable_pop;
    uint16_t  p_cnt, q_cnt, len, col;

    for (int n = 0; n < this->num_loci; n++) {
	locus_id = this->rev_locus_index(n);
	loc      = catalog[locus_id];
	s        = this->locus(locus_id);
	len      = strlen(loc->con);

	ltally = new LocTally(len);
	this->loc_tally[n] = ltally;

	// for (uint i = 0; i < loc->snps.size(); i++) {
	//     uint col = loc->snps[i]->col;
	for (col = 0; col < len; col++) {

	    ltally->nucs[col].col       = col;
	    ltally->nucs[col].bp        = loc->sort_bp(col);
	    ltally->nucs[col].loc_id    = locus_id;

	    this->tally_ref_alleles(s, col, 
				    ltally->nucs[col].allele_cnt, 
				    ltally->nucs[col].p_allele, 
				    ltally->nucs[col].q_allele,
				    p_cnt, q_cnt);

	    //
	    // Is this site variable?
	    //
	    if (ltally->nucs[col].allele_cnt > 1)
		ltally->nucs[col].fixed = false;
	    
	    for (int j = 0; j < this->num_pops; j++) {
		//
		// Sum the number of individuals examined at this locus across populations.
		//
		ltally->nucs[col].num_indv += s[j]->nucs[col].num_indv;
		ltally->nucs[col].pop_cnt  += s[j]->nucs[col].num_indv > 0 ? 1 : 0;
	    }

	    for (int j = 0; j < this->num_pops; j++) {
		//
		// Sum the most frequent allele across populations.
		//
		if (s[j]->nucs[col].p_nuc == ltally->nucs[col].p_allele)
		    ltally->nucs[col].p_freq += 
			s[j]->nucs[col].p * (s[j]->nucs[col].num_indv / (double) ltally->nucs[col].num_indv);
		else 
		    ltally->nucs[col].p_freq += 
			(1 - s[j]->nucs[col].p) * (s[j]->nucs[col].num_indv / (double) ltally->nucs[col].num_indv);
		//
		// Sum observed heterozygosity across populations.
		//
		ltally->nucs[col].obs_het += 
		    s[j]->nucs[col].obs_het * (s[j]->nucs[col].num_indv / (double) ltally->nucs[col].num_indv);
	    }

	    //
	    // We want to report the most frequent allele as the P allele. Reorder the alleles 
	    // if necessary.
	    //
	    if (ltally->nucs[col].p_freq < 0.5) {
		char a = ltally->nucs[col].p_allele;
		ltally->nucs[col].p_allele = ltally->nucs[col].q_allele;
		ltally->nucs[col].q_allele = a;
		ltally->nucs[col].p_freq   = 1 - ltally->nucs[col].p_freq;
		uint b = p_cnt;
		p_cnt = q_cnt;
		q_cnt = b;
	    }

	    //
	    // Check if this is a private allele. Either the site is variable and
	    // the allele exists in one population, or the site is fixed and one
	    // population is homozygous for the private allele.
	    //
	    variable_pop = -1;

	    if (p_cnt == 1 && q_cnt > 1) {
	    	for (int j = 0; j < this->num_pops; j++)
	    	    if (s[j]->nucs[col].p_nuc == ltally->nucs[col].p_allele ||
			s[j]->nucs[col].q_nuc == ltally->nucs[col].p_allele)
			variable_pop = j;
	    } else if (p_cnt > 1 && q_cnt == 1) {
	    	for (int j = 0; j < this->num_pops; j++)
	    	    if (s[j]->nucs[col].p_nuc == ltally->nucs[col].q_allele ||
			s[j]->nucs[col].q_nuc == ltally->nucs[col].q_allele)
			variable_pop = j;
	    }
	    ltally->nucs[col].priv_allele = variable_pop;
	}
    }

    return 0;
}

template<class LocusT>
int PopSum<LocusT>::tally_ref_alleles(LocSum **s, int snp_index, 
				      short unsigned int &allele_cnt, 
				      char &p_allele, char &q_allele, 
				      short unsigned int &p_cnt, short unsigned int &q_cnt) 
{
    int  nucs[4] = {0};
    char nuc[2];

    p_allele   = 0;
    q_allele   = 0;
    allele_cnt = 0;

    for (int j = 0; j < this->num_pops; j++) {
	nuc[0] = 0;
	nuc[1] = 0;
        nuc[0] = s[j]->nucs[snp_index].p_nuc;
        nuc[1] = s[j]->nucs[snp_index].q_nuc;

	for (uint k = 0; k < 2; k++) 
	    switch(nuc[k]) {
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

    //
    // Determine how many alleles are present at this position in this population.
    // We cannot deal with more than two alternative alleles, if there are more than two
    // in a single population, print a warning and exclude this nucleotide position.
    //
    int i;
    for (i = 0; i < 4; i++)
	if (nucs[i] > 0) allele_cnt++;

    if (allele_cnt > 2) {
	p_allele = 0;
	q_allele = 0;
	return 0;
    }

    //
    // Record which nucleotide is the P allele and which is the Q allele.
    //
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

    //
    // Tabulate the number of populations the p_allele and the q_allele occur in.
    // 
    p_cnt = 0;
    q_cnt = 0;

    for (int j = 0; j < this->num_pops; j++) {
	nuc[0] = 0;
	nuc[1] = 0;
        nuc[0] = s[j]->nucs[snp_index].p_nuc;
        nuc[1] = s[j]->nucs[snp_index].q_nuc;

	for (uint k = 0; k < 2; k++) 
	    if (nuc[k] != 0 && nuc[k] == p_allele)
		p_cnt++;
	    else if (nuc[k] != 0 && nuc[k] == q_allele)
		q_cnt++;
    }

    return 1;
}

template<class LocusT>
PopPair *PopSum<LocusT>::Fst(int locus, int pop_1, int pop_2, int pos) 
{
    LocSum  *s_1  = this->pop(locus, pop_1);  /////// SLOW!
    LocSum  *s_2  = this->pop(locus, pop_2);
    PopPair *pair = new PopPair();

    //
    // If this locus only appears in one population do not calculate Fst.
    //
    if (s_1->nucs[pos].num_indv == 0 || s_2->nucs[pos].num_indv == 0) 
	return pair;

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
	return pair;

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

    if (allele_cnt > 2)
	return NULL;

    double tot_alleles = n_1 + n_2;
    double p_1 = round(n_1 * s_1->nucs[pos].p);
    double q_1 = n_1 - p_1;
    double p_2 = 
	s_1->nucs[pos].p_nuc == s_2->nucs[pos].p_nuc ? 
	s_2->nucs[pos].p : (1 - s_2->nucs[pos].p);
    p_2 = round(n_2 * p_2);
    double q_2 = n_2 - p_2;

    double pi_all = this->pi(tot_alleles, p_1 + p_2, q_1 + q_2);

    double bcoeff_1 = this->binomial_coeff(n_1, 2);
    double bcoeff_2 = this->binomial_coeff(n_2, 2);

    double num = (bcoeff_1 * pi_1) + (bcoeff_2 * pi_2);
    double den = pi_all * (bcoeff_1 + bcoeff_2);

    double Fst = 1 - (num / den);

    pair->alleles = tot_alleles;
    pair->fst     = Fst;
    pair->pi      = pi_all;

    this->fishers_exact_test(pair, p_1, q_1, p_2, q_2);

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

    //
    // Calculate Fst (corrected for different samples sizes) using an AMOVA method,
    // correcting for unequal sample sizes.
    // Derived from Weir, _Genetic Data Analysis II_, chapter 5, "F Statistics,", pp166-167.
    //
    double p_1_freq = s_1->nucs[pos].p;
    double q_1_freq = 1 - p_1_freq;
    double p_2_freq = 
	s_1->nucs[pos].p_nuc == s_2->nucs[pos].p_nuc ? 
	s_2->nucs[pos].p : (1 - s_2->nucs[pos].p);
    double q_2_freq = 1 - p_2_freq;

    double p_avg_cor = 
	( (s_1->nucs[pos].num_indv * p_1_freq) + (s_2->nucs[pos].num_indv * p_2_freq) ) / 
	( s_1->nucs[pos].num_indv + s_2->nucs[pos].num_indv );
    double n_avg_cor = 2 * ((s_1->nucs[pos].num_indv / 2) + (s_2->nucs[pos].num_indv / 2));

    pair->amova_fst =
	(
	 (s_1->nucs[pos].num_indv * pow((p_1_freq - p_avg_cor), 2) + 
	  s_2->nucs[pos].num_indv * pow((p_2_freq - p_avg_cor), 2))
	 / 
	 n_avg_cor 
	 )
	/ 
	(p_avg_cor * (1 - p_avg_cor));

    if (log_fst_comp) {
	pair->comp     = new double[18];
	pair->comp[0]  = n_1;
	pair->comp[1]  = n_2;
	pair->comp[2]  = tot_alleles;
	pair->comp[3]  = p_1;
	pair->comp[4]  = q_1;
	pair->comp[5]  = p_2;
	pair->comp[6]  = q_2;
	pair->comp[7]  = pi_1;
	pair->comp[8]  = pi_2;
	pair->comp[9]  = pi_all;
	pair->comp[10] = bcoeff_1;
	pair->comp[11] = bcoeff_2;
	pair->comp[12] = p_1_freq;
	pair->comp[13] = q_1_freq;
	pair->comp[14] = p_2_freq;
	pair->comp[15] = q_2_freq;
	pair->comp[16] = p_avg_cor;
	pair->comp[17] = n_avg_cor;
    }

    // //
    // // Calculate Fst using a pure parametric method (assumes allele counts are real, not 
    // // samples). Jakobsson, Edge, and Rosenberg. "The Relationship Between Fst and the 
    // // Frequency of the Most Frequent Allele." Genetics 193:515-528. Equation 4.
    // //
    // double sigma_1 = p_1_freq + q_1_freq;
    // double sigma_2 = p_2_freq + q_2_freq;
    // double delta_1 = fabs(p_1_freq - p_2_freq);
    // double delta_2 = fabs(q_1_freq - q_2_freq);

    // pair->jakob_fst = (pow(delta_1, 2) + pow(delta_2, 2)) / ( 4 - (pow(sigma_1, 2) + pow(sigma_2, 2)) );

    return pair;
}

template<class LocusT>
int PopSum<LocusT>::tally_fixed_pos(LocusT *locus, Datum **d, LocSum *s, int pos, uint start, uint end) 
{
    double num_indv = 0.0;
    char   p_nuc = 0;

    for (uint i = start; i <= end; i++) {
	if (d[i] == NULL || pos >= d[i]->len) continue;
	//
	// Before counting this individual, make sure the model definitively called this 
	// position as hEterozygous or hOmozygous.
	//
	if (d[i]->model[pos] == 'E') {
	    cerr << "Warning: heterozygous model call at fixed nucleotide position: " 
		 << "locus " << locus->id << " individual " << d[i]->id << "; position: " << pos << "\n";
	}
	num_indv++;
	p_nuc = locus->con[pos];
    }
    //
    // Record the results in the PopSum object.
    //
    s->nucs[pos].loc_id   = locus->id;
    s->nucs[pos].bp       = locus->sort_bp(pos);
    s->nucs[pos].fixed    = true;
    s->nucs[pos].num_indv = num_indv;
    s->nucs[pos].alleles  = 2 * num_indv;

    if (num_indv > 0) {
	s->nucs[pos].p        =  1.0;
	s->nucs[pos].p_nuc    =  p_nuc;
	s->nucs[pos].obs_hom  =  1.0;
	s->nucs[pos].obs_het  =  0.0;
	s->nucs[pos].exp_hom  =  1.0;
	s->nucs[pos].exp_het  =  0.0;
	s->nucs[pos].stat[0]  =  0.0; // pi
	s->nucs[pos].stat[1]  = -7.0; // fis
    }

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
	if (d[i] == NULL || pos >= d[i]->len || d[i]->model[pos] == 'U') continue;

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

    if (allele_cnt > 2)
	return -1;

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
    // Calculate observed genotype frequencies.
    //
    double num_indv = 0.0;
    double obs_het  = 0.0;
    double obs_p    = 0.0;
    double obs_q    = 0.0;

    for (i = start; i <= end; i++) {
	if (d[i] == NULL || pos >= d[i]->len) continue;
	//
	// Before counting this individual, make sure the model definitively called this 
	// position as hEterozygous or hOmozygous.
	//
	if (d[i]->model[pos] == 'E' || d[i]->model[pos] == 'O')
	    num_indv++;
	else
	    continue;

	if (d[i]->obshap.size() > 1 &&
	    this->tally_observed_haplotypes(d[i]->obshap, snp_index) == 2)
		obs_het++;
	else if (d[i]->obshap[0][snp_index] == p_allele)
	    obs_p++;
	else if (d[i]->obshap[0][snp_index] == q_allele)
	    obs_q++;
    }
    //cerr << "  Num Individuals: " << num_indv << "; Obs Hets: " << obs_het << "; Obs P: " << obs_p << "; Obs Q: " << obs_q << "\n";

    if (num_indv == 0) return 0;

    //
    // Calculate total number of alleles
    //
    double tot_alleles = num_indv * 2;
    double allele_p    = obs_het + (2 * obs_p);
    double allele_q    = obs_het + (2 * obs_q);

    //
    // Calculate Pi, equivalent to expected heterozygosity (exp_het)
    //
    s->nucs[pos].stat[0] = this->pi(tot_alleles, allele_p, allele_q);

    if (s->nucs[pos].stat[0] == 0.0)
	s->nucs[pos].fixed = true;

    //
    // Convert to allele frequencies
    //
    allele_p = allele_p / tot_alleles;
    allele_q = allele_q / tot_alleles;

    //cerr << "  P allele frequency: " << allele_p << "; Q allele frequency: " << allele_q << "\n";

    // //
    // // If the minor allele frequency is below the cutoff, set it to zero.
    // //
    // if (minor_allele_freq > 0) {
    // 	if (allele_p < allele_q) {
    // 	    if (allele_p < minor_allele_freq) {
    // 		s->nucs[pos].pi            = 0.0;
    // 		s->nucs[pos].fixed         = true;
    // 		s->nucs[pos].filtered_site = true;
    // 		return 0;
    // 	    }
    // 	} else {
    // 	    if (allele_q < minor_allele_freq) {
    // 		s->nucs[pos].pi            = 0.0;
    // 		s->nucs[pos].fixed         = true;
    // 		s->nucs[pos].filtered_site = true;
    // 		return 0;
    // 	    }
    // 	}
    // }

    //
    // Calculate expected genotype frequencies.
    //
    double exp_het = 2 * allele_p * allele_q; // 2pq
    // double exp_p   = allele_p * allele_p;     // p^2
    // double exp_q   = allele_q * allele_q;     // q^2

    //cerr << "  Expected Het: " << exp_het << "; Expected P: " << exp_p << "; Expected Q: " << exp_q << "\n";

    obs_het = obs_het / num_indv;
    obs_p   = obs_p   / num_indv;
    obs_q   = obs_q   / num_indv;

    //cerr << "  Obs Hets Freq: " << obs_het << "; Obs P Freq: " << obs_p << "; Obs Q Freq: " << obs_q << "\n";

    //
    // Record the results in the PopSum object.
    //
    s->nucs[pos].loc_id   = locus->id;
    s->nucs[pos].bp       = locus->sort_bp(pos);
    s->nucs[pos].num_indv = num_indv;
    s->nucs[pos].alleles  = tot_alleles;
    s->nucs[pos].p        = allele_p > allele_q ? allele_p : allele_q;
    s->nucs[pos].p_nuc    = allele_p > allele_q ? p_allele : q_allele;
    s->nucs[pos].q_nuc    = allele_p > allele_q ? q_allele : p_allele;
    s->nucs[pos].obs_hom  = 1 - obs_het;
    s->nucs[pos].obs_het  = obs_het;
    s->nucs[pos].exp_hom  = 1 - exp_het;
    s->nucs[pos].exp_het  = exp_het;

    //
    // Calculate F_is, the inbreeding coefficient of an individual (I) relative to the subpopulation (S):
    //   Fis = (exp_het - obs_het) / exp_het
    //
    double fis = s->nucs[pos].pi == 0 ? -7 : (s->nucs[pos].pi - obs_het) / s->nucs[pos].pi;

    s->nucs[pos].stat[1] = fis;

    return 0;
}

template<class LocusT>
int PopSum<LocusT>::tally_observed_haplotypes(vector<char *> &obshap, int snp_index) 
{
    int  nucs[4] = {0};
    char nuc;

    //
    // Pull each allele for this SNP from the observed haplotype.
    //
    for (uint j = 0; j < obshap.size(); j++) {
        nuc = obshap[j][snp_index];

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

    int allele_cnt = 0;
    for (int i = 0; i < 4; i++)
	if (nucs[i] > 0) allele_cnt++;

    return allele_cnt;
}

template<class LocusT>
int PopSum<LocusT>::fishers_exact_test(PopPair *pair, double p_1, double q_1, double p_2, double q_2)
{
    //                            | Allele1 | Allele2 | 
    // Fisher's Exact Test:  -----+---------+---------+
    //                       Pop1 |   p_1   |   q_1   |
    //                       Pop2 |   p_2   |   q_2   |
    // Probability p:
    // p = ((p_1 + q_1)!(p_2 + q_2)!(p_1 + p_2)!(q_1 + q_2)!) / (n! p_1! q_1! p_2! q_2!)
    //
    // According to:
    //   Jerrold H. Zar, "A fast and efficient algorithm for the Fisher exact test."
    //   Behavior Research Methods, Instruments, & Computers 1987, 19(4): 43-44
    //
    // Probability p can be calculated as three binomial coefficients:
    //   Let p_1 + q_1 = r_1; p_2 + q_2 = r_2; p_1 + p_2 = c_1; q_1 + q_2 = c_2
    //
    //   p = (r_1 choose p_1)(r_2 choose p_2) / (n choose c_1)
    //
    // Fisher's Exact test algorithm implemented according to Sokal and Rohlf, _Biometry_, section 17.4.
    //

    double r_1  = p_1 + q_1;
    double r_2  = p_2 + q_2;
    double c_1  = p_1 + p_2;
    double d_1  = p_1 * q_2;
    double d_2  = p_2 * q_1;
    double n    = r_1 + r_2;
    double p    = 0.0;

    // char p1_str[32], p2_str[32], q1_str[32], q2_str[32];
    // sprintf(p1_str, "% 3.0f", p_1);
    // sprintf(q1_str, "% 3.0f", q_1);
    // sprintf(p2_str, "% 3.0f", p_2);
    // sprintf(q2_str, "% 3.0f", q_2);
    //
    // cerr 
    // 	<< "     | Allele1 | Allele2 | " << "\n"
    // 	<< "-----+---------+---------+" << "\n"
    // 	<< "Pop1 |   " << p1_str << "   |   " << q1_str << "   |" << "\n"
    // 	<< "Pop2 |   " << p2_str << "   |   " << q2_str << "   |" << "\n\n";

    //
    // Compute the first tail.
    //
    double p1     = p_1;
    double q1     = q_1;
    double p2     = p_2;
    double q2     = q_2;
    double tail_1 = 0.0;
    double den    = this->binomial_coeff(n, c_1);

    //
    // If (p_1*q_2 - p_2*q_1) < 0 decrease cells p_1 and q_2 by one and add one to p_2 and q_1. 
    // Compute p and repeat until one or more cells equal 0.
    //
    if (d_1 - d_2 < 0) {
	do {
	    p = (this->binomial_coeff(r_1, p1) * this->binomial_coeff(r_2, p2)) / den;

	    tail_1 += p;
	    p1--;
	    q2--;
	    p2++;
	    q1++;
	} while (p1 >= 0 && q2 >= 0);

    } else {
	//
	// Else, if (p_1*q_2 - p_2*q_1) > 0 decrease cells p_2 and q_1 by one and add one to p_1 and q_2. 
	// Compute p and repeat until one or more cells equal 0.
	//
	do {
	    p = (this->binomial_coeff(r_1, p1) * this->binomial_coeff(r_2, p2)) / den;

	    tail_1 += p;

	    p2--;
	    q1--;
	    p1++;
	    q2++;
	} while (p2 >= 0 && q1 >= 0);
    }

    //
    // Compute the second tail.
    //
    double tail_2 = 0.0;
    p = 0;

    //
    // If (p_1*q_2 - p_2*q_1) < 0, set to zero the smaller of the two frequencies, adjusting the other values 
    // to keep the marginals the same.
    //
    if (d_1 - d_2 < 0) {
	if (p2 < q1) {
	    q2 += p2;
	    p1 += p2;
	    q1 -= p2;
	    p2  = 0;
	} else {
	    p1 += q1;
	    q2 += q1;
	    p2 -= q1;
	    q1  = 0;
	}
    } else {
	if (p1 < q2) {
	    q1 += p1;
	    p2 += p1;
	    q2 -= p1;
	    p1  = 0;
	} else {
	    p2 += q2;
	    q1 += q2;
	    p1 -= q2;
	    q2  = 0;
	}
    }

    //
    // If (p_1*q_2 - p_2*q_1) < 0 decrease cells p_1 and q_2 by one and add one to p_2 and q_1. 
    // Compute p and repeat until tail_2 > tail_1.
    //
    if (d_1 - d_2 < 0) {
	do {
	    p = (this->binomial_coeff(r_1, p1) * this->binomial_coeff(r_2, p2)) / den;

	    tail_2 += p;

	    p1--;
	    q2--;
	    p2++;
	    q1++;
	} while (tail_2 < tail_1 && p1 >= 0 && q2 >= 0);

	tail_2 -= p;

    } else {
	//
	// Else, if (p_1*q_2 - p_2*q_1) > 0 decrease cells p_2 and q_1 by one and add one to p_1 and q_2. 
	// Compute p and repeat until one or more cells equal 0.
	//
	do {
	    p = (this->binomial_coeff(r_1, p1) * this->binomial_coeff(r_2, p2)) / den;

	    tail_2 += p;

	    p2--;
	    q1--;
	    p1++;
	    q2++;
	} while (tail_2 < tail_1 && p2 >= 0 && q1 >= 0);

	tail_2 -= p;
    }

    pair->fet_p = tail_1 + tail_2;

    if (pair->fet_p > 1.0) pair->fet_p = 1.0;

    //
    // Calculate the odds ratio. To account for possible cases were one allele frequency is
    // zero, we will increment all allele frequencies by one.
    //
    if (p_1 == 0 || q_1 == 0 || p_2 == 0 || q_2 == 0) {
	p_1++;
	q_1++;
	p_2++;
	q_2++;
    }
    pair->fet_or = (p_1 * q_2) / (q_1 * p_2);

    double ln_fet_or = pair->fet_or > 0 ? log(pair->fet_or) : 0.0;

    //
    // Calculate the standard error of the natural log of the odds ratio
    //
    double se = pair->fet_or > 0 ? sqrt((1 / p_1) + (1 / q_1) + (1 / p_2) + (1 / q_2)) : 0.0;

    //
    // Calculate the confidence intervals of the natural log of the odds ratio.
    //
    double ln_ci_low  = pair->fet_or > 0 ? ln_fet_or - (1.96 * se) : 0;
    double ln_ci_high = pair->fet_or > 0 ? ln_fet_or + (1.96 * se) : 0;

    //
    // Convert the confidence intervals out of natural log space
    //
    pair->ci_low  = pair->fet_or > 0 ? exp(ln_ci_low)  : 0;
    pair->ci_high = pair->fet_or > 0 ? exp(ln_ci_high) : 0;
    pair->lod     = fabs(log10(pair->fet_or));

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

template<class LocusT>
LocTally *PopSum<LocusT>::locus_tally(int locus) 
{
    return this->loc_tally[this->locus_order[locus]];
}

#endif // __POPSUM_H__
