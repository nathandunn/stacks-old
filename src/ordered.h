// -*-mode:c++; c-style:k&r; c-basic-offset:4;-*-
//
// Copyright 2014, Julian Catchen <jcatchen@uoregon.edu>
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

#ifndef __ORDERED_H__
#define __ORDERED_H__

#include <fstream>
using std::ifstream;
using std::ofstream;
using std::cin;
using std::cout;
using std::cerr;
#include <vector>
using std::vector;
#include <map>
using std::map;
#include <set>
using std::set;

#include "PopSum.h"

enum loc_type {haplotype, snp};

template<class StatT>
class Ordered {
public:
    ofstream        *log_fh;
    PopSum<CSLocus> *psum;

    int incompatible_loci;
    int multiple_loci;

    Ordered()  { }
    virtual ~Ordered() { }

    int init_sites(vector<StatT *> &, map<uint, uint> &, vector<CSLocus *> &);
    int init_sites(vector<StatT *> &, map<uint, uint> &, vector<CSLocus *> &, int);
    int init_sites(vector<StatT *> &, map<uint, uint> &, vector<CSLocus *> &, int, int);
    int init_haplotypes(vector<StatT *> &, map<uint, uint> &, vector<CSLocus *> &);
};

template<class StatT>
int 
Ordered<StatT>::init_sites(vector<StatT *> &sites, map<uint, uint> &sites_key, vector<CSLocus *> &sorted_loci)
{
    CSLocus   *loc;
    LocTally  *ltally;
    int        len;
    set<int>   bps;

    //
    // We need to create an array to store all the SNPs for exporting. We must
    // account for positions in the genome that are covered by more than one RAD tag.
    //
    for (uint pos = 0; pos < sorted_loci.size(); pos++) {
	loc    = sorted_loci[pos];
	len    = strlen(loc->con);
	ltally = this->psum->locus_tally(loc->id);

	for (int k = 0; k < len; k++) {
	    if (ltally->nucs[k].allele_cnt == 2)
		bps.insert(ltally->nucs[k].bp);
	}
    }

    sites.resize(bps.size(), NULL);

    //
    // Create a key describing where in the sites array to find each basepair coordinate.
    //
    set<int>::iterator it;
    int i = 0;
    for (it = bps.begin(); it != bps.end(); it++) {
	sites_key[*it] = i;
	i++;
    }

    return 0;
}

template<class StatT>
int 
Ordered<StatT>::init_sites(vector<StatT *> &sites, map<uint, uint> &sites_key, vector<CSLocus *> &sorted_loci, int pop_id)
{
    CSLocus *loc;
    LocSum  *lsum;
    int      len;
    set<int> bps;

    //
    // We need to create an array to store all the summary statistics for smoothing. We must
    // account for positions in the genome that are covered by more than one RAD tag.
    //
    for (uint pos = 0; pos < sorted_loci.size(); pos++) {
	loc  = sorted_loci[pos];
	len  = strlen(loc->con);
	lsum = this->psum->pop(loc->id, pop_id);

	for (int k = 0; k < len; k++) {
	    if (lsum->nucs[k].num_indv > 0)
		bps.insert(lsum->nucs[k].bp);
	}
    }

    sites.resize(bps.size(), NULL);

    //
    // Create a key describing where in the sites array to find each basepair coordinate.
    //
    set<int>::iterator it;
    int i = 0;
    for (it = bps.begin(); it != bps.end(); it++) {
	sites_key[*it] = i;
	i++;
    }

    return 0;
}

template<class StatT>
int 
Ordered<StatT>::init_sites(vector<StatT *> &sites, map<uint, uint> &sites_key, vector<CSLocus *> &sorted_loci, int pop_id_1, int pop_id_2)
{
    CSLocus *loc;
    LocSum  *lsum_1, *lsum_2;
    int      len;
    set<int> bps;

    //
    // We need to create an array to store all the pair values for computing smoothed Fst. We must
    // account for positions in the genome that are covered by more than one RAD tag.
    //
    for (uint pos = 0; pos < sorted_loci.size(); pos++) {
	loc    = sorted_loci[pos];
	len    = strlen(loc->con);
	lsum_1 = this->psum->pop(loc->id, pop_id_1);
	lsum_2 = this->psum->pop(loc->id, pop_id_2);

	for (int k = 0; k < len; k++) {
	    if (lsum_1->nucs[k].num_indv > 0 &&
		lsum_2->nucs[k].num_indv > 0)
		bps.insert(lsum_1->nucs[k].bp); // slow
	}
    }

    sites.resize(bps.size(), NULL);

    //
    // Create a key describing where in the sites array to find each basepair coordinate.
    //
    set<int>::iterator it;
    int i = 0;
    for (it = bps.begin(); it != bps.end(); it++) {
	sites_key[*it] = i; // slow
	i++;
    }

    return 0;
}

template<class StatT>
int 
Ordered<StatT>::init_haplotypes(vector<StatT *> &sites, map<uint, uint> &sites_key, vector<CSLocus *> &sorted_loci) 
{
    CSLocus *loc;
    int      bp;
    set<int> bps;

    for (uint pos = 0; pos < sorted_loci.size(); pos++) {
	loc = sorted_loci[pos];
	bp  = loc->sort_bp();

	bps.insert(bp);
    }

    sites.resize(bps.size(), NULL);

    //
    // Create a key describing where in the sites array to find each basepair coordinate.
    //
    set<int>::iterator it;
    int i = 0;
    for (it = bps.begin(); it != bps.end(); it++) {
	sites_key[*it] = i;
	i++;
    }

    return 0;
}

template<class StatT>
class OHaplotypes: public Ordered<StatT> {
public:
    OHaplotypes(): Ordered<StatT>() { }

    int order(vector<StatT *> &, map<uint, uint> &, vector<CSLocus *> &);
};

template<class StatT>
int 
OHaplotypes<StatT>::order(vector<StatT *> &sites, map<uint, uint> &sites_key, vector<CSLocus *> &sorted_loci) 
{
    this->init_haplotypes(sites, sites_key, sorted_loci);

    return 0;
};

template<class StatT>
class OPopPair: public Ordered<StatT> {
public:
    OPopPair(PopSum<CSLocus> *psum, ofstream &log_fh): Ordered<StatT>() { 
	this->log_fh = &log_fh;
	this->psum   = psum;
    }

    int order(vector<StatT *> &, map<uint, uint> &, vector<CSLocus *> &, int, int);
};

template<class StatT>
int 
OPopPair<StatT>::order(vector<StatT *> &sites, map<uint, uint> &sites_key, vector<CSLocus *> &sorted_loci, int pop_1, int pop_2) 
{
    CSLocus *loc;
    StatT   *pair;
    int      len;

    this->incompatible_loci = 0;
    this->multiple_loci     = 0;

    this->init_sites(sites, sites_key, sorted_loci, pop_1, pop_2);

    for (uint pos = 0; pos < sorted_loci.size(); pos++) {
	loc = sorted_loci[pos];
	len = strlen(loc->con);

	for (int k = 0; k < len; k++) {

	    pair = this->psum->Fst(loc->id, pop_1, pop_2, k);

	    //
	    // Locus is incompatible, log this position.
	    //
	    if (pair == NULL) {
		this->incompatible_loci++;
		*(this->log_fh) << "between_population\t"
				<< "incompatible_locus\t"
				<< loc->id << "\t"
				<< loc->loc.chr << "\t"
				<< loc->sort_bp(k) << "\t"
				<< k << "\t" 
				<< pop_1 << "\t" 
				<< pop_2 << "\n";
		delete pair;
		continue;
	    }

	    pair->loc_id = loc->id;
	    pair->bp     = loc->sort_bp(k);
	    pair->col    = k;

	    //
	    // Locus is fixed in both populations, or was only found in one population.
	    //
	    if (pair->pi == 0) {
		delete pair;
		continue;
	    }

	    //
	    // Check if this basepair position is already covered by a RAD site.
	    //
	    if (sites[sites_key[pair->bp]] != NULL) {
		this->multiple_loci++;
		*(this->log_fh) << "between_population\t"
				<< "multiple_locus\t"
				<< loc->id << "\t"
				<< loc->loc.chr << "\t"
				<< pair->bp << "\t"
				<< k << "\t" 
				<< pop_key[pop_1] << "\t" 
				<< pop_key[pop_2] << "\n";
		delete pair;
		continue;
	    }

	    sites[sites_key[pair->bp]] = pair;
	}
    }

    return 0;
};

template<class StatT>
class OSumStat: public Ordered<StatT> {
public:
    OSumStat(PopSum<CSLocus> *psum, ofstream &log_fh): Ordered<StatT>() { 
	this->log_fh = &log_fh;
	this->psum   = psum;
    }

    int order(vector<StatT *> &, vector<CSLocus *> &, int);
};

template<class StatT>
int 
OSumStat<StatT>::order(vector<StatT *> &sites, vector<CSLocus *> &sorted_loci, int pop_id) 
{
    this->incompatible_loci = 0;
    this->multiple_loci     = 0;

    map<uint, uint> sites_key;

    this->init_sites(sites, sites_key, sorted_loci, pop_id);

    CSLocus *loc;
    LocSum  *lsum;
    int      len;

    //
    // Assign nucleotides to their proper, ordered location in the genome, 
    // checking that a site hasn't already been covered by another RAD locus.
    //
    for (uint pos = 0; pos < sorted_loci.size(); pos++) {
	loc  = sorted_loci[pos];
	len  = strlen(loc->con);
	lsum = this->psum->pop(loc->id, pop_id);

	for (int k = 0; k < len; k++) {
	    if (lsum->nucs[k].num_indv == 0) continue;

	    if (sites_key.count(lsum->nucs[k].bp) == 0) {
		cerr << "Error: locus " << lsum->nucs[k].loc_id << " at " << lsum->nucs[k].bp << "bp is not defined in the sites map.\n";

	    } else if (sites[sites_key[lsum->nucs[k].bp]] == NULL) {
		sites[sites_key[lsum->nucs[k].bp]] = &(lsum->nucs[k]);

	    } else {
		this->multiple_loci++;
		*(this->log_fh) << "within_population\t"
				<< "multiple_locus\t"
				<< loc->id << "\t"
				<< loc->loc.chr << "\t"
				<< lsum->nucs[k].bp << "\t"
				<< k << "\t" 
				<< pop_id << "\t"
				<< "conflicts with locus " << sites[sites_key[lsum->nucs[k].bp]]->loc_id << "\n";
	    }
	}
    }

    return 0;
};

template<class StatT>
class OLocTally: public Ordered<StatT> {
public:
    OLocTally(PopSum<CSLocus> *psum, ofstream &log_fh): Ordered<StatT>() { 
	this->log_fh = &log_fh;
	this->psum   = psum;
    }

    int order(vector<StatT *> &, vector<CSLocus *> &);
};

template<class StatT>
int 
OLocTally<StatT>::order(vector<StatT *> &sites, vector<CSLocus *> &sorted_loci) 
{
    this->incompatible_loci = 0;
    this->multiple_loci     = 0;

    map<uint, uint> sites_key;

    this->init_sites(sites, sites_key, sorted_loci);

    CSLocus   *loc;
    LocTally  *ltally;
    int        len;

    //
    // Assign nucleotides to their proper, ordered location in the genome, 
    // checking that a site hasn't already been covered by another RAD locus.
    //
    for (uint pos = 0; pos < sorted_loci.size(); pos++) {
	loc    = sorted_loci[pos];
	len    = strlen(loc->con);
	ltally = this->psum->locus_tally(loc->id);

	for (int k = 0; k < len; k++) {
	    if (ltally->nucs[k].allele_cnt != 2) continue;

	    if (sites_key.count(ltally->nucs[k].bp) == 0) {
		cerr << "Error: locus " << ltally->nucs[k].loc_id << " at " << ltally->nucs[k].bp << "bp is not defined in the sites map.\n";

	    } else if (sites[sites_key[ltally->nucs[k].bp]] == NULL) {
		sites[sites_key[ltally->nucs[k].bp]] = &(ltally->nucs[k]);

	    } else {
		this->multiple_loci++;
		*(this->log_fh) << "within_population\t"
				<< "multiple_locus\t"
				<< loc->id << "\t"
				<< loc->loc.chr << "\t"
				<< ltally->nucs[k].bp << "\t"
				<< k << "\t" 
				<< "conflicts with locus " << sites[sites_key[ltally->nucs[k].bp]]->loc_id << "\n";
	    }
	}
    }

    return 0;
};

#endif // __ORDERED_H__
