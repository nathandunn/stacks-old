// -*-mode:c++; c-style:k&r; c-basic-offset:4;-*-
//
// Copyright 2014-2017, Julian Catchen <jcatchen@illinois.edu>
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

#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include <set>

#include "MetaPopInfo.h"
#include "PopSum.h"

extern MetaPopInfo mpopi;

enum loc_type {haplotype, snp};

template<class StatT>
class Ordered {
public:
    ofstream *log_fh;

    int incompatible_loci;
    int multiple_loci;

    Ordered()  {}
    virtual ~Ordered() {}

    int init_sites(vector<const StatT *> &, map<uint, uint> &, const vector<LocBin *> &);
    int init_sites(vector<const StatT *> &, map<uint, uint> &, const vector<LocBin *> &, uint);
    int init_sites(vector<const StatT *> &, map<uint, uint> &, const vector<LocBin *> &, uint, uint);
    int init_haplotypes(vector<const StatT *> &, map<uint, uint> &, const vector<LocBin *> &);
    int init_haplotypes(vector<const StatT *> &, map<uint, uint> &, const vector<LocBin *> &, uint, uint);
};

template<class StatT>
int
Ordered<StatT>::init_sites(vector<const StatT *> &sites, map<uint, uint> &sites_key, const vector<LocBin *> &sorted_loci)
{
    const CSLocus  *loc;
    const LocTally *ltally;
    int      len;
    set<int> bps;

    //
    // We need to create an array to store all the SNPs for exporting. We must
    // account for positions in the genome that are covered by more than one RAD tag.
    //
    for (uint pos = 0; pos < sorted_loci.size(); pos++) {
        loc    = sorted_loci[pos]->cloc;
        len    = loc->len;
        ltally = sorted_loci[pos]->s->meta_pop();

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
Ordered<StatT>::init_sites(vector<const StatT *> &sites, map<uint, uint> &sites_key, const vector<LocBin *> &sorted_loci, uint pop_id)
{
    const CSLocus *loc;
    const LocSum  *lsum;
    int      len;
    set<int> bps;

    //
    // We need to create an array to store all the summary statistics for smoothing. We must
    // account for positions in the genome that are covered by more than one RAD tag.
    //
    for (uint pos = 0; pos < sorted_loci.size(); pos++) {
        loc  = sorted_loci[pos]->cloc;
        len  = strlen(loc->con);
        lsum = sorted_loci[pos]->s->per_pop(pop_id);

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
Ordered<StatT>::init_sites(vector<const StatT *> &sites, map<uint, uint> &sites_key, const vector<LocBin *> &sorted_loci, uint pop_id_1, uint pop_id_2)
{
    const CSLocus *loc;
    const LocSum  *lsum_1, *lsum_2;
    int      len;
    set<int> bps;

    //
    // We need to create an array to store all the pair values for computing smoothed Fst. We must
    // account for positions in the genome that are covered by more than one RAD tag.
    //
    for (uint pos = 0; pos < sorted_loci.size(); pos++) {
        loc    = sorted_loci[pos]->cloc;
        len    = strlen(loc->con);
        lsum_1 = sorted_loci[pos]->s->per_pop(pop_id_1);
        lsum_2 = sorted_loci[pos]->s->per_pop(pop_id_2);

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
Ordered<StatT>::init_haplotypes(vector<const StatT *> &sites, map<uint, uint> &sites_key, const vector<LocBin *> &sorted_loci)
{
    const CSLocus *loc;
    int      bp;
    set<int> bps;

    for (uint pos = 0; pos < sorted_loci.size(); pos++) {
        loc = sorted_loci[pos]->cloc;
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

    int order(vector<const StatT *> &, map<uint, uint> &, const vector<LocBin *> &);
    int order(vector<const StatT *> &, map<uint, uint> &, const vector<LocBin *> &, const vector<StatT *> &);
};

template<class StatT>
int
OHaplotypes<StatT>::order(vector<const StatT *> &sites, map<uint, uint> &sites_key, const vector<LocBin *> &sorted_loci)
{
    this->init_haplotypes(sites, sites_key, sorted_loci);

    return 0;
};

template<class StatT>
int
OHaplotypes<StatT>::order(vector<const StatT *> &sites, map<uint, uint> &sites_key, const vector<LocBin *> &sorted_loci, const vector<StatT *> &div)
{
    StatT *pair;

    this->init_haplotypes(sites, sites_key, sorted_loci);

    for (uint i = 0; i < div.size(); i++) {
        pair = div[i];

        if (pair == NULL)
            continue;
        
        sites[sites_key[pair->bp]] = pair;
    }

    return 0;
};

template<class StatT>
class OPopPair: public Ordered<StatT> {
public:
    OPopPair(ofstream &log_fh): Ordered<StatT>() {
        this->log_fh = &log_fh;
    }

    int order(vector<const StatT *> &, map<uint, uint> &, const vector<LocBin *> &, const vector<StatT **> &);
};

template<class StatT>
int
OPopPair<StatT>::order(vector<const StatT *> &sites, map<uint, uint> &sites_key, const vector<LocBin *> &sorted_loci, const vector<StatT **> &div)
{
    CSLocus *loc;
    StatT  **pair;
    uint     cloc_len;
    bool     found = false;

    this->incompatible_loci = 0;
    this->multiple_loci     = 0;

    uint pop_1=UINT_MAX, pop_2=UINT_MAX;
    for (uint i = 0; i < div.size(); i++) {
        cloc_len = strlen(sorted_loci[i]->cloc->con);
        pair     = div[i];
        
        for (uint j = 0; j < cloc_len; j++) {
            if (pair[j] != NULL) {
                pop_1 = pair[j]->pop_1;
                pop_2 = pair[j]->pop_2;
                found = true;
                break;
            }
        }
        if (found == true) break;
    }
    assert(found);

    this->init_sites(sites, sites_key, sorted_loci, pop_1, pop_2);

    for (uint i = 0; i < div.size(); i++) {
        loc      = sorted_loci[i]->cloc;
        cloc_len = strlen(loc->con);
        pair     = div[i];

        for (uint pos = 0; pos < cloc_len; pos++) {
            if (pair[pos] == NULL)
                continue;

            //
            // Check if this basepair position is already covered by a RAD site.
            //
            if (sites[sites_key[pair[pos]->bp]] != NULL) {
                this->multiple_loci++;
                *(this->log_fh) << "between_population\t"
                                << "multiple_locus\t"
                                << loc->id << "\t"
                                << loc->loc.chr() << "\t"
                                << pair[pos]->bp + 1   << "\t"
                                << pos << "\t"
                                << mpopi.pops()[pair[pos]->pop_1].name << "\t"
                                << mpopi.pops()[pair[pos]->pop_2].name << "\n";
                continue;
            }

            sites[sites_key[pair[pos]->bp]] = pair[pos];
        }
    }

    return 0;
};

template<class StatT>
class OSumStat: public Ordered<StatT> {
public:
    OSumStat(ofstream &log_fh): Ordered<StatT>() {
        this->log_fh = &log_fh;
    }

    int order(vector<const StatT *> &, const vector<LocBin *> &, uint);
};

template<class StatT>
int
OSumStat<StatT>::order(vector<const StatT *> &sites, const vector<LocBin *> &sorted_loci, uint pop_id)
{
    this->incompatible_loci = 0;
    this->multiple_loci     = 0;

    map<uint, uint> sites_key;

    this->init_sites(sites, sites_key, sorted_loci, pop_id);

    const CSLocus *loc;
    const LocSum  *lsum;
    int            len;

    //
    // Assign nucleotides to their proper, ordered location in the genome,
    // checking that a site hasn't already been covered by another RAD locus.
    //
    for (uint pos = 0; pos < sorted_loci.size(); pos++) {
        loc  = sorted_loci[pos]->cloc;
        len  = strlen(loc->con);
        lsum = sorted_loci[pos]->s->per_pop(pop_id);

        for (int k = 0; k < len; k++) {
            if (lsum->nucs[k].num_indv == 0) continue;

            if (sites_key.count(lsum->nucs[k].bp) == 0) {
                cerr << "Error: locus " << lsum->nucs[k].loc_id << " at " << lsum->nucs[k].bp +1 << "bp is not defined in the sites map.\n";

            } else if (sites[sites_key[lsum->nucs[k].bp]] == NULL) {
                sites[sites_key[lsum->nucs[k].bp]] = &(lsum->nucs[k]);

            } else {
                this->multiple_loci++;
                *(this->log_fh) << "within_population\t"
                                << "multiple_locus\t"
                                << loc->id << "\t"
                                << loc->loc.chr() << "\t"
                                << lsum->nucs[k].bp +1 << "\t"
                                << k << "\t"
                                << mpopi.pops()[pop_id].name << "\t"
                                << "conflicts with locus " << sites[sites_key[lsum->nucs[k].bp]]->loc_id << "\n";
            }
        }
    }

    return 0;
};

template<class StatT>
class OLocTally: public Ordered<StatT> {
public:
    OLocTally(ofstream &log_fh): Ordered<StatT>() {
        this->log_fh = &log_fh;
    }
    OLocTally(): Ordered<StatT>() {
        this->log_fh = NULL;
    }

    int order(vector<const StatT *> &, const vector<LocBin *> &);
};

template<class StatT>
int
OLocTally<StatT>::order(vector<const StatT *> &sites, const vector<LocBin *> &sorted_loci)
{
    this->incompatible_loci = 0;
    this->multiple_loci     = 0;

    map<uint, uint> sites_key;

    this->init_sites(sites, sites_key, sorted_loci);

    const CSLocus  *loc;
    const LocTally *ltally;
    int len;

    //
    // Assign nucleotides to their proper, ordered location in the genome,
    // checking that a site hasn't already been covered by another RAD locus.
    //
    for (uint pos = 0; pos < sorted_loci.size(); pos++) {
        loc    = sorted_loci[pos]->cloc;
        len    = loc->len;
        ltally = sorted_loci[pos]->s->meta_pop();

        for (int k = 0; k < len; k++) {
            if (ltally->nucs[k].allele_cnt != 2) continue;

            if (sites_key.count(ltally->nucs[k].bp) == 0) {
                cerr << "Error: locus " << ltally->nucs[k].loc_id << " at " << ltally->nucs[k].bp +1 << "bp is not defined in the sites map.\n";

            } else if (sites[sites_key[ltally->nucs[k].bp]] == NULL) {
                sites[sites_key[ltally->nucs[k].bp]] = &(ltally->nucs[k]);

            } else {
                this->multiple_loci++;
                if (this->log_fh != NULL)
                    *(this->log_fh) << "within_population\t"
                                    << "multiple_locus\t"
                                    << loc->id << "\t"
                                    << loc->loc.chr() << "\t"
                                    << ltally->nucs[k].bp +1 << "\t"
                                    << k << "\t"
                                    << "conflicts with locus " << sites[sites_key[ltally->nucs[k].bp]]->loc_id << "\n";
            }
        }
    }

    return 0;
};

#endif // __ORDERED_H__
