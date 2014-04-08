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
    ofstream        log_fh;
    loc_type        type;
    map<uint, uint> sites_key;

public:
    Ordered(loc_type type)  { 
	this->type   = type;
    }
    virtual ~Ordered() { 
    }

    int init_sites(vector<StatT> &, int, vector<CSLocus *> &, PopSum<CSLocus> *);
    int init_haplotypes(vector<StatT> &, int, vector<CSLocus *> &, PopSum<CSLocus> *);
};

template<class StatT>
int 
Ordered<StatT>::init_sites(vector<StatT> &sites, int pop_id, 
			   vector<CSLocus *> &sorted_loci, PopSum<CSLocus> *psum) 
{
    CSLocus *loc;
    LocSum  *lsum;
    int      len;
    set<int> bps;

    for (uint pos = 0; pos < sorted_loci.size(); pos++) {
	loc  = sorted_loci[pos];
	len  = strlen(loc->con);
	lsum = psum->pop(loc->id, pop_id);

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
	this->sites_key[*it] = i;
	i++;
    }

    return 0;
}

template<class StatT>
int 
Ordered<StatT>::init_haplotypes(vector<StatT> &sites, int pop_id, 
				vector<CSLocus *> &sorted_loci, PopSum<CSLocus> *psum) 
{
    CSLocus *loc;
    LocSum  *lsum;
    set<int> bps;

    for (uint pos = 0; pos < sorted_loci.size(); pos++) {
	loc      = sorted_loci[pos];
	lsum     = psum->pop(loc->id, pop_id);
	lsum->bp = loc->sort_bp();

	if (lsum->alleles > 0) 
	    bps.insert(lsum->bp);
    }

    sites.resize(bps.size(), NULL);

    //
    // Create a key describing where in the sites array to find each basepair coordinate.
    //
    set<int>::iterator it;
    int i = 0;
    for (it = bps.begin(); it != bps.end(); it++) {
	this->sites_key[*it] = i;
	i++;
    }

    return 0;
}

template<class StatT>
class OHaplotypes: public Ordered<StatT> {
public:
    OHaplotypes(loc_type type): Ordered<StatT>(type) { }

    int init_sites(vector<StatT> &);
    int init_haplotypes(vector<StatT> &);
    int order(vector<StatT> &, int, vector<CSLocus *> &, PopSum<CSLocus> *);
};

template<class StatT>
int 
OHaplotypes<StatT>::order(vector<StatT> &sites, int pop_id, 
			  vector<CSLocus *> &sorted_loci, PopSum<CSLocus> *psum) 
{
    return 0;
};

#endif // __ORDERED_H__
