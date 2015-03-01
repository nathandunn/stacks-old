// -*-mode:c++; c-style:k&r; c-basic-offset:4;-*-
//
// Copyright 2013, Julian Catchen <jcatchen@uoregon.edu>
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
// catalog_utils.cc -- common routines for manipulating catalog objects.
//
// Julian Catchen
// jcatchen@uoregon.edu
// University of Oregon
//
#include "catalog_utils.h"

int 
reduce_catalog(map<int, CSLocus *> &catalog, set<int> &whitelist, set<int> &blacklist) 
{
    map<int, CSLocus *> list;
    map<int, CSLocus *>::iterator it;
    CSLocus *loc;

    if (whitelist.size() == 0 && blacklist.size() == 0) 
	return 0;
 
    int i = 0;
    for (it = catalog.begin(); it != catalog.end(); it++) {
	loc = it->second;

	if (whitelist.size() > 0 && whitelist.count(loc->id) == 0) continue;
	if (blacklist.count(loc->id)) continue;

	list[it->first] = it->second;
	i++;
    }

    catalog = list;

    return i;
}

int 
implement_single_snp_whitelist(map<int, CSLocus *> &catalog, map<int, set<int> > &whitelist) 
{
    map<int, set<int> > new_wl;
    CSLocus *loc;

    if (whitelist.size() > 0) {
	map<int, set<int> >::iterator it;

	for (it = whitelist.begin(); it != whitelist.end(); it++) {
	    loc = catalog[it->first];

	    if (it->second.size() == 0)
		new_wl[loc->id].insert(loc->snps[0]->col);
	    else
		for (uint i = 0; i < loc->snps.size(); i++) {
		    if (it->second.count(loc->snps[i]->col) == 0)
			continue;	
		    new_wl[loc->id].insert(loc->snps[i]->col);
		    break;
		}
	}
    } else {
	map<int, CSLocus *>::iterator it;

	for (it = catalog.begin(); it != catalog.end(); it++) {
	    loc = it->second;

	    if (loc->snps.size() > 0) 
		new_wl[loc->id].insert(loc->snps[0]->col);
	}
    }

    whitelist = new_wl;

    return 0;
}

int 
implement_random_snp_whitelist(map<int, CSLocus *> &catalog, map<int, set<int> > &whitelist) 
{
    map<int, set<int> > new_wl;
    CSLocus *loc;
    uint index;

    if (whitelist.size() > 0) {
	map<int, set<int> >::iterator it;

	for (it = whitelist.begin(); it != whitelist.end(); it++) {
	    loc = catalog[it->first];

	    if (it->second.size() == 0) {
		index = rand() % loc->snps.size();
		new_wl[loc->id].insert(loc->snps[index]->col);
	    } else {
		do {
		    index = rand() % loc->snps.size();
		} while (it->second.count(loc->snps[index]->col) == 0);
		new_wl[loc->id].insert(loc->snps[index]->col);
	    }
	}
    } else {
	map<int, CSLocus *>::iterator it;

	for (it = catalog.begin(); it != catalog.end(); it++) {
	    loc = it->second;

	    if (loc->snps.size() > 0) {
		index = rand() % loc->snps.size();
		new_wl[loc->id].insert(loc->snps[index]->col);
	    }
	}
    }

    whitelist = new_wl;

    return 0;
}

int
check_whitelist_integrity(map<int, CSLocus *> &catalog, map<int, set<int> > &whitelist)
{
    if (whitelist.size() == 0) return 0;

    int rm_snps = 0;
    int rm_loci = 0;

    CSLocus *loc;
    map<int, set<int> >::iterator it;
    set<int>::iterator sit;
    map<int, set<int> > new_wl;

    cerr << "Checking the integrity of the whitelist...";

    for (it = whitelist.begin(); it != whitelist.end(); it++) {
	if (catalog.count(it->first) == 0) {
	    rm_loci++;
	    cerr << "\n  Removing locus " << it->first << " from whitelist as it does not exist in the catalog.";
	} else {
	    loc = catalog[it->first];

	    if (it->second.size() == 0) {
		new_wl.insert(make_pair(it->first, std::set<int>()));
		continue;
	    }

	    set<int> cat_snps;
	    for (uint i = 0; i < loc->snps.size(); i++)
		cat_snps.insert(loc->snps[i]->col);

	    for (sit = it->second.begin(); sit != it->second.end(); sit++)
		if (cat_snps.count(*sit)) {
		    new_wl[it->first].insert(*sit);
		} else {
		    rm_snps++;
		    cerr << "\n  Removing SNP at column " << *sit << " in locus " << it->first << " from whitelist as it does not exist in the catalog.";
		}
	}
    }

    whitelist = new_wl;

    if (rm_loci > 0 || rm_snps > 0) cerr << "\n";

    cerr << "done.\n"
	 << "Removed " << rm_loci << " loci and " << rm_snps << " SNPs from the whitelist that were not found in the catalog.\n";

    return 0;
}

int 
reduce_catalog(map<int, CSLocus *> &catalog, map<int, set<int> > &whitelist, set<int> &blacklist) 
{
    map<int, CSLocus *> list;
    map<int, CSLocus *>::iterator it;
    CSLocus *loc;

    if (whitelist.size() == 0 && blacklist.size() == 0) 
	return 0;
 
    int i = 0;
    for (it = catalog.begin(); it != catalog.end(); it++) {
	loc = it->second;

	if (whitelist.size() > 0 && whitelist.count(loc->id) == 0) continue;
	if (blacklist.count(loc->id)) continue;

	list[it->first] = it->second;
	i++;
    }

    catalog = list;

    return i;
}

int 
reduce_catalog_snps(map<int, CSLocus *> &catalog, map<int, set<int> > &whitelist, PopMap<CSLocus> *pmap) 
{
    map<int, CSLocus *>::iterator it;
    CSLocus *loc;
    Datum  **d;

    if (whitelist.size() == 0) 
	return 0;
 
    //
    // We want to prune out SNP objects that are not in the whitelist.
    //
    int           pos;
    vector<SNP *> tmp;
    vector<uint>  cols;
    set<string>   obshaps;
    set<string>::iterator sit;
    for (it = catalog.begin(); it != catalog.end(); it++) {
	loc = it->second;

	if (whitelist[loc->id].size() == 0)
	    continue;

	tmp.clear();
	cols.clear();

	d = pmap->locus(loc->id);

	for (uint i = 0; i < loc->snps.size(); i++) {
	    if (whitelist[loc->id].count(loc->snps[i]->col) > 0) {
		tmp.push_back(loc->snps[i]);
		cols.push_back(i);
	    } else {
		//
		// Change the model calls in the samples to no longer contain this SNP.
		//
		pos = loc->snps[i]->col;
		for (int j = 0; j < pmap->sample_cnt(); j++) {
		    if (d[j] == NULL || pos >= d[j]->len) 
			continue;
		    if (d[j]->model != NULL) {
			d[j]->model[pos] = 'U';
		    }
		}

		delete loc->snps[i];
	    }
	}
	loc->snps.clear();
	for (uint i = 0; i < tmp.size(); i++)
	    loc->snps.push_back(tmp[i]);

	//
	// Now we need to adjust the matched haplotypes to sync to 
	// the SNPs left in the catalog.
	//
	// Reducing the lengths of the haplotypes  may create 
	// redundant (shorter) haplotypes, we need to remove these.
	//
	for (int i = 0; i < pmap->sample_cnt(); i++) {
	    if (d[i] == NULL) continue;

	    for (uint j = 0; j < d[i]->obshap.size(); j++) {
		for (uint k = 0; k < cols.size(); k++)
		    d[i]->obshap[j][k] = d[i]->obshap[j][cols[k]];
		d[i]->obshap[j][cols.size()] = '\0';
		obshaps.insert(d[i]->obshap[j]);
	    }
	    uint j = 0;
	    for (sit = obshaps.begin(); sit != obshaps.end(); sit++) {
		strcpy(d[i]->obshap[j], (*sit).c_str());
		j++;
	    }
	    while (j < d[i]->obshap.size()) {
		delete [] d[i]->obshap[j];
		j++;
	    }
	    d[i]->obshap.resize(obshaps.size());
	    obshaps.clear();
	}
    }

    return 0;
}
