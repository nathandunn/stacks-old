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

#ifndef __POPMAP_H__
#define __POPMAP_H__

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

class Datum {
public:
    int            id;           // Stack ID
    int            tot_depth;    // Stack depth
    vector<int>    depth;        // Stack depth of each matching allele
    bool           corrected;    // Has this genotype call been corrected
    char          *model;        // String representing SNP model output for each nucleotide at this locus.
    char          *gtype;        // Genotype
    char          *trans_gtype;  // Translated Genotype
    vector<char *> obshap;       // Observed Haplotypes
    vector<SNP *>  snps;
    Datum()  { corrected = false; gtype = NULL; trans_gtype = NULL; model = NULL; tot_depth = 0; }
    ~Datum() {
    	for (uint i = 0; i < this->obshap.size(); i++)
    	    delete [] this->obshap[i];
    	for (uint i = 0; i < this->snps.size(); i++)
	    delete this->snps[i];
    	delete [] this->gtype;
	delete [] this->trans_gtype;
	delete [] this->model;
    }
};

template<class LocusT=Locus>
class PopMap {
    set<pair<int, int> > blacklist;
    int      num_loci;
    int      num_samples;
    Datum ***data;
    map<int, int> locus_order;  // LocusID => ArrayIndex; map catalog IDs to their first dimension 
                                // position in the Datum array.
    map<int, int> rev_locus_order;
    map<int, int> sample_order; // SampleID => ArrayIndex; map defining at what position in 
                                // the second dimension of the datum array each sample is stored.
    map<int, int> rev_sample_order;

public:
    map<string, vector<LocusT *> > ordered_loci; // Loci ordered by genomic position

    PopMap(int, int);
    ~PopMap();

    int populate(vector<int> &, map<int, LocusT*> &, vector<vector<CatMatch *> > &, int);

    int loci_cnt() { return this->num_loci; }
    int rev_locus_index(int index) { return this->rev_locus_order[index]; }
    int sample_cnt() { return this->num_samples; }
    int rev_sample_index(int index) { return this->rev_sample_order[index]; }

    Datum **locus(int);
    Datum  *datum(int, int);
    bool    blacklisted(int, int);
};

template<class LocusT>
PopMap<LocusT>::PopMap(int num_samples, int num_loci) {
    this->data = new Datum **[num_loci];

    for (int i = 0; i < num_loci; i++) {
	this->data[i] = new Datum *[num_samples];

	for (int j = 0; j < num_samples; j++)
	    this->data[i][j] = NULL;
    }

    this->num_samples = num_samples;
    this->num_loci    = num_loci;
}

template<class LocusT>
PopMap<LocusT>::~PopMap() {
    for (int i = 0; i < this->num_loci; i++) {
	for (int j = 0; j < this->num_samples; j++)
	    delete this->data[i][j];
	delete [] this->data[i];
    }
    delete [] this->data;
}

template<class LocusT>
int PopMap<LocusT>::populate(vector<int> &sample_ids,
			     map<int, LocusT*> &catalog,
			     vector<vector<CatMatch *> > &matches,
			     int min_stack_depth=0) {
    //
    // Record the array position of each sample that we will load.
    //
    for (uint i = 0; i < sample_ids.size(); i++) {
	this->sample_order[sample_ids[i]] = i;
	this->rev_sample_order[i]         = sample_ids[i];
    }

    //
    // Create an index showing what position each catalog locus is stored at in the datum 
    // array. Create a second index allowing ordering of Loci by genomic position.
    //
    typename std::map<int, LocusT*>::iterator it;
    uint i = 0;
    for (it = catalog.begin(); it != catalog.end(); it++) {
	this->locus_order[it->first] = i;
	this->rev_locus_order[i]     = it->first;

	if (strlen(it->second->loc.chr) > 0)
	    this->ordered_loci[it->second->loc.chr].push_back(it->second);
	i++;
    }
    //
    // Sort the catalog loci on each chromosome according to base pair.
    //
    typename map<string, vector<LocusT*> >::iterator cit;
    for (cit = this->ordered_loci.begin(); cit != this->ordered_loci.end(); cit++)
	sort(cit->second.begin(), cit->second.end(), bp_compare);

    //
    // Populate the datum array
    //
    Datum *d;
    int    locus, sample;

    for (i = 0; i < matches.size(); i++) {
	for (uint j = 0; j < matches[i].size(); j++) {
	    sample = this->sample_order[matches[i][j]->sample_id];

	    if (this->locus_order.count(matches[i][j]->cat_id) == 0)
		continue;

	    locus  = this->locus_order[matches[i][j]->cat_id];

	    if (matches[i][j]->depth < min_stack_depth)
		continue;

	    // cerr << "Translating sample id: " << matches[i][j]->sample_id << " to index " << sample << "\n";
	    // cerr << "Translating locus id: " << matches[i][j]->cat_id << " to index " << locus << "\n";

	    if (this->data[locus][sample] == NULL) {

		if (this->blacklist.count(make_pair(matches[i][j]->sample_id, matches[i][j]->cat_id)) == 0) {
		    // cerr << "Creating new datum for tag ID: " << matches[i][j]->tag_id << "\n";
		    d = new Datum;
		    d->id = matches[i][j]->tag_id;
		    char *h = new char[strlen(matches[i][j]->haplotype) + 1];
		    strcpy(h, matches[i][j]->haplotype);
		    d->obshap.push_back(h);
		    d->depth.push_back(matches[i][j]->depth);
		    d->tot_depth += matches[i][j]->depth;
		    this->data[locus][sample] = d;

		    catalog[matches[i][j]->cat_id]->hcnt++;
		}
	    } else {
		// cerr << "  Adding haplotype to existing datum: sample: " << matches[i][j]->sample_id << ". tag: " << matches[i][j]->tag_id << "\n";
		//
		// Check that the IDs of the two matches are the same. If not, then two tags 
		// match this locus and the locus is invalid, set back to NULL.
		//
		if (matches[i][j]->tag_id == this->data[locus][sample]->id) {
		    char *h = new char[strlen(matches[i][j]->haplotype) + 1];
		    strcpy(h, matches[i][j]->haplotype);
		    this->data[locus][sample]->obshap.push_back(h);
		    this->data[locus][sample]->depth.push_back(matches[i][j]->depth);
		    this->data[locus][sample]->tot_depth += matches[i][j]->depth;
		} else {
		    //cerr << "    Deleting sample, multiple tag matches\n";
		    delete this->data[locus][sample];
		    this->data[locus][sample] = NULL;
		    this->blacklist.insert(make_pair(matches[i][j]->sample_id, matches[i][j]->cat_id));
		    catalog[matches[i][j]->cat_id]->hcnt--;
		}
	    }
	}
    }

    return 0;
}

template<class LocusT>
Datum **PopMap<LocusT>::locus(int locus) {
    return this->data[this->locus_order[locus]];
}

template<class LocusT>
Datum  *PopMap<LocusT>::datum(int locus, int sample) {
    return this->data[this->locus_order[locus]][this->sample_order[sample]];
}

template<class LocusT>
bool PopMap<LocusT>::blacklisted(int locus, int sample) {
    if (this->blacklist.count(make_pair(sample, locus)) > 0)
	return true;
    else
	return false;
}

#endif // __POPMAP_H__
