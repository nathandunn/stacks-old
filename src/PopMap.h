// -*-mode:c++; c-style:k&r; c-basic-offset:4;-*-
//
// Copyright 2011-2016, Julian Catchen <jcatchen@illinois.edu>
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

#include <exception>
using std::exception;
#include <cstring>
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
#include "locus.h"
#include "aln_utils.h"
#include "MetaPopInfo.h"
#include "Vcf.h"

class Datum {
public:
    int            id;            // Stack ID
    int            merge_partner; // Stack ID of merged datum, if this datum was merged/phased from two, overlapping datums.
    int            len;           // Length of locus
    int            tot_depth;     // Stack depth
    vector<int>    depth;         // Stack depth of each matching allele
    bool           corrected;     // Has this genotype call been corrected
    char          *model;         // String representing SNP model output for each nucleotide at this locus.
    char          *gtype;         // Genotype
    char          *trans_gtype;   // Translated Genotype
    char          *cigar;         // CIGAR string describing how the datum aligns to the catalog locus.
    double         lnl;           // Log likelihood of this locus.
    vector<char *> obshap;        // Observed Haplotypes
    vector<SNP *>  snps;          // All calls for this sample for this locus. size() is [len].
    Datum()  {
        this->id            = -1;
        this->corrected     = false;
        this->gtype         = NULL;
        this->trans_gtype   = NULL;
        this->model         = NULL;
        this->cigar         = NULL;
        this->tot_depth     = 0;
        this->len           = 0;
        this->lnl           = 0.0;
        this->merge_partner = 0;
    }
    ~Datum() {
        for (uint i = 0; i < this->obshap.size(); i++)
            delete [] this->obshap[i];
        for (uint i = 0; i < this->snps.size(); i++)
            delete this->snps[i];
        delete [] this->gtype;
        delete [] this->trans_gtype;
        delete [] this->model;
        delete [] this->cigar;
    }
    int add_model(const char *model)
    {
        if (this->cigar == NULL) {
            this->len   = strlen(model);
            this->model = new char[this->len + 1];
            strcpy(this->model, model);

        } else {
            vector<pair<char, uint> > c;
            this->len   = parse_cigar(this->cigar, c);
            this->model = new char[this->len + 1];
            apply_cigar_to_model_seq(this->model, this->len, model, c);
            // cerr << "Cigar: " << this->cigar << "\n"
            //      << "Old model:    " << model << "\n"
            //      << "Gapped model: " << this->model << "\n";
        }
        return 0;
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

    // Sets [sample_order] and [rev_sample_order]. The array indexes
    // will be the same as those of [mpopi.samples_].
    void add_metapop_info(const MetaPopInfo& mpopi);

public:
    map<string, vector<LocusT *> > ordered_loci; // Loci ordered by genomic position

    PopMap(int n_samples, int n_loci);
    ~PopMap();

    // Populates the PopMap based on sstack matches files.
    // The catalog is modified (LocusT must be CSLocus, and
    // members [cnt, hcnt, confounded_cnt] are modified).
    int populate(const vector<int>& sample_ids, map<int, LocusT*>& catalog, const vector<vector<CatMatch *> >& matches);

    // Populates the PopMap based on VCF records.
    // The catalog is modified (LocusT must be CSLocus, and
    // members [cnt, hcnt] are modified).
    // N.B. The IDs of the loci in the catalog MUST be the same
    // as the indexes in the records vector.
    int populate(const MetaPopInfo& mpopi, map<int, LocusT*>& catalog, const vector<VcfRecord>& records, const VcfHeader& header);

    int order_loci(const map<int, LocusT*>& catalog);
    int prune(set<int>& loc_ids);

    int loci_cnt() { return this->num_loci; }
    int locus_index(int id) {return this->locus_order.at(id);}
    int rev_locus_index(int index) { if (this->rev_locus_order.count(index) == 0) return -1; return this->rev_locus_order[index]; }

    int sample_cnt() { return this->num_samples; }
    int sample_index(int id) { if (this->sample_order.count(id) == 0) return -1; return this->sample_order[id]; }
    int rev_sample_index(int index) { if (this->rev_sample_order.count(index) == 0) return -1; return this->rev_sample_order[index]; }

    Datum **locus(int id);
    Datum  *datum(int loc_id, int sample_id);

    bool    blacklisted(int loc_id, int sample_id);
};

template<class LocusT>
void PopMap<LocusT>::add_metapop_info(const MetaPopInfo& mpopi) {
    for(size_t i = 0; i < mpopi.samples().size(); ++i) {
        sample_order[mpopi.samples()[i].id] = i;
        rev_sample_order[i] = mpopi.samples()[i].id;
    }
}

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
int PopMap<LocusT>::populate(const vector<int> &sample_ids,
                             map<int, LocusT*> &catalog,
                             const vector<vector<CatMatch *> > &matches) {
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
    typename std::map<int, LocusT*>::const_iterator it;
    uint i = 0;
    for (it = catalog.begin(); it != catalog.end(); it++) {
        this->locus_order[it->first] = i;
        this->rev_locus_order[i]     = it->first;
        i++;
    }

    //
    // Sort the catalog loci on each chromosome according to base pair.
    //
    this->order_loci(catalog);

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
                    d->lnl        = matches[i][j]->lnl;

                    if (matches[i][j]->cigar != NULL) {
                        d->cigar = new char[strlen(matches[i][j]->cigar) + 1];
                        strcpy(d->cigar, matches[i][j]->cigar);
                    }

                    this->data[locus][sample] = d;

                    catalog[matches[i][j]->cat_id]->hcnt++;
                    catalog[matches[i][j]->cat_id]->cnt++;
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

                    if (matches[i][j]->cigar != NULL && strcmp(this->data[locus][sample]->cigar, matches[i][j]->cigar) != 0)
                        cerr << "Warning: disparate CIGAR strings, catalog locus " << matches[i][j]->cat_id
                             << "; sample ID: " << matches[i][j]->sample_id << "; sample locus ID: " << matches[i][j]->tag_id
                             << "; datum cigar: " << this->data[locus][sample]->cigar << "; matches cigar: " << matches[i][j]->cigar << "\n";
                    
                    this->data[locus][sample]->obshap.push_back(h);
                    this->data[locus][sample]->depth.push_back(matches[i][j]->depth);
                    this->data[locus][sample]->tot_depth += matches[i][j]->depth;
                    this->data[locus][sample]->lnl        = matches[i][j]->lnl;

                } else {
                    //cerr << "    Deleting sample, multiple tag matches\n";
                    delete this->data[locus][sample];
                    this->data[locus][sample] = NULL;
                    this->blacklist.insert(make_pair(matches[i][j]->sample_id, matches[i][j]->cat_id));
                    catalog[matches[i][j]->cat_id]->hcnt--;
                    catalog[matches[i][j]->cat_id]->confounded_cnt++;
                }
            }
        }
    }

    return 0;
}

template<class LocusT>
int PopMap<LocusT>::populate(const MetaPopInfo& mpopi,
             map<int, LocusT*>& catalog,
             const vector<VcfRecord>& records,
             const VcfHeader& header) {

    // Initialize [sample_order], [rev_sample_order].
    add_metapop_info(mpopi);

    // Initalize [locus_order], [rev_locus_order].
    size_t loc_index = 0;
    for (typename map<int, LocusT*>::iterator
            l = catalog.begin();
            l != catalog.end();
            ++l) {
        locus_order[l->first] = loc_index;
        rev_locus_order[loc_index] = l->first;
        ++loc_index;
    }

    // Initialize [ordered_loci].
    order_loci(catalog);

    /*
     * Fill the PopMap.
     *
     * We observe the following rules to create the Datums :
     * [id] is the locus id.
     * [len] is the locus length (expected to be one, for a SNP)
     * [model] "E" or "O" according to the SAMPLE/GT field, or
     *     "U" if the GT field is absent.
     * [obshap] the nucleotide(s) observed for this SNP for this individual.
     *
     * When no depth information is available, [tot_depth] and the [depths]
     * of all alleles are set to 0.
     * (n.b. the parsing of depth information in VCF is not implemented as
     * of Mar 21, 2016.)
     *
     * When no likelihood information is available, [lnl] is set to 0.
     * (n.b. the parsing of depth information in VCF is not implemented as
     * of Mar 21, 2016.)
     *
     * The following members are left unset, on the premise that
     * "populations" does not use them :
     * corrected, genotype, trans_genotype
     *
     * [merge_partner] is set later by [merge_datums()] (in populations.cc).
     * [snps] is only used by [write_vcf()] and [write_vcf_strict()], which
     *     first call [populate_snp_calls()] then use the SNP data in the
     *     datum.
     */

    loc_index = 0;
    for (typename map<int, LocusT*>::iterator
            l = catalog.begin();
            l != catalog.end();
            ++l) {
        LocusT* loc = l->second;
        const VcfRecord& rec = records[loc->id]; // n.b. assumes locus ID == record index.

        for (size_t s = 0; s < mpopi.samples().size(); ++s) {
            size_t vcf_index = header.sample_indexes().at(mpopi.samples()[s].name);
            pair<int, int> gt = rec.parse_genotype(rec.samples.at(vcf_index));
            if (gt == pair<int,int>(-1,-1))
                continue;

            Datum* d = new Datum();
            data[loc_index][s] = d;
            ++loc->cnt;
            ++loc->hcnt;

            // id, len, tot_depth, lnl
            d->id = loc->id;
            d->len = loc->len;
            d->tot_depth = 0;
            d->lnl = 0;

            // model, obshap, depth
            d->model = new char[2];
            if (gt.first == gt.second) {
                strcpy(d->model, "O");
                const string& allele = rec.allele(gt.first);
                d->obshap.push_back(new char[allele.size()+1]);
                strcpy(d->obshap[0], allele.c_str());
                d->depth.push_back(0);
            } else {
                strcpy(d->model, "E");
                const string& allele1 = rec.allele(gt.first);
                const string& allele2 = rec.allele(gt.second);
                d->obshap.push_back(new char[allele1.size()+1]);
                d->obshap.push_back(new char[allele2.size()+1]);
                strcpy(d->obshap[0], allele1.c_str());
                strcpy(d->obshap[1], allele2.c_str());
                d->depth.push_back(0);
                d->depth.push_back(0);
            }
        }
        ++loc_index;
    }

    return 0;
}

template<class LocusT>
int PopMap<LocusT>::order_loci(const map<int, LocusT*> &catalog)
{
    this->ordered_loci.clear();

    typename std::map<int, LocusT*>::const_iterator it;

    for (it = catalog.begin(); it != catalog.end(); it++) {
        if (strlen(it->second->loc.chr) > 0)
            this->ordered_loci[it->second->loc.chr].push_back(it->second);
    }

    //
    // Sort the catalog loci on each chromosome according to base pair.
    //
    typename map<string, vector<LocusT*> >::iterator cit;
    for (cit = this->ordered_loci.begin(); cit != this->ordered_loci.end(); cit++)
        sort(cit->second.begin(), cit->second.end(), bp_compare);

    return 0;
}

template<class LocusT>
int PopMap<LocusT>::prune(set<int> &remove_ids) {
    uint new_size = this->num_loci - remove_ids.size();
    uint loc_id;
    map<int, int> new_loc_order, new_rev_loc_order;

    Datum ***d = new Datum **[new_size];

    int j = 0;
    for (int i = 0; i < this->num_loci; i++) {

        loc_id = this->rev_locus_order[i];

        if (remove_ids.count(loc_id) == 0) {
            // Keep this locus.
            d[j] = this->data[i];
            new_loc_order[loc_id] = j;
            new_rev_loc_order[j] = loc_id;
            j++;

        } else {
            // Remove this locus.
            for (int k = 0; k < this->num_samples; k++)
                delete this->data[i][k];
            delete [] this->data[i];
        }
    }

    delete [] this->data;

    this->data = d;
    this->locus_order.clear();
    this->locus_order = new_loc_order;
    this->rev_locus_order.clear();
    this->rev_locus_order = new_rev_loc_order;
    this->num_loci = new_size;

    //
    // Re-sort the catalog loci on each chromosome according to base pair.
    //
    map<string, vector<LocusT *> > new_ordered_loci;
    typename map<string, vector<LocusT*> >::iterator cit;

    for (cit = this->ordered_loci.begin(); cit != this->ordered_loci.end(); cit++) {
        for (uint k = 0; k < cit->second.size(); k++) {
            if (remove_ids.count(cit->second[k]->id) == 0)
                new_ordered_loci[cit->first].push_back(cit->second[k]);
        }
    }

    this->ordered_loci.clear();
    this->ordered_loci = new_ordered_loci;

    for (cit = this->ordered_loci.begin(); cit != this->ordered_loci.end(); cit++)
        sort(cit->second.begin(), cit->second.end(), bp_compare);

    return new_size;
}

template<class LocusT>
Datum **PopMap<LocusT>::locus(int id) {
    return this->data[this->locus_order[id]];
}

template<class LocusT>
Datum  *PopMap<LocusT>::datum(int loc_id, int sample_id) {
    return this->data[this->locus_order[loc_id]][this->sample_order[sample_id]];
}

template<class LocusT>
bool PopMap<LocusT>::blacklisted(int locus, int sample) {
    if (this->blacklist.count(make_pair(sample, locus)) > 0)
        return true;
    else
        return false;
}

#endif // __POPMAP_H__
