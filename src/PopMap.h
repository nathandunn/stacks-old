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
#include <cstring>
#include <string>
#include <vector>
#include <map>
#include <set>
#include <numeric>
using std::accumulate;
#include <algorithm>
#include <utility>

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
    const MetaPopInfo& metapopinfo;
    int      num_loci;
    set<pair<int, int> > blacklist;
    Datum ***data;
    map<int, int> locus_order;  // LocusID => ArrayIndex; map catalog IDs to their first dimension
                                // position in the Datum array.
    map<int, int> rev_locus_order;
    map<string, vector<LocusT *> > ordered_loci_; // Loci ordered by genomic position

public:
    //
    // Initialize & amend the PopMap.
    //

    PopMap(const MetaPopInfo& mpopi, int n_loci);
    ~PopMap();

    // Populates the PopMap based on sstack matches files.
    // The catalog is modified (LocusT must be CSLocus, and
    // members [cnt, hcnt, confounded_cnt] are modified).
    int populate(map<int, LocusT*>& catalog, const vector<vector<CatMatch *> >& matches);

    // Populates the PopMap based on Stacks (v2) files.
    void populate(map<int, LocusT*>& catalog, const unordered_map<int,vector<VcfRecord>>& cloci_records, const VcfHeader& header);
    static void populate_locus(Datum **locdata, const LocusT& cloc, const vector<VcfRecord>& cloc_records, const VcfHeader& header, const MetaPopInfo& mpopi);
    // Populates one locus based on external VCF files.
    static void populate_locus(Datum **locdata, LocusT& cloc, const VcfRecord record, const VcfHeader& header, const MetaPopInfo &mpopi);

    // Populates the PopMap based on VCF (SNP) records.
    // The catalog is modified (LocusT must be CSLocus, and
    // members [cnt, hcnt] are modified).
    // N.B. The IDs of the loci in the catalog MUST be the same
    // as the indexes in the records vector.
    int populate(map<int, LocusT*>& catalog, const vector<VcfRecord>& records, const VcfHeader& header);

    int order_loci(const map<int, LocusT*>& catalog);
    int prune(set<int>& loc_ids);
    map<string, vector<LocusT *> >& ordered_loci_nconst() {return ordered_loci_;}

    //
    // Obtain indexes, etc.
    //

    int loci_cnt() const { return this->num_loci; }
    int locus_index(int id) const {return locus_order.at(id);}
    int rev_locus_index(int index) const {try {return rev_locus_order.at(index);} catch(exception&) {return -1;}}

    int sample_cnt() const { return metapopinfo.samples().size(); }
    int sample_index(int id) const {try {return metapopinfo.get_sample_index(id);} catch (exception&) {return -1;}}
    int rev_sample_index(int index) const {return metapopinfo.samples().at(index).id;}

    const map<string,vector<LocusT*>>& ordered_loci() const {return ordered_loci_;}

    bool blacklisted(int loc_id, int sample_id) const {return blacklist.count({sample_id, loc_id});}

    //
    // Access the Datums.
    //

    Datum **locus(int id);
    Datum  *datum(int loc_id, int sample_id);

};

template<class LocusT>
PopMap<LocusT>::PopMap(const MetaPopInfo& mpopi, int num_loci): metapopinfo(mpopi)
{
    this->data = new Datum **[num_loci];

    for (int i = 0; i < num_loci; i++) {
        this->data[i] = new Datum *[metapopinfo.samples().size()];

        for (size_t j = 0; j < metapopinfo.samples().size(); j++)
            this->data[i][j] = NULL;
    }

    this->num_loci    = num_loci;
}

template<class LocusT>
PopMap<LocusT>::~PopMap() {
    for (int i = 0; i < this->num_loci; i++) {
        for (size_t j = 0; j < metapopinfo.samples().size(); j++)
            delete this->data[i][j];
        delete [] this->data[i];
    }
    delete [] this->data;
}

template<class LocusT>
void PopMap<LocusT>::populate(map<int, LocusT*>& catalog,
                              const unordered_map<int,vector<VcfRecord>>& cloci_records,
                              const VcfHeader& header)
{
    // Initalize [locus_order], [rev_locus_order].
    size_t i = 0;
    for (auto& cloc : catalog) {
        locus_order[cloc.second->id] = i;
        rev_locus_order[i] = cloc.second->id;
        ++i;
    }

    // Initialize [ordered_loci].
    order_loci(catalog);

    for (auto& cloc_pair : catalog) {
        LocusT& cloc = *cloc_pair.second;
        const vector<VcfRecord>& cloc_records = cloci_records.at(cloc.id);
        Datum** locdata = this->data[this->locus_index(cloc.id)];

        populate_locus(locdata, cloc, cloc_records, header, this->metapopinfo);
        for (int s = 0; s < sample_cnt(); ++s) {
            if (locdata[s] != NULL) {
                ++cloc.cnt;
                ++cloc.hcnt;
            }
        }
    }
}

template<class LocusT>
void PopMap<LocusT>::populate_locus(Datum** locdata,
                                    const LocusT& cloc,
                                    const vector<VcfRecord>& cloc_records,
                                    const VcfHeader& header,
                                    const MetaPopInfo& mpopi)
{
    /*
     * Fill the PopMap.
     *
     * We observe the following rules to create the Datums (@ when the value
     * is obvious) :
     * [id] the c-locus ID.
     * [len] the full length of the c-locus
     * [model] @
     * [cigar] is left to NULL.
     * [obshap] haplotypes, relative to [model].
     * [tot_depth], [depths] set to 0 as they don't make sense for
     *     variable-coverage data.
     * [lnl] is set to 0.
     *
     * Other members need not be set, c.f. `populate(vector<VcfRecord>&)`.
     */
    string model;
    model.resize(cloc.len, 'U');
    pair<string,string> obshaps;

    for (size_t sample = 0; sample < mpopi.samples().size(); ++sample) {

        size_t sample_vcf_i = header.sample_index(mpopi.samples()[sample].name);
        bool   no_data      = true;

        // Get the sample's model string.
        auto rec = cloc_records.begin();
        for (size_t col = 0; col < cloc.len; ++col) {
            if (rec == cloc_records.end()) {
                // No more VCF records.
                break;
            } else if (rec->pos() != col) {
                // No VCF record for this column, skip it.
                assert(model[col] == 'U'); // Position never changes as it's missing for all samples.
                continue;
            }

            const char* gt_str = rec->find_sample(sample_vcf_i);
            if (rec->is_monomorphic()) {
                // Monomorphic site.
                assert(rec->count_formats() == 1 && strcmp(rec->format0(),"DP")==0); // Only the samples overall depths are given.
                if (gt_str[0] == '.') {
                    model[col] = 'U';
                } else {
                    if (no_data)
                        no_data = false;
                    model[col] = 'O';
                }
            } else {
                // Polymorphic site.
                pair<int,int> gt = rec->parse_genotype_nochecks(gt_str);
                if (gt.first == -1) {
                    model[col] = 'U';
                } else {
                    if (no_data)
                        no_data = false;
                    model[col] = gt.first == gt.second ? 'O' : 'E';
                }
            }
            ++rec;
        }
        assert(rec == cloc_records.end());

        if (!no_data) {
            // Create the Datum (see rules above).
            Datum* d = new Datum();
            locdata[sample] = d;

            d->id = cloc.id;
            d->len = cloc.len;
            d->model = new char[cloc.len+1];
            strncpy(d->model, model.c_str(), cloc.len+1);
            d->cigar = NULL;
            if (cloc.snps.empty()) {
                d->obshap.push_back(new char[10]);
                strncpy(d->obshap[0], "consensus", 10);
            } else {
                // Build the haplotypes from the records. (Note: We can't get
                // the indexes of the SNP records directly from `cloc.snps`
                // as the records vector skips columns without any calls.
                vector<const VcfRecord*> snp_records;
                for (auto& rec : cloc_records) {
                    if (!rec.is_monomorphic()) {
                        assert(rec.pos() == cloc.snps.at(snp_records.size())->col);
                        snp_records.push_back(&rec);
                    }
                }
                assert(snp_records.size() == cloc.snps.size());
                VcfRecord::util::build_haps(obshaps, snp_records, sample_vcf_i);

                // Record the haplotypes.
                for (size_t i=0; i<2; ++i)
                    d->obshap.push_back(new char[snp_records.size()+1]);
                strncpy(d->obshap[0], obshaps.first.c_str(), snp_records.size()+1);
                strncpy(d->obshap[1], obshaps.second.c_str(), snp_records.size()+1);

                // Discard het calls that aren't phased, so that the model is
                // always 'U' (rather than 'E') where the obshaps are 'N'.
                for (size_t snp_i=0; snp_i<cloc.snps.size(); ++snp_i) {
                    if (d->obshap[0][snp_i] == 'N') {
                        size_t col = cloc.snps[snp_i]->col;
                        assert(d->model[col] == 'U' || d->model[col] == 'E');
                        if (d->model[col] == 'E')
                            d->model[col] = 'U';
                    }
                }
            }
            for (size_t i=0; i<d->obshap.size(); ++i)
                d->depth.push_back(0);
            d->tot_depth = 0;
            d->lnl = 0;
        }
    }
}

template<class LocusT>
int PopMap<LocusT>::populate(map<int, LocusT*> &catalog,
                             const vector<vector<CatMatch *> > &matches)
{
    //
    // Create an index showing what position each catalog locus is stored at in the datum
    // array. Create a second index allowing ordering of Loci by genomic position.
    //
    typename map<int, LocusT*>::const_iterator it;
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
            sample = metapopinfo.get_sample_index(matches[i][j]->sample_id);

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
int PopMap<LocusT>::populate(map<int, LocusT*>& catalog,
                             const vector<VcfRecord>& records,
                             const VcfHeader& header)
{
    // Initalize [locus_order], [rev_locus_order].
    size_t loc_index = 0;
    for (typename map<int, LocusT*>::iterator l = catalog.begin(); l != catalog.end(); ++l) {
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
     *     If one of a sample's VCF alleles is missing ('.') or has an index
     *     corresponding to the special '*' allele, the Datum* is left to NULL.
     *
     * [cigar] is left to NULL. (Updated May 24, 2017. It used not to exist.)
     *
     * When no depth information is available, [tot_depth] and the [depths]
     * of all alleles are set to 0.
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
     */

    loc_index = 0;
    for (typename map<int, LocusT*>::iterator l = catalog.begin(); l != catalog.end(); ++l) {

        LocusT* loc = l->second;

        const VcfRecord& rec = records[loc->id]; // n.b. assumes locus ID == record index.
        array<const char*,5> alleles; // Max. "ACGT*"
        {
            size_t i=0;
            for (auto a=rec.begin_alleles(); a!=rec.end_alleles(); ++a) {
                alleles.at(i) = *a;
                ++i;
            }
        }
        size_t ad_index;
        ad_index = rec.index_of_gt_subfield("AD");
        size_t n_alleles = rec.count_alleles();

        vector<int> ad;
        for (size_t s = 0; s < metapopinfo.samples().size(); ++s) {
            size_t vcf_index = header.sample_index(metapopinfo.samples()[s].name);
            const char* sample = rec.find_sample(vcf_index);

            pair<int, int> gt = rec.parse_genotype(sample);
            if (gt.first < 0
                    || gt.second < 0
                    || strcmp(alleles.at(gt.first),"*")==0
                    || strcmp(alleles.at(gt.second),"*")==0)
                // Missing or incomplete genotype.
                continue;

            if (ad_index != SIZE_MAX) {
                ad.clear();
                string ad_str = rec.parse_gt_subfield(sample, ad_index);
                size_t start = 0;
                size_t coma;
                try {
                     while ((coma = ad_str.find(',', start)) != string::npos) {
                         ad.push_back(stoi(ad_str.substr(start,coma)));
                         if (ad.back() < 0)
                             throw exception();
                         start=coma+1;
                     }
                     ad.push_back(stoi(ad_str.substr(start)));
                     if (ad.back() < 0)
                         throw exception();
                     if (ad.size() != n_alleles)
                             throw exception();
                 } catch (exception& e) {
                     cerr << "Warning: Badly formatted AD string '" << ad_str
                          << "' in VCF record '" << rec.chrom() << ":" << rec.pos() << "'.\n";
                     ad = vector<int>(n_alleles, 0);
                 }
            }

            Datum* d = new Datum();
            data[loc_index][s] = d;
            ++loc->cnt;
            ++loc->hcnt;

            // id, len, lnl
            d->id = loc->id;
            d->len = loc->len;
            d->lnl = 0;

            // model, obshap, depth
            d->model = new char[2];
            if (gt.first == gt.second) {
                strcpy(d->model, "O");
                const char* allele = alleles.at(gt.first);
                size_t len = strlen(allele);
                d->obshap.push_back(new char[len+1]);
                memcpy(d->obshap[0], allele, len+1);
                if (ad_index != SIZE_MAX)
                    d->depth = { ad[gt.first] };
                else
                    d->depth = {0};
            } else {
                strcpy(d->model, "E");
                const char* allele0 = alleles.at(gt.first);
                const char* allele1 = alleles.at(gt.second);
                size_t len0 = strlen(allele0);
                size_t len1 = strlen(allele1);
                d->obshap.push_back(new char[len0+1]);
                d->obshap.push_back(new char[len1+1]);
                memcpy(d->obshap[0], allele0, len0+1);
                memcpy(d->obshap[1], allele1, len1+1);
                if (ad_index != SIZE_MAX)
                    d->depth = {ad[gt.first], ad[gt.second]};
                else
                    d->depth = {0, 0};
            }

            // tot_depth
            d->tot_depth = std::accumulate(d->depth.begin(), d->depth.end(), 0);
        }
        ++loc_index;
    }

    return 0;
}

template<class LocusT> void
PopMap<LocusT>::populate_locus(Datum **locdata, LocusT &cloc, const VcfRecord rec, const VcfHeader& header, const MetaPopInfo &mpopi)
{
    /*
     * Fill the PopMap.
     *
     * We observe the following rules to create the Datums :
     * [id] is the locus id.
     * [len] is the locus length (expected to be one, for a SNP)
     * [model] "E" or "O" according to the SAMPLE/GT field, or
     *     "U" if the GT field is absent.
     * [obshap] the nucleotide(s) observed for this SNP for this individual.
     *     If one of a sample's VCF alleles is missing ('.') or has an index
     *     corresponding to the special '*' allele, the Datum* is left to NULL.
     *
     * [cigar] is left to NULL. (Updated May 24, 2017. It used not to exist.)
     *
     * When no depth information is available, [tot_depth] and the [depths]
     * of all alleles are set to 0.
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
     */
    array<const char*,5> alleles; // Max. "ACGT*"
    {
        size_t i=0;
        for (auto a=rec.begin_alleles(); a!=rec.end_alleles(); ++a) {
            alleles.at(i) = *a;
            ++i;
        }
    }
    size_t ad_index;
    ad_index = rec.index_of_gt_subfield("AD");
    size_t n_alleles = rec.count_alleles();

    vector<int> ad;
    for (size_t s = 0; s < mpopi.samples().size(); ++s) {
        size_t vcf_index = header.sample_index(mpopi.samples()[s].name);
        const char* sample = rec.find_sample(vcf_index);

        pair<int, int> gt = rec.parse_genotype(sample);
        if (gt.first < 0
            || gt.second < 0
            || strcmp(alleles.at(gt.first),"*")==0
            || strcmp(alleles.at(gt.second),"*")==0)
            // Missing or incomplete genotype.
            continue;

        if (ad_index != SIZE_MAX) {
            ad.clear();
            string ad_str = rec.parse_gt_subfield(sample, ad_index);
            size_t start = 0;
            size_t coma;
            try {
                while ((coma = ad_str.find(',', start)) != string::npos) {
                    ad.push_back(stoi(ad_str.substr(start,coma)));
                    if (ad.back() < 0)
                        throw exception();
                    start=coma+1;
                }
                ad.push_back(stoi(ad_str.substr(start)));
                if (ad.back() < 0)
                    throw exception();
                if (ad.size() != n_alleles)
                    throw exception();
            } catch (exception& e) {
                cerr << "Warning: Badly formatted AD string '" << ad_str
                     << "' in VCF record '" << rec.chrom() << ":" << rec.pos() << "'.\n";
                ad = vector<int>(n_alleles, 0);
            }
        }

        Datum* d = new Datum();
        locdata[s] = d;
        ++cloc.cnt;
        ++cloc.hcnt;

        // id, len, lnl
        d->id  = cloc.id;
        d->len = cloc.len;
        d->lnl = 0;

        // model, obshap, depth
        d->model = new char[2];
        if (gt.first == gt.second) {
            strcpy(d->model, "O");
            const char* allele = alleles.at(gt.first);
            size_t len = strlen(allele);
            d->obshap.push_back(new char[len+1]);
            memcpy(d->obshap[0], allele, len+1);
            if (ad_index != SIZE_MAX)
                d->depth = { ad[gt.first] };
            else
                d->depth = {0};
        } else {
            strcpy(d->model, "E");
            const char* allele0 = alleles.at(gt.first);
            const char* allele1 = alleles.at(gt.second);
            size_t len0 = strlen(allele0);
            size_t len1 = strlen(allele1);
            d->obshap.push_back(new char[len0+1]);
            d->obshap.push_back(new char[len1+1]);
            memcpy(d->obshap[0], allele0, len0+1);
            memcpy(d->obshap[1], allele1, len1+1);
            if (ad_index != SIZE_MAX)
                d->depth = {ad[gt.first], ad[gt.second]};
            else
                d->depth = {0, 0};
        }

        // tot_depth
        d->tot_depth = std::accumulate(d->depth.begin(), d->depth.end(), 0);
    }
}

template<class LocusT>
int PopMap<LocusT>::order_loci(const map<int, LocusT*> &catalog)
{
    this->ordered_loci_.clear();

    typename map<int, LocusT*>::const_iterator it;

    for (it = catalog.begin(); it != catalog.end(); it++) {
        if (!it->second->loc.empty())
            this->ordered_loci_[it->second->loc.chr()].push_back(it->second);
    }

    //
    // Sort the catalog loci on each chromosome according to base pair.
    //
    typename map<string, vector<LocusT*> >::iterator cit;
    for (cit = this->ordered_loci_.begin(); cit != this->ordered_loci_.end(); cit++)
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
            for (size_t k = 0; k < metapopinfo.samples().size(); k++)
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

    for (cit = this->ordered_loci_.begin(); cit != this->ordered_loci_.end(); cit++) {
        for (uint k = 0; k < cit->second.size(); k++) {
            if (remove_ids.count(cit->second[k]->id) == 0)
                new_ordered_loci[cit->first].push_back(cit->second[k]);
        }
    }

    this->ordered_loci_.clear();
    this->ordered_loci_ = new_ordered_loci;

    for (cit = this->ordered_loci_.begin(); cit != this->ordered_loci_.end(); cit++)
        sort(cit->second.begin(), cit->second.end(), bp_compare);

    return new_size;
}

template<class LocusT>
Datum **PopMap<LocusT>::locus(int id) {
    return this->data[this->locus_order[id]];
}

template<class LocusT>
Datum  *PopMap<LocusT>::datum(int loc_id, int sample_id) {
    return this->data[this->locus_order[loc_id]][metapopinfo.get_sample_index(sample_id)];
}

#endif // __POPMAP_H__
