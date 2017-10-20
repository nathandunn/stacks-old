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
    vector<int>    depth;         // Stack depth of each matching allele
    bool           corrected;     // Has this genotype call been corrected
    char          *model;         // String representing SNP model output for each nucleotide at this locus.
    char          *gtype;         // Genotype
    char          *trans_gtype;   // Translated Genotype
    char          *cigar;         // CIGAR string describing how the datum aligns to the catalog locus.
    vector<char *> obshap;        // Observed Haplotypes

    Datum()  {
        this->id            = -1;
        this->corrected     = false;
        this->gtype         = NULL;
        this->trans_gtype   = NULL;
        this->model         = NULL;
        this->cigar         = NULL;
        this->len           = 0;
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

    size_t tot_depth(size_t snp_index) const; // TODO
    size_t allele_depth(size_t snp_index, char allele) const; //TODO
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

    // Populates the PopMap based on Stacks (v2) files.
    static void populate_locus(Datum **locdata, const LocusT& cloc, const vector<VcfRecord>& cloc_records, const VcfHeader& header, const MetaPopInfo& mpopi);
    // Populates one locus based on external VCF files.
    static void populate_locus(Datum **locdata, LocusT& cloc, const VcfRecord record, const VcfHeader& header, const MetaPopInfo &mpopi);

    //
    // Obtain indexes, etc.
    //

    int loci_cnt() const { return this->num_loci; }
    int locus_index(int id) const {return locus_order.at(id);}
    int rev_locus_index(int index) const {try {return rev_locus_order.at(index);} catch(exception&) {return -1;}}

    int sample_cnt() const { return metapopinfo.samples().size(); }
    int sample_index(int id) const {try {return metapopinfo.get_sample_index(id);} catch (exception&) {return -1;}}
    int rev_sample_index(int index) const {return metapopinfo.samples().at(index).id;}

    bool blacklisted(int loc_id, int sample_id) const {return blacklist.count({sample_id, loc_id});}

    //
    // Access the Datums.
    //

    Datum **locus(int id) { return this->data[this->locus_order[id]]; }
    Datum  *datum(int loc_id, int sample_id) { return this->locus(loc_id)[metapopinfo.get_sample_index(sample_id)]; }

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
        }
    }
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
    }
}

#endif // __POPMAP_H__
