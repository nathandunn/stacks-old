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
    struct SNPData {
        size_t tot_depth;
        Counts<Nt2> nt_depths;

        SNPData() : tot_depth(0) {}
    };

    int            id;            // Stack ID
    int            merge_partner; // Stack ID of merged datum, if this datum was merged/phased from two, overlapping datums.
    int            len;           // Length of locus
    bool           corrected;     // Has this genotype call been corrected
    char          *model;         // String representing SNP model output for each nucleotide at this locus.
    char          *gtype;         // Genotype
    char          *trans_gtype;   // Translated Genotype
    char          *cigar;         // CIGAR string describing how the datum aligns to the catalog locus.
    vector<char *> obshap;        // Observed Haplotypes
    vector<SNPData> snpdata;

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
     * [len] the length of the c-locus
     * [model] @
     * [cigar] is left to NULL.
     * [obshap] haplotypes, relative to [model].
     * [tot_depth], [nt_depths] @
     *
     * Other members need not be set, c.f. `populate(VcfRecord>&)`.
     */

    string model;
    model.resize(cloc.len, 'U');
    pair<string,string> obshaps;

    for (size_t sample = 0; sample < mpopi.samples().size(); ++sample) {
        size_t sample_vcf_i = header.sample_index(mpopi.samples()[sample].name);
    try {
        bool   no_data      = true;

        //
        // Get the sample's model string.
        //

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
            //
            // Create the Datum (see rules above).
            //

            Datum* d = new Datum();
            locdata[sample] = d;

            d->id = cloc.id;
            d->len = cloc.len;
            d->model = new char[cloc.len+1];
            strncpy(d->model, model.c_str(), cloc.len+1);
            d->cigar = NULL;
            if (cloc.snps.empty()) {
                //
                // Monomorphic locus.
                //

                d->obshap.push_back(new char[10]);
                strncpy(d->obshap[0], "consensus", 10);

            } else {
                //
                // Polymorphic locus.
                //

                // Get the list of SNP records.
                vector<const VcfRecord*> snp_records;
                for (auto& rec : cloc_records) {
                    // N.B. Columns and record indexes can be out of phase as
                    // columns without coverage don't have a record.
                    if (!rec.is_monomorphic()) {
                        assert(rec.pos() == cloc.snps.at(snp_records.size())->col);
                        snp_records.push_back(&rec);
                    }
                }
                assert(snp_records.size() == cloc.snps.size());

                // Save extra information for the SNP records.
                for (const VcfRecord* rec : snp_records) {
                try {
                    d->snpdata.push_back(Datum::SNPData());
                    const char* sample_gt = rec->find_sample(sample_vcf_i);
                    if (d->model[rec->pos()] == 'U')
                        continue;

                    // DP.
                    size_t dp_index = rec->index_of_gt_subfield("DP");
                    const char* dp_str = VcfRecord::util::find_gt_subfield(sample_gt, dp_index);
                    if (dp_str == NULL) {
                        cerr << "Error: DP field is missing.\n";
                        throw exception();
                    }
                    long tot_depth = strtol(dp_str, NULL, 10);
                    if (tot_depth < 0) {
                        cerr << "Error: Bad DP field.\n";
                        throw exception();
                    }
                    d->snpdata.back().tot_depth = tot_depth;

                    // AD.
                    size_t ad_index = rec->index_of_gt_subfield("AD");
                    const char* ad_str = VcfRecord::util::find_gt_subfield(sample_gt, ad_index);
                    if (ad_str == NULL) {
                        cerr << "Error: AD field is missing.\n";
                        throw exception();
                    }
                    size_t i = 0;
                    --ad_str;
                    do {
                        ++ad_str;
                        long ad = strtol(ad_str, const_cast<char**>(&ad_str), 10);
                        if (ad < 0)
                            throw exception();
                        d->snpdata.back().nt_depths.increment(Nt2(*rec->find_allele(i)), ad);
                        ++i;
                    } while (*ad_str == ',');
                    if (i != rec->count_alleles()) {
                        cerr << "Error: Bad AD field.\n";
                        throw exception();
                    }

                } catch (exception&) {
                    cerr << "Error: In VCF record '" << rec->chrom() << ":" << rec->pos()+1 << "'.\n";
                    throw;
                }
                }

                // Build the haplotypes from the records.
                VcfRecord::util::build_haps(obshaps, snp_records, sample_vcf_i);
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
    } catch (exception&) {
        cerr << "Error: At the " << (sample_vcf_i+1) << "th sample.\n"
             << "Error: (Locus " << cloc.id << ")\n"
             << "Error: Bad GStacksVCF file.\n";
        throw;
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
     * When no likelihood information is available, the [gt+liks] are not set.
     * xxx (n.b. the parsing of depth information in VCF is not implemented as
     * of Oct, 2017.)
     *
     * The following members are left unset, on the premise that
     * "populations" does not use them :
     * corrected, genotype, trans_genotype
     *
     * [merge_partner] is set later by [merge_datums()] (in populations.cc).
     */
    array<const char*,5> alleles; // Max. "ACGT*"
    size_t n_alleles = 0;
    try {
        for (auto a = rec.begin_alleles(); a != rec.end_alleles(); ++a) {
            alleles.at(n_alleles) = *a;
            ++n_alleles;
        }
    } catch (std::out_of_range&) {
        // Max. "ACGT*"
        cerr << "Error: Expected at most 5 SNP alleles; in VCF record '"
             << rec.chrom() << ":" << rec.pos() << "'.\n";
        throw exception();
    }

    size_t dp_index = rec.index_of_gt_subfield("DP");
    size_t ad_index = rec.index_of_gt_subfield("AD");

    for (size_t s = 0; s < mpopi.samples().size(); ++s) {
        size_t vcf_index = header.sample_index(mpopi.samples()[s].name);
        const char* sample = rec.find_sample(vcf_index);

        pair<int, int> gt = rec.parse_genotype(sample);
        if (gt.first < 0
            || gt.second < 0
            || *alleles.at(gt.first) == '*'
            || *alleles.at(gt.second) == '*')
            // Missing or incomplete genotype.
            continue;

        long tot_depth = 0;
        if (dp_index != SIZE_MAX) {
            const char* dp_str = VcfRecord::util::find_gt_subfield(sample, dp_index);
            tot_depth = strtol(dp_str, NULL, 10);
            if (tot_depth < 0) {
                cerr << "Warning: Bad DP field; in VCF record '"
                     << rec.chrom() << ":" << rec.pos() << "'.\n";
                tot_depth = 0;
            }
        }

        array<long,5> ad;
        ad.fill(0);
        if (ad_index != SIZE_MAX) {
            const char* ad_str = VcfRecord::util::find_gt_subfield(sample, ad_index);
            try {
                size_t i = 0;
                const char* start = ad_str;
                --start;
                do {
                    ++start;
                    ad.at(i) = strtol(start, const_cast<char**>(&start), 10);
                    if (ad[i] < 0)
                        throw exception();
                    ++i;
                } while (*start == ',');
                if (i != n_alleles)
                    throw exception();
            } catch (exception& e) {
                cerr << "Warning: Badly formatted AD field '" << ad_str
                     << "' in VCF record '" << rec.chrom() << ":" << rec.pos() << "'.\n";
                ad.fill(0);
            }
        }

        Datum* d = new Datum();
        locdata[s] = d;
        ++cloc.cnt;
        ++cloc.hcnt;

        // id, len, lnl
        d->id  = cloc.id;
        d->len = cloc.len;
        d->model = new char[2];
        d->model[1] = '\0';

        // depth
        d->snpdata.push_back(Datum::SNPData());
        d->snpdata[0].tot_depth = tot_depth;
        for (size_t i = 0; i < n_alleles; ++i)
            d->snpdata[0].nt_depths.increment(Nt2(*alleles.at(i)), ad[i]);

        // model, obshap
        if (gt.first == gt.second) {
            d->model[0] = 'O';
            char* a = new char[2];
            a[1] = '\0';
            a[0] = *alleles.at(gt.first);
            d->obshap.push_back(a);
        } else {
            d->model[0] = 'E';
            char* a0 = new char[2];
            char* a1 = new char[2];
            a0[0] = *alleles.at(gt.first);
            a1[0] = *alleles.at(gt.second);
            a0[1] = '\0';
            a1[1] = '\0';
            d->obshap.push_back(a0);
            d->obshap.push_back(a1);
        }
    }
}

#endif // __POPMAP_H__
