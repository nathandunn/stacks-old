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
        uint8_t gq;
        GtLiks gtliks;

        SNPData() : tot_depth(0), gq(UINT8_MAX) {}
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
    PopMap(const MetaPopInfo& mpopi, int n_loci);
    ~PopMap();

    // Populates the Locus & PopMap based on Stacks (v2) files.
    static void populate_internal(
        LocusT* cloc,
        Datum** locdata,
        const Seq& fasta_record,
        const vector<VcfRecord>& vcf_records,
        const VcfHeader& vcf_header,
        const MetaPopInfo* mpopi
    );
    // Populates the Locus & PopMap based on external VCF files.
    static void populate_external(
        LocusT* cloc,
        Datum** locdata,
        int cloc_id,
        const VcfRecord& record,
        const VcfHeader& header,
        const MetaPopInfo* mpopi
    );

    // TODO {
    static void populate_locus(
        Datum **locdata,
        LocusT& cloc,
        const VcfRecord record,
        const VcfHeader& header,
        const MetaPopInfo &mpopi
    );
    //TODO }

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
void PopMap<LocusT>::populate_internal(
    LocusT* cloc,
    Datum** locdata,
    const Seq& fasta_record,
    const vector<VcfRecord>& vcf_records,
    const VcfHeader& vcf_header,
    const MetaPopInfo* mpopi
) { try {
    assert(fasta_record.id != NULL);
    assert(fasta_record.seq != NULL);
    assert(!fasta_record.comment.empty());
    assert(!vcf_records.empty());

    // Parse the FASTA record.
    // ==========
    // Locus ID.
    cloc->id = is_integer(fasta_record.id);
    if (cloc->id < 0) {
        cerr << "Error: Unable to parse locus ID.\n";
        throw exception();
    }
    // Consensus sequence.
    cloc->add_consensus(fasta_record.seq);
    if (cloc->len < vcf_records.back().pos() + 1) {
        cerr << "Error: Bad consensus length.\n";
        throw exception();
    }
    // Locus position:
    // If the analysis is reference-based, there will be a ' pos=...' field on
    // the fasta ID line, placed as a comment adjacent to the locus ID.
    // Overlap: if this data is de novo, there will potentially be an overlap
    // length between the single and paired-end contigs.
    assert(cloc->loc.empty());
    const char *p, *q;
    p = fasta_record.comment.c_str();
    do {
        q = p;
        while (*q != '\0' && *q != ' ' && *q != '\t')
            ++q;
        if (strncmp(p, "pos=", 4) == 0) {
            p += 4;
            cloc->loc = PhyLoc(string(p, q));
        } else if (strncmp(p, "contig=overlapped:", 18) == 0) {
            p += 18;
            cloc->overlap = is_integer(p);
            cloc->pe_ctg = true;
        } else if (strncmp(p, "contig=separate", 15) == 0) {
            cloc->pe_ctg = true;
        }
        p = q;
        while(*p != '\0' && (*p == ' ' || *p == '\t'))
            ++p;
    } while (*p != '\0');
    if (cloc->loc.empty())
        cloc->loc = PhyLoc("", 0, strand_plus);

    // Parse the VCF records.
    // ==========

    // CSLocus.
    // ----------
    for (auto& rec : vcf_records) {
        if (!rec.is_monomorphic()) {
            SNP* snp = new SNP();
            cloc->snps.push_back(snp);
            snp->type = snp_type_het;
            snp->col = rec.pos();
            auto a = rec.begin_alleles();
            snp->rank_1 = **a;
            ++a;
            assert(a != rec.end_alleles());
            snp->rank_2 = **a;
            ++a;
            if (a != rec.end_alleles()) {
                snp->rank_3 = **a;
                ++a;
                if (a != rec.end_alleles())
                    snp->rank_4 = **a;
            }
            snp->lratio = 0.0;
        }
    }

    // Genotypes.
    // ----------
    size_t snp_i = 0;
    --snp_i;
    vector<Nt2> rec_alleles;
    size_t dp_index = SIZE_MAX;
    size_t ad_index = SIZE_MAX;
    size_t gq_index = SIZE_MAX;
    size_t gl_index = SIZE_MAX;
    for (auto& rec : vcf_records)
    { try {
        size_t col = rec.pos();
        bool snp_rec = !rec.is_monomorphic();
        if (snp_rec) {
            if (!rec.is_snp()) {
                cerr << "Error: Not a SNP.\n";
                throw exception();
            }
            ++snp_i;
            rec_alleles.clear();
            for(auto a = rec.begin_alleles(); a != rec.end_alleles(); ++a)
                rec_alleles.push_back(Nt2(**a));
            if (rec.index_of_gt_subfield("PS") != 1)
                throw exception();
            dp_index = rec.index_of_gt_subfield("DP");
            ad_index = rec.index_of_gt_subfield("AD");
            gq_index = rec.index_of_gt_subfield("GQ");
            gl_index = rec.index_of_gt_subfield("GL");
        } else {
            assert(rec.count_formats() == 1 && strcmp(rec.format0(),"DP")==0);
        }
        VcfRecord::iterator gt_itr = rec.begin_samples();
        for (size_t sample_vcf_i = 0;
            sample_vcf_i < vcf_header.samples().size();
            ++gt_itr, ++sample_vcf_i
        ) { try {
            assert(gt_itr != rec.end_samples());
            // Check if the sample has data.
            const char* gt_str = *gt_itr;
            if (gt_str[0] == '.')
                continue;
            size_t sample_mpopi_i = mpopi->get_sample_index(vcf_header.samples()[sample_vcf_i]); //TODO inefficient
            Datum* d = locdata[sample_mpopi_i];
            if (snp_rec) {
                // Check that this isn't an unphased HET.
                // N.B.: GStacks writes a single phase set per sample.
                pair<int,int> gt = rec.parse_genotype_nochecks(gt_str);
                assert(gt.first >= 0);
                if (gt.first != gt.second) {
                    const char* ps = strchr(gt_str, ':');
                    if (ps == NULL)
                        throw exception();
                    ++ps;
                    if (ps[0] == '.')
                        continue;
                }
            }
            // This sample has data for this locus for this column.
            if (d == NULL) {
                ++cloc->cnt;
                d = new Datum();
                locdata[sample_mpopi_i] = d;
                d->id = cloc->id;
                d->len = cloc->len;
                d->model = new char[cloc->len+1];
                memset(d->model, 'U', cloc->len);
                d->model[cloc->len] = '\0';
                if (!cloc->snps.empty()) {
                    d->snpdata = vector<Datum::SNPData>(cloc->snps.size());
                    for (size_t i=0; i<2; ++i) {
                        char* hap = new char[cloc->snps.size()+1];
                        d->obshap.push_back(hap);
                        memset(hap, 'N', cloc->snps.size());
                        hap[cloc->snps.size()] = '\0';
                    }
                } else {
                    d->obshap.push_back(new char[10]);
                    strncpy(d->obshap[0], "consensus", 10);
                }
            }
            // Handle the non-SNP-record case.
            if (!snp_rec) {
                d->model[col] = 'O';
                continue;
            }
            // SNP column.
            pair<int,int> gt = rec.parse_genotype_nochecks(gt_str);
            d->model[col] = (gt.first == gt.second ? 'O' : 'E');
            d->obshap[0][snp_i] = char(rec_alleles.at(gt.first));
            d->obshap[1][snp_i] = char(rec_alleles.at(gt.second));
            // Record additional information on the call.
            Datum::SNPData& s = d->snpdata[snp_i];
            s.tot_depth = VcfRecord::util::parse_gt_dp(gt_str, dp_index);
            s.nt_depths = VcfRecord::util::parse_gt_ad(gt_str, ad_index, rec_alleles);
            s.gq = VcfRecord::util::parse_gt_gq(gt_str, gq_index);
            s.gtliks = VcfRecord::util::parse_gt_gl(gt_str, gl_index, rec_alleles);
        } catch (exception&) {
            cerr << "Error: At the " << (sample_vcf_i+1) << "th sample.\n";
            throw;
        }}
    } catch (exception&) {
        cerr << "Error: In VCF record '" << rec.chrom() << ":" << rec.pos()+1 << "'.\n";
        throw;
    }}

    // Finalize the CSLocus (??).
    // ==========
    if (!cloc->snps.empty()) {
        string hap;
        for (size_t i = 0; i < mpopi->samples().size(); ++i) {
            Datum* d = locdata[i];
            if (d == NULL)
                continue;
            assert(d->obshap.size() == 2);
            if (strchr(d->obshap[0], 'N') == NULL) {
                assert(strchr(d->obshap[1], 'N') == NULL);
                hap = d->obshap[0];
                ++cloc->alleles[move(hap)];
                hap = d->obshap[1];
                ++cloc->alleles[move(hap)];
            }
        }
        cloc->populate_alleles();
    }

} catch (exception&) {
    cerr << "Error: Locus " << fasta_record.id << "\n"
         << "Error: Bad GStacks files.\n";
    throw;
}}

template<class LocusT> void
PopMap<LocusT>::populate_external(
    LocusT* cloc,
    Datum** locdata,
    int cloc_id,
    const VcfRecord& rec,
    const VcfHeader& header,
    const MetaPopInfo* mpopi
) {
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

    for (size_t s = 0; s < mpopi->samples().size(); ++s) {
        size_t vcf_index = header.sample_index(mpopi->samples()[s].name);
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
        ++cloc->cnt;
        ++cloc->hcnt;

        // id, len, lnl
        d->id  = cloc->id;
        d->len = cloc->len;
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
