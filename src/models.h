// -*-mode:c++; c-style:k&r; c-basic-offset:4;-*-
//
// Copyright 2010-2012, Julian Catchen <jcatchen@uoregon.edu>
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

#ifndef __MODELS_H__
#define __MODELS_H__

#include "constants.h"
#include "utils.h"
#include "DNASeq4.h"
#include "stacks.h"
#include "locus.h"
#include "mstack.h"

//
// Possible models for calling nucleotide positions as fixed or variable
//
enum modelt {fixed, snp, bounded};

//
// For use with the multinomial model to call fixed nucleotides.
//
extern int barcode_size;
extern double barcode_err_freq;

extern double heterozygote_limit;
extern double homozygote_limit;
extern double bound_low;  // For the bounded-snp model.
extern double bound_high; // For the bounded-snp model.
extern double p_freq;     // For the fixed model.

bool lrtest(double lnl_althyp, double lnl_nullhyp, double threshold); // Where threshold is a value in the Chi2 distribution.

void     set_model_thresholds (double alpha);
snp_type call_snp (double l_ratio);
snp_type call_snp (double lnl_hom, double lnl_het);

double   lnl_multinomial_model_hom (double total, double n1);
double   lnl_multinomial_model_het (double total, double n1n2);
double   lr_multinomial_model         (double nuc_1, double nuc_2, double nuc_3, double nuc_4);
double   lr_bounded_multinomial_model (double nuc_1, double nuc_2, double nuc_3, double nuc_4);

void call_bounded_multinomial_snp(MergedStack *, int, map<char, int> &, bool);
void call_bounded_multinomial_snp(Locus *, int, map<char, int> &);
void call_multinomial_snp(MergedStack *, int, map<char, int> &, bool);
void call_multinomial_snp(Locus *, int, map<char, int> &);
void call_multinomial_fixed(MergedStack *, int, map<char, int> &);

double   heterozygous_likelihood(int, map<char, int> &);
double   homozygous_likelihood(int, map<char, int> &);

class SampleSiteData {
    Counts<Nt2> depths_;
    GtLiks lnls_;
    // The genotype call and the corresponding nucleotides.
    // hom {nt, Nt4::n} | het {min_nt, max_nt} | unk {Nt4::n, Nt4::n}
    // For hets, the two nucleotides are sorted lexically (A<C<G<T).
    snp_type call_;
    array<Nt4, 2> nts_;
public:
    SampleSiteData() : depths_(), lnls_(), call_(snp_type_unk), nts_{Nt4::n,Nt4::n} {}

    const Counts<Nt2>& depths() const {return depths_;}
          Counts<Nt2>& depths()       {return depths_;}
    const GtLiks& lnls() const {return lnls_;}
          GtLiks& lnls()       {return lnls_;}

    bool has_coverage() const {return depths_.sum() != 0;}

    snp_type call() const {return call_;}
    Nt4 nt0() const {assert(call_==snp_type_hom || call_==snp_type_het); return nts_[0];}
    Nt4 nt1() const {assert(call_==snp_type_het); return nts_[1];}

    void add_call(snp_type c, Nt4 rank0_nt, Nt4 rank1_nt);
};

class SiteCall {
    Counts<Nt2> tot_depths_;
    map<Nt4, size_t> alleles_;
    vector<SampleSiteData> sample_data_;
public:
    SiteCall(const Counts<Nt2>& tot_depths, map<Nt4, size_t>&& alleles, vector<SampleSiteData>&& sample_data)
        : tot_depths_(tot_depths), alleles_(move(alleles)), sample_data_(move(sample_data))
        {}
    const Counts<Nt2>& tot_depths() const {return tot_depths_;}
    size_t tot_depth() const {return tot_depths_.sum();}
    const map<Nt4, size_t>& alleles() const {return alleles_;}
    const vector<SampleSiteData>& sample_data() const {return sample_data_;}

    static map<Nt4,size_t> tally_allele_freqs(const vector<SampleSiteData>& spldata);
};

class Model {
public:
    virtual ~Model() {}
    virtual SiteCall call(const CLocAlnSet::site_iterator& site) const = 0;
};

//
// MultinomialModel: the standard Stacks v.1 model described in Hohenloe2010.
//
class MultinomialModel : public Model {
public:
    SiteCall call(const CLocAlnSet::site_iterator& site) const;
};

//
// MarukiHighModel: the model of Maruki & Lynch (2017) for high-coverage data.
//
class MarukiHighModel : public Model {
    double calc_hom_lnl(double n, double n1) const;
    double calc_het_lnl(double n, double n1n2) const;
public:
    SiteCall call(const CLocAlnSet::site_iterator& site) const;
};

//
// ==================
// Inline definitions
// ==================
//

inline
bool lrtest(double lnl_althyp, double lnl_nullhyp, double threshold) {
    return 2.0 * (lnl_althyp - lnl_nullhyp) > threshold;
}

inline
snp_type call_snp (double l_ratio) {
    if (l_ratio <= heterozygote_limit)
        return snp_type_het;
    else if (l_ratio >= homozygote_limit)
        return snp_type_hom;
    else
        return snp_type_unk;
}

inline
snp_type call_snp (double lnl_hom, double lnl_het) {
    return call_snp(2.0 * (lnl_hom - lnl_het));
}

inline
double lnl_multinomial_model_hom (double total, double n1) {
    if (n1 == total)
        return 0.0;
    else if (n1 < 0.25 * total)
        return total * log(0.25); // With epsilon estimate bounded at 1.0
    else
        return n1 * log(n1/total) + (total-n1) * log( (total-n1)/(3.0*total) );
}

inline
double lnl_multinomial_model_het (double total, double n1n2) {
    if (n1n2 == total)
        return total * log(0.5);
    else if (n1n2 < 0.5 * total)
        return total * log(0.25); // With epsilon estimate bounded at 1.0
    else
        return n1n2 * log( n1n2/(2.0*total) ) + (total-n1n2) * log( (total-n1n2)/(2.0*total) );
}

inline
double lr_multinomial_model_legacy (double nuc_1, double nuc_2, double nuc_3, double nuc_4) {
    //
    // This function is to check that the refactored function gives the same
    // results as the original code (i.e. this code).
    //

    double total = nuc_1 + nuc_2 + nuc_3 + nuc_4;
    assert(total > 0.0);

    double l_ratio = (nuc_1 * log(nuc_1 / total));

    if (total - nuc_1 > 0.0)
        l_ratio += ((total - nuc_1) * log((total - nuc_1) / (3.0 * total)));

    if (nuc_1 + nuc_2 > 0.0)
        l_ratio -= ((nuc_1 + nuc_2) * log((nuc_1 + nuc_2) / (2.0 * total)));

    if (nuc_3 + nuc_4 > 0.0)
        l_ratio -= ((nuc_3 + nuc_4) * log((nuc_3 + nuc_4) / (2.0 * total)));

    l_ratio *= 2.0;

    return l_ratio;
}

inline
double lr_multinomial_model (double nuc_1, double nuc_2, double nuc_3, double nuc_4) {
    //
    // Method of Paul Hohenlohe <hohenlohe@uidaho.edu>, personal communication.
    //
    // For a diploid individual, there are ten possible genotypes
    // (four homozygous and six heterozygous genotypes).  We calculate
    // the likelihood of each possible genotype by using a multinomial
    // sampling distribution, which gives the probability of observing
    // a set of read counts (n1,n2,n3,n4) given a particular genotype.
    //

    double total = nuc_1 + nuc_2 + nuc_3 + nuc_4;
    assert(total > 0.0);

    double l_ratio = 2.0 * (lnl_multinomial_model_hom(total, nuc_1) - lnl_multinomial_model_het(total, nuc_1+nuc_2));

    assert(nuc_1+nuc_2 == 0 || almost_equal(l_ratio, lr_multinomial_model_legacy(nuc_1,nuc_2,nuc_3,nuc_4)));
    return l_ratio;
}

inline
double lr_bounded_multinomial_model (double nuc_1, double nuc_2, double nuc_3, double nuc_4) {

    //
    // Method of Paul Hohenlohe <hohenlohe@uidaho.edu>, personal communication.
    //

    double total = nuc_1 + nuc_2 + nuc_3 + nuc_4;
    assert(total > 0.0);

    //
    // Calculate the site specific error rate for homozygous and heterozygous genotypes.
    //
    double epsilon_hom  = (4.0 / 3.0) * ((total - nuc_1) / total);
    double epsilon_het  = 2.0 * ((nuc_3 + nuc_4) / total);

    //
    // Check if the error rate is above or below the specified bound.
    //
    if (epsilon_hom < bound_low)
        epsilon_hom = bound_low;
    else if (epsilon_hom > bound_high)
        epsilon_hom = bound_high;

    if (epsilon_het < bound_low)
        epsilon_het = bound_low;
    else if (epsilon_het > bound_high)
        epsilon_het = bound_high;

    //
    // Calculate the log likelihood for the homozygous and heterozygous genotypes.
    //
    double ln_L_hom = nuc_1 * log(1 - ((3.0/4.0) * epsilon_hom));
    ln_L_hom += epsilon_hom > 0.0 ? ((nuc_2 + nuc_3 + nuc_4) * log(epsilon_hom / 4.0)) : 0.0;

    double ln_L_het = (nuc_1 + nuc_2) * log(0.5 - (epsilon_het / 4.0));
    ln_L_het += epsilon_het > 0.0 ? ((nuc_3 + nuc_4) * log(epsilon_het / 4.0)) : 0.0;

    //
    // Calculate the likelihood ratio.
    //
    double l_ratio  = 2.0 * (ln_L_hom - ln_L_het);

    // cerr << "  Nuc_1: " << nuc_1 << " Nuc_2: " << nuc_2 << " Nuc_3: " << nuc_3 << " Nuc_4: " << nuc_4
    //   << " epsilon homozygote: " << epsilon_hom
    //   << " epsilon heterozygote: " << epsilon_het
    //   << " Log likelihood hom: " << ln_L_hom
    //   << " Log likelihood het: " << ln_L_het
    //   << " Likelihood ratio: " << l_ratio << "\n";

    return l_ratio;
}

inline
void SampleSiteData::add_call(snp_type c, Nt4 rank0_nt, Nt4 rank1_nt) {
    call_ = c;
    if (call_ == snp_type_hom) {
        nts_[0] = rank0_nt;
    } else if (call_ == snp_type_het) {
        if (rank0_nt < rank1_nt)
            nts_ = {rank0_nt, rank1_nt};
        else
            nts_ = {rank1_nt, rank0_nt};
    }
}

#endif // __MODELS_H__
