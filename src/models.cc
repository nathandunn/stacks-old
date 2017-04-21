// -*-mode:c++; c-style:k&r; c-basic-offset:4;-*-
//
// Copyright 2010 - 2012, Julian Catchen <jcatchen@uoregon.edu>
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
// models.cc -- routines to detect polymorphism (snp) and detect a lack of polymorphism (fixed).
//
// Julian Catchen
// jcatchen@uoregon.edu
// University of Oregon
//
// $Id$
//
#include "models.h"

int    barcode_size     = 5;
double barcode_err_freq = 0.0;

double heterozygote_limit = -3.84;
double homozygote_limit   = 3.84;
double bound_low          = 0.0;
double bound_high         = 1.0;
double p_freq             = 0.5;

vector<pair<char, int>> sort_acgt(const map<char, int>&);
void record_snp(SNP& snp, snp_type type, uint col, double l_ratio, const vector<pair<char, int>>& nuc);
void record_dummy_snp(SNP& snp, uint col);

bool set_model_type(modelt& model, const string& arg) {
    if (arg == "snp")
        model = ::snp;
    else if (arg == "fixed")
        model = ::fixed;
    else if (arg == "bounded")
        model = ::bounded;
    else if (arg == "marukihigh")
        model = ::marukihigh;
    else if (arg == "marukilow")
        model = ::marukilow;
    else
        return false;

    return true;
}

void set_model_thresholds(double alpha) {
    if (alpha == 0.1) {
        heterozygote_limit = -2.71;
        homozygote_limit   =  2.71;
    } else if (alpha == 0.05) {
        heterozygote_limit = -3.84;
        homozygote_limit   =  3.84;
    } else if (alpha == 0.01) {
        heterozygote_limit = -6.64;
        homozygote_limit   =  6.64;
    } else if (alpha == 0.001) {
        heterozygote_limit = -10.83;
        homozygote_limit   =  10.83;
    } else {
        cerr << "Error: Unsupported alpha value '" << alpha << "'.\n";
        throw exception();
    }
}

vector<pair<char, int>> sort_acgt(const map<char, int>& counts) {
    vector<pair<char, int>> nuc;
    for (map<char, int>::const_iterator i=counts.begin(); i!=counts.end(); i++)
        if (i->first != 'N')
            nuc.push_back(make_pair(i->first, i->second));
    sort(nuc.begin(), nuc.end(), compare_pair);
    return nuc;
}

void record_snp(SNP& snp, snp_type type, uint col, double l_ratio, const vector<pair<char, int>>& nuc) {

    snp.type   = type;
    snp.col    = col;
    snp.lratio = l_ratio;
    snp.rank_1 = nuc[0].first;

    switch (type) {
    case snp_type_het:
        snp.rank_2 = nuc[1].first;
        break;
    case snp_type_hom:
        snp.rank_2 = '-';
        break;
    default:
        // snp_type_unk, with at least one observation (otherwise record_dummy_snp
        // would have been called)
        // A rank_2 nucleotide is set only if at least two nucleotides were observed.
        snp.rank_2 = nuc[1].second > 0 ? nuc[1].first : '-';
        break;
    }

    snp.rank_3 = 0;
    snp.rank_4 = 0;
}

void record_dummy_snp(SNP& snp, uint col) {
    snp.type   = snp_type_unk;
    snp.col    = col;
    snp.lratio = 0.0;
    snp.rank_1 = 'N';
    snp.rank_2 = '-';
    snp.rank_3 = 0;
    snp.rank_4 = 0;
}

void call_multinomial_snp(MergedStack *tag, int col, map<char, int> &n, bool record_snps) {
    vector<pair<char, int> > nuc = sort_acgt(n);
    if (nuc[0].second == 0) {
        if (record_snps) {
            tag->snps.push_back(new SNP());
            record_dummy_snp(*tag->snps.back(), col);
        }
        return;
    }

    double l_ratio = lr_multinomial_model(nuc[0].second, nuc[1].second, nuc[2].second, nuc[3].second);
    snp_type type = call_snp(l_ratio);
    if (record_snps) {
        tag->snps.push_back(new SNP());
        record_snp(*tag->snps.back(), type, col, l_ratio, nuc);
    }
}

void call_bounded_multinomial_snp(MergedStack *tag, int col, map<char, int> &n, bool record_snps) {
    vector<pair<char, int> > nuc = sort_acgt(n);
    if (nuc[0].second == 0) {
        if (record_snps) {
            tag->snps.push_back(new SNP());
            record_dummy_snp(*tag->snps.back(), col);
        }
        return;
    }

    double l_ratio = lr_bounded_multinomial_model(nuc[0].second, nuc[1].second, nuc[2].second, nuc[3].second);
    snp_type type = call_snp(l_ratio);
    if (record_snps) {
        tag->snps.push_back(new SNP());
        record_snp(*tag->snps.back(), type, col, l_ratio, nuc);
    }
}

void call_multinomial_snp(Locus *tag, int col, map<char, int> &n) {
    vector<pair<char, int> > nuc = sort_acgt(n);
    if (nuc[0].second == 0) {
        record_dummy_snp(*tag->snps[col], col);
        return;
    }

    double l_ratio = lr_multinomial_model(nuc[0].second, nuc[1].second, nuc[2].second, nuc[3].second);
    snp_type type = call_snp(l_ratio);
    record_snp(*tag->snps[col], type, col, l_ratio, nuc);
}

void call_bounded_multinomial_snp(Locus *tag, int col, map<char, int> &n) {
    vector<pair<char, int> > nuc = sort_acgt(n);
    if (nuc[0].second == 0) {
        record_dummy_snp(*tag->snps[col], col);
        return;
    }

    double l_ratio = lr_bounded_multinomial_model(nuc[0].second, nuc[1].second, nuc[2].second, nuc[3].second);
    snp_type type = call_snp(l_ratio);
    record_snp(*tag->snps[col], type, col, l_ratio, nuc);
}

void call_multinomial_fixed (MergedStack *tag, int col, map<char, int> &n) {
    const double nucleotide_fixed_limit = 1.92;

    vector<pair<char, int> > nuc = sort_acgt(n);
    if (nuc[0].second == 0) {
        tag->snps.push_back(new SNP());
        record_dummy_snp(*tag->snps.back(), col);
        return;
    }

    double l_ratio;
    double nuc_1 = nuc[0].second;
    double nuc_2 = nuc[1].second;
    {
        //
        // Method of Paul Hohenlohe <hohenlo@uoregon.edu>, personal communication.
        //
        // Each population sample contains DNA from 6 individuals, so a
        // sample of 12 alleles from the population. We want to assign a
        // nucleotide (A,C,G,T) to each position where the population is
        // fixed or nearly so, and N to each position that is either
        // polymorphic within the population or has insufficient coverage
        // depth to make a call. We can do this with a likelihood ratio
        // test of the read counts, testing whether the allele frequency
        // of the dominant allele is significantly larger than some
        // threshold p) , stepping through each nucleotide position across
        // RAD tags.
        //

        double epsilon = -1 * (log(1 - barcode_err_freq) / barcode_size);

        l_ratio  =
            nuc_1 * log( ((4 * nuc_1 * (1 - epsilon)) + ((nuc_1 + nuc_2) * epsilon)) /
                         ((4 * p_freq * (nuc_1 + nuc_2) * (1 - epsilon)) + ((nuc_1 + nuc_2) * epsilon)) );

        l_ratio +=
            nuc_2 * log( ((4 * nuc_2 * (1 - epsilon)) + ((nuc_1 + nuc_2) * epsilon)) /
                         ((4 * (1 - p_freq) * (nuc_1 + nuc_2) * (1 - epsilon)) + ((nuc_1 + nuc_2) * epsilon)) );

        //cerr << "Nuc_1: " << nuc_1 << " Nuc_2: " << nuc_2 << " Likelihood ratio: " << l_ratio << "\n";
    }

    double n_ratio = nuc_1 / (nuc_1 + nuc_2);

    snp_type type = n_ratio < p_freq || l_ratio < nucleotide_fixed_limit ? snp_type_unk : snp_type_hom;

    tag->snps.push_back(new SNP());
    record_snp(*tag->snps.back(), type, col, l_ratio, nuc);
}

//
// ln L(1/2) = ln(n! / n_1!n_2!n_3!n_4!) +
//               (n_1 + n_2) * ln(n_1 + n_2 / 2n) +
//               (n_3 + n_4) * ln(n_3 + n_4 / 2n)
//
double
heterozygous_likelihood(int col, map<char, int> &nuc)
{
    vector<pair<char, int> > cnts;
    map<char, int>::iterator i;

    double n = 0;
    for (i = nuc.begin(); i != nuc.end(); i++) {
        n += i->second;
        cnts.push_back(make_pair(i->first, i->second));
    }

    sort(cnts.begin(), cnts.end(), compare_pair);

    double n_1 = cnts[0].second;
    double n_2 = cnts[1].second;
    double n_3 = cnts[2].second;
    double n_4 = cnts[3].second;

    double term_1 =
        reduced_log_factorial(n, n_1) -
        (log_factorial(n_2) + log_factorial(n_3) + log_factorial(n_4));

    double term_3 = (n_3 + n_4 > 0) ? log((n_3 + n_4) / (2 * n)) : 0;

    double lnl =
        term_1 +
        ((n_1 + n_2) * log((n_1 + n_2) / (2 * n))) +
        ((n_3 + n_4) * term_3);

    return lnl;
}

//
// ln L(1/1) = ln(n! / n_1!n_2!n_3!n_4!) +
//               n_1 * ln(n_1 / n) +
//               (n - n_1) * ln(n - n_1 / 3n)
//
double
homozygous_likelihood(int col, map<char, int> &nuc)
{
    vector<pair<char, int> > cnts;
    map<char, int>::iterator i;

    double n = 0;
    for (i = nuc.begin(); i != nuc.end(); i++) {
        n += i->second;
        cnts.push_back(make_pair(i->first, i->second));
    }

    sort(cnts.begin(), cnts.end(), compare_pair);

    double n_1 = cnts[0].second;
    double n_2 = cnts[1].second;
    double n_3 = cnts[2].second;
    double n_4 = cnts[3].second;

    double term_1 =
        reduced_log_factorial(n, n_1) -
        (log_factorial(n_2) + log_factorial(n_3) + log_factorial(n_4));

    double term_3 = n - n_1 > 0 ? log((n - n_1) / (3 * n)) : 0;

    double lnl =
        term_1 +
        (n_1 * log(n_1 / n)) +
        ((n - n_1) * term_3);

    return lnl;
}

Nt2 SiteCall::most_frequent_allele() const {
    assert(!alleles().empty());
    auto a = alleles().begin();
    auto best = a;
    for(; a!=alleles().end(); ++a)
        if (a->second > best->second)
            best = a;
    return best->first;
}

map<Nt2,size_t> SiteCall::tally_allele_freqs(const vector<SampleCall>& spldata) {
    //
    // Iterate over the SampleCall's & record the existing alleles and their
    // frequencies.
    //

    map<Nt2,size_t> allele_freqs;

    Counts<Nt2> counts;
    for (const SampleCall& sd : spldata) {
        switch (sd.call()) {
        case snp_type_hom :
            counts.increment(sd.nt0());
            counts.increment(sd.nt0());
            break;
        case snp_type_het :
            counts.increment(sd.nt0());
            counts.increment(sd.nt1());
            break;
        default:
            // snp_type_unk
            break;
        }
    }

    array<pair<size_t,Nt2>,4> sorted_alleles = counts.sorted();
    size_t i = 0;
    while (i<4 && sorted_alleles[i].first != 0) {
        allele_freqs.insert({sorted_alleles[i].second, sorted_alleles[i].first});
        ++i;
    }

    return allele_freqs;
}

SiteCall MultinomialModel::call(vector<Counts<Nt2>>&& sample_depths) const {

    size_t n_samples = sample_depths.size();

    Counts<Nt2> tot_depths;
    for (const Counts<Nt2>& depths : sample_depths)
        tot_depths += depths;

    //
    // Make genotype calls.
    //
    vector<SampleCall> sample_calls (n_samples);
    array<pair<size_t,Nt2>,4> sorted;
    for (size_t sample=0; sample<n_samples; ++sample) {
        const Counts<Nt2>& depths = sample_depths[sample];
        size_t dp = depths.sum();
        if (dp == 0)
            continue;

        sorted = depths.sorted();
        SampleCall& c = sample_calls[sample];
        double lnl_hom = lnl_multinomial_model_hom(dp, sorted[0].first);
        double lnl_het = lnl_multinomial_model_het(dp, sorted[0].first+sorted[1].first);
        c.lnls().set(sorted[0].second, sorted[0].second, lnl_hom);
        c.lnls().set(sorted[0].second, sorted[1].second, lnl_het);

        snp_type gt_call = call_snp(lnl_hom, lnl_het);
        c.set_call(gt_call, sorted[0].second, sorted[1].second);
    }

    //
    // Record the existing alleles and their frequencies.
    //
    map<Nt2,size_t> allele_freqs = SiteCall::tally_allele_freqs(sample_calls);

    //
    // Compute the missing genotype likelihoods, if any.
    //
    if (allele_freqs.size() < 2) {
        sample_calls = vector<SampleCall>();
    } else {
        for (size_t sample=0; sample<n_samples; ++sample) {
            const Counts<Nt2>& depths = sample_depths[sample];
            size_t dp = depths.sum();
            if (dp == 0)
                continue;

            SampleCall& c = sample_calls[sample];
            for (auto nt1_it=allele_freqs.begin(); nt1_it!=allele_freqs.end(); ++nt1_it) {
                Nt2 nt1 (nt1_it->first);
                // Homozygote.
                if (!c.lnls().has_lik(nt1, nt1)) {
                    double lnl = lnl_multinomial_model_hom(dp, depths[nt1]);
                    c.lnls().set(nt1, nt1, lnl);
                }
                // Heterozygote(s).
                auto nt2_it = nt1_it;
                ++nt2_it;
                for (; nt2_it!=allele_freqs.end(); ++nt2_it) {
                    Nt2 nt2 (nt2_it->first);
                    if (!c.lnls().has_lik(nt1, nt2)) {
                        double lnl = lnl_multinomial_model_het(dp, depths[nt1]+depths[nt2]);
                        c.lnls().set(nt1, nt2, lnl);
                    }
                }
            }
        }
    }

    return SiteCall(tot_depths, move(sample_depths), move(allele_freqs), move(sample_calls));
}

double MarukiHighModel::calc_hom_lnl(double n, double n1) const {
    // This returns the same value as the Hohenlohe ('snp/binomial') model except
    // when the error rate estimate is bounded at 1.0 (i.e. `n1 < 0.25*n` c.f.
    // Hohenlohe equations).
    if (n1 == n)
        return 0.0;
    else if (n1 == 0.0)
        return n * log(1.0/3.0);
    else
        return n1 * log(n1/n) + (n-n1) * log( (n-n1)/(3.0*n) );
}

double MarukiHighModel::calc_het_lnl(double n, double n1n2) const {
    // This returns the same value as the Hohenlohe ('snp/binomial') model except
    // when the error rate estimate is bounded at 1.0 (i.e. `n1n2 < 0.5*n` c.f.
    // Hohenlohe equations).
    if (n1n2 == n)
        return n * log(0.5);
    else if (n1n2 < (1.0/3.0) * n)
        return n1n2 * log(1.0/6.0) + (n-n1n2) * log(1.0/3.0);
    else
        return n1n2 * log( n1n2/(2.0*n) ) + (n-n1n2) * log(1.0 - n1n2/(2.0*n) );
}

SiteCall MarukiHighModel::call(vector<Counts<Nt2>>&& sample_depths) const {

    /*
     * For this model the procedure is:
     * I. Obtain the most commonly seen nucleotide M.
     * II. Look for alternative alleles, if any:
     *     For each sample:
     *         If there is a genotype significantly better than MM:
     *             (This genotype can be Mm, mm or mn.)
     *             The site is polymorphic.
     *             Record m as an alternative allele.
     *             If the genotype is mn AND is significantly better than mm:
     *                 Record n as an allele.
     * III. Given the known alleles, compute the likelihoods for all possible
     *      genotypes (n.b. most of the non-trivial ones have already been
     *      computed), and call genotypes.
     */

    const size_t n_samples = sample_depths.size();

    Counts<Nt2> tot_depths;
    for (const Counts<Nt2>& depths : sample_depths)
        tot_depths += depths;

    if (tot_depths.sum() == 0)
        return SiteCall(tot_depths, move(sample_depths), map<Nt2,size_t>(), vector<SampleCall>());

    //
    // I.
    // Count the observed nucleotides of the site for all samples; set
    // `SampleCall::depths_`.
    // Then find the most common nucleotide across the population.
    //
    Nt2 nt_ref = tot_depths.sorted()[0].second;

    //
    // II.
    // Look for alternative alleles; start filling SampleCall::lnls_.
    //
    set<Nt2> alleles;
    vector<SampleCall> sample_calls (n_samples);
    alleles.insert(nt_ref); // Note: There are limit cases (esp. low coverage over many
                            // individuals) where ref_nt would not appear in any of the
                            // significant genotypes; we record it as an allele anyway.
    array<pair<size_t,Nt2>,4> sorted;
    for (size_t sample=0; sample<n_samples; ++sample) {
        const Counts<Nt2>& depths = sample_depths[sample];
        size_t dp = depths.sum();
        if (dp == 0)
            continue;

        // Find the best genotype for the sample -- this is either the homzygote
        // for the rank0 nucleotide or the heterozygote for the rank0 and rank1
        // nucleotides.
        sorted = depths.sorted();
        Nt2 nt0 = sorted[0].second;
        Nt2 nt1 = sorted[1].second;
        size_t dp0 = sorted[0].first;
        size_t dp1 = sorted[1].first;

        if (dp1 == 0 && nt0 == nt_ref)
            // We can wait until we know whether the site is fixed.
            continue;

        SampleCall& c = sample_calls[sample];
        double lnl_hom = calc_hom_lnl(dp, dp0);
        double lnl_het = calc_het_lnl(dp, dp0+dp1);
        c.lnls().set(nt0, nt0, lnl_hom);
        c.lnls().set(nt0, nt1, lnl_het);

        // Compare this likelihood to that of the ref,ref homozygote.
        static const double& polymorphism_signif_thr = homozygote_limit; //xxx Hack.
        if (nt0 == nt_ref) {
            if (lrtest(lnl_het, lnl_hom, polymorphism_signif_thr))
                // Record the alternative allele.
                alleles.insert(nt1);
        } else {
            double lnl_ref = calc_hom_lnl(dp, depths[nt_ref]);
            c.lnls().set(nt_ref, nt_ref, lnl_ref);
            double lnl_best = std::max(lnl_hom, lnl_het);
            if (lrtest(lnl_best, lnl_ref, polymorphism_signif_thr)) {
                // Record one alternative allele.
                alleles.insert(nt0);
                if (nt1!=nt_ref && lrtest(lnl_het, lnl_hom, polymorphism_signif_thr))
                    // Record a second alternative allele (the SNP is at least ternary).
                    alleles.insert(nt1);
            }
        }
    }

    //
    // III.
    // Compute the likelihoods for all possible genotypes & call genotypes.
    //
    if (alleles.size() == 1) {
        sample_calls = vector<SampleCall>();
    } else {
        for (size_t sample=0; sample<n_samples; ++sample) {
            const Counts<Nt2>& depths = sample_depths[sample];
            size_t dp = depths.sum();
            if (dp == 0)
                continue;

            SampleCall& c = sample_calls[sample];
            for (auto nt1=alleles.begin(); nt1!=alleles.end(); ++nt1) {
                // Homozygote.
                if (!c.lnls().has_lik(*nt1, *nt1))
                    c.lnls().set(*nt1, *nt1, calc_hom_lnl(dp, depths[*nt1]));
                // Heterozygote(s).
                auto nt2 = nt1;
                ++nt2;
                for (; nt2!=alleles.end(); ++nt2)
                    if (!c.lnls().has_lik(*nt1, *nt2))
                        c.lnls().set(*nt1, *nt2, calc_het_lnl(dp, depths[*nt1]+depths[*nt2]));
            }

            // Call the genotype -- skiping ignored alleles.
            sorted = depths.sorted();
            auto nt0 = sorted.begin();
            while(nt0 != sorted.end() && !alleles.count(nt0->second))
                ++nt0;
            auto nt1 = nt0;
            ++nt1;
            while(nt1 != sorted.end() && !alleles.count(nt1->second))
                ++nt1;
            if (nt1 != sorted.end()) {
                double lnl_hom = c.lnls().at(nt0->second, nt0->second);
                double lnl_het = c.lnls().at(nt0->second, nt1->second);
                c.set_call(call_snp(lnl_hom, lnl_het), nt0->second, nt1->second);
            }
        }
    }

    //
    // Finally, compute allele frequencies.
    //
    map<Nt2,size_t> allele_freqs;
    if (alleles.size() == 1)
        allele_freqs.insert({*alleles.begin(),-1});
    else
        allele_freqs = SiteCall::tally_allele_freqs(sample_calls);
    assert(allele_freqs.size() == alleles.size());

    return SiteCall(tot_depths, move(sample_depths), move(allele_freqs), move(sample_calls));
}
