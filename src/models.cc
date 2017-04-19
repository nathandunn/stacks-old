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

SiteCall MultinomialModel::call(const CLocAlnSet::site_iterator& site) const {

    //
    // Look at this site in each sample; make genotype calls.
    //
    size_t tot_depth = 0;
    vector<SampleCall> sample_calls;
    sample_calls.reserve(site.mpopi().samples().size());
    Counts<Nt4> counts;
    array<pair<size_t,Nt4>,4> sorted;
    for (size_t s=0; s<site.mpopi().samples().size(); ++s) {
        site.counts(counts, s);
        sorted = counts.sorted();

        if (sorted[0].first == 0) {
            sample_calls.push_back(SampleCall());
        } else {
            snp_type gt_call = call_snp(lr_multinomial_model(
                    sorted[0].first,
                    sorted[1].first,
                    sorted[2].first,
                    sorted[3].first
                    ));
            sample_calls.push_back(SampleCall(
                    counts,
                    gt_call,
                    sorted[0].second,
                    sorted[1].second
                    ));
            tot_depth += sample_calls.back().depths().sum();
        }
    }

    //
    // Iterate over the SampleCalls & record the genotypes that were found.
    //
    map<Nt4, size_t> alleles;
    counts.clear();
    for (const SampleCall& sc : sample_calls) {
        switch (sc.call()) {
        case snp_type_hom :
            counts.increment(sc.nt0());
            counts.increment(sc.nt0());
            break;
        case snp_type_het :
            counts.increment(sc.nt0());
            counts.increment(sc.nt1());
            break;
        default:
            // snp_type_unk
            break;
        }
    }
    sorted = counts.sorted();
    if (sorted[0].first > 0) {
        // At least one allele was observed.
        alleles.insert({sorted[0].second, sorted[0].first});

        if (sorted[1].first > 0) {
            // SNP with at least two alleles.
            alleles.insert({sorted[0].second, sorted[0].first});

            if (sorted[2].first > 0) {
                alleles.insert({sorted[0].second, sorted[0].first});

                if (sorted[3].first > 0) {
                    // Quaternary SNP.
                    alleles.insert({sorted[0].second, sorted[0].first});
                }
            }
        }
    }

    return SiteCall(tot_depth, move(alleles), move(sample_calls));
}
