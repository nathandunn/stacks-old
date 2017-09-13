// -*-mode:c++; c-style:k&r; c-basic-offset:4;-*-
//
// Copyright 2011-2017, Julian Catchen <jcatchen@illinois.edu>
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

#include "PopSum.h"

LocPopSum::LocPopSum(size_t cloc_len, const MetaPopInfo& mpopi)
{
    this->_pop_cnt  = mpopi.pops().size();
    this->_meta_pop = new LocTally(cloc_len);
    this->_per_pop  = new LocSum * [this->_pop_cnt];

    for (size_t i = 0; i < mpopi.pops().size(); i++)
        this->_per_pop[i] = new LocSum(cloc_len);
}

LocPopSum::~LocPopSum()
{
    delete [] this->_per_pop;
    delete this->_meta_pop;
}

int
LocPopSum::sum_pops(const CSLocus *cloc, const Datum **d, const MetaPopInfo &mpopi,
                    bool verbose, ofstream &log_fh)
{
    uint len = strlen(cloc->con);
    int res;
    set<int> snp_cols;
    int incompatible_loci = 0;

    for (uint i = 0; i < this->_pop_cnt; i++) {
        LocSum *s = this->_per_pop[i];

        const Pop& pop = mpopi.pops().at(i);
        //
        // Check if this locus has already been filtered and is NULL in all individuals.
        //
        bool filtered = true;
        for (uint k = pop.first_sample; k <= pop.last_sample; k++) {
            if (d[k] != NULL)
                filtered = false;
        }
        if (filtered == true) {
            for (uint k = 0; k < len; k++) {
                s->nucs[k].filtered_site = true;
            }
            continue;
        }

        //
        // The catalog records which nucleotides are heterozygous. For these nucleotides we will
        // calculate observed genotype frequencies, allele frequencies, and expected genotype frequencies.
        //
        for (uint k = 0; k < cloc->snps.size(); k++) {
            res = this->tally_heterozygous_pos(cloc, d, s,
                                               cloc->snps[k]->col, k, pop.first_sample, pop.last_sample);
            //
            // If site is incompatible (too many alleles present), log it.
            //
            if (res < 0) {
                s->nucs[cloc->snps[k]->col].incompatible_site = true;

                incompatible_loci++;
                if (verbose)
                    log_fh << "within_population\t"
                           << "incompatible_locus\t"
                           << cloc->id << "\t"
                           << cloc->loc.chr() << "\t"
                           << cloc->sort_bp(cloc->snps[k]->col) +1 << "\t"
                           << cloc->snps[k]->col << "\t"
                           << pop.name << "\n";
            }

            snp_cols.insert(cloc->snps[k]->col);
        }
        //
        // For all other fixed sites, we just need to record them.
        //
        for (uint k = 0; k < len; k++) {
            if (snp_cols.count(k)) continue;
            this->tally_fixed_pos(cloc, d, s, k, pop.first_sample, pop.last_sample);
        }

        snp_cols.clear();
    }

    return 0;
}

int
LocPopSum::tally_fixed_pos(const CSLocus *cloc, const Datum **d, LocSum *s,
                           int pos, uint start, uint end)
{
    double num_indv = 0.0;
    char   p_nuc    = 0;

    for (uint i = start; i <= end; i++) {
        if (d[i] == NULL || pos >= d[i]->len) continue;
        //
        // Before counting this individual, make sure the model definitively called this
        // position as hEterozygous or hOmozygous.
        //
        if (d[i]->model[pos] == 'E') {
            cerr << "Model: " << d[i]->model << "\n";
            cerr << "Warning: heterozygous model call at fixed nucleotide position: "
                 << "locus " << cloc->id << " individual " << d[i]->id << "; position: " << pos << "\n";
        }
        num_indv++;
        p_nuc = cloc->con[pos];
    }
    //
    // Record the results in the PopSum object.
    //
    s->nucs[pos].loc_id   = cloc->id;
    s->nucs[pos].bp       = cloc->sort_bp(pos);
    s->nucs[pos].fixed    = true;
    s->nucs[pos].num_indv = num_indv;
    s->nucs[pos].alleles  = 2 * num_indv;

    if (num_indv > 0) {
        s->nucs[pos].p        =  1.0;
        s->nucs[pos].p_nuc    =  p_nuc;
        s->nucs[pos].obs_hom  =  1.0;
        s->nucs[pos].obs_het  =  0.0;
        s->nucs[pos].exp_hom  =  1.0;
        s->nucs[pos].exp_het  =  0.0;
        s->nucs[pos].stat[0]  =  0.0; // pi
        s->nucs[pos].stat[1]  = -7.0; // fis
    }

    return 0;
}

int
LocPopSum::tally_heterozygous_pos(const CSLocus *cloc, const Datum **d, LocSum *s,
                                  int pos, int snp_index, uint start, uint end)
{
    //
    // Tally up the genotype frequencies.
    //
    int  nucs[4] = {0};
    uint i;
    char nuc;

    //cerr << "  Calculating summary stats at het locus " << cloc->id << " position " << pos << "; snp_index: " << snp_index << "\n";

    //
    // Iterate over each individual in this sub-population.
    //
    for (i = start; i <= end; i++) {
        if (d[i] == NULL || pos >= d[i]->len || d[i]->model[pos] == 'U') continue;

        //
        // Pull each allele for this SNP from the observed haplotype.
        //
        for (uint j = 0; j < d[i]->obshap.size(); j++) {
            nuc = d[i]->obshap[j][snp_index];

            switch(nuc) {
            case 'A':
            case 'a':
                nucs[0]++;
                break;
            case 'C':
            case 'c':
                nucs[1]++;
                break;
            case 'G':
            case 'g':
                nucs[2]++;
                break;
            case 'T':
            case 't':
                nucs[3]++;
                break;
            }
        }
    }

    //
    // Determine how many alleles are present at this position in this population.
    // We cannot deal with more than two alternative alleles, if there are more than two
    // in a single population, print a warning and exclude this nucleotide position.
    //
    int allele_cnt = 0;
    for (i = 0; i < 4; i++)
        if (nucs[i] > 0) allele_cnt++;

    if (allele_cnt > 2)
        return -1;

    //
    // Record which nucleotide is the P allele and which is the Q allele.
    //
    char p_allele = 0;
    char q_allele = 0;

    i = 0;
    while (p_allele == 0 && i < 4) {
        if (nucs[i] > 0) {
            switch(i) {
            case 0:
                p_allele = 'A';
                break;
            case 1:
                p_allele = 'C';
                break;
            case 2:
                p_allele = 'G';
                break;
            case 3:
                p_allele = 'T';
                break;
            }
        }
        i++;
    }
    while (q_allele == 0 && i < 4) {
        if (nucs[i] > 0) {
            switch(i) {
            case 1:
                q_allele = 'C';
                break;
            case 2:
                q_allele = 'G';
                break;
            case 3:
                q_allele = 'T';
                break;
            }
        }
        i++;
    }
    //cerr << "  P Allele: " << p_allele << "; Q Allele: " << q_allele << "\n";

    //
    // Calculate observed genotype frequencies.
    //
    double num_indv = 0.0;
    double obs_het  = 0.0;
    double obs_p    = 0.0;
    double obs_q    = 0.0;

    for (i = start; i <= end; i++) {
        if (d[i] == NULL || pos >= d[i]->len) continue;
        //
        // Before counting this individual, make sure the model definitively called this
        // position as hEterozygous or hOmozygous.
        //
        if (d[i]->model[pos] == 'E' || d[i]->model[pos] == 'O')
            num_indv++;
        else
            continue;

        if (d[i]->obshap.size() > 1 &&
            this->tally_observed_haplotypes(d[i]->obshap, snp_index) == 2)
                obs_het++;
        else if (d[i]->obshap[0][snp_index] == p_allele)
            obs_p++;
        else if (d[i]->obshap[0][snp_index] == q_allele)
            obs_q++;
    }
    //cerr << "  Num Individuals: " << num_indv << "; Obs Hets: " << obs_het << "; Obs P: " << obs_p << "; Obs Q: " << obs_q << "\n";

    if (num_indv == 0) return 0;

    //
    // Calculate total number of alleles
    //
    double tot_alleles = num_indv * 2;
    double allele_p    = obs_het + (2 * obs_p);
    double allele_q    = obs_het + (2 * obs_q);

    //
    // Calculate Pi, equivalent to expected heterozygosity (exp_het)
    //
    ////s->nucs[pos].stat[0] = this->pi(tot_alleles, allele_p, allele_q);

    if (s->nucs[pos].stat[0] == 0.0)
        s->nucs[pos].fixed = true;

    //
    // Convert to allele frequencies
    //
    allele_p = allele_p / tot_alleles;
    allele_q = allele_q / tot_alleles;

    //cerr << "  P allele frequency: " << allele_p << "; Q allele frequency: " << allele_q << "\n";

    // //
    // // If the minor allele frequency is below the cutoff, set it to zero.
    // //
    // if (minor_allele_freq > 0) {
    //  if (allele_p < allele_q) {
    //      if (allele_p < minor_allele_freq) {
    //          s->nucs[pos].pi            = 0.0;
    //          s->nucs[pos].fixed         = true;
    //          s->nucs[pos].filtered_site = true;
    //          return 0;
    //      }
    //  } else {
    //      if (allele_q < minor_allele_freq) {
    //          s->nucs[pos].pi            = 0.0;
    //          s->nucs[pos].fixed         = true;
    //          s->nucs[pos].filtered_site = true;
    //          return 0;
    //      }
    //  }
    // }

    //
    // Calculate expected genotype frequencies.
    //
    double exp_het = 2 * allele_p * allele_q; // 2pq
    // double exp_p   = allele_p * allele_p;     // p^2
    // double exp_q   = allele_q * allele_q;     // q^2

    //cerr << "  Expected Het: " << exp_het << "; Expected P: " << exp_p << "; Expected Q: " << exp_q << "\n";

    obs_het = obs_het / num_indv;
    obs_p   = obs_p   / num_indv;
    obs_q   = obs_q   / num_indv;

    //cerr << "  Obs Hets Freq: " << obs_het << "; Obs P Freq: " << obs_p << "; Obs Q Freq: " << obs_q << "\n";

    //
    // Record the results in the PopSum object.
    //
    s->nucs[pos].loc_id   = cloc->id;
    s->nucs[pos].bp       = cloc->sort_bp(pos);
    s->nucs[pos].num_indv = num_indv;
    s->nucs[pos].alleles  = tot_alleles;
    s->nucs[pos].p        = allele_p > allele_q ? allele_p : allele_q;
    s->nucs[pos].p_nuc    = allele_p > allele_q ? p_allele : q_allele;
    s->nucs[pos].q_nuc    = allele_p > allele_q ? q_allele : p_allele;
    s->nucs[pos].obs_hom  = 1 - obs_het;
    s->nucs[pos].obs_het  = obs_het;
    s->nucs[pos].exp_hom  = 1 - exp_het;
    s->nucs[pos].exp_het  = exp_het;

    //
    // Calculate F_is, the inbreeding coefficient of an individual (I) relative to the subpopulation (S):
    //   Fis = (exp_het - obs_het) / exp_het
    //
    double fis = s->nucs[pos].pi == 0 ? -7 : (s->nucs[pos].pi - obs_het) / s->nucs[pos].pi;

    s->nucs[pos].stat[1] = fis;

    return 0;
}

int
LocPopSum::tally_metapop(const CSLocus *cloc)
{
    int       variable_pop;
    uint16_t  p_cnt, q_cnt, col;
    LocSum   **s   = this->_per_pop;
    LocTally  *mp  = this->_meta_pop;
    uint16_t   len = strlen(cloc->con);

    for (col = 0; col < len; col++) {

        mp->nucs[col].col    = col;
        mp->nucs[col].bp     = cloc->sort_bp(col);
        mp->nucs[col].loc_id = cloc->id;

        this->tally_ref_alleles(col,
                                mp->nucs[col].allele_cnt,
                                mp->nucs[col].p_allele,
                                mp->nucs[col].q_allele,
                                p_cnt, q_cnt);

        //
        // Is this site variable?
        //
        if (mp->nucs[col].allele_cnt > 1)
            mp->nucs[col].fixed = false;

        for (uint j = 0; j < this->_pop_cnt; j++) {
            //
            // Sum the number of individuals examined at this locus across populations.
            //
            mp->nucs[col].num_indv += s[j]->nucs[col].num_indv;
            mp->nucs[col].pop_cnt  += s[j]->nucs[col].num_indv > 0 ? 1 : 0;
        }

        for (uint j = 0; j < this->_pop_cnt; j++) {
            //
            // Sum the most frequent allele across populations.
            //
            if (s[j]->nucs[col].p_nuc == mp->nucs[col].p_allele)
                mp->nucs[col].p_freq +=
                    s[j]->nucs[col].p * (s[j]->nucs[col].num_indv / (double) mp->nucs[col].num_indv);
            else
                mp->nucs[col].p_freq +=
                    (1 - s[j]->nucs[col].p) * (s[j]->nucs[col].num_indv / (double) mp->nucs[col].num_indv);
            //
            // Sum observed heterozygosity across populations.
            //
            mp->nucs[col].obs_het +=
                s[j]->nucs[col].obs_het * (s[j]->nucs[col].num_indv / (double) mp->nucs[col].num_indv);
        }

        //
        // We want to report the most frequent allele as the P allele. Reorder the alleles
        // if necessary.
        // XXX Possibly unstable for p_freq ~ 0.5. @Nick (July 2016)
        //
        if (mp->nucs[col].p_freq < 0.5) {
            char a = mp->nucs[col].p_allele;
            mp->nucs[col].p_allele = mp->nucs[col].q_allele;
            mp->nucs[col].q_allele = a;
            mp->nucs[col].p_freq   = 1 - mp->nucs[col].p_freq;
            uint b = p_cnt;
            p_cnt = q_cnt;
            q_cnt = b;
        }

        //
        // Check if this is a private allele. Either the site is variable and
        // the allele exists in one population, or the site is fixed and one
        // population is homozygous for the private allele.
        //
        variable_pop = -1;

        if (p_cnt == 1 && q_cnt > 1) {
            for (uint j = 0; j < this->_pop_cnt; j++)
                if (s[j]->nucs[col].p_nuc == mp->nucs[col].p_allele ||
                    s[j]->nucs[col].q_nuc == mp->nucs[col].p_allele)
                    variable_pop = j;
        } else if (p_cnt > 1 && q_cnt == 1) {
            for (uint j = 0; j < this->_pop_cnt; j++)
                if (s[j]->nucs[col].p_nuc == mp->nucs[col].q_allele ||
                    s[j]->nucs[col].q_nuc == mp->nucs[col].q_allele)
                    variable_pop = j;
        }
        mp->nucs[col].priv_allele = variable_pop;
    }

    return 0;
}

int
LocPopSum::tally_ref_alleles(int snp_index, uint16_t &allele_cnt,
                             char &p_allele, char &q_allele,
                             uint16_t &p_cnt, uint16_t &q_cnt)
{
    int  nucs[4] = {0};
    char nuc[2];

    p_allele   = 0;
    q_allele   = 0;
    allele_cnt = 0;

    for (uint j = 0; j < this->_pop_cnt; j++) {
        nuc[0] = 0;
        nuc[1] = 0;
        nuc[0] = this->_per_pop[j]->nucs[snp_index].p_nuc;
        nuc[1] = this->_per_pop[j]->nucs[snp_index].q_nuc;

        for (uint k = 0; k < 2; k++)
            switch(nuc[k]) {
            case 'A':
            case 'a':
                nucs[0]++;
                break;
            case 'C':
            case 'c':
                nucs[1]++;
                break;
            case 'G':
            case 'g':
                nucs[2]++;
                break;
            case 'T':
            case 't':
                nucs[3]++;
                break;
            }
    }

    //
    // Determine how many alleles are present at this position in this population.
    // We cannot deal with more than two alternative alleles, if there are more than two
    // in a single population, print a warning and exclude this nucleotide position.
    //
    int i;
    for (i = 0; i < 4; i++)
        if (nucs[i] > 0) allele_cnt++;

    if (allele_cnt > 2) {
        p_allele = 0;
        q_allele = 0;
        return 0;
    }

    //
    // Record which nucleotide is the P allele and which is the Q allele.
    //
    i = 0;
    while (p_allele == 0 && i < 4) {
        if (nucs[i] > 0) {
            switch(i) {
            case 0:
                p_allele = 'A';
                break;
            case 1:
                p_allele = 'C';
                break;
            case 2:
                p_allele = 'G';
                break;
            case 3:
                p_allele = 'T';
                break;
            }
        }
        i++;
    }
    while (q_allele == 0 && i < 4) {
        if (nucs[i] > 0) {
            switch(i) {
            case 1:
                q_allele = 'C';
                break;
            case 2:
                q_allele = 'G';
                break;
            case 3:
                q_allele = 'T';
                break;
            }
        }
        i++;
    }

    //
    // Tabulate the number of populations the p_allele and the q_allele occur in.
    //
    p_cnt = 0;
    q_cnt = 0;

    for (uint j = 0; j < this->_pop_cnt; j++) {
        nuc[0] = 0;
        nuc[1] = 0;
        nuc[0] = this->_per_pop[j]->nucs[snp_index].p_nuc;
        nuc[1] = this->_per_pop[j]->nucs[snp_index].q_nuc;

        for (uint k = 0; k < 2; k++)
            if (nuc[k] != 0 && nuc[k] == p_allele)
                p_cnt++;
            else if (nuc[k] != 0 && nuc[k] == q_allele)
                q_cnt++;
    }

    return 1;
}

int
LocPopSum::tally_observed_haplotypes(const vector<char *> &obshap, int snp_index)
{
    int  nucs[4] = {0};
    char nuc;

    //
    // Pull each allele for this SNP from the observed haplotype.
    //
    for (uint j = 0; j < obshap.size(); j++) {
        nuc = obshap[j][snp_index];

        switch(nuc) {
        case 'A':
        case 'a':
            nucs[0]++;
            break;
        case 'C':
        case 'c':
            nucs[1]++;
            break;
        case 'G':
        case 'g':
            nucs[2]++;
            break;
        case 'T':
        case 't':
            nucs[3]++;
            break;
        }
    }

    int allele_cnt = 0;
    for (int i = 0; i < 4; i++)
        if (nucs[i] > 0) allele_cnt++;

    return allele_cnt;
}
