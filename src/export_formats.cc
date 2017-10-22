
// -*-mode:c++; c-style:k&r; c-basic-offset:4;-*-
//
// Copyright 2012-2017, Julian Catchen <jcatchen@illinois.edu>
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
#include <algorithm>
#include <vector>

#include "ordered.h"
#include "sql_utilities.h"
#include "MetaPopInfo.h"
#include "nucleotides.h"
#include "Vcf.h"

#include "export_formats.h"

using namespace std;

extern InputMode   input_mode;
extern int         batch_id;
extern string      in_path;
extern string      out_path;
extern string      out_prefix;
extern bool        phylip_var;
extern bool        loci_ordered;
extern bool        ordered_export;
extern bool        merge_sites;
extern bool        write_gtypes;
extern bool        expand_id;
extern string      enz;
extern set<string> debug_flags;

extern MetaPopInfo      mpopi;
extern map<string, int> renz_olap;

int
Export::transpose(ifstream &ifh, vector<string> &transposed)
{
    vector<string> fields;
    string buf;
    size_t len;
    char   line[max_len];

    transposed.clear();

    do {
        buf.clear();

        //
        // Read the one line from the file.
        //
        ifh.getline(line, max_len);

        if (!ifh.good()) {
            return 0;
        }

        //
        // Check if we read a full line.
        //
        while (ifh.fail() && !ifh.eof()) {
            buf += line;
            ifh.clear();
            ifh.getline(line, max_len);
        }

        len = strlen(line);
        if (len > 0 && line[len - 1] == '\r') line[len - 1] = '\0';

        buf += line;

        //
        // Break the line up by tabs.
        //
        parse_tsv(buf.c_str(), fields);

        if (transposed.size() == 0) {
            transposed.push_back("");
            for (uint i = 0; i < fields.size(); i++)
                transposed.push_back(fields[i]);

        } else {
            assert(fields.size() == transposed.size());
            for (uint i = 0; i < fields.size(); i++)
                transposed.at(i) += "\t" + fields[i];
        }
    } while (!ifh.eof());

    return 1;
}

SumstatsExport::SumstatsExport() : Export(ExportType::sumstats)
{
    this->_path = out_path + out_prefix + ".sumstats.tsv";
}

int
SumstatsExport::open(const MetaPopInfo *mpopi)
{
    this->_mpopi   = mpopi;
    this->_pop_cnt = this->_mpopi->pops().size();

    this->_fh.open(this->_path.c_str(), ofstream::out);
    if (this->_fh.fail()) {
        cerr << "Error opening sumstats file '" << this->_path << "'\n";
        exit(1);
    }
    this->_fh.precision(fieldw);
    this->_fh.setf(std::ios::fixed);

    cerr << "Population-level summary statistics will be written to '" << this->_path << "'\n";

    return 0;
}

int
SumstatsExport::write_header()
{
    //
    // Write the population members.
    //
    for (auto& pop : mpopi.pops()) {
        this->_fh << "# " << pop.name << "\t";
        for (size_t i = pop.first_sample; i <= pop.last_sample; i++) {
            this->_fh << mpopi.samples()[i].name;
            if (i < pop.last_sample)
                this->_fh << ",";
        }
        this->_fh << "\n";
    }

    this->_fh << "# Locus ID" << "\t"
              << "Chr"      << "\t"
              << "BP"       << "\t"
              << "Col"      << "\t"
              << "Pop ID"   << "\t"
              << "P Nuc"    << "\t"
              << "Q Nuc"    << "\t"
              << "N"        << "\t"
              << "P"        << "\t"
              << "Obs Het"  << "\t"
              << "Obs Hom"  << "\t"
              << "Exp Het"  << "\t"
              << "Exp Hom"  << "\t"
              << "Pi"       << "\t"
              << "Smoothed Pi"  << "\t"
              << "Smoothed Pi P-value"  << "\t"
              << "Fis"          << "\t"
              << "Smoothed Fis" << "\t"
              << "Smoothed Fis P-value" << "\t"
              << "Private"      << "\n";
    return 0;
}

int
SumstatsExport::write_batch(const vector<LocBin *> &loci)
{
    CSLocus        *cloc;
    const LocSum   *s;
    const LocTally *t;
    double          p_freq;

    for (uint i = 0; i < loci.size(); i++) {
        cloc = loci[i]->cloc;
        t    = loci[i]->s->meta_pop();

        for (uint pos = 0; pos < cloc->len; pos++) {

            if (!t->nucs[pos].fixed) {
                for (uint pop = 0; pop < this->_pop_cnt; pop++) {

                    s = loci[i]->s->per_pop(pop);

                    if (s->nucs[pos].num_indv == 0)
                        continue;

                    this->_fh
                       << cloc->id << "\t"
                       << cloc->loc.chr() << "\t"
                       << cloc->sort_bp(pos) + 1 << "\t"
                       << pos << "\t"
                       << this->_mpopi->pops()[pop].name << "\t";

                    //
                    // Output the p and q alleles in the same order in each population.
                    //
                    if (t->nucs[pos].p_allele == s->nucs[pos].p_nuc) {
                        if (s->nucs[pos].q_nuc == 0)
                            this->_fh << s->nucs[pos].p_nuc << "\t" << "-";
                        else
                            this->_fh << s->nucs[pos].p_nuc << "\t" << s->nucs[pos].q_nuc;
                        p_freq = s->nucs[pos].p;

                    } else {
                        if (s->nucs[pos].q_nuc == 0)
                            this->_fh << "-\t" << s->nucs[pos].p_nuc;
                        else
                            this->_fh << s->nucs[pos].q_nuc << "\t" << s->nucs[pos].p_nuc;
                        p_freq = 1 - s->nucs[pos].p;
                    }

                    this->_fh << "\t" << (int) s->nucs[pos].num_indv << "\t"
                              << std::setprecision(8)      << p_freq << "\t"
                              << std::setprecision(fieldw) << s->nucs[pos].obs_het << "\t"
                              << s->nucs[pos].obs_hom   << "\t"
                              << s->nucs[pos].exp_het   << "\t"
                              << s->nucs[pos].exp_hom   << "\t"
                              << s->nucs[pos].stat[0]   << "\t" // Pi
                              << s->nucs[pos].smoothed[0] << "\t"  // Smoothed Pi
                              << s->nucs[pos].bs[0]       << "\t"  // Pi bootstrapped p-value
                              << (s->nucs[pos].stat[1] == -7.0 ? 0.0 : s->nucs[pos].stat[1]) << "\t"  // Fis
                              << s->nucs[pos].smoothed[1] << "\t"  // Smoothed Fis
                              << s->nucs[pos].bs[1]       << "\t"; // Fis bootstrapped p-value.
                    (t->nucs[pos].priv_allele == (int) pop) ? this->_fh << "1\n" : this->_fh << "0\n";
                }
            }
        }
    }

    return 0;
}

HapstatsExport::HapstatsExport() : Export(ExportType::hapstats)
{
    this->_path = out_path + out_prefix + ".hapstats.tsv";
}

int
HapstatsExport::open(const MetaPopInfo *mpopi)
{
    this->_mpopi   = mpopi;
    this->_pop_cnt = this->_mpopi->pops().size();

    this->_fh.open(this->_path.c_str(), ofstream::out);
    if (this->_fh.fail()) {
        cerr << "Error opening sumstats file '" << this->_path << "'\n";
        exit(1);
    }
    this->_fh.precision(fieldw);
    this->_fh.setf(std::ios::fixed);

    cerr << "Population-level haplotype summary statistics will be written to '" << this->_path << "'\n";

    return 0;
}

int
HapstatsExport::write_header()
{
    //
    // Write the population members.
    //
    for (auto& pop : this->_mpopi->pops()) {
        this->_fh << "# " << pop.name << "\t";

        for (size_t i = pop.first_sample; i <= pop.last_sample; i++) {
            this->_fh << this->_mpopi->samples()[i].name;
            if (i < pop.last_sample)
                this->_fh << ",";
        }
        this->_fh << "\n";
    }

    this->_fh
        << "# Locus ID"     << "\t"
        << "Chr"            << "\t"
        << "BP"             << "\t"
        << "Pop ID"         << "\t"
        << "N"              << "\t"
        << "Haplotype Cnt"  << "\t"
        << "Gene Diversity" << "\t"
        << "Smoothed Gene Diversity"      << "\t"
        << "Smoothed Gene Diversity P-value"      << "\t"
        << "Haplotype Diversity"          << "\t"
        << "Smoothed Haplotype Diversity" << "\t"
        << "Smoothed Haplotype Diversity P-value" << "\t"
        << "Haplotypes"                   << "\n";

    return 0;
}

int
HapstatsExport::write_batch(const vector<LocBin *> &loci)
{
    const LocStat *l;

    for (uint i = 0; i < loci.size(); i++) {
        const vector<Pop> &pops = this->_mpopi->pops();

        for (uint pop = 0; pop < this->_pop_cnt; pop++) {

            l = loci[i]->s->hapstats_per_pop(pop);

            if (l == NULL)
                continue;

            this->_fh
                << loci[i]->cloc->id        << "\t"
                << loci[i]->cloc->loc.chr() << "\t"
                << l->bp + 1        << "\t"
                << pops[pop].name   << "\t"
                << (int) l->alleles << "\t"
                << l->hap_cnt       << "\t"
                << l->stat[0]       << "\t"
                << l->smoothed[0]   << "\t"
                << l->bs[0]         << "\t"
                << l->stat[1]       << "\t"
                << l->smoothed[1]   << "\t"
                << l->bs[1]         << "\t"
                << l->hap_str       << "\n";
        }
    }

    return 0;
}

SnpDivergenceExport::SnpDivergenceExport() : Export(ExportType::snpdivergence)
{
    return;
}

int
SnpDivergenceExport::open(const MetaPopInfo *mpopi)
{
    this->_mpopi = mpopi;

    //
    // Open an SNP Divergence output file for each pair of populations.
    //
    string    path;
    ofstream *fh;

    for (uint pop_1 = 0; pop_1 < this->_mpopi->pops().size(); pop_1++) {
        const Pop& pop_1p = this->_mpopi->pops()[pop_1];

        for (uint pop_2 = pop_1 + 1; pop_2 < this->_mpopi->pops().size(); pop_2++) {
            const Pop& pop_2p = this->_mpopi->pops()[pop_2];

            path = out_path + out_prefix + ".fst_" + pop_1p.name + "-" + pop_2p.name + ".tsv";
            fh = new ofstream(path.c_str(), ofstream::out);
            if (fh->fail()) {
                cerr << "Error opening Fst output file '" << path << "'\n";
                exit(1);
            }
            fh->precision(fieldw);
            fh->setf(std::ios::fixed);

            this->_fhs.push_back(fh);
        }
    }

    return 0;
}

int
SnpDivergenceExport::write_header()
{
    //
    // Write a header into each SNP Divergence output file (one for each pair of populations).
    //
    ofstream *fh;

    uint i = 0;
    for (uint pop_1 = 0; pop_1 < this->_mpopi->pops().size(); pop_1++) {
        for (uint pop_2 = pop_1 + 1; pop_2 < this->_mpopi->pops().size(); pop_2++) {
            fh = this->_fhs[i];

            *fh << "# Locus ID" << "\t"
                << "Pop 1 ID"   << "\t"
                << "Pop 2 ID"   << "\t"
                << "Chr"        << "\t"
                << "BP"         << "\t"
                << "Column"     << "\t"
                << "Overall Pi" << "\t"
                << "AMOVA Fst"  << "\t"
                << "Fisher's P" << "\t"
                << "Odds Ratio" << "\t"
                << "CI Low"     << "\t"
                << "CI High"    << "\t"
                << "LOD"        << "\t"
                << "Corrected AMOVA Fst"        << "\t"
                << "Smoothed AMOVA Fst"         << "\t"
                << "Smoothed AMOVA Fst P-value" << "\t"
                << "Window SNP Count";

            //
            // If requested, log Fst component calculations to a file.
            //
            if (log_fst_comp) {
                *fh << "\t"
                    << "n_1" << "\t"
                    << "n_2" << "\t"
                    << "tot_alleles" << "\t"
                    << "p_1" << "\t"
                    << "q_1" << "\t"
                    << "p_2" << "\t"
                    << "q_2" << "\t"
                    << "pi_1" << "\t"
                    << "pi_2" << "\t"
                    << "pi_all" << "\t"
                    << "bcoeff_1" << "\t"
                    << "bcoeff_2" << "\t"
                    << "binomial_fst" << "\t"
                    << "p_1_freq" << "\t"
                    << "q_1_freq" << "\t"
                    << "p_2_freq" << "\t"
                    << "q_2_freq" << "\t"
                    << "p_avg_cor" << "\t"
                    << "n_avg_cor" << "\t"
                    << "amova_fst" << "\n";
            } else {
                *fh << "\n";
            }

            i++;
        }
    }
    return 0;
}

int
SnpDivergenceExport::write_batch_pairwise(const vector<LocBin *> &loci, const vector<vector<PopPair **>> &div)
{
    ofstream *fh;
    PopPair **pp;
    LocBin   *loc;
    size_t    cloc_len;

    for (uint i = 0; i < div.size(); i++) {

        fh = this->_fhs[i];
        assert(div[i].size() == loci.size());

        for (uint j = 0; j < div[i].size(); j++) {

            loc      = loci[j];
            cloc_len = loc->cloc->len;
            pp       = div[i][j];

            for (uint pos = 0; pos < cloc_len; pos++) {

                if (pp[pos] == NULL)
                        continue;

                *fh << pp[pos]->loc_id      << "\t"
                    << this->_mpopi->pops()[pp[pos]->pop_1].name << "\t"
                    << this->_mpopi->pops()[pp[pos]->pop_2].name << "\t"
                    << loc->cloc->loc.chr() << "\t"
                    << pp[pos]->bp + 1      << "\t"
                    << pp[pos]->col         << "\t"
                    << pp[pos]->pi          << "\t"
                    << pp[pos]->amova_fst   << "\t"
                    << std::setprecision(9)      << pp[pos]->fet_p  << "\t"
                    << std::setprecision(fieldw) << pp[pos]->fet_or << "\t"
                    << pp[pos]->ci_low      << "\t"
                    << pp[pos]->ci_high     << "\t"
                    << pp[pos]->lod         << "\t"
                    << pp[pos]->stat[1]     << "\t"
                    << pp[pos]->smoothed[1] << "\t"
                    << pp[pos]->bs[1]       << "\t"
                    << pp[pos]->snp_cnt;

                if (log_fst_comp) {
                    *fh << "\t"
                        << pp[pos]->comp[0]   << "\t"
                        << pp[pos]->comp[1]   << "\t"
                        << pp[pos]->comp[2]   << "\t"
                        << pp[pos]->comp[3]   << "\t"
                        << pp[pos]->comp[4]   << "\t"
                        << pp[pos]->comp[5]   << "\t"
                        << pp[pos]->comp[6]   << "\t"
                        << pp[pos]->comp[7]   << "\t"
                        << pp[pos]->comp[8]   << "\t"
                        << pp[pos]->comp[9]   << "\t"
                        << pp[pos]->comp[10]  << "\t"
                        << pp[pos]->comp[11]  << "\t"
                        << pp[pos]->fst       << "\t"
                        << pp[pos]->comp[12]  << "\t"
                        << pp[pos]->comp[13]  << "\t"
                        << pp[pos]->comp[14]  << "\t"
                        << pp[pos]->comp[15]  << "\t"
                        << pp[pos]->comp[16]  << "\t"
                        << pp[pos]->comp[17]  << "\t"
                        << pp[pos]->amova_fst << "\n";
                } else {
                    *fh << "\n";
                }
            }
        }
    }
    return 0;
}

HapDivergenceExport::HapDivergenceExport() : Export(ExportType::hapdivergence)
{
    return;
}

int
HapDivergenceExport::open(const MetaPopInfo *mpopi)
{
    this->_mpopi = mpopi;

    //
    // Open an Haplotype Divergence output file for each pair of populations.
    //
    string    path;
    ofstream *fh;

    for (uint pop_1 = 0; pop_1 < this->_mpopi->pops().size(); pop_1++) {
        const Pop& pop_1p = this->_mpopi->pops()[pop_1];

        for (uint pop_2 = pop_1 + 1; pop_2 < this->_mpopi->pops().size(); pop_2++) {
            const Pop& pop_2p = this->_mpopi->pops()[pop_2];

            path = out_path + out_prefix + ".phistats_" + pop_1p.name + "-" + pop_2p.name + ".tsv";
            fh = new ofstream(path.c_str(), ofstream::out);
            if (fh->fail()) {
                cerr << "Error opening Fst output file '" << path << "'\n";
                exit(1);
            }
            fh->precision(fieldw);
            fh->setf(std::ios::fixed);

            this->_fhs.push_back(fh);
        }
    }

    //
    // Open an single Haplotype Divergence output file for the full meta population.
    //
    path = out_path + out_prefix + ".phistats.tsv";
    fh = new ofstream(path.c_str(), ofstream::out);
    if (fh->fail()) {
        cerr << "Error opening Fst output file '" << path << "'\n";
        exit(1);
    }
    fh->precision(fieldw);
    fh->setf(std::ios::fixed);

    this->_metapop_fh = fh;

    return 0;
}

int
HapDivergenceExport::write_header()
{
    //
    // Write a header into each haplotype divergence output file (one for each pair of populations).
    //
    ofstream *fh;

    uint i = 0;
    for (uint pop_1 = 0; pop_1 < this->_mpopi->pops().size(); pop_1++) {
        for (uint pop_2 = pop_1 + 1; pop_2 < this->_mpopi->pops().size(); pop_2++) {
            vector<int> subpop_ids;
            subpop_ids.push_back(pop_1);
            subpop_ids.push_back(pop_2);
            fh = this->_fhs[i];

            //
            // Write the population members.
            //
            for (int k : subpop_ids) {
                const Pop& pop_k = this->_mpopi->pops()[k]; // This is [pop_i], then [pop_j].
                *fh << "# Population " << pop_k.name << "\t";
                for (size_t n = pop_k.first_sample; n <= pop_k.last_sample; n++) {
                    *fh << this->_mpopi->samples()[n].name;
                    if (n < pop_k.last_sample)
                        *fh << ",";
                }
                *fh << "\n";
            }

            *fh << "# Locus ID" << "\t"
                << "Pop 1 ID"   << "\t"
                << "Pop 2 ID"   << "\t"
                << "Chr"        << "\t"
                << "BP"         << "\t";
            if (log_fst_comp)
                *fh << "SSD(WP)"     << "\t"
                    << "SSD(AP/WG)"  << "\t"
                    << "SSD(AG)"     << "\t"
                    << "SSD(TOTAL)"  << "\t"
                    << "MSD(WP)"     << "\t"
                    << "MSD(AP/WG)"  << "\t"
                    << "MSD(AG)"     << "\t"
                    << "MSD(TOTAL)"  << "\t"
                    << "n"           << "\t"
                    << "n'"          << "\t"
                    << "n''"         << "\t"
                    << "Sigma2_a"    << "\t"
                    << "Sigma2_b"    << "\t"
                    << "Sigma2_c"    << "\t"
                    << "Sigma_Total" << "\t";
            *fh << "phi_st"          << "\t"
                << "Smoothed Phi_st" << "\t"
                << "Smoothed Phi_st P-value" << "\t"
                << "Fst'"            << "\t"
                << "Smoothed Fst'"   << "\t"
                << "Smoothed Fst' P-value"   << "\t"
                << "D_est"          << "\t"
                << "Smoothed D_est" << "\t"
                << "Smoothed D_est P-value" << "\n";

            i++;
        }
    }

    //
    // Write a header into the meta population haplotype divergence output file.
    //
    fh = this->_metapop_fh;

    //
    // Write the population members.
    //
    for (auto& pop : this->_mpopi->pops()) {
        *fh << "# Population " << pop.name << "\t";
        for (size_t k = pop.first_sample; k <= pop.last_sample; k++) {
            *fh << this->_mpopi->samples()[k].name;
            if (k < pop.last_sample)
                *fh << ",";
        }
        *fh << "\n";
    }

    //
    // Write the group members.
    //
    for (auto& group : this->_mpopi->groups()) {
        *fh << "# Group " << group.name << "\t";
        for (size_t i_pop : group.pops) {
            *fh << this->_mpopi->pops()[i_pop].name;
            if (i_pop != group.pops.back())
                *fh << ",";
        }
        *fh << "\n";
    }

    *fh << "# Locus ID" << "\t"
        << "Chr"        << "\t"
        << "BP"         << "\t";
    if (log_fst_comp)
        *fh << "SSD(WP)"     << "\t"
            << "SSD(AP/WG)"  << "\t"
            << "SSD(AG)"     << "\t"
            << "SSD(TOTAL)"  << "\t"
            << "MSD(WP)"     << "\t"
            << "MSD(AP/WG)"  << "\t"
            << "MSD(AG)"     << "\t"
            << "MSD(TOTAL)"  << "\t"
            << "n"           << "\t"
            << "n'"          << "\t"
            << "n''"         << "\t"
            << "Sigma2_a"    << "\t"
            << "Sigma2_b"    << "\t"
            << "Sigma2_c"    << "\t"
            << "Sigma_Total" << "\t";
    *fh << "phi_st"          << "\t"
        << "Smoothed Phi_st" << "\t"
        << "Smoothed Phi_st P-value" << "\t"
        << "Fst'"            << "\t"
        << "Smoothed Fst'"   << "\t"
        << "Smoothed Fst' P-value"   << "\t"
        << "D_est"          << "\t"
        << "Smoothed D_est" << "\t"
        << "Smoothed D_est P-value" << "\n";

    return 0;
}

int
HapDivergenceExport::write_batch_pairwise(const vector<LocBin *> &loci,
                                          const vector<vector<HapStat *>> &div,
                                          const vector<HapStat *> &metapop_div)
{
    ofstream *fh;
    HapStat  *h;
    LocBin   *loc;

    for (uint i = 0; i < div.size(); i++) {

        fh = this->_fhs[i];

        assert(div[i].size() == loci.size());

        for (uint j = 0; j < div[i].size(); j++) {
            loc = loci[j];
            h   = div[i][j];

            if (h == NULL) continue;

            *fh << h->loc_id    << "\t"
                << this->_mpopi->pops()[h->pop_1].name << "\t"
                << this->_mpopi->pops()[h->pop_2].name << "\t"
                << loc->cloc->loc.chr() << "\t"
                << h->bp +1     << "\t";
            if (log_fst_comp)
                *fh << h->comp[0]  << "\t"
                    << h->comp[1]  << "\t"
                    << h->comp[2]  << "\t"
                    << h->comp[3]  << "\t"
                    << h->comp[4]  << "\t"
                    << h->comp[5]  << "\t"
                    << h->comp[6]  << "\t"
                    << h->comp[7]  << "\t"
                    << h->comp[8]  << "\t"
                    << h->comp[9]  << "\t"
                    << h->comp[10] << "\t"
                    << h->comp[11] << "\t"
                    << h->comp[12] << "\t"
                    << h->comp[13] << "\t"
                    << h->comp[14] << "\t";
            *fh << h->stat[0]     << "\t"
                << h->smoothed[0] << "\t"
                << h->bs[0]       << "\t"
                << h->stat[3]     << "\t"
                << h->smoothed[3] << "\t"
                << h->bs[3]       << "\t"
                << h->stat[4]     << "\t"
                << h->smoothed[4] << "\t"
                << h->bs[4]       << "\n";
        }
    }

    fh = this->_metapop_fh;

    assert(metapop_div.size() == loci.size());

    for (uint j = 0; j < metapop_div.size(); j++) {
        loc = loci[j];
        h   = metapop_div[j];

        if (h == NULL) continue;

        *fh << h->loc_id    << "\t"
            << loc->cloc->loc.chr() << "\t"
            << h->bp +1     << "\t";
        if (log_fst_comp)
            *fh << h->comp[0]  << "\t"
                << h->comp[1]  << "\t"
                << h->comp[2]  << "\t"
                << h->comp[3]  << "\t"
                << h->comp[4]  << "\t"
                << h->comp[5]  << "\t"
                << h->comp[6]  << "\t"
                << h->comp[7]  << "\t"
                << h->comp[8]  << "\t"
                << h->comp[9]  << "\t"
                << h->comp[10] << "\t"
                << h->comp[11] << "\t"
                << h->comp[12] << "\t"
                << h->comp[13] << "\t"
                << h->comp[14] << "\t";
        *fh << h->stat[0]     << "\t"
            << h->smoothed[0] << "\t"
            << h->bs[0]       << "\t"
            << h->stat[3]     << "\t"
            << h->smoothed[3] << "\t"
            << h->bs[3]       << "\t"
            << h->stat[4]     << "\t"
            << h->smoothed[4] << "\t"
            << h->bs[4]       << "\n";
    }

    return 0;
}

GenotypesExport::GenotypesExport() : Export(ExportType::genotypes)
{
    this->_path = out_path + out_prefix + (write_gtypes ? ".genotypes.tsv" : ".haplotypes.tsv");
}

int
GenotypesExport::open(const MetaPopInfo *mpopi)
{
    this->_mpopi = mpopi;

    this->_fh.open(this->_path.c_str(), ofstream::out);
    if (this->_fh.fail()) {
        cerr << "Error opening markers file '" << this->_path << "'\n";
        exit(1);
    }

    cerr << "Raw Genotypes/Haplotypes will be written to '" << this->_path << "'\n";

    return 0;
}

int
GenotypesExport::write_header()
{
    this->_fh << "# Catalog Locus ID"    << "\t";
    if (expand_id)
        this->_fh << "\t";
    if (write_gtypes)
        this->_fh << "Marker\t";
    this->_fh << "Cnt";

    for (size_t i = 0; i < this->_mpopi->samples().size(); i++) {
        this->_fh << "\t" << this->_mpopi->samples()[i].name;
    }
    this->_fh << "\n";

    return 0;
}

int
GenotypesExport::write_batch(const vector<LocBin *> &loci)
{
    CSLocus *cloc;

    //
    // Output each locus.
    //
    for (uint i = 0; i < loci.size(); i++) {
        cloc = loci[i]->cloc;

        stringstream id;
        cloc->annotation.length() > 0 ?
            id << cloc->id << "|" << cloc->annotation : id << cloc->id;

        this->_fh << id.str();

        if (expand_id) {
            if (cloc->annotation.length() > 0)
                id << "\t" << cloc->id << "\t" << cloc->annotation;
            else if (strlen(cloc->loc.chr()) > 0)
                id << "\t" << cloc->id << "\t" << cloc->loc.chr() << "_" << cloc->loc.bp +1;
            else
                id << "\t" << cloc->id << "\t";
        }

        if (write_gtypes)
            this->_fh << "\t" << cloc->marker;

        write_gtypes ? this->_fh << "\t" << cloc->gcnt : this->_fh << "\t" << cloc->hcnt;

        Datum **d = loci[i]->d;
        string  obshap;

        for (size_t i = 0; i < this->_mpopi->samples().size(); i++) {
            this->_fh << "\t";

            if (d[i] == NULL)
                this->_fh << "-";
            else {
                if (write_gtypes) {
                    this->_fh << d[i]->gtype;
                } else {
                    obshap = "";
                    for (uint j = 0; j < d[i]->obshap.size(); j++)
                        obshap += string(d[i]->obshap[j]) + "/";
                    obshap = obshap.substr(0, obshap.length()-1);
                    this->_fh << obshap;
                }
            }
        }

        this->_fh << "\n";
    }

    return 0;
}

MarkersExport::MarkersExport() : Export(ExportType::markers)
{
    this->_path = out_path + out_prefix + ".markers.tsv";
}

int
MarkersExport::open(const MetaPopInfo *mpopi)
{
    this->_mpopi = mpopi;

    this->_fh.open(this->_path.c_str(), ofstream::out);
    if (this->_fh.fail()) {
        cerr << "Error opening markers file '" << this->_path << "'\n";
        exit(1);
    }
    this->_fh.precision(fieldw);
    this->_fh.setf(std::ios::fixed);

    cerr << "Genotyping markers will be written to '" << this->_path << "'\n";

    return 0;
}

int
MarkersExport::write_header()
{
    this->_fh << "# Catalog Locus ID"    << "\t"
              << "\t"
              << "Total Genotypes"     << "\t"
              << "Max"                 << "\t"
              << "Genotype Freqs"      << "\t"
              << "F"                   << "\t"
              << "Mean Log Likelihood" << "\t"
              << "Genotype Map"        << "\t"
              << "\n";
    return 0;
}

int
MarkersExport::write_batch(const vector<LocBin *> &loci)
{
    stringstream gtype_map;

    for (uint i = 0; i < loci.size(); i++) {
        string freq  = "";
        double max   = 0.0;
        int    total = 0;
        gtype_map.str("");

        if (loci[i]->cloc->marker.length() > 0) {
            tally_haplotype_freq(loci[i]->cloc, loci[i]->d, this->_mpopi->samples().size(), total, max, freq);

            //
            // Record the haplotype to genotype map.
            //
            map<string, string>::iterator j;
            for (j = loci[i]->cloc->gmap.begin(); j != loci[i]->cloc->gmap.end(); j++)
                gtype_map << j->first << ":" << j->second << ";";
        }

        this->_fh
            << loci[i]->cloc->id  << "\t"
            << "\t"                       // Marker
            << total              << "\t"
            << max                << "\t"
            << freq               << "\t"
            << loci[i]->cloc->f   << "\t"
            << 1.0                << "\t" //(Formerly `cloc->lnl`.)
            << gtype_map.str()    << "\t"
            << "\n";
    }

    return 0;
}

FastaLociExport::FastaLociExport() : Export(ExportType::fasta_loci)
{
    this->_path = out_path + out_prefix + ".loci.fa";
}

int
FastaLociExport::open(const MetaPopInfo *mpopi)
{
    this->_mpopi = mpopi;

    this->_fh.open(this->_path.c_str(), ofstream::out);
    if (this->_fh.fail()) {
        cerr << "Error opening FASTA Loci file '" << this->_path << "'\n";
        exit(1);
    }

    cerr << "FASTA consensus sequences for each locus in the metapopulation  will be written to '" << this->_path << "'\n";

    return 0;
}

int
FastaLociExport::write_header()
{
    //
    // Obtain the current date.
    //
    time_t     rawtime;
    struct tm *timeinfo;
    char       date[32];
    time(&rawtime);
    timeinfo = localtime(&rawtime);
    strftime(date, 32, "%B %d, %Y", timeinfo);

    this->_fh << "Stacks version " << VERSION << "; " << date << "\n";
    return 0;
}

int
FastaLociExport::write_batch(const vector<LocBin *> &loci)
{
    for (uint i = 0; i < loci.size(); i++) {
        LocBin* loc = loci[i];

        this->_fh << ">CLocus_" << loc->cloc->id;
        if (strcmp(loc->cloc->loc.chr(), "un") != 0)
            this->_fh << " [" << loc->cloc->loc.chr() << ", " << loc->cloc->sort_bp() + 1 << ", " << (loc->cloc->loc.strand == strand_plus ? "+" : "-");
        this->_fh << '\n' << loc->cloc->con << '\n';

#ifdef DEBUG
        bool no_samples = true;
        Datum** d = loc->d;
        for (uint j = 0; j < this->_mpopi->samples().size(); j++)
            if (d[j] != NULL)
                no_samples = false;
        assert(!no_samples);
#endif
    }

    return 0;
}

FastaRawExport::FastaRawExport() : Export(ExportType::fasta_raw)
{
    this->_path = out_path + out_prefix + ".samples-raw.fa";
}

int
FastaRawExport::write_header()
{
    //
    // Obtain the current date.
    //
    time_t     rawtime;
    struct tm *timeinfo;
    char       date[32];
    time(&rawtime);
    timeinfo = localtime(&rawtime);
    strftime(date, 32, "%B %d, %Y", timeinfo);

    this->_fh << "Stacks version " << VERSION << "; " << date << "\n";
    return 0;
}

int
FastaRawExport::open(const MetaPopInfo *mpopi)
{
    this->_mpopi = mpopi;

    this->_fh.open(this->_path.c_str(), ofstream::out);
    if (this->_fh.fail()) {
        cerr << "Error opening FASTA raw samples file '" << this->_path << "'\n";
        exit(1);
    }

    cerr << "Raw FASTA consensus sequences for each sample will be written to '" << this->_path << "'\n";

    return 0;
}

int
FastaRawExport::write_batch(const vector<LocBin *> &loci)
{
    //
    // Write a FASTA file containing each allele from each locus from
    // each sample in the population.
    //
    LocBin *loc;
    Datum **d;
    char   *seq;

    for (uint i = 0; i < loci.size(); i++) {
        loc = loci[i];
        d   = loc->d;
        seq = new char[loc->cloc->len + 1];
        strcpy(seq, loc->cloc->con);

        for (uint j = 0; j < this->_mpopi->samples().size(); j++) {
            if (d[j] == NULL)
                continue;

            for (uint k = 0; k < d[j]->obshap.size(); k++) {

                for (uint i = 0; i < loc->cloc->snps.size(); i++) {
                    uint col = loc->cloc->snps[i]->col;
                    seq[col] = col < loc->cloc->len ? d[j]->obshap[k][i] : loc->cloc->con[col];
                }

                this->_fh << ">CLocus_" << loc->cloc->id
                          << "_Sample_" << this->_mpopi->samples()[j].id
                          << "_Locus_"  << d[j]->id
                          << "_Allele_" << k
                          << " ["       << this->_mpopi->samples()[j].name;

                if (strcmp(loc->cloc->loc.chr(), "un") != 0)
                    this->_fh << "; " << loc->cloc->loc.chr() << ", " << loc->cloc->sort_bp() + 1 << ", " << (loc->cloc->loc.strand == strand_plus ? "+" : "-");
                this->_fh << "]\n"
                          << seq << "\n";
            }
        }
        delete [] seq;
    }

    return 0;
}

FastaSamplesExport::FastaSamplesExport() : Export(ExportType::fasta_samples)
{
    this->_path = out_path + out_prefix + ".samples.fa";
}

int
FastaSamplesExport::write_header()
{
    //
    // Obtain the current date.
    //
    time_t     rawtime;
    struct tm *timeinfo;
    char       date[32];
    time(&rawtime);
    timeinfo = localtime(&rawtime);
    strftime(date, 32, "%B %d, %Y", timeinfo);

    this->_fh << "Stacks version " << VERSION << "; " << date << "\n";
    return 0;
}

int
FastaSamplesExport::open(const MetaPopInfo *mpopi)
{
    this->_mpopi = mpopi;

    this->_fh.open(this->_path.c_str(), ofstream::out);
    if (this->_fh.fail()) {
        cerr << "Error opening FASTA samples file '" << this->_path << "'\n";
        exit(1);
    }

    cerr << "FASTA consensus sequences for each sample will be written to '" << this->_path << "'\n";

    return 0;
}

int
FastaSamplesExport::write_batch(const vector<LocBin *> &loci)
{
    LocBin *loc;
    Datum **d;
    char   *seq;

    for (uint i = 0; i < loci.size(); i++) {
        loc = loci[i];
        d   = loc->d;
        seq = new char[loc->cloc->len + 1];
        strcpy(seq, loc->cloc->con);

        for (uint j = 0; j < this->_mpopi->samples().size(); j++) {
            if (d[j] == NULL)
                continue;
            if (d[j]->obshap.size() > 2)
                continue;

            if (d[j]->obshap.size() == 1) {

                for (uint i = 0; i < loc->cloc->snps.size(); i++) {
                    uint col = loc->cloc->snps[i]->col;
                    seq[col] = col < loc->cloc->len ? d[j]->obshap[0][i] : loc->cloc->con[col];
                }

                this->_fh << ">CLocus_" << loc->cloc->id
                          << "_Sample_" << this->_mpopi->samples()[j].id
                          << "_Locus_"  << d[j]->id
                          << "_Allele_" << 0
                          << " ["       << this->_mpopi->samples()[j].name;
                if (strcmp(loc->cloc->loc.chr(), "un") != 0)
                    this->_fh << "; " << loc->cloc->loc.chr() << ", " << loc->cloc->sort_bp() + 1 << ", " << (loc->cloc->loc.strand == strand_plus ? "+" : "-");
                this->_fh << "]\n"
                          << seq << "\n";

                this->_fh << ">CLocus_" << loc->cloc->id
                          << "_Sample_" << this->_mpopi->samples()[j].id
                          << "_Locus_"  << d[j]->id
                          << "_Allele_" << 1
                          << " ["       << mpopi.samples()[j].name;
                if (strcmp(loc->cloc->loc.chr(), "un") != 0)
                    this->_fh << "; " << loc->cloc->loc.chr() << ", " << loc->cloc->sort_bp() + 1 << ", " << (loc->cloc->loc.strand == strand_plus ? "+" : "-");
                this->_fh << "]\n"
                          << seq << "\n";

            } else {
                for (uint k = 0; k < d[j]->obshap.size(); k++) {
                    for (uint i = 0; i < loc->cloc->snps.size(); i++) {
                        uint col = loc->cloc->snps[i]->col;
                        seq[col] = col < loc->cloc->len ? d[j]->obshap[k][i] : loc->cloc->con[col];
                    }

                    this->_fh << ">CLocus_" << loc->cloc->id
                              << "_Sample_" <<  this->_mpopi->samples()[j].id
                              << "_Locus_"  << d[j]->id
                              << "_Allele_" << k
                              << " ["       <<  this->_mpopi->samples()[j].name;
                    if (strcmp(loc->cloc->loc.chr(), "un") != 0)
                        this->_fh << "; " << loc->cloc->loc.chr() << ", " << loc->cloc->sort_bp() + 1 << ", " << (loc->cloc->loc.strand == strand_plus ? "+" : "-");
                    this->_fh << "]\n"
                              << seq << "\n";
                }
            }
        }

        delete [] seq;
    }

    return 0;
}
/*
int
write_vcf_ordered(map<int, CSLocus *> &catalog,
                  PopMap<CSLocus> *pmap, PopSum<CSLocus> *psum,
                  map<int, pair<merget, int> > &merge_map, ofstream &log_fh)
{
    //
    // Write a VCF file as defined here: http://www.1000genomes.org/node/101
    // Order SNPs by genomic position (and handle overlapping loci).
    //

    string path = out_path + out_prefix + ".vcf";
    cerr << "Writing ordered population data to VCF file '" << path << "'\n";
    log_fh << "\n#\n# Generating SNP-based VCF export.\n#\n";

    VcfHeader header;
    header.add_std_meta();
    for(auto& s : mpopi.samples())
        header.add_sample(s.name);
    VcfWriter writer (path, move(header));

    // We need to order the SNPs taking into account overlapping loci.
    OLocTally<NucTally> *ord = new OLocTally<NucTally>(log_fh);

    for (auto it = pmap->ordered_loci().begin(); it != pmap->ordered_loci().end(); it++) {
        vector<NucTally *> sites;
        // ord->order(sites, it->second);

        for (uint pos = 0; pos < sites.size(); pos++) {
            if (catalog.count(sites[pos]->loc_id) == 0) {
                cerr << "Unable to find locus id " << sites[pos]->loc_id << "\n";
                continue;
            }
            CSLocus* loc = catalog[sites[pos]->loc_id];
            uint16_t col = sites[pos]->col;
            int snp_index = loc->snp_index(col);
            if (snp_index < 0) {
                cerr << "Warning: unable to locate SNP call in column " << col << " for catalog locus " << loc->id << "\n";
                continue;
            }
            Datum** d = pmap->locus(loc->id);


            const char ref = sites[pos]->p_allele;
            const char alt = sites[pos]->q_allele;
            char freq_alt[32];
            sprintf(freq_alt, "%0.3f", 1 - sites[pos]->p_freq);

            VcfRecord rec;
            rec.append_chrom(string(loc->loc.chr()));
            rec.append_pos(loc->sort_bp(col) + 1);
            rec.append_id(to_string(loc->id) + "_" + to_string(col));
            rec.append_allele(Nt2(loc->loc.strand == strand_plus ? ref : reverse(ref)));
            rec.append_allele(Nt2(loc->loc.strand == strand_plus ? alt : reverse(alt)));
            rec.append_qual(".");
            rec.append_filters("PASS");
            rec.append_info(string("NS=") + to_string(sites[pos]->num_indv));
            rec.append_info(string("AF=") + freq_alt);
            rec.append_format("GT");
            rec.append_format("DP");
            rec.append_format("AD");

            for (int j = 0; j < pmap->sample_cnt(); j++) {
                stringstream sample;

                if (d[j] == NULL || col >= d[j]->len) {
                    // Data does not exist.
                    sample << "./.:0:.,.";
                } else if (d[j]->model[col] == 'U') {
                    // Data exists, but the model call was uncertain.
                    sample << "./.:" << d[j]->tot_depth << ":.,.";
                } else {
                    char allele1, allele2;
                    tally_observed_haplotypes(d[j]->obshap, snp_index, allele1, allele2);

                    if (allele1 == 0) {
                        // More than two alleles in this sample
                        sample << "./.:" << d[j]->tot_depth << ":.,.";
                    } else {
                        // Write the genotype.

                        int dp1, dp2;
                        find_datum_allele_depths(d[j], snp_index, allele1, allele2, dp1, dp2);

                        if (allele2 == 0) {
                            // homozygote
                            if (allele1 == ref) {
                                sample << "0/0:" << d[j]->tot_depth << ":" << dp1 << "," << dp2;
                            } else {
                                sample << "1/1:" << d[j]->tot_depth << ":" << dp2 << "," << dp1;
                            }
                        } else {
                            // heterozygote
                            sample << "0/1:" << d[j]->tot_depth
                                   << ":" << (allele1 == ref ? dp1 : dp2) << "," << (allele1 == ref ? dp2 : dp1);
                        }
                    }
                }
                rec.append_sample(sample.str());
            }
            writer.write_record(rec);
        }
    }

    return 0;
}

int
write_vcf(map<int, CSLocus *> &catalog,
          PopMap<CSLocus> *pmap, PopSum<CSLocus> *psum,
          map<int, pair<merget, int> > &merge_map)
{
    //
    // Write a VCF file as defined here: http://www.1000genomes.org/node/101
    //

    string path = out_path + out_prefix + ".vcf";
    cerr << "Writing population data to VCF file '" << path << "'\n";

    VcfHeader header;
    header.add_std_meta();
    header.add_meta(VcfMeta::predefs::info_locori);
    for(auto& s : mpopi.samples())
        header.add_sample(s.name);
    VcfWriter writer (path, move(header));

    map<string, vector<CSLocus *> >::const_iterator it;
    CSLocus *loc;
    Datum   **d;
    LocTally *t;

    for (it = pmap->ordered_loci().begin(); it != pmap->ordered_loci().end(); it++) {

        // We need to order the SNPs so negative and positive strand SNPs are properly ordered.
        vector<GenPos> ordered_snps;
        for (uint pos = 0; pos < it->second.size(); pos++) {
            loc = it->second[pos];
            t   = psum->locus_tally(loc->id);

            for (uint i = 0; i < loc->snps.size(); i++) {
                uint16_t col = loc->snps[i]->col;
                if (t->nucs[col].allele_cnt == 2)
                    ordered_snps.push_back(GenPos(loc->id, i, loc->sort_bp(col)));
            }
        }
        sort(ordered_snps.begin(), ordered_snps.end());

        for (uint pos = 0; pos < ordered_snps.size(); pos++) {
            loc = catalog[ordered_snps[pos].id];
            uint16_t col = loc->snps[ordered_snps[pos].snp_index]->col;
            int snp_index = loc->snp_index(col);
            if (snp_index < 0) {
                cerr << "Warning: unable to locate SNP call in column " << col << " for locus #" << loc->id << "\n";
                continue;
            }
            t   = psum->locus_tally(loc->id);
            d = pmap->locus(loc->id);

            const char ref = t->nucs[col].p_allele;
            const char alt = t->nucs[col].q_allele;
            char freq_alt[32];
            sprintf(freq_alt, "%0.3f", 1 - t->nucs[col].p_freq);

            VcfRecord rec;
            rec.append_chrom(string(loc->loc.chr()));
            rec.append_pos(loc->sort_bp(col) + 1);
            rec.append_id(to_string(loc->id) + "_" + to_string(col));
            rec.append_allele(Nt2(loc->loc.strand == strand_plus ? ref : reverse(ref)));
            rec.append_allele(Nt2(loc->loc.strand == strand_plus ? alt : reverse(alt)));
            rec.append_qual(".");
            rec.append_filters("PASS");
            rec.append_info(string("NS=") + to_string(t->nucs[col].num_indv));
            rec.append_info(string("AF=") + freq_alt);
            rec.append_info(string("locori=") + (loc->loc.strand == strand_plus ? "p" : "m"));
            rec.append_format("GT");
            rec.append_format("DP");
            rec.append_format("AD");

            for (int j = 0; j < pmap->sample_cnt(); j++) {
                stringstream sample;

                if (d[j] == NULL || col >= uint(d[j]->len)) {
                    // Data does not exist.
                    sample << "./.:0:.,.";
                } else if (d[j]->model[col] == 'U') {
                    // Data exists, but the model call was uncertain.
                    sample << "./.:" << d[j]->tot_depth << ":.,.";
                } else {
                    char allele1, allele2;
                    tally_observed_haplotypes(d[j]->obshap, ordered_snps[pos].snp_index, allele1, allele2);

                    if (allele1 == 0) {
                        // More than two alleles in this sample
                        sample << "./.:" << d[j]->tot_depth << ":.,.";
                    } else {
                        // Write the genotype.

                        int dp1, dp2;
                        find_datum_allele_depths(d[j], snp_index, allele1, allele2, dp1, dp2);

                        if (allele2 == 0) {
                            // homozygote
                            if (allele1 == ref) {
                                sample << "0/0:" << d[j]->tot_depth << ":" << dp1 << "," << dp2;
                            } else {
                                sample << "1/1:" << d[j]->tot_depth << ":" << dp2 << "," << dp1;
                            }
                        } else {
                            // heterozygote
                            sample << "0/1:" << d[j]->tot_depth
                                   << ":" << (allele1 == ref ? dp1 : dp2) << "," << (allele1 == ref ? dp2 : dp1);
                        }
                    }
                }
                rec.append_sample(sample.str());
            }
            writer.write_record(rec);
        }
    }

    return 0;
}

int
write_vcf_haplotypes(map<int, CSLocus *> &catalog,
                     PopMap<CSLocus> *pmap,
                     PopSum<CSLocus> *psum)
{
    //
    // Write a VCF file as defined here: http://samtools.github.io/hts-specs/
    //XXX Datum::obshap *is not* ordered, I think. See Bitbucket. @Nick (June 2016)
    //

    string path = out_path + out_prefix + ".haplotypes.vcf";
    cerr << "Writing population data haplotypes to VCF file '" << path << "'\n";

    VcfHeader header;
    header.add_std_meta();
    for(auto& s : mpopi.samples())
        header.add_sample(s.name);
    VcfWriter writer (path, move(header));

    CSLocus  *loc;
    Datum   **d;
    char      allele[id_len];

    for (auto chr = pmap->ordered_loci().begin(); chr != pmap->ordered_loci().end(); chr++) {
        for (uint pos = 0; pos < chr->second.size(); pos++) {
            loc = chr->second[pos];
            d   = pmap->locus(loc->id);

            map<string, double> hap_freq;
            const double n_alleles = count_haplotypes_at_locus(0, pmap->sample_cnt() - 1, (const Datum **) d, hap_freq);
            if (hap_freq.size() <= 1)
                // Monomorphic locus.
                // XXX What does [hap_freq.size()==1] mean ? @Nick (July 2016)
                continue;
            for(auto& h : hap_freq)
                // Convert counts to freqs.
                h.second /= n_alleles;

            //
            // Order the haplotypes according to most frequent. Record the ordered position or each
            // haplotype and convert them from counts to frequencies.
            //

            VcfRecord rec;
            rec.append_chrom(string(loc->loc.chr()));
            rec.append_pos(loc->sort_bp() + 1);
            rec.append_id(to_string(loc->id));

            //alleles
            vector<pair<string, double> > ordered_hap (hap_freq.begin(), hap_freq.end());
            sort(ordered_hap.begin(), ordered_hap.end(), compare_pair_haplotype);
            map<string, int> hap_index;
            for (size_t i = 0; i < ordered_hap.size(); i++) {
                string h = ordered_hap[i].first;
                rec.append_allele(loc->loc.strand == strand_plus ? h : rev_comp(h));
                hap_index[h] = i;
            }

            rec.append_qual(".");
            rec.append_filters("PASS");

            //info
            rec.append_info(string("locori=") + (loc->loc.strand == strand_plus ? "p" : "m"));
            stringstream ss;
            ss << n_alleles/2;
            rec.append_info(string("NS=") + ss.str());
            string af = "AF=";
            sprintf(allele, "%0.3f", ordered_hap[1].second); //NB. hap_freq.size() >= 2
            af += allele;
            for (auto h=ordered_hap.begin()+2; h!=ordered_hap.end(); ++h) {
                sprintf(allele, "%0.3f", h->second);
                af += string(",") + allele;
            }
            rec.append_info(af);

            //format
            rec.append_format("GT");
            rec.append_format("DP");

            for (int j = 0; j < pmap->sample_cnt(); j++) {
                stringstream sample;

                if (d[j] == NULL) {
                    // Data does not exist.
                    sample << "./.:0";
                } else if (d[j]->obshap.size() > 2) {
                    // More than two alleles in this sample
                    sample << "./.:" << d[j]->tot_depth;
                } else if (d[j]->obshap.size() == 1) {
                    // Homozygote.
                    char* h = d[j]->obshap[0];
                    int i = uncalled_haplotype(h) ? -1 : hap_index.at(h);
                    if(i >= 0)
                        sample << i << "/" << i << ":" << d[j]->tot_depth;
                    else
                        sample << "./.:" << d[j]->tot_depth;
                } else {
                    // Heterozygote.
                    char* h1 = d[j]->obshap[0];
                    char* h2 = d[j]->obshap[1];
                    int i1 = uncalled_haplotype(h1) ? -1 : hap_index.at(h1);
                    int i2 = uncalled_haplotype(h2) ? -1 : hap_index.at(h2);
                    if(i1 >= 0 && i2 >= 0)
                        sample << (i1 < i2 ? i1 : i2) << "/" << (i1 < i2 ? i2 : i1) << ":" << d[j]->tot_depth;
                    else if (i1 >= 0)
                        sample << i1 << "/.:" << d[j]->tot_depth;
                    else if (i2 >= 0)
                        sample << i2 << "/.:" << d[j]->tot_depth;
                }
                rec.append_sample(sample.str());
            }
            writer.write_record(rec);
        }
    }
    return 0;
}
*/

GenePopExport::GenePopExport() : OrderableExport(ExportType::genepop) {}

int
GenePopExport::open(const MetaPopInfo *mpopi)
{
    this->_mpopi = mpopi;

    //
    // Write a GenePop file as defined here: http://kimura.univ-montp2.fr/~rousset/Genepop.htm
    //
    this->_path = out_path + out_prefix + ".genepop";

    //
    // Open a temporary file.
    //
    this->_tmp_path = out_path + out_prefix + ".genepop.part";
    this->_tmpfh.open(this->_tmp_path);

    cerr << "Polymorphic sites in GenePop format will be written to '" << this->_path << "'\n";

    return 0;
}

int
GenePopExport::write_header()
{
    for (size_t p = 0; p < this->_mpopi->pops().size(); ++p) {
        const Pop& pop = this->_mpopi->pops()[p];
        for (size_t j = pop.first_sample; j <= pop.last_sample; j++) {
            this->_tmpfh << mpopi.samples()[j].name;
            if (j < this->_mpopi->samples().size() - 1)
                this->_tmpfh << "\t";
        }
    }
    this->_tmpfh << "\n";

    return 0;
}

int
OrderableExport::write_batch(const vector<LocBin*> &loci)
{

    if (ordered_export) {
        //
        // We need to order the SNPs to take into account overlapping loci.
        //
        vector<const NucTally *> sites;
        OLocTally<NucTally> ord;
        ord.order(sites, loci);

        uint k = 0;
        for (uint pos = 0; pos < sites.size(); pos++) {
            int loc_id = sites[pos]->loc_id;

            //
            // Advance through the loci list until we find the proper locus.
            //
            while (k < loci.size() && loci[k]->cloc->id != loc_id)
                k++;
            assert(loci[k]->cloc->id == loc_id);

            const LocBin* loc  = loci[k];
            size_t col       = sites[pos]->col;
            size_t snp_index = loc->cloc->snp_index(col);

            this->write_site(loc->cloc, loc->s, loc->d, col, snp_index);
        }

    } else {
        for (uint k = 0; k < loci.size(); k++) {
            const LocBin* loc = loci[k];
            const CSLocus* cloc = loc->cloc;
            Datum const*const* d = loc->d;
            const LocTally* t = loc->s->meta_pop();

            for (uint snp_index = 0; snp_index < cloc->snps.size(); snp_index++) {
                uint col = cloc->snps[snp_index]->col;

                if (t->nucs[col].allele_cnt != 2)
                    continue;

                this->write_site(cloc, loc->s, d, col, snp_index);
            }
        }
    }

    return 0;
}

int
GenePopExport::write_site(const CSLocus *loc, const LocPopSum *lps, Datum const*const* d, size_t col, size_t snp_index)
{
    map<char, string> nuc_map;
    nuc_map['A'] = "01";
    nuc_map['C'] = "02";
    nuc_map['G'] = "03";
    nuc_map['T'] = "04";

    const LocSum *s;
    char p_allele, q_allele;

    this->_tmpfh << loc->id << "_" << col;

    for (size_t p = 0; p < this->_mpopi->pops().size(); ++p) {
        const Pop& pop = this->_mpopi->pops()[p];
        s = lps->per_pop(p);

        for (size_t j = pop.first_sample; j <= pop.last_sample; j++) {

            if (s->nucs[col].incompatible_site ||
                s->nucs[col].filtered_site) {
                //
                // This site contains more than two alleles in this population or was filtered
                // due to a minor allele frequency that is too low.
                //
                this->_tmpfh << "\t0000";
            } else if (d[j] == NULL || col >= uint(d[j]->len)) {
                //
                // Data does not exist.
                //
                this->_tmpfh << "\t0000";
            } else if (d[j]->model[col] == 'U') {
                //
                // Data exists, but the model call was uncertain.
                //
                this->_tmpfh << "\t0000";
            } else {
                //
                // Tally up the nucleotide calls.
                //
                tally_observed_haplotypes(d[j]->obshap, snp_index, p_allele, q_allele);

                if (p_allele == 0 && q_allele == 0) {
                    // More than two potential alleles.
                    this->_tmpfh << "\t0000";
                } else if (p_allele == 0) {
                    this->_tmpfh << "\t" << nuc_map[q_allele] << nuc_map[q_allele];

                } else if (q_allele == 0) {
                    this->_tmpfh << "\t" << nuc_map[p_allele] << nuc_map[p_allele];

                } else {
                    this->_tmpfh << "\t" << nuc_map[p_allele] << nuc_map[q_allele];
                }
            }
        }
    }
    this->_tmpfh << "\n";
    return 0;
}

int
GenePopExport::post_processing()
{
    //
    // Obtain the current date.
    //
    time_t     rawtime;
    struct tm *timeinfo;
    char       date[32];
    time(&rawtime);
    timeinfo = localtime(&rawtime);
    strftime(date, 32, "%B %d, %Y", timeinfo);

    this->_fh.open(this->_path.c_str(), ofstream::out);
    if (this->_fh.fail()) {
        cerr << "Error opening GenePop file '" << this->_path << "'\n";
        return(0);
    }

    //
    // Output the header line.
    //
    this->_fh << "# Stacks version " << VERSION << "; Genepop version 4.1.3; " << date << "\n";

    this->_intmpfh.open(this->_tmp_path.c_str(), ofstream::in);
    if (this->_intmpfh.fail()) {
        cerr << "Error opening GenePop file '" << this->_tmp_path << "'\n";
        return(0);
    }

    vector<string> transposed_lines;

    this->transpose(this->_intmpfh, transposed_lines);

    assert(transposed_lines.size() == this->_mpopi->samples().size() + 1);

    size_t line_cnt = 0;
    size_t pos;
    //
    // The first line has a list of locus IDs, convert these to comma-separated.
    //
    for (uint i = 1; i < transposed_lines[line_cnt].size(); i++)
        if (transposed_lines[line_cnt][i] == '\t')
            transposed_lines[line_cnt][i] = ',';

    this->_fh << transposed_lines[line_cnt].substr(1) << "\n";
    line_cnt++;

    for (size_t p = 0; p < this->_mpopi->pops().size(); ++p) {
        const Pop& pop = this->_mpopi->pops()[p];
        this->_fh << "pop\n";

        for (size_t j = pop.first_sample; j <= pop.last_sample; j++) {
            pos = transposed_lines[line_cnt].find_first_of('\t');
            transposed_lines[line_cnt][pos] = ',';
            this->_fh
                << transposed_lines[line_cnt].substr(0, pos + 1) << "\t"
                << transposed_lines[line_cnt].substr(pos + 1) << "\n";
            line_cnt++;
        }
    }

    return 1;
}

void
GenePopExport::close()
{
    //
    // Close and delete the temporary files.
    //
    this->_intmpfh.close();

    this->_tmpfh.close();
    remove(this->_tmp_path.c_str());

    this->_fh.close();
    return;
}
int
VcfExport::open(const MetaPopInfo *mpopi)
{
    this->_path = out_path + out_prefix + ".vcf";
    cerr << "SNPs and calls will be written in VCF format to '" << this->_path << "'\n";

    this->_mpopi = mpopi;

    VcfHeader header;
    header.add_std_meta();
    header.add_meta(VcfMeta::predefs::info_locori);
    for(auto& s : this->_mpopi->samples())
        header.add_sample(s.name);

    this->_writer = new VcfWriter(this->_path, move(header));

    return 0;
}

int VcfExport::write_site(
        const CSLocus* cloc,
        const LocPopSum* psum,
        Datum const*const* d,
        size_t col,
        size_t index
) {
    const LocTally* t = psum->meta_pop();

    const char ref = t->nucs[col].p_allele;
    const char alt = t->nucs[col].q_allele;
    char freq_alt[32];
    sprintf(freq_alt, "%0.3f", 1 - t->nucs[col].p_freq);

    VcfRecord rec;
    rec.append_chrom(string(cloc->loc.chr()));
    rec.append_pos(cloc->sort_bp(col) + 1);
    rec.append_id(to_string(cloc->id) + "_" + to_string(col));
    rec.append_allele(Nt2(cloc->loc.strand == strand_plus ? ref : reverse(ref)));
    rec.append_allele(Nt2(cloc->loc.strand == strand_plus ? alt : reverse(alt)));
    rec.append_qual(".");
    rec.append_filters("PASS");
    rec.append_info(string("NS=") + to_string(t->nucs[col].num_indv));
    rec.append_info(string("AF=") + freq_alt);
    rec.append_info(string("locori=") + (cloc->loc.strand == strand_plus ? "p" : "m"));
    rec.append_format("GT");
    rec.append_format("DP");
    rec.append_format("AD");

    for (size_t s=0; s<this->_mpopi->samples().size(); s++) {
        stringstream sample;

        if (d[s] == NULL || col >= uint(d[s]->len) || d[s]->model[col] == 'U') {
            // Data does not exist.
            sample << ".";
        } else {
            if (d[s]->model[col] == 'O') {
                assert(d[s]->obshap.size() == 1
                       || (d[s]->obshap.size() == 2 && d[s]->obshap[0][index] == d[s]->obshap[1][index]));
                sample << (d[s]->obshap[0][index] == ref ? "0/0" : "1/1");
            } else {
                assert(d[s]->model[col] == 'E');
                assert(d[s]->obshap.size() == 2 && d[s]->obshap[0][index] != d[s]->obshap[1][index]);
                sample << "0/1";
            }

            // DP.
            sample << ":" << d[s]->snpdata[index].tot_depth;
            // AD.
            if (d[s]->snpdata[index].nt_depths.sum() > 0)
                sample << ":" << d[s]->snpdata[index].nt_depths[Nt2(ref)]
                       << "," << d[s]->snpdata[index].nt_depths[Nt2(alt)];
            else
                sample << ":.";
        }
        rec.append_sample(sample.str());
    }
    this->_writer->write_record(rec);

    return 0;
}

/*
int
write_structure(map<int, CSLocus *> &catalog,
                PopMap<CSLocus> *pmap,
                PopSum<CSLocus> *psum)
{
    //
    // Write a Structure file as defined here: http://pritch.bsd.uchicago.edu/structure.html
    //
    // To avoid linked SNPs (which Structure can't handle), we will only output the first
    // SNP from each variable locus.
    //
    string file = out_path + out_prefix + ".structure.tsv";

    cerr << "Writing population data to Structure file '" << file << "'...";

    ofstream fh(file.c_str(), ofstream::out);

    if (fh.fail()) {
        cerr << "Error opening Structure file '" << file << "'\n";
        exit(1);
    }

    //
    // Obtain the current date.
    //
    time_t     rawtime;
    struct tm *timeinfo;
    char       date[32];
    time(&rawtime);
    timeinfo = localtime(&rawtime);
    strftime(date, 32, "%B %d, %Y", timeinfo);

    //
    // Output the header.
    //
    fh << "# Stacks v" << VERSION << "; " << " Structure v2.3; " << date << "\n";

    map<string, vector<CSLocus *> >::const_iterator it;
    CSLocus  *loc;
    Datum   **d;
    LocSum  **s;
    LocTally *t;

    for (it = pmap->ordered_loci().begin(); it != pmap->ordered_loci().end(); it++) {
        for (uint pos = 0; pos < it->second.size(); pos++) {
            loc = it->second[pos];
            t   = psum->locus_tally(loc->id);

            for (uint i = 0; i < loc->snps.size(); i++) {
                uint col = loc->snps[i]->col;
                if (t->nucs[col].allele_cnt == 2)
                    fh << "\t" << loc->id << "_" << col;
            }
        }
    }
    fh << "\n";

    map<char, string> nuc_map;
    nuc_map['A'] = "1";
    nuc_map['C'] = "2";
    nuc_map['G'] = "3";
    nuc_map['T'] = "4";

    char      p_allele, q_allele;

    for (size_t p=0; p<mpopi.pops().size(); ++p) {
        const Pop& pop = mpopi.pops()[p];

        for (size_t j = pop.first_sample; j <= pop.last_sample; j++) {
            //
            // Output all the loci for this sample, printing only the p allele
            //
            fh << mpopi.samples()[j].name << "\t" << pop.name;

            for (it = pmap->ordered_loci().begin(); it != pmap->ordered_loci().end(); it++) {
                for (uint pos = 0; pos < it->second.size(); pos++) {
                    loc = it->second[pos];

                    s = psum->locus(loc->id);
                    d = pmap->locus(loc->id);
                    t = psum->locus_tally(loc->id);

                    for (uint i = 0; i < loc->snps.size(); i++) {
                        uint col = loc->snps[i]->col;

                        //
                        // If this site is fixed in all populations or has too many alleles don't output it.
                        //
                        if (t->nucs[col].allele_cnt != 2)
                            continue;

                        if (s[p]->nucs[col].incompatible_site ||
                            s[p]->nucs[col].filtered_site) {
                            //
                            // This site contains more than two alleles in this population or was filtered
                            // due to a minor allele frequency that is too low.
                            //
                            fh << "\t" << "0";
                        } else if (d[j] == NULL || col >= uint(d[j]->len)) {
                            //
                            // Data does not exist.
                            //
                            fh << "\t" << "0";
                        } else if (d[j]->model[col] == 'U') {
                            //
                            // Data exists, but the model call was uncertain.
                            //
                            fh << "\t" << "0";
                        } else {
                            //
                            // Tally up the nucleotide calls.
                            //
                            tally_observed_haplotypes(d[j]->obshap, i, p_allele, q_allele);

                            if (p_allele == 0 && q_allele == 0)
                                fh << "\t" << "0";
                            else if (p_allele == 0)
                                fh << "\t" << nuc_map[q_allele];
                            else
                                fh << "\t" << nuc_map[p_allele];
                        }
                    }
                }
            }
            fh << "\n";

            //
            // Output all the loci for this sample again, now for the q allele
            //
            fh << mpopi.samples()[j].name << "\t" << pop.name;

            for (it = pmap->ordered_loci().begin(); it != pmap->ordered_loci().end(); it++) {
                for (uint pos = 0; pos < it->second.size(); pos++) {
                    loc = it->second[pos];

                    s = psum->locus(loc->id);
                    d = pmap->locus(loc->id);
                    t = psum->locus_tally(loc->id);

                    for (uint i = 0; i < loc->snps.size(); i++) {
                        uint col = loc->snps[i]->col;

                        if (t->nucs[col].allele_cnt != 2)
                            continue;

                        if (s[p]->nucs[col].incompatible_site ||
                            s[p]->nucs[col].filtered_site) {
                            fh << "\t" << "0";
                        } else if (d[j] == NULL || col >= uint(d[j]->len)) {
                            fh << "\t" << "0";
                        } else if (d[j]->model[col] == 'U') {
                            fh << "\t" << "0";
                        } else {
                            tally_observed_haplotypes(d[j]->obshap, i, p_allele, q_allele);

                            if (p_allele == 0 && q_allele == 0)
                                fh << "\t" << "0";
                            else if (q_allele == 0)
                                fh << "\t" << nuc_map[p_allele];
                            else
                                fh << "\t" << nuc_map[q_allele];
                        }
                    }
                }
            }
            fh << "\n";
        }
    }

    fh.close();

    cerr << "done.\n";

    return 0;
}

int
write_structure_ordered(map<int, CSLocus *> &catalog,
                        PopMap<CSLocus> *pmap,
                        PopSum<CSLocus> *psum,
                        ofstream &log_fh)
{
    //
    // Write a Structure file as defined here: http://pritch.bsd.uchicago.edu/structure.html
    //
    // To avoid linked SNPs (which Structure can't handle), we will only output the first
    // SNP from each variable locus.
    //
    string file = out_path + out_prefix + ".structure.tsv";

    cerr << "Writing population data to Structure file '" << file << "'...";

    ofstream fh(file.c_str(), ofstream::out);

    if (fh.fail()) {
        cerr << "Error opening Structure file '" << file << "'\n";
        exit(1);
    }

    //
    // Obtain the current date.
    //
    time_t     rawtime;
    struct tm *timeinfo;
    char       date[32];
    time(&rawtime);
    timeinfo = localtime(&rawtime);
    strftime(date, 32, "%B %d, %Y", timeinfo);

    //
    // Output the header.
    //
    fh << "# Stacks v" << VERSION << "; " << " Structure v2.3; " << date << "\n";

    map<string, vector<NucTally *> > genome_sites;
    map<string, vector<CSLocus *> >::const_iterator it;
    CSLocus  *loc;
    Datum   **d;
    LocSum  **s;

    //
    // We need to order the SNPs to take into account overlapping loci.
    //
    OLocTally<NucTally> *ord = new OLocTally<NucTally>(log_fh);

    for (it = pmap->ordered_loci().begin(); it != pmap->ordered_loci().end(); it++) {
        vector<NucTally *> &sites = genome_sites[it->first];
        //ord->order(sites, it->second);

        for (uint pos = 0; pos < sites.size(); pos++)
            fh << "\t" << sites[pos]->loc_id << "_" << sites[pos]->col;
    }
    fh << "\n";

    map<char, string> nuc_map;
    nuc_map['A'] = "1";
    nuc_map['C'] = "2";
    nuc_map['G'] = "3";
    nuc_map['T'] = "4";

    char      p_allele, q_allele;
    uint      col, snp_index;

    for (size_t p=0; p<mpopi.pops().size(); ++p) {
        const Pop& pop = mpopi.pops()[p];

        for (size_t j = pop.first_sample; j <= pop.last_sample; j++) {
            //
            // Output all the loci for this sample, printing only the p allele
            //
            fh << mpopi.samples()[j].name << "\t" << pop.name;

            for (it = pmap->ordered_loci().begin(); it != pmap->ordered_loci().end(); it++) {
                vector<NucTally *> &sites = genome_sites[it->first];

                for (uint pos = 0; pos < sites.size(); pos++) {
                    loc = catalog[sites[pos]->loc_id];
                    s   = psum->locus(loc->id);
                    d   = pmap->locus(loc->id);
                    col = sites[pos]->col;

                    if (s[p]->nucs[col].incompatible_site ||
                        s[p]->nucs[col].filtered_site) {
                        //
                        // This site contains more than two alleles in this population or was filtered
                        // due to a minor allele frequency that is too low.
                        //
                        fh << "\t" << "0";
                    } else if (d[j] == NULL || col >= uint(d[j]->len)) {
                        //
                        // Data does not exist.
                        //
                        fh << "\t" << "0";
                    } else if (d[j]->model[col] == 'U') {
                        //
                        // Data exists, but the model call was uncertain.
                        //
                        fh << "\t" << "0";
                    } else {
                        snp_index = loc->snp_index(col);
                        //
                        // Tally up the nucleotide calls.
                        //
                        tally_observed_haplotypes(d[j]->obshap, snp_index, p_allele, q_allele);

                        if (p_allele == 0 && q_allele == 0)
                            fh << "\t" << "0";
                        else if (p_allele == 0)
                            fh << "\t" << nuc_map[q_allele];
                        else
                            fh << "\t" << nuc_map[p_allele];
                    }
                }
            }
            fh << "\n";

            //
            // Output all the loci for this sample again, now for the q allele
            //
            fh << mpopi.samples()[j].name << "\t" << pop.name;

            for (it = pmap->ordered_loci().begin(); it != pmap->ordered_loci().end(); it++) {
                vector<NucTally *> &sites = genome_sites[it->first];

                for (uint pos = 0; pos < sites.size(); pos++) {
                    loc = catalog[sites[pos]->loc_id];
                    s   = psum->locus(loc->id);
                    d   = pmap->locus(loc->id);
                    col = sites[pos]->col;

                    if (s[p]->nucs[col].incompatible_site ||
                        s[p]->nucs[col].filtered_site) {
                        fh << "\t" << "0";
                    } else if (d[j] == NULL || col >= uint(d[j]->len)) {
                        fh << "\t" << "0";
                    } else if (d[j]->model[col] == 'U') {
                        fh << "\t" << "0";
                    } else {
                        snp_index = loc->snp_index(col);
                        tally_observed_haplotypes(d[j]->obshap, snp_index, p_allele, q_allele);

                        if (p_allele == 0 && q_allele == 0)
                            fh << "\t" << "0";
                        else if (q_allele == 0)
                            fh << "\t" << nuc_map[p_allele];
                        else
                            fh << "\t" << nuc_map[q_allele];
                    }
                }
            }
            fh << "\n";
        }
    }

    fh.close();

    cerr << "done.\n";

    return 0;
}

int
write_hzar(map<int, CSLocus *> &catalog,
           PopMap<CSLocus> *pmap,
           PopSum<CSLocus> *psum)
{
    //
    // Write a Hybrid Zone Analysis using R (HZAR) file as defined here:
    //    http://cran.r-project.org/web/packages/hzar/hzar.pdf
    //
    string file = out_path + out_prefix + ".hzar.csv";

    cerr << "Writing population data to HZAR file '" << file << "'...";

    ofstream fh(file.c_str(), ofstream::out);

    if (fh.fail()) {
        cerr << "Error opening HZAR file '" << file << "'\n";
        exit(1);
    }

    //
    // Obtain the current date.
    //
    time_t     rawtime;
    struct tm *timeinfo;
    char       date[32];
    time(&rawtime);
    timeinfo = localtime(&rawtime);
    strftime(date, 32, "%B %d, %Y", timeinfo);

    //
    // Output the header.
    //
    fh << "# Stacks v" << VERSION << "; " << " HZAR v0.2-5; " << date << "\n"
       << "Population" << "," << "Distance";

    map<string, vector<CSLocus *> >::const_iterator it;
    CSLocus  *loc;
    LocSum  **s;
    LocTally *t;

    for (it = pmap->ordered_loci().begin(); it != pmap->ordered_loci().end(); it++) {
        for (uint pos = 0; pos < it->second.size(); pos++) {
            loc = it->second[pos];
            s   = psum->locus(loc->id);
            t   = psum->locus_tally(loc->id);

            for (uint i = 0; i < loc->snps.size(); i++) {
                uint col = loc->snps[i]->col;
                if (t->nucs[col].allele_cnt == 2) {
                    fh << "," << loc->id << "_" << col << ".A"
                       << "," << loc->id << "_" << col << ".B"
                       << "," << loc->id << "_" << col << ".N";
                }
            }
        }
    }
    fh << "\n";

    for (size_t p=0; p<mpopi.pops().size(); ++p) {
        const Pop& pop = mpopi.pops()[p];

        fh << pop.name << ",";

        for (it = pmap->ordered_loci().begin(); it != pmap->ordered_loci().end(); it++) {
            for (uint pos = 0; pos < it->second.size(); pos++) {
                loc = it->second[pos];

                s = psum->locus(loc->id);
                t = psum->locus_tally(loc->id);

                for (uint i = 0; i < loc->snps.size(); i++) {
                    uint col = loc->snps[i]->col;

                    //
                    // If this site is fixed in all populations or has too many alleles don't output it.
                    //
                    if (t->nucs[col].allele_cnt != 2)
                        continue;

                    if (s[p]->nucs[col].num_indv == 0 ||
                        s[p]->nucs[col].incompatible_site ||
                        s[p]->nucs[col].filtered_site) {
                        fh << ",0,0,0";
                        continue;
                    }

                    if (t->nucs[col].p_allele == s[p]->nucs[col].p_nuc)
                        fh << "," << s[p]->nucs[col].p << "," << 1 - s[p]->nucs[col].p << ",";
                    else
                        fh << "," << 1 - s[p]->nucs[col].p << "," << s[p]->nucs[col].p << ",";

                    fh << s[p]->nucs[col].num_indv * 2;
                }
            }
        }
        fh << "\n";
    }

    fh.close();

    cerr << "done.\n";

    return 0;
}

int
write_treemix(map<int, CSLocus *> &catalog,
              PopMap<CSLocus> *pmap,
              PopSum<CSLocus> *psum)
{
    //
    // Write a TreeMix file (Pickrell and Pritchard, 2012 PLoS Genetics)
    //    https://bitbucket.org/nygcresearch/treemix/wiki/Home
    //
    string file = out_path + out_prefix + ".treemix";

    cerr << "Writing population data to TreeMix file '" << file << "'; ";

    ofstream fh(file.c_str(), ofstream::out);

    if (fh.fail()) {
        cerr << "Error opening TreeMix file '" << file << "'\n";
        exit(1);
    }

    file += ".log";

    cerr << "logging nucleotide positions to '" << file << "'...";

    ofstream log_fh(file.c_str(), ofstream::out);

    if (log_fh.fail()) {
        cerr << "Error opening Phylip Log file '" << file << "'\n";
        exit(1);
    }

    //
    // Obtain the current date.
    //
    time_t     rawtime;
    struct tm *timeinfo;
    char       date[32];
    time(&rawtime);
    timeinfo = localtime(&rawtime);
    strftime(date, 32, "%B %d, %Y", timeinfo);

    log_fh << "# Stacks v" << VERSION << "; " << " TreeMix v1.1; " << date << "\n"
           << "# Line\tLocus ID\tColumn\tChr\tBasepair\n";

    //
    // Output the header.
    //
    fh << "# Stacks v" << VERSION << "; " << " TreeMix v1.1; " << date << "\n";

    map<string, vector<CSLocus *> >::const_iterator it;
    CSLocus  *loc;
    LocSum  **s;
    LocTally *t;

    //
    // Output a space-separated list of the populations on the first line.
    //
    stringstream sstr;
    for (auto& pop : mpopi.pops())
        sstr << pop.name << " ";

    fh << sstr.str().substr(0, sstr.str().length() - 1) << "\n";

    double p_freq, p_cnt, q_cnt, allele_cnt;
    long int line = 1;

    for (it = pmap->ordered_loci().begin(); it != pmap->ordered_loci().end(); it++) {
        for (uint pos = 0; pos < it->second.size(); pos++) {
            loc = it->second[pos];

            s = psum->locus(loc->id);
            t = psum->locus_tally(loc->id);

            for (uint i = 0; i < loc->snps.size(); i++) {
                uint col = loc->snps[i]->col;

                sstr.str("");

                //
                // If this site is fixed in all populations or has too many alleles don't output it.
                //
                if (t->nucs[col].allele_cnt != 2)
                    continue;

                for (size_t p=0; p<mpopi.pops().size(); ++p) {

                    if (s[p]->nucs[col].num_indv == 0 ||
                        s[p]->nucs[col].incompatible_site ||
                        s[p]->nucs[col].filtered_site) {
                        sstr << "0,0 ";
                        continue;
                    }

                    p_freq = (t->nucs[col].p_allele == s[p]->nucs[col].p_nuc) ?
                        s[p]->nucs[col].p :
                        1 - s[p]->nucs[col].p;

                    allele_cnt = s[p]->nucs[col].num_indv * 2;
                    p_cnt      = round(allele_cnt * p_freq);
                    q_cnt      = allele_cnt - p_cnt;
                    sstr << (int) p_cnt << "," << (int) q_cnt << " ";
                }

                if (sstr.str().length() == 0)
                    continue;

                fh << sstr.str().substr(0, sstr.str().length() - 1) << "\n";
                log_fh << line << "\t" << loc->id << "\t" << col << "\t" << loc->loc.chr() << "\t" << loc->sort_bp(col) + 1 << "\n";
                line++;
            }
        }
    }

    fh.close();
    log_fh.close();

    cerr << "done.\n";

    return 0;
}

int
write_fastphase(map<int, CSLocus *> &catalog,
                PopMap<CSLocus> *pmap,
                PopSum<CSLocus> *psum)
{
    //
    // Write a fastPHASE file as defined here: http://stephenslab.uchicago.edu/software.html
    //
    // Data will be written as independent, bi-allelic SNPs. We will write one file per chromosome.
    //
    cerr << "Writing population data to fastPHASE files...";

    map<string, vector<CSLocus *> >::const_iterator it;
    CSLocus  *loc;
    Datum   **d;
    LocSum  **s;
    LocTally *t;

    for (it = pmap->ordered_loci().begin(); it != pmap->ordered_loci().end(); it++) {

        string file = out_path + out_prefix + "." + it->first + ".fastphase.inp";

        ofstream fh(file.c_str(), ofstream::out);

        if (fh.fail()) {
            cerr << "Error opening fastPHASE file '" << file << "'\n";
            exit(1);
        }

        //
        // Tally up the number of sites
        //
        int  total_sites = 0;
        uint col;
        for (uint pos = 0; pos < it->second.size(); pos++) {
            loc = it->second[pos];
            t   = psum->locus_tally(loc->id);

            for (uint i = 0; i < loc->snps.size(); i++) {
                col = loc->snps[i]->col;
                if (t->nucs[col].allele_cnt == 2)
                    total_sites++;
            }
        }

        //
        // Output the total number of SNP sites and the number of individuals.
        //
        fh << mpopi.samples().size() << "\n"
           << total_sites    << "\n";

        //
        // We need to determine an ordering that can take into account overlapping RAD sites.
        //
        vector<GenPos> ordered_snps;
        for (uint pos = 0; pos < it->second.size(); pos++) {
            loc = it->second[pos];
            t   = psum->locus_tally(loc->id);

            for (uint i = 0; i < loc->snps.size(); i++) {
                col = loc->snps[i]->col;
                if (t->nucs[col].allele_cnt == 2)
                    ordered_snps.push_back(GenPos(loc->id, i, loc->sort_bp(col)));
            }
        }
        sort(ordered_snps.begin(), ordered_snps.end());

        //
        // Output the position of each site according to its basepair.
        //
        fh << "P";
        for (uint pos = 0; pos < ordered_snps.size(); pos++) {
            loc = catalog[ordered_snps[pos].id];
            col = loc->snps[ordered_snps[pos].snp_index]->col;
            fh << " " << ordered_snps[pos].bp +1;
        }
        fh << "\n";

        //
        // Output a line of 'S' characters, one per site, indicating that these are SNP markers.
        //
        string snp_markers, gtypes_str;
        snp_markers.assign(total_sites, 'S');
        fh << snp_markers << '\n';

        //
        // Now output each sample name followed by a new line, then all of the genotypes for that sample
        // on two lines.
        //

        char         p_allele, q_allele;
        stringstream gtypes;

        for (size_t p=0; p<mpopi.pops().size(); ++p) {
            const Pop& pop = mpopi.pops()[p];

            for (size_t j = pop.first_sample; j <= pop.last_sample; j++) {
                //
                // Output all the loci for this sample, printing only the p allele
                //
                fh << mpopi.samples()[j].name << "\n";

                gtypes.str("");
                for (uint pos = 0; pos < ordered_snps.size(); pos++) {
                    loc = catalog[ordered_snps[pos].id];
                    col = loc->snps[ordered_snps[pos].snp_index]->col;

                    s = psum->locus(loc->id);
                    d = pmap->locus(loc->id);
                    t = psum->locus_tally(loc->id);

                    //
                    // If this site is fixed in all populations or has too many alleles don't output it.
                    //
                    if (t->nucs[col].allele_cnt != 2)
                        continue;

                    if (s[p]->nucs[col].incompatible_site ||
                        s[p]->nucs[col].filtered_site) {
                        //
                        // This site contains more than two alleles in this population or was filtered
                        // due to a minor allele frequency that is too low.
                        //
                        gtypes << "? ";

                    } else if (d[j] == NULL || col >= uint(d[j]->len)) {
                        //
                        // Data does not exist.
                        //
                        gtypes << "? ";
                    } else if (d[j]->model[col] == 'U') {
                        //
                        // Data exists, but the model call was uncertain.
                        //
                        gtypes << "? ";
                    } else {
                        //
                        // Tally up the nucleotide calls.
                        //
                        tally_observed_haplotypes(d[j]->obshap, ordered_snps[pos].snp_index, p_allele, q_allele);

                        if (p_allele == 0 && q_allele == 0)
                            gtypes << "? ";
                        else if (p_allele == 0)
                            gtypes << q_allele << " ";
                        else
                            gtypes << p_allele << " ";
                    }
                }
                gtypes_str = gtypes.str();
                fh << gtypes_str.substr(0, gtypes_str.length() - 1) << "\n";

                //
                // Output all the loci for this sample again, now for the q allele
                //
                gtypes.str("");
                for (uint pos = 0; pos < ordered_snps.size(); pos++) {
                    loc = catalog[ordered_snps[pos].id];
                    col = loc->snps[ordered_snps[pos].snp_index]->col;


                    s = psum->locus(loc->id);
                    d = pmap->locus(loc->id);
                    t = psum->locus_tally(loc->id);

                    if (t->nucs[col].allele_cnt != 2)
                        continue;

                    if (s[p]->nucs[col].incompatible_site ||
                        s[p]->nucs[col].filtered_site) {
                        gtypes << "? ";

                    } else if (d[j] == NULL || col >= uint(d[j]->len)) {
                        gtypes << "? ";

                    } else if (d[j]->model[col] == 'U') {
                        gtypes << "? ";

                    } else {
                        tally_observed_haplotypes(d[j]->obshap, ordered_snps[pos].snp_index, p_allele, q_allele);

                        if (p_allele == 0 && q_allele == 0)
                            gtypes << "? ";
                        else if (q_allele == 0)
                            gtypes << p_allele << " ";
                        else
                            gtypes << q_allele << " ";
                    }
                }
                gtypes_str = gtypes.str();
                fh << gtypes_str.substr(0, gtypes_str.length() - 1) << "\n";
            }
        }

        fh.close();
    }

    cerr << "done.\n";

    return 0;
}

int
write_phase(map<int, CSLocus *> &catalog,
            PopMap<CSLocus> *pmap,
            PopSum<CSLocus> *psum)
{
    //
    // Write a PHASE file as defined here: http://stephenslab.uchicago.edu/software.html
    //
    // Data will be written as mixture of multiple allele, linked RAD sites
    // (SNPs within a single RAD locus are already phased), and bi-allelic SNPs. We
    // will write one file per chromosome.
    //
    cerr << "Writing population data to PHASE files...";

    map<string, vector<CSLocus *> >::const_iterator it;
    CSLocus  *loc;
    Datum   **d;
    LocSum  **s;
    LocTally *t;

    for (it = pmap->ordered_loci().begin(); it != pmap->ordered_loci().end(); it++) {

        string file = out_path + out_prefix + "." + it->first + ".phase.inp";

        ofstream fh(file.c_str(), ofstream::out);

        if (fh.fail()) {
            cerr << "Error opening PHASE file '" << file << "'\n";
            exit(1);
        }

        //
        // We need to determine an ordering for all legitimate loci/SNPs.
        //
        uint           col;
        vector<GenPos> ordered_loci;
        for (uint pos = 0; pos < it->second.size(); pos++) {
            loc = it->second[pos];
            t   = psum->locus_tally(loc->id);

            if (loc->snps.size() == 0) continue;

            //
            // Will we output this locus as a haplotype or as a SNP?
            //
            if (loc->snps.size() > 1) {
                //
                // Check that there aren't too many haplotypes (PHASE has a max of 50).
                //
                if (loc->alleles.size() > 40) continue;

                //
                // Iterate over the population to determine that this subset of the population
                // has data at this locus.
                //
                d = pmap->locus(loc->id);
                for (int j = 0; j < pmap->sample_cnt(); j++) {
                    if (d[j] != NULL &&
                        d[j]->obshap.size() > 0 &&
                        d[j]->obshap.size() <= 2) {
                        //
                        // Data exists, and there are the correct number of haplotypes.
                        //
                        ordered_loci.push_back(GenPos(loc->id, 0, loc->sort_bp(), haplotype));
                        break;
                    }
                }
            } else {
                col = loc->snps[0]->col;

                if (t->nucs[col].allele_cnt == 2)
                    ordered_loci.push_back(GenPos(loc->id, 0, loc->sort_bp(col), snp));
            }
        }
        sort(ordered_loci.begin(), ordered_loci.end());

        //
        // Output the total number of SNP sites and the number of individuals.
        //
        fh << mpopi.samples().size()      << "\n"
           << ordered_loci.size() << "\n";

        //
        // Output the position of each site according to its basepair.
        //
        fh << "P";
        for (uint pos = 0; pos < ordered_loci.size(); pos++)
            fh << " " << ordered_loci[pos].bp +1;
        fh << "\n";

        //
        // Output a line of 'S' characters for SNP markers, 'M' characters for multiallelic haplotypes.
        //
        for (uint pos = 0; pos < ordered_loci.size(); pos++) {
            if (pos > 0) fh << " ";
            fh << (ordered_loci[pos].type == snp ? "S" : "M");
        }
        fh << "\n";

        //
        // Now output each sample name followed by a new line, then all of the genotypes for that sample
        // on two lines.
        //

        string       gtypes_str;
        bool         found;
        char         p_allele, q_allele;
        stringstream gtypes;

        for (size_t p=0; p<mpopi.pops().size(); ++p) {
            const Pop& pop = mpopi.pops()[p];

            for (size_t j = pop.first_sample; j <= pop.last_sample; j++) {
                //
                // Output all the loci for this sample, printing only the p allele
                //
                fh << mpopi.samples()[j].name << "\n";

                gtypes.str("");
                for (uint pos = 0; pos < ordered_loci.size(); pos++) {
                    loc = catalog[ordered_loci[pos].id];

                    s = psum->locus(loc->id);
                    d = pmap->locus(loc->id);
                    t = psum->locus_tally(loc->id);

                    //
                    // Will we output this locus as a haplotype or as a SNP?
                    //
                    if (ordered_loci[pos].type == haplotype) {
                        if (d[j] == NULL) {
                            //
                            // Data does not exist.
                            //
                            gtypes << "-1 ";
                        } else {
                            //
                            // Data exists, output the first haplotype. We will assume the haplotypes are
                            // numbered by their position in the loc->strings vector.
                            //
                            if (d[j]->obshap.size() > 2)  {
                                // cerr << "Warning: too many haplotypes, catalog locus: " << loc->id << "\n";
                                gtypes << "-1 ";
                            } else {
                                found = false;
                                for (uint k = 0; k < loc->strings.size(); k++)
                                    if (d[j]->obshap[0] == loc->strings[k].first) {
                                        found = true;
                                        gtypes << k + 1 << " ";
                                    }
                                if (found == false)
                                    cerr << "Unable to find haplotype " << d[j]->obshap[0] << " from individual "
                                         << mpopi.samples()[j].name << "; catalog locus: " << loc->id << "\n";
                            }
                        }
                    } else {
                        col = loc->snps[ordered_loci[pos].snp_index]->col;

                        if (s[p]->nucs[col].incompatible_site ||
                            s[p]->nucs[col].filtered_site) {
                            //
                            // This site contains more than two alleles in this population or was filtered
                            // due to a minor allele frequency that is too low.
                            //
                            gtypes << "? ";

                        } else if (d[j] == NULL || col >= uint(d[j]->len)) {
                            //
                            // Data does not exist.
                            //
                            gtypes << "? ";
                        } else if (d[j]->model[col] == 'U') {
                            //
                            // Data exists, but the model call was uncertain.
                            //
                            gtypes << "? ";
                        } else {
                            //
                            // Tally up the nucleotide calls.
                            //
                            tally_observed_haplotypes(d[j]->obshap, ordered_loci[pos].snp_index, p_allele, q_allele);

                            if (p_allele == 0 && q_allele == 0)
                                gtypes << "? ";
                            else if (p_allele == 0)
                                gtypes << q_allele << " ";
                            else
                                gtypes << p_allele << " ";
                        }
                    }
                }
                gtypes_str = gtypes.str();
                fh << gtypes_str.substr(0, gtypes_str.length() - 1) << "\n";

                //
                // Output all the loci for this sample again, now for the q allele
                //
                gtypes.str("");
                for (uint pos = 0; pos < ordered_loci.size(); pos++) {
                    loc = catalog[ordered_loci[pos].id];

                    s = psum->locus(loc->id);
                    d = pmap->locus(loc->id);
                    t = psum->locus_tally(loc->id);

                    //
                    // Will we output this locus as a haplotype or as a SNP?
                    //
                    if (ordered_loci[pos].type == haplotype) {
                        if (d[j] == NULL) {
                            //
                            // Data does not exist.
                            //
                            gtypes << "-1 ";
                        } else {
                            //
                            // Data exists, output the second haplotype. We will assume the haplotypes are
                            // numbered by their position in the loc->strings vector.
                            //
                            if (d[j]->obshap.size() > 2)  {
                                // cerr << "Warning: too many haplotypes, catalog locus: " << loc->id << "\n";
                                gtypes << "-1 ";
                            } else if (d[j]->obshap.size() > 1)  {
                                found = false;
                                for (uint k = 0; k < loc->strings.size(); k++)
                                    if (d[j]->obshap[1] == loc->strings[k].first) {
                                        found = true;
                                        gtypes << k + 1 << " ";
                                    }
                                if (found == false)
                                    cerr << "Unable to find haplotype " << d[j]->obshap[1] << " from individual "
                                         << mpopi.samples()[j].name << "; catalog locus: " << loc->id << "\n";
                            } else {
                                found = false;
                                for (uint k = 0; k < loc->strings.size(); k++)
                                    if (d[j]->obshap[0] == loc->strings[k].first) {
                                        found = true;
                                        gtypes << k + 1 << " ";
                                    }
                                if (found == false)
                                    cerr << "Unable to find haplotype " << d[j]->obshap[0] << " from individual "
                                         << mpopi.samples()[j].name << "; catalog locus: " << loc->id << "\n";
                            }
                        }
                    } else {
                        col = loc->snps[ordered_loci[pos].snp_index]->col;

                        if (s[p]->nucs[col].incompatible_site ||
                            s[p]->nucs[col].filtered_site) {
                            gtypes << "? ";

                        } else if (d[j] == NULL || col >= uint(d[j]->len)) {
                            gtypes << "? ";

                        } else if (d[j]->model[col] == 'U') {
                            gtypes << "? ";

                        } else {
                            tally_observed_haplotypes(d[j]->obshap, ordered_loci[pos].snp_index, p_allele, q_allele);

                            if (p_allele == 0 && q_allele == 0)
                                gtypes << "? ";
                            else if (q_allele == 0)
                                gtypes << p_allele << " ";
                            else
                                gtypes << q_allele << " ";
                        }
                    }
                }
                gtypes_str = gtypes.str();
                fh << gtypes_str.substr(0, gtypes_str.length() - 1) << "\n";
            }
        }

        fh.close();
    }

    cerr << "done.\n";

    return 0;
}

int
write_plink(map<int, CSLocus *> &catalog,
            PopMap<CSLocus> *pmap,
            PopSum<CSLocus> *psum)
{
    //
    // Write a PLINK file as defined here: http://pngu.mgh.harvard.edu/~purcell/plink/data.shtml
    //
    // We will write one file per chromosome.
    //
    cerr << "Writing population data to PLINK files...";

    //
    // Obtain the current date.
    //
    time_t     rawtime;
    struct tm *timeinfo;
    char       date[32];
    time(&rawtime);
    timeinfo = localtime(&rawtime);
    strftime(date, 32, "%B %d, %Y", timeinfo);

    map<string, vector<CSLocus *> >::const_iterator it;
    CSLocus  *loc;
    Datum   **d;
    LocSum  **s;
    LocTally *t;
    string    chr;

    //
    // First, write a markers file containing each marker, the chromosome it falls on,
    // an empty centiMorgan field, and finally its genomic position in basepairs.
    //
    string file = out_path + out_prefix + ".plink.map";

    ofstream fh(file.c_str(), ofstream::out);

    if (fh.fail()) {
        cerr << "Error opening PLINK markers file '" << file << "'\n";
        exit(1);
    }

    //
    // Output the header.
    //
    fh << "# Stacks v" << VERSION << "; " << " PLINK v1.07; " << date << "\n";

    for (it = pmap->ordered_loci().begin(); it != pmap->ordered_loci().end(); it++) {
        chr = it->first;

        for (uint pos = 0; pos < it->second.size(); pos++) {
            loc = it->second[pos];
            t   = psum->locus_tally(loc->id);

            for (uint i = 0; i < loc->snps.size(); i++) {
                uint col = loc->snps[i]->col;
                if (t->nucs[col].allele_cnt == 2)
                    fh << chr << "\t"
                       << loc->id << "_" << col << "\t"
                       << "0\t"
                       << loc->sort_bp(col) +1 << "\n";
            }
        }
    }
    fh.close();

    //
    // Now output the genotypes in a separate file.
    //
    file = out_path + out_prefix + ".plink.ped";

    fh.open(file.c_str(), ofstream::out);

    if (fh.fail()) {
        cerr << "Error opening PLINK markers file '" << file << "'\n";
        exit(1);
    }

    fh << "# Stacks v" << VERSION << "; " << " PLINK v1.07; " << date << "\n";

    char p_allele, q_allele;

    //
    //  marker, output the genotypes for each sample in two successive columns.
    //
    for (size_t p=0; p<mpopi.pops().size(); ++p) {
        const Pop& pop = mpopi.pops()[p];

        for (size_t j = pop.first_sample; j <= pop.last_sample; j++) {

            fh << pop.name << "\t"
               << mpopi.samples()[j].name << "\t"
               << "0\t"  // Paternal ID
               << "0\t"  // Maternal ID
               << "0\t"  // Sex
               << "0";   // Phenotype

            for (it = pmap->ordered_loci().begin(); it != pmap->ordered_loci().end(); it++) {
                for (uint pos = 0; pos < it->second.size(); pos++) {
                    loc = it->second[pos];

                    s = psum->locus(loc->id);
                    d = pmap->locus(loc->id);
                    t = psum->locus_tally(loc->id);

                    for (uint i = 0; i < loc->snps.size(); i++) {
                        uint col = loc->snps[i]->col;

                        //
                        // If this site is fixed in all populations or has too many alleles don't output it.
                        //
                        if (t->nucs[col].allele_cnt != 2)
                            continue;
                        //
                        // Output the p and q alleles
                        //
                        if (s[p]->nucs[col].incompatible_site ||
                            s[p]->nucs[col].filtered_site) {
                            //
                            // This site contains more than two alleles in this population or was filtered
                            // due to a minor allele frequency that is too low.
                            //
                            fh << "\t" << "0" << "\t" << "0";
                        } else if (d[j] == NULL || col >= uint(d[j]->len)) {
                            //
                            // Data does not exist.
                            //
                            fh << "\t" << "0" << "\t" << "0";
                        } else if (d[j]->model[col] == 'U') {
                            //
                            // Data exists, but the model call was uncertain.
                            //
                            fh << "\t" << "0" << "\t" << "0";
                        } else {
                            //
                            // Tally up the nucleotide calls.
                            //
                            tally_observed_haplotypes(d[j]->obshap, i, p_allele, q_allele);

                            if (p_allele == 0 && q_allele == 0)
                                fh << "\t" << "0" << "\t" << "0";
                            else if (p_allele == 0)
                                fh << "\t" << q_allele << "\t" << q_allele;
                            else if (q_allele == 0)
                                fh << "\t" << p_allele << "\t" << p_allele;
                            else
                                fh << "\t" << p_allele << "\t" << q_allele;
                        }
                    }
                }
            }
            fh << "\n";
        }
    }

    fh.close();

    cerr << "done.\n";

    return 0;
}

int
write_beagle(map<int, CSLocus *> &catalog,
            PopMap<CSLocus> *pmap,
            PopSum<CSLocus> *psum)
{
    //
    // Write a Beagle file as defined here: http://faculty.washington.edu/browning/beagle/beagle.html
    //
    // We will write one file per chromosome, per population.
    //
    cerr << "Writing population data to unphased Beagle files...";

    //
    // Obtain the current date.
    //
    time_t     rawtime;
    struct tm *timeinfo;
    char       date[32];
    time(&rawtime);
    timeinfo = localtime(&rawtime);
    strftime(date, 32, "%B %d, %Y", timeinfo);

    map<string, vector<CSLocus *> >::const_iterator it;
    CSLocus  *loc;
    Datum   **d;
    LocSum  **s;
    LocTally *t;
    uint      col;

    stringstream pop_name;
    string       file;

    for (it = pmap->ordered_loci().begin(); it != pmap->ordered_loci().end(); it++) {

        //
        // We need to determine an ordering that can take into account overlapping RAD sites.
        //
        vector<GenPos> ordered_snps;
        for (uint pos = 0; pos < it->second.size(); pos++) {
            loc = it->second[pos];
            t   = psum->locus_tally(loc->id);

            for (uint i = 0; i < loc->snps.size(); i++) {
                col = loc->snps[i]->col;
                if (t->nucs[col].allele_cnt == 2)
                    ordered_snps.push_back(GenPos(loc->id, i, loc->sort_bp(col)));
            }
        }
        sort(ordered_snps.begin(), ordered_snps.end());

        //
        // Now output the genotypes in a separate file for each population.
        //
        for (size_t p=0; p<mpopi.pops().size(); ++p) {
            const Pop& pop = mpopi.pops()[p];

            //
            // Open a markers file containing each marker, its genomic position in basepairs
            // and the two alternative alleles at this position.
            //
            file = out_path + out_prefix + "." + pop.name + "-" + it->first + ".unphased.bgl.markers";

            ofstream mfh(file.c_str(), ofstream::out);
            if (mfh.fail()) {
                cerr << "Error opening Beagle markers file '" << file << "'\n";
                exit(1);
            }
            mfh << "# Stacks v" << VERSION << "; " << " Beagle v3.3; " << date << "\n";

            //
            // Open the genotypes file.
            //
            file = out_path + out_prefix + "." + pop.name + "-" + it->first + ".unphased.bgl";

            ofstream fh(file.c_str(), ofstream::out);
            if (fh.fail()) {
                cerr << "Error opening Beagle genotypes file '" << file << "'\n";
                exit(1);
            }
            fh << "# Stacks v" << VERSION << "; " << " Beagle v3.3; " << date << "\n";

            char p_allele, q_allele;
            //
            // Output a list of all the samples in this population.
            //
            fh << "I\tid";
            for (size_t j = pop.first_sample; j <= pop.last_sample; j++)
                fh << "\t" << mpopi.samples()[j].name << "\t" << mpopi.samples()[j].name;
            fh << "\n";

            //
            // Output population IDs for each sample.
            //
            fh << "S\tid";
            for (size_t j = pop.first_sample; j <= pop.last_sample; j++)
                fh << "\t" << pop.name << "\t" << pop.name;
            fh << "\n";

            //
            // For each marker, output the genotypes for each sample in two successive columns.
            //
            for (uint pos = 0; pos < ordered_snps.size(); pos++) {
                loc = catalog[ordered_snps[pos].id];

                s   = psum->locus(loc->id);
                d   = pmap->locus(loc->id);
                t   = psum->locus_tally(loc->id);
                col = loc->snps[ordered_snps[pos].snp_index]->col;

                //
                // If this site is fixed in all populations or has too many alleles don't output it.
                //
                if (t->nucs[col].allele_cnt != 2)
                    continue;

                //
                // If this site is monomorphic in this population don't output it.
                //
                if (s[p]->nucs[col].pi == 0.0)
                    continue;

                //
                // Output this locus to the markers file.
                //
                mfh << loc->id << "_" << col << "\t"
                    << loc->sort_bp(col) +1  << "\t"
                    << t->nucs[col].p_allele << "\t"
                    << t->nucs[col].q_allele << "\n";

                fh << "M" << "\t" << loc->id << "_" << col;

                for (size_t j = pop.first_sample; j <= pop.last_sample; j++) {
                    //
                    // Output the p allele
                    //
                    if (s[p]->nucs[col].incompatible_site ||
                        s[p]->nucs[col].filtered_site) {
                        //
                        // This site contains more than two alleles in this population or was filtered
                        // due to a minor allele frequency that is too low.
                        //
                        fh << "\t" << "?";

                    } else if (d[j] == NULL || col >= uint(d[j]->len)) {
                        //
                        // Data does not exist.
                        //
                        fh << "\t" << "?";
                    } else if (d[j]->model[col] == 'U') {
                        //
                        // Data exists, but the model call was uncertain.
                        //
                        fh << "\t" << "?";
                    } else {
                        //
                        // Tally up the nucleotide calls.
                        //
                        tally_observed_haplotypes(d[j]->obshap, ordered_snps[pos].snp_index, p_allele, q_allele);

                        if (p_allele == 0 && q_allele == 0)
                            fh << "\t" << "?";
                        else if (p_allele == 0)
                            fh << "\t" << q_allele;
                        else
                            fh << "\t" << p_allele;
                    }

                    //
                    // Now output the q allele
                    //
                    if (s[p]->nucs[col].incompatible_site ||
                        s[p]->nucs[col].filtered_site) {
                        fh << "\t" << "?";

                    } else if (d[j] == NULL || col >= uint(d[j]->len)) {
                        fh << "\t" << "?";

                    } else if (d[j]->model[col] == 'U') {
                        fh << "\t" << "?";

                    } else {
                        if (p_allele == 0 && q_allele == 0)
                            fh << "\t" << "?";
                        else if (q_allele == 0)
                            fh << "\t" << p_allele;
                        else
                            fh << "\t" << q_allele;
                    }
                }
                fh << "\n";
            }

            fh.close();
            mfh.close();
        }
    }

    cerr << "done.\n";

    return 0;
}

int
write_beagle_phased(map<int, CSLocus *> &catalog,
                    PopMap<CSLocus> *pmap,
                    PopSum<CSLocus> *psum)
{
    //
    // Write a Beagle file as a set of haplotpyes as defined here:
    //   http://faculty.washington.edu/browning/beagle/beagle.html
    //
    // We will write one file per chromosome.
    //
    cerr << "Writing population data to phased Beagle files...";

    //
    // Obtain the current date.
    //
    time_t     rawtime;
    struct tm *timeinfo;
    char       date[32];
    time(&rawtime);
    timeinfo = localtime(&rawtime);
    strftime(date, 32, "%B %d, %Y", timeinfo);

    map<string, vector<CSLocus *> >::const_iterator it;
    CSLocus  *loc;
    Datum   **d;

    stringstream pop_name;
    string       file;

    for (it = pmap->ordered_loci().begin(); it != pmap->ordered_loci().end(); it++) {

        //
        // We need to determine an ordering for all legitimate loci/SNPs.
        //
        vector<GenPos> ordered_loci;
        for (uint pos = 0; pos < it->second.size(); pos++) {
            loc = it->second[pos];

            if (loc->snps.size() == 0) continue;

            //
            // Check that there aren't too many haplotypes (PHASE has a max of 50).
            //
            if (loc->alleles.size() > 40) continue;

            //
            // Iterate over the population to determine that this subset of the population
            // has data at this locus.
            //
            d = pmap->locus(loc->id);
            for (int j = 0; j < pmap->sample_cnt(); j++) {
                if (d[j] != NULL &&
                    d[j]->obshap.size() > 0 &&
                    d[j]->obshap.size() <= 2) {
                    //
                    // Data exists, and their are the corrent number of haplotypes.
                    //
                    ordered_loci.push_back(GenPos(loc->id, 0, loc->sort_bp(), haplotype));
                    break;
                }
            }
        }
        sort(ordered_loci.begin(), ordered_loci.end());

        //
        // Now output the genotypes in a separate file for each population.
        //

        for (size_t i_pop=0; i_pop<mpopi.pops().size(); ++i_pop) {
            const Pop& pop = mpopi.pops()[i_pop];

            //
            // Open a file for writing the markers: their genomic position in basepairs
            // and the two alternative alleles at this position.
            //
            file = out_path + out_prefix + "." + pop.name + "-" + it->first + ".phased.bgl.markers";

            ofstream mfh(file.c_str(), ofstream::out);
            if (mfh.fail()) {
                cerr << "Error opening Beagle markers file '" << file << "'\n";
                exit(1);
            }
            mfh << "# Stacks v" << VERSION << "; " << " Beagle v3.3; " << date << "\n";

            //
            // Now output the haplotypes in a separate file.
            //
            file = out_path + out_prefix + "." + pop.name + "-" + it->first + ".phased.bgl";

            ofstream fh(file.c_str(), ofstream::out);
            if (fh.fail()) {
                cerr << "Error opening Beagle markers file '" << file << "'\n";
                exit(1);
            }
            fh << "# Stacks v" << VERSION << "; " << " Beagle v3.3; " << date << "\n";

            //
            // Output a list of all the samples in the data set.
            //
            fh << "I\tid";
            for (size_t j = pop.first_sample; j <= pop.last_sample; j++)
                fh << "\t" << mpopi.samples()[j].name << "\t" << mpopi.samples()[j].name;
            fh << "\n";

            //
            // Output population IDs for each sample.
            //
            fh << "S\tid";
            for (size_t j = pop.first_sample; j <= pop.last_sample; j++)
                fh << "\t" << pop.name << "\t" << pop.name;
            fh << "\n";

            for (uint pos = 0; pos < ordered_loci.size(); pos++) {
                loc = catalog[ordered_loci[pos].id];
                d   = pmap->locus(loc->id);

                //
                // If this locus is monomorphic in this population don't output it.
                //
                set<string> haplotypes;
                for (size_t j = pop.first_sample; j <= pop.last_sample; j++) {
                    if (d[j] == NULL) continue;

                    if (d[j]->obshap.size() == 2) {
                        haplotypes.insert(d[j]->obshap[0]);
                        haplotypes.insert(d[j]->obshap[1]);
                    } else {
                        haplotypes.insert(d[j]->obshap[0]);
                    }
                }
                if (haplotypes.size() == 1) continue;

                //
                // Output this locus to the markers file.
                //
                mfh << loc->id << "\t"
                    << loc->sort_bp() +1;
                for (uint j = 0; j < loc->strings.size(); j++)
                    mfh << "\t" << loc->strings[j].first;
                mfh << "\n";

                //
                // For each marker, output the genotypes for each sample in two successive columns.
                //
                fh << "M" << "\t" << loc->id;

                for (size_t j = pop.first_sample; j <= pop.last_sample; j++) {
                    //
                    // Output the p and the q haplotype
                    //
                    if (d[j] == NULL) {
                        //
                        // Data does not exist.
                        //
                        fh << "\t" << "?" << "\t" << "?";
                    } else {
                        //
                        // Data exists, output the first haplotype. We will assume the haplotypes are
                        // numbered by their position in the loc->strings vector.
                        //
                        if (d[j]->obshap.size() > 2)
                            fh << "\t" << "?" << "\t" << "?";
                        else if (d[j]->obshap.size() == 2)
                            fh << "\t" << d[j]->obshap[0] << "\t" << d[j]->obshap[1];
                        else
                            fh << "\t" << d[j]->obshap[0] << "\t" << d[j]->obshap[0];
                    }
                }
                fh << "\n";
            }
            fh.close();
            mfh.close();
        }
    }

    cerr << "done.\n";

    return 0;
}

int
write_phylip(map<int, CSLocus *> &catalog,
             PopMap<CSLocus> *pmap,
             PopSum<CSLocus> *psum)
{
    //
    // We want to find loci where each locus is fixed within a population but variable between populations.
    //
    // We will write those loci to a Phylip file as defined here:
    //     http://evolution.genetics.washington.edu/phylip/doc/main.html#inputfiles
    //
    string file = out_path + out_prefix + ".phylip";

    cerr << "Writing population data to Phylip file '" << file << "'; ";

    ofstream fh(file.c_str(), ofstream::out);

    if (fh.fail()) {
        cerr << "Error opening Phylip file '" << file << "'\n";
        exit(1);
    }

    file += ".log";

    cerr << "logging nucleotide positions to '" << file << "'...";

    ofstream log_fh(file.c_str(), ofstream::out);

    if (log_fh.fail()) {
        cerr << "Error opening Phylip Log file '" << file << "'\n";
        exit(1);
    }

    //
    // Obtain the current date.
    //
    time_t     rawtime;
    struct tm *timeinfo;
    char       date[32];
    time(&rawtime);
    timeinfo = localtime(&rawtime);
    strftime(date, 32, "%B %d, %Y", timeinfo);

    log_fh << "# Stacks v" << VERSION << "; " << " Phylip sequential; " << date << "\n"
           << "# Seq Pos\tLocus ID\tColumn\tPopulation\n";

    map<string, vector<CSLocus *> >::const_iterator it;
    CSLocus  *loc;
    LocSum  **s;
    LocTally *t;

    int  pop_cnt = psum->pop_cnt();
    char nuc = '?';

    //
    // A map storing, for each population, the concatenated list of interspecific nucleotides.
    //
    map<int, string> interspecific_nucs;

    int index = 0;
    for (it = pmap->ordered_loci().begin(); it != pmap->ordered_loci().end(); it++) {
        for (uint pos = 0; pos < it->second.size(); pos++) {
            loc = it->second[pos];

            s = psum->locus(loc->id);
            t = psum->locus_tally(loc->id);

            for (uint i = 0; i < loc->snps.size(); i++) {
                uint col = loc->snps[i]->col;

                if (phylip_var == false) {
                    //
                    // We are looking for loci that are fixed within each population, but are
                    // variable between one or more populations.
                    //
                    if (t->nucs[col].fixed == true || t->nucs[col].allele_cnt != 2 || t->nucs[col].pop_cnt < 2)
                        continue;

                    bool fixed_within = true;
                    for (int j = 0; j < pop_cnt; j++) {
                        if (s[j]->nucs[col].num_indv == 0)
                            continue;
                        if (s[j]->nucs[col].fixed == false) {
                            fixed_within = false;
                            break;
                        }
                    }
                    if (fixed_within == false) continue;

                    log_fh << index << "\t" << loc->id << "\t" << col << "\t";

                    for (int p=0; p<pop_cnt; p++) {
                        if (s[p]->nucs[col].num_indv > 0) {
                            interspecific_nucs[p] += s[p]->nucs[col].p_nuc;
                            log_fh << mpopi.pops()[p].name << ":" << s[p]->nucs[col].p_nuc << ",";
                        } else {
                            interspecific_nucs[p] += 'N';
                            log_fh << mpopi.pops()[p].name << ":N" << ",";
                        }
                    }
                    log_fh << "\n";
                    index++;

                } else {
                    //
                    // Encode SNPs that are variable within a population as well, using IUPAC notation:
                    //     http://en.wikipedia.org/wiki/Nucleic_acid_notation#IUPAC_notation
                    //
                    if (t->nucs[col].allele_cnt != 2)
                        continue;

                    log_fh << index << "\t" << loc->id << "\t" << col << "\t";

                    for (int j=0; j<pop_cnt; j++) {

                        switch(s[j]->nucs[col].p_nuc) {
                        case 0:
                            nuc = 'N';
                            break;
                        case 'A':
                            switch(s[j]->nucs[col].q_nuc) {
                            case 'C':
                                nuc = 'M';
                                break;
                            case 'G':
                                nuc = 'R';
                                break;
                            case 'T':
                                nuc = 'W';
                                break;
                            case 0:
                                nuc = 'A';
                                break;
                            }
                            break;
                        case 'C':
                            switch(s[j]->nucs[col].q_nuc) {
                            case 'A':
                                nuc = 'M';
                                break;
                            case 'G':
                                nuc = 'S';
                                break;
                            case 'T':
                                nuc = 'Y';
                                break;
                            case 0:
                                nuc = 'C';
                                break;
                            }
                            break;
                        case 'G':
                            switch(s[j]->nucs[col].q_nuc) {
                            case 'A':
                                nuc = 'R';
                                break;
                            case 'C':
                                nuc = 'S';
                                break;
                            case 'T':
                                nuc = 'K';
                                break;
                            case 0:
                                nuc = 'G';
                                break;
                            }
                            break;
                        case 'T':
                            switch(s[j]->nucs[col].q_nuc) {
                            case 'A':
                                nuc = 'W';
                                break;
                            case 'C':
                                nuc = 'Y';
                                break;
                            case 'G':
                                nuc = 'K';
                                break;
                            case 0:
                                nuc = 'T';
                                break;
                            }
                            break;
                        }
                        interspecific_nucs[j] += nuc;
                        log_fh << mpopi.pops()[j].name << ":" << nuc << ",";

                    }
                    log_fh << "\n";
                    index++;
                }
            }
        }
    }

    if (interspecific_nucs.size() == 0) {
        cerr << "  No data is available to write to the Phylip file.\n";
        return 0;
    }

    fh << mpopi.pops().size() << "    " << interspecific_nucs.begin()->second.length() << "\n";
    for (size_t i_pop=0; i_pop<mpopi.pops().size(); ++i_pop) {
        const Pop& pop = mpopi.pops()[i_pop];

        char id_str[id_len];
        sprintf(id_str, "%s", pop.name.c_str());
        uint len = strlen(id_str);
        for (uint j = len; j < 10; j++)
            id_str[j] = ' ';
        id_str[9] = '\0';

        fh << id_str << " " << interspecific_nucs[i_pop] << "\n";
    }

    //
    // Output the header.
    //
    fh << "# Stacks v" << VERSION << "; " << " Phylip sequential; " << date << "\n";

    fh.close();
    log_fh.close();

    cerr << "done.\n";

    return 0;
}

int
write_fullseq_phylip(map<int, CSLocus *> &catalog,
                     PopMap<CSLocus> *pmap,
                     PopSum<CSLocus> *psum)
{
    //
    // We want to write all variable loci in Phylip interleaved format. Polymorphic positions
    // will be encoded using IUPAC notation.
    //
    // We will write those loci to a Phylip file as defined here:
    //     http://evolution.genetics.washington.edu/phylip/doc/main.html#inputfiles
    //
    string file = out_path + out_prefix + ".fullseq.phylip";

    cerr << "Writing full sequence population data to Phylip file '" << file << "'; ";

    ofstream fh(file.c_str(), ofstream::out);

    if (fh.fail()) {
        cerr << "Error opening Phylip file '" << file << "'\n";
        exit(1);
    }

    //
    // We will also write a file that allows us to specify each RAD locus as a separate partition
    // for use in phylogenetics programs.
    //
    file = out_path + out_prefix + ".fullseq.partitions.phylip";

    ofstream par_fh(file.c_str(), ofstream::out);

    if (par_fh.fail()) {
        cerr << "Error opening Phylip partitions file '" << file << "'\n";
        exit(1);
    }

    file = out_path + out_prefix + ".fullseq.phylip.log";
    cerr << "logging nucleotide positions to '" << file << "'...";

    ofstream log_fh(file.c_str(), ofstream::out);

    if (log_fh.fail()) {
        cerr << "Error opening Phylip Log file '" << file << "'\n";
        exit(1);
    }

    //
    // Obtain the current date.
    //
    time_t     rawtime;
    struct tm *timeinfo;
    char       date[32];
    time(&rawtime);
    timeinfo = localtime(&rawtime);
    strftime(date, 32, "%B %d, %Y", timeinfo);

    log_fh << "# Stacks v" << VERSION << "; " << " Phylip interleaved; " << date << "\n"
           << "# Locus ID\tLine Number";
    if (loci_ordered) log_fh << "\tChr\tBasepair";
    log_fh << "\n";

    map<string, vector<CSLocus *> >::const_iterator it;
    CSLocus  *loc;
    LocSum  **s;
    LocTally *t;

    int  pop_cnt = psum->pop_cnt();
    char nuc = '\0';

    bool include;
    uint len = 0;

    //
    // Determine the length of sequence we will output.
    //
    for (it = pmap->ordered_loci().begin(); it != pmap->ordered_loci().end(); it++) {
        for (uint pos = 0; pos < it->second.size(); pos++) {
            loc = it->second[pos];

            t = psum->locus_tally(loc->id);

            include = true;
            for (uint i = 0; i < loc->snps.size(); i++) {
                uint col = loc->snps[i]->col;

                if (t->nucs[col].allele_cnt != 2)
                    include = false;
            }

            if (include)
                len += strlen(loc->con);
        }
    }

    map<int, string> outstrs;
    fh << mpopi.pops().size() << "    " << len << "\n";
    for (size_t i_pop=0; i_pop<mpopi.pops().size(); ++i_pop) {
        const Pop& pop = mpopi.pops()[i_pop];

        char id_str[id_len];
        sprintf(id_str, "%s", pop.name.c_str());
        len = strlen(id_str);
        for (uint j = len; j < 10; j++)
            id_str[j] = ' ';
        id_str[9] = '\0';

        outstrs[i_pop] = string(id_str) + " ";
    }

    char *seq;
    int   line  = 1;
    int   index = 1;
    int   cnt   = 1;

    for (it = pmap->ordered_loci().begin(); it != pmap->ordered_loci().end(); it++) {
        for (uint pos = 0; pos < it->second.size(); pos++) {
            loc = it->second[pos];

            s = psum->locus(loc->id);
            t = psum->locus_tally(loc->id);

            include = true;
            for (uint i = 0; i < loc->snps.size(); i++) {
                uint col = loc->snps[i]->col;

                if (t->nucs[col].allele_cnt != 2)
                    include = false;
            }

            if (!include)
                continue;

            seq = new char[loc->len + 1];
            strcpy(seq, loc->con);

            for (int j = 0; j < pop_cnt; j++) {
                for (uint i = 0; i < loc->snps.size(); i++) {
                    uint col = loc->snps[i]->col;

                    //
                    // Encode SNPs that are variable within a population using IUPAC notation:
                    //     http://en.wikipedia.org/wiki/Nucleic_acid_notation#IUPAC_notation
                    //
                    switch(s[j]->nucs[col].p_nuc) {
                    case 0:
                        nuc = 'N';
                        break;
                    case 'A':
                        switch(s[j]->nucs[col].q_nuc) {
                        case 'C':
                            nuc = 'M';
                            break;
                        case 'G':
                            nuc = 'R';
                            break;
                        case 'T':
                            nuc = 'W';
                            break;
                        case 0:
                            nuc = 'A';
                            break;
                        }
                        break;
                    case 'C':
                        switch(s[j]->nucs[col].q_nuc) {
                        case 'A':
                            nuc = 'M';
                            break;
                        case 'G':
                            nuc = 'S';
                            break;
                        case 'T':
                            nuc = 'Y';
                            break;
                        case 0:
                            nuc = 'C';
                            break;
                        }
                        break;
                    case 'G':
                        switch(s[j]->nucs[col].q_nuc) {
                        case 'A':
                            nuc = 'R';
                            break;
                        case 'C':
                            nuc = 'S';
                            break;
                        case 'T':
                            nuc = 'K';
                            break;
                        case 0:
                            nuc = 'G';
                            break;
                        }
                        break;
                    case 'T':
                        switch(s[j]->nucs[col].q_nuc) {
                        case 'A':
                            nuc = 'W';
                            break;
                        case 'C':
                            nuc = 'Y';
                            break;
                        case 'G':
                            nuc = 'K';
                            break;
                        case 0:
                            nuc = 'T';
                            break;
                        }
                        break;
                    }

                    seq[col] = nuc;
                }

                outstrs[j] += string(seq);
            }
            delete [] seq;

            log_fh << line << "\t" << loc->id;
            if (loci_ordered) log_fh << "\t" << loc->loc.chr() << "\t" << loc->sort_bp() + 1;
            log_fh << "\n";

            for (size_t i_pop=0; i_pop<mpopi.pops().size(); ++i_pop) {
                fh << outstrs[i_pop] << "\n";
                outstrs[i_pop] = "";
                line++;
            }
            fh << "\n";
            line++;

            par_fh << "DNA, p" << cnt << "=" << index << "-" << index + loc->len - 1 << "\n";
            index += loc->len;
            cnt++;
        }
    }

    //
    // Output the header.
    //
    fh << "# Stacks v" << VERSION << "; " << " Phylip interleaved; " << date << "\n";

    fh.close();
    par_fh.close();
    log_fh.close();

    cerr << "done.\n";

    return 0;
}
*/

int
tally_observed_haplotypes(const vector<char *> &obshap, int snp_index, char &p_allele, char &q_allele)
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

    //
    // Determine how many alleles are present at this position in this population.
    // We cannot deal with more than two alternative alleles, if there are more than two
    // in a single population, print a warning and exclude this nucleotide position.
    //
    int i;
    int allele_cnt = 0;
    for (i = 0; i < 4; i++)
        if (nucs[i] > 0) allele_cnt++;

    if (allele_cnt > 2) {
        p_allele = 0;
        q_allele = 0;
        return -1;
    }

    //
    // Record which nucleotide is the P allele and which is the Q allele.
    // (The P allele is the first one alphabetically, and the Q allele the second
    // one, if any.)
    //
    p_allele = 0;
    q_allele = 0;

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

    return 0;
}

int
tally_haplotype_freq(CSLocus *locus, Datum **d, uint sample_cnt,
                     int &total, double &max, string &freq_str)
{
    map<string, double> freq;

    total = 0;
    max   = 0;

    for (uint i = 0; i < sample_cnt; i++) {
        if (d[i] == NULL)
            continue;

        if (d[i]->gtype[0] != '-') {
            freq[d[i]->gtype]++;
            total++;
        }
    }

    if (total == 0)
        return 0;

    double frac;
    stringstream s;
    char   f[id_len];
    map<string, double>::iterator it;
    for (it = freq.begin(); it != freq.end(); it++) {
        frac = (double) it->second / (double) total * 100;
        if (frac > max) max = frac;
        sprintf(f, "(%0.1f%%);", frac);
        s << it->first << ":" << it->second << f;
    }

    freq_str = s.str();

    return 0;
}
