// -*-mode:c++; c-style:k&r; c-basic-offset:4;-*-
//
// Copyright 2012-2015, Julian Catchen <jcatchen@illinois.edu>
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

#ifndef __POPULATIONS_H__
#define __POPULATIONS_H__

#ifdef _OPENMP
#include <omp.h>    // OpenMP library
#endif
#include <getopt.h> // Process command-line options
#include <dirent.h> // Open/Read contents of a directory
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <utility>
using std::pair;
using std::make_pair;
#include <string>
using std::string;
#include <iostream>
#include <fstream>
using std::ifstream;
using std::ofstream;
using std::cin;
using std::cout;
using std::cerr;
using std::endl;
#include <iomanip>
using std::setw;
using std::setprecision;
using std::fixed;
#include <sstream>
using std::stringstream;
#include <vector>
using std::vector;
#include <map>
using std::map;
#include <set>
using std::set;

#include "constants.h"
#include "stacks.h"
#include "locus.h"
#include "renz.h"
#include "PopMap.h"
#include "PopSum.h"
#include "utils.h"
#include "log_utils.h"
#include "catalog_utils.h"
#include "sql_utilities.h"
#include "genotype_dictionaries.h"
#include "ordered.h"
#include "smoothing.h"
#include "bootstrap.h"

enum corr_type {p_value, bonferroni_win, bonferroni_gen, no_correction};
enum bs_type   {bs_exact, bs_approx, bs_none};
enum merget    {merge_sink, merge_src};
enum phaset    {merge_failure, simple_merge, complex_phase, nomapping_fail, multimapping_fail, multiple_fails};
enum class InputMode {stacks, vcf};

const int max_snp_dist = 500;

void    help( void );
void    version( void );
int     parse_command_line(int, char**);
int     build_file_list();
int     load_marker_list(string, set<int> &);
int     load_marker_column_list(string, map<int, set<int> > &);
int     apply_locus_constraints(map<int, CSLocus *> &, PopMap<CSLocus> *, ofstream &);
int     prune_polymorphic_sites(map<int, CSLocus *> &, PopMap<CSLocus> *, PopSum<CSLocus> *, map<int, set<int> > &, set<int> &, ofstream &);
int     log_haplotype_cnts(map<int, CSLocus *> &, ofstream &);
bool    order_unordered_loci(map<int, CSLocus *> &);
int     merge_shared_cutsite_loci(map<int, CSLocus *> &, PopMap<CSLocus> *, PopSum<CSLocus> *, map<int, pair<merget, int> > &, ofstream &);
phaset  merge_and_phase_loci(PopMap<CSLocus> *, CSLocus *, CSLocus *, set<int> &, ofstream &);
int     merge_datums(int, int, Datum **, Datum **, set<string> &, int);
int     merge_csloci(CSLocus *, CSLocus *, set<string> &);
int     tabulate_haplotypes(map<int, CSLocus *> &, PopMap<CSLocus> *);
int     create_genotype_map(CSLocus *, PopMap<CSLocus> *);
int     call_population_genotypes(CSLocus *, PopMap<CSLocus> *);
int     translate_genotypes(map<string, string> &, map<string, map<string, string> > &, map<int, CSLocus *> &, PopMap<CSLocus> *, map<int, string> &, set<int> &); // This function doesn't exist (March 24, 2016)
int     correct_fst_bonferroni_win(vector<PopPair *> &);
int     bootstrap_fst_approximate_dist(vector<double> &, vector<int>  &, double *, int *, map<int, vector<double> > &); // not used (March 23, 2016)
int     kernel_smoothed_popstats(map<int, CSLocus *> &, PopMap<CSLocus> *, PopSum<CSLocus> *, int, ofstream &);
int     bootstrap_popstats_approximate_dist(vector<double> &, vector<double> &, vector<int>  &, double *, int *, int, map<int, vector<double> > &, map<int, vector<double> > &); // not used (March 23, 2016)
double  bootstrap_approximate_pval(int, double, map<int, vector<double> > &);

int      calculate_summary_stats(map<int, CSLocus *> &, PopMap<CSLocus> *, PopSum<CSLocus> *);
int      calculate_haplotype_stats(map<int, CSLocus *> &, PopMap<CSLocus> *, PopSum<CSLocus> *);
int      kernel_smoothed_hapstats(vector<CSLocus *> &, PopSum<CSLocus> *, int, double *);
int      calculate_haplotype_divergence(map<int, CSLocus *> &, PopMap<CSLocus> *, PopSum<CSLocus> *);
int      calculate_haplotype_divergence_pairwise(map<int, CSLocus *> &, PopMap<CSLocus> *, PopSum<CSLocus> *);
bool     fixed_locus(Datum **, vector<int> &);

int      nuc_substitution_dist(map<string, int> &, double **);
int      nuc_substitution_identity(map<string, int> &, double **);
int      nuc_substitution_identity_max(map<string, int> &, double **);

HapStat *haplotype_amova(map<int, int> &, Datum **, LocSum **, vector<int> &);
double   amova_ssd_total(vector<string> &, map<string, int> &, double **);
double   amova_ssd_wp(vector<int> &, map<int, vector<int> > &, map<string, int> &, map<int, vector<string> > &, double **);
double   amova_ssd_ap_wg(vector<int> &, map<int, vector<int> > &, map<string, int> &, map<int, vector<string> > &, double **, double **);
double   amova_ssd_ag(vector<int> &, map<int, vector<int> > &, map<string, int> &, map<int, vector<string> > &, double **, double);

double   haplotype_d_est(Datum **, LocSum **, vector<int> &);
LocStat *haplotype_diversity(int, int, Datum **);
double   count_haplotypes_at_locus(int, int, Datum**, map<string, double>&);

//int  tally_ref_alleles(LocSum **, int, int, char &, char &); //unused; also commented out in the .cc
//int  load_snp_calls(string,  PopMap<CSLocus> *); //no implementation

//bool compare_pop_map(pair<int, string>, pair<int, string>); //no implementation; the function is in [sql_utilities.h]
bool hap_compare(pair<string, int>, pair<string, int>);

inline
bool uncalled_haplotype(const char *haplotype)
{
    for (const char *p = haplotype; *p != '\0'; p++)
        if (*p == 'N' || *p == 'n')
            return true;
    return false;
}

inline
double count_haplotypes_at_locus(int start, int end, Datum **d, map<string, double> &hap_cnts)
{
    double n = 0.0; // n.b. Unless hap_cnts was not empty, this is [hap_cnts.size()].

    for (int i = start; i <= end; i++) {
        if (d[i] == NULL)
            continue;

        const vector<char*>& haps = d[i]->obshap;
        if (haps.size() > 2) {
            // Too many haplotypes, ignored.
            continue;

        } else if (haps.size() == 1) {
            // Homozygote.
            if(!uncalled_haplotype(d[i]->obshap[0])) {
                n += 2;
                hap_cnts[d[i]->obshap[0]] += 2;
            }
        } else {
            // Heterozygote.
            for (uint j = 0; j < d[i]->obshap.size(); j++) {
                if(!uncalled_haplotype(d[i]->obshap[0])) {
                    n++;
                    hap_cnts[d[i]->obshap[j]]++;
                }
            }
        }
    }

    return n;
}

#endif // __POPULATIONS_H__
