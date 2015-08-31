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

const int max_snp_dist = 500;

class GenPos {
public:
    uint     id;
    uint     bp;
    uint     snp_index;
    loc_type type;

    GenPos(int id, int snp_index, int bp) {
	this->id        = id;
	this->snp_index = snp_index;
	this->bp        = bp;
	this->type      = snp;
    }
    GenPos(int id, int snp_index, int bp, loc_type type) {
	this->id        = id;
	this->snp_index = snp_index;
	this->bp        = bp;
	this->type      = type;
    }
};

void    help( void );
void    version( void );
int     parse_command_line(int, char**);
int     build_file_list(vector<pair<int, string> > &, map<int, pair<int, int> > &, map<int, vector<int> > &);
int     load_marker_list(string, set<int> &);
int     load_marker_column_list(string, map<int, set<int> > &);
int     apply_locus_constraints(map<int, CSLocus *> &, PopMap<CSLocus> *, map<int, pair<int, int> > &, ofstream &);
int     prune_polymorphic_sites(map<int, CSLocus *> &, PopMap<CSLocus> *, PopSum<CSLocus> *, map<int, pair<int, int> > &, map<int, set<int> > &, set<int> &, ofstream &);
int     log_haplotype_cnts(map<int, CSLocus *> &, ofstream &);
bool    order_unordered_loci(map<int, CSLocus *> &);
int     merge_shared_cutsite_loci(map<int, CSLocus *> &, PopMap<CSLocus> *, PopSum<CSLocus> *, map<int, pair<merget, int> > &, ofstream &);
phaset  merge_and_phase_loci(PopMap<CSLocus> *, CSLocus *, CSLocus *, set<int> &, ofstream &);
int     merge_datums(int, int, Datum **, Datum **, set<string> &, int);
int     merge_csloci(CSLocus *, CSLocus *, set<string> &);
int     datum_adjust_snp_positions(map<int, pair<merget, int> > &, CSLocus *, Datum *, map<int, SNPRes *> &);
int     tabulate_haplotypes(map<int, CSLocus *> &, PopMap<CSLocus> *);
int     create_genotype_map(CSLocus *, PopMap<CSLocus> *);
int     call_population_genotypes(CSLocus *, PopMap<CSLocus> *);
int     tally_haplotype_freq(CSLocus *, PopMap<CSLocus> *, int &, double &, string &);
int     translate_genotypes(map<string, string> &, map<string, map<string, string> > &, map<int, CSLocus *> &, PopMap<CSLocus> *, map<int, string> &, set<int> &);
int     correct_fst_bonferroni_win(vector<PopPair *> &);
int     bootstrap_fst_approximate_dist(vector<double> &, vector<int>  &, double *, int *, map<int, vector<double> > &);
int     kernel_smoothed_popstats(map<int, CSLocus *> &, PopMap<CSLocus> *, PopSum<CSLocus> *, int, ofstream &);
int     bootstrap_popstats_approximate_dist(vector<double> &, vector<double> &, vector<int>  &, double *, int *, int, map<int, vector<double> > &, map<int, vector<double> > &);
double  bootstrap_approximate_pval(int, double, map<int, vector<double> > &);

int      calculate_summary_stats(vector<pair<int, string> > &, map<int, pair<int, int> > &, map<int, CSLocus *> &, PopMap<CSLocus> *, PopSum<CSLocus> *);
int      calculate_haplotype_stats(vector<pair<int, string> > &, map<int, pair<int, int> > &, map<int, CSLocus *> &, PopMap<CSLocus> *, PopSum<CSLocus> *);
int      kernel_smoothed_hapstats(vector<CSLocus *> &, PopSum<CSLocus> *, int, double *);
int      calculate_haplotype_divergence(vector<pair<int, string> > &, map<int, pair<int, int> > &, map<int, vector<int> > &, map<int, CSLocus *> &, PopMap<CSLocus> *, PopSum<CSLocus> *);
int      calculate_haplotype_divergence_pairwise(vector<pair<int, string> > &, map<int, pair<int, int> > &, map<int, vector<int> > &, map<int, CSLocus *> &, PopMap<CSLocus> *, PopSum<CSLocus> *);
double   count_haplotypes_at_locus(int, int, Datum **, map<string, double> &);
bool     fixed_locus(map<int, pair<int, int> > &, Datum **, vector<int> &);
bool     uncalled_haplotype(const char *);

int      nuc_substitution_dist(map<string, int> &, double **);
int      nuc_substitution_identity(map<string, int> &, double **);
int      nuc_substitution_identity_max(map<string, int> &, double **);

HapStat *haplotype_amova(map<int, int> &, map<int, pair<int, int> > &, Datum **, LocSum **, vector<int> &);
double   amova_ssd_total(vector<string> &, map<string, int> &, double **);
double   amova_ssd_wp(vector<int> &, map<int, vector<int> > &, map<string, int> &, map<int, vector<string> > &, double **);
double   amova_ssd_ap_wg(vector<int> &, map<int, vector<int> > &, map<string, int> &, map<int, vector<string> > &, double **, double **);
double   amova_ssd_ag(vector<int> &, map<int, vector<int> > &, map<string, int> &, map<int, vector<string> > &, double **, double);

double   haplotype_d_est(map<int, pair<int, int> > &, Datum **, LocSum **, vector<int> &);
LocStat *haplotype_diversity(int, int, Datum **);

int  write_sql(map<int, CSLocus *> &, PopMap<CSLocus> *);
int  write_fst_stats(vector<pair<int, string> > &, map<int, pair<int, int> > &, map<int, CSLocus *> &, PopMap<CSLocus> *, PopSum<CSLocus> *, ofstream &);
int  write_generic(map<int, CSLocus *> &, PopMap<CSLocus> *, map<int, string> &, bool);
int  write_genomic(map<int, CSLocus *> &, PopMap<CSLocus> *);
int  write_fasta(map<int, CSLocus *> &, PopMap<CSLocus> *, map<int, string> &, vector<int> &);
int  write_strict_fasta(map<int, CSLocus *> &, PopMap<CSLocus> *, map<int, string> &, vector<int> &);
int  write_vcf(map<int, CSLocus *> &, PopMap<CSLocus> *, PopSum<CSLocus> *, map<int, string> &, vector<int> &, map<int, pair<merget, int> > &);
int  write_vcf_ordered(map<int, CSLocus *> &, PopMap<CSLocus> *, PopSum<CSLocus> *, map<int, string> &, vector<int> &, map<int, pair<merget, int> > &, ofstream &);
int  write_vcf_haplotypes(map<int, CSLocus *> &, PopMap<CSLocus> *, PopSum<CSLocus> *, map<int, string> &, vector<int> &);
int  populate_snp_calls(map<int, CSLocus *> &, PopMap<CSLocus> *, map<int, string> &, vector<int> &, map<int, pair<merget, int> > &);
int  find_datum_allele_depths(Datum *, int, char, char, int, int &, int &);
int  write_genepop(map<int, CSLocus *> &, PopMap<CSLocus> *, PopSum<CSLocus> *, map<int, pair<int, int> > &, map<int, string> &);
int  write_genepop_ordered(map<int, CSLocus *> &, PopMap<CSLocus> *, PopSum<CSLocus> *, map<int, pair<int, int> > &, map<int, string> &, ofstream &);
int  write_structure(map<int, CSLocus *> &, PopMap<CSLocus> *, PopSum<CSLocus> *, map<int, pair<int, int> > &, map<int, string> &);
int  write_structure_ordered(map<int, CSLocus *> &, PopMap<CSLocus> *, PopSum<CSLocus> *, map<int, pair<int, int> > &, map<int, string> &, ofstream &);
int  write_phase(map<int, CSLocus *> &, PopMap<CSLocus> *, PopSum<CSLocus> *, map<int, pair<int, int> > &, map<int, string> &);
int  write_fastphase(map<int, CSLocus *> &, PopMap<CSLocus> *, PopSum<CSLocus> *, map<int, pair<int, int> > &, map<int, string> &);
int  write_beagle(map<int, CSLocus *> &, PopMap<CSLocus> *, PopSum<CSLocus> *, map<int, pair<int, int> > &, map<int, string> &);
int  write_beagle_phased(map<int, CSLocus *> &, PopMap<CSLocus> *, PopSum<CSLocus> *, map<int, pair<int, int> > &, map<int, string> &);
int  write_plink(map<int, CSLocus *> &, PopMap<CSLocus> *, PopSum<CSLocus> *, map<int, pair<int, int> > &, map<int, string> &);
int  write_hzar(map<int, CSLocus *> &, PopMap<CSLocus> *, PopSum<CSLocus> *, map<int, pair<int, int> > &, map<int, string> &);
int  write_treemix(map<int, CSLocus *> &, PopMap<CSLocus> *, PopSum<CSLocus> *, map<int, pair<int, int> > &, map<int, string> &);
int  write_phylip(map<int, CSLocus *> &, PopMap<CSLocus> *, PopSum<CSLocus> *, map<int, pair<int, int> > &, map<int, string> &);
int  write_fullseq_phylip(map<int, CSLocus *> &, PopMap<CSLocus> *, PopSum<CSLocus> *, map<int, pair<int, int> > &, map<int, string> &);
int  tally_observed_haplotypes(vector<char *> &, int, char &, char &);
int  tally_ref_alleles(LocSum **, int, int, char &, char &);
int  load_snp_calls(string,  PopMap<CSLocus> *);

bool compare_pop_map(pair<int, string>, pair<int, string>);
bool hap_compare(pair<string, int>, pair<string, int>);
bool compare_genpos(GenPos, GenPos);

#endif // __POPULATIONS_H__
