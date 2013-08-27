// -*-mode:c++; c-style:k&r; c-basic-offset:4;-*-
//
// Copyright 2012, Julian Catchen <jcatchen@uoregon.edu>
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
#include "sql_utilities.h"
#include "genotype_dictionaries.h"

enum map_types {unk, none, gen, dh, cp, bc1, f2};
enum out_types {rqtl, joinmap, genomic};
enum corr_type {p_value, bonferroni_win, bonferroni_gen, no_correction};
enum bs_type   {bs_exact, bs_approx, bs_none};

const int max_snp_dist = 500;

//
// Bootstrap resamplign structure.
//
class BSample {
public:
    int    bp;
    int    alleles;
    double f;
    double pi;

    BSample() {
	this->bp      = 0;
	this->alleles = 0;
	this->f       = 0.0;
	this->pi      = -1.0;
    }
};

void    help( void );
void    version( void );
int     parse_command_line(int, char**);
int     build_file_list(vector<pair<int, string> > &, map<int, pair<int, int> > &);
int     load_marker_list(string, set<int> &);
int     reduce_catalog(map<int, CSLocus *> &, set<int> &, set<int> &);
int     apply_locus_constraints(map<int, CSLocus *> &, PopMap<CSLocus> *, map<int, pair<int, int> > &);
bool    order_unordered_loci(map<int, CSLocus *> &);
int     tabulate_haplotypes(map<int, CSLocus *> &, PopMap<CSLocus> *);
int     create_genotype_map(CSLocus *, PopMap<CSLocus> *);
int     call_population_genotypes(CSLocus *, PopMap<CSLocus> *);
int     tally_haplotype_freq(CSLocus *, PopMap<CSLocus> *, int &, double &, string &);
int     translate_genotypes(map<string, string> &, map<string, map<string, string> > &, map<int, CSLocus *> &, PopMap<CSLocus> *, map<int, string> &, set<int> &);
int     correct_fst_bonferroni_win(vector<PopPair *> &);
int     init_chr_pairs(map<string, vector<PopPair *> > &, string, map<uint, uint> &, vector<CSLocus *> &);
int     kernel_smoothed_fst(vector<PopPair *> &, double *, int *);
int     bootstrap_fst(vector<double> &, vector<PopPair *> &, double *);
int     bootstrap_fst_approximate_dist(vector<double> &, vector<int>  &, double *, int *, map<int, vector<double> > &);
int     init_chr_sites(vector<SumStat *> &, int, vector<CSLocus *> &, PopSum<CSLocus> *, uint &, ofstream &);
int     kernel_smoothed_popstats(map<int, CSLocus *> &, PopMap<CSLocus> *, PopSum<CSLocus> *, int, ofstream &);
int     bootstrap_popstats(vector<double> &, vector<double> &, vector<SumStat *> &, int, int, double *, SumStat *); 
int     bootstrap_popstats_approximate_dist(vector<double> &, vector<double> &, vector<int>  &, double *, int *, int, map<int, vector<double> > &, map<int, vector<double> > &);
double  bootstrap_pval(double, vector<double> &);
double  bootstrap_approximate_pval(int, double, map<int, vector<double> > &);
double *calculate_weights(void);

int  write_sql(map<int, CSLocus *> &, PopMap<CSLocus> *);
int  write_summary_stats(vector<pair<int, string> > &, map<int, pair<int, int> > &, map<int, CSLocus *> &, PopMap<CSLocus> *, PopSum<CSLocus> *);
int  write_fst_stats(vector<pair<int, string> > &, map<int, pair<int, int> > &, map<int, CSLocus *> &, PopMap<CSLocus> *, PopSum<CSLocus> *, ofstream &);
int  write_linkage_stats(map<int, pair<int, int> > &, map<int, CSLocus *> &, PopMap<CSLocus> *, PopSum<CSLocus> *);
int  write_generic(map<int, CSLocus *> &, PopMap<CSLocus> *, map<int, string> &, bool);
int  write_genomic(map<int, CSLocus *> &, PopMap<CSLocus> *);
int  write_vcf(map<int, CSLocus *> &, PopMap<CSLocus> *, PopSum<CSLocus> *, map<int, string> &, vector<int> &);
int  write_genepop(map<int, CSLocus *> &, PopMap<CSLocus> *, PopSum<CSLocus> *, map<int, pair<int, int> > &, map<int, string> &);
int  write_structure(map<int, CSLocus *> &, PopMap<CSLocus> *, PopSum<CSLocus> *, map<int, pair<int, int> > &, map<int, string> &);
int  write_phase(map<int, CSLocus *> &, PopMap<CSLocus> *, PopSum<CSLocus> *, map<int, pair<int, int> > &, map<int, string> &);
int  write_beagle(map<int, CSLocus *> &, PopMap<CSLocus> *, PopSum<CSLocus> *, map<int, pair<int, int> > &, map<int, string> &);
int  write_plink(map<int, CSLocus *> &, PopMap<CSLocus> *, PopSum<CSLocus> *, map<int, pair<int, int> > &, map<int, string> &);
int  write_phylip(map<int, CSLocus *> &, PopMap<CSLocus> *, PopSum<CSLocus> *, map<int, pair<int, int> > &, map<int, string> &);
int  tally_observed_haplotypes(vector<char *> &, int, char &, char &);
int  tally_ref_alleles(LocSum **, int, int, char &, char &);
int  load_snp_calls(string,  PopMap<CSLocus> *);

bool compare_pop_map(pair<int, string>, pair<int, string>);
bool hap_compare(pair<string, int>, pair<string, int>);

#endif // __POPULATIONS_H__
