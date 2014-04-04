// -*-mode:c++; c-style:k&r; c-basic-offset:4;-*-
//
// Copyright 2010, Julian Catchen <jcatchen@uoregon.edu>
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

#ifndef __GENOTYPES_H__
#define __GENOTYPES_H__

#ifdef _OPENMP
#include <omp.h>    // OpenMP library
#endif
#include <getopt.h> // Process command-line options
#include <dirent.h> // Open/Read contents of a directory
#include <math.h>
#include <stdlib.h>
#include <string.h>
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
#include "input.h"
#include "stacks.h"
#include "locus.h"
#include "renz.h"
#include "PopMap.h"
#include "sql_utilities.h"
#include "catalog_utils.h"
#include "utils.h"
#include "genotype_dictionaries.h"

//
// Chi-squared distribution critical values
//   [0] => p-values
//   [1] => one degree of freedom
//   [2] => two degrees of freedom
//   [3] => three degrees of freedom
//
const int chisq_crit_values_size = 10;
const double chisq_crit_values[4][10] = {
    {0.50, 0.25, 0.20, 0.15, 0.10, 0.05,  0.01,  0.005,  0.001,  0.0005},
    {0.46, 1.32, 1.64, 2.07, 2.71, 3.84,  6.63,  7.880, 10.830, 12.1200},
    {1.39, 2.77, 3.22, 3.79, 4.61, 5.99,  9.21, 10.600, 13.820, 15.2000},
    {2.37, 4.11, 4.64, 5.32, 6.25, 7.81, 11.34, 12.840, 16.270, 17.7300}
};

void   help( void );
void   version( void );
int    parse_command_line(int, char**);
int    build_file_list(vector<string> &);
int    load_marker_list(string, set<int> &);
int    identify_parental_ids(map<int, CSLocus *> &, vector<int> &, set<int> &);
int    find_markers(map<int, CSLocus *> &, PopMap<CSLocus> *, set<int> &);
int    calculate_f(map<int, CSLocus *> &, PopMap<CSLocus> *, set<int> &);
int    create_genotype_map(CSLocus *, PopMap<CSLocus> *, set<int> &);
int    apply_locus_constraints(map<int, CSLocus *> &, PopMap<CSLocus> *);
int    call_population_genotypes(CSLocus *, PopMap<CSLocus> *, map<string, map<string, string> > &);
int    tally_progeny_haplotypes(CSLocus *, PopMap<CSLocus> *, set<int> &, int &, double &, string &);
int    map_specific_genotypes(map<int, CSLocus *> &, PopMap<CSLocus> *, set<int> &);
int    translate_genotypes(map<string, string> &, map<string, map<string, string> > &, map<int, CSLocus *> &, PopMap<CSLocus> *, map<int, string> &, set<int> &);
int    automated_corrections(map<int, string> &, set<int> &, map<int, CSLocus *> &, vector<vector<CatMatch *> > &, PopMap<CSLocus> *);
int    check_uncalled_snps(CSLocus *, Locus *, Datum *);
int    call_alleles(vector<SNP *> &, vector<char *> &, vector<string> &);
int    check_homozygosity(vector<char *> &, int, char, char, string &);
int    manual_corrections(string, PopMap<CSLocus> *);

int    correct_cp_markers_missing_alleles(set<int> &, map<int, CSLocus *> &, PopMap<CSLocus> *);
int    calc_segregation_distortion(map<string, map<string, double> > &, map<int, CSLocus *> &, PopMap<CSLocus> *, set<int> &); 
double tally_generic_gtypes(int, PopMap<CSLocus> *, set<int> &, map<string, int> &);
double tally_translated_gtypes(int, PopMap<CSLocus> *, set<int> &, map<string, string> &, map<string, int> &);
double chisq_test(map<string, map<string, double> > &, map<string, int> &, string, double);
double chisq_pvalue(int, double);

int  export_gen_map(map<int, CSLocus *> &, PopMap<CSLocus> *, set<int> &, map<int, string> &);
int  export_cp_map(map<int, CSLocus *> &,  PopMap<CSLocus> *, set<int> &, map<int, string> &);
int  export_bc1_map(map<int, CSLocus *> &, PopMap<CSLocus> *, set<int> &, map<int, string> &);
int  export_dh_map(map<int, CSLocus *> &,  PopMap<CSLocus> *, set<int> &, map<int, string> &);
int  export_f2_map(map<int, CSLocus *> &,  PopMap<CSLocus> *, set<int> &, map<int, string> &);

int  write_generic(map<int, CSLocus *> &, PopMap<CSLocus> *, map<int, string> &, set<int> &, bool);
int  write_sql(map<int, CSLocus *> &,     PopMap<CSLocus> *, set<int> &);
int  write_joinmap(map<int, CSLocus *> &, PopMap<CSLocus> *, map<string, string> &, map<int, string> &, set<int> &);
int  write_onemap(map<int, CSLocus *> &,  PopMap<CSLocus> *, map<string, string> &, map<int, string> &, set<int> &);
int  write_onemap_mapmaker(map<int, CSLocus *> &,  PopMap<CSLocus> *, map<string, string> &, map<int, string> &, set<int> &);
int  write_rqtl(map<int, CSLocus *> &,    PopMap<CSLocus> *, map<string, string> &, map<int, string> &, set<int> &);
int  write_genomic(map<int, CSLocus *> &, PopMap<CSLocus> *);

bool hap_compare(pair<string, int>, pair<string, int>);

#endif // __GENOTYPES_H__
