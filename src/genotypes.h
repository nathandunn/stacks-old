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
#include "stacks.h"
#include "PopMap.h"
#include "sql_utilities.h"
#include "genotype_dictionaries.h"

enum map_types {unk, none, gen, dh, cp, bc1, f2};
enum out_types {rqtl, joinmap, genomic};

//
// Catalog Locus Class
//
class CLocus : public Locus {
public:
    CLocus() : Locus() { this->f = 0.0; this->hcnt = 0; this->gcnt = 0; this->trans_gcnt = 0; };
    string annotation;
    string marker;
    double f;                 // Inbreeder's coefficient
    map<string, string> gmap; // Observed haplotype to genotype map for this locus.
    int hcnt;                 // Number of progeny containing a haplotype for this locus.
    int gcnt;                 // Number of progeny containing a valid genotype.
    int trans_gcnt;           // Number of progeny containing a valid 
                              // genotype, translated for a particular map type.
};

void help( void );
void version( void );
int  parse_command_line(int, char**);
int  build_file_list(vector<string> &);
int  load_marker_list(string, set<int> &);
int  identify_parental_ids(map<int, CLocus *> &, set<int> &);
int  reduce_catalog(map<int, CLocus *> &, set<int> &, set<int> &);
int  find_markers(map<int, CLocus *> &, PopMap<CLocus> *, set<int> &);
int  calculate_f(map<int, CLocus *> &, PopMap<CLocus> *, set<int> &);
int  create_genotype_map(CLocus *, PopMap<CLocus> *, set<int> &);
int  call_population_genotypes(CLocus *, PopMap<CLocus> *, map<string, map<string, string> > &);
int  tally_progeny_haplotypes(CLocus *, PopMap<CLocus> *, set<int> &, int &, double &, string &);
int  translate_genotypes(map<string, string> &, map<string, map<string, string> > &, map<int, CLocus *> &, PopMap<CLocus> *, map<int, string> &, set<int> &);

int  automated_corrections(map<int, string> &, set<int> &, map<int, CLocus *> &, vector<vector<CatMatch *> > &, PopMap<CLocus> *);
int  check_uncalled_snps(CLocus *, Locus *, Datum *);
int  call_alleles(vector<SNP *> &, vector<char *> &, vector<string> &);
int  check_homozygosity(vector<char *> &, int, char, char, string &);

int  export_gen_map(map<int, CLocus *> &, PopMap<CLocus> *, set<int> &, map<int, string> &);
int  export_cp_map(map<int, CLocus *> &, PopMap<CLocus> *, set<int> &, map<int, string> &);
int  export_bc1_map(map<int, CLocus *> &, PopMap<CLocus> *, set<int> &, map<int, string> &);
int  export_dh_map(map<int, CLocus *> &, PopMap<CLocus> *, set<int> &, map<int, string> &);
int  export_f2_map(map<int, CLocus *> &, PopMap<CLocus> *, set<int> &, map<int, string> &);

int  write_generic(map<int, CLocus *> &, PopMap<CLocus> *, map<int, string> &, set<int> &, bool);
int  write_sql(map<int, CLocus *> &,     PopMap<CLocus> *, set<int> &);
int  write_joinmap(map<int, CLocus *> &, PopMap<CLocus> *, map<string, string> &, map<int, string> &, set<int> &);
int  write_rqtl(map<int, CLocus *> &,    PopMap<CLocus> *, map<string, string> &, map<int, string> &, set<int> &);
int  write_genomic(map<int, CLocus *> &, PopMap<CLocus> *);

bool hap_compare(pair<string, int>, pair<string, int>);
bool compare_pair(pair<char, int>, pair<char, int>);

#endif // __GENOTYPES_H__
