// -*-mode:c++; c-style:k&r; c-basic-offset:4;-*-
//
// Copyright 2013, Julian Catchen <jcatchen@uoregon.edu>
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

#ifndef __RXSTACKS_H__
#define __RXSTACKS_H__

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
#include <iomanip> // std::setprecision
#include <sstream>
using std::stringstream;
#include <vector>
using std::vector;
#include <queue>
using std::queue;
#include <map>
using std::map;
#include <set>
using std::set;

#include "constants.h"
#include "stacks.h"
#include "locus.h"
#include "PopMap.h"
#include "PopSum.h"
#include "utils.h"
#include "sql_utilities.h"
#include "models.h"
#include "mst.h"

void   help( void );
void   version( void );
int    parse_command_line(int, char**);
int    build_file_list(vector<pair<int, string> > &);
int    init_log(int, char **, ofstream &, ofstream &, ofstream &);
int    sum_haplotype_counts(map<int, CSLocus *> &, PopMap<CSLocus> *);
int    prune_mst_haplotypes(CSLocus *, Datum *, Locus *, unsigned long int &, ofstream &);
int    prune_locus_haplotypes(CSLocus *, Datum *, Locus *, unsigned long int &, ofstream &);
string convert_catalog_haplotype_to_sample(string, CSLocus *, Locus *);
int    remove_haplotype(CSLocus *, Locus *, string, unsigned long &, ofstream &, string);
int    dist(string, string);
int    measure_error(CSLocus *, Locus *, Datum *, ofstream &);
int    calc_lnl_means(map<int, CSLocus *> &, PopMap<CSLocus> *);
int    prune_nucleotides(CSLocus *, Locus *, Datum *, ofstream &, unsigned long int &,
			 unsigned long int &, unsigned long int &, unsigned long int &,
			 unsigned long int &, unsigned long int &, unsigned long int &);
int    invoke_model(Locus *, int, map<char, int> &);
int    call_alleles(Locus *, set<int> &); 
int    generate_matched_haplotypes(CSLocus *, Locus *, Datum *);
int    fill_catalog_snps(map<int, CSLocus *> &);
int    log_model_calls(Locus *, ofstream &,
		       unsigned long int &, unsigned long int &, unsigned long int &,
		       unsigned long int &, unsigned long int &, unsigned long int &);
int    write_results(string, map<int, Locus *> &);

#endif // __RXSTACKS_H__
