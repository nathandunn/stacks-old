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

#ifndef __PHASEDSTACKS_H__
#define __PHASEDSTACKS_H__

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
#include <iomanip>
#include <vector>
using std::vector;
#include <map>
using std::map;
#include <set>
using std::set;

#include "constants.h"
#include "utils.h"

class Sample {
public:
    string name;
    int    id;
    int    size;
    char  *nucs_1;
    char  *nucs_2;

    Sample() {
	this->id   = 0;
	this->size = 0;
	this->nucs_1 = NULL;
	this->nucs_2 = NULL;
    }
    ~Sample() {
	if (this->nucs_1 != NULL)
	    delete [] this->nucs_1;
	if (this->nucs_2 != NULL)
	    delete [] this->nucs_2;
    }
};

class NucSum {
public:
    uint  bp;
    uint  clocus;
    uint  nuc[4];
    float freq;

    // nuc[0] == A
    // nuc[1] == C
    // nuc[2] == G
    // nuc[3] == T

    NucSum() {
	this->freq   = 1.0;
	this->bp     = 0;
	this->clocus = 0;
	for (uint i = 0; i < 4; i++)
	    this->nuc[i] = 0;
    }
};

class PhasedSummary {
    map<string, int> sample_map;

public:
    uint     size;
    uint     sample_cnt;
    NucSum  *nucs;
    Sample  *samples;
    float  **dprime;

    PhasedSummary(uint num_samples, uint num_genotypes) {
	this->sample_cnt = num_samples;
	this->samples    = new Sample[this->sample_cnt];
	this->size       = num_genotypes;
	this->nucs       = new NucSum[this->size];
	this->dprime     = new float *[this->size];
	for (uint i = 0; i < this->size; i++) {
	    this->dprime[i] = new float[this->size];
	    memset(this->dprime[i], 0, this->size);
	}
    }
    ~PhasedSummary() {
	if (this->nucs != NULL)
	    delete [] this->nucs;
	if (this->dprime != NULL) {
	    for (uint i = 0; i < this->size; i++)
		delete [] this->dprime[i];
	    delete [] this->dprime;
	}
	if (this->samples != NULL)
	    delete [] this->samples;
    }
    int add_sample(string name) {
	uint i = this->sample_map.size();
	this->sample_map[name] = i;
	this->samples[i].name = name;
	return i;
    }
};

void  help( void );
void  version( void );
int   parse_command_line(int, char**);
int   build_file_list(vector<pair<int, string> > &);
PhasedSummary *parse_phase(string);
int   summarize_phased_genotypes(PhasedSummary *);
int   calc_dprime(PhasedSummary *);
int   assign_alleles(NucSum, char &, char &, float &, float &);
int   write_dprime(string, PhasedSummary *);

#endif // __PHASEDSTACKS_H__
