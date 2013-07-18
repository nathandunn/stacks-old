// -*-mode:c++; c-style:k&r; c-basic-offset:4;-*-
//
// Copyright 2010-2012, Julian Catchen <jcatchen@uoregon.edu>
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

#ifndef __MODELS_H__
#define __MODELS_H__

#include <string.h>
#include <math.h>
#include <string>
using std::string;
#include <vector>
using std::vector;
#include <map>
using std::map;
#include <utility>
using std::pair;
using std::make_pair;
#include <algorithm>
#include <iostream>
using std::ifstream;
using std::cin;
using std::cout;
using std::cerr;
using std::endl;

#include "constants.h"
#include "utils.h"
#include "mstack.h"
#include "locus.h"

//
// Possible models for calling nucleotide positions as fixed or variable
//
enum modelt {fixed, snp, bounded};

//
// For use with the multinomial model to call fixed nucleotides.
//
extern const int barcode_size;
extern modelt model_type;
extern double p_freq;
extern double barcode_err_freq;
extern double bound_low;
extern double bound_high;
extern double heterozygote_limit;
extern double homozygote_limit;

snp_type call_bounded_multinomial_snp(MergedStack *, int, map<char, int> &, bool);
snp_type call_bounded_multinomial_snp(Locus *, int, map<char, int> &);
snp_type call_multinomial_snp(MergedStack *, int, map<char, int> &, bool);
snp_type call_multinomial_snp(Locus *, int, map<char, int> &);
int      call_multinomial_fixed(MergedStack *, int, map<char, int> &);
double   heterozygous_likelihood(int, map<char, int> &);
double   homozygous_likelihood(int, map<char, int> &);

#endif // __MODELS_H__
