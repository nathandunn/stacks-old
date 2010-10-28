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

#ifndef __CMB_H__
#define __CMB_H__

#include <math.h>
#include <iostream>
#include <sstream>
using std::stringstream;
using std::cin;
using std::cout;
using std::cerr;
#include <vector>
using std::vector;
#include <string>
using std::string;
#include <set>
using std::set;
#include <utility>
using std::pair;
using std::make_pair;

#include "constants.h"
#include "utils.h"

typedef unsigned int uint;

void   write_cmb(int *, int);

typedef struct cmb {
    int  size;
    int *elem;
} Cmb;

class CombSet {
    //
    // Given these two variables, we will select N choose K combinations. 
    // This combination will be stored in sets, and we will then decrement K by 1
    // and continue to generate sets.
    //
    // Once we have generated all the combinations of a particular size, K, we 
    // will generate all combinations between the different sized sets, creating
    // compound sets.
    //
    int num_elements;  // N elements from which we wish to produce combinations
    int max_set_size;  // maximum set size, K, the largest subset we wish to select.

    int            index;
    vector<int>    size;
    vector<int>    lens;
    vector<int **> sets;
    vector<pair<int, int> > compound_set;
    vector<Cmb *>           compound_comb;

    int      count_valid_comb(long, int, int);
    bool     valid(int *, int);
    int      make_compound_set();
    int    **generate_combinations(int, int, long &, bool);
    int      next_combination(int *, int, int);
    long int num_combinations(int, int);

 public:
    CombSet(int, int);
    ~CombSet();

    Cmb **next(int map[] = NULL);
    void  reset();
    void  destroy(Cmb **);
};

#endif // __CMB_H__
