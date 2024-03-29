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

#ifndef __EXSTACKS_H__
#define __EXSTACKS_H__

#ifdef _OPENMP
#include <omp.h>    // OpenMP library
#endif
#include <getopt.h> // Process command-line options
#include <cmath>
#include <cstdlib>
#include <utility>

#include <string>

#include <iostream>
#include <fstream>
#include <sstream>

#include <vector>

#include <map>

#include <set>

#include <queue>
using std::queue;

#include "constants.h"
#include "stacks.h"
#include "sql_utilities.h"

void help( void );
void version( void );
int  parse_command_line(int, char**);
int  parse_tsv(const char *, vector<string> &);
bool compare_dist(pair<int, int>, pair<int, int>);
int  write_simple_output(Locus *);
bool compare_loci(Locus *, Locus *);

#endif // __EXSTACKS_H__
