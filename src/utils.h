// -*-mode:c++; c-style:k&r; c-basic-offset:4;-*-
//
// Copyright 2010-2015, Julian Catchen <jcatchen@illinois.edu>
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

#ifndef __UTILS_H__
#define __UTILS_H__

#include <stdlib.h>
#include <errno.h>
#include <limits.h>
#include <math.h>
#include <iostream>
using std::cerr;
using std::endl;
#include <utility>
using std::pair;
using std::make_pair;

#include <dirent.h>

#include "constants.h"
#include "stacks.h"

char   reverse(char);
char  *rev_comp(const char *);
void   reverse_string(char *);
int    is_integer(const char *);
double is_double(const char *);

double factorial(double);
double reduced_factorial(double, double);

double log_factorial(double);
double reduced_log_factorial(double, double);
//
// Comparison functions for the STL sort routine
//
bool compare_ints(int, int);
bool compare_pair(pair<char, int>, pair<char, int>);
bool compare_pair_intint(pair<int, int>, pair<int, int>);
bool compare_pair_intdouble(pair<int, double>, pair<int, double>);
bool compare_pair_stringint(pair<string, int>, pair<string, int>);
bool compare_pair_snp(pair<string, SNP *>, pair<string, SNP *>);
bool compare_pair_haplotype(pair<string, double>, pair<string, double>);
bool compare_pair_haplotype_rev(pair<string, double>, pair<string, double>);

//
// Comparison classes for STL sets
//
struct int_increasing {
    bool operator() (const int& lhs, const int& rhs) const {
        return lhs < rhs;
    }
};

struct int_decreasing {
    bool operator() (const int& lhs, const int& rhs) const {
        return lhs > rhs;
    }
};

//
// Wrapper for directory parsing functions.
// e.g. for(DirIterator e (path); e; ++e) {...}
//
class DirIterator {
    DIR* dir;
    struct dirent* entry;
public:
    DirIterator(const std::string& dir_path) : dir(NULL), entry(NULL) {
        dir = opendir(dir_path.c_str());
        if (dir == NULL) {
            cerr << "Error: Unable to open directory '" << dir_path << "' for reading.\n";
            throw std::exception();
        }
        entry = readdir(dir);
    }
    ~DirIterator() {closedir(dir);}

    const char* name() const {return entry->d_name;}

    operator bool() const {return entry!=NULL;}
    DirIterator& operator++() {entry = readdir(dir); return *this;}
    dirent* operator*() {return entry;}
};

#endif // __UTILS_H__
