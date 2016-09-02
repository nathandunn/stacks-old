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

#include "stacks.h"

const char nt_A = 'A';
const char nt_C = 'C';
const char nt_G = 'G';
const char nt_T = 'T';
const char nt_N = 'N';

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

/* NtCounts
 * ========== */
class NtCounts {
public:
    NtCounts() : ar_{0,0,0,0,0} {}

    void add(char c);
    uint count(char nt) const;
    uint coverage() const {return nA() + nC() + nG() + nT();}
    char consensus() const;

    uint nA() const {return ar_[0];}
    uint nC() const {return ar_[1];}
    uint nG() const {return ar_[2];}
    uint nT() const {return ar_[3];}
    uint nN() const {return ar_[4];}

    const uint* array() const {return ar_;}

protected:
    uint ar_[5];

    uint& nA() {return ar_[0];}
    uint& nC() {return ar_[1];}
    uint& nG() {return ar_[2];}
    uint& nT() {return ar_[3];}
    uint& nN() {return ar_[4];}
};

inline
void NtCounts::add(char c) {
    if (c == nt_A)
        ++nA();
    else if (c == nt_C)
        ++nC();
    else if (c == nt_G)
        ++nG();
    else if (c == nt_T)
        ++nT();
    else if (c == nt_N)
        ++nN();
    else
        throw std::exception();
}

inline
uint NtCounts::count(char c) const {
    if (c == nt_A)
        return nA();
    else if (c == nt_C)
        return nC();
    else if (c == nt_G)
        return nG();
    else if (c == nt_T)
        return nT();
    else if (c == nt_N)
        return nN();
    throw std::exception();
}

inline
char NtCounts::consensus() const {
    char c = nt_A;
    uint n = nA();
    if (n < nC()) {
        c = nt_C;
        n = nC();
    }
    if (n < nG()) {
        c = nt_G;
        n = nG();
    }
    if (n < nT()) {
        c = nt_T;
    }
    return c;
}

#endif // __UTILS_H__
