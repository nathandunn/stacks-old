// -*-mode:c++; c-style:k&r; c-basic-offset:4;-*-
//
// Copyright (c) 2014 University of Oregon
// Created by Julian Catchen <jcatchen@uoregon.edu>
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

char   reverse(char);
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
bool compare_pair_intdouble(pair<int, double>, pair<int, double>);
bool compare_pair_snp(pair<string, SNP *>, pair<string, SNP *>);


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

#endif // __UTILS_H__
