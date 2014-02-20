// -*-mode:c++; c-style:k&r; c-basic-offset:4;-*-
//
// Copyright (c) 2014 University of Oregon
// Created by Julian Catchen <jcatchen@uoregon.edu>
//

//
// utils.cc -- common routines needed in multiple object files.
//
// Julian Catchen
// jcatchen@uoregon.edu
// University of Oregon
//
#include "utils.h"

char reverse(char c) {
    switch (c) {
    case 'A':
    case 'a':
	return 'T';
	break;
    case 'C':
    case 'c':
	return 'G';
	break;
    case 'G':
    case 'g':
	return 'C';
	break;
    case 'T':
    case 't':
	return 'A';
	break;
    case 'N':
    case 'n':
    case '.':
	return 'N';
	break;
    }

    return 'N';
}

int is_integer(const char *str) {
    //
    // Adapted from the strtol manpage.
    //
    char *endptr;

    // To distinguish success/failure after call
    errno = 0;
    long val = strtol(str, &endptr, 10);

    //
    // Check for various possible errors
    //
    if ((errno == ERANGE && (val == LONG_MAX || val == LONG_MIN))
	|| (errno != 0 && val == 0)) {
	return -1;
    }

    if (endptr == str || *endptr != '\0')
	return -1;

    return (int) val;
}

double is_double(const char *str) {
    //
    // Adapted from the strtol manpage.
    //
    char *endptr;

    // To distinguish success/failure after call
    errno = 0;
    double val = strtod(str, &endptr);

    //
    // Check for various possible errors
    //
    if ((errno == ERANGE && (val == LONG_MAX || val == LONG_MIN))
	|| (errno != 0 && val == 0)) {
	return -1;
    }

    if (endptr == str || *endptr != '\0')
	return -1;

    return val;
}

double factorial(double n) {
    double fact = 1;

    for (double i = n; i > 1; i--)
        fact *= i;

    return fact;
}

double reduced_factorial(double n, double d) {
    double f = n - d;

    if (f < 0) 
        return 0;
    else if (f == 0) 
        return 1;
    else if (f == 1)
        return n;

    f = n;
    n--;
    while (n > d) {
        f *= n;
        n--;
    }

    return f;
}

double log_factorial(double n) {
    double fact = 0;

    for (double i = n; i > 1; i--)
        fact += log(i);

    return fact;
}

double reduced_log_factorial(double n, double d) {
    double f = n - d;

    if (f < 0) 
        return 0;
    else if (f == 0) 
        return 0;
    else if (f == 1)
        return log(n);

    f = log(n);
    n--;
    while (n > d) {
        f += log(n);
        n--;
    }

    return f;
}

bool compare_pair(pair<char, int> a, pair<char, int> b) {
    return (a.second > b.second);
}

bool compare_pair_intdouble(pair<int, double> a, pair<int, double> b) {
    return (a.second < b.second);
}

bool compare_ints(int a, int b) {
    return (a > b);
}

bool compare_pair_snp(pair<string, SNP *> a, pair<string, SNP *> b) {
    return (a.second->col < b.second->col);
}
