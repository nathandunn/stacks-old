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

//
// utils.cc -- common routines needed in multiple object files.
//
// Julian Catchen
// jcatchen@uoregon.edu
// University of Oregon
//
// $Id$
//
#include "utils.h"

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

bool compare_pair(pair<char, int> a, pair<char, int> b) {
    return (a.second > b.second);
}

