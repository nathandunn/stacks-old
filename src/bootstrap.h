// -*-mode:c++; c-style:k&r; c-basic-offset:4;-*-
//
// Copyright 2014, Julian Catchen <jcatchen@uoregon.edu>
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

#ifndef __BOOTSTRAP_H__
#define __BOOTSTRAP_H__

#include <math.h>

extern double   sigma;
extern int      bootstrap_reps;
extern set<int> bootstraplist;

template<class StatT=PopStat>
class Bootstrap {
    double *weights; // Weight matrix to apply while smoothing.

    int calc_weights();

public:
    Bootstrap(int size)  { 
	this->calc_weights();

    }
    ~Bootstrap() { 
	delete [] this->weights;
    }

    int initialize();
    int execute();
};

template<class StatT>
int 
Bootstrap<StatT>::calc_weights() 
{
    int limit = 3 * sigma;
    //
    // Calculate weights for window smoothing operations.
    //
    // For each genomic region centered on a nucleotide position c, the contribution of the population 
    // genetic statistic at position p to the region average was weighted by the Gaussian function:
    //   exp( (-1 * (p - c)^2) / (2 * sigma^2))
    //
    this->weights = new double[limit + 1];

    for (int i = 0; i <= limit; i++)
	this->weights[i] = exp((-1 * pow(i, 2)) / (2 * pow(sigma, 2)));

    return 0;
}

template<class StatT>
int
Bootstrap<StatT>::initialize()
{

    return 0;
}

template<class StatT>
int
Bootstrap<StatT>::execute()
{

    return 0;
}

#endif // __BOOTSTRAP_H__
