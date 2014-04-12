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
#include <vector>
using std::vector;

extern double   sigma;
extern int      bootstrap_reps;
extern bool     bootstrap_wl;
extern set<int> bootstraplist;

//
// Bootstrap resamplign structure.
//
class BSample {
public:
    int    bp;
    int    alleles;
    double stat[PopStatSize];

    BSample() {
	this->bp      = 0;
	this->alleles = 0;
	for (uint i = 0; i < PopStatSize; i++)
	    this->stat[i] = 0.0;

    }
};

template<class StatT=PopStat>
class Bootstrap {
    double                 *weights; // Weight matrix to apply while smoothing.
    vector<vector<double> > stats;
    int                     num_stats;

    int calc_weights();

public:
    Bootstrap(int size)  { 
	this->num_stats = size;
	this->calc_weights();

    }
    ~Bootstrap() { 
	delete [] this->weights;
    }

    int    add_data(vector<StatT *> &);
    int    execute(vector<StatT *> &);
    double pval(double, vector<double> &);
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
Bootstrap<StatT>::add_data(vector<StatT *> &sites)
{
    for (uint i = 0; i < sites.size(); i++) {
	if (sites[i] != NULL)
	    for (uint j = 0; j < this->num_stats; j++)
		this->stats[j].push_back(sites[i]->stat[j]);
    }

    return 0;
}

template<class StatT>
int
Bootstrap<StatT>::execute(vector<StatT *> &sites)
{
    #pragma omp parallel
    { 
	PopStat *c;
	double final_weight, sum, weighted_stat[PopStatSize];
	int  dist, index, limit_l, limit_u;
	int  limit = 3 * sigma;
	uint pos_l = 0;
	uint pos_u = 0;

        #pragma omp for schedule(dynamic, 1)  
	for (uint pos_c = 0; pos_c < sites.size(); pos_c++) {
	    c = sites[pos_c];

	    if (c == NULL)
		continue;

	    if (bootstrap_wl && bootstraplist.count(c->loc_id) == 0)
		continue;

	    limit_l = c->bp - limit > 0 ? c->bp - limit : 0;
	    limit_u = c->bp + limit;

	    while (pos_l < sites.size()) {
		if (sites[pos_l] == NULL) {
		    pos_l++;
		} else {
		    if (sites[pos_l]->bp < limit_l) 
			pos_l++;
		    else
			break;
		}
	    }
	    while (pos_u < sites.size()) {
		if (sites[pos_u] == NULL) {
		    pos_u++;
		} else {
		    if (sites[pos_u]->bp < limit_u)
			pos_u++;
		    else
			break;
		}
	    }
	    if (pos_u < sites.size() && sites[pos_u]->bp > limit_u) pos_u--;

	    int size = 0;
	    for (uint i = pos_l; i < pos_u;  i++) {
		if (sites[i] == NULL) continue;
		size++;
	    }

	    //
	    // Allocate an array of bootstrap resampling objects.
	    //
	    BSample *bs = new BSample[size];

	    //
	    // Populate the BSample objects.
	    //
	    int j = 0;
	    for (uint i = pos_l; i < pos_u;  i++) {
		if (sites[i] == NULL) continue;

		bs[j].bp      = sites[i]->bp;
		bs[j].alleles = sites[i]->alleles;
		j++;
	    }

	    vector<vector<double> > resampled_stats;
	    for (uint i = 0; i < this->num_stats; i++)
		resampled_stats[i].reserve(bootstrap_reps);

	    // cerr << "Window starts at " << bs[0].bp << "; centered on " << c->bp << "\n";

	    //
	    // Bootstrap this bitch.
	    //
	    for (int i = 0; i < bootstrap_reps; i++) {
		// cerr << "  Bootsrap rep " << i << "\n";

		for (uint k = 0; k < this->num_stats; k++)
		    weighted_stat[k] = 0.0;
		sum = 0.0;

		for (j = 0; j < size; j++) {

		    dist = bs[j].bp > c->bp ? bs[j].bp - c->bp : c->bp - bs[j].bp;

		    //
		    // Resample for this round of bootstrapping.
		    //
		    for (uint k = 0; k < this->num_stats; k++) {
			index         = (int) (this->stats[k].size() * (random() / (RAND_MAX + 1.0)));
			bs[j].stat[k] = this->stats[k][index];
			// cerr << "      WinPos: " << j << "; Randomly selecting " << index << " out of " << fst_samples.size() << " possible values giving Fst value: " << bs[j].f << "\n";
		    }

		    final_weight = (bs[j].alleles - 1) * weights[dist];
		    for (uint k = 0; k < this->num_stats; k++)
			weighted_stat[k] += bs[j].stat[k] * final_weight;
		    sum += final_weight;
		}

		// cerr << "    New weighted Fst value: " << weighted_fst / sum << "\n";
		for (uint k = 0; k < this->num_stats; k++)
		    resampled_stats[k].push_back(weighted_stat[k] / sum);
	    }

	    //
	    // Cacluate the p-value for this window based on the empirical Fst distribution.
	    //
	    for (uint k = 0; k < this->num_stats; k++) {
		sort(resampled_stats[k].begin(), resampled_stats[k].end());
		c->bs[k] = this->pval(c->smoothed[k], resampled_stats[k]);
	    }

	    delete [] bs;
	}
    }

    return 0;
}

template<class StatT>
double
Bootstrap<StatT>::pval(double stat, vector<double> &dist)
{
    vector<double>::iterator up;
    double pos;

    up = upper_bound(dist.begin(), dist.end(), stat);

    if (up == dist.begin())
	pos = 1;
    else if (up == dist.end())
	pos = dist.size();
    else 
	pos = up - dist.begin() + 1;

    double res = 1.0 - (pos / (double) dist.size());

    // cerr << "Generated Smoothed Fst Distribution:\n";
    // for (uint n = 0; n < dist.size(); n++)
    // 	cerr << "  n: " << n << "; Fst: " << dist[n] << "\n";

    // cerr << "Comparing Fst value: " << stat 
    // 	 << " at position " << (up - dist.begin()) << " out of " 
    // 	 << dist.size() << " positions (converted position: " << pos << "); pvalue: " << res << ".\n";

    return res;
}


#endif // __BOOTSTRAP_H__
