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

#ifndef __SMOOTHING_H__
#define __SMOOTHING_H__

#include <math.h>
#include "smoothing_utils.h"

template<class StatT=PopStat>
class KSmooth {
    uint    size;    // Number of elements expected in the StatT class to smooth.
    double *weights; // Weight matrix to apply while smoothing.

public:
    KSmooth(int size)  { 
	this->size    = size;
	this->weights = calc_weights();

    }
    ~KSmooth() { 
	delete [] this->weights;
    }

    int smooth(vector<StatT *> &popstats);
};

template<class StatT>
int
KSmooth<StatT>::smooth(vector<StatT *> &popstats)
{
    //
    // To generate smooth genome-wide distributions of Fst, we calculate a kernel-smoothing 
    // moving average of Fst values along each ordered chromosome.
    //
    // For each genomic region centered on a nucleotide position c, the contribution of the population 
    // genetic statistic at position p to the region average was weighted by the Gaussian function:
    //   exp( (-1 * (p - c)^2) / (2 * sigma^2))
    // 
    // In addition, we weight each position according to (n_k - 1), where n_k is the number of alleles
    // sampled at that location.
    //
    // By default, sigma = 150Kb, for computational efficiency, only calculate average out to 3sigma.
    //
    #pragma omp parallel
    { 
	int      dist;
	uint     pos_l, pos_u;
	double   sum, final_weight;
	PopStat *c, *p;

	pos_l = 0;
	pos_u = 0;

        #pragma omp for schedule(dynamic, 1)
	for (uint pos_c = 0; pos_c < popstats.size(); pos_c++) {
	    c = popstats[pos_c];

	    if (c == NULL)
		continue;

	    for (uint i = 0; i < this->size; i++)
		c->smoothed[i] = 0.0;
	    sum = 0.0;

	    determine_window_limits(popstats, c->bp, pos_l, pos_u);

	    for (uint pos_p = pos_l; pos_p < pos_u; pos_p++) {
		p = popstats[pos_p];

		if (p == NULL)
		    continue;

		dist = p->bp > c->bp ? p->bp - c->bp : c->bp - p->bp;

		if (dist > limit || dist < 0) {
		    #pragma omp critical
		    {
			cerr << "ERROR: current basepair is out of the sliding window.\n"
			     << "  Calculating sliding window; start position: " << pos_l << ", " << (popstats[pos_l] == NULL ? -1 : popstats[pos_l]->bp) << "bp; end position: " 
			     << pos_u << ", " << (popstats[pos_u] == NULL ? -1 : popstats[pos_u]->bp) << "bp; center: " 
			     << pos_c << ", " << popstats[pos_c]->bp << "bp\n"
			     << "  Current position: " << pos_p << ", " << popstats[pos_p]->bp << "; Dist: " << dist << "\n"
			     << "  Window positions:\n";

			for (uint j = pos_l; j < pos_u; j++) {
			    p = popstats[j];
			    if (p == NULL) continue;
			    cerr << "    Position: " << j << "; " << p->bp << "bp\n";
			}
			//exit(0);
		    }
		    continue;
		}

		// sites_cnt++;

		final_weight = (p->alleles - 1) * this->weights[dist];
		for (uint i = 0; i < this->size; i++)
		    c->smoothed[i] += p->stat[i] * final_weight;
		sum += final_weight;

		// if (c->loc_id == 9314) {
		//     cerr << "    id: " << p->loc_id 
		// 	 << "; dist: " << dist 
		// 	 << "; weight: " << weights[dist] 
		// 	 << "; final_weight: " << final_weight 
		// 	 << "; fst': " << p->stat[3] 
		// 	 << "; sum: " << sum 
		// 	 << "; smoothed: " << c->smoothed[3] << "\n";
		// }
	    }

	    // sites_per_snp += (sites_cnt / snp_cnt);
	    // tot_windows++;
	    //
	    // if (snp_cnt < max_snp_dist) {
	    //     #pragma omp atomic
	    //     snp_dist[snp_cnt]++;
	    // }
	    // c->snp_cnt = snp_cnt;

	    for (uint i = 0; i < this->size; i++)
		c->smoothed[i] /= sum;
	}
    }

    return 0;
}

#endif // __SMOOTHING_H__
