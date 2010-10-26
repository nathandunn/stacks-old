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
// cmb.cc -- routines to implement the Combination generating class: CombSet.
//
// Julian Catchen
// jcatchen@uoregon.edu
// University of Oregon
//
// $Id$
//
#include "cmb.h"

CombSet::~CombSet() {
    int num_sets = this->sets.size();
    int set, i;

    for (set = 0; set < num_sets; set++) {
        for (i = 0; i < this->lens[set]; i++)
            delete [] sets[set][i];
        delete [] sets[set];
    }

    num_sets = this->compound_comb.size();

    for (set = 0; set < num_sets; set++)
        delete [] this->compound_comb[set]->elem;
    delete this->compound_comb[set];
}

CombSet::CombSet(int n, int k) {
    this->max_set_size = n >= k ? k : n;
    this->num_elements = n;
    this->index        = 0;

    int   i;
    int   set_size = this->max_set_size;
    int   size;
    int **comb;
    Cmb  *new_comb;

    while (set_size > 0) {
        //
        // How many combinations will we make?
        //
        size = (int) num_combinations(this->num_elements, set_size);

        //
        // Generate all combinations, N choose K; N=num_elements, K=set_size
        //
        comb = this->generate_combinations(this->num_elements, set_size, size, false);

        this->lens.push_back(size);
        this->size.push_back(set_size);
        this->sets.push_back(comb);

        set_size--;
    }

    //
    // Create a single list of all the combinations generated above.
    //
    this->make_compound_set();

    //
    // Finally, generate all combinations of the compound set.
    //
    for (set_size = 1; set_size < this->num_elements; set_size++) {

        size = (int) num_combinations(this->compound_set.size(), set_size);

        cerr << "Num elements: " << this->compound_set.size() << "; subsets of size K: " << set_size << "; total: " << size;
        comb = this->generate_combinations(this->compound_set.size(), set_size, size, true);
        cerr << "; Valid sets: " << size << "\n";

        for (i = 0; i < size; i++) {
            new_comb = new Cmb;
            new_comb->size = set_size;
            new_comb->elem = comb[i];

            this->compound_comb.push_back(new_comb);
        }

        delete [] comb;
    }

    //
    // We can add the final combination by hand, when set_size == this->num_elements, since it generates
    // the most combinations and always results in all but 1 combination being removed.
    //
    new_comb = new Cmb;
    new_comb->size = this->num_elements;
    new_comb->elem = new int[this->num_elements];
    int offset = 0;
    for (i = 0; i < this->max_set_size - 1; i++)
        offset += this->lens[i];
    for (i = 0; i < this->num_elements; i++)
        new_comb->elem[i] = offset + i;

    this->compound_comb.push_back(new_comb);

    cerr << "Total compound combinations: " << this->compound_comb.size() << "\n";
}

int CombSet::make_compound_set() {

    int num_sets = this->lens.size();

    for (int i = 0; i < num_sets; i++)
        for (int j = 0; j < this->lens[i]; j++)
            this->compound_set.push_back(make_pair(i, j));

    return 0;
}

bool CombSet::valid(int *comb, int comb_size) {
    int count;
    int index;
    int k, n, i, j;

    count = 0;

    for (i = 0; i < comb_size; i++) {
        index = comb[i];

        // sets vector index number
        k = this->compound_set[index].first;

        count += this->size[k];
    }

    //
    // Combination not valid
    //
    if (count != this->num_elements)
        return false;

    //
    // Check that each of the original elements exists only once in all subsets.
    //
    pair<set<int>::iterator,bool> ret;
    set<int> s;

    for (i = 0; i < comb_size; i++) {
        index = comb[i];

        // sets vector index number
        k = this->compound_set[index].first;
        // combination number 
        n = this->compound_set[index].second;

        //
        // Items must be present only once in the final combination.
        //
        for (j = 0; j < this->size[k]; j++) {
            ret = s.insert(this->sets[k][n][j]);
            if (ret.second == false)
                return false;
        }
    }

    return true;
}

int CombSet::count_valid_comb(int total, int n, int k) {
    int *curr = new int[k];

    //
    // Setup the initial combination
    //
    int comb_num = 0;

    for (int i = 0; i < k; i++)
        curr[i] = i;

    if (this->valid(curr, k))
        comb_num++;

    //
    // Generate each successive combination
    //
    int j = 1;
    while (j < total) {
        next_combination(curr, n, k);
        if (this->valid(curr, k))
            comb_num++;
        j++;
    }
    delete [] curr;

    return comb_num;
}

int **CombSet::generate_combinations(int n, int k, int &total, bool validate) {
    int **comb;
    int  *curr;

    if (validate)
        total = this->count_valid_comb(total, n, k);

    comb = new int * [total];
    for (int i = 0; i < total; i++)
        comb[i] = new int[k];
    curr = new int[k];

    //
    // Setup the initial combination
    //
    int comb_num = 0;

    for (int i = 0; i < k; i++)
        curr[i] = i;

    if (validate == false || this->valid(curr, k)) {
        for (int i = 0; i < k; i++)
            comb[comb_num][i] = curr[i];
        comb_num++;
    }

    //
    // Generate each successive combination
    //
    while (comb_num < total) {
        next_combination(curr, n, k);

        if (validate == false || this->valid(curr, k)) {
            for (int i = 0; i < k; i++)
                comb[comb_num][i] = curr[i];
            comb_num++;
        }
    }

    delete [] curr;

    return comb;
}

int CombSet::next_combination(int *comb, int n, int k) {
    int i;

    // 
    // The zero'th position has been incremented to its maximal value,
    // it's not possible to further increment values in the set.
    //
    if (comb[0] > n - k)
        return 0;

    //
    // Increment the last position in the set.
    //
    i = k - 1;
    comb[i]++;

    //
    // Check if the last position has reached its maximal possible value, 
    // if so, move back one position, and increment it.
    //
    while ((i > 0) && (comb[i] >= n - k + 1 + i)) {
        i--;
        comb[i]++;
    }

    //
    // Move from the position we incremented above back out to the final position
    //
    for (i = i + 1; i < k; i++)
        comb[i] = comb[i - 1] + 1;

    return 1;
}

double CombSet::num_combinations(int n, int k) {
    double r = 1;

    for (int i = n; i >= (n - k + 1); i--)
        r *= i;
    double s = factorial(k);

    float num_comb = r / s;

    return lroundf(num_comb);
}

//
// Return a variable length array of Cmb objects, terminated by a NULL pointer.
//
Cmb **CombSet::next(int map[]) {

    if (this->index >= (int) this->compound_comb.size())
        return NULL;

    int  index, i, j, k, n;
    int  size = this->compound_comb[this->index]->size;
    int *e    = this->compound_comb[this->index]->elem;

    Cmb **c = new Cmb * [size + 1];

    for (i = 0; i < size; i++) {
        index = e[i];
        // sets vector index number
        k = this->compound_set[index].first;
        // combination number 
        n = this->compound_set[index].second;

        c[i] = new Cmb;
        c[i]->size = this->size[k];
        c[i]->elem = new int[this->size[k]];

        for (j = 0; j < this->size[k]; j++)
            c[i]->elem[j] = (map == NULL) ? 
                this->sets[k][n][j] : 
                map[this->sets[k][n][j]];
    }

    c[size] = NULL;

    this->index++;

    return c;
}

void CombSet::reset() {
    this->index = 0;
}

void CombSet::destroy(Cmb **cmb) {

    for (uint j = 0; cmb[j] != NULL; j++) {
        delete [] cmb[j]->elem;
        delete cmb[j];
    }
    delete [] cmb;
}

void write_cmb(int *comb, int size) {
    stringstream s;
    string t;

    s << "{";

    for (int i = 0; i < size; i++)
        s << comb[i] << ", ";
    t = s.str().substr(0, s.str().length() - 2);
    t += "}";

    cerr << t << "\n";
}
