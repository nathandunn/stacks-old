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
        cerr << "Num elements: " << this->num_elements << "; subsets of size K: " << set_size << "\n";
        comb = this->generate_combinations(this->num_elements, set_size, size);

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

        cerr << "Num elements: " << this->compound_set.size() << "; subsets of size K: " << set_size << "; total: " << size << "\n";
        comb = this->generate_combinations(this->compound_set.size(), set_size, size);

        //
        // Remove combinations that repeat elements
        //
        this->remove_duplicates(comb, size, set_size);

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

int CombSet::remove_duplicates(int **&comb, int &size, int comb_size) {
    int **clean      = NULL;
    int   clean_size = 0;

    pair<set<int>::iterator,bool> ret;
    set<int> s;
    int count;
    int index;
    int k, n, i, j, r;

    //
    // Initialize key to keep track of which permutations we keep/discard
    //
    int  *key = new int[size];
    for (i = 0; i < size; i++)
        key[i] = 1;

    for (i = 0; i < size; i++) {
        count = 0;

        //cerr << "Looking at combination " << i << "\n";

        for (j = 0; j < comb_size; j++) {
            index = comb[i][j];

            // sets vector index number
            k = this->compound_set[index].first;
            // combination number 
            n = this->compound_set[index].second;

            count += this->size[k];
        }

        if (count != this->num_elements) {
            key[i] = 0; // Discard this combination
            continue;
        }

        //
        // Check that each of the original elements exists only once in all subsets.
        //
        for (j = 0; j < comb_size; j++) {
            index = comb[i][j];

            // sets vector index number
            k = this->compound_set[index].first;
            // combination number 
            n = this->compound_set[index].second;

            //cerr << "    Set Num: " << k << "; Offset: " << n << "; Size: " << this->size[k] << "\n";

            //
            // Items must be present only once in the final combination.
            //
            r = 0;
            while (key[i] == 1 && r < this->size[k]) {
                ret = s.insert(this->sets[k][n][r]);
                if (ret.second == false) {
                    key[i] = 0; // Discard this combination
                }
                r++;
            }

        }

        s.clear();

        if (key[i] == 1) { 
            clean_size++;
        }
    }

    if (clean_size > 0)
        clean = new int * [clean_size];

    j = 0;
    for (i = 0; i < size; i++) {
        if (key[i] == 1) {
            clean[j] = comb[i];
            j++;
        }
        else 
            delete [] comb[i];
    }
     
    delete [] comb;
    delete [] key;

    comb = clean;
    size = clean_size;

    return 0;
}

int **CombSet::generate_combinations(int n, int k, int total) {
    int **comb;

    cerr << "Number of combinations given a set of size " << n << " choosing subsets of " << k << ": " << total << "\n";

    comb = new int * [total];
    for (int i = 0; i < total; i++)
        comb[i] = new int[k];

    //
    // Setup the initial combination
    //
    int comb_num = 0;

    for (int i = 0; i < k; i++)
        comb[comb_num][i] = i;
    comb_num++;

    //
    // Generate and print all the other combinations
    //
    while (comb_num < total) {
        next_combination(comb[comb_num - 1], comb[comb_num], n, k);
        comb_num++;
    }

    return comb;
}

int CombSet::next_combination(int *prev_comb, int *comb, int n, int k) {
    int i;

    // 
    // The zero'th position has been incremented to its maximal value,
    // it's not possible to further increment values in the set.
    //
    if (prev_comb[0] > n - k)
        return 0;

    //
    // Copy the previous combination before incrementing it to the next combination
    //
    for (i = 0; i < k; i++)
        comb[i] = prev_comb[i];

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
Cmb **CombSet::next() {

    if (this->index >= (int) this->compound_comb.size())
        return NULL;

    int  index, i, k, n;
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
        c[i]->elem = this->sets[k][n];
    }

    c[size] = NULL;

    this->index++;

    return c;
}

void CombSet::reset() {
    this->index = 0;
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

double factorial(int n) {
    double fact = 1;

    for (int i = n; i > 0; i--)
        fact *= i;

    return fact;
}
