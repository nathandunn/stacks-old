// -*-mode:c++; c-style:k&r; c-basic-offset:4;-*-
//
// Copyright 2016, Julian Catchen <jcatchen@illinois.edu>
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

#ifndef __ALN_UTILS_H__
#define __ALN_UTILS_H__

#include <utility>
#include <string>
#include <vector>
#include <tuple>

#include "constants.h"
#include "utils.h"

string invert_cigar(string);
int    parse_cigar(const char *, vector<pair<char, uint> > &, bool check_correctness = false);
string apply_cigar_to_seq(const char *, vector<pair<char, uint> > &);
string remove_cigar_from_seq(const char *, vector<pair<char, uint> > &);
string apply_cigar_to_model_seq(const char *, vector<pair<char, uint> > &);
int    apply_cigar_to_seq(char *, uint, const char *, vector<pair<char, uint> > &);
int    apply_cigar_to_model_seq(char *, uint, const char *, vector<pair<char, uint> > &);
std::tuple<uint,uint,uint> cigar_lengths(const vector<pair<char, uint>>&);
void cigar_extend_right(vector<pair<char, uint>>&, size_t);
void cigar_extend_left(vector<pair<char, uint>>&, size_t);

#include "locus.h"
int    adjust_snps_for_gaps(vector<pair<char, uint> > &, Locus *);
int    adjust_and_add_snps_for_gaps(vector<pair<char, uint> > &, Locus *);
int    remove_snps_from_gaps(vector<pair<char, uint> > &, Locus *);

//
// Inline definitions.
// ==========
//

inline
void cigar_extend_right(vector<pair<char, uint>>& cig, size_t len) {
    if (cig.back().first == 'D')
        cig.back().second += len;
    else
        cig.push_back({'D', len});
}

inline
void cigar_extend_left(vector<pair<char, uint>>& cig, size_t len) {
    if (cig.front().first == 'D')
        cig.front().second += len;
    else
        cig.insert(cig.begin(), {'D', len});
}

#endif  // __ALN_UTILS_H__
