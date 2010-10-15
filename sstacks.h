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

// -*-mode:c++; c-style:k&r; c-basic-offset:4;-*-

#ifndef __SSTACKS_H__
#define __SSTACKS_H__

#include <omp.h>    // OpenMP library
#include <getopt.h> // Process command-line options
#include <math.h>
#include <stdlib.h>
#include <utility>
using std::pair;
using std::make_pair;

#include <string>
using std::string;

#include <iostream>
#include <fstream>
using std::ifstream;
using std::ofstream;
using std::cin;
using std::cout;
using std::cerr;
using std::endl;

#include <vector>
using std::vector;

#include <map>
using std::map;

#include <set>
using std::set;

#ifdef __GNUC__
#include <ext/hash_map>
#include <ext/hash_fun.h>
using __gnu_cxx::hash_map;
using __gnu_cxx::hash;
#else
#include <hash_map>
#endif

#include "constants.h"
#include "stacks.h"
#include "sql_utilities.h"

//
// Query Locus Class
//
class QLocus : public Locus {
 public:
    vector<pair<int, allele_type> > matches;   // Matching tags found for the catalog. Stored as catalog ID/allele_type pair.

    int add_match(int, allele_type);
};

int QLocus::add_match(int id, allele_type type) {
    this->matches.push_back(make_pair(id, type));

    return 0;
}

struct eqstr {
    bool operator()(const char* s1, const char* s2) const {
	return strcmp(s1, s2) == 0;
    }
};

typedef hash_map<const char *, vector<pair<int, allele_type> >, hash<const char *>, eqstr> HashMap;

void help( void );
void version( void );
int  parse_command_line(int, char**);
int  populate_hash(map<int, Locus *> &, HashMap &);
int  find_matches(map<int, Locus *> &, map<int, QLocus *> &);
int  verify_match(Locus *, QLocus *, uint);
bool compare_dist(pair<int, int>, pair<int, int>);
int  write_matches(map<int, QLocus *> &);
bool compare_pair(pair<string, SNP *>, pair<string, SNP *>);

#endif // __SSTACKS_H__
