// -*-mode:c++; c-style:k&r; c-basic-offset:4;-*-
//
// Copyright 2011, Julian Catchen <jcatchen@uoregon.edu>
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

#ifndef __RENZ_H__
#define __RENZ_H__

#include <map>
using std::map;
#include <string>
using std::string;

const char *sbfI[]    = {"TGCAGG",            // CCTGCA/GG, SbfI
			 "CCTGCA"};
const char *pstI[]    = {"TGCAG",             // CTGCA/G, PstI
			 "CTGCA"};
const char *ecoRI[]   = {"AATTC",             // G/AATTC, EcoRI
			 "GAATT"};
const char *sgrAI[]   = {"CCGGCG", "CCGGTG",  // CR/CCGGYG, SgrAI; R=A or G; Y=C or T
			 "CGCCGG", "CACCGG"};
const char *notI[]    = {"GGCCGC",            // GC/GGCCGC, NotI
			 "GCGGCC"};
const char *haeIII[]  = {"GGCC",
			 "GGCC"};
const char *aluI[]    = {"AGCT",
			 "AGCT"};
const char *apeKI[]   = {"CAGC", "CTGC",    // G/CWGC, ApeKI; W=A or T
 			 "GTCG", "GACG"};
const char *mseI[]    = {"TTAA",
			 "TTAA"};
const char *hindIII[] = {"AAGCT",
			 "TTCGA"};
const char *dpnII[]   = {"GATC",
			 "GATC"};

void 
initialize_renz(map<string, const char **> &renz, map<string, int> &renz_cnt, map<string, int> &renz_len) {

    renz["sbfI"]    = sbfI;    // CCTGCA/GG, SbfI
    renz["pstI"]    = pstI;    // CTGCA/G, PstI
    renz["notI"]    = notI;    // GC/GGCCGC, NotI
    renz["ecoRI"]   = ecoRI;   // G/AATTC, EcoRI
    renz["sgrAI"]   = sgrAI;   // CR/CCGGYG, SgrAI; R=A or G; Y=C or T
    renz["apeKI"]   = apeKI;   // G/CWGC, ApeKI; W=A or T
    renz["hindIII"] = hindIII; // A/AGCTT, HindIII
    renz["dpnII"]   = dpnII;   // GATC, DpnII

    renz_cnt["sbfI"]    = 1;
    renz_cnt["pstI"]    = 1;
    renz_cnt["notI"]    = 1;
    renz_cnt["ecoRI"]   = 1;
    renz_cnt["sgrAI"]   = 2;
    renz_cnt["apeKI"]   = 2;
    renz_cnt["hindIII"] = 1;
    renz_cnt["dpnII"]   = 1;

    renz_len["sbfI"]    = 6;
    renz_len["pstI"]    = 5;
    renz_len["notI"]    = 6;
    renz_len["ecoRI"]   = 5;
    renz_len["sgrAI"]   = 6;
    renz_len["apeKI"]   = 4;
    renz_len["hindIII"] = 5;
    renz_len["dpnII"]   = 4;
}

#endif // __RENZ_H__
