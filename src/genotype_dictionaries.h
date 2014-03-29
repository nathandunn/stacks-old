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

#ifndef __GENOTYPE_DICTIONARIES_H__
#define __GENOTYPE_DICTIONARIES_H__

void 
initialize_dictionaries(map<string, map<string, string> > &global_dictionary) 
{
    global_dictionary["ab/--"]["a"]  = "aa";
    global_dictionary["ab/--"]["b"]  = "bb";

    global_dictionary["--/ab"]["a"]  = "aa";
    global_dictionary["--/ab"]["b"]  = "bb";

    global_dictionary["aa/bb"]["a"]  = "aa";
    global_dictionary["aa/bb"]["ab"] = "ab";
    global_dictionary["aa/bb"]["b"]  = "bb";

    global_dictionary["ab/ac"]["a"]  = "aa";
    global_dictionary["ab/ac"]["ab"] = "ab";
    global_dictionary["ab/ac"]["b"]  = "bb";
    global_dictionary["ab/ac"]["ac"] = "ac";
    global_dictionary["ab/ac"]["c"]  = "cc";
    global_dictionary["ab/ac"]["bc"] = "bc";

    global_dictionary["ab/cd"]["a"]  = "aa";
    global_dictionary["ab/cd"]["ab"] = "ab";
    global_dictionary["ab/cd"]["b"]  = "bb";
    global_dictionary["ab/cd"]["c"]  = "cc";
    global_dictionary["ab/cd"]["cd"] = "cd";
    global_dictionary["ab/cd"]["d"]  = "dd";
    global_dictionary["ab/cd"]["ac"] = "ac";
    global_dictionary["ab/cd"]["ad"] = "ad";
    global_dictionary["ab/cd"]["bc"] = "bc";
    global_dictionary["ab/cd"]["bd"] = "bd";

    global_dictionary["ab/aa"]["a"]  = "aa";
    global_dictionary["ab/aa"]["ab"] = "ab";
    global_dictionary["ab/aa"]["b"]  = "bb";

    global_dictionary["aa/ab"]["a"]  = "aa";
    global_dictionary["aa/ab"]["ab"] = "ab";
    global_dictionary["aa/ab"]["b"]  = "bb";

    global_dictionary["ab/cc"]["a"]  = "aa";
    global_dictionary["ab/cc"]["ab"] = "ab";
    global_dictionary["ab/cc"]["b"]  = "bb";
    global_dictionary["ab/cc"]["c"]  = "cc";
    global_dictionary["ab/cc"]["ac"] = "ac";
    global_dictionary["ab/cc"]["bc"] = "bc";

    global_dictionary["cc/ab"]["a"]  = "aa";
    global_dictionary["cc/ab"]["ab"] = "ab";
    global_dictionary["cc/ab"]["b"]  = "bb";
    global_dictionary["cc/ab"]["c"]  = "cc";
    global_dictionary["cc/ab"]["ac"] = "ac";
    global_dictionary["cc/ab"]["bc"] = "bc";

    global_dictionary["ab/ab"]["a"]  = "aa";
    global_dictionary["ab/ab"]["b"]  = "bb";
    global_dictionary["ab/ab"]["ab"] = "ab";
}

void 
load_joinmap_cp_dictionary(map<string, string> &types, 
			   map<string, map<string, string> > &dictionary,
			   map<string, map<string, double> > &segregation_ratios)
{
    types["ab/--"] = "lmx--";
    types["--/ab"] = "--xnp";
    types["ab/aa"] = "lmxll";
    types["aa/ab"] = "nnxnp";
    types["ab/a-"] = "lmx--";
    types["-a/ab"] = "--xnp";
    types["ab/c-"] = "lmx--";
    types["-c/ab"] = "--xnp";
    types["ab/cc"] = "lmx--";
    types["cc/ab"] = "--xnp";
    types["ab/ab"] = "hkxhk";
    types["ab/ac"] = "efxeg";
    types["ab/cd"] = "abxcd";

    dictionary["lmx--"]["-"]  = "--";
    dictionary["lmx--"]["aa"] = "ll";
    dictionary["lmx--"]["ab"] = "lm";
    dictionary["lmx--"]["ac"] = "ll";
    dictionary["lmx--"]["bc"] = "lm";

    dictionary["--xnp"]["-"]  = "--";
    dictionary["--xnp"]["aa"] = "nn";
    dictionary["--xnp"]["ab"] = "np";
    dictionary["--xnp"]["ac"] = "nn";
    dictionary["--xnp"]["bc"] = "np";

    dictionary["lmxll"]["-"]  = "--";
    dictionary["lmxll"]["aa"] = "ll";
    dictionary["lmxll"]["ab"] = "lm";

    dictionary["nnxnp"]["-"]  = "--";
    dictionary["nnxnp"]["aa"] = "nn";
    dictionary["nnxnp"]["ab"] = "np";

    dictionary["hkxhk"]["-"]  = "--";
    dictionary["hkxhk"]["ab"] = "hk";
    dictionary["hkxhk"]["aa"] = "hh";
    dictionary["hkxhk"]["bb"] = "kk";

    dictionary["efxeg"]["-"]  = "--";
    dictionary["efxeg"]["ab"] = "ef";
    dictionary["efxeg"]["ac"] = "eg";
    dictionary["efxeg"]["bc"] = "fg";
    dictionary["efxeg"]["aa"] = "ee";

    dictionary["abxcd"]["-"]  = "--";
    dictionary["abxcd"]["ac"] = "ac";
    dictionary["abxcd"]["ad"] = "ad";
    dictionary["abxcd"]["bc"] = "bc";
    dictionary["abxcd"]["bd"] = "bd";

    segregation_ratios["lmx--"]["ll"] = 0.50;
    segregation_ratios["lmx--"]["lm"] = 0.50;

    segregation_ratios["--xnp"]["nn"] = 0.50;
    segregation_ratios["--xnp"]["np"] = 0.50;

    segregation_ratios["lmxll"]["ll"] = 0.50;
    segregation_ratios["lmxll"]["lm"] = 0.50;

    segregation_ratios["nnxnp"]["nn"] = 0.50;
    segregation_ratios["nnxnp"]["np"] = 0.50;

    segregation_ratios["hkxhk"]["hk"] = 0.50;
    segregation_ratios["hkxhk"]["hh"] = 0.25;
    segregation_ratios["hkxhk"]["kk"] = 0.25;

    segregation_ratios["efxeg"]["ef"] = 0.25;
    segregation_ratios["efxeg"]["eg"] = 0.25;
    segregation_ratios["efxeg"]["fg"] = 0.25;
    segregation_ratios["efxeg"]["ee"] = 0.25;

    segregation_ratios["abxcd"]["ac"] = 0.25;
    segregation_ratios["abxcd"]["ad"] = 0.25;
    segregation_ratios["abxcd"]["bc"] = 0.25;
    segregation_ratios["abxcd"]["bd"] = 0.25;

    return;
}

void 
load_onemap_cp_dictionary(map<string, string> &types, map<string, map<string, string> > &dictionary) 
{
    types["ab/--"] = "abxoo";
    types["--/ab"] = "ooxab";
    types["ab/aa"] = "abxaa";
    types["aa/ab"] = "aaxab";
    types["ab/ab"] = "abxab";
    types["ab/ac"] = "abxac";
    types["ab/cd"] = "abxcd";

    // D1.11
    dictionary["abxoo"]["-"]  = "-";
    dictionary["abxoo"]["aa"] = "a";
    dictionary["abxoo"]["bb"] = "b";

    // D2.16
    dictionary["ooxab"]["-"]  = "-";
    dictionary["ooxab"]["aa"] = "a";
    dictionary["ooxab"]["bb"] = "b";

    // D1.10
    dictionary["abxaa"]["-"]  = "-";
    dictionary["abxaa"]["aa"] = "a";
    dictionary["abxaa"]["ab"] = "ab";

    // D2.15
    dictionary["aaxab"]["-"]  = "-";
    dictionary["aaxab"]["aa"] = "a";
    dictionary["aaxab"]["ab"] = "ab";

    // B3.7
    dictionary["abxab"]["-"]  = "-";
    dictionary["abxab"]["ab"] = "ab";
    dictionary["abxab"]["aa"] = "a";
    dictionary["abxab"]["bb"] = "b";

    // A.2
    dictionary["abxac"]["-"]  = "-";
    dictionary["abxac"]["ab"] = "ba";
    dictionary["abxac"]["ac"] = "ac";
    dictionary["abxac"]["bc"] = "bc";
    dictionary["abxac"]["aa"] = "a";

    // A.1
    dictionary["abxcd"]["-"]  = "-";
    dictionary["abxcd"]["ac"] = "ac";
    dictionary["abxcd"]["ad"] = "ad";
    dictionary["abxcd"]["bc"] = "bc";
    dictionary["abxcd"]["bd"] = "bd";

    return;
}

void 
load_generic_dictionaries(map<string, map<string, string> > &global_dictionary) 
{
    global_dictionary["ab/--"]["a"]  = "aa";
    global_dictionary["ab/--"]["b"]  = "bb";

    global_dictionary["--/ab"]["a"]  = "aa";
    global_dictionary["--/ab"]["b"]  = "bb";

    global_dictionary["aa/bb"]["a"]  = "-";
    global_dictionary["aa/bb"]["ab"] = "ab";
    global_dictionary["aa/bb"]["b"]  = "-";

    global_dictionary["ab/ac"]["a"]  = "aa";
    global_dictionary["ab/ac"]["ab"] = "ab";
    global_dictionary["ab/ac"]["b"]  = "-";
    global_dictionary["ab/ac"]["ac"] = "ac";
    global_dictionary["ab/ac"]["c"]  = "-";
    global_dictionary["ab/ac"]["bc"] = "bc";

    global_dictionary["ab/cd"]["a"]  = "-";
    global_dictionary["ab/cd"]["ab"] = "-";
    global_dictionary["ab/cd"]["b"]  = "-";
    global_dictionary["ab/cd"]["c"]  = "-";
    global_dictionary["ab/cd"]["cd"] = "-";
    global_dictionary["ab/cd"]["d"]  = "-";
    global_dictionary["ab/cd"]["ac"] = "ac";
    global_dictionary["ab/cd"]["ad"] = "ad";
    global_dictionary["ab/cd"]["bc"] = "bc";
    global_dictionary["ab/cd"]["bd"] = "bd";

    global_dictionary["ab/aa"]["a"]  = "aa";
    global_dictionary["ab/aa"]["ab"] = "ab";
    global_dictionary["ab/aa"]["b"]  = "-";

    global_dictionary["aa/ab"]["a"]  = "aa";
    global_dictionary["aa/ab"]["ab"] = "ab";
    global_dictionary["aa/ab"]["b"]  = "-";

    global_dictionary["ab/cc"]["a"]  = "-";
    global_dictionary["ab/cc"]["ab"] = "-";
    global_dictionary["ab/cc"]["b"]  = "-";
    global_dictionary["ab/cc"]["c"]  = "-";
    global_dictionary["ab/cc"]["ac"] = "ac";
    global_dictionary["ab/cc"]["bc"] = "bc";

    global_dictionary["cc/ab"]["a"]  = "-";
    global_dictionary["cc/ab"]["ab"] = "-";
    global_dictionary["cc/ab"]["b"]  = "-";
    global_dictionary["cc/ab"]["c"]  = "-";
    global_dictionary["cc/ab"]["ac"] = "ac";
    global_dictionary["cc/ab"]["bc"] = "bc";

    global_dictionary["ab/ab"]["a"]  = "aa";
    global_dictionary["ab/ab"]["b"]  = "bb";
    global_dictionary["ab/ab"]["ab"] = "ab";
}

inline
int 
encode_gtype(char a) 
{ 
    switch (a) {
    case 'A':
	return 0;
    case 'C':
	return 1;
    case 'G':
	return 2;
    case 'T':
	return 3;
    }

    return -1;
}

int 
encoded_gtypes[4][4] = 
{
  // A  C  G   T
    {1, 2, 3,  4}, // A
    {2, 5, 6,  7}, // C
    {3, 6, 8,  9}, // G
    {4, 7, 9, 10}  // T
};

#endif // __GENOTYPE_DICTIONARIES_H__
