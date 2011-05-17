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

void initialize_dictionaries(map<string, map<string, string> > &global_dictionary) {
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
    global_dictionary["ab/ac"]["ab"] = "ab";
    global_dictionary["ab/ac"]["c"]  = "cc";
    global_dictionary["ab/ac"]["ac"] = "ac";

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
    global_dictionary["aa/ab"]["ab"] = "ab";
    global_dictionary["aa/ab"]["bb"] = "bb";
    global_dictionary["aa/ab"]["c"]  = "cc";
    global_dictionary["aa/ab"]["ac"] = "ac";
    global_dictionary["aa/ab"]["bc"] = "bc";

    global_dictionary["cc/ab"]["aa"] = "aa";
    global_dictionary["cc/ab"]["ab"] = "ab";
    global_dictionary["cc/ab"]["bb"] = "bb";
    global_dictionary["cc/ab"]["c"]  = "cc";
    global_dictionary["cc/ab"]["ac"] = "ac";
    global_dictionary["cc/ab"]["bc"] = "bc";

    global_dictionary["ab/ab"]["a"]  = "aa";
    global_dictionary["ab/ab"]["b"]  = "bb";
    global_dictionary["ab/ab"]["ab"] = "ab";
}

#endif // __GENOTYPE_DICTIONARIES_H__
