// -*-mode:c++; c-style:k&r; c-basic-offset:4;-*-
//
// Copyright 2013, Julian Catchen <jcatchen@uoregon.edu>
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
// catalog_utils.cc -- common routines for manipulating catalog objects.
//
// Julian Catchen
// jcatchen@uoregon.edu
// University of Oregon
//
#include "catalog_utils.h"

int 
reduce_catalog(map<int, CSLocus *> &catalog, set<int> &whitelist, set<int> &blacklist) 
{
    map<int, CSLocus *> list;
    map<int, CSLocus *>::iterator it;
    CSLocus *loc;

    if (whitelist.size() == 0 && blacklist.size() == 0) 
	return 0;
 
    int i = 0;
    for (it = catalog.begin(); it != catalog.end(); it++) {
	loc = it->second;

	if (whitelist.size() > 0 && whitelist.count(loc->id) == 0) continue;
	if (blacklist.count(loc->id)) continue;

	list[it->first] = it->second;
	i++;
    }

    catalog = list;

    return i;
}

