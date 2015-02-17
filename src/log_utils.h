// -*-mode:c++; c-style:k&r; c-basic-offset:4;-*-
//
// Copyright 2013-2015, Julian Catchen <jcatchen@illinois.edu>
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

#ifndef __LOG_UTILS_H__
#define __LOG_UTILS_H__

#include <time.h>
#include <iostream>
#include <fstream>
using std::ifstream;
using std::ofstream;
using std::cin;
using std::cout;
using std::cerr;
using std::endl;
#include <sstream>
using std::stringstream;

#include "constants.h"

int init_log(ofstream &, int, char **);

#endif // __LOG_UTILS_H__
