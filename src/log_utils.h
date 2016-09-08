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

// TeeBuf
// ==========
// From http://wordaligned.org/articles/cpp-streambufs */
class TeeBuf: public std::streambuf {
public:
    // Construct a streambuf which tees output to both input
    // streambufs.
    TeeBuf(std::streambuf* sb1, std::streambuf* sb2)
        : sb1(sb1) , sb2(sb2)
        {}

private:
    std::streambuf* sb1;
    std::streambuf* sb2;
    // This tee buffer has no buffer. So every character "overflows"
    // and can be put directly into the teed buffers.
    virtual int overflow(int c) {
        if (c == EOF) {
            return !EOF;
        } else {
            int r1 = sb1->sputc(c);
            int r2 = sb2->sputc(c);
            return r1 == EOF || r2 == EOF ? EOF : c;
        }
    }

    // Sync both teed buffers.
    virtual int sync() {
        int r1 = sb1->pubsync();
        int r2 = sb2->pubsync();
        return r1 == 0 && r2 == 0 ? 0 : -1;
    }
};

// LogAlterator
// ==========
struct LogAlterator {
    std::ofstream logfile;
    std::ostream true_cout;
    std::ostream true_cerr;
    TeeBuf lo;
    TeeBuf le;

    // Constructor
    // ----------
    // Construct the log fstream;
    // Keep track of the streambufs of cout and cerr;
    // Construst the teeing streambufs and let cout and cerr use them.
    LogAlterator(std::ofstream&& logf)
            : logfile(std::move(logf))
            , true_cout(std::cout.rdbuf())
            , true_cerr(std::cerr.rdbuf())
            , lo(std::cout.rdbuf(), logfile.rdbuf())
            , le(std::cerr.rdbuf(), logfile.rdbuf())
            {
        std::cout.rdbuf(&lo);
        std::cerr.rdbuf(&le);
    }

    // Destructor
    // ---------
    // Restore cout and cerr.
    ~LogAlterator() {
        std::cout.rdbuf(true_cout.rdbuf());
        std::cerr.rdbuf(true_cerr.rdbuf());
    }
};

#endif // __LOG_UTILS_H__
