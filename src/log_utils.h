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
#include <sstream>

#include "config.h"
#include "constants.h"

inline
int init_log(std::ostream &fh, int argc, char **argv) {
    //
    // Obtain the current date.
    //
    time_t     rawtime;
    struct tm *timeinfo;
    char       date[32];
    time(&rawtime);
    timeinfo = localtime(&rawtime);
    strftime(date, 32, "%F %T", timeinfo);

    //
    // Write the command line that was executed.
    //
    for (int i = 0; i < argc; i++) {
        fh << argv[i];
        if (i < argc - 1) fh << " ";
    }
    fh << "\n" << argv[0] << " version " << VERSION << " executed " << date << "\n\n";

    return 0;
}

// TeeBuf
// ==========
// From http://wordaligned.org/articles/cpp-streambufs */
// A streambuf which tees to two streambufs.
// This tee buffer has no actual buffer, so every character "overflows"
// directly into the teed buffers.
class TeeBuf: public std::streambuf {
public:

    TeeBuf(std::streambuf* sb1, std::streambuf* sb2)
        : sb1(sb1) , sb2(sb2)
        {}

private:
    std::streambuf* sb1;
    std::streambuf* sb2;

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
    std::ofstream l; // The actual log file
    std::ostream o; // Just stdout
    std::ostream e; // Just stderr
    TeeBuf lo_buf;
    TeeBuf le_buf;

    // Constructor
    // ----------
    // Construct the log fstream;
    // Keep track of the streambufs of cout and cerr;
    // Construst the teeing streambufs and let cout and cerr use them.
    LogAlterator(const std::string& log_path, bool quiet = false)
            : l(log_path)
            , o(std::cout.rdbuf())
            , e(std::cerr.rdbuf())
            , lo_buf(std::cout.rdbuf(), l.rdbuf())
            , le_buf(std::cerr.rdbuf(), l.rdbuf())
            {
        if (quiet) {
            // Use the fstream buffer only, and not at all stdout and stderr.
            std::cout.rdbuf(l.rdbuf());
            std::cerr.rdbuf(l.rdbuf());
        } else {
            std::cout.rdbuf(&lo_buf);
            std::cerr.rdbuf(&le_buf);
        }
    }

    // Destructor
    // ---------
    // Restore cout and cerr.
    ~LogAlterator() {
        std::cout.rdbuf(o.rdbuf());
        std::cerr.rdbuf(e.rdbuf());
    }
};

#endif // __LOG_UTILS_H__
