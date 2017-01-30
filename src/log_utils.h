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

#include <string>
#include <iostream>
#include <fstream>

#include "constants.h"

int init_log(std::ostream &fh, int argc, char **argv);

// Returns e.g. "23.2%".
std::string as_percentage(double d);

std::string to_string(const FileT& ft);

// TeeBuf
// ==========
// A streambuf which tees to two streambufs.
// From http://wordaligned.org/articles/cpp-streambufs */
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
class TeeBuf;
class LogAlterator {
public:
    std::ofstream l; // The actual log file
    std::ostream o;  // Just stdout
    std::ostream e;  // Just stderr

    // Construct the log alterator.
    // std::cout and std::cerr will also write to the log file (if quiet
	// is true, output to stdout and stderr is suppressed i.e. they only
    // write to the log file).
    LogAlterator(const std::string& log_path, bool quiet = false);

    // Upon destruction, restore cout and cerr.
    ~LogAlterator() {
        std::cout.rdbuf(o.rdbuf());
        std::cerr.rdbuf(e.rdbuf());
    }

private:
    TeeBuf lo_buf;
    TeeBuf le_buf;
};

#endif // __LOG_UTILS_H__
