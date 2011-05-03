// -*-mode:c++; c-style:k&r; c-basic-offset:4;-*-
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

#ifndef __SAM_H__
#define __SAM_H__

//
// Code to parse Sam format. This format is created for
// reads that have been aligned to a reference genome. It takes the tab-separated form:
//
// <query> <strand> <chromosome> <base pair> ... <sequence> <phred quality score> ...
//
// One record per line.
//

#include "input.h"

class Sam: public Input {

 public:
    Sam(const char *path) : Input(path) {};
    ~Sam() {};
    Seq *next_seq();
    Seq *next_seq(Seq *) {};
};

Seq *Sam::next_seq() {
    vector<string> parts;

    //
    // Read a record from the file and place it in a Seq object, skipping header definitions.
    //
    do {
        this->fh.getline(this->line, max_len);

        if (!this->fh.good())
            return NULL;

        parse_tsv(this->line, parts);

    } while (parts[0][0] == '@');

    Seq *s = new Seq(parts[0].c_str(), parts[9].c_str(), parts[10].c_str(), parts[2].c_str(), atoi(parts[3].c_str()));

    return s;
}

#endif // __SAM_H__
