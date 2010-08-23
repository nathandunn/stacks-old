// -*-mode:c++; c-style:k&r; c-basic-offset:4;-*-
//
// Code to parse Sam format. This format is created for
// reads that have been aligned to a reference genome. It takes the tab-separated form:
//
// <query> <strand> <chromosome> <base pair> ... <sequence> <phred quality score> ...
//
// One record per line.
//
#ifndef __SAM_H__
#define __SAM_H__

#include "input.h"

class Sam: public Input {

 public:
    Sam(const char *path) : Input(path) {};
    ~Sam() {};
    Seq *next_seq();
};

Seq *Sam::next_seq() {
    vector<string> parts;

    //
    // Read a record from the file and place it in a Seq object
    //
    this->fh.getline(this->line, max_len);

    if (!this->fh.good()) {
	return NULL;
    }

    parse_tsv(this->line, parts);

    Seq *s = new Seq(parts[0].c_str(), parts[9].c_str(), parts[10].c_str(), parts[2].c_str(), atoi(parts[3].c_str()));

    return s;
}

#endif // __SAM_H__
