// -*-mode:c++; c-style:k&r; c-basic-offset:4;-*-
//
// Code to parse the internal (and tempoary) data format. This format is created for
// reads that have been aligned to a reference genome. It takes the tab-separated form:
//
// <chromosome> <base pair> <sequence> <phred quality score>
//
// One record per line.
//
#ifndef __TSV_H__
#define __TSV_H__

#include "input.h"

class Tsv: public Input {

 public:
    Tsv(const char *path) : Input(path) {};
    ~Tsv() {};
    Seq *next_seq();
};

Seq *Tsv::next_seq() {
    vector<string> parts;

    //
    // Read a record from the file and place it in a Seq object
    //
    this->fh.getline(this->line, max_len);

    if (!this->fh.good()) {
	return NULL;
    }

    parse_tsv(this->line, parts);

    string id = parts[0] + "_" + parts[1];

    Seq *s = new Seq(id.c_str(), parts[2].c_str(), parts[3].c_str(), parts[0].c_str(), atoi(parts[1].c_str()));

    return s;
}

#endif // __TSV_H__
