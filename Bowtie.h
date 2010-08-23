// -*-mode:c++; c-style:k&r; c-basic-offset:4;-*-
//
// Code to parse Bowtie's alignment format. This format is created for
// reads that have been aligned to a reference genome. It takes the tab-separated form:
//
// <query> <strand> <chromosome> <base pair> <sequence> <phred quality score> <flag> <mismatches>
//
// One record per line.
//
#ifndef __BOWTIE_H__
#define __BOWTIE_H__

#include "input.h"

class Bowtie: public Input {

 public:
    Bowtie(const char *path) : Input(path) {};
    ~Bowtie() {};
    Seq *next_seq();
};

Seq *Bowtie::next_seq() {
    vector<string> parts;

    //
    // Read a record from the file and place it in a Seq object
    //
    this->fh.getline(this->line, max_len);

    if (!this->fh.good()) {
	return NULL;
    }

    parse_tsv(this->line, parts);

    Seq *s = new Seq(parts[0].c_str(), parts[4].c_str(), parts[5].c_str(), parts[2].c_str(), atoi(parts[3].c_str()));

    return s;
}

#endif // __BOWTIE_H__
