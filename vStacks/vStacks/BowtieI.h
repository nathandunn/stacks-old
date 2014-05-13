// -*-mode:c++; c-style:k&r; c-basic-offset:4;-*-
//
// Copyright (c) 2014 University of Oregon
// Created by Julian Catchen <jcatchen@uoregon.edu>
//
//

#ifndef __BOWTIEI_H__
#define __BOWTIEI_H__

//
// Code to parse Bowtie's alignment format. This format is created for
// reads that have been aligned to a reference genome. It takes the tab-separated form:
//
// <query> <strand> <chromosome> <base pair> <sequence> <phred quality score> <flag> <mismatches>
//
// One record per line.
//

#include "input.h"

class Bowtie: public Input {

 public:
    Bowtie(const char *path) : Input(path) {};
    ~Bowtie() {};
    Seq *next_seq();
    int  next_seq(Seq &) { return 0; };
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

    strand_type strand = parts[1] == "+" ? plus : minus;

    //
    // If the read was aligned on the reverse strand (and is therefore reverse complemented)
    // alter the start point of the alignment to reflect the right-side of the read, at the
    // end of the RAD cut site.
    //
    int bp = strand == plus ? atoi(parts[3].c_str()) : atoi(parts[3].c_str()) + parts[4].length();

    Seq *s = new Seq(parts[0].c_str(), parts[4].c_str(), parts[5].c_str(), 
		     parts[2].c_str(), bp, strand);

    return s;
}

#endif // __BOWTIEI_H__
