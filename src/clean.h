// -*-mode:c++; c-style:k&r; c-basic-offset:4;-*-
//
// Copyright 2011, Julian Catchen <jcatchen@uoregon.edu>
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

#ifndef __CLEAN_H__
#define __CLEAN_H__

#include <string>
using std::string;
#include <map>
using std::map;
#include <iostream>
#include <fstream>
using std::ofstream;
using std::cin;
using std::cout;
using std::cerr;
using std::endl;

#include "input.h"

typedef struct read {
    char   *barcode;
    char   *machine;
    int     lane;
    int     tile;
    int     x;
    int     y;
    int     index;
    int     read;
    char   *seq;
    char   *phred;
    int    *int_scores;
    int     filter;
    int     retain;
    unsigned int len;
    double  win_len;
    double  stop_pos;
} Read;

extern int barcode_size;

int  parse_illumina_v1(const char *);
int  parse_illumina_v2(const char *);
int  parse_input_record(Seq *, Read *);
int  write_fastq(map<string, ofstream *> &, Read *, bool, bool);
int  write_fastq(ofstream *, Seq *);
int  write_fastq(ofstream *, Seq *, string);
int  write_fasta(map<string, ofstream *> &, Read *, bool, bool);
int  write_fasta(ofstream *, Seq *);
int  write_fasta(ofstream *, Seq *, string);
int  rev_complement(char *, bool, bool);
int  reverse_qual(char *, bool, bool);

#endif // __CLEAN_H__
