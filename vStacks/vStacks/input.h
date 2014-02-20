// -*-mode:c++; c-style:k&r; c-basic-offset:4;-*-
//
// Copyright (c) 2014 University of Oregon
// Created by Julian Catchen <jcatchen@uoregon.edu>
//

#ifndef __INPUT_H__
#define __INPUT_H__

#include <stdlib.h>
#include <string.h>
#include <string>
#include <vector>
#include <iostream>
#include <fstream>
using std::ifstream;
using std::string;
using std::vector;
using std::cin;
using std::cout;
using std::cerr;
using std::endl;

#include "constants.h"
#include "stacks.h"

typedef unsigned int uint;

class Seq {
 public:
    char *id;
    char *seq;
    char *qual;

    // Location information for a mapped sequence
    char  *loc_str;
    PhyLoc loc;

    Seq( void );
    Seq(const char *, const char *);
    Seq(const char *, const char *, const char *);
    Seq(const char *, const char *, const char *, const char *, uint, strand_type);
    ~Seq( void ) { 
	delete[] id; 
	delete[] seq; 
	delete[] qual; 
	delete[] loc_str; 
    }
};

//
// The base class for all of our Input classes, such as Tsv, Fastq, Fasta, etc.
//
class Input {
 public:
    string   path;
    ifstream fh;
    char     line[max_len];

    Input();
    Input(const char *path);
    virtual ~Input();
    virtual Seq *next_seq() = 0;
    virtual int  next_seq(Seq &) = 0;
};

char *rev_comp(const char *);
int   parse_tsv(const char *, vector<string> &);
int   parse_ssv(const char *, vector<string> &);
int   read_line(ifstream &, char **, int *);

#endif // __INPUT_H__
