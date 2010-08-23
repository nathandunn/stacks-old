// -*-mode:c++; c-style:k&r; c-basic-offset:4;-*-

#ifndef __INPUT_H__
#define __INPUT_H__

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

class Seq {
 public:
    char *id;
    char *seq;
    char *qual;

    // Location information for a mapped sequence
    char *loc_str;
    char *chr;
    uint  bp;

    Seq( void );
    Seq(const char *, const char *);
    Seq(const char *, const char *, const char *);
    Seq(const char *, const char *, const char *, const char *, uint);
    ~Seq( void ) { delete[] id; delete[] seq; delete[] qual; delete[] chr; delete[] loc_str; }
};

//
// The base class for all of our Input classes, such as Tsv, Fastq, Fasta, etc.
//
class Input {
 public:
    ifstream fh;
    char     line[max_len];

    Input(const char *path);
    virtual ~Input();
    virtual Seq *next_seq() = 0;
};


int parse_tsv(const char *, vector<string> &);



#endif // __INPUT_H__
