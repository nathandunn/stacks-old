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
#include <tr1/unordered_map>
using std::tr1::unordered_map;

#include "input.h"
#include "kmers.h"

enum fastqt   {generic_fastq, illv1_fastq, illv2_fastq};

enum barcodet {null_null,
	       inline_null,   index_null, 
	       inline_inline, index_index, 
	       inline_index,  index_inline};

extern int      bc_size_1, bc_size_2;
extern barcodet barcode_type;
extern int      truncate_seq;
extern double   win_size;
extern bool     paired;
extern bool     recover;
extern int      barcode_dist;

class BarcodePair {
public:
    string se;    // Single-end barcode
    string pe;    // Paired-end barcode

    BarcodePair()
    {
	this->se = "";
	this->pe = "";
    }
    BarcodePair(char *p, char *q)
    {
	this->se = string(p);
	this->pe = string(q);
    }
    BarcodePair(string se, string pe)
    {
	this->se = se;
	this->pe = pe;
    }
    void set(char *p, char *q)
    {
	this->se = string(p);
	this->pe = string(q);
    }
    void set(char *p)
    {
	this->se = string(p);
	this->pe = "";
    }
    string str() 
    {
	if (this->pe.length() > 0)
	    return string(this->se + "-" + this->pe);
	else
	    return this->se;
    }
    friend bool operator<(const BarcodePair &lhs, const BarcodePair &rhs)
    {
	if (lhs.se < rhs.se)
	    return true;
	else if (lhs.se == rhs.se && lhs.pe < rhs.pe)
	    return true;
	else
	    return false;
    }
    friend bool operator==(const BarcodePair &lhs, const BarcodePair &rhs)
    {
	return (lhs.se == rhs.se && lhs.pe == rhs.pe);
    }
    friend ofstream& operator<<(ofstream &out, const BarcodePair &bp)
    {
	if (bp.pe.length() > 0)
	    out << bp.se << "-" << bp.pe;
	else
	    out << bp.se;
	return out;
    }
};

class Read {
public:
    fastqt  fastq_type;
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
    bool    filter;
    int     inline_bc_len;
    int     retain;
    unsigned int size;
    unsigned int len;
    double  win_len;
    double  stop_pos;

    Read(uint buf_len, int read, int barcode_size, double win_size) {
	this->barcode       = new char[id_len  + 1];
	this->machine       = new char[id_len  + 1];
	this->seq           = new char[buf_len + 1];
	this->phred         = new char[buf_len + 1];
	this->int_scores    = new  int[buf_len];
	this->size          = buf_len + 1;
	this->inline_bc_len = barcode_size;
	this->read          = read;

	this->retain        = 1;
	this->tile          = 0;
	this->lane          = 0;
	this->x             = 0;
	this->y             = 0;
	this->index         = 0;

	this->barcode[0] = '\0';
	this->machine[0] = '\0';
	this->seq[0]     = '\0';
	this->phred[0]   = '\0';

	this->set_len(buf_len);
    }
    ~Read() {
	delete [] this->barcode;
	delete [] this->machine;
	delete [] this->seq;
	delete [] this->phred;
	delete [] this->int_scores;
    }
    int resize(int size) {
	delete [] this->seq;
	delete [] this->phred;
	delete [] this->int_scores;
	this->size  = size;
	this->seq   = new char[this->size];
	this->phred = new char[this->size];
	this->int_scores = new int[this->len];

	this->set_len(size - 1);

	return 0;
    }
    int set_len(uint buf_len) {
	if (buf_len == this->len)
	    return 0;

	//
	// Set the parameters for checking read quality later in processing.
	// Window length is 15% (rounded) of the sequence length.
	//
	this->len      = buf_len - this->inline_bc_len;
	this->win_len  = round((double) this->len * win_size);

	if (this->win_len < 1) 
	    this->win_len = 1;

	this->len     += this->inline_bc_len;
	this->stop_pos = this->len - this->win_len;

	return 0;
    }
};

typedef unordered_map<const char *, vector<int>, std::tr1::hash<const char *>, eqstr> AdapterHash;

int  parse_illumina_v1(const char *);
int  parse_illumina_v2(const char *);
int  parse_input_record(Seq *, Read *);
int  rev_complement(char *, int, bool);
int  reverse_qual(char *, int, bool);

int  process_barcode(Read *, Read *, BarcodePair &, 
		     map<BarcodePair, ofstream *> &,
		     set<string> &, set<string> &, 
		     map<BarcodePair, map<string, long> > &, map<string, long> &); 
bool correct_barcode(set<string> &, Read *);

int  filter_adapter_seq(Read *, char *, int, AdapterHash &, int, int, int);
int  init_adapter_seq(int, char *, int &, AdapterHash &, vector<char *> &);
int  free_adapter_seq(vector<char *> &);

int  check_quality_scores(Read *, int, int, int, int);


#endif // __CLEAN_H__
