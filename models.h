// -*-mode:c++; c-style:k&r; c-basic-offset:4;-*-

#ifndef __MODELS_H__
#define __MODELS_H__

#include <math.h>
#include <string>
using std::string;
#include <vector>
using std::vector;
#include <map>
using std::map;
#include <utility>
using std::pair;
using std::make_pair;
#include <iostream>
using std::ifstream;
using std::cin;
using std::cout;
using std::cerr;
using std::endl;

#include "stacks.h"

//
// Possible models for calling nucleotide positions as fixed or variable
//
enum modelt {fixed, snp};

//
// For use with the multinomial model to call fixed nucleotides.
//
extern const int barcode_size;
extern const int snp_min;
extern modelt model_type;
extern double p_freq;
extern double barcode_err_freq;

int  call_multinomial_snp(MergedStack *, int, map<char, int> &);
int  call_multinomial_fixed(MergedStack *, int, map<char, int> &);
bool compare_pair(pair<char, int>, pair<char, int>);

#endif // __MODELS_H__
