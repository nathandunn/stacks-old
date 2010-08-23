// -*-mode:c++; c-style:k&r; c-basic-offset:4;-*-

#ifndef __EXSTACKS_H__
#define __EXSTACKS_H__

#define VERSION 0.10

#include <omp.h>    // OpenMP library
#include <getopt.h> // Process command-line options
#include <math.h>
#include <stdlib.h>
#include <utility>
using std::pair;
using std::make_pair;

#include <string>
using std::string;

#include <iostream>
#include <fstream>
#include <sstream>
using std::ifstream;
using std::ofstream;
using std::stringstream;
using std::cin;
using std::cout;
using std::cerr;
using std::endl;

#include <vector>
using std::vector;

#include <map>
using std::map;

#include <set>
using std::set;

#include <queue>
using std::queue;

#include "constants.h"
#include "stacks.h"
#include "sql_utilities.h"

void help( void );
void version( void );
int  parse_command_line(int, char**);
int  parse_tsv(const char *, vector<string> &);
bool compare_dist(pair<int, int>, pair<int, int>);
int  write_simple_output(Locus *);
bool compare_loci(Locus *, Locus *);

#endif // __EXSTACKS_H__
