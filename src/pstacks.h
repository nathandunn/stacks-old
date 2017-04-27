// -*-mode:c++; c-style:k&r; c-basic-offset:4;-*-
//
// Copyright 2010-2015, Julian Catchen <jcatchen@illinois.edu>
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

#ifndef __PSTACKS_H__
#define __PSTACKS_H__

#ifdef _OPENMP
#include <omp.h>    // OpenMP library
#endif

#include <cerrno>
#include <zlib.h>   // Support for gzipped output files.

#include <getopt.h> // Process command-line options
#include <cstring>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip> // std::setprecision

#include <unordered_map>
#include <vector>
#include <map>
#include <set>
#include <utility>

#include "config.h"
#include "constants.h"
#include "stacks.h"     // Major data structures for holding stacks
#include "mstack.h"
#include "kmers.h"
#include "utils.h"
#include "models.h"     // Contains maximum likelihood statistical models.
#include "Tsv.h"        // Reading input files in Tab-separated values format
#include "BowtieI.h"    // Reading input files in Bowtie format
#include "SamI.h"       // Reading input files in SAM format
#include "BamI.h"       // Reading input files in BAM format
#include "DNANSeq.h"

#ifdef HAVE_SPARSEHASH
typedef sparse_hash_map<DNANSeq, vector<Seq*> > HashMap;
#else
typedef unordered_map<DNANSeq, vector<Seq*> > HashMap;
#endif

void   help( void );
void   version( void );
int    parse_command_line(int, char**);
void   load_radtags(string, HashMap &);
int    reduce_radtags(HashMap &, map<int, PStack *> &);
void   populate_merged_tags(map<int, PStack *>& unique, map<int, MergedStack *>& merged);
void   delete_low_cov_loci(map<int, MergedStack *>&, const map<int, PStack *>&);

void report_options(ostream& os);

//
// Debugging
//
int  dump_stacks(map<int, PStack *> &);
int  dump_merged_stacks(map<int, MergedStack *> &);

#endif // __PSTACKS_H__
