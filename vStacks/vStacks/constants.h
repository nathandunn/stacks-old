// -*-mode:c++; c-style:k&r; c-basic-offset:4;-*-
//
// Copyright (c) 2014 University of Oregon
// Created by Julian Catchen <jcatchen@uoregon.edu>
//

#ifndef __CONSTANTS_H__
#define __CONSTANTS_H__

//
// Pull in the configuration variables from the configure script
//
#if HAVE_CONFIG_H
#include "config.h"
#endif

//
// Maximum line length for parsing input files.
//
const int max_len = 1024;

//
// Maximum length of idetifiers, such as sequence IDs and chromosome names.
//
const int id_len = 255;

//
// Supported file types
//
enum file_type {unknown, sql, fasta, fastq, gzfasta, gzfastq, bowtie, sam, bam, tsv, bustard, phase, fastphase, beagle};

#endif
