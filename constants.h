// -*-mode:c++; c-style:k&r; c-basic-offset:4;-*-

#ifndef __CONSTANTS_H__
#define __CONSTANTS_H__

//
// Stacks version number
//
const double stacks_version = 0.90;

//
// Maximum line length for parsing input files.
//
const int max_len = 1024;

//
// Maximum length of idetifiers, such as sequence IDs and chromosome names.
//
const int id_len = 64;

//
// Supported file types
//
enum file_type {unknown, sql, fasta, fastq, bowtie, sam, tsv};

#endif
