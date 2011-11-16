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

//
// clean.cc -- common routines for processing and cleaning raw seqeunce data.
//
// Julian Catchen
// jcatchen@uoregon.edu
// University of Oregon
//
// $Id: clean.cc 2146 2011-08-02 22:11:50Z catchen $
//

#include "clean.h"

int parse_illumina_v1(const char *file) {
    const char *p, *q;

    //
    // Parse a file name that looks like: s_7_1_0001_qseq.txt ... s_7_1_0120_qseq.txt
    // but exclude the paired-end files:  s_7_2_0001_qseq.txt ... s_7_2_0120_qseq.txt
    //
    if (file[0] != 's') 
	return 0;

    int underscore_cnt = 0;
    for (p = file; *p != '\0'; p++) {
	if (*p == '_') {
	    underscore_cnt++;
	    q = p;
	}
    }

    if (underscore_cnt != 4)
	return 0;

    // Check the file suffix.
    if (strncmp(q, "_qseq.txt", 8) != 0)
	return 0;

    // Make sure it is not the paired-end file
    p  = file;
    p += 3;
    if (strncmp(p, "_1_", 3) != 0)
	return 0;

    //
    // Return the position of the paired-end number, so the other file name can be generated.
    //
    return (p + 1 - file);
}

int parse_illumina_v2(const char *file) {
    const char *p, *q;

    //
    // Parse a file name that looks like: lane6_NoIndex_L006_R1_003.fastq
    // but exclude the paired-end files:  lane6_NoIndex_L006_R2_003.fastq
    //
    if (strncmp(file, "lane", 4) != 0) 
	return 0;

    //
    // Make sure it ends in "fastq"
    //
    for (q = file; *q != '\0'; q++);
    for (p = q; *p != '.' && p > file; p--);
    if (strncmp(p, ".fastq", 6) != 0)
	return 0;

    int underscore_cnt = 0;
    for (p = file; *p != '\0'; p++) {
	if (*p == '_') underscore_cnt++;
	q = p;
    }

    if (underscore_cnt != 4)
	return 0;

    // Check the lane exists
    for (p = file; *p != '_' && *p != '\0'; p++);
    p++;
    for (; *p != '_' && *p != '\0'; p++);
    p++;

    if (*p != 'L')
	return 0;

    // Make sure it is not the paired-end file
    for (; *p != '_' && *p != '\0'; p++);
    if (strncmp(p, "_R1_", 4) != 0)
	return 0;

    //
    // Return the position of the paired-end number, so the other file name can be generated.
    //
    return (p + 2 - file); 
}

int write_fasta(map<string, ofstream *> &fhs, Read *href, bool paired_end) {
    char tile[id_len];
    sprintf(tile, "%04d", href->tile);

    int offset = paired_end ? 0 : barcode_size;

    *(fhs[href->barcode]) <<
	">" << href->barcode <<
	"_" << href->lane <<
	"_" << tile << 
	"_" << href->x <<
	"_" << href->y <<
	"_" << href->read << "\n" <<
	href->seq + offset << "\n";

    return 0;
}

int write_fasta(ofstream *fh, Seq *href) {
    *fh <<
	">" << 
	href->id  << "\n" <<
	href->seq << "\n";

    return 0;
}

int write_fastq(map<string, ofstream *> &fhs, Read *href, bool paired_end) {
    //
    // Write the sequence and quality scores in FASTQ format. 
    //
    char tile[id_len];
    sprintf(tile, "%04d", href->tile);

    int offset = paired_end ? 0 : barcode_size;

    if (fhs.count(href->barcode) == 0)
	cerr << "Writing to unknown barcode: '" << href->barcode << "'\n";

    *(fhs[href->barcode]) <<
	"@" << href->barcode <<
	"_" << href->lane << 
	"_" << tile << 
	"_" << href->x << 
	"_" << href->y << 
	"_" << href->read << "\n" <<
	href->seq + offset << "\n" <<
	"+\n" <<
	href->phred + offset << "\n";

    return 0;
}

int write_fastq(ofstream *fh, Seq *href) {
    *fh <<
	"@" << href->id << "\n" <<
	href->seq << "\n" <<
	"+\n" <<
	href->qual << "\n";

    return 0;
}

int write_fastq(ofstream *fh, Seq *href, string msg) {
    *fh <<
	"@" << href->id << "|" << msg << "\n" <<
	href->seq << "\n" <<
	"+\n" <<
	href->qual << "\n";

    return 0;
}

int write_fasta(ofstream *fh, Seq *href, string msg) {
    *fh <<
	">" << 
	href->id  << "|" << msg << "\n" <<
	href->seq << "\n";

    return 0;
}
