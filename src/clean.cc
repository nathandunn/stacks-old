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

int parse_input_record(Seq *s, Read *r) {
    char *p, *q;
    //
    // Count the number of colons to differentiate Illumina version.
    // CASAVA 1.8+ has a FASTQ header like this:
    //  @HWI-ST0747:155:C01WHABXX:8:1101:6455:26332 1:N:0:
    //
    //
    // Or, parse FASTQ header from previous versions that looks like this:
    //  @HWI-ST0747_0141:4:1101:1240:2199#0/1
    //  @HWI-ST0747_0143:2:2208:21290:200914#0/1
    //
    char *stop = s->id + strlen(s->id);
    int   cnt  = 0;

    for (p = s->id, q = p; q < stop; q++)
	cnt += *q == ':'? 1 : 0;

    if (cnt == 9) {
	//
	// According to Illumina manual, "CASAVA v1.8 User Guide" page 41:
	// @<instrument>:<run number>:<flowcell ID>:<lane>:<tile>:<x-pos>:<y-pos> <read>:<is filtered>:<control number>:<index sequence>
	//
	for (p = s->id, q = p; *q != ':' && q < stop; q++);
	*q = '\0';
	strcpy(r->machine, p);
	*q = ':';

	// Run number.
	for (p = q+1, q = p; *q != ':' && q < stop; q++);
	//*q = '\0';

	// Flowcell ID.
	for (p = q+1, q = p; *q != ':' && q < stop; q++);
	//*q = '\0';

	for (p = q+1, q = p; *q != ':' && q < stop; q++);
	*q = '\0';
	r->lane = atoi(p);
	*q = ':';

	for (p = q+1, q = p; *q != ':' && q < stop; q++);
	*q = '\0';
	r->tile = atoi(p);
	*q = ':';

	for (p = q+1, q = p; *q != ':' && q < stop; q++);
	*q = '\0';
	r->x = atoi(p);
	*q = ':';

	for (p = q+1, q = p; *q != ' ' && q < stop; q++);
	*q = '\0';
	r->y = atoi(p);
	*q = ' ';

	for (p = q+1, q = p; *q != ':' && q < stop; q++);
	*q = '\0';
	r->read = atoi(p);
	*q = ':';

	for (p = q+1, q = p; *q != ':' && q < stop; q++);
	*q = '\0';
	r->filter = *p == 'Y' ? true : false;
	*q = ':';

    } else {
	for (p = s->id, q = p; *q != ':' && q < stop; q++);
	*q = '\0';
	strcpy(r->machine, p);
	*q = ':';

	for (p = q+1, q = p; *q != ':' && q < stop; q++);
	*q = '\0';
	r->lane = atoi(p);
	*q = ':';

	for (p = q+1, q = p; *q != ':' && q < stop; q++);
	*q = '\0';
	r->tile = atoi(p);
	*q = ':';

	for (p = q+1, q = p; *q != ':' && q < stop; q++);
	*q = '\0';
	r->x = atoi(p);
	*q = ':';

	for (p = q+1, q = p; *q != '#' && q < stop; q++);
	*q = '\0';
	r->y = atoi(p);
	*q = '#';

	for (p = q+1, q = p; *q != '/' && q < stop; q++);
	*q = '\0';
	r->index = atoi(p);
	*q = '/';

	for (p = q+1, q = p; *q != '\0' && q < stop; q++);
	r->read = atoi(p);
    }

    strncpy(r->seq,     s->seq,  r->len); 
    r->seq[r->len]   = '\0';
    strncpy(r->phred,   s->qual, r->len);
    r->phred[r->len] = '\0';

    if (barcode_size > 0) {
	strncpy(r->barcode, r->seq,  barcode_size);
	r->barcode[barcode_size] = '\0';
    }

    r->retain = 1;

    return 0;
}

int write_fasta(map<string, ofstream *> &fhs, Read *href, bool barcode, bool overhang) {
    char tile[id_len];
    sprintf(tile, "%04d", href->tile);

    int offset;
    offset  = barcode  ? barcode_size : 0;
    offset += overhang ? 1 : 0;

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

int write_fastq(map<string, ofstream *> &fhs, Read *href, bool barcode, bool overhang) {
    //
    // Write the sequence and quality scores in FASTQ format. 
    //
    char tile[id_len];
    sprintf(tile, "%04d", href->tile);

    int offset;
    offset  = barcode  ? barcode_size : 0;
    offset += overhang ? 1 : 0;

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

int rev_complement(char *seq, bool barcode, bool overhang) {
    char *p, *q;
    int offset;
    offset  = barcode  ? barcode_size : 0;
    offset += overhang ? 1 : 0;
    q       = seq + offset;

    int len   = strlen(q);
    int j     = 0;
    char *com = new char[len + 1]; 
   
    for (p = q + len - 1; p >= q; p--) {
        switch (*p) {
        case 'A':
        case 'a':
            com[j] = 'T';
            break;
        case 'C':
        case 'c':
            com[j] = 'G';
            break;
        case 'G':
        case 'g':
            com[j] = 'C';
            break;
        case 'T':
        case 't':
            com[j] = 'A';
            break;
        }
        j++;
    }
    com[len] = '\0';

    for (j = 0; j < len; j++)
	q[j] = com[j];

    delete [] com;

    return 0;
}

int reverse_qual(char *qual, bool barcode, bool overhang) {
    char *p, *q;
    int offset;
    offset  = barcode  ? barcode_size : 0;
    offset += overhang ? 1 : 0;
    q       = qual + offset;

    int len   = strlen(q);
    int j     = 0;
    char *com = new char[len + 1]; 
   
    for (p = q + len - 1; p >= q; p--) {
	com[j] = *p;
        j++;
    }
    com[len] = '\0';

    for (j = 0; j < len; j++)
	q[j] = com[j];

    delete [] com;

    return 0;
}
