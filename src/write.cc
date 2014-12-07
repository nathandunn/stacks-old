// -*-mode:c++; c-style:k&r; c-basic-offset:4;-*-
//
// Copyright 2013-2014, Julian Catchen <jcatchen@uoregon.edu>
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
// write.cc -- common routines for writing FASTA/FASTQ records to a file..
//
// Julian Catchen
// jcatchen@uoregon.edu
// University of Oregon
//

#include "write.h"

int 
write_fasta(ofstream *fh, Read *href, bool overhang) {
    char tile[id_len];
    sprintf(tile, "%04d", href->tile);

    int offset = href->inline_bc_len;
    offset += overhang ? 1 : 0;

    if (href->fastq_type != generic_fastq)
    	*fh <<
    	    ">" << href->lane <<
    	    "_" << tile << 
    	    "_" << href->x <<
    	    "_" << href->y <<
    	    "_" << href->read << "\n" <<
    	    href->seq + offset << "\n";
    else 
    	*fh <<
    	    ">" << href->machine <<
    	    "_" << href->read << "\n" <<
    	    href->seq + offset << "\n";

    return 0;
}

int 
write_fasta(gzFile *fh, Read *href, bool overhang) {
    stringstream sstr;
    char tile[id_len];
    sprintf(tile, "%04d", href->tile);

    int offset = href->inline_bc_len;
    offset += overhang ? 1 : 0;

    if (href->fastq_type != generic_fastq)
    	sstr <<
    	    ">" << href->lane <<
    	    "_" << tile << 
    	    "_" << href->x <<
    	    "_" << href->y <<
    	    "_" << href->read << "\n" <<
    	    href->seq + offset << "\n";
    else 
    	sstr <<
    	    ">" << href->machine <<
    	    "_" << href->read << "\n" <<
    	    href->seq + offset << "\n";

    gzputs(*fh, sstr.str().c_str());

    return 0;
}

int 
write_fasta(ofstream *fh, Seq *href) {
    *fh <<
	">" << 
	href->id  << "\n" <<
	href->seq << "\n";

    return 0;
}

int 
write_fasta(gzFile *fh, Seq *href) {
    stringstream sstr;

    sstr <<
	">" << 
	href->id  << "\n" <<
	href->seq << "\n";

    gzputs(*fh, sstr.str().c_str());

    return 0;
}

int 
write_fastq(ofstream *fh, Read *href, bool overhang) {
    //
    // Write the sequence and quality scores in FASTQ format. 
    //
    char tile[id_len];
    sprintf(tile, "%04d", href->tile);

    int offset = href->inline_bc_len;
    offset += overhang ? 1 : 0;

    if (href->fastq_type != generic_fastq)
    	*fh <<
    	    "@" << href->lane << 
    	    "_" << tile << 
    	    "_" << href->x << 
    	    "_" << href->y << 
    	    "_" << href->read << "\n" <<
    	    href->seq + offset << "\n" <<
    	    "+\n" <<
    	    href->phred + offset << "\n";
    else
    	*fh <<
    	    "@" << href->machine << 
    	    "_" << href->read << "\n" <<
    	    href->seq + offset << "\n" <<
    	    "+\n" <<
    	    href->phred + offset << "\n";

    return 0;
}

int 
write_fastq(gzFile *fh, Read *href, bool overhang) {
    //
    // Write the sequence and quality scores in FASTQ format. 
    //
    stringstream sstr;
    char tile[id_len];
    sprintf(tile, "%04d", href->tile);

    int offset = href->inline_bc_len;
    offset += overhang ? 1 : 0;

    if (href->fastq_type != generic_fastq)
    	sstr <<
    	    "@" << href->lane << 
    	    "_" << tile << 
    	    "_" << href->x << 
    	    "_" << href->y << 
    	    "_" << href->read << "\n" <<
    	    href->seq + offset << "\n" <<
    	    "+\n" <<
    	    href->phred + offset << "\n";
    else
    	sstr <<
    	    "@" << href->machine << 
    	    "_" << href->read << "\n" <<
    	    href->seq + offset << "\n" <<
    	    "+\n" <<
    	    href->phred + offset << "\n";

    gzputs(*fh, sstr.str().c_str());

    return 0;
}

int 
write_fastq(ofstream *fh, Seq *href) {
    *fh <<
	"@" << href->id << "\n" <<
	href->seq << "\n" <<
	"+\n" <<
	href->qual << "\n";

    return 0;
}

int 
write_fastq(gzFile *fh, Seq *href) {
    stringstream sstr;
    sstr <<
	"@" << href->id << "\n" <<
	href->seq << "\n" <<
	"+\n" <<
	href->qual << "\n";

    gzputs(*fh, sstr.str().c_str());

    return 0;
}

int 
write_fastq(ofstream *fh, Seq *href, string msg) {
    *fh <<
	"@" << href->id << "|" << msg << "\n" <<
	href->seq << "\n" <<
	"+\n" <<
	href->qual << "\n";

    return 0;
}

int 
write_fastq(gzFile *fh, Seq *href, string msg) {
    stringstream sstr;
    sstr <<
	"@" << href->id << "|" << msg << "\n" <<
	href->seq << "\n" <<
	"+\n" <<
	href->qual << "\n";

    gzputs(*fh, sstr.str().c_str());

    return 0;
}

int 
write_fasta(ofstream *fh, Seq *href, string msg) {
    *fh <<
	">" << 
	href->id  << "|" << msg << "\n" <<
	href->seq << "\n";

    return 0;
}

int 
write_fasta(gzFile *fh, Seq *href, string msg) {
    stringstream sstr;
    sstr <<
	">" << 
	href->id  << "|" << msg << "\n" <<
	href->seq << "\n";

    gzputs(*fh, sstr.str().c_str());

    return 0;
}
