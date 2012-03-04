// -*-mode:c++; c-style:k&r; c-basic-offset:4;-*-
//
// Copyright 2010, Julian Catchen <jcatchen@uoregon.edu>
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

#ifndef __SAM_H__
#define __SAM_H__

//
// Code to parse Sam format. This format is created for
// reads that have been aligned to a reference genome. It takes the tab-separated form:
//
// <query> <strand> <chromosome> <base pair> ... <sequence> <phred quality score> ...
//
// One record per line.
//

#include "input.h"

class Sam: public Input {
    int parse_cigar(int, const char *);

 public:
    Sam(const char *path) : Input(path) {};
    ~Sam() {};
    Seq *next_seq();
    int  next_seq(Seq &) { return 0; };
};

Seq *
Sam::next_seq() 
{
    vector<string> parts;
    int flag;

    //
    // Read a record from the file and place it in a Seq object, skipping header 
    // definitions and unaligned sequences.
    //
    do {
        this->fh.getline(this->line, max_len);

        if (!this->fh.good())
            return NULL;

        parse_tsv(this->line, parts);

	//
	// According to SAM spec FLAGs are the second field, 
	// if FLAG bit 0x4 is set, sequence is not mapped.
	//
	flag = atoi(parts[1].c_str());
	flag = flag  & 4;
	flag = flag >> 2;

    } while (parts[0][0] == '@' || flag == 1);

    //
    // Check which strand this is aligned to: 
    //   SAM reference: FLAG bit 0x10 - sequence is reverse complemented
    //
    flag = atoi(parts[1].c_str());
    flag = flag & 16;
    flag = flag >> 4;

    //
    // If the read was aligned on the reverse strand (and is therefore reverse complemented)
    // alter the start point of the alignment to reflect the right-side of the read, at the
    // end of the RAD cut site.
    //
    // To accomplish this, we must parse the alignment CIGAR string
    //
    int bp = flag ? this->parse_cigar(atoi(parts[3].c_str()), parts[5].c_str()) : atoi(parts[3].c_str());

    //
    // Sam format has a 1-based offset for chrmosome/basepair positions, adjust it to match
    // the Stacks, 0-based offset.
    //
    bp--;

    Seq *s = new Seq(parts[0].c_str(), parts[9].c_str(), parts[10].c_str(), // Read ID, Sequence, Quality
		     parts[2].c_str(), bp, flag ? minus : plus);            // Chr, BasePair, Strand

    return s;
}

int 
Sam::parse_cigar(int aln_bp, const char *cigar)
{
    char buf[id_len];
    int  dist;
    const char *p, *q;

    p = cigar;

    //cerr << "Starting bp: " << aln_bp << "\n";

    while (*p != '\0') {
	q = p + 1;

	while (*q != '\0' && isdigit(*q))
	    q++;
	strncpy(buf, p, q - p);
	buf[q-p+1] = '\0';
	dist = atoi(buf);

	switch(*q) {
	case 'S':
	case 'D':
	    break;
	case 'M':
	case 'I':
	    aln_bp += dist;
	    break;
	}
	p = q + 1;

	//cerr << "  CIGAR: " << cigar << "; Dist: " << dist << "; Type: " << *q << "; aln_bp: " << aln_bp << " [" << *p << "]\n";
    }

    return aln_bp;
}

#endif // __SAM_H__
