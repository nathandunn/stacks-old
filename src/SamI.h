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

#ifndef __SAMI_H__
#define __SAMI_H__

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
    int parse_cigar(const char *, vector<pair<char, uint> > &, bool);
    int find_start_bp_neg(int, vector<pair<char, uint> > &);
    int find_start_bp_pos(int, vector<pair<char, uint> > &);
    int edit_gaps(vector<pair<char, uint> > &, char *);

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
    int  flag;
    uint len;

    //
    // Read a record from the file and place it in a Seq object, skipping header 
    // definitions and unaligned sequences.
    //
    do {
        this->fh.getline(this->line, max_len);

        if (!this->fh.good())
            return NULL;

	len = strlen(this->line);
	if (this->line[len - 1] == '\r') this->line[len - 1] = '\0';

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
    vector<pair<char, uint> > cigar;
    this->parse_cigar(parts[5].c_str(), cigar, flag);

    int bp = flag ? 
	this->find_start_bp_neg(atoi(parts[3].c_str()), cigar) : 
	this->find_start_bp_pos(atoi(parts[3].c_str()), cigar);

    //
    // Sam format has a 1-based offset for chrmosome/basepair positions, adjust it to match
    // the Stacks, 0-based offset.
    //
    bp--;

    Seq *s = new Seq(parts[0].c_str(), parts[9].c_str(), parts[10].c_str(), // Read ID, Sequence, Quality
		     parts[2].c_str(), bp, flag ? minus : plus);            // Chr, BasePair, Strand

    if (cigar.size() > 0)
	this->edit_gaps(cigar, s->seq);

    return s;
}

int 
Sam::parse_cigar(const char *cigar_str, vector<pair<char, uint> > &cigar, bool orientation)
{
    char buf[id_len];
    int  dist;
    const char *p, *q;

    p = cigar_str;

    if (*p == '*') return 0;

    while (*p != '\0') {
	q = p + 1;

	while (*q != '\0' && isdigit(*q))
	    q++;
	strncpy(buf, p, q - p);
	buf[q-p] = '\0';
	dist = atoi(buf);

	//
	// If aligned to the negative strand, sequence has been reverse complemented and 
	// CIGAR string should be interpreted in reverse.
	//
	if (orientation == plus)
	    cigar.push_back(make_pair(*q, dist));
	else
	    cigar.insert(cigar.begin(), make_pair(*q, dist));

	p = q + 1;
    }

    return 0;
}

int 
Sam::find_start_bp_neg(int aln_bp, vector<pair<char, uint> > &cigar)
{
    uint size = cigar.size();
    char op;
    uint dist;

    for (uint i = 0; i < size; i++)  {
	op   = cigar[i].first;
	dist = cigar[i].second;

	switch(op) {
	case 'I':
	    break;
	case 'S':
	    if (i < size - 1)
		aln_bp += dist;
	    break;
	case 'M':
	case 'D':
	    aln_bp += dist;
	    break;
	}
    }

    return aln_bp - 1;
}

int 
Sam::find_start_bp_pos(int aln_bp, vector<pair<char, uint> > &cigar)
{
    char op;
    uint dist;

    op   = cigar[0].first;
    dist = cigar[0].second;

    if (op == 'S')
	aln_bp -= dist;

    return aln_bp;
}

int 
Sam::edit_gaps(vector<pair<char, uint> > &cigar, char *seq)
{
    char buf[id_len];
    uint size = cigar.size();
    char op;
    uint dist, bp, len, buf_len, j, k, stop;

    len = strlen(seq);
    bp  = 0;

    for (uint i = 0; i < size; i++)  {
	op   = cigar[i].first;
	dist = cigar[i].second;

	switch(op) {
	case 'S':
	    stop = bp + dist;
	    stop = stop > len ? len : stop;
	    while (bp < stop) {
		seq[bp] = 'N';
		bp++;
	    }
	    break;
	case 'D':
	    //
	    // A deletion has occured in the read relative to the reference genome.
	    // Pad the read with sufficent Ns to match the deletion, shifting the existing
	    // sequence down. Trim the final length to keep the read length consistent.
	    //
	    strncpy(buf, seq + bp, id_len - 1);
	    buf[id_len - 1] = '\0';
	    buf_len         = strlen(buf);

	    stop = bp + dist;
	    stop = stop > len ? len : stop;
	    while (bp < stop) {
		seq[bp] = 'N';
		bp++;
	    }

	    j = bp;
	    k = 0;
	    while (j < len && k < buf_len) {
		seq[j] = buf[k];
		k++;
		j++;
	    }
	    break;
	case 'I':
	    //
	    // An insertion has occurred in the read relative to the reference genome. Delete the
	    // inserted bases and pad the end of the read with Ns.
	    //
	    k = bp + dist;
	    strncpy(buf, seq + k, id_len - 1);
	    buf[id_len - 1] = '\0';
	    buf_len         = strlen(buf);

	    j = bp;
	    k = 0;
	    while (j < len && k < buf_len) {
		seq[j] = buf[k];
		k++;
		j++;
	    }

	    stop = j + dist;
	    stop = stop > len ? len : stop;
	    while (j < stop) {
		seq[j] = 'N';
		j++;
	    }
	    break;
	case 'M':
	    bp += dist;
	    break;
	}
    }

    return 0;
}

#endif // __SAMI_H__
