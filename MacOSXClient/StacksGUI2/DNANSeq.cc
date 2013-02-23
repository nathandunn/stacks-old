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
// DNANSeq.cc
//
// Julian Catchen
// jcatchen@uoregon.edu
// University of Oregon
//
// $Id: DNANSeq.cc 2133 2011-06-07 04:07:41Z catchen $
//

#include "DNANSeq.h"

DNANSeq::DNANSeq(int size) {
    int bytes;

    this->bits = size * bits_per_nuc;

    bytes = BITNSLOTS(this->bits);

    this->s = new unsigned char[bytes];
    memset(this->s, 0, bytes);
}

DNANSeq::DNANSeq(int size, unsigned char *seq) {
    unsigned int bytes;

    this->bits = size * bits_per_nuc;

    bytes = BITNSLOTS(this->bits);

    this->s = new unsigned char[bytes];
    for (unsigned int i = 0; i < bytes; i++)
	this->s[i] = seq[i];
}

DNANSeq::DNANSeq(int size, const char *seq) {

    this->bits = size * bits_per_nuc;
    int bytes  = BITNSLOTS(this->bits);
    this->s    = new unsigned char[bytes];

    memset(this->s, 0, bytes);

    int bit = 0;

    for (int i = 0; i < size; i++) {
	switch (seq[i]) {
	case 'A':
	case 'a':
	    // A == 000
	    bit += 3;
	    break;
	case 'C':
	case 'c':
	    // C == 001
	    bit += 2;
	    BITSET(this->s, bit);
	    bit++;
	    break;
	case 'G':
	case 'g':
	    // G == 010
	    bit++;
	    BITSET(this->s, bit);
	    bit++;
	    bit++;
	    break;
	case 'T':
	case 't':
	    // T == 011
	    bit++;
	    BITSET(this->s, bit);
	    bit++;
	    BITSET(this->s, bit);
	    bit++;
	    break;
	case 'N':
	case 'n':
	case '.':
	    // N == 100
	    BITSET(this->s, bit);
	    bit += 3;
	    break;
	}
    }
}

DNANSeq::~DNANSeq() {
    delete [] this->s;
}

char DNANSeq::operator[](int pos) {
    unsigned char c, base;
    int bit;

    if (pos > ((this->bits / bits_per_nuc) - 1)) return '\0';

    bit = pos * bits_per_nuc;

    c    = 0;
    base = 'X';

    for (int i = bits_per_nuc - 1; i >= 0; i--) {
	if (BITTEST(this->s, bit))
	    c |= 1 << i;
	bit++;
    }

    switch (c) {
    case 0:
	base = 'A';
	break;
    case 1:
	base = 'C';
	break;
    case 2:
	base = 'G';
	break;
    case 3:
	base = 'T';
	break;
    case 4:
	base = 'N';
	break;
    default:
	cerr << "Unknown character " << (int) c << "\n";
	break;
    }
    //cerr << "  Decoding character " << pos << ", '" << base << "'\n";

    return base;
}

int DNANSeq::size() {
    return this->bits / bits_per_nuc;
}

char *DNANSeq::subseq(char *seq, int start, int end) {
    int i;

    for (i = start; i <= end; i++)
	seq[i - start] = this->operator[](i);

    seq[i - start] = '\0';

    return seq;
}

char *DNANSeq::seq(char *seq) {
    int i;
    int end = this->bits / bits_per_nuc;

    for (i = 0; i < end; i++)
	seq[i] = this->operator[](i);

    seq[i] = '\0';

    return seq;
}

char *DNANSeq::seq() {
    int i;
    int  size = this->bits / bits_per_nuc;
    char *seq = new char[size + 1];

    for (i = 0; i < size; i++)
	seq[i] = this->operator[](i);

    seq[i] = '\0';

    return seq;
}
