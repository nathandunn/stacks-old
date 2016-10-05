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

#include <iostream>

#include "DNANSeq.h"

using namespace std;

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

DNANSeq::DNANSeq(const DNANSeq& other) : bits(other.bits) {
    const int n_bytes = BITNSLOTS(bits);
    s = new unsigned char[n_bytes];
    memcpy(s, other.s, n_bytes);
}

DNANSeq::~DNANSeq() {
    delete [] this->s;
}

char DNANSeq::operator[](int pos) const {
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

int DNANSeq::size() const {
    return this->bits / bits_per_nuc;
}

void DNANSeq::seq(char* str) const {
    for (int i = 0; i < size(); i++)
        str[i] = (*this)[i];
    str[size()] = '\0';
}

string DNANSeq::seq() const {
    string str;
    str.reserve(size());
    for (int i = 0; i < size(); i++)
        str.push_back((*this)[i]);
    return str;
}

void DNANSeq::extend(int before, int after) {
    if (before == 0 && after == 0)
        return;

    int old_bits = bits;
    unsigned char* old_s = s;

    bits += (before + after) * bits_per_nuc;
    s = new unsigned char[BITNSLOTS(bits)];
    memset(s, 0, BITNSLOTS(bits));

    int bit = 0;

    // Prepend N's
    for (int i = 0; i < before; ++i) {
        BITSET(s, bit);
        bit += 3;
    }

    // Copy the old sequence.
    for (int old_bit=0; old_bit < old_bits; ++old_bit) {
        if (BITTEST(old_s, old_bit))
            BITSET(s, bit);
        ++bit;
    }

    // Append N's.
    for (int i = 0; i < after; ++i) {
        BITSET(s, bit);
        bit += 3;
    }

    delete[] old_s;
}
