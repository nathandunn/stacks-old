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

#ifndef __DNANSeq_H__
#define __DNANSeq_H__

#include <string.h>
#include <limits.h>

#define BITMASK(b)     (1 << ((b) % CHAR_BIT))
#define BITSLOT(b)     ((b) / CHAR_BIT)
#define BITSET(a, b)   ((a)[BITSLOT(b)] |= BITMASK(b))
#define BITCLEAR(a, b) ((a)[BITSLOT(b)] &= ~BITMASK(b))
#define BITTEST(a, b)  ((a)[BITSLOT(b)] & BITMASK(b))
#define BITNSLOTS(nb)  ((nb + CHAR_BIT - 1) / CHAR_BIT)

//
// We expect (and C++ defines) an unsigned char as 8 bits.
//
const unsigned short int bits_per_nuc = 3;
const unsigned short int byte_size    = 8;

//
// DNA Sequence Storage Class
//
// Two-bit compression, four bases per byte of storage:
//    A == 000
//    C == 001
//    G == 010
//    T == 011
//    N == 100
//
class DNANSeq {
public:
    //
    // The number of bits required to store string of DNA string
    //
    unsigned short int bits;
    //
    // Array of bytes to store DNA sequence.
    //
    unsigned char *s;

    DNANSeq(int);
    DNANSeq(int, const char *);
    DNANSeq(int, unsigned char *);
    ~DNANSeq();

    char  operator[](int);
    int   size();
    char *seq(char *);
    char *seq();
    char *subseq(char *, int, int);
};

#include <iostream>
#include <fstream>
#include <sstream>
using std::stringstream;
using std::cin;
using std::cout;
using std::cerr;

// namespace __gnu_cxx {
//     template<>
//     struct hash<DNANSeq *>
//     {
// 	size_t
// 	operator()(DNANSeq *__s) const {
// 	    unsigned long   __h = 0;
// 	    unsigned int  bytes = BITNSLOTS(__s->bits);
// 	    for (unsigned int i = 0; i < bytes; i++)
// 		__h = 5 * __h + __s->s[i];
// 	    return size_t(__h);
// 	}
//     };
// }

struct hash_dnanseq {
    size_t operator()(DNANSeq *__s) const
    {
	size_t __result = static_cast<size_t>(14695981039346656037ULL);
	unsigned short int __bytes  = BITNSLOTS(__s->bits);
	for (unsigned short int i = 0; i < __bytes; i++) {
	    __result ^= static_cast<size_t>(__s->s[i]);
	    __result *= static_cast<size_t>(1099511628211ULL);
	}

	return __result;
    }
};

struct dnanseq_eqstr {
    bool operator()(DNANSeq *s1, DNANSeq *s2) const {
	unsigned int bytes = BITNSLOTS(s1->bits);
	for (unsigned int i = 0; i < bytes; i++)
	    if (s1->s[i] != s2->s[i]) return false;
	return true;
    }
};

#endif // __DNANSeq_H__
