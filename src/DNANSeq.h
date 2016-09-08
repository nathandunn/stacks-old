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

#include <cstring>
#include <climits>

#define BITMASK(b)     (1 << ((b) % CHAR_BIT))
#define BITSLOT(b)     ((b) / CHAR_BIT)
#define BITSET(a, b)   ((a)[BITSLOT(b)] |= BITMASK(b))
#define BITCLEAR(a, b) ((a)[BITSLOT(b)] &= ~BITMASK(b))
#define BITTEST(a, b)  ((a)[BITSLOT(b)] & BITMASK(b))
#define BITNSLOTS(nb)  ((nb + CHAR_BIT - 1) / CHAR_BIT)

//
// DNA Sequence Storage Class
//
// Three-bit compression, 2.667 bases per byte of storage:
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
    DNANSeq(const char* s) : DNANSeq(strlen(s), s) {}
    DNANSeq(int, unsigned char *);
    DNANSeq(const DNANSeq&);
    DNANSeq& operator= (const DNANSeq& other) =delete;
    ~DNANSeq();

    char  operator[](int);
    int   size() const;
    char *seq(char *);
    char *seq();
    char *subseq(char *, int, int);

    void extend(int before, int after);

    bool operator== (const DNANSeq& other) const {
        unsigned int bytes = BITNSLOTS(bits);
        for (unsigned int i = 0; i < bytes; i++)
            if (s[i] != other.s[i])
                return false;
        return true;
    }

};

#include <iostream>
#include <fstream>
#include <sstream>
using std::stringstream;
using std::cin;
using std::cout;
using std::cerr;

// Specializations for std::hash
// Based on GCC
namespace std {
template<>
struct hash<DNANSeq> {
    size_t operator()(const DNANSeq& seq) const {
        size_t __result = static_cast<size_t>(14695981039346656037ULL);
        unsigned short int __bytes  = BITNSLOTS(seq.bits);
        for (unsigned short int i = 0; i < __bytes; i++) {
            __result ^= static_cast<size_t>(seq.s[i]);
            __result *= static_cast<size_t>(1099511628211ULL);
        }
        return __result;
    }
};

template<>
struct hash<const DNANSeq*> {
    size_t operator()(const DNANSeq* seq) const {return hash<DNANSeq>()(*seq);}
};
}

// namespace __gnu_cxx {
//     template<>
//     struct hash<DNANSeq *>
//     {
//      size_t
//      operator()(DNANSeq *__s) const {
//          unsigned long   __h = 0;
//          unsigned int  bytes = BITNSLOTS(__s->bits);
//          for (unsigned int i = 0; i < bytes; i++)
//              __h = 5 * __h + __s->s[i];
//          return size_t(__h);
//      }
//     };
// }

/* struct hash_dnanseq {
    // Deprecated
    size_t operator()(const DNANSeq* seq) const {return std::hash<DNANSeq>{}(*seq);}
};

struct dnanseq_eqstr {
    // Deprecated
    bool operator()(const DNANSeq *s1, const DNANSeq *s2) const {return *s1 == *s2;}
}; */

#endif // __DNANSeq_H__
