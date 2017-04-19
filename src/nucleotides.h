#ifndef NUCLEOTIDES_H
#define NUCLEOTIDES_H

#include "constants.h"

class Nt2;

//
// Nt4: A nucleotide coded on 4 bits.
// These definitions are compatible with those of htslib--htslib supports
// partially ambiguous nucleotides ('R', etc.) but we convert everything to 15
// (i.e. 0xF, 'N').
//
class Nt4 {
    size_t nt_;

public:
    Nt4() : nt_(Nt4::n) {}
    Nt4(char c) : nt_(from_ch[size_t(c)]) {}
    Nt4(size_t i) : nt_(i) {}
    Nt4(int i) : nt_(i) {}
    Nt4(const Nt4& other) : nt_(other.nt_) {}
    Nt4(const Nt2 nt2);
    Nt4& operator=(const Nt4& other) {nt_ = other.nt_; return *this;}

    Nt4 rev_compl() const {return rev_compl_[size_t(nt_)];}
    size_t index() const {return to_index[size_t(nt_)];}
    explicit operator size_t () const {return nt_;}
    explicit operator int () const {return nt_;}
    explicit operator char () const {return to_ch[nt_];}
    bool operator== (Nt4 other) const {return nt_ == other.nt_;}
    bool operator< (Nt4 other) const {return nt_ < other.nt_;}

    static constexpr size_t nbits = 4;
    static const Nt4 $; // 0000 (0)
    static const Nt4 a; // 0001 (1)
    static const Nt4 c; // 0010 (2)
    static const Nt4 g; // 0100 (4)
    static const Nt4 t; // 1000 (8)
    static const Nt4 n; // 1111 (15)
    static const array<Nt4,5> all; // All of the above.

    static constexpr size_t max() {return (1 << nbits) - 1;}

private:
    // Trivial ASCII-like hash table giving the 4-bits value of a nucleotide letter.
    // e.g. ch_to_nt [ (int)'G' ] == 4
    // Adapted from `htslib::seq_nt16_table` (hts.cc).
    static const Nt4 from_ch[256];

    // Trivial hash table to convert Nt2 nucleotides to Nt4.
    static const Nt4 from_nt2[4];

    // Trivial hash table giving the nucleotide letter of a 4-bits value.
    // e.g. nt_to_ch[4] == 'C'
    static const char to_ch[16];

    // Table giving the reverse complement of the nucleotide.
    static const Nt4 rev_compl_[16];

    // Convert $,a,c,g,t,n into a 0-5 index for arrays of objects indexed by possible Nt4 characters.
    static const size_t to_index[16];
};

//
// Nt2: A nucleotide coded on 2 bits.
//
class Nt2 {
    size_t nt_;

public:
    Nt2() : nt_(Nt2::a) {}
    Nt2(char c) : nt_(from_ch[size_t(c)]) {}
    Nt2(Nt4 nt4) : nt_(from_nt4[size_t(nt4)]) {assert(!(nt4==Nt4::n));}
    Nt2(size_t i) : nt_(i) {}
    Nt2(int i) : nt_(i) {}
    Nt2(const Nt2& other) : nt_(other.nt_) {}
    Nt2& operator=(const Nt2& other) {nt_ = other.nt_; return *this;}

    Nt2 rev_compl() const {return rev_compl_[nt_];}

    explicit operator size_t () const {return nt_;}
    explicit operator int () const {return nt_;}
    explicit operator char () const {return to_ch[nt_];}
    bool operator== (Nt2 other) const {return nt_ == other.nt_;}
    bool operator< (Nt2 other) const {return nt_ < other.nt_;}

    static constexpr size_t nbits = 2;
    static const Nt2 a;
    static const Nt2 c;
    static const Nt2 g;
    static const Nt2 t;
    static const array<Nt2,4> all;

    static constexpr size_t max() {return (1 << nbits) - 1;}

private:
    static const Nt2 from_ch[256];
    static const Nt2 from_nt4[16];
    static const char to_ch[4];

    static const Nt2 rev_compl_[4];
};


//
// Inline definitions
// ==========
//

inline // (Note: Trivial, but this can't be defined in-class.)
Nt4::Nt4(const Nt2 nt2) : nt_(from_nt2[size_t(nt2)]) {}

#endif
