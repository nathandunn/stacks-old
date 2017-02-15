#ifndef DNASEQ4_H
#define DNASEQ4_H

#include <cstring>
#include <iostream>
#include <vector>
#include <string>
#include <functional>

#include "constants.h"

// Definitions for nucleotides coded on 4 bits.
// These definitions are compatible with those of htslib--htslib supports
// partially ambiguous nucleotides ('R', etc.) but we convert everything to 15
// (i.e. 0xF, 'N').
struct Nt4 {
    static const size_t nbits = 4;

    static const size_t a = 1; //0001
    static const size_t c = 2; //0010
    static const size_t g = 4; //0100
    static const size_t t = 8; //1000
    static const size_t n = 15;//1111

    // Trivial ASCII-like hash table giving the 4-bits value of a nucleotide letter.
    // e.g. ch_to_nt [ (int)'G' ] == 4
    // Adapted from `htslib::seq_nt16_table` (hts.cc).
    static const size_t ch_to_nt[256];

    // Trivial hash table giving the nucleotide letter of a 4-bits value.
    // e.g. nt_to_ch[4] == 'C'
    static const char nt_to_ch[16];
};

// Definitions for nucleotides coded on 2 bits.
struct Nt2 {
    static const size_t nbits = 2;

    static const size_t a = 0;
    static const size_t c = 1;
    static const size_t g = 2;
    static const size_t t = 3;

    static const size_t ch_to_nt[256];
    static const char nt_to_ch[4];

    static const size_t nt4_to_nt[16];
};

// A dinucleotide coded on one byte, where the first nucleotide uses the high bits.
// This is compatible with BAM/HTSLIB.
class DiNuc {
    uchar x_;

public:
    DiNuc() : x_(0) {}
    DiNuc(const DiNuc& other) : x_(other.x_) {}
    DiNuc(uchar x) : x_(x) {}
    DiNuc(size_t u1, size_t u2) : x_(0) {set_first(u1); set_second(u2);}
    DiNuc(char c1, char c2) : DiNuc(Nt4::ch_to_nt[(size_t) c1], Nt4::ch_to_nt[(size_t) c2]) {}
    DiNuc& operator= (const DiNuc& other) {x_ = other.x_; return *this;}

    size_t first() const {return x_ >>4;}
    size_t second() const {return x_ & 15;}

    bool operator== (const DiNuc& other) const {return x_ == other.x_;}
    bool operator<  (const DiNuc& other) const {return x_ < other.x_;}

private:
    void set_first(size_t u) {x_ |= u <<4;}
    void set_second(size_t u) {x_ |= u;}

    void clear_first() {x_ &= 15;}
    void clear_second() {x_ &= ~15;}

public:
    uchar x() const {return x_;} // c.f. hash<DNASeq4>
};

// A sequence of 4-bits A, C, G, T and N.
// Uses DiNuc and is thus compatible with BAM/HTSLIB.
class DNASeq4 {
    size_t l_;
    std::vector<DiNuc> v_;

public:
    DNASeq4() : l_(0), v_() {}
    DNASeq4(const DNASeq4& other) : l_(other.l_), v_(other.v_) {}
    DNASeq4(DNASeq4&& other) : l_(other.l_), v_(std::move(other.v_)) {}
    DNASeq4(const char* s, size_t len);
    DNASeq4(const uchar* arr, size_t len) : l_(len), v_(arr, arr+len/2+len%2) {} // `len` is the length of the sequence (not of the array)
    DNASeq4(const std::string& s) : DNASeq4(s.c_str(), s.size()) {}
    DNASeq4& operator= (const DNASeq4& other) {l_ = other.l_; v_ = other.v_; return *this;}
    DNASeq4& operator= (DNASeq4&& other) {l_ = other.l_; v_ = std::move(other.v_); return *this;}

    size_t length() const {return l_;}

    uchar operator[] (size_t i) const {return i%2==0 ? v_[i/2].first() : v_[i/2].second();}
    bool  operator== (const DNASeq4& other) const {return l_ == other.l_ && v_ == other.v_;}
    bool  operator<  (const DNASeq4& other) const {return l_ < other.l_ ? true : v_ < other.v_;}
    friend class std::hash<DNASeq4>;

    // Iterator.
    class const_iterator {
        std::vector<DiNuc>::const_iterator vi_;
        bool first_;

    public:
        const_iterator(std::vector<DiNuc>::const_iterator vi, bool f) : vi_(vi), first_(f) {}
        bool operator!= (const_iterator other) {return ! (vi_ == other.vi_? first_ == other.first_ : false);}
        const_iterator& operator++ () {if (first_) {first_ = false;} else {++vi_; first_ = true;} return *this; }

        // Get the nucleotide.
        uchar operator* () {return first_ ? vi_->first() : vi_->second();}
    };
    const_iterator begin() const {return const_iterator(v_.begin(), true);}
    const_iterator end()   const {return length()%2==0 ? const_iterator(v_.end(), true) : const_iterator(--v_.end(), false);}

    // Methods to allow to memcpy into a htslib bam1_t.
    size_t nbytes() const {return v_.size() * sizeof (DiNuc);}
    const uchar* vdata() const {return (uchar*) v_.data();}
};

namespace std {
template<>
struct hash<DNASeq4> {
    size_t operator() (const DNASeq4& s) const {
        size_t x = static_cast<size_t>(14695981039346656037ULL);
        for (const DiNuc& d : s.v_) {
            x ^= static_cast<size_t>(d.x());
            x *= static_cast<size_t>(1099511628211ULL);
        }
        return x;
    }
};
}

#endif //DNASEQ4_h
