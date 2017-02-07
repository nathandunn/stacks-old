#ifndef DNASEQ4_H
#define DNASEQ4_H

#include <cstring>
#include <iostream>
#include <vector>
#include <string>
#include <functional>

#include "constants.h"

class DiNuc;
class DNASeq4;

namespace nt4 {

// Definitions for nucleotides coded on 4 bits.

// These definitions are compatible with those of htslib--htslib supports
// partially ambiguous nucleotides ('R', etc.) but we convert everything to 15
// (i.e. 0xF, 'N').

// Trivial ASCII-like hash table giving the 4-bits value of a nucleotide letter.
// Adapted from `htslib::seq_nt16_table` (hts.cc).
extern const uchar c2u[256];

extern const uchar a;
extern const uchar c;
extern const uchar g;
extern const uchar t;
extern const uchar n;

// Trivial hash table giving the nucleotide letter of a 4-bits value.
extern const char u2c[16];

} //namespace nt4

// A dinucleotide, coded on one byte.
// The first nucleotide uses the high bits.
class DiNuc {
    uchar x_;

public:
    DiNuc() : x_(0) {}
    DiNuc(const DiNuc& other) : x_(other.x_) {}
    DiNuc(uchar x) : x_(x) {}
    DiNuc(uchar u1, uchar u2) : x_(0) {set_first(u1); set_second(u2);}
    DiNuc(char c1, char c2) : DiNuc(nt4::c2u[(int) c1], nt4::c2u[(int) c2]) {}
    DiNuc& operator= (const DiNuc& other) {x_ = other.x_; return *this;}

    uchar first() const {return x_ >>4;}
    uchar second() const {return x_ & 15;}

    bool operator== (const DiNuc& other) const {return x_ == other.x_;}
    bool operator<  (const DiNuc& other) const {return x_ < other.x_;}

private:
    void set_first(uchar u) {x_ |= u <<4;}
    void set_second(uchar u) {x_ |= u;}

    void clear_first() {x_ &= 15;}
    void clear_second() {x_ &= ~15;}

    friend class std::hash<DNASeq4>;
};

// A sequence of 4-bits A, C, G, T and N.
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
            x ^= static_cast<size_t>(d.x_);
            x *= static_cast<size_t>(1099511628211ULL);
        }
        return x;
    }
};
}

#endif //DNASEQ4_h
