#include <cstring>
#include <iostream>
#include <vector>
#include <string>
#include <functional>

#include "hts.h" // For `seq_nt16_table`, `seq_nt16_str`.

#include "constants.h"

class DiNuc;
class DNASeq4;

// A dinucleotide, coded on one byte (4 bits per nucleotide).
// A,C,T,G,N are 1,2,4,8,15 like in Htslib, and the first
// nucleotide is coded with the high bits.
// Unlike in Htslib, '\0' (rather than '=') is 0, and other
// letters are converted to Ns.
class DiNuc {
    uchar u_;

public:
    DiNuc() : u_(0) {}
    DiNuc(const DiNuc& other) : u_(other.u_) {}
    DiNuc(char c1, char c2) : u_(0) {set_first(c2u[(uchar) c1]); set_second(c2u[(uchar) c2]);}

    uchar first() const {return u_ >>4;}
    uchar second() const {return u_ & 15;}
    char cfirst() const {return u2c[first()];}
    char csecond() const {return u2c[second()];}

    DiNuc& operator=(const DiNuc& other) {u_ = other.u_; return *this;}
    bool operator==(const DiNuc& other) const {return u_ == other.u_;}
    bool operator<(const DiNuc& other) const {return u_ < other.u_;}

private:
    void set_first(uchar u) {u_ |= u <<4;}
    void set_second(uchar u) {u_ |= u;}

    void clear_first() {u_ &= 15;}
    void clear_second() {u_ &= ~15;}

    // Trivial ASCII-like hash table for characters, adapted from `htslib::seq_nt16_table` (hts.h).
    // Value is 1 for 'A' and 'a', etc., 0 for '\0' (htslib uses '='), 15 otherwise.
    static const uchar c2u[256];
    static const char u2c[16];

public:
    uchar value() const {return u_;} // for hashing.
};

class DNASeq4 {
    std::vector<DiNuc> v_;

public:
    DNASeq4() : v_() {}
    DNASeq4(const DNASeq4& other) : v_(other.v_) {}
    DNASeq4(DNASeq4&& other) : v_(std::move(other.v_)) {}
    DNASeq4(const char* s, size_t len);
    DNASeq4(const std::string& s) : DNASeq4(s.c_str(), s.size()) {}

    size_t length() const {return v_.size()*2 - size_t(v_.back().second()==0);}
    char nuc(size_t i) const {return i%2==0 ? v_[i].cfirst() : v_[i].csecond();}

    size_t vsize() const {return v_.size();}
    DiNuc operator[](size_t i) const {return v_[i];}
    const DiNuc* vdata() const {return v_.data();}

    DNASeq4& operator=(const DNASeq4& other) {v_ = other.v_; return *this;}
    DNASeq4& operator=(DNASeq4&& other) {v_ = std::move(other.v_); return *this;}
    bool operator==(const DNASeq4& other) const {return v_ == other.v_;}
    bool operator<(const DNASeq4& other) const {return v_ < other.v_;}
};

// Specialization for std::hash
// Based on GCC
namespace std {
template<>
struct hash<DNASeq4> {
    size_t operator()(const DNASeq4& s) const {
        size_t x = static_cast<size_t>(14695981039346656037ULL);
        for (size_t i = 0; i < s.vsize(); i++) {
            x ^= static_cast<size_t>(s[i].value());
            x *= static_cast<size_t>(1099511628211ULL);
        }
        return x;
    }
};
}
