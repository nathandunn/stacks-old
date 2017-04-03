#ifndef DNASEQ4_H
#define DNASEQ4_H

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
    static const vector<size_t> all; // All of the above.

    static size_t from(char c) {return ch_to_nt4[size_t(c)];}
    static char to_ch(size_t nt4) {return nt4_to_ch[nt4];}
    static size_t rev_compl(size_t nt4) {return rev_compl_nt4[nt4];}

private:
    // Trivial ASCII-like hash table giving the 4-bits value of a nucleotide letter.
    // e.g. ch_to_nt [ (int)'G' ] == 4
    // Adapted from `htslib::seq_nt16_table` (hts.cc).
    static const size_t ch_to_nt4[256];

    // Trivial hash table giving the nucleotide letter of a 4-bits value.
    // e.g. nt_to_ch[4] == 'C'
    static const char nt4_to_ch[16];

    // Table giving the reverse complement of the nucleotide.
    static const size_t rev_compl_nt4[16];
};

// Definitions for nucleotides coded on 2 bits.
struct Nt2 {
    static const size_t nbits = 2;

    static const size_t a = 0;
    static const size_t c = 1;
    static const size_t g = 2;
    static const size_t t = 3;
    static const vector<size_t> all;

    static size_t from(char c) {return ch_to_nt2[size_t(c)];}
    static size_t from_nt4(size_t nt4) {return nt4_to_nt2[nt4];}
    static char to_ch(size_t nt2) {return nt2_to_ch[nt2];}
    static size_t rev_compl(size_t nt2) {return rev_compl_nt2[nt2];}

private:
    static const size_t ch_to_nt2[256];
    static const size_t nt4_to_nt2[16];
    static const char nt2_to_ch[4];

    static const size_t rev_compl_nt2[4];
};

class Nt4Counts {
    // Array of counts, containing the count of A's at index Nt4::a, of C's at
    // Nt4::c, G's at Nt4::g, T's at Nt4::t and N's at Nt4::n.
    size_t counts_[16];
    size_t* sorted_[4]; // Pointers to the counts of A, C, G, and T above.

public:
    Nt4Counts()
        : sorted_{counts_+Nt4::a, counts_+Nt4::c, counts_+Nt4::g, counts_+Nt4::t}
        {memset(counts_, 0xFF, 16 * sizeof(size_t)); reset();}

    void reset() {for (size_t nt4 : Nt4::all) counts_[nt4]=0;}
    void increment(size_t nt4) {++counts_[nt4];}
    void sort();

    size_t count(size_t nt4) const {return counts_[nt4];}
    const size_t* rank1() const {return sorted_[3];}
    const size_t* rank2() const {return sorted_[2];}
    const size_t* rank3() const {return sorted_[1];}
    const size_t* rank4() const {return sorted_[0];}
    size_t nt4_of(const size_t* count_ptr) const {return count_ptr-counts_;} // Returns e.g. Nt4::a.

    friend ostream& operator<< (ostream& os, const Nt4Counts& loc);
};

// A sequence of nucleotides on a uint64_t, where the first nucleotide uses the
// low bits.
template<class Nt>
class NtArray {
    uint64_t a_;

public:
    NtArray() : a_(0) {}
    NtArray(int) : a_(-1) {}
    NtArray(const NtArray<Nt>& other) : a_(other.a_) {}
    NtArray<Nt>& operator= (const NtArray<Nt>& other) {a_ = other.a_; return *this;}

    void set(size_t i, size_t nt) {a_ |= nt << (i*Nt::nbits);}
    void clear(size_t i) {a_ &= ~(lowbits << (i*Nt::nbits));}

    size_t operator[] (size_t i) const {return (a_ >> (i*Nt::nbits)) & lowbits;}
    bool operator== (const NtArray<Nt>& other) const {return a_ == other.a_;}
    bool operator<  (const NtArray<Nt>& other) const {return a_ < other.a_;}
    friend class std::hash<NtArray>;

    // Methods for kmers
    void push_front(size_t nt) {a_ <<= Nt::nbits; a_ |= nt;}
    void pop_front() {a_ >>= Nt::nbits;}

    static const size_t n_nts = sizeof a_ * 8 / Nt::nbits;

private:
    static const uint64_t lowbits = (1 << Nt::nbits) - 1;
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
    DiNuc(char c1, char c2) : DiNuc(Nt4::from(c1), Nt4::from(c2)) {}
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
    vector<DiNuc> v_;

public:
    DNASeq4() : l_(0), v_() {}
    DNASeq4(const DNASeq4& other) : l_(other.l_), v_(other.v_) {}
    DNASeq4(DNASeq4&& other) : l_(other.l_), v_(move(other.v_)) {}
    DNASeq4(const char* s, size_t len);
    DNASeq4(const uchar* arr, size_t len) : l_(len), v_(arr, arr+len/2+len%2) {} // `len` is the length of the sequence (not of the array)
    DNASeq4(const string& s) : DNASeq4(s.c_str(), s.size()) {}
    DNASeq4& operator= (const DNASeq4& other) {l_ = other.l_; v_ = other.v_; return *this;}
    DNASeq4& operator= (DNASeq4&& other) {l_ = other.l_; v_ = move(other.v_); return *this;}

    size_t length() const {return l_;}
    string str() const;
    DNASeq4 rev_compl() const;

    size_t operator[] (size_t i) const {return i%2==0 ? v_[i/2].first() : v_[i/2].second();}
    bool  operator== (const DNASeq4& other) const {return l_ == other.l_ && v_ == other.v_;}
    bool  operator<  (const DNASeq4& other) const {return l_ < other.l_ ? true : v_ < other.v_;}
    friend class std::hash<DNASeq4>;

    // Iterator.
    class iterator {
        vector<DiNuc>::const_iterator vi_;
        bool first_;

    public:
        iterator(vector<DiNuc>::const_iterator vi, bool f) : vi_(vi), first_(f) {}
        bool operator!= (iterator other) const {return ! (vi_ == other.vi_? first_ == other.first_ : false);}
        iterator& operator++ () {if (first_) {first_ = false;} else {++vi_; first_ = true;} return *this; }
        iterator& operator-- () {if (first_) {--vi_; first_ = false;} else {first_ = true;} return *this; }

        // Get the (Nt4) nucleotide.
        size_t operator* () const {return first_ ? vi_->first() : vi_->second();}
    };
    iterator begin() const {return iterator(v_.begin(), true);}
    iterator end()   const {return length()%2==0 ? iterator(v_.end(), true) : iterator(--v_.end(), false);}

    // Methods to allow to memcpy into a htslib bam1_t.
    size_t nbytes() const {return v_.size() * sizeof (DiNuc);}
    const uchar* vdata() const {return (uchar*) v_.data();}
};

inline
void Nt4Counts::sort() {
    std::sort(
        sorted_, sorted_+4,
        [] (const size_t* cnt1, const size_t* cnt2) {
            // Primarily on the count.
            // Secondarily on the pointer (i.e. array index/Nt4 value).
            return std::tie(*cnt1, cnt1) < std::tie(*cnt2, cnt2);
        }
    );
}

inline
ostream& operator<< (ostream& os, const Nt4Counts& cnts) {
    const size_t* r1 = cnts.rank1();
    const size_t* r2 = cnts.rank2();
    const size_t* r3 = cnts.rank3();
    const size_t* r4 = cnts.rank4();
    os << Nt4::to_ch(cnts.nt4_of(r1)) << ":" << *r1 << " "
       << Nt4::to_ch(cnts.nt4_of(r2)) << ":" << *r2 << " "
       << Nt4::to_ch(cnts.nt4_of(r3)) << ":" << *r4 << " "
       << Nt4::to_ch(cnts.nt4_of(r4)) << ":" << *r4;
    return os;
}

namespace std { template<class Nt>
struct hash<NtArray<Nt>> { size_t operator() (const NtArray<Nt>& a) const {
    return hash<uint64_t>()(a.a_);
}};}

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
