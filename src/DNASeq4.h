#ifndef DNASEQ4_H
#define DNASEQ4_H

#include "constants.h"

class Nt2;

// Definitions for nucleotides coded on 4 bits.
// These definitions are compatible with those of htslib--htslib supports
// partially ambiguous nucleotides ('R', etc.) but we convert everything to 15
// (i.e. 0xF, 'N').
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

    explicit operator size_t () const {return nt_;}
    explicit operator char () const {return to_ch[nt_];}
    bool operator== (Nt4 other) const {return nt_ == other.nt_;}
    bool operator< (Nt4 other) const {return nt_ < other.nt_;}

    static constexpr size_t nbits = 4;
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
};

// Definitions for nucleotides coded on 2 bits.
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

template<typename Nt>
class Counts {
    // Array of counts, containing the count of A's at index Nt::a,of C's at
    // index Nt::c, etc.
    array<size_t,Nt::max()+1> counts_;

public:
    Counts() {
        for (size_t& c : counts_)
            c=-1;
        for (Nt nt : Nt::all)
            counts_[size_t(nt)] = 0;
    }
    Counts(const Counts<Nt4>& nt4counts);
    Counts(const Counts<Nt2>& nt2counts);

    void clear() {for (Nt nt : Nt::all) counts_[size_t(nt)]=0;}
    void increment(Nt nt) {++counts_[size_t(nt)];}

    // Get the count for a given nucleotide.
    size_t operator[] (Nt nt) const {return counts_[size_t(nt)];}

    size_t sum() const {return (*this)[Nt::a] + (*this)[Nt::c] + (*this)[Nt::g] + (*this)[Nt::t];}
    array<pair<size_t,Nt>,4> sorted();

    Counts& operator+= (const Counts& other)
        {for (Nt nt : Nt::all) counts_[size_t(nt)] += other.counts_[size_t(nt)]; return *this;}

    // Print the counts.
    friend ostream& operator<< (ostream& os, const Counts& cnts)
        {for (Nt nt : Nt::all) os << char(nt) << ":" << cnts[nt] << " "; return os;}
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

    void set(size_t i, Nt nt) {a_ |= uint64_t(size_t(nt)) << (i*Nt::nbits);}
    void clear(size_t i) {a_ &= ~(lowbits << (i*Nt::nbits));}

    Nt operator[] (size_t i) const {return Nt(size_t((a_ >> (i*Nt::nbits)) & lowbits));}
    bool operator== (const NtArray<Nt>& other) const {return a_ == other.a_;}
    bool operator<  (const NtArray<Nt>& other) const {return a_ < other.a_;}
    friend class std::hash<NtArray>;

    // Methods for kmers
    void push_front(Nt nt) {a_ <<= Nt::nbits; a_ |= uint64_t(size_t(nt));}
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
    DiNuc(Nt4 nt1, Nt4 nt2) : x_(0) {set_first(nt1); set_second(nt2);}
    DiNuc(char nt1, char nt2) : DiNuc(Nt4(nt1), Nt4(nt2)) {}
    DiNuc& operator= (const DiNuc& other) {x_ = other.x_; return *this;}

    Nt4 first() const {return Nt4(size_t(x_ >>4));}
    Nt4 second() const {return Nt4(size_t(x_ & 15));}
    void first(Nt4 nt) {clear_first(); set_first(nt);}
    void second(Nt4 nt) {clear_second(); set_second(nt);}

    bool operator== (const DiNuc& other) const {return x_ == other.x_;}
    bool operator<  (const DiNuc& other) const {return x_ < other.x_;}

private:
    void set_first(Nt4 nt) {x_ |= uchar(size_t(nt)) <<4;}
    void set_second(Nt4 nt) {x_ |= uchar(size_t(nt));}

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
    DNASeq4(DNASeq4&& other) : l_(other.l_), v_(move(other.v_)) {other.clear();}
    DNASeq4(const char* s, size_t len);
    DNASeq4(const uchar* arr, size_t len) : l_(len), v_(arr, arr+len/2+len%2) {} // `len` is the length of the sequence (not of the array)
    DNASeq4(const string& s) : DNASeq4(s.c_str(), s.size()) {}
    DNASeq4& operator= (const DNASeq4& other) {l_ = other.l_; v_ = other.v_; return *this;}
    DNASeq4& operator= (DNASeq4&& other) {l_ = other.l_; v_ = move(other.v_); other.clear(); return *this;}

    size_t length() const {return l_;}
    string str() const;
    DNASeq4 rev_compl() const;
    void set(size_t i, Nt4 nt) {i%2==0 ? v_[i/2].first(nt) : v_[i/2].second(nt);}
    void clear() {l_ = 0; v_ = vector<DiNuc>();}
    void reserve(size_t len) {v_.reserve(len/2+len%2);}
    void append(const DNASeq4& other);

    Nt4 operator[] (size_t i) const {return i%2==0 ? v_[i/2].first() : v_[i/2].second();}
    bool  operator== (const DNASeq4& other) const {return l_ == other.l_ && v_ == other.v_;}
    bool  operator<  (const DNASeq4& other) const {return l_ < other.l_ ? true : v_ < other.v_;}
    friend class std::hash<DNASeq4>;
    friend ostream& operator<< (ostream& os, const DNASeq4& seq) {for (Nt4 nt : seq) os << char(nt); return os;}

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
        Nt4 operator* () const {return first_ ? vi_->first() : vi_->second();}
    };
    iterator begin() const {return iterator(v_.begin(), true);}
    iterator end()   const {return length()%2==0 ? iterator(v_.end(), true) : iterator(--v_.end(), false);}

    // Methods to allow to memcpy into a htslib bam1_t.
    size_t nbytes() const {return v_.size() * sizeof (DiNuc);}
    const uchar* vdata() const {return (uchar*) v_.data();}
};

//
// Inline definitions.
// ==========
//

inline
Nt4::Nt4(const Nt2 nt2) : nt_(from_nt2[size_t(nt2)]) {}

template<> inline
Counts<Nt2>::Counts(const Counts<Nt2>& nt2counts) : counts_(nt2counts.counts_) {}
template<> inline
Counts<Nt2>::Counts(const Counts<Nt4>& nt4counts) : Counts() {
    for (Nt2 nt2 : Nt2::all) // We thus ignore Nt4::n.
        counts_[size_t(nt2)] = nt4counts[Nt4(nt2)];
}
template<> inline
Counts<Nt4>::Counts(const Counts<Nt4>& nt4counts) : counts_(nt4counts.counts_) {}
template<> inline
Counts<Nt4>::Counts(const Counts<Nt2>& nt2counts) : Counts() {
    for (Nt2 nt2 : Nt2::all)
        counts_[size_t(Nt4(nt2))] = nt2counts[nt2];
}

template<typename Nt>
array<pair<size_t,Nt>,4> Counts<Nt>::sorted() {
    array<pair<size_t,Nt>,4> arr {{
        {(*this)[Nt::t], Nt::t},
        {(*this)[Nt::g], Nt::g},
        {(*this)[Nt::c], Nt::c},
        {(*this)[Nt::a], Nt::a}
    }};
    // Sort by decreasing value. Primarily on the count, secondarily on
    // the Nt4 value.
    std::sort(arr.rbegin(), arr.rend());
    return arr;
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
