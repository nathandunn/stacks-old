#ifndef DEBRUIJN_H
#define DEBRUIJN_H

#include <vector>
#include <array>
#include <unordered_map>

#include "constants.h"
#include "DNASeq4.h"

class Kmer {
    NtArray<Nt2> a_;

public:
    // Build a kmer from a sequence iterator. (No bounds checking, sequence must
    // be long enough.)
    Kmer(size_t km_len, DNASeq4::iterator si) : a_() {
        for (size_t i=0; i<km_len; ++i) {
            a_.set(i, Nt2::nt4_to_nt[si.nt()]);
            ++si;
        }
    }

    // Create the predecessor/successor kmer given an edge/nucleotide
    Kmer pred(size_t km_len, size_t nt) const
        {Kmer k (*this); k.a_.clear(km_len-1); k.a_.push_front(nt); return k;}
    Kmer succ(size_t km_len, size_t nt) const
        {Kmer k (*this); k.a_.pop_front(); k.a_.set(km_len-1, nt); return k;}

    bool operator==(const Kmer& other) const {return a_ == other.a_;}
    friend class std::hash<Kmer>;
};

namespace std { template<>
struct hash<Kmer> { size_t operator() (const Kmer& km) const {
    return hash()(km.a_);
}};}

class Node {
    // NodeInfo info; // sequence of the kmer, etc.

    std::array<size_t, 4> pred_; // predecessors, packed to the left
    std::array<size_t, 4> succ_; // successors, packed to the left

public:
    size_t n_pred() const;
    size_t n_succ() const;
    size_t pred(size_t i) const;
    size_t succ(size_t i) const;

    //
    // Code for simple paths ('sp')
    // ----------
    //
private:
    size_t sp_first; // First node of the path
    size_t sp_last;  // Last node of the path

    size_t sp_n_pred() const;
    size_t sp_n_succ() const;
    size_t sp_pred(size_t i) const;
    size_t sp_succ(size_t i) const;
};

class Graph {
    size_t km_len;
    std::unordered_map<Kmer, size_t> map;
    std::vector<Node> nodes;

public:
    void create(const CLocReadSet& reads);
};

//
// ==================
// Inline definitions
// ==================
//

#endif
