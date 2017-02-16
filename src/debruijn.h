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
    // Given a sequence iterator, builds the first valid kmer (i.e. skipping Ns).
    // If no kmer can be found, an empty object is returned.
    // After the call `first` points to the nucleotide immediately past the kmer.
    Kmer(size_t km_len, DNASeq4::iterator& first, DNASeq4::iterator past) : a_() {
        // Find a series of km_len good nucleotides.
        DNASeq4::iterator km_start = first;
        size_t n_good = 0;
        while(first != past && n_good != km_len) {
            if (first.nt() == Nt4::n) {
                // start again
                km_start = first;
                n_good = 0;
                ++first;
            } else {
                ++n_good;
                ++first;
            }
        }
        // Build the kmer.
        if (n_good == km_len) {
            for (size_t i=0; i<km_len; ++i) {
                a_.set(i, Nt2::nt4_to_nt[km_start.nt()]);
                ++km_start;
            }
        }
    }

    // Create the predecessor/successor kmer given an edge/nucleotide
    Kmer pred(size_t km_len, size_t nt) const
        {Kmer k (*this); k.a_.clear(km_len-1); k.a_.push_front(nt); return k;}
    Kmer succ(size_t km_len, size_t nt) const
        {Kmer k (*this); k.a_.pop_front(); k.a_.set(km_len-1, nt); return k;}

    bool empty() const {return a_ == NtArray<Nt2>();}

    bool operator==(const Kmer& other) const {return a_ == other.a_;}
    friend class std::hash<Kmer>;
};

namespace std { template<>
struct hash<Kmer> { size_t operator() (const Kmer& km) const {
    return hash<NtArray<Nt2>>()(km.a_);
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
    void create(const CLocReadSet& readset);
};

//
// ==================
// Inline definitions
// ==================
//

void Graph::create(const CLocReadSet& readset) {

    //
    // Count all kmers.
    //
    for (const Read& r : readset.reads()) {

        // Build the first kmer.
        DNASeq4::iterator next_nt = r.seq.begin();
        Kmer km = Kmer(km_len, next_nt, r.seq.end());
        if (km.empty()) {
            cerr << "Oops, no " << km_len << "-mers in " << r.seq.str() << "\n"; //
            continue;
        }

        // Record it.
        ++map[km];

        // Walk the sequence.
        while (next_nt != r.seq.end()) {
            size_t nt4 = next_nt.nt();
            if (nt4 == Nt4::n) {
                km = Kmer(km_len, next_nt, r.seq.end());
                if (km.empty())
                    // Not enough sequence remaining to make another kmer.
                    break;
            } else {
                km = km.succ(km_len, Nt2::nt4_to_nt[nt4]);
            }
            ++map[km];
        }
    }

    // TODO Build the graph
}

#endif
