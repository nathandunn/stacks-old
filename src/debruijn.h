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

struct NodeData {
    Kmer km;
    size_t count;

    NodeData(const Kmer& k, size_t c) : km(k), count(c) {}
};

class Node {
    NodeData d_;
    Node* pred_[4];
    Node* succ_[4];

public:
    Node(const NodeData& d) : d_(d), pred_(), succ_(), sp_first(), sp_last() {}

    size_t n_pred() const;
    size_t n_succ() const;
    Node& pred(size_t i) const;
    Node& succ(size_t i) const;
    void pred(size_t i, Node* n) {pred_[i] = n;}
    void succ(size_t i, Node* n) {succ_[i] = n;}

    // Data
    const Kmer& km() const {return d_.km;}
    size_t count() const {return d_.count;}

    //
    // Code for simple paths ('sp')
    // ----------
    //
private:
    Node* sp_first; // First node of the path
    Node* sp_last;  // Last node of the path

public:
    size_t sp_n_pred() const;
    size_t sp_n_succ() const;
    Node& sp_pred(size_t i) const;
    Node& sp_succ(size_t i) const;
};

struct KmMapValue {
    union {
        size_t count;
        Node* node;
    };
    KmMapValue() : count(0) {}
};

class Graph {
    size_t km_len;
    std::unordered_map<Kmer, KmMapValue> map;
    std::vector<Node> nodes;

public:
    void create(const CLocReadSet& readset, size_t min_kmer_count);
};

//
// ==================
// Inline definitions
// ==================
//

void Graph::create(const CLocReadSet& readset, size_t min_kmer_count) {

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
        ++map[km].count;

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
            ++map[km].count;
        }
    }

    //
    // Build the standalone nodes.
    //
    for(auto km=map.begin(); km!=map.end();) {
        if (km->second.count < min_kmer_count) {
            map.erase(km++);
        } else {
            nodes.push_back(Node(NodeData(km->first, km->second.count)));
            // Replace the count of the kmer with a pointer to its node.
            km->second.node = &nodes.back();
            ++km;
        }
    }

    //
    // Build the edges.
    //
    for (Node& n : nodes) {
        // Check each possible successor.
        for (size_t nt2=0; nt2<4; ++nt2) {
            auto km = map.find(n.km().succ(km_len, nt2));
            if (km != map.end()) {
                n.succ(nt2, km->second.node);
                km->second.node->pred(nt2, &n);
            }
        }
    }
}

#endif
