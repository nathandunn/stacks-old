#ifndef DEBRUIJN_H
#define DEBRUIJN_H

#include <vector>
#include <array>
#include <unordered_map>

#include "constants.h"

class Kmer {
    uint64_t kmer;

    bool operator==(const Kmer& other) const {return kmer == other.kmer;}
    friend class std::hash<Kmer>;
};

class Node {
    // NodeInfo info; // sequence of the kmer, etc.

    std::array<size_t, 4> pred; // predecessors, packed to the left
    std::array<size_t, 4> succ; // successors, packed to the left

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
    size_t kmlen;
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

namespace std { template<>
class std::hash<Kmer> { size_t operator() (const Kmer& km) const {
    return std::hash(km.kmer);
}};}

#endif
