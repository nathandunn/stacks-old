#ifndef DEBRUIJN_H
#define DEBRUIJN_H

#include <iostream>
#include <vector>
#include <list>
#include <unordered_map>
#include <unordered_set>

#include "constants.h"
#include "DNASeq4.h"
#include "locus.h"

//
// ==========
// Kmer
// ==========
//

class Kmer {
    NtArray<Nt2> a_;

public:
    Kmer() : a_() {}

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

    // The first nucleotide.
    size_t front() const {return a_[0];}
    size_t back(size_t km_len) const {return a_[km_len-1];}

    // Create the predecessor/successor kmer given an edge/nucleotide
    Kmer pred(size_t km_len, size_t nt) const
        {Kmer k (*this); k.a_.clear(km_len-1); k.a_.push_front(nt); return k;}
    Kmer succ(size_t km_len, size_t nt) const
        {Kmer k (*this); k.a_.pop_front(); k.a_.set(km_len-1, nt); return k;}

    bool empty() const {return a_ == NtArray<Nt2>();} //TODO This is true for A homopolymers! Use uint64_t(-1) as the default value, safe for kmers <=31bp
    std::string str(size_t km_len) const {
        std::string s;
        s.reserve(km_len);
        for (size_t i=0; i<km_len; ++i)
            s.push_back(Nt2::nt_to_ch[a_[i]]);
        return s;
    }

    bool operator==(const Kmer& other) const {return a_ == other.a_;}
    friend class std::hash<Kmer>;
};

namespace std { template<>
struct hash<Kmer> { size_t operator() (const Kmer& km) const {
    return hash<NtArray<Nt2>>()(km.a_);
}};}

//
// ==========
// Node
// ==========
//

class SPath;

struct NodeData {
    Kmer km;
    size_t count;

    NodeData() : km(), count(0) {}
    NodeData(const Kmer& k, size_t c) : km(k), count(c) {}
};

class Node {
    NodeData d_;
    Node* pred_[4];
    Node* succ_[4];

public:
    Node() : d_(), pred_(), succ_(), sp_() {}
    Node(const NodeData& d) : d_(d), pred_(), succ_(), sp_() {}
    ~Node() {}

    void set_pred(size_t nt2, Node* n) {pred_[nt2] = n;}
    void set_succ(size_t nt2, Node* n) {succ_[nt2] = n;}

    size_t n_pred() const {return size_t(pred_[0]!=NULL) + size_t(pred_[1]!=NULL) + size_t(pred_[2]!=NULL) + size_t(pred_[3]!=NULL);}
    size_t n_succ() const {return size_t(succ_[0]!=NULL) + size_t(succ_[1]!=NULL) + size_t(succ_[2]!=NULL) + size_t(succ_[3]!=NULL);}
    Node* pred(size_t nt2) {return pred_[nt2];}
    Node* succ(size_t nt2) {return succ_[nt2];}

    Node* first_pred() {Node* s = pred_[0]; for (size_t nt2=1; nt2<4; ++nt2) { if (s != NULL) break; s = pred(nt2);} return s;}
    Node* first_succ() {Node* s = succ_[0]; for (size_t nt2=1; nt2<4; ++nt2) { if (s != NULL) break; s = succ(nt2);} return s;}

    const Kmer& km() const {return d_.km;}
    size_t count() const {return d_.count;}

private:
    SPath* sp_;
    friend class SPath;
};

//
// ==========
// SPath (Simple Path)
// ==========
//
// Each node includes a SPath pointer so as to not have to maintain
// separate edge information. This pointer is handled by the methods in
// SPath, not by those of Node. It is null except for the first and last
// nodes of the path.
//

struct SPathData {
    size_t n_nodes;
    size_t km_cumcount;

    SPathData() : n_nodes(0), km_cumcount(0) {}
};

class SPath {
    Node* first_;
    Node* last_;
    SPathData d_;

public:
    bool being_visited;
    bool was_visited;

    SPath() : first_(), last_(), d_(), being_visited(), was_visited() {}
    SPath(Node* first);
    void update_ptrs() {first_->sp_ = this; last_->sp_ = this;}

    size_t n_pred() const {return n_pred();}
    size_t n_succ() const {return last_->n_succ();}
    SPath* pred(size_t nt2) {Node* n = first_->pred(nt2); return n == NULL ? NULL : n->sp_;}
    SPath* succ(size_t nt2) {Node* n = last_->succ(nt2); return n == NULL ? NULL : n->sp_;}

    SPath* first_pred() {Node* n = first_->first_pred(); return n == NULL ? NULL : n->sp_;}
    SPath* first_succ() {Node* n = last_->first_succ(); return n == NULL ? NULL : n->sp_;}

    size_t n_nodes() const {return d_.n_nodes;}
    size_t km_cumcount() const {return d_.km_cumcount;}
    std::string contig_str(size_t km_len);

    const Node* first() const {return first_;}
};

//
// ==========
// Graph
// ==========
//

struct KmMapValue {
    // Holds a kmer count before nodes are built, and a node index after that
    // (and the kmer count is copied to NodeData::count).
    union {
        size_t count;
        size_t node;
    };
    KmMapValue() : count(0) {}
};

class Graph {
    const size_t km_len_;

    std::unordered_map<Kmer, KmMapValue> map_;
    std::vector<Node> nodes_;
    std::vector<SPath> simple_paths_;
    std::vector<SPath*> sorted_spaths_; // The simple paths, sorted topologically, with the terminal (no successors) ones first.

public:
    Graph(size_t km_length) : km_len_(km_length) {}
    void rebuild(const CLocReadSet& readset, size_t min_kmer_count);

    size_t n_simple_paths() const {return simple_paths_.size();}
    bool topo_sort();

    void dump_gfa(const std::string& path);

private:
    // Resets the object.
    void clear() {nodes_.resize(0); map_.clear(); simple_paths_.resize(0); sorted_spaths_.resize(0);}
    size_t index_of(const Node* n) const {return n - nodes_.data();}

    bool topo_sort(SPath* p);
};

#endif
