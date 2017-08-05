#ifndef DEBRUIJN_H
#define DEBRUIJN_H

#include <cstdint>
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
// Up to 31-mers. 32 is possible in principle but `empty()` would be true for a
// series of 32 T's.
//

class Kmer {
    NtArray<Nt2> a_;

public:
    Kmer() : a_(-1) {}

    // Given a sequence iterator, builds the first valid kmer (i.e. skipping Ns).
    // If no kmer can be found, an empty object is returned.
    // After the call `first` points to the nucleotide immediately past the kmer.
    Kmer(size_t km_len, DNASeq4::iterator& first, DNASeq4::iterator past);

    // The first nucleotide.
    Nt2 front() const {return a_[0];}
    Nt2 back(size_t km_len) const {return a_[km_len-1];}

    // Create the predecessor/successor kmer given an edge/nucleotide.
    Kmer pred(size_t km_len, Nt2 nt) const {Kmer k (*this); k.a_.clear(km_len-1); k.a_.push_front(nt); return k;}
    Kmer succ(size_t km_len, Nt2 nt) const {Kmer k (*this); k.a_.pop_front(); k.a_.set(km_len-1, nt); return k;}

    bool empty() const {return *this == Kmer();}
    string str(size_t km_len) const {
        string s;
        s.reserve(km_len);
        for (size_t i=0; i<km_len; ++i)
            s.push_back(char(a_[i]));
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
    uint32_t count;

    NodeData() : km(), count(0) {}
    NodeData(const Kmer& k, uint32_t c) : km(k), count(c) {}
};

class Node {
    NodeData d_;
    Node* pred_[4];
    Node* succ_[4];

public:
    Node() : d_(), pred_(), succ_(), sp_() {}
    Node(const NodeData& d) : d_(d), pred_(), succ_(), sp_() {}
    ~Node() {}

    void set_succ(size_t nt2, Node* n) {succ_[nt2] = n; n->pred_[size_t(d_.km.front())] = this;}

    size_t n_pred() const {return size_t(pred_[0]!=NULL) + size_t(pred_[1]!=NULL) + size_t(pred_[2]!=NULL) + size_t(pred_[3]!=NULL);}
    size_t n_succ() const {return size_t(succ_[0]!=NULL) + size_t(succ_[1]!=NULL) + size_t(succ_[2]!=NULL) + size_t(succ_[3]!=NULL);}
    const Node* pred(size_t nt2) const {return pred_[nt2];}
          Node* pred(size_t nt2)       {return pred_[nt2];} // (non-const)
    const Node* succ(size_t nt2) const {return succ_[nt2];}
          Node* succ(size_t nt2)       {return succ_[nt2];}
    const Node* first_succ() const {const Node* s = succ_[0]; for (size_t nt2=1; nt2<4; ++nt2) { if (s != NULL) break; s = succ(nt2);} return s;}
          Node* first_succ()       {return (Node*) ((const Node*)this)->first_succ();}

    const Kmer& km() const {return d_.km;}
    size_t count() const {return d_.count;}

private:
    mutable SPath* sp_;
    friend class SPath;
};

//
// ==========
// SPath
// ==========
//
// Node-like object for a 'simple path'.
//
// Each node includes a SPath pointer so as to not have to maintain
// separate edge information. This pointer is handled by the methods in
// SPath, not by those of Node. It is null except for the first and last
// nodes of the path.
//
// Includes an externally usable pointer that is used by graph algorithms as
// a 'visited' flag and to insert information relevant to the algorithm into
// the object.
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
    mutable void* visitdata; // For graph algorithms.

    SPath() : first_(), last_(), d_(), visitdata() {}
    SPath(Node* first);
    void update_ptrs() {first_->sp_ = this; last_->sp_ = this;}

    size_t n_pred() const {return first_->n_pred();}
    size_t n_succ() const {return last_->n_succ();}
    const SPath* pred(size_t nt2) const {const Node* n = first_->pred(nt2); return n == NULL ? NULL : n->sp_;}
          SPath* pred(size_t nt2)       {return (SPath*) ((const SPath*)this)->pred(nt2);}
    const SPath* succ(size_t nt2) const {const Node* n = last_->succ(nt2); return n == NULL ? NULL : n->sp_;}
          SPath* succ(size_t nt2)       {return (SPath*) ((const SPath*)this)->succ(nt2);}
    const SPath* first_succ() const {const Node* n = last_->first_succ(); return n == NULL ? NULL : n->sp_;}
          SPath* first_succ()       {return (SPath*) ((const SPath*)this)->first_succ();}

    size_t n_nodes() const {return d_.n_nodes;}
    size_t km_cumcount() const {return d_.km_cumcount;}

    const Node* first() const {return first_;}

    template<typename SPathIt>
    static string contig_str(SPathIt first, SPathIt past, size_t km_len);
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

    unordered_map<Kmer, KmMapValue> map_;
    vector<Node> nodes_;
    vector<SPath> simple_paths_;

    vector<SPath*> sorted_spaths_; // The simple paths, sorted topologically, with the terminal (no successors) ones first.
    map<const void*,set<SPath*>> components_; // The simple paths, by component;
                                              // one simple path address serves as an ID for each component.
    map<const SPath*,const void*> sp_to_component_;

public:
    Graph(size_t km_length) : km_len_(km_length) {}
    void rebuild(const vector<const DNASeq4*>& reads, size_t min_kmer_count);

    size_t empty() const {return simple_paths_.empty();}

    // Find all connected components. This uses an undirected depth first search.
    void compute_components();

    // Finds the best path in the graph.
    // Return false if the graph is not a DAG.
    bool find_best_path(vector<const SPath*>& best_path);

    void dump_gfa(const string& path) const;

private:
    // Clears all members (except km_len_).
    void clear();

    size_t index_of(const Node* n) const {return n - nodes_.data();}

    // Sort topologically. Returns false if the graph is not a DAG.
    // The nodes record whether they are a parent in the current recursion
    // (and we use uchars instead of bools because vector<bool> is specialized).
    bool topo_sort();
    bool topo_sort(SPath* p, vector<uchar>& visitdata);

    // c.f. find_components()
    void propagate_component_id(const SPath* p, void* id);
};

//
// ==========
// Inline definitions
// ==========
//

inline
void Graph::clear() {
    nodes_.resize(0);
    map_.clear();
    simple_paths_.resize(0);
    sorted_spaths_.resize(0);
    components_.clear();
    sp_to_component_.clear();
}

inline
Kmer::Kmer(size_t km_len, DNASeq4::iterator& first, DNASeq4::iterator past) : a_() {
    // Find a series of km_len good nucleotides.
    DNASeq4::iterator km_start = first;
    size_t n_good = 0;
    while(first != past && n_good != km_len) {
        if (*first == Nt4::n) {
            // start again
            ++first;
            km_start = first;
            n_good = 0;
        } else {
            ++n_good;
            ++first;
        }
    }
    // Build the kmer.
    if (n_good == km_len) {
        for (size_t i=0; i<km_len; ++i) {
            a_.set(i, Nt2(*km_start));
            ++km_start;
        }
    }
}

template<typename SPathIt>
string SPath::contig_str(SPathIt first, SPathIt past, size_t km_len) {

    // Compute the final size.
    size_t ctg_len = km_len - 1;
    for (SPathIt sp=first; sp!=past; ++sp)
        ctg_len += (*sp)->n_nodes();

    // Initialize the contig with the contents of the first kmer minus its last
    // nucleotide.
    string ctg = (*first)->first_->km().str(km_len);
    ctg.pop_back();
    ctg.reserve(ctg_len);

    // Extend the contig; loop over every Node of every SPath.
    while (first != past) {
        Node* n = (*first)->first_;
        while (n != (*first)->last_) {
            ctg.push_back(char(n->km().back(km_len)));
            n = n->first_succ();
        }
        ctg.push_back(char(n->km().back(km_len)));
        ++first;
    }

    return ctg;
}

#endif
