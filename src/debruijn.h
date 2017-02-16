#ifndef DEBRUIJN_H
#define DEBRUIJN_H

#include <vector>
#include <list>
#include <unordered_map>
#include <unordered_set>

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
    Node(const NodeData& d) : d_(d), pred_(), succ_(), sp_last_(), sp_first_() {}
    void set_pred(size_t nt2, Node* n) {pred_[nt2] = n;}
    void set_succ(size_t nt2, Node* n) {succ_[nt2] = n;}

    size_t n_pred() const {return size_t(pred_[0]!=NULL) + size_t(pred_[1]!=NULL) + size_t(pred_[2]!=NULL) + size_t(pred_[3]!=NULL);}
    size_t n_succ() const {return size_t(succ_[0]!=NULL) + size_t(succ_[1]!=NULL) + size_t(succ_[2]!=NULL) + size_t(succ_[3]!=NULL);}
    Node* pred(size_t nt2) {return pred_[nt2];}
    Node* succ(size_t nt2) {return succ_[nt2];}

    Node* first_succ() {Node* s = succ_[0]; for (size_t nt2=1; nt2<4; ++nt2) { if (s != NULL) break; s = succ(nt2);} return s;}

    const Kmer& km() const {return d_.km;}
    size_t count() const {return d_.count;}

    //
    // Code for simple paths ('sp')
    // A simple path is identified by its first node.
    //
private:
    Node* sp_last_;  // Last node of the path. Set for the first path of the simple path.
    Node* sp_first_; // First node of the path. Set for the last node of the simple path.

public:
    void set_sp_last(Node* n) {sp_last_ = n;}
    void set_sp_first(Node* n) {sp_first_ = n;}

    size_t sp_n_pred() const {is_spfirst(this); return n_pred();}
    size_t sp_n_succ() const {is_spfirst(this); return sp_last_->n_succ();}
    Node* sp_pred(size_t nt2) {is_spfirst(this); is_splast(pred(nt2)); return is_set(pred(nt2)->sp_first_);}
    Node* sp_succ(size_t nt2) {is_spfirst(this); is_splast(sp_last_); return is_set(sp_last_->succ(nt2));}

private:
    //xxx debug
    static Node* is_set(Node* n) {assert(n!=NULL); return n;}
    static void is_spfirst(const Node* n) {assert(n->sp_last_!=NULL);}
    static void is_splast(const Node* n) {assert(n->sp_first_!=NULL);}
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
    std::list<Node*> nodes_wo_preds;

public:
    void create(const CLocReadSet& readset, size_t min_kmer_count);

private:
    std::unordered_set<Node*> sp_visited;
    void clear_sp_visited() {sp_visited.erase(sp_visited.begin(), sp_visited.end());}

    // Recursively builds the simple paths, starting at `start`. Updates sp_visited.
    void build_simple_paths(Node* first);
};

//
// ==================
// Inline definitions
// ==================
//

inline
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
                n.set_succ(nt2, km->second.node);
                km->second.node->set_pred(nt2, &n);
            }
        }
    }

    //
    // Record nodes that don't have predecessors.
    //
    for (Node& n : nodes)
        if (n.n_pred() == 0)
            nodes_wo_preds.push_back(&n);

    //
    // Build the simple paths.
    //
    clear_sp_visited();
    for (Node* n : nodes_wo_preds)
        build_simple_paths(n);
}

inline
void Graph::build_simple_paths(Node* first) {
    if (sp_visited.count(first))
        return;
    sp_visited.insert(first);

    Node* n = first;
    Node* s = n->first_succ();
    while (n->n_succ() == 1 && s->n_pred() == 1) {
        // Extend the simple path.
        n = s;
        s = n->first_succ();
    }
    first->set_sp_last(n);
    n->set_sp_first(first);

    // Check why the simple path ended.
    if (s == NULL) {
        // i.e. `n->n_succ() == 0`
        // No successors. End the recursion.

    } else if (n->n_succ() > 1) {
        // Several successors.
        for (size_t nt2=0; nt2<4; ++nt2)
            if (n->succ(nt2) != NULL)
                build_simple_paths(n->succ(nt2));
    } else {
        // i.e. `s->n_pred() > 1` (as `s` at least has one predecessor, `n`)
        // The next node has several predecessors.
        assert(s->n_pred() != 0); //xxx debug
        first->set_sp_last(n);
        n->set_sp_first(first);

        build_simple_paths(s);
    }
}

#endif
