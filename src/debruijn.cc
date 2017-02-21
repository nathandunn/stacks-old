#include <iostream>
#include <iomanip>
#include <sstream>

#include "constants.h"
#include "debruijn.h"

using namespace std;

void Node::sp_build() {
    Node* n = this;
    Node* s = n->first_succ();
    while (n->n_succ() == 1 && s->n_pred() == 1) {
        // Extend the simple path.
        n = s;
        s = n->first_succ();
    }

    // Record the start & end.
    this->sp_last_ = n;
    n->sp_first_ = this;
}

string Node::sp_path_str(size_t km_len) {
    is_spfirst(this);

    // Seeding the sequence with the entire first kmer gives contigs with
    // overlapping ends. However this is not practical for visualization purposes
    // as the lengths of the sequences do not correspond to the number of nodes;
    // instead the lengths are `n_nodes + km_len - 1`.
    // Thus here we use the last nucleotide of the node/kmer to represent the
    // graph in GFA.
    //string s = d_.km.str(km_len); // e.g. same as Minia
    string s;

    Node* n=this;
    while (n != sp_last_) {
        s.push_back(Nt2::nt_to_ch[n->d_.km.back(km_len)]);
        n = n->first_succ();
    }
    s.push_back(Nt2::nt_to_ch[n->d_.km.back(km_len)]);

    return s;
}

size_t Node::sp_cum_count() {
    is_spfirst(this);

    size_t cumcount = 0;
    Node* n=this;
    while (n != sp_last_) {
        cumcount += n->d_.count;
        n = n->first_succ();
    }
    // n == sp_last_
    cumcount += n->d_.count;

    return cumcount;
}
double Node::sp_mean_count() {
    is_spfirst(this);

    size_t cumcount = 0;
    size_t n_nodes = 0;

    Node* n=this;
    while (n != sp_last_) {
        ++n_nodes;
        cumcount += n->d_.count;
        n = n->first_succ();
    }
    // n == sp_last_
    ++n_nodes;
    cumcount += n->d_.count;

    return (double) cumcount / n_nodes;
}

void Graph::create(const CLocReadSet& readset, size_t min_kmer_count) {

    //cerr << "Building graph...\n"; //debug

    clear();

    //
    // Count all kmers.
    //
    for (const Read& r : readset.reads()) {

        // Build the first kmer.
        DNASeq4::iterator next_nt = r.seq.begin();
        Kmer km = Kmer(km_len_, next_nt, r.seq.end());
        if (km.empty()) {
            cerr << "Oops, no " << km_len_ << "-mers in " << r.seq.str() << "\n"; //
            continue;
        }

        // Record it.
        ++map_[km].count;

        // Walk the sequence.
        while (next_nt != r.seq.end()) {
            size_t nt4 = next_nt.nt();
            if (nt4 == Nt4::n) {
                ++next_nt;
                km = Kmer(km_len_, next_nt, r.seq.end());
                if (km.empty())
                    // Not enough sequence remaining to make another kmer.
                    break;
            } else {
                km = km.succ(km_len_, Nt2::nt4_to_nt[nt4]);
                ++next_nt;
            }
            ++map_[km].count;
        }
    }
    //cerr << "Found " << map_.size() << " kmers in " << readset.reads().size() <<" reads.\n"; //debug

    //
    // Build the standalone nodes.
    //
    for(auto km=map_.begin(); km!=map_.end();) {
        if (km->second.count < min_kmer_count) {
            map_.erase(km++);
        } else {
            nodes_.push_back(Node(NodeData(km->first, km->second.count)));
            // Replace the count of the kmer with the index of the corresponding node.
            km->second.node = nodes_.size() - 1;
            ++km;
        }
    }
    assert(map_.size() == nodes_.size());
    //cerr << "Built " << nodes_.size() << " nodes.\n"; //debug

    //
    // Build the edges.
    //
    for (Node& n : nodes_) {
        // Check each possible successor kmer.
        for (size_t nt2=0; nt2<4; ++nt2) {
            auto km = map_.find(n.km().succ(km_len_, nt2));
            if (km != map_.end()) {
                if (&nodes_[km->second.node] == &n)
                    // homopolymer, omit the edge
                    continue;
                n.set_succ(nt2, &nodes_[km->second.node]);
                nodes_[km->second.node].set_pred(n.km().front(), &n);
            }
        }
    }

    //
    // Build the simple paths.
    // There are three types of simple path starts:
    // * convergence nodes (several predecessors)
    // * nodes without predecessors
    // * successors of divergence nodes (one predecessor that has several successors)
    //
    for (Node& n : nodes_) {
        if (n.n_pred() != 1) {
            n.sp_build();
            simple_paths_.insert(&n);
        }

        if (n.n_succ() > 1) {
            for (size_t nt2=0; nt2<4; ++nt2) {
                Node* s = n.succ(nt2);
                if (s != NULL) {
                    s->sp_build();
                    simple_paths_.insert(s);
                }
            }
        }
    }
    //cerr << "Built " << simple_paths_.size() << " simple paths.\n"; //debug
}

void Graph::clear() {
    nodes_.resize(0);
    map_.clear();
    simple_paths_.clear();
}

void Graph::dump_gfa(const std::string& path) {
    ofstream ofs (path);
    if (!ofs) {
        cerr << "Error: Failed to open '" << path << "' for writing.\n";
        throw exception();
    }

    // Write the header.
    ofs << "H\tVN:Z:1.0\n";

    // Write the contigs.
    for (Node* sp : simple_paths_)
        ofs << "S\t" << index_of(sp) << "\t" << sp->sp_path_str(km_len_) << "\tKC:i:" << sp->sp_cum_count() << "\n";

    // Write the edges.
    for (Node* sp : simple_paths_) {
        for (size_t nt2=0; nt2<4; ++nt2) {
            Node* succ = sp->sp_succ(nt2);
            if (succ != NULL)
                ofs << "L\t" << index_of(sp) << "\t+\t" << index_of(succ) << "\t+\tM" << (km_len_-1) << "\n";
        }
    }
}
