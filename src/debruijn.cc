#include <iostream>
#include <iomanip>
#include <sstream>

#include "constants.h"
#include "debruijn.h"

using namespace std;

SPath::SPath(Node* first) : first_(first), last_(NULL), d_() {
    d_.n_nodes = 1;
    d_.km_cumcount = first->count();

    Node* n = first_;
    Node* s = n->first_succ();
    while (n->n_succ() == 1 && s->n_pred() == 1) {
        // Extend the simple path.
        n = s;
        ++d_.n_nodes;
        d_.km_cumcount += n->count();

        s = n->first_succ();
    }
    last_ = n;
}

string SPath::contig_str(size_t km_len) {

    string s = first_->km().str(km_len);
    Node* n = first_;
    while (n != last_) {
        s.push_back(Nt2::nt_to_ch[n->km().back(km_len)]);
        n = n->first_succ();
    }
    s.push_back(Nt2::nt_to_ch[n->km().back(km_len)]);

    return s;
}

void Graph::rebuild(const CLocReadSet& readset, size_t min_kmer_count) {

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
        if (n.n_pred() != 1)
            simple_paths_.push_back(SPath(&n));

        if (n.n_succ() > 1) {
            for (size_t nt2=0; nt2<4; ++nt2) {
                Node* s = n.succ(nt2);
                if (s != NULL) {
                    simple_paths_.push_back(SPath(s));
                }
            }
        }
    }
    for (SPath& p : simple_paths_)
        // A separate loop is required because the vector might have resized.
        p.update_ptrs();
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

    // Write the simple paths.
    for (SPath& p : simple_paths_)
        // n.b. In principle the length of the contigs should be (n_nodes+km_len-1).
        // However for visualization purposes we use n_nodes (for now at least).
        ofs << "S\t" << index_of(p.first()) << "\t*\tLN:i:" << p.n_nodes() << "\tKC:i:" << p.km_cumcount() << "\n";

    // Write the edges.
    for (SPath& p : simple_paths_) {
        for (size_t nt2=0; nt2<4; ++nt2) {
            SPath* succ = p.succ(nt2);
            if (succ != NULL)
                ofs << "L\t" << index_of(p.first()) << "\t+\t" << index_of(succ->first()) << "\t+\tM" << (km_len_-1) << "\n";
        }
    }
}
