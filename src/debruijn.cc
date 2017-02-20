#include <iostream>
#include <iomanip>
#include <sstream>

#include "constants.h"
#include "debruijn.h"

using namespace std;

string Node::sp_path_str(size_t km_len) {
    is_spfirst(this);

    // Seeding the sequence with the entire first kmer gives the same sequences
    // as those of Minia. However this is not practical for visualization purposes
    // as the lengths of the sequences do not correspond to the number of nodes;
    // instead the lengths are `n_nodes + km_len - 1`.
    // Thus here we use the last nucleotide of the node/kmer to represent the
    // graph in FastG.
    //string s = d_.km.str(km_len); // Same as Minia
    string s;

    Node* n=this;
    while (n != sp_last_) {
        s.push_back(Nt2::nt_to_ch[n->d_.km.back(km_len)]);
        n = n->first_succ();
    }
    s.push_back(Nt2::nt_to_ch[n->d_.km.back(km_len)]);

    return s;
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
                ++next_nt;
                km = Kmer(km_len, next_nt, r.seq.end());
                if (km.empty())
                    // Not enough sequence remaining to make another kmer.
                    break;
            } else {
                km = km.succ(km_len, Nt2::nt4_to_nt[nt4]);
                ++next_nt;
            }
            ++map[km].count;
        }
    }
    //cerr << "Found " << map.size() << " kmers in " << readset.reads().size() <<" reads.\n"; //debug

    //
    // Build the standalone nodes.
    //
    for(auto km=map.begin(); km!=map.end();) {
        if (km->second.count < min_kmer_count) {
            map.erase(km++);
        } else {
            nodes.push_back(Node(NodeData(km->first, km->second.count)));
            // Replace the count of the kmer with the index of the corresponding node.
            km->second.node = nodes.size() - 1;
            ++km;
        }
    }
    assert(map.size() == nodes.size());
    //cerr << "Built " << nodes.size() << " nodes.\n"; //debug

    //
    // Build the edges.
    //
    for (Node& n : nodes) {
        // Check each possible successor kmer.
        for (size_t nt2=0; nt2<4; ++nt2) {
            auto km = map.find(n.km().succ(km_len, nt2));
            if (km != map.end()) {
                if (&nodes[km->second.node] == &n)
                    // homopolymer, omit the edge
                    continue;
                n.set_succ(nt2, &nodes[km->second.node]);
                nodes[km->second.node].set_pred(n.km().front(), &n);
            }
        }
    }

    //
    // Record nodes that don't have predecessors.
    //
    for (Node& n : nodes)
        if (n.n_pred() == 0)
            nodes_wo_preds.push_back(&n);
    //cerr << "Found " << nodes_wo_preds.size() << " nodes without predecessors.\n"; //debug

    //
    // Build the simple paths.
    //
    sp_visited.clear();
    for (Node* n : nodes_wo_preds)
        build_simple_paths(n);
    //cerr << "Built " << sp_visited.size() << " simple paths.\n"; //debug
}

void Graph::clear() {
    nodes.resize(0);
    map.clear();
    nodes_wo_preds.clear();
    sp_visited.clear();
}

void Graph::build_simple_paths(Node* sp_first) {
    if (sp_visited.count(sp_first))
        return;
    sp_visited.insert(sp_first);

    Node* n = sp_first;
    Node* s = n->first_succ();
    while (n->n_succ() == 1 && s->n_pred() == 1) {
        // Extend the simple path.
        n = s;
        s = n->first_succ();
    }
    sp_first->set_sp_last(n);
    n->set_sp_first(sp_first);

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
        // i.e. `s->n_pred() != 1`
        // The next node has several predecessors (as `s` at least has one
        // predecessor, `n`).
        assert(s->n_pred() != 0);
        build_simple_paths(s);
    }
}

void Graph::dump_fg(const string& fastg_path) {
    ofstream ofs (fastg_path);
    if (!ofs) {
        cerr << "Error: Failed to open '" << fastg_path << "' for writing.\n";
        throw exception();
    }

    sp_visited.clear();
    for (Node* n : nodes_wo_preds)
        dump_fg(n, ofs);
}

void Graph::dump_fg(Node* sp, ostream& os) {
    if (sp_visited.count(sp))
        return;
    sp_visited.insert(sp);

    // Write the header.
    os << ">" << fg_header(sp);

    // Write the neighboring unitigs, if any.
    if (sp->sp_n_succ() > 0) {
        os << ":";
        stringstream ss;
        for (size_t nt2=0; nt2<4; ++nt2) {
            Node* succ = sp->sp_succ(nt2);
            if (succ != NULL)
                ss << "," << fg_header(succ);
        }
        os << ss.str().substr(1) << ";\n";
    } else {
        os << "\n";
    }

    // Write the sequence.
    os << sp->sp_path_str(km_len) << "\n";

    // Recurse.
    for (size_t nt2=0; nt2<4; ++nt2) {
        Node* n = sp->sp_succ(nt2);
        if (n != NULL)
            dump_fg(n, os);
    }
}

string Graph::fg_header(Node* sp) {
    stringstream ss;
    ss << std::fixed << setprecision(1);
    size_t id = sp - nodes.data();
    ss << "NODE_" << id << "_length_" << sp->sp_n_nodes() << "_cov_" << sp->sp_mean_count() << "_ID_" << id;
    return ss.str();
}
