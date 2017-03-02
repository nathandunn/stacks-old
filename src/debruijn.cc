#include <iostream>
#include <iomanip>
#include <sstream>

#include "constants.h"
#include "debruijn.h"

using namespace std;

SPath::SPath(Node* first) : first_(first), last_(NULL), d_(), visitdata(NULL) {
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
                km = km.succ(km_len_, Nt2::from_nt4(nt4));
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

bool Graph::topo_sort() {
    if (!sorted_spaths_.empty())
        // Already sorted.
        return true;

    vector<uchar> visitdata; // 0/1s; whether each spath is a parent in the recursion
    visitdata.reserve(simple_paths_.size());

    // Note: no need to reset the SPath::visitdata's to NULL as this is
    // the first algo to run.

    for (SPath& p : simple_paths_) {
        if(!topo_sort(&p, visitdata)) {
            sorted_spaths_.resize(0);
            return false;
        }
    }

    return true;
}

bool Graph::topo_sort(SPath* p, vector<uchar>& visitdata) {
    if (p->visitdata != NULL) {
        if (*(uchar*) p->visitdata)
            // The recursion looped; not a DAG.
            return false;
        else
            // Joining a known path from a different root.
            return true;
    } else {
        visitdata.push_back(true);
        p->visitdata = (void*)&visitdata.back();
        for (size_t nt2=0; nt2<4; ++nt2) {
            SPath* s = p->succ(nt2);
            if (s != NULL)
                if (!topo_sort(s, visitdata))
                    return false;
        }
        *(uchar*)p->visitdata = false;

        sorted_spaths_.push_back(p);
    }

    return true;
}

bool Graph::find_best_path(vector<const SPath*>& best_path) {
    best_path.resize(0);

    assert(!empty());
    if(!topo_sort())
        // Not a DAG.
        return false;

    #ifdef DEBUG
    // This is unnecessary as we iterate on a sort.
    for (const SPath& p : simple_paths_)
        p.visitdata = NULL;
    #endif

    // Compute the best score at each node.
    vector<size_t> scores;
    scores.reserve(sorted_spaths_.size());

    for (const SPath* p : sorted_spaths_) {
        //n.b. Terminal nodes were added first.
        size_t succ_scores[4];
        for (size_t nt2=0; nt2<4; ++nt2) {
            const SPath* succ = p->succ(nt2);
            if (succ != NULL)
                //n.b. as the graph is sorted, succ->visitdata has been set.
                succ_scores[nt2] = *(size_t*)succ->visitdata;
            else
                succ_scores[nt2] = 0;
        }
        scores.push_back(*std::max_element(succ_scores, succ_scores+4) + p->km_cumcount());
        p->visitdata = &scores.back();
    }

    // Find the best starting node.
    auto p = sorted_spaths_.rbegin();
    const SPath* best_start = *p;
    size_t best_score = *(size_t*)best_start->visitdata;
    ++p;
    while(p != sorted_spaths_.rend()) {
        if (*(size_t*)(*p)->visitdata > best_score) {
            best_start = *p;
            best_score = *(size_t*)best_start->visitdata;
        }
        ++p;
    }

    // Find the best path.
    const SPath* curr = best_start;
    while(curr != NULL) {
        best_path.push_back(curr);

        size_t succ_scores[4];
        for (size_t nt2=0; nt2<4; ++nt2) {
            const SPath* succ = curr->succ(nt2);
            if (succ != NULL)
                succ_scores[nt2] = *(size_t*)succ->visitdata;
            else
                succ_scores[nt2] = 0;
        }

        size_t* best_succ_score = std::max_element(succ_scores, succ_scores+4);
        curr = curr->succ(best_succ_score - succ_scores);
        // if succ_scores was {0,0,0,0}, curr is null
    }

    return true;
}

void Graph::dump_gfa(const string& path) const {
    ofstream ofs (path);
    if (!ofs) {
        cerr << "Error: Failed to open '" << path << "' for writing.\n";
        throw exception();
    }

    // Write the header.
    ofs << "H\tVN:Z:1.0\n";

    // Write the simple paths.
    for (const SPath& p : simple_paths_)
        // n.b. In principle the length of the contigs should be (n_nodes+km_len-1).
        // However for visualization purposes we use n_nodes (for now at least).
        ofs << "S\t" << index_of(p.first()) << "\t*\tLN:i:" << p.n_nodes() << "\tKC:i:" << p.km_cumcount() << "\n";

    // Write the edges.
    for (const SPath& p : simple_paths_) {
        for (size_t nt2=0; nt2<4; ++nt2) {
            const SPath* succ = p.succ(nt2);
            if (succ != NULL)
                ofs << "L\t" << index_of(p.first()) << "\t+\t" << index_of(succ->first()) << "\t+\tM" << (km_len_-1) << "\n";
        }
    }
}
