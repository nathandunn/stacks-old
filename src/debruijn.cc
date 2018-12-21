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

void SPath::erase(size_t km_len) {
    if (first_ != last_) {
        assert(first_->n_succ() == 1);
        assert(last_->n_pred() == 1);
    }
    // Disconnect the SPath from the graph.
    for (size_t nt2=0; nt2<4; ++nt2) {
        SPath* p = pred(nt2);
        SPath* s = succ(nt2);
        if (p != NULL && p != this)
            p->last_->rm_succ(size_t(first_->km().back(km_len)), first_);
        if (s != NULL && s != this)
            last_->rm_succ(nt2, s->first_);
    }
    // Clear the SPath.
    Node* next = NULL;
    for (Node* n = first_; ; n = next) {
        if (n != last_) {
            assert(n->n_succ() == 1);
            next = n->first_succ();
        }
        *n = Node();
        if (n == last_)
            break;
    }
    *this = SPath();
}

void SPath::merge_forward() {
    assert(n_succ() == 1);
    SPath& s = *first_succ();
    assert(s.n_pred() == 1);
    assert(s.pred(size_t(this->last()->km().front())) == this);
    // Rebuild the path.
    *this = SPath(first_);
    assert(last_ == s.last_);
    // Clear the successor SPath.
    s = SPath();
    // Update the Node::sp_'s.
    for (Node* n = first_; ; n = n->first_succ()) {
        n->sp_ = this;
        if (n == last_)
            break;
    }
}

void Graph::rebuild(const vector<const DNASeq4*>& reads, size_t min_kmer_count) {

    //cerr << "Building graph...\n"; //debug

    clear();

    //
    // Fill the kmer map & count all kmers.
    //
    Kmer km;
    for (const DNASeq4* s : reads) {
        Kmerizer kmers {km_len_, *s};
        Kmer km;
        while ((km = kmers.next()))
            ++map_[km].count;
    }
    //cerr << "Found " << map_.size() << " kmers in " << readset.reads().size() <<" reads.\n"; //debug
    // Remove homopolymers.
    for (Nt2 nt : Nt2::all)
        map_.erase(Kmer::homopolymer(km_len_, nt));

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
        assert(!n.empty());
        // Check each possible successor kmer.
        for (size_t nt2=0; nt2<4; ++nt2) {
            auto km = map_.find(n.km().succ(km_len_, nt2));
            if (km != map_.end()) {
                assert(&nodes_[km->second.node] != &n); // Homopolymers were removed.
                assert(km->first.back(km_len_) == Nt2(nt2));
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
        assert(!n.empty());
        if (n.n_pred() != 1)
            simple_paths_.push_back(SPath(&n));

        if (n.n_succ() > 1) {
            for (size_t nt2=0; nt2<4; ++nt2) {
                Node* s = n.succ(nt2);
                if (s != NULL && s->n_pred() == 1) {
                    simple_paths_.push_back(SPath(s));
                }
            }
        }
    }
    for (SPath& p : simple_paths_) {
        // A separate loop is required because the vector might have resized.
        assert(!p.empty()); // We've just built the graph; no pieces should have been deleted yet!
        p.update_ptrs();
    }
    //cerr << "Built " << simple_paths_.size() << " simple paths.\n"; //debug
}

bool Graph::remove_cycles() {
    if (!sorted_spaths_.empty())
        // The graph is already acyclic.
        return false;
    return remove_microsat_dimer_cycles();
}

bool Graph::remove_microsat_dimer_cycles() {
    bool edited = false;
    for (SPath& p : simple_paths_) {
        if (p.empty())
            continue;
        if (p.n_nodes() == 1) {
            for (size_t nt2=0; nt2<4; ++nt2) {
                SPath* s = p.succ(nt2);
                if (s == NULL)
                    continue;
                if (s->n_nodes() == 1 && s->succ(size_t(p.last()->km().back(km_len_))) == &p) {
                    if (remove_microsat_dimer_cycle(p, *s))
                        edited = true;
                    break;
                }
            }
        } else if (p.n_nodes() == 2) {
            // This is actually a degenerate case of the above. If no additional paths break the
            // simple path, the two dimer nodes may, depending on the oddness of the kmer
            // size and of the microsatellite tract length, collapse in one single simple
            // path that loops on itself. (e.g. GATATG/k=3: GAT-ATA-TAT-ATG.)
            for (size_t nt2=0; nt2<4; ++nt2) {
                if (p.succ(nt2) == &p) {
                    p.erase(km_len_);
                    edited = true;
                    break;
                }
            }
        }
    }
    return edited;
}

bool Graph::remove_microsat_dimer_cycle(SPath& p, SPath& q) {
    // Make sure we god this right.
    size_t p_front = size_t(p.first()->km().front());
    size_t q_front = size_t(q.first()->km().front());
    size_t p_back = size_t(p.last()->km().back(km_len_));
    size_t q_back = size_t(q.last()->km().back(km_len_));
    assert(p.succ(q_back) == &q && p.pred(q_front) == &q);
    assert(q.succ(p_back) == &p && q.pred(p_front) == &p);
    // Detect the configuration we're in.
    SPath* rm;
    SPath* keep;
    assert(p.n_succ() >= 1 && p.n_pred() >= 1); // c.f. above.
    assert(q.n_succ() >= 1 && q.n_pred() >= 1);
    size_t p_external = bool(p.n_succ() - 1) + bool(p.n_pred() - 1);
    size_t q_external = bool(q.n_succ() - 1) + bool(q.n_pred() - 1);
    if (p_external == 2 && q_external == 2) {
        // Most general case; we ignore it as there's no clear solution and it's
        // rare anyway.
        #ifdef DEBUG
        cout << "DEBUG: Ignored a microsat dimer that is in the most general configuration.\n";
        #endif
        return false;
    } else if (p_external < 2 && q_external < 2) {
        DOES_NOT_HAPPEN;
        // Because this is the degenerate case when the two nodes of the cycle are
        // in the same SPath, and it is handled as we detect the cycle.
    } else if (p_external < 2) {
        assert(q_external == 2);
        rm = &p;
        keep = &q;
    } else if (q_external < 2) {
        assert(p_external == 2);
        rm = &q;
        keep = &p;
    } else {
        DOES_NOT_HAPPEN;
    }
    // Erase one of the two SPaths.
    rm->erase(km_len_);
    // Merge the other SPath as necessary.
    if (keep->n_succ() == 1 && keep->first_succ()->n_pred() == 1)
        keep->merge_forward();
    if (keep->n_pred() == 1) {
        SPath* pred = NULL;
        for (size_t nt2=0; nt2<4; ++nt2) {
            if (keep->pred(nt2) != NULL) {
                pred = keep->pred(nt2);
                break;
            }
        }
        assert(pred != NULL);
        if (pred->n_succ() == 1) {
            pred->merge_forward();
            assert(keep->empty());
        }
    }
    return true;
}

bool Graph::topo_sort() const {
    if (!sorted_spaths_.empty())
        // Already sorted.
        return true;

    vector<uchar> visitdata; // 0/1s; whether each spath is a parent in the recursion
    visitdata.reserve(simple_paths_.size());
    for (const SPath& p : simple_paths_)
        if (!p.empty())
            p.visitdata = NULL;

    for (const SPath& p : simple_paths_) {
        if (p.empty())
            continue;
        if(!topo_sort(&p, visitdata)) {
            sorted_spaths_.resize(0);
            return false;
        }
    }

    return true;
}

bool Graph::topo_sort(const SPath* p, vector<uchar>& visitdata) const {
    if (p->visitdata != NULL) {
        if (*(uchar*) p->visitdata)
            // The recursion looped; not a DAG.
            return false;
        else
            // Joining a known path from a different root.
            return true;
    } else {
        visitdata.push_back(true); // n.b. Enough memory was reserved.
        p->visitdata = (void*)&visitdata.back();
        for (size_t nt2=0; nt2<4; ++nt2) {
            const SPath* s = p->succ(nt2);
            if (s != NULL)
                if (!topo_sort(s, visitdata))
                    return false;
        }
        *(uchar*)p->visitdata = false;

        sorted_spaths_.push_back(p);
    }

    return true;
}

bool Graph::find_best_path(vector<const SPath*>& best_path) const {
    best_path.resize(0);

    assert(!empty());
    if(!topo_sort())
        // Not a DAG.
        return false;

    #ifdef DEBUG
    // This is unnecessary as we iterate on a sort.
    for (const SPath& p : simple_paths_)
        if (!p.empty())
            p.visitdata = NULL;
    #endif

    // Compute the best score at each node.
    vector<size_t> scores;
    scores.reserve(sorted_spaths_.size());

    for (const SPath* p : sorted_spaths_) {
        assert(!p->empty());
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

void Graph::dump_gfa(const string& path, bool individual_nodes) const {
    ofstream ofs (path);
    if (!ofs) {
        cerr << "Error: Failed to open '" << path << "' for writing.\n";
        throw exception();
    }

    // Write the header.
    ofs << "H\tVN:Z:1.0\n";

    if (!individual_nodes) {
        // Write the vertices.
        for (const SPath& p : simple_paths_)
            if (!p.empty())
                // n.b. In principle the length of the contigs should be (n_nodes+km_len-1).
                // However for visualization purposes we use n_nodes (for now at least).
                ofs << "S\t" << index_of(p.first()) << "\t*\tLN:i:" << p.n_nodes() << "\tKC:i:" << p.km_cumcount() << "\n";

        // Write the edges.
        for (const SPath& p : simple_paths_) {
            if (p.empty())
                continue;
            for (size_t nt2=0; nt2<4; ++nt2) {
                const SPath* succ = p.succ(nt2);
                if (succ != NULL)
                    ofs << "L\t" << index_of(p.first()) << "\t+\t" << index_of(succ->first()) << "\t+\tM" << (km_len_-1) << "\n";
            }
        }
    } else {
        // Write the vertices.
        for (const Node& n : nodes_) {
            if (n.empty())
                continue;
            ofs << "S\t" << index_of(&n) << "\t*\tLN:i:1\tKC:i:" << n.count() << "\tseq:" << n.km().str(km_len_) << "\n";
        }
        // Write the edges.
        for (const Node& n : nodes_) {
            if (n.empty())
                continue;
            for (size_t nt2=0; nt2<4; ++nt2) {
                const Node* succ = n.succ(nt2);
                if (succ != NULL)
                    ofs << "L\t" << index_of(&n) << "\t+\t" << index_of(succ) << "\t+\tM" << (km_len_-1) << "\n";
            }
        }
    }
}

vector<vector<const SPath*>> Graph::components() {
    map<const void*, vector<const SPath*>> components_map;
    for (SPath& p : simple_paths_)
        if (!p.empty())
            p.visitdata = NULL;
    for (SPath& p : simple_paths_) {
        if (p.empty())
            continue;
        propagate_component_id(&p, NULL);
        components_map[p.visitdata].push_back(&p);
    }
    vector<vector<const SPath*>> components;
    for (auto& c : components_map)
        components.push_back(move(c.second));
    return components;
}

void Graph::propagate_component_id(const SPath* p, void* id) {
    if (p->visitdata == NULL) {
        if (id == NULL)
            // (One simple path address serves as an ID for each component.)
            id = const_cast<SPath*>(p);
        p->visitdata = id;
        const SPath* neighbor;
        for (size_t nt2=0; nt2<4; ++nt2) {
            neighbor = p->pred(nt2);
            if (neighbor != NULL)
                propagate_component_id(neighbor, id);
            neighbor = p->succ(nt2);
            if (neighbor != NULL)
                propagate_component_id(neighbor, id);
        }
    }
}
