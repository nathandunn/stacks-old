// -*-mode:c++; c-style:k&r; c-basic-offset:4;-*-
//
// Copyright 2017, Julian Catchen <jcatchen@illinois.edu>
//
// This file is part of Stacks.
//
// Stacks is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// Stacks is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with Stacks.  If not, see <http://www.gnu.org/licenses/>.
//

#include "SuffixTree.h"

STNode::~STNode()
{
    for (uint i = 0; i < NT4cnt; i++) {
	if (this->edges_[i] != NULL)
	    delete this->edges_[i];
    }
    this->suffix_link_ = NULL;
}

STEdge *
STNode::add_edge(Nt4 nuc, int pos)
{
    size_t n = nuc.index();

    if (this->edges_[n] != NULL)
        delete this->edges_[n];

    this->edges_[n] = new STEdge(pos);

    return this->edges_[n];
}

STNode *
STNode::add_suffix_link(STNode *node)
{
    this->suffix_link_ = node;
    return this->suffix_link_;
}

size_t
SuffixTree::align(DNASeq4 query, vector<pair<size_t, size_t> > &alns)
{
    int qcnt     = 0;
    int q_len    = query.length();
    int q_stop   = q_len - 1;
    int seq_stop = this->seq_.length() - 1;
    int next_pos = 0;

    STNode *active_node = this->root;
    Nt4     active_edge = active_node->edge(query[0]) != NULL ? query[0] : Nt4::$;
    int     node_pos    = active_node->edge(active_edge)->start();
    int     node_stop   = active_node->edge(active_edge)->end() == -1 ? seq_stop : active_node->edge(active_edge)->end();
    
    if (active_edge == Nt4::$)
        return 0;
    else {
        qcnt++;
        // cerr << "  Exiting the root on active edge " << char(active_edge) << "; qcnt: " << qcnt << "\n";
    }

    //
    // Walk the tree until we reach a non-matching nucleotide or hit a leaf node.
    //
    do {
        next_pos = node_pos + 1;

        if (next_pos > node_stop && active_node->edge(active_edge)->succ() == NULL) {
            //
            // We have reached a leaf node.
            //
            break;

        } else if (next_pos > node_stop && active_node->edge(active_edge)->succ() != NULL) {
            //
            // Traverse to the next node.
            //
            active_node = active_node->edge(active_edge)->succ();
            active_edge = active_node->edge(query[qcnt]) == NULL ? Nt4::$ : query[qcnt];
            //
            // The next required edge is not available in the successor node.
            //
            if (active_edge == Nt4::$)
                break;

            node_pos  = active_node->edge(active_edge)->start();
            node_stop = active_node->edge(active_edge)->end() == -1 ? seq_stop : active_node->edge(active_edge)->end();
            qcnt++;
            // cerr << "  Traversing to node " << active_node->id() << "; active edge: " << char(active_edge) << "; qcnt: " << qcnt << "\n";

        } else if (qcnt <= q_stop && this->seq_[next_pos] != query[qcnt]) {
            //
            // Next nucleotide does not match.
            //
            break;

        } else {
            node_pos++;
            qcnt++;
            // cerr << "    Matching character '" << char(this->seq_[next_pos]) << "' on active edge " << char(active_edge) << "; qcnt: " << qcnt << "\n";
        }
    } while (qcnt <= q_stop);

    if (qcnt < this->min_align)
        return 0;

    if (active_node->edge(active_edge)->succ() == NULL) {
        //
        // Are we at a leaf node?
        //
        size_t aln_pos = node_pos - qcnt + 1;
        alns.push_back(make_pair(aln_pos, qcnt));

    } else if (node_pos < node_stop) {
        //
        // Otherwise, traverse this path to the first leaf node to determine the alignment position.
        //
        size_t aln_pos = find_leaf_dist(active_node->edge(active_edge)->succ()) - qcnt + 1;
        alns.push_back(make_pair(aln_pos, qcnt));

    } else {
        //
        // Traverse the remaining paths out of active_node to determine all the alignments for this fragment.
        //
        vector<size_t> dists;

        find_all_leaf_dists(active_node->edge(active_edge)->succ(), dists);

        for (uint i = 0; i < dists.size(); i++)
            alns.push_back(make_pair(dists[i] - qcnt + 1, qcnt));
    }
    
    return 0;
}

size_t
SuffixTree::align(const char *query, vector<pair<size_t, size_t> > &alns)
{
    int qcnt     = 0;
    int q_len    = strlen(query);
    int q_stop   = q_len - 1;
    int seq_stop = this->seq_.length() - 1;
    int next_pos = 0;

    STNode *active_node = this->root;
    Nt4     active_edge = active_node->edge(Nt4(query[0])) != NULL ? Nt4(query[0]) : Nt4::$;
    int     node_pos    = active_node->edge(active_edge)->start();
    int     node_stop   = active_node->edge(active_edge)->end() == -1 ? seq_stop : active_node->edge(active_edge)->end();
    
    if (active_edge == Nt4::$)
        return 0;
    else
        qcnt++;

    //
    // Walk the tree until we reach a non-matching nucleotide or hit a leaf node.
    //
    do {
        next_pos = node_pos + 1;

        if (next_pos > node_stop && active_node->edge(active_edge)->succ() == NULL) {
            //
            // We have reached a leaf node.
            //
            break;

        } else if (next_pos > node_stop && active_node->edge(active_edge)->succ() != NULL) {
            //
            // Traverse to the next node.
            //
            active_node = active_node->edge(active_edge)->succ();
            active_edge = active_node->edge(Nt4(query[qcnt])) == NULL ? Nt4::$ : Nt4(query[qcnt]);
            //
            // The next required edge is not available in the successor node.
            //
            if (active_edge == Nt4::$)
                break;

            node_pos  = active_node->edge(active_edge)->start();
            node_stop = active_node->edge(active_edge)->end() == -1 ? seq_stop : active_node->edge(active_edge)->end();
            qcnt++;

        } else if (qcnt <= q_stop && this->seq_[next_pos] != Nt4(query[qcnt])) {
            //
            // Next nucleotide does not match.
            //
            break;

        } else {
            node_pos++;
            qcnt++;
        }
    } while (qcnt <= q_stop);

    if (qcnt < this->min_align)
        return 0;

    if (active_edge == Nt4::$) {
        //
        // Traverse the remaining paths out of active_node to determine all the alignments for this fragment.
        //
        vector<size_t> dists;

        find_all_leaf_dists(active_node, dists);

        for (uint i = 0; i < dists.size(); i++)
            alns.push_back(make_pair(dists[i] - qcnt + 1, qcnt));

    } else if (active_node->edge(active_edge)->succ() == NULL) {
        //
        // Are we at a leaf node?
        //
        size_t aln_pos = node_pos - qcnt + 1;
        alns.push_back(make_pair(aln_pos, qcnt));

    } else if (node_pos < node_stop) {
        //
        // Otherwise, traverse this path to the first leaf node to determine the alignment position.
        //
        size_t aln_pos = find_leaf_dist(active_node->edge(active_edge)->succ()) - qcnt + 1;
        alns.push_back(make_pair(aln_pos, qcnt));

    } else {
        //
        // Traverse the remaining paths out of active_node's successor to determine all the alignments for this fragment.
        //
        vector<size_t> dists;

        find_all_leaf_dists(active_node->edge(active_edge)->succ(), dists);

        for (uint i = 0; i < dists.size(); i++)
            alns.push_back(make_pair(dists[i] - qcnt + 1, qcnt));
    }
    
    return 0;
}

size_t
SuffixTree::find_leaf_dist(STNode *node)
{
    uint dist = 0;

    for (uint i = 0; i < NT4cnt; i++) {
        if (node->edge(i) != NULL) {
            uint end = node->edge(i)->end() == -1 ? this->seq_.length() - 1 : node->edge(i)->end();
            uint len = end - node->edge(i)->start() + 1;

            dist = len;

            if (node->edge(i)->succ() != NULL)
                dist += this->find_leaf_dist(node->edge(i)->succ());
            break;
        }
    }

    return dist;
}

size_t
SuffixTree::find_all_leaf_dists(STNode *node, vector<size_t> &dists)
{
    for (uint i = 0; i < NT4cnt; i++) {
        if (node->edge(i) != NULL) {
            uint end = node->edge(i)->end() == -1 ? this->seq_.length() - 1 : node->edge(i)->end();
            uint len = end - node->edge(i)->start() + 1;

            if (node->edge(i)->succ() != NULL) {
                this->find_leaf_dist(node->edge(i)->succ());

                for (uint j = 0; j < dists.size(); j++)
                    dists[j] += len;
            } else {
                dists.push_back(this->seq_.length() - 1 - len);
            }
        }
    }

    return 0;
}

size_t
SuffixTree::build_tree()
{
    //
    // We will build the tree in O(n) time using Ukkonen's algorithm.
    //   Algorithm guided by:
    //    1) Gusfield, D. Algorithms on Strings, Trees, and Sequences: Computer Science and Computational Biology. 1997. Chapters 5-6.
    //    2) Goller, J. (jogojapan). Accessed April 27, 2017. http://stackoverflow.com/questions/9452701/ukkonens-suffix-tree-algorithm-in-plain-english
    //

    uint    slen        = this->seq_.length();
    STNode *active_node = this->root;
    Nt4     active_edge = Nt4::$;
    int     seq_index   = 0;
    int     active_len  = 0;
    int     remainder   = 0;
    int     next_pos    = 0;
    int     end_pos     = -1;
    int     forward_len = 0;
    size_t  id          = 2;

    STNode *split_node, *prev_node;
    STEdge *old_edge;

    for (uint i = 0; i < slen; i++) {
        seq_index = i - remainder;

	//
        // With each letter in the string, we have one insertion remaining to make.
        //
        remainder++;

        //
        // Increment the end position for all prefixes in the suffix tree.
        //
	end_pos += 1;

	// cerr << "beginning step i: " << i
        //      << ";\n  active node: " << active_node->id()
        //      << ", active edge: " << char(active_edge)
        //      << ", active len: " << active_len
        //      << "; remainder: " << remainder
        //      << "; seq index: " << seq_index << "\n";
	// cerr << "  Adding '" << char(this->seq_[i]) << "' to the tree.\n";

	//
	// Set the active_edge if we can. Grab the next position for insertion in the suffix tree.
	//
	if (active_edge == Nt4::$ && active_node->edge(this->seq_[i]) != NULL)
	    active_edge = this->seq_[i];
	next_pos = active_edge == Nt4::$ ? -1 : active_node->edge(active_edge)->start() + active_len;
        
	//
	// Compare the next character in the suffix tree to the current character for insertion.
	//
	if (next_pos > -1 && this->seq_[i] == this->seq_[next_pos]) {
	    //
	    // If the character is already in the tree, do not insert it, but set the active point to it.
	    //
	    // cerr << "    Not inserting '" << char(this->seq_[i]) << "', it is in the tree.\n";
            active_len++;

	    //
	    // If, after incrementing the active_len, we have reached an internal node, move the active point to this node.
	    //
            forward_len += this->forward_nodes(&active_node, active_edge, active_len, seq_index, remainder, end_pos);

	} else {
	    //
	    // Otherwise, we have reached a suffix that is not yet present in the tree.
	    //   - We currently have remainder insertions queued up to add to the tree, which we will iterate over now.
	    //   - For each insertion: split the active edge exiting the active node at position active_len.
	    //       - Adjust the active point.
	    //   - If this is the second or later insertion of this round, connect the latter node to the former inserted node with a suffix link.
	    //   - We stop inserting suffixes if we find a suffix is already in the tree.

	    bool stop_insertion  = false;
            bool add_suffix_link = false;

	    while (stop_insertion == false && remainder > 0) {
                string suf = this->seq_.str().substr(seq_index, remainder);
		// cerr << "  Inserting suffix '" << suf << "' into the tree.\n";

                if (active_edge == Nt4::$ && active_node->edge(this->seq_[i]) != NULL)
                    active_edge = this->seq_[i];

		//
		// Is an edge, representing the next character to insert, not already present in the tree?
		//
		if (active_edge == Nt4::$) {
                    // cerr << "    Adding edge " << char(this->seq_[i]) << " to node " << active_node->id() << "\n";
		    active_node->add_edge(this->seq_[i], i);
                    if (active_node == this->root)
                        active_len++;
                    remainder--;

		} else if (active_len == 0 && active_node->edge(active_edge) != NULL) {
		    active_len++;
                    forward_len += this->forward_nodes(&active_node, active_edge, active_len, seq_index, remainder, end_pos);
		    stop_insertion = true;
                    // cerr << "  Stopping, active_edge " << char(active_edge) << " is already in the tree.\n";
		    continue;

		} else {
		    //
		    // Split the node necessary to insert the next suffix.
		    //
		    // cerr << "      Splitting node " << active_node->id() << " at active edge: " << char(active_edge) << "; active len: " << active_len << " (created node " << id << ")\n";
                    next_pos   = active_node->edge(active_edge)->start() + active_len;
		    split_node = new STNode(id);
		    old_edge   = active_node->edge(active_edge);
		    split_node->add_edge(this->seq_[next_pos], next_pos);
		    split_node->add_edge(this->seq_[i], i);

		    //
		    // Does the edge being split already point to an internal node? If so, reconnect the new node to the existing nodes.
		    //
		    if (old_edge->succ() != NULL) {
			split_node->edge(this->seq_[next_pos])->succ(old_edge->succ());
			split_node->edge(this->seq_[next_pos])->end(old_edge->end());
		    }

		    old_edge->succ(split_node);
		    old_edge->end(next_pos - 1);

		    if (add_suffix_link == true) {
			prev_node->add_suffix_link(split_node);
			// cerr << "      Adding suffix link between node " << prev_node->id() << " and " << split_node->id() << "\n";
		    }

                    add_suffix_link = true;
                    prev_node       = split_node;
                    remainder--;
		}

		id++;

		//
		// Update the active point. 
		//
		if (active_node != this->root) {
		    //
		    // If there is a suffix link, follow it.
		    //
		    if (active_node->suffix_link() != NULL) {
			// cerr << "      Following suffix link to node " << active_node->suffix_link()->id() << "\n";
			active_node = active_node->suffix_link();
                        seq_index++;
                        forward_len--;
			// cerr << "      Resetting the active node to " << active_node->id() << ", active edge: " << char(active_edge) << ".\n";
                        forward_len += this->forward_nodes(&active_node, active_edge, active_len, seq_index + forward_len, remainder, end_pos);
			// cerr << "          Active node: " << active_node->id() << ", active len: " << active_len << ", actvie edge: " << char(active_edge) << ".\n";

		    } else {
			//
			// Otherwise, reset the active node to the root.
			//
                        seq_index++;
                        active_node = this->root;
                        forward_len = 0;
                        if (active_node->edge(this->seq_[seq_index]) == NULL) {
                            active_edge = Nt4::$;
                        } else {
                            active_edge   = this->seq_[seq_index];
                            active_len    = end_pos - seq_index;
                        }
			// cerr << "      Resetting the active node to the root, active edge: " << char(active_edge) << ".\n";
                        forward_len += this->forward_nodes(&active_node, active_edge, active_len, seq_index, remainder, end_pos);
			// cerr << "          Active node: " << active_node->id() << ", active len: " << active_len << ", actvie edge: " << char(active_edge) << ".\n";
		    }

		} else {
		    //
		    // Rule 1: after insertion at root, active_node remains root, active len is reduced by one, active edge is advanced to the next suffix.
		    //
		    active_len--;
		    seq_index++;

                    if (active_len > 0) {
                        // cerr << "    Changing the active edge from: " << char(active_edge) << " to " << char(this->seq_[seq_index]) << ", active len: " << active_len << "\n";
                        active_edge = this->seq_[seq_index];
                        this->forward_nodes(&active_node, active_edge, active_len, seq_index, remainder, end_pos);
                    } else {
                        // cerr << "    Completed suffix insertions.\n";
                        active_edge = Nt4::$;
                    }
                }

		// cerr << "      Remainder: " << remainder << "\n";
	    }
	}

        // cerr << "  ending step i: " << i
        //      << "; active node: " << active_node->id()
        //      << ", active edge: " << char(active_edge)
        //      << ", active len: " << active_len
        //      << "; remainder: " << remainder
        //      << ", forward len: " << forward_len << "\n";
    }

    return 0;
}

inline int
SuffixTree::forward_nodes(STNode **active_node, Nt4 &active_edge, int &active_len, int seq_index, int &remainder, int end_pos)
{
    // cerr << "      Adjusting the active point (node " << (*active_node)->id() << ", active len: " << active_len << ") to the proper internal node.\n";

    if ((*active_node)->edge(active_edge) == NULL)
        return active_len;

    int local_len     = (*active_node)->edge(active_edge)->end() - (*active_node)->edge(active_edge)->start() + 1;
    int forwarded_len = 0;
    
    while (active_len >= local_len && (*active_node)->edge(active_edge)->succ() != NULL) {
	active_len    -= local_len;
	seq_index     += local_len;
        forwarded_len += local_len;
	*active_node = (*active_node)->edge(active_edge)->succ();
	// cerr << "      Forwarding active node to node " << (*active_node)->id() << ", active len: " << active_len << ", active edge: " << char(active_edge) << ", seq_index: " << seq_index << "\n";
	if (active_len == 0) {
	    active_edge = Nt4::$;
	} else {
	    active_edge = this->seq_[seq_index];
	    local_len   = (*active_node)->edge(active_edge)->end() == -1 ? end_pos : (*active_node)->edge(active_edge)->end();
	    local_len   = local_len - (*active_node)->edge(active_edge)->start() + 1;
	}
        // cerr << "        Reset active edge: " << char(active_edge) << "\n";
    }

    return forwarded_len;
}

size_t
SuffixTree::write_dot(ofstream &fh)
{
    string s = this->seq_.str();

    fh << "digraph G {\n";
    
    queue<STNode *> q;
    STNode *node;
    string  label;
    stringstream range;
    int     leaf_id = 1;

    q.push(this->root);

    while (q.size() > 0) {
        node = q.front();
        q.pop();

	fh << "  " << node->id() << "\n";

        for (uint i = 0; i < NT4cnt; i++) {
            if (node->edge(i) != NULL) {

                uint end = node->edge(i)->end() == -1 ? this->seq_.length() - 1 : node->edge(i)->end();
                uint len = end - node->edge(i)->start() + 1;
                
                label = s.substr(node->edge(i)->start(), len);
		range.str("");
		range << "/* s: " << node->edge(i)->start() << ", e: " << end << " */";

                if (node->edge(i)->succ() != NULL) {
                    q.push(node->edge(i)->succ());
                    fh << "  "
                       << node->id() << " -> "
                       << node->edge(i)->succ()->id()
                       << " [label=\"" << label << "\"]; " << range.str() << "\n";
                } else {
                    fh << "  leaf" << leaf_id << " [shape=point]\n"
		       << "  "
                       << node->id() << " -> leaf" << leaf_id
                       << " [label=\"" << label << "\"]; " << range.str() << "\n";
                    leaf_id++;
                }
            }
        }

        if (node->suffix_link() != NULL) {
            fh << "  "
               << node->id() << " -> "
               << node->suffix_link()->id() << " [style=dashed]\n";
        }
    }

    fh << "}\n";

    return 0;
}

size_t
SuffixTree::write_suffixes(ostream &fh)
{
    string s = this->seq_.str();

    vector<string>  suffixes;
    string          suffix;

    this->write_suffix(suffixes, suffix, this->root);

    sort(suffixes.begin(), suffixes.end(), compare_str_len);

    for (uint i = 0; i < suffixes.size(); i++)
        fh << suffixes[i] << "\n";
    
    return 0;
}

size_t
SuffixTree::write_suffix(vector<string> &suffixes, string suffix, STNode *node)
{
    for (uint i = 0; i < NT4cnt; i++) {
        if (node->edge(i) != NULL) {
            uint end = node->edge(i)->end() == -1 ? this->seq_.length() - 1 : node->edge(i)->end();
            uint len = end - node->edge(i)->start() + 1;
            suffix += this->seq_.str().substr(node->edge(i)->start(), len);

            if (node->edge(i)->succ() != NULL) {
                this->write_suffix(suffixes, suffix, node->edge(i)->succ());
            } else {
                suffixes.push_back(suffix);
            }
            suffix = suffix.substr(0, suffix.length() - len);
        }
    }

    return 0;
}

bool
compare_staln(STAln a, STAln b)
{
    return a.subj_pos == b.subj_pos;
}
