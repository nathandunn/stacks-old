// -*-mode:c++; c-style:k&r; c-basic-offset:4;-*-
//
// Copyright 2010, Julian Catchen <jcatchen@uoregon.edu>
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

//
// ustacks -- build denovo stacks
//
// Julian Catchen
// jcatchen@uoregon.edu
// University of Oregon
//
// $Id$
//
#include "ustacks.h"

//
// Global variables to hold command-line options.
//
file_type in_file_type;
string    in_file;
string    out_path;
int       num_threads       = 1;
int       sql_id            = 0;
bool      set_kmer_len      = true;
int       kmer_len          = 0;
int       min_merge_cov     = 2;
int       dump_graph        = 0;
int       retain_rem_reads  = false;
int       deleverage_stacks = 0;
int       remove_rep_stacks = 0;
int       max_utag_dist     = 2;
int       max_rem_dist      = 4;
double    cov_mean          = 0.0;
double    cov_stdev         = 0.0;
double    cov_scale         = 1;
int       deleverage_trigger;
int       removal_trigger;
//
// For use with the multinomial model to call fixed nucleotides.
//
modelt model_type        = snp;
double p_freq            = 0.5;
double barcode_err_freq  = 0.0;

int main (int argc, char* argv[]) {

    parse_command_line(argc, argv);

    //
    // Set the max remainder distance to be greater than the max_utag_dist.
    //
    max_rem_dist = max_utag_dist + 2;

    //
    // Set the number of OpenMP parallel threads to execute.
    //
    #ifdef _OPENMP
    omp_set_num_threads(num_threads);
    #endif

    HashMap           radtags;
    map<int, Seq *>   remainders;
    set<int>          merge_map;
    map<int, Stack *> unique;

    load_radtags(in_file, radtags);

    reduce_radtags(radtags, unique, remainders);

    if (cov_mean == 0 || cov_stdev == 0)
	calc_coverage_distribution(unique, cov_mean, cov_stdev);

    cerr << "Coverage mean: " << cov_mean << "; stdev: " << cov_stdev << "\n";

    calc_triggers(cov_mean, cov_stdev, deleverage_trigger, removal_trigger);

    cerr << "Deleveraging trigger: " << deleverage_trigger << "; Removal trigger: " << removal_trigger << "\n";

    map<int, MergedStack *> merged;

    populate_merged_tags(unique, merged);

    if (remove_rep_stacks) {
	cerr << "Calculating distance for removing repetitive stacks.\n";
	calc_kmer_distance(merged, 1);
	cerr << "Removing repetitive stacks.\n";
	remove_repetitive_stacks(unique, merged);
    }

    for (int i = 1; i <= max_utag_dist; i++) {
	cerr << "Calculating distance, round " << i << "\n";
	calc_kmer_distance(merged, i);

	cerr << "  Merging radtags, round " << i << "\n";
	merge_radtags(unique, merged, merge_map, i);

	call_consensus(merged, unique, remainders, false);

	merge_map.clear();
    }

    calc_merged_coverage_distribution(unique, merged);

    //dump_merged_tags(merged);

    cerr << "Merging remainder radtags\n";
    merge_remainders(merged, remainders);

    // Call the consensus sequence again, now that remainder tags have been merged.
    call_consensus(merged, unique, remainders, true);

    count_raw_reads(unique, merged);

    cerr << "Writing results\n";
    write_results(merged, unique, remainders);

    return 0;
}

int merge_remainders(map<int, MergedStack *> &merged, map<int, Seq *> &rem) {
    map<int, Seq *>::iterator it;
    int j, k;

    // OpenMP can't parallelize random access iterators, so we convert
    // our map to a vector of integer keys.
    vector<int> keys;
    for (it = rem.begin(); it != rem.end(); it++) 
	keys.push_back(it->first);

    cerr << "  " << keys.size() << " remainder sequences left to merge.\n";

    //
    // Calculate the number of k-mers we will generate. If kmer_len == 0,
    // determine the optimal length for k-mers.
    //
    int con_len   = strlen(merged.begin()->second->con);
    if (set_kmer_len) kmer_len = determine_kmer_length(con_len, max_rem_dist);
    int num_kmers = con_len - kmer_len + 1;

    cerr << "  Number of kmers per sequence: " << num_kmers << "\n";

    //
    // Calculate the minimum number of matching k-mers required for a possible sequence match.
    //
    int min_hits = calc_min_kmer_matches(kmer_len, max_rem_dist, con_len, true);

    KmerHashMap kmer_map;
    populate_kmer_hash(merged, kmer_map, kmer_len);
    int utilized = 0;

    #pragma omp parallel private(it, k)
    { 
        #pragma omp for schedule(dynamic) 
	for (j = 0; j < (int) keys.size(); j++) {
	    it = rem.find(keys[j]);
	    Seq *r = it->second;

            //
            // Generate the k-mers for this remainder sequence
            //
            vector<char *> rem_kmers;
            generate_kmers(r->seq, kmer_len, num_kmers, rem_kmers);

            map<int, int> hits;
            vector<int>::iterator map_it;
            //
            // Lookup the occurances of each remainder k-mer in the MergedStack k-mer map
            //
            for (k = 0; k < num_kmers; k++) {
                if (kmer_map.find(rem_kmers[k]) != kmer_map.end())
                    for (map_it  = kmer_map[rem_kmers[k]].begin();
                         map_it != kmer_map[rem_kmers[k]].end();
                         map_it++)
                        hits[*map_it]++;
            }

            //
            // Iterate through the list of hits. For each hit that has more than min_hits
            // check its full length to verify a match.
            //
	    map<int, int> dists;
            map<int, int>::iterator hit_it;
            for (hit_it = hits.begin(); hit_it != hits.end(); hit_it++) {
                if (hit_it->second < min_hits) continue;

                int d = dist(merged[hit_it->first], r);
                //
                // Store the distance between these two sequences if it is
                // below the maximum distance
                //
                if (d <= max_rem_dist) {
                    dists[hit_it->first] = d;
                }
            }

            //
            // Free the k-mers we generated for this remainder read
            //
            for (k = 0; k < num_kmers; k++)
                delete rem_kmers[k];

	    // Check to see if there is a uniquely low distance, if so,
	    // merge this remainder tag. If not, discard it, since we
	    // can't locate a single best-fitting Stack to merge it into.
	    map<int, int>::iterator s;
	    int min_id = -1;
	    int count  =  0;
	    int dist   =  max_rem_dist + 1;

	    for (s = dists.begin(); s != dists.end(); s++) {
		if ((*s).second < dist) {
		    min_id = (*s).first;
		    count  = 1;
		    dist   = (*s).second;
		} else if ((*s).second == dist) {
		    count++;
		}
	    }

	    // Found a merge partner.
	    if (min_id >= 0 && count == 1) {
		r->utilized = true;
                #pragma omp critical
		{
		    merged[min_id]->remtags.push_back(it->first);
		    utilized++;
		}
	    }
	}
    }

    cerr << "  Matched " << utilized << " remainder reads; unable to match " << keys.size() - utilized << " remainder reads.\n";

    return 0;
}

int call_alleles(MergedStack *mtag, vector<char *> &reads, vector<read_type> &read_types) {
    int     row;
    int     height = reads.size();
    string  allele;
    char   *base;
    vector<SNP *>::iterator snp;

    if (mtag->snps.size() == 0)
	return 0;

    for (row = 0; row < height; row++) {
	allele.clear();

	//
	// Only call a haplotype from primary reads.
	//
	if (read_types[row] == secondary) continue;

	for (snp = mtag->snps.begin(); snp != mtag->snps.end(); snp++) {
	    base = reads[row];
	    base = base + (*snp)->col;

	    //
	    // Check to make sure the nucleotide at the location of this SNP is
	    // of one of the two possible states the multinomial model called.
	    //
	    if (*base == (*snp)->rank_1 || *base == (*snp)->rank_2) 
		allele += *base;
	    else
		break;
	}

	if (allele.size() == mtag->snps.size())
	    mtag->alleles[allele]++;
    }

    return 0;
}

int call_consensus(map<int, MergedStack *> &merged, map<int, Stack *> &unique, map<int, Seq *> &rem, bool invoke_model) {
    //
    // OpenMP can't parallelize random access iterators, so we convert
    // our map to a vector of integer keys.
    //
    map<int, MergedStack *>::iterator it;
    vector<int> keys;
    for (it = merged.begin(); it != merged.end(); it++) 
	keys.push_back(it->first);

    int i;
    #pragma omp parallel private(i)
    { 
        #pragma omp for schedule(dynamic) 
	for (i = 0; i < (int) keys.size(); i++) {
	    MergedStack *mtag;
	    Stack *utag;

	    mtag = merged[keys[i]];

	    //
	    // Create a two-dimensional array, each row containing one read. For
	    // each unique tag that has been merged together, add the sequence for
	    // that tag into our array as many times as it originally occurred. 
	    //
	    vector<int>::iterator j;
	    vector<char *> reads;
	    vector<read_type> read_types;

	    for (j = mtag->utags.begin(); j != mtag->utags.end(); j++) {
		utag = unique[*j];

		for (uint k = 0; k < utag->count; k++) {
		    reads.push_back(utag->seq);
		    read_types.push_back(primary);
		}
	    }

	    // For each remainder tag that has been merged into this Stack, add the sequence. 
	    for (j = mtag->remtags.begin(); j != mtag->remtags.end(); j++) {
		reads.push_back(rem[*j]->seq);
		read_types.push_back(secondary);
	    }

	    //
	    // Iterate over each column of the array and call the consensus base.
	    //
	    int row, col;
	    int length = strlen(reads[0]);
	    int height = reads.size();
	    string con;
	    map<char, int> nuc;
	    map<char, int>::iterator max, n;
	    char *base;

	    for (col = 0; col < length; col++) {
		nuc['A'] = 0; 
		nuc['C'] = 0;
		nuc['G'] = 0;
		nuc['T'] = 0;

		for (row = 0; row < height; row++) {
		    base = reads[row];
		    base = base + col;
		    //cerr << "    Row: " << row << " Col: " << col << " Base: " << *base << "\n";
		    nuc[*base]++;
		}

		//
		// Find the base with a plurality of occurances and call it.
		//
		max = nuc.end();

		for (n = nuc.begin(); n != nuc.end(); n++) {

		    if (max == nuc.end() || n->second > max->second)
			max = n;
		}
		con += max->first;

		// Search this column for the presence of a SNP
		if (invoke_model) 
		    model_type == snp ? 
                        call_multinomial_snp(mtag, col, nuc, false) :
                        call_multinomial_fixed(mtag, col, nuc);
	    }

	    if (invoke_model) {
		call_alleles(mtag, reads, read_types);

                if (model_type == fixed) {
                    //
                    // Mask nucleotides that are not fixed.
                    //
                    vector<SNP *>::iterator s;
                    for (s = mtag->snps.begin(); s != mtag->snps.end(); s++) {
                        con.replace((*s)->col, 1, "N");
                    }
                }
            }

	    mtag->add_consensus(con.c_str());
	}
    }

    return 0;
}

int populate_merged_tags(map<int, Stack *> &unique, map<int, MergedStack *> &merged) {
    map<int, Stack *>::iterator i;
    map<int, MergedStack *>::iterator it_new, it_old;
    Stack      *utag;
    MergedStack *mtag;
    int        k = 0;

    it_old = merged.begin();

    for (i = unique.begin(); i != unique.end(); i++) {
	utag = (*i).second;
	mtag = new MergedStack;

	mtag->id    = k;
	mtag->count = utag->count;
	mtag->utags.push_back(utag->id);
	mtag->add_consensus(utag->seq);

	// Insert the new MergedStack giving a hint as to which position
	// to insert it at.
	it_new = merged.insert(it_old, pair<int, MergedStack *>(k, mtag));
	it_old = it_new;
	k++;
    }

    return 0;
}

int merge_radtags(map<int, Stack *> &unique, map<int, MergedStack *> &merged, set<int> &merge_map, int round) {
    map<int, MergedStack *> new_merged;
    map<int, MergedStack *>::iterator i;
    MergedStack *tag_1, *tag_2;

    //cerr << "Tags to merge: " << merged.size() << "\n";
    int id = 0;

    for (i = merged.begin(); i != merged.end(); i++) {
	tag_1 = i->second;

	//
	// This tag may already have been merged by an earlier operation.
	//
	if (merge_map.count(tag_1->id) > 0)
	    continue;

	set<int>                          unique_merge_list;
	queue<int>                        merge_list;
	pair<set<int>::iterator,bool>     ret;
	vector<pair<int, int> >::iterator k;

	//
	// If this tag is masked, do not try to merge it.
	//
	if (tag_1->masked) {
	    unique_merge_list.insert(tag_1->id);
	    tag_2 = merge_tags(merged, unique_merge_list, id);
	    new_merged.insert(pair<int, MergedStack *>(id, tag_2));
	    id++;
	    continue;
	}

	//
	// Construct a list of MergedStacks that are within a particular distance
	// of this tag.
	//
	unique_merge_list.insert(tag_1->id);
	merge_list.push(tag_1->id);

	while (!merge_list.empty()) {
	    tag_2 = merged[merge_list.front()];
	    merge_list.pop();

	    for (k = tag_2->dist.begin(); k != tag_2->dist.end(); k++) {
                ret = unique_merge_list.insert(k->first);

		//
		// If this Tag has not already been added to the merge list (i.e. we were able
		// to insert it in to our unique_merge_list, which is a set), add it for consideration
		// later in the loop.
		//
		if (ret.second == true)
		    merge_list.push((*k).first);
	    }
	}

	//
	// Merge these tags together into a new MergedStack object.
	//
	tag_1 = merge_tags(merged, unique_merge_list, id);

	//
	// Record the nodes that have been merged in this round.
	//
	set<int>::iterator j;
	for (j = unique_merge_list.begin(); j != unique_merge_list.end(); j++)
	    merge_map.insert(*j);

	//
	// If the depth of coverage of the merged tag is greater than the deleverage trigger
	// then execute the deleveraging algorithm.
	//
	if (deleverage_stacks && 
	    (tag_1->count > deleverage_trigger ||
	     (int) unique_merge_list.size() > max_delv_stacks)) {
	    vector<MergedStack *> tags;
	    int num_clusters = determine_single_linkage_clusters(merged, unique_merge_list, tag_1->count);

	    //set<int>::iterator uit;
	    //cerr << "Unique merge list: ";
	    //for (uit = unique_merge_list.begin(); uit != unique_merge_list.end(); uit++) {
	    //cerr << "  " << *uit << "\n";
	    //}

	    deleverage(unique, merged, unique_merge_list, num_clusters, round, tags);

	    for (uint t = 0; t < tags.size(); t++) {
		tags[t]->id = id;
		new_merged.insert(pair<int, MergedStack *>(id, tags[t]));
		id++;
	    }
	    delete tag_1;

	} else {
	    new_merged.insert(pair<int, MergedStack *>(id, tag_1));
	    id++;
	}
    }

    merged = new_merged;

    return 0;
}

MergedStack *merge_tags(map<int, MergedStack *> &merged, set<int> &merge_list, int id) {
    set<int>::iterator    i;
    vector<int>::iterator j;
    MergedStack *tag_1, *tag_2;

    tag_1     = new MergedStack;
    tag_1->id = id;

    for (i = merge_list.begin(); i != merge_list.end(); i++) {
	tag_2 = merged[(*i)];

	tag_1->deleveraged     = tag_2->deleveraged     ? true : tag_1->deleveraged;
	tag_1->masked          = tag_2->masked          ? true : tag_1->masked;
	tag_1->blacklisted     = tag_2->blacklisted     ? true : tag_1->blacklisted;
	tag_1->lumberjackstack = tag_2->lumberjackstack ? true : tag_1->lumberjackstack;

	for (j = tag_2->utags.begin(); j != tag_2->utags.end(); j++)
	    tag_1->utags.push_back(*j);

        tag_1->count += tag_2->count;
    }

    return tag_1;
}

int remove_repetitive_stacks(map<int, Stack *> &unique, map<int, MergedStack *> &merged) {
    //
    // If enabled, check the depth of coverage of each unique tag, and remove 
    // from consideration any tags with depths greater than removal_trigger. These tags
    // are likely to be multiple repetitive sites that have been merged together. 
    // Because large stacks of unique tags are likely to also generate many one-off
    // sequencing error reads, remove all seqeunces that are a distance of one from
    // the RAD-Tag with high depth of coverage.
    //
    map<int, MergedStack *>::iterator i;
    vector<pair<int, int> >::iterator k;
    map<int, MergedStack *> new_merged;
    MergedStack *tag_1, *tag_2;
    set<int> already_merged;

    //
    // First, iterate through the stacks and populate a list of tags that will be removed
    // (those above the removal_trigger and those 1 nucleotide away). If we don't construct
    // this list first, we will inadvertantly merge short stacks that end up being a
    // single nucleotide away from one of the lumberjack stacks found later in the process.
    //

    int id = 0;

    //
    // Merge all stacks that are over the removal trigger with their nearest neighbors and
    // mask them so they are not further processed by the program.
    //
    for (i = merged.begin(); i != merged.end(); i++) {
	tag_1 = i->second;

	//
	// Don't process a tag that has already been merged.
	//
        if (already_merged.count(tag_1->id) > 0)
	    continue;


	if (tag_1->count > removal_trigger) {
            set<int> unique_merge_list;
            unique_merge_list.insert(tag_1->id);
	    already_merged.insert(tag_1->id);

	    for (k = tag_1->dist.begin(); k != tag_1->dist.end(); k++) {
                if (already_merged.count(k->first) == 0) {
                    already_merged.insert(k->first);
                    unique_merge_list.insert(k->first);
                }
	    }

	    tag_1->lumberjackstack = true;
	    tag_1->masked  = true;

            //
            // Merge these tags together into a new MergedStack object.
            //
            tag_2 = merge_tags(merged, unique_merge_list, id);
            tag_2->add_consensus(tag_1->con);

            new_merged.insert(make_pair(id, tag_2));
            id++;
	}
    }

    //
    // Move the non-lumberjack stacks, unmodified, into the new merged map.
    //
    for (i = merged.begin(); i != merged.end(); i++) {
	tag_1 = i->second;

	if (already_merged.count(tag_1->id) > 0)
	    continue;

	set<int> unique_merge_list;
	unique_merge_list.insert(tag_1->id);

	tag_2 = merge_tags(merged, unique_merge_list, id);
	tag_2->add_consensus(tag_1->con);

	new_merged.insert(make_pair(id, tag_2));
	id++;
    }

    cerr << "  Removed " << already_merged.size() << " stacks.\n";

    merged = new_merged;

    cerr << "  " << merged.size() << " stacks remain for merging.\n";

    return 0;
}

int determine_single_linkage_clusters(map<int, MergedStack *> &merged, set<int> &merge_list, int tag_depth) {
    set<int>::iterator i;
    MergedStack *tag;

    return 2;

    // 
    // We want to determine how many clusters to break this lumberjack stack into. We start
    // by dividing the total depth of the tag by the expected coverage. This assumes that every
    // tag has the same depth. However, it is quite common for one or two clusters within the
    // stack to have very large depths, so we want to adjust the total number of clusters we
    // expect there to be by accounting for very large clusters within the lumberjackstack.
    //
    //  Lumberjack stack:
    //  |--------------------------------------------------------------------------------|
    //  Expected stacks as determined by cov_mean parameter:
    //  |--------|--------|--------|--------|--------|--------|--------|--------|--------|
    //  Actual components of stack:
    //  |---------------------------------------|----------------|----|---|--|--|--|--|--|
    //
    double num_clusters = round(tag_depth / cov_mean);
    double cluster_cov;
    int    cluster_count = 0;

    if (num_clusters <= 1) return 2;

    for (i = merge_list.begin(); i != merge_list.end(); i++) {
	tag           = merged[*i];
	cluster_cov   = floor(tag->count / cov_mean); 

	//cerr << "Cluster coverage: " << cluster_cov << "; Tag depth: " << tag->count << "; Expected cov: " << cov_mean << "\n";

	if (cluster_cov > 1) {
	    num_clusters -= cluster_cov;
	    cluster_count++;
	}
    }
    num_clusters = cluster_count + num_clusters;

    //
    // Make sure there aren't more clusters than individual stacks
    //
    num_clusters = num_clusters > merge_list.size() ? merge_list.size() : num_clusters;

    //cerr << "Number of clusters: " << num_clusters 
    //     << "; unique_merge_list size: " << merge_list.size() 
    //     << "; expected coverage: " << cov_mean
    //     << "; Lumberjackstack depth: " << tag_depth << "\n";
    
    return (int) num_clusters;
}

int deleverage(map<int, Stack *> &unique, 
	       map<int, MergedStack *> &merged, 
	       set<int> &merge_list, 
	       int num_clusters, 
	       int round, 
	       vector<MergedStack *> &deleveraged_tags) {
    set<int>::iterator i;
    vector<pair<int, int> >::iterator j;
    MergedStack *tag, *tag_1, *tag_2;
    vector<int> keys;

    //
    // Create a two-dimensional map to hold distances between nodes.
    //
    map<int, map<int, double> > dist_map, scaled_dist_map;

    for (i = merge_list.begin(); i != merge_list.end(); i++)
	keys.push_back(*i);

    //cerr << "    Deleveraging: " << merged[keys[0]]->con << "\n";

    //
    // There must be at least two tags before we can cluster the data.
    //
    if (keys.size() <= 2) {
	set<int> e;
	for (uint k = 0; k < keys.size(); k++)
	    e.insert(keys[k]);

	tag = merge_tags(merged, e, 0);
	deleveraged_tags.push_back(tag);

	return 0;
    }

    //
    // We want to measure the distance between all nodes in the merge_list.
    //
    for (uint k = 0; k < keys.size(); k++) {
	for (uint l = 0; l < keys.size(); l++) {
	    tag_1 = merged[keys[k]];
	    tag_2 = merged[keys[l]];

	    int d = dist(tag_1, tag_2);
	    dist_map[keys[k]][keys[l]] = d;
	    dist_map[keys[l]][keys[k]] = d;
	}
    }
    
    keys.clear();

    uint s, t, depth_1, depth_2, id;
    map<int, map<int, double> >::iterator q;
    map<int, double>::iterator r;
    double scale, d;

    // 
    // Record the IDs for the stacks in the distance map
    //
    q = dist_map.begin();
    for (r = (*q).second.begin(); r != (*q).second.end(); r++) {
	//cerr << "  Key: " << (*r).first << "; Number: " << (*q).second.count((*r).first) << "\n";
	keys.push_back((*r).first);
    }

    //
    // Scale the distances between nodes
    //
    for (s = 0; s < keys.size(); s++) {
	for (t = 0; t < keys.size(); t++) {
	    if (s == t) continue;

	    depth_1 = merged[keys[s]]->count;
	    depth_2 = merged[keys[t]]->count;
	    d       = abs(depth_2 - depth_1);

	    scale = d == 0 ? 1 : log(d);
	    //cerr << 
	    //"Distance between " << keys[s] << " (" << depth_1 << ") and " << keys[t] << " (" << depth_2 << "): " << 
	    //d << " (" << dist_map[keys[s]][keys[t]] << 
	    //"; scale: " << scale << 
	    //"; total: " << dist_map[keys[s]][keys[t]] * scale << ")\n";
	    scaled_dist_map[keys[s]][keys[t]] = dist_map[keys[s]][keys[t]] * scale;
	}
    }

    map<int, set<int> > cluster_map;
    map<int, set<int> >::iterator c;
    set<int>::iterator it;
    stringstream gout_file;

    hclust(keys, scaled_dist_map, num_clusters, cluster_map);

    if (dump_graph) {
	gout_file  << out_path.c_str() << "pptagu_" << keys[0] << ".dot";
	dump_stack_graph(gout_file.str(), unique, merged, keys, scaled_dist_map, cluster_map);
	gout_file.str("");
	gout_file  << out_path.c_str() << "pptagu_unscaled_" << keys[0] << ".dot";
	cerr << gout_file.str() << "\n";
	dump_stack_graph(gout_file.str(), unique, merged, keys, dist_map, cluster_map);
    }

    //
    // Merge the clustered tags together
    //
    id = 0;
    for (c = cluster_map.begin(); c != cluster_map.end(); c++) {

	tag = merge_tags(merged, c->second, id);
	tag->deleveraged = true;
	tag->masked      = true;
	//
	// Check to make sure the resulting tags aren't too far apart to merge. If so,
	// mask this tag off, so it is not considered in later merging steps.
	//
	if (!check_deleveraged_dist(dist_map, c->second, round))
	    tag->blacklisted = true;

	deleveraged_tags.push_back(tag);
	id++;
    }

    return 0;
}

int check_deleveraged_dist(map<int, map<int, double> > &dist_map, set<int> &merge_list, int max_dist) {
    set<int>::iterator i, j;
    double dist = 0;

    if (merge_list.size() == 1)
	return 1;

    for (i = merge_list.begin(); i != merge_list.end(); i++) {
	for (j = i, j++; j != merge_list.end(); j++) {
	    dist = dist_map[*i][*j];
	    //cerr << "Distance is: " << dist << " (between " << *i << " and " << *j << ")\n";
	    if (dist > max_dist)
		return 0;
	}
    }

    return 1;
}

int calc_kmer_distance(map<int, MergedStack *> &merged, int utag_dist) {
    //
    // Calculate the distance (number of mismatches) between each pair
    // of Radtags. We expect all radtags to be the same length;
    //
    KmerHashMap kmer_map;
    map<int, MergedStack *>::iterator it;
    MergedStack *tag_1, *tag_2;
    int i, j;

    // OpenMP can't parallelize random access iterators, so we convert
    // our map to a vector of integer keys.
    vector<int> keys;
    for (it = merged.begin(); it != merged.end(); it++) 
	keys.push_back(it->first);

    //
    // Calculate the number of k-mers we will generate. If kmer_len == 0,
    // determine the optimal length for k-mers.
    //
    int con_len   = strlen(merged[keys[0]]->con);
    if (set_kmer_len)
        kmer_len = determine_kmer_length(con_len, utag_dist);
    int num_kmers = con_len - kmer_len + 1;

    cerr << "  Number of kmers per sequence: " << num_kmers << "\n";

    //
    // Calculate the minimum number of matching k-mers required for a possible sequence match.
    //
    int min_hits = calc_min_kmer_matches(kmer_len, utag_dist, con_len, true);

    populate_kmer_hash(merged, kmer_map, kmer_len);
 
    #pragma omp parallel private(i, j, tag_1, tag_2)
    { 
        #pragma omp for schedule(dynamic) 
        for (i = 0; i < (int) keys.size(); i++) {
            tag_1 = merged[keys[i]];

            // Don't compute distances for masked tags
            if (tag_1->masked) continue;

            map<int, int> hits;
            vector<int>::iterator map_it;
            int d;
            //
            // Lookup the occurances of each k-mer in the kmer_map
            //
            for (j = 0; j < num_kmers; j++) {

                for (map_it  = kmer_map[tag_1->kmers[j]].begin(); 
                     map_it != kmer_map[tag_1->kmers[j]].end(); 
                     map_it++)
                    hits[*map_it]++;
            }

            //cerr << "  Tag " << tag_1->id << " hit " << hits.size() << " kmers.\n";

            //
            // Iterate through the list of hits. For each hit that has more than min_hits
            // check its full length to verify a match.
            //
            map<int, int>::iterator hit_it;
            for (hit_it = hits.begin(); hit_it != hits.end(); hit_it++) {
                //cerr << "  Tag " << hit_it->first << " has " << hit_it->second << " hits (min hits: " << min_hits << ")\n";

                if (hit_it->second < min_hits) continue;

                //cerr << "  Match found, checking full-length match\n";

                tag_2 = merged[hit_it->first];

                // Don't compute distances for masked tags
                if (tag_2->masked) continue;

                // Don't compare tag_1 against itself.
                if (tag_1 == tag_2) continue;

                d = dist(tag_1, tag_2);
                //cerr << "    Distance: " << d << "\n";

                //
                // Store the distance between these two sequences if it is
                // below the maximum distance (which governs those
                // sequences to be merged in the following step of the
                // algorithm.)
                //
                if (d == utag_dist)
                    tag_1->add_dist(tag_2->id, d);
            }
            // Sort the vector of distances.
            sort(tag_1->dist.begin(), tag_1->dist.end(), compare_dist);
        }
    }

    return 0;
}

int calc_distance(map<int, MergedStack *> &merged, int utag_dist) {
    //
    // Calculate the distance (number of mismatches) between each pair
    // of Radtags. We expect all radtags to be the same length;
    //
    map<int, MergedStack *>::iterator it;
    MergedStack *tag_1, *tag_2;
    int i, j;

    // OpenMP can't parallelize random access iterators, so we convert
    // our map to a vector of integer keys.
    vector<int> keys;
    for (it = merged.begin(); it != merged.end(); it++) 
	keys.push_back(it->first);

    #pragma omp parallel private(i, j, tag_1, tag_2)
    { 
        #pragma omp for schedule(dynamic) 
	for (i = 0; i < (int) keys.size(); i++) {

	    tag_1 = merged[keys[i]];

	    // Don't compute distances for masked tags
	    if (tag_1->masked) continue;

	    int d;

	    for (j = 0; j < (int) keys.size(); j++) {
		tag_2 = merged[keys[j]];

		// Don't compute distances for masked tags
		if (tag_2->masked) continue;

		// Don't compare tag_1 against itself.
		if (tag_1 == tag_2) continue;

		d = dist(tag_1, tag_2);
                //cerr << "  Distance: " << d << "\n";

		//
		// Store the distance between these two sequences if it is
		// below the maximum distance (which governs those
		// sequences to be merged in the following step of the
		// algorithm.)
		//
		if (d == utag_dist) {
		    tag_1->add_dist(tag_2->id, d);
                    //cerr << "  HIT.\n";
		}
	    }

	    // Sort the vector of distances.
	    sort(tag_1->dist.begin(), tag_1->dist.end(), compare_dist);
	}
    }

    return 0;
}

int reduce_radtags(HashMap &radtags, map<int, Stack *> &unique, map<int, Seq *> &rem) {
    HashMap::iterator it;
    vector<SeqId *>::iterator fit;

    Seq  *s;
    Stack *u;
    int   global_id = 1;

    for (it = radtags.begin(); it != radtags.end(); it++) {
	if (it->second.count < min_merge_cov) {
	    //
	    // Don't record this unique RAD-Tag if its coverage is below
	    // the specified cutoff. However, add the reads to the remainder
	    // vector for later processing.
	    //
	    for (fit = (*it).second.id.begin(); fit != (*it).second.id.end(); fit++) {
		s = new Seq((*fit)->id, (*it).first);
		rem[global_id] = s;
		global_id++;
	    }
	} else if (it->second.count > 1) {
	    //
	    // Populate a Stack object for this unique radtag. Create a
	    // map of the IDs for the sequences that have been
	    // collapsed into this radtag.
	    //
	    u         = new Stack;
	    u->id    = global_id;
	    u->count  = (*it).second.count;
	    u->add_seq((*it).first);

	    // Copy the original Fastq IDs from which this unique radtag was built.
	    for (fit = (*it).second.id.begin(); fit != (*it).second.id.end(); fit++) {
		u->add_id((*fit)->id);
	    }

	    unique[u->id] = u;
	    global_id++;
	}
    }

    if (unique.size() == 0) {
        cerr << "Error: Unable to form any stacks, data appear to be unique.\n";
        exit(1);
    }

    return 0;
}

int calc_coverage_distribution(map<int, Stack *> &unique, double &mean, double &stdev) {
    map<int, Stack *>::iterator i;
    double m     = 0.0;
    double s     = 0.0;
    double sum   = 0.0;
    uint   max   = 0;
    double total = 0.0;

    map<int, int> depth_dist;
    map<int, int>::iterator j;

    for (i = unique.begin(); i != unique.end(); i++) {
	m += i->second->count;
	total++;

        depth_dist[i->second->count]++;

	if (i->second->count > max)
	    max = i->second->count;
    }

    mean = round(m / total);

    //
    // Calculate the standard deviation
    //
    total = 0.0;

    for (i = unique.begin(); i != unique.end(); i++) {
	total++;
	s = i->second->count;
	sum += pow((s - mean), 2);
    }

    stdev = sqrt(sum / (total - 1));

    cerr << "  Mean coverage depth is " << mean << "; Std Dev: " << stdev << " Max: " << max << "\n";

    //
    // Output the distribution of stack depths
    //
    //for (j = depth_dist.begin(); j != depth_dist.end(); j++)
    //    cerr << j->first << "\t" << j->second << "\n";

    return 0;
}

double calc_merged_coverage_distribution(map<int, Stack *> &unique, map<int, MergedStack *> &merged) {
    map<int, MergedStack *>::iterator it;
    vector<int>::iterator            k;
    Stack *tag;
    double m     = 0.0;
    double s     = 0.0;
    double sum   = 0.0;
    double mean  = 0.0;
    double max   = 0.0;
    double stdev = 0.0;

    for (it = merged.begin(); it != merged.end(); it++) {
	m = 0.0;
	for (k = it->second->utags.begin(); k != it->second->utags.end(); k++) {
	    tag  = unique[*k];
	    m   += tag->count;
	}
	if (m > max) max = m;

	sum += m;
    }

    mean = sum / (double) merged.size();

    //
    // Calculate the standard deviation
    //
    for (it = merged.begin(); it != merged.end(); it++) {
	s = 0.0;
	for (k = it->second->utags.begin(); k != it->second->utags.end(); k++) {
	    tag  = unique[*k];
	    s   += tag->count;
	}
	sum += pow((s - mean), 2);
    }

    stdev = sqrt(sum / (merged.size() - 1));

    cerr << "  Mean merged coverage depth is " << mean << "; Std Dev: " << stdev << "; Max: " << max << "\n";

    return m;
}

int count_raw_reads(map<int, Stack *> &unique, map<int, MergedStack *> &merged) {
    map<int, MergedStack *>::iterator it;
    map<int, Stack *>::iterator sit;
    vector<int>::iterator k;
    Stack *tag;
    long int m = 0;

    map<int, int> uniq_ids;
    map<int, int>::iterator uit;

    for (it = merged.begin(); it != merged.end(); it++) {
	for (k = it->second->utags.begin(); k != it->second->utags.end(); k++) {
	    tag  = unique[*k];
	    m   += tag->count;

            if (uniq_ids.count(*k) == 0)
               uniq_ids[*k] = 0;
            uniq_ids[*k]++;
	}
        m += it->second->remtags.size();
    }

    for (uit = uniq_ids.begin(); uit != uniq_ids.end(); uit++)
       if (uit->second > 1)
           cerr << "  Unique stack #" << uit->first << " appears in " << uit->second << " merged stacks.\n";

    cerr << "Number of utilized reads: " << m << "\n";

    //for (sit = unique.begin(); sit != unique.end(); sit++)
    //   if (uniq_ids.count(sit->first) == 0)
    //       cerr << "  Stack " << sit->first << ": '" << sit->second->seq << "' unused.\n";

    return 0;
}

int write_results(map<int, MergedStack *> &m, map<int, Stack *> &u, map<int, Seq *> &r) {
    map<int, MergedStack *>::iterator i;
    vector<SeqId *>::iterator  j;
    vector<int>::iterator      k;
    vector<SNP *>::iterator    s;
    map<string, int>::iterator t;
    MergedStack *tag_1;
    Stack *tag_2;
    Seq  *rem;

    //
    // Parse the input file name to create the output files
    //
    size_t pos_1 = in_file.find_last_of("/");
    size_t pos_2 = in_file.find_last_of(".");
    string tag_file = out_path + in_file.substr(pos_1 + 1, (pos_2 - pos_1 - 1)) + ".tags.tsv";
    string snp_file = out_path + in_file.substr(pos_1 + 1, (pos_2 - pos_1 - 1)) + ".snps.tsv";
    string all_file = out_path + in_file.substr(pos_1 + 1, (pos_2 - pos_1 - 1)) + ".alleles.tsv";

    // Open the output files for writing.
    std::ofstream tags(tag_file.c_str());
    std::ofstream snps(snp_file.c_str());
    std::ofstream alle(all_file.c_str());
    int id;

    for (i = m.begin(); i != m.end(); i++) {
	float total = 0;
	tag_1 = i->second;

	// First write the consensus sequence
	tags << "0" << "\t" 
	     << sql_id << "\t" 
	     << tag_1->id << "\t" 
             << tag_1->loc.chr << "\t"
             << tag_1->loc.bp << "\t"
	     << "consensus\t" << "\t\t" 
	     << tag_1->con << "\t" 
	     << tag_1->deleveraged << "\t" 
	     << tag_1->blacklisted << "\t"
	     << tag_1->lumberjackstack << "\n";

	// Now write out the components of each unique tag merged into this one.
	id = 0;
	for (k = tag_1->utags.begin(); k != tag_1->utags.end(); k++) {
	    tag_2  = u[*k];
	    total += tag_2->count;

	    for (j = tag_2->map.begin(); j != tag_2->map.end(); j++) {
		tags << "0" << "\t" << sql_id << "\t" << tag_1->id << "\t\t\t" << "primary\t" << id << "\t" << (*j)->id << "\t" << tag_2->seq << "\t\t\t\n";
	    }

	    id++;
	}

	//
	// Write out the remainder tags merged into this unique tag.
	//
	total += tag_1->remtags.size();
	for (k = tag_1->remtags.begin(); k != tag_1->remtags.end(); k++) {
	    rem = r[*k];
	    tags << "0" << "\t" << sql_id << "\t" << tag_1->id << "\t\t\t" << "secondary\t\t" << rem->id << "\t" << rem->seq << "\t\t\t\n";
	}

	// Write out any SNPs detected in this unique tag.
	for (s = tag_1->snps.begin(); s != tag_1->snps.end(); s++) {
	    snps << "0" << "\t" << sql_id << "\t" << tag_1->id << "\t" << (*s)->col << "\t" << (*s)->lratio << "\t" << (*s)->rank_1 << "\t" << (*s)->rank_2 << "\n";
	}

	// Write the expressed alleles seen for the recorded SNPs and
	// the percentage of tags a particular allele occupies.
	for (t = tag_1->alleles.begin(); t != tag_1->alleles.end(); t++) {
	    alle << "0" << "\t" << sql_id << "\t" << tag_1->id << "\t" << (*t).first << "\t" << (((*t).second/total) * 100) << "\t" << (*t).second << "\n";
	}
    }

    tags.close();
    snps.close();
    alle.close();

    //
    // If specified, output reads not utilized in any stacks.
    //
    if (retain_rem_reads) {
	string unused_file = out_path + in_file.substr(pos_1 + 1, (pos_2 - pos_1 - 1)) + ".unused.fa";

	std::ofstream unused(unused_file.c_str());

	map<int, Seq *>::iterator r_it;
	for (r_it = r.begin(); r_it != r.end(); r_it++)
	    if (r_it->second->utilized == false)
		unused << ">" << r_it->second->id << "\n" << r_it->second->seq << "\n";

	unused.close();
    }

    return 0;
}

int hclust(vector<int> &keys, map<int, map<int, double> > &dist_map, int num_clusters, map<int, set<int> > &cluster_map) {
    uint s, t;

    //for (s = 0; s < keys.size(); s++) {
    //cerr << "Key: " << keys[s] << "\n";
    //}

    //
    // Generate a two-dimensional, ragged array containing the distances
    // between stacks.
    // From the cluster.h software manual:
    //  This matrix contains all the distances between the items that
    //  are being clustered.  The distance matrix can therefore be
    //  stored as a ragged array, with the number of columns in each
    //  row equal to the (zero-offset) row number. The distance
    //  between items i and j is stored in location [i][j] if j < i,
    //  in [j][i] if j > i, while it is zero if j=i. Note that the
    //  first row of the distance matrix is empty. It is included for
    //  computational convenience, as including an empty row requires
    //  minimal storage.
    //
    map<int, int> key_map;
    double **datamatrix = new double*[keys.size()];
    for (s = 1; s < keys.size(); s++)
	datamatrix[s] = new double[s];

    for (s = 0; s < keys.size(); s++) {
	key_map[s] = keys[s];

	for (t = 0; t < s; t++) {
	    datamatrix[s][t] = dist_map[keys[s]][keys[t]];
	    //cerr << datamatrix[s][t] << "\t";
	}
	//cerr << "\n";
    }

    //
    // Hierarchically cluster the data points, using single-linkage clustering.
    //
    Node *tree = treecluster(keys.size(), keys.size(), NULL, NULL, NULL, 0, 0, 's', datamatrix);

    // Free the data matrix as it has been modified by the clustering function.
    for (s = 1; s < keys.size(); s++) 
	delete datamatrix[s];
    delete datamatrix;

    // Indication that the treecluster routine failed
    if (tree == NULL) {
	cerr << "Tree clustering algorithm failed.\n";
	return 0;
    }

    int *clusterid = new int[keys.size()];

    cuttree(keys.size(), tree, num_clusters, clusterid);

    //
    // Record the stack clusters
    //
    for (s = 0; s < keys.size(); s++) {
	cluster_map[clusterid[s]].insert(key_map[s]);
	//cerr << "Stack " << s << " [uid: " << key_map[s] << "] is in cluster: " << clusterid[s] << "\n";
    }

    delete tree;
    delete clusterid;

    return 0;
}

int dump_stack_graph(string data_file, 
		     map<int, Stack *> &unique, 
		     map<int, MergedStack *> &merged, 
		     vector<int> &keys, 
		     map<int, map<int, double> > &dist_map, 
		     map<int, set<int> > &cluster_map) {
    uint s, t;
    double d, scale, scaled_d;
    char label[32];
    vector<string> colors;
    std::ofstream data(data_file.c_str());

    size_t pos_1 = data_file.find_last_of("/");
    size_t pos_2 = data_file.find_last_of(".");
    string title = data_file.substr(pos_1 + 1, pos_2 - pos_1 - 1);

    //
    // Output a list of IDs so we can locate these stacks in the final results.
    //
    for (s = 0; s < keys.size(); s++)
	data << "/* " << keys[s] << ": " << unique[merged[keys[s]]->utags[0]]->map[0]->id << "; depth: " << merged[keys[s]]->count << " */\n";

    //
    // Output a specification to visualize the stack graph using graphviz:
    //   http://www.graphviz.org/
    //
    data << "graph " << title.c_str() << " {\n"
	 << "rankdir=LR\n"
	 << "size=\"20!\"\n"
	 << "overlap=false\n"
	 << "node [shape=circle style=filled fillcolor=\"#3875d7\" fontname=\"Arial\"];\n"
	 << "edge [fontsize=8.0 fontname=\"Arial\" color=\"#aaaaaa\"];\n";

    colors.push_back("red");
    colors.push_back("blue"); 
    colors.push_back("green"); 
    colors.push_back("brown"); 
    colors.push_back("purple");

    map<int, set<int> >::iterator c;
    set<int>::iterator it;
    int color_index = 0;
    string color;

    // Write out the clusters created by R, prior to writing all the nodes and connections.
    s = 0;
    for (c = cluster_map.begin(); c != cluster_map.end(); c++) {
	data << "subgraph " << s << " {\n"
	     << "    edge [penwidth=5 fontsize=12.0 fontcolor=\"black\" color=\"black\"]\n";

	if ((*c).second.size() == 1) {
	    color = "white";
	    data << "    node [fillcolor=" << color.c_str() << " fontcolor=\"black\"]\n";
	} else {
	    color = colors[color_index % colors.size()];
	    data << "    node [fillcolor=" << color.c_str() << " fontcolor=\"white\"]\n";
	    color_index++;
	}

	for (it = (*c).second.begin(); it != (*c).second.end(); it++) {
	    data << "    " << *it << "\n";
	}

	if ((*c).second.size() > 1) {
	    uint j = 0;
	    for (it = (*c).second.begin(); it != (*c).second.end(); it++) {
		data << *it;
		if (j < (*c).second.size() - 1)
		    data << " -- ";
		j++;
	    }
	}

	data << "}\n";
	s++;
    }

    //
    // Scale the graph to display on a 10 inch canvas. Find the largest edge weight
    // and scale the edge lengths to fit the canvas.
    //
    for (s = 0; s < keys.size(); s++)
	for (t = s+1; t < keys.size(); t++)
	    scale = dist_map[keys[s]][keys[t]] > scale ? dist_map[keys[s]][keys[t]] : scale;
    scale = scale / 20;

    for (s = 0; s < keys.size(); s++) {
	for (t = s+1; t < keys.size(); t++) {
	    d = dist_map[keys[s]][keys[t]];
	    scaled_d = d / scale;
	    scaled_d = scaled_d < 0.75 ? 0.75 : scaled_d;
	    sprintf(label, "%.1f", d);
	    data << keys[s] << " -- " << keys[t] << " [len=" << scaled_d << ", label=" << label << "];\n";
	}
    }

    data << "}\n";

    data.close();

    return 0;
}

int dump_unique_tags(map<int, Stack *> &u) {
    map<int, Stack *>::iterator it;
    vector<SeqId *>::iterator fit;
    vector<pair<int, int> >::iterator pit;
    vector<int>::iterator mit;

    for (it = u.begin(); it != u.end(); it++) {

	cerr << "UniqueTag UID: " << (*it).second->id << "\n"
	     << "  Seq:       "   << (*it).second->seq << "\n"
	     << "  IDs:       "; 

	for (fit = (*it).second->map.begin(); fit != (*it).second->map.end(); fit++)
	    cerr << (*fit)->id << " ";

	cerr << "\n\n";
    }

    return 0;
}

int dump_merged_tags(map<int, MergedStack *> &m) {
    map<int, MergedStack *>::iterator it;
    vector<pair<int, int> >::iterator pit;
    vector<int>::iterator fit;

    for (it = m.begin(); it != m.end(); it++) {

	cerr << "MergedStack ID: " << it->second->id << "\n"
	     << "  Consensus:  ";
	if (it->second->con != NULL)
	    cerr << it->second->con << "\n";
	else 
	    cerr << "\n";
	cerr << "  IDs:        "; 

	for (fit = it->second->utags.begin(); fit != it->second->utags.end(); fit++)
	    cerr << (*fit) << " ";

	cerr << "\n"
	     << "  Distances: ";

	for (pit = it->second->dist.begin(); pit != it->second->dist.end(); pit++)
	    cerr << (*pit).first << ": " << (*pit).second << ", ";

	cerr << "\n\n";
    }

    return 0;
}

int load_radtags(string in_file, HashMap &radtags) {
    Input *fh;
    Seq   *c;

    if (in_file_type == fasta)
        fh = new Fasta(in_file.c_str());
    else if (in_file_type == fastq)
        fh = new Fastq(in_file.c_str());

    cerr << "  Parsing " << in_file.c_str() << "\n";
    long int corrected = 0;
    long int i = 0;
    while ((c = fh->next_seq()) != NULL) {
        if (i % 10000 == 0) cerr << "Loading RAD-Tag " << i << "       \r";

	for (char *p = c->seq; *p != '\0'; p++)
	    switch (*p) {
	    case 'N':
	    case 'n':
	    case '.':
		*p = 'A';
		corrected++;
	    }

	radtags[c->seq].add_id(c->id);
	radtags[c->seq].count++;
        i++;
    }
    cerr << "Loaded " << i << " RAD-Tags; inserted " << radtags.size() << " elements into the RAD-Tags hash map.\n";

    if (i == 0) {
        cerr << "Error: Unable to load data from '" << in_file.c_str() << "'.\n";
        exit(1);
    }

    cerr << corrected << " reads contained uncalled nucleotides that were modified.\n";

    //
    // Check to make sure all the reads are of the same length.
    //
    HashMap::iterator it;
    uint len = strlen(radtags.begin()->first);
    for (it = radtags.begin(); it != radtags.end(); it++)
        if (strlen((*it).first) != len)
            cerr << "  Warning: '" << (*it).second.id[0]->id << "' has a different length.\n";

    //
    // Close the file and delete the Input object.
    //
    delete fh;

    return 0;
}

int calc_triggers(double cov_mean, double cov_stdev, int &deleverage_trigger, int &removal_trigger) {

    deleverage_trigger = (int) round(cov_mean + cov_stdev * cov_scale);
    removal_trigger    = (int) round(cov_mean + (cov_stdev * 2) * cov_scale);

    return 0;

//     //
//     // Calculate the deleverage trigger. Assume RAD-Tags are selected from
//     // the sample for sequencing randomly, forming a poisson distribution
//     // representing the depths of coverage of RAD-Tags in the sample. Calculate 
//     // the trigger value that is larger than the depth of coverage of 99.9999% of stacks.
//     //
//     long double lambda = cov_mean;
//     int k         = 0;
//     long double d = 0.0;
//     long double e = 0.0;
//     long double f = 0.0;
//     long double g = 0.0;
//     long double h = 0.0;
//     long double i = 0.0;

//     do {
// 	e  = exp(-1 * lambda);
// 	g  = pow(lambda, k);
// 	f  = factorial(k);
// 	h  = (e * g);
// 	i  = h / f;
// 	d += i;

// 	//cerr << "iteration " << k << "; e: " << e << " h: " << h << " g: " << g << " F: " << f << " i: " << i << " D: " << d << "\n";
// 	k++;
//     } while (d < 0.999999);

//     return k - 1;
}

long double factorial(int i) {
    long double f = 1;

    if (i == 0) return 1;

    do {
	f = f * i;
	i--;
    } while (i > 0);

    return f;
}

int parse_command_line(int argc, char* argv[]) {
    int c;
     
    while (1) {
	static struct option long_options[] = {
	    {"help",        no_argument,       NULL, 'h'},
            {"version",     no_argument,       NULL, 'v'},
	    {"infile_type", required_argument, NULL, 't'},
	    {"file",        required_argument, NULL, 'f'},
	    {"outpath",     required_argument, NULL, 'o'},
	    {"id",          required_argument, NULL, 'i'},
	    {"min_cov",     required_argument, NULL, 'm'},
	    {"max_dist",    required_argument, NULL, 'M'},
	    {"num_threads", required_argument, NULL, 'p'},
	    {"deleverage",  no_argument,       NULL, 'd'},
	    {"remove_rep",  no_argument,       NULL, 'r'},
	    {"retain_rem",  no_argument,       NULL, 'R'},
	    {"graph",       no_argument,       NULL, 'g'},
	    {"exp_cov",     no_argument,       NULL, 'E'},
	    {"cov_stdev",   no_argument,       NULL, 's'},
	    {"cov_scale",   no_argument,       NULL, 'S'},
	    {"model_type",  required_argument, NULL, 'T'},
	    {"bc_err_freq", required_argument, NULL, 'e'},
	    {0, 0, 0, 0}
	};
	
	// getopt_long stores the option index here.
	int option_index = 0;
     
	c = getopt_long(argc, argv, "hRvdrgf:o:i:m:e:E:s:S:p:t:M:T:", long_options, &option_index);
     
	// Detect the end of the options.
	if (c == -1)
	    break;
     
	switch (c) {
	case 'h':
	    help();
	    break;
     	case 't':
            if (strcmp(optarg, "tsv") == 0)
                in_file_type = tsv;
            else if (strcmp(optarg, "bowtie") == 0)
                in_file_type = bowtie;
            else if (strcmp(optarg, "sam") == 0)
                in_file_type = sam;
            else if (strcmp(optarg, "fasta") == 0)
                in_file_type = fasta;
            else if (strcmp(optarg, "fastq") == 0)
                in_file_type = fastq;
            else
                in_file_type = unknown;
	    break;
     	case 'f':
	    in_file = optarg;
	    break;
	case 'o':
	    out_path = optarg;
	    break;
	case 'i':
	    sql_id = atoi(optarg);
	    break;
	case 'm':
	    min_merge_cov = atoi(optarg);
	    break;
	case 'M':
	    max_utag_dist = atoi(optarg);
	    break;
	case 'd':
	    deleverage_stacks++;
	    break;
	case 'r':
	    remove_rep_stacks++;
	    break;
	case 'R':
	    retain_rem_reads = true;
	    break;
	case 'g':
	    dump_graph++;
	    break;
	case 'E':
	    cov_mean = atof(optarg);
	    break;
	case 's':
	    cov_stdev = atof(optarg);
	    break;
	case 'S':
	    cov_scale = atof(optarg);
	    break;
     	case 'T':
            if (strcmp(optarg, "snp") == 0) {
                model_type = snp;
            } else if (strcmp(optarg, "fixed") == 0) {
                model_type = fixed;
            } else {
                cerr << "Unknown model type specified '" << optarg << "'\n";
                help();
            }
	case 'e':
	    barcode_err_freq = atof(optarg);
	    break;
	case 'p':
	    num_threads = atoi(optarg);
	    break;
        case 'v':
            version();
            break;
	case '?':
	    // getopt_long already printed an error message.
	    help();
	    break;
     
	default:
	    cerr << "Unknown command line option '" << (char) c << "'\n";
	    help();
	    abort();
	}
    }

    if (in_file.length() == 0 || in_file_type == unknown) {
	cerr << "You must specify an input file of a supported type.\n";
	help();
    }

    if (out_path.length() == 0) 
	out_path = ".";

    if (out_path.at(out_path.length() - 1) != '/') 
	out_path += "/";

    if (model_type == fixed && barcode_err_freq == 0) {
	cerr << "You must specify the barcode error frequency.\n";
	help();
    }

    return 0;
}

void version() {
    std::cerr << "ustacks " << VERSION << "\n\n";

    exit(0);
}

void help() {
    std::cerr << "ustacks " << VERSION << "\n"
              << "ustacks -t file_type -f file_path [-d] [-r] [-o path] [-i id] [-e errfreq] [-m min_cov] [-M max_dist] [-p num_threads] [-R] [-h]" << "\n"
              << "  p: enable parallel execution with num_threads threads.\n"
	      << "  t: input file Type. Supported types: fasta, fastq, bowtie, sam, tsv.\n"
              << "  f: input file path.\n"
	      << "  o: output path to write results." << "\n"
	      << "  i: SQL ID to insert into the output to identify this sample." << "\n"
	      << "  m: Minimum depth of coverage required to create a stack." << "\n"
	      << "  M: Maximum distance (in nucleotides) allowed between stacks." << "\n"
	      << "  d: enable the Deleveraging algorithm, used for resolving over merged tags." << "\n"
	      << "  r: enable the Removal algorithm, to drop highly-repetitive stacks (and nearby errors) from the algorithm." << "\n"
	      << "  e: specify the barcode error frequency (0 < e < 1) if using the 'fixed' model.\n"
	      << "  R: retain unused reads.\n"
	      << "  h: display this help messsage." << "\n\n";

    exit(0);
}
