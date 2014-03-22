// -*-mode:c++; c-style:k&r; c-basic-offset:4;-*-
//
// Copyright 2010-2013, Julian Catchen <jcatchen@uoregon.edu>
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
bool      call_sec_hapl     = true;
bool      set_kmer_len      = true;
int       kmer_len          = 0;
int       max_kmer_len      = 19;
int       min_merge_cov     = 3;
uint      max_subgraph      = 3;
int       dump_graph        = 0;
int       retain_rem_reads  = false;
int       deleverage_stacks = 0;
int       remove_rep_stacks = 0;
int       max_utag_dist     = 2;
int       max_rem_dist      = -1;
double    cov_mean          = 0.0;
double    cov_stdev         = 0.0;
double    cov_scale         = 1;
int       deleverage_trigger;
int       removal_trigger;
//
// For use with the multinomial model to call fixed nucleotides.
//
modelt model_type         = snp;
double alpha              = 0.05;
double bound_low          = 0.0;
double bound_high         = 1.0;
double p_freq             = 0.5;
double barcode_err_freq   = 0.0;
double heterozygote_limit = -3.84;
double homozygote_limit   =  3.84;

int main (int argc, char* argv[]) {

    parse_command_line(argc, argv);

    //
    // Set the max remainder distance to be greater than the max_utag_dist, if it is not
    // specified on the command line.
    //
    if (max_rem_dist == -1) max_rem_dist = max_utag_dist + 2;

    cerr << "Min depth of coverage to create a stack: " << min_merge_cov << "\n"
	 << "Max distance allowed between stacks: " << max_utag_dist << "\n"
	 << "Max distance allowed to align secondary reads: " << max_rem_dist << "\n"
	 << "Max number of stacks allowed per de novo locus: " << max_subgraph << "\n"
	 << "Deleveraging algorithm: " << (deleverage_stacks ? "enabled" : "disabled") << "\n"
	 << "Removal algorithm: " << (remove_rep_stacks ? "enabled" : "disabled") << "\n"
	 << "Model type: "; 
    switch (model_type) {
    case snp:
	cerr << "SNP\n";
	break;
    case fixed: 
	cerr << "Fixed\n";
	break;
    case bounded:
	cerr << "Bounded; lower epsilon bound: " << bound_low << "; upper bound: " << bound_high << "\n";
	break;
    }
    cerr << "Alpha significance level for model: " << alpha << "\n";

    //
    // Set limits to call het or homozygote according to chi-square distribution with one 
    // degree of freedom:
    //   http://en.wikipedia.org/wiki/Chi-squared_distribution#Table_of_.CF.872_value_vs_p-value
    //
    if (alpha == 0.1) {
	heterozygote_limit = -2.71;
	homozygote_limit   =  2.71;
    } else if (alpha == 0.05) {
	heterozygote_limit = -3.84;
	homozygote_limit   =  3.84;
    } else if (alpha == 0.01) {
	heterozygote_limit = -6.64;
	homozygote_limit   =  6.64;
    } else if (alpha == 0.001) {
	heterozygote_limit = -10.83;
	homozygote_limit   =  10.83;
    }

    //
    // Set the number of OpenMP parallel threads to execute.
    //
    #ifdef _OPENMP
    omp_set_num_threads(num_threads);
    #endif

    DNASeqHashMap     radtags;
    vector<DNASeq *>  radtags_keys;
    map<int, Rem *>   remainders;
    set<int>          merge_map;
    map<int, Stack *> unique;

    load_radtags(in_file, radtags, radtags_keys);

    reduce_radtags(radtags, unique, remainders);

    free_radtags_hash(radtags, radtags_keys);

    // dump_unique_tags(unique);

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

    cerr << "Calculating distance between stacks...\n";
    calc_kmer_distance(merged, max_utag_dist);

    cerr << "Merging stacks, maximum allowed distance: " << max_utag_dist << " nucleotide(s)\n";
    merge_stacks(unique, remainders, merged, merge_map, max_utag_dist);

    call_consensus(merged, unique, remainders, false);

    calc_merged_coverage_distribution(unique, merged);

    //dump_merged_tags(merged);

    cerr << "Merging remainder radtags\n";
    merge_remainders(merged, remainders);

    // Call the consensus sequence again, now that remainder tags have been merged.
    call_consensus(merged, unique, remainders, true);

    count_raw_reads(unique, remainders, merged);

    cerr << "Writing results\n";
    write_results(merged, unique, remainders);

    return 0;
}

int merge_remainders(map<int, MergedStack *> &merged, map<int, Rem *> &rem) {
    map<int, Rem *>::iterator it;
    int j, k;

    // OpenMP can't parallelize random access iterators, so we convert
    // our map to a vector of integer keys.
    vector<int> keys;
    uint tot = 0;
    for (it = rem.begin(); it != rem.end(); it++) {
    	keys.push_back(it->first);
	tot += it->second->count();
    }

    cerr << "  " << tot << " remainder sequences left to merge.\n";

    if (max_rem_dist <= 0) {
	cerr << "  Matched 0 remainder reads; unable to match " << tot << " remainder reads.\n";
	return 0;
    }

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

    KmerHashMap    kmer_map;
    vector<char *> kmer_map_keys;
    populate_kmer_hash(merged, kmer_map, kmer_map_keys, kmer_len);
    int utilized = 0;

    //
    // Create a character buffer to hold the Rem sequence, this is faster
    // than repeatedly decoding the DNASeq buffers.
    // 
    //it = rem.find(keys[0]);
    //char *buf = new char[it->second->seq->size + 1];

    #pragma omp parallel private(it, k)
    { 
        #pragma omp for schedule(dynamic) 
    	for (j = 0; j < (int) keys.size(); j++) {
    	    it = rem.find(keys[j]);
    	    Rem  *r   = it->second;
	    char *buf = new char[r->seq->size + 1];

            //
            // Generate the k-mers for this remainder sequence
            //
            vector<char *> rem_kmers;
            buf = r->seq->seq(buf);
            generate_kmers(buf, kmer_len, num_kmers, rem_kmers);

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

                int d = dist(merged[hit_it->first], buf);
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
                delete [] rem_kmers[k];

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

	    delete [] buf;

    	    // Found a merge partner.
    	    if (min_id >= 0 && count == 1) {
    		r->utilized = true;
                #pragma omp critical
    		{
    		    merged[min_id]->remtags.push_back(it->first);
    		    utilized += it->second->count();
    		}
    	    }
    	}
    }

    free_kmer_hash(kmer_map, kmer_map_keys);
    //delete [] buf;

    cerr << "  Matched " << utilized << " remainder reads; unable to match " << tot - utilized << " remainder reads.\n";

    return 0;
}

int call_alleles(MergedStack *mtag, vector<DNASeq *> &reads, vector<read_type> &read_types) {
    int     row;
    int     height = reads.size();
    string  allele;
    DNASeq *d;
    char    base;
    vector<SNP *>::iterator snp;

    for (row = 0; row < height; row++) {
	allele.clear();

	uint snp_cnt = 0;

	//
	// Only call a haplotype from primary reads.
	//
	if (!call_sec_hapl && read_types[row] == secondary) continue;

	for (snp = mtag->snps.begin(); snp != mtag->snps.end(); snp++) {
	    if ((*snp)->type != snp_type_het) continue;

	    snp_cnt++;

            d    = reads[row];
	    base = (*d)[(*snp)->col];

	    //
	    // Check to make sure the nucleotide at the location of this SNP is
	    // of one of the two possible states the multinomial model called.
	    //
	    if (base == (*snp)->rank_1 || base == (*snp)->rank_2) 
		allele += base;
	    else
		break;
	}

	if (snp_cnt > 0 && allele.length() == snp_cnt)
	    mtag->alleles[allele]++;
    }

    return 0;
}

int call_consensus(map<int, MergedStack *> &merged, map<int, Stack *> &unique, map<int, Rem *> &rem, bool invoke_model) {
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
	MergedStack *mtag;
	Stack       *utag;
	Rem         *r;

        #pragma omp for schedule(dynamic) 
    	for (i = 0; i < (int) keys.size(); i++) {
    	    mtag = merged[keys[i]];

    	    //
    	    // Create a two-dimensional array, each row containing one read. For
    	    // each unique tag that has been merged together, add the sequence for
    	    // that tag into our array as many times as it originally occurred. 
    	    //
    	    vector<int>::iterator j;
    	    vector<DNASeq *>  reads;
    	    vector<read_type> read_types;

    	    for (j = mtag->utags.begin(); j != mtag->utags.end(); j++) {
    		utag = unique[*j];

    		for (uint k = 0; k < utag->count(); k++) {
    		    reads.push_back(utag->seq);
    		    read_types.push_back(primary);
    		}
    	    }

    	    // For each remainder tag that has been merged into this Stack, add the sequence. 
    	    for (j = mtag->remtags.begin(); j != mtag->remtags.end(); j++) {
		r = rem[*j];

    		for (uint k = 0; k < r->count(); k++) {
		    reads.push_back(r->seq);
		    read_types.push_back(secondary);
		}
    	    }

    	    //
    	    // Iterate over each column of the array and call the consensus base.
    	    //
    	    int row, col;
    	    int length = reads[0]->size;
    	    int height = reads.size();
    	    string con;
    	    map<char, int> nuc;
    	    map<char, int>::iterator max, n;
            DNASeq *d;

    	    for (col = 0; col < length; col++) {
    		nuc['A'] = 0; 
    		nuc['C'] = 0;
    		nuc['G'] = 0;
    		nuc['T'] = 0;

    		for (row = 0; row < height; row++) {
    		    d = reads[row];
		    if (nuc.count((*d)[col]))
			nuc[(*d)[col]]++;
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

		//
    		// Search this column for the presence of a SNP
		//
    		if (invoke_model) 
		    switch(model_type) {
		    case snp:
                        call_multinomial_snp(mtag, col, nuc, true);
			break;
		    case bounded:
			call_bounded_multinomial_snp(mtag, col, nuc, true);
			break;
		    case fixed:
                        call_multinomial_fixed(mtag, col, nuc);
			break;
		    }
    	    }

    	    if (invoke_model) {
    		call_alleles(mtag, reads, read_types);

                if (model_type == fixed) {
                    //
                    // Mask nucleotides that are not fixed.
                    //
                    vector<SNP *>::iterator s;
                    for (s = mtag->snps.begin(); s != mtag->snps.end(); s++) {
			if ((*s)->type == snp_type_unk)
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
    Stack       *utag;
    MergedStack *mtag;
    int          k = 0;

    it_old = merged.begin();

    for (i = unique.begin(); i != unique.end(); i++) {
	utag = (*i).second;
	mtag = new MergedStack;

	mtag->id    = k;
	mtag->count = utag->count();
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

int merge_stacks(map<int, Stack *> &unique, map<int, Rem *> &rem, map<int, MergedStack *> &merged, set<int> &merge_map, int round) {
    map<int, MergedStack *> new_merged;
    map<int, MergedStack *>::iterator it, it_old, it_new;
    MergedStack *tag_1, *tag_2;
    vector<set<int> > merge_lists;

    uint index     = 0;
    int  cohort_id  = 0;
    int  id         = 1;
    uint delev_cnt = 0;
    uint blist_cnt = 0;

    for (it = merged.begin(); it != merged.end(); it++) {
	tag_1 = it->second;

	//
	// This tag may already have been merged by an earlier operation.
	//
	if (merge_map.count(tag_1->id) > 0)
	    continue;

	queue<int>                        merge_list;
	pair<set<int>::iterator,bool>     ret;
	vector<pair<int, int> >::iterator k;

	merge_lists.push_back(set<int>());

	if (tag_1->masked) {
	    merge_lists[index].insert(tag_1->id);
	    index++;
	    continue;
	}

	//
	// Construct a list of MergedStacks that are within a particular distance
	// of this tag.
	//
	merge_lists[index].insert(tag_1->id);
	merge_list.push(tag_1->id);

	while (!merge_list.empty()) {
	    tag_2 = merged[merge_list.front()];
	    merge_list.pop();

	    for (k = tag_2->dist.begin(); k != tag_2->dist.end(); k++) {
		ret = merge_lists[index].insert(k->first);

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
	// Record the nodes that have been merged in this round.
	//
	set<int>::iterator j;
	for (j = merge_lists[index].begin(); j != merge_lists[index].end(); j++)
	    merge_map.insert(*j);

	index++;
    }

    #pragma omp parallel private(tag_1, tag_2)
    {
	vector<MergedStack *> merged_tags;

        #pragma omp for reduction(+:delev_cnt) reduction(+:blist_cnt)
	for (uint index = 0; index < merge_lists.size(); index++) {
	    //
	    // Deal with the simple case of a single locus that does not need to be merged.
	    //
	    if (merge_lists[index].size() == 1) {
		tag_1 = merged[*(merge_lists[index].begin())];
		tag_2 = merge_tags(merged, merge_lists[index], 0);

		//
		// If this tag is masked, keep the old cohort_id.
		//
		if (tag_1->masked) {
		    tag_2->cohort_id = tag_1->cohort_id;
		} else {
		    tag_2->cohort_id = cohort_id;
		    #pragma omp atomic
		    cohort_id++;
		}
		merged_tags.push_back(tag_2);
		continue;
	    }

	    //
	    // Break large loci down by constructing a minimum
	    // spanning tree and severing long distance edges.
	    //
	    if (deleverage_stacks) {
		vector<MergedStack *> tags;
		bool delev;

		deleverage(unique, rem, merged, merge_lists[index], cohort_id, tags);

		if (tags.size() == 1) {
		    delev = false;
		} else {
		    delev_cnt++;
		    delev = true;
		}

		for (uint t = 0; t < tags.size(); t++) {
		    //tags[t]->id = id;
		    tags[t]->deleveraged = delev;

		    if (tags[t]->utags.size() > max_subgraph) {
			tags[t]->masked      = true;
			tags[t]->blacklisted = true;
			blist_cnt++;
		    }

		    //new_merged.insert(pair<int, MergedStack *>(id, tags[t]));
		    merged_tags.push_back(tags[t]);
		    //id++;
		}

		#pragma omp atomic
		cohort_id++;

	    } else {
		//
		// If not deleveraging, merge these tags together into a new MergedStack object.
		//
		tag_2 = merge_tags(merged, merge_lists[index], 0);
		tag_2->cohort_id = cohort_id;

		if (tag_2->utags.size() > max_subgraph) {
		    tag_2->masked      = true;
		    tag_2->blacklisted = true;
		    blist_cnt++;
		}

		//new_merged.insert(pair<int, MergedStack *>(id, tag_2));
		merged_tags.push_back(tag_2);

                #pragma omp atomic
		cohort_id++;
		//id++;
	    }
	}

	//
	// Merge the accumulated tags into the new_merged map.
	//
        #pragma omp critical
	{
	    it_old = merged.begin();
	    for (uint j = 0; j < merged_tags.size(); j++) {
		merged_tags[j]->id = id;
		it_new = new_merged.insert(it_old, pair<int, MergedStack *>(id, merged_tags[j]));
		it_old = it_new;
		id++;
	    }
	}
    }

    uint new_cnt = new_merged.size();
    uint old_cnt = merged.size();

    //
    // Free the memory from the old map of merged tags.
    //
    for (it = merged.begin(); it != merged.end(); it++)
	delete it->second;

    merged = new_merged;

    cerr << "  " << old_cnt << " stacks merged into " << new_cnt 
	 << " stacks; deleveraged " << delev_cnt 
	 << " stacks; removed " << blist_cnt << " stacks.\n";

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

	for (j = tag_2->remtags.begin(); j != tag_2->remtags.end(); j++)
	    tag_1->remtags.push_back(*j);

        tag_1->count += tag_2->count;
    }

    return tag_1;
}

MergedStack *merge_tags(map<int, MergedStack *> &merged, int *merge_list, int merge_list_size, int id) {
    int                   i;
    vector<int>::iterator j;
    MergedStack *tag_1, *tag_2;

    tag_1     = new MergedStack;
    tag_1->id = id;

    for (i = 0; i < merge_list_size; i++) {
	tag_2 = merged[merge_list[i]];

	tag_1->deleveraged     = tag_2->deleveraged     ? true : tag_1->deleveraged;
	tag_1->masked          = tag_2->masked          ? true : tag_1->masked;
	tag_1->blacklisted     = tag_2->blacklisted     ? true : tag_1->blacklisted;
	tag_1->lumberjackstack = tag_2->lumberjackstack ? true : tag_1->lumberjackstack;

	for (j = tag_2->utags.begin(); j != tag_2->utags.end(); j++)
	    tag_1->utags.push_back(*j);

	for (j = tag_2->remtags.begin(); j != tag_2->remtags.end(); j++)
	    tag_1->remtags.push_back(*j);

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
	    tag_1->masked          = true;
	    tag_1->blacklisted     = true;

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

    //
    // Free the memory from the old map of merged tags.
    //
    map<int, MergedStack *>::iterator it;
    for (it = merged.begin(); it != merged.end(); it++)
	delete it->second;

    merged = new_merged;

    cerr << "  " << merged.size() << " stacks remain for merging.\n";

    return 0;
}

int deleverage(map<int, Stack *> &unique, 
               map<int, Rem *> &rem,
	       map<int, MergedStack *> &merged, 
	       set<int> &merge_list, 
	       int cohort_id, 
	       vector<MergedStack *> &deleveraged_tags) {
    set<int>::iterator it;
    vector<pair<int, int> >::iterator j;
    MergedStack *tag_1, *tag_2;
    uint k, l;

    //
    // Create a minimum spanning tree in order to determine the minimum distance
    // between each node in the list.
    //
    MinSpanTree *mst = new MinSpanTree;
    vector<int>  keys;

    for (it = merge_list.begin(); it != merge_list.end(); it++) {
	keys.push_back(*it);
        mst->add_node(*it);
        tag_1 = merged[*it];

	// cerr << "  " << *it << " -> " << tag_1->utags[0] << "\n";
    }

    //
    // Measure the distance between each pair of nodes and add edges to our
    // minimum spanning tree.
    //
    Node *n_1, *n_2;
    for (k = 0; k < keys.size(); k++) {
        tag_1 = merged[keys[k]];
        n_1   = mst->node(keys[k]);

	for (l = k+1; l < keys.size(); l++) {
	    tag_2 = merged[keys[l]];
            n_2   = mst->node(keys[l]);

	    int d = dist(tag_1, tag_2);

            n_1->add_edge(mst->node(keys[l]), d);
            n_2->add_edge(mst->node(keys[k]), d);
	}
    }

    //
    // Build the minimum spanning tree.
    //
    mst->build_tree();

    //
    // Visualize the MST
    //
    if (dump_graph) {
        stringstream gout_file;
        size_t pos_1 = in_file.find_last_of("/");
        size_t pos_2 = in_file.find_last_of(".");
        gout_file << out_path << in_file.substr(pos_1 + 1, (pos_2 - pos_1 - 1)) << "_" << keys[0] << ".dot";
        string vis = mst->vis(true);
        ofstream gvis(gout_file.str().c_str());
        gvis << vis;
        gvis.close();
    }

    set<int> visited;
    set<int, int_increasing> dists;
    queue<Node *> q;

    Node *n = mst->head();
    q.push(n);

    while (!q.empty()) {
        n = q.front();
        q.pop();
        visited.insert(n->id);

        for (uint i = 0; i < n->min_adj_list.size(); i++) {
            if (visited.count(n->min_adj_list[i]->id) == 0) {
                q.push(n->min_adj_list[i]);
		// cerr << n->id << " -> " << n->min_adj_list[i]->id << ": ";

		//
		// Find the edge distance.
		//
		for (uint j = 0; j < n->edges.size(); j++)
		    if (n->edges[j]->child == n->min_adj_list[i]) {
			// cerr << n->edges[j]->dist << "\n";
			dists.insert(n->edges[j]->dist);
		    }
	    }
        }
    }

    //
    // This set is sorted by definition. Check if there is more than a single 
    // distance separating stacks.
    //
    if (dists.size() == 1) {
	tag_1 = merge_tags(merged, merge_list, 0);
	deleveraged_tags.push_back(tag_1);
	delete mst;
	return 0;
    }

    uint min_dist = *(dists.begin());

    //
    // If there is more than a single distance, split the minimum spanning tree
    // into separate loci, by cutting the tree at the larger distances.
    //
    set<int> uniq_merge_list;
    visited.clear();
    n = mst->head();
    q.push(n);
    int id = 0;

    uniq_merge_list.insert(n->id);
    while (!q.empty()) {
        n = q.front();
        q.pop();
        visited.insert(n->id);

        for (uint i = 0; i < n->min_adj_list.size(); i++) {
            if (visited.count(n->min_adj_list[i]->id) == 0) {
                q.push(n->min_adj_list[i]);

		for (uint j = 0; j < n->edges.size(); j++) {
		    if (n->edges[j]->child == n->min_adj_list[i])
			if (n->edges[j]->dist > min_dist) {

			    // cerr << "Merging the following stacks into a locus:\n";
			    for (it = uniq_merge_list.begin(); it != uniq_merge_list.end(); it++) {
				tag_1 = merged[*it];
				// cerr << "  " << *it << " -> " << tag_1->utags[0] << "\n";
			    }

			    tag_1 = merge_tags(merged, uniq_merge_list, id);
			    tag_1->cohort_id = cohort_id;
			    deleveraged_tags.push_back(tag_1);
			    uniq_merge_list.clear();
			    id++;
			}
		}

		uniq_merge_list.insert(n->min_adj_list[i]->id);
	    }
        }
    }

    // cerr << "Merging the following stacks into a locus:\n";
    for (it = uniq_merge_list.begin(); it != uniq_merge_list.end(); it++) {
	tag_1 = merged[*it];
	// cerr << "  " << *it << " -> " << tag_1->utags[0] << "\n";
    }

    tag_1 = merge_tags(merged, uniq_merge_list, id);
    tag_1->cohort_id = cohort_id;
    deleveraged_tags.push_back(tag_1);
    uniq_merge_list.clear();

    delete mst;

    return 0;
}

int calc_kmer_distance(map<int, MergedStack *> &merged, int utag_dist) {
    //
    // Calculate the distance (number of mismatches) between each pair
    // of Radtags. We expect all radtags to be the same length;
    //
    KmerHashMap    kmer_map;
    vector<char *> kmer_map_keys;
    MergedStack   *tag_1, *tag_2;
    map<int, MergedStack *>::iterator it;

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

    populate_kmer_hash(merged, kmer_map, kmer_map_keys, kmer_len);
 
    #pragma omp parallel private(tag_1, tag_2)
    { 
        #pragma omp for schedule(dynamic) 
        for (uint i = 0; i < keys.size(); i++) {
            tag_1 = merged[keys[i]];

            // Don't compute distances for masked tags
            if (tag_1->masked) continue;

            vector<char *> query_kmers;
            generate_kmers(tag_1->con, kmer_len, num_kmers, query_kmers);

            map<int, int> hits;
            int d;
            //
            // Lookup the occurances of each k-mer in the kmer_map
            //
            for (int j = 0; j < num_kmers; j++) {
                for (uint k = 0; k < kmer_map[query_kmers[j]].size(); k++)
                    hits[kmer_map[query_kmers[j]][k]]++;
            }

            //
            // Free the k-mers we generated for this query
            //
            for (int j = 0; j < num_kmers; j++)
                delete [] query_kmers[j];

            // cerr << "  Tag " << tag_1->id << " hit " << hits.size() << " kmers.\n";

            //
            // Iterate through the list of hits. For each hit that has more than min_hits
            // check its full length to verify a match.
            //
            map<int, int>::iterator hit_it;
            for (hit_it = hits.begin(); hit_it != hits.end(); hit_it++) {
                // cerr << "  Tag " << hit_it->first << " has " << hit_it->second << " hits (min hits: " << min_hits << ")\n";

                if (hit_it->second < min_hits) continue;

                // cerr << "  Match found, checking full-length match\n";

                tag_2 = merged[hit_it->first];

                // Don't compute distances for masked tags
                if (tag_2->masked) continue;

                // Don't compare tag_1 against itself.
                if (tag_1 == tag_2) continue;

                d = dist(tag_1, tag_2);
                // cerr << "    Distance: " << d << "\n";

                //
                // Store the distance between these two sequences if it is
                // below the maximum distance (which governs those
                // sequences to be merged in the following step of the
                // algorithm.)
                //
                if (d <= utag_dist)
                    tag_1->add_dist(tag_2->id, d);
            }
            // Sort the vector of distances.
            sort(tag_1->dist.begin(), tag_1->dist.end(), compare_dist);
        }
    }

    free_kmer_hash(kmer_map, kmer_map_keys);

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

int reduce_radtags(DNASeqHashMap &radtags, map<int, Stack *> &unique, map<int, Rem *> &rem) {
    DNASeqHashMap::iterator it;

    Rem   *r;
    Stack *u;
    int   global_id = 1;

    for (it = radtags.begin(); it != radtags.end(); it++) {
    	if (it->second.count() < min_merge_cov) {
    	    //
    	    // Don't record this unique RAD-Tag if its coverage is below
    	    // the specified cutoff. However, add the reads to the remainder
    	    // vector for later processing.
    	    //
	    r     = new Rem;
    	    r->id = global_id;
    	    r->add_seq(it->first);

    	    for (uint i = 0; i < it->second.ids.size(); i++)
    		r->add_id(it->second.ids[i]);

	    rem[r->id] = r;
	    global_id++;

    	} else {
    	    //
    	    // Populate a Stack object for this unique radtag. Create a
    	    // map of the IDs for the sequences that have been
    	    // collapsed into this radtag.
    	    //
    	    u     = new Stack;
    	    u->id = global_id;
    	    u->add_seq(it->first);

    	    // Copy the original Fastq IDs from which this unique radtag was built.
    	    for (uint i = 0; i < it->second.ids.size(); i++)
    		u->add_id(it->second.ids[i]);

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

int
free_radtags_hash(DNASeqHashMap &radtags, vector<DNASeq *> &radtags_keys)
{
    for (uint i = 0; i < radtags_keys.size(); i++)
	delete radtags_keys[i];

    radtags.clear();

    return 0;
}

int calc_coverage_distribution(map<int, Stack *> &unique, double &mean, double &stdev) {
    map<int, Stack *>::iterator i;
    double m     = 0.0;
    double s     = 0.0;
    double sum   = 0.0;
    uint   max   = 0;
    uint   cnt   = 0;
    double total = 0.0;

    map<int, int> depth_dist;
    map<int, int>::iterator j;

    for (i = unique.begin(); i != unique.end(); i++) {
	cnt = i->second->count();
	m += cnt;
	total++;

        depth_dist[cnt]++;

	if (cnt > max)
	    max = cnt;
    }

    mean = round(m / total);

    //
    // Calculate the standard deviation
    //
    total = 0.0;

    for (i = unique.begin(); i != unique.end(); i++) {
	total++;
	s = i->second->count();
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
	    m   += tag->count();
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
	    s   += tag->count();
	}
	sum += pow((s - mean), 2);
    }

    stdev = sqrt(sum / (merged.size() - 1));

    cerr << "  Mean merged coverage depth is " << mean << "; Std Dev: " << stdev << "; Max: " << max << "\n";

    return m;
}

int count_raw_reads(map<int, Stack *> &unique, map<int, Rem *> &rem, map<int, MergedStack *> &merged) {
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
	    m   += tag->count();

            if (uniq_ids.count(*k) == 0)
               uniq_ids[*k] = 0;
            uniq_ids[*k]++;
	}
	for (uint j = 0; j < it->second->remtags.size(); j++)
	    m += rem[it->second->remtags[j]]->count();
	    //m += it->second->remtags.size();
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

int 
write_results(map<int, MergedStack *> &m, map<int, Stack *> &u, map<int, Rem *> &r) 
{
    map<int, MergedStack *>::iterator i;
    vector<int>::iterator      k;
    vector<SNP *>::iterator    s;
    map<string, int>::iterator t;
    MergedStack  *tag_1;
    Stack        *tag_2;
    Rem          *rem;
    stringstream  sstr;

    bool gzip = (in_file_type == gzfastq || in_file_type == gzfasta) ? true : false;

    //
    // Read in the set of sequencing IDs so they can be included in the output.
    //
    vector<char *> seq_ids;
    load_seq_ids(seq_ids);

    //
    // Parse the input file name to create the output files
    //
    size_t pos_1 = in_file.find_last_of("/");
    size_t pos_2 = in_file.find_last_of(".");

    if (in_file.substr(pos_2) == ".gz") {
	in_file = in_file.substr(0, pos_2);
	pos_2   = in_file.find_last_of(".");
    }

    string tag_file = out_path + in_file.substr(pos_1 + 1, (pos_2 - pos_1 - 1)) + ".tags.tsv";
    string snp_file = out_path + in_file.substr(pos_1 + 1, (pos_2 - pos_1 - 1)) + ".snps.tsv";
    string all_file = out_path + in_file.substr(pos_1 + 1, (pos_2 - pos_1 - 1)) + ".alleles.tsv";

    if (gzip) {
	tag_file += ".gz";
	snp_file += ".gz";
	all_file += ".gz";
    }

    //
    // Open the output files for writing.
    //
    gzFile   gz_tags, gz_snps, gz_alle;
    ofstream tags, snps, alle;
    if (gzip) {
	gz_tags = gzopen(tag_file.c_str(), "wb");
	if (!gz_tags) {
	    cerr << "Error: Unable to open gzipped tag file '" << tag_file << "': " << strerror(errno) << ".\n";
	    exit(1);
	}
	gzbuffer(gz_tags, libz_buffer_size);
	gz_snps = gzopen(snp_file.c_str(), "wb");
	if (!gz_snps) {
	    cerr << "Error: Unable to open gzipped snps file '" << snp_file << "': " << strerror(errno) << ".\n";
	    exit(1);
	}
	gzbuffer(gz_snps, libz_buffer_size);
	gz_alle = gzopen(all_file.c_str(), "wb");
	if (!gz_alle) {
	    cerr << "Error: Unable to open gzipped alleles file '" << all_file << "': " << strerror(errno) << ".\n";
	    exit(1);
	}
	gzbuffer(gz_alle, libz_buffer_size);
    } else {
	tags.open(tag_file.c_str());
	if (tags.fail()) {
	    cerr << "Error: Unable to open tag file for writing.\n";
	    exit(1);
	}
	snps.open(snp_file.c_str());
	if (snps.fail()) {
	    cerr << "Error: Unable to open SNPs file for writing.\n";
	    exit(1);
	}
	alle.open(all_file.c_str());
	if (alle.fail()) {
	    cerr << "Error: Unable to open allele file for writing.\n";
	    exit(1);
	}
    }

    int id;

    char *buf = new char[m.begin()->second->len + 1];

    for (i = m.begin(); i != m.end(); i++) {
	float total = 0;
	tag_1 = i->second;

	// First write the consensus sequence
	sstr << "0"              << "\t" 
	     << sql_id           << "\t" 
	     << tag_1->id        << "\t"
	    //<< tag_1->cohort_id << "\t"
	     << ""               << "\t" // chr
	     << 0                << "\t" // bp
	     << "+"              << "\t" // strand
	     << "consensus\t"    << "\t" 
	     << "\t" 
	     << tag_1->con         << "\t" 
	     << tag_1->deleveraged << "\t" 
	     << tag_1->blacklisted << "\t"
	     << tag_1->lumberjackstack << "\n";

	//
	// Write a sequence recording the output of the SNP model for each nucleotide.
	//
	sstr << "0" << "\t" 
	     << sql_id << "\t" 
	     << tag_1->id << "\t" 
	    //<< "\t" // cohort_id
	     << "\t"  // chr
	     << "\t"  // bp
	     << "\t"  // strand
	     << "model\t" << "\t"
	     << "\t";
	for (s = tag_1->snps.begin(); s != tag_1->snps.end(); s++) {
	    if ((*s)->type == snp_type_het)
		sstr << "E";
	    else if ((*s)->type == snp_type_hom)
		sstr << "O";
	    else
		sstr << "U";
	}
	sstr << "\t" 
	     << "\t"
	     << "\t"
	     << "\n";

	if (gzip) gzputs(gz_tags, sstr.str().c_str()); else tags << sstr.str();
	sstr.str("");

	// Now write out the components of each unique tag merged into this one.
	id = 0;
	for (k = tag_1->utags.begin(); k != tag_1->utags.end(); k++) {
	    tag_2  = u[*k];
	    total += tag_2->count();

	    for (uint j = 0; j < tag_2->map.size(); j++) {
		sstr << "0"       << "\t"
		     << sql_id    << "\t"
		     << tag_1->id << "\t"
		    //<< "\t" // cohort_id
		     << "\t" // chr
		     << "\t" // bp
		     << "\t" // strand
		     << "primary\t" 
		     << id << "\t" 
		     << seq_ids[tag_2->map[j]] << "\t" 
		     << tag_2->seq->seq(buf) 
		     << "\t\t\t\n";

		if (gzip) gzputs(gz_tags, sstr.str().c_str()); else tags << sstr.str();
		sstr.str("");
	    }

	    id++;
	}

	//
	// Write out the remainder tags merged into this unique tag.
	//
	for (k = tag_1->remtags.begin(); k != tag_1->remtags.end(); k++) {
	    rem    = r[*k];
	    total += rem->map.size();

	    for (uint j = 0; j < rem->map.size(); j++)
		sstr << "0"       << "\t" 
		     << sql_id    << "\t" 
		     << tag_1->id << "\t"
		    //<< "\t" // cohort_id
		     << "\t" // chr
		     << "\t" // bp
		     << "\t" // strand
		     << "secondary\t"
		     << "\t" 
		     << seq_ids[rem->map[j]] << "\t" 
		     << rem->seq->seq(buf) 
		     << "\t\t\t\n";

	    if (gzip) gzputs(gz_tags, sstr.str().c_str()); else tags << sstr.str();
	    sstr.str("");
	}

	//
	// Write out any SNPs detected in this unique tag.
	//
	for (s = tag_1->snps.begin(); s != tag_1->snps.end(); s++) {
	    if ((*s)->type == snp_type_het)
		sstr << "0" << "\t" << sql_id << "\t" << tag_1->id << "\t" 
		     << (*s)->col << "\t" << (*s)->lratio << "\t" 
		     << (*s)->rank_1 << "\t" << (*s)->rank_2 << "\t\t\n";
	}

	if (gzip) gzputs(gz_snps, sstr.str().c_str()); else snps << sstr.str();
	sstr.str("");

	//
	// Write the expressed alleles seen for the recorded SNPs and
	// the percentage of tags a particular allele occupies.
	//
	for (t = tag_1->alleles.begin(); t != tag_1->alleles.end(); t++) {
	    sstr << "0" << "\t" << sql_id << "\t" << tag_1->id << "\t" << (*t).first << "\t" << (((*t).second/total) * 100) << "\t" << (*t).second << "\n";
	}
	if (gzip) gzputs(gz_alle, sstr.str().c_str()); else alle << sstr.str();
	sstr.str("");

    }

    if (gzip) {
	gzclose(gz_tags);
	gzclose(gz_snps);
	gzclose(gz_alle);
    } else {
	tags.close();
	snps.close();
	alle.close();
    }

    //
    // Free sequence IDs.
    //
    for (uint i = 0; i < seq_ids.size(); i++)
	delete [] seq_ids[i];

    //
    // If specified, output reads not utilized in any stacks.
    //
    if (retain_rem_reads) {
	string unused_file = out_path + in_file.substr(pos_1 + 1, (pos_2 - pos_1 - 1)) + ".unused.fa";

	gzFile   gz_unused;
	ofstream unused;

	if (gzip) {
	    unused_file += ".gz";
	    gz_unused = gzopen(unused_file.c_str(), "wb");
	    if (!gz_unused) {
		cerr << "Error: Unable to open gzipped discard file '" << unused_file << "': " << strerror(errno) << ".\n";
		exit(1);
	    }
	    gzbuffer(gz_unused, libz_buffer_size);
	} else {
	    unused.open(unused_file.c_str());
	    if (unused.fail()) {
		cerr << "Error: Unable to open discard file for writing.\n";
		exit(1);
	    }
	}

 	map<int, Rem *>::iterator r_it;
	for (r_it = r.begin(); r_it != r.end(); r_it++) {
	    if (r_it->second->utilized == false)
		sstr << ">" << r_it->second->id << "\n" << r_it->second->seq->seq(buf) << "\n";
	    if (gzip) gzputs(gz_unused, sstr.str().c_str()); else unused << sstr.str();
	    sstr.str("");
	}

	if (gzip) gzclose(gz_unused); else unused.close();
    }

    delete [] buf;

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
	data << "/* " << keys[s] << ": " << unique[merged[keys[s]]->utags[0]]->map[0] << "; depth: " << merged[keys[s]]->count << " */\n";

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
    vector<pair<int, int> >::iterator pit;
    vector<int>::iterator mit;
    char *c;

    for (it = u.begin(); it != u.end(); it++) {
	c = (*it).second->seq->seq();

	cerr << "UniqueTag UID: " << (*it).second->id << "\n"
	     << "  Seq:       "   << c << "\n"
	     << "  IDs:       "; 

	for (uint j = 0; j < it->second->map.size(); j++)
	    cerr << it->second->map[j] << " ";

	cerr << "\n\n";

	delete [] c;
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

int load_radtags(string in_file, DNASeqHashMap &radtags, vector<DNASeq *> &radtags_keys) {
    Input *fh;
    DNASeq *d;

    if (in_file_type == fasta)
        fh = new Fasta(in_file.c_str());
    else if (in_file_type == fastq)
        fh = new Fastq(in_file.c_str());
    else if (in_file_type == gzfasta)
        fh = new GzFasta(in_file.c_str());
    else if (in_file_type == gzfastq)
        fh = new GzFastq(in_file.c_str());

    cerr << "Parsing " << in_file.c_str() << "\n";
    long  int corrected = 0;
    long  int i         = 0;
    short int seql      = 0;
    short int prev_seql = 0;
    bool  len_mismatch  = false;

    Seq c;
    c.id   = new char[id_len];
    c.seq  = new char[max_len];
    c.qual = new char[max_len];

    while ((fh->next_seq(c)) != 0) {
        if (i % 10000 == 0) cerr << "  Loading RAD-Tag " << i << "       \r";

	prev_seql = seql;
	seql      = 0;

	for (char *p = c.seq; *p != '\0'; p++, seql++)
	    switch (*p) {
	    case 'N':
	    case 'n':
	    case '.':
		*p = 'A';
		corrected++;
	    }

	if (seql != prev_seql && prev_seql > 0) len_mismatch = true;

	d = new DNASeq(seql, c.seq);

	pair<DNASeqHashMap::iterator, bool> r;

	r = radtags.insert(make_pair(d, HVal()));
	(*r.first).second.add_id(i);
	radtags_keys.push_back(d);
        i++;
    }
    cerr << "Loaded " << i << " RAD-Tags; inserted " << radtags.size() << " elements into the RAD-Tags hash map.\n";

    if (i == 0) {
        cerr << "Error: Unable to load data from '" << in_file.c_str() << "'.\n";
        exit(1);
    }

    cerr << "  " << corrected << " reads contained uncalled nucleotides that were modified.\n";

    if (len_mismatch)
	cerr << "  Warning: different sequence lengths detected, this will interfere with Stacks algorithms.\n";

    //
    // Close the file and delete the Input object.
    //
    delete fh;

    return 0;
}

int
load_seq_ids(vector<char *> &seq_ids)
{
    Input *fh;

    if (in_file_type == fasta)
        fh = new Fasta(in_file.c_str());
    else if (in_file_type == fastq)
        fh = new Fastq(in_file.c_str());
    else if (in_file_type == gzfasta)
        fh = new GzFasta(in_file.c_str());
    else if (in_file_type == gzfastq)
        fh = new GzFastq(in_file.c_str());

    cerr << "  Refetching sequencing IDs from " << in_file.c_str() << "... ";    

    char *id;
    Seq c;
    c.id   = new char[id_len];
    c.seq  = new char[max_len];
    c.qual = new char[max_len];

    while ((fh->next_seq(c)) != 0) {
	id = new char[strlen(c.id) + 1];
	strcpy(id, c.id);
	seq_ids.push_back(id);
    }
    cerr << "read " << seq_ids.size() << " sequence IDs.\n";

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
	    {"help",         no_argument,       NULL, 'h'},
            {"version",      no_argument,       NULL, 'v'},
	    {"infile_type",  required_argument, NULL, 't'},
	    {"file",         required_argument, NULL, 'f'},
	    {"outpath",      required_argument, NULL, 'o'},
	    {"id",           required_argument, NULL, 'i'},
	    {"min_cov",      required_argument, NULL, 'm'},
	    {"max_dist",     required_argument, NULL, 'M'},
	    {"max_sec_dist", required_argument, NULL, 'N'},
	    {"max_locus_stacks", required_argument, NULL, 'K'},
	    {"num_threads",  required_argument, NULL, 'p'},
	    {"deleverage",   no_argument,       NULL, 'd'},
	    {"remove_rep",   no_argument,       NULL, 'r'},
	    {"retain_rem",   no_argument,       NULL, 'R'},
	    {"graph",        no_argument,       NULL, 'g'},
	    {"exp_cov",      no_argument,       NULL, 'E'},
	    {"cov_stdev",    no_argument,       NULL, 's'},
	    {"cov_scale",    no_argument,       NULL, 'S'},
	    {"sec_hapl",     no_argument,       NULL, 'H'},
	    {"model_type",   required_argument, NULL, 'T'},
	    {"bc_err_freq",  required_argument, NULL, 'e'},
	    {"bound_low",    required_argument, NULL, 'L'},
	    {"bound_high",   required_argument, NULL, 'U'},
	    {"alpha",        required_argument, NULL, 'A'},
	    {0, 0, 0, 0}
	};
	
	// getopt_long stores the option index here.
	int option_index = 0;
     
	c = getopt_long(argc, argv, "hHvdrgA:L:U:f:o:i:m:e:E:s:S:p:t:M:N:K:T:", long_options, &option_index);
     
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
            else if (strcmp(optarg, "fasta") == 0)
                in_file_type = fasta;
            else if (strcmp(optarg, "fastq") == 0)
                in_file_type = fastq;
	    else if (strcasecmp(optarg, "gzfasta") == 0)
                in_file_type = gzfasta;
	    else if (strcasecmp(optarg, "gzfastq") == 0)
                in_file_type = gzfastq;
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
	    sql_id = is_integer(optarg);
	    if (sql_id < 0) {
		cerr << "SQL ID (-i) must be an integer, e.g. 1, 2, 3\n";
		help();
	    }
	    break;
	case 'm':
	    min_merge_cov = atoi(optarg);
	    break;
	case 'M':
	    max_utag_dist = atoi(optarg);
	    break;
	case 'N':
	    max_rem_dist = atoi(optarg);
	    break;
	case 'd':
	    deleverage_stacks++;
	    break;
	case 'r':
	    remove_rep_stacks++;
	    break;
	case 'K':
	    max_subgraph = atoi(optarg);
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
            } else if (strcmp(optarg, "bounded") == 0) {
                model_type = bounded;
            } else {
                cerr << "Unknown model type specified '" << optarg << "'\n";
                help();
            }
	case 'e':
	    barcode_err_freq = atof(optarg);
	    break;
	case 'L':
	    bound_low  = atof(optarg);
	    break;
	case 'U':
	    bound_high = atof(optarg);
	    break;
	case 'A':
	    alpha = atof(optarg);
	    break;
	case 'H':
	    call_sec_hapl = false;
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

    if (alpha != 0.1 && alpha != 0.05 && alpha != 0.01 && alpha != 0.001) {
	cerr << "SNP model alpha significance level must be either 0.1, 0.05, 0.01, or 0.001.\n";
	help();
    }

    if (bound_low != 0 && (bound_low < 0 || bound_low >= 1.0)) {
	cerr << "SNP model lower bound must be between 0.0 and 1.0.\n";
	help();
    }

    if (bound_high != 1 && (bound_high <= 0 || bound_high > 1.0)) {
	cerr << "SNP model upper bound must be between 0.0 and 1.0.\n";
	help();
    }

    if (bound_low > 0 || bound_high < 1.0) {
	model_type = bounded;
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
              << "ustacks -t file_type -f file_path [-d] [-r] [-o path] [-i id] [-m min_cov] [-M max_dist] [-p num_threads] [-R] [-H] [-h]" << "\n"
	      << "  t: input file Type. Supported types: fasta, fastq, gzfasta, or gzfastq.\n"
              << "  f: input file path.\n"
	      << "  o: output path to write results." << "\n"
	      << "  i: SQL ID to insert into the output to identify this sample." << "\n"
	      << "  m: Minimum depth of coverage required to create a stack (default 3)." << "\n"
	      << "  M: Maximum distance (in nucleotides) allowed between stacks (default 2)." << "\n"
	      << "  N: Maximum distance allowed to align secondary reads to primary stacks (default: M + 2).\n"
	      << "  R: retain unused reads.\n"
	      << "  H: disable calling haplotypes from secondary reads.\n"
              << "  p: enable parallel execution with num_threads threads.\n"
	      << "  h: display this help messsage.\n\n"
	      << "  Stack assembly options:\n"
	      << "    r: enable the Removal algorithm, to drop highly-repetitive stacks (and nearby errors) from the algorithm." << "\n"
	      << "    d: enable the Deleveraging algorithm, used for resolving over merged tags." << "\n"
	      << "    --max_locus_stacks <num>: maximum number of stacks at a single de novo locus (default 3).\n"
	      << "  Model options:\n" 
	      << "    --model_type: either 'snp' (default), 'bounded', or 'fixed'\n"
	      << "    For the SNP or Bounded SNP model:\n"
	      << "      --alpha <num>: chi square significance level required to call a heterozygote or homozygote, either 0.1, 0.05 (default), 0.01, or 0.001.\n"
	      << "    For the Bounded SNP model:\n"
	      << "      --bound_low <num>: lower bound for epsilon, the error rate, between 0 and 1.0 (default 0).\n"
	      << "      --bound_high <num>: upper bound for epsilon, the error rate, between 0 and 1.0 (default 1).\n"
	      << "    For the Fixed model:\n"
	      << "      --bc_err_freq <num>: specify the barcode error frequency, between 0 and 1.0.\n";

    exit(0);
}
