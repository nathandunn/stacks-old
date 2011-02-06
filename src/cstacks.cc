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
// cstacks -- Create a catalog of Stacks.
//
// Julian Catchen
// jcatchen@uoregon.edu
// University of Oregon
//
// $Id$
//

#include "cstacks.h"

// Global variables to hold command-line options.
queue<pair<int, string> > samples;
string  out_path;
int     batch_id     = 0;
int     mult_matches = 0;
int     ctag_dist    = 0;
searcht search_type = sequence;
int     num_threads  = 1;

int main (int argc, char* argv[]) {

    parse_command_line(argc, argv);

    //
    // Set the number of OpenMP parallel threads to execute.
    //
    #ifdef _OPENMP
    omp_set_num_threads(num_threads);
    #endif

    map<int, CLocus *> catalog;
    map<int, CLocus *>::iterator cat_it;

    pair<int, string> s = samples.front();
    samples.pop();

    cerr << "Initializing catalog...\n";
    if (!initialize_catalog(s, catalog)) {
        cerr << "Failed to initialize the catalog.\n";
        return 1;
    }

    int i = 2;
    while (!samples.empty()) {
        map<int, QLocus *> sample;

	cerr << "Processing sample " << i << "\n";

	s = samples.front();
	samples.pop();

	if (!load_loci(s.second, sample)) {
            cerr << "Failed to load sample " << i << "\n";
            continue;
        }
	//dump_loci(sample);

        if (search_type == sequence) {
            cerr << "Searching for sequence matches...\n";
            find_kmer_matches_by_sequence(catalog, sample, ctag_dist);
        } else if (search_type == genomic_loc) {
            cerr << "Searching for matches by genomic location...\n";
            find_matches_by_genomic_loc(catalog, sample);
        }

	cerr << "Merging matches into catalog...\n";
	merge_matches(catalog, sample, s, ctag_dist);

        //
        // Regenerate the alleles for the catalog tags after merging the new sample into the catalog.
	//
        for (cat_it = catalog.begin(); cat_it != catalog.end(); cat_it++)
            cat_it->second->populate_alleles();

	i++;
    }

    cerr << "Writing catalog...\n";
    write_catalog(catalog);

    return 0;
}

int characterize_mismatch_snps(CLocus *catalog_tag, QLocus *query_tag) {
    vector<Match *>::iterator mat_it;
    vector<pair<allele_type, string> >::iterator all_it;
    
    mat_it = query_tag->matches.begin();

    const char *c = NULL;
    //
    // Fetch the proper catalog allele
    //
    for (all_it = catalog_tag->strings.begin(); all_it != catalog_tag->strings.end(); all_it++)
        if (all_it->first == (*mat_it)->cat_type) 
            c = all_it->second.c_str();
    if (c == NULL) return 0;

    const char *q = NULL;
    //
    // Fetch the proper query allele
    //
    for (all_it = query_tag->strings.begin(); all_it != query_tag->strings.end(); all_it++)
        if (all_it->first == (*mat_it)->query_type) 
            q = all_it->second.c_str();
    if (q == NULL) return 0;

    //
    // For each mismatch found, create a SNP object
    //
    const char *beg = c;
    const char *end = c + strlen(c);

    while (c < end) {
	if (*c != *q) {
            SNP *s = new SNP;
            s->col    =  c - beg;
            s->lratio = 0;
            s->rank_1 = *c;
            s->rank_2 = *q;

            merge_allele(catalog_tag, s);
            merge_allele(query_tag, s);

            catalog_tag->snps.push_back(s);

            s = new SNP;
            s->col    =  c - beg;
            s->lratio = 0;
            s->rank_1 = *q;
            s->rank_2 = *c;

            query_tag->snps.push_back(s);
        }
	c++;
	q++;
    }

    return 1;
}

int merge_matches(map<int, CLocus *> &catalog, map<int, QLocus *> &sample, pair<int, string> &sample_file, int ctag_dist) {
    map<int, QLocus *>::iterator i;
    vector<Match *>::iterator mat_it;
    CLocus *ctag;
    QLocus *qtag;

    for (i = sample.begin(); i != sample.end(); i++) {
	qtag = i->second;

        //
        // If this stack didn't match an existing catalog stack, add this stack to the 
        // catalog as a new stack.
        //
	if (qtag->matches.size() == 0) {
	    add_unique_tag(sample_file, catalog, qtag);
	    continue;
	}

        //
        // Check for multiple matches. We will reduce the list of Match objects, which
        // contain matches to multiple alleles for a single locus, to the smallest distance
        // for a locus.
        //
	map<int, uint> local_matches;
	map<int, uint>::iterator j;
	for (mat_it = qtag->matches.begin(); mat_it != qtag->matches.end(); mat_it++) {
            if (local_matches.count((*mat_it)->cat_id) == 0)
                local_matches[(*mat_it)->cat_id] = (*mat_it)->dist;
            else if ((*mat_it)->dist < local_matches[(*mat_it)->cat_id])
                local_matches[(*mat_it)->cat_id] = (*mat_it)->dist;
        }

        uint min_dist    =  1000;
        uint num_matches =  0;
        int  min_cat_id  = -1;
        //
        // Find the minimum distance and then check how many matches have that distance.
        //
        for (j = local_matches.begin(); j != local_matches.end(); j++)
            min_dist = j->second < min_dist ? j->second : min_dist;

        for (j = local_matches.begin(); j != local_matches.end(); j++)
            if (j->second == min_dist) {
                num_matches++;
                min_cat_id = j->first;
            }

	//
	// Emit a warning if the query tag matches more than one tag in the catalog.
	//
	if (num_matches > 1) {
	    cerr << 
		"  Warning: sample " << sample_file.second << ", tag " << qtag->id << 
		", matches more than one tag in the catalog: ";
	    for (j = local_matches.begin(); j != local_matches.end(); j++)
		cerr << j->first << " ";
	    cerr << "\n";

	    //
	    // Don't record matches to multiple catalog entries unless instructed
	    // to do so by the command line option.
	    //
	    if (!mult_matches) continue;
	}

        ctag = catalog[min_cat_id];

        if (ctag == NULL) 
            cerr << "  Unable to locate catalog tag " << min_cat_id << "\n";

        //
        // If mismatches are allowed between query and catalog tags, identify the 
        // mismatches and convert them into SNP objects to be merged into the catalog tag.
        //
        if (ctag_dist > 0 && !characterize_mismatch_snps(ctag, qtag))
            cerr 
                << "  Error characterizing mismatch SNPs " 
                << sample_file.second << ", tag " << qtag->id 
                << " with catalog tag " << ctag->id << "\n";

	//
	// Merge the SNPs and alleles from the sample into the catalog tag.
	//
	if (!ctag->merge_snps(qtag)) {
	    cerr << "Error merging " << sample_file.second << ", tag " << qtag->id <<
		" with catalog tag " << ctag->id << "\n";
	}

	ctag->sources.push_back(make_pair(sample_file.first, qtag->id));
    }

    return 0;
}

int add_unique_tag(pair<int, string> &sample_file, map<int, CLocus *> &catalog, QLocus *qloc) {
    vector<SNP *>::iterator i;
    map<string, int>::iterator j;

    int cid = catalog.size();

    CLocus *c = new CLocus;
    c->id = cid + 1;
    c->add_consensus(qloc->con);
    c->sources.push_back(make_pair(sample_file.first, qloc->id));
    catalog[c->id] = c;

    //
    // Add the physical genome location of this locus.
    //
    strncpy(c->loc.chr, qloc->loc.chr, id_len);
    c->loc.chr[id_len - 1] = '\0';
    c->loc.bp = qloc->loc.bp;

    // cerr << "Adding sample: " << qloc->id << " to the catalog as ID: " << c->id << "\n";

    for (i = qloc->snps.begin(); i != qloc->snps.end(); i++) {
	SNP *snp    = new SNP;
	snp->col    = (*i)->col;
	snp->lratio = (*i)->lratio;
	snp->rank_1 = (*i)->rank_1;
	snp->rank_2 = (*i)->rank_2;

	c->snps.push_back(snp);
    }

    for (j = qloc->alleles.begin(); j != qloc->alleles.end(); j++) {
	c->alleles[j->first] = j->second;
    }

    c->populate_alleles();

    return 0;
}

int find_kmer_matches_by_sequence(map<int, CLocus *> &catalog, map<int, QLocus *> &sample, int ctag_dist) {
    //
    // Calculate the distance (number of mismatches) between each pair
    // of Radtags. We expect all radtags to be the same length;
    //
    CatKmerHashMap kmer_map;
    map<int, QLocus *>::iterator it;
    vector<pair<allele_type, string> >::iterator allele;
    QLocus *tag_1;
    CLocus *tag_2;
    int i, j;

    // OpenMP can't parallelize random access iterators, so we convert
    // our map to a vector of integer keys.
    vector<int> keys;
    for (it = sample.begin(); it != sample.end(); it++) 
	keys.push_back(it->first);

    //
    // Calculate the number of k-mers we will generate. If kmer_len == 0,
    // determine the optimal length for k-mers.
    //
    int con_len   = strlen(sample[keys[0]]->con);
    int kmer_len  = determine_kmer_length(con_len, ctag_dist);
    int num_kmers = con_len - kmer_len + 1;

    cerr << "  Number of kmers per sequence: " << num_kmers << "\n";

    //
    // Calculate the minimum number of matching k-mers required for a possible sequence match.
    //
    int min_hits = calc_min_kmer_matches(kmer_len, ctag_dist, con_len, true);

    populate_kmer_hash(catalog, kmer_map, kmer_len);

    cerr << "  " << catalog.size() << " items in the catalog, " << kmer_map.size() << " elements in the catalog hash.\n";
 
    #pragma omp parallel private(i, j, tag_1, tag_2, allele)
    { 
        #pragma omp for schedule(dynamic) 
        for (i = 0; i < (int) keys.size(); i++) {
            tag_1 = sample[keys[i]];

            for (allele = tag_1->strings.begin(); allele != tag_1->strings.end(); allele++) {

                vector<char *> kmers;
                generate_kmers(allele->second.c_str(), kmer_len, num_kmers, kmers);

                map<int, vector<allele_type> > hits;
                vector<pair<allele_type, int> >::iterator map_it;
                int d;
                //
                // Lookup the occurances of each k-mer in the kmer_map
                //
                for (j = 0; j < num_kmers; j++) {

                    if (kmer_map.count(kmers[j]) > 0)
                        for (map_it  = kmer_map[kmers[j]].begin();
                             map_it != kmer_map[kmers[j]].end();
                             map_it++)
                            hits[map_it->second].push_back(map_it->first);
                }

                //
                // Free the allocated k-mers.
                //
                for (j = 0; j < num_kmers; j++)
                    delete [] kmers[j];
                kmers.clear();

                //cerr << "  Tag " << tag_1->id << " hit " << hits.size() << " kmers.\n";

                map<int, vector<allele_type> >::iterator hit_it;
                vector<allele_type>::iterator            all_it;

                //
                // Iterate through the list of hits. For each hit, total up the hits to the various alleles.
                // Any allele that has more than min_hits check its full length to verify a match.
                //
                for (hit_it = hits.begin(); hit_it != hits.end(); hit_it++) {
                    //cerr << "  Tag " << hit_it->first << " has " << hit_it->second << " hits (min hits: " << min_hits << ")\n";

                    map<allele_type, int>           allele_cnts;
                    map<allele_type, int>::iterator cnt_it;

                    for (all_it = hit_it->second.begin(); all_it != hit_it->second.end(); all_it++)
                        allele_cnts[*all_it]++;

                    for (cnt_it = allele_cnts.begin(); cnt_it != allele_cnts.end(); cnt_it++) {
                        //cerr << "      allele " << cnt_it->first << " has " << cnt_it->second << " hits\n"; 

                        if (cnt_it->second < min_hits) continue;

                        //cerr << "  Match found, checking full-length match\n";

                        tag_2 = catalog[hit_it->first];

                        d = dist(allele->second.c_str(), tag_2, cnt_it->first);

                        if (d < 0)
                            cerr << 
                                "Unknown error calculating distance between " << 
                                tag_1->id << " and " << tag_2->id << "; query allele: " << allele->first << "\n";

                        //cerr << "    Distance: " << d << " CTAG_DIST: " << ctag_dist << "\n";

                        //
                        // Add a match to the query sequence: catalog ID, catalog allele, query allele, distance
                        //
                        if (d <= ctag_dist)
                            tag_1->add_match(tag_2->id, cnt_it->first, allele->first, d);
                    }
                }
            }

            // Sort the vector of distances.
            sort(tag_1->matches.begin(), tag_1->matches.end(), compare_matches);
        }
    }

    return 0;
}

bool compare_matches(Match *a, Match *b) {
    return (a->dist < b->dist);
}

int find_matches_by_sequence(map<int, CLocus *> &catalog, map<int, QLocus *> &sample) {
    //
    // Calculate the distance (number of mismatches) between each pair
    // of Radtags. We expect all radtags to be the same length;
    //
    map<int, QLocus *>::iterator i;
    map<int, CLocus *>::iterator j;
    int k;

    // OpenMP can't parallelize random access iterators, so we convert
    // our map to a vector of integer keys.
    vector<int> keys;
    for (i = sample.begin(); i != sample.end(); i++) 
	keys.push_back(i->first);

    #pragma omp parallel private(i, j, k)
    {
        #pragma omp for schedule(dynamic) 
	for (k = 0; k < (int) keys.size(); k++) {

	    i = sample.find(keys[k]);

	    vector<pair<allele_type, string> >::iterator r, s;

	    //
	    // Iterate through the possible SAMPLE alleles
	    //
	    for (r = i->second->strings.begin(); r != i->second->strings.end(); r++) {

		for (j = catalog.begin(); j != catalog.end(); j++) {
		    //
		    // Iterate through the possible CATALOG alleles
		    //
		    for (s = j->second->strings.begin(); s != j->second->strings.end(); s++) {
			if (r->second == s->second) {
			    //cerr << "Found a match between " << i->first << " (" << r->first << ") and " << j->first << " (" << s->first << ")\n";

			    i->second->add_match(j->second->id, s->first, r->first, 0);
			}
		    }
		}
	    }
	}
    }

    return 0;
}

int find_matches_by_genomic_loc(map<int, CLocus *> &catalog, map<int, QLocus *> &sample) {
    //
    // Calculate the distance (number of mismatches) between each pair
    // of Radtags. We expect all radtags to be the same length;
    //
    map<int, QLocus *>::iterator i;
    map<int, CLocus *>::iterator j;
    int k;

    // OpenMP can't parallelize random access iterators, so we convert
    // our map to a vector of integer keys.
    vector<int> keys;
    for (i = sample.begin(); i != sample.end(); i++) 
	keys.push_back(i->first);

    #pragma omp parallel private(i, j, k)
    {
        #pragma omp for schedule(dynamic) 
	for (k = 0; k < (int) keys.size(); k++) {

	    i = sample.find(keys[k]);

            for (j = catalog.begin(); j != catalog.end(); j++) {
                
                if (strcmp(i->second->loc.chr, j->second->loc.chr) == 0 && 
                    i->second->loc.bp == j->second->loc.bp) {

                    i->second->add_match(j->second->id, "", "", 0);
                }
            }
        }
    }

    return 0;
}

int write_catalog(map<int, CLocus *> &catalog) {
    map<int, CLocus *>::iterator i;
    CLocus  *tag;
    set<int> matches;

    //
    // Parse the input file names to create the output file
    //
    stringstream prefix; 
    prefix << out_path << "batch_" << batch_id;

    //
    // Output the tags
    //
    string out_file = prefix.str() + ".catalog.tags.tsv";
    ofstream cat_file(out_file.c_str());
    out_file = prefix.str() + ".catalog.snps.tsv";
    ofstream snp_file(out_file.c_str());
    out_file = prefix.str() + ".catalog.alleles.tsv";
    ofstream all_file(out_file.c_str());

    for (i = catalog.begin(); i != catalog.end(); i++) {
	tag = i->second;

	write_simple_output(tag, cat_file, snp_file, all_file);

// 	//
// 	// Debugging code.
// 	//
//  	if (tag->matches.size() > 0) {
//  	    cerr << "Found match: " << tag->matches.size() << " matches.\n";
//  	}

//  	for (mat_it = tag->matches.begin(); mat_it != tag->matches.end(); mat_it++) {
//  	    matched_tag = parent_2[(*mat_it).first];
//  	    if (strcmp(matched_tag->con, tag->con) == 0)
//  		cerr << "  Match " << tag->id << " -> " << matched_tag->id << ": consensus matches.\n";
//  	    else 
//  		cerr << "  Match " << tag->id << " -> " << matched_tag->id << ": consensus does not match;\n" 
//  		     << "    " << matched_tag->con << "\n" 
//  		     << "    " << tag->con << "\n";
//  	}

    }

    cat_file.close();
    snp_file.close();
    all_file.close();

    return 0;
}

int merge_allele(Locus *locus, SNP *snp) {
    map<int, pair<string, SNP *> > columns;
    map<int, pair<string, SNP *> >::iterator c;
    vector<SNP *>::iterator i;

    for (i = locus->snps.begin(); i != locus->snps.end(); i++)
	columns[(*i)->col] = make_pair("sample", *i);

    //
    // Is this column already represented from the previous sample?
    //
    if (columns.count(snp->col)) {
        //
        // Check to make sure we aren't combining inconsistent SNPS.
        //
        if (columns[snp->col].second->rank_1 != snp->rank_1 &&
            columns[snp->col].second->rank_1 != snp->rank_2)
            cerr << "  Warning: inconsistent merging of SNPs for catalog ID " << locus->id << "\n";
        else
            columns[snp->col] = make_pair("both", snp);
    }
    else
        columns[snp->col] = make_pair("merge", snp);

    vector<pair<string, SNP *> > merged_snps;

    for (c = columns.begin(); c != columns.end(); c++) 
	merged_snps.push_back((*c).second);

    //
    // Sort the SNPs by column
    //
    sort(merged_snps.begin(), merged_snps.end(), compare_pair);

    //
    // Modify any existing alleles to account for this new SNP. If there are not any alleles, 
    // create new ones.
    //
    stringstream sallele;
    set<string> merged_alleles;
    string allele, new_allele;
    int pos;

    if (locus->alleles.size() == 0) {
        sallele << locus->con[snp->col];
        merged_alleles.insert(sallele.str());
    }

    map<string, int>::iterator j;
    vector<pair<string, SNP *> >::iterator k;    

    for (j = locus->alleles.begin(); j != locus->alleles.end(); j++) {
	allele     = j->first;
	new_allele = "";
	pos        = 0;

        //cerr << "Allele length: " << allele.size() << "\n";

	for (k = merged_snps.begin(); k != merged_snps.end(); k++) {
	    //
	    // If we inserted a SNP from the sample, add the proper nucleotide from the consensus
	    // sequence to account for it in the allele string.
	    //
	    if ((*k).first == "merge") {
		new_allele += locus->con[(*k).second->col];
                //cerr << "  Adding char from consensus position " << (*k).second->col << "\n";
	    } else {
		new_allele += allele[pos];
                //cerr << "  Adding char from allele position " << pos << "\n";
		pos++;
	    }
	}

	merged_alleles.insert(new_allele);
    }

    set<string>::iterator s;

    locus->alleles.clear();
    for (s = merged_alleles.begin(); s != merged_alleles.end(); s++) {
	locus->alleles[*s] = 0;
    }

    return 1;
}

bool compare_pair(pair<string, SNP *> a, pair<string, SNP *> b) {
    return (a.second->col < b.second->col);
}

int CLocus::merge_snps(QLocus *matched_tag) {
    vector<SNP *>::iterator i;
    map<string, int>::iterator j;
    vector<pair<string, SNP *> >::iterator k;
    map<int, pair<string, SNP *> > columns;
    map<int, pair<string, SNP *> >::iterator c;

    vector<pair<string, SNP *> > merged_snps;
    set<string> merged_alleles;
    set<string>::iterator s;

    for (i = this->snps.begin(); i != this->snps.end(); i++)
	columns[(*i)->col] = make_pair("catalog", *i);

    for (i = matched_tag->snps.begin(); i != matched_tag->snps.end(); i++) {
	//
	// Is this column already represented from the previous sample?
	//
	if (columns.count((*i)->col)) {
	    //
	    // Check to make sure we aren't combining inconsistent SNPS.
	    //
	    if (columns[(*i)->col].second->rank_1 != (*i)->rank_1 &&
		columns[(*i)->col].second->rank_1 != (*i)->rank_2)
		cerr << "  Warning: inconsistent merging of SNPs for catalog ID " << this->id << "\n";
	    else
		columns[(*i)->col] = make_pair("both", *i);
	}
	else
	    columns[(*i)->col] = make_pair("sample", *i);
    }

    for (c = columns.begin(); c != columns.end(); c++) 
	merged_snps.push_back((*c).second);

    //
    // Sort the SNPs by column
    //
    sort(merged_snps.begin(), merged_snps.end(), compare_pair);

    //
    // Merge the alleles accounting for any SNPs added from either of the two samples.
    //
    string allele, new_allele;
    int pos;

    for (j = this->alleles.begin(); j != this->alleles.end(); j++) {
	allele     = j->first;
	new_allele = "";
	pos        = 0;

	for (k = merged_snps.begin(); k != merged_snps.end(); k++) {
	    //
	    // If we inserted a SNP from the sample, add the proper nucleotide from the consensus
	    // sequence to account for it in the allele string.
	    //
	    if ((*k).first == "sample") {
		new_allele += this->con[(*k).second->col];
	    } else {
		new_allele += allele[pos];
		pos++;
	    }
	}

	merged_alleles.insert(new_allele);
    }

    for (j = matched_tag->alleles.begin(); j != matched_tag->alleles.end(); j++) {
	allele     = j->first;
	new_allele = "";
	pos        = 0;

	for (k = merged_snps.begin(); k != merged_snps.end(); k++) {
	    if ((*k).first == "catalog") {
		new_allele += matched_tag->con[(*k).second->col];
	    } else {
		new_allele += allele[pos];
		pos++;
	    }
	}

	merged_alleles.insert(new_allele);
    }

    //
    // Update the catalog entry's list of SNPs and alleles
    //
    this->snps.clear();

    for (k = merged_snps.begin(); k != merged_snps.end(); k++) {
	SNP *snp    = new SNP;
	snp->col    = (*k).second->col;
	snp->lratio = (*k).second->lratio;
	snp->rank_1 = (*k).second->rank_1;
	snp->rank_2 = (*k).second->rank_2;

	this->snps.push_back(snp);
    }

    for (k = merged_snps.begin(); k != merged_snps.end(); k++)
        delete k->second;

    this->alleles.clear();
    for (s = merged_alleles.begin(); s != merged_alleles.end(); s++) {
	this->alleles[*s] = 0;
    }

    return 1;
}

int populate_kmer_hash(map<int, CLocus *> &catalog, CatKmerHashMap &kmer_map, int kmer_len) {
    map<int, CLocus *>::iterator it;
    vector<pair<allele_type, string> >::iterator allele;
    vector<char *> kmers;
    CLocus        *tag;
    char          *hash_key;
    bool           exists;
    int            j;

    //
    // Break each stack down into k-mers and create a hash map of those k-mers
    // recording in which sequences they occur.
    //
    int num_kmers = strlen(catalog.begin()->second->con) - kmer_len + 1;

    for (it = catalog.begin(); it != catalog.end(); it++) {
        tag = it->second;

        //
        // Iterate through the possible Catalog alleles
        //
        for (allele = tag->strings.begin(); allele != tag->strings.end(); allele++) {
            //
            // Generate and hash the kmers for this allele string
            //
            generate_kmers(allele->second.c_str(), kmer_len, num_kmers, kmers);

            for (j = 0; j < num_kmers; j++) {

                exists = kmer_map.count(kmers[j]) == 0 ? false : true;

                if (exists) {
                    hash_key = kmers[j];
                } else {
                    hash_key = new char [strlen(kmers[j]) + 1];
                    strcpy(hash_key, kmers[j]);
                }

                kmer_map[hash_key].push_back(make_pair(allele->first, tag->id));
            }

            for (j = 0; j < num_kmers; j++)
                delete [] kmers[j];
            kmers.clear();
        }
    }

    //dump_kmer_map(kmer_map);

    return 0;
}

int write_simple_output(CLocus *tag, ofstream &cat_file, ofstream &snp_file, ofstream &all_file) {
    vector<SNP *>::iterator           snp_it;
    map<string, int>::iterator        all_it;
    vector<pair<int, int> >::iterator src_it;
    string sources;

    for (src_it = tag->sources.begin(); src_it != tag->sources.end(); src_it++) {
	stringstream s; 
	s << (*src_it).first << "_" << (*src_it).second << ",";
	sources += s.str();
    }
    sources = sources.substr(0, sources.length() - 1);

    cat_file << 
	"0"          << "\t" << 
	batch_id     << "\t" <<
	tag->id      << "\t" <<
        tag->loc.chr << "\t" <<
        tag->loc.bp  << "\t" <<
	"consensus"  << "\t" <<
	"0"          << "\t" <<
	sources      << "\t" <<
	tag->con     << "\t" << 
        0            << "\t" <<  // These flags are unused in cstacks, but important in ustacks
        0            << "\t" <<
        0            << "\n";

    //
    // Output the SNPs associated with the catalog tag
    //
    for (snp_it = tag->snps.begin(); snp_it != tag->snps.end(); snp_it++)
	snp_file << "0\t" << 
	    batch_id          << "\t" <<
	    tag->id          << "\t" << 
	    (*snp_it)->col    << "\t" << 
	    (*snp_it)->lratio << "\t" << 
	    (*snp_it)->rank_1 << "\t" << 
	    (*snp_it)->rank_2 << "\n";

    //
    // Output the alleles associated with the two matched tags
    //
    for (all_it = tag->alleles.begin(); all_it != tag->alleles.end(); all_it++)
	all_file << "0\t" << 
	    batch_id  << "\t" <<
	    tag->id  << "\t" << 
	    all_it->first << "\t" <<
            0 << "\t" <<              // These two fields are used in the pstacks output, not in
            0 << "\n";

    return 0;
}

int initialize_catalog(pair<int, string> &sample, map<int, CLocus *> &catalog) {
    map<int, CLocus *> tmp_catalog;

    //
    // Parse the input files.
    //
    if (!load_loci(sample.second, tmp_catalog))
        return 0;

    //
    // Iterate over the catalog entires and renumber them after recording the source of
    // locus.
    //
    map<int, CLocus *>::iterator j;
    int k = 1;
    for (j = tmp_catalog.begin(); j != tmp_catalog.end(); j++) {
	j->second->sources.push_back(make_pair(sample.first, j->second->id));
        j->second->id = k;

        catalog[k] = j->second;

	k++;        
    }

    return 1;
}

int parse_command_line(int argc, char* argv[]) {
    int c, sid;
    string sstr;

    while (1) {
	static struct option long_options[] = {
	    {"help",         no_argument,       NULL, 'h'},
            {"version",      no_argument,       NULL, 'v'},
	    {"mmatches",     no_argument,       NULL, 'm'},
	    {"genomic_loc",  no_argument,       NULL, 'g'},
	    {"batch_id",     required_argument, NULL, 'b'},
	    {"ctag_dist",    required_argument, NULL, 'n'},
	    {"sample",       required_argument, NULL, 's'},
	    {"sample_id",    required_argument, NULL, 'S'},
	    {"outpath",      required_argument, NULL, 'o'},
	    {"num_threads",  required_argument, NULL, 'p'},
	    {0, 0, 0, 0}
	};
	
	// getopt_long stores the option index here.
	int option_index = 0;
     
	c = getopt_long(argc, argv, "hgvmo:s:S:b:p:n:", long_options, &option_index);
     
	// Detect the end of the options.
	if (c == -1)
	    break;
     
	switch (c) {
	case 'h':
	    help();
	    break;
	case 'b':
	    batch_id = atoi(optarg);
	    break;
	case 'n':
	    ctag_dist = atoi(optarg);
	    break;
	case 'm':
	    mult_matches++;
	    break;
	case 'g':
	    search_type = genomic_loc;
	    break;
	case 's':
	    sstr = optarg;

	    c = getopt_long(argc, argv, "ho:s:S:b:", long_options, &option_index);
	    if (c != 'S') {
		cerr << "You must specify a sample ID after each sample file.\n";
		help();
		abort();
	    }

	    sid = atoi(optarg);
	    samples.push(make_pair(sid, sstr));
	    break;
     	case 'o':
	    out_path = optarg;
	    break;
        case 'v':
            version();
            break;
	case 'p':
	    num_threads = atoi(optarg);
	    break;
	case '?':
	    // getopt_long already printed an error message.
	    help();
	    break;
	default:
	    help();
	    abort();
	}
    }

    if (samples.size() == 0) {
	cerr << "You must specify at least one sample file.\n";
	help();
    }

    if (out_path.length() == 0) 
	out_path = ".";

    if (out_path.at(out_path.length() - 1) != '/') 
	out_path += "/";

    return 0;
}

void version() {
    std::cerr << "cstacks " << VERSION << "\n\n";

    exit(0);
}

void help() {
    std::cerr << "cstacks " << VERSION << "\n"
              << "cstacks -b batch_id -s sample_file -S id [-s sample_file_2 -S id_2 ...] [-o path] [-p num_threads] [-n num] [-g] [-h]" << "\n"
              << "  p: enable parallel execution with num_threads threads.\n"
	      << "  b: MySQL ID of this batch." << "\n"
	      << "  s: TSV file from which to load radtags." << "\n"
	      << "  S: MySQL ID of the sample." << "\n"
	      << "  o: output path to write results." << "\n"
	      << "  m: include tags in the catalog that match to more than one entry." << "\n"
              << "  n: number of mismatches allowed between sample tags when generating the catalog." << "\n"
              << "  g: base catalog matching on genomic location, not sequence identity." << "\n"
	      << "  h: display this help messsage." << "\n\n";

    exit(0);
}
