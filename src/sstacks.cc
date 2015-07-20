// -*-mode:c++; c-style:k&r; c-basic-offset:4;-*-
//
// Copyright 2010-2015, Julian Catchen <jcatchen@illinois.edu>
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
// sstacks -- search for occurances of stacks in a catalog of stacks.
//

#include "sstacks.h"

// Global variables to hold command-line options.
queue<string> samples;
string  catalog_path;
string  out_path;
FileT   in_file_type = FileT::sql;
int     num_threads  = 1;
int     batch_id     = 0;
int     samp_id      = 0; 
int     catalog      = 0;
bool    verify_haplotypes       = true;
bool    impute_haplotypes       = true;
bool    require_uniq_haplotypes = false;
searcht search_type             = sequence;

int main (int argc, char* argv[]) {

    parse_command_line(argc, argv);

    uint sample_cnt = samples.size();

    //
    // Set the number of OpenMP parallel threads to execute.
    //
    #ifdef _OPENMP
    omp_set_num_threads(num_threads);
    #endif

    map<int, Locus *>  catalog;
    bool compressed = false;
    int  res;

    if (search_type == sequence) 
 	cerr << "Searching for matches by sequence identity...\n";
    else if (search_type == genomic_loc)
	cerr << "Searching for matches by genomic location...\n";

    catalog_path += ".catalog";
    res = load_loci(catalog_path, catalog, false, false, compressed);

    if (res == 0) {
	cerr << "Unable to parse catalog, '" << catalog_path << "'\n";
	return 0;
    }

    string sample_path;
    int    i = 1;

    while (!samples.empty()) {
        map<int, QLocus *> sample;

	sample_path = samples.front();
	samples.pop();

	cerr << "Processing sample '" << sample_path << "' [" << i << " of " << sample_cnt << "]\n";

	res = load_loci(sample_path, sample, false, false, compressed);

	if (res == 0) {
	    cerr << "Unable to parse '" << sample_path << "'\n";
	    return 0;
	}

	in_file_type = compressed == true ? FileT::gzsql : FileT::sql;

	//
	// Assign the ID for this sample data.
	//
	samp_id = sample.begin()->second->sample_id;

	//dump_loci(catalog);
	//dump_loci(sample);

	if (search_type == sequence) {
	    cerr << "Searching for sequence matches...\n";
	    find_matches_by_sequence(catalog, sample);

	} else if (search_type == genomic_loc) {
	    cerr << "Searching for matches by genomic location...\n";
	    find_matches_by_genomic_loc(catalog, sample);
	}

	write_matches(sample_path, sample);
	i++;
    }

    return 0;
}

int find_matches_by_genomic_loc(map<int, Locus *> &sample_1, map<int, QLocus *> &sample_2) {
    //
    // Calculate the distance (number of mismatches) between each pair
    // of Radtags. We expect all radtags to be the same length;
    //
    map<int, QLocus *>::iterator i;
    map<int, Locus *>::iterator j;
    int  k;
    char id[id_len];

    //
    // Build a hash map out of the first sample (usually the catalog)
    //
    //
    // Create a map of the genomic locations of stacks in sample_1
    //
    cerr << "  Creating map of genomic locations...";

    map<string, set<int> > locations;
    for (j = sample_1.begin(); j != sample_1.end(); j++) {
        snprintf(id, id_len - 1, "%s|%d|%c", 
		 j->second->loc.chr, 
		 j->second->loc.bp, 
		 j->second->loc.strand == plus ? '+' : '-');
        locations[id].insert(j->second->id);
    }

    cerr << "done.\n";

    // OpenMP can't parallelize random access iterators, so we convert
    // our map to a vector of integer keys.
    vector<int> keys;
    for (i = sample_2.begin(); i != sample_2.end(); i++) 
	keys.push_back(i->first);

    //
    // Initialize some counters
    //
    unsigned long matches = 0;
    unsigned long nomatch = 0;
    unsigned long nosnps  = 0;
    unsigned long tot_hap = 0;
    unsigned long ver_hap = 0;

    #pragma omp parallel private(i, j, k, id)
    {
	unsigned long verified;
        #pragma omp for reduction(+:matches) reduction(+:tot_hap) reduction(+:ver_hap) reduction(+:nomatch) reduction(+:nosnps)
	for (k = 0; k < (int) keys.size(); k++) {

	    i = sample_2.find(keys[k]);
	    snprintf(id, id_len - 1, "%s|%d|%c", 
		     i->second->loc.chr, 
		     i->second->loc.bp, 
		     i->second->loc.strand == plus ? '+' : '-');

	    if (locations.count(id) > 0) {
		Locus *tag;
	        set<int>::iterator loc_it;
		vector<pair<allele_type, string> >::iterator q;

		matches++;

		for (loc_it = locations[id].begin(); loc_it != locations[id].end(); loc_it++) {
		    tag = sample_1[*loc_it];

		    //
		    // Generate haplotypes for query tag relative to the catalog tag.
		    //
		    set<string> query_haplotypes;
		    generate_query_haplotypes(tag, i->second, query_haplotypes);
		    tot_hap += query_haplotypes.size() > 0 ? query_haplotypes.size() : 1;

		    if (verify_haplotypes) {
			verified = verify_genomic_loc_match(tag, i->second, query_haplotypes, nosnps);
			ver_hap += verified;
			if (verified == 0) nomatch++;
		    } else {
			i->second->add_match(tag->id, tag->strings.begin()->first);
		    }
		}
	    }
	}
    }

    cerr << keys.size() << " stacks matched against the catalog containing " << sample_1.size() << " loci.\n" 
	 << "  " << matches << " matching loci, " << nomatch << " contained no verified haplotypes.\n"
	 << "  " << nosnps  << " loci contained SNPs unaccounted for in the catalog and were excluded.\n"
	 << "  " << tot_hap << " total haplotypes examined from matching loci, " << ver_hap << " verified.\n";

    return 0;
}

int verify_genomic_loc_match(Locus *s1_tag, QLocus *s2_tag, set<string> &query_haplotypes, unsigned long &nosnps) {
    vector<SNP *>::iterator i, j;

    //
    // We have found a match between the genomic location of s1 and s2. We now want
    // to verify that the haplotypes are consistent between the tags, i.e. they 
    // have the same number and types of SNPs.
    //

    //
    // 1. First, if there are no SNPs present in either the query or catalog, just
    //    check that the strings match.
    //
    uint min_len = s1_tag->len > s2_tag->len ? s2_tag->len : s1_tag->len;

    if (s1_tag->snps.size() == 0 && 
	s2_tag->snps.size() == 0 && 
	strncmp(s1_tag->con, s2_tag->con, min_len) == 0) {
	s2_tag->add_match(s1_tag->id, "consensus");
	return 1;
    }

    //
    // 2. Second, we will check that the query locus (s2_tag) does not have any SNPs 
    //    lacking in the catalog tag (s1_tag).
    //
    bool found;
    for (j = s2_tag->snps.begin(); j != s2_tag->snps.end(); j++) {
	found = false;
	//
	// SNP occurs in a column that is beyond the length of the catalog
	//
	if ((*j)->col > min_len - 1)
	    continue;

	for (i = s1_tag->snps.begin(); i != s1_tag->snps.end(); i++) {
	    if ((*i)->col == (*j)->col)
		found = true;
	}
	//
	// Query locus posses a SNP not present in the catalog.
	//
	if (found == false) {
	    nosnps++;
	    return 0;
	}
    }

    //
    // Finally, check that one of the constructed alleles matches the allele
    // passed in on the stack.
    //
    string cat_haplotype;
    vector<pair<allele_type, string> >::iterator c;
    set<string>::iterator a;

    uint matches = 0;
    for (a = query_haplotypes.begin(); a != query_haplotypes.end(); a++) {

	if (impute_haplotypes) {
	    int res = impute_haplotype(*a, s1_tag->strings, cat_haplotype);

	    if (res > 0) {
		//
		// If the matching haplotype was imputed, record the depths of the query alleles
		// under the new, imputed alleles.
		//
		if (s2_tag->alleles.count(cat_haplotype) == 0) {
		    if (s2_tag->alleles.count(*a) > 0)
			s2_tag->alleles[cat_haplotype] = s2_tag->alleles[*a];
		    else
			s2_tag->alleles[cat_haplotype] = s2_tag->depth;
		}
		//cerr << s2_tag->id << "; Adding cat haplotype: " << cat_haplotype << " based on depth of " << *a << ", " << s2_tag->alleles[cat_haplotype] << "\n";
		s2_tag->add_match(s1_tag->id, cat_haplotype);
		matches++;
	    } else if (res < 0) {
		cerr << "  Failure imputing haplotype for catalog locus: " << s1_tag->id << " and query tag: " << s2_tag->id << "\n";
	    }
	} else {
	    for (c = s1_tag->strings.begin(); c != s1_tag->strings.end(); c++)
		if (*a == c->first) {
		    //cerr << "  Adding match between " << s1_tag->id << " and " << c->first << "\n";
		    s2_tag->add_match(s1_tag->id, c->first);
		    matches++;
		}
	}
    }

    return matches;
}

// int impute_haplotype(string query_haplotype, 
// 		     vector<pair<allele_type, string> > &cat_haplotypes, 
// 		     string &match) {

//     uint max_len = query_haplotype.length() > cat_haplotypes[0].first.length() ?
// 	query_haplotype.length() : 
// 	cat_haplotypes[0].first.length();

//     //cerr << "Query len: " << query_haplotype.length() << "; Max length: " << max_len << "\n";

//     vector<string> cur, next;

//     for (uint i = 0; i < cat_haplotypes.size(); i++)
// 	cur.push_back(cat_haplotypes[i].first);
//     match = "";

//     //
//     // Examine the haplotypes one SNP at a time. If we are able to uniquely 
//     // determine the catalog haplotype that the query haplotype corresponds 
//     // to, return it.
//     //
//     uint j = 0;
//     while (cur.size() > 1 && j < max_len) {

// 	for (uint i = 0; i < cur.size(); i++) {
// 	    //cerr << "Comparing query[" << j << "]: '" << query_haplotype[j] << "' to catalog '" << cur[i][j] << "'\n"; 
// 	    if (query_haplotype[j] == cur[i][j]) {
// 		//cerr << "  Keeping this haplotype.\n";
// 		next.push_back(cur[i]);
// 	    }
// 	}
// 	cur = next;
// 	next.clear();
// 	j++;
//     }

//     //
//     // If there is only one left, make sure what we have of the haplotype does match
//     // and its not simply an erroneously called haplotype.
//     //
//     if (cur.size() == 1 && 
// 	strncmp(cur[0].c_str(), query_haplotype.c_str(), max_len) == 0) {
// 	match = cur[0];
// 	return 1;
//     }

//     //
//     // If, after examining all the available SNPs in the query haplotype, there is
//     // still more than a single possible catalog haplotype, then we can't impute it.
//     //
//     return 0;
// }

int impute_haplotype(string query_haplotype, 
		     vector<pair<allele_type, string> > &cat_haplotypes, 
		     string &match) {

    if (cat_haplotypes.size() == 0) {
	cerr << "Warning: malformed catalog tag: missing haplotype information.\n";
	return -1;
    }

    //cerr << "Examining " << query_haplotype << "\n";
    uint max_len = query_haplotype.length() > cat_haplotypes[0].first.length() ?
	query_haplotype.length() : 
	cat_haplotypes[0].first.length();

    //cerr << "Query len: " << query_haplotype.length() << "; Max length: " << max_len << "\n";

    vector<string> cur, next;
    uint match_cnt, no_n_cnt;

    for (uint i = 0; i < cat_haplotypes.size(); i++)
	cur.push_back(cat_haplotypes[i].first);
    match = "";

    //
    // Examine the haplotypes one SNP at a time. If we are able to uniquely 
    // determine the catalog haplotype that the query haplotype corresponds 
    // to, return it.
    //
    uint j = 0;
    while (cur.size() > 1 && j < max_len) {

	for (uint i = 0; i < cur.size(); i++) {
	    //cerr << "Comparing query[" << j << "]: '" << query_haplotype << "' to catalog '" << cur[i] << "'\n"; 
	    if (require_uniq_haplotypes && (query_haplotype[j] == cur[i][j] || query_haplotype[j] == 'N')) {
		//cerr << "  Keeping this haplotype.\n";
		next.push_back(cur[i]);
	    } else if (query_haplotype[j] == cur[i][j]) {
		//cerr << "  Keeping this haplotype.\n";
		next.push_back(cur[i]);
	    } //else {
		//cerr << "  Discarding this haplotype.\n";
	    //}
	}
	cur = next;
	next.clear();
	j++;
    }

    //
    // If there is only one left, make sure what we have of the haplotype does match
    // and its not simply an erroneously called haplotype.
    //
    no_n_cnt  = 0;
    match_cnt = 0;
    if (cur.size() == 1) {
	if (require_uniq_haplotypes) {
	    for (uint k = 0; k < max_len; k++)
		if (query_haplotype[k] != 'N') no_n_cnt++;
	    for (uint k = 0; k < max_len; k++)
		if (cur[0][k] == query_haplotype[k]) match_cnt++;

	    if (match_cnt == no_n_cnt) {
		//cerr << "Keeping " << query_haplotype << "\n";
		match = cur[0];
		return 1;
	    }
	} else {
	    if (strncmp(cur[0].c_str(), query_haplotype.c_str(), max_len) == 0) {
		match = cur[0];
		return 1;
	    }
	}
    }

    //
    // If, after examining all the available SNPs in the query haplotype, there is
    // still more than a single possible catalog haplotype, then we can't impute it.
    //
    return 0;
}

int 
generate_query_haplotypes(Locus *s1_tag, QLocus *s2_tag, set<string> &query_haplotypes)
{
    //
    // Construct a set of haplotypes from the query locus relative to the catalog locus.
    // (The query locus already has a set of haplotypes, however, they don't necessarily 
    //  account for all the SNPs in the catalog, so we will augment them with sequence
    //  from the consensus.)
    //
    if (s1_tag->snps.size() == 0 && s2_tag->snps.size() == 0)
	return 0;

    vector<pair<string, SNP *> >   merged_snps;
    map<int, pair<string, SNP *> > columns;
    map<int, pair<string, SNP *> >::iterator c;
    vector<pair<string, SNP *> >::iterator   k;
    vector<SNP *>::iterator i;

    for (i = s1_tag->snps.begin(); i != s1_tag->snps.end(); i++)
	columns[(*i)->col] = make_pair("catalog", *i);

    for (i = s2_tag->snps.begin(); i != s2_tag->snps.end(); i++) {
	//
	// Is this column already represented in the catalog?
	//
	if (columns.count((*i)->col))
	    columns[(*i)->col] = make_pair("both", *i);
	else
	    columns[(*i)->col] = make_pair("query", *i);
    }

    for (c = columns.begin(); c != columns.end(); c++) 
	merged_snps.push_back((*c).second);

    //
    // Sort the SNPs by column
    //
    sort(merged_snps.begin(), merged_snps.end(), compare_pair_snp);

    map<string, int> converted_alleles;
    map<string, int>::iterator b;
    string old_allele, new_allele;
    int    pos;

    for (b = s2_tag->alleles.begin(); b != s2_tag->alleles.end(); b++) {
	old_allele = b->first;
	new_allele = "";
	pos        = 0;

	for (k = merged_snps.begin(); k != merged_snps.end(); k++) {
	    //
	    // If the SNPs from the catalog haplotype beyond the length of the query, add Ns
	    //
	    if (k->first == "catalog") {
		new_allele += (k->second->col > s2_tag->len - 1) ? 'N' : s2_tag->con[k->second->col];
	    } else {
		new_allele += old_allele[pos];
		pos++;
	    }
	}
	query_haplotypes.insert(new_allele);
	converted_alleles[new_allele] = b->second;

	// cerr << "Adding haplotype: " << new_allele << " [" << b->first << "]\n";
    }

    if (s2_tag->alleles.size() == 0) {
	new_allele = "";
	for (k = merged_snps.begin(); k != merged_snps.end(); k++) {
	    new_allele += (k->second->col > s2_tag->len - 1) ? 'N' : s2_tag->con[k->second->col];
	}
	query_haplotypes.insert(new_allele);
	// cerr << "Adding haplotype 2: " << new_allele << "\n";
    } else {
	s2_tag->alleles.clear();
	for (b = converted_alleles.begin(); b != converted_alleles.end(); b++)
	    s2_tag->alleles[b->first] = b->second;
    }

    return 0;
}

int 
find_matches_by_sequence(map<int, Locus *> &sample_1, map<int, QLocus *> &sample_2)
{
    map<int, QLocus *>::iterator i;
    uint min_tag_len;

    //
    // We don't assume all radtags will be the same length.
    //
    min_tag_len = 
	sample_1.begin()->second->len > sample_2.begin()->second->len ?
	sample_2.begin()->second->len : sample_1.begin()->second->len; 

    //
    // Build a hash map out of the first sample (usually the catalog),
    // using only the minimum length substring of the longest reads;
    //
    HashMap sample_1_map;
    populate_hash(sample_1, sample_1_map, min_tag_len);

    // OpenMP can't parallelize random access iterators, so we convert
    // our map to a vector of integer keys.
    vector<int> keys;
    for (i = sample_2.begin(); i != sample_2.end(); i++) 
	keys.push_back(i->first);

    //
    // Initialize some counters
    //
    unsigned long matches = 0;
    unsigned long mmatch  = 0;
    unsigned long nosnps  = 0;
    unsigned long nomatch = 0;
    unsigned long tot_hap = 0;
    unsigned long ver_hap = 0;

    #pragma omp parallel
    {
        #pragma omp for reduction(+:matches) reduction(+:tot_hap) reduction(+:ver_hap) reduction(+:nomatch) reduction(+:mmatch)
 	for (uint k = 0; k < keys.size(); k++) {
	    QLocus *query = sample_2[keys[k]];	    

            //
            // Iterate through the haplotypes for this tag in sample_2
            //
            HashMap::iterator hit;
            vector<pair<allele_type, string> >::iterator q;  // Query records allele_type/search string pairs
            vector<pair<int, allele_type> >::iterator c;     // Hash map records id/allele_type pairs
	    map<string, vector<string> > haplo_hits;
	    set<int> loci_hit;

            for (q = query->strings.begin(); q != query->strings.end(); q++) {
                // cerr << "  Looking for haplotype: " << q->first << " with sequence " << q->second.substr(0, min_tag_len) << "\n";

                hit = sample_1_map.find(q->second.substr(0, min_tag_len).c_str());

                if (hit != sample_1_map.end()) {
		    tot_hap++;
                    // cerr << "    Found a match for " << hit->first << "\n";

                    for (c = hit->second.begin(); c != hit->second.end(); c++) {
			//
			// Record the catalog loci hit by the haplotypes of this query locus.
			//
			loci_hit.insert(c->first);

			//
			// Record the haplotypes hit between the query and catalog loci.
			//
			haplo_hits[q->first].push_back(c->second);

			if (verify_haplotypes == false)
			    query->add_match(c->first, c->second);
		    }
                }
            }

	    if (loci_hit.size() > 0) matches++;

	    if (verify_haplotypes && loci_hit.size() > 0) {
		uint verified = verify_sequence_match(sample_1, query, loci_hit, haplo_hits, 
						      min_tag_len, mmatch, nosnps);
		ver_hap += verified;
		if (verified == 0) nomatch++;
	    }
        }
    }

    cerr << keys.size() << " stacks compared against the catalog containing " << sample_1.size() << " loci.\n" 
	 << "  " << matches << " matching loci, " << nomatch << " contained no verified haplotypes.\n"
	 << "  " << mmatch  << " loci matched more than one catalog locus and were excluded.\n"
	 << "  " << nosnps  << " loci contained SNPs unaccounted for in the catalog and were excluded.\n"
	 << "  " << tot_hap << " total haplotypes examined from matching loci, " << ver_hap << " verified.\n";

    return 0;
}

int verify_sequence_match(map<int, Locus *> &sample_1, QLocus *query, 
			  set<int> &loci_hit, map<string, vector<string> > &haplo_hits, 
			  uint min_tag_len, unsigned long &mmatch, unsigned long &nosnps) {
    //
    // 1. Check that this query locus matches just a single catalog locus.
    //
    if (loci_hit.size() > 1) {
        mmatch++;
        return 0;
    }

    Locus *cat = sample_1[*(loci_hit.begin())];

    //
    // 2. Make sure the query has no SNPs unaccounted for in the catalog.
    //
    vector<SNP *>::iterator i, j;
    bool found;

    for (i = query->snps.begin(); i != query->snps.end(); i++) {
	found = false;
	//
	// SNP occurs in a column that is beyond the length of the catalog
	//
	if ((*i)->col > min_tag_len - 1)
	    continue;

	for (j = cat->snps.begin(); j != cat->snps.end(); j++) {
	    if ((*i)->col == (*j)->col)
		found = true;
	}
	//
	// Query locus posses a SNP not present in the catalog.
	//
	if (found == false) {
	    nosnps++;
	    return 0;
	}
    }

    //
    // 3. We want a one-to-one correspondance between a query haplotype and a 
    //    catalog haplotype. This relationship fails when the catalog and query seqeunces
    //    are different lengths and the full length haplotype can not be determined.
    //
    map<string, vector<string> >::iterator it;
    map<string, int> cat_hap, query_hap;
    
    for (it = haplo_hits.begin(); it != haplo_hits.end(); it++) {
	query_hap[it->first] = it->second.size();
	for (uint j = 0; j < it->second.size(); j++)
	    cat_hap[it->second[j]]++;
    }

    uint verified = 0;
    for (it = haplo_hits.begin(); it != haplo_hits.end(); it++)
	for (uint j = 0; j < it->second.size(); j++) {
	    if (cat_hap[it->second[j]] == 1 &&
		query_hap[it->first] == 1) {
		verified++;
		query->add_match(cat->id, it->second[j]);
		//
		// If the matching haplotype was imputed, record the depths of the query alleles
		// under the new, imputed alleles.
		//
		if (query->alleles.count(it->second[j]) == 0) {
		    if (query->alleles.count(it->first) > 0)
			query->alleles[it->second[j]] = query->alleles[it->first];
		    else
			query->alleles[it->second[j]] = query->depth;
		}
	    }
	}

    return verified;
}

int populate_hash(map<int, Locus *> &sample, HashMap &hash_map, int min_tag_len) {
    map<int, Locus *>::iterator it;
    vector<pair<allele_type, string> >::iterator all_it;
    Locus *tag;
    char  *key;

    //
    // Create a hash map out of the set of alleles for each Locus.
    //
    for (it = sample.begin(); it != sample.end(); it++) {
        tag = it->second;

        for (all_it = tag->strings.begin(); all_it != tag->strings.end(); all_it++) {
	    key = new char[min_tag_len + 1];
	    strncpy(key, all_it->second.c_str(), min_tag_len);
	    key[min_tag_len] = '\0';

	    hash_map[key].push_back(make_pair(tag->id, all_it->first));
	}
    }

    //dump_kmer_map(kmer_map);

    return 0;
}

int 
write_matches(string sample_path, map<int, QLocus *> &sample) 
{
    map<int, QLocus *>::iterator i;

    //
    // Parse the input file names to create the output file
    //
    size_t pos_1    = sample_path.find_last_of("/");
    string out_file = out_path + sample_path.substr(pos_1 + 1)  + ".matches.tsv";

    if (in_file_type == FileT::gzsql)
	out_file += ".gz";

    //
    // Open the output files for writing.
    //
    gzFile   gz_matches;
    ofstream matches;
    if (in_file_type == FileT::gzsql) {
	gz_matches = gzopen(out_file.c_str(), "wb");
	if (!gz_matches) {
	    cerr << "Error: Unable to open gzipped matches file '" << out_file << "': " << strerror(errno) << ".\n";
	    exit(1);
	}
        #if ZLIB_VERNUM >= 0x1240
	gzbuffer(gz_matches, libz_buffer_size);
	#endif
    } else {
	matches.open(out_file.c_str());
	if (matches.fail()) {
	    cerr << "Error: Unable to open matches file for writing.\n";
	    exit(1);
	}
    }

    //
    // Record the version of Stacks used and the date generated as a comment in the catalog.
    //
    // Obtain the current date.
    //
    stringstream log;
    time_t       rawtime;
    struct tm   *timeinfo;
    char         date[32];
    time(&rawtime);
    timeinfo = localtime(&rawtime);
    strftime(date, 32, "%F %T", timeinfo);
    log << "# sstacks version " << VERSION << "; generated on " << date << "\n"; 
    if (in_file_type == FileT::gzsql) 
        gzputs(gz_matches, log.str().c_str());
    else
        matches << log.str();

    QLocus *qloc;
    string       type;
    uint         match_depth;
    stringstream sstr;

    cerr << "Outputing to file " << out_file << "\n";

    for (i = sample.begin(); i != sample.end(); i++) {
	qloc = i->second;

	for (uint j = 0; j < qloc->matches.size(); j++) {
	    if (verify_haplotypes == false && search_type == genomic_loc)
		match_depth = qloc->depth;
	    else
		match_depth = 
		    qloc->alleles.count(qloc->matches[j]->cat_type) > 0 ? 
		    qloc->alleles[qloc->matches[j]->cat_type] : qloc->depth;

	    sstr << 
		"0"            << "\t" <<
		batch_id       << "\t" <<
		qloc->matches[j]->cat_id   << "\t" <<
		samp_id        << "\t" <<
		qloc->id  << "\t" << 
		qloc->matches[j]->cat_type << "\t" <<
		match_depth    << "\t" <<
		qloc->lnl << "\n";
	}

	if (in_file_type == FileT::gzsql) gzputs(gz_matches, sstr.str().c_str()); else matches << sstr.str();
	sstr.str("");
    }

        if (in_file_type == FileT::gzsql)
	    gzclose(gz_matches);
	else
	    matches.close();

    return 0;
}

int parse_command_line(int argc, char* argv[]) {
	string sample_file;
	int c;
     
    while (1) {
	static struct option long_options[] = {
	    {"help",              no_argument, NULL, 'h'},
            {"version",           no_argument, NULL, 'v'},
	    {"genomic_loc",       no_argument, NULL, 'g'},
	    {"verify_hap",        no_argument, NULL, 'x'},
	    {"uniq_haplotypes",   no_argument, NULL, 'u'},
	    {"num_threads", required_argument, NULL, 'p'},
	    {"batch_id",    required_argument, NULL, 'b'},
	    {"catalog",     required_argument, NULL, 'c'},
	    {"sample_2",    required_argument, NULL, 's'},
	    {"outpath",     required_argument, NULL, 'o'},
	    {0, 0, 0, 0}
	};
	
	// getopt_long stores the option index here.
	int option_index = 0;
     
	c = getopt_long(argc, argv, "hgxuvs:c:o:b:p:", long_options, &option_index);
     
	// Detect the end of the options.
	if (c == -1)
	    break;
     
	switch (c) {
	case 'h':
	    help();
	    break;
	case 'p':
	    num_threads = atoi(optarg);
	    break;
	case 'b':
	    batch_id = is_integer(optarg);
	    if (batch_id < 0) {
		cerr << "Batch ID (-b) must be an integer, e.g. 1, 2, 3\n";
		help();
	    }
	    break;
	case 's':
	    sample_file = optarg;
	    samples.push(sample_file);
	    break;
	case 'g':
	    search_type = genomic_loc;
	    break;
     	case 'o':
	    out_path = optarg;
	    break;
	case 'c': 
	    catalog_path = optarg;
	    break;
	case 'x': 
	    verify_haplotypes = false;
	    break;
	case 'u':
	    require_uniq_haplotypes = true;
	    break;
        case 'v':
            version();
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

    if (catalog_path.length() == 0) {
	cerr << "You must specify the prefix path to the catalog.\n";
	help();
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
    std::cerr << "sstacks " << VERSION << "\n\n";

    exit(0);
}

void help() {
    std::cerr << "sstacks " << VERSION << "\n"
              << "sstacks -b batch_id -c catalog_file -s sample_file [-s sample_file_2 ...] [-o path] [-p num_threads] [-g] [-x] [-v] [-h]" << "\n"
              << "  p: enable parallel execution with num_threads threads.\n"
	      << "  b: MySQL ID of this batch." << "\n"
	      << "  c: TSV file from which to load the catalog loci." << "\n"
	      << "  s: filename prefix from which to load sample loci." << "\n"
	      << "  o: output path to write results." << "\n"
              << "  g: base matching on genomic location, not sequence identity." << "\n"
	      << "  x: don't verify haplotype of matching locus." << "\n"
	      << "  v: print program version." << "\n"
	      << "  h: display this help messsage." << "\n\n";

    exit(0);
}
