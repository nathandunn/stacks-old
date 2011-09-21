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
// sstacks -- search for occurances of stacks in a catalog of stacks.
//
// Julian Catchen
// jcatchen@uoregon.edu
// University of Oregon
//
// $Id$
//

#include "sstacks.h"

// Global variables to hold command-line options.
string  sample_1_file;
string  sample_2_file;
string  out_path;
int     num_threads = 1;
int     batch_id    = 0;
int     samp_id     = 0; 
int     catalog     = 0;
bool    verify_haplotypes = true;
bool    impute_haplotypes = true;
searcht search_type = sequence;

int main (int argc, char* argv[]) {

    parse_command_line(argc, argv);

    //
    // Set the number of OpenMP parallel threads to execute.
    //
    #ifdef _OPENMP
    omp_set_num_threads(num_threads);
    #endif

    map<int, Locus *>  sample_1;
    map<int, QLocus *> sample_2;
    int res;

    if (catalog) {
        sample_1_file += ".catalog";
	res = load_loci(sample_1_file, sample_1, false);
    } else {
	res = load_loci(sample_1_file, sample_1, false);
    }

    if (res == 0) {
	cerr << "Unable to parse '" << sample_1_file << "'\n";
	return 0;
    }

    res = load_loci(sample_2_file, sample_2, false);

    if (res == 0) {
	cerr << "Unable to parse '" << sample_2_file << "'\n";
	return 0;
    }

    //dump_loci(sample_1);
    //dump_loci(sample_2);

    if (search_type == sequence) {
	cerr << "Searching for sequence matches...\n";
	find_matches_by_sequence(sample_1, sample_2);

    } else if (search_type == genomic_loc) {
	cerr << "Searching for matches by genomic location...\n";
	find_matches_by_genomic_loc(sample_1, sample_2);
    }

    write_matches(sample_2);

    return 0;
}

int find_matches_by_genomic_loc(map<int, Locus *> &sample_1, map<int, QLocus *> &sample_2) {
    //
    // Calculate the distance (number of mismatches) between each pair
    // of Radtags. We expect all radtags to be the same length;
    //
    map<int, QLocus *>::iterator i;
    map<int, Locus *>::iterator j;
    int k, min_tag_len;
    char id[id_len];

    min_tag_len = strlen(sample_1.begin()->second->con); 

    //
    // Build a hash map out of the first sample (usually the catalog)
    //
    //
    // Create a map of the genomic locations of stacks in sample_1
    //
    cerr << "  Creating map of genomic locations...";

    map<string, set<int> > locations;
    for (j = sample_1.begin(); j != sample_1.end(); j++) {
        snprintf(id, id_len - 1, "%s_%d", j->second->loc.chr, j->second->loc.bp);
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
    unsigned long tot_hap = 0;
    unsigned long ver_hap = 0;

    #pragma omp parallel private(i, j, k, id)
    {
	unsigned long verified;
        #pragma omp for reduction(+:matches) reduction(+:tot_hap) reduction(+:ver_hap) reduction(+:nomatch)
	for (k = 0; k < (int) keys.size(); k++) {

	    i = sample_2.find(keys[k]);
	    snprintf(id, id_len - 1, "%s_%d", i->second->loc.chr, i->second->loc.bp);

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
			verified = verify_genomic_loc_match(tag, i->second, query_haplotypes);
			ver_hap += verified;
			if (verified == 0) nomatch++;
		    } else {
			for (q = tag->strings.begin(); q != tag->strings.end(); q++) {
			    i->second->add_match(tag->id, q->first);
			    ver_hap++;
			}
		    }
		}
	    }
	}
    }

    cerr << keys.size() << " stacks matched against the catalog containing " << sample_1.size() << " loci.\n" 
	 << "  " << matches << " matching loci, " << nomatch << " contained no verified haplotypes.\n"
	 << "  " << tot_hap << " total haplotypes examined from matching loci, " << ver_hap << " verified.\n";

    return 0;
}

int verify_genomic_loc_match(Locus *s1_tag, QLocus *s2_tag, set<string> &query_haplotypes) {
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

	for (i = s1_tag->snps.begin(); i != s1_tag->snps.end(); i++) {
	    if ((*i)->col == (*j)->col)
		found = true;
	}
	//
	// Query locus posses a SNP not present in the catalog.
	//
	if (found == false)
	    return 0;
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

	    if (impute_haplotype(*a, s1_tag->strings, cat_haplotype)) {
		if (s2_tag->alleles.count((*a).c_str()))
		    s2_tag->alleles[cat_haplotype] = s2_tag->alleles[(*a).c_str()];
		else
		    s2_tag->alleles[cat_haplotype] = s2_tag->depth;

		//cerr << s2_tag->id << "; Adding cat haplotype: " << cat_haplotype << " based on depth of " << *a << ", " << s2_tag->alleles[cat_haplotype] << "\n";
		s2_tag->add_match(s1_tag->id, cat_haplotype);
		matches++;
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

int impute_haplotype(string query_haplotype, 
		     vector<pair<allele_type, string> > &cat_haplotypes, 
		     string &match) {

    uint max_len = query_haplotype.length() > cat_haplotypes[0].first.length() ?
	query_haplotype.length() : 
	cat_haplotypes[0].first.length();

    //cerr << "Query len: " << query_haplotype.length() << "; Max length: " << max_len << "\n";

    vector<string> cur, next;

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
	    //cerr << "Comparing query[" << j << "]: '" << query_haplotype[j] << "' to catalog '" << cur[i][j] << "'\n"; 
	    if (query_haplotype[j] == cur[i][j]) {
		//cerr << "  Keeping this haplotype.\n";
		next.push_back(cur[i]);
	    }
	}
	cur = next;
	j++;
    }

    //
    // If, after examining all the available SNPs in the query haplotype, there is
    // still more than a single possible catalog haplotype, then we can't impute it.
    //
    if (cur.size() == 1) {
	match = cur[0];
	return 1;
    }

    return 0;
}

int generate_query_haplotypes(Locus *s1_tag, QLocus *s2_tag, set<string> &query_haplotypes) {
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
    sort(merged_snps.begin(), merged_snps.end(), compare_pair);

    map<string, int>::iterator b;
    string old_allele, new_allele;
    int    pos;

    for (b = s2_tag->alleles.begin(); b != s2_tag->alleles.end(); b++) {
	old_allele = b->first;
	new_allele = "";
	pos        = 0;

	for (k = merged_snps.begin(); k != merged_snps.end(); k++) {
	    if ((*k).first == "catalog") {
		new_allele += s2_tag->con[(*k).second->col];
	    } else {
		new_allele += old_allele[pos];
		pos++;
	    }
	}
	query_haplotypes.insert(new_allele);
	//cerr << "Adding haplotype: " << new_allele << "\n";
    }

    if (s2_tag->alleles.size() == 0) {
	new_allele = "";
	for (k = merged_snps.begin(); k != merged_snps.end(); k++)
	    new_allele += s2_tag->con[(*k).second->col];
	query_haplotypes.insert(new_allele);
	//cerr << "Adding haplotype: " << new_allele << "\n";
    }

    return 0;
}

int find_matches_by_sequence(map<int, Locus *> &sample_1, map<int, QLocus *> &sample_2) {
    map<int, QLocus *>::iterator i;
    int k, min_tag_len;
    QLocus *tag;

    //
    // We don't assume all radtags will be the same length, we allow the sample_2 tags
    // to be shorter than the sample_1 tags (generally the catalog).
    //
    min_tag_len = strlen(sample_1.begin()->second->con); 
    //strlen(i->second->con) > strlen(j->second->con) ?
    //  strlen(j->second->con) : strlen(i->second->con); 

    //
    // Build a hash map out of the first sample (usually the catalog)
    //
    HashMap sample_1_map;
    populate_hash(sample_1, sample_1_map);

    // OpenMP can't parallelize random access iterators, so we convert
    // our map to a vector of integer keys.
    vector<int> keys;
    for (i = sample_2.begin(); i != sample_2.end(); i++) 
	keys.push_back(i->first);

    #pragma omp parallel private(tag, k)
    {
        #pragma omp for schedule(dynamic) 
	for (k = 0; k < (int) keys.size(); k++) {
            tag = sample_2[keys[k]];	    

            //
            // Iterate through the possible alleles for this tag in sample_2
            //
            HashMap::iterator hit;
            vector<pair<allele_type, string> >::iterator q;  // Query records allele_type/search string pairs
            vector<pair<int, allele_type> >::iterator c;     // Hash map records id/allele_type pairs

            for (q = tag->strings.begin(); q != tag->strings.end(); q++) {
                //cerr << "Looking for sequence " << s->second << "\n";
                hit = sample_1_map.find(q->second.c_str());

                if (hit != sample_1_map.end()) {
                    //cerr << "  Found a match for " << hit->first << "\n";

                    //
                    // Found one or more matches between a set of alleles. Verify it is legitamate before recording it.
                    //
                    for (c = hit->second.begin(); c != hit->second.end(); c++)
                        if (verify_sequence_match(sample_1[c->first], tag, min_tag_len)) {
                            tag->add_match(c->first, c->second);
                        }
                }
            }
        }
    }

    return 0;
}

int verify_sequence_match(Locus *s1_tag, QLocus *s2_tag, uint min_tag_len) {
    vector<SNP *>::iterator i, j;
    //
    // We have found a match between two alleles of s1 and s2. We now want
    // to verify that the alleles are consistent between the tags, i.e. they 
    // have the same number and types of SNPs.
    //
    // 1. Make sure s2_tag has the same number of SNPs as s1_tag (s1_tag 
    //    is generally the catalog, so s1_tag may have more SNPs than
    //    s2_tag, but not less).
    //
    if (s2_tag->snps.size() > s1_tag->snps.size()) {
	//cerr << "Match not verified, too many SNPs for " << s2_tag->id << "\n";
	return 0;
    }

    //
    // 2. Okay, s2_tag has the same or less number of SNPs as s1_tag. Now make
    //    sure they occur at the same columns and have the same alternate 
    //    nucleotides.
    //
    uint snp_count = 0;

    for (i = s2_tag->snps.begin(); i != s2_tag->snps.end(); i++) {
	for (j = s1_tag->snps.begin(); j != s1_tag->snps.end(); j++) {
            //
            // SNP occurs in a column that is beyond the length of the s2_tag
            //
            if ((*i)->col > min_tag_len - 1) {
                snp_count++;
                continue;
            }
	    //
	    // SNP occurs in the same column and has the same alternate nucleotides.
	    //
	    if ((*i)->col == (*j)->col) {
		if (((*i)->rank_1 == (*j)->rank_1 &&
		     (*i)->rank_2 == (*j)->rank_2) ||
		    ((*i)->rank_2 == (*j)->rank_1 &&
		     (*i)->rank_1 == (*j)->rank_2))
		    snp_count++;
	    }
	}
    }

    // We matched each SNP in s2_tag against one in s1_tag
    if (snp_count == s2_tag->snps.size())
	return 1;

    //cerr << "Match not verified, SNPs in S2 do not match catalog: " << snp_count << " vs " << s2_tag->snps.size() << "\n";

    return 0;
}

int populate_hash(map<int, Locus *> &sample, HashMap &hash_map) {
    map<int, Locus *>::iterator it;
    vector<pair<allele_type, string> >::iterator all_it;
    Locus *tag;

    //
    // Create a hash map out of the set of alleles for each Unique Tag.
    //
    for (it = sample.begin(); it != sample.end(); it++) {
        tag = it->second;

        for (all_it = tag->strings.begin(); all_it != tag->strings.end(); all_it++)
            hash_map[all_it->second.c_str()].push_back(make_pair(tag->id, all_it->first));
    }

    //dump_kmer_map(kmer_map);

    return 0;
}

int write_matches(map<int, QLocus *> &sample) {
    map<int, QLocus *>::iterator i;
    vector<pair<int, allele_type> >::iterator s;

    //
    // Parse the input file names to create the output file
    //
    size_t pos_1    = sample_2_file.find_last_of("/");
    string out_file = out_path + sample_2_file.substr(pos_1 + 1)  + ".matches.tsv";

    // Open the output files for writing.
    std::ofstream matches(out_file.c_str());
    string type;
    uint   match_depth;
    cerr << "Outputing to file " << out_file.c_str() << "\n";

    for (i = sample.begin(); i != sample.end(); i++) {

	for (s = (*i).second->matches.begin(); s != (*i).second->matches.end(); s++) {
	    match_depth = 
		i->second->alleles.count(s->second) > 0 ? 
		i->second->alleles[s->second] : i->second->depth;
	    matches << 
		"0"           << "\t" <<
		batch_id      << "\t" <<
		s->first      << "\t" <<
		samp_id       << "\t" <<
		i->second->id << "\t" << 
		s->second     << "\t" <<
		match_depth   << "\n";
	}
    }

    matches.close();

    return 0;
}

bool compare_pair(pair<string, SNP *> a, pair<string, SNP *> b) {
    return (a.second->col < b.second->col);
}

int parse_command_line(int argc, char* argv[]) {
    int c;
     
    while (1) {
	static struct option long_options[] = {
	    {"help",        no_argument,       NULL, 'h'},
            {"version",     no_argument,       NULL, 'v'},
	    {"num_threads", required_argument, NULL, 'p'},
	    {"batch_id",    required_argument, NULL, 'b'},
	    {"sample_1",    required_argument, NULL, 'r'},
	    {"sample_1_id", required_argument, NULL, 'R'},
	    {"sample_2",    required_argument, NULL, 's'},
	    {"sample_2_id", required_argument, NULL, 'S'},
	    {"genomic_loc", no_argument,       NULL, 'g'},
	    {"verify_hap",  no_argument,       NULL, 'x'},
	    {"outpath",     required_argument, NULL, 'o'},
	    {0, 0, 0, 0}
	};
	
	// getopt_long stores the option index here.
	int option_index = 0;
     
	c = getopt_long(argc, argv, "hgxvr:s:c:o:R:S:b:p:", long_options, &option_index);
     
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
	    batch_id = atoi(optarg);
	    break;
	case 'r':
	    sample_1_file = optarg;
	    catalog = 0;
	    break;
	case 's':
	    sample_2_file = optarg;
	    break;
	case 'S':
	    samp_id = atoi(optarg);
	    break;
	case 'g':
	    search_type = genomic_loc;
	    break;
     	case 'o':
	    out_path = optarg;
	    break;
	case 'c': 
	    sample_1_file = optarg;
	    catalog++;
	    break;
	case 'x': 
	    verify_haplotypes = false;
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

    if (sample_1_file.length() == 0) {
	cerr << "You must specify two sample files.\n";
	help();
    }

    if (sample_2_file.length() == 0) {
	cerr << "You must specify two sample files.\n";
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
              << "sstacks -b batch_id -c catalog_file -s sample_file [-S id] [-r sample_file] [-o path] [-p num_threads] [-g] [-x] [-v] [-h]" << "\n"
              << "  p: enable parallel execution with num_threads threads.\n"
	      << "  b: MySQL ID of this batch." << "\n"
	      << "  c: TSV file from which to load the catalog RAD-Tags." << "\n"
	      << "  r: Load the TSV file of a single sample instead of a catalog." << "\n"
	      << "  s: TSV file from which to load sample RAD-Tags." << "\n"
	      << "  S: MySQL ID of the specified sample." << "\n"
	      << "  o: output path to write results." << "\n"
              << "  g: base matching on genomic location, not sequence identity." << "\n"
	      << "  x: don't verify haplotype of matching locus." << "\n"
	      << "  v: print program version." << "\n"
	      << "  h: display this help messsage." << "\n\n";

    exit(0);
}
