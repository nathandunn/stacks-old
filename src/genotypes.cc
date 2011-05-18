// -*-mode:c++; c-style:k&r; c-basic-offset:4;-*-
//
// Copyright 2011, Julian Catchen <jcatchen@uoregon.edu>
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
// genotypes -- genotype a mapping cross.
//
// Julian Catchen
// jcatchen@uoregon.edu
// University of Oregon
//
// $Id: genotypes.cc 2117 2011-05-07 00:54:22Z catchen $
//

#include "genotypes.h"

// Global variables to hold command-line options.
int       num_threads = 1;
int       batch_id    = 0;
map_types map_type    = none;
out_types out_type    = joinmap;
string    in_path;
string    out_path;
string    out_file;
string    bl_file;
string    wl_file;
int       progeny_limit   = 1;
bool      corrections     = false;
bool      expand_id       = false;
bool      sql_out         = false;
int       min_stack_depth = 0;
int       min_hom_seqs    = 5;
double    min_het_seqs    = 0.05;
double    max_het_seqs    = 0.1;

set<int> whitelist, blacklist;

//
// Dictionaries to hold legal genotypes for different map types.
//
map<string, map<string, string> > global_dictionary;

int main (int argc, char* argv[]) {

    parse_command_line(argc, argv);

    //
    // Set the number of OpenMP parallel threads to execute.
    //
    #ifdef _OPENMP
    omp_set_num_threads(num_threads);
    #endif

    initialize_dictionaries(global_dictionary);

    vector<string> files;
    build_file_list(files);

    if (wl_file.length() > 0) {
	load_marker_list(wl_file, whitelist);
	cerr << "Loaded " << whitelist.size() << " whitelisted markers.\n";
    }
    if (bl_file.length() > 0) {
	load_marker_list(bl_file, blacklist);
	cerr << "Loaded " << blacklist.size() << " blacklisted markers.\n";
    }

    //
    // Load the catalog
    //
    stringstream catalog_file;
    map<int, CLocus *> catalog;
    int res;
    catalog_file << in_path << "batch_" << batch_id << ".catalog";
    if ((res = load_loci(catalog_file.str(), catalog, false)) == 0) {
    	cerr << "Unable to load the catalog '" << catalog_file.str() << "'\n";
     	return 0;
    }

    set<int> parent_ids;
    identify_parental_ids(catalog, parent_ids);

    //
    // Implement the black/white list
    //
    reduce_catalog(catalog, whitelist, blacklist);

    //
    // Load matches to the catalog
    //
    vector<vector<CatMatch *> > catalog_matches;
    map<int, string>            samples;
    vector<int>                 sample_ids;
    for (uint i = 0; i < files.size(); i++) {
	vector<CatMatch *> m;
	load_catalog_matches(in_path + files[i], m);
	catalog_matches.push_back(m);
	samples[m[0]->sample_id] = files[i];
	sample_ids.push_back(m[0]->sample_id);
    }

    //
    // Create the population map
    //
    cerr << "Populating observed haplotypes for " << sample_ids.size() << " samples, " << catalog.size() << " loci.\n";
    PopMap<CLocus> *pmap = new PopMap<CLocus>(sample_ids.size(), catalog.size());
    pmap->populate(sample_ids, catalog, catalog_matches, min_stack_depth);

    //
    // Identify mappable markers in the parents
    //
    find_markers(catalog, pmap, parent_ids);

    //
    // Calculate F, inbreeding coefficient
    //
    calculate_f(catalog, pmap, parent_ids);

    //
    // Create genotypes maps
    //
    map<int, CLocus *>::iterator it;
    for (it = catalog.begin(); it != catalog.end(); it++) {
	if (it->second->marker.length() > 0) {
	    create_genotype_map(it->second, pmap, parent_ids);
	    call_population_genotypes(it->second, pmap, global_dictionary);
	}
    }

    //
    // Make automated corrections
    // 
    if (corrections)
	automated_corrections(samples, parent_ids, catalog, catalog_matches, pmap);

    //
    // If a mapping type was specified, output it.
    //
    switch(map_type) {
    case dh:
    	export_dh_map(catalog, pmap, parent_ids, samples);
    	break;
    case cp:
    	export_cp_map(catalog, pmap, parent_ids, samples);
    	break;
    case bc1:
    	export_bc1_map(catalog, pmap, parent_ids, samples);
    	break;
    case f2:
    	export_f2_map(catalog, pmap, parent_ids, samples);
    	break;
    case gen:
	export_gen_map(catalog, pmap, parent_ids, samples);
	break;
    case none:
    case unk:
	break;
    }

    //
    // Output the observe haplotypes.
    //
    write_generic(catalog, pmap, samples, parent_ids, false);

    return 0;
}

int reduce_catalog(map<int, CLocus *> &catalog, set<int> &whitelist, set<int> &blacklist) {
    map<int, CLocus *> list;
    map<int, CLocus *>::iterator it;
    CLocus *loc;

    if (whitelist.size() == 0 && blacklist.size() == 0) 
	return 0;
 
    for (it = catalog.begin(); it != catalog.end(); it++) {
	loc = it->second;

	if (whitelist.size() > 0 && whitelist.count(loc->id) == 0) continue;
	if (blacklist.count(loc->id)) continue;
	list[it->first] = it->second;
    }

    catalog = list;

    return 0;
}

int identify_parental_ids(map<int, CLocus *> &catalog, set<int> &parents) {
    map<int, CLocus *>::iterator it;
    CLocus *loc;
    int     sample_id;

    for (it = catalog.begin(); it != catalog.end(); it++) {
	loc = it->second;
	for (uint i = 0; i < loc->comp.size(); i++) {
	    sample_id = (int) strtol(loc->comp[i], NULL, 10);
	}
	parents.insert(sample_id);
    }

    set<int>::iterator sit;
    cerr << "Identified parent IDs: ";
    for (sit = parents.begin(); sit != parents.end(); sit++)
	cerr << *sit << " ";
    cerr << "\n";

    return 0;
}

int find_markers(map<int, CLocus *> &catalog, PopMap<CLocus> *pmap, set<int> &parent_ids) {
    map<int, CLocus *>::iterator it;
    vector<char *>::iterator hit;
    set<int>::iterator p, q;
    int pid_1, pid_2, parent_count, allele_cnt_1, allele_cnt_2;
    Datum  *d_1, *d_2;
    CLocus *loc;

    for (it = catalog.begin(); it != catalog.end(); it++) {
	loc = it->second;
	//
	// Count the number of parental tags matching this catalog tag. A proper marker should
	// contain a single representative from each parent; multiple alleles must be called from 
	// a single tag from a single parent.
	//
	if (parent_ids.size() == 1) {
	    p = parent_ids.begin();
	    pid_1 = *p;
	    pid_2 = -1;
	    d_1   = pmap->blacklisted(loc->id, pid_1) ? NULL : pmap->datum(loc->id, pid_1);
	    d_2   = NULL;
	} else {
	    p = parent_ids.begin();
	    q = p++;

	    pid_1 = *p < *q ? *p : *q;
	    pid_2 = *p < *q ? *q : *p;
	    if (pmap->blacklisted(loc->id, pid_1) ||
		pmap->blacklisted(loc->id, pid_2)) {
		d_1 = NULL;
		d_2 = NULL;
	    } else {
		d_1   = pmap->datum(loc->id, pid_1);
		d_2   = pmap->datum(loc->id, pid_2);
	    }
	}

	parent_count = 0;
	if (d_1 != NULL) parent_count++;
	if (d_2 != NULL) parent_count++;

	//
	// Locus is present in both parents.
	//
	if (parent_count == 2) {
	    allele_cnt_1 = d_1->obshap.size();
	    allele_cnt_2 = d_2->obshap.size();

            //
            // Determine the number of unique alleles
            //
	    set<string> unique_alleles;

	    for (hit = d_1->obshap.begin(); hit != d_1->obshap.end(); hit++)
                unique_alleles.insert(*hit);
	    for (hit = d_2->obshap.begin(); hit != d_2->obshap.end(); hit++)
                unique_alleles.insert(*hit);
            int num_unique_alleles = unique_alleles.size();

	    //
	    // Locus is heterozygous in both parents. However, the number of alleles present distinguishes 
            // what type of marker it is. Four unique alleles requries an ab/cd marker, while four 
            // alleles that are the same in both parents requires an ab/ab marker. Finally, three unique 
            // alleles requires either an ab/ac marker.
	    //
	    if (allele_cnt_1 == 2 && allele_cnt_2 == 2) {
                if (num_unique_alleles == 3)
                    loc->marker = "ab/ac";
		else if (num_unique_alleles == 2)
                    loc->marker = "ab/ab";
		else
                    loc->marker = "ab/cd";
	    //
	    // Locus is homozygous in one parent and heterozygous in the other.
	    //
	    } else if (allele_cnt_1 == 2 && allele_cnt_2 == 1) {

                if (num_unique_alleles == 3)
                    loc->marker = "ab/cc";
                else if (num_unique_alleles == 2)
                    loc->marker = "ab/aa";
	    //
	    // Locus is homozygous in one parent and heterozygous in the other.
	    //
	    } else if (allele_cnt_1 == 1 && allele_cnt_2 == 2) {

                if (num_unique_alleles == 3)
                    loc->marker = "cc/ab";
                else if (num_unique_alleles == 2)
                    loc->marker = "aa/ab";
            //
            // Locus is homozygous in both parents, but heterozygous between parents.
            //
	    } else if (allele_cnt_1 == 1 && allele_cnt_2 == 1) {

		if (strcmp(d_1->obshap[0], "consensus") != 0 &&
		    strcmp(d_2->obshap[0], "consensus") != 0)
                    loc->marker = "aa/bb";
	    }
        //
        // Locus only exists in one parent.
	//
	} else if (parent_count == 1) {
	    if (d_1 != NULL && d_1->obshap.size() == 2)
		loc->marker = "ab/--";
	    else if (d_2 != NULL && d_2->obshap.size() == 2)
		loc->marker = "--/ab";
	}
    }

    return 0;
}

int calculate_f(map<int, CLocus *> &catalog, PopMap<CLocus> *pmap, set<int> &parent_ids) {
    map<int, CLocus *>::iterator it;
    map<char, int>::iterator j;
    Datum **d;
    CLocus *loc;

    for (it = catalog.begin(); it != catalog.end(); it++) {
	loc = it->second;
	d   = pmap->locus(loc->id);

	if (loc->snps.size() == 0) continue;

	double tot  = 0.0;
	double hets = 0;
	double p, q, h, h0;
	map<char, int> alle;

	for (int i = 0; i < pmap->sample_cnt(); i++) {
	    if (d[i] == NULL) 
		continue;
	    if (parent_ids.count(pmap->rev_sample_index(i))) 
		continue;

	    tot++;

	    if (d[i]->obshap.size() > 1) hets++;

	    //
	    // We are measuring the first SNP in the haplotype
	    //
	    for (uint j = 0; j < d[i]->obshap.size(); j++)
		alle[d[i]->obshap[j][0]]++;
	}

	if (alle.size() > 2 || tot == 0)
	    continue;

	j  = alle.begin();
	p  = j->second;
	j++;
	q  = j->second;
	h  = hets / tot;            // Observed frequency of heterozygotes in the population
	h0 = 2 * (p/tot) * (q/tot); // 2PQ, expected frequency of hets under Hardy-Weinberg 
	if (h0 > 0)
	    loc->f = (h0 - h) / h0;

	//cerr << "P: " << p << "; Q: " << q << "; Hets: " << hets << "; total: " << tot << "; f: " << loc->f << "\n";
    }

    return 0;
}

int create_genotype_map(CLocus *locus, PopMap<CLocus> *pmap, set<int> &parent_ids) {
    //
    // Create a genotype map. For any set of alleles, this routine will
    // assign each allele to one of the constituent genotypes, e.g. given the 
    // marker type 'aaxbb' and the alleles 'A' from the male, and 'G'
    // from the female, will assign 'G' == 'bb' and 'A'== 'aa'. It assumes that 
    // recombination may have occurred as with an F2, F3 or later cross.
    //
    //cerr << "Creating genotype map for catalog ID " << locus->id  << ", marker: " << locus->marker << ".\n";

    //
    // Get the parent IDs ordered
    //
    set<int>::iterator p = parent_ids.begin();
    set<int>::iterator q = p++;
    int pid_1 = *p < *q ? *p : *q;
    int pid_2 = *p < *q ? *q : *p;

    set<char> p1_gtypes, p2_gtypes;
    set<char>::iterator i;
    map<char, int> legal_gtypes, com_gtypes;
    //
    // First, identify any alleles that are common between the two parents.
    //
    p1_gtypes.insert(locus->marker[0]);
    p1_gtypes.insert(locus->marker[1]);
    p2_gtypes.insert(locus->marker[3]);
    p2_gtypes.insert(locus->marker[4]);
    for (i = p1_gtypes.begin(); i != p1_gtypes.end(); i++)
	if (*i != '-') legal_gtypes[*i]++;
    for (i = p2_gtypes.begin(); i != p2_gtypes.end(); i++)
	if (*i != '-') legal_gtypes[*i]++;
    //
    // Find the common genotypes
    //
    vector<char> types;
    map<char, int>::iterator j;

    for (j = legal_gtypes.begin(); j != legal_gtypes.end(); j++)
	if (j->second > 1) types.push_back(j->first);
    sort(types.begin(), types.end());

    Datum *d_1, *d_2;
    map<string, int> haplotypes;
    map<string, int>::iterator k;
    vector<pair<string, int> > sorted_haplotypes;
    d_1 = pmap->datum(locus->id, pid_1);
    d_2 = pmap->datum(locus->id, pid_2);

    if (d_1 != NULL) {
	for (uint n = 0; n < d_1->obshap.size(); n++)
	    haplotypes[d_1->obshap[n]]++;
    }
    if (d_2 != NULL) {
	for (uint n = 0; n < d_2->obshap.size(); n++)
	    haplotypes[d_2->obshap[n]]++;
    }

    // 
    // Sort the haplotypes map by value
    //
    for (k = haplotypes.begin(); k != haplotypes.end(); k++)
	sorted_haplotypes.push_back(*k);
    sort(sorted_haplotypes.begin(), sorted_haplotypes.end(), hap_compare);

    for (uint n = 0, index = 0; n < sorted_haplotypes.size() && index < types.size(); n++, index++) {
	if (sorted_haplotypes[n].second > 1) {
	    locus->gmap[sorted_haplotypes[n].first] = types[index];
	    com_gtypes[types[index]]++;
	    //cerr << "  Assinging common allele " << sorted_haplotypes[n].first << " to genotype '" << locus->gmap[sorted_haplotypes[n].first] << "'\n";
	}
    }

    //
    // Now, examine the remaining first parent alleles.
    //
    if (d_1 != NULL) {
	legal_gtypes.clear();
	for (i = p1_gtypes.begin(); i != p1_gtypes.end(); i++)
	    if (*i != '-' && com_gtypes.count(*i) == 0) {
		//cerr << "  Adding " << *i << " to first parent genotypes\n";
		legal_gtypes[*i]++;
	    }
	types.clear();
	for (j = legal_gtypes.begin(); j != legal_gtypes.end(); j++)
	    types.push_back(j->first);
	sort(types.begin(), types.end());

	for (uint n = 0, index = 0; n < d_1->obshap.size() && index < types.size(); n++, index++) {
	    if (locus->gmap.count(d_1->obshap[n])) {
		index--;
		continue;
	    }
	    locus->gmap[d_1->obshap[n]] = types[index];
	    //cerr << "  Assinging '" << d_1->obshap[n] << "' to first parent genotype '" << locus->gmap[d_1->obshap[n]] << "'\n";
	}
    }

    //
    // Finally, repeat in the second parent.
    //
    if (d_2 != NULL) {
	legal_gtypes.clear();
	for (i = p2_gtypes.begin(); i != p2_gtypes.end(); i++)
	    if (*i != '-' && com_gtypes.count(*i) == 0) {
		//cerr << "  Adding " << *i << " to second genotypes\n";
		legal_gtypes[*i]++;
	    }
	types.clear();
	for (j = legal_gtypes.begin(); j != legal_gtypes.end(); j++)
	    types.push_back(j->first);
	sort(types.begin(), types.end());

	for (uint n = 0, index = 0; n < d_2->obshap.size() && index < types.size(); n++, index++) {
	    if (locus->gmap.count(d_2->obshap[n])) {
		index--;
		continue;
	    }
	    locus->gmap[d_2->obshap[n]] = types[index];
	    //cerr << "  Assinging '" << d_2->obshap[n] << "' to second parent genotype '" << locus->gmap[d_2->obshap[n]] << "'\n";
	}
    }

    return 0;
}

int call_population_genotypes(CLocus *locus, 
			      PopMap<CLocus> *pmap, 
			      map<string, map<string, string> > &dictionary) {
    //
    // Fetch the array of observed haplotypes from the population
    //
    Datum **d = pmap->locus(locus->id);

    for (int i = 0; i < pmap->sample_cnt(); i++) {
	if (d[i] == NULL) 
	    continue;

	vector<string> gtypes;
	string gtype;

	for (uint j = 0; j < d[i]->obshap.size(); j++) {
	    //
	    // Impossible allele encountered.
	    //
	    if (locus->gmap.count(d[i]->obshap[j]) == 0) {
		gtypes.clear();
		gtypes.push_back("-");
		goto impossible;
	    }

	    gtypes.push_back(locus->gmap[d[i]->obshap[j]]);
	    //cerr << "Genotype: " << locus->gmap[d[i]->obshap[j]] << "\n";
	}

    impossible:
	sort(gtypes.begin(), gtypes.end());
 	for (uint j = 0; j < gtypes.size(); j++) {
	    gtype += gtypes[j];
	    //cerr << "  Adding genotype to string: " << gtypes[j] << "; " << gtype << "\n";
	}

 	string m = dictionary[locus->marker].count(gtype) ? 
	    dictionary[locus->marker][gtype] : 
	    "-";

	d[i]->gtype = new char[m.length() + 1];
	strcpy(d[i]->gtype, m.c_str());

	if (m != "-")
	    locus->gcnt++;

	//cerr << "Assigning datum, marker: " << locus->marker << ", string: " << m << ", haplotype: " << d[i]->obshap[0] << ", gtype: " << gtype << "\n";
     }

    return 0;
}

int automated_corrections(map<int, string> &samples, set<int> &parent_ids, map<int, CLocus *> &catalog,
			  vector<vector<CatMatch *> > &matches, PopMap<CLocus> *pmap) {
    int sample_id, catalog_id, tag_id;
    Datum *d;
    Locus *s;
    string file;

    cerr << "Performing automated corrections...\n";

    for (uint i = 0; i < matches.size(); i++) {
	sample_id = matches[i][0]->sample_id;
	file      = samples[sample_id];

	if (parent_ids.count(sample_id)) continue;

	map<int, Locus *> stacks;
	int res;
	if ((res = load_loci(in_path + file, stacks, true)) == 0) {
	    cerr << "Unable to load sample file '" << file << "'\n";
	    return 0;
	}

	set<pair<int, int> > processed;

	for (uint j = 0; j < matches[i].size(); j++) {
	     catalog_id = matches[i][j]->cat_id;
	     sample_id  = matches[i][j]->sample_id;
	     tag_id     = matches[i][j]->tag_id;

	     if (catalog.count(catalog_id) == 0) continue;

	     //
	     // There are multiple matches per stack, but we only need to process
	     // each stack once to make corrections.
	     //
	     if (processed.count(make_pair(catalog_id, tag_id)) == 0 &&
		 catalog[catalog_id]->marker.length() > 0) {
		 d = pmap->datum(catalog_id, sample_id);

		 if (d != NULL && strcmp(d->gtype, "-") != 0) {
		     s = stacks[tag_id];
		     check_uncalled_snps(catalog[catalog_id], s, d);
		 }

		 processed.insert(make_pair(catalog_id, tag_id));
	     }
	}

	//
	// Free up memory
	//
	map<int, Locus *>::iterator it;
	for (it = stacks.begin(); it != stacks.end(); it++)
	    delete it->second;
    }

    //
    // Summarize correction results
    //
    long pot_gen = 0;
    long tot_gen = 0;
    long tot_cor = 0;
    long het_cor = 0;
    long rem_cor = 0;
    int  markers = 0;
    map<int, CLocus *>::iterator it;
    for (it = catalog.begin(); it != catalog.end(); it++) {
	if (it->second->marker.length() == 0) continue;
	markers++;

	for (int j = 0; j < pmap->sample_cnt(); j++)  {
	    sample_id = pmap->rev_sample_index(j);
	    if (parent_ids.count(sample_id)) continue;

	    d = pmap->datum(it->first, sample_id);

	    pot_gen++;

	    if (d == NULL) continue;

	    tot_gen++;

	    if (d->corrected == true) {
		tot_cor++;

		if (d->gtype[0] == '-')
		    rem_cor++;
		else 
		    het_cor++;
	    }
	}
    }
    cerr << pot_gen << " potential genotypes in " << markers << " markers, " 
	 << tot_gen << " populated; " 
	 << tot_cor << " corrected, " 
	 << het_cor << " converted to heterozygotes, " 
	 << rem_cor << " unsupported homozygotes removed.\n";

    return 0;
}

int check_uncalled_snps(CLocus *clocus, Locus *stack, Datum *d) {
    //
    // If this locus is already known to be multi-allele, return, we only want 
    // to verify uncalled SNPs.
    //
    if (strlen(d->gtype) > 1 && d->gtype[0] != d->gtype[1])
	return 0;

    //cerr << "Catalog locus: " << clocus->id << ", marker: " << clocus->marker << ", tag_id: " << stack->id << "; Starting gtype: " << d->gtype << "\n";

    vector<SNP *> verified_snps;
    string status = "false";
    string homozygous;

    for (uint i = 0; i < clocus->snps.size(); i++) {
	check_homozygosity(stack->reads, 
			   clocus->snps[i]->col, clocus->snps[i]->rank_1, clocus->snps[i]->rank_2, 
			   homozygous);

 	if (homozygous == "unknown")
 	    status = "true";
	else if (homozygous == "false")
 	    verified_snps.push_back(clocus->snps[i]);
     }

    if (status == "true") {
	d->corrected = true;
	delete [] d->gtype;
	d->gtype = new char[2];
	strcpy(d->gtype, "-");
	return 0;
    } else if (verified_snps.size() < clocus->snps.size()) {
	return 0;
    }

    //
    // Detect the alleles present from the verified SNPs
    //
    vector<string> haplotypes;
    call_alleles(verified_snps, stack->reads, haplotypes);

    vector<string> types;
    for (uint i = 0; i < haplotypes.size(); i++) {
	if (clocus->gmap.count(haplotypes[i])) {
	    //cerr << "    Adding haplotype '" << haplotypes[i] << "', which maps to '" << clocus->gmap[haplotypes[i]] << "' to the genotype\n";
	    types.push_back(clocus->gmap[haplotypes[i]]);
	} else {
	    types.push_back("-");
	}
    }

    sort(types.begin(), types.end());

    string genotype;
    for (uint i = 0; i < types.size(); i++)
	genotype += types[i];

    genotype = 
	global_dictionary[clocus->marker].count(genotype) ? 
	global_dictionary[clocus->marker][genotype] : 
	"-";

    if (genotype != "-" && strcmp(genotype.c_str(), d->gtype) != 0) {
	d->corrected = true;
	delete [] d->gtype;
	d->gtype = new char[genotype.length() + 1];
	strcpy(d->gtype, genotype.c_str());
    }
    //cerr << "  Catalog: " << clocus->id << ", stack: " << stack->id << ", Ending Genotype: " << d->gtype << "\n\n";

    return 0;
}

int call_alleles(vector<SNP *> &snps, vector<char *> &reads, vector<string> &haplotypes) {
    map<string, int> a;
    int  height = reads.size();
    char base;

    for (int i = 0; i < height; i++) {
	string haplotype;

	for (uint j = 0; j < snps.size(); j++) {
	    base = reads[i][snps[j]->col];

	    //
	    // Check to make sure the nucleotide at the location of this SNP is
	    // of one of the two possible states the multinomial model called.
	    //
	    if (base == snps[j]->rank_1 || base == snps[j]->rank_2) 
		haplotype += base;
	    else
		break;
	}

	if (haplotype.length() == snps.size())
	    a[haplotype]++;
    }

    map<string, int>::iterator it;
    for (it = a.begin(); it != a.end(); it++)
	haplotypes.push_back(it->first);

    return 0;
}

int check_homozygosity(vector<char *> &reads, int col, char rank_1, char rank_2, string &homozygous) {

    //cerr << "  Examining col " << col << ", rank 1: " << rank_1 << "; rank 2: " << rank_2 << "\n";

    int height = reads.size();
    homozygous = "true";

    if (height < min_hom_seqs) {
	homozygous = "unknown";
	return 0;
    }

    map<char, int> nuc;
    vector<pair<char, int> > sorted_nuc;

    nuc['A'] = 0;
    nuc['C'] = 0;
    nuc['G'] = 0;
    nuc['T'] = 0;

    for (int j = 0; j < height; j++)
	nuc[reads[j][col]]++;

    map<char, int>::iterator i;
    for (i = nuc.begin(); i != nuc.end(); i++)
	sorted_nuc.push_back(make_pair(i->first, i->second));

    sort(sorted_nuc.begin(), sorted_nuc.end(), compare_pair);
    
    //
    // Check if more than a single nucleotide occurs in this column. Only 
    // count nucleotides that are part of the called SNP, do not count 
    // error-generated nucleotides.
    //
    if ((sorted_nuc[0].second > 0) && 
	(sorted_nuc[0].first == rank_1 || sorted_nuc[0].first == rank_2) &&
	(sorted_nuc[1].second > 0) && 
	(sorted_nuc[1].first == rank_1 || sorted_nuc[1].first == rank_2)) {
	homozygous = "false";
    }

    //
    // If we find a second nucleotide present, check its prevelence. If it is 
    // less than 1/20 of the total reads, don't count a heterozygote. If it is 
    // less than 1/10 report that we can't tell if its homozygous or not. Otherwise,
    // report this tag as a heterozygote.
    //
    double frac = (double) sorted_nuc[1].second / (double) height;
    //cerr << "    Frac: " << frac << "; Second-most Prominent Nuc count: " << sorted_nuc[1].second << "; Depth: " << height << "\n";

    if (homozygous == "false" && frac < min_het_seqs)
	homozygous = "true";
    else if (homozygous == "false" && frac < max_het_seqs)
	homozygous = "unknown";

    return 0;
}

int export_gen_map(map<int, CLocus *> &catalog, PopMap<CLocus> *pmap, set<int> &parent_ids, map<int, string> &samples) {
    //
    // We wish to export, a set of generic genotypes, not specific to any mapping type.
    // 

    //
    // Output the results
    //
    write_generic(catalog, pmap, samples, parent_ids, true);

    if (sql_out)
	write_sql(catalog, pmap, parent_ids);

    return 0;
}

int export_f2_map(map<int, CLocus *> &catalog, PopMap<CLocus> *pmap, set<int> &parent_ids, map<int, string> &samples) {
    //
    // We wish to export, according to the JoinMap manual, a locus genotype file (loc-file), 
    // which contains the information of all the loci for a single segregating population.
    // 
    // We are exporting an F2 population type:
    // The result of selfing the F1 of a cross between two fully homozygous diploid parents.
    //
    // Genotype codes for an F2 population, depending on the locus segregation type.
    //
    // Seg. type   Possible genotypes
    // ---------   ------------------
    // <aaxbb>     a, b, h, –
    // <abxcd>     a, b, h, –
    // <abxaa>     a, –
    // <aaxab>     b, –
    // <abxcc>     a, b, – 
    // <ccxab>     b, –
    //
    map<string, string> types;
    types["aa/bb"] = "aaxbb";
    types["ab/cd"] = "abxcd";
    types["ab/aa"] = "abxaa";
    types["aa/ab"] = "aaxab";
    types["ab/cc"] = "abxcc";
    types["cc/ab"] = "ccxab";

    map<string, map<string, string> > dictionary;
    dictionary["aaxbb"]["aa"] = "a";
    dictionary["aaxbb"]["ab"] = "h";
    dictionary["aaxbb"]["bb"] = "b";
    dictionary["aaxbb"]["-"]  = "-";

    dictionary["abxcd"]["aa"] = "a";
    dictionary["abxcd"]["ab"] = "a";
    dictionary["abxcd"]["bb"] = "a";
    dictionary["abxcd"]["cc"] = "b";
    dictionary["abxcd"]["cd"] = "b";
    dictionary["abxcd"]["dd"] = "b";
    dictionary["abxcd"]["ac"] = "h";
    dictionary["abxcd"]["ad"] = "h";
    dictionary["abxcd"]["bc"] = "h";
    dictionary["abxcd"]["bd"] = "h";
    dictionary["abxcd"]["-"]  = "-";

    dictionary["abxaa"]["aa"] = "-";
    dictionary["abxaa"]["ab"] = "-";
    dictionary["abxaa"]["bb"] = "a";
    dictionary["abxaa"]["-"]  = "-";

    dictionary["aaxab"]["aa"] = "-";
    dictionary["aaxab"]["ab"] = "-";
    dictionary["aaxab"]["bb"] = "b";
    dictionary["aaxab"]["-"]  = "-";

    dictionary["abxcc"]["a"]  = "a";
    dictionary["abxcc"]["ab"] = "a";
    dictionary["abxcc"]["bb"] = "a";
    dictionary["abxcc"]["cc"] = "b";
    dictionary["abxcc"]["ac"] = "-";
    dictionary["abxcc"]["bc"] = "-";
    dictionary["abxcc"]["-"]  = "-";

    dictionary["ccxab"]["aa"] = "b";
    dictionary["ccxab"]["ab"] = "b";
    dictionary["ccxab"]["bb"] = "b";
    dictionary["ccxab"]["cc"] = "a";
    dictionary["ccxab"]["ac"] = "-";
    dictionary["ccxab"]["bc"] = "-";
    dictionary["ccxab"]["-"]  = "-";

    //
    // Translate the genotypes for this particular map type.
    //
    translate_genotypes(types, dictionary, catalog, pmap, samples, parent_ids);

    //
    // Output the results
    //
    switch(out_type) {
    case joinmap:
	write_joinmap(catalog, pmap, types, samples, parent_ids);
	break;
    case rqtl:
	write_rqtl(catalog, pmap, types, samples, parent_ids);
	break;
    }

    if (sql_out)
	write_sql(catalog, pmap, parent_ids);

    return 0;
}

int export_dh_map(map<int, CLocus *> &catalog, PopMap<CLocus> *pmap, set<int> &parent_ids, map<int, string> &samples) {
    //
    // We wish to export, according to the JoinMap manual, a locus genotype file (loc-file), 
    // which contains the information of all the loci for a single segregating population.
    // 
    // We are exporting a DH population type:
    //   a doubled haploid population: the result of doubling the gametes of a single heterozygous 
    //   diploid individual.
    //
    // Segregation type codes for population type DH, from Joinmap manual:
    //
    // Code      Description
    // -------   -----------
    // <aaxbb>   locus homozygous in both parents, heterozygous between the parents
    //
    // Genotype codes for a CP population, depending on the locus segregation type.
    //
    // Seg. type   Possible genotypes
    // ---------   ------------------
    //     a       the one genotype
    //     b       the other genotype
    //
    map<string, string> types;
    types["ab/--"] = "abx--";
    types["--/ab"] = "--xab";

    map<string, map<string, string> > dictionary;
    dictionary["abx--"]["aa"] = "a";
    dictionary["abx--"]["bb"] = "b";
    dictionary["abx--"]["-"]  = "-";

    dictionary["--xab"]["aa"] = "a";
    dictionary["--xab"]["bb"] = "b";
    dictionary["--xab"]["-"]  = "-";

    //
    // Translate the genotypes for this particular map type.
    //
    translate_genotypes(types, dictionary, catalog, pmap, samples, parent_ids);

    //
    // Output the results
    //
    switch(out_type) {
    case joinmap:
	write_joinmap(catalog, pmap, types, samples, parent_ids);
	break;
    case rqtl:
	write_rqtl(catalog, pmap, types, samples, parent_ids);
	break;
    }

    if (sql_out)
	write_sql(catalog, pmap, parent_ids);

    return 0;
}

int export_bc1_map(map<int, CLocus *> &catalog, PopMap<CLocus> *pmap, set<int> &parent_ids, map<int, string> &samples) {
    //
    // We wish to export, according to the JoinMap manual, a locus genotype file (loc-file), 
    // which contains the information of all the loci for a single segregating population.
    // 
    // We are exporting a BC1 population type:
    //   a first generation backcross population: the result of crossing the F1 of a cross between 
    //   two fully homozygous diploid parents to one of the parents.
    //
    // Segregation type codes for population type BC1, from Joinmap manual:
    //
    // Code      Description
    // -------   -----------
    // <aaxbb>   locus homozygous in both parents, heterozygous between the parents
    //
    // Genotype codes for a CP population, depending on the locus segregation type.
    //
    // Seg. type   Possible genotypes
    // ---------   ------------------
    //     a       homozygote or haploid as the first parent
    //     b       homozygote or haploid as the second parent
    //     h       heterozygote (as the F1)
    //
    map<string, string> types;
    types["aa/bb"] = "aaxbb";
    types["bb/aa"] = "bbxaa";
    types["ab/cc"] = "abxcc";
    types["cc/ab"] = "ccxab";

    map<string, map<string, string> > dictionary;
    dictionary["aaxbb"]["-"]  = "-";
    dictionary["aaxbb"]["aa"] = "b";
    dictionary["aaxbb"]["ab"] = "h";
    dictionary["aaxbb"]["bb"] = "h";

    dictionary["bbxaa"]["-"]  = "-";
    dictionary["bbxaa"]["aa"] = "h";
    dictionary["bbxaa"]["ab"] = "h";
    dictionary["bbxaa"]["bb"] = "a";

    dictionary["abxcc"]["-"]  = "-";
    dictionary["abxcc"]["ac"] = "h";
    dictionary["abxcc"]["bc"] = "h";
    dictionary["abxcc"]["ab"] = "b";
    dictionary["abxcc"]["aa"] = "b";
    dictionary["abxcc"]["bb"] = "b";

    dictionary["ccxab"]["-"]  = "-";
    dictionary["ccxab"]["ac"] = "h";
    dictionary["ccxab"]["bc"] = "h";
    dictionary["ccxab"]["ab"] = "a";
    dictionary["ccxab"]["aa"] = "a";
    dictionary["ccxab"]["bb"] = "a";

    //
    // Translate the genotypes for this particular map type.
    //
    translate_genotypes(types, dictionary, catalog, pmap, samples, parent_ids);

    //
    // Output the results
    //
    switch(out_type) {
    case joinmap:
	write_joinmap(catalog, pmap, types, samples, parent_ids);
	break;
    case rqtl:
	write_rqtl(catalog, pmap, types, samples, parent_ids);
	break;
    }

    if (sql_out)
	write_sql(catalog, pmap, parent_ids);

    return 0;
}

int export_cp_map(map<int, CLocus *> &catalog, PopMap<CLocus> *pmap, set<int> &parent_ids, map<int, string> &samples) {
    //
    // We wish to export, according to the JoinMap manual, a locus genotype file (loc-file), 
    // which contains the information of all the loci for a single segregating population.
    // 
    // We are exporting a CP population type:
    //   a population resulting from a cross between two heterogeneously 
    //   heterozygous and homozygous diploid parents, linkage phases originally 
    //   (possibly) unknown.
    //
    // Segregation type codes for population type CP, from Joinmap manual:
    //
    // Code      Description
    // -------   -----------
    // <abxcd>   locus heterozygous in both parents, four alleles
    // <efxeg>   locus heterozygous in both parents, three alleles
    // <hkxhk>   locus heterozygous in both parents, two alleles
    // <lmxll>   locus heterozygous in the first parent
    // <nnxnp>   locus heterozygous in the second parent
    //
    // Genotype codes for a CP population, depending on the locus segregation type.
    //
    // Seg. type   Possible genotypes
    // ---------   ------------------
    // <abxcd>     ac, ad, bc, bd, ––
    // <efxeg>     ee, ef, eg, fg, ––
    // <hkxhk>     hh, hk, kk, h-, k-, ––
    // <lmxll>     ll, lm, –– 
    // <nnxnp>     nn, np, ––
    //
    map<string, string> types;
    types["ab/--"] = "lmx--";
    types["--/ab"] = "--xnp";
    types["ab/aa"] = "lmxll";
    types["aa/ab"] = "nnxnp";
    types["ab/ab"] = "hkxhk";
    types["ab/ac"] = "efxeg";
    types["ab/cd"] = "abxcd";

    map<string, map<string, string> > dictionary;
    dictionary["lmx--"]["-"]  = "--";
    dictionary["lmx--"]["aa"] = "ll";
    dictionary["lmx--"]["bb"] = "lm";

    dictionary["--xnp"]["-"]  = "--";
    dictionary["--xnp"]["aa"] = "nn";
    dictionary["--xnp"]["bb"] = "np";

    dictionary["lmxll"]["-"]  = "--";
    dictionary["lmxll"]["aa"] = "ll";
    dictionary["lmxll"]["ab"] = "lm";

    dictionary["nnxnp"]["-"]  = "--";
    dictionary["nnxnp"]["aa"] = "nn";
    dictionary["nnxnp"]["ab"] = "np";

    dictionary["hkxhk"]["-"]  = "--";
    dictionary["hkxhk"]["ab"] = "hk";
    dictionary["hkxhk"]["aa"] = "hh";
    dictionary["hkxhk"]["bb"] = "kk";

    dictionary["efxeg"]["-"]  = "--";
    dictionary["efxeg"]["ab"] = "ef";
    dictionary["efxeg"]["ac"] = "eg";
    dictionary["efxeg"]["bc"] = "fg";
    dictionary["efxeg"]["aa"] = "ee";

    dictionary["abxcd"]["-"]  = "--";
    dictionary["abxcd"]["ac"] = "ac";
    dictionary["abxcd"]["ad"] = "ad";
    dictionary["abxcd"]["bc"] = "bc";
    dictionary["abxcd"]["bd"] = "bd";

    //
    // Translate the genotypes for this particular map type.
    //
    translate_genotypes(types, dictionary, catalog, pmap, samples, parent_ids);

    //
    // Output the results
    //
    switch(out_type) {
    case joinmap:
	write_joinmap(catalog, pmap, types, samples, parent_ids);
	break;
    case rqtl:
	write_rqtl(catalog, pmap, types, samples, parent_ids);
	break;
    }

    if (sql_out)
	write_sql(catalog, pmap, parent_ids);

    return 0;
}

int translate_genotypes(map<string, string> &types, map<string, map<string, string> > &dictionary, 
			map<int, CLocus *> &catalog, PopMap<CLocus> *pmap, map<int, string> &samples,
			set<int> &parent_ids) {
    map<int, CLocus *>::iterator it;
    CLocus *loc;

    for (it = catalog.begin(); it != catalog.end(); it++) {
	loc = it->second;

	string marker = types.count(loc->marker) ? types[loc->marker] : "";
	Datum **d = pmap->locus(loc->id);

	for (int i = 0; i < pmap->sample_cnt(); i++) {
	    if (d[i] == NULL) continue;

	    if (parent_ids.count(pmap->rev_sample_index(i))) continue;

	    //cerr << "Examining progeny " << samples[pmap->rev_sample_index(i)] << "; marker: " << loc->marker << "\n";

	    string m;

	    if (marker.length() == 0) {
		m = dictionary[marker]["-"];
	    } else {
		m = dictionary[marker].count(d[i]->gtype) ? 
		    dictionary[marker][d[i]->gtype] : 
		    dictionary[marker]["-"];
	    }
	    d[i]->trans_gtype = new char[m.length() + 1];

	    //
	    // If the genotype was corrected, output it in uppercase letteres.
	    //
	    if (d[i]->corrected) {
		for (uint k = 0; k < m.length(); k++)
		    d[i]->trans_gtype[k] = toupper(m[k]);
		d[i]->trans_gtype[m.length()] = '\0';
	    } else {
		strcpy(d[i]->trans_gtype, m.c_str());
	    }
	    if (m != dictionary[marker]["-"])
		loc->trans_gcnt++;
	    //cerr << "  allele: " << d[i]->trans_gtype << "; trans_gcnt: " << loc->trans_gcnt << "\n";
	}
    }

    return 0;
}

// sub apply_corrected_genotypes {
//     my ($sth, $loci) = @_;

//     my (%corrections, $locus, $row, $sample);

//     print STDERR "Applying manually corrected genotypes to export data...\n";

//     //
//     // Fetch the manual corrections from the database.
//     //
//     $sth->{'corr'}->execute($batch_id)
// 	or die("Unable to select results from $db.\n");

//     while ($row = $sth->{'corr'}->fetchrow_hashref()) {
//         if (!defined($corrections{$row->{'catalog_id'}})) {
//             $corrections{$row->{'catalog_id'}} = {};
//         }
//         $corrections{$row->{'catalog_id'}}->{$row->{'file'}} = $row->{'genotype'};
//     }

//     foreach $locus (@{$loci}) {
//         next if (!defined($corrections{$locus->{'id'}}));

//         foreach $sample (keys %{$corrections{$locus->{'id'}}}) {
//             $locus->{'progeny'}->{$sample} = $corrections{$locus->{'id'}}->{$sample};
//         }
//     }
// }

int tally_progeny_haplotypes(CLocus *locus, PopMap<CLocus> *pmap, set<int> &parent_ids, 
			     int &total, double &max, string &freq_str) {

    map<string, double> freq;
    Datum **d = pmap->locus(locus->id);

    total = 0;
    max   = 0;

    //cerr << "Examining marker: " << locus->id << "\n";

    for (int i = 0; i < pmap->sample_cnt(); i++) {
	if (parent_ids.count(pmap->rev_sample_index(i))) continue;
	if (d[i] == NULL) continue;

	//cerr << "  Sample: " << i << "; Haplotype: " << d[i]->obshap[0] << "; Genotype: " << d[i]->gtype << "\n";
	if (d[i]->gtype[0] != '-') {
	    freq[d[i]->gtype]++;
	    total++;
	}
    }

    if (total == 0)
	return 0;

    double frac;
    stringstream s;
    char   f[id_len];
    map<string, double>::iterator it;
    for (it = freq.begin(); it != freq.end(); it++) {
	frac = (double) it->second / (double) total * 100;
	if (frac > max) max = frac;
	sprintf(f, "(%0.1f%%);", frac);
	s << it->first << ":" << it->second << f;
    }

    freq_str = s.str();

    return 0;
}

int write_sql(map<int, CLocus *> &catalog, PopMap<CLocus> *pmap, set<int> &parent_ids) {

    if (map_type == none)
	return 0;

    stringstream pop_name;
    pop_name << "batch_" << batch_id << ".markers.tsv";
    string file = in_path + pop_name.str();

    cerr << "Writing SQL markers file to '" << file << "'\n";

    ofstream fh(file.c_str(), ofstream::out);

    if (fh.fail()) {
        cerr << "Error opening markers SQL file '" << file << "'\n";
	exit(1);
    }

    map<int, CLocus *>::iterator it;
    CLocus *loc;
    char    f[id_len], g[id_len];

    for (it = catalog.begin(); it != catalog.end(); it++) {
	loc = it->second;

	if (loc->marker.length() == 0) continue;

	string freq  = "";
	double max   = 0.0;
	int    total = 0;
	tally_progeny_haplotypes(loc, pmap, parent_ids, total, max, freq);

	sprintf(f, "%0.1f", max);
	sprintf(g, "%0.2f", loc->f);

	fh << 0 << "\t" 
	   << batch_id << "\t" 
	   << loc->id << "\t" 
	   << loc->marker << "\t"
	   << total << "\t"
	   << f << "\t"
	   << freq << "\t"
	   << g << "\n";
    }

    fh.close();

    pop_name.str("");
    pop_name << "batch_" << batch_id << ".genotypes_" << progeny_limit << ".txt";
    file = in_path + pop_name.str();

    cerr << "Writing SQL genotypes file to '" << file << "'\n";

    fh.open(file.c_str(), ofstream::out);

    if (fh.fail()) {
        cerr << "Error opening genotypes SQL file '" << file << "'\n";
	exit(1);
    }

    for (it = catalog.begin(); it != catalog.end(); it++) {
	loc = it->second;

	if (map_type == gen && loc->gcnt < progeny_limit) 
	    continue;
	else if (map_type != gen && loc->trans_gcnt < progeny_limit) 
	    continue;

	Datum **d = pmap->locus(loc->id);

	for (int i = 0; i < pmap->sample_cnt(); i++) {
	    if (parent_ids.count(pmap->rev_sample_index(i))) continue;

	    fh << 0 << "\t"
	       << batch_id << "\t"
	       << loc->id << "\t"
	       << pmap->rev_sample_index(i) << "\t";

	    if (d[i] == NULL) 
		map_type == cp ? fh << "--\n" : fh << "-\n";
	    else if (d[i]->trans_gtype != NULL)
		fh << d[i]->trans_gtype << "\n";
	    else
		fh << d[i]->gtype << "\n";
	}
    }

    fh.close();

    return 0;
}

int write_generic(map<int, CLocus *> &catalog, PopMap<CLocus> *pmap, map<int, string> &samples, set<int> &parent_ids, bool write_gtypes) {

    stringstream pop_name;
    pop_name << "batch_" << batch_id;
    if (write_gtypes)
	pop_name << ".genotypes_" << progeny_limit << ".tsv";
    else 
	pop_name << ".haplotypes_" << progeny_limit << ".tsv";

    string file = in_path + pop_name.str();

    ofstream fh(file.c_str(), ofstream::out);

    if (fh.fail()) {
        cerr << "Error opening generic output file '" << file << "'\n";
	exit(1);
    }

    //
    // Count the number of markers that have enough samples to output.
    //
    map<int, CLocus *>::iterator it;
    CLocus *loc;
    int num_loci = 0;

    for (it = catalog.begin(); it != catalog.end(); it++) {
	loc = it->second;
	if (write_gtypes == false && loc->hcnt < progeny_limit) continue;
	if (write_gtypes == true  && loc->gcnt < progeny_limit) continue;

	num_loci++;
    }
    cerr << "Writing " << num_loci << " loci to " << (write_gtypes ? "genotype" : "observed haplotype") << " file, '" << file << "'\n";

    //
    // Write the header
    //
    fh << "Catalog ID\t";
    if (expand_id)
	fh << "\t";
    if (write_gtypes)
	fh << "Marker\t";
    fh << "Cnt\t"
       << "F\t";

    map<int, string>::iterator s;
    for (int i = 0; i < pmap->sample_cnt(); i++) {
	if (write_gtypes && parent_ids.count(pmap->rev_sample_index(i))) 
	    continue;
	fh << samples[pmap->rev_sample_index(i)];
	if (i < pmap->sample_cnt() - 1) 
	    fh << "\t";
    }
    fh << "\n";

    //
    // Output each locus.
    //
    for (it = catalog.begin(); it != catalog.end(); it++) {
	loc = it->second;

	if (write_gtypes == false && loc->hcnt < progeny_limit) continue;
	if (write_gtypes == true  && loc->gcnt < progeny_limit) continue;

	stringstream id;
        loc->annotation.length() > 0 ? 
            id << loc->id << "|" << loc->annotation : id << loc->id;

	fh << id.str();

        if (expand_id) {
            if (loc->annotation.length() > 0)
                id << "\t" << loc->id << "\t" << loc->annotation;
	    else if (strlen(loc->loc.chr) > 0)
		id << "\t" << loc->id << "\t" << loc->loc.chr << "_" << loc->loc.bp;
	    else
                id << "\t" << loc->id << "\t";
        }

	if (write_gtypes)
	    fh << "\t" << loc->marker;

	write_gtypes ? fh << "\t" << loc->gcnt : fh << "\t" << loc->hcnt;
	fh << "\t" << loc->f;

	Datum **d = pmap->locus(loc->id);
	string  obshap;

	for (int i = 0; i < pmap->sample_cnt(); i++) {
	    if (write_gtypes && parent_ids.count(pmap->rev_sample_index(i))) 
		continue;
	    fh << "\t";

	    if (d[i] == NULL)
		fh << "-";
	    else
		if (write_gtypes) {
		    fh << d[i]->gtype;
		} else {
		    obshap = "";
		    for (uint j = 0; j < d[i]->obshap.size(); j++)
			obshap += string(d[i]->obshap[j]) + "/";
		    obshap = obshap.substr(0, obshap.length()-1);
		    fh << obshap;
		}
	}

	fh << "\n";
    }

    fh.close();

    return 0;
}

int write_joinmap(map<int, CLocus *> &catalog, PopMap<CLocus> *pmap, map<string, string> &types, map<int, string> &samples, set<int> &parent_ids) {

    stringstream pop_name;
    pop_name << "batch_" << batch_id << ".genotypes_" << progeny_limit;
    string file = in_path + pop_name.str() + ".loc";

    ofstream fh(file.c_str(), ofstream::out);

    if (fh.fail()) {
        cerr << "Error opening joinmap output file '" << file << "'\n";
	exit(1);
    }

    //
    // Count the number of mappable progeny
    //
    map<int, CLocus *>::iterator it;
    CLocus *loc;
    int num_loci = 0;

    for (it = catalog.begin(); it != catalog.end(); it++) {
	loc = it->second;
	if (loc->trans_gcnt < progeny_limit) continue;

	num_loci++;
    }
    cerr << "Writing " << num_loci << " loci to JoinMap file, '" << file << "'\n";

    map<int, string> map_types;
    map_types[cp]  = "CP";
    map_types[dh]  = "DH";
    map_types[bc1] = "BC1";
    map_types[f2]  = "F2";

    //
    // Output the header of the file
    //
    fh << "name = " << pop_name.str() << "\n"
       << "popt = " << map_types[map_type] << "\n"
       << "nloc = " << num_loci << "\n"
       << "nind = " << pmap->sample_cnt() - parent_ids.size() << "\n\n";

    //
    // Output each locus.
    //
    for (it = catalog.begin(); it != catalog.end(); it++) {
	loc = it->second;

	if (loc->trans_gcnt < progeny_limit) continue;

	stringstream id;
        loc->annotation.length() > 0 ? 
            id << loc->id << "|" << loc->annotation : id << loc->id;

	fh << id.str() << "\t";

        if (expand_id) {
	    id.str("");

            if (loc->annotation.length() > 0)
                id << loc->id << "\t" << loc->annotation;
	    else if (strlen(loc->loc.chr) > 0)
		id << loc->id << "\t" << loc->loc.chr << "_" << loc->loc.bp;
	    else
                id << loc->id << "\t";

            fh << id.str() << "\t";
        }

	if (types[loc->marker] == "lmx--")
	    fh << "<lmxll>";
	else if (types[loc->marker] == "--xnp")
	    fh << "<nnxnp>";
	else
	    fh << "<" << types[loc->marker] << ">";

	Datum **d = pmap->locus(loc->id);

	for (int i = 0; i < pmap->sample_cnt(); i++) {
	    if (parent_ids.count(pmap->rev_sample_index(i))) continue;
	    fh << "\t";

	    if (d[i] == NULL) 
		map_type == cp ? fh << "--" : fh << "-";
	    else
		fh << d[i]->trans_gtype;
	}

	fh << "\n";
    }

    fh << "\nindividual names:\n";

    map<int, string>::iterator s;
    for (s = samples.begin(); s != samples.end(); s++) {
	if (parent_ids.count(s->first)) continue;
	fh << s->second << "\n";
    }

    fh.close();

    return 0;
}

int write_rqtl(map<int, CLocus *> &catalog, PopMap<CLocus> *pmap, map<string, string> &types, map<int, string> &samples, set<int> &parent_ids) {
    stringstream pop_name;
    pop_name << "batch_" << batch_id << ".genotypes_" << progeny_limit;
    string file = in_path + pop_name.str() + ".loc";

    cerr << "Writing R/QTL output file to '" << file << "'\n";

    ofstream fh(file.c_str(), ofstream::out);

    if (fh.fail()) {
        cerr << "Error opening R/QTL output file '" << file << "'\n";
	exit(1);
    }

    //
    // Count the number of mappable progeny
    //
    map<int, CLocus *>::iterator it;
    CLocus *loc;
    int num_loci = 0;

    for (it = catalog.begin(); it != catalog.end(); it++) {
	loc = it->second;
	if (loc->trans_gcnt < progeny_limit) continue;

	num_loci++;
    }
    cerr<< "found " << num_loci << " loci\n";

    map<int, string> map_types;
    map_types[cp]  = "CP";
    map_types[dh]  = "DH";
    map_types[bc1] = "BC1";
    map_types[f2]  = "F2";

     //
     // Output the header of the file, followed by the list of markers, one per column
     //
    fh << "# Exported: "   << pop_name.str() << "\n"
       << "# Map Type: "   << map_types[map_type] << "\n"
       << "# Num Loci: "   << num_loci << "\n"
       << "# Num Samples " << pmap->sample_cnt() - parent_ids.size() << "\n";

    for (it = catalog.begin(); it != catalog.end(); it++) {
	loc = it->second;
	if (loc->gcnt < progeny_limit) continue;

	fh << ",";

	stringstream id;
        loc->annotation.length() > 0 ? 
            id << loc->id << "|" << loc->annotation : id << loc->id;
	fh << id.str();
    }
    fh << "\n";

    //
    // Output the chromosome (if available) for each marker and then the location
    //
    for (it = catalog.begin(); it != catalog.end(); it++) {
	loc = it->second;
	if (loc->gcnt < progeny_limit) continue;

	fh << ",";

	string chr;
        chr = strlen(loc->loc.chr) > 0 ? loc->loc.chr : "1";

	fh << chr;
    }
    fh << "\n";

    int i = 1;
    for (it = catalog.begin(); it != catalog.end(); it++) {
	loc = it->second;
	if (loc->gcnt < progeny_limit) continue;

	fh << ",";

        int bp = loc->loc.bp > 0 ? loc->loc.bp : i;

	fh << bp;
	i++;
    }
    fh << "\n";

    //
    // For each sample, print out the genotypes for each marker
    //
    Datum *d;
    for (int i = 0; i < pmap->sample_cnt(); i++) {
	if (parent_ids.count(pmap->rev_sample_index(i))) continue;
	fh << samples[pmap->rev_sample_index(i)];

	for (it = catalog.begin(); it != catalog.end(); it++) {
	    loc = it->second;
	    if (loc->gcnt < progeny_limit) continue;

	    d = pmap->datum(loc->id, i);
	    fh << ",";

	    if (d == NULL) 
		map_type == cp ? fh << "--" : fh << "-";
	    else
		fh << d->trans_gtype;
	}
	fh << "\n";
    }

    fh.close();

    return 0;
}

// sub create_imputed_genotype_map {
//     my ($order, $marker, $tag_id, $parents, $progeny, $map) = @_;

//     my (@keys, $key, $m, $type, $allele, $uall);

//     my (%gtypes, %allcnt, %uniqall);

//     //
//     // Count up the number of each type of observed allele in the progeny,
//     // record whether those genotypes are heterozygous or homozygous.
//     //
//     foreach $key (keys %{$progeny}) {
// 	my $alleles;

// 	print STDERR "Examining progeny $key\n" if ($debug);
// 	//
// 	// Discard progeny with more than one locus matched to this catalog tag.
// 	//
// 	@keys = keys %{$progeny->{$key}};
// 	next if (scalar(@keys) > 1);

// 	$alleles = join("|", sort @{$progeny->{$key}->{$keys[0]}});

// 	if (!defined($allcnt{$alleles})) {
// 	    $allcnt{$alleles} = scalar(@{$progeny->{$key}->{$keys[0]}});
// 	}
// 	//print STDERR "Adding genotype $alleles\n";

// 	$gtypes{$alleles}++;

// 	foreach $allele (@{$progeny->{$key}->{$keys[0]}}) {
// 	    $uniqall{$allele}++;
// 	}
//     }
 
//     //
//     // Examine the first parent alleles (the only alleles we have, since 
//     // we are imputing the second parent.
//     //
//     my @parents = keys %{$parents};
//     my %legal_genotypes = ();
//     $key = $order->{$parents[0]} eq "first"  ? $parents[0] : $parents[1];
//     $m   = substr($marker, 0, 2);

//     foreach $type (split(//, $m)) {
// 	//print STDERR "  Adding $type to genotypes\n" if ($debug);
//         $legal_genotypes{$type}++;
//     }
//     my @types = sort keys %legal_genotypes;

//     if ($marker eq "lmxll") {
// 	@keys = sort {$gtypes{$b} <=> $gtypes{$a}} keys %gtypes;
// 	//
// 	// Discard heterozygous alleles and find the highest frequency homozygote, 
// 	// this is the "l" in the "lmxll" marker.
// 	//
// 	while ($allcnt{$keys[0]} == 2) {
// 	    shift @keys;
// 	}
// 	$map->{$keys[0]} = shift @types;
// 	print STDERR "  Assinging '$keys[0]' to first parent genotype '", $map->{$keys[0]}, "'\n" if ($debug);

// 	foreach $uall (sort {$uniqall{$b} <=> $uniqall{$a}} keys %uniqall) {
// 	    if ($uall ne $keys[0]) {
// 		$allele = $uall;
// 		last;
// 	    }
// 	}
// 	$map->{$allele} = shift @types;
// 	print STDERR "  Assinging '$allele' to first parent genotype '", $map->{$allele}, "'\n" if ($debug);
//     }

// }

int load_marker_list(string path, set<int> &list) {
    char     line[id_len];
    ifstream fh(path.c_str(), ifstream::in);

    if (fh.fail()) {
        cerr << "Error opening white/black list file '" << path << "'\n";
	exit(1);
    }

    int   marker;
    char *e;

    while (fh.good()) {
	fh.getline(line, id_len);

	if (strlen(line) == 0) continue;

	marker = (int) strtol(line, &e, 10);

	if (*e == '\0')
	    list.insert(marker);
    }

    fh.close();

    if (list.size() == 0) {
 	cerr << "Unable to load any markers from '" << path << "'\n";
	help();
    }

    return 0;
}

int build_file_list(vector<string> &files) {
    int    pos;
    string file;
    struct dirent *direntry;

    DIR *dir = opendir(in_path.c_str());

    if (dir == NULL) {
	cerr << "Unable to open directory '" << in_path << "' for reading.\n";
	exit(1);
    }

    while ((direntry = readdir(dir)) != NULL) {
	file = direntry->d_name;

	if (file == "." || file == "..")
	    continue;

	if (file.substr(0, 6) == "batch_")
	    continue;

	pos = file.find_first_of(".");
	if (file.substr(pos+1) == "tags.tsv")
	    files.push_back(file.substr(0, pos));
    }

    if (files.size() == 0) {
	cerr << "Unable to locate any input files to process within '" << in_path << "'\n";
    }

    cerr << "Found " << files.size() << " input file(s).\n";

    return 0;
}

bool hap_compare(pair<string, int> a, pair<string, int> b) {
    return (a.second > b.second);
}

bool compare_pair(pair<char, int> a, pair<char, int> b) {
    return (a.second > b.second);
}

int parse_command_line(int argc, char* argv[]) {
    int c;
     
    while (1) {
	static struct option long_options[] = {
	    {"help",        no_argument,       NULL, 'h'},
            {"version",     no_argument,       NULL, 'v'},
            {"corr",        no_argument,       NULL, 'c'},
            {"sql",         no_argument,       NULL, 's'},
	    {"num_threads", required_argument, NULL, 'p'},
	    {"batch_id",    required_argument, NULL, 'b'},
	    {"in_path",     required_argument, NULL, 'P'},
	    {"map_type",    required_argument, NULL, 't'},
	    {"out_type",    required_argument, NULL, 'o'},
	    {"progeny",     required_argument, NULL, 'r'},
	    {"min_depth",   required_argument, NULL, 'm'},
	    {"whitelist",   required_argument, NULL, 'W'},
	    {"blacklist",   required_argument, NULL, 'B'},
	    {0, 0, 0, 0}
	};
	
	// getopt_long stores the option index here.
	int option_index = 0;
     
	c = getopt_long(argc, argv, "hvcsib:p:t:o:r:P:m:W:B:", long_options, &option_index);
     
	// Detect the end of the options.
	if (c == -1)
	    break;
     
	switch (c) {
	case 'h':
	    help();
	    break;
	case 'P':
	    in_path = optarg;
	    break;
	case 'b':
	    batch_id = atoi(optarg);
	    break;
	case 't':
	    if (strcasecmp(optarg, "cp") == 0)
                map_type = cp;
	    else if (strcasecmp(optarg, "bc1") == 0)
		map_type = bc1;
	    else if (strcasecmp(optarg, "f2") == 0)
		map_type = f2;
	    else if (strcasecmp(optarg, "dh") == 0)
		map_type = dh;
	    else if (strcasecmp(optarg, "gen") == 0)
		map_type = gen;
	    else
		map_type = unk;
	    break;
	case 'o':
	    if (strcasecmp(optarg, "joinmap") == 0)
                out_type = joinmap;
	    else if (strcasecmp(optarg, "rqtl") == 0)
		out_type = rqtl;
	    break;
	case 'r':
	    progeny_limit = atoi(optarg);
	    break;
	case 'c':
	    corrections = true;
	    break;
	case 'i':
	    expand_id = true;
	    break;
	case 's':
	    sql_out = true;
	    break;
	case 'W':
	    wl_file = optarg;
	    break;
	case 'B':
	    bl_file = optarg;
	    break;
	case 'm':
	    min_stack_depth = atoi(optarg);
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

    if (in_path.length() == 0) {
	cerr << "You must specify a path to the directory containing Stacks output files.\n";
	help();
    }

    if (in_path.at(in_path.length() - 1) != '/') 
	in_path += "/";

    if (batch_id == 0) {
	cerr << "You must specify a batch ID.\n";
	help();
    }

    if (map_type != cp  &&
	map_type != dh  &&
	map_type != bc1 &&
	map_type != f2  &&
	map_type != gen &&
	map_type != none) {
        cerr << "You must specify a valid map type. 'CP', 'DH', 'F2', 'BC1' and 'GEN' are the currently supported map types.\n";
        help();
    }

    if (map_type == none && min_stack_depth > 0)
	cerr << "Warning: using a minimum stack depth when building genetic markers is not recommended.\n";

    return 0;
}

void version() {
    std::cerr << "genotypes " << VERSION << "\n\n";

    exit(0);
}

void help() {
    std::cerr << "genotypes " << VERSION << "\n"
              << "genotypes -b batch_id -P path [-r min] [-t map_type -o type] [-B blacklist] [-W whitelist] [-c] [-C] [-s] [-v] [-h]" << "\n"
	      << "  b: Batch ID to examine when exporting from the catalog.\n"
	      << "  r: minimum number of progeny required to print a marker.\n"
	      << "  c: make automated corrections to the data.\n"
	      << "  P: path to the Stacks output files.\n"
	      << "  t: map type to write. 'CP', 'DH', 'F2' and 'BC1' are the currently supported map types.\n"
	      << "  o: output file type to write, 'joinmap' and 'rqtl' are currently supported.\n"
	      << "  C: apply manual corrections (that were made via the web interface) to the data.\n"
	      << "  s: output a file to import results into an SQL database.\n"
	      << "  B: specify a file containing Blacklisted markers to be excluded from the export.\n"
	      << "  W: specify a file containign Whitelisted markers to include in the export.\n"
	      << "  v: print program version." << "\n"
	      << "  h: display this help messsage." << "\n\n";

    exit(0);
}
