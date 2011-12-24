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
// populations -- generate population genetic statistics and output 
// haplotypes in a population context.
//
// Julian Catchen
// jcatchen@uoregon.edu
// University of Oregon
//
// $Id: populations.cc 2117 2011-05-07 00:54:22Z catchen $
//

#include "populations.h"

// Global variables to hold command-line options.
int       num_threads = 1;
int       batch_id    = 0;
string    in_path;
string    out_path;
string    out_file;
string    pmap_path;
string    bl_file;
string    wl_file;
string    enz;
int       progeny_limit   = 1;
bool      corrections     = false;
bool      expand_id       = false;
bool      sql_out         = false;
int       min_stack_depth = 0;

map<int, pair<int, int> > pop_indexes;
set<int> whitelist, blacklist;

//
// Hold information about restriction enzymes
//
map<string, const char **> renz;
map<string, int>           renz_cnt;
map<string, int>           renz_len;

int main (int argc, char* argv[]) {

    initialize_renz(renz, renz_cnt, renz_len);

    parse_command_line(argc, argv);

    //
    // Set the number of OpenMP parallel threads to execute.
    //
    #ifdef _OPENMP
    omp_set_num_threads(num_threads);
    #endif

    vector<pair<int, string> > files;
    build_file_list(files, pop_indexes);

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
	load_catalog_matches(in_path + files[i].second, m);

	if (m.size() == 0) {
	    cerr << "Warning: unable to find any matches in file '" << files[i].second << "', excluding this sample from genotypes analysis.\n";
	    continue;
	}

	catalog_matches.push_back(m);
 	samples[m[0]->sample_id] = files[i].second;
	sample_ids.push_back(m[0]->sample_id);
    }

    //
    // Create the population map
    //
    cerr << "Populating observed haplotypes for " << sample_ids.size() << " samples, " << catalog.size() << " loci.\n";
    PopMap<CLocus> *pmap = new PopMap<CLocus>(sample_ids.size(), catalog.size());
    pmap->populate(sample_ids, catalog, catalog_matches, min_stack_depth);

    cerr << "Loading model outputs for " << sample_ids.size() << " samples, " << catalog.size() << " loci.\n";
    map<int, CLocus *>::iterator it;
    map<int, ModRes *>::iterator mit;
    Datum  *d;
    CLocus *loc;
    //
    // Load the output from the SNP calling model for each individual at each locus. This
    // model output string looks like this:
    //   OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOEOOOOOOEOOOOOOOOOOOOOOOOOOOOOOOOOOOOOUOOOOUOOOOOO
    // and records model calls for each nucleotide: O (hOmozygous), E (hEterozygous), U (Unknown)
    //
    for (uint i = 0; i < files.size(); i++) {
	map<int, ModRes *> modres;
	load_model_results(in_path + files[i].second, modres);

	if (modres.size() == 0) {
	    cerr << "Warning: unable to find any model results in file '" << files[i].second << "', excluding this sample from population analysis.\n";
	    continue;
	}

	for (it = catalog.begin(); it != catalog.end(); it++) {
	    loc = it->second;
	    d = pmap->datum(loc->id, sample_ids[i]);

	    if (d != NULL) {
		d->model = new char[strlen(modres[d->id]->model) + 1];
		strcpy(d->model, modres[d->id]->model);
	    }
	}

	for (mit = modres.begin(); mit != modres.end(); mit++)
	    delete mit->second;
	modres.clear();
    }

    uint pop_id, start_index, end_index;
    map<int, pair<int, int> >::iterator pit;

    PopSum<CLocus> *psum = new PopSum<CLocus>(pmap->loci_cnt(), pop_indexes.size());
    psum->initialize(pmap);

    for (pit = pop_indexes.begin(); pit != pop_indexes.end(); pit++) {
	start_index = pit->second.first;
	end_index   = pit->second.second;
	pop_id      = pit->first;
	cerr << "Generating nucleotide-level summary statistics for population " << pop_id << "\n";
	psum->add_population(catalog, pmap, pop_id, start_index, end_index);
    }

    //
    // Output the locus-level summary statistics.
    //
    write_summary_stats(files, pop_indexes, catalog, psum);

    //
    // Output the observed haplotypes.
    //
    write_generic(catalog, pmap, samples, parent_ids, false);

    //
    // Output nucleotide-level genotype calls for each individual.
    //
    write_genomic(catalog, pmap);

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

	//cerr << "Sample Id: " << pmap->rev_sample_index(i) << "\n";

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
	    //cerr << "  Observed Haplotype: " << d[i]->obshap[j] << ", Genotype: " << locus->gmap[d[i]->obshap[j]] << "\n";
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

int write_genomic(map<int, CLocus *> &catalog, PopMap<CLocus> *pmap) {
    stringstream pop_name;
    pop_name << "batch_" << batch_id << ".genomic_" << progeny_limit << ".tsv";

    string file = in_path + pop_name.str();

    ofstream fh(file.c_str(), ofstream::out);

    if (fh.fail()) {
        cerr << "Error opening genomic output file '" << file << "'\n";
	exit(1);
    }

    //
    // Count the number of markers that have enough samples to output.
    //
    map<int, CLocus *>::iterator cit;
    CLocus *loc;
    int num_loci = 0;

    for (cit = catalog.begin(); cit != catalog.end(); cit++) {
	loc = cit->second;
	if (loc->hcnt < progeny_limit) continue;

	num_loci += loc->len - renz_len[enz];
    }
    cerr << "Writing " << num_loci << " nucleotide positions to genomic file, '" << file << "'\n";

    //
    // Write the header
    //
    fh << num_loci << "\t" << pmap->sample_cnt() << "\n";

    //
    // Output each locus.
    //
    map<string, vector<CLocus *> >::iterator it;
    int  a, b;

    uint  rcnt = renz_cnt[enz];
    uint  rlen = renz_len[enz];
    char *p;

    for (it = pmap->ordered_loci.begin(); it != pmap->ordered_loci.end(); it++) {
	for (uint i = 0; i < it->second.size(); i++) {
	    loc = it->second[i];

	    if (loc->hcnt < progeny_limit) continue;

	    Datum **d = pmap->locus(loc->id);
	    set<int> snp_locs;
	    string   obshap;

	    for (uint i = 0; i < loc->snps.size(); i++)
		snp_locs.insert(loc->snps[i]->col);

	    uint start = 0;
	    uint end   = loc->len;
	    //
	    // Check for the existence of the restriction enzyme cut site, mask off 
	    // its output.
	    //
	    for (uint n = 0; n < rcnt; n++)
	    	if (strncmp(loc->con, renz[enz][n], rlen) == 0)
	    	    start += renz_len[enz];
	    if (start == 0) {
	    	p = loc->con + (loc->len - rlen);
	    	for (uint n = rcnt; n < rcnt + rcnt; n++)
	    	    if (strncmp(p, renz[enz][n], rlen) == 0)
	    		end -= renz_len[enz];
	    }

	    uint k = 0;
	    for (uint n = start; n < end; n++) {
		fh << loc->id << "\t" << loc->loc.chr << "\t" << loc->loc.bp + n;

 		if (snp_locs.count(n) == 0) {
		    for (int j = 0; j < pmap->sample_cnt(); j++) {
			a = encode_gtype(loc->con[n]);
			fh << "\t" << encoded_gtypes[a][a];
		    }
		} else {
		    for (int j = 0; j < pmap->sample_cnt(); j++) {
			fh << "\t";

			if (d[j] == NULL)
			    fh << "0";
			else 
			    switch (d[j]->obshap.size()) {
			    case 1:
				a = encode_gtype(d[j]->obshap[0][k]);
				fh << encoded_gtypes[a][a];
				break;
			    case 2:
				a = encode_gtype(d[j]->obshap[0][k]);
				b = encode_gtype(d[j]->obshap[1][k]);
				fh << encoded_gtypes[a][b];
				break;
			    default:
				fh << "0";
				break;
			    }
		    }
		    k++;
		}
		fh << "\n";
	    }
	}
    }

    fh.close();

    return 0;
}

int 
write_summary_stats(vector<pair<int, string> > &files, map<int, pair<int, int> > &pop_indexes, 
		    map<int, CLocus *> &catalog, PopSum<CLocus> *psum) 
{
    stringstream pop_name;
    pop_name << "batch_" << batch_id << ".sumstats_" << progeny_limit << ".tsv";

    string file = in_path + pop_name.str();

    ofstream fh(file.c_str(), ofstream::out);

    if (fh.fail()) {
        cerr << "Error opening generic output file '" << file << "'\n";
	exit(1);
    }

    int start, end;
    //
    // Write the population members.
    //
    map<int, pair<int, int> >::iterator pit;
    for (pit = pop_indexes.begin(); pit != pop_indexes.end(); pit++) {
	start = pit->second.first;
	end   = pit->second.second;
	fh << pit->first << "\t";
	for (int i = start; i <= end; i++) {
	    fh << files[i].second;
	    if (i < end) fh << ",";
	}
	fh << "\n";
    }

    cerr << "Writing " << catalog.size() << " loci to summary statistics file, '" << file << "'\n";

    map<int, CLocus *>::iterator it;
    CLocus  *loc;
    LocSum **s;
    int      len;
    int      pop_cnt = psum->pop_cnt();

    fh << "Locus ID" << "\t"
       << "Chr" << "\t"
       << "BP" << "\t"
       << "Pop ID" << "\t"
       << "N" << "\t"
       << "P" << "\t"
       << "Obs Het" << "\t"
       << "Obs Hom" << "\t"
       << "Exp Het" << "\t"
       << "Exp Hom" << "\n";

    for (it = catalog.begin(); it != catalog.end(); it++) {
	loc = it->second;

	s = psum->locus(loc->id);
	len = strlen(loc->con);

	for (int i = 0; i < len; i++) {

	    for (int j = 0; j < pop_cnt; j++) {

		if (s[j]->nucs[i].num_indv < progeny_limit) continue;

		fh << loc->id << "\t"
		   << loc->loc.chr << "\t"
		   << loc->loc.bp + i << "\t"
		   << psum->rev_pop_index(j) << "\t"
		   << s[j]->nucs[i].num_indv << "\t"
		   << s[j]->nucs[i].p << "\t"
		   << s[j]->nucs[i].obs_het << "\t"
		   << s[j]->nucs[i].obs_hom << "\t"
		   << s[j]->nucs[i].exp_het << "\t"
		   << s[j]->nucs[i].exp_hom << "\n";
	    }
	}
    }

    return 0;
}

int
write_generic(map<int, CLocus *> &catalog, PopMap<CLocus> *pmap, 
	      map<int, string> &samples, set<int> &parent_ids, bool write_gtypes)
{
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

int build_file_list(vector<pair<int, string> > &files, map<int, pair<int, int> > &pop_indexes) {
    char   line[max_len];
    vector<string> parts;

    ifstream fh(pmap_path.c_str(), ifstream::in);

    if (fh.fail()) {
        cerr << "Error opening population map '" << pmap_path << "'\n";
	exit(1);
    }

    while (fh.good()) {
	fh.getline(line, max_len);

	if (strlen(line) == 0) continue;

	//
	// Parse the population map, we expect:
	// <file name> <tab> <population ID>
	//
	parse_tsv(line, parts);

	files.push_back(make_pair(atoi(parts[1].c_str()), parts[0]));
    }

    fh.close();

    //
    // Sort the files according to population ID.
    //
    sort(files.begin(), files.end(), compare_pop_map);

    if (files.size() == 0) {
	cerr << "Unable to locate any input files to process within '" << in_path << "'\n";
    }

    cerr << "Found " << files.size() << " input file(s).\n";

    //
    // Determine the start/end index for each population in the files array.
    //
    int start  = 0;
    int end    = 0;
    int pop_id = files[0].first;

    do {
	end++;
	if (pop_id != files[end].first) {
	    pop_indexes[pop_id] = make_pair(start, end - 1);
	    start  = end;
	    pop_id = files[end].first;
	}
    } while (end < (int) files.size());

    cerr << "  " << pop_indexes.size() << " populations found\n";

    map<int, pair<int, int> >::iterator it;
    for (it = pop_indexes.begin(); it != pop_indexes.end(); it++) {
	start = it->second.first;
	end   = it->second.second;
	cerr << "    " << it->first << ": ";
	for (int i = start; i <= end; i++) {
	    cerr << files[i].second;
	    if (i < end) cerr << ", ";
	}
	cerr << "\n";
    }

    return 0;
}

bool compare_pop_map(pair<int, string> a, pair<int, string> b) {
    return (a.first < b.first);
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
	    {"progeny",     required_argument, NULL, 'r'},
	    {"min_depth",   required_argument, NULL, 'm'},
	    {"renz",        required_argument, NULL, 'e'},
	    {"pop_map",     required_argument, NULL, 'M'},
	    {"whitelist",   required_argument, NULL, 'W'},
	    {"blacklist",   required_argument, NULL, 'B'},
	    {0, 0, 0, 0}
	};
	
	// getopt_long stores the option index here.
	int option_index = 0;
     
	c = getopt_long(argc, argv, "hvcsib:p:t:o:r:M:P:m:e:W:B:", long_options, &option_index);
     
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
	case 'M':
	    pmap_path = optarg;
	    break;
	case 'b':
	    batch_id = atoi(optarg);
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
	case 'e':
	    enz = optarg;
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

    if (pmap_path.length() == 0) {
	cerr << "You must specify a path to the population map.\n";
	help();
    }

    if (batch_id == 0) {
	cerr << "You must specify a batch ID.\n";
	help();
    }

    if (enz.length() == 0) {
	cerr << "You must specify the restriction enzyme used with 'genomic' output.\n";
	help();
    }

    if (renz.count(enz) == 0) {
	cerr << "Unrecognized restriction enzyme specified: '" << enz.c_str() << "'.\n";
	help();
    }

    return 0;
}

void version() {
    std::cerr << "populations " << VERSION << "\n\n";

    exit(0);
}

void help() {
    std::cerr << "populations " << VERSION << "\n"
              << "populations -b batch_id -P path -M path [-r min] [-m min] [-t map_type -o type] [-B blacklist] [-W whitelist] [-c] [-s] [-e renz] [-v] [-h]" << "\n"
	      << "  b: Batch ID to examine when exporting from the catalog.\n"
	      << "  P: path to the Stacks output files.\n"
	      << "  M: path to the population map, a tab separated file describing which individuals belong in which population.\n"
	      << "  r: minimum number of individuals required to process a locus.\n"
	      << "  m: specify a minimum stack depth required before exporting a locus in a particular individual.\n"
	      << "  s: output a file to import results into an SQL database.\n"
	      << "  B: specify a file containing Blacklisted markers to be excluded from the export.\n"
	      << "  W: specify a file containign Whitelisted markers to include in the export.\n"
	      << "  e: restriction enzyme, required if generating 'genomic' output.\n"
	      << "  v: print program version." << "\n"
	      << "  h: display this help messsage." << "\n\n";

    exit(0);
}
