// -*-mode:c++; c-style:k&r; c-basic-offset:4;-*-
//
// Copyright 2013, Julian Catchen <jcatchen@uoregon.edu>
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
// rxstacks -- make model call corrections and haplotype corrections
// across a population of samples.
//
// Julian Catchen
// jcatchen@uoregon.edu
// University of Oregon
//

#include "rxstacks.h"

// Global variables to hold command-line options.
int    num_threads = 1;
int    batch_id    = 0;
string in_path;
string out_path;
double confounded_limit  = 0.75;
bool   filter_confounded = false;
bool   prune_haplotypes  = false;
int    max_haplotype_cnt = 0;
bool   lnl_dist          = false;
bool   filter_lnl        = false;
double lnl_limit         = 0.0;
bool   verbose           = false;

//
// For use with the multinomial model to call fixed nucleotides.
//
modelt model_type         = snp;
double alpha              = 0.1;
double bound_low          = 0.0;
double bound_high         = 1.0;
double p_freq             = 0.5;
double barcode_err_freq   = 0.0;
double heterozygote_limit = -2.71;
double homozygote_limit   =  2.71;
const int barcode_size    = 5;

int main (int argc, char* argv[]) {

    parse_command_line(argc, argv);

    cerr
	<< "Log liklihood filtering: " << (filter_lnl        == true ? "on" : "off") << "; threshold: " << lnl_limit << "\n"
	<< "Prune haplotypes: "        << (prune_haplotypes  == true ? "on" : "off") << "\n"
	<< "Filter confounded loci: "  << (filter_confounded == true ? "on" : "off") << "\n";

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
    if (verbose) num_threads = 1;
    #ifdef _OPENMP
    omp_set_num_threads(num_threads);
    #endif

    vector<pair<int, string> > files;
    if (!build_file_list(files))
	exit(1);

    //
    // Open and initialize the log file.
    //
    ofstream log_fh;
    init_log(argc, argv, log_fh);

    //
    // Load the catalog
    //
    stringstream catalog_file;
    map<int, CSLocus *> catalog;
    int res;
    catalog_file << in_path << "batch_" << batch_id << ".catalog";
    if ((res = load_loci(catalog_file.str(), catalog, false, false)) == 0) {
    	cerr << "Unable to load the catalog '" << catalog_file.str() << "'\n";
     	return 0;
    }

    //
    // Let's fill in the SNP model calls to include both hets and homozygotes to 
    // make it easier to iterate over them later.
    //
    fill_catalog_snps(catalog);

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
	    cerr << "Warning: unable to find any matches in file '" << files[i].second << "', excluding this sample from population analysis.\n";
	    continue;
	}

	catalog_matches.push_back(m);
	if (samples.count(m[0]->sample_id) == 0) {
	    samples[m[0]->sample_id] = files[i].second;
	    sample_ids.push_back(m[0]->sample_id);
	} else {
	    cerr << "Fatal error: sample ID " << m[0]->sample_id << " occurs twice in this data set, likely the pipeline was run incorrectly.\n";
	    exit(0);
	}
    }

    //
    // Create the population map
    // 
    cerr << "Populating observed haplotypes for " << sample_ids.size() << " samples, " << catalog.size() << " loci.\n";
    PopMap<CSLocus> *pmap = new PopMap<CSLocus>(sample_ids.size(), catalog.size());
    pmap->populate(sample_ids, catalog, catalog_matches);

    //
    // Sum haplotype counts across the population for each catalog locus.
    //
    sum_haplotype_counts(catalog, pmap);

    //
    // Calculate mean log likelihood across the population for each catalog locus.
    //
    calc_lnl_means(catalog, pmap);

    int    catalog_id, sample_id, tag_id;
    string file;

    //
    // Process samples matched to the catalog, one by one.
    //
    for (uint i = 0; i < catalog_matches.size(); i++) {
	sample_id = catalog_matches[i][0]->sample_id;
	file      = samples[sample_id];

	cerr << "Loading stacks from sample " << file << " [" << i+1 << " of " << catalog_matches.size() << "]...\n";

	//// if (sample_id != 160 && sample_id != 147) continue;

	map<int, Locus *> stacks;
	int res;
	if ((res = load_loci(in_path + file, stacks, true, true)) == 0) {
	    cerr << "Unable to load sample file '" << file << "'\n";
	    continue;
	}

	cerr << "Making corrections to sample " << file << "...";
	if (verbose)
	    log_fh << "\n# Sample " << file << ", Sample ID " << sample_id << "\n"
		   << "# Sample Id\t" 
		   << "Locus ID\t" 
		   << "SNP Col\t" 
		   << "Orig Value\t" 
		   << "Corr Value\n";
	
	set<pair<int, int> >           uniq_matches;
	set<pair<int, int> >::iterator it;
	vector<pair<int, int> >        matches;

	//
	// There are multiple matches per stack, but we only need to process
	// each stack once to make corrections.
	//
	for (uint j = 0; j < catalog_matches[i].size(); j++) {
	    catalog_id = catalog_matches[i][j]->cat_id;
	    tag_id     = catalog_matches[i][j]->tag_id;

	    uniq_matches.insert(make_pair(catalog_id, tag_id));
	}

	//
	// Put the catalog/tag ID pairs into a vector for parallel processing.
	//
	for (it = uniq_matches.begin(); it != uniq_matches.end(); it++)
	    matches.push_back(*it);

	unsigned long int nuc_cnt        = 0;
	unsigned long int unk_hom_cnt    = 0;
	unsigned long int unk_het_cnt    = 0;
	unsigned long int het_unk_cnt    = 0;
	unsigned long int hom_unk_cnt    = 0;
	unsigned long int het_hom_cnt    = 0;
	unsigned long int hom_het_cnt    = 0;
	unsigned long int conf_loci_cnt  = 0;
	unsigned long int pruned_hap_cnt = 0;
	unsigned long int blacklist_cnt  = 0;
	unsigned long int lnl_cnt        = 0;

        #pragma omp parallel private(catalog_id, tag_id)
	{ 
	    Datum   *d;
	    Locus   *loc;
	    CSLocus *cloc;

            #pragma omp for schedule(dynamic, 1) reduction(+:nuc_cnt) reduction(+:unk_hom_cnt) reduction(+:unk_het_cnt) \
		reduction(+:hom_unk_cnt) reduction(+:het_unk_cnt) reduction(+:hom_het_cnt) reduction(+:het_hom_cnt) \
		reduction(+:conf_loci_cnt) reduction(+:pruned_hap_cnt) reduction(+:blacklist_cnt) reduction(+:lnl_cnt)
	    for (uint j = 0; j < matches.size(); j++) {
		catalog_id = matches[j].first;
		tag_id     = matches[j].second;

		//// if (tag_id == 8334) {
		////     cerr << "Hit the tag.\n";
		//// }

		if (catalog.count(catalog_id) == 0) continue;

		cloc = catalog[catalog_id];
		loc  = stacks[tag_id];

		if (filter_confounded && 
		    ((double) cloc->confounded_cnt / (double) cloc->cnt > confounded_limit)) {
		    // cerr << "Catalog locus " << cloc->id << " is confounded; confounded cnt: " 
		    // 	 << cloc->confounded_cnt << "; total: " << cloc->cnt 
		    // 	 << "; freq: " << (double) cloc->confounded_cnt / (double)cloc->cnt << "\n";
		    loc->blacklisted = true;
		    conf_loci_cnt++;
		}

		d = pmap->datum(catalog_id, sample_id);

		if (d == NULL) continue;

		if (filter_lnl && cloc->lnl < lnl_limit) {
		    loc->blacklisted = true;
		    lnl_cnt++;
		}

		prune_nucleotides(cloc, loc, log_fh,
				  nuc_cnt,
				  unk_hom_cnt, unk_het_cnt, 
				  hom_unk_cnt, het_unk_cnt, 
				  hom_het_cnt, het_hom_cnt);

		//
		// Prune haplotypes from this locus.
		//
		if (prune_haplotypes) {
		    prune_locus_haplotypes(cloc, d, loc, pruned_hap_cnt);

		    if (loc->blacklisted) {
			if (verbose)
			    log_fh << "Blacklisted sample " << sample_id << ", locus " << loc->id << " due to inability to call haplotypes.\n";
			blacklist_cnt++;
		    }
		}
	    }
        }

	cerr << "done.\n";

	unsigned long int total = unk_hom_cnt + unk_het_cnt + hom_unk_cnt + het_unk_cnt + hom_het_cnt + het_hom_cnt;

	cerr << conf_loci_cnt << " confounded loci were blacklisted and not processed.\n"
	     << "Total nucleotides processed: " << nuc_cnt << "\n"
	     << "  Total nucleotides converted: " << total << "\n"
	     << "    Converted from unknown to homozygous:      " << unk_hom_cnt << " nucleotides.\n"
	     << "    Converted from unknown to heterozygous:    " << unk_het_cnt << " nucleotides.\n"
	     << "    Converted from homozygous to unknown:      " << hom_unk_cnt << " nucleotides.\n"
	     << "    Converted from heterozygous to unknown:    " << het_unk_cnt << " nucleotides.\n"
	     << "    Converted from homozygous to heterozygous: " << hom_het_cnt << " nucleotides.\n"
	     << "    Converted from heterozygous to homozygous: " << het_hom_cnt << " nucleotides.\n"
	     << "    Blacklisted: " << blacklist_cnt << " loci due to inability to call haplotypes.\n"
	     << "    Blacklisted: " << lnl_cnt << " loci due to log likelihoods below threshold.\n"
	     << pruned_hap_cnt << " haplotypes were pruned.\n";

	if (verbose)
	    log_fh << "# Sample\t"
		   << "Confounded loci\t"
		   << "Total nucs\t"
		   << "Total nucs converted\t"
		   << "Unk to Hom\t"
		   << "Unk to Het\t"
		   << "Hom to Unk\t"
		   << "Het to Unk\t"
		   << "Hom to Het\t"
		   << "Het to Hom\t"
		   << "Blacklisted\t"
		   << "Lnl Blacklisted\t"
		   << "Pruned Haplotypes\n";
	log_fh << file           << "\t"
	       << conf_loci_cnt  << "\t"
	       << nuc_cnt        << "\t"
	       << total          << "\t"
	       << unk_hom_cnt    << "\t"
	       << unk_het_cnt    << "\t"
	       << hom_unk_cnt    << "\t"
	       << het_unk_cnt    << "\t"
	       << hom_het_cnt    << "\t"
	       << het_hom_cnt    << "\t"
	       << blacklist_cnt  << "\t"
	       << lnl_cnt        << "\t"
	       << pruned_hap_cnt << "\n";

	cerr << "Writing modified stacks, SNPs, alleles to '" << out_path << "'...";

	//
	// Rewrite stacks, model outputs, and haplotypes.
	//
	write_results(file, stacks);

	//
	// Free up memory
	//
	map<int, Locus *>::iterator stack_it;
	for (stack_it = stacks.begin(); stack_it != stacks.end(); stack_it++)
	    delete stack_it->second;
    }

    log_fh.close();

    return 0;
}

int
calc_lnl_means(map<int, CSLocus *> &catalog, PopMap<CSLocus> *pmap)
{
    map<int, CSLocus *>::iterator it;
    CSLocus *cloc;
    Datum  **d;
    uint     cnt, mid;
    double   median, mean;
    vector<double> lnls, means;
    ofstream log_fh;

    if (lnl_dist) {
	stringstream log;
	log << "batch_" << batch_id << ".rxstacks_lnls.tsv";
	string log_path = out_path + log.str();
	log_fh.open(log_path.c_str(), ofstream::out);

	if (log_fh.fail()) {
	    cerr << "Error opening log file '" << log_path << "'\n";
	    exit(1);
	}
	log_fh << "# Catalog Locus\tMean\tMedian\n";
    }

    for (it = catalog.begin(); it != catalog.end(); it++) {
	cloc = it->second;

	d    = pmap->locus(cloc->id);
	cnt  = pmap->sample_cnt();
	mean     = 0.0;
	lnls.clear();

	for (uint i = 0; i < cnt; i++) {
	    if (d[i] == NULL) continue;

	    lnls.push_back(d[i]->lnl);
	    mean += d[i]->lnl;
	}

	sort(lnls.begin(), lnls.end());

	mid    = lnls.size() / 2;
	median = lnls.size() % 2 == 0 ? lnls[mid] + lnls[mid+1] / 2.0 : lnls[mid+1];
	mean   = mean / (double) lnls.size(); 

	cloc->lnl = mean;

	means.push_back(mean);

	if (lnl_dist)
	    log_fh << cloc->id << "\t" 
		   << mean << "\t"
		   << median << "\n";
    }

    if (lnl_dist)
	log_fh.close();

    return 0;
}

int
sum_haplotype_counts(map<int, CSLocus *> &catalog, PopMap<CSLocus> *pmap)
{
    map<int, CSLocus *>::iterator it;
    CSLocus *cloc;
    Datum  **d;
    uint     cnt;

    for (it = catalog.begin(); it != catalog.end(); it++) {
	cloc = it->second;

	d   = pmap->locus(cloc->id);
	cnt = pmap->sample_cnt();

	for (uint i = 0; i < cnt; i++) {
	    if (d[i] == NULL) continue;

	    if (d[i]->obshap.size() == 1) {
		if (cloc->hap_cnts.count(d[i]->obshap[0]) == 0)
		    cloc->hap_cnts[d[i]->obshap[0]]  = 2;
		else
		    cloc->hap_cnts[d[i]->obshap[0]] += 2;
	    } else {
		for (uint j = 0; j < d[i]->obshap.size(); j++)
		    if (cloc->hap_cnts.count(d[i]->obshap[j]) == 0)
			cloc->hap_cnts[d[i]->obshap[j]]  = 1;
		    else
			cloc->hap_cnts[d[i]->obshap[j]] += 1;
	    }
	}
    }

    return 0;
}

int
prune_locus_haplotypes(CSLocus *cloc, Datum *d, Locus *loc, unsigned long &pruned_hap_cnt)
{
    if (d->obshap.size() < 2) return 0;

    //
    // Identify the two most frequent haplotypes in this sample.
    //
    vector<pair<string, double> > haplotypes;
    double weighted_hap;

    for (uint i = 0; i < d->obshap.size(); i++) {
	//
	// Lookup the number of occurrences of this haplotype in the 
	// population as well as the depth of the haplotype in this indiviudal. 
	// We will weight the occurrences of the haplotype in the population by the natural log
	// of the read depth of the haplotype in this individual, storing the result.
	//
	weighted_hap = (double) cloc->hap_cnts[d->obshap[i]] * log((double) d->depth[i]);
	haplotypes.push_back(make_pair(string(d->obshap[i]), weighted_hap));
    }

    //
    // Sort according to haplotype frequency.
    //
    sort(haplotypes.begin(), haplotypes.end(), compare_pair_haplotype);

    if (haplotypes.size() == 0) {
	cerr << "Error processing catalog locus " << cloc->id << "\n";
	return -1;
    }

    //
    // Prune out excess haplotypes.
    //
    map<string, int>::iterator it;
    uint   j, k, cat_snp, cat_idx, loc_snp, loc_idx;
    string hap;

    for (uint i = 2; i < haplotypes.size(); i++) {
	//
	// Make sure that those haplotypes we want to discard occur at a frequency lower
	// than the second most frequent haplotype, instead of being tied for second.
	//
	if (haplotypes[i].second >= haplotypes[1].second ||
	    (max_haplotype_cnt > 0 && haplotypes[i].second > max_haplotype_cnt))
	    continue;

	cat_snp =  0;
	cat_idx = -1;
	loc_snp =  0;
	loc_idx = -1;
	k   = -1;
	j   = -1;
	hap = "";

	do {
	    j++;
	    loc_idx++;
	    //
	    // Advance to a het in the sample locus.
	    //
	    while (j < loc->snps.size() && loc->snps[j]->type != snp_type_het) j++;
	    if (j >= loc->snps.size()) break;
	    loc_snp = loc->snps[j]->col;

	    do {
		k++;
		cat_idx++;
		//
		// Advance to the het in the catalog locus that corresponds to the sample locus.
		//
		while (k < cloc->snps.size() && cloc->snps[k]->type != snp_type_het) k++;
		if (k >= cloc->snps.size()) break;
		cat_snp = cloc->snps[k]->col;

	    } while (cat_snp < loc_snp);

	    //
	    // Extract out the nucleotide from the catalog haplotype that matches the sample 
	    // haplotype. For example, catalog haplotype may be 'ACGTG' while sample haplotype 
	    // is 'CT'.
	    //
	    if (j < loc->snps.size() && k < cloc->snps.size() && cat_snp == loc_snp) {
		hap += haplotypes[i].first.at(cat_idx);
	    } else {
		cerr << "Error processing catalog locus " << cloc->id << "\n";
		return -1;
	    }

	} while (j < loc->snps.size());

	//
	// Remove the haplotype.
	//
	it = loc->alleles.find(hap);

	if (it != loc->alleles.end()) {
	    loc->alleles.erase(it);
	    pruned_hap_cnt++;
	} 
	// else {
	//     cerr << "Error erasing haplotype '" << hap << "' in catalog locus " << cloc->id << ", sample locus " << loc->id << "\n";
	// }
    }

    return 0;
}

int
prune_nucleotides(CSLocus *cloc, Locus *loc, ofstream &log_fh, unsigned long int &nuc_cnt, 
		  unsigned long int &unk_hom_cnt, unsigned long int &unk_het_cnt, 
		  unsigned long int &hom_unk_cnt, unsigned long int &het_unk_cnt, 
		  unsigned long int &hom_het_cnt, unsigned long int &het_hom_cnt)
{
    map<char, int>      nucs;
    set<char>           cnucs;
    set<char>::iterator it;
    set<int>            rows;

    nuc_cnt += loc->len;

    for (uint i = 0; i < loc->snps.size() && i < cloc->snps.size(); i++) {
	//
	// Either their is an unknown call in locus, or, there is a snp in the catalog and any state in the locus.
	//
	if ((loc->snps[i]->type == snp_type_unk) || 
	    (cloc->snps[i]->type == snp_type_het && loc->snps[i]->type == snp_type_hom)) {

	    // cerr << "  Looking at SNP call in tag " << loc->id << " at position " << i << "; col: " << loc->snps[i]->col << "\n"
	    // 	 << "    Catalog column: " << cloc->snps[i]->col << " (" << i << "); Sample column: " << loc->snps[i]->col << " (" << i << ")\n"
	    // 	 << "    Sample has model call type: " << (loc->snps[i]->type == snp_type_unk ? "Unknown" : "Homozygous") << "; nucleotides: '" 
	    // 	 << loc->snps[i]->rank_1 << "' and '" << loc->snps[i]->rank_2 << "'\n"
	    // 	 << "    Catalog has model call type: " << (cloc->snps[i]->type == snp_type_het ? "Heterozygous" : "Homozygous") << "; nucleotides: '" 
	    // 	 << cloc->snps[i]->rank_1 << "' and '" << cloc->snps[i]->rank_2 << "'\n";

	    if (loc->snps[i]->rank_1 == 'N' || cloc->snps[i]->rank_1 == 'N') continue;

	    cnucs.insert(cloc->snps[i]->rank_1);
	    if (cloc->snps[i]->rank_2 != 0) cnucs.insert(cloc->snps[i]->rank_2);
	    if (cloc->snps[i]->rank_3 != 0) cnucs.insert(cloc->snps[i]->rank_3);
	    if (cloc->snps[i]->rank_4 != 0) cnucs.insert(cloc->snps[i]->rank_4);

	    // cerr << "    Catalog has nucleotides: ";
	    // for (it = cnucs.begin(); it != cnucs.end(); it++)
	    // 	cerr << *it << ", ";
	    // cerr << "\n";

	    //
	    // Tally the number of occurances of each nucleotide also present in the
	    // catalog in order to fuel the snp calling model.
	    //
	    // Note reads that contain nucleotides not present in the catalog so they
	    // can be excluded when calling haplotypes from the read.
	    //
	    nucs['A'] = 0;
	    nucs['C'] = 0;
	    nucs['G'] = 0;
	    nucs['T'] = 0;
	    nucs['N'] = 0;

	    for (uint k = 0; k < loc->reads.size(); k++) {
	    	if (cnucs.count(loc->reads[k][i]) > 0)
	    	    nucs[loc->reads[k][i]]++;
		else if (loc->reads[k][i] != 'N')
		    rows.insert(k);
	    }

	    //
	    // Test pruned data for homozygosity or heterozygosity.
	    //
	    invoke_model(loc, i, nucs);

	    log_model_calls(loc, log_fh,
			    unk_hom_cnt, unk_het_cnt, 
			    hom_unk_cnt, het_unk_cnt, 
			    hom_het_cnt, het_hom_cnt);
	}

	nucs.clear();
	cnucs.clear();
    }

    //
    // Re-call alleles.
    //
    loc->alleles.clear();

    call_alleles(loc, rows);

    //
    // If SNPs were called at this locus but no alleles could be determined,
    // blacklist this tag. This can occur if a fixed nucleotide isn't captured in
    // the catalog and all the reads are removed for the purpose of reading haplotypes.
    //
    if (loc->alleles.size() <= 1)
	for (uint j = 0; j < loc->snps.size(); j++)
	    if (loc->snps[j]->type == snp_type_het) {
		loc->blacklisted = 1;
		break;
	    }

    return 0;
}

int 
invoke_model(Locus *loc, int col, map<char, int> &nucs) 
{
    //
    // Search this column for the presence of a SNP
    //
    switch(model_type) {
    case snp:
	call_multinomial_snp(loc, col, nucs);
	break;
    case bounded:
	call_bounded_multinomial_snp(loc, col, nucs);
	break;
    default: 
	break;
    }

    return 0;
}

int 
call_alleles(Locus *loc, set<int> &rows) 
{
    int      row;
    int      height = loc->reads.size();
    string   allele;
    char     base;
    vector<SNP *>::iterator snp;

    for (row = 0; row < height; row++) {
	//
	// If a read had a nucleotide not present in the catalog, do not call
	// a haplotype from it.
	//
	if (rows.count(row) > 0)
	    continue;

	allele.clear();

	uint snp_cnt = 0;

	for (snp = loc->snps.begin(); snp != loc->snps.end(); snp++) {
	    if ((*snp)->type != snp_type_het) continue;

	    snp_cnt++;

	    base = loc->reads[row][(*snp)->col];

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
	    loc->alleles[allele]++;
    }

    return 0;
}

int
log_model_calls(Locus *loc, ofstream &log_fh,
		unsigned long int &unk_hom_cnt, unsigned long int &unk_het_cnt, 
		unsigned long int &hom_unk_cnt, unsigned long int &het_unk_cnt, 
		unsigned long int &hom_het_cnt, unsigned long int &het_hom_cnt)
{
    //
    // Log model call changes
    //
    for (uint j = 0; j < loc->snps.size(); j++) {
	switch(loc->model[j]) {
	case 'U':
	    switch(loc->snps[j]->type) {
	    case snp_type_het:
		if (verbose)
		    log_fh << loc->sample_id << "\t" << loc->id << "\t" << loc->snps[j]->col << "\t" << 'U' << "\t" << 'E' << "\n";
		unk_het_cnt++;
		break;
	    case snp_type_hom:
		if (verbose)
		    log_fh << loc->sample_id << "\t" << loc->id << "\t" << loc->snps[j]->col << "\t" << 'U' << "\t" << 'O' << "\n";
		unk_hom_cnt++;
		break;
	    case snp_type_unk:
	    default:
		break;
	    }
	    break;
	case 'E':
	    switch(loc->snps[j]->type) {
	    case snp_type_het:
		break;
	    case snp_type_hom:
		if (verbose)
		    log_fh << loc->sample_id << "\t" << loc->id << "\t" << loc->snps[j]->col << "\t" << 'E' << "\t" << 'O' << "\n";
		het_hom_cnt++;
		break;
	    case snp_type_unk:
	    default:
	        if (verbose)
		    log_fh << loc->sample_id << "\t" << loc->id << "\t" << loc->snps[j]->col << "\t" << 'E' << "\t" << 'U' << "\n";
		het_unk_cnt++;
		break;
	    }
	    break;
	case 'O':
	default:
	    switch(loc->snps[j]->type) {
	    case snp_type_het:
		if (verbose)
		    log_fh << loc->sample_id << "\t" << loc->id << "\t" << loc->snps[j]->col << "\t" << 'O' << "\t" << 'E' << "\n";
		hom_het_cnt++;
		break;
	    case snp_type_hom:
		break;
	    case snp_type_unk:
	    default:
		if (verbose)
		    log_fh << loc->sample_id << "\t" << loc->id << "\t" << loc->snps[j]->col << "\t" << 'O' << "\t" << 'U' << "\n";
		hom_unk_cnt++;
		break;
	    }
	    break;
	}
    }

    return 0;
}

int 
write_results(string file, map<int, Locus *> &m) 
{
    map<int, Locus *>::iterator i;
    vector<char *>::iterator    j;
    vector<int>::iterator       k;
    map<string, int>::iterator  t;
    Locus *tag_1;

    //
    // Parse the input file name to create the output files
    //
    string tag_file = out_path + file + ".tags.tsv";
    string snp_file = out_path + file + ".snps.tsv";
    string all_file = out_path + file + ".alleles.tsv";

    // Open the output files for writing.
    std::ofstream tags(tag_file.c_str());
    std::ofstream snps(snp_file.c_str());
    std::ofstream alle(all_file.c_str());

    int wrote = 0;

    for (i = m.begin(); i != m.end(); i++) {
	tag_1 = i->second;

	wrote++;

	// First write the consensus sequence
	tags << "0" << "\t" 
	     << tag_1->sample_id << "\t" 
	     << tag_1->id << "\t" 
             << tag_1->loc.chr << "\t"
             << tag_1->loc.bp << "\t"
             << (tag_1->loc.strand == plus ? "+" : "-") << "\t"
	     << "consensus" << "\t" 
	     << "\t"
	     << "\t" 
	     << tag_1->con << "\t" 
	     << (tag_1->deleveraged == true ? "1" : "0")     << "\t"
	     << (tag_1->blacklisted == true ? "1" : "0")     << "\t"
	     << (tag_1->lumberjackstack == true ? "1" : "0") << "\t"
	     << tag_1->lnl << "\n";

	//
	// Write a sequence recording the output of the SNP model for each nucleotide.
	//
	tags << "0"              << "\t" 
	     << tag_1->sample_id << "\t" 
	     << tag_1->id        << "\t" 
             << "\t"
             << "\t"
             << "\t"
	     << "model"          << "\t"
	     << "\t"
	     << "\t";
	for (uint j = 0; j < tag_1->snps.size(); j++) {
	    switch(tag_1->snps[j]->type) {
	    case snp_type_het:
		tags << "E";
		break;
	    case snp_type_hom:
		tags << "O";
		break;
	    default:
		tags << "U";
		break;
	    }
	}
	tags << "\t"
	     << "\t"
	     << "\t"
	     << "\t"
	     << "\n";
	
	//
	// Now write out each read from this locus.
	//
	for (uint j = 0; j < tag_1->reads.size(); j++) {
	    tags << "0"                 << "\t" 
		 << tag_1->sample_id    << "\t" 
		 << tag_1->id           << "\t\t\t\t";

	    if (tag_1->comp_type[j] == primary)
		tags << "primary" << "\t";
	    else
		tags << "secondary" << "\t";

	    tags << tag_1->comp_cnt[j]  << "\t" 
		 << tag_1->comp[j]      << "\t" 
		 << tag_1->reads[j]     << "\t\t\t\t\n";
	}

	//
	// Write out the model calls for each nucleotide in this locus.
	//
	for (uint j = 0; j < tag_1->snps.size(); j++) {
	    snps << "0"                 << "\t"
		 << tag_1->sample_id    << "\t" 
		 << tag_1->id           << "\t" 
		 << tag_1->snps[j]->col << "\t";

	    switch(tag_1->snps[j]->type) {
	    case snp_type_het:
		snps << "E\t";
		break;
	    case snp_type_hom:
		snps << "O\t";
		break;
	    default:
		snps << "U\t";
		break;
	    }

	    snps << std::fixed   << std::setprecision(2)
		 << tag_1->snps[j]->lratio << "\t" 
		 << tag_1->snps[j]->rank_1 << "\t" 
		 << (tag_1->snps[j]->rank_2 == 0 ? '-' : tag_1->snps[j]->rank_2) << "\t\t\n";
	}

	//
	// Write the expressed alleles seen for the recorded SNPs and
	// the percentage of tags a particular allele occupies.
	//
        char pct[id_len];
	for (t = tag_1->alleles.begin(); t != tag_1->alleles.end(); t++) {
            sprintf(pct, "%.2f", ((t->second/double(tag_1->reads.size())) * 100));
	    alle << "0"              << "\t"
		 << tag_1->sample_id << "\t"
		 << tag_1->id        << "\t"
		 << t->first         << "\t"
		 << pct              << "\t"
		 << t->second        << "\n";
	}
    }

    tags.close();
    snps.close();
    alle.close();

    cerr << "wrote " << wrote << " loci.\n";

    return 0;
}

int build_file_list(vector<pair<int, string> > &files) {
    vector<string> parts;
    string f;

    //
    // Read all the files from the Stacks directory.
    //
    uint   pos;
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

	pos = file.rfind(".tags.tsv");
	if (pos < file.length())
	    files.push_back(make_pair(1, file.substr(0, pos)));
    }

    closedir(dir);

    if (files.size() == 0) {
	cerr << "Unable to locate any input files to process within '" << in_path << "'\n";
	return 0;
    }

    cerr << "Found " << files.size() << " input file(s).\n";

    return 1;
}

int
fill_catalog_snps(map<int, CSLocus *> &catalog)
{
    map<int, CSLocus *>::iterator it;
    CSLocus *cloc;

    for (it = catalog.begin(); it != catalog.end(); it++) {
	cloc = it->second;

	queue<SNP *> snps;
	for (uint j = 0; j < cloc->snps.size(); j++)
	    snps.push(cloc->snps[j]);

	cloc->snps.clear();

	for (uint j = 0; j < cloc->len; j++) {
	    if (snps.size() > 0 && snps.front()->col == j) {
		cloc->snps.push_back(snps.front());
		snps.pop();
	    } else {
		SNP *snp = new SNP;
		snp->type   = snp_type_hom;
		snp->col    = j;
		snp->lratio = 0;
		snp->rank_1 = cloc->con[j];
		snp->rank_2 = 0;
		cloc->snps.push_back(snp);
	    }
	}
    }

    return 0;
}

int
init_log(int argc, char **argv, ofstream &log_fh)
{
    //
    // Open the log file.
    //
    stringstream log;
    log << "batch_" << batch_id << ".rxstacks.log";
    string log_path = out_path + log.str();
    log_fh.open(log_path.c_str(), ofstream::out);

    if (log_fh.fail()) {
        cerr << "Error opening log file '" << log_path << "'\n";
	exit(1);
    }

    //
    // Obtain the current date.
    //
    time_t     rawtime;
    struct tm *timeinfo;
    char       date[32];
    time(&rawtime);
    timeinfo = localtime(&rawtime);
    strftime(date, 32, "%F %T", timeinfo);

    log_fh << "#";
    for (int i = 0; i < argc; i++)
	log_fh << " " << argv[i]; 
    log_fh << "\n" << "# rxstacks executed " << date;

    if (!verbose)
	log_fh << "\n" 
	       << "# Sample\t"
	       << "Confounded loci\t"
	       << "Total nucs\t"
	       << "Total nucs converted\t"
	       << "Unk to Hom\t"
	       << "Unk to Het\t"
	       << "Hom to Unk\t"
	       << "Het to Unk\t"
	       << "Hom to Het\t"
	       << "Het to Hom\t"
	       << "Pruned Haplotypes\n";

    return 0;
}

int 
parse_command_line(int argc, char* argv[]) 
{
    int c;
     
    while (1) {
	static struct option long_options[] = {
	    {"help",         no_argument,       NULL, 'h'},
            {"version",      no_argument,       NULL, 'v'},
	    {"conf_filter",  no_argument,       NULL, 'F'},
	    {"prune_haplo",  no_argument,       NULL, 'H'},
	    {"lnl_filter",   no_argument,       NULL, 'G'},
	    {"lnl_dist",     no_argument,       NULL, 'D'},
	    {"verbose",      no_argument,       NULL, 'V'},
	    {"num_threads",  required_argument, NULL, 't'},
	    {"batch_id",     required_argument, NULL, 'b'},
	    {"in_path",      required_argument, NULL, 'P'},
	    {"outpath",      required_argument, NULL, 'o'},
	    {"model_type",   required_argument, NULL, 'T'},
	    {"bound_low",    required_argument, NULL, 'L'},
	    {"bound_high",   required_argument, NULL, 'U'},
	    {"alpha",        required_argument, NULL, 'A'},
	    {"conf_lim",     required_argument, NULL, 'C'},
	    {"max_haplo",    required_argument, NULL, 'M'},
	    {"lnl_lim",      required_argument, NULL, 'I'},
	    {0, 0, 0, 0}
	};
	
	// getopt_long stores the option index here.
	int option_index = 0;
     
	c = getopt_long(argc, argv, "hvVFGDHo:t:b:P:T:L:U:A:C:I:", long_options, &option_index);
     
	// Detect the end of the options.
	if (c == -1)
	    break;
     
	switch (c) {
	case 'h':
	    help();
	    break;
	case 't':
	    num_threads = atoi(optarg);
	    break;
	case 'P':
	    in_path = optarg;
	    break;
	case 'b':
	    batch_id = is_integer(optarg);
	    if (batch_id < 0) {
		cerr << "Batch ID (-b) must be an integer, e.g. 1, 2, 3\n";
		help();
	    }
	    break;
	case 'o':
	    out_path = optarg;
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
	case 'L':
	    bound_low  = atof(optarg);
	    break;
	case 'U':
	    bound_high = atof(optarg);
	    break;
	case 'A':
	    alpha = atof(optarg);
	    break;
	case 'F':
	    filter_confounded = true;
	    break;
	case 'C':
	    confounded_limit  = is_double(optarg);
	    filter_confounded = true;
	    break;
	case 'H':
	    prune_haplotypes = true;
	    break;
	case 'M':
	    max_haplotype_cnt = is_integer(optarg);
	    break;
	case 'G':
	    filter_lnl = true;
	    break;
	case 'I':
	    lnl_limit  = is_double(optarg);
	    filter_lnl = true;
	    break;
	case 'D':
	    lnl_dist = true;
	    break;
	case 'V':
	    verbose = true;
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

    if (out_path.length() == 0) 
	out_path = in_path;

    if (out_path.at(out_path.length() - 1) != '/') 
	out_path += "/";

    if (batch_id == 0) {
	cerr << "You must specify a batch ID.\n";
	help();
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

    if (filter_confounded == true && (confounded_limit < 0 || confounded_limit > 1.0)) {
	cerr << "Confounded locus limit is a percentage and must be between 0.0 and 1.0.\n";
	help();
    }

    return 0;
}

void version() {
    std::cerr << "rxstacks " << VERSION << "\n\n";

    exit(0);
}

void help() {
    std::cerr << "rxstacks " << VERSION << "\n"
              << "rxstacks -b batch_id -P path [-o path] [-t threads] [-v] [-h]" << "\n"
	      << "  b: Batch ID to examine when exporting from the catalog.\n"
	      << "  P: path to the Stacks output files.\n"
	      << "  o: output path to write results.\n"
	      << "  t: number of threads to run in parallel sections of code.\n"
	      << "  v: print program version." << "\n"
	      << "  h: display this help messsage." << "\n\n"
	      << "  Filtering options:\n"
	      << "    --lnl_filter: filter catalog loci based on the mean log likelihood of the catalog locus in the population.\n"
	      << "      --lnl_lim <limit>: minimum log likelihood required to keep a catalog locus.\n"
	      << "      --lnl_dist: print distribution of mean log likelihoods for catalog loci.\n"
	      << "    --conf_filter: filter confounded loci.\n"
	      << "    --conf_lim <limit>: between 0.0 and 1.0 (default 0.75), proportion of loci in population that must be confounded relative to the catalog locus.\n"
	      << "    --prune_haplo: prune out non-biological haplotypes unlikely to occur in the population.\n"
	      << "    --max_haplo_cnt <limit>: only consider haplotypes for pruning if they occur in fewer than max_haplo_cnt samples.\n"
	      << "  Model options:\n" 
	      << "    --model_type <type>: either 'snp' (default), 'bounded', or 'fixed'\n"
	      << "    For the SNP or Bounded SNP model:\n"
	      << "      --alpha <num>: chi square significance level required to call a heterozygote or homozygote, either 0.1 (default), 0.05, 0.01, or 0.001.\n"
	      << "    For the Bounded SNP model:\n"
	      << "      --bound_low <num>: lower bound for epsilon, the error rate, between 0 and 1.0 (default 0).\n"
	      << "      --bound_high <num>: upper bound for epsilon, the error rate, between 0 and 1.0 (default 1).\n"
	      << "  Logging Options:\n"
	      << "      --verbose: extended logging, including coordinates of all changed nucleotides (forces single-threaded execution).\n";
    exit(0);
}
