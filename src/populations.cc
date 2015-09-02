// -*-mode:c++; c-style:k&r; c-basic-offset:4;-*-
//
// Copyright 2012-2015, Julian Catchen <jcatchen@illinois.edu>
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

#include "populations.h"

// Global variables to hold command-line options.
int       num_threads =  1;
int       batch_id    = -1;
string    in_path;
string    out_path;
string    out_file;
string    pmap_path;
string    bl_file;
string    wl_file;
string    bs_wl_file;
string    enz;
double    sigma             = 150000.0;
double    sample_limit      = 0.0;
int       population_limit  = 1;
bool      calc_fstats       = false;
bool      bootstrap         = false;
bool      bootstrap_fst     = false;
bool      bootstrap_pifis   = false;
bool      bootstrap_phist   = false;
bool      bootstrap_div     = false;
bs_type   bootstrap_type    = bs_exact;
int       bootstrap_reps    = 100;
bool      bootstrap_wl      = false;
bool      write_single_snp  = false;
bool      write_random_snp  = false;
bool      merge_sites       = false;
bool      expand_id         = false;
bool      sql_out           = false;
bool      vcf_out           = false;
bool      vcf_haplo_out     = false;
bool      fasta_out         = false;
bool      fasta_strict_out  = false;
bool      genepop_out       = false;
bool      genomic_out       = false;
bool      structure_out     = false;
bool      phase_out         = false;
bool      fastphase_out     = false;
bool      beagle_out        = false;
bool      beagle_phased_out = false;
bool      plink_out         = false;
bool      hzar_out          = false;
bool      treemix_out       = false;
bool      phylip_out        = false;
bool      phylip_var        = false;
bool      phylip_var_all    = false;
bool      ordered_export    = false;
bool      kernel_smoothed   = false;
bool      loci_ordered      = false;
bool      log_fst_comp      = false;
bool      verbose           = false;
bool      filter_lnl        = false;
double    lnl_limit         = 0.0;
int       min_stack_depth   = 0;
double    merge_prune_lim   = 1.0;
double    minor_allele_freq = 0.0;
double    max_obs_het       = 1.0;
double    p_value_cutoff    = 0.05;
corr_type fst_correction    = no_correction;

map<int, string>          pop_key, grp_key;
map<int, pair<int, int> > pop_indexes;
map<int, vector<int> >    grp_members;

set<int> blacklist, bootstraplist;
map<int, set<int> > whitelist;

//
// Hold information about restriction enzymes
//
map<string, const char **> renz;
map<string, int>           renz_cnt;
map<string, int>           renz_len;
map<string, int>           renz_olap;

int main (int argc, char* argv[]) {

    initialize_renz(renz, renz_cnt, renz_len);
    initialize_renz_olap(renz_olap);

    parse_command_line(argc, argv);

    cerr
	<< "Fst kernel smoothing: " << (kernel_smoothed == true ? "on" : "off") << "\n"
	<< "Bootstrap resampling: ";
    if (bootstrap)
	cerr << "on, " << (bootstrap_type == bs_exact ? "exact; " : "approximate; ") << bootstrap_reps << " reptitions\n";
    else
	cerr << "off\n";
    cerr
	<< "Percent samples limit per population: " << sample_limit << "\n"
	<< "Locus Population limit: " << population_limit << "\n"
	<< "Minimum stack depth: " << min_stack_depth << "\n"
	<< "Log liklihood filtering: " << (filter_lnl == true ? "on"  : "off") << "; threshold: " << lnl_limit << "\n"
	<< "Minor allele frequency cutoff: " << minor_allele_freq << "\n"
	<< "Maximum observed heterozygosity cutoff: " << max_obs_het << "\n"
	<< "Applying Fst correction: ";
    switch(fst_correction) {
    case p_value:
	cerr << "P-value correction.\n";
	break;
    case bonferroni_win:
	cerr << "Bonferroni correction within sliding window.\n";
	break;
    case bonferroni_gen:
	cerr << "Bonferroni correction across genome wide sites.\n";
	break;
    case no_correction:
	cerr << "none.\n";
	break;
    }

    //
    // Set the number of OpenMP parallel threads to execute.
    //
    #ifdef _OPENMP
    omp_set_num_threads(num_threads);
    #endif

    //
    // Seed the random number generator
    //
    srandom(time(NULL));

    vector<pair<int, string> > files;
    if (!build_file_list(files, pop_indexes, grp_members))
	exit(1);

    if (wl_file.length() > 0) {
	load_marker_column_list(wl_file, whitelist);
	cerr << "Loaded " << whitelist.size() << " whitelisted markers.\n";
    }
    if (bl_file.length() > 0) {
	load_marker_list(bl_file, blacklist);
	cerr << "Loaded " << blacklist.size() << " blacklisted markers.\n";
    }
    if (bs_wl_file.length() > 0) {
	load_marker_list(bs_wl_file, bootstraplist);
	cerr << "Loaded " << bootstraplist.size() << " markers to include when bootstrapping.\n";
    }

    //
    // Open the log file.
    //
    stringstream log;
    log << "batch_" << batch_id << ".populations.log";
    string log_path = in_path + log.str();
    ofstream log_fh(log_path.c_str(), ofstream::out);
    if (log_fh.fail()) {
        cerr << "Error opening log file '" << log_path << "'\n";
	exit(1);
    }
    init_log(log_fh, argc, argv);

    //
    // Load the catalog
    //
    stringstream catalog_file;
    map<int, CSLocus *> catalog;
    bool compressed = false;
    int  res;
    catalog_file << in_path << "batch_" << batch_id << ".catalog";
    if ((res = load_loci(catalog_file.str(), catalog, false, false, compressed)) == 0) {
    	cerr << "Unable to load the catalog '" << catalog_file.str() << "'\n";
     	return 0;
    }

    //
    // Check the whitelist.
    //
    check_whitelist_integrity(catalog, whitelist);

    //
    // Implement the black/white list
    //
    reduce_catalog(catalog, whitelist, blacklist);

    //
    // If the catalog is not reference aligned, assign an arbitrary ordering to catalog loci.
    //
    loci_ordered = order_unordered_loci(catalog);

    //
    // Load matches to the catalog
    //
    vector<vector<CatMatch *> > catalog_matches;
    map<int, string>            samples;
    vector<int>                 sample_ids;
    for (int i = 0; i < (int) files.size(); i++) {
	vector<CatMatch *> m;
	load_catalog_matches(in_path + files[i].second, m);

	if (m.size() == 0) {
	    cerr << "Warning: unable to find any matches in file '" << files[i].second << "', excluding this sample from population analysis.\n";
	    //
	    // This case is generated by an existing, but empty file.
	    // Remove this sample from the population index which was built from 
	    // existing files, but we couldn't yet check for empty files.
	    //
	    map<int, pair<int, int> >::iterator pit;
	    for (pit = pop_indexes.begin(); pit != pop_indexes.end(); pit++)
		if (i >= pit->second.first && i <= pit->second.second) {
		    pit->second.second--;
		    pit++;
		    while (pit != pop_indexes.end()) {
			pit->second.first--;
			pit->second.second--;
			pit++;
		    }
		    break;
		}

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
    // Tabulate haplotypes present and in what combinations.
    //
    tabulate_haplotypes(catalog, pmap);

    //
    // Output a list of heterozygous loci and the associate haplotype frequencies.
    //
    if (sql_out)
    	write_sql(catalog, pmap);

    log_fh << "# Distribution of population loci.\n";
    log_haplotype_cnts(catalog, log_fh);

    apply_locus_constraints(catalog, pmap, pop_indexes, log_fh);

    log_fh << "# Distribution of population loci after applying locus constraints.\n";
    log_haplotype_cnts(catalog, log_fh);

    cerr << "Loading model outputs for " << sample_ids.size() << " samples, " << catalog.size() << " loci.\n";
    map<int, CSLocus *>::iterator it;
    map<int, ModRes *>::iterator mit;
    Datum   *d;
    CSLocus *loc;

    //
    // Load the output from the SNP calling model for each individual at each locus. This
    // model output string looks like this:
    //   OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOEOOOOOOEOOOOOOOOOOOOOOOOOOOOOOOOOOOOOUOOOOUOOOOOO
    // and records model calls for each nucleotide: O (hOmozygous), E (hEterozygous), U (Unknown)
    //
    for (uint i = 0; i < sample_ids.size(); i++) {
    	map<int, ModRes *> modres;
    	load_model_results(in_path + samples[sample_ids[i]], modres);

    	if (modres.size() == 0) {
    	    cerr << "Warning: unable to find any model results in file '" << samples[sample_ids[i]] << "', excluding this sample from population analysis.\n";
    	    continue;
    	}

    	for (it = catalog.begin(); it != catalog.end(); it++) {
    	    loc = it->second;
    	    d = pmap->datum(loc->id, sample_ids[i]);

    	    if (d != NULL) {
		if (modres.count(d->id) == 0) {
		    cerr << "Fatal error: Unable to find model data for catalog locus " << loc->id 
			 << ", sample ID " << sample_ids[i] << ", sample locus " << d->id 
			 << "; likely IDs were mismatched when running pipeline.\n";
		    exit(0);
		}
		d->len   = strlen(modres[d->id]->model);
    		d->model = new char[d->len + 1];
    		strcpy(d->model, modres[d->id]->model);
    	    }
    	}

    	for (mit = modres.begin(); mit != modres.end(); mit++)
    	    delete mit->second;
    	modres.clear();
    }

    uint pop_id, start_index, end_index;
    map<int, pair<int, int> >::iterator pit;

    PopSum<CSLocus> *psum = new PopSum<CSLocus>(pmap->loci_cnt(), pop_indexes.size());
    psum->initialize(pmap);

    for (pit = pop_indexes.begin(); pit != pop_indexes.end(); pit++) {
    	start_index = pit->second.first;
    	end_index   = pit->second.second;
    	pop_id      = pit->first;
    	cerr << "Generating nucleotide-level summary statistics for population '" << pop_key[pop_id] << "'\n";
    	psum->add_population(catalog, pmap, pop_id, start_index, end_index, verbose, log_fh);
    }

    cerr << "Tallying loci across populations...";
    psum->tally(catalog);
    cerr << "done.\n";

    //
    // We have removed loci that were below the -r and -p thresholds. Now we need to
    // identify individual SNPs that are below the -r threshold or the minor allele
    // frequency threshold (-a). In these cases we will remove the SNP, but keep the locus.
    //
    blacklist.clear();
    int pruned_snps = prune_polymorphic_sites(catalog, pmap, psum, pop_indexes, whitelist, blacklist, log_fh);
    cerr << "Pruned " << pruned_snps << " variant sites due to filter constraints.\n";

    if (!verbose)
	cerr << "  (enable the --verbose flag to record the reason why each site was filtered in the batch_X.populations.log file.)\n";
    
    //
    // Create an artificial whitelist if the user requested only the first or a random SNP per locus.
    //
    if (write_single_snp)
	implement_single_snp_whitelist(catalog, psum, whitelist);
    else if (write_random_snp)
	implement_random_snp_whitelist(catalog, psum, whitelist);

    //
    // Remove the accumulated SNPs 
    //
    cerr << "Removing " << blacklist.size() << " additional loci for which all variant sites were filtered...";
    set<int> empty_list;
    reduce_catalog(catalog, empty_list, blacklist);
    reduce_catalog_snps(catalog, whitelist, pmap);
    int retained = pmap->prune(blacklist);
    cerr << " retained " << retained << " loci.\n";

    //
    // Merge loci that overlap on a common restriction enzyme cut site.
    //
    map<int, pair<merget, int> > merge_map;
    if (merge_sites && loci_ordered)
	merge_shared_cutsite_loci(catalog, pmap, psum, merge_map, log_fh);

    //
    // Regenerate summary statistics after pruning SNPs and  merging loci.
    //
    delete psum;
    psum = new PopSum<CSLocus>(pmap->loci_cnt(), pop_indexes.size());
    psum->initialize(pmap);
    for (pit = pop_indexes.begin(); pit != pop_indexes.end(); pit++) {
	start_index = pit->second.first;
	end_index   = pit->second.second;
	pop_id      = pit->first;
	cerr << "Regenerating nucleotide-level summary statistics for population '" << pop_key[pop_id] << "'\n";
	psum->add_population(catalog, pmap, pop_id, start_index, end_index, verbose, log_fh);
    }
    cerr << "Re-tallying loci across populations...";
    psum->tally(catalog);
    cerr << "done.\n";

    for (pit = pop_indexes.begin(); pit != pop_indexes.end(); pit++) {
	pop_id = pit->first;
	if (kernel_smoothed && loci_ordered) {
	    cerr << "  Generating kernel-smoothed population statistics...\n";
	    kernel_smoothed_popstats(catalog, pmap, psum, pop_id, log_fh);
	}
    }

    calculate_haplotype_stats(files, pop_indexes, catalog, pmap, psum);

    if (calc_fstats) {
    	calculate_haplotype_divergence(files, pop_indexes, grp_members, catalog, pmap, psum);

    	calculate_haplotype_divergence_pairwise(files, pop_indexes, grp_members, catalog, pmap, psum);
    }

    //
    // Calculate and output the locus-level summary statistics.
    //
    calculate_summary_stats(files, pop_indexes, catalog, pmap, psum);

    //
    // Output the observed haplotypes.
    //
    write_generic(catalog, pmap, samples, false);

    //
    // Output data in requested formats
    //
    if (fasta_out)
    	write_fasta(catalog, pmap, samples, sample_ids);

    if (fasta_strict_out)
	write_strict_fasta(catalog, pmap, samples, sample_ids);

    if (genepop_out && ordered_export)
    	write_genepop_ordered(catalog, pmap, psum, pop_indexes, samples, log_fh);
    else if (genepop_out)
    	write_genepop(catalog, pmap, psum, pop_indexes, samples);

    if (structure_out && ordered_export)
    	write_structure_ordered(catalog, pmap, psum, pop_indexes, samples, log_fh);
    else if (structure_out) 
    	write_structure(catalog, pmap, psum, pop_indexes, samples);

    if (fastphase_out)
    	write_fastphase(catalog, pmap, psum, pop_indexes, samples);

    if (phase_out)
    	write_phase(catalog, pmap, psum, pop_indexes, samples);

    if (beagle_out)
    	write_beagle(catalog, pmap, psum, pop_indexes, samples);

    if (beagle_phased_out)
    	write_beagle_phased(catalog, pmap, psum, pop_indexes, samples);

    if (plink_out)
    	write_plink(catalog, pmap, psum, pop_indexes, samples);

    if (hzar_out)
    	write_hzar(catalog, pmap, psum, pop_indexes, samples);

    if (treemix_out)
    	write_treemix(catalog, pmap, psum, pop_indexes, samples);
    
    if (phylip_out || phylip_var)
    	write_phylip(catalog, pmap, psum, pop_indexes, samples);

    if (phylip_var_all)
    	write_fullseq_phylip(catalog, pmap, psum, pop_indexes, samples);

    if (vcf_haplo_out)
    	write_vcf_haplotypes(catalog, pmap, psum, samples, sample_ids);

    if (vcf_out && ordered_export)
    	write_vcf_ordered(catalog, pmap, psum, samples, sample_ids, merge_map, log_fh);
    else if (vcf_out)
    	write_vcf(catalog, pmap, psum, samples, sample_ids, merge_map);

    //
    // Calculate and write Fst.
    //
    if (calc_fstats)
    	write_fst_stats(files, pop_indexes, catalog, pmap, psum, log_fh);

    //
    // Output nucleotide-level genotype calls for each individual.
    //
    if (genomic_out)
    	write_genomic(catalog, pmap);

    log_fh.close();

    return 0;
}

int
apply_locus_constraints(map<int, CSLocus *> &catalog, 
			PopMap<CSLocus> *pmap, 
			map<int, pair<int, int> > &pop_indexes, 
			ofstream &log_fh)
{
    uint pop_id, start_index, end_index;
    CSLocus *loc;
    Datum  **d;

    if (sample_limit == 0 && population_limit == 0 && min_stack_depth == 0) return 0;

    if (verbose)
	log_fh << "\n#\n# List of loci removed by first filtering stage of sample and population constraints\n#\n"
	       << "# Action\tLocus ID\tChr\tBP\tColumn\tReason\n";

    map<int, CSLocus *>::iterator it;
    map<int, pair<int, int> >::iterator pit;

    uint pop_cnt   = pop_indexes.size();
    int *pop_order = new int [pop_cnt];

    // Which population each sample belongs to.
    int *samples   = new int [pmap->sample_cnt()];

    // For the current locus, how many samples in each population.
    int *pop_cnts  = new int [pop_cnt];

    // The total number of samples in each population.
    int *pop_tot   = new int [pop_cnt];

    pop_id = 0;
    for (pit = pop_indexes.begin(); pit != pop_indexes.end(); pit++) {
	start_index = pit->second.first;
	end_index   = pit->second.second;
	pop_tot[pop_id]  = 0;

	for (uint i = start_index; i <= end_index; i++) {
	    samples[i] = pop_id;
	    pop_tot[pop_id]++;
	}
	pop_order[pop_id] = pit->first;
	pop_id++;
    }

    for (uint i = 0; i < pop_cnt; i++)
	pop_cnts[i] = 0;

    double pct       = 0.0;
    bool   pop_limit = false;
    int    pops      = 0;
    int    below_stack_dep  = 0;
    uint   below_lnl_thresh = 0;
    set<int> blacklist;

    for (it = catalog.begin(); it != catalog.end(); it++) {
	loc = it->second;
	d   = pmap->locus(loc->id);

	for (int i = 0; i < pmap->sample_cnt(); i++) {
	    //
	    // Check that each sample is over the minimum stack depth for this locus.
	    //
	    if (d[i] != NULL && 
		min_stack_depth > 0 && 
		d[i]->tot_depth < min_stack_depth) {
		below_stack_dep++;
		delete d[i];
		d[i] = NULL;
		loc->hcnt--;
	    }

	    //
	    // Check that each sample is over the log likelihood threshold.
	    //
	    if (d[i] != NULL && 
		filter_lnl   && 
		d[i]->lnl < lnl_limit) {
		below_lnl_thresh++;
		delete d[i];
		d[i] = NULL;
		loc->hcnt--;
	    }
	}

	//
	// Tally up the count of samples in this population.
	//
	for (int i = 0; i < pmap->sample_cnt(); i++) {
	    if (d[i] != NULL)
		pop_cnts[samples[i]]++;
	}

	//
	// Check that the counts for each population are over sample_limit. If not, zero out 
	// the members of that population.
	//
	for (uint i = 0; i < pop_cnt; i++) {
	    pct = (double) pop_cnts[i] / (double) pop_tot[i];

	    if (pop_cnts[i] > 0 && pct < sample_limit) {
		//cerr << "Removing population " << pop_order[i] << " at locus: " << loc->id << "; below sample limit: " << pct << "\n";
		start_index = pop_indexes[pop_order[i]].first;
		end_index   = pop_indexes[pop_order[i]].second;

		for (uint j  = start_index; j <= end_index; j++) {
		    if (d[j] != NULL) {
			delete d[j];
			d[j] = NULL;
			loc->hcnt--;
		    }
		}
		pop_cnts[i] = 0;
	    }
	}

	//
	// Check that this locus is present in enough populations.
	//
	for (uint i = 0; i < pop_cnt; i++)
	    if (pop_cnts[i] > 0) pops++;
	if (pops < population_limit) {
	    //cerr << "Removing locus: " << loc->id << "; below population limit: " << pops << "\n";
	    pop_limit = true;
	}

	if (pop_limit) {
	    blacklist.insert(loc->id);

	    if (verbose)
		log_fh << "removed_locus\t"
		       << loc->id << "\t"
		       << loc->loc.chr << "\t"
		       << loc->sort_bp() << "\t"
		       << 0 << "\tfailed_population_limit\n";
	}

	for (uint i = 0; i < pop_cnt; i++)
	    pop_cnts[i] = 0;
	pop_limit = false;
	pops      = 0;
    }

    //
    // Remove loci
    //
    if (min_stack_depth > 0) 
	cerr << "Removed " << below_stack_dep << " samples from loci that are below the minimum stack depth of " << min_stack_depth << "x\n";
    if (filter_lnl)
	cerr << "Removed " << below_lnl_thresh << " samples from loci that are below the log likelihood threshold of " << lnl_limit << "\n";
    cerr << "Removing " << blacklist.size() << " loci that did not pass sample/population constraints...";
    set<int> whitelist;
    reduce_catalog(catalog, whitelist, blacklist);
    int retained = pmap->prune(blacklist);
    cerr << " retained " << retained << " loci.\n";

    delete [] pop_cnts;
    delete [] pop_tot;
    delete [] pop_order;
    delete [] samples;

    if (retained == 0)
	exit(0);

    return 0;
}

int
prune_polymorphic_sites(map<int, CSLocus *> &catalog, 
			PopMap<CSLocus> *pmap,
			PopSum<CSLocus> *psum,
			map<int, pair<int, int> > &pop_indexes, 
			map<int, set<int> > &whitelist, set<int> &blacklist,
			ofstream &log_fh)
{
    map<int, set<int> > new_wl;
    vector<int> pop_prune_list;
    CSLocus  *loc;
    LocTally *t;
    LocSum  **s;
    Datum   **d;
    bool      sample_prune, maf_prune, het_prune, inc_prune;
    int       size, pruned = 0;
    uint      pop_id, start_index, end_index;
    
    if (verbose)
	log_fh << "\n#\n# List of pruned nucleotide sites\n#\n"
	       << "# Action\tLocus ID\tChr\tBP\tColumn\tReason\n";

    //
    // If the whitelist is populated, use it as a guide for what loci to consider.
    // 
    // Construct a new whitelist along the way, that is a subset of the existing list.
    //
    if (whitelist.size() > 0) {
	map<int, set<int> >::iterator it;

	for (it = whitelist.begin(); it != whitelist.end(); it++) {
	    //
	    // A locus on the whitelist may have already been filtered out.
	    //
	    if (catalog.count(it->first) == 0)
		continue;

	    loc = catalog[it->first];
	    t   = psum->locus_tally(loc->id);
	    s   = psum->locus(loc->id);

	    //
	    // Check that each SNP in this locus is above the sample_limit and that 
	    // each SNP is above the minor allele frequency. If so, add it back to
	    // the whiteliest.
	    //
	    size = it->second.size();
	    for (uint i = 0; i < loc->snps.size(); i++) {

		//
		// If it is not already in the whitelist, ignore it.
		//
		if (size > 0 && it->second.count(loc->snps[i]->col) == 0)
		    continue;

		//
		// If the site is fixed, ignore it.
		//
		if (t->nucs[loc->snps[i]->col].fixed == true)
		    continue;

		sample_prune = false;
		maf_prune    = false;
		het_prune    = false;
		inc_prune    = false;
		pop_prune_list.clear();
		
		for (int j = 0; j < psum->pop_cnt(); j++) {
		    pop_id = psum->rev_pop_index(j);

		    if (s[j]->nucs[loc->snps[i]->col].incompatible_site)
			inc_prune = true;

		    else if (s[j]->nucs[loc->snps[i]->col].num_indv == 0 ||
			     (double) s[j]->nucs[loc->snps[i]->col].num_indv / (double) psum->pop_size(pop_id) < sample_limit)
			pop_prune_list.push_back(pop_id);
		}

		//
		// Check how many populations have to be pruned out due to sample limit. If less than
		// population limit, prune them; if more than population limit, mark locus for deletion.
		//
		if ((psum->pop_cnt() - pop_prune_list.size()) < (uint) population_limit) {
		    sample_prune = true;
		} else {
		    for (uint j = 0; j < pop_prune_list.size(); j++) {
			if (s[psum->pop_index(pop_prune_list[j])]->nucs[loc->snps[i]->col].num_indv == 0) continue;

		    	start_index = pop_indexes[pop_prune_list[j]].first;
		    	end_index   = pop_indexes[pop_prune_list[j]].second;
		    	d           = pmap->locus(loc->id);

			for (uint k = start_index; k <= end_index; k++) {
			    if (d[k] == NULL || loc->snps[i]->col >= (uint) d[k]->len) 
				continue;
			    if (d[k]->model != NULL) {
				d[k]->model[loc->snps[i]->col] = 'U';
			    }
			}
		    }
		}
		
		if (t->nucs[loc->snps[i]->col].allele_cnt > 1) {
		    //
		    // Test for minor allele frequency.
		    //
		    if ((1 - t->nucs[loc->snps[i]->col].p_freq) < minor_allele_freq)
			maf_prune = true;
		    //
		    // Test for observed heterozygosity.
		    //
		    if (t->nucs[loc->snps[i]->col].obs_het > max_obs_het)
		    	het_prune = true;
		}

		if (maf_prune == false && het_prune == false && sample_prune == false && inc_prune == false) {
		    new_wl[loc->id].insert(loc->snps[i]->col);
		} else {
		    pruned++;
		    if (verbose) {
			log_fh << "pruned_polymorphic_site\t"
			       << loc->id << "\t"
			       << loc->loc.chr << "\t"
			       << loc->sort_bp(loc->snps[i]->col) << "\t"
			       << loc->snps[i]->col << "\t"; 
			if (inc_prune)
			    log_fh << "incompatible_site\n";
			else if (sample_prune)
			    log_fh << "sample_limit\n";
			else if (maf_prune)
			    log_fh << "maf_limit\n";
			else if (het_prune)
			    log_fh << "obshet_limit\n";
			else
			    log_fh << "unknown_reason\n";
		    }
		}
	    }

	    //
	    // If no SNPs were retained for this locus, then mark it to be removed entirely.
	    //
	    if (new_wl.count(loc->id) == 0) {
		if (verbose)
		    log_fh << "removed_locus\t"
			   << loc->id << "\t"
			   << loc->loc.chr << "\t"
			   << loc->sort_bp() << "\t"
			   << 0 << "\tno_snps_remaining\n";
		blacklist.insert(loc->id);
	    }
	}

    } else {
	//
	// Otherwise, just iterate over the catalog.
	//
	map<int, CSLocus *>::iterator it;
	for (it = catalog.begin(); it != catalog.end(); it++) {
	    loc = it->second;

	    //
	    // If this locus is fixed, don't try to filter it out.
	    //
	    if (loc->snps.size() == 0) {
		new_wl.insert(make_pair(loc->id, std::set<int>()));
		continue;
	    }

	    t = psum->locus_tally(loc->id);
	    s = psum->locus(loc->id);

	    for (uint i = 0; i < loc->snps.size(); i++) {

		//
		// If the site is fixed, ignore it.
		//
		if (t->nucs[loc->snps[i]->col].fixed == true)
		    continue;

		sample_prune = false;
		maf_prune    = false;
		het_prune    = false;
		inc_prune    = false;
		pop_prune_list.clear();
		
		for (int j = 0; j < psum->pop_cnt(); j++) {
		    pop_id = psum->rev_pop_index(j);

		    if (s[j]->nucs[loc->snps[i]->col].incompatible_site)
			inc_prune = true;
		    else if (s[j]->nucs[loc->snps[i]->col].num_indv == 0 ||
			     (double) s[j]->nucs[loc->snps[i]->col].num_indv / (double) psum->pop_size(pop_id) < sample_limit)
			pop_prune_list.push_back(pop_id);
		}

		//
		// Check how many populations have to be pruned out due to sample limit. If less than
		// population limit, prune them; if more than population limit, mark locus for deletion.
		//
		if ((psum->pop_cnt() - pop_prune_list.size()) < (uint) population_limit) {
		    sample_prune = true;
		} else {
		    for (uint j = 0; j < pop_prune_list.size(); j++) {
			if (s[psum->pop_index(pop_prune_list[j])]->nucs[loc->snps[i]->col].num_indv == 0) continue;
			
		    	start_index = pop_indexes[pop_prune_list[j]].first;
		    	end_index   = pop_indexes[pop_prune_list[j]].second;
		    	d           = pmap->locus(loc->id);

			for (uint k = start_index; k <= end_index; k++) {
			    if (d[k] == NULL || loc->snps[i]->col >= (uint) d[k]->len) 
				continue;
			    if (d[k]->model != NULL) {
				d[k]->model[loc->snps[i]->col] = 'U';
			    }
			}
		    }
		}
		
		if (t->nucs[loc->snps[i]->col].allele_cnt > 1) {
		    //
		    // Test for minor allele frequency.
		    //
		    if ((1 - t->nucs[loc->snps[i]->col].p_freq) < minor_allele_freq)
			maf_prune = true;
		    //
		    // Test for observed heterozygosity.
		    //
		    if (t->nucs[loc->snps[i]->col].obs_het > max_obs_het)
		    	het_prune = true;
		}

		if (maf_prune == false && het_prune == false && sample_prune == false && inc_prune == false) {
		    new_wl[loc->id].insert(loc->snps[i]->col);
		} else {
		    pruned++;
		    if (verbose) {
			log_fh << "pruned_polymorphic_site\t"
			       << loc->id << "\t"
			       << loc->loc.chr << "\t"
			       << loc->sort_bp(loc->snps[i]->col) << "\t"
			       << loc->snps[i]->col << "\t";
			if (inc_prune)
			    log_fh << "incompatible_site\n";
			else if (sample_prune)
			    log_fh << "sample_limit\n";
			else if (maf_prune)
			    log_fh << "maf_limit\n";
			else if (het_prune)
			    log_fh << "obshet_limit\n";
			else
			    log_fh << "unknown_reason\n";
		    }
		}
	    }

	    //
	    // If no SNPs were retained for this locus, then mark it to be removed entirely.
	    //
	    if (new_wl.count(loc->id) == 0) {
		if (verbose)
		    log_fh << "removed_locus\t"
			   << loc->id << "\t"
			   << loc->loc.chr << "\t"
			   << loc->sort_bp() << "\t"
			   << 0 << "\tno_snps_remaining\n";
		blacklist.insert(loc->id);
	    }
	}
    }

    whitelist = new_wl;

    return pruned;
}

bool 
order_unordered_loci(map<int, CSLocus *> &catalog) 
{
    map<int, CSLocus *>::iterator it;
    CSLocus *loc;
    set<string> chrs;

    for (it = catalog.begin(); it != catalog.end(); it++) {
	loc = it->second;
	if (strlen(loc->loc.chr) > 0) 
	    chrs.insert(loc->loc.chr);
    }

    //
    // This data is already reference aligned.
    //
    if (chrs.size() > 0)
	return true;

    cerr << "Catalog is not reference aligned, arbitrarily ordering catalog loci.\n";

    uint bp = 1;
    for (it = catalog.begin(); it != catalog.end(); it++) {
	loc = it->second;
	loc->loc.chr = new char[3];
	strcpy(loc->loc.chr, "un");
	loc->loc.bp  = bp;

	bp += strlen(loc->con);
    }

    return false;
}

int
log_haplotype_cnts(map<int, CSLocus *> &catalog, ofstream &log_fh)
{
    map<int, CSLocus *>::iterator it;
    map<int, int> valid, absent, confounded;

    CSLocus *loc;
    int missing;

    for (it = catalog.begin(); it != catalog.end(); it++) {
	loc = it->second;

	if (valid.count(loc->hcnt) == 0)
	    valid[loc->hcnt] = 1;
	else
	    valid[loc->hcnt]++;

	if (confounded.count(loc->confounded_cnt) == 0)
	    confounded[loc->confounded_cnt] = 1;
	else
	    confounded[loc->confounded_cnt]++;

	missing = loc->cnt - loc->hcnt;

	if (absent.count(missing) == 0)
	    absent[missing] = 1;
	else
	    absent[missing]++;
    }

    map<int, int>::iterator cnt_it;

    log_fh << "# Distribution of valid loci matched to catalog locus.\n"
	   << "# Valid samples at locus\tCount\n";
    for (cnt_it = valid.begin(); cnt_it != valid.end(); cnt_it++)
	log_fh << cnt_it->first << "\t" << cnt_it->second << "\n";

    log_fh << "# Distribution of confounded loci at catalog locus.\n"
	   << "# Confounded samples at locus\tCount\n";
    for (cnt_it = confounded.begin(); cnt_it != confounded.end(); cnt_it++)
	log_fh << cnt_it->first << "\t" << cnt_it->second << "\n";

    log_fh << "# Distribution of missing loci at catalog loci.\n"
	   << "# Absent samples at locus\tCount\n";
    for (cnt_it = absent.begin(); cnt_it != absent.end(); cnt_it++)
	log_fh << cnt_it->first << "\t" << cnt_it->second << "\n";

    return 0;
}

int 
tabulate_haplotypes(map<int, CSLocus *> &catalog, PopMap<CSLocus> *pmap) 
{
    map<int, CSLocus *>::iterator it;
    vector<char *>::iterator hit;
    Datum  **d;
    CSLocus *loc;
    double   mean, cnt;

    for (it = catalog.begin(); it != catalog.end(); it++) {
	loc = it->second;
	d   = pmap->locus(loc->id);

	mean = 0.0;
	cnt  = 0.0;

	for (int i = 0; i < pmap->sample_cnt(); i++) {
	    if (d[i] == NULL) 
		continue;

	    if (d[i]->obshap.size() > 1)
		loc->marker = "heterozygous";

	    mean += d[i]->lnl;
	    cnt++;
	}

	if (loc->marker.length() > 0) {
     	    create_genotype_map(loc, pmap);
	    call_population_genotypes(loc, pmap);
	}

	loc->lnl = mean / cnt;
    }

    return 0;
}

int
merge_shared_cutsite_loci(map<int, CSLocus *> &catalog, 
			  PopMap<CSLocus> *pmap, PopSum<CSLocus> *psum,
			  map<int, pair<merget, int> > &merge_map,
			  ofstream &log_fh)
{
    map<string, vector<CSLocus *> >::iterator it;
    CSLocus *cur, *next;
    Datum  **d_1, **d_2;
    double   prune_pct;
    uint unmergable, tot_loci, tot_samp;
    uint success           = 0;
    uint failure           = 0;
    uint overlap           = 0;
    uint simple_merge_cnt  = 0;
    uint complex_merge_cnt = 0;
    uint missing_samps_cnt = 0;
    uint phase_fail_cnt    = 0;
    uint nomapping_cnt     = 0;
    uint multimapping_cnt  = 0;
    uint multifails_cnt    = 0;

    tot_loci = pmap->loci_cnt();

    set<int> loci_to_destroy;
    map<int, int> missing_samps_dist;

    cerr << "To merge adjacent loci at least " << merge_prune_lim * 100 << "% of samples must have both adjacent loci;"
	 << " the remaining " << 100 - (merge_prune_lim * 100) << "% of individuals will be pruned.\n"
	 << "Attempting to merge adjacent loci that share a cutsite...";
	 
    if (verbose)
	log_fh << "\n#\n# List of locus pairs that share a cutsite that failed to merge because they could not be phased.\n#\n";

    //
    // Iterate over each chromosome.
    //
    for (it = pmap->ordered_loci.begin(); it != pmap->ordered_loci.end(); it++) {
	//
	// Iterate over each ordered locus on this chromosome.
	//
	next = it->second[0];
	for (uint pos = 1; pos < it->second.size(); pos++) {
	    cur  = next;
	    next = it->second[pos];

	    //
	    // Do these two loci overlap?
	    //   +Must occur on opposite strands
	    //   +Must overlap according to the length of the cutsite.
	    //
	    if (((cur->loc.strand == minus && next->loc.strand == plus) &&
		 ((int) (cur->loc.bp  - next->loc.bp + 1) == renz_olap[enz])) ||
		((cur->loc.strand == plus  && next->loc.strand == minus) &&
		 ((int) (next->loc.bp - cur->loc.bp  + 1) == renz_olap[enz]))) {
		overlap++;

		d_1        = pmap->locus(cur->id);
		d_2        = pmap->locus(next->id);
		unmergable = 0;
		tot_samp   = 0;

		//
		// Check if all members of the population contain these two loci (or are missing both).
		//
		for (int i = 0; i < pmap->sample_cnt(); i++) {
		    if (d_1[i] != NULL || d_2[i] != NULL)
			tot_samp++;
		    if ((d_1[i] != NULL && d_2[i] == NULL) ||
			(d_1[i] == NULL && d_2[i] != NULL))
			unmergable++;
		}

		prune_pct = (double) (tot_samp - unmergable) / (double) tot_samp;

		//
		// If some of the individuals only have one locus and not the other, prune them out.
		//
		if (prune_pct < 1.0 && prune_pct >= merge_prune_lim) {
		    for (int i = 0; i < pmap->sample_cnt(); i++)
			if (d_1[i] != NULL && d_2[i] == NULL) {
			    delete d_1[i];
			    d_1[i] = NULL;
			} else if (d_1[i] == NULL && d_2[i] != NULL) {
			    delete d_2[i];
			    d_2[i] = NULL;
			}
		}

		//
		// If possible, merge the two loci together.
		//
		if (prune_pct < merge_prune_lim) {
		    int pct = (int) (prune_pct * 100);
		    missing_samps_dist[pct]++;
		    if (verbose) log_fh << "Missing samples, Sink Locus: " << cur->id << "; Source Locus: " << next->id << "; " 
					<< pct << "% present (" << 100 - pct << "% missing)\n";
		    missing_samps_cnt++;
		    failure++;
		    continue;
		} 

		phaset res = merge_and_phase_loci(pmap, cur, next, loci_to_destroy, log_fh);
		switch(res) {
		case multiple_fails:
		    if (verbose) log_fh << "Failed to phase, Sink Locus: " << cur->id << "; Source Locus: " << next->id << "; " 
					<< "multiple failures\n";
		    multifails_cnt++;
		    phase_fail_cnt++;
		    failure++;
		    break;
		case multimapping_fail:
		    if (verbose) log_fh << "Failed to phase, Sink Locus: " << cur->id << "; Source Locus: " << next->id << "; " 
					<< "multimapping in one or more individuals\n";
		    multimapping_cnt++;
		    phase_fail_cnt++;
		    failure++;
		    break;
		case nomapping_fail:
		    if (verbose) log_fh << "Failed to phase, Sink Locus: " << cur->id << "; Source Locus: " << next->id << "; " 
					<< "no mapping in one or more individuals\n";
		    nomapping_cnt++;
		    phase_fail_cnt++;
		    failure++;
		    break;
		case complex_phase:
		    if (verbose) log_fh << "Phased Sink Locus: " << cur->id << "; Source Locus: " << next->id << "; "
					<< "a complex phasing operation.\n";
		    complex_merge_cnt++;
		    success++;
		    merge_map[cur->id]  = make_pair(merge_sink, next->id);
		    merge_map[next->id] = make_pair(merge_src, cur->id);
		    break;
		case simple_merge:
		    if (verbose) log_fh << "Phased Sink Locus: " << cur->id << "; Source Locus: " << next->id << "; "
					<< "a simple merge operation.\n";
		    simple_merge_cnt++;
		    success++;
		    merge_map[cur->id]  = make_pair(merge_sink, next->id);
		    merge_map[next->id] = make_pair(merge_src, cur->id);
		    break;
		default:
		    cerr << "Warning: Merge failure.\n";
		    break;
		}
	    }
	}
    }

    //
    // Remove those loci that have been merged from both the popualtion map and catalog.
    //
    set<int> emptyset;
    pmap->prune(loci_to_destroy);
    reduce_catalog(catalog, emptyset, loci_to_destroy);

    cerr << "done.\n"
	 << "Of " << tot_loci << " loci, "
	 << overlap << " pairs share a cutsite; "
	 << success << " pairs were merged; "
	 << failure << " pairs failed to merge; "
	 << pmap->loci_cnt() << " remaining loci.\n"
	 << "  Of those merged, " << simple_merge_cnt << " required only a simple merge without phasing; "
	 << "while " << complex_merge_cnt << " required phasing.\n"
	 << "  Of those that failed to merge, " << missing_samps_cnt << " were missing one of the two haplotypes in one or more samples; "
	 << "while " << phase_fail_cnt << " failed to be phased.\n" 
	 << "    Of those that failed to phase, " << nomapping_cnt << " failed due to a lack of haplotype mappings; " 
	 << multimapping_cnt << " failed due to multiple haplotype mappings; " << multifails_cnt << " failed due to both.\n";

    log_fh << "\n#\n# Merging adjacent loci with a shared restriction enzyme cutsite\n#\n"
	   << "Of " << tot_loci << " loci, "
	   << overlap << " pairs share a cutsite; "
	   << success << " pairs were merged; "
	   << failure << " pairs failed to merge; "
	   << pmap->loci_cnt() << " remaining loci.\n"
	   << "  Of those merged, " << simple_merge_cnt << " required only a simple merge without phasing; "
	   << "while " << complex_merge_cnt << " required phasing.\n"
	   << "  Of those that failed to merge, " << missing_samps_cnt << " were missing one of the two haplotypes in one or more samples; "
	   << "while " << phase_fail_cnt << " failed to be phased.\n" 
	   << "    Of those that failed to phase, " << nomapping_cnt << " failed due to a lack of haplotype mappings; " 
	   << multimapping_cnt << " failed due to multiple haplotype mappings; " << multifails_cnt << " failed due to both.\n";
    log_fh << "#\n# Distribution of loci with samples missing one of two loci to be merged\n"
	   << "# Percent samples with both loci present\tNumber of cases\n";
    map<int, int>::iterator mit;
    for (mit = missing_samps_dist.begin(); mit != missing_samps_dist.end(); mit++)
	log_fh << mit->first << "\t" << mit->second << "\n";
    log_fh << "\n";

    return 0;
}

phaset
merge_and_phase_loci(PopMap<CSLocus> *pmap, CSLocus *cur, CSLocus *next, 
		     set<int> &loci_to_destroy,
		     ofstream &log_fh)
{
    Datum **d_1 = pmap->locus(cur->id);
    Datum **d_2 = pmap->locus(next->id);

    set<int>    phased_results;
    set<string> phased_haplotypes;
    string      merged_hap;
    char       *h_1, *h_2;
    int         merge_type;

    if (verbose) log_fh << "Attempting to phase source locus " << cur->id << " with sink locus " << next->id << "\n";

    int sample_cnt        = 0;
    int phased_sample_cnt = 0;
    //
    // Take a census of the already phased haplotypes. We have phased haplotypes 
    // if for individual i:
    //   1. d_1 has a single haplotype and d_2 has a single haplotype
    //   2. d_1 has a single haplotpye and d_2 has multiple haplotypes
    //   3. d_1 has multiple haplotpyes and d_2 has a single haplotype
    //
    // If one or both of the loci have no SNPs, then the haplotype is
    // recorded as "consensus." Check that condition before we start merging.
    //
    if (cur->snps.size() > 0 && next->snps.size() > 0)
	merge_type = 0;
    else if (cur->snps.size() == 0)
	merge_type = 1;
    else if (next->snps.size() == 0)
	merge_type = 2;
    else
	merge_type = 3;

    for (int i = 0; i < pmap->sample_cnt(); i++) {
	if (d_1[i] == NULL || d_2[i] == NULL)
	    continue;
	else if (d_1[i]->obshap.size() > 1 && d_2[i]->obshap.size() > 1)
	    continue;
	else {
	    for (uint j = 0; j < d_1[i]->obshap.size(); j++) {
		for (uint k = 0; k < d_2[i]->obshap.size(); k++) {
		    switch (merge_type) {
		    case 0:
			merged_hap = string(d_1[i]->obshap[j]) + string(d_2[i]->obshap[k]);
			break;
		    case 1:
			merged_hap = string(d_2[i]->obshap[k]);
			break;
		    case 2:
			merged_hap = string(d_1[i]->obshap[j]);
			break;
		    case 3:
		    default:
			merged_hap = "consensus";
			break;
		    }
		    phased_haplotypes.insert(merged_hap);
		    // cerr << "Phasing: '" << d_1[i]->obshap[j] << "' + '" << d_2[i]->obshap[k] << "' => '" << merged_hap << "'\n";
		}
	    }
	    phased_sample_cnt++;
	    sample_cnt++;
	}
    }

    //
    // Indicate that these two loci had a simple merge, with no phasing necessary.
    //
    phased_results.insert(simple_merge);

    //
    // Now we need to check if we can phase the remaining haplotypes.
    //
    for (int i = 0; i < pmap->sample_cnt(); i++) {
	if (d_1[i] == NULL || d_2[i] == NULL)
	    continue;
	else if (d_1[i]->obshap.size() > 1 && d_2[i]->obshap.size() > 1) {
	    // cerr << "Attempting to phase individual " << i << ": " << d_1[i]->id << " / " << d_2[i]->id << "\n";

	    sample_cnt++;
	    //
	    // We should be able to find a sinlge phasing mapping for each haplotype from d_1 to d_2 
	    // that includes all the haplotypes in these two loci.
	    //
	    vector<pair<char *, char *> > seen_phased;
	    uint tot_obshap = d_1[i]->obshap.size() + d_2[i]->obshap.size();
	    uint phased_cnt = 0;
	    for (uint j = 0; j < d_1[i]->obshap.size(); j++) {
		for (uint k = 0; k < d_2[i]->obshap.size(); k++) {
		    // cerr << "  " << d_1[i]->obshap[j] << " + " << d_2[i]->obshap[k];
		    //
		    // Record each pair of haplotypes that has been seen phased previously.
		    //
		    if (phased_haplotypes.count(string(d_1[i]->obshap[j]) + string(d_2[i]->obshap[k]))) {
			seen_phased.push_back(make_pair(d_1[i]->obshap[j], d_2[i]->obshap[k]));
			// cerr << " => " << d_1[i]->obshap[j] << d_2[i]->obshap[k];
		    }
		    // cerr << "\n";
		}
	    }
	    //
	    // Now, we will iterate over all sets of phased haplotypes and look
	    // for combinations that use all four individual haplotypes.
	    //
	    for (uint j = 0; j < seen_phased.size(); j++) {
		for (uint k = j; k < seen_phased.size(); k++) {
		    set<char *> incorporated_haplotypes;
		    //
		    // Count the number of distinct char pointers. If this combination
		    // of haplotypes includes all unphased haplotypes, count it.
		    //
		    incorporated_haplotypes.insert(seen_phased[j].first);
		    incorporated_haplotypes.insert(seen_phased[j].second);
		    incorporated_haplotypes.insert(seen_phased[k].first);
		    incorporated_haplotypes.insert(seen_phased[k].second);
		    if (incorporated_haplotypes.size() == tot_obshap)
			phased_cnt++;
		}
	    }

	    //
	    // If one pair of haplotypes is mapped, but the other is not, assume the second pair or
	    // haplotypes must be phased by process of elimination.
	    //
	    if (phased_cnt == 0 && seen_phased.size() == 1) {
	    	h_1 = seen_phased[0].first  == d_1[i]->obshap[1] ?
	    	    d_1[i]->obshap[0] : d_1[i]->obshap[1];
	    	h_2 = seen_phased[0].second == d_2[i]->obshap[1] ?
	    	    d_2[i]->obshap[0] : d_2[i]->obshap[1];
	    	phased_haplotypes.insert(string(h_1) + string(h_2));
	    	phased_cnt++;
	    	// cerr << "  Phasing: '" << hap_1 << "' + '" << hap_2 << "' => '" << string(hap_1) + string(hap_2) << "'\n";
	    }

	    if (phased_cnt == 0) {
		phased_results.insert(nomapping_fail);
		if (verbose) log_fh << "    Locus NOT phased in individual " << i << "; loci " << d_1[i]->id << " / " << d_2[i]->id << " no mapping found.\n";
	    } else if (phased_cnt == 1) {
		phased_sample_cnt++;
		phased_results.insert(complex_phase);
	    } else {
		phased_results.insert(multimapping_fail);
		if (verbose) log_fh << "    Locus NOT phased in individual " << i << "; loci " << d_1[i]->id << " / " << d_2[i]->id << " multiple mappings found.\n";
	    }
	}
    }

    if (phased_sample_cnt != sample_cnt) {
	if (phased_results.count(nomapping_fail) > 0 && 
	    phased_results.count(multimapping_fail) > 0)
	    return multiple_fails;
	else if (phased_results.count(nomapping_fail) > 0)
	    return nomapping_fail;
	else if (phased_results.count(multimapping_fail) > 0)
	    return multimapping_fail;
	else {
	    cerr << "WE SHOULD NOT GET HERE\n";
	    return merge_failure;
	}
    }

    //
    // Okay, merge these two loci together.
    //
    if (!merge_datums(pmap->sample_cnt(), cur->len, d_1, d_2, phased_haplotypes, merge_type))
	return merge_failure;

    //
    // Merge the catalog entries together.
    //
    if (!merge_csloci(cur, next, phased_haplotypes))
    	return merge_failure;

    //
    // Mark the merged locus for destruction.
    //
    loci_to_destroy.insert(next->id);

    if (phased_results.count(complex_phase) > 0)
	return complex_phase;
    return simple_merge;
}

int
merge_csloci(CSLocus *sink, CSLocus *src, set<string> &phased_haplotypes)
{
    //
    // We assume that we are merging two loci: one on the negative strand, one on the
    // positive. We will keep the sink cslocus and delete the src cslocus.
    //   -> The sink cslocus is assumed to be on the negative strand.
    //

    //
    // 1. Reverse complement the SNP coordinates in the sink locus so that they are 
    //    enumerated on the positive strand. Complement the alleles as well.
    //
    for (uint j = 0; j < sink->snps.size(); j++) {
	sink->snps[j]->col    = sink->len - sink->snps[j]->col - 1;
	sink->snps[j]->rank_1 = reverse(sink->snps[j]->rank_1);
	sink->snps[j]->rank_2 = reverse(sink->snps[j]->rank_2);
	sink->snps[j]->rank_3 = reverse(sink->snps[j]->rank_3);
	sink->snps[j]->rank_4 = reverse(sink->snps[j]->rank_4);
    }

    //
    // 2. Adjust the SNP coordinates in the src locus to account for the now, longer length.
    //
    for (uint j = 0; j < src->snps.size(); j++)
	src->snps[j]->col = sink->len + src->snps[j]->col - renz_olap[enz];

    //
    // 3. Combine SNPs between the two catalog loci: add the SNPs from the sink (formerly on the 
    //    negative strand) in reverse order, followed by the SNPs from the src.
    //
    vector<SNP *> tmpsnp;
    for (int j = (int) sink->snps.size() - 1; j >= 0; j--)
	tmpsnp.push_back(sink->snps[j]);
    for (uint j = 0; j < src->snps.size(); j++)
	tmpsnp.push_back(src->snps[j]);
    sink->snps.clear();
    for (uint j = 0; j < tmpsnp.size(); j++)
	sink->snps.push_back(tmpsnp[j]);

    //
    // 4. Adjust the genomic location of the sink locus.
    //
    uint bp = sink->sort_bp();
    sink->loc.bp     = bp;
    sink->loc.strand = plus;

    //
    // 5. Adjust the length of the sequence.
    //
    sink->len += src->len - renz_olap[enz];

    //
    // 6. Merge the consensus sequence together.
    //
    char *new_con = rev_comp(sink->con);
    delete [] sink->con;
    sink->con = new_con;
    new_con   = new char[sink->len + 1];
    strcpy(new_con, sink->con);
    delete [] sink->con;
    sink->con = new_con;
    new_con  += src->len - renz_olap[enz];
    strcpy(new_con, src->con);

    //
    // 7. Record the now phased haplotypes.
    //
    sink->alleles.clear();
    set<string>::iterator it;
    for (it = phased_haplotypes.begin(); it != phased_haplotypes.end(); it++)
	sink->alleles[*it] = 0;

    // cerr << "CSLocus " << sink->id << ":\n"
    // 	 << "Length: " << sink->len << "; Chr: " << sink->loc.chr << "; BP: " << sink->sort_bp() << "; strand: " << (sink->loc.strand == plus ? "+" : "-") << "\n"
    // 	 << "  SNPs:\n";
    // for (uint j = 0; j < sink->snps.size(); j++) 
    // 	cerr << "    Col: " << sink->snps[j]->col 
    // 	     << "    Rank 1: " << sink->snps[j]->rank_1
    // 	     << "    Rank 2: " << sink->snps[j]->rank_2 << "\n";
    // cerr << "  Alleles:\n";
    // map<string, int>::iterator ait;
    // for (ait = sink->alleles.begin(); ait != sink->alleles.end(); ait++) 
    // 	cerr << "    " << ait->first << "\n";

    return 1;
}

int
merge_datums(int sample_cnt, 
	     int sink_locus_len,
	     Datum **sink, Datum **src, 
	     set<string> &phased_haplotypes, 
	     int merge_type)
{
    char           tmphap[id_len], *new_hap;
    uint           haplen, model_len, offset;
    vector<SNP *>  tmpsnp;
    vector<string> tmpobshap;
    vector<int>    tmpobsdep;

    //
    // We assume that we are merging two loci: one on the negative strand, one on the
    // positive. We will keep the sink datum and delete the src datum.
    //   -The sink datum is assumed to be on the negative strand.
    //
    for (int i = 0; i < sample_cnt; i++) {
	if (sink[i] == NULL && src[i] == NULL)
	    continue;
	else if (sink[i] == NULL || src[i] == NULL)
	    cerr << "Unexpected condition in merging datums: one datum is NULL while the other is not.\n";

	//
	// 1. Reverse complement the SNP coordinates in the sink locus so that they are 
	//    enumerated on the positive strand. Complement the alleles as well.
	//
	for (uint j = 0; j < sink[i]->snps.size(); j++) {
	    sink[i]->snps[j]->col    = sink[i]->len - sink[i]->snps[j]->col - 1;
	    sink[i]->snps[j]->rank_1 = reverse(sink[i]->snps[j]->rank_1);
	    sink[i]->snps[j]->rank_2 = reverse(sink[i]->snps[j]->rank_2);
	    sink[i]->snps[j]->rank_3 = reverse(sink[i]->snps[j]->rank_3);
	    sink[i]->snps[j]->rank_4 = reverse(sink[i]->snps[j]->rank_4);
	}

	//
	// 2. Adjust the SNP coordinates in the src locus to account for the now, longer length.
	//
	for (uint j = 0; j < src[i]->snps.size(); j++)
	    src[i]->snps[j]->col = sink[i]->len + src[i]->snps[j]->col - renz_olap[enz];

	//
	// 3. Reverse complement the observed haplotypes in the sink locus.
	//
	haplen = strlen(sink[i]->obshap[0]);
	for (uint j = 0; j < sink[i]->obshap.size(); j++) {
	    for (uint k = 0; k < haplen; k++)
		tmphap[k] = reverse(sink[i]->obshap[j][haplen - k - 1]);
	    tmphap[haplen] = '\0';
	    strcpy(sink[i]->obshap[j], tmphap);
	}

	//
	// 4. Combine SNPs between the two datums: add the SNPs from the sink (formerly on the 
	//    negative strand) in reverse order, followed by the SNPs from the src.
	//
	tmpsnp.clear();
	for (int j = (int) sink[i]->snps.size() - 1; j >= 0; j--)
	    tmpsnp.push_back(sink[i]->snps[j]);
	for (uint j = 0; j < src[i]->snps.size(); j++)
	    tmpsnp.push_back(src[i]->snps[j]);
	sink[i]->snps.clear();
	for (uint j = 0; j < tmpsnp.size(); j++)
	    sink[i]->snps.push_back(tmpsnp[j]);
    }

    //
    // 5. Combine observed haplotypes between the two datums while phasing them.
    //    5.1 First combine the haplotypes from samples that are already in phase.
    //
    string      merged_hap;
    vector<int> to_be_phased;
    phased_haplotypes.clear();
    for (int i = 0; i < sample_cnt; i++) {
	if (sink[i] == NULL && src[i] == NULL)
	    continue;

	if (sink[i]->obshap.size() > 1 && src[i]->obshap.size() > 1) {
	    to_be_phased.push_back(i);
	    continue;
	} else {
	    tmpobshap.clear();
	    tmpobsdep.clear();
	    for (uint j = 0; j < sink[i]->obshap.size(); j++) {
		for (uint k = 0; k < src[i]->obshap.size(); k++) {
		    switch (merge_type) {
		    case 0:
			merged_hap = string(sink[i]->obshap[j]) + string(src[i]->obshap[k]);
			break;
		    case 1:
			merged_hap = string(src[i]->obshap[j]);
			break;
		    case 2:
			merged_hap = string(sink[i]->obshap[j]);
			break;
		    case 3:
		    default:
			merged_hap = "consensus";
			break;
		    }
		    phased_haplotypes.insert(merged_hap);
		    tmpobshap.push_back(merged_hap);
		    tmpobsdep.push_back((sink[i]->depth[j] + src[i]->depth[k]) / 2);
		}
	    }
	    sink[i]->depth.clear();
	    for (uint j = 0; j < sink[i]->obshap.size(); j++)
		delete [] sink[i]->obshap[j];
	    sink[i]->obshap.clear();
	    for (uint j = 0; j < tmpobshap.size(); j++) {
		new_hap = new char[tmpobshap[j].length() + 1];
		strcpy(new_hap, tmpobshap[j].c_str());
		sink[i]->obshap.push_back(new_hap);
		sink[i]->depth.push_back(tmpobsdep[j]);
	    }
	}
    }
    //
    //    5.2 Phase and combine the haplotypes from the remaining samples.
    //
    int index;
    for (uint i = 0; i < to_be_phased.size(); i++) {
	index = to_be_phased[i];
	tmpobshap.clear();
	tmpobsdep.clear();

	vector<pair<char *, char *> > seen_phased;
	uint tot_obshap = sink[index]->obshap.size() + src[index]->obshap.size();

	for (uint j = 0; j < sink[index]->obshap.size(); j++) {
	    for (uint k = 0; k < src[index]->obshap.size(); k++) {
		if (phased_haplotypes.count(string(sink[index]->obshap[j]) + string(src[index]->obshap[k])))
		    seen_phased.push_back(make_pair(sink[index]->obshap[j], src[index]->obshap[k]));
	    }
	}

	for (uint j = 0; j < seen_phased.size(); j++) {
	    for (uint k = j; k < seen_phased.size(); k++) {
		set<char *> incorporated_haplotypes;
		incorporated_haplotypes.insert(seen_phased[j].first);
		incorporated_haplotypes.insert(seen_phased[j].second);
		incorporated_haplotypes.insert(seen_phased[k].first);
		incorporated_haplotypes.insert(seen_phased[k].second);
		if (incorporated_haplotypes.size() == tot_obshap) {
		    tmpobshap.push_back(string(seen_phased[j].first) + string(seen_phased[j].second));
		    tmpobshap.push_back(string(seen_phased[k].first) + string(seen_phased[k].second));
		    //tmpobsdep.push_back((sink[index]->depth[j] + src[index]->depth[k]) / 2);
		}
	    }
	}

	sink[index]->depth.clear();
	for (uint j = 0; j < sink[index]->obshap.size(); j++)
	    delete [] sink[index]->obshap[j];
	sink[index]->obshap.clear();
	for (uint j = 0; j < tmpobshap.size(); j++) {
	    new_hap = new char[tmpobshap[j].length() + 1];
	    strcpy(new_hap, tmpobshap[j].c_str());
	    sink[index]->obshap.push_back(new_hap);
	    // sink[index]->depth.push_back(tmpobsdep[j]);
	}
    }

    //
    // 6. Merge model calls; Set the length; combine the two depth and lnl measures together.
    //
    string model_calls;
    char  *p;

    for (int i = 0; i < sample_cnt; i++) {
	if (sink[i] == NULL && src[i] == NULL)
	    continue;

	//
	// Merge the two strings of model calls together.
	// We need to check if the locus for this individual is shorter than the catalog
	// locus. If so, we need to expand out the model call array to be the proper length.
	//
	reverse_string(sink[i]->model);
	offset = 0;
	model_calls.clear();
	if (sink_locus_len > sink[i]->len) {
	    offset = sink_locus_len - sink[i]->len;
	    model_calls.assign(offset, 'N');
	}
	model_len = offset + sink[i]->len + src[i]->len - renz_olap[enz];
	model_calls.append(sink[i]->model);
	delete [] sink[i]->model;
	sink[i]->model = new char[model_len + 1];
	strcpy(sink[i]->model, model_calls.c_str());
	p  = sink[i]->model;
	p += offset + sink[i]->len - renz_olap[enz];
	strcpy(p, src[i]->model);
	
	sink[i]->len       = model_len;
	sink[i]->tot_depth = (sink[i]->tot_depth + src[i]->tot_depth) / 2;
	sink[i]->lnl       = (sink[i]->lnl + src[i]->lnl) / 2.0;

	//
	// Record which datum was merged into this one.
	//
	sink[i]->merge_partner = src[i]->id;
    }

    return 1;
}

int
datum_adjust_snp_positions(map<int, pair<merget, int> > &merge_map, 
			   CSLocus *loc, Datum *datum, 
			   map<int, SNPRes *> &snpres)
{
    //
    // We will start with the 'sink' locus, which was originally on the negative strand:
    //   1. If the locus was shorter than the catalog locus, pad the difference.
    //   2. Convert to positive strand: Reverse the order, complement the alleles,
    //      alter the internal column position.
    //
    SNP    *snp;
    SNPRes *snpr     = snpres[datum->id];
    int     index    = 0;
    int     stop_pos = renz_olap[enz] - 1;

    //
    // We know the catalog was padded since we already padded hte model call string
    // if it was necessary when originally merging.
    //
    while (datum->model[index] == 'N') {
	snp         = new SNP;
	snp->col    = index;
	snp->lratio = 0.0;
	snp->rank_1 = 'N';
	snp->type   = snp_type_unk;
	datum->snps.push_back(snp);
	index++;
    }

    for (int j = snpr->snps.size() - 1; j > stop_pos; j--) {
	snp         = new SNP;
	snp->col    = index;
	snp->lratio = snpr->snps[j]->lratio;
	snp->rank_1 = reverse(snpr->snps[j]->rank_1);
	snp->rank_2 = reverse(snpr->snps[j]->rank_2);
	snp->rank_3 = reverse(snpr->snps[j]->rank_3);
	snp->rank_4 = reverse(snpr->snps[j]->rank_4);
	datum->snps.push_back(snp);
	index++;
    }

    //
    // Now we fetch the former locus, the 'src', which was originally on the positive strand.
    // All we have to do is adjust the column position of each SNP.
    //
    snpr = snpres[datum->merge_partner];

    for (uint j = 0; j < snpres[datum->id]->snps.size(); j++) {
	snp         = new SNP;
	snp->col    = index;
	snp->lratio = snpr->snps[j]->lratio;
	snp->rank_1 = snpr->snps[j]->rank_1;
	snp->rank_2 = snpr->snps[j]->rank_2;
	snp->rank_3 = snpr->snps[j]->rank_3;
	snp->rank_4 = snpr->snps[j]->rank_4;
	datum->snps.push_back(snp);
	index++;
    }

    return 0;
}

int 
create_genotype_map(CSLocus *locus, PopMap<CSLocus> *pmap) 
{
    //
    // Create a genotype map. For any set of haplotypes, this routine will
    // assign each haplotype to a genotype, e.g. given the haplotypes 
    // 'AC' and 'GT' in the population, this routine will assign 'AC' == 'a' 
    // and 'GT' == 'b'. If an individual is homozygous for 'AC', they will be 
    // assigned an 'aa' genotype.
    //
    //cerr << "Creating genotype map for catalog ID " << locus->id  << ", marker: " << locus->marker << ".\n";

    char gtypes[26] ={'a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j',
		      'k', 'l', 'm', 'n', 'o', 'p', 'q', 'r', 's', 't',
		      'u', 'v', 'w', 'x', 'y', 'z'};

    Datum **d;
    map<string, int> haplotypes;
    map<string, int>::iterator k;
    vector<pair<string, int> > sorted_haplotypes;

    d = pmap->locus(locus->id);

    for (int i = 0; i < pmap->sample_cnt(); i++) {

	if (d[i] != NULL)
	    for (uint n = 0; n < d[i]->obshap.size(); n++)
		haplotypes[d[i]->obshap[n]]++;
    }

    //
    // Check that there are not more haplotypes than we have encodings.
    //
    if (haplotypes.size() > 26) return 0;

    // 
    // Sort the haplotypes map by value
    //
    for (k = haplotypes.begin(); k != haplotypes.end(); k++)
	sorted_haplotypes.push_back(*k);
    sort(sorted_haplotypes.begin(), sorted_haplotypes.end(), hap_compare);

    for (uint n = 0, index = 0; n < sorted_haplotypes.size() && index <= 26; n++, index++) {
	locus->gmap[sorted_haplotypes[n].first] = gtypes[index];
	//cerr << "GMAP: " << sorted_haplotypes[n].first << " == " << gtypes[index] << "\n";
    }

    return 0;
}

int call_population_genotypes(CSLocus *locus, 
			      PopMap<CSLocus> *pmap) {
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

 	string m = gtype.length() == 1 ? 
	    gtype + gtype : gtype;

	d[i]->gtype = new char[m.length() + 1];
	strcpy(d[i]->gtype, m.c_str());

	if (m != "-")
	    locus->gcnt++;

	//cerr << "Assigning datum, marker: " << locus->marker << ", string: " << m << ", haplotype: " << d[i]->obshap[0] << ", gtype: " << gtype << "\n";
     }

    return 0;
}

int tally_haplotype_freq(CSLocus *locus, PopMap<CSLocus> *pmap,
			 int &total, double &max, string &freq_str) {

    map<string, double> freq;
    Datum **d = pmap->locus(locus->id);

    total = 0;
    max   = 0;

    //cerr << "Examining marker: " << locus->id << "\n";

    for (int i = 0; i < pmap->sample_cnt(); i++) {
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

int write_genomic(map<int, CSLocus *> &catalog, PopMap<CSLocus> *pmap) {
    stringstream pop_name;
    pop_name << "batch_" << batch_id << ".genomic.tsv";

    string file = in_path + pop_name.str();

    ofstream fh(file.c_str(), ofstream::out);

    if (fh.fail()) {
        cerr << "Error opening genomic output file '" << file << "'\n";
	exit(1);
    }

    //
    // Count the number of markers that have enough samples to output.
    //
    map<int, CSLocus *>::iterator cit;
    CSLocus *loc;
    int num_loci = 0;

    for (cit = catalog.begin(); cit != catalog.end(); cit++) {
	loc = cit->second;

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
    map<string, vector<CSLocus *> >::iterator it;
    int  a, b;

    uint  rcnt = enz.length() ? renz_cnt[enz] : 0;
    uint  rlen = enz.length() ? renz_len[enz] : 0;
    char *p;

    for (it = pmap->ordered_loci.begin(); it != pmap->ordered_loci.end(); it++) {
	for (uint i = 0; i < it->second.size(); i++) {
	    loc = it->second[i];

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
		fh << loc->id << "\t" << loc->loc.chr << "\t" << loc->sort_bp(n);

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
calculate_haplotype_stats(vector<pair<int, string> > &files, map<int, pair<int, int> > &pop_indexes, 
			  map<int, CSLocus *> &catalog, PopMap<CSLocus> *pmap, PopSum<CSLocus> *psum) 
{
    map<string, vector<CSLocus *> >::iterator it;
    CSLocus  *loc;
    Datum   **d;
    LocStat  *l;

    //
    // Instantiate the kernel smoothing and bootstrap objects if requested.
    //
    KSmooth<LocStat>     *ks;
    OHaplotypes<LocStat> *ord;
    Bootstrap<LocStat>   *bs;
    if (kernel_smoothed && loci_ordered) {
	ks  = new KSmooth<LocStat>(2);
	ord = new OHaplotypes<LocStat>();
    }

    //
    // Open output file and print header.
    //
    stringstream pop_name;
    pop_name << "batch_" << batch_id << ".hapstats" << ".tsv";
    string file = in_path + pop_name.str();

    ofstream fh(file.c_str(), ofstream::out);
    if (fh.fail()) {
        cerr << "Error opening haplotype stats file '" << file << "'\n";
	exit(1);
    }
    fh.precision(fieldw);
    fh.setf(std::ios::fixed);

    map<int, pair<int, int> >::iterator pit;
    int start, end, pop_id;

    //
    // Write the population members.
    //
    for (pit = pop_indexes.begin(); pit != pop_indexes.end(); pit++) {
	start = pit->second.first;
	end   = pit->second.second;
	fh << "# " << pop_key[pit->first] << "\t";
	for (int i = start; i <= end; i++) {
	    fh << files[i].second;
	    if (i < end) fh << ",";
	}
	fh << "\n";
    }

    fh << "# Batch ID "    << "\t"
       << "Locus ID"       << "\t"
       << "Chr"            << "\t"
       << "BP"             << "\t"
       << "Pop ID"         << "\t"
       << "N"              << "\t"
       << "Haplotype Cnt"  << "\t"
       << "Gene Diversity" << "\t"
       << "Smoothed Gene Diversity"      << "\t"
       << "Smoothed Gene Diversity P-value"      << "\t"
       << "Haplotype Diversity"          << "\t"
       << "Smoothed Haplotype Diversity" << "\t"
       << "Smoothed Haplotype Diversity P-value" << "\t"
       << "Haplotypes"                   << "\n";

    //
    // Iterate over the members of each population.
    //
    for (pit = pop_indexes.begin(); pit != pop_indexes.end(); pit++) {
	start     = pit->second.first;
	end       = pit->second.second;
	pop_id    = pit->first;

    	cerr << "Generating haplotype-level summary statistics for population '" << pop_key[pop_id] << "'\n";
	map<string, vector<LocStat *> > genome_locstats;

	for (it = pmap->ordered_loci.begin(); it != pmap->ordered_loci.end(); it++) {

	    if (bootstrap_div)
		bs = new Bootstrap<LocStat>(2);

	    vector<LocStat *> &locstats = genome_locstats[it->first];
	    map<uint, uint>    locstats_key;
	    ord->order(locstats, locstats_key, it->second);

	    for (uint pos = 0; pos < it->second.size(); pos++) {
		loc = it->second[pos];
		d   = pmap->locus(loc->id);

		if (loc->snps.size() == 0) 
		    continue;

		// cerr << "Looking at locus " << loc->id << "\n";

		l = haplotype_diversity(start, end, d);

		if (l != NULL) {
		    l->loc_id = loc->id;
		    l->bp     = loc->sort_bp();
		    locstats[locstats_key[l->bp]] = l;
		}
	    }

	    if (kernel_smoothed && loci_ordered) {
		cerr << "    Generating kernel-smoothed statistics on chromosome " << it->first << "\n";
		ks->smooth(locstats);
	    }

	    if (bootstrap_div) 
		bs->add_data(locstats);
	}

	for (it = pmap->ordered_loci.begin(); it != pmap->ordered_loci.end(); it++) {
	    vector<LocStat *> &locstats = genome_locstats[it->first];

	    if (bootstrap_div)
		bs->execute(locstats);

	    //
	    // Write results.
	    //
	    for (uint k = 0; k < locstats.size(); k++) {
		l = locstats[k];
		if (l == NULL) continue;

		fh << batch_id         << "\t"
		   << l->loc_id        << "\t"
		   << it->first        << "\t"
		   << l->bp + 1        << "\t"
		   << pop_key[pop_id]  << "\t"
		   << (int) l->alleles << "\t"
		   << l->hap_cnt       << "\t"
		   << l->stat[0]       << "\t"
		   << l->smoothed[0]   << "\t"
		   << l->bs[0]         << "\t"
		   << l->stat[1]       << "\t"
		   << l->smoothed[1]   << "\t"
		   << l->bs[1]         << "\t"
		   << l->hap_str       << "\n";
	    }

	    for (uint k = 0; k < locstats.size(); k++) 
		delete locstats[k];
	}

	if (bootstrap_div) 
	    delete bs;
    }

    if (kernel_smoothed && loci_ordered) {
	delete ks;
	delete ord;
    }

    fh.close();

    return 0;
}

int
nuc_substitution_dist(map<string, int> &hap_index, double **hdists) 
{
    vector<string> haplotypes;
    map<string, int>::iterator it;
    uint i, j;

    for (it = hap_index.begin(); it != hap_index.end(); it++)
	haplotypes.push_back(it->first);

    const char *p, *q;
    double dist;

    for (i = 0; i < haplotypes.size(); i++) {
	for (j = i; j < haplotypes.size(); j++) {

	    dist = 0.0;
	    p    = haplotypes[i].c_str();
	    q    = haplotypes[j].c_str();

	    while (*p != '\0' && *q != '\0') {
		if (*p != *q) dist++;
		p++;
		q++;
	    }

	    hdists[i][j] = dist;
	    hdists[j][i] = dist;
	}
    }

    // //
    // // Print the distance matrix.
    // //
    // cerr << "  ";
    // for (hit = loc_hap_index.begin(); hit != loc_hap_index.end(); hit++)
    // 	cerr << "\t" << hit->first;
    // cerr << "\n";
    // for (hit = loc_hap_index.begin(); hit != loc_hap_index.end(); hit++) {
    // 	cerr << "  " << hit->first;
    // 	for (hit_2 = loc_hap_index.begin(); hit_2 != loc_hap_index.end(); hit_2++)
    // 	    cerr << "\t" << hdists[hit->second][hit_2->second];
    // 	cerr << "\n";
    // }
    // cerr << "\n";

    return 0;
}

int
nuc_substitution_identity(map<string, int> &hap_index, double **hdists) 
{
    vector<string> haplotypes;
    map<string, int>::iterator it;
    uint i, j;

    for (it = hap_index.begin(); it != hap_index.end(); it++)
	haplotypes.push_back(it->first);

    double dist;

    for (i = 0; i < haplotypes.size(); i++) {
	for (j = i; j < haplotypes.size(); j++) {

	    if (haplotypes[i] == haplotypes[j])
		dist = 0.0;
	    else
		dist = 1.0;

	    hdists[i][j] = dist;
	    hdists[j][i] = dist;
	}
    }

    return 0;
}

int
nuc_substitution_identity_max(map<string, int> &hap_index, double **hdists) 
{
    vector<string> haplotypes;
    map<string, int>::iterator it;
    uint i, j;

    for (it = hap_index.begin(); it != hap_index.end(); it++)
	haplotypes.push_back(it->first);

    for (i = 0; i < haplotypes.size(); i++) {
	for (j = i; j < haplotypes.size(); j++) {
	    hdists[i][j] = 1.0;
	    hdists[j][i] = 1.0;
	}
    }

    return 0;
}

int 
calculate_haplotype_divergence(vector<pair<int, string> > &files, 
			       map<int, pair<int, int> > &pop_indexes, 
			       map<int, vector<int> > &master_grp_members,
			       map<int, CSLocus *> &catalog, PopMap<CSLocus> *pmap, PopSum<CSLocus> *psum) 
{
    map<string, vector<CSLocus *> >::iterator it;

    if (bootstrap_phist)
	cerr << "Calculating halotype F statistics across all populations/groups and bootstrap resampling...\n";
    else
	cerr << "Calculating haplotype F statistics across all populations/groups...\n";

    //
    // Create a list of all the groups we have.
    //
    map<int, vector<int> >::iterator git;
    map<int, int> pop_grp_key;
    for (git = master_grp_members.begin(); git != master_grp_members.end(); git++)
	for (uint i = 0; i < git->second.size(); i++)
	    pop_grp_key[git->second[i]] = git->first;

    //
    // Create a list of all the populations we have.
    //
    vector<int> pop_ids;

    map<int, pair<int, int> >::iterator pit;
    for (pit = pop_indexes.begin(); pit != pop_indexes.end(); pit++)
	pop_ids.push_back(pit->first);

    //
    // Instantiate the kernel smoothing object and associated ordering object if requested.
    //
    KSmooth<HapStat>     *ks;
    OHaplotypes<HapStat> *ord;
    Bootstrap<HapStat>   *bs;
    if (kernel_smoothed && loci_ordered) {
	ks  = new KSmooth<HapStat>(5);
	ord = new OHaplotypes<HapStat>();
    }

    if (bootstrap_phist)
	bs = new Bootstrap<HapStat>(5);

    map<string, vector<HapStat *> > genome_hapstats;

    uint cnt = 0;
    for (it = pmap->ordered_loci.begin(); it != pmap->ordered_loci.end(); it++) {
	string chr = it->first;

	cerr << "  Generating haplotype F statistics for " << chr << "...";

	map<uint, uint>    hapstats_key;
	vector<HapStat *> &hapstats = genome_hapstats[chr];
	ord->order(hapstats, hapstats_key, it->second);

        #pragma omp parallel
	{ 
	    CSLocus  *loc;
	    LocSum  **s;
	    Datum   **d;
	    HapStat  *h;

            #pragma omp for schedule(dynamic, 1) reduction(+:cnt)
	    for (uint pos = 0; pos < it->second.size(); pos++) {
		loc = it->second[pos];
		s   = psum->locus(loc->id);
		d   = pmap->locus(loc->id);

		if (loc->snps.size() == 0)
		    continue;

		//
		// If this locus only appears in one population or there is only a single haplotype,
		// do not calculate haplotype F stats.
		//
		if (fixed_locus(pop_indexes, d, pop_ids))
		    continue;

		cnt++;
		// cerr << "Processing locus " << loc->id << "\n";

		h = haplotype_amova(pop_grp_key, pop_indexes, d, s, pop_ids);

		if (h != NULL) {
		    h->stat[4] = haplotype_d_est(pop_indexes, d, s, pop_ids);

		    h->loc_id = loc->id;
		    h->bp     = loc->sort_bp();
		    hapstats[hapstats_key[h->bp]] = h;
		}
	    }
	}

	if (bootstrap_phist)
	    bs->add_data(hapstats);

	cerr << "done.\n";

	//
	// Calculate kernel-smoothed Fst values.
	//
	if (kernel_smoothed && loci_ordered) {
	    cerr << "  Generating kernel-smoothed haplotype F statistics for " << it->first << "...";
	    ks->smooth(hapstats);
	    cerr << "done.\n";
	}
    }

    if (bootstrap_phist) {
	for (it = pmap->ordered_loci.begin(); it != pmap->ordered_loci.end(); it++)
	    bs->execute(genome_hapstats[it->first]);
    }

    cerr << "done.\n";

    if (kernel_smoothed && loci_ordered) {
	delete ks;
	delete ord;
    }

    if (bootstrap_phist)
	delete bs;

    cerr << "Writing haplotype F statistics... ";

    stringstream pop_name;
    pop_name << "batch_" << batch_id << ".phistats" << ".tsv";

    string file = in_path + pop_name.str();

    ofstream fh(file.c_str(), ofstream::out);
    if (fh.fail()) {
        cerr << "Error opening haplotype Phi_st file '" << file << "'\n";
	exit(1);
    }
    fh.precision(fieldw);
    fh.setf(std::ios::fixed);

    //
    // Write the population members.
    //
    int start, end;
    for (pit = pop_indexes.begin(); pit != pop_indexes.end(); pit++) {
	start = pit->second.first;
	end   = pit->second.second;
	fh << "# Population " << pop_key[pit->first] << "\t";
	for (int k = start; k <= end; k++) {
	    fh << files[k].second;
	    if (k < end) fh << ",";
	}
	fh << "\n";
    }

    //
    // Write the group members.
    //
    for (git = grp_members.begin(); git != grp_members.end(); git++) {
	end = git->second.size();
	fh << "# Group " << grp_key[git->first] << "\t";
	for (int k = 0; k < end; k++) {
	    fh << pop_key[git->second[k]];
	    if (k < end - 1) fh << ",";
	}
	fh << "\n";
    }

    fh << "# Batch ID " << "\t"
       << "Locus ID"    << "\t"
       << "Chr"         << "\t"
       << "BP"          << "\t"
       << "PopCnt"      << "\t";
    if (log_fst_comp)
	fh << "SSD(WP)"     << "\t"
	   << "SSD(AP/WG)"  << "\t"
	   << "SSD(AG)"     << "\t"
	   << "SSD(TOTAL)"  << "\t"
	   << "MSD(WP)"     << "\t" 
	   << "MSD(AP/WG)"  << "\t" 
	   << "MSD(AG)"     << "\t"
	   << "MSD(TOTAL)"  << "\t"
	   << "n"           << "\t"
	   << "n'"          << "\t"
	   << "n''"         << "\t"
	   << "Sigma2_a"    << "\t" 
	   << "Sigma2_b"    << "\t"
	   << "Sigma2_c"    << "\t"
	   << "Sigma_Total" << "\t";
    fh << "phi_st"          << "\t"
       << "Smoothed Phi_st" << "\t"
       << "Smoothed Phi_st P-value" << "\t"
       << "Phi_ct"          << "\t"
       << "Smoothed Phi_ct" << "\t"
       << "Smoothed Phi_ct P-value" << "\t"
       << "Phi_sc"          << "\t"
       << "Smoothed Phi_sc" << "\t"
       << "Smoothed Phi_sc P-value" << "\t"
       << "Fst'"            << "\t"
       << "Smoothed Fst'"   << "\t"
       << "Smoothed Fst' P-value"   << "\t"
       << "D_est"           << "\t"
       << "Smoothed D_est"  << "\t"
       << "Smoothed D_est P-value"  << "\n";

    for (it = pmap->ordered_loci.begin(); it != pmap->ordered_loci.end(); it++) {
	string chr = it->first;

	vector<HapStat *> &hapstats = genome_hapstats[chr];

	for (uint k = 0; k < hapstats.size(); k++) {
	    if (hapstats[k] == NULL) continue;

	    fh << batch_id            << "\t"
	       << hapstats[k]->loc_id << "\t"
	       << chr                 << "\t"
	       << hapstats[k]->bp     << "\t"
	       << hapstats[k]->popcnt << "\t";
	    if (log_fst_comp)
		fh << hapstats[k]->comp[0]  << "\t"
		   << hapstats[k]->comp[1]  << "\t"
		   << hapstats[k]->comp[2]  << "\t"
		   << hapstats[k]->comp[3]  << "\t"
		   << hapstats[k]->comp[4]  << "\t" 
		   << hapstats[k]->comp[5]  << "\t" 
		   << hapstats[k]->comp[6]  << "\t"
		   << hapstats[k]->comp[7]  << "\t"
		   << hapstats[k]->comp[8]  << "\t"
		   << hapstats[k]->comp[9]  << "\t"
		   << hapstats[k]->comp[10] << "\t"
		   << hapstats[k]->comp[11] << "\t" 
		   << hapstats[k]->comp[12] << "\t"
		   << hapstats[k]->comp[13] << "\t"
		   << hapstats[k]->comp[14] << "\t";
	    fh << hapstats[k]->stat[0]     << "\t"
	       << hapstats[k]->smoothed[0] << "\t"
	       << hapstats[k]->bs[0]       << "\t"
	       << hapstats[k]->stat[1]     << "\t"
	       << hapstats[k]->smoothed[1] << "\t"
	       << hapstats[k]->bs[1]       << "\t"
	       << hapstats[k]->stat[2]     << "\t"
	       << hapstats[k]->smoothed[2] << "\t"
	       << hapstats[k]->bs[2]       << "\t"
	       << hapstats[k]->stat[3]     << "\t"
	       << hapstats[k]->smoothed[3] << "\t"
	       << hapstats[k]->bs[3]       << "\t"
	       << hapstats[k]->stat[4]     << "\t"
	       << hapstats[k]->smoothed[4] << "\t"
	       << hapstats[k]->bs[4]       << "\n";

	    delete hapstats[k];
	}
    }

    fh.close();

    cerr << "wrote " << cnt << " loci to haplotype Phi_st file, '" << file << "'\n";

    return 0;
}

int 
calculate_haplotype_divergence_pairwise(vector<pair<int, string> > &files, 
					map<int, pair<int, int> > &pop_indexes, 
					map<int, vector<int> > &master_grp_members,
					map<int, CSLocus *> &catalog, PopMap<CSLocus> *pmap, PopSum<CSLocus> *psum) 
{
    map<string, vector<CSLocus *> >::iterator it;

    if (bootstrap_phist)
	cerr << "Calculating pairwise halotype F statistics and bootstrap resampling...\n";
    else
	cerr << "Calculating pairwise haplotype F statistics...\n";

    //
    // Assign all individuals to one group for the pairwise calculations.
    //
    map<int, vector<int> >::iterator git;
    map<int, int> pop_grp_key;
    for (git = master_grp_members.begin(); git != master_grp_members.end(); git++)
	for (uint i = 0; i < git->second.size(); i++)
	    pop_grp_key[git->second[i]] = 1;

    map<int, pair<int, int> >::iterator pit;
    vector<int> pop_ids;
    for (pit = pop_indexes.begin(); pit != pop_indexes.end(); pit++)
	pop_ids.push_back(pit->first);

    //
    // Instantiate the kernel smoothing object if requested.
    //
    KSmooth<HapStat>     *ks;
    OHaplotypes<HapStat> *ord;
    Bootstrap<HapStat>   *bs;
    if (kernel_smoothed && loci_ordered) {
	ks  = new KSmooth<HapStat>(5);
	ord = new OHaplotypes<HapStat>();
    }

    for (uint i = 0; i < pop_ids.size(); i++) {
	for (uint j = i + 1; j < pop_ids.size(); j++) {

	    if (bootstrap_phist)
		bs = new Bootstrap<HapStat>(5);

	    map<string, vector<HapStat *> > genome_hapstats;
	    vector<int> subpop_ids;

	    subpop_ids.push_back(pop_ids[i]);
	    subpop_ids.push_back(pop_ids[j]);

	    cerr << "  Processing populations '" << pop_key[pop_ids[i]] << "' and '" << pop_key[pop_ids[j]] << "'\n";

	    uint cnt = 0;
	    for (it = pmap->ordered_loci.begin(); it != pmap->ordered_loci.end(); it++) {
		string chr = it->first;

		cerr << "    Generating pairwise haplotype F statistics for " << chr << "...";

		map<uint, uint>    hapstats_key;
		vector<HapStat *> &hapstats = genome_hapstats[chr];
		ord->order(hapstats, hapstats_key, it->second);

                #pragma omp parallel
		{ 
		    CSLocus  *loc;
		    LocSum  **s;
		    Datum   **d;
		    HapStat  *h;

                    #pragma omp for schedule(dynamic, 1) reduction(+:cnt)
		    for (uint pos = 0; pos < it->second.size(); pos++) {
			loc = it->second[pos];
			s   = psum->locus(loc->id);
			d   = pmap->locus(loc->id);

			if (loc->snps.size() == 0)
			    continue;

			//
			// If this locus only appears in one population or there is only a single haplotype,
			// do not calculate haplotype F stats.
			//
			if (fixed_locus(pop_indexes, d, subpop_ids))
			    continue;

			cnt++;
			// cerr << "Processing locus " << loc->id << "\n";

			h = haplotype_amova(pop_grp_key, pop_indexes, d, s, subpop_ids);

			if (h != NULL) {
			    h->stat[4] = haplotype_d_est(pop_indexes, d, s, subpop_ids);

			    h->loc_id = loc->id;
			    h->bp     = loc->sort_bp();
			    hapstats[hapstats_key[h->bp]] = h;
			}
		    }
		}

		if (bootstrap_phist)
		    bs->add_data(hapstats);

		cerr << "done.\n";

		//
		// Calculate kernel-smoothed Fst values.
		//
		if (kernel_smoothed && loci_ordered) {
		    cerr << "    Generating kernel-smoothed Phi_st for " << it->first << "...";
		    ks->smooth(hapstats);
		    cerr << "done.\n";
		}
	    }

	    if (bootstrap_phist) {
		for (it = pmap->ordered_loci.begin(); it != pmap->ordered_loci.end(); it++)
		    bs->execute(genome_hapstats[it->first]);
	    }

	    cerr << "done.\n";

	    if (bootstrap_phist)
		delete bs;

	    cerr << "Writing haplotype F statistics... ";

	    stringstream pop_name;
	    pop_name << "batch_" << batch_id << ".phistats_" << pop_key[pop_ids[i]] << "-" << pop_key[pop_ids[j]] << ".tsv";

	    string file = in_path + pop_name.str();

	    ofstream fh(file.c_str(), ofstream::out);
	    if (fh.fail()) {
		cerr << "Error opening haplotype Phi_st file '" << file << "'\n";
		exit(1);
	    }
	    fh.precision(fieldw);
	    fh.setf(std::ios::fixed);

	    //
	    // Write the population members.
	    //
	    int start, end;
	    for (uint k = 0; k < subpop_ids.size(); k++) {
		start = pop_indexes[subpop_ids[k]].first;
		end   = pop_indexes[subpop_ids[k]].second;
		fh << "# Population " << pop_key[subpop_ids[k]] << "\t";
		for (int n = start; n <= end; n++) {
		    fh << files[n].second;
		    if (n < end) fh << ",";
		}
		fh << "\n";
	    }

	    fh << "# Batch ID " << "\t"
	       << "Locus ID"    << "\t"
	       << "Pop 1 ID"   << "\t"
	       << "Pop 2 ID"   << "\t"
	       << "Chr"         << "\t"
	       << "BP"          << "\t";
	    if (log_fst_comp)
		fh << "SSD(WP)"     << "\t"
		   << "SSD(AP/WG)"  << "\t"
		   << "SSD(AG)"     << "\t"
		   << "SSD(TOTAL)"  << "\t"
		   << "MSD(WP)"     << "\t" 
		   << "MSD(AP/WG)"  << "\t" 
		   << "MSD(AG)"     << "\t"
		   << "MSD(TOTAL)"  << "\t"
		   << "n"           << "\t"
		   << "n'"          << "\t"
		   << "n''"         << "\t"
		   << "Sigma2_a"    << "\t" 
		   << "Sigma2_b"    << "\t"
		   << "Sigma2_c"    << "\t"
		   << "Sigma_Total" << "\t";
	    fh << "phi_st"          << "\t"
	       << "Smoothed Phi_st" << "\t"
	       << "Smoothed Phi_st P-value" << "\t"
	       << "Fst'"            << "\t"
	       << "Smoothed Fst'"   << "\t"
	       << "Smoothed Fst' P-value"   << "\t"
	       << "D_est"          << "\t"
	       << "Smoothed D_est" << "\t"
	       << "Smoothed D_est P-value" << "\n";

	    for (it = pmap->ordered_loci.begin(); it != pmap->ordered_loci.end(); it++) {
		string chr = it->first;

		vector<HapStat *> &hapstats = genome_hapstats[chr];

		for (uint k = 0; k < hapstats.size(); k++) {
		    if (hapstats[k] == NULL) continue;

		    fh << batch_id            << "\t"
		       << hapstats[k]->loc_id << "\t"
		       << pop_key[pop_ids[i]] << "\t"
		       << pop_key[pop_ids[j]] << "\t"
		       << chr                 << "\t"
		       << hapstats[k]->bp     << "\t";
		    if (log_fst_comp)
			fh << hapstats[k]->comp[0]  << "\t"
			   << hapstats[k]->comp[1]  << "\t"
			   << hapstats[k]->comp[2]  << "\t"
			   << hapstats[k]->comp[3]  << "\t"
			   << hapstats[k]->comp[4]  << "\t" 
			   << hapstats[k]->comp[5]  << "\t" 
			   << hapstats[k]->comp[6]  << "\t"
			   << hapstats[k]->comp[7]  << "\t"
			   << hapstats[k]->comp[8]  << "\t"
			   << hapstats[k]->comp[9]  << "\t"
			   << hapstats[k]->comp[10] << "\t"
			   << hapstats[k]->comp[11] << "\t" 
			   << hapstats[k]->comp[12] << "\t"
			   << hapstats[k]->comp[13] << "\t"
			   << hapstats[k]->comp[14] << "\t";
		    fh << hapstats[k]->stat[0]     << "\t"
		       << hapstats[k]->smoothed[0] << "\t"
		       << hapstats[k]->bs[0]       << "\t"
		       << hapstats[k]->stat[3]     << "\t"
		       << hapstats[k]->smoothed[3] << "\t"
		       << hapstats[k]->bs[3]       << "\t"
		       << hapstats[k]->stat[4]     << "\t"
		       << hapstats[k]->smoothed[4] << "\t"
		       << hapstats[k]->bs[4]       << "\n";

		    delete hapstats[k];
		}
	    }

	    fh.close();

	    cerr << "wrote " << cnt << " loci to pairwise haplotype file, '" << file << "'\n";
	}
    }

    if (kernel_smoothed && loci_ordered) {
	delete ks;
	delete ord;
    }

    return 0;
}

bool
fixed_locus(map<int, pair<int, int> > &pop_indexes, Datum **d, vector<int> &pop_ids)
{
    set<string>               loc_haplotypes;
    map<int, vector<string> > pop_haplotypes;
    int start, end, pop_id;
    int pop_cnt = pop_ids.size();

    for (int p = 0; p < pop_cnt; p++) {
	start  = pop_indexes[pop_ids[p]].first;
	end    = pop_indexes[pop_ids[p]].second;
	pop_id = pop_ids[p];

	for (int i = start; i <= end; i++) {
	    if (d[i] == NULL) continue;

	    if (d[i]->obshap.size() > 2) { 
		continue;

	    } else if (d[i]->obshap.size() == 1) {
		if (!uncalled_haplotype(d[i]->obshap[0])) {
		    loc_haplotypes.insert(d[i]->obshap[0]);
		    pop_haplotypes[pop_id].push_back(d[i]->obshap[0]);
		    pop_haplotypes[pop_id].push_back(d[i]->obshap[0]);
		}
	    } else {
		for (uint j = 0; j < d[i]->obshap.size(); j++) {
		    if (!uncalled_haplotype(d[i]->obshap[0])) {
		        loc_haplotypes.insert(d[i]->obshap[j]);
			pop_haplotypes[pop_id].push_back(d[i]->obshap[j]);
		    }
		}
	    }
	}
    }

    uint valid_pops = 0;

    for (int p = 0; p < pop_cnt; p++) {
	pop_id = pop_ids[p];

	if (pop_haplotypes[pop_id].size() > 0) 
	    valid_pops++;
    }

    //
    // Check that more than one population has data for this locus.
    //
    if (valid_pops <= 1) 
	return true;

    //
    // Check that there is more than one haplotype at this locus.
    //
    if (loc_haplotypes.size() == 1)
	return true;

    return false;
}

inline bool
uncalled_haplotype(const char *haplotype)
{
    for (const char *p = haplotype; *p != '\0'; p++)
	if (*p == 'N' || *p == 'n')
	    return true;
    return false;
}

inline double
count_haplotypes_at_locus(int start, int end, Datum **d, map<string, double> &hap_cnts)
{
    double n = 0.0;

    for (int i = start; i <= end; i++) {
	if (d[i] == NULL) continue;

	if (d[i]->obshap.size() > 2) { 
	    continue;

	} else if (d[i]->obshap.size() == 1) {
	    if(!uncalled_haplotype(d[i]->obshap[0])) {
		n += 2;
		hap_cnts[d[i]->obshap[0]] += 2;
	    }
	} else {
	    for (uint j = 0; j < d[i]->obshap.size(); j++) {
		if(!uncalled_haplotype(d[i]->obshap[0])) {
		    n++;
		    hap_cnts[d[i]->obshap[j]]++;
		}
	    }
	}
    }
    
    return n;
}

LocStat *
haplotype_diversity(int start, int end, Datum **d)
{
    map<string, double>::iterator hit;
    vector<string>      haplotypes;
    map<string, double> hap_freq;
    map<string, int>    hap_index;
    double   n              = 0.0;
    double   gene_diversity = 0.0;
    double   hapl_diversity = 0.0;
    LocStat *lstat;

    //
    // Tabulate the haplotypes in this population.
    //
    n = count_haplotypes_at_locus(start, end, d, hap_freq);

    // cerr << "  " << n << " total haplotypes observed.\n";

    //
    // If this haplotype is fixed, don't calculate any statistics.
    //
    if (n == 0) 
	return NULL;

    lstat = new LocStat;

    //
    // Store a summary of the haplotype counts to output below.
    //
    stringstream sstr;
    for (hit = hap_freq.begin(); hit != hap_freq.end(); hit++)
	sstr << hit->first << ":" << hit->second  << ";";
    lstat->hap_str = sstr.str().substr(0, sstr.str().length() - 1);

    //
    // Determine an ordering for the haplotypes. Convert haplotype counts into frequencies.
    //
    uint k = 0;
    for (hit = hap_freq.begin(); hit != hap_freq.end(); hit++) {
	hap_index[hit->first] = k;
	haplotypes.push_back(hit->first);
	k++;

	// cerr << "  Haplotype '" << hit->first << "' occured " << hit->second  << " times; ";

	hit->second = hit->second / n;

	// cerr << " frequency of " << hit->second << "%\n";
    }
    //
    // Initialize a two-dimensional array to hold distances between haplotyes.
    //
    double **hdists = new double *[hap_index.size()];
    for (k = 0; k < hap_index.size(); k++) {
	hdists[k] = new double[hap_index.size()];
	memset(hdists[k], 0, hap_index.size());
    }

    //
    // Calculate the distances between haplotypes.
    //
    nuc_substitution_dist(hap_index, hdists);

    //
    // Calculate haplotype diversity, Pi.
    //
    for (uint i = 0; i < haplotypes.size(); i++) {
	for (uint j = 0; j < haplotypes.size(); j++) {
	    hapl_diversity += 
		hap_freq[haplotypes[i]] * 
		hap_freq[haplotypes[j]] * 
		hdists[hap_index[haplotypes[i]]][hap_index[haplotypes[j]]];
	}
    }
    hapl_diversity = (n / (n-1)) * hapl_diversity;

    //
    // Calculate gene diversity. 
    //
    for (uint i = 0; i < haplotypes.size(); i++) {
	gene_diversity += hap_freq[haplotypes[i]] * hap_freq[haplotypes[i]];
    }
    gene_diversity = (n / (n - 1)) * (1 - gene_diversity);

    lstat->alleles = n;
    lstat->stat[0] = gene_diversity;
    lstat->stat[1] = hapl_diversity;
    lstat->hap_cnt = haplotypes.size();

    // cerr << "  Population " << pop_id << " has haplotype diversity (pi) of " << s[pop_index]->pi << "\n";

    for (k = 0; k < hap_index.size(); k++)
	delete hdists[k];
    delete hdists;

    return lstat;
}

HapStat *
haplotype_amova(map<int, int> &pop_grp_key, map<int, pair<int, int> > &pop_indexes, 
		Datum **d, LocSum **s, vector<int> &pop_ids)
{
    map<string, int>          loc_hap_index;
    vector<string>            loc_haplotypes;
    map<int, vector<string> > pop_haplotypes;
    map<int, vector<int> >    grp_members;
    vector<int>               grps;

    map<string, int>::iterator hit, hit_2;
    map<int, pair<int, int> >::iterator pit;
    int start, end, pop_id, pop_id_1;

    HapStat  *h;
    int pop_cnt = pop_ids.size();

    //
    // Tabulate the occurences of haplotypes at this locus.
    //
    for (int p = 0; p < pop_cnt; p++) {
	start  = pop_indexes[pop_ids[p]].first;
	end    = pop_indexes[pop_ids[p]].second;
	pop_id = pop_ids[p];

	for (int i = start; i <= end; i++) {
	    if (d[i] == NULL) continue;

	    if (d[i]->obshap.size() > 2) { 
		continue;

	    } else if (d[i]->obshap.size() == 1) {
		if(!uncalled_haplotype(d[i]->obshap[0])) {
		    loc_hap_index[d[i]->obshap[0]]++;
		    loc_haplotypes.push_back(d[i]->obshap[0]);
		    loc_haplotypes.push_back(d[i]->obshap[0]);
		    pop_haplotypes[pop_id].push_back(d[i]->obshap[0]);
		    pop_haplotypes[pop_id].push_back(d[i]->obshap[0]);
		}
	    } else {
		for (uint j = 0; j < d[i]->obshap.size(); j++) {
		    if(!uncalled_haplotype(d[i]->obshap[0])) {
			loc_hap_index[d[i]->obshap[j]]++;
			loc_haplotypes.push_back(d[i]->obshap[j]);
			pop_haplotypes[pop_id].push_back(d[i]->obshap[j]);
		    }
		}
	    }
	}
    }

    //
    // What is the total number of populations that had valid haplotypes.
    //
    double valid_pop_cnt = 0.0;
    for (int p = 0; p < pop_cnt; p++) {
	pop_id = pop_ids[p];
	if (pop_haplotypes[pop_id].size() > 0)
	    valid_pop_cnt++;
    }

    //
    // If we filtered a population out at this locus make sure that we still have at least one
    // representative present in each group.
    //
    set<int> uniq_grps;
    for (pit = pop_indexes.begin(); pit != pop_indexes.end(); pit++) {
	pop_id = pit->first;

	if (pop_haplotypes.count(pop_id) > 0) {
	    uniq_grps.insert(pop_grp_key[pop_id]);
	    grp_members[pop_grp_key[pop_id]].push_back(pop_id);
	}
    }
    set<int>::iterator uit;
    for (uit = uniq_grps.begin(); uit != uniq_grps.end(); uit++)
	grps.push_back(*uit);

    if (grps.size() == 0)
	return NULL;

    // cerr << "Groups: ";
    // for (uint i = 0; i < grps.size(); i++)
    //     cerr << grps[i] << ", ";
    // cerr << "\n";
    // for (git = grp_members.begin(); git != grp_members.end(); git++) {
    //     cerr << "Group " << git->first << ": ";
    //     for (uint i = 0; i < git->second.size(); i++)
    // 	cerr << git->second[i] << ", ";
    //     cerr << "\n";
    // }

    //
    // Determine an ordering for the haplotypes.
    //
    uint m = 0;
    for (hit = loc_hap_index.begin(); hit != loc_hap_index.end(); hit++) {
	loc_hap_index[hit->first] = m;
	m++;
    }

    //
    // Initialize a two-dimensional array to hold distances between haplotyes.
    //
    double **hdists     = new double *[loc_hap_index.size()];
    double **hdists_max = new double *[loc_hap_index.size()];
    for (uint k = 0; k < loc_hap_index.size(); k++) {
	hdists[k] = new double[loc_hap_index.size()];
	memset(hdists[k], 0, loc_hap_index.size());
	hdists_max[k] = new double[loc_hap_index.size()];
	memset(hdists_max[k], 0, loc_hap_index.size());
    }

    //
    // Calculate the distances between haplotypes.
    //
    nuc_substitution_dist(loc_hap_index, hdists);

    //
    // Calculate the sum of squared distances in each subset: total, within populations, across populations 
    // and withing groups, and across groups.
    //
    double ssd_total = amova_ssd_total(loc_haplotypes, loc_hap_index, hdists);
    double ssd_wp    = amova_ssd_wp(grps, grp_members, loc_hap_index, pop_haplotypes, hdists);
    double ssd_ap_wg = amova_ssd_ap_wg(grps, grp_members, loc_hap_index, pop_haplotypes, hdists, hdists);
    double ssd_ag    = grps.size() > 1 ? amova_ssd_ag(grps, grp_members, loc_hap_index, pop_haplotypes, hdists, ssd_total) : 0.0;

    //
    // Calculate n
    //
    double n        = 0.0;
    double n_1      = 0.0;
    double n_2      = 0.0;
    double s_g      = 0.0;
    double tot_cnt  = 0.0;
    double grp_cnt  = 0.0;
    double num_grps = grps.size();
    double a        = 0.0;
    double b        = 0.0;

    for (uint g = 0; g < num_grps; g++) {
	for (uint r = 0; r < grp_members[grps[g]].size(); r++) {
	    pop_id_1 = grp_members[grps[g]][r];
	    tot_cnt += (double) pop_haplotypes[pop_id_1].size();
	}
    }
    for (uint g = 0; g < num_grps; g++) {
	grp_cnt = 0.0;
	for (uint r = 0; r < grp_members[grps[g]].size(); r++) {
	    pop_id_1 = grp_members[grps[g]][r];
	    grp_cnt += (double) pop_haplotypes[pop_id_1].size();
	}

	a = 0.0;
	for (uint r = 0; r < grp_members[grps[g]].size(); r++) {
	    pop_id_1 = grp_members[grps[g]][r];
	    a += (double) (pop_haplotypes[pop_id_1].size() * pop_haplotypes[pop_id_1].size()) / grp_cnt;
	}
	s_g += a;
    }
    n = (tot_cnt - s_g) / (double) (valid_pop_cnt - num_grps);

    // cerr << "  n: "<< n << "\n";

    if (num_grps > 1) {
	//
	// Calculate n'
	//
	a = 0.0;
	for (uint g = 0; g < num_grps; g++) {
	    for (uint r = 0; r < grp_members[grps[g]].size(); r++) {
		pop_id_1 = grp_members[grps[g]][r];
		a += ((double) (pop_haplotypes[pop_id_1].size() * pop_haplotypes[pop_id_1].size()) / tot_cnt);
	    }
	}
	n_1 = (s_g - a) / (double) (num_grps - 1.0);

	// cerr << "  n': "<< n_1 << "\n";

	//
	// Calculate n''
	//
	for (uint g = 0; g < num_grps; g++) {
	    a = 0.0;
	    for (uint r = 0; r < grp_members[grps[g]].size(); r++) {
		pop_id_1 = grp_members[grps[g]][r];
		a += pop_haplotypes[pop_id_1].size();
	    }
	    b += ((a * a) / tot_cnt);
	}
	n_2 = (tot_cnt - b) / (double) (num_grps - 1);

	// cerr << "  n'': "<< n_2 << "\n";
    }

    //
    // Calculate the mean square deviations, equal to SSD divided by degrees of freedom.
    //
    double msd_ag    = num_grps > 1 ? ssd_ag / (double) (num_grps - 1) : 0.0;
    double msd_ap_wg = ssd_ap_wg / ((double) (valid_pop_cnt - num_grps));
    double msd_wp    = ssd_wp    / ((double) (loc_haplotypes.size() - valid_pop_cnt));
    double msd_total = ssd_total / ((double) (loc_haplotypes.size() - 1));

    double sigma_c     = msd_wp;
    double sigma_b     = n > 0 ? (msd_ap_wg - sigma_c) / n : 0.0;
    double sigma_a     = 0.0;

    if (grps.size() > 1)
	sigma_a = (msd_ag - sigma_c - (n_1 * sigma_b)) / n_2;

    // Arlequin seems to sum the variance components instead of independently calculating sigma_total: MSD(total) = SSD(total)/degrees.of.freedom
    double sigma_total = sigma_a + sigma_b + sigma_c; // msd_total; 

    double phi_st = 0.0;
    double phi_ct = 0.0;
    double phi_sc = 0.0;

    if (grps.size() > 1) {
	phi_st = sigma_total > 0.0 ? (sigma_a + sigma_b) / sigma_total : 0.0;
	phi_ct = sigma_total > 0.0 ?  sigma_a / sigma_total : 0.0;
	phi_sc = (sigma_a + sigma_b) > 0.0 ?  sigma_b / (sigma_b + sigma_c) : 0.0;
    } else {
	phi_st = sigma_total > 0.0 ? sigma_b / sigma_total : 0.0;
    }

    // cerr << "  MSD(AG): " << msd_ag  << "; MSD(AP/WG): " << msd_ap_wg << "; MSD(WP): " << msd_wp  << "; MSD(TOTAL): "  << msd_total   << "\n"
    //      << "  Sigma_a: " << sigma_a << "; Sigma_b: "    << sigma_b   << "; Sigma_c: " << sigma_c << "; Sigma_Total: " << sigma_total << "\n"
    //      << "  Phi_st: "  << phi_st  << "; Phi_ct: "     << phi_ct    << "; Phi_sc: "  << phi_sc  << "\n";


    //
    // Calculate Fst' = Fst / Fst_max
    //
    // First calculate Fst.
    //
    // To calculate Fst instead of Phi_st, we need to reset our distance matrix to return 1 if haplotypes are different, 0 otherwise.
    //
    nuc_substitution_identity(loc_hap_index, hdists);
    ssd_wp    = amova_ssd_wp(grps, grp_members, loc_hap_index, pop_haplotypes, hdists);
    ssd_ap_wg = amova_ssd_ap_wg(grps, grp_members, loc_hap_index, pop_haplotypes, hdists, hdists);
    //
    // Calculate the mean square deviations, equal to SSD divided by degrees of freedom.
    //
    msd_ap_wg   = ssd_ap_wg / ((double) (valid_pop_cnt - num_grps));
    msd_wp      = ssd_wp    / ((double) (loc_haplotypes.size() - valid_pop_cnt));
    sigma_c     = msd_wp;
    sigma_b     = n > 0 ? (msd_ap_wg - sigma_c) / n : 0.0;
    sigma_total = sigma_b + sigma_c;

    double fst = sigma_total > 0.0 ? sigma_b / sigma_total : 0.0;

    //
    // Now calculate Fst_max.
    //
    // Reset our distance matrix to give maximum possible distance between haplotypes
    // and recalculate sum of squared deviations across groups.
    //
    nuc_substitution_identity_max(loc_hap_index, hdists_max);
    ssd_ap_wg = amova_ssd_ap_wg(grps, grp_members, loc_hap_index, pop_haplotypes, hdists, hdists_max);

    //
    // Recalculate the mean square deviations, given maximum divergence between populations.
    //
    msd_ap_wg = ssd_ap_wg / ((double) (valid_pop_cnt - num_grps));
    sigma_b   = n > 0 ? (msd_ap_wg - sigma_c) / n : 0.0;

    double fst_max = sigma_total > 0.0 ? sigma_b / sigma_total : 0.0;
    double fst_1   = fst_max > 0.0     ? fst / fst_max         : 0.0;

    //
    // Cache the results so we can print them in order below, once the parallel code has executed.
    //
    h = new HapStat;
    h->alleles = tot_cnt;
    h->popcnt  = valid_pop_cnt;

    if (log_fst_comp) {
	h->comp = new double[15];
	h->comp[0]  = ssd_wp;
	h->comp[1]  = ssd_ap_wg;
	h->comp[2]  = ssd_ag;
	h->comp[3]  = ssd_total;
	h->comp[4]  = msd_wp;
	h->comp[5]  = msd_ap_wg;
	h->comp[6]  = msd_ag;
	h->comp[7]  = msd_total;
	h->comp[8]  = n;
	h->comp[9]  = n_1;
	h->comp[10] = n_2;
	h->comp[11] = sigma_a;
	h->comp[12] = sigma_b;
	h->comp[13] = sigma_c;
	h->comp[14] = sigma_total;
    }

    h->stat[0] = phi_st;
    h->stat[1] = phi_ct;
    h->stat[2] = phi_sc;
    h->stat[3] = fst_1;

    for (uint k = 0; k < loc_hap_index.size(); k++) {
	delete [] hdists[k];
	delete [] hdists_max[k];
    }
    delete [] hdists;
    delete [] hdists_max;

    return h;
}

double
amova_ssd_total(vector<string> &loc_haplotypes, map<string, int> &loc_hap_index, double **hdists)
{
    //
    // Calculate sum of squared deviations for the total sample, SSD(Total)
    //
    double ssd_total = 0.0;

    for (uint j = 0; j < loc_haplotypes.size(); j++) {
	for (uint k = 0; k < loc_haplotypes.size(); k++) {
	    ssd_total += hdists[loc_hap_index[loc_haplotypes[j]]][loc_hap_index[loc_haplotypes[k]]];
	    // cerr << j << "\t" 
	    // 	 << k << "\t" 
	    // 	 << loc_haplotypes[j] << "\t" 
	    // 	 << loc_haplotypes[k] << "\t" 
	    // 	 << hdists[loc_hap_index[loc_haplotypes[j]]][loc_hap_index[loc_haplotypes[k]]] << "\n";
	}
    }
    ssd_total = (1.0 / (double) (2*loc_haplotypes.size())) * ssd_total;
    // cerr << "  ssd_total: "<< ssd_total << "\n";

    return ssd_total;
}

double
amova_ssd_wp(vector<int> &grps, map<int, vector<int> > &grp_members, 
	     map<string, int> &loc_hap_index, map<int, vector<string> > &pop_haplotypes, 
	     double **hdists)
{
    //
    // Calculate the sum of squared deviations within populations, SSD(WP)
    //
    double ssd_wp = 0.0;
    double ssd    = 0.0;
    int    pop_id;

    for (uint g = 0; g < grps.size(); g++) {
	for (uint i = 0; i < grp_members[grps[g]].size(); i++) {
	    pop_id = grp_members[grps[g]][i];
	    ssd = 0.0;

	    for (uint j = 0; j < pop_haplotypes[pop_id].size(); j++) {
		for (uint k = 0; k < pop_haplotypes[pop_id].size(); k++) {
		    ssd += hdists[loc_hap_index[pop_haplotypes[pop_id][j]]][loc_hap_index[pop_haplotypes[pop_id][k]]];
		    // cerr << pop_id << "\t"
		    // 	 << j << "\t"
		    // 	 << k << "\t" 
		    // 	 << loc_haplotypes[j] << "\t" 
		    // 	 << loc_haplotypes[k] << "\t" 
		    // 	 << hdists[loc_hap_index[loc_haplotypes[j]]][loc_hap_index[loc_haplotypes[k]]] << "\n";
		}
	    }

	    if (pop_haplotypes[pop_id].size() > 0)
		ssd_wp += (1.0 / (double) (2*pop_haplotypes[pop_id].size())) * ssd;
	}
    }
    // cerr << "  ssd_wp: "<< ssd_wp << "\n";

    return ssd_wp;
}

double
amova_ssd_ap_wg(vector<int> &grps, map<int, vector<int> > &grp_members, 
		map<string, int> &loc_hap_index, map<int, vector<string> > &pop_haplotypes, 
		double **hdists_1, double **hdists_2)
{
    //
    // Calculate the sum of squared deviations across populations and within groups, SSD(AP/WG)
    //
    double ssd_ap_wg = 0.0;
    double ssd       = 0.0;
    double ssd_1     = 0.0;
    double ssd_2     = 0.0;
    double den       = 0.0;
    int    pop_id, pop_id_1, pop_id_2;

    for (uint g = 0; g < grps.size(); g++) {

	ssd_1 = 0.0;
	for (uint r = 0; r < grp_members[grps[g]].size(); r++) {
	    pop_id_1 = grp_members[grps[g]][r];

	    for (uint j = 0; j < pop_haplotypes[pop_id_1].size(); j++) {

		for (uint s = 0; s < grp_members[grps[g]].size(); s++) {
		    pop_id_2 = grp_members[grps[g]][s];

		    for (uint k = 0; k < pop_haplotypes[pop_id_2].size(); k++) {
			if (pop_id_1 == pop_id_2)
			    ssd_1 += hdists_1[loc_hap_index[pop_haplotypes[pop_id_1][j]]][loc_hap_index[pop_haplotypes[pop_id_2][k]]];
			else
			    ssd_1 += hdists_2[loc_hap_index[pop_haplotypes[pop_id_1][j]]][loc_hap_index[pop_haplotypes[pop_id_2][k]]];
		    }
		}
	    }
	}

	den = 0.0;
	for (uint r = 0; r < grp_members[grps[g]].size(); r++) {
	    pop_id_1 = grp_members[grps[g]][r];
	    den += 2 * pop_haplotypes[pop_id_1].size();
	}

	ssd_1 = ssd_1 / den;

	ssd_2 = 0.0;
	for (uint r = 0; r < grp_members[grps[g]].size(); r++) {
	    pop_id = grp_members[grps[g]][r];
	    ssd = 0.0;

	    for (uint j = 0; j < pop_haplotypes[pop_id].size(); j++) {
		for (uint k = 0; k < pop_haplotypes[pop_id].size(); k++) {
		    ssd += hdists_1[loc_hap_index[pop_haplotypes[pop_id][j]]][loc_hap_index[pop_haplotypes[pop_id][k]]];
		}
	    }

	    if (pop_haplotypes[pop_id].size() > 0)
		ssd_2 += (1.0 / (double) (2*pop_haplotypes[pop_id].size())) * ssd;
	}

	ssd_ap_wg += ssd_1 - ssd_2;
    }
    // cerr << "  ssd_ap_wg: "<< ssd_ap_wg << "\n";

    return ssd_ap_wg;
}

double
amova_ssd_ag(vector<int> &grps, map<int, vector<int> > &grp_members, 
	     map<string, int> &loc_hap_index, map<int, vector<string> > &pop_haplotypes, 
	     double **hdists, double ssd_total)
{
    //
    // Calculate the sum of squared deviations across groups, SSD(AG)
    //
    int    pop_id_1, pop_id_2;
    double ssd_ag = 0.0;
    double ssd    = 0.0;
    double ssd_1  = 0.0;
    double den    = 0.0;

    for (uint g = 0; g < grps.size(); g++) {
	ssd_1 = 0.0;

	for (uint r = 0; r < grp_members[grps[g]].size(); r++) {
	    pop_id_1 = grp_members[grps[g]][r];

	    for (uint j = 0; j < pop_haplotypes[pop_id_1].size(); j++) {

		for (uint s = 0; s < grp_members[grps[g]].size(); s++) {
		    pop_id_2 = grp_members[grps[g]][s];

		    for (uint k = 0; k < pop_haplotypes[pop_id_2].size(); k++) {
			ssd_1 += hdists[loc_hap_index[pop_haplotypes[pop_id_1][j]]][loc_hap_index[pop_haplotypes[pop_id_2][k]]];
		    }
		}
	    }
	}

	den = 0.0;
	for (uint r = 0; r < grp_members[grps[g]].size(); r++) {
	    pop_id_1 = grp_members[grps[g]][r];
	    den += 2 * pop_haplotypes[pop_id_1].size();
	}

	ssd += ssd_1 / den;
    }

    ssd_ag = ssd_total - ssd;

    // cerr << "  ssd_ag: "<< ssd_ag << "\n";

    return ssd_ag;
}

double
haplotype_d_est(map<int, pair<int, int> > &pop_indexes, Datum **d, LocSum **s, vector<int> &pop_ids)
{
    //
    // Calculate D_est, fixation index, as described by 
    //   Bird, et al., 2011, Detecting and measuring genetic differentiation
    //     +-Equation 11
    // and
    //   Jost, 2008, GST and its relatives do not measure differentiation, Molecular Ecology
    //     +- Equation 13, D_est_chao
    //
    map<string, double>            loc_haplotypes;
    map<int, map<string, double> > pop_haplotypes;
    map<int, double>               pop_totals;

    map<string, double>::iterator it;
    int start, end, pop_id;

    uint pop_cnt = pop_ids.size();

    //
    // Tabulate the occurences of haplotypes at this locus.
    //
    for (uint p = 0; p < pop_cnt; p++) {
    	start  = pop_indexes[pop_ids[p]].first;
    	end    = pop_indexes[pop_ids[p]].second;
    	pop_id = pop_ids[p];

    	for (int i = start; i <= end; i++) {
    	    if (d[i] == NULL) continue;

    	    if (d[i]->obshap.size() > 2) { 
    		continue;

    	    } else if (d[i]->obshap.size() == 1) {
    		loc_haplotypes[d[i]->obshap[0]]         += 2;
    		pop_haplotypes[pop_id][d[i]->obshap[0]] += 2;

    	    } else {
    		for (uint j = 0; j < d[i]->obshap.size(); j++) {
		    loc_haplotypes[d[i]->obshap[j]]++;
		    pop_haplotypes[pop_id][d[i]->obshap[j]]++;
    		}
    	    }
    	}

	for (it = pop_haplotypes[pop_id].begin(); it != pop_haplotypes[pop_id].end(); it++)
	    pop_totals[pop_id] += it->second;
    }

    double x = 0.0;

    for (it = loc_haplotypes.begin(); it != loc_haplotypes.end(); it++) {

	double freq_sum_sq = 0.0;
	double freq_sq_sum = 0.0;
	for (uint p = 0; p < pop_cnt; p++) {
	    pop_id = pop_ids[p];
	    freq_sum_sq += (pop_haplotypes[pop_id][it->first] / pop_totals[pop_id]);
	    freq_sq_sum += pow((pop_haplotypes[pop_id][it->first] / pop_totals[pop_id]), 2);
	}
	freq_sum_sq = pow(freq_sum_sq, 2);

	x += (freq_sum_sq - freq_sq_sum) / (pop_cnt - 1);
    }

    double y = 0.0;

    for (it = loc_haplotypes.begin(); it != loc_haplotypes.end(); it++) {
	for (uint p = 0; p < pop_cnt; p++) {
	    pop_id = pop_ids[p];

	    y += (pop_haplotypes[pop_id][it->first] * (pop_haplotypes[pop_id][it->first] - 1)) /
		(pop_totals[pop_id] * (pop_totals[pop_id] - 1));
	}
    }

    double d_est = 1.0 - (x / y);

    return d_est;
}

int 
calculate_summary_stats(vector<pair<int, string> > &files, map<int, pair<int, int> > &pop_indexes, 
			map<int, CSLocus *> &catalog, PopMap<CSLocus> *pmap, PopSum<CSLocus> *psum) 
{
    map<string, vector<CSLocus *> >::iterator it;
    CSLocus  *loc;
    LocSum  **s;
    LocTally *t;
    int       len;
    int       pop_cnt = psum->pop_cnt();

    //
    // Calculate the means for each summary statistic.
    //
    int    *private_cnt;
    double *num_indv_mean, *p_mean, *obs_het_mean, *obs_hom_mean, *exp_het_mean, *exp_hom_mean, *pi_mean, *fis_mean;
    double *num_indv_var,  *p_var,  *obs_het_var,  *obs_hom_var,  *exp_het_var,  *exp_hom_var,  *pi_var,  *fis_var;
    double *num_indv_mean_all, *p_mean_all, *obs_het_mean_all, *obs_hom_mean_all, *exp_het_mean_all, *exp_hom_mean_all, *pi_mean_all, *fis_mean_all;
    double *num_indv_var_all,  *p_var_all,  *obs_het_var_all,  *obs_hom_var_all,  *exp_het_var_all,  *exp_hom_var_all,  *pi_var_all,  *fis_var_all;
    double *n, *n_all, *var_sites;
    private_cnt       = new int[pop_cnt];
    n                 = new double[pop_cnt];
    var_sites         = new double[pop_cnt];
    num_indv_mean     = new double[pop_cnt];
    num_indv_var      = new double[pop_cnt];
    p_mean            = new double[pop_cnt];
    p_var             = new double[pop_cnt];
    obs_het_mean      = new double[pop_cnt];
    obs_het_var       = new double[pop_cnt];
    obs_hom_mean      = new double[pop_cnt];
    obs_hom_var       = new double[pop_cnt];
    exp_het_mean      = new double[pop_cnt];
    exp_het_var       = new double[pop_cnt];
    exp_hom_mean      = new double[pop_cnt];
    exp_hom_var       = new double[pop_cnt];
    pi_mean           = new double[pop_cnt];
    pi_var            = new double[pop_cnt];
    fis_mean          = new double[pop_cnt];
    fis_var           = new double[pop_cnt];

    n_all             = new double[pop_cnt];
    num_indv_mean_all = new double[pop_cnt];
    num_indv_var_all  = new double[pop_cnt];
    p_mean_all        = new double[pop_cnt];
    p_var_all         = new double[pop_cnt];
    obs_het_mean_all  = new double[pop_cnt];
    obs_het_var_all   = new double[pop_cnt];
    obs_hom_mean_all  = new double[pop_cnt];
    obs_hom_var_all   = new double[pop_cnt];
    exp_het_mean_all  = new double[pop_cnt];
    exp_het_var_all   = new double[pop_cnt];
    exp_hom_mean_all  = new double[pop_cnt];
    exp_hom_var_all   = new double[pop_cnt];
    pi_mean_all       = new double[pop_cnt];
    pi_var_all        = new double[pop_cnt];
    fis_mean_all      = new double[pop_cnt];
    fis_var_all       = new double[pop_cnt];

    for (int j = 0; j < pop_cnt; j++) {
	private_cnt[j]   = 0;
	n[j]             = 0.0;
	var_sites[j]     = 0.0;
	num_indv_mean[j] = 0.0;
	num_indv_var[j]  = 0.0;
	p_mean[j]        = 0.0;
	p_var[j]         = 0.0;
	obs_het_mean[j]  = 0.0;
	obs_het_var[j]   = 0.0;
	obs_hom_mean[j]  = 0.0;
	obs_hom_var[j]   = 0.0;
	exp_het_mean[j]  = 0.0;
	exp_het_var[j]   = 0.0;
	exp_hom_mean[j]  = 0.0;
	exp_hom_var[j]   = 0.0;
	pi_mean[j]       = 0.0;
	pi_var[j]        = 0.0;
	fis_mean[j]      = 0.0;
	fis_var[j]       = 0.0;

	n_all[j]             = 0.0;
	num_indv_mean_all[j] = 0.0;
	num_indv_var_all[j]  = 0.0;
	p_mean_all[j]        = 0.0;
	p_var_all[j]         = 0.0;
	obs_het_mean_all[j]  = 0.0;
	obs_het_var_all[j]   = 0.0;
	obs_hom_mean_all[j]  = 0.0;
	obs_hom_var_all[j]   = 0.0;
	exp_het_mean_all[j]  = 0.0;
	exp_het_var_all[j]   = 0.0;
	exp_hom_mean_all[j]  = 0.0;
	exp_hom_var_all[j]   = 0.0;
	pi_mean_all[j]       = 0.0;
	pi_var_all[j]        = 0.0;
	fis_mean_all[j]      = 0.0;
	fis_var_all[j]       = 0.0;
    }

    for (it = pmap->ordered_loci.begin(); it != pmap->ordered_loci.end(); it++) {
	for (uint pos = 0; pos < it->second.size(); pos++) {
	    loc = it->second[pos];

	    s = psum->locus(loc->id);
	    t = psum->locus_tally(loc->id);
	    len = strlen(loc->con);

	    for (int i = 0; i < len; i++) {
		//
		// Compile private alleles
		//
		if (t->nucs[i].priv_allele >= 0)
		    private_cnt[t->nucs[i].priv_allele]++;

		if (t->nucs[i].allele_cnt == 2) {

		    for (int j = 0; j < pop_cnt; j++) {

			if (s[j]->nucs[i].num_indv == 0) continue;

			n[j]++;

			if (s[j]->nucs[i].pi > 0) var_sites[j]++;

			num_indv_mean[j] += s[j]->nucs[i].num_indv;
			p_mean[j]        += s[j]->nucs[i].p;
			obs_het_mean[j]  += s[j]->nucs[i].obs_het;
			obs_hom_mean[j]  += s[j]->nucs[i].obs_hom;
			exp_het_mean[j]  += s[j]->nucs[i].exp_het;
			exp_hom_mean[j]  += s[j]->nucs[i].exp_hom;
			pi_mean[j]       += s[j]->nucs[i].stat[0];
			fis_mean[j]      += s[j]->nucs[i].stat[1] != -7.0 ? s[j]->nucs[i].stat[1] : 0.0;

			n_all[j]++;
			num_indv_mean_all[j] += s[j]->nucs[i].num_indv;
			p_mean_all[j]        += s[j]->nucs[i].p;
			obs_het_mean_all[j]  += s[j]->nucs[i].obs_het;
			obs_hom_mean_all[j]  += s[j]->nucs[i].obs_hom;
			exp_het_mean_all[j]  += s[j]->nucs[i].exp_het;
			exp_hom_mean_all[j]  += s[j]->nucs[i].exp_hom;
			pi_mean_all[j]       += s[j]->nucs[i].stat[0];
			fis_mean_all[j]      += s[j]->nucs[i].stat[1] != -7.0 ? s[j]->nucs[i].stat[1] : 0.0;
		    }

		} else if (t->nucs[i].allele_cnt == 1) {
		    for (int j = 0; j < pop_cnt; j++) {
			if (s[j]->nucs[i].num_indv == 0) continue;

			n_all[j]++;
			num_indv_mean_all[j] += s[j]->nucs[i].num_indv;
			p_mean_all[j]        += s[j]->nucs[i].p;
			obs_het_mean_all[j]  += s[j]->nucs[i].obs_het;
			obs_hom_mean_all[j]  += s[j]->nucs[i].obs_hom;
			exp_het_mean_all[j]  += s[j]->nucs[i].exp_het;
			exp_hom_mean_all[j]  += s[j]->nucs[i].exp_hom;
			pi_mean_all[j]       += s[j]->nucs[i].stat[0];
			fis_mean_all[j]      += s[j]->nucs[i].stat[1] != -7.0 ? s[j]->nucs[i].stat[1] : 0.0;
		    }
		}
	    }
	}
    }

    for (int j = 0; j < pop_cnt; j++) {
	num_indv_mean[j] = num_indv_mean[j] / n[j];
	p_mean[j]        = p_mean[j]        / n[j];
	obs_het_mean[j]  = obs_het_mean[j]  / n[j];
	obs_hom_mean[j]  = obs_hom_mean[j]  / n[j];
	exp_het_mean[j]  = exp_het_mean[j]  / n[j];
	exp_hom_mean[j]  = exp_hom_mean[j]  / n[j];
	pi_mean[j]       = pi_mean[j]       / n[j];
	fis_mean[j]      = fis_mean[j]      / n[j];

	num_indv_mean_all[j] = num_indv_mean_all[j] / n_all[j];
	p_mean_all[j]        = p_mean_all[j]        / n_all[j];
	obs_het_mean_all[j]  = obs_het_mean_all[j]  / n_all[j];
	obs_hom_mean_all[j]  = obs_hom_mean_all[j]  / n_all[j];
	exp_het_mean_all[j]  = exp_het_mean_all[j]  / n_all[j];
	exp_hom_mean_all[j]  = exp_hom_mean_all[j]  / n_all[j];
	pi_mean_all[j]       = pi_mean_all[j]       / n_all[j];
	fis_mean_all[j]      = fis_mean_all[j]      / n_all[j];
    }

    stringstream pop_name;
    pop_name << "batch_" << batch_id << ".sumstats" << ".tsv";

    string file = in_path + pop_name.str();

    ofstream fh(file.c_str(), ofstream::out);
    if (fh.fail()) {
        cerr << "Error opening sumstats file '" << file << "'\n";
	exit(1);
    }
    fh.precision(fieldw);
    fh.setf(std::ios::fixed);

    double p_freq;
    int start, end;
    //
    // Write the population members.
    //
    map<int, pair<int, int> >::iterator pit;
    for (pit = pop_indexes.begin(); pit != pop_indexes.end(); pit++) {
	start = pit->second.first;
	end   = pit->second.second;
	fh << "# " << pit->first << "\t";
	for (int i = start; i <= end; i++) {
	    fh << files[i].second;
	    if (i < end) fh << ",";
	}
	fh << "\n";
    }

    cerr << "Writing " << catalog.size() << " loci to summary statistics file, '" << file << "'\n";

    fh << "# Batch ID " << "\t"
       << "Locus ID" << "\t"
       << "Chr"      << "\t"
       << "BP"       << "\t"
       << "Col"      << "\t"
       << "Pop ID"   << "\t"
       << "P Nuc"    << "\t"
       << "Q Nuc"    << "\t"
       << "N"        << "\t"
       << "P"        << "\t"
       << "Obs Het"  << "\t"
       << "Obs Hom"  << "\t"
       << "Exp Het"  << "\t"
       << "Exp Hom"  << "\t"
       << "Pi"       << "\t"
       << "Smoothed Pi"  << "\t"
       << "Smoothed Pi P-value"  << "\t"
       << "Fis"          << "\t"
       << "Smoothed Fis" << "\t"
       << "Smoothed Fis P-value" << "\t"
       << "Private"      << "\n";

    for (it = pmap->ordered_loci.begin(); it != pmap->ordered_loci.end(); it++) {
	for (uint pos = 0; pos < it->second.size(); pos++) {
	    loc = it->second[pos];

	    s = psum->locus(loc->id);
	    t = psum->locus_tally(loc->id);
	    len = strlen(loc->con);

	    for (int i = 0; i < len; i++) {

		// 
		// If this site is fixed in all populations, DON'T output it. If it is variable,
		// or fixed within populations but variable among, DO output it.
		//
		if (t->nucs[i].allele_cnt == 2) {

		    for (int j = 0; j < pop_cnt; j++) {

			if (s[j]->nucs[i].num_indv == 0) continue;

			fh << batch_id << "\t"
			   << loc->id << "\t"
			   << loc->loc.chr << "\t"
			   << loc->sort_bp(i) + 1 << "\t"
			   << i << "\t"
			   << pop_key[psum->rev_pop_index(j)] << "\t";

			//
			// Output the p and q alleles in the same order in each population.
			//
			if (t->nucs[i].p_allele == s[j]->nucs[i].p_nuc) {
			    if (s[j]->nucs[i].q_nuc == 0)
				fh << s[j]->nucs[i].p_nuc << "\t" << "-";
			    else
				fh << s[j]->nucs[i].p_nuc << "\t" << s[j]->nucs[i].q_nuc;
			    p_freq = s[j]->nucs[i].p;

			} else {
			    if (s[j]->nucs[i].q_nuc == 0)
				fh << "-\t" << s[j]->nucs[i].p_nuc;
			    else
				fh << s[j]->nucs[i].q_nuc << "\t" << s[j]->nucs[i].p_nuc;
			    p_freq = 1 - s[j]->nucs[i].p;
			}

			fh << "\t" << (int) s[j]->nucs[i].num_indv << "\t"
			   << std::setprecision(8)      << p_freq << "\t"
			   << std::setprecision(fieldw) << s[j]->nucs[i].obs_het << "\t"
			   << s[j]->nucs[i].obs_hom   << "\t"
			   << s[j]->nucs[i].exp_het   << "\t"
			   << s[j]->nucs[i].exp_hom   << "\t"
			   << s[j]->nucs[i].stat[0]   << "\t" // Pi
			   << s[j]->nucs[i].smoothed[0] << "\t"  // Smoothed Pi
			   << s[j]->nucs[i].bs[0]       << "\t"  // Pi bootstrapped p-value
			   << (s[j]->nucs[i].stat[1] == -7.0 ? 0.0 : s[j]->nucs[i].stat[1]) << "\t"  // Fis
			   << s[j]->nucs[i].smoothed[1] << "\t"  // Smoothed Fis
			   << s[j]->nucs[i].bs[1]       << "\t"; // Fis bootstrapped p-value.
			(t->nucs[i].priv_allele == j) ? fh << "1\n" : fh << "0\n";

			//
			// Tabulate the residuals to calculate the variance.
			//
			num_indv_var[j] += pow((s[j]->nucs[i].num_indv - num_indv_mean[j]), 2);
			p_var[j]        += pow((s[j]->nucs[i].p        - p_mean[j]),        2);
			obs_het_var[j]  += pow((s[j]->nucs[i].obs_het  - obs_het_mean[j]),  2);
			obs_hom_var[j]  += pow((s[j]->nucs[i].obs_hom  - obs_hom_mean[j]),  2);
			exp_het_var[j]  += pow((s[j]->nucs[i].exp_het  - exp_het_mean[j]),  2);
			exp_hom_var[j]  += pow((s[j]->nucs[i].exp_hom  - exp_hom_mean[j]),  2);
			pi_var[j]       += pow((s[j]->nucs[i].stat[0]  - pi_mean[j]),       2);
			fis_var[j]      += pow((s[j]->nucs[i].stat[1]  - fis_mean[j]),      2);

			num_indv_var_all[j] += pow((s[j]->nucs[i].num_indv - num_indv_mean_all[j]), 2);
			p_var_all[j]        += pow((s[j]->nucs[i].p        - p_mean_all[j]),        2);
			obs_het_var_all[j]  += pow((s[j]->nucs[i].obs_het  - obs_het_mean_all[j]),  2);
			obs_hom_var_all[j]  += pow((s[j]->nucs[i].obs_hom  - obs_hom_mean_all[j]),  2);
			exp_het_var_all[j]  += pow((s[j]->nucs[i].exp_het  - exp_het_mean_all[j]),  2);
			exp_hom_var_all[j]  += pow((s[j]->nucs[i].exp_hom  - exp_hom_mean_all[j]),  2);
			pi_var_all[j]       += pow((s[j]->nucs[i].stat[0]  - pi_mean_all[j]),       2);
			fis_var_all[j]      += pow((s[j]->nucs[i].stat[1]  - fis_mean_all[j]),      2);
		    }
		} else if (t->nucs[i].allele_cnt == 1) {
		    for (int j = 0; j < pop_cnt; j++) {
			if (s[j]->nucs[i].num_indv == 0) continue;

			num_indv_var_all[j] += pow((s[j]->nucs[i].num_indv - num_indv_mean_all[j]), 2);
			p_var_all[j]        += pow((s[j]->nucs[i].p        - p_mean_all[j]),        2);
			obs_het_var_all[j]  += pow((s[j]->nucs[i].obs_het  - obs_het_mean_all[j]),  2);
			obs_hom_var_all[j]  += pow((s[j]->nucs[i].obs_hom  - obs_hom_mean_all[j]),  2);
			exp_het_var_all[j]  += pow((s[j]->nucs[i].exp_het  - exp_het_mean_all[j]),  2);
			exp_hom_var_all[j]  += pow((s[j]->nucs[i].exp_hom  - exp_hom_mean_all[j]),  2);
			pi_var_all[j]       += pow((s[j]->nucs[i].stat[0]  - pi_mean_all[j]),       2);
			fis_var_all[j]      += pow((s[j]->nucs[i].stat[1]  - fis_mean_all[j]),      2);
		    }
		}
	    }
	}
    }

    //
    // Calculate the variance.
    //
    for (int j = 0; j < pop_cnt; j++) {
	num_indv_var[j] = num_indv_var[j] / (n[j] - 1);
	p_var[j]        = p_var[j] / (n[j] - 1);
	obs_het_var[j]  = obs_het_var[j] / (n[j] - 1);
	obs_hom_var[j]  = obs_hom_var[j] / (n[j] - 1);
	exp_het_var[j]  = exp_het_var[j] / (n[j] - 1);
	exp_hom_var[j]  = exp_hom_var[j] / (n[j] - 1);
	pi_var[j]       = pi_var[j] / (n[j] - 1);
	fis_var[j]      = fis_var[j] / (n[j] - 1);

	num_indv_var_all[j] = num_indv_var_all[j] / (n_all[j] - 1);
	p_var_all[j]        = p_var_all[j] / (n_all[j] - 1);
	obs_het_var_all[j]  = obs_het_var_all[j] / (n_all[j] - 1);
	obs_hom_var_all[j]  = obs_hom_var_all[j] / (n_all[j] - 1);
	exp_het_var_all[j]  = exp_het_var_all[j] / (n_all[j] - 1);
	exp_hom_var_all[j]  = exp_hom_var_all[j] / (n_all[j] - 1);
	pi_var_all[j]       = pi_var_all[j] / (n_all[j] - 1);
	fis_var_all[j]      = fis_var_all[j] / (n_all[j] - 1);
    }

    fh.close();

    pop_name.str("");
    pop_name << "batch_" << batch_id << ".sumstats_summary" << ".tsv";

    file = in_path + pop_name.str();

    fh.open(file.c_str(), ofstream::out);

    if (fh.fail()) {
        cerr << "Error opening sumstats summary file '" << file << "'\n";
	exit(1);
    }

    //
    // Write out summary statistics of the summary statistics.
    //
    fh << "# Variant positions\n"
       << "# Pop ID\t"
       << "Private\t"
       << "Num Indv\t"
       << "Var\t"
       << "StdErr\t"
       << "P\t"
       << "Var\t"
       << "StdErr\t"
       << "Obs Het\t"
       << "Var\t"
       << "StdErr\t"
       << "Obs Hom\t"
       << "Var\t"
       << "StdErr\t"
       << "Exp Het\t"
       << "Var\t"
       << "StdErr\t"
       << "Exp Hom\t"
       << "Var\t"
       << "StdErr\t"
       << "Pi\t"
       << "Var\t"
       << "StdErr\t"
       << "Fis\t"
       << "Var\t"
       << "StdErr\n";

    double *sq_n     = new double[pop_cnt];
    double *sq_n_all = new double[pop_cnt];

    for (int j = 0; j < pop_cnt; j++) {
	sq_n[j]     = sqrt(n[j]);
	sq_n_all[j] = sqrt(n_all[j]);
    }

    for (int j = 0; j < pop_cnt; j++)
	fh << pop_key[psum->rev_pop_index(j)] << "\t" 
	   << private_cnt[j]         << "\t"
	   << num_indv_mean[j]       << "\t"
	   << num_indv_var[j]        << "\t"
	   << sqrt(num_indv_var[j]) / sq_n[j] << "\t"
	   << p_mean[j]              << "\t"
	   << p_var[j]               << "\t"
	   << sqrt(p_var[j]) / sq_n[j] << "\t"
	   << obs_het_mean[j]        << "\t"
	   << obs_het_var[j]         << "\t"
	   << sqrt(obs_het_var[j]) / sq_n[j] << "\t"
	   << obs_hom_mean[j]        << "\t"
	   << obs_hom_var[j]         << "\t"
	   << sqrt(obs_hom_var[j]) / sq_n[j] << "\t"
	   << exp_het_mean[j]        << "\t"
	   << exp_het_var[j]         << "\t"
	   << sqrt(exp_het_var[j]) / sq_n[j] << "\t"
	   << exp_hom_mean[j]        << "\t"
	   << exp_hom_var[j]         << "\t"
	   << sqrt(exp_hom_var[j]) / sq_n[j] << "\t"
	   << pi_mean[j]             << "\t"
	   << pi_var[j]              << "\t"
	   << sqrt(pi_var[j]) / sq_n[j] << "\t"
	   << fis_mean[j]            << "\t"
	   << fis_var[j]             << "\t"
	   << sqrt(num_indv_var[j]) / sq_n[j] << "\n";

    fh << "# All positions (variant and fixed)\n"
       << "# Pop ID\t"
       << "Private\t"
       << "Sites\t"
       << "Variant Sites\t"
       << "Polymorphic Sites\t"
       << "% Polymorphic Loci\t"
       << "Num Indv\t"
       << "Var\t"
       << "StdErr\t"
       << "P\t"
       << "Var\t"
       << "StdErr\t"
       << "Obs Het\t"
       << "Var\t"
       << "StdErr\t"
       << "Obs Hom\t"
       << "Var\t"
       << "StdErr\t"
       << "Exp Het\t"
       << "Var\t"
       << "StdErr\t"
       << "Exp Hom\t"
       << "Var\t"
       << "StdErr\t"
       << "Pi\t"
       << "Var\t"
       << "StdErr\t"
       << "Fis\t"
       << "Var\t"
       << "StdErr\n";

    for (int j = 0; j < pop_cnt; j++) {
	fh << pop_key[psum->rev_pop_index(j)] << "\t" 
	   << private_cnt[j]             << "\t"
	   << n_all[j]                   << "\t"
	   << n[j]                       << "\t"
	   << var_sites[j]               << "\t"
	   << var_sites[j] / n_all[j] * 100 << "\t"
	   << num_indv_mean_all[j]       << "\t"
	   << num_indv_var_all[j]        << "\t"
	   << sqrt(num_indv_var_all[j]) / sq_n_all[j] << "\t"
	   << p_mean_all[j]              << "\t"
	   << p_var_all[j]               << "\t"
	   << sqrt(p_var_all[j]) / sq_n_all[j] << "\t"
	   << obs_het_mean_all[j]        << "\t"
	   << obs_het_var_all[j]         << "\t"
	   << sqrt(obs_het_var_all[j]) / sq_n_all[j] << "\t"
	   << obs_hom_mean_all[j]        << "\t"
	   << obs_hom_var_all[j]         << "\t"
	   << sqrt(obs_hom_var_all[j]) / sq_n_all[j] << "\t"
	   << exp_het_mean_all[j]        << "\t"
	   << exp_het_var_all[j]         << "\t"
	   << sqrt(exp_het_var_all[j]) / sq_n_all[j] << "\t"
	   << exp_hom_mean_all[j]        << "\t"
	   << exp_hom_var_all[j]         << "\t"
	   << sqrt(exp_hom_var_all[j]) / sq_n_all[j] << "\t"
	   << pi_mean_all[j]             << "\t"
	   << pi_var_all[j]              << "\t"
	   << sqrt(pi_var_all[j]) / sq_n_all[j] << "\t"
	   << fis_mean_all[j]            << "\t"
	   << fis_var_all[j]             << "\t"
	   << sqrt(num_indv_var_all[j]) / sq_n_all[j] << "\n";
    }

    delete [] private_cnt;
    delete [] n;
    delete [] var_sites;
    delete [] sq_n;
    delete [] num_indv_mean;
    delete [] num_indv_var;
    delete [] p_mean;
    delete [] p_var;
    delete [] obs_het_mean;
    delete [] obs_het_var;
    delete [] obs_hom_mean;
    delete [] obs_hom_var;
    delete [] exp_het_mean;
    delete [] exp_het_var;
    delete [] exp_hom_mean;
    delete [] exp_hom_var;
    delete [] pi_mean;
    delete [] pi_var;
    delete [] fis_mean;
    delete [] fis_var;

    delete [] n_all;
    delete [] sq_n_all;
    delete [] num_indv_mean_all;
    delete [] num_indv_var_all;
    delete [] p_mean_all;
    delete [] p_var_all;
    delete [] obs_het_mean_all;
    delete [] obs_het_var_all;
    delete [] obs_hom_mean_all;
    delete [] obs_hom_var_all;
    delete [] exp_het_mean_all;
    delete [] exp_het_var_all;
    delete [] exp_hom_mean_all;
    delete [] exp_hom_var_all;
    delete [] pi_mean_all;
    delete [] pi_var_all;
    delete [] fis_mean_all;
    delete [] fis_var_all;

    fh.close();

    return 0;
}

int 
write_fst_stats(vector<pair<int, string> > &files, map<int, pair<int, int> > &pop_indexes, 
		map<int, CSLocus *> &catalog, PopMap<CSLocus> *pmap, PopSum<CSLocus> *psum, ofstream &log_fh) 
{
    //
    // We want to iterate over each pair of populations and calculate Fst at each 
    // nucleotide of each locus.
    //
    vector<double> means;
    vector<int>    pops;
    map<int, pair<int, int> >::iterator pit;
    for (pit = pop_indexes.begin(); pit != pop_indexes.end(); pit++)
	pops.push_back(pit->first);

    if (pops.size() == 1) return 0;

    //
    // Instantiate the kernel smoothing object if requested.
    //
    OPopPair<PopPair>  *ord = new OPopPair<PopPair>(psum, log_fh);
    KSmooth<PopPair>   *ks;
    Bootstrap<PopPair> *bs;
    if (kernel_smoothed && loci_ordered)
	ks  = new KSmooth<PopPair>(2);

    for (uint i = 0; i < pops.size(); i++) {

	for (uint j = i + 1; j < pops.size(); j++) {
	    int pop_1 = pops[i];
	    int pop_2 = pops[j];

	    double sum = 0.0;
	    double cnt = 0.0;

	    stringstream pop_name;
	    pop_name << "batch_" << batch_id << ".fst_" << pop_key[pop_1] << "-" << pop_key[pop_2] << ".tsv";

	    string   file = in_path + pop_name.str();
	    ofstream fh(file.c_str(), ofstream::out);
	    if (fh.fail()) {
		cerr << "Error opening Fst output file '" << file << "'\n";
		exit(1);
	    }
	    fh.precision(fieldw);
	    fh.setf(std::ios::fixed);

	    cerr << "Calculating Fst for populations '" << pop_key[pop_1] << "' and '" << pop_key[pop_2] << "' and writing it to file, '" << file << "'\n";

	    fh << "# Batch ID" << "\t"
	       << "Locus ID"   << "\t"
	       << "Pop 1 ID"   << "\t"
	       << "Pop 2 ID"   << "\t"
	       << "Chr"        << "\t"
	       << "BP"         << "\t"
	       << "Column"     << "\t"
	       << "Overall Pi" << "\t"
	       << "Fst"        << "\t"
	       << "Fisher's P" << "\t"
	       << "Odds Ratio" << "\t"
	       << "CI Low"     << "\t"
	       << "CI High"    << "\t"
	       << "LOD"        << "\t"
	       << "Corrected Fst" << "\t"
	       << "Smoothed Fst"  << "\t"
	       << "AMOVA Fst" << "\t"
	       << "Corrected AMOVA Fst" << "\t"
	       << "Smoothed AMOVA Fst" << "\t"
	       << "Smoothed AMOVA Fst P-value" << "\t"
	       << "Window SNP Count";

	    //
	    // If requested, log Fst component calculations to a file.
	    //
	    if (log_fst_comp) {
		fh << "\t"
		   << "n_1" << "\t"
		   << "n_2" << "\t"
		   << "tot_alleles" << "\t"
		   << "p_1" << "\t"
		   << "q_1" << "\t"
		   << "p_2" << "\t"
		   << "q_2" << "\t"
		   << "pi_1" << "\t"
		   << "pi_2" << "\t"
		   << "pi_all" << "\t"
		   << "bcoeff_1" << "\t"
		   << "bcoeff_2" << "\t"
		   << "binomial_fst" << "\t"
		   << "p_1_freq" << "\t"
		   << "q_1_freq" << "\t"
		   << "p_2_freq" << "\t"
		   << "q_2_freq" << "\t"
		   << "p_avg_cor" << "\t"
		   << "n_avg_cor" << "\t"
		   << "amova_fst" << "\n";
	    } else {
		fh << "\n";
	    }

	    if (bootstrap_fst)
		bs = new Bootstrap<PopPair>(2);

	    map<string, vector<CSLocus *> >::iterator it;
	    map<string, vector<PopPair *> > genome_pairs;
	    // int snp_dist[max_snp_dist] = {0};

	    for (it = pmap->ordered_loci.begin(); it != pmap->ordered_loci.end(); it++) {
		string chr = it->first;

		map<uint, uint>    pairs_key;
		vector<PopPair *> &pairs = genome_pairs[chr];

		//
		// Order loci between the two populations and calculate Fst
		//
		ord->order(pairs, pairs_key, it->second, pop_1, pop_2);

		//
		// Apply user-selected correction to the Fst values.
		//
		double correction;
		switch(fst_correction) {
		case p_value:
		    for (uint i = 0; i < pairs.size(); i++) {
			if (pairs[i] != NULL) {
			    pairs[i]->stat[0] = pairs[i]->fet_p < p_value_cutoff ? pairs[i]->fst : 0;
			    pairs[i]->stat[1] = pairs[i]->fet_p < p_value_cutoff ? pairs[i]->amova_fst : 0;
			}
		    }
		    break;
		case bonferroni_win:
		    correct_fst_bonferroni_win(pairs);
		    break;
		case bonferroni_gen:
		    correction = p_value_cutoff / catalog.size();
		    for (uint i = 0; i < pairs.size(); i++) {
			if (pairs[i] != NULL) {
			    pairs[i]->stat[0] = pairs[i]->fet_p < correction ? pairs[i]->fst : 0;
			    pairs[i]->stat[1] = pairs[i]->fet_p < correction ? pairs[i]->amova_fst : 0;
			}
		    }
		    break;
		case no_correction:
		    for (uint i = 0; i < pairs.size(); i++) {
			if (pairs[i] != NULL) {
			    pairs[i]->stat[0] = pairs[i]->fst;
			    pairs[i]->stat[1] = pairs[i]->amova_fst;
			}
		    }
		    break;
		}

		//
		// If bootstrapping is enabled, record all Fst values.
		//
		if (bootstrap_fst)
		    bs->add_data(pairs);

		//
		// Calculate kernel-smoothed Fst values.
		//
		if (kernel_smoothed && loci_ordered) {
		    cerr << "  Generating kernel-smoothed Fst for " << it->first << ".\n";
		    ks->smooth(pairs);
		}
	    }

	    //
	    // If bootstrap resampling method is approximate, generate our single, empirical distribution.
	    //
	    map<int, vector<double> > approx_fst_dist;
	    // if (bootstrap_fst && bootstrap_type == bs_approx) 
	    // 	bootstrap_fst_approximate_dist(fst_samples, allele_depth_samples, weights, snp_dist, approx_fst_dist);

	    for (it = pmap->ordered_loci.begin(); it != pmap->ordered_loci.end(); it++) {
		string chr = it->first;
		vector<PopPair *> &pairs = genome_pairs[chr];

		//
		// Bootstrap resample this chromosome.
		//
		if (bootstrap_fst && bootstrap_type == bs_exact) {
		    cerr << "  Bootstrap resampling kernel-smoothed Fst for " << it->first << ".\n";
		    bs->execute(pairs);
		}

		for (uint i = 0; i < pairs.size(); i++) {

		    if (pairs[i] == NULL)
			continue;

		    //
		    // Calculate Fst P-value from approximate distribution.
		    //
		    // if (bootstrap_fst && bootstrap_type == bs_approx)
		    //     pairs[i]->bs[0] = bootstrap_approximate_pval(pairs[i]->snp_cnt, pairs[i]->stat[0], approx_fst_dist);

		    cnt++;
		    sum += pairs[i]->stat[1]; // Corrected AMOVA Fst

		    fh << batch_id          << "\t"
		       << pairs[i]->loc_id  << "\t"
		       << pop_key[pop_1]    << "\t"
		       << pop_key[pop_2]    << "\t"
		       << chr               << "\t"
		       << pairs[i]->bp      << "\t"
		       << pairs[i]->col     << "\t"
		       << pairs[i]->pi      << "\t"
		       << pairs[i]->fst     << "\t"
		       << std::setprecision(9) << pairs[i]->fet_p << "\t"
		       << pairs[i]->fet_or  << "\t"
		       << pairs[i]->ci_low  << "\t"
		       << pairs[i]->ci_high << "\t"
		       << pairs[i]->lod     << "\t"
		       << pairs[i]->stat[0]     << "\t"
		       << pairs[i]->smoothed[0] << "\t"
		       << pairs[i]->amova_fst   << "\t"
		       << pairs[i]->stat[1]     << "\t"
		       << pairs[i]->smoothed[1] << "\t"
		       << pairs[i]->bs[1] << "\t"
		       << pairs[i]->snp_cnt;

		    if (log_fst_comp) {
			fh << "\t"
			   << pairs[i]->comp[0]   << "\t"
			   << pairs[i]->comp[1]   << "\t"
			   << pairs[i]->comp[2]   << "\t"
			   << pairs[i]->comp[3]   << "\t"
			   << pairs[i]->comp[4]   << "\t"
			   << pairs[i]->comp[5]   << "\t"
			   << pairs[i]->comp[6]   << "\t"
			   << pairs[i]->comp[7]   << "\t"
			   << pairs[i]->comp[8]   << "\t"
			   << pairs[i]->comp[9]   << "\t"
			   << pairs[i]->comp[10]  << "\t"
			   << pairs[i]->comp[11]  << "\t"
			   << pairs[i]->fst       << "\t"
			   << pairs[i]->comp[12]  << "\t"
			   << pairs[i]->comp[13]  << "\t"
			   << pairs[i]->comp[14]  << "\t"
			   << pairs[i]->comp[15]  << "\t"
			   << pairs[i]->comp[16]  << "\t"
			   << pairs[i]->comp[17]  << "\t"
			   << pairs[i]->amova_fst << "\n";
		    } else {
			fh << "\n";
		    }

		    delete pairs[i];
		}
	    }
	    cerr << "Pop 1: " << pop_key[pop_1] << "; Pop 2: " << pop_key[pop_2] << "; mean Fst: " << (sum / cnt) << "\n";
	    means.push_back(sum / cnt);

	    cerr << "Pooled populations '" << pop_key[pop_1] << "' and '" << pop_key[pop_2] << "' contained: " << ord->incompatible_loci << " incompatible loci; " 
		 << ord->multiple_loci << " nucleotides covered by more than one RAD locus.\n";
	    fh.close();

	    if (bootstrap_fst)
		delete bs;
	}
    }

    //
    // Write out the mean Fst measure of each pair of populations.
    //
    stringstream pop_name;
    pop_name << "batch_" << batch_id << ".fst_summary.tsv";

    string file = in_path + pop_name.str();
    ofstream fh(file.c_str(), ofstream::out);

    if (fh.fail()) {
	cerr << "Error opening generic output file '" << file << "'\n";
	exit(1);
    }

    //
    // Write out X-axis header.
    //
    for (uint i = 0; i < pops.size(); i++)
	fh << "\t" << pop_key[pops[i]];
    fh << "\n";

    uint n = 0;
    for (uint i = 0; i < pops.size() - 1; i++) {
	fh << pop_key[pops[i]];

	for (uint k = 0; k <= i; k++)
	    fh << "\t";

	for (uint j = i + 1; j < pops.size(); j++) {
	    fh << "\t" << means[n];
	    n++;
	}
	fh << "\n";
    }

    fh.close();

    delete ord;
    if (kernel_smoothed && loci_ordered) {
	delete ks;
    }

    return 0;
}

int
correct_fst_bonferroni_win(vector<PopPair *> &pairs)
{
    int      limit = 3 * sigma;
    int      limit_l, limit_u;
    double   correction;
    uint     cnt, pos_l, pos_u;

    pos_l = 0;
    pos_u = 0;

    for (uint pos_c = 0; pos_c < pairs.size(); pos_c++) {
	if (pairs[pos_c] == NULL) continue;

	limit_l = pairs[pos_c]->bp - limit > 0 ? pairs[pos_c]->bp - limit : 0;
	limit_u = pairs[pos_c]->bp + limit;

	while (pos_l <  pairs.size()) {
	    if (pairs[pos_l] == NULL) {
		pos_l++;
	    } else {
		if (pairs[pos_l]->bp < limit_l) 
		    pos_l++;
		else
		    break;
	    }
	}    
	while (pos_u < pairs.size()) {
	    if (pairs[pos_u] == NULL) {
		pos_u++;
	    } else {
		if (pairs[pos_u]->bp < limit_u)
		    pos_u++;
		else
		    break;
	    }
	}

	cnt = 0;
	for (uint i = pos_l; i < pos_u; i++) {
	    if (pairs[i] == NULL) continue;
	    cnt++;
	}

	correction = p_value_cutoff / cnt;
	pairs[pos_c]->stat[0] = pairs[pos_c]->fet_p < correction ? pairs[pos_c]->fst : 0;
    }

    return 0;
}

int 
kernel_smoothed_popstats(map<int, CSLocus *> &catalog, PopMap<CSLocus> *pmap, PopSum<CSLocus> *psum, int pop_id, ofstream &log_fh) 
{
    // int snp_dist[max_snp_dist] = {0};
    // int sites_per_snp = 0;
    // int tot_windows = 0;
    map<string, vector<CSLocus *> >::iterator it;
    map<string, vector<SumStat *> > genome_sites;

    //
    // Instantiate the kernel smoothing object if requested.
    //
    KSmooth<SumStat>   *ks  = new KSmooth<SumStat>(2);
    OSumStat<SumStat>  *ord = new OSumStat<SumStat>(psum, log_fh);
    Bootstrap<SumStat> *bs;

    if (bootstrap_pifis)
	bs = new Bootstrap<SumStat>(2);

    for (it = pmap->ordered_loci.begin(); it != pmap->ordered_loci.end(); it++) {
	vector<SumStat *> &sites = genome_sites[it->first];

	ord->order(sites, it->second, pop_id);
	if (bootstrap_pifis) bs->add_data(sites);
    }

    cerr << "    Population '" << pop_key[pop_id] << "' contained " << ord->multiple_loci << " nucleotides covered by more than one RAD locus.\n";

    for (it = pmap->ordered_loci.begin(); it != pmap->ordered_loci.end(); it++) {
	if (bootstrap_pifis)
	    cerr << "    Smoothing and bootstrapping chromosome " << it->first << "\n";
	else
	    cerr << "    Smoothing chromosome " << it->first << "\n";

	vector<SumStat *> &sites = genome_sites[it->first];

	ks->smooth(sites);

	if (bootstrap_pifis && bootstrap_type == bs_exact)
	    bs->execute_mixed(sites);
    }

    delete ks;
    delete ord;
    if (bootstrap_pifis) delete bs;

//     //
//     // If bootstrap resampling method is approximate, generate our single, empirical distribution.
//     //
//     map<int, vector<double> > approx_fis_dist;
//     map<int, vector<double> > approx_pi_dist;
//     if (bootstrap && bootstrap_type == bs_approx) {
//     	sites_per_snp = sites_per_snp / tot_windows;

//     	// cerr << "Sites per snp: " << sites_per_snp << "\n";

//     	bootstrap_popstats_approximate_dist(fis_samples, pi_samples, allele_depth_samples, 
//     					    weights, snp_dist, sites_per_snp, 
//     					    approx_fis_dist, approx_pi_dist);

//     	for (it = pmap->ordered_loci.begin(); it != pmap->ordered_loci.end(); it++) {

//     	    for (uint pos = 0; pos < it->second.size(); pos++) {
//     		loc  = it->second[pos];
//     		len  = strlen(loc->con);
//     		lsum = psum->pop(loc->id, pop_id);

//     		for (int k = 0; k < len; k++)
//     		    if (lsum->nucs[k].num_indv > 0 && bootstrap && lsum->nucs[k].pi > 0) {
//     			//
//     			// Calculate Fis/Pi p-values from approximate distribution.
//     			//
//     			lsum->nucs[k].wFis_pval = bootstrap_approximate_pval(lsum->nucs[k].snp_cnt, lsum->nucs[k].wFis, approx_fis_dist);
//     			lsum->nucs[k].wPi_pval  = bootstrap_approximate_pval(lsum->nucs[k].snp_cnt, lsum->nucs[k].wPi, approx_pi_dist);
//     		    }
//     	    }
//     	}
//     }

    return 0;
}

int
bootstrap_popstats_approximate_dist(vector<double> &fis_samples,
				    vector<double> &pi_samples,
				    vector<int>  &allele_samples,
				    double *weights, int *snp_dist, int sites_per_snp,
				    map<int, vector<double> > &approx_fis_dist,
				    map<int, vector<double> > &approx_pi_dist)
{
    //
    // Allocate an array of bootstrap resampling objects.
    //
    int win_size = 6 * sigma + 1;
    int win_cntr = win_size / 2;

    //
    // Initialize the Fst distribution map.
    //
    for (int i = 0; i < max_snp_dist; i++) {
	if (snp_dist[i] == 0.0) continue;

	// cerr << "SNP Dist: " << i << " snps occurred " << snp_dist[i] << "\n";
	approx_fis_dist[i] = vector<double> ();
	approx_fis_dist[i].reserve(bootstrap_reps);

	approx_pi_dist[i] = vector<double> ();
	approx_pi_dist[i].reserve(bootstrap_reps);
    }

    vector<int> poss;
    poss.reserve(max_snp_dist);
    double weighted_fis, weighted_pi, sum_fis, sum_pi, final_weight_fis, final_weight_pi;
    // int    index_1, index_2;
    int    pos, index_3, dist, start, end;
    int    half = sites_per_snp / 2;

    for (int i = 0; i < max_snp_dist; i++) {
	if (snp_dist[i] == 0.0) continue;

	cerr << "  Generating NULL distribution for " << i << " SNPs...\n";

	// #pragma omp parallel private(poss, pos, index_1, index_2, index_3, dist, sum_fis, sum_pi, weighted_fis, weighted_pi, final_weight_fis, final_weight_pi)
        #pragma omp parallel private(poss, pos, index_3, dist, sum_fis, sum_pi, weighted_fis, weighted_pi, final_weight_fis, final_weight_pi)
	{ 
	    BSample *bs  = new BSample[win_size];

	    //
	    // Populate the BSample objects.
	    //
	    for (int n = 0; n < win_size;  n++)
		bs[n].bp = n + 1;

	    vector<double> fiss, pis;

	    //
	    // Bootstrap this bitch.
	    //
            #pragma omp for schedule(dynamic, 1)
	    for (int j = 0; j < bootstrap_reps; j++) {
		// cerr << "    Bootsrap rep " << j << "\n";

		//
		// First SNP is always placed at the center of the window.
		//
		pos     = win_cntr;
		// index_1 = (int) (fis_samples.size()    * (random() / (RAND_MAX + 1.0)));
		// index_2 = (int) (pi_samples.size()     * (random() / (RAND_MAX + 1.0)));
		index_3 = (int) (allele_samples.size() * (random() / (RAND_MAX + 1.0)));
		//
		// Fill in the area around the SNP with fixed sites.
		//
		start = pos - half > 0 ? pos - half : 0;
		end   = pos + half < win_size ? pos + half : win_size;
		for (int n = start; n < end; n++) {
		    // bs[n].f       = 0;
		    // bs[n].pi      = 0;
		    bs[n].alleles = bs[pos].alleles;
		    poss.push_back(n);
		}
		// bs[pos].f       = fis_samples[index_1];
		// bs[pos].pi      = pi_samples[index_2];
		bs[pos].alleles = allele_samples[index_3];
		// cerr << "      Placing SNP at position: " << pos << "; with data from " << index_1 << " filling area from " << start << " to " << end << "\n";

		//
		// Randomly select the positions and values for each SNP to populate the window
		// 
		for (int k = 0; k < i - 1; k++) {
		    pos     = (int) (win_size * (random() / (RAND_MAX + 1.0)));
		    // index_1 = (int) (fis_samples.size()    * (random() / (RAND_MAX + 1.0)));
		    // index_2 = (int) (pi_samples.size()     * (random() / (RAND_MAX + 1.0)));
		    index_3 = (int) (allele_samples.size() * (random() / (RAND_MAX + 1.0)));

		    poss.push_back(pos);
		    //
		    // Fill in the area around the SNP with fixed sites.
		    //
		    start = pos - half > 0 ? pos - half : 0;
		    end   = pos + half < win_size ? pos + half : win_size;
		    for (int n = start; n < end; n++) {
			// bs[n].f       = 0;
			// bs[n].pi      = 0;
			bs[n].alleles = bs[pos].alleles;
			poss.push_back(n);
		    }
		    // bs[pos].f       = fis_samples[index_1];
		    // bs[pos].pi      = pi_samples[index_2];
		    bs[pos].alleles = allele_samples[index_3];
		    // cerr << "      Placing SNP at position: " << pos << "; with data from " << index_1 << " filling area from " << start << " to " << end << "\n";
		}

		weighted_fis = 0.0;
		sum_fis      = 0.0;
		weighted_pi  = 0.0;
		sum_pi       = 0.0;

		for (int n = 0; n < win_size; n++) {
		    // if (bs[n].pi < 0.0)
		    // continue;
		    //
		    // Calculate weighted Fst at this position.
		    //
		    dist = bs[n].bp > bs[win_cntr].bp ? bs[n].bp - bs[win_cntr].bp : bs[win_cntr].bp - bs[n].bp;

		    final_weight_fis = (bs[n].alleles - 1) * weights[dist];
		    // weighted_fis    += bs[n].f * final_weight_fis;
		    sum_fis         += final_weight_fis;

		    final_weight_pi  = (bs[n].alleles - 1) * weights[dist];
		    // weighted_pi     += bs[n].pi * final_weight_pi;
		    sum_pi          += final_weight_pi;
		}

		fiss.push_back(weighted_fis / sum_fis);
		pis.push_back(weighted_pi  / sum_pi);
		// cerr << "      New weighted fis value: " << weighted_fis / sum_fis << "; size: " << fiss.size() << "\n";

		for (uint n = 0; n < poss.size(); n++) {
		    // bs[poss[n]].f  = 0.0;
		    // bs[poss[n]].pi = -1.0;
		}
		poss.clear();
	    }

// 	    #pragma omp critical
// 	    {
// 		vector<double> &f = approx_fis_dist[i];
// 		for (uint n = 0; n < fiss.size(); n++)
// 		    f.push_back(fiss[n]);
// 		vector<double> &p = approx_pi_dist[i];
// 		for (uint n = 0; n < pis.size(); n++)
// 		    p.push_back(pis[n]);
// 	    }

	    delete [] bs;
	}

	sort(approx_fis_dist[i].begin(), approx_fis_dist[i].end());
	sort(approx_pi_dist[i].begin(),  approx_pi_dist[i].end());
    }

    return 0;
}

int
bootstrap_fst_approximate_dist(vector<double> &fst_samples,
			       vector<int>  &allele_samples,
			       double *weights, int *snp_dist,
			       map<int, vector<double> > &approx_fst_dist)
{
    //
    // Allocate an array of bootstrap resampling objects.
    //
    int win_size = 6 * sigma + 1;
    int win_cntr = win_size / 2;

    //
    // Initialize the Fst distribution map.
    //
    for (int i = 0; i < max_snp_dist; i++) {
	if (snp_dist[i] == 0.0) continue;

	// cerr << "SNP Dist: " << i << " snps occurred " << snp_dist[i] << "\n";
	approx_fst_dist[i] = vector<double> ();
	approx_fst_dist[i].reserve(bootstrap_reps);
    }

    vector<int> poss;
    poss.reserve(max_snp_dist);
    double weighted_fst, sum, final_weight;
    //int    index_1;
    int    pos, index_2, dist;

    for (int i = 0; i < max_snp_dist; i++) {
	if (snp_dist[i] == 0.0) continue;

	cerr << "  Generating NULL distribution for " << i << " SNPs...\n";

	// #pragma omp parallel private(poss, pos, index_1, index_2, dist, sum, weighted_fst, final_weight)
        #pragma omp parallel private(poss, pos, index_2, dist, sum, weighted_fst, final_weight)
	{ 
	    BSample *bs  = new BSample[win_size];

	    //
	    // Populate the BSample objects.
	    //
	    for (int n = 0; n < win_size;  n++)
		bs[n].bp = n + 1;

	    vector<double> fsts;

	    //
	    // Bootstrap this bitch.
	    //
            #pragma omp for schedule(dynamic, 1)
	    for (int j = 0; j < bootstrap_reps; j++) {
		// cerr << "Bootsrap rep " << j << "\n";

		//
		// First SNP is always placed at the center of the window.
		//
		pos     = win_cntr;
		// index_1 = (int) (fst_samples.size() * (random() / (RAND_MAX + 1.0)));
		index_2 = (int) (allele_samples.size() * (random() / (RAND_MAX + 1.0)));
		// bs[pos].f       = fst_samples[index_1];
		bs[pos].alleles = allele_samples[index_2];

		//
		// Randomly select the positions and values for each SNP to populate the window
		// 
		for (int k = 0; k < i - 1; k++) {
		    pos     = (int) (win_size * (random() / (RAND_MAX + 1.0)));
		    // index_1 = (int) (fst_samples.size() * (random() / (RAND_MAX + 1.0)));
		    index_2 = (int) (allele_samples.size() * (random() / (RAND_MAX + 1.0)));
		    // bs[pos].f       = fst_samples[index_1];
		    // bs[pos].alleles = allele_samples[index_2];
		    // cerr << "  " << j << ": Placing SNP at position: " << pos << " with data from index " << index_1 << "\n";

		    poss.push_back(pos);
		}

		weighted_fst = 0.0;
		sum          = 0.0;

		for (int n = 0; n < win_size; n++) {
		    // if (bs[n].f == 0.0)
		    // continue;
		    //
		    // Calculate weighted Fst at this position.
		    //
		    dist = bs[n].bp > bs[win_cntr].bp ? bs[n].bp - bs[win_cntr].bp : bs[win_cntr].bp - bs[n].bp;

		    final_weight  = (bs[n].alleles - 1) * weights[dist];
		    // weighted_fst += bs[n].f * final_weight;
		    sum          += final_weight;
		}

		fsts.push_back(weighted_fst / sum);
		// cerr << "    New weighted Fst value: " << weighted_fst / sum << "; size: " << fsts.size() << "\n";

		// for (uint n = 0; n < poss.size(); n++)
		// bs[poss[n]].f = 0.0;
		poss.clear();
	    }

// 	    #pragma omp critical
// 	    {
// 		vector<double> &f = approx_fst_dist[i];
// 		for (uint n = 0; n < fsts.size(); n++)
// 		    f.push_back(fsts[n]);
// 	    }

	    delete [] bs;
	}

	sort(approx_fst_dist[i].begin(), approx_fst_dist[i].end());
    }

    return 0;
}

double
bootstrap_approximate_pval(int snp_cnt, double stat, map<int, vector<double> > &approx_dist)
{
    if (approx_dist.count(snp_cnt) == 0)
	return 1.0;

    vector<double>::iterator up;
    vector<double> &dist = approx_dist[snp_cnt];
    double pos;

    up  = upper_bound(dist.begin(), dist.end(), stat);

    if (up == dist.begin())
	pos = 1;
    else if (up == dist.end())
	pos = dist.size();
    else 
	pos = up - dist.begin() + 1;

    double res = 1.0 - (pos / (double) dist.size());

    // cerr << "Generated Approx Smoothed Fst Distribution:\n";
    // for (uint n = 0; n < dist.size(); n++)
    // 	cerr << "  n: " << n << "; Fst: " << dist[n] << "\n";

    // cerr << "Comparing Fst value: " << stat 
    // 	 << " at position " << (up - dist.begin()) << " out of " 
    // 	 << dist.size() << " positions (converted position: " << pos << "); pvalue: " << res << ".\n";

    return res;
}

int
write_generic(map<int, CSLocus *> &catalog, PopMap<CSLocus> *pmap, 
	      map<int, string> &samples, bool write_gtypes)
{
    stringstream pop_name;
    pop_name << "batch_" << batch_id;
    if (write_gtypes)
	pop_name << ".genotypes.tsv";
    else 
	pop_name << ".haplotypes.tsv";

    string file = in_path + pop_name.str();

    ofstream fh(file.c_str(), ofstream::out);

    if (fh.fail()) {
        cerr << "Error opening generic output file '" << file << "'\n";
	exit(1);
    }

    //
    // Count the number of markers that have enough samples to output.
    //
    map<int, CSLocus *>::iterator it;
    CSLocus *loc;
    int num_loci = catalog.size();

    cerr << "Writing " << num_loci << " loci to " << (write_gtypes ? "genotype" : "observed haplotype") << " file, '" << file << "'\n";

    //
    // Write the header
    //
    fh << "Catalog ID\t";
    if (expand_id)
	fh << "\t";
    if (write_gtypes)
	fh << "Marker\t";
    fh << "Cnt\t";

    for (int i = 0; i < pmap->sample_cnt(); i++) {
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

	Datum **d = pmap->locus(loc->id);
	string  obshap;

	for (int i = 0; i < pmap->sample_cnt(); i++) {
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

int 
write_sql(map<int, CSLocus *> &catalog, PopMap<CSLocus> *pmap) 
{
    stringstream pop_name;
    pop_name << "batch_" << batch_id << ".markers.tsv";
    string file = in_path + pop_name.str();

    cerr << "Writing SQL markers file to '" << file << "'\n";

    ofstream fh(file.c_str(), ofstream::out);
    if (fh.fail()) {
        cerr << "Error opening markers SQL file '" << file << "'\n";
	exit(1);
    }
    fh.precision(fieldw);
    fh.setf(std::ios::fixed);

    fh << "# SQL ID"            << "\t" 
       << "Batch ID"            << "\t" 
       << "Catalog Locus ID"    << "\t" 
       << "\t"
       << "Total Genotypes"     << "\t"
       << "Max"                 << "\t"
       << "Genotype Freqs"      << "\t"
       << "F"                   << "\t"
       << "Mean Log Likelihood" << "\t"
       << "Genotype Map"        << "\t"
       << "\n";

    map<int, CSLocus *>::iterator it;
    CSLocus *loc;
    stringstream gtype_map;

    for (it = catalog.begin(); it != catalog.end(); it++) {
	loc = it->second;

	string freq  = "";
	double max   = 0.0;
	int    total = 0;
	gtype_map.str("");

	if (loc->marker.length() > 0) {
	    tally_haplotype_freq(loc, pmap, total, max, freq);

	    //
	    // Record the haplotype to genotype map.
	    //
	    map<string, string>::iterator j;
	    for (j = loc->gmap.begin(); j != loc->gmap.end(); j++)
		gtype_map << j->first << ":" << j->second << ";";
	}

	fh << 0 << "\t" 
	   << batch_id << "\t" 
	   << loc->id  << "\t" 
	   << "\t"              // Marker
	   << total    << "\t"
	   << max      << "\t"
	   << freq     << "\t"
	   << loc->f   << "\t"
	   << loc->lnl << "\t"
           << gtype_map.str() << "\t"
	   << "\n";
    }

    fh.close();

    return 0;
}

int 
write_fasta(map<int, CSLocus *> &catalog, PopMap<CSLocus> *pmap, map<int, string> &samples, vector<int> &sample_ids) 
{
    //
    // Write a FASTA file containing each allele from each locus from 
    // each sample in the population.
    //
    stringstream pop_name;
    pop_name << "batch_" << batch_id << ".fa";
    string file = in_path + pop_name.str();

    cerr << "Writing population alleles to FASTA file '" << file << "'\n";

    ofstream fh(file.c_str(), ofstream::out);

    if (fh.fail()) {
        cerr << "Error opening FASTA file '" << file << "'\n";
	exit(1);
    }

    map<string, vector<CSLocus *> >::iterator it;
    CSLocus *loc;
    Datum  **d;
    char    *seq;

    for (it = pmap->ordered_loci.begin(); it != pmap->ordered_loci.end(); it++) {
	for (uint pos = 0; pos < it->second.size(); pos++) {
	    loc = it->second[pos];
	    d   = pmap->locus(loc->id);
	    seq = new char[loc->len + 1];
	    strcpy(seq, loc->con);

	    for (int j = 0; j < pmap->sample_cnt(); j++) {
		if (d[j] == NULL) 
		    continue;

		for (uint k = 0; k < d[j]->obshap.size(); k++) {

		    for (uint i = 0; i < loc->snps.size(); i++) {
			uint col = loc->snps[i]->col;
			seq[col] = col < loc->len ? d[j]->obshap[k][i] : loc->con[col];
		    }

		    fh << ">CLocus_" << loc->id 
		       << "_Sample_" << pmap->rev_sample_index(j)
		       << "_Locus_"  << d[j]->id 
		       << "_Allele_" << k
		       << " ["       << samples[pmap->rev_sample_index(j)];

		    if (strcmp(loc->loc.chr, "un") != 0)
			fh << "; " << loc->loc.chr << ", " << loc->sort_bp() + 1 << ", " << (loc->loc.strand == plus ? "+" : "-");
		    fh << "]\n"
		       << seq << "\n";
		}
	    }
	    delete [] seq;
	}
    }

    fh.close();

    return 0;
}

int 
write_strict_fasta(map<int, CSLocus *> &catalog, PopMap<CSLocus> *pmap, map<int, string> &samples, vector<int> &sample_ids) 
{
    //
    // Write a FASTA file containing each allele from each locus from 
    // each sample in the population.
    //
    stringstream pop_name;
    pop_name << "batch_" << batch_id << ".strict.fa";
    string file = in_path + pop_name.str();

    cerr << "Writing strict population alleles to FASTA file '" << file << "'\n";

    ofstream fh(file.c_str(), ofstream::out);

    if (fh.fail()) {
        cerr << "Error opening strict FASTA file '" << file << "'\n";
	exit(1);
    }

    map<string, vector<CSLocus *> >::iterator it;
    CSLocus *loc;
    Datum  **d;
    char    *seq;

    for (it = pmap->ordered_loci.begin(); it != pmap->ordered_loci.end(); it++) {
	for (uint pos = 0; pos < it->second.size(); pos++) {
	    loc = it->second[pos];
	    d   = pmap->locus(loc->id);
	    seq = new char[loc->len + 1];
	    strcpy(seq, loc->con);

	    for (int j = 0; j < pmap->sample_cnt(); j++) {
		if (d[j] == NULL) 
		    continue;
		if (d[j]->obshap.size() > 2) 
		    continue;

		if (d[j]->obshap.size() == 1) {

		    for (uint i = 0; i < loc->snps.size(); i++) {
			uint col = loc->snps[i]->col;
			seq[col] = col < loc->len ? d[j]->obshap[0][i] : loc->con[col];
		    }

		    fh << ">CLocus_" << loc->id 
		       << "_Sample_" << pmap->rev_sample_index(j)
		       << "_Locus_"  << d[j]->id 
		       << "_Allele_" << 0
		       << " ["       << samples[pmap->rev_sample_index(j)];
		    if (strcmp(loc->loc.chr, "un") != 0)
			fh << "; " << loc->loc.chr << ", " << loc->sort_bp() + 1 << ", " << (loc->loc.strand == plus ? "+" : "-");
		    fh << "]\n"
		       << seq << "\n";

		    fh << ">CLocus_" << loc->id 
		       << "_Sample_" << pmap->rev_sample_index(j)
		       << "_Locus_"  << d[j]->id 
		       << "_Allele_" << 1
		       << " ["       << samples[pmap->rev_sample_index(j)];
		    if (strcmp(loc->loc.chr, "un") != 0)
			fh << "; " << loc->loc.chr << ", " << loc->sort_bp() + 1 << ", " << (loc->loc.strand == plus ? "+" : "-");
		    fh << "]\n"
		       << seq << "\n";

		} else {
		    for (uint k = 0; k < d[j]->obshap.size(); k++) {
			for (uint i = 0; i < loc->snps.size(); i++) {
			    uint col = loc->snps[i]->col;
			    seq[col] = col < loc->len ? d[j]->obshap[k][i] : loc->con[col];
			}

			fh << ">CLocus_" << loc->id 
			   << "_Sample_" << pmap->rev_sample_index(j)
			   << "_Locus_"  << d[j]->id 
			   << "_Allele_" << k
			   << " ["       << samples[pmap->rev_sample_index(j)];
			if (strcmp(loc->loc.chr, "un") != 0)
			    fh << "; " << loc->loc.chr << ", " << loc->sort_bp() + 1 << ", " << (loc->loc.strand == plus ? "+" : "-");
			fh << "]\n"
			   << seq << "\n";
		    }
		}
	    }

	    delete [] seq;
	}
    }

    fh.close();

    return 0;
}

int 
write_vcf_ordered(map<int, CSLocus *> &catalog, 
		  PopMap<CSLocus> *pmap, PopSum<CSLocus> *psum, 
		  map<int, string> &samples, vector<int> &sample_ids,
		  map<int, pair<merget, int> > &merge_map, ofstream &log_fh) 
{
    //
    // Write a VCF file as defined here: http://www.1000genomes.org/node/101
    //
    stringstream pop_name;
    pop_name << "batch_" << batch_id << ".vcf";
    string file = in_path + pop_name.str();

    ofstream fh(file.c_str(), ofstream::out);

    if (fh.fail()) {
        cerr << "Error opening VCF file '" << file << "'\n";
	exit(1);
    }

    //
    // Load SNP data so that model likelihoods can be output to VCF file.
    //
    cerr << "In preparation for VCF export, loading SNP data for " << samples.size() << " samples.\n";

    populate_snp_calls(catalog, pmap, samples, sample_ids, merge_map);

    cerr << "Writing population data to VCF file '" << file << "'\n";

    log_fh << "\n#\n# Generating SNP-based VCF export.\n#\n";

    //
    // Obtain the current date.
    //
    time_t     rawtime;
    struct tm *timeinfo;
    char       date[32];
    time(&rawtime);
    timeinfo = localtime(&rawtime);
    strftime(date, 32, "%Y%m%d", timeinfo);

    //
    // Output the header.
    //
    fh << "##fileformat=VCFv4.0\n"
       << "##fileDate=" << date << "\n"
       << "##source=\"Stacks v" << VERSION << "\"\n"
       << "##INFO=<ID=NS,Number=1,Type=Integer,Description=\"Number of Samples With Data\">\n"
       << "##INFO=<ID=AF,Number=.,Type=Float,Description=\"Allele Frequency\">\n"
       << "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n"
       << "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Read Depth\">\n"
       << "##FORMAT=<ID=AD,Number=1,Type=Integer,Description=\"Allele Depth\">\n"
       << "##FORMAT=<ID=GL,Number=.,Type=Float,Description=\"Genotype Likelihood\">\n"
       << "#CHROM" << "\t" << "POS" << "\t" << "ID" << "\t" << "REF" << "\t" << "ALT" << "\t" 
       << "QUAL" << "\t" << "FILTER" << "\t" << "INFO" << "\t" << "FORMAT";

    for (int i = 0; i < pmap->sample_cnt(); i++)
	fh << "\t" << samples[pmap->rev_sample_index(i)];
    fh << "\n";    

    map<string, vector<CSLocus *> >::iterator it;
    CSLocus *loc;
    Datum  **d;
    int      gt_1, gt_2, dp_1, dp_2;
    char     p_allele, q_allele, p_str[32], q_str[32];
    uint16_t col;
    int      snp_index;

    //
    // We need to order the SNPs taking into account overlapping loci.
    //
    OLocTally<NucTally> *ord = new OLocTally<NucTally>(psum, log_fh);

    for (it = pmap->ordered_loci.begin(); it != pmap->ordered_loci.end(); it++) {
	vector<NucTally *> sites;
	ord->order(sites, it->second);

     	for (uint pos = 0; pos < sites.size(); pos++) {
	    if (catalog.count(sites[pos]->loc_id) == 0) {
		cerr << "Unable to find locus id " << sites[pos]->loc_id << "\n";
		continue;
	    }
    	    loc = catalog[sites[pos]->loc_id];
	    col = sites[pos]->col;

	    sprintf(p_str, "%0.3f", sites[pos]->p_freq);
	    sprintf(q_str, "%0.3f", 1 - sites[pos]->p_freq);

	    //
	    // If on the negative strand, complement the alleles.
	    //
	    p_allele = loc->loc.strand == minus ? reverse(sites[pos]->p_allele) : sites[pos]->p_allele;
	    q_allele = loc->loc.strand == minus ? reverse(sites[pos]->q_allele) : sites[pos]->q_allele;

	    fh << loc->loc.chr << "\t" 
	       << loc->sort_bp(col) + 1 << "\t" 
	       << loc->id    << "\t"
	       << p_allele   << "\t"            // REFerence allele
	       << q_allele   << "\t"            // ALTernate allele
	       << "."        << "\t"            // QUAL
	       << "PASS"     << "\t"            // FILTER
	       << "NS="      << sites[pos]->num_indv << ";"   // INFO
	       << "AF="      << p_str << "," << q_str << "\t" // INFO
	       << "GT:DP:AD:GL";                // FORMAT

	    snp_index = loc->snp_index(col);
	    if (snp_index < 0) {
		cerr << "Warning, unable to locate SNP call in column " << col << " for catalog locus " << loc->id << "\n";
		fh << "\n";
		continue;
	    }

	    d = pmap->locus(loc->id);

	    for (int j = 0; j < pmap->sample_cnt(); j++) {
		fh << "\t";

		if (d[j] == NULL) {
		    //
		    // Data does not exist.
		    //
		    fh << "./.:0:.,.:.,.,.";
		} else if (d[j]->model[col] == 'U') {
		    //
		    // Data exists, but the model call was uncertain.
		    //
		    fh << "./.:" << d[j]->tot_depth << ":.,.:.,.,.";
		} else {
		    //
		    // Tally up the nucleotide calls.
		    //
		    tally_observed_haplotypes(d[j]->obshap, snp_index, p_allele, q_allele);

		    if (p_allele == 0 && q_allele == 0) {
			// More than two potential alleles.
			fh << "./.:" << d[j]->tot_depth << ":.,.:.,.,.";
		    } else {
			find_datum_allele_depths(d[j], snp_index, sites[pos]->p_allele, sites[pos]->q_allele, p_allele+q_allele, dp_1, dp_2);

			if (p_allele == 0) {
			    gt_1 = q_allele == sites[pos]->p_allele ? 0 : 1;
			    fh << gt_1 << "/" << gt_1 << ":" << d[j]->tot_depth << ":" << dp_1 << "," << dp_2;
			} else if (q_allele == 0) {
			    gt_1 = p_allele == sites[pos]->p_allele ? 0 : 1;
			    fh << gt_1 << "/" << gt_1 << ":" << d[j]->tot_depth << ":" << dp_1 << "," << dp_2;
			} else {
			    gt_1 = p_allele == sites[pos]->p_allele ? 0 : 1;
			    gt_2 = q_allele == sites[pos]->p_allele ? 0 : 1;
			    fh << gt_1 << "/" << gt_2 << ":" << d[j]->tot_depth << ":" << dp_1 << "," << dp_2;
			}
			//
			// Output the likelihood for this model call.
			//
			if (col < d[j]->snps.size()) {
			    fh << ":.," << d[j]->snps[col]->lratio << ",.";
			} else {
			    cerr << "Warning, unable to locate SNP call in column " << col << " for catalog locus " << loc->id << ", tag ID " << d[j]->id << "\n";
			    fh << ":.,.,.";
			}
		    }
		}
	    }
	    fh << "\n";
	}
    }
    fh.close();

    return 0;
}

int 
write_vcf(map<int, CSLocus *> &catalog, 
	  PopMap<CSLocus> *pmap, PopSum<CSLocus> *psum, 
	  map<int, string> &samples, vector<int> &sample_ids,
	  map<int, pair<merget, int> > &merge_map) 
{
    //
    // Write a VCF file as defined here: http://www.1000genomes.org/node/101
    //
    stringstream pop_name;
    pop_name << "batch_" << batch_id << ".vcf";
    string file = in_path + pop_name.str();

    ofstream fh(file.c_str(), ofstream::out);

    if (fh.fail()) {
        cerr << "Error opening VCF file '" << file << "'\n";
	exit(1);
    }

    cerr << "In preparation for VCF export, loading SNP data for " << samples.size() << " samples.\n";
    //
    // Load SNP data so that model likelihoods can be output to VCF file.
    //
    populate_snp_calls(catalog, pmap, samples, sample_ids, merge_map);

    //
    // Obtain the current date.
    //
    time_t     rawtime;
    struct tm *timeinfo;
    char       date[32];
    time(&rawtime);
    timeinfo = localtime(&rawtime);
    strftime(date, 32, "%Y%m%d", timeinfo);

    cerr << "Writing population data to VCF file '" << file << "'\n";

    //
    // Output the header.
    //
    fh << "##fileformat=VCFv4.0\n"
       << "##fileDate=" << date << "\n"
       << "##source=\"Stacks v" << VERSION << "\"\n"
       << "##INFO=<ID=NS,Number=1,Type=Integer,Description=\"Number of Samples With Data\">\n"
       << "##INFO=<ID=AF,Number=.,Type=Float,Description=\"Allele Frequency\">\n"
       << "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n"
       << "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Read Depth\">\n"
       << "##FORMAT=<ID=AD,Number=1,Type=Integer,Description=\"Allele Depth\">\n"
       << "##FORMAT=<ID=GL,Number=.,Type=Float,Description=\"Genotype Likelihood\">\n"
       << "#CHROM" << "\t" << "POS" << "\t" << "ID" << "\t" << "REF" << "\t" << "ALT" << "\t" 
       << "QUAL" << "\t" << "FILTER" << "\t" << "INFO" << "\t" << "FORMAT";

    for (int i = 0; i < pmap->sample_cnt(); i++)
	fh << "\t" << samples[pmap->rev_sample_index(i)];
    fh << "\n";    

    map<string, vector<CSLocus *> >::iterator it;
    CSLocus *loc;
    Datum   **d;
    LocTally *t;
    int       gt_1, gt_2, dp_1, dp_2;
    double    num_indv;
    char      p_allele, q_allele, p_str[32], q_str[32];
    int       snp_index;

    for (it = pmap->ordered_loci.begin(); it != pmap->ordered_loci.end(); it++) {

	//
	// We need to order the SNPs so negative and positive strand SNPs are properly ordered.
	//
	vector<GenPos> ordered_loci;
	uint           col;
    	for (uint pos = 0; pos < it->second.size(); pos++) {
    	    loc = it->second[pos];
	    t   = psum->locus_tally(loc->id);

    	    for (uint i = 0; i < loc->snps.size(); i++) {
    		col = loc->snps[i]->col;
		if (t->nucs[col].allele_cnt == 2)
		    ordered_loci.push_back(GenPos(loc->id, i, loc->sort_bp(col)));
	    }
	}
	sort(ordered_loci.begin(), ordered_loci.end(), compare_genpos);

	for (uint pos = 0; pos < ordered_loci.size(); pos++) {
    	    loc = catalog[ordered_loci[pos].id];
	    col = loc->snps[ordered_loci[pos].snp_index]->col;
	    t   = psum->locus_tally(loc->id);

	    num_indv = (double) t->nucs[col].num_indv;

	    sprintf(p_str, "%0.3f", t->nucs[col].p_freq);
	    sprintf(q_str, "%0.3f", 1 - t->nucs[col].p_freq);

	    //
	    // If on the negative strand, complement the alleles.
	    //
	    p_allele = loc->loc.strand == minus ? reverse(t->nucs[col].p_allele) : t->nucs[col].p_allele;
	    q_allele = loc->loc.strand == minus ? reverse(t->nucs[col].q_allele) : t->nucs[col].q_allele;

	    fh << loc->loc.chr << "\t" 
	       << loc->sort_bp(col) + 1 << "\t" 
	       << loc->id << "\t"
	       << p_allele << "\t"              // REFerence allele
	       << q_allele << "\t"              // ALTernate allele
	       << "."        << "\t"            // QUAL
	       << "PASS"     << "\t"            // FILTER
	       << "NS="      << num_indv << ";" // INFO
	       << "AF="      << p_str << "," << q_str << "\t" // INFO
	       << "GT:DP:AD:GL";                // FORMAT

	    snp_index = loc->snp_index(col);
	    if (snp_index < 0) {
		cerr << "Warning, unable to locate SNP call in column " << col << " for catalog locus " << loc->id << "\n";
		fh << "\n";
		continue;
	    }

	    d = pmap->locus(loc->id);

	    for (int j = 0; j < pmap->sample_cnt(); j++) {
		fh << "\t";

		if (d[j] == NULL) {
		    //
		    // Data does not exist.
		    //
		    fh << "./.:0:.,.:.,.,.";
		} else if (d[j]->model[col] == 'U') {
		    //
		    // Data exists, but the model call was uncertain.
		    //
		    fh << "./.:" << d[j]->tot_depth << ":.,.:.,.,.";
		} else {
		    //
		    // Tally up the nucleotide calls.
		    //
		    tally_observed_haplotypes(d[j]->obshap, ordered_loci[pos].snp_index, p_allele, q_allele);

		    if (p_allele == 0 && q_allele == 0) {
			// More than two potential alleles.
			fh << "./.:" << d[j]->tot_depth << ":.,.:.,.,.";
		    } else {
			find_datum_allele_depths(d[j], snp_index, t->nucs[col].p_allele, t->nucs[col].q_allele, p_allele+q_allele, dp_1, dp_2);

			if (p_allele == 0) {
			    gt_1 = q_allele == t->nucs[col].p_allele ? 0 : 1;
			    fh << gt_1 << "/" << gt_1 << ":" << d[j]->tot_depth << ":" << dp_1 << "," << dp_2;
			} else if (q_allele == 0) {
			    gt_1 = p_allele == t->nucs[col].p_allele ? 0 : 1;
			    fh << gt_1 << "/" << gt_1 << ":" << d[j]->tot_depth << ":" << dp_1 << "," << dp_2;
			} else {
			    gt_1 = p_allele == t->nucs[col].p_allele ? 0 : 1;
			    gt_2 = q_allele == t->nucs[col].p_allele ? 0 : 1;
			    fh << gt_1 << "/" << gt_2 << ":" << d[j]->tot_depth << ":" << dp_1 << "," << dp_2;
			}
			//
			// Output the likelihood measure for this model call.
			//
			if (snp_index >= 0) {
			    fh << ":.," << d[j]->snps[snp_index]->lratio << ",.";
			} else {
			    cerr << "Warning, unable to locate SNP call in column " << col << " for catalog locus " << loc->id << ", tag ID " << d[j]->id << "\n";
			    fh << ":.,.,.";
			}
		    }
		}
	    }
	    fh << "\n";
	}
    }
    fh.close();

    return 0;
}

int
populate_snp_calls(map<int, CSLocus *> &catalog, PopMap<CSLocus> *pmap,
		   map<int, string> &samples, vector<int> &sample_ids,
		   map<int, pair<merget, int> > &merge_map)
{
    map<int, CSLocus *>::iterator cit;
    map<int, SNPRes *>::iterator sit;
    CSLocus *loc;
    Datum   *datum;
    SNPRes  *snpr;
    SNP     *snp;

    for (uint i = 0; i < sample_ids.size(); i++) {
	map<int, SNPRes *> snpres;
	load_snp_calls(in_path + samples[sample_ids[i]], snpres);

    	for (cit = catalog.begin(); cit != catalog.end(); cit++) {
    	    loc = cit->second;
    	    datum = pmap->datum(loc->id, sample_ids[i]);

    	    if (datum != NULL && snpres.count(datum->id)) {

		if (merge_sites && merge_map.count(loc->id)) { 
		    datum_adjust_snp_positions(merge_map, loc, datum, snpres);
		} else {
		    //
		    // Deep copy the SNP objects.
		    //
		    snpr = snpres[datum->id];
		    for (uint j = 0; j < snpr->snps.size(); j++) {
			snp         = new SNP;
			snp->col    = snpr->snps[j]->col;
			snp->lratio = snpr->snps[j]->lratio;
			snp->rank_1 = snpr->snps[j]->rank_1;
			snp->rank_2 = snpr->snps[j]->rank_2;
			snp->rank_3 = snpr->snps[j]->rank_3;
			snp->rank_4 = snpr->snps[j]->rank_4;

			datum->snps.push_back(snp);
		    }
		}
    	    }
    	}

	for (sit = snpres.begin(); sit != snpres.end(); sit++)
	    delete sit->second;
    }

    return 0;
}

int
find_datum_allele_depths(Datum *d, int snp_index, char p_allele, char q_allele, int allele_cnt, int &dp_1, int &dp_2)
{
    dp_1 = 0;
    dp_2 = 0;

    if (allele_cnt == 1) {

	//
	// There is a single observed haplotype for this locus, e.g. GA.
	//
	if (d->obshap.size() == 1) {
	    if (d->obshap[0][snp_index] == p_allele) {
		dp_1 = d->depth[0];
		dp_2 = 0;
	    } else {
		dp_1 = 0;
		dp_2 = d->depth[0];
	    }
	} else {
	    //
	    // This SNP position is homozygous, but the locus is heterozygous, so there is more
	    // than one observed haplotype, e.g. GA / TA.
	    //
	    if (d->obshap[0][snp_index] == p_allele) {
		dp_1 = d->tot_depth;
		dp_2 = 0;
	    } else  {
		dp_1 = 0;
		dp_2 = d->tot_depth;
	    }
	}

    } else {
	//
	// This SNP position is heterozygous.
	//
	for (uint i = 0; i < d->obshap.size(); i++) {
	    if (d->obshap[i][snp_index] == p_allele)
		dp_1 = d->depth[i];
	    else if (d->obshap[i][snp_index] == q_allele)
		dp_2 = d->depth[i];
	}
    }

    if (dp_1 == 0 && dp_2 == 0)
	cerr << "Warning: Unable to find allele depths for datum " << d->id << "\n";

    return 0;
}

int 
write_vcf_haplotypes(map<int, CSLocus *> &catalog, 
		     PopMap<CSLocus> *pmap, PopSum<CSLocus> *psum, 
		     map<int, string> &samples, vector<int> &sample_ids) 
{
    //
    // Write a VCF file as defined here: http://samtools.github.io/hts-specs/
    //
    stringstream pop_name;
    pop_name << "batch_" << batch_id << ".haplotypes.vcf";
    string file = in_path + pop_name.str();

    cerr << "Writing population data haplotypes to VCF file '" << file << "'\n";

    ofstream fh(file.c_str(), ofstream::out);

    if (fh.fail()) {
        cerr << "Error opening VCF file '" << file << "'\n";
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
    strftime(date, 32, "%Y%m%d", timeinfo);

    //
    // Output the header.
    //
    fh << "##fileformat=VCFv4.2\n"
       << "##fileDate=" << date << "\n"
       << "##source=\"Stacks v" << VERSION << "\"\n"
       << "##INFO=<ID=NS,Number=1,Type=Integer,Description=\"Number of Samples With Data\">\n"
       << "##INFO=<ID=AF,Number=.,Type=Float,Description=\"Allele Frequency\">\n"
       << "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n"
       << "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Read Depth\">\n"
       << "#CHROM" << "\t" << "POS" << "\t" << "ID" << "\t" << "REF" << "\t" << "ALT" << "\t" 
       << "QUAL" << "\t" << "FILTER" << "\t" << "INFO" << "\t" << "FORMAT";

    for (int i = 0; i < pmap->sample_cnt(); i++)
	fh << "\t" << samples[pmap->rev_sample_index(i)];
    fh << "\n";    

    map<string, vector<CSLocus *> >::iterator it;
    map<string, double>::iterator hit;
    map<string, double> hap_freq;
    map<string, int>    hap_index;
    vector<pair<string, int> > ordered_hap;
    CSLocus  *loc;
    Datum   **d;
    double    num_indv, num_hap;
    char      allele[id_len];

    for (it = pmap->ordered_loci.begin(); it != pmap->ordered_loci.end(); it++) {
    	for (uint pos = 0; pos < it->second.size(); pos++) {
    	    loc = it->second[pos];
	    d   = pmap->locus(loc->id);

	    hap_freq.clear();
	    hap_index.clear();
	    ordered_hap.clear();

	    num_hap = count_haplotypes_at_locus(0, pmap->sample_cnt() - 1, d, hap_freq);

	    if (num_hap == 0 || hap_freq.size() == 1) 
		continue;

	    num_indv = num_hap / 2.0;

	    //
	    // Order the haplotypes according to most frequent. Record the ordered position or each 
	    // haplotype and convert them from counts to frequencies.
	    //
	    for (hit = hap_freq.begin(); hit != hap_freq.end(); hit++) {
		ordered_hap.push_back(make_pair(hit->first, hit->second));
		hit->second = hit->second / num_hap;
	    }
	    sort(ordered_hap.begin(), ordered_hap.end(), compare_pair_haplotype);
	    for (uint i = 0; i < ordered_hap.size(); i++)
		hap_index[ordered_hap[i].first] = i;

	    string alt_str, freq_str;
	    for (uint i = 1; i < ordered_hap.size(); i++) {
		alt_str  += ordered_hap[i].first;
		sprintf(allele, "%0.3f", hap_freq[ordered_hap[i].first]);
		freq_str += allele;
		if (i < ordered_hap.size() - 1) {
		    alt_str  += ",";
		    freq_str += ",";
		}
	    }

	    fh << loc->loc.chr         << "\t" 
	       << loc->sort_bp() + 1   << "\t" 
	       << loc->id              << "\t"
	       << ordered_hap[0].first << "\t"       // REFerence haplotypes
	       << alt_str        << "\t"             // ALTernate haplotypes
	       << "."            << "\t"             // QUAL
	       << "PASS"         << "\t"             // FILTER
	       << "NS="          << num_indv << ";"  // INFO
	       << "AF="          << freq_str << "\t" // INFO
	       << "GT:DP";                           // FORMAT

	    for (int j = 0; j < pmap->sample_cnt(); j++) {
		fh << "\t";

		if (d[j] == NULL) {
		    //
		    // Data does not exist.
		    //
		    fh << "./.:0";

		} else if (d[j]->obshap.size() > 2) { 
		    fh << "./.:" << d[j]->tot_depth;

		} else if (d[j]->obshap.size() == 1) {
		    if(uncalled_haplotype(d[j]->obshap[0]))
		        fh << "./.:" << d[j]->tot_depth;
		    else
			fh << hap_index[d[j]->obshap[0]] << "/" << hap_index[d[j]->obshap[0]] << ":" << d[j]->tot_depth;
		} else {
		    if(!uncalled_haplotype(d[j]->obshap[0]) &&
		       !uncalled_haplotype(d[j]->obshap[1]))
		        fh << hap_index[d[j]->obshap[0]] << "/" << hap_index[d[j]->obshap[1]] << ":" << d[j]->tot_depth;
		    else if (!uncalled_haplotype(d[j]->obshap[0]))
		        fh << hap_index[d[j]->obshap[0]] << "/" << "." << ":" << d[j]->tot_depth;
		    else if (!uncalled_haplotype(d[j]->obshap[1]))
		        fh << "." << "/" << hap_index[d[j]->obshap[1]] << ":" << d[j]->tot_depth;
	        }
	    }
	    fh << "\n";
	}
    }

    fh.close();

    return 0;
}

int 
write_genepop(map<int, CSLocus *> &catalog, 
	      PopMap<CSLocus> *pmap, 
	      PopSum<CSLocus> *psum, 
	      map<int, pair<int, int> > &pop_indexes, 
	      map<int, string> &samples) 
{
    //
    // Write a GenePop file as defined here: http://kimura.univ-montp2.fr/~rousset/Genepop.htm
    //
    stringstream pop_name;
    pop_name << "batch_" << batch_id << ".genepop";
    string file = in_path + pop_name.str();

    cerr << "Writing population data to GenePop file '" << file << "'\n";

    ofstream fh(file.c_str(), ofstream::out);

    if (fh.fail()) {
        cerr << "Error opening GenePop file '" << file << "'\n";
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
    strftime(date, 32, "%B %d, %Y", timeinfo);

    //
    // Output the header line.
    //
    fh << "Stacks version " << VERSION << "; Genepop version 4.1.3; " << date << "\n";

    map<int, pair<int, int> >::iterator pit;
    map<int, CSLocus *>::iterator it;
    CSLocus  *loc;
    Datum   **d;
    LocSum  **s;
    LocTally *t;
    int      start_index, end_index, col, pop_id;
    char     p_allele, q_allele;

    //
    // Determine how many loci will be output, then output all the loci on the second line, comma-separated.
    //
    uint cnt = 0;
    for (it = catalog.begin(); it != catalog.end(); it++) {
	loc = it->second;
	for (uint j = 0; j < loc->snps.size(); j++) {
	    col = loc->snps[j]->col;
	    t   = psum->locus_tally(loc->id);

	    if (t->nucs[col].allele_cnt != 2) 
		continue;
	    cnt++;
	}
    }

    uint i = 0;
    for (it = catalog.begin(); it != catalog.end(); it++) {
	loc = it->second;
	for (uint j = 0; j < loc->snps.size(); j++) {
	    col = loc->snps[j]->col;
	    t   = psum->locus_tally(loc->id);

	    // 
	    // If this site is fixed in all populations or has too many alleles don't output it.
	    //
	    if (t->nucs[col].allele_cnt != 2) 
		continue;
	    i++;
	    fh << loc->id << "_" << col;
	    if (i <  cnt) fh << ",";
	}
    }
    fh << "\n";

    map<char, string> nuc_map;
    nuc_map['A'] = "01";
    nuc_map['C'] = "02";
    nuc_map['G'] = "03";
    nuc_map['T'] = "04";

    for (pit = pop_indexes.begin(); pit != pop_indexes.end(); pit++) {
	pop_id      = psum->pop_index(pit->first);
	start_index = pit->second.first;
	end_index   = pit->second.second;

	fh << "pop\n";

	for (int j = start_index; j <= end_index; j++) {

	    fh << samples[pmap->rev_sample_index(j)] << ",";

	    for (it = catalog.begin(); it != catalog.end(); it++) {
		loc = it->second;
		d   = pmap->locus(loc->id);
		s   = psum->locus(loc->id);
		t   = psum->locus_tally(loc->id);

		for (i = 0; i < loc->snps.size(); i++) {
		    uint col = loc->snps[i]->col;

		    if (t->nucs[col].allele_cnt != 2) 
			continue;

		    if (s[pop_id]->nucs[col].incompatible_site ||
			s[pop_id]->nucs[col].filtered_site) {
			//
			// This site contains more than two alleles in this population or was filtered
			// due to a minor allele frequency that is too low.
			//
			fh << "\t0000";
		    } else if (d[j] == NULL) {
			//
			// Data does not exist.
			//
			fh << "\t0000";
		    } else if (d[j]->model[col] == 'U') {
			//
			// Data exists, but the model call was uncertain.
			//
			fh << "\t0000";
		    } else {
			//
			// Tally up the nucleotide calls.
			//
			tally_observed_haplotypes(d[j]->obshap, i, p_allele, q_allele);

			if (p_allele == 0 && q_allele == 0) {
			    // More than two potential alleles.
			    fh << "\t0000";
			} else if (p_allele == 0) {
			    fh << "\t" << nuc_map[q_allele] << nuc_map[q_allele];

			} else if (q_allele == 0) {
			    fh << "\t" << nuc_map[p_allele] << nuc_map[p_allele];

			} else {
			    fh << "\t" << nuc_map[p_allele] << nuc_map[q_allele];
			}
		    }
		}
	    }
	    fh << "\n";
	}
    }

    fh.close();

    return 0;
}

int 
write_genepop_ordered(map<int, CSLocus *> &catalog, 
		      PopMap<CSLocus> *pmap, 
		      PopSum<CSLocus> *psum, 
		      map<int, pair<int, int> > &pop_indexes, 
		      map<int, string> &samples, ofstream &log_fh) 
{
    //
    // Write a GenePop file as defined here: http://kimura.univ-montp2.fr/~rousset/Genepop.htm
    //
    stringstream pop_name;
    pop_name << "batch_" << batch_id << ".genepop";
    string file = in_path + pop_name.str();

    cerr << "Writing population data to GenePop file '" << file << "'\n";

    ofstream fh(file.c_str(), ofstream::out);

    if (fh.fail()) {
        cerr << "Error opening GenePop file '" << file << "'\n";
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
    strftime(date, 32, "%B %d, %Y", timeinfo);

    //
    // Output the header line.
    //
    fh << "Stacks version " << VERSION << "; Genepop version 4.1.3; " << date << "\n";

    map<string, vector<NucTally *> > genome_sites;
    map<int, pair<int, int> >::iterator pit;
    map<string, vector<CSLocus *> >::iterator it;
    CSLocus  *loc;
    Datum   **d;
    LocSum  **s;
    int      start_index, end_index, pop_id;
    uint     col, snp_index;
    char     p_allele, q_allele;

    //
    // We need to order the SNPs to take into account overlapping loci.
    //
    OLocTally<NucTally> *ord = new OLocTally<NucTally>(psum, log_fh);

    //
    // Output all the loci on the second line, comma-separated.
    //
    int chrs = pmap->ordered_loci.size();
    int cnt  = 0;
    for (it = pmap->ordered_loci.begin(); it != pmap->ordered_loci.end(); it++) {
	vector<NucTally *> &sites = genome_sites[it->first];
	ord->order(sites, it->second);
	cnt++;

    	for (uint pos = 0; pos < sites.size(); pos++) {
	    fh << sites[pos]->loc_id << "_" << sites[pos]->col;
	    if (cnt < chrs || pos < sites.size() - 1) fh << ",";
	}
    }
    fh << "\n";

    map<char, string> nuc_map;
    nuc_map['A'] = "01";
    nuc_map['C'] = "02";
    nuc_map['G'] = "03";
    nuc_map['T'] = "04";

    for (pit = pop_indexes.begin(); pit != pop_indexes.end(); pit++) {
	pop_id      = psum->pop_index(pit->first);
	start_index = pit->second.first;
	end_index   = pit->second.second;

	fh << "pop\n";

	for (int j = start_index; j <= end_index; j++) {

	    fh << samples[pmap->rev_sample_index(j)] << ",";

	    for (it = pmap->ordered_loci.begin(); it != pmap->ordered_loci.end(); it++) {
		vector<NucTally *> &sites = genome_sites[it->first];

		for (uint pos = 0; pos < sites.size(); pos++) {
		    loc = catalog[sites[pos]->loc_id];
		    s   = psum->locus(loc->id);
		    d   = pmap->locus(loc->id);
		    col = sites[pos]->col;

		    if (s[pop_id]->nucs[col].incompatible_site ||
			s[pop_id]->nucs[col].filtered_site) {
			//
			// This site contains more than two alleles in this population or was filtered
			// due to a minor allele frequency that is too low.
			//
			fh << "\t0000";
		    } else if (d[j] == NULL) {
			//
			// Data does not exist.
			//
			fh << "\t0000";
		    } else if (d[j]->model[col] == 'U') {
			//
			// Data exists, but the model call was uncertain.
			//
			fh << "\t0000";
		    } else {
			snp_index = loc->snp_index(col);
			//
			// Tally up the nucleotide calls.
			//
			tally_observed_haplotypes(d[j]->obshap, snp_index, p_allele, q_allele);

			if (p_allele == 0 && q_allele == 0) {
			    // More than two potential alleles.
			    fh << "\t0000";
			} else if (p_allele == 0) {
			    fh << "\t" << nuc_map[q_allele] << nuc_map[q_allele];

			} else if (q_allele == 0) {
			    fh << "\t" << nuc_map[p_allele] << nuc_map[p_allele];

			} else {
			    fh << "\t" << nuc_map[p_allele] << nuc_map[q_allele];
			}
		    }
		}
	    }
	    fh << "\n";
	}
    }

    fh.close();

    return 0;
}

int 
write_structure(map<int, CSLocus *> &catalog, 
		PopMap<CSLocus> *pmap, 
		PopSum<CSLocus> *psum, 
		map<int, pair<int, int> > &pop_indexes, 
		map<int, string> &samples) 
{
    //
    // Write a Structure file as defined here: http://pritch.bsd.uchicago.edu/structure.html
    //
    // To avoid linked SNPs (which Structure can't handle), we will only output the first
    // SNP from each variable locus.
    //
    stringstream pop_name;
    pop_name << "batch_" << batch_id << ".structure.tsv";
    string file = in_path + pop_name.str();

    cerr << "Writing population data to Structure file '" << file << "'...";

    ofstream fh(file.c_str(), ofstream::out);

    if (fh.fail()) {
        cerr << "Error opening Structure file '" << file << "'\n";
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
    strftime(date, 32, "%B %d, %Y", timeinfo);

    //
    // Output the header.
    //
    fh << "# Stacks v" << VERSION << "; " << " Structure v2.3; " << date << "\n";

    map<string, vector<CSLocus *> >::iterator it;
    CSLocus  *loc;
    Datum   **d;
    LocSum  **s;
    LocTally *t;

    for (it = pmap->ordered_loci.begin(); it != pmap->ordered_loci.end(); it++) {
    	for (uint pos = 0; pos < it->second.size(); pos++) {
    	    loc = it->second[pos];
	    t   = psum->locus_tally(loc->id);

    	    for (uint i = 0; i < loc->snps.size(); i++) {
    		uint col = loc->snps[i]->col;
		if (t->nucs[col].allele_cnt == 2)
		    fh << "\t" << loc->id << "_" << col;
	    }
	}
    }
    fh << "\n";

    map<char, string> nuc_map;
    nuc_map['A'] = "1";
    nuc_map['C'] = "2";
    nuc_map['G'] = "3";
    nuc_map['T'] = "4";

    map<int, pair<int, int> >::iterator pit;
    int       start_index, end_index, pop_id, p;
    char      p_allele, q_allele;

    for (pit = pop_indexes.begin(); pit != pop_indexes.end(); pit++) {
	p           = psum->pop_index(pit->first);
	pop_id      = pit->first;
	start_index = pit->second.first;
	end_index   = pit->second.second;

	for (int j = start_index; j <= end_index; j++) {
	    //
	    // Output all the loci for this sample, printing only the p allele
	    //
	    fh << samples[pmap->rev_sample_index(j)] << "\t" << pop_id;

	    for (it = pmap->ordered_loci.begin(); it != pmap->ordered_loci.end(); it++) {
		for (uint pos = 0; pos < it->second.size(); pos++) {
		    loc = it->second[pos];

		    s = psum->locus(loc->id);
		    d = pmap->locus(loc->id);
		    t = psum->locus_tally(loc->id);

		    for (uint i = 0; i < loc->snps.size(); i++) {
			uint col = loc->snps[i]->col;

			// 
			// If this site is fixed in all populations or has too many alleles don't output it.
			//
			if (t->nucs[col].allele_cnt != 2) 
			    continue;

			if (s[p]->nucs[col].incompatible_site ||
			    s[p]->nucs[col].filtered_site) {
			    //
			    // This site contains more than two alleles in this population or was filtered
			    // due to a minor allele frequency that is too low.
			    //
			    fh << "\t" << "0";
			} else if (d[j] == NULL) {
			    //
			    // Data does not exist.
			    //
			    fh << "\t" << "0";
			} else if (d[j]->model[col] == 'U') {
			    //
			    // Data exists, but the model call was uncertain.
			    //
			    fh << "\t" << "0";
			} else {
			    //
			    // Tally up the nucleotide calls.
			    //
			    tally_observed_haplotypes(d[j]->obshap, i, p_allele, q_allele);

			    if (p_allele == 0 && q_allele == 0)
				fh << "\t" << "0";
			    else if (p_allele == 0)
				fh << "\t" << nuc_map[q_allele];
			    else
				fh << "\t" << nuc_map[p_allele];
			}
		    }
		}
    	    }
	    fh << "\n";

	    //
	    // Output all the loci for this sample again, now for the q allele
	    //
	    fh << samples[pmap->rev_sample_index(j)] << "\t" << pop_id;

	    for (it = pmap->ordered_loci.begin(); it != pmap->ordered_loci.end(); it++) {
		for (uint pos = 0; pos < it->second.size(); pos++) {
		    loc = it->second[pos];

		    s = psum->locus(loc->id);
		    d = pmap->locus(loc->id);
		    t = psum->locus_tally(loc->id);

		    for (uint i = 0; i < loc->snps.size(); i++) {
			uint col = loc->snps[i]->col;

			if (t->nucs[col].allele_cnt != 2) 
			    continue;

			if (s[p]->nucs[col].incompatible_site ||
			    s[p]->nucs[col].filtered_site) {
			    fh << "\t" << "0";
			} else if (d[j] == NULL) {
			    fh << "\t" << "0";
			} else if (d[j]->model[col] == 'U') {
			    fh << "\t" << "0";
			} else {
			    tally_observed_haplotypes(d[j]->obshap, i, p_allele, q_allele);

			    if (p_allele == 0 && q_allele == 0)
				fh << "\t" << "0";
			    else if (q_allele == 0)
				fh << "\t" << nuc_map[p_allele];
			    else
				fh << "\t" << nuc_map[q_allele];
			}
		    }
		}
    	    }
	    fh << "\n";
    	}
    }

    fh.close();

    cerr << "done.\n";

    return 0;
}

int 
write_structure_ordered(map<int, CSLocus *> &catalog, 
			PopMap<CSLocus> *pmap, 
			PopSum<CSLocus> *psum, 
			map<int, pair<int, int> > &pop_indexes, 
			map<int, string> &samples, ofstream &log_fh) 
{
    //
    // Write a Structure file as defined here: http://pritch.bsd.uchicago.edu/structure.html
    //
    // To avoid linked SNPs (which Structure can't handle), we will only output the first
    // SNP from each variable locus.
    //
    stringstream pop_name;
    pop_name << "batch_" << batch_id << ".structure.tsv";
    string file = in_path + pop_name.str();

    cerr << "Writing population data to Structure file '" << file << "'...";

    ofstream fh(file.c_str(), ofstream::out);

    if (fh.fail()) {
        cerr << "Error opening Structure file '" << file << "'\n";
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
    strftime(date, 32, "%B %d, %Y", timeinfo);

    //
    // Output the header.
    //
    fh << "# Stacks v" << VERSION << "; " << " Structure v2.3; " << date << "\n";

    map<string, vector<NucTally *> > genome_sites;
    map<string, vector<CSLocus *> >::iterator it;
    CSLocus  *loc;
    Datum   **d;
    LocSum  **s;

    //
    // We need to order the SNPs to take into account overlapping loci.
    //
    OLocTally<NucTally> *ord = new OLocTally<NucTally>(psum, log_fh);

    for (it = pmap->ordered_loci.begin(); it != pmap->ordered_loci.end(); it++) {
	vector<NucTally *> &sites = genome_sites[it->first];
	ord->order(sites, it->second);

    	for (uint pos = 0; pos < sites.size(); pos++)
	    fh << "\t" << sites[pos]->loc_id << "_" << sites[pos]->col;
    }
    fh << "\n";

    map<char, string> nuc_map;
    nuc_map['A'] = "1";
    nuc_map['C'] = "2";
    nuc_map['G'] = "3";
    nuc_map['T'] = "4";

    map<int, pair<int, int> >::iterator pit;
    int       start_index, end_index, pop_id, p;
    char      p_allele, q_allele;
    uint      col, snp_index;

    for (pit = pop_indexes.begin(); pit != pop_indexes.end(); pit++) {
	p           = psum->pop_index(pit->first);
	pop_id      = pit->first;
	start_index = pit->second.first;
	end_index   = pit->second.second;

	for (int j = start_index; j <= end_index; j++) {
	    //
	    // Output all the loci for this sample, printing only the p allele
	    //
	    fh << samples[pmap->rev_sample_index(j)] << "\t" << pop_id;

	    for (it = pmap->ordered_loci.begin(); it != pmap->ordered_loci.end(); it++) {
		vector<NucTally *> &sites = genome_sites[it->first];

		for (uint pos = 0; pos < sites.size(); pos++) {
		    loc = catalog[sites[pos]->loc_id];
		    s   = psum->locus(loc->id);
		    d   = pmap->locus(loc->id);
		    col = sites[pos]->col;

		    if (s[p]->nucs[col].incompatible_site ||
			s[p]->nucs[col].filtered_site) {
			//
			// This site contains more than two alleles in this population or was filtered
			// due to a minor allele frequency that is too low.
			//
			fh << "\t" << "0";
		    } else if (d[j] == NULL) {
			//
			// Data does not exist.
			//
			fh << "\t" << "0";
		    } else if (d[j]->model[col] == 'U') {
			//
			// Data exists, but the model call was uncertain.
			//
			fh << "\t" << "0";
		    } else {
			snp_index = loc->snp_index(col);
			//
			// Tally up the nucleotide calls.
			//
			tally_observed_haplotypes(d[j]->obshap, snp_index, p_allele, q_allele);

			if (p_allele == 0 && q_allele == 0)
			    fh << "\t" << "0";
			else if (p_allele == 0)
			    fh << "\t" << nuc_map[q_allele];
			else
			    fh << "\t" << nuc_map[p_allele];
		    }
		}
	    }
	    fh << "\n";

	    //
	    // Output all the loci for this sample again, now for the q allele
	    //
	    fh << samples[pmap->rev_sample_index(j)] << "\t" << pop_id;

	    for (it = pmap->ordered_loci.begin(); it != pmap->ordered_loci.end(); it++) {
		vector<NucTally *> &sites = genome_sites[it->first];

		for (uint pos = 0; pos < sites.size(); pos++) {
		    loc = catalog[sites[pos]->loc_id];
		    s   = psum->locus(loc->id);
		    d   = pmap->locus(loc->id);
		    col = sites[pos]->col;

		    if (s[p]->nucs[col].incompatible_site ||
			s[p]->nucs[col].filtered_site) {
			fh << "\t" << "0";
		    } else if (d[j] == NULL) {
			fh << "\t" << "0";
		    } else if (d[j]->model[col] == 'U') {
			fh << "\t" << "0";
		    } else {
			snp_index = loc->snp_index(col);
			tally_observed_haplotypes(d[j]->obshap, snp_index, p_allele, q_allele);

			if (p_allele == 0 && q_allele == 0)
			    fh << "\t" << "0";
			else if (q_allele == 0)
			    fh << "\t" << nuc_map[p_allele];
			else
			    fh << "\t" << nuc_map[q_allele];
		    }
		}
	    }
	    fh << "\n";
	}
    }

    fh.close();

    cerr << "done.\n";

    return 0;
}

int 
write_hzar(map<int, CSLocus *> &catalog, 
	   PopMap<CSLocus> *pmap, 
	   PopSum<CSLocus> *psum, 
	   map<int, pair<int, int> > &pop_indexes, 
	   map<int, string> &samples) 
{
    //
    // Write a Hybrid Zone Analysis using R (HZAR) file as defined here: 
    //    http://cran.r-project.org/web/packages/hzar/hzar.pdf
    //
    stringstream pop_name;
    pop_name << "batch_" << batch_id << ".hzar.csv";
    string file = in_path + pop_name.str();

    cerr << "Writing population data to HZAR file '" << file << "'...";

    ofstream fh(file.c_str(), ofstream::out);

    if (fh.fail()) {
        cerr << "Error opening HZAR file '" << file << "'\n";
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
    strftime(date, 32, "%B %d, %Y", timeinfo);

    //
    // Output the header.
    //
    fh << "# Stacks v" << VERSION << "; " << " HZAR v0.2-5; " << date << "\n"
       << "Population" << "," << "Distance";

    map<string, vector<CSLocus *> >::iterator it;
    CSLocus  *loc;
    LocSum  **s;
    LocTally *t;

    for (it = pmap->ordered_loci.begin(); it != pmap->ordered_loci.end(); it++) {
    	for (uint pos = 0; pos < it->second.size(); pos++) {
    	    loc = it->second[pos];
    	    s   = psum->locus(loc->id);
	    t   = psum->locus_tally(loc->id);

    	    for (uint i = 0; i < loc->snps.size(); i++) {
    		uint col = loc->snps[i]->col;
		if (t->nucs[col].allele_cnt == 2) {
		    fh << "," << loc->id << "_" << col << ".A"
		       << "," << loc->id << "_" << col << ".B"
		       << "," << loc->id << "_" << col << ".N";
		}
	    }
	}
    }
    fh << "\n";

    map<int, pair<int, int> >::iterator pit;
    int pop_id, p;

    for (pit = pop_indexes.begin(); pit != pop_indexes.end(); pit++) {
	p      = psum->pop_index(pit->first);
	pop_id = pit->first;

	fh << pop_key[pop_id] << ",";

	for (it = pmap->ordered_loci.begin(); it != pmap->ordered_loci.end(); it++) {
	    for (uint pos = 0; pos < it->second.size(); pos++) {
		loc = it->second[pos];

		s = psum->locus(loc->id);
		t = psum->locus_tally(loc->id);

		for (uint i = 0; i < loc->snps.size(); i++) {
		    uint col = loc->snps[i]->col;

		    // 
		    // If this site is fixed in all populations or has too many alleles don't output it.
		    //
		    if (t->nucs[col].allele_cnt != 2) 
			continue;

		    if (s[p]->nucs[col].num_indv == 0 ||
			s[p]->nucs[col].incompatible_site ||
			s[p]->nucs[col].filtered_site) {
			fh << ",0,0,0";
			continue;
		    }
		    
		    if (t->nucs[col].p_allele == s[p]->nucs[col].p_nuc)
			fh << "," << s[p]->nucs[col].p << "," << 1 - s[p]->nucs[col].p << ",";
		    else
			fh << "," << 1 - s[p]->nucs[col].p << "," << s[p]->nucs[col].p << ",";

		    fh << s[p]->nucs[col].num_indv * 2;
		}
	    }
    	}
	fh << "\n";
    }

    fh.close();

    cerr << "done.\n";

    return 0;
}

int 
write_treemix(map<int, CSLocus *> &catalog, 
	      PopMap<CSLocus> *pmap, 
	      PopSum<CSLocus> *psum, 
	      map<int, pair<int, int> > &pop_indexes, 
	      map<int, string> &samples) 
{
    //
    // Write a TreeMix file (Pickrell and Pritchard, 2012 PLoS Genetics)
    //    https://bitbucket.org/nygcresearch/treemix/wiki/Home
    //
    stringstream pop_name;
    pop_name << "batch_" << batch_id << ".treemix";
    string file = in_path + pop_name.str();

    cerr << "Writing population data to TreeMix file '" << file << "'; ";

    ofstream fh(file.c_str(), ofstream::out);

    if (fh.fail()) {
        cerr << "Error opening TreeMix file '" << file << "'\n";
    	exit(1);
    }

    pop_name << ".log";
    file = in_path + pop_name.str();

    cerr << "logging nucleotide positions to '" << file << "'...";

    ofstream log_fh(file.c_str(), ofstream::out);

    if (log_fh.fail()) {
        cerr << "Error opening Phylip Log file '" << file << "'\n";
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
    strftime(date, 32, "%B %d, %Y", timeinfo);

    log_fh << "# Stacks v" << VERSION << "; " << " TreeMix v1.1; " << date << "\n"
	   << "# Line\tLocus ID\tColumn\tChr\tBasepair\n";

    //
    // Output the header.
    //
    fh << "# Stacks v" << VERSION << "; " << " TreeMix v1.1; " << date << "\n";

    map<string, vector<CSLocus *> >::iterator it;
    map<int, pair<int, int> >::iterator pit;
    CSLocus  *loc;
    LocSum  **s;
    LocTally *t;
    int       p;

    //
    // Output a space-separated list of the populations on the first line.
    //
    stringstream sstr;
    for (pit = pop_indexes.begin(); pit != pop_indexes.end(); pit++)
	sstr << pop_key[pit->first] << " ";
    
    fh << sstr.str().substr(0, sstr.str().length() - 1) << "\n";

    double p_freq, p_cnt, q_cnt, allele_cnt;
    long int line = 1;
    
    for (it = pmap->ordered_loci.begin(); it != pmap->ordered_loci.end(); it++) {
	for (uint pos = 0; pos < it->second.size(); pos++) {
	    loc = it->second[pos];

	    s = psum->locus(loc->id);
	    t = psum->locus_tally(loc->id);
		
	    for (uint i = 0; i < loc->snps.size(); i++) {
		uint col = loc->snps[i]->col;

		sstr.str("");

		// 
		// If this site is fixed in all populations or has too many alleles don't output it.
		//
		if (t->nucs[col].allele_cnt != 2) 
		    continue;

		for (pit = pop_indexes.begin(); pit != pop_indexes.end(); pit++) {
		    p = psum->pop_index(pit->first);

		    if (s[p]->nucs[col].num_indv == 0 ||
			s[p]->nucs[col].incompatible_site ||
			s[p]->nucs[col].filtered_site) {
			sstr << "0,0 ";
			continue;
		    }
		    
		    p_freq = (t->nucs[col].p_allele == s[p]->nucs[col].p_nuc) ? 
			s[p]->nucs[col].p :
			1 - s[p]->nucs[col].p;
		    
		    allele_cnt = s[p]->nucs[col].num_indv * 2;
		    p_cnt      = round(allele_cnt * p_freq);
		    q_cnt      = allele_cnt - p_cnt;
		    sstr << (int) p_cnt << "," << (int) q_cnt << " ";
		}

		if (sstr.str().length() == 0)
		    continue;

		fh << sstr.str().substr(0, sstr.str().length() - 1) << "\n";
		log_fh << line << "\t" << loc->id << "\t" << col << "\t" << loc->loc.chr << "\t" << loc->sort_bp(col) + 1 << "\n";
		line++;
	    }
    	}
    }

    fh.close();
    log_fh.close();
    
    cerr << "done.\n";

    return 0;
}

int 
write_fastphase(map<int, CSLocus *> &catalog, 
		PopMap<CSLocus> *pmap, 
		PopSum<CSLocus> *psum, 
		map<int, pair<int, int> > &pop_indexes, 
		map<int, string> &samples) 
{
    //
    // Write a fastPHASE file as defined here: http://stephenslab.uchicago.edu/software.html
    //
    // Data will be written as independent, bi-allelic SNPs. We will write one file per chromosome.
    //
    cerr << "Writing population data to fastPHASE files...";

    map<string, vector<CSLocus *> >::iterator it;
    CSLocus  *loc;
    Datum   **d;
    LocSum  **s;
    LocTally *t;

    for (it = pmap->ordered_loci.begin(); it != pmap->ordered_loci.end(); it++) {

	stringstream pop_name;
	pop_name << "batch_" << batch_id << "." << it->first << ".fastphase.inp";
	string file = in_path + pop_name.str();

	ofstream fh(file.c_str(), ofstream::out);

	if (fh.fail()) {
	    cerr << "Error opening fastPHASE file '" << file << "'\n";
	    exit(1);
	}

	//
	// Tally up the number of sites
	//
	int  total_sites = 0;
	uint col;
    	for (uint pos = 0; pos < it->second.size(); pos++) {
    	    loc = it->second[pos];
	    t   = psum->locus_tally(loc->id);

    	    for (uint i = 0; i < loc->snps.size(); i++) {
    		col = loc->snps[i]->col;
		if (t->nucs[col].allele_cnt == 2)
		    total_sites++;
	    }
	}

	//
	// Output the total number of SNP sites and the number of individuals.
	//
	fh << samples.size() << "\n"
	   << total_sites    << "\n";

	//
	// We need to determine an ordering that can take into account overlapping RAD sites.
	//
	vector<GenPos> ordered_loci;
    	for (uint pos = 0; pos < it->second.size(); pos++) {
    	    loc = it->second[pos];
	    t   = psum->locus_tally(loc->id);

    	    for (uint i = 0; i < loc->snps.size(); i++) {
    		col = loc->snps[i]->col;
		if (t->nucs[col].allele_cnt == 2)
		    ordered_loci.push_back(GenPos(loc->id, i, loc->sort_bp(col)));
	    }
	}
	sort(ordered_loci.begin(), ordered_loci.end(), compare_genpos);

	//
	// Output the position of each site according to its basepair.
	//
	fh << "P";
    	for (uint pos = 0; pos < ordered_loci.size(); pos++) {
    	    loc = catalog[ordered_loci[pos].id];
	    col = loc->snps[ordered_loci[pos].snp_index]->col;
	    fh << " " << ordered_loci[pos].bp;
	}
	fh << "\n";

	//
	// Output a line of 'S' characters, one per site, indicating that these are SNP markers.
	//
	string snp_markers, gtypes_str;
	snp_markers.assign(total_sites, 'S');
	fh << snp_markers << '\n';	    

	//
	// Now output each sample name followed by a new line, then all of the genotypes for that sample
	// on two lines.
	//

	map<int, pair<int, int> >::iterator pit;
	int          start_index, end_index, pop_id;
	char         p_allele, q_allele;
	stringstream gtypes;

	for (pit = pop_indexes.begin(); pit != pop_indexes.end(); pit++) {
	    pop_id      = psum->pop_index(pit->first);
	    start_index = pit->second.first;
	    end_index   = pit->second.second;

	    for (int j = start_index; j <= end_index; j++) {
		//
		// Output all the loci for this sample, printing only the p allele
		//
		fh << samples[pmap->rev_sample_index(j)] << "\n";
		
		gtypes.str("");
		for (uint pos = 0; pos < ordered_loci.size(); pos++) {
		    loc = catalog[ordered_loci[pos].id];
		    col = loc->snps[ordered_loci[pos].snp_index]->col;

		    s = psum->locus(loc->id);
		    d = pmap->locus(loc->id);
		    t = psum->locus_tally(loc->id);

		    // 
		    // If this site is fixed in all populations or has too many alleles don't output it.
		    //
		    if (t->nucs[col].allele_cnt != 2) 
			continue;

		    if (s[pop_id]->nucs[col].incompatible_site ||
			s[pop_id]->nucs[col].filtered_site) {
			//
			// This site contains more than two alleles in this population or was filtered
			// due to a minor allele frequency that is too low.
			//
			gtypes << "? ";

		    } else if (d[j] == NULL) {
			//
			// Data does not exist.
			//
			gtypes << "? ";
		    } else if (d[j]->model[col] == 'U') {
			//
			// Data exists, but the model call was uncertain.
			//
			gtypes << "? ";
		    } else {
			//
			// Tally up the nucleotide calls.
			//
			tally_observed_haplotypes(d[j]->obshap, ordered_loci[pos].snp_index, p_allele, q_allele);

			if (p_allele == 0 && q_allele == 0)
			    gtypes << "? ";
			else if (p_allele == 0)
			    gtypes << q_allele << " ";
			else
			    gtypes << p_allele << " ";
		    }
		}
		gtypes_str = gtypes.str();
		fh << gtypes_str.substr(0, gtypes_str.length() - 1) << "\n";

		//
		// Output all the loci for this sample again, now for the q allele
		//
		gtypes.str("");
		for (uint pos = 0; pos < ordered_loci.size(); pos++) {
		    loc = catalog[ordered_loci[pos].id];
		    col = loc->snps[ordered_loci[pos].snp_index]->col;


		    s = psum->locus(loc->id);
		    d = pmap->locus(loc->id);
		    t = psum->locus_tally(loc->id);

		    if (t->nucs[col].allele_cnt != 2) 
			continue;

		    if (s[pop_id]->nucs[col].incompatible_site ||
			s[pop_id]->nucs[col].filtered_site) {
			gtypes << "? ";

		    } else if (d[j] == NULL) {
			gtypes << "? ";

		    } else if (d[j]->model[col] == 'U') {
			gtypes << "? ";

		    } else {
			tally_observed_haplotypes(d[j]->obshap, ordered_loci[pos].snp_index, p_allele, q_allele);

			if (p_allele == 0 && q_allele == 0)
			    gtypes << "? ";
			else if (q_allele == 0)
			    gtypes << p_allele << " ";
			else
			    gtypes << q_allele << " ";
		    }
		}
		gtypes_str = gtypes.str();
		fh << gtypes_str.substr(0, gtypes_str.length() - 1) << "\n";
	    }
	}

	fh.close();
    }

    cerr << "done.\n";

    return 0;
}

int 
write_phase(map<int, CSLocus *> &catalog, 
	    PopMap<CSLocus> *pmap, 
	    PopSum<CSLocus> *psum, 
	    map<int, pair<int, int> > &pop_indexes, 
	    map<int, string> &samples) 
{
    //
    // Write a PHASE file as defined here: http://stephenslab.uchicago.edu/software.html
    //
    // Data will be written as mixture of multiple allele, linked RAD sites 
    // (SNPs within a single RAD locus are already phased), and bi-allelic SNPs. We 
    // will write one file per chromosome.
    //
    cerr << "Writing population data to PHASE files...";

    map<string, vector<CSLocus *> >::iterator it;
    CSLocus  *loc;
    Datum   **d;
    LocSum  **s;
    LocTally *t;

    for (it = pmap->ordered_loci.begin(); it != pmap->ordered_loci.end(); it++) {

	stringstream pop_name;
	pop_name << "batch_" << batch_id << "." << it->first << ".phase.inp";
	string file = in_path + pop_name.str();

	ofstream fh(file.c_str(), ofstream::out);

	if (fh.fail()) {
	    cerr << "Error opening PHASE file '" << file << "'\n";
	    exit(1);
	}

	//
	// We need to determine an ordering for all legitimate loci/SNPs.
	//
	uint           col;
	vector<GenPos> ordered_loci;
    	for (uint pos = 0; pos < it->second.size(); pos++) {
    	    loc = it->second[pos];
	    t   = psum->locus_tally(loc->id);

	    if (loc->snps.size() == 0) continue;

	    //
	    // Will we output this locus as a haplotype or as a SNP?
	    //
	    if (loc->snps.size() > 1) {
		//
		// Check that there aren't too many haplotypes (PHASE has a max of 50).
		//
		if (loc->alleles.size() > 40) continue;

		//
		// Iterate over the population to determine that this subset of the population
		// has data at this locus.
		//
		d = pmap->locus(loc->id);
		for (int j = 0; j < pmap->sample_cnt(); j++) {
		    if (d[j] != NULL &&
			d[j]->obshap.size() > 0 && 
			d[j]->obshap.size() <= 2) {
			//
			// Data exists, and there are the correct number of haplotypes.
			//
			ordered_loci.push_back(GenPos(loc->id, 0, loc->sort_bp(), haplotype));
			break;
		    }
		}
	    } else {
		col = loc->snps[0]->col;

		if (t->nucs[col].allele_cnt == 2)
		    ordered_loci.push_back(GenPos(loc->id, 0, loc->sort_bp(col), snp));
	    }
	}
	sort(ordered_loci.begin(), ordered_loci.end(), compare_genpos);

	//
	// Output the total number of SNP sites and the number of individuals.
	//
	fh << samples.size()      << "\n"
	   << ordered_loci.size() << "\n";

	//
	// Output the position of each site according to its basepair.
	//
	fh << "P";
    	for (uint pos = 0; pos < ordered_loci.size(); pos++)
	    fh << " " << ordered_loci[pos].bp;
	fh << "\n";

	//
	// Output a line of 'S' characters for SNP markers, 'M' characters for multiallelic haplotypes.
	//
    	for (uint pos = 0; pos < ordered_loci.size(); pos++) {
	    if (pos > 0) fh << " ";
	    fh << (ordered_loci[pos].type == snp ? "S" : "M");
	}
	fh << "\n";

	//
	// Now output each sample name followed by a new line, then all of the genotypes for that sample
	// on two lines.
	//

	map<int, pair<int, int> >::iterator pit;
	string       gtypes_str;
	bool         found;
	int          start_index, end_index, pop_id;
	char         p_allele, q_allele;
	stringstream gtypes;

	for (pit = pop_indexes.begin(); pit != pop_indexes.end(); pit++) {
	    pop_id      = psum->pop_index(pit->first);
	    start_index = pit->second.first;
	    end_index   = pit->second.second;

	    for (int j = start_index; j <= end_index; j++) {
		//
		// Output all the loci for this sample, printing only the p allele
		//
		fh << samples[pmap->rev_sample_index(j)] << "\n";
		
		gtypes.str("");
		for (uint pos = 0; pos < ordered_loci.size(); pos++) {
		    loc = catalog[ordered_loci[pos].id];

		    s = psum->locus(loc->id);
		    d = pmap->locus(loc->id);
		    t = psum->locus_tally(loc->id);

		    //
		    // Will we output this locus as a haplotype or as a SNP?
		    //
		    if (ordered_loci[pos].type == haplotype) {
			if (d[j] == NULL) {
			    //
			    // Data does not exist.
			    //
			    gtypes << "-1 ";
			} else {
			    //
			    // Data exists, output the first haplotype. We will assume the haplotypes are 
			    // numbered by their position in the loc->strings vector.
			    //
			    if (d[j]->obshap.size() > 2)  {
				// cerr << "Warning: too many haplotypes, catalog locus: " << loc->id << "\n";
				gtypes << "-1 ";
			    } else {
				found = false;
				for (uint k = 0; k < loc->strings.size(); k++)
				    if (d[j]->obshap[0] == loc->strings[k].first) {
					found = true;
					gtypes << k + 1 << " ";
				    }
				if (found == false)
				    cerr << "Unable to find haplotype " << d[j]->obshap[0] << " from individual " 
					 << samples[pmap->rev_sample_index(j)] << "; catalog locus: " << loc->id << "\n";
			    }
			}
		    } else {
			col = loc->snps[ordered_loci[pos].snp_index]->col;

			if (s[pop_id]->nucs[col].incompatible_site ||
			    s[pop_id]->nucs[col].filtered_site) {
			    //
			    // This site contains more than two alleles in this population or was filtered
			    // due to a minor allele frequency that is too low.
			    //
			    gtypes << "? ";

			} else if (d[j] == NULL) {
			    //
			    // Data does not exist.
			    //
			    gtypes << "? ";
			} else if (d[j]->model[col] == 'U') {
			    //
			    // Data exists, but the model call was uncertain.
			    //
			    gtypes << "? ";
			} else {
			    //
			    // Tally up the nucleotide calls.
			    //
			    tally_observed_haplotypes(d[j]->obshap, ordered_loci[pos].snp_index, p_allele, q_allele);

			    if (p_allele == 0 && q_allele == 0)
				gtypes << "? ";
			    else if (p_allele == 0)
				gtypes << q_allele << " ";
			    else
				gtypes << p_allele << " ";
			}
		    }
		}
		gtypes_str = gtypes.str();
		fh << gtypes_str.substr(0, gtypes_str.length() - 1) << "\n";

		//
		// Output all the loci for this sample again, now for the q allele
		//
		gtypes.str("");
		for (uint pos = 0; pos < ordered_loci.size(); pos++) {
		    loc = catalog[ordered_loci[pos].id];

		    s = psum->locus(loc->id);
		    d = pmap->locus(loc->id);
		    t = psum->locus_tally(loc->id);

		    //
		    // Will we output this locus as a haplotype or as a SNP?
		    //
		    if (ordered_loci[pos].type == haplotype) {
			if (d[j] == NULL) {
			    //
			    // Data does not exist.
			    //
			    gtypes << "-1 ";
			} else {
			    //
			    // Data exists, output the second haplotype. We will assume the haplotypes are 
			    // numbered by their position in the loc->strings vector.
			    //
			    if (d[j]->obshap.size() > 2)  {
				// cerr << "Warning: too many haplotypes, catalog locus: " << loc->id << "\n";
				gtypes << "-1 ";
			    } else if (d[j]->obshap.size() > 1)  {
				found = false;
				for (uint k = 0; k < loc->strings.size(); k++)
				    if (d[j]->obshap[1] == loc->strings[k].first) {
					found = true;
					gtypes << k + 1 << " ";
				    }
				if (found == false)
				    cerr << "Unable to find haplotype " << d[j]->obshap[1] << " from individual " 
					 << samples[pmap->rev_sample_index(j)] << "; catalog locus: " << loc->id << "\n";
			    } else {
				found = false;
				for (uint k = 0; k < loc->strings.size(); k++)
				    if (d[j]->obshap[0] == loc->strings[k].first) {
					found = true;
					gtypes << k + 1 << " ";
				    }
				if (found == false)
				    cerr << "Unable to find haplotype " << d[j]->obshap[0] << " from individual " 
					 << samples[pmap->rev_sample_index(j)] << "; catalog locus: " << loc->id << "\n";
			    }
			}
		    } else {
			col = loc->snps[ordered_loci[pos].snp_index]->col;

			if (s[pop_id]->nucs[col].incompatible_site ||
			    s[pop_id]->nucs[col].filtered_site) {
			    gtypes << "? ";

			} else if (d[j] == NULL) {
			    gtypes << "? ";

			} else if (d[j]->model[col] == 'U') {
			    gtypes << "? ";

			} else {
			    tally_observed_haplotypes(d[j]->obshap, ordered_loci[pos].snp_index, p_allele, q_allele);

			    if (p_allele == 0 && q_allele == 0)
				gtypes << "? ";
			    else if (q_allele == 0)
				gtypes << p_allele << " ";
			    else
				gtypes << q_allele << " ";
			}
		    }
		}
		gtypes_str = gtypes.str();
		fh << gtypes_str.substr(0, gtypes_str.length() - 1) << "\n";
	    }
	}

	fh.close();
    }

    cerr << "done.\n";

    return 0;
}

int 
write_plink(map<int, CSLocus *> &catalog, 
	    PopMap<CSLocus> *pmap, 
	    PopSum<CSLocus> *psum, 
	    map<int, pair<int, int> > &pop_indexes, 
	    map<int, string> &samples) 
{
    //
    // Write a PLINK file as defined here: http://pngu.mgh.harvard.edu/~purcell/plink/data.shtml
    //
    // We will write one file per chromosome.
    //
    cerr << "Writing population data to PLINK files...";

    //
    // Obtain the current date.
    //
    time_t     rawtime;
    struct tm *timeinfo;
    char       date[32];
    time(&rawtime);
    timeinfo = localtime(&rawtime);
    strftime(date, 32, "%B %d, %Y", timeinfo);

    map<string, vector<CSLocus *> >::iterator it;
    CSLocus  *loc;
    Datum   **d;
    LocSum  **s;
    LocTally *t;
    string    chr;

    //
    // First, write a markers file containing each marker, the chromosome it falls on,
    // an empty centiMorgan field, and finally its genomic position in basepairs.
    //
    stringstream pop_name;
    pop_name << "batch_" << batch_id << ".plink.map";
    string file = in_path + pop_name.str();

    ofstream fh(file.c_str(), ofstream::out);

    if (fh.fail()) {
	cerr << "Error opening PLINK markers file '" << file << "'\n";
	exit(1);
    }

    //
    // Output the header.
    //
    fh << "# Stacks v" << VERSION << "; " << " PLINK v1.07; " << date << "\n";

    for (it = pmap->ordered_loci.begin(); it != pmap->ordered_loci.end(); it++) {
	chr = it->first;

    	for (uint pos = 0; pos < it->second.size(); pos++) {
    	    loc = it->second[pos];
	    t   = psum->locus_tally(loc->id);

    	    for (uint i = 0; i < loc->snps.size(); i++) {
    		uint col = loc->snps[i]->col;
		if (t->nucs[col].allele_cnt == 2)
		    fh << chr << "\t"
		       << loc->id << "_" << col << "\t"
		       << "0\t" 
		       << loc->sort_bp(col) << "\n";
	    }
	}
    }
    fh.close();

    //
    // Now output the genotypes in a separate file.
    //
    pop_name.str("");
    pop_name << "batch_" << batch_id << ".plink.ped";
    file = in_path + pop_name.str();

    fh.open(file.c_str(), ofstream::out);

    if (fh.fail()) {
	cerr << "Error opening PLINK markers file '" << file << "'\n";
	exit(1);
    }

    fh << "# Stacks v" << VERSION << "; " << " PLINK v1.07; " << date << "\n";

    map<int, pair<int, int> >::iterator pit;
    int  start_index, end_index, pop_id;
    char p_allele, q_allele;

    //
    //  marker, output the genotypes for each sample in two successive columns.
    //
    for (pit = pop_indexes.begin(); pit != pop_indexes.end(); pit++) {
	pop_id      = psum->pop_index(pit->first);
	start_index = pit->second.first;
	end_index   = pit->second.second;

	for (int j = start_index; j <= end_index; j++) {

	    fh << pit->first << "\t"
	       << samples[pmap->rev_sample_index(j)] << "\t"
	       << "0\t"  // Paternal ID
	       << "0\t"  // Maternal ID
	       << "0\t"  // Sex
	       << "0";   // Phenotype

	    for (it = pmap->ordered_loci.begin(); it != pmap->ordered_loci.end(); it++) {
		for (uint pos = 0; pos < it->second.size(); pos++) {
		    loc = it->second[pos];

		    s = psum->locus(loc->id);
		    d = pmap->locus(loc->id);
		    t = psum->locus_tally(loc->id);
 
		    for (uint i = 0; i < loc->snps.size(); i++) {
			uint col = loc->snps[i]->col;

			// 
			// If this site is fixed in all populations or has too many alleles don't output it.
			//
			if (t->nucs[col].allele_cnt != 2) 
			    continue;
			//
			// Output the p and q alleles
			//
			if (s[pop_id]->nucs[col].incompatible_site ||
			    s[pop_id]->nucs[col].filtered_site) {
			    //
			    // This site contains more than two alleles in this population or was filtered
			    // due to a minor allele frequency that is too low.
			    //
			    fh << "\t" << "0" << "\t" << "0";
			} else if (d[j] == NULL) {
			    //
			    // Data does not exist.
			    //
			    fh << "\t" << "0" << "\t" << "0";
			} else if (d[j]->model[col] == 'U') {
			    //
			    // Data exists, but the model call was uncertain.
			    //
			    fh << "\t" << "0" << "\t" << "0";
			} else {
			    //
			    // Tally up the nucleotide calls.
			    //
			    tally_observed_haplotypes(d[j]->obshap, i, p_allele, q_allele);

			    if (p_allele == 0 && q_allele == 0)
				fh << "\t" << "0" << "\t" << "0";
			    else if (p_allele == 0)
				fh << "\t" << q_allele << "\t" << q_allele;
			    else if (q_allele == 0)
				fh << "\t" << p_allele << "\t" << p_allele;
			    else
				fh << "\t" << p_allele << "\t" << q_allele;
			}
		    }
		}
	    }
	    fh << "\n";
	}
    }

    fh.close();

    cerr << "done.\n";

    return 0;
}

int 
write_beagle(map<int, CSLocus *> &catalog, 
	    PopMap<CSLocus> *pmap, 
	    PopSum<CSLocus> *psum, 
	    map<int, pair<int, int> > &pop_indexes, 
	    map<int, string> &samples) 
{
    //
    // Write a Beagle file as defined here: http://faculty.washington.edu/browning/beagle/beagle.html
    //
    // We will write one file per chromosome, per population.
    //
    cerr << "Writing population data to unphased Beagle files...";

    //
    // Obtain the current date.
    //
    time_t     rawtime;
    struct tm *timeinfo;
    char       date[32];
    time(&rawtime);
    timeinfo = localtime(&rawtime);
    strftime(date, 32, "%B %d, %Y", timeinfo);

    map<string, vector<CSLocus *> >::iterator it;
    CSLocus  *loc;
    Datum   **d;
    LocSum  **s;
    LocTally *t;
    uint      col;

    stringstream pop_name;
    string       file;

    for (it = pmap->ordered_loci.begin(); it != pmap->ordered_loci.end(); it++) {

	//
	// We need to determine an ordering that can take into account overlapping RAD sites.
	//
	vector<GenPos> ordered_loci;
	for (uint pos = 0; pos < it->second.size(); pos++) {
	    loc = it->second[pos];
	    t   = psum->locus_tally(loc->id);

	    for (uint i = 0; i < loc->snps.size(); i++) {
		col = loc->snps[i]->col;
		if (t->nucs[col].allele_cnt == 2)
		    ordered_loci.push_back(GenPos(loc->id, i, loc->sort_bp(col)));
	    }
	}
	sort(ordered_loci.begin(), ordered_loci.end(), compare_genpos);

	//
	// Now output the genotypes in a separate file for each population.
	//
	map<int, pair<int, int> >::iterator pit;
	int  start_index, end_index, pop_id;

	for (pit = pop_indexes.begin(); pit != pop_indexes.end(); pit++) {
	    pop_id      = psum->pop_index(pit->first);
	    start_index = pit->second.first;
	    end_index   = pit->second.second;

	    //
	    // Open a markers file containing each marker, its genomic position in basepairs
	    // and the two alternative alleles at this position.
	    //
	    pop_name.str("");
	    pop_name << "batch_" << batch_id << "." << pop_key[pit->first] << "-" << it->first << ".unphased.bgl.markers";
	    file = in_path + pop_name.str();

	    ofstream mfh(file.c_str(), ofstream::out);
	    if (mfh.fail()) {
		cerr << "Error opening Beagle markers file '" << file << "'\n";
		exit(1);
	    }
	    mfh << "# Stacks v" << VERSION << "; " << " Beagle v3.3; " << date << "\n";

	    //
	    // Open the genotypes file.
	    //
	    pop_name.str("");
	    pop_name << "batch_" << batch_id << "." << pop_key[pit->first] << "-" << it->first << ".unphased.bgl";
	    file = in_path + pop_name.str();

	    ofstream fh(file.c_str(), ofstream::out);
	    if (fh.fail()) {
		cerr << "Error opening Beagle genotypes file '" << file << "'\n";
		exit(1);
	    }
	    fh << "# Stacks v" << VERSION << "; " << " Beagle v3.3; " << date << "\n";

	    char p_allele, q_allele;
	    //
	    // Output a list of all the samples in this population.
	    //
	    fh << "I\tid";
	    for (int j = start_index; j <= end_index; j++)
		fh << "\t" << samples[pmap->rev_sample_index(j)] << "\t" << samples[pmap->rev_sample_index(j)];
	    fh << "\n";

	    //
	    // Output population IDs for each sample.
	    //
	    fh << "S\tid";
	    for (int j = start_index; j <= end_index; j++)
		fh << "\t" << pit->first << "\t" << pit->first;
	    fh << "\n";

	    //
	    // For each marker, output the genotypes for each sample in two successive columns.
	    //
	    for (uint pos = 0; pos < ordered_loci.size(); pos++) {
		loc = catalog[ordered_loci[pos].id];

		s   = psum->locus(loc->id);
		d   = pmap->locus(loc->id);
		t   = psum->locus_tally(loc->id);
		col = loc->snps[ordered_loci[pos].snp_index]->col;

		// 
		// If this site is fixed in all populations or has too many alleles don't output it.
		//
		if (t->nucs[col].allele_cnt != 2) 
		    continue;

		//
		// If this site is monomorphic in this population don't output it.
		//
		if (s[pop_id]->nucs[col].pi == 0.0)
		    continue;

		//
		// Output this locus to the markers file.
		//
		mfh << loc->id << "_" << col << "\t" 
		    << loc->sort_bp(col)     << "\t" 
		    << t->nucs[col].p_allele << "\t" 
		    << t->nucs[col].q_allele << "\n";

		fh << "M" << "\t" << loc->id << "_" << col;

		for (int j = start_index; j <= end_index; j++) {
		    //
		    // Output the p allele
		    //
		    if (s[pop_id]->nucs[col].incompatible_site ||
			s[pop_id]->nucs[col].filtered_site) {
			//
			// This site contains more than two alleles in this population or was filtered
			// due to a minor allele frequency that is too low.
			//
			fh << "\t" << "?";

		    } else if (d[j] == NULL) {
			//
			// Data does not exist.
			//
			fh << "\t" << "?";
		    } else if (d[j]->model[col] == 'U') {
			//
			// Data exists, but the model call was uncertain.
			//
			fh << "\t" << "?";
		    } else {
			//
			// Tally up the nucleotide calls.
			//
			tally_observed_haplotypes(d[j]->obshap, ordered_loci[pos].snp_index, p_allele, q_allele);

			if (p_allele == 0 && q_allele == 0)
			    fh << "\t" << "?";
			else if (p_allele == 0)
			    fh << "\t" << q_allele;
			else
			    fh << "\t" << p_allele;
		    }

		    //
		    // Now output the q allele
		    //
		    if (s[pop_id]->nucs[col].incompatible_site ||
			s[pop_id]->nucs[col].filtered_site) {
			fh << "\t" << "?";

		    } else if (d[j] == NULL) {
			fh << "\t" << "?";

		    } else if (d[j]->model[col] == 'U') {
			fh << "\t" << "?";

		    } else {
			if (p_allele == 0 && q_allele == 0)
			    fh << "\t" << "?";
			else if (q_allele == 0)
			    fh << "\t" << p_allele;
			else
			    fh << "\t" << q_allele;
		    }
		}
		fh << "\n";
	    }

	    fh.close();
	    mfh.close();
	}
    }

    cerr << "done.\n";

    return 0;
}

int 
write_beagle_phased(map<int, CSLocus *> &catalog, 
		    PopMap<CSLocus> *pmap, 
		    PopSum<CSLocus> *psum, 
		    map<int, pair<int, int> > &pop_indexes, 
		    map<int, string> &samples) 
{
    //
    // Write a Beagle file as a set of haplotpyes as defined here: 
    //   http://faculty.washington.edu/browning/beagle/beagle.html
    //
    // We will write one file per chromosome.
    //
    cerr << "Writing population data to phased Beagle files...";

    //
    // Obtain the current date.
    //
    time_t     rawtime;
    struct tm *timeinfo;
    char       date[32];
    time(&rawtime);
    timeinfo = localtime(&rawtime);
    strftime(date, 32, "%B %d, %Y", timeinfo);

    map<string, vector<CSLocus *> >::iterator it;
    CSLocus  *loc;
    Datum   **d;

    stringstream pop_name;
    string       file;

    for (it = pmap->ordered_loci.begin(); it != pmap->ordered_loci.end(); it++) {

	//
	// We need to determine an ordering for all legitimate loci/SNPs.
	//
	vector<GenPos> ordered_loci;
    	for (uint pos = 0; pos < it->second.size(); pos++) {
    	    loc = it->second[pos];

	    if (loc->snps.size() == 0) continue;

	    //
	    // Check that there aren't too many haplotypes (PHASE has a max of 50).
	    //
	    if (loc->alleles.size() > 40) continue;

	    //
	    // Iterate over the population to determine that this subset of the population
	    // has data at this locus.
	    //
	    d = pmap->locus(loc->id);
	    for (int j = 0; j < pmap->sample_cnt(); j++) {
		if (d[j] != NULL &&
		    d[j]->obshap.size() > 0 && 
		    d[j]->obshap.size() <= 2) {
		    //
		    // Data exists, and their are the corrent number of haplotypes.
		    //
		    ordered_loci.push_back(GenPos(loc->id, 0, loc->sort_bp(), haplotype));
		    break;
		}
	    }
	}
	sort(ordered_loci.begin(), ordered_loci.end(), compare_genpos);

	//
	// Now output the genotypes in a separate file for each population.
	//
	map<int, pair<int, int> >::iterator pit;
	int  start_index, end_index, pop_id;

	for (pit = pop_indexes.begin(); pit != pop_indexes.end(); pit++) {
	    pop_id      = psum->pop_index(pit->first);
	    start_index = pit->second.first;
	    end_index   = pit->second.second;

	    //
	    // Open a file for writing the markers: their genomic position in basepairs
	    // and the two alternative alleles at this position.
	    //
	    pop_name.str("");
	    pop_name << "batch_" << batch_id << "." << pop_key[pit->first] << "-" << it->first << ".phased.bgl.markers";
	    file = in_path + pop_name.str();

	    ofstream mfh(file.c_str(), ofstream::out);
	    if (mfh.fail()) {
		cerr << "Error opening Beagle markers file '" << file << "'\n";
		exit(1);
	    }
	    mfh << "# Stacks v" << VERSION << "; " << " Beagle v3.3; " << date << "\n";

	    //
	    // Now output the haplotypes in a separate file.
	    //
	    pop_name.str("");
	    pop_name << "batch_" << batch_id << "." << pop_key[pit->first] << "-" << it->first << ".phased.bgl";
	    file = in_path + pop_name.str();

	    ofstream fh(file.c_str(), ofstream::out);
	    if (fh.fail()) {
		cerr << "Error opening Beagle markers file '" << file << "'\n";
		exit(1);
	    }
	    fh << "# Stacks v" << VERSION << "; " << " Beagle v3.3; " << date << "\n";

	    //
	    // Output a list of all the samples in the data set.
	    //
	    fh << "I\tid";
	    for (int j = start_index; j <= end_index; j++)
		fh << "\t" << samples[pmap->rev_sample_index(j)] << "\t" << samples[pmap->rev_sample_index(j)];
	    fh << "\n";

	    //
	    // Output population IDs for each sample.
	    //
	    fh << "S\tid";
	    for (int j = start_index; j <= end_index; j++)
		fh << "\t" << pop_id << "\t" << pop_id;
	    fh << "\n";

	    for (uint pos = 0; pos < ordered_loci.size(); pos++) {
		loc = catalog[ordered_loci[pos].id];
		d   = pmap->locus(loc->id);

		//
		// If this locus is monomorphic in this population don't output it.
		//
		set<string> haplotypes;
		for (int j = start_index; j <= end_index; j++) {
		    if (d[j] == NULL) continue;

		    if (d[j]->obshap.size() == 2) {
			haplotypes.insert(d[j]->obshap[0]);
			haplotypes.insert(d[j]->obshap[1]);
		    } else {
			haplotypes.insert(d[j]->obshap[0]);
		    }
		}
		if (haplotypes.size() == 1) continue;

		//
		// Output this locus to the markers file.
		//
		mfh << loc->id << "\t" 
		    << loc->sort_bp();
		for (uint j = 0; j < loc->strings.size(); j++)
		    mfh << "\t" << loc->strings[j].first;
		mfh << "\n";

		//
		// For each marker, output the genotypes for each sample in two successive columns.
		//
		fh << "M" << "\t" << loc->id;

		for (int j = start_index; j <= end_index; j++) {
		    //
		    // Output the p and the q haplotype
		    //
		    if (d[j] == NULL) {
			//
			// Data does not exist.
			//
			fh << "\t" << "?" << "\t" << "?";
		    } else {
			//
			// Data exists, output the first haplotype. We will assume the haplotypes are 
			// numbered by their position in the loc->strings vector.
			//
			if (d[j]->obshap.size() > 2)
			    fh << "\t" << "?" << "\t" << "?";
			else if (d[j]->obshap.size() == 2)
			    fh << "\t" << d[j]->obshap[0] << "\t" << d[j]->obshap[1];
			else
			    fh << "\t" << d[j]->obshap[0] << "\t" << d[j]->obshap[0];
		    }
		}
		fh << "\n";
	    }
	    fh.close();
	    mfh.close();
	}
    }   

    cerr << "done.\n";

    return 0;
}

int 
write_phylip(map<int, CSLocus *> &catalog, 
	     PopMap<CSLocus> *pmap, 
	     PopSum<CSLocus> *psum, 
	     map<int, pair<int, int> > &pop_indexes, 
	     map<int, string> &samples) 
{
    //
    // We want to find loci where each locus is fixed within a population but variable between populations.
    //
    // We will write those loci to a Phylip file as defined here: 
    //     http://evolution.genetics.washington.edu/phylip/doc/main.html#inputfiles
    //
    stringstream pop_name;
    pop_name << "batch_" << batch_id << ".phylip";
    string file = in_path + pop_name.str();

    cerr << "Writing population data to Phylip file '" << file << "'; ";

    ofstream fh(file.c_str(), ofstream::out);

    if (fh.fail()) {
        cerr << "Error opening Phylip file '" << file << "'\n";
    	exit(1);
    }

    pop_name << ".log";
    file = in_path + pop_name.str();

    cerr << "logging nucleotide positions to '" << file << "'...";

    ofstream log_fh(file.c_str(), ofstream::out);

    if (log_fh.fail()) {
        cerr << "Error opening Phylip Log file '" << file << "'\n";
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
    strftime(date, 32, "%B %d, %Y", timeinfo);

    log_fh << "# Stacks v" << VERSION << "; " << " Phylip sequential; " << date << "\n"
	   << "# Seq Pos\tLocus ID\tColumn\tPopulation\n";

    map<string, vector<CSLocus *> >::iterator it;
    CSLocus  *loc;
    LocSum  **s;
    LocTally *t;

    map<int, pair<int, int> >::iterator pit;
    int  pop_cnt = psum->pop_cnt();
    int  pop_id;
    char nuc;

    //
    // A map storing, for each population, the concatenated list of interspecific nucleotides.
    //
    map<int, string> interspecific_nucs;

    int index = 0;
    for (it = pmap->ordered_loci.begin(); it != pmap->ordered_loci.end(); it++) {
	for (uint pos = 0; pos < it->second.size(); pos++) {
	    loc = it->second[pos];

 	    s = psum->locus(loc->id);
	    t = psum->locus_tally(loc->id);

	    for (uint i = 0; i < loc->snps.size(); i++) {
		uint col = loc->snps[i]->col;

		if (phylip_var == false) {
		    //
		    // We are looking for loci that are fixed within each population, but are 
		    // variable between one or more populations.
		    //
		    if (t->nucs[col].fixed == true || t->nucs[col].allele_cnt != 2 || t->nucs[col].pop_cnt < 2)
			continue;

		    bool fixed_within = true;
		    for (int j = 0; j < pop_cnt; j++) {
			if (s[j]->nucs[col].num_indv == 0)
			    continue;
			if (s[j]->nucs[col].fixed == false) {
			    fixed_within = false;
			    break;
			}
		    }
		    if (fixed_within == false) continue;

		    log_fh << index << "\t" << loc->id << "\t" << col << "\t";

		    for (int j = 0; j < pop_cnt; j++) {
			pop_id = psum->rev_pop_index(j);

			if (s[j]->nucs[col].num_indv > 0) {
			    interspecific_nucs[pop_id] += s[j]->nucs[col].p_nuc;
			    log_fh << pop_key[pop_id] << ":" << s[j]->nucs[col].p_nuc << ",";
			} else {
			    interspecific_nucs[pop_id] += 'N';
			    log_fh << pop_key[pop_id] << ":N" << ",";
			}
		    }
		    log_fh << "\n";
		    index++;

		} else {
		    //
		    // Encode SNPs that are variable within a population as well, using IUPAC notation:
		    //     http://en.wikipedia.org/wiki/Nucleic_acid_notation#IUPAC_notation
		    //
		    if (t->nucs[col].allele_cnt != 2)
			continue;

		    log_fh << index << "\t" << loc->id << "\t" << col << "\t";

		    for (int j = 0; j < pop_cnt; j++) {
			pop_id = psum->rev_pop_index(j);

			switch(s[j]->nucs[col].p_nuc) {
			case 0:
			    nuc = 'N';
			    break;
			case 'A':
			    switch(s[j]->nucs[col].q_nuc) {
			    case 'C':
				nuc = 'M';
				break;
			    case 'G':
				nuc = 'R';
				break;
			    case 'T':
				nuc = 'W';
				break;
			    case 0:
				nuc = 'A';
				break;
			    }
			    break;
			case 'C':
			    switch(s[j]->nucs[col].q_nuc) {
			    case 'A':
				nuc = 'M';
				break;
			    case 'G':
				nuc = 'S';
				break;
			    case 'T':
				nuc = 'Y';
				break;
			    case 0:
				nuc = 'C';
				break;
			    }
			    break;
			case 'G':
			    switch(s[j]->nucs[col].q_nuc) {
			    case 'A':
				nuc = 'R';
				break;
			    case 'C':
				nuc = 'S';
				break;
			    case 'T':
				nuc = 'K';
				break;
			    case 0:
				nuc = 'G';
				break;
			    }
			    break;
			case 'T':
			    switch(s[j]->nucs[col].q_nuc) {
			    case 'A':
				nuc = 'W';
				break;
			    case 'C':
				nuc = 'Y';
				break;
			    case 'G':
				nuc = 'K';
				break;
			    case 0:
				nuc = 'T';
				break;
			    }
			    break;
			}
			interspecific_nucs[pop_id] += nuc;
			log_fh << pop_key[pop_id] << ":" << nuc << ",";

		    }
		    log_fh << "\n";
		    index++;
		}
	    }
	}
    }

    if (interspecific_nucs.size() == 0) {
	cerr << "  No data is available to write to the Phylip file.\n";
	return 0;
    }

    char id_str[id_len];
    uint len;

    fh << pop_indexes.size() << "    " << interspecific_nucs.begin()->second.length() << "\n";
    for (pit = pop_indexes.begin(); pit != pop_indexes.end(); pit++) {
	pop_id = pit->first;

	sprintf(id_str, "%s", pop_key[pop_id].c_str());
	len = strlen(id_str);
	for (uint j = len; j < 10; j++)
	    id_str[j] = ' ';
	id_str[9] = '\0'; 

	fh << id_str << " " << interspecific_nucs[pop_id] << "\n";
    }

    //
    // Output the header.
    //
    fh << "# Stacks v" << VERSION << "; " << " Phylip sequential; " << date << "\n";

    fh.close();
    log_fh.close();

    cerr << "done.\n";

    return 0;
}

int 
write_fullseq_phylip(map<int, CSLocus *> &catalog, 
		     PopMap<CSLocus> *pmap, 
		     PopSum<CSLocus> *psum, 
		     map<int, pair<int, int> > &pop_indexes, 
		     map<int, string> &samples) 
{
    //
    // We want to write all variable loci in Phylip interleaved format. Polymorphic positions
    // will be encoded using IUPAC notation.
    //
    // We will write those loci to a Phylip file as defined here: 
    //     http://evolution.genetics.washington.edu/phylip/doc/main.html#inputfiles
    //
    stringstream pop_name;
    pop_name << "batch_" << batch_id << ".fullseq.phylip";
    string file = in_path + pop_name.str();

    cerr << "Writing full sequence population data to Phylip file '" << file << "'; ";

    ofstream fh(file.c_str(), ofstream::out);

    if (fh.fail()) {
        cerr << "Error opening Phylip file '" << file << "'\n";
    	exit(1);
    }

    //
    // We will also write a file that allows us to specify each RAD locus as a separate partition
    // for use in phylogenetics programs.
    //
    pop_name.str("");
    pop_name << "batch_" << batch_id << ".fullseq.partitions.phylip";
    file = in_path + pop_name.str();

    ofstream par_fh(file.c_str(), ofstream::out);

    if (par_fh.fail()) {
        cerr << "Error opening Phylip partitions file '" << file << "'\n";
    	exit(1);
    }

    pop_name.str("");
    pop_name << "batch_" << batch_id << "fullseq.phylip.log";
    file = in_path + pop_name.str();

    cerr << "logging nucleotide positions to '" << file << "'...";

    ofstream log_fh(file.c_str(), ofstream::out);

    if (log_fh.fail()) {
        cerr << "Error opening Phylip Log file '" << file << "'\n";
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
    strftime(date, 32, "%B %d, %Y", timeinfo);

    log_fh << "# Stacks v" << VERSION << "; " << " Phylip interleaved; " << date << "\n"
	   << "# Locus ID\tLine Number";
    if (loci_ordered) log_fh << "\tChr\tBasepair";
    log_fh << "\n";

    map<string, vector<CSLocus *> >::iterator it;
    CSLocus  *loc;
    LocSum  **s;
    LocTally *t;

    map<int, pair<int, int> >::iterator pit;
    int  pop_cnt = psum->pop_cnt();
    int  pop_id;
    char nuc;

    bool include;
    char id_str[id_len];
    uint len = 0;

    //
    // Determine the length of sequence we will output.
    //
    for (it = pmap->ordered_loci.begin(); it != pmap->ordered_loci.end(); it++) {
	for (uint pos = 0; pos < it->second.size(); pos++) {
	    loc = it->second[pos];

 	    t = psum->locus_tally(loc->id);

	    include = true;
	    for (uint i = 0; i < loc->snps.size(); i++) {
		uint col = loc->snps[i]->col;

		if (t->nucs[col].allele_cnt != 2)
		    include = false;
	    }

	    if (include)
		len += strlen(loc->con);
	}
    }

    map<int, string> outstrs;
    
    fh << pop_indexes.size() << "    " << len << "\n";

    for (pit = pop_indexes.begin(); pit != pop_indexes.end(); pit++) {
	pop_id = pit->first;

	outstrs[pop_id] = "";
	sprintf(id_str, "%s", pop_key[pop_id].c_str());
	len = strlen(id_str);
	for (uint j = len; j < 10; j++)
	    id_str[j] = ' ';
	id_str[9] = '\0'; 

	outstrs[pop_id] += string(id_str) + " ";
    }

    char *seq;
    int   line  = 1;
    int   index = 1;
    int   cnt   = 1;
    
    for (it = pmap->ordered_loci.begin(); it != pmap->ordered_loci.end(); it++) {
	for (uint pos = 0; pos < it->second.size(); pos++) {
	    loc = it->second[pos];

 	    s = psum->locus(loc->id);
	    t = psum->locus_tally(loc->id);

	    include = true;
	    for (uint i = 0; i < loc->snps.size(); i++) {
		uint col = loc->snps[i]->col;

		if (t->nucs[col].allele_cnt != 2)
		    include = false;
	    }

	    if (!include)
		continue;

	    seq = new char[loc->len + 1];
	    strcpy(seq, loc->con);

	    for (int j = 0; j < pop_cnt; j++) {
		pop_id = psum->rev_pop_index(j);

		for (uint i = 0; i < loc->snps.size(); i++) {
		    uint col = loc->snps[i]->col;

		    //
		    // Encode SNPs that are variable within a population using IUPAC notation:
		    //     http://en.wikipedia.org/wiki/Nucleic_acid_notation#IUPAC_notation
		    //
		    switch(s[j]->nucs[col].p_nuc) {
		    case 0:
			nuc = 'N';
			break;
		    case 'A':
			switch(s[j]->nucs[col].q_nuc) {
			case 'C':
			    nuc = 'M';
			    break;
			case 'G':
			    nuc = 'R';
			    break;
			case 'T':
			    nuc = 'W';
			    break;
			case 0:
			    nuc = 'A';
			    break;
			}
			break;
		    case 'C':
			switch(s[j]->nucs[col].q_nuc) {
			case 'A':
			    nuc = 'M';
			    break;
			case 'G':
			    nuc = 'S';
			    break;
			case 'T':
			    nuc = 'Y';
			    break;
			case 0:
			    nuc = 'C';
			    break;
			}
			break;
		    case 'G':
			switch(s[j]->nucs[col].q_nuc) {
			case 'A':
			    nuc = 'R';
			    break;
			case 'C':
			    nuc = 'S';
			    break;
			case 'T':
			    nuc = 'K';
			    break;
			case 0:
			    nuc = 'G';
			    break;
			}
			break;
		    case 'T':
			switch(s[j]->nucs[col].q_nuc) {
			case 'A':
			    nuc = 'W';
			    break;
			case 'C':
			    nuc = 'Y';
			    break;
			case 'G':
			    nuc = 'K';
			    break;
			case 0:
			    nuc = 'T';
			    break;
			}
			break;
		    }

		    seq[col] = nuc;
		}
		
		outstrs[pop_id] += string(seq);
	    }
	    delete [] seq;

	    log_fh << line << "\t" << loc->id;
	    if (loci_ordered) log_fh << "\t" << loc->loc.chr << "\t" << loc->sort_bp() + 1;
	    log_fh << "\n";
	    
	    for (pit = pop_indexes.begin(); pit != pop_indexes.end(); pit++) {
		pop_id = pit->first;
		fh << outstrs[pop_id] << "\n";
		outstrs[pop_id] = "";
		line++;
	    }
	    fh << "\n";
	    line++;

	    par_fh << "DNA, p" << cnt << "=" << index << "-" << index + loc->len - 1 << "\n";
	    index += loc->len;
	    cnt++;
	}
    }

    //
    // Output the header.
    //
    fh << "# Stacks v" << VERSION << "; " << " Phylip interleaved; " << date << "\n";

    fh.close();
    par_fh.close();
    log_fh.close();

    cerr << "done.\n";

    return 0;
}

int 
tally_ref_alleles(LocSum **s, int pop_cnt, int snp_index, char &p_allele, char &q_allele) 
{
    int  nucs[4] = {0};
    char nuc[2];

    for (int j = 0; j < pop_cnt; j++) {
	nuc[0] = 0;
	nuc[1] = 0;
        nuc[0] = s[j]->nucs[snp_index].p_nuc;
        nuc[1] = s[j]->nucs[snp_index].q_nuc;

	for (uint k = 0; k < 2; k++) 
	    switch(nuc[k]) {
	    case 'A':
	    case 'a':
		nucs[0]++;
		break;
	    case 'C':
	    case 'c':
		nucs[1]++;
		break;
	    case 'G':
	    case 'g':
		nucs[2]++;
		break;
	    case 'T':
	    case 't':
		nucs[3]++;
		break;
	    }
    }

    //
    // Determine how many alleles are present at this position in this population.
    // We cannot deal with more than two alternative alleles, if there are more than two
    // in a single population, print a warning and exclude this nucleotide position.
    //
    int i;
    int allele_cnt = 0;
    for (i = 0; i < 4; i++)
	if (nucs[i] > 0) allele_cnt++;

    if (allele_cnt > 2) {
	p_allele = 0;
	q_allele = 0;
	return 0;
    }

    //
    // Record which nucleotide is the P allele and which is the Q allele.
    //
    p_allele = 0;
    q_allele = 0;

    i = 0;
    while (p_allele == 0 && i < 4) {
	if (nucs[i] > 0) {
	    switch(i) {
	    case 0:
		p_allele = 'A';
		break;
	    case 1:
		p_allele = 'C';
		break;
	    case 2:
		p_allele = 'G';
		break;
	    case 3:
		p_allele = 'T';
		break;
	    }
	}
	i++;
    }
    while (q_allele == 0 && i < 4) {
	if (nucs[i] > 0) {
	    switch(i) {
	    case 1:
		q_allele = 'C';
		break;
	    case 2:
		q_allele = 'G';
		break;
	    case 3:
		q_allele = 'T';
		break;
	    }
	}
	i++;
    }

    return 1;
}

int 
tally_observed_haplotypes(vector<char *> &obshap, int snp_index, char &p_allele, char &q_allele) 
{
    int  nucs[4] = {0};
    char nuc;

    //
    // Pull each allele for this SNP from the observed haplotype.
    //
    for (uint j = 0; j < obshap.size(); j++) {
        nuc = obshap[j][snp_index];

	switch(nuc) {
	case 'A':
	case 'a':
	    nucs[0]++;
	    break;
	case 'C':
	case 'c':
	    nucs[1]++;
	    break;
	case 'G':
	case 'g':
	    nucs[2]++;
	    break;
	case 'T':
	case 't':
	    nucs[3]++;
	    break;
	}
    }

    //
    // Determine how many alleles are present at this position in this population.
    // We cannot deal with more than two alternative alleles, if there are more than two
    // in a single population, print a warning and exclude this nucleotide position.
    //
    int i;
    int allele_cnt = 0;
    for (i = 0; i < 4; i++)
	if (nucs[i] > 0) allele_cnt++;

    if (allele_cnt > 2) {
	p_allele = 0;
	q_allele = 0;
	return -1;
    }

    //
    // Record which nucleotide is the P allele and which is the Q allele.
    //
    p_allele = 0;
    q_allele = 0;

    i = 0;
    while (p_allele == 0 && i < 4) {
	if (nucs[i] > 0) {
	    switch(i) {
	    case 0:
		p_allele = 'A';
		break;
	    case 1:
		p_allele = 'C';
		break;
	    case 2:
		p_allele = 'G';
		break;
	    case 3:
		p_allele = 'T';
		break;
	    }
	}
	i++;
    }
    while (q_allele == 0 && i < 4) {
	if (nucs[i] > 0) {
	    switch(i) {
	    case 1:
		q_allele = 'C';
		break;
	    case 2:
		q_allele = 'G';
		break;
	    case 3:
		q_allele = 'T';
		break;
	    }
	}
	i++;
    }

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
    char *p, *e;

    while (fh.good()) {
	fh.getline(line, id_len);

	if (strlen(line) == 0) continue;

	//
	// Skip commented lines.
	//
	for (p = line; isspace(*p) && *p != '\0'; p++);
	if (*p == '#') continue;

	marker = (int) strtol(line, &e, 10);

	if (*e == '\0')
	    list.insert(marker);
    }

    fh.close();

    if (list.size() == 0) {
 	cerr << "Unable to load any markers from '" << path << "'\n";
	exit(1);
    }

    return 0;
}

int load_marker_column_list(string path, map<int, set<int> > &list) {
    char     line[id_len];
    ifstream fh(path.c_str(), ifstream::in);

    if (fh.fail()) {
        cerr << "Error opening white/black list file '" << path << "'\n";
	exit(1);
    }

    vector<string> parts;
    uint  marker, col;
    char *p, *e;

    uint line_num = 1;
    while (fh.good()) {
	fh.getline(line, id_len);

	if (strlen(line) == 0) continue;

	//
	// Skip commented lines.
	//
	for (p = line; isspace(*p) && *p != '\0'; p++);
	if (*p == '#') continue;

	//
	// Parse the whitelist, we expect:
	// <marker>[<tab><snp column>]
	//
	parse_tsv(line, parts);

	if (parts.size() > 2) {
	    cerr << "Too many columns in whitelist " << path << "' at line " << line_num << "\n";
	    exit(1);

	} else if (parts.size() == 2) {
	    marker = (int) strtol(parts[0].c_str(), &e, 10);
	    if (*e != '\0') {
		cerr << "Unable to parse whitelist, '" << path << "' at line " << line_num << "\n";
		exit(1);
	    }
	    col = (int) strtol(parts[1].c_str(), &e, 10);
	    if (*e != '\0') {
		cerr << "Unable to parse whitelist, '" << path << "' at line " << line_num << "\n";
		exit(1);
	    }
	    list[marker].insert(col);

	} else {
	    marker = (int) strtol(parts[0].c_str(), &e, 10);
	    if (*e != '\0') {
		cerr << "Unable to parse whitelist, '" << path << "' at line " << line << "\n";
		exit(1);
	    }
	    list.insert(make_pair(marker, std::set<int>()));
	}

	line_num++;
    }

    fh.close();

    if (list.size() == 0) {
 	cerr << "Unable to load any markers from '" << path << "'\n";
	help();
    }

    return 0;
}

int 
build_file_list(vector<pair<int, string> > &files, 
		map<int, pair<int, int> > &pop_indexes, 
		map<int, vector<int> > &grp_members) 
{
    char             line[max_len];
    vector<string>   parts;
    map<string, int> pop_key_rev, grp_key_rev;
    set<string>      pop_names, grp_names;
    string f;
    uint   len;

    if (pmap_path.length() > 0) {
	cerr << "Parsing population map.\n";

	ifstream fh(pmap_path.c_str(), ifstream::in);

	if (fh.fail()) {
	    cerr << "Error opening population map '" << pmap_path << "'\n";
	    return 0;
	}

	uint pop_id = 0;
	uint grp_id = 0;

	while (fh.good()) {
	    fh.getline(line, max_len);

	    len = strlen(line);
	    if (len == 0) continue;

	    //
	    // Check that there is no carraige return in the buffer.
	    //
	    if (line[len - 1] == '\r') line[len - 1] = '\0';

	    //
	    // Ignore comments
	    //
	    if (line[0] == '#') continue;

	    //
	    // Parse the population map, we expect:
	    // <file name><tab><population string>[<tab><group string>]
	    //
	    parse_tsv(line, parts);

	    if (parts.size() < 2 || parts.size() > 3) {
		cerr << "Population map is not formated correctly: expecting two or three, tab separated columns, found " << parts.size() << ".\n";
		return 0;
	    }

	    //
	    // Have we seen this population or group before?
	    //
	    if (pop_names.count(parts[1]) == 0) {
		pop_names.insert(parts[1]);
		pop_id++;
		pop_key[pop_id]       = parts[1];
		pop_key_rev[parts[1]] = pop_id;

		//
		// If this is the first time we have seen this population, but not the
		// first time we have seen this group, add the population to the group list.
		//
		if (parts.size() == 3 && grp_key_rev.count(parts[2]) > 0)
		    grp_members[grp_key_rev[parts[2]]].push_back(pop_id);
	    }
	    if (parts.size() == 3 && grp_names.count(parts[2]) == 0) {
		grp_names.insert(parts[2]);
		grp_id++;
		grp_key[grp_id] = parts[2];
		grp_key_rev[parts[2]] = grp_id;

		//
		// Associate the current population with the group.
		//
		grp_members[grp_id].push_back(pop_id);
	    }

	    //
	    // Test that file exists before adding to list.
	    //
	    ifstream test_fh;
	    gzFile   gz_test_fh;

	    f = in_path.c_str() + parts[0] + ".matches.tsv";
	    test_fh.open(f.c_str());

	    if (test_fh.fail()) {
		//
		// Test for a gzipped file.
		//
		f = in_path.c_str() + parts[0] + ".matches.tsv.gz";
		gz_test_fh = gzopen(f.c_str(), "rb");
		if (!gz_test_fh) {
		    cerr << " Unable to find " << f.c_str() << ", excluding it from the analysis.\n";
		} else {
		    gzclose(gz_test_fh);
		    files.push_back(make_pair(pop_key_rev[parts[1]], parts[0]));
		}
	    } else {
		test_fh.close();
		files.push_back(make_pair(pop_key_rev[parts[1]], parts[0]));
	    }
	}

	fh.close();
    } else {
	cerr << "No population map specified, building file list.\n";

	//
	// If no population map is specified, read all the files from the Stacks directory.
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
	    if (pos < file.length()) {
		files.push_back(make_pair(1, file.substr(0, pos)));
	    } else {
		pos = file.rfind(".tags.tsv.gz");
		if (pos < file.length())
		    files.push_back(make_pair(1, file.substr(0, pos)));
	    }
	}

	pop_key[1] = "1";

	closedir(dir);
    }

    if (files.size() == 0) {
	cerr << "Unable to locate any input files to process within '" << in_path << "'\n";
	return 0;
    }

    //
    // Sort the files according to population ID.
    //
    sort(files.begin(), files.end(), compare_pop_map);

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

    pop_indexes.size() == 1 ?
	cerr << "  " << pop_indexes.size() << " population found\n" :
	cerr << "  " << pop_indexes.size() << " populations found\n";

    if (population_limit > (int) pop_indexes.size()) {
	cerr << "Population limit (" 
	     << population_limit 
	     << ") larger than number of popualtions present, adjusting parameter to " 
	     << pop_indexes.size() << "\n";
	population_limit = pop_indexes.size();
    }

    map<int, pair<int, int> >::iterator it;
    for (it = pop_indexes.begin(); it != pop_indexes.end(); it++) {
	start = it->second.first;
	end   = it->second.second;
	cerr << "    " << pop_key[it->first] << ": ";
	for (int i = start; i <= end; i++) {
	    cerr << files[i].second;
	    if (i < end) cerr << ", ";
	}
	cerr << "\n";
    }

    //
    // If no group membership is specified in the population map, create a default 
    // group with each population ID as a member.
    //
    if (grp_members.size() == 0) {
	for (it = pop_indexes.begin(); it != pop_indexes.end(); it++)
	    grp_members[1].push_back(it->first);
	grp_key[1] = "1";
    }

    grp_members.size() == 1 ?
	cerr << "  " << grp_members.size() << " group of populations found\n" :
	cerr << "  " << grp_members.size() << " groups of populations found\n";
    map<int, vector<int> >::iterator git;
    for (git = grp_members.begin(); git != grp_members.end(); git++) {
	cerr << "    " << grp_key[git->first] << ": ";
	for (uint i = 0; i < git->second.size(); i++) {
	    cerr << pop_key[git->second[i]];
	    if (i < git->second.size() - 1) cerr << ", ";
	}
	cerr << "\n";
    }

    return 1;
}

bool compare_pop_map(pair<int, string> a, pair<int, string> b) {
    if (a.first == b.first)
	return (a.second < b.second);
    return (a.first < b.first);
}

bool hap_compare(pair<string, int> a, pair<string, int> b) {
    return (a.second > b.second);
}

bool compare_genpos(GenPos a, GenPos b) {
    return (a.bp < b.bp);
}

int parse_command_line(int argc, char* argv[]) {
    int c;
     
    while (1) {
	static struct option long_options[] = {
	    {"help",           no_argument,       NULL, 'h'},
            {"version",        no_argument,       NULL, 'v'},
            {"verbose",        no_argument,       NULL, 'd'},
            {"sql",            no_argument,       NULL, 's'},
            {"vcf",            no_argument,       NULL, 'V'},
            {"vcf_haplotypes", no_argument,       NULL, 'n'},
            {"fasta",          no_argument,       NULL, 'F'},
            {"fasta_strict",   no_argument,       NULL, 'J'},
            {"structure",      no_argument,       NULL, 'S'},
            {"fastphase",      no_argument,       NULL, 'A'},
            {"phase",          no_argument,       NULL, 'C'},
            {"beagle",         no_argument,       NULL, 'E'},
            {"beagle_phased",  no_argument,       NULL, 'H'},
            {"plink",          no_argument,       NULL, 'K'},
            {"genomic",        no_argument,       NULL, 'g'},
	    {"genepop",        no_argument,       NULL, 'G'},
	    {"phylip",         no_argument,       NULL, 'Y'},
	    {"phylip_var",     no_argument,       NULL, 'L'},
	    {"phylip_var_all", no_argument,       NULL, 'T'},
	    {"hzar",           no_argument,       NULL, 'Z'},
	    {"treemix",        no_argument,       NULL, 'U'},
	    {"merge_sites",    no_argument,       NULL, 'D'},
	    {"window_size",    required_argument, NULL, 'w'},
	    {"num_threads",    required_argument, NULL, 't'},
	    {"batch_id",       required_argument, NULL, 'b'},
	    {"in_path",        required_argument, NULL, 'P'},
	    {"progeny",        required_argument, NULL, 'r'},
	    {"min_depth",      required_argument, NULL, 'm'},
	    {"renz",           required_argument, NULL, 'e'},
	    {"pop_map",        required_argument, NULL, 'M'},
	    {"whitelist",      required_argument, NULL, 'W'},
	    {"blacklist",      required_argument, NULL, 'B'},
	    {"write_single_snp",  no_argument,       NULL, 'I'},
	    {"write_random_snp",  no_argument,       NULL, 'j'},
	    {"ordered_export",    no_argument,       NULL, 'N'},
            {"kernel_smoothed",   no_argument,       NULL, 'k'},
            {"fstats",            no_argument,       NULL, '6'},
            {"log_fst_comp",      no_argument,       NULL, 'l'},
            {"bootstrap_type",    required_argument, NULL, 'O'},
	    {"bootstrap_reps",    required_argument, NULL, 'R'},
	    {"bootstrap_wl",      required_argument, NULL, 'Q'},
            {"bootstrap",         no_argument,       NULL, '1'},
            {"bootstrap_fst",     no_argument,       NULL, '2'},
            {"bootstrap_phist",   no_argument,       NULL, '3'},
            {"bootstrap_div",     no_argument,       NULL, '4'},
            {"bootstrap_pifis",   no_argument,       NULL, '5'},
	    {"min_populations",   required_argument, NULL, 'p'},
	    {"min_maf",           required_argument, NULL, 'a'},
	    {"max_obs_het",       required_argument, NULL, 'q'},
	    {"lnl_lim",           required_argument, NULL, 'c'},
	    {"merge_prune_lim",   required_argument, NULL, 'i'},
	    {"fst_correction",    required_argument, NULL, 'f'},
	    {"p_value_cutoff",    required_argument, NULL, 'u'},
	    {0, 0, 0, 0}
	};
	
	// getopt_long stores the option index here.
	int option_index = 0;
     
	c = getopt_long(argc, argv, "ACDEFGHJKLNSTUVYZ123456dghjklnsva:b:c:e:f:i:m:o:p:q:r:t:u:w:B:I:M:O:P:R:Q:W:", long_options, &option_index);

	// Detect the end of the options.
	if (c == -1)
	    break;
     
	switch (c) {
	case 'h':
	    help();
	    break;
	case 'd':
	    verbose = true;
	    break;
	case 't':
	    num_threads = atoi(optarg);
	    break;
	case 'P':
	    in_path = optarg;
	    break;
	case 'M':
	    pmap_path = optarg;
	    break;
	case 'D':
	    merge_sites = true;
	    break;
	case 'i':
	    merge_prune_lim = is_double(optarg);
	    break;
	case 'q':
	    max_obs_het = is_double(optarg);
	    break;
	case 'b':
	    batch_id = is_integer(optarg);
	    if (batch_id < 0) {
		cerr << "Batch ID (-b) must be an integer, e.g. 1, 2, 3\n";
		help();
	    }
	    break;
	case 'r':
	    sample_limit = atof(optarg);
	    break;
	case 'p':
	    population_limit = atoi(optarg);
	    break;
	case 'k':
	    kernel_smoothed = true;
	    calc_fstats     = true;
	    break;
	case '6':
	    calc_fstats = true;
	    break;
 	case 'l':
	    log_fst_comp = true;
	    break;
	case '1':
	    bootstrap       = true;
	    bootstrap_fst   = true;
	    bootstrap_phist = true;
	    bootstrap_pifis = true;
	    bootstrap_div   = true;
	    break;
	case '2':
	    bootstrap_fst   = true;
	    break;
	case '3':
	    bootstrap_phist = true;
	    break;
	case '4':
	    bootstrap_div = true;
	    break;
	case '5':
	    bootstrap_pifis = true;
	    break;
	case 'O':
	    if (strcasecmp(optarg, "exact") == 0)
		bootstrap_type = bs_exact;
	    else if (strcasecmp(optarg, "approx") == 0)
		bootstrap_type = bs_approx;
	    else {
		cerr << "Unknown bootstrap type specified '" << optarg << "'\n";
		help();
	    }
	    break;
	case 'R':
	    bootstrap_reps = atoi(optarg);
	    break;
	case 'Q':
	    bs_wl_file = optarg;
	    bootstrap_wl = true;
	    break;
	case 'c':
	    lnl_limit  = is_double(optarg);
	    filter_lnl = true;
	    break;
	case 'I':
	    write_single_snp = true;
	    break;
	case 'j':
	    write_random_snp = true;
	    break;
	case 'N':
	    ordered_export = true;
	    break;
	case 's':
	    sql_out = true;
	    break;
	case 'V':
	    vcf_out = true;
	    break;
	case 'n':
	    vcf_haplo_out = true;
	    break;
	case 'F':
	    fasta_out = true;
	    break;
	case 'J':
	    fasta_strict_out = true;
	    break;
	case 'G':
	    genepop_out = true;
	    break;
	case 'S':
	    structure_out = true;
	    break;
	case 'A':
	    fastphase_out = true;
	    break;
	case 'C':
	    phase_out = true;
	    break;
	case 'E':
	    beagle_out = true;
	    break;
	case 'H':
	    beagle_phased_out = true;
	    break;
	case 'K':
	    plink_out = true;
	    break;
	case 'Z':
	    hzar_out = true;
	    break;
	case 'Y':
	    phylip_out = true;
	    break;
	case 'L':
	    phylip_var = true;
	    break;
	case 'T':
	    phylip_var_all = true;
	    break;
	case 'U':
	    treemix_out = true;
	    break;
	case 'g':
	    genomic_out = true;
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
	case 'a':
	    minor_allele_freq = atof(optarg);
	    break;
	case 'f':
	    if (strcasecmp(optarg, "p_value") == 0)
		fst_correction = p_value;
	    else if (strcasecmp(optarg, "bonferroni_win") == 0)
		fst_correction = bonferroni_win;
	    else if (strcasecmp(optarg, "bonferroni_gen") == 0)
		fst_correction = bonferroni_gen;
	    else {
		cerr << "Unknown Fst correction specified '" << optarg << "'\n";
		help();
	    }
	    break;
	case 'u':
	    p_value_cutoff = atof(optarg);
	    break;
	case 'e':
	    enz = optarg;
	    break;
	case 'w':
	    sigma = atof(optarg);
	    break;
        case 'v':
            version();
            break;
	case '?':
	    // getopt_long already printed an error message.
	    help();
	    break;
	default:
	    cerr << "Unknown command line option: '" << (char) c << "'\n";
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
	cerr << "A population map was not specified, all samples will be read from '" << in_path << "' as a single popultaion.\n";
    }

    if (batch_id < 0) {
	cerr << "You must specify a batch ID.\n";
	help();
    }

    if (enz.length() > 0 && renz.count(enz) == 0) {
	cerr << "Unrecognized restriction enzyme specified: '" << enz.c_str() << "'.\n";
	help();
    }

    if (merge_prune_lim != 1.0) {
	if (merge_prune_lim > 1.0)
	    merge_prune_lim = merge_prune_lim / 100;

	if (merge_prune_lim < 0 || merge_prune_lim > 1.0) {
	    cerr << "Unable to parse the merge sites pruning limit.\n";
	    help();
	}
    }

    if (minor_allele_freq > 0) {
	if (minor_allele_freq > 1)
	    minor_allele_freq = minor_allele_freq / 100;

	if (minor_allele_freq > 0.5) {
	    cerr << "Unable to parse the minor allele frequency.\n";
	    help();
	}
    }

    if (max_obs_het != 1.0) {
	if (max_obs_het > 1)
	    max_obs_het = max_obs_het / 100;

	if (max_obs_het < 0 || max_obs_het > 1.0) {
	    cerr << "Unable to parse the maximum observed heterozygosity.\n";
	    help();
	}
    }

    if (sample_limit > 0) {
	if (sample_limit > 1)
	    sample_limit = sample_limit / 100;

	if (sample_limit > 1.0) {
	    cerr << "Unable to parse the sample limit frequency\n";
	    help();
	}
    }

    if (write_single_snp && write_random_snp) {
	cerr << "Please specify either '--write_single_snp' or '--write_random_snp', not both.\n";
	help();
    }

    if (merge_sites == true && enz.length() == 0) {
	cerr << "You must specify the restriction enzyme associated with this data set to merge overlaping cutsites.\n";
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
              << "populations -b batch_id -P path -M path [-r min] [-m min] [-B blacklist] [-W whitelist] [-s] [-e renz] [-t threads] [-v] [-h]" << "\n"
	      << "  b: Batch ID to examine when exporting from the catalog.\n"
	      << "  P: path to the Stacks output files.\n"
	      << "  M: path to the population map, a tab separated file describing which individuals belong in which population.\n"
	      << "  s: output a file to import results into an SQL database.\n"
	      << "  B: specify a file containing Blacklisted markers to be excluded from the export.\n"
	      << "  W: specify a file containing Whitelisted markers to include in the export.\n"
	      << "  e: restriction enzyme, required if generating 'genomic' output.\n"
	      << "  t: number of threads to run in parallel sections of code.\n"
	      << "  v: print program version." << "\n"
	      << "  h: display this help messsage." << "\n\n"
	      << "  Merging and Phasing:\n"
	      << "    --merge_sites: merge loci that were produced from the same restriction enzyme cutsite (requires reference-aligned data).\n"
	      << "    --merge_prune_lim: when merging adjacent loci, if at least X% samples posses both loci prune the remaining samples out of the analysis.\n"
	      << "  Data Filtering:\n"
	      << "    r: minimum percentage of individuals in a population required to process a locus for that population.\n"
	      << "    p: minimum number of populations a locus must be present in to process a locus.\n"
	      << "    m: specify a minimum stack depth required for individuals at a locus.\n"
	      << "    f: specify a correction to be applied to Fst values: 'p_value', 'bonferroni_win', or 'bonferroni_gen'.\n"
	      << "    --min_maf: specify a minimum minor allele frequency required to process a nucleotide site at a locus (0 < min_maf < 0.5).\n"
	      << "    --max_obs_het: specify a maximum observed heterozygosity required to process a nucleotide site at a locus.\n"
	      << "    --p_value_cutoff [num]: required p-value to keep an Fst measurement (0.05 by default). Also used as base for Bonferroni correction.\n"
	      << "    --lnl_lim [num]: filter loci with log likelihood values below this threshold.\n"
	      << "    --write_single_snp: restrict data analysis to only the first SNP per locus.\n"
	      << "    --write_random_snp: restrict data analysis to one random SNP per locus.\n\n"
	      << "  Fstats:\n"
	      << "    --fstats: enable SNP and haplotype-based F statistics.\n\n"
	      << "  Kernel-smoothing algorithm:\n" 
	      << "    k: enable kernel-smoothed Pi, Fis, Fst, Fst', and Phi_st calculations.\n"
	      << "    --window_size [num]: distance over which to average values (sigma, default 150,000bp; window is 3sigma in length).\n\n"
	      << "  Bootstrap Resampling:\n" 
	      << "    --bootstrap: turn on boostrap resampling for all smoothed statistics.\n"
	      << "    --bootstrap_pifis: turn on boostrap resampling for smoothed SNP-based Pi and Fis calculations.\n"
	      << "    --bootstrap_fst: turn on boostrap resampling for smoothed Fst calculations based on pairwise population comparison of SNPs.\n"
	      << "    --bootstrap_div: turn on boostrap resampling for smoothed haplotype diveristy and gene diversity calculations based on haplotypes.\n"
	      << "    --bootstrap_phist: turn on boostrap resampling for smoothed Phi_st calculations based on haplotypes.\n"
	      << "    --bootstrap_reps [num]: number of bootstrap resamplings to calculate (default 100).\n"
	      << "    --bootstrap_wl [path]: only bootstrap loci contained in this whitelist.\n\n"
	      << "  File ouput options:\n"
	      << "    --ordered_export: if data is reference aligned, exports will be ordered; only a single representative of each overlapping site.\n"
	      << "    --genomic: output each nucleotide position (fixed or polymorphic) in all population members to a file.\n"
	      << "    --fasta: output full sequence for each unique haplotype, from each sample locus in FASTA format, regardless of plausibility.\n"
	      << "    --fasta_strict: output full sequence for each haplotype, from each sample locus in FASTA format, only for biologically plausible loci.\n"
	      << "    --vcf: output SNPs in Variant Call Format (VCF).\n"
	      << "    --vcf_haplotypes: output haplotypes in Variant Call Format (VCF).\n"
	      << "    --genepop: output results in GenePop format.\n"
	      << "    --structure: output results in Structure format.\n"
	      << "    --phase: output genotypes in PHASE format.\n"
	      << "    --fastphase: output genotypes in fastPHASE format.\n"
	      << "    --beagle: output genotypes in Beagle format.\n"
	      << "    --beagle_phased: output haplotypes in Beagle format.\n"
	      << "    --plink: output genotypes in PLINK format.\n"
	      << "    --hzar: output genotypes in Hybrid Zone Analysis using R (HZAR) format.\n"
	      << "    --phylip: output nucleotides that are fixed-within, and variant among populations in Phylip format for phylogenetic tree construction.\n"
	      << "    --phylip_var: include variable sites in the phylip output encoded using IUPAC notation.\n"
	      << "    --phylip_var_all: include all sequence as well as variable sites in the phylip output encoded using IUPAC notation.\n"
	      << "    --treemix: output SNPs in a format useable for the TreeMix program (Pickrell and Pritchard).\n\n"
	      << "  Debugging:\n"
	      << "    --verbose: turn on additional logging.\n"
	      << "    --log_fst_comp: log components of Fst/Phi_st calculations to a file.\n";
    // << "    --bootstrap_type [exact|approx]: enable bootstrap resampling for population statistics (reference genome required).\n"

    exit(0);
}
