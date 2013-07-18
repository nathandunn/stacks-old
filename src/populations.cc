// -*-mode:c++; c-style:k&r; c-basic-offset:4;-*-
//
// Copyright 2012-2013, Julian Catchen <jcatchen@uoregon.edu>
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
double    sigma             = 150000;
double    sample_limit      = 0;
int       progeny_limit     = 0;
int       population_limit  = 1;
bool      bootstrap         = false;
bs_type   bootstrap_type    = bs_none;
int       bootstrap_reps    = 100;
bool      corrections       = false;
bool      write_single_snp  = false;
bool      expand_id         = false;
bool      sql_out           = false;
bool      vcf_out           = false;
bool      genepop_out       = false;
bool      genomic_out       = false;
bool      structure_out     = false;
bool      phase_out         = false;
bool      beagle_out        = false;
bool      plink_out         = false;
bool      phylip_out        = false;
bool      phylip_var        = false;
bool      kernel_smoothed   = false;
bool      linkage_stats     = false;
bool      loci_ordered      = false;
bool      log_fst_comp      = false;
int       min_stack_depth   = 0;
double    minor_allele_freq = 0;
double    p_value_cutoff    = 0.05;
corr_type fst_correction    = no_correction;

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
	<< "Minor allele frequency cutoff: " << minor_allele_freq << "\n"
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
    if (!build_file_list(files, pop_indexes))
	exit(1);

    if (wl_file.length() > 0) {
	load_marker_list(wl_file, whitelist);
	cerr << "Loaded " << whitelist.size() << " whitelisted markers.\n";
    }
    if (bl_file.length() > 0) {
	load_marker_list(bl_file, blacklist);
	cerr << "Loaded " << blacklist.size() << " blacklisted markers.\n";
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

    apply_locus_constraints(catalog, pmap, pop_indexes);

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
    	cerr << "Generating nucleotide-level summary statistics for population " << pop_id << "\n";
    	psum->add_population(catalog, pmap, pop_id, start_index, end_index, log_fh);

	if (kernel_smoothed && loci_ordered) {
	    cerr << "  Generating kernel-smoothed population statistics";
	    if (bootstrap) cerr << " and bootstrap resampling";
	    cerr << "...\n";
	    kernel_smoothed_popstats(catalog, pmap, psum, pop_id, log_fh);
	}
    }

    cerr << "Tallying loci across populations...";
    psum->tally(catalog);
    cerr << "done.\n";

    //
    // Idenitfy polymorphic loci, tabulate haplotypes present.
    //
    tabulate_haplotypes(catalog, pmap);

    //
    // Output a list of heterozygous loci and the associate haplotype frequencies.
    //
    if (sql_out)
	write_sql(catalog, pmap);

    //
    // Output the locus-level summary statistics.
    //
    write_summary_stats(files, pop_indexes, catalog, pmap, psum);

    //
    // Output data in requested formats
    //
    if (vcf_out)
	write_vcf(catalog, pmap, psum, samples, sample_ids);

    if (genepop_out)
	write_genepop(catalog, pmap, psum, pop_indexes, samples);

    if (structure_out)
	write_structure(catalog, pmap, psum, pop_indexes, samples);

    if (phase_out)
	write_phase(catalog, pmap, psum, pop_indexes, samples);

    if (beagle_out)
	write_beagle(catalog, pmap, psum, pop_indexes, samples);

    if (plink_out)
	write_plink(catalog, pmap, psum, pop_indexes, samples);

    if (phylip_out)
	write_phylip(catalog, pmap, psum, pop_indexes, samples);

    //
    // Output the observed haplotypes.
    //
    write_generic(catalog, pmap, samples, false);

    //
    // Calculate and write Fst.
    //
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
			map<int, pair<int, int> > &pop_indexes)
{
    uint pop_id, start_index, end_index;
    CSLocus *loc;
    Datum  **d;

    if (sample_limit == 0 && population_limit == 0 && min_stack_depth == 0) return 0;

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
    bool   pro_limit = false;
    bool   pop_limit = false;
    int    pops      = 0;
    int    below_stack_dep = 0;
    set<int> blacklist;

    for (it = catalog.begin(); it != catalog.end(); it++) {
	loc = it->second;
	d   = pmap->locus(loc->id);

	//
	// Check that each sample is over the minimum stack depth for this locus.
	//
	if (min_stack_depth > 0) 
	    for (int i = 0; i < pmap->sample_cnt(); i++) {
		if (d[i] != NULL && d[i]->tot_depth < min_stack_depth) {
		    below_stack_dep++;
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
	// Check that the counts for each population are over progeny_limit. If not, zero out 
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

	if (pop_limit)
	    blacklist.insert(loc->id);

	for (uint i = 0; i < pop_cnt; i++)
	    pop_cnts[i] = 0;
	pro_limit = false;
	pop_limit = false;
	pops      = 0;
    }

    //
    // Remove loci
    //
    if (min_stack_depth > 0) 
	cerr << "Removed " << below_stack_dep << " samples from loci that are below the minimum stack depth of " << min_stack_depth << "x\n";
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
reduce_catalog(map<int, CSLocus *> &catalog, set<int> &whitelist, set<int> &blacklist) 
{
    map<int, CSLocus *> list;
    map<int, CSLocus *>::iterator it;
    CSLocus *loc;

    if (whitelist.size() == 0 && blacklist.size() == 0) 
	return 0;
 
    int i = 0;
    for (it = catalog.begin(); it != catalog.end(); it++) {
	loc = it->second;

	if (whitelist.size() > 0 && whitelist.count(loc->id) == 0) continue;
	if (blacklist.count(loc->id)) continue;

	list[it->first] = it->second;
	i++;
    }

    catalog = list;

    return i;
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

int tabulate_haplotypes(map<int, CSLocus *> &catalog, PopMap<CSLocus> *pmap) {
    map<int, CSLocus *>::iterator it;
    vector<char *>::iterator hit;
    Datum  **d;
    CSLocus *loc;

    for (it = catalog.begin(); it != catalog.end(); it++) {
	loc = it->second;
	d   = pmap->locus(loc->id);

	for (int i = 0; i < pmap->sample_cnt(); i++) {
	    if (d[i] == NULL) 
		continue;

	    if (d[i]->obshap.size() > 1)
		loc->marker = "heterozygous";
	}

	if (loc->marker.length() > 0) {
     	    create_genotype_map(loc, pmap);
	    call_population_genotypes(loc, pmap);
	}
    }

    return 0;
}

int create_genotype_map(CSLocus *locus, PopMap<CSLocus> *pmap) {
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
write_linkage_stats(map<int, pair<int, int> > &pop_indexes, map<int, CSLocus *> &catalog, 
		    PopMap<CSLocus> *pmap, PopSum<CSLocus> *psum) 
{
    //
    // We want to iterate over each population and calculate D', r^2 and chi-square significance
    // for the first nucleotide of each locus.
    //
    // map<string, vector<CSLocus *> >::iterator it;

    // vector<int> pops;
    // map<int, pair<int, int> >::iterator pit;
    // for (pit = pop_indexes.begin(); pit != pop_indexes.end(); pit++)
    // 	pops.push_back(pit->first);

    // int      num_loci = pmap->loci_cnt();
    // double **dprime   = new double *[num_loci];
    // double **rsq      = new double *[num_loci];
    // double **chisq    = new double *[num_loci];

    // for (uint i = 0; i < num_loci; i++) {
    // 	dprime[i] = new double[num_loci];
    // 	rsq[i]    = new double[num_loci];
    // 	chisq[i]  = new double[num_loci];
    // }

    // for (uint i = 0; i < pops.size(); i++) {

    // 	int start = pop_indexes[pops[i]].first;
    // 	int end   = pop_indexes[pops[i]].second;

    // 	//
    // 	// Clear the two-dimensional arrays for the next population.
    // 	//
    // 	for (uint j = 0; j < num_loci; j++) {
    // 	    memset(dprime[j], 0, num_loci);
    // 	    memset(rsq[j],    0, num_loci);
    // 	    memset(chisq[j],  0, num_loci);
    // 	}

    // 	CSLocus  *loc_1;
    // 	LocSum  **s_1;
    // 	Datum   **d_1;

    // 	uint col_1, col_2, tot;
    // 	char p_nuc_1, q_nuc_1, p_nuc_2, q_nuc_2;

    // 	for (it = pmap->ordered_loci.begin(); it != pmap->ordered_loci.end(); it++) {
    // 	    string chr = it->first;
    // 	    for (uint pos_1 = 0; pos_1 < it->second.size(); pos_1++) {
    // 		loc_1 = it->second[pos_1];
    // 		s_1   = psum->locus(loc_1->id);
    // 		d_1   = pmap->locus(loc_1->id);

    // 		for (uint k = 0; k < loc_1->snps.size(); k++) {
    // 		    col_1   = loc_1->snps[k]->col;
    // 		    p_nuc_1 = s_1->nucs[col_1].p_nuc;
    // 		    q_nuc_1 = s_1->nucs[col_1].q_nuc;

    // 		    //
    // 		    // Compare this locus against every other locus in the data set.
    // 		    //
    // 		    CSLocus  *loc_2;
    // 		    LocSum  **s_2;
    // 		    Datum   **d_2;

    // 		    map<string, vector<CSLocus *> >::iterator chr_it;
    // 		    for (uint pos_2 = 0; pos_2 < chr_it->second.size(); pos_2++) {
    // 			loc_2 = chr_it->second[pos_2];
    // 			s_2   = psum->locus(loc_2->id);
    // 			d_2   = pmap->locus(loc_2->id);

    // 			for (uint n = 0; n < loc_2->snps.size(); n++) {
    // 			    col_2   = loc_2->snps[n]->col;
    // 			    p_nuc_2 = s_2->nucs[col_2].p_nuc;
    // 			    q_nuc_2 = s_2->nucs[col_2].q_nuc;

    // 			    tot = 0;

    // 			    for (uint m = start; m <= end; m++) {
    // 				if (d_1[m] == NULL || col_1 >= d_1[m]->len ||
    // 				    d_2[m] == NULL || col_2 >= d_2[m]->len) 
    // 				    continue;
    // 				if (d_1[m]->obshap.size() < 2 ||
    // 				    d_2[m]->obshap.size() < 2 ||
    // 				    psum->tally_observed_haplotypes(d_1[m]->obshap, k) != 2 ||
    // 				    psum->tally_observed_haplotypes(d_2[m]->obshap, n) != 2)
    // 				    continue;

    // 				tot++;

    // 				if (d_1[m]->obshap[0] == q_nuc_1) {
    // 				    if (d_2[m]->obshap[0] == q_nuc_2)
    // 					q1q2++;
    // 				    else if (d_2[m]->obshap[0] == p_nuc_2)
    // 					q1p2++;
    // 				} else if (d_1[m]->obshap[0] == p_nuc_1) {
    // 				    if (d_2[m]->obshap[0] == q_nuc_2)
    // 					p1q2++;
    // 				    else if (d_2[m]->obshap[0] == p_nuc_2)
    // 					p1p2++;
    // 				} else if (d_1[m]->obshap[1] == q_nuc_1) {
    // 				    if (d_2[m]->obshap[0] == q_nuc_2)
    // 					q1q2++;
    // 				    else if (d_2[m]->obshap[0] == p_nuc_2)
    // 					q1p2++;
    // 				} else if (d_1[m]->obshap[1] == p_nuc_1) {
    // 				    if (d_2[m]->obshap[0] == q_nuc_2)
    // 					p1q2++;
    // 				    else if (d_2[m]->obshap[0] == p_nuc_2)
    // 					p1p2++;
    // 				}
    // 			    }

    // 			    if (write_single_snp) 
    // 				break;
    // 			}
    // 		    }
		    
    // 		    if (write_single_snp) 
    // 			break;
    // 		}
    // 	    }
    // 	}
    // }

    // //
    // // Free the two-dimensional arrays.
    // //
    // for (uint i = 0; i < num_loci; i++) {
    // 	delete [] dprime[i];
    //     delete [] rsq[i];
    // 	delete [] chisq[i];
    // }
    // delete [] dprime;
    // delete [] rsq;
    // delete [] chisq;

    return 0;
}

int 
write_summary_stats(vector<pair<int, string> > &files, map<int, pair<int, int> > &pop_indexes, 
		    map<int, CSLocus *> &catalog, PopMap<CSLocus> *pmap, PopSum<CSLocus> *psum) 
{
    map<string, vector<CSLocus *> >::iterator it;
    CSLocus  *loc;
    LocSum  **s;
    LocTally *t;
    int       len;
    int       pop_cnt = psum->pop_cnt();
    char      fisstr[32], wfisstr[32], wpistr[32], sitesstr[32];

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
			pi_mean[j]       += s[j]->nucs[i].pi;
			fis_mean[j]      += s[j]->nucs[i].Fis;

			n_all[j]++;
			num_indv_mean_all[j] += s[j]->nucs[i].num_indv;
			p_mean_all[j]        += s[j]->nucs[i].p;
			obs_het_mean_all[j]  += s[j]->nucs[i].obs_het;
			obs_hom_mean_all[j]  += s[j]->nucs[i].obs_hom;
			exp_het_mean_all[j]  += s[j]->nucs[i].exp_het;
			exp_hom_mean_all[j]  += s[j]->nucs[i].exp_hom;
			pi_mean_all[j]       += s[j]->nucs[i].pi;
			fis_mean_all[j]      += s[j]->nucs[i].Fis;
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
			pi_mean_all[j]       += s[j]->nucs[i].pi;
			fis_mean_all[j]      += s[j]->nucs[i].Fis;
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
		// If this site is fixed in all populations, DON'T output it. If it is variable
		// or fixed within populations but variable among, DO output it.
		//
		if (t->nucs[i].allele_cnt == 2) {

		    for (int j = 0; j < pop_cnt; j++) {

			if (s[j]->nucs[i].num_indv == 0) continue;

			sprintf(fisstr,  "%0.10f", s[j]->nucs[i].Fis);
			sprintf(wfisstr, "%0.10f", s[j]->nucs[i].wFis);
			sprintf(wpistr,  "%0.10f", s[j]->nucs[i].wPi);

			fh << batch_id << "\t"
			   << loc->id << "\t"
			   << loc->loc.chr << "\t"
			   << loc->sort_bp(i) << "\t"
			   << i << "\t"
			   << psum->rev_pop_index(j) << "\t"
			   << s[j]->nucs[i].p_nuc << "\t";

			if (s[j]->nucs[i].q_nuc != 0) 
			    fh << s[j]->nucs[i].q_nuc;

			fh << "\t" << s[j]->nucs[i].num_indv << "\t"
			   << s[j]->nucs[i].p         << "\t"
			   << s[j]->nucs[i].obs_het   << "\t"
			   << s[j]->nucs[i].obs_hom   << "\t"
			   << s[j]->nucs[i].exp_het   << "\t"
			   << s[j]->nucs[i].exp_hom   << "\t"
			   << s[j]->nucs[i].pi        << "\t"
			   << wpistr                  << "\t"
			   << s[j]->nucs[i].wPi_pval  << "\t"
			   << fisstr                  << "\t"
			   << wfisstr                 << "\t"
			   << s[j]->nucs[i].wFis_pval << "\t";
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
			pi_var[j]       += pow((s[j]->nucs[i].pi       - pi_mean[j]),       2);
			fis_var[j]      += pow((s[j]->nucs[i].Fis      - fis_mean[j]),      2);

			num_indv_var_all[j] += pow((s[j]->nucs[i].num_indv - num_indv_mean_all[j]), 2);
			p_var_all[j]        += pow((s[j]->nucs[i].p        - p_mean_all[j]),        2);
			obs_het_var_all[j]  += pow((s[j]->nucs[i].obs_het  - obs_het_mean_all[j]),  2);
			obs_hom_var_all[j]  += pow((s[j]->nucs[i].obs_hom  - obs_hom_mean_all[j]),  2);
			exp_het_var_all[j]  += pow((s[j]->nucs[i].exp_het  - exp_het_mean_all[j]),  2);
			exp_hom_var_all[j]  += pow((s[j]->nucs[i].exp_hom  - exp_hom_mean_all[j]),  2);
			pi_var_all[j]       += pow((s[j]->nucs[i].pi       - pi_mean_all[j]),       2);
			fis_var_all[j]      += pow((s[j]->nucs[i].Fis      - fis_mean_all[j]),      2);
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
			pi_var_all[j]       += pow((s[j]->nucs[i].pi       - pi_mean_all[j]),       2);
			fis_var_all[j]      += pow((s[j]->nucs[i].Fis      - fis_mean_all[j]),      2);
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
	fh << psum->rev_pop_index(j) << "\t" 
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
	sprintf(sitesstr, "%.0f", n_all[j]);

	fh << psum->rev_pop_index(j)     << "\t" 
	   << private_cnt[j]             << "\t"
	   << sitesstr                   << "\t"
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

    double *weights = calculate_weights();

    for (uint i = 0; i < pops.size(); i++) {

	for (uint j = i + 1; j < pops.size(); j++) {
	    int pop_1 = pops[i];
	    int pop_2 = pops[j];

	    double sum = 0.0;
	    double cnt = 0.0;

	    int incompatible_loci = 0;
	    int multiple_loci     = 0;

	    stringstream pop_name;
	    pop_name << "batch_" << batch_id << ".fst_" << pop_1 << "-" << pop_2 << ".tsv";

	    string file = in_path + pop_name.str();
	    ofstream fh(file.c_str(), ofstream::out);

	    if (fh.fail()) {
		cerr << "Error opening Fst output file '" << file << "'\n";
		exit(1);
	    }

	    //
	    // If requested, log Fst component calculations to a file.
	    //
	    ofstream *comp_fh = NULL;
	    if (log_fst_comp) {
		pop_name.str("");
		pop_name << "batch_" << batch_id << ".fst-comp_" << pop_1 << "-" << pop_2 << ".tsv";

		file = in_path + pop_name.str();
		comp_fh = new ofstream(file.c_str(), ofstream::out);

		if (comp_fh->fail()) {
		    cerr << "Error opening Fst component output file '" << file << "'\n";
		    exit(1);
		}

		*comp_fh << "# " << "n_1" << "\t"
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
			 << "binomial_fst" << "\t\t"
			 << "p_1_freq" << "\t"
			 << "q_1_freq" << "\t"
			 << "p_2_freq" << "\t"
			 << "q_2_freq" << "\t"
			 << "p_avg_cor" << "\t"
			 << "n_avg_cor" << "\t"
			 << "amova_fst" << "\n";
	    }

	    cerr << "Calculating Fst for populations " << pop_1 << " and " << pop_2 << " and writing it to file, '" << file << "'\n";

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
	       << "Smoothed Fst P-value" << "\t"
	       << "AMOVA Fst" << "\t"
	       << "Corrected AMOVA Fst" << "\t"
	       << "Smoothed AMOVA Fst" << "\t"
	       << "Window SNP Count" << "\n";

	    map<string, vector<CSLocus *> >::iterator it;
	    CSLocus *loc;
	    PopPair *pair;
	    int      len;
	    char     fst_str[32], wfst_str[32], cfst_str[32], afst_str[32], cafst_str[32], wafst_str[32];

	    map<string, vector<PopPair *> > genome_pairs;
	    vector<double> fst_samples;
	    vector<int>    allele_depth_samples;
	    int            snp_dist[max_snp_dist] = {0};

	    for (it = pmap->ordered_loci.begin(); it != pmap->ordered_loci.end(); it++) {
		string chr = it->first;
		map<uint, uint> pairs_key;

		init_chr_pairs(genome_pairs, chr, pairs_key, it->second);

		vector<PopPair *> &pairs = genome_pairs[chr];

		for (uint pos = 0; pos < it->second.size(); pos++) {
		    loc = it->second[pos];
		    len = strlen(loc->con);

		    for (int k = 0; k < len; k++) {

			pair = psum->Fst(loc->id, pop_1, pop_2, k, comp_fh);

			//
			// Locus is incompatible, log this position.
			//
			if (pair == NULL) {
			    incompatible_loci++;
			    log_fh << "between_population\t"
				   << "incompatible_locus\t"
				   << loc->id << "\t"
				   << loc->loc.chr << "\t"
				   << loc->sort_bp(k) << "\t"
				   << k << "\t" 
				   << pop_1 << "\t" 
				   << pop_2 << "\n";
			    delete pair;
			    continue;
			}

			pair->loc_id = loc->id;
			pair->bp     = loc->sort_bp(k);
			pair->col    = k;

			//
			// Locus is fixed in both populations, or was only found in one population.
			//
			if (pair->pi == 0) {
			    delete pair;
			    continue;
			}

			//
			// Check if this basepair position is already covered by a RAD site.
			//
			if (pairs[pairs_key[pair->bp]] != NULL) {
			    multiple_loci++;
			    log_fh << "between_population\t"
				   << "multiple_locus\t"
				   << loc->id << "\t"
				   << loc->loc.chr << "\t"
				   << pair->bp << "\t"
				   << k << "\t" 
				   << pop_1 << "\t" 
				   << pop_2 << "\n";
			    delete pair;
			    continue;
			}

			pairs[pairs_key[pair->bp]] = pair;
		    }
		}

		//
		// Apply user-selected correction to the Fst values.
		//
		double correction;
		switch(fst_correction) {
		case p_value:
		    for (uint i = 0; i < pairs.size(); i++) {
			if (pairs[i] != NULL) {
			    pairs[i]->cfst       = pairs[i]->fet_p < p_value_cutoff ? pairs[i]->fst : 0;
			    pairs[i]->camova_fst = pairs[i]->fet_p < p_value_cutoff ? pairs[i]->amova_fst : 0;
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
			    pairs[i]->cfst       = pairs[i]->fet_p < correction ? pairs[i]->fst : 0;
			    pairs[i]->camova_fst = pairs[i]->fet_p < correction ? pairs[i]->amova_fst : 0;
			}
		    }
		    break;
		case no_correction:
		    for (uint i = 0; i < pairs.size(); i++) {
			if (pairs[i] != NULL) {
			    pairs[i]->cfst = pairs[i]->fst;
			    pairs[i]->camova_fst = pairs[i]->amova_fst;
			}
		    }
		    break;
		}

		//
		// If bootstrapping is enabled, record all Fst values.
		//
		if (bootstrap)
		    for (uint i = 0; i < pairs.size(); i++) {
			if (pairs[i] != NULL) {
			    fst_samples.push_back(pairs[i]->cfst);
			    allele_depth_samples.push_back(pairs[i]->alleles);
			}
		    }

		//
		// Calculate kernel-smoothed Fst values.
		//
		if (kernel_smoothed && loci_ordered) {
		    cerr << "  Generating kernel-smoothed Fst for " << it->first << ".\n";
		    kernel_smoothed_fst(pairs, weights, snp_dist);
		}
	    }

	    //
	    // If bootstrap resampling method is approximate, generate our single, empirical distribution.
	    //
	    map<int, vector<double> > approx_fst_dist;
	    if (bootstrap && bootstrap_type == bs_approx) 
		bootstrap_fst_approximate_dist(fst_samples, allele_depth_samples, weights, snp_dist, approx_fst_dist);

	    for (it = pmap->ordered_loci.begin(); it != pmap->ordered_loci.end(); it++) {
		string chr = it->first;
		vector<PopPair *> &pairs = genome_pairs[chr];

		//
		// Bootstrap resample this chromosome.
		//
		if (bootstrap && bootstrap_type == bs_exact) {
		    cerr << "  Bootstrap resampling kernel-smoothed Fst for " << it->first << ".\n";
		    bootstrap_fst(fst_samples, pairs, weights);
		}

		for (uint i = 0; i < pairs.size(); i++) {

		    if (pairs[i] == NULL)
			continue;

		    // if (pairs[i]->bp < 2110730 || pairs[i]->bp > 2110750) {
		    //     continue;
		    // }

		    //
		    // Calculate Fst P-value from approximate distribution.
		    //
		    if (bootstrap && bootstrap_type == bs_approx)
			pairs[i]->wfst_pval = bootstrap_approximate_pval(pairs[i]->snp_cnt, pairs[i]->wfst, approx_fst_dist);

		    cnt++;
		    sum += pairs[i]->camova_fst;

		    sprintf(fst_str,   "%0.10f", pairs[i]->fst);
		    sprintf(cfst_str,  "%0.10f", pairs[i]->cfst);
		    sprintf(wfst_str,  "%0.10f", pairs[i]->wfst);
		    sprintf(afst_str,  "%0.10f", pairs[i]->amova_fst);
		    sprintf(cafst_str, "%0.10f", pairs[i]->camova_fst);
		    sprintf(wafst_str, "%0.10f", pairs[i]->wamova_fst);

		    fh << batch_id          << "\t"
		       << pairs[i]->loc_id  << "\t"
		       << pop_1             << "\t"
		       << pop_2             << "\t"
		       << chr               << "\t"
		       << pairs[i]->bp      << "\t"
		       << pairs[i]->col     << "\t"
		       << pairs[i]->pi      << "\t"
		       << fst_str           << "\t"
		       << pairs[i]->fet_p   << "\t"
		       << pairs[i]->fet_or  << "\t"
		       << pairs[i]->ci_low  << "\t"
		       << pairs[i]->ci_high << "\t"
		       << pairs[i]->lod     << "\t"
		       << cfst_str          << "\t"
		       << wfst_str          << "\t"
		       << pairs[i]->wfst_pval << "\t"
		       << afst_str            << "\t"
		       << cafst_str           << "\t"
		       << wafst_str           << "\t"
		       << pairs[i]->snp_cnt   << "\n";

		    delete pairs[i];
		}
	    }
	    cerr << "Pop 1: " << pop_1 << "; Pop 2: " << pop_2 << "; mean Fst: " << (sum / cnt) << "\n";
	    means.push_back(sum / cnt);

	    cerr << "Pooled populations " << pop_1 << " and " << pop_2 << " contained: " << incompatible_loci << " incompatible loci; " 
		 << multiple_loci << " nucleotides covered by more than one RAD locus.\n";
	    fh.close();

	    if (log_fst_comp) 
		comp_fh->close();
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
	fh << "\t" << pops[i];
    fh << "\n";

    uint n = 0;
    for (uint i = 0; i < pops.size() - 1; i++) {
	fh << pops[i];

	for (uint k = 0; k <= i; k++)
	    fh << "\t";

	for (uint j = i + 1; j < pops.size(); j++) {
	    fh << "\t" << means[n];
	    n++;
	}
	fh << "\n";
    }

    fh.close();

    delete [] weights;

    return 0;
}

int
init_chr_pairs(map<string, vector<PopPair *> > &genome_chrs, string chr, map<uint, uint> &pairs_key, vector<CSLocus *> &sorted_loci)
{
    CSLocus *loc;
    int      len, bp;

    //
    // We need to create an array to store all the pair values for computing Fst. We must
    // account for positions in the genome that are covered by more than one RAD tag.
    //
    set<int> bps;

    for (uint pos = 0; pos < sorted_loci.size(); pos++) {
	loc = sorted_loci[pos];
	len = strlen(loc->con);

	for (int k = 0; k < len; k++) {
	    bp  = loc->sort_bp(k);
	    bps.insert(bp);
	}
    }

    genome_chrs[chr].resize(bps.size(), NULL);

    set<int>::iterator it;
    int i = 0;
    for (it = bps.begin(); it != bps.end(); it++) {
	pairs_key[*it] = i;
	i++;
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
	pairs[pos_c]->cfst = pairs[pos_c]->fet_p < correction ? pairs[pos_c]->fst : 0;
    }

    return 0;
}

int 
kernel_smoothed_popstats(map<int, CSLocus *> &catalog, PopMap<CSLocus> *pmap, PopSum<CSLocus> *psum, int pop_id, ofstream &log_fh) 
{
    //
    // We calculate a kernel-smoothing moving average of Pi and Fis values along each ordered chromosome.
    //
    // For each genomic region centered on a nucleotide position c, the contribution of the population 
    // genetic statistic at position p to the region average was weighted by the Gaussian function:
    //   exp( (-1 * (p - c)^2) / (2 * sigma^2))
    // 
    // In addition, we weight each position according to (n_k - 1), where n_k is the number of alleles
    // sampled at that location.
    //
    // By default, sigma = 150Kb, for computational efficiency, only calculate average out to 3sigma.
    //
    int limit = 3 * sigma;
    int alleles;

    //
    // Precalculate weights.
    //
    double *weights = calculate_weights();

    //
    // Create a list of all nucleotide positions in this population to make it easy to iterate over
    // with the sliding window.
    //
    map<string, vector<CSLocus *> >::iterator it;
    CSLocus *loc;
    LocSum  *lsum;
    int      len;

    vector<double> pi_samples;
    vector<double> fis_samples;
    vector<int>    allele_depth_samples;
    int            snp_dist[max_snp_dist] = {0};
    int            sites_per_snp = 0;
    int            tot_windows = 0;
    
    for (it = pmap->ordered_loci.begin(); it != pmap->ordered_loci.end(); it++) {

	for (uint pos = 0; pos < it->second.size(); pos++) {
	    loc  = it->second[pos];
	    len  = strlen(loc->con);
	    lsum = psum->pop(loc->id, pop_id);

	    for (int k = 0; k < len; k++)
		if (lsum->nucs[k].num_indv > 0 && bootstrap && lsum->nucs[k].pi > 0) {
		    fis_samples.push_back(lsum->nucs[k].Fis);
		    pi_samples.push_back(lsum->nucs[k].pi);
		    allele_depth_samples.push_back(lsum->nucs[k].num_indv * 2);
		}
	}
    }

    uint multiple_loci = 0;

    for (it = pmap->ordered_loci.begin(); it != pmap->ordered_loci.end(); it++) {
	vector<SumStat *> sites;

	// if (it->first != "groupI") continue;

	init_chr_sites(sites, pop_id, it->second, psum, multiple_loci, log_fh);

	cerr << "    Processing chromosome " << it->first << "\n";
	#pragma omp parallel private(alleles)
	{
	    int      dist, limit_l, limit_u;
	    uint     pos_l, pos_u;
	    double   weighted_fis, weighted_pi, sum_fis, sum_pi, final_weight_fis, final_weight_pi;
	    SumStat *c, *p;

	    pos_l = 0;
	    pos_u = 0;

	    //
	    // Center the window on each variable nucleotide position.
	    //
            #pragma omp for schedule(dynamic, 1) reduction(+:tot_windows) reduction(+:sites_per_snp)
	    for (uint pos_c = 0; pos_c < sites.size(); pos_c++) {
		c = sites[pos_c];

		if (c->pi == 0)
		    continue;

		////
		////if (c->bp != 1190640) continue;

		weighted_fis = 0.0;
		sum_fis      = 0.0;
		weighted_pi  = 0.0;
		sum_pi       = 0.0;

		limit_l = c->bp - limit > 0 ? c->bp - limit : 0;
		limit_u = c->bp + limit;

		while (pos_l < sites.size()) {
		    if (sites[pos_l] == NULL) {
			pos_l++;
		    } else {
			if (sites[pos_l]->bp < limit_l) 
			    pos_l++;
			else
			    break;
		    }
		}
		while (pos_u < sites.size()) {
		    if (sites[pos_u] == NULL) {
			pos_u++;
		    } else {
			if (sites[pos_u]->bp < limit_u)
			    pos_u++;
			else
			    break;
		    }
		}
		//if (pos_u < sites.size() && sites[pos_u]->bp > limit_u) pos_u--;

		// if (pos_u < sites.size())
		// 	cerr << "Calculating sliding window; start position: " << pos_l << ", " << sites[pos_l]->bp << "bp; end position: " 
		// 	     << pos_u << ", " << sites[pos_u]->bp << "bp; center: " 
		// 	     << pos_c << ", " << sites[pos_c]->bp << "bp; total size: " << sites[pos_u]->bp - sites[pos_l]->bp << "\n";
		int snp_cnt   = 0;
		int sites_cnt = 0;

		for (uint pos_p = pos_l; pos_p < pos_u; pos_p++) {
		    p = sites[pos_p];

		    if (p == NULL)
			continue;

		    sites_cnt++;

		    alleles = p->num_indv * 2;
		    dist    = p->bp > c->bp ? p->bp - c->bp : c->bp - p->bp;

		    if (dist > limit || dist < 0) {
                        #pragma omp critical
			{
			    cerr << "ERROR: current basepair is out of the sliding window.\n"
				 << "  Calculating sliding window; start position: " << pos_l << ", " << (sites[pos_l] == NULL ? -1 : sites[pos_l]->bp) << "bp; end position: " 
				 << pos_u << ", " << (sites[pos_u] == NULL ? -1 : sites[pos_u]->bp) << "bp; center: " 
				 << pos_c << ", " << sites[pos_c]->bp << "bp\n"
				 << "  Current position: " << pos_p << ", " << sites[pos_p]->bp << "; Dist: " << dist << "\n"
				 << "  Window positions:\n";

			    for (uint j = pos_l; j < pos_u; j++) {
				p = sites[j];
				if (p == NULL) continue;
				cerr << "    Position: " << j << "; " << p->bp << "bp\n";
			    }
			}
			continue;
		    }
		    //cerr << "Window centered on: " << c->bp << "; Examining bp " << p->bp << "; distance from center: " << dist << "\n";

		    if (p->pi > 0) {
			snp_cnt++;

			final_weight_fis = (alleles - 1) * weights[dist];
			weighted_fis    += p->Fis * final_weight_fis;
			sum_fis         += final_weight_fis;
		    }

		    final_weight_pi = (alleles - 1) * weights[dist];
		    weighted_pi    += p->pi * final_weight_pi;
		    sum_pi         += final_weight_pi;
		}

		sites_per_snp += (sites_cnt / snp_cnt);
		tot_windows++;

		if (snp_cnt < max_snp_dist) {
                    #pragma omp atomic
		    snp_dist[snp_cnt]++;
		}

		c->snp_cnt = snp_cnt;
		c->wFis    = weighted_fis / sum_fis;
		c->wPi     = weighted_pi  / sum_pi;

		if (bootstrap && bootstrap_type == bs_exact)
		    bootstrap_popstats(fis_samples, pi_samples, sites, pos_l, pos_u, weights, c);
	    }
	}
	sites.clear();
    }
    cerr << "Population " << pop_id << " contained " << multiple_loci << " nucleotides covered by more than one RAD locus.\n";

    //
    // If bootstrap resampling method is approximate, generate our single, empirical distribution.
    //
    map<int, vector<double> > approx_fis_dist;
    map<int, vector<double> > approx_pi_dist;
    if (bootstrap && bootstrap_type == bs_approx) {
    	sites_per_snp = sites_per_snp / tot_windows;

    	// cerr << "Sites per snp: " << sites_per_snp << "\n";

    	bootstrap_popstats_approximate_dist(fis_samples, pi_samples, allele_depth_samples, 
    					    weights, snp_dist, sites_per_snp, 
    					    approx_fis_dist, approx_pi_dist);

    	for (it = pmap->ordered_loci.begin(); it != pmap->ordered_loci.end(); it++) {

    	    for (uint pos = 0; pos < it->second.size(); pos++) {
    		loc  = it->second[pos];
    		len  = strlen(loc->con);
    		lsum = psum->pop(loc->id, pop_id);

    		for (int k = 0; k < len; k++)
    		    if (lsum->nucs[k].num_indv > 0 && bootstrap && lsum->nucs[k].pi > 0) {
    			//
    			// Calculate Fis/Pi p-values from approximate distribution.
    			//
    			lsum->nucs[k].wFis_pval = bootstrap_approximate_pval(lsum->nucs[k].snp_cnt, lsum->nucs[k].wFis, approx_fis_dist);
    			lsum->nucs[k].wPi_pval  = bootstrap_approximate_pval(lsum->nucs[k].snp_cnt, lsum->nucs[k].wPi, approx_pi_dist);
    		    }
    	    }
    	}
    }

    delete [] weights;

    return 0;
}

int
init_chr_sites(vector<SumStat *> &sites, 
	       int pop_id, vector<CSLocus *> &sorted_loci, PopSum<CSLocus> *psum, 
	       uint &multiple_loci, ofstream &log_fh)
{
    CSLocus *loc;
    LocSum  *lsum;
    int      len;

    //
    // We need to create an array to store all the nucleotide values for computing kernel-smoothed
    // statistics for this population. We must account for positions in the genome that are covered 
    // by more than one RAD tag.
    //
    set<int> bps;

    for (uint pos = 0; pos < sorted_loci.size(); pos++) {
	loc  = sorted_loci[pos];
	len  = strlen(loc->con);
	lsum = psum->pop(loc->id, pop_id);

	for (int k = 0; k < len; k++) {
	    if (lsum->nucs[k].num_indv > 0) 
		bps.insert(lsum->nucs[k].bp);
	}
    }

    sites.resize(bps.size(), NULL);

    //
    // Create a key describing where in the sites array to find each basepair coordinate.
    //
    map<uint, uint>    sites_key;
    set<int>::iterator it;
    int i = 0;
    for (it = bps.begin(); it != bps.end(); it++) {
	sites_key[*it] = i;
	i++;
    }

    //
    // Assign nucleotides to their proper, ordered location in the genome, 
    // checking that a site hasn't already been covered by another RAD locus.
    //
    for (uint pos = 0; pos < sorted_loci.size(); pos++) {
	loc  = sorted_loci[pos];
	len  = strlen(loc->con);
	lsum = psum->pop(loc->id, pop_id);

	for (int k = 0; k < len; k++) {
	    if (lsum->nucs[k].num_indv == 0) continue;

	    if (sites_key.count(lsum->nucs[k].bp) == 0) {
		cerr << "Error: locus " << lsum->nucs[k].loc_id << " at " << lsum->nucs[k].bp << "bp is not defined in the sites map.\n";

	    } else if (sites[sites_key[lsum->nucs[k].bp]] == NULL) {
		sites[sites_key[lsum->nucs[k].bp]] = &(lsum->nucs[k]);

	    } else {
		multiple_loci++;
		log_fh << "within_population\t"
		       << "multiple_locus\t"
		       << loc->id << "\t"
		       << loc->loc.chr << "\t"
		       << lsum->nucs[k].bp << "\t"
		       << k << "\t" 
		       << pop_id << "\t"
		       << "conflicts with locus " << sites[sites_key[lsum->nucs[k].bp]]->loc_id << "\n";
	    }
	}
    }

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
    int    pos, index_1, index_2, index_3, dist, start, end;
    int    half = sites_per_snp / 2;

    for (int i = 0; i < max_snp_dist; i++) {
	if (snp_dist[i] == 0.0) continue;

	cerr << "  Generating NULL distribution for " << i << " SNPs...\n";

        #pragma omp parallel private(poss, pos, index_1, index_2, index_3, dist, sum_fis, sum_pi, weighted_fis, weighted_pi, final_weight_fis, final_weight_pi)
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
		index_1 = (int) (fis_samples.size()    * (random() / (RAND_MAX + 1.0)));
		index_2 = (int) (pi_samples.size()     * (random() / (RAND_MAX + 1.0)));
		index_3 = (int) (allele_samples.size() * (random() / (RAND_MAX + 1.0)));
		//
		// Fill in the area around the SNP with fixed sites.
		//
		start = pos - half > 0 ? pos - half : 0;
		end   = pos + half < win_size ? pos + half : win_size;
		for (int n = start; n < end; n++) {
		    bs[n].f       = 0;
		    bs[n].pi      = 0;
		    bs[n].alleles = bs[pos].alleles;
		    poss.push_back(n);
		}
		bs[pos].f       = fis_samples[index_1];
		bs[pos].pi      = pi_samples[index_2];
		bs[pos].alleles = allele_samples[index_3];
		// cerr << "      Placing SNP at position: " << pos << "; with data from " << index_1 << " filling area from " << start << " to " << end << "\n";

		//
		// Randomly select the positions and values for each SNP to populate the window
		// 
		for (int k = 0; k < i - 1; k++) {
		    pos     = (int) (win_size * (random() / (RAND_MAX + 1.0)));
		    index_1 = (int) (fis_samples.size()    * (random() / (RAND_MAX + 1.0)));
		    index_2 = (int) (pi_samples.size()     * (random() / (RAND_MAX + 1.0)));
		    index_3 = (int) (allele_samples.size() * (random() / (RAND_MAX + 1.0)));

		    poss.push_back(pos);
		    //
		    // Fill in the area around the SNP with fixed sites.
		    //
		    start = pos - half > 0 ? pos - half : 0;
		    end   = pos + half < win_size ? pos + half : win_size;
		    for (int n = start; n < end; n++) {
			bs[n].f       = 0;
			bs[n].pi      = 0;
			bs[n].alleles = bs[pos].alleles;
			poss.push_back(n);
		    }
		    bs[pos].f       = fis_samples[index_1];
		    bs[pos].pi      = pi_samples[index_2];
		    bs[pos].alleles = allele_samples[index_3];
		    // cerr << "      Placing SNP at position: " << pos << "; with data from " << index_1 << " filling area from " << start << " to " << end << "\n";
		}

		weighted_fis = 0.0;
		sum_fis      = 0.0;
		weighted_pi  = 0.0;
		sum_pi       = 0.0;

		for (int n = 0; n < win_size; n++) {
		    if (bs[n].pi < 0.0)
			continue;
		    //
		    // Calculate weighted Fst at this position.
		    //
		    dist = bs[n].bp > bs[win_cntr].bp ? bs[n].bp - bs[win_cntr].bp : bs[win_cntr].bp - bs[n].bp;

		    final_weight_fis = (bs[n].alleles - 1) * weights[dist];
		    weighted_fis    += bs[n].f * final_weight_fis;
		    sum_fis         += final_weight_fis;

		    final_weight_pi  = (bs[n].alleles - 1) * weights[dist];
		    weighted_pi     += bs[n].pi * final_weight_pi;
		    sum_pi          += final_weight_pi;
		}

		fiss.push_back(weighted_fis / sum_fis);
		pis.push_back(weighted_pi  / sum_pi);
		// cerr << "      New weighted fis value: " << weighted_fis / sum_fis << "; size: " << fiss.size() << "\n";

		for (uint n = 0; n < poss.size(); n++) {
		    bs[poss[n]].f  = 0.0;
		    bs[poss[n]].pi = -1.0;
		}
		poss.clear();
	    }

	    #pragma omp critical
	    {
		vector<double> &f = approx_fis_dist[i];
		for (uint n = 0; n < fiss.size(); n++)
		    f.push_back(fiss[n]);
		vector<double> &p = approx_pi_dist[i];
		for (uint n = 0; n < pis.size(); n++)
		    p.push_back(pis[n]);
	    }

	    delete [] bs;
	}

	sort(approx_fis_dist[i].begin(), approx_fis_dist[i].end());
	sort(approx_pi_dist[i].begin(),  approx_pi_dist[i].end());
    }

    return 0;
}

int
bootstrap_popstats(vector<double> &fis_samples, vector<double> &pi_samples, 
		   vector<SumStat *> &sites, int pos_l, int pos_u, 
		   double *weights,
		   SumStat *c)
{
    int    size = pos_u - pos_l;
    int    dist, index;
    double weighted_fis, weighted_pi, partial_weighted_fis, partial_weighted_pi; 
    double partial_sum_fis, partial_sum_pi, sum_fis, sum_pi;
    double final_weight_fis, final_weight_pi;

    //
    // Allocate an array of bootstrap resampling objects.
    //
    BSample *bs = new BSample[size];

    //
    // Populate the BSample objects.
    //
    int j = 0;
    for (int i = pos_l; i < pos_u;  i++) {
	bs[j].f       = sites[i]->Fis;
 	bs[j].pi      = sites[i]->pi;
	bs[j].bp      = sites[i]->bp;
	bs[j].alleles = sites[i]->num_indv * 2;
	j++;
    }

    vector<double> fiss, pis;
    fiss.reserve(bootstrap_reps);
    pis.reserve(bootstrap_reps);

    // cerr << "Window starts at " << bs[0].bp << "; Window center is " << c->bp << "\n";

    //
    // Precompute the fraction of the window that will not change during resampling.
    //
    partial_weighted_fis = 0.0;
    partial_sum_fis      = 0.0;
    partial_weighted_pi  = 0.0;
    partial_sum_pi       = 0.0;

    for (j = 0; j < size; j++) {
	if (bs[j].pi > 0) continue;

	dist = bs[j].bp > c->bp ? bs[j].bp - c->bp : c->bp - bs[j].bp;

	final_weight_fis      = (bs[j].alleles - 1.0) * weights[dist];
	partial_weighted_fis += bs[j].f * final_weight_fis;
	partial_sum_fis      += final_weight_fis;

	final_weight_pi      = (bs[j].alleles - 1.0) * weights[dist];
	partial_weighted_pi += bs[j].pi * final_weight_pi;
	partial_sum_pi      += final_weight_pi;
    }

    //
    // Bootstrap this bitch.
    //
    for (int i = 0; i < bootstrap_reps; i++) {
	// cerr << "  Bootsrap rep " << i << "\n";

	weighted_fis = partial_weighted_fis;
	weighted_pi  = partial_weighted_pi;
	sum_fis      = partial_sum_fis;
	sum_pi       = partial_sum_pi;

	for (j = 0; j < size; j++) {
	    //
	    // Resample for this round of bootstrapping.
	    //
	    if (bs[j].pi == 0) 
		continue;

	    index = (int) (fis_samples.size() * (random() / (RAND_MAX + 1.0)));
	    // cerr << "      WinPos: " << j << "; Randomly selecting " << index << " out of " << pi_samples.size() << " possible values giving Pi value: " << pi_samples[index] << "\n";

	    //
	    // Calculate weighted Fis/Pi at this position.
	    //
	    dist = bs[j].bp > c->bp ? bs[j].bp - c->bp : c->bp - bs[j].bp;

	    final_weight_fis = (bs[j].alleles - 1.0) * weights[dist];
	    weighted_fis    += fis_samples[index] * final_weight_fis;
	    sum_fis         += final_weight_fis;

	    final_weight_pi = (bs[j].alleles - 1.0) * weights[dist];
	    weighted_pi    += pi_samples[index] * final_weight_pi;
	    sum_pi         += final_weight_pi;
	}

	// cerr << "    New weighted Pi value: " << weighted_pi / sum_pi << "\n";
	fiss.push_back(weighted_fis / sum_fis);
	pis.push_back(weighted_pi / sum_pi);
    }

    delete [] bs;

    //
    // Cacluate the p-value for this window based on the empirical Fst distribution.
    //
    sort(fiss.begin(), fiss.end());
    c->wFis_pval = bootstrap_pval(c->wFis, fiss);

    sort(pis.begin(), pis.end());
    c->wPi_pval = bootstrap_pval(c->wPi, pis);

    return 0;
}

int
kernel_smoothed_fst(vector<PopPair *> &pairs, double *weights, int *snp_dist) 
{
    //
    // To generate smooth genome-wide distributions of Fst, we calculate a kernel-smoothing 
    // moving average of Fst values along each ordered chromosome.
    //
    // For each genomic region centered on a nucleotide position c, the contribution of the population 
    // genetic statistic at position p to the region average was weighted by the Gaussian function:
    //   exp( (-1 * (p - c)^2) / (2 * sigma^2))
    // 
    // In addition, we weight each position according to (n_k - 1), where n_k is the number of alleles
    // sampled at that location.
    //
    // By default, sigma = 150Kb, for computational efficiency, only calculate average out to 3sigma.
    //
    #pragma omp parallel
    { 
	int      limit = 3 * sigma;
	int      dist, limit_l, limit_u;
	uint     pos_l, pos_u;
	double   weighted_fst, weighted_amova_fst, sum, final_weight;
	PopPair *c, *p;

	pos_l = 0;
	pos_u = 0;

        #pragma omp for schedule(dynamic, 1)
	for (uint pos_c = 0; pos_c < pairs.size(); pos_c++) {
	    c = pairs[pos_c];

	    if (c == NULL)
		continue;

	    weighted_fst       = 0.0;
	    weighted_amova_fst = 0.0;
	    sum                = 0.0;

	    limit_l = c->bp - limit > 0 ? c->bp - limit : 0;
	    limit_u = c->bp + limit;

	    while (pos_l < pairs.size()) {
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
	    //if (pos_u < pairs.size() && pairs[pos_u]->bp > limit_u) 
		//do { pos_u--; } while (pairs[pos_u] == NULL);

	    // cerr << "Calculating sliding window; start position: " << pos_l << ", " << pairs[pos_l]->bp << "bp; end position: " 
	    //      << pos_u << ", " << pairs[pos_u]->bp << "bp; center: " 
	    //      << pos_c << ", " << pairs[pos_c]->bp << "bp\n";
	    int snp_cnt = 0;

	    for (uint pos_p = pos_l; pos_p < pos_u; pos_p++) {
		p = pairs[pos_p];

		if (p == NULL)
		    continue;

		snp_cnt++;

		dist = p->bp > c->bp ? p->bp - c->bp : c->bp - p->bp;

		if (dist > limit || dist < 0) {
		    #pragma omp critical
		    {
			cerr << "ERROR: current basepair is out of the sliding window.\n"
			     << "  Calculating sliding window; start position: " << pos_l << ", " << (pairs[pos_l] == NULL ? -1 : pairs[pos_l]->bp) << "bp; end position: " 
			     << pos_u << ", " << (pairs[pos_u] == NULL ? -1 : pairs[pos_u]->bp) << "bp; center: " 
			     << pos_c << ", " << pairs[pos_c]->bp << "bp\n"
			     << "  Current position: " << pos_p << ", " << pairs[pos_p]->bp << "; Dist: " << dist << "\n"
			     << "  Window positions:\n";

			for (uint j = pos_l; j < pos_u; j++) {
			    p = pairs[j];
			    if (p == NULL) continue;
			    cerr << "    Position: " << j << "; " << p->bp << "bp\n";
			}
			//exit(0);
		    }
		    continue;
		}

		final_weight        = (p->alleles - 1) * weights[dist];
		weighted_fst       += p->cfst * final_weight;
		weighted_amova_fst += p->camova_fst * final_weight;
		sum                += final_weight;
	    }

	    // cerr << "Fst measure at " << c->bp << "bp has " << snp_cnt << " snps.\n";

	    if (snp_cnt < max_snp_dist) {
		#pragma omp atomic
		snp_dist[snp_cnt]++;
	    }

	    c->snp_cnt    = snp_cnt;
	    c->wfst       = weighted_fst / sum;
	    c->wamova_fst = weighted_amova_fst / sum;
	}
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
    int    pos, index_1, index_2, dist;

    for (int i = 0; i < max_snp_dist; i++) {
	if (snp_dist[i] == 0.0) continue;

	cerr << "  Generating NULL distribution for " << i << " SNPs...\n";

        #pragma omp parallel private(poss, pos, index_1, index_2, dist, sum, weighted_fst, final_weight)
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
		index_1 = (int) (fst_samples.size() * (random() / (RAND_MAX + 1.0)));
		index_2 = (int) (allele_samples.size() * (random() / (RAND_MAX + 1.0)));
		bs[pos].f       = fst_samples[index_1];
		bs[pos].alleles = allele_samples[index_2];

		//
		// Randomly select the positions and values for each SNP to populate the window
		// 
		for (int k = 0; k < i - 1; k++) {
		    pos     = (int) (win_size * (random() / (RAND_MAX + 1.0)));
		    index_1 = (int) (fst_samples.size() * (random() / (RAND_MAX + 1.0)));
		    index_2 = (int) (allele_samples.size() * (random() / (RAND_MAX + 1.0)));
		    bs[pos].f       = fst_samples[index_1];
		    bs[pos].alleles = allele_samples[index_2];
		    // cerr << "  " << j << ": Placing SNP at position: " << pos << " with data from index " << index_1 << "\n";

		    poss.push_back(pos);
		}

		weighted_fst = 0.0;
		sum          = 0.0;

		for (int n = 0; n < win_size; n++) {
		    if (bs[n].f == 0.0)
			continue;
		    //
		    // Calculate weighted Fst at this position.
		    //
		    dist = bs[n].bp > bs[win_cntr].bp ? bs[n].bp - bs[win_cntr].bp : bs[win_cntr].bp - bs[n].bp;

		    final_weight  = (bs[n].alleles - 1) * weights[dist];
		    weighted_fst += bs[n].f * final_weight;
		    sum          += final_weight;
		}

		fsts.push_back(weighted_fst / sum);
		// cerr << "    New weighted Fst value: " << weighted_fst / sum << "; size: " << fsts.size() << "\n";

		for (uint n = 0; n < poss.size(); n++)
		    bs[poss[n]].f = 0.0;
		poss.clear();
	    }

	    #pragma omp critical
	    {
		vector<double> &f = approx_fst_dist[i];
		for (uint n = 0; n < fsts.size(); n++)
		    f.push_back(fsts[n]);
	    }

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

double
bootstrap_pval(double stat, vector<double> &dist)
{
    vector<double>::iterator up;
    double pos;

    up = upper_bound(dist.begin(), dist.end(), stat);

    if (up == dist.begin())
	pos = 1;
    else if (up == dist.end())
	pos = dist.size();
    else 
	pos = up - dist.begin() + 1;

    double res = 1.0 - (pos / (double) dist.size());

    // cerr << "Generated Smoothed Fst Distribution:\n";
    // for (uint n = 0; n < dist.size(); n++)
    // 	cerr << "  n: " << n << "; Fst: " << dist[n] << "\n";

    // cerr << "Comparing Fst value: " << stat 
    // 	 << " at position " << (up - dist.begin()) << " out of " 
    // 	 << dist.size() << " positions (converted position: " << pos << "); pvalue: " << res << ".\n";

    return res;
}

int
bootstrap_fst(vector<double> &fst_samples, vector<PopPair *> &pairs, double *weights) 
{
    #pragma omp parallel
    { 
	PopPair *c;
	double final_weight, weighted_fst, sum;
	int  dist, index, limit_l, limit_u;
	int  limit = 3 * sigma;
	uint pos_l = 0;
	uint pos_u = 0;

        #pragma omp for schedule(dynamic, 1)  
	for (uint pos_c = 0; pos_c < pairs.size(); pos_c++) {
	    c = pairs[pos_c];

	    if (c == NULL)
		continue;

	    limit_l = c->bp - limit > 0 ? c->bp - limit : 0;
	    limit_u = c->bp + limit;

	    while (pos_l < pairs.size()) {
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
	    if (pos_u < pairs.size() && pairs[pos_u]->bp > limit_u) pos_u--;

	    int size = 0;
	    for (uint i = pos_l; i < pos_u;  i++) {
		if (pairs[i] == NULL) continue;
		size++;
	    }

	    //
	    // Allocate an array of bootstrap resampling objects.
	    //
	    BSample *bs = new BSample[size];

	    //
	    // Populate the BSample objects.
	    //
	    int j = 0;
	    for (uint i = pos_l; i < pos_u;  i++) {
		if (pairs[i] == NULL) continue;

		bs[j].bp      = pairs[i]->bp;
		bs[j].alleles = pairs[i]->alleles;
		j++;
	    }

	    vector<double> fsts;
	    fsts.reserve(bootstrap_reps);
	    vector<double>::iterator up;

	    // cerr << "Window starts at " << bs[0].bp << "; centered on " << c->bp << "\n";

	    //
	    // Bootstrap this bitch.
	    //
	    for (int i = 0; i < bootstrap_reps; i++) {
		// cerr << "  Bootsrap rep " << i << "\n";

		weighted_fst = 0.0;
		sum          = 0.0;

		for (j = 0; j < size; j++) {
		    //
		    // Resample for this round of bootstrapping.
		    //
		    index   = (int) (fst_samples.size() * (random() / (RAND_MAX + 1.0)));
		    bs[j].f = fst_samples[index];
		    // cerr << "      WinPos: " << j << "; Randomly selecting " << index << " out of " << fst_samples.size() << " possible values giving Fst value: " << bs[j].f << "\n";

		    //
		    // Calculate weighted Fst at this position.
		    //
		    dist = bs[j].bp > c->bp ? bs[j].bp - c->bp : c->bp - bs[j].bp;

		    final_weight  = (bs[j].alleles - 1) * weights[dist];
		    weighted_fst += bs[j].f * final_weight;
		    sum          += final_weight;
		}

		// cerr << "    New weighted Fst value: " << weighted_fst / sum << "\n";
		fsts.push_back(weighted_fst / sum);
	    }

	    //
	    // Cacluate the p-value for this window based on the empirical Fst distribution.
	    //
	    sort(fsts.begin(), fsts.end());
 	    c->wfst_pval = bootstrap_pval(c->wfst, fsts);

	    delete [] bs;
	}
    }

    return 0;
}

double *
calculate_weights() 
{
    int limit = 3 * sigma;
    //
    // Calculate weights for window smoothing operations.
    //
    // For each genomic region centered on a nucleotide position c, the contribution of the population 
    // genetic statistic at position p to the region average was weighted by the Gaussian function:
    //   exp( (-1 * (p - c)^2) / (2 * sigma^2))
    //
    double *weights = new double[limit + 1];

    for (int i = 0; i <= limit; i++)
	weights[i] = exp((-1 * pow(i, 2)) / (2 * pow(sigma, 2)));

    return weights;
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

    map<int, CSLocus *>::iterator it;
    CSLocus *loc;
    char    f[id_len], g[id_len];
    stringstream gtype_map;

    for (it = catalog.begin(); it != catalog.end(); it++) {
	loc = it->second;

	if (loc->marker.length() == 0) continue;

	string freq  = "";
	double max   = 0.0;
	int    total = 0;
	tally_haplotype_freq(loc, pmap, total, max, freq);

	sprintf(f, "%0.1f", max);
	sprintf(g, "%0.2f", loc->f);

	//
	// Record the haplotype to genotype map.
	//
	map<string, string>::iterator j;
	gtype_map.str("");
	for (j = loc->gmap.begin(); j != loc->gmap.end(); j++)
	    gtype_map << j->first << ":" << j->second << ";";

	fh << 0 << "\t" 
	   << batch_id << "\t" 
	   << loc->id << "\t" 
	   << "\t"              // Marker
	   << total << "\t"
	   << f << "\t"
	   << freq << "\t"
	   << g << "\t"
           << gtype_map.str() <<"\n";
    }

    fh.close();

    return 0;
}

int 
write_vcf(map<int, CSLocus *> &catalog, PopMap<CSLocus> *pmap, PopSum<CSLocus> *psum, map<int, string> &samples, vector<int> &sample_ids) 
{
    //
    // Write a VCF file as defined here: http://www.1000genomes.org/node/101
    //
    stringstream pop_name;
    pop_name << "batch_" << batch_id << ".vcf";
    string file = in_path + pop_name.str();

    cerr << "Writing population data to VCF file '" << file << "'\n";

    ofstream fh(file.c_str(), ofstream::out);

    if (fh.fail()) {
        cerr << "Error opening VCF file '" << file << "'\n";
	exit(1);
    }

    //
    // Load SNP data so that model likelihoods can be output to VCF file.
    //
    cerr << "Loading SNP data for " << samples.size() << " samples.\n";
    map<int, CSLocus *>::iterator cit;
    map<int, SNPRes *>::iterator sit;
    CSLocus *loc;
    Datum   *datum;

    for (uint i = 0; i < sample_ids.size(); i++) {
	map<int, SNPRes *> snpres;
	load_snp_calls(in_path + samples[sample_ids[i]], snpres);

    	for (cit = catalog.begin(); cit != catalog.end(); cit++) {
    	    loc = cit->second;
    	    datum = pmap->datum(loc->id, sample_ids[i]);

    	    if (datum != NULL && snpres.count(datum->id)) {
		for (uint j = 0; j < snpres[datum->id]->snps.size(); j++)
		    datum->snps.push_back(snpres[datum->id]->snps[j]);
		snpres[datum->id]->snps.clear();
    	    }
    	}

	for (sit = snpres.begin(); sit != snpres.end(); sit++)
	    delete sit->second;
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
    fh << "##fileformat=VCFv4.0\n"
       << "##fileDate=" << date << "\n"
       << "##source=\"Stacks v" << VERSION << "\"\n"
       << "##INFO=<ID=NS,Number=1,Type=Integer,Description=\"Number of Samples With Data\">\n"
       << "##INFO=<ID=AF,Number=.,Type=Float,Description=\"Allele Frequency\">\n"
       << "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n"
       << "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Read Depth\">\n"
       << "##FORMAT=<ID=GL,Number=.,Type=Float,Description=\"Genotype Likelihood\">\n"
       << "#CHROM" << "\t" << "POS" << "\t" << "ID" << "\t" << "REF" << "\t" << "ALT" << "\t" 
       << "QUAL" << "\t" << "FILTER" << "\t" << "INFO" << "\t" << "FORMAT";

    for (int i = 0; i < pmap->sample_cnt(); i++)
	fh << "\t" << samples[pmap->rev_sample_index(i)];
    fh << "\n";    

    map<string, vector<CSLocus *> >::iterator it;
    Datum   **d;
    LocSum  **s;
    LocTally *t;
    int       len, gt_1, gt_2;
    double    p_freq, num_indv;
    int       pop_cnt = psum->pop_cnt();
    char      p_allele, q_allele, p_str[32], q_str[32];

    for (it = pmap->ordered_loci.begin(); it != pmap->ordered_loci.end(); it++) {
	for (uint pos = 0; pos < it->second.size(); pos++) {
	    loc = it->second[pos];

	    s   = psum->locus(loc->id);
	    t   = psum->locus_tally(loc->id);
	    len = strlen(loc->con);

	    for (uint i = 0; i < loc->snps.size(); i++) {
		uint col = loc->snps[i]->col;

		num_indv = (double) t->nucs[col].num_indv;
		p_freq   = 0.0;
		for (int j = 0; j < pop_cnt; j++) {
		    //
		    // Sum the most frequent allele across populations.
		    //
		    if (s[j]->nucs[col].p_nuc == t->nucs[col].p_allele)
			p_freq += s[j]->nucs[col].p * (s[j]->nucs[col].num_indv / num_indv);
		    else 
			p_freq += (1 - s[j]->nucs[col].p) * (s[j]->nucs[col].num_indv / num_indv);
		}

		// 
		// If this site is fixed in all populations or has too many alleles don't output it.
		//
		if (t->nucs[col].allele_cnt != 2) 
		    continue;

		sprintf(p_str, "%0.3f", p_freq);
		sprintf(q_str, "%0.3f", 1 - p_freq);

		fh << loc->loc.chr << "\t" 
		   << loc->sort_bp(col) << "\t" 
		   << loc->id << "\t"
		   << t->nucs[col].p_allele << "\t"            // REFerence allele
		   << t->nucs[col].q_allele << "\t"            // ALTernate allele
		   << "."        << "\t"                       // QUAL
		   << "PASS"     << "\t"                       // FILTER
		   << "NS="      << num_indv << ";"            // INFO
		   << "AF="      << p_str << ":" << q_str << ";" << "\t" // INFO
		   << "GT:DP:GL";                              // FORMAT

		d = pmap->locus(loc->id);

		for (int j = 0; j < pmap->sample_cnt(); j++) {
		    fh << "\t";

		    if (d[j] == NULL) {
			//
			// Data does not exist.
			//
			fh << ".:0:.,.,.";
		    } else if (d[j]->model[col] == 'U') {
			//
			// Data exists, but the model call was uncertain.
			//
			fh << ".:" << d[j]->tot_depth << ":.,.,.";
		    } else {
			//
			// Tally up the nucleotide calls.
			//
			tally_observed_haplotypes(d[j]->obshap, i, p_allele, q_allele);

			if (p_allele == 0 && q_allele == 0) {
			    // More than two potential alleles.
			    fh << ".:" << d[j]->tot_depth << ":.,.,.";
			} else if (p_allele == 0) {
			    gt_1 = q_allele == t->nucs[col].p_allele ? 0 : 1;
			    fh << gt_1 << "/" << gt_1 << ":" << d[j]->tot_depth << ":.,.,.";
			} else if (q_allele == 0) {
			    gt_1 = p_allele == t->nucs[col].p_allele ? 0 : 1;
			    fh << gt_1 << "/" << gt_1 << ":" << d[j]->tot_depth << ":.,.,.";
			} else {
			    gt_1 = p_allele == t->nucs[col].p_allele ? 0 : 1;
			    gt_2 = q_allele == t->nucs[col].p_allele ? 0 : 1;
			    fh << gt_1 << "/" << gt_2 << ":" << d[j]->tot_depth;
			    //
			    // Find the heterozygous SNP call for this column and output it.
			    //
			    int  snp_index = -1;
			    uint k;
			    for (k = 0; k < d[j]->snps.size(); k++)
				if (d[j]->snps[k]->col == col) {
				    snp_index = k;
				    break;
				}
			    if (snp_index >= 0) {
				fh << ":.," << d[j]->snps[k]->lratio << ",.";
			    } else {
				cerr << "Warning, unable to locate SNP call for catalog locus " << loc->id << ", tag ID " << d[j]->id << "\n";
				fh << ":.,.,.";
			    }
			}
		    }
		}
		fh << "\n";
	    }
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
    int      len, start_index, end_index, col, pop_id;
    char     p_allele, q_allele;

    //
    // Output all the loci on the second line, comma-separated.
    //
    uint i   = 0;
    uint cnt = catalog.size();
    uint locus_index = 0;

    for (it = catalog.begin(); it != catalog.end(); it++) {
	loc = it->second;
	for (uint j = 0; j < loc->snps.size(); j++) {
	    col = loc->snps[j]->col;
	    s   = psum->locus(loc->id);
	    t   = psum->locus_tally(loc->id);

	    // 
	    // If this site is fixed in all populations or has too many alleles don't output it.
	    //
	    if (t->nucs[col].allele_cnt != 2) 
		continue;

	    fh << loc->id << "_" << col;
	    if (i <  cnt - 1) fh << ",";
	    //if (i == cnt - 1 && j < loc->snps.size() - 1) fh << ",";
	    if (write_single_snp) break;
	}
	i++;
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

	    locus_index = 0;
	    for (it = catalog.begin(); it != catalog.end(); it++) {
		loc = it->second;

		len = strlen(loc->con);
		d   = pmap->locus(loc->id);
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
		    if (write_single_snp) break;
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
    fh << "# Stacks v" << VERSION << "; " << " Structure v2.3; " << date << "\n"
       << "\t";

    map<string, vector<CSLocus *> >::iterator it;
    CSLocus  *loc;
    Datum   **d;
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
		    fh << "\t" << loc->id;
		    if (write_single_snp) 
			break;
		    else 
			fh << "_" << col;
		}
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
			if (write_single_snp) break;
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
			if (write_single_snp) break;
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
write_phase(map<int, CSLocus *> &catalog, 
	    PopMap<CSLocus> *pmap, 
	    PopSum<CSLocus> *psum, 
	    map<int, pair<int, int> > &pop_indexes, 
	    map<int, string> &samples) 
{
    //
    // Write a PHASE/fastPHASE file as defined here: http://stephenslab.uchicago.edu/software.html
    //
    // We will write one file per chromosome.
    //
    cerr << "Writing population data to PHASE/fastPHASE files...";

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
	// Tally up the number of sites
	//
	int  total_sites = 0;
	uint col;
    	for (uint pos = 0; pos < it->second.size(); pos++) {
    	    loc = it->second[pos];
	    t   = psum->locus_tally(loc->id);

    	    for (uint i = 0; i < loc->snps.size(); i++) {
    		col = loc->snps[i]->col;
		if (t->nucs[col].allele_cnt == 2) {
		    total_sites++;
		    if (write_single_snp) 
			break;
		}
	    }
	}

	//
	// Output the total number of SNP sites and the number of individuals.
	//
	fh << samples.size() << "\n"
	   << total_sites << "\n";

	//
	// Output the position of each site according to its basepair.
	//
	fh << "P";
    	for (uint pos = 0; pos < it->second.size(); pos++) {
    	    loc = it->second[pos];
	    t   = psum->locus_tally(loc->id);

    	    for (uint i = 0; i < loc->snps.size(); i++) {
    		col = loc->snps[i]->col;
		if (t->nucs[col].allele_cnt == 2) {
		    fh << " " << loc->sort_bp(col);
		    if (write_single_snp) 
			break;
		}
	    }
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
		for (uint pos = 0; pos < it->second.size(); pos++) {
		    loc = it->second[pos];

		    s = psum->locus(loc->id);
		    d = pmap->locus(loc->id);
		    t = psum->locus_tally(loc->id);

		    for (uint i = 0; i < loc->snps.size(); i++) {
			col = loc->snps[i]->col;

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
			    tally_observed_haplotypes(d[j]->obshap, i, p_allele, q_allele);

			    if (p_allele == 0 && q_allele == 0)
				gtypes << "? ";
			    else if (p_allele == 0)
				gtypes << q_allele << " ";
			    else
				gtypes << p_allele << " ";
			}
			if (write_single_snp) break;
		    }
		}
		gtypes_str = gtypes.str();
		fh << gtypes_str.substr(0, gtypes_str.length() - 1) << "\n";

		//
		// Output all the loci for this sample again, now for the q allele
		//
		gtypes.str("");
		for (uint pos = 0; pos < it->second.size(); pos++) {
		    loc = it->second[pos];

		    s = psum->locus(loc->id);
		    d = pmap->locus(loc->id);
		    t = psum->locus_tally(loc->id);

		    for (uint i = 0; i < loc->snps.size(); i++) {
			col = loc->snps[i]->col;

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
			    tally_observed_haplotypes(d[j]->obshap, i, p_allele, q_allele);

			    if (p_allele == 0 && q_allele == 0)
				gtypes << "? ";
			    else if (q_allele == 0)
				gtypes << p_allele << " ";
			    else
				gtypes << q_allele << " ";
			}
			if (write_single_snp) break;
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
		if (t->nucs[col].allele_cnt == 2) {
		    fh << chr << "\t"
		       << loc->id;
		    if (!write_single_snp) 
			fh << "_" << col;
		    fh << "\t0\t" << loc->sort_bp(col) << "\n";
		    if (write_single_snp)
			break;
		}
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

			if (write_single_snp) break;
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
    // We will write one file per chromosome.
    //
    cerr << "Writing population data to Beagle files...";

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

    for (it = pmap->ordered_loci.begin(); it != pmap->ordered_loci.end(); it++) {

	//
	// First, write a markers file containing each marker, its genomic position in basepairs
	// and the two alternative alleles at this position.
	//
	stringstream pop_name;
	pop_name << "batch_" << batch_id << "." << it->first << ".markers";
	string file = in_path + pop_name.str();

	ofstream fh(file.c_str(), ofstream::out);

	if (fh.fail()) {
	    cerr << "Error opening Beagle markers file '" << file << "'\n";
	    exit(1);
	}

	//
	// Output the header.
	//
	fh << "# Stacks v" << VERSION << "; " << " Beagle v3.3; " << date << "\n";

    	for (uint pos = 0; pos < it->second.size(); pos++) {
    	    loc = it->second[pos];
	    t   = psum->locus_tally(loc->id);

    	    for (uint i = 0; i < loc->snps.size(); i++) {
    		uint col = loc->snps[i]->col;
		if (t->nucs[col].allele_cnt == 2) {
		    fh << loc->id;
		    if (!write_single_snp) 
			fh << "_" << col;
		    fh << "\t" << loc->sort_bp(col) << "\t" 
		       << t->nucs[col].p_allele     << "\t" 
		       << t->nucs[col].q_allele     << "\n";
		    if (write_single_snp)
			break;
		}
	    }
	}

	fh.close();

	//
	// Now output the genotypes in a separate file.
	//
	pop_name.str("");
	pop_name << "batch_" << batch_id << "." << it->first << ".unphased.bgl";
	file = in_path + pop_name.str();

	fh.open(file.c_str(), ofstream::out);

	if (fh.fail()) {
	    cerr << "Error opening Beagle markers file '" << file << "'\n";
	    exit(1);
	}

	fh << "# Stacks v" << VERSION << "; " << " Beagle v3.3; " << date << "\n";

	map<int, pair<int, int> >::iterator pit;
	int  start_index, end_index, pop_id;
	char p_allele, q_allele;

	//
	// Output a list of all the samples in the data set.
	//
	for (pit = pop_indexes.begin(); pit != pop_indexes.end(); pit++) {
	    pop_id      = psum->pop_index(pit->first);
	    start_index = pit->second.first;
	    end_index   = pit->second.second;

	    fh << "I\tid";
	    for (int j = start_index; j <= end_index; j++)
		fh << "\t" << samples[pmap->rev_sample_index(j)] << "\t" << samples[pmap->rev_sample_index(j)];
	}
	fh << "\n";

	//
	// Output population IDs for each sample.
	//
	for (pit = pop_indexes.begin(); pit != pop_indexes.end(); pit++) {
	    pop_id      = pit->first;
	    start_index = pit->second.first;
	    end_index   = pit->second.second;

	    fh << "S\tid";
	    for (int j = start_index; j <= end_index; j++)
		fh << "\t" << pop_id << "\t" << pop_id;
	}
	fh << "\n";

	//
	// For each marker, output the genotypes for each sample in two successive columns.
	//
	for (pit = pop_indexes.begin(); pit != pop_indexes.end(); pit++) {
	    pop_id      = psum->pop_index(pit->first);
	    start_index = pit->second.first;
	    end_index   = pit->second.second;

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

		    fh << "M" << "\t" << loc->id;
		    if (!write_single_snp)
			fh << "_" << loc->snps[i]->col;

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
			    tally_observed_haplotypes(d[j]->obshap, i, p_allele, q_allele);

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
			    tally_observed_haplotypes(d[j]->obshap, i, p_allele, q_allele);

			    if (p_allele == 0 && q_allele == 0)
				fh << "\t" << "?";
			    else if (q_allele == 0)
				fh << "\t" << p_allele;
			    else
				fh << "\t" << q_allele;
			}
		    }
		    fh << "\n";
		    if (write_single_snp) break;
		}
	    }
	}

	fh.close();
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

    log_fh << "# Seq Pos\tLocus ID\tColumn\tPopulation\n";

    map<string, vector<CSLocus *> >::iterator it;
    CSLocus  *loc;
    Datum   **d;
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
	    d = pmap->locus(loc->id);
	    t = psum->locus_tally(loc->id);

	    for (uint i = 0; i < loc->snps.size(); i++) {
		uint col = loc->snps[i]->col;

		if (phylip_var == false) {
		    //
		    // We are looking for loci that are fixed within each population, but are 
		    // variable between one or more populations.
		    //
		    if (t->nucs[col].fixed == false || t->nucs[col].allele_cnt == 1 || t->nucs[col].pop_cnt < 2)
			continue;

		    log_fh << index << "\t" << loc->id << "\t" << col << "\t";

		    for (int j = 0; j < pop_cnt; j++) {
			pop_id = psum->rev_pop_index(j);

			if (s[j]->nucs[col].num_indv > 0) {
			    interspecific_nucs[pop_id] += s[j]->nucs[col].p_nuc;
			    log_fh << pop_id << ":" << s[j]->nucs[col].p_nuc << ",";
			} else {
			    interspecific_nucs[pop_id] += 'N';
			    log_fh << pop_id << ":N" << ",";
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
			log_fh << pop_id << ":" << nuc << ",";

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

	sprintf(id_str, "%d", pop_id);
	len = strlen(id_str);
	for (uint j = len; j < 10; j++)
	    id_str[j] = ' ';
	id_str[9] = '\0'; 

	fh << id_str << " " << interspecific_nucs[pop_id] << "\n";
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
    fh << "# Stacks v" << VERSION << "; " << " Phylip sequential; " << date << "\n";

    fh.close();
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
    char   pop_id_str[id_len];
    vector<string> parts;
    string f;
    uint   len;

    if (pmap_path.length() > 0) {
	cerr << "Parsing population map.\n";

	ifstream fh(pmap_path.c_str(), ifstream::in);

	if (fh.fail()) {
	    cerr << "Error opening population map '" << pmap_path << "'\n";
	    return 0;
	}

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
	    // <file name> <tab> <population ID>
	    //
	    parse_tsv(line, parts);

	    if (parts.size() != 2) {
		cerr << "Population map is not formated correctly: expecting two, tab separated columns, found " << parts.size() << ".\n";
		return 0;
	    }

	    strncpy(pop_id_str, parts[1].c_str(), id_len);
	    for (int i = 0; i < id_len && pop_id_str[i] != '\0'; i++)
		if (!isdigit(pop_id_str[i])) {
		    cerr << "Population map is not formated correctly: expecting numerical ID in second column, found '" << parts[1] << "'.\n";
		    return 0;
		}

	    //
	    // Test that file exists before adding to list.
	    //
	    f = in_path.c_str() + parts[0] + ".matches.tsv";
	    ifstream test_fh(f.c_str(), ifstream::in);

	    if (test_fh.fail()) {
		cerr << " Unable to find " << f.c_str() << ", excluding it from the analysis.\n";
	    } else {
		test_fh.close();
		files.push_back(make_pair(atoi(parts[1].c_str()), parts[0]));
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
	    if (pos < file.length())
		files.push_back(make_pair(1, file.substr(0, pos)));
	}

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
	cerr << "    " << it->first << ": ";
	for (int i = start; i <= end; i++) {
	    cerr << files[i].second;
	    if (i < end) cerr << ", ";
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

int parse_command_line(int argc, char* argv[]) {
    int c;
     
    while (1) {
	static struct option long_options[] = {
	    {"help",        no_argument,       NULL, 'h'},
            {"version",     no_argument,       NULL, 'v'},
            {"corr",        no_argument,       NULL, 'c'},
            {"sql",         no_argument,       NULL, 's'},
            {"vcf",         no_argument,       NULL, 'V'},
            {"structure",   no_argument,       NULL, 'S'},
            {"phase",       no_argument,       NULL, 'A'},
            {"beagle",      no_argument,       NULL, 'E'},
            {"plink",       no_argument,       NULL, 'K'},
            {"genomic",     no_argument,       NULL, 'g'},
	    {"genepop",     no_argument,       NULL, 'G'},
	    {"phylip",      no_argument,       NULL, 'Y'},
	    {"phylip_var",  no_argument,       NULL, 'L'},
	    {"window_size", required_argument, NULL, 'w'},
	    {"num_threads", required_argument, NULL, 't'},
	    {"batch_id",    required_argument, NULL, 'b'},
	    {"in_path",     required_argument, NULL, 'P'},
	    {"progeny",     required_argument, NULL, 'r'},
	    {"min_depth",   required_argument, NULL, 'm'},
	    {"renz",        required_argument, NULL, 'e'},
	    {"pop_map",     required_argument, NULL, 'M'},
	    {"whitelist",   required_argument, NULL, 'W'},
	    {"blacklist",   required_argument, NULL, 'B'},
	    {"write_single_snp",  no_argument,       NULL, 'I'},
            {"kernel_smoothed",   no_argument,       NULL, 'k'},
            {"log_fst_comp",      no_argument,       NULL, 'l'},
            {"bootstrap",         required_argument, NULL, 'O'},
	    {"bootstrap_reps",    required_argument, NULL, 'R'},
	    {"min_populations",   required_argument, NULL, 'p'},
	    {"minor_allele_freq", required_argument, NULL, 'a'},
	    {"fst_correction",    required_argument, NULL, 'f'},
	    {"p_value_cutoff",    required_argument, NULL, 'u'},
	    {0, 0, 0, 0}
	};
	
	// getopt_long stores the option index here.
	int option_index = 0;
     
	c = getopt_long(argc, argv, "hlkKSALEYVGgvcsib:p:t:o:r:M:P:m:e:W:B:I:w:a:f:p:u:R:O:", long_options, &option_index);
     
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
	case 'M':
	    pmap_path = optarg;
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
	    break;
	case 'l':
	    log_fst_comp = true;
	    break;
	case 'O':
	    bootstrap = true;
	    if (strcmp(optarg, "exact") == 0)
		bootstrap_type = bs_exact;
	    else if (strcmp(optarg, "approx") == 0)
		bootstrap_type = bs_approx;
	    else {
		cerr << "Unknown bootstrap type specified '" << optarg << "'\n";
		help();
	    }
	    break;
	case 'R':
	    bootstrap_reps = atoi(optarg);
	    break;
	case 'c':
	    corrections = true;
	    break;
	case 'i':
	    expand_id = true;
	    break;
	case 'I':
	    write_single_snp = true;
	    break;
	case 's':
	    sql_out = true;
	    break;
	case 'V':
	    vcf_out = true;
	    break;
	case 'G':
	    genepop_out = true;
	    break;
	case 'S':
	    structure_out = true;
	    break;
	case 'A':
	    phase_out = true;
	    break;
	case 'E':
	    beagle_out = true;
	    break;
	case 'K':
	    plink_out = true;
	    break;
	case 'Y':
	    phylip_out = true;
	    break;
	case 'L':
	    phylip_var = true;
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
	    if (strcmp(optarg, "p_value") == 0)
		fst_correction = p_value;
	    else if (strcmp(optarg, "bonferroni_win") == 0)
		fst_correction = bonferroni_win;
	    else if (strcmp(optarg, "bonferroni_gen") == 0)
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

    if (batch_id == 0) {
	cerr << "You must specify a batch ID.\n";
	help();
    }

    if (enz.length() > 0 && renz.count(enz) == 0) {
	cerr << "Unrecognized restriction enzyme specified: '" << enz.c_str() << "'.\n";
	help();
    }

    if (minor_allele_freq > 0) {
	if (minor_allele_freq > 1)
	    minor_allele_freq = minor_allele_freq / 100;

	if (minor_allele_freq > 0.5) {
	    cerr << "Unable to parse the minor allele frequency\n";
	    help();
	}
    }

    if (progeny_limit > 0) {
	if (progeny_limit > 1)
	    progeny_limit = progeny_limit / 100;

	if (progeny_limit > 1.0) {
	    cerr << "Unable to parse the progeny limit frequency\n";
	    help();
	}
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
	      << "  Data Filtering:\n"
	      << "    r: minimum percentage of individuals in a population required to process a locus for that population.\n"
	      << "    p: minimum number of populations a locus must be present in to process a locus.\n"
	      << "    m: specify a minimum stack depth required for individuals at a locus.\n"
	      << "    a: specify a minimum minor allele frequency required to process a nucleotide site at a locus (0 < a < 0.5).\n"
	      << "    f: specify a correction to be applied to Fst values: 'p_value', 'bonferroni_win', or 'bonferroni_gen'.\n"
	      << "    --p_value_cutoff [num]: required p-value to keep an Fst measurement (0.05 by default). Also used as base for Bonferroni correction.\n\n"
	      << "  Kernel-smoothing algorithm:\n" 
	      << "    k: enable kernel-smoothed Pi, Fis, and Fst calculations.\n"
	      << "    --window_size [num]: distance over which to average values (sigma, default 150Kb)\n\n"
	      << "  Bootstrap Resampling:\n" 
	      << "    --bootstrap [exact|approx]: enable bootstrap resampling for population statistics (reference genome required).\n"
	      << "    --bootstrap_reps [num]: number of bootstrap resamplings to calculate (default 100).\n\n"
	      << "  File ouput options:\n"
	      << "    --genomic: output each nucleotide position (fixed or polymorphic) in all population members to a file.\n"
	      << "    --vcf: output results in Variant Call Format (VCF).\n"
	      << "    --genepop: output results in GenePop format.\n"
	      << "    --structure: output results in Structure format.\n"
	      << "    --phase: output genotypes in PHASE/fastPHASE format.\n"
	      << "    --beagle: output genotypes in Beagle format.\n"
	      << "    --plink: output genotypes in PLINK format.\n"
	      << "    --phylip: output nucleotides that are fixed-within, and variant among populations in Phylip format for phylogenetic tree construction.\n"
	      << "      --phylip_var: include variable sites in the phylip output encoded using IUPAC notation.\n"
	      << "    --write_single_snp: write only the first SNP per locus in Genepop and Structure outputs.\n\n"
	      << "  Debugging:\n"
	      << "    --log_fst_comp: log components of Fst calculations to a file.\n";

    exit(0);
}
