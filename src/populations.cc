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
double    sigma           = 150000;
int       progeny_limit   = 1;
bool      corrections     = false;
bool      expand_id       = false;
bool      sql_out         = false;
bool      genomic_out     = false;
bool      kernel_fst      = false;
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
    map<int, CLocus *> catalog;
    int res;
    catalog_file << in_path << "batch_" << batch_id << ".catalog";
    if ((res = load_loci(catalog_file.str(), catalog, false)) == 0) {
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
    order_unordered_loci(catalog);

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
    	psum->add_population(catalog, pmap, pop_id, start_index, end_index, log_fh);
    }

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
    // Output the observed haplotypes.
    //
    write_generic(catalog, pmap, samples, false);

    //
    // Output the locus-level summary statistics.
    //
    write_summary_stats(files, pop_indexes, catalog, pmap, psum);

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

int order_unordered_loci(map<int, CLocus *> &catalog) {
    map<int, CLocus *>::iterator it;
    CLocus *loc;
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
	return 0;

    cerr << "Catalog is not reference aligned, arbitrarily ordering catalog loci.\n";

    uint bp = 1;
    for (it = catalog.begin(); it != catalog.end(); it++) {
	loc = it->second;
	loc->loc.chr = new char[3];
	strcpy(loc->loc.chr, "un");
	loc->loc.bp  = bp;

	bp += strlen(loc->con);
    }
    

    return 0;
}

int tabulate_haplotypes(map<int, CLocus *> &catalog, PopMap<CLocus> *pmap) {
    map<int, CLocus *>::iterator it;
    vector<char *>::iterator hit;
    Datum  **d;
    CLocus  *loc;

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

int create_genotype_map(CLocus *locus, PopMap<CLocus> *pmap) {
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

int call_population_genotypes(CLocus *locus, 
			      PopMap<CLocus> *pmap) {
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

int tally_haplotype_freq(CLocus *locus, PopMap<CLocus> *pmap,
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

    uint  rcnt = enz.length() ? renz_cnt[enz] : 0;
    uint  rlen = enz.length() ? renz_len[enz] : 0;
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
		    map<int, CLocus *> &catalog, PopMap<CLocus> *pmap, PopSum<CLocus> *psum) 
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
	fh << "# " << pit->first << "\t";
	for (int i = start; i <= end; i++) {
	    fh << files[i].second;
	    if (i < end) fh << ",";
	}
	fh << "\n";
    }

    cerr << "Writing " << catalog.size() << " loci to summary statistics file, '" << file << "'\n";

    map<string, vector<CLocus *> >::iterator it;
    CLocus  *loc;
    LocSum **s;
    int      len;
    int      pop_cnt = psum->pop_cnt();
    char     fisstr[32];
    bool     fixed, plimit;

    fh << "# Locus ID" << "\t"
       << "Chr"      << "\t"
       << "BP"       << "\t"
       << "Pop ID"   << "\t"
       << "N"        << "\t"
       << "P"        << "\t"
       << "Obs Het"  << "\t"
       << "Obs Hom"  << "\t"
       << "Exp Het"  << "\t"
       << "Exp Hom"  << "\t"
       << "Pi"       << "\t"
       << "Fis"      << "\n";

    for (it = pmap->ordered_loci.begin(); it != pmap->ordered_loci.end(); it++) {
	for (uint pos = 0; pos < it->second.size(); pos++) {
	    loc = it->second[pos];

	    s = psum->locus(loc->id);
	    len = strlen(loc->con);

	    for (int i = 0; i < len; i++) {

		// 
		// If this site is fixed in all populations, don't output it.
		//
		fixed  = true;
		plimit = true;
		for (int j = 0; j < pop_cnt; j++) {
		    if (s[j]->nucs[i].num_indv < progeny_limit) 
			plimit = false;
		    if (s[j]->nucs[i].pi != 0) 
			fixed = false;
		}

		if (!fixed && plimit) {
		    for (int j = 0; j < pop_cnt; j++) {

			sprintf(fisstr, "%0.10f", s[j]->nucs[i].Fis);

			fh << loc->id << "\t"
			   << loc->loc.chr << "\t"
			   << loc->loc.bp + i << "\t"
			   << psum->rev_pop_index(j) << "\t"
			   << s[j]->nucs[i].num_indv << "\t"
			   << s[j]->nucs[i].p << "\t"
			   << s[j]->nucs[i].obs_het << "\t"
			   << s[j]->nucs[i].obs_hom << "\t"
			   << s[j]->nucs[i].exp_het << "\t"
			   << s[j]->nucs[i].exp_hom << "\t"
			   << s[j]->nucs[i].pi      << "\t"
			   << fisstr << "\n";
		    }
		}
	    }
	}
    }

    return 0;
}

int 
write_fst_stats(vector<pair<int, string> > &files, map<int, pair<int, int> > &pop_indexes, 
		map<int, CLocus *> &catalog, PopMap<CLocus> *pmap, PopSum<CLocus> *psum, ofstream &log_fh) 
{
    //
    // We want to iterate over each pair of populations and calculate Fst at each 
    // nucleotide of each locus.
    //
    vector<int> pops;
    map<int, pair<int, int> >::iterator pit;
    for (pit = pop_indexes.begin(); pit != pop_indexes.end(); pit++)
	pops.push_back(pit->first);

    if (pops.size() == 1) return 0;

    for (uint i = 0; i < pops.size(); i++) {
	for (uint j = i + 1; j < pops.size(); j++) {
	    int pop_1 = pops[i];
	    int pop_2 = pops[j];

	    int incompatible_loci = 0;

	    stringstream pop_name;
	    pop_name << "batch_" << batch_id << ".fst_" << pop_1 << "-" << pop_2 << ".tsv";

	    string file = in_path + pop_name.str();
	    ofstream fh(file.c_str(), ofstream::out);

	    if (fh.fail()) {
		cerr << "Error opening generic output file '" << file << "'\n";
		exit(1);
	    }

	    cerr << "Calculating Fst for populations " << pop_1 << " and " << pop_2 << " and writing it to file, '" << file << "'\n";

	    fh << "# Locus ID"   << "\t"
	       << "Chr"          << "\t"
	       << "BP"           << "\t"
	       << "Column"       << "\t"
	       << "Overall Pi"   << "\t"
	       << "Fst"          << "\t"
	       << "Smoothed Fst" << "\n";

	    map<string, vector<CLocus *> >::iterator it;
	    CLocus  *loc;
	    PopPair *pair;
	    int      len;
	    char     fst_str[32], wfst_str[32];

	    for (it = pmap->ordered_loci.begin(); it != pmap->ordered_loci.end(); it++) {

		vector<PopPair *> pairs;

		for (uint pos = 0; pos < it->second.size(); pos++) {
		    loc = it->second[pos];
		    len = strlen(loc->con);

		    for (int k = 0; k < len; k++) {

			pair = psum->Fst(loc->id, pop_1, pop_2, k, progeny_limit);

			//
			// Locus is incompatible, log this position.
			//
			if (pair == NULL) {
			    pairs.push_back(NULL);
			    incompatible_loci++;
			    log_fh << loc->id << "\t"
				   << loc->loc.chr << "\t"
				   << loc->loc.bp + k << "\t"
				   << k << "\t" 
				   << pop_1 << "\t" 
				   << pop_2 << "\n";
			    continue;
			}

			//
			// Locus is fixed in both populations, or was only found in one population.
			//
			if (pair->pi == 0) {
			    delete pair;
			    pairs.push_back(NULL);
			    continue;
			}

			pair->bp = loc->loc.bp + k;
			pairs.push_back(pair);
		    }
		}

		//
		// Calculate kernel-smoothed Fst values.
		//
		if (kernel_fst) {
		    cerr << "  Generating kernel-smoothed Fst for " << it->first << ".\n";
		    kernel_smoothed_fst(pairs);
		}

		uint i = 0;
		for (uint pos = 0; pos < it->second.size(); pos++) {
		    loc = it->second[pos];
		    len = strlen(loc->con);

		    for (int k = 0; k < len; k++) {

			if (pairs[i] == NULL) {
			    i++;
			    continue;
			}

			sprintf(fst_str,  "%0.10f", pairs[i]->fst);
			sprintf(wfst_str, "%0.10f", pairs[i]->wfst);

			fh << loc->id << "\t"
			   << loc->loc.chr << "\t"
			   << loc->loc.bp + k << "\t"
			   << k << "\t"
			   << pairs[i]->pi << "\t"
			   << fst_str << "\t"
			   << wfst_str << "\n";

			delete pairs[i];
			i++;
		    }
		}
	    }
	    cerr << "Pooled populations " << pop_1 << " and " << pop_2 << " contained " << incompatible_loci << " incompatible loci.\n";
	    fh.close();
	}
    }

    return 0;
}

int
kernel_smoothed_fst(vector<PopPair *> &pairs) {
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
    int      limit = 3 * sigma;
    int      dist;
    double   weighted_fst, sum, final_weight;
    PopPair *c, *p;

    //
    // Precalculate weights.
    //
    double *weights = new double[limit + 1];
    for (int i = 0; i <= limit; i++)
	weights[i] = exp((-1 * pow(i, 2)) / (2 * pow(sigma, 2)));

    for (uint pos_c = 0; pos_c < pairs.size(); pos_c++) {
	c = pairs[pos_c];

	if (c == NULL)
	    continue;

	weighted_fst = 0.0;
	sum          = 0.0;

	for (uint pos_p = 0; pos_p < pairs.size(); pos_p++) {
	    p = pairs[pos_p];

	    if (p == NULL)
		continue;

	    dist = p->bp > c->bp ? p->bp - c->bp : c->bp - p->bp;
	    if (dist > limit)
		continue;

	    final_weight  = (p->alleles - 1) * weights[dist];
	    weighted_fst += p->fst * final_weight;
	    sum          += final_weight;
	}

	c->wfst = weighted_fst / sum;
    }

    return 0;
}

int
write_generic(map<int, CLocus *> &catalog, PopMap<CLocus> *pmap, 
	      map<int, string> &samples, bool write_gtypes)
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
    fh << "Cnt\t";

    map<int, string>::iterator s;
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

int write_sql(map<int, CLocus *> &catalog, PopMap<CLocus> *pmap) {

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

    if (pmap_path.length() > 0) {
	ifstream fh(pmap_path.c_str(), ifstream::in);

	if (fh.fail()) {
	    cerr << "Error opening population map '" << pmap_path << "'\n";
	    exit(1);
	}

	while (fh.good()) {
	    fh.getline(line, max_len);

	    if (strlen(line) == 0) continue;
	    if (line[0] == '#') continue;

	    //
	    // Parse the population map, we expect:
	    // <file name> <tab> <population ID>
	    //
	    parse_tsv(line, parts);

	    files.push_back(make_pair(atoi(parts[1].c_str()), parts[0]));
	}

	fh.close();
    } else {
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
    if (a.first == b.first)
	return (a.second < b.second);
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
            {"genomic",     no_argument,       NULL, 'g'},
            {"kernel_fst",  no_argument,       NULL, 'k'},
	    {"window_size", required_argument, NULL, 'w'},
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
     
	c = getopt_long(argc, argv, "hkgvcsib:p:t:o:r:M:P:m:e:W:B:w:", long_options, &option_index);
     
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
	case 'k':
	    kernel_fst = true;
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

    return 0;
}

void version() {
    std::cerr << "populations " << VERSION << "\n\n";

    exit(0);
}

void help() {
    std::cerr << "populations " << VERSION << "\n"
              << "populations -b batch_id -P path -M path [-r min] [-m min] [-g] [-B blacklist] [-W whitelist] [-s] [-e renz] [-v] [-h]" << "\n"
	      << "  b: Batch ID to examine when exporting from the catalog.\n"
	      << "  P: path to the Stacks output files.\n"
	      << "  M: path to the population map, a tab separated file describing which individuals belong in which population.\n"
	      << "  r: minimum number of individuals required to process a locus.\n"
	      << "  m: specify a minimum stack depth required before exporting a locus in a particular individual.\n"
	      << "  s: output a file to import results into an SQL database.\n"
	      << "  g: output each nucleotide position in all population members to a file.\n"
	      << "  B: specify a file containing Blacklisted markers to be excluded from the export.\n"
	      << "  W: specify a file containign Whitelisted markers to include in the export.\n"
	      << "  e: restriction enzyme, required if generating 'genomic' output.\n"
	      << "  v: print program version." << "\n"
	      << "  h: display this help messsage." << "\n\n"
	      << "  Kernel-smoothing algorithm:\n" 
	      << "   k: enable kernel-smoothed Fst calculation.\n"
	      << "   --window_size: distance over which to average values (sigma, default 150Kb)\n";

    exit(0);
}
