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
int       num_threads = 1;
int       batch_id    = 0;
bool      sql_out     = false;
string    in_path;
string    out_path;
string    out_file;

int main (int argc, char* argv[]) {

    parse_command_line(argc, argv);

    //
    // Set the number of OpenMP parallel threads to execute.
    //
    #ifdef _OPENMP
    omp_set_num_threads(num_threads);
    #endif

    vector<pair<int, string> > files;
    if (!build_file_list(files))
	exit(1);

    //
    // Open the log file.
    //
    stringstream log;
    log << "batch_" << batch_id << ".rxstacks.log";
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
    if ((res = load_loci(catalog_file.str(), catalog, false)) == 0) {
    	cerr << "Unable to load the catalog '" << catalog_file.str() << "'\n";
     	return 0;
    }

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

    cerr << "Loading model outputs for " << sample_ids.size() << " samples, " << catalog.size() << " loci.\n";
    map<int, CSLocus *>::iterator it;
    map<int, ModRes *>::iterator mit;
    Datum   *d;
    CSLocus *loc;

    //
    // Load the output from the SNP calling model for each individual at each locus.
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

    //
    // Idenitfy polymorphic loci, tabulate haplotypes present.
    //
    tabulate_haplotypes(catalog, pmap);

    //
    // Output a list of heterozygous loci and the associate haplotype frequencies.
    //
    if (sql_out)
	write_sql(catalog, pmap);

    log_fh.close();

    return 0;
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

int build_file_list(vector<pair<int, string> > &files) {
    char   line[max_len];
    vector<string> parts;
    string f;
    uint   len;

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

bool hap_compare(pair<string, int> a, pair<string, int> b) {
    return (a.second > b.second);
}

int parse_command_line(int argc, char* argv[]) {
    int c;
     
    while (1) {
	static struct option long_options[] = {
	    {"help",        no_argument,       NULL, 'h'},
            {"version",     no_argument,       NULL, 'v'},
            {"sql",         no_argument,       NULL, 's'},
	    {"num_threads", required_argument, NULL, 't'},
	    {"batch_id",    required_argument, NULL, 'b'},
	    {"in_path",     required_argument, NULL, 'P'},
	    {0, 0, 0, 0}
	};
	
	// getopt_long stores the option index here.
	int option_index = 0;
     
	c = getopt_long(argc, argv, "hvst:b:P:", long_options, &option_index);
     
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
	case 's':
	    sql_out = true;
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

    return 0;
}

void version() {
    std::cerr << "rxstacks " << VERSION << "\n\n";

    exit(0);
}

void help() {
    std::cerr << "rxstacks " << VERSION << "\n"
              << "rxstackss -b batch_id -P path [-t threads] [-v] [-h]" << "\n"
	      << "  b: Batch ID to examine when exporting from the catalog.\n"
	      << "  P: path to the Stacks output files.\n"
	      << "  t: number of threads to run in parallel sections of code.\n"
	      << "  v: print program version." << "\n"
	      << "  h: display this help messsage." << "\n\n";

    exit(0);
}
