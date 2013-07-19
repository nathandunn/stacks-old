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
    string log_path = out_path + log.str();
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

    int      catalog_id, sample_id, tag_id;
    string   file;
    Datum   *d;
    Locus   *loc;
    CSLocus *cloc;

    //
    // Process samples matched to the catalog, one by one.
    //
    for (uint i = 0; i < catalog_matches.size(); i++) {
	sample_id = catalog_matches[i][0]->sample_id;
	file      = samples[sample_id];

	cerr << "Making corrections to sample " << file << "...\n";

	map<int, Locus *> stacks;
	int res;
	if ((res = load_loci(in_path + file, stacks, true, true)) == 0) {
	    cerr << "Unable to load sample file '" << file << "'\n";
	    return 0;
	}

	set<pair<int, int> > processed;

	for (uint j = 0; j < catalog_matches[i].size(); j++) {
	     catalog_id = catalog_matches[i][j]->cat_id;
	     sample_id  = catalog_matches[i][j]->sample_id;
	     tag_id     = catalog_matches[i][j]->tag_id;

	     if (catalog.count(catalog_id) == 0) continue;

	     //
	     // There are multiple matches per stack, but we only need to process
	     // each stack once to make corrections.
	     //
	     if (processed.count(make_pair(catalog_id, tag_id)) == 0) {
		 processed.insert(make_pair(catalog_id, tag_id));

		 d = pmap->datum(catalog_id, sample_id);

		 if (d == NULL) continue;

		 cloc = catalog[catalog_id];
		 loc  = stacks[tag_id];

		 prune_nucleotides(cloc, loc, log_fh);
	     }
	}

	//
	// Rewrite stacks, model outputs, and haplotypes.
	//
	write_results(file, stacks);

	//
	// Free up memory
	//
	map<int, Locus *>::iterator it;
	for (it = stacks.begin(); it != stacks.end(); it++)
	    delete it->second;
    }

    log_fh.close();

    return 0;
}


int
prune_nucleotides(CSLocus *cloc, Locus *loc, ofstream &log_fh)
{
    map<char, int>      nucs;
    set<char>           cnucs;
    set<char>::iterator it;
    set<int>            rows;

    for (uint i = 0; i < loc->snps.size() && i < cloc->snps.size(); i++) {
	//
	// replace this with two cases: either their is an unknown call in locus, and no data in catalog, allowing us 
	// to look for a homozygote, or, there is a snp in the catalog and any state in the locus.
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
	    invoke_model(loc, i, nucs, log_fh);
	}

	nucs.clear();
	cnucs.clear();
    }

    //
    // Re-call alleles.
    //
    loc->alleles.clear();

    call_alleles(loc, rows);

    return 0;
}

int 
invoke_model(Locus *loc, int col, map<char, int> &nucs, ofstream &log_fh) 
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

    log_model_calls(log_fh, loc);

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
	// If a read had a nucleotide not present in the catalot, do not call
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
log_model_calls(ofstream &log_fh, Locus *loc)
{
    //
    // Log model call changes
    //
    for (uint j = 0; j < loc->snps.size(); j++) {
	switch(loc->model[j]) {
	case 'U':
	    switch(loc->snps[j]->type) {
	    case snp_type_het:
		log_fh << loc->sample_id << "\t" << loc->id << "\t" << loc->snps[j]->col << "\t" << 'U' << "\t" << 'E' << "\n";
		break;
	    case snp_type_hom:
		log_fh << loc->sample_id << "\t" << loc->id << "\t" << loc->snps[j]->col << "\t" << 'U' << "\t" << 'O' << "\n";
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
		log_fh << loc->sample_id << "\t" << loc->id << "\t" << loc->snps[j]->col << "\t" << 'E' << "\t" << 'O' << "\n";
		break;
	    case snp_type_unk:
	    default:
		log_fh << loc->sample_id << "\t" << loc->id << "\t" << loc->snps[j]->col << "\t" << 'E' << "\t" << 'U' << "\n";
		break;
	    }
	    break;
	case 'O':
	default:
	    switch(loc->snps[j]->type) {
	    case snp_type_het:
		log_fh << loc->sample_id << "\t" << loc->id << "\t" << loc->snps[j]->col << "\t" << 'O' << "\t" << 'E' << "\n";
		break;
	    case snp_type_hom:
		break;
	    case snp_type_unk:
	    default:
		log_fh << loc->sample_id << "\t" << loc->id << "\t" << loc->snps[j]->col << "\t" << 'O' << "\t" << 'U' << "\n";
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
	     << "consensus\t" << "\t\t" 
	     << tag_1->con << "\t" 
	     << "0" << "\t" 
	     << "0" << "\t"
	     << "0" << "\n";

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
		 << tag_1->reads[j]     << "\t\t\t\n";
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

    cerr << "  Wrote " << wrote << " loci.\n";

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
parse_command_line(int argc, char* argv[]) 
{
    int c;
     
    while (1) {
	static struct option long_options[] = {
	    {"help",        no_argument,       NULL, 'h'},
            {"version",     no_argument,       NULL, 'v'},
	    {"num_threads", required_argument, NULL, 't'},
	    {"batch_id",    required_argument, NULL, 'b'},
	    {"in_path",     required_argument, NULL, 'P'},
	    {"outpath",     required_argument, NULL, 'o'},
	    {"model_type",   required_argument, NULL, 'T'},
	    {"bound_low",    required_argument, NULL, 'L'},
	    {"bound_high",   required_argument, NULL, 'U'},
	    {"alpha",        required_argument, NULL, 'A'},
	    {0, 0, 0, 0}
	};
	
	// getopt_long stores the option index here.
	int option_index = 0;
     
	c = getopt_long(argc, argv, "hvo:t:b:P:T:L:U:A:", long_options, &option_index);
     
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

    return 0;
}

void version() {
    std::cerr << "rxstacks " << VERSION << "\n\n";

    exit(0);
}

void help() {
    std::cerr << "rxstacks " << VERSION << "\n"
              << "rxstackss -b batch_id -P path [-o path] [-t threads] [-v] [-h]" << "\n"
	      << "  b: Batch ID to examine when exporting from the catalog.\n"
	      << "  P: path to the Stacks output files.\n"
	      << "  o: output path to write results.\n"
	      << "  t: number of threads to run in parallel sections of code.\n"
	      << "  v: print program version." << "\n"
	      << "  h: display this help messsage." << "\n\n"
	      << "  Model options:\n" 
	      << "    --model_type <type>: either 'snp' (default), 'bounded', or 'fixed'\n"
	      << "    For the SNP or Bounded SNP model:\n"
	      << "      --alpha <num>: chi square significance level required to call a heterozygote or homozygote, either 0.1 (default), 0.05, 0.01, or 0.001.\n"
	      << "    For the Bounded SNP model:\n"
	      << "      --bound_low <num>: lower bound for epsilon, the error rate, between 0 and 1.0 (default 0).\n"
	      << "      --bound_high <num>: upper bound for epsilon, the error rate, between 0 and 1.0 (default 1).\n";
    exit(0);
}
