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
    if ((res = load_loci(catalog_file.str(), catalog, false, false)) == 0) {
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

		 prune_reads(cloc, loc);
	     }
	}

	//
	// Re-execute model across locus nucleotides.
	//

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

	break;
    }

    // //
    // // Idenitfy polymorphic loci, tabulate haplotypes present.
    // //
    // tabulate_haplotypes(catalog, pmap);

    log_fh.close();

    return 0;
}


int
prune_reads(CSLocus *cloc, Locus *loc)
{
    //
    // 1. Identify nucleotides in the locus for which model calls could not 
    //    be made.
    // 2. Check these nucleotides in the catalog to identify which alleles 
    //    are present in the wider population. 
    // 3. Remove reads from locus that contain alleles that are not present 
    //    in the wider population.
    //

    set<char> nucs;
    set<char>::iterator it;

    set<int>  read_blacklist;
    set<int>::iterator  it_bl;

    int  sample_col  = 0;
    int  catalog_col = 0;
    uint j           = 0;

    for (uint i = 0; i < loc->snps.size(); i++) {
	if (loc->snps[i]->type == snp_type_unk && 
	    loc->snps[i]->rank_1 != 'N') {
	    cerr << "  We have found an unknown SNP call in tag " << loc->id << " at position " << i << "; col: " << loc->snps[i]->col << "\n";

	    cerr << "    Sample has nucleotides: '" << loc->snps[i]->rank_1 << "' and '" << loc->snps[i]->rank_2 << "'\n";

	    //
	    // Advance to the appropriate catalog locus column.
	    //
	    while (j < cloc->snps.size() && catalog_col < sample_col) {
		catalog_col = cloc->snps[j]->col;

		if (catalog_col < sample_col)
		    j++;
	    }

	    if (catalog_col != sample_col) {
		//
		// If there is no SNP object in the catalog locus to represent this site
		// then we assume it is a homozygous site in the population.
		//
		nucs.insert(cloc->con[sample_col]);
	    } else {
		//
		// Otherwise it is a heterozygous site in the catalog and we can check the 
		// alleles in the catalog against those nucleotides at this site.
		//
		nucs.insert(cloc->snps[j]->rank_1);
		if (cloc->snps[j]->rank_2 != 0) nucs.insert(cloc->snps[j]->rank_2);
		if (cloc->snps[j]->rank_3 != 0) nucs.insert(cloc->snps[j]->rank_3);
		if (cloc->snps[j]->rank_4 != 0) nucs.insert(cloc->snps[j]->rank_4);
	    }

	    cerr << "    Catalog has nucleotides: '";
	    for (it = nucs.begin(); it != nucs.end(); it++)
		cerr << *it << ", ";
	    cerr << "\n";

	    //
	    // Mark reads for removal that contain alleles not present in the rest of 
	    // the population.
	    //
	    for (uint k = 0; k < loc->reads.size(); k++) {
		if (nucs.count(loc->reads[k][sample_col]) == 0 &&
		    loc->reads[k][sample_col] != 'N') 
		    read_blacklist.insert(k);
	    }
	}

	sample_col++;
	nucs.clear();
    }

    if (read_blacklist.size() > 0) {
	cerr << "    Reads marked for removal: ";
	for (it_bl = read_blacklist.begin(); it_bl != read_blacklist.end(); it_bl++)
	    cerr << *it_bl << ", ";
	cerr << "\n";

	//
	// Remove the reads.
	//
	vector<char *> reads, read_ids;
	for (uint k = 0; k < loc->reads.size(); k++) {
	    if (read_blacklist.count(k) > 0) {
		delete [] loc->reads[k];
		delete [] loc->comp[k];
	    } else {
		reads.push_back(loc->reads[k]);
		read_ids.push_back(loc->comp[k]);
	    }
	}

	loc->reads = reads;
	loc->comp  = read_ids;
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
		 << tag_1->snps[j]->rank_2 << "\t\t\n";
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

bool hap_compare(pair<string, int> a, pair<string, int> b) {
    return (a.second > b.second);
}

int parse_command_line(int argc, char* argv[]) {
    int c;
     
    while (1) {
	static struct option long_options[] = {
	    {"help",        no_argument,       NULL, 'h'},
            {"version",     no_argument,       NULL, 'v'},
	    {"num_threads", required_argument, NULL, 't'},
	    {"batch_id",    required_argument, NULL, 'b'},
	    {"in_path",     required_argument, NULL, 'P'},
	    {"outpath",     required_argument, NULL, 'o'},
	    {0, 0, 0, 0}
	};
	
	// getopt_long stores the option index here.
	int option_index = 0;
     
	c = getopt_long(argc, argv, "hvo:t:b:P:", long_options, &option_index);
     
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
	      << "  h: display this help messsage." << "\n\n";

    exit(0);
}
