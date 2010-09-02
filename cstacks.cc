// -*-mode:c++; c-style:k&r; c-basic-offset:4;-*-
//
// cstacks -- Create a catalog of Stacks.
//
// Julian Catchen
// jcatchen@uoregon.edu
// University of Oregon
//

#include "cstacks.h"

// Global variables to hold command-line options.
queue<pair<int, string> > samples;
string  out_path;
int     batch_id     = 0;
int     mult_matches = 0;
searcht search_type = sequence;
int     num_threads  = 1;

int main (int argc, char* argv[]) {

    parse_command_line(argc, argv);

    //
    // Set the number of OpenMP parallel threads to execute.
    //
    omp_set_num_threads(num_threads);

    map<int, CLocus *> catalog;
    map<int, QLocus *> sample;

    pair<int, string> s = samples.front();
    samples.pop();

    cerr << "Initializing catalog...\n";
    if (!initialize_catalog(s, catalog)) {
        cerr << "Failed to initialize the catalog.\n";
        return 1;
    }

    int i = 2;
    while (!samples.empty()) {
	cerr << "Processing sample " << i << "\n";

	s = samples.front();
	samples.pop();

	if (!load_loci(s.second, sample)) {
            cerr << "Failed to load sample " << i << "\n";
            continue;
        }
	//dump_loci(sample);

        if (search_type == sequence) {
            cerr << "Searching for sequence matches...\n";
            find_matches_by_sequence(catalog, sample);
        } else if (search_type == genomic_loc) {
            cerr << "Searching for matches by genomic location...\n";
            find_matches_by_genomic_loc(catalog, sample);
        }

	cerr << "Merging matches into catalog...\n";
	merge_matches(catalog, sample, s);

	i++;
    }

    cerr << "Writing catalog...\n";
    write_catalog(catalog);

    return 0;
}

int merge_matches(map<int, CLocus *> &catalog, map<int, QLocus *> &sample, pair<int, string> &sample_file) {
    map<int, QLocus *>::iterator i;
    vector<pair<int, allele_type> >::iterator mat_it;
    CLocus *ctag;
    QLocus *utag;

    for (i = sample.begin(); i != sample.end(); i++) {
	utag = i->second;

	//
	// Emit a warning if the sample 1 tag matches more than one tag in sample 2.
	//
	set<int> local_matches;
	set<int>::iterator j;
	for (mat_it = utag->matches.begin(); mat_it != utag->matches.end(); mat_it++)
	    local_matches.insert(mat_it->first);

        //
        // If this stack didn't match an existing catalog stack, add this stack to the 
        // catalog as a new stack.
        //
	if (local_matches.size() == 0) {
	    add_unique_tag(sample_file, catalog, utag);
	    continue;
	}

	if (local_matches.size() > 1) {
	    cerr << 
		"  Warning: sample " << sample_file.second << ", tag " << utag->id << 
		", matches more than one tag in the catalog: ";
	    for (j = local_matches.begin(); j != local_matches.end(); j++)
		cerr << *j << " ";
	    cerr << "\n";

	    //
	    // Don't record matches to multiple catalog entries unless instructed
	    // to do so by the command line option.
	    //
	    if (!mult_matches) continue;
	}

	//
	// Record the tag in the second sample that has been matched to one in the first sample
	// so it will not need to be output again below.
	//
	j = local_matches.begin();
        ctag = catalog[(*j)];

        if (ctag == NULL) 
            cerr << "  Unable to locate catalog tag " << *j << "\n";

	//
	// Merge the SNPs and alleles assigned to the two matched tags
	//
	if (!ctag->merge_snps(utag)) {
	    cerr << "Error merging " << sample_file.second << ", tag " << utag->id <<
		" with catalog tag " << ctag->id << "\n";
	}

	ctag->sources.push_back(make_pair(sample_file.first, utag->id));
	ctag->populate_alleles();
    }

    return 0;
}

int add_unique_tag(pair<int, string> &sample_file, map<int, CLocus *> &catalog, QLocus *qloc) {
    vector<SNP *>::iterator i;
    map<string, int>::iterator j;

    int cid = catalog.size();

    CLocus *c = new CLocus;
    c->id = cid + 1;
    c->add_consensus(qloc->con);
    c->sources.push_back(make_pair(sample_file.first, qloc->id));
    catalog[c->id] = c;

    // cerr << "Adding sample: " << qloc->id << " to the catalog as ID: " << c->id << "\n";

    for (i = qloc->snps.begin(); i != qloc->snps.end(); i++) {
	SNP *snp    = new SNP;
	snp->col    = (*i)->col;
	snp->lratio = (*i)->lratio;
	snp->rank_1 = (*i)->rank_1;
	snp->rank_2 = (*i)->rank_2;

	c->snps.push_back(snp);
    }

    for (j = qloc->alleles.begin(); j != qloc->alleles.end(); j++) {
	c->alleles[j->first] = j->second;
    }

    c->populate_alleles();

    return 0;
}

int find_matches_by_sequence(map<int, CLocus *> &catalog, map<int, QLocus *> &sample) {
    //
    // Calculate the distance (number of mismatches) between each pair
    // of Radtags. We expect all radtags to be the same length;
    //
    map<int, QLocus *>::iterator i;
    map<int, CLocus *>::iterator j;
    int k;

    // OpenMP can't parallelize random access iterators, so we convert
    // our map to a vector of integer keys.
    vector<int> keys;
    for (i = sample.begin(); i != sample.end(); i++) 
	keys.push_back(i->first);

    #pragma omp parallel private(i, j, k)
    {
        #pragma omp for schedule(dynamic) 
	for (k = 0; k < (int) keys.size(); k++) {

	    i = sample.find(keys[k]);

	    vector<pair<allele_type, string> >::iterator r, s;

	    //
	    // Iterate through the possible SAMPLE alleles
	    //
	    for (r = i->second->strings.begin(); r != i->second->strings.end(); r++) {

		for (j = catalog.begin(); j != catalog.end(); j++) {
		    //
		    // Iterate through the possible CATALOG alleles
		    //
		    for (s = j->second->strings.begin(); s != j->second->strings.end(); s++) {
			if (r->second == s->second) {
			    //cerr << "Found a match between " << i->first << " (" << r->first << ") and " << j->first << " (" << s->first << ")\n";

			    i->second->add_match(j->second->id, r->first);
			}
		    }
		}
	    }
	}
    }

    return 0;
}

int find_matches_by_genomic_loc(map<int, CLocus *> &catalog, map<int, QLocus *> &sample) {
    //
    // Calculate the distance (number of mismatches) between each pair
    // of Radtags. We expect all radtags to be the same length;
    //
    map<int, QLocus *>::iterator i;
    map<int, CLocus *>::iterator j;
    int k;

    // OpenMP can't parallelize random access iterators, so we convert
    // our map to a vector of integer keys.
    vector<int> keys;
    for (i = sample.begin(); i != sample.end(); i++) 
	keys.push_back(i->first);

    #pragma omp parallel private(i, j, k)
    {
        #pragma omp for schedule(dynamic) 
	for (k = 0; k < (int) keys.size(); k++) {

	    i = sample.find(keys[k]);

            for (j = catalog.begin(); j != catalog.end(); j++) {
                
                if (strcmp(i->second->loc.chr, j->second->loc.chr) == 0 && 
                    i->second->loc.bp == j->second->loc.bp) {

                    i->second->add_match(j->second->id, "");
                }
            }
        }
    }

    return 0;
}

int write_catalog(map<int, CLocus *> &catalog) {
    map<int, CLocus *>::iterator i;
    CLocus  *tag;
    set<int> matches;

    //
    // Parse the input file names to create the output file
    //
    stringstream prefix; 
    prefix << out_path << "batch_" << batch_id;

    //
    // Output the tags
    //
    string out_file = prefix.str() + ".catalog.tags.tsv";
    ofstream cat_file(out_file.c_str());
    out_file = prefix.str() + ".catalog.snps.tsv";
    ofstream snp_file(out_file.c_str());
    out_file = prefix.str() + ".catalog.alleles.tsv";
    ofstream all_file(out_file.c_str());

    for (i = catalog.begin(); i != catalog.end(); i++) {
	tag = i->second;

	write_simple_output(tag, cat_file, snp_file, all_file);

// 	//
// 	// Debugging code.
// 	//
//  	if (tag->matches.size() > 0) {
//  	    cerr << "Found match: " << tag->matches.size() << " matches.\n";
//  	}

//  	for (mat_it = tag->matches.begin(); mat_it != tag->matches.end(); mat_it++) {
//  	    matched_tag = parent_2[(*mat_it).first];
//  	    if (strcmp(matched_tag->con, tag->con) == 0)
//  		cerr << "  Match " << tag->id << " -> " << matched_tag->id << ": consensus matches.\n";
//  	    else 
//  		cerr << "  Match " << tag->id << " -> " << matched_tag->id << ": consensus does not match;\n" 
//  		     << "    " << matched_tag->con << "\n" 
//  		     << "    " << tag->con << "\n";
//  	}

    }

    cat_file.close();
    snp_file.close();
    all_file.close();

    return 0;
}

int CLocus::merge_snps(QLocus *matched_tag) {
    vector<SNP *>::iterator i;
    map<string, int>::iterator j;
    vector<pair<string, SNP *> >::iterator k;
    map<int, pair<string, SNP *> > columns;
    map<int, pair<string, SNP *> >::iterator c;

    vector<pair<string, SNP *> > merged_snps;
    set<string> merged_alleles;
    set<string>::iterator s;

    for (i = this->snps.begin(); i != this->snps.end(); i++)
	columns[(*i)->col] = make_pair("catalog", *i);

    for (i = matched_tag->snps.begin(); i != matched_tag->snps.end(); i++) {
	//
	// Is this column already represented from the previous sample?
	//
	if (columns.count((*i)->col)) {
	    //
	    // Check to make sure we aren't combining inconsistent SNPS.
	    //
	    if (columns[(*i)->col].second->rank_1 != (*i)->rank_1 &&
		columns[(*i)->col].second->rank_1 != (*i)->rank_2)
		cerr << "  Warning: inconsistent merging of SNPs for catalog ID " << this->id << "\n";
	    else
		columns[(*i)->col] = make_pair("both", *i);
	}
	else
	    columns[(*i)->col] = make_pair("sample", *i);
    }

    for (c = columns.begin(); c != columns.end(); c++) 
	merged_snps.push_back((*c).second);

    //
    // Sort the SNPs by column
    //
    sort(merged_snps.begin(), merged_snps.end(), compare_pair);

    //
    // Merge the alleles accounting for any SNPs added from either of the two samples.
    //
    string allele, new_allele;
    int pos;

    for (j = this->alleles.begin(); j != this->alleles.end(); j++) {
	allele     = j->first;
	new_allele = "";
	pos        = 0;

	for (k = merged_snps.begin(); k != merged_snps.end(); k++) {
	    //
	    // If we inserted a SNP from the sample, add the proper nucleotide from the consensus
	    // sequence to account for it in the allele string.
	    //
	    if ((*k).first == "sample") {
		new_allele += this->con[(*k).second->col];
	    } else {
		new_allele += allele[pos];
		pos++;
	    }
	}

	merged_alleles.insert(new_allele);
    }

    for (j = matched_tag->alleles.begin(); j != matched_tag->alleles.end(); j++) {
	allele     = j->first;
	new_allele = "";
	pos        = 0;

	for (k = merged_snps.begin(); k != merged_snps.end(); k++) {
	    if ((*k).first == "catalog") {
		new_allele += matched_tag->con[(*k).second->col];
	    } else {
		new_allele += allele[pos];
		pos++;
	    }
	}

	merged_alleles.insert(new_allele);
    }

    //
    // Update the catalog entry's list of SNPs and alleles
    //
    this->snps.clear();
    for (k = merged_snps.begin(); k != merged_snps.end(); k++) {
	SNP *snp    = new SNP;
	snp->col    = (*k).second->col;
	snp->lratio = (*k).second->lratio;
	snp->rank_1 = (*k).second->rank_1;
	snp->rank_2 = (*k).second->rank_2;

	this->snps.push_back(snp);
    }

    this->alleles.clear();
    for (s = merged_alleles.begin(); s != merged_alleles.end(); s++) {
	this->alleles[*s] = 0;
    }


    return 1;
}

bool compare_pair(pair<string, SNP *> a, pair<string, SNP *> b) {
    return (a.second->col < b.second->col);
}

int write_simple_output(CLocus *tag, ofstream &cat_file, ofstream &snp_file, ofstream &all_file) {
    vector<SNP *>::iterator           snp_it;
    map<string, int>::iterator        all_it;
    vector<pair<int, int> >::iterator src_it;
    string sources;

    for (src_it = tag->sources.begin(); src_it != tag->sources.end(); src_it++) {
	stringstream s; 
	s << (*src_it).first << "_" << (*src_it).second << ",";
	sources += s.str();
    }
    sources = sources.substr(0, sources.length() - 1);

    cat_file << 
	"0"          << "\t" << 
	batch_id     << "\t" <<
	tag->id      << "\t" <<
        tag->loc.chr << "\t" <<
        tag->loc.bp  << "\t" <<
	"consensus"  << "\t" <<
	"0"          << "\t" <<
	sources      << "\t" <<
	tag->con     << "\t" << 
        0            << "\t" <<  // These flags are unused in cstacks, but important in ustacks
        0            << "\t" <<
        0            << "\n";

    //
    // Output the SNPs associated with the catalog tag
    //
    for (snp_it = tag->snps.begin(); snp_it != tag->snps.end(); snp_it++)
	snp_file << "0\t" << 
	    batch_id          << "\t" <<
	    tag->id          << "\t" << 
	    (*snp_it)->col    << "\t" << 
	    (*snp_it)->lratio << "\t" << 
	    (*snp_it)->rank_1 << "\t" << 
	    (*snp_it)->rank_2 << "\n";

    //
    // Output the alleles associated with the two matched tags
    //
    for (all_it = tag->alleles.begin(); all_it != tag->alleles.end(); all_it++)
	all_file << "0\t" << 
	    batch_id  << "\t" <<
	    tag->id  << "\t" << 
	    all_it->first << "\t" <<
            0 << "\t" <<              // These two fields are used in the pstacks output, not in
            0 << "\n";

    return 0;
}

int initialize_catalog(pair<int, string> &sample, map<int, CLocus *> &catalog) {
    map<int, CLocus *> tmp_catalog;

    //
    // Parse the input files.
    //
    if (!load_loci(sample.second, tmp_catalog))
        return 0;

    //
    // Iterate over the catalog entires and renumber them after recording the source of
    // locus.
    //
    map<int, CLocus *>::iterator j;
    int k = 1;
    for (j = tmp_catalog.begin(); j != tmp_catalog.end(); j++) {
	j->second->sources.push_back(make_pair(sample.first, j->second->id));
        j->second->id = k;

        catalog[k] = j->second;

	k++;        
    }

    return 1;
}

int parse_command_line(int argc, char* argv[]) {
    int c, sid;
    string sstr;

    while (1) {
	static struct option long_options[] = {
	    {"help",         no_argument,       NULL, 'h'},
            {"version",      no_argument,       NULL, 'v'},
	    {"mmatches",     no_argument,       NULL, 'm'},
	    {"genomic_loc",  no_argument,       NULL, 'g'},
	    {"batch_id",     required_argument, NULL, 'b'},
	    {"sample",       required_argument, NULL, 's'},
	    {"sample_id",    required_argument, NULL, 'S'},
	    {"outpath",      required_argument, NULL, 'o'},
	    {"num_threads",  required_argument, NULL, 'p'},
	    {0, 0, 0, 0}
	};
	
	// getopt_long stores the option index here.
	int option_index = 0;
     
	c = getopt_long(argc, argv, "hgvmo:s:S:b:p:", long_options, &option_index);
     
	// Detect the end of the options.
	if (c == -1)
	    break;
     
	switch (c) {
	case 'h':
	    help();
	    break;
	case 'b':
	    batch_id = atoi(optarg);
	    break;
	case 'm':
	    mult_matches++;
	    break;
	case 'g':
	    search_type = genomic_loc;
	    break;
	case 's':
	    sstr = optarg;

	    c = getopt_long(argc, argv, "ho:s:S:b:", long_options, &option_index);
	    if (c != 'S') {
		cerr << "You must specify a sample ID after each sample file.\n";
		help();
		abort();
	    }

	    sid = atoi(optarg);
	    samples.push(make_pair(sid, sstr));
	    break;
     	case 'o':
	    out_path = optarg;
	    break;
        case 'v':
            version();
            break;
	case 'p':
	    num_threads = atoi(optarg);
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

    if (samples.size() == 0) {
	cerr << "You must specify at least one sample file.\n";
	help();
    }

    if (out_path.length() == 0) 
	out_path = ".";

    if (out_path.at(out_path.length() - 1) != '/') 
	out_path += "/";

    return 0;
}

void version() {
    std::cerr << "cstacks " << VERSION << "\n\n";

    exit(0);
}

void help() {
    std::cerr << "cstacks " << VERSION << "\n"
              << "cstacks -b batch_id -s sample_file -S id [-s sample_file_2 -S id_2 ...] [-o path] [-p num_threads] [-g] [-h]" << "\n"
              << "  p: enable parallel execution with num_threads threads.\n"
	      << "  b: MySQL ID of this batch." << "\n"
	      << "  s: TSV file from which to load radtags." << "\n"
	      << "  S: MySQL ID of the sample." << "\n"
	      << "  o: output path to write results." << "\n"
	      << "  m: include tags in the catalog that match to more than one entry." << "\n"
              << "  g: base catalog matching on genomic location, not sequence identity." << "\n"
	      << "  h: display this help messsage." << "\n\n";

    exit(0);
}
