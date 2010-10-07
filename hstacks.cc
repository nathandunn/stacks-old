// -*-mode:c++; c-style:k&r; c-basic-offset:4;-*-
//
// hstacks -- find homologous stacks among a set of samples
//
// Match stacks between samples to identify homologous loci. Stacks
// may contain masked sites (N's), resulting from using the fixed-model in
// the ustacks program, or an explicit number of mismatches per tag may be
// specified, allowing non-exact matching homologus sites to be identified
// across a set of samples.
//
// Julian Catchen
// jcatchen@uoregon.edu
// University of Oregon
//

#include "hstacks.h"

// Global variables to hold command-line options.
string in_path;
string out_path;
int    batch_id        = 0;
int    num_threads     = 1;
int    stack_depth_min = 1;

int main (int argc, char* argv[]) {

    parse_command_line(argc, argv);

    //
    // Set the number of OpenMP parallel threads to execute.
    //
    omp_set_num_threads(num_threads);

    vector<string>           input_files;
    vector<string>::iterator in_file;
    map<int, HLocus *>       samples;

    build_file_list(in_path, input_files);

    int id = 0;

    for (in_file = input_files.begin(); in_file != input_files.end(); in_file++) {
        map<int, HLocus *> sample;
        map<int, HLocus *>::iterator it;

	size_t pos_1     = (*in_file).find_last_of("/");
	size_t pos_2     = (*in_file).find_last_of(".");
	string sample_id = (*in_file).substr(pos_1 + 1, (pos_2 - pos_1 - 1));

	load_loci(*in_file, sample);

        //
        // Give each locus a unique ID among all samples
        //
        for (it = sample.begin(); it != sample.end(); it++) {
            it->second->uniq_id = id;
            samples[id]         = (*it).second;
            id++;
        }
    }

    //
    // Calculate distance between tags in different samples.
    //
    calc_distance(samples, 1);

    //
    // Write out tags matched between samples.
    //
    write_homologous_loci(samples);

    return 0;
}

int calc_distance(map<int, HLocus *> &loci, int utag_dist) {
    //
    // Calculate the distance (number of mismatches) between each pair
    // of Radtags. We expect all radtags to be the same length;
    //
    map<int, HLocus *>::iterator it;
    HLocus *tag_1, *tag_2;
    int i, j;

    cerr << "Calculating distance between stacks...\n";

    // OpenMP can't parallelize random access iterators, so we convert
    // our map to a vector of integer keys.
    vector<int> keys;
    for (it = loci.begin(); it != loci.end(); it++) 
        keys.push_back(it->first);

    #pragma omp parallel private(i, j, tag_1, tag_2)
    { 
        #pragma omp for schedule(dynamic) 
	for (i = 0; i < (int) keys.size(); i++) {

	    tag_1 = loci[keys[i]];

	    int d;

	    for (j = 0; j < (int) keys.size(); j++) {
		tag_2 = loci[keys[j]];

		// Don't compare tag_1 against itself.
		if (tag_1 == tag_2)
		    continue;

		d = dist(tag_1, tag_2);

		//
		// Store the distance between these two sequences if it is
		// below the maximum distance.
		//
		if (d == utag_dist) {
                    if (tag_1->depth < stack_depth_min ||
                        tag_2->depth < stack_depth_min)
                        continue;

		    tag_1->add_match(tag_2->uniq_id, d);
		}
	    }

	    // Sort the vector of distances.
	    sort(tag_1->matches.begin(), tag_1->matches.end(), compare_dist);
	}
    }

    return 0;
}

int dist(HLocus *tag_1, HLocus *tag_2) {
    int   dist = 0;
    char *p    = tag_1->con;
    char *q    = tag_2->con;
    char *end  = p + strlen(p);

    // Count the number of characters that are different
    // between the two sequences. Don't count wildcard 'N'
    // nucleotides.
    while (p < end) {
	dist += ((*p == *q) || (*q == 'N' || *p == 'N')) ? 0 : 1;
	p++;
	q++;
    }

    return dist;
}

bool compare_dist(Match *a, Match *b) {
    return (a->dist < b->dist);
}

int write_homologous_loci(map<int, HLocus *> &samples) {
    map<int, HLocus *>::iterator i;
    vector<int>::iterator k;
    set<int>   write_map;
    HLocus *tag_1, *tag_2;

    cerr << "Writing homologous stacks...\n";

    //
    // Parse the input file name to create the output files
    //
    stringstream prefix; 
    prefix << out_path << "batch_" << batch_id;

    string tag_file = prefix.str() + ".homologous.tags.tsv";
    string nuc_file = prefix.str() + ".homologous.nucs.tsv";

    // Open the output files for writing.
    std::ofstream tags(tag_file.c_str());
    std::ofstream nucs(nuc_file.c_str());
    int id = 0;

    for (i = samples.begin(); i != samples.end(); i++) {
	tag_1 = i->second;

	//
	// This tag may already have been merged by an earlier operation.
	//
	if (write_map.find(tag_1->uniq_id) != write_map.end())
	    continue;

	set<int> unique_merge_list;
	set<int>::iterator it;

	trace_stack_graph(tag_1, samples, unique_merge_list);

	for (it = unique_merge_list.begin(); it != unique_merge_list.end(); it++) {
	    tag_2 = samples[(*it)];

	    // Record the nodes that have been merged in this round.
	    write_map.insert(tag_2->uniq_id);

            //
            // For each tag we are outputting, output the depth of coverage for each
            // nucleotide in that stack separately (in order to calculate correlations
            // between depth of coverage and fixed/non-fixed nucleotides.
            //
            if (unique_merge_list.size() > 1) {
                char *p, *end;
                end = tag_2->con + strlen(tag_2->con);
                for (p = tag_2->con; p < end; p++)
                    nucs 
                        << tag_2->sample_id << "_" << tag_2->id << "\t" 
                        << *p << "\t" 
                        << tag_2->depth << "\n";
            }

            //
            // Output the consensus tag for this population.
            //
	    tags 
		<< "0" << "\t" 
		<< batch_id << "\t" 
		<< id << "\t" 
		<< tag_2->uniq_id << "\t" 
		<< tag_2->sample_id << "_" <<  tag_2->id << "\t" 
		<< tag_2->con << "\n";
	}

	id++;
    }

    tags.close();
    nucs.close();

    cerr << "  Wrote " << id - 1 << " tags to " << tag_file << "\n";

    return 0;
}

int trace_stack_graph(HLocus *tag_1, map<int, HLocus *> &loci, set<int> &unique_merge_list) {
    queue<int>                        merge_list;
    pair<set<int>::iterator,bool>     ret;
    vector<Match *>::iterator k;
    HLocus *tag_2;

    unique_merge_list.insert(tag_1->uniq_id);
    merge_list.push(tag_1->uniq_id);

    while (!merge_list.empty()) {
	tag_2 = loci[merge_list.front()];
	merge_list.pop();

	for (k = tag_2->matches.begin(); k != tag_2->matches.end(); k++) {
	    ret = unique_merge_list.insert((*k)->id);

	    //
	    // If this Tag has not already been added to the merge list (i.e. we were able
	    // to insert it in to our unique_merge_list, which is a set), add it for consideration
	    // later in the loop.
	    //
	    if (ret.second == true)
		merge_list.push((*k)->id);
	}
    }

    return 0;
}

int build_file_list(string in_path, vector<string> &sql_files) {
    struct dirent *dirp;
    string d;
    size_t found;

    DIR *dp = opendir(in_path.c_str());

    if (dp == NULL) {
        cerr << "Error (" << errno << ") opening " << in_path << "\n";
        return errno;
    }

    while ((dirp = readdir(dp)) != NULL) {
	d     = string(dirp->d_name);
	found = d.find("tags.tsv");

	if (found != string::npos) {
            size_t pos = d.find(".tags.tsv");            
            d = in_path + d.substr(0, pos); 
	    sql_files.push_back(d);
	}
    }

    closedir(dp);

    cerr << "Identified " << sql_files.size() << " samples.\n";

    return 0;
}

int parse_command_line(int argc, char* argv[]) {
    int c;
     
    while (1) {
	static struct option long_options[] = {
	    {"help",        no_argument,       NULL, 'h'},
            {"version",     no_argument,       NULL, 'v'},
	    {"inpath",      required_argument, NULL, 'p'},
	    {"outpath",     required_argument, NULL, 'o'},
	    {"batch_id",    required_argument, NULL, 'b'},
	    {0, 0, 0, 0}
	};
	
	// getopt_long stores the option index here.
	int option_index = 0;
     
	c = getopt_long(argc, argv, "hvi:p:o:b:e:m:", long_options, &option_index);
     
	// Detect the end of the options.
	if (c == -1)
	    break;
     
	switch (c) {
	case 'h':
	    help();
	    break;
     	case 'v':
	    version();
	    break;
	case 'i':
	    in_path = optarg;
	    break;
	case 'o':
	    out_path = optarg;
	    break;
	case 'b':
	    batch_id = atoi(optarg);
	    break;
	case 'm':
	    stack_depth_min = atoi(optarg);
	    break;
	case 'p':
	    num_threads = atoi(optarg);
	    break;
	case '?':
	    // getopt_long already printed an error message.
	    help();
	    break;
	default:
	    cerr << "Unknown command line option '" << (char) c << "'\n";
	    help();
	    abort();
	}
    }

    if (in_path.length() == 0) {
	cerr << "You must specify a path to a set of input files.\n";
	help();
    }

    if (in_path.at(in_path.length() - 1) != '/') 
	in_path += "/";

    if (out_path.length() == 0) 
	out_path = ".";

    if (out_path.at(out_path.length() - 1) != '/') 
	out_path += "/";

    return 0;
}

void version() {
    std::cerr << "hstacks " << stacks_version << "\n\n";

    exit(0);
}

void help() {
    std::cerr << "hstacks " << stacks_version << "\n"
              << "hstacks -i path [-o path] [-b batch_id] [-m min] [-p min_threads] [-h]" << "\n"
	      << "  i: path to the set of SQL files from which to load loci." << "\n"
	      << "  o: output path to write results." << "\n"
	      << "  b: SQL Batch ID to insert into the output to identify a group of samples." << "\n"
              << "  m: minimum stack depth required for a locus to be included in the search." << "\n"
              << "  p: enable parallel execution with num_threads threads.\n"
	      << "  h: display this help messsage." << "\n\n";

    exit(0);
}
