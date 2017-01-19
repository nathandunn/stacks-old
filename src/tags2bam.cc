#include <iostream>
#include <fstream>

#include <getopt.h>

#include "constants.h"
#include "sql_utilities.h"
#include "log_utils.h"

#include "tags2bam.h"

using namespace std;

//
// Argument globs.
//
bool quiet = false;
string prefix_path;

//
// Extra globs.
//
LogAlterator* lg = NULL;

// main()
// ==========
int main(int argc, char* argv[]) {

    // Parse arguments
    parse_command_line(argc, argv);

    // Open the log
    string lg_path = prefix_path + ".t2bam.log";
    if(!quiet)
        cout << "Logging to '" << lg_path << "'." << endl;
    lg = new LogAlterator(lg_path, quiet);
    init_log(lg->l, argc, argv);
    report_options(cout);
    cout << "\n";

    // Read the sample's files.
    map<int, Locus*> sloci;
    unordered_map<int, int> sloc_to_cloc;
    read_sample_files(sloci, sloc_to_cloc);

    // Write the BAM file.
    write_bam_file(sloci, sloc_to_cloc);

    cout << "tags2bam is done." << endl;
    return 0;
}

void read_sample_files(map<int, Locus*>& sloci, unordered_map<int, int>& sloc_to_cloc) {
    // Read the matches file.
    cout << "Reading the matches file..." << endl;
    vector<CatMatch*> matches;
    load_catalog_matches(prefix_path, matches);
    if (matches.empty()) {
        cerr << "Error: Unable to load matches from '"
             << prefix_path + ".matches.tsv(.gz)'.\n";
        throw exception();
    }
    int sample_id = matches[0]->sample_id;

    // Retrieve the list of bijective loci
    vector<pair<int, int> > bij_loci = retrieve_bijective_loci(matches);
    sloc_to_cloc.reserve(bij_loci.size());
    sloc_to_cloc.insert(bij_loci.begin(), bij_loci.end());

    // Read the tags file
    cout << "Reading the tags file..." << endl;
    bool tmpbool;
    if(load_loci(prefix_path, sloci, 2, false, tmpbool) != 1) {
        cerr << "Error: could not find stacks files '" << prefix_path << ".*' (tags, snps and/or alleles).\n";
        throw exception();
    }

    // Discard loci that aren't bijective with the catalog
    for (auto sloc=sloci.begin(); sloc!=sloci.end();) {
        if (sloc_to_cloc.count(sloc->second->id)) {
            sloc++;
        } else {
            delete sloc->second;
            sloci.erase(sloc++);
        }
    }

    for (CatMatch* m : matches)
        delete m;
}

void write_bam_file(const std::map<int, Locus*>& sloci, const std::unordered_map<int, int>& sloc_to_cloc) {
    //TODO
}

const string help_string = string() +
        "tags2bam " + VERSION  + "\n"
        "tags2bam -s sample_prefix\n"
        "\n"
        "  -s,--prefix: prefix path for the sample.\n"
        "\n"
        ;

void bad_args() {
    cerr << help_string;
    exit(-1);
}

void parse_command_line(int argc, char* argv[]) {

    static const option long_options[] = {
        {"help",         no_argument,       NULL, 'h'},
        {"quiet",        no_argument,       NULL, 'q'},
        {"version",      no_argument,       NULL,  1000},
        {"prefix",       required_argument, NULL, 's'},
        {0, 0, 0, 0}
    };

    int c;
    int long_options_i;
    while (true) {

        c = getopt_long(argc, argv, "hqs:", long_options, &long_options_i);

        // Detect the end of the options.
        if (c == -1)
            break;

        switch (c) {
        case 'h':
            cout << help_string;
            exit(0);
            break;
        case 'q':
            quiet = true;
            break;
        case 1000: //version
            cout << "pstacks_pe " << VERSION << "\n";
            exit(0);
            break;
        case 's':
            prefix_path = optarg;
            break;
        case '?':
            bad_args();
            break;
        }
    }

    //
    // Check command consistency.
    //

    if (prefix_path.empty()) {
        cerr << "Error: A sample prefix path is required (-s).\n";
        bad_args();
    }
}

void report_options(ostream& os) {
    os << "Configuration for this run:\n";
    os << "  Sample prefix: '" << prefix_path << "'\n";
}
