#include <cmath>
#include <cstring>
#include <iostream>
#include <fstream>

#include <getopt.h>

#include "constants.h"
#include "sql_utilities.h"
#include "log_utils.h"
#include "BamI.h"
#include "DNASeq4.h"

void parse_command_line(int argc, char* argv[]);
void report_options(ostream& os);
void run();

//
// Argument globals.
//
bool quiet = false;
string prefix_path;
string pe_reads_path;

//
// Extra globals.
//
const string prog_name = "tsv2bam";
LogAlterator* lg = NULL;

// main()
// ==========
int main(int argc, char* argv[]) {
    IF_NDEBUG_TRY

    // Parse arguments
    parse_command_line(argc, argv);

    // Open the log
    string lg_path = prefix_path + "." + prog_name + ".log";
    if(!quiet)
        cout << "Logging to '" << lg_path << "'." << endl;
    lg = new LogAlterator(lg_path, quiet);
    init_log(lg->l, argc, argv);
    report_options(cout);
    cout << "\n";

    // Run the program.
    run();

    // Cleanup.
    cout << prog_name << " is done.\n";
    delete lg;
    return 0;
    IF_NDEBUG_CATCH_ALL_EXCEPTIONS
}

void run() {
}

void parse_command_line(int argc, char* argv[]) {

    const string help_string = string() +
            prog_name + " " + VERSION  + "\n" +
            prog_name + " -s sample_prefix [-f paired_reads]\n"
            "\n"
            "  -s,--prefix: prefix path for the sample.\n"
            "  -f,--pe-reads: path to the sample's paired-end read sequences (if any).\n"
            "\n"
            ;

    const auto bad_args = [&help_string](){
        cerr << help_string;
        exit(13);
    };

    static const option long_options[] = {
        {"help",         no_argument,       NULL, 'h'},
        {"quiet",        no_argument,       NULL, 'q'},
        {"version",      no_argument,       NULL,  1000},
        {"prefix",       required_argument, NULL, 's'},
        {"pe-reads",     required_argument, NULL, 'f'},
        {0, 0, 0, 0}
    };

    int c;
    int long_options_i;
    while (true) {

        c = getopt_long(argc, argv, "hqs:f:", long_options, &long_options_i);

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
            cout << prog_name << " " << VERSION << "\n";
            exit(0);
            break;
        case 's':
            prefix_path = optarg;
            break;
        case 'f':
            pe_reads_path = optarg;
            break;
        default:
            bad_args();
            break;
        }
    }

    //
    // Check command consistency.
    //
    if (optind < argc) {
        cerr << "Error: Failed to parse command line: '" << argv[optind] << "' is seen as a positional argument. Expected no positional arguments.\n";
        bad_args();
    }

    if (prefix_path.empty()) {
        cerr << "Error: A sample prefix path is required (-s).\n";
        bad_args();
    }

    if (!pe_reads_path.empty() && guess_file_type(pe_reads_path) == FileT::unknown) {
        cerr << "Error: Failed to recognize the format of '" << pe_reads_path << "'.\n";
        bad_args();
    }
}

void report_options(ostream& os) {
    os << "Configuration for this run:\n";
    os << "  Sample prefix: '" << prefix_path << "'\n";
    if (!pe_reads_path.empty())
        os << "  Paired-end reads path: '" << pe_reads_path << "'\n";
}
