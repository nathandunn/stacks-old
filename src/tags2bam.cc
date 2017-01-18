#include <iostream>
#include <fstream>

#include <getopt.h>

#include "constants.h"
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

    // TODO

    cout << "tags2bam is done." << endl;
    return 0;
}

const string help_string = string() +
        "pstacks_pe " + VERSION  + "\n"
        "pstacks_pe -s sample_prefix -f paired_reads_file [-t type] [-p n_threads] [--model type [...]]\n"
        "\n"
        "  -s,--prefix: prefix path for the sample.\n"
        "  -f,--paired_reads: path to the paired reads alignments.\n"
        "  -t,--filetype: alignment file type. Supported types: bowtie, sam, or bam. (Default: guess)\n"
        "  -p,--num_threads: number of threads to use for parallel execution. (Default: 1)\n"
        "\n"
        "Model:\n"
        "  --model: Model for calling SNPs. One of 'snp', 'bounded', or 'fixed'. (Default: snp)\n"
        "  For the SNP and Bounded SNP models:\n"
        "    --alpha: required confidence levels for SNP calls. One of 0.1, 0.05 (default), 0.01, or 0.001.\n"
        "  For the Bounded SNP model:\n"
        "    --bound_low <float>: lower bound for epsilon, the error rate. Between 0 and 1 (default 0).\n"
        "    --bound_high <num>: upper bound for epsilon, the error rate, between 0 and 1.0 (default 1).\n"
        "  For the Fixed model:\n"
        "    --bc_err_freq <num>: specify the barcode error frequency, between 0 and 1.0.\n"
        "\n"
        "Misc.: -q,--quiet; -h,--help; --version.\n"
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
        cerr << "You must specify the prefix path of the sample (-s).\n";
        bad_args();
    }
}

void report_options(ostream& os) {
    os << "Configuration for this run:\n";
    os << "  Sample prefix: '" << prefix_path << "'\n";
}
