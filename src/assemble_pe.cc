#include <string>

#include <getopt.h>

#include "constants.h"
#include "log_utils.h"
#include "MetaPopInfo.h"
#include "locus.h"
#include "BamCLocReader.h"
#include "debruijn.h"

using namespace std;

//
// Argument globs.
//
bool quiet = false;
string bam_path;
string out_dir;
size_t km_length = -1;
size_t min_km_count = 2;
set<int> locus_wl;
bool gfa_out = false;
bool fasta_out = false;

//
// Extra globs.
//
LogAlterator* lg = NULL;

//
// Function declarations.
//
void process_one_locus(const CLocReadSet& loc, Graph& graph);

void parse_command_line(int argc, char* argv[]);
void report_options(ostream& os);

int main(int argc, char** argv) {

    // Parse arguments.
    parse_command_line(argc, argv);

    // Open the log.
    string lg_path = out_dir + "assemble_pe.log";
    if(!quiet)
        cout << "Logging to '" << lg_path << "'." << endl;
    lg = new LogAlterator(lg_path, quiet);
    init_log(lg->l, argc, argv);
    report_options(cout);
    cout << "\n" << flush;

    // Open the BAM file and parse the header.
    MetaPopInfo mpopi;
    BamCLocReader bam_fh = BamCLocReader(bam_path, mpopi);

    // Process every locus
    CLocReadSet loc (mpopi);
    Graph graph (km_length);
    size_t n_loci = 0;
    if (locus_wl.empty()) {
        // No whitelist.
        while (bam_fh.read_one_locus(loc)) {
            process_one_locus(loc, graph);
            ++n_loci;
        }

    } else {
        while (bam_fh.read_one_locus(loc) && !locus_wl.empty()) {
            if (locus_wl.count(loc.id())) {
                process_one_locus(loc, graph);
                ++n_loci;
                locus_wl.erase(loc.id());
            }
        }
    }
    cout << "Processed " << n_loci << " loci.\n";

    cout << "assemble_pe is done.\n";
    return 0;
}

void process_one_locus(const CLocReadSet& loc, Graph& graph) {
    graph.rebuild(loc, min_km_count);
    if (graph.n_simple_paths() == 0)
        // Graph is empty.
        return;

    if (!graph.topo_sort())
        // Not a DAG.
        return;

    vector<const SPath*> best_path = graph.find_best_path();

    if (gfa_out)
        graph.dump_gfa(out_dir + to_string(loc.id()) + ".gfa");

    if (fasta_out) {
        ofstream fasta (out_dir + to_string(loc.id()) + ".fa");
        for (const Read& r : loc.reads())
            fasta << ">" << r.name << "\n" << r.seq.str() << "\n";
    }
}

const string help_string = string() +
        "assemble_pe " + VERSION  + "\n"
        "assemble_pe -b bam_file -O out_dir\n"
        "\n"
        "  -b: path to the input sorted read, in BAM format\n"
        "  -O: path to an output directory\n"
        "  -k: kmer length\n"
        "  --min-cov: minimum coverage to consider a kmer\n"
        "  -W,--whitelist: a whitelist of locus IDs\n"
        "  --gfa: output a GFA file for each locus\n"
        "  --fasta: output a Fasta file for each locus\n"
        "\n"
        ;

void bad_args() {
    cerr << help_string;
    exit(13);
}

void parse_command_line(int argc, char* argv[]) {

    static const option long_options[] = {
        {"help",         no_argument,       NULL, 'h'},
        {"quiet",        no_argument,       NULL, 'q'},
        {"version",      no_argument,       NULL,  1000},
        {"bam-file",     required_argument, NULL, 'b'},
        {"out-dir",      required_argument, NULL, 'O'},
        {"kmer-length",  required_argument, NULL, 'k'},
        {"min-cov",      required_argument, NULL,  1003},
        {"whitelist",    required_argument, NULL,  'W'},
        {"gfa",          no_argument,       NULL,  1004},
        {"fasta",        no_argument,       NULL,  1002},
        {0, 0, 0, 0}
    };

    string wl_path;

    int c;
    int long_options_i;
    while (true) {

        c = getopt_long(argc, argv, "hqb:O:k:W:", long_options, &long_options_i);

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
        case 'b':
            bam_path = optarg;
            break;
        case 'O':
            out_dir = optarg;
            break;
        case 'k':
            km_length = atoi(optarg);
            break;
        case 1003://min-cov
            min_km_count = atoi(optarg);
            break;
        case 'W':
            wl_path = optarg;
            break;
        case 1004://gfa
            gfa_out = true;
            break;
        case 1002://fasta
            fasta_out = true;
            break;
        case '?':
            bad_args();
            break;
        }
    }

    // Check command consistency.
    if (optind < argc) {
        cerr << "Error: Failed to parse command line: '" << argv[optind] << "' is seen as a positional argument. Expected no positional arguments.\n";
        bad_args();
    }
    if (bam_path.empty()) {
        cerr << "Error: An input BAM file must be provided (-b).\n";
        bad_args();
    }
    if (out_dir.empty()) {
        cerr << "Error: An output directory must be provided (-O).\n";
        bad_args();
    }
    if (km_length == size_t(-1)) {
        cerr << "Error: A kmer length must be provided (-k).\n";
        bad_args();
    }

    // Process arguments.
    if (out_dir.back() != '/')
        out_dir += '/';

    if (!wl_path.empty()) {
        ifstream wl_fh (wl_path);
        if (!wl_fh) {
            cerr << "Error: Failed to open '" << wl_path << "' for reading.\n";
            throw exception();
        }
        int id;
        while (wl_fh >> id)
            locus_wl.insert(id);
        if (locus_wl.empty()) {
            cerr << "Error: Whitelist '" << wl_path << "' appears empty.\n";
            throw exception();
        }
    }
}

void report_options(ostream& os) {
    os << "Configuration for this run:\n"
       << "  Sorted reads file: '" << bam_path << "'\n"
       << "  Output directory: '" << out_dir << "'\n"
       << "  Kmer length: " << km_length << "\n"
       << "  Min coverage: " << min_km_count << "\n"
       ;

    if (!locus_wl.empty())
        os << "  Whitelist of " << locus_wl.size() << " loci.\n";
}
