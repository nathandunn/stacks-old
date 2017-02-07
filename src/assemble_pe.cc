#include <string>

#include <getopt.h>

#include "constants.h"
#include "log_utils.h"
#include "MetaPopInfo.h"
#include "locus.h"
#include "BamI.h"

using namespace std;

//
// Argument globs.
//
bool quiet = false;
string bam_path;
string out_dir;

//
// Extra globs.
//
LogAlterator* lg = NULL;

//
// Function declarations.
//
bool read_one_locus(CLocReadSet& loc, Bam* bam_f, const map<string, size_t>& rg_to_sample);
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
    cout << "\n";

    // Open the BAM file.
    Bam* bam_f = new Bam(bam_path.c_str());

    // Get the list of samples from the BAM header.
    MetaPopInfo mpopi;
    map<string, size_t> rg_to_sample; // (Read group ID, sample index)
    {
        BamHeader::ReadGroups read_groups = bam_f->h().read_groups();

        // Create the MetaPopInfo object
        vector<string> samples;
        for (auto rg : read_groups)
            samples.push_back(rg.second.at("SM"));
        mpopi.init_names(samples);

        // Get the (read group : sample) map
        for (auto rg : read_groups)
            rg_to_sample.insert({rg.first, mpopi.get_sample_index(rg.second.at("SM"))});
    }

    // Process every locus
    CLocReadSet loc (mpopi);
    bool eof;
    if (!bam_f->next_record()) {
        cerr << "Error: Failed to read records from BAM file '" << bam_path << "'.\n";
        throw exception();
    }
    do {
        eof = !read_one_locus(loc, bam_f, rg_to_sample);

        cout << "Locus #" << loc.id() << "\n";
        for (const Read& r : loc.reads()) {
            cout << r.name << " (" << mpopi.samples()[loc.sample_of(r)].name << ") ";
            for (auto nt = r.s.begin(); nt != r.s.end(); ++nt)
                cout << nt4::u2c[*nt];
            cout << "\n";
        }

        // TODO Process the locus.

    } while (!eof);

    return 0;
}

bool read_one_locus(CLocReadSet& loc, Bam* bam_f, const map<string, size_t>& rg_to_sample) {
    loc.clear();
    const BamRecord& rec = bam_f->r(); // the current record

    // Parse the locus ID.
    int32_t curr_chrom = rec.chrom();
    loc.id(stoi(bam_f->h().chrom_str(curr_chrom)));

    // Read all the reads of the locus, and one more.
    do {
        loc.add(Read(rec.seq(), rec.qname()), rg_to_sample.at(rec.read_group()));
        if(!bam_f->next_record())
            // EOF.
            return false;
    } while (rec.chrom() == curr_chrom);

    return true;
}

const string help_string = string() +
        "assemble_pe " + VERSION  + "\n"
        "assemble_pe -b bam_file -O out_dir\n"
        "\n"
        "  -b: path to the input sorted read, in BAM format\n"
        "  -b: path to an output directory\n"
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
        {"bam_file",     required_argument, NULL, 'b'},
        {"out_dir",      required_argument, NULL, 'O'},
        {0, 0, 0, 0}
    };

    int c;
    int long_options_i;
    while (true) {

        c = getopt_long(argc, argv, "hqb:O:", long_options, &long_options_i);

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

    // Process arguments.
    if (out_dir.back() != '/')
        out_dir += '/';
}

void report_options(ostream& os) {
    os << "Configuration for this run:\n"
       << "  Sorted reads file: '" << bam_path << "'\n"
       << "  Output directory: '" << out_dir << "'\n"
       ;
}
