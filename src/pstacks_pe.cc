#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <list>
#include <map>
#include <unordered_set>
#include <unordered_map>

#include <getopt.h>
#ifdef _OPENMP
#include <omp.h>
#endif

#include "constants.h"
#include "log_utils.h"
#include "stacks.h"
#include "models.h"
#include "input.h"
#include "BamI.h"
#include "SamI.h"
#include "BowtieI.h"
#include "Tsv.h"
#include "sql_utilities.h"

#include "pstacks_base.h"
#include "pstacks_pe.h"

#ifdef DEBUG
#define IF_NDEBUG_TRY
#define IF_NDEBUG_CATCH_ALL_EXCEPTIONS
#else
#define IF_NDEBUG_TRY \
    try {
#define IF_NDEBUG_CATCH_ALL_EXCEPTIONS \
    } catch (const std::exception& e) { \
        std::cerr << "Aborted. (" << e.what() << ").\n"; \
        return -1; \
    }
#endif

using namespace std;

//
// Argument globs.
//
string prefix_path;
string paired_alns_path;
FileT  in_file_type = FileT::unknown;
int    num_threads        = 1;
bool quiet = false;

modelt model_type         = snp;
double alpha              = 0.05;
double bound_low          = 0.0;
double bound_high         = 1.0;
double barcode_err_freq   = 0.0;

//
// Extra globs.
//
LogAlterator* lg = NULL;
int    sql_id             = -1;
double heterozygote_limit;
double homozygote_limit;

const int barcode_size    = 5;
double p_freq             = 0.5; // const

set<string> debug_flags;
#define DEBUG_FWREADS "FWREADS"
#define DEBUG_READNAMES "READNAMES"
const set<string> known_debug_flags {DEBUG_FWREADS, DEBUG_READNAMES};

// main()
// ==========
int main(int argc, char* argv[]) {
    IF_NDEBUG_TRY

    // Parse arguments
    parse_command_line(argc, argv);

    // Open the log
    string lg_path = prefix_path + ".pstacks_pe.log";
    if(!quiet)
        cout << "Logging to '" << lg_path << "'.\n";
    lg = new LogAlterator(lg_path, quiet);
    init_log(lg->l, argc, argv);
    report_options(cout);
    cout << endl;

    // Initialize stuff.
    set_model_thresholds(alpha);
#ifdef _OPENMP
    omp_set_num_threads(num_threads);
#endif

    // Retrieve bijective matches to the catalog.
    // ----------

    // Parse the matches file.
    cout << "Loading matches to the catalog..." << endl;
    vector<CatMatch*> matches;
    load_catalog_matches(prefix_path, matches);
    if (matches.empty()) {
        cerr << "Error: Unable to load matches from '"
             << prefix_path + ".matches.tsv(.gz)'.\n";
        throw exception();
    }
    sql_id = matches[0]->sample_id;

    unordered_set<int> bij_sloci = retrieve_bijective_sloci(matches);

    for (CatMatch* m : matches)
        delete m;
    matches.clear();

     // Parse the tags files; sort the reads per locus.
     // ----------
    cout << "Reading read-to-locus information from the tags file..." << endl;
    unordered_map<string, size_t> read_name_to_loc; // map of (read name, loc index)
    vector<FwLocInfo> fwloci_info;
    bool is_input_gzipped;
    link_reads_to_loci(bij_sloci, read_name_to_loc, fwloci_info, is_input_gzipped);
    const size_t n_loci = fwloci_info.size();
    cout << "This sample covers " << n_loci
         << " catalog loci with " << read_name_to_loc.size() << " reads." << endl;

    // Load the paired-ends into PStacks.
    // ----------
    cout << "Loading the paired-end sequences..." << endl;
    vector<vector<PStack> > stacks_per_loc = load_aligned_reads(fwloci_info, read_name_to_loc);
    read_name_to_loc.clear();

    // Merge PStacks into MergedStacks
    // ----------
    // We check that there are paired-end reads for each locus.
    //
    // This is not necessarily the case if forward and paired reads were
    // aligned independently : the forward read may have been kept
    // while the paired-end read was discarded. Loci without any paired-end
    // reads just do not have entries in the `tags_pe.tsv` file.
    cout << "Merging stacks..." << endl;

    vector<MergedStack> pe_loci;
    pe_loci.reserve(n_loci);
    for (size_t i=0; i<n_loci; ++i) {
        vector<PStack>& stacks = stacks_per_loc[i];
        pe_loci.push_back(stacks.empty() ? MergedStack() : merge_pstacks(stacks, fwloci_info[i].id));
    }

    // Report results.
    size_t n_pe_loci = 0;
    size_t n_stacks = 0;
    for (const MergedStack& l : pe_loci) {
        if (l.utags.empty())
            continue;
        ++n_pe_loci;
        n_stacks+=l.utags.size();
    }
    cout << "Created " << n_pe_loci << " paired loci with " << n_stacks << " stacks." << endl;

    // Call SNPs and alleles.
    // ----------
    cout << "Calling SNPs..." << endl;

    // Create the maps to pass to `call_consensus` and `write_results`.
    map<int, MergedStack*> loci_map;
    map<int, PStack*> stacks_map;
    for (size_t i=0; i<n_loci; ++i) {
        MergedStack& loc = pe_loci[i];
        if(loc.utags.empty())
            continue;
        loci_map.insert({loc.id, &loc});

        for (PStack& stack : stacks_per_loc[i])
            stacks_map.insert({stack.id, &stack});
    }

    // Call the variants.
    call_consensus(loci_map, stacks_map, true);

    // Write results.
    // ----------
    cout << "Writing results..." << endl;
    write_results(loci_map, stacks_map, is_input_gzipped, true);

    cout << "pstacks_pe is done." << endl;
    return 0;

    IF_NDEBUG_CATCH_ALL_EXCEPTIONS
}

void convert_fw_read_name_to_paired(string& read_name) {

    // Check the format.
    if (read_name.length() < 2
            || (read_name.substr(read_name.length()-2) != "/1"
                    && read_name.substr(read_name.length()-2) != "_1")
            ){
        cerr << "Error: Unrecognized read name format; expected '"
             << read_name << "' to end with '/1' or '_1'.\n";
        throw exception();
    }

    // Change the 1 into a 2.
    read_name.back() = '2';
}

void link_reads_to_loci(
        const unordered_set<int>& bij_sloci,
        unordered_map<string, size_t>& pread_name_to_loc,
        vector<FwLocInfo>& sloci_info,
        bool& is_input_gzipped
        ) {

    // Parse the sloci in the tags file; guess the names of the associated
    // paired-end reads.
    map<int, Locus*> sloci;
    if(load_loci(prefix_path, sloci, 1, false, is_input_gzipped) != 1) {
        cerr << "Error: could not find stacks files '" << prefix_path << ".*' (tags, snps and/or alleles).\n";
        throw exception();
    }

    const bool process_names = debug_flags.count(DEBUG_FWREADS) ? false : true;
    for (const auto& element : sloci) {
        const Locus& sloc = *element.second;
        if (not bij_sloci.count(sloc.id))
            continue;

        sloci_info.push_back(FwLocInfo(sloc.id, sloc.loc));

        // For each first read,
        for (const char* fread_name : sloc.comp) {
            string pread_name (fread_name);
            if (process_names)
                convert_fw_read_name_to_paired(pread_name);
            pread_name_to_loc.insert( {pread_name, sloci_info.size()-1} );
        }
    }

    for (const auto& element : sloci)
        delete element.second;

    // If required, save the list of reads.
    if (debug_flags.count(DEBUG_READNAMES)) {
        ofstream readnames_f ("readnames");
        for (auto& r : pread_name_to_loc)
            readnames_f << r.first << "\n";
    }
}

vector<vector<PStack> > load_aligned_reads(
        const vector<FwLocInfo>& fwloci_info,
        const std::unordered_map<std::string, size_t>& read_name_to_loc
        ) {
    vector<vector<PStack> > stacks_per_loc;

    // Open the alignment file.
    Input* pe_reads_f;
    if (in_file_type == FileT::bam)
        pe_reads_f = new Bam(paired_alns_path.c_str());
    else if (in_file_type == FileT::sam)
        pe_reads_f = new Sam(paired_alns_path.c_str());
    else if (in_file_type == FileT::bowtie)
        pe_reads_f = new Bowtie(paired_alns_path.c_str());
    else if (in_file_type == FileT::tsv)
        pe_reads_f = new Tsv(paired_alns_path.c_str());
    else
        assert(NEVER);

    // Read and stack the alignments, per locus.
    vector<set<PStack> > stack_sets_per_loc (fwloci_info.size());
    size_t n_used_reads = 0;
    size_t next_stack_id = 0;
    Seq seq;
    while(pe_reads_f->next_seq(seq)) {
        auto it = read_name_to_loc.find(seq.id);
        if (it == read_name_to_loc.end())
            continue;
        if (! fwloci_info[it->second].is_upstream_of(seq))
            continue;

        ++n_used_reads;
        set<PStack>& loc = stack_sets_per_loc[it->second];

        PStack key;
        key.loc = seq.loc;
        key.add_seq(seq.seq);

        auto ins = loc.insert(move(key));
        set<PStack>::iterator stack = ins.first;
        if (ins.second) {
            PStack::set_id_of(stack, next_stack_id);
            ++next_stack_id;
        }
        PStack::add_read_to(stack, seq.id);
    }
    if (n_used_reads == 0) {
        cerr << "Error: Failed to find any matching paired-end reads in '" << paired_alns_path << "'." << endl;
        throw exception();
    }

    // Get the PStacks out of their set.
    stacks_per_loc.reserve(fwloci_info.size());
    for (set<PStack>& loc : stack_sets_per_loc) {
        stacks_per_loc.push_back(vector<PStack>(loc.begin(), loc.end())); // But see set::extract in c++17
        loc.clear();
    }

    cout << "Found " << n_used_reads << " aligned paired-end reads, created " << next_stack_id << " stacks." << endl;
    return stacks_per_loc;
}

MergedStack merge_pstacks(vector<PStack>& pstacks, int loc_id) {
    MergedStack locus;
    assert(!pstacks.empty());
    locus.id = loc_id;

    // Determine the positions spanned by the PStacks.
    // ----------
    list<Contig> contigs;
    for (const PStack& p : pstacks) {
        auto c = contigs.begin();
        for (; c!=contigs.end(); ++c)
            if (c->add(p))
                break;
        if (c == contigs.end())
            // This PStack doesn't overlap any of the existing contigs.
            contigs.push_back(Contig(p));
    }

    // Some contigs may have secondarily become overlapping, merge them.
    for (auto c=contigs.begin(); c!=contigs.end(); ++c) {
        auto other = c;
        ++other;
        while (other != contigs.end()) {
            if (c->add(*other)) {
                contigs.erase(other);
                other=c; // start from ctg again
            }
            ++other;
        }
    }

    // Find the contig with the most reads.
    // ----------
    const Contig* best_c = NULL;
    size_t best_n_reads = 0;
    for(Contig& c : contigs) {
        if (c.count() > best_n_reads) {
            best_c = &c;
            best_n_reads = best_c->count();
        }
    }

    // Build the MergedStack.
    // ----------
    locus.loc = best_c->loc;
    locus.add_consensus(string('N', best_c->len).c_str());
    for (const PStack* s : best_c->stacks)
        locus.utags.push_back(s->id);
    locus.count = best_n_reads;

    // Extend the stacks so that they all have the same span, and destroy the
    // PStacks that we couldn't use.
    // ----------
    for (PStack& p : pstacks) {
        if (best_c->stacks.count(&p))
            p.extend(best_c->loc, best_c->len);
        else
            p.clear();
    }
    pstacks.erase(remove_if( // note: The `Contig::stack`s are invalidated.
            pstacks.begin(),
            pstacks.end(),
            [] (PStack& s) {return s.seq == NULL;}
            ), pstacks.end());

    assert(pstacks.size() == best_c->stacks.size());

    return locus;
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
        {"version",      no_argument,       NULL,  1000},
        {"quiet",        no_argument, NULL, 'q'},
        {"num_threads",  required_argument, NULL, 'p'},
        {"prefix",       required_argument, NULL, 's'},
        {"paired_reads", required_argument, NULL, 'f'},
        {"filetype",     required_argument, NULL, 't'},
        {"model",        required_argument, NULL,  1001},
        {"alpha",        required_argument, NULL,  1002},
        {"bound_low",    required_argument, NULL,  1003},
        {"bound_high",   required_argument, NULL,  1004},
        {"bc_err_freq",  required_argument, NULL,  1005},
        {"bc_err_freq",  required_argument, NULL,  1005},
        {"debug_flags",  required_argument, NULL,  999},
        {0, 0, 0, 0}
    };

    int c;
    int long_options_i;
    while (true) {

        c = getopt_long(argc, argv, "hvs:f:t:p:q", long_options, &long_options_i);

        // Detect the end of the options.
        if (c == -1)
            break;

        switch (c) {
        case 'h':
            cout << help_string;
            exit(0);
            break;
        case 1000: //version
            cout << "pstacks_pe " << VERSION << "\n";
            exit(0);
            break;
        case 'q':
            quiet = true;
            break;
        case 'p':
            num_threads = atoi(optarg);
            break;
        case 's':
            prefix_path = optarg;
            break;
        case 'f':
            paired_alns_path = optarg;
            break;
        case 't':
            if (strcmp(optarg, "bowtie") == 0)
                in_file_type = FileT::bowtie;
            else if (strcmp(optarg, "sam") == 0)
                in_file_type = FileT::sam;
            else if (strcmp(optarg, "bam") == 0)
                in_file_type = FileT::bam;
            else if (strcmp(optarg, "tsv") == 0)
                in_file_type = FileT::tsv;
            else {
                cerr << "Error: Unknown file type '" << optarg << "'.\n";
                bad_args();
            }
            break;
        case 1001: //model
            if (strcmp(optarg, "snp") == 0) {
                model_type = snp;
            } else if (strcmp(optarg, "fixed") == 0) {
                model_type = ::fixed;
            } else if (strcmp(optarg, "bounded") == 0) {
                model_type = bounded;
            } else {
                cerr << "Error: Unknown model type '" << optarg << "'.\n";
                bad_args();
            }
            break;
        case 1002: //alpha
            alpha = atof(optarg);
            if (alpha != 0.1 && alpha != 0.05 && alpha != 0.01 && alpha != 0.001) {
                cerr << "SNP model alpha significance level must be either 0.1, 0.05, 0.01, or 0.001.\n";
                bad_args();
            }
            break;
        case 1003: //bound_low
            bound_low  = atof(optarg);
            if (bound_low < 0. || bound_low > 1.) {
                cerr << "SNP model lower bound must be between 0.0 and 1.0.\n";
                bad_args();
            }
            break;
        case 1004: //bound_high
            bound_high = atof(optarg);
            if (bound_high < 0. || bound_high > 1.) {
                cerr << "SNP model upper bound must be between 0.0 and 1.0.\n";
                bad_args();
            }
            break;
        case 1005: //bc_err_freq
            barcode_err_freq = atof(optarg);
            break;
        case 999: //debug_flags
        {
            stringstream ss (optarg);
            string s;
            while (getline(ss, s, ',')) {
                if (known_debug_flags.count(s)) {
                    debug_flags.insert(s);
                } else {
                    cerr << "? Error: Unknown error flag '" << s << "'.\n";
                    bad_args();
                }
            }
            break;
        }
        case '?':
            bad_args();
            break;
        default:
            cerr << "Option '--" << long_options[long_options_i].name << "' is deprecated.\n";
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

    if (paired_alns_path.empty()) {
        cerr << "You must specify the path to the paired reads alignments.\n";
        bad_args();
    }

    if (in_file_type == FileT::unknown) {
        if (paired_alns_path.length() >= 4 && paired_alns_path.substr(paired_alns_path.length() - 4) == ".bam") {
            in_file_type = FileT::bam;
        } else if (paired_alns_path.length() >= 4 && paired_alns_path.substr(paired_alns_path.length() - 4) == ".sam") {
            in_file_type = FileT::sam;
        } else {
            cerr << "Unable to guess the type of the input file, please specify it.\n";
            bad_args();
        }
    }

    if (model_type == ::fixed && barcode_err_freq == 0) {
        cerr << "You must specify the barcode error frequency.\n";
        bad_args();
    }
}

void report_options(ostream& os) {
    os << "Configuration for this run:\n";
    os << "  Sample prefix: '" << prefix_path << "'\n";
    os << "  Paired-end reads alignments: '" << paired_alns_path << "'\n";
    os << "  Paired-end reads alignment type: " << to_string(in_file_type) << "\n";
    os << "  Number of threads: " << num_threads << "\n";

    // Model.
    if (model_type == snp) {
        os << "  Model: snp\n"
           << "  Model alpha: " << alpha << "\n";
    } else if (model_type == bounded) {
        os << "  Model: snp\n"
           << "  Model alpha: " << alpha << "\n"
           << "  Model lower bound: " << bound_low << "\n"
           << "  Model higher bound: " << bound_high << "\n";
    } else if (model_type == ::fixed) {
        os << "  Model: fixed\n"
           << "  Model barcode err. prob.: " << barcode_err_freq << "\n";
    }

    if (!debug_flags.empty()) {
        os << "  DEBUG FLAGS:";
        for (const string& flag : debug_flags)
            os << " " << flag;
        os << "\n";
    }
}
