#include <cmath>
#include <cstring>
#include <iostream>
#include <fstream>

#include <getopt.h>

#include "constants.h"
#include "sql_utilities.h"
#include "log_utils.h"
#include "locus.h"
#include "gzFastq.h"
#include "gzFasta.h"
#include "FastqI.h"
#include "FastaI.h"
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
bool dbg_reversed_pe_reads = false;

//
// Extra globals.
//
const string prog_name = "tsv2bam";
LogAlterator* logger = NULL;

// main()
// ==========
int main(int argc, char* argv[]) {
    IF_NDEBUG_TRY

    // Parse arguments
    parse_command_line(argc, argv);

    // Open the log
    string lg_path = prefix_path + "." + prog_name + ".log";
    logger = new LogAlterator(lg_path, quiet, argc, argv);
    report_options(cout);
    cout << "\n";

    // Run the program.
    run();

    // Cleanup.
    cout << prog_name << " is done.\n";
    delete logger;
    return 0;
    IF_NDEBUG_CATCH_ALL_EXCEPTIONS
}

void run() {

    //
    // Read the sample's tsv files.
    //
    map<int, Locus*> sloci;
    vector<CatMatch*> matches;
    int sample_id = -1;
    {
        cout << "Reading the matches file..." << endl;
        load_catalog_matches(prefix_path, matches);
        if (matches.empty()) {
            cerr << "Error: Unable to load matches from '"
                 << prefix_path + ".matches.tsv(.gz)'.\n";
            throw exception();
        }
        sample_id = matches[0]->sample_id;

        cout << "Reading the tags file..." << endl;
        bool tmpbool;
        int rv = load_loci(prefix_path, sloci, 2, false, tmpbool);
        if(rv != 1) {
            cerr << "Error: Could not find stacks files '" << prefix_path << ".*' (tags, snps and/or alleles).\n";
            throw exception();
        }
    }

    //
    // Retrieve the list of bijective loci.
    //
    unordered_map<int, int> sloc_to_cloc;
    {
        vector<pair<int, int> > bij_loci = retrieve_bijective_loci(matches);
        cerr << bij_loci.size() << " sample loci ("
             << as_percentage((double) bij_loci.size() / sloci.size())
             << " of " << sloci.size()
             << ") had a one-to-one relationship with the catalog.\n";
        if (bij_loci.empty()) {
            cerr << "Error: No usable matches to the catalog.\n";
            throw exception();
        }
        sloc_to_cloc.reserve(bij_loci.size());
        sloc_to_cloc.insert(bij_loci.begin(), bij_loci.end());
    }

    //
    // Discard loci that aren't bijective with the catalog.
    //
    for (auto sloc=sloci.begin(); sloc!=sloci.end();) {
        if (sloc_to_cloc.count(sloc->second->id)) {
            assert(!sloc->second->blacklisted);
            ++sloc;
        } else {
            delete sloc->second;
            sloci.erase(sloc++);
        }
    }

    //
    // Sort the loci by catalog ID & prepare the loading of paired-end reads.
    //
    struct Loc {
        int cloc_id;
        Locus* sloc;
        map<DNASeq4, vector<string>>* pe_reads;

        bool operator< (const Loc& other) const {return cloc_id < other.cloc_id;}
    };

    vector<Loc> sorted_loci;
    {
        for (auto& sloc : sloci)
            sorted_loci.push_back({sloc_to_cloc.at(sloc.second->id), sloc.second, NULL});
        std::sort(sorted_loci.begin(), sorted_loci.end());
    }

    //
    // Load the paired-end reads, if any.
    //
    if(!pe_reads_path.empty()) {
        cerr << "Loading paired-end reads...\n";

        // Initialize the paired-end read structures.
        for (Loc& loc : sorted_loci)
            loc.pe_reads = new map<DNASeq4, vector<string>>();

        // Open the paired-end reads file.
        Input* pe_reads_f = NULL;
        switch (guess_file_type(pe_reads_path)) {
        case FileT::gzfastq: pe_reads_f = new GzFastq(pe_reads_path.c_str()); break;
        case FileT::gzfasta: pe_reads_f = new GzFasta(pe_reads_path.c_str()); break;
        case FileT::fastq: pe_reads_f = new Fastq(pe_reads_path.c_str()); break;
        case FileT::fasta: pe_reads_f = new Fasta(pe_reads_path.c_str()); break;
        case FileT::bam: pe_reads_f = new Bam(pe_reads_path.c_str()); break;
        default:
            cerr << "Error: '" << pe_reads_path << "' has unexpected format '"
                 << to_string(guess_file_type(pe_reads_path)) << "'.\n";
            throw exception();
            break;
        }

        // Get the read-name-to-locus map.
        unordered_map<string, Loc*> readname_to_loc;
        for (Loc& loc : sorted_loci) {
            for (const char* read_name : loc.sloc->comp) {
                string read_name2 (read_name);
                strip_read_number(read_name2);
                readname_to_loc.insert({move(read_name2), &loc});
            }
        }

        // Load the reads.
        size_t n_pe_reads = 0;
        size_t n_used_reads = 0;
        Seq seq;
        seq.id   = new char[id_len];
        seq.seq  = new char[max_len];
        seq.qual = new char[max_len];
        while(pe_reads_f->next_seq(seq)) {
            ++n_pe_reads;
            string id (seq.id);
            strip_read_number(id);
            auto loc_it = readname_to_loc.find(id);
            if (loc_it == readname_to_loc.end())
                continue;

            ++n_used_reads;
            map<DNASeq4, vector<string> >& loc_stacks = *loc_it->second->pe_reads;
            pair<const DNASeq4,vector<string>>& stack = *loc_stacks.insert({DNASeq4(seq.seq), vector<string>()}).first;
            stack.second.push_back(move(id));
        }
        if (n_used_reads == 0) {
            cerr << "Error: Failed to find any matching paired-end reads in '" << pe_reads_path << "'.\n";
            throw exception();
        }
        cerr << "Found a paired-end read for " << n_used_reads
             << " (" << as_percentage((double) n_used_reads / readname_to_loc.size())
             << ") of the assembled forward reads." << endl;

        delete pe_reads_f;
    }

    //
    // Open the BAM file.
    //
    htsFile* bam_f = NULL;
    {
        string matches_bam_path = prefix_path + ".matches.bam";
        cerr << "Outputing to '" << matches_bam_path << "'\n";
        bam_f = hts_open(matches_bam_path.c_str(), "wb");
    }

    //
    // Write the BAM header.
    //
    vector<int> targets; // Catalog loci.
    {
        cerr << "Writing the header...\n";
        string sample_name = prefix_path.substr(prefix_path.find_last_of('/')+1);
        stringstream header_text;
        header_text << "@HD\tVN:1.5\tSO:coordinate\n"
                    << "@RG\tID:" << sample_id << "\tSM:" << sample_name << "\tid:" << sample_id << "\n";

        // For samtools-merge to be happy, we must provide the complete list of
        // catalog loci. So we actually have to load the catalog.
        {
            int batch_id = matches[0]->batch_id;
            string dir_path = prefix_path.substr(0, prefix_path.find_last_of('/')+1);
            if (dir_path.empty())
                dir_path = "./";
            string catalog_prefix = dir_path + "batch_" + to_string(batch_id) + ".catalog";

            map<int, Locus*> catalog;
            bool tmp;
            int rv = load_loci(catalog_prefix, catalog, 0, false, tmp);
            if (rv != 1) {
                cerr << "Error: Unable to load catalog '" << catalog_prefix << "'.\n";
                throw exception();
            }

            targets.reserve(catalog.size());
            for (auto& cloc : catalog) {
                // Add the @SQ line. Length is mandatory; all loci declare to be 10000bp.
                header_text << "@SQ\tSN:" << cloc.first << "\tLN:10000";
                if (cloc.second->loc.chr != NULL) {
                    const PhyLoc& pos = cloc.second->loc;
                    header_text << "\tpos1:" << pos.chr
                                << ':' << (pos.bp+1)
                                << ':' << (pos.strand == strand_plus ? '+' : '-');
                }
                header_text << "\n";
                targets.push_back(cloc.first);
                delete cloc.second;
            }
        }

        // Write it.
        write_bam_header(bam_f, header_text.str());
    }

    //
    // Write the BAM records.
    //
    {
        cerr << "Writing the records...\n";
        size_t target_i = 0;
        BamRecord rec;
        for (auto& loc : sorted_loci) {
            while (targets[target_i] != loc.cloc_id)
                ++target_i;

            if(!pe_reads_path.empty()) {
                // Write the assembled reads.
                for (size_t j=0; j<loc.sloc->comp.size();++j) {
                    string name (loc.sloc->comp[j]);
                    strip_read_number(name);

                    const char* seq = loc.sloc->reads[j];
                    size_t seq_len = strlen(seq);
                    rec.assign(
                            name,
                            BAM_FREAD1,
                            target_i,
                            0,
                            {{'M', seq_len}},
                            DNASeq4(seq, seq_len),
                            sample_id
                            );
                    rec.write_to(bam_f);
                }

                // Write the paired-end reads.
                for (const pair<DNASeq4, vector<string>>& stack : *loc.pe_reads) {
                    for (const string& name : stack.second) {
                        rec.assign(
                                name,
                                BAM_FREAD2 | BAM_FREVERSE,
                                target_i,
                                0,
                                {{'S', stack.first.length()}},
                                dbg_reversed_pe_reads ? stack.first : stack.first.rev_compl(),
                                sample_id
                                );
                        rec.write_to(bam_f);
                    }
                }

            } else {
                // Write the assembled reads.
                // But in this case we do not expect the read names to have a
                // specific suffix, and we do not set the FREAD1 flag.
                for (size_t j=0; j<loc.sloc->comp.size();++j) {
                    const char* name = loc.sloc->comp[j];
                    const char* seq = loc.sloc->reads[j];
                    size_t seq_len = strlen(seq);
                    rec.assign(
                            string(name),
                            0,
                            target_i,
                            0,
                            {{'M', seq_len}},
                            DNASeq4(seq, seq_len),
                            sample_id
                            );
                    rec.write_to(bam_f);
                }
            }
        }
    }

    //
    // Cleanup.
    //
    for (CatMatch* m : matches)
        delete m;
    for (auto& sloc : sloci)
        delete sloc.second;
    if(!pe_reads_path.empty())
        for (Loc& loc : sorted_loci)
            delete loc.pe_reads;
    hts_close(bam_f);
}

void parse_command_line(int argc, char* argv[]) {

    const string help_string = string() +
            prog_name + " " + VERSION  + "\n" +
            prog_name + " -s sample_prefix [-f paired_reads]\n"
            "\n"
            "  -s,--prefix: prefix path for the sample. (Required.)\n"
            "  -f,--pe-reads: path to the file containing the sample's paired-end read sequences, if any.\n"
            "\n"
#ifdef DEBUG
            "Debug options:\n"
            "  --reversed-pe-reads: for simulated paired-end reads written on the same\n"
            "                       strand as the forward ones\n"
            "\n"
#endif
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
        {"reversed-pe-reads", no_argument,  NULL, 2000},
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
        case 2000: // reversed-pe-reads
            dbg_reversed_pe_reads = true;
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

    const set<FileT> accepted_formats {FileT::gzfastq, FileT::fastq, FileT::gzfasta, FileT::fasta, FileT::bam};
    if (!pe_reads_path.empty() && !accepted_formats.count(guess_file_type(pe_reads_path))) {
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
