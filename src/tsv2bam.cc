#include <cmath>
#include <cstring>
#include <iostream>
#include <fstream>

#include <getopt.h>
#include <unistd.h>

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
#include "catalog_utils.cc"

void parse_command_line(int argc, char* argv[]);
void report_options(ostream& os);
int  run();
void run(const vector<int>& cloc_ids,
         const string& header_sq_lines,
         const unordered_map<int,size_t>& cloc_lengths,
         size_t sample_i, ostream& os
         );
void cigar_apply_to_locus(Locus* l, const Cigar& c);

//
// Argument globals.
//
bool quiet = false;
string in_dir;
int batch_id = -1;
vector<string> samples;
string popmap_path;
vector<string> pe_reads_paths;
FileT pe_reads_format;
int num_threads = 1;

bool dbg_reversed_pe_reads = false;

//
// Extra globals.
//
const string prog_name = "tsv2bam";
unique_ptr<LogAlterator> logger = NULL;

// main()
// ==========
int main(int argc, char* argv[]) {
    IF_NDEBUG_TRY

    // Parse arguments
    parse_command_line(argc, argv);

    // Open the log
    string lg_path = in_dir + prog_name + ".log";
    logger.reset(new LogAlterator(lg_path, quiet, argc, argv));
    report_options(cout);
    cout << "\n";

    if (!pe_reads_paths.empty())
        cout << "Paired-end reads files found, e.g. '" << pe_reads_paths.at(0) << "'.\n";

    // Initialize OPENMP.
    #ifdef _OPENMP
    omp_set_num_threads(num_threads);
    #endif

    // Run the program.
    int rv = run();
    if (rv != 0)
        return rv;

    // Cleanup.
    cout << "\n" << prog_name << " is done.\n";
    return 0;
    IF_NDEBUG_CATCH_ALL_EXCEPTIONS
}

int run() {

    //
    // Read the catalog.
    //
    cout << "Loading the catalog..." << endl;
    map<int, Locus*> catalog;
    string catalog_prefix = in_dir + "batch_" + to_string(batch_id) + ".catalog";
    bool dummy;
    int rv = load_loci(catalog_prefix, catalog, 0, false, dummy, false);
    if (rv != 1) {
        cerr << "Error: Unable to load catalog '" << catalog_prefix << "'.\n";
        throw exception();
    }

    //
    // Create the BAM target list.
    //
    stringstream header_sq_lines;
    vector<int> cloc_ids;
    unordered_map<int,size_t> cloc_lengths;
    cloc_ids.reserve(catalog.size());
    for (auto& cloc : catalog) {
        // Add the @SQ line.
        header_sq_lines << "@SQ\tSN:" << cloc.first
                        << "\tLN:" << cloc.second->len
                        << "\tnc:" << cloc.second->comp.size();

        if (!cloc.second->loc.empty()) {
            const PhyLoc& pos = cloc.second->loc;
            header_sq_lines << "\tpos:" << pos.chr()
                        << ':' << (pos.bp+1)
                        << ':' << (pos.strand == strand_plus ? '+' : '-');
        }
        header_sq_lines << "\n";
        cloc_ids.push_back(cloc.first);
        cloc_lengths[cloc.first] = cloc.second->len;
        delete cloc.second;
        cloc.second = NULL;
    }
    catalog.clear();

    //
    // Process the samples.
    //
    vector<string> outputs (samples.size());
    int omp_return = 0;
    #pragma omp parallel
    {
        #pragma omp for schedule(dynamic)
        for (size_t i=0; i<samples.size(); ++i) {
            if (omp_return != 0)
                continue;
            try {
                #pragma omp critical(cout)
                cout << "Processing sample '" << samples[i] << "'...\n" << flush;

                stringstream ss;
                run(cloc_ids, header_sq_lines.str(), cloc_lengths, i, ss);
                outputs[i] = ss.str();

            } catch (exception& e) {
                #pragma omp critical(exception)
                omp_return = stacks_handle_exceptions(e);
            }
        }
    }
    if (omp_return != 0)
        return omp_return;

    //
    // Write the saved logs.
    //
    for (const string& o : outputs)
        cout << o;

    return 0;
}

void run(const vector<int>& cloc_ids,
    const string& header_sq_lines,
    const unordered_map<int,size_t>& cloc_lengths,
    size_t sample_i,
    ostream& os
) {
    ostream& cout = os;
    cout << "\nSample '" << samples.at(sample_i) << "':\n";

    string prefix_path = in_dir + samples.at(sample_i);

    //
    // Read the sample's tsv files.
    //
    map<int, Locus*> sloci;
    vector<CatMatch*> matches;
    int sample_id = -1;
    {
        load_catalog_matches(prefix_path, matches, false);
        if (matches.empty()) {
            cerr << "Error: Unable to load matches from '"
                 << prefix_path + ".matches.tsv(.gz)'.\n";
            throw exception();
        } else if (matches[0]->batch_id != batch_id) {
            cerr << "Error: The catalog batch ID is " << batch_id
                 << " but matches file '" << prefix_path
                 << ".matches.tsv.gz' corresponds to batch ID "
                 << matches[0]->batch_id << ".\n";
            throw exception();
        }
        sample_id = matches[0]->sample_id;

        bool dummy;
        int rv = load_loci(prefix_path, sloci, 2, false, dummy, false);
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
        cout << bij_loci.size() << " sample loci ("
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
    // Align the sample loci to the catalog.
    //
    int prev_sloc = -1;
    Cigar c;
    for (const CatMatch* m : matches) {
        assert(m->tag_id >= prev_sloc); // Matches files are sorted by sloc_id
        if (m->tag_id == prev_sloc)
            // We have already aligned this sample locus.
            continue;
        prev_sloc = m->tag_id;

        auto itr = sloci.find(m->tag_id);
        if (itr == sloci.end())
            // Loci that aren't bijective with the catalog have been removed.
            continue;
        Locus* l = itr->second;

        if (cloc_lengths.at(m->cat_id) == l->len)
            // No alignment to do.
            continue;
        if (m->cigar == NULL || m->cigar[0] == '\0') {
            c.clear();
            c.push_back({'M', l->len});
        } else {
            parse_cigar(m->cigar, c);
            // `apply_cigar_to_seq` returns a padded sequence, so insertions would be
            // a problem. However, "match" CIGARs should in principle not contain them.
            assert(cigar_is_MDI(c) && strchr(m->cigar, 'I') == NULL);
        }
        size_t cig_len = cigar_length_ref(c);
        if (cig_len < cloc_lengths.at(m->cat_id))
            cigar_extend_right(c, cloc_lengths.at(m->cat_id) - cig_len);
        cigar_apply_to_locus(l, c);
        assert(l->len == cloc_lengths.at(m->cat_id));
    }

    //
    // Sort the loci by catalog ID & prepare the loading of paired-end reads.
    //
    struct Loc {
        int cloc_id;
        Locus* sloc; // 'Sample' locus (includes the fw reads in Locus::comp and Locus::reads)
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
    if(!pe_reads_paths.empty()) {
        string pe_reads_path = pe_reads_paths.at(sample_i);

        // Initialize the paired-end read structures.
        for (Loc& loc : sorted_loci)
            loc.pe_reads = new map<DNASeq4, vector<string>>();

        // Open the paired-end reads file.
        Input* pe_reads_f = NULL;
        switch (pe_reads_format) {
        case FileT::gzfastq: pe_reads_f = new GzFastq(pe_reads_path.c_str()); break;
        case FileT::gzfasta: pe_reads_f = new GzFasta(pe_reads_path.c_str()); break;
        case FileT::fastq: pe_reads_f = new Fastq(pe_reads_path.c_str()); break;
        case FileT::fasta: pe_reads_f = new Fasta(pe_reads_path.c_str()); break;
        case FileT::bam: pe_reads_f = new Bam(pe_reads_path.c_str()); break;
        default: DOES_NOT_HAPPEN; break;
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
        cout << "Found a paired-end read for " << n_used_reads
             << " (" << as_percentage((double) n_used_reads / readname_to_loc.size())
             << ") of the assembled forward reads." << endl;

        delete pe_reads_f;
    }

    //
    // Write the BAM header.
    //
    string header;
    header += "@HD\tVN:1.5\tSO:coordinate\n";
    header += "@RG\tID:" + to_string(sample_id) + "\tSM:" + samples.at(sample_i) + "\tid:" + to_string(sample_id) + '\n';
    header += header_sq_lines;

    //
    // Open the BAM file.
    //
    string matches_bam_path = prefix_path + ".matches.bam";
    Bam bam_f (matches_bam_path, BamHeader(header));

    //
    // Write the BAM records.
    //
    {
        const vector<int>& targets = cloc_ids;
        size_t target_i = 0;
        BamRecord rec;
        Cigar cig;
        for (auto& loc : sorted_loci) {
            while (targets[target_i] != loc.cloc_id)
                ++target_i;

            if(!pe_reads_paths.empty()) {
                // Write the forward reads.
                for (size_t j=0; j<loc.sloc->comp.size();++j) {
                    string name (loc.sloc->comp[j]);
                    strip_read_number(name);

                    const char* seq = loc.sloc->reads[j];
                    size_t seq_len = strlen(seq);
                    cig.assign({{'M', seq_len}});
                    rec.assign(
                            name,
                            BAM_FPAIRED | BAM_FREAD1,
                            target_i,
                            0,
                            cig,
                            DNASeq4(seq, seq_len),
                            sample_id
                            );
                    bam_f.write(rec);
                }

                // Write the paired-end reads.
                for (const pair<DNASeq4, vector<string>>& stack : *loc.pe_reads) {
                    for (const string& name : stack.second) {
                        cig.assign({{'X', 1}, {'S', stack.first.length()-1}});
                        rec.assign(
                                name,
                                BAM_FPAIRED | BAM_FREAD2 | BAM_FREVERSE,
                                target_i,
                                0,
                                cig,
                                dbg_reversed_pe_reads ? stack.first : stack.first.rev_compl(),
                                sample_id
                                );
                        bam_f.write(rec);
                    }
                }

            } else {
                // Write the (possibly pooled) reads.
                // In this case we do not expect the read names to have a
                // specific suffix, and we do not set the FREAD1 flag.
                for (size_t j=0; j<loc.sloc->comp.size();++j) {
                    const char* name = loc.sloc->comp[j];
                    const char* seq = loc.sloc->reads[j];
                    size_t seq_len = strlen(seq);
                    cig.assign({{'M', seq_len}});
                    rec.assign(
                            string(name),
                            0,
                            target_i,
                            0,
                            cig,
                            DNASeq4(seq, seq_len),
                            sample_id
                            );
                    bam_f.write(rec);
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
    if(!pe_reads_paths.empty())
        for (Loc& loc : sorted_loci)
            delete loc.pe_reads;
}

void cigar_apply_to_locus(Locus* l, const Cigar& c) {
    // Align the consensus & model.
    // (@Nick Sep2017: Changing `con` & `model` is for consistency, I don't
    // think these variables are actually used later.)
    size_t cloc_len = cigar_length_ref(c);
    string new_con = apply_cigar_to_seq(l->con, c);
    string new_model = apply_cigar_to_model_seq(l->model, c);
    assert(new_con.length() == cloc_len);
    assert(new_model.length() == cloc_len);

    if (l->len != cloc_len) {
        assert(cloc_len > l->len); // As matches files don't contain I operations.
        delete[] l->con;
        delete[] l->model;
        l->con = new char[cloc_len+1];
        l->model = new char[cloc_len+1];
    }
    strncpy(l->con, new_con.c_str(), cloc_len+1);
    strncpy(l->model, new_model.c_str(), cloc_len+1);

    // Align the reads.
    for (char*& r : l->reads) {
        string new_r = apply_cigar_to_seq(r, c);
        assert(new_r.length() == cloc_len);
        if (l->len != cloc_len) {
            delete[] r;
            r = new char[cloc_len+1];
        }
        strncpy(r, new_r.c_str(), cloc_len+1);
    }
    l->len = cloc_len;
}

void parse_command_line(int argc, char* argv[]) {

    const string help_string = string() +
            prog_name + " " + VERSION  + "\n" +
            prog_name + " -P stacks_dir -M popmap [-b batch_id] [-R paired_reads_dir]\n" +
            prog_name + " -P stacks_dir -s sample [-s sample ...] [-b batch_id] [-R paired_reads_dir]\n" +
            "\n"
            "  -P,--in-dir: input directory.\n"
            "  -M,--popmap: population map.\n"
            "  -s,--sample: name of one sample.\n"
            "  -R,--pe-reads-dir: directory where to find the paired-end reads files (in fastq/fasta/bam (gz) format).\n"
            "  -b: catalog batch ID (default: guess).\n"
            "  -t: number of threads to use (default: 1).\n"
            "\n"
#ifdef DEBUG
            "Debug options:\n"
            "  --dbg-reversed-pe-reads: for simulated paired-end reads written on the same\n"
            "                           strand as the forward ones\n"
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
        {"in-dir",       required_argument, NULL, 'P'},
        {"popmap",       required_argument, NULL, 'M'},
        {"sample",       required_argument, NULL, 's'},
        {"pe-reads-dir", required_argument, NULL, 'R'},
        {"batch-id",     required_argument, NULL, 'b'},
        {"threads",      required_argument, NULL, 't'},
        {"dbg-reversed-pe-reads", no_argument, NULL, 2000},
        {0, 0, 0, 0}
    };

    string pe_reads_dir;

    int c;
    int long_options_i;
    while (true) {

        c = getopt_long(argc, argv, "hqt:P:b:M:s:R:", long_options, &long_options_i);

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
        case 'P':
            in_dir = optarg;
            if (in_dir.back() != '/')
                in_dir += '/';
            break;
        case 'M':
            popmap_path = optarg;
            break;
        case 's':
            samples.push_back(optarg);
            break;
        case 'R':
            pe_reads_dir = optarg;
            if (pe_reads_dir.back() != '/')
                pe_reads_dir += '/';
            break;
        case 'b':
            batch_id = is_integer(optarg);
            if (batch_id < 0) {
                cerr << "Error: Illegal -b option value '" << optarg << "'.\n";
                bad_args();
            }
            break;
        case 't':
            num_threads = is_integer(optarg);
            if (num_threads < 0) {
                cerr << "Error: Illegal -t option value '" << optarg << "'.\n";
                bad_args();
            }
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

    // -P
    if (in_dir.empty()) {
        cerr << "Error: An input directory is required (-P).\n";
        bad_args();
    }

    // -M & -s
    if (popmap_path.empty() && samples.empty()) {
        cerr << "Error: One of -M or -s should be provided.\n";
        bad_args();
    } else if (!popmap_path.empty() && !samples.empty()) {
        cerr << "Error: only one of -M or -s should be provided.\n";
        bad_args();
    }

    // -b
    if (batch_id < 0) {
        vector<int> cat_ids = find_catalogs(in_dir);
        if (cat_ids.empty()) {
            cerr << "Error: Unable to find a catalog in '" << in_dir << "'.\n";
            bad_args();
        } else if (cat_ids.size() == 1) {
            batch_id = cat_ids[0];
        }  else {
            cerr << "Error: Input directory contains several catalogs, please specify -b.\n";
            bad_args();
        }
    }

    //
    // Process arguments.
    //
    if (!popmap_path.empty()) {
           MetaPopInfo mpopi;
           mpopi.init_popmap(popmap_path);
           for (const Sample& s : mpopi.samples())
               samples.push_back(s.name);
    }

    if (!pe_reads_dir.empty()) {
        // Find the first paired-end read file.
        string ext;
        for (const string& e : {".fq.gz", ".fastq.gz", ".fq", ".fastq", ".fa.gz", ".fasta.gz", ".fa", ".fasta", ".bam"}) {
            if (access((pe_reads_dir + samples.at(0) + ".2" + e).c_str(), F_OK) == 0) {
                // File exists.
                ext = e;
                break;
            }
        }
        if (ext.empty()) {
            cerr << "Error: Unable to find the first paired-end reads file at '" << pe_reads_dir+samples.at(0) << ".2.*'\n";
            throw exception();
        }
        // Record the paired-end reads paths.
        for (const string& sample : samples)
            pe_reads_paths.push_back(pe_reads_dir + sample + ".2" + ext);
        pe_reads_format = guess_file_type(pe_reads_paths.at(0));
        if (pe_reads_format == FileT::unknown)
            DOES_NOT_HAPPEN;
    }
}

void report_options(ostream& os) {
    os << "Configuration for this run:\n";
    os << "  Stacks directory: '" << in_dir << "'\n"
       << "  Batch ID: " << batch_id << "\n";
    if (!popmap_path.empty())
        os << "  Population map: '" << popmap_path << "'\n";
    os << "  Num. samples: " << samples.size() << "\n";
    if (!pe_reads_paths.empty())
        os << "  Paired-end reads directory: '"
           << pe_reads_paths[0].substr(0, pe_reads_paths[0].find_last_of('/')) << "/'\n";
    if (num_threads > 1)
        os << "  Multithreaded.\n";
}
