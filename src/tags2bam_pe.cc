#include <cmath>
#include <cstring>
#include <iostream>
#include <fstream>

#include <getopt.h>

#include "sam.h"

#include "constants.h"
#include "sql_utilities.h"
#include "log_utils.h"
#include "gzFastq.h"
#include "BamI.h"
#include "DNASeq4.h"

#include "tags2bam_pe.h"

using namespace std;

//
// Argument globs.
//
bool quiet = false;
string prefix_path;
string pe_reads_path;
const FileT pe_reads_format = FileT::gzfastq; // xxx Const, because I don't have guess_file_type() in branch nick yet.

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
    string lg_path = prefix_path + ".tags2bam_pe.log";
    if(!quiet)
        cout << "Logging to '" << lg_path << "'." << endl;
    lg = new LogAlterator(lg_path, quiet);
    init_log(lg->l, argc, argv);
    report_options(cout);
    cout << "\n";

    // Read the sample's stacks files.
    vector<int> cloci; // catalog locus IDs
    unordered_map<string, size_t> read_name_to_loc; // map of {read name, loc index}
    int sample_id;
    read_sample_files(cloci, read_name_to_loc, sample_id);

    // Read the paired-end reads.
    vector<map<DNASeq4, vector<string> > > pe_reads_by_loc (cloci.size());
    load_pe_reads(pe_reads_by_loc, read_name_to_loc);

    // Sort the loci by catalog ID (they are sorted by sample ID).
    map<int, size_t> sorted_loci; // (cloc_id, index)
    for (size_t i=0; i<cloci.size(); ++i)
        sorted_loci.insert({cloci[i], i});

    // Write the BAM file.
    write_bam_file(sorted_loci, pe_reads_by_loc, sample_id);

    cout << "tags2bam_pe is done." << endl;
    return 0;
}

void read_sample_files(
        std::vector<int>& cloci,
        std::unordered_map<string, size_t>& read_name_to_loc,
        int& sample_id
        ){

    // Read the matches file.
    cout << "Reading the matches file..." << endl;
    vector<CatMatch*> matches;
    load_catalog_matches(prefix_path, matches);
    if (matches.empty()) {
        cerr << "Error: Unable to load matches from '"
             << prefix_path + ".matches.tsv(.gz)'.\n";
        throw exception();
    }

    // Retrieve the sample ID.
    sample_id = matches[0]->sample_id;

    // Retrieve the list of bijective loci.
    unordered_map<int, int> sloc_to_cloc;
    vector<pair<int, int> > bij_loci = retrieve_bijective_loci(matches);
    sloc_to_cloc.reserve(bij_loci.size());
    sloc_to_cloc.insert(bij_loci.begin(), bij_loci.end());

    // Read the tags file.
    cout << "Reading the tags file..." << endl;
    map<int, Locus*> sloci;
    bool tmpbool;
    int rv = load_loci(prefix_path, sloci, 1, false, tmpbool);
    if(rv != 1) {
        cerr << "Error: could not find stacks files '" << prefix_path << ".*' (tags, snps and/or alleles).\n";
        throw exception();
    }

    // Discard loci that aren't bijective with the catalog.
    for (auto sloc=sloci.begin(); sloc!=sloci.end();) {
        if (sloc_to_cloc.count(sloc->second->id)) {
            sloc++;
        } else {
            delete sloc->second;
            sloci.erase(sloc++);
        }
    }

    // Link paired reads to catalog loci.
    cloci.reserve(sloci.size());
    size_t i=0;
    for (auto& sloc : sloci) {
        const Locus& sloc2 = *sloc.second;

        for (const char* read_name : sloc2.comp) {
            string pe_read (read_name);
            convert_fw_read_name_to_paired(pe_read);
            read_name_to_loc.insert({move(pe_read), i});
        }

        cloci.push_back(sloc_to_cloc.at(sloc2.id));
        i++;
    }

    // Clean up.
    for (CatMatch* m : matches)
        delete m;

    for (auto& sloc : sloci)
        delete sloc.second;
}

void load_pe_reads(
        vector<map<DNASeq4, vector<string> > >& pe_reads_by_loc,
        const unordered_map<string, size_t>& read_name_to_loc
        ) {

    cout << "Reading the paired-end reads..." << endl;

    Input* pe_reads_f;
    if (pe_reads_format == FileT::gzfastq)
        pe_reads_f = new GzFastq(pe_reads_path.c_str());
    else
        throw exception(); // xxx

    size_t n_used_reads = 0;
    Seq seq;
    seq.id   = new char[id_len]; // Necessary or GzFastq will segfault.
    seq.seq  = new char[max_len];
    seq.qual = new char[max_len];
    while(pe_reads_f->next_seq(seq)) {
        string id (seq.id);
        DNASeq4 seq4 (seq.seq); // xxx could assign

        auto loc = read_name_to_loc.find(id);
        if (loc == read_name_to_loc.end())
            continue;

        ++n_used_reads;

        map<DNASeq4, vector<string> >& loc_stacks = pe_reads_by_loc.at(loc->second);
        auto& stack = *loc_stacks.insert({move(seq4), vector<string>()}).first;
        stack.second.push_back(move(id));
    }

    if (n_used_reads == 0) {
        cerr << "Error: Failed to find any matching paired-end reads in '" << pe_reads_path << "'.\n";
        throw exception();
    }

    cout << "Assigned " << n_used_reads << " paired-end reads to catalog loci." << endl;
}

void write_bam_file(
        const map<int, size_t>& sorted_loci,
        const vector<map<DNASeq4, vector<string> > >& pe_reads_by_loc,
        int sample_id
        ){

    cout << "Writing the bam file..." << endl;

    string bam_path = prefix_path + ".pe_reads.bam";
    htsFile* bam_f = hts_open(bam_path.c_str(), "wb");

    // Write the header.
    // ----------

    // We want a header with the following structure:
    // @HD VN:1.5 SO:coordinate
    // @RG ID:sample_id SM:sample_name
    // @SQ SN:cloc_id LN:(2^31-1) -- No length.
    // @SQ etc.

    string sample_name = prefix_path.substr(prefix_path.find_last_of('/')+1);
    string header_text = string() +
            "@HD\tVN:1.5\tSO:coordinate\n"
            "@RG\tID:" + to_string(sample_id) + "\tSM:" + sample_name + "\n";

    vector<pair<string, uint32_t> > chrs;
    for (auto& loc : sorted_loci)
        chrs.push_back({to_string(loc.first), pow<size_t>(2,31)-1});

    write_bam_header(bam_f, header_text, chrs);

    // Write the records.
    // ----------

    size_t loc_i = 0; // Locus index; we use the same loop as for the header.
    bam1_t* r = bam_init1(); // As we're going to reuse the exact same fields for all the reads.

    for (auto& loc : sorted_loci) {
        for (auto& stack : pe_reads_by_loc.at(loc.second)) {
            for (const string& read_name : stack.second) {

                const DNASeq4& seq = stack.first;

                // Prepare the `aux` array.
                size_t l_aux = 3 + sizeof(int32_t);
                uchar aux[l_aux];
                aux[0] = 'R';
                aux[1] = 'G';
                aux[2] = 'i';
                *(int32_t*)(aux+3) = sample_id;

                // bam1_t::core
                bam1_core_t& c = r->core;
                c.tid = loc_i;
                c.pos = 0;
                c.qual = 255;
                c.l_qname = read_name.length() + 1; //n.b. `l_qname` includes the null character.
                c.flag = 0;
                c.n_cigar = 1; //An aligned sequence should have at least one cigar op.
                c.l_qseq = seq.length();
                c.mtid = -1;
                c.mpos = -1;

                // bam1_t::data
                // Htslib says: "bam1_t::data -- all variable-length data, concatenated;
                // structure: qname-cigar-seq-qual-aux, concatenated".

                r->l_data = r->core.l_qname + r->core.n_cigar*sizeof(uint32_t) + seq.vsize() + seq.length() + l_aux;
                r->m_data = r->l_data;
                r->data = new uchar[r->m_data];

                uchar* p = r->data;
                strcpy((char*)p, read_name.c_str());
                p += r->core.l_qname;
                *(uint32_t*)p = seq.length() <<BAM_CIGAR_SHIFT; // Barcodes have their length on the 24 high bits & op on the low 8 bits.
                *(uint32_t*)p |= 4; // S, c.f. Spec.
                p += r->core.n_cigar*sizeof(uint32_t);
                memcpy(p, seq.vdata(), seq.vsize()); // n.b. `sizeof(DiNuc)==1`
                p += seq.vsize();
                memset(p, 0xFF, seq.length());
                p += seq.length();
                memcpy(p, aux, l_aux);

                // bam1_t::core.bin
                // c.f. `sam_parse1()`; I have no idea what this is.
                uint32_t* cigar = (uint32_t*)(r->data + r->core.l_qname);
                c.bin = hts_reg2bin(c.pos, c.pos + bam_cigar2rlen(c.n_cigar, cigar), 14, 5);

                int rv = bam_write1(bam_f->fp.bgzf, r);
                if (rv == -1) {
                    cerr << "Error: Writing of BAM record failed.\n";
                    throw exception();
                }
            }
        }
        loc_i++;
    }

    hts_close(bam_f);
}

const string help_string = string() +
        "tags2bam_pe " + VERSION  + "\n"
        "tags2bam_pe -s sample_prefix -f pe_reads_path \n"
        "\n"
        "  -s,--prefix: prefix path for the sample.\n"
        "  -f,--pe_reads: path to paired-end reads.\n"
        "\n"
        ;

void bad_args() {
    cerr << help_string;
    exit(-1);
}

void report_options(ostream& os) {
    os << "Configuration for this run:\n"
       << "  Sample prefix: '" << prefix_path << "'\n"
       << "  Paired-end reads path: '" << pe_reads_path << "'\n"
       << "  Paired-end reads format: " + to_string(pe_reads_format) + "\n";
}

void parse_command_line(int argc, char* argv[]) {

    static const option long_options[] = {
        {"help",         no_argument,       NULL, 'h'},
        {"quiet",        no_argument,       NULL, 'q'},
        {"version",      no_argument,       NULL,  1000},
        {"prefix",       required_argument, NULL, 's'},
        {"pe_reads",     required_argument, NULL, 'f'},
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
            cout << "pstacks_pe " << VERSION << "\n";
            exit(0);
            break;
        case 's':
            prefix_path = optarg;
            break;
        case 'f':
            pe_reads_path = optarg;
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

    if (pe_reads_path.empty()) {
        cerr << "Error: A path to paired-end reads is required (-f).\n";
        bad_args();
    }
}
