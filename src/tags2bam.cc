#include <cmath>
#include <cstring>
#include <iostream>
#include <fstream>

#include <getopt.h>

#include "sam.h"

#include "constants.h"
#include "sql_utilities.h"
#include "log_utils.h"
#include "BamI.h"
#include "DNASeq4.h"

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
    string lg_path = prefix_path + ".tags2bam.log";
    if(!quiet)
        cout << "Logging to '" << lg_path << "'." << endl;
    lg = new LogAlterator(lg_path, quiet);
    init_log(lg->l, argc, argv);
    report_options(cout);
    cout << "\n";

    // Read the sample's files.
    map<int, Locus*> sloci;
    unordered_map<int, int> sloc_to_cloc;
    int sample_id;
    read_sample_files(sloci, sloc_to_cloc, sample_id);

    // Sort the loci by catalog ID.
    map<int, Locus*> sorted_loci; // (cloc_id, sloc)
    for (auto sloc : sloci)
        sorted_loci.insert({sloc_to_cloc.at(sloc.second->id), sloc.second});

    // Write the BAM file.
    write_bam_file(sloci, sample_id);

    cout << "tags2bam is done." << endl;
    return 0;
}

void read_sample_files(map<int, Locus*>& sloci, unordered_map<int, int>& sloc_to_cloc, int& sample_id) {
    // Read the matches file.
    cout << "Reading the matches file..." << endl;
    vector<CatMatch*> matches;
    load_catalog_matches(prefix_path, matches);
    if (matches.empty()) {
        cerr << "Error: Unable to load matches from '"
             << prefix_path + ".matches.tsv(.gz)'.\n";
        throw exception();
    }

    sample_id = matches[0]->sample_id;

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

void write_bam_file(const map<int, Locus*>& sorted_loci, int sample_id) {

    cout << "Writing the bam file..." << endl;

    string bam_path = prefix_path + ".matches.bam";
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
        const Locus* sloc = loc.second;
        for (size_t j=0; j<sloc->comp.size();++j) {
            const char* name = sloc->comp[j];
            DNASeq4 seq (sloc->reads[j], strlen(sloc->reads[j]));

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
            c.l_qname = strlen(name) + 1; //n.b. `l_qname` includes the null character.
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
            strcpy((char*)p, name);
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
        loc_i++;
    }

    hts_close(bam_f);
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
