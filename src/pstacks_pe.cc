#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <map>
#include <unordered_set>
#include <unordered_map>

#include "input.h"
#include "BamI.h"
#include "sql_utilities.h"

#include "stacks.h"
#include "models.h"

#include "pstacks_base.h"
//#include "pstacks_pe.h"

using namespace std;

// From pstacks {
const int barcode_size = 5;

FileT  in_file_type;
//string in_file;
//FileT  out_file_type;
//string out_path;
int    sql_id        = 0;
int    min_stack_cov = 3;
int    num_threads   = 1;

modelt model_type         = snp;
double alpha              = 0.05;
double bound_low          = 0.0;
double bound_high         = 1.0;
double p_freq             = 0.5;
double barcode_err_freq   = 0.0;
double heterozygote_limit = -3.84;
double homozygote_limit   =  3.84;
// end pstacks }

/*
 * pstacks_pe
 * **********
 */

int batch_id;
string prefix_path;
string first_reads_path;
string paired_reads_path;
string paired_alns_path;

/* CLocReadSet
 * =========== */
struct ReadsByCLoc {
    typedef map<const DNANSeq*, vector<Seq> > CLocReadSet;

    unordered_set<const DNANSeq*> unique_seqs;
    vector<CLocReadSet> readsets;

    size_t n_used_reads;

    // Read the input file, saving the reads that belong to one of the catalog loci.
    ReadsByCLoc(
            Input* pe_reads_f,
            size_t n_cloci,
            const unordered_map<string, size_t> read_name_to_cloc
            )
            : unique_seqs()
            , readsets(n_cloci)
            , n_used_reads(0) {
        Seq seq;
        while(pe_reads_f->next_seq(seq)) {
            auto cloc_it = read_name_to_cloc.find(seq.id);
            if (cloc_it != read_name_to_cloc.end()) {
                ++n_used_reads;
                add_seq_to_cloc(cloc_it->second, seq);
            }
        }
    }

    // add_seq_to_cloc(): Adds a read to the given clocus.
    void add_seq_to_cloc(size_t cloc, const Seq& seq) {
        const DNANSeq* seq_ptr = *unique_seqs.insert(new DNANSeq(seq.seq)).first;
        vector<Seq>& stack = readsets.at(cloc)[seq_ptr]; // First call creates the element.
        stack.push_back(seq);
        stack.back().drop_seq();
    }

    void clear() {
        readsets.clear();
        for (auto seq_ptr : unique_seqs)
            delete seq_ptr;
        unique_seqs.clear();
    }
};

/* Function declarations.
 * ========== */

// link_reads_to_cloci()
// ----------
// Parse the matches and tags (and fastq) files to link the paired reads to catalog loci.
// Uses globals `prefix_path`, `first_reads_path` and `second_reads_path`.
// Returns the number of reads in the fastq file.
void link_reads_to_cloci(unordered_map<string, size_t>& read_name_to_cloc, vector<int>& cloc_to_cloc_id);

/* main()
 * ========== */
int main(int argc, char** argv) {

    // Fix options
    prefix_path = "./s13_an_01";
    first_reads_path = "./s13_an_01.1.fq";
    paired_reads_path = "./s13_an_01.2.fq";
    paired_alns_path = "./s13_an_01.2.bam";
    in_file_type = FileT::bam;

    /*
     * Parse the matches, tags and fastq files to assign catalog loci to the reads.
     * ----------
     */
    cerr << "Reading read-to-locus information from the matches and tags files...\n";
    vector<int> cloc_to_cloc_id;
    unordered_map<string, size_t> read_name_to_cloc; // [(read index, cloc index)]
    link_reads_to_cloci(read_name_to_cloc, cloc_to_cloc_id);

    /*
     * Loading the paired-ends.
     * ----------
     */
    cerr << "Loading the paired-end sequences...\n";
    Input* pe_reads_f;
    if (in_file_type == FileT::bam)
        pe_reads_f = new Bam(paired_alns_path.c_str());
    const ReadsByCLoc reads_by_cloc (pe_reads_f, cloc_to_cloc_id.size(), read_name_to_cloc);
    delete pe_reads_f;
    cerr << "Used " << reads_by_cloc.n_used_reads << " aligned paired-end reads.\n";

    //todo {
    /*cerr << "Calculating the average coverage...\n";
    double avgcov = 0;
    for (const auto& c : contigs) {
        double c_avgcov = 0;
        for (const auto& col : c.cols)
            c_avgcov += col.coverage();
        c_avgcov /= c.cols.size();
        avgcov += c_avgcov;
    }
    avgcov /= contigs.size();
    cerr << "Average coverage is " << avgcov << "\n";*/
    //todo }

    // Convert the contigs to locus objects.
    map<int, MergedStack*> loci;

    return 0;
}

void link_reads_to_cloci(unordered_map<string, size_t>& pread_name_to_cloc, vector<int>& cloc_to_cloc_id) {

    // Load the matches.
    // Look for sample & catalog loci in a bijective relationship.
    // Assign indexes to cloci.
    unordered_map<int, size_t> sloc_id_to_cloc;

    vector<CatMatch*> matches;
    load_catalog_matches(prefix_path, matches);

    unordered_map<int, set<int> > cloc_id_to_sloc_ids;
    unordered_map<int, set<int> > sloc_id_to_cloc_ids;
    for (const CatMatch* m : matches) {
        cloc_id_to_sloc_ids[m->cat_id].insert(m->tag_id);
        sloc_id_to_cloc_ids[m->tag_id].insert(m->cat_id);
    }
    for (const auto& sloc : sloc_id_to_cloc_ids) {
        if (sloc.second.size() == 1) {
            const int cloc_id = *sloc.second.begin();
            if (cloc_id_to_sloc_ids.at(*sloc.second.begin()).size() == 1) {
                // Bijective, keep it.
                sloc_id_to_cloc.insert({sloc.first, cloc_to_cloc_id.size()});
                cloc_to_cloc_id.push_back(cloc_id);
            }
        }
    }

    cloc_id_to_sloc_ids.clear();
    sloc_id_to_cloc_ids.clear();
    for (const CatMatch* m : matches)
        delete m;
    matches.clear();

    // Read the tags file.
    // Link the first reads to cloci.
    unordered_map<string, size_t> fread_name_to_cloc;

    map<int, Locus*> sloci;
    bool tmp;
    if(load_loci(prefix_path, sloci, true, false, tmp) != 1) {
        cerr << "Error: could not find stacks files '" << prefix_path << ".*' (tags, snps and/or alleles).\n";
        throw exception();
    }

    for (const auto& element : sloci) {
        const Locus& sloc = *element.second;
        for (const char* fread_name : sloc.comp) {
            fread_name_to_cloc.insert({string(fread_name), sloc_id_to_cloc.at(sloc.id)});
        }
    }
    for (const auto& element : sloci)
        delete element.second;
    sloci.clear();
    sloc_id_to_cloc.clear();

    //xxx For now, we work with the first reads only...
    swap(pread_name_to_cloc, fread_name_to_cloc);

    // Attempt to guess the format of the read names;
    // otherwise read the fastq file indexes to link the first & paired names.
    // ...
}
