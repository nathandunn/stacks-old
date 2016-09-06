#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <deque>
#include <set>
#include <map>
#include <unordered_map>

#include "utils.h"
#include "models.h"
#include "input.h"
#include "BamI.h"
#include "stacks.h"
#include "sql_utilities.h"

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

/* AlnContig
 * ========== */
struct AlnContig {
    PhyLoc loc;
    deque<NtCounts> cols;

    size_t skipped_records; //xxx default 0

    void add(const Seq& alnseq);
    string consensus() const;

private:
    void init(const PhyLoc& l, size_t n) {loc=l; cols=deque<NtCounts>(n);}
};

string AlnContig::consensus() const {
    string s;
    s.reserve(cols.size());
    for (const NtCounts& c : cols)
        s.push_back(c.consensus());
    return s;
}

void AlnContig::add(const Seq& alnseq) {
    size_t seqlen = strlen(alnseq.seq);

    if (loc.chr == NULL)
        init(alnseq.loc, seqlen);

    //xxx If the record does not overlap the existing contig, skip it.
    if (strcmp(alnseq.loc.chr, loc.chr) != 0
            || alnseq.loc.strand != loc.strand
            || alnseq.loc.bp > loc.bp + cols.size() - 1
            || alnseq.loc.bp + seqlen -1 < loc.bp
            ) {
        ++skipped_records;
        return;
    }

    // Extend the contig, if necessary.
    while (loc.bp > alnseq.loc.bp) {
        --loc.bp;
        cols.push_front(NtCounts());
    }
    while (loc.bp + cols.size() < alnseq.loc.bp + seqlen)
        cols.push_back(NtCounts());

    // Record the nucleotides of the sequence.
    const char* nt = alnseq.seq;
    deque<NtCounts>::iterator col = cols.begin();
    for (size_t i=0; i < alnseq.loc.bp - loc.bp; ++i)
        ++col;
    for (size_t i=0; i<seqlen; ++i) {
        col->add(*nt);
        ++nt;
        ++col;
    }
}

/* Function declarations.
 * ========== */

// link_reads_to_cloci()
// ----------
// Parse the matches, tags and fastq files to link the paired reads to catalog loci.
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

    // Parse the matches, tags and fastq files to assign catalog loci to the reads.
    cerr << "Reading read-to-locus information from the matches and tags files...\n";
    vector<int> cloc_to_cloc_id;
    unordered_map<string, size_t> read_name_to_cloc; // [(read index, cloc index)]
    link_reads_to_cloci(read_name_to_cloc, cloc_to_cloc_id);

    // Read the paired-ends, build the contig objects.
    cerr << "Reading the alignments & building the contigs...\n";
    Input* pe_reads_f;
    if (in_file_type == FileT::bam)
        pe_reads_f = new Bam(paired_alns_path.c_str());
    vector<AlnContig> contigs = vector<AlnContig>(cloc_to_cloc_id.size());
    size_t n_preads_used;
    Seq alnseq;
    while(pe_reads_f->next_seq(alnseq)) {
        auto cloc = read_name_to_cloc.find(alnseq.id);
        if (cloc != read_name_to_cloc.end()) {
            ++n_preads_used;
            AlnContig& c = contigs[cloc->second];
            c.add(alnseq);
        }
    }
    cerr << "Used " << n_preads_used << " aligned paired-end reads.\n";

    //todo
    cerr << "Calculating the average coverage...\n";
    double avgcov = 0;
    for (const auto& c : contigs) {
        double c_avgcov = 0;
        for (const auto& col : c.cols)
            c_avgcov += col.coverage();
        c_avgcov /= c.cols.size();
        avgcov += c_avgcov;
    }
    avgcov /= contigs.size();
    cerr << "Average coverage is " << avgcov << "\n";

    // Convert the contigs to locus objects.
    map<int, MergedStack*> loci;

    return 0;
}

void link_reads_to_cloci(unordered_map<string, size_t>& pread_name_to_cloc, vector<int>& cloc_to_cloc_id) {

    // Load the matches.
    // Look for sample & catalog loci in a bijective relationship.
    // Assign indexes to cloci.
    map<int, size_t> sloc_id_to_cloc;

    vector<CatMatch*> matches;
    load_catalog_matches(prefix_path, matches);

    map<int, set<int> > cloc_id_to_sloc_ids;
    map<int, set<int> > sloc_id_to_cloc_ids;
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
