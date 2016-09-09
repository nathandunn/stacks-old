#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <map>
#include <unordered_set>
#include <unordered_map>
#include <cassert>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "constants.h"
#include "log_utils.h"
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
string paired_alns_path;
LogAlterator* log_alt = NULL;

/* CLocReadSet
 * =========== */
struct ReadsByCLoc {
public:
    size_t n_used_reads;

    // Read the input file, saving the reads that belong to one of the catalog loci.
    ReadsByCLoc(Input* pe_reads_f,
                size_t n_cloci,
                const unordered_map<string, size_t>& read_name_to_cloc);

    // Obtain the MergedStack's and PStack's.
    void convert_to_pmstacks(const vector<int>& cloc_to_cloc_id,
                             map<int, PStack*>& pstacks,
                             map<int, MergedStack*>& mstacks
                             );

private:
    // For each locus, we group reads by sequence.
    typedef map<const DNANSeq*, vector<Seq> > CLocReadSet;

    unordered_set<const DNANSeq*> unique_seqs;
    vector<CLocReadSet> readsets;

    // Add a read to the given clocus.
    void add_seq_to_cloc(size_t cloc, const Seq& seq) {
        const DNANSeq* seq_ptr = *unique_seqs.insert(new DNANSeq(seq.seq)).first;
        vector<Seq>& stack = readsets.at(cloc)[seq_ptr]; // First call creates the element.
        stack.push_back(seq);
        stack.back().drop_seq();
    }
};

ReadsByCLoc::ReadsByCLoc(
        Input* pe_reads_f,
        size_t n_cloci,
        const unordered_map<string, size_t>& read_name_to_cloc
        )
: n_used_reads(0)
, unique_seqs()
, readsets(n_cloci)
{
    Seq seq;
    while(pe_reads_f->next_seq(seq)) {
        auto cloc_it = read_name_to_cloc.find(seq.id);
        if (cloc_it != read_name_to_cloc.end()) {
            ++n_used_reads;
            add_seq_to_cloc(cloc_it->second, seq);
        }
    }
}

void ReadsByCLoc::convert_to_pmstacks(
        const vector<int>& cloc_to_cloc_id,
        map<int, PStack*>& pstacks,
        map<int, MergedStack*>& mstacks
        ) {
    int pstack_id = 0;

    for (size_t cloc=0; cloc<readsets.size(); ++cloc) {

        //
        // First, obtain the raw PStacks of this c-locus.
        // The PStacks share their sequence and location.
        //
        vector<PStack*> cloc_pstacks;
        for (auto& identical_reads : readsets[cloc]) {
            map<PhyLoc, PStack*> cloc_pstacks_by_loc;

            for (const Seq& s : identical_reads.second) {
                PStack*& p = cloc_pstacks_by_loc.insert({s.loc, NULL}).first->second; // n.b. ref to pointer
                if (p == NULL) {
                    // This element was just created; this is a novel (location, seq) combination.
                    p = new PStack();
                    p->id = pstack_id;
                    ++pstack_id;
                    p->loc = s.loc;
                    p->add_seq(identical_reads.first);

                    cloc_pstacks.push_back(p);
                }
                ++p->count;
                p->add_id(s.id);
            }
        }

        //
        // Determine the positions spanned by the PStacks.
        //
        assert(!cloc_pstacks.empty()); // The sample has "matches".
        const PStack* first_pstack = *cloc_pstacks.begin();
        PhyLoc loc = first_pstack->loc;
        size_t len = first_pstack->seq->size();
        set<const PStack*> non_olap; //xxx For now, ignore non-overlapping pstacks.
        for (const PStack* p : cloc_pstacks) {
            // Check that the PStack is on the same chromosome and strand.
            if (strcmp(p->loc.chr, loc.chr) != 0
                    || p->loc.strand != loc.strand) {
                non_olap.insert(p);
                continue;
            }

            if (loc.strand == strand_plus) {
                // Check that the PStack overlaps.
                if (p->loc.bp > loc.bp + len - 1 || p->loc.bp +p->seq->size() - 1 < loc.bp) {
                    non_olap.insert(p);
                    continue;
                }
                if (loc.bp > p->loc.bp) {
                    // Extend left.
                    len += p->loc.bp - loc.bp;
                    loc.bp = p->loc.bp;
                }
                if (loc.bp + len < p->loc.bp +p->seq->size()) {
                    // Extend right.
                    len += p->loc.bp +p->seq->size() - loc.bp + len;
                }
            } else {
                // loc.strand == strand_minus

                // Check that the PStack overlaps.
                if (p->loc.bp -p->seq->size() + 1 > loc.bp || p->loc.bp < loc.bp - len + 1) {
                    non_olap.insert(p);
                    continue;
                }
                if (p->loc.bp -p->seq->size() > p->loc.bp -p->seq->size()) {
                    // Extend left.
                    len += p->loc.bp -p->seq->size() - p->loc.bp -p->seq->size();
                }
                if (loc.bp < p->loc.bp) {
                    // Extend right.
                    len += p->loc.bp - loc.bp;
                    loc.bp = p->loc.bp;
                }
            }
        }

        //
        // Extend the PStacks so that they all have the same location and length.
        //
        for (PStack* p : cloc_pstacks) {
            if (non_olap.count(p))
                continue;
            p->extend(loc, len);
        }

        //
        // Create the MergedStack of the clocus
        //
        MergedStack* m = new MergedStack();
        m->id = cloc_to_cloc_id[cloc];
        m->add_consensus(first_pstack->seq);
        assert(m->len == len);
        m->loc = loc;
        for (const PStack* p : cloc_pstacks) {
            if (non_olap.count(p))
                continue;
            m->count += p->count;
            m->utags.push_back(p->id);
        }

        //
        // Insert the MergedStack and the PStack's of the c-locus in
        // in the global objects.
        //
        mstacks.insert({m->id, m});
        for (PStack* p : cloc_pstacks) {
            if (non_olap.count(p))
                continue;
            pstacks.insert({p->id, p});
        }

    } // for(c-locus)

    return;
}

/* link_reads_to_cloci()
 * ==========
 * Parses the matches and tags (and fastq) files to link the paired reads to catalog loci.
 * Also, sets `gzipped_input` to the appropriate value.
 * Uses globals `prefix_path`, `first_reads_path` and `second_reads_path`.
 */
void link_reads_to_cloci(unordered_map<string, size_t>& read_name_to_cloc, vector<int>& cloc_to_cloc_id, bool& gzipped_input);

/* main()
 * ========== */
int main(int argc, char** argv) {
#ifndef DEBUG
try {
#endif

    // Fix options
    prefix_path = "./s13_an_01";
    sql_id        = 1;
    //paired_alns_path = "/projects/catchenlab/rochette/sbk/scan/samples/s13_an_01.bam";
    paired_alns_path = "/home/rochette/src/stacks/n_tmp/20160908.1sample_red/s13_an_01.bijective.bam";
    in_file_type = FileT::bam;

    string log_path = prefix_path + ".pstacks_pe.log";
    log_alt = new LogAlterator(ofstream(log_path));

#ifdef _OPENMP
    omp_set_num_threads(num_threads);
#endif

    /*
     * Parse the matches, tags and fastq files to assign catalog loci to the reads.
     * ----------
     */
    cerr << "Reading read-to-locus information from the matches and tags files...\n";
    vector<int> cloc_to_cloc_id;
    unordered_map<string, size_t> read_name_to_cloc; // [(read index, cloc index)]
    bool is_input_gzipped;
    link_reads_to_cloci(read_name_to_cloc, cloc_to_cloc_id, is_input_gzipped);
    cerr << "Found " << read_name_to_cloc.size() << " paired-end reads spanning "
           << cloc_to_cloc_id.size() << " catalog loci.\n";

    /*
     * Load the paired-ends.
     * ----------
     */
    cerr << "Loading the paired-end sequences...\n";
    Input* pe_reads_f;
    if (in_file_type == FileT::bam)
        pe_reads_f = new Bam(paired_alns_path.c_str());
    ReadsByCLoc reads_by_cloc (pe_reads_f, cloc_to_cloc_id.size(), read_name_to_cloc);
    delete pe_reads_f;
    read_name_to_cloc.clear();
    cerr << "Used " << reads_by_cloc.n_used_reads << " aligned paired-end reads.\n"; //xxx "used"

    /*
     * Convert the data to PStack's and MStack's.
     */
    cerr << "Assembling the paired-end sequences...\n";
    map<int, MergedStack*> loci;
    map<int, PStack*> stacks;
    reads_by_cloc.convert_to_pmstacks(cloc_to_cloc_id, stacks, loci);
    cerr << "Created " << loci.size() << " loci, made of " << stacks.size() << " stacks.\n";

    /*
     * Call SNPs and alleles.
     */
    call_consensus(loci, stacks, true);

    write_results(loci, stacks, is_input_gzipped, true);

    return 0;

#ifndef DEBUG
} catch (const exception& e) {
    cerr << "Terminated after an error occurred (" << e.what() << ").\n";
    return -1;
}
#endif
}

void link_reads_to_cloci(unordered_map<string, size_t>& pread_name_to_cloc, vector<int>& cloc_to_cloc_id, bool& is_input_gzipped) {

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
    map<int, Locus*> sloci;
    if(load_loci(prefix_path, sloci, true, false, is_input_gzipped) != 1) {
        cerr << "Error: could not find stacks files '" << prefix_path << ".*' (tags, snps and/or alleles).\n";
        throw exception();
    }

    // For each first read, guess the name of the paired-end
    // read (names are expected to end in '/1' or '_1'),
    // and link them to cloci.
    for (const auto& element : sloci) {
        const Locus& sloc = *element.second;
        for (const char* fread_name : sloc.comp) {
            string pread_name (fread_name);
            if (pread_name.length() < 2
                    || (pread_name.substr(pread_name.length()-2) != "/1"
                            && pread_name.substr(pread_name.length()-2) != "_1")
                    ) {
                cerr << "Error: Unrecognized read name format; expected '"
                     << pread_name << "' to end with '/1' or '_1'.\n";
                throw exception();
            }
            //pread_name.at(pread_name.length()-1) = '2';
            pread_name_to_cloc.insert({pread_name, sloc_id_to_cloc.at(sloc.id)});
        }
    }

    // Delete what `load_loci()` allocated.
    for (const auto& element : sloci)
        delete element.second;

    //xxx Save list of reads.
    /*ofstream readnames_f ("readnames");
    for (auto& r : pread_name_to_cloc)
        readnames_f << r.first << "\n";
    */

    // Attempt to guess the format of the read names;
    // otherwise read the fastq file indexes to link the first & paired names.
    // ...

}
