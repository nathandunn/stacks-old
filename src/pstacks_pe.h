#include <ostream>
#include <string>
#include <vector>
#include <unordered_set>
#include <unordered_map>

void parse_command_line(int argc, char* argv[]);
void report_options(std::ostream& fh);

// retrieve_bijective_sloci()
// ----------
// Parse the matches file and return the ids of the sample loci that are in a
// bijective relationship with the catalog.
// Uses global `prefix_path`.
std::unordered_set<int> retrieve_bijective_sloci();

// convert_fw_read_name_to_paired()
// ----------
// Given a forward read name, guess the paired-end read name. The forward read
// name is expected to end in '/1' or '_1'.
void convert_fw_read_name_to_paired(std::string& read_name);

// link_reads_to_cloci()
// ----------
// Parses the matches and tags (and fastq) files to link the paired reads to catalog loci.
// Also, sets `gzipped_input` to the appropriate value.
// Uses globals `prefix_path` and `paired_alns_path`.
void link_reads_to_loci(
        std::unordered_map<std::string, size_t>& read_name_to_loc,
        std::vector<int>& sloc_ids,
        bool& gzipped_input
        );

// CLocReadSet
// ----------
struct ReadsByCLoc {
public:
    size_t n_used_reads;

    // Read the input file, saving the reads that belong to one of the catalog loci.
    ReadsByCLoc(Input* pe_reads_f,
                size_t n_cloci,
                const std::unordered_map<string, size_t>& read_name_to_cloc);

    // Obtain the MergedStack's and PStack's.
    // (This progressively clears the CLocReadSet's.)
    void convert_to_pmstacks(const vector<int>& sloc_to_sloc_id,
                             map<int, PStack*>& pstacks,
                             map<int, MergedStack*>& mstacks
                             );

private:
    // For each locus, we group reads by sequence.
    typedef std::map<const DNANSeq*, std::vector<Seq> > CLocReadSet;

    std::unordered_set<const DNANSeq*> unique_seqs;
    std::vector<CLocReadSet> readsets;

    // Add a read to the given clocus.
    void add_seq_to_cloc(std::size_t cloc, Seq& seq) {
        const DNANSeq* key = new DNANSeq(seq.seq);
        auto insertion = unique_seqs.insert(key);
        if (!insertion.second)
            delete key;
        key = *insertion.first;
        std::vector<Seq>& stack = readsets.at(cloc)[key]; // First call constructs the vector<Seq>.
        seq.delete_seq(); // Now stored in `unique_seqs`.
        stack.push_back(seq);
    }
};
