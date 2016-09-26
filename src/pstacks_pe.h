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
// Uses glob `prefix_path`.
std::unordered_set<int> retrieve_bijective_sloci();

// convert_fw_read_name_to_paired()
// ----------
// Given a forward read name, guess the paired-end read name. The forward read
// name is expected to end in '/1' or '_1'.
void convert_fw_read_name_to_paired(std::string& read_name);

// link_reads_to_cloci()
// ----------
// Extracts the reads composing each locus from the tags file, and guesses the
// names of the paired-end reads (via `convert_fw_read_name_to_paired`).
// Uses glob `prefix_path`.
void link_reads_to_loci(
        const std::unordered_set<int>& bij_sloci,
        std::unordered_map<std::string, size_t>& read_name_to_loc,
        std::vector<int>& sloc_ids,
        bool& gzipped_input
        );

// load_aligned_reads()
// ----------
// Collapses reads into PStacks according to their sequence and PhyLoc.
// Uses glob `paired_alns_path`.
vector<vector<PStack> > load_aligned_reads(
        size_t n_loci,
        const std::unordered_map<std::string, size_t>& read_name_to_loc
        );

// merge_pstacks()
// ----------
// Creates a MergedStack from a set of PStack's. The PStacks are "extended" to
// the length of the contig.
MergedStack merge_pstacks(std::vector<PStack>& pstacks, int loc_id);
