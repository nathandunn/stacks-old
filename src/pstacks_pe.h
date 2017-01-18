#include <ostream>
#include <string>
#include <vector>
#include <unordered_set>
#include <unordered_map>

#include "stacks.h"

void parse_command_line(int argc, char* argv[]);
void report_options(std::ostream& fh);

// FwLocInfo
// ----------
// Keeps track of the IDs and genomic positions of the forward loci.
struct FwLocInfo {
    int id;
    PhyLoc loc;

    FwLocInfo(int i, const PhyLoc& l) : id(i), loc(l) {}

    bool is_upstream_of(const Seq& s) const;
};

// Contig
// ----------
// Keeps track of the positions spanned by a set of PStacks.
struct Contig {
    PhyLoc loc;
    int len;
    std::set<const PStack*> stacks;

    Contig(const PStack& s) : loc(s.loc), len(s.seq->size()), stacks({&s}) {}

    bool overlaps(const PhyLoc& other_loc, int other_len) const;
    // broaden(): same return value as overlaps(); updates `loc` and `len`.
    bool broaden(const PhyLoc& other_loc, int other_len);
    // add(): same return value as broaden(); updates `stacks`.
    bool add(const PStack& s);
    bool add(const Contig& s);
    // count(): number of reads in the contig.
    size_t count() const;
};

// link_reads_to_cloci()
// ----------
// Extracts the reads composing each locus from the tags file, and guesses the
// names of the paired-end reads (via `convert_fw_read_name_to_paired`).
// Uses glob `prefix_path`.
void link_reads_to_loci(
        std::vector<FwLocInfo>& loc_info,
        std::unordered_map<std::string, size_t>& read_name_to_loc
        );

// load_aligned_reads()
// ----------
// Collapses reads into PStacks according to their sequence and PhyLoc.
// Uses glob `paired_alns_path`.
std::vector<std::vector<PStack> > load_aligned_reads(
        const std::vector<FwLocInfo>& fwloci_info,
        const std::unordered_map<std::string, size_t>& read_name_to_loc
        );

// merge_pstacks()
// ----------
// Creates a MergedStack from a set of PStack's. The PStacks are "extended" to
// the length of the contig.
MergedStack merge_pstacks(std::vector<PStack>& pstacks, int loc_id);

//
// Inline definitions.
// ----------
//

inline
bool FwLocInfo::is_upstream_of(const Seq& s) const {
    static const int max_allowed_distance = 5000;

    if (strcmp(loc.chr, s.loc.chr) != 0)
        return false;

    // note: It is unneeded to check for the read's strand.

    if (loc.strand == strand_plus)
        return s.loc.bp >= loc.bp
                && s.loc.bp < loc.bp + max_allowed_distance;
    else
        return s.loc.bp <= loc.bp
                && int(s.loc.bp) > int(loc.bp) - max_allowed_distance;
}


inline
bool Contig::overlaps(const PhyLoc& oloc, int olen) const {
    // Check chromosome & strand.
    if (strcmp(oloc.chr, loc.chr) != 0
            || oloc.strand != loc.strand)
        return false;

    // The given locus does not overlap if it starts after the end of the contig
    // or ends before its start.
    if (loc.strand == strand_plus)
        return ! (
                oloc.bp > loc.bp + len - 1
                || oloc.bp + olen - 1 < loc.bp
                );
    else
        return ! (
                oloc.bp - olen + 1 > loc.bp
                || oloc.bp < loc.bp - len + 1
                );
}

inline
bool Contig::broaden(const PhyLoc& oloc, int olen) {
    if (! overlaps(oloc, olen))
        return false;

    if (loc.strand == strand_plus) {
        if (loc.bp > oloc.bp) {
            // Extend left.
            len += loc.bp - oloc.bp;
            loc.bp = oloc.bp;
        }
        if (loc.bp + len < oloc.bp + olen) {
            // Extend right.
            len += oloc.bp + olen - (loc.bp + len);
        }
    } else {
        if (loc.bp < oloc.bp) {
            // Extend right.
            len += oloc.bp - loc.bp;
            loc.bp = oloc.bp;
        }
        if (int(loc.bp) - len > int(oloc.bp) - olen) {
            // Extend left.
            len += int(loc.bp) - len - (int(oloc.bp) - olen);
        }
    }

    return true;
}

inline
bool Contig::add(const PStack& s) {
    if (broaden(s.loc, s.seq->size())) {
        stacks.insert(&s);
        return true;
    } else {
        return false;
    }
}

inline
bool Contig::add(const Contig& other) {
    if (broaden(other.loc, other.len)) {
        stacks.insert(other.stacks.begin(), other.stacks.end());
        return true;
    } else {
        return false;
    }
}

inline
size_t Contig::count() const {
    size_t n_reads = 0;
    for (const PStack* s : stacks)
        n_reads += s->count;
    return n_reads;
}
