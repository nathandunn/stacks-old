#ifndef RYSTACKS_H
#define RYSTACKS_H

#include "constants.h"
#include "nucleotides.h"
#include "models.h"
#include "SuffixTree.h"
#include "GappedAln.h"

void parse_command_line(int argc, char* argv[]);
void report_options(ostream& os);

//
// PhasedHet
// ----------
// Phasing data in a phase-set-based format so that it is eaasy to write a
// phased VCF.
//
struct PhasedHet {
    size_t phase_set; // N.B. The convention in VCF is to use the column of the first phased SNP for this.
    Nt2 left_allele;
    Nt2 right_allele;
};

//
// SnpAlleleCooccurrenceCounter
// ----------
// Convenience matrix class to record SNP allele co-occurences in
// reads/read pairs.
//
// The matrix is n_snps*n_snps but we only use the i<j half.
//
class SnpAlleleCooccurrenceCounter {
    size_t n_snps_;
    vector<array<array<size_t,4>,4>> cooccurences_;
public:
    SnpAlleleCooccurrenceCounter(size_t n_snps) : n_snps_(n_snps), cooccurences_(n_snps_*n_snps_) {};
    void clear();

    size_t& at(size_t snp_i1, Nt2 snp1_allele, size_t snp_i2, Nt2 snp2_allele);
};

//
// ProcessingStats
// ----------
// Statistics produced by `LocusProcessor::process()`.
//
struct ProcessingStats {
    size_t n_nonempty_loci;
    size_t n_loci_w_pe_reads;
    size_t n_loci_almost_no_pe_reads;
    size_t n_loci_pe_graph_not_dag;

    map<pair<size_t,size_t>,size_t> n_badly_phased_samples; // { {n_bad_samples, n_tot_samples} : count }

    size_t n_loci_phasing_issues()  const {size_t n=0; for(auto& e: n_badly_phased_samples) n+=e.second; return n;}
    size_t n_loci_no_pe_reads()     const {return n_nonempty_loci - n_loci_w_pe_reads;}
    size_t n_loci_usable_pe_reads() const {return n_loci_w_pe_reads - n_loci_almost_no_pe_reads - n_loci_pe_graph_not_dag;}

    ProcessingStats& operator+= (const ProcessingStats& other);
};

//
// LocusProcessor
// ----------
// Functor for processing loci. Thread-specific.
//
class LocusProcessor {
public:
    LocusProcessor() : stats_(), loc_id_(-1), loc_pos_(), mpopi_(NULL), o_vcf_(), o_fa_() {}

    // Process a locus.
    void process(CLocReadSet&& loc);

    // Access the output. Statistics & movable fasta/vcf per-locus text outputs.
    const ProcessingStats& stats() const {return stats_;}
    string& vcf_out() {return o_vcf_;}
    string& fasta_out() {return o_fa_;}

private:
    ProcessingStats stats_;
    int loc_id_;
    PhyLoc loc_pos_;
    const MetaPopInfo* mpopi_;
    string o_vcf_;
    string o_fa_;

    string assemble_contig(const vector<const DNASeq4*>& seqs);

    int align_reads_to_contig(SuffixTree *st, GappedAln *g_aln, DNASeq4 query, AlignRes &aln_res);

    // For each sample, phase heterozygous SNPs.
    vector<map<size_t,PhasedHet>> phase_hets (
            const vector<SiteCall>& calls,
            const CLocAlnSet& aln_loc,
            set<size_t>& inconsistent_samples
    ) const;

    // Prune the phasing output to keep only the best/largest phase set.
    void rm_supernumerary_phase_sets (vector<map<size_t,PhasedHet>>& phase_data) const;

    // Create the fasta/vcf text outputs.
    void write_one_locus (
            const CLocAlnSet& aln_loc,
            const vector<SiteCounts>& depths,
            const vector<SiteCall>& calls,
            const vector<map<size_t,PhasedHet>>& phase_data // {col : phasedhet} maps, for all samples
    );
};

//
// For debugging and benchmarking.
// ==========
//

//
// dbg_extract_cigar
// ----------
// Extracts a true alignment from a read ID. The ID is expected to contain
// `cig1=<CIGAR1>` and possibly `cig2=<CIGAR2>` fields, e.g.:
// `t00800n:msp_00:a2:r0:cig1=77M1I72M851H:cig2=211H150M639H/1`
//
Cigar dbg_extract_cigar(const string& read_id);

//
// from_true_alignments
// ----------
// Creates a CLocAlnSet from a CLocReadSet, using true alignments.
//
void from_true_alignments(CLocAlnSet& aln_loc, CLocReadSet&& loc);

//
// Clocks
// ----------
// For Benchmarking. Thread specific.
//
struct Clocks {
    double clocking;
    double reading;
    double processing;
    double writing_fa;
    double writing_vcf;
    double sum() const {return clocking+reading+processing+writing_fa+writing_vcf;}

    Clocks& operator+= (const Clocks& other);
    Clocks& operator/= (double d);
    friend Clocks operator+ (const Clocks& lhs, const Clocks& rhs) {Clocks sum (lhs); sum+=rhs; return sum;}
};

#endif
