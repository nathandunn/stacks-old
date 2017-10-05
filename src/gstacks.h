// -*-mode:c++; c-style:k&r; c-basic-offset:4;-*-
//
// Copyright 2017, Julian Catchen <jcatchen@illinois.edu>
//
// This file is part of Stacks.
//
// Stacks is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// Stacks is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with Stacks.  If not, see <http://www.gnu.org/licenses/>.
//
#ifndef RYSTACKS_H
#define RYSTACKS_H

#include "constants.h"
#include "nucleotides.h"
#include "models.h"
#include "SuffixTree.h"
#include "GappedAln.h"

void parse_command_line(int argc, char* argv[]);
void report_options(ostream& os);

enum class GStacksInputT {unknown, denovo, refbased};

 //
// PhasedHet & PhaseSet
// ----------
//
struct PhasedHet {
    size_t phase_set;//TODO // N.B. The convention in VCF is to use the column of the first phased SNP for this.
    Nt4 left_allele;
    Nt4 right_allele;

    bool operator== (const PhasedHet& other) const
        {return phase_set == other.phase_set && left_allele == other.left_allele && right_allele == other.right_allele;}
};
class PhaseSet {
    vector<PhasedHet> phased_hets_;
public:
    PhaseSet() : phased_hets_{} {}
    void clear() {phased_hets_.clear();}
    bool empty() const {return phased_hets_.empty();}

    PhaseSet(size_t n_hets) : phased_hets_(n_hets, {SIZE_MAX, Nt4::n, Nt4::n}) {}
    size_t size() const {return phased_hets_.size();}
    const PhasedHet& het(size_t het_i) const {return phased_hets_[het_i];}

    // Add the first node to the phase set.
    void add_het(size_t het_i, array<Nt2,2> alleles);
    // Add a singleton node to the phase set, according to the given edge.
    void add_het(size_t het_i, array<Nt2,2> alleles, Nt2 nt_i, size_t het_j, Nt2 nt_j);
    // Merge with another phase set, according to the given edge.
    void merge_with(const PhaseSet& other, size_t het_i, Nt2 nt_i, size_t het_j, Nt2 nt_j);
    // Check that a new edge between two nodes within the phase set is consistent.
    bool is_edge_consistent(size_t het_i, Nt2 nt_i, size_t het_j, Nt2 nt_j) const;
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

    const size_t& at(size_t snp_i1, Nt2 snp1_allele, size_t snp_i2, Nt2 snp2_allele) const;
          size_t& at(size_t snp_i1, Nt2 snp1_allele, size_t snp_i2, Nt2 snp2_allele)
              {return (size_t&)((const SnpAlleleCooccurrenceCounter&)*this).at(snp_i1, snp1_allele, snp_i2, snp2_allele);}
};

//
// ProcessingStats
// ----------
// Statistics produced by `LocusProcessor::process()`.
//
class ProcessingStats {
public:
    size_t n_nonempty_loci;
    size_t n_loci_w_pe_reads;
    size_t n_loci_almost_no_pe_reads;
    size_t n_loci_pe_graph_not_dag;
    size_t n_aln_reads;
    size_t n_tot_reads;
    size_t n_se_pe_loc_overlaps;
    size_t mean_se_pe_loc_overlap;

    map<pair<size_t,size_t>,size_t> n_badly_phased_samples; // { {n_bad_samples, n_tot_samples} : count }

    size_t n_loci_phasing_issues()  const { size_t n = 0; for (auto& e: n_badly_phased_samples) n+=e.second; return n; }
    size_t n_loci_no_pe_reads()     const { return n_nonempty_loci - n_loci_w_pe_reads; }
    size_t n_loci_usable_pe_reads() const { return n_loci_w_pe_reads - n_loci_almost_no_pe_reads - n_loci_pe_graph_not_dag; }

    ProcessingStats& operator+= (const ProcessingStats& other);
};

//
// LocData
// ----------
// Data regarding the locus LocusProcessor currently works on.
//
class LocData {
public:
    int id;
    PhyLoc pos;
    const MetaPopInfo* mpopi;

    bool overlapped;

    string o_vcf;
    string o_fa;
    string o_details;
    stringstream details_ss;
    
    void clear();
};

//
// LocusProcessor
// ----------
// Functor for processing loci. Thread-specific.
//
class LocusProcessor {
public:
    LocusProcessor() : stats_(), loc_() {}

    // Process a locus.
    void process(CLocReadSet& loc);
    void process(CLocAlnSet& aln_loc);

    // Access the output. Statistics & movable fasta/vcf per-locus text outputs.
    const ProcessingStats& stats() const {return stats_;}
    string& vcf_out() {return loc_.o_vcf;}
    string& fasta_out() {return loc_.o_fa;}
    string& details_out() {return loc_.o_details;}

private:
    ProcessingStats stats_;
    mutable LocData loc_;

    string assemble_contig(const vector<const DNASeq4*>& seqs);
    bool add_read_to_aln(
            CLocAlnSet& aln_loc,
            AlignRes& aln_res,
            SRead&& read,
            GappedAln* aligner,
            SuffixTree* stree
            ) const;

    int align_reads_to_contig(SuffixTree *st, GappedAln *g_aln, DNASeq4 query, AlignRes &aln_res) const;
    int find_locus_overlap(SuffixTree *st, DNASeq4 se_consensus) const;

    DNASeq4 build_consensus(
            const vector<SiteCounts>& depths,
            const vector<SiteCall>& calls
            ) const;

    // For each sample, phase heterozygous SNPs.
    vector<map<size_t,PhasedHet>> phase_hets (
            const vector<SiteCall>& calls,
            const CLocAlnSet& aln_loc,
            set<size_t>& inconsistent_samples
            ) const;

    void count_pairwise_cooccurrences(
            SnpAlleleCooccurrenceCounter& cooccurrences,
            const CLocAlnSet& aln_loc,
            size_t sample,
            const vector<size_t>& snp_cols,
            const vector<size_t>& het_snps,
            const vector<const SampleCall*>& sample_het_calls
            ) const;

    // Assemble phase sets.
    // This is based on the graph of cooccurrences between the het alleles.
    // We build a graph in which the nodes are phased pairs of hets alleles
    // (`PhasedHet`s) and phase sets subgraphs of such nodes. The edges between
    // pairs of alleles are based on those in the graph of cooccurrences, but
    // in a phased context they appear as one of two types:
    // * parallel (+): left-allele-to-right-allele (left-to-left) or right-to-right.
    // * crossed (-): left-to-right, or right-to-left.
    // Allele linkage is transitive (through composition of +/- edges).
    // The assembly fails if there are inconsistent edges, i.e. edges that link
    // (transitively) one allele of one het SNP to both alleles of another het
    // SNP.
    // In practice, we start with an empty graph and iterate on edges, adding
    // nodes/PhasedHets when we first see them and checking that new edges are
    // consistent with the growing graph.
    bool assemble_phase_sets(
            vector<PhaseSet>& phase_sets,
            const vector<size_t>& het_snps,
            const vector<const SampleCall*>& sample_het_calls,
            const SnpAlleleCooccurrenceCounter& cooccurrences
            ) const;

    // Create the fasta/vcf text outputs.
    void write_one_locus (
            const CLocAlnSet& aln_loc,
            const vector<SiteCounts>& depths,
            const vector<SiteCall>& calls,
            const vector<map<size_t,PhasedHet>>& phase_data // {col : phasedhet} maps, for all samples
    );

    // (debug) Write a sample's haplotype graph.
    void write_sample_hapgraph(
            ostream& os,
            size_t sample,
            const vector<size_t>& het_snps,
            const vector<size_t>& snp_cols,
            const vector<const SampleCall*>& sample_het_calls,
            const SnpAlleleCooccurrenceCounter& cooccurences
            ) const;

    //
    // (debug) Creates a CLocAlnSet from a CLocReadSet, using the true
    // reference but computing paired-end reads alignments. Same requirements as
    // `from_true_alignments()`.
    //
    void using_true_reference(CLocAlnSet& aln_loc, CLocReadSet&& loc);

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
void from_true_alignments(CLocAlnSet& aln_loc, CLocReadSet&& loc, bool merge_reads);

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
    double writing_details;
    double sum() const {return clocking + reading + processing + writing_fa + writing_vcf + writing_details;}

    Clocks& operator+= (const Clocks& other);
    Clocks& operator/= (double d);
    friend Clocks operator+ (const Clocks& lhs, const Clocks& rhs) {Clocks sum (lhs); sum+=rhs; return sum;}
};

#endif
