#include <getopt.h>

#include "zlib.h"

#include "constants.h"
#include "utils.h"
#include "log_utils.h"
#include "catalog_utils.h"
#include "locus.h"
#include "locus_readers.h"
#include "debruijn.h"
#include "GappedAln.h"
#include "aln_utils.h"
#include "Alignment.h"
#include "models.h"

using namespace std;

int stacks_handle_exceptions(const exception& e) {
    std::cerr << "Aborted.";
    if (typeid(e) != typeid(std::exception))
        std::cerr << " (" << e.what() << ")";
    std::cerr << "\n";
    return 13;
}

struct PhasedHet {
    size_t phase_set; // N.B. The convention in VCF is to use the column of the first phased SNP for this.
    Nt2 left_allele;
    Nt2 right_allele;
};

class SnpAlleleCooccurrenceCounter {
    size_t n_snps_;
    vector<array<array<size_t,4>,4>> cooccurences_;
public:
    SnpAlleleCooccurrenceCounter(size_t n_snps);
    size_t& at(size_t snp_i1, Nt2 snp1_allele, size_t snp_i2, Nt2 snp2_allele);
    void clear();
};

struct ProcessingStats {
    size_t n_nonempty_loci;
    size_t n_loci_w_pe_reads;
    size_t n_loci_almost_no_pe_reads;
    size_t n_loci_pe_graph_not_dag;

    size_t n_loci_phasing_issues() const;
    size_t n_loci_no_pe_reads() const {return n_nonempty_loci - n_loci_w_pe_reads;}
    size_t n_loci_usable_pe_reads() const {return n_loci_w_pe_reads - n_loci_almost_no_pe_reads - n_loci_pe_graph_not_dag;}

    map<pair<size_t,size_t>,size_t> n_badly_phased_samples; // { {n_bad_samples, n_tot_samples} : count }

    ProcessingStats& operator+= (const ProcessingStats& other);
};

struct ProcessingOutput {
    string vcf;
    string fa;

    void write(VcfWriter& vcf_f, gzFile gzfasta_f) const {
        vcf_f.file() << vcf;
        gzwrite(gzfasta_f, fa.c_str(), fa.length());
    }
};

class LocusProcessor {
public:
    LocusProcessor() : stats_() {}
    void process(CLocReadSet&& loc);

    const ProcessingStats& stats() const {return stats_;}
    const ProcessingOutput& out() const {return out_;}

private:
    ProcessingStats stats_;
    ProcessingOutput out_;

    vector<map<size_t,PhasedHet>> phase_hets (
            const vector<SiteCall>& calls,
            const CLocAlnSet& aln_loc,
            set<size_t>& inconsistent_samples
    ) const;
    void rm_supernumerary_phase_sets (
            vector<map<size_t,PhasedHet>>& phase_data
    ) const;
    void write_one_locus (
            const CLocAlnSet& aln_loc,
            const vector<SiteCounts>& depths,
            const vector<SiteCall>& calls,
            const vector<map<size_t,PhasedHet>>& phase_data // {col : phasedhet} maps, for all samples
    );
};

void parse_command_line(int argc, char* argv[]);
void report_options(ostream& os);

Cigar dbg_extract_cigar(const string& read_id);

//
// Argument globals.
//
bool quiet = false;
int num_threads = 1;
string in_dir;
int batch_id = -1;
set<int> locus_wl;
bool ignore_pe_reads = false;

modelt model_type = marukilow;
unique_ptr<const Model> model;

size_t km_length = 31;
size_t min_km_count = 2;

bool dbg_no_haplotypes = false;
bool dbg_write_gfa = false;
bool dbg_write_alns = false;
bool dbg_write_hapgraphs = false;
bool dbg_write_nt_depths = false;
bool dbg_true_alns = false;

//
// Additional globals.
//
const string prog_name = "rystacks";
unique_ptr<LogAlterator> logger;
gzFile o_gzfasta_f = NULL;
unique_ptr<VcfWriter> o_vcf_f;
ofstream o_aln_f;
ofstream o_hapgraphs_f;

int main(int argc, char** argv) {
try {

    //
    // Parse arguments.
    //
    parse_command_line(argc, argv);

    //
    // Open the log.
    //
    string lg_path = in_dir + "batch_" + to_string(batch_id) + "." + prog_name + ".log";
    logger.reset(new LogAlterator(lg_path, quiet, argc, argv));
    report_options(cout);
    cout << "\n" << flush;

    //
    // Initialize OPENMP.
    //
    #ifdef _OPENMP
    omp_set_num_threads(num_threads);
    #endif

    //
    // Open the BAM file and parse the header.
    //
    BamCLocReader bam_fh (in_dir + "batch_" + to_string(batch_id) + ".catalog.bam");

    //
    // Open the output files.
    //
    string o_gzfasta_path = in_dir + "batch_" + to_string(batch_id) + "." + prog_name + ".fa.gz";
    o_gzfasta_f = gzopen(o_gzfasta_path.c_str(), "wb");
    check_open(o_gzfasta_f, o_gzfasta_path);

    string o_vcf_path = in_dir + "batch_" + to_string(batch_id) + "." + prog_name + ".vcf.gz";
    VcfHeader vcf_header;
    vcf_header.add_std_meta();
    for(auto& s : bam_fh.mpopi().samples())
        vcf_header.add_sample(s.name);
    o_vcf_f.reset(new VcfWriter(o_vcf_path, move(vcf_header)));

    if (dbg_write_alns) {
        string o_aln_path = in_dir + "batch_" + to_string(batch_id) + "." + prog_name + ".alns";
        o_aln_f.open(o_aln_path);
        check_open(o_aln_f, o_aln_path);
        o_aln_f <<
            "# This prints observed read haplotypes:\n"
            "# show_loc() { loc=$1; sed -n \"/^END $loc\\b/ q; /^BEGIN $loc\\b/,$ p\" ${RY_DIR:-.}/batch_1.rystacks.alns | tail -n+2; }\n"
            "# snp_cols() { loc=$1; awk \"\\$1==$loc; \\$1>$loc {exit}\" ${RY_DIR:-.}/batch_1.rystacks.vcf | awk '$5!=\".\"' | cut -f2 | paste -sd ','; }\n"
            "# show_haps() { loc=$1; cols=$2; spl=$3; show_loc $loc | grep \"\\b$spl\\b\" | cut -f3 | cut -c \"$cols\" | sort; }\n"
            "# true_loci() { loc=$1; spl=$2; show_loc $loc | grep \"\\b$spl\\b\" | grep -v ref | cut -d: -f1 | sort -u; }\n"
            ;
    }

    if (dbg_write_hapgraphs) {
        string o_hapgraphs_path = in_dir + "batch_" + to_string(batch_id) + "." + prog_name + ".hapgraphs.dot";
        o_hapgraphs_f.open(o_hapgraphs_path);
        check_open(o_hapgraphs_f, o_hapgraphs_path);
        o_hapgraphs_f << "# dot -Tpdf -O batch_1.rystacks.hapgraphs.dot\n"
                      << "# loc=371\n"
                      << "# { g=batch_1.rystacks.hapgraphs.dot; sed -n '0,/^subgraph/p' $g | head -n-1; sed -n \"/^subgraph cluster_loc$loc\\b/,/^}/p\" $g; echo \\}; } | dot -Tpdf -o haps.$loc.pdf\n"
                      << "graph {\n"
                      << "edge[color=\"grey60\",fontsize=12,labeljust=\"l\"];\n";
    }

    //
    // Process every locus
    //
    ProcessingStats stats {};
    cout << "Processing " << (locus_wl.empty() ? "all" : "whitelisted") << " loci...\n" << flush;
    const size_t n_loci = locus_wl.empty() ? bam_fh.n_loci() : locus_wl.size();
    ProgressMeter progress (cout, n_loci);

    // For parallelization.
    int omp_return = 0;
    size_t next_to_write = 0;
    int next_id_to_write = bam_fh.target2id(next_to_write);
    if (!locus_wl.empty())
        // We don't use the `next_to_write` BAM index; instead we get the IDs from the whitelist.
        next_id_to_write = *locus_wl.begin();
    map<int,ProcessingOutput> outputs; // map of {loc_id, output}
    #pragma omp parallel
    {
        LocusProcessor loc_proc;
        CLocReadSet loc (bam_fh.mpopi());
        #pragma omp for schedule(dynamic)
        for (size_t i=0; i<n_loci; ++i) {
            if (omp_return != 0)
                continue;
            try {
                if (locus_wl.empty()) {
                    #pragma omp critical
                    bam_fh.read_one_locus(loc);
                } else {
                    #pragma omp critical
                    {
                        do {
                            if (!bam_fh.read_one_locus(loc)) {
                                cerr << "Error: Some whitelisted loci weren't found in the BAM file.\n";
                                omp_return = stacks_handle_exceptions(exception());
                                break;
                            }
                        } while (!locus_wl.count(loc.id()));
                    }
                    if (omp_return != 0)
                        continue;
                }

                int loc_id = loc.id();
                loc_proc.process(move(loc));

                #pragma omp critical
                {
                    if (loc_id == next_id_to_write) {
                        // Write it.
                        loc_proc.out().write(*o_vcf_f, o_gzfasta_f);
                        ++progress;
                        if (locus_wl.empty()) {
                            ++next_to_write;
                            if (next_to_write < bam_fh.n_loci())
                                next_id_to_write = bam_fh.target2id(next_to_write);
                            else
                                // We just wrote the last locus.
                                assert(outputs.empty()); // (Thus the next outputs.find() will fail.)
                        } else {
                            locus_wl.erase(loc_id);
                            if (!locus_wl.empty())
                                next_id_to_write = *locus_wl.begin();
                            else
                                assert(outputs.empty());
                        }
                        // Write stored output, if any.
                        map<int,ProcessingOutput>::iterator output;
                        while ((output = outputs.find(next_id_to_write)) != outputs.end()) {
                            output->second.write(*o_vcf_f, o_gzfasta_f);
                            outputs.erase(output);
                            ++progress;
                            if (locus_wl.empty()) {
                                ++next_to_write;
                                if (next_to_write < bam_fh.n_loci())
                                    next_id_to_write = bam_fh.target2id(next_to_write);
                                else
                                    // We just wrote the last locus.
                                    assert(outputs.empty());
                            } else {
                                locus_wl.erase(loc_id);
                                if (!locus_wl.empty())
                                    next_id_to_write = *locus_wl.begin();
                                else
                                    assert(outputs.empty());
                            }
                        }
                    } else {
                        // Store output for later.
                        outputs.insert( {loc_id, loc_proc.out()} );
                    }
                }

            } catch (exception& e) {
                #pragma omp critical
                omp_return = stacks_handle_exceptions(e);
            }
        }

        #pragma omp critical
        stats += loc_proc.stats();
    }
    if (omp_return != 0)
         return omp_return;
    progress.done();
    assert(locus_wl.empty());

    //
    // Report statistics on the analysis.
    //
    {
        size_t tot = stats.n_nonempty_loci;
        size_t ph = stats.n_loci_phasing_issues();
        auto pct = [tot](size_t n) {return as_percentage((double) n / tot) ;};
        cout << "\n"
             << "Processed " << tot << " loci:\n"
             << "  (All loci are always conserved)\n"
             << "  " << ph << " loci had phasing issues (" << pct(ph) << ")\n"
             ;
        if (stats.n_loci_w_pe_reads > 0) {
            assert(!ignore_pe_reads);
            size_t no_pe = stats.n_loci_no_pe_reads() + stats.n_loci_almost_no_pe_reads;
            size_t pe_dag = stats.n_loci_pe_graph_not_dag;
            size_t pe_good = stats.n_loci_usable_pe_reads();
            cout << "  " << no_pe << " loci had no or almost no paired-end reads (" << pct(no_pe) << ")\n"
                 << "  " << pe_dag
                 << " loci had paired-end reads that couldn't be assembled into a contig (" << pct(pe_dag) << ")\n"
                 << "  " << pe_good << " loci had usable paired-end reads (" << pct(pe_good) << ")\n"
                 ;
        }
        cout << "\n";

        logger->l << "BEGIN badly_phased\n"
                  << "n_tot_samples\tn_bad_samples\tn_loci\n";
        for (auto& elem : stats.n_badly_phased_samples)
            logger->l << elem.first.second << '\t' << elem.first.first << '\t' << elem.second << '\n';
        logger->l << "END badly_phased\n\n";
    }

    // Cleanup & return.
    gzclose(o_gzfasta_f);
    model.reset();
    if (dbg_write_hapgraphs)
        o_hapgraphs_f << "}\n";
    cout << prog_name << " is done.\n";
    return 0;

} catch (exception& e) {
    return stacks_handle_exceptions(e);
}
}

SnpAlleleCooccurrenceCounter::SnpAlleleCooccurrenceCounter(size_t n_snps)
    : n_snps_(n_snps),
      cooccurences_(n_snps_*n_snps_) // n*n matrix, athough we only use the i<j half.
    {}

size_t& SnpAlleleCooccurrenceCounter::at(size_t snp_i1, Nt2 snp1_allele, size_t snp_i2, Nt2 snp2_allele) {
    assert(snp_i1 < snp_i2);
    return cooccurences_[snp_i1*n_snps_+snp_i2][size_t(snp1_allele)][size_t(snp2_allele)];
}

void SnpAlleleCooccurrenceCounter::clear() {
    for(size_t i=0; i<n_snps_; ++i)
        for(size_t j=i+1; j<n_snps_; ++j)
            cooccurences_[i*n_snps_+j] = array<array<size_t,4>,4>();
}

size_t ProcessingStats::n_loci_phasing_issues() const {
    size_t n = 0;
    for (auto& elem : n_badly_phased_samples)
        n += elem.second;
    return n;
}

ProcessingStats& ProcessingStats::operator+= (const ProcessingStats& other) {
    n_nonempty_loci += other.n_nonempty_loci;
    n_loci_w_pe_reads += other.n_loci_w_pe_reads;
    n_loci_almost_no_pe_reads += other.n_loci_almost_no_pe_reads;
    n_loci_pe_graph_not_dag += other.n_loci_pe_graph_not_dag;

    for (auto count : other.n_badly_phased_samples)
        n_badly_phased_samples[count.first] += count.second;

    return *this;
}

void LocusProcessor::process(CLocReadSet&& loc) {
    if (loc.reads().empty())
        return;
    ++stats_.n_nonempty_loci;

    //
    // Process the paired-end reads.
    //
    CLocAlnSet pe_aln_loc (loc.mpopi());
    pe_aln_loc.id(loc.id());
    pe_aln_loc.pos(loc.pos());
    if (!ignore_pe_reads) {
        if (!dbg_true_alns) {
            do { // (Avoiding nested ifs.)
                if (loc.pe_reads().empty())
                    break;
                ++stats_.n_loci_w_pe_reads;

                //
                // Assemble the reads.
                //
                string pe_contig;
                vector<const DNASeq4*> seqs_to_assemble;
                for (const Read& r : loc.pe_reads())
                    seqs_to_assemble.push_back(&r.seq);

                Graph graph (km_length);
                graph.rebuild(seqs_to_assemble, min_km_count);
                if (graph.empty()) {
                    ++stats_.n_loci_almost_no_pe_reads;
                    break;
                }

                if (dbg_write_gfa) {
                    graph.dump_gfa(in_dir + to_string(loc.id()) + ".spaths.gfa");
                    graph.dump_gfa(in_dir + to_string(loc.id()) + ".nodes.gfa", true);
                }

                vector<const SPath*> best_path;
                if (!graph.find_best_path(best_path)) {
                    // Not a DAG.
                    ++stats_.n_loci_pe_graph_not_dag;
                    break;
                }

                //
                // Align each read to the contig.
                //
                string ctg = SPath::contig_str(best_path.begin(), best_path.end(), km_length);
                pe_aln_loc.ref(DNASeq4(ctg));

                GappedAln aligner;
                for (SRead& r : loc.pe_reads()) {
                    string seq = r.seq.str();
                    aligner.init(r.seq.length(), ctg.length());
                    aligner.align(seq, ctg);
                    Cigar cigar;
                    parse_cigar(aligner.result().cigar.c_str(), cigar);
                    if (cigar.size() > 10)
                        // Read didn't align, discard it. xxx Refine this.
                        continue;
                    pe_aln_loc.add(SAlnRead(move((Read&)r), move(cigar), r.sample));
                }
            } while (false);

        } else {
            // Use true alignments. Expect "cig2=STRING" in the read IDs.
            if (!loc.pe_reads().empty()) {
                for (SRead& r : loc.pe_reads()) {
                    Cigar cigar = dbg_extract_cigar(r.name);
                    simplify_cigar_to_MDI(cigar);
                    if (pe_aln_loc.ref().empty())
                        pe_aln_loc.ref(DNASeq4(string(cigar_length_ref(cigar), 'N')));
                    if (cigar_length_ref(cigar) != pe_aln_loc.ref().length()) {
                        cerr << "Error: DEBUG: ref-length of cig2 in read '" << r.name
                             << "' was expected to be " << pe_aln_loc.ref().length()
                             << ", not " << cigar_length_ref(cigar) << ".\n";
                        throw exception();
                    }
                    pe_aln_loc.add(SAlnRead(move((Read&)r), move(cigar), r.sample));
                }
            }
        }
    }

    //
    // Build the foward-reads object.
    //
    CLocAlnSet fw_aln_loc (loc.mpopi());
    fw_aln_loc.id(loc.id());
    fw_aln_loc.pos(loc.pos());
    if (!dbg_true_alns) {
        fw_aln_loc.ref(DNASeq4(loc.reads().at(0).seq));
        for (SRead& r : loc.reads()) {
            if (r.seq.length() != fw_aln_loc.ref().length()) {
                cerr << "DEBUG: Error: Can't handle reads of different legnths.\n"; //xxx
                throw exception();
            }
            fw_aln_loc.add(SAlnRead(move((Read&)r), {{'M',r.seq.length()}}, r.sample));
        }
    } else {
        // Use true alignments. Expect "cig1=STRING" in the read IDs.
        for (SRead& r : loc.reads()) {
            Cigar cigar = dbg_extract_cigar(r.name);
            simplify_cigar_to_MDI(cigar);
            if (fw_aln_loc.ref().empty())
                fw_aln_loc.ref(DNASeq4(string(cigar_length_ref(cigar), 'N')));
            if (cigar_length_ref(cigar) != pe_aln_loc.ref().length()) {
                cerr << "Error: DEBUG: ref-length of cig1 in read '" << r.name
                     << "' was expected to be " << pe_aln_loc.ref().length()
                     << ", not " << cigar_length_ref(cigar) << ".\n";
                throw exception();
            }
            fw_aln_loc.add(SAlnRead(move((Read&)r), move(cigar), r.sample));
        }
    }
    loc.clear();

    //
    // Merge the forward & paired-end contigs.
    //
    CLocAlnSet aln_loc (fw_aln_loc.mpopi());
    aln_loc.id(fw_aln_loc.id());
    aln_loc.pos(fw_aln_loc.pos());
    if (!dbg_true_alns) {
        if (pe_aln_loc.reads().empty()) {
            aln_loc = move(fw_aln_loc);
        } else {
            CLocAlnSet dummy (fw_aln_loc.mpopi());
            dummy.id(fw_aln_loc.id());
            dummy.ref(DNASeq4(string(10, 'N')));
            aln_loc = CLocAlnSet::juxtapose(
                    move(fw_aln_loc),
                    CLocAlnSet::juxtapose(move(dummy), move(pe_aln_loc))
                    );
            aln_loc.merge_paired_reads();
        }

    } else {
        // With true alignments, everything should be in the same reference already.
        aln_loc = move(fw_aln_loc);
        if (!pe_aln_loc.reads().empty()) {
            if (pe_aln_loc.ref().length() != aln_loc.ref().length()) {
                cerr << "Error: DEBUG: locus length is " << aln_loc.ref().length()
                     << " for the forward locus half and " << pe_aln_loc.ref().length()
                     << " for the reverse half.\n";
                throw exception();
            }
            for (SAlnRead& r : pe_aln_loc.reads())
                aln_loc.add(move(r));
            aln_loc.merge_paired_reads();
        }
    }
    fw_aln_loc.clear();
    pe_aln_loc.clear();

    if (dbg_write_alns)
        #pragma omp critical
        o_aln_f << "BEGIN " << aln_loc.id() << "\n"
                << aln_loc
                << "\nEND " << aln_loc.id() << "\n";

    //
    // Call SNPs.
    //
    vector<SiteCounts> depths;
    vector<SiteCall> calls;
    calls.reserve(aln_loc.ref().length());
    for(CLocAlnSet::site_iterator site (aln_loc); bool(site); ++site) {
        depths.push_back(site.counts());
        calls.push_back(model->call(depths.back()));
    }

    // Update the consensus sequence.
    DNASeq4 new_ref = aln_loc.ref();
    assert(calls.size() == new_ref.length());
    for (size_t i=0; i<calls.size(); ++i) {
        if (calls[i].alleles().empty())
            new_ref.set(i, Nt4::n);
        else
            new_ref.set(i, calls[i].most_frequent_allele());
    }
    aln_loc.ref(move(new_ref));

    // Call haplotypes.
    vector<map<size_t,PhasedHet>> phase_data;
    if (!dbg_no_haplotypes) {
        set<size_t> inconsistent_samples;
        phase_data = phase_hets(calls, aln_loc, inconsistent_samples);
        if (!inconsistent_samples.empty())
            ++stats_.n_badly_phased_samples[ {inconsistent_samples.size(), aln_loc.n_samples()} ];

        rm_supernumerary_phase_sets(phase_data);
    }

    write_one_locus(aln_loc, depths, calls, phase_data);
}

vector<map<size_t,PhasedHet>> LocusProcessor::phase_hets (
        const vector<SiteCall>& calls,
        const CLocAlnSet& aln_loc,
        set<size_t>& inconsistent_samples
) const {
    vector<map<size_t,PhasedHet>> phased_samples (aln_loc.mpopi().samples().size());

    vector<size_t> snp_cols; // The SNPs of this locus.
    for (size_t i=0; i<aln_loc.ref().length(); ++i)
        if (calls[i].alleles().size() > 1)
            snp_cols.push_back(i);

    if (snp_cols.empty())
        return phased_samples;

    stringstream o_hapgraph_ss;
    if (dbg_write_hapgraphs) {
        o_hapgraph_ss << "subgraph cluster_loc" << aln_loc.id() << " {\n"
                      << "\tlabel=\"locus " << aln_loc.id() << "\";\n"
                      << "\t# snp columns: ";
        join(snp_cols, ',', o_hapgraph_ss);
        o_hapgraph_ss << "\n";
    }

    SnpAlleleCooccurrenceCounter cooccurences (snp_cols.size());
    for (size_t sample=0; sample<aln_loc.mpopi().samples().size(); ++sample) {
        if (aln_loc.sample_reads(sample).empty())
            continue;

        vector<size_t> het_snps; // The heterozygote SNPs of this sample.
        for(size_t snp_i=0; snp_i<snp_cols.size(); ++snp_i)
            if (calls[snp_cols[snp_i]].sample_calls()[sample].call() == snp_type_het)
                het_snps.push_back(snp_i);
        // Check that we have >= 2 hets.
        if (het_snps.size() == 0) {
            continue;
        } else if (het_snps.size() == 1) {
            // Sample has trivial 1nt-long haplotypes.
            size_t col = snp_cols[het_snps[0]];
            const SampleCall& c = calls[col].sample_calls()[sample];
            phased_samples[sample].insert({col, {col, c.nt0(), c.nt1()}});
            continue;
        }

        vector<const SampleCall*> sample_het_calls;
        for (size_t het_i=0; het_i<het_snps.size(); ++het_i)
            sample_het_calls.push_back(&calls[snp_cols[het_snps[het_i]]].sample_calls()[sample]);

        // Iterate over reads, record seen haplotypes (as pairwise cooccurrences).
        cooccurences.clear();
        vector<Nt4> read_hap (het_snps.size());
        for (size_t read_i : aln_loc.sample_reads(sample)) {
            // Build the haplotype.
            size_t curr_col = 0;
            Alignment::iterator read_itr = aln_loc.reads()[read_i].aln();
            for (size_t het_i=0; het_i<het_snps.size(); ++het_i) {
                size_t col = snp_cols[het_snps[het_i]];
                read_itr += col - curr_col;
                curr_col = col;
                Nt4 nt = *read_itr;
                if (nt == Nt4::n) {
                    read_hap[het_i] = Nt4::n;
                } else {
                    const SampleCall& c = *sample_het_calls[het_i];
                    Nt2 nt2 = Nt2(nt);
                    if (nt2 == c.nt0() || nt2 == c.nt1())
                        read_hap[het_i] = nt;
                    else
                        read_hap[het_i] = Nt4::n;
                }
            }

            // Record the pairwise cooccurrences.
            for (size_t i=0; i<het_snps.size(); ++i) {
                Nt4 nti = read_hap[i];
                if (nti == Nt4::n)
                    continue;
                for (size_t j=i+1; j<het_snps.size(); ++j) {
                    Nt4 ntj = read_hap[j];
                    if (ntj == Nt4::n)
                        continue;
                    ++cooccurences.at(het_snps[i], Nt2(nti), het_snps[j], Nt2(ntj));
                }
            }
        }

        // Write the dot graph, if required.
        if (dbg_write_hapgraphs) {
            auto nodeid = [&aln_loc,&sample](size_t col, Nt2 allele)
                    {return string("l")+to_string(aln_loc.id())+"s"+to_string(sample)+"c"+to_string(col)+char(allele);};

            // Initialize the subgraph.
            size_t n_reads = aln_loc.sample_reads(sample).size();
            size_t n_merged_reads = 0;
            for (size_t read_i : aln_loc.sample_reads(sample))
                if (aln_loc.reads()[read_i].name.back() == 'm')
                    ++n_merged_reads;
            const string& sample_name = aln_loc.mpopi().samples()[sample].name;
            o_hapgraph_ss << "\tsubgraph cluster_sample" << sample << " {\n"
                          << "\t\tlabel=\""
                          << "i" << sample << " '" << sample_name << "'\\n"
                          << "nreads=" << n_reads << ",merged=" << n_merged_reads
                          << "\";\n"
                          << "\t\tstyle=dashed;\n"
                          << "\t\t# heterozygous columns: ";
            vector<size_t> het_cols;
            for (size_t snp_i : het_snps)
                het_cols.push_back(snp_cols[snp_i]);
            join(het_cols, ',', o_hapgraph_ss);
            o_hapgraph_ss << "\n";

            // Write the node labels.
            for (size_t i=0; i<het_snps.size(); ++i) {
                array<Nt2,2> alleles = sample_het_calls[i]->nts();
                size_t col = snp_cols[het_snps[i]];
                for (Nt2 allele : alleles)
                    o_hapgraph_ss << "\t\t" << nodeid(col, allele)
                                  << " [label=<"
                                  << "<sup><font point-size=\"10\">" << col+1 << "</font></sup>" << allele
                                  << ">];\n";
            }

            // Write the edges.
            for (size_t het_i=0; het_i<het_snps.size(); ++het_i) {
                array<Nt2,2> alleles_i = sample_het_calls[het_i]->nts();
                size_t snp_i = het_snps[het_i];
                for (size_t het_j=het_i+1; het_j<het_snps.size(); ++het_j) {
                    array<Nt2,2> alleles_j = sample_het_calls[het_j]->nts();
                    size_t snp_j = het_snps[het_j];
                    for (Nt2 nti : alleles_i) {
                        for (Nt2 ntj : alleles_j) {
                            size_t n = cooccurences.at(snp_i, nti, snp_j, ntj);
                            if (n == 0)
                                continue;
                            o_hapgraph_ss << "\t\t" << nodeid(snp_cols[snp_i],nti)
                                          << " -- " << nodeid(snp_cols[snp_j],ntj) << " [";
                            if (n==1)
                                o_hapgraph_ss << "style=dotted";
                            else
                                o_hapgraph_ss << "label=\"" << n << "\",penwidth=" << n;
                            o_hapgraph_ss << "];\n";
                        }
                    }
                }
            }
            o_hapgraph_ss << "\t}\n";
        }

        // Call haplotypes.
        // This is based on a graph of cooccurrences in which nodes are the het
        // alleles. Subgraphs represent haplotypes, and may not include more than
        // one node for each SNP.
        vector<vector<Nt4>> haps;
        bool inconsistent = false;
        // We keep track of which haplotype each allele is currently part of.
        vector<array<size_t,4>> allele_to_hap (het_snps.size(), {SIZE_MAX,SIZE_MAX,SIZE_MAX,SIZE_MAX});
        for (size_t het_i=0; het_i<het_snps.size(); ++het_i) {
            array<Nt2,2> allelesi = sample_het_calls[het_i]->nts();
            size_t snp_i = het_snps[het_i];
            for (size_t het_j=het_i+1; het_j<het_snps.size(); ++het_j) {
                array<Nt2,2> allelesj = sample_het_calls[het_j]->nts();
                size_t snp_j = het_snps[het_j];
                for (Nt2 nti : allelesi) {
                    for (Nt2 ntj : allelesj) {
                        size_t n = cooccurences.at(snp_i, nti, snp_j, ntj);
                        if (n == 0)
                            continue;

                        // Discard low-coverage edges.
                        const size_t min_n_cooccurrences = 2;
                        if (n < min_n_cooccurrences)
                            continue;

                        size_t& hap_i = allele_to_hap[het_i][size_t(nti)];
                        size_t& hap_j = allele_to_hap[het_j][size_t(ntj)];
                        if (hap_i == SIZE_MAX && hap_j == SIZE_MAX) {
                            // Both nodes are singletons. Start a new haplotype.
                            hap_i = haps.size();
                            hap_j = haps.size();
                            haps.push_back(vector<Nt4>(het_snps.size(), Nt4::n));
                            haps.back()[het_i] = Nt4(nti);
                            haps.back()[het_j] = Nt4(ntj);
                        } else if (hap_i == hap_j) {
                            // Nodes are already in the same haplotype/subgraph.
                            // Nothing to do.
                        } else if (hap_j == SIZE_MAX) {
                            // Add ntj to nti's haplotype (i.e. subgraph).
                            vector<Nt4>& hap = haps[hap_i];
                            assert(hap[het_i] == Nt4(nti));
                            if (hap[het_j] == Nt4::n) {
                                hap[het_j] = Nt4(ntj);
                                hap_j = hap_i;
                            } else {
                                assert(hap[het_j] != Nt4(ntj)); // As `cooccurrences` is a half-matrix.
                                // allele nti at position het_i (i.e. column `snp_cols[het_snps[het_i]]`)
                                // is already phased with another allele at position het_j.
                                inconsistent = true;
                                break;
                            }
                        } else if (hap_i == SIZE_MAX) {
                            // Same as immediately above, reversed (add nti to ntj's haplotype).
                            vector<Nt4>& hap = haps[hap_j];
                            assert(hap[het_j] == Nt4(ntj));
                            if (hap[het_i] == Nt4::n) {
                                hap[het_i] = Nt4(nti);
                                hap_i = hap_j;
                            } else {
                                assert(hap[het_i] != Nt4(nti));
                                inconsistent = true;
                                break;
                            }
                        } else {
                            // Both nodes are already in a haplotype. The haplotype graph
                            // is consistent only if the positions spanned by the two
                            // haplotypes are disjoint.
                            assert(haps[hap_i][het_i] == Nt4(nti));
                            assert(haps[hap_j][het_j] == Nt4(ntj));
                            for (size_t k=0; k<het_snps.size(); ++k) {
                                if (haps[hap_i][k] != Nt4::n && haps[hap_j][k] != Nt4::n) {
                                    assert(haps[hap_i][k] != haps[hap_j][k]);
                                    inconsistent = true;
                                    break;
                                }
                            }
                            if (inconsistent) {
                                break;
                            } else {
                                // Merge the two haplotypes/subgraphs.
                                vector<Nt4>& hap = haps[hap_i];
                                vector<Nt4>& rm_hap = haps[hap_j];
                                assert(haps[hap_i][het_j] == Nt4::n);
                                assert(haps[hap_j][het_i] == Nt4::n);
                                for (size_t k=0; k<het_snps.size(); ++k) {
                                    if (rm_hap[k] != Nt4::n) {
                                        // Transfer the allele to the kept haplotype.
                                        hap[k] = rm_hap[k];
                                        allele_to_hap[k][size_t(Nt2(hap[k]))] = hap_i;
                                    }
                                }
                                rm_hap = vector<Nt4>();
                                #ifdef DEBUG
                                // Check that the discarded haplotype has become inaccessible.
                                size_t rm_hap_i = &rm_hap - haps.data();
                                for (auto& het : allele_to_hap)
                                    for (auto& nt : het)
                                        assert(nt != rm_hap_i);
                                #endif
                            }
                        }
                    }
                    if (inconsistent)
                        break;
                }
                if (inconsistent)
                    break;
            }
            if (inconsistent)
                break;
        }

        if (inconsistent) {
            inconsistent_samples.insert(sample);
        } else {
            // Record haplotypes.
            // Iterate over pairs of haplotypes; overlaps become 'phase sets'.
            // There may be three or more haplotypes but in this case, because
            // these haplotypes are compatible by construction, phased/overlapping
            // columns are disjoint across all pairs.
            for (size_t i=0; i<haps.size(); ++i) {
                if (haps[i].empty())
                    // Deleted remnant of a merger.
                    continue;

                for (size_t j=i+1; j<haps.size(); ++j) {
                    if (haps[j].empty())
                        continue;

                    size_t phase_set = SIZE_MAX;
                    for (size_t het_i=0; het_i<het_snps.size(); ++het_i) {
                        Nt4 hapi_nt = haps[i][het_i];
                        Nt4 hapj_nt = haps[j][het_i];
                        if (hapi_nt == Nt4::n || hapj_nt == Nt4::n)
                            continue;

                        size_t col = snp_cols[het_snps[het_i]];
                        if (phase_set == SIZE_MAX)
                            // This is the first observation for this pair of haplotypes.
                            phase_set = col;

                        assert(!phased_samples[sample].count(col)); // Disjoint resolved columns across pairs.
                        phased_samples[sample][col] = PhasedHet({phase_set, hapi_nt, hapj_nt});
                    }
                }
            }
            // Record singleton nodes (that are implicit in our representation).
            for (size_t het_i=0; het_i<het_snps.size(); ++het_i) {
                size_t col = snp_cols[het_snps[het_i]];
                if (!phased_samples[sample].count(col)) {
                    array<Nt2,2> alleles = sample_het_calls[het_i]->nts();
                    phased_samples[sample][col] = PhasedHet({col, alleles[0], alleles[1]});
                }
            }
            assert(!phased_samples[sample].empty());
        }
    }
    if (dbg_write_hapgraphs) {
        o_hapgraph_ss << "}\n";
        #pragma omp critical
        o_hapgraphs_f << o_hapgraph_ss.str();
    }

    return phased_samples;
}

void LocusProcessor::rm_supernumerary_phase_sets (
        vector<map<size_t,PhasedHet>>& phase_data
) const {
    // We only keep one phase set per sample.
    for (auto& sample : phase_data) {
        if (sample.empty())
            continue;
        // Tally phase sets & find the one with the most SNPs.
        map<size_t,size_t> phase_set_sizes;
        for (auto& phasedhet : sample)
            ++phase_set_sizes[phasedhet.second.phase_set];
        auto best_ps = phase_set_sizes.begin();
        for (auto ps=++phase_set_sizes.begin(); ps!=phase_set_sizes.end(); ++ps)
            if (ps->second > best_ps->second)
                best_ps = ps;
        // Remove all phase sets but one.
        for (auto phasedhet=sample.begin(); phasedhet!=sample.end();) {
            if (phasedhet->second.phase_set == best_ps->first)
                ++phasedhet;
            else
                sample.erase(phasedhet++);
        }
    }
}

void LocusProcessor::write_one_locus (
        const CLocAlnSet& aln_loc,
        const vector<SiteCounts>& depths,
        const vector<SiteCall>& calls,
        const vector<map<size_t,PhasedHet>>& phase_data
) {
    char loc_id[16];
    sprintf(loc_id, "%d", aln_loc.id());

    const DNASeq4& ref = aln_loc.ref();
    const MetaPopInfo& mpopi = aln_loc.mpopi();

    //
    // Vcf output.
    //
    assert(depths.size() == ref.length());
    assert(calls.size() == ref.length());
    vector<size_t> sample_sites_w_data (mpopi.samples().size(), 0);
    stringstream vcf_records;
    VcfRecord rec;
    for (size_t i=0; i<ref.length(); ++i) {
        const SiteCounts& sitedepths = depths[i];
        const SiteCall& sitecall = calls[i];
        if (sitecall.alleles().empty())
            // No useful data at this site.
            continue;

        if (!ref[i].is_acgt())
            continue;

        // Determine which alleles exist, and their order.
        // (n.b. As of Apr4,2017 the ref allele might not be the most frequent one.)
        vector<Nt2> vcf_alleles;
        map<Nt2, size_t> vcf_allele_indexes;
        {
            vcf_alleles.push_back(ref[i]);
            vcf_allele_indexes.insert({Nt2(ref[i]), 0});

            // Sort the alleles by frequency.
            vector<pair<double, Nt2>> sorted_alleles;
            for (auto& a : sitecall.alleles())
                sorted_alleles.push_back({a.second, a.first});
            sort(sorted_alleles.rbegin(), sorted_alleles.rend()); // (decreasing)

            // The reference allele has already been added to vcf_alleles; exclude it.
            for (auto iter=sorted_alleles.begin(); iter!=sorted_alleles.end(); ++iter) {
                if (iter->second == vcf_alleles[0]) {
                    sorted_alleles.erase(iter);
                    break;
                }
            }

            // Record the alternative alleles.
            for (auto& alt_allele : sorted_alleles) {
                vcf_allele_indexes.insert({alt_allele.second, vcf_alleles.size()});
                vcf_alleles.push_back(alt_allele.second);
            }
        }

        //
        // Create the VCF record.
        //

        // Chrom & pos.
        rec.clear();
        rec.append_chrom(loc_id);
        rec.append_pos(i+1);
        rec.append_id();

        // Alleles.
        for (Nt2 nt : vcf_alleles)
            rec.append_allele(nt);

        rec.append_qual();
        rec.append_filters();

        if(vcf_alleles.size() == 1) {
            // Fixed site.

            // Info/DP.
            rec.append_info(string("DP=") + to_string(sitedepths.tot.sum()));
            // Info/AD.
            Nt2 ref_nt = sitecall.alleles().begin()->first;
            rec.append_info(string("AD=") + to_string(sitedepths.tot[ref_nt]));
            if (dbg_write_nt_depths) {
                // Info/cnts.
                stringstream cnts;
                join(sitedepths.tot.arr(), ',', cnts);
                rec.append_info(string("cnts=") + cnts.str());
            }
            // Format.
            rec.append_format("DP");
            // Genotypes.
            for (size_t sample=0; sample<mpopi.samples().size(); ++sample) {
                size_t dp = sitedepths.samples[sample].sum();
                if (dp == 0) {
                    rec.append_sample(".");
                    continue;
                }
                ++sample_sites_w_data[sample];
                rec.append_sample(to_string(dp));
            }

        } else {
            // Polymorphic site.

            // Info/DP.
            rec.append_info(string("DP=") + to_string(sitedepths.tot.sum()));
            // Info/AD.
            vector<size_t> ad;
            for (Nt2 nt : vcf_alleles)
                ad.push_back(sitedepths.tot[nt]);
            stringstream ss;
            join(ad, ',', ss);
            rec.append_info(string("AD=") + ss.str());
            // Info/AF.
            vector<double> alt_freqs;
            for (auto nt=++vcf_alleles.begin(); nt!=vcf_alleles.end(); ++nt) // rem. always >1 alleles.
                alt_freqs.push_back(sitecall.alleles().at(*nt));
            rec.append_info(VcfRecord::util::fmt_info_af(alt_freqs));
            if (dbg_write_nt_depths) {
                // Info/cnts.
                stringstream cnts;
                join(sitedepths.tot.arr(), ',', cnts);
                rec.append_info(string("cnts=") + cnts.str());
            }

            // Format.
            rec.append_format("GT");
            if (!dbg_no_haplotypes)
                rec.append_format("PS"); // Phase set.
            rec.append_format("DP");
            rec.append_format("AD");
            rec.append_format("GL");

            // Genotypes.
            for (size_t sample=0; sample<mpopi.samples().size(); ++sample) {
                const Counts<Nt2>& sdepths = sitedepths.samples[sample];
                const SampleCall& scall = sitecall.sample_calls()[sample];
                if (sdepths.sum() == 0) {
                    // No data for this sample.
                    rec.append_sample(".");
                    continue;
                }
                ++sample_sites_w_data[sample];

                stringstream genotype;
                // GT field.
                array<size_t,2> gt;
                switch (scall.call()) {
                case snp_type_hom:
                    gt[0] = vcf_allele_indexes.at(scall.nt0());
                    genotype << gt[0] << '/' << gt[0];
                    if (!dbg_no_haplotypes)
                        genotype << ":.";
                    break;
                case snp_type_het:
                    if (!dbg_no_haplotypes) {
                        assert(!phase_data.empty());
                        auto itr = phase_data[sample].find(i);
                        if (itr != phase_data[sample].end()) {
                            // There is phase data for this heterozygous site.
                            const PhasedHet& p = itr->second;
                            genotype << vcf_allele_indexes.at(p.left_allele)
                                     << '|'
                                     << vcf_allele_indexes.at(p.right_allele)
                                     << ':' << (p.phase_set + 1);
                        } else {
                            // No phase data.
                            gt[0] = vcf_allele_indexes.at(scall.nt0());
                            gt[1] = vcf_allele_indexes.at(scall.nt1());
                            sort(gt.begin(), gt.end()); // (Prevents '1/0'.)
                            genotype << gt[0] << '/' << gt[1];
                            genotype << ":.";
                        }
                    } else {
                        gt[0] = vcf_allele_indexes.at(scall.nt0());
                        gt[1] = vcf_allele_indexes.at(scall.nt1());
                        sort(gt.begin(), gt.end()); // (Prevents '1/0'.)
                        genotype << gt[0] << '/' << gt[1];
                    }
                    break;
                default:
                    genotype << '.';
                    if (!dbg_no_haplotypes)
                        genotype << ":.";
                    break;
                }
                // DP field.
                genotype << ':' << sdepths.sum();
                // AD field.
                vector<size_t> ad;
                ad.reserve(vcf_alleles.size());
                for (Nt2 nt : vcf_alleles)
                    ad.push_back(sdepths[nt]);
                genotype << ':';
                join(ad, ',', genotype);
                // GL field.
                genotype << ':' << VcfRecord::util::fmt_gt_gl(vcf_alleles, scall.lnls());
                // cnts field.
                if (dbg_write_nt_depths) {
                    genotype << ":";
                    join(sdepths.arr(), ',', genotype);
                }
                // Push it.
                rec.append_sample(genotype.str());
            }
        }
        assert(rec.count_samples() == o_vcf_f->header().samples().size());

        vcf_records << rec;
    }
    out_.vcf = vcf_records.str();

    //
    // Fasta output.
    //

    // Determine the number of samples for this locus. Some samples may have
    // been discarded (as of May 26, 2017, this would be because their haplotypes
    // were inconsistent).
    set<size_t> samples_w_reads;
    for (const SAlnRead& r : aln_loc.reads())
        samples_w_reads.insert(r.sample);
    size_t n_remaining_samples = 0;
    for (size_t sample_n_sites : sample_sites_w_data)
        if (sample_n_sites > 0)
            ++n_remaining_samples;

    // Write the record.
    out_.fa.clear();
    out_.fa += '>';
    out_.fa += loc_id;
    if (!aln_loc.pos().empty()) {
        const PhyLoc& p = aln_loc.pos();
        char pos[16];
        sprintf(pos, "%u", p.bp+1);
        out_.fa += " pos=";
        out_.fa += p.chr();
        out_.fa += ':';
        out_.fa += pos;
        out_.fa += ':';
        out_.fa += (p.strand == strand_plus ? '+' : '-');
    }
    char n_spls[32];
    sprintf(n_spls, "%zu", n_remaining_samples);
    out_.fa += " NS=";
    out_.fa += n_spls;
    if (n_remaining_samples != samples_w_reads.size()) {
        assert(n_remaining_samples < samples_w_reads.size());
        sprintf(n_spls, "%zu", samples_w_reads.size() - n_remaining_samples);
        out_.fa += " n_discarded_samples=";
        out_.fa += n_spls;
    }
    out_.fa += '\n';
    out_.fa += ref.str();
    out_.fa += '\n';
}

Cigar dbg_extract_cigar(const string& read_id) {
    static const char keyword1[] = "cig1=";
    static const char keyword2[] = "cig2=";
    static const size_t kw_len = sizeof(keyword1) - 1;
    static_assert(sizeof(keyword1) == sizeof(keyword2), "");

    Cigar cigar;

    // Find the start.
    const char* cig_start = read_id.c_str();
    const char* kw = keyword1;
    if (read_id.back() == '2' && read_id.at(read_id.length()-2) == '/')
        kw = keyword2;
    while (*cig_start != '\0'
            && ! (*cig_start == kw[0] && strncmp(cig_start, kw, kw_len) == 0))
        ++cig_start;
    if (*cig_start == '\0') {
        cerr << "Error: DEBUG: Coulnd't find cigar in read ID '" << read_id
             << "'; expected 'cigar=.+[^A-Z0-9]'.\n";
        throw exception();
    }
    cig_start += kw_len;

    // Find the end.
    const char* cig_past = cig_start;
    while(*cig_past != '\0' && is_cigar_char[uint(*cig_past)])
        ++cig_past;

    // Extract the cigar.
    parse_cigar(string(cig_start, cig_past).c_str(), cigar, true);

    return cigar;
}

const string help_string = string() +
        prog_name + " " + VERSION  + "\n" +
        prog_name + " -P in_dir\n"
        "\n"
        "  -P: input directory\n"
        "  -b: batch ID (default: guess)\n"
        "  The input directory must contain a 'batch_X.catalog.bam' file, that\n"
        "  the user should generate after running tsv2bam with e.g.:\n"
        "  \"samtools merge ./batch_1.catalog.bam ./*.matches.bam\"\n"
        "\n"
        "  --ignore-pe-reads: ignore paired-end reads even if present in the input\n"
        "  -W,--whitelist: path to a whitelist of locus IDs\n"
        "  -t,--threads: number of threads to use (default: 1)\n"
        "\n"
        "Model options:\n"
        "  --model: model to use to call variants and genotypes; one of\n"
        "           marukilow (default), marukihigh, or snp\n"
        "  --var-alpha: alpha threshold for discovering SNPs (default: 0.05)\n"
        "  --gt-alpha: alpha threshold for calling genotypes (default: 0.05)\n"
        "\n"
        "Expert options:\n"
        "  --kmer-length: kmer length for the de Bruijn graph (default: 31)\n"
        "  --min-kmer-cov: minimum coverage to consider a kmer (default: 2)\n"
        "\n"
#ifdef DEBUG
        "Debug options:\n"
        "  --no-haps: disable phasing\n"
        "  --gfa: output a GFA file for each locus\n"
        "  --alns: output a file showing the contigs & alignments\n"
        "  --hap-graphs: output a dot graph file showing phasing information\n"
        "  --depths: write detailed depth data in the output VCF\n"
        "  --true-alns: use true alignments (for simulated data; read IDs must\n"
        "               include 'cig1=...' and 'cig2=...' fields.\n"
        "\n"
#endif
        ;

void parse_command_line(int argc, char* argv[]) {

    auto bad_args = [](){
        cerr << help_string;
        exit(13);
    };

try {

    static const option long_options[] = {
        {"version",      no_argument,       NULL,  1000},
        {"help",         no_argument,       NULL,  'h'},
        {"quiet",        no_argument,       NULL,  'q'},
        {"in-dir",       required_argument, NULL,  'P'},
        {"batch-id",     required_argument, NULL,  'b'},
        {"whitelist",    required_argument, NULL,  'W'},
        {"threads",      required_argument, NULL,  't'},
        {"ignore-pe-reads", no_argument,    NULL,  1012},
        {"model",        required_argument, NULL,  1006},
        {"gt-alpha",     required_argument, NULL,  1005},
        {"var-alpha",    required_argument, NULL,  1008},
        {"kmer-length",  required_argument, NULL,  1001},
        {"min-kmer-cov", required_argument, NULL,  1002},
        {"true-alns",    no_argument,       NULL,  1011},
        {"no-haps",      no_argument,       NULL,  1009},
        {"gfa",          no_argument,       NULL,  1003},
        {"alns",         no_argument,       NULL,  1004},
        {"hap-graphs",   no_argument,       NULL,  1010},
        {"depths",       no_argument,       NULL,  1007},
        {0, 0, 0, 0}
    };

    double gt_alpha = 0.05;
    double var_alpha = 0.05;
    string wl_path;

    int c;
    int long_options_i;
    while (true) {

        c = getopt_long(argc, argv, "hqP:b:W:t:", long_options, &long_options_i);

        if (c == -1)
            break;

        switch (c) {
        case 1000: //version
            cout << prog_name << " " << VERSION << "\n";
            exit(0);
            break;
        case 'h':
            cout << help_string;
            exit(0);
            break;
        case 'q':
            quiet = true;
            break;
        case 'P':
            in_dir = optarg;
            break;
        case 'b':
            batch_id = atoi(optarg);
            break;
        case 1006: //model
            model_type = parse_model_type(optarg);
            break;
        case 1005: //gt-alpha
            gt_alpha = atof(optarg);
            break;
        case 1008: //var-alpha
            var_alpha = atof(optarg);
            break;
        case 'W':
            wl_path = optarg;
            break;
        case 't':
            num_threads = is_integer(optarg);
            if (num_threads < 0) {
                cerr << "Error: Illegal -t option value '" << optarg << "'.\n";
                bad_args();
            }
            break;
        case 1012: //ignore-pe-reads
            ignore_pe_reads = true;
            break;
        case 1001://kmer-length
            km_length = atoi(optarg);
            break;
        case 1002://min-cov
            min_km_count = atoi(optarg);
            break;
        case 1011://true-alns
            dbg_true_alns = true;
            break;
        case 1009://no-haps
            dbg_no_haplotypes = true;
            break;
        case 1003://gfa
            dbg_write_gfa = true;
            break;
        case 1004://aln
            dbg_write_alns = true;
            break;
        case 1010://hap-graphs
            dbg_write_hapgraphs = true;
            break;
        case 1007://depths
            dbg_write_nt_depths = true;
            break;
        case '?':
            bad_args();
            break;
        }
    }

    // Check command consistency.
    if (optind < argc) {
        cerr << "Error: Failed to parse command line: '" << argv[optind]
             << "' is seen as a positional argument. Expected no positional arguments.\n";
        bad_args();
    }
    if (in_dir.empty()) {
        cerr << "Error: An input directory must be provided (-P).\n";
        bad_args();
    }

    switch (model_type) {
    case snp:        model.reset(new MultinomialModel(gt_alpha)); break;
    case marukihigh: model.reset(new MarukiHighModel(gt_alpha, var_alpha));  break;
    case marukilow:  model.reset(new MarukiLowModel(gt_alpha, var_alpha));   break;
    default:
        cerr << "Error: Model choice '" << to_string(model_type) << "' is not supported.\n";
        bad_args();
        break;
    }

    // Process arguments.
    if (in_dir.back() != '/')
        in_dir += '/';

    if (batch_id < 0) {
        vector<int> cat_ids = find_catalogs(in_dir);
        if (cat_ids.size() == 1) {
            batch_id = cat_ids[0];
        } else if (cat_ids.empty()) {
            cerr << "Error: Unable to find a catalog in '" << in_dir << "'.\n";
            bad_args();
        } else {
            cerr << "Error: Input directory contains several catalogs, please specify -b.\n";
            bad_args();
        }
    }

    if (!wl_path.empty()) {
        ifstream wl_fh (wl_path);
        check_open(wl_fh, wl_path);
        int id;
        while (wl_fh >> id)
            locus_wl.insert(id);
        if (locus_wl.empty()) {
            cerr << "Error: Whitelist '" << wl_path << "' appears empty.\n";
            throw exception();
        }
    }

} catch (std::invalid_argument&) {
    bad_args();
}
}

void report_options(ostream& os) {
    os << "Configuration for this run:\n"
       << "  Input directory: '" << in_dir << "'\n"
       << "  Batch ID: " << batch_id << "\n"
       << "  Model: " << *model << "\n";

    if (!locus_wl.empty())
        os << "  Whitelist of " << locus_wl.size() << " loci.\n";
    if (ignore_pe_reads)
        os << "  Ignoring paired-end reads.\n";
    if (km_length != 31)
        os << "  Kmer length: " << km_length << "\n";
    if (min_km_count != 2)
        os << "  Min coverage: " << min_km_count << "\n";
}