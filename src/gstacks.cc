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
#include <getopt.h>
#include <zlib.h>

#include "gstacks.h"

#include "constants.h"
#include "utils.h"
#include "log_utils.h"
#include "catalog_utils.h"
#include "locus.h"
#include "locus_readers.h"
#include "debruijn.h"
#include "aln_utils.h"
#include "Alignment.h"
#include "models.h"

//
// Argument globals.
//
GStacksInputT  input_type = GStacksInputT::unknown;
string         popmap_path;  // Set if --popmap is given. For logging purposes.
vector<string> sample_names; // Set if --popmap is given.
vector<string> in_bams;      // Set if -B is given, or if -P is given without -M.

string o_prefix;

BamCLocBuilder::Config refbased_cfg {false, 1000, 10, 0.20, 1, false};

int    num_threads       = 1;
bool   quiet             = false;
bool   ignore_pe_reads   = false;
double min_aln_cov       = 0.75;
int    min_se_pe_overlap = 5;
bool   detailed_output   = false;

modelt model_type = marukilow;
unique_ptr<const Model> model;

size_t km_length         = 31;
size_t min_km_count      = 2;
size_t max_fragment_alns = 2;

bool   dbg_no_overlaps     = false;
bool   dbg_no_haplotypes   = false;
bool   dbg_write_gfa       = false;
bool   dbg_write_alns      = false;
bool   dbg_write_hapgraphs = false;
bool   dbg_write_nt_depths = false;
bool   dbg_true_reference  = false;
bool   dbg_true_alns       = false;

//
// Additional globals.
//
const string prog_name = "gstacks";
unique_ptr<LogAlterator> logger;
gzFile o_gzfasta_f = NULL;
unique_ptr<VcfWriter> o_vcf_f;
unique_ptr<VersatileWriter> o_details_f;

ofstream o_aln_f;
ofstream o_hapgraphs_f;
const char o_aln_header[] =
    "# This prints observed read haplotypes:\n"
    "# show_loc() { loc=$1; cat ${RY_DIR:-.}/gstacks.alns | sed -n \"/^END $loc\\b/ q; /^BEGIN $loc\\b/,$ p\" | tail -n+2; }\n"
    "# snp_cols() { loc=$1; zcat ${RY_DIR:-.}/gstacks.vcf.gz | awk \"\\$1==$loc; \\$1>$loc {exit}\" | awk '$5!=\".\"' | cut -f2 | paste -sd ','; }\n"
    "# show_haps() { loc=$1; cols=$2; spl=$3; show_loc $loc | grep \"\\b$spl\\b\" | cut -f3 | cut -c \"$cols\" | sort; }\n"
    "# true_loci() { loc=$1; spl=$2; show_loc $loc | grep \"\\b$spl\\b\" | grep -v ref | cut -d: -f1 | sort -u; }\n"
    ;
const char o_hapgraphs_header[] =
    "# dot -Tpdf -O gstacks.hapgraphs.dot\n"
    "# loc=371\n"
    "# { g=gstacks.hapgraphs.dot; sed -n '0,/^subgraph/p' $g | head -n-1; sed -n \"/^subgraph cluster_loc$loc\\b/,/^}/p\" $g; echo \\}; } | dot -Tpdf -o haps.$loc.pdf\n"
    "graph {\n"
    "edge[color=\"grey60\",fontsize=12,labeljust=\"l\"];\n";
    ;

//
// main
// ==========
//

int
main(int argc, char** argv)
{
try {

    //
    // Parse arguments.
    //
    parse_command_line(argc, argv);

    //
    // Open the input BAM file(s).
    //
    vector<Bam*> bam_f_ptrs;
    for (const string& in_bam : in_bams)
        bam_f_ptrs.push_back(new Bam(in_bam.c_str()));

    unique_ptr<BamCLocReader> bam_cloc_reader;
    unique_ptr<BamCLocBuilder> bam_cloc_builder;
    const MetaPopInfo* bam_mpopi;
    if (input_type == GStacksInputT::denovo_popmap || input_type == GStacksInputT::denovo_merger) {
        bam_cloc_reader.reset(new BamCLocReader(move(bam_f_ptrs)));
        // (Note: For de novo (Stacks) input we use read groups, even when a
        // popmap was specified, since the input files ought to have them and
        // this allows to detect corruption.)
        bam_mpopi = &bam_cloc_reader->mpopi();
    } else if (input_type == GStacksInputT::refbased_popmap || input_type == GStacksInputT::refbased_list) {
        bam_cloc_builder.reset(new BamCLocBuilder(move(bam_f_ptrs), refbased_cfg, sample_names));
        bam_mpopi = &bam_cloc_builder->mpopi();
    }

    //
    // Open the log.
    //
    logger.reset(new LogAlterator(o_prefix + ".log", quiet, argc, argv));
    report_options(cout);
    cout << "\n" << flush;

    //
    // Initialize OPENMP.
    //
    #ifdef _OPENMP
    omp_set_num_threads(num_threads);
    #endif

    //
    // Open the output files.
    //
    string o_gzfasta_path = o_prefix + ".fa.gz";
    o_gzfasta_f = gzopen(o_gzfasta_path.c_str(), "wb");
    check_open(o_gzfasta_f, o_gzfasta_path);

    VcfHeader vcf_header;
    vcf_header.add_std_meta();
    for(auto& s : bam_mpopi->samples())
        vcf_header.add_sample(s.name);
    o_vcf_f.reset(new VcfWriter(o_prefix + ".vcf.gz", move(vcf_header)));

    if (detailed_output) {
        o_details_f.reset(new VersatileWriter(o_prefix + ".details.gz"));
        *o_details_f << "# show_loc() { loc=$1; zcat ./gstacks.details.gz | sed -rn \"/^BEGIN locus $loc\\b/,\\$ p; /^END locus $loc\\b/ q;\"; }\n";
    }

    if (dbg_write_alns) {
        string o_aln_path = o_prefix + ".alns";
        o_aln_f.open(o_aln_path);
        check_open(o_aln_f, o_aln_path);
        o_aln_f << o_aln_header;
    }

    if (dbg_write_hapgraphs) {
        string o_hapgraphs_path = o_prefix + ".hapgraphs.dot";
        o_hapgraphs_f.open(o_hapgraphs_path);
        check_open(o_hapgraphs_f, o_hapgraphs_path);
        o_hapgraphs_f << o_hapgraphs_header;
    }

    //
    // Process every locus
    //
    GenotypeStats gt_stats {};

    // For clocking.
    Timers t_threads_totals;
    Timer t_parallel;
    Timer t_writing_vcf;
    size_t n_writes = 0;
    size_t max_size_before_write = 0;

    // For parallelization.
    int omp_return = 0;
    std::deque<pair<bool,string>> fa_outputs;
    std::deque<pair<bool,string>> vcf_outputs;
    std::deque<pair<bool,string>> det_outputs; // (details)
    size_t next_fa_to_write = 0; // locus index, in the input BAM file.
    size_t next_vcf_to_write = 0;
    size_t next_det_to_write = 0;

    if (input_type == GStacksInputT::denovo_popmap || input_type == GStacksInputT::denovo_merger) {

        const size_t n_loci = bam_cloc_reader->n_loci();
        ContigStats ctg_stats {};

        cout << "Processing all loci...\n" << flush;
        ProgressMeter progress (cout, true, n_loci);

        t_parallel.restart();
        #pragma omp parallel
        {
            LocusProcessor loc_proc;
            Timers& t = loc_proc.timers();
            CLocReadSet loc (*bam_mpopi);

            #pragma omp for schedule(dynamic)
            for (size_t i=0; i<n_loci; ++i) {
                if (omp_return != 0)
                    continue;
                try {
                    // Read a locus from the BAM file.
                    t.reading.restart();
                    #pragma omp critical(read)
                    try {
                        if (omp_return == 0)
                            bam_cloc_reader->read_one_locus(loc);
                    } catch (exception& e) {
                        omp_return = stacks_handle_exceptions(e);
                    }
                if (omp_return == 0) {
                    size_t loc_i = loc.bam_i();
                    t.reading.stop();

                    // Process it.
                    t.processing.restart();
                    loc_proc.process(loc);
                    t.processing.stop();

                    // Write the FASTA output.
                    t.writing_fa.restart();
                    #pragma omp critical(write_fa)
                    {
                        for (size_t i=next_fa_to_write+fa_outputs.size(); i<=loc_i; ++i)
                            fa_outputs.push_back( {false, string()} );
                        fa_outputs[loc_i - next_fa_to_write] = {true, move(loc_proc.fasta_out())};

                        while (!fa_outputs.empty() && fa_outputs.front().first) {
                            const string& fa = fa_outputs.front().second;
                            if (!fa.empty() && gzwrite(o_gzfasta_f, fa.c_str(), fa.length()) <= 0)
                                throw std::ios::failure("gzwrite");
                            fa_outputs.pop_front();
                            ++next_fa_to_write;
                        }
                    }
                    t.writing_fa.stop();

                    // Write the VCF output.
                    t.writing_vcf.restart();
                    #pragma omp critical(write_vcf)
                    {
                        for (size_t i=next_vcf_to_write+vcf_outputs.size(); i<=loc_i; ++i)
                            vcf_outputs.push_back( {false, string()} );
                        vcf_outputs[loc_i - next_vcf_to_write] = {true, move(loc_proc.vcf_out()) };

                        if (vcf_outputs.front().first) {
                            ++n_writes;
                            if (vcf_outputs.size() > max_size_before_write)
                                max_size_before_write = vcf_outputs.size();
                            t_writing_vcf.restart();
                            do {
                                o_vcf_f->file() << vcf_outputs.front().second;
                                vcf_outputs.pop_front();
                                ++next_vcf_to_write;
                                ++progress;
                            } while (!vcf_outputs.empty() && vcf_outputs.front().first);
                            t_writing_vcf.stop();
                        }
                    }
                    t.writing_vcf.stop();

                    // Write the detailed output.
                    if (detailed_output) {
                        t.writing_details.restart();
                        #pragma omp critical(write_details)
                        {
                            for (size_t i=next_det_to_write+det_outputs.size(); i<=loc_i; ++i)
                                det_outputs.push_back( {false, string()} );
                            det_outputs[loc_i - next_det_to_write] = {true, move(loc_proc.details_out())};

                            while (!det_outputs.empty() && det_outputs.front().first) {
                                *o_details_f << det_outputs.front().second;
                                det_outputs.pop_front();
                                ++next_det_to_write;
                            }
                        }
                        t.writing_details.stop();
                    }
                }} catch (exception& e) {
                    #pragma omp critical(exc)
                    omp_return = stacks_handle_exceptions(e);
                }
            }

            // Tally the thread statistics.
            #pragma omp critical(stats)
            {
                ctg_stats += loc_proc.ctg_stats();
                gt_stats  += loc_proc.gt_stats();
                t_threads_totals += loc_proc.timers();
            }
        }
        if (omp_return != 0)
            return omp_return;
        t_parallel.stop();
        progress.done();
        cout << '\n';
        assert(gt_stats.n_genotyped_loci == ctg_stats.n_nonempty_loci);

        // Report assembly statistics.
        if (ctg_stats.n_loci_w_pe_reads > 0) {
            assert(!ignore_pe_reads);
            size_t no_pe   = ctg_stats.n_loci_no_pe_reads() + ctg_stats.n_loci_almost_no_pe_reads;
            size_t pe_ndag  = ctg_stats.n_loci_pe_graph_not_dag;
            size_t pe_ctg = ctg_stats.n_loci_ctg();
            auto pct = [&ctg_stats](size_t n) { return as_percentage((double) n / ctg_stats.n_nonempty_loci); };

            ostream os (cout.rdbuf());
            os << std::fixed << std::setprecision(1);
            os << "Attempted to assemble and align paired-end reads for " << ctg_stats.n_nonempty_loci << " loci:\n"
               << "  " << no_pe << " loci had no or almost no paired-end reads (" << pct(no_pe) << ");\n"
               << "  " << pe_ndag << " loci had paired-end reads that couldn't be assembled into a contig ("
               << pct(pe_ndag) << ");\n"
               << "  For the remaining " << pe_ctg << " loci (" << pct(pe_ctg) << "), a paired-end contig was assembled;\n"
               << "  Average contig size was " << ctg_stats.ctg_avg_length() << " bp;\n"
               << "  Out of " << ctg_stats.n_tot_reads << " paired-end reads in these loci (mean "
               << (double) ctg_stats.n_aln_reads / pe_ctg << " reads per locus),\n"
               << "    " << ctg_stats.n_aln_reads << " were successfuly aligned ("
               << as_percentage((double) ctg_stats.n_aln_reads / ctg_stats.n_tot_reads) << ");\n"
               << "  " << ctg_stats.n_overlaps << " paired-end contigs overlapped the forward region ("
               << as_percentage((double) ctg_stats.n_overlaps / ctg_stats.n_loci_ctg()) << "; mean overlap: "
               << ctg_stats.mean_olap_length() << "bp).\n"
               << "\n";
        }

    } else if (input_type == GStacksInputT::refbased_popmap || input_type == GStacksInputT::refbased_list) {

        cout << "Processing all loci...\n" << flush;
        ProgressMeter progress (cout, false, 1000);
        bool eof = false;

        t_parallel.restart();
        #pragma omp parallel
        try {
            LocusProcessor loc_proc;
            Timers& t = loc_proc.timers();
            CLocAlnSet aln_loc;

            bool thread_eof;
            #pragma omp atomic read
            thread_eof = eof;
            while(!thread_eof && omp_return == 0) {
                // Read a locus from the BAM file.
                t.reading.restart();
                #pragma omp critical(read)
                try {
                    if (!(thread_eof = eof) && omp_return == 0)
                        thread_eof = eof = !bam_cloc_builder->build_one_locus(aln_loc);
                } catch (exception& e) {
                    omp_return = stacks_handle_exceptions(e);
                }
                t.reading.stop();
                if (thread_eof || omp_return != 0)
                    break;

                // Process it.
                t.processing.restart();
                assert(aln_loc.id() >= 1);
                size_t loc_i = aln_loc.id() - 1;
                aln_loc.merge_paired_reads();
                loc_proc.process(aln_loc);
                t.processing.stop();

                // Write the FASTA output.
                t.writing_fa.restart();
                #pragma omp critical(write_fa)
                {
                    for (size_t i=next_fa_to_write+fa_outputs.size(); i<=loc_i; ++i)
                        fa_outputs.push_back( {false, string()} );
                    fa_outputs[loc_i - next_fa_to_write] = {true, move(loc_proc.fasta_out())};

                    while (!fa_outputs.empty() && fa_outputs.front().first) {
                        const string& fa = fa_outputs.front().second;
                        if (!fa.empty() && gzwrite(o_gzfasta_f, fa.c_str(), fa.length()) <= 0)
                            throw std::ios::failure("gzwrite");
                        fa_outputs.pop_front();
                        ++next_fa_to_write;
                    }
                }
                t.writing_fa.stop();

                // Write the VCF output.
                t.writing_vcf.restart();
                #pragma omp critical(write_vcf)
                {
                    for (size_t i=next_vcf_to_write+vcf_outputs.size(); i<=loc_i; ++i)
                        vcf_outputs.push_back( {false, string()} );
                    vcf_outputs[loc_i - next_vcf_to_write] = {true, move(loc_proc.vcf_out()) };

                    if (vcf_outputs.front().first) {
                        ++n_writes;
                        if (vcf_outputs.size() > max_size_before_write)
                            max_size_before_write = vcf_outputs.size();
                        t_writing_vcf.restart();
                        do {
                            o_vcf_f->file() << vcf_outputs.front().second;
                            vcf_outputs.pop_front();
                            ++next_vcf_to_write;
                            ++progress;
                        } while (!vcf_outputs.empty() && vcf_outputs.front().first);
                        t_writing_vcf.stop();
                    }
                }
                t.writing_vcf.stop();

                // Write the detailed output.
                if (detailed_output) {
                    t.writing_details.restart();
                    #pragma omp critical(write_details)
                    {
                        for (size_t i=next_det_to_write+det_outputs.size(); i<=loc_i; ++i)
                            det_outputs.push_back( {false, string()} );
                        det_outputs[loc_i - next_det_to_write] = {true, move(loc_proc.details_out())};

                        while (!det_outputs.empty() && det_outputs.front().first) {
                            *o_details_f << det_outputs.front().second;
                            det_outputs.pop_front();
                            ++next_det_to_write;
                        }
                    }
                    t.writing_details.stop();
                }

                #pragma omp atomic read
                thread_eof = eof;
            }

            gt_stats += loc_proc.gt_stats();
            t_threads_totals += loc_proc.timers();
        } catch (exception& e) {
            #pragma omp critical (exc)
            omp_return = stacks_handle_exceptions(e);
        }
        if (omp_return != 0)
            return omp_return;
        t_parallel.stop();
        progress.done();
        cout << '\n';

        // Report statistics on the input BAM.
        const BamCLocBuilder::BamStats& bam_stats = bam_cloc_builder->stats();
        cout << "Read " << bam_stats.n_records << " BAM records:\n"
             << "  kept " << bam_stats.n_primary_kept() << " primary alignments\n"
             << "  skipped " << bam_stats.n_primary_mapq << " primary alignments with insufficient mapping qualities ("
             << as_percentage((double)bam_stats.n_primary_mapq / bam_stats.n_primary) << ")\n"
             << "  skipped " << bam_stats.n_primary_softclipped << " excessively soft-clipped primary alignments ("
             << as_percentage((double)bam_stats.n_primary_softclipped / bam_stats.n_primary) << ")\n"
             << "  Skipped " << bam_stats.n_secondary << " secondary alignments\n"
             << "  Skipped " << bam_stats.n_supplementary << " supplementary alignments\n"
             << "  Skipped " << bam_stats.n_unmapped << " unmapped reads\n"
             << "\n"
             ;

    } else DOES_NOT_HAPPEN; //input_type

    // Report statistics on genotyping and haplotyping.
    {
        size_t tot = gt_stats.n_genotyped_loci;
        size_t ph  = gt_stats.n_loci_phasing_issues();
        auto   pct = [tot](size_t n) { return as_percentage((double) n / tot); };

        cout << "Built and genotyped " << tot << " loci:\n"
             << "  (All loci are always conserved)\n"
             << "  one or more samples were excluded in " << ph << " loci (" << pct(ph) << ") because of phasing issues\n\n";

        logger->l << "BEGIN badly_phased\n"
                  << "n_tot_samples\tn_bad_samples\tn_loci\n";
        for (auto& elem : gt_stats.n_badly_phased_samples)
            logger->l << elem.first.second << '\t' << elem.first.first << '\t' << elem.second << '\n';
        logger->l << "END badly_phased\n\n";
    }

    #ifdef DEBUG
    if (typeid(*model) == typeid(MarukiLowModel&))
        cerr << "DEBUG: marukilow: " << ((const MarukiLowModel&)*model).n_underflows() << " underflows occurred.\n\n";
    #endif

    //
    // Report clockings.
    //
    {
        double ll  = t_parallel.elapsed();
        double v   = t_writing_vcf.elapsed();

        double r   = t_threads_totals.reading.elapsed() / num_threads;
        double p   = t_threads_totals.processing.elapsed() / num_threads;
        double w_f = t_threads_totals.writing_fa.elapsed() / num_threads;
        double w_v = t_threads_totals.writing_vcf.elapsed() / num_threads;
        double w_d = t_threads_totals.writing_details.elapsed() / num_threads;

        double a   = t_threads_totals.assembling.elapsed() / num_threads;
        double o   = t_threads_totals.olap_aligning.elapsed() / num_threads;
        double g   = t_threads_totals.geno_haplotyping.elapsed() / num_threads;
        double b_v = t_threads_totals.building_vcf.elapsed() / num_threads;

        double c = t_parallel.consumed()
                 + t_threads_totals.reading.consumed() / num_threads
                 + t_threads_totals.processing.consumed() / num_threads
                 + t_threads_totals.writing_fa.consumed() / num_threads
                 + t_threads_totals.writing_vcf.consumed() / num_threads
                 + t_threads_totals.writing_details.consumed() / num_threads
                 + t_threads_totals.assembling.consumed() / num_threads
                 + t_threads_totals.olap_aligning.consumed() / num_threads
                 + t_threads_totals.geno_haplotyping.consumed() / num_threads
                 + t_threads_totals.building_vcf.consumed() / num_threads
                 + t_writing_vcf.consumed()
                 ;

        ostream os (logger->l.rdbuf());
        os << std::fixed << std::setprecision(2);
        os << "BEGIN clockings\n"
           << "Num. threads: " << num_threads << "\n"
           << "Parallel time: " << ll << "\n"
           << "Average thread time spent:\n"
           << std::setw(8) << r  << "  reading (" << as_percentage(r / ll) << ")\n"
           << std::setw(8) << p << "  processing (" << as_percentage(p / ll) << ")\n";
        if (a != 0.0)
            // De novo mode & paired-ends.
            os << std::setw(16) << a << "  assembling (" << as_percentage(a / ll) << ")\n"
               << std::setw(16) << o << "  aligning/overlapping (" << as_percentage(o / ll) << ")\n";
        os << std::setw(16) << g << "  genotyping/haplotyping (" << as_percentage(g / ll) << ")\n"
           << std::setw(16) << b_v << "  building_vcf (" << as_percentage(b_v / ll) << ")\n"
           << std::setw(8) << w_f << "  writing_fa (" << as_percentage(w_f / ll) << ")\n"
           << std::setw(8) << w_v << "  writing_vcf (" << as_percentage(w_v / ll) << ")\n";
        if (detailed_output)
            os << std::setw(8) << w_d << "  writing_details (" << as_percentage(w_d / ll) << ")\n";
        os << std::setw(8) << c << "  clocking (" << as_percentage(c / ll) << ")\n"
           << "Total time spent writing vcf: " << v << " (" << as_percentage(v / ll) << ")\n"
           << "VCFwrite block size: mean=" << (double) gt_stats.n_genotyped_loci / n_writes
               << "(n=" << n_writes << "); max=" << max_size_before_write << "\n"
           << "END clockings\n\n"
           ;
    }

    //
    // Cleanup & return.
    //
    gzclose(o_gzfasta_f);
    model.reset();
    if (dbg_write_hapgraphs)
        o_hapgraphs_f << "}\n";
    cout << "gstacks is done.\n";
    return 0;

} catch (exception& e) {
    model.reset();
    return stacks_handle_exceptions(e);
}
}


//
// SnpAlleleCooccurrenceCounter
// ============================
//

const size_t& SnpAlleleCooccurrenceCounter::at(size_t snp_i1, Nt2 snp1_allele, size_t snp_i2, Nt2 snp2_allele) const {
    assert(snp_i1 < snp_i2);
    return cooccurences_[snp_i1*n_snps_+snp_i2][size_t(snp1_allele)][size_t(snp2_allele)];
}

void SnpAlleleCooccurrenceCounter::clear() {
    for(size_t i=0; i<n_snps_; ++i)
        for(size_t j=i+1; j<n_snps_; ++j)
            cooccurences_[i*n_snps_+j] = array<array<size_t,4>,4>();
}

//
// GenotypeStats & ContigStats
// ===============
//

GenotypeStats& GenotypeStats::operator+= (const GenotypeStats& other) {
    this->n_genotyped_loci           += other.n_genotyped_loci;
    for (auto count : other.n_badly_phased_samples)
        this->n_badly_phased_samples[count.first] += count.second;

    return *this;
}

ContigStats& ContigStats::operator+= (const ContigStats& other) {
    this->n_nonempty_loci           += other.n_nonempty_loci;
    this->n_loci_w_pe_reads         += other.n_loci_w_pe_reads;
    this->n_loci_almost_no_pe_reads += other.n_loci_almost_no_pe_reads;
    this->n_loci_pe_graph_not_dag   += other.n_loci_pe_graph_not_dag;
    this->length_ctg_tot            += other.length_ctg_tot;
    this->n_aln_reads               += other.n_aln_reads;
    this->n_tot_reads               += other.n_tot_reads;
    this->n_overlaps                += other.n_overlaps;
    this->length_overlap_tot        += other.length_overlap_tot;

    return *this;
}

//
// LocData
// =======
//
void LocData::clear() {
    id = -1;
    pos.clear();
    mpopi = NULL;
    ctg_status = unknown;
    o_vcf.clear();
    o_fa.clear();
    o_details.clear();
    details_ss.clear();
    details_ss.str(string());
}

//
// LocusProcessor
// ==============
//

void
LocusProcessor::process(CLocReadSet& loc)
{
    if (loc.reads().empty())
        return;
    ++ctg_stats_.n_nonempty_loci;

    this->loc_.clear();
    this->loc_.id  =  loc.id();
    this->loc_.pos =  loc.pos();
    this->loc_.mpopi   = &loc.mpopi();
    if (detailed_output)
        loc_.details_ss << "BEGIN locus " << loc_.id << "\n";

    //
    // Build the alignment matrix.
    //
    CLocAlnSet aln_loc;
    aln_loc.reinit(loc_.id, loc_.pos, loc_.mpopi);

    if (dbg_true_alns) {
        from_true_alignments(aln_loc, move(loc), true);
    } else if (dbg_true_reference) {
        using_true_reference(aln_loc, move(loc));
    } else {
        //
        // Transfer the already aligned foward-reads.
        //
        aln_loc.ref(DNASeq4(loc.reads().at(0).seq.length())); // Just N's.
        for (SRead& r : loc.reads())
            aln_loc.add(SAlnRead(Read(move(r.seq), move(r.name)), {{'M',r.seq.length()}}, r.sample));
        aln_loc.recompute_consensus();

        //
        // Sometimes the consensus in the catalog.tags.tsv file extends further
        // than any of the validly matching loci. In this case there's no
        // coverage at the end of the 'forward' alignment matrix, but the reads
        // are still N-padded. If this happens, hard-clip the Ns.
        //
        if (*--aln_loc.ref().end() == Nt4::n)
            aln_loc.hard_clip_right_Ns();

        //
        // Process the paired-end reads, if any.
        //
        do { // Avoiding nested IFs.
            if (ignore_pe_reads)
                break;
            if (loc.pe_reads().empty())
                break;
            ++this->ctg_stats_.n_loci_w_pe_reads;
            if (detailed_output) {
                loc_.details_ss << "BEGIN pe_reads\n";
                for (const SRead& r : loc.pe_reads())
                    loc_.details_ss << "pe_read"
                                    << '\t' << r.name
                                    << '\t' << loc_.mpopi->samples()[r.sample].name
                                    << '\t' << r.seq << '\n';
                loc_.details_ss << "END pe_reads\n";
            }

            // Assemble a contig.
            timers_.assembling.restart();
            vector<const DNASeq4*> seqs_to_assemble;
            for (const Read& r : loc.pe_reads())
                seqs_to_assemble.push_back(&r.seq);
            DNASeq4 ctg = DNASeq4(assemble_contig(seqs_to_assemble));
            timers_.assembling.stop();
            if (ctg.empty())
                break;
            ctg_stats_.length_ctg_tot += ctg.length();

            //
            // Build a SuffixTree of the contig.
            //
            timers_.olap_aligning.restart();
            SuffixTree *stree   = new SuffixTree(ctg);
            stree->build_tree();

            //
            // Initialize a gapped alignment obect for use in SE/PE contig overlap and read alignmnets.
            //
            GappedAln aligner;
            AlignRes  aln_res;

            //
            // Determine if there is overlap -- and how much -- between the SE and PE contigs.
            //   We will query the PE contig suffix tree using the SE consensus sequence.
            //
            int    overlap;
            string overlap_cigar;
            if (dbg_no_overlaps)
                overlap = 0;
            else
                overlap = this->find_locus_overlap(stree, &aligner, aln_loc.ref(), overlap_cigar);

            if (overlap > 0) {
                this->loc_.ctg_status = LocData::overlapped;
                this->ctg_stats_.n_overlaps++;
                this->ctg_stats_.length_overlap_tot += overlap;
            } else {
                this->loc_.ctg_status = LocData::separate;
            }

            if (detailed_output)
                loc_.details_ss
                        << "BEGIN contig\n"
                        << "pe_contig\t"     << ctg << '\n'
                        << "fw_consensus\t"  << aln_loc.ref() << '\n'
                        << "overlap\t"       << overlap << '\t' << overlap_cigar << "\n"
                        ;

            //
            // Extend the reference sequence & recompute the suffix tree.
            //
            CLocAlnSet tmp_loc;
            tmp_loc.reinit(loc_.id, loc_.pos, loc_.mpopi);
            tmp_loc.ref(move(ctg));

            aln_loc = CLocAlnSet::juxtapose(move(aln_loc), move(tmp_loc), (overlap > 0 ? -long(overlap) : +10));
            if (detailed_output)
                loc_.details_ss << "merger\t" << aln_loc.ref() << "\nEND contig\n";

            delete stree;
            stree   = new SuffixTree(aln_loc.ref());
            stree->build_tree();

            //
            // Align the paired-end reads.
            //
            this->ctg_stats_.n_tot_reads += loc.pe_reads().size();
            if (detailed_output)
                loc_.details_ss << "BEGIN pe_alns\n";
            for (SRead& r : loc.pe_reads()) {
                if (add_read_to_aln(aln_loc, aln_res, move(r), &aligner, stree)) {
                    this->ctg_stats_.n_aln_reads++;
                    if (detailed_output)
                        loc_.details_ss << "pe_aln_local"
                                        << '\t' << aln_loc.reads().back().name
                                        << '\t' << aln_res.subj_pos + 1 << ':' << aln_res.cigar
                                        << '\n';
                } else {
                    if (detailed_output)
                        loc_.details_ss << "pe_aln_fail\t" << aln_loc.reads().back().name << '\n';
                }
            }
            if (detailed_output)
                loc_.details_ss << "END pe_alns\n";
            delete stree;
            timers_.olap_aligning.stop();

            aln_loc.merge_paired_reads();

        } while (false);
    }
    loc.clear();

    process(aln_loc);
}

void
LocusProcessor::process(CLocAlnSet& aln_loc)
{
    assert(!aln_loc.reads().empty());
    ++gt_stats_.n_genotyped_loci;
    if (input_type == GStacksInputT::denovo_popmap || input_type == GStacksInputT::denovo_merger) {
        // Called from process(CLocReadSet&).
        assert(this->loc_.id == aln_loc.id());
    } else {
        this->loc_.clear();
        this->loc_.id = aln_loc.id();
        this->loc_.pos = aln_loc.pos();
        this->loc_.mpopi = &aln_loc.mpopi();
        if (detailed_output)
            this->loc_.details_ss << "BEGIN locus " << loc_.id << "\n";
    }

    aln_loc.recompute_consensus();

    if (detailed_output) {
        loc_.details_ss << "BEGIN aln_matrix\n";
        for (const SAlnRead& read : aln_loc.reads())
            loc_.details_ss << read.name << '\t' << loc_.mpopi->samples()[read.sample].name << '\t' << read.cigar << '\n';
        loc_.details_ss << "END aln_matrix\n";
    }

    if (dbg_write_alns)
        #pragma omp critical
        o_aln_f << "BEGIN " << loc_.id << "\n"
                << aln_loc
                << "\nEND " << loc_.id << "\n";

    //
    // Call SNPs. Determine the consensus sequence.
    //
    timers_.geno_haplotyping.restart();
    vector<SiteCounts> depths;
    vector<SiteCall> calls;
    depths.reserve(aln_loc.ref().length());
    calls.reserve(aln_loc.ref().length());
    DNASeq4 new_consensus = aln_loc.ref();
    for (CLocAlnSet::site_iterator site (aln_loc); bool(site); ++site) {
        depths.push_back(site.counts());
        calls.push_back(model->call(depths.back()));
        if (!calls.back().alleles().empty())
            new_consensus.set(site.col(), calls.back().most_frequent_allele());
        else
            // For the high/low Maruki" models this actually only happens when
            // there is no coverage; for the Hohenlohe model it may also happen
            // when there isn't a single significant call.
            // Keep the already computed majority-rule nucleotide/don't do anything.
            ;
    }
    aln_loc.ref(move(new_consensus));

    // Call haplotypes.
    vector<map<size_t,PhasedHet>> phase_data;
    if (!dbg_no_haplotypes) {
        set<size_t> inconsistent_samples;
        phase_data = phase_hets(calls, aln_loc, inconsistent_samples);
        if (!inconsistent_samples.empty())
            ++gt_stats_.n_badly_phased_samples[ {inconsistent_samples.size(), aln_loc.n_samples()} ];
    }
    timers_.geno_haplotyping.stop();

    write_one_locus(aln_loc, depths, calls, phase_data);

    if (detailed_output) {
        loc_.details_ss << "END locus " << loc_.id << "\n";
        loc_.o_details = loc_.details_ss.str();
    }
}

bool
LocusProcessor::align_reads_to_contig(SuffixTree *st, GappedAln *g_aln, DNASeq4 &enc_query, AlignRes &aln_res) const
{
    vector<STAln> alns, final_alns;
    vector<pair<size_t, size_t> > step_alns;
    string      query  = enc_query.str();
    const char *q      = query.c_str();
    const char *q_stop = q + query.length();
    size_t      q_pos  = 0;
    size_t      id     = 0;
    char c[id_len];

    do {
        step_alns.clear();

        q_pos = q - query.c_str();

        st->align(q, step_alns);

        if (step_alns.size() == 0 || step_alns.size() > max_fragment_alns) {
            q++;
        } else {
            for (uint i = 0; i < step_alns.size(); i++)
                alns.push_back(STAln(id, q_pos, step_alns[i].first, step_alns[i].second));
            q += step_alns[0].second + 1;
            id++;
        }
    } while (q < q_stop);

    //
    // No alignments to the suffix tree were found.
    //
    if (alns.size() == 0)
        return false;

    //
    // Perfect alignmnet to the suffix tree. Return result.
    //
    if (alns.size() == 1 && alns[0].aln_len == query.length()) {
        snprintf(c, id_len, "%luM", query.length());
        aln_res.cigar    = c;
        aln_res.subj_pos = alns[0].subj_pos;
        return true;
    }

    //
    // Find a consistent set of suffix tree alignments, ordered in a directed, acyclic grapgh (DAG).
    //
    this->suffix_tree_hits_to_dag(query.length(), alns, final_alns);

    g_aln->init(query.length(), st->seq_len(), true);
    if (g_aln->align_constrained(query, st->seq_str(), final_alns)) {
        aln_res = g_aln->result();
    }

    return true;
}

int
LocusProcessor::suffix_tree_hits_to_dag(size_t query_len, vector<STAln> &alns, vector<STAln> &final_alns) const
{
    //
    // 1. Sort the alignment fragments so they are primarily ordered by subject position, secondarily by query position.
    //
    sort(alns.begin(), alns.end(), [](STAln a, STAln b)
         {
            if (a.subj_pos == b.subj_pos)
                return a.query_pos < b.query_pos;
            else
                return a.subj_pos < b.subj_pos;
         });

    int    gap_len, end_pos;
    double score, scale;
    bool   term;
    vector<size_t> term_nodes;
    //
    // 2. Traverse the list of fragments and add links to reachable nodes, recording the score.
    //    2.1. A successor node is reachable from the predecessor node if both the query and subject
    //         positions are advanced relative to the predecessor, and there is no sequence overlap
    //         between the two nodes.
    //    2.2. We weight the score by the number of 'gap' nucleotides between the two fragments to
    //         penalize fragments that are further away.
    //    2.3  Record the maximal link in the graph.
    //    2.4. Determine if a node is a terminal node.
    //
    for (uint i = 0; i < alns.size(); i++) {
        end_pos = alns[i].subj_pos + alns[i].aln_len - 1;
        term    = true;

        assert(end_pos > 0);

        for (uint j = i + 1; j < alns.size(); j++) {
            if (alns[j].query_pos > alns[i].query_pos &&
                alns[j].subj_pos  > alns[i].subj_pos  &&
                alns[j].subj_pos  > (uint) end_pos) {

                term    = false;
                gap_len = alns[j].subj_pos - end_pos - 1;
                scale   = gap_len > (int) query_len ? 1 : alns[i].aln_len * ((double) gap_len / (double) query_len);
                score   = alns[i].aln_len - scale; // Raw score.
                score  += alns[i].max._score;      // Score adjusted for the highest incoming node.

                assert(gap_len >= 0);
                assert(score   >= 0);

                alns[j].links.push_back(STLink(i, score));
                if (alns[j].max._index == j || score > alns[j].max._score)
                    alns[j].max = STLink(i, score);
            }
        }
        if (term)
            term_nodes.push_back(i);
    }

    //
    // 3. Find the terminal node with the highest score. 
    //
    double max_score = alns[term_nodes.front()].max._score;
    size_t max_index = term_nodes.front(); 
    for (uint i = 1; i < term_nodes.size(); i++) {
        if (alns[term_nodes[i]].max._score > max_score) {
            max_index = term_nodes[i];
            max_score = alns[term_nodes[i]].max._score;
        }
    }

    //
    // 4. Backtrack to get the optimal path.
    //
    vector<size_t> optimal;

    if (max_score == 0 && max_index == 0) {
        //
        // None of the fragments were connected in the DAG, select the largest fragment.
        //
        max_score = alns[term_nodes.front()].aln_len;
        max_index = term_nodes.front();
        for (uint i = 1; i < term_nodes.size(); i++)
            if (alns[i].aln_len > max_score) {
                max_index = term_nodes[i];
                max_score = alns[term_nodes[i]].aln_len;
            }
        optimal.push_back(max_index);
                
    } else {

        uint n = max_index;
        while (alns[n].links.size() > 0) {
            optimal.push_back(n);
            n = alns[n].max._index;
        }
        optimal.push_back(n);
    }    
    assert(optimal.size() > 0);

    final_alns.clear();
    for (int i = optimal.size() - 1; i >= 0; i--)
        final_alns.push_back(alns[optimal[i]]);

    return 0;
}

size_t
LocusProcessor::find_locus_overlap(SuffixTree *stree, GappedAln *g_aln, const DNASeq4 &se_consensus, string &overlap_cigar) const
{
    vector<STAln> alns, final_alns;
    vector<pair<size_t, size_t> > step_alns;
    AlignRes aln_res;

    string      query  = se_consensus.str();
    const char *q      = query.c_str();
    const char *q_stop = q + query.length();
    size_t      q_pos  = 0;
    size_t      id     = 1;
    const char *p;

    do {
        step_alns.clear();

        q_pos = q - query.c_str();

        stree->align(q, step_alns);

        if (step_alns.size() == 0 || step_alns.size() > max_fragment_alns) {
            q++;
        } else {
            for (uint i = 0; i < step_alns.size(); i++)
                alns.push_back(STAln(id, q_pos, step_alns[i].first, step_alns[i].second));
            q += step_alns[0].second + 1;
            id++;
        }
    } while (q < q_stop);

    if (alns.size() == 0) {
        //
        // If no alignments have been found, search the tails of the query and subject for any overlap
        // that is too small to be picked up by the SuffixTree.
        //
        int min_olap = (int) stree->min_aln() > min_se_pe_overlap ? stree->min_aln() : min_se_pe_overlap;
        string pe_ctg = stree->seq_str().substr(0, min_olap);
        p = pe_ctg.c_str();
        q = query.c_str() + (query.length() - min_olap);

        while (min_olap >= min_se_pe_overlap && *q != '\0') {
            if (strncmp(q, p, min_olap) == 0)
                return min_olap;
            min_olap--;
            q++;
        }

        return 0;
    }

    //
    // Perfect alignmnet to the suffix tree that occupies the end of the
    // query and the beginning of the subject. Return result.
    //
    if (alns.size() == 1) {
        size_t query_stop = alns.front().query_pos + alns.front().aln_len - 1;

        if (query_stop == (query.length() - 1) && alns.front().subj_pos == 0)
            return alns.front().aln_len;
    }

    //
    // Otherwise, find a consistent set of suffix tree alignments, ordered in a directed, acyclic grapgh (DAG).
    //
    this->suffix_tree_hits_to_dag(query.length(), alns, final_alns);

    //
    // Create a gapped alignment to determine the exact overlap.
    //
    g_aln->init(query.length(), stree->seq_len(), true);
    if (g_aln->align_constrained(query, stree->seq_str(), final_alns)) {
        aln_res       = g_aln->result();
        overlap_cigar = aln_res.cigar;
    }

    Cigar cigar;
    parse_cigar(aln_res.cigar.c_str(), cigar, true);

    //
    // If the first element in the CIGAR is a soft-masked region, ignore it.
    //
    size_t offset = 0;
    if (cigar.front().first != 'S') {
        offset = aln_res.subj_pos;
    }

    return offset + cigar_length_ref(cigar);
}

string
LocusProcessor::assemble_contig(const vector<const DNASeq4*>& seqs)
{
    Graph graph (km_length);

    graph.rebuild(seqs, min_km_count);
    if (graph.empty()) {
        ++ctg_stats_.n_loci_almost_no_pe_reads;
        return string();
    }

    if (dbg_write_gfa) {
        graph.dump_gfa(o_prefix + "." + to_string(loc_.id) + ".spaths.gfa");
        graph.dump_gfa(o_prefix + "." + to_string(loc_.id) + ".nodes.gfa", true);
    }

    vector<const SPath*> best_path;
    if (!graph.find_best_path(best_path)) {
        // Not a DAG.
        ++ctg_stats_.n_loci_pe_graph_not_dag;
        return string();
    }

    return SPath::contig_str(best_path.begin(), best_path.end(), km_length);
}

bool LocusProcessor::add_read_to_aln(
        CLocAlnSet& aln_loc,
        AlignRes& aln_res,
        SRead&& r,
        GappedAln* aligner,
        SuffixTree* stree
) const {

    if (!this->align_reads_to_contig(stree, aligner, r.seq, aln_res))
        return false;

    if (aln_res.pct_id < min_aln_cov)
        return false;

    Cigar cigar;
    parse_cigar(aln_res.cigar.c_str(), cigar);
    simplify_cigar_to_MDI(cigar);
    cigar_extend_left(cigar, aln_res.subj_pos);
    assert(cigar_length_ref(cigar) <= aln_loc.ref().length());
    cigar_extend_right(cigar, aln_loc.ref().length() - cigar_length_ref(cigar));

    aln_loc.add(SAlnRead(move((Read&)r), move(cigar), r.sample));
    return true;
}

vector<map<size_t,PhasedHet>> LocusProcessor::phase_hets (
        const vector<SiteCall>& calls,
        const CLocAlnSet& aln_loc,
        set<size_t>& inconsistent_samples
) const {
    inconsistent_samples.clear();

    vector<map<size_t,PhasedHet>> phased_samples (loc_.mpopi->samples().size());
    vector<size_t> snp_cols; // The SNPs of this locus.
    for (size_t i=0; i<aln_loc.ref().length(); ++i)
        if (calls[i].alleles().size() > 1)
            snp_cols.push_back(i);

    if (snp_cols.empty())
        return phased_samples;

    stringstream o_hapgraph_ss;
    if (dbg_write_hapgraphs) {
        o_hapgraph_ss << "\n"
                      << "subgraph cluster_loc" << loc_.id << " {\n"
                      << "\tlabel=\"locus " << loc_.id << "\";\n"
                      << "\t# snp columns: ";
        join(snp_cols, ',', o_hapgraph_ss);
        o_hapgraph_ss << "\n";
    }

    SnpAlleleCooccurrenceCounter cooccurrences (snp_cols.size());
    for (size_t sample=0; sample<loc_.mpopi->samples().size(); ++sample) {
        if (aln_loc.sample_reads(sample).empty())
            continue;

        //
        // Find heterozygous positions & check that we have >= 2 hets.
        //
        vector<size_t> het_snps;
        for(size_t snp_i=0; snp_i<snp_cols.size(); ++snp_i)
            if (calls[snp_cols[snp_i]].sample_calls()[sample].call() == snp_type_het)
                het_snps.push_back(snp_i);
        if (het_snps.size() == 0) {
            // Sample is homozygous.
            continue;
        } else if (het_snps.size() == 1) {
            // One heterozygous SNP; sample has trivial 1nt-long haplotypes.
            size_t col = snp_cols[het_snps[0]];
            const SampleCall& c = calls[col].sample_calls()[sample];
            phased_samples[sample].insert({col, {col, c.nt0(), c.nt1()}});
            continue;
        }

        vector<const SampleCall*> sample_het_calls;
        for (size_t het_i=0; het_i<het_snps.size(); ++het_i)
            sample_het_calls.push_back(&calls[snp_cols[het_snps[het_i]]].sample_calls()[sample]);

        //
        // Iterate over reads, record seen haplotypes (as pairwise cooccurrences).
        //
        count_pairwise_cooccurrences(cooccurrences, aln_loc, sample, snp_cols, het_snps, sample_het_calls);

        if (dbg_write_hapgraphs)
            write_sample_hapgraph(o_hapgraph_ss, sample, het_snps, snp_cols, sample_het_calls, cooccurrences);

        //
        // Assemble phase sets.
        //
        vector<PhaseSet> phase_sets;
        if (!assemble_phase_sets(phase_sets, het_snps, sample_het_calls, cooccurrences)) {
            inconsistent_samples.insert(sample);
            continue;
        }

        //
        // Record phase sets.
        //
        map<size_t,PhasedHet>& sample_phase_data = phased_samples[sample];
        for (const PhaseSet& ps : phase_sets) {
            assert(ps.size() == het_snps.size());
            size_t phase_set_id = SIZE_MAX;
            for (size_t het_i=0; het_i<het_snps.size(); ++het_i) {
                if (ps.het(het_i).left_allele == Nt4::n)
                    continue;

                size_t col = snp_cols[het_snps[het_i]];
                if (phase_set_id == SIZE_MAX)
                    phase_set_id = col; // Recommended value, c.f. VCF specification.
                PhasedHet ph = ps.het(het_i);
                ph.phase_set = phase_set_id;
                sample_phase_data[col] = ph;
            }
        }

        // Record singleton nodes (that are implicit in our representation).
        for (size_t het_i=0; het_i<het_snps.size(); ++het_i) {
            size_t col = snp_cols[het_snps[het_i]];
            if (!sample_phase_data.count(col)) {
                array<Nt2,2> alleles = sample_het_calls[het_i]->nts();
                sample_phase_data[col] = PhasedHet({col, alleles[0], alleles[1]});
            }
        }
        assert(!sample_phase_data.empty());

        //
        // Remove all phase sets but the largest one.
        //
        map<size_t,size_t> phase_set_sizes;
        for (auto& phasedhet : sample_phase_data)
            ++phase_set_sizes[phasedhet.second.phase_set];
        auto best_ps = phase_set_sizes.begin();
        for (auto ps=++phase_set_sizes.begin(); ps!=phase_set_sizes.end(); ++ps)
            if (ps->second > best_ps->second)
                best_ps = ps;
        for (auto phasedhet=sample_phase_data.begin(); phasedhet!=sample_phase_data.end();) {
            if (phasedhet->second.phase_set == best_ps->first)
                ++phasedhet;
            else
                sample_phase_data.erase(phasedhet++);
        }
    }

    if (dbg_write_hapgraphs) {
        o_hapgraph_ss << "}\n";
        #pragma omp critical
        o_hapgraphs_f << o_hapgraph_ss.str();
    }

    return phased_samples;
}

void LocusProcessor::count_pairwise_cooccurrences(
        SnpAlleleCooccurrenceCounter& cooccurrences,
        const CLocAlnSet& aln_loc,
        size_t sample,
        const vector<size_t>& snp_cols,
        const vector<size_t>& het_snps,
        const vector<const SampleCall*>& sample_het_calls
        ) const {

    cooccurrences.clear();
    vector<Nt4> read_hap (het_snps.size());
    for (size_t read_i : aln_loc.sample_reads(sample)) {
        // Build the haplotype.
        size_t curr_col = 0;
        Alignment::iterator read_itr (aln_loc.reads()[read_i].aln());
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
                ++cooccurrences.at(het_snps[i], Nt2(nti), het_snps[j], Nt2(ntj));
            }
        }
    }
}

bool LocusProcessor::assemble_phase_sets(
        vector<PhaseSet>& phase_sets,
        const vector<size_t>& het_snps,
        const vector<const SampleCall*>& sample_het_calls,
        const SnpAlleleCooccurrenceCounter& cooccurrences
) const {

    static const size_t min_n_cooccurrences = 2;

    phase_sets.clear();

    // We keep track of which phase set each het is currently part of.
    vector<size_t> allele_to_ps (het_snps.size(), SIZE_MAX);

    for (size_t het_i=0; het_i<het_snps.size(); ++het_i) {
        size_t snp_i = het_snps[het_i];
        size_t& ps_i = allele_to_ps[het_i]; //(by reference; always up to date)
        array<Nt2,2> allelesi = sample_het_calls[het_i]->nts();
        for (size_t het_j=het_i+1; het_j<het_snps.size(); ++het_j) {
            size_t snp_j = het_snps[het_j];
            size_t& ps_j = allele_to_ps[het_j];
            array<Nt2,2> allelesj = sample_het_calls[het_j]->nts();

            for (Nt2 nti : allelesi) {
                for (Nt2 ntj : allelesj) {
                    size_t n = cooccurrences.at(snp_i, nti, snp_j, ntj);
                    if (n == 0)
                        continue;

                    if (n < min_n_cooccurrences)
                        // Discard low-coverage edges.
                        continue;

                    if (ps_i == SIZE_MAX && ps_j == SIZE_MAX) {
                        // Both nodes are singletons. Start a new phase set.
                        ps_i = phase_sets.size();
                        ps_j = phase_sets.size();
                        phase_sets.push_back(PhaseSet(het_snps.size()));
                        phase_sets.back().add_het(het_i, allelesi);
                        phase_sets.back().add_het(het_j, allelesj, ntj, het_i, nti);
                    } else if (ps_i == SIZE_MAX) {
                        // Add `het_i` to `ps_j`.
                        phase_sets[ps_j].add_het(het_i, allelesi, nti, het_j, ntj);
                        ps_i = ps_j;
                    } else if (ps_j == SIZE_MAX) {
                        // Add `het_j` to `ps_i`.
                        phase_sets[ps_i].add_het(het_j, allelesj, ntj, het_i, nti);
                        ps_j = ps_i;
                    } else if (ps_i != ps_j) {
                        // Merge `ps_j` into `ps_i`.
                        phase_sets[ps_i].merge_with(phase_sets[ps_j], het_i, nti, het_j, ntj);
                        phase_sets[ps_j].clear();
                        size_t merged_ps = ps_j;
                        for (size_t& allele_ps : allele_to_ps)
                            if (allele_ps == merged_ps)
                                allele_ps = ps_i;
                    } else {
                        assert(ps_i == ps_j);
                        // Check that the edge is consistent.
                        if (!phase_sets[ps_i].is_edge_consistent(het_i, nti, het_j, ntj))
                            return false;
                    }
                }
            }
        }
    }

    // Purge empty phase sets.
    phase_sets.erase(std::remove_if(
            phase_sets.begin(), phase_sets.end(),
            [] (const PhaseSet& ps) {return ps.empty();}
            ),phase_sets.end());
    return true;
}

void PhaseSet::add_het(size_t het_i, array<Nt2,2> alleles) {
    assert(phased_hets_.size() > het_i);
    assert(phased_hets_ == PhaseSet(phased_hets_.size()).phased_hets_); // Uninitialized.

    phased_hets_[het_i].left_allele  = Nt4(alleles[0]);
    phased_hets_[het_i].right_allele = Nt4(alleles[1]);
}

void PhaseSet::add_het(size_t het_i, array<Nt2,2> alleles, Nt2 nt_i, size_t het_j, Nt2 nt_j) {
    assert(phased_hets_.size() > std::max(het_i, het_j));
    assert(phased_hets_[het_i].left_allele == Nt4::n); // `i` node shouldn't exist yet.
    assert(phased_hets_[het_j].left_allele != Nt4::n); //

    bool crossed_edge =
            (nt_j == phased_hets_[het_j].left_allele)
            ^ (nt_i == alleles[0]);

    if (crossed_edge) {
        // Flip the alleles so that the all the edges within the phase set are parallel.
        phased_hets_[het_i].left_allele  = Nt4(alleles[1]);
        phased_hets_[het_i].right_allele = Nt4(alleles[0]);
    } else {
        phased_hets_[het_i].left_allele  = Nt4(alleles[0]);
        phased_hets_[het_i].right_allele = Nt4(alleles[1]);
    }
}

void PhaseSet::merge_with(const PhaseSet& other, size_t het_i, Nt2 nt_i, size_t het_j, Nt2 nt_j) {
    assert(phased_hets_.size() > std::max(het_i, het_j));
    assert(phased_hets_.size() == other.phased_hets_.size());
    assert(phased_hets_[het_i].left_allele != Nt4::n);
    assert(other.phased_hets_[het_j].left_allele != Nt4::n);

    bool crossed_edge =
            (nt_i == phased_hets_[het_i].left_allele)
            ^ (nt_j == other.phased_hets_[het_j].left_allele);

    for (size_t i=0; i<phased_hets_.size(); ++i) {
        if (other.phased_hets_[i].left_allele != Nt4::n) {
            assert(other.phased_hets_[i].right_allele != Nt4::n); // Both alleles should be set/not set together.
            assert(phased_hets_[i].left_allele == Nt4::n); // Phase sets should be non-overlapping.
            if (crossed_edge) {
                phased_hets_[i].left_allele  = other.phased_hets_[i].right_allele;
                phased_hets_[i].right_allele = other.phased_hets_[i].left_allele;
            } else {
                phased_hets_[i].left_allele  = other.phased_hets_[i].left_allele;
                phased_hets_[i].right_allele = other.phased_hets_[i].right_allele;
            }
        }
    }
}

bool PhaseSet::is_edge_consistent(size_t het_i, Nt2 nt_i, size_t het_j, Nt2 nt_j) const {
    assert(phased_hets_[het_i].left_allele != Nt4::n);
    assert(phased_hets_[het_j].left_allele != Nt4::n);

    bool crossed_edge =
            (nt_i == phased_hets_[het_i].left_allele)
            ^ (nt_j == phased_hets_[het_j].left_allele);

    // We build phase sets so that all inside edges are parallel
    // (by flipping the alleles when necessary), so the proposed
    // edge is inconsistent if it is crossed.
    return !crossed_edge;
}

void LocusProcessor::write_one_locus (
        const CLocAlnSet& aln_loc,
        const vector<SiteCounts>& depths,
        const vector<SiteCall>& calls,
        const vector<map<size_t,PhasedHet>>& phase_data
) {
    char loc_id[16];
    sprintf(loc_id, "%d", loc_.id);

    const DNASeq4& ref = aln_loc.ref();
    const MetaPopInfo& mpopi = aln_loc.mpopi();

    //
    // Vcf output.
    //
    timers_.building_vcf.restart();
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
                            std::sort(gt.begin(), gt.end()); // (Prevents '1/0'.)
                            genotype << gt[0] << '/' << gt[1];
                            genotype << ":.";
                        }
                    } else {
                        gt[0] = vcf_allele_indexes.at(scall.nt0());
                        gt[1] = vcf_allele_indexes.at(scall.nt1());
                        std::sort(gt.begin(), gt.end()); // (Prevents '1/0'.)
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
    loc_.o_vcf = vcf_records.str();
    timers_.building_vcf.stop();

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
    string& fa = loc_.o_fa;
    assert(fa.empty());
    fa += '>';
    fa += loc_id;
    if (!aln_loc.pos().empty()) {
        const PhyLoc& p = aln_loc.pos();
        char pos[16];
        sprintf(pos, "%u", p.bp+1);
        fa += " pos=";
        fa += p.chr();
        fa += ':';
        fa += pos;
        fa += ':';
        fa += (p.strand == strand_plus ? '+' : '-');
    }
    char n_spls[32];
    sprintf(n_spls, "%zu", n_remaining_samples);
    fa += " NS=";
    fa += n_spls;
    if (n_remaining_samples != samples_w_reads.size()) {
        assert(n_remaining_samples < samples_w_reads.size());
        sprintf(n_spls, "%zu", samples_w_reads.size() - n_remaining_samples);
        fa += " n_discarded_samples=";
        fa += n_spls;
    }
    switch (loc_.ctg_status) {
    case LocData::overlapped: fa += " contig=overlapped"; break;
    case LocData::separate:   fa += " contig=separate";   break;
    default: break;
    }
    fa += '\n';
    fa += ref.str();
    fa += '\n';
}

void LocusProcessor::write_sample_hapgraph(
        ostream& os,
        size_t sample,
        const vector<size_t>& het_snps,
        const vector<size_t>& snp_cols,
        const vector<const SampleCall*>& sample_het_calls,
        const SnpAlleleCooccurrenceCounter& cooccurrences
        ) const {

    auto nodeid = [this,&sample](size_t col, Nt2 allele)
            {return string("l")+to_string(loc_.id)+"s"+to_string(sample)+"c"+to_string(col)+char(allele);};

    // Initialize the subgraph.
    os << "\tsubgraph cluster_sample" << sample << " {\n"
                  << "\t\tlabel=\""
                  << "i" << sample << " '" << loc_.mpopi->samples()[sample].name << "'\\n"
                  << "\";\n"
                  << "\t\tstyle=dashed;\n"
                  << "\t\t# heterozygous columns: ";
    vector<size_t> het_cols;
    for (size_t snp_i : het_snps)
        het_cols.push_back(snp_cols[snp_i]);
    join(het_cols, ',', os);
    os << "\n";

    // Write the node labels.
    for (size_t i=0; i<het_snps.size(); ++i) {
        array<Nt2,2> alleles = sample_het_calls[i]->nts();
        size_t col = snp_cols[het_snps[i]];
        for (Nt2 allele : alleles)
            os << "\t\t" << nodeid(col, allele)
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
                    size_t n = cooccurrences.at(snp_i, nti, snp_j, ntj);
                    if (n == 0)
                        continue;
                    os << "\t\t" << nodeid(snp_cols[snp_i],nti)
                                  << " -- " << nodeid(snp_cols[snp_j],ntj) << " [";
                    if (n==1)
                        os << "style=dotted";
                    else
                        os << "label=\"" << n << "\",penwidth=" << n;
                    os << "];\n";
                }
            }
        }
    }
    os << "\t}\n";
}

//
// Debugging & benchmarking
// ========================
//

Cigar dbg_extract_cigar(const string& read_id) {
    static const char keyword1[] = "cig1=";
    static const char keyword2[] = "cig2=";
    static const size_t kw_len = sizeof(keyword1) - 1;
    static_assert(sizeof(keyword1) == sizeof(keyword2), "");

    Cigar cigar;

    bool paired_end = read_id.back() == '2' && read_id.at(read_id.length()-2) == '/';

    // Find the start.
    const char* cig_start = read_id.c_str();
    const char* kw = keyword1;
    if (paired_end)
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

    // For paired-end reads, reverse complement (i.e. just reverse) the cigar.
    if (paired_end)
        std::reverse(cigar.begin(), cigar.end());

    return cigar;
}

void from_true_alignments(CLocAlnSet& aln_loc, CLocReadSet&& loc, bool merge_reads) {
    //
    // Add forward reads.
    //
    for (SRead& r : loc.reads()) {
        Cigar cigar = dbg_extract_cigar(r.name);
        simplify_cigar_to_MDI(cigar);
        if (aln_loc.ref().empty()) {
            aln_loc.ref(DNASeq4(string(cigar_length_ref(cigar), 'N')));
        } else if (cigar_length_ref(cigar) != aln_loc.ref().length()) {
            cerr << "Error: DEBUG: ref-length of cig1 in read '" << r.name
                 << "' was expected to be " << aln_loc.ref().length()
                 << ", not " << cigar_length_ref(cigar) << ".\n";
            throw exception();
        }
        // Undo the alignments. We remove the Ns added by tsv2bam (that make
        // the sequence longer than the read), but not the Ns added by
        // ustacks/pstacks (as the sequence was trimmed accordingly).
        r.seq.shift_Ns_towards_the_end();
        r.seq.resize(cigar_length_query(cigar));
        aln_loc.add(SAlnRead(move((Read&)r), move(cigar), r.sample));
    }

    //
    // Add paired-end reads, if any.
    //
    if (!ignore_pe_reads && !loc.pe_reads().empty()) {
        for (SRead& r : loc.pe_reads()) {
            Cigar cigar = dbg_extract_cigar(r.name);
            simplify_cigar_to_MDI(cigar);
            if (cigar_length_ref(cigar) != aln_loc.ref().length()) {
                cerr << "Error: DEBUG: ref-length of cig2 in read '" << r.name
                     << "' was expected to be " << aln_loc.ref().length()
                     << ", not " << cigar_length_ref(cigar) << ".\n";
                throw exception();
            }
            aln_loc.add(SAlnRead(move((Read&)r), move(cigar), r.sample));
        }
        if (merge_reads)
            aln_loc.merge_paired_reads();
    }
}

void LocusProcessor::using_true_reference(CLocAlnSet& aln_loc, CLocReadSet&& loc) {
    // Reconstruct the reference sequence, based on true alignments.
    from_true_alignments(aln_loc, move(loc), false);
    DNASeq4 true_ref (aln_loc.ref().length());
    size_t i = 0;
    Counts<Nt4> cnts;
    for (CLocAlnSet::site_iterator site (aln_loc); bool(site); ++site, ++i) {
        site.counts(cnts);
        pair<size_t,Nt4> best_nt = cnts.sorted()[0];
        true_ref.set(i, best_nt.first > 0 ? best_nt.second : Nt4::n);
    }
    aln_loc.ref(move(true_ref));

    // Undo the paired-end read (true) alignments.
    // N.B. `from_true_alignments()` adds the paired-end reads at the end.
    vector<SRead> pe_reads;
    while(!aln_loc.reads().empty() && aln_loc.reads().back().name.back() == '2') {
        SAlnRead& r = aln_loc.reads().back();
        pe_reads.push_back(SRead(Read(move(r.seq), move(r.name)), r.sample));
        aln_loc.reads().pop_back();
    }

    if (!pe_reads.empty()) {
        // Compute the paired-end reads alignments.
        SuffixTree stree (aln_loc.ref());
        stree.build_tree();
        GappedAln aligner (pe_reads.front().seq.length(), aln_loc.ref().length(), true);
        AlignRes aln_res;
        for (SRead& r : pe_reads)
            if(add_read_to_aln(aln_loc, aln_res, move(r), &aligner, &stree))
                ++ctg_stats_.n_aln_reads;
        if (detailed_output)
            for (auto& r : aln_loc.reads())
                if (r.is_read2())
                    loc_.details_ss << "pe_aln_global"
                                << '\t' << r.name
                                << '\t' << r.cigar
                                << '\n';

        // Merge forward & paired-end reads.
        aln_loc.merge_paired_reads();
    }
}

Timers& Timers::operator+= (const Timers& other) {
    reading += other.reading;
    processing += other.processing;
    writing_fa += other.writing_fa;
    writing_vcf += other.writing_vcf;
    writing_details += other.writing_details;

    assembling += other.assembling;
    olap_aligning += other.olap_aligning;
    geno_haplotyping += other.geno_haplotyping;
    building_vcf += other.building_vcf;

    return *this;
}

//
// Arguments
// ==========
//

const string help_string = string() +
        "gstacks " + VERSION  + "\n" +
        "\n"
        "De novo mode:\n"
        "  gstacks -P stacks_dir -M popmap\n"
        "\n"
        "  -P: input directory\n"
        "\n"
        "Reference-based mode:\n"
        "  gstacks -I in_dir -M popmap -O out_dir [--paired]\n"
        "  gstacks -B bam_file [-B ...] -O out_dir [--paired]\n"
        "\n"
        "  -I: input directory containing '[SAMPLE_NAME].bam' files (or such-named\n"
        "      symbolic links).\n"
        "  -B: input BAM file(s)\n"
        "  --paired: reads are paired (RAD loci will be defined by READ1 alignments)\n"
        "\n"
        "  The input BAM file(s) must be sorted by coordinate. With -B, records\n"
        "  must be assigned to samples using BAM \"reads groups\" (gstacks uses"
        "  the ID \"identifier\" and SM \"sample name\" fields). Read groups,\n"
        "  if repeated in several files, must be consistent. With -I, read groups\n"
        "  are ignored.\n"
        "\n"
        "Shared options:\n"
        "  -M: path to a population map giving the list of samples\n"
        "  -O: output directory (default: none with -B; with -P same as the input\n"
        "      directory)\n"
        "  -t,--threads: number of threads to use (default: 1)\n"
        "  --details: write a more detailed output\n"
        "  --ignore-pe-reads: ignore paired-end reads even if present in the input\n"
        "                     (in reference-based mode, this implies --paired)\n"
        "\n"
        "Model options:\n"
        "  --model: model to use to call variants and genotypes; one of\n"
        "           marukilow (default), marukihigh, or snp\n"
        "  --var-alpha: alpha threshold for discovering SNPs (default: 0.05 for marukilow)\n"
        "  --gt-alpha: alpha threshold for calling genotypes (default: 0.05)\n"
        "\n"
        "Advanced options:\n"
        "  (De novo mode)\n"
        "  --kmer-length: kmer length for the de Bruijn graph (default: 31, max. 31)\n"
        "  --min-kmer-cov: minimum coverage to consider a kmer (default: 2)\n"
        "\n"
        "  (Reference-based mode)\n"
        "  --min-mapq: minimum PHRED-scaled mapping quality to consider a read (default: 10)\n"
        "  --max-clipped: maximum soft-clipping level, in fraction of read length (default: 0.20)\n"
        "  -m,--min-spl-reads: minimum number of reads for a sample to be considered at a locus (default: 1)\n"
        "  --max-insert-len: maximum allowed sequencing insert length (for --paired; default: 1000)\n"
        "\n"
#ifdef DEBUG
        "Debug options:\n"
        "  --dbg-no-overlaps: disable overlapping\n"
        "  --dbg-no-haps: disable phasing\n"
        "  --dbg-gfa: output a GFA file for each locus\n"
        "  --dbg-alns: output a file showing the contigs & alignments\n"
        "  --hap-graphs: output a dot graph file showing phasing information\n"
        "  --dbg-depths: write detailed depth data in the output VCF\n"
        "  --dbg-true-alns: use true alignments (for simulated data; read IDs must\n"
        "               include 'cig1=...' and 'cig2=...' fields.\n"
        "  --dbg-true-reference: align paired-end reads to the true reference\n"
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
        {"stacks-dir",   required_argument, NULL,  'P'},
        {"in-dir",       required_argument, NULL,  'I'},
        {"in-bam",       required_argument, NULL,  'B'},
        {"popmap",       required_argument, NULL,  'M'},
        {"out-dir",      required_argument, NULL,  'O'},
        {"paired",       no_argument,       NULL,  1007},
        {"threads",      required_argument, NULL,  't'},
        {"model",        required_argument, NULL,  1006},
        {"gt-alpha",     required_argument, NULL,  1005},
        {"var-alpha",    required_argument, NULL,  1008},
        {"kmer-length",  required_argument, NULL,  1001},
        {"min-kmer-cov", required_argument, NULL,  1002},
        {"ignore-pe-reads", no_argument,    NULL,  1012},
        {"details",      no_argument,       NULL,  1013},
        {"min-mapq",     required_argument, NULL,  1014},
        {"max-clipped",  required_argument, NULL,  1015},
        {"min-spl-reads", required_argument, NULL, 'm'},
        {"max-insert-len", required_argument, NULL,  1016},
        //debug options
        {"dbg-gfa",      no_argument,       NULL,  2003},
        {"dbg-alns",     no_argument,       NULL,  2004}, {"alns", no_argument, NULL, 3004},
        {"dbg-depths",   no_argument,       NULL,  2007},
        {"dbg-hap-graphs", no_argument,     NULL,  2010},
        {"dbg-true-reference", no_argument, NULL,  2012},
        {"dbg-true-alns", no_argument,      NULL,  2011}, {"true-alns", no_argument, NULL, 3011},
        {"dbg-no-overlaps", no_argument,    NULL,  2008},
        {"dbg-no-haps",  no_argument,       NULL,  2009},
        {0, 0, 0, 0}
    };

    string stacks_dir;
    string in_dir;
    string out_dir;
    double gt_alpha = 0.05;
    double var_alpha = 0.05;

    int c;
    int long_options_i;
    while (true) {

        c = getopt_long(argc, argv, "hqP:I:B:M:O:W:t:m:", long_options, &long_options_i);

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
            stacks_dir = optarg;
            break;
        case 'I':
            in_dir = optarg;
            break;
        case 'B':
            in_bams.push_back(optarg);
            break;
        case 'M':
            popmap_path = optarg;
            break;
        case 'O':
            out_dir = optarg;
            break;
        case 1007: //paired
            refbased_cfg.paired = true;
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
        case 't':
            num_threads = is_integer(optarg);
            if (num_threads < 0) {
                cerr << "Error: Illegal -t option value '" << optarg << "'.\n";
                bad_args();
            }
            break;
        case 1012: //ignore-pe-reads
            ignore_pe_reads = true;
            refbased_cfg.ign_pe_reads = true;
            break;
        case 1001://kmer-length
            km_length = atoi(optarg);
            break;
        case 1002://min-cov
            min_km_count = atoi(optarg);
            break;
        case 1013://details
            detailed_output = true;
            break;
        case 1014://min-mapq
            refbased_cfg.min_mapq = stoi(optarg);
            if (refbased_cfg.min_mapq > 255) {
                cerr << "Error: Illegal --min-mapq value '" << optarg << "'.\n";
                bad_args();
            }
            break;
        case 1015://max-clipped
            refbased_cfg.max_clipped = atof(optarg);
            if (refbased_cfg.max_clipped < 0.0 || refbased_cfg.max_clipped > 1.0) {
                cerr << "Error: Illegal --max-clipped value '" << optarg << "'.\n";
                bad_args();
            }
            break;
        case 'm'://min-spl-reads
            refbased_cfg.min_reads_per_sample = stoi(optarg);
            break;
        case 1016://max-insert
            refbased_cfg.max_insert_refsize = stoi(optarg);
            break;
    
        //
        // Debug options
        //
        case 2012://dbg-true-alns
            dbg_true_reference = true;
            break;
        case 3011:
        case 2011://dbg-true-alns
            dbg_true_alns = true;
            break;
        case 2003://dbg-gfa
            dbg_write_gfa = true;
            break;
        case 3004:
        case 2004://dbg-alns
            dbg_write_alns = true;
            break;
        case 2010://dbg-hap-graphs
            dbg_write_hapgraphs = true;
            break;
        case 2007://dbg-depths
            dbg_write_nt_depths = true;
            break;
        case 2008://dbg-no-haps
            dbg_no_overlaps = true;
            break;
        case 2009://dbg-no-haps
            dbg_no_haplotypes = true;
            break;
        case '?':
            bad_args();
            break;
        }
    }

    //
    // Check command consistency.
    //

    if (optind < argc) {
        cerr << "Error: Failed to parse command line: '" << argv[optind]
             << "' is seen as a positional argument. Expected no positional arguments.\n";
        bad_args();
    }

    size_t n_modes_given = !stacks_dir.empty() + !in_dir.empty() + !in_bams.empty();
    if (n_modes_given != 1) {
        cerr << "Error: Please specify exactly one of -P, -I or -B.\n";
        bad_args();
    }

    typedef GStacksInputT In;
    if (!stacks_dir.empty())
        input_type = popmap_path.empty()
            ? In::denovo_merger
            : In::denovo_popmap;
    else if (!in_dir.empty())
        input_type = In::refbased_popmap;
    else if (!in_bams.empty())
        input_type = In::refbased_list;
    else
        DOES_NOT_HAPPEN;

    if (input_type == In::refbased_popmap
        && popmap_path.empty()
    ) {
        cerr << "Error: Please specify a population map (-M).\n";
        bad_args();
    } else if (input_type == In::refbased_list
        && !popmap_path.empty()
    ) {
        cerr << "Error: Please specify -I/-M or -B, not both.\n";
        bad_args();
    }

    if ((input_type == In::refbased_popmap || input_type == In::refbased_list)
        && out_dir.empty()
    ) {
        cerr << "Error: Please specify an output directory (-O).\n";
        bad_args();
    }

    if (refbased_cfg.paired
        && (input_type == In::denovo_popmap || input_type == In::denovo_merger)
    ) {
        cerr << "Error: --paired is for the reference-based mode.\n";
        bad_args();
    }

    //
    // Process arguments.
    //

    for (string* dir : {&stacks_dir, &in_dir, &out_dir})
        if (!dir->empty() && dir->back() != '/')
            *dir += '/';

    switch (model_type) {
    case snp:        model.reset(new MultinomialModel(gt_alpha)); break;
    case marukihigh: model.reset(new MarukiHighModel(gt_alpha, var_alpha));  break;
    case marukilow:  model.reset(new MarukiLowModel(gt_alpha, var_alpha));   break;
    default:
        cerr << "Error: Model choice '" << to_string(model_type) << "' is not supported.\n";
        bad_args();
        break;
    }

    if (!popmap_path.empty()) {
        MetaPopInfo m;
        m.init_popmap(popmap_path);
        for (auto& s : m.samples())
            sample_names.push_back(s.name);
    }

    if (input_type == In::denovo_popmap) {
        for (const string& s : sample_names)
            in_bams.push_back(stacks_dir + s + ".matches.bam");
        if (out_dir.empty())
            out_dir = stacks_dir;
    } else if (input_type == In::denovo_merger) {
        in_bams.push_back(stacks_dir + "catalog.bam");
        if (out_dir.empty())
            out_dir = stacks_dir;
    } else if (input_type == In::refbased_popmap) {
        for (const string& s : sample_names)
            in_bams.push_back(in_dir + s + ".bam");
    } else {
        assert(input_type == In::refbased_list);
    }

    o_prefix = out_dir + prog_name;
    check_or_mk_dir(out_dir);

    if (in_bams.empty())
        DOES_NOT_HAPPEN;

} catch (std::invalid_argument&) {
    bad_args();
}
}

void report_options(ostream& os) {
    os << "Configuration for this run:\n";

    switch (input_type) {
    case GStacksInputT::denovo_popmap:
    case GStacksInputT::denovo_merger:
        os << "  Input mode: denovo\n";
        break;
    case GStacksInputT::refbased_popmap:
    case GStacksInputT::refbased_list:
        os << "  Input mode: reference-based";
        if (refbased_cfg.paired)
            os << ", paired-end\n";
        else
            os << ", single-end\n";
        break;
    default:
        DOES_NOT_HAPPEN;
        break;
    }

    if (!popmap_path.empty())
        os << "  Population map: '" << popmap_path << "'\n";

    os << "  Input files: " << in_bams.size() << ", e.g. '" << in_bams.front()<< "'\n"
       << "  Output to: '" << o_prefix << ".*'\n"
       << "  Model: " << *model << "\n";

    if (ignore_pe_reads)
        os << "  Ignoring paired-end reads.\n";
    if (km_length != 31)
        os << "  Kmer length: " << km_length << "\n";
    if (min_km_count != 2)
        os << "  Min coverage: " << min_km_count << "\n";
}
