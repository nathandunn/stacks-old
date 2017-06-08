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

struct PhasedHet {
    size_t phase_set; // N.B. The convention in VCF is to use the column of the first phased SNP for this.
    Nt2 left_allele;
    Nt2 right_allele;
};

class SnpAlleleCooccurrenceCounter {
    size_t n_snps_;
    vector<array<array<size_t,4>,4>> cooccurences_;
public:
    SnpAlleleCooccurrenceCounter(size_t n_snps)
        : n_snps_(n_snps),
          cooccurences_(n_snps_*n_snps_) // n*n matrix, athough we only use the i<j half.
        {}
    size_t& at(size_t snp_i1, Nt2 snp1_allele, size_t snp_i2, Nt2 snp2_allele)
        {assert(snp_i1 < snp_i2); return cooccurences_[snp_i1*n_snps_+snp_i2][size_t(snp1_allele)][size_t(snp2_allele)];}
    void clear()
        {for(size_t i=0; i<n_snps_; ++i) for(size_t j=i+1; j<n_snps_; ++j) cooccurences_[i*n_snps_+j] = array<array<size_t,4>,4>();}
};

void parse_command_line(int argc, char* argv[]);
void report_options(ostream& os);
bool process_one_locus(CLocReadSet&& loc);
vector<map<size_t,PhasedHet>> phase_hets(const vector<SiteCall>& calls,
                                         const CLocAlnSet& aln_loc,
                                         set<size_t>& inconsistent_samples
                                         );
void write_one_locus(const CLocAlnSet& aln_loc,
                     const vector<SiteCounts>& depths,
                     const vector<SiteCall>& calls,
                     const vector<map<size_t,PhasedHet>>& phase_data); // {col : phasedhet} maps, for all samples

//
// Argument globals.
//
bool quiet = false;
string in_dir;
int batch_id = -1;
modelt model_type = snp;
const Model* model = NULL;
set<int> locus_wl;
size_t km_length = 31;
size_t min_km_count = 2;
bool write_haplotypes = false;
bool write_gfa = false;
bool write_alns = false;
bool write_hapgraphs = false;
bool vcf_write_depths = false;

//
// Extra globals.
//
const string prog_name = "rystacks";
LogAlterator* logger = NULL;
gzFile o_gzfasta_f = NULL;
VcfWriter* o_vcf_f = NULL;
ofstream o_models_f;
ofstream o_aln_f;
ofstream o_hapgraphs_f;

int main(int argc, char** argv) {

    // Parse arguments.
    parse_command_line(argc, argv);

    // Open the log.
    string lg_path = in_dir + prog_name + ".log";
    logger = new LogAlterator(lg_path, quiet, argc, argv);
    report_options(cout);
    cout << "\n" << flush;

    // Open the BAM file and parse the header.
    BamCLocReader bam_fh (in_dir + "batch_" + to_string(batch_id) + ".catalog.bam");

    // Open the output files.
    string o_gzfasta_path = in_dir + "batch_" + to_string(batch_id) + "." + prog_name + ".fa.gz";
    o_gzfasta_f = gzopen(o_gzfasta_path.c_str(), "wb");
    check_open(o_gzfasta_f, o_gzfasta_path);

    string o_vcf_path = in_dir + "batch_" + to_string(batch_id) + "." + prog_name + ".vcf";
    VcfHeader vcf_header;
    vcf_header.add_meta(VcfMeta::predefs::info_DP);
    vcf_header.add_meta(VcfMeta::predefs::info_AF);
    vcf_header.add_meta(VcfMeta::predefs::info_AD);
    vcf_header.add_meta(VcfMeta::predefs::format_GT);
    vcf_header.add_meta(VcfMeta::predefs::format_DP);
    vcf_header.add_meta(VcfMeta::predefs::format_AD);
    vcf_header.add_meta(VcfMeta::predefs::format_GL);
    for(auto& s : bam_fh.mpopi().samples())
        vcf_header.add_sample(s.name);
    o_vcf_f = new VcfWriter(o_vcf_path, move(vcf_header));

    string o_models_path = in_dir + "batch_" + to_string(batch_id) + "." + prog_name + ".tsv";
    o_models_f.open(o_models_path);
    check_open(o_models_f, o_models_path);

    if (write_alns) {
        string o_aln_path = in_dir + "batch_" + to_string(batch_id) + "." + prog_name + ".alns";
        o_aln_f.open(o_aln_path);
        check_open(o_aln_f, o_aln_path);
        o_aln_f <<
            "# This prints observed read haplotypes:\n"
            "# loc=39\n"
            "# sample=BT_2827.13\n"
            "# cols=$(grep -E \"^$loc\\b\" batch_1.rystacks.vcf | awk '$5!=\".\"' | cut -f2 | paste -sd ',') # (SNPs.)\n"
            "# sed -n \"/^BEGIN $loc\\b/,/^END $loc\\b/ p\" batch_1.rystacks.alns | grep \"$sample\" | cut -f3 | cut -c$cols | sort | uniq -c | sort -nr\n";
    }

    if (write_hapgraphs) {
        string o_hapgraphs_path = in_dir + "batch_" + to_string(batch_id) + "." + prog_name + ".hapgraphs.dot";
        o_hapgraphs_f.open(o_hapgraphs_path);
        check_open(o_hapgraphs_f, o_hapgraphs_path);
        o_hapgraphs_f << "# dot -Tsvg -O batch_1.rystacks.hapgraphs.dot\n"
                      << "graph {\n"
                      << "edge[color=\"grey60\",arrowsize=0.8,fontsize=12,labeljust=\"l\"];\n";
    }

    // Process every locus
    cout << "Processing all loci...\n" << flush;
    CLocReadSet loc (bam_fh.mpopi());
    size_t n_loci = 0;
    size_t n_discarded = 0;
    if (locus_wl.empty()) {
        // No whitelist.
        while (bam_fh.read_one_locus(loc)) {
            ++n_loci;
            if (!process_one_locus(move(loc)))
                ++n_discarded;
        }

    } else {
        while (bam_fh.read_one_locus(loc) && !locus_wl.empty()) {
            if (locus_wl.count(loc.id())) {
                ++n_loci;
                if (!process_one_locus(move(loc)))
                    ++n_discarded;
                locus_wl.erase(loc.id());
            }
        }
    }
    cout << "Processed " << n_loci << " loci; retained " << (n_loci-n_discarded) << " of them.\n";

    gzclose(o_gzfasta_f);
    delete o_vcf_f;
    delete model;
    if (write_hapgraphs)
        o_hapgraphs_f << "}\n";

    cout << prog_name << " is done.\n";
    delete logger;
    return 0;
}

bool process_one_locus(CLocReadSet&& loc) {
    assert(!loc.reads().empty());

    //
    // Process the paired-end reads.
    //
    CLocAlnSet pe_aln_loc (loc.mpopi());
    pe_aln_loc.id(loc.id());
    do { // (Avoiding nested ifs.)
        if (loc.pe_reads().empty())
            break;

        //
        // Assemble the reads.
        //
        string pe_contig;
        vector<const DNASeq4*> seqs_to_assemble;
        for (const Read& r : loc.pe_reads())
            seqs_to_assemble.push_back(&r.seq);

        Graph graph (km_length);
        graph.rebuild(seqs_to_assemble, min_km_count);
        if (graph.empty())
            break;

        if (write_gfa)
            graph.dump_gfa(in_dir + to_string(loc.id()) + ".gfa");

        vector<const SPath*> best_path;
        if (!graph.find_best_path(best_path))
            // Not a DAG.
            break;

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

    // Build the foward-reads object.

    CLocAlnSet fw_aln_loc (loc.mpopi(), loc.id());
    fw_aln_loc.ref(DNASeq4(loc.reads().at(0).seq));
    for (SRead& r : loc.reads()) {
        if (r.seq.length() != fw_aln_loc.ref().length()) {
            cerr << "DEBUG: Error: Can't handle reads of different legnths.\n"; //xxx
            throw exception();
        }
        fw_aln_loc.add(SAlnRead(move((Read&)r), {{'M',r.seq.length()}}, r.sample));
    }

    //
    // Merge the forward & paired-end contigs.
    //
    CLocAlnSet aln_loc (loc.mpopi(), loc.id());
    if (pe_aln_loc.reads().empty()) {
        aln_loc = move(fw_aln_loc);
    } else {
        CLocAlnSet dummy (loc.mpopi(), loc.id());
        dummy.ref(DNASeq4(string(10, 'N')));
        aln_loc = CLocAlnSet::juxtapose(
                move(fw_aln_loc),
                CLocAlnSet::juxtapose(move(dummy), move(pe_aln_loc))
                );
        aln_loc.merge_paired_reads();
    }

    if (write_alns)
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
    if (write_haplotypes) {
        set<size_t> inconsistent_samples;
        phase_data = phase_hets(calls, aln_loc, inconsistent_samples);
        /* TODO {
        for (size_t sample : inconsistent_samples) {
            // Observed haplotypes are inconsistent given the sample's diploidy.
            for (size_t i=0; i<aln_loc.ref().length(); ++i) {
                depths[i].samples[sample] = Counts<Nt2>();
                calls[i].discard_sample(sample);
            }
        }
        // TODO } */
    }

    write_one_locus(aln_loc, depths, calls, phase_data);

    return true;
}

vector<map<size_t,PhasedHet>> phase_hets(const vector<SiteCall>& calls,
                                         const CLocAlnSet& aln_loc,
                                         set<size_t>& inconsistent_samples
                                         ){
    vector<map<size_t,PhasedHet>> phased_samples (aln_loc.mpopi().samples().size());

    vector<size_t> snp_cols;
    for (size_t i=0; i<aln_loc.ref().length(); ++i)
        if (calls[i].alleles().size() > 1)
            snp_cols.push_back(i);

    if (snp_cols.empty())
        return phased_samples;

    if (write_hapgraphs) {
        o_hapgraphs_f << "subgraph cluster_loc" << aln_loc.id() << " {\n"
                      << "\tlabel=\"locus " << aln_loc.id() << "\";\n"
                      << "\t# snp columns: ";
        join(snp_cols, ',', o_hapgraphs_f);
        o_hapgraphs_f << "\n";
    }

    SnpAlleleCooccurrenceCounter cooccurences (snp_cols.size());
    for (size_t sample=0; sample<aln_loc.mpopi().samples().size(); ++sample) {
        if (aln_loc.sample_reads(sample).empty())
            continue;

        vector<size_t> het_snps;
        for(size_t snp_i=0; snp_i<snp_cols.size(); ++snp_i)
            if (calls[snp_cols[snp_i]].sample_calls()[sample].call() == snp_type_het)
                het_snps.push_back(snp_i);

        if (het_snps.size() == 0) {
            continue;
        } else if (het_snps.size() == 1) {
            // Sample has trivial 1nt-long haplotypes.
            size_t col = snp_cols[het_snps[0]];
            const SampleCall& c = calls[col].sample_calls()[sample];
            phased_samples[sample].insert({col, {col, c.nt0(), c.nt1()}});
            continue;
        }

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
                const SampleCall& c = calls[col].sample_calls()[sample];
                read_hap[het_i] = (nt == c.nt0() || nt == c.nt1()) ? nt : Nt4::n;
            }

            // Record the pairwise cooccurrences.
            for (size_t i=0; i<het_snps.size(); ++i) {
                for (size_t j=i+1; j<het_snps.size(); ++j) {
                    Nt4 nti = read_hap[i];
                    Nt4 ntj = read_hap[j];
                    if (nti == Nt4::n || ntj == Nt4::n)
                        continue;
                    ++cooccurences.at(het_snps[i], Nt2(nti), het_snps[j], Nt2(ntj));
                }
            }
        }

        // Initialize the dot graph, if required.
        auto nodeid = [&aln_loc,&sample](size_t col, Nt2 allele)
                {return string("l")+to_string(aln_loc.id())+"s"+to_string(sample)+"c"+to_string(col)+char(allele);};
        if (write_hapgraphs) {
            o_hapgraphs_f << "\tsubgraph cluster_sample" << sample << " {\n"
                          << "\t\tlabel=\"sample" << sample << "\";\n"
                          << "\t\tstyle=dashed;\n"
                          << "\t\t# heterozygous columns: ";
            vector<size_t> het_cols;
            for (size_t snp_i : het_snps)
                het_cols.push_back(snp_cols[snp_i]);
            join(het_cols, ',', o_hapgraphs_f);
            o_hapgraphs_f << "\n";

            for (size_t snp_i : het_snps) {
                size_t col = snp_cols[snp_i];
                const SampleCall& c = calls[col].sample_calls()[sample];
                for (Nt2 allele : {c.nt0(), c.nt1()})
                    o_hapgraphs_f << "\t\t" << nodeid(col, allele)
                                  << " [label=<"
                                  << "<sup><font point-size=\"10\">" << col << "</font></sup>" << allele
                                  << ">];\n";
            }
        }

        // Iterate over pairwise cooccurrences.
        for (size_t het_i=0; het_i<het_snps.size(); ++het_i) {
            size_t snp_i = het_snps[het_i];
            size_t coli = snp_cols[snp_i];
            for (size_t het_j=het_i+1; het_j<het_snps.size(); ++het_j) {
                size_t snp_j = het_snps[het_j];
                size_t colj = snp_cols[snp_j];
                for (auto& nti : calls[coli].alleles()) {
                    for (auto& ntj : calls[colj].alleles()) {
                        size_t n = cooccurences.at(snp_i, nti.first, snp_j, ntj.first);
                        if (n == 0)
                            continue;
                        if (write_hapgraphs) {
                            o_hapgraphs_f << "\t\t" << nodeid(coli,nti.first) << " -- " << nodeid(colj,ntj.first) << " [";
                            if (n==1)
                                o_hapgraphs_f << "style=dotted";
                            else
                                o_hapgraphs_f << "label=\"" << n << "\",penwidth=" << n;
                            o_hapgraphs_f << "];\n";
                        }
                        // ...
                    }
                }
            }
        }

        if (write_hapgraphs)
            o_hapgraphs_f << "\t}\n";
    }
    if (write_hapgraphs)
        o_hapgraphs_f << "}\n";

    return phased_samples;
}

void write_one_locus(
        const CLocAlnSet& aln_loc,
        const vector<SiteCounts>& depths,
        const vector<SiteCall>& calls,
        const vector<map<size_t,PhasedHet>>& phase_data
        ){

    size_t loc_id = aln_loc.id();
    const DNASeq4& ref = aln_loc.ref();
    const MetaPopInfo& mpopi = aln_loc.mpopi();

    //
    // Vcf output.
    //
    assert(depths.size() == ref.length());
    assert(calls.size() == ref.length());
    vector<size_t> sample_sites_w_data (mpopi.samples().size(), 0);
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

        // Create the VCF record.
        VcfRecord rec;
        rec.type_m() = Vcf::RType::expl;
        rec.chrom_m() = to_string(loc_id);
        rec.pos_m() = i+1;

        // Alleles.
        for (Nt2 nt : vcf_alleles)
            rec.alleles_m().push_back(string(1, char(nt)));

        if(rec.n_alleles() == 1) {
            // Fixed site.

            // Info/DP.
            rec.info_m().push_back(string("DP=") + to_string(sitedepths.tot.sum()));
            // Info/AD.
            Nt2 ref_nt = sitecall.alleles().begin()->first;
            rec.info_m().push_back(string("AD=") + to_string(sitedepths.tot[ref_nt]));
            if (vcf_write_depths) {
                // Info/cnts.
                stringstream cnts;
                join(sitedepths.tot.arr(), ',', cnts);
                rec.info_m().push_back(string("cnts=") + cnts.str());
            }
            // Format.
            rec.format_m().push_back("DP");
            // Genotypes.
            for (size_t sample=0; sample<mpopi.samples().size(); ++sample) {
                size_t dp = sitedepths.samples[sample].sum();
                if (dp == 0) {
                    rec.samples_m().push_back(".");
                    continue;
                }
                ++sample_sites_w_data[sample];
                rec.samples_m().push_back(to_string(dp));
            }

        } else {
            // Polymorphic site.

            // Info/DP.
            rec.info_m().push_back(string("DP=") + to_string(sitedepths.tot.sum()));
            // Info/AD.
            vector<size_t> ad;
            for (Nt2 nt : vcf_alleles)
                ad.push_back(sitedepths.tot[nt]);
            stringstream ss;
            join(ad, ',', ss);
            rec.info_m().push_back(string("AD=") + ss.str());
            // Info/AF.
            vector<double> alt_freqs;
            for (auto nt=++vcf_alleles.begin(); nt!=vcf_alleles.end(); ++nt) // rem. always >1 alleles.
                alt_freqs.push_back(sitecall.alleles().at(*nt));
            rec.info_m().push_back(VcfRecord::util::fmt_info_af(alt_freqs));
            if (vcf_write_depths) {
                // Info/cnts.
                stringstream cnts;
                join(sitedepths.tot.arr(), ',', cnts);
                rec.info_m().push_back(string("cnts=") + cnts.str());
            }

            // Format.
            rec.format_m().push_back("GT");
            if (write_haplotypes)
                rec.format_m().push_back("PS"); // Phase set.
            rec.format_m().push_back("DP");
            rec.format_m().push_back("AD");
            rec.format_m().push_back("GL");

            // Genotypes.
            rec.samples_m().reserve(mpopi.samples().size());
            for (size_t sample=0; sample<mpopi.samples().size(); ++sample) {
                const Counts<Nt2>& sdepths = sitedepths.samples[sample];
                const SampleCall& scall = sitecall.sample_calls()[sample];
                if (sdepths.sum() == 0) {
                    // No data for this sample.
                    rec.samples_m().push_back(".");
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
                    if (write_haplotypes)
                        genotype << ":.";
                    break;
                case snp_type_het:
                    if (write_haplotypes) {
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
                    if (write_haplotypes)
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
                if (vcf_write_depths) {
                    genotype << ":";
                    join(sdepths.arr(), ',', genotype);
                }
                // Push it.
                rec.samples_m().push_back(genotype.str());
            }
        }

        // Write the record.
        o_vcf_f->write_record(rec);
    }

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

    // Write the fasta record.
    gzputs(o_gzfasta_f, ">");
    gzputs(o_gzfasta_f, to_string(loc_id).c_str());
    gzputs(o_gzfasta_f, " NS=");
    gzputs(o_gzfasta_f, to_string(n_remaining_samples).c_str());
    if (n_remaining_samples != samples_w_reads.size()) {
        assert(n_remaining_samples < samples_w_reads.size());
        gzputs(o_gzfasta_f, " n_discarded_samples=");
        gzputs(o_gzfasta_f, to_string(samples_w_reads.size() - n_remaining_samples).c_str());
    }
    gzputs(o_gzfasta_f, "\n");
    gzputs(o_gzfasta_f, ref.str().c_str());
    gzputs(o_gzfasta_f, "\n");

    //
    // Models/tsv output.
    // LOCID \t LINETYPE \t SAMPLEID \t CONTENTS
    //

    // Consensus.
    o_models_f << loc_id << "\tconsensus\t\t" << ref << "\n";

    // Model.
    o_models_f << loc_id << "\tmodel\t\t";
    for (auto& c : calls)
        o_models_f << c.alleles().size();
    o_models_f << "\n";

    // Depth.
    // One two-digit hex number per position (max 0xFF).
    o_models_f << loc_id << "\tdepth\t\t" << std::hex;
    for (size_t i=0; i<ref.length(); ++i) {
        size_t dp = depths[i].tot.sum();
        if (dp <= 0xF)
            o_models_f << "0" << dp;
        else if (dp <= 0xFF)
            o_models_f << dp;
        else
            o_models_f << 0xFF;
    }
    o_models_f << std::dec << "\n";

    // For each sample.
    for (size_t s=0; s<mpopi.samples().size(); ++s) {
        int sample_id = mpopi.samples()[s].id;

        // Model.
        o_models_f << loc_id << "\ts_model\t" << sample_id << "\t";
        for (size_t i=0; i<ref.length(); ++i) {
            const SiteCall& c = calls[i];
            if (c.alleles().size() == 0) {
                o_models_f << "U";
            } else if (c.alleles().size() == 1) {
                o_models_f << "O";
            } else {
                switch (c.sample_calls()[s].call()) {
                case snp_type_hom: o_models_f << "O"; break;
                case snp_type_het: o_models_f << "E"; break;
                case snp_type_unk: o_models_f << "U"; break;
                }
            }
        }
        o_models_f << "\n";

        // Depths.
        // Four two-digit hex numbers (A,C,T,G) per position.
        o_models_f << loc_id << "\ts_depths\t" << sample_id << "\t" << std::hex;
        for (size_t i=0; i<ref.length(); ++i) {
            // For each site...
            const SiteCounts& sitedepths = depths[i];
            for (Nt2 nt : Nt2::all) {
                // For each of A, C, G and T...
                size_t dp = sitedepths.samples[s][nt];
                if (dp <= 0xF)
                    o_models_f << "0" << dp;
                else if (dp <= 0xFF)
                    o_models_f << dp;
                else
                    o_models_f << 0xFF;
            }
        }
        o_models_f << std::dec << "\n";
    }
}

const string help_string = string() +
        prog_name + " " + VERSION  + "\n" +
        prog_name + " -P in_dir\n"
        "\n"
        "  -P: input directory (must contain a batch_X.catalog.bam file)\n"
        "  -b: batch ID (default: guess)\n"
        "  -W,--whitelist: a whitelist of locus IDs\n"
        "  --model: model to use to call variants and genotypes;\n"
        "           one of snp (default), marukihigh, or marukilow\n"
        "  --gt-alpha: alpha threshold for calling genotypes (default: 0.05)\n"
        "  --var-alpha: alpha threshold for discovering variants (default: 0.05)\n"
        "\n"
        "Debug options:\n"
        "  --kmer-length: kmer length (default: 31)\n"
        "  --min-cov: minimum coverage to consider a kmer (default: 2)\n"
        "  --gfa: output a GFA file for each locus\n"
        "  --alns: output a file showing the contigs & alignments\n"
        "  --hap-graphs: output a dot graph file showing phasing information\n"
        "  --depths: write detailed depth data in the output VCF\n"
        "  --haps: output a phased VCF /!\\ conflicts with paired-end input\n"
        "\n"
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
        {"model",        required_argument, NULL,  1006},
        {"gt-alpha",     required_argument, NULL,  1005},
        {"var-alpha",    required_argument, NULL,  1008},
        {"whitelist",    required_argument, NULL,  'W'},
        {"kmer-length",  required_argument, NULL,  1001},
        {"min-cov",      required_argument, NULL,  1002},
        {"haps",         no_argument,       NULL,  1009},
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

        c = getopt_long(argc, argv, "hqP:b:W:", long_options, &long_options_i);

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
        case 1001://kmer-length
            km_length = atoi(optarg);
            break;
        case 1002://min-cov
            min_km_count = atoi(optarg);
            break;
        case 1009://haps
            write_haplotypes = true;
            break;
        case 1003://gfa
            write_gfa = true;
            break;
        case 1004://aln
            write_alns = true;
            break;
        case 1010://hap-graphs
            write_hapgraphs = true;
            break;
        case 1007://depths
            vcf_write_depths = true;
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
    case snp:        model = new MultinomialModel(gt_alpha); break;
    case marukihigh: model = new MarukiHighModel(gt_alpha, var_alpha);  break;
    case marukilow:  model = new MarukiLowModel(gt_alpha, var_alpha);   break;
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
        if (!wl_fh) {
            cerr << "Error: Failed to open '" << wl_path << "' for reading.\n";
            throw exception();
        }
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
    if (km_length != 31)
        os << "  Kmer length: " << km_length << "\n";
    if (min_km_count != 2)
        os << "  Min coverage: " << min_km_count << "\n";
}
