#include <getopt.h>

#include "zlib.h"

#include "constants.h"
#include "utils.h"
#include "log_utils.h"
#include "catalog_utils.h"
#include "locus.h"
#include "BamCLocReader.h"
#include "debruijn.h"
#include "GappedAln.h"
#include "aln_utils.h"
#include "Alignment.h"
#include "models.h"

using namespace std;

class SiteCall;
struct SampleCall;
void parse_command_line(int argc, char* argv[]);
void report_options(ostream& os);
bool process_one_locus(CLocReadSet&& loc);
void write_one_locus(const CLocAlnSet& aln_loc, const vector<SiteCall>& calls);

class SiteCall {
    size_t tot_depth_;
    map<Nt4, size_t> alleles_;
    vector<SampleCall> sample_calls_;

public:
    SiteCall(const CLocAlnSet::site_iterator& site);

    size_t tot_depth() const {return tot_depth_;}
    const map<Nt4, size_t>& alleles() const {return alleles_;}
    const vector<SampleCall>& sample_calls() const {return sample_calls_;}

    //Nt4 most_frequent_nt() const;
};

struct SampleCall {
    array<size_t, 4> depths;
    snp_type call;
    array<Nt4, 2> nts; // hom {nt, Nt4::n} | het {min_nt, max_nt} | unk {Nt4::n, Nt4::n}

    SampleCall() : depths{0, 0, 0, 0}, call(snp_type_unk), nts{0, 0} {}

    size_t tot_depth() const {return depths[0]+depths[1]+depths[2]+depths[3];}
};

/*
Nt4 SiteCall::most_frequent_nt() const {
    assert(!alleles.empty());

    auto iter = alleles.begin();
    Nt4 nt = iter->first;
    size_t count = iter->second;
    ++iter;

    while (iter != alleles.end()) {
        if (iter->second > count) {
            nt = iter->first;
            count = iter->second;
        }
        ++iter;
    }

    return nt;
}
*/

//
// Argument globs.
//
bool quiet = false;
string in_dir;
int batch_id = -1;
double gt_alpha = 0.05;
set<int> locus_wl;
size_t km_length = 31;
size_t min_km_count = 2;
bool gfa_out = false;
bool aln_out = false;

//
// Extra globs.
//
const string prog_name = "rystacks";
LogAlterator* lg = NULL;
gzFile o_gzfasta_f = NULL;
VcfWriter* o_vcf_f = NULL;
ofstream o_models_f;
ofstream o_aln_f;

int main(int argc, char** argv) {

    // Parse arguments.
    parse_command_line(argc, argv);

    // Open the log.
    string lg_path = in_dir + prog_name + ".log";
    if(!quiet)
        cout << "Logging to '" << lg_path << "'." << endl;
    lg = new LogAlterator(lg_path, quiet);
    init_log(lg->l, argc, argv);
    report_options(cout);
    cout << "\n" << flush;

    // Initialize the model globs.
    set_model_thresholds(gt_alpha);

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
    vcf_header.add_meta(VcfMeta::predefs::format_GT);
    vcf_header.add_meta(VcfMeta::predefs::format_DP);
    vcf_header.add_meta(VcfMeta::predefs::format_AD);
    for(auto& s : bam_fh.mpopi().samples())
        vcf_header.add_sample(s.name);
    o_vcf_f = new VcfWriter(o_vcf_path, move(vcf_header));

    string o_models_path = in_dir + "batch_" + to_string(batch_id) + "." + prog_name + ".tsv";
    o_models_f.open(o_models_path);
    check_open(o_models_f, o_models_path);

    if (aln_out) {
        string o_aln_path = in_dir + "batch_" + to_string(batch_id) + "." + prog_name + ".aln";
        o_aln_f.open(o_aln_path);
        check_open(o_aln_f, o_aln_path);
        o_aln_f << "# id=123; sed -n \"/^BEGIN $id$/,/^END $id$/p\" batch_" << batch_id << ".rystacks.aln | less -x25\n";
    }

    // Process every locus
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

    cout << prog_name << " is done.\n";
    delete lg;
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

        if (gfa_out)
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
    for (SRead& r : loc.reads())
        fw_aln_loc.add(SAlnRead(move((Read&)r), {{'M',r.seq.length()}}, r.sample));

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
    }

    if (aln_out)
        o_aln_f << "BEGIN " << aln_loc.id() << "\n"
                << aln_loc
                << "\nEND " << aln_loc.id() << "\n";

    //
    // Call SNPs.
    //
    vector<SiteCall> calls;
    CLocAlnSet::site_iterator site (aln_loc);
    while(site) {
        calls.push_back(SiteCall(site));
        ++site;
    }

    /*
    // Update the consensus sequence, if necessary.
    // TODO Do the cases below happen? Add print statements...?
    DNASeq4 new_ref = aln_loc.ref();
    assert(calls.size() == new_ref.length());
    for (size_t i=0; i<calls.size(); ++i) {
        const SiteCall& c = calls[i];
        if (c.alleles.size() > 0) {
            Nt4 most_frequent = c.most_frequent_nt();
            if (new_ref[i] != most_frequent)
                new_ref.set(i, most_frequent);
        }
        // else if (new_ref[i] != Nt4::n) {}
    }
    aln_loc.ref(move(new_ref));
    */

    write_one_locus(aln_loc, calls);

    return true;
}

SiteCall::SiteCall(const CLocAlnSet::site_iterator& site)
        : tot_depth_(0), alleles_(), sample_calls_()
        {

    // N.B. For now we use the old binomial model.

    //
    // Look at this site in each sample; make genotype calls.
    //
    Nt4Counts counts;
    for (size_t s=0; s<site.mpopi().samples().size(); ++s) {
        SampleCall s_call;

        site.counts(counts, s);
        if (counts.at_rank(0) == 0) {
            s_call.depths = {0, 0, 0, 0};
            s_call.call = snp_type_unk;
            s_call.nts = {Nt4::n, Nt4::n};
            sample_calls_.push_back(s_call);
            continue;
        }

        s_call.depths = {counts.at(Nt4::a), counts.at(Nt4::c), counts.at(Nt4::g), counts.at(Nt4::t)};
        tot_depth_ += s_call.tot_depth();

        s_call.call = call_snp(lr_multinomial_model(counts.at_rank(0), counts.at_rank(1), counts.at_rank(2), counts.at_rank(3)));
        switch (s_call.call) {
        case snp_type_hom:
            s_call.nts = {counts.nt_of_rank(0), Nt4::n};
            break;
        case snp_type_het:
            s_call.nts = {counts.nt_of_rank(0), counts.nt_of_rank(1)};
            sort(s_call.nts.begin(), s_call.nts.end());
            break;
        default:
            s_call.nts = {Nt4::n, Nt4::n};
            break;
        }
        sample_calls_.push_back(s_call);

    }

    //
    // Iterate over the SampleCalls & record the genotypes that were found.
    //
    counts.reset();
    for (const SampleCall& sc : sample_calls_) {
        switch (sc.call) {
        case snp_type_hom : counts.increment(sc.nts[0]); counts.increment(sc.nts[0]); break;
        case snp_type_het : counts.increment(sc.nts[0]); counts.increment(sc.nts[1]); break;
        default: break;
        }
    }
    counts.sort();
    if (counts.at_rank(0) > 0) {
        // At least one allele was observed.
        alleles_.insert({counts.nt_of_rank(0), counts.at_rank(0)});

        if (counts.at_rank(1) > 0) {
            // SNP with at least two alleles.
            alleles_.insert({counts.nt_of_rank(1), counts.at_rank(1)});

            if (counts.at_rank(2) > 0) {
                alleles_.insert({counts.nt_of_rank(2), counts.at_rank(2)});

                if (counts.at_rank(3) > 0) {
                    // Quaternary SNP.
                    alleles_.insert({counts.nt_of_rank(3), counts.at_rank(3)});
                }
            }
        }
    }
}

void write_one_locus(const CLocAlnSet& aln_loc, const vector<SiteCall>& calls) {

    size_t loc_id = aln_loc.id();
    const DNASeq4& ref = aln_loc.ref();
    const MetaPopInfo& mpopi = aln_loc.mpopi();

    //
    // Fasta output.
    //

    // Determine the number of samples that have reads for this locus.
    set<size_t> loc_samples;
    for (const SAlnRead& r : aln_loc.reads())
        loc_samples.insert(r.sample);

    // Write the fasta record.
    gzputs(o_gzfasta_f, ">");
    gzputs(o_gzfasta_f, (to_string(loc_id) + " NS=" + to_string(loc_samples.size())).c_str());
    gzputs(o_gzfasta_f, "\n");
    gzputs(o_gzfasta_f, ref.str().c_str());
    gzputs(o_gzfasta_f, "\n");

    //
    // Vcf output.
    //
    assert(calls.size() == ref.length());
    for (size_t i=0; i<ref.length(); ++i) {
        const SiteCall& sitecall = calls[i];
        if (sitecall.alleles().empty())
            // No data at this site. xxx Reconsider?
            continue;

        // Determine which alleles exist, and their order.
        // (n.b. As of Apr4,2017 the ref allele might not be the most frequent one.)
        vector<Nt4> vcf_alleles;
        map<Nt4, size_t> vcf_allele_indexes;
        {
            vcf_alleles.push_back(ref[i]);
            vcf_allele_indexes.insert({ref[i], 0});

            // Sort the alleles by frequency.
            vector<pair<size_t, Nt4>> sorted_alleles;
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
        rec.type = Vcf::RType::expl;
        rec.chrom = to_string(loc_id);
        rec.pos = i+1;

        // Alleles.
        for (Nt4 nt : vcf_alleles)
            rec.alleles.push_back(string(1, char(nt)));

        // INFO/DP.
        rec.info.push_back({"DP", to_string(sitecall.tot_depth())});

        // INFO/AF.
        if (vcf_alleles.size() > 1) {
            size_t tot_count = 0;
            for (auto& a : sitecall.alleles())
                tot_count += a.second;

            vector<double> alt_freqs;
            for (auto nt=++vcf_alleles.begin(); nt!=vcf_alleles.end(); ++nt) // rem. always >1 alleles.
                alt_freqs.push_back((double)sitecall.alleles().at(*nt) / tot_count);

            rec.info.push_back(VcfRecord::util::fmt_info_af(alt_freqs));
        }

        // Genotypes.
        rec.format.push_back("GT");
        rec.format.push_back("DP");
        rec.format.push_back("AD");
        rec.samples.reserve(mpopi.samples().size());
        for (size_t s=0; s<mpopi.samples().size(); ++s) {
            const SampleCall& s_call = sitecall.sample_calls()[s];

            size_t dp = s_call.depths[0] + s_call.depths[1] + s_call.depths[2] + s_call.depths[3];
            if (dp == 0) {
                // No data for this sample.
                rec.samples.push_back(".");
                continue;
            }

            stringstream genotype;

            // GT field.
            vector<size_t> gt;
            switch (s_call.call) {
            case snp_type_hom:
                gt.push_back(vcf_allele_indexes.at(s_call.nts[0]));
                genotype << gt[0] << '/' << gt[0];
                break;
            case snp_type_het:
                gt.push_back(vcf_allele_indexes.at(s_call.nts[0]));
                gt.push_back(vcf_allele_indexes.at(s_call.nts[1]));
                sort(gt.begin(), gt.end()); // (Prevents '1/0'.)
                genotype << gt[0] << '/' << gt[1];
                break;
            default:
                genotype << '.';
                break;
            }

            // DP field.
            genotype << ':' << dp;

            // AD field.
            vector<size_t> ad;
            ad.reserve(vcf_alleles.size());
            for (Nt4 nt : vcf_alleles)
                ad.push_back(s_call.depths[size_t(Nt2(nt))]);
            genotype << ':';
            join(ad, ',', genotype);

            rec.samples.push_back(genotype.str());
        }

        // Write the record.
        o_vcf_f->write_record(rec);
    }

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
    for (auto& c : calls) {
        size_t dp = c.tot_depth();
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
        for (auto& c : calls) {
            switch (c.sample_calls()[s].call) {
            case snp_type_hom: o_models_f << "O"; break;
            case snp_type_het: o_models_f << "E"; break;
            case snp_type_unk: o_models_f << "U"; break;
            }
        }
        o_models_f << "\n";

        // Depths.
        // Four two-digit hex numbers (A,C,T,G) per position.
        o_models_f << loc_id << "\ts_depths\t" << sample_id << "\t" << std::hex;
        for (auto& c : calls) {
            // For each site/position.
            for (size_t nt=0; nt<4; ++nt) {
                size_t dp = c.sample_calls()[s].depths[nt];
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
        "  --gt-alpha: alpha threshold for calling genotypes\n"
        "\n"
        "Alignment options:\n"
        "  --kmer-length: kmer length (default: 31)\n"
        "  --min-cov: minimum coverage to consider a kmer (default: 2)\n"
        "  --gfa: output a GFA file for each locus\n"
        "  --aln: output a file showing the contigs & alignments\n"
        "\n"
        ;

void parse_command_line(int argc, char* argv[]) {

    auto bad_args = [](){
        cerr << help_string;
        exit(13);
    };

    static const option long_options[] = {
        {"version",      no_argument,       NULL,  1000},
        {"help",         no_argument,       NULL,  'h'},
        {"quiet",        no_argument,       NULL,  'q'},
        {"in-dir",       required_argument, NULL,  'P'},
        {"batch-id",     required_argument, NULL,  'b'},
        {"gt-alpha",     required_argument, NULL,  1005},
        {"whitelist",    required_argument, NULL,  'W'},
        {"kmer-length",  required_argument, NULL,  1001},
        {"min-cov",      required_argument, NULL,  1002},
        {"gfa",          no_argument,       NULL,  1003},
        {"aln",          no_argument,       NULL,  1004},
        {0, 0, 0, 0}
    };

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
        case 1005: //gt-alpha
            gt_alpha = atof(optarg);
            if (gt_alpha != 0.1 && gt_alpha != 0.05 && gt_alpha != 0.01 && gt_alpha != 0.001) {
                cerr << "Error: Illegal --gt-alpha value; pick one of {0.1, 0.05, 0.01, 0.001}.\n";
                bad_args();
            }
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
        case 1003://gfa
            gfa_out = true;
            break;
        case 1004://aln
            aln_out = true;
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
}

void report_options(ostream& os) {
    os << "Configuration for this run:\n"
       << "  Input directory: '" << in_dir << "'\n"
       << "  Batch ID: " << batch_id << "\n"
       ;
    if (!locus_wl.empty())
        os << "  Whitelist of " << locus_wl.size() << " loci.\n";
    if (km_length != 31)
        os << "  Kmer length: " << km_length << "\n";
    if (min_km_count != 2)
        os << "  Min coverage: " << min_km_count << "\n";
}
