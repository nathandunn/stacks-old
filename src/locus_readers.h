#ifndef LOCUS_READERS_H
#define LOCUS_READERS_H

#include "constants.h"
#include "MetaPopInfo.h"
#include "BamI.h"
#include "locus.h"
#include "Vcf.h"

void build_mpopi(MetaPopInfo& mpopi, map<string, size_t>& rg_to_sample, const Bam& bam_f);

class BamCLocReader {
    Bam* bam_f_;
    bool eof_;
    MetaPopInfo mpopi_;
    map<string, size_t> rg_to_sample_;
    vector<const char*> loci_aln_positions_;

    int32_t loc_i_; // Index of the next locus (chromosome). Incremented by read_one_locus().

public:
    BamCLocReader(Bam** bam_f);
    ~BamCLocReader() {if (bam_f_) delete bam_f_;}

    const Bam* bam_f() const {return bam_f_;}
    size_t n_loci() const {return bam_f_->h().n_ref_chroms();}
    int target2id(size_t target_i) const {return atoi(bam_f_->h().chrom_str(target_i));}
    const MetaPopInfo& mpopi() const {return mpopi_;}

    // Reads one locus. Returns false on EOF.
    bool read_one_locus(CLocReadSet& readset);
};

class BamCLocBuilder {
public:
    BamCLocBuilder(Bam** bam_f);
    ~BamCLocBuilder() {if (bam_f_) delete bam_f_;}

    // Reads one locus. Returns false on EOF.
    bool build_one_locus(CLocAlnSet& readset);    
    
    const MetaPopInfo& mpopi() const {return mpopi_;}
    const Bam* bam_f() const {return bam_f_;}

private:
    static const size_t max_insert_refsize_ = 2000;
    static const size_t min_reads_per_sample_ = 3;

    Bam* bam_f_;
    MetaPopInfo mpopi_;
    map<string, size_t> rg_to_sample_;
    size_t n_loci_built_;

    bool eof_;
    pair<PhyLoc,SAlnRead> next_record_;
    bool read_and_parse_next_record(); // Calls `bam_f_->next_record()`, then sets `next_record_`. Returns false on EOF.

    // Buffers to store the reads of the sliding window.
    bool fill_window();
    map<PhyLoc,vector<SAlnRead>> fw_reads_by_5prime_pos;
    //std::multimap<PhyLoc,SAlnRead> pe_reads_by_5prime_pos;
    //map<const char*,std::multimap<PhyLoc,SAlnRead>::iterator, LessCStrs> pe_reads_by_name;
};

class VcfCLocReader {
    VcfParser vcf_f_;
    VcfRecord next_rec_;
    bool eof_;
 public:
    VcfCLocReader(): vcf_f_(), next_rec_(), eof_(false) {};
    VcfCLocReader(const string& vcf_path);

    int  open(string &vcf_path);
    const VcfHeader& header() const {return vcf_f_.header();}
    void set_sample_ids(MetaPopInfo& mpopi) const;

    // Reads one locus. Returns false on EOF.
    bool read_one_locus(vector<VcfRecord>& records);
};

//
// Inline definitions
// ==================
//

inline
void
build_mpopi(
        MetaPopInfo& mpopi,
        map<string, size_t>& rg_to_sample,
        const Bam& bam_f
) {
    assert(mpopi.samples().empty());
    assert(rg_to_sample.empty());

    vector<string> rg_ids;
    vector<string> samples;
    vector<int> sample_ids;

    // Parse the @RG header lines.
    const char* p = bam_f.h().text();
    size_t line = 1;
    try {
        while (true) {
            if (strncmp(p, "@RG\t", 4) == 0) {
                // Sample line.
                rg_ids.push_back(string());
                samples.push_back(string());
                sample_ids.push_back(-1);

                p += 4;
                while (*p && *p != '\n') {
                    const char* q = p;
                    while (*p && *p != '\n' && *p != '\t')
                        ++p;

                    if (strncmp(q, "ID:", 3) == 0) {
                        rg_ids.back() = string(q+3, p);
                    } else if (strncmp(q, "SM:", 3) == 0) {
                        samples.back() = string(q+3, p);
                    } else if (strncmp(q, "id:", 3) == 0) {
                        char* end;
                        sample_ids.back() = strtol(q+3, &end, 10);
                        if (end != p)
                            throw exception();
                    }
                    if (*p == '\t')
                        ++p;
                }

                if (rg_ids.back().empty() || samples.back().empty() || sample_ids.back() == -1)
                    throw exception();
            }

            p = strchr(p, '\n');
            if (p == NULL)
                break;
            ++p;
            ++line;
        }
    } catch (exception&) {
        cerr << "Error: Malformed BAM header, at line " << line << ".\n";
        throw;
    }

    // Initialize the MetaPopInfo.
    mpopi.init_names(samples);

    // Set the sample IDs.
    assert(sample_ids.size() == samples.size());
    for (size_t s=0; s<samples.size(); ++s)
        mpopi.set_sample_id(mpopi.get_sample_index(samples[s]),sample_ids[s]);

    // Initialize `rg_to_sample_`.
    assert(rg_ids.size() == samples.size());
    for (size_t s=0; s<mpopi.samples().size(); ++s)
        rg_to_sample[rg_ids[s]] = mpopi.get_sample_index(samples[s]);
}

inline
BamCLocReader::BamCLocReader(Bam** bam_f)
        : bam_f_(*bam_f),
          eof_(false),
          mpopi_(),
          rg_to_sample_(),
          loc_i_(-1)
        {

    *bam_f = NULL;

    //
    // Create the MetaPopInfo object from the header.
    //
    build_mpopi(mpopi_, rg_to_sample_, *bam_f_);

    //
    // If the first locus has alignment information: this is a ref-based analysis,
    // record all the alignment positions. Otherwise, assume it's a de novo analysis.
    //
    {
        const char* p = bam_f_->h().text();
        size_t line = 1;
        while (true) {
            if (strncmp(p, "@SQ\t", 4) == 0) {
                loci_aln_positions_.push_back(NULL);
                while (*p && *p != '\n') {
                    if (strncmp(p, "pos:", 4) == 0) {
                        loci_aln_positions_.back() = p+4;
                        break; // Skip the rest of the line.
                    }
                    while (*p && *p != '\n' && *p != '\t')
                        ++p;
                    if (*p == '\t')
                        ++p;
                }
                if (loci_aln_positions_.back() == NULL) {
                    if (loci_aln_positions_.size() == 1) {
                        // de novo.
                        loci_aln_positions_ = vector<const char*>();
                        break;
                    } else {
                        cerr << "Error: In BAM header, at line " << line
                             << ": alignment information is missing.\n";
                        throw exception();
                    }
                }
            }

            p = strchr(p, '\n');
            if (p == NULL)
                break;
            ++p;
            ++line;
        }
        assert(loci_aln_positions_.empty() || loci_aln_positions_.size() == bam_f_->h().n_ref_chroms());
    }

    //
    // Read the very first record.
    //
    if (!bam_f_->next_record()) {
        cerr << "Error: No records in BAM file '" << bam_f_->path << "'.\n";
        throw exception();
    } else if (bam_f_->r().is_unmapped()) {
        cerr << "Error: BAM file '" << bam_f_->path << "' unexpectedly contains unmapped records.\n";
        throw exception();
    }
}

inline
bool BamCLocReader::read_one_locus(CLocReadSet& readset) {
    assert(&readset.mpopi() == &mpopi_); // Otherwise sample indexes may be misleading.

    ++loc_i_;
    if (loc_i_ == int32_t(n_loci())) {
        assert(eof_);
        return false;
    }

    readset.clear();
    readset.bam_i(loc_i_);
    readset.id(atoi(bam_f_->h().chrom_str(loc_i_)));
    if (!loci_aln_positions_.empty()) {
        const char* p = loci_aln_positions_.at(loc_i_);
        const char* q = p;
        while (*q && *q != '\n' && *q != '\t')
            ++q;
        readset.pos(PhyLoc(string(p, q)));
    }

    if (!eof_) {
        // Read all the reads of the locus, and one more.
        if (bam_f_->r().chrom() < loc_i_) {
            cerr << "Error: BAM file isn't properly sorted.\n";
            throw exception();
        }
        while (bam_f_->r().chrom() == loc_i_) {
            const BamRecord& rec = bam_f_->r();
            const char* rg = rec.read_group();
            if (rg == NULL) {
                cerr << "Error: Corrupted BAM file: missing read group.\n";
                throw exception();
            }
            if (rec.is_read2())
                readset.add_pe(SRead(Read(rec.seq(), rec.qname()+"/2"), rg_to_sample_.at(rg)));
            else if (rec.is_read1())
                readset.add(SRead(Read(rec.seq(), rec.qname()+"/1"), rg_to_sample_.at(rg)));
            else
                // If tsv2bam wasn't given paired-end reads, no flag was set and the
                // read names were left unchanged, so we also don't touch them.
                readset.add(SRead(Read(rec.seq(), rec.qname()), rg_to_sample_.at(rg)));

            if (!bam_f_->next_record()) {
                eof_ = true;
                break;
            }
        }
    }

    if (readset.reads().empty()) {
        // A catalog locus might not have reads, e.g. if some samples were
        // omitted in the `samtools merge`, but this would be weird.
        static bool warn = false;
        if (warn) {
            warn = false;
            cerr << "Warning: Some catalog loci do not have any reads in the BAM file.\n";
        }
    }
    return true;
}

inline
BamCLocBuilder::BamCLocBuilder(Bam** bam_f)
:
    bam_f_(*bam_f),
    n_loci_built_(0),
    eof_(false),
    next_record_(PhyLoc(), SAlnRead(AlnRead(Read(DNASeq4(), string()), Cigar()), SIZE_MAX))
{
    *bam_f = NULL;

    // Create the MetaPopInfo object from the header.
    build_mpopi(mpopi_, rg_to_sample_, *bam_f_);

    // Read the very first record.
    if (!read_and_parse_next_record()) {
        cerr << "Error: No records in BAM file '" << bam_f_->path << "'.\n";
        throw exception();
    }
}

inline
bool
BamCLocBuilder::read_and_parse_next_record()
{
    if(!bam_f_->next_record())
        return false;

    const BamRecord& r = bam_f_->r();

    int32_t chrom = r.chrom();
    int32_t pos = r.pos();
    strand_type strand = strand_plus;
    string name = r.qname();
    Cigar cigar = r.cigar();
    DNASeq4 seq = r.seq();
    size_t sample = rg_to_sample_.at(string(r.read_group()));

    if (r.is_rev_compl()) {
        pos += cigar_length_ref(cigar) - 1;
        strand = strand_minus;
        std::reverse(cigar.begin(), cigar.end());
        seq = seq.rev_compl();
    }

    next_record_ = {
        PhyLoc(bam_f_->h().chrom_str(chrom), pos, strand),
        SAlnRead(AlnRead(Read(move(seq), move(name)), move(cigar)), sample)
    };
    return true;
}

inline
bool
BamCLocBuilder::fill_window()
{
    if (!eof_) {
        while (fw_reads_by_5prime_pos.empty()
                || (strcmp(next_record_.first.chr(), fw_reads_by_5prime_pos.begin()->first.chr()) == 0
                    && next_record_.first.bp <= fw_reads_by_5prime_pos.begin()->first.bp + max_insert_refsize_)
                ) {

            vector<SAlnRead>& v = fw_reads_by_5prime_pos.insert(
                    std::make_pair(next_record_.first, vector<SAlnRead>()) // (`make_pair` is necessary here until c++17.)
                    ).first->second;
            v.push_back(move(next_record_.second));

            if (!read_and_parse_next_record()) {
                eof_ = true;
                break;
            }
        }
    }

    return !fw_reads_by_5prime_pos.empty();
}

inline
bool
BamCLocBuilder::build_one_locus(CLocAlnSet& aln_loc)
{
    aln_loc.clear();

    //
    // Find the next locus.
    //
    while(true) {
        // Fill the window.
        if (!fill_window())
            return false;

        //
        // Apply filters to the putative locus.
        //
        assert(!fw_reads_by_5prime_pos.empty());
        vector<SAlnRead>& loc_reads = fw_reads_by_5prime_pos.begin()->second;

        // Remove reads from samples that have less that `min_reads_per_sample_` reads.
        vector<size_t> n_reads_per_sample (mpopi_.samples().size(), 0);
        for (const SAlnRead& read : loc_reads)
            ++n_reads_per_sample[read.sample];
        loc_reads.erase(std::remove_if(
            loc_reads.begin(), loc_reads.end(),
            [&n_reads_per_sample] (const SAlnRead& read) {return n_reads_per_sample[read.sample] < min_reads_per_sample_;}
            ), loc_reads.end());

        if (loc_reads.empty()) {
            // Discard the locus; regenerate the window and retry.
            fw_reads_by_5prime_pos.erase(fw_reads_by_5prime_pos.begin());
            continue;
        }
        break;
    }

    //
    // Build the locus object.
    //
    aln_loc.reinit(n_loci_built_+1, fw_reads_by_5prime_pos.begin()->first, &mpopi_);
    for (SAlnRead& read : fw_reads_by_5prime_pos.begin()->second)
        aln_loc.add(move(read));
    fw_reads_by_5prime_pos.erase(fw_reads_by_5prime_pos.begin());
    ++n_loci_built_;
    return true;
}

inline
VcfCLocReader::VcfCLocReader(const string& vcf_path)
        : vcf_f_(vcf_path),
          next_rec_(),
          eof_(false)
{
    // Read the very first record.
    if(!vcf_f_.next_record(next_rec_))
        eof_ = true;
}

inline int
VcfCLocReader::open(string &vcf_path)
{
    this->vcf_f_.open(vcf_path);
    
    // Read the very first record.
    if(!this->vcf_f_.next_record(next_rec_))
        eof_ = true;

    return 0;
}

inline
void VcfCLocReader::set_sample_ids(MetaPopInfo& mpopi) const {
    //TODO Write & retrieve actual sample IDs using the VCF header.
    for (size_t i = 0; i < mpopi.samples().size(); ++i)
        mpopi.set_sample_id(i, i+1); //id=i+1
}

inline
bool VcfCLocReader::read_one_locus(vector<VcfRecord>& records) {
    records.clear();
    if (eof_)
        return false;

    // Parse the locus ID.
    static string curr_chrom;
    curr_chrom.assign(next_rec_.chrom());
    records.push_back(move(next_rec_));

    // Read all the records of the locus, and one more.
    while (strcmp(records.back().chrom(), curr_chrom.c_str()) == 0) {
        records.push_back(VcfRecord());
        if (!vcf_f_.next_record(records.back())) {
            eof_ = true;
            break;
        }
    }

    // Remove the first record of the next locus (or the record for which EOF
    // was encountered) from the vector.
    next_rec_ = move(records.back());
    records.pop_back();

    return true;
}

#endif
