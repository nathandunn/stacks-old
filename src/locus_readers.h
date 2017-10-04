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
    BamCLocBuilder(Bam** bam_f, size_t min_reads_per_sample_, bool paired_mode, size_t max_insert_refsize=0);
    ~BamCLocBuilder() {if (bam_f_) delete bam_f_;}

    // Reads one locus. Returns false on EOF.
    bool build_one_locus(CLocAlnSet& readset);    
    
    const MetaPopInfo& mpopi() const {return mpopi_;}
    const Bam* bam_f() const {return bam_f_;}

private:
    size_t min_reads_per_sample_ = 3;
    bool paired_mode_;
    size_t max_insert_refsize_;

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
    std::multimap<PhyLoc,SAlnRead> pe_reads_by_5prime_pos;
    map<const char*,std::multimap<PhyLoc,SAlnRead>::iterator, LessCStrs> pe_reads_by_name;
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
BamCLocBuilder::BamCLocBuilder(
        Bam** bam_f,
        size_t min_reads_per_sample,
        bool paired_mode,
        size_t max_insert_refsize
) :
    min_reads_per_sample_(min_reads_per_sample),
    paired_mode_(paired_mode),
    max_insert_refsize_(paired_mode ? max_insert_refsize_ : 0),
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

    if (r.is_read1())
        name += "/1";
    else if (r.is_read2())
        name += "/2";
    cigar_simplify_to_MDI(cigar);
        
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
                || (
                    // This defines the sliding window. On the right side we use the left boundary
                    // of the alignment (i.e. the BAM "pos" field) rather than the 5prime bp of the
                    // next record--this is less stringent, and prevents abberant CIGARs to cause
                    // issues.
                    strcmp(next_record_.first.chr(), fw_reads_by_5prime_pos.begin()->first.chr()) == 0
                    && bam_f_->r().pos() <= fw_reads_by_5prime_pos.begin()->first.bp + max_insert_refsize_
                )) {
            
            bool treat_as_fw;
            const string& name = next_record_.second.name;
            if (!paired_mode_) {
                treat_as_fw = true;
            } else if (name.back() == '1' && (*++name.rbegin() == '/' || *++name.rbegin() == '_')) {
                // n.b. The SAlnRead name has been edited according to the BAM flags.
                // We don't discriminate the case when the READ1/READ2 flags are set
                // and the case when these flags are not set but the read names end
                // with /1, /2.
                treat_as_fw = true;
            } else if (name.back() == '2' && (*++name.rbegin() == '/' || *++name.rbegin() == '_')) {
                treat_as_fw = false;
            } else {
                cerr << "Error: --paired was specified but it is unclear whether BAM record '"
                     << bam_f_->r().qname() << "' corresponds to a forward or reverse read.\n";
                throw exception();
            }

            if (treat_as_fw) {
                vector<SAlnRead>& v = fw_reads_by_5prime_pos.insert(
                        std::make_pair(next_record_.first, vector<SAlnRead>()) // (`make_pair` is necessary here until c++17.)
                        ).first->second;
                v.push_back(move(next_record_.second));                
            } else {
                auto itr = pe_reads_by_5prime_pos.insert(move(next_record_));
                pe_reads_by_name.insert( {itr->second.name.c_str(), itr} );
            }

            if (!read_and_parse_next_record()) {
                eof_ = true;
                break;
            }
        }
    }

    // Cleanup paired-end reads that are out of the window (for paired-end reads,
    // the windows extends in both directions).
    if (!fw_reads_by_5prime_pos.empty()) {
        const PhyLoc& window_center = fw_reads_by_5prime_pos.begin()->first;
        while (!pe_reads_by_5prime_pos.empty()
                && (strcmp(pe_reads_by_5prime_pos.begin()->first.chr(), window_center.chr()) != 0
                    || pe_reads_by_5prime_pos.begin()->first.bp + max_insert_refsize_ < window_center.bp + 1)
                ) {
            pe_reads_by_name.erase(pe_reads_by_5prime_pos.begin()->second.name.c_str());
            pe_reads_by_5prime_pos.erase(pe_reads_by_5prime_pos.begin());
        }
        // Also check at the end of the multimap, as the BAM chromosome order isn't
        // necessarily alphabetic (while the PhyLoc order is).
        while (!pe_reads_by_5prime_pos.empty()
                && strcmp(pe_reads_by_5prime_pos.rbegin()->first.chr(), window_center.chr()) != 0
                ) {
            pe_reads_by_name.erase(pe_reads_by_5prime_pos.rbegin()->second.name.c_str());
            pe_reads_by_5prime_pos.erase(--pe_reads_by_5prime_pos.end());
        }
    } else {
        pe_reads_by_5prime_pos.clear();
        pe_reads_by_name.clear();
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

        // Apply filters to the putative locus.
        // i.e. Remove reads from samples that have less that `min_reads_per_sample_` reads.
        assert(!fw_reads_by_5prime_pos.empty());
        vector<SAlnRead>& loc_reads = fw_reads_by_5prime_pos.begin()->second;

        vector<size_t> n_reads_per_sample (mpopi_.samples().size(), 0);
        for (const SAlnRead& read : loc_reads)
            ++n_reads_per_sample[read.sample];
        loc_reads.erase(std::remove_if(
            loc_reads.begin(), loc_reads.end(),
            [&n_reads_per_sample,this] (const SAlnRead& read) {return n_reads_per_sample[read.sample] < min_reads_per_sample_;}
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

    // Forward reads.
    vector<SAlnRead> fw_reads = move(fw_reads_by_5prime_pos.begin()->second);

    // Paired-end reads.
    vector<SAlnRead> pe_reads;
    if (!pe_reads_by_5prime_pos.empty()) {
        const PhyLoc& fw_5prime = fw_reads_by_5prime_pos.begin()->first;
        string pe_name;
        for (const SAlnRead& fw_read : fw_reads) {
            assert(fw_read.name.back() == '1');
            pe_name = fw_read.name;
            pe_name.back() = '2';
            auto pe_name_itr = pe_reads_by_name.find(pe_name.c_str());
            if (pe_name_itr != pe_reads_by_name.end()) {
                const PhyLoc& pe_5prime = pe_name_itr->second->first;
                SAlnRead& pe_read = pe_name_itr->second->second;

                size_t fw_len = cigar_length_ref(fw_read.cigar);
                size_t pe_len = cigar_length_ref(pe_read.cigar);

                // Check that the relative positions of the reads are compatible.
                assert(strcmp(fw_5prime.chr(), pe_5prime.chr()) == 0); // The window.
                if (
                        (fw_5prime.strand == strand_plus && pe_5prime.strand == strand_minus
                            && fw_5prime.bp + std::max(fw_len, pe_len) <= pe_5prime.bp + 1
                            && fw_5prime.bp + max_insert_refsize_ >= pe_5prime.bp + 1)
                        ||  
                        (fw_5prime.strand == strand_minus && pe_5prime.strand == strand_plus
                            && pe_5prime.bp + std::max(fw_len, pe_len) <= fw_5prime.bp + 1
                            && pe_5prime.bp + max_insert_refsize_ >= fw_5prime.bp + 1)
                        ) {
                    // Put the paired read in the reference of the forward reads.
                    pe_read.seq = pe_read.seq.rev_compl();
                    std::reverse(pe_read.cigar.begin(), pe_read.cigar.end());

                    size_t extend;
                    if (fw_5prime.strand == strand_plus)
                        extend = pe_5prime.bp - fw_5prime.bp + 1 - pe_len;
                    else
                        extend = fw_5prime.bp - pe_5prime.bp + 1 - pe_len;
                    cigar_extend_left(pe_read.cigar, extend);
                    
                    pe_reads.push_back(move(pe_read));
                }

                pe_reads_by_5prime_pos.erase(pe_name_itr->second);
                pe_reads_by_name.erase(pe_name_itr);
            }
        }
    }

    // Determine the locus length.
    size_t max_cig_length = 0;
    for (const SAlnRead& read : fw_reads) {
        size_t l = cigar_length_ref(read.cigar);
        if (l > max_cig_length)
            max_cig_length = l;
    }
    for (const SAlnRead& read : pe_reads) {
        size_t l = cigar_length_ref(read.cigar);
        if (l > max_cig_length)
            max_cig_length = l;
    }

    // Build the actual object.
    aln_loc.reinit(n_loci_built_+1, fw_reads_by_5prime_pos.begin()->first, &mpopi_);
    aln_loc.ref(DNASeq4(max_cig_length));
    for (SAlnRead& read : fw_reads) {
        cigar_extend_right(read.cigar, max_cig_length - cigar_length_ref(read.cigar));
        aln_loc.add(move(read));
    }
    for (SAlnRead& read : pe_reads) {
        cigar_extend_right(read.cigar, max_cig_length - cigar_length_ref(read.cigar));
        aln_loc.add(move(read));
    }

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
