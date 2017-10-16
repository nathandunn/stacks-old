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
    struct Config {
        bool paired;
        size_t max_insert_refsize;
        size_t min_mapq;
        double max_clipped;
        size_t min_reads_per_sample;
        bool ign_pe_reads;
    };

    struct BamStats {
        size_t n_records;
        size_t n_primary;
        size_t n_primary_mapq;
        size_t n_primary_softclipped;
        size_t n_secondary;
        size_t n_supplementary;
        size_t n_unmapped;

        size_t n_primary_kept() const {return n_primary - n_primary_mapq - n_primary_softclipped;}
    };

    BamCLocBuilder(Bam** bam_f, const Config& cfg);
    ~BamCLocBuilder() {if (bam_f_) delete bam_f_;}

    // Reads one locus. Returns false on EOF.
    bool build_one_locus(CLocAlnSet& readset);    

    const MetaPopInfo& mpopi() const {return mpopi_;}
    const Bam* bam_f() const {return bam_f_;}
    const BamStats& stats() const {return stats_;}

private:
    // Structure to store the reads' 5prime positions.
    struct Pos5 {
        int32_t chrom;
        int32_t bp;
        bool    strand; // Plus is true, minus is false.

        bool operator<(const Pos5& other) const {
            if (chrom == other.chrom) {
                if (bp == other.bp) {
                    if (strand == other.strand)
                        return false; // Equal.
                    return strand == false; // Minus strand first.
                }
                return bp < other.bp;
            }
            return chrom < other.chrom;
        }
    };

    Config cfg_;

    Bam* bam_f_;
    MetaPopInfo mpopi_;
    map<string, size_t> rg_to_sample_;

    BamStats stats_;
    size_t n_loci_built_;

    BamRecord next_record_;
    bool treat_next_record_as_fw_;
    bool next_record();

    // Buffers to store the reads of the sliding window.
    bool fill_window();
    bool eof_;
    map<Pos5,vector<BamRecord>> fw_reads_by_5prime_pos_;
    std::multimap<Pos5,BamRecord> pe_reads_by_5prime_pos_;
    map<const char*,std::multimap<Pos5,BamRecord>::iterator, LessCStrs> pe_reads_by_name_;
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

                if (rg_ids.back().empty() || samples.back().empty())
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

    if (samples.empty()) {
        cerr << "Error: Found no @RG lines (read groups/sample identifiers) in the header of BAM file '"
             << bam_f.path << "'.\n";
        throw exception();
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
    for (const Sample& s : mpopi_.samples()) {
        if (s.id == -1) {
            cerr << "Error: Sample '" << s.name << "' is missing its (integer) ID; was tsv2bam run correctly?\n";
            throw exception();
        }
    }

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
                cerr << "Error: Corrupted BAM file: missing read group for record '" << rec.qname() << "'.\n";
                throw exception();
            }
            if (rec.is_read2())
                readset.add_pe(SRead(Read(rec.seq(), string(rec.qname())+"/2"), rg_to_sample_.at(rg)));
            else if (rec.is_read1())
                readset.add(SRead(Read(rec.seq(), string(rec.qname())+"/1"), rg_to_sample_.at(rg)));
            else
                // If tsv2bam wasn't given paired-end reads, no flag was set and the
                // read names were left unchanged, so we also don't touch them.
                readset.add(SRead(Read(rec.seq(), string(rec.qname())), rg_to_sample_.at(rg)));

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
        const Config& cfg
) :
    cfg_(cfg),
    bam_f_(*bam_f),
    stats_(),
    n_loci_built_(0),
    treat_next_record_as_fw_(false),
    eof_(false)
    {
    *bam_f = NULL;

    if (!cfg_.paired)
        cfg_.max_insert_refsize = 0;

    // Create the MetaPopInfo object from the header. Assign sample IDs.
    build_mpopi(mpopi_, rg_to_sample_, *bam_f_);
    for (size_t i=0; i<mpopi_.samples().size(); ++i)
        mpopi_.set_sample_id(i, i+1);

    // Read the very first record.
    if (!next_record()) {
        cerr << "Error: No usable records in BAM file '" << bam_f_->path << "'.\n";
        throw exception();
    }
}

inline
bool
BamCLocBuilder::next_record()
{
    while (true) {
        if(!bam_f_->next_record())
            return false;

        const BamRecord& r = bam_f_->r();

        if (cfg_.ign_pe_reads && r.is_read2())
            continue;

        ++stats_.n_records;

        // Check if the record is primary.
        if (r.is_unmapped()) {
            ++stats_.n_unmapped;
            continue;
        } else if (r.is_secondary()) {
            ++stats_.n_secondary;
            continue;
        } else if (r.is_supplementary()) {
            ++stats_.n_supplementary;
            continue;
        }
        assert(r.is_primary());
        ++stats_.n_primary;

        // Check the MAPQ
        if (r.mapq() < cfg_.min_mapq) {
            ++stats_.n_primary_mapq;
            continue;
        }

        // Assess whether this alignment should be treated as a forward read.
        // N.B. We don't discriminate the case when the READ1/READ2 flags are
        // set and the case when these flags are not set but the read names end
        // with /1, /2.
        if (!cfg_.paired) {
            treat_next_record_as_fw_ = true;
        } else if (r.is_read1()) {
            treat_next_record_as_fw_ = true;
        } else if (r.is_read2()) {
            treat_next_record_as_fw_ = false;
        } else {
            cerr << "Error: BAM record '" << r.qname()
                 << "' is not flagged as READ1 nor READ2 (while reading BAM file '"
                 << bam_f_->path << "'in paired mode).\n";
            throw exception();
        }

        // Check that there's a sequence.
        if (r.hts_l_seq() <= 0) {
            cerr << "Error: No sequence at BAM primary record '" << r.qname() << "'.\n";
            throw exception();
        }

        // Check the CIGAR.
        if (r.hts_n_cigar() == 0) {
            cerr << "Error: Empty CIGAR at BAM primary record '" << r.qname() << "'.\n";
            throw exception();
        }
        if (treat_next_record_as_fw_) {
            if (r.is_rev_compl()) {
                if(BamRecord::cig_op_t(r.hts_cigar()[r.hts_n_cigar()-1]) != BAM_CMATCH) {
                    // A cutsite is expected, and it is not aligned.
                    ++stats_.n_primary_softclipped;
                    continue;
                }
            } else {
                if(BamRecord::cig_op_t(r.hts_cigar()[0]) != BAM_CMATCH) {
                    ++stats_.n_primary_softclipped;
                    continue;
                }
            }
        }
        size_t softclipped = 0;
        for (size_t i=0; i<r.hts_n_cigar(); ++i)
            if (BamRecord::cig_op_t(r.hts_cigar()[i]) == BAM_CSOFT_CLIP)
                softclipped += BamRecord::cig_op_len(r.hts_cigar()[i]);
        // xxx Terminal insertions should count as soft-clipping.
        if ((double) softclipped / r.hts_l_seq() > cfg_.max_clipped) {
            // Too much soft-clipping.
            ++stats_.n_primary_softclipped;
            continue;
        }

        // Assign to `next_record_`.
        next_record_ = move(bam_f_->r());

        break;
    }
    return true;
}

inline
bool
BamCLocBuilder::fill_window()
{
    /*
     * This defines the behavior of the sliding window.
     *
     * We guarantee that any usable alignment of which the leftmost base
     * is within `max_insert_size` bases of the first base of the leftmost
     * cutsite is loaded.
     *
     * `fw_reads_by_5prime_pos_` is a map of {[first cutsite base pos] : BamRecord}
     * Cutsites can be on the left or right of the read alignment depending on the
     * strand, but this is easy to handle since `BamCLocBuilder::next_record()`
     * ensures that the BAM record's TLEN field is properly set.
     *
     * `pe_reads_by_5prime_pos_` also uses the 5prime position (i.e. that has been
     * tlen-corrected if the strand is minus) so that we can compute insert lengths
     * properly later. Paired-end reads need to be sorted by position so that we
     * can discard those that were not matched and have become too far behind to
     * find one at the given `max_insert_size`.
     */

    if (!eof_) {
        const BamRecord* previous = NULL;
        while (fw_reads_by_5prime_pos_.empty()
                || (
                    next_record_.chrom() == fw_reads_by_5prime_pos_.begin()->first.chrom
                    && size_t(next_record_.pos()) <= fw_reads_by_5prime_pos_.begin()->first.bp + cfg_.max_insert_refsize
                )
        ) {
            // Determine the 5prime position.
            Pos5 rec_5prime_pos;
            rec_5prime_pos.chrom = next_record_.chrom();
            if (next_record_.is_rev_compl()) {
                rec_5prime_pos.bp = next_record_.pos() + BamRecord::cig_ref_len(next_record_) - 1;
                rec_5prime_pos.strand = false;
            } else {
                rec_5prime_pos.bp = next_record_.pos();
                rec_5prime_pos.strand = true;
            }

            // Save the record.
            if (treat_next_record_as_fw_) {
                auto locus_itr = fw_reads_by_5prime_pos_.insert(
                        // (Cannot use {first, second} here until c++17.)
                        std::make_pair(move(rec_5prime_pos), vector<BamRecord>())
                        ).first;
                locus_itr->second.push_back(move(next_record_));
                previous = &locus_itr->second.back(); // At least stable until we push_back again.
            } else {
                auto pe_read_itr = pe_reads_by_5prime_pos_.insert(
                        std::make_pair(move(rec_5prime_pos), move(next_record_))
                        );
                pe_reads_by_name_.insert( {pe_read_itr->second.qname(), pe_read_itr} );
                previous = &pe_read_itr->second;
            }

            if (!next_record()) {
                eof_ = true;
                break;
            }
            if (next_record_.chrom() < previous->chrom()
                    || (next_record_.pos() < previous->pos() && next_record_.chrom() == previous->chrom())
            ) {
                cerr << "Error: BAM file is not properly sorted (record '" << next_record_.qname()
                     << "' at " << bam_f_->h().chrom_str(next_record_.chrom()) << ':' << next_record_.pos()
                     << " should come before record '" << previous->qname() << "' at "
                     << bam_f_->h().chrom_str(previous->chrom()) << ':' << previous->pos() << ").\n";
                throw exception();
            }
        }
    }
    if (fw_reads_by_5prime_pos_.empty())
        // No more loci.
        return false;

    if (cfg_.paired) {
        // Clean up paired-end reads that are too far behind.
        Pos5 leftmost_cutsite = fw_reads_by_5prime_pos_.begin()->first;
        while (!pe_reads_by_5prime_pos_.empty()
                && (
                    pe_reads_by_5prime_pos_.begin()->first.chrom != leftmost_cutsite.chrom
                    || pe_reads_by_5prime_pos_.begin()->first.bp + cfg_.max_insert_refsize < size_t(leftmost_cutsite.bp) + 1
                )
        ) {
            pe_reads_by_name_.erase(pe_reads_by_5prime_pos_.begin()->second.qname());
            pe_reads_by_5prime_pos_.erase(pe_reads_by_5prime_pos_.begin());
        }
        // Assuming proper sorting (which we enforce) there shouldn't be reads
        // from further-ranked chromosomes.
        assert(pe_reads_by_5prime_pos_.empty() || pe_reads_by_5prime_pos_.rbegin()->first.chrom == leftmost_cutsite.chrom);
    } else {
        assert(pe_reads_by_5prime_pos_.empty());
    }

    return true;
}

inline
bool
BamCLocBuilder::build_one_locus(CLocAlnSet& aln_loc)
{
    aln_loc.clear();

    //
    // Find the next locus.
    //

    vector<size_t> fw_samples;
    while(true) {
        //
        // Fill the window.
        //
        if (!fill_window())
            return false;
        assert(!fw_reads_by_5prime_pos_.empty());

        //
        // Apply filters to the putative locus.
        //

        vector<BamRecord>& loc_records = fw_reads_by_5prime_pos_.begin()->second;

        // Get the samples of the records, and tally sample occurrences.
        fw_samples.clear();
        vector<size_t> n_reads_per_sample (mpopi_.samples().size(), 0);
        for (const BamRecord& r : loc_records) {
            const char* rg = r.read_group();
            if (rg == NULL) {
                cerr << "Error: Corrupted BAM file: missing read group for record '" << r.qname() << "'.\n";
                throw exception();
            }
            size_t sample = rg_to_sample_.at(string(rg));
            fw_samples.push_back(sample);
            ++n_reads_per_sample[sample];
        }

        // Remove reads from samples that have less than `cfg_.min_reads_per_sample` reads.
        assert(fw_samples.size() == loc_records.size());
        for (size_t i=0; i<loc_records.size(); ++i) {
            if (n_reads_per_sample[fw_samples[i]] < cfg_.min_reads_per_sample) {
                // Mark it for deletion.
                loc_records[i].destroy();
                assert(loc_records[i].empty());
                fw_samples[i] = SIZE_MAX;
            }
        }
        // Splice deleted records out of the vectors.
        loc_records.erase(std::remove_if(
                loc_records.begin(), loc_records.end(),
                [] (const BamRecord& r) {return r.empty();}
                ), loc_records.end());
        fw_samples.erase(std::remove_if(
                fw_samples.begin(), fw_samples.end(),
                [] (size_t s) {return s == SIZE_MAX;}
                ), fw_samples.end());

        if (loc_records.empty()) {
            // Discard the locus; regenerate the window and retry.
            fw_reads_by_5prime_pos_.erase(fw_reads_by_5prime_pos_.begin());
            continue;
        }
        break;
    }

    //
    // Build the locus object.
    //

    // Position.
    Pos5 cutsite = fw_reads_by_5prime_pos_.begin()->first;
    PhyLoc loc_pos (bam_f_->h().chrom_str(cutsite.chrom), cutsite.bp, cutsite.strand ? strand_plus : strand_minus);
    size_t max_ref_span = 0;

    // Forward reads.
    const vector<BamRecord>& fw_records = fw_reads_by_5prime_pos_.begin()->second;
    vector<SAlnRead> fw_reads;
    assert(fw_samples.size() == fw_records.size());
    for (size_t i=0; i<fw_records.size(); ++i) {
        const BamRecord& r = fw_records[i];
        size_t sample = fw_samples[i];

        string name (r.qname());
        if (r.is_read1())
            name += "/1";
        else if (r.is_read2())
            // (Happens when --paired is off.)
            name += "/2";
        DNASeq4 seq = r.seq();
        Cigar cigar = r.cigar();
        cigar_simplify_to_MDI(cigar);

        size_t ref_span;
        if (cutsite.strand == false) {
            assert(r.is_rev_compl());
            std::reverse(cigar.begin(), cigar.end());
            seq = seq.rev_compl();

            assert(r.pos() <= cutsite.bp);
            ref_span = cutsite.bp - r.pos() + 1;
        } else {
            ref_span = BamRecord::cig_ref_len(r);
        }
        if (ref_span > max_ref_span)
            max_ref_span = ref_span;

        fw_reads.push_back(SAlnRead(AlnRead(Read(move(seq), move(name)), move(cigar)), sample));
    }

    // Paired-end reads.
    vector<SAlnRead> pe_reads;
    if (!pe_reads_by_5prime_pos_.empty()) {
        for (size_t i=0; i<fw_records.size(); ++i) {
            const BamRecord& fw_rec = fw_records[i];

            // Find a mate.
            auto pe_name_itr = pe_reads_by_name_.find(fw_rec.qname());
            if (pe_name_itr == pe_reads_by_name_.end())
                continue;

            std::multimap<Pos5,BamRecord>::iterator pe_rec_itr = pe_name_itr->second;
            Pos5 pe_pos5 = pe_rec_itr->first;
            const BamRecord& pe_rec = pe_rec_itr->second;

            // Check that the mate is in the window and that the reads don't
            // extend past one another. //xxx Oct 2017 @Nick: This may actually
            // be okay (c.f. the "clocaln::juxtapose" assert failures upon assembly
            // of the barcode behind READ1), but should we allow it?
            assert(pe_pos5.chrom == cutsite.chrom); // By construction of the window.
            if (pe_pos5.strand == cutsite.strand)
                continue;
            size_t pe_ref_len;
            if (cutsite.strand == true) {
                assert(pe_pos5.bp >= pe_rec.pos()); // The pe read is reverse complemented, so its 5prime is on the right.
                pe_ref_len = pe_pos5.bp - pe_rec.pos() + 1;
                if (pe_rec.pos() < cutsite.bp
                        || size_t(pe_pos5.bp - cutsite.bp + 1) > cfg_.max_insert_refsize)
                    continue;
            } else {
                assert(pe_pos5.bp == pe_rec.pos());
                pe_ref_len = BamRecord::cig_ref_len(pe_rec);
                if (pe_rec.pos() + pe_ref_len - 1 > size_t(cutsite.bp)
                        || size_t(cutsite.bp - pe_pos5.bp + 1) > cfg_.max_insert_refsize)
                    continue;
            }
            size_t insert_length = std::labs(pe_pos5.bp -cutsite.bp) + 1;
            if (insert_length > max_ref_span)
                max_ref_span = insert_length;

            // Check that the read groups are consistent.
            const char* rg = pe_rec.read_group();
            if (rg == NULL) {
                cerr << "Error: Corrupted BAM file: missing read group for record '" << pe_rec.qname() << "'.\n";
                throw exception();
            }
            size_t pe_sample = rg_to_sample_.at(string(rg));
            if (pe_sample != fw_samples[i]) {
                cerr << "Warning: Paired reads '" << fw_rec.qname() << "' and '"
                     << pe_rec.qname() << "' belong to different samples ('"
                     << mpopi_.samples()[fw_samples[i]].name << "' and '"
                     << mpopi_.samples()[pe_sample].name << "').\n";
                continue;
            }

            // Save the paired-end read.
            // We need to put the paired-end read in the same reference as the
            // forward reads. (1) BAM records that are part of a locus on the
            // minus strand should be reverse complemented.
            // (2) The first reference base in the cigar should be `cutsite.bp`,
            // regardless of the direction of the locus.
            string name (pe_rec.qname());
            assert(pe_rec.is_read2());
            name += "/2";
            DNASeq4 seq = pe_rec.seq();
            Cigar cigar = pe_rec.cigar();
            cigar_simplify_to_MDI(cigar);
            if (cutsite.strand == true) {
                assert(pe_rec.pos() >= cutsite.bp);
                cigar_extend_left(cigar, pe_rec.pos() - cutsite.bp);
            } else {
                assert(size_t(cutsite.bp) >= pe_rec.pos() + pe_ref_len - 1);
                cigar_extend_right(cigar, cutsite.bp - (pe_rec.pos() + pe_ref_len - 1));
                std::reverse(cigar.begin(), cigar.end());
                seq = seq.rev_compl();
            }
            pe_reads.push_back(SAlnRead(AlnRead(Read(move(seq), move(name)), move(cigar)), pe_sample));

            pe_reads_by_name_.erase(pe_name_itr);
            pe_reads_by_5prime_pos_.erase(pe_rec_itr);
        }
    }

    // Build the actual object.
    aln_loc.reinit(n_loci_built_++, loc_pos, &mpopi_);
    aln_loc.ref(DNASeq4(max_ref_span));
    for (SAlnRead& read : fw_reads) {
        cigar_extend_right(read.cigar, max_ref_span - cigar_length_ref(read.cigar));
        aln_loc.add(move(read));
    }
    for (SAlnRead& read : pe_reads) {
        cigar_extend_right(read.cigar, max_ref_span - cigar_length_ref(read.cigar));
        aln_loc.add(move(read));
    }

    fw_reads_by_5prime_pos_.erase(fw_reads_by_5prime_pos_.begin());
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
    //xxx We could write & retrieve actual sample IDs using the VCF header.
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
