#ifndef LOCUS_READERS_H
#define LOCUS_READERS_H

#include "constants.h"
#include "MetaPopInfo.h"
#include "BamI.h"
#include "locus.h"
#include "Vcf.h"

// BamPopInfo : class to handle BAM read groups and interface
// them with MetaPopInfo Samples.
class BamPopInfo {
    MetaPopInfo mpopi_;
    typedef map<string,size_t> RGMap;
    vector<RGMap> read_groups_;
    vector<map<const char*, RGMap::iterator, LessCStrs>> rg_to_sample_;

public:
    BamPopInfo(const vector<Bam*>& bam_fs); //, const vector<size_t>& samples={});
    void reassign_ids();

    const MetaPopInfo& mpopi() const {return mpopi_;}
    size_t sample_of(const BamRecord& rec, size_t bam_f_i=0);
};

// BamCLocReader: Class to read catalog loci from files by the de novo pipeline.
class BamCLocReader {
    vector<Bam*> bam_fs_;

    BamPopInfo bpopi_;
    vector<const char*> loci_aln_positions_;

    int32_t loc_i_; // Index of the next locus (chromosome). Incremented by read_one_locus().

    vector<uchar> eofs_; // For each BAM file. (Vector of bools; `uchar`s because std specializes `vector<bool>`.)
    vector<BamRecord> next_records_; // The next record, for each file.

public:
    BamCLocReader(vector<Bam*>&& bam_fs);
    ~BamCLocReader() {for (Bam* bam_f : bam_fs_) delete bam_f;}

    const vector<Bam*>& bam_fs() const {return bam_fs_;}
    size_t n_loci() const {return bam_fs_[0]->h().n_ref_chroms();}
    int target2id(size_t target_i) const {return atoi(bam_fs_[0]->h().chrom_str(target_i));}
    const MetaPopInfo& mpopi() const {return bpopi_.mpopi();}

    // Reads one locus. Returns false on EOF.
    bool read_one_locus(CLocReadSet& readset);
};

// BamCLocBuilder: Class to build catalog loci from third party alignment files.
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

        size_t n_loci_built;

        size_t n_primary_kept() const {return n_primary - n_primary_mapq - n_primary_softclipped;}
    };

    BamCLocBuilder(vector<Bam*>&& bam_fs, const Config& cfg);
    ~BamCLocBuilder() {for(Bam* bam_f : bam_fs_) delete bam_f;}

    // Reads one locus. Returns false on EOF.
    bool build_one_locus(CLocAlnSet& readset);    

    const MetaPopInfo& mpopi() const {return bpopi_.mpopi();}
    const vector<Bam*>& bam_fs() const {return bam_fs_;}
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
    BamStats stats_;

    Bam* bam_f_;
    vector<Bam*> bam_fs_;
    BamPopInfo bpopi_;

    bool next_record(size_t bam_f_i);
    bool next_record() {return next_record(0);}
    vector<uchar> eofs_; // For each BAM file. (Vector of bools; `uchar`s because std specializes `vector<bool>`.)
    vector<BamRecord> next_records_; // The next record, for each file.
    vector<uchar> treat_next_records_as_fw_; // Whether to treat the above as forward reads. c.f. `cfg_.paired`/--paired.

    // Buffers to store the reads of the sliding window.
    bool fill_window();
    void add_next_record_to_the_window(size_t bam_f_i);
    typedef size_t SampleIdx;
    map<Pos5,vector<pair<BamRecord,SampleIdx>>> fw_reads_by_5prime_pos_;
    std::multimap<Pos5,pair<BamRecord,SampleIdx>> pe_reads_by_5prime_pos_;
    map<const char*,
        std::multimap<Pos5,pair<BamRecord,SampleIdx>>::iterator,
        LessCStrs> pe_reads_by_name_;
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
BamPopInfo::BamPopInfo(const vector<Bam*>& bam_fs)
:
    read_groups_(bam_fs.size()),
    rg_to_sample_(bam_fs.size())
{
    if (bam_fs.empty())
        DOES_NOT_HAPPEN;

    struct SampleData {
        string rg;
        string name;
        int stacks_id;
    };
    vector<vector<SampleData>> readgroups_per_f  (bam_fs.size());

    for (size_t bam_f_i=0; bam_f_i<bam_fs.size(); ++bam_f_i) {
        Bam* bam_f = bam_fs[bam_f_i];
        size_t line = 1;
    try {
        // Parse the @RG header lines.
        const char* p = bam_f->h().text();
        string rg;
        string name;
        int stacks_id = -1;
        while (true) {
            if (strncmp(p, "@RG\t", 4) == 0) {
                rg.clear();
                name.clear();
                stacks_id = -1;
                // Parse the line.
                p += 4;
                while (*p && *p != '\n') {
                    const char* q = p;
                    while (*p && *p != '\n' && *p != '\t')
                        ++p;

                    if (strncmp(q, "ID:", 3) == 0) {
                        rg.assign(q+3, p);
                    } else if (strncmp(q, "SM:", 3) == 0) {
                        name.assign(q+3, p);
                    } else if (strncmp(q, "id:", 3) == 0) {
                        char* end;
                        stacks_id = strtol(q+3, &end, 10);
                        if (end != p)
                            throw exception();
                    }
                    if (*p == '\t')
                        ++p;
                }
                if (rg.empty()) {
                    cerr << "Error: @RG line has no ID (identifier) field.\n";
                    throw exception();
                } else if (name.empty()) {
                    cerr << "Error: @RG line has no SM (sample name) field.\n";
                    throw exception();
                }
                // Record the read group.
                readgroups_per_f[bam_f_i].push_back({move(rg), move(name), stacks_id});
            }
            p = strchr(p, '\n');
            if (p == NULL)
                break;
            ++p;
            ++line;
        }
    } catch (exception&) {
        cerr << "Error: Malformed BAM header.\n"
             << "Error: At header line " << line << " in file '" << bam_f->path << "'.\n";
        throw;
    }}

    // Check that all the files have read groups.
    bool no_rg = false;
    for (size_t bam_f_i=0; bam_f_i<bam_fs.size(); ++bam_f_i) {
        if (readgroups_per_f[bam_f_i].empty()) {
            no_rg = true;
            cerr << "Error: No read group(s) in BAM file '"
                 << bam_fs[bam_f_i]->path << "'.\n";
        }
    }
    if (no_rg)
        throw exception();

    // Check that duplicated readgroup IDs and sample names appear in a single
    // (ID,NAME) combination.
    map<string,pair<string,size_t>> name_to_rg;
    map<string,pair<string,size_t>> rg_to_name;
    map<string,int> name_to_stacks_id;
    for (size_t bam_f_i=0; bam_f_i<bam_fs.size(); ++bam_f_i) {
        for (const SampleData& d : readgroups_per_f[bam_f_i]) {
            auto rv1 = name_to_rg.insert({d.name, {d.rg, bam_f_i}});
            auto rv2 = rg_to_name.insert({d.rg, {d.name, bam_f_i}});
            if (!rv1.second && rv1.first->second.first != d.rg) {
                const pair<string,pair<string,size_t>>& n2r = *rv1.first;
                cerr << "Error: Incompatible read groups (ID:" << d.rg << ", SM:"
                     << d.name << ") and (ID:" << n2r.second.first << ", SM:"
                     << n2r.first << ").\n"
                     << "Error: In BAM files '" << bam_fs[bam_f_i]->path
                     << "' and '" << bam_fs[n2r.second.second]->path << ".\n";
                throw exception();
            }
            if (!rv2.second && rv2.first->second.first != d.name) {
                const pair<string,pair<string,size_t>>& r2n = *rv2.first;
                cerr << "Error: Incompatible read groups (ID:" << d.rg << ", SM:"
                     << d.name << ") and (ID:" << r2n.first << ", SM:"
                     << r2n.second.first << ").\n"
                     << "Error: In BAM files '" << bam_fs[bam_f_i]->path
                     << "' and '" << bam_fs[r2n.second.second]->path << ".\n";
                throw exception();
            }
            name_to_stacks_id[d.name] = d.stacks_id; // (Overwriting may happen.)
        }
    }

    // Initialize the MetaPopInfo.
    vector<string> names_all;
    for (auto& elem : name_to_rg)
        names_all.push_back(elem.first);
    mpopi_.init_names(names_all);

    // Set the sample IDs.
    for (auto& elem : name_to_stacks_id)
        mpopi_.set_sample_id(mpopi_.get_sample_index(elem.first), elem.second);

    // Initialize `read_groups_`, `rg_to_sample_`.
    for (size_t bam_f_i=0; bam_f_i<bam_fs.size(); ++bam_f_i) {
        for (SampleData& d : readgroups_per_f[bam_f_i]) {
            auto rv = read_groups_[bam_f_i].emplace(move(d.rg), mpopi_.get_sample_index(d.name));
            if (rv.second)
                rg_to_sample_[bam_f_i].emplace(rv.first->first.c_str(), rv.first);
        }
    }
}

inline
void BamPopInfo::reassign_ids() {
    for (size_t i=0; i<mpopi_.samples().size(); ++i)
        mpopi_.set_sample_id(i, i+1);
}

inline
size_t BamPopInfo::sample_of(const BamRecord& rec, size_t bam_f_i) {
    const char* rg = rec.read_group();
    if (rg == NULL) {
        cerr << "Error: Missing read group for BAM record '" << rec.qname() << "'.\n";
        throw exception();
    }
    size_t sample;
    try {
        sample = rg_to_sample_.at(bam_f_i).at(rg)->second;
    } catch(std::out_of_range&) {
        cerr << "Error: BAM record '" << rec.qname() << "' has read group '" << rg
             << "', which was not declared in the header.\n";
        throw exception();
    }
    return sample;
}

inline
BamCLocReader::BamCLocReader(vector<Bam*>&& bam_fs)
:
    bam_fs_(move(bam_fs)),
    bpopi_(bam_fs_),
    loc_i_(-1),
    eofs_(bam_fs_.size(), false),
    next_records_(bam_fs_.size())
{
    // Check that headers are consistent.
    size_t bam_f_i=1;
    if (bam_fs_.empty())
        DOES_NOT_HAPPEN;
    try {
        for (; bam_f_i<bam_fs_.size(); ++bam_f_i)
            BamHeader::check_same_ref_chroms(bam_fs_[0]->h(), bam_fs_[bam_f_i]->h());
    } catch (exception& e) {
        cerr << "Error: Inconsistent BAM headers; in files '" << bam_fs_[0]->path
             << "' and '" << bam_fs_[bam_f_i]->path << "'.\n";
        throw;
    }

    // Check that Stacks sample IDs are present.
    for (const Sample& s : mpopi().samples())
        if (s.id == -1)
            cerr << "WARNING: Sample '" << s.name << "' is missing its (integer) ID; was tsv2bam run correctly?!\n";

    //
    // Load genomic alignment information from the header.
    // If the first locus has alignment information: this is a ref-based analysis,
    // record all the alignment positions. Otherwise, assume it's a de novo analysis.
    //
    const char* p = bam_fs_[0]->h().text();
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
    if(!loci_aln_positions_.empty() && loci_aln_positions_.size() != bam_fs_[0]->h().n_ref_chroms())
        DOES_NOT_HAPPEN;
}

inline
bool BamCLocReader::read_one_locus(CLocReadSet& readset) {
    assert(&readset.mpopi() == &mpopi()); // Otherwise sample indexes may be misleading.

    ++loc_i_;
    if (loc_i_ == 0) {
        // Read the first records.
        for (size_t bam_f_i=0; bam_f_i<bam_fs_.size(); ++bam_f_i) {
            Bam* bam_f = bam_fs_[bam_f_i];
            BamRecord& rec = next_records_[bam_f_i];
            uchar& eof = eofs_[bam_f_i];
            if ((eof = !bam_f->next_record(rec, true))) {
                cerr << "Error: No records in BAM file '" << bam_f->path << "'.\n";
                throw exception();
            } else if (!rec.is_primary()) {
                cerr << "Error: BAM file '" << bam_f->path
                     << "' unexpectedly contains non-primary records.\n";
                throw exception();
            }
        }
    } else if (loc_i_ == int32_t(n_loci())) {
        assert(std::accumulate(eofs_.begin(), eofs_.end(), size_t(0)) == eofs_.size());
        return false;
    }

    readset.clear();
    readset.bam_i(loc_i_);
    readset.id(atoi(bam_fs_[0]->h().chrom_str(loc_i_)));
    if (!loci_aln_positions_.empty()) {
        const char* p = loci_aln_positions_.at(loc_i_);
        const char* q = p;
        while (*q && *q != '\n' && *q != '\t')
            ++q;
        readset.pos(PhyLoc(string(p, q)));
    }

    for (size_t bam_f_i=0; bam_f_i<bam_fs_.size(); ++bam_f_i) {
        Bam* bam_f = bam_fs_[bam_f_i];
        BamRecord& rec = next_records_[bam_f_i];
        uchar& eof = eofs_[bam_f_i];

        // Read all the reads of the locus, and one more.
        while (!eof && rec.chrom() == loc_i_) {
            if (rec.is_read2())
                readset.add_pe(SRead(Read(rec.seq(), string(rec.qname())+"/2"), bpopi_.sample_of(rec, bam_f_i)));
            else if (rec.is_read1())
                readset.add(SRead(Read(rec.seq(), string(rec.qname())+"/1"), bpopi_.sample_of(rec, bam_f_i)));
            else
                // If tsv2bam wasn't given paired-end reads, no flag was set and the
                // read names were left unchanged, so we also don't touch them.
                readset.add(SRead(Read(rec.seq(), string(rec.qname())), bpopi_.sample_of(rec, bam_f_i)));

            eof = !bam_f->next_record(rec, true);
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
        vector<Bam*>&& bam_fs,
        const Config& cfg
) :
    cfg_(cfg),
    stats_(),
    bam_fs_(move(bam_fs)),
    bpopi_(bam_fs_),
    eofs_(bam_fs_.size(), false),
    next_records_(bam_fs_.size()),
    treat_next_records_as_fw_(bam_fs_.size())
{
    if (!cfg_.paired)
        cfg_.max_insert_refsize = 0;
    bpopi_.reassign_ids();

    size_t bam_f_i=1;
    if (bam_fs_.empty())
        DOES_NOT_HAPPEN;
    try {
        for (; bam_f_i<bam_fs_.size(); ++bam_f_i)
            BamHeader::check_same_ref_chroms(bam_fs_[0]->h(), bam_fs_[bam_f_i]->h());
    } catch (exception& e) {
        cerr << "Error: Inconsistent BAM headers; in files '" << bam_fs_[0]->path
             << "' and '" << bam_fs_[bam_f_i]->path << "'.\n";
        throw;
    }
}

inline
bool BamCLocBuilder::next_record(size_t bam_f_i)
{
    Bam* bam_f = bam_fs_[bam_f_i];
    BamRecord& r = next_records_[bam_f_i];
    uchar& treat_as_fw = treat_next_records_as_fw_[bam_f_i];

    while (true) {
        if(!bam_f->next_record(r, true))
            return false;

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
            treat_as_fw = true;
        } else if (r.is_read1()) {
            treat_as_fw = true;
        } else if (r.is_read2()) {
            treat_as_fw = false;
        } else {
            cerr << "Error: --paired is on but BAM record '" << r.qname()
                 << "' is not flagged as READ1 nor READ2.\n";
            throw exception();
        }

        // Check that there's a sequence.
        if (r.hts_l_seq() <= 0) {
            cerr << "Error: No sequence at BAM record '" << r.qname() << "'.\n";
            throw exception();
        }

        // Check the CIGAR.
        if (r.hts_n_cigar() == 0) {
            cerr << "Error: Empty CIGAR at BAM record '" << r.qname() << "'.\n";
            throw exception();
        }
        if (treat_as_fw) {
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
        break;
    }
    return true;
}

inline
void BamCLocBuilder::add_next_record_to_the_window(size_t bam_f_i)
{
    //
    // See `fill_window()`.
    // This inserts `next_record_` into `fw_reads_by_5prime_pos_` (if it is a
    // forward read or treated as such) or `pe_reads_by_5prime_pos_` (if it is
    // a paired-end read).
    //

    BamRecord& rec = next_records_[bam_f_i];

    // Determine the 5prime position.
    Pos5 rec_5prime_pos;
    rec_5prime_pos.chrom = rec.chrom();
    if (rec.is_rev_compl()) {
        rec_5prime_pos.strand = false;
        rec_5prime_pos.bp = rec.pos() + BamRecord::cig_ref_len(rec) - 1;
    } else {
        rec_5prime_pos.strand = true;
        rec_5prime_pos.bp = rec.pos();
    }

    // Determine the sample.
    size_t sample = bpopi_.sample_of(rec, bam_f_i);

    // Save the record.
    if (treat_next_records_as_fw_[bam_f_i]) {
        auto locus_itr = fw_reads_by_5prime_pos_.insert(
                std::make_pair(move(rec_5prime_pos), vector<pair<BamRecord,SampleIdx>>())
                ).first;
        locus_itr->second.push_back({move(rec), sample});
    } else {
        auto pe_read_itr = pe_reads_by_5prime_pos_.insert(
                std::make_pair(move(rec_5prime_pos), std::make_pair(move(rec), move(sample)))
                );
        // Also add an entry in the {pe_read_name : pe_record} map.
        pe_reads_by_name_.insert( {pe_read_itr->second.first.qname(), pe_read_itr} );
    }
}

inline
bool BamCLocBuilder::fill_window()
{
    /*
     * This defines the behavior of the sliding window.
     *
     * This guarantees that any usable alignment of which the leftmost base
     * is within `max_insert_size` bases of the first base of the leftmost
     * cutsite is loaded.
     *
     * `fw_reads_by_5prime_pos_` is a map of {[first cutsite base pos] : BamRecord}
     * Cutsites can be on the left or right of the read alignment depending on the
     * strand, but this is easy to handle since `BamCLocBuilder::next_record()`
     * ensures that the BAM record's TLEN field is properly set.
     *
     * `pe_reads_by_5prime_pos_` also uses the 5prime position (i.e. that has been
     * tlen-corrected if the strand is minus) to ease insert lengths computations.
     * Paired-end reads need to be sorted by position so that we can discard those
     * that were not matched and have become too far behind to find one at the given
     * `max_insert_size`. Paired-end records can also be accessed by name (for
     * mate matching) via `pe_reads_by_name_`.
     */

    for (size_t bam_f_i=0; bam_f_i<bam_fs_.size(); ++bam_f_i) {
    try {
        uchar& eof = eofs_[bam_f_i];
        BamRecord& rec = next_records_[bam_f_i];
        if (!eof) {
            if (rec.empty()) {
                if ((eof = !next_record(bam_f_i))) {
                    if (bam_fs_[bam_f_i]->n_records_read() == 0)
                        cerr << "Error: No BAM records.\n";
                    else
                        cerr << "Error: All records were discarded.\n";
                    throw exception();
                }
                assert(!rec.empty());
            }
            while (!eof && fw_reads_by_5prime_pos_.empty()) {
                add_next_record_to_the_window(bam_f_i); // (Note: `next_record_` may be a paired-end read.)
                eof = !next_record(bam_f_i);
            }
            while (!eof
                    && rec.chrom() == fw_reads_by_5prime_pos_.begin()->first.chrom
                    && size_t(rec.pos()) <= fw_reads_by_5prime_pos_.begin()->first.bp + cfg_.max_insert_refsize
            ) {
                add_next_record_to_the_window(bam_f_i);
                eof = !next_record(bam_f_i);
            }
        }
    } catch (exception& e) {
        cerr << "Error: (At the " << bam_fs_[bam_f_i]->n_records_read()
             << "th record in file '" << bam_fs_[bam_f_i]->path << "'.)\n";
        throw;
    }}

    if (fw_reads_by_5prime_pos_.empty()) {
        assert(std::accumulate(eofs_.begin(), eofs_.end(), size_t(0)) == eofs_.size());
        return false;
    }

    if (cfg_.paired) {
        // Clean up paired-end reads that are too far behind.
        Pos5 leftmost_cutsite = fw_reads_by_5prime_pos_.begin()->first;
        while (!pe_reads_by_5prime_pos_.empty()
                && (
                    pe_reads_by_5prime_pos_.begin()->first.chrom != leftmost_cutsite.chrom
                    || pe_reads_by_5prime_pos_.begin()->first.bp + cfg_.max_insert_refsize < size_t(leftmost_cutsite.bp) + 1
                )
        ) {
            pe_reads_by_name_.erase(pe_reads_by_5prime_pos_.begin()->second.first.qname());
            pe_reads_by_5prime_pos_.erase(pe_reads_by_5prime_pos_.begin());
        }
        // As we enforce sorting there shouldn't be reads from further-ranked chromosomes.
        assert(pe_reads_by_5prime_pos_.empty() || pe_reads_by_5prime_pos_.rbegin()->first.chrom == leftmost_cutsite.chrom);
    } else {
        assert(pe_reads_by_5prime_pos_.empty());
    }

    return true;
}

inline
bool BamCLocBuilder::build_one_locus(CLocAlnSet& aln_loc)
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

        vector<pair<BamRecord,SampleIdx>>& loc_records = fw_reads_by_5prime_pos_.begin()->second;

        // Get the samples of the records, and tally sample occurrences.
        vector<size_t> n_reads_per_sample (mpopi().samples().size(), 0);
        for (const pair<BamRecord,SampleIdx>& r : loc_records)
            ++n_reads_per_sample[r.second];

        // Remove reads from samples that have less than `cfg_.min_reads_per_sample` reads.
        for (pair<BamRecord,SampleIdx>& r : loc_records)
            if (n_reads_per_sample[r.second] < cfg_.min_reads_per_sample)
                // Bad sample, mark the record for deletion.
                r.first.destroy();
        loc_records.erase(std::remove_if(
                loc_records.begin(), loc_records.end(),
                [] (const pair<BamRecord,SampleIdx>& r) { return r.first.empty(); }
                ), loc_records.end());

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
    PhyLoc loc_pos (bam_fs_.at(0)->h().chrom_str(cutsite.chrom),
                    cutsite.bp,
                    cutsite.strand ? strand_plus : strand_minus);
    size_t max_ref_span = 0;

    // Forward reads.
    const vector<pair<BamRecord,SampleIdx>>& fw_records = fw_reads_by_5prime_pos_.begin()->second;
    vector<SAlnRead> fw_reads;
    for (const pair<BamRecord,SampleIdx>& pair : fw_records) {
        const BamRecord& r = pair.first;
        size_t sample = pair.second;

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
        for (auto& fw_pair : fw_records) {
            const BamRecord& fw_rec = fw_pair.first;

            // Find a mate.
            auto pe_name_itr = pe_reads_by_name_.find(fw_rec.qname());
            if (pe_name_itr == pe_reads_by_name_.end())
                continue;
            std::multimap<Pos5,pair<BamRecord,SampleIdx>>::iterator pe_rec_itr = pe_name_itr->second;

            Pos5 pe_pos5 = pe_rec_itr->first;
            const BamRecord& pe_rec = pe_rec_itr->second.first;
            size_t pe_sample = pe_rec_itr->second.second;

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
            if (pe_sample != fw_pair.second) {
                cerr << "Warning: Paired reads '" << fw_rec.qname() << "' and '"
                     << pe_rec.qname() << "' belong to different samples ('"
                     << mpopi().samples()[fw_pair.second].name << "' and '"
                     << mpopi().samples()[pe_sample].name << "').\n";
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
    aln_loc.reinit(++stats_.n_loci_built, loc_pos, &mpopi());
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
