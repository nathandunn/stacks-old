#ifndef LOCUS_READERS_H
#define LOCUS_READERS_H

#include "constants.h"
#include "MetaPopInfo.h"
#include "BamI.h"
#include "locus.h"
#include "Vcf.h"

class BamCLocReader {
    Bam* bam_f_;
    MetaPopInfo mpopi_;
    map<string, size_t> rg_to_sample_;
    vector<const char*> loci_aln_positions_;

    int32_t loc_i_; // Index of the next locus (chromosome). Incremented by read_one_locus().

public:
    BamCLocReader(const string& bam_path);
    ~BamCLocReader() {if (bam_f_) delete bam_f_;}

    size_t n_loci() const {return bam_f_->h().n_ref_chroms();}
    const MetaPopInfo& mpopi() const {return mpopi_;}
    int target2id(size_t target_i) const {return atoi(bam_f_->h().chrom_str(target_i));}

    // Reads one locus. Returns false on EOF.
    bool read_one_locus(CLocReadSet& readset);
};

class VcfCLocReader {
    VcfParser vcf_f_;
    VcfRecord next_rec_;
    bool eof_;
public:
    VcfCLocReader(const string& vcf_path);

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
BamCLocReader::BamCLocReader(const string& bam_path)
        : bam_f_(NULL),
          mpopi_(),
          rg_to_sample_(),
          loc_i_(-1)
        {

    bam_f_ = new Bam(bam_path.c_str());

    //
    // Create the MetaPopInfo object from the header.
    //
    {
        vector<string> rg_ids;
        vector<string> samples;
        vector<int> sample_ids;

        // Parse the @RG header lines.
        const char* p = bam_f_->h().text();
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
        mpopi_.init_names(samples);

        // Set the sample IDs.
        assert(sample_ids.size() == samples.size());
        for (size_t s=0; s<samples.size(); ++s)
            mpopi_.set_sample_id(mpopi_.get_sample_index(samples[s]),sample_ids[s]);

        // Initialize `rg_to_sample_`.
        assert(rg_ids.size() == samples.size());
        for (size_t s=0; s<mpopi_.samples().size(); ++s)
            rg_to_sample_[rg_ids[s]] = mpopi_.get_sample_index(samples[s]);
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
        cerr << "Error: No records in BAM file '" << bam_path << "'.\n";
        throw exception();
    } else if (bam_f_->r().is_unmapped()) {
        cerr << "Error: BAM file '" << bam_path << "' unexpectedly contains unmapped records.\n";
        throw exception();
    }
}

inline
bool BamCLocReader::read_one_locus(CLocReadSet& readset) {
    assert(&readset.mpopi() == &mpopi_); // Otherwise sample indexes may be misleading.

    ++loc_i_;
    if (size_t(loc_i_) == n_loci()) {
        assert(bam_f_->eof());
        return false;
    }

    readset.clear();
    readset.id(atoi(bam_f_->h().chrom_str(loc_i_)));
    if (!loci_aln_positions_.empty()) {
        const char* p = loci_aln_positions_.at(loc_i_);
        const char* q = p;
        while (*q && *q != '\n' && *q != '\t')
            ++q;
        readset.pos(PhyLoc(string(p, q)));
    }

    if (!bam_f_->eof()) {
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

            if (!bam_f_->next_record())
                break;
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
VcfCLocReader::VcfCLocReader(const string& vcf_path)
        : vcf_f_(vcf_path),
          next_rec_(),
          eof_(false)
{
    // Read the very first record.
    if(!vcf_f_.next_record(next_rec_))
        eof_ = true;
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
