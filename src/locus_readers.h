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

public:
    BamCLocReader(const string& bam_path);
    ~BamCLocReader() {if (bam_f_) delete bam_f_;}

    const MetaPopInfo& mpopi() const {return mpopi_;}

    // Reads one locus. Returns false on EOF.
    bool read_one_locus(CLocReadSet& readset);
};

class VcfCLocReader {
    unique_ptr<VcfAbstractParser> vcf_f_;
    VcfRecord next_rec_;
    bool eof_;
public:
    VcfCLocReader(const string& vcf_path);

    const VcfHeader& header() const {return vcf_f_->header();}
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
          rg_to_sample_()
        {

    bam_f_ = new Bam(bam_path.c_str());

    // Get the list of samples from the header.
    BamHeader::ReadGroups read_groups = bam_f_->h().read_groups();

    // Create the MetaPopInfo object.
    vector<string> samples;
    for (auto& rg : read_groups)
        samples.push_back(rg.second.at("SM"));
    mpopi_.init_names(samples);

    // Parse sample IDs, if any.
    for (auto& rg : read_groups) {
        auto id = rg.second.find("id");
        if (id != rg.second.end())
            mpopi_.set_sample_id(mpopi_.get_sample_index(rg.second.at("SM")), stoi(id->second));
    }

    // Fill the (read group : sample) map.
    for (auto& rg : read_groups)
        rg_to_sample_.insert({rg.first, mpopi_.get_sample_index(rg.second.at("SM"))});

    // Read the very first record.
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

    readset.clear();
    if (bam_f_ == NULL)
        // EOF was hit at the end of the previous locus.
        return false;

    const BamRecord& rec = bam_f_->r(); //the current record

    // Parse the locus ID.
    int32_t curr_chrom = rec.chrom();
    static int32_t last_chrom = -1;
    if (curr_chrom < last_chrom) {
        cerr << "Error: BAM file isn't properly sorted.\n";
        throw exception();
    }
    last_chrom = curr_chrom;
    readset.id(atoi(bam_f_->h().chrom_str(curr_chrom)));

    // Read all the reads of the locus, and one more.
    do {
        if (rec.is_read2())
            readset.add_pe(SRead(Read(rec.seq(), rec.qname()+"/2"), rg_to_sample_.at(rec.read_group())));
        else if (rec.is_read1())
            readset.add(SRead(Read(rec.seq(), rec.qname()+"/1"), rg_to_sample_.at(rec.read_group())));
        else
            // If tsv2bam wasn't given paired-end reads, no flag was set and the
            // read names were left unchanged, so we also don't touch them.
            readset.add(SRead(Read(rec.seq(), rec.qname()), rg_to_sample_.at(rec.read_group())));

        if(!bam_f_->next_record()) {
            // EOF
            delete bam_f_;
            bam_f_ = NULL;
            break;
        }
    } while (rec.chrom() == curr_chrom);

    return true;
}

inline
VcfCLocReader::VcfCLocReader(const string& vcf_path)
        : vcf_f_(Vcf::adaptive_open(vcf_path)),
          next_rec_(),
          eof_(false)
{
    vcf_f_->read_header();

    // Read the very first record.
    if(!vcf_f_->next_record(next_rec_))
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
    string curr_chrom = next_rec_.chrom;
    records.push_back(move(next_rec_));

    // Read all the records of the locus, and one more.
    while (records.back().chrom == curr_chrom) {
        records.push_back(VcfRecord());
        if (!vcf_f_->next_record(records.back())) {
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
