#ifndef BAMCLOCREADER_H
#define BAMCLOCREADER_H

#include <vector>
#include <string>
#include <map>

#include "constants.h"
#include "MetaPopInfo.h"
#include "BamI.h"
#include "locus.h"

class BamCLocReader {
    Bam* bam_f_;
    const MetaPopInfo& mpopi_;
    map<string, size_t> rg_to_sample_;

public:
    BamCLocReader(const string& bam_path, MetaPopInfo& mpopi);
    ~BamCLocReader() {if (bam_f_) delete bam_f_;}

    // Reads one locus. Returns false on EOF.
    bool read_one_locus(CLocReadSet& readset);
};

inline
BamCLocReader::BamCLocReader(const string& bam_path, MetaPopInfo& mpopi)
        : bam_f_(NULL),
          mpopi_(mpopi),
          rg_to_sample_()
        {

    bam_f_ = new Bam(bam_path.c_str());

    // Get the list of samples from the header.
    BamHeader::ReadGroups read_groups = bam_f_->h().read_groups();

    // Create the MetaPopInfo object.
    vector<string> samples;
    for (auto& rg : read_groups)
        samples.push_back(rg.second.at("SM"));
    mpopi.init_names(samples);

    // Fill the (read group : sample) map.
    for (auto& rg : read_groups)
        rg_to_sample_.insert({rg.first, mpopi.get_sample_index(rg.second.at("SM"))});

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

#endif
