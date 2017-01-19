#include <cstring>

#include "BamI.h"

using namespace std;

void write_bam_header(
        htsFile* bam_f,
        const string& header_text,
        const vector<pair<string, uint32_t> >& chrs
        ) {

    bam_hdr_t* hdr = bam_hdr_init();

    // Unstructured text part.
    hdr->l_text = header_text.length()+1; // null-terminated
    hdr->text = new char[hdr->l_text];
    strcpy(hdr->text, header_text.c_str());

    // Reference sequences (catalog loci).
    hdr->n_targets = chrs.size();
    hdr->target_len = new uint32_t[hdr->n_targets];
    hdr->target_name = new char*[hdr->n_targets];
    size_t i=0;
    for (auto& c : chrs) {
        hdr->target_len[i] = c.second;
        hdr->target_name[i] = new char[c.first.length()+1];
        strcpy(hdr->target_name[i], c.first.c_str());
        i++;
    }

    // Write the header.
    int r = bam_hdr_write(bam_f->fp.bgzf, hdr);
    if (r != 0) {
        cerr << "Error: Writing of BAM header failed (`bam_hdr_write()`returned " << r << ").\n";
        throw exception();
    }

    // Clean up.
    bam_hdr_destroy(hdr);

}
