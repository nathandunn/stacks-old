#include <cstring>

#include "constants.h"
#include "BamI.h"

using namespace std;

void write_bam_header(
        htsFile* bam_f,
        const string& header
        ) {

    bam_hdr_t* hdr = sam_hdr_parse(header.length()+1, header.c_str());
    hdr->l_text = header.length()+1; // null-terminated
    hdr->text = new char[hdr->l_text];
    strcpy(hdr->text, header.c_str());

    // Write the header.
    int r = bam_hdr_write(bam_f->fp.bgzf, hdr);
    if (r != 0) {
        cerr << "Error: Writing of BAM header failed (`bam_hdr_write()`returned " << r << ").\n";
        throw exception();
    }

    // Clean up.
    bam_hdr_destroy(hdr);
}
