#include <cstring>
#include <string>
#include <map>

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

BamHeader::ReadGroups BamHeader::read_groups() const {
    ReadGroups rgs;

    // xxx (Feb2017) I currently store the sample ID in the read group ID, but
    // BAM read group IDs are unstable. It would be better to use a custom field
    // such as 'si' (sample id).

    const char* end = h_-> text + h_->l_text;
    const char* p = h_->text;
    while (*p) {
        // Process a line.

        if (p + 3 < end && strncmp(p, "@RG\t", 4) == 0) {
            // This is a @RG line.
            map<string, string> rg;

            p+=4;
            while (*p && *p != '\n') {
                // Process a tag.

                string tag;
                while (*p && *p != '\n' && *p != '\t' && *p != ':') {
                    tag += *p;
                    ++p;
                }
                if (*p == ':')
                    ++p;

                string value;
                while (*p && *p != '\n' && *p != '\t') {
                    value += *p;
                    ++p;
                }

                if (tag.size() == 2 && !value.empty())
                    rg.insert({tag, value});

                if (*p == '\t')
                    ++p;
            }

            try {
                string id = rg.at("ID"); // might throw
                rgs.insert({move(id), move(rg)});
            } catch (out_of_range&) {} // ignore
        }

        while(*p && *p != '\n')
            ++p;

        if (*p == '\n')
            ++p;
    }

    return rgs;
}
