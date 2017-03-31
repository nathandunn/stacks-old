#include <cstring>
#include <string>
#include <map>

#include "constants.h"
#include "BamI.h"

using namespace std;

const size_t o = -1;
const size_t cigar_c2i[128] = {
    o,o,o,o, o,o,o,o, o,o,o,o, o,o,o,o,
    o,o,o,o, o,o,o,o, o,o,o,o, o,o,o,o,
    o,o,o,o, o,o,o,o, o,o,o,o, o,o,o,o,
    o,o,o,o, o,o,o,o, o,o,o,o, o,7,o,o,

    o,o,9,o, 2,o,o,o, 5,1,o,o, o,0,3,o,
    6,o,o,4, o,o,o,o, 8,o,o,o, o,o,o,o,
    o,o,o,o, o,o,o,o, o,o,o,o, o,o,o,o,
    o,o,o,o, o,o,o,o, o,o,o,o, o,o,o,o,
};

void write_bam_header(
        htsFile* bam_f,
        const string& header
        ) {

    bam_hdr_t* hdr = sam_hdr_parse(header.length()+1, header.c_str());
    hdr->l_text = header.length()+1; // null-terminated
    hdr->text = (char*) malloc(hdr->l_text);
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

void BamRecord::assign(
        const string& name,
        uint16_t flg,
        int32_t chr_index,
        int32_t aln_pos,
        const vector<pair<char,uint>>& cig,
        const DNASeq4& seq,
        size_t read_group
        ) {

    // bam1_t::core
    c_.tid = chr_index;
    c_.pos = aln_pos;
    c_.bin = 0; // No idea
    c_.qual = 255;
    c_.l_qname = name.length() + 1; // `l_qname` includes the trailing \0.
    c_.flag = flg;
    c_.n_cigar = cig.size();
    c_.l_qseq = seq.length();
    c_.mtid = -1;
    c_.mpos = -1;
    c_.isize = -1; // No idea

    // bam1_t::data
    // Htslib says: "bam1_t::data -- all variable-length data, concatenated;
    // structure: qname-cigar-seq-qual-aux, concatenated".

    // Prepare the `aux` data.
    string rg = string() + "RG" + "Z" + to_string(read_group);

    // Determine the length of `data`.
    size_t l_aux = rg.length() + 1;
    r_->l_data = c_.l_qname + c_.n_cigar*sizeof(uint32_t) + seq.nbytes() + seq.length() + l_aux;
    if (r_->l_data > r_->m_data) {
        if (r_->data != NULL)
            free(r_->data);
        r_->m_data = r_->l_data;
        r_->data = (uchar*) malloc(r_->m_data);
    }

    // Fill the data array.
    uchar* p = r_->data;
    //qname
    strcpy((char*)p, name.c_str());
    p += c_.l_qname;
    //cigar
    for (const pair<char, uint>& op : cig) {
        // Cigars are uint32_t's with the length on the 32 high bits & op on the low 4 bits.
        *(uint32_t*)p = uint32_t(op.second) <<BAM_CIGAR_SHIFT;
        *(uint32_t*)p |= cigar_c2i[size_t(op.first)];
        p += c_.n_cigar*sizeof(uint32_t);
    }
    //seq & qual
    memcpy(p, seq.vdata(), seq.nbytes());
    p += seq.nbytes();
    memset(p, 0xFF, seq.length());
    p += seq.length();
    //aux
    memcpy(p, rg.c_str(), rg.length()+1);

    // bam1_t::core.bin
    // c.f. `sam_parse1()`; I have no idea what this is.
    uint32_t* cigar = (uint32_t*)(r_->data + c_.l_qname);
    c_.bin = hts_reg2bin(c_.pos, c_.pos + bam_cigar2rlen(c_.n_cigar, cigar), 14, 5);
}

void BamRecord::write_to(htsFile* bam_f) const {
    if (bam_write1(bam_f->fp.bgzf, r_) < 0) {
        cerr << "Error: Writing of BAM record failed.\n";
        throw exception();
    }
}
