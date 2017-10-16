#include <cstring>
#include <string>
#include <map>

#include "constants.h"
#include "BamI.h"

using namespace std;

const uint32_t o = UINT32_MAX;
const uint32_t cigar_c2i[128] = {
    o,o,o,o, o,o,o,o, o,o,o,o, o,o,o,o,
    o,o,o,o, o,o,o,o, o,o,o,o, o,o,o,o,
    o,o,o,o, o,o,o,o, o,o,o,o, o,o,o,o,
    o,o,o,o, o,o,o,o, o,o,o,o, o,7,o,o,

    o,o,9,o, 2,o,o,o, 5,1,o,o, o,0,3,o,
    6,o,o,4, o,o,o,o, 8,o,o,o, o,o,o,o,
    o,o,o,o, o,o,o,o, o,o,o,o, o,o,o,o,
    o,o,o,o, o,o,o,o, o,o,o,o, o,o,o,o,
};

void check_open(const htsFile* bam_f, const string& path) {
    if (bam_f == NULL) {
        cerr << "Error: Failed to open BAM file '" << path << "'.\n";
        throw exception();
    }
}

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
    r_->core.tid = chr_index;
    r_->core.pos = aln_pos;
    r_->core.bin = 0; // No idea
    r_->core.qual = 255;
    r_->core.l_qname = name.length() + 1; // `l_qname` includes the trailing \0.
    r_->core.flag = flg;
    r_->core.n_cigar = cig.size();
    r_->core.l_qseq = seq.length();
    r_->core.mtid = -1;
    r_->core.mpos = -1;
    r_->core.isize = -1; // No idea

    // bam1_t::data
    // Htslib says: "bam1_t::data -- all variable-length data, concatenated;
    // structure: qname-cigar-seq-qual-aux, concatenated".

    // Prepare the `aux` data.
    string rg = string() + "RG" + "Z" + to_string(read_group);

    // Determine the length of `data`.
    size_t l_aux = rg.length() + 1;
    r_->l_data = r_->core.l_qname + r_->core.n_cigar*sizeof(uint32_t) + seq.nbytes() + seq.length() + l_aux;
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
    p += r_->core.l_qname;
    //cigar
    for (const pair<char, uint>& op : cig) {
        // Cigars are uint32_t's with the length on the 28 high bits & op on the low 4 bits.
        *(uint32_t*)p = (uint32_t(op.second) <<BAM_CIGAR_SHIFT) | cigar_c2i[size_t(op.first)];
        p += sizeof(uint32_t);
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
    uint32_t* cigar = (uint32_t*)(r_->data + r_->core.l_qname);
    r_->core.bin = hts_reg2bin(r_->core.pos, r_->core.pos + bam_cigar2rlen(r_->core.n_cigar, cigar), 14, 5);
}

Bam::Bam(const char *path) : Input(), bam_fh(NULL), n_records_read_(0), hdr() {
    this->path   = string(path);
    bam_fh = hts_open(path, "r");
    if (bam_fh == NULL) {
        cerr << "Error: Failed to open BAM file '" << path << "'.\n";
        throw exception();
    } else if (bam_fh->format.format != bam) {
        cerr << "Error: '" << path << "':";
        if (bam_fh->format.format == sam)
            cerr << " this is a SAM file (and BAM was specified).\n";
        else
            cerr << " not a BAM file.\n";
        throw exception();
    }
    hdr.reinit(bam_fh);
};
