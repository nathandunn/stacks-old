// -*-mode:c++; c-style:k&r; c-basic-offset:4;-*-
//
// Copyright 2013-2016, Julian Catchen <jcatchen@illinois.edu>
//
// This file is part of Stacks.
//
// Stacks is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// Stacks is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with Stacks.  If not, see <http://www.gnu.org/licenses/>.
//

#ifndef __BAMI_H__
#define __BAMI_H__
//#ifdef HAVE_BAM

//
// Code to parse binary BAM format. This format is created for
// reads that have been aligned to a reference genome.
//

#include <string>
#include <vector>
#include <utility>

#include "input.h"
#include "sam.h"

// Write a header to a BAM file.
// The header text should include all tags except @SQ, and seems to be written
// in plain text in the BAM file.
void write_bam_header(htsFile* bam_f, const std::string& header_text, const std::vector<std::pair<std::string, uint32_t> >& chrs);

class Bam: public Input {
    htsFile   *bam_fh;
    bam_hdr_t *bamh;
    bam1_t    *aln;

    map<uint, string> chrs;

    int parse_header();
    int parse_cigar(const char *, vector<pair<char, uint> > &, bool);
    int parse_bam_cigar(vector<pair<char, uint> > &, bool);
    int find_start_bp_pos(int, vector<pair<char, uint> > &);
    int find_start_bp_neg(int, vector<pair<char, uint> > &);
    int edit_gaps(vector<pair<char, uint> > &, char *);

 public:
    Bam(const char *path) : Input() {
        this->path   = string(path);
        this->bam_fh = hts_open(path, "r");
        this->aln    = bam_init1();

        this->parse_header();
    };
    ~Bam() {
        hts_close(this->bam_fh);
        bam_hdr_destroy(this->bamh);
        bam_destroy1(this->aln);
    };
    Seq *next_seq();
    int  next_seq(Seq&);
};

int
Bam::parse_header()
{
    this->bamh = bam_hdr_init();
    this->bamh = sam_hdr_read(this->bam_fh);

    for (uint j = 0; j < (uint) this->bamh->n_targets; j++) {
        //
        // Record the mapping from integer ID to chromosome name that we will see in BAM records.
        //
        this->chrs[j] = string(this->bamh->target_name[j]);
    }

    return 0;
}

Seq *
Bam::next_seq()
{
    Seq* s = new Seq();
    if(next_seq(*s) != 1) {
        delete s;
        s = NULL;
    }
    return s;
}

int
Bam::next_seq(Seq& s)
{
    int bytes_read = 0;
    int sflag      = 0;
    int flag       = 0;

    //
    // Read a record from the file, skipping unmapped reads,  and place it in a Seq object.
    //
    do {
        bytes_read = sam_read1(this->bam_fh, this->bamh, this->aln);

        if (bytes_read <= 0)
            return 0;

        flag = ((this->aln->core.flag & BAM_FUNMAP) != 0);

    } while (flag == 1);

    //
    // Check which strand this is aligned to:
    //   SAM reference: FLAG bit 0x10 - sequence is reverse complemented
    //
    sflag = ((this->aln->core.flag & BAM_FREVERSE) != 0);

    //
    // If the read was aligned on the reverse strand (and is therefore reverse complemented)
    // alter the start point of the alignment to reflect the right-side of the read, at the
    // end of the RAD cut site.
    //
    // To accomplish this, we must parse the alignment CIGAR string
    //
    vector<pair<char, uint> > cigar;
    this->parse_bam_cigar(cigar, sflag);

    uint bp = sflag ?
        this->find_start_bp_neg(this->aln->core.pos, cigar) :
        this->find_start_bp_pos(this->aln->core.pos, cigar);

    //
    // Check if this is the primary or secondary alignment.
    //
    alnt aln_type = pri_aln;
    flag = ((this->aln->core.flag & BAM_FSECONDARY) != 0);
    if (flag)
	aln_type = sec_aln;
    //
    // Check if this is a supplemenatry (chimeric) alignment (not yet defined in Bam.h).
    //
    flag = ((this->aln->core.flag & 2048) != 0);
    if (flag)
	aln_type = sup_aln;

    //
    // Fetch the sequence.
    //
    string  seq;
    uint8_t j;

    seq.reserve(this->aln->core.l_qseq);

    for (int i = 0; i < this->aln->core.l_qseq; i++) {
        j = bam_seqi(bam_get_seq(this->aln), i);
        switch(j) {
        case 1:
            seq += 'A';
            break;
        case 2:
            seq += 'C';
            break;
        case 4:
            seq += 'G';
            break;
        case 8:
            seq += 'T';
            break;
        case 15:
            seq += 'N';
            break;
        }
    }

    //
    // Fetch the quality score.
    //
    string   qual;
    uint8_t *q = bam_get_qual(this->aln);
    for (int i = 0; i < this->aln->core.l_qseq; i++) {
        qual += char(int(q[i]) + 33);
    }

    string chr = this->chrs[this->aln->core.tid];

    //
    // Calculate the percentage of the sequence that was aligned to the reference.
    //
    double len = 0.0;
    for (uint i = 0; i < cigar.size(); i++)
        switch (cigar[i].first) {
        case 'M':
        case 'I':
        case '=':
            len += cigar[i].second;
        }
    double pct_aln = len / double(seq.length());

    s = Seq((const char *) bam_get_qname(this->aln), seq.c_str(), qual.c_str(),
            chr.c_str(), bp, sflag ? strand_minus : strand_plus, 
            aln_type, pct_aln);

    if (cigar.size() > 0)
        this->edit_gaps(cigar, s.seq);

    return 1;
}

int
Bam::parse_bam_cigar(vector<pair<char, uint> > &cigar, bool orientation)
{
    int  op, len;
    char c;
    uint32_t *cgr = bam_get_cigar(this->aln);

    for (int k = 0; k < this->aln->core.n_cigar; k++) {
        op  = cgr[k] &  BAM_CIGAR_MASK;
        len = cgr[k] >> BAM_CIGAR_SHIFT;

        switch(op) {
        case BAM_CMATCH:
            c = 'M';
            break;
        case BAM_CINS:
            c = 'I';
            break;
        case BAM_CDEL:
            c = 'D';
            break;
        case BAM_CREF_SKIP:
            c = 'N';
            break;
        case BAM_CSOFT_CLIP:
            c = 'S';
            break;
        case BAM_CHARD_CLIP:
            c = 'H';
            break;
        case BAM_CPAD:
            c = 'P';
            break;
        }

        //
        // If aligned to the negative strand, sequence has been reverse complemented and
        // CIGAR string should be interpreted in reverse.
        //
        if (orientation == strand_plus)
            cigar.push_back(make_pair(c, len));
        else
            cigar.insert(cigar.begin(), make_pair(c, len));
    }

    return 0;
}

int
Bam::parse_cigar(const char *cigar_str, vector<pair<char, uint> > &cigar, bool orientation)
{
    char buf[id_len];
    int  dist;
    const char *p, *q;

    p = cigar_str;

    if (*p == '*') return 0;

    while (*p != '\0') {
        q = p + 1;

        while (*q != '\0' && isdigit(*q))
            q++;
        strncpy(buf, p, q - p);
        buf[q-p] = '\0';
        dist = atoi(buf);

        //
        // If aligned to the negative strand, sequence has been reverse complemented and
        // CIGAR string should be interpreted in reverse.
        //
        if (orientation == strand_plus)
            cigar.push_back(make_pair(*q, dist));
        else
            cigar.insert(cigar.begin(), make_pair(*q, dist));

        p = q + 1;
    }

    return 0;
}

int
Bam::find_start_bp_neg(int aln_bp, vector<pair<char, uint> > &cigar)
{
    uint size = cigar.size();
    char op;
    uint dist;

    for (uint i = 0; i < size; i++)  {
        op   = cigar[i].first;
        dist = cigar[i].second;

        switch(op) {
        case 'I':
            break;
        case 'S':
            if (i < size - 1)
                aln_bp += dist;
            break;
        case 'M':
        case 'D':
            aln_bp += dist;
            break;
        }
    }

    return aln_bp - 1;
}

int
Bam::find_start_bp_pos(int aln_bp, vector<pair<char, uint> > &cigar)
{
    char op;
    uint dist;

    op   = cigar[0].first;
    dist = cigar[0].second;

    if (op == 'S')
        aln_bp -= dist;

    return aln_bp;
}

int
Bam::edit_gaps(vector<pair<char, uint> > &cigar, char *seq)
{
    char *buf;
    uint  size = cigar.size();
    char  op;
    uint  dist, bp, len, buf_len, buf_size, j, k, stop;

    len = strlen(seq);
    bp  = 0;

    buf      = new char[len + 1];
    buf_size = len + 1;

    for (uint i = 0; i < size; i++)  {
        op   = cigar[i].first;
        dist = cigar[i].second;

        switch(op) {
        case 'S':
            stop = bp + dist;
            stop = stop > len ? len : stop;
            while (bp < stop) {
                seq[bp] = 'N';
                bp++;
            }
            break;
        case 'D':
            //
            // A deletion has occured in the read relative to the reference genome.
            // Pad the read with sufficent Ns to match the deletion, shifting the existing
            // sequence down. Trim the final length to keep the read length consistent.
            //
            k = bp >= len ? len : bp;

            strncpy(buf, seq + k, buf_size - 1);
            buf[buf_size - 1] = '\0';
            buf_len         = strlen(buf);

            stop = bp + dist;
            stop = stop > len ? len : stop;
            while (bp < stop) {
                seq[bp] = 'N';
                bp++;
            }

            j = bp;
            k = 0;
            while (j < len && k < buf_len) {
                seq[j] = buf[k];
                k++;
                j++;
            }
            break;
        case 'I':
            //
            // An insertion has occurred in the read relative to the reference genome. Delete the
            // inserted bases and pad the end of the read with Ns.
            //
            if (bp >= len) break;

            k = bp + dist > len ? len : bp + dist;
            strncpy(buf, seq + k, buf_size - 1);
            buf[buf_size - 1] = '\0';
            buf_len           = strlen(buf);

            j = bp;
            k = 0;
            while (j < len && k < buf_len) {
                seq[j] = buf[k];
                k++;
                j++;
            }

            stop = j + dist;
            stop = stop > len ? len : stop;
            while (j < stop) {
                seq[j] = 'N';
                j++;
            }
            break;
        case 'M':
            bp += dist;
            break;
        default:
            break;
        }
    }

    delete [] buf;

    return 0;
}

//#endif // HAVE_BAM
#endif // __BAMI_H__
