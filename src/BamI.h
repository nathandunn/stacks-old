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

//
// Code to parse binary BAM format. This format is created for
// reads that have been aligned to a reference genome.
//

#ifdef HAVE_BAM

#include <algorithm>

#include "sam.h" //htslib

#include "stacks.h"
#include "input.h"

int bam_find_start_bp(int, strand_type, const vector<pair<char, uint> > &);
int bam_edit_gaps(vector<pair<char, uint> > &, char *);

class Bam: public Input {
private:
    htsFile   *bam_fh;
    bam_hdr_t *bamh;
    bam1_t    *aln;

    map<uint, string> chrs;

    int parse_header();
    void parse_cigar(vector<pair<char, uint> > &);

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

class BamRecord {
    bam1_t* r_;

public:
    BamRecord() : r_(bam_init1()) {}
    ~BamRecord() {bam_destroy1(r_);}

    string qname() const;
    uint16_t flag() const;
    int32_t chrom() const;
    int32_t pos() const;
    uint8_t mapq() const;
    vector<pair<char, uint>> cigar() const;
    // (rnext)
    // (pnext)
    // (tlen)
    //DNASeq4  seq() const; // xxx (Feb2017) Import DNASeq4 from n_pe
    // (qual)
    vector<string> aux() const;

    AlnT aln_type();
    bool is_rev_compl();
    bool is_paired();
    bool is_fw_read();
    bool is_rev_read();

    //friend bool Bam::next();
};

//
// Inline definitions
// ----------
//

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

    //
    // Read a record from the file, skipping unmapped reads,  and place it in a Seq object.
    //
    if (sam_read1(bam_fh, bamh, aln) < 0)
        return false;

    //
    // Parse the type of the record.
    //
    AlnT aln_type;
    if (aln->core.flag & BAM_FUNMAP)
        aln_type = AlnT::null;
    else if (aln->core.flag & BAM_FSECONDARY)
        aln_type = AlnT::secondary;
    else if (aln->core.flag & BAM_FSUPPLEMENTARY)
        aln_type = AlnT::supplementary;
    else
        aln_type = AlnT::primary;

    //
    // Fetch the sequence.
    //
    string  seq;
    seq.reserve(aln->core.l_qseq);
    for (int i = 0; i < aln->core.l_qseq; i++) {
        uint8_t j = bam_seqi(bam_get_seq(aln), i);
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
        default:
            // assert(false);
            break;
        }
    }

    //
    // Fetch the quality score.
    //
    string   qual;
    uint8_t *q = bam_get_qual(aln);
    for (int i = 0; i < this->aln->core.l_qseq; i++) {
        qual += char(int(q[i]) + 33);
    }

    if (aln_type == AlnT::null) {
        s = Seq(bam_get_qname(aln), seq.c_str(), qual.c_str());
    } else {
        // Fetch the chromosome.
        string chr = chrs[aln->core.tid];

        //
        // Check which strand this is aligned to:
        //   SAM reference: FLAG bit 0x10 - sequence is reverse complemented
        //
        strand_type strand = aln->core.flag & BAM_FREVERSE ? strand_minus : strand_plus;

        //
        // Parse the alignment CIGAR string.
        // If aligned to the negative strand, sequence has been reverse complemented and
        // CIGAR string should be interpreted in reverse.
        //
        vector<pair<char, uint> > cigar;
        parse_cigar(cigar);
        if (strand == strand_minus)
            std::reverse(cigar.begin(), cigar.end());

        //
        // If the read was aligned on the reverse strand (and is therefore reverse complemented)
        // alter the start point of the alignment to reflect the right-side of the read, at the
        // end of the RAD cut site.
        //

        uint bp = bam_find_start_bp(aln->core.pos, strand, cigar);

        //
        // Calculate the percentage of the sequence that was aligned to the reference.
        //
        uint clipped = 0;
        for (auto& op : cigar)
            if (op.first == 'S')
                clipped += op.second;
        double pct_clipped = (double) clipped / seq.length();

        s = Seq(bam_get_qname(aln), seq.c_str(), qual.c_str(),
                chr.c_str(), bp, strand, aln_type, pct_clipped, aln->core.qual);

        if (cigar.size() > 0)
            bam_edit_gaps(cigar, s.seq);
    }

    return true;
}

void
Bam::parse_cigar(vector<pair<char, uint> > &cigar)
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
        case BAM_CEQUAL:
            c = '=';
            break;
        case BAM_CDIFF:
            c = 'X';
            break;
        case BAM_CINS:
            c = 'I';
            break;
        case BAM_CDEL:
            c = 'D';
            break;
        case BAM_CSOFT_CLIP:
            c = 'S';
            break;
        case BAM_CREF_SKIP:
            c = 'N';
            break;
        case BAM_CHARD_CLIP:
            c = 'H';
            break;
        case BAM_CPAD:
            c = 'P';
            break;
        default:
            cerr << "Warning: Unknown CIGAR operation (current read name is '" << bam_get_qname(aln) << "').\n";
            break;
        }

        cigar.push_back(make_pair(c, len));
    }
}

int
bam_find_start_bp(int aln_bp, strand_type strand, const vector<pair<char, uint> > &cigar)
{
    if (strand == strand_plus) {
        if (cigar.at(0).first == 'S')
            aln_bp -= cigar.at(0).second;
    } else {
        // assert(strand == strand_minus);
        for (uint i = 0; i < cigar.size(); i++)  {
            char op   = cigar[i].first;
            uint dist = cigar[i].second;

            switch(op) {
            case 'I':
            case 'H':
                break;
            case 'S':
                if (i < cigar.size() - 1)
                    aln_bp += dist;
                break;
            case 'M':
            case '=':
            case 'X':
            case 'D':
            case 'N':
                aln_bp += dist;
                break;
            default:
                break;
            }
        }
        aln_bp -= 1;
    }

    return aln_bp;
}

int
bam_edit_gaps(vector<pair<char, uint> > &cigar, char *seq)
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
        case '=':
        case 'X':
            bp += dist;
            break;
        default:
            break;
        }
    }

    delete [] buf;

    return 0;
}

#else  // If HAVE_BAM is undefined and BAM library is not present.

#include "input.h"

class Bam: public Input {
 public:
    Bam(const char *path) : Input() { cerr << "BAM support was not enabled when Stacks was compiled.\n"; };
    ~Bam() {};
    Seq *next_seq()      { return NULL; };
    int  next_seq(Seq &) { return 0; };
};

#endif // HAVE_BAM

#endif // __BAMI_H__
