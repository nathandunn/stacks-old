// -*-mode:c++; c-style:k&r; c-basic-offset:4;-*-
//
// Copyright 2013, Julian Catchen <jcatchen@uoregon.edu>
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

#include "input.h"
#include "bam.h"

class Bam: public Input {
    bamFile  bam_fh;
    bam1_t  *aln;

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
	this->bam_fh = bam_open(path, "r");
	this->aln    = bam_init1();

	this->parse_header();
    };
    ~Bam() {
	bam_close(this->bam_fh);
	bam_destroy1(this->aln);
    };
    Seq *next_seq();
    int  next_seq(Seq &) { return 0; };
};

int
Bam::parse_header()
{
    bam_header_t *bamh = bam_header_init();
    bamh = bam_header_read(this->bam_fh);

    for (uint j = 0; j < (uint) bamh->n_targets; j++) {
	//
	// Record the mapping from integer ID to chromosome name that we will see in BAM records.
	//
	this->chrs[j] = string(bamh->target_name[j]);
    }

    bam_header_destroy(bamh);

    return 0;
}

Seq *
Bam::next_seq() 
{
    int bytes_read = 0;
    int flag       = 0;

    //
    // Read a record from the file, skipping unmapped reads,  and place it in a Seq object.
    //
    do {
	bytes_read = bam_read1(this->bam_fh, this->aln);

        if (bytes_read <= 0)
            return NULL;

	flag = ((this->aln->core.flag & BAM_FUNMAP) != 0);

    } while (flag == 1);

    //
    // Check which strand this is aligned to: 
    //   SAM reference: FLAG bit 0x10 - sequence is reverse complemented
    //
    flag = ((this->aln->core.flag & BAM_FREVERSE) != 0);

    //
    // If the read was aligned on the reverse strand (and is therefore reverse complemented)
    // alter the start point of the alignment to reflect the right-side of the read, at the
    // end of the RAD cut site.
    //
    // To accomplish this, we must parse the alignment CIGAR string
    //
    vector<pair<char, uint> > cigar;
    this->parse_bam_cigar(cigar, flag);

    uint bp = flag ? 
	this->find_start_bp_neg(this->aln->core.pos, cigar) : 
	this->find_start_bp_pos(this->aln->core.pos, cigar);

    //
    // Fetch the sequence.
    //
    string  seq;
    uint8_t j;
    for (int i = 0; i < this->aln->core.l_qseq; i++) {
	j = bam1_seqi(bam1_seq(this->aln), i);
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
    uint8_t *q = bam1_qual(this->aln);
    for (int i = 0; i < this->aln->core.l_qseq; i++) {
	qual += char(int(q[i]) + 33);
    }

    string chr = this->chrs[this->aln->core.tid];

    Seq *s = new Seq((const char *) bam1_qname(this->aln), seq.c_str(), qual.c_str(),
		     chr.c_str(), bp, flag ? minus : plus);

    if (cigar.size() > 0)
    	this->edit_gaps(cigar, s->seq);

    return s;
}

int 
Bam::parse_bam_cigar(vector<pair<char, uint> > &cigar, bool orientation)
{
    int  op, len;
    char c;
    uint32_t *cgr = bam1_cigar(this->aln);

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
	if (orientation == plus)
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
	if (orientation == plus)
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
    char buf[id_len];
    uint size = cigar.size();
    char op;
    uint dist, bp, len, buf_len, j, k, stop;

    len = strlen(seq);
    bp  = 0;

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
	    strncpy(buf, seq + bp, id_len - 1);
	    buf[id_len - 1] = '\0';
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
	    // inserted bases and padd the end of the read with Ns.
	    //
	    k = bp + dist;
	    strncpy(buf, seq + k, id_len - 1);
	    buf[id_len - 1] = '\0';
	    buf_len         = strlen(buf);

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
	}
    }

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
