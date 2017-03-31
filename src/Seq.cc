#include <utility>

#include "constants.h"
#include "utils.h"
#include "Seq.h"

Seq::Seq() {
    this->id       = NULL;
    this->seq      = NULL;
    this->qual     = NULL;
    this->loc_str  = NULL;
    this->aln_type = AlnT::null;
    this->pct_clipped  = 0.0;
    this->map_qual = 255;
}

Seq::Seq(const Seq& other)
    : loc(other.loc) {
    if (other.id != NULL) {
        id = new char[strlen(other.id)+1];
        strcpy(id, other.id);
    } else {
        id = NULL;
    }
    if (other.seq != NULL) {
        seq = new char[strlen(other.seq)+1];
        strcpy(seq, other.seq);
    } else {
        seq = NULL;
    }
    if (other.qual != NULL) {
        qual = new char[strlen(other.qual)+1];
        strcpy(qual, other.qual);
    } else {
        qual = NULL;
    }
    if (other.loc_str != NULL) {
        loc_str = new char[strlen(other.loc_str)+1];
        strcpy(loc_str, other.loc_str);
    } else {
        loc_str = NULL;
    }

    pct_clipped  = other.pct_clipped;
    aln_type = other.aln_type;
    map_qual = other.map_qual;
}

Seq::Seq(const char *id, const char *seq) {
    this->id       = new char[strlen(id)   + 1];
    this->seq      = new char[strlen(seq)  + 1];
    this->qual     = NULL;
    this->loc_str  = NULL;

    strcpy(this->id,   id);
    strcpy(this->seq,  seq);

    this->aln_type = AlnT::null;
    this->pct_clipped  = 0.0;
    this->map_qual = 255;
}

Seq::Seq(const char *id, const char *seq, const char *qual)  {
    this->id       = new char[strlen(id)   + 1];
    this->seq      = new char[strlen(seq)  + 1];
    this->qual     = new char[strlen(qual) + 1];
    this->loc_str  = NULL;

    strcpy(this->id,   id);
    strcpy(this->seq,  seq);
    strcpy(this->qual, qual);

    this->aln_type = AlnT::null;
    this->pct_clipped  = 0.0;
    this->map_qual = 255;
}

Seq::Seq(const char *id, const char *seq, const char *qual, const char *chr, uint bp, strand_type strand)  {
    this->id      = new char[strlen(id)   + 1];
    this->qual    = new char[strlen(qual) + 1];
    this->loc_str = new char[strlen(chr)  + 15];

    strcpy(this->id,   id);
    strcpy(this->qual, qual);
    this->loc.set(chr, bp, strand);

    sprintf(this->loc_str, "%s|%d|%c", chr, bp, strand == strand_plus ? '+' : '-');

    //
    // Reverse complement sequences from the negative strand
    //
    if (strand == strand_plus) {
        this->seq = new char[strlen(seq)  + 1];
        strcpy(this->seq, seq);
    } else {
        this->seq = rev_comp(seq);
    }

    this->aln_type = AlnT::primary;
    this->pct_clipped  = 0.0;
    this->map_qual = 255;
}

Seq::Seq(const char *id, const char *seq, const char *qual, const char *chr, uint bp, strand_type strand, AlnT aln_type, double pct_clipped, int map_qual)  {
    this->id      = new char[strlen(id)   + 1];
    this->qual    = new char[strlen(qual) + 1];
    this->loc_str = new char[strlen(chr)  + 15];

    strcpy(this->id,   id);
    strcpy(this->qual, qual);
    this->loc.set(chr, bp, strand);

    sprintf(this->loc_str, "%s|%d|%c", chr, bp, strand == strand_plus ? '+' : '-');

    //
    // Reverse complement sequences from the negative strand
    //
    if (strand == strand_plus) {
        this->seq = new char[strlen(seq)  + 1];
        strcpy(this->seq, seq);
    } else {
        this->seq = rev_comp(seq);
    }

    this->aln_type = aln_type;
    this->pct_clipped  = pct_clipped;
    this->map_qual = map_qual;
}

void swap(Seq& s1, Seq& s2) {
    char  *ptr;
    AlnT   a;
    double p;

    ptr = s1.id;
    s1.id = s2.id;
    s2.id = ptr;

    ptr = s1.seq;
    s1.seq = s2.seq;
    s2.seq = ptr;

    ptr = s1.qual;
    s1.qual = s2.qual;
    s2.qual = ptr;

    ptr = s1.loc_str;
    s1.loc_str = s2.loc_str;
    s2.loc_str = ptr;

    a = s1.aln_type;
    s1.aln_type = s2.aln_type;
    s2.aln_type = a;

    p = s1.pct_clipped;
    s1.pct_clipped = s2.pct_clipped;
    s2.pct_clipped = p;

    swap(s1.loc, s2.loc);

    std::swap(s1.map_qual, s2.map_qual);
}
