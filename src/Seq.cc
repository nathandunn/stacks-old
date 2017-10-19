#include <utility>

#include "constants.h"
#include "utils.h"
#include "Seq.h"

const char PhyLoc::empty_str[1] = {'\0'};

PhyLoc::PhyLoc(const string& s) : PhyLoc() {
    try {
        // Chromosome.
        const char* p = strchr(s.c_str(), ':');
        if (p == NULL || p == s.c_str())
            throw exception();
        chr_ = new char[p-s.c_str()+1];
        strncpy(chr_, s.c_str(), p-s.c_str());
        chr_[p-s.c_str()] = '\0';

        // BP.
        ++p;
        char* end;
        bp = strtol(p, &end, 10) - 1;
        if (end == p)
            throw exception();

        // Strand.
        if (*end != '\0') {
            if (*end != ':')
                throw exception();
            p = end + 1;
            if (!( (*p == '+' || *p == '-') && *(p+1) == '\0' ) )
                throw exception();
            strand = (*p == '+') ? strand_plus : strand_minus;
        }
    } catch (exception&) {
        cerr << "Error: Malformed genomic position '" << s << "'.\n";
        throw;
    }
}

Seq::Seq() {
    this->id       = NULL;
    this->comment  = NULL;
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
    if (other.comment != NULL) {
        comment = new char[strlen(other.comment)+1];
        strcpy(comment, other.comment);
    } else {
        comment = NULL;
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
    this->comment  = NULL;
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
    this->comment  = NULL;
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
    this->comment  = NULL;
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
    this->comment = NULL;

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
