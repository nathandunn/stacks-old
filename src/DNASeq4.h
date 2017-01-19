#include "hts.h" // For `seq_nt16_table`, `seq_nt16_str`.

#include "constants.h"

// A dinucleotide, coded on one byte (4 bits per nucleotide).
// 1,2,4,8,15 are A,C,T,G,N like in Htslib, and the high-bits
// correspond to the first nucleotide.
class DiNuc {
    uchar u_;

public:
    DiNuc() : u_(0) {}
    DiNuc(const DiNuc& other) : u_(other.u_) {}
    DiNuc(char c1, char c2) : u_(0) {set_first(c2u(c1)); set_first(c2u(c2));}

    uchar ufirst() {return u_ & lowbits;}
    uchar usecond() {return u_ >>4;}

    char first() {return u2c(ufirst());}
    char second() {return u2c(usecond());}

    static uchar c2u(char c) {return c2u_table[c];}
    static char u2c(uchar u) {return u2c_table[u];}

private:
    static const uchar lowbits = 15;
    static const uchar highbits = ~15;
    static const uchar DiNuc::c2u_table[256];
    static const char DiNuc::u2c_table[16];

    void set_first(uchar u) {u_ |= u;}
    void set_second(uchar u) {u_ |= u <<4;}

    void clear_first() {u_ &= lowbits;}
    void clear_second() {u_ &= highbits;}
};
