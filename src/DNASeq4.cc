#include "DNASeq4.h"

using namespace std;

const size_t Nt4::ch_to_nt[256] = {
    0,n,n,n, n,n,n,n, n,n,n,n, n,n,n,n, // 0x
    n,n,n,n, n,n,n,n, n,n,n,n, n,n,n,n, // 1x
    n,n,n,n, n,n,n,n, n,n,n,n, n,n,n,n, // 2x
    n,n,n,n, n,n,n,n, n,n,n,n, n,0,n,n, // 3x rem. 0 for 0x3D ('=') because htslib does it...
    n,a,n,c, n,n,n,g, n,n,n,n, n,n,n,n, // 4x
    n,n,n,n, t,n,n,n, n,n,n,n, n,n,n,n, // 5x
    n,a,n,c, n,n,n,g, n,n,n,n, n,n,n,n, // 6x
    n,n,n,n, t,n,n,n, n,n,n,n, n,n,n,n, // 7x

    n,n,n,n, n,n,n,n, n,n,n,n, n,n,n,n,
    n,n,n,n, n,n,n,n, n,n,n,n, n,n,n,n,
    n,n,n,n, n,n,n,n, n,n,n,n, n,n,n,n,
    n,n,n,n, n,n,n,n, n,n,n,n, n,n,n,n,
    n,n,n,n, n,n,n,n, n,n,n,n, n,n,n,n,
    n,n,n,n, n,n,n,n, n,n,n,n, n,n,n,n,
    n,n,n,n, n,n,n,n, n,n,n,n, n,n,n,n,
    n,n,n,n, n,n,n,n, n,n,n,n, n,n,n,n
};

const char Nt4::nt_to_ch[16] {
    '=','A','C','?','G','?','?','?','T','?','?','?','?','?','?','N'
};

DNASeq4::DNASeq4(const char* s, size_t len) : l_(len), v_() {
    v_.reserve(len/2 + len%2);
    for (size_t i=0; i<len; i+=2)
        v_.push_back(DiNuc(s[i], s[i+1])); //n.b. `s` is null-terminated
}
