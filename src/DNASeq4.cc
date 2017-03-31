#include "DNASeq4.h"

using namespace std;

const size_t Nt4::ch_to_nt4[256] = {
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

const char Nt4::nt4_to_ch[16] {
    '=','A','C','?','G','?','?','?','T','?','?','?','?','?','?','N'
};

const size_t Nt2::ch_to_nt2[256] = {
    0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0, // 0x
    0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0, // 1x
    0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0, // 2x
    0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0, // 3x
    0,a,0,c, 0,0,0,g, 0,0,0,0, 0,0,0,0, // 4x
    0,0,0,0, t,0,0,0, 0,0,0,0, 0,0,0,0, // 5x
    0,a,0,c, 0,0,0,g, 0,0,0,0, 0,0,0,0, // 6x
    0,0,0,0, t,0,0,0, 0,0,0,0, 0,0,0,0, // 7x

    0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0,
    0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0,
    0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0,
    0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0,
    0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0,
    0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0,
    0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0,
    0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0
};

const char Nt2::nt2_to_ch[4] {
    'A','C','G','T'
};

const size_t Nt2::nt4_to_nt2[16] = {
    0,a,c,0,g,0,0,0,t,0,0,0,0,0,0,0
};

DNASeq4::DNASeq4(const char* s, size_t len) : l_(len), v_() {
    v_.reserve(len/2 + len%2);
    for (size_t i=0; i<len; i+=2)
        v_.push_back(DiNuc(s[i], s[i+1])); //n.b. `s` is null-terminated
}

string DNASeq4::str() const {
    string s;
    s.reserve(l_);
    for (auto nt=begin(); nt!=end(); ++nt)
        s.push_back(*nt);
    return s;
}
