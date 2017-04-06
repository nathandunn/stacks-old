#include "DNASeq4.h"

using namespace std;

const Nt4 Nt4::a (1);  // 0001 (1)
const Nt4 Nt4::c (2);  // 0010 (2)
const Nt4 Nt4::g (4);  // 0100 (4)
const Nt4 Nt4::t (8);  // 1000 (8)
const Nt4 Nt4::n (15); // 1111 (15)
const vector<Nt4> Nt4::all = {a, c, g, t, n};

const Nt4 Nt4::from_ch[256] = {
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

const char Nt4::to_ch[16] {
    '=','A','C','?','G','?','?','?','T','?','?','?','?','?','?','N'
};

const Nt4 Nt4::rev_compl_[16] {
    n,t,g,n,c,n,n,n,a,n,n,n,n,n,n,n
};

const Nt2 Nt2::a (0);
const Nt2 Nt2::c (1);
const Nt2 Nt2::g (2);
const Nt2 Nt2::t (3);
const vector<Nt2> Nt2::all = {a, c, g, t};

const Nt2 Nt2::from_ch[256] = {
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

const Nt2 Nt2::from_nt4[16] = {
    0,a,c,0,g,0,0,0,t,0,0,0,0,0,0,0
};

const char Nt2::to_ch[4] {
    'A','C','G','T'
};

const Nt2 Nt2::rev_compl_[4] {
    t,g,c,a
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
        s.push_back(char(*nt));
    return s;
}

DNASeq4 DNASeq4::rev_compl() const {
    assert(v_.size() == (l_+1)/2);

    DNASeq4 rev;
    rev.l_ = l_;
    rev.v_.reserve(v_.size());

    iterator nt = end();
    for (size_t i=0; i<l_/2; ++i)
        // Push two nucleotides, l_/2 times.
        rev.v_.push_back(
                DiNuc((*--nt).rev_compl(), (*--nt).rev_compl())
                );

    if (l_ % 2 == 1)
        rev.v_.push_back(DiNuc((*--nt).rev_compl(), Nt4(0)));

    assert(rev.v_.size() == v_.size() && !(nt != begin()));
    return rev;
}
