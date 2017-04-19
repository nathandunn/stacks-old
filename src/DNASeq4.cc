#include "DNASeq4.h"

DNASeq4::DNASeq4(const char* s, size_t len) : l_(len), v_() {
    reserve(l_);
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
    rev.reserve(rev.l_);

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

void DNASeq4::append(const DNASeq4& other) {

    if (other.l_ == 0)
        return;

    auto nt = other.begin();
    if (l_%2==1) {
        v_.back().second(*nt);
        ++nt;
    }

    l_ += other.l_;
    reserve(l_);

    while(nt != other.end()) {
        Nt4 first = *nt;
        ++nt;
        if (nt != other.end()) {
            v_.push_back(DiNuc(first, *nt));
            ++nt;
        } else {
            v_.push_back(DiNuc(first, Nt4(0)));
            break;
        }
    }
}
