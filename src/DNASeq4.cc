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
    DNASeq4 rev;
    rev.l_ = l_;
    rev.reserve(rev.l_);

    iterator nt = end();
    if (nt != begin()) {
        do {
            --nt;
            Nt4 first = (*nt).rev_compl();
            if (nt != begin()) {
                --nt;
                Nt4 second = (*nt).rev_compl();
                rev.v_.push_back(DiNuc(first, second));
            } else {
                rev.v_.push_back(DiNuc((*nt).rev_compl(), Nt4(0)));
                break;
            }
        } while(nt != begin());
    }

    return rev;
}

void DNASeq4::append(iterator first, iterator past) {

    if (!(first != past))
        return;

    if (l_%2==1) {
        v_.back().second(*first);
        ++l_;
        ++first;
    }

    l_ += past - first;
    reserve(l_);

    while(first != past) {
        Nt4 prev = *first;
        ++first;
        if (first != past) {
            v_.push_back(DiNuc(prev, *first));
            ++first;
        } else {
            v_.push_back(DiNuc(prev, Nt4(0)));
            break;
        }
    }
}
