#ifndef ALIGNMENT_H
#define ALIGNMENT_H

#include "constants.h"
#include "DNASeq4.h"
#include "stacks.h"

typedef vector<pair<char, uint>> Cigar;

class Alignment {
    const DNASeq4* seq_;
    const Cigar* cig_;

public:
    Alignment(const DNASeq4& seq, const Cigar& cigar)
        : seq_(&seq), cig_(&cigar)
        {assert(check_cigar());}

    // N.B. Inefficient; prefer iteration.
    Nt4 operator[] (size_t ref_i) const;

    string str() const {string s; for(iterator it (*this); it; ++it) s.push_back(char(*it)); return s;}

private:
    bool check_cigar() const;

public:
    // Iterator.
    // We have to use a range-style iterator to skip insertions (as we can't
    // peek at the next CIGAR operation if we don't know if we have reached the
    // end or not / we can't dereference cig_.end()).
    class iterator {
        Cigar::const_iterator cig_it_;
        Cigar::const_iterator cig_past_;
        size_t pos_; // Position in the current cigar op.
        DNASeq4::iterator seq_it_;
        DNASeq4::iterator seq_past_; // For debugging purposes; only used in assert's.

    public:
        iterator(const Alignment& a)
            : cig_it_(a.cig_->begin()), cig_past_(a.cig_->end()), pos_(0), seq_it_(a.seq_->begin()), seq_past_(a.seq_->end())
            {skip_insertion();}
        iterator& operator++ ();
        operator bool() const {return cig_it_ != cig_past_;}

        Nt4 operator* () const {if (cig_it_->first=='M') return *seq_it_; else {assert(cig_it_->first=='D'); return Nt4::n;}}

    private:
        void skip_insertion();
    };
};

struct AlnRead : Read {
    Cigar cigar;
    AlnRead(Read&& r, Cigar&& c) : Read(move(r)), cigar(move(c)) {}

    Alignment aln() const {return Alignment(seq, cigar);}
};

//
// ==================
// Inline definitions
// ==================
//

inline
Nt4 Alignment::operator[] (size_t ref_i) const {

    size_t seq_i = 0;
    for (auto op=cig_->begin(); op!=cig_->end(); ++op) {
        if (op->first == 'M') {
            if (ref_i < op->second)
                // This is the relevant cigar operation.
                return (*seq_)[seq_i+ref_i];

            // Consume ref & seq.
            seq_i += op->second;
            ref_i -= op->second;

        } else if (op->first == 'D') {
            if (ref_i < op->second)
                // This is the relevant cigar operation.
                return Nt4::n;

            // Consume ref.
            ref_i -= op->second;
        } else if (op->first == 'I') {
            // Consume seq.
            seq_i += op->second;
        } else {
            assert(false);
        }
    }
    // `ref_i` wasn't entirely consumed.
    throw std::out_of_range("Alignment::op[]: out_of_range");
    return Nt4::n;
}

inline
Alignment::iterator& Alignment::iterator::operator++ () {
    assert(cig_it_ != cig_past_);

    //
    // If the current cigar operation is M, advance in the sequence.
    // (rem. The current operation is never I, as they're skipped.)
    //
    if (cig_it_->first == 'M') {
        assert(seq_it_ != seq_past_);
        ++seq_it_;
    }

    //
    // Advance in the cigar.
    //
    ++pos_;
    if (pos_ == cig_it_->second) {
        // Enter the next CIGAR operation.
        pos_ = 0;
        ++cig_it_;
        skip_insertion();
    }
    // Upon reaching the end of the CIGAR, check that the entire sequence was
    // also consumed.
    assert(cig_it_ == cig_past_ ? !(seq_it_ != seq_past_) : true);

    return *this;
}

inline
void Alignment::iterator::skip_insertion() {
    if (cig_it_ != cig_past_ && cig_it_->first == 'I') {
        // Op is I; skip this insertion.
        for (size_t i=0; i<cig_it_->second; ++i) {
            assert(seq_it_ != seq_past_);
            ++seq_it_;
        }
        ++cig_it_;
    }
}

inline
bool Alignment::check_cigar() const {
    size_t seq_i = 0;
    for (auto& op : *cig_) {
        if (op.first == 'M' || op.first == 'I')
            // M and I consume the sequence.
            seq_i += op.second;
        else if (op.first != 'D')
            // Oops, the class only knows M, I and D.
            return false;
    }

    return seq_i == seq_->length();
}

#endif
