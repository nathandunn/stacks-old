#ifndef ALIGNMENT_H
#define ALIGNMENT_H

#include "constants.h"
#include "DNASeq4.h"
#include "stacks.h"

typedef vector<pair<char, uint>> Cigar;

class Alignment {
    const DNASeq4* seq_;
    Cigar cig_;

public:
    Alignment(const DNASeq4& seq, Cigar&& cigar)
        : seq_(&seq), cig_(move(cigar))
        {assert(check_cigar());}

    // N.B. Copying of course copies the pointer as well; this is not always
    // what we want, c.f. AlnRead
    Alignment(Alignment&&) = default;
    Alignment& operator= (Alignment&&) = default;
    Cigar&& move_cigar() {return move(cig_);}
    void assign(const DNASeq4& seq, Cigar&& cigar) {seq_ = &seq; cig_ = move(cigar);}

    // N.B. Inefficient; prefer iteration.
    Nt4 operator[] (size_t ref_i) const;

    const Cigar& cigar() const {return cig_;}
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
            : cig_it_(a.cig_.begin()), cig_past_(a.cig_.end()), pos_(0), seq_it_(a.seq_->begin()), seq_past_(a.seq_->end())
            {skip_insertion();}
        iterator& operator++ ();
        operator bool() const {return cig_it_ != cig_past_;}

        Nt4 operator* () const {if (cig_it_->first=='M') return *seq_it_; else {assert(cig_it_->first=='D'); return Nt4::n;}}

    private:
        void skip_insertion();
    };
};

struct AlnRead : Read {
    Alignment aln;
    AlnRead(Read&& r, Cigar&& c) : Read(move(r)), aln(seq, move(c)) {}

    AlnRead(AlnRead&& other)
        : Read(move(other)), aln(seq, other.aln.move_cigar()) //n.b. `seq` is `this->seq`.
        {}
    AlnRead& operator= (AlnRead&& other)
        {Read::operator=(move(other)); aln.assign(seq, other.aln.move_cigar()); return *this;}
};

//
// ==================
// Inline definitions
// ==================
//

inline
Nt4 Alignment::operator[] (size_t ref_i) const {
    size_t seq_i = 0;
    auto op = cig_.begin();
    while (ref_i >= op->second || op->first == 'I') {
        if(op->first == 'M') {
            // Consumes ref & seq.
            seq_i += op->second;
            ref_i -= op->second;
        } else if (op->first == 'D') {
            // Consumes ref.
            ref_i -= op->second;
        } else if (op->first == 'I') {
            // Consumes seq.
            seq_i += op->second;
        }
        ++op;
        if (op == cig_.end())
            throw std::out_of_range(string("out_of_range in Alignment::op[]: +")+to_string(ref_i));
    }

    seq_i += ref_i;
    return (*seq_)[seq_i];
}

inline
Alignment::iterator& Alignment::iterator::operator++ () {
    assert(cig_it_ != cig_past_);

    if (cig_it_->first == 'M') {
        assert(seq_it_ != seq_past_);
        ++seq_it_;
    }

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
    for (auto& op : cig_) {
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
