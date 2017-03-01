#ifndef ALIGNMENT_H
#define ALIGNMENT_H

#include "constants.h"
#include "DNASeq4.h"
#include "stacks.h"

class Alignment {
public:
    typedef std::vector<std::pair<char, size_t>> Cigar;

private:
    Cigar cig_;
    const DNASeq4& seq_;

public:
    Alignment(const DNASeq4& seq, Cigar&& cigar)
        : cig_(std::move(cigar)), seq_(seq)
        {}
    Alignment(Alignment&&) = default; // Trivial move constructor.
    Alignment(const Alignment&) = delete; // Forbid copying.
    Alignment& operator= (Alignment&&) = default;
    Alignment& operator= (const Alignment&) = delete;

    // N.B. Inefficient; prefer iteration.
    size_t operator[] (size_t ref_i) const;

    // Iterator.
    // We have to use a range-style iterator to let operator++ skip insertions
    // (as we can't peek at the next CIGAR operation if we don't know if we
    // have reached the end or not / we can't dereference cig_.end()).
    class range_iterator {
        Cigar::const_iterator cig_it_;
        Cigar::const_iterator cig_past_;
        size_t pos_; // Position in the current cigar op.
        DNASeq4::iterator seq_it_;
        DNASeq4::iterator seq_past_; // For debugging purposes; only used in assert's.

    public:
        range_iterator(const Alignment& a)
            : cig_it_(a.cig_.begin()), cig_past_(a.cig_.end()), pos_(0), seq_it_(a.seq_.begin()), seq_past_(a.seq_.end())
            {}
        range_iterator& operator++ ();
        operator bool() const {return cig_it_ != cig_past_;}

        size_t nt() const {if (cig_it_->first=='M') return seq_it_.nt(); else {assert(cig_it_->first=='D'); return Nt4::n;}}
        char operator* () const {return Nt4::nt_to_ch[nt()];}
    };
};

struct AlnRead : Read {
    Alignment aln;
    AlnRead(Read&& r, Alignment&& a) : Read(std::move(r)), aln(std::move(a)) {}
};

// AlnSite
// A site in a multiple alignment.
class AlnSite {
    const vector<Alignment::range_iterator>& col_; // One iterator per read in the alignment.

public:
    AlnSite(const vector<Alignment::range_iterator>& column) : col_(column) {}
    void counts(Nt4Counts& counts) const {counts.reset(); for (auto& read: col_) counts.increment(read.nt()); counts.sort();}
};


//
// ==================
// Inline definitions
// ==================
//

inline
size_t Alignment::operator[] (size_t ref_i) const {
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
            throw std::out_of_range(std::string("out_of_range in Alignment::op[]: +")+std::to_string(ref_i));
    }

    seq_i += ref_i;
    return seq_[seq_i];
}

inline
Alignment::range_iterator& Alignment::range_iterator::operator++ () {
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
        if (cig_it_ != cig_past_) {
            if (cig_it_->first == 'I') {
                // Op is I; skip this insertion.
                for (size_t i=0; i<cig_it_->second; ++i) {
                    assert(seq_it_ != seq_past_);
                    ++seq_it_;
                }
                ++cig_it_;
            }
        } else {
            // End of the CIGAR; check that the entire sequence was consumed.
            assert(!(seq_it_ != seq_past_));
        }
    }

    return *this;
}

#endif
