// -*-mode:c++; c-style:k&r; c-basic-offset:4;-*-
//
// Copyright 2010, Julian Catchen <jcatchen@uoregon.edu>
//
// This file is part of Stacks.
//
// Stacks is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// Stacks is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with Stacks.  If not, see <http://www.gnu.org/licenses/>.
//

#ifndef __STACKS_H__
#define __STACKS_H__

#include <cstring>
#include <string>
#include <vector>
#include <map>
#include <set>
#include <algorithm>
#include <utility>
#include<iostream>
#include<sstream>

#include "constants.h"
#include "nucleotides.h"
#include "Seq.h"
#include "DNASeq.h"
#include "DNANSeq.h"
#include "DNASeq4.h"
#include "MetaPopInfo.h"

typedef unsigned int uint;
typedef string allele_type;

enum snp_type    {snp_type_het, snp_type_hom, snp_type_unk};
enum read_type   {primary, secondary};
enum searcht     {sequence, genomic_loc};
enum corr_type   {p_value, bonferroni_win, bonferroni_gen, no_correction};

class SNP {
 public:
    snp_type type;   // Heterozygous or homozygous
    uint     col;
    float    lratio;
    char     rank_1;
    char     rank_2;
    char     rank_3;
    char     rank_4;

    SNP() {
        col    = 0;
        lratio = 0.0;
        rank_1 = 0;
        rank_2 = 0;
        rank_3 = 0;
        rank_4 = 0;
    }
};

class Gap {
public:
    uint start;
    uint end;

    Gap(uint s, uint e) {
        start = s;
        end   = e;
    }
};

class Aln {
public:
    uint   id;
    uint   gap_cnt;
    double pct_id;
    string cigar;
    Aln() {
        this->id      = 0;
        this->pct_id  = 0.0;
        this->gap_cnt = 0;
    }
    Aln(uint id, string cigar, double pct_id, uint gap_cnt) {
        this->id      = id;
        this->cigar   = cigar;
        this->pct_id  = pct_id;
        this->gap_cnt = gap_cnt;
    }
};

class PStack {
 public:
    uint            id;
    uint         count; // Number of identical reads forming this stack
    DNANSeq       *seq; // Sequence read
    vector<char *> map; // List of sequence read IDs merged into this stack
    PhyLoc         loc; // Physical genome location of this stack.

    PStack()  {
        id     = 0;
        count  = 0;
        seq    = NULL;
    }
    PStack(const PStack& other);
    PStack& operator= (PStack&& other);
    PStack& operator= (const PStack& other) = delete;

    ~PStack() {
        if (seq!=NULL)
            delete seq;
        for (unsigned int i = 0; i < this->map.size(); i++)
            delete [] this->map[i];
    }

    int  add_id(const char *);
    int  add_seq(const char *);
    int  add_seq(const DNANSeq *);
    void add_read(const char* read_name) {
        char* copy = new char[strlen(read_name)+1];
        strcpy(copy, read_name);
        map.push_back(copy);
        ++count;
    }

    // extend(): Extends the PStack to the desired span.
    void extend(const PhyLoc& phyloc, int length);

    void clear();
    bool operator< (const PStack& other) const;

    static void set_id_of(set<PStack>::iterator pstack, int id) {
        const_cast<PStack&>(*pstack).id = id;
    }
    static void add_read_to(set<PStack>::iterator pstack, const char* read_name) {
        const_cast<PStack&>(*pstack).add_read(read_name);
    }

};

class Stack {
 public:
    uint         id;
    DNANSeq     *seq; // Sequence read
    vector<uint> map; // List of sequence read IDs merged into this stack

    Stack()  {
        id  = 0;
        seq = NULL;
    }
    ~Stack() {
        delete this->seq;
    }
    uint count() { return this->map.size(); }
    int  add_id(uint);
    int  add_seq(const char *);
    int  add_seq(const DNANSeq *);
};

class Rem {
 public:
    uint          id;
    vector<uint> map; // List of sequence read IDs merged into this stack
    DNANSeq     *seq; // Sequence read
    bool    utilized;

    Rem();
    Rem(int, uint, DNANSeq *);
    ~Rem() {
        delete this->seq;
    }
    uint count() { return this->map.size(); }
    int  add_id(uint);
    int  add_seq(const char *);
    int  add_seq(const DNANSeq *);
};

class CatMatch {
public:
    int    batch_id;
    int    cat_id; // c-locus ID
    int    sample_id; // sample ID
    int    tag_id; // s-locus ID
    int    depth;
    double lnl;
    char  *haplotype;
    char  *cigar;

    CatMatch() {
        batch_id  = 0;
        cat_id    = 0;
        sample_id = 0;
        tag_id    = 0;
        depth     = 0;
        lnl       = 0.0;
        haplotype = NULL;
        cigar     = NULL;
    }
    ~CatMatch() {
        delete [] haplotype;
        delete [] cigar;
    }
};

class ModRes {
public:
    int   sample_id;
    int   tag_id;
    char *model;

    ModRes(int samp_id, int tag_id, const char *model) {
        this->sample_id = samp_id;
        this->tag_id    = tag_id;
        this->model     = new char[strlen(model) + 1];
        strcpy(this->model, model);
    }
    ~ModRes() {
        delete [] this->model;
    }
};

class SNPRes {
public:
    int   sample_id;
    int   tag_id;
    vector<SNP *> snps;

    SNPRes(int samp_id, int tag_id) {
        this->sample_id = samp_id;
        this->tag_id    = tag_id;
    }
    ~SNPRes() {
        for (uint i = 0; i < this->snps.size(); i++)
            delete this->snps[i];
        this->snps.clear();
    }
};

struct Read {
    DNASeq4 seq;
    string name;

    Read(DNASeq4&& s, string&& n)
        : seq(move(s)), name(move(n))
        {}
    Read(Read&&) = default;
    Read& operator= (Read&&) = default;

    bool is_read2() const;
};

//
// Counts: A class to store nucleotide counts.
// e.g. Counts<Nt2> is a std::array of {n_A, n_C, n_G, n_T}.
//
template<typename Nt>
class Counts {
    // Array of counts, containing the count of A's at index Nt::a,of C's at
    // index Nt::c, etc.
    array<size_t,Nt::max()+1> counts_;

public:
    Counts() {
        for (size_t& c : counts_)
            c=-1;
        for (Nt nt : Nt::all)
            counts_[size_t(nt)] = 0;
    }
    Counts(const Counts<Nt4>& nt4counts);
    Counts(const Counts<Nt2>& nt2counts);

    void clear() {for (Nt nt : Nt::all) counts_[size_t(nt)]=0;}
    void increment(Nt nt) {++counts_[size_t(nt)];}
    void increment(Nt nt, size_t cnt) {counts_[size_t(nt)] += cnt;}

    // Get the count for a given nucleotide.
    size_t operator[] (Nt nt) const {return counts_[size_t(nt)];}
    const array<size_t,Nt::max()+1>& arr() const {return counts_;}

    size_t sum() const {return (*this)[Nt::a] + (*this)[Nt::c] + (*this)[Nt::g] + (*this)[Nt::t];}
    array<pair<size_t,Nt>,4> sorted() const;

    Counts& operator+= (const Counts& other)
        {for (Nt nt : Nt::all) counts_[size_t(nt)] += other.counts_[size_t(nt)]; return *this;}

    // Print the counts.
    template<typename Nt_> friend ostream& operator<< (ostream& os, const Counts<Nt_>& cnts);
};

struct SiteCounts {
    Counts<Nt2> tot; // The sum over all samples.
    vector<Counts<Nt2>> samples; // With size() == mpopi->samples().size()/
    const MetaPopInfo* mpopi;
};

//
// GtLiks: A class to store the likelihoods of SNP genotypes.
//
class GtLiks {
    array<double,10> lnliks_; // {AA,AC,CC,AG,CG,GG,AT,CT,GT,TT} similar to VCF.
public:
    GtLiks() : lnliks_{{1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0}} {}
    double at(Nt2 n1, Nt2 n2) const {return at(gt_index(n1,n2));}
    bool has_lik(Nt2 n1, Nt2 n2) const {return has_lik(gt_index(n1,n2));}
    void set(Nt2 n1, Nt2 n2, double lnl) {set(gt_index(n1, n2), lnl);}

    double at(size_t gt) const {assert(has_lik(gt)); return lnliks_[gt];}
    bool has_lik(size_t gt) const {return lnliks_[gt] != 1.0;}
    void set(size_t gt, double lnl) {assert(!std::isnan(lnl)); assert(lnl<=0.0); assert(!has_lik(gt)); lnliks_[gt] = lnl;}

    static size_t gt_index(Nt2 n1, Nt2 n2) {
        if(n1<n2)
            return size_t(n1) + (size_t(n2)*(size_t(n2)+1)) / 2;
        else
            return size_t(n2) + (size_t(n1)*(size_t(n1)+1)) / 2;
    }

    // For debugging.
    friend ostream& operator<<(ostream& os, const GtLiks& liks);
};

//
// Inline definitions
// ----------
//

inline
PStack::PStack(const PStack& other)
        : id(other.id)
        , count (other.count)
        , seq (new DNANSeq(*other.seq))
        , map ()
        , loc (other.loc)
        {
    map.reserve(other.map.size());
    for (const char* readname : other.map) {
        char* copy = new char[strlen(readname)+1];
        strcpy(copy, readname);
        map.push_back(copy);
    }
}

inline
PStack& PStack::operator= (PStack&& other) {
    id = other.id;
    count = other.count;
    seq = other.seq;
    other.seq = NULL;
    swap(map, other.map);
    swap(loc, other.loc);
    return *this;
}

inline
void PStack::clear() {
    id = 0;
    count = 0;
    if(seq!=NULL) {
        delete seq;
        seq = NULL;
    }
    map.clear();
    loc.clear();
}

inline
bool PStack::operator< (const PStack& other) const {
    if (loc < other.loc)
        return true;
    else if (other.loc < loc)
        return false;
    else
        // Same genomic loci, compare sequences.
        return *seq < *other.seq;
}

template<> inline
Counts<Nt2>::Counts(const Counts<Nt2>& nt2counts) : counts_(nt2counts.counts_) {}
template<> inline
Counts<Nt2>::Counts(const Counts<Nt4>& nt4counts) : Counts() {
    for (Nt2 nt2 : Nt2::all) // We thus ignore Nt4::n.
        counts_[size_t(nt2)] = nt4counts[Nt4(nt2)];
}
template<> inline
Counts<Nt4>::Counts(const Counts<Nt4>& nt4counts) : counts_(nt4counts.counts_) {}
template<> inline
Counts<Nt4>::Counts(const Counts<Nt2>& nt2counts) : Counts() {
    for (Nt2 nt2 : Nt2::all)
        counts_[size_t(Nt4(nt2))] = nt2counts[nt2];
}

template<typename Nt>
array<pair<size_t,Nt>,4> Counts<Nt>::sorted() const {
    array<pair<size_t,Nt>,4> arr {{
        {(*this)[Nt::t], Nt::t},
        {(*this)[Nt::g], Nt::g},
        {(*this)[Nt::c], Nt::c},
        {(*this)[Nt::a], Nt::a}
    }};
    // Sort by decreasing value. Primarily on the count, secondarily on
    // the Nt4 value.
    std::sort(arr.rbegin(), arr.rend());
    return arr;
}

template<typename Nt>
ostream& operator<< (ostream& os, const Counts<Nt>& cnts) {
    bool first=true;
    for (Nt nt : Nt::all) {
        if (first)
            first = false;
        else
            os << " ";
        os << nt << ":" << cnts[nt];
    }
    return os;
}

inline
bool Read::is_read2() const {
    return name.length() >= 2
            && name.back() == '2'
            && name[name.length()-2] == '/';
}

#endif // __STACKS_H__
