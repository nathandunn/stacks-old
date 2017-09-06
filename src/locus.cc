// -*-mode:c++; c-style:k&r; c-basic-offset:4;-*-OA
//
// Copyright 2013-2015, Julian Catchen <jcatchen@illinois.edu>
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

//
// locus.cc -- routines for the Locus class and its derivatives.
//
#include "locus.h"

uint
Locus::sort_bp(uint k) const
{
    if (this->loc.strand == strand_plus)
        return this->loc.bp + k;
    else
        return (k == 0 ? this->loc.bp - this->len + 1 : this->loc.bp - k);
}

int
Locus::snp_index(uint col) const
{
    for (uint i = 0; i < this->snps.size(); i++)
        if (this->snps[i]->col == col)
            return i;
    return -1;
}

int
Locus::add_consensus(const char *seq)
{
    if (this->con != NULL)
        delete [] this->con;

    this->len = strlen(seq);
    this->con = new char[this->len + 1];
    strcpy(this->con, seq);

    return 0;
}

int
Locus::add_model(const char *seq)
{
    if (this->model != NULL)
        delete [] this->model;

    this->model = new char[this->len + 1];
    strncpy(this->model, seq, this->len);
    this->model[this->len] = '\0';

    return 0;
}

int
Locus::populate_alleles()
{
    vector<SNP *>::iterator  i;
    map<string, int>::iterator j;
    string s;
    uint   k;

    if (this->len > strlen(this->con))
	cerr << "Recorded locus->len: " << this->len << "; consensus length: " << strlen(this->con) << "\n";
    
    //
    // Is this effective?
    //
    for (uint n = 0; n < this->strings.size(); n++) {
        this->strings[n].first.clear();
        this->strings[n].second.clear();
    }
    this->strings.clear();

    if (this->snps.size() == 0) {
        this->strings.push_back(make_pair("consensus", this->con));
        return 0;
    }

    for (j = this->alleles.begin(); j != this->alleles.end(); j++) {
        s = this->con;
        k = 0;

        for (i = this->snps.begin(); i != this->snps.end(); i++) {
            if ((*i)->type == snp_type_het && (*i)->col < this->len && k < j->first.length()) {
                s.replace((*i)->col, 1, 1, j->first[k]);
                k++;
            }
        }

        this->strings.push_back(make_pair(j->first, s));
    }

    return 0;
}

bool
bp_compare(Locus *a, Locus *b)
{
    return (a->sort_bp() < b->sort_bp());
}

int
adjust_snps_for_gaps(Cigar &cigar, Locus *loc)
{
    uint   size = cigar.size();
    char   op;
    uint   dist, bp, stop, offset, snp_index;

    bp        = 0;
    offset    = 0;
    snp_index = 0;

    for (uint i = 0; i < size; i++)  {
        op   = cigar[i].first;
        dist = cigar[i].second;

        switch(op) {
        case 'D':
            offset += dist;
            break;
        case 'I':
        case 'M':
        case 'S':
            stop = bp + dist;
            while (bp < stop && snp_index < loc->snps.size()) {
                if (loc->snps[snp_index]->col == bp) {
                    loc->snps[snp_index]->col += offset;
                    snp_index++;
                }
                bp++;
            }
            break;
        default:
            break;
        }
    }

    return 0;
}

int
adjust_and_add_snps_for_gaps(Cigar &cigar, Locus *loc)
{
    uint   size = cigar.size();
    char   op;
    uint   dist, bp, new_bp, stop, snp_cnt;
    SNP   *s;

    bp      = 0;
    new_bp  = 0;
    snp_cnt = loc->snps.size();

    vector<SNP *> snps;

    for (uint i = 0; i < size; i++)  {
        op   = cigar[i].first;
        dist = cigar[i].second;

        switch(op) {
        case 'D':
            stop = new_bp + dist;
            while (new_bp < stop) {
                s = new SNP;
                s->col    = new_bp;
                s->type   = snp_type_unk;
                s->rank_1 = 'N';
                snps.push_back(s);
                new_bp++;
            }
            break;
        case 'I':
        case 'M':
        case 'S':
            stop = bp + dist > snp_cnt ? snp_cnt : bp + dist;
            while (bp < stop) {
                loc->snps[bp]->col = new_bp;
                snps.push_back(loc->snps[bp]);
                bp++;
                new_bp++;
            }
            break;
        default:
            break;
        }
    }

    loc->snps.clear();

    for (uint i = 0; i < snps.size(); i++)
        loc->snps.push_back(snps[i]);

    return 0;
}

int
remove_snps_from_gaps(Cigar &cigar, Locus *loc)
{
    uint   size = cigar.size();
    char   op;
    uint   dist, bp, new_bp, stop, snp_cnt;

    bp      = 0;
    new_bp  = 0;
    snp_cnt = loc->snps.size();

    vector<SNP *> snps;

    for (uint i = 0; i < size; i++)  {
        op   = cigar[i].first;
        dist = cigar[i].second;

        switch(op) {
        case 'D':
            stop = bp + dist;
            while (bp < stop) {
                delete loc->snps[bp];
                bp++;
            }
            break;
        case 'I':
        case 'M':
        case 'S':
            stop = bp + dist > snp_cnt ? snp_cnt : bp + dist;
            while (bp < stop) {
                loc->snps[bp]->col = new_bp;
                snps.push_back(loc->snps[bp]);
                bp++;
                new_bp++;
            }
            break;
        default:
            break;
        }
    }

    loc->snps.clear();

    for (uint i = 0; i < snps.size(); i++)
        loc->snps.push_back(snps[i]);

    return 0;
}

QLocus::~QLocus()
{
    vector<Match *>::iterator it;

    for (it = this->matches.begin(); it != this->matches.end(); it++)
        delete *it;
}

int
QLocus::add_match(int catalog_id, allele_type cat_type, allele_type query_type, int distance)
{
    Match *m = new Match;

    m->cat_id     = catalog_id;
    m->cat_type   = cat_type;
    m->query_type = query_type;
    m->dist       = distance;

    this->matches.push_back(m);

    return 0;
}

int
QLocus::add_match(int catalog_id, allele_type cat_type, allele_type query_type, int distance, string cigar)
{
    Match *m = new Match;

    m->cat_id     = catalog_id;
    m->cat_type   = cat_type;
    m->query_type = query_type;
    m->dist       = distance;
    m->cigar      = cigar;

    this->matches.push_back(m);

    return 0;
}

int
QLocus::add_match(int catalog_id, allele_type cat_type)
{
    Match *m = new Match;

    m->cat_id     = catalog_id;
    m->cat_type   = cat_type;
    m->query_type = "";
    m->dist       = 0;

    this->matches.push_back(m);

    return 0;
}

int
QLocus::clear_matches()
{
    vector<Match *>::iterator it;

    for (it = this->matches.begin(); it != this->matches.end(); it++)
        delete *it;
    this->matches.clear();

    return 0;
}

void
CLocAlnSet::clear() {
    id_= -1;
    aln_pos_.clear();
    ref_ = DNASeq4();
    mpopi_ = NULL;
    reads_.clear();
    reads_per_sample_.clear();
}


void
CLocAlnSet::merge_paired_reads()
{
    //
    // Sort reads by name. Paired reads should have the same name but end with
    // respectively "/1" and "/2".
    //
    sort(this->reads_.begin(), this->reads_.end(),
         [](const SAlnRead& r1, const SAlnRead& r2) { return r1.name < r2.name; }
         );

    // Merge paired reads.
    for (auto r1 = this->reads_.begin(); r1 != this->reads_.end(); ++r1) {
        auto  r2 = r1;
        ++r2;

        if (r2 == this->reads_.end())
            break;
        
        const string& n1 = r1->name;
        const string& n2 = r2->name;
        const size_t   l = n1.length();

        if (n2.length() == l && l >= 2 &&
            n1[l-2] == '/' && n1[l-1] == '1' &&
            n2[l-2] == '/' && n2[l-1] == '2' && 
            n1.substr(0, l-2) == n2.substr(0, l-2)) {

            // r1 and r2 are paired, merge them.
            assert(r1->sample == r2->sample);
            *r1 = SAlnRead(AlnRead::merger_of(move(*r1), move(*r2)), r1->sample);
            
            assert(cigar_length_ref(r1->cigar) == ref_.length());

            // Mark r2 for removal and skip it.
            r2->seq.clear();
            ++r1;
            if (r1 == this->reads_.end())
                break;
        }
    }

    // Remove emptied reads.
    reads_.erase(std::remove_if(reads_.begin(),
                                reads_.end(),
                                [](const Read& r) { return r.seq.empty(); }
                                ),
                 reads_.end());

    // Refresh `reads_per_sample_`.
    reads_per_sample_ = vector<vector<size_t>>(mpopi().samples().size());
    for (size_t i=0; i<reads_.size(); ++i)
        reads_per_sample_[reads_[i].sample].push_back(i);
}

ostream& operator<< (ostream& os, const CLocAlnSet& loc) {
    os << "ref\t.\t" << loc.ref().str();
    for (auto& r : loc.reads_)
        os << "\n" << r.name << "\t" << loc.mpopi().samples()[r.sample].name << "\t" << r.aln();
    return os;
}

CLocAlnSet
CLocAlnSet::juxtapose(CLocAlnSet&& left, CLocAlnSet&& right, long offset)
{
    assert(offset == +10); // xxx Overlapping contigs aren't actually implemented.

    assert(left.id() == right.id());
    assert(left.pos() == right.pos());
    assert(&left.mpopi() == &right.mpopi());

    CLocAlnSet merged (move(left));
    size_t left_ref_len = merged.ref().length();

    // Extend the reference sequence.
    for (long i=0; i<offset; ++i)
        merged.ref_.push_back(Nt4::n);
    merged.ref_.append(right.ref().begin(), right.ref().end());

    // Extend the left reads.
    for (SAlnRead& r : merged.reads_) {
        assert(offset >= 0);
        cigar_extend_right(r.cigar, offset + right.ref().length());
        assert(cigar_length_ref(r.cigar) == merged.ref().length());
    }

    // Extend & add the right reads.
    for (SAlnRead& r : right.reads_) {
        assert(offset >= 0);
        cigar_extend_left(r.cigar, offset + left_ref_len);
        merged.add(move(r));
    }
    right.reads_ = vector<SAlnRead>();
    right.reads_per_sample_ = vector<vector<size_t>>(right.mpopi().samples().size());

    return merged;
}
