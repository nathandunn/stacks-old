// -*-mode:c++; c-style:k&r; c-basic-offset:4;-*-
//
// Copyright 2016, Julian Catchen <jcatchen@illinois.edu>
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
// aln_utils.cc -- common routines needed for dealing with gapped alignments.
//

#include "aln_utils.h"

string
invert_cigar(string cigar)
{
    for (uint i = 0; i < cigar.length(); i++) {
        if (cigar[i] == 'I')
            cigar[i] = 'D';
        else if (cigar[i] == 'D')
            cigar[i] = 'I';
    }

    return cigar;
}

int 
parse_cigar(const char *cigar_str, vector<pair<char, uint> > &cigar)
{
    char buf[id_len];
    int  dist;
    const char *p, *q;

    uint seqlen = 0;

    cigar.clear();

    p = cigar_str;

    while (*p != '\0') {
        q = p + 1;

        while (*q != '\0' && isdigit(*q))
            q++;
        strncpy(buf, p, q - p);
        buf[q-p] = '\0';
        dist = is_integer(buf);

        cigar.push_back(make_pair(*q, dist));

        seqlen += dist;

        p = q + 1;
    }

    return seqlen;
}

string
apply_cigar_to_seq(const char *seq, vector<pair<char, uint> > &cigar)
{
    uint   size = cigar.size();
    char   op;
    uint   dist, bp, edited_bp, stop;
    string edited_seq;

    //
    // Calculate the overall sequence length.
    //
    uint seqlen = 0;
    for (uint i = 0; i < size; i++)
        seqlen += cigar[i].second;

    bp  = 0;

    edited_seq.reserve(seqlen);

    for (uint i = 0; i < size; i++)  {
        op   = cigar[i].first;
        dist = cigar[i].second;

        switch(op) {
        case 'S':
            stop = bp + dist;
            while (bp < stop) {
                edited_seq.push_back('N');
                bp++;
            }
            break;
        case 'D':
	    edited_bp = 0;
            while (edited_bp < dist) {
                edited_seq.push_back('N');
		edited_bp++;
            }
            break;
        case 'I':
        case 'M':
            stop = bp + dist;
            while (bp < stop) {
                edited_seq.push_back(seq[bp]);
                bp++;
            }
            break;
        default:
            break;
        }
    }

    return edited_seq;
}

string
apply_cigar_to_model_seq(const char *seq, vector<pair<char, uint> > &cigar)
{
    uint   size = cigar.size();
    char   op;
    uint   dist, bp, edited_bp, stop;
    string edited_seq;

    //
    // Calculate the overall sequence length.
    //
    uint seqlen = 0;
    for (uint i = 0; i < size; i++)
        seqlen += cigar[i].second;

    bp  = 0;

    edited_seq.reserve(seqlen);

    for (uint i = 0; i < size; i++)  {
        op   = cigar[i].first;
        dist = cigar[i].second;

        switch(op) {
        case 'S':
            stop = bp + dist;
            while (bp < stop) {
                edited_seq.push_back('U');
                bp++;
            }
            break;
        case 'D':
	    edited_bp = 0;
            while (edited_bp < dist) {
                edited_seq.push_back('U');
		edited_bp++;
            }
            break;
        case 'I':
        case 'M':
            stop = bp + dist;
            while (bp < stop) {
                edited_seq.push_back(seq[bp]);
                bp++;
            }
            break;
        default:
            break;
        }
    }

    return edited_seq;
}

int
apply_cigar_to_seq(char *seq, uint seq_len, const char *old_seq, vector<pair<char, uint> > &cigar)
{
    uint   size = cigar.size();
    char   op;
    uint   dist, bp, seq_bp, stop;

    bp         = 0;
    seq_bp     = 0;

    for (uint i = 0; i < size; i++)  {
        op   = cigar[i].first;
        dist = cigar[i].second;

        switch(op) {
        case 'S':
            stop = seq_bp + dist;
            stop = stop > seq_len ? seq_len : stop;
            while (seq_bp < stop) {
                seq[seq_bp] = 'N';
                seq_bp++;
                bp++;
            }
            break;
        case 'D':
            stop = seq_bp + dist;
            stop = stop > seq_len ? seq_len : stop;
            while (seq_bp < stop) {
                seq[seq_bp] = 'N';
		seq_bp++;
            }
            break;
        case 'I':
        case 'M':
            stop = bp + dist;
            stop = stop > seq_len ? seq_len : stop;
            while (bp < stop) {
                seq[seq_bp] = old_seq[bp];
                seq_bp++;
                bp++;
            }
            break;
        default:
            break;
        }
    }

    seq[seq_len] = '\0';

    return 0;
}

int
apply_cigar_to_model_seq(char *seq, uint seq_len, const char *model, vector<pair<char, uint> > &cigar)
{
    uint   size = cigar.size();
    char   op;
    uint   dist, model_bp, seq_bp, stop;

    model_bp  = 0;
    seq_bp    = 0;

    for (uint i = 0; i < size; i++)  {
        op   = cigar[i].first;
        dist = cigar[i].second;

        switch(op) {
        case 'S':
            stop = seq_bp + dist;
            stop = stop > seq_len ? seq_len : stop;
            while (seq_bp < stop) {
                seq[seq_bp] = 'U';
                seq_bp++;
                model_bp++;
            }
            break;
        case 'D':
            stop = seq_bp + dist;
            stop = stop > seq_len ? seq_len : stop;
            while (seq_bp < stop) {
                seq[seq_bp] = 'U';
		seq_bp++;
            }
            break;
        case 'I':
        case 'M':
            stop = model_bp + dist;
            stop = stop > seq_len ? seq_len : stop;
            while (model_bp < stop) {
                seq[seq_bp] = model[model_bp];
                seq_bp++;
                model_bp++;
            }
            break;
        default:
            break;
        }
    }

    seq[seq_len] = '\0';

    return 0;
}

string
remove_cigar_from_seq(const char *seq, vector<pair<char, uint> > &cigar)
{
    uint   size = cigar.size();
    char   op;
    uint   dist, bp, edited_bp, stop;
    string edited_seq;

    //
    // Calculate the overall sequence length.
    //
    uint seqlen = 0;
    for (uint i = 0; i < size; i++)
        seqlen += cigar[i].first != 'D' ? cigar[i].second : 0; 

    bp  = 0;

    edited_seq.reserve(seqlen);

    for (uint i = 0; i < size; i++)  {
        op   = cigar[i].first;
        dist = cigar[i].second;

        switch(op) {
        case 'D':
	    edited_bp = 0;
            while (edited_bp < dist) {
                edited_bp++;
		bp++;
            }
            break;
        case 'I':
        case 'M':
            stop = bp + dist;
            while (bp < stop) {
                edited_seq.push_back(seq[bp]);
                bp++;
            }
            break;
        default:
            break;
        }
    }

    return edited_seq;
}

int
adjust_snps_for_gaps(vector<pair<char, uint> > &cigar, Locus *loc)
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
adjust_and_add_snps_for_gaps(vector<pair<char, uint> > &cigar, Locus *loc)
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
remove_snps_from_gaps(vector<pair<char, uint> > &cigar, Locus *loc)
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
