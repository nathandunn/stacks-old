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

// shell:
// echo "
// print(', '.join([ ('true' if chr(i) in set('0123456789=DHIMNPSX') else 'false') for i in range(256) ]))
// " | python3 | fold -s
const bool is_cigar_char[256] = {
    false, false, false, false, false, false, false, false, false, false, false,
    false, false, false, false, false, false, false, false, false, false, false,
    false, false, false, false, false, false, false, false, false, false, false,
    false, false, false, false, false, false, false, false, false, false, false,
    false, false, false, false, true, true, true, true, true, true, true, true,
    true, true, false, false, false, true, false, false, false, false, false,
    false, true, false, false, false, true, true, false, false, false, true, true,
    false, true, false, false, true, false, false, false, false, true, false,
    false, false, false, false, false, false, false, false, false, false, false,
    false, false, false, false, false, false, false, false, false, false, false,
    false, false, false, false, false, false, false, false, false, false, false,
    false, false, false, false, false, false, false, false, false, false, false,
    false, false, false, false, false, false, false, false, false, false, false,
    false, false, false, false, false, false, false, false, false, false, false,
    false, false, false, false, false, false, false, false, false, false, false,
    false, false, false, false, false, false, false, false, false, false, false,
    false, false, false, false, false, false, false, false, false, false, false,
    false, false, false, false, false, false, false, false, false, false, false,
    false, false, false, false, false, false, false, false, false, false, false,
    false, false, false, false, false, false, false, false, false, false, false,
    false, false, false, false, false, false, false, false, false, false, false,
    false, false, false, false, false, false, false, false, false, false, false,
    false, false, false, false, false, false, false, false, false, false, false,
    false
};

ostream& operator<< (ostream& os, const Cigar& cig) {
    for (auto op : cig)
        os << op.second << (op.first == '\0' ? '?' : op.first);
    return os;
}

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
parse_cigar(const char *cigar_str, Cigar &cigar, bool check_correctness)
{
    cigar.clear();
    const char* p      = cigar_str;
    uint        seqlen = 0;

    while(*p != '\0') {
        char* q;
        uint len = strtol(p, &q, 10);
        char c = *q;

        if (check_correctness) {
            if (q == p || c == '\0') {
                // No number or no qualifier, respectively.
                cerr << "Error: Malformed CIGAR string '" << cigar_str << "'.\n";
                throw exception();
            }
            switch (c) {
            case 'M':
            case '=':
            case 'X':
            case 'I':
            case 'D':
            case 'S':
            case 'N':
            case 'H':
            case 'P':
                break;
            default:
                cerr << "Warning: Unknown CIGAR operation '" << c << "' in '" << cigar_str << "'.\n";
                break;
            }
        }

        cigar.push_back({c, len});
        seqlen += len;
        p = q + 1;
    }

    return seqlen;
}

string
apply_cigar_to_seq(const char *seq, Cigar &cigar)
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
apply_cigar_to_model_seq(const char *seq, Cigar &cigar)
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
apply_cigar_to_seq(char *seq, uint seq_len, const char *old_seq, Cigar &cigar)
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
apply_cigar_to_model_seq(char *seq, uint seq_len, const char *model, Cigar &cigar)
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
remove_cigar_from_seq(const char *seq, Cigar &cigar)
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

std::tuple<uint,uint,uint> cigar_lengths(const Cigar& cigar) {
    uint padded_len = 0;
    uint ref_len = 0;
    uint seq_len = 0;
    for (auto& op : cigar) {
        padded_len += op.second;
        switch (op.first) {
        case 'M':
        case '=':
        case 'X':
            // Consume both ref & seq.
            ref_len += op.second;
            seq_len += op.second;
            break;
        case 'I':
            // Consume seq.
            seq_len += op.second;
            break;
        case 'D':
        case 'S':
        case 'N':
        case 'H':
        case 'P':
            // Consume ref.
            ref_len += op.second;
            break;
        default:
            DOES_NOT_HAPPEN;
            break;
        }
    }

    return std::make_tuple(padded_len, ref_len, seq_len);
}

void simplify_cigar_to_MDI(Cigar& cig) {
    if (cig.empty())
        return;

    // Replace operations with the relevant equivalent in "MDI".
    for (auto& op : cig) {
        switch (op.first) {
        case '=':
        case 'X': op.first = 'M'; break;
        case 'S': op.first = 'I'; break;
        case 'N':
        case 'H': op.first = 'D'; break;
        default: break;
        case 'P': op.first = '\0'; break;
        }
    }

    // Collapse identical successive operations.
    auto prev = cig.rbegin();
    auto op = ++cig.rbegin(); // n.b. `cig` isn't empty.
    while(op != cig.rend()) {
        if (op->first == prev->first) {
            op->second += prev->second;
            prev->first = '\0';
        }
        ++prev;
        ++op;
    }

    // Remove '\0' operations.
    cig.erase(std::remove_if(
            cig.begin(), cig.end(),
            [](const pair<char,uint>& op){return op.first == '\0';}
            ), cig.end());
}
