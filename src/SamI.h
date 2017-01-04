// -*-mode:c++; c-style:k&r; c-basic-offset:4;-*-
//
// Copyright 2010-2016, Julian Catchen <jcatchen@illinois.edu>
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

#ifndef __SAMI_H__
#define __SAMI_H__

//
// Code to parse Sam format. This format is created for
// reads that have been aligned to a reference genome. It takes the tab-separated form:
//
// <query> <strand> <chromosome> <base pair> ... <sequence> <phred quality score> ...
//
// One record per line.
//
#include "stacks.h"
#include "input.h"
#include "BamI.h"

class Sam: public Input {
    int parse_cigar(const char *, vector<pair<char, uint> > &, strand_type);

 public:
    Sam(const char *path);
    ~Sam() {}
    Seq *next_seq();
    int  next_seq(Seq& s);
};

Sam::Sam(const char *path)
: Input(path)
{

}

Seq* Sam::next_seq() {
    Seq* s = new Seq();
    if(next_seq(*s) != 1) {
        delete s;
        s = NULL;
    }
    return s;
}

int Sam::next_seq(Seq& s) {
    vector<string> parts;
    int  flag  = 0;
    uint len;

    //
    // Read a record from the file and place it in a Seq object, skipping header
    // definitions and unaligned sequences.
    //
    do {
        fh.getline(line, max_len);

        if (!fh.good())
            return 0;

        len = strlen(line);
        if (line[len - 1] == '\r') line[len - 1] = '\0';

        parse_tsv(line, parts);

        //
        // According to SAM spec FLAGs are the second field,
        // if FLAG bit 0x4 is set, sequence is not mapped.
        //
        flag = atoi(parts[1].c_str());

    } while (parts[0][0] == '@');

    //
    // Parse the type of the record.
    //
    AlnT aln_type;
    if (flag & BAM_FUNMAP)
        aln_type = AlnT::null;
    else if (flag & BAM_FSECONDARY)
        aln_type = AlnT::secondary;
    else if (flag & BAM_FSUPPLEMENTARY)
        aln_type = AlnT::supplementary;
    else
        aln_type = AlnT::primary;

    if (aln_type == AlnT::null) {
        s = Seq(parts[0].c_str(), parts[9].c_str(), parts[10].c_str()); // Read ID, Sequence, Quality
    } else {
        //
        // Check which strand this is aligned to:
        //
        strand_type strand = flag & BAM_FREVERSE ? strand_minus : strand_plus;

        //
        // Parse the alignment CIGAR string.
        //
        vector<pair<char, uint> > cigar;
        parse_cigar(parts[5].c_str(), cigar, strand);

        //
        // If the read was aligned on the reverse strand (and is therefore reverse complemented)
        // alter the start point of the alignment to reflect the right-side of the read, at the
        // end of the RAD cut site.
        //
        int bp = bam_find_start_bp(atoi(parts[3].c_str()), strand, cigar);
        bp--; // SAM uses 1-based genome positions.

        //
        // Calculate the percentage of the sequence that was aligned to the reference.
        //
        len = 0;
        for (uint i = 0; i < cigar.size(); i++)
            switch (cigar[i].first) {
            case 'M':
            case 'I':
            case '=':
            case 'X':
                len += cigar[i].second;
                break;
            case 'D':
            case 'S':
                break;
            case 'H':
            case 'N':
                static bool emitted_hn_warning = false;
                if(aln_type == AlnT::primary && !emitted_hn_warning) {
                    cerr << "Warning: Some CIGARs contained H and/or N operations.\n";
                    emitted_hn_warning = true;
                }
                break;
            default:
                cerr << "Error parsing CIGAR string '" << cigar[i].second << cigar[i].first << "'.\n";
                break;
            }
        double pct_aln = double(len) / double(parts[9].length());

        s = Seq(parts[0].c_str(), parts[9].c_str(), parts[10].c_str(), // Read ID, Sequence, Quality
                parts[2].c_str(), bp, strand, aln_type, pct_aln); // Chromosome, etc.

        if (cigar.size() > 0)
            bam_edit_gaps(cigar, s.seq);
    }
    return true;
}

int
Sam::parse_cigar(const char *cigar_str, vector<pair<char, uint> > &cigar, strand_type strand)
{
    char buf[id_len];
    int  dist;
    const char *p, *q;

    p = cigar_str;

    if (*p == '*') return 0;

    while (*p != '\0') {
        q = p + 1;

        while (*q != '\0' && isdigit(*q))
            q++;
        strncpy(buf, p, q - p);
        buf[q-p] = '\0';
        dist = atoi(buf);

        //
        // If aligned to the negative strand, sequence has been reverse complemented and
        // CIGAR string should be interpreted in reverse.
        //
        if (strand == strand_plus)
            cigar.push_back(make_pair(*q, dist));
        else
            cigar.insert(cigar.begin(), make_pair(*q, dist));

        p = q + 1;
    }

    return 0;
}

#endif // __SAMI_H__
