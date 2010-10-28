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

//
// kmers.cc -- routines to generate and hash K-mers
//
// Julian Catchen
// jcatchen@uoregon.edu
// University of Oregon
//
// $Id$
//
#include "kmers.h"

int determine_kmer_length(int read_len, int dist) {
    int kmer_len, span, min_matches;

    //
    // If distance allowed between sequences is 0, then k-mer length equals read length.
    //
    if (dist == 0) 
        return read_len;

    //
    // Longer k-mer lengths will provide a smaller hash, with better key placement.
    // Increase the kmer_len until we start to miss hits at the given distance. Then
    // back the kmer_len off one unit to get the final value.
    //
    for (kmer_len = 5; kmer_len < read_len; kmer_len += 2) {
        span = (kmer_len * (dist + 1)) - 1;

        min_matches = read_len - span;

        if (min_matches <= 0) break;
    }

    if (kmer_len >= read_len) {
        cerr << "Unable to find a suitable k-mer length for matching.\n";
        exit(1);
    }

    kmer_len -= 2;

    cerr << 
        "  Distance allowed between stacks: " << dist << "\n" <<
        "  Using a k-mer length of " << kmer_len << "\n";

    return kmer_len;
}

int calc_min_kmer_matches(int kmer_len, int dist, int read_len, bool exit_err) {
    int span, min_matches;

    span = (kmer_len * (dist + 1)) - 1;

    min_matches = read_len - span;

    cerr << "  Miniumum number of k-mers to define a match: " << min_matches << "\n";

    if (exit_err && min_matches <= 0) {
        cerr << 
            "Combination of k-mer length (" << kmer_len << ") and edit distance (" << dist << ") allows for " <<
            "sequences to be missed by the matching algorithm.\n";
        exit(1);
    }

    return min_matches;
}

int generate_kmers(const char *seq, int kmer_len, int num_kmers, vector<char *> &kmers) {
    char *kmer;
    const char *k = seq;

    for (int i = 0; i < num_kmers; i++) {
        kmer = new char[kmer_len + 1];
        strncpy(kmer, k, kmer_len);
        kmer[kmer_len] = '\0';
        kmers.push_back(kmer);
        k++;
    }

    return 0;
}

int generate_permutations(map<int, char **> &pstrings, int width) {
    int   i, j, rem, div, num;
    char *p;
    // 
    // Given a k-mer that allows wildcards -- 'N' characters, we need to generate all
    // possible k-mers. To do so, we will generate a range of numbers that we convert to
    // base 4, assuming that 0 = 'A', 1 = 'C', 2 = 'G', 3 = 'T'.
    //
    const int base = 4;
    int range      = (int) pow(4, width);

    //
    // Create an array of strings to hold the permuted nucleotides.
    //
    char **strings = new char * [range];
    for (i = 0; i < range; i++)
        strings[i] = new char[width + 1];

    for (i = 0; i < range; i++) {
        for (j = 0; j < width; j++) 
            strings[i][j] = 'A';
        strings[i][width] = '\0';
    }

    for (i = 0; i < range; i++) {
        //
        // Convert this number to base 4
        //
        p   = strings[i]; p += width - 1;
        num = i;
        do {
            div = (int) floor(num / base);
            rem = num % base;

            switch(rem) {
            case 0:
                *p = 'A';
                break;
            case 1:
                *p = 'C';
                break;
            case 2:
                *p = 'G';
                break;
            case 3:
                *p = 'T';
                break;
            }
            num = div;
            p--;
        } while (div > 0);
    }

    pstrings[width] = strings;

    return 0;
}

int populate_kmer_hash(map<int, MergedStack *> &merged, KmerHashMap &kmer_map, int kmer_len, bool exclude_masked_tags) {
    map<int, MergedStack *>::iterator it;
    MergedStack *tag;
    int j;
    //
    // Break each stack down into k-mers and create a hash map of those k-mers
    // recording in which sequences they occur.
    //
    int num_kmers = strlen(merged.begin()->second->con) - kmer_len + 1;

    for (it = merged.begin(); it != merged.end(); it++) {
        tag = it->second;

        // Don't compute distances for masked tags
        if (exclude_masked_tags && tag->masked) continue;

        generate_kmers(tag->con, kmer_len, num_kmers, tag->kmers);

        // Hash the kmers
        for (j = 0; j < num_kmers; j++)
            kmer_map[tag->kmers[j]].push_back(tag->id);
    }

    //dump_kmer_map(kmer_map);

    return 0;
}

int populate_kmer_hash(map<int, Locus *> &catalog, CatKmerHashMap &kmer_map, int kmer_len) {
    map<int, Locus *>::iterator it;
    vector<pair<allele_type, string> >::iterator allele;
    vector<char *> kmers;
    Locus         *tag;
    char          *hash_key;
    bool           exists;
    int            j;

    //
    // Break each stack down into k-mers and create a hash map of those k-mers
    // recording in which sequences they occur.
    //
    int num_kmers = strlen(catalog.begin()->second->con) - kmer_len + 1;

    for (it = catalog.begin(); it != catalog.end(); it++) {
        tag = it->second;

        //
        // Iterate through the possible Catalog alleles
        //
        for (allele = tag->strings.begin(); allele != tag->strings.end(); allele++) {
            //
            // Generate and hash the kmers for this allele string
            //
            generate_kmers(allele->second.c_str(), kmer_len, num_kmers, kmers);

            for (j = 0; j < num_kmers; j++) {

                exists = kmer_map.count(kmers[j]) == 0 ? false : true;

                if (exists) {
                    hash_key = kmers[j];
                } else {
                    hash_key = new char [strlen(kmers[j]) + 1];
                    strcpy(hash_key, kmers[j]);
                }

                kmer_map[hash_key].push_back(make_pair(allele->first, tag->id));
            }

            for (j = 0; j < num_kmers; j++)
                delete [] kmers[j];
            kmers.clear();
        }
    }

    //dump_kmer_map(kmer_map);

    return 0;
}

int dist(const char *tag_1, Locus *tag_2, allele_type allele) {
    int   dist = 0;
    const char *p    = tag_1;
    const char *q    = NULL;
    const char *end  = p + strlen(p);

    //
    // Identify which matching string has the proper allele
    //
    vector<pair<allele_type, string> >::iterator it;

    for (it = tag_2->strings.begin(); it != tag_2->strings.end(); it++)
        if (it->first == allele) 
            q = it->second.c_str();
    if (q == NULL) return -1;

    // Count the number of characters that are different
    // between the two sequences.
    while (p < end) {
	dist += (*p == *q) ? 0 : 1;
	p++; 
	q++;
    }

    return dist;
}

int dist(Locus *tag_1, Locus *tag_2) {
    int   dist = 0;
    char *p    = tag_1->con;
    char *q    = tag_2->con;
    char *end  = p + strlen(p);

    // Count the number of characters that are different
    // between the two sequences.
    while (p < end) {
	dist += (*p == *q) ? 0 : 1;
	p++; 
	q++;
    }

    return dist;
}

int dist(MergedStack *tag_1, MergedStack *tag_2) {
    int   dist = 0;
    char *p    = tag_1->con;
    char *q    = tag_2->con;
    char *end  = p + strlen(p);

    // Count the number of characters that are different
    // between the two sequences.
    while (p < end) {
	dist += (*p == *q) ? 0 : 1;
	p++; 
	q++;
    }

    return dist;
}

int dist(MergedStack *tag_1, Seq *rem) {
    int   dist = 0;
    char *p    = tag_1->con;
    char *q    = rem->seq;
    char *end  = p + strlen(p);

    // Count the number of characters that are different
    // between the two sequences.
    while (p < end) {
	dist += (*p == *q) ? 0 : 1;
	p++; 
	q++;
    }

    return dist;
}

bool compare_dist(pair<int, int> a, pair<int, int> b) {
    return (a.second < b.second);
}

int dump_kmer_map(KmerHashMap &kmer_map) {
    KmerHashMap::iterator kit;
    vector<int>::iterator vit;

    cerr << kmer_map.size() << " keys in the map.\n";

    int i = 1;
    for (kit = kmer_map.begin(); kit != kmer_map.end(); kit++) {
        cerr << "Key #" << i << " " << kit->first << ": ";
        for (vit = (kit->second).begin(); vit != (kit->second).end(); vit++) 
            cerr << " " << *vit;
        cerr << "\n";
        i++;

        if (i > 1000) break;
    }

    return 0;
}
