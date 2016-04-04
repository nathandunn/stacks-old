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

//
// ustacks -- build denovo stacks
//

#include "ustacks.h"

//
// Global variables to hold command-line options.
//
FileT   in_file_type;
string  in_file;
string  out_path;
int     num_threads       = 1;
int     sql_id            = 0;
bool    call_sec_hapl     = true;
bool    set_kmer_len      = true;
int     kmer_len          = 0;
int     max_kmer_len      = 19;
int     min_merge_cov     = 3;
uint    max_subgraph      = 3;
int     dump_graph        = 0;
int     retain_rem_reads  = false;
int     deleverage_stacks = 0;
int     remove_rep_stacks = 0;
int     max_utag_dist     = 2;
int     max_rem_dist      = -1;
bool    gapped_alignments = false;
double  min_match_len     = 0.80;
double  max_gaps          = 2;
int     deleverage_trigger;
int     removal_trigger;
//
// For use with the multinomial model to call fixed nucleotides.
//
modelt model_type         = snp;
double alpha              = 0.05;
double bound_low          = 0.0;
double bound_high         = 1.0;
double p_freq             = 0.5;
double barcode_err_freq   = 0.0;
double heterozygote_limit = -3.84;
double homozygote_limit   =  3.84;

int main (int argc, char* argv[]) {

    parse_command_line(argc, argv);

    //
    // Set the max remainder distance to be greater than the max_utag_dist, if it is not
    // specified on the command line.
    //
    if (max_rem_dist == -1) max_rem_dist = max_utag_dist + 2;

    cerr << "ustacks paramters selected:\n"
         << "  Min depth of coverage to create a stack: " << min_merge_cov << "\n"
         << "  Max distance allowed between stacks: " << max_utag_dist << "\n"
         << "  Max distance allowed to align secondary reads: " << max_rem_dist << "\n"
         << "  Max number of stacks allowed per de novo locus: " << max_subgraph << "\n"
         << "  Deleveraging algorithm: " << (deleverage_stacks ? "enabled" : "disabled") << "\n"
         << "  Removal algorithm: " << (remove_rep_stacks ? "enabled" : "disabled") << "\n"
         << "  Model type: "; 
    switch (model_type) {
    case snp:
        cerr << "SNP\n";
        break;
    case fixed: 
        cerr << "Fixed\n";
        break;
    case bounded:
        cerr << "Bounded; lower epsilon bound: " << bound_low << "; upper bound: " << bound_high << "\n";
        break;
    }
    cerr << "  Alpha significance level for model: " << alpha << "\n"
         << "  Gapped alignments: " << (gapped_alignments ? "enabled" : "disabled") << "\n";

    //
    // Set limits to call het or homozygote according to chi-square distribution with one 
    // degree of freedom:
    //   http://en.wikipedia.org/wiki/Chi-squared_distribution#Table_of_.CF.872_value_vs_p-value
    //
    if (alpha == 0.1) {
        heterozygote_limit = -2.71;
        homozygote_limit   =  2.71;
    } else if (alpha == 0.05) {
        heterozygote_limit = -3.84;
        homozygote_limit   =  3.84;
    } else if (alpha == 0.01) {
        heterozygote_limit = -6.64;
        homozygote_limit   =  6.64;
    } else if (alpha == 0.001) {
        heterozygote_limit = -10.83;
        homozygote_limit   =  10.83;
    }

    //
    // Set the number of OpenMP parallel threads to execute.
    //
    #ifdef _OPENMP
    omp_set_num_threads(num_threads);
    #endif

    DNASeqHashMap     radtags;
    vector<DNANSeq *> radtags_keys;
    map<int, Rem *>   remainders;
    set<int>          merge_map;
    map<int, Stack *> unique;

    load_radtags(in_file, radtags, radtags_keys);

    reduce_radtags(radtags, unique, remainders);

    free_radtags_hash(radtags, radtags_keys);

    // dump_unique_tags(unique);

    double cov_mean, cov_stdev, cov_max;
    
    calc_coverage_distribution(unique, cov_mean, cov_stdev, cov_max);
    cerr << "Initial coverage mean: " << cov_mean << "; Std Dev: " << cov_stdev << "; Max: " << cov_max << "\n";

    calc_triggers(cov_mean, cov_stdev, 1, deleverage_trigger, removal_trigger);

    cerr << "Deleveraging trigger: " << deleverage_trigger << "; Removal trigger: " << removal_trigger << "\n";

    map<int, MergedStack *> merged;

    populate_merged_tags(unique, merged);

    cerr << merged.size() << " initial stacks were populated; " << remainders.size() << " stacks were set aside as secondary reads.\n";

    if (remove_rep_stacks) {
        cerr << "Calculating distance for removing repetitive stacks.\n";
        calc_kmer_distance(merged, 1);
        cerr << "Removing repetitive stacks.\n";
        remove_repetitive_stacks(unique, merged);
    }

    calc_coverage_distribution(unique, merged, cov_mean, cov_stdev, cov_max);
    cerr << "Post-Repeat Removal, coverage depth Mean: " << cov_mean << "; Std Dev: " << cov_stdev << "; Max: " << cov_max << "\n";

    cerr << "Calculating distance between stacks...\n";
    calc_kmer_distance(merged, max_utag_dist);

    cerr << "Merging stacks, maximum allowed distance: " << max_utag_dist << " nucleotide(s)\n";
    merge_stacks(unique, remainders, merged, merge_map, max_utag_dist);

    call_consensus(merged, unique, remainders, false);

    calc_coverage_distribution(unique, merged, cov_mean, cov_stdev, cov_max);
    cerr << "After merging, coverage depth Mean: " << cov_mean << "; Std Dev: " << cov_stdev << "; Max: " << cov_max << "\n";

    //dump_merged_tags(merged);

    cerr << "Merging remainder radtags\n";
    merge_remainders(merged, remainders);

    calc_coverage_distribution(unique, remainders, merged, cov_mean, cov_stdev, cov_max);
    cerr << "After remainders merged, coverage depth Mean: " << cov_mean << "; Std Dev: " << cov_stdev << "; Max: " << cov_max << "\n";

    if (gapped_alignments) {
        call_consensus(merged, unique, remainders, false);

        cerr << "Searching for gaps between merged stacks...\n";
        search_for_gaps(merged, min_match_len);

        merge_gapped_alns(unique, remainders, merged);

        calc_coverage_distribution(unique, remainders, merged, cov_mean, cov_stdev, cov_max);
        cerr << "After gapped alignments, coverage depth Mean: " << cov_mean << "; Std Dev: " << cov_stdev << "; Max: " << cov_max << "\n";
    }
    
    //
    // Call the final consensus sequence and invoke the SNP model.
    //
    cerr << "Calling final consensus sequences, invoking SNP-calling model...\n";
    call_consensus(merged, unique, remainders, true);

    count_raw_reads(unique, remainders, merged);

    cerr << "Writing loci, SNPs, and alleles to '" << out_path << "'...\n";
    write_results(merged, unique, remainders);
    cerr << "done.\n";

    return 0;
}

int
merge_gapped_alns(map<int, Stack *> &unique, map<int, Rem *> &rem, map<int, MergedStack *> &merged)
{
    map<int, MergedStack *> new_merged;
    map<int, MergedStack *>::iterator it;
    MergedStack *tag_1, *tag_2, *merged_tag;

    int  id        = 1;
    uint merge_cnt = 0;

    set<int> processed;
    string   cigar_1, cigar_2;

    for (it = merged.begin(); it != merged.end(); it++) {
        if (processed.count(it->first) > 0)
            continue;

        tag_1 = it->second;
        //
        // No gapped alignments for this stack, or this stack has already been set aside.
        //
        if (tag_1->masked || tag_1->alns.size() != 1)
            continue;

        //
        // Found a gapped alignment. Make sure the alignments are the same.
        //
        tag_2   = merged[tag_1->alns[0].first];
        cigar_1 = tag_1->alns[0].second;
        cigar_2 = tag_2->alns.size() != 1 ? "" : tag_2->alns[0].second;

        for (uint i = 0; i < cigar_1.length(); i++) {
            if (cigar_1[i] == 'I')
                cigar_1[i] = 'D';
            else if (cigar_1[i] == 'D')
                cigar_1[i] = 'I';
        }

        if (cigar_1 == cigar_2) {
            //
            // Edit the sequences to accommodate any added gaps.
            //
            vector<pair<char, uint> > cigar;
            // cerr << "CIGAR to parse: " << tag_1->alns[0].second << "\n";
            parse_cigar(tag_1->alns[0].second.c_str(), cigar);

            uint   gap_cnt = 0;
            double aln_len = 0;
            char   op;
            for (uint j = 0; j < cigar.size(); j++) {
                op = cigar[j].first;
                switch (op) {
                case 'I':
                case 'D':
                    gap_cnt++;
                    break;
                case 'M':
                    aln_len += cigar[j].second;
                    break;
                }
            }

            //
            // Check that the alignment still contains fewer than 
            // max_utag_dist mismatches.
            //
            if (dist(tag_1, tag_2, cigar) > max_utag_dist)
                continue;
            //
            // If the alignment has too many gaps, skip it.
            //
            if (gap_cnt > (max_gaps + 1))
                continue;
            //
            // If the alignment doesn't span enough of the two sequences, skip it.
            //
            if ((aln_len / (double) tag_1->len) < min_match_len)
                continue;

            edit_gapped_seqs(unique, rem, tag_1, cigar);

            // cerr << "CIGAR to parse: " << tag_2->alns[0].second << "\n";
            cigar.clear();
            parse_cigar(tag_2->alns[0].second.c_str(), cigar);
            edit_gapped_seqs(unique, rem, tag_2, cigar);

            //
            // Merge the tags.
            //
            merged_tag     = merge_tags(tag_1, tag_2, id);
            new_merged[id] = merged_tag;
            id++;

            //
            // Record the gaps.
            //
            uint pos = 0;
            for (uint j = 0; j < cigar.size(); j++) {
                if (cigar[j].first == 'I' || cigar[j].first == 'D')
                    merged_tag->gaps.push_back(Gap(pos, pos + cigar[j].second));
                pos += cigar[j].second;
            }

            processed.insert(tag_1->id);
            processed.insert(tag_2->id);

            merge_cnt++;
        }
    }

    set<int> merge_set;
    for (it = merged.begin(); it != merged.end(); it++) {
        if (processed.count(it->first))
            continue;
        tag_1          = it->second;
        merge_set.insert(tag_1->id);
        tag_2          = merge_tags(merged, merge_set, id);
        new_merged[id] = tag_2;
        merge_set.clear();
        id++;
    }

    uint new_cnt = new_merged.size();
    uint old_cnt = merged.size();

    //
    // Free the memory from the old map of merged tags.
    //
    for (it = merged.begin(); it != merged.end(); it++)
        delete it->second;

    merged = new_merged;

    cerr << "  " << old_cnt << " stacks merged into " << new_cnt 
         << " stacks; merged " << merge_cnt 
         << " gapped alignments.\n";

    return 0;
}

int
edit_gapped_seqs(map<int, Stack *> &unique, map<int, Rem *> &rem, MergedStack *tag, vector<pair<char, uint> > &cigar)
{
    int    stack_id;
    Stack *s;
    Rem   *r;
    char  *buf = new char[tag->len + 1];

    for (uint i = 0; i < tag->utags.size(); i++) {
        stack_id = tag->utags[i];
        s = unique[stack_id];

        buf = s->seq->seq(buf);
        edit_gaps(cigar, buf);

        delete s->seq;
        s->seq = new DNANSeq(tag->len, buf);
    }

    for (uint i = 0; i < tag->remtags.size(); i++) {
        stack_id = tag->remtags[i];
        r = rem[stack_id];

        buf = r->seq->seq(buf);
        edit_gaps(cigar, buf);

        delete r->seq;
        r->seq = new DNANSeq(tag->len, buf);
    }

    delete [] buf;

    return 0;
}

int 
parse_cigar(const char *cigar_str, vector<pair<char, uint> > &cigar)
{
    char buf[id_len];
    int  dist;
    const char *p, *q;

    p = cigar_str;

    while (*p != '\0') {
        q = p + 1;

        while (*q != '\0' && isdigit(*q))
            q++;
        strncpy(buf, p, q - p);
        buf[q-p] = '\0';
        dist = atoi(buf);

        cigar.push_back(make_pair(*q, dist));

        p = q + 1;
    }

    return 0;
}

int 
edit_gaps(vector<pair<char, uint> > &cigar, char *seq)
{
    char *buf;
    uint  size = cigar.size();
    char  op;
    uint  dist, bp, len, buf_len, buf_size, j, k, stop;

    len = strlen(seq);
    bp  = 0;

    buf      = new char[len + 1];
    buf_size = len + 1;

    for (uint i = 0; i < size; i++)  {
        op   = cigar[i].first;
        dist = cigar[i].second;

        switch(op) {
        case 'S':
            stop = bp + dist;
            stop = stop > len ? len : stop;
            while (bp < stop) {
                seq[bp] = 'N';
                bp++;
            }
            break;
        case 'D':
            //
            // A deletion has occured in the read relative to the reference genome.
            // Pad the read with sufficent Ns to match the deletion, shifting the existing
            // sequence down. Trim the final length to keep the read length consistent.
            //
            k = bp >= len ? len : bp;
            
            strncpy(buf, seq + k, buf_size - 1);
            buf[buf_size - 1] = '\0';
            buf_len         = strlen(buf);

            stop = bp + dist;
            stop = stop > len ? len : stop;
            while (bp < stop) {
                seq[bp] = 'N';
                bp++;
            }

            j = bp;
            k = 0;
            while (j < len && k < buf_len) {
                seq[j] = buf[k];
                k++;
                j++;
            }
            break;
        case 'I':
        case 'M':
            bp += dist;
            break;
        default:
            break;
        }
    }

    delete [] buf;

    return 0;
}

int
dist(MergedStack *tag_1, MergedStack *tag_2, vector<pair<char, uint> > &cigar)
{
    uint  size = cigar.size();
    char  op;
    uint  dist, len, pos_1, pos_2, stop;
    int   mismatches = 0;

    len   = strlen(tag_1->con);
    pos_1 = 0;
    pos_2 = 0;

    for (uint i = 0; i < size; i++)  {
        op   = cigar[i].first;
        dist = cigar[i].second;

        switch(op) {
        case 'D':
            //
            // A deletion has occured in tag_1 relative to tag_2.
            //
            pos_2 += dist;
            break;
        case 'I':
            //
            // An insertion has occured in tag_1 relative to tag_2.
            //
            pos_1 += dist;
            break;
        case 'M':
            stop = pos_1 + dist;
            while (pos_1 < stop && pos_1 < len && pos_2 < len) {
                if (tag_1->con[pos_1] != tag_2->con[pos_2])
                    mismatches++;
                pos_1++;
                pos_2++;
            }
            break;
        default:
            break;
        }
    }

    return mismatches;
}

int
search_for_gaps(map<int, MergedStack *> &merged, double min_match_len)
{
    //
    // Search for loci that can be merged with a gapped alignment.
    //
    KmerHashMap    kmer_map;
    vector<char *> kmer_map_keys;
    MergedStack   *tag_1, *tag_2;
    map<int, MergedStack *>::iterator it;

    //
    // OpenMP can't parallelize random access iterators, so we convert
    // our map to a vector of integer keys.
    //
    vector<int> keys;
    for (it = merged.begin(); it != merged.end(); it++) 
        keys.push_back(it->first);

    //
    // Calculate the number of k-mers we will generate. If kmer_len == 0,
    // determine the optimal length for k-mers.
    //
    int con_len   = strlen(merged[keys[0]]->con);
    int kmer_len  = 19;
    int num_kmers = con_len - kmer_len + 1;

    //
    // Calculate the minimum number of matching k-mers required for a possible sequence match.
    //
    int min_hits = (round((double) con_len * min_match_len) - (kmer_len * max_gaps)) - kmer_len + 1;

    cerr << "  Searching with a k-mer length of " << kmer_len << " (" << num_kmers << " k-mers per read); " << min_hits << " k-mer hits required.\n";

    populate_kmer_hash(merged, kmer_map, kmer_map_keys, kmer_len);
 
    #pragma omp parallel private(tag_1, tag_2)
    {
	GappedAln *aln = new GappedAln(con_len);
        
        #pragma omp for schedule(dynamic) 
        for (uint i = 0; i < keys.size(); i++) {
            tag_1 = merged[keys[i]];

            // Don't compute distances for masked tags
            if (tag_1->masked) continue;

            vector<char *> query_kmers;
            generate_kmers(tag_1->con, kmer_len, num_kmers, query_kmers);

            map<int, int> hits;
            //
            // Lookup the occurances of each k-mer in the kmer_map
            //
            for (int j = 0; j < num_kmers; j++) {
                for (uint k = 0; k < kmer_map[query_kmers[j]].size(); k++)
                    hits[kmer_map[query_kmers[j]][k]]++;
            }

            //
            // Free the k-mers we generated for this query
            //
            for (int j = 0; j < num_kmers; j++)
                delete [] query_kmers[j];

            //
            // Iterate through the list of hits. For each hit that has more than min_hits
            // check its full length to verify a match.
            //
            map<int, int>::iterator hit_it;
            for (hit_it = hits.begin(); hit_it != hits.end(); hit_it++) {

                if (hit_it->second < min_hits) continue;

                tag_2 = merged[hit_it->first];

                // Don't compute distances for masked tags
                if (tag_2->masked) continue;

                // Don't compare tag_1 against itself.
                if (tag_1 == tag_2) continue;

                if (aln->align(tag_1->con, tag_2->con))
		    tag_1->alns.push_back(make_pair(tag_2->id, aln->cigar));
            }
        }

	delete aln;
    }

    free_kmer_hash(kmer_map, kmer_map_keys);

    return 0;
}

// int
// init_alignment(int len, double ***matrix, AlignPath ***path)
// {
//     uint m = len + 1;
//     uint n = len + 1;

//     *matrix = new double * [m];
//     for (uint i = 0; i < m; i++)
//         (*matrix)[i] = new double [n];

//     *path = new AlignPath * [m];
//     for (uint i = 0; i < m; i++)
//         (*path)[i] = new AlignPath [n];

//     return 0;
// }

// int
// free_alignment(int m, double **matrix, AlignPath **path)
// {
//     for (int i = 0; i < m; i++) {
//         delete [] matrix[i];
//         delete [] path[i];
//     }
//     delete [] matrix;
//     delete [] path;

//     return 0;
// }

// int
// align(MergedStack *tag_1, MergedStack *tag_2, double **matrix, AlignPath **path)
// {
//     //         j---->
//     //        [0][1][2][3]...[n-1]
//     //       +--------------------
//     // i [0] | [i][j]
//     // | [1] |
//     // | [2] |
//     // v [3] |
//     //   ... |
//     // [m-1] |
//     // 
//     uint m = tag_1->len + 1;
//     uint n = tag_2->len + 1;
    
//     //
//     // Initialize the first column and row of the dynamic programming
//     // matrix and the path array.
//     //
//     path[0][0].diag = false;
//     path[0][0].up   = false;
//     path[0][0].left = false;
//     matrix[0][0]    = 0.0;
//     for (uint i = 1; i < m; i++) {
//         matrix[i][0]    = path[i - 1][0].up ? matrix[i - 1][0] + gapext_score : matrix[i - 1][0] + gapopen_score;
//         path[i][0].diag = false;
//         path[i][0].up   = true;
//         path[i][0].left = false;
//     }
//     for (uint j = 1; j < n; j++) {
//         matrix[0][j]    = path[0][j - 1].left ? matrix[0][j - 1] + gapext_score : matrix[0][j - 1] + gapopen_score;
//         path[0][j].diag = false;
//         path[0][j].up   = false;
//         path[0][j].left = true;
//     }

//     double  score_down, score_diag, score_right;
//     double  scores[3];
//     dynprog direction[3];
    
//     for (uint i = 1; i < m; i++) {
//         for (uint j = 1; j < n; j++) {
//             // Calculate the score:
//             //   1) If we were to move down from the above cell.
//             score_down   = matrix[i - 1][j];
//             score_down  += path[i - 1][j].up ?  gapext_score : gapopen_score;
//             //   2) If we were to move diagonally from the above and left cell.
//             score_diag   = matrix[i - 1][j - 1] + (tag_1->con[i - 1] == tag_2->con[j - 1] ? match_score : mismatch_score);
//             //   3) If we were to move over from the cell left of us.
//             score_right  = matrix[i][j - 1];
//             score_right += path[i][j - 1].left ? gapext_score : gapopen_score;

//             //
//             // Sort the scores, highest to lowest.
//             //
//             scores[0]    = score_down;
//             direction[0] = dynp_down;
//             scores[1]    = score_diag;
//             direction[1] = dynp_diag;
//             scores[2]    = score_right;
//             direction[2] = dynp_right;

//             if (scores[0] < scores[1])
//                 swap(scores, direction, 0, 1);
//             if (scores[1] < scores[2])
//                 swap(scores, direction, 1, 2);
//             if (scores[0] < scores[1])
//                 swap(scores, direction, 0, 1);

//             matrix[i][j] = scores[0];

//             if (scores[0] > scores[1]) {
//                 //
//                 // One path is best.
//                 //
//                 switch (direction[0]) {
//                 case dynp_diag:
//                     path[i][j].diag = true;
//                     path[i][j].up   = false;
//                     path[i][j].left = false;
//                     break;
//                 case dynp_down:
//                     path[i][j].diag = false;
//                     path[i][j].up   = true;
//                     path[i][j].left = false;
//                     break;
//                 case dynp_right:
//                 default:
//                     path[i][j].diag = false;
//                     path[i][j].up   = false;
//                     path[i][j].left = true;
//                 }
                
//             } else if (scores[0] == scores[1]) {
//                 //
//                 // Two of the paths are equivalent.
//                 //
//                 switch (direction[0]) {
//                 case dynp_diag:
//                     path[i][j].diag = true;
                    
//                     switch (direction[1]) {
//                     case dynp_down:
//                         path[i][j].up   = true;
//                         path[i][j].left = false;
//                         break;
//                     default:
//                     case dynp_right:
//                         path[i][j].up   = false;
//                         path[i][j].left = true;
//                         break;
//                     }
//                     break;
//                 case dynp_down:
//                     path[i][j].up = true;
                    
//                     switch (direction[1]) {
//                     case dynp_right:
//                         path[i][j].diag  = false;
//                         path[i][j].left = true;
//                         break;
//                     default:
//                     case dynp_diag:
//                         path[i][j].diag  = true;
//                         path[i][j].left = false;
//                         break;
//                     }
//                     break;
//                 default:
//                 case dynp_right:
//                     path[i][j].left = true;
                    
//                     switch (direction[1]) {
//                     case dynp_diag:
//                         path[i][j].diag = true;
//                         path[i][j].up   = false;
//                         break;
//                     default:
//                     case dynp_down:
//                         path[i][j].diag = false;
//                         path[i][j].up   = true;
//                         break;
//                     }
//                     break;
//                 }
                
//             } else {
//                 //
//                 // All paths equivalent.
//                 //
//                 path[i][j].diag = true;
//                 path[i][j].up   = true;
//                 path[i][j].left = true;
//             }
//         }
//     }

//     // dump_alignment(tag_1, tag_2, matrix, path);

//     trace_alignment(tag_1, tag_2, path);
    
//     return 0;
// }

// inline int
// swap(double *scores, dynprog *direction, int index_1, int index_2)
// {
//     double swap        = scores[index_1];
//     scores[index_1]    = scores[index_2];
//     scores[index_2]    = swap;
//     dynprog swapdir    = direction[index_1];
//     direction[index_1] = direction[index_2];
//     direction[index_2] = swapdir;

//     return 0;
// }

// int
// trace_alignment(MergedStack *tag_1, MergedStack *tag_2, AlignPath **path)
// {
//     //         j---->
//     //        [0][1][2][3]...[n-1]
//     //       +--------------------
//     // i [0] | [i][j]
//     // | [1] |
//     // | [2] |
//     // v [3] |
//     //   ... |
//     // [m-1] |
//     // 
//     int    m = tag_1->len + 1;
//     int    n = tag_2->len + 1;
//     int    i, j, cnt, len, gaps;
//     string cigar;
//     char   buf[id_len];

//     vector<pair<string, int> > alns;
//     bool more_paths = true;
    
//     do {
//         more_paths = false;

//         i = m - 1;
//         j = n - 1;

//         string aln_1, aln_2;

//         while (i > 0 || j > 0) {
//             cnt  = path[i][j].count();

//             if (cnt > 1) more_paths = true;

//             if (path[i][j].diag) {
//                 aln_1 += tag_1->con[i - 1];
//                 aln_2 += tag_2->con[j - 1];
//                 if (cnt > 1) path[i][j].diag = false;
//                 i--;
//                 j--;
//             } else if (path[i][j].up) {
//                 aln_1 += tag_1->con[i - 1];
//                 aln_2 += "-";
//                 if (cnt > 1) path[i][j].up = false;
//                 i--;
//             } else if (path[i][j].left) {
//                 aln_1 += "-";
//                 aln_2 += tag_2->con[j - 1];
//                 if (cnt > 1) path[i][j].left = false;
//                 j--;
//             }
//         }

//         reverse(aln_1.begin(), aln_1.end());
//         reverse(aln_2.begin(), aln_2.end());

//         //
//         // Convert to CIGAR strings.
//         //
//         cigar = "";
//         len   = aln_1.length();
//         gaps  = 0;
//         i     = 0;
//         while (i < len) {
//             if (aln_1[i] != '-' && aln_2[i] != '-') {
//                 cnt = 0;
//                 do {
//                     cnt++;
//                     i++;
//                 } while (i < len && aln_1[i] != '-' && aln_2[i] != '-');
//                 sprintf(buf, "%dM", cnt);

//             } else if (aln_1[i] == '-') {
//                 cnt = 0;
//                 do {
//                     cnt++;
//                     i++;
//                 } while (i < len && aln_1[i] == '-');
//                 sprintf(buf, "%dD", cnt);
//                 gaps++;

//             } else {
//                 cnt = 0;
//                 do {
//                     cnt++;
//                     i++;
//                 } while (i < len && aln_2[i] == '-');
//                 sprintf(buf, "%dI", cnt);
//                 gaps++;
//             }

//             cigar += buf;
//         }

//         alns.push_back(make_pair(cigar, gaps));
        
//         // cerr << aln_1 << " [" << cigar << ", gaps: " << gaps << "]\n"
//         //      << aln_2 << "\n";

//     } while (more_paths);

//     cigar = "";

//     if (alns.size() == 1) {
//         cigar = alns[0].first;
//         // cerr << "Final alignment: " << cigar << "; gaps: " << alns[0].second << "\n";

//     } else {
//         sort(alns.begin(), alns.end(), compare_pair_stringint);
//         if (alns[0].second < alns[1].second) {
//             cigar = alns[0].first;
//             // cerr << "Final alignment: " << cigar << "; gaps: " << alns[0].second << "\n";
//         }
//     }

//     if (cigar.length() > 0) {
//         tag_1->alns.push_back(make_pair(tag_2->id, cigar));

//         return 1;
//     }

//     return 0;
// }

// int
// dump_alignment(MergedStack *tag_1, MergedStack *tag_2, double **matrix, AlignPath **path)
// {
//     //         j---->
//     //        [0][1][2][3]...[n-1]
//     //       +--------------------
//     // i [0] | [i][j]
//     // | [1] |
//     // | [2] |
//     // v [3] |
//     //   ... |
//     // [m-1] |
//     // 
//     uint m = tag_1->len + 1;
//     uint n = tag_2->len + 1;

//     //
//     // Output the score matrix.
//     //
//     cout << "         ";
//     for (uint j = 0; j < tag_2->len; j++)
//         cout << "   " << tag_2->con[j] << "  |";
//     cout << "\n";

//     cout << "  ";
//     for (uint j = 0; j < n; j++)
//         printf("% 6.1f|", matrix[0][j]);
//     cout << "\n";

//     for (uint i = 1; i < m; i++) {
//         cout << tag_1->con[i - 1] << " ";
//         for (uint j = 0; j < n; j++)
//             printf("% 6.1f|", matrix[i][j]);
//         cout << "\n";
//     }

//     cout << "\n";

//     //
//     // Output the path matrix.
//     //
//     cout << "       ";
//     for (uint j = 0; j < tag_2->len; j++)
//         cout << "  " << tag_2->con[j] << " |";
//     cout << "\n";

//     cout << "  ";
//     for (uint j = 0; j < n; j++) {
//         cout << " ";
//         path[0][j].diag ? cout << "d" : cout << " ";
//         path[0][j].up   ? cout << "u" : cout << " ";
//         path[0][j].left ? cout << "l" : cout << " ";
//         cout << "|";
//     }
//     cout << "\n";

//     for (uint i = 1; i < m; i++) {
//         cout << tag_1->con[i - 1] << " ";
//         for (uint j = 0; j < n; j++) {
//             cout << " ";
//             path[i][j].diag ? cout << "d" : cout << " ";
//             path[i][j].up   ? cout << "u" : cout << " ";
//             path[i][j].left ? cout << "l" : cout << " ";
//             cout << "|";
//         }
//         cout << "\n";
//     }

//     cout << "\n";
    
//     return 0;
// }

int
merge_remainders(map<int, MergedStack *> &merged, map<int, Rem *> &rem)
{
    map<int, Rem *>::iterator it;
    int j, k;

    // OpenMP can't parallelize random access iterators, so we convert
    // our map to a vector of integer keys.
    vector<int> keys;
    uint tot = 0;
    for (it = rem.begin(); it != rem.end(); it++) {
        keys.push_back(it->first);
        tot += it->second->count();
    }

    cerr << "  " << tot << " remainder sequences left to merge.\n";

    if (max_rem_dist <= 0) {
        cerr << "  Matched 0 remainder reads; unable to match " << tot << " remainder reads.\n";
        return 0;
    }

    //
    // Calculate the number of k-mers we will generate. If kmer_len == 0,
    // determine the optimal length for k-mers.
    //
    int con_len   = strlen(merged.begin()->second->con);
    if (set_kmer_len) kmer_len = determine_kmer_length(con_len, max_rem_dist);
    int num_kmers = con_len - kmer_len + 1;

    //
    // Calculate the minimum number of matching k-mers required for a possible sequence match.
    //
    int min_hits = calc_min_kmer_matches(kmer_len, max_rem_dist, con_len, set_kmer_len ? true : false);

    cerr << "  Distance allowed between stacks: " << max_rem_dist << "; searching with a k-mer length of " << kmer_len << " (" << num_kmers << " k-mers per read); " << min_hits << " k-mer hits required.\n";
    
    KmerHashMap    kmer_map;
    vector<char *> kmer_map_keys;
    populate_kmer_hash(merged, kmer_map, kmer_map_keys, kmer_len);
    int utilized = 0;

    //
    // Create a character buffer to hold the Rem sequence, this is faster
    // than repeatedly decoding the DNASeq buffers.
    // 
    //it = rem.find(keys[0]);
    //char *buf = new char[it->second->seq->size + 1];

    #pragma omp parallel private(it, k)
    { 
        #pragma omp for schedule(dynamic) 
        for (j = 0; j < (int) keys.size(); j++) {
            it = rem.find(keys[j]);
            Rem  *r   = it->second;
            char *buf = new char[r->seq->size() + 1];

            //
            // Generate the k-mers for this remainder sequence
            //
            vector<char *> rem_kmers;
            buf = r->seq->seq(buf);
            generate_kmers(buf, kmer_len, num_kmers, rem_kmers);

            map<int, int> hits;
            vector<int>::iterator map_it;
            //
            // Lookup the occurances of each remainder k-mer in the MergedStack k-mer map
            //
            for (k = 0; k < num_kmers; k++) {
                if (kmer_map.find(rem_kmers[k]) != kmer_map.end())
                    for (map_it  = kmer_map[rem_kmers[k]].begin();
                         map_it != kmer_map[rem_kmers[k]].end();
                         map_it++)
                        hits[*map_it]++;
            }

            //
            // Iterate through the list of hits. For each hit that has more than min_hits
            // check its full length to verify a match.
            //
            map<int, int> dists;
            map<int, int>::iterator hit_it;
            for (hit_it = hits.begin(); hit_it != hits.end(); hit_it++) {
                if (hit_it->second < min_hits) continue;

                int d = dist(merged[hit_it->first], buf);
                //
                // Store the distance between these two sequences if it is
                // below the maximum distance
                //
                if (d <= max_rem_dist) {
                    dists[hit_it->first] = d;
                }
            }

            //
            // Free the k-mers we generated for this remainder read
            //
            for (k = 0; k < num_kmers; k++)
                delete [] rem_kmers[k];

            // Check to see if there is a uniquely low distance, if so,
            // merge this remainder tag. If not, discard it, since we
            // can't locate a single best-fitting Stack to merge it into.
            map<int, int>::iterator s;
            int min_id = -1;
            int count  =  0;
            int dist   =  max_rem_dist + 1;

            for (s = dists.begin(); s != dists.end(); s++) {
                if ((*s).second < dist) {
                    min_id = (*s).first;
                    count  = 1;
                    dist   = (*s).second;
                } else if ((*s).second == dist) {
                    count++;
                }
            }

            delete [] buf;

            // Found a merge partner.
            if (min_id >= 0 && count == 1) {
                r->utilized = true;
                #pragma omp critical
                {
                    merged[min_id]->remtags.push_back(it->first);
                    utilized += it->second->count();
                }
            }
        }
    }

    free_kmer_hash(kmer_map, kmer_map_keys);
    //delete [] buf;

    cerr << "  Matched " << utilized << " remainder reads; unable to match " << tot - utilized << " remainder reads.\n";

    return 0;
}

int
call_alleles(MergedStack *mtag, vector<DNANSeq *> &reads, vector<read_type> &read_types)
{
    int     row;
    int     height = reads.size();
    string  allele;
    DNANSeq *d;
    char    base;
    vector<SNP *>::iterator snp;

    for (row = 0; row < height; row++) {
        allele.clear();

        uint snp_cnt = 0;

        //
        // Only call a haplotype from primary reads.
        //
        if (!call_sec_hapl && read_types[row] == secondary) continue;

        for (snp = mtag->snps.begin(); snp != mtag->snps.end(); snp++) {
            if ((*snp)->type != snp_type_het) continue;

            snp_cnt++;

            d    = reads[row];
            base = (*d)[(*snp)->col];

            //
            // Check to make sure the nucleotide at the location of this SNP is
            // of one of the two possible states the multinomial model called.
            //
            if (base == (*snp)->rank_1 || base == (*snp)->rank_2) 
                allele += base;
            else
                break;
        }

        if (snp_cnt > 0 && allele.length() == snp_cnt)
            mtag->alleles[allele]++;
    }

    return 0;
}

int
call_consensus(map<int, MergedStack *> &merged, map<int, Stack *> &unique, map<int, Rem *> &rem, bool invoke_model)
{
    //
    // OpenMP can't parallelize random access iterators, so we convert
    // our map to a vector of integer keys.
    //
    map<int, MergedStack *>::iterator it;
    vector<int> keys;
    for (it = merged.begin(); it != merged.end(); it++) 
        keys.push_back(it->first);

    int i;
    #pragma omp parallel private(i)
    { 
        MergedStack *mtag;
        Stack       *utag;
        Rem         *r;

        #pragma omp for schedule(dynamic) 
        for (i = 0; i < (int) keys.size(); i++) {
            mtag = merged[keys[i]];
            
            //
            // Create a two-dimensional array, each row containing one read. For
            // each unique tag that has been merged together, add the sequence for
            // that tag into our array as many times as it originally occurred. 
            //
            vector<int>::iterator j;
            vector<DNANSeq *>  reads;
            vector<read_type> read_types;

            for (j = mtag->utags.begin(); j != mtag->utags.end(); j++) {
                utag = unique[*j];

                for (uint k = 0; k < utag->count(); k++) {
                    reads.push_back(utag->seq);
                    read_types.push_back(primary);
                }
            }

            // For each remainder tag that has been merged into this Stack, add the sequence. 
            for (j = mtag->remtags.begin(); j != mtag->remtags.end(); j++) {
                r = rem[*j];

                for (uint k = 0; k < r->count(); k++) {
                    reads.push_back(r->seq);
                    read_types.push_back(secondary);
                }
            }

            //
            // Iterate over each column of the array and call the consensus base.
            //
            uint row, col;
            uint length = reads[0]->size();
            uint height = reads.size();
            string con;
            map<char, int> nuc;
            map<char, int>::iterator max, n;
            DNANSeq *d;

            uint cur_gap = mtag->gaps.size() > 0 ? 0 : 1;

            for (col = 0; col < length; col++) {
                //
                // Don't invoke the model within gaps.
                //
                if (cur_gap < mtag->gaps.size() && col == mtag->gaps[cur_gap].start) {
                    do {
                        con += 'N';
                        SNP *snp    = new SNP;
                        snp->type   = snp_type_unk;
                        snp->col    = col;
			snp->rank_1 = '-';
			snp->rank_2 = '-';
                        mtag->snps.push_back(snp);
                        col++;
                    } while (col < mtag->gaps[cur_gap].end && col < length);
                    col--;
                    cur_gap++;
                    continue;
                }

                nuc['A'] = 0; 
                nuc['G'] = 0;
                nuc['C'] = 0;
                nuc['T'] = 0;

                for (row = 0; row < height; row++) {
                    d = reads[row];
                    if (nuc.count((*d)[col]))
                        nuc[(*d)[col]]++;
                }

                //
                // Find the base with a plurality of occurances and call it.
                //
                max = nuc.end();

                for (n = nuc.begin(); n != nuc.end(); n++) {

                    if (max == nuc.end() || n->second > max->second)
                        max = n;
                }
                con += max->first;

                //
                // Search this column for the presence of a SNP
                //
                if (invoke_model) 
                    switch(model_type) {
                    case snp:
                        call_multinomial_snp(mtag, col, nuc, true);
                        break;
                    case bounded:
                        call_bounded_multinomial_snp(mtag, col, nuc, true);
                        break;
                    case fixed:
                        call_multinomial_fixed(mtag, col, nuc);
                        break;
                    }
            }

            if (invoke_model) {
                call_alleles(mtag, reads, read_types);

                if (model_type == fixed) {
                    //
                    // Mask nucleotides that are not fixed.
                    //
                    vector<SNP *>::iterator s;
                    for (s = mtag->snps.begin(); s != mtag->snps.end(); s++) {
                        if ((*s)->type == snp_type_unk)
                            con.replace((*s)->col, 1, "N");
                    }
                }
            }

            mtag->add_consensus(con.c_str());
        }
    }

    return 0;
}

int
populate_merged_tags(map<int, Stack *> &unique, map<int, MergedStack *> &merged)
{
    map<int, Stack *>::iterator i;
    map<int, MergedStack *>::iterator it_new, it_old;
    Stack       *utag;
    MergedStack *mtag;
    int          k = 0;

    it_old = merged.begin();

    for (i = unique.begin(); i != unique.end(); i++) {
        utag = (*i).second;
        mtag = new MergedStack;

        mtag->id    = k;
        mtag->count = utag->count();
        mtag->utags.push_back(utag->id);
        mtag->add_consensus(utag->seq);

        // Insert the new MergedStack giving a hint as to which position
        // to insert it at.
        it_new = merged.insert(it_old, pair<int, MergedStack *>(k, mtag));
        it_old = it_new;
        k++;
    }

    return 0;
}

int
merge_stacks(map<int, Stack *> &unique, map<int, Rem *> &rem, map<int, MergedStack *> &merged, set<int> &merge_map, int round)
{
    map<int, MergedStack *> new_merged;
    map<int, MergedStack *>::iterator it, it_old, it_new;
    MergedStack *tag_1, *tag_2;
    vector<set<int> > merge_lists;

    uint index     = 0;
    int  cohort_id  = 0;
    int  id         = 1;
    uint delev_cnt = 0;
    uint blist_cnt = 0;

    for (it = merged.begin(); it != merged.end(); it++) {
        tag_1 = it->second;

        //
        // This tag may already have been merged by an earlier operation.
        //
        if (merge_map.count(tag_1->id) > 0)
            continue;

        queue<int>                        merge_list;
        pair<set<int>::iterator,bool>     ret;
        vector<pair<int, int> >::iterator k;

        merge_lists.push_back(set<int>());

        if (tag_1->masked) {
            merge_lists[index].insert(tag_1->id);
            index++;
            continue;
        }

        //
        // Construct a list of MergedStacks that are within a particular distance
        // of this tag.
        //
        merge_lists[index].insert(tag_1->id);
        merge_list.push(tag_1->id);

        while (!merge_list.empty()) {
            tag_2 = merged[merge_list.front()];
            merge_list.pop();

            for (k = tag_2->dist.begin(); k != tag_2->dist.end(); k++) {
                ret = merge_lists[index].insert(k->first);

                //
                // If this Tag has not already been added to the merge list (i.e. we were able
                // to insert it in to our unique_merge_list, which is a set), add it for consideration
                // later in the loop.
                //
                if (ret.second == true)
                    merge_list.push((*k).first);
            }
        }

        //
        // Record the nodes that have been merged in this round.
        //
        set<int>::iterator j;
        for (j = merge_lists[index].begin(); j != merge_lists[index].end(); j++)
            merge_map.insert(*j);

        index++;
    }

    #pragma omp parallel private(tag_1, tag_2)
    {
        vector<MergedStack *> merged_tags;

        #pragma omp for reduction(+:delev_cnt) reduction(+:blist_cnt)
        for (uint index = 0; index < merge_lists.size(); index++) {
            //
            // Deal with the simple case of a single locus that does not need to be merged.
            //
            if (merge_lists[index].size() == 1) {
                tag_1 = merged[*(merge_lists[index].begin())];
                tag_2 = merge_tags(merged, merge_lists[index], 0);

                //
                // If this tag is masked, keep the old cohort_id.
                //
                if (tag_1->masked) {
                    tag_2->cohort_id = tag_1->cohort_id;
                } else {
                    tag_2->cohort_id = cohort_id;
                    #pragma omp atomic
                    cohort_id++;
                }
                merged_tags.push_back(tag_2);
                continue;
            }

            //
            // Break large loci down by constructing a minimum
            // spanning tree and severing long distance edges.
            //
            if (deleverage_stacks) {
                vector<MergedStack *> tags;
                bool delev;

                deleverage(unique, rem, merged, merge_lists[index], cohort_id, tags);

                if (tags.size() == 1) {
                    delev = false;
                } else {
                    delev_cnt++;
                    delev = true;
                }

                for (uint t = 0; t < tags.size(); t++) {
                    //tags[t]->id = id;
                    tags[t]->deleveraged = delev;

                    if (tags[t]->utags.size() > max_subgraph) {
                        tags[t]->masked      = true;
                        tags[t]->blacklisted = true;
                        blist_cnt++;
                    }

                    //new_merged.insert(pair<int, MergedStack *>(id, tags[t]));
                    merged_tags.push_back(tags[t]);
                    //id++;
                }

                #pragma omp atomic
                cohort_id++;

            } else {
                //
                // If not deleveraging, merge these tags together into a new MergedStack object.
                //
                tag_2 = merge_tags(merged, merge_lists[index], 0);
                tag_2->cohort_id = cohort_id;

                if (tag_2->utags.size() > max_subgraph) {
                    tag_2->masked      = true;
                    tag_2->blacklisted = true;
                    blist_cnt++;
                }

                //new_merged.insert(pair<int, MergedStack *>(id, tag_2));
                merged_tags.push_back(tag_2);

                #pragma omp atomic
                cohort_id++;
                //id++;
            }
        }

        //
        // Merge the accumulated tags into the new_merged map.
        //
        #pragma omp critical
        {
            it_old = merged.begin();
            for (uint j = 0; j < merged_tags.size(); j++) {
                merged_tags[j]->id = id;
                it_new = new_merged.insert(it_old, pair<int, MergedStack *>(id, merged_tags[j]));
                it_old = it_new;
                id++;
            }
        }
    }

    uint new_cnt = new_merged.size();
    uint old_cnt = merged.size();

    //
    // Free the memory from the old map of merged tags.
    //
    for (it = merged.begin(); it != merged.end(); it++)
        delete it->second;

    merged = new_merged;

    cerr << "  " << old_cnt << " stacks merged into " << new_cnt 
         << " stacks; deleveraged " << delev_cnt 
         << " stacks; removed " << blist_cnt << " stacks.\n";

    return 0;
}

MergedStack *
merge_tags(MergedStack *tag_1, MergedStack *tag_2, int id)
{
    MergedStack *new_tag;

    new_tag     = new MergedStack;
    new_tag->id = id;

    new_tag->deleveraged      = tag_2->deleveraged      || tag_1->deleveraged;
    new_tag->masked           = tag_2->masked           || tag_1->masked;
    new_tag->blacklisted      = tag_2->blacklisted      || tag_1->blacklisted;
    new_tag->gappedlumberjack = tag_2->gappedlumberjack || tag_1->gappedlumberjack;
    new_tag->lumberjackstack  = tag_2->lumberjackstack  || tag_1->lumberjackstack;

    for (uint i = 0; i < tag_1->utags.size(); i++)
        new_tag->utags.push_back(tag_1->utags[i]);

    for (uint i = 0; i < tag_1->remtags.size(); i++)
        new_tag->remtags.push_back(tag_1->remtags[i]);

    new_tag->count = tag_1->count;

    for (uint i = 0; i < tag_2->utags.size(); i++)
        new_tag->utags.push_back(tag_2->utags[i]);

    for (uint i = 0; i < tag_2->remtags.size(); i++)
        new_tag->remtags.push_back(tag_2->remtags[i]);

    new_tag->count += tag_2->count;

    return new_tag;
}

MergedStack *merge_tags(map<int, MergedStack *> &merged, set<int> &merge_list, int id) {
    set<int>::iterator    i;
    vector<int>::iterator j;
    MergedStack *tag_1, *tag_2;

    tag_1     = new MergedStack;
    tag_1->id = id;

    for (i = merge_list.begin(); i != merge_list.end(); i++) {
        tag_2 = merged[(*i)];

        tag_1->deleveraged     = tag_2->deleveraged     ? true : tag_1->deleveraged;
        tag_1->masked          = tag_2->masked          ? true : tag_1->masked;
        tag_1->blacklisted     = tag_2->blacklisted     ? true : tag_1->blacklisted;
        tag_1->lumberjackstack = tag_2->lumberjackstack ? true : tag_1->lumberjackstack;

        for (j = tag_2->utags.begin(); j != tag_2->utags.end(); j++)
            tag_1->utags.push_back(*j);

        for (j = tag_2->remtags.begin(); j != tag_2->remtags.end(); j++)
            tag_1->remtags.push_back(*j);

        tag_1->count += tag_2->count;
    }

    return tag_1;
}

MergedStack *merge_tags(map<int, MergedStack *> &merged, int *merge_list, int merge_list_size, int id) {
    int                   i;
    vector<int>::iterator j;
    MergedStack *tag_1, *tag_2;

    tag_1     = new MergedStack;
    tag_1->id = id;

    for (i = 0; i < merge_list_size; i++) {
        tag_2 = merged[merge_list[i]];

        tag_1->deleveraged     = tag_2->deleveraged     ? true : tag_1->deleveraged;
        tag_1->masked          = tag_2->masked          ? true : tag_1->masked;
        tag_1->blacklisted     = tag_2->blacklisted     ? true : tag_1->blacklisted;
        tag_1->lumberjackstack = tag_2->lumberjackstack ? true : tag_1->lumberjackstack;

        for (j = tag_2->utags.begin(); j != tag_2->utags.end(); j++)
            tag_1->utags.push_back(*j);

        for (j = tag_2->remtags.begin(); j != tag_2->remtags.end(); j++)
            tag_1->remtags.push_back(*j);

        tag_1->count += tag_2->count;
    }

    return tag_1;
}

int
remove_repetitive_stacks(map<int, Stack *> &unique, map<int, MergedStack *> &merged)
{
    //
    // If enabled, check the depth of coverage of each unique tag, and remove 
    // from consideration any tags with depths greater than removal_trigger. These tags
    // are likely to be multiple repetitive sites that have been merged together. 
    // Because large stacks of unique tags are likely to also generate many one-off
    // sequencing error reads, remove all seqeunces that are a distance of one from
    // the RAD-Tag with high depth of coverage.
    //
    map<int, MergedStack *>::iterator i;
    vector<pair<int, int> >::iterator k;
    map<int, MergedStack *> new_merged;
    MergedStack *tag_1, *tag_2;
    set<int> already_merged;

    //
    // First, iterate through the stacks and populate a list of tags that will be removed
    // (those above the removal_trigger and those 1 nucleotide away). Sort the list of
    // stacks so that we process them from largest depth to shortest so the same stacks
    // are always merged/removed..
    //
    vector<pair<int, int> > ordered_tags;
    for (i = merged.begin(); i != merged.end(); i++) {
        if (i->second->count > removal_trigger)
            ordered_tags.push_back(make_pair(i->second->id, i->second->count));
    }
    sort(ordered_tags.begin(), ordered_tags.end(), compare_pair_intint);

    pair<set<int>::iterator,bool> ret;
    int id = 0;

    //
    // Merge all stacks that are over the removal trigger with their nearest neighbors and
    // mask them so they are not further processed by the program.
    //
    for (uint j = 0; j < ordered_tags.size(); j++) {
        tag_1 = merged[ordered_tags[j].first];

        //
        // Don't process a tag that has already been merged.
        //
        if (already_merged.count(tag_1->id) > 0)
            continue;

        //
        // Construct a list of MergedStacks that are either:
        //   within a distance of 1 nucleotide of this tag, or
        //   are they themselves above the lumberjack stacks limit.
        //
        queue<int> merge_queue;
        set<int>   merge_list;
        merge_queue.push(tag_1->id);
        merge_list.insert(tag_1->id);
        already_merged.insert(tag_1->id);

        while (!merge_queue.empty()) {
            tag_2 = merged[merge_queue.front()];
            merge_queue.pop();

            if (tag_2->count < removal_trigger)
                continue;

            for (k = tag_2->dist.begin(); k != tag_2->dist.end(); k++) {
                ret = already_merged.insert(k->first);

                if (ret.second == true) {
                    merge_queue.push(k->first);
                    merge_list.insert(k->first);
                }
            }
        }
        
        //
        // Merge these tags together into a new MergedStack object.
        //
        tag_2 = merge_tags(merged, merge_list, id);
        tag_2->add_consensus(tag_1->con);

        tag_2->lumberjackstack = true;
        tag_2->masked          = true;
        tag_2->blacklisted     = true;

        new_merged.insert(make_pair(id, tag_2));
        id++;
    }

    //
    // Move the non-lumberjack stacks, unmodified, into the new merged map.
    //
    for (i = merged.begin(); i != merged.end(); i++) {
        tag_1 = i->second;

        if (already_merged.count(tag_1->id) > 0)
            continue;

        set<int> merge_list;
        merge_list.insert(tag_1->id);

        tag_2 = merge_tags(merged, merge_list, id);
        tag_2->add_consensus(tag_1->con);

        new_merged.insert(make_pair(id, tag_2));
        id++;
    }

    cerr << "  Removed " << already_merged.size() << " stacks.\n";

    //
    // Free the memory from the old map of merged tags.
    //
    map<int, MergedStack *>::iterator it;
    for (it = merged.begin(); it != merged.end(); it++)
        delete it->second;

    merged = new_merged;

    cerr << "  " << merged.size() << " stacks remain for merging.\n";

    return 0;
}

int deleverage(map<int, Stack *> &unique, 
               map<int, Rem *> &rem,
               map<int, MergedStack *> &merged, 
               set<int> &merge_list, 
               int cohort_id, 
               vector<MergedStack *> &deleveraged_tags) {
    set<int>::iterator it;
    vector<pair<int, int> >::iterator j;
    MergedStack *tag_1, *tag_2;
    uint k, l;

    //
    // Create a minimum spanning tree in order to determine the minimum distance
    // between each node in the list.
    //
    MinSpanTree *mst = new MinSpanTree;
    vector<int>  keys;

    for (it = merge_list.begin(); it != merge_list.end(); it++) {
        keys.push_back(*it);
        mst->add_node(*it);
        tag_1 = merged[*it];

        // cerr << "  " << *it << " -> " << tag_1->utags[0] << "\n";
    }

    //
    // Measure the distance between each pair of nodes and add edges to our
    // minimum spanning tree.
    //
    Node *n_1, *n_2;
    for (k = 0; k < keys.size(); k++) {
        tag_1 = merged[keys[k]];
        n_1   = mst->node(keys[k]);

        for (l = k+1; l < keys.size(); l++) {
            tag_2 = merged[keys[l]];
            n_2   = mst->node(keys[l]);

            int d = dist(tag_1, tag_2);

            n_1->add_edge(mst->node(keys[l]), d);
            n_2->add_edge(mst->node(keys[k]), d);
        }
    }

    //
    // Build the minimum spanning tree.
    //
    mst->build_tree();

    //
    // Visualize the MST
    //
    if (dump_graph) {
        stringstream gout_file;
        size_t pos_1 = in_file.find_last_of("/");
        size_t pos_2 = in_file.find_last_of(".");
        gout_file << out_path << in_file.substr(pos_1 + 1, (pos_2 - pos_1 - 1)) << "_" << keys[0] << ".dot";
        string vis = mst->vis(true);
        ofstream gvis(gout_file.str().c_str());
        gvis << vis;
        gvis.close();
    }

    set<int> visited;
    set<int, int_increasing> dists;
    queue<Node *> q;

    Node *n = mst->head();
    q.push(n);

    while (!q.empty()) {
        n = q.front();
        q.pop();
        visited.insert(n->id);

        for (uint i = 0; i < n->min_adj_list.size(); i++) {
            if (visited.count(n->min_adj_list[i]->id) == 0) {
                q.push(n->min_adj_list[i]);
                // cerr << n->id << " -> " << n->min_adj_list[i]->id << ": ";

                //
                // Find the edge distance.
                //
                for (uint j = 0; j < n->edges.size(); j++)
                    if (n->edges[j]->child == n->min_adj_list[i]) {
                        // cerr << n->edges[j]->dist << "\n";
                        dists.insert(n->edges[j]->dist);
                    }
            }
        }
    }

    //
    // This set is sorted by definition. Check if there is more than a single 
    // distance separating stacks.
    //
    if (dists.size() == 1) {
        tag_1 = merge_tags(merged, merge_list, 0);
        deleveraged_tags.push_back(tag_1);
        delete mst;
        return 0;
    }

    uint min_dist = *(dists.begin());

    //
    // If there is more than a single distance, split the minimum spanning tree
    // into separate loci, by cutting the tree at the larger distances.
    //
    set<int> uniq_merge_list;
    visited.clear();
    n = mst->head();
    q.push(n);
    int id = 0;

    uniq_merge_list.insert(n->id);
    while (!q.empty()) {
        n = q.front();
        q.pop();
        visited.insert(n->id);

        for (uint i = 0; i < n->min_adj_list.size(); i++) {
            if (visited.count(n->min_adj_list[i]->id) == 0) {
                q.push(n->min_adj_list[i]);

                for (uint j = 0; j < n->edges.size(); j++) {
                    if (n->edges[j]->child == n->min_adj_list[i])
                        if (n->edges[j]->dist > min_dist) {

                            // cerr << "Merging the following stacks into a locus:\n";
                            for (it = uniq_merge_list.begin(); it != uniq_merge_list.end(); it++) {
                                tag_1 = merged[*it];
                                // cerr << "  " << *it << " -> " << tag_1->utags[0] << "\n";
                            }

                            tag_1 = merge_tags(merged, uniq_merge_list, id);
                            tag_1->cohort_id = cohort_id;
                            deleveraged_tags.push_back(tag_1);
                            uniq_merge_list.clear();
                            id++;
                        }
                }

                uniq_merge_list.insert(n->min_adj_list[i]->id);
            }
        }
    }

    // cerr << "Merging the following stacks into a locus:\n";
    for (it = uniq_merge_list.begin(); it != uniq_merge_list.end(); it++) {
        tag_1 = merged[*it];
        // cerr << "  " << *it << " -> " << tag_1->utags[0] << "\n";
    }

    tag_1 = merge_tags(merged, uniq_merge_list, id);
    tag_1->cohort_id = cohort_id;
    deleveraged_tags.push_back(tag_1);
    uniq_merge_list.clear();

    delete mst;

    return 0;
}

int calc_kmer_distance(map<int, MergedStack *> &merged, int utag_dist) {
    //
    // Calculate the distance (number of mismatches) between each pair
    // of Radtags. We expect all radtags to be the same length;
    //
    KmerHashMap    kmer_map;
    vector<char *> kmer_map_keys;
    MergedStack   *tag_1, *tag_2;
    map<int, MergedStack *>::iterator it;

    // OpenMP can't parallelize random access iterators, so we convert
    // our map to a vector of integer keys.
    vector<int> keys;
    for (it = merged.begin(); it != merged.end(); it++) 
        keys.push_back(it->first);

    //
    // Calculate the number of k-mers we will generate. If kmer_len == 0,
    // determine the optimal length for k-mers.
    //
    int con_len   = strlen(merged[keys[0]]->con);
    if (set_kmer_len)
        kmer_len = determine_kmer_length(con_len, utag_dist);
    int num_kmers = con_len - kmer_len + 1;

    //
    // Calculate the minimum number of matching k-mers required for a possible sequence match.
    //
    int min_hits = calc_min_kmer_matches(kmer_len, utag_dist, con_len, set_kmer_len ? true : false);

    cerr << "  Distance allowed between stacks: " << utag_dist << "; searching with a k-mer length of " << kmer_len << " (" << num_kmers << " k-mers per read); " << min_hits << " k-mer hits required.\n";

    populate_kmer_hash(merged, kmer_map, kmer_map_keys, kmer_len);
 
    #pragma omp parallel private(tag_1, tag_2)
    { 
        #pragma omp for schedule(dynamic) 
        for (uint i = 0; i < keys.size(); i++) {
            tag_1 = merged[keys[i]];

            // Don't compute distances for masked tags
            if (tag_1->masked) continue;

            vector<char *> query_kmers;
            generate_kmers(tag_1->con, kmer_len, num_kmers, query_kmers);

            map<int, int> hits;
            int d;
            //
            // Lookup the occurances of each k-mer in the kmer_map
            //
            for (int j = 0; j < num_kmers; j++) {
                for (uint k = 0; k < kmer_map[query_kmers[j]].size(); k++)
                    hits[kmer_map[query_kmers[j]][k]]++;
            }

            //
            // Free the k-mers we generated for this query
            //
            for (int j = 0; j < num_kmers; j++)
                delete [] query_kmers[j];

            // cerr << "  Tag " << tag_1->id << " hit " << hits.size() << " kmers.\n";

            //
            // Iterate through the list of hits. For each hit that has more than min_hits
            // check its full length to verify a match.
            //
            map<int, int>::iterator hit_it;
            for (hit_it = hits.begin(); hit_it != hits.end(); hit_it++) {
                // cerr << "  Tag " << hit_it->first << " has " << hit_it->second << " hits (min hits: " << min_hits << ")\n";

                if (hit_it->second < min_hits) continue;

                // cerr << "  Match found, checking full-length match\n";

                tag_2 = merged[hit_it->first];

                // Don't compute distances for masked tags
                if (tag_2->masked) continue;

                // Don't compare tag_1 against itself.
                if (tag_1 == tag_2) continue;

                d = dist(tag_1, tag_2);
                // cerr << "    Distance: " << d << "\n";

                //
                // Store the distance between these two sequences if it is
                // below the maximum distance (which governs those
                // sequences to be merged in the following step of the
                // algorithm.)
                //
                if (d <= utag_dist)
                    tag_1->add_dist(tag_2->id, d);
            }
            // Sort the vector of distances.
            sort(tag_1->dist.begin(), tag_1->dist.end(), compare_dist);
        }
    }

    free_kmer_hash(kmer_map, kmer_map_keys);

    return 0;
}

int calc_distance(map<int, MergedStack *> &merged, int utag_dist) {
    //
    // Calculate the distance (number of mismatches) between each pair
    // of Radtags. We expect all radtags to be the same length;
    //
    map<int, MergedStack *>::iterator it;
    MergedStack *tag_1, *tag_2;
    int i, j;

    // OpenMP can't parallelize random access iterators, so we convert
    // our map to a vector of integer keys.
    vector<int> keys;
    for (it = merged.begin(); it != merged.end(); it++) 
        keys.push_back(it->first);

    #pragma omp parallel private(i, j, tag_1, tag_2)
    { 
        #pragma omp for schedule(dynamic) 
        for (i = 0; i < (int) keys.size(); i++) {

            tag_1 = merged[keys[i]];

            // Don't compute distances for masked tags
            if (tag_1->masked) continue;

            int d;

            for (j = 0; j < (int) keys.size(); j++) {
                tag_2 = merged[keys[j]];

                // Don't compute distances for masked tags
                if (tag_2->masked) continue;

                // Don't compare tag_1 against itself.
                if (tag_1 == tag_2) continue;

                d = dist(tag_1, tag_2);
                //cerr << "  Distance: " << d << "\n";

                //
                // Store the distance between these two sequences if it is
                // below the maximum distance (which governs those
                // sequences to be merged in the following step of the
                // algorithm.)
                //
                if (d == utag_dist) {
                    tag_1->add_dist(tag_2->id, d);
                    //cerr << "  HIT.\n";
                }
            }

            // Sort the vector of distances.
            sort(tag_1->dist.begin(), tag_1->dist.end(), compare_dist);
        }
    }

    return 0;
}

int reduce_radtags(DNASeqHashMap &radtags, map<int, Stack *> &unique, map<int, Rem *> &rem) {
    DNASeqHashMap::iterator it;

    Rem   *r;
    Stack *u;
    int   global_id = 1;
    
    for (it = radtags.begin(); it != radtags.end(); it++) {
        if (it->second.count() < min_merge_cov) {
            //
            // Don't record this unique RAD-Tag if its coverage is below
            // the specified cutoff. However, add the reads to the remainder
            // vector for later processing.
            //
            r     = new Rem;
            r->id = global_id;
            r->add_seq(it->first);

            for (uint i = 0; i < it->second.ids.size(); i++)
                r->add_id(it->second.ids[i]);

            rem[r->id] = r;
            global_id++;

        } else {
            //
            // Populate a Stack object for this unique radtag. Create a
            // map of the IDs for the sequences that have been
            // collapsed into this radtag.
            //
            u     = new Stack;
            u->id = global_id;
            u->add_seq(it->first);

            // Copy the original Fastq IDs from which this unique radtag was built.
            for (uint i = 0; i < it->second.ids.size(); i++)
                u->add_id(it->second.ids[i]);

            unique[u->id] = u;
            global_id++;
        }
    }

    if (unique.size() == 0) {
        cerr << "Error: Unable to form any stacks, data appear to be unique.\n";
        exit(1);
    }

    return 0;
}

int
free_radtags_hash(DNASeqHashMap &radtags, vector<DNANSeq *> &radtags_keys)
{
    for (uint i = 0; i < radtags_keys.size(); i++)
        delete radtags_keys[i];

    radtags.clear();

    return 0;
}

int
calc_coverage_distribution(map<int, Stack *> &unique,
                           double &mean, double &stdev, double &max)
{
    map<int, Stack *>::iterator i;
    double m     = 0.0;
    double s     = 0.0;
    double sum   = 0.0;
    uint   cnt   = 0;
    double total = 0.0;

    mean  = 0.0;
    max   = 0.0;
    stdev = 0.0;

    map<int, int> depth_dist;
    map<int, int>::iterator j;

    for (i = unique.begin(); i != unique.end(); i++) {
        cnt = i->second->count();
        m += cnt;
        total++;

        depth_dist[cnt]++;

        if (cnt > max)
            max = cnt;
    }

    mean = m / total;

    //
    // Calculate the standard deviation
    //
    total = 0.0;

    for (i = unique.begin(); i != unique.end(); i++) {
        total++;
        s = i->second->count();
        sum += pow((s - mean), 2);
    }

    stdev = sqrt(sum / (total - 1));

    return 0;
}

int
calc_coverage_distribution(map<int, Stack *> &unique,
                           map<int, MergedStack *> &merged,
                           double &mean, double &stdev, double &max)
{
    map<int, MergedStack *>::iterator it;
    vector<int>::iterator             k;
    Stack *tag;
    double m   = 0.0;
    double s   = 0.0;
    double sum = 0.0;
    double cnt = 0.0;

    mean  = 0.0;
    max   = 0.0;
    stdev = 0.0;

    for (it = merged.begin(); it != merged.end(); it++) {
        if (it->second->blacklisted) continue;

        cnt++;
        m = 0.0;
        for (k = it->second->utags.begin(); k != it->second->utags.end(); k++) {
            tag  = unique[*k];
            m   += tag->count();
        }
        if (m > max) max = m;

        sum += m;
    }

    mean = sum / cnt;

    //
    // Calculate the standard deviation
    //
    for (it = merged.begin(); it != merged.end(); it++) {
        s = 0.0;
        for (k = it->second->utags.begin(); k != it->second->utags.end(); k++) {
            tag  = unique[*k];
            s   += tag->count();
        }
        sum += pow((s - mean), 2);
    }

    stdev = sqrt(sum / (cnt - 1));

    return 0;
}

int
calc_coverage_distribution(map<int, Stack *> &unique,
                           map<int, Rem *> &rem, 
                           map<int, MergedStack *> &merged,
                           double &mean, double &stdev, double &max)
{
    map<int, MergedStack *>::iterator it;
    vector<int>::iterator             k;
    Stack *tag;
    double m    = 0.0;
    double s    = 0.0;
    double sum  = 0.0;
    double cnt  = 0.0;
    
    mean  = 0.0;
    max   = 0.0;
    stdev = 0.0;

    for (it = merged.begin(); it != merged.end(); it++) {
        if (it->second->blacklisted) continue;

        cnt++;
        m = 0.0;
        for (k = it->second->utags.begin(); k != it->second->utags.end(); k++) {
            tag  = unique[*k];
            m   += tag->count();
        }
        for (uint j = 0; j < it->second->remtags.size(); j++)
            m += rem[it->second->remtags[j]]->count();

        if (m > max) max = m;

        sum += m;
    }

    mean = sum / cnt;

    //
    // Calculate the standard deviation
    //
    for (it = merged.begin(); it != merged.end(); it++) {
        s = 0.0;
        for (k = it->second->utags.begin(); k != it->second->utags.end(); k++) {
            tag  = unique[*k];
            s   += tag->count();
        }
        for (uint j = 0; j < it->second->remtags.size(); j++)
            s += rem[it->second->remtags[j]]->count();
        sum += pow((s - mean), 2);
    }

    stdev = sqrt(sum / (cnt - 1));

    return 0;
}

int count_raw_reads(map<int, Stack *> &unique, map<int, Rem *> &rem, map<int, MergedStack *> &merged) {
    map<int, MergedStack *>::iterator it;
    map<int, Stack *>::iterator sit;
    vector<int>::iterator k;
    Stack *tag;
    long int m = 0;

    map<int, int> uniq_ids;
    map<int, int>::iterator uit;

    for (it = merged.begin(); it != merged.end(); it++) {
        for (k = it->second->utags.begin(); k != it->second->utags.end(); k++) {
            tag  = unique[*k];
            m   += tag->count();

            if (uniq_ids.count(*k) == 0)
               uniq_ids[*k] = 0;
            uniq_ids[*k]++;
        }
        for (uint j = 0; j < it->second->remtags.size(); j++)
            m += rem[it->second->remtags[j]]->count();
            //m += it->second->remtags.size();
    }

    for (uit = uniq_ids.begin(); uit != uniq_ids.end(); uit++)
       if (uit->second > 1)
           cerr << "  Unique stack #" << uit->first << " appears in " << uit->second << " merged stacks.\n";

    cerr << "Number of utilized reads: " << m << "\n";

    //for (sit = unique.begin(); sit != unique.end(); sit++)
    //   if (uniq_ids.count(sit->first) == 0)
    //       cerr << "  Stack " << sit->first << ": '" << sit->second->seq << "' unused.\n";

    return 0;
}

int 
write_results(map<int, MergedStack *> &m, map<int, Stack *> &u, map<int, Rem *> &r) 
{
    map<int, MergedStack *>::iterator i;
    vector<int>::iterator      k;
    vector<SNP *>::iterator    s;
    map<string, int>::iterator t;
    MergedStack  *tag_1;
    Stack        *tag_2;
    Rem          *rem;
    stringstream  sstr;

    bool gzip = (in_file_type == FileT::gzfastq || in_file_type == FileT::gzfasta) ? true : false;

    //
    // Read in the set of sequencing IDs so they can be included in the output.
    //
    vector<char *> seq_ids;
    load_seq_ids(seq_ids);

    //
    // Parse the input file name to create the output files
    //
    size_t pos_1 = in_file.find_last_of("/");
    size_t pos_2 = in_file.find_last_of(".");

    if (in_file.substr(pos_2) == ".gz") {
        in_file = in_file.substr(0, pos_2);
        pos_2   = in_file.find_last_of(".");
    }

    string tag_file = out_path + in_file.substr(pos_1 + 1, (pos_2 - pos_1 - 1)) + ".tags.tsv";
    string snp_file = out_path + in_file.substr(pos_1 + 1, (pos_2 - pos_1 - 1)) + ".snps.tsv";
    string all_file = out_path + in_file.substr(pos_1 + 1, (pos_2 - pos_1 - 1)) + ".alleles.tsv";

    if (gzip) {
        tag_file += ".gz";
        snp_file += ".gz";
        all_file += ".gz";
    }

    //
    // Open the output files for writing.
    //
    gzFile   gz_tags, gz_snps, gz_alle;
    ofstream tags, snps, alle;
    if (gzip) {
        gz_tags = gzopen(tag_file.c_str(), "wb");
        if (!gz_tags) {
            cerr << "Error: Unable to open gzipped tag file '" << tag_file << "': " << strerror(errno) << ".\n";
            exit(1);
        }
        #if ZLIB_VERNUM >= 0x1240
        gzbuffer(gz_tags, libz_buffer_size);
        #endif
        gz_snps = gzopen(snp_file.c_str(), "wb");
        if (!gz_snps) {
            cerr << "Error: Unable to open gzipped snps file '" << snp_file << "': " << strerror(errno) << ".\n";
            exit(1);
        }
        #if ZLIB_VERNUM >= 0x1240
        gzbuffer(gz_snps, libz_buffer_size);
        #endif
        gz_alle = gzopen(all_file.c_str(), "wb");
        if (!gz_alle) {
            cerr << "Error: Unable to open gzipped alleles file '" << all_file << "': " << strerror(errno) << ".\n";
            exit(1);
        }
        #if ZLIB_VERNUM >= 0x1240
        gzbuffer(gz_alle, libz_buffer_size);
        #endif
    } else {
        tags.open(tag_file.c_str());
        if (tags.fail()) {
            cerr << "Error: Unable to open tag file for writing.\n";
            exit(1);
        }
        snps.open(snp_file.c_str());
        if (snps.fail()) {
            cerr << "Error: Unable to open SNPs file for writing.\n";
            exit(1);
        }
        alle.open(all_file.c_str());
        if (alle.fail()) {
            cerr << "Error: Unable to open allele file for writing.\n";
            exit(1);
        }
    }

    //
    // Record the version of Stacks used and the date generated as a comment in the catalog.
    //
    // Obtain the current date.
    //
    stringstream log;
    time_t       rawtime;
    struct tm   *timeinfo;
    char         date[32];
    time(&rawtime);
    timeinfo = localtime(&rawtime);
    strftime(date, 32, "%F %T", timeinfo);
    log << "# ustacks version " << VERSION << "; generated on " << date << "\n"; 
    if (gzip) {
        gzputs(gz_tags, log.str().c_str());
        gzputs(gz_snps, log.str().c_str());
        gzputs(gz_alle, log.str().c_str());
    } else {
        tags << log.str();
        snps << log.str();
        alle << log.str();
    }

    int id;

    char *buf = new char[m.begin()->second->len + 1];

    for (i = m.begin(); i != m.end(); i++) {
        float total = 0;
        tag_1 = i->second;

        //
        // Calculate the log likelihood of this merged stack.
        //
        tag_1->gen_matrix(u, r);
        tag_1->calc_likelihood();

        // First write the consensus sequence
        sstr << "0"              << "\t" 
             << sql_id           << "\t" 
             << tag_1->id        << "\t"
            //<< tag_1->cohort_id << "\t"
             << ""               << "\t" // chr
             << 0                << "\t" // bp
             << "+"              << "\t" // strand
             << "consensus\t"    << "\t" 
             << "\t" 
             << tag_1->con         << "\t" 
             << tag_1->deleveraged << "\t" 
             << tag_1->blacklisted << "\t"
             << tag_1->lumberjackstack << "\t"
             << tag_1->lnl << "\n";

        //
        // Write a sequence recording the output of the SNP model for each nucleotide.
        //
        sstr << "0" << "\t" 
             << sql_id << "\t" 
             << tag_1->id << "\t" 
            //<< "\t" // cohort_id
             << "\t"  // chr
             << "\t"  // bp
             << "\t"  // strand
             << "model\t" << "\t"
             << "\t";
        for (s = tag_1->snps.begin(); s != tag_1->snps.end(); s++) {
            switch((*s)->type) {
            case snp_type_het:
                sstr << "E";
                break;
            case snp_type_hom:
                sstr << "O";
                break;
            default:
                sstr << "U";
                break;
            }
        }
        sstr << "\t" 
             << "\t"
             << "\t"
             << "\t"
             << "\n";

        if (gzip) gzputs(gz_tags, sstr.str().c_str()); else tags << sstr.str();
        sstr.str("");

        //
        // Now write out the components of each unique tag merged into this locus.
        //
        id = 0;
        for (k = tag_1->utags.begin(); k != tag_1->utags.end(); k++) {
            tag_2  = u[*k];
            total += tag_2->count();

            for (uint j = 0; j < tag_2->map.size(); j++) {
                sstr << "0"       << "\t"
                     << sql_id    << "\t"
                     << tag_1->id << "\t"
                    //<< "\t" // cohort_id
                     << "\t" // chr
                     << "\t" // bp
                     << "\t" // strand
                     << "primary\t" 
                     << id << "\t" 
                     << seq_ids[tag_2->map[j]] << "\t" 
                     << tag_2->seq->seq(buf) 
                     << "\t\t\t\t\n";

                if (gzip) gzputs(gz_tags, sstr.str().c_str()); else tags << sstr.str();
                sstr.str("");
            }

            id++;
        }

        //
        // Write out the remainder tags merged into this unique tag.
        //
        for (k = tag_1->remtags.begin(); k != tag_1->remtags.end(); k++) {
            rem    = r[*k];
            total += rem->map.size();

            for (uint j = 0; j < rem->map.size(); j++)
                sstr << "0"       << "\t" 
                     << sql_id    << "\t" 
                     << tag_1->id << "\t"
                    //<< "\t" // cohort_id
                     << "\t" // chr
                     << "\t" // bp
                     << "\t" // strand
                     << "secondary\t"
                     << "\t" 
                     << seq_ids[rem->map[j]] << "\t" 
                     << rem->seq->seq(buf) 
                     << "\t\t\t\t\n";

            if (gzip) gzputs(gz_tags, sstr.str().c_str()); else tags << sstr.str();
            sstr.str("");
        }

        //
        // Write out the model calls for each nucleotide in this locus.
        //
        for (s = tag_1->snps.begin(); s != tag_1->snps.end(); s++) {
            sstr << "0"          << "\t" 
                 << sql_id       << "\t" 
                 << tag_1->id    << "\t" 
                 << (*s)->col    << "\t";

            switch((*s)->type) {
            case snp_type_het:
                sstr << "E\t";
                break;
            case snp_type_hom:
                sstr << "O\t";
                break;
            default:
                sstr << "U\t";
                break;
            }

            sstr << std::fixed   << std::setprecision(2)
                 << (*s)->lratio << "\t" 
                 << (*s)->rank_1 << "\t" 
                 << (*s)->rank_2 << "\t\t\n";
        }

        if (gzip) gzputs(gz_snps, sstr.str().c_str()); else snps << sstr.str();
        sstr.str("");

        //
        // Write the expressed alleles seen for the recorded SNPs and
        // the percentage of tags a particular allele occupies.
        //
        for (t = tag_1->alleles.begin(); t != tag_1->alleles.end(); t++) {
            sstr << "0"         << "\t" 
                 << sql_id      << "\t" 
                 << tag_1->id   << "\t" 
                 << (*t).first  << "\t" 
                 << (((*t).second/total) * 100) << "\t" 
                 << (*t).second << "\n";
        }
        if (gzip) gzputs(gz_alle, sstr.str().c_str()); else alle << sstr.str();
        sstr.str("");

    }

    if (gzip) {
        gzclose(gz_tags);
        gzclose(gz_snps);
        gzclose(gz_alle);
    } else {
        tags.close();
        snps.close();
        alle.close();
    }

    //
    // Free sequence IDs.
    //
    for (uint i = 0; i < seq_ids.size(); i++)
        delete [] seq_ids[i];

    //
    // If specified, output reads not utilized in any stacks.
    //
    if (retain_rem_reads) {
        string unused_file = out_path + in_file.substr(pos_1 + 1, (pos_2 - pos_1 - 1)) + ".unused.fa";

        gzFile   gz_unused;
        ofstream unused;

        if (gzip) {
            unused_file += ".gz";
            gz_unused = gzopen(unused_file.c_str(), "wb");
            if (!gz_unused) {
                cerr << "Error: Unable to open gzipped discard file '" << unused_file << "': " << strerror(errno) << ".\n";
                exit(1);
            }
            #if ZLIB_VERNUM >= 0x1240
            gzbuffer(gz_unused, libz_buffer_size);
            #endif
        } else {
            unused.open(unused_file.c_str());
            if (unused.fail()) {
                cerr << "Error: Unable to open discard file for writing.\n";
                exit(1);
            }
        }

        map<int, Rem *>::iterator r_it;
        for (r_it = r.begin(); r_it != r.end(); r_it++) {
            if (r_it->second->utilized == false)
                sstr << ">" << r_it->second->id << "\n" << r_it->second->seq->seq(buf) << "\n";
            if (gzip) gzputs(gz_unused, sstr.str().c_str()); else unused << sstr.str();
            sstr.str("");
        }

        if (gzip) gzclose(gz_unused); else unused.close();
    }

    delete [] buf;

    return 0;
}

int dump_stack_graph(string data_file, 
                     map<int, Stack *> &unique, 
                     map<int, MergedStack *> &merged, 
                     vector<int> &keys, 
                     map<int, map<int, double> > &dist_map, 
                     map<int, set<int> > &cluster_map) {
    uint s, t;
    double d, scale, scaled_d;
    char label[32];
    vector<string> colors;
    std::ofstream data(data_file.c_str());

    size_t pos_1 = data_file.find_last_of("/");
    size_t pos_2 = data_file.find_last_of(".");
    string title = data_file.substr(pos_1 + 1, pos_2 - pos_1 - 1);

    //
    // Output a list of IDs so we can locate these stacks in the final results.
    //
    for (s = 0; s < keys.size(); s++)
        data << "/* " << keys[s] << ": " << unique[merged[keys[s]]->utags[0]]->map[0] << "; depth: " << merged[keys[s]]->count << " */\n";

    //
    // Output a specification to visualize the stack graph using graphviz:
    //   http://www.graphviz.org/
    //
    data << "graph " << title.c_str() << " {\n"
         << "rankdir=LR\n"
         << "size=\"20!\"\n"
         << "overlap=false\n"
         << "node [shape=circle style=filled fillcolor=\"#3875d7\" fontname=\"Arial\"];\n"
         << "edge [fontsize=8.0 fontname=\"Arial\" color=\"#aaaaaa\"];\n";

    colors.push_back("red");
    colors.push_back("blue"); 
    colors.push_back("green"); 
    colors.push_back("brown"); 
    colors.push_back("purple");

    map<int, set<int> >::iterator c;
    set<int>::iterator it;
    int color_index = 0;
    string color;

    // Write out the clusters created by R, prior to writing all the nodes and connections.
    s = 0;
    for (c = cluster_map.begin(); c != cluster_map.end(); c++) {
        data << "subgraph " << s << " {\n"
             << "    edge [penwidth=5 fontsize=12.0 fontcolor=\"black\" color=\"black\"]\n";

        if ((*c).second.size() == 1) {
            color = "white";
            data << "    node [fillcolor=" << color.c_str() << " fontcolor=\"black\"]\n";
        } else {
            color = colors[color_index % colors.size()];
            data << "    node [fillcolor=" << color.c_str() << " fontcolor=\"white\"]\n";
            color_index++;
        }

        for (it = (*c).second.begin(); it != (*c).second.end(); it++) {
            data << "    " << *it << "\n";
        }

        if ((*c).second.size() > 1) {
            uint j = 0;
            for (it = (*c).second.begin(); it != (*c).second.end(); it++) {
                data << *it;
                if (j < (*c).second.size() - 1)
                    data << " -- ";
                j++;
            }
        }

        data << "}\n";
        s++;
    }

    //
    // Scale the graph to display on a 10 inch canvas. Find the largest edge weight
    // and scale the edge lengths to fit the canvas.
    //
    for (s = 0; s < keys.size(); s++)
        for (t = s+1; t < keys.size(); t++)
            scale = dist_map[keys[s]][keys[t]] > scale ? dist_map[keys[s]][keys[t]] : scale;
    scale = scale / 20;

    for (s = 0; s < keys.size(); s++) {
        for (t = s+1; t < keys.size(); t++) {
            d = dist_map[keys[s]][keys[t]];
            scaled_d = d / scale;
            scaled_d = scaled_d < 0.75 ? 0.75 : scaled_d;
            sprintf(label, "%.1f", d);
            data << keys[s] << " -- " << keys[t] << " [len=" << scaled_d << ", label=" << label << "];\n";
        }
    }

    data << "}\n";

    data.close();

    return 0;
}

int dump_unique_tags(map<int, Stack *> &u) {
    map<int, Stack *>::iterator it;
    vector<pair<int, int> >::iterator pit;
    vector<int>::iterator mit;
    char *c;

    for (it = u.begin(); it != u.end(); it++) {
        c = (*it).second->seq->seq();

        cerr << "UniqueTag UID: " << (*it).second->id << "\n"
             << "  Seq:       "   << c << "\n"
             << "  IDs:       "; 

        for (uint j = 0; j < it->second->map.size(); j++)
            cerr << it->second->map[j] << " ";

        cerr << "\n\n";

        delete [] c;
    }

    return 0;
}

int dump_merged_tags(map<int, MergedStack *> &m) {
    map<int, MergedStack *>::iterator it;
    vector<pair<int, int> >::iterator pit;
    vector<int>::iterator fit;

    for (it = m.begin(); it != m.end(); it++) {

        cerr << "MergedStack ID: " << it->second->id << "\n"
             << "  Consensus:  ";
        if (it->second->con != NULL)
            cerr << it->second->con << "\n";
        else 
            cerr << "\n";
        cerr << "  IDs:        "; 

        for (fit = it->second->utags.begin(); fit != it->second->utags.end(); fit++)
            cerr << (*fit) << " ";

        cerr << "\n"
             << "  Distances: ";

        for (pit = it->second->dist.begin(); pit != it->second->dist.end(); pit++)
            cerr << (*pit).first << ": " << (*pit).second << ", ";

        cerr << "\n\n";
    }

    return 0;
}

int load_radtags(string in_file, DNASeqHashMap &radtags, vector<DNANSeq *> &radtags_keys) {
    Input *fh = NULL;
    DNANSeq *d;

    if (in_file_type == FileT::fasta)
        fh = new Fasta(in_file.c_str());
    else if (in_file_type == FileT::fastq)
        fh = new Fastq(in_file.c_str());
    else if (in_file_type == FileT::gzfasta)
        fh = new GzFasta(in_file.c_str());
    else if (in_file_type == FileT::gzfastq)
        fh = new GzFastq(in_file.c_str());

    cerr << "Parsing " << in_file.c_str() << "\n";
    long  int corrected = 0;
    long  int i         = 0;
    short int seql      = 0;
    short int prev_seql = 0;
    bool  len_mismatch  = false;

    Seq c;
    c.id   = new char[id_len];
    c.seq  = new char[max_len];
    c.qual = new char[max_len];

    while ((fh->next_seq(c)) != 0) {
        if (i % 10000 == 0) cerr << "  Loading RAD-Tag " << i << "       \r";

        prev_seql = seql;
        seql      = 0;

        for (char *p = c.seq; *p != '\0'; p++, seql++)
            switch (*p) {
            case 'N':
            case 'n':
            case '.':
                *p = 'A';
                corrected++;
            }

        if (seql != prev_seql && prev_seql > 0) len_mismatch = true;

        d = new DNANSeq(seql, c.seq);

        pair<DNASeqHashMap::iterator, bool> r;

        r = radtags.insert(make_pair(d, HVal()));
        (*r.first).second.add_id(i);
        radtags_keys.push_back(d);
        i++;
    }
    cerr << "Loaded " << i << " RAD-Tags; inserted " << radtags.size() << " elements into the RAD-Tags hash map.\n";

    if (i == 0) {
        cerr << "Error: Unable to load data from '" << in_file.c_str() << "'.\n";
        exit(1);
    }

    cerr << "  " << corrected << " reads contained uncalled nucleotides that were modified.\n";

    if (len_mismatch)
        cerr << "  Warning: different sequence lengths detected, this will interfere with Stacks algorithms.\n";

    //
    // Close the file and delete the Input object.
    //
    delete fh;

    return 0;
}

int
load_seq_ids(vector<char *> &seq_ids)
{
    Input *fh = NULL;

    if (in_file_type == FileT::fasta)
        fh = new Fasta(in_file.c_str());
    else if (in_file_type == FileT::fastq)
        fh = new Fastq(in_file.c_str());
    else if (in_file_type == FileT::gzfasta)
        fh = new GzFasta(in_file.c_str());
    else if (in_file_type == FileT::gzfastq)
        fh = new GzFastq(in_file.c_str());

    cerr << "  Refetching sequencing IDs from " << in_file.c_str() << "... ";    

    char *id;
    Seq c;
    c.id   = new char[id_len];
    c.seq  = new char[max_len];
    c.qual = new char[max_len];

    while ((fh->next_seq(c)) != 0) {
        id = new char[strlen(c.id) + 1];
        strcpy(id, c.id);
        seq_ids.push_back(id);
    }
    cerr << "read " << seq_ids.size() << " sequence IDs.\n";

    delete fh;

    return 0;
}

int
calc_triggers(double cov_mean,
              double cov_stdev,
              double cov_scale,
              int &deleverage_trigger, int &removal_trigger)
{

    deleverage_trigger = (int) round(cov_mean + cov_stdev * cov_scale);
    removal_trigger    = (int) round(cov_mean + (cov_stdev * 2) * cov_scale);

    return 0;

//     //
//     // Calculate the deleverage trigger. Assume RAD-Tags are selected from
//     // the sample for sequencing randomly, forming a poisson distribution
//     // representing the depths of coverage of RAD-Tags in the sample. Calculate 
//     // the trigger value that is larger than the depth of coverage of 99.9999% of stacks.
//     //
//     long double lambda = cov_mean;
//     int k         = 0;
//     long double d = 0.0;
//     long double e = 0.0;
//     long double f = 0.0;
//     long double g = 0.0;
//     long double h = 0.0;
//     long double i = 0.0;

//     do {
//      e  = exp(-1 * lambda);
//      g  = pow(lambda, k);
//      f  = factorial(k);
//      h  = (e * g);
//      i  = h / f;
//      d += i;

//      //cerr << "iteration " << k << "; e: " << e << " h: " << h << " g: " << g << " F: " << f << " i: " << i << " D: " << d << "\n";
//      k++;
//     } while (d < 0.999999);

//     return k - 1;
}

long double factorial(int i) {
    long double f = 1;

    if (i == 0) return 1;

    do {
        f = f * i;
        i--;
    } while (i > 0);

    return f;
}

int parse_command_line(int argc, char* argv[]) {
    int c;
     
    while (1) {
        static struct option long_options[] = {
            {"help",             no_argument,       NULL, 'h'},
            {"version",          no_argument,       NULL, 'v'},
            {"infile_type",      required_argument, NULL, 't'},
            {"file",             required_argument, NULL, 'f'},
            {"outpath",          required_argument, NULL, 'o'},
            {"id",               required_argument, NULL, 'i'},
            {"min_cov",          required_argument, NULL, 'm'},
            {"max_dist",         required_argument, NULL, 'M'},
            {"max_sec_dist",     required_argument, NULL, 'N'},
            {"max_locus_stacks", required_argument, NULL, 'K'},
            {"k_len",            required_argument, NULL, 'k'},
            {"num_threads",      required_argument, NULL, 'p'},
            {"deleverage",       no_argument,       NULL, 'd'},
            {"remove_rep",       no_argument,       NULL, 'r'},
            {"retain_rem",       no_argument,       NULL, 'R'},
            {"graph",            no_argument,       NULL, 'g'},
            {"sec_hapl",         no_argument,       NULL, 'H'},
            {"gapped",           no_argument,       NULL, 'G'},
            {"max_gaps",         required_argument, NULL, 'X'},
            {"min_aln_len",      required_argument, NULL, 'x'},
            {"model_type",       required_argument, NULL, 'T'},
            {"bc_err_freq",      required_argument, NULL, 'e'},
            {"bound_low",        required_argument, NULL, 'L'},
            {"bound_high",       required_argument, NULL, 'U'},
            {"alpha",            required_argument, NULL, 'A'},
            {0, 0, 0, 0}
        };

        // getopt_long stores the option index here.
        int option_index = 0;
     
        c = getopt_long(argc, argv, "GhHvdrgRA:L:U:f:o:i:m:e:p:t:M:N:K:k:T:X:x:", long_options, &option_index);
     
        // Detect the end of the options.
        if (c == -1)
            break;
     
        switch (c) {
        case 'h':
            help();
            break;
        case 't':
            if (strcmp(optarg, "tsv") == 0)
                in_file_type = FileT::tsv;
            else if (strcmp(optarg, "fasta") == 0)
                in_file_type = FileT::fasta;
            else if (strcmp(optarg, "fastq") == 0)
                in_file_type = FileT::fastq;
            else if (strcasecmp(optarg, "gzfasta") == 0)
                in_file_type = FileT::gzfasta;
            else if (strcasecmp(optarg, "gzfastq") == 0)
                in_file_type = FileT::gzfastq;
            else
                in_file_type = FileT::unknown;
            break;
        case 'f':
            in_file = optarg;
            break;
        case 'o':
            out_path = optarg;
            break;
        case 'i':
            sql_id = is_integer(optarg);
            if (sql_id < 0) {
                cerr << "SQL ID (-i) must be an integer, e.g. 1, 2, 3\n";
                help();
            }
            break;
        case 'm':
            min_merge_cov = is_integer(optarg);
            break;
        case 'M':
            max_utag_dist = is_integer(optarg);
            break;
        case 'N':
            max_rem_dist = is_integer(optarg);
            break;
        case 'd':
            deleverage_stacks++;
            break;
        case 'r':
            remove_rep_stacks++;
            break;
        case 'K':
            max_subgraph = is_integer(optarg);
            break;
        case 'k':
            set_kmer_len = false;
            kmer_len     = is_integer(optarg);
            break;
        case 'R':
            retain_rem_reads = true;
            break;
        case 'g':
            dump_graph++;
            break;
        case 'G':
            gapped_alignments = true;
            break;
        case 'X':
            max_gaps = is_integer(optarg);
            break;
        case 'x':
            min_match_len = is_double(optarg);
            break;
        case 'T':
            if (strcmp(optarg, "snp") == 0) {
                model_type = snp;
            } else if (strcmp(optarg, "fixed") == 0) {
                model_type = fixed;
            } else if (strcmp(optarg, "bounded") == 0) {
                model_type = bounded;
            } else {
                cerr << "Unknown model type specified '" << optarg << "'\n";
                help();
            }
        case 'e':
            barcode_err_freq = is_double(optarg);
            break;
        case 'L':
            bound_low  = is_double(optarg);
            break;
        case 'U':
            bound_high = is_double(optarg);
            break;
        case 'A':
            alpha = is_double(optarg);
            break;
        case 'H':
            call_sec_hapl = false;
            break;
        case 'p':
            num_threads = is_integer(optarg);
            break;
        case 'v':
            version();
            break;
        case '?':
            // getopt_long already printed an error message.
            help();
            break;
     
        default:
            cerr << "Unknown command line option '" << (char) c << "'\n";
            help();
            abort();
        }
    }

    if (set_kmer_len == false && (kmer_len < 5 || kmer_len > 31)) {
        cerr << "Kmer length must be between 5 and 31bp.\n";
        help();
    }

    if (alpha != 0.1 && alpha != 0.05 && alpha != 0.01 && alpha != 0.001) {
        cerr << "SNP model alpha significance level must be either 0.1, 0.05, 0.01, or 0.001.\n";
        help();
    }

    if (bound_low != 0 && (bound_low < 0 || bound_low >= 1.0)) {
        cerr << "SNP model lower bound must be between 0.0 and 1.0.\n";
        help();
    }

    if (bound_high != 1 && (bound_high <= 0 || bound_high > 1.0)) {
        cerr << "SNP model upper bound must be between 0.0 and 1.0.\n";
        help();
    }

    if (bound_low > 0 || bound_high < 1.0) {
        model_type = bounded;
    }

    if (in_file.length() == 0 || in_file_type == FileT::unknown) {
        cerr << "You must specify an input file of a supported type.\n";
        help();
    }

    if (out_path.length() == 0) 
        out_path = ".";

    if (out_path.at(out_path.length() - 1) != '/') 
        out_path += "/";

    if (model_type == fixed && barcode_err_freq == 0) {
        cerr << "You must specify the barcode error frequency.\n";
        help();
    }

    return 0;
}

void version() {
    std::cerr << "ustacks " << VERSION << "\n\n";

    exit(0);
}

void help() {
    std::cerr << "ustacks " << VERSION << "\n"
              << "ustacks -t file_type -f file_path [-d] [-r] [-o path] [-i id] [-m min_cov] [-M max_dist] [-p num_threads] [-R] [-H] [-h]" << "\n"
              << "  t: input file Type. Supported types: fasta, fastq, gzfasta, or gzfastq.\n"
              << "  f: input file path.\n"
              << "  o: output path to write results.\n"
              << "  i: SQL ID to insert into the output to identify this sample.\n"
              << "  m: Minimum depth of coverage required to create a stack (default 3).\n"
              << "  M: Maximum distance (in nucleotides) allowed between stacks (default 2).\n"
              << "  N: Maximum distance allowed to align secondary reads to primary stacks (default: M + 2).\n"
              << "  R: retain unused reads.\n"
              << "  H: disable calling haplotypes from secondary reads.\n"
              << "  p: enable parallel execution with num_threads threads.\n"
              << "  h: display this help messsage.\n\n"
              << "  Stack assembly options:\n"
              << "    r: enable the Removal algorithm, to drop highly-repetitive stacks (and nearby errors) from the algorithm.\n"
              << "    d: enable the Deleveraging algorithm, used for resolving over merged tags.\n"
              << "    --max_locus_stacks <num>: maximum number of stacks at a single de novo locus (default 3).\n"
              << "     --k_len <len>: specify k-mer size for matching between alleles and loci (automatically calculated by default).\n\n"
              << "  Gapped assembly options:\n"
              << "    --gapped: preform gapped alignments between stacks.\n"
              << "    --max_gaps: number of gaps allowed between stacks before merging (default: 2).\n"
              << "    --min_aln_len: minimum length of aligned sequence in a gapped alignment (default: 0.80).\n\n"
              << "  Model options:\n" 
              << "    --model_type: either 'snp' (default), 'bounded', or 'fixed'\n"
              << "    For the SNP or Bounded SNP model:\n"
              << "      --alpha <num>: chi square significance level required to call a heterozygote or homozygote, either 0.1, 0.05 (default), 0.01, or 0.001.\n"
              << "    For the Bounded SNP model:\n"
              << "      --bound_low <num>: lower bound for epsilon, the error rate, between 0 and 1.0 (default 0).\n"
              << "      --bound_high <num>: upper bound for epsilon, the error rate, between 0 and 1.0 (default 1).\n"
              << "    For the Fixed model:\n"
              << "      --bc_err_freq <num>: specify the barcode error frequency, between 0 and 1.0.\n";

    exit(0);
}
