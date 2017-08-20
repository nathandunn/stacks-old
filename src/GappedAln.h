// -*-mode:c++; c-style:k&r; c-basic-offset:4;-*-
//
// Copyright 2016 - 2017, Julian Catchen <jcatchen@illinois.edu>
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

#ifndef __GAPPEDALN_H__
#define __GAPPEDALN_H__

#include "SuffixTree.h"

enum dynprog {dynp_down, dynp_right, dynp_diag};

class AlignRes {
public:
    string cigar;
    uint   gap_cnt;
    uint   contiguity;
    double pct_id;
    AlignRes() {
        this->gap_cnt    = 0;
        this->contiguity = 0;
        this->pct_id     = 0.0;
    }
    AlignRes(string&& cigar, uint gapcnt, uint contiguity, double pct_id) {
        this->cigar      = move(cigar);
        this->gap_cnt    = gapcnt;
        this->contiguity = contiguity;
        this->pct_id     = pct_id;
    }
};

class AlignPath {
public:
    bool diag;
    bool left;
    bool up;

    AlignPath() {
        diag = false;
        left = false;
        up   = false;
    }
    int count() {
        int cnt;
        cnt  = this->up   ? 1 : 0;
        cnt += this->left ? 1 : 0;
        cnt += this->diag ? 1 : 0;

        return cnt;
    }
};

//
// Needleman-Wunsch Alignment
//
const double gapopen_score  = -10;
const double gapext_score   = -0.5;
const double mismatch_score = -4;
const double match_score    =  5;

class GappedAln {
    uint        _m;
    uint        _n;
    uint        _m_size;
    uint        _n_size;
    double    **matrix;
    AlignPath **path;
    AlignRes    _aln;

    inline int swap(double *, dynprog *, int, int);
    int        trace_alignment(const string&, const string&);
    int        trace_constrained_alignment(const string&, const string&, int, int, int, int);
 public:
    GappedAln();
    GappedAln(int i) : GappedAln(i, i) {};
    GappedAln(int, int);
    ~GappedAln();

    int init(int, int);
    int align(const string&, const string&);
    int align_constrained(const string&, const string&, const vector<STAln> &);
    const AlignRes& result();

    int parse_cigar(vector<pair<char, uint> > &);
    int dump_alignment(const string&, const string&);
};

GappedAln::GappedAln()
{
    this->_m      = 0;
    this->_n      = 0;
    this->_m_size = this->_m;
    this->_n_size = this->_n;
    this->matrix  = NULL;
    this->path    = NULL;
}

GappedAln::GappedAln(int len_1, int len_2)
{
    this->_m      = len_1 + 1;
    this->_n      = len_2 + 1;
    this->_m_size = this->_m;
    this->_n_size = this->_n;

    this->matrix = new double * [this->_m];
    for (uint i = 0; i < this->_m; i++) {
        this->matrix[i] = new double [this->_n];
        memset(this->matrix[i], 0, sizeof(double) * this->_n);
    }

    this->path = new AlignPath * [this->_m];
    for (uint i = 0; i < this->_m; i++)
        this->path[i] = new AlignPath [this->_n];
}

GappedAln::~GappedAln()
{
    for (uint i = 0; i < this->_m_size; i++) {
        delete [] this->matrix[i];
        delete [] this->path[i];
    }
    delete [] this->matrix;
    delete [] this->path;
}

int
GappedAln::init(int size_1, int size_2)
{
    //
    // Resize the underlying matrix and path arrays, if necessary.
    //
    if ((size_1 + 1) > (int)this->_m_size || (size_2 + 1) > (int)this->_n_size) {
        for (uint i = 0; i < this->_m_size; i++) {
            delete [] this->matrix[i];
            delete [] this->path[i];
        }
        if (this->_m_size > 0) {
            delete [] this->matrix;
            delete [] this->path;
        }

        this->_m_size = size_1 + 1;
        this->_n_size = size_2 + 1;
        //
        // Resize the arrays to be 25% larger than the requested size.
        //
        int new_size = this->_m_size > this->_n_size ? this->_m_size : this->_n_size;
        new_size += int((double) new_size * 0.25);
        this->_m_size = new_size;
        this->_n_size = new_size;

        this->matrix = new double * [this->_m_size];
        for (uint i = 0; i < this->_m_size; i++) {
            this->matrix[i] = new double [this->_n_size];
            memset(this->matrix[i], 0, this->_n_size);
        }

        this->path = new AlignPath * [this->_m_size];
        for (uint i = 0; i < this->_m_size; i++)
            this->path[i] = new AlignPath [this->_n_size];
    }

    //
    // Otherwise, set the dimensions of the matrix and path arrays to be the sequence lengths.
    //
    this->_m = size_1 + 1;
    this->_n = size_2 + 1;

    return 0;
}

int
GappedAln::align(const string& tag_1, const string& tag_2)
{
    //         j---->
    //        [0][1][2][3]...[n-1]
    //       +--------------------
    // i [0] | [i][j]
    // | [1] |
    // | [2] |
    // v [3] |
    //   ... |
    // [m-1] |
    //

    //
    // Initialize the first column and row of the dynamic programming
    // matrix and the path array.
    //
    path[0][0].diag = false;
    path[0][0].up   = false;
    path[0][0].left = false;
    matrix[0][0]    = 0.0;
    for (uint i = 1; i < this->_m; i++) {
        this->matrix[i][0]    = this->path[i - 1][0].up ? this->matrix[i - 1][0] + gapext_score : this->matrix[i - 1][0] + gapopen_score;
        this->path[i][0].diag = false;
        this->path[i][0].up   = true;
        this->path[i][0].left = false;
    }
    for (uint j = 1; j < this->_n; j++) {
        this->matrix[0][j]    = this->path[0][j - 1].left ? this->matrix[0][j - 1] + gapext_score : this->matrix[0][j - 1] + gapopen_score;
        this->path[0][j].diag = false;
        this->path[0][j].up   = false;
        this->path[0][j].left = true;
    }

    double  score_down, score_diag, score_right;
    double  scores[3];
    dynprog direction[3];

    for (uint i = 1; i < this->_m; i++) {
        for (uint j = 1; j < this->_n; j++) {
            // Calculate the score:
            //   1) If we were to move down from the above cell.
            score_down   = this->matrix[i - 1][j];
            score_down  += this->path[i - 1][j].up ?  gapext_score : gapopen_score;
            //   2) If we were to move diagonally from the above and left cell.
            score_diag   = this->matrix[i - 1][j - 1] + (tag_1[i - 1] == tag_2[j - 1] ? match_score : mismatch_score);
            //   3) If we were to move over from the cell left of us.
            score_right  = this->matrix[i][j - 1];
            score_right += this->path[i][j - 1].left ? gapext_score : gapopen_score;

            //
            // Sort the scores, highest to lowest.
            //
            scores[0]    = score_down;
            direction[0] = dynp_down;
            scores[1]    = score_diag;
            direction[1] = dynp_diag;
            scores[2]    = score_right;
            direction[2] = dynp_right;

            if (scores[0] < scores[1])
                this->swap(scores, direction, 0, 1);
            if (scores[1] < scores[2])
                this->swap(scores, direction, 1, 2);
            if (scores[0] < scores[1])
                this->swap(scores, direction, 0, 1);

            this->matrix[i][j] = scores[0];

            if (scores[0] > scores[1]) {
                //
                // One path is best.
                //
                switch (direction[0]) {
                case dynp_diag:
                    this->path[i][j].diag = true;
                    this->path[i][j].up   = false;
                    this->path[i][j].left = false;
                    break;
                case dynp_down:
                    this->path[i][j].diag = false;
                    this->path[i][j].up   = true;
                    this->path[i][j].left = false;
                    break;
                case dynp_right:
                default:
                    this->path[i][j].diag = false;
                    this->path[i][j].up   = false;
                    this->path[i][j].left = true;
                }

            } else if (scores[0] == scores[1]) {
                //
                // Two of the paths are equivalent.
                //
                switch (direction[0]) {
                case dynp_diag:
                    this->path[i][j].diag = true;

                    switch (direction[1]) {
                    case dynp_down:
                        this->path[i][j].up   = true;
                        this->path[i][j].left = false;
                        break;
                    default:
                    case dynp_right:
                        this->path[i][j].up   = false;
                        this->path[i][j].left = true;
                        break;
                    }
                    break;
                case dynp_down:
                    this->path[i][j].up = true;

                    switch (direction[1]) {
                    case dynp_right:
                        this->path[i][j].diag  = false;
                        this->path[i][j].left = true;
                        break;
                    default:
                    case dynp_diag:
                        this->path[i][j].diag  = true;
                        this->path[i][j].left = false;
                        break;
                    }
                    break;
                default:
                case dynp_right:
                    this->path[i][j].left = true;

                    switch (direction[1]) {
                    case dynp_diag:
                        this->path[i][j].diag = true;
                        this->path[i][j].up   = false;
                        break;
                    default:
                    case dynp_down:
                        this->path[i][j].diag = false;
                        this->path[i][j].up   = true;
                        break;
                    }
                    break;
                }

            } else {
                //
                // All paths equivalent.
                //
                this->path[i][j].diag = true;
                this->path[i][j].up   = true;
                this->path[i][j].left = true;
            }
        }
    }

    // dump_alignment(tag_1, tag_2, matrix, path);

    if (this->trace_alignment(tag_1, tag_2))
        return 1;

    return 0;
}

int
GappedAln::align_constrained(const string& query, const string& subj, const vector<STAln> &alns)
{
    //         j----> subject
    //        [0][1][2][3]...[n-1]
    //       +--------------------
    // i [0] | [i][j]
    // | [1] |
    // | [2] |
    // v [3] |
    //   ... |
    // [m-1] |
    // query

    //
    // Initialize the first column and row of the dynamic programming
    // matrix and the path array.
    //
    path[0][0].diag = false;
    path[0][0].up   = false;
    path[0][0].left = false;
    matrix[0][0]    = 0.0;
    // for (uint i = 1; i < this->_m; i++) {
    //     this->matrix[i][0]    = this->path[i - 1][0].up ? this->matrix[i - 1][0] + gapext_score : this->matrix[i - 1][0] + gapopen_score;
    //     this->path[i][0].diag = false;
    //     this->path[i][0].up   = true;
    //     this->path[i][0].left = false;
    // }
    // for (uint j = 1; j < this->_n; j++) {
    //     this->matrix[0][j]    = this->path[0][j - 1].left ? this->matrix[0][j - 1] + gapext_score : this->matrix[0][j - 1] + gapopen_score;
    //     this->path[0][j].diag = false;
    //     this->path[0][j].up   = false;
    //     this->path[0][j].left = true;
    // }

    double  score_down, score_diag, score_right;
    double  scores[3];
    dynprog direction[3];
    uint    q, s;
    uint    q_min = this->_m;
    uint    q_max = 0;
    uint    s_min = this->_n;
    uint    s_max = 0;

    //
    // Fill in each pre-aligned region of the sequences starting from the first fragment as ordered by
    // the query sequence.
    // Next fill in the following connector sequence span the region between two pre-aligned regions.
    //
    for (uint n = 0; n < alns.size(); n++) {
        q_min = alns[n].query_pos < q_min ? alns[n].query_pos : q_min;
        s_min = alns[n].subj_pos  < s_min ? alns[n].subj_pos  : s_min;
        q_max = alns[n].query_pos + alns[n].aln_len > q_max ? alns[n].query_pos + alns[n].aln_len : q_max;
        s_max = alns[n].subj_pos  + alns[n].aln_len > s_max ? alns[n].subj_pos  + alns[n].aln_len : s_max;
 
        for (uint i = alns[n].query_pos + 1; i <= alns[n].query_pos + alns[n].aln_len; i++) {
            for (uint j = alns[n].subj_pos + 1; j <= alns[n].subj_pos + alns[n].aln_len; j++) {
                q = i - alns[n].query_pos;
                s = j - alns[n].subj_pos;
                if (q == s) {
                    this->matrix[i][j]    = 1000;
                    this->path[i][j].diag = true;
                    this->path[i][j].up   = false;
                    this->path[i][j].left = false;
                } else if (q < s) {
                    this->matrix[i][j]    = -1000;
                    this->path[i][j].diag = false;
                    this->path[i][j].up   = false;
                    this->path[i][j].left = true;
                } else {
                    this->matrix[i][j]    = -1000;
                    this->path[i][j].diag = false;
                    this->path[i][j].up   = true;
                    this->path[i][j].left = false;
                }                
            }
        }
    }

    q_min++;
    q_max++;
    s_min++;
    s_max++;

    //
    // Now, fill in the regions of the matrix between the predefined matching regions.
    //
    int q_start, q_end, s_start, s_end;

    for (uint m = 0, n = 1; m < alns.size() - 1; m++, n++) {
        q_start = alns[n].query_pos - (alns[m].query_pos + alns[m].aln_len);

        if (q_start <= 0)
            continue;

        q_start = 1 + alns[m].query_pos + alns[m].aln_len;
        q_end   = 1 + alns[n].query_pos - 1;
        s_start = 1 + alns[m].subj_pos + alns[m].aln_len;
        s_end   = 1 + alns[n].subj_pos - 1;

        cerr << "q_start: " << q_start << ", q_end: " << q_end << "; s_start: " << s_start << ", s_end: " << s_end << "\n";
        //
        // Bound the region we are about to score.
        //
        for (int j = s_start; j <= s_end + 1; j++) {
            this->matrix[q_start - 1][j]    = -1000;
            this->path[q_start - 1][j].diag = false;
            this->path[q_start - 1][j].up   = false;
            this->path[q_start - 1][j].left = true;
        }
        for (int j = s_start - 1; j <= s_end; j++) {
            this->matrix[q_end + 1][j]    = -1000;
            this->path[q_end + 1][j].diag = false;
            this->path[q_end + 1][j].up   = true;
            this->path[q_end + 1][j].left = false;
        }
        for (int i = q_start; i <= q_end + 1; i++) {
            this->matrix[i][s_start - 1]    = -1000;
            this->path[i][s_start - 1].diag = false;
            this->path[i][s_start - 1].up   = true;
            this->path[i][s_start - 1].left = false;
        }
        for (int i = q_start - 1; i <= q_end; i++) {
            this->matrix[i][s_end + 1]    = -1000;
            this->path[i][s_end + 1].diag = false;
            this->path[i][s_end + 1].up   = false;
            this->path[i][s_end + 1].left = true;
        }

        //
        // Score the bounded region.
        //
        for (int i = q_start; i <= q_end; i++) {
            for (int j = s_start; j <= s_end; j++) {
                // Calculate the score:
                //   1) If we were to move down from the above cell.
                score_down   = this->matrix[i - 1][j];
                score_down  += this->path[i - 1][j].up ?  gapext_score : gapopen_score;
                //   2) If we were to move diagonally from the above and left cell.
                score_diag   = this->matrix[i - 1][j - 1] + (query[i - 1] == subj[j - 1] ? match_score : mismatch_score);
                //   3) If we were to move over from the cell left of us.
                score_right  = this->matrix[i][j - 1];
                score_right += this->path[i][j - 1].left ? gapext_score : gapopen_score;

                //
                // Sort the scores, highest to lowest.
                //
                scores[0]    = score_down;
                direction[0] = dynp_down;
                scores[1]    = score_diag;
                direction[1] = dynp_diag;
                scores[2]    = score_right;
                direction[2] = dynp_right;

                if (scores[0] < scores[1])
                    this->swap(scores, direction, 0, 1);
                if (scores[1] < scores[2])
                    this->swap(scores, direction, 1, 2);
                if (scores[0] < scores[1])
                    this->swap(scores, direction, 0, 1);

                this->matrix[i][j] = scores[0];

                if (scores[0] > scores[1]) {
                    //
                    // One path is best.
                    //
                    switch (direction[0]) {
                    case dynp_diag:
                        this->path[i][j].diag = true;
                        this->path[i][j].up   = false;
                        this->path[i][j].left = false;
                        break;
                    case dynp_down:
                        this->path[i][j].diag = false;
                        this->path[i][j].up   = true;
                        this->path[i][j].left = false;
                        break;
                    case dynp_right:
                    default:
                        this->path[i][j].diag = false;
                        this->path[i][j].up   = false;
                        this->path[i][j].left = true;
                    }

                } else if (scores[0] == scores[1]) {
                    //
                    // Two of the paths are equivalent.
                    //
                    switch (direction[0]) {
                    case dynp_diag:
                        this->path[i][j].diag = true;

                        switch (direction[1]) {
                        case dynp_down:
                            this->path[i][j].up   = true;
                            this->path[i][j].left = false;
                            break;
                        default:
                        case dynp_right:
                            this->path[i][j].up   = false;
                            this->path[i][j].left = true;
                            break;
                        }
                        break;
                    case dynp_down:
                        this->path[i][j].up = true;

                        switch (direction[1]) {
                        case dynp_right:
                            this->path[i][j].diag  = false;
                            this->path[i][j].left = true;
                            break;
                        default:
                        case dynp_diag:
                            this->path[i][j].diag  = true;
                            this->path[i][j].left = false;
                            break;
                        }
                        break;
                    default:
                    case dynp_right:
                        this->path[i][j].left = true;

                        switch (direction[1]) {
                        case dynp_diag:
                            this->path[i][j].diag = true;
                            this->path[i][j].up   = false;
                            break;
                        default:
                        case dynp_down:
                            this->path[i][j].diag = false;
                            this->path[i][j].up   = true;
                            break;
                        }
                        break;
                    }

                } else {
                    //
                    // All paths equivalent.
                    //
                    this->path[i][j].diag = true;
                    this->path[i][j].up   = true;
                    this->path[i][j].left = true;
                }
            }
        }
    }

    this->dump_alignment(query, subj);

    if (this->trace_constrained_alignment(query, subj, q_min, q_max, s_min, s_max))
        return 1;

    return 0;
}

inline int
GappedAln::swap(double *scores, dynprog *direction, int index_1, int index_2)
{
    double swap        = scores[index_1];
    scores[index_1]    = scores[index_2];
    scores[index_2]    = swap;
    dynprog swapdir    = direction[index_1];
    direction[index_1] = direction[index_2];
    direction[index_2] = swapdir;

    return 0;
}

bool
compare_alignres(const AlignRes& a, const AlignRes& b)
{
    if (a.gap_cnt == b.gap_cnt) {

        if (a.pct_id == b.pct_id)
            return (a.contiguity > b.contiguity);
        else
            return (a.pct_id > b.pct_id);

    } else {
        return (a.gap_cnt < b.gap_cnt);
    }
}

int
GappedAln::trace_alignment(const string& tag_1, const string& tag_2)
{
    //         j---->
    //        [0][1][2][3]...[n-1]
    //       +--------------------
    // i [0] | [i][j]
    // | [1] |
    // | [2] |
    // v [3] |
    //   ... |
    // [m-1] |
    //
    int    i, j, cnt, len, gaps, contiguity;
    double ident;
    string cigar;
    char   buf[id_len];

    vector<AlignRes> alns;
    bool more_paths = true;
    bool seq_break  = false;

    do {
        more_paths = false;

        i = this->_m - 1;
        j = this->_n - 1;

        string aln_1, aln_2;

        while (i > 0 || j > 0) {
            cnt  = this->path[i][j].count();

            if (cnt > 1) more_paths = true;

            if (this->path[i][j].diag) {
                aln_1 += tag_1[i - 1];
                aln_2 += tag_2[j - 1];
                if (cnt > 1) this->path[i][j].diag = false;
                i--;
                j--;
            } else if (this->path[i][j].up) {
                aln_1 += tag_1[i - 1];
                aln_2 += "-";
                if (cnt > 1) this->path[i][j].up = false;
                i--;
            } else if (this->path[i][j].left) {
                aln_1 += "-";
                aln_2 += tag_2[j - 1];
                if (cnt > 1) this->path[i][j].left = false;
                j--;
            }
        }

        reverse(aln_1.begin(), aln_1.end());
        reverse(aln_2.begin(), aln_2.end());

        //
        // Convert to CIGAR strings.
        //
        cigar      = "";
        len        = aln_1.length();
        gaps       = 0;
        contiguity = 0;
        seq_break  = false;
        ident      = 0.0;
        i          = 0;
        while (i < len) {
            if (aln_1[i] != '-' && aln_2[i] != '-') {
                cnt = 0;
                do {
                    if (aln_1[i] == aln_2[i]) ident++;
                    cnt++;
                    i++;
                    if (seq_break == false) contiguity++;
                } while (i < len && aln_1[i] != '-' && aln_2[i] != '-');
                sprintf(buf, "%dM", cnt);

            } else if (aln_1[i] == '-') {
                cnt = 0;
                do {
                    cnt++;
                    i++;
                } while (i < len && aln_1[i] == '-');
                sprintf(buf, "%dD", cnt);
                gaps++;
                seq_break = true;

            } else {
                cnt = 0;
                do {
                    cnt++;
                    i++;
                } while (i < len && aln_2[i] == '-');
                sprintf(buf, "%dI", cnt);
                gaps++;
                seq_break = true;
            }

            cigar += buf;
        }

        alns.push_back(AlignRes(move(cigar), gaps, contiguity, (ident / (double) len)));

        // cerr << aln_1 << " [" << cigar << ", contiguity: " << contiguity << ", gaps: " << gaps << "]\n"
        //      << aln_2 << "\n";

    } while (more_paths);

    sort(alns.begin(), alns.end(), compare_alignres);
    this->_aln = alns[0];
    // cerr << "Final alignment: " << this->_aln.cigar << "; contiguity: " << contiguity << "; gaps: " << this->_aln.gap_cnt << "\n";

    return 1;
}

int
GappedAln::trace_constrained_alignment(const string& query, const string& subj, int q_start, int q_end, int s_start, int s_end)
{
    //         j---->      subject
    //        [0][1][2][3]...[n-1]
    //       +--------------------
    // i [0] | [i][j]
    // | [1] |
    // | [2] |
    // v [3] |
    //   ... |
    // [m-1] |
    // query
    //
    int    i, j, cnt, len, gaps, contiguity;
    double ident;
    string cigar;
    char   buf[id_len];

    vector<AlignRes> alns;
    bool more_paths = true;
    bool seq_break  = false;

    do {
        more_paths = false;

        i = q_end - 1;
        j = s_end - 1;

        string aln_1, aln_2;

        while (i >= q_start || j >= s_start) {
            cnt  = this->path[i][j].count();

            if (cnt > 1) more_paths = true;

            if (this->path[i][j].diag) {
                aln_1 += query[i - 1];
                aln_2 += subj[j - 1];
                if (cnt > 1) this->path[i][j].diag = false;
                i--;
                j--;
            } else if (this->path[i][j].up) {
                aln_1 += query[i - 1];
                aln_2 += "-";
                if (cnt > 1) this->path[i][j].up = false;
                i--;
            } else if (this->path[i][j].left) {
                aln_1 += "-";
                aln_2 += subj[j - 1];
                if (cnt > 1) this->path[i][j].left = false;
                j--;
            }
        }

        reverse(aln_1.begin(), aln_1.end());
        reverse(aln_2.begin(), aln_2.end());

        //
        // Convert to CIGAR strings.
        //
        cigar      = "";
        len        = aln_1.length();
        gaps       = 0;
        contiguity = 0;
        seq_break  = false;
        ident      = 0.0;
        i          = 0;
        while (i < len) {
            if (aln_1[i] != '-' && aln_2[i] != '-') {
                cnt = 0;
                do {
                    if (aln_1[i] == aln_2[i]) ident++;
                    cnt++;
                    i++;
                    if (seq_break == false) contiguity++;
                } while (i < len && aln_1[i] != '-' && aln_2[i] != '-');
                sprintf(buf, "%dM", cnt);

            } else if (aln_1[i] == '-') {
                cnt = 0;
                do {
                    cnt++;
                    i++;
                } while (i < len && aln_1[i] == '-');
                sprintf(buf, "%dD", cnt);
                gaps++;
                seq_break = true;

            } else {
                cnt = 0;
                do {
                    cnt++;
                    i++;
                } while (i < len && aln_2[i] == '-');
                sprintf(buf, "%dI", cnt);
                gaps++;
                seq_break = true;
            }

            cigar += buf;
        }

        alns.push_back(AlignRes(move(cigar), gaps, contiguity, (ident / (double) len)));

        // cerr << aln_1 << " [" << cigar << ", contiguity: " << contiguity << ", gaps: " << gaps << "]\n"
        //      << aln_2 << "\n";

    } while (more_paths);

    sort(alns.begin(), alns.end(), compare_alignres);
    this->_aln = alns[0];
    // cerr << "Final alignment: " << this->_aln.cigar << "; contiguity: " << contiguity << "; gaps: " << this->_aln.gap_cnt << "\n";

    return 1;
}

const AlignRes&
GappedAln::result()
{
    return this->_aln;
}

int
GappedAln::parse_cigar(vector<pair<char, uint> > &cigar)
{
    char buf[id_len];
    int  dist;
    const char *p, *q;

    p = this->_aln.cigar.c_str();

    cigar.clear();

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
GappedAln::dump_alignment(const string& tag_1, const string& tag_2)
{
    //         j---->
    //        [0][1][2][3]...[n-1]
    //       +--------------------
    // i [0] | [i][j]
    // | [1] |
    // | [2] |
    // v [3] |
    //   ... |
    // [m-1] |
    //

    //
    // Output the score matrix.
    //
    cout << "         ";
    for (uint j = 0; j < this->_n; j++)
        cout << "   " << tag_2[j] << "  |";
    cout << "\n";

    cout << "  ";
    for (uint j = 0; j < this->_n; j++)
        printf("% 6.1f|", this->matrix[0][j]);
    cout << "\n";

    for (uint i = 1; i < this->_m; i++) {
        cout << tag_1[i - 1] << " ";
        for (uint j = 0; j < this->_n; j++)
            printf("% 6.1f|", this->matrix[i][j]);
        cout << "\n";
    }

    cout << "\n";

    //
    // Output the path matrix.
    //
    cout << "       ";
    for (uint j = 0; j < this->_n; j++)
        cout << "  " << tag_2[j] << " |";
    cout << "\n";

    cout << "  ";
    for (uint j = 0; j < this->_n; j++) {
        cout << " ";
        this->path[0][j].diag ? cout << "d" : cout << " ";
        this->path[0][j].up   ? cout << "u" : cout << " ";
        this->path[0][j].left ? cout << "l" : cout << " ";
        cout << "|";
    }
    cout << "\n";

    for (uint i = 1; i < this->_m; i++) {
        cout << tag_1[i - 1] << " ";
        for (uint j = 0; j < this->_n; j++) {
            cout << " ";
            this->path[i][j].diag ? cout << "d" : cout << " ";
            this->path[i][j].up   ? cout << "u" : cout << " ";
            this->path[i][j].left ? cout << "l" : cout << " ";
            cout << "|";
        }
        cout << "\n";
    }

    cout << "\n";

    return 0;
}

#endif // __GAPPEDALN_H__
