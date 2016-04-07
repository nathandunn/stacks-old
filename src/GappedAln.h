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

#ifndef __GAPPEDALN_H__
#define __GAPPEDALN_H__

enum dynprog {dynp_down, dynp_right, dynp_diag};

class AlignRes {
public:
    string cigar;
    uint   gap_cnt;
    double pct_id;
    AlignRes() {
	this->gap_cnt = 0;
	this->pct_id  = 0.0;
    }
    AlignRes(string c, uint g, double p) {
	this->cigar   = c;
	this->gap_cnt = g;
	this->pct_id  = p;
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
    double    **matrix;
    AlignPath **path;
    AlignRes    _aln;

 public:    
    GappedAln(int);
    GappedAln(int, int);
    ~GappedAln();

    int        align(string, string);
    inline int swap(double *, dynprog *, int, int);
    int        trace_alignment(string, string);
    AlignRes   result();

    int parse_cigar(vector<pair<char, uint> > &);
    int dump_alignment(string, string);
};


GappedAln::GappedAln(int len)
{
    this->_m = len + 1;
    this->_n = len + 1;

    this->matrix = new double * [this->_m];
    for (uint i = 0; i < this->_m; i++)
        this->matrix[i] = new double [this->_n];

    this->path = new AlignPath * [this->_m];
    for (uint i = 0; i < this->_m; i++)
        this->path[i] = new AlignPath [this->_n];
}

GappedAln::GappedAln(int len_1, int len_2)
{
    this->_m = len_1 + 1;
    this->_n = len_2 + 1;

    this->matrix = new double * [this->_m];
    for (uint i = 0; i < this->_m; i++)
        this->matrix[i] = new double [this->_n];

    this->path = new AlignPath * [this->_m];
    for (uint i = 0; i < this->_m; i++)
        this->path[i] = new AlignPath [this->_n];
}

GappedAln::~GappedAln()
{
    for (int i = 0; i < this->_m; i++) {
        delete [] this->matrix[i];
        delete [] this->path[i];
    }
    delete [] this->matrix;
    delete [] this->path;
}

int
GappedAln::align(string tag_1, string tag_2)
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
compare_alignres(AlignRes a, AlignRes b)
{
    if (a.gap_cnt == b.gap_cnt)
	return (a.pct_id > b.pct_id);
    else
	return (a.gap_cnt < b.gap_cnt);
}

int
GappedAln::trace_alignment(string tag_1, string tag_2)
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
    int    i, j, cnt, len, gaps;
    double ident;
    string cigar;
    char   buf[id_len];

    vector<AlignRes> alns;
    bool more_paths = true;
    
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
        cigar = "";
        len   = aln_1.length();
        gaps  = 0;
	ident = 0.0;
        i     = 0;
        while (i < len) {
            if (aln_1[i] != '-' && aln_2[i] != '-') {
                cnt = 0;
                do {
		    if (aln_1[i] == aln_2[2]) ident++;
                    cnt++;
                    i++;
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

            } else {
                cnt = 0;
                do {
                    cnt++;
                    i++;
                } while (i < len && aln_2[i] == '-');
                sprintf(buf, "%dI", cnt);
                gaps++;
            }

            cigar += buf;
        }

        alns.push_back(AlignRes(cigar, gaps, (ident / (double) len)));
        
	// cerr << aln_1 << " [" << cigar << ", gaps: " << gaps << "]\n"
        //      << aln_2 << "\n";

    } while (more_paths);


    sort(alns.begin(), alns.end(), compare_alignres);
    this->_aln = alns[0];
    // cerr << "Final alignment: " << cigar << "; gaps: " << alns[0].second << "\n";

    return 1;
}

AlignRes
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
GappedAln::dump_alignment(string tag_1, string tag_2)
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

    cigar.clear();

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

string
apply_cigar_to_seq(const char *seq, vector<pair<char, uint> > &cigar)
{
    uint   size = cigar.size();
    char   op;
    uint   dist, bp, edited_bp, len, stop;
    string edited_seq;

    //
    // Calculate the overall sequence length.
    //
    uint seqlen = 0;
    for (uint i = 0; i < size; i++)
        seqlen += cigar[i].second;

    len = strlen(seq);
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

// int
// rank_alignments(vector<pair<int, string> > &alns)
// {
//     sort(alns.begin(), alns.end(), compare_alns);
    
//     return 0;
// }

// int
// compare_alns(pair<int, string> a, pair<int, string> b)
// {
//     char   op;
//     uint   dist, sum_a, sum_b;

//     vector<pair<char, uint> > cigar;

//     parse_cigar(a.second, cigar);

//     sum_a = 0;
//     sum_b = 0;

//     for (uint i = 0; i < cigar.size(); i++) {
// 	op   = cigar[i].first;
// 	dist = cigar[i].second;
// 	if (op == 'M')
// 	    sum_a += dist;
//     }

//     cigar.clear();

//     parse_cigar(b.second, cigar);

//     for (uint i = 0; i < cigar.size(); i++) {
// 	op   = cigar[i].first;
// 	dist = cigar[i].second;
// 	if (op == 'M')
// 	    sum_b += dist;
//     }

//     return (sum_a < sum_b);
// }

#endif // __GAPPEDALN_H__
