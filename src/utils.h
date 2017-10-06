// -*-mode:c++; c-style:k&r; c-basic-offset:4;-*-
//
// Copyright 2010-2015, Julian Catchen <jcatchen@illinois.edu>
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

#ifndef __UTILS_H__
#define __UTILS_H__

#include <cstdlib>
#include <cerrno>
#include <climits>
#include <cmath>
#include <iostream>
#include <utility>
#include <string>

#include <unistd.h>
#include <dirent.h>

#include <zlib.h>

#include "constants.h"
#include "stacks.h"

char   reverse(char);
char  *rev_comp(const char *);
string rev_comp(const string&);
void   reverse_string(char *);
int    is_integer(const char *);
double is_double(const char *);
int    tokenize_string(const char *, vector<string> &);

double factorial(double);
double reduced_factorial(double, double);

double log_factorial(double);
double reduced_log_factorial(double, double);

inline
bool almost_equal(double x, double y) {
    const double precision = 1.0e-9;
    if (x == 0.0 && y == 0.0)
        return true;
    if (!std::isnormal(x) || !std::isnormal(y)) {
        stringstream ss;
        ss << "almost_equal: x=" << x << ", y=" << y;
        throw std::domain_error(ss.str());
    }
    return std::abs(x-y) <= precision * std::abs(std::min(x,y));
}

//
// Comparison functions for the STL sort routine
//
bool compare_ints(int, int);
bool compare_pair(pair<char, int>, pair<char, int>);
bool compare_pair_intint(pair<int, int>, pair<int, int>);
bool compare_pair_intdouble(pair<int, double>, pair<int, double>);
bool compare_pair_stringint(pair<string, int>, pair<string, int>);
bool compare_pair_snp(pair<string, SNP *>, pair<string, SNP *>);
bool compare_pair_haplotype(pair<string, double>, pair<string, double>);
bool compare_pair_haplotype_rev(pair<string, double>, pair<string, double>);
bool compare_str_len(string, string);
struct LessCStrs
{
   bool operator() (const char* str1, const char* str2) const {return strcmp(str1, str2) < 0 ;}
} ;

//
// Comparison classes for STL sets
//
struct int_increasing {
    bool operator() (const int& lhs, const int& rhs) const {
        return lhs < rhs;
    }
};

struct int_decreasing {
    bool operator() (const int& lhs, const int& rhs) const {
        return lhs > rhs;
    }
};

//
// Join a range of elements into a stream.
//
template<typename IterableT, typename SepT>
void join(IterableT elements, const SepT& sep, ostream& os) {
    auto first = elements.begin();
    if (first != elements.end()) {
        os << *first;
        ++first;
        while (first != elements.end()) {
            os << sep << *first;
            ++first;
        }
    }
}

//
// Routines to check that files are open.
//
inline
void check_open (std::ifstream& fs, const string& path) {
    if (!fs.is_open()) {
        cerr << "Error: Failed to open '" << path << "' for reading.\n";
        throw exception();
    }
    fs.exceptions(fs.exceptions() | ios::badbit);
}
inline
void check_open (std::ofstream& fs, const string& path) {
    if (!fs.is_open()) {
        cerr << "Error: Failed to open '" << path << "' for writing.\n";
        throw exception();
    }
    fs.exceptions(fs.exceptions() | ios::badbit);
}
inline
void check_open (const gzFile fs, const string& path)
    {if (fs == NULL) {cerr << "Error: Failed to gz-open file '" << path << "'.\n"; throw exception();}}

//
// Check that a directory exists or try to create it.
//
void check_or_mk_dir(const string& path);

//
// Timing routine "gettime".
//
inline
double gettm() {
#if _POSIX_MONOTONIC_CLOCU
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return ts.tv_sec + ts.tv_nsec / 1.0e9;
#else
    return 0.0;
#endif
}

//
// Class to read lines from a plain text or compressed file indifferently.
//
class VersatileLineReader {
    string path_;
    size_t line_number_;
    bool is_gzipped_;

    ifstream ifs_;
    string ifsbuffer_;

    gzFile gzfile_;
    char* gzbuffer_;
    size_t gzbuffer_size_;
    size_t gzline_len_;
    static const size_t gzbuffer_init_size = 65536;

public:
    VersatileLineReader(const string& path);
    VersatileLineReader();
    ~VersatileLineReader();

    int open(string &);

    //
    // Reads one line from the file, removing the trailing '\n' (and '\r', if any).
    // Returns false on EOF, or throws an exception if the file doesn't end with a newline.
    // e.g.:
    // const char* line; size_t len; while (file.getline(line, len)) {...}
    //
    bool getline(const char*& line, size_t& len);

    const string& path() const {return path_;}
    size_t line_number() const {return line_number_;} // 1-based.
};

class VersatileWriter {
    const string path_;
    bool is_gzipped_;
    ofstream ofs_;
    gzFile gzfile_;

    void gzputs_(const char* s);

public:
    VersatileWriter(const string& path);
    ~VersatileWriter() {if(is_gzipped_) gzclose(gzfile_);}

    const string& path() const {return path_;}

    friend VersatileWriter& operator<< (VersatileWriter& w, char c);
    friend VersatileWriter& operator<< (VersatileWriter& w, const char* s);
    friend VersatileWriter& operator<< (VersatileWriter& w, const string& s);
    friend VersatileWriter& operator<< (VersatileWriter& w, int i);
    friend VersatileWriter& operator<< (VersatileWriter& w, size_t i);
};

//
// Wrapper for directory parsing functions.
// e.g. for(DirIterator e (path); e; ++e) {...}
//
class DirIterator {
    DIR* dir;
    struct dirent* entry;
public:
    DirIterator(const string& dir_path) : dir(NULL), entry(NULL) {
        dir = opendir(dir_path.c_str());
        if (dir == NULL) {
            cerr << "Error: Unable to open directory '" << dir_path << "' for reading.\n";
            throw exception();
        }
        entry = readdir(dir);
    }
    ~DirIterator() {closedir(dir);}

    const char* name() const {return entry->d_name;}

    operator bool() const {return entry!=NULL;}
    DirIterator& operator++() {entry = readdir(dir); return *this;}
    dirent* operator*() {return entry;}
};

// strip_read_number
// ----------
// Given a read name, removes the trailing /1, /2, _1 or _1.
inline
void strip_read_number(string& read_name) {
    if (read_name.size() >= 2) {
        char last = read_name[read_name.size()-1];
        char ante = read_name[read_name.size()-2];
        if ((last == '1' || last == '2') && (ante == '/' || ante == '_')) {
            // Remove the suffix & return.
            read_name.resize(read_name.size()-2);
            return;
        }
    }

    // Unexpected suffix.
    cerr << "Error: Unrecognized read name format: expected '"
         << read_name << "' to end with /1, /2, _1 or _2.\n";
    throw exception();
}

inline
VersatileWriter& operator<< (VersatileWriter& w, char c) {
    if (w.is_gzipped_) {
        if (gzputc(w.gzfile_, c) == -1)
            throw ios::failure("gzputc");
    } else {
        w.ofs_ << c;
    }
    return w;
}

inline
VersatileWriter& operator<< (VersatileWriter& w, const char* s) {
    if (w.is_gzipped_)
        w.gzputs_(s);
    else
        w.ofs_ << s;
    return w;
}

inline
void VersatileWriter::gzputs_(const char* s) {
    if (gzputs(gzfile_, s) == -1)
        throw ios::failure("gzputs");
}


inline
VersatileWriter& operator<< (VersatileWriter& w, const string& s) {
    if (w.is_gzipped_) {
        if (!s.empty() && gzwrite(w.gzfile_, s.c_str(), s.length()) <= 0)
            throw ios::failure("gzwrite");
    } else {
        w.ofs_ << s;
    }
    return w;
}

inline
VersatileWriter& operator<< (VersatileWriter& w, int i) {
    if (w.is_gzipped_) {
        char buf[16];
        sprintf(buf, "%d", i);
        w.gzputs_(buf);
    } else {
        w.ofs_ << i;
    }
    return w;
}

inline
VersatileWriter& operator<< (VersatileWriter& w, size_t i) {
    if (w.is_gzipped_) {
        char buf[32];
        sprintf(buf, "%zu", i);
        w.gzputs_(buf);
    } else {
        w.ofs_ << i;
    }
    return w;
}

#endif // __UTILS_H__
