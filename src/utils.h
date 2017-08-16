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

double factorial(double);
double reduced_factorial(double, double);

double log_factorial(double);
double reduced_log_factorial(double, double);

inline
bool almost_equal(double x, double y) {
    if (!std::isnormal(x) || !std::isnormal(y))
        throw std::domain_error("almost_equal");
    return std::abs(x-y) <= 1e-12 * std::abs(std::min(x,y));
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
void check_open (const std::ifstream& fs, const string& path)
    {if (!fs.is_open()) {cerr << "Error: Failed to open '" << path << "' for reading.\n"; throw exception();}}
inline
void check_open (const std::ofstream& fs, const string& path)
    {if (!fs.is_open()) {cerr << "Error: Failed to open '" << path << "' for writing.\n"; throw exception();}}
inline
void check_open (const gzFile fs, const string& path)
    {if (fs == NULL) {cerr << "Error: Failed to gz-open file '" << path << "'.\n"; throw exception();}}

//
// Class to read lines from a plain text or compressed file indifferently.
//
class VersatileLineReader {
    const string path_;
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
    ~VersatileLineReader();

    // Reads one line from the file, removing the trailing '\n' (and '\r', if any).
    // Returns false on EOF, or throws an exception if the file doesn't end with a newline.
    // e.g.:
    // const char* line; size_t len; while (file.getline(line, len)) {...}
    bool getline(const char*& line, size_t& len);

    const string& path() const {return path_;}
    size_t line_number() const {return line_number_;} // 1-based.
};

class VersatileWriter {
    const string path_;
    bool is_gzipped_;
    ofstream ofs_;
    gzFile gzfile_;

public:
    VersatileWriter(const string& path);
    ~VersatileWriter() {if(is_gzipped_) gzclose(gzfile_);}

    const string& path() const {return path_;}

    friend VersatileWriter& operator<< (VersatileWriter& w, char c)
        {if (w.is_gzipped_) gzputc(w.gzfile_, c); else w.ofs_ << c; return w;}
    friend VersatileWriter& operator<< (VersatileWriter& w, const char* s)
        {if (w.is_gzipped_) gzputs(w.gzfile_, s); else w.ofs_ << s; return w;}
    friend VersatileWriter& operator<< (VersatileWriter& w, const string& s)
        {if (w.is_gzipped_) gzwrite(w.gzfile_, s.c_str(), s.length()); else w.ofs_ << s; return w;}
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
VersatileWriter& operator<< (VersatileWriter& w, int i) {
    if (w.is_gzipped_) {
        char buf[16];
        sprintf(buf, "%d", i);
        gzputs(w.gzfile_, buf);
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
        gzputs(w.gzfile_, buf);
    } else {
        w.ofs_ << i;
    }
    return w;
}

#endif // __UTILS_H__
