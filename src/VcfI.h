// -*-mode:c++; c-style:k&r; c-basic-offset:4;-*-
//
// Copyright 2010-2015, Julian Catchen <jcatchen@uoregon.edu>
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

#ifndef __VCFI_H__
#define __VCFI_H__

#include <fstream>
#include <exception>
#include <utility>
#include <string>
#include <vector>
#include <set>

#include "config.h"

#ifdef HAVE_LIBZ
#include <zlib.h>
#endif

#include "constants.h"

using namespace std;

namespace Vcf {
    const size_t base_fields_no = 8; // CHROM POS ID REF ALT QUAL FILTER INFO [FORMAT SAMPLE ...]]
    const size_t chrom = 0;
    const size_t pos = 1;
    const size_t id = 2;
    const size_t ref = 3;
    const size_t alt = 4;
    const size_t qual = 5;
    const size_t filter = 6;
    const size_t info = 7;
    const size_t format = 8;
    const size_t first_sample = 9;

    enum Record_type {
        null,
        expl, // explicit : marker with explicitly written alleles (e.g. SNPs, small indels or entirely defined haplotypes)
        symbolic,
        breakend
    };

    const size_t line_buf_size = 4096;
    const set<string> common_format_fields ({"GT","AD","DP","GL"});
}

class Vcf_abstractparser;

/*
 * Vcf_basicrecord
 *
 * Fast datastructure to store VCF records
 */
class Vcf_basicrecord {
public:
    string chrom_; // required
    size_t pos_; // required
    string id_;
    string ref_; //required, case insensitive
    Vcf::Record_type type_;
    vector<string> alt_; // case insensitive
    string qual_;
    vector<string> filter_;
    vector<pair<string, string> > info_;
    vector<string> format_;
    vector<vector<string> > samples_;

    Vcf_basicrecord() : chrom_(), pos_(-1), id_(), ref_(), type_(Vcf::null),  alt_(), qual_(), filter_(), info_(), format_(), samples_() {}

    void clear();
};

/*
 * Vcf_header
 *
 * Stores the contents of a VCF header.
 */
class Vcf_header {
//public:
    //struct Vcf_info {};
    //struct Vcf_format {};

protected:
    string version_;
    string date_;
    string source_;
    string reference_;
    //map<string, string> contigs_; // (id,desc) items
    //map<string, string> alt_keys_; // same
    //map<string, string> filter_tags_; // same
    //map<string, map<string, string> > info_keys_;
    //map<string, map<string, string> > format_keys_;
    vector<string> samples_;
public:
    Vcf_header() : version_(), date_(), source_(), reference_(), samples_() {}

    //getters
    const string& version() const {return version_;}
    const string& date() const {return date_;}
    const string& source() const {return source_;}
    const string& reference() const {return reference_;}
    const vector<string>& samples() const {return samples_;}

    //setters
    void version(const string& v) {version_ = v;}
    void date(const string& d) {date_ = d;}
    void source(const string& s) {source_ = s;}
    void reference(const string& r) {reference_ = r;}

    void add_sample(const string& s) {samples_.push_back(s);}
};

/*
 * Vcf_abstractparser
 *
 * Main class for parsing VCF. The derived non-abstract classes
 * Vcf_parser and Vcf_gzparser only add the file_ attribute and
 * implement the getline(), eof() and check_eol() methods.
 *
 * At present, the parser does not handle :
 * -- records for which the fixed fields (CHROM to FORMAT) do not
 *    fit in the buffer. Such records usually correspond to large
 *    indels with several alleles. In these cases, the record
 *    passed to next_record() will be empty and have the default
 *    type 'null'.
 * -- records describing symbolic or breakends alleles. In these
 *    cases the record passed to next_record() will be of type
 *    'symbolic' or 'breakend', with only the fields CHROM, POS,
 *    ID and REF filled in.
 */
class Vcf_abstractparser {
protected:
    const string path_;
    Vcf_header header_;
    size_t header_lines_; // Number of header lines there were (set by the constructor).

    set<string> format_fields_to_keep_; // If set, the parser will return records with only these fields included.

    size_t line_number_;
    char line_[Vcf::line_buf_size];
    bool eol_; // Set by check_eol(). True if the buffer reaches the end of the currently parsed line.
    vector<char*> tabs_; // Vector of pointers to {first_tab, second_tab, ..., \0 } in line_
    vector<char*> bounds_; // vector of pointers to {leading_tab, first_sep, second_sep, ..., (trailing \t or \0) } in line_
    vector<bool> kept_format_fields_; // Keeps track of which FORMAT subfields were kept for this record.

    void read_header(); // Called by the constructors of the non-abstract derived classes (Vcf_parser and Vcf_gzparser).

    virtual void getline(char* ptr, size_t n) =0;
    virtual bool eof() const =0;
    virtual void check_eol() =0; // n.b. The implementation in gzparser relies on tabs_ to access the end of the string in line_.
    void read_while_not_eol(); // Read while eol_ is false.

    void add_sample(Vcf_basicrecord& record, char* tab1, char* tab2); // Basically 'record.samples_.push_back(parsed_value)'.

public:
    Vcf_abstractparser(const string& path); // Initialize the VCF parser for the given file, i.e. open the file and read the header.
    virtual ~Vcf_abstractparser() {}

    // Getters.
    const string& path() const {return path_;};
    const Vcf_header& header() const {return header_;};
    size_t header_lines() const {return header_lines_;}
    size_t line_number() const {return line_number_;}

    // Setters.
    void format_fields_to_keep(const set<string>& fields) {format_fields_to_keep_ = fields;}

    // Read a record. Returns 0 on EOF, 1 otherwise.
    int next_record(Vcf_basicrecord& record);
};

/*
 * Vcf_parser
 *
 * Implements Vcf_abstractparser for plain text files.
 */
class Vcf_parser : public Vcf_abstractparser {
    ifstream file_;
    void getline(char* ptr, size_t n) {file_.getline(ptr, n);}
    bool eof() const {return file_.eof();}
    void check_eol() {eol_ = ! file_.fail(); file_.clear();}
public:
    Vcf_parser(const string& path);
};

/*
 * Vcf_gzparser
 *
 * Implements Vcf_abstractparser for gzipped files.
 */
#ifdef HAVE_LIBZ
class Vcf_gzparser : public Vcf_abstractparser {
    gzFile file_;
    void getline(char* ptr, size_t n) {gzgets(file_, ptr, n);}
    bool eof() const {return gzeof(file_);}
    void check_eol();

    Vcf_gzparser(Vcf_gzparser& p) = delete;
    Vcf_gzparser& operator=(Vcf_gzparser& p) = delete;
public:
    Vcf_gzparser(const string& path);
    ~Vcf_gzparser() {gzclose(file_);}
};
#else
class Vcf_gzparser : public Vcf_abstractparser {
public:
    Vcf_gzparser(const string& path) : Vcf_abstractparser(path) {cerr << "Trying to read a gzipped file, but Gzip support was not enabled when Stacks was compiled." << endl; throw exception();}

    // Dummy definitions for the pure virtual methods.
    void getline(char* ptr, size_t n) {}
    bool eof() const {return true;}
    void check_eol() {}
};
#endif // HAVE_LIBZ

#endif // __VCFI_H__
