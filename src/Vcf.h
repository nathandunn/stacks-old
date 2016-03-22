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

#ifndef __VCF_H__
#define __VCF_H__

#include "config.h"

#include <fstream>
#include <exception>
#include <utility>
#include <string>
#include <vector>
#include <set>
#include <map>

#ifdef HAVE_LIBZ
#include <zlib.h>
#endif

using std::size_t;
using std::pair;
using std::vector;
using std::string;
using std::set;
using std::map;
using std::ifstream;

#include "constants.h"

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
        expl, // record with explicitly written alleles, e.g. SNPs, small indels or entirely defined haplotypes.
        symbolic, // e.g. ALT is "<DUP:TANDEM>".
        breakend // e.g. ALT is "G]17:198982]".
    };

    // Constants for the parser.
    const size_t line_buf_size = 4096;
    const set<string> common_format_fields ({"GT","AD","DP","GL"});
}

class VcfAbstractParser;

/*
 * VcfRecord
 * Datastructure to store VCF records
 */
struct VcfRecord {

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

    VcfRecord()
    : chrom_(), pos_(-1), id_(), ref_(), type_(Vcf::null),  alt_(), qual_(), filter_(), info_(), format_(), samples_()
    {}

private:
    //map<string,size_t> format_subfields_indexes_; //todo

public:
    void clear(); // Clears all the members.

    // Methods for further parsing
    //pair<size_t, size_t> parse_genotype(size_t sample_index) const; // Returns (first allele, second allele)
    //const string& allele(size_t index) const {if (index == 0) return ref_; else return alt_.at(index-1);}
};

/*
 * VcfHeader
 * Stores the contents of a VCF header.
 */
class VcfHeader {
//public:
    //struct Vcf_info {}; // Not implemented yet.
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
    VcfHeader() : version_(), date_(), source_(), reference_(), samples_() {}

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
 * VcfAbstractParser
 *
 * Main class for parsing VCF. The derived non-abstract classes
 * VcfParser and VcfGzParser only add the file_ attribute and
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
class VcfAbstractParser {
protected:
    const string path_;
    VcfHeader header_;
    size_t header_lines_; // Number of header lines there were (set by the constructor).

    set<string> format_fields_to_keep_; // If set, the parser will return records with only these fields included.

    size_t line_number_;
    char line_[Vcf::line_buf_size];
    bool eol_; // Set by check_eol(). True if the buffer reaches the end of the currently parsed line.
    bool eof_; // Set by getline(). True if EOF was reached.
    vector<char*> tabs_; // Vector of pointers to {first_tab, second_tab, ..., \0 } in line_
    vector<char*> bounds_; // vector of pointers to {leading_tab, first_sep, second_sep, ..., (trailing \t or \0) } in line_
    vector<bool> kept_format_fields_; // Keeps track of which FORMAT subfields were kept for this record.

    // Parses the header.
    // Called by the open method of the non-abstract derived classes (VcfParser and VcfGzParser).
    // Throws an exception if the header is malformed.
    void read_header();

    virtual void getline(char* ptr, size_t n) =0;
    virtual void check_eol() =0; // n.b. The implementation in gzparser relies on tabs_ to access the end of the string in line_.
    void read_while_not_eol(); // Read while eol_ is false.

    void add_sample(VcfRecord& record, char* tab1, char* tab2); // Basically 'record.samples_.push_back(parsed_value)'.

public:
    VcfAbstractParser();
    virtual ~VcfAbstractParser() {}

    // Opens the VCF file and reads the header.
    // Returns 0 if successful, 1 if the file could not be opened.
    virtual int open(const string& path) =0;

    // Getters.
    const string& path() const {return path_;};
    const VcfHeader& header() const {return header_;};
    size_t header_lines() const {return header_lines_;}
    size_t line_number() const {return line_number_;}

    // Setters.
    void format_fields_to_keep(const set<string>& fields) {format_fields_to_keep_ = fields;}

    // Reads a record. Returns false on EOF, true otherwise.
    bool next_record(VcfRecord& record);
};

/*
 * VcfParser
 * Implements VcfAbstractParser for plain text files.
 */
class VcfParser : public VcfAbstractParser {
    ifstream file_;
    void getline(char* ptr, size_t n) {file_.getline(ptr, n); eof_ = file_.eof();}
    void check_eol() {eol_ = ! file_.fail(); file_.clear();}
public:
    int open(const string& path);
};

#ifdef HAVE_LIBZ
/*
 * VcfGzParser
 * Implements VcfAbstractParser for gzipped files.
 */
class VcfGzParser : public VcfAbstractParser {
    gzFile file_;
    void getline(char* ptr, size_t n) {gzgets(file_, ptr, n); eof_ = gzeof(file_);}
    void check_eol();

    VcfGzParser(VcfGzParser& p) = delete; // No copy constructor.
    VcfGzParser& operator=(VcfGzParser& p) = delete;
public:
    VcfGzParser() : file_(NULL) {}
    ~VcfGzParser() {if(file_) gzclose(file_);}

    int open(const string& path);
};
#endif // HAVE_LIBZ

#endif // __VCF_H__
