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
#include <iostream>
#include <exception>
#include <stdexcept>
#include <utility>
#include <cstring>
#include <cstdlib>
#include <string>
#include <vector>
#include <set>
#include <map>

#ifdef HAVE_LIBZ
#include <zlib.h>
#endif

#include "constants.h"

using std::size_t;
using std::pair;
using std::vector;
using std::string;
using std::set;
using std::map;
using std::ifstream;
using std::cerr;
using std::exception;
using std::out_of_range;

class VcfAbstractParser;

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

const size_t gt_subfield = 0;

// enum for record types
enum class RType {
    null,
    expl,      // Record with explicitly written alleles, e.g. SNPs, small indels or entirely defined haplotypes
    invariant, // ALT is empty
    symbolic,  // ALT is e.g. "<DUP:TANDEM>"
    breakend   // ALT is e.g. "G]17:198982]"
};

// Constants for the parser.
const size_t line_buf_size = 4096;

// Open the given VCF file using VcfParser or VcfGzParser, depending on the
// suffix of the file. Return NULL if the opening failed.
// (The pointee is dynamically allocated and should be deleted.)
VcfAbstractParser* adaptive_open(const string& path);

} //namespace Vcf

/*
 * VcfRecord
 * Datastructure to store VCF records
 */
struct VcfRecord {

    string chrom; // required
    size_t pos; // required
    string id;
    vector<string> alleles; // allele 0 is REF and is required ; case insensitive
    string qual;
    vector<string> filter;
    vector<pair<string, string> > info;
    vector<string> format;
    vector<vector<string> > samples;

    Vcf::RType type;
    bool no_gt; // format[Vcf::gt_subfield_index] != "GT". Trusted by parse_genotype().

    VcfRecord()
    : chrom(), pos(-1), id(), alleles(), qual(), filter(), info(), format(),
      samples(), type(Vcf::RType::null), no_gt(true)
    {}

    inline void clear(); // Clears all the members.

    // Returns (first allele, second allele), or (-1,-1) if the record doesn't have a genotype.
    inline pair<int, int> parse_genotype(const vector<string>& sample) const;
    pair<int, int> parse_genotype(size_t sample_index) const {return parse_genotype(samples[sample_index]);}

    // Returns alleles[i], provided it actually exists. (Otherwise it is easy to write a VCF that causes a segfault.)
    inline const string& allele(size_t index) const;
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
    set<string> format_fields_to_keep_; // If set, the parser will return records that include these fields only.
    vector<bool> samples_to_keep_; // If set, the parser will return records that include these samples only.

    size_t line_number_;
    char line_[Vcf::line_buf_size];
    bool eol_; // Set by check_eol(). True if the buffer reaches the end of the currently parsed line.
    bool eof_; // Set by getline(). True if EOF was reached.
    vector<char*> tabs_; // Vector of pointers to {first_tab, second_tab, ..., \0 } in [line_].
    vector<char*> bounds_; // Vector of pointers to {leading_tab, first_sep, second_sep, ..., (trailing \t or \0) } in [line_].
    vector<bool> kept_format_fields_; // Keeps track of which FORMAT subfields were kept for this record.
    size_t sample_index_; // Index of the current (next) sample.

    // Parses the header.
    // Called by the open method of the non-abstract derived classes (VcfParser and VcfGzParser).
    // Throws an exception if the header is malformed.
    void read_header();

    virtual void getline(char* ptr, size_t n) =0;
    virtual void check_eol() =0; // rem. The implementation in gzparser relies on [tabs_] to access the end of the string in [line_].
    inline void read_to_eol(); // Reads while [eol_] is false.

    // Adds a sample to [record.samples_] if [samples_to_keep_.at(sample_index_)] is true.
    inline void add_sample(VcfRecord& record, char* tab1, char* tab2);

public:
    VcfAbstractParser();
    virtual ~VcfAbstractParser() {}

    // Opens the VCF file and reads the header.
    // Returns true if successful, false if the file could not be opened.
    virtual bool open(const string& path) =0;

    // Sets the format fields to keep.
    void format_fields_to_keep(const set<string>& fields) {format_fields_to_keep_ = fields;}
    // Sets the samples to keep.
    void samples_to_keep(const set<string>& samples);

    // Getters.
    const string& path() const {return path_;};
    const VcfHeader& header() const {return header_;};
    size_t header_lines() const {return header_lines_;}
    size_t line_number() const {return line_number_;}

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
    VcfParser() : file_() {}
    bool open(const string& path);
};

#ifdef HAVE_LIBZ
/*
 * VcfGzParser
 * Implements VcfAbstractParser for gzipped files.
 */
class VcfGzParser : public VcfAbstractParser {
    gzFile file_;
    void getline(char* ptr, size_t n) {if (gzgets(file_, ptr, n) == NULL) eof_ = true; else eof_ = false;}
    inline void check_eol();

    VcfGzParser(VcfGzParser& p) = delete; // No copy constructor.
    VcfGzParser& operator=(VcfGzParser& p) = delete;
public:
    VcfGzParser() : file_(NULL) {}
    ~VcfGzParser() {if(file_) gzclose(file_);}

    bool open(const string& path);
};
#endif // HAVE_LIBZ

inline
void VcfRecord::clear() {
    pos = -1;
    chrom.clear();
    id.clear();
    alleles.clear();
    qual.clear();
    filter.clear();
    info.clear();
    format.clear();
    samples.clear();
    type = Vcf::RType::null;
    no_gt = true;
}

inline
pair<int, int> VcfRecord::parse_genotype(const vector<string>& sample) const {
    pair<int, int> genotype;

    if (no_gt
            || sample[Vcf::gt_subfield].empty()
            || sample[Vcf::gt_subfield].front() == '.') {
        genotype = {-1,-1};
        return genotype;
    }

    const char* first = sample[0].c_str();
    const char* sep = strchr(first, '/');
    if (sep == NULL)
        sep = strchr(first, '|');

    if (sep == NULL // no separator
            || sep == first // first field is empty
            || *(sep+1) == '\0') { // second field is empty
        cerr << "Error: Malformed VCF genotype field '" << first
             << "', at marker '" << chrom << ":" << pos
             << "'.\n";
        throw exception();
    }
    genotype = {atoi(first), atoi(sep+1)};
    return genotype;
}

inline
const string& VcfRecord::allele(size_t index) const {
    try {
        return alleles.at(index);
    } catch (out_of_range& e) {
        cerr << "Error: Malformed VCF genotype for marker '"
             << chrom << ":" << pos
             << "' (allele index does not exist).\n";
        throw e;
    }
}

#endif // __VCF_H__
