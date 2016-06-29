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
using std::ofstream;
using std::cerr;

using std::exception;

class VcfAbstractParser;
class VcfHeader;

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

const string base_header = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO";

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
 * ==========
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
    vector<string> samples;

    Vcf::RType type;

    //map<string, size_t> allele_indexes;
    //void refresh_allele_indexes();

    VcfRecord()
    : chrom(), pos(-1), id(), alleles(), qual(), filter(), info(), format(),
      samples(), type(Vcf::RType::null)
    {}

    // Clears all the members.
    inline void clear();

    // Returns (first allele, second allele), '-1' meaning no data.
    inline pair<int, int> parse_genotype(const string& sample) const;

    // Returns alleles[i], provided it actually exists. (Otherwise it is easy to write a VCF that causes a segfault.)
    inline const string& allele(size_t index) const;
};

/*
 * VcfMeta
 * ==========
 * Represents one line of VCF metainformation.
 */
class VcfMeta {
    string key_;
    string value_;
public:
    VcfMeta(const string& k, const string& p) : key_(k), value_(p) {}

    const string& key() const {return key_;}
    const string& value() const {return value_;}

    static const map<string, VcfMeta> predefined;
};
//class Info : public Meta {...};

/*
 * VcfHeader
 * ==========
 * Stores the contents of a VCF header.
 */
class VcfHeader {
    vector<string> samples_;
    vector<VcfMeta> meta_;

    // map<string, size_t> sample_indexes;
    // map<string, vector<size_t> > meta_indexes;
public:
    VcfHeader() : samples_(), meta_() {}

    // Adds the meta lines VERSION, FILEDATE and SOURCE
    void init_meta(const string& fileformat = "VCFv4.2");

    void add_sample(const string& s) {samples_.push_back(s);}
    const vector<string>& samples() const {return samples_;}

    void add_meta(const VcfMeta& m) {meta_.push_back(m);}
    const vector<VcfMeta>& meta() const {return meta_;}
};

/*
 * VcfAbstractParser
 * ==========
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
    vector<bool> samples_to_keep_; // If set, the parser will return records that include these samples only.

    size_t line_number_;
    char line_[Vcf::line_buf_size];
    bool eol_; // Set by check_eol(). True if the buffer reaches the end of the currently parsed line.
    bool eof_; // Set by getline(). True if EOF was reached.
    vector<char*> tabs_; // Vector of pointers to {first_tab, second_tab, ..., \0 } in [line_].
    vector<char*> bounds_; // Vector of pointers to {leading_tab, first_sep, second_sep, ..., (trailing \t or \0) } in [line_].
    size_t sample_index_; // Index of the current (next) sample.

    virtual void getline(char* ptr, size_t n) =0;
    virtual void check_eol() =0; // rem. The implementation in gzparser relies on [tabs_] to access the end of the string in [line_].
    inline void read_to_eol(); // Reads while [eol_] is false.

    // Add a sample to [record.samples_] if [samples_to_keep_.at(sample_index_)] is true.
    inline void add_sample(VcfRecord& record, char* tab1, char* tab2);

public:
    VcfAbstractParser(const string& path);
    virtual ~VcfAbstractParser() {}

    // Assess the state of the underlying file.
    virtual bool fail() =0;

    // Getters.
    const string& path() const {return path_;};
    const VcfHeader& header() const {return header_;};
    size_t line_number() const {return line_number_;}

    // Parse the header.
    void read_header();

    // Read a record. Returns false on EOF, true otherwise.
    bool next_record(VcfRecord& record);

    // Set the samples to keep.
    void samples_to_keep(const set<string>& samples);
};

/*
 * VcfParser
 * =========
 * Implements VcfAbstractParser for plain text files.
 */
class VcfParser : public VcfAbstractParser {
    ifstream file_;
    void getline(char* ptr, size_t n) {file_.getline(ptr, n); eof_ = file_.eof();}
    void check_eol() {eol_ = ! file_.fail(); file_.clear();}
public:
    VcfParser(const string& path) : VcfAbstractParser(path), file_(path) {}
    bool fail() {return file_.fail();}
};

#ifdef HAVE_LIBZ
/*
 * VcfGzParser
 * ==========
 * Implements VcfAbstractParser for gzipped files.
 */
class VcfGzParser : public VcfAbstractParser {
    gzFile file_;
    void getline(char* ptr, size_t n) {if (gzgets(file_, ptr, n) == NULL) eof_ = true; else eof_ = false;}
    inline void check_eol();

    VcfGzParser(VcfGzParser& p) = delete; // No copy constructor.
    VcfGzParser& operator=(VcfGzParser& p) = delete;
public:
    VcfGzParser(const string& path);
    ~VcfGzParser() {if(file_) gzclose(file_);}
    bool fail() {return file_ == NULL;}
};
#endif // HAVE_LIBZ

/*
 * VcfWriter
 * ==========
 */
class VcfWriter {
private:
    const string path_;
    ofstream file_;

public:
    VcfWriter(const string& path) : path_(path), file_(path) {}
    bool fail() {return file_.fail();}

    void write_header(const VcfHeader& h);
    void write_record(const VcfRecord& r, const VcfHeader& h);
};

/*
 * Inline methods.
 * ==========
 */

inline void
get_bounds(vector<char*>& bounds, char* tab1, char* tab2, char sep)
{
    bounds.clear();
    bounds.push_back(tab1);
    do {
        bounds.push_back(strchr(bounds.back()+1, sep)); // this is last+1 if tab1==tab2
    } while (bounds.back());
    bounds.pop_back();
    bounds.push_back(tab2);
}

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
}

inline
pair<int, int> VcfRecord::parse_genotype(const string& sample) const {

    pair<int, int> genotype = {-1,-1};

    if (format.empty()
            || format[0] != "GT"
            || sample.empty()
            || sample[0] == '.') {
        return genotype;
    }
    const char* first = sample.c_str();
    const char* slash = strchr(first, '/');
    if (slash == NULL) {
        slash = strchr(first, '|');
        if (slash == NULL) {
            cerr << "Error: Malformed VCF genotype field '" << sample
                    << "', at marker '" << chrom << ":" << pos
                    << "'.\n";
            throw exception();
        }
    }
    if (*(slash+1) == '.') {
        static bool printed = false;
        if (not printed) {
            // Print the warning once.
            cerr << "Notice: Treating incomplete genotypes (e.g. '1/.') as missing.\n";
            printed = true;
        }
        return genotype;
    }

    const char* colon = strchr(slash, ':');
    try {
        genotype.first = std::stoi(string(first, slash));
        genotype.second = std::stoi(colon==NULL ? string(slash+1) : string(slash+1, colon));
    } catch (std::invalid_argument& e) {
        cerr << "Error: Malformed VCF genotype '" << sample
             << "', at marker '" << chrom << ":" << pos
             << "'.\n";
        throw e;
    }

    return genotype;
}

inline
const string& VcfRecord::allele(size_t index) const {
    try {
        return alleles.at(index);
    } catch (std::out_of_range& e) {
        cerr << "Error: Malformed VCF genotype for marker '"
             << chrom << ":" << pos
             << "' (allele index does not exist).\n";
        throw e;
    }
}


inline
void
VcfAbstractParser::read_to_eol()
{
    while(!eol_) {
        getline(line_, Vcf::line_buf_size);
        if(eof_) {
            cerr << "Error: VCF file " << path_ << " does not end with a newline.\n";
            throw exception();
        }
        tabs_.clear();
        tabs_.push_back(line_+strlen(line_));
        check_eol();
    }
}

inline void
VcfAbstractParser::add_sample(VcfRecord& record, char* tab1, char* tab2)
{
    if(tab2 == tab1+1)
        cerr << "Warning: malformed VCF record line (empty SAMPLE field should be marked by a dot)."
             << " Line " << line_number_ << " in file " << path_ << "'.\n";

    if (samples_to_keep_.empty() || samples_to_keep_.at(sample_index_))
        record.samples.push_back(string(tab1+1, tab2));

    ++sample_index_;
}

#ifdef HAVE_LIBZ
inline
void VcfGzParser::check_eol() {
    if(*(tabs_.back()-1) == '\n') { // rem. safe, gzgets never returns a null string
        eol_ = true;
        tabs_.back() = tabs_.back()-1;
        *tabs_.back() = '\0';
    } else {
        eol_= false;
    }
}
#endif // HAVE_LIBZ

#endif // __VCF_H__
