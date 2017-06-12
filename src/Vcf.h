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
#include "nucleotides.h"
#include "utils.h"
#include "stacks.h"

using std::out_of_range;

class VcfAbstractParser;
class VcfHeader;

namespace Vcf {

const size_t base_fields_no = 8; // CHROM POS ID REF ALT QUAL FILTER INFO [FORMAT SAMPLE ...]]
const string base_fields = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO";
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
// suffix of the file. Throws on failure.
// (The pointee is dynamically allocated and should be deleted.)
unique_ptr<VcfAbstractParser> adaptive_open(const string& path);

} //namespace Vcf

/*
 * VcfRecord
 * ==========
 * Datastructure to store VCF records
 */
class VcfRecord {
    vector<char> buffer_;
    vector<size_t> strings_; // Positions of all independent c-strings in the buffer.

    size_t pos_i_() const {return 0;} // Index in strings_ of the position c-string.
    size_t id_i_() const {return 1;}
    size_t allele0_i_() const {return 2;}
    size_t filters_i_;
    size_t qual_i_() const {return filters_i_+1;}
    size_t info0_i_() const {return filters_i_+2;}
    size_t format0_i_;
    size_t sample0_i_;

    const char* i2str(size_t i) const {return (char*) &buffer_[strings_.at(i)];}

    void append_str(const string& s) {append_str(s.c_str(), s.length());}
    void append_str(const char* str, size_t len);

public:
    VcfRecord()
    : buffer_(), strings_(), filters_i_(SIZE_MAX), format0_i_(SIZE_MAX), sample0_i_(SIZE_MAX)
    {}

    void assign(const char* rec, size_t len, const VcfHeader& header);
    bool check_record() const; // Not implemented. Check that empty fields are '.', etc.

    const char* chrom()          const {return buffer_.data();}
         size_t pos()            const {return atoi(i2str(pos_i_()));}
    const char* id()             const {return i2str(id_i_());}
    const char* allele(size_t i) const {return i2str(allele0_i_() + i);}
    const char* filters()        const {return i2str(filters_i_);}
    const char* qual()           const {return i2str(qual_i_());}
    const char* info(size_t i)   const {return i2str(info0_i_() + i);}
    const char* format(size_t i) const {return i2str(format0_i_ + i);} // (Will throw if the field doesn't exit.)
    const char* sample(size_t i) const {return i2str(sample0_i_ + i);}

    size_t n_alleles() const {return filters_i_ - allele0_i_();}
    size_t n_infos()   const {return format0_i_ - info0_i_();}
    size_t n_formats() const {return sample0_i_ - format0_i_;}
    size_t n_samples() const {return strings_.size() - sample0_i_;}

    Vcf::RType type() const;

    // Record creation functions.
    void clear();
    void append_chrom(const string& s);
    void append_pos(size_t pos);
    void append_id(const string& s);
    void append_allele(Nt2 nt);
    void append_allele(const string& s);
    void append_filters(const string& s);
    void append_qual(const string& s);
    void append_info(const string& s);
    void append_format(const string& s);
    void append_sample(const string& s);

    // TODO Deprecated functions.
    Vcf::RType& type_m();
    string& chrom_m();
    size_t& pos_m();
    string& id_m();
    vector<string>& alleles_m();
    string& qual_m();
    vector<string>& filter_m();
    vector<string>& info_m();
    vector<string>& format_m();
    vector<string>& samples_m();

    inline size_t index_of_gt_subfield(const string& key) const; // SIZE_MAX if not found.
    inline string parse_gt_subfield(const char* sample, size_t index) const;

    // Returns (first allele, second allele), '-1' meaning no data.
    // (The second version is for use with internal, stacks-generated files.)
    inline pair<int, int> parse_genotype(const char* sample) const;
    inline pair<int, int> parse_genotype_nochecks(const char* sample) const;

    inline bool is_snp() const;

    struct util {
        static string fmt_info_af(const vector<double>& alt_freqs);
        static string fmt_gt_gl(const vector<Nt2>& alleles, const GtLiks& liks);
        static GtLiks parse_gt_gl(const vector<Nt2>& alleles, const string& gl);
        static size_t n_genotypes(size_t n_alleles) {return (n_alleles*(n_alleles+1))/2;}

        // Builds the haplotypes of a sample over a set of (phased) records.
        // (At most one phase set is expected.)
        // Returns false if the haplotypes were incomplete (N's).
        static bool build_haps(pair<string,string>& haplotypes,
                                  const vector<const VcfRecord*>& snp_records,
                                  size_t sample_index);
    };
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

    struct predefs {
        static const VcfMeta info_AD;
        static const VcfMeta info_AF;
        static const VcfMeta info_DP;
        static const VcfMeta info_NS;

        static const VcfMeta format_AD;
        static const VcfMeta format_DP;
        static const VcfMeta format_GL;
        static const VcfMeta format_GQ;
        static const VcfMeta format_GT;
        static const VcfMeta format_HQ;

        // Custom.
        static const VcfMeta info_locori;
    };
};

/*
 * VcfHeader
 * ==========
 * Stores the contents of a VCF header.
 */
class VcfHeader {
    vector<string> samples_;
    vector<VcfMeta> meta_;

    map<string, size_t> sample_indexes_;
    // map<string, vector<size_t> > meta_indexes;
public:
    VcfHeader(const string& version = "VCFv4.2")
        : samples_(), meta_()
        {init_meta(version);}

    const vector<VcfMeta>& meta() const {return meta_;}
    const vector<string>& samples() const {return samples_;}
    const map<string, size_t>& sample_indexes() const {return sample_indexes_;}

    // Adds the meta lines VERSION, FILEDATE and SOURCE
    void add_meta(const VcfMeta& m) {meta_.push_back(m);}
    void add_sample(const string& s) {samples_.push_back(s); sample_indexes_.insert({s, samples_.size()-1});}

private:
    void init_meta(const string& version);
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
};

/*
 * VcfParser
 * =========
 * Implements VcfAbstractParser for plain text files.
 */
class VcfParser : public VcfAbstractParser {
    ifstream file_;
    inline void getline(char* ptr, size_t n);
    void check_eol() {eol_ = ! file_.fail(); file_.clear();}
public:
    VcfParser(const string& path) : VcfAbstractParser(path), file_(path) {check_open(file_, path);}
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
    inline void getline(char* ptr, size_t n);
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
    const VcfHeader header_;

public:
    VcfWriter(const string& path, VcfHeader&& header)
        : path_(path), file_(path), header_(header)
        {check_open(file_, path_); write_header();}

    const VcfHeader& header() const {return header_;}
    void write_record(const VcfRecord& r);

private:
    void write_header();
};

/*
 * Inline methods.
 * ==========
 */

inline
void VcfRecord::append_str(const char* str, size_t len) {
    strings_.push_back(buffer_.size());
    buffer_.resize(buffer_.size()+len+1);
    memcpy(buffer_.data(), str, len+1);
}

inline
void VcfRecord::append_chrom(const string& s) {
    assert(buffer_.empty());
    append_str(s);
}

inline
void VcfRecord::append_pos(size_t pos) {
    assert(strings_.empty());
    append_str(to_string(pos));
}

inline
void VcfRecord::append_id(const string& s) {
    assert(strings_.size()==1);
    append_str(s);
}

inline
void VcfRecord::append_allele(Nt2 nt) {
    assert(strings_.size() > 1);
    assert(filters_i_ == SIZE_MAX);
    char a[2] {char(nt), '\0'};
    append_str(a, 1);
}

inline
void VcfRecord::append_allele(const string& s) {
    assert(strings_.size() > 1);
    assert(filters_i_ == SIZE_MAX);
    append_str(s);
}

inline
void VcfRecord::append_filters(const string& s) {
    assert(strings_.size() > 2);
    assert(filters_i_ == SIZE_MAX);
    filters_i_ = strings_.size();
    append_str(s);
}

inline
void VcfRecord::append_qual(const string& s) {
    assert(filters_i_ != SIZE_MAX);
    assert(strings_.size() == filters_i_ + 1);
    append_str(s);
}

inline
void VcfRecord::append_info(const string& s) {
    assert(strings_.size() > filters_i_ + 1);
    assert(format0_i_ == SIZE_MAX);
    append_str(s);
}

inline
void VcfRecord::append_format(const string& s) {
    assert(strings_.size() > filters_i_ + 2);
    assert(sample0_i_ == SIZE_MAX);
    if (format0_i_ == SIZE_MAX)
        format0_i_ = strings_.size();
    append_str(s);
}

inline
void VcfRecord::append_sample(const string& s) {
    assert(format0_i_ != SIZE_MAX);
    if (sample0_i_ == SIZE_MAX)
        sample0_i_ = strings_.size();
    append_str(s);
}

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
    buffer_.clear();
    strings_.clear();
    filters_i_ = SIZE_MAX;
    format0_i_ = SIZE_MAX;
    sample0_i_ = SIZE_MAX;
}

inline
size_t VcfRecord::index_of_gt_subfield(const string& key) const {
    for (size_t i=0; i<n_formats(); ++i)
        if (strcmp(key.c_str(), format(i)) == 0)
            return i;
    return SIZE_MAX;
}

inline
string VcfRecord::parse_gt_subfield(const char* sample, size_t index) const {
    string subf;

    // Skip the first [index] colons.
    const char* first = sample;
    for(size_t i=0; i<index; ++i) {
        first = strchr(first, ':');
        if (first == NULL)
            // The requested field is not explicitly written, return the empty string.
            return subf;
        else
            first += 1;
    }

    const char* last = strchr(first, ':');
    subf = last == NULL ? string(first) : string(first, last);
    return subf;
}

inline
pair<int, int> VcfRecord::parse_genotype(const char* sample) const {

    pair<int, int> genotype = {-1,-1};

    assert(sample != NULL && sample[0] != '\0');
    if (n_formats() == 0
            || strcmp(format(0), "GT") == 0
            || sample[0] == '.') {
        return genotype;
    }
    const char* first = sample;
    const char* slash = strchr(first, '/');
    if (slash == NULL) {
        slash = strchr(first, '|');
        if (slash == NULL) {
            cerr << "Error: Malformed VCF genotype field '" << sample
                    << "', at marker '" << chrom() << ":" << pos()
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
        genotype.first = stoi(string(first, slash));
        genotype.second = stoi(colon==NULL ? string(slash+1) : string(slash+1, colon));
        if (genotype.first < 0
            || genotype.first >= int(n_alleles())
            || genotype.second < 0
            || genotype.second >= int(n_alleles()))
            throw exception();
    } catch (exception& e) {
        cerr << "Error: Malformed VCF genotype '" << sample
             << "', at marker '" << chrom() << ":" << pos()
             << "'.\n";
        throw e;
    }

    return genotype;
}

inline pair<int, int> VcfRecord::parse_genotype_nochecks(const char* sample) const {
    assert(n_formats() > 0 && strcmp(format(0),"GT")==0);
    assert(sample != NULL && sample[0] != '\0');

    pair<int, int> genotype = {-1,-1};
    if (sample[0] == '.')
        return genotype;

    const char* start = sample;
    char* end;

    // First allele.
    genotype.first = strtol(start, &end, 10);
    assert(end != start);
    assert(*end == '/' || *end == '|');

    // Second allele.
    start = end;
    ++start;
    genotype.second = strtol(start, &end, 10);
    assert(end != start);
    assert(*end == ':' || *end == '\0');

    return genotype;
}

inline
bool VcfRecord::is_snp() const {
    for (size_t i=0; i<n_alleles(); ++i)
        if (strlen(allele(i)) > 1)
            return false;
    return true;
}

inline
void VcfAbstractParser::read_to_eol()
{
    while(!eol_) {
        getline(line_, Vcf::line_buf_size);
        if(eof_) {
            eol_=true;
        } else {
            tabs_.clear();
            tabs_.push_back(line_+strlen(line_));
            check_eol();
        }
    }
}

inline
void VcfAbstractParser::add_sample(VcfRecord& record, char* tab1, char* tab2)
{
    if(tab2 == tab1+1)
        cerr << "Warning: malformed VCF record line (empty SAMPLE field should be marked by a dot)."
             << " Line " << line_number_ << " in file " << path_ << "'.\n";

    record.samples_m().push_back(string(tab1+1, tab2));

    ++sample_index_;
}

inline
void VcfParser::getline(char* ptr, size_t n) {
    file_.getline(ptr, n);
    if(file_.eof()) {
        eof_=file_.fail();
        if(!eof_)
            cerr << "Notice: File '" << path_ << "' does not end with a newline.\n";
    }
}

#ifdef HAVE_LIBZ
inline
void VcfGzParser::getline(char* buf, size_t n) {
    if (gzgets(file_, buf, n) == NULL)
        eof_ = true;
}

inline
void VcfGzParser::check_eol() {
    if(*(tabs_.back()-1) == '\n') { // rem. safe, gzgets returns NULL on EOF
        eol_ = true;
        tabs_.back() = tabs_.back()-1;
        *tabs_.back() = '\0';
    } else {
        eol_= false;
    }
}
#endif // HAVE_LIBZ

#endif // __VCF_H__
