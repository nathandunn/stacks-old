#include <iostream>
#include <cstring>
#include <cstdlib>
#include <exception>
#include <algorithm>

#include "Vcf.h"

using namespace std;

inline void get_bounds(vector<char*>& bounds, char* tab1, char* tab2, char sep) {
    bounds.clear();
    bounds.push_back(tab1);
    do {
        bounds.push_back(strchr(bounds.back()+1, sep)); // this is last+1 if tab1==tab2
    } while (bounds.back());
    bounds.pop_back();
    bounds.push_back(tab2);
}

void malformed_header(const string& path, const size_t line_no) {
    cerr << "Error: Malformed VCF header."
         << " Line " << line_no << " in file " << path
         << endl;
    throw exception();
}

inline void VcfRecord::clear() {
    pos_ = -1;
    chrom_.clear();
    id_.clear();
    ref_.clear();
    type_ = Vcf::null;
    alt_.clear();
    qual_.clear();
    filter_.clear();
    info_.clear();
    format_.clear();
    samples_.clear();
}

VcfAbstractParser::VcfAbstractParser()
: line_number_(0), header_lines_(0), eol_(true), eof_(false) {
    memset(line_, '\0', Vcf::line_buf_size);
}

void VcfAbstractParser::read_header() {
    while(true) {
        getline(line_, Vcf::line_buf_size);
        ++line_number_;
        if(eof_)
            malformed_header(path_, line_number_);
        if (line_[0] != '#')
            malformed_header(path_, line_number_);

        if(line_[1] == '#') {
            // Meta header line.

            char* equal = strchr(line_, '=');
            if(equal == NULL)
                malformed_header(path_, line_number_);
            char* end = equal + strlen(equal);

            // Check the length of the line (discard lines that don't fit in the buffer).
            tabs_.clear();
            tabs_.push_back(end);
            check_eol();
            if(!eol_) {
                read_while_not_eol();
                continue;
            }

            string key = string(line_+2, equal);
            transform(key.begin(), key.end(), key.begin(), ::toupper);
            if(key == "VERSION") {
                header_.version(string(equal+1, end));
            } else if(key == "FILEDATE") {
                header_.date(string(equal+1, end));
            } else if(key == "SOURCE") {
                header_.source(string(equal+1, end));
            } else if(key == "REFERENCE") {
                header_.reference(string(equal+1, end));
            } else if(key == "CONTIG") {
                ; // not implemented yet
            } else if(key == "ALT") {
                ; // not implemented yet
            } else if(key == "INFO") {
                ; // not implemented yet
            } else if(key == "FORMAT") {
                ; // not implemented yet
            } else {
                ; // not implemented yet ; skip or map<key,string>
            }

            continue;
        }
        break;
    }

    // Final header line.

    if(strncmp(line_, "#CHROM\t", 7) != 0) {
        cerr << "strncmp" << endl;
        malformed_header(path_, line_number_);
    }

    // Get tabs.
    tabs_.clear();
    tabs_.push_back(strchr(line_, '\t')); //rem. one tab guaranteed
    while (tabs_.back())
        tabs_.push_back(strchr(tabs_.back()+1, '\t'));
    tabs_.pop_back();
    tabs_.push_back(tabs_.back() + strlen(tabs_.back())); // final '\0'
    check_eol();

    // Check the number of tabs.
    if(!(tabs_.size() == Vcf::base_fields_no or tabs_.size() >= Vcf::base_fields_no+2))
        malformed_header(path_, line_number_);

    // Check windows line endings.
    if(*(tabs_.back()-1) == '\r') {
        tabs_.back() = tabs_.back()-1;
        *tabs_.back() = '\0';
    }

    // Parse sample names.
    if(tabs_.size() >= Vcf::base_fields_no+2) {
        for(size_t i=Vcf::first_sample-1; i < tabs_.size()-2; ++i)
            header_.add_sample(string(tabs_[i]+1, tabs_[i+1]));

        while(!eol_){
            // Buffer wasn't long enough : parse until eol is found.
            // Copy the truncated name at the beginning of the buffer.
            size_t last_field_len = tabs_.back() - *(tabs_.end()-2);
            *line_ = '\t';
            strcpy(line_+1, *(tabs_.end()-2)+1);

            getline(line_+last_field_len, Vcf::line_buf_size-last_field_len);
            if(eof_) {
                cerr << "Error: VCF file " << path_ << " does not end with a newline." << endl;
                throw exception();
            }

            // Get tabs.
            tabs_.clear();
            tabs_.push_back(line_);
            while (tabs_.back())
                tabs_.push_back(strchr(tabs_.back()+1, '\t'));
            tabs_.pop_back();
            tabs_.push_back(tabs_.back() + strlen(tabs_.back()));
            check_eol();

            // Check windows line endings.
            if(*(tabs_.back()-1) == '\r') {
                tabs_.back() = tabs_.back()-1;
                *tabs_.back() = '\0';
            }

            for(size_t i = 0; i<tabs_.size()-2;++i)
                header_.add_sample(string(tabs_[i]+1, tabs_[i+1]));
        }
        header_.add_sample(string(*(tabs_.end()-2)+1, tabs_.back()));

        if(find(header_.samples().begin(),header_.samples().end(), string()) != header_.samples().end())
            // empty field
            malformed_header(path_, line_number_);
    }

    header_lines_ = line_number_;
    return;
}

inline void VcfAbstractParser::read_while_not_eol() {
    while(!eol_) {
        getline(line_, Vcf::line_buf_size);
        if(eof_) {
            cerr << "Error: VCF file " << path_ << " does not end with a newline." << endl;
            throw exception();
        }
        tabs_.clear();
        tabs_.push_back(line_+strlen(line_));
        check_eol();
    }
}

inline void VcfAbstractParser::add_sample(VcfRecord& record, char* tab1, char* tab2) {
    if(tab2 == tab1+1)
        cerr << "Warning: malformed VCF record line (empty SAMPLE field should be marked by a dot)."
             << " Line " << line_number_ << " in file " << path_
             << endl;

    record.samples_.push_back(vector<string>());
    vector<string>& sample = record.samples_.back();
    get_bounds(bounds_, tab1, tab2, ':');
    if(bounds_.size()-1 > kept_format_fields_.size()) {
        cerr << "Error: malformed VCF genotype : '" << string(tab1+1, tab2)
                << "' (at most " << kept_format_fields_.size() << " subfields exepcted)."
                << " Line " << line_number_ << " in file " << path_
                << endl;
        throw exception ();
    }
    for(size_t i = 0; i<bounds_.size()-1; ++i)
        if(kept_format_fields_[i])
            sample.push_back(string(bounds_[i]+1, bounds_.at(i+1)));

    if(sample.size() < record.format_.size()) {
        while(sample.size() < record.format_.size())
            sample.push_back(string()); // We add missing values for the remaining subfields.
    }
}

int VcfParser::open(const string& path) {
    file_.exceptions(ifstream::badbit);

    file_.open(path);
    if(!file_.good()) {
        return 1;
    }

    read_header();
    return 0;
}

#ifdef HAVE_LIBZ
int VcfGzParser::open(const string& path) {
    file_ = gzopen(path.c_str(), "rb");
    if(!file_) {
        return 1;
    }
#if ZLIB_VERNUM >= 0x1240
    gzbuffer(file_, libz_buffer_size);
#endif

    read_header();
    return 0;
}

inline void VcfGzParser::check_eol() {
    if(*(tabs_.back()-1) == '\n') { // n.b. safe; gzgets never returns a null string
        eol_ = true;
        tabs_.back() = tabs_.back()-1;
        *tabs_.back() = '\0';
    } else {
        eol_= false;
    }
}
#endif

bool VcfAbstractParser::next_record(VcfRecord& record) {
    getline(line_, Vcf::line_buf_size);
    ++line_number_;
    if(eof_)
        return false;

    record.clear();

    // Get tabs.
    tabs_.clear();
    tabs_.push_back(strchr(line_, '\t'));
    while(tabs_.back())
        tabs_.push_back(strchr(tabs_.back()+1, '\t'));
    tabs_.pop_back();
    if(tabs_.size() > 0)
        tabs_.push_back(tabs_.back() + strlen(tabs_.back()));
    else
        tabs_.push_back(line_ + strlen(line_));
    check_eol();

    /* Check the number of fields (should be ==8 or >=10 depending on the
     * presence of samples) :
     * If the line fits in the buffer but the number of fields is wrong,
     * raise an error.
     * If the fixed fields don't fit in the buffer, return a null record. */
    if(eol_) {
        if(tabs_.size() < Vcf::base_fields_no + (header_.samples().empty() ? 0 : 2)) {
            cerr << "Error: malformed VCF record line (at least "
                    << Vcf::base_fields_no + (header_.samples().empty() ? 0 : 2)
                    << " fields required).\nLine " << line_number_ << " in file " << path_
                    << endl;
            throw exception ();
        }
    } else {
        if(header_.samples().empty() || tabs_.size() < Vcf::base_fields_no + 2) {
            read_while_not_eol();
            return true;
        }
    }

    // Check windows line endings.
    if(*(tabs_.back()-1) == '\r') {
        tabs_.back() = tabs_.back()-1;
        *tabs_.back() = '\0';
    }

    // Separate all substrings (replace \t's with \0's).
    for(size_t i = 0; i<tabs_.size()-1;++i)
        *tabs_[i] = '\0';

    //chrom
    if(line_ == tabs_[Vcf::chrom]) {
        cerr << "Error: malformed VCF record line (CHROM field is required)."
                << " Line " << line_number_ << " in file " << path_
             << endl;
        throw exception();
    }
    record.chrom_.assign(line_, tabs_[Vcf::chrom]);

    //pos
    if(tabs_[Vcf::pos-1]+1 == tabs_[Vcf::pos]) {
        cerr << "Warning: malformed VCF record line (POS field is required)."
                << " Line " << line_number_ << " in file " << path_
             << endl;
        throw exception();
    }
    record.pos_ = strtol(tabs_[Vcf::pos-1]+1, NULL, 10);

    //id
    if(tabs_[Vcf::id-1]+1 == tabs_[Vcf::id])
        cerr << "Warning: malformed VCF record line (empty ID field should be marked by a dot)."
                << " Line " << line_number_ << " in file " << path_
             << endl;
    if(*(tabs_[Vcf::id-1]+1) != '.')
        record.id_.assign(tabs_[Vcf::id-1]+1, tabs_[Vcf::id]);

    //ref
    if(tabs_[Vcf::ref-1]+1 == tabs_[Vcf::ref]) {
        cerr << "Warning: malformed VCF record line (REF field is required)."
                << " Line " << line_number_ << " in file " << path_
             << endl;
        throw exception();
    }
    record.ref_.assign(tabs_[Vcf::ref-1]+1, tabs_[Vcf::ref]);

    //alt
    if(tabs_[Vcf::alt-1]+1 == tabs_[Vcf::alt])
        cerr << "Warning: malformed VCF record line (empty ALT field should be marked by a dot)."
                << " Line " << line_number_ << " in file " << path_
             << endl;
    if(*(tabs_[Vcf::alt-1]+1) == '<') {
        record.type_ = Vcf::symbolic;
        read_while_not_eol();
        return true; // Do not parse records describing symbolic alleles.
    }
    if (strchr(tabs_[Vcf::alt-1]+1, '[') || strchr(tabs_[Vcf::alt-1]+1, ']')) {
        record.type_ = Vcf::breakend;
        read_while_not_eol();
        return true; // Do not parse records describing breakend alleles.
    }
    record.type_ = Vcf::expl;
    if(*(tabs_[Vcf::alt-1]+1) != '.') {
        get_bounds(bounds_, tabs_[Vcf::alt-1], tabs_[Vcf::alt], ',');
        for(size_t i = 0; i < bounds_.size()-1; ++i )
            record.alt_.push_back(string(bounds_.at(i)+1, bounds_.at(i+1)));
    }

    //qual
    if(tabs_[Vcf::qual-1]+1 == tabs_[Vcf::qual])
        cerr << "Warning: malformed VCF record line (empty QUAL field should be marked by a dot)."
                << " Line " << line_number_ << " in file " << path_
             << endl;
    if(*(tabs_[Vcf::qual-1]+1) != '.')
        record.qual_.assign(tabs_[Vcf::qual-1]+1, tabs_[Vcf::qual]);

    //filter
    if(tabs_[Vcf::filter-1]+1 == tabs_[Vcf::filter])
        cerr << "Warning: malformed VCF record line (empty FILTER field should be marked by a dot)."
                << " Line " << line_number_ << " in file " << path_
             << endl;
    if(*(tabs_[Vcf::filter-1]+1) != '.') {
        get_bounds(bounds_, tabs_[Vcf::filter-1], tabs_[Vcf::filter],';');
        for(size_t i = 0; i < bounds_.size()-1; ++i )
            record.filter_.push_back(string(bounds_.at(i)+1, bounds_.at(i+1)));
    }

    //info
    if(tabs_[Vcf::info-1]+1 == tabs_[Vcf::info])
        cerr << "Warning: malformed VCF record line (empty INFO field should be marked by a dot)."
        << " Line " << line_number_ << " in file " << path_
        << endl;
    if(*(tabs_[Vcf::info-1]+1) != '.') {
        get_bounds(bounds_, tabs_[Vcf::info-1], tabs_[Vcf::info], ';');
        for(size_t i = 0; i < bounds_.size()-1; ++i ) {
            char* equal = strchr(bounds_[i]+1, '=');
            if(equal == NULL or bounds_[i+1] - equal < 0)
                record.info_.push_back(pair<string,string>(string(bounds_[i]+1, bounds_[i+1]),string()));
            else
                record.info_.push_back(pair<string,string>(string(bounds_[i]+1, equal),string(equal+1,bounds_[i+1])));
        }
    }

    //format
    if (!header_.samples().empty() && *(tabs_[Vcf::format-1]+1) != '.') {
        if(tabs_[Vcf::format-1]+1 == tabs_[Vcf::format])
            cerr << "Warning: malformed VCF record line (empty FORMAT field should be marked by a dot)."
                    << " Line " << line_number_ << " in file " << path_
                 << endl;

        get_bounds(bounds_, tabs_[Vcf::format-1], tabs_[Vcf::format],':');
        kept_format_fields_.clear();
        for(size_t i = 0; i < bounds_.size()-1; ++i ) {
            string format = string(bounds_.at(i)+1, bounds_.at(i+1));
            if(format_fields_to_keep_.size() == 0 || format_fields_to_keep_.count(format)) {
                kept_format_fields_.push_back(true);
                record.format_.push_back(format);
            } else {
                kept_format_fields_.push_back(false);
            }
        }

        //samples
        for(size_t s = Vcf::base_fields_no; s < tabs_.size()-2; ++s)
            add_sample(record, tabs_[s], tabs_[s+1]);

        while(!eol_) {
            // Copy the trucated sample field to the beggining of the buffer and read some more.
            size_t lastfieldlen = tabs_.back() - *(tabs_.end()-2);
            *line_ = '\t';
            strcpy(line_+1, *(tabs_.end()-2)+1);

            getline(line_+lastfieldlen, Vcf::line_buf_size-lastfieldlen);
            if(eof_) {
                cerr << "Error: VCF file " << path_ << " does not end with a newline." << endl;
                throw exception();
            }
            // Get tabs.
            tabs_.clear();
            tabs_.push_back(line_);
            while (tabs_.back())
                tabs_.push_back(strchr(tabs_.back()+1, '\t'));
            tabs_.pop_back();
            tabs_.push_back(tabs_.back() + strlen(tabs_.back()));
            check_eol();

            // Check windows line endings.
            if(*(tabs_.back()-1) == '\r') {
                tabs_.back() = tabs_.back()-1;
                *tabs_.back() = '\0';
            }
            for(size_t i = 0; i<tabs_.size()-2;++i) {
                *tabs_[i+1] = '\0';
                add_sample(record, tabs_[i], tabs_[i+1]);
            }
        }
        add_sample(record, *(tabs_.end()-2), tabs_.back()); // last sample

        if(record.samples_.size() != header_.samples().size()) {
            cerr << "Error: malformed VCF record line ("
                 << header_.samples().size() << " SAMPLE fields expected, " << record.samples_.size() <<  " found)."
                 << " Line " << line_number_ << " in file " << path_
                 << endl;
            throw exception();
        }
    }

    return true;
}
