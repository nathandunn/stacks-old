#include <algorithm>

#include "Vcf.h"

using namespace std;

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

void
malformed_header(const string& path, const size_t line_no)
{
    cerr << "Error: Malformed VCF header. Line " << line_no << " in file '" << path << "'.\n";
    throw exception();
}

VcfAbstractParser::VcfAbstractParser()
: path_(), header_(), header_lines_(0), format_fields_to_keep_(), samples_to_keep_(),
  line_number_(0), eol_(true), eof_(false), tabs_(), bounds_(), kept_format_fields_(),
  sample_index_(-1)
{
    memset(line_, '\0', Vcf::line_buf_size);
}

void
VcfAbstractParser::read_header()
{
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
                read_to_eol();
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

    if(strncmp(line_, "#CHROM\t", 7) != 0)
        malformed_header(path_, line_number_);

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
                cerr << "Error: VCF file " << path_ << " does not end with a newline.\n";
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

inline void
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

    if (samples_to_keep_.empty() || samples_to_keep_.at(sample_index_)) {
        record.samples.push_back(vector<string>());
        vector<string>& sample = record.samples.back();
        get_bounds(bounds_, tab1, tab2, ':');
        if(bounds_.size()-1 > kept_format_fields_.size()) {
            // n.b. kept_format_fields_.size() is the total number of format fields
            // for this record, whereas record.format.size() is the number of fields
            // that were retained.
            cerr << "Error: malformed VCF genotype : '" << string(tab1+1, tab2)
                 << "' (at most " << kept_format_fields_.size() << " subfields exepcted)."
                 << " Line " << line_number_ << " in file " << path_ << "'.\n";
            throw exception ();
        }
        for(size_t i = 0; i<bounds_.size()-1; ++i)
            if(kept_format_fields_[i])
                sample.push_back(string(bounds_[i]+1, bounds_.at(i+1)));

        if(sample.size() < record.format.size()) {
            while(sample.size() < record.format.size())
                sample.push_back(string()); // We add missing values for the remaining subfields.
        }
    }
    ++sample_index_;
}

bool
VcfParser::open(const string& path)
{
    file_.exceptions(ifstream::badbit);

    file_.open(path);
    if(file_.fail()) {
        return false;
    }

    read_header();
    return true;
}

#ifdef HAVE_LIBZ
bool
VcfGzParser::open(const string& path)
{
    file_ = gzopen(path.c_str(), "rb");
    if(!file_) {
        return false;
    }
#if ZLIB_VERNUM >= 0x1240
    gzbuffer(file_, libz_buffer_size);
#endif

    read_header();
    return true;
}

void VcfAbstractParser::samples_to_keep(const set<string>& samples) {
    samples_to_keep_.clear();
    samples_to_keep_.reserve(header_.samples().size());
    for (vector<string>::const_iterator
            s=header_.samples().begin();
            s!=header_.samples().end();
            ++s)
        samples_to_keep_.push_back(bool(samples.count(*s)));
}

inline void
VcfGzParser::check_eol()
{
    if(*(tabs_.back()-1) == '\n') { // n.b. safe; gzgets never returns a null string
        eol_ = true;
        tabs_.back() = tabs_.back()-1;
        *tabs_.back() = '\0';
    } else {
        eol_= false;
    }
}
#endif

VcfAbstractParser*
Vcf::adaptive_open(const string& path)
{
    VcfAbstractParser* parser = NULL;
    if (path.length() >= 4 && path.substr(path.length()-4) == ".vcf") {
        parser = new VcfParser();
        if (not parser->open(path)) {
            // Opening failed
            delete parser;
            parser = NULL;
        }
#ifdef HAVE_LIBZ
    } else if (path.length() >= 7 && path.substr(path.length()-7) == ".vcf.gz") {
        parser = new VcfGzParser();
        if (not parser->open(path)) {
            delete parser;
            parser = NULL;
        }
#endif
    } else {
        cerr << "Error: File '" << path << "' : expected '.vcf(.gz)' suffix.";
        throw exception();
    }

    return parser;
}

bool
VcfAbstractParser::next_record(VcfRecord& record)
{
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
                    << " fields required).\nLine " << line_number_ << " in file '" << path_ << "'.\n";
            throw exception ();
        }
    } else {
        if(header_.samples().empty() || tabs_.size() < Vcf::base_fields_no + 2) {
            read_to_eol();
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

    // CHROM
    if(line_ == tabs_[Vcf::chrom]) {
        cerr << "Error: malformed VCF record line (CHROM field is required)."
                << " Line " << line_number_ << " in file '" << path_ << "'.\n";
        throw exception();
    }
    record.chrom.assign(line_, tabs_[Vcf::chrom]);

    // POS
    if(tabs_[Vcf::pos-1]+1 == tabs_[Vcf::pos]) {
        cerr << "Error: malformed VCF record line (POS field is required)."
                << " Line " << line_number_ << " in file '" << path_ << "'.\n";
        throw exception();
    }
    record.pos = atoi(tabs_[Vcf::pos-1]+1) -1 ; // VCF is 1-based

    // ID
    if(tabs_[Vcf::id-1]+1 == tabs_[Vcf::id])
        cerr << "Warning: malformed VCF record line (empty ID field should be marked by a dot)."
                << " Line " << line_number_ << " in file '" << path_ << "'.\n";
    if(*(tabs_[Vcf::id-1]+1) != '.')
        record.id.assign(tabs_[Vcf::id-1]+1, tabs_[Vcf::id]);

    // REF
    if(tabs_[Vcf::ref-1]+1 == tabs_[Vcf::ref]) {
        cerr << "Error: malformed VCF record line (REF field is required)."
                << " Line " << line_number_ << " in file '" << path_ << "'.\n";
        throw exception();
    }
    record.alleles.push_back(string(tabs_[Vcf::ref-1]+1, tabs_[Vcf::ref]));

    // ALT & determine the type of the record
    if(*(tabs_[Vcf::alt-1]+1) == '\0') {
        cerr << "Warning: malformed VCF record line (empty ALT field should be marked by a dot)."
                << " Line " << line_number_ << " in file '" << path_ << "'.\n";
        record.type = Vcf::RType::invariant;
    } else if(*(tabs_[Vcf::alt-1]+1) == '.') {
        record.type = Vcf::RType::invariant;
    } else {
        get_bounds(bounds_, tabs_[Vcf::alt-1], tabs_[Vcf::alt], ',');
        for (size_t i = 0; i < bounds_.size()-1; ++i )
            record.alleles.push_back(string(bounds_.at(i)+1, bounds_.at(i+1)));

        if (strchr(tabs_[Vcf::alt-1]+1, '[') || strchr(tabs_[Vcf::alt-1]+1, ']'))
            record.type = Vcf::RType::breakend;
        else if (strchr(tabs_[Vcf::alt-1]+1, '<'))
            record.type = Vcf::RType::symbolic;
        else
            record.type = Vcf::RType::expl;
    }

    // Do not parse symbolic & breakend records.
    if (record.type == Vcf::RType::symbolic || record.type == Vcf::RType::breakend) {
        read_to_eol();
        return true;
    }

    // QUAL
    if(tabs_[Vcf::qual-1]+1 == tabs_[Vcf::qual])
        cerr << "Warning: malformed VCF record line (empty QUAL field should be marked by a dot)."
                << " Line " << line_number_ << " in file '" << path_ << "'.\n";
    if(*(tabs_[Vcf::qual-1]+1) != '.')
        record.qual.assign(tabs_[Vcf::qual-1]+1, tabs_[Vcf::qual]);

    // FILTER
    if(tabs_[Vcf::filter-1]+1 == tabs_[Vcf::filter])
        cerr << "Warning: malformed VCF record line (empty FILTER field should be marked by a dot)."
                << " Line " << line_number_ << " in file '" << path_ << "'.\n";
    if(*(tabs_[Vcf::filter-1]+1) != '.') {
        get_bounds(bounds_, tabs_[Vcf::filter-1], tabs_[Vcf::filter],';');
        for(size_t i = 0; i < bounds_.size()-1; ++i )
            record.filter.push_back(string(bounds_.at(i)+1, bounds_.at(i+1)));
    }

    // INFO
    if(tabs_[Vcf::info-1]+1 == tabs_[Vcf::info])
        cerr << "Warning: malformed VCF record line (empty INFO field should be marked by a dot)."
        << " Line " << line_number_ << " in file '" << path_ << "'.\n";
    if(*(tabs_[Vcf::info-1]+1) != '.') {
        get_bounds(bounds_, tabs_[Vcf::info-1], tabs_[Vcf::info], ';');
        for(size_t i = 0; i < bounds_.size()-1; ++i ) {
            char* equal = strchr(bounds_[i]+1, '=');
            if(equal == NULL or bounds_[i+1] - equal < 0)
                record.info.push_back(pair<string,string>(string(bounds_[i]+1, bounds_[i+1]),string()));
            else
                record.info.push_back(pair<string,string>(string(bounds_[i]+1, equal),string(equal+1,bounds_[i+1])));
        }
    }

    // FORMAT
    if (!header_.samples().empty() && *(tabs_[Vcf::format-1]+1) != '.') {
        if(tabs_[Vcf::format-1]+1 == tabs_[Vcf::format])
            cerr << "Warning: malformed VCF record line (empty FORMAT field should be marked by a dot)."
                    << " Line " << line_number_ << " in file '" << path_ << "'.\n";

        get_bounds(bounds_, tabs_[Vcf::format-1], tabs_[Vcf::format],':');
        kept_format_fields_.clear();
        for(size_t i = 0; i < bounds_.size()-1; ++i ) {
            string format = string(bounds_.at(i)+1, bounds_.at(i+1));
            if(format_fields_to_keep_.empty() || format_fields_to_keep_.count(format)) {
                kept_format_fields_.push_back(true);
                record.format.push_back(format);
            } else {
                kept_format_fields_.push_back(false);
            }
        }
        if(record.format[Vcf::gt_subfield] == "GT")
            // There are genotypes for this record.
            record.no_gt = false;

        // SAMPLES
        sample_index_ = 0;
        for (size_t s = Vcf::base_fields_no; s < tabs_.size()-2; ++s) //n.b. -2
            add_sample(record, tabs_[s], tabs_[s+1]);

        // If the record line was not read entirely, copy the trucated sample field
        // to the beggining of the buffer and read some more.
        while (!eol_) {
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

            // Windows line endings.
            if(*(tabs_.back()-1) == '\r') {
                tabs_.back() = tabs_.back()-1;
                *tabs_.back() = '\0';
            }
            for(size_t i = 0; i<tabs_.size()-2;++i) {
                *tabs_[i+1] = '\0';
                add_sample(record, tabs_[i], tabs_[i+1]);
            }
        }
        //the actual last sample
        add_sample(record, *(tabs_.end()-2), tabs_.back());

        if(sample_index_ != header_.samples().size()) {
            cerr << "Error: malformed VCF record line ("
                 << header_.samples().size() << " SAMPLE fields expected, "
                 << sample_index_ <<  " found). File '" << path_ << "', line "
                 << line_number_ << " :\n" << line_ << "\n";
            throw exception();
        }
    }

    return true;
}
