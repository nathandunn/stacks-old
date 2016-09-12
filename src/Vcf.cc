#include <algorithm>
#include <ctime>
#include <iomanip>
#include <sstream>

#include "Vcf.h"

using namespace std;

const map<string, VcfMeta> VcfMeta::predefined = {
        {"INFO/NS", VcfMeta("INFO","<ID=NS,Number=1,Type=Integer,Description=\"Number of Samples With Data\">")},
        {"INFO/AF", VcfMeta("INFO","<ID=AF,Number=.,Type=Float,Description=\"Allele Frequency\">")},
        {"FORMAT/GT", VcfMeta("FORMAT","<ID=GT,Number=1,Type=String,Description=\"Genotype\">")},
        {"FORMAT/DP", VcfMeta("FORMAT","<ID=DP,Number=1,Type=Integer,Description=\"Read Depth\">")},
        {"FORMAT/AD", VcfMeta("FORMAT","<ID=AD,Number=1,Type=Integer,Description=\"Allele Depth\">")},
        {"FORMAT/GL", VcfMeta("FORMAT","<ID=GL,Number=.,Type=Float,Description=\"Genotype Likelihood\">")},
    };

void VcfHeader::init_meta(const string& fileformat) {
    add_meta(VcfMeta("fileformat", fileformat));

    time_t t;
    time(&t);
    char date[9];
    strftime(date, 9, "%Y%m%d", localtime(&t));
    add_meta(VcfMeta("fileDate", date));

    add_meta(VcfMeta("source", string("\"Stacks v") + VERSION + "\""));
}

VcfAbstractParser::VcfAbstractParser(const string& path)
: path_(path), header_(), line_number_(0), eol_(true), eof_(false), tabs_(), bounds_(), sample_index_(-1)
{
    memset(line_, '\0', Vcf::line_buf_size);
}

void
VcfAbstractParser::read_header()
{
    auto malformed = [this] () {
        cerr << "Error: Malformed header."
             << " At line " << line_number_ << " in file '" << path_ << "'.\n";
        throw exception();
    };

    while(true) {
        getline(line_, Vcf::line_buf_size);
        ++line_number_;
        if(eof_)
            malformed();

        if (line_[0] != '#')
            malformed();

        if(line_[1] == '#') {
            // Meta header line.

            char* equal = strchr(line_, '=');
            if(equal == NULL) {
                // (Notice: Skipping missing '=')
                tabs_.clear();
                tabs_.push_back(line_+strlen(line_));
                check_eol();
                read_to_eol();
                continue;
            }
            char* end = equal + strlen(equal);

            // Check the length of the line (discard lines that don't fit in the buffer).
            tabs_.clear();
            tabs_.push_back(end);
            check_eol();
            if(!eol_) {
                // (Notice: Skipping long line)
                read_to_eol();
                continue;
            }

            string key = string(line_+2, equal);
            transform(key.begin(), key.end(), key.begin(), ::toupper);
            header_.add_meta(VcfMeta(key, string(equal+1, end)));

            continue;
        }
        break;
    }

    // Final header line.

    if(strncmp(line_, Vcf::base_fields.c_str(), Vcf::base_fields.length()) != 0)
        malformed();

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
        malformed();

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
            if(eof_)
                malformed();

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
    }

    return;
}

#ifdef HAVE_LIBZ
VcfGzParser::VcfGzParser(const string& path)
: VcfAbstractParser(path), file_(NULL)
{
    file_ = gzopen(path_.c_str(), "rb");
#if ZLIB_VERNUM >= 0x1240
    gzbuffer(file_, libz_buffer_size);
#endif
}

#endif // HAVE_LIBZ

VcfAbstractParser*
Vcf::adaptive_open(const string& path)
{
    VcfAbstractParser* parser = NULL;
    if (path.length() >= 4 && path.substr(path.length()-4) == ".vcf") {
        parser = new VcfParser(path);
        if (parser->fail()) {
            // Opening failed
            delete parser;
            parser = NULL;
        }
#ifdef HAVE_LIBZ
    } else if (path.length() >= 7 && path.substr(path.length()-7) == ".vcf.gz") {
        parser = new VcfGzParser(path);
        if (parser->fail()) {
            delete parser;
            parser = NULL;
        }
#endif
    } else {
        cerr << "Error: File '" << path << "' : expected '.vcf(.gz)' suffix.\n";
        throw exception();
    }

    return parser;
}

bool
VcfAbstractParser::next_record(VcfRecord& record)
{
    getline(line_, Vcf::line_buf_size);
    if(eof_)
        return false;
    ++line_number_;

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
            cerr << "Warning: In file '" << path_ << "': skipping the very long record at line "
                 << line_number_ << ".\n";
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
        cerr << "Warning: Skipping malformed VCF record (missing CHROM value)."
                << " Line " << line_number_ << " in file '" << path_ << "'.\n";
        // Skip this record.
        // ([record.type] is still [null].)
        return true;
    }
    record.chrom.assign(line_, tabs_[Vcf::chrom]);

    // POS
    if(tabs_[Vcf::pos-1]+1 == tabs_[Vcf::pos]
        || *(tabs_[Vcf::pos-1]+1) == '.' ){
        cerr << "Warning: Skipping malformed VCF record (missing POS value)."
                << " Line " << line_number_ << " in file '" << path_ << "'.\n";
        // Skip this record.
        // ([record.type] is still [null].)
        return true;
    }
    record.pos = stoi(tabs_[Vcf::pos-1]+1) -1 ; // VCF is 1-based

    // ID
    if(tabs_[Vcf::id-1]+1 == tabs_[Vcf::id])
        cerr << "Notice: Empty ID field should be marked by a dot."
                << " Line " << line_number_ << " in file '" << path_ << "'.\n";
    if(*(tabs_[Vcf::id-1]+1) != '.')
        record.id.assign(tabs_[Vcf::id-1]+1, tabs_[Vcf::id]);

    // REF
    if(tabs_[Vcf::ref-1]+1 == tabs_[Vcf::ref]
        || *(tabs_[Vcf::ref-1]+1) == '.' ) {
        cerr << "Warning: malformed VCF record (missing REF value)."
                << " Line " << line_number_ << " in file '" << path_ << "'.\n";
        // Skip this record.
        // ([record.type] is still [null].)
        return true;
    }
    record.alleles.push_back(string(tabs_[Vcf::ref-1]+1, tabs_[Vcf::ref]));

    // ALT & determine the type of the record
    if(tabs_[Vcf::alt-1]+1 == tabs_[Vcf::alt]) {
        cerr << "Warning: Skipping malformed VCF record (expected ALT field to be marked by a dot)."
                << " Line " << line_number_ << " in file '" << path_ << "'.\n";
        // Skip this record.
        // ([record.type] is still [null].)
        return true;
    } else if(*(tabs_[Vcf::alt-1]+1) == '.') {
        record.type = Vcf::RType::invariant;
    } else {
        get_bounds(bounds_, tabs_[Vcf::alt-1], tabs_[Vcf::alt], ',');
        for (size_t i = 0; i < bounds_.size()-1; ++i ) {
            record.alleles.push_back(string(bounds_.at(i)+1, bounds_.at(i+1)));
            if (record.alleles.back().empty()) {
                record.alleles.pop_back();
                cerr << "Warning: Skipping malformed VCF record (malformed ALT field)."
                     << " Line " << line_number_ << " in file '" << path_ << "'.\n";
                // Skip this record.
                // ([record.type] is still [null].)
                return true;
            }
        }

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
        cerr << "Notice: Empty QUAL field should be marked by a dot."
             << " Line " << line_number_ << " in file '" << path_ << "'.\n";
    if(*(tabs_[Vcf::qual-1]+1) != '.')
        record.qual.assign(tabs_[Vcf::qual-1]+1, tabs_[Vcf::qual]);

    // FILTER
    if(tabs_[Vcf::filter-1]+1 == tabs_[Vcf::filter])
        cerr << "Notice: Empty FILTER field should be marked by a dot."
             << " Line " << line_number_ << " in file '" << path_ << "'.\n";
    if(*(tabs_[Vcf::filter-1]+1) != '.') {
        get_bounds(bounds_, tabs_[Vcf::filter-1], tabs_[Vcf::filter],';');
        for(size_t i = 0; i < bounds_.size()-1; ++i )
            record.filter.push_back(string(bounds_.at(i)+1, bounds_.at(i+1)));
    }

    // INFO
    if(tabs_[Vcf::info-1]+1 == tabs_[Vcf::info])
        cerr << "Notice: Empty INFO field should be marked by a dot."
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
            cerr << "Notice: Empty FORMAT field should be marked by a dot."
                 << " Line " << line_number_ << " in file '" << path_ << "'.\n";

        get_bounds(bounds_, tabs_[Vcf::format-1], tabs_[Vcf::format],':');
        for(size_t i = 0; i < bounds_.size()-1; ++i )
            record.format.push_back(string(bounds_.at(i)+1, bounds_.at(i+1)));

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
                cerr << "Warning: Encountered end of file while reading a record.\n";
                record.type = Vcf::RType::null;
                return false;
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
            cerr << "Error: malformed VCF record ("
                 << header_.samples().size() << " SAMPLE fields expected, "
                 << sample_index_ <<  " found). File '" << path_ << "', line "
                 << line_number_ << ".\n";
            throw exception();
        }
    }

    return true;
}

void VcfWriter::write_header(const VcfHeader& header) {
    for(const VcfMeta& m : header.meta())
        file_ << "##" << m.key() << "=" << m.value() << "\n";

    file_ << Vcf::base_header;
    if(not header.samples().empty())
        file_ << "\tFORMAT";
    for(const string& s : header.samples())
        file_ << "\t" << s;
    file_ << "\n";
}

void VcfWriter::write_record(const VcfRecord& r, const VcfHeader& h) {

    if (r.type != Vcf::RType::expl)
        // This is not implemented.
        throw exception();

    file_ << r.chrom
          << "\t" << r.pos
          << "\t" << (r.id.empty() ? "." : r.id)
          << "\t" << r.alleles.at(0);

    //ALT
    file_ << "\t";
    if (r.alleles.size() == 1) {
        file_ << ".";
    } else {
        auto allele = r.alleles.begin()+1;
        file_ << *allele;
        ++allele;
        while(allele != r.alleles.end()) {
            file_ << "," << *allele;
            ++allele;
        }
    }

    file_ << "\t" << (r.qual.empty() ? "." : r.qual);

    //FILTER
    file_ << "\t";
    if (r.filter.empty()) {
        file_ << "\t.";
    } else {
        auto filter = r.filter.begin();
        file_ << *filter;
        ++filter;
        while(filter != r.filter.end()) {
            file_ << ";" << *filter;
            ++filter;
        }
    }

    //INFO
    file_ << "\t";
    if (r.info.empty()) {
        file_ << ".";
    } else {
        auto i = r.info.begin();
        file_ << i->first << "=" << i->second;
        ++i;
        while(i != r.info.end()) {
            file_ << ";" << i->first << "=" << i->second;
            ++i;
        }
    }

    if (not h.samples().empty()) {
        //FORMAT
        file_ << "\t";
        if (r.format.empty()) {
            file_ << ".";
        } else {
            auto f = r.format.begin();
            file_ << *f;
            ++f;
            while(f != r.format.end()) {
                file_ << ":" << *f;
                ++f;
            }
        }

        //SAMPLES
        if (r.samples.size() != h.samples().size())
            throw exception();

        for (const string& s : r.samples)
            file_ << "\t" << s;
    }

    file_ << "\n";
}
