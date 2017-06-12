#include <algorithm>
#include <ctime>
#include <iomanip>
#include <sstream>

#include "Vcf.h"

using namespace std;

void VcfRecord::assign(const char* rec, size_t len, const VcfHeader& header) {
    buffer_.resize(len+1);
    memcpy(buffer_.data(), rec, len+1);
    if (buffer_[len-1] == '\n') {
        buffer_[len-1] = '\0';
        if (buffer_[len-2] == '\r')
            buffer_[len-2] = '\0';
    }

    // Parse the record.
    size_t n_exp_fields = header.samples().empty() ? 8 : 9 + header.samples().size();
    auto wrong_n_fields = [&rec, &len, &n_exp_fields](size_t obs){
        cerr << "Error: Expected VCF record to have "
             << n_exp_fields << " fields, not " << obs << "."
             << "Offending record (" << len << " characters) is:\n"
             << rec << '\n';
        throw exception();
    };

    size_t n_fields = 1;
    char* p = buffer_.data();
    char* q;
    auto next_field = [&p, &q, &n_fields, &wrong_n_fields]() {
        q = p;
        p = strchr(p, '\t');
        if (!p)
            wrong_n_fields(n_fields);
        ++n_fields;
        *p = '\0';
        ++p;
    };

    // chrom
    next_field();

    // pos
    next_field();
    strings_.push_back(q - buffer_.data());

    // id
    next_field();
    strings_.push_back(q - buffer_.data());

    // alleles (ref + alt)
    next_field();
    strings_.push_back(q - buffer_.data());
    next_field();
    if (*q != '.') { // If ALT is '.', we don't record the pointer to it.
        strings_.push_back(q - buffer_.data());
        while(q = strchr(q, ',')) {
            *q = '\0';
            ++q;
            strings_.push_back(q - buffer_.data());
        }
    }

    // filters
    filters_i_ = strings_.size();
    next_field();
    strings_.push_back(q - buffer_.data());

    // qual
    next_field();
    strings_.push_back(q - buffer_.data());

    if (header.samples().empty()) {
        // info (last field)
        if (n_fields != n_exp_fields)
            wrong_n_fields(n_fields);
        strings_.push_back(q - buffer_.data());
        while(p = strchr(p, ':')) {
            *p = '\0';
            ++p;
            strings_.push_back(p - buffer_.data());
        }
        format0_i_ = SIZE_MAX;
        sample0_i_ = SIZE_MAX;
    } else {
        // info
        next_field();
        strings_.push_back(q - buffer_.data());
        while(q = strchr(q, ';')) {
            *q = '\0';
            ++q;
            strings_.push_back(q - buffer_.data());
        }

        // format
        format0_i_ = strings_.size();
        next_field();
        strings_.push_back(q - buffer_.data());
        while(q = strchr(q, ':')) {
            *q = '\0';
            ++q;
            strings_.push_back(q - buffer_.data());
        }

        // samples (last fields)
        sample0_i_ = strings_.size();
        strings_.push_back(p - buffer_.data());
        while(p = strchr(p, '\t')) {
            *p = '\0';
            ++p;
            ++n_fields;
            strings_.push_back(p - buffer_.data());
        }
        if (n_fields != n_exp_fields)
            wrong_n_fields(n_fields);
    }
}

Vcf::RType VcfRecord::type() const {
    if (n_alleles() == 1)
        return Vcf::RType::invariant;
    for (size_t i=0; i<n_alleles(); ++i) {
        if (strchr(allele(i), '<'))
            return Vcf::RType::symbolic;
        else if (strchr(allele(i), '[') || strchr(allele(i), ']'))
            return Vcf::RType::breakend;
    }
    return Vcf::RType::expl;
}

string VcfRecord::util::fmt_info_af(const vector<double>& alt_freqs) {
    stringstream ss;
    ss << std::setprecision(3);
    join(alt_freqs, ',', ss);
    return string("AF=") + ss.str();
}

string VcfRecord::util::fmt_gt_gl(const vector<Nt2>& alleles, const GtLiks& liks) {
    static const size_t n_digits = 5;

    assert(!alleles.empty());
    vector<double> v;
    v.reserve(n_genotypes(alleles.size()));
    for (size_t a2=0; a2<alleles.size(); ++a2) {
        Nt2 a2nt = alleles[a2];
        for (size_t a1=0; a1<=a2; ++a1) {
            Nt2 a1nt = alleles[a1];
            v.push_back(liks.at(a1nt, a2nt) / log(10));
        }
    }
    stringstream ss;
    ss << std::setprecision(n_digits);
    join(v, ',', ss);
    return ss.str();
}

GtLiks VcfRecord::util::parse_gt_gl(const vector<Nt2>& alleles, const string& gl) {

    auto illegal = [&gl](const string& msg = ""){
        cerr << "Error: Illegal VCF genotype GL field: '" << gl << "'" << msg << ".\n";
        throw exception();
    };

    GtLiks liks;

    if (gl == ".")
        return liks;
    else if (gl.empty())
        illegal(" (empty field)");

    vector<double> v;
    double d = 0;
    const char* p = gl.c_str();
    char* end = NULL;
    d = std::strtod(p, &end);
    if (end == p)
        illegal();
    v.push_back(d);
    p = end;
    while(*p != '\0') {
        if (*p != ',')
            illegal(" (expected a comma)");
        ++p;
        d = std::strtod(p, &end);
        if (end == p)
            illegal();
        v.push_back(d);
        p = end;
    }
    if (v.size() != n_genotypes(alleles.size()))
        illegal(string(" (expected ") + to_string(n_genotypes(alleles.size())) + " values)");

    size_t gt_i = 0;
    for (size_t a2=0; a2<alleles.size(); ++a2) {
        Nt2 a2nt = alleles[a2];
        for (size_t a1=0; a1<=a2; ++a1) {
            Nt2 a1nt = alleles[a1];
            liks.set(a1nt, a2nt, v[gt_i] * log(10)); // Back to base e.
            ++gt_i;
        }
    }
    assert(gt_i == v.size());

    return liks;
}

bool VcfRecord::util::build_haps(pair<string,string>& haplotypes, const vector<const VcfRecord*>& snp_records, size_t sample_index) {
    haplotypes.first.resize(snp_records.size(), 'N');
    haplotypes.second.resize(snp_records.size(), 'N');

    #ifdef DEBUG
    size_t phase_set (-1);
    #endif
    bool complete = true;
    for (size_t i=0; i<snp_records.size(); ++i) {
        // For each SNP...
        // Parse the genotype.
        const VcfRecord& rec = *snp_records[i];
        pair<int,int> gt = rec.parse_genotype_nochecks(rec.sample(sample_index));

        // Record the (phased) alleles.
        if (gt.first == -1) {
            complete = false;
            haplotypes.first[i] = 'N';
            haplotypes.second[i] = 'N';
            continue;
        }
        haplotypes.first[i] = rec.allele(gt.first)[0];
        haplotypes.second[i] = rec.allele(gt.second)[0];

        #ifdef DEBUG
        // Check that there is only one phase set.
        assert(rec.n_formats() >= 2 && strcmp(rec.format(1), "PS") == 0);
        if (gt.first != gt.second) {
            size_t ps = stoi(rec.parse_gt_subfield(rec.sample(sample_index), 1));
            if (phase_set == size_t(-1))
                phase_set = ps;
            else
                assert(ps == phase_set);
        }
        #endif
    }

    return complete;
}

const VcfMeta VcfMeta::predefs::info_AD ("INFO","<ID=AD,Number=R,Type=Integer,Description=\"Total Depth for Each Allele\">");
const VcfMeta VcfMeta::predefs::info_AF ("INFO","<ID=AF,Number=A,Type=Float,Description=\"Allele Frequency\">");
const VcfMeta VcfMeta::predefs::info_DP ("INFO","<ID=DP,Number=1,Type=Integer,Description=\"Total Depth\">");
const VcfMeta VcfMeta::predefs::info_NS ("INFO","<ID=NS,Number=1,Type=Integer,Description=\"Number of Samples With Data\">");

const VcfMeta VcfMeta::predefs::format_AD ("FORMAT","<ID=AD,Number=R,Type=Integer,Description=\"Allele Depth\">");
const VcfMeta VcfMeta::predefs::format_DP ("FORMAT","<ID=DP,Number=1,Type=Integer,Description=\"Read Depth\">");
const VcfMeta VcfMeta::predefs::format_HQ ("FORMAT","<ID=HQ,Number=2,Type=Integer,Description=\"Haplotype Quality\">");
const VcfMeta VcfMeta::predefs::format_GL ("FORMAT","<ID=GL,Number=G,Type=Float,Description=\"Genotype Likelihood\">");
const VcfMeta VcfMeta::predefs::format_GQ ("FORMAT","<ID=GQ,Number=1,Type=Integer,Description=\"Genotype Quality\">");
const VcfMeta VcfMeta::predefs::format_GT ("FORMAT","<ID=GT,Number=1,Type=String,Description=\"Genotype\">");

const VcfMeta VcfMeta::predefs::info_locori ("INFO","<ID=locori,Number=1,Type=Character,Description=\"Orientation the corresponding Stacks locus aligns in\">");

void VcfHeader::init_meta(const string& version) {
    add_meta(VcfMeta("fileformat", version));

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
    check_open(file_, path);
#if ZLIB_VERNUM >= 0x1240
    gzbuffer(file_, libz_buffer_size);
#endif
}

#endif // HAVE_LIBZ

unique_ptr<VcfAbstractParser>
Vcf::adaptive_open(const string& path)
{
    FileT filet = guess_file_type(path);
    if (filet == FileT::vcf) {
        return unique_ptr<VcfAbstractParser>(new VcfParser(path));
#ifdef HAVE_LIBZ
    } else if (filet == FileT::gzvcf) {
        return unique_ptr<VcfAbstractParser>(new VcfGzParser(path));
#endif
    } else {
        cerr << "Error: File '" << path << "' : expected '.vcf(.gz)' suffix.\n";
        throw exception();
    }
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
    record.chrom_m().assign(line_, tabs_[Vcf::chrom]);

    // POS
    if(tabs_[Vcf::pos-1]+1 == tabs_[Vcf::pos]
        || *(tabs_[Vcf::pos-1]+1) == '.' ){
        cerr << "Warning: Skipping malformed VCF record (missing POS value)."
                << " Line " << line_number_ << " in file '" << path_ << "'.\n";
        // Skip this record.
        // ([record.type] is still [null].)
        return true;
    }
    record.pos_m() = stoi(tabs_[Vcf::pos-1]+1) -1 ; // VCF is 1-based

    // ID
    if(tabs_[Vcf::id-1]+1 == tabs_[Vcf::id])
        cerr << "Notice: Empty ID field should be marked by a dot."
                << " Line " << line_number_ << " in file '" << path_ << "'.\n";
    if(*(tabs_[Vcf::id-1]+1) != '.')
        record.id_m().assign(tabs_[Vcf::id-1]+1, tabs_[Vcf::id]);

    // REF
    if(tabs_[Vcf::ref-1]+1 == tabs_[Vcf::ref]
        || *(tabs_[Vcf::ref-1]+1) == '.' ) {
        cerr << "Warning: malformed VCF record (missing REF value)."
                << " Line " << line_number_ << " in file '" << path_ << "'.\n";
        // Skip this record.
        // ([record.type] is still [null].)
        return true;
    }
    record.alleles_m().push_back(string(tabs_[Vcf::ref-1]+1, tabs_[Vcf::ref]));

    // ALT & determine the type of the record
    if(tabs_[Vcf::alt-1]+1 == tabs_[Vcf::alt]) {
        cerr << "Warning: Skipping malformed VCF record (expected ALT field to be marked by a dot)."
                << " Line " << line_number_ << " in file '" << path_ << "'.\n";
        // Skip this record.
        // ([record.type] is still [null].)
        return true;
    } else if(*(tabs_[Vcf::alt-1]+1) == '.') {
        record.type_m() = Vcf::RType::invariant;
    } else {
        get_bounds(bounds_, tabs_[Vcf::alt-1], tabs_[Vcf::alt], ',');
        for (size_t i = 0; i < bounds_.size()-1; ++i ) {
            record.alleles_m().push_back(string(bounds_.at(i)+1, bounds_.at(i+1)));
            if (record.alleles_m().back().empty()) {
                record.alleles_m().pop_back();
                cerr << "Warning: Skipping malformed VCF record (malformed ALT field)."
                     << " Line " << line_number_ << " in file '" << path_ << "'.\n";
                // Skip this record.
                // ([record.type] is still [null].)
                return true;
            }
        }

        if (strchr(tabs_[Vcf::alt-1]+1, '[') || strchr(tabs_[Vcf::alt-1]+1, ']'))
            record.type_m() = Vcf::RType::breakend;
        else if (strchr(tabs_[Vcf::alt-1]+1, '<'))
            record.type_m() = Vcf::RType::symbolic;
        else
            record.type_m() = Vcf::RType::expl;
    }

    // Do not parse symbolic & breakend records.
    if (record.type() == Vcf::RType::symbolic || record.type() == Vcf::RType::breakend) {
        read_to_eol();
        return true;
    }

    // QUAL
    if(tabs_[Vcf::qual-1]+1 == tabs_[Vcf::qual])
        cerr << "Notice: Empty QUAL field should be marked by a dot."
             << " Line " << line_number_ << " in file '" << path_ << "'.\n";
    if(*(tabs_[Vcf::qual-1]+1) != '.')
        record.qual_m().assign(tabs_[Vcf::qual-1]+1, tabs_[Vcf::qual]);

    // FILTER
    if(tabs_[Vcf::filter-1]+1 == tabs_[Vcf::filter])
        cerr << "Notice: Empty FILTER field should be marked by a dot."
             << " Line " << line_number_ << " in file '" << path_ << "'.\n";
    if(*(tabs_[Vcf::filter-1]+1) != '.') {
        get_bounds(bounds_, tabs_[Vcf::filter-1], tabs_[Vcf::filter],';');
        for(size_t i = 0; i < bounds_.size()-1; ++i )
            record.filter_m().push_back(string(bounds_.at(i)+1, bounds_.at(i+1)));
    }

    // INFO
    if(tabs_[Vcf::info-1]+1 == tabs_[Vcf::info])
        cerr << "Notice: Empty INFO field should be marked by a dot."
             << " Line " << line_number_ << " in file '" << path_ << "'.\n";
    if(*(tabs_[Vcf::info-1]+1) != '.') {
        get_bounds(bounds_, tabs_[Vcf::info-1], tabs_[Vcf::info], ';');
        for(size_t i = 0; i < bounds_.size()-1; ++i )
                record.info_m().push_back(string(bounds_[i]+1, bounds_[i+1]));
    }

    // FORMAT
    if (!header_.samples().empty() && *(tabs_[Vcf::format-1]+1) != '.') {
        if(tabs_[Vcf::format-1]+1 == tabs_[Vcf::format])
            cerr << "Notice: Empty FORMAT field should be marked by a dot."
                 << " Line " << line_number_ << " in file '" << path_ << "'.\n";

        get_bounds(bounds_, tabs_[Vcf::format-1], tabs_[Vcf::format],':');
        for(size_t i = 0; i < bounds_.size()-1; ++i )
            record.format_m().push_back(string(bounds_.at(i)+1, bounds_.at(i+1)));

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
                record.type_m() = Vcf::RType::null;
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

void VcfWriter::write_header() {
    for(const VcfMeta& m : header_.meta())
        file_ << "##" << m.key() << "=" << m.value() << "\n";

    file_ << Vcf::base_header;
    if(not header_.samples().empty())
        file_ << "\tFORMAT";
    for(const string& s : header_.samples())
        file_ << "\t" << s;
    file_ << "\n";
}

void VcfWriter::write_record(const VcfRecord& r) {

    if (r.type() != Vcf::RType::expl)
        // This is not implemented.
        throw exception();

    file_ << r.chrom()
          << "\t" << r.pos()
          << "\t" << r.id()
          << "\t" << r.allele(0);

    //ALT
    file_ << "\t";
    if (r.n_alleles() == 1) {
        file_ << ".";
    } else {
        file_ << r.allele(1);
        for (size_t i=2; i<r.n_alleles(); ++i)
            file_ << "," << r.allele(i);
    }

    //QUAL
    file_ << "\t" << r.qual();

    //FILTER
    file_ << "\t" << r.filters();

    //INFO
    file_ << "\t";
    if (r.n_infos() == 0) {
        file_ << ".";
    } else {
        file_ << r.info(0);
        for (size_t i=1; i<r.n_infos(); ++i)
            file_ << ";" << r.info(i);
    }

    if (not header_.samples().empty()) {
        //FORMAT
        file_ << "\t";
        if (r.n_formats() == 0) {
            file_ << ".";
        } else {
            file_ << r.format(0);
            for (size_t i=1; i<r.n_formats(); ++i)
                file_ << ":" << r.format(i);
        }

        //SAMPLES
        assert(r.n_samples() != header_.samples().size());
        for (size_t i=0; i<r.n_samples(); ++i)
            file_ << "\t" << r.sample(i);
    }

    file_ << "\n";
}
