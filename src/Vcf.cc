#include <algorithm>
#include <ctime>
#include <iomanip>
#include <sstream>

#include "Vcf.h"

using namespace std;

const string VcfHeader::std_fields = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO";

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
        while((q = strchr(q, ','))) {
            *q = '\0'; // ALT fields become null-separated.
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
        while((p = strchr(p, ';'))) {
            *p = '\0'; // INFO fields become null-separated.
            ++p;
            strings_.push_back(p - buffer_.data());
        }
        format0_i_ = SIZE_MAX;
        sample0_i_ = SIZE_MAX;
    } else {
        // info
        next_field();
        strings_.push_back(q - buffer_.data());
        while((q = strchr(q, ';'))) {
            *q = '\0';
            ++q;
            strings_.push_back(q - buffer_.data());
        }

        // format
        format0_i_ = strings_.size();
        next_field();
        strings_.push_back(q - buffer_.data());
        while((q = strchr(q, ':'))) {
            *q = '\0'; // FORMAT fields become null-separated.
            ++q;
            strings_.push_back(q - buffer_.data());
        }

        // samples (last fields)
        sample0_i_ = strings_.size();
        strings_.push_back(p - buffer_.data());
        while((p = strchr(p, '\t'))) {
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

VcfParser::VcfParser(const string& path)
: path_(path), header_(), line_number_(0),
  is_gzipped_(false), ifs_(), ifsbuffer_(),
  gzfile_(NULL), gzbuffer_(NULL), gzbuffer_size_(0), gzline_len_(0)
{
    FileT ftype = guess_file_type(path_);
    if (ftype == FileT::vcf) {
        is_gzipped_ = false;
        ifs_.open(path_);
        check_open(ifs_, path_);
    } else if (ftype == FileT::gzvcf) {
        is_gzipped_ = true;
        gzfile_ = gzopen(path_.c_str(), "rb");
        check_open(gzfile_, path_);
        gzbuffer_size_ = gzbuffer_init_size;
        gzbuffer_ = new char[gzbuffer_size_];
    } else {
        cerr << "Error: File '" << path_ << "' : expected '.vcf(.gz)' suffix.\n";
        throw exception();
    }
    read_header();
}

void
VcfParser::read_header()
{
    auto malformed = [this] () {
        cerr << "Error: Malformed header."
             << " At line " << line_number_ << " in file '" << path_ << "'.\n";
        throw exception();
    };

    // Meta header lines.
    while(getline() && strncmp(line(), "##", 2) == 0) {
        const char* equal = strchr(line(), '=');
        if(!equal) {
            // Skip header lines missing the '='.
            continue;
        }
        header_.add_meta(VcfMeta(string(line()+2, equal), string(equal+1)));
    }

    // Final header line.
    if(strncmp(line(), VcfHeader::std_fields.c_str(), VcfHeader::std_fields.length()) != 0)
        malformed();

    // Parse sample names, if any.
    if(line_len() > VcfHeader::std_fields.length()) {
        const char* p = line() + VcfHeader::std_fields.length();
        const char* end = line() + line_len();

        // Check that FORMAT is present.
        const char format[] = "\tformat";
        const size_t format_len = sizeof(format) - 1;
        if(strncmp(p, format, format_len) != 0)
            malformed();
        p += format_len;
        if (*p != '\t')
            malformed();

        do {
            ++p;
            const char* next = strchr(p, '\t');
            if (!next)
                next = end;
            header_.add_sample(string(p, next));
            p = next;
        } while (p != end);
    }
}

bool VcfParser::getline() {
    auto truncated = [this]() {
        cerr << "Error: While reading '" << path_ << "' (file may be truncated).\n";
        throw exception();
    };

    if (!is_gzipped_) {
        if (!std::getline(ifs_, ifsbuffer_))
            return false;
        if (ifsbuffer_.back() != '\n')
            truncated();
        // Remove the '\n'.
        ifsbuffer_.pop_back();
        if (ifsbuffer_.back() == '\r')
            // Remove the '\r'.
            ifsbuffer_.pop_back();

    } else {
        if (!gzgets(gzfile_, gzbuffer_, gzbuffer_size_))
            return false;
        gzline_len_ = strlen(gzbuffer_);

        // Check the contents of the buffer.
        while (gzbuffer_[gzline_len_ - 1] != '\n') {
            if (gzline_len_ < gzbuffer_size_ - 1)
                // The buffer was long enough but file didn't end with a newline.
                truncated();

            assert(gzline_len_ == gzbuffer_size_ - 1);

            // The buffer wasn't long enough.
            // Get a new, wider buffer & copy what we've already read.
            char* old_buffer = gzbuffer_;
            gzbuffer_size_ *= 2;
            gzbuffer_ = new char[gzbuffer_size_];
            memcpy(gzbuffer_, old_buffer, gzline_len_ + 1);
            delete[] old_buffer;

            // Continue reading the line.
            char* start = gzbuffer_ + gzline_len_;
            assert(*start == '\0');
            if (!gzgets(gzfile_, start, gzbuffer_size_ - gzline_len_))
                // EOF; this means that the last character read by the previous
                // call was the last of the file. And it wasn't a newline.
                truncated();
            gzline_len_ += strlen(start);
        }
        // Remove the '\n'.
        gzbuffer_[gzline_len_ - 1] = '\0';
        --gzline_len_;
        if (gzbuffer_[gzline_len_ - 1] == '\r') {
            // Remove the '\r'.
            gzbuffer_[gzline_len_ - 1] = '\0';
            --gzline_len_;
        }
    }

    ++line_number_;
    return true;
}

void VcfWriter::write_header() {
    for(const VcfMeta& m : header_.meta())
        file_ << "##" << m.key() << "=" << m.value() << "\n";

    file_ << VcfHeader::std_fields;
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
        file_ << r.format(0);
        for (size_t i=1; i<r.n_formats(); ++i)
            file_ << ":" << r.format(i);

        //SAMPLES
        assert(r.n_samples() != header_.samples().size());
        for (size_t i=0; i<r.n_samples(); ++i)
            file_ << "\t" << r.sample(i);
    }

    file_ << "\n";
}

