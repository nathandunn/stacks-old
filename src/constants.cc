#include <regex>
#include <cassert>

#include "constants.h"

using namespace std;

const
map<string, FileT> known_extensions = {
    {".fa", FileT::fasta},
    {".fasta", FileT::fasta},
    {".fa.gz", FileT::gzfasta},
    {".fasta.gz", FileT::gzfasta},
    {".fq", FileT::fastq},
    {".fastq", FileT::fastq},
    {".fq.gz", FileT::gzfastq},
    {".fastq.gz", FileT::gzfastq},
    {".sam", FileT::sam},
    {".bam", FileT::bam}
};

string remove_suffix(FileT type, const string& orig) {
    string file (orig);

    int pos = file.find_last_of(".");

    if ((type == FileT::gzfastq || type == FileT::gzfasta) && file.substr(pos) == ".gz")
        file = file.substr(0, pos);

    pos = file.find_last_of(".");

    if (type == FileT::gzfastq || type == FileT::fastq) {

        if (file.substr(pos) == ".fastq" || file.substr(pos) == ".fq")
            file = file.substr(0, pos);

    } else if (type == FileT::gzfasta || type == FileT::fasta) {

        if (file.substr(pos) == ".fasta" || file.substr(pos) == ".fa")
            file = file.substr(0, pos);
    }

    return file;
}

regex init_file_ext_regex () {
    string s = "(";

    auto i = known_extensions.begin();
    assert(!known_extensions.empty());
    string ext = i->first;
    escape_char('.', ext);
    s += ext;
    ++i;
    while(i != known_extensions.end()) {
        ext = i->first;
        escape_char('.', ext);
        s += "|" + ext;
        ++i;
    }

    s += ")$";
    return regex(s);
}

FileT guess_file_type (const string& path) {

    static const regex reg = init_file_ext_regex();

    smatch m;
    regex_search(path, m, reg);

    if (m.empty())
        return FileT::unknown;
    else
        return known_extensions.at(m.str());
}

void escape_char(char c, string& s) {
    vector<size_t> dots;
    size_t i = -1;
    while ((i = s.find(c, i+1)) != string::npos)
        dots.push_back(i);

    for(auto j=dots.rbegin(); j!=dots.rend(); ++j)
        s.insert(*j, 1, '\\');
}
