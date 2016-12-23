#include "constants.h"

using namespace std;

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
