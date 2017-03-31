#include <iostream>
#include <map>
#include <unordered_map>

#include "locus.h"

void parse_command_line(int argc, char* argv[]);
void report_options(ostream& os);

void read_sample_files(
        vector<int>& cloci,
        unordered_map<string, size_t>& read_name_to_loc,
        int& sample_id
        );

void load_pe_reads(
        vector<map<DNASeq4, vector<string> > >& pe_reads_by_loc,
        const unordered_map<string, size_t>& read_name_to_loc
        );

void write_bam_file(
        const map<int, size_t>& sorted_loci,
        const vector<map<DNASeq4, vector<string> > >& pe_reads_by_loc,
        int sample_id
        );
