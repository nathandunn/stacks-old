#include <iostream>
#include <map>
#include <unordered_map>

#include "locus.h"

void parse_command_line(int argc, char* argv[]);
void report_options(std::ostream& os);

void read_sample_files(
        std::vector<int>& cloci,
        std::unordered_map<string, size_t>& read_name_to_loc,
        int& sample_id
        );

void load_pe_reads(
        std::vector<std::map<DNASeq4, std::vector<std::string> > >& pe_reads_by_loc,
        const std::unordered_map<std::string, size_t>& read_name_to_loc
        );

void write_bam_file(
        const std::map<int, size_t>& sorted_loci,
        const std::vector<std::map<DNASeq4, std::vector<std::string> > >& pe_reads_by_loc,
        int sample_id
        );
