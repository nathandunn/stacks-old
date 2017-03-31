#include <iostream>
#include <map>
#include <unordered_map>

#include "locus.h"

void parse_command_line(int argc, char* argv[]);
void report_options(ostream& os);

void read_sample_files(map<int, Locus*>& sloci, unordered_map<int, int>& sloc_to_cloc, int& sample_id);
void write_bam_file(const map<int, Locus*>& sorted_loci, int sample_id);
