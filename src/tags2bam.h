#include <iostream>
#include <map>
#include <unordered_map>

#include "locus.h"

void parse_command_line(int argc, char* argv[]);
void report_options(std::ostream& os);

void read_sample_files(std::map<int, Locus*>& sloci, std::unordered_map<int, int>& sloc_to_cloc);
void write_bam_file(const std::map<int, Locus*>& sloci, const std::unordered_map<int, int>& sloc_to_cloc);
