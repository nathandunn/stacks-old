// -*-mode:c++; c-style:k&r; c-basic-offset:4;-*-
//
// Copyright 2013-2015, Julian Catchen <jcatchen@illinois.edu>
//
// This file is part of Stacks.
//
// Stacks is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// Stacks is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with Stacks.  If not, see <http://www.gnu.org/licenses/>.
//

#include <ctime>
#include <fstream>
#include <sstream>
#include <iomanip>

#include "constants.h"
#include "log_utils.h"

using namespace std;

int
init_log(ostream &fh, int argc, char **argv)
{
    //
    // Obtain the current date.
    //
    time_t     rawtime;
    struct tm *timeinfo;
    char       date[32];
    time(&rawtime);
    timeinfo = localtime(&rawtime);
    strftime(date, 32, "%F %T", timeinfo);

    //
    // Write the command line that was executed.
    //
    for (int i = 0; i < argc; i++) {
        fh << argv[i];
        if (i < argc - 1) fh << " ";
    }
    fh << "\n" << argv[0] << " version " << VERSION << " executed " << date << "\n\n";

    return 0;
}

string as_percentage(double d) {
    stringstream ss;
    ss << std::fixed << setprecision(1);
    ss << d*100 << "%";
    return ss.str();
}

string to_string(const FileT& ft) {
    if (ft == FileT::unknown)
        return "unknown";
    if (ft == FileT::sql)
        return "sql";
    if (ft == FileT::gzsql)
        return "gzsql";
    if (ft == FileT::fasta)
        return "fasta";
    if (ft == FileT::gzfasta)
        return "gzfasta";
    if (ft == FileT::fastq)
        return "fastq";
    if (ft == FileT::gzfastq)
        return "gzfastq";
    if (ft == FileT::bowtie)
        return "bowtie";
    if (ft == FileT::sam)
        return "sam";
    if (ft == FileT::bam)
        return "bam";
    if (ft == FileT::tsv)
        return "tsv";
    if (ft == FileT::bustard)
        return "bustard";
    if (ft == FileT::phase)
        return "phase";
    if (ft == FileT::fastphase)
        return "fastphase";
    if (ft == FileT::beagle)
        return "beagle";
    assert(false);
    return "?!";
}

LogAlterator::LogAlterator(const std::string& log_path, bool quiet)
        : l(log_path)
        , o(std::cout.rdbuf())
        , e(std::cerr.rdbuf())
        , lo_buf(std::cout.rdbuf(), l.rdbuf())
        , le_buf(std::cerr.rdbuf(), l.rdbuf())
        {
    if (quiet) {
        // Use the fstream buffer only, and not at all stdout and stderr.
        std::cout.rdbuf(l.rdbuf());
        std::cerr.rdbuf(l.rdbuf());
    } else {
        std::cout.rdbuf(&lo_buf);
        std::cerr.rdbuf(&le_buf);
    }
}
