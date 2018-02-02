// -*-mode:c++; c-style:k&r; c-basic-offset:4;-*-
//
// Copyright 2010, Julian Catchen <jcatchen@uoregon.edu>
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

//
// exstacks -- an Example Stacks program. Read in a set of reference-aligned stacks
// that have had SNPs identified by the pstacks program, sort the loci by genomic
// position and print them out.
//
// Julian Catchen
// jcatchen@uoregon.edu
// University of Oregon
//

#include "exstacks.h"

// Global variables to hold command-line options.
queue<string> samples;
string        out_path;
int           num_threads = 1;

int main (int argc, char* argv[]) {

    parse_command_line(argc, argv);

    //
    // Set the number of OpenMP parallel threads to execute.
    //
    //omp_set_num_threads(num_threads);

    map<int, Locus *> sample;

    string s = samples.front();
    samples.pop();

    bool tmp;
    if (!load_loci(s, sample, 0, false, tmp)) {
        cerr << "Failed to load sample " << s.c_str() << "\n";
        return 1;
    }

    // dump_loci(sample);

    //
    // Sort the Loci by genomic location.
    //
    map<string, vector<Locus *> > sorted;
    set<string> chrs;
    map<int, Locus *>::iterator i;
    set<string>::iterator j;

    //
    // First, bin loci according to chromosome and record a list of all chromosomes.
    //
    for (i = sample.begin(); i != sample.end(); i++) {
        sorted[i->second->loc.chr()].push_back(i->second);
        chrs.insert(i->second->loc.chr());
    }

    //
    // Second, order the loci on each chromosome
    //
    for (j = chrs.begin(); j != chrs.end(); j++)
        sort(sorted[*j].begin(), sorted[*j].end(), compare_loci);

    //
    // Now print things out in order.
    //
    vector<Locus *>::iterator k;

    for (j = chrs.begin(); j != chrs.end(); j++) {
        cerr << "Examining chromosome " << *j << "\n";

        for (k = sorted[*j].begin(); k != sorted[*j].end(); k++) {
            write_simple_output(*k);
        }
    }

    return 0;
}

bool compare_loci(Locus *a, Locus *b) {
    return (a->loc.bp < b->loc.bp);
}

int write_simple_output(Locus *tag) {
    vector<SNP *>::iterator           snp_it;
    map<string, int>::iterator        all_it;
    vector<pair<int, int> >::iterator src_it;
    int i;

    cout <<
        "  Locus: " <<
        tag->id    << ", chr: " << tag->loc.chr() << " " << tag->loc.bp << "bp; " <<
        tag->con    << "\n";

    //
    // Output the SNPs associated with the catalog tag
    //
    i = 1;
    for (snp_it = tag->snps.begin(); snp_it != tag->snps.end(); snp_it++) {
        cout << "    SNP " << i << ": " <<
            "col: "    << (*snp_it)->col    << " " <<
            "lratio: " << (*snp_it)->lratio << " " <<
            "rank_1: " << (*snp_it)->rank_1 << " " <<
            "rank_2: " << (*snp_it)->rank_2 << "\n";
        i++;
    }

    //
    // Output the alleles associated with the two matched tags
    //
    i = 1;
    for (all_it = tag->alleles.begin(); all_it != tag->alleles.end(); all_it++) {
        cout << "    Allele " << i << ": " <<
            all_it->first << "\n";
        i++;
    }

    return 0;
}

int parse_command_line(int argc, char* argv[]) {
    int c;
    string sstr;

    while (1) {
        static struct option long_options[] = {
            {"help",         no_argument,       NULL, 'h'},
            {"version",      no_argument,       NULL, 'v'},
            {"sample",       required_argument, NULL, 's'},
            {"outpath",      required_argument, NULL, 'o'},
            {"num_threads",  required_argument, NULL, 'p'},
            {0, 0, 0, 0}
        };

        // getopt_long stores the option index here.
        int option_index = 0;

        c = getopt_long(argc, argv, "hvo:s:p:", long_options, &option_index);

        // Detect the end of the options.
        if (c == -1)
            break;

        switch (c) {
        case 'h':
            help();
            break;
        case 's':
            sstr = optarg;
            samples.push(sstr);
            break;
        case 'o':
            out_path = optarg;
            break;
        case 'v':
            version();
            break;
        case 'p':
            num_threads = atoi(optarg);
            break;
        case '?':
            // getopt_long already printed an error message.
            help();
            break;
        default:
            help();
            exit(1);
        }
    }

    if (samples.size() == 0) {
        cerr << "You must specify at least one sample file.\n";
        help();
    }

    if (out_path.length() == 0)
        out_path = ".";

    if (out_path.at(out_path.length() - 1) != '/')
        out_path += "/";

    return 0;
}

void version() {
    cerr << "exstacks " << stacks_version << "\n\n";

    exit(1);
}

void help() {
    cerr << "exstacks " << stacks_version << "\n"
              << "exstacks -s sample_file [-o path] [-p num_threads] [-h]" << "\n"
              << "  p: enable parallel execution with num_threads threads.\n"
              << "  s: TSV file from which to load radtags." << "\n"
              << "  o: output path to write results." << "\n"
              << "  m: include tags in the catalog that match to more than one entry." << "\n"
              << "  h: display this help messsage." << "\n\n";

    exit(1);
}
