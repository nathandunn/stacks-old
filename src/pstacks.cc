// -*-mode:c++; c-style:k&r; c-basic-offset:4;-*-
//
// Copyright 2010-2016, Julian Catchen <jcatchen@illinois.edu>
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
// pstacks -- search an existing set of stacks for polymorphisms
//

#include "log_utils.h"

#include "pstacks_base.h"
#include "pstacks.h"

using namespace std;
//
// Global variables to hold command-line options.
//
FileT  in_file_type;
string in_file;
string prefix_path;
int    sql_id        = -1;
int    min_stack_cov = 3;
double max_clipped   = 0.15;
int    min_mapping_qual = 10;
bool   keep_sec_alns = false;
int    num_threads   = 1;

//
// For use with the multinomial model to call fixed nucleotides.
//
modelt model_type         = snp;
double alpha              = 0.05;

int main (int argc, char* argv[]) {
    IF_NDEBUG_TRY

    parse_command_line(argc, argv);

    report_options(cerr);
    cerr << std::fixed << std::setprecision(2);

    //
    // Set limits to call het or homozygote according to chi-square distribution with one
    // degree of freedom:
    //   http://en.wikipedia.org/wiki/Chi-squared_distribution#Table_of_.CF.872_value_vs_p-value
    //
    if (alpha == 0.1) {
        heterozygote_limit = -2.71;
        homozygote_limit   =  2.71;
    } else if (alpha == 0.05) {
        heterozygote_limit = -3.84;
        homozygote_limit   =  3.84;
    } else if (alpha == 0.01) {
        heterozygote_limit = -6.64;
        homozygote_limit   =  6.64;
    } else if (alpha == 0.001) {
        heterozygote_limit = -10.83;
        homozygote_limit   =  10.83;
    }

    //
    // Set the number of OpenMP parallel threads to execute.
    //
    #ifdef _OPENMP
    omp_set_num_threads(num_threads);
    #endif

    map<int, PStack *> unique; // Unique {sequence, alignment position} combinations.
    HashMap*           radtags = new HashMap();
    load_radtags(in_file, *radtags);
    reduce_radtags(*radtags, unique);
    for (auto& stack : *radtags)
        for (Seq* read : stack.second)
            delete read;
    delete radtags;

    map<int, MergedStack *> merged;
    populate_merged_tags(unique, merged);

    delete_low_cov_loci(merged, unique);

    // Call the consensus sequence again, now that remainder tags have been merged.
    call_consensus(merged, unique, true);

    write_results(merged, unique, true, false);

    cerr << "pstacks is done.\n";

    return 0;
    IF_NDEBUG_CATCH_ALL_EXCEPTIONS
}

void populate_merged_tags(map<int, PStack *> &unique, map<int, MergedStack *> &merged) {
    map<int, PStack *>::iterator i;
    map<int, MergedStack *>::iterator it_new, it_old;
    map<string, set<int> > locations;
    map<string, set<int> >::iterator k;
    set<int>::iterator s;
    char         id[id_len];
    PStack      *u;
    MergedStack *m;
    int global_id = 1;

    //
    // Create a map of each unique Stack that has been aligned to the same genomic location.
    //
    for (i = unique.begin(); i != unique.end(); i++) {
        snprintf(id, id_len - 1, "%s|%d|%s",
                 i->second->loc.chr(),
                 i->second->loc.bp,
                 i->second->loc.strand == strand_plus ? "+" : "-");
        locations[id].insert(i->second->id);
    }

    it_old = merged.begin();

    for (k = locations.begin(); k != locations.end(); k++) {
        m = new MergedStack;
        m->id = global_id;

        //
        // Record the consensus and physical location for this stack.
        //
        s = k->second.begin();
        m->add_consensus(unique[*s]->seq);
        m->loc.set(unique[*s]->loc.chr(), unique[*s]->loc.bp, unique[*s]->loc.strand);

        //
        // Record the individual stacks that were aligned together.
        //
        for (; s != k->second.end(); s++) {
            u = unique[*s];
            m->count += u->count;
            m->utags.push_back(u->id);
        }

        //
        // Insert the new MergedStack giving a hint as to which position
        // to insert it at.
        //
        it_new = merged.insert(it_old, pair<int, MergedStack *>(global_id, m));
        it_old = it_new;
        global_id++;
    }

    double mean;
    double stdev;
    double max;
    calc_coverage_distribution(unique, merged, mean, stdev, max);

    cerr << "Created " << merged.size() << " loci; mean coverage is " << mean << " (stdev: " << stdev << ", max: " << size_t(max) << ").\n";
}

void delete_low_cov_loci(map<int, MergedStack *>& merged, const map<int, PStack*>& unique) {

    size_t n_deleted = 0;
    size_t n_reads = 0;
    size_t n_rm_reads = 0;
    vector<int> to_erase;

    for (auto& mtag : merged) {
        int depth = 0;
        for (int utag_id : mtag.second->utags)
            depth += unique.at(utag_id)->count;

        n_reads += depth;
        if (depth < min_stack_cov) {
            delete mtag.second;
            to_erase.push_back(mtag.first);

            n_deleted++;
            n_rm_reads += depth;
        }
    }

    for (int id : to_erase)
        merged.erase(id);

    double mean;
    double stdev;
    double max;
    calc_coverage_distribution(unique, merged, mean, stdev, max);

    cerr << "Discarded " << n_deleted << " low coverage loci comprising " << n_rm_reads
         << " (" << as_percentage((double) n_rm_reads/ n_reads) << ") reads.\n"
         << "Kept " << merged.size() << " loci; mean coverage is "
         << mean << " (stdev: " << stdev << ", max: " << size_t(max) << ").\n";
}

//
// This function assumes that there may be identical reads, mapped to multiple
// places in the genome. In this case, reads are broken down by read ID
// and split into different Stack objects.
//
int reduce_radtags(HashMap &radtags, map<int, PStack *> &unique) {
    HashMap::iterator it;
    vector<Seq *>::iterator sit;

    PStack *u;
    int    global_id = 1;

    for (it = radtags.begin(); it != radtags.end(); it++) {
        //
        // Make sure there aren't any reads of identical sequence that have been mapped to
        // different genomic locations.
        //
        map<string, int> locations;
        map<string, int>::iterator lit;
        for (sit = (*it).second.begin(); sit != (*it).second.end(); sit++)
            locations[(*sit)->loc_str]++;

        for (lit = locations.begin(); lit != locations.end(); lit++) {
            //
            // Populate a Stack object for this unique radtag.
            //
            u        = new PStack;
            u->id    = global_id;
            u->count = lit->second;
            u->add_seq(&it->first);

            //
            // Record the physical location of this stack.
            //
            for (sit = (*it).second.begin(); sit != (*it).second.end(); sit++) {
                if (strcmp((*sit)->loc_str, lit->first.c_str()) == 0) {
                    u->add_id((*sit)->id);
                    u->loc.set((*sit)->loc.chr(), (*sit)->loc.bp, (*sit)->loc.strand);
                }
            }

            unique[u->id] = u;
            global_id++;
        }
    }

    return 0;
}

//
// We expect tags to have already been aligned to a reference genome. Therefore, the tags
// are identified by their chromosome and basepair location.
//
void load_radtags(string in_file, HashMap &radtags) {
    Input *fh = NULL;
    Seq c;

    if (in_file_type == FileT::bowtie)
        fh = new Bowtie(in_file.c_str());
    else if (in_file_type == FileT::sam)
        fh = new Sam(in_file.c_str());
    else if (in_file_type == FileT::bam)
        fh = new Bam(in_file.c_str());
    else if (in_file_type == FileT::tsv)
        fh = new Tsv(in_file.c_str());

    cerr << "Reading alignments...\n";

    int primary_kept   = 0;
    int primary_qual   = 0;
    int primary_clipped = 0;
    int secondary_kept = 0;
    int secondary_disc = 0;
    int supplementary  = 0;
    int unmapped       = 0;
    while ((fh->next_seq(c)) != 0) {

        switch (c.aln_type) {
        case AlnT::null:
            unmapped++;
            continue;
            break;
        case AlnT::primary:
            if (c.map_qual < min_mapping_qual) {
                primary_qual++;
                continue;
            } else if (c.pct_clipped > max_clipped) {
                primary_clipped++;
                continue;
            } else {
                primary_kept++;
            }
            break;
        case AlnT::secondary:
            if (keep_sec_alns && c.pct_clipped <= max_clipped) {
                secondary_kept++;
            } else {
                secondary_disc++;
                continue;
            }
            break;
        case AlnT::supplementary:
            supplementary++;
            continue;
            break;
        default:
            DOES_NOT_HAPPEN;
            break;
        }

        HashMap::iterator element = radtags.insert({DNANSeq(strlen(c.seq), c.seq), vector<Seq*>()}).first;
        element->second.push_back(new Seq(c));
        Seq& the_seq = *element->second.back();
        if (the_seq.seq != NULL) {
            delete[] the_seq.seq;
            the_seq.seq = NULL;
        }
        if (the_seq.qual != NULL) {
            delete[] the_seq.qual;
            the_seq.qual = NULL;
        }
    }

    int n_primary = primary_kept+primary_qual+primary_clipped;
    cerr << "Done reading alignment records:\n"
         << "  Kept " << primary_kept << " primary alignments\n"
         << "  Skipped " << primary_qual << " (" << as_percentage((double) primary_qual / n_primary)
         << ") primary alignments with insufficient mapping qualities\n"
         << "  Skipped " << primary_clipped << " (" << as_percentage((double) primary_clipped / n_primary)
         << ") excessively soft-clipped primary alignments\n"
         << "  Skipped " << secondary_disc << " secondary alignments\n"
         << "  Skipped " << supplementary << " supplementary alignments\n"
         << "  Skipped " << unmapped << " unmapped reads\n";

    if(keep_sec_alns)
        cerr << "  Kept " << secondary_kept << " secondary alignments\n";

    if (radtags.empty()) {
        cerr << "Error: No input.\n";
        throw exception();
    }

    cerr << "Collapsed reads into " << radtags.size() << " stacks.\n";

    delete fh;
}

int dump_stacks(map<int, PStack *> &u) {
    map<int, PStack *>::iterator it;
    vector<char *>::iterator fit;
    vector<pair<int, int> >::iterator pit;
    vector<int>::iterator mit;

    for (it = u.begin(); it != u.end(); it++) {

        cerr << "Stack ID: " << (*it).second->id << "\n"
             << "  Seq:    " << (*it).second->seq->seq() << "\n"
             << "  IDs:    ";

        for (fit = (*it).second->map.begin(); fit != (*it).second->map.end(); fit++)
            cerr << *fit << " ";

        cerr << "\n\n";
    }

    return 0;
}

int dump_merged_stacks(map<int, MergedStack *> &m) {
    map<int, MergedStack *>::iterator it;
    vector<pair<int, int> >::iterator pit;
    vector<int>::iterator fit;

    for (it = m.begin(); it != m.end(); it++) {

        cerr << "MergedStack ID: " << it->second->id << "\n"
             << "  Consensus:  ";
        if (it->second->con != NULL)
            cerr << it->second->con << "\n";
        else
            cerr << "\n";
        cerr << "  IDs:        ";

        for (fit = it->second->utags.begin(); fit != it->second->utags.end(); fit++)
            cerr << (*fit) << " ";

        cerr << "\n"
             << "  Distances: ";

        for (pit = it->second->dist.begin(); pit != it->second->dist.end(); pit++)
            cerr << (*pit).first << ": " << (*pit).second << ", ";

        cerr << "\n\n";
    }

    return 0;
}

int parse_command_line(int argc, char* argv[]) {
    string out_path;
    int c;

    while (1) {
        static struct option long_options[] = {
            {"help",         no_argument,       NULL, 'h'},
            {"version",      no_argument,       NULL, 1000},
            {"infile_type",  required_argument, NULL, 't'},
            {"file",         required_argument, NULL, 'f'},
            {"outpath",      required_argument, NULL, 'o'},
            {"id",           required_argument, NULL, 'i'},
            {"min_cov",      required_argument, NULL, 'm'},
            {"max_clipped",  required_argument, NULL, 1001},
            {"min_mapq",     required_argument, NULL, 1002},
            {"keep_sec_aln", required_argument, NULL, 'k'},
            {"num_threads",  required_argument, NULL, 'p'},
            {"bc_err_freq",  required_argument, NULL, 'e'},
            {"model_type",   required_argument, NULL, 'T'},
            {"bound_low",    required_argument, NULL, 'L'},
            {"bound_high",   required_argument, NULL, 'U'},
            {"alpha",        required_argument, NULL, 'A'},
            {0, 0, 0, 0}
        };

        // getopt_long stores the option index here.
        int option_index = 0;

        c = getopt_long(argc, argv, "hkvOT:a:A:L:U:f:o:i:e:p:m:s:f:t:", long_options, &option_index);

        // Detect the end of the options.
        if (c == -1)
            break;

        switch (c) {
        case 'h':
            help();
            break;
        case 't':
            if (strcmp(optarg, "bowtie") == 0)
                in_file_type = FileT::bowtie;
            else if (strcmp(optarg, "sam") == 0)
                in_file_type = FileT::sam;
            else if (strcmp(optarg, "bam") == 0)
                in_file_type = FileT::bam;
            else if (strcmp(optarg, "tsv") == 0)
                in_file_type = FileT::tsv;
            else
                in_file_type = FileT::unknown;
            break;
        case 'f':
            in_file = optarg;
            break;
        case 'o':
            out_path = optarg;
            break;
        case 'i':
            sql_id = is_integer(optarg);
            break;
        case 'm':
            min_stack_cov = atoi(optarg);
            break;
        case 1001: //max_clipped
            max_clipped = is_double(optarg);
            if (max_clipped > 1)
                max_clipped = max_clipped / 100;

            if (max_clipped < 0 || max_clipped > 1.0) {
                cerr << "Unable to parse the maximum clipped proportion.\n";
                help();
            }
            break;
        case 1002:
            min_mapping_qual = is_integer(optarg);
            if (min_mapping_qual < 0) {
                cerr << "Unable to parse the minimum mapping quality.\n";
                help();
            }
            break;
        case 'k':
            keep_sec_alns = true;
            break;
        case 'e':
            barcode_err_freq = atof(optarg);
            break;
        case 'T':
            if (strcmp(optarg, "snp") == 0) {
                model_type = snp;
            } else if (strcmp(optarg, "fixed") == 0) {
                model_type = ::fixed;
            } else if (strcmp(optarg, "bounded") == 0) {
                model_type = bounded;
            } else {
                cerr << "Unknown model type specified '" << optarg << "'\n";
                help();
            }
            break;
        case 'L':
            bound_low  = atof(optarg);
            break;
        case 'U':
            bound_high = atof(optarg);
            break;
        case 'A':
            alpha = atof(optarg);
            break;
        case 'p':
            num_threads = atoi(optarg);
            break;
        case 1000:
            version();
            break;
        case '?':
            // getopt_long already printed an error message.
            help();
            break;

        default:
            cerr << "Unknown command line option '" << (char) c << "'\n";
            help();
            exit(1);
        }
    }

    if (optind < argc) {
        cerr << "Error: Failed to parse command line: '" << argv[optind] << "' is seen as a positional argument. Expected no positional arguments.\n";
        help();
    }

    if (alpha != 0.1 && alpha != 0.05 && alpha != 0.01 && alpha != 0.001) {
        cerr << "SNP model alpha significance level must be either 0.1, 0.05, 0.01, or 0.001.\n";
        help();
    }

    if (bound_low != 0 && (bound_low < 0 || bound_low >= 1.0)) {
        cerr << "SNP model lower bound must be between 0.0 and 1.0.\n";
        help();
    }

    if (bound_high != 1 && (bound_high <= 0 || bound_high > 1.0)) {
        cerr << "SNP model upper bound must be between 0.0 and 1.0.\n";
        help();
    }

    if (bound_low > 0 || bound_high < 1.0) {
        model_type = bounded;
    }

    if (model_type == ::fixed && barcode_err_freq == 0) {
        cerr << "You must specify the barcode error frequency.\n";
        help();
    }
    if (in_file.empty()) {
        cerr << "You must specify an input file.\n";
        help();
    }

    if (in_file_type == FileT::unknown) {
        in_file_type = guess_file_type(in_file);
        if (in_file_type == FileT::unknown) {
            cerr << "Unable to recongnize the extention of file '" << in_file << "'.\n";
            help();
        }
    }

    if (sql_id < 0) {
        cerr << "A sample ID must be provided.\n";
        help();
    }

    if (out_path.empty())
        out_path = ".";

    if (out_path.at(out_path.length() - 1) != '/')
        out_path += "/";

    string s = in_file.rfind('/') == string::npos ? in_file : in_file.substr(in_file.rfind('/') + 1);
    s = s.substr(0, s.rfind('.'));
    prefix_path = out_path + s;

    return 0;
}

void version() {
    cerr << "pstacks " << VERSION << "\n\n";

    exit(1);
}

void help() {
    cerr << "pstacks " << VERSION << "\n"
              << "pstacks -f file_path -i id -o path [-m min_cov] [-p num_threads]" << "\n"
              << "  f: input file path.\n"
              << "  i: a unique integer ID for this sample.\n"
              << "  o: output directory.\n"
              << "  m: minimum depth of coverage to report a stack (default 3).\n"
              << "  p: enable parallel execution with num_threads threads.\n"
              << "  t: input file Type. Supported types: bam, sam, bowtie (default: guess).\n"
              << "  --max_clipped <float>: alignments with more than this fraction of soft-clipped bases are discarded (default 15%).\n"
              << "  --min_mapq <int>: minimum required quality (default 10).\n"
              << "  --keep_sec_alns: keep secondary alignments (default: false, only keep primary alignments).\n"
              << "\n"
              << "  Model options:\n"
              << "    --model_type <type>: either 'snp' (default), 'bounded', or 'fixed'\n"
              << "       For the SNP or Bounded SNP model:\n"
              << "       --alpha <num>: chi square significance level required to call a heterozygote or homozygote, either 0.1, 0.05 (default), 0.01, or 0.001.\n"
              << "       For the Bounded SNP model:\n"
              << "       --bound_low <num>: lower bound for epsilon, the error rate, between 0 and 1.0 (default 0).\n"
              << "       --bound_high <num>: upper bound for epsilon, the error rate, between 0 and 1.0 (default 1).\n"
              << "       For the Fixed model:\n"
              << "       --bc_err_freq <num>: specify the barcode error frequency, between 0 and 1.0.\n";

    exit(1);
}

void report_options(ostream& os) {
    os << "pstacks parameters selected:\n"
       << "  Alignments file: " << in_file << "\n"
       << "  Output prefix: " << prefix_path << "\n"
       << "  Sample ID: " << sql_id << "\n"
       << "  Min locus depth: " << min_stack_cov << "\n"
       << "  Max clipped proportion: " << max_clipped << "\n"
       << "  Min mapping quality: " << min_mapping_qual << "\n";

    // Model.
    if (model_type == snp) {
        os << "  Model: snp\n"
           << "    Model alpha: " << alpha << "\n";
    } else if (model_type == bounded) {
        os << "  Model: snp\n"
           << "    Model alpha: " << alpha << "\n"
           << "    Model lower bound: " << bound_low << "\n"
           << "    Model higher bound: " << bound_high << "\n";
    } else if (model_type == ::fixed) {
        os << "  Model: fixed\n"
           << "    Model barcode err. prob.: " << barcode_err_freq << "\n";
    }
}
