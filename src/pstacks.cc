// -*-mode:c++; c-style:k&r; c-basic-offset:4;-*-
//
// Copyright 2010-2015, Julian Catchen <jcatchen@illinois.edu>
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

#include "pstacks.h"
#include "pstacks_base.h"

//
// Global variables to hold command-line options.
//
FileT  in_file_type;
string in_file;
FileT  out_file_type;
string out_path;
string prefix_path;
int    sql_id        = 0;
int    min_stack_cov = 3;
int    num_threads   = 1;

//
// For use with the multinomial model to call fixed nucleotides.
//
modelt model_type         = snp;
double alpha              = 0.05;
double bound_low          = 0.0;
double bound_high         = 1.0;
double p_freq             = 0.5;
double barcode_err_freq   = 0.0;
double heterozygote_limit = -3.84;
double homozygote_limit   =  3.84;

int main (int argc, char* argv[]) {

    parse_command_line(argc, argv);

    cerr << "Min depth of coverage to report a stack: " << min_stack_cov << "\n"
         << "Model type: ";
    switch (model_type) {
    case snp:
        cerr << "SNP\n";
        break;
    case fixed:
        cerr << "Fixed\n";
        break;
    case bounded:
        cerr << "Bounded; lower epsilon bound: " << bound_low << "; upper bound: " << bound_high << "\n";
        break;
    }
    cerr << "Alpha significance level for model: " << alpha << "\n";

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

    HashMap*           radtags = new HashMap();
    map<int, PStack *> unique;

    load_radtags(in_file, *radtags);

    reduce_radtags(*radtags, unique);

    for (auto& stack : *radtags)
        for (Seq* read : stack.second)
            delete read;
    delete radtags;

    //dump_stacks(unique);

    map<int, MergedStack *> merged;

    populate_merged_tags(unique, merged);

    //dump_merged_stacks(merged);

    size_t merged_size_old = merged.size();
    prune_low_coverage_loci(merged, unique);
    cerr << "Excluded " << merged_size_old - merged.size()
         << " loci due to insuffient depth of coverage; now working with "
         << merged.size() << " loci.\n";

    // Call the consensus sequence again, now that remainder tags have been merged.
    cerr << "Identifying polymorphic sites and calling consensus sequences...";
    call_consensus(merged, unique, true);
    cerr << "done.\n";

    count_raw_reads(unique, merged);

    calc_coverage_distribution(unique, merged);

    cerr << "Writing loci, SNPs, alleles to '" << prefix_path << ".*'...\n";
    const bool gzip = in_file_type == FileT::bam;
    write_results(merged, unique, gzip, false);

    return 0;
}

double calc_coverage_distribution(map<int, PStack *> &unique, map<int, MergedStack *> &merged) {
    map<int, MergedStack *>::iterator it;
    vector<int>::iterator             k;
    PStack *tag;
    double  depth = 0.0;
    double  total = 0.0;
    double  sum   = 0.0;
    double  mean  = 0.0;
    double  max   = 0.0;
    double  stdev = 0.0;

    for (it = merged.begin(); it != merged.end(); it++) {
        depth = 0.0;
        for (k = it->second->utags.begin(); k != it->second->utags.end(); k++) {
            tag    = unique[*k];
            depth += tag->count;
        }

        if (depth > max)
            max = depth;

        sum += depth;
        total++;
    }

    mean = sum / total;

    //
    // Calculate the standard deviation
    //
    for (it = merged.begin(); it != merged.end(); it++) {
        depth = 0.0;
        for (k = it->second->utags.begin(); k != it->second->utags.end(); k++) {
            tag    = unique[*k];
            depth += tag->count;
        }

        sum += pow((depth - mean), 2);
    }

    stdev = sqrt(sum / (total - 1));

    cerr << "  Mean coverage depth is " << mean << "; Std Dev: " << stdev << "; Max: " << max << "\n";

    return mean;
}

int count_raw_reads(map<int, PStack *> &unique, map<int, MergedStack *> &merged) {
    map<int, MergedStack *>::iterator it;
    vector<int>::iterator k;
    PStack *tag;
    long int m = 0;

    for (it = merged.begin(); it != merged.end(); it++) {
        for (k = it->second->utags.begin(); k != it->second->utags.end(); k++) {
            tag  = unique[*k];
            m   += tag->count;
        }
        m += it->second->remtags.size();
    }

    cerr << "  Number of utilized reads " << m << "\n";

    return 0;
}

void prune_low_coverage_loci(map<int, MergedStack *>& merged, const map<int, PStack *>& unique) {
    for(auto m_it = merged.begin(); m_it != merged.end();) {
        uint tot_depth = 0;
        for (int u : m_it->second->utags)
            tot_depth += unique.at(u)->count;

        if (tot_depth < min_stack_cov) {
            auto m_it_copy = m_it;
            ++m_it;
            delete m_it_copy->second;
            merged.erase(m_it_copy);
        } else {
            ++m_it;
        }
    }
}

int populate_merged_tags(map<int, PStack *> &unique, map<int, MergedStack *> &merged) {
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
                 i->second->loc.chr,
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
        m->loc.set(unique[*s]->loc.chr, unique[*s]->loc.bp, unique[*s]->loc.strand);

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

    cerr << "  Merged " << unique.size() << " unique Stacks into " << merged.size() << " loci.\n";

    return 0;
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
                    u->loc.set((*sit)->loc.chr, (*sit)->loc.bp, (*sit)->loc.strand);
                }
            }

            unique[u->id] = u;
            global_id++;
        }
    }

    cerr << "  " << radtags.size() << " unique stacks were aligned to " << unique.size() << " genomic locations.\n";

    return 0;
}

//
// We expect tags to have already been aligned to a reference genome. Therefore, the tags
// are identified by their chromosome and basepair location.
//
int load_radtags(string in_file, HashMap &radtags) {
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

    cerr << "Parsing " << in_file.c_str() << "\n";

    int i = 0;
    cerr << "Loading aligned sequences...";
    while ((fh->next_seq(c)) != 0) {
        if (i % 1000000 == 0 && i>0)
            cerr << i/1000000 << "M...";

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

        i++;
    }
    cerr << "done\n";

    if (i == 0) {
        cerr << "Error: Unable to load data from '" << in_file.c_str() << "'.\n";
        exit(-1);
    }

    cerr << "Loaded " << i << " sequence reads; " <<
        "identified " << radtags.size() << " unique stacks from those reads.\n";

    //
    // Close the file and delete the Input object.
    //
    delete fh;

    return 0;
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
            {"version",      no_argument,       NULL, 'v'},
            {"infile_type",  required_argument, NULL, 't'},
            {"file",         required_argument, NULL, 'f'},
            {"outpath",      required_argument, NULL, 'o'},
            {"id",           required_argument, NULL, 'i'},
            {"min_cov",      required_argument, NULL, 'm'},
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

        c = getopt_long(argc, argv, "hvOT:A:L:U:f:o:i:e:p:m:s:f:t:y:", long_options, &option_index);

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
            if (sql_id < 0) {
                cerr << "SQL ID (-i) must be an integer, e.g. 1, 2, 3\n";
                help();
            }
            break;
        case 'm':
            min_stack_cov = atoi(optarg);
            break;
        case 'e':
            barcode_err_freq = atof(optarg);
            break;
        case 'T':
            if (strcmp(optarg, "snp") == 0) {
                model_type = snp;
            } else if (strcmp(optarg, "fixed") == 0) {
                model_type = fixed;
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
        case 'v':
            version();
            break;
        case '?':
            // getopt_long already printed an error message.
            help();
            break;

        default:
            cerr << "Unknown command line option '" << (char) c << "'\n";
            help();
            abort();
        }
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

    if (in_file.length() == 0) {
        cerr << "You must specify an input file.\n";
        help();
    }

    if (in_file_type == FileT::unknown) {
        if (in_file.length() >= 4 && in_file.substr(in_file.length() - 4) == ".bam") {
            in_file_type = FileT::bam;
        } else if (in_file.length() >= 4 && in_file.substr(in_file.length() - 4) == ".sam") {
            in_file_type = FileT::sam;
        } else {
            cerr << "Unable to guess the type of the input file, please specify it.\n";
            help();
        }
    }

    if (model_type == fixed && barcode_err_freq == 0) {
        cerr << "You must specify the barcode error frequency.\n";
        help();
    }

    // Set `prefix_path`.
    if (out_path.length() == 0)
        out_path = ".";
    if (out_path.at(out_path.length() - 1) != '/')
        out_path += "/";
    string s = in_file.rfind('/') == string::npos ? in_file : in_file.substr(in_file.rfind('/') + 1);
    s = s.substr(0, s.rfind('.'));
    prefix_path = out_path + s;

    return 0;
}

void version() {
    std::cerr << "pstacks " << VERSION << "\n\n";

    exit(0);
}

void help() {
    std::cerr << "pstacks " << VERSION << "\n"
              << "pstacks -t file_type -f file_path [-o path] [-i id] [-m min_cov] [-p num_threads] [-h]" << "\n"
              << "  t: input file type (optional). Supported types: bowtie, sam, or bam.\n"
              << "  f: input file path.\n"
              << "  o: output path to write results.\n"
              << "  i: SQL ID to insert into the output to identify this sample.\n"
              << "  m: minimum depth of coverage to report a stack (default 1).\n"
              << "  p: enable parallel execution with num_threads threads.\n"
              << "  h: display this help messsage.\n"
              << "  Model options:\n"
              << "    --model_type <type>: either 'snp' (default), 'bounded', or 'fixed'\n"
              << "    For the SNP or Bounded SNP model:\n"
              << "      --alpha <num>: chi square significance level required to call a heterozygote or homozygote, either 0.1, 0.05 (default), 0.01, or 0.001.\n"
              << "    For the Bounded SNP model:\n"
              << "      --bound_low <num>: lower bound for epsilon, the error rate, between 0 and 1.0 (default 0).\n"
              << "      --bound_high <num>: upper bound for epsilon, the error rate, between 0 and 1.0 (default 1).\n"
              << "    For the Fixed model:\n"
              << "      --bc_err_freq <num>: specify the barcode error frequency, between 0 and 1.0.\n";

    exit(0);
}
