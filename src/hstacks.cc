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
// hstacks -- find homologous stacks among a set of samples
//
// Match stacks between samples to identify homologous loci. Stacks
// may contain masked sites (N's), resulting from using the fixed-model in
// the ustacks program, or an explicit number of mismatches per tag may be
// specified, allowing non-exact matching homologus sites to be identified
// across a set of samples.
//
// Julian Catchen
// jcatchen@uoregon.edu
// University of Oregon
//
// $Id$
//

#include "hstacks.h"

// Global variables to hold command-line options.
string in_path;
string out_path;
int    batch_id        = 0;
int    num_threads     = 1;
int    stack_depth_min = 1;
int    stack_dist      = 0;
int    n_limit         = 4;

int main (int argc, char* argv[]) {
    IF_NDEBUG_TRY

    parse_command_line(argc, argv);

    //
    // Set the number of OpenMP parallel threads to execute.
    //
    #ifdef _OPENMP
    omp_set_num_threads(num_threads);
    #endif

    vector<string>           input_files;
    vector<string>::iterator in_file;
    map<int, HLocus *>       samples;

    build_file_list(in_path, input_files);

    int id = 0;

    for (in_file = input_files.begin(); in_file != input_files.end(); in_file++) {
        map<int, HLocus *> sample;
        map<int, HLocus *>::iterator it;

        size_t pos_1     = (*in_file).find_last_of("/");
        size_t pos_2     = (*in_file).find_last_of(".");
        string sample_id = (*in_file).substr(pos_1 + 1, (pos_2 - pos_1 - 1));

        bool compressed = false;
        load_loci(*in_file, sample, 0, false, compressed);

        //
        // Give each locus a unique ID among all samples
        //
        for (it = sample.begin(); it != sample.end(); it++) {
            it->second->uniq_id = id;
            samples[id]         = (*it).second;
            id++;
        }
    }

    //
    // Calculate distance between tags in different samples.
    //
    cerr << "Calculating distance between stacks...\n";

    calc_kmer_distance(samples, stack_dist);

    //
    // Write out tags matched between samples.
    //
    write_homologous_loci(samples);

    return 0;
    IF_NDEBUG_CATCH_ALL_EXCEPTIONS
}

int calc_kmer_distance(map<int, HLocus *> &loci, int stack_dist) {
    //
    // Calculate the distance (number of mismatches) between each pair
    // of Radtags. We expect all radtags to be the same length;
    //
    CatKmerHashMap kmer_map;

    map<int, HLocus *>::iterator it;
    vector<pair<allele_type, string> >::iterator allele;
    HLocus *tag_1, *tag_2;
    int i, j;

    //
    // Calculate the number of k-mers we will generate. If kmer_len == 0,
    // determine the optimal length for k-mers.
    //
    int con_len   = strlen(loci.begin()->second->con);
    int kmer_len  = determine_kmer_length(con_len, stack_dist);
    int num_kmers = con_len - kmer_len + 1;

    //
    // Calculate the minimum number of matching k-mers required for a possible sequence match.
    //
    int min_hits = calc_min_kmer_matches(kmer_len, stack_dist, con_len, false);

    //
    // If more mismatches are allowed than can be handled by the k-mer algorithm, revert
    // to the simple, slow matching algorithm.
    //
    if (min_hits <= 0) {
        cerr << "  Unable to use k-mer matching due to the number of allowed mismatches. Switching to slower algorithm...\n";
        calc_distance(loci, stack_dist);
        return 0;
    }

    cerr << "  Number of kmers per sequence: " << num_kmers << "\n";

    populate_kmer_hash(loci, kmer_map, kmer_len);

    cerr << "  " << loci.size() << " loci, " << kmer_map.size() << " elements in the kmer hash.\n";

    //
    // OpenMP can't parallelize random access iterators, so we convert
    // our map to a vector of integer keys.
    //
    vector<int> keys;
    for (it = loci.begin(); it != loci.end(); it++)
        keys.push_back(it->first);

    #pragma omp parallel private(i, j, tag_1, tag_2, allele)
    {
        #pragma omp for schedule(dynamic)
        for (i = 0; i < (int) keys.size(); i++) {
            tag_1 = loci[keys[i]];

            for (allele = tag_1->strings.begin(); allele != tag_1->strings.end(); allele++) {

                vector<char *> kmers;
                generate_kmers(allele->second.c_str(), kmer_len, num_kmers, kmers);

                map<int, vector<allele_type> > hits;
                vector<pair<allele_type, int> >::iterator map_it;
                int d;
                //
                // Lookup the occurances of each k-mer in the kmer_map
                //
                for (j = 0; j < num_kmers; j++) {

                    if (kmer_map.count(kmers[j]) > 0)
                        for (map_it  = kmer_map[kmers[j]].begin();
                             map_it != kmer_map[kmers[j]].end();
                             map_it++)
                            hits[map_it->second].push_back(map_it->first);
                }

                //
                // Free the allocated k-mers.
                //
                for (j = 0; j < num_kmers; j++)
                    delete [] kmers[j];
                kmers.clear();

                //cerr << "  Tag " << tag_1->id << " hit " << hits.size() << " kmers.\n";

                map<int, vector<allele_type> >::iterator hit_it;
                vector<allele_type>::iterator            all_it;

                //
                // Iterate through the list of hits. For each hit, total up the hits to the various alleles.
                // Any allele that has more than min_hits check its full length to verify a match.
                //
                for (hit_it = hits.begin(); hit_it != hits.end(); hit_it++) {
                    //cerr << "  Tag " << hit_it->first << " has " << hit_it->second << " hits (min hits: " << min_hits << ")\n";

                    map<allele_type, int>           allele_cnts;
                    map<allele_type, int>::iterator cnt_it;

                    for (all_it = hit_it->second.begin(); all_it != hit_it->second.end(); all_it++)
                        allele_cnts[*all_it]++;

                    for (cnt_it = allele_cnts.begin(); cnt_it != allele_cnts.end(); cnt_it++) {
                        //cerr << "      allele " << cnt_it->first << " has " << cnt_it->second << " hits\n";

                        if (cnt_it->second < min_hits) continue;

                        //cerr << "  Match found, checking full-length match\n";

                        tag_2 = loci[hit_it->first];

                        d = dist(allele->second.c_str(), tag_2, cnt_it->first);

                        if (d < 0)
                            cerr <<
                                "Unknown error calculating distance between " <<
                                tag_1->id << " and " << tag_2->id << "; query allele: " << allele->first << "\n";

                        //cerr << "    Distance: " << d << " CTAG_DIST: " << ctag_dist << "\n";

                        //
                        // Add a match to the query sequence: catalog ID, catalog allele, query allele, distance
                        //
                        if (d <= stack_dist) {
                            if (tag_1->depth < stack_depth_min ||
                                tag_2->depth < stack_depth_min)
                                continue;

                            tag_1->add_match(tag_2->uniq_id, d);
                        }
                    }
                }
            }

            // Sort the vector of distances.
            sort(tag_1->matches.begin(), tag_1->matches.end(), compare_mdist);
        }
    }

    return 0;
}

int populate_kmer_hash(map<int, HLocus *> &loci, CatKmerHashMap &kmer_map, int kmer_len) {
    map<int, HLocus *>::iterator it;
    vector<pair<allele_type, string> >::iterator allele;
    vector<char *> kmers;
    HLocus *tag;
    char   *hash_key;
    bool    exists;
    int     j;

    //
    // Break each stack down into k-mers and create a hash map of those k-mers
    // recording in which sequences they occur.
    //
    int num_kmers = strlen(loci.begin()->second->con) - kmer_len + 1;

    for (it = loci.begin(); it != loci.end(); it++) {
        tag = it->second;

        //
        // Iterate through the possible Loci alleles
        //
        for (allele = tag->strings.begin(); allele != tag->strings.end(); allele++) {
            //
            // Generate and Hash the kmers for this allele string
            //
            generate_kmers(allele->second.c_str(), kmer_len, num_kmers, kmers);

            for (j = 0; j < num_kmers; j++) {

                exists = kmer_map.count(kmers[j]) == 0 ? false : true;

                if (exists) {
                    hash_key = kmers[j];
                } else {
                    hash_key = new char [strlen(kmers[j]) + 1];
                    strcpy(hash_key, kmers[j]);
                }

                kmer_map[hash_key].push_back(make_pair(allele->first, tag->uniq_id));
            }

            for (j = 0; j < num_kmers; j++)
                delete [] kmers[j];
            kmers.clear();
        }
    }

    //dump_kmer_map(kmer_map);

    return 0;
}

int calc_distance(map<int, HLocus *> &loci, int utag_dist) {
    //
    // Calculate the distance (number of mismatches) between each pair
    // of Radtags. We expect all radtags to be the same length;
    //
    map<int, HLocus *>::iterator it;
    HLocus *tag_1, *tag_2;
    int i, j;

    cerr << "Calculating distance between stacks...\n";

    // OpenMP can't parallelize random access iterators, so we convert
    // our map to a vector of integer keys.
    vector<int> keys;
    for (it = loci.begin(); it != loci.end(); it++)
        keys.push_back(it->first);

    #pragma omp parallel private(i, j, tag_1, tag_2)
    {
        #pragma omp for schedule(dynamic)
        for (i = 0; i < (int) keys.size(); i++) {

            tag_1 = loci[keys[i]];

            int d;

            for (j = 0; j < (int) keys.size(); j++) {
                tag_2 = loci[keys[j]];

                // Don't compare tag_1 against itself.
                if (tag_1 == tag_2)
                    continue;

                d = dist(tag_1, tag_2);

                //
                // Store the distance between these two sequences if it is
                // below the maximum distance.
                //
                if (d == utag_dist) {
                    if (tag_1->depth < stack_depth_min ||
                        tag_2->depth < stack_depth_min)
                        continue;

                    tag_1->add_match(tag_2->uniq_id, d);
                }
            }

            // Sort the vector of distances.
            sort(tag_1->matches.begin(), tag_1->matches.end(), compare_mdist);
        }
    }

    return 0;
}

int dist(HLocus *tag_1, HLocus *tag_2) {
    int   dist = 0;
    char *p    = tag_1->con;
    char *q    = tag_2->con;
    char *end  = p + strlen(p);

    // Count the number of characters that are different
    // between the two sequences. Don't count wildcard 'N'
    // nucleotides.
    while (p < end) {
        dist += ((*p == *q) || (*q == 'N' || *p == 'N')) ? 0 : 1;
        p++;
        q++;
    }

    return dist;
}

bool compare_mdist(Match *a, Match *b) {
    return (a->dist < b->dist);
}

int call_consensus(map<int, HLocus *> &loci, set<int> &merge_list,
                   string &consensus, vector<SNP *> &snps, vector<string> &alleles) {
    //
    // Create a two-dimensional array, each row containing one read. For
    // each unique tag that has been merged together, add the sequence for
    // that tag into our array as many times as it originally occurred.
    //
    HLocus *tag;
    set<int>::iterator j;
    vector<char *> reads;

    for (j = merge_list.begin(); j != merge_list.end(); j++) {
        tag = loci[*j];

        reads.push_back(tag->con);
    }

    //
    // Iterate over each column of the array and call the consensus base.
    //
    int row, col;
    int length = strlen(reads[0]);
    int height = reads.size();

    char *base;

    for (col = 0; col < length; col++) {
        vector<pair<char, int> > nuc;
        nuc.push_back(make_pair('A', 0));
        nuc.push_back(make_pair('C', 0));
        nuc.push_back(make_pair('G', 0));
        nuc.push_back(make_pair('T', 0));

        for (row = 0; row < height; row++) {
            base = reads[row];
            base = base + col;
            //cerr << "    Row: " << row << " Col: " << col << " Base: " << *base << "\n";

            switch(*base) {
            case 'A':
                nuc[0].second++;
                break;
            case 'C':
                nuc[1].second++;
                break;
            case 'G':
                nuc[2].second++;
                break;
            case 'T':
                nuc[3].second++;
                break;
            default:
                break;
            }
        }

        //cerr << "A: " << nuc[0].second << " C: " << nuc[1].second << " G: " << nuc[2].second << " T: " << nuc[3].second << "\n";

        //
        // Find the base with a plurality of occurances and call it.
        //
        sort(nuc.begin(), nuc.end(), compare_pair);

        consensus += nuc[0].second > 0 ? nuc[0].first : 'N';

        //
        // If the nucleotides are not fixed record a SNP.
        //
        if (nuc[0].second > 0 && nuc[1].second > 0) {
            SNP *s = new SNP;
            s->col    = col;
            s->lratio = 0;
            s->rank_1 = nuc[0].first;
            s->rank_2 = nuc[1].first;

            snps.push_back(s);
        }
    }

    if (!call_alleles(reads, snps, alleles)) {
        cerr << "Error calling alleles.\n";
        exit(1);
    }

    return 0;
}

int call_alleles(vector<char *> &reads, vector<SNP *> &snps, vector<string> &alleles) {
    int     row;
    int     height = reads.size();
    string  allele;
    char   *base;
    vector<SNP *>::iterator snp;

    if (snps.size() == 0)
        return 1;

    for (row = 0; row < height; row++) {
        allele.clear();

        for (snp = snps.begin(); snp != snps.end(); snp++) {
            base    = reads[row];
            base    = base + (*snp)->col;
            allele += *base;
        }

        if (allele.size() == snps.size())
            alleles.push_back(allele);
        else
            return 0;
    }

    return 1;
}

int write_homologous_loci(map<int, HLocus *> &samples) {
    map<int, HLocus *>::iterator i;
    vector<int>::iterator k;
    set<int>   write_map;
    HLocus *tag_1, *tag_2;

    cerr << "Writing homologous stacks...\n";

    //
    // Parse the input file name to create the output files
    //
    stringstream prefix;
    string out_file;
    prefix << out_path << "batch_" << batch_id;

    // Open the output files for writing.
    out_file = prefix.str() + ".homologous.nucs.tsv";
    ofstream nuc_file(out_file.c_str());
    out_file = prefix.str() + ".homologous.snps.tsv";
    ofstream snp_file(out_file.c_str());
    out_file = prefix.str() + ".homologous.alleles.tsv";
    ofstream all_file(out_file.c_str());
    out_file = prefix.str() + ".homologous.matches.tsv";
    ofstream mat_file(out_file.c_str());
    out_file = prefix.str() + ".homologous.tags.tsv";
    ofstream tag_file(out_file.c_str());

    int id = 1;

    for (i = samples.begin(); i != samples.end(); i++) {
        tag_1 = i->second;

        //
        // This tag may already have been merged by an earlier operation.
        //
        if (write_map.find(tag_1->uniq_id) != write_map.end())
            continue;

        set<int> unique_merge_list;
        set<string> unique_alleles;
        set<int>::iterator it;

        trace_stack_graph(tag_1, samples, unique_merge_list);

        //
        // Call the consensus for this locus and identify SNPs and associated alleles.
        //
        string       consensus;
        vector<SNP *>     snps;
        vector<string> alleles;

        call_consensus(samples, unique_merge_list, consensus, snps, alleles);

        //
        // Output the consensus tag for a locus in this sample.
        //
        tag_file <<
            "0"              << "\t" <<
            batch_id         << "\t" <<
            id               << "\t" <<
            tag_1->loc.chr() << "\t" <<
            tag_1->loc.bp    << "\t" <<
            "consensus"      << "\t" <<
            0                << "\t" <<
            ""               << "\t" <<
            consensus        << "\t" <<
            0                << "\t" <<  // These flags are unused in hstacks, but important in ustacks
            0                << "\t" <<
            0                << "\n";

        //
        // Output the SNPs and alleles
        //
        string allele;
        vector<SNP *>::iterator  s;
        set<string>::iterator    u;

        for (s = snps.begin(); s != snps.end(); s++)
            snp_file <<
                "0"          << "\t" <<
                batch_id     << "\t" <<
                id           << "\t" <<
                (*s)->col    << "\t" <<
                (*s)->lratio << "\t" <<
                (*s)->rank_1 << "\t" <<
                (*s)->rank_2 << "\n";

        for (uint a = 0; a < alleles.size(); a++)
            unique_alleles.insert(alleles[a]);

        for (u = unique_alleles.begin(); u != unique_alleles.end(); u++)
            all_file <<
                "0"        << "\t" <<
                batch_id   << "\t" <<
                id         << "\t" <<
                *u         << "\t" <<
                0          << "\t" <<
                0          << "\n";

        unique_alleles.clear();

        int sub_id = 0;
        int a      = 0;

        for (it = unique_merge_list.begin(); it != unique_merge_list.end(); it++) {
            tag_2 = samples[(*it)];

            // Record the nodes that have been merged in this round.
            write_map.insert(tag_2->uniq_id);

            //
            // For each tag we are outputting, output the depth of coverage for each
            // nucleotide in that stack separately (in order to calculate correlations
            // between depth of coverage and fixed/non-fixed nucleotides.
            //
            if (unique_merge_list.size() > 1) {
                char *p, *end;
                end = tag_2->con + strlen(tag_2->con);
                for (p = tag_2->con; p < end; p++)
                    nuc_file
                        << tag_2->sample_id << "_" << tag_2->id << "\t"
                        << *p << "\t"
                        << tag_2->depth << "\n";
            }

            //
            // Output the consensus sequenes for all homologous loci.
            //
            tag_file <<
                "0"              << "\t" <<
                batch_id         << "\t" <<
                id               << "\t" <<
                tag_2->loc.chr() << "\t" <<
                tag_2->loc.bp    << "\t" <<
                "primary"        << "\t" <<
                sub_id           << "\t" <<
                tag_2->sample_id << "_"  <<  tag_2->id << "\t" <<
                tag_2->con       << "\t" <<
                ""               << "\t" <<  // These flags are unused in hstacks, but important in ustacks
                ""               << "\t" <<
                ""               << "\n";

            allele = (alleles.size() == 0) ? "consensus" : alleles[a];

            mat_file <<
                "0"              << "\t" <<
                batch_id         << "\t" <<
                id               << "\t" <<
                tag_2->sample_id << "\t" <<
                tag_2->uniq_id   << "\t" <<
                allele           << "\n";

            sub_id++;
            a++;
        }

        id++;
    }

    tag_file.close();
    nuc_file.close();
    snp_file.close();
    all_file.close();
    mat_file.close();

    cerr << "  Wrote " << id - 1 << " tags to " << out_file << "\n";

    return 0;
}

int trace_stack_graph(HLocus *tag_1, map<int, HLocus *> &loci, set<int> &unique_merge_list) {
    queue<int>                    merge_list;
    pair<set<int>::iterator,bool> ret;
    vector<Match *>::iterator     k;
    HLocus *tag_2;

    unique_merge_list.insert(tag_1->uniq_id);
    merge_list.push(tag_1->uniq_id);

    while (!merge_list.empty()) {
        tag_2 = loci[merge_list.front()];
        merge_list.pop();

        for (k = tag_2->matches.begin(); k != tag_2->matches.end(); k++) {
            ret = unique_merge_list.insert((*k)->cat_id);

            //
            // If this Tag has not already been added to the merge list (i.e. we were able
            // to insert it in to our unique_merge_list, which is a set), add it for consideration
            // later in the loop.
            //
            if (ret.second == true)
                merge_list.push((*k)->cat_id);
        }
    }

    return 0;
}

int build_file_list(string in_path, vector<string> &sql_files) {
    struct dirent *dirp;
    string d;

    DIR *dp = opendir(in_path.c_str());

    if (dp == NULL) {
        cerr << "Error (" << errno << ") opening " << in_path << "\n";
        return errno;
    }

    while ((dirp = readdir(dp)) != NULL) {
        d     = string(dirp->d_name);

        if (d.find("tags.tsv") != string::npos &&
            d.find("batch") == string::npos) {
            size_t pos = d.find(".tags.tsv");
            d = in_path + d.substr(0, pos);
            sql_files.push_back(d);
        }
    }

    closedir(dp);

    cerr << "Identified " << sql_files.size() << " samples.\n";

    return 0;
}

HLocus::~HLocus()
{
    vector<Match *>::iterator it;

    for (it = this->matches.begin(); it != this->matches.end(); it++)
        delete *it;
}

int
HLocus::add_match(int id, int distance)
{
    Match *m = new Match;
    m->cat_id = id;
    m->dist   = distance;

    this->matches.push_back(m);

    return 0;
}

int
HLocus::populate_alleles()
{
    this->strings.clear();

    string s;
    int    j;
    vector<SNP *>::iterator    s_it;
    map<string, int>::iterator a;

    if (this->strings.size() == 0)
        this->strings.push_back(make_pair("consensus", this->con));

    else
        for (a = this->alleles.begin(); a != this->alleles.end(); a++) {
            s = this->con;
            j = 0;

            for (s_it = this->snps.begin(); s_it != this->snps.end(); s_it++) {
                s.replace((*s_it)->col, 1, 1, a->first[j]);
                j++;
            }

            this->strings.push_back(make_pair(a->first, s));
        }

    //
    // Count the number of N's
    //
    vector<int> col;
    char *p     = this->con;
    int   i     = 0;
    while (*p != '\0') {
        if (*p == 'N')
            col.push_back(i);
        i++;
        p++;
    }

    int n_cnt = col.size();

    if (n_cnt == 0) return 0;

    //
    // If there are too many Ns in this stack, do not include it in the
    // search.
    //
    if (n_cnt > n_limit) {
        this->strings.clear();
        return 0;
    }

    //
    // Generate all permutations of strings for n_cnt N's
    //
    if (pstrings.count(n_cnt) == 0)
        generate_permutations(pstrings, n_cnt);

    vector<pair<allele_type, string> > new_strings;
    vector<pair<allele_type, string> >::iterator k;
    vector<int>::iterator c;
    char  *q, **r;
    int    n = (int) pow(4, n_cnt);

    for (k = this->strings.begin(); k != this->strings.end(); k++) {

        for (i = 0; i < n; i++) {
            r = pstrings[n_cnt];
            q = r[i];

            s = k->second;
            j = 0;
            for (c = col.begin(); c != col.end(); c++) {
                //cerr << "Str: " << s << "; rep str: " << q << "; replacing col: " << *c << " with '" << q[j] << "'\n";
                s.replace(*c, 1, 1, q[j]);
                j++;
            }

            new_strings.push_back(make_pair(k->first, s));
        }
    }

    this->strings.clear();
    for (k = new_strings.begin(); k != new_strings.end(); k++)
        this->strings.push_back(*k);

    return 0;
}

int parse_command_line(int argc, char* argv[]) {
    int c;

    while (1) {
        static struct option long_options[] = {
            {"help",        no_argument,       NULL, 'h'},
            {"version",     no_argument,       NULL, 'v'},
            {"stack_dist",  required_argument, NULL, 'n'},
            {"depth_min",   required_argument, NULL, 'm'},
            {"inpath",      required_argument, NULL, 'p'},
            {"outpath",     required_argument, NULL, 'o'},
            {"n_limit",     required_argument, NULL, 'N'},
            {"batch_id",    required_argument, NULL, 'b'},
            {0, 0, 0, 0}
        };

        // getopt_long stores the option index here.
        int option_index = 0;

        c = getopt_long(argc, argv, "hvi:p:o:b:e:m:n:N:", long_options, &option_index);

        // Detect the end of the options.
        if (c == -1)
            break;

        switch (c) {
        case 'h':
            help();
            break;
        case 'v':
            version();
            break;
        case 'i':
            in_path = optarg;
            break;
        case 'o':
            out_path = optarg;
            break;
        case 'b':
            batch_id = atoi(optarg);
            break;
        case 'N':
            n_limit = atoi(optarg);
            break;
        case 'm':
            stack_depth_min = atoi(optarg);
            break;
        case 'n':
            stack_dist = atoi(optarg);
            break;
        case 'p':
            num_threads = atoi(optarg);
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

    if (in_path.length() == 0) {
        cerr << "You must specify a path to a set of input files.\n";
        help();
    }

    if (in_path.at(in_path.length() - 1) != '/')
        in_path += "/";

    if (out_path.length() == 0)
        out_path = ".";

    if (out_path.at(out_path.length() - 1) != '/')
        out_path += "/";

    return 0;
}

void version() {
    cerr << "hstacks " << VERSION << "\n\n";

    exit(1);
}

void help() {
    cerr << "hstacks " << VERSION << "\n"
              << "hstacks -i path [-o path] [-b batch_id] [-n mismatches] [-m min] [-p min_threads] [-N limit] [-h]" << "\n"
              << "  i: path to the set of SQL files from which to load loci." << "\n"
              << "  o: output path to write results." << "\n"
              << "  b: SQL Batch ID to insert into the output to identify a group of samples." << "\n"
              << "  m: minimum stack depth required for a locus to be included in the search." << "\n"
              << "  n: number of mismatches to allow between stacks." << "\n"
              << "  N: number of 'N' characters to allow in a stack (default: 4)." << "\n"
              << "  p: enable parallel execution with num_threads threads.\n"
              << "  h: display this help messsage." << "\n\n";

    exit(1);
}
