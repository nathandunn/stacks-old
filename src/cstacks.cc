// -*-mode:c++; c-style:k&r; c-basic-offset:4;-*-
//
// Copyright 2010-2018, Julian Catchen <jcatchen@illinois.edu>
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
// cstacks -- Create a catalog of Stacks.
//

#include "MetaPopInfo.h"
#include "cstacks.h"

// Global variables to hold command-line options.
queue<pair<int, string> > samples;
string  out_path;
string  catalog_path;
FileT   in_file_type      = FileT::sql;
int     ctag_dist         = 1;
bool    set_kmer_len      = true;
int     kmer_len          = 0;
searcht search_type       = sequence;
int     num_threads       = 1;
bool    mult_matches      = false;
bool    report_mmatches   = false;
bool    require_uniq_haplotypes = false;
bool    gapped_alignments = true;
double  min_match_len     = 0.80;
double  max_gaps          = 2.0;

using namespace std;

int main (int argc, char* argv[]) {
    IF_NDEBUG_TRY

    parse_command_line(argc, argv);

    uint sample_cnt = samples.size();

    cerr << "cstacks parameters selected:\n"
         << "  Loci matched based on sequence identity.\n"
         << "  Number of mismatches allowed between stacks: " << ctag_dist << "\n"
         << "  Gapped alignments: " << (gapped_alignments ? "enabled" : "disabled") << "\n"
         << "Constructing catalog from " << sample_cnt << " samples.\n";

    //
    // Set the number of OpenMP parallel threads to execute.
    //
    #ifdef _OPENMP
    omp_set_num_threads(num_threads);
    #endif

    map<int, CLocus *> catalog;
    map<int, CLocus *>::iterator cat_it;
    map<int, QLocus *>::iterator query_it;
    pair<int, string> s;
    bool compressed = false;
    int  i;

    set<int> seen_sample_ids; // For checking sample ID unicity.

    if (catalog_path.length() > 0) {
        cerr << "\nInitializing existing catalog...\n";
        if (!initialize_existing_catalog(catalog_path, catalog)) {
            cerr << "Error: Failed to initialize the catalog.\n";
            throw exception();
        }
        i = 1;

    } else {
        s = samples.front();
        samples.pop();

        cerr << "\nInitializing new catalog...\n";
        if (!initialize_new_catalog(s, catalog)) {
            cerr << "Error: Failed to initialize the catalog.\n";
            throw exception();
        }
        i = 2;
        seen_sample_ids.insert(catalog.begin()->second->sample_id);
    }

    //
    // Build an index of the catalog
    //
    map<string, int> cat_index;
    if (search_type == genomic_loc) {
        cerr << "Building an index of the catalog.\n";
        update_catalog_index(catalog, cat_index);
    }

    while (!samples.empty()) {
        map<int, QLocus *> sample;

        s = samples.front();
        samples.pop();

        cerr << "\nProcessing sample " << s.second << " [" << i << " of " << sample_cnt << "]\n";

        if (!load_loci(s.second, sample, 0, false, compressed)) {
            cerr << "Failed to load sample " << i << "\n";
            continue;
        }
        //
        // Assign the ID for this sample data.
        //
        s.first = sample.begin()->second->sample_id;

        // Check for unique sample IDs.
        if (!seen_sample_ids.insert(s.first).second) {
            // Insert failed.
            cerr << "Error: Sample ID '" << s.first << "' occurs more than once. Sample IDs must be unique." << endl;
            throw exception();
        }

        //dump_loci(sample);

        cerr << "Searching for sequence matches...\n";
        find_kmer_matches_by_sequence(catalog, sample, ctag_dist);

        if (gapped_alignments) {
            cerr << "Searching for gapped alignments...\n";
            search_for_gaps(catalog, sample, min_match_len, ctag_dist);
        }

        cerr << "Merging matches into catalog...\n";
        uint mmatches = 0;
        uint gmatches = 0;
        uint umatches = 0;
        uint nmatches = 0;
        merge_matches(catalog, sample, s, ctag_dist, nmatches, umatches, gmatches, mmatches);
        cerr << "  " << umatches << " loci were matched to a catalog locus.\n"
             << "  " << gmatches << " loci were matched to a catalog locus using gapped alignments.\n"
             << "  " << nmatches << " loci were newly added to the catalog.\n"
             << "  " << mmatches << " loci matched more than one catalog locus and were excluded.\n";

        //
        // Regenerate the alleles for the catalog tags after merging the new sample into the catalog.
        //
        for (cat_it = catalog.begin(); cat_it != catalog.end(); cat_it++) {
            cat_it->second->populate_alleles();
            cat_it->second->match_cnt = 0;
        }

        i++;

        for (query_it = sample.begin(); query_it != sample.end(); query_it++)
            delete query_it->second;
        sample.clear();
    }

    cerr << "\nWriting catalog in directory '" << out_path << "'.\n";
    write_catalog(catalog);
    cerr << "cstacks is done.\n";

    //
    // Free memory associated with the catalog.
    //
    for (cat_it = catalog.begin(); cat_it != catalog.end(); cat_it++)
        delete cat_it->second;
    catalog.clear();

    return 0;
    IF_NDEBUG_CATCH_ALL_EXCEPTIONS
}

int update_catalog_index(map<int, CLocus *> &catalog, map<string, int> &cat_index) {
    map<int, CLocus *>::iterator j;
    char id[id_len];

    for (j = catalog.begin(); j != catalog.end(); j++) {
        snprintf(id, id_len - 1, "%s|%d|%c",
                 j->second->loc.chr(),
                 j->second->loc.bp,
                 j->second->loc.strand == strand_plus ? '+' : '-');

        if (cat_index.count(id) == 0) {
            cat_index[id] = j->first;
        } else {
            if (cat_index[id] != j->first)
                cerr << "Error: Catalog index mismatch, key: '" << id << "'.\n";
        }
    }

    return 0;
}

int
characterize_mismatch_snps(CLocus *catalog_tag, QLocus *query_tag)
{
    set<int> snp_cols;
    uint i;
    for (i = 0; i < catalog_tag->snps.size(); i++)
        snp_cols.insert(catalog_tag->snps[i]->col);
    for (i = 0; i < query_tag->snps.size(); i++)
        snp_cols.insert(query_tag->snps[i]->col);

    //
    // For each mismatch found, create a SNP object
    //
    const char *c        = catalog_tag->con;
    const char *c_beg    = c;
    const char *c_end    = c + strlen(c);
    const char *q        = query_tag->con;
    const char *q_beg    = q;
    const char *q_end    = q + strlen(q);

    i = 0;
    while (c < c_end && q < q_end) {
        if (snp_cols.count(i) == 0 &&
            (*c != *q) && (*c != 'N' && *q != 'N')) {

            cerr << "Adding a new SNP at position " << c - c_beg << ", " << *c << "/" << *q << "\n";
            SNP *s = new SNP;
            s->type   = snp_type_het;
            s->col    = c - c_beg;
            s->lratio = 0;
            s->rank_1 = *c;
            s->rank_2 = *q;

            merge_allele(catalog_tag, s);
            merge_allele(query_tag, s);

            catalog_tag->snps.push_back(s);

            s = new SNP;
            s->type   = snp_type_het;
            s->col    = q - q_beg;
            s->lratio = 0;
            s->rank_1 = *q;
            s->rank_2 = *c;

            query_tag->snps.push_back(s);
        }
        c++;
        q++;
        i++;
    }

    return 1;
}

int
merge_matches(map<int, CLocus *> &catalog, map<int, QLocus *> &sample, pair<int, string> &sample_file, int ctag_dist,
              uint &new_matches, uint &unique_matches, uint &gapped_matches, uint &multiple_matches)
{
    map<int, QLocus *>::iterator i;
    CLocus *ctag;
    QLocus *qtag;
    string  cseq, qseq, cigar_str;
    int     cseq_len, match_index=-1;
    vector<pair<char, uint> > cigar;

    GappedAln *aln = new GappedAln();

    for (i = sample.begin(); i != sample.end(); i++) {
        qtag = i->second;

        //
        // If this stack didn't match an existing catalog stack, add this stack to the
        // catalog as a new stack.
        //
        if (qtag->matches.size() == 0) {
            add_unique_tag(sample_file, catalog, qtag);
            new_matches++;
            continue;
        }

        //
        // Check for multiple matches. We will reduce the list of Match objects, which
        // contain matches to multiple alleles for a single locus, to the smallest distance
        // for a locus.
        //
        map<int, uint> local_matches;
        for (uint k = 0; k < qtag->matches.size(); k++) {
            if (local_matches.count(qtag->matches[k]->cat_id) == 0)
                local_matches[qtag->matches[k]->cat_id] = qtag->matches[k]->dist;
            else if (qtag->matches[k]->dist < local_matches[qtag->matches[k]->cat_id])
                local_matches[qtag->matches[k]->cat_id] = qtag->matches[k]->dist;
        }

        uint min_dist    =  1000;
        uint num_matches =  0;
        int  min_cat_id  = -1;

        //
        // Find the minimum distance and then check how many matches have that distance.
        //
        for (map<int, uint>::iterator j = local_matches.begin(); j != local_matches.end(); j++)
            min_dist = j->second < min_dist ? j->second : min_dist;

        for (map<int, uint>::iterator j = local_matches.begin(); j != local_matches.end(); j++)
            if (j->second == min_dist) {
                num_matches++;
                min_cat_id = j->first;
            }

        //
        // Emit a warning if the query tag matches more than one tag in the catalog.
        //
        if (num_matches > 1) {
            multiple_matches++;
            if (report_mmatches) {
                cerr <<
                    "  Warning: sample " << sample_file.second << ", tag " << qtag->id <<
                    ", matches more than one tag in the catalog and was excluded: ";
                for (map<int, uint>::iterator j = local_matches.begin(); j != local_matches.end(); j++)
                    cerr << j->first << " ";
                cerr << "\n";
            }
            //
            // Don't record matches to multiple catalog entries unless instructed
            // to do so by the command line option.
            //
            if (!mult_matches) continue;
        }

        ctag = catalog[min_cat_id];

        if (ctag == NULL)
            cerr << "  Unable to locate catalog tag " << min_cat_id << "\n";

        cigar_str = "";

        for (uint k = 0; k < qtag->matches.size(); k++)
            if ((int)qtag->matches[k]->cat_id == min_cat_id) {
                cigar_str   = qtag->matches[k]->cigar;
                match_index = k;
                break;
            }

        //
        // If we have already matched a query locus to this catalog locus for the current
        // sample, we must re-align the sequences in case changes have been made to the
        // sequence by the previous matching sequence.
        //
        if (gapped_alignments && ctag->match_cnt > 0) {
            string query_allele, query_seq, cat_allele, cat_seq;
            // cerr << "    Warning: Catalog locus " << ctag->id
            //      << ", Sample " << qtag->sample_id << ", locus " << qtag->id
            //      << "; sequence length has changed since original alignment: "
            //      << qseq_len << " <-> " << ctag->len
            //      << "; re-aligning.\n";

            //
            // Find the proper query allele to align against the catalog. We can align
            // against the catalog consensus because the allele strings may have changed.
            //
            query_allele = qtag->matches[match_index]->query_type;
            for (uint k = 0; k < qtag->strings.size(); k++)
                if (qtag->strings[k].first == query_allele) {
                    query_seq = qtag->strings[k].second;
                    break;
                }
            aln->init(ctag->len, qtag->len);
            aln->align(ctag->con, query_seq);
            cigar_str = invert_cigar(aln->result().cigar);
        }

        bool gapped_aln = false;
        if (cigar_str.length() > 0)
            gapped_aln = true;

        //
        // If the match was a gapped alignment, adjust the lengths of the consensus sequences.
        // Adjust the postition of any SNPs that were shifted down sequence due to a gap.
        //
        if (gapped_aln) {
            cseq_len  = parse_cigar(cigar_str.c_str(), cigar);
            qseq      = apply_cigar_to_seq(qtag->con, cigar);
            adjust_snps_for_gaps(cigar, qtag);

            cigar_str = invert_cigar(cigar_str);
            cseq_len  = parse_cigar(cigar_str.c_str(), cigar);
            cseq      = apply_cigar_to_seq(ctag->con, cigar);
            adjust_snps_for_gaps(cigar, ctag);

            if (qseq.length() != cseq.length())
                cerr << "Sample ID: " << qtag->sample_id << "; Query ID: " << qtag->id << "; catalog ID: " << ctag->id << ";\n"
                     << "qloc: " << qtag->con << "\n"
                     << "qseq: " << qseq << "\n"
                     << "cloc: " << ctag->con << " [" << cigar_str << "]\n"
                     << "cseq: " << cseq << "\n";

            //
            // If the alignment modified the catalog locus, record it so we can re-align
            // any other matching sequences from this sample.
            //
            if ((uint)cseq_len > ctag->len)
                ctag->match_cnt++;

            //
            // Adjust the consensus sequences for both loci.
            //
            ctag->add_consensus(cseq.c_str());
            qtag->add_consensus(qseq.c_str());

            gapped_matches++;

        } else {
            unique_matches++;
            ctag->match_cnt++;
        }

        //
        // If mismatches are allowed between query and catalog tags, identify the
        // mismatches and convert them into SNP objects to be merged into the catalog tag.
        //
        if ((ctag_dist > 0 || search_type == genomic_loc) && !characterize_mismatch_snps(ctag, qtag))
            cerr
                << "  Error characterizing mismatch SNPs "
                << sample_file.second << ", tag " << qtag->id
                << " with catalog tag " << ctag->id << "\n";

        //
        // Merge the SNPs and alleles from the sample into the catalog tag.
        //
        if (!ctag->merge_snps(qtag)) {
            cerr << "Error merging " << sample_file.second << ", tag " << qtag->id
                 << " with catalog tag " << ctag->id << "\n";
        }

        //
        // Add any new sequence information into the catalog consensus.
        //
        if (gapped_aln) {
            for (uint k = 0; k < ctag->len && k < qtag->len; k++)
                if (qtag->con[k] != 'N' && ctag->con[k] == 'N')
                    ctag->con[k] = qtag->con[k];

        } else if (strlen(ctag->con) < strlen(qtag->con)) {
            ctag->add_consensus(qtag->con);
        }

        ctag->sources.push_back(make_pair(sample_file.first, qtag->id));
    }

    delete aln;

    return 0;
}

int add_unique_tag(pair<int, string> &sample_file, map<int, CLocus *> &catalog, QLocus *qloc) {
    vector<SNP *>::iterator i;
    map<string, int>::iterator j;

    int cid = catalog.size();

    CLocus *c = new CLocus;
    c->id = cid + 1;
    c->add_consensus(qloc->con);
    //
    // Record the source of this catalog tag.
    //
    c->sources.push_back(make_pair(sample_file.first, qloc->id));
    //
    // Add the physical genome location of this locus.
    //
    c->loc.set(qloc->loc.chr(), qloc->loc.bp, qloc->loc.strand);

    catalog[c->id] = c;

    // cerr << "Adding sample: " << qloc->id << " to the catalog as ID: " << c->id << "\n";

    for (i = qloc->snps.begin(); i != qloc->snps.end(); i++) {
        SNP *snp    = new SNP;
        snp->col    = (*i)->col;
        snp->type   = (*i)->type;
        snp->lratio = (*i)->lratio;
        snp->rank_1 = (*i)->rank_1;
        snp->rank_2 = (*i)->rank_2;

        c->snps.push_back(snp);
    }

    for (j = qloc->alleles.begin(); j != qloc->alleles.end(); j++) {
        c->alleles[j->first] = j->second;
    }

    c->populate_alleles();

    return 0;
}

int find_kmer_matches_by_sequence(map<int, CLocus *> &catalog, map<int, QLocus *> &sample, int ctag_dist) {
    //
    // Calculate the distance (number of mismatches) between each pair
    // of Radtags. We expect all radtags to be the same length;
    //
    KmerHashMap                       kmer_map;
    map<int, pair<allele_type, int> > allele_map;
    vector<char *>                    kmer_map_keys;
    map<int, QLocus *>::iterator      it;
    vector<pair<allele_type, string> >::iterator allele;
    QLocus *tag_1;
    CLocus *tag_2;

    // OpenMP can't parallelize random access iterators, so we convert
    // our map to a vector of integer keys.
    vector<int> keys;
    for (it = sample.begin(); it != sample.end(); it++)
        keys.push_back(it->first);

    //
    // Calculate the number of k-mers we will generate. If kmer_len == 0,
    // determine the optimal length for k-mers.
    //
    int con_len   = strlen(sample[keys[0]]->con);
    if (set_kmer_len) kmer_len = determine_kmer_length(con_len, ctag_dist);
    int num_kmers = con_len - kmer_len + 1;

    //
    // Calculate the minimum number of matching k-mers required for a possible sequence match.
    //
    int min_hits = calc_min_kmer_matches(kmer_len, ctag_dist, con_len, set_kmer_len ? true : false);

    cerr << "  Distance allowed between stacks: " << ctag_dist
         << "; searching with a k-mer length of " << kmer_len << " (" << num_kmers << " k-mers per read); "
         << min_hits << " k-mer hits required.\n";

    // clock_t time_1, time_2, time_3, time_4;
    // double  per_locus = 0.0;

    // time_1 = clock();
    populate_kmer_hash(catalog, kmer_map, kmer_map_keys, allele_map, kmer_len);
    // time_2 = clock();

    cerr << "  " << catalog.size() << " loci in the catalog, " << kmer_map.size() << " kmers in the catalog hash.\n";

    #pragma omp parallel private(tag_1, tag_2, allele)
    {
        KmerHashMap::iterator   h;
        vector<char *>          kmers;
        set<string>             uniq_kmers;
        vector<int>             hits;
        vector<pair<int, int> > ordered_hits;
        uint                    hit_cnt, index, prev_id, allele_id, hits_size;
        int                     d;
        pair<allele_type, int>  cat_hit;

        initialize_kmers(kmer_len, num_kmers, kmers);

        #pragma omp for
        for (uint i = 0; i < keys.size(); i++) {
            tag_1 = sample[keys[i]];

            // time_3 = clock();

            for (allele = tag_1->strings.begin(); allele != tag_1->strings.end(); allele++) {

                generate_kmers_lazily(allele->second.c_str(), kmer_len, num_kmers, kmers);

                //
                // We want to create a list of unique kmers to search with; otherwise, repetitive kmers will
                // generate, multiple, spurious hits in sequences with multiple copies of the same kmer.
                //
                uniq_kmers.clear();
                for (int j = 0; j < num_kmers; j++)
                    uniq_kmers.insert(kmers[j]);

                hits.clear();
                ordered_hits.clear();

                //
                // Lookup the occurances of each k-mer in the kmer_map
                //
                for (set<string>::iterator j = uniq_kmers.begin(); j != uniq_kmers.end(); j++) {

                    h = kmer_map.find(j->c_str());

                    if (h != kmer_map.end())
                        for (uint k = 0; k <  h->second.size(); k++)
                            hits.push_back(h->second[k]);
                }

                //
                // Sort the vector of indexes; provides the number of hits to each allele/locus
                // and orders them largest to smallest.
                //
                sort(hits.begin(), hits.end());

                //
                // Iterate through the list of hits and collapse them down by number of kmer hits per allele.
                //
                hits_size = hits.size();

                if (hits_size == 0)
                    continue;

                prev_id   = hits[0];
                index     = 0;

                do {
                    hit_cnt   = 0;
                    allele_id = prev_id;

                    while (index < hits_size && (uint) hits[index] == prev_id) {
                        hit_cnt++;
                        index++;
                    }

                    if (index < hits_size)
                        prev_id = hits[index];

                    if (hit_cnt >= (uint) min_hits)
                        ordered_hits.push_back(make_pair(allele_id, hit_cnt));

                } while (index < hits_size);

                for (uint j = 0; j < ordered_hits.size(); j++) {
                    cat_hit = allele_map.at(ordered_hits[j].first);
                    hit_cnt = ordered_hits[j].second;

                    tag_2 = catalog[cat_hit.second];

                    d = dist(allele->second.c_str(), tag_2, cat_hit.first);

                    assert(d >= 0);

                    //
                    // Check if any of the mismatches occur at the 3' end of the read. If they
                    // do, they may indicate a frameshift is present at the 3' end of the read.
                    // If found, do not merge these tags, leave them for the gapped alignmnet
                    // algorithm.
                    //
                    if (d <= ctag_dist && check_frameshift(allele->second.c_str(), tag_2, cat_hit.first, (size_t) ctag_dist))
                        continue;

                    //
                    // Add a match to the query sequence: catalog ID, catalog allele, query allele, distance
                    //
                    if (d <= ctag_dist)
                        tag_1->add_match(tag_2->id, cat_hit.first, allele->first, d);
                }
            }

            // Sort the vector of distances.
            sort(tag_1->matches.begin(), tag_1->matches.end(), compare_matches);

            // time_4 = clock();
            // per_locus += (time_4 - time_3);
        }

        //
        // Free the allocated k-mers.
        //
        for (uint j = 0; j < kmers.size(); j++)
            delete [] kmers[j];
        kmers.clear();
    }

    // cerr << "Time to kmerize catalog: " << time_2 - time_1 << "\n"
    //      << "Average time per locus:  " << per_locus / (double) keys.size() << "\n";

    free_kmer_hash(kmer_map, kmer_map_keys);

    return 0;
}

int
search_for_gaps(map<int, CLocus *> &catalog, map<int, QLocus *> &sample, double min_match_len, double ctag_dist)
{
    //
    // Search for loci that can be merged with a gapped alignment.
    //
    KmerHashMap                       kmer_map;
    map<int, pair<allele_type, int> > allele_map;
    vector<char *>                    kmer_map_keys;
    map<int, QLocus *>::iterator      it;
    QLocus *tag_1;
    CLocus *tag_2;

    //
    // OpenMP can't parallelize random access iterators, so we convert
    // our map to a vector of integer keys.
    //
    vector<int> keys;
    for (it = sample.begin(); it != sample.end(); it++)
        keys.push_back(it->first);

    //
    // Calculate the number of k-mers we will generate. If kmer_len == 0,
    // determine the optimal length for k-mers.
    //
    int con_len   = strlen(sample[keys[0]]->con);
    int kmer_len  = 19;
    int num_kmers = con_len - kmer_len + 1;

    //
    // Calculate the minimum number of matching k-mers required for a possible sequence match.
    //
    int min_hits = (round((double) con_len * min_match_len) - (kmer_len * max_gaps)) - kmer_len + 1;

    cerr << "  Searching with a k-mer length of " << kmer_len << " (" << num_kmers << " k-mers per read); " << min_hits << " k-mer hits required.\n";

    // clock_t time_1, time_2, time_3, time_4;
    // double  per_locus = 0.0;

    // time_1 = clock();
    populate_kmer_hash(catalog, kmer_map, kmer_map_keys, allele_map, kmer_len);
    // time_2 = clock();

    #pragma omp parallel private(tag_1, tag_2)
    {
        KmerHashMap::iterator     h;
        AlignRes                  aln_res;
        vector<char *>            kmers;
        set<string>               uniq_kmers;
        vector<int>               hits;
        vector<pair<int, int> >   ordered_hits;
        uint                      hit_cnt, index, prev_id, allele_id, hits_size, stop, top_hit;
        int                       d;
        vector<pair<char, uint> > cigar;
        pair<allele_type, int>    cat_hit;
        string                    cat_seq;

        GappedAln *aln = new GappedAln();

        initialize_kmers(kmer_len, num_kmers, kmers);

        #pragma omp for schedule(dynamic)
        for (uint i = 0; i < keys.size(); i++) {
            tag_1 = sample[keys[i]];

            //
            // If we already matched this locus to the catalog without using gapped alignments, skip it now.
            //
            if (tag_1->matches.size() > 0)
                continue;

            // time_3 = clock();

            for (vector<pair<allele_type, string> >::iterator allele = tag_1->strings.begin(); allele != tag_1->strings.end(); allele++) {

                generate_kmers_lazily(allele->second.c_str(), kmer_len, num_kmers, kmers);

                //
                // We want to create a list of unique kmers to search with; otherwise, repetitive kmers will
                // generate, multiple, spurious hits in sequences with multiple copies of the same kmer.
                //
                uniq_kmers.clear();
                for (int j = 0; j < num_kmers; j++)
                    uniq_kmers.insert(kmers[j]);

                hits.clear();
                ordered_hits.clear();

                //
                // Lookup the occurances of each k-mer in the kmer_map
                //
                for (set<string>::iterator j = uniq_kmers.begin(); j != uniq_kmers.end(); j++) {

                    h = kmer_map.find(j->c_str());

                    if (h != kmer_map.end())
                        for (uint k = 0; k <  h->second.size(); k++)
                            hits.push_back(h->second[k]);
                }

                //
                // Sort the vector of indexes; provides the number of hits to each allele/locus
                // and orders them largest to smallest.
                //
                sort(hits.begin(), hits.end());

                //
                // Iterate through the list of hits and collapse them down by number of kmer hits per allele.
                //
                hits_size = hits.size();

                if (hits_size == 0)
                    continue;

                prev_id   = hits[0];
                index     = 0;

                do {
                    hit_cnt   = 0;
                    allele_id = prev_id;

                    while (index < hits_size && (uint) hits[index] == prev_id) {
                        hit_cnt++;
                        index++;
                    }

                    if (index < hits_size)
                        prev_id = hits[index];

                    if (hit_cnt >= (uint) min_hits)
                        ordered_hits.push_back(make_pair(allele_id, hit_cnt));

                } while (index < hits_size);

                if (ordered_hits.size() == 0)
                    continue;

                //
                // Process the hits from most kmer hits to least kmer hits.
                //
                sort(ordered_hits.begin(), ordered_hits.end(), compare_pair_intint);

                //
                // Only try to align the sequences with the most kmers in common.
                //
                top_hit = ordered_hits[0].second;
                stop    = 1;
                for (uint j = 1; j < ordered_hits.size(); j++)
                    if ((uint)ordered_hits[j].second < top_hit) {
                        stop = j;
                        break;
                    }

                for (uint j = 0; j < stop; j++) {
                    cat_hit = allele_map.at(ordered_hits[j].first);
                    hit_cnt = ordered_hits[j].second;

                    tag_2 = catalog[cat_hit.second];

                    cat_seq = "";
                    for (uint k = 0; k < tag_2->strings.size(); k++)
                        if (tag_2->strings[k].first == cat_hit.first) {
                            cat_seq = tag_2->strings[k].second;
                            break;
                        }

                    aln->init(tag_2->len, tag_1->len);

                    if (aln->align(cat_seq, allele->second)) {
                        cigar.clear();
                        aln->parse_cigar(cigar);
                        aln_res = aln->result();
                        d       = dist(cat_seq.c_str(), allele->second.c_str(), cigar);

                        //
                        // If the alignment has too many gaps, skip it.
                        //
                        if (aln_res.gap_cnt <= (max_gaps + 1)) {
                            //
                            // If the alignment doesn't span enough of the two sequences, skip it.
                            //
                            if (aln_res.pct_id >= min_match_len) {

                                if (d <= ctag_dist)
                                    tag_1->add_match(tag_2->id, cat_hit.first, allele->first, d, invert_cigar(aln_res.cigar));
                            }
                        }
                    }
                }
            }

            // time_4 = clock();
            // per_locus += (time_4 - time_3);
        }

        //
        // Free the k-mers we generated for this query
        //
        for (uint j = 0; j < kmers.size(); j++)
            delete [] kmers[j];
        kmers.clear();

        delete aln;
    }

    // cerr << "Time to kmerize catalog: " << time_2 - time_1 << "\n"
    //      << "Average time per locus:  " << per_locus / (double) keys.size() << "\n";

    free_kmer_hash(kmer_map, kmer_map_keys);

    return 0;
}

bool compare_matches(Match *a, Match *b) {
    return (a->dist < b->dist);
}

int find_matches_by_sequence(map<int, CLocus *> &catalog, map<int, QLocus *> &sample) {
    //
    // Calculate the distance (number of mismatches) between each pair
    // of Radtags. We expect all radtags to be the same length;
    //
    map<int, QLocus *>::iterator i;
    map<int, CLocus *>::iterator j;
    int k;

    // OpenMP can't parallelize random access iterators, so we convert
    // our map to a vector of integer keys.
    vector<int> keys;
    for (i = sample.begin(); i != sample.end(); i++)
        keys.push_back(i->first);

    #pragma omp parallel private(i, j, k)
    {
        #pragma omp for schedule(dynamic)
        for (k = 0; k < (int) keys.size(); k++) {

            i = sample.find(keys[k]);

            vector<pair<allele_type, string> >::iterator r, s;

            //
            // Iterate through the possible SAMPLE alleles
            //
            for (r = i->second->strings.begin(); r != i->second->strings.end(); r++) {

                for (j = catalog.begin(); j != catalog.end(); j++) {
                    //
                    // Iterate through the possible CATALOG alleles
                    //
                    for (s = j->second->strings.begin(); s != j->second->strings.end(); s++) {
                        if (r->second == s->second) {
                            //cerr << "Found a match between " << i->first << " (" << r->first << ") and " << j->first << " (" << s->first << ")\n";

                            i->second->add_match(j->second->id, s->first, r->first, 0);
                        }
                    }
                }
            }
        }
    }

    return 0;
}

int find_matches_by_genomic_loc(map<string, int> &cat_index, map<int, QLocus *> &sample) {
    map<int, QLocus *>::iterator i;
    map<int, CLocus *>::iterator j;

    //
    // OpenMP can't parallelize random access iterators, so we convert
    // our map to a vector of integer keys.
    //
    vector<int> keys;
    for (i = sample.begin(); i != sample.end(); i++)
        keys.push_back(i->first);

    #pragma omp parallel private(i, j)
    {
        char id[id_len];

        #pragma omp for
        for (int k = 0; k < (int) keys.size(); k++) {

            i = sample.find(keys[k]);

            snprintf(id, id_len - 1, "%s|%d|%c",
                     i->second->loc.chr(),
                     i->second->loc.bp,
                     i->second->loc.strand == strand_plus ? '+' : '-');

            if (cat_index.count(id) > 0)
                i->second->add_match(cat_index[id], "", "", 0);
        }
    }

    return 0;
}

int write_catalog(map<int, CLocus *> &catalog) {

    map<int, CLocus *>::iterator i;
    CLocus  *tag;
    set<int> matches;

    bool gzip = (in_file_type == FileT::gzsql) ? true : false;

    string tag_file = out_path + "catalog.tags.tsv";
    string snp_file = out_path + "catalog.snps.tsv";
    string all_file = out_path + "catalog.alleles.tsv";

    if (gzip) {
        tag_file += ".gz";
        snp_file += ".gz";
        all_file += ".gz";
    }

    //
    // Open the output files for writing.
    //
    gzFile   gz_tags, gz_snps, gz_alle;
    ofstream tags, snps, alle;
    if (gzip) {
        gz_tags = gzopen(tag_file.c_str(), "wb");
        if (!gz_tags) {
            cerr << "Error: Unable to open gzipped catalog tag file '" << tag_file << "': " << strerror(errno) << ".\n";
            exit(1);
        }
        #if ZLIB_VERNUM >= 0x1240
        gzbuffer(gz_tags, libz_buffer_size);
        #endif
        gz_snps = gzopen(snp_file.c_str(), "wb");
        if (!gz_snps) {
            cerr << "Error: Unable to open gzipped catalog snps file '" << snp_file << "': " << strerror(errno) << ".\n";
            exit(1);
        }
        #if ZLIB_VERNUM >= 0x1240
        gzbuffer(gz_snps, libz_buffer_size);
        #endif
        gz_alle = gzopen(all_file.c_str(), "wb");
        if (!gz_alle) {
            cerr << "Error: Unable to open gzipped catalog alleles file '" << all_file << "': " << strerror(errno) << ".\n";
            exit(1);
        }
        #if ZLIB_VERNUM >= 0x1240
        gzbuffer(gz_alle, libz_buffer_size);
        #endif
    } else {
        tags.open(tag_file.c_str());
        if (tags.fail()) {
            cerr << "Error: Unable to open catalog tag file for writing.\n";
            exit(1);
        }
        snps.open(snp_file.c_str());
        if (snps.fail()) {
            cerr << "Error: Unable to open catalog SNPs file for writing.\n";
            exit(1);
        }
        alle.open(all_file.c_str());
        if (alle.fail()) {
            cerr << "Error: Unable to open catalog alleles file for writing.\n";
            exit(1);
        }
    }

    //
    // Record the version of Stacks used and the date generated as a comment in the catalog.
    //
    // Obtain the current date.
    //
    stringstream log;
    time_t       rawtime;
    struct tm   *timeinfo;
    char         date[32];
    time(&rawtime);
    timeinfo = localtime(&rawtime);
    strftime(date, 32, "%F %T", timeinfo);
    log << "# cstacks version " << VERSION << "; catalog generated on " << date << "\n";
    if (gzip) {
        gzputs(gz_tags, log.str().c_str());
        gzputs(gz_snps, log.str().c_str());
        gzputs(gz_alle, log.str().c_str());
    } else {
        tags << log.str();
        snps << log.str();
        alle << log.str();
    }

    for (i = catalog.begin(); i != catalog.end(); i++) {
        tag = i->second;

        if (gzip)
            write_gzip_output(tag, gz_tags, gz_snps, gz_alle);
        else
            write_simple_output(tag, tags, snps, alle);
    }

    if (gzip) {
        gzclose(gz_tags);
        gzclose(gz_snps);
        gzclose(gz_alle);
    } else {
        tags.close();
        snps.close();
        alle.close();
    }

    return 0;
}

int merge_allele(Locus *locus, SNP *snp) {
    map<int, pair<string, SNP *> > columns;
    map<int, pair<string, SNP *> >::iterator c;
    vector<SNP *>::iterator i;
    SNP *lsnp;

    for (i = locus->snps.begin(); i != locus->snps.end(); i++)
        columns[(*i)->col] = make_pair("sample", *i);

    if (columns.count(snp->col)) {
        lsnp = columns[snp->col].second;

        //
        // If this is a new allele for this nucleotide, add it to the catalog SNP.
        //
        bool rank_1_exists = false;
        bool rank_2_exists = false;

        if (snp->rank_1 == lsnp->rank_1 ||
            snp->rank_1 == lsnp->rank_2 ||
            snp->rank_1 == lsnp->rank_3 ||
            snp->rank_1 == lsnp->rank_4) {
            rank_1_exists = true;
        }
        if (snp->rank_2 == lsnp->rank_1 ||
            snp->rank_2 == lsnp->rank_2 ||
            snp->rank_2 == lsnp->rank_3 ||
            snp->rank_2 == lsnp->rank_4) {
            rank_2_exists = true;
        }

        if (rank_1_exists == false) {
            if (lsnp->rank_3 == 0)
                lsnp->rank_3 = snp->rank_1;
            else
                lsnp->rank_4 = snp->rank_1;
        }
        if (rank_2_exists == false) {
            if (lsnp->rank_3 == 0)
                lsnp->rank_3 = snp->rank_2;
            else
                lsnp->rank_4 = snp->rank_2;
        }

        columns[snp->col] = make_pair("both", lsnp);
    } else {
        columns[snp->col] = make_pair("merge", snp);
    }

    vector<pair<string, SNP *> > merged_snps;

    for (c = columns.begin(); c != columns.end(); c++)
        merged_snps.push_back((*c).second);

    //
    // Sort the SNPs by column
    //
    sort(merged_snps.begin(), merged_snps.end(), compare_pair_snp);

    //
    // Modify any existing alleles to account for this new SNP. If there are not any alleles,
    // create new ones.
    //
    stringstream sallele;
    set<string> merged_alleles;
    string allele, new_allele;
    int pos;

    if (locus->alleles.size() == 0) {
        sallele << locus->con[snp->col];
        merged_alleles.insert(sallele.str());
    }

    map<string, int>::iterator j;
    vector<pair<string, SNP *> >::iterator k;

    for (j = locus->alleles.begin(); j != locus->alleles.end(); j++) {
        allele     = j->first;
        new_allele = "";
        pos        = 0;

        // cerr << "Allele length: " << allele.size() << "\n";

        for (k = merged_snps.begin(); k != merged_snps.end(); k++) {
            //
            // If we inserted a SNP from the sample, add the proper nucleotide from the consensus
            // sequence to account for it in the allele string.
            //
            if ((*k).first == "merge") {
                new_allele += locus->con[(*k).second->col];
                // cerr << "  Adding char '" << locus->con[k->second->col] << "' from consensus position " << (*k).second->col << "\n";
            } else {
                new_allele += allele[pos];
                // cerr << "  Adding char '" << allele[pos] << "' from allele position " << pos << "\n";
                pos++;
            }
        }

        merged_alleles.insert(new_allele);
    }

    set<string>::iterator s;

    locus->alleles.clear();
    for (s = merged_alleles.begin(); s != merged_alleles.end(); s++) {
        locus->alleles[*s] = 0;
    }

    return 1;
}

int CLocus::merge_snps(QLocus *matched_tag) {
    vector<SNP *>::iterator i;
    map<string, int>::iterator j;
    vector<pair<string, SNP *> >::iterator k;
    map<int, pair<string, SNP *> > columns;
    map<int, pair<string, SNP *> >::iterator c;

    vector<pair<string, SNP *> > merged_snps;
    set<string> merged_alleles;
    set<string>::iterator s;
    SNP *csnp;

    for (i = this->snps.begin(); i != this->snps.end(); i++)
        columns[(*i)->col] = make_pair("catalog", *i);

    for (i = matched_tag->snps.begin(); i != matched_tag->snps.end(); i++) {
        //
        // Is this column already represented from the previous sample?
        //
        if (columns.count((*i)->col)) {
            csnp = columns[(*i)->col].second;

            //
            // If this is a new allele for this nucleotide, add it to the catalog SNP.
            //
            bool rank_1_exists = false;
            bool rank_2_exists = false;

            if ((*i)->rank_1 == csnp->rank_1 ||
                (*i)->rank_1 == csnp->rank_2 ||
                (*i)->rank_1 == csnp->rank_3 ||
                (*i)->rank_1 == csnp->rank_4) {
                rank_1_exists = true;
            }
            if ((*i)->rank_2 == csnp->rank_1 ||
                (*i)->rank_2 == csnp->rank_2 ||
                (*i)->rank_2 == csnp->rank_3 ||
                (*i)->rank_2 == csnp->rank_4) {
                rank_2_exists = true;
            }

            if (rank_1_exists == false) {
                if (csnp->rank_3 == 0)
                    csnp->rank_3 = (*i)->rank_1;
                else
                    csnp->rank_4 = (*i)->rank_1;
            }
            if (rank_2_exists == false) {
                if (csnp->rank_3 == 0)
                    csnp->rank_3 = (*i)->rank_2;
                else
                    csnp->rank_4 = (*i)->rank_2;
            }

            columns[(*i)->col] = make_pair("both", csnp);
        } else {
            columns[(*i)->col] = make_pair("sample", *i);
        }
    }

    for (c = columns.begin(); c != columns.end(); c++)
        merged_snps.push_back((*c).second);

    //
    // Sort the SNPs by column
    //
    sort(merged_snps.begin(), merged_snps.end(), compare_pair_snp);

    //
    // If the catalog tag has no defined alleles, create a matching haplotype
    // from the consensus sequence before merging in the new alleles.
    //
    string allele, new_allele;
    int    pos;

    if (this->alleles.size() == 0) {
        char c;
        new_allele = "";
        for (k = merged_snps.begin(); k != merged_snps.end(); k++) {
            csnp = k->second;
            c    = this->con[k->second->col];

            new_allele += (csnp->col > this->len - 1) ? 'N' : c;

            if (csnp->col > this->len - 1) continue;

            if (c != csnp->rank_1 &&
                c != csnp->rank_2 &&
                c != csnp->rank_3 &&
                c != csnp->rank_4) {

                if (csnp->rank_3 == 0)
                    csnp->rank_3 = c;
                else
                    csnp->rank_4 = c;
            }
        }

        if (new_allele.length() > 0)
            merged_alleles.insert(new_allele);
    }

    //
    // Merge the alleles accounting for any SNPs added from either of the two samples.
    //
    for (j = this->alleles.begin(); j != this->alleles.end(); j++) {
        allele     = j->first;
        new_allele = "";
        pos        = 0;

        for (k = merged_snps.begin(); k != merged_snps.end(); k++) {
            //
            // If we inserted a SNP from the sample, add the proper nucleotide from the consensus
            // sequence to account for it in the allele string.
            //
            if (k->first == "sample") {
                new_allele += k->second->col > this->len - 1 ? 'N' : this->con[k->second->col];
            } else {
                new_allele += allele[pos];
                pos++;
            }
        }

        merged_alleles.insert(new_allele);
    }

    for (j = matched_tag->alleles.begin(); j != matched_tag->alleles.end(); j++) {
        allele     = j->first;
        new_allele = "";
        pos        = 0;

        for (k = merged_snps.begin(); k != merged_snps.end(); k++) {
            if (k->first == "catalog") {
                new_allele += k->second->col > matched_tag->len - 1 ? 'N' : matched_tag->con[k->second->col];
            } else {
                new_allele += allele[pos];
                pos++;
            }
        }

        merged_alleles.insert(new_allele);
    }

    //
    // If the matching tag being merged into the catalog had no called SNPs
    // create alleles from the consensus sequence and check that catalog SNP
    // objects contain all the nucleoties.
    //
    if (matched_tag->alleles.size() == 0) {
        char c;
        new_allele = "";
        for (k = merged_snps.begin(); k != merged_snps.end(); k++) {
            csnp = k->second;
            c    = matched_tag->con[k->second->col];

            new_allele += (csnp->col > matched_tag->len - 1) ? 'N' : c;

            if (csnp->col > matched_tag->len - 1) continue;

            if (c != csnp->rank_1 &&
                c != csnp->rank_2 &&
                c != csnp->rank_3 &&
                c != csnp->rank_4) {

                if (csnp->rank_3 == 0)
                    csnp->rank_3 = c;
                else
                    csnp->rank_4 = c;
            }
        }

        if (new_allele.length() > 0)
            merged_alleles.insert(new_allele);
    }

    // //
    // // If the newly merged alleles contain Ns due to different sequence lengths,
    // // check if we can reduce the alleles as one of the longer allele haplotypes
    // // may fully encompass a shorter allele haplotype that has been padded with Ns.
    // //
    // if (require_uniq_haplotypes) this->reduce_alleles(merged_alleles);

    //
    // Update the catalog entry's list of SNPs and alleles
    //
    this->snps.clear();

    for (k = merged_snps.begin(); k != merged_snps.end(); k++) {
        SNP *snp    = new SNP;
        snp->col    = (*k).second->col;
        snp->type   = (*k).second->type;
        snp->lratio = 0.0;
        snp->rank_1 = (*k).second->rank_1;
        snp->rank_2 = (*k).second->rank_2;
        snp->rank_3 = (*k).second->rank_3;
        snp->rank_4 = (*k).second->rank_4;

        this->snps.push_back(snp);

        if (k->first == "catalog" || k->first == "both")
            delete k->second;
    }

    this->alleles.clear();
    for (s = merged_alleles.begin(); s != merged_alleles.end(); s++) {
        this->alleles[*s] = 0;
    }

    return 1;
}

int
CLocus::reduce_alleles(set<string> &alleles)
{
    set<string>::iterator it;
    uint len, max_len, match, ncnt;
    vector<string> haplotypes, cur, next;

    max_len = 0;
    for (it = alleles.begin(); it != alleles.end(); it++) {
        max_len = it->length() > max_len ? it->length() : max_len;
        haplotypes.push_back(*it);
    }

    len = alleles.size();
    alleles.clear();

    for (uint i = 0; i < len; i++) {
        //cerr << "Looking at haplotype[" << i << "]: " << haplotypes[i] << "\n";
        //
        // We will only look at strings that contain Ns.
        //
        if (haplotypes[i].find('N') == string::npos) {
            alleles.insert(haplotypes[i]);
            //cerr << "  No Ns, skipping...\n";
            continue;
        }

        uint k = 0;
        uint j = i + 1;
        while (k < len - 1) {
            cur.push_back(haplotypes[j % len]);
            k++;
            j++;
        }

        //
        // Examine the haplotype alleles one SNP at a time. If we are able to uniquely
        // determine a second haplotype that encompasses the first
        // to, return it.
        //
        j = 0;
        while (cur.size() > 1 && j < max_len) {

            for (k = 0; k < cur.size(); k++) {
                cerr << "Comparing haplotypes[" << i << "]: '" << haplotypes[i] << "' to '" << cur[k] << " at position " << j << "'\n";
                if (haplotypes[i][j] == cur[k][j] || haplotypes[i][j] == 'N') {
                    cerr << "  Keeping this haplotype.\n";
                    next.push_back(cur[k]);
                } else {
                    cerr << "  Discarding this haplotype.\n";
                }
            }
            cur = next;
            next.clear();
            j++;
        }

        //
        // If there is only one left, make sure what we have of the haplotype does match
        // and its not simply an erroneously called haplotype. If so, then this haplotype
        // is encompassed by another, longer haplotype and we do not need to keep it.
        //
        ncnt  = 0;
        match = 0;
        if (cur.size() > 1) {
            cerr << "Discarding " << haplotypes[i] << "\n";
            continue;
        } else if (cur.size() == 1) {
            for (k = 0; k < max_len; k++)
                if (haplotypes[i][k] != 'N') ncnt++;
            for (k = 0; k < max_len; k++)
                if (cur[0][k] == haplotypes[i][k]) match++;
            if (match == ncnt) {
                cerr << "Discarding " << haplotypes[i] << "\n";
                continue;
            }
        }

        cerr << "Keeping " << haplotypes[i] << "\n";
        alleles.insert(haplotypes[i]);
    }

    return 0;
}

int
populate_kmer_hash(map<int, CLocus *> &catalog, CatKmerHashMap &kmer_map, vector<char *> &kmer_map_keys, int kmer_len)
{
    map<int, CLocus *>::iterator it;
    vector<pair<allele_type, string> >::iterator allele;
    vector<char *> kmers;
    CLocus        *tag;
    char          *hash_key;
    bool           exists;
    int            j;

    //
    // Break each stack down into k-mers and create a hash map of those k-mers
    // recording in which sequences they occur.
    //
    int num_kmers = strlen(catalog.begin()->second->con) - kmer_len + 1;

    for (it = catalog.begin(); it != catalog.end(); it++) {
        tag = it->second;

        //
        // Iterate through the possible Catalog alleles
        //
        for (allele = tag->strings.begin(); allele != tag->strings.end(); allele++) {
            //
            // Generate and hash the kmers for this allele string
            //
            generate_kmers(allele->second.c_str(), kmer_len, num_kmers, kmers);

            for (j = 0; j < num_kmers; j++) {
                hash_key = kmers[j];
                exists   = kmer_map.count(hash_key) == 0 ? false : true;

                kmer_map[hash_key].push_back(make_pair(allele->first, tag->id));

                if (exists)
                    delete [] kmers[j];
                else
                    kmer_map_keys.push_back(hash_key);
            }

            kmers.clear();
        }
    }

    //dump_kmer_map(kmer_map);

    return 0;
}

int
write_simple_output(CLocus *tag, ofstream &cat_file, ofstream &snp_file, ofstream &all_file)
{
    vector<SNP *>::iterator           snp_it;
    map<string, int>::iterator        all_it;
    vector<pair<int, int> >::iterator src_it;
    string sources;

    for (src_it = tag->sources.begin(); src_it != tag->sources.end(); src_it++) {
        stringstream s;
        s << (*src_it).first << "_" << (*src_it).second << ",";
        sources += s.str();
    }
    sources = sources.substr(0, sources.length() - 1);

    cat_file <<
        tag->id      << "\t" <<
        tag->loc.chr() << "\t" <<
        tag->loc.bp  << "\t" <<
        (tag->loc.strand == strand_plus ? "+" : "-") << "\t" <<
        "consensus"  << "\t" <<
        "0"          << "\t" <<
        sources      << "\t" <<
        tag->con     << "\t" <<
        0            << "\t" <<  // These flags are unused in cstacks, but important in ustacks
        0            << "\t" <<
        0            << "\t" <<
        0            << "\n";

    //
    // Output the SNPs associated with the catalog tag
    //
    for (snp_it = tag->snps.begin(); snp_it != tag->snps.end(); snp_it++) {
        snp_file <<   tag->id        << "\t"
                 << (*snp_it)->col << "\t";

        switch((*snp_it)->type) {
        case snp_type_het:
            snp_file << "E\t";
            break;
        case snp_type_hom:
            snp_file << "O\t";
            break;
        default:
            snp_file << "U\t";
            break;
        }

        snp_file <<
            (*snp_it)->lratio << "\t" <<
            (*snp_it)->rank_1 << "\t" <<
            (*snp_it)->rank_2 << "\t" <<
            ((*snp_it)->rank_3 == 0 ? '-' : (*snp_it)->rank_3) << "\t" <<
            ((*snp_it)->rank_4 == 0 ? '-' : (*snp_it)->rank_4) << "\n";
    }

    //
    // Output the alleles associated with the two matched tags
    //
    for (all_it = tag->alleles.begin(); all_it != tag->alleles.end(); all_it++)
        all_file <<
            tag->id       << "\t" <<
            all_it->first << "\t" <<
            "0"           << "\t" <<    // These two fields are used in the
            "0"           << "\n";      // ustacks/pstacks output, not in cstacks.

    return 0;
}

int
write_gzip_output(CLocus *tag, gzFile &cat_file, gzFile &snp_file, gzFile &all_file)
{
    vector<SNP *>::iterator           snp_it;
    map<string, int>::iterator        all_it;
    vector<pair<int, int> >::iterator src_it;
    string       sources;
    stringstream sstr;

    for (src_it = tag->sources.begin(); src_it != tag->sources.end(); src_it++) {
        sstr << (*src_it).first << "_" << (*src_it).second << ",";
    }
    sources = sstr.str();
    sources = sources.substr(0, sources.length() - 1);

    sstr.str("");

    sstr << "0"         << "\t"
         << tag->id     << "\t"
         << "consensus" << "\t"
         << "0"         << "\t"
         << sources     << "\t"
         << tag->con    << "\t"
         << 0           << "\t" // These flags are unused in cstacks, but important in ustacks
         << 0           << "\t"
         << 0           << "\n";

    gzputs(cat_file, sstr.str().c_str());
    sstr.str("");

    //
    // Output the SNPs associated with the catalog tag
    //
    for (snp_it = tag->snps.begin(); snp_it != tag->snps.end(); snp_it++) {
        sstr << "0"            << "\t"
             <<   tag->id      << "\t"
             << (*snp_it)->col << "\t";

        switch((*snp_it)->type) {
        case snp_type_het:
            sstr << "E\t";
            break;
        case snp_type_hom:
            sstr << "O\t";
            break;
        default:
            sstr << "U\t";
            break;
        }

        sstr <<
            (*snp_it)->lratio << "\t" <<
            (*snp_it)->rank_1 << "\t" <<
            (*snp_it)->rank_2 << "\t" <<
            ((*snp_it)->rank_3 == 0 ? '-' : (*snp_it)->rank_3) << "\t" <<
            ((*snp_it)->rank_4 == 0 ? '-' : (*snp_it)->rank_4) << "\n";
    }

    gzputs(snp_file, sstr.str().c_str());
    sstr.str("");

    //
    // Output the alleles associated with the two matched tags
    //
    for (all_it = tag->alleles.begin(); all_it != tag->alleles.end(); all_it++)
        sstr << "0"     << "\t"
             << tag->id << "\t"
             << all_it->first << "\t"
             << 0       << "\t"
             << 0       << "\n";

    gzputs(all_file, sstr.str().c_str());

    return 0;
}

int
initialize_new_catalog(pair<int, string> &sample, map<int, CLocus *> &catalog)
{
    map<int, CLocus *> tmp_catalog;
    bool compressed = false;

    //
    // Parse the input files.
    //
    if (!load_loci(sample.second, tmp_catalog, 0, false, compressed))
        return 0;

    in_file_type = compressed == true ? FileT::gzsql : FileT::sql;

    sample.first = tmp_catalog.begin()->second->sample_id;

    //
    // Iterate over the catalog entires and renumber them after recording the source of
    // locus.
    //
    map<int, CLocus *>::iterator j;
    int k = 1;
    for (j = tmp_catalog.begin(); j != tmp_catalog.end(); j++) {
        j->second->sources.push_back(make_pair(sample.first, j->second->id));
        j->second->id = k;

        catalog[k] = j->second;

        k++;
    }

    cerr << "  " << catalog.size() << " loci were newly added to the catalog.\n";

    return 1;
}

int
initialize_existing_catalog(string catalog_path, map<int, CLocus *> &catalog)
{
    bool compressed;

    //
    // Parse the input files.
    //
    catalog_path += ".catalog";
    if (!load_loci(catalog_path, catalog, 0, false, compressed))
        return 0;

    in_file_type = compressed == true ? FileT::gzsql : FileT::sql;

    //
    // Iterate over the catalog entires and convert the stack components
    // into source objects, to record what samples each locus came from.
    //
    map<int, CLocus *>::iterator j;
    CLocus *loc;
    char   *p, *q;
    int     sample_id, locus_id;

    for (j = catalog.begin(); j != catalog.end(); j++) {
        loc = j->second;

        for (uint i = 0; i < loc->comp.size(); i++) {
            //
            // Parse the ID into sample ID / locus ID, given 43_1356, parse into
            // sample ID 43 and locus ID 1356.
            //
            for (p = loc->comp[i]; *p != '_' && *p != '\0'; p++);
            if (*p != '_')
                return 0;
            p++;
            sample_id = strtol(loc->comp[i], &q, 10);
            if (*q != '_')
                return 0;

            locus_id = strtol(p, &q, 10);

            if (*q != '\0')
                return 0;

            loc->sources.push_back(make_pair(sample_id, locus_id));
        }
    }

    return 1;
}

int parse_command_line(int argc, char* argv[]) {
    string in_dir;
    string popmap_path;

    while (1) {
        static struct option long_options[] = {
            {"help",            no_argument, NULL, 'h'},
            {"version",         no_argument, NULL, 1000},
            {"mmatches",        no_argument, NULL, 'm'},
            {"aligned",         no_argument, NULL, 'g'},
            {"uniq_haplotypes", no_argument, NULL, 'u'},
            {"report_mmatches", no_argument, NULL, 'R'},
            {"disable_gapped",  no_argument, NULL, 'G'},
            {"max_gaps",        required_argument, NULL, 'X'},
            {"min_aln_len",     required_argument, NULL, 'x'},
            {"ctag_dist",       required_argument, NULL, 'n'},
            {"k_len",           required_argument, NULL, 'k'},
            {"catalog",         required_argument, NULL, 'c'},
            {"sample",          required_argument, NULL, 's'},
            {"outpath",         required_argument, NULL, 'o'},
            {"num_threads",     required_argument, NULL, 'p'},
            {0, 0, 0, 0}
        };

        int option_index = 0;
        int c = getopt_long(argc, argv, "hgvuRmGX:x:o:s:c:p:n:k:P:M:", long_options, &option_index);

        // Detect the end of the options.
        if (c == -1)
            break;

        switch (c) {
        case 'h':
            help();
            break;
        case 'n':
            ctag_dist = is_integer(optarg);
            break;
        case 'k':
            set_kmer_len = false;
            kmer_len     = is_integer(optarg);
            break;
        case 'm':
            mult_matches = true;
            break;
        case 'R':
            report_mmatches = true;
            break;
        case 'g':
            search_type = genomic_loc;
            break;
        case 's':
            samples.push(make_pair(0, optarg));
            break;
        case 'c':
            catalog_path = optarg;
            break;
        case 'o':
            out_path = optarg;
            break;
        case 'u':
            require_uniq_haplotypes = true;
            break;
        case 'G':
            gapped_alignments = false;
            break;
        case 'X':
            max_gaps = is_double(optarg);
            break;
        case 'x':
            min_match_len = is_double(optarg);
            break;
        case 1000:
            version();
            break;
        case 'p':
            num_threads = is_integer(optarg);
            break;
        case 'P':
            in_dir = optarg;
            break;
        case 'M':
            popmap_path = optarg;
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

    if (optind < argc) {
        cerr << "Error: Failed to parse command line: '" << argv[optind] << "' is seen as a positional argument. Expected no positional arguments.\n";
        help();
    }

    if (in_dir.empty() && samples.empty()) {
        cerr << "Error: You must specify one of -P or -s.\n";
        help();
    } else if ((!in_dir.empty() || !popmap_path.empty())
            && (!samples.empty() || !out_path.empty())) {
        cerr << "Error: Please use options -P/-M or -s/-o, not both.\n";
        help();
    }

    if (!in_dir.empty()) {
        if (popmap_path.empty()) {
            cerr << "Error: Please specify a population map (-M).\n";
            help();
        }

        if (in_dir.back() != '/')
            in_dir += "/";

        // Set `samples`.
        MetaPopInfo popmap;
        popmap.init_popmap(popmap_path);
        for (const Sample& s : popmap.samples())
            samples.push({0, in_dir + s.name});

        // Set `out_path`.
        out_path = in_dir;

    } else if (!samples.empty()) {
        if (out_path.empty())
            out_path = ".";

        if (out_path.back() != '/')
            out_path += "/";
    }

    if (set_kmer_len == false && (kmer_len < 5 || kmer_len > 31)) {
        cerr << "Kmer length must be between 5 and 31bp.\n";
        help();
    }

    return 0;
}

void version() {
    cerr << "cstacks " << VERSION << "\n";

    exit(1);
}

void help() {
    cerr << "cstacks " << VERSION << "\n"
              << "cstacks -P in_dir -M popmap [-n num_mismatches] [--gapped] [-p num_threads]" << "\n"
              << "cstacks --aligned -P in_dir -M popmap [-p num_threads]" << "\n"
              << "cstacks -s sample1_path [-s sample2_path ...] -o path [-n num_mismatches] [--gapped] [-p num_threads]" << "\n"
              << "cstacks --aligned -s sample1_path [-s sample2_path ...] -o path [-p num_threads]" << "\n"
              << "\n"
              << "  P: path to the directory containing Stacks files.\n"
              << "  M: path to a population map file.\n"
              << "  n: number of mismatches allowed between sample loci when build the catalog (default 1)." << "\n"
              << "  p: enable parallel execution with num_threads threads.\n"
              << "  s: sample prefix from which to load loci into the catalog." << "\n"
              << "  o: output path to write results." << "\n"
              << "  --catalog <path>: add to an existing catalog.\n"
              << "\n"
              << "Gapped assembly options:\n"
              << "  --max_gaps: number of gaps allowed between stacks before merging (default: 2).\n"
              << "  --min_aln_len: minimum length of aligned sequence in a gapped alignment (default: 0.80).\n"
              << "  --disable_gapped: disable gapped alignments between stacks (default: use gapped alignments).\n"
              << "\n"
              << "Advanced options:\n"
              << "  m: include tags in the catalog that match to more than one entry (default false)." << "\n"
              << "  --k_len <len>: specify k-mer size for matching between between catalog loci (automatically calculated by default).\n"
              << "  --report_mmatches: report query loci that match more than one catalog locus.\n";

    exit(1);
}
