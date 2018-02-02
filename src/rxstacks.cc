// -*-mode:c++; c-style:k&r; c-basic-offset:4;-*-
//
// Copyright 2013-2016, Julian Catchen <jcatchen@illinois.edu>
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
// rxstacks -- make model call corrections and haplotype corrections
// across a population of samples.
//

#include "MetaPopInfo.h"
#include "catalog_utils.h"

#include "rxstacks.h"

// Global variables to hold command-line options.
int    num_threads = 1;
int    batch_id    = -1;
string in_path;
string out_path;
FileT  in_file_type      = FileT::sql;
double confounded_limit  = 0.75;
bool   filter_confounded = false;
bool   prune_haplotypes  = false;
int    max_haplotype_cnt = 0;
bool   lnl_dist          = false;
bool   filter_lnl        = false;
double lnl_limit         = 0.0;
bool   verbose           = false;

//
// For use with the multinomial model to call fixed nucleotides.
//
modelt model_type         = snp;
double alpha              = 0.1;

int main (int argc, char* argv[]) {
    IF_NDEBUG_TRY

    parse_command_line(argc, argv);

    cerr
        << "Log liklihood filtering: " << (filter_lnl        == true ? "on"  : "off") << "; threshold: " << lnl_limit << "\n"
        << "Prune haplotypes: "        << (prune_haplotypes  == true ? "yes" : "no")  << "\n"
        << "Filter confounded loci: "  << (filter_confounded == true ? "yes" : "no")  << "\n";

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
    if (verbose) num_threads = 1;
    #ifdef _OPENMP
    omp_set_num_threads(num_threads);
    #endif

    MetaPopInfo mpopi;
    mpopi.init_directory(in_path);

    //
    // Open and initialize the log files.
    //
    ofstream log_fh, log_snp_fh, log_hap_fh;
    init_log(argc, argv, log_fh, log_snp_fh, log_hap_fh);

    //
    // Load the catalog
    //
    stringstream catalog_file;
    map<int, CSLocus *> catalog;
    bool compressed = false;
    catalog_file << in_path << "batch_" << batch_id << ".catalog";
    int res = load_loci(catalog_file.str(), catalog, 0, false, compressed);
    if (res  == 0) {
        cerr << "Error: Unable to load the catalog '" << catalog_file.str() << "'\n";
        throw exception();
    }

    in_file_type = compressed == true ? FileT::gzsql : FileT::sql;

    //
    // Let's fill in the SNP model calls to include both hets and homozygotes to
    // make it easier to iterate over them later.
    //
    fill_catalog_snps(catalog);

    //
    // Load matches to the catalog
    //
    vector<vector<CatMatch *> > catalog_matches;
    vector<size_t> samples_to_remove;
    set<size_t> seen_samples;
    for (size_t i=0; i<mpopi.samples().size(); ++i) {
        catalog_matches.push_back(vector<CatMatch*>());
        vector<CatMatch *>& m = catalog_matches.back();

        const Sample& sample = mpopi.samples()[i];
        load_catalog_matches(in_path + sample.name, m);

        if (m.size() == 0) {
            cerr << "Warning: Absent or malformed matches file '"
                 << sample.name << ".matches.tsv(.gz)"
                 <<"', excluding this sample from population analysis.\n";
            samples_to_remove.push_back(i);
            catalog_matches.pop_back(); // n.b. Index shift will be resolved by the call to MetaPopInfo::delete_samples().
            continue;
        }

        size_t sample_id = m[0]->sample_id;
        if (seen_samples.count(sample_id) > 0) {
            cerr << "Error: sample ID " << sample_id << " occurs twice in this data set, likely the pipeline was run incorrectly.\n";
            throw exception();
        }
        seen_samples.insert(sample_id);
        mpopi.set_sample_id(i, sample_id);
    }
    mpopi.delete_samples(samples_to_remove);
    if (mpopi.samples().size() == 0) {
        cerr << "Error: Couln't find any matches files.\n";
        throw exception();
    }
    // [mpopi] is definitive.
    cerr << "Working on " << mpopi.samples().size() << " samples.\n";

    //
    // Create the population map
    //
    cerr << "Populating observed haplotypes for " << mpopi.samples().size() << " samples, " << catalog.size() << " loci.\n";
    PopMap<CSLocus> *pmap = new PopMap<CSLocus>(mpopi, catalog.size());
    pmap->populate(catalog, catalog_matches);

    //
    // Sum haplotype counts across the population for each catalog locus.
    //
    sum_haplotype_counts(catalog, pmap);

    //
    // Calculate mean log likelihood across the population for each catalog locus.
    //
    calc_lnl_means(catalog, pmap);

    int    catalog_id, tag_id;

    //
    // Process samples matched to the catalog, one by one.
    //
    for (uint i = 0; i < mpopi.samples().size(); i++) {
        const Sample& sample = mpopi.samples()[i];

        cerr << "Loading stacks from sample " << sample.name << " [" << i+1 << " of " << mpopi.samples().size() << "]...\n";

        //////
        //////
        ////// if (sample_id != 176) continue;
        //////
        //////

        map<int, Locus *> stacks;
        int res;
        if ((res = load_loci(in_path + sample.name, stacks, 2, true, compressed)) == 0) {
            cerr << "Unable to load sample file '" << sample.name << "'\n";
            continue;
        }

        cerr << "Making corrections to sample " << sample.name << "...";

        set<pair<int, int> >           uniq_matches;
        set<pair<int, int> >::iterator it;
        vector<pair<int, int> >        matches;

        //
        // There are multiple matches per stack, but we only need to process
        // each stack once to make corrections.
        //
        for (uint j = 0; j < catalog_matches[i].size(); j++) {
            catalog_id = catalog_matches[i][j]->cat_id;
            tag_id     = catalog_matches[i][j]->tag_id;

            uniq_matches.insert(make_pair(catalog_id, tag_id));
        }

        //
        // Put the catalog/tag ID pairs into a vector for parallel processing.
        //
        for (it = uniq_matches.begin(); it != uniq_matches.end(); it++)
            matches.push_back(*it);

        unsigned long int nuc_cnt        = 0;
        unsigned long int unk_hom_cnt    = 0;
        unsigned long int unk_het_cnt    = 0;
        unsigned long int het_unk_cnt    = 0;
        unsigned long int hom_unk_cnt    = 0;
        unsigned long int het_hom_cnt    = 0;
        unsigned long int hom_het_cnt    = 0;
        unsigned long int conf_loci_cnt  = 0;
        unsigned long int pruned_hap_cnt = 0;
        unsigned long int pruned_mst_hap_cnt = 0;
        unsigned long int blacklist_cnt  = 0;
        unsigned long int lnl_cnt        = 0;

        #pragma omp parallel private(catalog_id, tag_id)
        {
            Datum   *d;
            Locus   *loc;
            CSLocus *cloc;
            uint     seq_len;
            string   seq, model;
            char    *adj_seq;
            vector<pair<char, uint> > cigar;
            vector<char *>            reads;

            #pragma omp for schedule(dynamic, 1) reduction(+:nuc_cnt) reduction(+:unk_hom_cnt) reduction(+:unk_het_cnt) \
                reduction(+:hom_unk_cnt) reduction(+:het_unk_cnt) reduction(+:hom_het_cnt) reduction(+:het_hom_cnt) \
                reduction(+:conf_loci_cnt) reduction(+:pruned_hap_cnt) reduction(+:pruned_mst_hap_cnt) reduction(+:blacklist_cnt) reduction(+:lnl_cnt)
            for (uint j = 0; j < matches.size(); j++) {
                catalog_id = matches[j].first;
                tag_id     = matches[j].second;

                //if (tag_id == 10970) {
                //    cerr << "Hit the tag.\n";
                //}

                //// if (catalog_id != 3080) continue;

                if (catalog.count(catalog_id) == 0) continue;

                cloc = catalog[catalog_id];
                loc  = stacks[tag_id];

                if (filter_confounded &&
                    ((double) cloc->confounded_cnt / (double) cloc->cnt > confounded_limit)) {
                    // cerr << "Catalog locus " << cloc->id << " is confounded; confounded cnt: "
                    //   << cloc->confounded_cnt << "; total: " << cloc->cnt
                    //   << "; freq: " << (double) cloc->confounded_cnt / (double)cloc->cnt << "\n";
                    loc->blacklisted = true;
                    conf_loci_cnt++;
                    continue;
                }

                d = pmap->datum(catalog_id, sample.id);

                if (d == NULL) continue;

                if (filter_lnl && cloc->lnl < lnl_limit) {
                    loc->blacklisted = true;
                    lnl_cnt++;
                    continue;
                }

                //
                // If this locus was matched to the catalog using a gapped alignment, adjust the sequence
                // and SNP positions.
                //
                if (d->cigar != NULL) {
                    seq_len = parse_cigar(d->cigar, cigar);
                    seq     = apply_cigar_to_seq(loc->con, cigar);
                    model   = apply_cigar_to_model_seq(loc->model, cigar);
                    loc->add_consensus(seq.c_str());
                    loc->add_model(model.c_str());
                    adjust_and_add_snps_for_gaps(cigar, loc);

                    for (uint k = 0; k < loc->reads.size(); k++) {
                        reads.push_back(loc->reads[k]);

                        adj_seq = new char[seq_len + 1];
                        apply_cigar_to_seq(adj_seq, seq_len, loc->reads[k], cigar);
                        loc->reads[k] = adj_seq;
                    }
                }

                prune_nucleotides(cloc, loc, d, log_snp_fh,
                                  nuc_cnt,
                                  unk_hom_cnt, unk_het_cnt,
                                  hom_unk_cnt, het_unk_cnt,
                                  hom_het_cnt, het_hom_cnt);

                //
                // Prune haplotypes from this locus.
                //
                if (prune_haplotypes) {
                    prune_mst_haplotypes(cloc, d, loc, pruned_mst_hap_cnt, log_hap_fh);

                    prune_locus_haplotypes(cloc, d, loc, pruned_hap_cnt, log_hap_fh);

                    if (loc->blacklisted)
                        blacklist_cnt++;
                }

                //
                // If this locus was matched to the catalog using a gapped alignment, de-adjust the sequence
                // and SNP positions for writing back to disk.
                //
                if (d->cigar != NULL) {
                    seq   = remove_cigar_from_seq(loc->con, cigar);
                    model = remove_cigar_from_seq(loc->model, cigar);
                    loc->add_consensus(seq.c_str());
                    loc->add_model(model.c_str());
                    remove_snps_from_gaps(cigar, loc);

                    for (uint k = 0; k < loc->reads.size(); k++) {
                        adj_seq = loc->reads[k];
                        loc->reads[k] = reads[k];
                        delete [] adj_seq;
                    }
                    reads.clear();
                }
            }
        }

        cerr << "done.\n";

        unsigned long int total = unk_hom_cnt + unk_het_cnt + hom_unk_cnt + het_unk_cnt + hom_het_cnt + het_hom_cnt;

        cerr << "Total nucleotides processed: " << nuc_cnt << "\n"
             << "  Total nucleotides converted: " << total << "\n"
             << "    Converted from unknown to homozygous:      " << unk_hom_cnt << " nucleotides.\n"
             << "    Converted from unknown to heterozygous:    " << unk_het_cnt << " nucleotides.\n"
             << "    Converted from homozygous to unknown:      " << hom_unk_cnt << " nucleotides.\n"
             << "    Converted from heterozygous to unknown:    " << het_unk_cnt << " nucleotides.\n"
             << "    Converted from homozygous to heterozygous: " << hom_het_cnt << " nucleotides.\n"
             << "    Converted from heterozygous to homozygous: " << het_hom_cnt << " nucleotides.\n"
             << "Pruned: "      << pruned_mst_hap_cnt << " haplotypes using a tree method.\n"
             << "Pruned: "      << pruned_hap_cnt     << " haplotypes using a rare haplotype method.\n"
             << "Blacklisted: " << blacklist_cnt      << " loci due to inability to call haplotypes.\n"
             << "Blacklisted: " << lnl_cnt            << " loci due to log likelihoods below threshold.\n"
             << "Blacklisted: " << conf_loci_cnt      << " confounded loci.\n";

        log_fh << sample.name           << "\t"
               << nuc_cnt        << "\t"
               << total          << "\t"
               << unk_hom_cnt    << "\t"
               << unk_het_cnt    << "\t"
               << hom_unk_cnt    << "\t"
               << het_unk_cnt    << "\t"
               << hom_het_cnt    << "\t"
               << het_hom_cnt    << "\t"
               << blacklist_cnt  << "\t"
               << conf_loci_cnt  << "\t"
               << lnl_cnt        << "\t"
               << pruned_hap_cnt << "\t"
               << pruned_mst_hap_cnt << "\n";

        cerr << "Writing modified stacks, SNPs, alleles to '" << out_path << "'...";

        //
        // Rewrite stacks, model outputs, and haplotypes.
        //
        write_results(sample.name, stacks);

        //
        // Free up memory
        //
        cerr << "Freeing memory...";
        map<int, Locus *>::iterator stack_it;
        for (stack_it = stacks.begin(); stack_it != stacks.end(); stack_it++)
            delete stack_it->second;
        stacks.clear();
        cerr << "done.\n";
    }

    log_fh.close();

    if (verbose) {
        log_snp_fh.close();
        log_hap_fh.close();
    }

    //
    // Free memory associated with the catalog and matches to the catalog.
    //
    for (map<int, CSLocus *>::iterator cat_it = catalog.begin(); cat_it != catalog.end(); cat_it++)
        delete cat_it->second;

    for (uint i = 0; i < catalog_matches.size(); i++)
        for (uint j = 0; j < catalog_matches[i].size(); j++)
            delete catalog_matches[i][j];

    cerr << "rxstacks is done.\n";
    return 0;
    IF_NDEBUG_CATCH_ALL_EXCEPTIONS
}

int
dist(string hap_1, string hap_2) {
    int   dist  = 0;
    const char *p     = hap_1.c_str();
    const char *q     = hap_2.c_str();
    const char *p_end = p + hap_1.length();
    const char *q_end = q + hap_2.length();

    //
    // Count the number of characters that are different
    // between the two sequences.
    //
    while (p < p_end && q < q_end) {
        dist += (*p == *q) ? 0 : 1;
        p++;
        q++;
    }

    return dist;
}

int
calc_lnl_means(map<int, CSLocus *> &catalog, PopMap<CSLocus> *pmap)
{
    map<int, CSLocus *>::iterator it;
    CSLocus *cloc;
    Datum  **d;
    uint     cnt, mid, tot;
    double   median, mean;
    vector<double> lnls;
    ofstream log_fh;

    if (lnl_dist) {
        stringstream log;
        log << "batch_" << batch_id << ".rxstacks_lnls.tsv";
        string log_path = out_path + log.str();
        log_fh.open(log_path.c_str(), ofstream::out);

        if (log_fh.fail()) {
            cerr << "Error opening log file '" << log_path << "'\n";
            exit(1);
        }
        log_fh << "# Catalog Locus\tMean\tMedian\n";
    }

    tot = 0;

    for (it = catalog.begin(); it != catalog.end(); it++) {
        cloc = it->second;

        d    = pmap->locus(cloc->id);
        cnt  = pmap->sample_cnt();
        mean = 0.0;
        lnls.clear();

        for (uint i = 0; i < cnt; i++) {
            if (d[i] == NULL) continue;

            lnls.push_back(d[i]->lnl);
            mean += d[i]->lnl;
        }

        if (lnls.size() == 0)
            continue;

        sort(lnls.begin(), lnls.end());
        mid    = lnls.size() / 2;
        median = lnls.size() % 2 == 0 ? lnls[mid-1] + lnls[mid] / 2.0 : lnls[mid];
        mean   = mean / (double) lnls.size();

        cloc->lnl = mean;

        //
        // If the mean log likelihood for this catalog locus is below the threshold, count it as
        // its constituent components will be filtered as encountered later.
        //
        if (filter_lnl && cloc->lnl < lnl_limit)
            tot++;

        if (lnl_dist)
            log_fh << cloc->id << "\t"
                   << mean << "\t"
                   << median << "\n";
    }

    if (lnl_dist)
        log_fh.close();

    //
    // Print number of catalog loci that are confounded and will be removed.
    //
    cerr << tot << " catalog loci will be removed from the analysis due to log likelihoods below the threshold.\n";

    return 0;
}

int
sum_haplotype_counts(map<int, CSLocus *> &catalog, PopMap<CSLocus> *pmap)
{
    map<int, CSLocus *>::iterator it;
    CSLocus *cloc;
    Datum  **d;
    uint     cnt;

    for (it = catalog.begin(); it != catalog.end(); it++) {
        cloc = it->second;

        d   = pmap->locus(cloc->id);
        cnt = pmap->sample_cnt();

        for (uint i = 0; i < cnt; i++) {
            if (d[i] == NULL) continue;

            if (d[i]->obshap.size() == 1) {
                if (cloc->hap_cnts.count(d[i]->obshap[0]) == 0)
                    cloc->hap_cnts[d[i]->obshap[0]]  = 2;
                else
                    cloc->hap_cnts[d[i]->obshap[0]] += 2;
            } else {
                for (uint j = 0; j < d[i]->obshap.size(); j++)
                    if (cloc->hap_cnts.count(d[i]->obshap[j]) == 0)
                        cloc->hap_cnts[d[i]->obshap[j]]  = 1;
                    else
                        cloc->hap_cnts[d[i]->obshap[j]] += 1;
            }
        }
    }

    return 0;
}

int
prune_mst_haplotypes(CSLocus *cloc, Datum *d, Locus *loc, unsigned long &pruned_hap_cnt, ofstream &log_fh)
{
    //
    // Create a minimum spanning tree in order to determine the minimum distance
    // between each node in the list.
    //
    mst::MinSpanTree *mst = new mst::MinSpanTree();

    map<string, int>::iterator it;
    vector<string> haps;
    for (it = cloc->hap_cnts.begin(); it != cloc->hap_cnts.end(); it++) {
        mst->add_node(it->first);
        haps.push_back(it->first);
    }

    //
    // We are going to connect nodes in the graph when a SNP occurs in one
    // of the positions of the haplotype.
    //
    mst::Node *n_1, *n_2;

    uint snp_pos = 0;
    for (uint i = 0; i < cloc->snps.size(); i++) {
        if (cloc->snps[i]->type != snp_type_het)
            continue;

        for (uint j = 0; j < haps.size(); j++) {
            for (uint k = j + 1; k < haps.size(); k++) {
                //
                // If these two haplotypes differ by this SNP (and only this SNP), connect them in the graph.
                //
                if (haps[j].at(snp_pos) != haps[k].at(snp_pos) &&
                    dist(haps[j], haps[k]) == 1) {
                    n_1 = mst->node(haps[j]);
                    n_2 = mst->node(haps[k]);
                    n_1->add_edge(n_2, 1);
                    n_2->add_edge(n_1, 1);
                }
            }
        }
        snp_pos++;
    }

    //
    // Build the minimum spanning tree.
    //
    mst->build_tree();

    //
    // Sort the haplotypes by read depth in this sample
    //
    vector<pair<string, double> > haplotypes;

    for (uint i = 0; i < d->obshap.size(); i++)
        haplotypes.push_back(make_pair(string(d->obshap[i]), (double) d->depth[i]));

    uint size = haplotypes.size();

    //
    // Sort according to haplotype frequency.
    //
    sort(haplotypes.begin(), haplotypes.end(), compare_pair_haplotype_rev);

    if (size <= 2) {
        delete mst;
        return 0;
    }

    //
    // Pull out the two most frequently occuring haplotypes.
    //
    string hap_1, hap_2;
    double hap_1_depth, hap_2_depth;
    hap_1       = haplotypes[size - 1].first;
    hap_1_depth = haplotypes[size - 1].second;
    haplotypes.pop_back();

    if (haplotypes[size - 2].second > haplotypes[size - 3].second) {
        hap_2       = haplotypes[size - 2].first;
        hap_2_depth = haplotypes[size - 2].second;
        haplotypes.pop_back();
    } else {
        hap_2       = "";
        hap_2_depth = 0.0;
    }

    //
    // For each remaining haplotpye, check if it can be merged into a node (haplotype) no
    // more than one nucleotide apart. If there is more than one, merge it into the more
    // frequently occuring haplotype.
    //
    string hap, src_hap, dest_hap, label;
    double max, weighted;

    for (uint i = 0; i < haplotypes.size(); i++) {

        //
        // Find the current haplotype in the MST.
        //
        n_1 = mst->node(haplotypes[i].first);

        max      = 0.0;
        hap      = "";
        weighted = 0.0;
        //
        // Check any potential edges in the graph for merging.
        //
        for (uint j = 0; j < n_1->edges.size(); j++) {
            label = n_1->edges[j]->child->label;

            if (label == hap_1) {
                weighted = (double) cloc->hap_cnts[label] * log(hap_1_depth);
                // cerr << "Cloc hap: " << label << "; popcnt: " << cloc->hap_cnts[label] << "; hap depth: " << hap_1_depth << "; weighted: " << weighted << "\n";
            } else if (label == hap_2) {
                weighted = (double) cloc->hap_cnts[label] * log(hap_2_depth);
                // cerr << "Cloc hap: " << label << "; popcnt: " << cloc->hap_cnts[label] << "; hap depth: " << hap_2_depth << "; weighted: " << weighted << "\n";
            } else
                continue;

            if (weighted == max) {
                //
                // There is more than one identical possibility, we can do no more.
                //
                hap = "";
                break;

            } else if (weighted > max) {
                max = weighted;
                hap = label;
            }
        }

        if (hap.length() == 0)
            continue;

        src_hap  = convert_catalog_haplotype_to_sample(haplotypes[i].first, cloc, loc);
        dest_hap = convert_catalog_haplotype_to_sample(hap, cloc, loc);

        if (verbose) {
            #pragma omp critical
            log_fh << cloc->id            << "\t"
                   << loc->sample_id      << "\t"
                   << loc->id             << "\t"
                   << src_hap             << "\t"
                   << haplotypes[i].first << "\t"
                   << dest_hap            << "\t"
                   << hap                 << "\t"
                   << "mst"               << "\n";
        }
        pruned_hap_cnt++;

        //
        // Remove the haplotype.
        //
        it = loc->alleles.find(src_hap);

        if (it != loc->alleles.end()) {
            loc->alleles.erase(it);
        }

        //
        // Add to the count of the merged-to haplotype.
        //
        if (loc->alleles.count(dest_hap) > 0) {
            loc->alleles[dest_hap]++;
        } else {
            cerr << "Error finding allele\n";
        }
    }

    //
    // Update the matched haplotypes in the Datum object, so the haplotype pruner can
    // operate on newly generated, spurious haplotypes.
    //
    generate_matched_haplotypes(cloc, loc, d);
    delete mst;

    return 0;
}

int
prune_locus_haplotypes(CSLocus *cloc, Datum *d, Locus *loc, unsigned long &pruned_hap_cnt, ofstream &log_fh)
{
    if (d->obshap.size() < 2) return 0;

    //
    // Identify the two most frequent haplotypes in this sample.
    //
    vector<pair<string, double> > haplotypes;
    double weighted_hap;

    for (uint i = 0; i < d->obshap.size(); i++) {
        //
        // Lookup the number of occurrences of this haplotype in the
        // population as well as the depth of the haplotype in this indiviudal.
        // We will weight the occurrences of the haplotype in the population by the natural log
        // of the read depth of the haplotype in this individual, storing the result.
        //
        weighted_hap = (double) cloc->hap_cnts[d->obshap[i]] * log((double) d->depth[i]);
        haplotypes.push_back(make_pair(string(d->obshap[i]), weighted_hap));
    }

    //
    // Sort according to haplotype frequency.
    //
    sort(haplotypes.begin(), haplotypes.end(), compare_pair_haplotype);

    if (haplotypes.size() == 0) {
        cerr << "Error processing catalog locus " << cloc->id << "\n";
        return -1;
    }

    //
    // Prune out excess haplotypes.
    //
    for (uint i = 2; i < haplotypes.size(); i++) {
        //
        // Make sure that those haplotypes we want to discard occur at a frequency lower
        // than the second most frequent haplotype, instead of being tied for second.
        //
        if (haplotypes[i].second >= haplotypes[1].second ||
            (max_haplotype_cnt > 0 && haplotypes[i].second > max_haplotype_cnt))
            continue;

        remove_haplotype(cloc, loc, haplotypes[i].first, pruned_hap_cnt, log_fh, "rare_step_1");
        haplotypes.erase(haplotypes.begin() + i);
    }

    //
    // If there are more than two haplotypes remaining and the second, third, etc
    // haplotype exist only in this individual, prune them out.
    //

    if (haplotypes.size() > 2) {
        int    stop_pos  = haplotypes.size() - 1;
        int    start_pos = stop_pos;
        double score     = haplotypes[stop_pos].second;
        while (start_pos > 1) {
            if (cloc->hap_cnts[haplotypes[start_pos].first] == 1 &&
                haplotypes[start_pos].second == score)
                start_pos--;
            else
                break;
        }

        if (start_pos < stop_pos) {
            for (int i = start_pos; i <= stop_pos; i++)
                remove_haplotype(cloc, loc, haplotypes[i].first, pruned_hap_cnt, log_fh, "rare_step_1");
        }
    }

    //
    // Update the matched haplotypes in the Datum object, so the haplotype pruner can
    // operate on newly generated, spurious haplotypes.
    //
    generate_matched_haplotypes(cloc, loc, d);

    return 0;
}


string
convert_catalog_haplotype_to_sample(string cat_haplotype, CSLocus *cloc, Locus *loc)
{
    int cat_snp =  0;
    int cat_idx = -1;
    int loc_snp =  0;
    int loc_idx = -1;
    int k       = -1;
    int j       = -1;
    string hap;

    do {
        j++;
        loc_idx++;
        //
        // Advance to a het in the sample locus.
        //
        while (j < (int) loc->snps.size() && loc->snps[j]->type != snp_type_het) j++;
        if (j >= (int) loc->snps.size()) break;
        loc_snp = loc->snps[j]->col;

        do {
            k++;
            cat_idx++;
            //
            // Advance to the het in the catalog locus that corresponds to the sample locus.
            //
            while (k < (int) cloc->snps.size() && cloc->snps[k]->type != snp_type_het) k++;
            if (k >= (int) cloc->snps.size()) break;
            cat_snp = cloc->snps[k]->col;

        } while (cat_snp < loc_snp);

        //
        // Extract out the nucleotide from the catalog haplotype that matches the sample
        // haplotype. For example, catalog haplotype may be 'ACGTG' while sample haplotype
        // is 'CT'.
        //
        if (j < (int) loc->snps.size() && k < (int) cloc->snps.size() && cat_snp == loc_snp) {
            hap += cat_haplotype.at(cat_idx);
        } else {
            cerr << "Error processing catalog locus " << cloc->id << "\n";
            return "";
        }

    } while (j < (int) loc->snps.size());

    return hap;
}

int
remove_haplotype(CSLocus *cloc, Locus *loc, string haplotype,
                 unsigned long &pruned_hap_cnt, ofstream &log_fh, string alg_type)
{
    map<string, int>::iterator it;
    string hap = "";

    hap = convert_catalog_haplotype_to_sample(haplotype, cloc, loc);

    if (verbose) {
        #pragma omp critical
        log_fh << cloc->id       << "\t"
               << loc->sample_id << "\t"
               << loc->id        << "\t"
               << hap            << "\t"
               << haplotype      << "\t"
               << "\t"
               << "\t"
               << alg_type       << "\n";
    }

    //
    // Remove the haplotype.
    //
    it = loc->alleles.find(hap);

    if (it != loc->alleles.end()) {
        loc->alleles.erase(it);
        pruned_hap_cnt++;
    }

    //
    // Decrement the count for this haplotype in the catalog locus.
    //
    if (cloc->hap_cnts.count(haplotype) > 0)
        cloc->hap_cnts[haplotype]--;

    return 0;
}

int
prune_nucleotides(CSLocus *cloc, Locus *loc, Datum *d, ofstream &log_fh, unsigned long int &nuc_cnt,
                  unsigned long int &unk_hom_cnt, unsigned long int &unk_het_cnt,
                  unsigned long int &hom_unk_cnt, unsigned long int &het_unk_cnt,
                  unsigned long int &hom_het_cnt, unsigned long int &het_hom_cnt)
{
    map<char, int>      nucs;
    set<char>           cnucs;
    set<char>::iterator it;
    set<int>            rows;

    nuc_cnt += loc->len;

    for (uint i = 0; i < loc->snps.size() && i < cloc->snps.size(); i++) {

        //
        // Safety checks.
        //
        if (loc->snps[i]->col != cloc->snps[i]->col)
            cerr << "Warning: comparing mismatched SNPs in catalog locus " << cloc->id
                 << " and sample " << loc->sample_id << ", locus " << loc->id << "; col: " << loc->snps[i]->col << "\n";
        if (loc->snps[i]->type == snp_type_het && cloc->snps[i]->type == snp_type_hom)
            cerr << "Warning: sample locus is variable while catalog locus is fixed; catalog locus " << cloc->id
                 << " and sample " << loc->sample_id << ", locus " << loc->id << "; col: " << loc->snps[i]->col << "\n";

        //
        // Either there is an unknown call in locus, or, there is a snp in the catalog and any state in the locus.
        //
        if ((loc->snps[i]->type == snp_type_unk) ||
            (cloc->snps[i]->type == snp_type_het && loc->snps[i]->type == snp_type_hom)) {

            // cerr << "  Looking at SNP call in tag " << loc->id << " at position " << i << "; col: " << loc->snps[i]->col << "\n"
            //   << "    Catalog column: " << cloc->snps[i]->col << " (" << i << "); Sample column: " << loc->snps[i]->col << " (" << i << ")\n"
            //   << "    Sample has model call type: " << (loc->snps[i]->type == snp_type_unk ? "Unknown" : "Homozygous") << "; nucleotides: '"
            //   << loc->snps[i]->rank_1 << "' and '" << loc->snps[i]->rank_2 << "'; model call: " << loc->model[loc->snps[i]->col] << "\n"
            //   << "    Catalog has model call type: " << (cloc->snps[i]->type == snp_type_het ? "Heterozygous" : "Homozygous") << "; nucleotides: '"
            //   << cloc->snps[i]->rank_1 << "' and '" << cloc->snps[i]->rank_2 << "'\n";

            if (loc->snps[i]->rank_1 == 'N' || cloc->snps[i]->rank_1 == 'N' ||
                loc->snps[i]->rank_1 == '-' || cloc->snps[i]->rank_1 == '-')
                continue;

            cnucs.insert(cloc->snps[i]->rank_1);
            if (cloc->snps[i]->rank_2 != 0) cnucs.insert(cloc->snps[i]->rank_2);
            if (cloc->snps[i]->rank_3 != 0) cnucs.insert(cloc->snps[i]->rank_3);
            if (cloc->snps[i]->rank_4 != 0) cnucs.insert(cloc->snps[i]->rank_4);

            // cerr << "    Catalog has nucleotides: ";
            // for (it = cnucs.begin(); it != cnucs.end(); it++)
            //  cerr << *it << ", ";
            // cerr << "\n";

            //
            // Tally the number of occurances of each nucleotide also present in the
            // catalog in order to fuel the snp calling model.
            //
            // Note reads that contain nucleotides not present in the catalog so they
            // can be excluded when calling haplotypes from the read.
            //
            nucs['A'] = 0;
            nucs['C'] = 0;
            nucs['G'] = 0;
            nucs['T'] = 0;
            nucs['N'] = 0;

            for (uint k = 0; k < loc->reads.size(); k++) {
                if (cnucs.count(loc->reads[k][i]) > 0)
                    nucs[loc->reads[k][i]]++;
                else if (loc->reads[k][i] != 'N')
                    rows.insert(k);
            }

            //
            // Test pruned data for homozygosity or heterozygosity.
            //
            invoke_model(loc, i, nucs);

            cnucs.clear();
        }
    }

    if (verbose) {
        #pragma omp critical
        log_model_calls(loc, log_fh,
                        unk_hom_cnt, unk_het_cnt,
                        hom_unk_cnt, het_unk_cnt,
                        hom_het_cnt, het_hom_cnt);
    } else {
        log_model_calls(loc, log_fh,
                        unk_hom_cnt, unk_het_cnt,
                        hom_unk_cnt, het_unk_cnt,
                        hom_het_cnt, het_hom_cnt);
    }

    //
    // Re-call alleles.
    //
    loc->alleles.clear();

    call_alleles(loc, rows);

    //
    // If SNPs were called at this locus but no alleles could be determined,
    // blacklist this tag. This can occur if a fixed nucleotide isn't captured in
    // the catalog and all the reads are removed for the purpose of reading haplotypes.
    //
    if (loc->alleles.size() <= 1)
        for (uint j = 0; j < loc->snps.size(); j++)
            if (loc->snps[j]->type == snp_type_het) {
                loc->blacklisted = 1;
                break;
            }

    //
    // Update the matched haplotypes in the Datum object, so the haplotype pruner can
    // operate on newly generated, spurious haplotypes.
    //
    generate_matched_haplotypes(cloc, loc, d);

    return 0;
}

int
invoke_model(Locus *loc, int col, map<char, int> &nucs)
{
    //
    // Search this column for the presence of a SNP
    //
    switch(model_type) {
    case snp:
        call_multinomial_snp(loc, col, nucs);
        break;
    case bounded:
        call_bounded_multinomial_snp(loc, col, nucs);
        break;
    default:
        break;
    }

    return 0;
}

int
call_alleles(Locus *loc, set<int> &rows)
{
    int      row;
    int      height = loc->reads.size();
    string   allele;
    char     base;
    vector<SNP *>::iterator snp;

    for (row = 0; row < height; row++) {
        //
        // If a read had a nucleotide not present in the catalog, do not call
        // a haplotype from it.
        //
        if (rows.count(row) > 0)
            continue;

        allele.clear();

        uint snp_cnt = 0;

        for (snp = loc->snps.begin(); snp != loc->snps.end(); snp++) {
            if ((*snp)->type != snp_type_het) continue;

            snp_cnt++;

            base = loc->reads[row][(*snp)->col];

            //
            // Check to make sure the nucleotide at the location of this SNP is
            // of one of the two possible states the multinomial model called.
            //
            if (base == (*snp)->rank_1 || base == (*snp)->rank_2)
                allele += base;
            else
                break;
        }

        if (snp_cnt > 0 && allele.length() == snp_cnt)
            loc->alleles[allele]++;
    }

    return 0;
}

int
generate_matched_haplotypes(CSLocus *cloc, Locus *loc, Datum *d)
{
    //
    // Free the existing matched haplotypes.
    //
    for (uint i = 0; i < d->obshap.size(); i++)
        delete [] d->obshap[i];
    d->obshap.clear();
    d->depth.clear();

    //
    // Construct a set of haplotypes from the locus relative to the catalog locus.
    // (The locus already has a set of haplotypes, however, they don't necessarily
    //  account for all the SNPs in the catalog, so we will augment them with sequence
    //  from the consensus.)
    //
    vector<pair<string, SNP *> >   merged_snps;
    map<int, pair<string, SNP *> > columns;
    map<int, pair<string, SNP *> >::iterator c;
    vector<pair<string, SNP *> >::iterator   k;

    for (uint i = 0; i < cloc->snps.size(); i++) {
        if (cloc->snps[i]->type != snp_type_het)
            continue;
        columns[cloc->snps[i]->col] = make_pair("catalog", cloc->snps[i]);
    }

    for (uint i = 0; i < loc->snps.size();  i++) {
        if (loc->snps[i]->type != snp_type_het)
            continue;

        //
        // Is this column already represented in the catalog?
        //
        if (columns.count(loc->snps[i]->col))
            columns[loc->snps[i]->col] = make_pair("both", loc->snps[i]);
        else
            columns[loc->snps[i]->col] = make_pair("query", loc->snps[i]);
    }

    for (c = columns.begin(); c != columns.end(); c++)
        merged_snps.push_back((*c).second);

    //
    // Sort the SNPs by column
    //
    sort(merged_snps.begin(), merged_snps.end(), compare_pair_snp);

    map<string, int>::iterator b;
    string old_allele, new_allele;
    int    pos;

    for (b = loc->alleles.begin(); b != loc->alleles.end(); b++) {
        old_allele = b->first;
        new_allele = "";
        pos        = 0;

        for (k = merged_snps.begin(); k != merged_snps.end(); k++) {
            //
            // If the SNPs from the catalog haplotype beyond the length of the query, add Ns
            //
            if (k->first == "catalog") {
                new_allele += (k->second->col > loc->len - 1) ? 'N' : loc->con[k->second->col];
            } else {
                new_allele += old_allele[pos];
                pos++;
            }
        }

        char *h = new char[new_allele.length() + 1];
        strcpy(h, new_allele.c_str());
        d->obshap.push_back(h);
        d->depth.push_back(b->second);

        // loc->alleles[new_allele] = b->second;
        // cerr << "Adding haplotype: " << new_allele << " [" << b->first << "]\n";
    }

    return 0;
}

int
log_model_calls(Locus *loc, ofstream &log_fh,
                unsigned long int &unk_hom_cnt, unsigned long int &unk_het_cnt,
                unsigned long int &hom_unk_cnt, unsigned long int &het_unk_cnt,
                unsigned long int &hom_het_cnt, unsigned long int &het_hom_cnt)
{
    //
    // Log model call changes
    //
    for (uint j = 0; j < loc->snps.size(); j++) {
        switch(loc->model[j]) {
        case 'U':
            switch(loc->snps[j]->type) {
            case snp_type_het:
                if (verbose)
                    log_fh << loc->sample_id << "\t" << loc->id << "\t" << loc->snps[j]->col << "\t" << 'U' << "\t" << 'E' << "\n";
                unk_het_cnt++;
                break;
            case snp_type_hom:
                if (verbose)
                    log_fh << loc->sample_id << "\t" << loc->id << "\t" << loc->snps[j]->col << "\t" << 'U' << "\t" << 'O' << "\n";
                unk_hom_cnt++;
                break;
            case snp_type_unk:
            default:
                break;
            }
            break;
        case 'E':
            switch(loc->snps[j]->type) {
            case snp_type_het:
                break;
            case snp_type_hom:
                if (verbose)
                    log_fh << loc->sample_id << "\t" << loc->id << "\t" << loc->snps[j]->col << "\t" << 'E' << "\t" << 'O' << "\n";
                het_hom_cnt++;
                break;
            case snp_type_unk:
            default:
                if (verbose)
                    log_fh << loc->sample_id << "\t" << loc->id << "\t" << loc->snps[j]->col << "\t" << 'E' << "\t" << 'U' << "\n";
                het_unk_cnt++;
                break;
            }
            break;
        case 'O':
        default:
            switch(loc->snps[j]->type) {
            case snp_type_het:
                if (verbose)
                    log_fh << loc->sample_id << "\t" << loc->id << "\t" << loc->snps[j]->col << "\t" << 'O' << "\t" << 'E' << "\n";
                hom_het_cnt++;
                break;
            case snp_type_hom:
                break;
            case snp_type_unk:
            default:
                if (verbose)
                    log_fh << loc->sample_id << "\t" << loc->id << "\t" << loc->snps[j]->col << "\t" << 'O' << "\t" << 'U' << "\n";
                hom_unk_cnt++;
                break;
            }
            break;
        }
    }

    return 0;
}

int
write_results(string file, map<int, Locus *> &m)
{
    map<int, Locus *>::iterator i;
    vector<char *>::iterator    j;
    vector<int>::iterator       k;
    map<string, int>::iterator  t;
    Locus                      *tag_1;
    stringstream                sstr;

    bool gzip = (in_file_type == FileT::gzsql) ? true : false;

    //
    // Parse the input file name to create the output files
    //
    string tag_file = out_path + file + ".tags.tsv";
    string snp_file = out_path + file + ".snps.tsv";
    string all_file = out_path + file + ".alleles.tsv";
    string mod_file = out_path + file + ".models.tsv";

    if (gzip) {
        tag_file += ".gz";
        snp_file += ".gz";
        all_file += ".gz";
        mod_file += ".gz";
    }

    //
    // Open the output files for writing.
    //
    gzFile   gz_tags=NULL, gz_snps=NULL, gz_alle=NULL, gz_mods=NULL;
    ofstream tags, snps, alle, mods;
    if (gzip) {
        gz_tags = gzopen(tag_file.c_str(), "wb");
        if (!gz_tags) {
            cerr << "Error: Unable to open gzipped tag tag file '" << tag_file << "': " << strerror(errno) << ".\n";
            exit(1);
        }
        #if ZLIB_VERNUM >= 0x1240
        gzbuffer(gz_tags, libz_buffer_size);
        #endif
        gz_mods = gzopen(mod_file.c_str(), "wb");
        if (!gz_mods) {
            cerr << "Error: Unable to open gzipped model file '" << mod_file << "': " << strerror(errno) << ".\n";
            exit(1);
        }
        #if ZLIB_VERNUM >= 0x1240
        gzbuffer(gz_mods, libz_buffer_size);
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
            cerr << "Error: Unable to open tag file for writing.\n";
            exit(1);
        }
        mods.open(mod_file.c_str());
        if (mods.fail()) {
            cerr << "Error: Unable to open model file for writing.\n";
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
    log << "# rxstacks version " << VERSION << "; generated on " << date << "\n";
    if (gzip) {
        gzputs(gz_tags, log.str().c_str());
        gzputs(gz_mods, log.str().c_str());
        gzputs(gz_snps, log.str().c_str());
        gzputs(gz_alle, log.str().c_str());
    } else {
        tags << log.str();
        mods << log.str();
        snps << log.str();
        alle << log.str();
    }

    int wrote = 0;

    for (i = m.begin(); i != m.end(); i++) {
        tag_1 = i->second;

        wrote++;

        // First write the consensus sequence
        sstr << "0" << "\t"
             << tag_1->sample_id << "\t"
             << tag_1->id << "\t"
             << tag_1->loc.chr() << "\t"
             << tag_1->loc.bp << "\t"
             << (tag_1->loc.strand == strand_plus ? "+" : "-") << "\t"
             << "consensus" << "\t"
             << "\t"
             << "\t"
             << tag_1->con << "\t"
             << (tag_1->deleveraged == true ? "1" : "0")     << "\t"
             << (tag_1->blacklisted == true ? "1" : "0")     << "\t"
             << (tag_1->lumberjackstack == true ? "1" : "0") << "\t"
             << tag_1->lnl << "\n";

        //
        // Write a sequence recording the output of the SNP model for each nucleotide.
        //
        sstr << "0"              << "\t"
             << tag_1->sample_id << "\t"
             << tag_1->id        << "\t"
             << "\t"
             << "\t"
             << "\t"
             << "model"          << "\t"
             << "\t"
             << "\t";
        for (uint j = 0; j < tag_1->snps.size(); j++) {
            switch(tag_1->snps[j]->type) {
            case snp_type_het:
                sstr << "E";
                break;
            case snp_type_hom:
                sstr << "O";
                break;
            default:
                sstr << "U";
                break;
            }
        }
        sstr << "\t"
             << "\t"
             << "\t"
             << "\t"
             << "\n";

        if (gzip) gzputs(gz_tags, sstr.str().c_str()); else tags << sstr.str();
        if (gzip) gzputs(gz_mods, sstr.str().c_str()); else mods << sstr.str();
        sstr.str("");

        //
        // Now write out each read from this locus.
        //
        for (uint j = 0; j < tag_1->reads.size(); j++) {
            sstr << "0"                 << "\t"
                 << tag_1->sample_id    << "\t"
                 << tag_1->id           << "\t\t\t\t";

            if (tag_1->comp_type[j] == primary) {
                sstr << "primary" << "\t"
                     << tag_1->comp_cnt[j]  << "\t";
            } else {
                sstr << "secondary" << "\t\t";
            }

            sstr << tag_1->comp[j]      << "\t"
                 << tag_1->reads[j]     << "\t\t\t\t\n";
        }

        if (gzip) gzputs(gz_tags, sstr.str().c_str()); else tags << sstr.str();
        sstr.str("");

        //
        // Write out the model calls for each nucleotide in this locus.
        //
        for (uint j = 0; j < tag_1->snps.size(); j++) {
            sstr << "0"                 << "\t"
                 << tag_1->sample_id    << "\t"
                 << tag_1->id           << "\t"
                 << tag_1->snps[j]->col << "\t";

            switch(tag_1->snps[j]->type) {
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

            sstr << std::fixed   << std::setprecision(3)
                 << tag_1->snps[j]->lratio << "\t"
                 << tag_1->snps[j]->rank_1 << "\t"
                 << (tag_1->snps[j]->rank_2 == 0 ? '-' : tag_1->snps[j]->rank_2) << "\t\t\n";
        }

        if (gzip) gzputs(gz_snps, sstr.str().c_str()); else snps << sstr.str();
        sstr.str("");

        //
        // Write the expressed alleles seen for the recorded SNPs and
        // the percentage of tags a particular allele occupies.
        //
        char pct[id_len];
        for (t = tag_1->alleles.begin(); t != tag_1->alleles.end(); t++) {
            sprintf(pct, "%.2f", ((t->second/double(tag_1->reads.size())) * 100));
            sstr << "0"              << "\t"
                 << tag_1->sample_id << "\t"
                 << tag_1->id        << "\t"
                 << t->first         << "\t"
                 << pct              << "\t"
                 << t->second        << "\n";
        }

        if (gzip) gzputs(gz_alle, sstr.str().c_str()); else alle << sstr.str();
        sstr.str("");
    }

    if (gzip) {
        gzclose(gz_tags);
        gzclose(gz_mods);
        gzclose(gz_snps);
        gzclose(gz_alle);
    } else {
        tags.close();
        mods.close();
        snps.close();
        alle.close();
    }

    cerr << "wrote " << wrote << " loci.\n";

    return 0;
}

int
fill_catalog_snps(map<int, CSLocus *> &catalog)
{
    map<int, CSLocus *>::iterator it;
    CSLocus *cloc;

    for (it = catalog.begin(); it != catalog.end(); it++) {
        cloc = it->second;

        queue<SNP *> snps;
        for (uint j = 0; j < cloc->snps.size(); j++)
            snps.push(cloc->snps[j]);

        cloc->snps.clear();

        for (uint j = 0; j < cloc->len; j++) {
            if (snps.size() > 0 && snps.front()->col == j) {
                cloc->snps.push_back(snps.front());
                snps.pop();
            } else {
                SNP *snp = new SNP;
                snp->type   = snp_type_hom;
                snp->col    = j;
                snp->lratio = 0;
                snp->rank_1 = cloc->con[j];
                snp->rank_2 = 0;
                cloc->snps.push_back(snp);
            }
        }
    }

    return 0;
}

int
init_log(int argc, char **argv, ofstream &log_fh, ofstream &log_snp_fh, ofstream &log_hap_fh)
{
    stringstream log;
    stringstream sstr;

    //
    // Open the log files.
    //
    log << out_path << "batch_" << batch_id << ".rxstacks.log";
    log_fh.open(log.str().c_str(), ofstream::out);

    if (log_fh.fail()) {
        cerr << "Error opening log file '" << log.str() << "'\n";
        exit(1);
    }

    if (verbose) {
        log.str("");
        log << out_path << "batch_" << batch_id << ".rxstacks.snps.log";
        log_snp_fh.open(log.str().c_str(), ofstream::out);

        if (log_snp_fh.fail()) {
            cerr << "Error opening log file '" << log.str() << "'\n";
            exit(1);
        }

        log.str("");
        log << out_path << "batch_" << batch_id << ".rxstacks.haplotypes.log";
        log_hap_fh.open(log.str().c_str(), ofstream::out);

        if (log_hap_fh.fail()) {
            cerr << "Error opening log file '" << log.str() << "'\n";
            exit(1);
        }
    }

    //
    // Obtain the current date.
    //
    time_t     rawtime;
    struct tm *timeinfo;
    char       date[32];
    time(&rawtime);
    timeinfo = localtime(&rawtime);
    strftime(date, 32, "%F %T", timeinfo);

    sstr << "#";
    for (int i = 0; i < argc; i++)
        sstr << " " << argv[i];
    sstr << "\n" << "# rxstacks executed " << date;

    log_fh << sstr.str() << "\n"
           << "# Sample\t"
           << "Total nucs\t"
           << "Total nucs converted\t"
           << "Unk to Hom\t"
           << "Unk to Het\t"
           << "Hom to Unk\t"
           << "Het to Unk\t"
           << "Hom to Het\t"
           << "Het to Hom\t"
           << "Confounded loci\t"
           << "Lnl Filtered loci\t"
           << "Pruned Haplotypes\t"
           << "MST-Pruned Haplotypes\n";

    if (verbose) {
        log_snp_fh << sstr.str() << "\n"
                   << "# Sample Id\t"
                   << "Locus ID\t"
                   << "SNP Col\t"
                   << "Orig Value\t"
                   << "Corr Value\n";
        log_hap_fh << sstr.str() << "\n"
                   << "# Catalog Locus\t"
                   << "Sample\t"
                   << "Sample Locus\t"
                   << "Sample Haplotype\t"
                   << "Catalog Haplotype\t"
                   << "Corrected Sample Haplotype\t"
                   << "Corrected Catalog Haplotype\t"
                   << "Algorithm\n";
    }

    return 0;
}

int
parse_command_line(int argc, char* argv[])
{
    int c;

    while (1) {
        static struct option long_options[] = {
            {"help",         no_argument,       NULL, 'h'},
            {"version",      no_argument,       NULL, 'v'},
            {"conf_filter",  no_argument,       NULL, 'F'},
            {"prune_haplo",  no_argument,       NULL, 'H'},
            {"lnl_filter",   no_argument,       NULL, 'G'},
            {"lnl_dist",     no_argument,       NULL, 'D'},
            {"verbose",      no_argument,       NULL, 'V'},
            {"num_threads",  required_argument, NULL, 't'},
            {"batch_id",     required_argument, NULL, 'b'},
            {"in_path",      required_argument, NULL, 'P'},
            {"outpath",      required_argument, NULL, 'o'},
            {"model_type",   required_argument, NULL, 'T'},
            {"bound_low",    required_argument, NULL, 'L'},
            {"bound_high",   required_argument, NULL, 'U'},
            {"alpha",        required_argument, NULL, 'A'},
            {"conf_lim",     required_argument, NULL, 'C'},
            {"max_haplo",    required_argument, NULL, 'M'},
            {"lnl_lim",      required_argument, NULL, 'I'},
            {0, 0, 0, 0}
        };

        // getopt_long stores the option index here.
        int option_index = 0;

        c = getopt_long(argc, argv, "hvVFGDHo:t:b:P:T:L:U:A:C:I:", long_options, &option_index);

        // Detect the end of the options.
        if (c == -1)
            break;

        switch (c) {
        case 'h':
            help();
            break;
        case 't':
            num_threads = atoi(optarg);
            break;
        case 'P':
            in_path = optarg;
            break;
        case 'b':
            batch_id = is_integer(optarg);
            break;
        case 'o':
            out_path = optarg;
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
        case 'L':
            bound_low  = atof(optarg);
            break;
        case 'U':
            bound_high = atof(optarg);
            break;
        case 'A':
            alpha = atof(optarg);
            break;
        case 'F':
            filter_confounded = true;
            break;
        case 'C':
            confounded_limit  = is_double(optarg);
            filter_confounded = true;
            break;
        case 'H':
            prune_haplotypes = true;
            break;
        case 'M':
            max_haplotype_cnt = is_integer(optarg);
            break;
        case 'G':
            filter_lnl = true;
            break;
        case 'I':
            lnl_limit  = is_double(optarg);
            filter_lnl = true;
            break;
        case 'D':
            lnl_dist = true;
            break;
        case 'V':
            verbose = true;
            break;
        case 'v':
            version();
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

    if (in_path.length() == 0) {
        cerr << "You must specify a path to the directory containing Stacks output files.\n";
        help();
    }

    if (out_path.length() == 0) {
        cerr << "You must specify a path to a directory where to write output files.\n";
        help();
    }

    if (in_path.at(in_path.length() - 1) != '/')
        in_path += "/";

    if (out_path.at(out_path.length() - 1) != '/')
        out_path += "/";

    if (batch_id < 0) {
        vector<int> cat_ids = find_catalogs(in_path);
        if (cat_ids.size() == 1) {
            batch_id = cat_ids[0];
        } else if (cat_ids.empty()) {
            cerr << "Error: Unable to find a catalog in '" << in_path << "'.\n";
            help();
        } else {
            cerr << "Error: Input directory contains several catalogs, please specify -b.\n";
            help();
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

    if (filter_confounded == true && (confounded_limit < 0 || confounded_limit > 1.0)) {
        cerr << "Confounded locus limit is a percentage and must be between 0.0 and 1.0.\n";
        help();
    }

    return 0;
}

void version() {
    cerr << "rxstacks " << VERSION << "\n\n";

    exit(1);
}

void help() {
    cerr << "rxstacks " << VERSION << "\n"
              << "rxstacks -P path -o path [-t threads] [-b batch_id]" << "\n"
              << "  P: path to the Stacks output files.\n"
              << "  o: output path to write results ('.' to override the current files).\n"
              << "  t: number of threads to run in parallel sections of code.\n"
              << "  b: database/batch ID of the input catalog to consider (default: guess).\n"
              << "  Filtering options:\n"
              << "    --lnl_filter: filter catalog loci based on the mean log likelihood of the catalog locus in the population.\n"
              << "      --lnl_lim <limit>: minimum log likelihood required to keep a catalog locus.\n"
              << "      --lnl_dist: print distribution of mean log likelihoods for catalog loci.\n"
              << "    --conf_filter: filter confounded loci.\n"
              << "    --conf_lim <limit>: between 0.0 and 1.0 (default 0.75), proportion of loci in population that must be confounded relative to the catalog locus.\n"
              << "    --prune_haplo: prune out non-biological haplotypes unlikely to occur in the population.\n"
              << "    --max_haplo <limit>: only consider haplotypes for pruning if they occur in fewer than max_haplo_cnt samples.\n"
              << "  Model options:\n"
              << "    --model_type <type>: either 'snp' (default), 'bounded', or 'fixed'\n"
              << "    For the SNP or Bounded SNP model:\n"
              << "      --alpha <num>: chi square significance level required to call a heterozygote or homozygote, either 0.1 (default), 0.05, 0.01, or 0.001.\n"
              << "    For the Bounded SNP model:\n"
              << "      --bound_low <num>: lower bound for epsilon, the error rate, between 0 and 1.0 (default 0).\n"
              << "      --bound_high <num>: upper bound for epsilon, the error rate, between 0 and 1.0 (default 1).\n"
              << "  Logging Options:\n"
              << "      --verbose: extended logging, including coordinates of all changed nucleotides (forces single-threaded execution).\n";
    exit(1);
}
