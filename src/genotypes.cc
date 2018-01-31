// -*-mode:c++; c-style:k&r; c-basic-offset:4;-*-
//
// Copyright 2011-2014, Julian Catchen <jcatchen@uoregon.edu>
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
// genotypes -- genotype a mapping cross.
//
// Julian Catchen
// jcatchen@uoregon.edu
// University of Oregon
//

#include "MetaPopInfo.h"

#include "genotypes.h"

// Global variables to hold command-line options.
int       num_threads =  1;
int       batch_id    = -1;
map_types map_type    = gen;
out_types out_type    = joinmap;
string    in_path;
string    out_path;
string    cor_path;
string    out_file;
string    bl_file;
string    wl_file;
string    enz;
int       progeny_limit    = 1;
bool      man_corrections  = false;
bool      corrections      = false;
bool      expand_id        = false;
bool      sql_out          = false;
bool      filter_lnl       = false;
double    lnl_limit        = 0.0;
double    chisq_pval_limit = 0.05;
int       min_stack_depth  = 0;
int       min_hom_seqs     = 5;
double    min_het_seqs     = 0.05;
double    max_het_seqs     = 0.1;

set<int> whitelist, blacklist;

//
// Hold information about restriction enzymes
//
map<string, const char **> renz;
map<string, int>           renz_cnt;
map<string, int>           renz_len;

//
// Dictionaries to hold legal genotypes for different map types.
//
map<string, map<string, string> > global_dictionary;

int main (int argc, char* argv[]) {
    IF_NDEBUG_TRY

    initialize_renz(renz, renz_cnt, renz_len);

    parse_command_line(argc, argv);

    //
    // Set the number of OpenMP parallel threads to execute.
    //
    #ifdef _OPENMP
    omp_set_num_threads(num_threads);
    #endif

    initialize_dictionaries(global_dictionary);

    MetaPopInfo mpopi;
    mpopi.init_directory(in_path);

    if (wl_file.length() > 0) {
        load_marker_list(wl_file, whitelist);
        cerr << "Loaded " << whitelist.size() << " whitelisted markers.\n";
    }
    if (bl_file.length() > 0) {
        load_marker_list(bl_file, blacklist);
        cerr << "Loaded " << blacklist.size() << " blacklisted markers.\n";
    }

    //
    // Load the catalog
    //
    stringstream catalog_file;
    map<int, CSLocus *> catalog;
    bool compressed = false;
    catalog_file << in_path << "batch_" << batch_id << ".catalog";
    int res = load_loci(catalog_file.str(), catalog, 0, false, compressed);
    if (res == 0) {
        cerr << "Error: Unable to load the catalog '" << catalog_file.str() << "'\n";
        throw exception();
    }

    //
    // Implement the black/white list
    //
    reduce_catalog(catalog, whitelist, blacklist);

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

    // Set the deprecated [samples] variable to the value it should have
    map<int, string> samples; // map of {sample_id, sample_name}
    mpopi.fill_samples(samples);

    // Retrieve the IDs of the parents
    vector<int> sample_ids;
    for (auto& sample : mpopi.samples())
        sample_ids.push_back(sample.id);
    sort(sample_ids.begin(), sample_ids.end());
    set<int> parent_ids;
    identify_parental_ids(catalog, sample_ids, parent_ids);

    //
    // Create the population map
    //
    cerr << "Populating observed haplotypes for " << mpopi.samples().size() << " samples, " << catalog.size() << " loci.\n";
    PopMap<CSLocus> *pmap = new PopMap<CSLocus>(mpopi, catalog.size());
    pmap->populate(catalog, catalog_matches);

    apply_locus_constraints(catalog, pmap);

    //
    // Identify mappable markers in the parents
    //
    find_markers(catalog, pmap, parent_ids);

    //
    // Calculate F, inbreeding coefficient
    //
    calculate_f(catalog, pmap, parent_ids);

    //
    // Create genotypes maps, calculate mean log likelihood for each locus.
    //
    map<int, CSLocus *>::iterator it;
    Datum  **d;
    CSLocus *loc;
    double   mean, cnt;
    for (it = catalog.begin(); it != catalog.end(); it++) {
        loc = it->second;
        d   = pmap->locus(loc->id);

        if (loc->marker.length() > 0) {
            create_genotype_map(loc, pmap, parent_ids);
            call_population_genotypes(loc, pmap, global_dictionary);
        }

        mean = 0.0;
        cnt  = 0.0;
        for (int i = 0; i < pmap->sample_cnt(); i++) {
            if (d[i] == NULL)
                continue;
            mean += d[i]->lnl;
            cnt++;
        }
        loc->lnl = mean / cnt;
    }

    //
    // Make automated corrections
    //
    if (corrections)
        automated_corrections(samples, parent_ids, catalog, catalog_matches, pmap);

    //
    // Check markers for potentially missing alleles.
    //
    switch(map_type) {
    case cp:
        correct_cp_markers_missing_alleles(parent_ids, catalog, pmap);
        break;
    case dh:
    case bc1:
    case f2:
    case gen:
    case none:
    case unk:
        break;
    }

    //
    // Reassign genotypes according to specific map type, record any
    // marker corrections made by detecting missing alleles.
    //
    if (map_type != gen)
        map_specific_genotypes(catalog, pmap, parent_ids);

    //
    // Incorporate manual corrections exported from a Stacks SQL database.
    //
    if (man_corrections)
        manual_corrections(cor_path, pmap);

    //
    // Calculate segregation distortion using chi-squared test.
    //
    map<string, map<string, double> > seg_ratios;

    if (map_type != gen) {
        load_segregation_ratios(map_type, seg_ratios);
        calc_segregation_distortion(seg_ratios, catalog, pmap, parent_ids);
    }

    //
    // If a mapping type was specified, output it.
    //
    switch(map_type) {
    case dh:
        export_dh_map(catalog, pmap, parent_ids, samples);
        break;
    case cp:
        export_cp_map(catalog, pmap, parent_ids, samples);
        break;
    case bc1:
        export_bc1_map(catalog, pmap, parent_ids, samples);
        break;
    case f2:
        export_f2_map(catalog, pmap, parent_ids, samples);
        break;
    case gen:
        export_gen_map(catalog, pmap, parent_ids, samples);
        break;
    case none:
    case unk:
        break;
    }

    if (sql_out)
        write_sql(catalog, pmap, parent_ids);

    //
    // Output the observed haplotypes.
    //
    write_generic(catalog, pmap, samples, parent_ids, false);

    if (out_type == genomic)
        write_genomic(catalog, pmap);

    cerr << "genotypes is done.\n";
    return 0;
    IF_NDEBUG_CATCH_ALL_EXCEPTIONS
}

int
apply_locus_constraints(map<int, CSLocus *> &catalog,
                        PopMap<CSLocus> *pmap)
{
    CSLocus *loc;
    Datum  **d;
    uint below_stack_dep  = 0;
    uint below_lnl_thresh = 0;

    if (min_stack_depth == 0) return 0;

    map<int, CSLocus *>::iterator it;

    for (it = catalog.begin(); it != catalog.end(); it++) {
        loc = it->second;
        d   = pmap->locus(loc->id);

        for (int i = 0; i < pmap->sample_cnt(); i++) {
            //
            // Check that each sample is over the minimum stack depth for this locus.
            //
            if (d[i] != NULL &&
                d[i]->tot_depth < min_stack_depth) {
                below_stack_dep++;
                delete d[i];
                d[i] = NULL;
            }

            //
            // Check that each sample is over the log likelihood threshold.
            //
            if (d[i] != NULL &&
                filter_lnl   &&
                d[i]->lnl < lnl_limit) {
                below_lnl_thresh++;
                delete d[i];
                d[i] = NULL;
            }
        }
    }

    if (min_stack_depth > 0)
        cerr << "Removed " << below_stack_dep << " samples from loci that are below the minimum stack depth of " << min_stack_depth << "x\n";
    if (filter_lnl)
        cerr << "Removed " << below_lnl_thresh << " samples from loci that are below the log likelihood threshold of " << lnl_limit << "\n";

    return 0;
}

int identify_parental_ids(map<int, CSLocus *> &catalog, vector<int> &sample_ids, set<int> &parents) {
    set<int> catalog_parents;
    map<int, CSLocus *>::iterator it;
    CSLocus *loc;
    int      sample_id;

    //
    // We assume the catalog was constructed from one or more parents of one
    // or more crosses. These are listed in the catalog.tags.tsv file, column 8.
    //
    for (it = catalog.begin(); it != catalog.end(); it++) {
        loc = it->second;
        for (uint i = 0; i < loc->comp.size(); i++) {
            sample_id = (int) strtol(loc->comp[i], NULL, 10);
            catalog_parents.insert(sample_id);
        }
    }

    //
    // Now we want to iterate over those individuals genotypes found when
    // searching the Stacks directory and crosscheck those found in the catalog.
    //
    for (uint i = 0; i < sample_ids.size(); i++) {
        if (catalog_parents.count(sample_ids[i]) > 0)
            parents.insert(sample_ids[i]);
    }

    set<int>::iterator sit;
    cerr << "Identified parent IDs: ";
    for (sit = parents.begin(); sit != parents.end(); sit++)
        cerr << *sit << " ";
    cerr << "\n";

    return 0;
}

int find_markers(map<int, CSLocus *> &catalog, PopMap<CSLocus> *pmap, set<int> &parent_ids) {
    map<int, CSLocus *>::iterator it;
    vector<char *>::iterator hit;
    set<int>::iterator p, q;
    int pid_1, pid_2, parent_count, allele_cnt_1, allele_cnt_2;
    Datum   *d_1, *d_2;
    CSLocus *loc;

    if (parent_ids.size() > 2) return 0;

    for (it = catalog.begin(); it != catalog.end(); it++) {
        loc = it->second;
        //
        // Count the number of parental tags matching this catalog tag. A proper marker should
        // contain a single representative from each parent; multiple alleles must be called from
        // a single tag from a single parent.
        //
        if (parent_ids.size() == 1) {
            p = parent_ids.begin();
            pid_1 = *p;
            pid_2 = -1;
            d_1   = pmap->blacklisted(loc->id, pid_1) ? NULL : pmap->datum(loc->id, pid_1);
            d_2   = NULL;
        } else {
            p = parent_ids.begin();
            q = p++;

            pid_1 = *p < *q ? *p : *q;
            pid_2 = *p < *q ? *q : *p;
            if (pmap->blacklisted(loc->id, pid_1) ||
                pmap->blacklisted(loc->id, pid_2)) {
                d_1 = NULL;
                d_2 = NULL;
            } else {
                d_1   = pmap->datum(loc->id, pid_1);
                d_2   = pmap->datum(loc->id, pid_2);
            }
        }

        parent_count = 0;
        if (d_1 != NULL) parent_count++;
        if (d_2 != NULL) parent_count++;

        //
        // Locus is present in both parents.
        //
        if (parent_count == 2) {
            allele_cnt_1 = d_1->obshap.size();
            allele_cnt_2 = d_2->obshap.size();

            //
            // Determine the number of unique alleles
            //
            set<string> unique_alleles;

            for (hit = d_1->obshap.begin(); hit != d_1->obshap.end(); hit++)
                unique_alleles.insert(*hit);
            for (hit = d_2->obshap.begin(); hit != d_2->obshap.end(); hit++)
                unique_alleles.insert(*hit);
            int num_unique_alleles = unique_alleles.size();

            //
            // Locus is heterozygous in both parents. However, the number of alleles present distinguishes
            // what type of marker it is. Four unique alleles requries an ab/cd marker, while four
            // alleles that are the same in both parents requires an ab/ab marker. Finally, three unique
            // alleles requires either an ab/ac marker.
            //
            if (allele_cnt_1 == 2 && allele_cnt_2 == 2) {
                if (num_unique_alleles == 3)
                    loc->marker = "ab/ac";
                else if (num_unique_alleles == 2)
                    loc->marker = "ab/ab";
                else
                    loc->marker = "ab/cd";
            //
            // Locus is homozygous in one parent and heterozygous in the other.
            //
            } else if (allele_cnt_1 == 2 && allele_cnt_2 == 1) {

                if (num_unique_alleles == 3)
                    loc->marker = "ab/cc";
                else if (num_unique_alleles == 2)
                    loc->marker = "ab/aa";
            //
            // Locus is homozygous in one parent and heterozygous in the other.
            //
            } else if (allele_cnt_1 == 1 && allele_cnt_2 == 2) {

                if (num_unique_alleles == 3)
                    loc->marker = "cc/ab";
                else if (num_unique_alleles == 2)
                    loc->marker = "aa/ab";
            //
            // Locus is homozygous in both parents, but heterozygous between parents.
            //
            } else if (allele_cnt_1 == 1 && allele_cnt_2 == 1) {

                if (strcmp(d_1->obshap[0], d_2->obshap[0]) != 0)
                    loc->marker = "aa/bb";
            }
        //
        // Locus only exists in one parent.
        //
        } else if (parent_count == 1) {
            if (d_1 != NULL && d_1->obshap.size() == 2)
                loc->marker = "ab/--";
            else if (d_2 != NULL && d_2->obshap.size() == 2)
                loc->marker = "--/ab";
        }
    }

    return 0;
}

int calculate_f(map<int, CSLocus *> &catalog, PopMap<CSLocus> *pmap, set<int> &parent_ids) {
    map<int, CSLocus *>::iterator it;
    map<char, int>::iterator j;
    Datum  **d;
    CSLocus *loc;

    for (it = catalog.begin(); it != catalog.end(); it++) {
        loc = it->second;
        d   = pmap->locus(loc->id);

        if (loc->snps.size() == 0) continue;

        double tot  = 0.0;
        double hets = 0;
        double p, q, h, h0;
        map<char, int> alle;

        for (int i = 0; i < pmap->sample_cnt(); i++) {
            if (d[i] == NULL)
                continue;
            if (parent_ids.count(pmap->rev_sample_index(i)))
                continue;

            tot++;

            if (d[i]->obshap.size() > 1) hets++;

            //
            // We are measuring the first SNP in the haplotype
            //
            for (uint j = 0; j < d[i]->obshap.size(); j++)
                alle[d[i]->obshap[j][0]]++;
        }

        if (alle.size() > 2 || tot == 0)
            continue;

        j  = alle.begin();
        p  = j->second;
        j++;
        q  = j->second;
        h  = hets / tot;            // Observed frequency of heterozygotes in the population
        h0 = 2 * (p/tot) * (q/tot); // 2PQ, expected frequency of hets under Hardy-Weinberg
        if (h0 > 0)
            loc->f = (h0 - h) / h0;

        //cerr << "P: " << p << "; Q: " << q << "; Hets: " << hets << "; total: " << tot << "; f: " << loc->f << "\n";
    }

    return 0;
}

int create_genotype_map(CSLocus *locus, PopMap<CSLocus> *pmap, set<int> &parent_ids) {
    //
    // Create a genotype map. For any set of alleles, this routine will
    // assign each allele to one of the constituent genotypes, e.g. given the
    // marker type 'aaxbb' and the alleles 'A' from the male, and 'G'
    // from the female, will assign 'G' == 'bb' and 'A'== 'aa'. It assumes that
    // recombination may have occurred as with an F2, F3 or later cross.
    //
    //cerr << "Creating genotype map for catalog ID " << locus->id  << ", marker: " << locus->marker << ".\n";

    //
    // Get the parent IDs ordered
    //
    set<int>::iterator p = parent_ids.begin();
    set<int>::iterator q = p++;
    int pid_1 = *p < *q ? *p : *q;
    int pid_2 = *p < *q ? *q : *p;

    set<char> p1_gtypes, p2_gtypes;
    set<char>::iterator i;
    map<char, int> legal_gtypes, com_gtypes;
    //
    // First, identify any alleles that are common between the two parents.
    //
    p1_gtypes.insert(locus->marker[0]);
    p1_gtypes.insert(locus->marker[1]);
    p2_gtypes.insert(locus->marker[3]);
    p2_gtypes.insert(locus->marker[4]);
    for (i = p1_gtypes.begin(); i != p1_gtypes.end(); i++)
        if (*i != '-') legal_gtypes[*i]++;
    for (i = p2_gtypes.begin(); i != p2_gtypes.end(); i++)
        if (*i != '-') legal_gtypes[*i]++;
    //
    // Find the common genotypes
    //
    vector<char> types;
    map<char, int>::iterator j;

    for (j = legal_gtypes.begin(); j != legal_gtypes.end(); j++)
        if (j->second > 1) types.push_back(j->first);
    sort(types.begin(), types.end());

    Datum *d_1, *d_2;
    map<string, int> haplotypes;
    map<string, int>::iterator k;
    vector<pair<string, int> > sorted_haplotypes;
    d_1 = pmap->datum(locus->id, pid_1);
    d_2 = pmap->datum(locus->id, pid_2);

    if (d_1 != NULL) {
        for (uint n = 0; n < d_1->obshap.size(); n++)
            haplotypes[d_1->obshap[n]]++;
    }
    if (d_2 != NULL) {
        for (uint n = 0; n < d_2->obshap.size(); n++)
            haplotypes[d_2->obshap[n]]++;
    }

    //
    // Sort the haplotypes map by value
    //
    for (k = haplotypes.begin(); k != haplotypes.end(); k++)
        sorted_haplotypes.push_back(*k);
    sort(sorted_haplotypes.begin(), sorted_haplotypes.end(), hap_compare);

    for (uint n = 0, index = 0; n < sorted_haplotypes.size() && index < types.size(); n++, index++) {
        if (sorted_haplotypes[n].second > 1) {
            locus->gmap[sorted_haplotypes[n].first] = types[index];
            com_gtypes[types[index]]++;
            // cerr << "  Assigning common allele " << sorted_haplotypes[n].first << " to genotype '" << locus->gmap[sorted_haplotypes[n].first] << "'\n";
        }
    }

    //
    // Now, examine the remaining first parent alleles.
    //
    if (d_1 != NULL) {
        legal_gtypes.clear();
        for (i = p1_gtypes.begin(); i != p1_gtypes.end(); i++)
            if (*i != '-' && com_gtypes.count(*i) == 0) {
                // cerr << "  Adding " << *i << " to first parent genotypes\n";
                legal_gtypes[*i]++;
            }
        types.clear();
        for (j = legal_gtypes.begin(); j != legal_gtypes.end(); j++)
            types.push_back(j->first);
        sort(types.begin(), types.end());

        for (uint n = 0, index = 0; n < d_1->obshap.size() && index < types.size(); n++, index++) {
            if (locus->gmap.count(d_1->obshap[n])) {
                index--;
                continue;
            }
            locus->gmap[d_1->obshap[n]] = types[index];
            // cerr << "  Assinging '" << d_1->obshap[n] << "' to first parent genotype '" << locus->gmap[d_1->obshap[n]] << "'\n";
        }
    }

    //
    // Finally, repeat in the second parent.
    //
    if (d_2 != NULL) {
        legal_gtypes.clear();
        for (i = p2_gtypes.begin(); i != p2_gtypes.end(); i++)
            if (*i != '-' && com_gtypes.count(*i) == 0) {
                // cerr << "  Adding " << *i << " to second genotypes\n";
                legal_gtypes[*i]++;
            }
        types.clear();
        for (j = legal_gtypes.begin(); j != legal_gtypes.end(); j++)
            types.push_back(j->first);
        sort(types.begin(), types.end());

        for (uint n = 0, index = 0; n < d_2->obshap.size() && index < types.size(); n++, index++) {
            if (locus->gmap.count(d_2->obshap[n])) {
                index--;
                continue;
            }
            locus->gmap[d_2->obshap[n]] = types[index];
            // cerr << "  Assinging '" << d_2->obshap[n] << "' to second parent genotype '" << locus->gmap[d_2->obshap[n]] << "'\n";
        }
    }

    return 0;
}

int call_population_genotypes(CSLocus *locus,
                              PopMap<CSLocus> *pmap,
                              map<string, map<string, string> > &dictionary) {
    //
    // Fetch the array of observed haplotypes from the population
    //
    Datum **d = pmap->locus(locus->id);

    for (int i = 0; i < pmap->sample_cnt(); i++) {
        if (d[i] == NULL)
            continue;

        vector<string> gtypes;
        string gtype;

        //cerr << "Sample Id: " << pmap->rev_sample_index(i) << "\n";

        for (uint j = 0; j < d[i]->obshap.size(); j++) {
            //
            // Impossible allele encountered.
            //
            if (locus->gmap.count(d[i]->obshap[j]) == 0) {
                gtypes.clear();
                gtypes.push_back("-");
                goto impossible;
            }

            gtypes.push_back(locus->gmap[d[i]->obshap[j]]);
            //cerr << "  Observed Haplotype: " << d[i]->obshap[j] << ", Genotype: " << locus->gmap[d[i]->obshap[j]] << "\n";
        }

    impossible:
        sort(gtypes.begin(), gtypes.end());
        for (uint j = 0; j < gtypes.size(); j++) {
            gtype += gtypes[j];
            //cerr << "  Adding genotype to string: " << gtypes[j] << "; " << gtype << "\n";
        }

        string m = dictionary[locus->marker].count(gtype) ?
            dictionary[locus->marker][gtype] :
            "--";

        if (d[i]->gtype != NULL)
            delete d[i]->gtype;

        d[i]->gtype = new char[m.length() + 1];
        strcpy(d[i]->gtype, m.c_str());

        if (m != "--")
            locus->gcnt++;

        //cerr << "Assigning datum, marker: " << locus->marker << ", string: " << m << ", haplotype: " << d[i]->obshap[0] << ", gtype: " << gtype << "\n";
     }

    return 0;
}

int
correct_cp_markers_missing_alleles(set<int> &parent_ids, map<int, CSLocus *> &catalog, PopMap<CSLocus> *pmap)
{
    map<int, CSLocus *>::iterator it;
    CSLocus *loc;
    Datum  **d;
    map<string, map<string, double> > seg_ratios;

    //
    // The following segregation ratios will occur when one of the parents
    // in the cross is missing an allele. We will not see these ratios
    // in one of these markers with no missing alleles.
    //
    seg_ratios["ab/aa"]["aa"] = 0.50;
    seg_ratios["ab/aa"]["ab"] = 0.25;
    seg_ratios["ab/aa"]["bb"] = 0.25;

    seg_ratios["aa/ab"]["aa"] = 0.50;
    seg_ratios["aa/ab"]["ab"] = 0.25;
    seg_ratios["aa/ab"]["bb"] = 0.25;

    seg_ratios["ab/cc"]["ac"] = 0.25;
    seg_ratios["ab/cc"]["bc"] = 0.25;
    seg_ratios["ab/cc"]["aa"] = 0.25;
    seg_ratios["ab/cc"]["bb"] = 0.25;

    seg_ratios["cc/ab"]["ac"] = 0.25;
    seg_ratios["cc/ab"]["bc"] = 0.25;
    seg_ratios["cc/ab"]["aa"] = 0.25;
    seg_ratios["cc/ab"]["bb"] = 0.25;

    cerr << "Testing catalog loci for mapping parents with missing alleles...";
    int corrected = 0;

    for (it = catalog.begin(); it != catalog.end(); it++) {
        loc = it->second;

        //
        // We only want to examine markers where one parent is homozygous.
        //
        if (loc->marker != "ab/aa" &&
            loc->marker != "aa/ab" &&
            loc->marker != "ab/cc" &&
            loc->marker != "cc/ab") continue;

        map<string, int> cnts;

        //
        // Calculate initial segregation distortion.
        //
        double n          = tally_generic_gtypes(loc->id, pmap, parent_ids, cnts);
        double chisq_pval = chisq_test(seg_ratios, cnts, loc->marker, n);

        //
        // Check if our genotype ratios match the segregation ratios specified above. If so,
        // we have a dropped allele in one of the parents.
        //
        if (n == 0 || chisq_pval < chisq_pval_limit)
            continue;

        corrected++;

        if (loc->marker == "ab/aa")
            loc->marker = "ab/a-";
        else if (loc->marker == "aa/ab")
            loc->marker = "-a/ab";
        else if (loc->marker == "ab/cc")
            loc->marker = "ab/c-";
        else if (loc->marker == "cc/ab")
            loc->marker = "-c/ab";

        d = pmap->locus(loc->id);

        if (loc->marker == "ab/a-" || loc->marker == "-a/ab") {

            for (int i = 0; i < pmap->sample_cnt(); i++) {
                if (d[i] == NULL) continue;

                if (parent_ids.count(pmap->rev_sample_index(i))) continue;

                if (strcmp(d[i]->gtype, "bb") == 0)
                    strcpy(d[i]->gtype, "ab");
            }

        } else if (loc->marker == "ab/c-") {

            for (int i = 0; i < pmap->sample_cnt(); i++) {
                if (d[i] == NULL) continue;

                if (parent_ids.count(pmap->rev_sample_index(i))) continue;

                if (strcmp(d[i]->gtype, "bb") == 0)
                    strcpy(d[i]->gtype, "bd");
                else if (strcmp(d[i]->gtype, "aa") == 0)
                    strcpy(d[i]->gtype, "ad");
            }

        } else if (loc->marker == "-c/ab") {

            for (int i = 0; i < pmap->sample_cnt(); i++) {
                if (d[i] == NULL) continue;

                if (parent_ids.count(pmap->rev_sample_index(i))) continue;

                if (strcmp(d[i]->gtype, "bb") == 0)
                    strcpy(d[i]->gtype, "ad");
                else if (strcmp(d[i]->gtype, "aa") == 0)
                    strcpy(d[i]->gtype, "ac");
                else if (strcmp(d[i]->gtype, "bc") == 0)
                    strcpy(d[i]->gtype, "bd");
                else if (strcmp(d[i]->gtype, "ac") == 0)
                    strcpy(d[i]->gtype, "bc");
            }
        }
    }

    //
    // Now we will deal with aa/bb markers separately, since there can be three possible
    // missing allele situations:
    //   aa: 50%, ab: 50% - we have an aa/b- marker, which should be mapped as an --/ab
    //   bb: 50%, ab: 50% - we have an -a/bb marker, which should be mapped as an ab/--
    //   aa: 33%, ab: 33%, bb: 33% - we have an -a/b- maker, which should be mapped as an ab/ab, but
    //                               we can't disambiguate the aa bb genotypes so it can't be mapped.
    //
    map<string, map<string, double> > seg_ratio_1, seg_ratio_2;
    seg_ratio_1["aa/bb"]["aa"] = 0.50;
    seg_ratio_1["aa/bb"]["ab"] = 0.50;

    seg_ratio_2["aa/bb"]["bb"] = 0.50;
    seg_ratio_2["aa/bb"]["ab"] = 0.50;

    for (it = catalog.begin(); it != catalog.end(); it++) {
        loc = it->second;

        if (loc->marker != "aa/bb") continue;

        map<string, int> cnts;

        double n          = tally_generic_gtypes(loc->id, pmap, parent_ids, cnts);
        double chisq_pval = chisq_test(seg_ratio_1, cnts, loc->marker, n);

        if (n == 0) continue;

        if (chisq_pval >= chisq_pval_limit) {
            corrected++;

            loc->marker = "aa/b-";

            d = pmap->locus(loc->id);
            for (int i = 0; i < pmap->sample_cnt(); i++) {
                if (d[i] == NULL) continue;

                if (parent_ids.count(pmap->rev_sample_index(i))) continue;

                if (strcmp(d[i]->gtype, "ab") == 0)
                    strcpy(d[i]->gtype, "bb");
            }

        } else {
            chisq_pval = chisq_test(seg_ratio_2, cnts, loc->marker, n);

            if (chisq_pval >= chisq_pval_limit) {
                corrected++;
                loc->marker = "-a/bb";

                d = pmap->locus(loc->id);
                for (int i = 0; i < pmap->sample_cnt(); i++) {
                    if (d[i] == NULL) continue;

                    if (parent_ids.count(pmap->rev_sample_index(i))) continue;

                    if (strcmp(d[i]->gtype, "ab") == 0)
                        strcpy(d[i]->gtype, "aa");
                }

            }
        }
    }

    cerr << "corrected " << corrected << " catalog loci.\n";

    return 0;
}

int
automated_corrections(map<int, string> &samples, set<int> &parent_ids, map<int, CSLocus *> &catalog,
                      vector<vector<CatMatch *> > &matches, PopMap<CSLocus> *pmap)
{
    int sample_id, catalog_id, tag_id;
    Datum *d;
    Locus *s;
    string file;

    cerr << "Performing automated corrections...\n";

    for (uint i = 0; i < matches.size(); i++) {
        sample_id = matches[i][0]->sample_id;
        file      = samples[sample_id];

        //if (sample_id != 29) continue;
        if (parent_ids.count(sample_id)) continue;

        map<int, Locus *> stacks;
        bool compressed = false;
        int  res;
        if ((res = load_loci(in_path + file, stacks, 2, false, compressed)) == 0) {
            cerr << "Unable to load sample file '" << file << "'\n";
            return 0;
        }

        set<pair<int, int> > processed;

        for (uint j = 0; j < matches[i].size(); j++) {
             catalog_id = matches[i][j]->cat_id;
             sample_id  = matches[i][j]->sample_id;
             tag_id     = matches[i][j]->tag_id;

             if (catalog.count(catalog_id) == 0) continue;

             //
             // There are multiple matches per stack, but we only need to process
             // each stack once to make corrections.
             //
             if (processed.count(make_pair(catalog_id, tag_id)) == 0 &&
                 catalog[catalog_id]->marker.length() > 0) {
                 d = pmap->datum(catalog_id, sample_id);

                 //cerr << "Accessing catalog ID " << catalog_id << "; sample: " << sample_id << "; marker: " << catalog[catalog_id]->marker << ": d: " << d << "; gtype: " << d->gtype << "\n";

                 if (d != NULL && strcmp(d->gtype, "--") != 0) {
                     s = stacks[tag_id];
                     check_uncalled_snps(catalog[catalog_id], s, d);
                 }

                 processed.insert(make_pair(catalog_id, tag_id));
             }
        }

        //
        // Free up memory
        //
        map<int, Locus *>::iterator it;
        for (it = stacks.begin(); it != stacks.end(); it++)
            delete it->second;
    }

    //
    // Summarize correction results
    //
    long pot_gen = 0;
    long tot_gen = 0;
    long tot_cor = 0;
    long het_cor = 0;
    long rem_cor = 0;
    int  markers = 0;
    map<int, CSLocus *>::iterator it;
    for (it = catalog.begin(); it != catalog.end(); it++) {
        if (it->second->marker.length() == 0) continue;
        markers++;

        for (int j = 0; j < pmap->sample_cnt(); j++)  {
            sample_id = pmap->rev_sample_index(j);
            if (parent_ids.count(sample_id)) continue;

            d = pmap->datum(it->first, sample_id);

            pot_gen++;

            if (d == NULL) continue;

            tot_gen++;

            if (d->corrected == true) {
                tot_cor++;

                if (strcmp(d->gtype, "--") == 0)
                    rem_cor++;
                else
                    het_cor++;
            }
        }
    }
    cerr << pot_gen << " potential genotypes in " << markers << " markers, "
         << tot_gen << " populated; "
         << tot_cor << " corrected, "
         << het_cor << " converted to heterozygotes, "
         << rem_cor << " unsupported homozygotes removed.\n";

    return 0;
}

int check_uncalled_snps(CSLocus *clocus, Locus *stack, Datum *d) {
    //
    // If this locus is already known to be multi-allele, return, we only want
    // to verify uncalled SNPs.
    //
    if (strlen(d->gtype) > 1 && d->gtype[0] != d->gtype[1])
        return 0;

    //cerr << "Catalog locus: " << clocus->id << ", marker: " << clocus->marker << ", tag_id: " << stack->id << "; Starting gtype: " << d->gtype << "\n";

    vector<SNP *> verified_snps;
    string status = "false";
    string homozygous;

    for (uint i = 0; i < clocus->snps.size(); i++) {
        check_homozygosity(stack->reads,
                           clocus->snps[i]->col, clocus->snps[i]->rank_1, clocus->snps[i]->rank_2,
                           homozygous);

        if (homozygous == "unknown")
            status = "true";
        else if (homozygous == "false")
            verified_snps.push_back(clocus->snps[i]);
     }

    if (status == "true") {
        d->corrected = true;
        delete [] d->gtype;
        d->gtype = new char[2];
        strcpy(d->gtype, "-");
        return 0;
    } else if (verified_snps.size() < clocus->snps.size()) {
        return 0;
    }

    //
    // Detect the alleles present from the verified SNPs
    //
    vector<string> haplotypes;
    call_alleles(verified_snps, stack->reads, haplotypes);

    vector<string> types;
    for (uint i = 0; i < haplotypes.size(); i++) {
        if (clocus->gmap.count(haplotypes[i])) {
            //cerr << "    Adding haplotype '" << haplotypes[i] << "', which maps to '" << clocus->gmap[haplotypes[i]] << "' to the genotype\n";
            types.push_back(clocus->gmap[haplotypes[i]]);
        } else {
            //cerr << "    Adding haplotype '-' for " << haplotypes[i] << "\n";
            types.push_back("-");
        }
    }

    sort(types.begin(), types.end());

    string genotype;
    for (uint i = 0; i < types.size(); i++)
        genotype += types[i];

    //cerr << "Final genotype: " << genotype << "\n";

    genotype =
        global_dictionary[clocus->marker].count(genotype) ?
        global_dictionary[clocus->marker][genotype] :
        "--";

    //cerr << "Final translated genotype: " << genotype << "\n";

    if (strcmp(genotype.c_str(), d->gtype) != 0) {
        d->corrected = true;
        delete [] d->gtype;
        d->gtype = new char[genotype.length() + 1];
        strcpy(d->gtype, genotype.c_str());
    }
    //cerr << "  Catalog: " << clocus->id << ", stack: " << stack->id << ", Ending Genotype: " << d->gtype << "\n\n";

    return 0;
}

int call_alleles(vector<SNP *> &snps, vector<char *> &reads, vector<string> &haplotypes) {
    map<string, int> a;
    int  height = reads.size();
    char base;

    for (int i = 0; i < height; i++) {
        string haplotype;

        for (uint j = 0; j < snps.size(); j++) {
            base = reads[i][snps[j]->col];

            //
            // Check to make sure the nucleotide at the location of this SNP is
            // of one of the two possible states the multinomial model called.
            //
            if (base == snps[j]->rank_1 || base == snps[j]->rank_2)
                haplotype += base;
            else
                break;
        }

        if (haplotype.length() == snps.size())
            a[haplotype]++;
    }

    map<string, int>::iterator it;
    for (it = a.begin(); it != a.end(); it++) {
        //cerr << "    Calling haplotype: " << it->first << "\n";
        haplotypes.push_back(it->first);
    }

    return 0;
}

int check_homozygosity(vector<char *> &reads, int col, char rank_1, char rank_2, string &homozygous) {

    //cerr << "  Examining col " << col << ", rank 1: " << rank_1 << "; rank 2: " << rank_2 << "\n";

    int height = reads.size();
    homozygous = "true";

    if (height < min_hom_seqs) {
        homozygous = "unknown";
        return 0;
    }

    map<char, int> nuc;
    vector<pair<char, int> > sorted_nuc;

    nuc['A'] = 0;
    nuc['C'] = 0;
    nuc['G'] = 0;
    nuc['T'] = 0;

    for (int j = 0; j < height; j++)
        nuc[reads[j][col]]++;

    map<char, int>::iterator i;
    for (i = nuc.begin(); i != nuc.end(); i++)
        sorted_nuc.push_back(make_pair(i->first, i->second));

    sort(sorted_nuc.begin(), sorted_nuc.end(), compare_pair);

    //
    // Check if more than a single nucleotide occurs in this column. Only
    // count nucleotides that are part of the called SNP, do not count
    // error-generated nucleotides. Also, check that the sorting was successful
    // by ensuring that sorted_nuc[0] > sorted_nuc[1] > sorted_nuc[2].
    //
    if (sorted_nuc[2].second > 0 && sorted_nuc[1].second <= sorted_nuc[2].second) {
        homozygous = "unknown";
        return 0;
    }

    // cerr << "Sorted_nuc[0], '" << sorted_nuc[0].first << "', count: " << sorted_nuc[0].second
    //          << "; Sorted_nuc[1], '" << sorted_nuc[1].first << "', count: " << sorted_nuc[1].second
    //          << "; Sorted_nuc[2], '" << sorted_nuc[2].first << "', count: " << sorted_nuc[2].second << "\n";

    if ((sorted_nuc[0].second > 0) &&
        (sorted_nuc[0].first == rank_1 || sorted_nuc[0].first == rank_2) &&
        (sorted_nuc[1].second > 0) &&
        (sorted_nuc[1].first == rank_1 || sorted_nuc[1].first == rank_2)) {
        homozygous = "false";
    }

    //
    // If we find a second nucleotide present, check its prevelence. If it is
    // less than 1/20 of the total reads, don't count a heterozygote. If it is
    // less than 1/10 report that we can't tell if its homozygous or not. Otherwise,
    // report this tag as a heterozygote.
    //
    double frac = (double) sorted_nuc[1].second / (double) height;
    //cerr << "    Frac: " << frac << "; Second-most Prominent Nuc count: " << sorted_nuc[1].second << "; Depth: " << height << "\n";

    if (homozygous == "false" && frac < min_het_seqs)
        homozygous = "true";
    else if (homozygous == "false" && frac < max_het_seqs)
        homozygous = "unknown";

    //cerr << "      Homozygous: " << homozygous << "\n";

    return 0;
}

int
manual_corrections(string cor_path, PopMap<CSLocus> *pmap)
{
    //
    // Load manual corrections from a tab-seprated file, as exported from a Stacks SQL
    // dataabse. Has the format:
    //   id<tab>batch_id<tab>catalog_id<tab>sample_id<tab>genotype
    //
    char     line[max_len];
    ifstream fh(cor_path.c_str(), ifstream::in);

    if (fh.fail()) {
        cerr << "Error opening manual corrections file '" << cor_path << "'\n";
        exit(1);
    }

    vector<string> parts;
    int   catalog_id, sample_id, len;
    char  gtype[id_len];
    char *e;

    int   line_num = 0;
    int   total    = 0;
    int   skipped  = 0;
    int   success  = 0;
    int   i;

    while (fh.good()) {
        fh.getline(line, id_len);
        line_num++;

        len = strlen(line);
        if (len == 0) continue;

        //
        // Check that there is no carraige return in the buffer.
        //
        if (line[len - 1] == '\r') line[len - 1] = '\0';

        //
        // Ignore comments
        //
        if (line[0] == '#') continue;

        parse_tsv(line, parts);

        if (parts.size() != 5) {
            cerr << "Error parsing '" << line << "' at line: " << line_num << ". (" << parts.size() << " fields).\n";
            return 0;
        }

        catalog_id = (int) strtol(parts[2].c_str(), &e, 10);
        if (*e != '\0') {
            cerr << "Error parsing '" << parts[2].c_str() << "' at line: " << line_num << ".\n";
            return 0;
        }
        sample_id = (int) strtol(parts[3].c_str(), &e, 10);
        if (*e != '\0') {
            cerr << "Error parsing '" << parts[3].c_str() << "' at line: " << line_num << ".\n";
            return 0;
        }
        strcpy(gtype, parts[4].c_str());

        //
        // Overwrite the genotype in the PopMap.
        //
        Datum **d = pmap->locus(catalog_id);
        total++;

        if (d == NULL) {
            skipped++;
            continue;
        }

        if ((i = pmap->sample_index(sample_id)) < 0) {
            skipped++;
            continue;
        }

        if (d[i] == NULL) {
            skipped++;
            continue;
        }

        for (uint k = 0; k < strlen(gtype); k++)
            gtype[k] = tolower(gtype[k]);

        if (strcmp(gtype, "--") == 0)
            strcpy(gtype, "-");

        if (d[i]->gtype != NULL)
            delete [] d[i]->gtype;

        d[i]->gtype = new char[strlen(gtype) + 1];
        strcpy(d[i]->gtype, gtype);

        d[i]->corrected = true;

        success++;
    }

    fh.close();

    cerr << "Successfully imported "
         << success << " manually corrected genotypes. Skipped "
         << skipped << " genotypes due to invalid catalog or sample IDs, "
         << total << " genotypes read from file.\n";

    return 0;
}

int export_gen_map(map<int, CSLocus *> &catalog, PopMap<CSLocus> *pmap, set<int> &parent_ids, map<int, string> &samples) {
    //
    // We wish to export, a set of generic genotypes, not specific to any mapping type.
    //

    //
    // Mark those genotypes that have been corrected in uppercase letters.
    //
    map<int, CSLocus *>::iterator it;
    CSLocus *loc;
    uint     len;

    for (it = catalog.begin(); it != catalog.end(); it++) {
        loc = it->second;

        Datum **d = pmap->locus(loc->id);

        for (int i = 0; i < pmap->sample_cnt(); i++) {
            if (d[i] == NULL) continue;

            if (parent_ids.count(pmap->rev_sample_index(i))) continue;

            if (d[i]->corrected) {
                len = strlen(d[i]->gtype);
                for (uint k = 0; k < len; k++)
                    d[i]->gtype[k] = toupper(d[i]->gtype[k]);
            }
        }
    }

    //
    // Output the results
    //
    write_generic(catalog, pmap, samples, parent_ids, true);

    return 0;
}

int export_f2_map(map<int, CSLocus *> &catalog, PopMap<CSLocus> *pmap, set<int> &parent_ids, map<int, string> &samples) {
    //
    // We wish to export, according to the JoinMap manual, a locus genotype file (loc-file),
    // which contains the information of all the loci for a single segregating population.
    //
    // We are exporting an F2 population type:
    // The result of selfing the F1 of a cross between two fully homozygous diploid parents.
    //
    // Genotype codes for an F2 population, depending on the locus segregation type.
    //
    // Seg. type   Possible genotypes
    // ---------   ------------------
    // <aaxbb>     a, b, h, –
    // <abxcd>     a, b, h, –
    // <abxaa>     a, –
    // <aaxab>     b, –
    // <abxcc>     a, b, –
    // <ccxab>     b, –
    //
    map<string, string> types;
    map<string, map<string, string> > dictionary;

    load_mm_f2_dictionary(types, dictionary);

    //
    // Translate the genotypes for this particular map type.
    //
    translate_genotypes(types, dictionary, catalog, pmap, samples, parent_ids);

    //
    // Output the results
    //
    switch(out_type) {
    case joinmap:
        write_joinmap(catalog, pmap, types, samples, parent_ids);
        break;
    case rqtl:
        write_rqtl(catalog, pmap, types, samples, parent_ids);
        break;
    case onemap:
        write_onemap_mapmaker(catalog, pmap, types, samples, parent_ids);
        break;
    default:
        break;
    }

    return 0;
}

int export_dh_map(map<int, CSLocus *> &catalog, PopMap<CSLocus> *pmap, set<int> &parent_ids, map<int, string> &samples) {
    //
    // We wish to export, according to the JoinMap manual, a locus genotype file (loc-file),
    // which contains the information of all the loci for a single segregating population.
    //
    // We are exporting a DH population type:
    //   a doubled haploid population: the result of doubling the gametes of a single heterozygous
    //   diploid individual.
    //
    // Segregation type codes for population type DH, from Joinmap manual:
    //
    // Code      Description
    // -------   -----------
    // <aaxbb>   locus homozygous in both parents, heterozygous between the parents
    //
    // Genotype codes for a CP population, depending on the locus segregation type.
    //
    // Seg. type   Possible genotypes
    // ---------   ------------------
    //     a       the one genotype
    //     b       the other genotype
    //
    map<string, string> types;
    map<string, map<string, string> > dictionary;

    load_mm_dh_dictionary(types, dictionary);

    //
    // Translate the genotypes for this particular map type.
    //
    translate_genotypes(types, dictionary, catalog, pmap, samples, parent_ids);

    //
    // Output the results
    //
    switch(out_type) {
    case joinmap:
        write_joinmap(catalog, pmap, types, samples, parent_ids);
        break;
    case rqtl:
        write_rqtl(catalog, pmap, types, samples, parent_ids);
        break;
    default:
        break;
    }

    return 0;
}

int export_bc1_map(map<int, CSLocus *> &catalog, PopMap<CSLocus> *pmap, set<int> &parent_ids, map<int, string> &samples) {
    //
    // We wish to export, according to the JoinMap manual, a locus genotype file (loc-file),
    // which contains the information of all the loci for a single segregating population.
    //
    // We are exporting a BC1 population type:
    //   a first generation backcross population: the result of crossing the F1 of a cross between
    //   two fully homozygous diploid parents to one of the parents.
    //
    // Segregation type codes for population type BC1, from Joinmap manual:
    //
    // Code      Description
    // -------   -----------
    // <aaxbb>   locus homozygous in both parents, heterozygous between the parents
    //
    // Genotype codes for a BC1 population, depending on the locus segregation type.
    //
    // Seg. type   Possible genotypes
    // ---------   ------------------
    //     a       homozygote or haploid as the first parent
    //     b       homozygote or haploid as the second parent
    //     h       heterozygote (as the F1)
    //
    map<string, string> types;
    map<string, map<string, string> > dictionary;

    load_mm_bc_dictionary(types, dictionary);

    //
    // Translate the genotypes for this particular map type.
    //
    translate_genotypes(types, dictionary, catalog, pmap, samples, parent_ids);

    //
    // Output the results
    //
    switch(out_type) {
    case joinmap:
        write_joinmap(catalog, pmap, types, samples, parent_ids);
        break;
    case rqtl:
        write_rqtl(catalog, pmap, types, samples, parent_ids);
        break;
    case onemap:
        write_onemap_mapmaker(catalog, pmap, types, samples, parent_ids);
        break;
    default:
        break;
    }

    return 0;
}

int export_cp_map(map<int, CSLocus *> &catalog, PopMap<CSLocus> *pmap, set<int> &parent_ids, map<int, string> &samples) {
    //
    // We wish to export, according to the JoinMap manual, a locus genotype file (loc-file),
    // which contains the information of all the loci for a single segregating population.
    //
    // We are exporting a CP population type:
    //   a population resulting from a cross between two heterogeneously
    //   heterozygous and homozygous diploid parents, linkage phases originally
    //   (possibly) unknown.
    //
    // Segregation type codes for population type CP, from Joinmap manual:
    //
    // Code      Description
    // -------   -----------
    // <abxcd>   locus heterozygous in both parents, four alleles
    // <efxeg>   locus heterozygous in both parents, three alleles
    // <hkxhk>   locus heterozygous in both parents, two alleles
    // <lmxll>   locus heterozygous in the first parent
    // <nnxnp>   locus heterozygous in the second parent
    //
    // Genotype codes for a CP population, depending on the locus segregation type.
    //
    // Seg. type   Possible genotypes
    // ---------   ------------------
    // <abxcd>     ac, ad, bc, bd, ––
    // <efxeg>     ee, ef, eg, fg, ––
    // <hkxhk>     hh, hk, kk, h-, k-, ––
    // <lmxll>     ll, lm, ––
    // <nnxnp>     nn, np, ––
    //
    map<string, string> types;
    map<string, map<string, string> > dictionary;

    switch(out_type) {
    case joinmap:
        load_joinmap_cp_dictionary(types, dictionary);
        break;
    case onemap:
        load_onemap_cp_dictionary(types, dictionary);
        break;
    default:
        break;
    }

    //
    // Translate the genotypes for this particular map type.
    //
    translate_genotypes(types, dictionary, catalog, pmap, samples, parent_ids);

    //
    // Output the results
    //
    switch(out_type) {
    case joinmap:
        write_joinmap(catalog, pmap, types, samples, parent_ids);
        break;
    case onemap:
        write_onemap(catalog, pmap, types, samples, parent_ids);
        break;
    default:
        break;
    }

    return 0;
}

int
calc_segregation_distortion(map<string, map<string, double> > &seg_ratios, map<int, CSLocus *> &catalog,
                            PopMap<CSLocus> *pmap, set<int> &parent_ids)
{
    map<string, string>               types;
    map<string, map<string, string> > dictionary;

    switch(map_type) {
    case dh:
        load_dh_dictionary(types, dictionary);
        break;
    case cp:
        load_cp_dictionary(types, dictionary);
        break;
    case bc1:
        load_mm_bc_dictionary(types, dictionary);
        break;
    case f2:
        load_mm_f2_dictionary(types, dictionary);
        break;
    case gen:
    case none:
    case unk:
        break;
    }

    map<int, CSLocus *>::iterator it;
    CSLocus *loc;

    for (it = catalog.begin(); it != catalog.end(); it++) {
        loc = it->second;

        if (seg_ratios.count(loc->marker) == 0) continue;

        map<string, int> cnts;

        double n = tally_translated_gtypes(loc->id, pmap, parent_ids, dictionary[loc->marker], cnts);

        if (n == 0) continue;

        // cerr << "ID: " << loc->id << "; marker: " << loc->marker << "\n";

        loc->chisq = chisq_test(seg_ratios, cnts, loc->marker, n);
    }

    return 0;
}

double
tally_translated_gtypes(int loc_id, PopMap<CSLocus> *pmap, set<int> &parent_ids,
                        map<string, string> &dictionary, map<string, int> &cnts)
{
    Datum **d = pmap->locus(loc_id);
    double  n = 0.0;

    for (int i = 0; i < pmap->sample_cnt(); i++) {
        if (d[i] == NULL) continue;

        if (parent_ids.count(pmap->rev_sample_index(i))) continue;

        if (strcmp(d[i]->gtype, "--") == 0)
            continue;

        n++;

        if (cnts.count(dictionary[d[i]->gtype]) > 0)
            cnts[dictionary[d[i]->gtype]]++;
        else
            cnts[dictionary[d[i]->gtype]] = 1;
    }

    return n;
}

double
tally_generic_gtypes(int loc_id, PopMap<CSLocus> *pmap, set<int> &parent_ids, map<string, int> &cnts)
{
    Datum **d = pmap->locus(loc_id);
    double  n = 0.0;

    for (int i = 0; i < pmap->sample_cnt(); i++) {
        if (d[i] == NULL) continue;

        if (parent_ids.count(pmap->rev_sample_index(i))) continue;

        if (strcmp(d[i]->gtype, "--") == 0)
            continue;

        n++;

        if (cnts.count(d[i]->gtype) > 0)
            cnts[d[i]->gtype]++;
        else
            cnts[d[i]->gtype] = 1;
    }

    return n;
}

double
chisq_test(map<string, map<string, double> > &seg_ratios, map<string, int> &cnts, string marker, double n)
{
    //
    // Calculate chi-square value.
    //   sit->second * n == the expected value for this genotype
    //
    double chisq  = 0.0;
    double exp    = 0.0;
    double obs    = 0.0;
    double df     = seg_ratios[marker].size() - 1;

    map<string, double>::iterator sit;

    for (sit = seg_ratios[marker].begin(); sit != seg_ratios[marker].end(); sit++) {
        obs = cnts.count(sit->first) == 0 ? 0 : cnts[sit->first];
        exp = sit->second * n;
        // cerr << "      category: " << sit->first << "; obs: " << obs << "; exp: " << exp << "\n";

        chisq += ((obs - exp) * (obs - exp)) / exp;
    }
    // cerr << "    df: " << df << "; Chisq value: " << chisq << "; pvalue: " << chisq_pvalue(df, chisq) << "\n";

    //
    // Determine p-value
    //
    return chisq_pvalue(df, chisq);
}

double
chisq_pvalue(int df, double chisq)
{
    int i = 0;
    while (chisq > chisq_crit_values[df][i] &&
           i < chisq_crit_values_size) {
        i++;
    }

    if (i == chisq_crit_values_size)
        return chisq_crit_values[0][chisq_crit_values_size - 1];

    return chisq_crit_values[0][i];
}

int
map_specific_genotypes(map<int, CSLocus *> &catalog, PopMap<CSLocus> *pmap, set<int> &parent_ids)
{
    map<string, string>               types;
    map<string, map<string, string> > dictionary;

    switch(map_type) {
    case dh:
        load_dh_dictionary(types, dictionary);
        break;
    case cp:
        load_cp_dictionary(types, dictionary);
        break;
    case bc1:
        load_bc_dictionary(types, dictionary);
        break;
    case f2:
        load_f2_dictionary(types, dictionary);
        break;
    case gen:
    case none:
    case unk:
        break;
    }

    map<int, CSLocus *>::iterator it;
    string   marker, m;
    Datum  **d;
    CSLocus *loc;

    for (it = catalog.begin(); it != catalog.end(); it++) {
        loc = it->second;
        loc->gcnt = 0;

        if (loc->marker.length() == 0) continue;

        if (types.count(loc->marker)) {
            loc->uncor_marker = loc->marker;
            loc->marker       = types[loc->marker];
            marker = loc->marker;
        } else {
            marker = "";
        }

        d = pmap->locus(loc->id);

        for (int i = 0; i < pmap->sample_cnt(); i++) {
            if (d[i] == NULL) continue;

            if (parent_ids.count(pmap->rev_sample_index(i))) continue;

            if (marker.length() == 0) {
                m = "--";
            } else {
                m = dictionary[marker].count(d[i]->gtype) ?
                    dictionary[marker][d[i]->gtype] :
                    dictionary[marker]["--"];
            }

            strcpy(d[i]->gtype, m.c_str());

            if (m != dictionary[marker]["--"])
                loc->gcnt++;
        }
    }

    return 0;
}

int
translate_genotypes(map<string, string> &types, map<string, map<string, string> > &dictionary,
                    map<int, CSLocus *> &catalog, PopMap<CSLocus> *pmap, map<int, string> &samples,
                    set<int> &parent_ids)
{
    map<int, CSLocus *>::iterator it;
    CSLocus *loc;

    for (it = catalog.begin(); it != catalog.end(); it++) {
        loc = it->second;

        string marker = types.count(loc->marker) ? types[loc->marker] : "";
        Datum **d = pmap->locus(loc->id);

        for (int i = 0; i < pmap->sample_cnt(); i++) {
            if (d[i] == NULL) continue;

            if (parent_ids.count(pmap->rev_sample_index(i))) continue;

            //cerr << "Examining progeny " << samples[pmap->rev_sample_index(i)] << "; marker: " << loc->marker << "\n";

            string m;

            if (marker.length() == 0) {
                m = dictionary[marker]["--"];
            } else {
                m = dictionary[marker].count(d[i]->gtype) ?
                    dictionary[marker][d[i]->gtype] :
                    dictionary[marker]["--"];
            }
            d[i]->trans_gtype = new char[m.length() + 1];

            //
            // If the genotype was corrected, output it in uppercase letters.
            //
            if (d[i]->corrected) {
                for (uint k = 0; k < m.length(); k++)
                    d[i]->trans_gtype[k] = toupper(m[k]);
                d[i]->trans_gtype[m.length()] = '\0';
            } else {
                strcpy(d[i]->trans_gtype, m.c_str());
            }
            if (m != dictionary[marker]["--"])
                loc->trans_gcnt++;
            //cerr << "  allele: " << d[i]->trans_gtype << "; trans_gcnt: " << loc->trans_gcnt << "\n";
        }
    }

    return 0;
}

int tally_progeny_haplotypes(CSLocus *locus, PopMap<CSLocus> *pmap, set<int> &parent_ids,
                             int &total, double &max, string &freq_str) {
    char gtype[id_len];
    map<string, double> freq;
    Datum **d = pmap->locus(locus->id);

    total = 0;
    max   = 0;

    //cerr << "Examining marker: " << locus->id << "\n";

    for (int i = 0; i < pmap->sample_cnt(); i++) {
        if (parent_ids.count(pmap->rev_sample_index(i))) continue;
        if (d[i] == NULL) continue;

        //cerr << "  Sample: " << i << "; Haplotype: " << d[i]->obshap[0] << "; Genotype: " << d[i]->gtype << "\n";
        if (strcmp(d[i]->gtype, "--") != 0) {
            //
            // Automated corrections will uppercase genotypes, convert them back to lowercase
            // in order to tally them properly.
            //
            int j = 0;
            while (d[i]->gtype[j] != '\0') {
                    gtype[j] = tolower(d[i]->gtype[j]);
                    j++;
            }
            gtype[j] = '\0';
            freq[gtype]++;
            total++;
        }
    }

    if (total == 0)
        return 0;

    double frac;
    stringstream s;
    char   f[id_len];
    map<string, double>::iterator it;
    for (it = freq.begin(); it != freq.end(); it++) {
        frac = (double) it->second / (double) total * 100;
        if (frac > max) max = frac;
        sprintf(f, "(%0.1f%%);", frac);
        s << it->first << ":" << it->second << f;
    }

    freq_str = s.str().substr(0, s.str().length() - 1);

    return 0;
}

int write_sql(map<int, CSLocus *> &catalog, PopMap<CSLocus> *pmap, set<int> &parent_ids) {

    if (map_type == none)
        return 0;

    stringstream pop_name;
    pop_name << "batch_" << batch_id << ".markers.tsv";
    string file = in_path + pop_name.str();

    cerr << "Writing SQL markers file to '" << file << "'\n";

    ofstream fh(file.c_str(), ofstream::out);

    if (fh.fail()) {
        cerr << "Error opening markers SQL file '" << file << "'\n";
        exit(1);
    }

    fh << "# SQL ID"         << "\t"
       << "Batch ID"         << "\t"
       << "Catalog Locus ID" << "\t"
       << "Marker Type"      << "\t"
       << "Total Genotypes"  << "\t"
       << "Max"              << "\t"
       << "Genotype Freqs"   << "\t"
       << "Segregation Distortion" << "\t"
       << "Mean Log Likelihood"    << "\t"
       << "Genotype Map"           << "\t"
       << "Uncorrected Marker"     << "\n";

    map<int, CSLocus *>::iterator it;
    map<string, string>::iterator j;
    CSLocus     *loc;
    char         max_str[id_len];
    stringstream gtype_map;

    for (it = catalog.begin(); it != catalog.end(); it++) {
        loc = it->second;

        if (loc->marker.length() == 0) continue;

        double max   = 0.0;
        int    total = 0;
        string freq, map;
        tally_progeny_haplotypes(loc, pmap, parent_ids, total, max, freq);

        sprintf(max_str, "%0.2f", max);

        //
        // Record the haplotype to genotype map.
        //
        gtype_map.str("");
        for (j = loc->gmap.begin(); j != loc->gmap.end(); j++)
            gtype_map << j->first << ":" << j->second << ";";
        map = gtype_map.str().substr(0, gtype_map.str().length() - 1);

        fh << 0           << "\t"
           << batch_id    << "\t"
           << loc->id     << "\t"
           << loc->marker << "\t"
           << total       << "\t"
           << max_str     << "\t"
           << freq        << "\t"
           << loc->chisq  << "\t"
           << loc->lnl    << "\t"
           << map         << "\t"
           << (loc->uncor_marker.length() == 0 ? loc->marker : loc->uncor_marker) << "\n";
    }

    fh.close();

    pop_name.str("");
    pop_name << "batch_" << batch_id << ".genotypes_" << progeny_limit << ".txt";
    file = in_path + pop_name.str();

    cerr << "Writing SQL genotypes file to '" << file << "'\n";

    fh.open(file.c_str(), ofstream::out);

    if (fh.fail()) {
        cerr << "Error opening genotypes SQL file '" << file << "'\n";
        exit(1);
    }

    fh << "# SQL ID"         << "\t"
       << "Batch ID"         << "\t"
       << "Catalog Locus ID" << "\t"
       << "Sample ID"        << "\t"
       << "Genotype"         << "\n";

    for (it = catalog.begin(); it != catalog.end(); it++) {
        loc = it->second;

        if (loc->gcnt < progeny_limit)
            continue;

        Datum **d = pmap->locus(loc->id);

        for (int i = 0; i < pmap->sample_cnt(); i++) {
            if (parent_ids.count(pmap->rev_sample_index(i))) continue;

            fh << 0 << "\t"
               << batch_id << "\t"
               << loc->id << "\t"
               << pmap->rev_sample_index(i) << "\t";

            if (d[i] == NULL)
                map_type == cp ? fh << "--\n" : fh << "-\n";
            else
                fh << d[i]->gtype << "\n";
        }
    }

    fh.close();

    return 0;
}

int write_genomic(map<int, CSLocus *> &catalog, PopMap<CSLocus> *pmap) {
    stringstream pop_name;
    pop_name << "batch_" << batch_id << ".genomic_" << progeny_limit << ".tsv";

    string file = in_path + pop_name.str();

    ofstream fh(file.c_str(), ofstream::out);

    if (fh.fail()) {
        cerr << "Error opening genomic output file '" << file << "'\n";
        exit(1);
    }

    //
    // Count the number of markers that have enough samples to output.
    //
    map<int, CSLocus *>::iterator cit;
    const CSLocus *loc;
    int num_loci = 0;

    for (cit = catalog.begin(); cit != catalog.end(); cit++) {
        loc = cit->second;
        if (loc->hcnt < progeny_limit) continue;

        num_loci += loc->len - renz_len[enz];
    }
    cerr << "Writing " << num_loci << " nucleotide positions to genomic file, '" << file << "'\n";

    //
    // Write the header
    //
    fh << num_loci << "\t" << pmap->sample_cnt() << "\n";

    //
    // Output each locus.
    //
    map<string, vector<CSLocus *> >::const_iterator it;
    int  a, b;

    uint  rcnt = renz_cnt[enz];
    uint  rlen = renz_len[enz];
    char *p;

    for (it = pmap->ordered_loci().begin(); it != pmap->ordered_loci().end(); it++) {
        for (uint i = 0; i < it->second.size(); i++) {
            loc = it->second[i];

            if (loc->hcnt < progeny_limit) continue;

            Datum **d = pmap->locus(loc->id);
            set<int> snp_locs;
            string   obshap;

            for (uint i = 0; i < loc->snps.size(); i++)
                snp_locs.insert(loc->snps[i]->col);

            uint start = 0;
            uint end   = loc->len;
            //
            // Check for the existence of the restriction enzyme cut site, mask off
            // its output.
            //
            for (uint n = 0; n < rcnt; n++)
                    if (strncmp(loc->con, renz[enz][n], rlen) == 0)
                        start += renz_len[enz];
            if (start == 0) {
                    p = loc->con + (loc->len - rlen);
                    for (uint n = rcnt; n < rcnt + rcnt; n++)
                        if (strncmp(p, renz[enz][n], rlen) == 0)
                            end -= renz_len[enz];
            }

            uint k = 0;
            for (uint n = start; n < end; n++) {
                fh << loc->id << "\t" << loc->loc.chr() << "\t" << loc->loc.bp + n;

                if (snp_locs.count(n) == 0) {
                    for (int j = 0; j < pmap->sample_cnt(); j++) {
                        a = encode_gtype(loc->con[n]);
                        fh << "\t" << encoded_gtypes[a][a];
                    }
                } else {
                    for (int j = 0; j < pmap->sample_cnt(); j++) {
                        fh << "\t";

                        if (d[j] == NULL)
                            fh << "0";
                        else
                            switch (d[j]->obshap.size()) {
                            case 1:
                                a = encode_gtype(d[j]->obshap[0][k]);
                                fh << encoded_gtypes[a][a];
                                break;
                            case 2:
                                a = encode_gtype(d[j]->obshap[0][k]);
                                b = encode_gtype(d[j]->obshap[1][k]);
                                fh << encoded_gtypes[a][b];
                                break;
                            default:
                                fh << "0";
                                break;
                            }
                    }
                    k++;
                }
                fh << "\n";
            }
        }
    }

    fh.close();

    return 0;
}

int write_generic(map<int, CSLocus *> &catalog, PopMap<CSLocus> *pmap, map<int, string> &samples, set<int> &parent_ids, bool write_gtypes) {

    stringstream pop_name;
    pop_name << "batch_" << batch_id;
    if (write_gtypes)
        pop_name << ".genotypes_" << progeny_limit << ".tsv";
    else
        pop_name << ".haplotypes_" << progeny_limit << ".tsv";

    string file = in_path + pop_name.str();

    ofstream fh(file.c_str(), ofstream::out);

    if (fh.fail()) {
        cerr << "Error opening generic output file '" << file << "'\n";
        exit(1);
    }

    //
    // Count the number of markers that have enough samples to output.
    //
    map<int, CSLocus *>::iterator it;
    CSLocus *loc;
    int num_loci = 0;

    for (it = catalog.begin(); it != catalog.end(); it++) {
        loc = it->second;
        if (write_gtypes == false && loc->hcnt < progeny_limit) continue;
        if (write_gtypes == true  && loc->gcnt < progeny_limit) continue;

        num_loci++;
    }
    cerr << "Writing " << num_loci << " loci to " << (write_gtypes ? "genotype" : "observed haplotype") << " file, '" << file << "'\n";

    //
    // Write the header
    //
    fh << "# Catalog ID\t";
    if (expand_id)
        fh << "\t";
    if (write_gtypes)
        fh << "Marker\t";
    fh << "Cnt\t"
       << "Seg Dist\t";

    map<int, string>::iterator s;
    for (int i = 0; i < pmap->sample_cnt(); i++) {
        if (write_gtypes && parent_ids.count(pmap->rev_sample_index(i)))
            continue;
        fh << samples[pmap->rev_sample_index(i)];
        if (i < pmap->sample_cnt() - 1)
            fh << "\t";
    }
    fh << "\n";

    //
    // Output each locus.
    //
    for (it = catalog.begin(); it != catalog.end(); it++) {
        loc = it->second;

        if (write_gtypes == false && loc->hcnt < progeny_limit) continue;
        if (write_gtypes == true  && loc->gcnt < progeny_limit) continue;

        stringstream id;
        loc->annotation.length() > 0 ?
            id << loc->id << "|" << loc->annotation : id << loc->id;

        fh << id.str();

        if (expand_id) {
            if (loc->annotation.length() > 0)
                id << "\t" << loc->id << "\t" << loc->annotation;
            else if (strlen(loc->loc.chr()) > 0)
                id << "\t" << loc->id << "\t" << loc->loc.chr() << "_" << loc->loc.bp;
            else
                id << "\t" << loc->id << "\t";
        }

        if (write_gtypes)
            fh << "\t" << loc->marker;

        write_gtypes ? fh << "\t" << loc->gcnt : fh << "\t" << loc->hcnt;
        fh << "\t" << loc->chisq;

        Datum **d = pmap->locus(loc->id);
        string  obshap;

        for (int i = 0; i < pmap->sample_cnt(); i++) {
            if (write_gtypes && parent_ids.count(pmap->rev_sample_index(i)))
                continue;
            fh << "\t";

            if (d[i] == NULL)
                fh << "-";
            else
                if (write_gtypes) {
                    fh << d[i]->gtype;
                } else {
                    obshap = "";
                    for (uint j = 0; j < d[i]->obshap.size(); j++)
                        obshap += string(d[i]->obshap[j]) + "/";
                    obshap = obshap.substr(0, obshap.length()-1);
                    fh << obshap;
                }
        }

        fh << "\n";
    }

    fh.close();

    return 0;
}

int
write_joinmap(map<int, CSLocus *> &catalog, PopMap<CSLocus> *pmap, map<string, string> &types, map<int, string> &samples, set<int> &parent_ids)
{
    stringstream pop_name;
    pop_name << "batch_" << batch_id << ".genotypes_" << progeny_limit;
    string file = in_path + pop_name.str() + ".loc";

    ofstream fh(file.c_str(), ofstream::out);

    if (fh.fail()) {
        cerr << "Error opening joinmap output file '" << file << "'\n";
        exit(1);
    }

    //
    // Count the number of mappable progeny
    //
    map<int, CSLocus *>::iterator it;
    CSLocus *loc;
    int num_loci = 0;

    for (it = catalog.begin(); it != catalog.end(); it++) {
        loc = it->second;
        if (loc->trans_gcnt < progeny_limit) continue;

        num_loci++;
    }
    cerr << "Writing " << num_loci << " loci to JoinMap file, '" << file << "'\n";

    map<int, string> map_types;
    map_types[cp]  = "CP";
    map_types[dh]  = "DH";
    map_types[bc1] = "BC1";
    map_types[f2]  = "F2";

    //
    // Output the header of the file
    //
    fh << "name = " << pop_name.str() << "\n"
       << "popt = " << map_types[map_type] << "\n"
       << "nloc = " << num_loci << "\n"
       << "nind = " << pmap->sample_cnt() - parent_ids.size() << "\n\n";

    //
    // Output each locus.
    //
    for (it = catalog.begin(); it != catalog.end(); it++) {
        loc = it->second;

        if (loc->trans_gcnt < progeny_limit) continue;

        stringstream id;
        loc->annotation.length() > 0 ?
            id << loc->id << "|" << loc->annotation : id << loc->id;

        fh << id.str() << "\t";

        if (expand_id) {
            id.str("");

            if (loc->annotation.length() > 0)
                id << loc->id << "\t" << loc->annotation;
            else if (strlen(loc->loc.chr()) > 0)
                id << loc->id << "\t" << loc->loc.chr() << "_" << loc->loc.bp;
            else
                id << loc->id << "\t";

            fh << id.str() << "\t";
        }

        if (types[loc->marker] == "lmx--")
            fh << "<lmxll>";
        else if (types[loc->marker] == "--xnp")
            fh << "<nnxnp>";
        else
            fh << "<" << types[loc->marker] << ">";

        Datum **d = pmap->locus(loc->id);

        for (int i = 0; i < pmap->sample_cnt(); i++) {
            if (parent_ids.count(pmap->rev_sample_index(i))) continue;
            fh << "\t";

            if (d[i] == NULL)
                map_type == cp ? fh << "--" : fh << "-";
            else
                fh << d[i]->trans_gtype;
        }

        fh << "\n";
    }

    fh << "\nindividual names:\n";

    map<int, string>::iterator s;
    for (s = samples.begin(); s != samples.end(); s++) {
        if (parent_ids.count(s->first)) continue;
        fh << s->second << "\n";
    }

    fh.close();

    return 0;
}

int
write_onemap_mapmaker(map<int, CSLocus *> &catalog, PopMap<CSLocus> *pmap, map<string, string> &types, map<int, string> &samples, set<int> &parent_ids)
{
    stringstream pop_name;
    pop_name << "batch_" << batch_id << ".genotypes_" << progeny_limit;
    string file = in_path + pop_name.str() + ".onemap.txt";

    ofstream fh(file.c_str(), ofstream::out);

    if (fh.fail()) {
        cerr << "Error opening joinmap output file '" << file << "'\n";
        exit(1);
    }

    //
    // Count the number of mappable progeny
    //
    map<int, CSLocus *>::iterator it;
    CSLocus *loc;
    int num_loci = 0;

    for (it = catalog.begin(); it != catalog.end(); it++) {
        loc = it->second;
        if (loc->trans_gcnt < progeny_limit) continue;

        num_loci++;
    }
    cerr << "Writing " << num_loci << " loci to OneMap file, '" << file << "'\n";

    //
    // Output map type.
    //
    if (map_type == f2 )
        fh << "data type f2 intercross\n";
    else if (map_type == bc1)
        fh << "data type f2 backcross\n";

    //
    // Output the header: number of individuals, number of markers, number of
    // quantitative traits (none).
    //
    fh << pmap->sample_cnt() - parent_ids.size() << " " << num_loci << " " << "0\n\n";

    //
    // Output each locus.
    //
    for (it = catalog.begin(); it != catalog.end(); it++) {
        loc = it->second;

        if (loc->trans_gcnt < progeny_limit) continue;

        fh << "*" << loc->id;

        Datum **d = pmap->locus(loc->id);

        for (int i = 0; i < pmap->sample_cnt(); i++) {
            if (parent_ids.count(pmap->rev_sample_index(i))) continue;
            fh << " ";

            if (d[i] == NULL)
                fh << "-";
            else
                fh << d[i]->trans_gtype;
        }
        fh << "\n";
    }

    fh.close();

    return 0;
}

int
write_onemap(map<int, CSLocus *> &catalog, PopMap<CSLocus> *pmap, map<string, string> &types, map<int, string> &samples, set<int> &parent_ids)
{
    stringstream pop_name;
    pop_name << "batch_" << batch_id << ".genotypes_" << progeny_limit;
    string file = in_path + pop_name.str() + "onemap.tsv";

    ofstream fh(file.c_str(), ofstream::out);

    if (fh.fail()) {
        cerr << "Error opening joinmap output file '" << file << "'\n";
        exit(1);
    }

    //
    // Count the number of mappable progeny
    //
    map<int, CSLocus *>::iterator it;
    CSLocus *loc;
    int num_loci = 0;

    for (it = catalog.begin(); it != catalog.end(); it++) {
        loc = it->second;
        if (loc->trans_gcnt < progeny_limit) continue;

        num_loci++;
    }
    cerr << "Writing " << num_loci << " loci to OneMap file, '" << file << "'\n";

    map<string, string> marker_types;
    marker_types["abxoo"] = "D1.11";
    marker_types["ooxab"] = "D2.16";
    marker_types["abxaa"] = "D1.10";
    marker_types["aaxab"] = "D2.15";
    marker_types["abxab"] = "B3.7";
    marker_types["abxac"] = "A.2";
    marker_types["abxcd"] = "A.1";


    //
    // Output the header: number of individuals followed by number of markers.
    //
    fh << pmap->sample_cnt() - parent_ids.size() << "\t" << num_loci << "\n";

    //
    // Output each locus.
    //
    for (it = catalog.begin(); it != catalog.end(); it++) {
        loc = it->second;

        if (loc->trans_gcnt < progeny_limit) continue;

        fh << "*" << loc->id << " "
           << marker_types[types[loc->marker]] << "\t";

        Datum **d = pmap->locus(loc->id);

        for (int i = 0; i < pmap->sample_cnt(); i++) {
            if (parent_ids.count(pmap->rev_sample_index(i))) continue;

            if (d[i] == NULL)
                fh << "-";
            else
                fh << d[i]->trans_gtype;

            if (i < pmap->sample_cnt() - 1)
                fh << ",";
        }

        fh << "\n";
    }

    fh.close();

    return 0;
}

int
write_rqtl(map<int, CSLocus *> &catalog, PopMap<CSLocus> *pmap, map<string, string> &types, map<int, string> &samples, set<int> &parent_ids)
{
    stringstream pop_name;
    pop_name << "batch_" << batch_id << ".genotypes_" << progeny_limit;
    string file = in_path + pop_name.str() + ".rqtl.tsv";

    ofstream fh(file.c_str(), ofstream::out);

    if (fh.fail()) {
        cerr << "Error opening R/QTL output file '" << file << "'\n";
        exit(1);
    }

    //
    // Count the number of mappable progeny
    //
    map<int, CSLocus *>::iterator it;
    CSLocus *loc;
    int num_loci = 0;

    for (it = catalog.begin(); it != catalog.end(); it++) {
        loc = it->second;
        if (loc->trans_gcnt < progeny_limit) continue;

        num_loci++;
    }
    cerr << "Writing " << num_loci << " loci to R/QTL file, '" << file << "'\n";

    map<int, string> map_types;
    map_types[cp]  = "CP";
    map_types[dh]  = "DH";
    map_types[bc1] = "BC1";
    map_types[f2]  = "F2";

     //
     // Output the header of the file, followed by the list of markers, one per column
     //
    fh << "# Exported: "   << pop_name.str() << "\n"
       << "# Map Type: "   << map_types[map_type] << "\n"
       << "# Num Loci: "   << num_loci << "\n"
       << "# Num Samples " << pmap->sample_cnt() - parent_ids.size() << "\n";

    for (it = catalog.begin(); it != catalog.end(); it++) {
        loc = it->second;
        if (loc->gcnt < progeny_limit) continue;

        fh << ",";

        stringstream id;
        loc->annotation.length() > 0 ?
            id << loc->id << "|" << loc->annotation : id << loc->id;
        fh << id.str();
    }
    fh << "\n";

    //
    // Output the chromosome (if available) for each marker and then the location
    //
    for (it = catalog.begin(); it != catalog.end(); it++) {
        loc = it->second;
        if (loc->gcnt < progeny_limit) continue;

        fh << ",";

        string chr;
        chr = strlen(loc->loc.chr()) > 0 ? loc->loc.chr() : "1";

        fh << chr;
    }
    fh << "\n";

    int i = 1;
    for (it = catalog.begin(); it != catalog.end(); it++) {
        loc = it->second;
        if (loc->gcnt < progeny_limit) continue;

        fh << ",";

        int bp = loc->loc.bp > 0 ? loc->loc.bp : i;

        fh << bp;
        i++;
    }
    fh << "\n";

    //
    // For each sample, print out the genotypes for each marker
    //
    Datum *d;
    for (int i = 0; i < pmap->sample_cnt(); i++) {
        if (parent_ids.count(pmap->rev_sample_index(i))) continue;

        fh << samples[pmap->rev_sample_index(i)];

        for (it = catalog.begin(); it != catalog.end(); it++) {
            loc = it->second;
            //if (loc->gcnt < progeny_limit) continue;

            d = pmap->datum(loc->id, pmap->rev_sample_index(i));
            fh << ",";

            if (d == NULL)
                map_type == cp ? fh << "--" : fh << "-";
            else
                fh << d->trans_gtype;
        }
        fh << "\n";
    }

    fh.close();

    return 0;
}

// sub create_imputed_genotype_map {
//     my ($order, $marker, $tag_id, $parents, $progeny, $map) = @_;

//     my (@keys, $key, $m, $type, $allele, $uall);

//     my (%gtypes, %allcnt, %uniqall);

//     //
//     // Count up the number of each type of observed allele in the progeny,
//     // record whether those genotypes are heterozygous or homozygous.
//     //
//     foreach $key (keys %{$progeny}) {
//         my $alleles;

//         print STDERR "Examining progeny $key\n" if ($debug);
//         //
//         // Discard progeny with more than one locus matched to this catalog tag.
//         //
//         @keys = keys %{$progeny->{$key}};
//         next if (scalar(@keys) > 1);

//         $alleles = join("|", sort @{$progeny->{$key}->{$keys[0]}});

//         if (!defined($allcnt{$alleles})) {
//             $allcnt{$alleles} = scalar(@{$progeny->{$key}->{$keys[0]}});
//         }
//         //print STDERR "Adding genotype $alleles\n";

//         $gtypes{$alleles}++;

//         foreach $allele (@{$progeny->{$key}->{$keys[0]}}) {
//             $uniqall{$allele}++;
//         }
//     }

//     //
//     // Examine the first parent alleles (the only alleles we have, since
//     // we are imputing the second parent.
//     //
//     my @parents = keys %{$parents};
//     my %legal_genotypes = ();
//     $key = $order->{$parents[0]} eq "first"  ? $parents[0] : $parents[1];
//     $m   = substr($marker, 0, 2);

//     foreach $type (split(//, $m)) {
//         //print STDERR "  Adding $type to genotypes\n" if ($debug);
//         $legal_genotypes{$type}++;
//     }
//     my @types = sort keys %legal_genotypes;

//     if ($marker eq "lmxll") {
//         @keys = sort {$gtypes{$b} <=> $gtypes{$a}} keys %gtypes;
//         //
//         // Discard heterozygous alleles and find the highest frequency homozygote,
//         // this is the "l" in the "lmxll" marker.
//         //
//         while ($allcnt{$keys[0]} == 2) {
//             shift @keys;
//         }
//         $map->{$keys[0]} = shift @types;
//         print STDERR "  Assinging '$keys[0]' to first parent genotype '", $map->{$keys[0]}, "'\n" if ($debug);

//         foreach $uall (sort {$uniqall{$b} <=> $uniqall{$a}} keys %uniqall) {
//             if ($uall ne $keys[0]) {
//                 $allele = $uall;
//                 last;
//             }
//         }
//         $map->{$allele} = shift @types;
//         print STDERR "  Assinging '$allele' to first parent genotype '", $map->{$allele}, "'\n" if ($debug);
//     }

// }

int load_marker_list(string path, set<int> &list) {
    char     line[id_len];
    ifstream fh(path.c_str(), ifstream::in);

    if (fh.fail()) {
        cerr << "Error opening white/black list file '" << path << "'\n";
        exit(1);
    }

    int   marker;
    char *e;

    while (fh.good()) {
        fh.getline(line, id_len);

        if (strlen(line) == 0) continue;

        marker = (int) strtol(line, &e, 10);

        if (*e == '\0')
            list.insert(marker);
    }

    fh.close();

    if (list.size() == 0) {
        cerr << "Unable to load any markers from '" << path << "'\n";
        help();
    }

    return 0;
}

bool hap_compare(pair<string, int> a, pair<string, int> b) {
    return (a.second > b.second);
}

int parse_command_line(int argc, char* argv[]) {
    int c;

    while (1) {
        static struct option long_options[] = {
            {"help",          no_argument,       NULL, 'h'},
            {"version",       no_argument,       NULL, 'v'},
            {"corr",          no_argument,       NULL, 'c'},
            {"sql",           no_argument,       NULL, 's'},
            {"num_threads",   required_argument, NULL, 'p'},
            {"batch_id",      required_argument, NULL, 'b'},
            {"in_path",       required_argument, NULL, 'P'},
            {"map_type",      required_argument, NULL, 't'},
            {"out_type",      required_argument, NULL, 'o'},
            {"progeny",       required_argument, NULL, 'r'},
            {"min_depth",     required_argument, NULL, 'm'},
            {"min_hom_seqs",  required_argument, NULL, 'H'},
            {"min_het_seqs",  required_argument, NULL, 'N'},
            {"max_het_seqs",  required_argument, NULL, 'X'},
            {"renz",          required_argument, NULL, 'e'},
            {"whitelist",     required_argument, NULL, 'W'},
            {"blacklist",     required_argument, NULL, 'B'},
            {"man_corr",      required_argument, NULL, 'C'},
            {"lnl_lim",       required_argument, NULL, 'L'},
            {0, 0, 0, 0}
        };

        // getopt_long stores the option index here.
        int option_index = 0;

        c = getopt_long(argc, argv, "hvcsib:p:t:o:r:P:m:e:H:N:X:W:B:C:L:", long_options, &option_index);

        // Detect the end of the options.
        if (c == -1)
            break;

        switch (c) {
        case 'h':
            help();
            break;
        case 'P':
            in_path = optarg;
            break;
        case 'b':
            batch_id = is_integer(optarg);
            if (batch_id < 0) {
                cerr << "Batch ID (-b) must be an integer, e.g. 1, 2, 3\n";
                help();
            }
            break;
        case 't':
            if (strcasecmp(optarg, "cp") == 0)
                map_type = cp;
            else if (strcasecmp(optarg, "bc1") == 0)
                map_type = bc1;
            else if (strcasecmp(optarg, "f2") == 0)
                map_type = f2;
            else if (strcasecmp(optarg, "dh") == 0)
                map_type = dh;
            else if (strcasecmp(optarg, "gen") == 0)
                map_type = gen;
            else
                map_type = unk;
            break;
        case 'o':
            if (strcasecmp(optarg, "joinmap") == 0)
                out_type = joinmap;
            else if (strcasecmp(optarg, "rqtl") == 0)
                out_type = rqtl;
            else if (strcasecmp(optarg, "onemap") == 0)
                out_type = onemap;
            else if (strcasecmp(optarg, "genomic") == 0)
                out_type = genomic;
            break;
        case 'r':
            progeny_limit = atoi(optarg);
            break;
        case 'c':
            corrections = true;
            break;
        case 'L':
            lnl_limit  = is_double(optarg);
            filter_lnl = true;
            break;
        case 'i':
            expand_id = true;
            break;
        case 's':
            sql_out = true;
            break;
        case 'W':
            wl_file = optarg;
            break;
        case 'B':
            bl_file = optarg;
            break;
        case 'C':
            man_corrections = true;
            cor_path = optarg;
            break;
        case 'm':
            min_stack_depth = is_integer(optarg);
            break;
        case 'H':
            min_hom_seqs = is_integer(optarg);
            break;
        case 'N':
            min_het_seqs = is_double(optarg);
            break;
        case 'X':
            max_het_seqs = is_double(optarg);
            break;
        case 'e':
            enz = optarg;
            enz.at(0) = tolower(enz.at(0));
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

    if (in_path.at(in_path.length() - 1) != '/')
        in_path += "/";

    if (batch_id < 0) {
        cerr << "You must specify a batch ID.\n";
        help();
    }

    if (map_type != cp  &&
        map_type != dh  &&
        map_type != bc1 &&
        map_type != f2  &&
        map_type != gen &&
        map_type != none) {
        cerr << "You must specify a valid map type. 'CP', 'DH', 'F2', 'BC1' and 'GEN' are the currently supported map types.\n";
        help();
    }

    if (map_type != none && min_stack_depth > 0)
        cerr << "Warning: using a minimum stack depth when building genetic markers is not recommended.\n";

    if (out_type == genomic && enz.length() == 0) {
        cerr << "You must specify the restriction enzyme used with 'genomic' output.\n";
        help();
    }

    if (out_type == genomic && renz.count(enz) == 0) {
        cerr << "Unrecognized restriction enzyme specified: '" << enz.c_str() << "'.\n";
        help();
    }

    return 0;
}

void version() {
    cerr << "genotypes " << VERSION << "\n\n";

    exit(1);
}

void help() {
    cerr << "genotypes " << VERSION << "\n"
              << "genotypes -b batch_id -P path [-r min] [-m min] [-t map_type -o type] [-B blacklist] [-W whitelist] [-c] [-s] [-e renz] [-v] [-h]" << "\n"
              << "  b: Batch ID to examine when exporting from the catalog.\n"
              << "  r: minimum number of progeny required to print a marker.\n"
              << "  c: make automated corrections to the data.\n"
              << "  P: path to the Stacks output files.\n"
              << "  t: map type to write. 'CP', 'DH', 'F2', 'BC1' and 'GEN' are the currently supported map types.\n"
              << "  o: output file type to write, 'joinmap', 'onemap', 'rqtl', and 'genomic' are currently supported.\n"
              << "  m: specify a minimum stack depth required before exporting a locus in a particular individual.\n"
              << "  s: output a file to import results into an SQL database.\n"
              << "  B: specify a file containing Blacklisted markers to be excluded from the export.\n"
              << "  W: specify a file containign Whitelisted markers to include in the export.\n"
              << "  e: restriction enzyme, required if generating 'genomic' output.\n"
              << "  v: print program version." << "\n"
              << "  h: display this help messsage." << "\n"
              << "  Filtering options:\n"
              << "    --lnl_lim [num]: filter loci with log likelihood values below this threshold.\n"
              << "  Automated corrections options:\n"
              << "    --min_hom_seqs: minimum number of reads required at a stack to call a homozygous genotype (default 5).\n"
              << "    --min_het_seqs: below this minor allele frequency a stack is called a homozygote, above it (but below --max_het_seqs) it is called unknown (default 0.05).\n"
              << "    --max_het_seqs: minimum frequency of minor allele to call a heterozygote (default 0.1).\n"
              << "  Manual corrections options:\n"
              << "    --cor_path <path>: path to file containing manual genotype corrections from a Stacks SQL database to incorporate into output.\n";

    exit(1);
}
