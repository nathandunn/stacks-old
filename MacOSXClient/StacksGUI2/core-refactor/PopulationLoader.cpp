//
// Created by ndunn on 3/1/13.
//
// To change the template use AppCode | Preferences | File Templates.
//


#include "PopulationLoader.hpp"
#include "CLocus.hpp"
#include "input.h"
#include "dirent.h"
#include "PopMap.h"

using std::pair ;

int PopulationLoader::reduce_catalog(map<int, CLocus *> catalog,set<int> whitelist,set<int> blacklist){
    map<int, CLocus *> list;
    map<int, CLocus *>::iterator it;
    CLocus *loc;

    if (whitelist.size() == 0 && blacklist.size() == 0)
        return 0;

    int i = 0;
    for (it = catalog.begin(); it != catalog.end(); it++) {
        loc = it->second;

        if (whitelist.size() > 0 && whitelist.count(loc->id) == 0) continue;
        if (blacklist.count(loc->id)) continue;

        list[it->first] = it->second;
        i++;
    }

    catalog = list;

    return i;
}

int PopulationLoader::reduce_catalog(map<int, CLocus *> catalog){

    set<int> whitelist, blacklist;
    return reduce_catalog(catalog,whitelist,blacklist);

}

bool PopulationLoader::order_unordered_loci(map<int, CLocus *> catalog) {
    map<int, CLocus *>::iterator it;
    CLocus *loc;
    set<string> chrs;

    for (it = catalog.begin(); it != catalog.end(); it++) {
        loc = it->second;
        if (strlen(loc->loc.chr) > 0)
            chrs.insert(loc->loc.chr);
    }

    //
    // This data is already reference aligned.
    //
    if (chrs.size() > 0)
        return true;

    cerr << "Catalog is not reference aligned, arbitrarily ordering catalog loci.\n";

    uint bp = 1;
    for (it = catalog.begin(); it != catalog.end(); it++) {
        loc = it->second;
        loc->loc.chr = new char[3];
        strcpy(loc->loc.chr, "un");
        loc->loc.bp  = bp;

        bp += strlen(loc->con);
    }

    return false;

}

int PopulationLoader::build_file_list(string pmap_path, vector<int, string> files , map<int, pair<int,int>> &pop_indexes) {
    char   line[max_len];
    char   pop_id_str[id_len];
    vector<string> parts;
    string f;
    uint   len;
    string in_path ;
    int population_limit = 1;


    if (pmap_path.length() > 0) {
        cerr << "Parsing population map.\n";

        ifstream fh(pmap_path.c_str(), ifstream::in);

        if (fh.fail()) {
            cerr << "Error opening population map '" << pmap_path << "'\n";
            return 0;
        }

        while (fh.good()) {
            fh.getline(line, max_len);

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

            //
            // Parse the population map, we expect:
            // <file name> <tab> <population ID>
            //
            parse_tsv(line, parts);

            if (parts.size() != 2) {
                cerr << "Population map is not formated correctly: expecting two, tab separated columns, found " << parts.size() << ".\n";
                return 0;
            }

            strncpy(pop_id_str, parts[1].c_str(), id_len);
            for (int i = 0; i < id_len && pop_id_str[i] != '\0'; i++)
                if (!isdigit(pop_id_str[i])) {
                    cerr << "Population map is not formated correctly: expecting numerical ID in second column, found '" << parts[1] << "'.\n";
                    return 0;
                }

            //
            // Test that file exists before adding to list.
            //
            f = in_path.c_str() + parts[0] + ".matches.tsv";
            ifstream test_fh(f.c_str(), ifstream::in);

            if (test_fh.fail()) {
                cerr << " Unable to find " << f.c_str() << ", excluding it from the analysis.\n";
            } else {
                test_fh.close();
                files.push_back(make_pair(atoi(parts[1].c_str()), parts[0]));
            }
        }

        fh.close();
    } else {
        cerr << "No population map specified, building file list.\n";

        //
        // If no population map is specified, read all the files from the Stacks directory.
        //
        uint   pos;
        string file;
        struct dirent *direntry;

        DIR *dir = opendir(in_path.c_str());

        if (dir == NULL) {
            cerr << "Unable to open directory '" << in_path << "' for reading.\n";
            exit(1);
        }

        while ((direntry = readdir(dir)) != NULL) {
            file = direntry->d_name;

            if (file == "." || file == "..")
                continue;

            if (file.substr(0, 6) == "batch_")
                continue;

            pos = file.rfind(".tags.tsv");
            if (pos < file.length())
                files.push_back(make_pair(1, file.substr(0, pos)));
        }

        closedir(dir);
    }

    if (files.size() == 0) {
        cerr << "Unable to locate any input files to process within '" << in_path << "'\n";
        return 0;
    }

    //
    // Sort the files according to population ID.
    //
    sort(files.begin(), files.end(), compare_pop_map);

    cerr << "Found " << files.size() << " input file(s).\n";

    //
    // Determine the start/end index for each population in the files array.
    //
    int start  = 0;
    int end    = 0;
    int pop_id = files[0].first;

    do {
        end++;
        if (pop_id != files[end].first) {
            pop_indexes[pop_id] = make_pair(start, end - 1);
            start  = end;
            pop_id = files[end].first;
        }
    } while (end < (int) files.size());

    cerr << "  " << pop_indexes.size() << " populations found\n";

    if (population_limit > (int) pop_indexes.size()) {
        cerr << "Population limit ("
                << population_limit
                << ") larger than number of popualtions present, adjusting parameter to "
                << pop_indexes.size() << "\n";
        population_limit = pop_indexes.size();
    }

    map<int, pair<int, int> >::iterator it;
    for (it = pop_indexes.begin(); it != pop_indexes.end(); it++) {
        start = it->second.first;
        end   = it->second.second;
        cerr << "    " << it->first << ": ";
        for (int i = start; i <= end; i++) {
            cerr << files[i].second;
            if (i < end) cerr << ", ";
        }
        cerr << "\n";
    }

    return 1;
}


bool PopulationLoader::compare_pop_map(pair<int, string> a, pair<int, string> b) {
    if (a.first == b.first)
        return (a.second < b.second);
    return (a.first < b.first);
}


int PopulationLoader::apply_locus_constraints(map<int, CLocus *> &catalog, PopMap<CLocus> *pmap, map<int, pair<int, int> > &pop_indexes){
    uint pop_id, start_index, end_index;
    CLocus *loc;
    Datum **d;

    double    sample_limit      = 0;
    int       progeny_limit     = 0;
    int       population_limit  = 1;
    int       min_stack_depth   = 0;


    if (sample_limit == 0 && population_limit == 0 && min_stack_depth == 0) return 0;

    map<int, CLocus *>::iterator it;
    map<int, pair<int, int> >::iterator pit;

    uint pop_cnt   = pop_indexes.size();
    int *pop_order = new int [pop_cnt];

    // Which population each sample belongs to.
    int *samples   = new int [pmap->sample_cnt()];

    // For the current locus, how many samples in each population.
    int *pop_cnts  = new int [pop_cnt];

    // The total number of samples in each population.
    int *pop_tot   = new int [pop_cnt];

    pop_id = 0;
    for (pit = pop_indexes.begin(); pit != pop_indexes.end(); pit++) {
        start_index = pit->second.first;
        end_index   = pit->second.second;
        pop_tot[pop_id]  = 0;

        for (uint i = start_index; i <= end_index; i++) {
            samples[i] = pop_id;
            pop_tot[pop_id]++;
        }
        pop_order[pop_id] = pit->first;
        pop_id++;
    }

    for (uint i = 0; i < pop_cnt; i++)
        pop_cnts[i] = 0;

    double pct       = 0.0;
    bool   pro_limit = false;
    bool   pop_limit = false;
    int    pops      = 0;
    int    below_stack_dep = 0;
    set<int> blacklist;

    for (it = catalog.begin(); it != catalog.end(); it++) {
        loc = it->second;
        d   = pmap->locus(loc->id);

        //
        // Check that each sample is over the minimum stack depth for this locus.
        //
        if (min_stack_depth > 0)
            for (int i = 0; i < pmap->sample_cnt(); i++) {
                if (d[i] != NULL && d[i]->tot_depth < min_stack_depth) {
                    below_stack_dep++;
                    delete d[i];
                    d[i] = NULL;
                    loc->hcnt--;
                }
            }

        //
        // Tally up the count of samples in this population.
        //
        for (int i = 0; i < pmap->sample_cnt(); i++) {
            if (d[i] != NULL)
                pop_cnts[samples[i]]++;
        }

        //
        // Check that the counts for each population are over progeny_limit. If not, zero out
        // the members of that population.
        //
        for (uint i = 0; i < pop_cnt; i++) {
            pct = (double) pop_cnts[i] / (double) pop_tot[i];

            if (pop_cnts[i] > 0 && pct < sample_limit) {
                //cerr << "Removing population " << pop_order[i] << " at locus: " << loc->id << "; below sample limit: " << pct << "\n";
                start_index = pop_indexes[pop_order[i]].first;
                end_index   = pop_indexes[pop_order[i]].second;

                for (uint j  = start_index; j <= end_index; j++) {
                    if (d[j] != NULL) {
                        delete d[j];
                        d[j] = NULL;
                        loc->hcnt--;
                    }
                }
                pop_cnts[i] = 0;
            }
        }

        //
        // Check that this locus is present in enough populations.
        //
        for (uint i = 0; i < pop_cnt; i++)
            if (pop_cnts[i] > 0) pops++;
        if (pops < population_limit) {
            //cerr << "Removing locus: " << loc->id << "; below population limit: " << pops << "\n";
            pop_limit = true;
        }

        if (pop_limit)
            blacklist.insert(loc->id);

        for (uint i = 0; i < pop_cnt; i++)
            pop_cnts[i] = 0;
        pro_limit = false;
        pop_limit = false;
        pops      = 0;
    }

    //
    // Remove loci
    //
    if (min_stack_depth > 0)
        cerr << "Removed " << below_stack_dep << " samples from loci that are below the minimum stack depth of " << min_stack_depth << "x\n";
    cerr << "Removing " << blacklist.size() << " loci that did not pass sample/population constraints...";
    set<int> whitelist;
    reduce_catalog(catalog, whitelist, blacklist);
    int retained = pmap->prune(blacklist);
    cerr << " retained " << retained << " loci.\n";

    delete [] pop_cnts;
    delete [] pop_tot;
    delete [] pop_order;
    delete [] samples;

    if (retained == 0)
        exit(0);

    return 0;

}

