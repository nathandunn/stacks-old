// -*-mode:c++; c-style:k&r; c-basic-offset:4;-*-
//
// Copyright 2013, Julian Catchen <jcatchen@uoregon.edu>
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
// catalog_utils.cc -- common routines for manipulating catalog objects.
//
// Julian Catchen
// jcatchen@uoregon.edu
// University of Oregon
//
#include <regex>

#include "utils.h"

#include "catalog_utils.h"

using namespace std;

const regex catalog_tags_regex ("^batch_([0-9]+).catalog.tags.tsv(.gz)?$");

vector<int> find_catalogs(const string& dir_path) {
    vector<int> ids;

    for (DirIterator e (dir_path); e; ++e) {
        smatch m;
        string name (e.name());
        regex_match(name, m, catalog_tags_regex);
        if (!m.empty())
            ids.push_back(stoi(m[1].str()));
    }

    return ids;
}

int
reduce_catalog(map<int, CSLocus *> &catalog, set<int> &whitelist, set<int> &blacklist)
{
    map<int, CSLocus *> list;
    map<int, CSLocus *>::iterator it;
    CSLocus *loc;

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

int
check_whitelist_integrity(map<int, CSLocus *> &catalog, map<int, set<int> > &whitelist)
{
    if (whitelist.size() == 0) return 0;

    int rm_snps = 0;
    int rm_loci = 0;

    CSLocus *loc;
    map<int, set<int> >::iterator it;
    set<int>::iterator sit;
    map<int, set<int> > new_wl;

    cerr << "Checking the integrity of the whitelist...";

    for (it = whitelist.begin(); it != whitelist.end(); it++) {
        if (catalog.count(it->first) == 0) {
            rm_loci++;
            cerr << "\n  Removing locus " << it->first << " from whitelist as it does not exist in the catalog.";
        } else {
            loc = catalog[it->first];

            if (it->second.size() == 0) {
                new_wl.insert(make_pair(it->first, set<int>()));
                continue;
            }

            set<int> cat_snps;
            for (uint i = 0; i < loc->snps.size(); i++)
                cat_snps.insert(loc->snps[i]->col);

            for (sit = it->second.begin(); sit != it->second.end(); sit++)
                if (cat_snps.count(*sit)) {
                    new_wl[it->first].insert(*sit);
                } else {
                    rm_snps++;
                    cerr << "\n  Removing SNP at column " << *sit << " in locus " << it->first << " from whitelist as it does not exist in the catalog.";
                }
        }
    }

    whitelist = new_wl;

    if (rm_loci > 0 || rm_snps > 0) cerr << "\n";

    cerr << "done.\n"
         << "Removed " << rm_loci << " loci and " << rm_snps << " SNPs from the whitelist that were not found in the catalog.\n";

    return 0;
}

int
reduce_catalog(map<int, CSLocus *> &catalog, map<int, set<int> > &whitelist, set<int> &blacklist)
{
    map<int, CSLocus *> list;
    map<int, CSLocus *>::iterator it;
    CSLocus *loc;

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

/*int
reduce_catalog_snps(map<int, CSLocus *> &catalog, map<int, set<int> > &whitelist, PopMap<CSLocus> *pmap)
{
    map<int, CSLocus *>::iterator it;
    CSLocus *loc;
    Datum  **d;

    if (whitelist.size() == 0)
        return 0;

    //
    // We want to prune out SNP objects that are not in the whitelist.
    //
    int              pos;
    vector<SNP *>    tmp;
    vector<uint>     cols;
    map<string, int> obshaps;
    map<string, int>::iterator sit;
    for (it = catalog.begin(); it != catalog.end(); it++) {
        loc = it->second;

        if (whitelist[loc->id].size() == 0)
            continue;

        tmp.clear();
        cols.clear();

        d = pmap->locus(loc->id);

        for (uint i = 0; i < loc->snps.size(); i++) {
            if (whitelist[loc->id].count(loc->snps[i]->col) > 0) {
                tmp.push_back(loc->snps[i]);
                cols.push_back(i);
            } else {
                //
                // Change the model calls in the samples to no longer contain this SNP.
                //
                pos = loc->snps[i]->col;
                for (int j = 0; j < pmap->sample_cnt(); j++) {
                    if (d[j] == NULL || pos >= d[j]->len)
                        continue;
                    if (d[j]->model != NULL) {
                        d[j]->model[pos] = 'U';
                    }
                }

                delete loc->snps[i];
            }
        }
        loc->snps.clear();
        for (uint i = 0; i < tmp.size(); i++)
            loc->snps.push_back(tmp[i]);

        map<string, int>::iterator it;
        char allele_old[id_len], allele_new[id_len];
        //
        // We need to adjust the catalog's list of haplotypes/alleles
        // for this locus to account for the pruned SNPs.
        //
        for (it = loc->alleles.begin(); it != loc->alleles.end(); it++) {
            strncpy(allele_old, it->first.c_str(), id_len - 1);
            allele_old[id_len - 1] = '\0';

            for (uint k = 0; k < cols.size(); k++)
                allele_new[k] = allele_old[cols[k]];
            allele_new[cols.size()] = '\0';
            obshaps[string(allele_new)] += it->second;
        }
        loc->alleles.clear();
        for (sit = obshaps.begin(); sit != obshaps.end(); sit++) {
            loc->alleles[sit->first] = sit->second;
        }
        obshaps.clear();

        loc->populate_alleles();

        //
        // Now we need to adjust the matched haplotypes to sync to
        // the SNPs left in the catalog.
        //
        // Reducing the lengths of the haplotypes  may create
        // redundant (shorter) haplotypes, we need to remove these.
        //
        for (int i = 0; i < pmap->sample_cnt(); i++) {
            if (d[i] == NULL) continue;

            for (uint j = 0; j < d[i]->obshap.size(); j++) {
                for (uint k = 0; k < cols.size(); k++)
                    d[i]->obshap[j][k] = d[i]->obshap[j][cols[k]];
                d[i]->obshap[j][cols.size()] = '\0';
                obshaps[d[i]->obshap[j]] += d[i]->depth[j];
            }
            uint j = 0;
            for (sit = obshaps.begin(); sit != obshaps.end(); sit++) {
                strcpy(d[i]->obshap[j], sit->first.c_str());
                d[i]->depth[j] = sit->second;
                j++;
            }
            while (j < d[i]->obshap.size()) {
                delete [] d[i]->obshap[j];
                j++;
            }
            d[i]->obshap.resize(obshaps.size());
            d[i]->depth.resize(obshaps.size());
            obshaps.clear();
        }
    }

    return 0;
}*/

CSLocus *
new_cslocus(const VcfRecord rec, int id)
{
    CSLocus* loc = new CSLocus();

    loc->sample_id = 0;
    loc->id        = id;
    loc->len       = 1;
    loc->con       = new char[2];
    strcpy(loc->con, rec.allele0());

    loc->loc.set(rec.chrom(), (uint)rec.pos(), strand_plus);

    for (auto a = rec.begin_alleles(); a != rec.end_alleles(); ++a) {
        if (strcmp(*a, "*")==0)
            continue;
        loc->alleles.insert({string(*a), 0});
    }

    loc->depth = 0;

    loc->snps.push_back(new SNP());
    SNP& snp = *loc->snps.back();
    snp.col = 0;
    snp.type = snp_type_unk;
    vector<char*> snp_alleles = {&snp.rank_1, &snp.rank_2, &snp.rank_3, &snp.rank_4};
    try {
        uint j = 0;
        for (auto a=rec.begin_alleles(); a!=rec.end_alleles(); ++a) {
            assert(strlen(*a) == 1);
            if (**a == '*')
                continue;
            *snp_alleles.at(j) = **a;
            j++;
        }
    } catch (out_of_range& e) {
        cerr << "Warning: Skipping malformed VCF SNP record '"
             << rec.chrom() << ":" << rec.pos() << "'."
             << " Alleles were:";
        for (auto a=rec.begin_alleles(); a!=rec.end_alleles(); ++a)
            cerr << " '" << *a << "';";
        cerr << ".\n";
        delete loc->snps[0];
        delete loc->con;
        delete loc;

        return NULL;
    }

    loc->populate_alleles();

    return loc;
}

CSLocus* new_cslocus(const Seq& consensus, const vector<VcfRecord>& records, int id) {
    CSLocus* loc = new CSLocus();

    // sample_id
    loc->sample_id = 0;

    // id
    loc->id = id;

    // con + len
    assert(consensus.seq != NULL);
    loc->add_consensus(consensus.seq);

    // loc
    // If the analysis is reference-based, there will be a ' pos=...' field on
    // the fasta ID line, placed as a comment adjacent to the locus ID..
    assert(loc->loc.empty());

    const char *p, *q;
    p = consensus.comment;

    do {
        q = p;
        while (*q != '\0' && *q != ' ' && *q != '\t')
            ++q;

        if (strncmp(p, "pos=", 4) == 0) {
            p += 4;
            loc->loc = PhyLoc(string(p, q));
            break;
        }

        p = *q == '\0' ? q : q + 1;

    } while (*p != '\0');

    if (loc->loc.empty())
        loc->loc = PhyLoc("", 0, strand_plus); // n.b. Not the same as PhyLoc(); with this `PhyLoc::chr != NULL`.

    // snps
    vector<const VcfRecord*> snp_records;
    for (auto& rec : records) {
        if (!rec.is_monomorphic()) {
            snp_records.push_back(&rec);
            assert(rec.is_snp());
        }
    }
    for (const VcfRecord* rec : snp_records) {
        SNP* snp = new SNP();
        loc->snps.push_back(snp);
        snp->type = snp_type_het;
        snp->col = rec->pos();
        auto a = rec->begin_alleles();
        snp->rank_1 = (*a)[0];
        ++a;
        assert(a != rec->end_alleles());
        snp->rank_2 = (*a)[0];
        ++a;
        if (a != rec->end_alleles()) {
            snp->rank_3 = (*a)[0];
            ++a;
            if (a != rec->end_alleles())
                snp->rank_4 = (*a)[0];
        }
        snp->lratio = 0.0;
    }

    // alleles
    if (!snp_records.empty()) {
        pair<string,string> haplotypes;
        for (size_t sample=0; sample<records.at(0).count_samples(); ++sample) {
            VcfRecord::util::build_haps(haplotypes, snp_records, sample);
            bool complete = true;
            for (size_t i=0; i<snp_records.size(); ++i) {
                if (haplotypes.first[i] == 'N' || haplotypes.second[i] == 'N') {
                    complete = false;
                    break;
                }
            }
            if (complete) {
                ++loc->alleles[haplotypes.first];
                ++loc->alleles[haplotypes.second];
            }
        }
    }

    // strings
    loc->populate_alleles();

    // depth
    loc->depth = 0;

    return loc;
}
