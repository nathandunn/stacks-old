//
// Created by ndunn on 3/4/13.
//
// To change the template use AppCode | Preferences | File Templates.
//


#include "LociLoader.hpp"

bool LociLoader::hap_compare(pair<string, int> a, pair<string, int> b) {
    return (a.second > b.second);
}

int LociLoader::tabulate_haplotypes(map<int, CSLocus *>  & catalog, PopMap<CSLocus> *pmap) {
    map<int, CSLocus *>::iterator it;
    vector<char *>::iterator hit;
    Datum  **d;
    CSLocus *loc;

    for (it = catalog.begin(); it != catalog.end(); it++) {
        loc = it->second;
        d   = pmap->locus(loc->id);

        for (int i = 0; i < pmap->sample_cnt(); i++) {
            if (d[i] == NULL)
                continue;

            if (d[i]->obshap.size() > 1)
                loc->marker = "heterozygous";
        }

        if (loc->marker.length() > 0) {
            create_genotype_map(loc, pmap);
            call_population_genotypes(loc, pmap);
        }
    }

    return 0;
}

int LociLoader::create_genotype_map(CSLocus *locus, PopMap<CSLocus> *pmap) {
    //
    // Create a genotype map. For any set of haplotypes, this routine will
    // assign each haplotype to a genotype, e.g. given the haplotypes
    // 'AC' and 'GT' in the population, this routine will assign 'AC' == 'a'
    // and 'GT' == 'b'. If an individual is homozygous for 'AC', they will be
    // assigned an 'aa' genotype.
    //
    //cerr << "Creating genotype map for catalog ID " << locus->id  << ", marker: " << locus->marker << ".\n";

    char gtypes[26] ={'a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j',
            'k', 'l', 'm', 'n', 'o', 'p', 'q', 'r', 's', 't',
            'u', 'v', 'w', 'x', 'y', 'z'};

    Datum **d;
    map<string, int> haplotypes;
    map<string, int>::iterator k;
    vector<pair<string, int> > sorted_haplotypes;

    d = pmap->locus(locus->id);

    for (int i = 0; i < pmap->sample_cnt(); i++) {

        if (d[i] != NULL)
            for (uint n = 0; n < d[i]->obshap.size(); n++)
                haplotypes[d[i]->obshap[n]]++;
    }

    //
    // Check that there are not more haplotypes than we have encodings.
    //
    if (haplotypes.size() > 26) return 0;

    //
    // Sort the haplotypes map by value
    //
    for (k = haplotypes.begin(); k != haplotypes.end(); k++)
        sorted_haplotypes.push_back(*k);
    sort(sorted_haplotypes.begin(), sorted_haplotypes.end(), hap_compare);

    for (uint n = 0, index = 0; n < sorted_haplotypes.size() && index <= 26; n++, index++) {
        locus->gmap[sorted_haplotypes[n].first] = gtypes[index];
        //cerr << "GMAP: " << sorted_haplotypes[n].first << " == " << gtypes[index] << "\n";
    }

    return 0;

}

int LociLoader::call_population_genotypes(CSLocus *locus, PopMap<CSLocus> *pmap) {
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

        string m = gtype.length() == 1 ?
                gtype + gtype : gtype;

        d[i]->gtype = new char[m.length() + 1];
        strcpy(d[i]->gtype, m.c_str());

        if (m != "-")
            locus->gcnt++;

        //cerr << "Assigning datum, marker: " << locus->marker << ", string: " << m << ", haplotype: " << d[i]->obshap[0] << ", gtype: " << gtype << "\n";
    }

    return 0;

}
