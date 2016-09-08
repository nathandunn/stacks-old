#include "models.h"
#include "pstacks_base.h"

using std::set;
using std::map;
using std::vector;

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

int call_consensus(map<int, MergedStack *> &merged, map<int, PStack *> &unique, bool invoke_model) {
    //
    // OpenMP can't parallelize random access iterators, so we convert
    // our map to a vector of integer keys.
    //
    map<int, MergedStack *>::iterator it;
    vector<int> keys;
    for (it = merged.begin(); it != merged.end(); it++)
        keys.push_back(it->first);

    int i;
    #pragma omp parallel private(i)
    {
        #pragma omp for schedule(dynamic)
        for (i = 0; i < (int) keys.size(); i++) {
            MergedStack *mtag;
            PStack *utag;

            mtag = merged[keys[i]];

            //
            // Create a two-dimensional array, each row containing one read. For
            // each unique tag that has been merged together, add the sequence for
            // that tag into our array as many times as it originally occurred.
            //
            vector<int>::iterator j;
            vector<DNANSeq *> reads;

            for (j = mtag->utags.begin(); j != mtag->utags.end(); j++) {
                utag = unique[*j];

                for (uint k = 0; k < utag->count; k++) {
                    reads.push_back(utag->seq);
                }
            }

            //
            // Iterate over each column of the array and call the consensus base.
            //
            int row, col;
            int length = reads[0]->size();
            int height = reads.size();
            string con;
            map<char, int> nuc;
            map<char, int>::iterator max, n;
            DNANSeq *d;

            for (col = 0; col < length; col++) {
                nuc['A'] = 0;
                nuc['C'] = 0;
                nuc['G'] = 0;
                nuc['T'] = 0;
                nuc['N'] = 0;

                for (row = 0; row < height; row++) {
                    d = reads[row];
                    if (nuc.count((*d)[col]))
                        nuc[(*d)[col]]++;
                }

                //
                // Find the base with a plurality of occurances and call it.
                //
                max = nuc.end();

                for (n = nuc.begin(); n != nuc.end(); n++) {
                    if (n->first == 'N')
                        continue;
                    if (max == nuc.end() || n->second > max->second)
                        max = n;
                }
                con += max->second == 0 ? 'N' : max->first;

                //
                // Search this column for the presence of a SNP
                //
                if (invoke_model)
                    switch(model_type) {
                    case snp:
                        call_multinomial_snp(mtag, col, nuc, true);
                        break;
                    case bounded:
                        call_bounded_multinomial_snp(mtag, col, nuc, true);
                        break;
                    case fixed:
                        call_multinomial_fixed(mtag, col, nuc);
                        break;
                    }
            }

            if (invoke_model) {
                call_alleles(mtag, reads);

                if (model_type == fixed) {
                    //
                    // Mask nucleotides that are not fixed.
                    //
                    vector<SNP *>::iterator s;
                    for (s = mtag->snps.begin(); s != mtag->snps.end(); s++) {
                        if ((*s)->type == snp_type_unk)
                            con.replace((*s)->col, 1, "N");
                    }
                }
            }

            mtag->add_consensus(con.c_str());

            //
            // If SNPs were called at this locus but no alleles could be determined,
            // blacklist this tag. This can occur if there are two many uncalled bases
            // in the locus (Ns), such that haplotypes can't be consistently read
            // due to the presence of the Ns in the reads.
            //
            if (mtag->alleles.empty())
                for (uint j = 0; j < mtag->snps.size(); j++)
                    if (mtag->snps[j]->type == snp_type_het) {
                        mtag->blacklisted = 1;
                        break;
                    }
        }
    }

    return 0;
}

int call_alleles(MergedStack *mtag, vector<DNANSeq *> &reads) {
    int      row;
    int      height = reads.size();
    string   allele;
    char     base;
    vector<SNP *>::iterator snp;
    DNANSeq *d;

    for (row = 0; row < height; row++) {
        allele.clear();

        uint snp_cnt = 0;

        for (snp = mtag->snps.begin(); snp != mtag->snps.end(); snp++) {
            if ((*snp)->type != snp_type_het) continue;

            snp_cnt++;

            d    = reads[row];
            base = (*d)[(*snp)->col];

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
            mtag->alleles[allele]++;
    }

    return 0;
}
