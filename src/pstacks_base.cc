#include <iostream>
#include <iomanip>
#include <sstream>
#include <fstream>
#include <string>
#include <vector>
#include <map>

#include <zlib.h>

#include "constants.h"
#include "models.h"

#include "pstacks_base.h"

using namespace std;

int call_consensus(map<int, MergedStack *> &merged, map<int, PStack *> &unique, bool invoke_model) {
    cerr << "Identifying polymorphic sites and calling consensus sequences...";
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
                    case ::fixed:
                        call_multinomial_fixed(mtag, col, nuc);
                        break;
                    default:
                        DOES_NOT_HAPPEN;
                        break;
                    }
            }

            if (invoke_model) {
                call_alleles(mtag, reads);

                if (model_type == ::fixed) {
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
        }
    }

    size_t n_blacklisted = 0;
    for (const auto& mtag : merged)
        if (mtag.second->blacklisted)
            ++n_blacklisted;

    cerr << "done. (" << n_blacklisted << " loci were blacklisted.)\n";
    return 0;
}

void call_alleles(MergedStack *mtag, vector<DNANSeq *> &reads) {

    vector<SNP*> het_snps;
    for (SNP* snp : mtag->snps)
        if (snp->type == snp_type_het)
            het_snps.push_back(snp);

    if (!het_snps.empty()) {
        string allele;
        allele.reserve(het_snps.size());
        for (DNANSeq* r : reads) {
            allele.clear();
            for(SNP* hsnp : het_snps) {
                //
                // Check to make sure the nucleotide at the location of this SNP is
                // of one of the two possible states the multinomial model called.
                //
                char base = (*r)[hsnp->col];
                if (base != hsnp->rank_1 && base != hsnp->rank_2)
                    break;
                allele.push_back(base);
            }

            if (allele.length() == het_snps.size())
                mtag->alleles[allele]++;
        }

        if (mtag->alleles.empty())
            mtag->blacklisted = true;
    }
}

void
calc_coverage_distribution(const map<int, PStack*> &unique,
                           const map<int, MergedStack*> &merged,
                           double &mean, double &stdev, double &max)
{
    max   = 0.0;
    mean  = 0.0;
    stdev = 0.0;

    size_t not_blacklisted = 0;
    for (const pair<int, MergedStack*>& mtag : merged) {
        if (mtag.second->blacklisted)
            continue;

        ++not_blacklisted;

        double depth = 0.0;
        for (int utag_id : mtag.second->utags)
            depth += unique.at(utag_id)->count;

        mean += depth;
        if (depth > max)
            max = depth;
    }
    mean /= not_blacklisted;

    for (const pair<int, MergedStack*>& mtag : merged) {
        if (mtag.second->blacklisted)
            continue;

        double depth = 0.0;
        for (int utag_id : mtag.second->utags)
            depth += unique.at(utag_id)->count;

        stdev += pow(depth - mean, 2);
    }

    stdev /= not_blacklisted;
    stdev = sqrt(stdev);
}

int write_results(map<int, MergedStack *> &m,
                  map<int, PStack *> &u,
                  bool gzip,
                  bool paired_end
                  ) {
    map<int, MergedStack *>::iterator i;
    vector<char *>::iterator   j;
    vector<int>::iterator      k;
    vector<SNP *>::iterator    s;
    map<string, int>::iterator t;
    MergedStack *tag_1;
    PStack      *tag_2;
    stringstream sstr;

    //
    // Determine the names of the output files
    //
    string gz = (gzip ? ".gz" : "");
    string pe = (paired_end ? "_pe" : "");
    string tag_file = prefix_path + ".tags" + pe + ".tsv" + gz;
    string snp_file = prefix_path + ".snps" + pe + ".tsv" + gz;
    string all_file = prefix_path + ".alleles" + pe + ".tsv" + gz;
    string mod_file = prefix_path + ".models" + pe + ".tsv" + gz;

    //
    // Open the output files for writing.
    //
    gzFile   gz_tags, gz_snps, gz_alle, gz_mods;
    ofstream tags, snps, alle, mods;
    if (gzip) {
        gz_tags = gzopen(tag_file.c_str(), "wb");
        if (!gz_tags) {
            cerr << "Error: Unable to open gzipped tag file '" << tag_file << "': " << strerror(errno) << ".\n";
            exit(1);
        }
        #if ZLIB_VERNUM >= 0x1240
        gzbuffer(gz_tags, libz_buffer_size);
        #endif
        gz_mods = gzopen(mod_file.c_str(), "wb");
        if (!gz_mods) {
            cerr << "Error: Unable to open gzipped tag file '" << tag_file << "': " << strerror(errno) << ".\n";
            exit(1);
        }
        #if ZLIB_VERNUM >= 0x1240
        gzbuffer(gz_mods, libz_buffer_size);
        #endif
        gz_snps = gzopen(snp_file.c_str(), "wb");
        if (!gz_snps) {
            cerr << "Error: Unable to open gzipped snps file '" << snp_file << "': " << strerror(errno) << ".\n";
            exit(1);
        }
        #if ZLIB_VERNUM >= 0x1240
        gzbuffer(gz_snps, libz_buffer_size);
        #endif
        gz_alle = gzopen(all_file.c_str(), "wb");
        if (!gz_alle) {
            cerr << "Error: Unable to open gzipped alleles file '" << all_file << "': " << strerror(errno) << ".\n";
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
            cerr << "Error: Unable to open tag file for writing.\n";
            exit(1);
        }
        snps.open(snp_file.c_str());
        if (snps.fail()) {
            cerr << "Error: Unable to open SNPs file for writing.\n";
            exit(1);
        }
        alle.open(all_file.c_str());
        if (alle.fail()) {
            cerr << "Error: Unable to open allele file for writing.\n";
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
    log << "# pstacks version " << VERSION << "; generated on " << date << "\n";
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

    int id;

    int   wrote       = 0;
    int   blacklisted = 0;

    for (i = m.begin(); i != m.end(); i++) {
        tag_1 = i->second;

        //
        // Calculate the log likelihood of this merged stack.
        //
        tag_1->gen_matrix(u);
        tag_1->calc_likelihood();

        wrote++;

        if (tag_1->blacklisted) blacklisted++;

        // First write the consensus sequence
        sstr << "0" << "\t"
             << sql_id << "\t"
             << tag_1->id << "\t"
             << tag_1->loc.chr << "\t"
             << tag_1->loc.bp << "\t"
             << (tag_1->loc.strand == strand_plus ? "+" : "-") << "\t"
             << "consensus\t" << "\t\t"
             << tag_1->con << "\t"
             << tag_1->deleveraged << "\t"
             << tag_1->blacklisted << "\t"
             << tag_1->lumberjackstack << "\t"
             << tag_1->lnl << "\n";

        //
        // Write a sequence recording the output of the SNP model for each nucleotide.
        //
        sstr << "0" << "\t"
             << sql_id << "\t"
             << tag_1->id << "\t"
             << "\t"
             << "\t"
             << "\t"
             << "model\t" << "\t"
             << "\t";
        for (s = tag_1->snps.begin(); s != tag_1->snps.end(); s++) {
            switch((*s)->type) {
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

        // Now write out the components of each unique tag merged into this one.
        id = 0;
        for (k = tag_1->utags.begin(); k != tag_1->utags.end(); k++) {
            tag_2  = u[*k];
            for (j = tag_2->map.begin(); j != tag_2->map.end(); j++) {
                sstr << "0" << "\t" << sql_id << "\t" << tag_1->id << "\t\t\t\t" << "primary\t" << id << "\t" << *j << "\t" << tag_2->seq->seq() << "\t\t\t\t\n";
                if (gzip) gzputs(gz_tags, sstr.str().c_str()); else tags << sstr.str();
                sstr.str("");
            }
            id++;
        }

        //
        // Write out the model calls for each nucleotide in this locus.
        //
        for (s = tag_1->snps.begin(); s != tag_1->snps.end(); s++) {
            sstr << "0"          << "\t"
                 << sql_id       << "\t"
                 << tag_1->id    << "\t"
                 << (*s)->col    << "\t";

            switch((*s)->type) {
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

            sstr << std::fixed   << std::setprecision(2)
                 << (*s)->lratio << "\t"
                 << (*s)->rank_1 << "\t"
                 << (*s)->rank_2 << "\t\t\n";
        }

        if (gzip) gzputs(gz_snps, sstr.str().c_str()); else snps << sstr.str();
        sstr.str("");

        // Write the expressed alleles seen for the recorded SNPs and
        // the percentage of tags a particular allele occupies.
        //
        uint total = 0;
        for (int stack : tag_1->utags)
            total += u.at(stack)->count;
        char pct[id_len];
        for (t = tag_1->alleles.begin(); t != tag_1->alleles.end(); t++) {
            sprintf(pct, "%.2f", (((double)t->second/total) * 100));
            sstr << "0"       << "\t"
                 << sql_id    << "\t"
                 << tag_1->id << "\t"
                 << t->first  << "\t"
                 << pct       << "\t"
                 << t->second << "\n";
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

    cerr << "  Wrote " << wrote << " loci (of which " << blacklisted << " are blacklisted)." << endl;

    return 0;
}
