// -*-mode:c++; c-style:k&r; c-basic-offset:4;-*-
//
// Copyright 2010-2016, Julian Catchen <jcatchen@illinois.edu>
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

#include "log_utils.h"

#include "pstacks.h"

using namespace std;
//
// Global variables to hold command-line options.
//
FileT  in_file_type;
string in_file;
string out_path;
int    sql_id        = -1;
int    min_stack_cov = 3;
double max_clipped   = 0.15;
bool   keep_sec_alns = false;
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

    report_options(cerr);

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

    map<int, PStack *> unique; // Unique {sequence, alignment position} combinations.
    HashMap*           radtags = new HashMap();
    load_radtags(in_file, *radtags);
    reduce_radtags(*radtags, unique);
    for (auto& stack : *radtags)
        for (Seq* read : stack.second)
            delete read;
    delete radtags;

    map<int, MergedStack *> merged;
    populate_merged_tags(unique, merged);

    delete_low_cov_loci(merged, unique);

    // Call the consensus sequence again, now that remainder tags have been merged.
    call_consensus(merged, unique, true);

    write_results(merged, unique);

    cerr << "pstacks is done.\n";

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

int write_results(map<int, MergedStack *> &m, map<int, PStack *> &u) {
    map<int, MergedStack *>::iterator i;
    vector<char *>::iterator   j;
    vector<int>::iterator      k;
    vector<SNP *>::iterator    s;
    map<string, int>::iterator t;
    MergedStack *tag_1;
    PStack      *tag_2;
    stringstream sstr;

    //
    // Parse the input file name to create the output files
    //
    size_t pos_1 = in_file.find_last_of("/");
    size_t pos_2 = in_file.find_last_of(".");
    string prefix = out_path + in_file.substr(pos_1 + 1, (pos_2 - pos_1 - 1));
    string tag_file = prefix + ".tags.tsv";
    string snp_file = prefix + ".snps.tsv";
    string all_file = prefix + ".alleles.tsv";
    string mod_file = prefix + ".models.tsv";

    cerr << "Writing tags, models, snps, and alleles files...\n";

    bool gzip = (in_file_type == FileT::bam) ? true : false;
    if (gzip) {
        tag_file += ".gz";
        snp_file += ".gz";
        all_file += ".gz";
        mod_file += ".gz";
    }

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

    char *buf; // = new char[m.begin()->second->len + 1];

    for (i = m.begin(); i != m.end(); i++) {
        tag_1 = i->second;

        float total = 0;
        for (k = tag_1->utags.begin(); k != tag_1->utags.end(); k++)
            total += u[*k]->count;

        //
        // Calculate the log likelihood of this merged stack.
        //
        tag_1->gen_matrix(u);
        tag_1->calc_likelihood();

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
            buf = tag_2->seq->seq();

            for (j = tag_2->map.begin(); j != tag_2->map.end(); j++) {
                sstr << "0" << "\t" << sql_id << "\t" << tag_1->id << "\t\t\t\t" << "primary\t" << id << "\t" << *j << "\t" << buf << "\t\t\t\t\n";
                if (gzip) gzputs(gz_tags, sstr.str().c_str()); else tags << sstr.str();
                sstr.str("");
            }
            id++;
            delete [] buf;
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
        char pct[id_len];
        for (t = tag_1->alleles.begin(); t != tag_1->alleles.end(); t++) {
            sprintf(pct, "%.2f", ((t->second/total) * 100));
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

    return 0;
}

void populate_merged_tags(map<int, PStack *> &unique, map<int, MergedStack *> &merged) {
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

    double mean;
    double stdev;
    double max;
    calc_coverage_distribution(unique, merged, mean, stdev, max);

    cerr << "Created " << merged.size() << " loci; mean coverage is " << mean << " (stdev: " << stdev << ", max: " << max << ").\n";
}

void delete_low_cov_loci(map<int, MergedStack *>& merged, const map<int, PStack*>& unique) {

    size_t n_deleted = 0;
    size_t n_reads = 0;
    vector<int> to_erase;

    for (auto& mtag : merged) {
        int depth = 0;
        for (int utag_id : mtag.second->utags)
            depth += unique.at(utag_id)->count;

        if (depth < min_stack_cov) {
            delete mtag.second;
            to_erase.push_back(mtag.first);

            n_deleted++;
            n_reads += depth;
        }
    }

    for (int id : to_erase)
        merged.erase(id);

    cerr << "Discarded " << n_deleted << " low coverage loci comprising " << n_reads << " reads.\n"
         << "Now working with " << merged.size() << " loci.\n";
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

    return 0;
}

//
// We expect tags to have already been aligned to a reference genome. Therefore, the tags
// are identified by their chromosome and basepair location.
//
void load_radtags(string in_file, HashMap &radtags) {
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

    cerr << "Reading alignments...\n";

    int primary_kept   = 0;
    int primary_disc   = 0;
    int secondary_kept = 0;
    int secondary_disc = 0;
    int supplementary  = 0;
    int unmapped       = 0;
    while ((fh->next_seq(c)) != 0) {

        switch (c.aln_type) {
        case AlnT::null:
            unmapped++;
            continue;
            break;
        case AlnT::primary:
            if (c.pct_clipped <= max_clipped) {
                primary_kept++;
            } else {
                primary_disc++;
                continue;
            }
            break;
        case AlnT::secondary:
            if (keep_sec_alns && c.pct_clipped <= max_clipped) {
                secondary_kept++;
            } else {
                secondary_disc++;
                continue;
            }
            break;
        case AlnT::supplementary:
            supplementary++;
            continue;
            break;
        default:
            // assert(false);
            break;
        }

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
    }

    cerr << "Done reading alignment records:\n"
         << "  Kept " << primary_kept << " primary alignments\n"
         << "  Skipped " << primary_disc << " (" << as_percentage((double) primary_disc / (primary_disc+primary_kept))
         << ") excessively soft-clipped primary alignments\n"
         << "  Skipped " << secondary_disc << " secondary alignments\n"
         << "  Skipped " << supplementary << " supplementary alignments\n"
         << "  Skipped " << unmapped << " unmapped reads\n";

    if(keep_sec_alns)
        cerr << "  Kept " << secondary_kept << " secondary alignments\n";

    if (radtags.empty()) {
        cerr << "Error: No input.\n";
        throw std::exception();
    }

    cerr << "Collapsed reads into " << radtags.size() << " stacks.\n";

    delete fh;
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
    int c;

    while (1) {
        static struct option long_options[] = {
            {"help",         no_argument,       NULL, 'h'},
            {"version",      no_argument,       NULL, 1000},
            {"infile_type",  required_argument, NULL, 't'},
            {"file",         required_argument, NULL, 'f'},
            {"outpath",      required_argument, NULL, 'o'},
            {"id",           required_argument, NULL, 'i'},
            {"min_cov",      required_argument, NULL, 'm'},
            {"max_clipped",  required_argument, NULL, 1001},
            {"keep_sec_aln", required_argument, NULL, 'k'},
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

        c = getopt_long(argc, argv, "hkvOT:a:A:L:U:f:o:i:e:p:m:s:f:t:", long_options, &option_index);

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
            break;
        case 'm':
            min_stack_cov = atoi(optarg);
            break;
        case 1001: //max_clipped
            max_clipped = is_double(optarg);
            if (max_clipped > 1)
                max_clipped = max_clipped / 100;

            if (max_clipped < 0 || max_clipped > 1.0) {
                cerr << "Unable to parse the maximum clipped proportion.\n";
                help();
            }
            break;
        case 'k':
            keep_sec_alns = true;
            break;
        case 'e':
            barcode_err_freq = atof(optarg);
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
        case 1000:
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

    if (in_file.empty()) {
        cerr << "You must specify an input file.\n";
        help();
    }

    if (in_file_type == FileT::unknown) {
        in_file_type = guess_file_type(in_file);
        if (in_file_type == FileT::unknown) {
            cerr << "Unable to recongnize the extention of file '" << in_file << "'.\n";
            help();
        }
    }

    if (sql_id < 0) {
        cerr << "A sample ID must be provided.\n";
        help();
    }

    if (out_path.length() == 0)
        out_path = ".";

    if (out_path.at(out_path.length() - 1) != '/')
        out_path += "/";

    if (model_type == ::fixed && barcode_err_freq == 0) {
        cerr << "You must specify the barcode error frequency.\n";
        help();
    }

    return 0;
}

void version() {
    std::cerr << "pstacks " << VERSION << "\n\n";

    exit(0);
}

void help() {
    std::cerr << "pstacks " << VERSION << "\n"
              << "pstacks -f file_path -i id -o path [-m min_cov] [-p num_threads]" << "\n"
              << "  f: input file path.\n"
              << "  i: a unique integer ID for this sample.\n"
              << "  o: output directory.\n"
              << "  m: minimum depth of coverage to report a stack (default 3).\n"
              << "  p: enable parallel execution with num_threads threads.\n"
              << "  t: input file Type. Supported types: bam, sam, bowtie (default: guess).\n"
              << "  --max_clipped <float>: alignments with more than this fraction of soft-clipped bases are discarded (default 15%).\n"
              << "  --keep_sec_alns: keep secondary alignments (default: false, only keep primary alignments).\n"
              << "\n"
              << "  Model options:\n"
              << "    --model_type <type>: either 'snp' (default), 'bounded', or 'fixed'\n"
              << "       For the SNP or Bounded SNP model:\n"
              << "       --alpha <num>: chi square significance level required to call a heterozygote or homozygote, either 0.1, 0.05 (default), 0.01, or 0.001.\n"
              << "       For the Bounded SNP model:\n"
              << "       --bound_low <num>: lower bound for epsilon, the error rate, between 0 and 1.0 (default 0).\n"
              << "       --bound_high <num>: upper bound for epsilon, the error rate, between 0 and 1.0 (default 1).\n"
              << "       For the Fixed model:\n"
              << "       --bc_err_freq <num>: specify the barcode error frequency, between 0 and 1.0.\n";

    exit(0);
}

void report_options(std::ostream& os) {
    os << "Alignments file: " << in_file << "\n"
         << "Output directory: " << out_path << "\n"
         << "Sample ID: " << sql_id << "\n"
         << "Min locus depth: " << min_stack_cov << "\n"
         << "Max clipped proportion: " << max_clipped << "\n";

    // Model.
    if (model_type == snp) {
        os << "Model: snp\n"
           << "Model alpha: " << alpha << "\n";
    } else if (model_type == bounded) {
        os << "Model: snp\n"
           << "Model alpha: " << alpha << "\n"
           << "Model lower bound: " << bound_low << "\n"
           << "Model higher bound: " << bound_high << "\n";
    } else if (model_type == ::fixed) {
        os << "Model: fixed\n"
           << "Model barcode err. prob.: " << barcode_err_freq << "\n";
    }
}
