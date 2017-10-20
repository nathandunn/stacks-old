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
// sql_utilities.h -- template routines to read and write Stacks SQL file formats.
//
#ifndef __SQL_UTILITIES_H__
#define __SQL_UTILITIES_H__

#include <unordered_set>

#include "input.h"
#include "utils.h"

//
// The expected number of tab-separated fields in our SQL input files.
//
const uint num_tags_fields    = 14;
const uint num_snps_fields    = 10;
const uint num_alleles_fields =  6;
const uint num_matches_fields =  9;

void load_catalog_matches(string sample,  vector<CatMatch *> &matches, bool verbose=true);
int load_model_results(string sample,  map<int, ModRes *> &modres);
int load_snp_calls(string sample,  map<int, SNPRes *> &snpres);

// retrieve_bijective_sloci()
// ----------
// Returns pairs of (sample locus ID, catalog locus ID) for which there is a bijective relation.
vector<pair<int, int> > retrieve_bijective_loci(const vector<CatMatch*>& matches);
vector<pair<int, int> > retrieve_bijective_loci(const vector<pair<int,int>>& sloc_cloc_id_pairs);

template <class LocusT>
int
load_loci(const string& sample,  map<int, LocusT *> &loci, int store_reads, bool load_all_model_calls, bool &compressed, bool verbose=true)
{
    using namespace std;

    LocusT        *c;
    SNP           *snp;
    string         f;
    char          *cmp;
    const char    *p, *q;
    int            len;
    vector<string> parts;
    set<int>       blacklisted;
    long int       line_num;
    ifstream       fh;
    gzFile         gz_fh;

    char *line      = new char[max_len];
    int   size      = max_len;
    bool  gzip      = false;
    bool  open_fail = true;
    int   fh_status = 1;

    // //
    // // First, try to parse the models file to pull in the consensus sequence and model string
    // // for each locus. If the models file is not available or we are requested to store the
    // // reads from each stack, fall back to the tags file.
    // //
    // if (!store_reads) {
    //     f = sample + ".models.tsv";
    //     fh.open(f.c_str(), ifstream::in);
    //     open_fail = fh.fail() ? true : false;
    // }

    // if (!store_reads && open_fail) {
    //     //
    //     // Test for a gzipped MODELs file.
    //     //
    //     f = sample + ".models.tsv.gz";
    //     gz_fh = gzopen(f.c_str(), "rb");
    //     if (!gz_fh) {
    //         open_fail = true;
    //     } else {
    //         open_fail = false;
    //         #if ZLIB_VERNUM >= 0x1240
    //         gzbuffer(gz_fh, libz_buffer_size);
    //         #endif
    //         gzip = true;
    //     }
    // }

    if (open_fail) {
        //
        // Test for a TAGs file.
        //
        f = sample + ".tags.tsv";
        fh.open(f.c_str(), ifstream::in);
        if (fh.fail())
            open_fail = true;
        else
            open_fail = false;
    }

    if (open_fail) {
        //
        // Test for a gzipped TAGs file.
        //
        f = sample + ".tags.tsv.gz";
        gz_fh = gzopen(f.c_str(), "rb");
        if (!gz_fh) {
            open_fail = true;
        } else {
            open_fail = false;
            #if ZLIB_VERNUM >= 0x1240
            gzbuffer(gz_fh, libz_buffer_size);
            #endif
            gzip = true;
        }
    }

    if (open_fail) {
        cerr << " Unable to open '" << sample << "'\n";
        return 0;
    }

    if (verbose)
        cerr << "  Parsing " << f.c_str() << "\n";

    uint id;

    line_num = 1;
    while (fh_status) {
        fh_status = (gzip == true) ? read_gzip_line(gz_fh, &line, &size) : read_line(fh, &line, &size);

        if (!fh_status && strlen(line) == 0)
            continue;

        if (is_comment(line)) continue;

        parse_tsv(line, parts);

        if (parts.size() != num_tags_fields) {
            cerr << "Error parsing " << f.c_str() << " at line: " << line_num << ". (" << parts.size() << " fields).\n";
            return 0;
        }

        id = atoi(parts[2].c_str());

        if (parts[6] != "consensus") {
            if (blacklisted.count(id)) continue;

            //
            // Make sure this locus has already been defined (consensus sequence SHOULD always
            // be specified first in the file for a particular locus).
            //
            if (loci.count(id) > 0) {
                //
                // Read the model sequence, a series of letters specifying if the model called a
                // homozygous base (O), a heterozygous base (E), or if the base type was unknown (U).
                //
                if (parts[6] == "model") {
                    loci[id]->model = new char[parts[9].length() + 1];
                    strcpy(loci[id]->model, parts[9].c_str());

                } else {
                    //
                    // Otherwise, we expect a primary or secondary read, record these if specified.
                    //
                    loci[id]->depth++;

                    if (store_reads >= 1) {
                        if (store_reads >= 2) {
                            // Load the actual sequences (otherwise don't).
                            char *read = new char[parts[9].length() + 1];
                            strcpy(read, parts[9].c_str());
                            loci[id]->reads.push_back(read);
                        }

                        char *read_id = new char[parts[8].length() + 1];
                        strcpy(read_id, parts[8].c_str());
                        loci[id]->comp.push_back(read_id);
                        //
                        // Store the internal stack number for this read.
                        //
                        loci[id]->comp_cnt.push_back(atoi(parts[7].c_str()));

                        //
                        // Store the read type.
                        //
                        if (parts[6] == "primary")
                            loci[id]->comp_type.push_back(primary);
                        else
                            loci[id]->comp_type.push_back(secondary);
                    }
                }

                continue;
            } else {
                cerr << "Error parsing " << f.c_str() << " at line: " << line_num << ". (stack " << id << " does not exist).\n";
                return 0;
            }
        }

        //
        // Do not include blacklisted tags in the catalog. They are tags that are composed
        // of noise and/or repetitive sequence.
        //
        if (parts[11] == "1") {
            blacklisted.insert(id);
            continue;
        }

        c = new LocusT;
        c->sample_id = atoi(parts[1].c_str());
        c->id        = id;
        c->add_consensus(parts[9].c_str());

        //
        // Read in the flags
        //
        c->deleveraged     = (parts[10] == "1" ? true : false);
        c->lumberjackstack = (parts[12] == "1" ? true : false);

        //
        // Parse the physical genome location of this locus.
        //
        c->loc.set(parts[3].c_str(), atoi(parts[4].c_str()), (parts[5] == "+" ? strand_plus : strand_minus));

        //
        // Parse the components of this stack (either the Illumina ID, or the catalog constituents)
        //
        q = parts[8].c_str();
        while (*q != '\0') {
            for (p = q; *q != ',' && *q != '\0'; q++);
            len = q - p;
            cmp = new char[len + 1];
            strncpy(cmp, p, len);
            cmp[len] = '\0';
            c->comp.push_back(cmp);
            if (*q != '\0') q++;
        }

        loci[c->id] = c;

        line_num++;
    }

    if (gzip)
        gzclose(gz_fh);
    else
        fh.close();

    //
    // Next, parse the SNP file and load model calls.
    //
    gzip      = false;
    fh_status = 1;
    line_num  = 1;

    f = sample + ".snps.tsv";
    fh.open(f.c_str(), ifstream::in);
    if (fh.fail()) {
        //
        // Test for a gzipped file.
        //
        f = sample + ".snps.tsv.gz";
        gz_fh = gzopen(f.c_str(), "rb");
        if (!gz_fh) {
            cerr << " Unable to open '" << sample << "'\n";
            return 0;
        }
        #if ZLIB_VERNUM >= 0x1240
        gzbuffer(gz_fh, libz_buffer_size);
        #endif
        gzip       = true;
        compressed = true;
    }
    if (verbose)
        cerr << "  Parsing " << f.c_str() << "\n";

    while (fh_status) {
        fh_status = (gzip == true) ? read_gzip_line(gz_fh, &line, &size) : read_line(fh, &line, &size);

        if (!fh_status && strlen(line) == 0)
            continue;

        if (is_comment(line)) continue;

        parse_tsv(line, parts);

        if (parts.size() != num_snps_fields && parts.size() != num_snps_fields - 2) {
            cerr << "Error parsing " << f.c_str() << " at line: " << line_num << ". (" << parts.size() << " fields).\n";
            return 0;
        }

        id = atoi(parts[2].c_str());

        if (blacklisted.count(id))
            continue;

        //
        // Only load heterozygous model calls.
        //
        if (load_all_model_calls == false && parts[4] != "E")
            continue;

        snp         = new SNP;
        snp->col    = atoi(parts[3].c_str());
        snp->lratio = atof(parts[5].c_str());
        snp->rank_1 = parts[6].at(0);
        snp->rank_2 = parts[7].at(0) == '-' ? 0 : parts[7].at(0);

        if (parts[4] == "E")
            snp->type = snp_type_het;
        else if (parts[4] == "O")
            snp->type = snp_type_hom;
        else
            snp->type = snp_type_unk;

        if (parts.size() == 10) {
            if (parts[8].length() == 0 || parts[8].at(0) == '-')
                snp->rank_3 = 0;
            else
                snp->rank_3 = parts[8].at(0);

            if (parts[9].length() == 0 || parts[9].at(0) == '-')
                snp->rank_4 = 0;
            else
                snp->rank_4 = parts[9].at(0);
        }

        if (loci.count(id) > 0) {
            loci[id]->snps.push_back(snp);
        } else {
            cerr << "Error parsing " << f.c_str() << " at line: " << line_num << ". SNP asks for nonexistent locus with ID: " << id << "\n";
            return 0;
        }

        line_num++;
    }

    if (gzip)
        gzclose(gz_fh);
    else
        fh.close();

    //
    // Finally, parse the Alleles file
    //
    gzip      = false;
    fh_status = 1;
    line_num  = 1;

    f = sample + ".alleles.tsv";
    fh.open(f.c_str(), ifstream::in);
    if (fh.fail()) {
        //
        // Test for a gzipped file.
        //
        f = sample + ".alleles.tsv.gz";
        gz_fh = gzopen(f.c_str(), "rb");
        if (!gz_fh) {
            cerr << " Unable to open '" << sample << "'\n";
            return 0;
        }
        #if ZLIB_VERNUM >= 0x1240
        gzbuffer(gz_fh, libz_buffer_size);
        #endif
        gzip       = true;
        compressed = true;
    }
    if (verbose)
        cerr << "  Parsing " << f.c_str() << "\n";

    while (fh_status) {
        fh_status = (gzip == true) ? read_gzip_line(gz_fh, &line, &size) : read_line(fh, &line, &size);

        if (!fh_status && strlen(line) == 0)
            continue;

        if (is_comment(line)) continue;

        parse_tsv(line, parts);

        if (parts.size() != num_alleles_fields) {
            cerr << "Error parsing " << f.c_str() << " at line: " << line_num << ". (" << parts.size() << " fields).\n";
            return 0;
        }

        id = atoi(parts[2].c_str());

        if (blacklisted.count(id))
            continue;

        if (loci.count(id) > 0) {
            loci[id]->alleles[parts[3]] = atoi(parts[5].c_str());
        } else {
            cerr << "Error parsing " << f.c_str() << " at line: " << line_num << ". SNP asks for nonexistent locus with ID: " << id << "\n";
            return 0;
        }

        line_num++;
    }

    if (gzip)
        gzclose(gz_fh);
    else
        fh.close();

    //
    // Populate the strings member with the sequence for each allele for each Locus.
    //
    typename map<int, LocusT *>::iterator i;
    for (i = loci.begin(); i != loci.end(); i++)
        i->second->populate_alleles();

    free(line);

    return 1;
}

template <class LocusT>
int dump_loci(map<int, LocusT *> &u) {
    typename map<int, LocusT *>::iterator i;
    vector<SNP *>::iterator      s;

    for (i = u.begin(); i != u.end(); i++) {

        cerr << "Locus ID:    " << i->second->id << "\n"
             << "  Consensus: " << i->second->con << "\n"
             << "  Genomic Location: " << i->second->loc.chr << "; " << i->second->loc.bp +1 << "bp\n"
             << "  SNPs:\n";

        for (s = i->second->snps.begin(); s != i->second->snps.end(); s++)
            cerr << "    Col: " << (*s)->col << " rank 1: " << (*s)->rank_1 << " rank 2: " << (*s)->rank_2 << "\n";

        cerr << "\n";
    }

    return 0;
}

#endif // __SQL_UTILITIES_H__
