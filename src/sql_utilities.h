// -*-mode:c++; c-style:k&r; c-basic-offset:4;-*-
//
// Copyright 2010, Julian Catchen <jcatchen@uoregon.edu>
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
// Julian Catchen
// jcatchen@uoregon.edu
// University of Oregon
//
#ifndef __SQL_UTILITIES_H__
#define __SQL_UTILITIES_H__

#include "input.h"

//
// The expected number of tab-separated fields in our SQL input files.
//
const uint num_tags_fields    = 13;
const uint num_snps_fields    =  9;
const uint num_alleles_fields =  6;
const uint num_matches_fields =  7;

template <class LocusT>
int load_loci(string sample,  map<int, LocusT *> &loci, bool store_reads) {
    LocusT        *c;
    SNP           *snp;
    string         f;
    char          *line, *cmp;
    const char    *p, *q;
    int            len, size;
    vector<string> parts;
    set<int>       blacklisted;
    long int       line_num;
    ifstream       fh;

    line = (char *) malloc(sizeof(char) * max_len);
    size = max_len;

    // 
    // First, parse the tag file and pull in the consensus sequence
    // for each Radtag.
    //
    f = sample + ".tags.tsv";
    fh.open(f.c_str(), ifstream::in);

    if (fh.fail()) {
        cerr << " Unable to open " << f.c_str() << "\n";
        return 0;
    } else {
        cerr << "  Parsing " << f.c_str() << "\n";
    }

    uint id;

    line_num = 0;
    while (fh.good()) {
        read_line(fh, &line, &size);

	if (!fh.good() && strlen(line) == 0)
	    continue;

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

		    if (store_reads) {
			char *read = new char[parts[9].length() + 1];
			strcpy(read, parts[9].c_str());
			loci[id]->reads.push_back(read);
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
        // Parse the physical genome location of this locus.
        //
	c->loc.set(parts[3].c_str(), atoi(parts[4].c_str()), (parts[5] == "+" ? plus : minus));

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

    fh.close();

    // 
    // Next, parse the SNP file
    //
    f = sample + ".snps.tsv";
    fh.open(f.c_str(), ifstream::in);

    if (fh.fail()) {
        cerr << " Unable to open " << f.c_str() << "\n";
        return 0;
    } else {
        cerr << "  Parsing " << f.c_str() << "\n";
    }

    line_num = 0;
    while (fh.good()) {
        read_line(fh, &line, &size);

	if (!fh.good() && strlen(line) == 0)
	    continue;

	parse_tsv(line, parts);

        if (parts.size() != num_snps_fields && parts.size() != num_snps_fields - 2) {
            cerr << "Error parsing " << f.c_str() << " at line: " << line_num << ". (" << parts.size() << " fields).\n";
            return 0;
        }

        id = atoi(parts[2].c_str());

	if (blacklisted.count(id))
	    continue;

	snp         = new SNP;
	snp->col    = atoi(parts[3].c_str());
	snp->lratio = atof(parts[4].c_str());
	snp->rank_1 = parts[5].at(0);
	snp->rank_2 = parts[6].at(0);

	if (parts.size() == 9) {
	    snp->rank_3 = parts[7].length() == 0 ? 0 : parts[7].at(0);
	    snp->rank_4 = parts[8].length() == 0 ? 0 : parts[8].at(0);
	}

        if (loci.count(id) > 0) {
            loci[id]->snps.push_back(snp);
        } else {
            cerr << "Error parsing " << f.c_str() << " at line: " << line_num << ". SNP asks for nonexistent locus with ID: " << id << "\n";
            return 0;
        }

        line_num++;
    }

    fh.close();

    // 
    // Finally, parse the Alleles file
    //
    f = sample + ".alleles.tsv";
    fh.open(f.c_str(), ifstream::in);

    if (fh.fail()) {
        cerr << " Unable to open " << f.c_str() << "\n";
        return 0;
    } else {
        cerr << "  Parsing " << f.c_str() << "\n";
    }

    line_num = 0;
    while (fh.good()) {
        read_line(fh, &line, &size);

	if (!fh.good() && strlen(line) == 0)
	    continue;

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

    //
    // Populate the strings member with the sequence for each allele for each Locus.
    //
    typename map<int, LocusT *>::iterator i;
    for (i = loci.begin(); i != loci.end(); i++)
        i->second->populate_alleles();


    fh.close();

    delete [] line;

    return 1;
}

template <class LocusT>
int dump_loci(map<int, LocusT *> &u) {
    typename map<int, LocusT *>::iterator i;
    vector<SNP *>::iterator      s;

    for (i = u.begin(); i != u.end(); i++) {

	cerr << "Locus ID:    " << i->second->id << "\n"
	     << "  Consensus: " << i->second->con << "\n"
             << "  Genomic Location: " << i->second->loc.chr << "; " << i->second->loc.bp << "bp\n"
	     << "  SNPs:\n";

	for (s = i->second->snps.begin(); s != i->second->snps.end(); s++)
	    cerr << "    Col: " << (*s)->col << " rank 1: " << (*s)->rank_1 << " rank 2: " << (*s)->rank_2 << "\n";

	cerr << "\n";
    }

    return 0;
}

int load_catalog_matches(string sample,  vector<CatMatch *> &matches) {
    CatMatch      *m;
    string         f;
    char           line[max_len];
    vector<string> parts;
    long int       line_num;
    ifstream       fh;

    f = sample + ".matches.tsv";
    fh.open(f.c_str(), ifstream::in);

    if (fh.fail()) {
        cerr << " Unable to open " << f.c_str() << "\n";
        return 0;
    } else {
        cerr << "  Parsing " << f.c_str() << "\n";
    }

    line_num = 0;
    while (fh.good()) {
	fh.getline(line, max_len);
        line_num++;

	if (!fh.good() && strlen(line) == 0)
	    continue;

	parse_tsv(line, parts);

        if (parts.size() != num_matches_fields) {
            cerr << "Error parsing " << f.c_str() << " at line: " << line_num << ". (" << parts.size() << " fields).\n";
            return 0;
        }

	m = new CatMatch;
	m->batch_id  = atoi(parts[1].c_str());
	m->cat_id    = atoi(parts[2].c_str());
        m->sample_id = atoi(parts[3].c_str());
	m->tag_id    = atoi(parts[4].c_str());
	m->haplotype = new char[parts[5].length() + 1];
	strcpy(m->haplotype, parts[5].c_str());
	m->depth     = atoi(parts[6].c_str());

	matches.push_back(m);
    }

    fh.close();

    return 0;
}

int load_model_results(string sample,  map<int, ModRes *> &modres) {
    string         f;
    char          *line;
    int            size;
    vector<string> parts;
    long int       line_num;
    ifstream       fh;

    line = (char *) malloc(sizeof(char) * max_len);
    size = max_len;

    // 
    // First, parse the tag file and pull in the consensus sequence
    // for each Radtag.
    //
    f = sample + ".tags.tsv";
    fh.open(f.c_str(), ifstream::in);

    if (fh.fail()) {
        cerr << " Unable to open " << f.c_str() << "\n";
        return 0;
    } else {
        cerr << "  Parsing " << f.c_str() << "\n";
    }

    ModRes *mod;
    uint    tag_id, samp_id;

    line_num = 0;
    while (fh.good()) {
        read_line(fh, &line, &size);
        line_num++;

	if (!fh.good() && strlen(line) == 0)
	    continue;

	parse_tsv(line, parts);

        if (parts.size() != num_tags_fields) {
            cerr << "Error parsing " << f.c_str() << " at line: " << line_num << ". (" << parts.size() << " fields).\n";
            return 0;
        }

	//
	// Read the model sequence, a series of letters specifying if the model called a
	// homozygous base (O), a heterozygous base (E), or if the base type was unknown (U).
	//
	if (parts[6] != "model") continue;

	samp_id = atoi(parts[1].c_str()); 
	tag_id  = atoi(parts[2].c_str());
	mod     = new ModRes(samp_id, tag_id, parts[9].c_str());

	modres[tag_id] = mod;
    }

    fh.close();

    delete [] line;

    return 1;
}

#endif // __SQL_UTILITIES_H__
