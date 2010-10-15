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

// -*-mode:c++; c-style:k&r; c-basic-offset:4;-*-
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
const uint num_tags_fields    = 12;
const uint num_snps_fields    =  7;
const uint num_alleles_fields =  6;

template <class LocusT>
int load_loci(string &sample,  map<int, LocusT *> &loci) {
    LocusT        *c;
    SNP           *snp;
    string         f;
    char           line[max_len];
    vector<string> parts;
    set<int>       blacklisted;
    long int       line_num;
    ifstream       fh;

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
	fh.getline(line, max_len);

	if (!fh.good() && strlen(line) == 0)
	    continue;

	parse_tsv(line, parts);

        if (parts.size() != num_tags_fields) {
            cerr << "Error parsing " << f.c_str() << " at line: " << line_num << ". (" << parts.size() << " fields).\n";
            return 0;
        }

        id = atoi(parts[2].c_str());

	if (parts[5] != "consensus") {
            if (blacklisted.count(id)) {
                continue;
            } else if (loci.count(id) > 0) {
                loci[id]->depth++;
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
	if (parts[10] == "1") {
	    blacklisted.insert(id);
	    continue;
	}

	c = new LocusT;
        c->sample_id = atoi(parts[1].c_str());
	c->id        = id;
	c->add_consensus(parts[8].c_str());

        //
        // Parse the physical genome location of this locus.
        //
        strncpy(c->loc.chr, parts[3].c_str(), id_len);
        c->loc.chr[id_len - 1] = '\0';
        c->loc.bp = atoi(parts[4].c_str());

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
	fh.getline(line, max_len);

	if (!fh.good() && strlen(line) == 0)
	    continue;

	parse_tsv(line, parts);

        if (parts.size() != num_snps_fields) {
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
	fh.getline(line, max_len);

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

#endif // __SQL_UTILITIES_H__
