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

#ifndef __LOCUS_H__
#define __LOCUS_H__

#include <string.h>
#include <string>
using std::string;
#include <vector>
using std::vector;
#include <map>
using std::map;
#include <algorithm>
#include <utility>
using std::pair;
using std::make_pair;

#include "constants.h"
#include "stacks.h"

typedef struct match {
    uint        cat_id;
    allele_type cat_type;
    allele_type query_type;
    uint        dist;
} Match;

class Locus {
 public:
    int          id; // Locus ID
    int   sample_id; // Sample ID
    int       depth; // Stack depth
    char       *con; // Consensus sequence
    char     *model; // Model calls for each nucleotide
    uint        len; // Sequence length
    double      lnl; // Log likelihood of this locus

    //
    // Flags
    //
    bool blacklisted;
    bool deleveraged;
    bool lumberjackstack;

    vector<char *>          comp;   // Raw components in this stack.
    vector<char *>         reads;   // Sequence reads contributing to this locus.
    vector<uint>        comp_cnt;   // Counter for internal stacks merged into this locus.
    vector<read_type>  comp_type;   // Read types for reads contributing to this locus.
    PhyLoc                   loc;   // Physical genome location of this stack.
    vector<SNP *>           snps;   // Single Nucleotide Polymorphisms in this stack.
    map<string, int>     alleles;   // Map of the allelic configuration of SNPs in this stack along with the count of each
    vector<pair<allele_type, string> > strings; // Strings for matching (representing the various allele combinations)

    Locus()  { 
	id              = 0; 
	sample_id       = 0; 
	depth           = 0; 
	model           = NULL;
	con             = NULL; 
	len             = 0;
	lnl             = 0.0;
	blacklisted     = false;
        deleveraged     = false;
	lumberjackstack = false;
    }
    virtual ~Locus() { 
        delete [] con; 
	delete [] model;
        for (uint i = 0; i < snps.size(); i++)
            delete snps[i];
        for (uint i = 0; i < comp.size(); i++)
            delete [] comp[i];
        for (uint i = 0; i < reads.size(); i++)
            delete [] reads[i];
    }
    uint sort_bp(uint k = 0);
    int snp_index(uint);
    int add_consensus(const char *);
    virtual int populate_alleles();
};

//
// Query Locus Class
//
class QLocus : public Locus {
 public:
    vector<Match *> matches;   // Matching tags found for the catalog.

    QLocus(): Locus() {}
    ~QLocus();

    int add_match(int, allele_type, allele_type, int);
    int add_match(int, allele_type);
};

//
// Catalog Locus Class, for use in cstacks, contains catalog loci and records the
// constiuent tags this locus was built from.
//
class CLocus : public Locus {
 public:
    vector<pair<int, int> > sources;   // Sample/ID pairs for the sources contributing to this catalog entry

    int merge_snps(QLocus *);
    int reduce_alleles(set<string> &);
};

//
// Catalog Summary Locus Class; used in genotypes and populations, records a catalog 
// locus with summary information derived from individuals in the population. 
//
class CSLocus : public Locus {
public:
    CSLocus() : Locus() { 
	this->f          = 0.0; 
	this->cnt  = 0; 
	this->hcnt       = 0; 
	this->gcnt       = 0; 
	this->trans_gcnt = 0; 
	this->chisq      = 1.0;
	this->confounded_cnt = 0;
    };
    string annotation;
    string marker;
    string uncor_marker;
    map<string, int> hap_cnts;    // Counts of each observed haplotype for this locus in the population.
    double f;                     // Inbreeder's coefficient
    map<string, string> gmap;     // Observed haplotype to genotype map for this locus.
    int confounded_cnt;           // Number of samples/progeny containing confounded loci (more than one 
                                  //   locus from an individual sample matches this catalog locus).
    int hcnt;                     // Number of samples/progeny containing a haplotype for this locus.
    int cnt;                      // Number of samples/progeny containing data for this locus.
    int gcnt;                     // Number of progeny containing a valid genotype.
    int trans_gcnt;               // Number of progeny containing a valid 
                                  //   genotype, translated for a particular map type.
    double chisq;             // Chi squared p-value testing the null hypothesis of no segregation distortion.
};

bool bp_compare(Locus *, Locus *);

#endif // __LOCUS_H__
