// -*-mode:c++; c-style:k&r; c-basic-offset:4;-*-
//
// Copyright 2013-2017, Julian Catchen <jcatchen@illinois.edu>
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

#ifndef __CATALOG_UTILS_H__
#define __CATALOG_UTILS_H__

#include <vector>
#include <string>
#include <map>
#include <set>

#include "constants.h"
#include "stacks.h"
#include "locus.h"
#include "PopMap.h"
#include "PopSum.h"
#include "Vcf.h"

// find_catalogs()
// Looks for catalog files in the given directory and returns the associated ID(s).
vector<int> find_catalogs(const string& dir_path);

int check_whitelist_integrity(map<int, CSLocus *> &, map<int, set<int> > &);
int reduce_catalog(map<int, CSLocus *> &, set<int> &, set<int> &);
int reduce_catalog(map<int, CSLocus *> &, map<int, set<int> > &, set<int> &);
int reduce_catalog_snps(map<int, CSLocus *> &, map<int, set<int> > &, PopMap<CSLocus> *);

/*
 * create_catalog(vector<VcfRecord>&):
 * Creates a catalog based on VCF SNP records.
 *
 * We observe the following rules to create the catalog loci :
 * [sample_id] (batch number) Always set to 0.
 * [id] VCF records do not intrinsically have locus ids; we use the SNP records indexes.
 * [len] Always set to 1.
 * [con] We use the reference nucleotide as the consensus.
 * [loc] Use the chromosome and position given by each record (n.b. the
 *       VCF format requires these field), and strand "strand_plus".
 * [snps] Use the ref+alt alleles.
 *     [col] Always set to 0 (first nucleotide in the consensus).
 *     [type] always 'unk', on the premise that it is not used by populations.
 *     [lratio] Always set to 0.
 *     [rank_1] The ref allele.
 *     [rank_2], [rank_3], [rank_4] The alt allele(s).
 * [alleles] Use the ref+alt alleles in the order they appear, skipping
 *     the special '*' allele ('variant is irrelevant in certain context
 *     because of an neighboring structural polymorphism') if present.
 * [strings] We fill this by calling Locus::populate_alleles().
 * [cnt] Set to the approriate value when filling the PopMap.
 * [hcnt] (Same as above.)
 * [confounded_cnt] (Same as above.)
 * [gmap] This is filled by tabulate_haplotypes() (in "populations.cc").
 * [gcnt] (Same as above.)
 * [marker] (Same as above.)
 *
 * When no depth information is available, the [depth] and the depths
 * of the [alleles] are set to 0. (n.b. the parsing of depth information in
 * VCF is not implemented as of Mar 21, 2016.)
 *
 * When no likelihood information is available, [lnl] is set to 0. (n.b.
 * the parsing of likelihood information in VCF is not implemented as of
 * Mar 21, 2016.)
 *
 * The following members are left unset, on the premise that
 * "populations" does not use them :
 * model, blacklisted, deleveraged, lumberjack, components, reads, comp_cnt,
 * comp_type, annotation, uncor_marker, hap_cnts, f, trans_gcnt, chisq.
 *
 * REMOVED 2017-10-20 (~v2.0Beta2b), prototype was:
 * ```
 * map<int, CSLocus*>* create_catalog(const vector<VcfRecord>& vcf_records);
 * ```
 *
 * new_cslocus
 * ===========
 *
 * Creates a CSLocus based on a consensus sequence and a set of VCF records for
 * the locus. The fasta identifier and all the VCF records should have the same
 * locus ID.
 *
 * We observe the following rules to create the locus (@ when the value is obvious) :
 * [sample_id] (batch number) Always set to 0.
 * [id] @
 * [con] the fasta sequence
 * [len] @
 * [loc] for denovo, this should be {"", 0, strand_plus} TODO ref-based.
 * [snps] polymorphic positions taken from the VCF; i.e. VCF records where
 *         ALT is "." are ignored.
 *     [col] @
 *     [type] always "snp_type_het"
 *     [lratio] always 0 (this isn't used by populations).
 *     [rank_1, 2, 3, 4] @
 * [alleles] We fill the haplotype frequency map by reconstituting each sample's
 *         pair of haplotypes from the phased SNP genotypes (in the VCF).
 * [strings] @
 * [depth] is always 0.
 * [lnl] is always 0.
 *
 * Other members need not be set, c.f. `create_catalog(vector<VcfRecord>&)`.
 */
CSLocus *new_cslocus(const Seq& consensus, const vector<VcfRecord>& records, int id);

//
// Create a single catalog locus based on an external VCF file.
//
CSLocus *new_cslocus(const VcfRecord rec, int id);

#endif // __CATALOG_UTILS_H__
