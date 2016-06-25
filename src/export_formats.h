#ifndef EXPORT_FORMATS_H
#define EXPORT_FORMATS_H

#include <iostream>
#include <utility>
#include <map>

#include "locus.h"
#include "PopMap.h"
#include "PopSum.h"
#include "ordered.h" // for "snp"
#include "populations.h" // for "merget", "uncalled_haplotype()", "count_haplotypes_at_locus()"

using std::ofstream;
using std::pair;
using std::map;

class GenPos {
public:
    uint     id;
    uint     bp;
    uint     snp_index;
    loc_type type;

    GenPos(int id, int snp_index, int bp) {
    this->id        = id;
    this->snp_index = snp_index;
    this->bp        = bp;
    this->type      = snp;
    }
    GenPos(int id, int snp_index, int bp, loc_type type) {
    this->id        = id;
    this->snp_index = snp_index;
    this->bp        = bp;
    this->type      = type;
    }
};

int write_sql(map<int, CSLocus *> &, PopMap<CSLocus> *);
int write_fst_stats(map<int, CSLocus *> &, PopMap<CSLocus> *, PopSum<CSLocus> *, ofstream &);
int write_generic(map<int, CSLocus *> &, PopMap<CSLocus> *, bool);
int write_genomic(map<int, CSLocus *> &, PopMap<CSLocus> *);
int write_fasta(map<int, CSLocus *> &, PopMap<CSLocus> *);
int write_strict_fasta(map<int, CSLocus *> &, PopMap<CSLocus> *);
int write_vcf(map<int, CSLocus *> &, PopMap<CSLocus> *, PopSum<CSLocus> *, map<int, pair<merget, int> > &);
int write_vcf_ordered(map<int, CSLocus *> &, PopMap<CSLocus> *, PopSum<CSLocus> *, map<int, pair<merget, int> > &, ofstream &);
int write_vcf_haplotypes(map<int, CSLocus *> &, PopMap<CSLocus> *, PopSum<CSLocus> *);
int write_genepop(map<int, CSLocus *> &, PopMap<CSLocus> *, PopSum<CSLocus> *);
int write_genepop_ordered(map<int, CSLocus *> &, PopMap<CSLocus> *, PopSum<CSLocus> *, ofstream &);
int write_structure(map<int, CSLocus *> &, PopMap<CSLocus> *, PopSum<CSLocus> *);
int write_structure_ordered(map<int, CSLocus *> &, PopMap<CSLocus> *, PopSum<CSLocus> *, ofstream &);
int write_phase(map<int, CSLocus *> &, PopMap<CSLocus> *, PopSum<CSLocus> *);
int write_fastphase(map<int, CSLocus *> &, PopMap<CSLocus> *, PopSum<CSLocus> *);
int write_beagle(map<int, CSLocus *> &, PopMap<CSLocus> *, PopSum<CSLocus> *);
int write_beagle_phased(map<int, CSLocus *> &, PopMap<CSLocus> *, PopSum<CSLocus> *);
int write_plink(map<int, CSLocus *> &, PopMap<CSLocus> *, PopSum<CSLocus> *);
int write_hzar(map<int, CSLocus *> &, PopMap<CSLocus> *, PopSum<CSLocus> *);
int write_treemix(map<int, CSLocus *> &, PopMap<CSLocus> *, PopSum<CSLocus> *);
int write_phylip(map<int, CSLocus *> &, PopMap<CSLocus> *, PopSum<CSLocus> *);
int write_fullseq_phylip(map<int, CSLocus *> &, PopMap<CSLocus> *, PopSum<CSLocus> *);

int populate_snp_calls(map<int, CSLocus *> &, PopMap<CSLocus> *, map<int, pair<merget, int> > &);
int find_datum_allele_depths(Datum *, int, char, char, int, int &, int &);
int tally_observed_haplotypes(vector<char *> &, int, char &, char &);
int tally_haplotype_freq(CSLocus *, PopMap<CSLocus> *, int &, double &, string &);
int datum_adjust_snp_positions(map<int, pair<merget, int> > &, CSLocus *, Datum *, map<int, SNPRes *> &);

#endif // EXPORT_FORMATS_H
