// -*-mode:c++; c-style:k&r; c-basic-offset:4;-*-
//
// Copyright 2012-2017, Julian Catchen <jcatchen@illinois.edu>
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
#ifndef EXPORT_FORMATS_H
#define EXPORT_FORMATS_H

#include <iostream>
#include <utility>
#include <map>

#include "locus.h"
#include "PopMap.h"
#include "PopSum.h"
#include "ordered.h" // for "snp"
#include "populations.h" // for "merget", "InputMode", "uncalled_haplotype()", "count_haplotypes_at_locus()"

enum class ExportType {markers, sumstats, hapstats, snpdivergence, hapdivergence,
        fasta_loci, fasta_raw, fasta_samples, structure, genepop, ordered_genepop, vcf, ordered_vcf};

class Export {
 protected:
    ExportType _type;
    string     _path;
    ofstream   _fh;

 public:
    Export(ExportType t): _type(t) {};
    virtual ~Export() {};
    virtual int  open(const MetaPopInfo *) = 0;
    virtual int  write_header()    = 0;
    virtual int  write_batch(const vector<LocBin *> &) = 0;
    virtual int  post_processing() = 0;
    virtual void close()           = 0;

    ExportType type() { return this->_type; }
    int transpose(ifstream &ifh, vector<string> &transposed);
};

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

    bool operator<(const GenPos& other) const {return bp < other.bp;}
};

class MarkersExport: public Export {
    //
    // Output a list of heterozygous loci and the associated haplotype frequencies.
    //
    const MetaPopInfo *_mpopi;
    
 public:
    MarkersExport();
    ~MarkersExport() {};
    int  open(const MetaPopInfo *mpopi);
    int  write_header();
    int  write_batch(const vector<LocBin *> &);
    int  post_processing() { return 0; }
    void close() {
        this->_fh.close();
        return;
    }
};

class SumstatsExport: public Export {
    //
    // Output the locus-level summary statistics.
    //
    const MetaPopInfo *_mpopi;
    uint  _pop_cnt;

 public:
    SumstatsExport();
    ~SumstatsExport() {};
    int  open(const MetaPopInfo *mpopi);
    int  write_header();
    int  write_batch(const vector<LocBin *> &);
    int  post_processing() { return 0; }
    void close() {
        this->_fh.close();
        return;
    }
};

class HapstatsExport: public Export {
    //
    // Output the locus-level haplotype statistics.
    //
    const MetaPopInfo *_mpopi;
    uint  _pop_cnt;

 public:
    HapstatsExport();
    ~HapstatsExport() {};
    int  open(const MetaPopInfo *mpopi);
    int  write_header();
    int  write_batch(const vector<LocBin *> &);
    int  post_processing() { return 0; }
    void close() {
        this->_fh.close();
        return;
    }
};

class SnpDivergenceExport: public Export {
    //
    // Output the SNP-level divergence statistics.
    //
    const MetaPopInfo *_mpopi;
    vector<ofstream *> _fhs;
    
 public:
    SnpDivergenceExport();
    ~SnpDivergenceExport() {
        for (uint i = 0; i < this->_fhs.size(); i++)
            delete this->_fhs[i];
    };
    int  open(const MetaPopInfo *mpopi);
    int  write_header();
    int  write_batch(const vector<LocBin *> &) { return 0; }
    int  write_batch_pairwise(const vector<LocBin *> &, const vector<vector<PopPair **>> &);
    int  post_processing() { return 0; }
    void close() {
        for (uint i = 0; i < this->_fhs.size(); i++)
            this->_fhs[i]->close();
        return;
    }
};

class HapDivergenceExport: public Export {
    //
    // Output the SNP-level divergence statistics.
    //
    const MetaPopInfo *_mpopi;
    vector<ofstream *> _fhs;
    ofstream *_metapop_fh;
    
 public:
    HapDivergenceExport();
    ~HapDivergenceExport() {
        for (uint i = 0; i < this->_fhs.size(); i++)
            delete this->_fhs[i];
    };
    int  open(const MetaPopInfo *mpopi);
    int  write_header();
    int  write_batch(const vector<LocBin *> &) { return 0; }
    int  write_batch_pairwise(const vector<LocBin *> &, const vector<vector<HapStat *>> &, const vector<HapStat *> &);
    int  post_processing() { return 0; }
    void close() {
        for (uint i = 0; i < this->_fhs.size(); i++)
            this->_fhs[i]->close();
        return;
    }
};

class GenePopExport: public Export {
    //
    // Output a list of heterozygous loci and the associated haplotype frequencies.
    //
    const MetaPopInfo *_mpopi;
    int      _fd;
    string   _tmp_path;
    ostream *_tmpfh;
    ifstream _intmpfh;

 public:
    GenePopExport();
    ~GenePopExport() {};
    int  open(const MetaPopInfo *mpopi);
    int  write_header();
    int  write_batch(const vector<LocBin *> &);
    int  post_processing();
    void close();

 private:
    int  write_site(const CSLocus *, const LocPopSum *, const Datum **, size_t, size_t);
};

class OrderedGenePopExport: public Export {
    //
    // Output a list of heterozygous loci and the associated haplotype frequencies.
    //
    const MetaPopInfo *_mpopi;
    string   _tmp_path;
    ofstream _tmpfh;

 public:
    OrderedGenePopExport();
    ~OrderedGenePopExport() {};
    int  open(const MetaPopInfo *mpopi);
    int  write_header() { return 0; }
    int  write_batch(const vector<LocBin *> &);
    int  post_processing();
    void close();
};

class FastaLociExport: public Export {
    //
    // Output a list of heterozygous loci and the associated haplotype frequencies.
    //
    const MetaPopInfo *_mpopi;
    
 public:
    FastaLociExport();
    ~FastaLociExport() {};
    int  open(const MetaPopInfo *mpopi);
    int  write_header();
    int  write_batch(const vector<LocBin *> &);
    int  post_processing() { return 0; }
    void close() {
        this->_fh.close();
        return;
    }
};

class FastaRawExport: public Export {
    //
    // Output a list of heterozygous loci and the associated haplotype frequencies.
    //
    const MetaPopInfo *_mpopi;
    
 public:
    FastaRawExport();
    ~FastaRawExport() {};
    int  open(const MetaPopInfo *mpopi);
    int  write_header();
    int  write_batch(const vector<LocBin *> &);
    int  post_processing() { return 0; }
    void close() {
        this->_fh.close();
        return;
    }
};

class FastaSamplesExport: public Export {
    //
    // Output a list of heterozygous loci and the associated haplotype frequencies.
    //
    const MetaPopInfo *_mpopi;
    
 public:
    FastaSamplesExport();
    ~FastaSamplesExport() {};
    int  open(const MetaPopInfo *mpopi);
    int  write_header();
    int  write_batch(const vector<LocBin *> &);
    int  post_processing() { return 0; }
    void close() {
        this->_fh.close();
        return;
    }
};

/*
int write_generic(map<int, CSLocus *> &, PopMap<CSLocus> *, bool);
int write_genomic(map<int, CSLocus *> &, PopMap<CSLocus> *);
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
*/

int find_datum_allele_depths(Datum *, int, char, char, int &, int &);
int tally_observed_haplotypes(const vector<char *> &, int, char &, char &);
int tally_haplotype_freq(CSLocus *, Datum **, uint, int &, double &, string &);

#endif // EXPORT_FORMATS_H
