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

#ifndef __POPULATIONS_H__
#define __POPULATIONS_H__

#ifdef _OPENMP
#include <omp.h>    // OpenMP library
#endif
#include <getopt.h> // Process command-line options
#include <dirent.h> // Open/Read contents of a directory
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <utility>
#include <string>
#include <iostream>
#include <fstream>
#include <iomanip>
using std::setw;
using std::setprecision;
using std::fixed;
#include <sstream>
#include <vector>
#include <map>
#include <set>
using std::unique_ptr;

#include "constants.h"
#include "stacks.h"
#include "locus.h"
#include "renz.h"
#include "PopMap.h"
#include "PopSum.h"
#include "utils.h"
#include "log_utils.h"
#include "catalog_utils.h"
#include "sql_utilities.h"
#include "genotype_dictionaries.h"
#include "ordered.h"
#include "smoothing.h"
#include "bootstrap.h"
#include "MetaPopInfo.h"
#include "gzFasta.h"
#include "locus_readers.h"

enum corr_type {p_value, bonferroni_win, bonferroni_gen, no_correction};
enum bs_type   {bs_exact, bs_approx, bs_none};
enum merget    {merge_sink, merge_src};
enum phaset    {merge_failure, simple_merge, complex_phase, nomapping_fail, multimapping_fail, multiple_fails};
enum class InputMode {stacks, stacks2, vcf};

const int max_snp_dist = 500;

struct LocBin {
    CSLocus   *cloc;
    Datum    **d;
    LocPopSum *s;

    LocBin(): cloc(NULL), d(NULL), s(NULL) {}
    ~LocBin() {
        if (this->cloc != NULL) delete cloc;
        if (this->d    != NULL) delete [] d;
        if (this->s    != NULL) delete s;
    }
};

//
// Class for accumulating summary statistics as we read each batch of data.
// After all loci are processed, output the summary statistics summary.
//
class SumStatsSummary {
    size_t  _pop_cnt;
    int    *_private_cnt;
    double *_num_indv_mean,         *_p_mean,         *_obs_het_mean,         *_obs_hom_mean,         *_exp_het_mean,         *_exp_hom_mean,         *_pi_mean,         *_fis_mean;
    double *_num_indv_acc_mean,     *_p_acc_mean,     *_obs_het_acc_mean,     *_obs_hom_acc_mean,     *_exp_het_acc_mean,     *_exp_hom_acc_mean,     *_pi_acc_mean,     *_fis_acc_mean;
    double *_num_indv_var,          *_p_var,          *_obs_het_var,          *_obs_hom_var,          *_exp_het_var,          *_exp_hom_var,          *_pi_var,          *_fis_var;
    double *_num_indv_mean_all,     *_p_mean_all,     *_obs_het_mean_all,     *_obs_hom_mean_all,     *_exp_het_mean_all,     *_exp_hom_mean_all,     *_pi_mean_all,     *_fis_mean_all;
    double *_num_indv_acc_mean_all, *_p_acc_mean_all, *_obs_het_acc_mean_all, *_obs_hom_acc_mean_all, *_exp_het_acc_mean_all, *_exp_hom_acc_mean_all, *_pi_acc_mean_all, *_fis_acc_mean_all;
    double *_num_indv_var_all,      *_p_var_all,      *_obs_het_var_all,      *_obs_hom_var_all,      *_exp_het_var_all,      *_exp_hom_var_all,      *_pi_var_all,      *_fis_var_all;
    double *_n, *_n_all, *_var_sites;
    double *_sq_n, *_sq_n_all;

public:
    SumStatsSummary(size_t pop_cnt);
    ~SumStatsSummary();
    int accumulate(const vector<LocBin *> &);
    int final_calculation();
    int write_results();

private:
    double online_variance(double x, double &acc_mean, double n);
};

//
// Class for storing distributions of data related to catalog loci.
//
class CatalogDists {
    map<size_t, size_t> _pre_valid,  _pre_absent,  _pre_confounded;  // Prior to any filtering. 
    map<size_t, size_t> _post_valid, _post_absent, _post_confounded; // After filtering.
    map<size_t, size_t> _pre_snps_per_loc, _post_snps_per_loc;
    
public:
    int accumulate_pre_filtering(const CSLocus *);
    int accumulate(const vector<LocBin *> &);
    int write_results(ostream &log_fh);
};

//
// Class for filtering whole loci based on the sample and population limits (-r, -p).
//
class LocusFilter {
public:
    LocusFilter() {
        this->_pop_cnt    = 0;
        this->_sample_cnt = 0;
        this->_pop_order  = NULL; // The array order of each population.
        this->_samples    = NULL; // Which population each sample belongs to.
        this->_pop_cnts   = NULL; // For a locus, how many samples are present in each population.
        this->_pop_tot    = NULL; // The total number of samples in each population.
        this->_filtered_loci  = 0;
        this->_total_loci     = 0;
        this->_filtered_sites = 0;
        this->_total_sites    = 0;
    }
    LocusFilter(MetaPopInfo *mpopi) {
        assert(mpopi != NULL);
        this->init(mpopi);
    }
    ~LocusFilter() {
        if (this->_pop_order != NULL)
            delete [] this->_pop_order;
        if (this->_samples != NULL)
            delete [] this->_samples;
        if (this->_pop_cnts != NULL)
            delete [] this->_pop_cnts;
        if (this->_pop_tot != NULL)
            delete [] this->_pop_tot;
    }
    // Remove the copy and move constructors.
    LocusFilter(const LocusFilter &) = delete; 
    LocusFilter(LocusFilter &&) = delete;

    void   init(MetaPopInfo *mpopi);
    bool   filter(MetaPopInfo *mpopi, Datum **d);
    bool   prune_sites(MetaPopInfo *mpopi, CSLocus *cloc, Datum **d, LocPopSum *s, ostream &log_fh);
    void   keep_single_snp(MetaPopInfo *mpopi, CSLocus *cloc, Datum **d, LocPopSum *s, map<int, set<int>> &wl);
    void   keep_random_snp(MetaPopInfo *mpopi, CSLocus *cloc, Datum **d, LocPopSum *s, map<int, set<int>> &wl);

    size_t filtered()       const { return this->_filtered_loci; }
    size_t total()          const { return this->_total_loci; }
    size_t filtered_sites() const { return this->_filtered_sites; }
    size_t total_sites()    const { return this->_total_sites; }
    
private:
    size_t  _pop_cnt;
    size_t  _sample_cnt;
    size_t *_pop_order;
    size_t *_samples;
    size_t *_pop_cnts;
    size_t *_pop_tot;
    size_t  _filtered_loci;
    size_t  _total_loci;
    size_t  _filtered_sites;
    size_t  _total_sites;
    
    void reset();
};

//
// BatchLocusProcessor
// ----------
// Class for processing loci in batches, or per chromosome.
//
class BatchLocusProcessor {
public:
    BatchLocusProcessor():
        _input_mode(InputMode::stacks2), _batch_size(0), _mpopi(NULL),
        _vcf_parser(), _cloc_reader(), _fasta_reader(),
        _vcf_header(NULL), _next_loc(NULL), _cloc_vcf_rec(NULL), _ext_vcf_rec(NULL), _catalog(NULL),
        _blacklist(), _whitelist(), _loc_filter(), _unordered_bp(1) {}
    BatchLocusProcessor(InputMode mode, size_t batch_size, MetaPopInfo *popi):
        _input_mode(mode), _batch_size(batch_size), _mpopi(popi),
        _vcf_parser(), _cloc_reader(), _fasta_reader(),
        _vcf_header(NULL), _next_loc(NULL), _cloc_vcf_rec(NULL), _ext_vcf_rec(NULL), _catalog(NULL),
        _blacklist(), _whitelist(), _loc_filter(), _unordered_bp(1) {}
    BatchLocusProcessor(InputMode mode, size_t batch_size): 
        _input_mode(mode), _batch_size(batch_size), _mpopi(NULL),
        _vcf_parser(), _cloc_reader(), _fasta_reader(),
        _vcf_header(NULL), _next_loc(NULL), _cloc_vcf_rec(NULL), _ext_vcf_rec(NULL), _catalog(NULL),
        _blacklist(), _whitelist(), _loc_filter(), _unordered_bp(1) {}
    ~BatchLocusProcessor() {
        if (this->_ext_vcf_rec != NULL)
            delete this->_ext_vcf_rec;
        if (this->_cloc_vcf_rec != NULL)
            delete this->_cloc_vcf_rec;
        if (this->_catalog != NULL)
            delete this->_catalog;
    };
    
    int            init(int, string, string);
    size_t         next_batch(ostream &);
    int            hapstats(ostream &);
    int            write_distributions(ostream &log_fh) { return this->_dists.write_results(log_fh); }
    
    MetaPopInfo*   pop_info()      { return this->_mpopi; }
    int            pop_info(MetaPopInfo *popi) { this->_mpopi = popi; return 0; }
    VcfParser&     vcf_reader()    { return this->_vcf_parser; }
    GzFasta&       fasta_reader()  { return this->_fasta_reader; }
    VcfCLocReader& cloc_reader()   { return this->_cloc_reader; }
    VcfHeader*     vcf_header()    { return this->_vcf_header; }
    size_t         batch_size()    { return this->_batch_size; }
    size_t         batch_size(size_t bsize) { this->_batch_size = bsize; return bsize; }
    
    const map<int, CSLocus *>& catalog()         { return *this->_catalog; }
    const unordered_map<int, vector<VcfRecord>>& cloc_vcf_records() { return *this->_cloc_vcf_rec; }
    const vector<VcfRecord>&   ext_vcf_records() { return *this->_ext_vcf_rec; }
    const set<int>&            blacklist()       { return this->_blacklist; }
    const map<int, set<int>>&  whitelist()       { return this->_whitelist; }
    const LocusFilter&         filter()          { return this->_loc_filter; }
    const vector<LocBin *>&    loci()            { return this->_loci; }
    const CatalogDists&        dists()           { return this->_dists; }

private:
    InputMode    _input_mode;
    size_t       _batch_size; // Number of loci to process at a time.
    MetaPopInfo *_mpopi;      // Population Map

    // Parsers
    VcfParser     _vcf_parser;
    VcfCLocReader _cloc_reader;
    GzFasta       _fasta_reader;

    // Data stores
    VcfHeader        *_vcf_header;
    vector<LocBin *>  _loci;
    LocBin           *_next_loc;
    CatalogDists      _dists;

    // Data stores for external VCF
    unordered_map<int, vector<VcfRecord>> *_cloc_vcf_rec;
    vector<VcfRecord>                     *_ext_vcf_rec;
    map<int, CSLocus *>                   *_catalog;

    // Controls for which loci are loaded
    set<int>           _blacklist;
    map<int, set<int>> _whitelist;
    LocusFilter        _loc_filter;
    size_t             _unordered_bp;
    
private:
    int    init_external_loci(string, string);
    int    init_stacks_loci(int, string, string);
    size_t next_batch_external_loci(ostream &);
    size_t next_batch_stacks_loci(ostream &);
};

void    help( void );
void    version( void );
int     parse_command_line(int, char**);
void    output_parameters(ostream &);
void    open_log(ofstream &);
int     build_file_list();
int     load_marker_list(string, set<int> &);
int     load_marker_column_list(string, map<int, set<int> > &);
int     apply_locus_constraints(map<int, CSLocus *> &, PopMap<CSLocus> *, ofstream &);
int     prune_polymorphic_sites(map<int, CSLocus *> &, PopMap<CSLocus> *, PopSum<CSLocus> *, map<int, set<int> > &, set<int> &, ofstream &);
int     log_haplotype_cnts(map<int, CSLocus *> &, ofstream &);
bool    order_unordered_loci(map<int, CSLocus *> &);
int     merge_shared_cutsite_loci(map<int, CSLocus *> &, PopMap<CSLocus> *, PopSum<CSLocus> *, map<int, pair<merget, int> > &, ofstream &);
phaset  merge_and_phase_loci(PopMap<CSLocus> *, CSLocus *, CSLocus *, set<int> &, ofstream &);
int     merge_datums(int, int, Datum **, Datum **, set<string> &, int);
int     merge_csloci(CSLocus *, CSLocus *, set<string> &);
int     tabulate_haplotypes(map<int, CSLocus *> &, PopMap<CSLocus> *);
int     tabulate_locus_haplotypes(CSLocus *, Datum **, int);
int     create_genotype_map(CSLocus *, Datum **, int);
int     call_population_genotypes(CSLocus *, Datum **, int);
int     translate_genotypes(map<string, string> &, map<string, map<string, string> > &, map<int, CSLocus *> &, PopMap<CSLocus> *, map<int, string> &, set<int> &); // This function doesn't exist (March 24, 2016)
int     correct_fst_bonferroni_win(vector<PopPair *> &);
int     bootstrap_fst_approximate_dist(vector<double> &, vector<int>  &, double *, int *, map<int, vector<double> > &); // not used (March 23, 2016)
int     kernel_smoothed_popstats(map<int, CSLocus *> &, PopMap<CSLocus> *, PopSum<CSLocus> *, int, ofstream &);
int     bootstrap_popstats_approximate_dist(vector<double> &, vector<double> &, vector<int>  &, double *, int *, int, map<int, vector<double> > &, map<int, vector<double> > &); // not used (March 23, 2016)
double  bootstrap_approximate_pval(int, double, map<int, vector<double> > &);

int      calculate_summary_stats(map<int, CSLocus *> &, PopMap<CSLocus> *, PopSum<CSLocus> *);
int      calculate_haplotype_stats(map<int, CSLocus *> &, PopMap<CSLocus> *, PopSum<CSLocus> *);
int      kernel_smoothed_hapstats(vector<CSLocus *> &, PopSum<CSLocus> *, int, double *);
int      calculate_haplotype_divergence(map<int, CSLocus *> &, PopMap<CSLocus> *, PopSum<CSLocus> *);
int      calculate_haplotype_divergence_pairwise(map<int, CSLocus *> &, PopMap<CSLocus> *, PopSum<CSLocus> *);
bool     fixed_locus(Datum **, vector<int> &);

int      nuc_substitution_identity(map<string, int> &, double **);
int      nuc_substitution_identity_max(map<string, int> &, double **);

HapStat *haplotype_amova(Datum **, LocSum **, vector<int> &);
double   amova_ssd_total(vector<string> &, map<string, int> &, double **);
double   amova_ssd_wp(vector<int> &, map<int, vector<int> > &, map<string, int> &, map<int, vector<string> > &, double **);
double   amova_ssd_ap_wg(vector<int> &, map<int, vector<int> > &, map<string, int> &, map<int, vector<string> > &, double **, double **);
double   amova_ssd_ag(vector<int> &, map<int, vector<int> > &, map<string, int> &, map<int, vector<string> > &, double **, double);

double   haplotype_d_est(Datum **, LocSum **, vector<int> &);
LocStat *haplotype_diversity(int, int, Datum **);

void log_snps_per_loc_distrib(ostream&, map<int, CSLocus*>&);

//int  tally_ref_alleles(LocSum **, int, int, char &, char &); //unused; also commented out in the .cc
//int  load_snp_calls(string,  PopMap<CSLocus> *); //no implementation

//bool compare_pop_map(pair<int, string>, pair<int, string>); //no implementation; the function is in [sql_utilities.h]
bool hap_compare(const pair<string,int>&, const pair<string,int>&);

void vcfcomp_simplify_pmap (map<int, CSLocus*>& catalog, PopMap<CSLocus>* pmap);

#endif // __POPULATIONS_H__
