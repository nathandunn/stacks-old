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

//
// populations -- generate population genetic statistics and output
// haplotypes in a population context.
//
#include <cctype>

#include "export_formats.h"
#include "populations.h"

using namespace std;

// Global variables to hold command-line options.
InputMode input_mode  = InputMode::stacks2;
int       num_threads =  1;
string    in_path;
string    in_vcf_path;
string    out_path;
string    pmap_path;
string    bl_file;
string    wl_file;
string    bs_wl_file;
string    enz;
double    sigma             = 150000.0;
double    sample_limit      = 0.0;
int       population_limit  = 1;
int       batch_size        = 10000;
bool      calc_fstats       = false;
bool      bootstrap         = false;
bool      bootstrap_fst     = false;
bool      bootstrap_pifis   = false;
bool      bootstrap_phist   = false;
bool      bootstrap_div     = false;
bs_type   bootstrap_type    = bs_exact;
int       bootstrap_reps    = 100;
bool      bootstrap_wl      = false;
bool      write_single_snp  = false;
bool      write_random_snp  = false;
bool      merge_sites       = false;
bool      expand_id         = false;
bool      write_gtypes      = false;
bool      vcf_haplo_out     = false;
bool      fasta_loci_out    = false;
bool      fasta_samples_out = false;
bool      fasta_samples_raw_out = false;
bool      genomic_out       = false;
bool      structure_out     = false;
bool      phase_out         = false;
bool      fastphase_out     = false;
bool      beagle_out        = false;
bool      beagle_phased_out = false;
bool      plink_out         = false;
bool      hzar_out          = false;
bool      treemix_out       = false;
bool      phylip_out        = false;
bool      phylip_var        = false;
bool      phylip_var_all    = false;
bool      ordered_export    = false;
bool      smooth_fstats     = false;
bool      smooth_popstats   = false;
bool      loci_ordered      = false;
bool      log_fst_comp      = false;
bool      verbose           = false;
bool      filter_lnl        = false;
double    lnl_limit         = 0.0;
int       min_stack_depth   = 0;
double    merge_prune_lim   = 1.0;
double    minor_allele_freq = 0.0;
double    max_obs_het       = 1.0;
double    p_value_cutoff    = 0.05;
corr_type fst_correction    = no_correction;
set<string> debug_flags;

string           out_prefix;
MetaPopInfo      mpopi;
vector<Export *> exports;
set<int>         bootstraplist;

//
// Hold information about restriction enzymes
//
map<string, const char **> renz;
map<string, int>           renz_cnt;
map<string, int>           renz_len;
map<string, int>           renz_olap;

int main (int argc, char* argv[]) {
    IF_NDEBUG_TRY

#ifndef HAVE_LIBZ
    cerr << "Stacks was compiled without zlib, and will refuse to parse compressed files.\n";
#endif

    //
    // Initialize the globals that need it.
    //
    initialize_renz(renz, renz_cnt, renz_len);
    initialize_renz_olap(renz_olap);
    srandom(time(NULL));

    //
    // Parse the command line.
    //
    parse_command_line(argc, argv);

    //
    // Open and initialize the log file.
    //
    LogAlterator *logger = new LogAlterator(out_path + out_prefix + ".log", false, argc, argv);

    output_parameters(cerr);

    //
    // Set the number of OpenMP parallel threads to execute.
    //
    #ifdef _OPENMP
    omp_set_num_threads(num_threads);
    #endif

    //
    // Read the population map file, if any.
    //
    if (not pmap_path.empty()) {
        cerr << "Parsing population map...\n";
        mpopi.init_popmap(pmap_path);
        cerr << "The population map contained " << mpopi.samples().size() << " samples, "
             << mpopi.pops().size() << " population(s), " << mpopi.groups().size() << " group(s).\n";
    }

    //
    // Locate and open input files, read VCF headers, parse population map, load black/white lists.
    //
    BatchLocusProcessor bloc(input_mode, batch_size, &mpopi);

    bloc.init(in_path, pmap_path);

    //
    // Report information on the structure of the populations specified.
    //
    mpopi.status();

    if (size_t(population_limit) > mpopi.pops().size()) {
        cerr << "Notice: Population limit (" << population_limit << ")"
             << " larger than number of popualtions present, adjusting parameter to "
             << mpopi.pops().size() << "\n";
        population_limit = mpopi.pops().size();
    }

    //
    // Setup the default data exports.
    //
    exports.push_back(new MarkersExport());
    exports.push_back(new GenotypesExport());
    exports.push_back(new SumstatsExport());
    exports.push_back(new HapstatsExport());

    SnpDivergenceExport *sdiv_exp = NULL;
    HapDivergenceExport *hdiv_exp = NULL;
    if (calc_fstats) {
        sdiv_exp = new SnpDivergenceExport();
        exports.push_back(sdiv_exp);
        hdiv_exp = new HapDivergenceExport();
        exports.push_back(hdiv_exp);
    }

    //
    // Setup the kernel smoothing apparatus.
    //
    LocusSmoothing *smooth = NULL;
    if (smooth_fstats || smooth_popstats)
        smooth = new LocusSmoothing(&mpopi, logger->l);

    //
    // Setup the divergence statistics calculator, if requested.
    //
    LocusDivergence *ldiv = NULL;
    if (calc_fstats)
        ldiv = new LocusDivergence(&mpopi);

    //
    // Open the export files and write any headers.
    //
    for (uint i = 0; i < exports.size(); i++) {
        exports[i]->open((const MetaPopInfo *) &mpopi);
        exports[i]->write_header();
    }

    //
    // Initialize the summary statistics object which will accumulate the metapopulation summary statistics
    // as the individual loci are read and processed.
    //
    SumStatsSummary sumstats(mpopi.pops().size());

    const LocusFilter &filter = bloc.filter();
    cerr << "\nProcessing data in batches;\n";
    int loc_cnt = 0;

    do {
        //
        // Read the next set of loci to process.
        // - If data are denovo, load blim._batch_size loci.
        // - If data are reference aligned, load one chromosome.
        // - Filter the loci according to command line parameters (-r, -p, --maf, --write_single_snp, etc.)
        // - Sort the loci by basepair if they are ordered.
        //
        cerr << "  Begin batch " << bloc.next_batch_number() << "...";
        loc_cnt  = bloc.next_batch(logger->l);

        cerr << "analyzed " << filter.batch_total() << " loci";
        if (loci_ordered && bloc.loci().size() > 0)
            cerr <<  " from " << bloc.loci().front()->cloc->loc.chr();
        cerr << "; filtered " << filter.batch_filtered() << " loci; " << filter.batch_seen() << " loci seen.\n";

        if (loc_cnt == 0 && filter.batch_seen() == 0) break;

        sumstats.accumulate(bloc.loci());

        //
        // Calculate haplotype and gene diversity per locus per population.
        //
        bloc.hapstats(logger->l);

        //
        // Export this subset of the loci.
        //
        for (uint i = 0; i < exports.size(); i++)
            exports[i]->write_batch(bloc.loci());

        //
        // Calculate and report the extent of overlap between different RAD loci.
        //
        if (loci_ordered)
            cerr << "    " << bloc.report_locus_overlap( (verbose ? &logger->l : NULL) ) << "\n";

        //
        // Calculate divergence statistics (Fst), if requested.
        //
        if (calc_fstats) {
            cerr << "    Calculating F statistics...";
            ldiv->snp_divergence(bloc.loci());
            ldiv->haplotype_divergence_pairwise(bloc.loci());
            ldiv->haplotype_divergence(bloc.loci());
            cerr << "done.\n";
        }

        //
        // Smooth population statistics across individual populations, and between populations.
        //
        if ( (smooth_fstats || smooth_popstats) && loci_ordered == false) {
            cerr << "    Notice: Smoothing was requested (-k), but will not be performed as the loci are not ordered.\n";
        } else if (smooth_fstats || smooth_popstats) {
            cerr << "    Generating kernel-smoothed population statistics...";
            if (smooth_popstats) {
                smooth->snpstats(bloc.loci(), logger->l);
                smooth->hapstats(bloc.loci(), logger->l);
            }
            if (smooth_fstats) {
                smooth->snp_divergence(bloc.loci(), ldiv->snp_values(), logger->l);
                smooth->hap_divergence(bloc.loci(), ldiv->haplotype_values(), ldiv->metapop_haplotype_values(), logger->l);
            }
            cerr << "done.\n";
        }

        if (calc_fstats) {
            sdiv_exp->write_batch_pairwise(bloc.loci(), ldiv->snp_values());
            hdiv_exp->write_batch_pairwise(bloc.loci(), ldiv->haplotype_values(), ldiv->metapop_haplotype_values());
            ldiv->clear(bloc.loci());
        }

    } while (loc_cnt > 0 || filter.batch_seen() > 0);

    //
    // Report what we read from the input files.
    //
    bloc.summarize(cerr);

    cerr << "\n"
         << "Removed " << filter.filtered() << " loci that did not pass sample/population constraints from " << filter.seen() << " loci.\n"
         << "Kept " << filter.total() << " loci, composed of " << filter.total_sites() << " sites; "
         << filter.filtered_sites() << " of those sites were filtered, " << filter.variant_sites() << " variant sites remained.\n";

    //
    // Do the final sumstats calculations and write the sumstats summary files.
    //
    sumstats.final_calculation();
    sumstats.write_results();

    if (calc_fstats)
        ldiv->write_summary(out_path + out_prefix);

    //
    // Write out the distributions of catalog loci.
    //
    bloc.write_distributions(logger->l);
    bloc.cleanup();

    if (smooth_fstats || smooth_popstats)
        delete smooth;
    if (calc_fstats)
        delete ldiv;

    //
    // Close the export files and do any required post processing.
    //
    for (uint i = 0; i < exports.size(); i++) {
        exports[i]->post_processing();
        exports[i]->close();
        delete exports[i];
    }

    delete logger;

    cerr << "Populations is done.\n";

    return 0;

    IF_NDEBUG_CATCH_ALL_EXCEPTIONS
}

int
BatchLocusProcessor::init(string in_path, string pmap_path)
{
    //
    // Read the blacklist and whitelist to control which loci we load.
    //
    int cnt;
    if (bl_file.length() > 0) {
        cnt = this->_loc_filter.load_blacklist(bl_file);
        cerr << "Loaded " << cnt << " blacklisted markers.\n";
    }
    if (wl_file.length() > 0) {
        cnt = this->_loc_filter.load_whitelist(wl_file);
        cerr << "Loaded " << cnt << " whitelisted markers.\n";
        this->_user_supplied_whitelist = true;
    }

    if (this->_input_mode == InputMode::vcf)
        this->init_external_loci(in_vcf_path, pmap_path);
    else
        this->init_stacks_loci(in_path, pmap_path);

    return 0;
}

size_t
BatchLocusProcessor::next_batch(ostream &log_fh)
{
    size_t loc_cnt;

    //
    // Clear out any loci from the previous batch.
    //
    for (uint i = 0; i < this->_loci.size(); i++)
        delete this->_loci[i];
    this->_loci.clear();

    this->_batch_num++;

    if (this->_input_mode == InputMode::vcf)
        loc_cnt = this->next_batch_external_loci(log_fh);
    else
        loc_cnt = this->next_batch_stacks_loci(log_fh);

    return loc_cnt;
}

int
BatchLocusProcessor::init_stacks_loci(string in_path, string pmap_path)
{
    //
    // Open the files.
    //
    string catalog_fa_path  = in_path + "gstacks.fa.gz";
    string catalog_vcf_path = in_path + "gstacks.vcf.gz";

    this->_fasta_reader.open(catalog_fa_path);
    this->_cloc_reader.open(catalog_vcf_path);

    // Create the population map or check that all samples have data.
    if (pmap_path.empty()) {
        cerr << "No population map specified, using all samples...\n";
        this->_mpopi->init_names(this->cloc_reader().header().samples());
    } else {
        size_t n_samples_before = this->_mpopi->samples().size();
        this->_mpopi->intersect_with(this->cloc_reader().header().samples());
        size_t n_rm_samples = n_samples_before - this->_mpopi->samples().size();

        if (n_rm_samples > 0) {
            cerr << "Warning: No genotype data exists for " << n_rm_samples
                 << " of the samples listed in the population map.\n";
            if (this->_mpopi->samples().empty()) {
                cerr << "Error: No more samples.\n";
                throw exception();
            }
        }
    }

    this->cloc_reader().set_sample_ids(*this->_mpopi);

    //
    // Initialize the locus filter after we have constructed the population map.
    //
    this->_loc_filter.init(this->pop_info());

    //
    // Store the VCF header data.
    //
    this->_vcf_header = new VcfHeader(this->cloc_reader().header());

    return 0;
}

size_t
BatchLocusProcessor::next_batch_stacks_loci(ostream &log_fh)
{
    vector<VcfRecord> records;
    Seq seq;
    seq.id      = new char[id_len];
    seq.comment = new char[id_len];
    seq.seq     = new char[max_len];

    LocBin *loc;
    int     cloc_id, rv;
    string  prev_chr, cur_chr;
    size_t  loc_cnt = 0;

    this->_loci.clear();
    this->_loc_filter.batch_clear();

    //
    // Check if we queued a LocBin object from the last round of reading.
    //
    if (this->_next_loc != NULL) {
        this->_loci.push_back(this->_next_loc);
        prev_chr        = this->_next_loc->cloc->loc.chr();
        this->_loc_filter.locus_seen();
        this->_loc_filter.keep_locus(this->_next_loc);
        this->_next_loc = NULL;
        loc_cnt++;
    }

    do {
        if (!this->_cloc_reader.read_one_locus(records)) break;

        //
        // Get the current locus ID.
        //
        assert(!records.empty());
        cloc_id = is_integer(records[0].chrom());
        assert(cloc_id >= 0);

        //
        // Find the corresponding fasta record. (Note: c-loci with very low coverage
        // might be entirely missing from the VCF; in this case ignore them.)
        //
        do {
            rv = this->_fasta_reader.next_seq(seq);
        } while (rv != 0 && is_integer(seq.id) != cloc_id);

        if (rv == 0) {
            cerr << "Error: catalog VCF and FASTA files are discordant, maybe trucated. rv: " << rv << "; cloc_id: " << cloc_id << "\n";
            throw exception();
        }

        //
        // Create and populate a new catalog locus.
        //
        loc = new LocBin(this->_mpopi->samples().size());
        loc->cloc = new_cslocus(seq, records, cloc_id);

        //
        // Set the current chromosome.
        //
        cur_chr = loc->cloc->loc.chr();
        if (prev_chr.length() == 0)
            prev_chr = cur_chr;

        //
        // If we are in ordered mode, and there were no loci analyzed on the previous chromosome,
        // keep processing the current chromosome without returning.
        //
        if (cur_chr != prev_chr && this->_loc_filter.batch_filtered() == this->_loc_filter.batch_seen())
            prev_chr = cur_chr;

        this->_loc_filter.locus_seen();

        //
        // Check if this locus is filtered.
        //
        if (this->_user_supplied_whitelist && this->_loc_filter.whitelist_filter(cloc_id)) {
            records.clear();
            delete loc;
            continue;
        }
        if (this->_loc_filter.blacklist_filter(cloc_id)) {
            records.clear();
            delete loc;
            continue;
        }

        //
        // Create and populate a map of the population genotypes.
        //
        loc->d = new Datum *[loc->sample_cnt];
        for (size_t i = 0; i < this->_mpopi->samples().size(); i++) loc->d[i] = NULL;
        PopMap<CSLocus>::populate_locus(loc->d, *loc->cloc, records, *this->_vcf_header, *this->_mpopi);
        records.clear();

        this->_dists.accumulate_pre_filtering(loc->cloc);

        //
        // Apply locus constraints to remove entire loci below the -r/-p thresholds.
        //
        if (this->_loc_filter.filter(this->_mpopi, loc->d)) {
            delete loc;
            continue;
        }

        //
        // Create the PopSum object and compute the summary statistics for this locus.
        //
        loc->s = new LocPopSum(strlen(loc->cloc->con), *this->_mpopi);
        loc->s->sum_pops(loc->cloc, (const Datum **) loc->d, *this->_mpopi, verbose, cerr);
        loc->s->tally_metapop(loc->cloc);

        //
        // If write_single_snp or write_random_snp has been specified, mark sites to be pruned using the whitelist.
        //
        if (write_single_snp)
            this->_loc_filter.keep_single_snp(loc->cloc, loc->s->meta_pop());
        else if (write_random_snp)
            this->_loc_filter.keep_random_snp(loc->cloc, loc->s->meta_pop());
        //
        // Prune the sites according to the whitelist.
        //
        this->_loc_filter.prune_sites_with_whitelist(this->_mpopi, loc->cloc, loc->d, this->_user_supplied_whitelist);

        //
        // Identify individual SNPs that are below the -r threshold or the minor allele
        // frequency threshold (-a). In these cases we will remove the SNP, but keep the locus.
        // If all SNPs are filtered, delete the locus.
        //
        if (this->_loc_filter.prune_sites_with_filters(this->_mpopi, loc->cloc, loc->d, loc->s, log_fh)) {
            delete loc;
            continue;
        }

        //
        // If these data are unordered, provide an arbitrary ordering.
        //
        if (loc->cloc->loc.empty()) {
            loc->cloc->loc.set("un", this->_unordered_bp, strand_plus);
            this->_unordered_bp += loc->cloc->len;
        } else {
            loci_ordered = true;
        }

        //
        // Regenerate summary statistics after pruning SNPs.
        //
        loc->s->sum_pops(loc->cloc, (const Datum **) loc->d, (const MetaPopInfo &) *this->_mpopi, verbose, cerr);
        loc->s->tally_metapop(loc->cloc);

        //
        // Tabulate haplotypes present and in what combinations.
        //
        tabulate_locus_haplotypes(loc->cloc, loc->d, this->_mpopi->samples().size());

        if (cur_chr == prev_chr) {
            this->_loci.push_back(loc);
            loc_cnt++;

            this->_loc_filter.keep_locus(loc);

        } else {
            this->_next_loc = loc;
            this->_loc_filter.locus_unsee();
        }

    } while (( loci_ordered && prev_chr == cur_chr) ||
             (!loci_ordered && loc_cnt   < this->_batch_size));

    //
    // Record the post-filtering distribution of catalog loci for this batch.
    //
    this->_dists.accumulate(this->_loci);

    //
    // Sort the catalog loci, if possible.
    //
    if (loci_ordered)
        sort(this->_loci.begin(), this->_loci.end(),
             [] (const LocBin *a, const LocBin *b) -> bool {
                 return a->cloc->loc.bp < b->cloc->loc.bp;
             });

    return loc_cnt;
}

int
BatchLocusProcessor::init_external_loci(string in_path, string pmap_path)
{
    //
    // Open the VCF file
    //
    cerr << "Opening the VCF file...\n";
    this->_vcf_parser.open(in_path);

    if (this->_vcf_parser.header().samples().empty()) {
        cerr << "Error: No samples in VCF file '" << in_path << "'.\n";
        throw exception();
    }

    // Reconsider the MetaPopInfo in light of the VCF header.
    if (pmap_path.empty()) {
        cerr << "No population map specified, creating one from the VCF header...\n";
        this->_mpopi->init_names(this->_vcf_parser.header().samples());

    } else {
        // Intersect the samples present in the population map and the VCF.
        size_t n_samples_before = this->_mpopi->samples().size();

        this->_mpopi->intersect_with(this->_vcf_parser.header().samples());

        size_t n_rm_samples = n_samples_before - this->_mpopi->samples().size();
        if (n_rm_samples > 0) {
            cerr << "Warning: Of the samples listed in the population map, "
                 << n_rm_samples << " could not be found in the VCF :";
            if (this->_mpopi->samples().empty()) {
                cerr << "Error: No more samples.\n";
                throw exception();
            }
        }
    }

    // Create arbitrary sample IDs.
    for (size_t i = 0; i < this->_mpopi->samples().size(); ++i)
        this->_mpopi->set_sample_id(i, i+1); //id=i+1

    //
    // Initialize the locus filter after we have constructed the population map.
    //
    this->_loc_filter.init(this->pop_info());

    //
    // Store the VCF header data.
    //
    this->_vcf_header = new VcfHeader(this->vcf_reader().header());

    this->_total_ext_vcf = 0;

    return 0;
}

size_t
BatchLocusProcessor::next_batch_external_loci(ostream &log_fh)
{
    //
    // VCF mode
    //
    LocBin   *loc;
    string    prev_chr, cur_chr;
    size_t    loc_cnt = 0;
    size_t    cloc_id = 1;
    VcfRecord rec;

    //
    // Check if we queued a LocBin object from the last round of reading.
    //
    if (this->_next_loc != NULL) {
        this->_loci.push_back(this->_next_loc);
        cloc_id         = this->_next_loc->cloc->id + 1;
        prev_chr        = this->_next_loc->cloc->loc.chr();
        this->_next_loc = NULL;
        loc_cnt++;
    }

    do {
        if (!this->_vcf_parser.next_record(rec)) break;

        this->_total_ext_vcf++;

        // Check for a SNP.
        if (not rec.is_snp()) {
            this->_skipped_notsnp.push_back(this->_vcf_parser.line_number());
            continue;
        }

        // Check for a filtered-out SNP
        if (strncmp(rec.filters(), ".", 2) != 0 && strncmp(rec.filters(), "PASS", 5) != 0) {
            this->_skipped_filter.push_back(this->_vcf_parser.line_number());
            continue;
        }

        //
        // Create and populate a new catalog locus.
        //
        loc = new LocBin(this->_mpopi->samples().size());
        loc->cloc = new_cslocus(rec, cloc_id);

        if (loc->cloc == NULL) {
            delete loc;
            continue;
        }

        //
        // Create and populate a map of the population genotypes.
        //
        loc->d = new Datum *[loc->sample_cnt];
        for (size_t i = 0; i < this->_mpopi->samples().size(); i++) loc->d[i] = NULL;
        PopMap<CSLocus>::populate_locus(loc->d, *loc->cloc, rec, (const VcfHeader &) *this->_vcf_header, (const MetaPopInfo &) *this->_mpopi);

        this->_dists.accumulate_pre_filtering(loc->cloc);

        //
        // Apply locus constraints to remove entire loci below the -r/-p thresholds.
        //
        if (this->_loc_filter.filter(this->_mpopi, loc->d)) {
            delete loc;
            continue;
        }

        //
        // Create the PopSum object and compute the summary statistics for this locus.
        //
        loc->s = new LocPopSum(strlen(loc->cloc->con), (const MetaPopInfo &) *this->_mpopi);
        loc->s->sum_pops(loc->cloc, (const Datum **) loc->d, (const MetaPopInfo &) *this->_mpopi, verbose, cerr);
        loc->s->tally_metapop(loc->cloc);

        //
        // If write_single_snp or write_random_snp has been specified, mark sites to be pruned using the whitelist.
        //
        if (write_single_snp)
            this->_loc_filter.keep_single_snp(loc->cloc, loc->s->meta_pop());
        else if (write_random_snp)
            this->_loc_filter.keep_random_snp(loc->cloc, loc->s->meta_pop());
        //
        // Prune the sites according to the whitelist.
        //
        this->_loc_filter.prune_sites_with_whitelist(this->_mpopi, loc->cloc, loc->d, this->_user_supplied_whitelist);

        //
        // Identify individual SNPs that are below the -r threshold or the minor allele
        // frequency threshold (-a). In these cases we will remove the SNP, but keep the locus.
        // If all SNPs are filtered, delete the locus.
        //
        if (this->_loc_filter.prune_sites_with_filters(this->_mpopi, loc->cloc, loc->d, loc->s, log_fh)) {
            delete loc;
            continue;
        }

        //
        // Regenerate summary statistics after pruning SNPs.
        //
        loc->s->sum_pops(loc->cloc, (const Datum **) loc->d, (const MetaPopInfo &) *this->_mpopi, verbose, cerr);
        loc->s->tally_metapop(loc->cloc);

        //
        // If these data are unordered, provide an arbitrary ordering.
        //
        if (loc->cloc->loc.empty()) {
            loc->cloc->loc.set("un", this->_unordered_bp, strand_plus);
            this->_unordered_bp += strlen(loc->cloc->con);
        } else {
            loci_ordered = true;
        }

        cur_chr = loc->cloc->loc.chr();
        if (prev_chr.length() == 0)
            prev_chr = cur_chr;

        //
        // Tabulate haplotypes present and in what combinations.
        //
        tabulate_locus_haplotypes(loc->cloc, loc->d, this->_mpopi->samples().size());

        if (cur_chr == prev_chr) {
            this->_loci.push_back(loc);
            loc_cnt++;
        } else {
            this->_next_loc = loc;
        }
        cloc_id++;

    } while (( loci_ordered && prev_chr == cur_chr) ||
             (!loci_ordered && loc_cnt   < this->_batch_size));

    //
    // Record the post-filtering distribution of catalog loci for this batch.
    //
    this->_dists.accumulate(this->_loci);

    //
    // Sort the catalog loci, if possible.
    //
    if (loci_ordered)
        sort(this->_loci.begin(), this->_loci.end(),
             [] (const LocBin *a, const LocBin *b) -> bool {
                 return a->cloc->loc.bp < b->cloc->loc.bp;
             });

    return loc_cnt;
}

int
BatchLocusProcessor::summarize(ostream &log_fh)
{
    if (this->_input_mode == InputMode::vcf) {
        log_fh << "Found " << this->_total_ext_vcf << " SNP records in file '" << in_vcf_path
             << "'. (Skipped " << this->_skipped_filter.size() << " already filtered-out SNPs and "
             << this->_skipped_notsnp.size() << " non-SNP records ; more with --verbose.)\n";
        if (verbose && not this->_skipped_notsnp.empty()) {
            log_fh << "The following VCF record lines were determined not to be SNPs and skipped :";
            for (vector<size_t>::const_iterator l = this->_skipped_notsnp.begin(); l != this->_skipped_notsnp.end(); ++l)
                log_fh << " " << *l;
            log_fh << "\n";
        }
    } else {
    }

    return 0;
}

int
BatchLocusProcessor::hapstats(ostream &log_fh)
{
    for (uint i = 0; i < this->_loci.size(); i++) {
        LocBin *loc = this->_loci[i];

        loc->s->calc_hapstats(loc->cloc, (const Datum **) loc->d, *this->_mpopi);
    }

    return 0;
}

int
BatchLocusProcessor::cleanup()
{
    if (this->_input_mode == InputMode::vcf) {
        delete this->_vcf_header;

    } else {
        delete this->_vcf_header;
    }

    return 0;
}

string
BatchLocusProcessor::report_locus_overlap(ofstream *log_fh)
{
    const CSLocus  *loc;
    const LocSum   *lsum;

    //
    // Create a vector of maps, to record for each population, which sites overlap between different loci.
    //
    vector<map<uint, set<uint>>> sites;
    for (uint i = 0; i < this->_mpopi->pops().size(); i++) {
        sites.push_back(map<uint, set<uint>>());
        map<uint, set<uint>> &pop_sites = sites.back();

        for (uint pos = 0; pos < this->_loci.size(); pos++) {
            loc  = this->_loci[pos]->cloc;
            lsum = this->_loci[pos]->s->per_pop(i);

            for (uint k = 0; k < loc->len; k++) {
                if (lsum->nucs[k].num_indv == 0)
                    continue;

                if (pop_sites.count(lsum->nucs[k].bp) == 0)
                    pop_sites[lsum->nucs[k].bp] = set<uint>();
                pop_sites[lsum->nucs[k].bp].insert(loc->id);
            }
        }
    }

    map<uint, set<uint>> final_map;

    size_t multiple_sites = 0;
    //
    // Iterate over each population and combine all sites that overlap in one or more populations.
    //
    for (uint i = 0; i < sites.size(); i++) {
        for (map<uint, set<uint>>::iterator it = sites[i].begin(); it != sites[i].end(); it++) {

            if (final_map.count(it->first) == 0)
                final_map[it->first] = set<uint>();

            for (set<uint>::iterator j = it->second.begin(); j != it->second.end(); j++)
                final_map[it->first].insert(*j);
        }
    }

    size_t uniq_sites = final_map.size();

    if (verbose && log_fh != NULL)
        *log_fh << "The following loci overlap on genomic segment " << this->_loci[0]->cloc->loc.chr() << ":\n";

    for (map<uint, set<uint>>::iterator it = final_map.begin(); it != final_map.end(); it++) {
        if (it->second.size() > 1) {
            multiple_sites++;

            if (verbose && log_fh != NULL) {
                for (set<uint>::iterator j = it->second.begin(); j != it->second.end(); j++)
                    *log_fh << "    " << *j << ", ";
                *log_fh << it->first + 1  << "bp\t" << "\n";
            }
        }
    }

    char   buf[max_len];
    double pct;

    pct = (double) multiple_sites / (double) uniq_sites * 100;

    snprintf(buf, max_len,
             "%ld unique genomic sites, %ld sites covered by multiple loci (%.03f%%)",
             uniq_sites, multiple_sites, pct);

    return string(buf);
}


void
LocusFilter::batch_clear()
{
    this->_batch_total_loci    = 0;
    this->_batch_filtered_loci = 0;
    this->_batch_seen_loci     = 0;
}

void
LocusFilter::locus_seen()
{
    this->_seen_loci++;
    this->_batch_seen_loci++;
}

void
LocusFilter::locus_unsee()
{
    this->_seen_loci--;
    this->_batch_seen_loci--;
}

void
LocusFilter::keep_locus(LocBin *loc)
{
    this->_total_loci++;
    this->_batch_total_loci++;

    //
    // Count up the number of sites and variable sites.
    //
    const LocTally *t = loc->s->meta_pop();

    for (uint i = 0; i < loc->cloc->len; i++) {
        if (t->nucs[i].fixed == false)
            this->_variant_sites++;
        this->_total_sites++;
    }
}

bool
LocusFilter::whitelist_filter(size_t locus_id)
{
    if (this->_whitelist.count(locus_id))
        return false;

    this->_filtered_loci++;
    this->_batch_filtered_loci++;

    return true;
}

bool
LocusFilter::blacklist_filter(size_t locus_id)
{
    if (this->_blacklist.count(locus_id)) {
        this->_filtered_loci++;
        this->_batch_filtered_loci++;
        return true;
    }

    return false;
}


bool
LocusFilter::filter(MetaPopInfo *mpopi, Datum **d)
{
    this->reset();

    for (size_t i = 0; i < this->_sample_cnt; i++) {
        //
        // Check that each sample is over the log likelihood threshold.
        //
        if (d[i] != NULL &&
            filter_lnl   &&
            d[i]->lnl < lnl_limit) {
            // below_lnl_thresh++;
            delete d[i];
            d[i] = NULL;
            // loc->hcnt--;
        }
    }

    //
    // Tally up the count of samples in this population.
    //
    for (size_t i = 0; i < this->_sample_cnt; i++) {
        if (d[i] != NULL)
            this->_pop_cnts[this->_samples[i]]++;
    }

    //
    // Check that the counts for each population are over sample_limit. If not, zero out
    // the members of that population.
    //
    double pct = 0.0;

    for (uint i = 0; i < this->_pop_cnt; i++) {
        const Pop& pop = mpopi->pops()[this->_pop_order[i]];

        pct = (double) this->_pop_cnts[i] / (double) this->_pop_tot[i];

        if (this->_pop_cnts[i] > 0 && pct < sample_limit) {
            for (uint j = pop.first_sample; j <= pop.last_sample; j++) {
                if (d[j] != NULL) {
                    delete d[j];
                    d[j] = NULL;
                    // loc->hcnt--;
                }
            }
            this->_pop_cnts[i] = 0;
        }
    }

    //
    // Check that this locus is present in enough populations.
    //
    bool pop_limit = false;
    int  pops      = 0;

    for (uint i = 0; i < this->_pop_cnt; i++)
        if (this->_pop_cnts[i] > 0) pops++;
    if (pops < population_limit)
        pop_limit = true;

    if (pop_limit) {
        this->_filtered_loci++;
        this->_batch_filtered_loci++;
    }

    return pop_limit;
}

void
LocusFilter::init(MetaPopInfo *mpopi)
{
    this->_pop_cnt    = mpopi->pops().size();
    this->_sample_cnt = mpopi->samples().size();

    assert(this->_pop_cnt > 0);
    assert(this->_sample_cnt > 0);

    if (this->_pop_order != NULL)
        delete [] this->_pop_order;
    if (this->_samples != NULL)
        delete [] this->_samples;
    if (this->_pop_cnts != NULL)
        delete [] this->_pop_cnts;
    if (this->_pop_tot != NULL)
        delete [] this->_pop_tot;

    this->_pop_order  = new size_t [this->_pop_cnt];
    this->_samples    = new size_t [this->_sample_cnt];
    this->_pop_cnts   = new size_t [this->_pop_cnt];
    this->_pop_tot    = new size_t [this->_pop_cnt];

    this->_filtered_loci  = 0;
    this->_total_loci     = 0;
    this->_filtered_sites = 0;
    this->_total_sites    = 0;

    size_t pop_sthg = 0;

    for (size_t i_pop = 0; i_pop < mpopi->pops().size(); ++i_pop) {
        const Pop& pop = mpopi->pops()[i_pop];
        this->_pop_tot[pop_sthg]  = 0;

        for (uint i = pop.first_sample; i <= pop.last_sample; i++) {
            this->_samples[i] = pop_sthg;
            this->_pop_tot[pop_sthg]++;
        }
        this->_pop_order[pop_sthg] = i_pop;
        pop_sthg++;
    }
}

void
LocusFilter::reset()
{
    memset(this->_pop_cnts, 0, this->_pop_cnt * sizeof(size_t));
}

int
LocusFilter::keep_single_snp(const CSLocus *cloc, const LocTally *t)
{
    set<int> new_wl;

    //
    // If this locus is not in the whitelist, or it is in the whitelist but no specific
    // SNPs are specified in the whitelist for this locus -- so all SNPs are included,
    // choose the first variant.
    //
    if (this->_whitelist.count(cloc->id) == 0 || this->_whitelist[cloc->id].size() == 0) {
        for (uint i = 0; i < cloc->snps.size(); i++)
            if (t->nucs[cloc->snps[i]->col].fixed == false) {
                new_wl.insert(cloc->snps[i]->col);
                break;
            }

    } else {
        //
        // Otherwise, choose the first SNP that is already in the whitelist.
        //
        for (uint i = 0; i < cloc->snps.size(); i++) {
            if (this->_whitelist[cloc->id].count(cloc->snps[i]->col) == 0 ||
                t->nucs[cloc->snps[i]->col].fixed == true)
                continue;
            new_wl.insert(cloc->snps[i]->col);
            break;
        }
    }

    this->_whitelist[cloc->id] = new_wl;

    return 0;
}

int
LocusFilter::keep_random_snp(const CSLocus *cloc, const LocTally *t)
{
    set<int> new_wl, seen;
    int      index;

    if (this->_whitelist.count(cloc->id) == 0 && cloc->snps.size() == 0) {
        this->_whitelist[cloc->id] = new_wl;
        return 0;
    }

    //
    // If no specific SNPs are specified in the whitelist for this locus,
    // then all SNPs are included, choose a random variant.
    //
    if (this->_whitelist.count(cloc->id) == 0  || this->_whitelist[cloc->id].size() == 0) {
        do {
            index = rand() % cloc->snps.size();
            seen.insert(index);
        }  while (t->nucs[cloc->snps[index]->col].fixed == true && seen.size() < cloc->snps.size());

        if (t->nucs[cloc->snps[index]->col].fixed == false)
            new_wl.insert(cloc->snps[index]->col);

    } else {
        //
        // Otherwise, choose a SNP that is already in the whitelist randomly.
        //
        do {
            index = rand() % cloc->snps.size();
            seen.insert(index);
        } while (this->_whitelist.count(cloc->snps[index]->col) == 0 && seen.size() < cloc->snps.size());

        if (t->nucs[cloc->snps[index]->col].fixed == false)
            new_wl.insert(cloc->snps[index]->col);
    }

    this->_whitelist[cloc->id] = new_wl;

    return 0;
}

int
LocusFilter::prune_sites_with_whitelist(MetaPopInfo *mpopi, CSLocus *cloc, Datum **d, bool user_wl)
{
    if (user_wl == false && write_single_snp == false && write_random_snp == false)
        return 0;

    assert(this->_whitelist.count(cloc->id));

    //
    // When a user whitelist is not supplied but write_single_snp or write_random_snp is,
    // there can be homozygous loci that end up in the whitelist.
    //
    if (cloc->snps.size() == 0)
        return 0;

    if (user_wl == true && this->_whitelist[cloc->id].size() == 0)
        for (uint i = 0; i < cloc->snps.size(); i++)
            this->_whitelist[cloc->id].insert(cloc->snps[i]->col);

    prune_sites(cloc, d, this->_whitelist[cloc->id]);

    return 0;
}

int
LocusFilter::prune_sites(CSLocus *cloc, Datum **d, set<int> &keep)
{
    //
    // We want to prune out SNP objects that are not in the whitelist.
    //
    int              pos;
    vector<SNP *>    tmp;
    vector<uint>     cols;
    map<string, int> obshaps;
    map<string, int>::iterator sit;

    for (uint i = 0; i < cloc->snps.size(); i++) {
        if (keep.count(cloc->snps[i]->col) > 0) {
            tmp.push_back(cloc->snps[i]);
            cols.push_back(i);
        } else {
            //
            // Change the model calls in the samples to no longer contain this SNP.
            //
            pos = cloc->snps[i]->col;
            for (uint j = 0; j < this->_sample_cnt; j++) {
                if (d[j] == NULL || pos >= d[j]->len)
                    continue;
                if (d[j]->model != NULL) {
                    d[j]->model[pos] = 'U';
                }
            }

            delete cloc->snps[i];
        }
    }
    cloc->snps.clear();
    for (uint i = 0; i < tmp.size(); i++)
        cloc->snps.push_back(tmp[i]);

    map<string, int>::iterator it;
    char allele_old[id_len], allele_new[id_len];
    //
    // We need to adjust the catalog's list of haplotypes/alleles
    // for this locus to account for the pruned SNPs.
    //
    for (it = cloc->alleles.begin(); it != cloc->alleles.end(); it++) {
        strncpy(allele_old, it->first.c_str(), id_len - 1);
        allele_old[id_len - 1] = '\0';

        for (uint k = 0; k < cols.size(); k++)
            allele_new[k] = allele_old[cols[k]];
        allele_new[cols.size()] = '\0';
        obshaps[string(allele_new)] += it->second;
    }
    cloc->alleles.clear();
    for (sit = obshaps.begin(); sit != obshaps.end(); sit++) {
        cloc->alleles[sit->first] = sit->second;
    }
    obshaps.clear();

    cloc->populate_alleles();

    //
    // Now we need to adjust the matched haplotypes to sync to
    // the SNPs left in the catalog.
    //
    // Reducing the lengths of the haplotypes  may create
    // redundant (shorter) haplotypes, we need to remove these.
    //
    for (uint i = 0; i < this->_sample_cnt; i++) {
        if (d[i] == NULL) continue;

        for (uint j = 0; j < d[i]->obshap.size(); j++) {
            for (uint k = 0; k < cols.size(); k++)
                d[i]->obshap[j][k] = d[i]->obshap[j][cols[k]];
            d[i]->obshap[j][cols.size()] = '\0';
            obshaps[d[i]->obshap[j]] += d[i]->depth[j];
        }
        uint j = 0;
        for (sit = obshaps.begin(); sit != obshaps.end(); sit++) {
            strcpy(d[i]->obshap[j], sit->first.c_str());
            d[i]->depth[j] = sit->second;
            j++;
        }
        while (j < d[i]->obshap.size()) {
            delete [] d[i]->obshap[j];
            j++;
        }
        d[i]->obshap.resize(obshaps.size());
        d[i]->depth.resize(obshaps.size());
        obshaps.clear();
    }

    return 0;
}

bool
LocusFilter::prune_sites_with_filters(MetaPopInfo *mpopi, CSLocus *cloc, Datum **d, LocPopSum *s, ostream &log_fh)
{
    set<int>    site_list;
    vector<int> pop_prune_list;
    bool        sample_prune, maf_prune, het_prune, inc_prune;

    //
    // If this locus is fixed, ignore it.
    //
    if (cloc->snps.size() == 0)
        return false;

    const LocSum   *sum;
    const LocTally *t;

    t = s->meta_pop();
    for (uint i = 0; i < cloc->snps.size(); i++) {
        //
        // If the site is fixed, ignore it.
        //
        if (t->nucs[cloc->snps[i]->col].fixed == true)
            continue;

        sample_prune = false;
        maf_prune    = false;
        het_prune    = false;
        inc_prune    = false;
        pop_prune_list.clear();

        for (size_t p = 0; p < s->pop_cnt(); ++p) {
            sum = s->per_pop(p);

            if (sum->nucs[cloc->snps[i]->col].incompatible_site)
                inc_prune = true;

            else if (sum->nucs[cloc->snps[i]->col].num_indv == 0 ||
                     (double) sum->nucs[cloc->snps[i]->col].num_indv / (double) this->_pop_tot[p] < sample_limit)
                pop_prune_list.push_back(p);
        }

        //
        // Check how many populations have to be pruned out due to sample limit. If less than
        // population limit, prune them; if more than population limit, mark locus for deletion.
        //
        if ((mpopi->pops().size() - pop_prune_list.size()) < (uint) population_limit) {
            sample_prune = true;
        } else {
            for (size_t p : pop_prune_list) {
                if (s->per_pop(p)->nucs[cloc->snps[i]->col].num_indv == 0)
                    continue;

                const Pop& pop = mpopi->pops()[p];
                for (uint k = pop.first_sample; k <= pop.last_sample; k++) {
                    if (d[k] == NULL || cloc->snps[i]->col >= (uint) d[k]->len)
                        continue;
                    if (d[k]->model != NULL) {
                        d[k]->model[cloc->snps[i]->col] = 'U';
                    }
                }
            }
        }

        if (t->nucs[cloc->snps[i]->col].allele_cnt > 1) {
            //
            // Test for minor allele frequency.
            //
            if ((1 - t->nucs[cloc->snps[i]->col].p_freq) < minor_allele_freq)
                maf_prune = true;
            //
            // Test for observed heterozygosity.
            //
            if (t->nucs[cloc->snps[i]->col].obs_het > max_obs_het)
                het_prune = true;
        }

        if (maf_prune == false && het_prune == false && sample_prune == false && inc_prune == false) {
            site_list.insert(cloc->snps[i]->col);
        } else {
            this->_filtered_sites++;

            if (verbose) {
                log_fh << "pruned_polymorphic_site\t"
                       << cloc->id << "\t"
                       << cloc->loc.chr() << "\t"
                       << cloc->sort_bp(cloc->snps[i]->col) +1 << "\t"
                       << cloc->snps[i]->col << "\t";
                if (inc_prune)
                    log_fh << "incompatible_site\n";
                else if (sample_prune)
                    log_fh << "sample_limit\n";
                else if (maf_prune)
                    log_fh << "maf_limit\n";
                else if (het_prune)
                    log_fh << "obshet_limit\n";
                else
                    log_fh << "unknown_reason\n";
            }
        }
    }

    //
    // If no SNPs were retained for this locus, then mark it to be removed entirely.
    //
    if (site_list.size() == 0) {
        if (verbose)
            log_fh << "removed_locus\t"
                   << cloc->id << "\t"
                   << cloc->loc.chr() << "\t"
                   << cloc->sort_bp() +1 << "\t"
                   << 0 << "\tno_snps_remaining\n";
        this->_filtered_loci++;
        this->_batch_filtered_loci++;
        return true;
    } else {
        //
        // Remove the filtered sites from this locus.
        //
        this->prune_sites(cloc, d, site_list);
    }

    return false;
}

int
CatalogDists::accumulate_pre_filtering(const CSLocus *loc)
{
    size_t missing;

    if (this->_pre_valid.count(loc->hcnt) == 0)
        this->_pre_valid[loc->hcnt] = 1;
    else
        this->_pre_valid[loc->hcnt]++;

    if (this->_pre_confounded.count(loc->confounded_cnt) == 0)
        this->_pre_confounded[loc->confounded_cnt] = 1;
    else
        this->_pre_confounded[loc->confounded_cnt]++;

    missing = loc->cnt - loc->hcnt;

    if (this->_pre_absent.count(missing) == 0)
        this->_pre_absent[missing] = 1;
    else
        this->_pre_absent[missing]++;

    if (this->_pre_snps_per_loc.count(loc->snps.size()) == 0)
        this->_pre_snps_per_loc[loc->snps.size()] = 1;
    else
        this->_pre_snps_per_loc[loc->snps.size()]++;

    return 0;
}

int
CatalogDists::accumulate(const vector<LocBin *> &loci)
{
    const CSLocus *loc;
    size_t missing;

    for (uint i = 0; i < loci.size(); i++) {
        loc = loci[i]->cloc;

        if (this->_post_valid.count(loc->hcnt) == 0)
            this->_post_valid[loc->hcnt] = 1;
        else
            this->_post_valid[loc->hcnt]++;

        if (this->_post_confounded.count(loc->confounded_cnt) == 0)
            this->_post_confounded[loc->confounded_cnt] = 1;
        else
            this->_post_confounded[loc->confounded_cnt]++;

        missing = loc->cnt - loc->hcnt;

        if (this->_post_absent.count(missing) == 0)
            this->_post_absent[missing] = 1;
        else
            this->_post_absent[missing]++;

        if (this->_post_snps_per_loc.count(loc->snps.size()) == 0)
            this->_post_snps_per_loc[loc->snps.size()] = 1;
        else
            this->_post_snps_per_loc[loc->snps.size()]++;

    }

    return 0;
}

int
CatalogDists::write_results(ostream &log_fh)
{
    map<size_t, size_t>::iterator cnt_it;

    log_fh << "# Distribution of valid samples matched to a catalog locus prior to filtering.\n"
           << "# Valid samples at locus\tCount\n";
    for (cnt_it = this->_pre_valid.begin(); cnt_it != this->_pre_valid.end(); cnt_it++)
        log_fh << cnt_it->first << "\t" << cnt_it->second << "\n";

    log_fh << "\n# Distribution of confounded samples for each catalog locus prior to filtering.\n"
           << "# Confounded samples at locus\tCount\n";
    for (cnt_it = this->_pre_confounded.begin(); cnt_it != this->_pre_confounded.end(); cnt_it++)
        log_fh << cnt_it->first << "\t" << cnt_it->second << "\n";

    log_fh << "\n# Distribution of missing samples for each catalog locus prior to filtering.\n"
           << "# Absent samples at locus\tCount\n";
    for (cnt_it = this->_pre_absent.begin(); cnt_it != this->_pre_absent.end(); cnt_it++)
        log_fh << cnt_it->first << "\t" << cnt_it->second << "\n";

    log_fh << "\n# Distribution of the number of SNPs per catalog locus prior to filtering.\n"
           << "# Number SNPs\tNumber loci\n";
    for (cnt_it = this->_pre_snps_per_loc.begin(); cnt_it != this->_pre_snps_per_loc.end(); cnt_it++)
        log_fh << cnt_it->first << "\t" << cnt_it->second << "\n";
    log_fh << "\n";

    log_fh << "# Distribution of valid samples matched to a catalog locus after filtering.\n"
           << "# Valid samples at locus\tCount\n";
    for (cnt_it = this->_post_valid.begin(); cnt_it != this->_post_valid.end(); cnt_it++)
        log_fh << cnt_it->first << "\t" << cnt_it->second << "\n";

    log_fh << "\n# Distribution of confounded samples for each catalog locus after filtering.\n"
           << "# Confounded samples at locus\tCount\n";
    for (cnt_it = this->_post_confounded.begin(); cnt_it != this->_post_confounded.end(); cnt_it++)
        log_fh << cnt_it->first << "\t" << cnt_it->second << "\n";

    log_fh << "\n# Distribution of missing samples for each catalog locus after filtering.\n"
           << "# Absent samples at locus\tCount\n";
    for (cnt_it = this->_post_absent.begin(); cnt_it != this->_post_absent.end(); cnt_it++)
        log_fh << cnt_it->first << "\t" << cnt_it->second << "\n";

    log_fh << "\n# Distribution of the number of SNPs per catalog locus after filtering.\n"
           << "# Number SNPs\tNumber loci\n";
    for (cnt_it = this->_post_snps_per_loc.begin(); cnt_it != this->_post_snps_per_loc.end(); cnt_it++)
        log_fh << cnt_it->first << "\t" << cnt_it->second << "\n";
    log_fh << "\n";

    return 0;
}

int
log_haplotype_cnts(map<int, CSLocus *> &catalog, ofstream &log_fh)
{
    map<int, CSLocus *>::iterator it;
    map<int, int> valid, absent, confounded;

    CSLocus *loc;
    int missing;

    for (it = catalog.begin(); it != catalog.end(); it++) {
        loc = it->second;

        if (valid.count(loc->hcnt) == 0)
            valid[loc->hcnt] = 1;
        else
            valid[loc->hcnt]++;

        if (confounded.count(loc->confounded_cnt) == 0)
            confounded[loc->confounded_cnt] = 1;
        else
            confounded[loc->confounded_cnt]++;

        missing = loc->cnt - loc->hcnt;

        if (absent.count(missing) == 0)
            absent[missing] = 1;
        else
            absent[missing]++;
    }

    map<int, int>::iterator cnt_it;

    log_fh << "# Distribution of valid loci matched to catalog locus.\n"
           << "# Valid samples at locus\tCount\n";
    for (cnt_it = valid.begin(); cnt_it != valid.end(); cnt_it++)
        log_fh << cnt_it->first << "\t" << cnt_it->second << "\n";

    log_fh << "# Distribution of confounded loci at catalog locus.\n"
           << "# Confounded samples at locus\tCount\n";
    for (cnt_it = confounded.begin(); cnt_it != confounded.end(); cnt_it++)
        log_fh << cnt_it->first << "\t" << cnt_it->second << "\n";

    log_fh << "# Distribution of missing loci at catalog loci.\n"
           << "# Absent samples at locus\tCount\n";
    for (cnt_it = absent.begin(); cnt_it != absent.end(); cnt_it++)
        log_fh << cnt_it->first << "\t" << cnt_it->second << "\n";

    return 0;
}

int
tabulate_locus_haplotypes(CSLocus *cloc, Datum **d, int sample_cnt)
{
    double mean = 0.0;
    double cnt  = 0.0;

    for (int i = 0; i < sample_cnt; i++) {
        if (d[i] == NULL)
            continue;

        if (d[i]->obshap.size() > 1)
            cloc->marker = "heterozygous";

        mean += d[i]->lnl;
        cnt++;
    }

    if (cloc->marker.length() > 0) {
        create_genotype_map(cloc, d, sample_cnt);
        call_population_genotypes(cloc, d, sample_cnt);
    }

    cloc->lnl = mean / cnt;

    return 0;
}

int
tabulate_haplotypes(map<int, CSLocus *> &catalog, PopMap<CSLocus> *pmap)
{
    map<int, CSLocus *>::iterator it;
    vector<char *>::iterator hit;
    Datum  **d;
    CSLocus *loc;
    double   mean, cnt;

    for (it = catalog.begin(); it != catalog.end(); it++) {
        loc = it->second;
        d   = pmap->locus(loc->id);

        mean = 0.0;
        cnt  = 0.0;

        for (int i = 0; i < pmap->sample_cnt(); i++) {
            if (d[i] == NULL)
                continue;

            if (d[i]->obshap.size() > 1)
                loc->marker = "heterozygous";

            mean += d[i]->lnl;
            cnt++;
        }

        if (loc->marker.length() > 0) {
            create_genotype_map(loc, d, pmap->sample_cnt());
            call_population_genotypes(loc, d, pmap->sample_cnt());
        }

        loc->lnl = mean / cnt;
    }

    return 0;
}

/*
int
merge_shared_cutsite_loci(map<int, CSLocus *> &catalog,
                          PopMap<CSLocus> *pmap, PopSum<CSLocus> *psum,
                          map<int, pair<merget, int> > &merge_map,
                          ofstream &log_fh)
{
    map<string, vector<CSLocus *> >::iterator it;
    CSLocus *cur, *next;
    Datum  **d_1, **d_2;
    double   prune_pct;
    uint unmergable, tot_loci, tot_samp;
    uint success           = 0;
    uint failure           = 0;
    uint overlap           = 0;
    uint simple_merge_cnt  = 0;
    uint complex_merge_cnt = 0;
    uint missing_samps_cnt = 0;
    uint phase_fail_cnt    = 0;
    uint nomapping_cnt     = 0;
    uint multimapping_cnt  = 0;
    uint multifails_cnt    = 0;

    tot_loci = pmap->loci_cnt();

    set<int> loci_to_destroy;
    map<int, int> missing_samps_dist;

    cerr << "To merge adjacent loci at least " << merge_prune_lim * 100 << "% of samples must have both adjacent loci;"
         << " the remaining " << 100 - (merge_prune_lim * 100) << "% of individuals will be pruned.\n"
         << "Attempting to merge adjacent loci that share a cutsite...";

    if (verbose)
        log_fh << "\n#\n# List of locus pairs that share a cutsite that failed to merge because they could not be phased.\n#\n";

    //
    // Iterate over each chromosome.
    //
    for (it = pmap->ordered_loci_nconst().begin(); it != pmap->ordered_loci_nconst().end(); it++) {
        //
        // Iterate over each ordered locus on this chromosome.
        //
        next = it->second[0];
        for (uint pos = 1; pos < it->second.size(); pos++) {
            cur  = next;
            next = it->second[pos];

            //
            // Do these two loci overlap?
            //   +Must occur on opposite strands
            //   +Must overlap according to the length of the cutsite.
            //
            if (((cur->loc.strand == strand_minus && next->loc.strand == strand_plus) &&
                 ((int) (cur->loc.bp  - next->loc.bp + 1) == renz_olap[enz])) ||
                ((cur->loc.strand == strand_plus  && next->loc.strand == strand_minus) &&
                 ((int) (next->loc.bp - cur->loc.bp  + 1) == renz_olap[enz]))) {
                overlap++;

                d_1        = pmap->locus(cur->id);
                d_2        = pmap->locus(next->id);
                unmergable = 0;
                tot_samp   = 0;

                //
                // Check if all members of the population contain these two loci (or are missing both).
                //
                for (int i = 0; i < pmap->sample_cnt(); i++) {
                    if (d_1[i] != NULL || d_2[i] != NULL)
                        tot_samp++;
                    if ((d_1[i] != NULL && d_2[i] == NULL) ||
                        (d_1[i] == NULL && d_2[i] != NULL))
                        unmergable++;
                }

                prune_pct = (double) (tot_samp - unmergable) / (double) tot_samp;

                //
                // If some of the individuals only have one locus and not the other, prune them out.
                //
                if (prune_pct < 1.0 && prune_pct >= merge_prune_lim) {
                    for (int i = 0; i < pmap->sample_cnt(); i++)
                        if (d_1[i] != NULL && d_2[i] == NULL) {
                            delete d_1[i];
                            d_1[i] = NULL;
                        } else if (d_1[i] == NULL && d_2[i] != NULL) {
                            delete d_2[i];
                            d_2[i] = NULL;
                        }
                }

                //
                // If possible, merge the two loci together.
                //
                if (prune_pct < merge_prune_lim) {
                    int pct = (int) (prune_pct * 100);
                    missing_samps_dist[pct]++;
                    if (verbose) log_fh << "Missing samples, Sink Locus: " << cur->id << "; Source Locus: " << next->id << "; "
                                        << pct << "% present (" << 100 - pct << "% missing)\n";
                    missing_samps_cnt++;
                    failure++;
                    continue;
                }

                phaset res = merge_and_phase_loci(pmap, cur, next, loci_to_destroy, log_fh);
                switch(res) {
                case multiple_fails:
                    if (verbose) log_fh << "Failed to phase, Sink Locus: " << cur->id << "; Source Locus: " << next->id << "; "
                                        << "multiple failures\n";
                    multifails_cnt++;
                    phase_fail_cnt++;
                    failure++;
                    break;
                case multimapping_fail:
                    if (verbose) log_fh << "Failed to phase, Sink Locus: " << cur->id << "; Source Locus: " << next->id << "; "
                                        << "multimapping in one or more individuals\n";
                    multimapping_cnt++;
                    phase_fail_cnt++;
                    failure++;
                    break;
                case nomapping_fail:
                    if (verbose) log_fh << "Failed to phase, Sink Locus: " << cur->id << "; Source Locus: " << next->id << "; "
                                        << "no mapping in one or more individuals\n";
                    nomapping_cnt++;
                    phase_fail_cnt++;
                    failure++;
                    break;
                case complex_phase:
                    if (verbose) log_fh << "Phased Sink Locus: " << cur->id << "; Source Locus: " << next->id << "; "
                                        << "a complex phasing operation.\n";
                    complex_merge_cnt++;
                    success++;
                    merge_map[cur->id]  = make_pair(merge_sink, next->id);
                    merge_map[next->id] = make_pair(merge_src, cur->id);
                    break;
                case simple_merge:
                    if (verbose) log_fh << "Phased Sink Locus: " << cur->id << "; Source Locus: " << next->id << "; "
                                        << "a simple merge operation.\n";
                    simple_merge_cnt++;
                    success++;
                    merge_map[cur->id]  = make_pair(merge_sink, next->id);
                    merge_map[next->id] = make_pair(merge_src, cur->id);
                    break;
                default:
                    cerr << "Warning: Merge failure.\n";
                    break;
                }
            }
        }
    }

    //
    // Remove those loci that have been merged from both the popualtion map and catalog.
    //
    set<int> emptyset;
    pmap->prune(loci_to_destroy);
    reduce_catalog(catalog, emptyset, loci_to_destroy);

    cerr << "done.\n"
         << "Of " << tot_loci << " loci, "
         << overlap << " pairs share a cutsite; "
         << success << " pairs were merged; "
         << failure << " pairs failed to merge; "
         << pmap->loci_cnt() << " remaining loci.\n"
         << "  Of those merged, " << simple_merge_cnt << " required only a simple merge without phasing; "
         << "while " << complex_merge_cnt << " required phasing.\n"
         << "  Of those that failed to merge, " << missing_samps_cnt << " were missing one of the two haplotypes in one or more samples; "
         << "while " << phase_fail_cnt << " failed to be phased.\n"
         << "    Of those that failed to phase, " << nomapping_cnt << " failed due to a lack of haplotype mappings; "
         << multimapping_cnt << " failed due to multiple haplotype mappings; " << multifails_cnt << " failed due to both.\n";

    log_fh << "\n#\n# Merging adjacent loci with a shared restriction enzyme cutsite\n#\n"
           << "Of " << tot_loci << " loci, "
           << overlap << " pairs share a cutsite; "
           << success << " pairs were merged; "
           << failure << " pairs failed to merge; "
           << pmap->loci_cnt() << " remaining loci.\n"
           << "  Of those merged, " << simple_merge_cnt << " required only a simple merge without phasing; "
           << "while " << complex_merge_cnt << " required phasing.\n"
           << "  Of those that failed to merge, " << missing_samps_cnt << " were missing one of the two haplotypes in one or more samples; "
           << "while " << phase_fail_cnt << " failed to be phased.\n"
           << "    Of those that failed to phase, " << nomapping_cnt << " failed due to a lack of haplotype mappings; "
           << multimapping_cnt << " failed due to multiple haplotype mappings; " << multifails_cnt << " failed due to both.\n";
    log_fh << "#\n# Distribution of loci with samples missing one of two loci to be merged\n"
           << "# Percent samples with both loci present\tNumber of cases\n";
    map<int, int>::iterator mit;
    for (mit = missing_samps_dist.begin(); mit != missing_samps_dist.end(); mit++)
        log_fh << mit->first << "\t" << mit->second << "\n";
    log_fh << "\n";

    return 0;
}
*/
phaset
merge_and_phase_loci(PopMap<CSLocus> *pmap, CSLocus *cur, CSLocus *next,
                     set<int> &loci_to_destroy,
                     ofstream &log_fh)
{
    Datum **d_1 = pmap->locus(cur->id);
    Datum **d_2 = pmap->locus(next->id);

    set<int>    phased_results;
    set<string> phased_haplotypes;
    string      merged_hap;
    char       *h_1, *h_2;
    int         merge_type;

    if (verbose) log_fh << "Attempting to phase source locus " << cur->id << " with sink locus " << next->id << "\n";

    int sample_cnt        = 0;
    int phased_sample_cnt = 0;
    //
    // Take a census of the already phased haplotypes. We have phased haplotypes
    // if for individual i:
    //   1. d_1 has a single haplotype and d_2 has a single haplotype
    //   2. d_1 has a single haplotpye and d_2 has multiple haplotypes
    //   3. d_1 has multiple haplotpyes and d_2 has a single haplotype
    //
    // If one or both of the loci have no SNPs, then the haplotype is
    // recorded as "consensus." Check that condition before we start merging.
    //
    if (cur->snps.size() > 0 && next->snps.size() > 0)
        merge_type = 0;
    else if (cur->snps.size() == 0)
        merge_type = 1;
    else if (next->snps.size() == 0)
        merge_type = 2;
    else
        merge_type = 3;

    for (int i = 0; i < pmap->sample_cnt(); i++) {
        if (d_1[i] == NULL || d_2[i] == NULL)
            continue;
        else if (d_1[i]->obshap.size() > 1 && d_2[i]->obshap.size() > 1)
            continue;
        else {
            for (uint j = 0; j < d_1[i]->obshap.size(); j++) {
                for (uint k = 0; k < d_2[i]->obshap.size(); k++) {
                    switch (merge_type) {
                    case 0:
                        merged_hap = string(d_1[i]->obshap[j]) + string(d_2[i]->obshap[k]);
                        break;
                    case 1:
                        merged_hap = string(d_2[i]->obshap[k]);
                        break;
                    case 2:
                        merged_hap = string(d_1[i]->obshap[j]);
                        break;
                    case 3:
                    default:
                        merged_hap = "consensus";
                        break;
                    }
                    phased_haplotypes.insert(merged_hap);
                    // cerr << "Phasing: '" << d_1[i]->obshap[j] << "' + '" << d_2[i]->obshap[k] << "' => '" << merged_hap << "'\n";
                }
            }
            phased_sample_cnt++;
            sample_cnt++;
        }
    }

    //
    // Indicate that these two loci had a simple merge, with no phasing necessary.
    //
    phased_results.insert(simple_merge);

    //
    // Now we need to check if we can phase the remaining haplotypes.
    //
    for (int i = 0; i < pmap->sample_cnt(); i++) {
        if (d_1[i] == NULL || d_2[i] == NULL)
            continue;
        else if (d_1[i]->obshap.size() > 1 && d_2[i]->obshap.size() > 1) {
            // cerr << "Attempting to phase individual " << i << ": " << d_1[i]->id << " / " << d_2[i]->id << "\n";

            sample_cnt++;
            //
            // We should be able to find a sinlge phasing mapping for each haplotype from d_1 to d_2
            // that includes all the haplotypes in these two loci.
            //
            vector<pair<char *, char *> > seen_phased;
            uint tot_obshap = d_1[i]->obshap.size() + d_2[i]->obshap.size();
            uint phased_cnt = 0;
            for (uint j = 0; j < d_1[i]->obshap.size(); j++) {
                for (uint k = 0; k < d_2[i]->obshap.size(); k++) {
                    // cerr << "  " << d_1[i]->obshap[j] << " + " << d_2[i]->obshap[k];
                    //
                    // Record each pair of haplotypes that has been seen phased previously.
                    //
                    if (phased_haplotypes.count(string(d_1[i]->obshap[j]) + string(d_2[i]->obshap[k]))) {
                        seen_phased.push_back(make_pair(d_1[i]->obshap[j], d_2[i]->obshap[k]));
                        // cerr << " => " << d_1[i]->obshap[j] << d_2[i]->obshap[k];
                    }
                    // cerr << "\n";
                }
            }
            //
            // Now, we will iterate over all sets of phased haplotypes and look
            // for combinations that use all four individual haplotypes.
            //
            for (uint j = 0; j < seen_phased.size(); j++) {
                for (uint k = j; k < seen_phased.size(); k++) {
                    set<char *> incorporated_haplotypes;
                    //
                    // Count the number of distinct char pointers. If this combination
                    // of haplotypes includes all unphased haplotypes, count it.
                    //
                    incorporated_haplotypes.insert(seen_phased[j].first);
                    incorporated_haplotypes.insert(seen_phased[j].second);
                    incorporated_haplotypes.insert(seen_phased[k].first);
                    incorporated_haplotypes.insert(seen_phased[k].second);
                    if (incorporated_haplotypes.size() == tot_obshap)
                        phased_cnt++;
                }
            }

            //
            // If one pair of haplotypes is mapped, but the other is not, assume the second pair or
            // haplotypes must be phased by process of elimination.
            //
            if (phased_cnt == 0 && seen_phased.size() == 1) {
                h_1 = seen_phased[0].first  == d_1[i]->obshap[1] ?
                    d_1[i]->obshap[0] : d_1[i]->obshap[1];
                h_2 = seen_phased[0].second == d_2[i]->obshap[1] ?
                    d_2[i]->obshap[0] : d_2[i]->obshap[1];
                phased_haplotypes.insert(string(h_1) + string(h_2));
                phased_cnt++;
                // cerr << "  Phasing: '" << hap_1 << "' + '" << hap_2 << "' => '" << string(hap_1) + string(hap_2) << "'\n";
            }

            if (phased_cnt == 0) {
                phased_results.insert(nomapping_fail);
                if (verbose) log_fh << "    Locus NOT phased in individual " << i << "; loci " << d_1[i]->id << " / " << d_2[i]->id << " no mapping found.\n";
            } else if (phased_cnt == 1) {
                phased_sample_cnt++;
                phased_results.insert(complex_phase);
            } else {
                phased_results.insert(multimapping_fail);
                if (verbose) log_fh << "    Locus NOT phased in individual " << i << "; loci " << d_1[i]->id << " / " << d_2[i]->id << " multiple mappings found.\n";
            }
        }
    }

    if (phased_sample_cnt != sample_cnt) {
        if (phased_results.count(nomapping_fail) > 0 &&
            phased_results.count(multimapping_fail) > 0)
            return multiple_fails;
        else if (phased_results.count(nomapping_fail) > 0)
            return nomapping_fail;
        else if (phased_results.count(multimapping_fail) > 0)
            return multimapping_fail;
        else {
            cerr << "WE SHOULD NOT GET HERE\n";
            return merge_failure;
        }
    }

    //
    // Okay, merge these two loci together.
    //
    if (!merge_datums(pmap->sample_cnt(), cur->len, d_1, d_2, phased_haplotypes, merge_type))
        return merge_failure;

    //
    // Merge the catalog entries together.
    //
    if (!merge_csloci(cur, next, phased_haplotypes))
        return merge_failure;

    //
    // Mark the merged locus for destruction.
    //
    loci_to_destroy.insert(next->id);

    if (phased_results.count(complex_phase) > 0)
        return complex_phase;
    return simple_merge;
}

int
merge_csloci(CSLocus *sink, CSLocus *src, set<string> &phased_haplotypes)
{
    //
    // We assume that we are merging two loci: one on the negative strand, one on the
    // positive. We will keep the sink cslocus and delete the src cslocus.
    //   -> The sink cslocus is assumed to be on the negative strand.
    //

    //
    // 1. Reverse complement the SNP coordinates in the sink locus so that they are
    //    enumerated on the positive strand. Complement the alleles as well.
    //
    for (uint j = 0; j < sink->snps.size(); j++) {
        sink->snps[j]->col    = sink->len - sink->snps[j]->col - 1;
        sink->snps[j]->rank_1 = reverse(sink->snps[j]->rank_1);
        sink->snps[j]->rank_2 = reverse(sink->snps[j]->rank_2);
        sink->snps[j]->rank_3 = reverse(sink->snps[j]->rank_3);
        sink->snps[j]->rank_4 = reverse(sink->snps[j]->rank_4);
    }

    //
    // 2. Adjust the SNP coordinates in the src locus to account for the now, longer length.
    //
    for (uint j = 0; j < src->snps.size(); j++)
        src->snps[j]->col = sink->len + src->snps[j]->col - renz_olap[enz];

    //
    // 3. Combine SNPs between the two catalog loci: add the SNPs from the sink (formerly on the
    //    negative strand) in reverse order, followed by the SNPs from the src.
    //
    vector<SNP *> tmpsnp;
    for (int j = (int) sink->snps.size() - 1; j >= 0; j--)
        tmpsnp.push_back(sink->snps[j]);
    for (uint j = 0; j < src->snps.size(); j++)
        tmpsnp.push_back(src->snps[j]);
    sink->snps.clear();
    for (uint j = 0; j < tmpsnp.size(); j++)
        sink->snps.push_back(tmpsnp[j]);

    //
    // 4. Adjust the genomic location of the sink locus.
    //
    uint bp = sink->sort_bp();
    sink->loc.bp     = bp;
    sink->loc.strand = strand_plus;

    //
    // 5. Adjust the length of the sequence.
    //
    sink->len += src->len - renz_olap[enz];

    //
    // 6. Merge the consensus sequence together.
    //
    char *new_con = rev_comp(sink->con);
    delete [] sink->con;
    sink->con = new_con;
    new_con   = new char[sink->len + 1];
    strcpy(new_con, sink->con);
    delete [] sink->con;
    sink->con = new_con;
    new_con  += src->len - renz_olap[enz];
    strcpy(new_con, src->con);

    //
    // 7. Record the now phased haplotypes.
    //
    sink->alleles.clear();
    set<string>::iterator it;
    for (it = phased_haplotypes.begin(); it != phased_haplotypes.end(); it++)
        sink->alleles[*it] = 0;

    // cerr << "CSLocus " << sink->id << ":\n"
    //   << "Length: " << sink->len << "; Chr: " << sink->loc.chr << "; BP: " << sink->sort_bp() << "; strand: " << (sink->loc.strand == strand_plus ? "+" : "-") << "\n"
    //   << "  SNPs:\n";
    // for (uint j = 0; j < sink->snps.size(); j++)
    //  cerr << "    Col: " << sink->snps[j]->col
    //       << "    Rank 1: " << sink->snps[j]->rank_1
    //       << "    Rank 2: " << sink->snps[j]->rank_2 << "\n";
    // cerr << "  Alleles:\n";
    // map<string, int>::iterator ait;
    // for (ait = sink->alleles.begin(); ait != sink->alleles.end(); ait++)
    //  cerr << "    " << ait->first << "\n";

    return 1;
}

int
merge_datums(int sample_cnt,
             int sink_locus_len,
             Datum **sink, Datum **src,
             set<string> &phased_haplotypes,
             int merge_type)
{
    char           tmphap[id_len], *new_hap;
    uint           haplen, model_len, offset;
    vector<SNP *>  tmpsnp;
    vector<string> tmpobshap;
    vector<int>    tmpobsdep;

    //
    // We assume that we are merging two loci: one on the negative strand, one on the
    // positive. We will keep the sink datum and delete the src datum.
    //   -The sink datum is assumed to be on the negative strand.
    //
    for (int i = 0; i < sample_cnt; i++) {
        if (sink[i] == NULL && src[i] == NULL)
            continue;
        else if (sink[i] == NULL || src[i] == NULL)
            cerr << "Unexpected condition in merging datums: one datum is NULL while the other is not.\n";

        //
        // 1. Reverse complement the observed haplotypes in the sink locus.
        //
        haplen = strlen(sink[i]->obshap[0]);
        for (uint j = 0; j < sink[i]->obshap.size(); j++) {
            for (uint k = 0; k < haplen; k++)
                tmphap[k] = reverse(sink[i]->obshap[j][haplen - k - 1]);
            tmphap[haplen] = '\0';
            strcpy(sink[i]->obshap[j], tmphap);
        }
    }

    //
    // 2. Combine observed haplotypes between the two datums while phasing them.
    //    2.1 First combine the haplotypes from samples that are already in phase.
    //
    string      merged_hap;
    vector<int> to_be_phased;
    phased_haplotypes.clear();
    for (int i = 0; i < sample_cnt; i++) {
        if (sink[i] == NULL && src[i] == NULL)
            continue;

        if (sink[i]->obshap.size() > 1 && src[i]->obshap.size() > 1) {
            to_be_phased.push_back(i);
            continue;
        } else {
            tmpobshap.clear();
            tmpobsdep.clear();
            for (uint j = 0; j < sink[i]->obshap.size(); j++) {
                for (uint k = 0; k < src[i]->obshap.size(); k++) {
                    switch (merge_type) {
                    case 0:
                        merged_hap = string(sink[i]->obshap[j]) + string(src[i]->obshap[k]);
                        break;
                    case 1:
                        merged_hap = string(src[i]->obshap[j]);
                        break;
                    case 2:
                        merged_hap = string(sink[i]->obshap[j]);
                        break;
                    case 3:
                    default:
                        merged_hap = "consensus";
                        break;
                    }
                    phased_haplotypes.insert(merged_hap);
                    tmpobshap.push_back(merged_hap);
                    tmpobsdep.push_back((sink[i]->depth[j] + src[i]->depth[k]) / 2);
                }
            }
            sink[i]->depth.clear();
            for (uint j = 0; j < sink[i]->obshap.size(); j++)
                delete [] sink[i]->obshap[j];
            sink[i]->obshap.clear();
            for (uint j = 0; j < tmpobshap.size(); j++) {
                new_hap = new char[tmpobshap[j].length() + 1];
                strcpy(new_hap, tmpobshap[j].c_str());
                sink[i]->obshap.push_back(new_hap);
                sink[i]->depth.push_back(tmpobsdep[j]);
            }
        }
    }
    //
    //    2.2 Phase and combine the haplotypes from the remaining samples.
    //
    int index;
    for (uint i = 0; i < to_be_phased.size(); i++) {
        index = to_be_phased[i];
        tmpobshap.clear();
        tmpobsdep.clear();

        vector<pair<char *, char *> > seen_phased;
        uint tot_obshap = sink[index]->obshap.size() + src[index]->obshap.size();

        for (uint j = 0; j < sink[index]->obshap.size(); j++) {
            for (uint k = 0; k < src[index]->obshap.size(); k++) {
                if (phased_haplotypes.count(string(sink[index]->obshap[j]) + string(src[index]->obshap[k])))
                    seen_phased.push_back(make_pair(sink[index]->obshap[j], src[index]->obshap[k]));
            }
        }

        for (uint j = 0; j < seen_phased.size(); j++) {
            for (uint k = j; k < seen_phased.size(); k++) {
                set<char *> incorporated_haplotypes;
                incorporated_haplotypes.insert(seen_phased[j].first);
                incorporated_haplotypes.insert(seen_phased[j].second);
                incorporated_haplotypes.insert(seen_phased[k].first);
                incorporated_haplotypes.insert(seen_phased[k].second);
                if (incorporated_haplotypes.size() == tot_obshap) {
                    tmpobshap.push_back(string(seen_phased[j].first) + string(seen_phased[j].second));
                    tmpobshap.push_back(string(seen_phased[k].first) + string(seen_phased[k].second));
                    //tmpobsdep.push_back((sink[index]->depth[j] + src[index]->depth[k]) / 2);
                }
            }
        }

        sink[index]->depth.clear();
        for (uint j = 0; j < sink[index]->obshap.size(); j++)
            delete [] sink[index]->obshap[j];
        sink[index]->obshap.clear();
        for (uint j = 0; j < tmpobshap.size(); j++) {
            new_hap = new char[tmpobshap[j].length() + 1];
            strcpy(new_hap, tmpobshap[j].c_str());
            sink[index]->obshap.push_back(new_hap);
            // sink[index]->depth.push_back(tmpobsdep[j]);
        }
    }

    //
    // 3. Merge model calls; Set the length; combine the two depth and lnl measures together.
    //
    string model_calls;
    char  *p;

    for (int i = 0; i < sample_cnt; i++) {
        if (sink[i] == NULL && src[i] == NULL)
            continue;

        //
        // Merge the two strings of model calls together.
        // We need to check if the locus for this individual is shorter than the catalog
        // locus. If so, we need to expand out the model call array to be the proper length.
        //
        reverse_string(sink[i]->model);
        offset = 0;
        model_calls.clear();
        if (sink_locus_len > sink[i]->len) {
            offset = sink_locus_len - sink[i]->len;
            model_calls.assign(offset, 'N');
        }
        model_len = offset + sink[i]->len + src[i]->len - renz_olap[enz];
        model_calls.append(sink[i]->model);
        delete [] sink[i]->model;
        sink[i]->model = new char[model_len + 1];
        strcpy(sink[i]->model, model_calls.c_str());
        p  = sink[i]->model;
        p += offset + sink[i]->len - renz_olap[enz];
        strcpy(p, src[i]->model);

        sink[i]->len       = model_len;
        sink[i]->tot_depth = (sink[i]->tot_depth + src[i]->tot_depth) / 2;
        sink[i]->lnl       = (sink[i]->lnl + src[i]->lnl) / 2.0;

        //
        // Record which datum was merged into this one.
        //
        sink[i]->merge_partner = src[i]->id;
    }

    return 1;
}

int
create_genotype_map(CSLocus *locus, Datum **d, int sample_cnt)
{
    //
    // Create a genotype map. For any set of haplotypes, this routine will
    // assign each haplotype to a genotype, e.g. given the haplotypes
    // 'AC' and 'GT' in the population, this routine will assign 'AC' == 'a'
    // and 'GT' == 'b'. If an individual is homozygous for 'AC', they will be
    // assigned an 'aa' genotype.
    //
    //cerr << "Creating genotype map for catalog ID " << locus->id  << ", marker: " << locus->marker << ".\n";

    char gtypes[26] ={'a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j',
                      'k', 'l', 'm', 'n', 'o', 'p', 'q', 'r', 's', 't',
                      'u', 'v', 'w', 'x', 'y', 'z'};

    map<string, int> haplotypes;
    map<string, int>::iterator k;
    vector<pair<string, int> > sorted_haplotypes;

    for (int i = 0; i < sample_cnt; i++) {

        if (d[i] != NULL)
            for (uint n = 0; n < d[i]->obshap.size(); n++)
                haplotypes[d[i]->obshap[n]]++;
    }

    //
    // Check that there are not more haplotypes than we have encodings.
    //
    if (haplotypes.size() > 26) return 0;

    //
    // Sort the haplotypes map by value
    //
    for (k = haplotypes.begin(); k != haplotypes.end(); k++)
        sorted_haplotypes.push_back(*k);
    sort(sorted_haplotypes.begin(), sorted_haplotypes.end(), hap_compare);

    for (uint n = 0, index = 0; n < sorted_haplotypes.size() && index <= 26; n++, index++) {
        locus->gmap[sorted_haplotypes[n].first] = gtypes[index];
        //cerr << "GMAP: " << sorted_haplotypes[n].first << " == " << gtypes[index] << "\n";
    }

    return 0;
}

int
call_population_genotypes(CSLocus *locus, Datum **d, int sample_cnt)
{
    //
    // Fetch the array of observed haplotypes from the population
    //
    for (int i = 0; i < sample_cnt; i++) {
        if (d[i] == NULL)
            continue;

        vector<string> gtypes;
        string gtype;

        //cerr << "Sample Id: " << pmap->rev_sample_index(i) << "\n";

        for (uint j = 0; j < d[i]->obshap.size(); j++) {
            //
            // Impossible allele encountered.
            //
            if (locus->gmap.count(d[i]->obshap[j]) == 0) {
                gtypes.clear();
                gtypes.push_back("-");
                goto impossible;
            }

            gtypes.push_back(locus->gmap[d[i]->obshap[j]]);
            //cerr << "  Observed Haplotype: " << d[i]->obshap[j] << ", Genotype: " << locus->gmap[d[i]->obshap[j]] << "\n";
        }

    impossible:
        sort(gtypes.begin(), gtypes.end());
        for (uint j = 0; j < gtypes.size(); j++) {
            gtype += gtypes[j];
            //cerr << "  Adding genotype to string: " << gtypes[j] << "; " << gtype << "\n";
        }

        string m = gtype.length() == 1 ?
            gtype + gtype : gtype;

        d[i]->gtype = new char[m.length() + 1];
        strcpy(d[i]->gtype, m.c_str());

        if (m != "-")
            locus->gcnt++;

        //cerr << "Assigning datum, marker: " << locus->marker << ", string: " << m << ", haplotype: " << d[i]->obshap[0] << ", gtype: " << gtype << "\n";
     }

    return 0;
}

SumStatsSummary::SumStatsSummary(size_t pop_cnt)
{
    this->_pop_cnt           = pop_cnt;
    this->_private_cnt       = new int[this->_pop_cnt];
    this->_n                 = new double[this->_pop_cnt];
    this->_var_sites         = new double[this->_pop_cnt];

    this->_num_indv_mean     = new double[this->_pop_cnt];
    this->_num_indv_acc_mean = new double[this->_pop_cnt];
    this->_num_indv_var      = new double[this->_pop_cnt];
    this->_p_mean            = new double[this->_pop_cnt];
    this->_p_acc_mean        = new double[this->_pop_cnt];
    this->_p_var             = new double[this->_pop_cnt];
    this->_obs_het_mean      = new double[this->_pop_cnt];
    this->_obs_het_acc_mean  = new double[this->_pop_cnt];
    this->_obs_het_var       = new double[this->_pop_cnt];
    this->_obs_hom_mean      = new double[this->_pop_cnt];
    this->_obs_hom_acc_mean  = new double[this->_pop_cnt];
    this->_obs_hom_var       = new double[this->_pop_cnt];
    this->_exp_het_mean      = new double[this->_pop_cnt];
    this->_exp_het_acc_mean  = new double[this->_pop_cnt];
    this->_exp_het_var       = new double[this->_pop_cnt];
    this->_exp_hom_mean      = new double[this->_pop_cnt];
    this->_exp_hom_acc_mean  = new double[this->_pop_cnt];
    this->_exp_hom_var       = new double[this->_pop_cnt];
    this->_pi_mean           = new double[this->_pop_cnt];
    this->_pi_acc_mean       = new double[this->_pop_cnt];
    this->_pi_var            = new double[this->_pop_cnt];
    this->_fis_mean          = new double[this->_pop_cnt];
    this->_fis_acc_mean      = new double[this->_pop_cnt];
    this->_fis_var           = new double[this->_pop_cnt];

    this->_n_all                 = new double[this->_pop_cnt];
    this->_num_indv_mean_all     = new double[this->_pop_cnt];
    this->_num_indv_acc_mean_all = new double[this->_pop_cnt];
    this->_num_indv_var_all      = new double[this->_pop_cnt];
    this->_p_mean_all            = new double[this->_pop_cnt];
    this->_p_acc_mean_all        = new double[this->_pop_cnt];
    this->_p_var_all             = new double[this->_pop_cnt];
    this->_obs_het_mean_all      = new double[this->_pop_cnt];
    this->_obs_het_acc_mean_all  = new double[this->_pop_cnt];
    this->_obs_het_var_all       = new double[this->_pop_cnt];
    this->_obs_hom_mean_all      = new double[this->_pop_cnt];
    this->_obs_hom_acc_mean_all  = new double[this->_pop_cnt];
    this->_obs_hom_var_all       = new double[this->_pop_cnt];
    this->_exp_het_mean_all      = new double[this->_pop_cnt];
    this->_exp_het_acc_mean_all  = new double[this->_pop_cnt];
    this->_exp_het_var_all       = new double[this->_pop_cnt];
    this->_exp_hom_mean_all      = new double[this->_pop_cnt];
    this->_exp_hom_acc_mean_all  = new double[this->_pop_cnt];
    this->_exp_hom_var_all       = new double[this->_pop_cnt];
    this->_pi_mean_all           = new double[this->_pop_cnt];
    this->_pi_acc_mean_all       = new double[this->_pop_cnt];
    this->_pi_var_all            = new double[this->_pop_cnt];
    this->_fis_mean_all          = new double[this->_pop_cnt];
    this->_fis_acc_mean_all      = new double[this->_pop_cnt];
    this->_fis_var_all           = new double[this->_pop_cnt];

    this->_sq_n     = new double[this->_pop_cnt];
    this->_sq_n_all = new double[this->_pop_cnt];

    for (uint j = 0; j < this->_pop_cnt; j++) {
        this->_private_cnt[j]       = 0;
        this->_n[j]                 = 0.0;
        this->_var_sites[j]         = 0.0;
        this->_num_indv_mean[j]     = 0.0;
        this->_num_indv_acc_mean[j] = 0.0;
        this->_num_indv_var[j]      = 0.0;
        this->_p_mean[j]            = 0.0;
        this->_p_acc_mean[j]        = 0.0;
        this->_p_var[j]             = 0.0;
        this->_obs_het_mean[j]      = 0.0;
        this->_obs_het_acc_mean[j]  = 0.0;
        this->_obs_het_var[j]       = 0.0;
        this->_obs_hom_mean[j]      = 0.0;
        this->_obs_hom_acc_mean[j]  = 0.0;
        this->_obs_hom_var[j]       = 0.0;
        this->_exp_het_mean[j]      = 0.0;
        this->_exp_het_acc_mean[j]  = 0.0;
        this->_exp_het_var[j]       = 0.0;
        this->_exp_hom_mean[j]      = 0.0;
        this->_exp_hom_acc_mean[j]  = 0.0;
        this->_exp_hom_var[j]       = 0.0;
        this->_pi_mean[j]           = 0.0;
        this->_pi_acc_mean[j]       = 0.0;
        this->_pi_var[j]            = 0.0;
        this->_fis_mean[j]          = 0.0;
        this->_fis_acc_mean[j]      = 0.0;
        this->_fis_var[j]           = 0.0;

        this->_n_all[j]                 = 0.0;
        this->_num_indv_mean_all[j]     = 0.0;
        this->_num_indv_acc_mean_all[j] = 0.0;
        this->_num_indv_var_all[j]      = 0.0;
        this->_p_mean_all[j]            = 0.0;
        this->_p_acc_mean_all[j]        = 0.0;
        this->_p_var_all[j]             = 0.0;
        this->_obs_het_mean_all[j]      = 0.0;
        this->_obs_het_acc_mean_all[j]  = 0.0;
        this->_obs_het_var_all[j]       = 0.0;
        this->_obs_hom_mean_all[j]      = 0.0;
        this->_obs_hom_acc_mean_all[j]  = 0.0;
        this->_obs_hom_var_all[j]       = 0.0;
        this->_exp_het_mean_all[j]      = 0.0;
        this->_exp_het_acc_mean_all[j]  = 0.0;
        this->_exp_het_var_all[j]       = 0.0;
        this->_exp_hom_mean_all[j]      = 0.0;
        this->_exp_hom_acc_mean_all[j]  = 0.0;
        this->_exp_hom_var_all[j]       = 0.0;
        this->_pi_mean_all[j]           = 0.0;
        this->_pi_acc_mean_all[j]       = 0.0;
        this->_pi_var_all[j]            = 0.0;
        this->_fis_mean_all[j]          = 0.0;
        this->_fis_acc_mean_all[j]      = 0.0;
        this->_fis_var_all[j]           = 0.0;
    }
}

int
SumStatsSummary::accumulate(const vector<LocBin *> &loci)
{
    //
    // We are calculating the mean, variance, and standard deviation for several variables.
    //   We will calculate them partially, for each set of loci input to the program using
    //   the algorithm described in:
    //     B. P. Welford. (1962) Note on a Method for Calculating Corrected Sums of Squares and
    //     Products. Technometrics: 4(3), pp. 419-420.
    //
    CSLocus        *cloc;
    const LocSum   *s;
    const LocTally *t;

    for (uint i = 0; i < loci.size(); i++) {
        cloc = loci[i]->cloc;
        t    = loci[i]->s->meta_pop();
        uint len = strlen(cloc->con);

        for (uint pos = 0; pos < len; pos++) {
            //
            // Compile private alleles
            //
            if (t->nucs[pos].priv_allele >= 0)
                _private_cnt[t->nucs[pos].priv_allele]++;

            if (t->nucs[pos].allele_cnt == 2) {

                for (uint pop = 0; pop < this->_pop_cnt; pop++) {

                    s = loci[i]->s->per_pop(pop);

                    if (s->nucs[pos].num_indv == 0) continue;

                    _n[pop]++;

                    if (s->nucs[pos].pi > 0) _var_sites[pop]++;

                    //
                    // Accumulate sums for each variable to calculate the means.
                    //
                    _num_indv_mean[pop] += s->nucs[pos].num_indv;
                    _p_mean[pop]        += s->nucs[pos].p;
                    _obs_het_mean[pop]  += s->nucs[pos].obs_het;
                    _obs_hom_mean[pop]  += s->nucs[pos].obs_hom;
                    _exp_het_mean[pop]  += s->nucs[pos].exp_het;
                    _exp_hom_mean[pop]  += s->nucs[pos].exp_hom;
                    _pi_mean[pop]       += s->nucs[pos].stat[0];
                    _fis_mean[pop]      += s->nucs[pos].stat[1] != -7.0 ? s->nucs[pos].stat[1] : 0.0;

                    _n_all[pop]++;
                    _num_indv_mean_all[pop] += s->nucs[pos].num_indv;
                    _p_mean_all[pop]        += s->nucs[pos].p;
                    _obs_het_mean_all[pop]  += s->nucs[pos].obs_het;
                    _obs_hom_mean_all[pop]  += s->nucs[pos].obs_hom;
                    _exp_het_mean_all[pop]  += s->nucs[pos].exp_het;
                    _exp_hom_mean_all[pop]  += s->nucs[pos].exp_hom;
                    _pi_mean_all[pop]       += s->nucs[pos].stat[0];
                    _fis_mean_all[pop]      += s->nucs[pos].stat[1] != -7.0 ? s->nucs[pos].stat[1] : 0.0;

                    //
                    // Accumulate a partial sum of squares to calculate the variance.
                    //
                    _num_indv_var[pop] += this->online_variance(s->nucs[pos].num_indv, _num_indv_acc_mean[pop], _n[pop]);
                    _p_var[pop]        += this->online_variance(s->nucs[pos].p,        _p_acc_mean[pop],        _n[pop]);
                    _obs_het_var[pop]  += this->online_variance(s->nucs[pos].obs_het,  _obs_het_acc_mean[pop],  _n[pop]);
                    _obs_hom_var[pop]  += this->online_variance(s->nucs[pos].obs_hom,  _obs_hom_acc_mean[pop],  _n[pop]);
                    _exp_het_var[pop]  += this->online_variance(s->nucs[pos].exp_het,  _exp_het_acc_mean[pop],  _n[pop]);
                    _exp_hom_var[pop]  += this->online_variance(s->nucs[pos].exp_hom,  _exp_hom_acc_mean[pop],  _n[pop]);
                    _pi_var[pop]       += this->online_variance(s->nucs[pos].stat[0],  _pi_acc_mean[pop],       _n[pop]);
                    _fis_var[pop]      += this->online_variance(s->nucs[pos].stat[1] != -7.0 ? s->nucs[pos].stat[1] : 0.0, _fis_acc_mean[pop], _n[pop]);

                    _num_indv_var_all[pop] += this->online_variance(s->nucs[pos].num_indv, _num_indv_acc_mean_all[pop], _n_all[pop]);
                    _p_var_all[pop]        += this->online_variance(s->nucs[pos].p,        _p_acc_mean_all[pop],        _n_all[pop]);
                    _obs_het_var_all[pop]  += this->online_variance(s->nucs[pos].obs_het,  _obs_het_acc_mean_all[pop],  _n_all[pop]);
                    _obs_hom_var_all[pop]  += this->online_variance(s->nucs[pos].obs_hom,  _obs_hom_acc_mean_all[pop],  _n_all[pop]);
                    _exp_het_var_all[pop]  += this->online_variance(s->nucs[pos].exp_het,  _exp_het_acc_mean_all[pop],  _n_all[pop]);
                    _exp_hom_var_all[pop]  += this->online_variance(s->nucs[pos].exp_hom,  _exp_hom_acc_mean_all[pop],  _n_all[pop]);
                    _pi_var_all[pop]       += this->online_variance(s->nucs[pos].stat[0],  _pi_acc_mean_all[pop],       _n_all[pop]);
                    _fis_var_all[pop]      += this->online_variance(s->nucs[pos].stat[1] != -7.0 ? s->nucs[pos].stat[1] : 0.0, _fis_acc_mean_all[pop], _n_all[pop]);
                }

            } else if (t->nucs[pos].allele_cnt == 1) {

                for (uint pop = 0; pop < this->_pop_cnt; pop++) {
                    s = loci[i]->s->per_pop(pop);

                    if (s->nucs[pos].num_indv == 0) continue;

                    _n_all[pop]++;
                    _num_indv_mean_all[pop] += s->nucs[pos].num_indv;
                    _p_mean_all[pop]        += s->nucs[pos].p;
                    _obs_het_mean_all[pop]  += s->nucs[pos].obs_het;
                    _obs_hom_mean_all[pop]  += s->nucs[pos].obs_hom;
                    _exp_het_mean_all[pop]  += s->nucs[pos].exp_het;
                    _exp_hom_mean_all[pop]  += s->nucs[pos].exp_hom;
                    _pi_mean_all[pop]       += s->nucs[pos].stat[0];
                    _fis_mean_all[pop]      += s->nucs[pos].stat[1] != -7.0 ? s->nucs[pos].stat[1] : 0.0;

                    _num_indv_var_all[pop] += this->online_variance(s->nucs[pos].num_indv, _num_indv_acc_mean_all[pop], _n_all[pop]);
                    _p_var_all[pop]        += this->online_variance(s->nucs[pos].p,        _p_acc_mean_all[pop],        _n_all[pop]);
                    _obs_het_var_all[pop]  += this->online_variance(s->nucs[pos].obs_het,  _obs_het_acc_mean_all[pop],  _n_all[pop]);
                    _obs_hom_var_all[pop]  += this->online_variance(s->nucs[pos].obs_hom,  _obs_hom_acc_mean_all[pop],  _n_all[pop]);
                    _exp_het_var_all[pop]  += this->online_variance(s->nucs[pos].exp_het,  _exp_het_acc_mean_all[pop],  _n_all[pop]);
                    _exp_hom_var_all[pop]  += this->online_variance(s->nucs[pos].exp_hom,  _exp_hom_acc_mean_all[pop],  _n_all[pop]);
                    _pi_var_all[pop]       += this->online_variance(s->nucs[pos].stat[0],  _pi_acc_mean_all[pop],       _n_all[pop]);
                    _fis_var_all[pop]      += this->online_variance(s->nucs[pos].stat[1] != -7.0 ? s->nucs[pos].stat[1] : 0.0, _fis_acc_mean_all[pop], _n_all[pop]);
                }
            }
        }
    }

    return 0;
}

inline double
SumStatsSummary::online_variance(double x, double &acc_mean, double n)
{
    double delta1, delta2;

    delta1    = x - acc_mean;
    acc_mean += delta1 / n;
    delta2    = x - acc_mean;
    return delta1 * delta2;
}

int
SumStatsSummary::final_calculation()
{
    //
    // Finish the mean calculation.
    //
    for (uint j = 0; j < this->_pop_cnt; j++) {
        _num_indv_mean[j] = _num_indv_mean[j] / _n[j];
        _p_mean[j]        = _p_mean[j]        / _n[j];
        _obs_het_mean[j]  = _obs_het_mean[j]  / _n[j];
        _obs_hom_mean[j]  = _obs_hom_mean[j]  / _n[j];
        _exp_het_mean[j]  = _exp_het_mean[j]  / _n[j];
        _exp_hom_mean[j]  = _exp_hom_mean[j]  / _n[j];
        _pi_mean[j]       = _pi_mean[j]       / _n[j];
        _fis_mean[j]      = _fis_mean[j]      / _n[j];

        _num_indv_mean_all[j] = _num_indv_mean_all[j] / _n_all[j];
        _p_mean_all[j]        = _p_mean_all[j]        / _n_all[j];
        _obs_het_mean_all[j]  = _obs_het_mean_all[j]  / _n_all[j];
        _obs_hom_mean_all[j]  = _obs_hom_mean_all[j]  / _n_all[j];
        _exp_het_mean_all[j]  = _exp_het_mean_all[j]  / _n_all[j];
        _exp_hom_mean_all[j]  = _exp_hom_mean_all[j]  / _n_all[j];
        _pi_mean_all[j]       = _pi_mean_all[j]       / _n_all[j];
        _fis_mean_all[j]      = _fis_mean_all[j]      / _n_all[j];
    }

    //
    // Finish the online variance calculation.
    //
    for (uint j = 0; j < this->_pop_cnt; j++) {
        _num_indv_var[j] = _num_indv_var[j] / (_n[j] - 1);
        _p_var[j]        = _p_var[j]        / (_n[j] - 1);
        _obs_het_var[j]  = _obs_het_var[j]  / (_n[j] - 1);
        _obs_hom_var[j]  = _obs_hom_var[j]  / (_n[j] - 1);
        _exp_het_var[j]  = _exp_het_var[j]  / (_n[j] - 1);
        _exp_hom_var[j]  = _exp_hom_var[j]  / (_n[j] - 1);
        _pi_var[j]       = _pi_var[j]       / (_n[j] - 1);
        _fis_var[j]      = _fis_var[j]      / (_n[j] - 1);

        _num_indv_var_all[j] = _num_indv_var_all[j] / (_n_all[j] - 1);
        _p_var_all[j]        = _p_var_all[j]        / (_n_all[j] - 1);
        _obs_het_var_all[j]  = _obs_het_var_all[j]  / (_n_all[j] - 1);
        _obs_hom_var_all[j]  = _obs_hom_var_all[j]  / (_n_all[j] - 1);
        _exp_het_var_all[j]  = _exp_het_var_all[j]  / (_n_all[j] - 1);
        _exp_hom_var_all[j]  = _exp_hom_var_all[j]  / (_n_all[j] - 1);
        _pi_var_all[j]       = _pi_var_all[j]       / (_n_all[j] - 1);
        _fis_var_all[j]      = _fis_var_all[j]      / (_n_all[j] - 1);
    }

    //
    // Calculate the first half of the standard deviation.
    //
    for (uint j = 0; j < this->_pop_cnt; j++) {
        _sq_n[j]     = sqrt(_n[j]);
        _sq_n_all[j] = sqrt(_n_all[j]);
    }

    return 0;
}

int
SumStatsSummary::write_results()
{
    string   path = out_path + out_prefix + ".sumstats_summary.tsv";
    ofstream fh(path.c_str(), ofstream::out);
    if (fh.fail()) {
        cerr << "Error opening sumstats summary file '" << path << "'\n";
        exit(1);
    }
    fh.precision(fieldw);
    fh.setf(std::ios::fixed);

    //
    // Write out summary statistics of the summary statistics.
    //
    fh << "# Variant positions\n"
       << "# Pop ID\t"
       << "Private\t"
       << "Num_Indv\t"
       << "Var\t"
       << "StdErr\t"
       << "P\t"
       << "Var\t"
       << "StdErr\t"
       << "Obs_Het\t"
       << "Var\t"
       << "StdErr\t"
       << "Obs_Hom\t"
       << "Var\t"
       << "StdErr\t"
       << "Exp_Het\t"
       << "Var\t"
       << "StdErr\t"
       << "Exp_Hom\t"
       << "Var\t"
       << "StdErr\t"
       << "Pi\t"
       << "Var\t"
       << "StdErr\t"
       << "Fis\t"
       << "Var\t"
       << "StdErr\n";

    for (uint j = 0; j < this->_pop_cnt; j++)
        fh << mpopi.pops()[j].name << "\t"
           << _private_cnt[j]         << "\t"
           << _num_indv_mean[j]       << "\t"
           << _num_indv_var[j]        << "\t"
           << sqrt(_num_indv_var[j]) / _sq_n[j] << "\t"
           << _p_mean[j]              << "\t"
           << _p_var[j]               << "\t"
           << sqrt(_p_var[j])         / _sq_n[j] << "\t"
           << _obs_het_mean[j]        << "\t"
           << _obs_het_var[j]         << "\t"
           << sqrt(_obs_het_var[j])   / _sq_n[j] << "\t"
           << _obs_hom_mean[j]        << "\t"
           << _obs_hom_var[j]         << "\t"
           << sqrt(_obs_hom_var[j])   / _sq_n[j] << "\t"
           << _exp_het_mean[j]        << "\t"
           << _exp_het_var[j]         << "\t"
           << sqrt(_exp_het_var[j])   / _sq_n[j] << "\t"
           << _exp_hom_mean[j]        << "\t"
           << _exp_hom_var[j]         << "\t"
           << sqrt(_exp_hom_var[j])   / _sq_n[j] << "\t"
           << _pi_mean[j]             << "\t"
           << _pi_var[j]              << "\t"
           << sqrt(_pi_var[j])        / _sq_n[j] << "\t"
           << _fis_mean[j]            << "\t"
           << _fis_var[j]             << "\t"
           << sqrt(_num_indv_var[j])  / _sq_n[j] << "\n";

    cerr << "\nPopulation summary statistics (more detail in populations.sumstats_summary.tsv):\n";

    for (uint j = 0; j < this->_pop_cnt; j++)
        cerr << "  " << mpopi.pops()[j].name << ": "
             << setprecision(fieldw) << _num_indv_mean[j] << " samples per locus; "
             << "pi: " << _pi_mean[j] << "; "
             << setprecision(10) << "all/variant/polymorphic sites: " << _n_all[j] << "/" << _n[j] << "/" << _var_sites[j] << "; "
             << setprecision(fieldw) << "private alleles: " << _private_cnt[j] << "\n";

    fh << "# All positions (variant and fixed)\n"
       << "# Pop ID\t"
       << "Private\t"
       << "Sites\t"
       << "Variant_Sites\t"
       << "Polymorphic_Sites\t"
       << "%Polymorphic_Loci\t"
       << "Num_Indv\t"
       << "Var\t"
       << "StdErr\t"
       << "P\t"
       << "Var\t"
       << "StdErr\t"
       << "Obs_Het\t"
       << "Var\t"
       << "StdErr\t"
       << "Obs_Hom\t"
       << "Var\t"
       << "StdErr\t"
       << "Exp_Het\t"
       << "Var\t"
       << "StdErr\t"
       << "Exp_Hom\t"
       << "Var\t"
       << "StdErr\t"
       << "Pi\t"
       << "Var\t"
       << "StdErr\t"
       << "Fis\t"
       << "Var\t"
       << "StdErr\n";

    for (uint j = 0; j < this->_pop_cnt; j++) {
        fh << mpopi.pops()[j].name << "\t"
           << _private_cnt[j]             << "\t" << setprecision(0)
           << _n_all[j]                   << "\t"
           << _n[j]                       << "\t"
           << _var_sites[j]               << "\t" << setprecision(fieldw)
           << _var_sites[j]             / _n_all[j] * 100 << "\t"
           << _num_indv_mean_all[j]       << "\t"
           << _num_indv_var_all[j]        << "\t"
           << sqrt(_num_indv_var_all[j]) / _sq_n_all[j] << "\t"
           << _p_mean_all[j]              << "\t"
           << _p_var_all[j]               << "\t"
           << sqrt(_p_var_all[j])        / _sq_n_all[j] << "\t"
           << _obs_het_mean_all[j]        << "\t"
           << _obs_het_var_all[j]         << "\t"
           << sqrt(_obs_het_var_all[j])  / _sq_n_all[j] << "\t"
           << _obs_hom_mean_all[j]        << "\t"
           << _obs_hom_var_all[j]         << "\t"
           << sqrt(_obs_hom_var_all[j])  / _sq_n_all[j] << "\t"
           << _exp_het_mean_all[j]        << "\t"
           << _exp_het_var_all[j]         << "\t"
           << sqrt(_exp_het_var_all[j])  / _sq_n_all[j] << "\t"
           << _exp_hom_mean_all[j]        << "\t"
           << _exp_hom_var_all[j]         << "\t"
           << sqrt(_exp_hom_var_all[j])  / _sq_n_all[j] << "\t"
           << _pi_mean_all[j]             << "\t"
           << _pi_var_all[j]              << "\t"
           << sqrt(_pi_var_all[j])       / _sq_n_all[j] << "\t"
           << _fis_mean_all[j]            << "\t"
           << _fis_var_all[j]             << "\t"
           << sqrt(_num_indv_var_all[j]) / _sq_n_all[j] << "\n";
    }

    fh.close();

    return 0;
}

SumStatsSummary::~SumStatsSummary()
{
    delete [] this->_private_cnt;
    delete [] this->_n;
    delete [] this->_var_sites;
    delete [] this->_sq_n;
    delete [] this->_num_indv_mean;
    delete [] this->_num_indv_acc_mean;
    delete [] this->_num_indv_var;
    delete [] this->_p_mean;
    delete [] this->_p_acc_mean;
    delete [] this->_p_var;
    delete [] this->_obs_het_mean;
    delete [] this->_obs_het_acc_mean;
    delete [] this->_obs_het_var;
    delete [] this->_obs_hom_mean;
    delete [] this->_obs_hom_acc_mean;
    delete [] this->_obs_hom_var;
    delete [] this->_exp_het_mean;
    delete [] this->_exp_het_acc_mean;
    delete [] this->_exp_het_var;
    delete [] this->_exp_hom_mean;
    delete [] this->_exp_hom_acc_mean;
    delete [] this->_exp_hom_var;
    delete [] this->_pi_mean;
    delete [] this->_pi_acc_mean;
    delete [] this->_pi_var;
    delete [] this->_fis_mean;
    delete [] this->_fis_acc_mean;
    delete [] this->_fis_var;

    delete [] this->_n_all;
    delete [] this->_sq_n_all;
    delete [] this->_num_indv_mean_all;
    delete [] this->_num_indv_acc_mean_all;
    delete [] this->_num_indv_var_all;
    delete [] this->_p_mean_all;
    delete [] this->_p_acc_mean_all;
    delete [] this->_p_var_all;
    delete [] this->_obs_het_mean_all;
    delete [] this->_obs_het_acc_mean_all;
    delete [] this->_obs_het_var_all;
    delete [] this->_obs_hom_mean_all;
    delete [] this->_obs_hom_acc_mean_all;
    delete [] this->_obs_hom_var_all;
    delete [] this->_exp_het_mean_all;
    delete [] this->_exp_het_acc_mean_all;
    delete [] this->_exp_het_var_all;
    delete [] this->_exp_hom_mean_all;
    delete [] this->_exp_hom_acc_mean_all;
    delete [] this->_exp_hom_var_all;
    delete [] this->_pi_mean_all;
    delete [] this->_pi_acc_mean_all;
    delete [] this->_pi_var_all;
    delete [] this->_fis_mean_all;
    delete [] this->_fis_acc_mean_all;
    delete [] this->_fis_var_all;
}

int
correct_fst_bonferroni_win(vector<PopPair *> &pairs)
{
    int      limit = 3 * sigma;
    int      limit_l, limit_u;
    double   correction;
    uint     cnt, pos_l, pos_u;

    pos_l = 0;
    pos_u = 0;

    for (uint pos_c = 0; pos_c < pairs.size(); pos_c++) {
        if (pairs[pos_c] == NULL) continue;

        limit_l = pairs[pos_c]->bp - limit > 0 ? pairs[pos_c]->bp - limit : 0;
        limit_u = pairs[pos_c]->bp + limit;

        while (pos_l <  pairs.size()) {
            if (pairs[pos_l] == NULL) {
                pos_l++;
            } else {
                if (pairs[pos_l]->bp < limit_l)
                    pos_l++;
                else
                    break;
            }
        }
        while (pos_u < pairs.size()) {
            if (pairs[pos_u] == NULL) {
                pos_u++;
            } else {
                if (pairs[pos_u]->bp < limit_u)
                    pos_u++;
                else
                    break;
            }
        }

        cnt = 0;
        for (uint i = pos_l; i < pos_u; i++) {
            if (pairs[i] == NULL) continue;
            cnt++;
        }

        correction = p_value_cutoff / cnt;
        pairs[pos_c]->stat[0] = pairs[pos_c]->fet_p < correction ? pairs[pos_c]->fst : 0;
    }

    return 0;
}

int
LocusSmoothing::snpstats(const vector<LocBin *> &loci, ofstream &log_fh)
{
    for (uint i = 0; i < this->_mpopi->pops().size(); i++) {

        vector<const SumStat *> sites;

        this->_ord_ss->order(sites, loci, i);
        this->_ks_ss->smooth(sites);
    }

    return 0;
}

int
LocusSmoothing::hapstats(const vector<LocBin *> &loci, ofstream &log_fh)
{

    for (uint i = 0; i < this->_mpopi->pops().size(); i++) {
        vector<const LocStat *> sites;

        this->_ord_ls->order(sites, loci);
        this->_ks_ls->smooth(sites);
    }

    return 0;
}

int
LocusSmoothing::snp_divergence(const vector<LocBin *> &loci, const vector<vector<PopPair **>> &div, ofstream &log_fh)
{
    for (uint i = 0; i < div.size(); i++) {
        assert(div[i].size() == loci.size());

        vector<const PopPair *> sites;

        if (this->_ord_pp->order(sites, loci, div[i]))
            this->_ks_pp->smooth(sites);
    }

    return 0;
}

int
LocusSmoothing::hap_divergence(const vector<LocBin *> &loci,
                               const vector<vector<HapStat *>> &div,
                               const vector<HapStat *> &metadiv,
                               ofstream &log_fh)
{
    for (uint i = 0; i < div.size(); i++) {
        assert(div[i].size() == loci.size());

        vector<const HapStat *> sites;

        this->_ord_hs->order(sites, loci, div[i]);
        this->_ks_hs->smooth(sites);
    }

    //
    // Kernel-smooth the haplotype divergence statistics for the metapopulation.
    //
    vector<const HapStat *> sites;

    this->_ord_hs->order(sites, loci, metadiv);
    this->_ks_hs->smooth(sites);

    return 0;
}

int
bootstrap_popstats_approximate_dist(vector<double> &fis_samples,
                                    vector<double> &pi_samples,
                                    vector<int>  &allele_samples,
                                    double *weights, int *snp_dist, int sites_per_snp,
                                    map<int, vector<double> > &approx_fis_dist,
                                    map<int, vector<double> > &approx_pi_dist)
{
    //
    // Allocate an array of bootstrap resampling objects.
    //
    int win_size = 6 * sigma + 1;
    int win_cntr = win_size / 2;

    //
    // Initialize the Fst distribution map.
    //
    for (int i = 0; i < max_snp_dist; i++) {
        if (snp_dist[i] == 0.0) continue;

        // cerr << "SNP Dist: " << i << " snps occurred " << snp_dist[i] << "\n";
        approx_fis_dist[i] = vector<double> ();
        approx_fis_dist[i].reserve(bootstrap_reps);

        approx_pi_dist[i] = vector<double> ();
        approx_pi_dist[i].reserve(bootstrap_reps);
    }

    vector<int> poss;
    poss.reserve(max_snp_dist);
    double weighted_fis, weighted_pi, sum_fis, sum_pi, final_weight_fis, final_weight_pi;
    // int    index_1, index_2;
    int    pos, index_3, dist, start, end;
    int    half = sites_per_snp / 2;

    for (int i = 0; i < max_snp_dist; i++) {
        if (snp_dist[i] == 0.0) continue;

        cerr << "  Generating NULL distribution for " << i << " SNPs...\n";

        // #pragma omp parallel private(poss, pos, index_1, index_2, index_3, dist, sum_fis, sum_pi, weighted_fis, weighted_pi, final_weight_fis, final_weight_pi)
        #pragma omp parallel private(poss, pos, index_3, dist, sum_fis, sum_pi, weighted_fis, weighted_pi, final_weight_fis, final_weight_pi)
        {
            BSample *bs  = new BSample[win_size];

            //
            // Populate the BSample objects.
            //
            for (int n = 0; n < win_size;  n++)
                bs[n].bp = n + 1;

            vector<double> fiss, pis;

            //
            // Bootstrap this bitch.
            //
            #pragma omp for schedule(dynamic, 1)
            for (int j = 0; j < bootstrap_reps; j++) {
                // cerr << "    Bootsrap rep " << j << "\n";

                //
                // First SNP is always placed at the center of the window.
                //
                pos     = win_cntr;
                // index_1 = (int) (fis_samples.size()    * (random() / (RAND_MAX + 1.0)));
                // index_2 = (int) (pi_samples.size()     * (random() / (RAND_MAX + 1.0)));
                index_3 = (int) (allele_samples.size() * (random() / (RAND_MAX + 1.0)));
                //
                // Fill in the area around the SNP with fixed sites.
                //
                start = pos - half > 0 ? pos - half : 0;
                end   = pos + half < win_size ? pos + half : win_size;
                for (int n = start; n < end; n++) {
                    // bs[n].f       = 0;
                    // bs[n].pi      = 0;
                    bs[n].alleles = bs[pos].alleles;
                    poss.push_back(n);
                }
                // bs[pos].f       = fis_samples[index_1];
                // bs[pos].pi      = pi_samples[index_2];
                bs[pos].alleles = allele_samples[index_3];
                // cerr << "      Placing SNP at position: " << pos << "; with data from " << index_1 << " filling area from " << start << " to " << end << "\n";

                //
                // Randomly select the positions and values for each SNP to populate the window
                //
                for (int k = 0; k < i - 1; k++) {
                    pos     = (int) (win_size * (random() / (RAND_MAX + 1.0)));
                    // index_1 = (int) (fis_samples.size()    * (random() / (RAND_MAX + 1.0)));
                    // index_2 = (int) (pi_samples.size()     * (random() / (RAND_MAX + 1.0)));
                    index_3 = (int) (allele_samples.size() * (random() / (RAND_MAX + 1.0)));

                    poss.push_back(pos);
                    //
                    // Fill in the area around the SNP with fixed sites.
                    //
                    start = pos - half > 0 ? pos - half : 0;
                    end   = pos + half < win_size ? pos + half : win_size;
                    for (int n = start; n < end; n++) {
                        // bs[n].f       = 0;
                        // bs[n].pi      = 0;
                        bs[n].alleles = bs[pos].alleles;
                        poss.push_back(n);
                    }
                    // bs[pos].f       = fis_samples[index_1];
                    // bs[pos].pi      = pi_samples[index_2];
                    bs[pos].alleles = allele_samples[index_3];
                    // cerr << "      Placing SNP at position: " << pos << "; with data from " << index_1 << " filling area from " << start << " to " << end << "\n";
                }

                weighted_fis = 0.0;
                sum_fis      = 0.0;
                weighted_pi  = 0.0;
                sum_pi       = 0.0;

                for (int n = 0; n < win_size; n++) {
                    // if (bs[n].pi < 0.0)
                    // continue;
                    //
                    // Calculate weighted Fst at this position.
                    //
                    dist = bs[n].bp > bs[win_cntr].bp ? bs[n].bp - bs[win_cntr].bp : bs[win_cntr].bp - bs[n].bp;

                    final_weight_fis = (bs[n].alleles - 1) * weights[dist];
                    // weighted_fis    += bs[n].f * final_weight_fis;
                    sum_fis         += final_weight_fis;

                    final_weight_pi  = (bs[n].alleles - 1) * weights[dist];
                    // weighted_pi     += bs[n].pi * final_weight_pi;
                    sum_pi          += final_weight_pi;
                }

                fiss.push_back(weighted_fis / sum_fis);
                pis.push_back(weighted_pi  / sum_pi);
                // cerr << "      New weighted fis value: " << weighted_fis / sum_fis << "; size: " << fiss.size() << "\n";

                for (uint n = 0; n < poss.size(); n++) {
                    // bs[poss[n]].f  = 0.0;
                    // bs[poss[n]].pi = -1.0;
                }
                poss.clear();
            }

//          #pragma omp critical
//          {
//              vector<double> &f = approx_fis_dist[i];
//              for (uint n = 0; n < fiss.size(); n++)
//                  f.push_back(fiss[n]);
//              vector<double> &p = approx_pi_dist[i];
//              for (uint n = 0; n < pis.size(); n++)
//                  p.push_back(pis[n]);
//          }

            delete [] bs;
        }

        sort(approx_fis_dist[i].begin(), approx_fis_dist[i].end());
        sort(approx_pi_dist[i].begin(),  approx_pi_dist[i].end());
    }

    return 0;
}

int
bootstrap_fst_approximate_dist(vector<double> &fst_samples,
                               vector<int>  &allele_samples,
                               double *weights, int *snp_dist,
                               map<int, vector<double> > &approx_fst_dist)
{
    //
    // Allocate an array of bootstrap resampling objects.
    //
    int win_size = 6 * sigma + 1;
    int win_cntr = win_size / 2;

    //
    // Initialize the Fst distribution map.
    //
    for (int i = 0; i < max_snp_dist; i++) {
        if (snp_dist[i] == 0.0) continue;

        // cerr << "SNP Dist: " << i << " snps occurred " << snp_dist[i] << "\n";
        approx_fst_dist[i] = vector<double> ();
        approx_fst_dist[i].reserve(bootstrap_reps);
    }

    vector<int> poss;
    poss.reserve(max_snp_dist);
    double weighted_fst, sum, final_weight;
    //int    index_1;
    int    pos, index_2, dist;

    for (int i = 0; i < max_snp_dist; i++) {
        if (snp_dist[i] == 0.0) continue;

        cerr << "  Generating NULL distribution for " << i << " SNPs...\n";

        // #pragma omp parallel private(poss, pos, index_1, index_2, dist, sum, weighted_fst, final_weight)
        #pragma omp parallel private(poss, pos, index_2, dist, sum, weighted_fst, final_weight)
        {
            BSample *bs  = new BSample[win_size];

            //
            // Populate the BSample objects.
            //
            for (int n = 0; n < win_size;  n++)
                bs[n].bp = n + 1;

            vector<double> fsts;

            //
            // Bootstrap this bitch.
            //
            #pragma omp for schedule(dynamic, 1)
            for (int j = 0; j < bootstrap_reps; j++) {
                // cerr << "Bootsrap rep " << j << "\n";

                //
                // First SNP is always placed at the center of the window.
                //
                pos     = win_cntr;
                // index_1 = (int) (fst_samples.size() * (random() / (RAND_MAX + 1.0)));
                index_2 = (int) (allele_samples.size() * (random() / (RAND_MAX + 1.0)));
                // bs[pos].f       = fst_samples[index_1];
                bs[pos].alleles = allele_samples[index_2];

                //
                // Randomly select the positions and values for each SNP to populate the window
                //
                for (int k = 0; k < i - 1; k++) {
                    pos     = (int) (win_size * (random() / (RAND_MAX + 1.0)));
                    // index_1 = (int) (fst_samples.size() * (random() / (RAND_MAX + 1.0)));
                    index_2 = (int) (allele_samples.size() * (random() / (RAND_MAX + 1.0)));
                    // bs[pos].f       = fst_samples[index_1];
                    // bs[pos].alleles = allele_samples[index_2];
                    // cerr << "  " << j << ": Placing SNP at position: " << pos << " with data from index " << index_1 << "\n";

                    poss.push_back(pos);
                }

                weighted_fst = 0.0;
                sum          = 0.0;

                for (int n = 0; n < win_size; n++) {
                    // if (bs[n].f == 0.0)
                    // continue;
                    //
                    // Calculate weighted Fst at this position.
                    //
                    dist = bs[n].bp > bs[win_cntr].bp ? bs[n].bp - bs[win_cntr].bp : bs[win_cntr].bp - bs[n].bp;

                    final_weight  = (bs[n].alleles - 1) * weights[dist];
                    // weighted_fst += bs[n].f * final_weight;
                    sum          += final_weight;
                }

                fsts.push_back(weighted_fst / sum);
                // cerr << "    New weighted Fst value: " << weighted_fst / sum << "; size: " << fsts.size() << "\n";

                // for (uint n = 0; n < poss.size(); n++)
                // bs[poss[n]].f = 0.0;
                poss.clear();
            }

//          #pragma omp critical
//          {
//              vector<double> &f = approx_fst_dist[i];
//              for (uint n = 0; n < fsts.size(); n++)
//                  f.push_back(fsts[n]);
//          }

            delete [] bs;
        }

        sort(approx_fst_dist[i].begin(), approx_fst_dist[i].end());
    }

    return 0;
}

double
bootstrap_approximate_pval(int snp_cnt, double stat, map<int, vector<double> > &approx_dist)
{
    if (approx_dist.count(snp_cnt) == 0)
        return 1.0;

    vector<double>::iterator up;
    vector<double> &dist = approx_dist[snp_cnt];
    double pos;

    up  = upper_bound(dist.begin(), dist.end(), stat);

    if (up == dist.begin())
        pos = 1;
    else if (up == dist.end())
        pos = dist.size();
    else
        pos = up - dist.begin() + 1;

    double res = 1.0 - (pos / (double) dist.size());

    // cerr << "Generated Approx Smoothed Fst Distribution:\n";
    // for (uint n = 0; n < dist.size(); n++)
    //  cerr << "  n: " << n << "; Fst: " << dist[n] << "\n";

    // cerr << "Comparing Fst value: " << stat
    //   << " at position " << (up - dist.begin()) << " out of "
    //   << dist.size() << " positions (converted position: " << pos << "); pvalue: " << res << ".\n";

    return res;
}

int
LocusFilter::load_blacklist(string path)
{
    char     line[id_len];
    ifstream fh(path.c_str(), ifstream::in);

    if (fh.fail()) {
        cerr << "Error opening white/black list file '" << path << "'\n";
        exit(1);
    }

    size_t line_num = 0;
    while (fh.getline(line, id_len)) {
        ++line_num;

        //
        // Skip blank & commented lines ; correct windows-style line ends.
        //
        size_t len = strlen(line);
        if (len == 0) {
            continue;
        } else if (line[len-1] == '\r') {
            line[len-1] = '\0';
            --len;
            if (len == 0)
                continue;
        }
        char* p = line;
        while (isspace(*p) && *p != '\0')
            ++p;
        if (*p == '#')
            continue;

        //
        // Parse the blacklist
        //
        char* e;
        int marker = (int) strtol(line, &e, 10);
        if (*e == '\0') {
            this->_blacklist.insert(marker);
        } else {
            cerr << "Error: Unable to parse blacklist '" << path << "' at line " << line_num << ".\n";
            throw exception();
        }
    }

    fh.close();

    if (this->_blacklist.size() == 0) {
        cerr << "Unable to load any markers from '" << path << "'\n";
        exit(1);
    }

    return (int) this->_blacklist.size();
}

int load_marker_list(string path, set<int> &list) {
    char     line[id_len];
    ifstream fh(path.c_str(), ifstream::in);

    if (fh.fail()) {
        cerr << "Error opening white/black list file '" << path << "'\n";
        exit(1);
    }

    size_t line_num = 0;
    while (fh.getline(line, id_len)) {
        ++line_num;

        //
        // Skip blank & commented lines ; correct windows-style line ends.
        //
        size_t len = strlen(line);
        if (len == 0) {
            continue;
        } else if (line[len-1] == '\r') {
            line[len-1] = '\0';
            --len;
            if (len == 0)
                continue;
        }
        char* p = line;
        while (isspace(*p) && *p != '\0')
            ++p;
        if (*p == '#')
            continue;

        //
        // Parse the blacklist
        //
        char* e;
        int marker = (int) strtol(line, &e, 10);
        if (*e == '\0') {
            list.insert(marker);
        } else {
            cerr << "Error: Unable to parse blacklist '" << path << "' at line " << line_num << ".\n";
            throw exception();
        }
    }

    fh.close();

    if (list.size() == 0) {
        cerr << "Unable to load any markers from '" << path << "'\n";
        exit(1);
    }

    return 0;
}

int
LocusFilter::load_whitelist(string path)
{
    char     line[id_len];
    ifstream fh(path.c_str(), ifstream::in);

    if (fh.fail()) {
        cerr << "Error opening white/black list file '" << path << "'\n";
        exit(1);
    }

    vector<string> parts;
    uint col;
    char *e;

    uint line_num = 1;
    while (fh.getline(line, id_len)) {

        //
        // Skip blank & commented lines ; correct windows-style line ends.
        //
        size_t len = strlen(line);
        if (len == 0) {
            continue;
        } else if (line[len-1] == '\r') {
            line[len-1] = '\0';
            --len;
            if (len == 0)
                continue;
        }
        char* p = line;
        while (isspace(*p) && *p != '\0')
            ++p;
        if (*p == '#')
            continue;

        //
        // Parse the whitelist, we expect:
        // <marker>[<tab><snp column>]
        //
        parse_tsv(line, parts);

        if (parts.size() > 2) {
            cerr << "Too many columns in whitelist " << path << "' at line " << line_num << "\n";
            exit(1);

        } else if (parts.size() == 2) {
            int marker = (int) strtol(parts[0].c_str(), &e, 10);
            if (*e != '\0') {
                cerr << "Unable to parse whitelist, '" << path << "' at line " << line_num << "\n";
                exit(1);
            }
            col = (int) strtol(parts[1].c_str(), &e, 10);
            if (*e != '\0') {
                cerr << "Unable to parse whitelist, '" << path << "' at line " << line_num << "\n";
                exit(1);
            }
            this->_whitelist[marker].insert(col);

        } else {
            int marker = (int) strtol(parts[0].c_str(), &e, 10);
            if (*e != '\0') {
                cerr << "Unable to parse whitelist, '" << path << "' at line " << line_num << "\n";
                exit(1);
            }
            this->_whitelist.insert(make_pair(marker, set<int>()));
        }

        line_num++;
    }

    fh.close();

    if (this->_whitelist.size() == 0) {
        cerr << "Unable to load any markers from '" << path << "'\n";
        help();
    }

    return (int) this->_whitelist.size();
}

int load_marker_column_list(string path, map<int, set<int> > &list) {
    char     line[id_len];
    ifstream fh(path.c_str(), ifstream::in);

    if (fh.fail()) {
        cerr << "Error opening white/black list file '" << path << "'\n";
        exit(1);
    }

    vector<string> parts;
    uint col;
    char *e;

    uint line_num = 1;
    while (fh.getline(line, id_len)) {

        //
        // Skip blank & commented lines ; correct windows-style line ends.
        //
        size_t len = strlen(line);
        if (len == 0) {
            continue;
        } else if (line[len-1] == '\r') {
            line[len-1] = '\0';
            --len;
            if (len == 0)
                continue;
        }
        char* p = line;
        while (isspace(*p) && *p != '\0')
            ++p;
        if (*p == '#')
            continue;

        //
        // Parse the whitelist, we expect:
        // <marker>[<tab><snp column>]
        //
        parse_tsv(line, parts);

        if (parts.size() > 2) {
            cerr << "Too many columns in whitelist " << path << "' at line " << line_num << "\n";
            exit(1);

        } else if (parts.size() == 2) {
            int marker = (int) strtol(parts[0].c_str(), &e, 10);
            if (*e != '\0') {
                cerr << "Unable to parse whitelist, '" << path << "' at line " << line_num << "\n";
                exit(1);
            }
            col = (int) strtol(parts[1].c_str(), &e, 10);
            if (*e != '\0') {
                cerr << "Unable to parse whitelist, '" << path << "' at line " << line_num << "\n";
                exit(1);
            }
            list[marker].insert(col);

        } else {
            int marker = (int) strtol(parts[0].c_str(), &e, 10);
            if (*e != '\0') {
                cerr << "Unable to parse whitelist, '" << path << "' at line " << line_num << "\n";
                exit(1);
            }
            list.insert(make_pair(marker, set<int>()));
        }

        line_num++;
    }

    fh.close();

    if (list.size() == 0) {
        cerr << "Unable to load any markers from '" << path << "'\n";
        help();
    }

    return 0;
}

bool
hap_compare(const pair<string, int>& a, const pair<string, int>& b)
{
    return (a.second > b.second);
}

void
output_parameters(ostream &fh)
{
    fh << "populations parameters selected:\n";
    if (input_mode == InputMode::vcf)
        fh << "  Input mode: VCF\n";
    fh
        << "  Percent samples limit per population: " << sample_limit << "\n"
        << "  Locus Population limit: " << population_limit << "\n"
        << "  Minimum stack depth: " << min_stack_depth << "\n"
        << "  Log liklihood filtering: " << (filter_lnl == true ? "on"  : "off") << "; threshold: " << lnl_limit << "\n"
        << "  Minor allele frequency cutoff: " << minor_allele_freq << "\n"
        << "  Maximum observed heterozygosity cutoff: " << max_obs_het << "\n"
        << "  Applying Fst correction: ";
    switch(fst_correction) {
    case p_value:
        fh << "P-value correction.\n";
        break;
    case bonferroni_win:
        fh << "Bonferroni correction within sliding window.\n";
        break;
    case bonferroni_gen:
        fh << "Bonferroni correction across genome wide sites.\n";
        break;
    case no_correction:
        fh << "none.\n";
        break;
    }
    fh << "  Pi/Fis kernel smoothing: " << (smooth_popstats == true ? "on" : "off") << "\n"
       << "  Fstats kernel smoothing: " << (smooth_fstats   == true ? "on" : "off") << "\n"
       << "  Bootstrap resampling: ";
    if (bootstrap)
        fh << "on, " << (bootstrap_type == bs_exact ? "exact; " : "approximate; ") << bootstrap_reps << " reptitions\n";
    else
        fh << "off\n";

    fh << "\n";
}

int
parse_command_line(int argc, char* argv[])
{
    Export *exp;

    while (1) {
        static struct option long_options[] = {
            {"help",           no_argument,       NULL, 'h'},
            {"version",        no_argument,       NULL, 'v'},
            {"verbose",        no_argument,       NULL, 'd'},
            {"vcf",            no_argument,       NULL, 1004},
            {"vcf_haplotypes", no_argument,       NULL, 'n'},
            {"fasta_loci",     no_argument,       NULL, 1006},
            {"fasta_samples",  no_argument,       NULL, 'J'}, {"fasta_strict", no_argument, NULL, 'J'},
            {"fasta_samples_raw", no_argument,    NULL, 'F'}, {"fasta", no_argument, NULL, 'F'},
            {"structure",      no_argument,       NULL, 'S'},
            {"fastphase",      no_argument,       NULL, 'A'},
            {"phase",          no_argument,       NULL, 'C'},
            {"beagle",         no_argument,       NULL, 'E'},
            {"beagle_phased",  no_argument,       NULL, 'H'},
            {"plink",          no_argument,       NULL, 'K'},
            {"genomic",        no_argument,       NULL, 'g'},
            {"genepop",        no_argument,       NULL, 'G'},
            {"phylip",         no_argument,       NULL, 'Y'},
            {"phylip_var",     no_argument,       NULL, 'L'},
            {"phylip_var_all", no_argument,       NULL, 'T'},
            {"hzar",           no_argument,       NULL, 'Z'},
            {"treemix",        no_argument,       NULL, 'U'},
            {"merge_sites",    no_argument,       NULL, 'D'},
            {"sigma",          required_argument, NULL, 1005},
            {"threads",        required_argument, NULL, 't'},
            {"in_path",        required_argument, NULL, 'P'},
            {"v1",             no_argument,       NULL, 2000},
            {"out_path",       required_argument, NULL, 'O'},
            {"in_vcf",         required_argument, NULL, 'V'},
            {"progeny",        required_argument, NULL, 'r'},
            {"min_depth",      required_argument, NULL, 'm'},
            {"renz",           required_argument, NULL, 'e'},
            {"popmap",         required_argument, NULL, 'M'},
            {"whitelist",      required_argument, NULL, 'W'},
            {"blacklist",      required_argument, NULL, 'B'},
            {"batch_size",     required_argument, NULL, 1999},
            {"write_single_snp",  no_argument,       NULL, 'I'},
            {"write_random_snp",  no_argument,       NULL, 'j'},
            {"ordered_export",    no_argument,       NULL, 1002},
            {"smooth",            no_argument,       NULL, 'k'},
            {"smooth_fstats",     no_argument,       NULL, 1007},
            {"smooth_popstats",   no_argument,       NULL, 1008},
            {"fstats",            no_argument,       NULL, '6'},
            {"log_fst_comp",      no_argument,       NULL, 'l'},
            {"bootstrap_type",    required_argument, NULL, 1001},
            {"bootstrap_reps",    required_argument, NULL, 1003},
            {"bootstrap_wl",      required_argument, NULL, 'Q'},
            {"bootstrap",         no_argument,       NULL, '1'},
            {"bootstrap_fst",     no_argument,       NULL, '2'},
            {"bootstrap_phist",   no_argument,       NULL, '3'},
            {"bootstrap_div",     no_argument,       NULL, '4'},
            {"bootstrap_pifis",   no_argument,       NULL, '5'},
            {"min_populations",   required_argument, NULL, 'p'},
            {"min_maf",           required_argument, NULL, 'a'},
            {"max_obs_het",       required_argument, NULL, 'q'},
            {"lnl_lim",           required_argument, NULL, 'c'},
            {"merge_prune_lim",   required_argument, NULL, 'i'},
            {"fst_correction",    required_argument, NULL, 'f'},
            {"p_value_cutoff",    required_argument, NULL, 'u'},
            {"debug_flags",       required_argument, NULL, 1000},
            {0, 0, 0, 0}
        };

        // getopt_long stores the option index here.
        int c = getopt_long(argc, argv, "ACDEFGHJKLNSTUV:YZ123456dghjklnva:c:e:f:i:m:o:p:q:r:t:u:w:B:I:M:O:P:R:Q:W:", long_options, NULL);

        // Detect the end of the options.
        if (c == -1)
            break;

        switch (c) {
        case 'h':
            help();
            break;
        case 'd':
            verbose = true;
            break;
        case 't':
            num_threads = atoi(optarg);
            break;
        case 'P':
            in_path = optarg;
            if (!in_path.empty() && in_path.back() != '/')
                in_path += "/";
            break;
        case 2000: //v1
            input_mode = InputMode::stacks;
            break;
        case 1999:
            batch_size = is_integer(optarg);
            break;
        case 'O':
            out_path = optarg;
            if (!out_path.empty() && out_path.back() != '/')
                out_path += "/";
            break;
        case 'V':
            in_vcf_path = optarg;
            break;
        case 'M':
            pmap_path = optarg;
            break;
        case 'D':
            merge_sites = true;
            break;
        case 'i':
            merge_prune_lim = is_double(optarg);
            if (merge_prune_lim > 1.0)
                merge_prune_lim = merge_prune_lim / 100;

            if (merge_prune_lim < 0 || merge_prune_lim > 1.0) {
                cerr << "Unable to parse the merge sites pruning limit.\n";
                help();
            }
            break;
        case 'q':
            max_obs_het = is_double(optarg);
            if (max_obs_het > 1)
                max_obs_het = max_obs_het / 100;

            if (max_obs_het < 0 || max_obs_het > 1.0) {
                cerr << "Unable to parse the maximum observed heterozygosity.\n";
                help();
            }
            break;
        case 'r':
            sample_limit = atof(optarg);
            if (sample_limit > 1)
                sample_limit = sample_limit / 100;

            if (sample_limit > 1.0) {
                cerr << "Unable to parse the sample limit frequency\n";
                help();
            }
            break;
        case 'p':
            population_limit = atoi(optarg);
            break;
        case 'k':
            smooth_popstats = true;
            smooth_fstats   = true;
            calc_fstats     = true;
            break;
        case 1007:
            smooth_fstats   = true;
            calc_fstats     = true;
            break;
        case 1008:
            smooth_popstats = true;
            break;
        case '6':
            calc_fstats = true;
            break;
        case 'l':
            log_fst_comp = true;
            break;
        case '1':
            bootstrap       = true;
            bootstrap_fst   = true;
            bootstrap_phist = true;
            bootstrap_pifis = true;
            bootstrap_div   = true;
            break;
        case '2':
            bootstrap_fst   = true;
            break;
        case '3':
            bootstrap_phist = true;
            break;
        case '4':
            bootstrap_div = true;
            break;
        case '5':
            bootstrap_pifis = true;
            break;
        case 1001:
            if (strcasecmp(optarg, "exact") == 0)
                bootstrap_type = bs_exact;
            else if (strcasecmp(optarg, "approx") == 0)
                bootstrap_type = bs_approx;
            else {
                cerr << "Unknown bootstrap type specified '" << optarg << "'\n";
                help();
            }
            break;
        case 1003:
            bootstrap_reps = atoi(optarg);
            break;
        case 'Q':
            bs_wl_file = optarg;
            bootstrap_wl = true;
            break;
        case 'c':
            lnl_limit  = is_double(optarg);
            filter_lnl = true;
            break;
        case 'I':
            write_single_snp = true;
            break;
        case 'j':
            write_random_snp = true;
            break;
        case 1002:
            ordered_export = true;
            break;
        case 1004:
            exports.push_back(new VcfExport());
            break;
        case 'n':
            vcf_haplo_out = true;
            break;
        case 1006:
            exp = new FastaLociExport();
            exports.push_back(exp);
            break;
        case 'F':
            exp = new FastaRawExport();
            exports.push_back(exp);
            break;
        case 'J':
            exp = new FastaSamplesExport();
            exports.push_back(exp);
            break;
        case 'G':
            exp = new GenePopExport();
            exports.push_back(exp);
            break;
        case 'S':
            structure_out = true;
            break;
        case 'A':
            fastphase_out = true;
            break;
        case 'C':
            phase_out = true;
            break;
        case 'E':
            beagle_out = true;
            break;
        case 'H':
            beagle_phased_out = true;
            break;
        case 'K':
            plink_out = true;
            break;
        case 'Z':
            hzar_out = true;
            break;
        case 'Y':
            phylip_out = true;
            break;
        case 'L':
            phylip_var = true;
            break;
        case 'T':
            phylip_var_all = true;
            break;
        case 'U':
            treemix_out = true;
            break;
        case 'g':
            genomic_out = true;
            break;
        case 'W':
            wl_file = optarg;
            break;
        case 'B':
            bl_file = optarg;
            break;
        case 'm':
            min_stack_depth = atoi(optarg);
            break;
        case 'a':
            minor_allele_freq = atof(optarg);
            if (minor_allele_freq > 1)
                minor_allele_freq = minor_allele_freq / 100;

            if (minor_allele_freq < 0 || minor_allele_freq > 0.5) {
                cerr << "Unable to parse the minor allele frequency.\n";
                help();
            }
            break;
        case 'f':
            if (strcasecmp(optarg, "p_value") == 0)
                fst_correction = p_value;
            else if (strcasecmp(optarg, "bonferroni_win") == 0)
                fst_correction = bonferroni_win;
            else if (strcasecmp(optarg, "bonferroni_gen") == 0)
                fst_correction = bonferroni_gen;
            else {
                cerr << "Unknown Fst correction specified '" << optarg << "'\n";
                help();
            }
            break;
        case 'u':
            p_value_cutoff = atof(optarg);
            break;
        case 'e':
            enz = optarg;
            enz.at(0) = tolower(enz.at(0));
            if (renz.count(enz) == 0) {
                cerr << "Unrecognized restriction enzyme specified: '" << enz.c_str() << "'.\n";
                help();
            }
            break;
        case 1005: //sigma
            sigma = atof(optarg);
            break;
        case 'v':
            version();
            exit(1);
            break;
        case '?':
            // getopt_long already printed an error message.
            help();
            break;
        case 1000:
        {
            static const set<string> known_debug_flags = {"VCFCOMP"};
            stringstream ss (optarg);
            string s;
            while (getline(ss, s, ',')) {
                if (known_debug_flags.count(s)) {
                    debug_flags.insert(s);
                } else {
                    cerr << "DEBUG> Error: Unknown error flag '" << s << "'.\n";
                    return -1;
                }
            }
            cerr << "DEBUG> Debug flag(s) : '" << optarg << "'.\n";

            if (debug_flags.count("VCFCOMP") && not write_random_snp) {
                write_single_snp = true;
                cerr << "DEBUG> Added --write_single_snp.\n";
            }

            break;
        }
        default:
            cerr << "Unknown command line option: '" << (char) c << "'\n";
            help();
            abort();
        }
    }

    if (optind < argc) {
        cerr << "Error: Failed to parse command line: '" << argv[optind] << "' is seen as a positional argument. Expected no positional arguments.\n";
        help();
    }

    //
    // Check argument constrains.
    //
    if (input_mode == InputMode::stacks && in_path.empty()) {
        cerr << "Error: Option --v1 requires -P to be given.\n";
        help();
    }

    if (!in_path.empty() && !in_vcf_path.empty()) {
        cerr << "Error: Please specify either '-P/--in_path' or '-V/--in_vcf', not both.\n";
        help();
    } else if (in_path.empty() && in_vcf_path.empty()) {
        cerr << "Error: One of '-P/--in_path' or '-V/--in_vcf' is required.\n";
        help();
    } else if (not in_vcf_path.empty()) {
        input_mode = InputMode::vcf;
    }

    if (input_mode == InputMode::stacks || input_mode == InputMode::stacks2) {

        if (pmap_path.empty())
            cerr << "A population map was not specified, all samples will be read from '" << in_path << "' as a single popultaion.\n";

        if (out_path.empty())
            out_path = in_path;

        out_prefix = "populations";

    } else if (input_mode == InputMode::vcf) {

        if (out_path.empty()) {
            cerr << "Error: Malformed arguments: input mode 'vcf' requires an output directory (--out_path).\n";
            help();
        }

        // Determine out_prefix
        string fname = in_vcf_path;
        if (in_vcf_path.find_last_of('/') != string::npos && in_vcf_path.back() != '/')
            fname = in_vcf_path.substr(in_vcf_path.find_last_of('/')+1);
        size_t trim = 0;
        if (fname.length() > 4 && fname.substr(fname.length()-4) == ".vcf")
            trim = 4;
        else if (fname.length() > 7 && fname.substr(fname.length()-7) == ".vcf.gz")
            trim = 7;
        out_prefix = fname.substr(0, fname.length()-trim);
        out_prefix += ".p";
    }

    // Other
    if (write_single_snp && write_random_snp) {
        cerr << "Please specify either '--write_single_snp' or '--write_random_snp', not both.\n";
        help();
    }

    if (merge_sites == true && enz.length() == 0) {
        cerr << "You must specify the restriction enzyme associated with this data set to merge overlaping cutsites.\n";
        help();
    }

    return 0;
}

void version() {
    cerr << "populations " << VERSION << "\n\n";
}

void help() {
    cerr << "populations " << VERSION << "\n"
         << "Usage:\n"
         << "populations -P dir [-O dir] [-M popmap] (filters) [--fstats] [-k [--sigma=150000] [--bootstrap [-N 100]]] (output formats)\n"
         << "populations -V vcf -O dir [-M popmap] (filters) [--fstats] [-k [--sigma=150000] [--bootstrap [-N 100]]] (output formats)\n"
         << "\n"
         << "  -P,--in_path: path to a directory containing GStacks ouput files.\n"
         << "  -V,--in_vcf: path to a standalone input VCF file.\n"
         << "  -O,--out_path: path to a directory where to write the output files. (Required by -V; otherwise defaults to value of -P.)\n"
         << "  -M,--popmap: path to a population map. (Format is 'SAMPLE1 \\t POP1 \\n SAMPLE2 ...'.)\n"
         << "  -t,--threads: number of threads to run in parallel sections of code.\n"
         << "  --batch_size [int]: the number of loci to process in a batch (default: 10,000 in de novo mode; in reference mode, one chromosome\n"
         << "                      per batch). Increase to speed analysis, uses more memory, decrease to save memory).\n"
         << "\n"
         << "Data Filtering:\n"
         << "  -p [int]: minimum number of populations a locus must be present in to process a locus.\n"
         << "  -r [float]: minimum percentage of individuals in a population required to process a locus for that population.\n"
         << "  --min_maf [float]: specify a minimum minor allele frequency required to process a nucleotide site at a locus (0 < min_maf < 0.5).\n"
         << "  --max_obs_het [float]: specify a maximum observed heterozygosity required to process a nucleotide site at a locus.\n"
         << "  -m [int]: specify a minimum stack depth required for individuals at a locus.\n"
         << "  --lnl_lim [float]: filter loci with log likelihood values below this threshold.\n"
         << "  --write_single_snp: restrict data analysis to only the first SNP per locus.\n"
         << "  --write_random_snp: restrict data analysis to one random SNP per locus.\n"
         << "  -B: path to a file containing Blacklisted markers to be excluded from the export.\n"
         << "  -W: path to a file containing Whitelisted markers to include in the export.\n"
         << "\n"
         << "Merging and Phasing:\n"
         << "  -e,--renz: restriction enzyme name.\n"
         << "  --merge_sites: merge loci that were produced from the same restriction enzyme cutsite (requires reference-aligned data).\n"
         << "  --merge_prune_lim: when merging adjacent loci, if at least X% samples posses both loci prune the remaining samples out of the analysis.\n"
         << "\n"
         << "Fstats:\n"
         << "  --fstats: enable SNP and haplotype-based F statistics.\n"
         << "  --fst_correction: specify a correction to be applied to Fst values: 'p_value', 'bonferroni_win', or 'bonferroni_gen'. Default: off.\n"
         << "  --p_value_cutoff [float]: maximum p-value to keep an Fst measurement. Default: 0.05. (Also used as base for Bonferroni correction.)\n"
         << "\n"
         << "Kernel-smoothing algorithm:\n"
         << "  -k,--smooth: enable kernel-smoothed Pi, Fis, Fst, Fst', and Phi_st calculations.\n"
         << "  --smooth_fstats: enable kernel-smoothed Fst, Fst', and Phi_st calculations.\n"
         << "  --smooth_popstats: enable kernel-smoothed Pi and Fis calculations.\n"
         << "  --sigma [int]: standard deviation of the kernel smoothing weight distribution. Default 150kb.\n"
         << "  --bootstrap: turn on boostrap resampling for all smoothed statistics.\n"
         << "  -N,--bootstrap_reps [int]: number of bootstrap resamplings to calculate (default 100).\n"
         << "  --bootstrap_pifis: turn on boostrap resampling for smoothed SNP-based Pi and Fis calculations.\n"
         << "  --bootstrap_fst: turn on boostrap resampling for smoothed Fst calculations based on pairwise population comparison of SNPs.\n"
         << "  --bootstrap_div: turn on boostrap resampling for smoothed haplotype diveristy and gene diversity calculations based on haplotypes.\n"
         << "  --bootstrap_phist: turn on boostrap resampling for smoothed Phi_st calculations based on haplotypes.\n"
         << "  --bootstrap_wl [path]: only bootstrap loci contained in this whitelist.\n"
         << "\n"
         << "File output options:\n"
         << "  --ordered_export: if data is reference aligned, exports will be ordered; only a single representative of each overlapping site.\n"
         << "  --genomic: output each nucleotide position (fixed or polymorphic) in all population members to a file (requires --renz).\n"
         << "  --fasta_samples: output the sequences of the two haplotypes of each (diploid) sample, for each locus, in FASTA format.\n"
         << "  --fasta_samples_raw: output all haplotypes observed in each sample, for each locus, in FASTA format.\n"
         << "  --fasta_loci: output consensus sequences of all loci, in FASTA format.\n"
         << "  --vcf: output SNPs in Variant Call Format (VCF).\n"
         << "  --vcf_haplotypes: output haplotypes in Variant Call Format (VCF).\n"
         << "  --genepop: output results in GenePop format.\n"
         << "  --structure: output results in Structure format.\n"
         << "  --phase: output genotypes in PHASE format.\n"
         << "  --fastphase: output genotypes in fastPHASE format.\n"
         << "  --beagle: output genotypes in Beagle format.\n"
         << "  --beagle_phased: output haplotypes in Beagle format.\n"
         << "  --plink: output genotypes in PLINK format.\n"
         << "  --hzar: output genotypes in Hybrid Zone Analysis using R (HZAR) format.\n"
         << "  --phylip: output nucleotides that are fixed-within, and variant among populations in Phylip format for phylogenetic tree construction.\n"
         << "  --phylip_var: include variable sites in the phylip output encoded using IUPAC notation.\n"
         << "  --phylip_var_all: include all sequence as well as variable sites in the phylip output encoded using IUPAC notation.\n"
         << "  --treemix: output SNPs in a format useable for the TreeMix program (Pickrell and Pritchard).\n"
         << "\n"
         << "Additional options:\n"
         << "  -h,--help: display this help messsage.\n"
         << "  -v,--version: print program version.\n"
         << "  --verbose: turn on additional logging.\n"
         << ("  --log_fst_comp: log components of Fst/Phi_st calculations to a file.\n");

              // << "    --bootstrap_type [exact|approx]: enable bootstrap resampling for population statistics (reference genome required).\n"

    exit(1);
}
