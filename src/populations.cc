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

#include <dirent.h>
#include <sys/types.h>
#include <sys/stat.h>

#include "export_formats.h"
#include "populations.h"

using namespace std;

// Global variables to hold command-line options.
InputMode input_mode  = InputMode::stacks2;
int       num_threads =  1;
int       batch_id    = -1;
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
bool      sql_out           = false;
bool      vcf_out           = false;
bool      vcf_haplo_out     = false;
bool      fasta_loci_out    = false;
bool      fasta_samples_out = false;
bool      fasta_samples_raw_out = false;
bool      genepop_out       = false;
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
bool      kernel_smoothed   = false;
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

string    out_prefix;

MetaPopInfo mpopi;

set<int> bootstraplist;

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
    output_parameters(cerr);

    //
    // Open and initialize the log file.
    //
    ofstream log_fh;
    open_log(log_fh);
    init_log(log_fh, argc, argv);

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

    BatchLocusProcessor bloc(input_mode, 100, &mpopi);

    //
    // Locate and open input files, read VCF headers, parse population map, load black/white lists.
    //
    bloc.init(batch_id, in_path, pmap_path);

    //
    // Read the next set of loci to process.
    // - If data are denovo, load blim._batch_size loci.
    // - If data are reference aligned, load one chromosome.
    //
    bloc.next_batch(log_fh);
    
    map<int, CSLocus *> catalog = bloc.catalog();

    //
    // Read the bootstrap-whitelist.
    //
    if (bs_wl_file.length() > 0) {
        load_marker_list(bs_wl_file, bootstraplist);
        cerr << "Loaded " << bootstraplist.size() << " markers to include when bootstrapping.\n";
    }
    set<int> blacklist(bloc.blacklist());
    map<int, set<int>> whitelist(bloc.whitelist());
    
    //
    // Retrieve the genomic order of loci.
    //
    loci_ordered = order_unordered_loci(catalog);

    // Report information on the MetaPopInfo.
    mpopi.status();

    if (size_t(population_limit) > mpopi.pops().size()) {
        cerr << "Notice: Population limit (" << population_limit << ")"
             << " larger than number of popualtions present, adjusting parameter to "
             << mpopi.pops().size() << "\n";
        population_limit = mpopi.pops().size();
    }

    //
    // Initialize the PopMap
    //
    cerr << "Populating observed haplotypes for " << mpopi.samples().size() << " samples, " << catalog.size() << " loci.\n";
    PopMap<CSLocus> *pmap = new PopMap<CSLocus>(mpopi, catalog.size());

    // Using Stacks v2 files.
    if (input_mode == InputMode::stacks2)
        pmap->populate(catalog, bloc.cloc_vcf_records(), *bloc.vcf_header());
    // ...or using VCF records.
    else if (input_mode == InputMode::vcf)
        pmap->populate(catalog, bloc.ext_vcf_records(), *bloc.vcf_header());

    //
    // Tabulate haplotypes present and in what combinations.
    //
    tabulate_haplotypes(catalog, pmap);

    //
    // Output a list of heterozygous loci and the associate haplotype frequencies.
    //
    if (sql_out)
        write_sql(catalog, pmap);

    log_fh << "# Distribution of population loci.\n";
    log_haplotype_cnts(catalog, log_fh);

    apply_locus_constraints(catalog, pmap, log_fh);
    if (pmap->loci_cnt() == 0) {
        cerr << "Error: All loci have been filtered out.\n";
        throw exception();
    }

    log_fh << "# Distribution of population loci after applying locus constraints.\n";
    log_haplotype_cnts(catalog, log_fh);

    //
    // Create the PopSum object and compute the summary statistics.
    //

    PopSum<CSLocus> *psum = new PopSum<CSLocus>(*pmap, mpopi);
    for (size_t i=0; i<mpopi.pops().size(); ++i) {
        cerr << "Generating nucleotide-level summary statistics for population '" << mpopi.pops()[i].name << "'\n";
        psum->add_population(catalog, pmap, i, verbose, log_fh);
    }

    cerr << "Tallying loci across populations...";
    psum->tally(catalog);
    cerr << "done.\n";

    //
    // We have removed loci that were below the -r and -p thresholds. Now we need to
    // identify individual SNPs that are below the -r threshold or the minor allele
    // frequency threshold (-a). In these cases we will remove the SNP, but keep the locus.
    //
    blacklist.clear();
    int pruned_snps = prune_polymorphic_sites(catalog, pmap, psum, whitelist, blacklist, log_fh);
    cerr << "Pruned " << pruned_snps << " variant sites due to filter constraints (more with --verbose).\n";

    //
    // Create an artificial whitelist if the user requested only the first or a random SNP per locus.
    //
    if (write_single_snp)
        implement_single_snp_whitelist(catalog, psum, whitelist);
    else if (write_random_snp)
        implement_random_snp_whitelist(catalog, psum, whitelist);

    //
    // Remove the accumulated SNPs
    //
    cerr << "Removing " << blacklist.size() << " additional loci for which all variant sites were filtered...";
    set<int> empty_list;
    reduce_catalog(catalog, empty_list, blacklist);
    reduce_catalog_snps(catalog, whitelist, pmap);
    int retained = pmap->prune(blacklist);
    cerr << " retained " << retained << " loci.\n";
    if (pmap->loci_cnt() == 0) {
        cerr << "Error: All loci have been filtered out.\n";
        throw exception();
    }

    //
    // Merge loci that overlap on a common restriction enzyme cut site.
    //
    map<int, pair<merget, int> > merge_map;
    if (merge_sites && loci_ordered)
        merge_shared_cutsite_loci(catalog, pmap, psum, merge_map, log_fh);

    //
    // Regenerate summary statistics after pruning SNPs and  merging loci.
    //

    if (debug_flags.count("VCFCOMP"))
        vcfcomp_simplify_pmap(catalog, pmap);

    delete psum;
    psum = new PopSum<CSLocus>(*pmap, mpopi);
    for (size_t i=0; i<mpopi.pops().size(); ++i) {
        cerr << "Regenerating nucleotide-level summary statistics for population '" << mpopi.pops()[i].name << "'\n";
        psum->add_population(catalog, pmap, i, verbose, log_fh);
    }
    cerr << "Re-tallying loci across populations...";
    psum->tally(catalog);
    cerr << "done.\n";

    if (kernel_smoothed) {
        if (loci_ordered) {
            for (size_t i=0; i<mpopi.pops().size(); ++i) {
                cerr << "  Generating kernel-smoothed population statistics for population '" << mpopi.pops()[i].name << "'...\n";
                kernel_smoothed_popstats(catalog, pmap, psum, i, log_fh);
            }
        } else {
            cerr << "Notice: Smoothing was requested (-k), but will not be performed as the loci are not ordered.\n";
        }
    }

    //
    // Log the SNPs per locus distribution.
    //
    log_snps_per_loc_distrib(log_fh, catalog);

    calculate_haplotype_stats(catalog, pmap, psum);

    if (calc_fstats) {
        calculate_haplotype_divergence(catalog, pmap, psum);
        calculate_haplotype_divergence_pairwise(catalog, pmap, psum);
    }

    //
    // Calculate and output the locus-level summary statistics.
    //
    calculate_summary_stats(catalog, pmap, psum);

    //
    // Output the observed haplotypes.
    //
    write_generic(catalog, pmap, false);

    //
    // Output data in requested formats
    //
    if (fasta_loci_out)
        write_fasta_loci(catalog, pmap);

    if (fasta_samples_out)
        write_fasta_samples(catalog, pmap);

    if (fasta_samples_raw_out)
        write_fasta_samples_raw(catalog, pmap);

    if (genepop_out && ordered_export)
        write_genepop_ordered(catalog, pmap, psum, log_fh);
    else if (genepop_out)
        write_genepop(catalog, pmap, psum);

    if (structure_out && ordered_export)
        write_structure_ordered(catalog, pmap, psum, log_fh);
    else if (structure_out)
        write_structure(catalog, pmap, psum);

    if (fastphase_out)
        write_fastphase(catalog, pmap, psum);

    if (phase_out)
        write_phase(catalog, pmap, psum);

    if (beagle_out)
        write_beagle(catalog, pmap, psum);

    if (beagle_phased_out)
        write_beagle_phased(catalog, pmap, psum);

    if (plink_out)
        write_plink(catalog, pmap, psum);

    if (hzar_out)
        write_hzar(catalog, pmap, psum);

    if (treemix_out)
        write_treemix(catalog, pmap, psum);

    if (phylip_out || phylip_var)
        write_phylip(catalog, pmap, psum);

    if (phylip_var_all)
        write_fullseq_phylip(catalog, pmap, psum);

    if (vcf_haplo_out)
        write_vcf_haplotypes(catalog, pmap, psum);

    if (vcf_out && ordered_export)
        write_vcf_ordered(catalog, pmap, psum, merge_map, log_fh);
    else if (vcf_out)
        write_vcf(catalog, pmap, psum, merge_map);

    //
    // Calculate and write Fst.
    //
    if (calc_fstats)
        write_fst_stats(catalog, pmap, psum, log_fh);

    //
    // Output nucleotide-level genotype calls for each individual.
    //
    if (genomic_out)
        write_genomic(catalog, pmap);

    for (auto& cloc : catalog)
        delete cloc.second;
    delete psum;
    delete pmap;
    log_fh.close();

    cerr << "Populations is done.\n";
    return 0;
    IF_NDEBUG_CATCH_ALL_EXCEPTIONS
}

int
BatchLocusProcessor::init(int batch_id, string in_path, string pmap_path)
{
    //
    // Read the blacklist and whitelist to control which loci we load..
    //
    if (bl_file.length() > 0) {
        load_marker_list(bl_file, this->_blacklist);
        cerr << "Loaded " << this->_blacklist.size() << " blacklisted markers.\n";
    }
    if (wl_file.length() > 0) {
        load_marker_column_list(wl_file, this->_whitelist);
        cerr << "Loaded " << this->_whitelist.size() << " whitelisted markers.\n";
        //// check_whitelist_integrity(catalog, whitelist);
    }
    
    if (this->_input_mode == InputMode::vcf)
        this->init_external_loci(in_path, pmap_path);
    else
        this->init_stacks_loci(batch_id, in_path, pmap_path);

    return 0;
}

size_t
BatchLocusProcessor::next_batch(ostream &log_fh)
{
    size_t loc_cnt;
    
    if (this->_input_mode == InputMode::vcf)
        loc_cnt = this->next_batch_external_loci(log_fh);
    else
        loc_cnt = this->next_batch_stacks_loci();

    return loc_cnt;
}

int
BatchLocusProcessor::init_stacks_loci(int batch_id, string in_path, string pmap_path)
{
    //
    // Open the files.
    //
    string catalog_fa_path  = in_path + "batch_" + to_string(batch_id) + ".gstacks.fa.gz";
    string catalog_vcf_path = in_path + "batch_" + to_string(batch_id) + ".gstacks.vcf.gz";

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
    // Store the VCF header data.
    //
    this->_vcf_header = new VcfHeader(this->cloc_reader().header());

    return 0;
}

size_t
BatchLocusProcessor::next_batch_stacks_loci()
{
    if (this->_catalog != NULL)
        delete this->_catalog;
    this->_catalog = new map<int, CSLocus *>();
    if (this->_cloc_vcf_rec != NULL)
        delete _cloc_vcf_rec;
    this->_cloc_vcf_rec = new unordered_map<int, vector<VcfRecord>>();

    // Read the files, create the loci.
    vector<VcfRecord> records;
    Seq seq;
    seq.id      = new char[id_len];
    seq.comment = new char[id_len];
    seq.seq     = new char[max_len];

    int    cloc_id, rv;
    size_t loc_cnt = 0;

    while (loc_cnt < this->_batch_size && this->_cloc_reader.read_one_locus(records)) {
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
        // Check if this locus is filtered.
        //
        if (this->_whitelist.size() > 0 && this->_whitelist.count(cloc_id) == 0) {
            records.clear();
            continue;
        }
        if (this->_blacklist.count(cloc_id)) {
            records.clear();
            continue;
        }

        // Create the CSLocus.
        this->_catalog->insert({cloc_id, new_cslocus(seq, records, cloc_id)});

        // Save the records (they are needed for the PopMap object).
        (*this->_cloc_vcf_rec)[cloc_id] = move(records);

        loc_cnt++;
    }

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
    // Store the VCF header data.
    //
    this->_vcf_header = new VcfHeader(this->vcf_reader().header());

    return 0;
}

size_t
BatchLocusProcessor::next_batch_external_loci(ostream &log_fh)
{
    //
    // VCF mode
    //
    
    if (this->_ext_vcf_rec != NULL)
        delete this->_ext_vcf_rec;
    this->_ext_vcf_rec = new vector<VcfRecord>();

    vector<size_t> skipped_notsnp;
    vector<size_t> skipped_filter;
    size_t loc_cnt = 0;

    this->_ext_vcf_rec->push_back(VcfRecord());
    VcfRecord* rec = & this->_ext_vcf_rec->back();

    while (loc_cnt < this->_batch_size && this->_vcf_parser.next_record(*rec)) {
        // Check for a SNP.
        if (not rec->is_snp()) {
            skipped_notsnp.push_back(this->_vcf_parser.line_number());
            continue;
        }

        // Check for a filtered-out SNP
        if (strncmp(rec->filters(), ".", 2) != 0 && strncmp(rec->filters(), "PASS", 5) != 0) {
            skipped_filter.push_back(this->_vcf_parser.line_number());
            continue;
        }

        // Save the SNP.
        this->_ext_vcf_rec->push_back(VcfRecord());
        rec = & this->_ext_vcf_rec->back();

        loc_cnt++;
    }

    this->_ext_vcf_rec->pop_back();

    cerr << "Found " << this->_ext_vcf_rec->size() << " SNP records in file '" << in_vcf_path
         << "'. (Skipped " << skipped_filter.size() << " already filtered-out SNPs and "
         << skipped_notsnp.size() << " non-SNP records ; more with --verbose.)\n";
    if (verbose && not skipped_notsnp.empty()) {
        log_fh << "The following VCF record lines were determined not to be SNPs and skipped :";
        for (vector<size_t>::const_iterator l=skipped_notsnp.begin(); l!=skipped_notsnp.end(); ++l)
            log_fh << " " << *l;
        log_fh << "\n";
    }
    if (this->_ext_vcf_rec->size() == 0) {
        cerr << "Error: No records.\n";
        throw exception();
    }

    if (this->_catalog != NULL)
        delete this->_catalog;
    this->_catalog = create_catalog(this->ext_vcf_records());

    return loc_cnt;
}

void
vcfcomp_simplify_pmap (map<int, CSLocus*>& catalog, PopMap<CSLocus>* pmap)
{
    cerr << "DEBUG Deleting information from the pmap & catalog so that they resemble what can be retrieved from a VCF.\n";
    // n.b. In this configuration we only have one SNP per locus so we don't have
    // to worry about what U's imply regarding haplotypes.
    size_t n_deleted = 0;
    size_t n_loci = pmap->loci_cnt();
    size_t n_samples = pmap->sample_cnt();
    for (size_t l=0; l<n_loci; ++l) {
        CSLocus* loc = catalog.at(pmap->rev_locus_index(l));
        if (loc->snps.empty())
            continue;
        if (loc->snps.size() > 1) {
            cerr << "Error: This requires --write_single_snp.\n";
            throw exception();
        }
        size_t col = loc->snps.at(0)->col;
        Datum** datums = pmap->locus(loc->id);
        for (size_t s=0; s<n_samples; ++s) {
            if (datums[s] == NULL)
                continue;
            if (size_t(datums[s]->len) <= col
                    || datums[s]->model[col] == 'U'
                    || datums[s]->obshap.empty()
                    || strcmp(datums[s]->obshap[0], "N") == 0
                    ) {
                delete datums[s];
                datums[s] = NULL;
                --loc->cnt;
                --loc->hcnt;
                ++n_deleted;
            }
        }
    }
    cerr << "? Deleted " << n_deleted << " 'Datums'.\n";

    set<int> myblacklist;
    // All NULL.
    for (auto& l : catalog) {
        Datum** data = pmap->locus(l.second->id);
        size_t non_null = 0;
        for (size_t i=0; i < size_t(pmap->sample_cnt()); ++i)
            if (data[i]!=NULL)
                ++non_null;
        if (non_null == 0)
            myblacklist.insert(l.second->id);
    }
    // Not two alleles.
    for (auto& l : catalog) {
        if (l.second->snps.empty()
                || l.second->snps[0]->rank_2 == 0
                || l.second->snps[0]->rank_3 != 0
                || l.second->alleles.size() != 2) {
            myblacklist.insert(l.second->id);
        } else {
            // Check the actual number of alleles
            Datum** data = pmap->locus(l.second->id);
            set<char> seen_alleles;
            for (size_t i=0; i<size_t(pmap->sample_cnt()); ++i) {
                Datum* d = data[i];
                if (d != NULL)
                    for(char* hapl : d->obshap)
                        seen_alleles.insert(hapl[0]);
            }
            if ((int((seen_alleles.count('A') || seen_alleles.count('a')))
                    + (seen_alleles.count('C') || seen_alleles.count('c'))
                    + (seen_alleles.count('T') || seen_alleles.count('t'))
                    + (seen_alleles.count('G') || seen_alleles.count('g')))
                    != 2)
                myblacklist.insert(l.second->id);
        }
    }
    // Same SNP in different loci.
    for (auto& chr : pmap->ordered_loci()) {
        map<size_t,vector<size_t> > seen_bp0; // (bp, [loc_id's])
        for (CSLocus* loc : chr.second)
            if (not loc->snps.empty())
                seen_bp0[loc->sort_bp(loc->snps[0]->col)].push_back(loc->id);
        for (auto& bp : seen_bp0)
            if (bp.second.size() > 1)
                for (size_t loc_id : bp.second)
                    myblacklist.insert(loc_id);
    }
    set<int> empty;
    reduce_catalog(catalog, empty, myblacklist);
    pmap->prune(myblacklist);
    cerr << "? Now working on " << catalog.size() << " loci (deleted " << myblacklist.size() << " loci).\n";
}

int
apply_locus_constraints(map<int, CSLocus *> &catalog,
                        PopMap<CSLocus> *pmap,
                        ofstream &log_fh)
{
    uint pop_sthg;
    CSLocus *loc;
    Datum  **d;

    if (sample_limit == 0 && population_limit == 0 && min_stack_depth == 0) return 0;

    if (verbose)
        log_fh << "\n#\n# List of loci removed by first filtering stage of sample and population constraints\n#\n"
               << "# Action\tLocus ID\tChr\tBP\tColumn\tReason\n";

    map<int, CSLocus *>::iterator it;

    uint pop_cnt   = mpopi.pops().size();
    int *pop_order = new int [pop_cnt];

    // Which population each sample belongs to.
    int *samples   = new int [pmap->sample_cnt()];

    // For the current locus, how many samples in each population.
    int *pop_cnts  = new int [pop_cnt];

    // The total number of samples in each population.
    int *pop_tot   = new int [pop_cnt];

    pop_sthg = 0;
    for (size_t i_pop=0; i_pop<mpopi.pops().size(); ++i_pop) {
        const Pop& pop = mpopi.pops()[i_pop];
        pop_tot[pop_sthg]  = 0;

        for (uint i = pop.first_sample; i <= pop.last_sample; i++) {
            samples[i] = pop_sthg;
            pop_tot[pop_sthg]++;
        }
        pop_order[pop_sthg] = i_pop;
        pop_sthg++;
    }

    for (uint i = 0; i < pop_cnt; i++)
        pop_cnts[i] = 0;

    double pct       = 0.0;
    bool   pop_limit = false;
    int    pops      = 0;
    int    below_stack_dep  = 0;
    uint   below_lnl_thresh = 0;
    set<int> blacklist;

    for (it = catalog.begin(); it != catalog.end(); it++) {
        loc = it->second;
        d   = pmap->locus(loc->id);

        for (int i = 0; i < pmap->sample_cnt(); i++) {
            //
            // Check that each sample is over the minimum stack depth for this locus.
            //
            if (d[i] != NULL &&
                min_stack_depth > 0 &&
                d[i]->tot_depth < min_stack_depth) {
                below_stack_dep++;
                delete d[i];
                d[i] = NULL;
                loc->hcnt--;
            }

            //
            // Check that each sample is over the log likelihood threshold.
            //
            if (d[i] != NULL &&
                filter_lnl   &&
                d[i]->lnl < lnl_limit) {
                below_lnl_thresh++;
                delete d[i];
                d[i] = NULL;
                loc->hcnt--;
            }
        }

        //
        // Tally up the count of samples in this population.
        //
        for (int i = 0; i < pmap->sample_cnt(); i++) {
            if (d[i] != NULL)
                pop_cnts[samples[i]]++;
        }

        //
        // Check that the counts for each population are over sample_limit. If not, zero out
        // the members of that population.
        //
        for (uint i = 0; i < pop_cnt; i++) {
            const Pop& pop = mpopi.pops()[pop_order[i]];

            pct = (double) pop_cnts[i] / (double) pop_tot[i];

            if (pop_cnts[i] > 0 && pct < sample_limit) {
                //cerr << "Removing population " << pop_order[i] << " at locus: " << loc->id << "; below sample limit: " << pct << "\n";
                for (uint j  = pop.first_sample; j <= pop.last_sample; j++) {
                    if (d[j] != NULL) {
                        delete d[j];
                        d[j] = NULL;
                        loc->hcnt--;
                    }
                }
                pop_cnts[i] = 0;
            }
        }

        //
        // Check that this locus is present in enough populations.
        //
        for (uint i = 0; i < pop_cnt; i++)
            if (pop_cnts[i] > 0) pops++;
        if (pops < population_limit) {
            //cerr << "Removing locus: " << loc->id << "; below population limit: " << pops << "\n";
            pop_limit = true;
        }

        if (pop_limit) {
            blacklist.insert(loc->id);

            if (verbose)
                log_fh << "removed_locus\t"
                       << loc->id << "\t"
                       << loc->loc.chr() << "\t"
                       << loc->sort_bp() +1 << "\t"
                       << 0 << "\tfailed_population_limit\n";
        }

        for (uint i = 0; i < pop_cnt; i++)
            pop_cnts[i] = 0;
        pop_limit = false;
        pops      = 0;
    }

    //
    // Remove loci
    //
    if (min_stack_depth > 0)
        cerr << "Removed " << below_stack_dep << " samples from loci that are below the minimum stack depth of " << min_stack_depth << "x\n";
    if (filter_lnl)
        cerr << "Removed " << below_lnl_thresh << " samples from loci that are below the log likelihood threshold of " << lnl_limit << "\n";
    cerr << "Removing " << blacklist.size() << " loci that did not pass sample/population constraints...";
    set<int> whitelist;
    reduce_catalog(catalog, whitelist, blacklist);
    int retained = pmap->prune(blacklist);
    cerr << " retained " << retained << " loci.\n";

    delete [] pop_cnts;
    delete [] pop_tot;
    delete [] pop_order;
    delete [] samples;

    return 0;
}

int
prune_polymorphic_sites(map<int, CSLocus *> &catalog,
                        PopMap<CSLocus> *pmap,
                        PopSum<CSLocus> *psum,
                        map<int, set<int> > &whitelist, set<int> &blacklist,
                        ofstream &log_fh)
{
    map<int, set<int> > new_wl;
    vector<int> pop_prune_list;
    CSLocus  *loc;
    LocTally *t;
    LocSum  **s;
    Datum   **d;
    bool      sample_prune, maf_prune, het_prune, inc_prune;
    int       size, pruned = 0;

    if (verbose)
        log_fh << "\n#\n# List of pruned nucleotide sites\n#\n"
               << "# Action\tLocus ID\tChr\tBP\tColumn\tReason\n";

    //
    // If the whitelist is populated, use it as a guide for what loci to consider.
    //
    // Construct a new whitelist along the way, that is a subset of the existing list.
    //
    if (whitelist.size() > 0) {
        map<int, set<int> >::iterator it;

        for (it = whitelist.begin(); it != whitelist.end(); it++) {
            //
            // A locus on the whitelist may have already been filtered out.
            //
            if (catalog.count(it->first) == 0)
                continue;

            loc = catalog[it->first];
            t   = psum->locus_tally(loc->id);
            s   = psum->locus(loc->id);

            //
            // Check that each SNP in this locus is above the sample_limit and that
            // each SNP is above the minor allele frequency. If so, add it back to
            // the whiteliest.
            //
            size = it->second.size();
            for (uint i = 0; i < loc->snps.size(); i++) {

                //
                // If it is not already in the whitelist, ignore it.
                //
                if (size > 0 && it->second.count(loc->snps[i]->col) == 0)
                    continue;

                //
                // If the site is fixed, ignore it.
                //
                if (t->nucs[loc->snps[i]->col].fixed == true)
                    continue;

                sample_prune = false;
                maf_prune    = false;
                het_prune    = false;
                inc_prune    = false;
                pop_prune_list.clear();

                for (size_t p=0; p<mpopi.pops().size(); ++p) {
                    if (s[p]->nucs[loc->snps[i]->col].incompatible_site)
                        inc_prune = true;

                    else if (s[p]->nucs[loc->snps[i]->col].num_indv == 0 ||
                             (double) s[p]->nucs[loc->snps[i]->col].num_indv / (double) psum->pop_size(p) < sample_limit)
                        pop_prune_list.push_back(p);
                }

                //
                // Check how many populations have to be pruned out due to sample limit. If less than
                // population limit, prune them; if more than population limit, mark locus for deletion.
                //
                if ((mpopi.pops().size() - pop_prune_list.size()) < (uint) population_limit) {
                    sample_prune = true;
                } else {
                    for (size_t p : pop_prune_list) {
                        if (s[p]->nucs[loc->snps[i]->col].num_indv == 0)
                            continue;

                        d = pmap->locus(loc->id);
                        const Pop& pop = mpopi.pops()[p];
                        for (uint k = pop.first_sample; k <= pop.last_sample; k++) {
                            if (d[k] == NULL || loc->snps[i]->col >= (uint) d[k]->len)
                                continue;
                            if (d[k]->model != NULL) {
                                d[k]->model[loc->snps[i]->col] = 'U';
                            }
                        }
                    }
                }

                if (t->nucs[loc->snps[i]->col].allele_cnt > 1) {
                    //
                    // Test for minor allele frequency.
                    //
                    if ((1 - t->nucs[loc->snps[i]->col].p_freq) < minor_allele_freq)
                        maf_prune = true;
                    //
                    // Test for observed heterozygosity.
                    //
                    if (t->nucs[loc->snps[i]->col].obs_het > max_obs_het)
                        het_prune = true;
                }

                if (maf_prune == false && het_prune == false && sample_prune == false && inc_prune == false) {
                    new_wl[loc->id].insert(loc->snps[i]->col);
                } else {
                    pruned++;
                    if (verbose) {
                        log_fh << "pruned_polymorphic_site\t"
                               << loc->id << "\t"
                               << loc->loc.chr() << "\t"
                               << loc->sort_bp(loc->snps[i]->col) +1 << "\t"
                               << loc->snps[i]->col << "\t";
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
            if (new_wl.count(loc->id) == 0) {
                if (verbose)
                    log_fh << "removed_locus\t"
                           << loc->id << "\t"
                           << loc->loc.chr() << "\t"
                           << loc->sort_bp() +1 << "\t"
                           << 0 << "\tno_snps_remaining\n";
                blacklist.insert(loc->id);
            }
        }

    } else {
        //
        // Otherwise, just iterate over the catalog.
        //
        map<int, CSLocus *>::iterator it;
        for (it = catalog.begin(); it != catalog.end(); it++) {
            loc = it->second;

            //
            // If this locus is fixed, don't try to filter it out.
            //
            if (loc->snps.size() == 0) {
                new_wl.insert(make_pair(loc->id, set<int>()));
                continue;
            }

            t = psum->locus_tally(loc->id);
            s = psum->locus(loc->id);

            for (uint i = 0; i < loc->snps.size(); i++) {

                //
                // If the site is fixed, ignore it.
                //
                if (t->nucs[loc->snps[i]->col].fixed == true)
                    continue;

                sample_prune = false;
                maf_prune    = false;
                het_prune    = false;
                inc_prune    = false;
                pop_prune_list.clear();

                for (size_t p = 0; p < mpopi.pops().size(); p++) {
                    if (s[p]->nucs[loc->snps[i]->col].incompatible_site)
                        inc_prune = true;
                    else if (s[p]->nucs[loc->snps[i]->col].num_indv == 0 ||
                             (double) s[p]->nucs[loc->snps[i]->col].num_indv / (double) psum->pop_size(p) < sample_limit)
                        pop_prune_list.push_back(p);
                }

                //
                // Check how many populations have to be pruned out due to sample limit. If less than
                // population limit, prune them; if more than population limit, mark locus for deletion.
                //
                if ((psum->pop_cnt() - pop_prune_list.size()) < (uint) population_limit) {
                    sample_prune = true;
                } else {
                    for (size_t p : pop_prune_list) {
                        if (s[p]->nucs[loc->snps[i]->col].num_indv == 0)
                            continue;

                        d = pmap->locus(loc->id);
                        const Pop& pop = mpopi.pops()[p];
                        for (uint k = pop.first_sample; k <= pop.last_sample; k++) {
                            if (d[k] == NULL || loc->snps[i]->col >= (uint) d[k]->len)
                                continue;
                            if (d[k]->model != NULL) {
                                d[k]->model[loc->snps[i]->col] = 'U';
                            }
                        }
                    }
                }

                if (t->nucs[loc->snps[i]->col].allele_cnt > 1) {
                    //
                    // Test for minor allele frequency.
                    //
                    if ((1 - t->nucs[loc->snps[i]->col].p_freq) < minor_allele_freq)
                        maf_prune = true;
                    //
                    // Test for observed heterozygosity.
                    //
                    if (t->nucs[loc->snps[i]->col].obs_het > max_obs_het)
                        het_prune = true;
                }

                if (maf_prune == false && het_prune == false && sample_prune == false && inc_prune == false) {
                    new_wl[loc->id].insert(loc->snps[i]->col);
                } else {
                    pruned++;
                    if (verbose) {
                        log_fh << "pruned_polymorphic_site\t"
                               << loc->id << "\t"
                               << loc->loc.chr() << "\t"
                               << loc->sort_bp(loc->snps[i]->col) +1 << "\t"
                               << loc->snps[i]->col << "\t";
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
            if (new_wl.count(loc->id) == 0) {
                if (verbose)
                    log_fh << "removed_locus\t"
                           << loc->id << "\t"
                           << loc->loc.chr() << "\t"
                           << loc->sort_bp() +1 << "\t"
                           << 0 << "\tno_snps_remaining\n";
                blacklist.insert(loc->id);
            }
        }
    }

    whitelist = new_wl;

    return pruned;
}

bool
order_unordered_loci(map<int, CSLocus *> &catalog)
{
    map<int, CSLocus *>::iterator it;
    CSLocus *loc;
    set<string> chrs;

    for (it = catalog.begin(); it != catalog.end(); it++) {
        loc = it->second;
        if (!loc->loc.empty())
            chrs.insert(loc->loc.chr());
    }

    //
    // This data is already reference aligned.
    //
    if (chrs.size() > 0)
        return true;

    cerr << "Catalog is not reference aligned, arbitrarily ordering catalog loci.\n";

    uint bp = 1;
    for (it = catalog.begin(); it != catalog.end(); it++) {
        loc = it->second;
        loc->loc = PhyLoc("un", bp);
        bp += strlen(loc->con);
    }

    return false;
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
            create_genotype_map(loc, pmap);
            call_population_genotypes(loc, pmap);
        }

        loc->lnl = mean / cnt;
    }

    return 0;
}

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
create_genotype_map(CSLocus *locus, PopMap<CSLocus> *pmap)
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

    Datum **d;
    map<string, int> haplotypes;
    map<string, int>::iterator k;
    vector<pair<string, int> > sorted_haplotypes;

    d = pmap->locus(locus->id);

    for (int i = 0; i < pmap->sample_cnt(); i++) {

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

int call_population_genotypes(CSLocus *locus,
                              PopMap<CSLocus> *pmap) {
    //
    // Fetch the array of observed haplotypes from the population
    //
    Datum **d = pmap->locus(locus->id);

    for (int i = 0; i < pmap->sample_cnt(); i++) {
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

int write_genomic(map<int, CSLocus *> &catalog, PopMap<CSLocus> *pmap) {
    string file = out_path + out_prefix + ".genomic.tsv";

    ofstream fh(file.c_str(), ofstream::out);

    if (fh.fail()) {
        cerr << "Error opening genomic output file '" << file << "'\n";
        exit(1);
    }

    uint rcnt = enz.length() ? renz_cnt[enz] : 0;
    uint rlen = enz.length() ? renz_len[enz] : 0;

    //
    // Count the number of markers that have enough samples to output.
    //
    map<int, CSLocus *>::iterator cit;
    const CSLocus *loc;
    int num_loci = 0;

    for (cit = catalog.begin(); cit != catalog.end(); cit++) {
        loc = cit->second;

        uint start = 0;
        uint end   = loc->len;
        if (end > rlen) {
            for (uint n = 0; n < rcnt; n++)
                if (strncmp(loc->con, renz[enz][n], rlen) == 0)
                    start += renz_len[enz];
        }
        num_loci += end - start;
    }
    cerr << "Writing " << num_loci << " nucleotide positions to genomic file, '" << file << "'\n";

    //
    // Write the header
    //
    fh << num_loci << "\t" << pmap->sample_cnt() << "\n";

    //
    // Output each locus.
    //
    map<string, vector<CSLocus *> >::const_iterator it;
    int  a, b;

    for (it = pmap->ordered_loci().begin(); it != pmap->ordered_loci().end(); it++) {
        for (uint i = 0; i < it->second.size(); i++) {
            loc = it->second[i];

            Datum **d = pmap->locus(loc->id);
            set<int> snp_locs;
            string   obshap;

            for (uint i = 0; i < loc->snps.size(); i++)
                snp_locs.insert(loc->snps[i]->col);

            uint start = 0;
            uint end   = loc->len;
            //
            // Check for the existence of the restriction enzyme cut site, mask off
            // its output.
            //
            if (end > rlen) {
                for (uint n = 0; n < rcnt; n++)
                    if (strncmp(loc->con, renz[enz][n], rlen) == 0)
                        start += renz_len[enz];
            }

            uint k = 0;
            for (uint n = start; n < end; n++) {
                fh << loc->id << "\t" << loc->loc.chr() << "\t" << loc->sort_bp(n) +1;

                if (snp_locs.count(n) == 0) {
                    for (int j = 0; j < pmap->sample_cnt(); j++) {
                        a = encode_gtype(loc->con[n]);
                        fh << "\t" << encoded_gtypes[a][a];
                    }
                } else {
                    for (int j = 0; j < pmap->sample_cnt(); j++) {
                        fh << "\t";

                        if (d[j] == NULL)
                            fh << "0";
                        else
                            switch (d[j]->obshap.size()) {
                            case 1:
                                a = encode_gtype(d[j]->obshap[0][k]);
                                fh << encoded_gtypes[a][a];
                                break;
                            case 2:
                                a = encode_gtype(d[j]->obshap[0][k]);
                                b = encode_gtype(d[j]->obshap[1][k]);
                                fh << encoded_gtypes[a][b];
                                break;
                            default:
                                fh << "0";
                                break;
                            }
                    }
                    k++;
                }
                fh << "\n";
            }
        }
    }

    fh.close();

    return 0;
}

int
calculate_haplotype_stats(map<int, CSLocus *> &catalog, PopMap<CSLocus> *pmap, PopSum<CSLocus> *psum)
{
    map<string, vector<CSLocus *> >::const_iterator it;
    const CSLocus  *loc;
    Datum   **d;
    LocStat  *l;

    //
    // Instantiate the kernel smoothing and bootstrap objects if requested.
    //
    KSmooth<LocStat>     *ks=NULL;
    OHaplotypes<LocStat> *ord=NULL;
    Bootstrap<LocStat>   *bs=NULL;
    if (kernel_smoothed && loci_ordered) {
        ks  = new KSmooth<LocStat>(2);
        ord = new OHaplotypes<LocStat>();
    }

    //
    // Open output file and print header.
    //
    string file = out_path + out_prefix + ".hapstats.tsv";

    ofstream fh(file.c_str(), ofstream::out);
    if (fh.fail()) {
        cerr << "Error opening haplotype stats file '" << file << "'\n";
        exit(1);
    }
    fh.precision(fieldw);
    fh.setf(std::ios::fixed);

    //
    // Write the population members.
    //
    for (auto& pop : mpopi.pops()) {
        fh << "# " << pop.name << "\t";
        for (size_t i = pop.first_sample; i <= pop.last_sample; i++) {
            fh << mpopi.samples()[i].name;
            if (i < pop.last_sample)
                fh << ",";
        }
        fh << "\n";
    }

    fh << "# Batch ID "    << "\t"
       << "Locus ID"       << "\t"
       << "Chr"            << "\t"
       << "BP"             << "\t"
       << "Pop ID"         << "\t"
       << "N"              << "\t"
       << "Haplotype Cnt"  << "\t"
       << "Gene Diversity" << "\t"
       << "Smoothed Gene Diversity"      << "\t"
       << "Smoothed Gene Diversity P-value"      << "\t"
       << "Haplotype Diversity"          << "\t"
       << "Smoothed Haplotype Diversity" << "\t"
       << "Smoothed Haplotype Diversity P-value" << "\t"
       << "Haplotypes"                   << "\n";

    //
    // Iterate over the members of each population.
    //
    for (auto& pop : mpopi.pops()) {

        cerr << "Generating haplotype-level summary statistics for population '" << pop.name << "'\n";
        map<string, vector<LocStat *> > genome_locstats;

        for (it = pmap->ordered_loci().begin(); it != pmap->ordered_loci().end(); it++) {

            if (bootstrap_div)
                bs = new Bootstrap<LocStat>(2);

            vector<LocStat *> &locstats = genome_locstats[it->first];
            map<uint, uint>    locstats_key;
            ord->order(locstats, locstats_key, it->second);

            for (uint pos = 0; pos < it->second.size(); pos++) {
                loc = it->second[pos];
                d   = pmap->locus(loc->id);

                if (loc->snps.size() == 0)
                    continue;

                // cerr << "Looking at locus " << loc->id << "\n";

                l = haplotype_diversity(pop.first_sample, pop.last_sample, d);

                if (l != NULL) {
                    l->loc_id = loc->id;
                    l->bp     = loc->sort_bp();
                    locstats[locstats_key[l->bp]] = l;
                }
            }

            if (kernel_smoothed && loci_ordered) {
                cerr << "    Generating kernel-smoothed statistics on chromosome " << it->first << "\n";
                ks->smooth(locstats);
            }

            if (bootstrap_div)
                bs->add_data(locstats);
        }

        for (it = pmap->ordered_loci().begin(); it != pmap->ordered_loci().end(); it++) {
            vector<LocStat *> &locstats = genome_locstats[it->first];

            if (bootstrap_div)
                bs->execute(locstats);

            //
            // Write results.
            //
            for (uint k = 0; k < locstats.size(); k++) {
                l = locstats[k];
                if (l == NULL) continue;

                fh << batch_id         << "\t"
                   << l->loc_id        << "\t"
                   << it->first        << "\t"
                   << l->bp + 1        << "\t"
                   << pop.name         << "\t"
                   << (int) l->alleles << "\t"
                   << l->hap_cnt       << "\t"
                   << l->stat[0]       << "\t"
                   << l->smoothed[0]   << "\t"
                   << l->bs[0]         << "\t"
                   << l->stat[1]       << "\t"
                   << l->smoothed[1]   << "\t"
                   << l->bs[1]         << "\t"
                   << l->hap_str       << "\n";
            }

            for (uint k = 0; k < locstats.size(); k++)
                delete locstats[k];
        }

        if (bootstrap_div)
            delete bs;
    }

    if (kernel_smoothed && loci_ordered) {
        delete ks;
        delete ord;
    }

    fh.close();

    return 0;
}

int
nuc_substitution_dist(map<string, int> &hap_index, double **hdists)
{
    vector<string> haplotypes;
    map<string, int>::iterator it;
    uint i, j;

    for (it = hap_index.begin(); it != hap_index.end(); it++)
        haplotypes.push_back(it->first);

    const char *p, *q;
    double dist;

    for (i = 0; i < haplotypes.size(); i++) {
        for (j = i; j < haplotypes.size(); j++) {

            dist = 0.0;
            p    = haplotypes[i].c_str();
            q    = haplotypes[j].c_str();

            while (*p != '\0' && *q != '\0') {
                if (*p != *q) dist++;
                p++;
                q++;
            }

            hdists[i][j] = dist;
            hdists[j][i] = dist;
        }
    }

    // //
    // // Print the distance matrix.
    // //
    // cerr << "  ";
    // for (hit = loc_hap_index.begin(); hit != loc_hap_index.end(); hit++)
    //  cerr << "\t" << hit->first;
    // cerr << "\n";
    // for (hit = loc_hap_index.begin(); hit != loc_hap_index.end(); hit++) {
    //  cerr << "  " << hit->first;
    //  for (hit_2 = loc_hap_index.begin(); hit_2 != loc_hap_index.end(); hit_2++)
    //      cerr << "\t" << hdists[hit->second][hit_2->second];
    //  cerr << "\n";
    // }
    // cerr << "\n";

    return 0;
}

int
nuc_substitution_identity(map<string, int> &hap_index, double **hdists)
{
    vector<string> haplotypes;
    map<string, int>::iterator it;
    uint i, j;

    for (it = hap_index.begin(); it != hap_index.end(); it++)
        haplotypes.push_back(it->first);

    double dist;

    for (i = 0; i < haplotypes.size(); i++) {
        for (j = i; j < haplotypes.size(); j++) {

            if (haplotypes[i] == haplotypes[j])
                dist = 0.0;
            else
                dist = 1.0;

            hdists[i][j] = dist;
            hdists[j][i] = dist;
        }
    }

    return 0;
}

int
nuc_substitution_identity_max(map<string, int> &hap_index, double **hdists)
{
    vector<string> haplotypes;
    map<string, int>::iterator it;
    uint i, j;

    for (it = hap_index.begin(); it != hap_index.end(); it++)
        haplotypes.push_back(it->first);

    for (i = 0; i < haplotypes.size(); i++) {
        for (j = i; j < haplotypes.size(); j++) {
            hdists[i][j] = 1.0;
            hdists[j][i] = 1.0;
        }
    }

    return 0;
}

int
calculate_haplotype_divergence(map<int, CSLocus *> &catalog, PopMap<CSLocus> *pmap, PopSum<CSLocus> *psum)
{
    map<string, vector<CSLocus *> >::const_iterator it;

    if (bootstrap_phist)
        cerr << "Calculating halotype F statistics across all populations/groups and bootstrap resampling...\n";
    else
        cerr << "Calculating haplotype F statistics across all populations/groups...\n";

    //
    // Create a list of all the populations we have.
    //
    vector<int> pop_ids;
    for (size_t i=0; i<mpopi.pops().size(); ++i)
        pop_ids.push_back(i);

    //
    // Instantiate the kernel smoothing object and associated ordering object if requested.
    //
    KSmooth<HapStat>     *ks=NULL;
    OHaplotypes<HapStat> *ord=NULL;
    Bootstrap<HapStat>   *bs=NULL;
    if (kernel_smoothed && loci_ordered) {
        ks  = new KSmooth<HapStat>(5);
        ord = new OHaplotypes<HapStat>();
    }

    if (bootstrap_phist)
        bs = new Bootstrap<HapStat>(5);

    map<string, vector<HapStat *> > genome_hapstats;

    uint cnt = 0;
    for (it = pmap->ordered_loci().begin(); it != pmap->ordered_loci().end(); it++) {
        string chr = it->first;

        cerr << "  Generating haplotype F statistics for " << chr << "...";

        map<uint, uint>    hapstats_key;
        vector<HapStat *> &hapstats = genome_hapstats[chr];
        ord->order(hapstats, hapstats_key, it->second);

        #pragma omp parallel
        {
            const CSLocus  *loc;
            LocSum  **s;
            Datum   **d;
            HapStat  *h;

            #pragma omp for schedule(dynamic, 1) reduction(+:cnt)
            for (uint pos = 0; pos < it->second.size(); pos++) {
                loc = it->second[pos];
                s   = psum->locus(loc->id);
                d   = pmap->locus(loc->id);

                if (loc->snps.size() == 0)
                    continue;

                //
                // If this locus only appears in one population or there is only a single haplotype,
                // do not calculate haplotype F stats.
                //
                if (fixed_locus(d, pop_ids))
                    continue;

                cnt++;
                // cerr << "Processing locus " << loc->id << "\n";

                h = haplotype_amova(d, s, pop_ids);

                if (h != NULL) {
                    h->stat[4] = haplotype_d_est(d, s, pop_ids);

                    h->loc_id = loc->id;
                    h->bp     = loc->sort_bp();
                    hapstats[hapstats_key[h->bp]] = h;
                }
            }
        }

        if (bootstrap_phist)
            bs->add_data(hapstats);

        cerr << "done.\n";

        //
        // Calculate kernel-smoothed Fst values.
        //
        if (kernel_smoothed && loci_ordered) {
            cerr << "  Generating kernel-smoothed haplotype F statistics for " << it->first << "...";
            ks->smooth(hapstats);
            cerr << "done.\n";
        }
    }

    if (bootstrap_phist) {
        for (it = pmap->ordered_loci().begin(); it != pmap->ordered_loci().end(); it++)
            bs->execute(genome_hapstats[it->first]);
    }

    cerr << "done.\n";

    if (kernel_smoothed && loci_ordered) {
        delete ks;
        delete ord;
    }

    if (bootstrap_phist)
        delete bs;

    cerr << "Writing haplotype F statistics... ";

    string file = out_path + out_prefix + ".phistats.tsv";

    ofstream fh(file.c_str(), ofstream::out);
    if (fh.fail()) {
        cerr << "Error opening haplotype Phi_st file '" << file << "'\n";
        exit(1);
    }
    fh.precision(fieldw);
    fh.setf(std::ios::fixed);

    //
    // Write the population members.
    //
    for (auto& pop : mpopi.pops()) {
        fh << "# Population " << pop.name << "\t";
        for (size_t k = pop.first_sample; k <= pop.last_sample; k++) {
            fh << mpopi.samples()[k].name;
            if (k < pop.last_sample)
                fh << ",";
        }
        fh << "\n";
    }

    //
    // Write the group members.
    //
    for (auto& group : mpopi.groups()) {
        fh << "# Group " << group.name << "\t";
        for (size_t i_pop : group.pops) {
            fh << mpopi.pops()[i_pop].name;
            if (i_pop != group.pops.back())
                fh << ",";
        }
        fh << "\n";
    }

    fh << "# Batch ID " << "\t"
       << "Locus ID"    << "\t"
       << "Chr"         << "\t"
       << "BP"          << "\t"
       << "PopCnt"      << "\t";
    if (log_fst_comp)
        fh << "SSD(WP)"     << "\t"
           << "SSD(AP/WG)"  << "\t"
           << "SSD(AG)"     << "\t"
           << "SSD(TOTAL)"  << "\t"
           << "MSD(WP)"     << "\t"
           << "MSD(AP/WG)"  << "\t"
           << "MSD(AG)"     << "\t"
           << "MSD(TOTAL)"  << "\t"
           << "n"           << "\t"
           << "n'"          << "\t"
           << "n''"         << "\t"
           << "Sigma2_a"    << "\t"
           << "Sigma2_b"    << "\t"
           << "Sigma2_c"    << "\t"
           << "Sigma_Total" << "\t";
    fh << "phi_st"          << "\t"
       << "Smoothed Phi_st" << "\t"
       << "Smoothed Phi_st P-value" << "\t"
       << "Phi_ct"          << "\t"
       << "Smoothed Phi_ct" << "\t"
       << "Smoothed Phi_ct P-value" << "\t"
       << "Phi_sc"          << "\t"
       << "Smoothed Phi_sc" << "\t"
       << "Smoothed Phi_sc P-value" << "\t"
       << "Fst'"            << "\t"
       << "Smoothed Fst'"   << "\t"
       << "Smoothed Fst' P-value"   << "\t"
       << "D_est"           << "\t"
       << "Smoothed D_est"  << "\t"
       << "Smoothed D_est P-value"  << "\n";

    for (it = pmap->ordered_loci().begin(); it != pmap->ordered_loci().end(); it++) {
        string chr = it->first;

        vector<HapStat *> &hapstats = genome_hapstats[chr];

        for (uint k = 0; k < hapstats.size(); k++) {
            if (hapstats[k] == NULL) continue;

            fh << batch_id            << "\t"
               << hapstats[k]->loc_id << "\t"
               << chr                 << "\t"
               << hapstats[k]->bp +1  << "\t"
               << hapstats[k]->popcnt << "\t";
            if (log_fst_comp)
                fh << hapstats[k]->comp[0]  << "\t"
                   << hapstats[k]->comp[1]  << "\t"
                   << hapstats[k]->comp[2]  << "\t"
                   << hapstats[k]->comp[3]  << "\t"
                   << hapstats[k]->comp[4]  << "\t"
                   << hapstats[k]->comp[5]  << "\t"
                   << hapstats[k]->comp[6]  << "\t"
                   << hapstats[k]->comp[7]  << "\t"
                   << hapstats[k]->comp[8]  << "\t"
                   << hapstats[k]->comp[9]  << "\t"
                   << hapstats[k]->comp[10] << "\t"
                   << hapstats[k]->comp[11] << "\t"
                   << hapstats[k]->comp[12] << "\t"
                   << hapstats[k]->comp[13] << "\t"
                   << hapstats[k]->comp[14] << "\t";
            fh << hapstats[k]->stat[0]     << "\t"
               << hapstats[k]->smoothed[0] << "\t"
               << hapstats[k]->bs[0]       << "\t"
               << hapstats[k]->stat[1]     << "\t"
               << hapstats[k]->smoothed[1] << "\t"
               << hapstats[k]->bs[1]       << "\t"
               << hapstats[k]->stat[2]     << "\t"
               << hapstats[k]->smoothed[2] << "\t"
               << hapstats[k]->bs[2]       << "\t"
               << hapstats[k]->stat[3]     << "\t"
               << hapstats[k]->smoothed[3] << "\t"
               << hapstats[k]->bs[3]       << "\t"
               << hapstats[k]->stat[4]     << "\t"
               << hapstats[k]->smoothed[4] << "\t"
               << hapstats[k]->bs[4]       << "\n";

            delete hapstats[k];
        }
    }

    fh.close();

    cerr << "wrote " << cnt << " loci to haplotype Phi_st file, '" << file << "'\n";

    return 0;
}

int
calculate_haplotype_divergence_pairwise(map<int, CSLocus *> &catalog, PopMap<CSLocus> *pmap, PopSum<CSLocus> *psum)
{
    map<string, vector<CSLocus *> >::const_iterator it;

    if (bootstrap_phist)
        cerr << "Calculating pairwise halotype F statistics and bootstrap resampling...\n";
    else
        cerr << "Calculating pairwise haplotype F statistics...\n";

    //
    // Assign all individuals to one group for the pairwise calculations.
    //
    vector<int> pop_ids;
    for (size_t i=0; i<mpopi.pops().size(); ++i)
        pop_ids.push_back(i);

    //
    // Instantiate the kernel smoothing object if requested.
    //
    KSmooth<HapStat>     *ks=NULL;
    OHaplotypes<HapStat> *ord=NULL;
    Bootstrap<HapStat>   *bs=NULL;
    if (kernel_smoothed && loci_ordered) {
        ks  = new KSmooth<HapStat>(5);
        ord = new OHaplotypes<HapStat>();
    }

    for (uint i = 0; i < mpopi.pops().size(); ++i) {
        const Pop& pop_i = mpopi.pops()[i];
        for (uint j = i + 1; j < mpopi.pops().size(); ++j) {
            const Pop& pop_j = mpopi.pops()[j];
            vector<int> subpop_ids;
            subpop_ids.push_back(i);
            subpop_ids.push_back(j);

            if (bootstrap_phist)
                bs = new Bootstrap<HapStat>(5);

            map<string, vector<HapStat *> > genome_hapstats;

            cerr << "  Processing populations '" << pop_i.name << "' and '" << pop_j.name << "'\n";

            uint cnt = 0;
            for (it = pmap->ordered_loci().begin(); it != pmap->ordered_loci().end(); it++) {
                string chr = it->first;

                cerr << "    Generating pairwise haplotype F statistics for " << chr << "...";

                map<uint, uint>    hapstats_key;
                vector<HapStat *> &hapstats = genome_hapstats[chr];
                ord->order(hapstats, hapstats_key, it->second);

                #pragma omp parallel
                {
                    const CSLocus  *loc;
                    LocSum  **s;
                    Datum   **d;
                    HapStat  *h;

                    #pragma omp for schedule(dynamic, 1) reduction(+:cnt)
                    for (uint pos = 0; pos < it->second.size(); pos++) {
                        loc = it->second[pos];
                        s   = psum->locus(loc->id);
                        d   = pmap->locus(loc->id);

                        if (loc->snps.size() == 0)
                            continue;

                        //
                        // If this locus only appears in one population or there is only a single haplotype,
                        // do not calculate haplotype F stats.
                        //
                        if (fixed_locus(d, subpop_ids))
                            continue;

                        cnt++;
                        // cerr << "Processing locus " << loc->id << "\n";

                        h = haplotype_amova(d, s, subpop_ids);

                        if (h != NULL) {
                            h->stat[4] = haplotype_d_est(d, s, subpop_ids);

                            h->loc_id = loc->id;
                            h->bp     = loc->sort_bp();
                            hapstats[hapstats_key[h->bp]] = h;
                        }
                    }
                }

                if (bootstrap_phist)
                    bs->add_data(hapstats);

                cerr << "done.\n";

                //
                // Calculate kernel-smoothed Fst values.
                //
                if (kernel_smoothed && loci_ordered) {
                    cerr << "    Generating kernel-smoothed Phi_st for " << it->first << "...";
                    ks->smooth(hapstats);
                    cerr << "done.\n";
                }
            }

            if (bootstrap_phist) {
                for (it = pmap->ordered_loci().begin(); it != pmap->ordered_loci().end(); it++)
                    bs->execute(genome_hapstats[it->first]);
            }

            cerr << "done.\n";

            if (bootstrap_phist)
                delete bs;

            cerr << "Writing haplotype F statistics... ";

            string file = out_path + out_prefix + ".phistats_" + pop_i.name + "-" + pop_j.name + ".tsv";

            ofstream fh(file.c_str(), ofstream::out);
            if (fh.fail()) {
                cerr << "Error opening haplotype Phi_st file '" << file << "'\n";
                exit(1);
            }
            fh.precision(fieldw);
            fh.setf(std::ios::fixed);

            //
            // Write the population members.
            //
            for (int k : subpop_ids) {
                const Pop& pop_k = mpopi.pops()[k]; // This is [pop_i], then [pop_j].
                fh << "# Population " << pop_k.name << "\t";
                for (size_t n = pop_k.first_sample; n <= pop_k.last_sample; n++) {
                    fh << mpopi.samples()[n].name;
                    if (n < pop_k.last_sample)
                        fh << ",";
                }
                fh << "\n";
            }

            fh << "# Batch ID " << "\t"
               << "Locus ID"    << "\t"
               << "Pop 1 ID"   << "\t"
               << "Pop 2 ID"   << "\t"
               << "Chr"         << "\t"
               << "BP"          << "\t";
            if (log_fst_comp)
                fh << "SSD(WP)"     << "\t"
                   << "SSD(AP/WG)"  << "\t"
                   << "SSD(AG)"     << "\t"
                   << "SSD(TOTAL)"  << "\t"
                   << "MSD(WP)"     << "\t"
                   << "MSD(AP/WG)"  << "\t"
                   << "MSD(AG)"     << "\t"
                   << "MSD(TOTAL)"  << "\t"
                   << "n"           << "\t"
                   << "n'"          << "\t"
                   << "n''"         << "\t"
                   << "Sigma2_a"    << "\t"
                   << "Sigma2_b"    << "\t"
                   << "Sigma2_c"    << "\t"
                   << "Sigma_Total" << "\t";
            fh << "phi_st"          << "\t"
               << "Smoothed Phi_st" << "\t"
               << "Smoothed Phi_st P-value" << "\t"
               << "Fst'"            << "\t"
               << "Smoothed Fst'"   << "\t"
               << "Smoothed Fst' P-value"   << "\t"
               << "D_est"          << "\t"
               << "Smoothed D_est" << "\t"
               << "Smoothed D_est P-value" << "\n";

            for (it = pmap->ordered_loci().begin(); it != pmap->ordered_loci().end(); it++) {
                string chr = it->first;

                vector<HapStat *> &hapstats = genome_hapstats[chr];

                for (uint k = 0; k < hapstats.size(); k++) {
                    if (hapstats[k] == NULL) continue;

                    fh << batch_id            << "\t"
                       << hapstats[k]->loc_id << "\t"
                       << pop_i.name          << "\t"
                       << pop_j.name          << "\t"
                       << chr                 << "\t"
                       << hapstats[k]->bp +1  << "\t";
                    if (log_fst_comp)
                        fh << hapstats[k]->comp[0]  << "\t"
                           << hapstats[k]->comp[1]  << "\t"
                           << hapstats[k]->comp[2]  << "\t"
                           << hapstats[k]->comp[3]  << "\t"
                           << hapstats[k]->comp[4]  << "\t"
                           << hapstats[k]->comp[5]  << "\t"
                           << hapstats[k]->comp[6]  << "\t"
                           << hapstats[k]->comp[7]  << "\t"
                           << hapstats[k]->comp[8]  << "\t"
                           << hapstats[k]->comp[9]  << "\t"
                           << hapstats[k]->comp[10] << "\t"
                           << hapstats[k]->comp[11] << "\t"
                           << hapstats[k]->comp[12] << "\t"
                           << hapstats[k]->comp[13] << "\t"
                           << hapstats[k]->comp[14] << "\t";
                    fh << hapstats[k]->stat[0]     << "\t"
                       << hapstats[k]->smoothed[0] << "\t"
                       << hapstats[k]->bs[0]       << "\t"
                       << hapstats[k]->stat[3]     << "\t"
                       << hapstats[k]->smoothed[3] << "\t"
                       << hapstats[k]->bs[3]       << "\t"
                       << hapstats[k]->stat[4]     << "\t"
                       << hapstats[k]->smoothed[4] << "\t"
                       << hapstats[k]->bs[4]       << "\n";

                    delete hapstats[k];
                }
            }

            fh.close();

            cerr << "wrote " << cnt << " loci to pairwise haplotype file, '" << file << "'\n";
        }
    }

    if (kernel_smoothed && loci_ordered) {
        delete ks;
        delete ord;
    }

    return 0;
}

bool
fixed_locus(Datum **d, vector<int> &pop_ids)
{
    set<string>               loc_haplotypes;
    map<int, vector<string> > pop_haplotypes;

    for (int pop_id : pop_ids) {
        const Pop& pop = mpopi.pops()[pop_id];
        for (size_t i = pop.first_sample; i <= pop.last_sample; i++) {
            if (d[i] == NULL) continue;

            if (d[i]->obshap.size() > 2) {
                continue;

            } else if (d[i]->obshap.size() == 1) {
                if (!uncalled_haplotype(d[i]->obshap[0])) {
                    loc_haplotypes.insert(d[i]->obshap[0]);
                    pop_haplotypes[pop_id].push_back(d[i]->obshap[0]);
                    pop_haplotypes[pop_id].push_back(d[i]->obshap[0]);
                }
            } else {
                for (uint j = 0; j < d[i]->obshap.size(); j++) {
                    if (!uncalled_haplotype(d[i]->obshap[0])) {
                        loc_haplotypes.insert(d[i]->obshap[j]);
                        pop_haplotypes[pop_id].push_back(d[i]->obshap[j]);
                    }
                }
            }
        }
    }

    uint valid_pops = 0;

    for (int pop_id : pop_ids) {
        if (pop_haplotypes[pop_id].size() > 0)
            valid_pops++;
    }

    //
    // Check that more than one population has data for this locus.
    //
    if (valid_pops <= 1)
        return true;

    //
    // Check that there is more than one haplotype at this locus.
    //
    if (loc_haplotypes.size() == 1)
        return true;

    return false;
}

LocStat *
haplotype_diversity(int start, int end, Datum **d)
{
    map<string, double>::iterator hit;
    vector<string>      haplotypes;
    map<string, double> hap_freq;
    map<string, int>    hap_index;
    double   n              = 0.0;
    double   gene_diversity = 0.0;
    double   hapl_diversity = 0.0;
    LocStat *lstat;

    //
    // Tabulate the haplotypes in this population.
    //
    n = count_haplotypes_at_locus(start, end, d, hap_freq);

    // cerr << "  " << n << " total haplotypes observed.\n";

    //
    // If this haplotype is fixed, don't calculate any statistics.
    //
    if (n == 0)
        return NULL;

    lstat = new LocStat;

    //
    // Store a summary of the haplotype counts to output below.
    //
    stringstream sstr;
    for (hit = hap_freq.begin(); hit != hap_freq.end(); hit++)
        sstr << hit->first << ":" << hit->second  << ";";
    lstat->hap_str = sstr.str().substr(0, sstr.str().length() - 1);

    //
    // Determine an ordering for the haplotypes. Convert haplotype counts into frequencies.
    //
    uint k = 0;
    for (hit = hap_freq.begin(); hit != hap_freq.end(); hit++) {
        hap_index[hit->first] = k;
        haplotypes.push_back(hit->first);
        k++;

        // cerr << "  Haplotype '" << hit->first << "' occured " << hit->second  << " times; ";

        hit->second = hit->second / n;

        // cerr << " frequency of " << hit->second << "%\n";
    }
    //
    // Initialize a two-dimensional array to hold distances between haplotyes.
    //
    double **hdists = new double *[hap_index.size()];
    for (k = 0; k < hap_index.size(); k++) {
        hdists[k] = new double[hap_index.size()];
        memset(hdists[k], 0, hap_index.size());
    }

    //
    // Calculate the distances between haplotypes.
    //
    nuc_substitution_dist(hap_index, hdists);

    //
    // Calculate haplotype diversity, Pi.
    //
    for (uint i = 0; i < haplotypes.size(); i++) {
        for (uint j = 0; j < haplotypes.size(); j++) {
            hapl_diversity +=
                hap_freq[haplotypes[i]] *
                hap_freq[haplotypes[j]] *
                hdists[hap_index[haplotypes[i]]][hap_index[haplotypes[j]]];
        }
    }
    hapl_diversity = (n / (n-1)) * hapl_diversity;

    //
    // Calculate gene diversity.
    //
    for (uint i = 0; i < haplotypes.size(); i++) {
        gene_diversity += hap_freq[haplotypes[i]] * hap_freq[haplotypes[i]];
    }
    gene_diversity = (n / (n - 1)) * (1 - gene_diversity);

    lstat->alleles = n;
    lstat->stat[0] = gene_diversity;
    lstat->stat[1] = hapl_diversity;
    lstat->hap_cnt = haplotypes.size();

    // cerr << "  Population " << pop_id << " has haplotype diversity (pi) of " << s[pop_index]->pi << "\n";

    for (k = 0; k < hap_index.size(); k++)
        delete hdists[k];
    delete [] hdists;

    return lstat;
}

HapStat *
haplotype_amova(Datum **d, LocSum **s, vector<int> &pop_ids)
{
    map<string, int>          loc_hap_index;
    vector<string>            loc_haplotypes;
    map<int, vector<string> > pop_haplotypes;
    map<int, vector<int> >    grp_members;
    vector<int>               grps;

    map<string, int>::iterator hit, hit_2;

    HapStat  *h;

    //
    // Tabulate the occurences of haplotypes at this locus.
    //
    for (int pop_id : pop_ids) {
        const Pop& pop = mpopi.pops()[pop_id];
        for (size_t i = pop.first_sample; i <= pop.last_sample; i++) {
            if (d[i] == NULL) continue;

            if (d[i]->obshap.size() > 2) {
                continue;

            } else if (d[i]->obshap.size() == 1) {
                if(!uncalled_haplotype(d[i]->obshap[0])) {
                    loc_hap_index[d[i]->obshap[0]]++;
                    loc_haplotypes.push_back(d[i]->obshap[0]);
                    loc_haplotypes.push_back(d[i]->obshap[0]);
                    pop_haplotypes[pop_id].push_back(d[i]->obshap[0]);
                    pop_haplotypes[pop_id].push_back(d[i]->obshap[0]);
                }
            } else {
                for (uint j = 0; j < d[i]->obshap.size(); j++) {
                    if(!uncalled_haplotype(d[i]->obshap[0])) {
                        loc_hap_index[d[i]->obshap[j]]++;
                        loc_haplotypes.push_back(d[i]->obshap[j]);
                        pop_haplotypes[pop_id].push_back(d[i]->obshap[j]);
                    }
                }
            }
        }
    }

    //
    // What is the total number of populations that had valid haplotypes.
    //
    double valid_pop_cnt = 0.0;
    for (int pop_id : pop_ids) {
        if (pop_haplotypes[pop_id].size() > 0)
            valid_pop_cnt++;
    }

    //
    // If we filtered a population out at this locus make sure that we still have at least one
    // representative present in each group.
    //
    set<int> uniq_grps;
    for (size_t pop_id=0; pop_id<mpopi.pops().size(); ++pop_id) {
        const Pop& pop = mpopi.pops()[pop_id];
        if (pop_haplotypes.count(pop_id) > 0) {
            uniq_grps.insert(pop.group);
            grp_members[pop.group].push_back(pop_id);
        }
    }
    set<int>::iterator uit;
    for (uit = uniq_grps.begin(); uit != uniq_grps.end(); uit++)
        grps.push_back(*uit);

    if (grps.size() == 0)
        return NULL;

    // cerr << "Groups: ";
    // for (uint i = 0; i < grps.size(); i++)
    //     cerr << grps[i] << ", ";
    // cerr << "\n";
    // for (git = grp_members.begin(); git != grp_members.end(); git++) {
    //     cerr << "Group " << git->first << ": ";
    //     for (uint i = 0; i < git->second.size(); i++)
    //  cerr << git->second[i] << ", ";
    //     cerr << "\n";
    // }

    //
    // Determine an ordering for the haplotypes.
    //
    uint m = 0;
    for (hit = loc_hap_index.begin(); hit != loc_hap_index.end(); hit++) {
        loc_hap_index[hit->first] = m;
        m++;
    }

    //
    // Initialize a two-dimensional array to hold distances between haplotyes.
    //
    double **hdists     = new double *[loc_hap_index.size()];
    double **hdists_max = new double *[loc_hap_index.size()];
    for (uint k = 0; k < loc_hap_index.size(); k++) {
        hdists[k] = new double[loc_hap_index.size()];
        memset(hdists[k], 0, loc_hap_index.size());
        hdists_max[k] = new double[loc_hap_index.size()];
        memset(hdists_max[k], 0, loc_hap_index.size());
    }

    //
    // Calculate the distances between haplotypes.
    //
    nuc_substitution_dist(loc_hap_index, hdists);

    //
    // Calculate the sum of squared distances in each subset: total, within populations, across populations
    // and withing groups, and across groups.
    //
    double ssd_total = amova_ssd_total(loc_haplotypes, loc_hap_index, hdists);
    double ssd_wp    = amova_ssd_wp(grps, grp_members, loc_hap_index, pop_haplotypes, hdists);
    double ssd_ap_wg = amova_ssd_ap_wg(grps, grp_members, loc_hap_index, pop_haplotypes, hdists, hdists);
    double ssd_ag    = grps.size() > 1 ? amova_ssd_ag(grps, grp_members, loc_hap_index, pop_haplotypes, hdists, ssd_total) : 0.0;

    //
    // Calculate n
    //
    double n        = 0.0;
    double n_1      = 0.0;
    double n_2      = 0.0;
    double s_g      = 0.0;
    double tot_cnt  = 0.0;
    double grp_cnt  = 0.0;
    double num_grps = grps.size();
    double a        = 0.0;
    double b        = 0.0;

    for (uint g = 0; g < num_grps; g++) {
        for (uint r = 0; r < grp_members[grps[g]].size(); r++) {
            int pop_id_1 = grp_members[grps[g]][r];
            tot_cnt += (double) pop_haplotypes[pop_id_1].size();
        }
    }
    for (uint g = 0; g < num_grps; g++) {
        grp_cnt = 0.0;
        for (uint r = 0; r < grp_members[grps[g]].size(); r++) {
            int pop_id_1 = grp_members[grps[g]][r];
            grp_cnt += (double) pop_haplotypes[pop_id_1].size();
        }

        a = 0.0;
        for (uint r = 0; r < grp_members[grps[g]].size(); r++) {
            int pop_id_1 = grp_members[grps[g]][r];
            a += (double) (pop_haplotypes[pop_id_1].size() * pop_haplotypes[pop_id_1].size()) / grp_cnt;
        }
        s_g += a;
    }
    n = (tot_cnt - s_g) / (double) (valid_pop_cnt - num_grps);

    // cerr << "  n: "<< n << "\n";

    if (num_grps > 1) {
        //
        // Calculate n'
        //
        a = 0.0;
        for (uint g = 0; g < num_grps; g++) {
            for (uint r = 0; r < grp_members[grps[g]].size(); r++) {
                int pop_id_1 = grp_members[grps[g]][r];
                a += ((double) (pop_haplotypes[pop_id_1].size() * pop_haplotypes[pop_id_1].size()) / tot_cnt);
            }
        }
        n_1 = (s_g - a) / (double) (num_grps - 1.0);

        // cerr << "  n': "<< n_1 << "\n";

        //
        // Calculate n''
        //
        for (uint g = 0; g < num_grps; g++) {
            a = 0.0;
            for (uint r = 0; r < grp_members[grps[g]].size(); r++) {
                int pop_id_1 = grp_members[grps[g]][r];
                a += pop_haplotypes[pop_id_1].size();
            }
            b += ((a * a) / tot_cnt);
        }
        n_2 = (tot_cnt - b) / (double) (num_grps - 1);

        // cerr << "  n'': "<< n_2 << "\n";
    }

    //
    // Calculate the mean square deviations, equal to SSD divided by degrees of freedom.
    //
    double msd_ag    = num_grps > 1 ? ssd_ag / (double) (num_grps - 1) : 0.0;
    double msd_ap_wg = ssd_ap_wg / ((double) (valid_pop_cnt - num_grps));
    double msd_wp    = ssd_wp    / ((double) (loc_haplotypes.size() - valid_pop_cnt));
    double msd_total = ssd_total / ((double) (loc_haplotypes.size() - 1));

    double sigma_c     = msd_wp;
    double sigma_b     = n > 0 ? (msd_ap_wg - sigma_c) / n : 0.0;
    double sigma_a     = 0.0;

    if (grps.size() > 1)
        sigma_a = (msd_ag - sigma_c - (n_1 * sigma_b)) / n_2;

    // Arlequin seems to sum the variance components instead of independently calculating sigma_total: MSD(total) = SSD(total)/degrees.of.freedom
    double sigma_total = sigma_a + sigma_b + sigma_c; // msd_total;

    double phi_st = 0.0;
    double phi_ct = 0.0;
    double phi_sc = 0.0;

    if (grps.size() > 1) {
        phi_st = sigma_total > 0.0 ? (sigma_a + sigma_b) / sigma_total : 0.0;
        phi_ct = sigma_total > 0.0 ?  sigma_a / sigma_total : 0.0;
        phi_sc = (sigma_a + sigma_b) > 0.0 ?  sigma_b / (sigma_b + sigma_c) : 0.0;
    } else {
        phi_st = sigma_total > 0.0 ? sigma_b / sigma_total : 0.0;
    }

    // cerr << "  MSD(AG): " << msd_ag  << "; MSD(AP/WG): " << msd_ap_wg << "; MSD(WP): " << msd_wp  << "; MSD(TOTAL): "  << msd_total   << "\n"
    //      << "  Sigma_a: " << sigma_a << "; Sigma_b: "    << sigma_b   << "; Sigma_c: " << sigma_c << "; Sigma_Total: " << sigma_total << "\n"
    //      << "  Phi_st: "  << phi_st  << "; Phi_ct: "     << phi_ct    << "; Phi_sc: "  << phi_sc  << "\n";

    //
    // Calculate Fst' = Fst / Fst_max
    //
    // First calculate Fst.
    //
    // To calculate Fst instead of Phi_st, we need to reset our distance matrix to return 1 if haplotypes are different, 0 otherwise.
    //
    nuc_substitution_identity(loc_hap_index, hdists);
    ssd_wp    = amova_ssd_wp(grps, grp_members, loc_hap_index, pop_haplotypes, hdists);
    ssd_ap_wg = amova_ssd_ap_wg(grps, grp_members, loc_hap_index, pop_haplotypes, hdists, hdists);
    //
    // Calculate the mean square deviations, equal to SSD divided by degrees of freedom.
    //
    msd_ap_wg   = ssd_ap_wg / ((double) (valid_pop_cnt - num_grps));
    msd_wp      = ssd_wp    / ((double) (loc_haplotypes.size() - valid_pop_cnt));
    sigma_c     = msd_wp;
    sigma_b     = n > 0 ? (msd_ap_wg - sigma_c) / n : 0.0;
    sigma_total = sigma_b + sigma_c;

    double fst = sigma_total > 0.0 ? sigma_b / sigma_total : 0.0;

    //
    // Now calculate Fst_max.
    //
    // Reset our distance matrix to give maximum possible distance between haplotypes
    // and recalculate sum of squared deviations across groups.
    //
    nuc_substitution_identity_max(loc_hap_index, hdists_max);
    ssd_ap_wg = amova_ssd_ap_wg(grps, grp_members, loc_hap_index, pop_haplotypes, hdists, hdists_max);

    //
    // Recalculate the mean square deviations, given maximum divergence between populations.
    //
    msd_ap_wg = ssd_ap_wg / ((double) (valid_pop_cnt - num_grps));
    sigma_b   = n > 0 ? (msd_ap_wg - sigma_c) / n : 0.0;

    double fst_max = sigma_total > 0.0 ? sigma_b / sigma_total : 0.0;
    double fst_1   = fst_max > 0.0     ? fst / fst_max         : 0.0;

    //
    // Cache the results so we can print them in order below, once the parallel code has executed.
    //
    h = new HapStat;
    h->alleles = tot_cnt;
    h->popcnt  = valid_pop_cnt;

    if (log_fst_comp) {
        h->comp = new double[15];
        h->comp[0]  = ssd_wp;
        h->comp[1]  = ssd_ap_wg;
        h->comp[2]  = ssd_ag;
        h->comp[3]  = ssd_total;
        h->comp[4]  = msd_wp;
        h->comp[5]  = msd_ap_wg;
        h->comp[6]  = msd_ag;
        h->comp[7]  = msd_total;
        h->comp[8]  = n;
        h->comp[9]  = n_1;
        h->comp[10] = n_2;
        h->comp[11] = sigma_a;
        h->comp[12] = sigma_b;
        h->comp[13] = sigma_c;
        h->comp[14] = sigma_total;
    }

    h->stat[0] = phi_st;
    h->stat[1] = phi_ct;
    h->stat[2] = phi_sc;
    h->stat[3] = fst_1;

    for (uint k = 0; k < loc_hap_index.size(); k++) {
        delete [] hdists[k];
        delete [] hdists_max[k];
    }
    delete [] hdists;
    delete [] hdists_max;

    return h;
}

double
amova_ssd_total(vector<string> &loc_haplotypes, map<string, int> &loc_hap_index, double **hdists)
{
    //
    // Calculate sum of squared deviations for the total sample, SSD(Total)
    //
    double ssd_total = 0.0;

    for (uint j = 0; j < loc_haplotypes.size(); j++) {
        for (uint k = 0; k < loc_haplotypes.size(); k++) {
            ssd_total += hdists[loc_hap_index[loc_haplotypes[j]]][loc_hap_index[loc_haplotypes[k]]];
            // cerr << j << "\t"
            //   << k << "\t"
            //   << loc_haplotypes[j] << "\t"
            //   << loc_haplotypes[k] << "\t"
            //   << hdists[loc_hap_index[loc_haplotypes[j]]][loc_hap_index[loc_haplotypes[k]]] << "\n";
        }
    }
    ssd_total = (1.0 / (double) (2*loc_haplotypes.size())) * ssd_total;
    // cerr << "  ssd_total: "<< ssd_total << "\n";

    return ssd_total;
}

double
amova_ssd_wp(vector<int> &grps, map<int, vector<int> > &grp_members,
             map<string, int> &loc_hap_index, map<int, vector<string> > &pop_haplotypes,
             double **hdists)
{
    //
    // Calculate the sum of squared deviations within populations, SSD(WP)
    //
    double ssd_wp = 0.0;
    double ssd    = 0.0;
    int    pop_id;

    for (uint g = 0; g < grps.size(); g++) {
        for (uint i = 0; i < grp_members[grps[g]].size(); i++) {
            pop_id = grp_members[grps[g]][i];
            ssd = 0.0;

            for (uint j = 0; j < pop_haplotypes[pop_id].size(); j++) {
                for (uint k = 0; k < pop_haplotypes[pop_id].size(); k++) {
                    ssd += hdists[loc_hap_index[pop_haplotypes[pop_id][j]]][loc_hap_index[pop_haplotypes[pop_id][k]]];
                    // cerr << pop_id << "\t"
                    //   << j << "\t"
                    //   << k << "\t"
                    //   << loc_haplotypes[j] << "\t"
                    //   << loc_haplotypes[k] << "\t"
                    //   << hdists[loc_hap_index[loc_haplotypes[j]]][loc_hap_index[loc_haplotypes[k]]] << "\n";
                }
            }

            if (pop_haplotypes[pop_id].size() > 0)
                ssd_wp += (1.0 / (double) (2*pop_haplotypes[pop_id].size())) * ssd;
        }
    }
    // cerr << "  ssd_wp: "<< ssd_wp << "\n";

    return ssd_wp;
}

double
amova_ssd_ap_wg(vector<int> &grps, map<int, vector<int> > &grp_members,
                map<string, int> &loc_hap_index, map<int, vector<string> > &pop_haplotypes,
                double **hdists_1, double **hdists_2)
{
    //
    // Calculate the sum of squared deviations across populations and within groups, SSD(AP/WG)
    //
    double ssd_ap_wg = 0.0;
    double ssd       = 0.0;
    double ssd_1     = 0.0;
    double ssd_2     = 0.0;
    double den       = 0.0;
    int    pop_id, pop_id_1, pop_id_2;

    for (uint g = 0; g < grps.size(); g++) {

        ssd_1 = 0.0;
        for (uint r = 0; r < grp_members[grps[g]].size(); r++) {
            pop_id_1 = grp_members[grps[g]][r];

            for (uint j = 0; j < pop_haplotypes[pop_id_1].size(); j++) {

                for (uint s = 0; s < grp_members[grps[g]].size(); s++) {
                    pop_id_2 = grp_members[grps[g]][s];

                    for (uint k = 0; k < pop_haplotypes[pop_id_2].size(); k++) {
                        if (pop_id_1 == pop_id_2)
                            ssd_1 += hdists_1[loc_hap_index[pop_haplotypes[pop_id_1][j]]][loc_hap_index[pop_haplotypes[pop_id_2][k]]];
                        else
                            ssd_1 += hdists_2[loc_hap_index[pop_haplotypes[pop_id_1][j]]][loc_hap_index[pop_haplotypes[pop_id_2][k]]];
                    }
                }
            }
        }

        den = 0.0;
        for (uint r = 0; r < grp_members[grps[g]].size(); r++) {
            pop_id_1 = grp_members[grps[g]][r];
            den += 2 * pop_haplotypes[pop_id_1].size();
        }

        ssd_1 = ssd_1 / den;

        ssd_2 = 0.0;
        for (uint r = 0; r < grp_members[grps[g]].size(); r++) {
            pop_id = grp_members[grps[g]][r];
            ssd = 0.0;

            for (uint j = 0; j < pop_haplotypes[pop_id].size(); j++) {
                for (uint k = 0; k < pop_haplotypes[pop_id].size(); k++) {
                    ssd += hdists_1[loc_hap_index[pop_haplotypes[pop_id][j]]][loc_hap_index[pop_haplotypes[pop_id][k]]];
                }
            }

            if (pop_haplotypes[pop_id].size() > 0)
                ssd_2 += (1.0 / (double) (2*pop_haplotypes[pop_id].size())) * ssd;
        }

        ssd_ap_wg += ssd_1 - ssd_2;
    }
    // cerr << "  ssd_ap_wg: "<< ssd_ap_wg << "\n";

    return ssd_ap_wg;
}

double
amova_ssd_ag(vector<int> &grps, map<int, vector<int> > &grp_members,
             map<string, int> &loc_hap_index, map<int, vector<string> > &pop_haplotypes,
             double **hdists, double ssd_total)
{
    //
    // Calculate the sum of squared deviations across groups, SSD(AG)
    //
    int    pop_id_1, pop_id_2;
    double ssd_ag = 0.0;
    double ssd    = 0.0;
    double ssd_1  = 0.0;
    double den    = 0.0;

    for (uint g = 0; g < grps.size(); g++) {
        ssd_1 = 0.0;

        for (uint r = 0; r < grp_members[grps[g]].size(); r++) {
            pop_id_1 = grp_members[grps[g]][r];

            for (uint j = 0; j < pop_haplotypes[pop_id_1].size(); j++) {

                for (uint s = 0; s < grp_members[grps[g]].size(); s++) {
                    pop_id_2 = grp_members[grps[g]][s];

                    for (uint k = 0; k < pop_haplotypes[pop_id_2].size(); k++) {
                        ssd_1 += hdists[loc_hap_index[pop_haplotypes[pop_id_1][j]]][loc_hap_index[pop_haplotypes[pop_id_2][k]]];
                    }
                }
            }
        }

        den = 0.0;
        for (uint r = 0; r < grp_members[grps[g]].size(); r++) {
            pop_id_1 = grp_members[grps[g]][r];
            den += 2 * pop_haplotypes[pop_id_1].size();
        }

        ssd += ssd_1 / den;
    }

    ssd_ag = ssd_total - ssd;

    // cerr << "  ssd_ag: "<< ssd_ag << "\n";

    return ssd_ag;
}

double
haplotype_d_est(Datum **d, LocSum **s, vector<int> &pop_ids)
{
    //
    // Calculate D_est, fixation index, as described by
    //   Bird, et al., 2011, Detecting and measuring genetic differentiation
    //     +-Equation 11
    // and
    //   Jost, 2008, GST and its relatives do not measure differentiation, Molecular Ecology
    //     +- Equation 13, D_est_chao
    //
    map<string, double>            loc_haplotypes;
    map<int, map<string, double> > pop_haplotypes;
    map<int, double>               pop_totals;

    map<string, double>::iterator it;

    uint pop_cnt = pop_ids.size();

    //
    // Tabulate the occurences of haplotypes at this locus.
    //
    for (int pop_id : pop_ids) {
        const Pop& pop = mpopi.pops()[pop_id];
        for (size_t i = pop.first_sample; i <= pop.last_sample; i++) {
            if (d[i] == NULL) {
                continue;
            } else if (d[i]->obshap.size() > 2) {
                continue;
            } else if (d[i]->obshap.size() == 1) {
                loc_haplotypes[d[i]->obshap[0]]         += 2;
                pop_haplotypes[pop_id][d[i]->obshap[0]] += 2;

            } else {
                for (uint j = 0; j < d[i]->obshap.size(); j++) {
                    loc_haplotypes[d[i]->obshap[j]]++;
                    pop_haplotypes[pop_id][d[i]->obshap[j]]++;
                }
            }
        }

        for (it = pop_haplotypes[pop_id].begin(); it != pop_haplotypes[pop_id].end(); it++)
            pop_totals[pop_id] += it->second;
    }

    double x = 0.0;

    for (it = loc_haplotypes.begin(); it != loc_haplotypes.end(); it++) {

        double freq_sum_sq = 0.0;
        double freq_sq_sum = 0.0;
        for (int pop_id : pop_ids) {
            freq_sum_sq += (pop_haplotypes[pop_id][it->first] / pop_totals[pop_id]);
            freq_sq_sum += pow((pop_haplotypes[pop_id][it->first] / pop_totals[pop_id]), 2);
        }
        freq_sum_sq = pow(freq_sum_sq, 2);

        x += (freq_sum_sq - freq_sq_sum) / (pop_cnt - 1);
    }

    double y = 0.0;

    for (it = loc_haplotypes.begin(); it != loc_haplotypes.end(); it++) {
        for (int pop_id : pop_ids) {
            y += (pop_haplotypes[pop_id][it->first] * (pop_haplotypes[pop_id][it->first] - 1)) /
                (pop_totals[pop_id] * (pop_totals[pop_id] - 1));
        }
    }

    double d_est = 1.0 - (x / y);

    return d_est;
}

void log_snps_per_loc_distrib(ostream& log_fh, map<int, CSLocus*>& catalog)
{

    // N.B. The method below gives the same numbers as by counting the SNPs that satisfy
    // `LocTally::nucs[col].allele_cnt >= 2`.
    // This is different than the sumstats file, that uses `LocTally::nucs[col].allele_cnt == 2`;
    // The difference seems to be about 1/1000.

    // Bin loci by number of SNPs.
    map<size_t, size_t> snps_per_loc_distrib;
    for (auto& cloc : catalog)
        ++snps_per_loc_distrib[cloc.second->snps.size()];

    // Fill the gaps in the distribution.
    size_t n_max = snps_per_loc_distrib.rbegin()->first;
    for(size_t i=0; i<n_max; ++i)
        snps_per_loc_distrib.insert({i, 0});

    // Write the distribution.
    log_fh << "# Distribution of the number of SNPs per locus.\n"
              "#n_snps\tn_loci\n";
    for(const auto& n_snps : snps_per_loc_distrib)
        log_fh << n_snps.first << "\t" << n_snps.second << "\n";
}

int
calculate_summary_stats(map<int, CSLocus *> &catalog, PopMap<CSLocus> *pmap, PopSum<CSLocus> *psum)
{
    map<string, vector<CSLocus *> >::const_iterator it;
    const CSLocus  *loc;
    LocSum  **s;
    LocTally *t;
    int       len;
    int       pop_cnt = psum->pop_cnt();

    //
    // Calculate the means for each summary statistic.
    //
    int    *private_cnt;
    double *num_indv_mean, *p_mean, *obs_het_mean, *obs_hom_mean, *exp_het_mean, *exp_hom_mean, *pi_mean, *fis_mean;
    double *num_indv_var,  *p_var,  *obs_het_var,  *obs_hom_var,  *exp_het_var,  *exp_hom_var,  *pi_var,  *fis_var;
    double *num_indv_mean_all, *p_mean_all, *obs_het_mean_all, *obs_hom_mean_all, *exp_het_mean_all, *exp_hom_mean_all, *pi_mean_all, *fis_mean_all;
    double *num_indv_var_all,  *p_var_all,  *obs_het_var_all,  *obs_hom_var_all,  *exp_het_var_all,  *exp_hom_var_all,  *pi_var_all,  *fis_var_all;
    double *n, *n_all, *var_sites;
    private_cnt       = new int[pop_cnt];
    n                 = new double[pop_cnt];
    var_sites         = new double[pop_cnt];
    num_indv_mean     = new double[pop_cnt];
    num_indv_var      = new double[pop_cnt];
    p_mean            = new double[pop_cnt];
    p_var             = new double[pop_cnt];
    obs_het_mean      = new double[pop_cnt];
    obs_het_var       = new double[pop_cnt];
    obs_hom_mean      = new double[pop_cnt];
    obs_hom_var       = new double[pop_cnt];
    exp_het_mean      = new double[pop_cnt];
    exp_het_var       = new double[pop_cnt];
    exp_hom_mean      = new double[pop_cnt];
    exp_hom_var       = new double[pop_cnt];
    pi_mean           = new double[pop_cnt];
    pi_var            = new double[pop_cnt];
    fis_mean          = new double[pop_cnt];
    fis_var           = new double[pop_cnt];

    n_all             = new double[pop_cnt];
    num_indv_mean_all = new double[pop_cnt];
    num_indv_var_all  = new double[pop_cnt];
    p_mean_all        = new double[pop_cnt];
    p_var_all         = new double[pop_cnt];
    obs_het_mean_all  = new double[pop_cnt];
    obs_het_var_all   = new double[pop_cnt];
    obs_hom_mean_all  = new double[pop_cnt];
    obs_hom_var_all   = new double[pop_cnt];
    exp_het_mean_all  = new double[pop_cnt];
    exp_het_var_all   = new double[pop_cnt];
    exp_hom_mean_all  = new double[pop_cnt];
    exp_hom_var_all   = new double[pop_cnt];
    pi_mean_all       = new double[pop_cnt];
    pi_var_all        = new double[pop_cnt];
    fis_mean_all      = new double[pop_cnt];
    fis_var_all       = new double[pop_cnt];

    for (int j = 0; j < pop_cnt; j++) {
        private_cnt[j]   = 0;
        n[j]             = 0.0;
        var_sites[j]     = 0.0;
        num_indv_mean[j] = 0.0;
        num_indv_var[j]  = 0.0;
        p_mean[j]        = 0.0;
        p_var[j]         = 0.0;
        obs_het_mean[j]  = 0.0;
        obs_het_var[j]   = 0.0;
        obs_hom_mean[j]  = 0.0;
        obs_hom_var[j]   = 0.0;
        exp_het_mean[j]  = 0.0;
        exp_het_var[j]   = 0.0;
        exp_hom_mean[j]  = 0.0;
        exp_hom_var[j]   = 0.0;
        pi_mean[j]       = 0.0;
        pi_var[j]        = 0.0;
        fis_mean[j]      = 0.0;
        fis_var[j]       = 0.0;

        n_all[j]             = 0.0;
        num_indv_mean_all[j] = 0.0;
        num_indv_var_all[j]  = 0.0;
        p_mean_all[j]        = 0.0;
        p_var_all[j]         = 0.0;
        obs_het_mean_all[j]  = 0.0;
        obs_het_var_all[j]   = 0.0;
        obs_hom_mean_all[j]  = 0.0;
        obs_hom_var_all[j]   = 0.0;
        exp_het_mean_all[j]  = 0.0;
        exp_het_var_all[j]   = 0.0;
        exp_hom_mean_all[j]  = 0.0;
        exp_hom_var_all[j]   = 0.0;
        pi_mean_all[j]       = 0.0;
        pi_var_all[j]        = 0.0;
        fis_mean_all[j]      = 0.0;
        fis_var_all[j]       = 0.0;
    }

    for (it = pmap->ordered_loci().begin(); it != pmap->ordered_loci().end(); it++) {
        for (uint pos = 0; pos < it->second.size(); pos++) {
            loc = it->second[pos];

            s = psum->locus(loc->id);
            t = psum->locus_tally(loc->id);
            len = strlen(loc->con);

            for (int i = 0; i < len; i++) {
                //
                // Compile private alleles
                //
                if (t->nucs[i].priv_allele >= 0)
                    private_cnt[t->nucs[i].priv_allele]++;

                if (t->nucs[i].allele_cnt == 2) {

                    for (int j = 0; j < pop_cnt; j++) {

                        if (s[j]->nucs[i].num_indv == 0) continue;

                        n[j]++;

                        if (s[j]->nucs[i].pi > 0) var_sites[j]++;

                        num_indv_mean[j] += s[j]->nucs[i].num_indv;
                        p_mean[j]        += s[j]->nucs[i].p;
                        obs_het_mean[j]  += s[j]->nucs[i].obs_het;
                        obs_hom_mean[j]  += s[j]->nucs[i].obs_hom;
                        exp_het_mean[j]  += s[j]->nucs[i].exp_het;
                        exp_hom_mean[j]  += s[j]->nucs[i].exp_hom;
                        pi_mean[j]       += s[j]->nucs[i].stat[0];
                        fis_mean[j]      += s[j]->nucs[i].stat[1] != -7.0 ? s[j]->nucs[i].stat[1] : 0.0;

                        n_all[j]++;
                        num_indv_mean_all[j] += s[j]->nucs[i].num_indv;
                        p_mean_all[j]        += s[j]->nucs[i].p;
                        obs_het_mean_all[j]  += s[j]->nucs[i].obs_het;
                        obs_hom_mean_all[j]  += s[j]->nucs[i].obs_hom;
                        exp_het_mean_all[j]  += s[j]->nucs[i].exp_het;
                        exp_hom_mean_all[j]  += s[j]->nucs[i].exp_hom;
                        pi_mean_all[j]       += s[j]->nucs[i].stat[0];
                        fis_mean_all[j]      += s[j]->nucs[i].stat[1] != -7.0 ? s[j]->nucs[i].stat[1] : 0.0;
                    }

                } else if (t->nucs[i].allele_cnt == 1) {
                    for (int j = 0; j < pop_cnt; j++) {
                        if (s[j]->nucs[i].num_indv == 0) continue;

                        n_all[j]++;
                        num_indv_mean_all[j] += s[j]->nucs[i].num_indv;
                        p_mean_all[j]        += s[j]->nucs[i].p;
                        obs_het_mean_all[j]  += s[j]->nucs[i].obs_het;
                        obs_hom_mean_all[j]  += s[j]->nucs[i].obs_hom;
                        exp_het_mean_all[j]  += s[j]->nucs[i].exp_het;
                        exp_hom_mean_all[j]  += s[j]->nucs[i].exp_hom;
                        pi_mean_all[j]       += s[j]->nucs[i].stat[0];
                        fis_mean_all[j]      += s[j]->nucs[i].stat[1] != -7.0 ? s[j]->nucs[i].stat[1] : 0.0;
                    }
                }
            }
        }
    }

    for (int j = 0; j < pop_cnt; j++) {
        num_indv_mean[j] = num_indv_mean[j] / n[j];
        p_mean[j]        = p_mean[j]        / n[j];
        obs_het_mean[j]  = obs_het_mean[j]  / n[j];
        obs_hom_mean[j]  = obs_hom_mean[j]  / n[j];
        exp_het_mean[j]  = exp_het_mean[j]  / n[j];
        exp_hom_mean[j]  = exp_hom_mean[j]  / n[j];
        pi_mean[j]       = pi_mean[j]       / n[j];
        fis_mean[j]      = fis_mean[j]      / n[j];

        num_indv_mean_all[j] = num_indv_mean_all[j] / n_all[j];
        p_mean_all[j]        = p_mean_all[j]        / n_all[j];
        obs_het_mean_all[j]  = obs_het_mean_all[j]  / n_all[j];
        obs_hom_mean_all[j]  = obs_hom_mean_all[j]  / n_all[j];
        exp_het_mean_all[j]  = exp_het_mean_all[j]  / n_all[j];
        exp_hom_mean_all[j]  = exp_hom_mean_all[j]  / n_all[j];
        pi_mean_all[j]       = pi_mean_all[j]       / n_all[j];
        fis_mean_all[j]      = fis_mean_all[j]      / n_all[j];
    }

    string file = out_path + out_prefix + ".sumstats.tsv";

    ofstream fh(file.c_str(), ofstream::out);
    if (fh.fail()) {
        cerr << "Error opening sumstats file '" << file << "'\n";
        exit(1);
    }
    fh.precision(fieldw);
    fh.setf(std::ios::fixed);

    double p_freq;
    //
    // Write the population members.
    //
    for (auto& pop : mpopi.pops()) {
        fh << "# " << pop.name << "\t";
        for (size_t i = pop.first_sample; i <= pop.last_sample; i++) {
            fh << mpopi.samples()[i].name;
            if (i < pop.last_sample)
                fh << ",";
        }
        fh << "\n";
    }

    cerr << "Writing " << catalog.size() << " loci to summary statistics file, '" << file << "'\n";

    fh << "# Batch ID " << "\t"
       << "Locus ID" << "\t"
       << "Chr"      << "\t"
       << "BP"       << "\t"
       << "Col"      << "\t"
       << "Pop ID"   << "\t"
       << "P Nuc"    << "\t"
       << "Q Nuc"    << "\t"
       << "N"        << "\t"
       << "P"        << "\t"
       << "Obs Het"  << "\t"
       << "Obs Hom"  << "\t"
       << "Exp Het"  << "\t"
       << "Exp Hom"  << "\t"
       << "Pi"       << "\t"
       << "Smoothed Pi"  << "\t"
       << "Smoothed Pi P-value"  << "\t"
       << "Fis"          << "\t"
       << "Smoothed Fis" << "\t"
       << "Smoothed Fis P-value" << "\t"
       << "Private"      << "\n";

    for (it = pmap->ordered_loci().begin(); it != pmap->ordered_loci().end(); it++) {
        for (uint pos = 0; pos < it->second.size(); pos++) {
            loc = it->second[pos];

            s = psum->locus(loc->id);
            t = psum->locus_tally(loc->id);
            len = strlen(loc->con);

            for (int i = 0; i < len; i++) {

                //
                // If this site is fixed in all populations, DON'T output it. If it is variable,
                // or fixed within populations but variable among, DO output it.
                //
                if (t->nucs[i].allele_cnt == 2) {

                    for (int j = 0; j < pop_cnt; j++) {

                        if (s[j]->nucs[i].num_indv == 0)
                            continue;

                        fh << batch_id << "\t"
                           << loc->id << "\t"
                           << loc->loc.chr() << "\t"
                           << loc->sort_bp(i) + 1 << "\t"
                           << i << "\t"
                           << mpopi.pops()[j].name << "\t";

                        //
                        // Output the p and q alleles in the same order in each population.
                        //
                        if (t->nucs[i].p_allele == s[j]->nucs[i].p_nuc) {
                            if (s[j]->nucs[i].q_nuc == 0)
                                fh << s[j]->nucs[i].p_nuc << "\t" << "-";
                            else
                                fh << s[j]->nucs[i].p_nuc << "\t" << s[j]->nucs[i].q_nuc;
                            p_freq = s[j]->nucs[i].p;

                        } else {
                            if (s[j]->nucs[i].q_nuc == 0)
                                fh << "-\t" << s[j]->nucs[i].p_nuc;
                            else
                                fh << s[j]->nucs[i].q_nuc << "\t" << s[j]->nucs[i].p_nuc;
                            p_freq = 1 - s[j]->nucs[i].p;
                        }

                        fh << "\t" << (int) s[j]->nucs[i].num_indv << "\t"
                           << std::setprecision(8)      << p_freq << "\t"
                           << std::setprecision(fieldw) << s[j]->nucs[i].obs_het << "\t"
                           << s[j]->nucs[i].obs_hom   << "\t"
                           << s[j]->nucs[i].exp_het   << "\t"
                           << s[j]->nucs[i].exp_hom   << "\t"
                           << s[j]->nucs[i].stat[0]   << "\t" // Pi
                           << s[j]->nucs[i].smoothed[0] << "\t"  // Smoothed Pi
                           << s[j]->nucs[i].bs[0]       << "\t"  // Pi bootstrapped p-value
                           << (s[j]->nucs[i].stat[1] == -7.0 ? 0.0 : s[j]->nucs[i].stat[1]) << "\t"  // Fis
                           << s[j]->nucs[i].smoothed[1] << "\t"  // Smoothed Fis
                           << s[j]->nucs[i].bs[1]       << "\t"; // Fis bootstrapped p-value.
                        (t->nucs[i].priv_allele == j) ? fh << "1\n" : fh << "0\n";

                        //
                        // Tabulate the residuals to calculate the variance.
                        //
                        num_indv_var[j] += pow((s[j]->nucs[i].num_indv - num_indv_mean[j]), 2);
                        p_var[j]        += pow((s[j]->nucs[i].p        - p_mean[j]),        2);
                        obs_het_var[j]  += pow((s[j]->nucs[i].obs_het  - obs_het_mean[j]),  2);
                        obs_hom_var[j]  += pow((s[j]->nucs[i].obs_hom  - obs_hom_mean[j]),  2);
                        exp_het_var[j]  += pow((s[j]->nucs[i].exp_het  - exp_het_mean[j]),  2);
                        exp_hom_var[j]  += pow((s[j]->nucs[i].exp_hom  - exp_hom_mean[j]),  2);
                        pi_var[j]       += pow((s[j]->nucs[i].stat[0]  - pi_mean[j]),       2);
                        fis_var[j]      += pow((s[j]->nucs[i].stat[1]  - fis_mean[j]),      2);

                        num_indv_var_all[j] += pow((s[j]->nucs[i].num_indv - num_indv_mean_all[j]), 2);
                        p_var_all[j]        += pow((s[j]->nucs[i].p        - p_mean_all[j]),        2);
                        obs_het_var_all[j]  += pow((s[j]->nucs[i].obs_het  - obs_het_mean_all[j]),  2);
                        obs_hom_var_all[j]  += pow((s[j]->nucs[i].obs_hom  - obs_hom_mean_all[j]),  2);
                        exp_het_var_all[j]  += pow((s[j]->nucs[i].exp_het  - exp_het_mean_all[j]),  2);
                        exp_hom_var_all[j]  += pow((s[j]->nucs[i].exp_hom  - exp_hom_mean_all[j]),  2);
                        pi_var_all[j]       += pow((s[j]->nucs[i].stat[0]  - pi_mean_all[j]),       2);
                        fis_var_all[j]      += pow((s[j]->nucs[i].stat[1]  - fis_mean_all[j]),      2);
                    }
                } else if (t->nucs[i].allele_cnt == 1) {
                    for (int j = 0; j < pop_cnt; j++) {
                        if (s[j]->nucs[i].num_indv == 0) continue;

                        num_indv_var_all[j] += pow((s[j]->nucs[i].num_indv - num_indv_mean_all[j]), 2);
                        p_var_all[j]        += pow((s[j]->nucs[i].p        - p_mean_all[j]),        2);
                        obs_het_var_all[j]  += pow((s[j]->nucs[i].obs_het  - obs_het_mean_all[j]),  2);
                        obs_hom_var_all[j]  += pow((s[j]->nucs[i].obs_hom  - obs_hom_mean_all[j]),  2);
                        exp_het_var_all[j]  += pow((s[j]->nucs[i].exp_het  - exp_het_mean_all[j]),  2);
                        exp_hom_var_all[j]  += pow((s[j]->nucs[i].exp_hom  - exp_hom_mean_all[j]),  2);
                        pi_var_all[j]       += pow((s[j]->nucs[i].stat[0]  - pi_mean_all[j]),       2);
                        fis_var_all[j]      += pow((s[j]->nucs[i].stat[1]  - fis_mean_all[j]),      2);
                    }
                }
            }
        }
    }

    //
    // Calculate the variance.
    //
    for (int j = 0; j < pop_cnt; j++) {
        num_indv_var[j] = num_indv_var[j] / (n[j] - 1);
        p_var[j]        = p_var[j] / (n[j] - 1);
        obs_het_var[j]  = obs_het_var[j] / (n[j] - 1);
        obs_hom_var[j]  = obs_hom_var[j] / (n[j] - 1);
        exp_het_var[j]  = exp_het_var[j] / (n[j] - 1);
        exp_hom_var[j]  = exp_hom_var[j] / (n[j] - 1);
        pi_var[j]       = pi_var[j] / (n[j] - 1);
        fis_var[j]      = fis_var[j] / (n[j] - 1);

        num_indv_var_all[j] = num_indv_var_all[j] / (n_all[j] - 1);
        p_var_all[j]        = p_var_all[j] / (n_all[j] - 1);
        obs_het_var_all[j]  = obs_het_var_all[j] / (n_all[j] - 1);
        obs_hom_var_all[j]  = obs_hom_var_all[j] / (n_all[j] - 1);
        exp_het_var_all[j]  = exp_het_var_all[j] / (n_all[j] - 1);
        exp_hom_var_all[j]  = exp_hom_var_all[j] / (n_all[j] - 1);
        pi_var_all[j]       = pi_var_all[j] / (n_all[j] - 1);
        fis_var_all[j]      = fis_var_all[j] / (n_all[j] - 1);
    }

    fh.close();

    file = out_path + out_prefix + ".sumstats_summary.tsv";

    fh.open(file.c_str(), ofstream::out);

    if (fh.fail()) {
        cerr << "Error opening sumstats summary file '" << file << "'\n";
        exit(1);
    }

    //
    // Write out summary statistics of the summary statistics.
    //
    fh << "# Variant positions\n"
       << "# Pop ID\t"
       << "Private\t"
       << "Num Indv\t"
       << "Var\t"
       << "StdErr\t"
       << "P\t"
       << "Var\t"
       << "StdErr\t"
       << "Obs Het\t"
       << "Var\t"
       << "StdErr\t"
       << "Obs Hom\t"
       << "Var\t"
       << "StdErr\t"
       << "Exp Het\t"
       << "Var\t"
       << "StdErr\t"
       << "Exp Hom\t"
       << "Var\t"
       << "StdErr\t"
       << "Pi\t"
       << "Var\t"
       << "StdErr\t"
       << "Fis\t"
       << "Var\t"
       << "StdErr\n";

    double *sq_n     = new double[pop_cnt];
    double *sq_n_all = new double[pop_cnt];

    for (int j = 0; j < pop_cnt; j++) {
        sq_n[j]     = sqrt(n[j]);
        sq_n_all[j] = sqrt(n_all[j]);
    }

    for (int j = 0; j < pop_cnt; j++)
        fh << mpopi.pops()[j].name << "\t"
           << private_cnt[j]         << "\t"
           << num_indv_mean[j]       << "\t"
           << num_indv_var[j]        << "\t"
           << sqrt(num_indv_var[j]) / sq_n[j] << "\t"
           << p_mean[j]              << "\t"
           << p_var[j]               << "\t"
           << sqrt(p_var[j]) / sq_n[j] << "\t"
           << obs_het_mean[j]        << "\t"
           << obs_het_var[j]         << "\t"
           << sqrt(obs_het_var[j]) / sq_n[j] << "\t"
           << obs_hom_mean[j]        << "\t"
           << obs_hom_var[j]         << "\t"
           << sqrt(obs_hom_var[j]) / sq_n[j] << "\t"
           << exp_het_mean[j]        << "\t"
           << exp_het_var[j]         << "\t"
           << sqrt(exp_het_var[j]) / sq_n[j] << "\t"
           << exp_hom_mean[j]        << "\t"
           << exp_hom_var[j]         << "\t"
           << sqrt(exp_hom_var[j]) / sq_n[j] << "\t"
           << pi_mean[j]             << "\t"
           << pi_var[j]              << "\t"
           << sqrt(pi_var[j]) / sq_n[j] << "\t"
           << fis_mean[j]            << "\t"
           << fis_var[j]             << "\t"
           << sqrt(num_indv_var[j]) / sq_n[j] << "\n";

    fh << "# All positions (variant and fixed)\n"
       << "# Pop ID\t"
       << "Private\t"
       << "Sites\t"
       << "Variant Sites\t"
       << "Polymorphic Sites\t"
       << "% Polymorphic Loci\t"
       << "Num Indv\t"
       << "Var\t"
       << "StdErr\t"
       << "P\t"
       << "Var\t"
       << "StdErr\t"
       << "Obs Het\t"
       << "Var\t"
       << "StdErr\t"
       << "Obs Hom\t"
       << "Var\t"
       << "StdErr\t"
       << "Exp Het\t"
       << "Var\t"
       << "StdErr\t"
       << "Exp Hom\t"
       << "Var\t"
       << "StdErr\t"
       << "Pi\t"
       << "Var\t"
       << "StdErr\t"
       << "Fis\t"
       << "Var\t"
       << "StdErr\n";

    for (int j = 0; j < pop_cnt; j++) {
        fh << mpopi.pops()[j].name << "\t"
           << private_cnt[j]             << "\t"
           << n_all[j]                   << "\t"
           << n[j]                       << "\t"
           << var_sites[j]               << "\t"
           << var_sites[j] / n_all[j] * 100 << "\t"
           << num_indv_mean_all[j]       << "\t"
           << num_indv_var_all[j]        << "\t"
           << sqrt(num_indv_var_all[j]) / sq_n_all[j] << "\t"
           << p_mean_all[j]              << "\t"
           << p_var_all[j]               << "\t"
           << sqrt(p_var_all[j]) / sq_n_all[j] << "\t"
           << obs_het_mean_all[j]        << "\t"
           << obs_het_var_all[j]         << "\t"
           << sqrt(obs_het_var_all[j]) / sq_n_all[j] << "\t"
           << obs_hom_mean_all[j]        << "\t"
           << obs_hom_var_all[j]         << "\t"
           << sqrt(obs_hom_var_all[j]) / sq_n_all[j] << "\t"
           << exp_het_mean_all[j]        << "\t"
           << exp_het_var_all[j]         << "\t"
           << sqrt(exp_het_var_all[j]) / sq_n_all[j] << "\t"
           << exp_hom_mean_all[j]        << "\t"
           << exp_hom_var_all[j]         << "\t"
           << sqrt(exp_hom_var_all[j]) / sq_n_all[j] << "\t"
           << pi_mean_all[j]             << "\t"
           << pi_var_all[j]              << "\t"
           << sqrt(pi_var_all[j]) / sq_n_all[j] << "\t"
           << fis_mean_all[j]            << "\t"
           << fis_var_all[j]             << "\t"
           << sqrt(num_indv_var_all[j]) / sq_n_all[j] << "\n";
    }

    delete [] private_cnt;
    delete [] n;
    delete [] var_sites;
    delete [] sq_n;
    delete [] num_indv_mean;
    delete [] num_indv_var;
    delete [] p_mean;
    delete [] p_var;
    delete [] obs_het_mean;
    delete [] obs_het_var;
    delete [] obs_hom_mean;
    delete [] obs_hom_var;
    delete [] exp_het_mean;
    delete [] exp_het_var;
    delete [] exp_hom_mean;
    delete [] exp_hom_var;
    delete [] pi_mean;
    delete [] pi_var;
    delete [] fis_mean;
    delete [] fis_var;

    delete [] n_all;
    delete [] sq_n_all;
    delete [] num_indv_mean_all;
    delete [] num_indv_var_all;
    delete [] p_mean_all;
    delete [] p_var_all;
    delete [] obs_het_mean_all;
    delete [] obs_het_var_all;
    delete [] obs_hom_mean_all;
    delete [] obs_hom_var_all;
    delete [] exp_het_mean_all;
    delete [] exp_het_var_all;
    delete [] exp_hom_mean_all;
    delete [] exp_hom_var_all;
    delete [] pi_mean_all;
    delete [] pi_var_all;
    delete [] fis_mean_all;
    delete [] fis_var_all;

    fh.close();

    return 0;
}

int
write_fst_stats(map<int, CSLocus *> &catalog, PopMap<CSLocus> *pmap, PopSum<CSLocus> *psum, ofstream &log_fh)
{
    //
    // We want to iterate over each pair of populations and calculate Fst at each
    // nucleotide of each locus.
    //
    if (mpopi.pops().size() == 1)
        return 0;

    vector<double> means;

    //
    // Instantiate the kernel smoothing object if requested.
    //
    OPopPair<PopPair>  *ord = new OPopPair<PopPair>(psum, log_fh);
    KSmooth<PopPair>   *ks=NULL;
    Bootstrap<PopPair> *bs=NULL;
    if (kernel_smoothed && loci_ordered) {
        cerr << "Instantiating the kernel smoothing window, using sigma = " << sigma << " with a sliding window size of " << 6 * sigma << "\n";
        ks  = new KSmooth<PopPair>(2);
    }

    for (uint pop_1 = 0; pop_1 < mpopi.pops().size(); pop_1++) {
        const Pop& pop_1p = mpopi.pops()[pop_1];
        for (uint pop_2 = pop_1 + 1; pop_2 < mpopi.pops().size(); pop_2++) {
            const Pop& pop_2p = mpopi.pops()[pop_2];

            double sum = 0.0;
            double cnt = 0.0;

            string file = out_path + out_prefix + ".fst_" + pop_1p.name + "-" + pop_2p.name + ".tsv";
            ofstream fh(file.c_str(), ofstream::out);
            if (fh.fail()) {
                cerr << "Error opening Fst output file '" << file << "'\n";
                exit(1);
            }
            fh.precision(fieldw);
            fh.setf(std::ios::fixed);

            cerr << "Calculating Fst for populations '" << pop_1p.name << "' and '" << pop_2p.name << "' and writing it to file, '" << file << "'\n";

            fh << "# Batch ID" << "\t"
               << "Locus ID"   << "\t"
               << "Pop 1 ID"   << "\t"
               << "Pop 2 ID"   << "\t"
               << "Chr"        << "\t"
               << "BP"         << "\t"
               << "Column"     << "\t"
               << "Overall Pi" << "\t"
               << "Fst"        << "\t"
               << "Fisher's P" << "\t"
               << "Odds Ratio" << "\t"
               << "CI Low"     << "\t"
               << "CI High"    << "\t"
               << "LOD"        << "\t"
               << "Corrected Fst" << "\t"
               << "Smoothed Fst"  << "\t"
               << "AMOVA Fst" << "\t"
               << "Corrected AMOVA Fst" << "\t"
               << "Smoothed AMOVA Fst" << "\t"
               << "Smoothed AMOVA Fst P-value" << "\t"
               << "Window SNP Count";

            //
            // If requested, log Fst component calculations to a file.
            //
            if (log_fst_comp) {
                fh << "\t"
                   << "n_1" << "\t"
                   << "n_2" << "\t"
                   << "tot_alleles" << "\t"
                   << "p_1" << "\t"
                   << "q_1" << "\t"
                   << "p_2" << "\t"
                   << "q_2" << "\t"
                   << "pi_1" << "\t"
                   << "pi_2" << "\t"
                   << "pi_all" << "\t"
                   << "bcoeff_1" << "\t"
                   << "bcoeff_2" << "\t"
                   << "binomial_fst" << "\t"
                   << "p_1_freq" << "\t"
                   << "q_1_freq" << "\t"
                   << "p_2_freq" << "\t"
                   << "q_2_freq" << "\t"
                   << "p_avg_cor" << "\t"
                   << "n_avg_cor" << "\t"
                   << "amova_fst" << "\n";
            } else {
                fh << "\n";
            }

            if (bootstrap_fst)
                bs = new Bootstrap<PopPair>(2);

            map<string, vector<CSLocus *> >::const_iterator it;
            map<string, vector<PopPair *> > genome_pairs;
            // int snp_dist[max_snp_dist] = {0};

            for (it = pmap->ordered_loci().begin(); it != pmap->ordered_loci().end(); it++) {
                string chr = it->first;

                map<uint, uint>    pairs_key;
                vector<PopPair *> &pairs = genome_pairs[chr];

                //
                // Order loci between the two populations and calculate Fst
                //
                ord->order(pairs, pairs_key, it->second, pop_1, pop_2);

                //
                // Apply user-selected correction to the Fst values.
                //
                double correction;
                switch(fst_correction) {
                case p_value:
                    for (uint i = 0; i < pairs.size(); i++) {
                        if (pairs[i] != NULL) {
                            pairs[i]->stat[0] = pairs[i]->fet_p < p_value_cutoff ? pairs[i]->fst : 0;
                            pairs[i]->stat[1] = pairs[i]->fet_p < p_value_cutoff ? pairs[i]->amova_fst : 0;
                        }
                    }
                    break;
                case bonferroni_win:
                    correct_fst_bonferroni_win(pairs);
                    break;
                case bonferroni_gen:
                    correction = p_value_cutoff / catalog.size();
                    for (uint i = 0; i < pairs.size(); i++) {
                        if (pairs[i] != NULL) {
                            pairs[i]->stat[0] = pairs[i]->fet_p < correction ? pairs[i]->fst : 0;
                            pairs[i]->stat[1] = pairs[i]->fet_p < correction ? pairs[i]->amova_fst : 0;
                        }
                    }
                    break;
                case no_correction:
                    for (uint i = 0; i < pairs.size(); i++) {
                        if (pairs[i] != NULL) {
                            pairs[i]->stat[0] = pairs[i]->fst;
                            pairs[i]->stat[1] = pairs[i]->amova_fst;
                        }
                    }
                    break;
                }

                //
                // If bootstrapping is enabled, record all Fst values.
                //
                if (bootstrap_fst)
                    bs->add_data(pairs);

                //
                // Calculate kernel-smoothed Fst values.
                //
                if (kernel_smoothed && loci_ordered) {
                    cerr << "  Generating kernel-smoothed Fst for " << it->first << ".\n";
                    ks->smooth(pairs);
                }
            }

            //
            // If bootstrap resampling method is approximate, generate our single, empirical distribution.
            //
            map<int, vector<double> > approx_fst_dist;
            // if (bootstrap_fst && bootstrap_type == bs_approx)
            //  bootstrap_fst_approximate_dist(fst_samples, allele_depth_samples, weights, snp_dist, approx_fst_dist);

            for (it = pmap->ordered_loci().begin(); it != pmap->ordered_loci().end(); it++) {
                string chr = it->first;
                vector<PopPair *> &pairs = genome_pairs[chr];

                //
                // Bootstrap resample this chromosome.
                //
                if (bootstrap_fst && bootstrap_type == bs_exact) {
                    cerr << "  Bootstrap resampling kernel-smoothed Fst for " << it->first << ".\n";
                    bs->execute(pairs);
                }

                for (uint i = 0; i < pairs.size(); i++) {

                    if (pairs[i] == NULL)
                        continue;

                    //
                    // Calculate Fst P-value from approximate distribution.
                    //
                    // if (bootstrap_fst && bootstrap_type == bs_approx)
                    //     pairs[i]->bs[0] = bootstrap_approximate_pval(pairs[i]->snp_cnt, pairs[i]->stat[0], approx_fst_dist);

                    cnt++;
                    sum += pairs[i]->stat[1]; // Corrected AMOVA Fst

                    fh << batch_id          << "\t"
                       << pairs[i]->loc_id  << "\t"
                       << pop_1p.name       << "\t"
                       << pop_2p.name       << "\t"
                       << chr               << "\t"
                       << pairs[i]->bp +1   << "\t"
                       << pairs[i]->col     << "\t"
                       << pairs[i]->pi      << "\t"
                       << pairs[i]->fst     << "\t"
                       << std::setprecision(9) << pairs[i]->fet_p << "\t"
                       << pairs[i]->fet_or  << "\t"
                       << pairs[i]->ci_low  << "\t"
                       << pairs[i]->ci_high << "\t"
                       << pairs[i]->lod     << "\t"
                       << pairs[i]->stat[0]     << "\t"
                       << pairs[i]->smoothed[0] << "\t"
                       << pairs[i]->amova_fst   << "\t"
                       << pairs[i]->stat[1]     << "\t"
                       << pairs[i]->smoothed[1] << "\t"
                       << pairs[i]->bs[1] << "\t"
                       << pairs[i]->snp_cnt;

                    if (log_fst_comp) {
                        fh << "\t"
                           << pairs[i]->comp[0]   << "\t"
                           << pairs[i]->comp[1]   << "\t"
                           << pairs[i]->comp[2]   << "\t"
                           << pairs[i]->comp[3]   << "\t"
                           << pairs[i]->comp[4]   << "\t"
                           << pairs[i]->comp[5]   << "\t"
                           << pairs[i]->comp[6]   << "\t"
                           << pairs[i]->comp[7]   << "\t"
                           << pairs[i]->comp[8]   << "\t"
                           << pairs[i]->comp[9]   << "\t"
                           << pairs[i]->comp[10]  << "\t"
                           << pairs[i]->comp[11]  << "\t"
                           << pairs[i]->fst       << "\t"
                           << pairs[i]->comp[12]  << "\t"
                           << pairs[i]->comp[13]  << "\t"
                           << pairs[i]->comp[14]  << "\t"
                           << pairs[i]->comp[15]  << "\t"
                           << pairs[i]->comp[16]  << "\t"
                           << pairs[i]->comp[17]  << "\t"
                           << pairs[i]->amova_fst << "\n";
                    } else {
                        fh << "\n";
                    }

                    delete pairs[i];
                }
            }
            cerr << "Pop 1: " << pop_1p.name << "; Pop 2: " << pop_2p.name << "; mean Fst: " << (sum / cnt) << "\n";
            means.push_back(sum / cnt);

            cerr << "Pooled populations '" << pop_1p.name << "' and '" << pop_2p.name << "' contained: " << ord->incompatible_loci << " incompatible loci; "
                 << ord->multiple_loci << " nucleotides covered by more than one RAD locus.\n";
            fh.close();

            if (bootstrap_fst)
                delete bs;
        }
    }

    //
    // Write out the mean Fst measure of each pair of populations.
    //
    string file = out_path + out_prefix + ".fst_summary.tsv";
    ofstream fh(file.c_str(), ofstream::out);

    if (fh.fail()) {
        cerr << "Error opening generic output file '" << file << "'\n";
        exit(1);
    }

    //
    // Write out X-axis header.
    //
    for (auto& pop : mpopi.pops())
        fh << "\t" << pop.name;
    fh << "\n";

    uint n = 0;
    for (uint i = 0; i < mpopi.pops().size() - 1; i++) {
        fh << mpopi.pops()[i].name;

        for (uint k = 0; k <= i; k++)
            fh << "\t";

        for (uint j = i + 1; j < mpopi.pops().size(); j++) {
            fh << "\t" << means[n];
            n++;
        }
        fh << "\n";
    }

    fh.close();

    delete ord;
    if (kernel_smoothed && loci_ordered) {
        delete ks;
    }

    return 0;
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
kernel_smoothed_popstats(map<int, CSLocus *> &catalog, PopMap<CSLocus> *pmap, PopSum<CSLocus> *psum, int pop_id, ofstream &log_fh)
{
    // int snp_dist[max_snp_dist] = {0};
    // int sites_per_snp = 0;
    // int tot_windows = 0;
    map<string, vector<CSLocus *> >::const_iterator it;
    map<string, vector<SumStat *> > genome_sites;

    //
    // Instantiate the kernel smoothing object if requested.
    //
    KSmooth<SumStat>   *ks  = new KSmooth<SumStat>(2);
    OSumStat<SumStat>  *ord = new OSumStat<SumStat>(psum, log_fh);
    Bootstrap<SumStat> *bs = NULL;

    if (bootstrap_pifis)
        bs = new Bootstrap<SumStat>(2);

    for (it = pmap->ordered_loci().begin(); it != pmap->ordered_loci().end(); it++) {
        vector<SumStat *> &sites = genome_sites[it->first];

        ord->order(sites, it->second, pop_id);
        if (bootstrap_pifis) bs->add_data(sites);
    }

    cerr << "    Population '" << mpopi.pops()[pop_id].name << "' contained " << ord->multiple_loci << " nucleotides covered by more than one RAD locus.\n";

    for (it = pmap->ordered_loci().begin(); it != pmap->ordered_loci().end(); it++) {
        if (bootstrap_pifis)
            cerr << "    Smoothing and bootstrapping chromosome " << it->first << "\n";
        else
            cerr << "    Smoothing chromosome " << it->first << "\n";

        vector<SumStat *> &sites = genome_sites[it->first];

        ks->smooth(sites);

        if (bootstrap_pifis && bootstrap_type == bs_exact)
            bs->execute_mixed(sites);
    }

    delete ks;
    delete ord;
    if (bootstrap_pifis) delete bs;

//     //
//     // If bootstrap resampling method is approximate, generate our single, empirical distribution.
//     //
//     map<int, vector<double> > approx_fis_dist;
//     map<int, vector<double> > approx_pi_dist;
//     if (bootstrap && bootstrap_type == bs_approx) {
//      sites_per_snp = sites_per_snp / tot_windows;

//      // cerr << "Sites per snp: " << sites_per_snp << "\n";

//      bootstrap_popstats_approximate_dist(fis_samples, pi_samples, allele_depth_samples,
//                                          weights, snp_dist, sites_per_snp,
//                                          approx_fis_dist, approx_pi_dist);

//      for (it = pmap->ordered_loci().begin(); it != pmap->ordered_loci().end(); it++) {

//          for (uint pos = 0; pos < it->second.size(); pos++) {
//              loc  = it->second[pos];
//              len  = strlen(loc->con);
//              lsum = psum->pop(loc->id, pop_id);

//              for (int k = 0; k < len; k++)
//                  if (lsum->nucs[k].num_indv > 0 && bootstrap && lsum->nucs[k].pi > 0) {
//                      //
//                      // Calculate Fis/Pi p-values from approximate distribution.
//                      //
//                      lsum->nucs[k].wFis_pval = bootstrap_approximate_pval(lsum->nucs[k].snp_cnt, lsum->nucs[k].wFis, approx_fis_dist);
//                      lsum->nucs[k].wPi_pval  = bootstrap_approximate_pval(lsum->nucs[k].snp_cnt, lsum->nucs[k].wPi, approx_pi_dist);
//                  }
//          }
//      }
//     }

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
write_generic(map<int, CSLocus *> &catalog, PopMap<CSLocus> *pmap, bool write_gtypes)
{
    string file = out_path + out_prefix + (write_gtypes ? ".genotypes.tsv" : ".haplotypes.tsv");

    ofstream fh(file.c_str(), ofstream::out);

    if (fh.fail()) {
        cerr << "Error opening generic output file '" << file << "'\n";
        exit(1);
    }

    //
    // Count the number of markers that have enough samples to output.
    //
    map<int, CSLocus *>::iterator it;
    CSLocus *loc;
    int num_loci = catalog.size();

    cerr << "Writing " << num_loci << " loci to " << (write_gtypes ? "genotype" : "observed haplotype") << " file, '" << file << "'\n";

    //
    // Write the header
    //
    fh << "Catalog ID\t";
    if (expand_id)
        fh << "\t";
    if (write_gtypes)
        fh << "Marker\t";
    fh << "Cnt\t";

    for (int i = 0; i < pmap->sample_cnt(); i++) {
        fh << mpopi.samples()[i].name;
        if (i < pmap->sample_cnt() - 1)
            fh << "\t";
    }
    fh << "\n";

    //
    // Output each locus.
    //
    for (it = catalog.begin(); it != catalog.end(); it++) {
        loc = it->second;

        stringstream id;
        loc->annotation.length() > 0 ?
            id << loc->id << "|" << loc->annotation : id << loc->id;

        fh << id.str();

        if (expand_id) {
            if (loc->annotation.length() > 0)
                id << "\t" << loc->id << "\t" << loc->annotation;
            else if (strlen(loc->loc.chr()) > 0)
                id << "\t" << loc->id << "\t" << loc->loc.chr() << "_" << loc->loc.bp +1;
            else
                id << "\t" << loc->id << "\t";
        }

        if (write_gtypes)
            fh << "\t" << loc->marker;

        write_gtypes ? fh << "\t" << loc->gcnt : fh << "\t" << loc->hcnt;

        Datum **d = pmap->locus(loc->id);
        string  obshap;

        for (int i = 0; i < pmap->sample_cnt(); i++) {
            fh << "\t";

            if (d[i] == NULL)
                fh << "-";
            else {
                if (write_gtypes) {
                    fh << d[i]->gtype;
                } else {
                    obshap = "";
                    for (uint j = 0; j < d[i]->obshap.size(); j++)
                        obshap += string(d[i]->obshap[j]) + "/";
                    obshap = obshap.substr(0, obshap.length()-1);
                    fh << obshap;
                }
            }
        }

        fh << "\n";
    }

    fh.close();

    return 0;
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
open_log(ofstream &log_fh)
{
    struct stat out_path_stat;

    if (stat(out_path.substr(0, out_path.length()-1).c_str(), &out_path_stat) == 0) {
        //
        // Path exists, check that it is a directory
        //
        if (!S_ISDIR(out_path_stat.st_mode)) {
            cerr << "Error: '" << out_path.substr(0, out_path.length()-1) << "' is not a directory.\n";
            throw exception();
        }

    } else if (mkdir(out_path.c_str(), ACCESSPERMS) != 0) {
        //
        // Failed to create the directory.
        //
        cerr << "Error: Failed to create directory '" << out_path << "'.\n";
        throw exception();
    }

    string log_path = out_path + out_prefix + ".populations.log";
    log_fh.open(log_path.c_str(), ofstream::out);
    if (log_fh.fail()) {
        cerr << "Error opening log file '" << log_path << "'\n";
        throw exception();
    }
}

void
output_parameters(ostream &fh)
{
    fh << "populations parameters selected:\n";
    if (input_mode == InputMode::vcf)
        fh << "  Input mode: VCF\n";
    fh << "  Fst kernel smoothing: " << (kernel_smoothed == true ? "on" : "off") << "\n"
         << "  Bootstrap resampling: ";
    if (bootstrap)
        fh << "on, " << (bootstrap_type == bs_exact ? "exact; " : "approximate; ") << bootstrap_reps << " reptitions\n";
    else
        fh << "off\n";
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
    fh << "\n";
}

int
parse_command_line(int argc, char* argv[])
{

    while (1) {
        static struct option long_options[] = {
            {"help",           no_argument,       NULL, 'h'},
            {"version",        no_argument,       NULL, 'v'},
            {"verbose",        no_argument,       NULL, 'd'},
            {"sql",            no_argument,       NULL, 's'},
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
            {"batch_id",       required_argument, NULL, 'b'},
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
            {"write_single_snp",  no_argument,       NULL, 'I'},
            {"write_random_snp",  no_argument,       NULL, 'j'},
            {"ordered_export",    no_argument,       NULL, 1002},
            {"kernel_smoothed",   no_argument,       NULL, 'k'},
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
            {"debug_flags",    required_argument, NULL, 1000},
            {0, 0, 0, 0}
        };

        // getopt_long stores the option index here.
        int c = getopt_long(argc, argv, "ACDEFGHJKLNSTUV:YZ123456dghjklnsva:b:c:e:f:i:m:o:p:q:r:t:u:w:B:I:M:O:P:R:Q:W:", long_options, NULL);

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
        case 'b':
            batch_id = is_integer(optarg);
            if (batch_id < 0) {
                cerr << "Batch ID (-b) must be an integer, e.g. 1, 2, 3\n";
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
            kernel_smoothed = true;
            calc_fstats     = true;
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
        case 's':
            sql_out = true;
            break;
        case 1004:
            vcf_out = true;
            break;
        case 'n':
            vcf_haplo_out = true;
            break;
        case 1006:
            fasta_loci_out = true;
            break;
        case 'F':
            fasta_samples_raw_out = true;
            break;
        case 'J':
            fasta_samples_out = true;
            break;
        case 'G':
            genepop_out = true;
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
    }

    if (input_mode == InputMode::stacks || input_mode == InputMode::stacks2) {

        if (pmap_path.empty())
            cerr << "A population map was not specified, all samples will be read from '" << in_path << "' as a single popultaion.\n";

        if (batch_id < 0) {
            vector<int> cat_ids = find_catalogs(in_path);
            if (cat_ids.size() == 1) {
                batch_id = cat_ids[0];
            } else if (cat_ids.empty()) {
                cerr << "Error: Unable to find a catalog in '" << in_path << "'.\n";
                help();
            } else {
                cerr << "Error: Input directory contains several catalogs, please specify -b.\n";
                help();
            }
        }

        if (out_path.empty())
            out_path = in_path;

        out_prefix = string("batch_") + to_string(batch_id);

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
              << "  -P,--in_path: path to the directory containing the Stacks files.\n"
              << "  -V,--in_vcf: path to an input VCF file.\n"
              << "  -O,--out_path: path to a directory where to write the output files. (Required by -V; otherwise defaults to value of -P.)\n"
              << "  -M,--popmap: path to a population map. (Format is 'SAMPLE1 \\t POP1 \\n SAMPLE2 ...'.)\n"
              << "  -t,--threads: number of threads to run in parallel sections of code.\n"
              << "  -b,--batch_id: ID of the catalog to consider (default: guess).\n"
              << "  -s,--sql_out: output a file to import results into an SQL database.\n"
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
              << "  -k,--kernel_smoothed: enable kernel-smoothed Pi, Fis, Fst, Fst', and Phi_st calculations.\n"
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
