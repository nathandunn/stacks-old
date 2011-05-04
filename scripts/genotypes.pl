#!/usr/bin/perl
#
# Copyright 2010, Julian Catchen <jcatchen@uoregon.edu>
#
# This file is part of Stacks.
#
# Stacks is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# Stacks is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with Stacks.  If not, see <http://www.gnu.org/licenses/>.
#

#
# Written by Julian Catchen <jcatchen@uoregon.edu>
#

use strict;
use DBI;

use constant min_hom_seqs   => 5;
use constant min_het_seqs   => 0.05;
use constant max_het_seqs   => 0.1;
use constant stacks_version => "_VERSION_";

use constant true    => 1;
use constant false   => 0;
use constant unknown => -1;

my $mysql_config    = "_PKGDATADIR_" . "sql/mysql.cnf";
my $debug           = 0;
my $db              = "";
my $map_type        = "";
my $out_type        = "joinmap";
my $cache_path      = "_LOCALSTATEDIR_" . "caches";
my $in_path        = "./";
my $progeny_limit   = 4;
my $batch_id        = 0;
my $category        = 0;
my $man_corrections = 0;
my $corrections     = 0;
my $blast_only      = 0;
my $exclude_blast   = 0;
my $expand_id       = 0;
my $impute          = 0;
my $out_file        = "";
my $bl_file         = "";
my $wl_file         = "";
my $sql             = 0;

my $map_types = {'DH'  => \&export_dh_map,
		 'CP'  => \&export_cp_map,
		 'BC1' => \&export_bc1_map,
		 'F2'  => \&export_f2_map};
my $out_types = {'rqtl'    => \&write_rqtl,
		 'joinmap' => \&write_joinmap};

parse_command_line();

#
# Conditionally include the caching module if user specifies automated corrections.
#
require CHI if ($corrections > 0);

#
# Make sure the cache location exists
#
if ($corrections > 0 && !-e $cache_path) {
    print STDERR "Unable to locate the path to store cache files '$cache_path'.\n";
    usage();
}

$cache_path .= "/" . $db . "_batch_" . $batch_id;

#
# Connect to the database and prepare our queries.
#
my (%sth);

prepare_sql_handles(\%sth);

if (!defined($map_types->{$map_type})) {
    print STDERR "Unknown map type.\n";
} else {
    $map_types->{$map_type}->(\%sth, $cache_path);
}

close_sql_handles(\%sth);

sub export_f2_map {
    my ($sth, $cache_path) = @_;

    #
    # We wish to export, according to the JoinMap manual, a locus genotype file (loc-file), 
    # which contains the information of all the loci for a single segregating population.
    # 
    # We are exporting an F2 population type:
    # The result of selfing the F1 of a cross between two fully homozygous diploid parents.
    #
    # Genotype codes for an F2 population, depending on the locus segregation type.
    #
    # Seg. type   Possible genotypes
    # ---------   ------------------
    # <aaxbb>     a, b, h, –
    # <abxcd>     a, b, h, –
    # <abxaa>     a, –
    # <aaxab>     b, –
    # <abxcc>     a, b, – 
    # <ccxab>     b, –

    my (%types, %genotypes, %order, @samples, %sample_ids, %marker_list, %snps, @loci, $locus, $mtype, $key, $file);

    $types{'aa/bb'} = "aaxbb";
    $types{'ab/cd'} = "abxcd";
    $types{'ab/aa'} = "abxaa";
    $types{'aa/ab'} = "aaxab";
    $types{'ab/cc'} = "abxcc";
    $types{'cc/ab'} = "ccxab";

    $genotypes{'aaxbb'} = {'a' => 0, 'b' => 0, 'h' => 0, '-' => 0};
    $genotypes{'abxcd'} = {'a' => 0, 'b' => 0, 'h' => 0, '-' => 0};
    $genotypes{'abxaa'} = {'a' => 0, '-' => 0};
    $genotypes{'aaxab'} = {'b' => 0, '-' => 0};
    $genotypes{'abxcc'} = {'a' => 0, 'b' => 0, '-' => 0};
    $genotypes{'ccxab'} = {'b' => 0, '-' => 0};

    #
    # Determine what order the parents occur in. Order is necessary to properly
    # assign genotypes to alleles in the progeny (JoinMap expects order to 
    # be preserved in markers where a SNP is sex specific, e.g. efxeg).
    #
    determine_parent_order($sth, \%order);

    #
    # Fetch a list of the samples we will output.
    #
    print STDERR "Fetching samples...\n";
    fetch_samples($sth, \@samples, \%sample_ids);

    #
    # Fetch a list of all the loci we will examine
    #
    print STDERR "Fetching loci...\n";
    fetch_loci($sth, \@loci, \%types);

    #
    # Cache the sequencing reads for markers in the parents/progeny.
    #
    my $cache;

    if ($corrections > 0) {
        print STDERR "Caching marker sequences...\n";
        $cache = cache_progeny_seqs(\%sth, $cache_path, \@loci, \%marker_list, \%sample_ids);

	#
	# Fetch the set of SNPs for use in making automated corrections
	#
	print STDERR "Fetching SNPs from the catalog...\n";
	fetch_snps(\%sth, \%snps);
    }

    #
    # Fetch the parents/progeny for each locus
    #
    my $i = 1;
    my $num_loci = scalar(@loci);
    foreach $locus (@loci) {
	#print STDERR "Processing locus $i of $num_loci          \r";

	$locus->{'progeny'} = {};

	call_f2_genotypes($sth, $cache, \%snps, \%sample_ids, \%genotypes, \%order, 
                          $locus->{'id'}, $types{$locus->{'marker'}}, $locus->{'progeny'});

	$i++;
    }
    print STDERR "\n";

    apply_corrected_genotypes($sth, \@loci) if ($man_corrections);

    #
    # Output the results
    #
    $out_types->{$out_type}->(\%types, \@samples, \@loci);

    write_loci_sql(\%types, \@samples, \%sample_ids, \@loci) if ($sql);
}

sub export_dh_map {
    my ($sth, $cache_path) = @_;

    #
    # We wish to export, according to the JoinMap manual, a locus genotype file (loc-file), 
    # which contains the information of all the loci for a single segregating population.
    # 
    # We are exporting a DH population type:
    #   a doubled haploid population: the result of doubling the gametes of a single heterozygous 
    #   diploid individual.
    #
    # Segregation type codes for population type DH, from Joinmap manual:
    #
    # Code      Description
    # -------   -----------
    # <aaxbb>   locus homozygous in both parents, heterozygous between the parents
    #
    # Genotype codes for a CP population, depending on the locus segregation type.
    #
    # Seg. type   Possible genotypes
    # ---------   ------------------
    #     a       the one genotype
    #     b       the other genotype
    #

    my (%types, %genotypes, %order, @samples, %sample_ids, %marker_list, %snps, @loci, $locus, $mtype, $key, $file);

    $types{'ab/--'} = "abx--";
    $types{'--/ab'} = "--xab";

    $genotypes{'abx--'} = {'a' => 0, 'b' => 0, '-' => 0};
    $genotypes{'--xab'} = {'a' => 0, 'b' => 0, '-' => 0};

    #
    # Determine what order the parents occur in. Order is necessary to properly
    # assign genotypes to alleles in the progeny.
    #
    determine_parent_order($sth, \%order);

    #
    # Fetch a list of the samples we will output.
    #
    print STDERR "Fetching samples...\n";
    fetch_samples($sth, \@samples, \%sample_ids);

    #
    # Fetch a list of all the loci we will examine
    #
    print STDERR "Fetching loci...\n";
    fetch_loci($sth, \@loci, \%types);

    #
    # Cache the sequencing reads for markers in the parents/progeny.
    #
    my $cache;

    if ($corrections > 0) {
        print STDERR "Caching marker sequences...\n";
        $cache = cache_progeny_seqs(\%sth, $cache_path, \@loci, \%marker_list, \%sample_ids);

	#
	# Fetch the set of SNPs for use in making automated corrections
	#
	print STDERR "Fetching SNPs from the catalog...\n";
	fetch_snps(\%sth, \%snps);
    }

    #
    # Fetch the parents/progeny for each locus
    #
    my $i = 1;
    my $num_loci = scalar(@loci);
    foreach $locus (@loci) {
	print STDERR "Processing locus $i of $num_loci          \r";

	$locus->{'progeny'} = {};

	call_dh_genotypes($sth, $cache, \%snps, \%sample_ids, \%genotypes, \%order, 
                          $locus->{'id'}, $types{$locus->{'marker'}}, $locus->{'progeny'});

	$i++;
    }
    print STDERR "\n";

    apply_corrected_genotypes($sth, \@loci) if ($man_corrections);

    #
    # Output the results
    #
    $out_types->{$out_type}->(\%types, \@samples, \@loci);

    write_loci_sql(\%types, \@samples, \%sample_ids, \@loci) if ($sql);
}

sub export_bc1_map {
    my ($sth, $cache_path) = @_;

    #
    # We wish to export, according to the JoinMap manual, a locus genotype file (loc-file), 
    # which contains the information of all the loci for a single segregating population.
    # 
    # We are exporting a BC1 population type:
    #   a first generation backcross population: the result of crossing the F1 of a cross between 
    #   two fully homozygous diploid parents to one of the parents.
    #
    # Segregation type codes for population type BC1, from Joinmap manual:
    #
    # Code      Description
    # -------   -----------
    # <aaxbb>   locus homozygous in both parents, heterozygous between the parents
    #
    # Genotype codes for a CP population, depending on the locus segregation type.
    #
    # Seg. type   Possible genotypes
    # ---------   ------------------
    #     a       homozygote or haploid as the first parent
    #     b       homozygote or haploid as the second parent
    #     h       heterozygote (as the F1)
    #

    my (%types, %genotypes, %order, @samples, %sample_ids, %marker_list, %snps, @loci, $locus, $mtype, $key, $file);

    $types{'aa/bb'} = "aaxbb";
    $types{'bb/aa'} = "bbxaa";
    $types{'ab/cc'} = "abxcc";
    $types{'cc/ab'} = "ccxab";

    $genotypes{'aaxbb'} = {'b' => 0, 'h' => 0, '-' => 0};
    $genotypes{'bbxaa'} = {'a' => 0, 'h' => 0, '-' => 0};
    $genotypes{'abxcc'} = {'b' => 0, 'h' => 0, '-' => 0};
    $genotypes{'ccxab'} = {'a' => 0, 'h' => 0, '-' => 0};

    #
    # Determine what order the parents occur in. Order is necessary to properly
    # assign genotypes to alleles in the progeny.
    #
    determine_parent_order($sth, \%order);

    #
    # Fetch a list of the samples we will output.
    #
    print STDERR "Fetching samples...\n";
    fetch_samples($sth, \@samples, \%sample_ids);

    #
    # Fetch a list of all the loci we will examine
    #
    print STDERR "Fetching loci...\n";
    fetch_loci($sth, \@loci, \%types);

    #
    # Cache the sequencing reads for markers in the parents/progeny.
    #
    my $cache;

    if ($corrections > 0) {
        print STDERR "Caching marker sequences...\n";
        $cache = cache_progeny_seqs(\%sth, $cache_path, \@loci, \%marker_list, \%sample_ids);

	#
	# Fetch the set of SNPs for use in making automated corrections
	#
	print STDERR "Fetching SNPs from the catalog...\n";
	fetch_snps(\%sth, \%snps);
    }

    #
    # Fetch the parents/progeny for each locus
    #
    my $i = 1;
    my $num_loci = scalar(@loci);
    foreach $locus (@loci) {
	print STDERR "Processing locus $i of $num_loci          \r";

	$locus->{'progeny'} = {};

	call_bc1_genotypes($sth, $cache, \%snps, \%sample_ids, \%genotypes, \%order, 
                           $locus->{'id'}, $types{$locus->{'marker'}}, $locus->{'progeny'});

	$i++;
    }
    print STDERR "\n";

    apply_corrected_genotypes($sth, \@loci) if ($man_corrections);

    #
    # Output the results
    #
    $out_types->{$out_type}->(\%types, \@samples, \@loci);

    write_loci_sql(\%types, \@samples, \%sample_ids, \@loci) if ($sql);
}

sub export_cp_map {
    my ($sth, $cache_path) = @_;

    #
    # We wish to export, according to the JoinMap manual, a locus genotype file (loc-file), 
    # which contains the information of all the loci for a single segregating population.
    # 
    # We are exporting a CP population type:
    #   a population resulting from a cross between two heterogeneously 
    #   heterozygous and homozygous diploid parents, linkage phases originally 
    #   (possibly) unknown.
    #
    # Segregation type codes for population type CP, from Joinmap manual:
    #
    # Code      Description
    # -------   -----------
    # <abxcd>   locus heterozygous in both parents, four alleles
    # <efxeg>   locus heterozygous in both parents, three alleles
    # <hkxhk>   locus heterozygous in both parents, two alleles
    # <lmxll>   locus heterozygous in the first parent
    # <nnxnp>   locus heterozygous in the second parent
    #
    # Genotype codes for a CP population, depending on the locus segregation type.
    #
    # Seg. type   Possible genotypes
    # ---------   ------------------
    # <abxcd>     ac, ad, bc, bd, ––
    # <efxeg>     ee, ef, eg, fg, ––
    # <hkxhk>     hh, hk, kk, h-, k-, ––
    # <lmxll>     ll, lm, –– 
    # <nnxnp>     nn, np, ––

    my (%types, %genotypes, %order, @samples, %sample_ids, %marker_list, %snps, @loci, $locus, $mtype, $key, $file);

    $types{'ab/--'} = "lmx--";
    $types{'--/ab'} = "--xnp";
    $types{'ab/aa'} = "lmxll";
    $types{'aa/ab'} = "nnxnp";
    $types{'ab/ab'} = "hkxhk";
    $types{'ab/ac'} = "efxeg";
    $types{'ab/cd'} = "abxcd";

    $genotypes{'lmx--'} = {'ll' => 0, 'lm' => 0, '--' => 0};
    $genotypes{'--xnp'} = {'nn' => 0, 'np' => 0, '--' => 0};
    $genotypes{'lmxll'} = {'ll' => 0, 'lm' => 0, '--' => 0};
    $genotypes{'nnxnp'} = {'nn' => 0, 'np' => 0, '--' => 0};
    $genotypes{'hkxhk'} = {'hh' => 0, 'hk' => 0, 'kk' => 0, '--' => 0};
    $genotypes{'efxeg'} = {'ee' => 0, 'ef' => 0, 'eg' => 0, 'fg' => 0, '--' => 0};
    $genotypes{'abxcd'} = {'ac' => 0, 'ad' => 0, 'bc' => 0, 'bd' => 0, '--' => 0};

    #
    # Determine what order the parents occur in. Order is necessary to properly
    # assign genotypes to alleles in the progeny (JoinMap expects order to 
    # be preserved in markers where a SNP is sex specific, e.g. efxeg).
    #
    determine_parent_order($sth, \%order);

    #
    # Fetch a list of the samples we will output.
    #
    print STDERR "Fetching samples...\n";
    fetch_samples($sth, \@samples, \%sample_ids);

    #
    # Fetch a list of all the loci we will examine
    #
    print STDERR "Fetching loci...\n";
    fetch_loci($sth, \@loci, \%types);

    #
    # Cache the sequencing reads for markers in the parents/progeny.
    #
    my $cache;

    if ($corrections > 0) {
        print STDERR "Caching marker sequences...\n";
        $cache = cache_progeny_seqs(\%sth, $cache_path, \@loci, \%marker_list, \%sample_ids);

	#
	# Fetch the set of SNPs for use in making automated corrections
	#
	print STDERR "Fetching SNPs from the catalog...\n";
	fetch_snps(\%sth, \%snps);
    }

    #
    # Fetch the parents/progeny for each locus
    #
    my $i = 1;
    my $num_loci = scalar(@loci);
    foreach $locus (@loci) {
	print STDERR "Processing locus $i of $num_loci          \r";

	$locus->{'progeny'} = {};

	call_cp_genotypes($sth, $cache, \%snps, \%sample_ids, \%genotypes, \%order, 
                          $locus->{'id'}, $types{$locus->{'marker'}}, $locus->{'progeny'});

	$i++;
    }
    print STDERR "\n";

    apply_corrected_genotypes($sth, \@loci) if ($man_corrections);

    #
    # Output the results
    #
    $out_types->{$out_type}->(\%types, \@samples, \@loci);

    write_loci_sql(\%types, \@samples, \%sample_ids, \@loci) if ($sql);
}

sub create_marker_list {
    my ($sth, $loci, $marker_list) = @_;

    my ($row, $locus, $i, $num_loci);

    $i        = 1;
    $num_loci = scalar(@{$loci});

    foreach $locus (@{$loci}) {
	print STDERR "  Processing locus $i of $num_loci          \r" if ($i % 100 == 0);

	$sth->{'marker_list'}->execute($batch_id, $locus->{'id'})
	    or die("Unable to select results from $db.\n");

	while ($row = $sth->{'marker_list'}->fetchrow_hashref()) {

	    if (!defined($marker_list->{$row->{'sample_id'}})) {
		$marker_list->{$row->{'sample_id'}} = {};
	    }

	    $marker_list->{$row->{'sample_id'}}->{$row->{'tag_id'}}++;
	}

	$i++;
    }

    print STDERR "\n";
}

sub cache_progeny_seqs {
    my ($sth, $cache_path, $loci, $marker_list, $samples) = @_;

    #
    # Check if the cache already exists, if so, open and return.
    #
    print STDERR "Checking for the existence of cache: $cache_path\n";
    if (-e $cache_path) {
	my $cache = CHI->new( driver => 'File', root_dir => $cache_path);
	return $cache;
    }

    #
    # Otherwise, populate the new cache.
    #
    my $cache = CHI->new(driver => 'File', root_dir => $cache_path);

    #
    # Create a list of parents/progeny associated with each catalog marker to cache
    #
    print STDERR "Creating marker list of parents/progeny...\n";
    create_marker_list($sth, $loci, $marker_list);

    my ($row, $num_samples, $i, $id, $sample, $key, $tags_fh, $line, $href, $path);

    $num_samples = scalar(keys %{$samples});

    $i = 1;
    foreach $sample (sort keys %{$samples}) {
	print STDERR "  Caching sample $i of $num_samples [$sample, $samples->{$sample}]\n";

	my %stacks;

	$path = $in_path . "/" . $sample . ".tags.tsv";

	open($tags_fh, "<$path") or die("Unable to open file '$path', $!\n");

	#$sth->{'reads'}->execute($samples->{$sample})
	#    or die("Unable to select results from $db.\n");

	#
	# Create a two-dimensional array, each row containing one read.
	#
	# while ($row = $sth->{'reads'}->fetchrow_hashref()) {
	while ($line = <$tags_fh>) {
	    $href = parse_record($line);

	    next if (!defined($marker_list->{$samples->{$sample}}->{$href->{'tag_id'}}));
	    next if ($href->{'rel'} eq "consensus");

	    if (!defined($stacks{$href->{'tag_id'}})) {
		$stacks{$href->{'tag_id'}} = [];
	    }

	    my $aref;
	    @{$aref} = split(//, $href->{'seq'});
	    push(@{$stacks{$href->{'tag_id'}}}, $aref);
	}

	#
	# Write out the stacks of seqeunces to the cache.
	#
	print STDERR "    Writing ", scalar(keys %stacks), " stacks to the file cache.\n";
	foreach $id (keys %stacks) {
	    $key = $samples->{$sample} . "_" . $id;
	    $cache->set($key, $stacks{$id}, "never");
	}

	close($tags_fh);

	$i++;
    }

    print STDERR "\n";

    return $cache;
}

sub parse_record {
    my ($line) = @_;

    my (@parts, $href);

    chomp $line;
    @parts = split(/\t/, $line);

    $href = {};
    $href->{'id'}        = $parts[0];
    $href->{'sample_id'} = $parts[1];
    $href->{'tag_id'}    = $parts[2];
    $href->{'chr'}       = $parts[3];
    $href->{'bp'}        = $parts[4];
    $href->{'rel'}       = $parts[5];
    $href->{'subseq_id'} = $parts[6];
    $href->{'read_id'}   = $parts[7];
    $href->{'seq'}       = $parts[8];
    $href->{'delev'}     = $parts[9];
    $href->{'black'}     = $parts[10];
    $href->{'removed'}   = $parts[11];

    return $href;
}

sub apply_corrected_genotypes {
    my ($sth, $loci) = @_;

    my (%corrections, $locus, $row, $sample);

    print STDERR "Applying manually corrected genotypes to export data...\n";

    #
    # Fetch the manual corrections from the database.
    #
    $sth->{'corr'}->execute($batch_id)
	or die("Unable to select results from $db.\n");

    while ($row = $sth->{'corr'}->fetchrow_hashref()) {
        if (!defined($corrections{$row->{'catalog_id'}})) {
            $corrections{$row->{'catalog_id'}} = {};
        }
        $corrections{$row->{'catalog_id'}}->{$row->{'file'}} = $row->{'genotype'};
    }

    foreach $locus (@{$loci}) {
        next if (!defined($corrections{$locus->{'id'}}));

        foreach $sample (keys %{$corrections{$locus->{'id'}}}) {
            $locus->{'progeny'}->{$sample} = $corrections{$locus->{'id'}}->{$sample};
        }
    }
}

sub write_loci_sql {
    my ($types, $samples, $sample_ids, $loci) = @_;

    my ($file, $out_fh, $locus, $num_loci, $sample);

    $file = $in_path . "/" . "batch_" . $batch_id . ".markers_" . $progeny_limit . ".txt";

    open($out_fh, ">$file") or die("Unable to open output file '$file'; $!\n");

    foreach $locus (@{$loci}) {
	next if (scalar(keys %{$locus->{'progeny'}}) < $progeny_limit);

	foreach $sample (@{$samples}) {
	    print $out_fh 
		"0\t",
		$batch_id, "\t",
		$locus->{'id'}, "\t",
		$sample_ids->{$sample}, "\t";

	    if (defined($locus->{'progeny'}->{$sample})) {
		print $out_fh $locus->{'progeny'}->{$sample}, "\n";
	    } else {
		$map_type eq "CP" ? print $out_fh "--", "\n" : print $out_fh "-", "\n";
	    }
	}
    }

    close($out_fh);
}

sub write_joinmap {
    my ($types, $samples, $loci) = @_;

    my ($out_fh, $bl_fh, $locus, $num_loci, $sample, $key, $progeny_cnt, $id);

    my %blacklist;
    #
    # Read in a blacklist
    #
    load_list($bl_file, \%blacklist);
    print STDERR "Read ", scalar(keys %blacklist), " blacklisted markers.\n" if (length($bl_file) > 0);

    my $pop_name = "batch_" . $batch_id . ".markers_" . $progeny_limit;
    my $file     = $in_path . "/" . $pop_name . ".loc";

    print STDERR "Writing Joinmap output file to '$file'\n";

    open($out_fh, ">$file") or die("Unable to open output file '$file'; $!\n");

    #
    # Count the number of mappable progeny
    #
    $num_loci = 0;

    foreach $locus (@{$loci}) {
	next if (defined($blacklist{$locus->{'id'}}));

        #
        # Count the number of progeny
        #
        $locus->{'progeny_cnt'} = 0;
	foreach $key (keys %{$locus->{'progeny'}}) {
            $locus->{'progeny_cnt'}++ if ($key ne "--" && $key ne "-");
        }
	next if ($locus->{'progeny_cnt'} < $progeny_limit);

	$num_loci++;
    }
    
    #
    # Output the header of the file
    #
    print $out_fh
	"name = $pop_name\n",
	"popt = $map_type\n",
	"nloc = ", $num_loci, "\n",
	"nind = ", scalar(@{$samples}), "\n\n";

    foreach $locus (@{$loci}) {
	next if (defined($blacklist{$locus->{'id'}}));
        next if ($locus->{'progeny_cnt'} < $progeny_limit);

        $id = length($locus->{'ext_id'}) > 0 ? 
            $locus->{'id'} . "|" . $locus->{'ext_id'} :
            $locus->{'id'};

	print $out_fh 
	    $id, "\t";

        if ($expand_id) {
            $id = length($locus->{'ext_id'}) > 0 ? 
                $locus->{'id'} . "\t" . $locus->{'ext_id'} :
                $locus->{'id'} . "\t";

            print $out_fh $id, "\t";
        }

	if ($types->{$locus->{'marker'}} eq "lmx--") {
	    print $out_fh "<lmxll>";
	} elsif ($types->{$locus->{'marker'}} eq "--xnp") {
	    print $out_fh "<nnxnp>";
	} else {
	    print $out_fh "<", $types->{$locus->{'marker'}}, ">";
	}

	foreach $sample (@{$samples}) {
	    print $out_fh "\t";

	    if (defined($locus->{'progeny'}->{$sample})) {
		print $out_fh $locus->{'progeny'}->{$sample};
	    } else {
		$map_type eq "CP" ? print $out_fh "--" : print $out_fh "-";
	    }
	}

	print $out_fh "\n";
    }

    print $out_fh "\nindividual names:\n";

    foreach $sample (@{$samples}) {
	print $out_fh $sample, "\n";
    }

    close($out_fh);
}

sub write_rqtl {
    my ($types, $samples, $loci) = @_;

    my ($out_fh, $bl_fh, $locus, $num_loci, $sample, $key, $progeny_cnt, $id);

    my %blacklist;
    #
    # Read in a blacklist
    #
    load_list($bl_file, \%blacklist);
    print STDERR "Read ", scalar(keys %blacklist), " blacklisted markers.\n" if (length($bl_file) > 0);

    my $pop_name = "batch_" . $batch_id . ".markers_" . $progeny_limit;
    my $file     = $in_path . "/" . $pop_name . ".loc";

    print STDERR "Writing R/QTL output file to '$file'\n";

    open($out_fh, ">$file") or die("Unable to open output file '$file'; $!\n");

    #
    # Count the number of mappable progeny
    #
    $num_loci = 0;

    foreach $locus (@{$loci}) {
	next if (defined($blacklist{$locus->{'id'}}));
        #
        # Count the number of progeny
        #
        $locus->{'progeny_cnt'} = 0;
	foreach $key (keys %{$locus->{'progeny'}}) {
            $locus->{'progeny_cnt'}++ if ($key ne "--" && $key ne "-");
        }
	next if ($locus->{'progeny_cnt'} < $progeny_limit);

	$num_loci++;
    }
    
    #
    # Output the header of the file, followed by the list of markers, one per column
    #
    print $out_fh
	"# Exported: $pop_name\n",
	"# Map Type: $map_type\n",
	"# Num Loci: ", $num_loci, "\n",
	"# Num Samples: ", scalar(@{$samples}), "\n",
	"Individual";

    foreach $locus (@{$loci}) {
	next if (defined($blacklist{$locus->{'id'}}));
        next if ($locus->{'progeny_cnt'} < $progeny_limit);

	print $out_fh ",";

        $id = length($locus->{'ext_id'}) > 0 ? 
            $locus->{'id'} . "|" . $locus->{'ext_id'} :
            $locus->{'id'};

	print $out_fh $id;
    }
    print $out_fh "\n";

    #
    # Output the chromosome (if available) for each marker and then the location
    #
    foreach $locus (@{$loci}) {
	next if (defined($blacklist{$locus->{'id'}}));
        next if ($locus->{'progeny_cnt'} < $progeny_limit);

	print $out_fh ",";

        $id = defined($locus->{'chr'}) ? $locus->{'chr'} : "1";

	print $out_fh $id;
    }
    print $out_fh "\n";

    my $i = 1;
    foreach $locus (@{$loci}) {
	next if (defined($blacklist{$locus->{'id'}}));
        next if ($locus->{'progeny_cnt'} < $progeny_limit);

	print $out_fh ",";

        $id = $locus->{'bp'} > 0 ? $locus->{'bp'} : $i;

	print $out_fh $id;
	$i++;
    }
    print $out_fh "\n";

    #
    # For each sample, print out the genotypes for each marker
    #
    foreach $sample (@{$samples}) {
	print $out_fh $sample;

	foreach $locus (@{$loci}) {
	    next if (defined($blacklist{$locus->{'id'}}));
	    next if ($locus->{'progeny_cnt'} < $progeny_limit);

	    print $out_fh ",";

	    if (defined($locus->{'progeny'}->{$sample})) {
		print $out_fh $locus->{'progeny'}->{$sample};
	    } else {
		$map_type eq "CP" ? print $out_fh "--" : print $out_fh "-";
	    }
	}

	print $out_fh "\n";
    }

    close($out_fh);
}

sub load_list {
    my ($file, $list) = @_;
    
    return if (length($file) == 0);

    my ($fh);

    open($fh, "<$file") or die("Unable to open list '$file'; $!\n");

    while (<$fh>) {
	chomp $_;
	$list->{$_}++;
    }

    close($fh);
}

sub fetch_loci {
    my ($sth, $loci, $genotypes) = @_;

    my ($row, $href, $wl_keys, $wl_fh);
    my %whitelist;

    #
    # Read in a whitelist
    #
    $wl_keys = 0;

    if (length($wl_file) > 0) {
        open($wl_fh, "<$wl_file") or die("Unable to open whitelist file '$wl_file'; $!\n");

        while (<$wl_fh>) {
            chomp $_;
            $whitelist{$_}++;
        }

        $wl_keys = scalar(keys %whitelist);
        print STDERR "Read ", $wl_keys, " whitelisted markers.\n";

        close($wl_fh);
    }

    $sth->{'loci'}->execute($progeny_limit, $batch_id)
	or die("Unable to select results from $db.\n");

    while ($row = $sth->{'loci'}->fetchrow_hashref()) {
        #next if ($row->{'tag_id'} != 3);

	$href = {};
	$href->{'id'}     = $row->{'tag_id'};
	$href->{'marker'} = $row->{'marker'};
	$href->{'ext_id'} = $row->{'external_id'};

        next if (!defined($genotypes->{$href->{'marker'}}));
	next if ($wl_keys && !defined($whitelist{$href->{'id'}}));

	push(@{$loci}, $href);
    }
}

sub fetch_samples {
    my ($sth, $samples, $sample_ids) = @_;

    my ($row);

    $sth->{'samp'}->execute($batch_id)
	or die("Unable to select results from $db.\n");

    while ($row = $sth->{'samp'}->fetchrow_hashref()) {
	if ($row->{'file'} !~ /f?e?male/) {
	    push(@{$samples}, $row->{'file'});
	}
	$sample_ids->{$row->{'file'}} = $row->{'id'};
    }
}

sub determine_parent_order {
    my ($sth, $order) = @_;

    my ($row);

    my @o = ('first', 'second');

    $sth->{'order'}->execute($batch_id)
	or die("Unable to select results from $db.\n");

    my $num_rows = $sth->{'order'}->rows();

    while (@o) {
	$row = $sth->{'order'}->fetchrow_hashref();
	$order->{$row->{'file'}} = shift(@o) if (defined($row));

	print STDERR "Order '$row->{'file'}': ", $order->{$row->{'file'}}, "\n" if ($debug);
	last if ($num_rows == 1);
    }
}

sub fetch_progeny {
    my ($sth, $order, $tag_id, $parents, $progeny) = @_;

    my ($key, $row);

    foreach $key (keys %{$order}) {
        $parents->{$key} = [];
    }

    $sth->{'match'}->execute($batch_id, $tag_id)
	or die("Unable to select results from $db.\n");

    while ($row = $sth->{'match'}->fetchrow_hashref()) {
	if ($row->{'type'} eq "parent") {
	    push(@{$parents->{$row->{'file'}}}, $row->{'allele'});

	} elsif ($row->{'type'} eq "progeny") {
	    if (!defined($progeny->{$row->{'file'}})) {
		$progeny->{$row->{'file'}} = {};
	    }

	    push(@{$progeny->{$row->{'file'}}->{$row->{'tag_id'}}}, $row->{'allele'});
	}
    }
}

sub call_f2_genotypes {
    my ($sth, $cache, $snps, $sample_ids, $genotypes, $order, $tag_id, $marker, $markers) = @_;

    my (@keys, $key, $allele, $m, $gtype, $corrected, %parents, %progeny, %genotype_map);

    my $check_genotypes = {'aaxbb' => \&check_uncalled_snps,
			   'abxcd' => \&check_uncalled_snps,
			   'abxaa' => \&check_uncalled_snps,
			   'aaxab' => \&check_uncalled_snps,
			   'abxcc' => \&check_uncalled_snps,
			   'ccxab' => \&check_uncalled_snps};

    my $create_genotype_map = {'aaxbb' => \&create_genotype_map,
                               'abxcd' => \&create_genotype_map,
                               'abxaa' => \&create_genotype_map,
                               'aaxab' => \&create_genotype_map,
                               'abxcc' => \&create_genotype_map,
                               'ccxab' => \&create_genotype_map};

    my $dictionary = {'aaxbb' => {'a'  => 'a',
				  'ab' => 'h',
				  'b'  => 'b'},
		      'abxcd' => {'a'  => 'a',
				  'ab' => 'a',
				  'b'  => 'a',
				  'c'  => 'b',
				  'cd' => 'b',
				  'd'  => 'b',
				  'ac' => 'h',
				  'ad' => 'h',
				  'bc' => 'h',
				  'bd' => 'h'},
		      'abxaa' => {'a'  => '-',
				  'ab' => '-',
				  'b'  => 'a'},
		      'aaxab' => {'a'  => '-',
				  'ab' => '-',
				  'b'  => 'b'},
		      'abxcc' => {'a'  => 'a',
				  'ab' => 'a',
				  'bb' => 'a',
				  'c'  => 'b',
				  'ac' => '-',
				  'bc' => '-'},
		      'ccxab' => {'aa' => 'b',
				  'ab' => 'b',
				  'bb' => 'b',
				  'c'  => 'a',
				  'ac' => '-',
				  'bc' => '-'}};

    fetch_progeny($sth, $order, $tag_id, \%parents, \%progeny);

    #
    # Create a map between alleles and genotype symbols
    #
    $create_genotype_map->{$marker}->($order, $marker, $tag_id, \%parents, \%progeny, \%genotype_map);

    foreach $key (keys %progeny) {
	my @alleles;

	print STDERR "Examining progeny $key; marker: $marker\n" if ($debug);

	#
	# If there is more than a single tag from a particular progeny, matching this tag in the 
	# catalog, then discard this progeny, since we don't know which tag is correct and the
	# catalog tag is probably overmerged.
	#
	@keys = keys %{$progeny{$key}};
	if (scalar(@keys) > 1) {
	    print STDERR 
		"Discarding progeny $key from catalog tag $tag_id with multiple tag matches: ", 
		join(" ", @keys), "\n" if ($debug);
	    next;
	}

	foreach $allele (@{$progeny{$key}->{$keys[0]}}) {
	    #
	    # Impossible allele encountered.
	    #
	    if (!defined($genotype_map{$allele})) {
		@alleles = ();
		push(@alleles, "-");
		last;
	    }

	    push(@alleles, $genotype_map{$allele});
	}
	$gtype = join("", sort @alleles);

	#
	# Check the genotype, correct any errors.
	#
	$corrected = false;

        if ($corrections > 0) {
            $m = $check_genotypes->{$marker}->($sth, $cache, $snps, $tag_id, $sample_ids->{$key}, $keys[0], \%genotype_map, $gtype);
	    $corrected = true if ($m ne $gtype);
        } else {
	    $m = $gtype;
	}

	$m = defined($dictionary->{$marker}->{$m}) ? $dictionary->{$marker}->{$m} : "-";

	# print STDERR "$key: GENOTYPE m: $m / gtype: $gtype\n";
	if (!defined($genotypes->{$marker}->{$m})) {
	    print STDERR "Warning: illegal genotype encountered in tag $tag_id, progeny $key, '$m'\n" if ($debug);
	} else  {
	    $markers->{$key} = ($corrected == true) ? uc($m) : $m;
	    print STDERR "  $key allele: $markers->{$key}\n" if ($debug);
	}
    }
}

sub call_dh_genotypes {
    my ($sth, $cache, $snps, $sample_ids, $genotypes, $order, $tag_id, $marker, $markers) = @_;

    my (@keys, $key, $allele, $gtype, $m, $corrected, %parents, %progeny, %genotype_map);

    my $check_genotypes = {'abx--' => \&check_uncalled_snps,
                           '--xab' => \&check_uncalled_snps
			   };

    my $create_genotype_map = {'abx--' => \&create_genotype_map,
                               '--xab' => \&create_genotype_map
                               };

    return if (!defined($create_genotype_map->{$marker}));

    my $dictionary = {'abx--' => {'a'  => 'a',
				  'b'  => 'b',
				  '-'  => '-'},
		      '--xab' => {'a'  => 'a',
				  'b'  => 'b',
				  '-'  => '-'}};

    fetch_progeny($sth, $order, $tag_id, \%parents, \%progeny);

    #
    # Create a map between alleles and genotype symbols
    #
    $create_genotype_map->{$marker}->($order, $marker, $tag_id, \%parents, \%progeny, \%genotype_map);

    foreach $key (keys %progeny) {
	my @alleles;

	print STDERR "Examining progeny $key; marker: $marker\n" if ($debug);
	#
	# If there is more than a single tag from a particular progeny, matching this tag in the 
	# catalog, then discard this progeny, since we don't know which tag is correct and the
	# catalog tag is probably overmerged.
	#
	@keys = keys %{$progeny{$key}};
	if (scalar(@keys) > 1) {
	    print STDERR 
		"Discarding progeny $key from catalog tag $tag_id with multiple tag matches: ",
		join(" ", @keys), "\n" if ($debug);
	    next;
	}

	foreach $allele (@{$progeny{$key}->{$keys[0]}}) {
	    #
	    # Impossible allele encountered.
	    #
	    if (!defined($genotype_map{$allele})) {
		@alleles = ();
		push(@alleles, "-");
		last;
	    }

	    push(@alleles, $genotype_map{$allele});
	}
	$gtype = join("", sort @alleles);

	#
	# Check the genotype, correct any errors.
	#
	$corrected = false;

        if ($corrections > 0) {
            $m = $check_genotypes->{$marker}->($sth, $cache, $snps, $tag_id, $sample_ids->{$key}, $keys[0], \%genotype_map, $gtype);
	    $corrected = true if ($m ne $gtype);
        } else {
	    $m = $gtype;
	}

	$m = defined($dictionary->{$marker}->{$m}) ? $dictionary->{$marker}->{$m} : "-";

	if (!defined($genotypes->{$marker}->{$m})) {
	    print STDERR "Warning: illegal genotype encountered in tag $tag_id, progeny $key, '$m'\n" if ($debug);
	} else  {
	    $markers->{$key} = ($corrected == true) ? uc($m) : $m;
	    print STDERR "  $key allele: $markers->{$key}\n" if ($debug);
	}
    }
}

sub call_bc1_genotypes {
    my ($sth, $cache, $snps, $sample_ids, $genotypes, $order, $tag_id, $marker, $markers) = @_;

    my (@keys, $key, $allele, $gtype, $m, $corrected, %parents, %progeny, %genotype_map);

    my $check_genotypes = {'aaxbb' => \&check_uncalled_snps,
                           'bbxaa' => \&check_uncalled_snps,
			   'abxcc' => \&check_uncalled_snps,
                           'ccxab' => \&check_uncalled_snps,
			   };

    my $create_genotype_map = {'aaxbb' => \&create_genotype_map,
                               'bbxaa' => \&create_genotype_map,
                               'abxcc' => \&create_genotype_map,
                               'ccxab' => \&create_genotype_map,
                               };

    return if (!defined($create_genotype_map->{$marker}));

    my $dictionary = {'aaxbb' => {'a'  => 'b',
				  'ab' => 'h',
				  'b'  => 'h'},
		      'bbxaa' => {'a'  => 'h',
				  'ab' => 'h',
				  'b'  => 'a'},
		      'abxcc' => {'ac' => 'h',
				  'bc' => 'h',
				  'ab' => 'b',
				  'a'  => 'b',
				  'b'  => 'b'},
		      'ccxab' => {'ac' => 'h',
				  'bc' => 'h',
				  'ab' => 'a',
				  'a'  => 'a',
				  'b'  => 'a'}};

    fetch_progeny($sth, $order, $tag_id, \%parents, \%progeny);

    #
    # Create a map between alleles and genotype symbols
    #
    $create_genotype_map->{$marker}->($order, $marker, $tag_id, \%parents, \%progeny, \%genotype_map);

    foreach $key (keys %progeny) {
	my @alleles;

	print STDERR "Examining progeny $key; marker: $marker\n" if ($debug);
	#
	# If there is more than a single tag from a particular progeny, matching this tag in the 
	# catalog, then discard this progeny, since we don't know which tag is correct and the
	# catalog tag is probably overmerged.
	#
	@keys = keys %{$progeny{$key}};
	if (scalar(@keys) > 1) {
	    print STDERR 
		"Discarding progeny $key from catalog tag $tag_id with multiple tag matches: ", 
		join(" ", @keys), "\n" if ($debug);
	    next;
	}

	foreach $allele (@{$progeny{$key}->{$keys[0]}}) {
	    #
	    # Impossible allele encountered.
	    #
	    if (!defined($genotype_map{$allele})) {
		@alleles = ();
		push(@alleles, "-");
		last;
	    }

	    push(@alleles, $genotype_map{$allele});
	}
	$gtype = join("", sort @alleles);

	#
	# Check the genotype, correct any errors.
	#
	$corrected = false;

        if ($corrections > 0) {
            $m = $check_genotypes->{$marker}->($sth, $cache, $snps, $tag_id, $sample_ids->{$key}, $keys[0], \%genotype_map, $gtype);
	    $corrected = true if ($m ne $gtype);
        } else {
	    $m = $gtype;
	}

	$m = defined($dictionary->{$marker}->{$m}) ? $dictionary->{$marker}->{$m} : "-";
	
	if (!defined($genotypes->{$marker}->{$m})) {
	    print STDERR "Warning: illegal genotype encountered in tag $tag_id, progeny $key, '$m'\n" if ($debug);
	} else  {
	    $markers->{$key} = ($corrected == true) ? uc($m) : $m;
	    print STDERR "  $key allele: $markers->{$key}\n" if ($debug);
	}
    }
}

sub call_cp_genotypes {
    my ($sth, $cache, $snps, $sample_ids, $genotypes, $order, $tag_id, $marker, $markers) = @_;

    my (@keys, $key, $allele, $gtype, $m, $corrected, %parents, %progeny, %genotype_map);

    my $check_genotypes = {'lmx--' => \&check_uncalled_snps,
			   '--xnp' => \&check_uncalled_snps,
			   'lmxll' => \&check_uncalled_snps,
			   'nnxnp' => \&check_uncalled_snps,
			   'hkxhk' => \&check_uncalled_snps,
			   'efxeg' => \&check_uncalled_snps,
			   'abxcd' => \&check_uncalled_snps
			   };

    my $create_genotype_map = {'lmx--' => \&create_genotype_map,
                               '--xnp' => \&create_genotype_map,
                               'lmxll' => \&create_genotype_map,
                               'nnxnp' => \&create_genotype_map,
                               'hkxhk' => \&create_genotype_map,
                               'efxeg' => \&create_genotype_map,
                               'abxcd' => \&create_genotype_map
                               };
    if ($impute) {
	$create_genotype_map->{'lmxll'} = \&create_imputed_genotype_map;
    }

    my $dictionary = {'lmx--' => {'l'  => 'll',
				  'm'  => 'lm'},
		      '--xnp' => {'n'  => 'nn',
				  'p'  => 'np'},
		      'lmxll' => {'l'  => 'll',
				  'lm' => 'lm'},
		      'nnxnp' => {'n'  => 'nn',
				  'np' => 'np'},
		      'hkxhk' => {'hk' => 'hk',
				  'k'  => 'kk',
				  'h'  => 'hh'},
		      'efxeg' => {'ef' => 'ef',
				  'eg' => 'eg',
				  'fg' => 'fg',
				  'e'  => 'ee'},
		      'abxcd' => {'ac' => 'ac',
				  'ad' => 'ad',
				  'bc' => 'bc',
				  'bd' => 'bd'}};

    fetch_progeny($sth, $order, $tag_id, \%parents, \%progeny);

    #
    # Create a map between alleles and genotype symbols
    #
    $create_genotype_map->{$marker}->($order, $marker, $tag_id, \%parents, \%progeny, \%genotype_map);

    foreach $key (keys %progeny) {
	my @alleles;

	print STDERR "Examining progeny $key; marker: $marker\n" if ($debug);

	#
	# If there is more than a single tag from a particular progeny, matching this tag in the 
	# catalog, then discard this progeny, since we don't know which tag is correct and the
	# catalog tag is probably overmerged.
	#
	@keys = keys %{$progeny{$key}};
	if (scalar(@keys) > 1) {
	    print STDERR 
		"Discarding progeny $key from catalog tag $tag_id with multiple tag matches: ", 
		join(" ", @keys), "\n" if ($debug);
	    next;
	}

	foreach $allele (@{$progeny{$key}->{$keys[0]}}) {
	    #
	    # Impossible allele encountered.
	    #
	    if (!defined($genotype_map{$allele})) {
		@alleles = ();
		push(@alleles, "-", "-");
		last;
	    }

	    push(@alleles, $genotype_map{$allele});
	}
	$gtype = join("", sort @alleles);

	#
	# Check the genotype, correct any errors.
	#
	$corrected = false;

        if ($corrections > 0) {
            $m = $check_genotypes->{$marker}->($sth, $cache, $snps, $tag_id, $sample_ids->{$key}, $keys[0], \%genotype_map, $gtype);
	    $corrected = true if ($m ne $gtype);
        } else {
	    $m = $gtype;
	}

	$m = defined($dictionary->{$marker}->{$m}) ? $dictionary->{$marker}->{$m} : "--";

	if (!defined($genotypes->{$marker}->{$m})) {
	    print STDERR "Warning: illegal genotype encountered in tag $tag_id, progeny $key, '$m'\n" if ($debug);
	} else  {
	    $markers->{$key} = ($corrected == true) ? uc($m) : $m;
	    print STDERR "  $key allele: $markers->{$key}\n" if ($debug);
	}
    }
}

sub create_genotype_map {
    my ($order, $marker, $tag_id, $parents, $progeny, $map) = @_;

    my (%genotypes, %legal_genotypes, %par_specific, @parents, @types, @keys, %alleles, %com_types);
    my ($type, $m, $key, $allele);
    #
    # Create a genotype map. For any set of alleles, this routine will
    # assign each allele to one of the constituent genotypes, e.g. given the 
    # marker type 'aaxbb' and the alleles 'A' from the male, and 'G'
    # from the female, will assign 'G' == 'bb' and 'A'== 'aa'. It assumes that 
    # recombination may have occurred as with an F2, F3 or later cross.
    #
    @parents = keys %{$parents};

    print STDERR "Creating genotype map.\n" if ($debug);

    #
    # First, identify any alleles that are common between the two parents.
    #
    # Record genotypes from first parent
    foreach $type (split(//, substr($marker, 0, 2))) {
	next if ($type eq "-");
        $par_specific{$type}++;
    }
    foreach $type (keys %par_specific) {
	$legal_genotypes{$type}++;
    }
    # Record genotypes from second parent
    foreach $type (split(//, substr($marker, 3, 2))) {
	next if ($type eq "-");
        $par_specific{$type}++;
    }
    foreach $type (keys %par_specific) {
	$legal_genotypes{$type}++;
    }

    # Find the common genotypes
    foreach $type (keys %legal_genotypes) {
	push(@types, $type) if ($legal_genotypes{$type} > 1);
    }
    @types = sort @types;

    foreach $allele (@{$parents->{$parents[0]}}, @{$parents->{$parents[1]}}) {
	$alleles{$allele}++;
    }
    @keys = sort {$alleles{$b} <=> $alleles{$a}} keys %alleles;

    foreach $allele (@keys) {
	if ($alleles{$allele} > 1) {
	    $map->{$allele} = shift @types;
	    $com_types{$map->{$allele}}++;
	    print STDERR "  Assinging common allele '$allele' to genotype '", $map->{$allele}, "'\n" if ($debug);
	}
    }

    #
    # Now, examine the remaining first parent alleles.
    #
    %legal_genotypes = ();
    $key = $order->{$parents[0]} eq "first"  ? $parents[0] : $parents[1];
    $m   = substr($marker, 0, 2);

    foreach $type (split(//, $m)) {
	next if ($type eq "-" || defined($com_types{$type}));
	print STDERR "  Adding $type to genotypes\n" if ($debug);
        $legal_genotypes{$type}++;
    }
    @types = reverse sort keys %legal_genotypes;

    if (scalar(@types) > 0) {
	foreach $allele (@{$parents->{$key}}) {
	    next if (defined($map->{$allele}));
	    $map->{$allele} = shift @types;
	    print STDERR "  Assinging '$allele' to first parent genotype '", $map->{$allele}, "'\n" if ($debug);
	}
    }

    #
    # Finally, repeat in the second parent.
    #
    %legal_genotypes = ();
    $key = $order->{$parents[0]} eq "second"  ? $parents[0] : $parents[1];
    $m   = substr($marker, 3, 2);

    foreach $type (split(//, $m)) {
	next if ($type eq "-" || defined($com_types{$type}));
	print STDERR "  Adding $type to genotypes\n" if ($debug);
        $legal_genotypes{$type}++;
    }
    @types = sort keys %legal_genotypes;

    if (scalar(@types) > 0) {
	foreach $allele (@{$parents->{$key}}) {
	    next if (defined($map->{$allele}));
	    $map->{$allele} = shift @types;
	    print STDERR "  Assinging '$allele' to second parent genotype '", $map->{$allele}, "'\n" if ($debug);
	}
    }
}

sub create_imputed_genotype_map {
    my ($order, $marker, $tag_id, $parents, $progeny, $map) = @_;

    my (@keys, $key, $m, $type, $allele, $uall);

    my (%gtypes, %allcnt, %uniqall);

    #
    # Count up the number of each type of observed allele in the progeny,
    # record whether those genotypes are heterozygous or homozygous.
    #
    foreach $key (keys %{$progeny}) {
	my $alleles;

	print STDERR "Examining progeny $key\n" if ($debug);
	#
	# Discard progeny with more than one locus matched to this catalog tag.
	#
	@keys = keys %{$progeny->{$key}};
	next if (scalar(@keys) > 1);

	$alleles = join("|", sort @{$progeny->{$key}->{$keys[0]}});

	if (!defined($allcnt{$alleles})) {
	    $allcnt{$alleles} = scalar(@{$progeny->{$key}->{$keys[0]}});
	}
	#print STDERR "Adding genotype $alleles\n";

	$gtypes{$alleles}++;

	foreach $allele (@{$progeny->{$key}->{$keys[0]}}) {
	    $uniqall{$allele}++;
	}
    }
 
    #
    # Examine the first parent alleles (the only alleles we have, since 
    # we are imputing the second parent.
    #
    my @parents = keys %{$parents};
    my %legal_genotypes = ();
    $key = $order->{$parents[0]} eq "first"  ? $parents[0] : $parents[1];
    $m   = substr($marker, 0, 2);

    foreach $type (split(//, $m)) {
	#print STDERR "  Adding $type to genotypes\n" if ($debug);
        $legal_genotypes{$type}++;
    }
    my @types = sort keys %legal_genotypes;

    if ($marker eq "lmxll") {
	@keys = sort {$gtypes{$b} <=> $gtypes{$a}} keys %gtypes;
	#
	# Discard heterozygous alleles and find the highest frequency homozygote, 
	# this is the "l" in the "lmxll" marker.
	#
	while ($allcnt{$keys[0]} == 2) {
	    shift @keys;
	}
	$map->{$keys[0]} = shift @types;
	print STDERR "  Assinging '$keys[0]' to first parent genotype '", $map->{$keys[0]}, "'\n" if ($debug);

	foreach $uall (sort {$uniqall{$b} <=> $uniqall{$a}} keys %uniqall) {
	    if ($uall ne $keys[0]) {
		$allele = $uall;
		last;
	    }
	}
	$map->{$allele} = shift @types;
	print STDERR "  Assinging '$allele' to first parent genotype '", $map->{$allele}, "'\n" if ($debug);
    }

}

sub check_uncalled_snps {
    my ($sth, $cache, $snps, $catalog_id, $sample_id, $tag_id, $map, $gtype) = @_;

    #
    # If this locus is already known to be multi-allele, return, we only want 
    # to verify uncalled SNPs.
    #
    return $gtype if (length($gtype) > 1);

    print STDERR "Sample: $sample_id; Tag: $tag_id: Genotype; $gtype\n" if ($debug);

    my ($reads, $row, $key, $homozygous, $rev_gtype, @snps, @nucs, $nuc, $allele, $orig_gtype);

    $orig_gtype = $gtype;
    $key        = $sample_id . "_" . $tag_id;
    $reads      = $cache->get($key);

    if (!defined($reads)) {
	print STDERR "Unable to locate cached object with key $key\n"; 
	return;
    }

    my (@verified_snps, @alleles, @types, $unk);

    $unk = false;

    foreach $row (@{$snps->{$catalog_id}}) {
        $homozygous = check_homozygosity($reads, $row->{'col'}, $row->{'rank_1'}, $row->{'rank_2'});

	if ($homozygous == unknown) {
	    $unk = true;
	} elsif ($homozygous == false) {
	    push(@verified_snps, $row);
	}
    }

    if ($unk == true) {
	return "-";
    } elsif (scalar(@verified_snps) < scalar(@{$snps->{$catalog_id}})) {
	return $gtype;
    }

    #
    # Detect the alleles present from the verified SNPs
    #
    call_alleles(\@verified_snps, $reads, \@alleles);

    $gtype = "";
    foreach $allele (@alleles) {
	print STDERR "    Adding allele '$allele', which maps to '", $map->{$allele}, "' to the genotype\n" if ($debug);
	push(@types, defined($map->{$allele}) ? $map->{$allele} : "-");
    }
    $gtype = join("", sort @types);

    print STDERR "  Ending Genotype: $gtype\n" if ($debug);

    return $gtype;
}

sub call_alleles {
    my ($snps, $reads, $alleles) = @_;

    my ($row, $height, $allele, $base, $snp);

    $height = scalar(@{$reads});

    my %a;

    foreach $row (0..$height-1) {
	$allele = "";

	foreach $snp (@{$snps}) {
	    $base = $reads->[$row][$snp->{'col'}];

	    #
	    # Check to make sure the nucleotide at the location of this SNP is
	    # of one of the two possible states the multinomial model called.
	    #
	    if ($base eq $snp->{'rank_1'} || $base eq $snp->{'rank_2'}) { 
		$allele .= $base;
	    } else {
		last;
	    }
	}

       $a{$allele}++ if (length($allele) == scalar(@{$snps}));
    }

    push(@{$alleles}, keys %a);
}

sub fetch_snps {
    my ($sth, $snps) = @_;

    my ($row);
    #
    # Create a two-dimensional array, each row containing one read.
    #
    $sth->{'snps'}->execute($batch_id)
	or die("Unable to select results from $db.\n");

    while ($row = $sth->{'snps'}->fetchrow_hashref()) {
	if (!defined($snps->{$row->{'tag_id'}})) {
	    $snps->{$row->{'tag_id'}} = [];
	}
	push(@{$snps->{$row->{'tag_id'}}}, $row);
    }
}

sub fetch_reads {
    my ($sth, $sample_id, $tag_id, $reads) = @_;

    my ($row);
    #
    # Create a two-dimensional array, each row containing one read.
    #
    $sth->{'reads'}->execute($sample_id, $tag_id)
	or die("Unable to select results from $db.\n");

    while ($row = $sth->{'reads'}->fetchrow_hashref()) {
	next if ($row->{'relationship'} eq "consensus");

	my $aref;
	@{$aref} = split(//, $row->{'seq'});
	push(@{$reads}, $aref);
    }
}

sub check_homozygosity {
    my ($reads, $col, $rank_1, $rank_2) = @_;

    my (@res, $height, $nuc, $homozygous, $row);

    print STDERR "  Examining col $col, rank 1: $rank_1; rank 2: $rank_2\n" if ($debug);

    $height     = scalar(@{$reads});
    $homozygous = true;

    return unknown if ($height < min_hom_seqs);

    $nuc = {'A' => 0, 'C' => 0, 'G' => 0, 'T' => 0};

    foreach $row (0..$height-1) {
	$nuc->{$reads->[$row][$col]}++;
    }

    @res  = sort {$nuc->{$b} <=> $nuc->{$a}} keys %{$nuc};

    #
    # Check if more than a single nucleotide occurs in this column. Only 
    # count nucleotides that are part of the called SNP, do not count 
    # error-generated nucleotides.
    #
    if (($nuc->{$res[0]} > 0) && 
	($res[0] eq $rank_1 || $res[0] eq $rank_2) &&
	$nuc->{$res[1]} > 0 && 
	($res[1] eq $rank_1 || $res[1] eq $rank_2)) {
	$homozygous = false;
    }

    #
    # If we find a second nucleotide present, check its prevelence. If it is 
    # less than 1/20 of the total reads, don't count a heterozygote. If it is 
    # less than 1/10 report that we can't tell if its homozygous or not. Otherwise,
    # report this tag as a heterozygote.
    #
    my $frac = $nuc->{$res[1]} / $height;
    print STDERR "    Frac: $frac; Second-most Prominent Nuc count: $nuc->{$res[1]}; Depth: $height\n" if ($debug);

    if ($homozygous == false && $frac < min_het_seqs) {
	$homozygous = true;
    } elsif ($homozygous == false && $frac < max_het_seqs) {
	$homozygous = unknown;
    }

    return $homozygous;
}

sub prepare_sql_handles {
    my ($sth, $outg) = @_;

    #
    # Connect to the database, check for the existence of a MySQL config file in the home
    # directory first, otherwise use the stacks-distributed one.
    #
    my $cnf = (-e $ENV{"HOME"} . "/.my.cnf") ? $ENV{"HOME"} . "/.my.cnf" : $mysql_config;
    $sth->{'dbh'} = DBI->connect("DBI:mysql:$db:mysql_read_default_file=$cnf")
	or die("Unable to connect to the $db MySQL Database!\n" . $DBI::errstr);

    my $query;

    $query = 
	"SELECT sample_id, tag_id " . 
	"FROM matches " . 
	"WHERE batch_id=? and catalog_id=?";
    $sth->{'marker_list'} = $sth->{'dbh'}->prepare($query) or die($sth->{'dbh'}->errstr());

    $query = 
	"SELECT samples.sample_id, tag_id, allele, file, type " . 
	"FROM matches " . 
	"JOIN samples ON (matches.sample_id=samples.id) " . 
	"WHERE samples.batch_id=? and catalog_id=?";
    $sth->{'match'} = $sth->{'dbh'}->prepare($query) or die($sth->{'dbh'}->errstr());

    $query = 
	"SELECT tag_id, marker, external_id, chr, bp FROM catalog_index " . 
        "LEFT JOIN catalog_annotations as ca ON (ca.batch_id=catalog_index.batch_id AND ca.catalog_id=catalog_index.tag_id) ".
	"WHERE marker!='' AND progeny>? AND catalog_index.batch_id=?";
    if ($blast_only) {
        $query .= " AND blast_hits>0";
    } elsif ($exclude_blast) {
        $query .= " AND blast_hits=0";
    }
    $sth->{'loci'} = $sth->{'dbh'}->prepare($query) or die($sth->{'dbh'}->errstr());

    $query =
	"SELECT id, file FROM samples WHERE batch_id=? ORDER BY id";
    $sth->{'samp'} = $sth->{'dbh'}->prepare($query) or die($sth->{'dbh'}->errstr());

    $query =
	"SELECT relationship, tag_id, seq FROM unique_tags WHERE sample_id=?";
    $sth->{'reads'} = $sth->{'dbh'}->prepare($query) or die($sth->{'dbh'}->errstr());

    $query =
	"SELECT tag_id, col, rank_1, rank_2 FROM catalog_snps WHERE batch_id=?";
    $sth->{'snps'} = $sth->{'dbh'}->prepare($query) or die($sth->{'dbh'}->errstr());

    $query =
	"SELECT file FROM samples WHERE type='parent' AND batch_id=? ORDER BY id";
    $sth->{'order'} = $sth->{'dbh'}->prepare($query) or die($sth->{'dbh'}->errstr());

    $query =
	"SELECT gc.catalog_id, gc.sample_id, gc.genotype, file " . 
        "FROM genotype_corrections as gc " . 
        "JOIN samples ON (gc.sample_id=samples.id) " .
        "WHERE gc.batch_id=?";
    $sth->{'corr'} = $sth->{'dbh'}->prepare($query) or die($sth->{'dbh'}->errstr());
}

sub close_sql_handles {
    my ($sth) = @_;

    my $key;

    foreach $key (keys %{$sth}) {
	next if ($key =~ /dbh/);

	$sth->{$key}->finish();
    }

    foreach $key (keys %{$sth}) {
	next if ($key !~ /dbh/);

	$sth->{$key}->disconnect();
    }
}

sub parse_command_line {
    while (@ARGV) {
	$_ = shift @ARGV;
	if    ($_ =~ /^-d$/) { $debug++; }
	elsif ($_ =~ /^-D$/) { $db            = shift @ARGV; }
	elsif ($_ =~ /^-P$/) { $in_path       = shift @ARGV; }
        elsif ($_ =~ /^-A$/) { $cache_path    = shift @ARGV; }
	elsif ($_ =~ /^-b$/) { $batch_id      = shift @ARGV; }
	elsif ($_ =~ /^-t$/) { $map_type      = shift @ARGV; }
	elsif ($_ =~ /^-o$/) { $out_type      = shift @ARGV; }
	elsif ($_ =~ /^-c$/) { $corrections++; }
	elsif ($_ =~ /^-C$/) { $man_corrections++; }
	elsif ($_ =~ /^-L$/) { $blast_only++; }
        elsif ($_ =~ /^-E$/) { $exclude_blast++; }
        elsif ($_ =~ /^-i$/) { $expand_id++; }
        elsif ($_ =~ /^-I$/) { $impute++; }
	elsif ($_ =~ /^-s$/) { $sql++; }
	elsif ($_ =~ /^-B$/) { $bl_file       = shift @ARGV; }
	elsif ($_ =~ /^-W$/) { $wl_file       = shift @ARGV; }
	elsif ($_ =~ /^-p$/) { $progeny_limit = shift @ARGV; }
	elsif ($_ =~ /^-v$/) { version(); exit(); }
	elsif ($_ =~ /^-h$/) { usage(); }
	else {
	    print STDERR "Unknown command line options received: $_\n";
	    usage();
	}
    }

    if (length($db) == 0) {
	print STDERR "You must specify a RAD-Tag database.\n";
	usage();
    }

    if ($batch_id == 0) {
	print STDERR "You must specify a batch ID.\n";
	usage();
    }

    if ($corrections && length($in_path) == 0) {
	print STDERR "You must specify a path to the directory containing Stacks output files.\n";
	usage();
    }
    $in_path = substr($in_path, 0, -1) if (substr($in_path, -1) eq "/");

    if ($blast_only > 0 && $exclude_blast > 0) {
        print STDERR "You can only select '-B' or '-E', not both.\n";
        usage();
    }

    if (length($wl_file) > 0 && !-e $wl_file) {
	print STDERR "Unable to locate the whitelist '$wl_file'.\n";
	usage();
    }

    if (length($bl_file) > 0 && !-e $bl_file) {
	print STDERR "Unable to locate the blacklist '$bl_file'.\n";
	usage();
    }

    if (!defined($map_types->{$map_type})) {
        print STDERR "You must specify a valid map type. 'CP', 'DH', 'F2' and 'BC1' are the currently supported map types.\n";
        usage();
    }

    if (!defined($out_types->{$out_type})) {
        print STDERR "You must specify a valid output file type. 'joinmap', and 'rqtl' are the currently supported output types.\n";
        usage();
    }
}

sub version {
    print STDERR "genotypes.pl ", stacks_version, "\n";

}

sub usage {
    version();

    print << "EOQ";
genotypes.pl -t map_type -D db -b id [-p min] [-c -P path -A path] [-o type] [-s] [-C] [-B blacklist] [-W whitelist] [-d] [-h]
  t: map type to write. 'CP', 'DH', 'F2' and 'BC1' are the currently supported map types.
  o: output file type to write, 'joinmap' and 'rqtl' are currently supported.
  D: Stacks database to examine.
  b: Batch ID to examine when exporting from the catalog.
  p: minimum number of progeny required to print a marker.
  c: make automated corrections to the data.
  P: path to the Stacks output files.
  A: specify a path to the data cache.
  C: apply manual corrections (that were made via the web interface) to the data.
  s: output a file to import results into an SQL database.
  B: specify a file containing Blacklisted markers to be excluded from the export.
  W: specify a file containign Whitelisted markers to include in the export.
  h: display this help message.
  d: turn on debug output.

EOQ
    exit(0);
}
