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
use CHI;

use constant true         => 1;
use constant false        => 0;
use constant unknown      => -1;
use constant min_hom_seqs => 5;
use constant min_het_seqs => 0.05;
use constant max_het_seqs => 0.1;

use constant stacks_version => "_VERSION_";

my $mysql_config    = "_PKGDATADIR_" . "sql/mysql.cnf";
my $debug           = 0;
my $db              = "";
my $map_type        = "";
my $cache_path      = "_LOCALSTATEDIR_" . "caches";
my $uni_path        = "";
my $progeny_limit   = 4;
my $batch_id        = 0;
my $category        = 0;
my $man_corrections = 0;
my $corrections     = 0;
my $blast_only      = 0;
my $exclude_blast   = 0;
my $expand_id       = 0;
my $out_file        = "";
my $bl_file         = "";
my $wl_file         = "";
my $sql             = 0;

parse_command_line();

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

if ($map_type eq "CP") {
    export_joinmap_cp_map(\%sth, $cache_path);

} elsif ($map_type eq "BC1") {
    export_joinmap_bc1_map(\%sth, $cache_path);

} elsif ($map_type eq "DH") {
    export_joinmap_dh_map(\%sth, $cache_path);

} else {
    print STDERR "Unknown map type.\n";
}

close_sql_handles(\%sth);

sub export_joinmap_dh_map {
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

    my (%types, %genotypes, %order, @samples, %sample_ids, %marker_list, @loci, $locus, $mtype, $key, $file);

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
    }

    #
    # Fetch the parents/progeny for each locus
    #
    my $i = 1;
    my $num_loci = scalar(@loci);
    foreach $locus (@loci) {
	print STDERR "Processing locus $i of $num_loci          \r";

	$locus->{'progeny'} = {};

	call_dh_genotypes($sth, $cache, \%sample_ids, \%genotypes, \%order, 
                          $locus->{'id'}, $types{$locus->{'marker'}}, $locus->{'progeny'});

	$i++;
    }
    print STDERR "\n";

    apply_corrected_genotypes($sth, \@loci) if ($man_corrections);

    if ($category) {
	foreach $key (keys %types) {
	    $mtype = $types{$key};

	    # Alter the filename to reflect the markers it contains
	    $file = $out_file;
	    $file =~ s/\.loc$//;
	    $file .= "." . $mtype;
	    $file .= ".loc";
	    print STDERR "Filename: $file\n";

	    # Escape characters that have special meaning in the regex.
	    $key =~ s/\-/\\\-/g;

	    write_loci($file, \%types, \@samples, \@loci, $key);
	}

    } else {
	write_loci($out_file, \%types, \@samples, \@loci, ".*");
    }

    write_loci_sql($out_file, \%types, \@samples, \%sample_ids, \@loci) if ($sql);
}

sub export_joinmap_bc1_map {
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

    my (%types, %genotypes, %order, @samples, %sample_ids, %marker_list, @loci, $locus, $mtype, $key, $file);

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
    }

    #
    # Fetch the parents/progeny for each locus
    #
    my $i = 1;
    my $num_loci = scalar(@loci);
    foreach $locus (@loci) {
	print STDERR "Processing locus $i of $num_loci          \r";

	$locus->{'progeny'} = {};

	call_bc1_genotypes($sth, $cache, \%sample_ids, \%genotypes, \%order, 
                           $locus->{'id'}, $types{$locus->{'marker'}}, $locus->{'progeny'});

	$i++;
    }
    print STDERR "\n";

    apply_corrected_genotypes($sth, \@loci) if ($man_corrections);

    if ($category) {
	foreach $key (keys %types) {
	    $mtype = $types{$key};

	    # Alter the filename to reflect the markers it contains
	    $file = $out_file;
	    $file =~ s/\.loc$//;
	    $file .= "." . $mtype;
	    $file .= ".loc";
	    print STDERR "Filename: $file\n";

	    # Escape characters that have special meaning in the regex.
	    $key =~ s/\-/\\\-/g;

	    write_loci($file, \%types, \@samples, \@loci, $key);
	}

    } else {
	write_loci($out_file, \%types, \@samples, \@loci, ".*");
    }

    write_loci_sql($out_file, \%types, \@samples, \%sample_ids, \@loci) if ($sql);
}

sub export_joinmap_cp_map {
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

    my (%types, %genotypes, %order, @samples, %sample_ids, %marker_list, @loci, $locus, $mtype, $key, $file);

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
    }

    #
    # Fetch the parents/progeny for each locus
    #
    my $i = 1;
    my $num_loci = scalar(@loci);
    foreach $locus (@loci) {
	print STDERR "Processing locus $i of $num_loci          \r";

	$locus->{'progeny'} = {};

	call_cp_genotypes($sth, $cache, \%sample_ids, \%genotypes, \%order, 
                          $locus->{'id'}, $types{$locus->{'marker'}}, $locus->{'progeny'});

	$i++;
    }
    print STDERR "\n";

    apply_corrected_genotypes($sth, \@loci) if ($man_corrections);

    if ($category) {
	foreach $key (keys %types) {
	    $mtype = $types{$key};

	    # Alter the filename to reflect the markers it contains
	    $file = $out_file;
	    $file =~ s/\.loc$//;
	    $file .= "." . $mtype;
	    $file .= ".loc";
	    print STDERR "Filename: $file\n";

	    # Escape characters that have special meaning in the regex.
	    $key =~ s/\-/\\\-/g;

	    write_loci($file, \%types, \@samples, \@loci, $key);
	}

    } else {
	write_loci($out_file, \%types, \@samples, \@loci, ".*");
    }

    write_loci_sql($out_file, \%types, \@samples, \%sample_ids, \@loci) if ($sql);
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

	$path = $uni_path . "/" . $sample . ".tags.tsv";

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
    my ($file, $types, $samples, $sample_ids, $loci) = @_;

    my ($row, $out_fh, $locus, $num_loci, $sample);

    $file =~ s/\.loc$/\.txt/;

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

sub write_loci {
    my ($file, $types, $samples, $loci, $marker_regex) = @_;

    my ($row, $out_fh, $bl_fh, $locus, $num_loci, $sample, $key, $progeny_cnt, $id);

    my %blacklist;
    #
    # Read in a blacklist
    #
    if (length($bl_file) > 0) {
        open($bl_fh, "<$bl_file") or die("Unable to open blacklist file '$file'; $!\n");

        while (<$bl_fh>) {
            chomp $_;
            $blacklist{$_}++;
        }

        print STDERR "Read ", scalar(keys %blacklist), " blacklisted markers.\n";

        close($bl_fh);
    }

    my ($pop_name) = ($out_file =~ /\/(.+)\..+$/);

    open($out_fh, ">$file") or die("Unable to open output file '$file'; $!\n");

    #
    # Count the number of mappable progeny
    #
    $num_loci = 0;

    foreach $locus (@{$loci}) {
	next if (defined($blacklist{$locus->{'id'}}));
	next if (scalar(keys %{$locus->{'progeny'}}) < $progeny_limit);
	next if ($locus->{'marker'} !~ /$marker_regex/);

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
	next if ($locus->{'marker'} !~ /$marker_regex/);

        #
        # Count the number of progeny
        #
        $progeny_cnt = 0;
	foreach $key (keys %{$locus->{'progeny'}}) {
            $progeny_cnt++ if ($key ne "--" && $key ne "-");
        }
        next if ($progeny_cnt < $progeny_limit);

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
        #next if ($row->{'tag_id'} != 4);

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

sub call_dh_genotypes {
    my ($sth, $cache, $sample_ids, $genotypes, $order, $tag_id, $marker, $markers) = @_;

    my (@keys, $key, $row, $allele, $m, $ptags, %parents, %progeny, %genotype_map, %rev_geno_map);

    my $check_genotypes = {'abx--' => \&check_lmx___genotype,
                           '--xab' => \&check___xnp_genotype
			   };

    my $create_genotype_map = {'abx--' => \&create_simple_genotype_map,
                               '--xab' => \&create_simple_genotype_map
                               };

    return if (!defined($create_genotype_map->{$marker}));

    fetch_progeny($sth, $order, $tag_id, \%parents, \%progeny);

    #
    # Create a map between alleles and genotype symbols
    #
    $create_genotype_map->{$marker}->($order, $marker, $tag_id, \%parents, \%genotype_map);

    #
    # Create a reverse map
    #
    foreach $key (keys %genotype_map) {
	print STDERR "Key: $key, Map: $genotype_map{$key}\n" if ($debug);
	$rev_geno_map{$genotype_map{$key}} = $key;
    }

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
	@alleles = sort @alleles;

        $m = $alleles[0];
        $m = "-" if (grep(/^\-$/, @alleles));

	#
	# Check the genotype, correct any errors.
	#
        if ($corrections > 0) {
            $m = $check_genotypes->{$marker}->($sth, $cache, $tag_id, $sample_ids->{$key}, $keys[0], \%genotype_map, \%rev_geno_map, $m);
        }

	if (!defined($genotypes->{$marker}->{lc($m)})) {
	    print STDERR "Warning: illegal genotype encountered in tag $tag_id, progeny $key, '$m'\n" if ($debug);
	} else  {
	    $markers->{$key} = $m;
	    print STDERR "  $key allele: $markers->{$key}\n" if ($debug);
	}
    }
}

sub call_bc1_genotypes {
    my ($sth, $cache, $sample_ids, $genotypes, $order, $tag_id, $marker, $markers) = @_;

    my (@keys, $key, $row, $allele, $m, $ptags, %parents, %progeny, %genotype_map, %rev_geno_map);

    my $check_genotypes = {'aaxbb' => \&check_aaxbb_genotype,
                           'bbxaa' => \&check_aaxbb_genotype,
			   'abxcc' => \&check_abxcc_genotype,
                           'ccxab' => \&check_ccxab_genotype
			   };

    my $create_genotype_map = {'aaxbb' => \&create_aaxbb_genotype_map,
                               'bbxaa' => \&create_bbxaa_genotype_map,
                               'abxcc' => \&create_abxcc_genotype_map,
                               'ccxab' => \&create_ccxab_genotype_map
                               };

    return if (!defined($create_genotype_map->{$marker}));

    fetch_progeny($sth, $order, $tag_id, \%parents, \%progeny);

    #
    # Create a map between alleles and genotype symbols
    #
    $create_genotype_map->{$marker}->($order, $marker, $tag_id, \%parents, \%genotype_map);

    #
    # Create a reverse map
    #
    foreach $key (keys %genotype_map) {
	print STDERR "Key: $key, Map: $genotype_map{$key}\n" if ($debug);
	$rev_geno_map{$genotype_map{$key}} = $key;
    }

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
	@alleles = sort @alleles;

        if (grep(/^h$/, @alleles)) {
            $m = "h";
        } else {
            $m = $alleles[0];
        }
        $m = "-" if (grep(/^\-$/, @alleles));

	#
	# Check the genotype, correct any errors.
	#
        if ($corrections > 0) {
            $m = $check_genotypes->{$marker}->($sth, $cache, $tag_id, $sample_ids->{$key}, $keys[0], \%genotype_map, \%rev_geno_map, $m);
        }

	if (!defined($genotypes->{$marker}->{lc($m)})) {
	    print STDERR "Warning: illegal genotype encountered in tag $tag_id, progeny $key, '$m'\n" if ($debug);
	} else  {
	    $markers->{$key} = $m;
	    print STDERR "  $key allele: $markers->{$key}\n" if ($debug);
	}
    }
}

sub call_cp_genotypes {
    my ($sth, $cache, $sample_ids, $genotypes, $order, $tag_id, $marker, $markers) = @_;

    my (@keys, $key, $row, $allele, $m, $ptags, %parents, %progeny, %genotype_map, %rev_geno_map);

    my $check_genotypes = {'lmx--' => \&check_lmx___genotype,
			   '--xnp' => \&check___xnp_genotype,
			   'lmxll' => \&check_lmxll_genotype,
			   'nnxnp' => \&check_nnxnp_genotype,
			   'hkxhk' => \&check_hkxhk_genotype,
			   'efxeg' => \&check_efxeg_genotype,
			   'abxcd' => \&check_abxcd_genotype
			   };

    my $create_genotype_map = {'lmx--' => \&create_simple_genotype_map,
                               '--xnp' => \&create_simple_genotype_map,
                               'lmxll' => \&create_simple_genotype_map,
                               'nnxnp' => \&create_simple_genotype_map,
                               'hkxhk' => \&create_simple_genotype_map,
                               'efxeg' => \&create_efxeg_genotype_map,
                               'abxcd' => \&create_abxcd_genotype_map
                               };

    fetch_progeny($sth, $order, $tag_id, \%parents, \%progeny);

    #
    # Create a map between alleles and genotype symbols
    #

    $create_genotype_map->{$marker}->($order, $marker, $tag_id, \%parents, \%genotype_map);

    #
    # Create a reverse map
    #
    foreach $key (keys %genotype_map) {
	print STDERR "Key: $key, Map: $genotype_map{$key}\n" if ($debug);
	$rev_geno_map{$genotype_map{$key}} = $key;
    }

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
	@alleles = sort @alleles;

	#
	# If the tag was male or female-only in the parents, than we can no have more 
	# than a single allele in any individual progeny.
	# 
	if (scalar(@alleles) == 2 &&
	    ($marker eq "lmx--" || $marker eq "--xnp")) {
	    # Illegal genotype
	    splice(@alleles, 0);

	#
	# In the cases where we have a common genotype between the parent, such as 'aa' or 'ee'
	# the pipeline will only report a single 'a' or 'e', since all the RAD-Tags for the two
	# were merged into a single allele by the pipeline.
	#
	} elsif (scalar(@alleles) == 1 && 
	    ($alleles[0] eq 'l' || $alleles[0] eq 'n' || 
	     $alleles[0] eq 'h' || $alleles[0] eq 'k')) {
	    push(@alleles, $alleles[0]);

	#
	# We are mapping ab/-- and --/ab tags as lmxll and nnxnp, respectively, since 
	# Joinmap does not have a specification for heterozygous, single parent tags.
	#
	} elsif (scalar(@alleles) == 1 && 
		 ($marker eq "lmx--" || $marker eq "--xnp")) {
	    unshift(@alleles, 'l') if ($alleles[0] eq 'm');
	    unshift(@alleles, 'n') if ($alleles[0] eq 'p');
	}

	$m = join("", @alleles);

	#
	# Check the genotype, correct any errors.
	#
        if ($corrections > 0) {
            $m = $check_genotypes->{$marker}->($sth, $cache, $tag_id, $sample_ids->{$key}, $keys[0], \%genotype_map, \%rev_geno_map, $m);
        }

	if (!defined($genotypes->{$marker}->{lc($m)})) {
	    print STDERR "Warning: illegal genotype encountered in tag $tag_id, progeny $key, '$m'\n" if ($debug);
	} else  {
	    $markers->{$key} = $m;
	    print STDERR "  $key allele: $markers->{$key}\n" if ($debug);
	}
    }
}

sub create_efxeg_genotype_map {
    my ($order, $marker, $tag_id, $parents, $map) = @_;

    my (@parents, @keys, %alleles, %rev_map);
    my ($allele, $key);

    #
    # First step, identify the allele that is common between the two parents.
    #
    @parents = keys %{$parents};
    foreach $allele (@{$parents->{$parents[0]}}, @{$parents->{$parents[1]}}) {
	$alleles{$allele}++;
    }
    @keys = sort {$alleles{$b} <=> $alleles{$a}} keys %alleles;

    if ($alleles{$keys[0]} != 2) {
        print STDERR "Warning: unable to determine genotype map for tag $tag_id with marker '$marker'\n";
        return;
    }

    $map->{$keys[0]} = 'e';
    $rev_map{'e'}    = $keys[0];

    print STDERR "Assigning '$keys[0]' to genotype 'e'\n" if ($debug);

    #
    # Now assign the remaining alleles, 'f' to the second allele in the first parent,
    # 'g' to the second allele in the second parent
    #
    $key = $order->{$parents[0]} eq "first" ? $parents[0] : $parents[1];

    $allele = $parents->{$key}->[0] eq $rev_map{'e'} ? 
        $parents->{$key}->[1] : $parents->{$key}->[0];
    $map->{$allele} = 'f';
    $rev_map{'f'} = $allele;
    print STDERR "Assigning '$rev_map{'f'}' to genotype 'f'\n" if ($debug);

    $key = $order->{$parents[0]} eq "second" ? $parents[0] : $parents[1];

    $allele = $parents->{$key}->[0] eq $rev_map{'e'} ? 
        $parents->{$key}->[1] : $parents->{$key}->[0];
    $map->{$allele} = 'g';
    $rev_map{'g'} = $allele;
    print STDERR "Assigning '$rev_map{'g'}' to genotype 'g'\n" if ($debug);
}

sub create_abxcd_genotype_map {
    my ($order, $marker, $tag_id, $parents, $map) = @_;

}

sub create_aaxbb_genotype_map {
    my ($order, $marker, $tag_id, $parents, $map) = @_;

    my (%genotypes, @parents, @types, %alleles);
    my ($key);

    @parents = keys %{$parents};

    $key = $order->{$parents[0]} eq "second" ? $parents[0] : $parents[1];
    $map->{$parents->{$key}->[0]} = "h";

    $key = $order->{$parents[0]} eq "first"  ? $parents[0] : $parents[1];
    $map->{$parents->{$key}->[0]} = "b";
}

sub create_bbxaa_genotype_map {
    my ($order, $marker, $tag_id, $parents, $map) = @_;

    my (%genotypes, @parents, @types, %alleles);
    my ($key);

    @parents = keys %{$parents};

    $key = $order->{$parents[0]} eq "first"  ? $parents[0] : $parents[1];
    $map->{$parents[$key]->[0]} = "h";

    $key = $order->{$parents[0]} eq "second" ? $parents[0] : $parents[1];
    $map->{$parents[$key]->[0]} = "a";
}

sub create_abxcc_genotype_map {
    my ($order, $marker, $tag_id, $parents, $map) = @_;

    my (%genotypes, @parents, @types, %alleles);
    my ($key, $allele);

    @parents = keys %{$parents};

    $key = $order->{$parents[0]} eq "second" ? $parents[0] : $parents[1];
    $map->{$parents->{$key}->[0]} = "h";

    $key = $order->{$parents[0]} eq "first"  ? $parents[0] : $parents[1];

    foreach $allele (@{$parents->{$parents[$key]}}) {
        $map->{$allele} = "b";
    }
}

sub create_ccxab_genotype_map {
    my ($order, $marker, $tag_id, $parents, $map) = @_;

    my (%genotypes, @parents, @types, %alleles);
    my ($key, $allele);

    @parents = keys %{$parents};

    $key = $order->{$parents[0]} eq "second" ? $parents[0] : $parents[1];
    foreach $allele (@{$parents->{$parents[$key]}}) {
        $map->{$allele} = "h";
    }

    $key = $order->{$parents[0]} eq "first"  ? $parents[0] : $parents[1];
    $map->{$parents->{$key}->[0]} = "b";
}

sub create_simple_genotype_map {
    my ($order, $marker, $tag_id, $parents, $map) = @_;

    my (%genotypes, @parents, @types, %alleles);
    my ($type, $m, $key, $allele);
    #
    # Create a genotype map. For any set of alleles, this routine will
    # assign each allele to one of the constituent genotypes, e.g. given the 
    # marker type 'abxaa' and the alleles 'A' and 'G' from the male, and 'G'
    # from the female, will assign 'G' == 'a' and 'A'== 'b'.
    #
    $m = substr($marker, 0, 2) . substr($marker, 3, 2);

    foreach $type (split(//, $m)) {
	next if ($type eq "-");
	$genotypes{$type}++;
    }

    #
    # Sort genotypes alphabetically.
    #
    @types = sort keys %genotypes;

    @parents = keys %{$parents};
    foreach $allele (@{$parents->{$parents[0]}}, @{$parents->{$parents[1]}}) {
	$alleles{$allele}++;
    }

    #
    # Order the alleles by the number of times each occurs. Given the lmxll or nnxnp tags,
    # the more numerous allele will always be assigned to the genotype that occurs
    # sooner in the alphabet ('l' before 'm', 'n' before 'p')
    #
    foreach $allele (sort {$alleles{$b} <=> $alleles{$a}} keys %alleles) {

	if (!defined($map->{$allele})) {
	    $key = shift @types;

            print STDERR "Assinging '$allele' to genotype '$key'\n" if ($debug);

	    if (!defined($key)) {
		print STDERR "Warning: impossible allele combination in parental genotypes, tag $tag_id\n";
		next;
	    }
	    $map->{$allele} = $key;
	}
    }
}

sub check_lmx___genotype {
    my ($sth, $cache, $catalog_id, $sample_id, $tag_id, $map, $rev_map, $gtype) = @_;

    return $gtype;
}

sub check___xnp_genotype {
    my ($sth, $cache, $catalog_id, $sample_id, $tag_id, $map, $rev_map, $gtype) = @_;

    return $gtype;
}

sub check_lmxll_genotype {
    my ($sth, $cache, $catalog_id, $sample_id, $tag_id, $map, $rev_map, $gtype) = @_;

    print STDERR "Sample: $sample_id; Tag: $tag_id: Genotype; $gtype\n" if ($debug);

    my ($reads, $row, $key, $homozygous);

    $key    = $sample_id . "_" . $tag_id;
    $reads  = $cache->get($key);

    if (!defined($reads)) {
	print STDERR "Unable to locate cached object with key $key\n"; 
	return;
    }

    $sth->{'snp'}->execute($batch_id, $catalog_id)
	or die("Unable to select results from $db.\n");
    $row = $sth->{'snp'}->fetchrow_hashref();

    $homozygous = check_homozygosity($reads, $row->{'col'}, $row->{'rank_1'}, $row->{'rank_2'});

    if ($gtype eq "ll") {
	if ($homozygous == true) {
	    $gtype = "ll"; 
	} elsif ($homozygous == false) {
	    $gtype = uc("lm");
	} else {
	    $gtype = "--";
	}
    } elsif ($gtype eq "m") {
	if ($homozygous) {
	    $gtype = "m"; 
	} elsif (!$homozygous) {
	    $gtype = uc("lm");
	} else {
	    $gtype = "--";
	}
    }

    print STDERR "  Ending Genotype: $gtype\n" if ($debug);

    return $gtype;
}

sub check_nnxnp_genotype {
    my ($sth, $cache, $catalog_id, $sample_id, $tag_id, $map, $rev_map, $gtype) = @_;

    print STDERR "Sample: $sample_id; Tag: $tag_id: Genotype; $gtype\n" if ($debug);

    my ($reads, $row, $key, $homozygous);

    $key   = $sample_id . "_" . $tag_id;
    $reads = $cache->get($key);

    if (!defined($reads)) {
	print STDERR "Unable to locate cached object with key $key\n"; 
	return;
    }

    $sth->{'snp'}->execute($batch_id, $catalog_id)
	or die("Unable to select results from $db.\n");
    $row = $sth->{'snp'}->fetchrow_hashref();

    $homozygous = check_homozygosity($reads, $row->{'col'}, $row->{'rank_1'}, $row->{'rank_2'});

    if ($gtype eq "nn") {
	if ($homozygous == true) {
	    $gtype = "nn"; 
	} elsif ($homozygous == false) {
	    $gtype = uc("np");
	} else {
	    $gtype = "--";
	}
    } elsif ($gtype eq "p") {
	if ($homozygous == true) {
	    $gtype = "p"; 
	} elsif ($homozygous == false) {
	    $gtype = uc("np");
	} else {
	    $gtype = "--";
	}
    }

    print STDERR "  Ending Genotype: $gtype\n" if ($debug);

    return $gtype;
}

sub check_hkxhk_genotype {
    my ($sth, $cache, $catalog_id, $sample_id, $tag_id, $map, $rev_map, $gtype) = @_;

    print STDERR "Sample: $sample_id; Tag: $tag_id: Genotype; $gtype\n" if ($debug);

    my ($reads, $row, $key, $homozygous);

    $key   = $sample_id . "_" . $tag_id;
    $reads = $cache->get($key);

    if (!defined($reads)) {
	print STDERR "Unable to locate cached object with key $key\n"; 
	return;
    }

    $sth->{'snp'}->execute($batch_id, $catalog_id)
	or die("Unable to select results from $db.\n");
    $row = $sth->{'snp'}->fetchrow_hashref();

    $homozygous = check_homozygosity($reads, $row->{'col'}, $row->{'rank_1'}, $row->{'rank_2'});

    if ($gtype eq "hh") {
	if ($homozygous == true) {
	    $gtype = "hh"; 
	} elsif ($homozygous == false) {
	    $gtype = uc("hk");
	} else {
	    $gtype = "--";
	}
    } elsif ($gtype eq "kk") {
	if ($homozygous == true) {
	    $gtype = "kk"; 
	} elsif ($homozygous == false) {
	    $gtype = uc("hk");
	} else {
	    $gtype = "--";
	}
    }

    print STDERR "  Ending Genotype: $gtype\n" if ($debug);

    return $gtype;
}

sub check_efxeg_genotype {
    my ($sth, $cache, $catalog_id, $sample_id, $tag_id, $map, $rev_map, $gtype) = @_;

    print STDERR "Sample: $sample_id; Tag: $tag_id: Genotype; $gtype\n" if ($debug);

    my ($reads, $row, @snps, @alleles, @nucs, $nuc, $allele, $rev_gtype, $key, $homozygous, $a, $orig_gtype);

    $orig_gtype = $gtype;
    $key        = $sample_id . "_" . $tag_id;
    $reads      = $cache->get($key);

    if (!defined($reads)) {
	print STDERR "Unable to locate cached object with key $key\n"; 
	return;
    }

    @alleles = split(//, $gtype);

    #
    # The ee genotype is the only genotype we can correct for an efxeg marker.
    #
    return $gtype if (scalar(@alleles) > 1);

    $sth->{'snp'}->execute($batch_id, $catalog_id)
	or die("Unable to select results from $db.\n");

    while ($row = $sth->{'snp'}->fetchrow_hashref()) {
        push(@snps, $row);
    }

    #
    # Iterate over each allele present in this progeny RAD-Tag. Check the column for each SNP in
    # the catalog tag this progeny tag matched. Check those SNP to ensure there is sufficient
    # depth to call a homozygote, or if there is sufficient reads to call a heterozygote.
    #
    foreach $a (@alleles) {

        $rev_gtype = $rev_map->{$a};
        print STDERR "  Allele: ", $a, "; reverse genotype: ", $rev_gtype, "\n" if ($debug);

        @nucs   = split(//, $rev_gtype);
        $allele = "";

        foreach $row (@snps) {
            $homozygous = check_homozygosity($reads, $row->{'col'}, $row->{'rank_1'}, $row->{'rank_2'});
            $nuc        = shift @nucs;

            if ($homozygous == unknown) {
                $allele .= "N";
            } else {
                $allele .= $nuc;
            }
        }

        print STDERR "    Adding allele '$allele', which maps to '", $map->{$allele}, "' to the genotype\n" if ($debug);
        $gtype .= defined($map->{$allele}) ? $map->{$allele} : "N";
    }

    $gtype = join("", sort split(//, $gtype));
    print STDERR 
	"  Allele: $allele; Gtype: $gtype\n",
	"  Ending Genotype: $gtype\n" if ($debug);

    if ($gtype eq $orig_gtype) {
        return $gtype;
    } else {
        return uc($gtype);
    }
}

sub check_abxcd_genotype {
    my ($sth, $cache, $catalog_id, $sample_id, $tag_id, $map, $rev_map, $gtype) = @_;

    return $gtype;
}

sub check_aaxbb_genotype {
    my ($sth, $cache, $catalog_id, $sample_id, $tag_id, $map, $rev_map, $gtype) = @_;

    return $gtype if ($gtype eq "h");

    print STDERR "Sample: $sample_id; Tag: $tag_id: Genotype; $gtype\n" if ($debug);

    my ($reads, $row, $key, $homozygous, $rev_gtype, @snps, @nucs, $nuc, $allele, $orig_gtype);

    $orig_gtype = $gtype;
    $key        = $sample_id . "_" . $tag_id;
    $reads      = $cache->get($key);

    if (!defined($reads)) {
	print STDERR "Unable to locate cached object with key $key\n"; 
	return;
    }

    $sth->{'snp'}->execute($batch_id, $catalog_id)
	or die("Unable to select results from $db.\n");

    while ($row = $sth->{'snp'}->fetchrow_hashref()) {
        push(@snps, $row);
    }

    $rev_gtype = $rev_map->{$gtype};
    print STDERR "  Allele: ", $gtype, "; reverse genotype: ", $rev_gtype, "\n" if ($debug);

    @nucs   = split(//, $rev_gtype);
    $allele = "";

    foreach $row (@snps) {
        $homozygous = check_homozygosity($reads, $row->{'col'}, $row->{'rank_1'}, $row->{'rank_2'});
        $nuc        = shift @nucs;

        if ($homozygous == true) {
            $allele .= $nuc;
        } elsif ($homozygous == false) {
            $allele .= ($nuc eq $row->{'rank_1'}) ? $row->{'rank_2'} : $row->{'rank_1'};
        } else {
            $allele .= "N";
        }
    }

    print STDERR "    Adding allele '$allele', which maps to '", $map->{$allele}, "' to the genotype\n" if ($debug);
    $gtype = defined($map->{$allele}) ? $map->{$allele} : $gtype;

    print STDERR "  Ending Genotype: $gtype\n" if ($debug);

    if ($gtype eq $orig_gtype) {
        return $gtype;
    } else {
        return uc($gtype);
    }
}

sub check_abxcc_genotype {
    my ($sth, $cache, $catalog_id, $sample_id, $tag_id, $map, $rev_map, $gtype) = @_;

    return $gtype;
}

sub check_ccxab_genotype {
    my ($sth, $cache, $catalog_id, $sample_id, $tag_id, $map, $rev_map, $gtype) = @_;

    return $gtype;
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
	"SELECT tag_id, marker, external_id FROM catalog_index " . 
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
	"SELECT col, rank_1, rank_2 FROM catalog_snps WHERE batch_id=? AND tag_id=? ORDER BY col";
    $sth->{'snp'} = $sth->{'dbh'}->prepare($query) or die($sth->{'dbh'}->errstr());

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
	elsif ($_ =~ /^-P$/) { $uni_path      = shift @ARGV; }
        elsif ($_ =~ /^-A$/) { $cache_path    = shift @ARGV; }
	elsif ($_ =~ /^-b$/) { $batch_id      = shift @ARGV; }
	elsif ($_ =~ /^-t$/) { $map_type      = shift @ARGV; }
	elsif ($_ =~ /^-c$/) { $corrections++; }
	elsif ($_ =~ /^-C$/) { $man_corrections++; }
	elsif ($_ =~ /^-B$/) { $blast_only++; }
        elsif ($_ =~ /^-E$/) { $exclude_blast++; }
        elsif ($_ =~ /^-i$/) { $expand_id++; }
	elsif ($_ =~ /^-s$/) { $sql++; }
	elsif ($_ =~ /^-l$/) { $bl_file       = shift @ARGV; }
	elsif ($_ =~ /^-w$/) { $wl_file       = shift @ARGV; }
	elsif ($_ =~ /^-o$/) { $out_file      = shift @ARGV; }
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

    if (length($uni_path) == 0) {
	print STDERR "You must specify a path to the directory containing Stacks output files.\n";
	usage();
    }

    if (length($out_file) == 0) {
	$out_file = "./batch_" . $batch_id . "_catalog.loc";
    }

    if ($blast_only > 0 && $exclude_blast > 0) {
        print STDERR "You can only select '-B' or '-E', not both.\n";
        usage();
    }

    if (length($map_type) == 0 || ($map_type ne "CP" && $map_type ne "BC1" && $map_type ne "DH")) {
        print STDERR "You must specify a map type. 'CP', 'DH' and 'BC1' are the currently supported map types.\n";
        usage();
    }
}

sub version {
    print STDERR "genotypes.pl ", stacks_version, "\n";

}

sub usage {
    version();

    print << "EOQ";
genotypes.pl -t map_type -D db -P path [-b id] [-p min] [-c] [-o path] [-s] [-C] [-B] [-E] [-A path] [-d] [-h]
  t: map type to write. 'CP', 'DH' and 'BC1' are the currently supported map types.
  P: path to the Stacks output files.
  D: radtag database to examine.
  b: Batch ID to examine when exporting from the catalog.
  o: output file.
  c: make automated corrections to the data.
  C: apply manual corrections (that were made via the web interface) to the data.
  p: minimum number of progeny required to print a marker.
  s: output a file to import results into an SQL database.
  E: exclude genotypes related to a marker with a BLAST hit.
  B: only export genotypes related to a marker with BLAST hits.
  A: specify a path to the data cache.
  h: display this help message.
  d: turn on debug output.

EOQ
    exit(0);
}
