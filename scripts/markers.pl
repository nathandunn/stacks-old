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
# Written by Julian Catchen <jatchen@uoregon.edu>
#

use strict;
use DBI;

use constant stacks_version => "_VERSION_";
use constant true  => 1;
use constant false => 0;

my $mysql_config   = "_PKGDATADIR_" . "/sql/mysql.cnf";
my $debug          = 0;
my $db             = "";
my $batch_id       = 0;
my $impute_markers = false;
my $out_path       = "./";

my %genotypes = (
		 'ab/cc' => {'ac' => 0, 'bc' => 0},
		 'cc/ab' => {'ac' => 0, 'bc' => 0},
		 'aa/bb' => {'aa' => 0, 'ab' => 0, 'bb' => 0},
		 'ab/--' => {'aa' => 0, 'bb' => 0, 'ab' => 0},
		 '--/ab' => {'aa' => 0, 'bb' => 0, 'ab' => 0},
		 'ab/aa' => {'aa' => 0, 'ab' => 0, 'bb' => 0},
		 'aa/ab' => {'aa' => 0, 'ab' => 0, 'bb' => 0},
		 'ab/ab' => {'aa' => 0, 'ab' => 0, 'bb' => 0},
		 'ab/ac' => {'aa' => 0, 'ab' => 0, 'ac' => 0, 'bb' => 0, 'bc' => 0},
		 'ab/cd' => {'aa' => 0, 'ab' => 0, 'ac' => 0, 'ad' => 0, 
			     'bb' => 0, 'bc' => 0, 'bd' => 0, 'cc' => 0, 'cd' => 0, 'dd' => 0}
		 );

parse_command_line();

#
# Connect to the database and prepare our queries.
#
my (%sth, %catalog, %order);

prepare_sql_handles(\%sth);

populate(\%sth, \%catalog, \%order);

find_markers(\%sth, \%catalog, \%order);

close_sql_handles(\%sth);

sub impute_marker {
    my ($parent, $progeny) = @_;

    my ($key, @keys, $cnt, $type, $err);

    my %mtypes = (
	'het' => {
	    'ab/cc' => {'ac' => 0, 'bc' => 0}, # 2 hets
	    'ab/--' => {'aa' => 0, 'bb' => 0}, # 2 hom
	    'ab/aa' => {'aa' => 0, 'ab' => 0}, # 1 hom, 1 het
	    'ab/ab' => {'aa' => 0, 'ab' => 0, 'bb' => 0}, # 2 hom, 1 het
	    'ab/ac' => {'aa' => 0, 'ac' => 0, 'bc' => 0, }, # 1 hom, 2 het
	    'ab/cd' => {'ac' => 0, 'ad' => 0, 'bc' => 0, 'bd' => 0} # 4 het
	},
	'hom' => {
	    'cc/ab' => {'ac' => 0, 'bc' => 0}, # 2 hets
	    'aa/bb' => {'ab' => 0}, # 1 het
	    'aa/ab' => {'aa' => 0, 'ab' => 0, 'bb' => 0} # 2 hom, 1 het
	});

    my %gtypes;
    my $het_cnt = 0; 
    my $hom_cnt = 0;
    my $par_cnt = 0;
    my $marker  = "unknown";

    foreach $key (keys %{$progeny}) {
	my $alleles;

	print STDERR "Examining progeny $key\n" if ($debug);
	#
	# Discard progeny with more than one locus matched to this catalog tag.
	#
	@keys = keys %{$progeny->{$key}};
	next if (scalar(@keys) > 1);

	$cnt     = scalar(@{$progeny->{$key}->{$keys[0]}});
	$alleles = join("|", sort @{$progeny->{$key}->{$keys[0]}});
	if (!defined($gtypes{$alleles})) {
	    if ($cnt == 1) {
		$hom_cnt++;
	    } elsif ($cnt == 2) {
		$het_cnt++;
	    }
	}    
	$gtypes{$alleles}++;
    }

    # print STDERR "Num Hom Gtypes: $hom_cnt; Num Het GTypes: $het_cnt\n";

    #
    # Is the parent heterozygous or homozygous
    #
    $par_cnt = scalar(@{$parent});

    if ($par_cnt == 1) {

	if ($hom_cnt == 0 && $het_cnt == 2) {
	    $marker = "cc/ab";
	} elsif ($hom_cnt == 0 && $het_cnt == 1) {
	    $marker = "aa/bb";
	} elsif ($hom_cnt == 2 && $het_cnt == 1) {
	    $marker = "aa/ab";
	}

    } elsif ($par_cnt == 2) {

	if ($hom_cnt == 0 && $het_cnt == 2) {
	    $marker = "ab/cc";
	} elsif ($hom_cnt == 2 && $het_cnt == 0) {
	    $marker = "ab/--";
	} elsif ($hom_cnt == 1 && $het_cnt == 1) {
	    $marker = "ab/aa";
	} elsif ($hom_cnt == 2 && $het_cnt == 1) {
	    $marker = "ab/ab";
	} elsif ($hom_cnt == 1 && $het_cnt == 2) {
	    $marker = "ab/ac";
	} elsif ($hom_cnt == 0 && $het_cnt == 4) {
	    $marker = "ab/cd";
	}
    }

    #
    # Check to make sure error alleles aren't interfering with the marker imputation
    #
    if ($marker eq "ab/ab") {
	my %types = ("ab/ab" => [.25, .25, .50], 
		     "ab/aa" => [.50, .50]);
	my $min_err  = 1000;
	my $min_type = "";

	foreach $type (keys %types) {
	    $err = cmp_hardy_weinberg($type, $types{$type}, \%gtypes);
	    # print STDERR "Type: $type, Fit: $err\n";
	    if ($err < $min_err) {
		$min_err = $err;
		$min_type = $type;
	    }
	}

	return $min_type;
    }

    return $marker;
}

sub cmp_hardy_weinberg {
    my ($type, $pcts, $gtypes) = @_;

    my ($gtype, $pct, $tot, $i, $j, $cnt, $obs_cnt, $diff, $key);
    my (@hw, @keys);

    $tot = 0;
    foreach $gtype (keys %{$gtypes}) {
	$tot += $gtypes->{$gtype};
    }

    #
    # Create a two-dimensional matrix comparing the expected Hardy-Weinberg frequencies
    # to the observed frequencies.
    #
    $cnt = scalar(@{$pcts}) - 1;
    foreach $i (0..$cnt) {
	$hw[$i] = [];

	foreach $gtype (keys %{$gtypes}) {
	    push(@{$hw[$i]}, abs($pcts->[$i] - $gtypes->{$gtype}/$tot));

	    # print STDERR 
	    # 	"Gtype: $gtype; ",
	    # 	"Cnt: ", $gtypes->{$gtype}/$tot, "; ",
	    # 	"Diff: ", abs($pcts->[$i] - $gtypes->{$gtype}/$tot), "\n";
	}
	# print STDERR "\n";
    }

    my (%hw_map, %hw_diff, %used);
    #
    # Select the best fit betwee each observed genotyped and hardy-weinberg frequency
    #
    @keys    = keys %{$gtypes};
    $obs_cnt = scalar(@keys) - 1;
    foreach $i (0..$cnt) {
	#print STDERR "Looking at HW Pct: $pcts->[$i]\n";

	foreach $j (0..$obs_cnt) {
	    $diff  = $hw[$i][$j];
	    $gtype = $keys[$j];

	    #print STDERR "  Diff for genotype $gtype: $diff\n";

	    if (!defined($used{$gtype}) &&
		(!defined($hw_diff{$i}) || $diff < $hw_diff{$i})) {
		#print STDERR "    Found a better fit\n";
		$hw_diff{$i} = $diff;
		$hw_map{$i}  = $gtype;
	    }
	}
	$used{$hw_map{$i}}++;
    }

    my $tot_err = 0;

    foreach $key (keys %hw_map) {
	#print STDERR "Mapped $key ($pcts->[$key]) to $hw_map{$key}\n";
	$tot_err += $hw_diff{$key};
    }

    return $tot_err;
}

sub find_markers {
    my ($sth, $catalog, $order) = @_;

    my ($key, $k, $allele_cnt_1, $allele_cnt_2, $annote, $sample_id, $tag_id, $parent_count, $tag_count, $ratio);

    open(OUT, ">$out_path/batch_" . $batch_id . ".markers.tsv") or 
	die("Unable to open output file: $out_path/batch_" . $batch_id . ".markers.tsv; $!\n");

    foreach $key (keys %{$catalog}) {
	#dump_tags($catalog->{$key}) if ($debug);

	my (%samples);
	#
	# Count the number of parental tags matching this catalog tag. A proper marker should
	# contain a single representative from each parent; multiple alleles must be called from 
	# a single tag from a single parent.
	#
	my @parents = keys %{$catalog->{$key}->{'par'}};

	foreach $k (@parents) {
	    ($sample_id, $tag_id) = split("_", $k);
	    $samples{$sample_id}++;
	}

	$parent_count = scalar(keys %samples);
	$tag_count    = scalar(@parents);

	if ($parent_count == 1 && $impute_markers == true) {
	    $annote = impute_marker($catalog->{$key}->{'par'}->{$parents[0]}, $catalog->{$key}->{'pro'});

	    if ($annote ne "unknown") {
		$ratio  = tally_progeny_alleles($order,
						$catalog->{$key}->{'par'}, 
						$catalog->{$key}->{'pro'},
						$key, $annote);
		print OUT
		    "0\t",
		    $batch_id, "\t",
		    $key, "\t",
		    $annote, "\t",
		    $ratio, "\n";
	    }
	    next;
	}

	#
	# Locus is present in both parents.
	#
	if ($parent_count == 2 && $tag_count == 2) {
	    $allele_cnt_1 = scalar(@{$catalog->{$key}->{'par'}->{$parents[0]}});
	    $allele_cnt_2 = scalar(@{$catalog->{$key}->{'par'}->{$parents[1]}});

            #
            # Determine the number of unique alleles
            #
            my (%unique_alleles, $num_unique_alleles, $allele);

            foreach $allele (@{$catalog->{$key}->{'par'}->{$parents[0]}}, 
                             @{$catalog->{$key}->{'par'}->{$parents[1]}}) {
                $unique_alleles{$allele}++;
            }
            $num_unique_alleles = scalar(keys %unique_alleles);

	    #
	    # Locus is heterozygous in both parents. However, the number of alleles present distinguishes 
            # what type of marker it is. Four unique alleles requries an ab/cd marker, while four 
            # alleles that are the same in both parents requires an ab/ab marker. Finally, three unique 
            # alleles requires either an ab/ac marker.
	    #
	    if ($allele_cnt_1 == 2 && $allele_cnt_2 == 2) {
                
                if ($num_unique_alleles == 3) {
                    $annote = "ab/ac";

                } elsif ($num_unique_alleles == 2) {
                    $annote = "ab/ab";

                } else {
                    $annote = "ab/cd";
                }

		$ratio  = tally_progeny_alleles($order,
						$catalog->{$key}->{'par'}, 
						$catalog->{$key}->{'pro'},
						$key, $annote);
		print OUT
		    "0\t",
		    $batch_id, "\t",
		    $key, "\t",
		    $annote, "\t",
		    $ratio, "\n";
	    #
	    # Rad-Tag is homozygous in one parent and heterozygous in the other.
	    #
	    } elsif ($allele_cnt_1 == 2 && $allele_cnt_2 == 1) {

                ($sample_id, $tag_id) = split("_", $parents[0]);

                if ($num_unique_alleles == 3) {
                    $annote = $order->{$sample_id} eq 'first' ? "ab/cc" : "cc/ab";

                } elsif ($num_unique_alleles == 2) {
                    $annote = $order->{$sample_id} eq 'first' ? "ab/aa" : "aa/ab";
                }

		$ratio = tally_progeny_alleles($order,
					       $catalog->{$key}->{'par'}, 
					       $catalog->{$key}->{'pro'},
					       $key, $annote);
		print OUT
		    "0\t",
		    $batch_id, "\t",
		    $key, "\t",
		    $annote, "\t",
		    $ratio, "\n";
	    #
	    # Rad-Tag is homozygous in one parent and heterozygous in the other.
	    #
	    } elsif ($allele_cnt_1 == 1 && $allele_cnt_2 == 2) {

                ($sample_id, $tag_id) = split("_", $parents[0]);

                if ($num_unique_alleles == 3) {
                    $annote = $order->{$sample_id} eq 'first' ? "cc/ab" : "ab/cc";

                } elsif ($num_unique_alleles == 2) {
                    $annote = $order->{$sample_id} eq 'first' ? "aa/ab" : "ab/aa";
                }

		$ratio = tally_progeny_alleles($order,
					       $catalog->{$key}->{'par'}, 
					       $catalog->{$key}->{'pro'},
					       $key, $annote);
		print OUT
		    "0\t",
		    $batch_id, "\t",
		    $key, "\t",
		    $annote, "\t",
		    $ratio, "\n";
            #
            # RAD-Tag is homozygous in both parents, but heterozygous between parents.
            #
	    } elsif ($allele_cnt_1 == 1 && $allele_cnt_2 == 1) {

                if ($catalog->{$key}->{'par'}->{$parents[0]}->[0] ne "consensus" &&
                    $catalog->{$key}->{'par'}->{$parents[1]}->[0] ne "consensus") {
                    $annote = "aa/bb";
                    $ratio  = tally_progeny_alleles($order,
						    $catalog->{$key}->{'par'}, 
                                                    $catalog->{$key}->{'pro'},
                                                    $key, $annote);
                    print OUT
                        "0\t",
                        $batch_id, "\t",
                        $key, "\t",
                        $annote, "\t",
                        $ratio, "\n";
                }
	    }
        #
        # Rad-Tag only exists in one parent.
	#
	} elsif ($parent_count == 1 && $tag_count == 1) {
	    $allele_cnt_1 = scalar(@{$catalog->{$key}->{'par'}->{$parents[0]}});

	    if ($allele_cnt_1 == 2) {
		($sample_id, $tag_id) = split("_", $parents[0]);
		$annote = $order->{$sample_id} eq 'first' ? "ab/--" : "--/ab";

		$ratio = tally_progeny_alleles($order,
					       $catalog->{$key}->{'par'}, 
					       $catalog->{$key}->{'pro'},
					       $key, $annote);

		print OUT
		    "0\t",
		    $batch_id, "\t",
		    $key, "\t",
		    $annote, "\t",
		    $ratio, "\n";
	    }
	}
    }

    close(OUT);
}

sub compare_snp_locations {
    my ($sth, $p_1, $p_2) = @_;

    my ($row, $sample_id, $tag_id, $col_1, $col_2);

    ($sample_id, $tag_id) = split("_", $p_1);

    $sth->{'snp'}->execute($sample_id, $tag_id)
	or die("Unable to select results from $db.\n");
    $row = $sth->{'snp'}->fetchrow_arrayref();
    $col_1 = $row->[0];

    ($sample_id, $tag_id) = split("_", $p_2);

    $sth->{'snp'}->execute($sample_id, $tag_id)
	or die("Unable to select results from $db.\n");
    $row = $sth->{'snp'}->fetchrow_arrayref();
    $col_2 = $row->[0];

    if ($col_1 == $col_2) {
	return 1;
    } else {
	return 0;
    }
}

sub tally_progeny_alleles {
    my ($order, $parents, $progeny, $tag_id, $marker) = @_;

    my (@keys, %markers, %genotype_map);
    my ($key, $allele, $m, $pct);
    
    #
    # Create a map between alleles and genotype symbols
    #
    create_genotype_map($marker, $tag_id, $order, $parents, \%genotype_map);

    my $dictionary = {'ab/--' => {'a'  => 'aa',
				  'b'  => 'bb'},
		      '--/ab' => {'a'  => 'aa',
				  'b'  => 'bb'},
		      'aa/bb' => {'a'  => 'aa',
				  'ab' => 'ab',
				  'b'  => 'bb'},
		      'ab/cd' => {'a'  => 'aa',
				  'ab' => 'ab',
				  'b'  => 'bb',
				  'c'  => 'cc',
				  'cd' => 'cd',
				  'd'  => 'dd',
				  'ac' => 'ac',
				  'ad' => 'ad',
				  'bc' => 'bc',
				  'bd' => 'bd'},
		      'ab/aa' => {'a'  => 'aa',
				  'ab' => 'ab',
				  'b'  => 'bb'},
		      'aa/ab' => {'a'  => 'aa',
				  'ab' => 'ab',
				  'b'  => 'bb'},
		      'ab/cc' => {'a'  => 'aa',
				  'ab' => 'ab',
				  'bb' => 'bb',
				  'c'  => 'cc',
				  'ac' => 'ac',
				  'bc' => 'bc'},
		      'cc/ab' => {'aa' => 'aa',
				  'ab' => 'ab',
				  'bb' => 'bb',
				  'c'  => 'cc',
				  'ac' => 'ac',
				  'bc' => 'bc'},
		      'ab/ab' => {'a'  => 'aa',
				  'b'  => 'bb',
				  'ab' => 'ab'}};

    foreach $key (keys %{$progeny}) {
	my @alleles;

	print STDERR "Examining progeny $key; marker: $marker\n" if ($debug);

	#
	# If there is more than a single tag from a particular progeny, matching this tag in the 
	# catalog, then discard this progeny, since we don't know which tag is correct and the
	# catalog tag is probably overmerged.
	#
	@keys = keys %{$progeny->{$key}};
	if (scalar(@keys) > 1) {
	    print STDERR 
		"Discarding progeny $key from catalog tag $tag_id with multiple tag matches: ", 
		join(" ", @keys), "\n" if ($debug);
	    next;
	}

	foreach $allele (@{$progeny->{$key}->{$keys[0]}}) {
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
	$m = join("", sort @alleles);
	print STDERR "Converting from '$m' to ", $dictionary->{$marker}->{$m}, "\n" if ($debug);

	$m = defined($dictionary->{$marker}->{$m}) ? $dictionary->{$marker}->{$m} : "-";


	if (!defined($genotypes{$marker}->{$m})) {
	    print STDERR "  Warning: illegal genotype encountered ('$m') in tag $tag_id, progeny $key\n" if ($debug);
	} else  {
	    $markers{$m}++;
	    print STDERR "  $m allele: $markers{$m}\n" if ($debug);
	}
    }

    my $ratio         = "";
    my $valid_progeny = 0;
    my $max_pct       = 0;

    foreach $m (sort keys %{$genotypes{$marker}}) {
	$valid_progeny += $markers{$m};
    }
    foreach $m (sort keys %{$genotypes{$marker}}) {
	$pct = $valid_progeny > 0 ? sprintf("%.1f", $markers{$m} / $valid_progeny * 100) : "0.0";
	if ($markers{$m} > 0) {
	    $ratio .= defined($markers{$m}) ? "$m:$markers{$m}($pct%);" : "$m:0(0%);";
	}
	$max_pct = $pct > $max_pct ? $pct : $max_pct;
    }

    return $valid_progeny . "\t" . $max_pct . "\t" . $ratio;
}

sub create_genotype_map {
    my ($marker, $tag_id, $order, $parents, $map) = @_;

    my (%genotypes, %legal_genotypes, @parents, @types, @keys, %alleles, %par_specific, %com_types);
    my ($type, $m, $key, $allele, $sample_id, $tag_id);
    #
    # Create a genotype map. For any set of alleles, this routine will
    # assign each allele to one of the constituent genotypes, e.g. given the 
    # marker type 'aaxbb' and the alleles 'A' from the male, and 'G'
    # from the female, will assign 'G' == 'bb' and 'A'== 'aa'. It assumes that 
    # recombination may have occurred as with an F2, F3 or later cross.
    #
    @parents = keys %{$parents};

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
    ($sample_id, $tag_id) = split("_", $parents[0]);
    $key = $order->{$sample_id} eq "first" ? $parents[0] : $parents[1];
    $m   = substr($marker, 0, 2);

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
	    print STDERR "  Assinging '$allele' to genotype '", $map->{$allele}, "'\n" if ($debug);
	}
    }

    #
    # Finally, repeat in the second parent.
    #
    %legal_genotypes = ();
    $key = $order->{$sample_id} eq "second"  ? $parents[0] : $parents[1];
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
	    print STDERR "  Assinging '$allele' to genotype '", $map->{$allele}, "'\n" if ($debug);
	}
    }

}

sub populate {
    my ($sth, $catalog, $order) = @_;

    my ($row, $key, $parents, $progeny);

    $sth->{'match'}->execute($batch_id)
	or die("Unable to select results from $db.\n");

    #
    # For each catalog tag, count the number of alleles present
    # in each of the parents and the progeny.
    #
    while ($row = $sth->{'match'}->fetchrow_hashref()) {

	if (!defined($catalog->{$row->{'catalog_id'}})) {
	    $catalog->{$row->{'catalog_id'}} = {};
	    $catalog->{$row->{'catalog_id'}}->{'par'} = {};
	    $catalog->{$row->{'catalog_id'}}->{'pro'} = {};
	}

	if ($row->{'type'} eq "parent") {

	    $parents = $catalog->{$row->{'catalog_id'}}->{'par'};
	    $key     = $row->{'sample_id'} . "_" . $row->{'tag_id'};

	    if (!defined($parents->{$key})) {
		$parents->{$key} = [];
	    }
	    push(@{$parents->{$key}}, $row->{'allele'});

	} elsif ($row->{'type'} eq "progeny") {
	    $progeny = $catalog->{$row->{'catalog_id'}}->{'pro'};

	    if (!defined($progeny->{$row->{'sample_id'}})) {
		$progeny->{$row->{'sample_id'}} = {};
	    }
	    if (!defined($progeny->{$row->{'sample_id'}}->{$row->{'tag_id'}})) {
		$progeny->{$row->{'sample_id'}}->{$row->{'tag_id'}} = [];
	    }
	    push(@{$progeny->{$row->{'sample_id'}}->{$row->{'tag_id'}}}, $row->{'allele'});

	} else {
	    print STDERR "Unknown catalog type '$row->{'type'}'\n";
	}
    }

    my @o = ('first', 'second');

    $sth->{'order'}->execute($batch_id)
	or die("Unable to select results from $db.\n");

    my $num_rows = $sth->{'order'}->rows();

    while (@o) {
	$row = $sth->{'order'}->fetchrow_hashref();
	$order->{$row->{'id'}} = shift(@o) if (defined($row));

	print STDERR "Order $row->{'id'}: ", $order->{$row->{'id'}}, "\n" if ($debug);
	last if ($num_rows == 1);
    }
}

sub dump_tags {
    my ($catalog) = @_;

    my ($key, $sample_id, $tag_id);

    print STDERR 
	"Examining catalog ID $key\n",
	"    ", scalar(keys %{$catalog->{'par'}}), " parents:\n";
    foreach $key (keys %{$catalog->{'par'}}) {
	print STDERR 
	    "      ", 
	    $key, ": ", 
	    scalar(@{$catalog->{'par'}->{$key}}), " alleles.\n";
    }

    print STDERR
	"    ", scalar(keys %{$catalog->{'pro'}}), " progeny:\n";
    foreach $sample_id (keys %{$catalog->{'pro'}}) {
	foreach $tag_id (keys %{$catalog->{'pro'}->{$sample_id}}) {
	    print STDERR 
		"      ", 
		$sample_id, "_", $tag_id, ": ", 
		scalar(@{$catalog->{'pro'}->{$sample_id}->{$tag_id}}), " alleles.\n";
	}
    }
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
	"SELECT catalog_id, matches.sample_id, file, type, tag_id, allele " . 
	"FROM matches " .
	"JOIN samples ON (matches.sample_id=samples.id) " .
	"WHERE matches.batch_id=?";
    $sth->{'match'} = $sth->{'dbh'}->prepare($query) or die($sth->{'dbh'}->errstr());

    $query =
	"SELECT col FROM snps WHERE sample_id=? AND tag_id=?";
    $sth->{'snp'} = $sth->{'dbh'}->prepare($query) or die($sth->{'dbh'}->errstr());

    $query =
	"SELECT id FROM samples WHERE type='parent' AND batch_id=? ORDER BY id";
    $sth->{'order'} = $sth->{'dbh'}->prepare($query) or die($sth->{'dbh'}->errstr());
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
	elsif ($_ =~ /^-D$/) { $db       = shift @ARGV; }
	elsif ($_ =~ /^-b$/) { $batch_id = shift @ARGV; }
	elsif ($_ =~ /^-o$/) { $out_path = shift @ARGV; }
	elsif ($_ =~ /^-i$/) { $impute_markers = true; }
	elsif ($_ =~ /^-v$/) { version(); exit(); }
	elsif ($_ =~ /^-h$/) { usage(); }
	else {
	    print STDERR "Unknown command line options received: $_\n";
	    usage();
	}
    }

    if ($batch_id == 0) {
	print STDERR "You must specify a batch ID.\n";
	usage();
    }

    if (length($db) == 0) {
        print STDERR "You must specify a database to access.\n";
        usage();
    }
}

sub version {
    print STDERR "markers.pl ", stacks_version, "\n";
}

sub usage {
    version();

    print << "EOQ";
markers.pl -b id [-D db] [-o path] [-d] [-h]
  D: radtag database to examine.
  b: Batch ID to examine.
  o: path to print output file.
  h: display this help message.
  d: turn on debug output.

EOQ
    exit(0);
}
