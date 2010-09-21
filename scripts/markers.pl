#!/usr/bin/perl
#
# Written by Julian Catchen <catchen@cs.uoregon.edu>
#

use strict;
use DBI;

my $debug     = 0;
my $db        = "";
my $batch_id  = 0;
my $out_path  = "./";

my %genotypes = (
		 'ab/cc' => {'ac' => 0, 'bc' => 0},
		 'cc/ab' => {'ac' => 0, 'bc' => 0},
		 'aa/bb' => {'aa' => 0, 'ab' => 0, 'bb' => 0},
		 'ab/--' => {'aa' => 0, 'ab' => 0},
		 '--/ab' => {'aa' => 0, 'ab' => 0},
		 'ab/aa' => {'aa' => 0, 'ab' => 0},
		 'aa/ab' => {'aa' => 0, 'ab' => 0},
		 'ab/ab' => {'aa' => 0, 'ab' => 0, 'bb' => 0},
		 'ab/ac' => {'aa' => 0, 'ab' => 0, 'ac' => 0, 'bc' => 0},
		 'ab/cd' => {'ac' => 0, 'ad' => 0, 'bc' => 0, 'bd' => 0}
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

sub find_markers {
    my ($sth, $catalog, $order) = @_;

    my ($key, $k, $allele_cnt_1, $allele_cnt_2, $annote, $sample_id, $tag_id, $parent_count, $tag_count, $ratio);

    open(OUT, ">$out_path/batch_" . $batch_id . ".markers.tsv") or 
	die("Unable to open output file: $out_path/batch_" . $batch_id . ".markers.tsv; $!\n");

    foreach $key (keys %{$catalog}) {
	dump_tags($catalog->{$key}) if ($debug);

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

	#
	# RAD-Tag is present in both parents.
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
	    # Rad-Tag is heterozygous in both parents. However, the number of alleles present distinguishes 
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

		$ratio  = tally_progeny_alleles($catalog->{$key}->{'par'}, 
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

                if ($num_unique_alleles == 3) {
                    $annote = "ab/cc";

                } elsif ($num_unique_alleles == 2) {
                    ($sample_id, $tag_id) = split("_", $parents[0]);
                    $annote = $order->{$sample_id} eq 'first' ? "ab/aa" : "aa/ab";
                }

		$ratio = tally_progeny_alleles($catalog->{$key}->{'par'}, 
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

                if ($num_unique_alleles == 3) {
                    $annote = "cc/ab";

                } elsif ($num_unique_alleles == 2) {
                    ($sample_id, $tag_id) = split("_", $parents[0]);
                    $annote = $order->{$sample_id} eq 'first' ? "ab/aa" : "aa/ab";
                }

		$ratio = tally_progeny_alleles($catalog->{$key}->{'par'}, 
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
                    $ratio  = tally_progeny_alleles($catalog->{$key}->{'par'}, 
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

		$ratio = tally_progeny_alleles($catalog->{$key}->{'par'}, 
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
    my ($parents, $progeny, $tag_id, $marker) = @_;

    my (@keys, %markers, %genotype_map);
    my ($key, $allele, $m, $pct);
    
    #
    # Create a map between alleles and genotype symbols
    #
    create_marker_map($marker, $tag_id, $parents, \%genotype_map);

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
	@alleles = sort @alleles;

	#
	# If the tag was male or female-only in the parents, than we can no have more 
	# than a single allele in any individual progeny.
	# 
	if (scalar(@alleles) == 2 &&
	    ($marker eq "ab/--" || $marker eq "--/ab")) {
	    # Illegal genotype
	    splice(@alleles, 0);

	#
	# In the cases where we have a common genotype between the parents, such as 'aa'
	# the pipeline will only report a single 'a', since all the RAD-Tags for the two
	# were merged into a single allele by the pipeline.
	#
	} elsif (scalar(@alleles) == 1 && $alleles[0] eq 'a') {
	    push(@alleles, $alleles[0]);

	#
	# We are mapping ab/ab tags.
	#
	} elsif (scalar(@alleles) == 1 && $marker eq "ab/ab" && $alleles[0] eq 'b') {
	    push(@alleles, $alleles[0]);

	#
	# We are mapping ab/-- and --/ab tags as lmxll and nnxnp, respectively, since 
	# Joinmap does not have a specification for heterozygous, single parent tags.
	#
	} elsif (scalar(@alleles) == 1 && 
		 ($marker eq "ab/--" || $marker eq "--/ab")) {
	    unshift(@alleles, 'a') if ($alleles[0] eq 'b');
	}

	$m = join("", @alleles);

	if (!defined($genotypes{$marker}->{$m})) {
	    print STDERR "  Warning: illegal genotype encountered in tag $tag_id, progeny $key\n" if ($debug);
	} else  {
	    $markers{$m}++;
	    print STDERR "  $key allele: $markers{$key}\n" if ($debug);
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
	$ratio .= defined($markers{$m}) ? "$m:$markers{$m}($pct%);" : "$m:0(0%);";
	$max_pct = $pct > $max_pct ? $pct : $max_pct;
    }

    return $valid_progeny . "\t" . $max_pct . "\t" . $ratio;
}

sub create_marker_map {
    my ($marker, $tag_id, $parents, $map) = @_;

    my (%genotypes, @nucs, %alleles);
    my ($nuc, $m, $key, $allele);

    #
    # Create a genotype map. For any set of alleles, this routine will
    # assign each allele to one of the constituent markers, e.g. given the 
    # marker type 'abxaa' and the alleles 'A' and 'G' from the male, and 'G'
    # from the female, will assign 'G' == 'a' and 'A' = 'b'.
    #
    $m = substr($marker, 0, 2) . substr($marker, 3, 2);

    foreach $nuc (split(//, $m)) {
	next if ($nuc eq "-");
	$genotypes{$nuc}++;
    }

    @nucs = sort keys %genotypes;

    foreach $key (keys %{$parents}) {
	foreach $allele (@{$parents->{$key}}) {
	    $alleles{$allele}++;
	}
    }

    foreach $allele (sort {$alleles{$b} <=> $alleles{$a}} keys %alleles) {

	if (!defined($map->{$allele})) {
	    $key = shift @nucs;

	    #die("Impossible allele combination!\n") if (!defined($key));
	    if (!defined($key)) {
		print STDERR "Warning: impossible allele combination in parental genotypes, tag $tag_id\n";
		next;
	    }   
	    $map->{$allele} = $key;
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

    # Connect to the database, assumes user has a MySQL ~/.my.cnf file to
    # specify the host, username and password
    $sth->{'dbh'} = DBI->connect("DBI:mysql:$db:mysql_read_default_file=" . $ENV{"HOME"} . "/.my.cnf")
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
}

sub usage {
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
