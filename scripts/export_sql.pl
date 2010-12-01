#!/usr/bin/perl -w
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
# Read in a set of filtering paramters, query a Stacks pipeline database based on
# those filters and write the results into a compact tab-separated values or excel file.
#
# Written by Julian Catchen <jcatchen@uoregon.edu>
#

use strict;
use DBI;
#use Excel::Writer::XLSX;
use Spreadsheet::WriteExcel;

use constant stacks_version => "_VERSION_";

my $mysql_config = "_PKGDATADIR_" . "sql/mysql.cnf";
my $out_file = "";
my $type     = "tsv";
my $batch_id = 0;
my $db       = "";
my $debug    = 0;

parse_command_line();

my (%sth, %loci, %samples, %filters);

prepare_sql_handles(\%sth, \%filters);

populate(\%sth, \%loci, \%samples, \%filters);

write_results(\%loci, \%samples);

print "Success\n";

sub populate {
    my ($sth, $loci, $samples, $filters) = @_;

    my (%delev, $row, $snp_row, $all_row, $gen_row);

    my @params;

    #
    # Cache all the stacks that were deleveraged.
    #
    $sth->{'delev'}->execute($batch_id);

    while ($row = $sth->{'delev'}->fetchrow_hashref()) {
        $delev{$row->{'sample_id'} . "_" . $row->{'tag_id'}}++;
    }

    #
    # Pull list of samples for this batch
    #
    $sth->{'samp'}->execute($batch_id);

    while ($row = $sth->{'samp'}->fetchrow_hashref()) {
        $samples->{$row->{'file'}} = $row->{'id'};
    }

    prepare_filter_parameters(\@params, $filters);

    #
    # Fetch the results and populate the array of groups.
    #
    $sth->{'tag'}->execute(@params);

    while ($row = $sth->{'tag'}->fetchrow_hashref()) {
        my $locus = {};
        $locus->{'id'}         = $row->{'tag_id'};
        $locus->{'annotation'} = defined($row->{'external_id'}) ? $row->{'external_id'} : "";
        $locus->{'chr'}        = $row->{'chr'};
        $locus->{'bp'}         = $row->{'bp'};
        $locus->{'delev'}      = 0;
        $locus->{'marker'}     = $row->{'marker'};
        $locus->{'seq'}        = $row->{'seq'};
        $locus->{'alleles'}    = "";
        $locus->{'snps'}       = "";
        $locus->{'gtypes'}     = {};

        $locus->{'num_alleles'}   = $row->{'alleles'};
        $locus->{'num_snps'}      = $row->{'snps'};
        $locus->{'num_parents'}   = $row->{'parents'};
        $locus->{'num_progeny'}   = $row->{'progeny'};
        $locus->{'valid_progeny'} = $row->{'valid_progeny'};
        $locus->{'num_ests'}      = $row->{'ests'};
        $locus->{'num_pe_tags'}   = $row->{'pe_radtags'};
        $locus->{'num_blast'}     = $row->{'blast_hits'};

        #
        # Fetch SNPs and Alleles
        #
        $sth->{'snp'}->execute($batch_id, $locus->{'id'});

        while ($snp_row = $sth->{'snp'}->fetchrow_hashref()) {
            $locus->{'snps'} .= $snp_row->{'col'} . "," . $snp_row->{'rank_1'} . ">" . $snp_row->{'rank_2'} . ";";
        }
        $locus->{'snps'} = substr($locus->{'snps'}, 0, -1) if (length($locus->{'snps'}) > 0);

        $sth->{'allele'}->execute($batch_id, $locus->{'id'});

        while ($all_row = $sth->{'allele'}->fetchrow_hashref()) {
            $locus->{'alleles'} .= $all_row->{'allele'} . ";";
        }
        $locus->{'alleles'} = substr($locus->{'alleles'}, 0, -1) if (length($locus->{'alleles'}) > 0);

        #
        # Add genotypes
        #
        $sth->{'mat'}->execute($batch_id, $locus->{'id'});

        while ($gen_row = $sth->{'mat'}->fetchrow_hashref()) {
            if (!defined($locus->{'gtypes'}->{$gen_row->{'file'}})) {
                $locus->{'gtypes'}->{$gen_row->{'file'}} = [];
            }

            push(@{$locus->{'gtypes'}->{$gen_row->{'file'}}}, 
                 {'file'   => $gen_row->{'file'},
                  'allele' => $gen_row->{'allele'},
                  'tag_id' => $gen_row->{'tag_id'}});

            #
            # Check if this particular sample was deleveraged
            #
            if (defined($delev{$gen_row->{'id'} . "_" . $gen_row->{'tag_id'}}) &&
                $delev{$gen_row->{'id'} . "_" . $gen_row->{'tag_id'}} >= 1) {
                $locus->{'delev'}++;
            }
        }

        $loci->{$row->{'tag_id'}} = $locus;
    }
}

sub prepare_filter_parameters {
    my ($params, $filters) = @_;

    my ($filter);

    push(@{$params}, $batch_id);

    foreach $filter (keys %{$filters}) {

        if ($filter eq "snps") {
            push(@{$params}, $filters->{'snps'});

        } elsif ($filter eq "alle") {
            push(@{$params}, $filters->{'alle'});

        } elsif ($filter eq "pare") {
            push(@{$params}, $filters->{'pare'});
	
        } elsif ($filter eq "prog") {
            push(@{$params}, $filters->{'prog'});

        } elsif ($filter eq "vprog") {
            push(@{$params}, $filters->{'vprog'});

        } elsif ($filter eq "cata") {
            push(@{$params}, $filters->{'cata'});

        } elsif ($filter eq "est") {
            push(@{$params}, 0);

        } elsif ($filter eq "pe") {
            push(@{$params}, 0);

        } elsif ($filter eq "blast") {
            push(@{$params}, 0);

        } elsif ($filter eq "mark") {
            if ($filters->{'mark'} eq "Any") {
                push(@{$params}, "%/%");
            } else {
                push(@{$params}, $filters->{'mark'});
            }
        }
    }
}

sub apply_query_filters {
    my ($filters) = @_;

    my ($query, $filter) = "";

    my %sql_filters = 
        ("cata"  => "(catalog_index.tag_id = ?)", 
         "alle"  => "(alleles >= ?)", 
         "snps"  => "(snps >= ?)",
         "pare"  => "(parents = ?)",
         "prog"  => "(progeny >= ?)",
         "vprog" => "(valid_progeny >= ?)",
         "mark"  => "(marker LIKE ?)", 
         "est"   => "(ests > ?)",
         "pe"    => "(pe_radtags > ?)",
         "blast" => "(blast_hits > ?)");
    
    if (scalar(keys %{$filters}) > 0) {

        foreach $filter (keys %{$filters}) {
            $query .= " AND ";
            $query .= $sql_filters{$filter};
        }
    }

    return $query;
}

sub write_results {
    my ($loci, $samples) = @_;

    my ($workbook, $worksheet);

    my ($out_fh, $str, $cat_id, $id, $locus, $gtypes, $types);

    if ($type eq "xls") {
        $workbook  = Spreadsheet::WriteExcel->new($out_file) or die("Unable to initiate excel spreadsheet.\n");
        $worksheet = $workbook->add_worksheet() or die("Unable to add a worksheet to our excel spreadsheet.\n");
    } else {
        open($out_fh, ">$out_file") or die("Unable to open output file '$out_file'\n");
    }

    #
    # Print the heading out for the spreadsheet
    #
    my $i = 0;

    $str = 
        "Catalog ID\t" . 
        "Annotation\t" . 
        "Chr\t" .
        "BP\t" .
        "Consensus Sequence\t" .
        "Num Parents\t" .
        "Num Progeny\t" . 
        "Num SNPs\t" .
        "SNPs\t" .
        "Num Alleles\t" .
        "Alleles\t" .
        "Deleveraged\t";

    foreach $id (keys %{$samples}) {
        $str .= $id . "\t";
    }
    $str  = substr($str, 0, -1);
    $str .= "\n";

    $type eq "xls" ? write_excel($worksheet, $i, $str) : print $out_fh $str;
    $i++;

    foreach $cat_id (keys %{$loci}) {
        $locus = $loci->{$cat_id};

        $str =
            $cat_id . "\t" .
            $locus->{'annotation'} . "\t" .
            $locus->{'chr'} . "\t" .
            $locus->{'bp'} . "\t" .
            $locus->{'seq'} . "\t" .
            $locus->{'num_parents'} . "\t" .
            $locus->{'num_progeny'} . "\t" .
            $locus->{'num_snps'} . "\t" .
            $locus->{'snps'} . "\t" . 
            $locus->{'num_alleles'} . "\t" .
            $locus->{'alleles'} . "\t" .
            $locus->{'delev'} . "\t";

        foreach $id (keys %{$samples}) {
            $types = $locus->{'gtypes'}->{$id};

            if (!defined($types)) {
                $str .= "\t";
                next;
            }

            foreach $type (@{$types}) {
                $str .= $type->{'allele'} . "/";
            }

            $str  = substr($str, 0, -1);
            $str .= "\t";
        }
        $str  = substr($str, 0, -1);
        $str .= "\n";

        $type eq "xls" ? write_excel($worksheet, $i, $str) : print $out_fh $str;

        $i++;
    }

    $str = "\n";
    $type eq "xls" ? write_excel($worksheet, $i, $str) : print $out_fh $str;
    $i++;

    foreach $id (keys %{$samples}) {
        $str = $samples->{$id} . "\t" . $id . "\n";
        $type eq "xls" ? write_excel($worksheet, $i, $str) : print $out_fh $str;
        $i++;
    }

    $type eq "xls" ? $workbook->close() : close($out_fh);
}

sub write_excel {
    my ($worksheet, $i, $str) = @_;

    chomp $str;
    my @row = split(/\t/, $str);

    my $j = 0;

    foreach my $r (@row) {
        $worksheet->write($i, $j, $r);
        $j++;
    }
}

sub prepare_sql_handles {
    my ($sth, $filters) = @_;

    #
    # Connect to the database, check for the existence of a MySQL config file in the home
    # directory first, otherwise use the stacks-distributed one.
    #
    my $cnf = (-e $ENV{"HOME"} . "/.my.cnf") ? $ENV{"HOME"} . "/.my.cnf" : $mysql_config;
    $sth->{'dbh'} = DBI->connect("DBI:mysql:$db:mysql_read_default_file=$cnf")
	or die("Unable to connect to the $db MySQL Database!\n" . $DBI::errstr);

    my $query;

    $query = 
        "SELECT sample_id, tag_id FROM tag_index " . 
        "WHERE batch_id=? AND deleveraged=true";
    $sth->{'delev'} = $sth->{'dbh'}->prepare($query) or die($sth->{'dbh'}->errstr());

    $query = 
        "SELECT id, file FROM samples " . 
        "WHERE batch_id=? ORDER BY id";
    $sth->{'samp'} = $sth->{'dbh'}->prepare($query) or die($sth->{'dbh'}->errstr());

    $query = 
        "SELECT allele FROM catalog_alleles " . 
        "WHERE batch_id=? AND tag_id=?";
    $sth->{'allele'} = $sth->{'dbh'}->prepare($query) or die($sth->{'dbh'}->errstr());

    $query = 
        "SELECT col, rank_1, rank_2 FROM catalog_snps " . 
        "WHERE batch_id=? AND tag_id=? ORDER BY col";
    $sth->{'snp'} = $sth->{'dbh'}->prepare($query) or die($sth->{'dbh'}->errstr());

    $query = 
        "SELECT samples.id as id, samples.sample_id, samples.type, file, tag_id, allele " . 
        "FROM matches " . 
        "JOIN samples ON (matches.sample_id=samples.id) " . 
        "WHERE matches.batch_id=? AND catalog_id=? ORDER BY samples.id";
    $sth->{'mat'} = $sth->{'dbh'}->prepare($query) or die($sth->{'dbh'}->errstr());

    $query = 
        "SELECT catalog_index.tag_id as tag_id, chr, bp, snps, alleles, parents, progeny, valid_progeny, " . 
        "seq, marker, max_pct, ratio, ests, pe_radtags, blast_hits, external_id " .
        "FROM catalog_index " .
        "JOIN catalog_tags ON (catalog_index.cat_id=catalog_tags.id) " . 
        "LEFT JOIN catalog_annotations ON " . 
        "(" . 
        "catalog_index.batch_id=catalog_annotations.batch_id AND " . 
        "catalog_index.tag_id=catalog_annotations.catalog_id" . 
        ") " .
        "WHERE catalog_index.batch_id=?";

    $query .= apply_query_filters($filters);

    $sth->{'tag'} = $sth->{'dbh'}->prepare($query) or die($sth->{'dbh'}->errstr());
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

    my ($filter, $name, $value);

    while (@ARGV) {
        $_ = shift @ARGV;
        if    ($_ =~ /^-f$/) { $out_file = shift @ARGV; }
        elsif ($_ =~ /^-o$/) { $type     = shift @ARGV; }
        elsif ($_ =~ /^-b$/) { $batch_id = shift @ARGV; }
        elsif ($_ =~ /^-D$/) { $db       = shift @ARGV; }
        elsif ($_ =~ /^-d$/) { $debug++; }
	elsif ($_ =~ /^-v$/) { version(); exit(); }
        elsif ($_ =~ /^-h$/) { usage(); }
        elsif ($_ =~ /^-F$/) {
            $filter = shift @ARGV;
            ($name, $value) = split(/=/, $filter);

            if (length($name) == 0 || length($value) == 0) {
                print STDERR "Error parsing filter '$filter'\n";
                usage();
            }
            $filters{$name} = $value;

        } else {
            print STDERR "Unknown command line option: '$_'\n";
            usage();
        }
    }

    if ($out_file eq "") {
	print STDERR "You must specify the file to write data to!\n";
	usage();
    }

    if ($type ne "tsv" && $type ne "xls") {
        print STDERR "Unknown output file type specified '$type'.\n";
        usage();
    }
}

sub version {
    print STDERR "export_sql.pl ", stacks_version, "\n";
}

sub usage {
    version();

    print STDERR 
        "export_sql.pl -D db -b batch_id -f file -o tsv|xls [-F filter=value ...] [-d] [-h]\n", 
        "    D: database to export from.\n",
        "    b: batch ID of the dataset to export.\n",
        "    f: file to output data.\n",
        "    o: type of data to export: 'tsv' or 'xls'.\n",
        "    F: one or more filters in the format name=value.\n",
        "    h: display this help message.\n",
        "    d: turn on debug output.\n\n";
    exit(0);
}
