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
# Process the data for a genetic map: build stacks in parents and progeny, 
# create a catalog from the parents, and match progeny against the catatlog.
# Call genotypes, and load all data into an MySQL database along the way.
#
# For the database interactions to work, the 'mysql' program is expected to be
# on the path and sufficient permissions set to access the specified database.
#
# By Julian Catchen <jcatchen@uoregon.edu>
#

use strict;
use constant stacks_version => "_VERSION_";

my $debug        = 0;
my $mysql_config = "_PKGDATADIR_" . "sql/mysql.cnf";
my $sql          = 1;
my $exe_path     = "_BINDIR_";
my $out_path     = "";
my $db           = "";
my $rep_tags     = 0;
my $min_cov      = 0;
my $min_rcov     = 0;
my $min_dist     = 0;
my $min_sdist    = 0;
my $dis_shapl    = 0;
my $fuzzy_match  = 0;
my $num_threads  = 0;
my $cov_scale    = 0;
my $batch_id     = -1;
my $sample_id    = 1;
my $desc         = ""; # Database description of this dataset
my $date         = ""; # Date relevent to this data, formatted for SQL: 2009-05-31

my @parents;
my @progeny;

parse_command_line();

my $cnf = (-e $ENV{"HOME"} . "/.my.cnf") ? $ENV{"HOME"} . "/.my.cnf" : $mysql_config;

#
# Check for the existence of the necessary pipeline programs
#
die ("Unable to find '" . $exe_path . "ustacks'.\n")   if (!-e $exe_path . "ustacks"   || !-x $exe_path . "ustacks");
die ("Unable to find '" . $exe_path . "cstacks'.\n")   if (!-e $exe_path . "cstacks"   || !-x $exe_path . "cstacks");
die ("Unable to find '" . $exe_path . "sstacks'.\n")   if (!-e $exe_path . "sstacks"   || !-x $exe_path . "sstacks");
die ("Unable to find '" . $exe_path . "genotypes'.\n") if (!-e $exe_path . "genotypes" || !-x $exe_path . "genotypes");
die ("Unable to find '" . $exe_path . "index_radtags.pl'.\n") if (!-e $exe_path . "index_radtags.pl" || !-x $exe_path . "index_radtags.pl");

my ($i, $log, $log_fh, $pfile, $file, $num_files, $parent, $sample, %map);

$i         = 1;
$num_files = scalar(@parents) + scalar(@progeny);

if ($sql) {
    #
    # SQL Batch ID for this set of Radtags, along with description and date of 
    # sequencing. Insert this batch data into the database.
    # 
    `mysql --defaults-file=$cnf $db -e "INSERT INTO batches SET id=$batch_id, description='$desc', date='$date'"`;
}

my (@types, $type);
foreach $parent (@parents) {
    push(@types, "parent");
}
foreach $parent (@progeny) {
    push(@types, "progeny");
}

my (@results, $minc, $minrc, $mind, $minsd, $rrep, $cmd, $cscale, $threads, $fuzzym, $dshapl);

$minc    = $min_cov     > 0 ? "-m $min_cov"     : "";
$minrc   = $min_rcov    > 0 ? "-m $min_rcov"    : $minc;
$mind    = $min_dist    > 0 ? "-M $min_dist"    : "";
$minsd   = $min_sdist   > 0 ? "-N $min_sdist"   : "";
$dshapl  = $dis_shapl   > 0 ? "-H"              : "";
$cscale  = $cov_scale   > 0 ? "-S $cov_scale"   : "";
$threads = $num_threads > 0 ? "-p $num_threads" : "";
$fuzzym  = $fuzzy_match > 0 ? "-n $fuzzy_match" : "";

#
# Open the log file
#
$log = "$out_path/denovo_map.log";
open($log_fh, ">$log") or die("Unable to open log file '$log'; $!\n");

foreach $sample (@parents, @progeny) {
    my ($ftype, $pfile) = "";

    my ($prefix, $suffix) = ($sample =~ /^(.+)\.(.+)$/);

    if ($prefix =~ /^.*\/.+$/) {
        ($pfile) = ($prefix =~ /^.*\/(.+)$/);
    } else {
        $pfile = $prefix;
    }

    if ($suffix =~ /^fa_?\d?$/ || $suffix =~ /^fasta_?\d?$/) {
        $ftype = "fasta";
    } elsif ($suffix =~ /^fq_?\d?$/ || $suffix =~ /^fastq_?\d?$/) {
        $ftype = "fastq";
    } else {
        die("Unknown input file type.\n");
    }

    $type = shift @types;

    printf("Identifying unique stacks; file % 3s of % 3s [%s]\n", $i, $num_files, $pfile);
    printf($log_fh "Identifying unique stacks; file % 3s of % 3s [%s]\n", $i, $num_files, $pfile);

    if ($sql) {
	`mysql --defaults-file=$cnf $db -e "INSERT INTO samples SET sample_id=$i, batch_id=$batch_id, type='$type', file='$pfile'"`;
	@results = `mysql --defaults-file=$cnf $db -N -B -e "SELECT id FROM samples WHERE sample_id=$i AND batch_id=$batch_id AND type='$type' AND file='$pfile'"`;
	chomp $results[0];
	$sample_id = $results[0];
    } else {
	print STDERR "mysql --defaults-file=$cnf $db -e \"INSERT INTO samples SET sample_id=$i, batch_id=$batch_id, type='$type', file='$pfile'\"\n";
    }

    $map{$pfile} = $sample_id;

    #
    # Calculate the expected coverage by dividing the number of reads for this sample
    # by the number of RAD sites in the genome.
    #
    if ($rep_tags) {
	$rrep = "-d -r";
    } else {
	$rrep = "";
    }

    if ($type eq "parent") {
	$cmd = $exe_path . "ustacks -t $ftype -f $sample -o $out_path -i $sample_id $rrep $minc $mind $minsd $dshapl $cscale $threads 2>&1";
    } elsif ($type eq "progeny") {
	$cmd = $exe_path . "ustacks -t $ftype -f $sample -o $out_path -i $sample_id $rrep $minrc $mind $minsd $dshapl $cscale $threads 2>&1";
    }
    print STDERR "$cmd\n";
    print $log_fh    "$cmd\n";
    @results = `$cmd`;
    print $log_fh @results;

    $file = "$out_path/$pfile" . ".tags.tsv";
    import_sql_file($file, "unique_tags");

    $file = "$out_path/$pfile" . ".snps.tsv";
    import_sql_file($file, "snps");

    $file = "$out_path/$pfile" . ".alleles.tsv";
    import_sql_file($file, "alleles");

    $i++;

    $sample_id++ if ($sql == 0);
}

my ($rid, $pfile, $parents, $cat_file);

#
# Generate catalog of RAD-Tags
#
print STDERR "Generating RAD-Tag catalog...\n";
foreach $sample (@parents) {
    my ($prefix, $suffix) = ($sample =~ /^(.+)\.(.+)$/);

    if ($prefix =~ /^.*\/.+$/) {
        ($pfile) = ($prefix =~ /^.*\/(.+)$/);
    } else {
        $pfile = $prefix;
    }

    $parents .= "-s $out_path/$pfile -S " . $map{$pfile} . " ";
}

$cat_file = "batch_" . $batch_id;
$cmd      = $exe_path . "cstacks -b $batch_id -o $out_path $parents $threads $fuzzym 2>&1";
print STDERR  "$cmd\n";
print $log_fh "$cmd\n";
@results =    `$cmd`;
print $log_fh @results;

print STDERR "Importing catalog to MySQL database\n";
$file = "$out_path/$cat_file" . ".catalog.tags.tsv";
import_sql_file($file, "catalog_tags");

$file = "$out_path/$cat_file" . ".catalog.snps.tsv";
import_sql_file($file, "catalog_snps");

$file = "$out_path/$cat_file" . ".catalog.alleles.tsv";
import_sql_file($file, "catalog_alleles");

#
# Match parents and progeny to the catalog
#
$i         = 1;
$num_files = scalar(@parents) + scalar(@progeny);

foreach $sample (@parents, @progeny) {

    my ($prefix, $suffix) = ($sample =~ /^(.+)\.(.+)$/);

    if ($prefix =~ /^.*\/.+$/) {
        ($pfile) = ($prefix =~ /^.*\/(.+)$/);
    } else {
        $pfile = $prefix;
    }

    printf(STDERR "Matching RAD-Tags to catalog; file % 3s of % 3s [%s]\n", $i, $num_files, $pfile);

    $rid = $map{$pfile};

    $cmd = $exe_path . "sstacks -b $batch_id -c $out_path/$cat_file -s $out_path/$pfile -S $rid -o $out_path $threads 2>&1";
    print STDERR  "$cmd\n";
    print $log_fh "$cmd\n";
    @results =    `$cmd`;
    print $log_fh @results;

    $file = "$out_path/" . $pfile . ".matches.tsv";
    import_sql_file($file, "matches");

    $i++;
}

#
# Generate a set of observed haplotypes and a set of markers and generic genotypes
#
$cmd = $exe_path . "genotypes -b $batch_id -P $out_path -t gen -r 1 -c -s  2>&1";
print STDERR  "$cmd\n";
print $log_fh "$cmd\n";
@results =    `$cmd`;
print $log_fh @results;

$file = "$out_path/batch_" . $batch_id . ".markers.tsv";
import_sql_file($file, "markers");

$file = "$out_path/batch_" . $batch_id . ".genotypes_1.txt";
import_sql_file($file, "catalog_genotypes");

if ($sql) {
    #
    # Index the radtags database
    #
    $cmd = $exe_path . "index_radtags.pl -D $db -t -c 2>&1";
    print STDERR  "$cmd\n";
    print $log_fh "$cmd\n";
    @results =    `$cmd`;
    print $log_fh @results;
}

close($log_fh);

sub import_sql_file {
    my ($file, $table) = @_;

    my (@results);

    @results = `mysql --defaults-file=$cnf $db -e "LOAD DATA LOCAL INFILE '$file' INTO TABLE $table"` if ($sql);
    print STDERR "mysql --defaults-file=$cnf $db -e \"LOAD DATA LOCAL INFILE '$file' INTO TABLE $table\"\n", @results, "\n";
}

sub parse_command_line {
    while (@ARGV) {
	$_ = shift @ARGV;
	if    ($_ =~ /^-p$/) { push(@parents, shift @ARGV); }
	elsif ($_ =~ /^-r$/) { push(@progeny, shift @ARGV); }
       	elsif ($_ =~ /^-t$/) { $rep_tags++; }
	elsif ($_ =~ /^-o$/) { $out_path    = shift @ARGV; }
	elsif ($_ =~ /^-m$/) { $min_cov     = shift @ARGV; }
	elsif ($_ =~ /^-P$/) { $min_rcov    = shift @ARGV; }
	elsif ($_ =~ /^-M$/) { $min_dist    = shift @ARGV; }
	elsif ($_ =~ /^-N$/) { $min_sdist   = shift @ARGV; }
        elsif ($_ =~ /^-n$/) { $fuzzy_match = shift @ARGV; }
        elsif ($_ =~ /^-T$/) { $num_threads = shift @ARGV; }
	elsif ($_ =~ /^-c$/) { $cov_scale   = shift @ARGV; }
	elsif ($_ =~ /^-D$/) { $desc        = shift @ARGV; }
	elsif ($_ =~ /^-e$/) { $exe_path    = shift @ARGV; }
	elsif ($_ =~ /^-b$/) { $batch_id    = shift @ARGV; }
	elsif ($_ =~ /^-s$/) { $sample_id   = shift @ARGV; }
	elsif ($_ =~ /^-a$/) { $date        = shift @ARGV; }
	elsif ($_ =~ /^-S$/) { $sql         = 0; }
	elsif ($_ =~ /^-B$/) { $db          = shift @ARGV; }
       	elsif ($_ =~ /^-H$/) { $dis_shapl++; }
	elsif ($_ =~ /^-d$/) { $debug++; }
	elsif ($_ =~ /^-v$/) { version(); exit(); }
	elsif ($_ =~ /^-h$/) { usage(); }
	else {
	    print STDERR "Unknown command line option: '$_'\n";
	    usage();
	}
    }

    $exe_path = $exe_path . "/"          if (substr($out_path, -1) ne "/");
    $out_path = substr($out_path, 0, -1) if (substr($out_path, -1) eq "/");

    if ($batch_id < 0) {
	print STDERR "You must specify a batch ID.\n";
	usage();
    }

    if (scalar(@parents) == 0) {
	print STDERR "You must specify at least one parent file.\n";
	usage();
    }
}

sub version {
    print STDERR "denovo_map.pl ", stacks_version, "\n";
}

sub usage {
    version();

    print STDERR <<EOQ; 
denovo_map.pl -p path -r path -o path [-t] [-m min_cov] [-M mismatches] [-n mismatches] [-T num_threads] [-b batch_id -D desc -a yyyy-mm-dd] [-S -s num] [-e path] [-d] [-h]
    p: path to a FASTQ/FASTA file containing parent sequences.
    r: path to a FASTQ/FASTA file containing progeny sequences.
    o: path to write pipeline output files.
    m: specify a minimum number of identical, raw reads required to create a stack.
    P: specify a minimum number of identical, raw reads required to create a stack in 'progeny' individuals.
    M: specify the number of mismatches allowed between loci when processing a single individual (default 2).
    N: specify the number of mismatches allowed when aligning secondary reads to primary stacks (default M+2).
    n: specify the number of mismatches allowed between loci when building the catalog (default 0).
    t: remove, or break up, highly repetitive RAD-Tags in the ustacks program.
    H: disable calling haplotypes from secondary reads.
    T: specify the number of threads to execute.
    B: specify a database to load data into.
    b: batch ID representing this dataset.
    D: batch description
    a: batch run date, yyyy-mm-dd
    S: disable recording SQL data in the database.
    s: starting sample_id, this is determined automatically if database interaction is enabled.
    e: executable path, location of pipeline programs.
    h: display this help message.
    d: turn on debug output.

EOQ

exit(0);
}
