#!/usr/bin/perl
#
# Copyright 2010-2013, Julian Catchen <jcatchen@uoregon.edu>
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
use POSIX;
use constant stacks_version => "_VERSION_";

my $mysql_config = "_PKGDATADIR_" . "sql/mysql.cnf";
my $dry_run      = 0;
my $sql          = 1;
my $exe_path     = "_BINDIR_";
my $out_path     = "";
my $popmap_path  = "";
my $db           = "";
my $data_type    = "map";
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
my $alpha        = 0;
my $bound_low    = 0;
my $bound_high   = 1;

my @parents;
my @progeny;
my @samples;

my $cmd_str = $0 . " " . join(" ", @ARGV);

parse_command_line();

check_input_files(\@parents, \@progeny, \@samples);

my $cnf = (-e $ENV{"HOME"} . "/.my.cnf") ? $ENV{"HOME"} . "/.my.cnf" : $mysql_config;

#
# Check for the existence of the necessary pipeline programs
#
die ("Unable to find '" . $exe_path . "ustacks'.\n")     if (!-e $exe_path . "ustacks"     || !-x $exe_path . "ustacks");
die ("Unable to find '" . $exe_path . "cstacks'.\n")     if (!-e $exe_path . "cstacks"     || !-x $exe_path . "cstacks");
die ("Unable to find '" . $exe_path . "sstacks'.\n")     if (!-e $exe_path . "sstacks"     || !-x $exe_path . "sstacks");
die ("Unable to find '" . $exe_path . "genotypes'.\n")   if (!-e $exe_path . "genotypes"   || !-x $exe_path . "genotypes");
die ("Unable to find '" . $exe_path . "populations'.\n") if (!-e $exe_path . "populations" || !-x $exe_path . "populations");
die ("Unable to find '" . $exe_path . "index_radtags.pl'.\n") if (!-e $exe_path . "index_radtags.pl" || !-x $exe_path . "index_radtags.pl");

my ($i, $log, $log_fh, $pfile, $file, $num_files, $parent, $sample, %map);

$i         = 1;
$num_files = scalar(@parents) + scalar(@progeny) + scalar(@samples);

my (@types, $type, @pop_ids, $pop, %pops);

parse_population_map(\@samples, \@pop_ids, \%pops) if ($data_type eq "population");

foreach $parent (@parents) {
    push(@types, "parent");
    push(@pop_ids, 1);
}
foreach $parent (@progeny) {
    push(@types, "progeny");
    push(@pop_ids, 1);
}
foreach $parent (@samples) {
    push(@types, "sample");
}

my (@results, $minc, $minrc, $mind, $minsd, $rrep, $cmd, $cscale, $threads, $fuzzym, $dshapl, $ppath, $pop_cnt, 
    $bl, $bh, $al, $model);

$pop_cnt = scalar(keys %pops);
$minc    = $min_cov     > 0 ? "-m $min_cov"     : "";
$minrc   = $min_rcov    > 0 ? "-m $min_rcov"    : $minc;
$mind    = $min_dist    > 0 ? "-M $min_dist"    : "";
$minsd   = $min_sdist   > 0 ? "-N $min_sdist"   : "";
$dshapl  = $dis_shapl   > 0 ? "-H"              : "";
$cscale  = $cov_scale   > 0 ? "-S $cov_scale"   : "";
$threads = $num_threads > 0 ? "-p $num_threads" : "";
$fuzzym  = $fuzzy_match > 0 ? "-n $fuzzy_match" : "";
$ppath   = length($popmap_path) > 0 ? "-M $popmap_path" : "";
$bl      = $bound_low   > 0 ? "--bound_low $bound_low" : "";
$bh      = $bound_high  < 1 ? "--bound_high $bound_high" : "";
$al      = $alpha       > 0 ? "--alpha $alpha" : "";
$model   = (length($bl) > 0 || length($bh) > 0) ? "--model_type bounded" : "";

#
# Open the log file
#
$log = "$out_path/denovo_map.log";
open($log_fh, ">$log") or die("Unable to open log file '$log'; $!\n");

print $log_fh 
    "ref_map.pl started at ", strftime("%Y-%m-%d %H:%M:%S",(localtime(time))), "\n",
    $cmd_str, "\n";

if ($sql == 1 && $dry_run == 0) {
    #
    # SQL Batch ID for this set of Radtags, along with description and date of 
    # sequencing. Insert this batch data into the database.
    # 
    `mysql --defaults-file=$cnf $db -e "INSERT INTO batches SET id=$batch_id, description='$desc', date='$date', type='$data_type'"`;
} else {
    print STDERR  "mysql --defaults-file=$cnf $db -e \"INSERT INTO batches SET id=$batch_id, description='$desc', date='$date', type='$data_type'\"\n";
    print $log_fh "mysql --defaults-file=$cnf $db -e \"INSERT INTO batches SET id=$batch_id, description='$desc', date='$date', type='$data_type'\"\n";
}

foreach $sample (@parents, @progeny, @samples) {
    my ($ftype, $pfile) = "";
    my $gzip = 0;

    my ($prefix, $suffix) = ($sample =~ /^(.+)\.(.+)$/);

    if ($suffix eq "gz") {
	$gzip = 1;
	($prefix, $suffix) = ($prefix =~ /^(.+)\.(.+)$/);
    }

    if ($prefix =~ /^.*\/.+$/) {
        ($pfile) = ($prefix =~ /^.*\/(.+)$/);
    } else {
        $pfile = $prefix;
    }

    if ($gzip == 1) {
	if ($suffix =~ /^fa_?\d?$/ || $suffix =~ /^fasta_?\d?$/) {
	    $ftype = "gzfasta";
	} elsif ($suffix =~ /^fq_?\d?$/ || $suffix =~ /^fastq_?\d?$/) {
	    $ftype = "gzfastq";
	} else {
	    die("Unknown input file type.\n");
	}
    } else {
	if ($suffix =~ /^fa_?\d?$/ || $suffix =~ /^fasta_?\d?$/) {
	    $ftype = "fasta";
	} elsif ($suffix =~ /^fq_?\d?$/ || $suffix =~ /^fastq_?\d?$/) {
	    $ftype = "fastq";
	} else {
	    die("Unknown input file type.\n");
	}
    }

    $type = shift @types;
    $pop  = shift @pop_ids;

    printf("Identifying unique stacks; file % 3s of % 3s [%s]\n", $i, $num_files, $pfile);
    printf($log_fh "Identifying unique stacks; file % 3s of % 3s [%s]\n", $i, $num_files, $pfile);

    if ($sql == 1 && $dry_run == 0) {
	`mysql --defaults-file=$cnf $db -e "INSERT INTO samples SET sample_id=$i, batch_id=$batch_id, type='$type', file='$pfile', pop_id=$pop"`;
	@results = `mysql --defaults-file=$cnf $db -N -B -e "SELECT id FROM samples WHERE sample_id=$i AND batch_id=$batch_id AND type='$type' AND file='$pfile'"`;
	chomp $results[0];
	$sample_id = $results[0];
    } else {
	print STDERR "mysql --defaults-file=$cnf $db -e \"INSERT INTO samples SET sample_id=$i, batch_id=$batch_id, type='$type', file='$pfile', pop_id=$pop\"\n";
	print $log_fh "mysql --defaults-file=$cnf $db -e \"INSERT INTO samples SET sample_id=$i, batch_id=$batch_id, type='$type', file='$pfile', pop_id=$pop\"\n";
    }

    $map{$pfile} = $sample_id;

    if ($rep_tags) {
	$rrep = "-d -r";
    } else {
	$rrep = "";
    }

    if ($type eq "parent" || $type eq "sample") {
	$cmd = $exe_path . "ustacks -t $ftype -f $sample -o $out_path -i $sample_id $rrep $minc $mind $minsd $dshapl $cscale $model $bl $bh $al $threads 2>&1";
    } elsif ($type eq "progeny") {
	$cmd = $exe_path . "ustacks -t $ftype -f $sample -o $out_path -i $sample_id $rrep $minrc $mind $minsd $dshapl $cscale $model $bl $bh $al $threads 2>&1";
    }
    print STDERR "$cmd\n";
    print $log_fh    "$cmd\n";
    @results = `$cmd` if ($dry_run == 0);
    write_results(\@results, $log_fh);

    $file = "$out_path/$pfile" . ".tags.tsv";
    import_sql_file($log_fh, $file, "unique_tags", 0);

    $file = "$out_path/$pfile" . ".snps.tsv";
    import_sql_file($log_fh, $file, "snps", 0);

    $file = "$out_path/$pfile" . ".alleles.tsv";
    import_sql_file($log_fh, $file, "alleles", 0);

    $i++;

    $sample_id++ if ($sql == 0);
}

my ($rid, $pfile, $parents, $cat_file);

#
# Generate catalog of RAD-Tags
#
print STDERR "Generating RAD-Tag catalog...\n";
foreach $sample (@parents, @samples) {
    my ($prefix, $suffix) = ($sample =~ /^(.+)\.(.+)$/);

    if ($suffix eq "gz") {
	($prefix, $suffix) = ($prefix =~ /^(.+)\.(.+)$/);
    }

    if ($prefix =~ /^.*\/.+$/) {
        ($pfile) = ($prefix =~ /^.*\/(.+)$/);
    } else {
        $pfile = $prefix;
    }

    $parents .= "-s $out_path/$pfile ";
}

$cat_file = "batch_" . $batch_id;
$cmd      = $exe_path . "cstacks -b $batch_id -o $out_path $parents $threads $fuzzym 2>&1";
print STDERR  "$cmd\n";
print $log_fh "$cmd\n";
@results =    `$cmd` if ($dry_run == 0);
print $log_fh @results;

print STDERR "Importing catalog to MySQL database\n";
$file = "$out_path/$cat_file" . ".catalog.tags.tsv";
import_sql_file($log_fh, $file, "catalog_tags", 0);

$file = "$out_path/$cat_file" . ".catalog.snps.tsv";
import_sql_file($log_fh, $file, "catalog_snps", 0);

$file = "$out_path/$cat_file" . ".catalog.alleles.tsv";
import_sql_file($log_fh, $file, "catalog_alleles", 0);

#
# Match parents and progeny to the catalog
#
$i         = 1;
$num_files = scalar(@parents) + scalar(@progeny) + scalar(@samples);

foreach $sample (@parents, @progeny, @samples) {

    my ($prefix, $suffix) = ($sample =~ /^(.+)\.(.+)$/);

    if ($suffix eq "gz") {
	($prefix, $suffix) = ($prefix =~ /^(.+)\.(.+)$/);
    }

    if ($prefix =~ /^.*\/.+$/) {
        ($pfile) = ($prefix =~ /^.*\/(.+)$/);
    } else {
        $pfile = $prefix;
    }

    printf(STDERR "Matching RAD-Tags to catalog; file % 3s of % 3s [%s]\n", $i, $num_files, $pfile);

    $cmd = $exe_path . "sstacks -b $batch_id -c $out_path/$cat_file -s $out_path/$pfile -o $out_path $threads 2>&1";
    print STDERR  "$cmd\n";
    print $log_fh "$cmd\n";
    @results =    `$cmd` if ($dry_run == 0);
    print $log_fh @results;

    $file = "$out_path/" . $pfile . ".matches.tsv";
    import_sql_file($log_fh, $file, "matches", 0);

    $i++;
}

if ($data_type eq "map") {
    #
    # Generate a set of observed haplotypes and a set of markers and generic genotypes
    #
    $cmd = $exe_path . "genotypes -b $batch_id -P $out_path -t gen -r 1 -c -s 2>&1";
    print STDERR  "$cmd\n";
    print $log_fh "$cmd\n";
    @results =    `$cmd` if ($dry_run == 0);
    print $log_fh @results;

    $file = "$out_path/batch_" . $batch_id . ".markers.tsv";
    import_sql_file($log_fh, $file, "markers", 0);

    $file = "$out_path/batch_" . $batch_id . ".genotypes_1.txt";
    import_sql_file($log_fh, $file, "catalog_genotypes", 0);
} else {
    $cmd = $exe_path . "populations -b $batch_id -P $out_path -s $ppath 2>&1";
    print STDERR  "$cmd\n";
    print $log_fh "$cmd\n";
    @results =    `$cmd` if ($dry_run == 0);
    print $log_fh @results;

    $file = "$out_path/batch_" . $batch_id . ".markers.tsv";
    import_sql_file($log_fh, $file, "markers", 0);

    $file = "$out_path/batch_" . $batch_id . ".sumstats.tsv";
    import_sql_file($log_fh, $file, "sumstats", $pop_cnt+1);

    #
    # Import the Fst files.
    #
    my (@keys, $m, $n);
    @keys = sort keys %pops;
    for ($m = 0; $m < scalar(@keys); $m++) {
	for ($n = $m+1; $n < scalar(@keys); $n++) {
	    $file = "$out_path/batch_" . $batch_id . ".fst_" . $keys[$m] . "-" . $keys[$n] . ".tsv";
	    import_sql_file($log_fh, $file, "fst", 1);
	}
    }
}

if ($sql) {
    #
    # Index the radtags database
    #
    print STDERR "Indexing the database...\n";
    $cmd = $exe_path . "index_radtags.pl -D $db -t -c 2>&1";
    print STDERR  "$cmd\n";
    print $log_fh "$cmd\n";
    @results =    `$cmd` if ($dry_run == 0);
    print $log_fh @results;
}

print $log_fh "denovo_map.pl completed at ", strftime("%Y-%m-%d %H:%M:%S",(localtime(time))), "\n";

close($log_fh);

sub parse_population_map {
    my ($samples, $pop_ids, $pops) = @_;

    my ($fh, @parts, $line, %ids, $file, $path);

    if (length($popmap_path) == 0) {
	foreach $path (@{$samples}) {
	    push(@{$pop_ids}, 1);
	    $pops->{1}++;
	}
	return;
    }

    open($fh, "<$popmap_path") or die("Unable to open population map, '$popmap_path', $!\n");

    while ($line = <$fh>) {
	chomp $line;
	@parts = split(/\t/, $line);

	if (scalar(@parts) > 2) {
	    die("Unable to parse population map, '$popmap_path' (map should contain no more than two columns).\n");
	}

	if ($parts[1] !~ /\d+/) {
	    die("Unable to parse population map, '$popmap_path' (population ID in second column should be an integer).\n");
	}

	$ids{$parts[0]} = $parts[1];
    }

    foreach $path (@{$samples}) {
	my ($prefix, $suffix) = ($path =~ /^(.+)\.(.+)$/);

	if ($prefix =~ /^.*\/.+$/) {
	    ($file) = ($prefix =~ /^.*\/(.+)$/);
	} else {
	    $file = $prefix;
	}

	if (!defined($ids{$file})) {
	    die("Unable to find '$file' in the population map, '$popmap_path'.\n");
	}

	push(@{$pop_ids}, $ids{$file});
	$pops->{$ids{$file}}++;
    }

    print STDERR "Parsed population map: ", scalar(@{$samples}), " files in ", scalar(keys %{$pops}), " populations.\n";

    close($fh);
}

sub check_input_files {
    my ($parents, $progeny, $samples) = @_;

    #
    # Check that no duplicate files were specified.
    #
    my (%files, $file);
    foreach $file (@{$parents}, @{$progeny}, @{$samples}) {
	$files{$file}++;
    }
    foreach $file (keys %files) {
	if ($files{$file} > 1) {
	    print STDERR "A duplicate file was specified which may create undefined results, '$file'\n";
	    usage();
	}
    }

    #
    # Check that all the files exist and are accessible.
    #
    foreach $file (@{$parents}) {
	if (!-e $file) {
	    print STDERR "Unable to locate parental file '$file'\n";
	    usage();
	}
    }
    print STDERR "Found ", scalar(@{$parents}), " parental file(s).\n" if (scalar(@{$parents}) > 0);

    foreach $file (@{$progeny}) {
	if (!-e $file) {
	    print STDERR "Unable to locate progeny file '$file'\n";
	    usage();
	}
    }
    print STDERR "Found ", scalar(@{$progeny}), " progeny file(s).\n" if (scalar(@{$progeny}) > 0);

    foreach $file (@{$samples}) {
	if (!-e $file) {
	    print STDERR "Unable to locate sample file '$file'\n";
	    usage();
	}
    }
    print STDERR "Found ", scalar(@{$samples}), " sample file(s).\n" if (scalar(@{$samples}) > 0);
}

sub write_results {
    my ($results, $log_fh) = @_;

    my $line;

    foreach $line (@{$results}) {
	if ($line =~ /\r/) { 
	    $line =~ s/^.+\r(.*\n)$/\1/; 
	}
	print $log_fh $line;
    }
}

sub import_sql_file {
    my ($log_fh, $file, $table, $skip_lines) = @_;

    my (@results, $ignore);

    $ignore = "IGNORE $skip_lines LINES" if ($skip_lines > 0);

    @results = `mysql --defaults-file=$cnf $db -e "LOAD DATA LOCAL INFILE '$file' INTO TABLE $table $ignore"` if ($sql == 1 && $dry_run == 0);
    print STDERR  "mysql --defaults-file=$cnf $db -e \"LOAD DATA LOCAL INFILE '$file' INTO TABLE $table $ignore\"\n", @results;
    print $log_fh "mysql --defaults-file=$cnf $db -e \"LOAD DATA LOCAL INFILE '$file' INTO TABLE $table $ignore\"\n", @results;
}

sub parse_command_line {
    while (@ARGV) {
	$_ = shift @ARGV;
	if    ($_ =~ /^-p$/) { push(@parents, shift @ARGV); }
	elsif ($_ =~ /^-r$/) { push(@progeny, shift @ARGV); }
	elsif ($_ =~ /^-s$/) { push(@samples, shift @ARGV); }
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
	elsif ($_ =~ /^-i$/) { $sample_id   = shift @ARGV; }
	elsif ($_ =~ /^-a$/) { $date        = shift @ARGV; }
	elsif ($_ =~ /^-S$/) { $sql         = 0; }
	elsif ($_ =~ /^-B$/) { $db          = shift @ARGV; }
	elsif ($_ =~ /^-O$/) { $popmap_path = shift @ARGV; }
       	elsif ($_ =~ /^-H$/) { $dis_shapl++; }
	elsif ($_ =~ /^--bound_low$/)  { $bound_low  = shift @ARGV; }
	elsif ($_ =~ /^--bound_high$/) { $bound_high = shift @ARGV; }
	elsif ($_ =~ /^--alpha$/)      { $alpha      = shift @ARGV; }
	elsif ($_ =~ /^-d$/) { $dry_run++; }
	elsif ($_ =~ /^-v$/) { version(); exit(); }
	elsif ($_ =~ /^-h$/) { usage(); }
	else {
	    print STDERR "Unknown command line option: '$_'\n";
	    usage();
	}
    }

    $exe_path = $exe_path . "/"          if (substr($out_path, -1) ne "/");
    $out_path = substr($out_path, 0, -1) if (substr($out_path, -1) eq "/");

    if ($batch_id !~ /^\d+$/ || $batch_id < 0) {
	print STDERR "You must specify a batch ID and it must be an integer (e.g. 1, 2, 3).\n";
	usage();
    }

    if ($sql > 0 && length($date) == 0) {
	$date = strftime("%Y-%m-%d", (localtime(time)));
    }

    if (scalar(@parents) > 0 && scalar(@samples) > 0) {
	print STDERR "You must specify either parent or sample files, but not both.\n";
	usage();
    }

    if (scalar(@parents) == 0 && scalar(@samples) == 0) {
	print STDERR "You must specify at least one parent or sample file.\n";
	usage();
    }

    if (scalar(@samples) > 0) {
	$data_type = "population";
    } else {
	$data_type = "map";
    }

    if ($alpha != 0 && $alpha != 0.1 && $alpha != 0.05 && $alpha != 0.01 && $alpha != 0.001) {
	print STDERR "The value of alpha must be 0.1, 0.05, 0.01, or 0.001.\n";
	usage();
    }

    if ($bound_low != 0 || $bound_high != 1) {
	if ($bound_low < 0 || $bound_low > 1 || $bound_high < 0 || $bound_high > 1) {
	    print STDERR "SNP model bounds must be between 0 and 1.0.\n";
	    usage();
	}
    }
}

sub version {
    print STDERR "denovo_map.pl ", stacks_version, "\n";
}

sub usage {
    version();

    print STDERR <<EOQ; 
denovo_map.pl -p path -r path [-s path] -o path [-t] [-m min_cov] [-M mismatches] [-n mismatches] [-T num_threads] [-O popmap] [-B db -b batch_id -D desc -a yyyy-mm-dd] [-S -i num] [-e path] [-d] [-h]
    p: path to a FASTQ/FASTA file containing parent sequences from a mapping cross.
    r: path to a FASTQ/FASTA file containing progeny sequences from a mapping cross.
    s: path to a FASTQ/FASTA file contiaining an individual sample from a population.
    o: path to write pipeline output files.
    O: if analyzing one or more populations, specify a pOpulation map.
    T: specify the number of threads to execute.
    e: executable path, location of pipeline programs.
    d: perform a dry run. Do not actually execute any programs, just print what would be executed.
    h: display this help message.

  Stack assembly options:
    m: specify a minimum number of identical, raw reads required to create a stack.
    P: specify a minimum number of identical, raw reads required to create a stack in 'progeny' individuals.
    M: specify the number of mismatches allowed between loci when processing a single individual (default 2).
    N: specify the number of mismatches allowed when aligning secondary reads to primary stacks (default M+2).
    n: specify the number of mismatches allowed between loci when building the catalog (default 0).
    t: remove, or break up, highly repetitive RAD-Tags in the ustacks program.
    H: disable calling haplotypes from secondary reads.

  Database options:
    b: batch ID representing this dataset.
    B: specify a database to load data into.
    D: batch description
    a: batch run date, yyyy-mm-dd
    S: disable recording SQL data in the database.
    i: starting sample_id, this is determined automatically if database interaction is enabled.

  SNP Model Options (these options are passed on to ustacks):
    --bound_low: lower bound for epsilon, the error rate, between 0 and 1.0 (default 0).
    --bound_high: upper bound for epsilon, the error rate, between 0 and 1.0 (default 1).
    --alpha: chi square significance level required to call a heterozygote or homozygote, either 0.1, 0.05 (default), 0.01, or 0.001.

EOQ

exit(0);
}
