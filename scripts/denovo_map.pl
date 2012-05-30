#!/usr/bin/perl
#
# Copyright 2010-2012, Julian Catchen <jcatchen@uoregon.edu>
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

if ($sql) {
    #
    # SQL Batch ID for this set of Radtags, along with description and date of 
    # sequencing. Insert this batch data into the database.
    # 
    `mysql --defaults-file=$cnf $db -e "INSERT INTO batches SET id=$batch_id, description='$desc', date='$date', type='$data_type'"`;
}

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

my (@results, $minc, $minrc, $mind, $minsd, $rrep, $cmd, $cscale, $threads, $fuzzym, $dshapl, $ppath, $pop_cnt);

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

#
# Open the log file
#
$log = "$out_path/denovo_map.log";
open($log_fh, ">$log") or die("Unable to open log file '$log'; $!\n");

print $log_fh 
    "ref_map.pl started at ", strftime("%Y-%m-%d %H:%M:%S",(localtime(time))), "\n",
    $cmd_str, "\n";

foreach $sample (@parents, @progeny, @samples) {
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
    $pop  = shift @pop_ids;

    printf("Identifying unique stacks; file % 3s of % 3s [%s]\n", $i, $num_files, $pfile);
    printf($log_fh "Identifying unique stacks; file % 3s of % 3s [%s]\n", $i, $num_files, $pfile);

    if ($sql) {
	`mysql --defaults-file=$cnf $db -e "INSERT INTO samples SET sample_id=$i, batch_id=$batch_id, type='$type', file='$pfile', pop_id=$pop"`;
	@results = `mysql --defaults-file=$cnf $db -N -B -e "SELECT id FROM samples WHERE sample_id=$i AND batch_id=$batch_id AND type='$type' AND file='$pfile'"`;
	chomp $results[0];
	$sample_id = $results[0];
    } else {
	print STDERR "mysql --defaults-file=$cnf $db -e \"INSERT INTO samples SET sample_id=$i, batch_id=$batch_id, type='$type', file='$pfile', pop_id=$pop\"\n";
	print $log_fh "mysql --defaults-file=$cnf $db -e \"INSERT INTO samples SET sample_id=$i, batch_id=$batch_id, type='$type', file='$pfile', pop_id=$pop\"\n";
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

    if ($type eq "parent" || $type eq "sample") {
	$cmd = $exe_path . "ustacks -t $ftype -f $sample -o $out_path -i $sample_id $rrep $minc $mind $minsd $dshapl $cscale $threads 2>&1";
    } elsif ($type eq "progeny") {
	$cmd = $exe_path . "ustacks -t $ftype -f $sample -o $out_path -i $sample_id $rrep $minrc $mind $minsd $dshapl $cscale $threads 2>&1";
    }
    print STDERR "$cmd\n";
    print $log_fh    "$cmd\n";
    @results = `$cmd`;
    print $log_fh @results;

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
    import_sql_file($log_fh, $file, "matches", 0);

    $i++;
}

if ($data_type eq "map") {
    #
    # Generate a set of observed haplotypes and a set of markers and generic genotypes
    #
    $cmd = $exe_path . "genotypes -b $batch_id -P $out_path -t gen -r 1 -c -s  2>&1";
    print STDERR  "$cmd\n";
    print $log_fh "$cmd\n";
    @results =    `$cmd`;
    print $log_fh @results;

    $file = "$out_path/batch_" . $batch_id . ".markers.tsv";
    import_sql_file($log_fh, $file, "markers", 0);

    $file = "$out_path/batch_" . $batch_id . ".genotypes_1.txt";
    import_sql_file($log_fh, $file, "catalog_genotypes", 0);
} else {
    $cmd = $exe_path . "populations -b $batch_id -P $out_path -r 1 -s $ppath 2>&1";
    print STDERR  "$cmd\n";
    print $log_fh "$cmd\n";
    @results =    `$cmd`;
    print $log_fh @results;

    $file = "$out_path/batch_" . $batch_id . ".markers.tsv";
    import_sql_file($log_fh, $file, "markers", 0);

    $file = "$out_path/batch_" . $batch_id . ".sumstats.tsv";
    import_sql_file($log_fh, $file, "sumstats", $pop_cnt+1);

    #
    # Import the Fst files.
    #
    my ($m, $n);
    foreach $m (sort keys %pops) {
	foreach $n (sort keys %pops) {
	    $file = "$out_path/batch_" . $batch_id . ".fst_" . $m . "-" . $n . ".tsv";
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
    @results =    `$cmd`;
    print $log_fh @results;
}

print $log_fh "denovo_map.pl completed at ", strftime("%Y-%m-%d %H:%M:%S",(localtime(time))), "\n";

close($log_fh);

sub parse_population_map {
    my ($samples, $pop_ids, $pops) = @_;

    my ($fh, @parts, $line, %ids, $file, $path);

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

sub import_sql_file {
    my ($log_fh, $file, $table, $skip_lines) = @_;

    my (@results, $ignore);

    $ignore = "IGNORE $skip_lines LINES" if ($skip_lines > 0);

    @results = `mysql --defaults-file=$cnf $db -e "LOAD DATA LOCAL INFILE '$file' INTO TABLE $table $ignore"` if ($sql);
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
}

sub version {
    print STDERR "denovo_map.pl ", stacks_version, "\n";
}

sub usage {
    version();

    print STDERR <<EOQ; 
denovo_map.pl -p path -r path [-s path] -o path [-t] [-m min_cov] [-M mismatches] [-n mismatches] [-T num_threads] [-O popmap] [-B db -b batch_id -D desc -a yyyy-mm-dd] [-S -i num] [-e path] [-h]
    p: path to a FASTQ/FASTA file containing parent sequences from a mapping cross.
    r: path to a FASTQ/FASTA file containing progeny sequences from a mapping cross.
    s: path to a Bowtie/SAM file contiaining an individual sample from a population.
    o: path to write pipeline output files.
    m: specify a minimum number of identical, raw reads required to create a stack.
    P: specify a minimum number of identical, raw reads required to create a stack in 'progeny' individuals.
    M: specify the number of mismatches allowed between loci when processing a single individual (default 2).
    N: specify the number of mismatches allowed when aligning secondary reads to primary stacks (default M+2).
    n: specify the number of mismatches allowed between loci when building the catalog (default 0).
    t: remove, or break up, highly repetitive RAD-Tags in the ustacks program.
    H: disable calling haplotypes from secondary reads.
    T: specify the number of threads to execute.
    O: if analyzing one or more populations, specify a pOpulation map
    B: specify a database to load data into.
    b: batch ID representing this dataset.
    D: batch description
    a: batch run date, yyyy-mm-dd
    S: disable recording SQL data in the database.
    i: starting sample_id, this is determined automatically if database interaction is enabled.
    e: executable path, location of pipeline programs.
    h: display this help message.

EOQ

exit(0);
}
