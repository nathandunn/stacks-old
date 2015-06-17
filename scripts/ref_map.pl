#!/usr/bin/env perl
#
# Copyright 2010-2015, Julian Catchen <jcatchen@illinois.edu>
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

use strict;
use POSIX;
use File::Temp qw/ mktemp /;
use File::Spec;
use constant stacks_version => "_VERSION_";

my $dry_run      = 0;
my $sql          = 1;
my $mysql_config = "_PKGDATADIR_" . "sql/mysql.cnf";
my $exe_path     = "_BINDIR_";
my $out_path     = "";
my $popmap_path  = "";
my $db           = "";
my $data_type    = "map";
my $batch_id     = -1;
my $sample_id    = 1;
my $desc         = ""; # Database description of this dataset
my $date         = ""; # Date relevent to this data, formatted for SQL: 2009-05-31

my @parents;
my @progeny;
my @samples;

my (@_pstacks, @_cstacks, @_sstacks, @_genotypes, @_populations);

my $cmd_str = $0 . " " . join(" ", @ARGV);

parse_command_line();

check_input_files(\@parents, \@progeny, \@samples);

my $cnf = (-e $ENV{"HOME"} . "/.my.cnf") ? $ENV{"HOME"} . "/.my.cnf" : $mysql_config;

#
# Check for the existence of the necessary pipeline programs
#
die ("Unable to find '" . $exe_path . "pstacks'.\n")     if (!-e $exe_path . "pstacks"     || !-x $exe_path . "pstacks");
die ("Unable to find '" . $exe_path . "cstacks'.\n")     if (!-e $exe_path . "cstacks"     || !-x $exe_path . "cstacks");
die ("Unable to find '" . $exe_path . "sstacks'.\n")     if (!-e $exe_path . "sstacks"     || !-x $exe_path . "sstacks");
die ("Unable to find '" . $exe_path . "genotypes'.\n")   if (!-e $exe_path . "genotypes"   || !-x $exe_path . "genotypes");
die ("Unable to find '" . $exe_path . "populations'.\n") if (!-e $exe_path . "populations" || !-x $exe_path . "populations");
die ("Unable to find '" . $exe_path . "index_radtags.pl'.\n") if (!-e $exe_path . "index_radtags.pl" || !-x $exe_path . "index_radtags.pl");

my ($i, $log, $log_fh, $pipe_fh, $pfile, $file, $num_files, $parent, $sample, %map);

$i         = 1;
$num_files = scalar(@parents) + scalar(@progeny) + scalar(@samples);

my (@types, $type, @pop_ids, $pop, %pops, @grp_ids, $grp, %grps);

parse_population_map(\@samples, \@pop_ids, \%pops, \@grp_ids, \%grps) if ($data_type eq "population");

foreach $parent (@parents) {
    push(@types, "parent");
    push(@pop_ids, "1");
    push(@grp_ids, "1");
}
foreach $parent (@progeny) {
    push(@types, "progeny");
    push(@pop_ids, "1");
    push(@grp_ids, "1");
}
foreach $parent (@samples) {
    push(@types, "sample");
}

my (@results, $cmd, $pop_cnt);

$pop_cnt = scalar(keys %pops);

#
# Open the log file
#
$log = "$out_path/ref_map.log";
open($log_fh, ">$log") or die("Unable to open log file '$log'; $!\n");

print $log_fh 
    "ref_map.pl version ", stacks_version, " started at ", strftime("%Y-%m-%d %H:%M:%S",(localtime(time))), "\n",
    $cmd_str, "\n";

if ($sql == 1) {
    #
    # SQL Batch ID for this set of Radtags, along with description and date of 
    # sequencing. Insert this batch data into the database.
    # 
    `mysql --defaults-file=$cnf $db -e "INSERT INTO batches SET id=$batch_id, description='$desc', date='$date', type='$data_type'"` if ($dry_run == 0);

    print $log_fh "mysql --defaults-file=$cnf $db -e \"INSERT INTO batches SET id=$batch_id, description='$desc', date='$date', type='$data_type'\"\n";
}

my $gzip = 0;

foreach $sample (@parents, @progeny, @samples) {
    my ($ftype, $pfile) = "";

    my ($prefix, $suffix) = ($sample =~ /^(.+)\.(.+)$/);

    if ($prefix =~ /^.*\/.+$/) {
        ($pfile) = ($prefix =~ /^.*\/(.+)$/);
    } else {
        $pfile = $prefix;
    }

    if ($suffix =~ /^bowtie$/) {
        $ftype = "bowtie";
    } elsif ($suffix =~ /^sam$/) {
        $ftype = "sam";
    } elsif ($suffix =~ /^bam$/) {
        $ftype = "bam";
	$gzip  = 1;
    } elsif ($suffix =~ /^map$/) {
        $ftype = "tsv";
    } else {
        die("Unknown input file type.\n");
    }

    $type = shift @types;
    $pop  = shift @pop_ids;
    $grp  = shift @grp_ids;

    printf("Identifying unique stacks; file % 3s of % 3s [%s]\n", $i, $num_files, $pfile);
    printf($log_fh "Identifying unique stacks; file % 3s of % 3s [%s]\n", $i, $num_files, $pfile);

    if ($sql == 1) {
	if ($dry_run == 0) {
	    `mysql --defaults-file=$cnf $db -e "INSERT INTO samples SET sample_id=$i, batch_id=$batch_id, type='$type', file='$pfile', pop_id='$pop', group_id='$grp'"`;
	    @results = `mysql --defaults-file=$cnf $db -N -B -e "SELECT id FROM samples WHERE sample_id=$i AND batch_id=$batch_id AND type='$type' AND file='$pfile'"`;
	    chomp $results[0];
	    $sample_id = $results[0];
	}
	print $log_fh "mysql --defaults-file=$cnf $db -e \"INSERT INTO samples SET sample_id=$i, batch_id=$batch_id, type='$type', file='$pfile', pop_id='$pop', group_id='$grp'\"\n";
    }

    $map{$pfile} = $sample_id;

    $cmd = $exe_path . "pstacks -t $ftype -f $sample -o $out_path -i $sample_id " . join(" ", @_pstacks) . " 2>&1";
    print STDERR  "  $cmd\n";
    print $log_fh "$cmd\n";
    @results = `$cmd` if ($dry_run == 0);
    write_results(\@results, $log_fh);

    print STDERR "  Loading pstacks output to $db..." if ($sql == 1);

    if ($gzip == 1) {
	$file = "$out_path/$pfile" . ".tags.tsv.gz";
	import_gzsql_file($log_fh, $file, "unique_tags", 1);

	$file = "$out_path/$pfile" . ".snps.tsv.gz";
	import_gzsql_file($log_fh, $file, "snps", 1);

	$file = "$out_path/$pfile" . ".alleles.tsv.gz";
	import_gzsql_file($log_fh, $file, "alleles", 1);

    } else {
	$file = "$out_path/$pfile" . ".tags.tsv";
	import_sql_file($log_fh, $file, "unique_tags", 1);

	$file = "$out_path/$pfile" . ".snps.tsv";
	import_sql_file($log_fh, $file, "snps", 1);

	$file = "$out_path/$pfile" . ".alleles.tsv";
	import_sql_file($log_fh, $file, "alleles", 1);
    }
    print STDERR "done.\n" if ($sql == 1);

    $i++;

    $sample_id++ if ($sql == 0);
}

my ($rid, $pfile, $parents, $cat_file);

#
# Generate catalog of RAD-Tags
#
print STDERR "Generating catalog...\n";
foreach $sample (@parents, @samples) {
    my ($prefix, $suffix) = ($sample =~ /^(.+)\.(.+)$/);

    if ($prefix =~ /^.*\/.+$/) {
        ($pfile) = ($prefix =~ /^.*\/(.+)$/);
    } else {
        $pfile = $prefix;
    }

    $parents .= "-s $out_path/$pfile ";
}

$cat_file = "batch_" . $batch_id;
$cmd      = $exe_path . "cstacks -g -b $batch_id -o $out_path $parents " . join(" ", @_cstacks) . " 2>&1";
print STDERR  "  $cmd\n";
print $log_fh "$cmd\n";

if ($dry_run == 0) {
    open($pipe_fh, "$cmd |");
    while (<$pipe_fh>) {
	print $log_fh $_;
	if ($_ =~ /failed/i) { print STDERR "Catalog construction failed.\n"; exit(1); }
    }
    close($pipe_fh);
}

print STDERR "  Importing catalog to MySQL database..." if ($sql == 1);

if ($gzip == 1) {
    $file = "$out_path/$cat_file" . ".catalog.tags.tsv.gz";
    import_gzsql_file($log_fh, $file, "catalog_tags", 1);

    $file = "$out_path/$cat_file" . ".catalog.snps.tsv.gz";
    import_gzsql_file($log_fh, $file, "catalog_snps", 1);

    $file = "$out_path/$cat_file" . ".catalog.alleles.tsv.gz";
    import_gzsql_file($log_fh, $file, "catalog_alleles", 1);

} else {
    $file = "$out_path/$cat_file" . ".catalog.tags.tsv";
    import_sql_file($log_fh, $file, "catalog_tags", 1);

    $file = "$out_path/$cat_file" . ".catalog.snps.tsv";
    import_sql_file($log_fh, $file, "catalog_snps", 1);

    $file = "$out_path/$cat_file" . ".catalog.alleles.tsv";
    import_sql_file($log_fh, $file, "catalog_alleles", 1);
}
print STDERR "done.\n" if ($sql == 1);

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

    printf(STDERR "Matching samples to catalog; file % 3s of % 3s [%s]\n", $i, $num_files, $pfile);

    $rid = $map{$pfile};

    $cmd = $exe_path . "sstacks -g -b $batch_id -c $out_path/$cat_file -s $out_path/$pfile -o $out_path " . join(" ", @_sstacks) . " 2>&1";
    print STDERR  "  $cmd\n";
    print $log_fh "$cmd\n";
    @results =    `$cmd` if ($dry_run == 0);
    print $log_fh @results;

    print STDERR "  Loading sstacks output to $db..." if ($sql == 1);

    if ($gzip == 1) {
	$file = "$out_path/" . $pfile . ".matches.tsv.gz";
	import_gzsql_file($log_fh, $file, "matches", 1);

    } else {
	$file = "$out_path/" . $pfile . ".matches.tsv";
	import_sql_file($log_fh, $file, "matches", 1);
    }
    print STDERR "done.\n" if ($sql == 1);

    $i++;
}

if ($data_type eq "map") {
    #
    # Generate a set of observed haplotypes and a set of markers and generic genotypes
    #
    printf(STDERR "Generating genotypes...\n");

    $cmd = $exe_path . "genotypes -b $batch_id -P $out_path -r 1 -c -s " . join(" ", @_genotypes) . " 2>&1";
    print STDERR  "$cmd\n";
    print $log_fh "$cmd\n";

    if ($dry_run == 0) {
	open($pipe_fh, "$cmd |");
	while (<$pipe_fh>) {
	    print $log_fh $_;
	}
	close($pipe_fh);
    }

    $file = "$out_path/batch_" . $batch_id . ".markers.tsv";
    import_sql_file($log_fh, $file, "markers", 1);

    $file = "$out_path/batch_" . $batch_id . ".genotypes_1.txt";
    import_sql_file($log_fh, $file, "catalog_genotypes", 1);

} else {
    printf(STDERR "Calculating population-level summary statistics\n");

    $cmd = $exe_path . "populations -b $batch_id -P $out_path -s " . join(" ", @_populations) . " 2>&1";
    print STDERR  "$cmd\n";
    print $log_fh "$cmd\n";

    if ($dry_run == 0) {
	open($pipe_fh, "$cmd |");
	while (<$pipe_fh>) {
	    print $log_fh $_;
	}
	close($pipe_fh);
    }

    $file = "$out_path/batch_" . $batch_id . ".markers.tsv";
    import_sql_file($log_fh, $file, "markers", 1);

    $file = "$out_path/batch_" . $batch_id . ".sumstats.tsv";
    import_sql_file($log_fh, $file, "sumstats", $pop_cnt+1);

    $file = "$out_path/batch_" . $batch_id . ".hapstats.tsv";
    import_sql_file($log_fh, $file, "hapstats", $pop_cnt+1);

    #
    # Import the Fst files.
    #
    my $fst_cnt = 0;
    my (@keys, $m, $n);
    @keys = sort keys %pops;
    for ($m = 0; $m < scalar(@keys); $m++) {
	for ($n = 0; $n < scalar(@keys); $n++) {
	    $file = "$out_path/batch_" . $batch_id . ".fst_" . $keys[$m] . "-" . $keys[$n] . ".tsv";

	    if (-e $file) {
		import_sql_file($log_fh, $file, "fst", 1);
		$fst_cnt++;
	    }
	}
    }
    print STDERR "Imported $fst_cnt SNP Fst file(s).\n";

    #
    # Import the Phi_st files.
    #
    $fst_cnt = 0;
    for ($m = 0; $m < scalar(@keys); $m++) {
	for ($n = 0; $n < scalar(@keys); $n++) {
	    $file = "$out_path/batch_" . $batch_id . ".phistats_" . $keys[$m] . "-" . $keys[$n] . ".tsv";

	    if (-e $file) {
		import_sql_file($log_fh, $file, "phist", 3);
		$fst_cnt++;
	    }
	}
    }
    print STDERR "Imported $fst_cnt Haplotype Fst file(s).\n";
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

print $log_fh "refmap_map.pl completed at ", strftime("%Y-%m-%d %H:%M:%S",(localtime(time))), "\n";

close($log_fh);

sub parse_population_map {
    my ($samples, $pop_ids, $pops, $grp_ids, $grps) = @_;

    my ($fh, @parts, $line, %ids, $file, $path);

    if (length($popmap_path) == 0) {
	foreach $path (@{$samples}) {
	    push(@{$pop_ids}, "1");
	    push(@{$grp_ids}, "1");
	    $pops->{"1"}++;
	    $grps->{"1"}++;
	}
	return;
    }

    open($fh, "<$popmap_path") or die("Unable to open population map, '$popmap_path', $!\n");

    while ($line = <$fh>) {
	chomp $line;
	@parts = split(/\t/, $line);

	if (scalar(@parts) > 3) {
	    die("Unable to parse population map, '$popmap_path' (map should contain no more than three columns).\n");
	}

	$ids{$parts[0]} = $parts[1];

	if (scalar(@parts) > 2) {
	    push(@{$grp_ids}, $parts[2]);
	    $grps->{$parts[2]}++;
	}
    }

    if (scalar(keys %{$grps}) == 0) {
	$grps->{"1"}++;
    }

    foreach $path (@{$samples}) {
	my ($prefix, $suffix);
	if ($path =~ /^.+\..+\.gz$/) {
	    ($prefix, $suffix) = ($path =~ /^(.+)\.(.+)\.gz$/);
	} else {
	    ($prefix, $suffix) = ($path =~ /^(.+)\.(.+)$/);
	}

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

    print STDERR "Parsed population map: ", scalar(@{$samples}), " files in ", scalar(keys %{$pops});
    scalar(keys %{$pops}) == 1 ?  print STDERR " population" : print STDERR " populations";
    print STDERR " and ", scalar(keys %{$grps});
    scalar(keys %{$grps}) == 1 ? print STDERR " group.\n" : print STDERR " groups.\n";

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

    if ($sql == 1) {
	print $log_fh "mysql --defaults-file=$cnf $db -e \"LOAD DATA LOCAL INFILE '$file' INTO TABLE $table $ignore\"\n", @results;
    }
}

sub import_gzsql_file {
    my ($log_fh, $file, $table, $skip_lines) = @_;

    my (@results, $ignore);

    $ignore = "IGNORE $skip_lines LINES" if ($skip_lines > 0);

    #
    # Get a temporary file name and create a named pipe.
    #
    my $tmpdir     = File::Spec->tmpdir();
    my $named_pipe = mktemp($tmpdir . "/denovo_map_XXXXXX");
    if ($sql == 1 && $dry_run == 0) {
	mkfifo($named_pipe, 0700) || die("Unable to create named pipe for loading gzipped data: $named_pipe, $!");
	print $log_fh "Streaming $file into named pipe $named_pipe.\n";
    }

    #
    # Dump our gzipped data onto the named pipe.
    #
    system("gunzip -c $file > $named_pipe &") if ($sql == 1 && $dry_run == 0);

    @results = `mysql --defaults-file=$cnf $db -e "LOAD DATA LOCAL INFILE '$named_pipe' INTO TABLE $table $ignore"` if ($sql == 1 && $dry_run == 0);

    if ($sql == 1) {
	print $log_fh "mysql --defaults-file=$cnf $db -e \"LOAD DATA LOCAL INFILE '$named_pipe' INTO TABLE $table $ignore\"\n", @results;
    }

    #
    # Remove the pipe.
    #
    unlink($named_pipe) if ($sql == 1 && $dry_run == 0);
}

sub parse_command_line {
    my $arg;

    while (@ARGV) {
	$_ = shift @ARGV;
	if    ($_ =~ /^-p$/) { push(@parents, shift @ARGV); }
	elsif ($_ =~ /^-r$/) { push(@progeny, shift @ARGV); }
	elsif ($_ =~ /^-s$/) { push(@samples, shift @ARGV); }
	elsif ($_ =~ /^-o$/) { $out_path    = shift @ARGV; }
	elsif ($_ =~ /^-D$/) { $desc        = shift @ARGV; }
	elsif ($_ =~ /^-e$/) { $exe_path    = shift @ARGV; }
	elsif ($_ =~ /^-b$/) { $batch_id    = shift @ARGV; }
	elsif ($_ =~ /^-i$/) { $sample_id   = shift @ARGV; }
	elsif ($_ =~ /^-a$/) { $date        = shift @ARGV; }
	elsif ($_ =~ /^-S$/) { $sql         = 0; }
	elsif ($_ =~ /^-B$/) { $db          = shift @ARGV; }
	elsif ($_ =~ /^-d$/) { $dry_run++; }
	elsif ($_ =~ /^-O$/) { 
	    $popmap_path = shift @ARGV;
	    push(@_populations, "-M " . $popmap_path); 

	} elsif ($_ =~ /^-A$/) { 
	    $arg = shift @ARGV;
	    push(@_genotypes, "-t " . $arg); 

	    $arg = lc($arg);
	    if ($arg ne "gen" && $arg ne "cp" && $arg ne "f2" && $arg ne "bc1" && $arg ne "dh") {
		print STDERR "Unknown genetic mapping cross specified: '$arg'\n";
		usage();
	    }

	} elsif ($_ =~ /^-T$/) {
	    $arg = shift @ARGV;
	    push(@_pstacks, "-p " . $arg); 
	    push(@_cstacks, "-p " . $arg); 
	    push(@_sstacks, "-p " . $arg); 
	    push(@_populations, "-t " . $arg); 

	} elsif ($_ =~ /^-m$/) { 
	    push(@_pstacks,   "-m " . shift @ARGV); 

	} elsif ($_ =~ /^-n$/) { 
	    push(@_cstacks, "-n " . shift @ARGV); 

	} elsif ($_ =~ /^--bound_low$/) { 
	    push(@_pstacks, "--bound_low " . shift @ARGV); 
	    push(@_pstacks, "--model_type bounded");

	} elsif ($_ =~ /^--bound_high$/) { 
	    push(@_pstacks, "--bound_high " . shift @ARGV); 
	    push(@_pstacks, "--model_type bounded");
	} elsif ($_ =~ /^--alpha$/) { 
	    push(@_pstacks, "--alpha " . shift @ARGV); 

	} elsif ($_ =~ /^-X$/) {
	    #
	    # Pass an arbitrary command-line option to a pipeline program.
	    #
	    # Command line option must be of the form '-X "program:option"'
	    #
	    $arg = shift @ARGV;
	    my ($prog, $opt) = ($arg =~ /^(\w+):(.+)$/);

	    if ($prog eq "pstacks") {
		push(@_pstacks, $opt); 

	    } elsif ($prog eq "cstacks") {
		push(@_cstacks, $opt); 

	    } elsif ($prog eq "sstacks") {
		push(@_sstacks, $opt); 

	    } elsif ($prog eq "genotypes") {
		push(@_genotypes, $opt); 

	    } elsif ($prog eq "populations") {
		push(@_populations, $opt); 
	    } else {
		print STDERR "Unknown pipeline program, '$arg'\n";
		usage();
	    }
	}
	elsif ($_ =~ /^-v$/) { version(); exit(); }
	elsif ($_ =~ /^-h$/) { usage(); }
	else {
	    print STDERR "Unknown command line option: '$_'\n";
	    usage();
	}
    }

    $exe_path = $exe_path . "/"          if (substr($exe_path, -1) ne "/");
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
}

sub version {
    print STDERR "ref_map.pl ", stacks_version, "\n";
}

sub usage {
    version();

    print STDERR <<EOQ; 
ref_map.pl -p path -r path [-s path] -o path [-n mismatches] [-m min_cov] [-T num_threads] [-A type] [-O popmap] [-B db -b batch_id -D "desc" -a yyyy-mm-dd] [-S -i id] [-e path] [-d] [-h]
    p: path to a Bowtie/SAM/BAM file containing parent sequences from a mapping cross.
    r: path to a Bowtie/SAM/BAM file containing progeny sequences from a mapping cross.
    s: path to a Bowtie/SAM/BAM file containing an individual sample from a population.
    o: path to write pipeline output files.
    A: if processing a genetic map, specify the cross type, 'CP', 'F2', 'BC1', 'DH', or 'GEN'.
    n: specify the number of mismatches allowed between loci when building the catalog (default 1).
    T: specify the number of threads to execute.
    m: specify the minimum depth of coverage to report a stack in pstacks (default 1).
    O: if analyzing one or more populations, specify a pOpulation map
    e: executable path, location of pipeline programs.
    d: perform a dry run. Do not actually execute any programs, just print what would be executed.
    h: display this help message.

   Database options:
    B: specify a database to load data into.
    b: batch ID representing this dataset in the database.
    D: batch description
    a: batch run date, yyyy-mm-dd
    S: disable recording SQL data in the database.
    i: starting sample_id, this is determined automatically if database interaction is enabled.

   SNP Model Options (these options are passed on to pstacks):
    --bound_low <num>: lower bound for epsilon, the error rate, between 0 and 1.0 (default 0).
    --bound_high <num>: upper bound for epsilon, the error rate, between 0 and 1.0 (default 1).
    --alpha <num>: chi square significance level required to call a heterozygote or homozygote, either 0.1, 0.05 (default), 0.01, or 0.001.

  Arbitrary command line options:
    -X "program:option": pass a command line option to one of the pipeline components, e.g.'-X "sstacks:-x"'.

EOQ

exit(0);
}
