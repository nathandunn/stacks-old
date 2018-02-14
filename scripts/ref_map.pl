#!/usr/bin/env perl
#
# Copyright 2010-2018, Julian Catchen <jcatchen@illinois.edu>
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

use strict;
use POSIX;
use File::Temp qw/ mktemp /;
use File::Spec;
use constant stacks_version => "_VERSION_";

use constant true  => 1;
use constant false => 0;

my $dry_run      = false;
my $sql          = false;
my $create_db    = false;
my $overw_db     = false;
my $mysql_config = "_PKGDATADIR_" . "sql/mysql.cnf";
my $mysql_tables = "_PKGDATADIR_" . "sql/stacks.sql";
my $exe_path     = "_BINDIR_";
my $out_path     = "";
my $popmap_path  = "";
my $sample_path  = "";
my $db           = "";
my $batch_id     = 1;
my $sample_id    = 1;
my $desc         = ""; # Database description of this dataset
my $date         = ""; # Date relevent to this data, formatted for SQL: 2009-05-31
my $gzip         = false;
my $time         = "";

my @parents;
my @progeny;
my @samples;

my (@_gstacks, @_populations);

my $cmd_str = $0 . " " . join(" ", @ARGV);

parse_command_line();

my $cnf = (-e $ENV{"HOME"} . "/.my.cnf") ? $ENV{"HOME"} . "/.my.cnf" : $mysql_config;

#
# Check for the existence of the necessary pipeline programs
#
foreach my $prog ("gstacks", "populations") {
    die "Unable to find '" . $exe_path . $prog . "'.\n" if (!-e $exe_path . $prog);
}

my ($log, $log_fh, $sample);

my (@sample_list, %pop_ids, %pops, %grp_ids, %grps, %sample_ids);

parse_population_map(\@sample_list, \%pop_ids, \%pops, \%grp_ids, \%grps);

initialize_samples(\@parents, \@progeny, \@samples, \@sample_list, \%pop_ids, \%grp_ids);

#
# Open the log file
#
$log = "$out_path/ref_map.log";
open($log_fh, ">$log") or die("Unable to open log file '$log'; $!\n");

print $log_fh
    "ref_map.pl version ", stacks_version, " started at ", strftime("%Y-%m-%d %H:%M:%S", (localtime(time))), "\n",
    $cmd_str, "\n";

initialize_database($log_fh, \@parents, \@progeny, \@samples, \%sample_ids) if ($sql == true);

execute_stacks($log_fh, $sample_id, \@parents, \@progeny, \@samples, \%sample_ids);

load_sql_data($log_fh, \%pops, \@parents, \@progeny, \@samples) if ($sql == true);

print $log_fh "\nref_map.pl completed at ", strftime("%Y-%m-%d %H:%M:%S", (localtime(time))), "\n";
close($log_fh);

sub check_return_value {
    # $? is a 16 bit int. Exit code is given by `$? & 255` if the process was
    # terminated by a signal, and by `$? >> 8` if it exited normally.
    my ($rv, $log_fh) = @_;
    if ($rv != 0) {
        my $code = ($rv >> 8) & 127;
        if ($rv & 255 || ($rv >> 8) > 127) {
            $code += 128;
        }
        my $msg = "\nref_map.pl: Aborted because the last command failed ($code";
        if ($code == 129 || $code == 130 || $code == 131) {
            $msg .= "/interrupted";
        } elsif ($code == 137 || $code == 143) {
            $msg .= "/killed";
        } elsif ($code == 134) {
            $msg .= "/SIGABRT";
        } elsif ($code == 139) {
            $msg .= "/segmentation fault";
        }
        $msg .= ")";
        print $log_fh ($msg . ".\n");
        print STDERR ($msg . "; see log file.\n");
        exit 1;
    }
}

sub execute_stacks {
    my ($log_fh, $sample_id, $parents, $progeny, $samples, $sample_ids) = @_;

    my (@results, @depths_of_cov);
    my ($pop_cnt, $sample, $num_files, $i, $cmd, $pipe_fh, $path, $cat_file);

    #
    # Call genotypes.
    #
    print STDERR "Calling variants, genotypes and haplotypes...\n";
    print $log_fh "\ngstacks\n==========\n";

    $cmd = $exe_path . "gstacks -I $sample_path -M $popmap_path -O $out_path";
    foreach (@_gstacks) {
        $cmd .= " " . $_;
    }
    print STDERR  "  $cmd\n\n";
    print $log_fh "$cmd\n\n";
    if (!$dry_run) {
        open($pipe_fh, "$time $cmd 2>&1 |");
        while (<$pipe_fh>) {
            print $log_fh $_;
        }
        close($pipe_fh);
        check_return_value($?, $log_fh);
    }

    printf(STDERR "Calculating population-level summary statistics\n");
    print $log_fh "\npopulations\n==========\n";

    $cmd = $exe_path . "populations" . " -P $out_path " . join(" ", @_populations);
    print STDERR  "  $cmd\n\n";
    print $log_fh "$cmd\n\n";

    if (!$dry_run) {
        open($pipe_fh, "$time $cmd 2>&1 |");
        while (<$pipe_fh>) {
            print $log_fh $_;
        }
        close($pipe_fh);
        check_return_value($?, $log_fh);
    }
    
    print STDERR  "ref_map.pl is done.\n";
    print $log_fh "ref_map.pl is done.\n";
}

sub parse_population_map {
    my ($sample_list, $pop_ids, $pops, $grp_ids, $grps) = @_;

    my ($fh, @parts, $line, $sample);

    return if (length($popmap_path) == 0);

    open($fh, "<$popmap_path") or die("Unable to open population map, '$popmap_path', $!\n");

    while ($line = <$fh>) {
        chomp $line;

        next if ($line =~ /^\s*#/);

        @parts = split(/\t/, $line);

        if (scalar(@parts) > 3) {
            die("Unable to parse population map, '$popmap_path' (map should contain no more than three columns).\n");
        }

        push(@{$sample_list}, $parts[0]);

        $pop_ids->{$parts[0]} = $parts[1];
        $pops->{$parts[1]}++;

        if (scalar(@parts) > 2) {
            $grp_ids->{$parts[0]} = $parts[2];
            $grps->{$parts[2]}++;
        }
    }

    if (scalar(keys %{$grps}) == 0) {
        $grps->{"1"}++;

        foreach $sample (@{$sample_list}) {
            $grp_ids->{$sample} = "1";
        }
    }

    print STDERR "Parsed population map: ", scalar(@{$sample_list}), " files in ", scalar(keys %{$pops});
    scalar(keys %{$pops}) == 1 ?  print STDERR " population" : print STDERR " populations";
    print STDERR " and ", scalar(keys %{$grps});
    scalar(keys %{$grps}) == 1 ? print STDERR " group.\n" : print STDERR " groups.\n";

    close($fh);
}

sub initialize_samples {
    my ($parents, $progeny, $samples, $sample_list, $pop_ids, $grp_ids) = @_;

    my ($local_gzip, $file, $prefix, $suffix, $path, $found, $i);

    if (scalar(@{$sample_list}) > 0 && scalar(@{$samples}) == 0) {
        my @suffixes = ("bam");
        my @fmts     = ("bam");

        #
        # If a population map was specified and no samples were provided on the command line.
        #
        foreach $sample (@{$sample_list}) {
            $found = false;

            for ($i = 0; $i < scalar(@suffixes); $i++) {
                $path = $sample_path . $sample . "." . $suffixes[$i];
                if (-e $path) {

                    $gzip = true if ($i == 1);

                    push(@{$samples}, {'path'   => $sample_path,
                                       'file'   => $sample,
                                       'suffix' => $suffixes[$i],
                                       'type'   => "sample",
                                       'fmt'    => $fmts[$i]});
                    $found = true;
                    last;
                }
            }

            if ($found == false) {
                die("Error: Failed to open '$sample_path$sample.bam'.\n");
            }
        }
    }

    #
    # If a population map was specified, make sure all samples in the list were found (and vice versa) and assign popualtion IDs.
    #
    if (scalar(@{$sample_list}) > 0) {

        my %sample_hash;

        foreach $sample (@{$samples}) {
            $sample_hash{$sample->{'file'}}++;

            if (!defined($pop_ids->{$sample->{'file'}})) {
                die("Unable to find an entry for '" . $sample->{'file'} . "' in the population map, '$popmap_path'.\n");
            } else {
                $sample->{'pop_id'} = $pop_ids->{$sample->{'file'}};
            }
            if (!defined($grp_ids->{$sample->{'file'}})) {
                die("Unable to find an entry for '" . $sample->{'file'} . "' in the population map, '$popmap_path'.\n");
            } else {
                $sample->{'grp_id'} = $grp_ids->{$sample->{'file'}};
            }
        }

        foreach $sample (@{$sample_list}) {
            if (!defined($sample_hash{$sample})) {
                die("Unable to find a file corresponding to the population map entry '" . $sample . "' in the population map, '$popmap_path'.\n");
            }
        }

    } else {
        foreach $sample (@{$parents}, @{$progeny}, @{$samples}) {
            $sample->{'pop_id'} = "1";
            $sample->{'grp_id'} = "1";
            $pop_ids->{$sample->{'file'}} = $sample->{'pop_id'};
            $grp_ids->{$sample->{'file'}} = $sample->{'grp_id'};
        }
    }

    #
    # Check that no duplicate files were specified.
    #
    my (%files, $file);
    foreach $file (@{$parents}, @{$progeny}, @{$samples}) {
        $files{$file}++;
    }
    foreach $file (keys %files) {
        if ($files{$file} > 1) {
            die("A duplicate file was specified which may create undefined results, '$file'\n");
        }
    }

    print STDERR "Found ", scalar(@{$parents}), " parental file(s).\n\n" if (scalar(@{$parents}) > 0);
    print STDERR "Found ", scalar(@{$progeny}), " progeny file(s).\n\n" if (scalar(@{$progeny}) > 0);
    print STDERR "Found ", scalar(@{$samples}), " sample file(s).\n\n" if (scalar(@{$samples}) > 0);
}

sub initialize_database {
    my ($log_fh, $parents, $progeny, $samples, $sample_ids) = @_;

    my (@results, $sample_id, $sample);

    print $log_fh "Initializing the database...\n";

    #
    # Create the database.
    #
    if ($create_db) {
        #
        # Check that the database doesn't already exist.
        #
        if ($dry_run == false) {
            @results = `mysql --defaults-file=$cnf -N -B -e "SHOW DATABASES LIKE '$db'"`;
            if (scalar(@results) > 0 && $overw_db == false) {
                die("Unable to create database '$db', it already exists.\n");
            }
        }

        if ($overw_db == true) {
            `mysql --defaults-file=$cnf -N -B -e "DROP DATABASE IF EXISTS $db"` if ($dry_run == false);
            print $log_fh "mysql --defaults-file=$cnf -N -B -e \"DROP DATABASE IF EXISTS '$db'\"\n";
        }

        `mysql --defaults-file=$cnf -e "CREATE DATABASE $db"` if ($dry_run == false);
        print $log_fh "mysql --defaults-file=$cnf $db -e \"CREATE DATABASE $db\"\n";
        `mysql --defaults-file=$cnf $db < $mysql_tables` if ($dry_run == false);
        print $log_fh "mysql --defaults-file=$cnf $db < $mysql_tables\n";
    }

    #
    # Set the SQL Batch ID for this set of loci, along with description and date of
    # sequencing. Insert this batch data into the database.
    #
    `mysql --defaults-file=$cnf $db -e "INSERT INTO batches SET id=$batch_id, description='$desc', date='$date', type='population'"` if ($dry_run == false);

    print $log_fh "mysql --defaults-file=$cnf $db -e \"INSERT INTO batches SET id=$batch_id, description='$desc', date='$date', type='population'\"\n";

    print $log_fh "Loading sample data into the MySQL database...\n";

    my $i = 1;

    foreach $sample (@{$parents}, @{$progeny}, @{$samples}) {

        if ($dry_run == false) {
            `mysql --defaults-file=$cnf $db -e "INSERT INTO samples SET sample_id=$i, batch_id=$batch_id, type='$sample->{'type'}', file='$sample->{'file'}', pop_id='$sample->{'pop_id'}', group_id='$sample->{'grp_id'}'"`;
            @results = `mysql --defaults-file=$cnf $db -N -B -e "SELECT id FROM samples WHERE sample_id=$i AND batch_id=$batch_id AND type='$sample->{'type'}' AND file='$sample->{'file'}'"`;
            chomp $results[0];
            $sample_id = $results[0];

            #
            # Save the sample ID to use when running pstacks.
            #
            $sample_ids->{$sample->{'file'}} = $sample_id;
        }
        print $log_fh "mysql --defaults-file=$cnf $db -e \"INSERT INTO samples SET sample_id=$i, batch_id=$batch_id, type='$sample->{'type'}', file='$sample->{'file'}', pop_id='$sample->{'pop_id'}', group_id='$sample->{'grp_id'}'\"\n";

        $i++;
    }

    print $log_fh "\n";
}

sub load_sql_data {
    my ($log_fh, $pops, $parents, $progeny, $samples) = @_;

    my ($pop_cnt, $sample, $num_files, $i, $file);

    print STDERR "\nComputation is complete, loading results to the database '$db'.\n";

    #
    # Load the outputs generated by populations to the database.
    #
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
    print STDERR "Imported $fst_cnt Fst file(s).\n";

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

    print $log_fh "\n";

    #
    # Index the radtags database
    #
    my ($cmd, @results);
    print STDERR "Indexing the database...\n";
    $cmd = $exe_path . "index_radtags.pl -D $db -t -c 2>&1";
    print STDERR  "$cmd\n";
    print $log_fh "$cmd\n";
    @results =    `$time $cmd` if ($dry_run == false);
    print $log_fh @results;
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

sub write_depths_of_cov {
    my ($depths, $log_fh) = @_;

    print STDERR "\nDepths of Coverage for Processed Samples:\n";
    print $log_fh "\nDepths of Coverage for Processed Samples:\n";

    foreach $a (@{$depths}) {
        print STDERR  $a->[0], ": ", $a->[1], "x\n";
        print $log_fh $a->[0], ": ", $a->[1], "x\n";
    }
}

sub import_sql_file {
    my ($log_fh, $file, $table, $skip_lines) = @_;

    my (@results, $ignore);

    $ignore = "IGNORE $skip_lines LINES" if ($skip_lines > 0);

    @results = `mysql --defaults-file=$cnf $db -e "LOAD DATA LOCAL INFILE '$file' INTO TABLE $table $ignore"` if ($sql == true && $dry_run == false);

    if ($sql == true) {
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
    if ($sql == true && $dry_run == false) {
        mkfifo($named_pipe, 0700) || die("Unable to create named pipe for loading gzipped data: $named_pipe, $!");
        print $log_fh "Streaming $file into named pipe $named_pipe.\n";
    }

    #
    # Dump our gzipped data onto the named pipe.
    #
    system("gunzip -c $file > $named_pipe &") if ($sql == true && $dry_run == false);

    @results = `mysql --defaults-file=$cnf $db -e "LOAD DATA LOCAL INFILE '$named_pipe' INTO TABLE $table $ignore"` if ($sql == true && $dry_run == false);

    if ($sql == 1) {
        print $log_fh "mysql --defaults-file=$cnf $db -e \"LOAD DATA LOCAL INFILE '$named_pipe' INTO TABLE $table $ignore\"\n", @results;
    }

    #
    # Remove the pipe.
    #
    unlink($named_pipe) if ($sql == true && $dry_run == false);
}

sub parse_command_line {
    my ($arg);

    while (@ARGV) {
        $_ = shift @ARGV;
        if    ($_ =~ /^-v$/) { version(); exit 1; }
        elsif ($_ =~ /^-h$/) { usage(); }
        elsif ($_ =~ /^-d$/ || $_ =~ /^--dry-run$/) { $dry_run = true; }
        elsif ($_ =~ /^-o$/) { $out_path  = shift @ARGV; }
        elsif ($_ =~ /^-e$/) { $exe_path  = shift @ARGV; }
        elsif ($_ =~ /^--samples$/) {
            $sample_path = shift @ARGV;
        } elsif ($_ =~ /^-O$/ || $_ =~ /^--popmap$/) {
            $popmap_path = shift @ARGV;
            push(@_populations, "-M " . $popmap_path);

        } elsif ($_ =~ /^--unpaired$/) {
            push(@_gstacks, "--unpaired");

        } elsif ($_ =~ /^--ignore-pe-reads$/) {
            push(@_gstacks, "--ignore-pe-reads");

        } elsif ($_ =~ /^-T$/) {
            $arg = shift @ARGV;
            push(@_gstacks,     "-t " . $arg);
            push(@_populations, "-t " . $arg);

        } elsif ($_ =~ /^--var-alpha$/) {
            $arg = shift @ARGV;
            push(@_gstacks, "--var_alpha " . $arg);

        } elsif ($_ =~ /^--gt-alpha$/) {
            $arg = shift @ARGV;
            push(@_gstacks, "--gt_alpha " . $arg);

        } elsif ($_ =~ /^-X$/) {
            #
            # Pass an arbitrary command-line option to a pipeline program.
            #
            # Command line option must be of the form '-X "program:option"'
            #
            $arg = shift @ARGV;
            my ($prog, $opt) = ($arg =~ /^(\w+):(.+)$/);

            if ($prog eq "gstacks") {
                push(@_gstacks, $opt);

            } elsif ($prog eq "populations") {
                push(@_populations, $opt);
            } else {
                print STDERR "Unknown pipeline program, '$arg'\n";
                usage();
            }
        } elsif ($_ =~ /^--time-components$/) {
            $time = '/usr/bin/time';
            if (! -e $time) {
                die "Error: '$time': No such file or directory.\n";
            }
        } else {
            print STDERR "Unknown command line option: '$_'\n";
            usage();
        }
    }

    $exe_path = $exe_path . "/"          if (substr($exe_path, -1) ne "/");
    $out_path = substr($out_path, 0, -1) if (substr($out_path, -1) eq "/");

    if ($sql == true && length($date) == 0) {
        $date = strftime("%Y-%m-%d", (localtime(time)));
    }

    if (length($popmap_path) == 0) {
        print STDERR "You must specify a population map that lists your sample names (--popmap).\n";
        usage();
    }

    if (length($sample_path) == 0) {
        print STDERR "You must specify the path to the directory containing the samples (--samples).\n";
        usage();
    }

    if (length($sample_path) > 0) {
        $sample_path .= "/" if (substr($sample_path, -1) ne "/");
    }
}

sub version {
    print STDERR "ref_map.pl ", stacks_version, "\n";
}

sub usage {
    version();

    print STDERR <<EOQ;
ref_map.pl --samples dir --popmap path [-s spacer] [--paired] -o dir [-X prog:"opts" ...]

  Input/Output files:
    --samples: path to the directory containing the samples BAM (or SAM) alignment files.
    --popmap: path to a population map file (format is "<name> TAB <pop>", one sample per line).
    s: spacer for file names: by default this is empty and the program looks for files
       named 'SAMPLE_NAME.bam'; if this option is given the program looks for files
       named 'SAMPLE_NAME.SPACER.bam'.
    o: path to an output directory.
    --unpaired: ignore read pairing (for ddRAD; treat READ2's as if they were READ1's)
    --ignore-pe-reads: ignore paired-end reads even if present in the input

  General options:
    X: additional options for specific pipeline components, e.g. -X "populations: -p 3 -r 0.50"
    T: the number of threads/CPUs to use (default: 1).
    d: Dry run. Do not actually execute anything, just print the individual pipeline commands 
       that would be executed.

  SNP model options:
    --var-alpha: significance level at which to call variant sites (for gstacks; default: 0.05).
    --gt-alpha: significance level at which to call genotypes (for gstacks; default: 0.05).

  Miscellaneous:
    --time-components (for benchmarking)
EOQ
    exit 1;
}
