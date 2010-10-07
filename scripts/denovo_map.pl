#!/usr/bin/perl
#
# By Julian Catchen <catchen@cs.uoregon.edu>
#

use strict;
use Bio::SeqIO;

my $debug       = 0;
my $sql         = 1;
my $exe_path    = $ENV{'HOME'} . "/research/solexa/radtags/bin";
my $out_path    = "";
my $white_list  = "";
my $db          = "";
my $rep_tags    = 0;
my $min_cov     = 0;
my $fuzzy_match = 0;
my $cov_scale   = 0;
my $batch_id    = 0;
my $sample_id   = 1;
my $desc        = ""; #"Lepisosteus oculatus RAD-Tag Samples";
my $date        = ""; #"2009-05-31";

my @parents;
my @progeny;

parse_command_line();

my ($i, $log, $log_fh, $pfile, $rfile, $file, $num_files, $parent, $sample, %map);

$i         = 1;
$num_files = scalar(@parents) + scalar(@progeny);

if ($sql) {
    #
    # SQL Batch ID for this set of Radtags, along with description and date of 
    # sequencing. Insert this batch data into the database.
    # 
    `mysql $db -e "INSERT INTO batches SET id=$batch_id, description='$desc', date='$date'"`;
}

my (@types, $type);
foreach $parent (@parents) {
    push(@types, "parent");
}
foreach $parent (@progeny) {
    push(@types, "progeny");
}

my (@results, $minc, $rrep, $cmd, $cscale, $threads, $fuzzym);

$minc    = $min_cov   > 0 ? "-m $min_cov"   : "";
$cscale  = $cov_scale > 0 ? "-S $cov_scale" : "";
$threads = "-p 15"; 
$fuzzym  = "-n $fuzzy_match";

#
# Open the log file
#
$log = "$out_path/denovo_map.log";
open($log_fh, ">$log") or die("Unable to open log file '$log'; $!\n");

foreach $sample (@parents, @progeny) {

    ($pfile) = ($sample =~ /^.*\/(.+)\.fastq_1$/);

    $type = shift @types;

    printf("Identifying unique radtags; file % 3s of % 3s [%s]\n", $i, $num_files, $pfile);
    printf($log_fh "Identifying unique radtags; file % 3s of % 3s [%s]\n", $i, $num_files, $pfile);

    if ($sql) {
	`mysql $db -e "INSERT INTO samples SET sample_id=$i, batch_id=$batch_id, type='$type', file='$pfile'"`;
	@results = `mysql $db -N -B -e "SELECT id FROM samples WHERE sample_id=$i AND batch_id=$batch_id AND type='$type' AND file='$pfile'"`;
	chomp $results[0];
	$sample_id = $results[0];
    } else {
	print STDERR "mysql $db -e \"INSERT INTO samples SET sample_id=$i, batch_id=$batch_id, type='$type', file='$pfile'\"\n";
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

    $cmd = "$exe_path/ustacks -t fastq -f $sample -o $out_path -b $batch_id -i $sample_id $rrep $minc $cscale $threads 2>&1";
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
    ($pfile)  = ($sample =~ /^.*\/(.+)\.fastq_1$/);
    $parents .= "-s $out_path/$pfile -S " . $map{$pfile} . " ";
}

$cat_file = "batch_" . $batch_id;
$cmd      = "$exe_path/cstacks -b $batch_id -o $out_path $parents $threads $fuzzym 2>&1";
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

    ($rfile) = ($sample =~ /^.*\/(.+)\.fastq_1$/);
    printf(STDERR "Matching RAD-Tags to catalog; file % 3s of % 3s [%s]\n", $i, $num_files, $rfile);

    $rid = $map{$rfile};

    $cmd = "$exe_path/sstacks -b $batch_id -c $out_path/$cat_file -s $out_path/$rfile -S $rid -o $out_path $threads 2>&1";
    print STDERR  "$cmd\n";
    print $log_fh "$cmd\n";
    @results =    `$cmd`;
    print $log_fh @results;

    $file = "$out_path/" . $rfile . ".matches.tsv";
    import_sql_file($file, "matches");

    $i++;
}

#
# Search for markers
#
print STDERR "$exe_path/markers.pl -b $batch_id -o $out_path\n";
`$exe_path/markers.pl -D $db -b $batch_id -o $out_path`;
$file = "$out_path/batch_" . $batch_id . ".markers.tsv";
import_sql_file($file, "markers");

#
# Index the radtags database
#
print STDERR "$exe_path/index-radtags.pl -D $db -t -c\n";
`$exe_path/index-radtags.pl -D $db -t -c`;

close($log_fh);

sub import_sql_file {
    my ($file, $table) = @_;

    my (@results);

    @results = `mysql $db -e "LOAD DATA LOCAL INFILE '$file' INTO TABLE $table"` if ($sql);
    print STDERR "mysql $db -e \"LOAD DATA LOCAL INFILE '$file' INTO TABLE $table\"\n", @results, "\n";
}

sub parse_command_line {
    while (@ARGV) {
	$_ = shift @ARGV;
	if    ($_ =~ /^-p$/) { push(@parents, shift @ARGV); }
	elsif ($_ =~ /^-r$/) { push(@progeny, shift @ARGV); }
       	elsif ($_ =~ /^-t$/) { $rep_tags++; }
	elsif ($_ =~ /^-o$/) { $out_path    = shift @ARGV; }
	elsif ($_ =~ /^-m$/) { $min_cov     = shift @ARGV; }
        elsif ($_ =~ /^-n$/) { $fuzzy_match = shift @ARGV; }
	elsif ($_ =~ /^-c$/) { $cov_scale   = shift @ARGV; }
	elsif ($_ =~ /^-D$/) { $desc        = shift @ARGV; }
	elsif ($_ =~ /^-b$/) { $batch_id    = shift @ARGV; }
	elsif ($_ =~ /^-s$/) { $sample_id   = shift @ARGV; }
	elsif ($_ =~ /^-a$/) { $date        = shift @ARGV; }
	elsif ($_ =~ /^-S$/) { $sql         = 0; }
	elsif ($_ =~ /^-B$/) { $db          = shift @ARGV; }
	elsif ($_ =~ /^-d$/) { $debug++; }
	elsif ($_ =~ /^-h$/) { usage(); }
	else {
	    print STDERR "Unknown command line option: '$_'\n";
	    usage();
	}
    }

    $out_path = substr($out_path, 0, -1) if (substr($out_path, -1) eq "/");

    if (scalar(@parents) == 0) {
	print STDERR "You must specify at least one parent file.\n";
	usage();
    }
}

sub usage {
    print STDERR <<EOQ; 
denovo_map.pl -p path -r path -o path [-t] [-m min_cov] [-n mismatches] [-c scale] [-D desc] [-b batch_id] [-s num] [-a yyyy-mm-dd] [-S] [-d] [-h]
    p: path to a FASTQ file containing parent sequences.
    r: path to a FASTQ file containing progeny sequences.
    o: path to output the cleaned files.
    b: batch ID
    s: starting sample_id, this is determined automatically if database interaction is enabled.
    S: disable recording SQL data in the database.
    D: batch description
    a: batch run date, yyyy-mm-dd
    m: specify a minimum depth of coverage to merge stacks.
    n: specify the number of mismatches allowed between loci when building the catalog.
    c: coverage scaling factor affecting when tags are deleveraged or removed (between 0 and 1).
    t: remove, or break up, highly repetitive RAD-Tags.
    B: specify a database to load data into.
    h: display this help message.
    d: turn on debug output.

EOQ

exit(0);
}
