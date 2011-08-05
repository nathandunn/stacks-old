#!/usr/bin/perl
#
# Copyright 2011, Julian Catchen <jcatchen@uoregon.edu>
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
use constant stacks_version => "_VERSION_";

my $debug        = 0;
my $paired       = 0;
my $amos         = 0;
my $ins_len_dist = 0;
my $single_path  = "";
my $paired_path  = "";
my $sanger_path  = "";
my $out_path     = ".";
my $white_list   = "";
my $hash_len     = 27;
my $insert_len   = 0;
my $exp_cov      = 0;
my $cov_cut      = 0.0;
my $min_len      = 0;
my $read_trk     = 0;
my $velveth      = "velveth";
my $velvetg      = "velvetg";
my $parse_afg    = "./parse-afg.pl";

parse_command_line();

my (@locus_files, $num_files, $file, $input_file, $output_file, $hres_file, $gres_file, $ace_file, $dist_file, $barcode);

build_file_list(\@locus_files);
$num_files = scalar(@locus_files);

my $i = 1;

my ($sing_data, $pair_data, $sang_data, $ins, $cov, $afg, $cut, $min, $read);

$ins  = $paired   ? "-ins_length $insert_len"   : "-ins_length auto";
$cov  = $paired   ? "-exp_cov $exp_cov"         : "-exp_cov auto";
$cut  = $cov_cut  ? "-cov_cutoff $cov_cut"      : "-cov_cutoff auto";
$min  = $min_len  ? "-min_contig_lgth $min_len" : "";
$read = $read_trk ? "-read_trkg yes"            : "";
$afg  = $amos || $ins_len_dist ? "-amos_file yes" : "";

foreach $file (@locus_files) {

    ($file) = ($file =~ /(.+)\.fa/); 
    $output_file = $out_path . "/" . $file;
    $hres_file   = $out_path . "/" . $file . "-h.output";
    $gres_file   = $out_path . "/" . $file . "-g.output";
    $ace_file    = $out_path . "/" . $file . "/" . $file . ".ace";
    $dist_file   = $out_path . "/" . $file . "-calc_ins_len.tsv";

    $sing_data = length($single_path) > 0 ? '-short '       . $single_path . "/" . $file . ".fa" : "";
    $pair_data = length($paired_path) > 0 ? '-shortPaired ' . $paired_path . "/" . $file . ".fa" : "";
    $sang_data = length($sanger_path) > 0 ? '-long '        . $sanger_path . "/" . $file . ".fa" : "";

    #
    # Write out the parameters for this assembly
    #
    open(PARAM, "> $out_path/parameters.txt") or die("Unable to open parameter file: $!\n");

    print PARAM
	"Single-end Data: ", $sing_data, "\n",
	"Paired-end Data: ", $pair_data, "\n",
	"Sanger Data:     ", $sang_data, "\n",
	"Hash Length:     ", $hash_len, "\n",
	"Insert Length:   ", $ins, "\n",
	"Coverage:        ", $cov, "\n",
	"Coverage Cutoff: ", $cut, "\n",
	"Miniumum contig length: ", $min, "\n",
	"Read Tracking: ", $read, "\n";

    close(PARAM);

    print STDERR "Assembling locus '$file'; run $i of $num_files\n";

    # Execute velveth to build hash table, then velvetg to assemble
    print STDERR "$velveth $output_file $hash_len -fasta $sing_data $pair_data $sang_data &> $hres_file\n" if ($debug);
    `$velveth $output_file $hash_len -fasta $sing_data $pair_data $sang_data &> $hres_file`;

    print STDERR "$velvetg $output_file $ins $cov $cut $min $read $afg &> $gres_file\n" if ($debug);
    `$velvetg $output_file $ins $cov $cut $min $read $afg &> $gres_file`;

    print STDERR "amos2ace $output_file/velvet_asm.afg -o $ace_file\n" if ($amos && $debug);
    #`amos2ace $output_file/velvet_asm.afg -o $ace_file` if ($amos);

    if ($paired && $ins_len_dist) {
	if ($paired_path =~ /.+\/\d+.?\d*\-\d+.?\d*\-\d+.?\d*\/01\/paired$/) {
	    print STDERR "Calculating insert length distribution\n";
	    print STDERR "$parse_afg -p $output_file/velvet_asm.afg > $dist_file\n" if ($debug);
	    `$parse_afg -i -p $output_file/velvet_asm.afg -c $output_file/pairs.fasta > $dist_file`;

	    unlink "$output_file/velvet_asm.afg" if (!$amos);
	} else {
	    print STDERR "Skipping insert length distribution calculation.\n";
	}
    }

    $i++;
}

sub build_file_list {
    my ($files) = @_;

    my (@ls, $line, $file, $path);

    # Load a white list of files to process if it is supplied.
    my @wl;
    if (length($white_list) > 0) {
	load_white_list(\@wl);
    }

    $path = length($paired_path) > 0 ? $paired_path : $single_path;

    @ls = `ls -1 $path/`;

    foreach $line (@ls) {
	chomp $line;

	next if (length($line) == 0);
        next if ($line !~ /.+\.fa$/ && $line !~ /.+\.fasta$/);	

	($file) = ($line =~ /(.+\.fas?t?a?)/); 

        if (scalar(@wl) > 0) {
	    next if (!grep(/$file/, @wl));
	}

	push(@{$files}, $file);
    }
}

sub load_white_list {
    my ($wl) = @_;

    open(WHITE, "<" . $white_list) 
	or die("Unable to open white list file '$white_list': $!\n");

    my $line   = "";

    while ($line = <WHITE>) {
	chomp $line;

	next if (length($line) == 0);

	push(@{$wl}, $line);
    }

    close(WHITE);
}

sub parse_command_line {
    while (@ARGV) {
	$_ = shift @ARGV;
	if    ($_ =~ /^-s$/) { $single_path = shift @ARGV; }
	elsif ($_ =~ /^-p$/) { $paired_path = shift @ARGV; }
	elsif ($_ =~ /^-l$/) { $sanger_path = shift @ARGV; }
	elsif ($_ =~ /^-o$/) { $out_path    = shift @ARGV; }
	elsif ($_ =~ /^-W$/) { $white_list  = shift @ARGV; }
	elsif ($_ =~ /^-I$/) { $insert_len  = shift @ARGV; }
	elsif ($_ =~ /^-C$/) { $exp_cov     = shift @ARGV; }
	elsif ($_ =~ /^-T$/) { $cov_cut     = shift @ARGV; }
	elsif ($_ =~ /^-R$/) { $read_trk    = shift @ARGV; }
	elsif ($_ =~ /^-M$/) { $min_len     = shift @ARGV; }
	elsif ($_ =~ /^-H$/) { $hash_len    = shift @ARGV; }
	elsif ($_ =~ /^-D$/) { $ins_len_dist++; }
	elsif ($_ =~ /^-P$/) { $paired++; }
	elsif ($_ =~ /^-A$/) { $amos++; }
	elsif ($_ =~ /^-v$/) { version(); exit(); }
	elsif ($_ =~ /^-d$/) { $debug++; }
	elsif ($_ =~ /^-h$/) { usage(); }
	else {
	    usage();
	}
    }

    $single_path = substr($single_path,  0, -1) if (substr($single_path, -1)  eq "/");
    $paired_path = substr($paired_path,  0, -1) if (substr($paired_path, -1)  eq "/");
    $out_path    = substr($out_path, 0, -1) if (substr($out_path, -1) eq "/");

    #if ($paired && $insert_len == 0) {
#	print STDERR "You must specify the insert length when processing paired-end reads.\n";
#	usage();
#    }

#    if ($paired && $exp_cov == 0) {
#	print STDERR "You must specify the expected coverage when processing paired-end reads.\n";
#	usage();
#    }
}

sub version {
    print STDERR "exec-velvet.pl ", stacks_version, "\n";
}

sub usage {
    version();

    print STDERR <<EOQ; 
exec-velvet.pl -p path -s path [-l path] -o path [-W white_list] [-P] [-I len] [-C exp_cov] [-H len] [-T cov_cut] [-A] [-d] [-h]
  p: path to the paired-end FASTA files to assemble.
  s: path to the single-end FASTA files to assemble.
  l: path to long, sanger-style reads to assemble.
  o: path to output the assembled files.
  W: a white list of files to process in the input path.
  H: length of overlap required for reads (hash length, default 27).
  P: process paired-end reads.
  I: insert length (for paired-end reads; see velvet documentation).
  C: expected coverage (for paired-end reads).
  T: coverage cutoff.
  M: minimum contig length, discard contigs shorter than this value.
  R: turn on velvet's read tracking (uses additional memory).
  D: calculate the observed insert length (for paired-end reads). 
  A: generate ACE files for viewing in Eagleview.
  h: display this help message.
  d: turn on debug output.

EOQ

  exit(0);
}
