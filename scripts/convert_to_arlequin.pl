#!/usr/bin/perl
#
# Copyright 2014, Julian Catchen <jcatchen@uoregon.edu>
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

use constant true  => 1;
use constant false => 0;

my $in_path  = "";
my $out_path = "";
my $locus_id = 0;

parse_command_line();

my %haplotypes;
my %hapmap;
my %pops;

my ($fh, $line, @parts, @haps, $hap, $pop, $h, $c);

open($fh, "<$in_path") or die("Unable to open tags output file '$in_path', $!\n");

while ($line = <$fh>) {
    chomp $line;
    @parts = split(/\t/, $line);

    next if ($parts[1] != $locus_id);

    print STDERR $line, "\n";

    $pops{$parts[4]} = {'cnt' => $parts[5], 
			'hap' => {}};

    #
    # Parse haplotype string.
    #
    @haps = split(/;/, $parts[13]);

    foreach $hap (@haps) {
	($h, $c) = ($hap =~ /(\w+):([0-9]+)/);
	$pops{$parts[4]}->{'hap'}->{$h} = $c;
	$haplotypes{$h}++;
    }
}

close($fh);

#
# Create a map of haplotypes.
#
my $i = 1;
foreach $hap (keys %haplotypes) {
    $hapmap{$hap} = "h" . $i;
    $i++;
}

my $num_samples = scalar(keys %pops);
#
# Write out the Arlequin output.
#

print <<EOF;
[Profile]

  NbSamples=$num_samples
  DataType=DNA
  GenotypicData=0
  GameticPhase=0
  LocusSeparator=TAB
  RecessiveData=0
  MissingData='?'

[Data]

[[HaplotypeDefinition]]

  HaplList={
EOF

foreach $hap (keys %hapmap) {
    print "    ", $hapmap{$hap}, " ", $hap, "\n";
}

print <<EOF;
  }

[[Samples]]

EOF

foreach $pop (keys %pops) {
    print 
	"  SampleName=\"", $pop, "\"\n",
	"  SampleSize=", $pops{$pop}->{'cnt'}, "\n",
	"  SampleData= {\n";

    foreach $hap (keys %{$pops{$pop}->{'hap'}}) {
	print "    ", $hapmap{$hap}, "  ", $pops{$pop}->{'hap'}->{$hap}, "\n";
    }
    print "  }\n";
}

print <<EOF;

[[Structure]] 

StructureName="One group"
NbGroups=1
Group={
EOF

foreach $pop (keys %pops) {
    print "    \"", $pop, "\"\n",
}

print "}\n";


sub parse_command_line {
    while (@ARGV) {
	$_ = shift @ARGV;
	if    ($_ =~ /^-p$/) { $in_path  = shift @ARGV; }
	elsif ($_ =~ /^-i$/) { $locus_id = shift @ARGV; }
	elsif ($_ =~ /^-v$/) { version(); exit(); }
	elsif ($_ =~ /^-h$/) { usage(); }
	else {
	    print STDERR "Unknown command line option: '$_'\n";
	    usage();
	}
    }

    if (! -e $in_path) {
	print STDERR "Unable to locate Stacks haplotype statistics file: $in_path\n";
	usage();
    }

    if ($locus_id == 0) {
	print STDERR "You must specify a locus ID to convert.\n";
	usage();
    }
}

sub version {
    print STDERR "convert_to_arlequin.pl ", stacks_version, "\n";
}

sub usage {
    version();

    print STDERR <<EOQ; 
convert_to_arlequin.pl -p path -i id [-d] [-h]
    p: path to the Stacks batch_X.hapstats.tsv file.
    i: locus ID to convert.
    h: display this help message.

EOQ

exit(0);
}
