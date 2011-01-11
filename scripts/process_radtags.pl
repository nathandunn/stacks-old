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
# Written by Julian Catchen <jcatchen@uoregon.edu>
#

use strict;

use constant true        => 0;
use constant false       => 1;
use constant score_limit => 10; # Quality Score >= 90%
use constant len_limit   => 31;

use constant stacks_version => "_VERSION_";

my $debug        = 0;
my $quality      = 0;
my $input_type   = "raw";
my $output_type  = "fasta";
my $interleaved  = false;
my $clean        = 0;
my $paired       = 0;
my $truncate     = 0;
my $recover      = 0;
my $in_path      = ".";
my $out_path     = ".";
my $barcode_list = "";
my $barcode_size = 5;

my %encoded;
my $tag  = 'sbfI';
my %tags = ('sbfI'  => ['TGCAGG'],  # CCTGCA/GG, SbfI
	    'pstI'  => ['TGCAG'],   # CTGCA/G, PstI
	    'ecoRI' => ['AATTC'],   # G/AATTC, EcoRI
	    'sgrAI' => ['CCGGCG', 'CCGGTG']); # CR/CCGGYG, SgrAI; R=A or G; Y=C or T

my $input_types  = {'fastq' => \&parse_fastq_record,
                    'raw'   => \&parse_raw_record
                    };
my $output_types = {'fasta' => \&write_fasta,
                    'fastq' => \&write_fastq,
                    'raw'   => \&write_raw
                    };

parse_command_line();

my (@files, @barcodes, $num_files, $prefix, $seq_fh, $i);

build_file_list(\@files);
load_barcode_list(\@barcodes) if (length($barcode_list) > 0);

$num_files = scalar(@files);

my (%pair_1_fhs, %pair_2_fhs, %barcode_log, %counters);

open_files(\@barcodes, \%pair_1_fhs, \%pair_2_fhs, \%counters);

$i = 1;

foreach $prefix (@files) {

    printf(STDERR "Processing file % 3s of % 3s [%s]\n", $i, $num_files, $prefix);

    $counters{$prefix} = {};
    $counters{$prefix}->{'total'}       = 0;
    $counters{$prefix}->{'low_quality'} = 0;
    $counters{$prefix}->{'noradtag'}    = 0;
    $counters{$prefix}->{'ambiguous'}   = 0;
    $counters{$prefix}->{'retained'}    = 0;
    $counters{$prefix}->{'orphaned'}    = 0;
    $counters{$prefix}->{'recovered'}   = 0;

    if ($paired) {
	process_paired_reads($prefix, \%pair_1_fhs, \%pair_2_fhs, $counters{$prefix}, \%barcode_log);
    } else {
	process_reads($prefix, \%pair_1_fhs, $counters{$prefix}, \%barcode_log);
    }

    print STDERR 
	"  ", 
        $counters{$prefix}->{'total'}, " total reads; ", 
	"-", $counters{$prefix}->{'ambiguous'}, " ambiguous barcodes; ",
	"-", $counters{$prefix}->{'noradtag'}, " ambiguous RAD-Tags; ",
	"+", $counters{$prefix}->{'recovered'}, " recovered; ",
	"-", $counters{$prefix}->{'low_quality'}, " low quality reads; ",
	"-", $counters{$prefix}->{'orphaned'}, " orphaned paired-end reads; ",
	$counters{$prefix}->{'retained'}, " retained reads.\n";

    $i++;
}

print_results(\@barcodes, \%counters, \%barcode_log);

close_file_handles(\%pair_1_fhs);
close_file_handles(\%pair_2_fhs);

sub process_reads {
    my ($prefix, $pair_1_fhs, $counter, $barcode_log) = @_;

    my ($line, $seq_fh, $suffix, $href);

    $suffix = ($input_type eq "raw") ? "_qseq.txt" : "_sequence.txt";

    open($seq_fh, "<" . $in_path . "/" . $prefix . $suffix) 
        or die("Unable to open sequence file '$prefix': $!\n");

    while (defined($href = $input_types->{$input_type}->($seq_fh))) {
	$counter->{'total'}++;

	process_singlet($href, $barcode_log, $counter, false);

	if ($href->{'retain'}) {
	    $output_types->{$output_type}->($pair_1_fhs, $href);
	}
    }

    close($seq_fh);
}

sub process_paired_reads {
    my ($prefix, $pair_1_fhs, $pair_2_fhs, $counter, $barcode_log) = @_;

    my ($suffix, $prefix_1, $prefix_2, $seq_fh_1, $seq_fh_2, $href_1, $href_2, $line);

    $prefix_1 = $prefix;

    if ($input_type eq "raw") {
        $prefix   =~ /^s_(\d)_1_(\d+)/;
        $prefix_2 = "s_" . $1 . "_2_" . $2;
    } else {
        $prefix   =~ /^s_(\d)_1/;
        $prefix_2 = "s_" . $1 . "_2";
    }

    $suffix = $input_type eq "raw" ? "_qseq.txt" : "_sequence.txt";

    open($seq_fh_1, "<" . $in_path . "/" . $prefix_1 . $suffix) 
        or die("Unable to open sequence file '$prefix_1': $!\n");
    open($seq_fh_2, "<" . $in_path . "/" . $prefix_2 . $suffix) 
        or die("Unable to open sequence file '$prefix_2': $!\n");

    while (defined(($href_1 = $input_types->{$input_type}->($seq_fh_1)))) {
	$href_2 = $input_types->{$input_type}->($seq_fh_2);

	$counter->{'total'} += 2;

	process_singlet($href_1, $barcode_log, $counter, false);

	if (!$href_1->{'retain'}) {
	    $href_2->{'retain'} = 0;
	    $counter->{'orphaned'}++;
	} else {
	    $href_2->{'barcode'} = $href_1->{'barcode'};
	    process_singlet($href_2, $barcode_log, $counter, true);
	}

	if ($href_1->{'retain'} && $href_2->{'retain'}) {
	    $output_types->{$output_type}->($pair_1_fhs, $href_1);
            $output_type eq "fasta" ?
                $output_types->{$output_type}->($pair_1_fhs, $href_2) :
                $output_types->{$output_type}->($pair_2_fhs, $href_2);

	} elsif ($href_1->{'retain'} && !$href_2->{'retain'}) {
	    $output_types->{$output_type}->($pair_1_fhs, $href_1);
	}
    }

    close($seq_fh_1);
    close($seq_fh_2);
}

sub process_singlet {
    my ($href, $barcode_log, $counter, $paired_end) = @_;

    my ($rad_cor, $t);

    #
    # If requested, truncate this read to $truncate nucleotides
    #
    if ($truncate) {
	$href->{'seq'}   = substr($href->{'seq'},   0, $truncate);
	$href->{'phred'} = substr($href->{'phred'}, 0, $truncate);
    }

    if ($paired_end == false) {
	# Pull off the barcode
	$href->{'barcode'} = substr($href->{'seq'}, 0, $barcode_size);

	#
	# Shorten the sequence by the length of the barcode. But, unlike de novo 
	# sequencing, there is no nucleotide overhang to remove.
	#
	$href->{'seq'}   = substr($href->{'seq'},   $barcode_size);
	$href->{'phred'} = substr($href->{'phred'}, $barcode_size);

	#
	# Log the barcodes we receive.
	#
	if (!defined($barcode_log->{$href->{'barcode'}})) {
	    $barcode_log->{$href->{'barcode'}} = {};
	    $barcode_log->{$href->{'barcode'}}->{'noradtag'} = 0;
	    $barcode_log->{$href->{'barcode'}}->{'total'}    = 0;
            $barcode_log->{$href->{'barcode'}}->{'retained'} = 0;
	}
	$barcode_log->{$href->{'barcode'}}->{'total'}++;
    
	#
	# Is this a legitimate barcode?
	#
	if (!defined($pair_1_fhs{$href->{'barcode'}})) {
	    #
	    # Try to correct the barcode.
	    #
	    if (!correct_barcode($href, $counter, $barcode_log)) {
		$counter->{'ambiguous'}++;
		$href->{'retain'} = 0;
		return;
	    }
	}

	#
	# Is the RADTAG intact?
	#
        $rad_cor = 0;
        foreach $t (@{$tags{$tag}}) {
            $rad_cor++ if (substr($href->{'seq'}, 0, length($tags{$tag})) ne $tags{$tag});
        }
        if ($rad_cor == 0) {
	    #
	    # Try to correct the RAD-Tag.
	    #
	    if (!correct_radtag($href, $counter)) {
		$barcode_log->{$href->{'barcode'}}->{'noradtag'}++;
		$counter->{'noradtag'}++;
		$href->{'retain'} = 0;
		return;
	    }
	}
    }

    # Drop this sequence if it has any uncalled nucleotides
    if ($clean && $href->{'seq'} =~ /\.|N/) {
	$counter->{'low_quality'}++;
	$href->{'retain'} = 0;
	return;
    }

    # Drop this sequence if it has low quality scores
    if($quality && !check_quality_scores($href)) {
	$counter->{'low_quality'}++;
	$href->{'retain'} = 0;
	return;
    }

    $barcode_log->{$href->{'barcode'}}->{'retained'}++;
    $counter->{'retained'}++;
}

sub parse_raw_record {
    my ($seq_fh) = @_;

    my $line = <$seq_fh>;

    return undef if (!defined($line));

    my (@parts, $href);

    chomp $line;
    @parts = split(/\t/, $line);

    $href = {};
    $href->{'machine'} = $parts[0];
    $href->{'run'}     = $parts[1];
    $href->{'lane'}    = $parts[2];
    $href->{'tile'}    = $parts[3];
    $href->{'x'}       = $parts[4];
    $href->{'y'}       = $parts[5];
    $href->{'index'}   = $parts[6];
    $href->{'read'}    = $parts[7];
    $href->{'seq'}     = $parts[8];
    $href->{'phred'}   = $parts[9];
    $href->{'filter'}  = $parts[10];
    $href->{'retain'}  = 1;

    return $href;
}

sub parse_fastq_record {
    my ($seq_fh) = @_;

    my (@parts, $line, $href);

    $line = <$seq_fh>;
    return undef if (!defined($line));

    chomp $line;
    $line  = substr($line, 1);
    @parts = split(/:/, $line);

    #print STDERR "ID: $line\n";

    ($parts[5], $parts[6], $parts[7]) = 
        ($parts[4] =~ /(\d+)\#(\d+)\/(\d+)/);

    $href = {};
    $href->{'machine'} = $parts[0];
    $href->{'run'}     = "";
    $href->{'lane'}    = $parts[1];
    $href->{'tile'}    = $parts[2];
    $href->{'x'}       = $parts[3];
    $href->{'y'}       = $parts[5];
    $href->{'index'}   = $parts[6];
    $href->{'read'}    = $parts[7];
    $href->{'retain'}  = 1;

    $line = <$seq_fh>;
    return undef if (!defined($line));

    chomp $line;
    $href->{'seq'} = $line;

    $line = <$seq_fh>;
    return undef if (!defined($line));
    $line = <$seq_fh>;
    return undef if (!defined($line));

    chomp $line;
    $href->{'phred'} = $line;

    return $href;
}

sub correct_barcode {
    my ($href, $counter, $barcode_log) = @_;

    return 0 if (!$recover);

    #
    # If the barcode sequence is off by no more than a single nucleotide, correct it.
    #
    my ($dist, $close, $barcode, $b, $old_barcode);

    $close = 0;

    foreach $barcode (@barcodes) {
	$dist = dist($barcode, $href->{'barcode'}); 

	if ($dist <= 1) {
	    $close++;
	    $b = $barcode;
	}
    }

    if ($close == 1) {
	#
	# Correct the barcode.
	#
	$old_barcode = $href->{'barcode'};
	$href->{'barcode'} = $b;
	$counter->{'recovered'}++;
	$barcode_log->{$old_barcode}->{'total'}--;

        if (!defined($barcode_log->{$href->{'barcode'}})) {
            $barcode_log->{$href->{'barcode'}}->{'total'}    = 0;
            $barcode_log->{$href->{'barcode'}}->{'retained'} = 0;
            $barcode_log->{$href->{'barcode'}}->{'noradtag'} = 0;
        }
	$barcode_log->{$href->{'barcode'}}->{'total'}++;
	return 1;
    }

    return 0;
}

sub correct_radtag {
    my ($href, $counter) = @_;

    return 0 if (!$recover);

    my ($t, $dist);

    #
    # If the RAD-Tag sequence is off by no more than a single nucleotide, correct it.
    #
    foreach $t (@{$tags{$tag}}) {

        $dist = dist($t, substr($href->{'seq'}, 0, length($t)));

        if ($dist <= 1) {
            #
            # Correct the read.
            #
            $href->{'seq'} = $t . substr($href->{'seq'}, length($t));
            $counter->{'recovered'}++;

            return 1;
        }
    }

    return 0;
}

sub dist {
    #
    # Calculate the number of mismatches between two strinsg. This is done by 
    # XORing the strings togeher, which will give a string of 1s and 0s corresponding
    # to the match or mismatch of the charaters. The tr operator then counts the 
    # number of 0s in the string (matches) and we subtract that from the length to
    # get the number of mismatches.
    #
    return length($_[0]) - (($_[0] ^ $_[1]) =~ tr/\0/\0/)
}

sub check_quality_scores {
    my ($href) = @_;

    #
    # Phred quality scores are discussed here:
    #  http://en.wikipedia.org/wiki/FASTQ_format
    #
    # Illumina 1.3+ encodes phred scores between ASCII values 64 (0 quality) and 104 (40 quality)
    #
    #   @ABCDEFGHIJKLMNOPQRSTUVWXYZ[\]^_`abcdefgh
    #   |         |         |         |         |
    #  64        74        84        94       104
    #   0        10(90%)   20(99%)   30(99.9%) 40(99.99%)
    #
    my ($j, $k, @scores, @int_scores, $score, $len, $win_len, $stop_pos, $trim, $mean, $int);

    #
    # Trim sequences if phred quality score drops below a threshold.
    #
    # Window length is 15% (rounded) of the sequence length.
    #
    $len      = length($href->{'seq'});
    $win_len  = int($len * .15 + 0.5);
    $stop_pos = $len - $win_len - 1;
    $trim     = 0;
    @scores   = split(//, $href->{'phred'});

    #
    # Convert the encoded quality scores to their integer values
    #
    foreach $j (0 .. $len - 1) {
        $int_scores[$j] = ord($scores[$j]) - 64;
    }

    #print STDERR 
    #   "PHRED: $href->{'phred'}\n",
    #	"  Len: $len; Window len: $win_len; Stop pos; $stop_pos\n";

    foreach $j (0 .. $stop_pos) {
	#print STDERR "  J: $j\n";
	$mean = 0;

	foreach $k ($j .. $j + $win_len - 1) {
	    $mean += $int_scores[$k];
	    #print STDERR "    $score: ord ", ord($score), "; int $int\n";
	}

	$mean = $mean / $win_len;

	#print STDERR "Index $j; mean: $mean\n";

	if ($mean < score_limit) {
	    $trim = $j;
	}

	last if ($trim);
    }

    if ($trim) {
	$href->{'retain'} = 0;
	return 0;
    }

    return 1;
}

sub write_fasta {
    my ($fhs, $href) = @_;

    print {$fhs->{$href->{'barcode'}}} 
    ">", $href->{'barcode'}, "_", $href->{'lane'}, "_", sprintf("%04d", $href->{'tile'}), "_", $href->{'x'}, "_", $href->{'y'}, "_", $href->{'read'}, "\n",
    $href->{'seq'}, "\n";
}

sub write_fastq {
    my ($fhs, $href) = @_;

    #
    # Write the sequence and quality scores in FASTQ format. 
    #
   print {$fhs->{$href->{'barcode'}}}
   "\@", $href->{'barcode'}, "_", $href->{'lane'}, "_", sprintf("%04d", $href->{'tile'}), "_", $href->{'x'}, "_", $href->{'y'}, "_", $href->{'read'}, "\n",
   $href->{'seq'}, "\n",
   "+\n",
   $href->{'phred'}, "\n";
}

sub write_raw {
    my ($fhs, $href) = @_;

    print {$fhs->{$href->{'barcode'}}}
    $href->{'machine'}, "\t",
    $href->{'run'}, "\t",
    $href->{'lane'}, "\t",
    $href->{'tile'}, "\t",
    $href->{'x'}, "\t",
    $href->{'y'}, "\t",
    $href->{'index'}, "\t",
    $href->{'read'}, "\t",
    $href->{'seq'}, "\t",
    $href->{'phred'}, "\t",
    $href->{'filter'}, "\n";
}

sub print_results {
    my ($barcodes, $counters, $barcode_log) = @_;

    my ($log, $file, $counter, $barcode);

    $log = $out_path . "/process_radtags.log";
    open(LOG, ">$log") or die("Unable to open log file: '$log'\n");

    print LOG 
	"File\t",
	"Retained Reads\t",
	"Low Quality\t",
	"Ambiguous Barcodes\t",
	"Ambiguous RAD-Tag\t",
	"Orphaned paired-end reads\t",
	"Total\n";

    foreach $file (sort keys %{$counters}) {
	print LOG 
	    $file, "\t",
	    $counters->{$file}->{'retained'},    "\t",
	    $counters->{$file}->{'low_quality'}, "\t",
	    $counters->{$file}->{'ambiguous'},   "\t",
	    $counters->{$file}->{'noradtag'},    "\t",
	    $counters->{$file}->{'orphaned'},    "\t",
	    $counters->{$file}->{'total'},       "\n";
    }

    my $c = {};
    $c->{'total'}       = 0;
    $c->{'low_quality'} = 0;
    $c->{'ambiguous'}   = 0;
    $c->{'noradtag'}    = 0;
    $c->{'orphaned'}    = 0;

    #
    # Total up the individual counters
    #
    foreach $file (sort keys %{$counters}) {
	$c->{'total'}         += $counters->{$file}->{'total'};
	$c->{'low_quality'}   += $counters->{$file}->{'low_quality'};
	$c->{'ambiguous'}     += $counters->{$file}->{'ambiguous'};
	$c->{'orphaned'}      += $counters->{$file}->{'orphaned'};
	$c->{'noradtag'}      += $counters->{$file}->{'noradtag'};
	$c->{'retained'}      += $counters->{$file}->{'retained'}; 
    }

    print STDERR 
	"$c->{'total'} total sequences;\n",
	"  $c->{'ambiguous'} uncalled nucleotide drops;\n",
	"  $c->{'low_quality'} low quality read drops;\n", 
	"  $c->{'noradtag'} ambiguous RAD-Tag drops;\n",
	"  $c->{'orphaned'} orphaned paired-end reads;\n",
	"$c->{'retained'} retained reads.\n";

    print LOG 
	"Total Sequences\t",      $c->{'total'},       "\n",
	"Ambiguous Barcodes\t",   $c->{'ambiguous'},   "\n",
	"Low Quality\t",          $c->{'low_quality'}, "\n",
	"Ambiguous RAD-Tag\t",    $c->{'noradtag'},    "\n",
	"Orphaned Paired-ends\t", $c->{'orphaned'},    "\n",
	"Retained Reads\t",       $c->{'retained'},    "\n";

    #
    # Print out barcode information.
    #
    print LOG 
	"\n",
	"Barcode\t",
	"Total\t",
	"No RadTag\t",
	"Retained\n";

    foreach $barcode (@{$barcodes}) {
        if (!defined($barcode_log->{$barcode})) {
            print LOG
                $barcode, "\t", "0\t", "0\t", "0\n";
        } else {
            print LOG 
                $barcode, "\t",
                $barcode_log->{$barcode}->{'total'}, "\t",
                $barcode_log->{$barcode}->{'noradtag'}, "\t",
                $barcode_log->{$barcode}->{'retained'}, "\n";
        }
    }

    print LOG 
	"\n",
	"Sequences not recorded\n",
	"Barcode\t",
	"Total\n";

    foreach $barcode (sort {$barcode_log->{$b}->{'total'} <=> $barcode_log->{$a}->{'total'}} keys %{$barcode_log}) {
	next if (defined($pair_1_fhs{$barcode}));
	next if ($barcode_log->{$barcode}->{'total'} == 0);

	print LOG 
	    $barcode, "\t",
	    $barcode_log->{$barcode}->{'total'}, "\n";
    }

    close(LOG);
}

sub open_files {
    my ($barcodes, $pair_1_fhs, $pair_2_fhs, $counters) = @_;

    my ($barcode, $suffix_1, $suffix_2);

    if ($interleaved) {
        if ($output_type eq "fastq") {
            $suffix_1 = ".fq";
        } elsif ($output_type eq "raw") {
            $suffix_1 = ".raw";
        } else {
            $suffix_1 = ".fa";
        }
    } else {
        if ($output_type eq "fastq") {
            $suffix_1 = ".fq_1";
            $suffix_2 = ".fq_2";
        } elsif ($output_type eq "raw") {
            $suffix_1 = ".raw_1";
            $suffix_2 = ".raw_2";
        } else {
            $suffix_1 = ".fa_1";
            $suffix_2 = ".fa_2";
        }
    }

    foreach $barcode (@{$barcodes}) {

        if ($interleaved == true) {
            open($pair_1_fhs->{$barcode}, ">" . $out_path . "/sample_" . $barcode . $suffix_1) 
                or die("Unable to open FASTA file for barcode '$barcode': $!\n");

        } else {
            open($pair_1_fhs->{$barcode}, ">" . $out_path . "/sample_" . $barcode . $suffix_1) 
                or die("Unable to open FASTA file for barcode '$barcode': $!\n");

            if ($paired) {
                open($pair_2_fhs->{$barcode}, ">" . $out_path . "/sample_" . $barcode . $suffix_2) 
                    or die("Unable to open quality score file for barcode '$barcode': $!\n");
            }
        }
    }
}

sub load_barcode_list {
    my ($bl) = @_;

    open(BC, "<" . $barcode_list) 
	or die("Unable to open barcode file '$barcode_list': $!\n");

    my ($line, $bc, $blen, $prev);

    while ($line = <BC>) {
	chomp $line;

	next if (length($line) == 0);

	($bc) = ($line =~ /^([ACGT]+)$/);

	push(@{$bl}, $bc);
    }

    close(BC);

    if (scalar(@{$bl}) == 0) {
	print STDERR "Unable to load any barcodes from '$barcode_list'\n";
	usage();
    }

    #
    # Determine the barcode length
    #
    $prev = length($bl->[0]);
    foreach $bc (@{$bl}) {
        $blen = length($bc);

        if ($prev != $blen) {
            print STDERR "Barcodes must all be the same length. Place different barcode lengths in separate runs.\n";
            usage();
        }
    }

    $barcode_size = $blen;
}

sub build_file_list {
    my ($files) = @_;

    my (@ls, $line, $prefix);

    @ls = $input_type eq "raw" ? 
        `ls -1 $in_path/s_*_1_*_qseq.txt 2> /dev/null` :
        `ls -1 $in_path/s_*_sequence.txt 2> /dev/null`;

    if (scalar(@ls) == 0) {
	print STDERR "Unable to locate any input files to process within '$in_path'\n";
	usage();
    }

    foreach $line (@ls) {
	chomp $line;

        if ($input_type eq "raw") {
            ($prefix) = ($line =~ /$in_path\/(s_\d_\d_\d{4})_qseq\.txt/);
        } else {
            ($prefix) = ($line =~ /$in_path\/(s_\d_?1?)_sequence\.txt/);
        }

	push(@{$files}, $prefix);
    }
}

sub close_file_handles {
    my ($fhs) = @_;

    foreach my $barcode (keys %{$fhs}) {
	close($fhs->{$barcode});
    }
}

sub parse_command_line {
    while (@ARGV) {
	$_ = shift @ARGV;
	if    ($_ =~ /^-p$/) { $in_path  = shift @ARGV; }
	elsif ($_ =~ /^-o$/) { $out_path = shift @ARGV; }
	elsif ($_ =~ /^-b$/) { $barcode_list = shift @ARGV; }
	elsif ($_ =~ /^-e$/) { $tag = shift @ARGV; }
	elsif ($_ =~ /^-t$/) { $truncate = shift @ARGV; }
	elsif ($_ =~ /^-F$/) { $output_type = "fastq"; }
	elsif ($_ =~ /^-R$/) { $output_type = "raw"; }
        elsif ($_ =~ /^-I$/) { $input_type = shift @ARGV; }
	elsif ($_ =~ /^-c$/) { $clean++; }
	elsif ($_ =~ /^-r$/) { $recover++; }
	elsif ($_ =~ /^-q$/) { $quality++; }
	elsif ($_ =~ /^-P$/) { $paired++; }
	elsif ($_ =~ /^-d$/) { $debug++; }
	elsif ($_ =~ /^-v$/) { version(); exit(); }
	elsif ($_ =~ /^-h$/) { usage(); }
	else {
	    usage();
	}
    }

    $in_path  = substr($in_path, 0, -1)  if (substr($in_path, -1)  eq "/");
    $out_path = substr($out_path, 0, -1) if (substr($out_path, -1) eq "/");

    if (length($barcode_list) == 0) {
	print STDERR "You must specify a list of barcodes.\n";
	usage();
    }

    if ($input_type ne "raw" && $input_type ne "fastq") {
        print STDERR "Input type must be either 'raw' or 'fastq'.\n";
        usage();
    }
}

sub version {
    print STDERR "process_radtags.pl ", stacks_version, "\n";
}

sub usage {
    version();

    print STDERR <<EOQ; 
process_radtags.pl -p path -o path -b file [-P] [-I input_type] [-e enzyme] [-c] [-q] [-r] [-t len] [-F|-R] [-d] [-h]
  p: path to the Solexa BUSTARD output files, or GERALD files (if -I is specified as 'fastq').
  o: path to output the processed files.
  P: input data contains paired-end reads.
  b: a list of barcodes for this run.
  I: input file type, either 'raw' for the raw BUSTARD output files, or 'fastq' for GERALD files (default 'raw').
  c: clean data, remove any read with an uncalled base.
  q: discard reads with low quality scores.
  r: rescue barcodes and RAD-Tags.
  e: specify the restriction enzyme to look for (either 'sbfI', 'pstI', 'ecoRI', or 'sgrAI').
  t: truncate final read length to this value.
  F: output a FASTQ file instead of a FASTA file.
  R: output data in the Raw BUSTARD format.
  h: display this help message.
  d: turn on debug output.

EOQ

  exit(0);
}
