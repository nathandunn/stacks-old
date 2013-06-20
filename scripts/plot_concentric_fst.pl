#!/usr/bin/perl
#
# Copyright 2013, Julian Catchen <jcatchen@uoregon.edu>
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
use GD;
use GD::Polyline;
use POSIX qw(floor ceil);
use PDF::API2;

use constant stacks_version => "_VERSION_";
use constant PI  => 4 * atan2(1, 1);
use constant in  => 1/72;
use constant px  => 1;
use constant pt  => 1;
use constant deg => 1;
use constant bp  => 1;
use constant truecolor => 1;
use constant true  => 1;
use constant false => 0;

my $debug     = 0;
my $org       = "gac";
my $out_path  = "";

#
# defines physical characteristics for each genome,
# including number of chromosomes and their names.
#
my %orgs = ('gac' => {'groupI'    => 28185914, 'groupII'   => 23295652, 'groupIII'   => 16798506, 'groupIV'  => 32632948,
		      'groupIX'   => 20249479, 'groupV'    => 12251397, 'groupVI'    => 17083675, 'groupVII' => 27937443,
		      'groupVIII' => 19368704, 'groupX'    => 15657440, 'groupXI'    => 16706052, 'groupXII' => 18401067,
		      'groupXIII' => 20083130, 'groupXIV'  => 15246461, 'groupXIX'   => 20240660, 'groupXV'  => 16198764,
		      'groupXVI'  => 18115788, 'groupXVII' => 14603141, 'groupXVIII' => 16282716, 'groupXX'  => 19732071,
		      'groupXXI'  => 11717487},
	    'mmu' => {'Y' => 15902555,   '19' => 61342430,  '18' => 90772031,  '17' => 95272651, '16' => 98319150,
		      '15' => 103494974, '13' => 120284312, '12' => 121257530, '11' => 121843856, '9' => 124076172,
		      '14' => 125194864, '10' => 129993255,  '8' => 131738871,  '6' => 149517037, '7' => 152524553,
		      '5' => 152537259,   '4' => 155630120,  '3' => 159599783,  'X' => 166650296, '2' => 181748087,
		      '1' => 197195432});

#
# Define the path to the font we plan to use in the image
#
my $font = "/usr/share/fonts/liberation/LiberationSans-Regular.ttf";

my (@chrs, @fst_files, %img, %fst);

$img{'type'}      = "pdf";
$img{'width'}     = 800/px;
$img{'min_dist'}  = 1000000/bp;
$img{'font_size'} = 16;

parse_command_line(\%img);

populate_fst(\@fst_files, \%fst);

if ($img{'type'} eq "gd") {
    generate_png(\@chrs, \%fst, \%img);

} else {
    generate_pdf(\@chrs, \%fst, \%img);
}

sub populate_fst {
    my ($fst_files, $fst) = @_;

    my ($fh, $line, @parts);

    foreach my $file (@{$fst_files}) {

	open($fh, "<$file") or die("Unable to open fst file '$file' ($!)\n");

	$fst->{$file} = {};

	while ($line = <$fh>) {
	    chomp $line;

	    next if (length($line) == 0 || $line =~ /^\s*#/);

	    @parts = split(/\t/, $line);

	    next if (!defined($orgs{$org}->{$parts[4]}));

	    push(@{$fst->{$file}->{$parts[4]}->{'bp'}},  $parts[5]);
	    push(@{$fst->{$file}->{$parts[4]}->{'fst'}}, $parts[19]);
	}

	close($fh);
    }
}

sub generate_pdf {
    my ($chrs, $fst, $img) = @_;

    my (%colors, $chr, $key);

    allocate_pdf_colors(\%colors);

    $img->{'width'}      = 800/px;
    $img->{'height'}     = 800/px;
    $img->{'margin_x'}   = 0/px;
    $img->{'margin_y'}   = 0/px;
    $img->{'arc_margin'} = 2/deg;

    #
    # Open the file and set the font
    #
    $img->{'pdf'}      = PDF::API2->new('-file' => $out_path);
    $img->{'pdf_font'} = $img->{'pdf'}->ttfont($font);
    
    # Get a new page
    $img->{'pdf_page'} = $img->{'pdf'}->page();

    # Set the size of the page
    $img->{'pdf_page'}->mediabox($img->{'width'}, $img->{'height'});

    # Determine the center of the image
    $img->{'center_x'}      = floor(($img->{'width'}  + $img->{'margin_x'}) / 2);
    $img->{'center_y'}      = floor(($img->{'height'} + $img->{'margin_y'}) / 2);
    $img->{'circle_width'}  = floor($img->{'width'}   * 0.8);
    $img->{'circle_height'} = floor($img->{'height'}  * 0.8);
    $img->{'min_radius'}    = floor($img->{'width'}   / 2 * 0.4);
    $img->{'max_radius'}    = floor($img->{'width'}   / 2 * 0.8);

    # Determine the length of each chromosome in degrees.
    determine_chromosome_lengths($img, $chrs);

    #
    # Determine the number of concentric circles to draw, and radial spacing for each.
    #
    my @fst_keys    = keys %{$fst};
    my $num_circles = scalar(@fst_keys); 
    my $radius_dist = ($img->{'max_radius'} - $img->{'min_radius'}) / ($num_circles - 1);

    #
    # Draw scale bars across the circles.
    #
    draw_scale_pdf($img, \%colors, \@chrs);

    #
    # Draw key for Fst coloration.
    #
    draw_key_pdf($img, \%colors);

    #
    # Draw the individual chromosomes as arcs along a circle
    #
    my $radius = $img->{'min_radius'};
    my $labels = false;
    my $fst_key;

    for $key (1..$num_circles) {
	$labels = $key == $num_circles ? true : false;
	$fst_key = shift @fst_keys;

	foreach $chr (@chrs) {
	    draw_chromosome_arc_pdf($img, \%colors, $chr, $radius, $labels);
	    draw_fst_pdf($img, $fst->{$fst_key}, $chr, $radius);
	}

	$radius += $radius_dist;
    }

    $img{'pdf'}->save();
}

sub draw_fst_pdf {
    my ($img, $fsts, $chr, $radius) = @_;

    my ($fst_cnt, $line, $bp, $fst, $start_point, $end_point, $x1, $y1, $x2, $y2, $c, $chr_len, $tics, $r, $g, $b, $color);

    my $min_rad = $radius - 6;
    my $max_rad = $radius + 6;

    $c       = $img->{'chrs'}->{$chr};
    $chr_len = $orgs{$org}->{$chr};
    $fst_cnt = scalar(@{$fsts->{$chr}->{'bp'}});

    foreach my $i (0..$fst_cnt - 1) {
	$bp  = $fsts->{$chr}->{'bp'}->[$i];
	$fst = $fsts->{$chr}->{'fst'}->[$i];

	($r, $g, $b) = hsv2rgb(scale_color($fst), 0.99, 0.99);
	$color       = color2hex($r, $g, $b);

	# print STDERR "Fst: $fst has color: $color\n";

	#
	# Scale the starting point for proper placement on the chromosome arc
	#
	$start_point = ($bp / $chr_len) * ($c->{'end_deg'} - $c->{'start_deg'});
	$start_point = $c->{'start_deg'} + $start_point;

	#
	# Convert the point to an x,y coordinate
	#
	$x1 = $min_rad * cos(deg2rad($start_point)) + $img->{'center_x'};
	$y1 = $min_rad * sin(deg2rad($start_point)) + $img->{'center_y'};

	$x2 = $max_rad * cos(deg2rad($start_point)) + $img->{'center_x'};
	$y2 = $max_rad * sin(deg2rad($start_point)) + $img->{'center_y'};

	$line = $img->{'pdf_page'}->gfx();

	$line->strokecolor($color);
	$line->linewidth(1/pt);

	$line->move($x1, $y1);
	$line->line($x2, $y2);
	$line->stroke();
    }
}

sub draw_key_pdf {
    my ($img, $colors) = @_;

    my ($x1, $y, $y1, $x2, $y2, $i, $r, $g, $b, $color);

    my $key = $img->{'pdf_page'}->gfx();

    $x1 = $img->{'width'}  - 50;
    $y1 = $img->{'height'} - 225;
    $x2 = $img->{'width'}  - 25;
    $y2 = $img->{'height'} - 25;

    $key->rectxy($x1, $y1, $x2, $y2);
    $key->strokecolor($colors->{'black'});
    $key->linewidth(2/pt);
    $key->stroke();
    $y = $y1;

    for ($i = 0; $i <= 1; $i += .005) {
	($r, $g, $b) = hsv2rgb(scale_color($i), 0.99, 0.99);
	$color       = color2hex($r, $g, $b);
	
	$key->strokecolor($color);
	$key->linewidth(1/pt);

	$key->move($x1, $y);
	$key->line($x2, $y);
	$key->stroke();

	$y++;
    }

    #
    # Draw the key title and scale bars
    #
    my $label = $img->{'pdf_page'}->text();
    $label->font($img->{'pdf_font'}, $img->{'font_size'} - 2);
    $label->fillcolor($colors->{'black'});

    # Get the size of the text we want to write
    my $text = "Fst";
    my $text_width  = $label->advancewidth($text);
    my $text_height = $img->{'font_size'} - 2;

    $label->translate(int(($x1 + $x2) / 2 + 0.5) - int($text_width / 2 + 0.5), $y2 + 5);
    $label->text($text);

    $label->font($img->{'pdf_font'}, $img->{'font_size'} - 6);
    $text_height = $img->{'font_size'} - 6;

    $text = "0.0";
    $text_width  = $label->advancewidth($text);
    $label->translate($x1 - $text_width - 3, $y1 - int($text_height / 2 + 0.5));
    $label->text($text);

    $text = "0.5";
    $text_width  = $label->advancewidth($text);
    $label->translate($x1 - $text_width - 3, int(($y1 + $y2) / 2 + 0.5) - int($text_height / 2 + 0.5));
    $label->text($text);

    $text = "1.0";
    $text_width  = $label->advancewidth($text);
    $label->translate($x1 - $text_width - 3, $y2 - int($text_height / 2 + 0.5));
    $label->text($text);
}

sub draw_scale_pdf {
    my ($img, $colors, $chrs) = @_;

    my ($chr, $line, $start_bp, $start_point, $end_point, $x1, $y1, $x2, $y2, $c, $chr_len, $tics, $text);

    my $min_rad = $img->{'min_radius'} - floor($img->{'width'} / 2 * 0.05);
    my $max_rad = $img->{'max_radius'} + floor($img->{'width'} / 2 * 0.05);

    foreach $chr (@{$chrs}) {
	$c       = $img->{'chrs'}->{$chr};
	$chr_len = $orgs{$org}->{$chr};

	#
	# Determine chromosome length rounded to nearest 10Mb.
	#
	$tics = int($chr_len / 5000000) % 5000000;

	foreach my $i (1..$tics) {
	    $start_bp = 5000000 * $i;

	    #
	    # Scale the starting point for proper placement on the chromosome arc
	    #
	    $start_point = ($start_bp / $chr_len) * ($c->{'end_deg'} - $c->{'start_deg'});
	    $start_point = $c->{'start_deg'} + $start_point;

	    #
	    # Convert the point to an x,y coordinate
	    #
	    $x1 = $min_rad * cos(deg2rad($start_point)) + $img->{'center_x'};
	    $y1 = $min_rad * sin(deg2rad($start_point)) + $img->{'center_y'};

	    $x2 = $max_rad * cos(deg2rad($start_point)) + $img->{'center_x'};
	    $y2 = $max_rad * sin(deg2rad($start_point)) + $img->{'center_y'};

	    $line = $img->{'pdf_page'}->gfx();
	    $line->linedash(6, 3);

	    if ($start_bp % 10000000 == 0) {
		$line->strokecolor($colors->{'black'});
		$line->linewidth(2/pt);
	    } else {
		$line->strokecolor($colors->{'grey'});
		$line->linewidth(1/pt);
	    }

	    $line->move($x1, $y1);
	    $line->line($x2, $y2);
	    $line->stroke();
 
	    #
	    # Draw a label every 10Mb on the interior of the circle.
	    #
	    if ($start_bp % 10000000 == 0) {
		my $text_loc  = deg2rad($start_point);

		$x1 = $min_rad * cos($text_loc) + $img->{'center_x'};
		$y1 = $min_rad * sin($text_loc) + $img->{'center_y'};

		my $label = $img->{'pdf_page'}->text();
		$label->font($img->{'pdf_font'}, $img->{'font_size'} - 6);
		$label->fillcolor($colors->{'black'});

		# Get the size of the text we want to write
		$text = int($start_bp / 1000000) . "Mb";
		my $text_width  = $label->advancewidth($text);
		my $text_height = $img->{'font_size'} - 6;

		# Adjust it, depending on the quadrant of the circle we are in
		my ($w, $h) = 0;

		if ($text_loc >= 0 && $text_loc <= 1.57) {
		    # Top right quadrant
		    $w -= ($text_width / 2 + 0.5);
		    $h -= $text_height;
		} elsif ($text_loc > 1.57 && $text_loc <= 3.14) {
		    # Top left quadrant
		    $w -= ($text_width / 2 + 0.5);
		    $h -= $text_height;
		} elsif ($text_loc > 3.14 && $text_loc <= 4.71) {
		    # Bottom left quadrant
		    $w -= ($text_width / 2 + 0.5);
		    $h += ($text_height / 2 + 0.5);
		} elsif ($text_loc > 4.71) {
		    # Bottom right quadrant
		    $w -= ($text_width / 2 + 0.5);
		    $h += int($text_height / 2 + 0.5);
		}

		$label->translate($x1+$w, $y1+$h);
		$label->text($text);
	    }
	}
    }

    $line->linedash();
}

sub draw_chromosome_arc_pdf {
    my ($img, $colors, $chr, $radius, $labels) = @_;

    my ($arc_len, $start_deg, $end_deg, $avg_arc) = 0;
    my ($text_width, $text_height, $text_loc, $text, $label, $arc, $x, $y);

    $start_deg = $img->{'chrs'}->{$chr}->{'start_deg'};
    $end_deg   = $img->{'chrs'}->{$chr}->{'end_deg'};

    $arc = $img->{'pdf_page'}->gfx();
    $arc->fillcolor($colors->{'grey'});
    $arc->strokecolor($colors->{'grey'});
    $arc->linewidth(12/pt);

    $arc->arc($img->{'center_x'}, $img->{'center_y'}, 
	      $radius, $radius, $start_deg, $end_deg, 1);
    $arc->stroke();

    #
    # Draw the chromosome name
    #
    if ($labels == true) {
	$text      = $chr;
	$text_loc  = deg2rad(int((($end_deg + $start_deg) / 2) + 0.5));

	$x = int(($img->{'max_radius'} + 20) * cos($text_loc) + $img->{'center_x'} + 0.5);
	$y = int(($img->{'max_radius'} + 20) * sin($text_loc) + $img->{'center_y'} + 0.5);

	$label = $img->{'pdf_page'}->text();
	$label->font($img->{'pdf_font'}, $img->{'font_size'});
	$label->fillcolor($colors->{'black'});

	# Get the size of the text we want to write
	$text_width  = $label->advancewidth($text);
	$text_height = $img->{'font_size'};

	# Adjust it, depending on the quadrant of the circle we are in
	my ($w, $h) = 0;

	if ($text_loc >= 0 && $text_loc <= 1.57) {
	    # Top right quadrant
	    $w  = 0;
	    $h -= int($text_height / 2 + 0.5);
	} elsif ($text_loc > 1.57 && $text_loc <= 3.14) {
	    # Top left quadrant
	    $w -= $text_width;
	    $h -= int($text_height / 2 + 0.5);
	} elsif ($text_loc > 3.14 && $text_loc <= 4.71) {
	    # Bottom left quadrant
	    $w -= $text_width;
	    $h -= int($text_height / 2 + 0.5);
	} elsif ($text_loc > 4.71) {
	    # Bottom right quadrant
	    $w  = 0;
	    $h -= int($text_height / 2 + 0.5);
	}

	$label->translate($x+$w, $y+$h);
	$label->text($text);
    }
    print STDERR
	"Drawing an arc centered at ($img->{'center_x'}, $img->{'center_y'}) with radius $img->{'radius'};\n",
	"    starting degree: $start_deg, end degree: $end_deg; length: ", 
	$end_deg - $start_deg, "; arc_len: $arc_len\n",
	"    text_loc: $text_loc, X: $x, Y: $y\n" if ($debug);
}

sub generate_png {
    my ($chrs, $fst, $img) = @_;

    my (%colors, $chr, $href);

    allocate_png_colors(\%colors);

    $img->{'width'}      = 800/px;
    $img->{'height'}     = 800/px;
    $img->{'margin_x'}   = 100/px;
    $img->{'margin_y'}   = 0/px;
    $img->{'arc_margin'} = 5/deg;

    # Create a new image.
    $img->{'gd'} = new GD::Image($img->{'width'}  + $img->{'margin_x'}, 
				 $img->{'height'} + $img->{'margin_y'}, truecolor);

    # Determine the center of the image.
    $img->{'center_x'}      = floor(($img->{'width'}  + $img->{'margin_x'}) / 2);
    $img->{'center_y'}      = floor(($img->{'height'} + $img->{'margin_y'}) / 2);
    $img->{'circle_width'}  = floor($img->{'width'}   * 0.8);
    $img->{'circle_height'} = floor($img->{'height'}  * 0.8);
    $img->{'radius'}        = floor($img->{'width'}   / 2 * 0.8);

    # Fill the background in white.
    $img->{'gd'}->filledRectangle(0, 0, 
				  $img->{'width'} + $img->{'margin_x'} - 1, 
				  $img->{'width'} + $img->{'margin_y'} - 1, 
				  $img->{'gd'}->colorAllocate(255, 255, 255));

    # Determine the length of each chromosome in degrees.
    determine_chromosome_lengths($img, \@chrs);

    # Draw the individual chromosomes as arcs along a circle
    foreach $chr (@chrs) {
	draw_chromosome_arc($img, \%colors, $chr);
    }

    open(OUT, ">$out_path") or die("Unable to open GD file '$out_path' ($!)\n");

    print STDERR "Output file: $out_path\n";

    # make sure we are writing to a binary stream
    binmode OUT;

    # Convert the image to PNG and print it
    print OUT $img->{'gd'}->png;

    close(OUT);
}

sub determine_chromosome_lengths {
    my ($img, $chrs) = @_;

    my ($total_chr_len, $total_arc_len, $chr, $num_chrs);

    # Determine the total length of the chromosomes
    $num_chrs      = scalar(@{$chrs});
    $total_chr_len = 0;

    foreach $chr (@{$chrs}) {
	$total_chr_len += $orgs{$org}->{$chr};
    }

    #
    # Determine the total size of the arc (in degrees) after 
    # removing space for margins around each chromosome.
    #
    $total_arc_len = (360 - ($num_chrs * $img->{'arc_margin'}));

    my ($arc_len, $start_deg, $end_deg, $degrees);

    #
    # Calculate arc length with each arc's size relative to the 
    # length of the other chromosomes. Round the answer and 
    # ensure no arcs are less than a single degree 
    #
    $degrees = 0;

    foreach $chr (@{$chrs}) {

	$arc_len   = int(($total_arc_len * ($orgs{$org}->{$chr} / $total_chr_len)) + 0.5);
	$arc_len   = $arc_len < 1 ? 1 : $arc_len;

	$start_deg = $degrees   + $img->{'arc_margin'};
	$end_deg   = $start_deg + $arc_len;

	$img->{'chrs'}->{$chr}->{'start_deg'} = $start_deg;
	$img->{'chrs'}->{$chr}->{'end_deg'}   = $end_deg;

	$degrees = $end_deg;
    }

    # Record the average arc length for later use in draw_chromosome_arc().
    $img->{'avg_arc'} = int($total_arc_len / $num_chrs + 0.5);
}

sub draw_chromosome_arc {
    my ($img, $colors, $chr) = @_;

    my ($arc_len, $start_deg, $mid_deg, $end_deg, $avg_arc) = 0;
    my (@bounds, $text, $text_loc, $x, $y, $polyline);

    $start_deg = $img->{'chrs'}->{$chr}->{'start_deg'};
    $end_deg   = $img->{'chrs'}->{$chr}->{'end_deg'};
    $arc_len   = $end_deg - $start_deg;

    # Create a new polyline
    $polyline = new GD::Polyline;

    $img->{'gd'}->setThickness(7);

    $x = int($img->{'radius'} * cos(deg2rad($start_deg)) + $img->{'center_x'} + 0.5);
    $y = int($img->{'radius'} * sin(deg2rad($start_deg)) + $img->{'center_y'} + 0.5);
    $polyline->addPt($x, $y);

    $mid_deg = $start_deg;

    #
    # A spline must have 3n+1 control points. We will use 16 control points when
    # the arcs are large, 7 control points when they are small to make sure we
    # get a smooth curve, and 3 control points when an arc is tiny. If we use
    # too many control points, they overlap, crashing the polyline module.
    #
    if ($arc_len < 10) {
	$img->{'control_pts'} = 3;
    } elsif ($arc_len < 45) {
	$img->{'control_pts'} = 7;
    } else {
	$img->{'control_pts'} = 16;
    }

    foreach (1 .. $img->{'control_pts'} - 2) {
	$mid_deg += int(($end_deg - $start_deg) / ($img->{'control_pts'} - 1) + 0.5);
	$x = int($img->{'radius'} * cos(deg2rad($mid_deg)) + $img->{'center_x'} + 0.5);
	$y = int($img->{'radius'} * sin(deg2rad($mid_deg)) + $img->{'center_y'} + 0.5);
	$polyline->addPt($x, $y);
    }

    $x = int($img->{'radius'} * cos(deg2rad($end_deg)) + $img->{'center_x'} + 0.5);
    $y = int($img->{'radius'} * sin(deg2rad($end_deg)) + $img->{'center_y'} + 0.5);
    $polyline->addPt($x, $y);

    $img->{'gd'}->polydraw($polyline->addControlPoints()->toSpline(), 
			   $colors->{'grey'});

    $img->{'gd'}->setThickness(1);

    #
    # Draw the chromosome name
    #
    $text      = $chr;
    $text_loc  = deg2rad(int((($end_deg + $start_deg) / 2) + 0.5));

    $x = int(($img->{'radius'} + 30) * cos($text_loc) + $img->{'center_x'} + 0.5);
    $y = int(($img->{'radius'} + 30) * sin($text_loc) + $img->{'center_y'} + 0.5);

    # Get the size of the text we want to write
    @bounds = GD::Image->stringFT($colors->{'green'}, $font, 16, 0, $x, $y, $text);
	
    # Adjust it, depending on the quadrant of the circle we are in
    my ($w, $h) = 0;
    if ($text_loc >= 0 && $text_loc <= 1.57) {
	$w  = 0;
	$h -= floor(($bounds[7] - $bounds[1]) / 2);
    } elsif ($text_loc > 1.57 && $text_loc <= 3.14) {
	$w -= $bounds[2] - $bounds[0];
	$h -= floor(($bounds[7] - $bounds[1]) / 2);
    } elsif ($text_loc > 3.14 && $text_loc <= 4.71) {
	$w -= $bounds[2] - $bounds[0];
	$h -= floor(($bounds[7] - $bounds[1]) / 2);
    } elsif ($text_loc > 4.71) {
	$w  = 0;
	$h -= floor(($bounds[7] - $bounds[1]) / 2);
    }

    $img->{'gd'}->stringFT($colors->{'black'}, $font, 16, 0, $x+$w, $y+$h, $text);

    print STDERR
	"Drawing an arc centered at ($img->{'center_x'}, $img->{'center_y'}) with radius $img->{'radius'};\n",
	"    starting degree: $start_deg, end degree: $end_deg; length: ", 
	$end_deg - $start_deg, "; arc_len: $arc_len\n",
	"    text_loc: $text_loc, X: $x, Y: $y\n" if ($debug);
}

sub deg2rad { 
    return PI * $_[0] / 180;
}

sub scale_color {
    my ($val) = @_;

    my $min_range = 1;
    my $max_range = 10;

    my $min_color = 0.0;
    my $max_color = 0.68;

    return $max_color if ($val == 0);

    $val = ($max_range - $min_range) * $val + 1;

    my $color = (log($val)/log(10)) * ($max_color - $min_color);
    #my $color = $val * ($max_color - $min_color);

    return $max_color - $color;

    #return log($val)/log(10);

    #return $val;
}

sub hsv2rgb {
    my ($H, $s, $v) = @_;

    #
    # Our value of $h is a fraction between 0 and 1, scale it to a number between 0 and 360.
    #
    my $h = int($H * 360 + 0.5);

    #
    # Conversion routine from Perl Monks: http://www.perlmonks.org?node_id=139485 
    #
    if ($s == 0) {
        return $v, $v, $v;
    }

    $h /= 60;
    my $i = floor($h);
    my $f = $h - $i;
    my $p = $v * (1 - $s);
    my $q = $v * (1 - ($s * $f));
    my $t = $v * (1 - ($s * (1 - $f)));

    if ( $i == 0 ) {
        return $v, $t, $p;
    }
    elsif ( $i == 1 ) {
        return $q, $v, $p;
    }
    elsif ( $i == 2 ) {
        return $p, $v, $t;
    }
    elsif ( $i == 3 ) {
        return $p, $q, $v;
    }
    elsif ( $i == 4 ) {
        return $t, $p, $v;
    }
    else {
        return $v, $p, $q;
    }
}

sub color2hex {
    my ($r, $g, $b) = @_;

    $r = int($r * 255 + 0.5);
    $g = int($g * 255 + 0.5);
    $b = int($b * 255 + 0.5);

    my $str = sprintf("#%02x%02x%02x", $r, $g, $b);

    return $str;
}

sub allocate_png_colors {
    my ($gd_img, $colors) = @_;

    $colors->{'black'} = $gd_img->colorAllocate(  0,   0,   0);
    $colors->{'white'} = $gd_img->colorAllocate(255, 255, 255);
    $colors->{'grey'}  = $gd_img->colorAllocate(127, 127, 127); #7f7f7f;
}

sub allocate_pdf_colors {
    my ($colors) = @_;

    $colors->{'black'}  = "#000000";
    $colors->{'white'}  = "#ffffff";
    $colors->{'grey'}   = "#7f7f7f";
}

sub parse_command_line {
    my ($img) = @_;

    my $chr_list;

    while (@ARGV) {
	$_ = shift @ARGV;
	if ($_ =~ /^-d$/)     { $debug++; }
	elsif ($_ =~ /^-o$/)  { $out_path  = shift @ARGV; }
	elsif ($_ =~ /^-O$/)  { $org       = shift @ARGV; }
	elsif ($_ =~ /^-c$/)  { $chr_list  = shift @ARGV; }
	elsif ($_ =~ /^-f$/)  { push(@fst_files, shift @ARGV); }
	elsif ($_ =~ /^-t$/)  { $img->{'type'}  = shift @ARGV; }
	elsif ($_ =~ /^-s$/)  { $img->{'width'} = shift @ARGV; }
	elsif ($_ =~ /^-h$/)  { usage(); }
	else {
	    print STDERR "Unknown command line options received: $_\n";
	    usage();
	}
    }

    if (length($out_path) == 0) {
	print STDERR "You must specify a path to write files to.\n";
	usage();
    }

    if ($img{'type'} ne "gd" && $img{'type'} ne "pdf") {
	print STDERR "Unknown image type specified: $img{'type'}\n";
	usage();
    }

    #
    # Generate a list of chromosomes to plot, either those the user specified, or
    # the full genome.
    #
    if (length($chr_list) > 0) {
	my @c = split(/,/, $chr_list);
	foreach my $chr (@c) {
	    if (defined($orgs{$org}->{$chr})) {
		push(@chrs, $chr);
	    }
	}
	if (scalar(@chrs) == 0) {
	    print STDERR "Unable to locate the requested chromosomes, '$chr_list', in $org\n";
	    help();
	}
    } else {
	@chrs = keys %{$orgs{$org}};
    }
}

sub usage {
	print << "EOQ";
plot_concentric_fst.pl -o path -f fst_path [-f fst_path...] [-c chr1,chr2] [-O org] [-h]
  o:  write output to this file.
  f:  Stacks Fst file.
  c:  comma-seperated list of chromosomes to plot.
  O:  organism to plot, either 'gac' or 'mmu'.
  t:  image type, either 'gd' or 'pdf'.
  h:  display this help message.
  d:  turn on debug output.


EOQ

    exit(0);
}
