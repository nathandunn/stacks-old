#!/usr/bin/env python2

import optparse
import sys
import math
import re
from operator import itemgetter
import cairo

#
# Global configuration variables.
#
outtype        = "pdf"
outpath        = ""
chrpath        = ""
files          = []
file_types     = []
draw_filenames = True
draw_plotkey   = True

#
# The img dictionary will hold image-specific dimensions.
#
img = {}
img['height']        = 800
img['width']         = 800
img['min_rad_pct']   = 0.55
img['max_rad_pct']   = 0.8
img['font_size']     = 16;
img['fst_scale_max'] = 1.0
img['div_scale_max'] = 1.0

colors = []
colors.append((255.0/255.0, 255.0/255.0, 255.0/255.0)) # ffffff White
colors.append((  0.0/255.0,  80.0/255.0, 242.0/255.0)) # 0050f2 Blue
colors.append((255.0/255.0,   0.0/255.0,   0.0/255.0)) # ff0000 Red

cols = {'fst'     : 18,
        'phist'   : 7,
        'fstp'    : 10,
        'hapdiv'  : 11,
        'hapavg'  : 11,
        'hapdiff' :  8}

#
# Global variables to hold plotting data.
#
chrs         = {}
chrs_default = []
chrs_sorted  = []
stats        = {}
means        = {}

def parse_file_type(option, opt, value, parser):
    global files
    global file_types

    files.append(value)
    if opt == "-f" or opt == "--fst":
        file_types.append('fst')
    elif opt == "-F" or opt == "--fstprine":
        file_types.append('fstp')
    elif opt == "-P" or opt == "--phist":
        file_types.append('phist')
    elif opt == "-H" or opt == "--hapdiv":
        file_types.append('hapdiv')
    elif opt == "-D" or opt == "--hapdiff":
        file_types.append('hapdiff')
    elif opt == "-A" or opt == "--hapavg":
        file_types.append('hapavg')

def parse_command_line(img):
    global outpath
    global outtype
    global chrpath
    global files
    global file_types
    global chrs_sorted
    global draw_filenames
    global draw_plotkey

    p = optparse.OptionParser()

    #
    # Add options.
    #
    p.add_option("-o", action="store", dest="outfile",
                 help="write output to this file.")
    p.add_option("-t", "--outtype", action="store", dest="outtype",
                 help="output type, either 'pdf' or 'svg'.")
    p.add_option("-c", "--chrs", action="store", dest="chrpath",
                 help="chromosome definition file: two column, chromosome and length, in plot order.")
    p.add_option("-C", action="store", dest="chr_list",
                 help="only plot this comma-seperated list of chromosomes.")
    p.add_option("-m", action="store", type="float", dest="min_rad",
                 help="minimum radius, between 0 and 1.0 (default 0.55).")
    p.add_option("-M", action="store", type="float", dest="max_rad",
                 help="maximum radius, between 0 and 1.0 (default 0.8).")
    p.add_option("-s", action="store", type="float", dest="img_width",
                 help="set height/width of image in points.")
    p.add_option("-n", "--no_filenames", action="store_false", dest="draw_filenames",
                 help="Suppress drawing of input filenames.")
    p.add_option("-k", "--no_key", action="store_false", dest="draw_key",
                 help="Suppress drawing of plot key.")
    p.add_option("-f", "--fst", action="callback", type="string", callback=parse_file_type,
                 help="Stacks Fst statistic from 'fst' file.")
    p.add_option("-F", "--fstprime", action="callback", type="string", callback=parse_file_type,
                 help="Stacks Fst' statistic from 'phistats' file.")
    p.add_option("-P", "--phist", action="callback", type="string", callback=parse_file_type,
                 help="Stacks Phi_st statistic from 'phistats' file.")
    p.add_option("-H", "--hapdiv", action="callback", type="string", callback=parse_file_type, 
                 help="Stacks haplotype statistics file, 'hapstats', will plot raw haplotype diversity.")
    p.add_option("-D", "--hapdiff", action="callback", type="string", callback=parse_file_type,
                 help="Stacks haplotype statistics difference file, 'hapstats_diff', will plot difference between two haplotype diversities.")
    p.add_option("-A", "--hapavg", action="callback", type="string", callback=parse_file_type,
                 help="Stacks haplotype statistics file, 'hapstats', will plot difference from haplotype diversity average.")
    p.add_option("--hapdiv_scale", action="store", type="float", dest="hapdiv_scale",
                 help="Maximum haplotype diversity value for drawing the color scale.")
    p.add_option("--fst_scale", action="store", type="float", dest="fst_scale",
                 help="Maximum Fst/Phist/Fst' value for drawing the color scale.")
    
    #
    # Parse the command line
    #
    (opts, args) = p.parse_args()

    if opts.outfile != None:
        outpath = opts.outfile
    if opts.outtype != None:
        outtype = opts.outtype
    if opts.chrpath != None:
        chrpath = opts.chrpath
    if opts.fst_scale != None:
        img['fst_scale_max'] = opts.fst_scale
    if opts.hapdiv_scale != None:
        img['div_scale_max'] = opts.hapdiv_scale
    if opts.min_rad != None:
        img['min_rad_pct'] = opts.min_rad
    if opts.max_rad != None:
        img['max_rad_pct'] = opts.max_rad
    if opts.img_width != None:
        img['width']  = opts.img_width
        img['height'] = opts.img_width
    if opts.draw_filenames != None:
        draw_filenames = opts.draw_filenames
    if opts.draw_key != None:
        draw_plotkey = opts.draw_key
        
    if len(outpath) == 0:
        print >> sys.stderr, "You must specify a path to write files to."
        p.print_help()
        sys.exit()

    if img['min_rad_pct'] < 0 or img['min_rad_pct'] > 1.0:
        print >> sys.stderr, "Minimum radius must be between 0 and 1.0."
        p.print_help()
        sys.exit()

    if img['max_rad_pct'] < 0 or img['max_rad_pct'] > 1.0:
        print >> sys.stderr, "Maximum radius must be between 0 and 1.0."
        p.print_help()
        sys.exit()

    if outtype != "svg" and outtype != "pdf":
        print >> sys.stderr, "Output type must be either 'pdf' or 'svg'."
        p.print_help()
        sys.exit()
        
    #
    # Generate a list of chromosomes to plot, either those the user specified, or
    # the full genome.
    #
    if opts.chr_list != None:
        c = opts.chr_list.split(",")
        for chr in c:
            if chr in chrs:
                chrs_sorted.append(chr)
        if len(chrs_sorted) == 0:
            print >> sys.stderr, "Unable to locate the requested chromosomes, '", opts.chr_list, "'"
            p.print_help()
            sys.exit()
    else:
        chrs_sorted = chrs_default

    #
    # Input files have to have unique names to be placed in the hash later.
    #
    for i in range(len(files)):
        files[i] = str(i + 1) + ". " + files[i]


def determine_chromosome_lengths(img, chrs, chrs_sorted):
    #
    # Determine the total length of the chromosomes
    #
    num_chrs      = len(chrs_sorted)
    total_chr_len = 0.0

    for chr in chrs_sorted:
        total_chr_len += chrs[chr]

    #
    # Determine the total size of the arc (in degrees) after 
    # removing space for margins around each chromosome.
    #
    total_arc_len = float(360 - (num_chrs * img['chr_margin']))

    #
    # Calculate arc length with each arc's size relative to the 
    # length of the other chromosomes. Ensure no arcs are less than a
    # single degree.
    #
    # In the cairo library, 0 degrees corresponds to 3 o'clock. We want the
    # chromosomes to draw from 12 o'clock, so we start at 270 degrees.
    #
    degrees = 270
    img['chrs'] = {}
    
    for chr in chrs_sorted:
        img['chrs'][chr] = {}

	arc_len   = total_arc_len * (chrs[chr] / total_chr_len)
	arc_len   = 1 if arc_len < 1 else arc_len

	start_deg = degrees   + img['chr_margin']
	end_deg   = start_deg + arc_len

	img['chrs'][chr]['start_deg'] = start_deg % 360
	img['chrs'][chr]['end_deg']   = end_deg % 360

	degrees = end_deg

    #
    # Record the average arc length for later use in draw_chromosome_arc().
    #
    img['avg_arc'] = total_arc_len / num_chrs

def get_x_y_coordinates(degree, radius):
    if degree <= 90:
        theta      = float(degree)
        opp_side   = radius * math.sin(math.radians(theta))
        adj_side   = radius * math.cos(math.radians(theta))
        x          = img['center_x'] + adj_side
        y          = img['center_y'] + opp_side
    elif degree <= 180:
        theta      = float(degree - 90.0)
        opp_side   = radius * math.sin(math.radians(theta))
        adj_side   = radius * math.cos(math.radians(theta))
        x          = img['center_x'] - opp_side
        y          = img['center_y'] + adj_side
    elif degree <= 270:
        theta      = float(degree - 180.0)
        opp_side   = radius * math.sin(math.radians(theta))
        adj_side   = radius * math.cos(math.radians(theta))
        x          = img['center_x'] - adj_side
        y          = img['center_y'] - opp_side
    else:
        theta      = float(degree - 270.0)
        opp_side   = radius * math.sin(math.radians(theta))
        adj_side   = radius * math.cos(math.radians(theta))
        x          = img['center_x'] + opp_side
        y          = img['center_y'] - adj_side

    return (x, y)

def draw_filenames(cr, img, files, file_types):
    textents    = cr.text_extents(files[0])
    text_height = textents[3]
    y = 5 + text_height
    i = len(files)
    j = i - 1;

    cr.select_font_face("Sans", cairo.FONT_SLANT_NORMAL, cairo.FONT_WEIGHT_NORMAL)
    cr.set_font_size(img['font_size'] - 4)
    cr.set_source_rgb(0, 0, 0) # Black
    
    while (j >= 0):
	f     = files[j];
	ftype = file_types[j];

        p = re.compile('.*\/(\w.+)$')
        m = p.match(f)
        file = m.group(1)
        
	if ftype == "fst":
	    file = "Fst; " + file
	elif ftype == "fstp":
	    file = "Fst'; " + file
	elif ftype == "phist":
	    file = "Phi_st; " + file
	elif ftype == "hapdiv":
	    file = "HapDiv; " + file
	elif ftype == "hapavg":
	    file = "HapAvg; " + file
	elif ftype == "hapdif":
	    file = "HapDif; " + file

        text = str(i) + ". " + file
        
        #
        # Get the size of the text we want to write, returns a tuple:
        #   (x, y, width, height, dx, dy)
        #
        textents = cr.text_extents(text)
        text_height = textents[3]

        cr.move_to(5, y);
        cr.show_text(text)

	y += text_height + 2
	i -= 1
        j -= 1

def draw_scale_bars(cr, img, chrs_sorted):

    min_rad = img['min_radius'] - math.floor(img['width'] / 2.0 * 0.03)
    max_rad = img['max_radius'] + math.floor(img['width'] / 2.0 * 0.08)

    cr.set_dash([4, 1])

    genome_len = 0.0
    tic_size   = 15000000

    for chr in chrs_sorted:
        genome_len += chrs[chr]
    if genome_len < 500000000:
        tic_size = 5000000

    for chr in chrs_sorted:
	c       = img['chrs'][chr]
	chr_len = chrs[chr]

	#
	# Determine chromosome length rounded to nearest 10Mb.
	#
	tics = int(chr_len / tic_size) % tic_size

	for i in range(1, tics + 1):
	    start_bp = tic_size * i

	    #
	    # Scale the starting point for proper placement on the chromosome arc
	    #
            if c['end_deg'] < c['start_deg']:
                start_degree = float((360 - c['start_deg']) + c['end_deg'])
            else:
	        start_degree = float(c['end_deg'] - c['start_deg'])
	    start_degree = (c['start_deg'] + ((float(start_bp) / float(chr_len)) * start_degree)) % 360

	    #
	    # Convert the point to an x,y coordinate
	    #
            (x1, y1) = get_x_y_coordinates(start_degree, min_rad)
            (x2, y2) = get_x_y_coordinates(start_degree, max_rad)

	    if start_bp % 10000000 == 0:
                cr.set_source_rgb(90.0/255.0, 90.0/255.0, 90.0/255.0) # Dark Grey
                cr.set_line_width(2)
	    else:
		cr.set_source_rgb(127.0/255.0, 127.0/255.0, 127.0/255.0) # Grey
                cr.set_line_width(1)

	    cr.move_to(x1, y1)
	    cr.line_to(x2, y2)
	    cr.stroke()
 
            #
            # Draw a label every 10Mb on the interior of the circle.
            #
            cr.set_font_size(img['font_size'] - 6)
            cr.set_source_rgb(0, 0, 0)

            if start_bp % 10000000 == 0:
                #
                # Get the size of the text we want to write, returns a tuple:
                #   (x, y, width, height, dx, dy)
                #
                text = str(start_bp / 1000000) + "Mb"
                textents = cr.text_extents(text)
                text_width  = textents[2]
                text_height = textents[3]

        	# Adjust it, depending on the quadrant of the circle we are in
        	w = 0
                h = 0
                if start_degree <= 90:
                    # Bottom right quadrant
                    w -= text_width
                    h -= text_height / 2.0
                elif start_degree <= 180:
                    # Bottom left quadrant
                    h += text_height / 2.0
                elif start_degree <= 270:
                    # Top left quadrant
                    h += text_height
                else:
                    # Top right quadrant
                    w -= text_width
                    h += text_height

                cr.move_to(x1+w, y1+h);
                cr.show_text(text)
    cr.set_dash([])

def draw_chromosome_arc(cr, img, chr, radius, labels):

    start_deg = img['chrs'][chr]['start_deg']
    end_deg   = img['chrs'][chr]['end_deg']

    # print >> sys.stderr, "Drawing arc: (", img['center_x'], ",", img['center_y'], "), radius: ", radius, "; start degree:", start_deg, "; end degree:", end_deg
    
    x = 0.0
    y = 0.0
    #
    # Draw the inside arc.
    #
    (x, y) = get_x_y_coordinates(start_deg, radius)
    cr.move_to(x, y)
    cr.set_line_width(img['chr_border'])
    cr.arc(img['center_x'], img['center_y'], radius, math.radians(start_deg), math.radians(end_deg))

    #
    # Determine the starting X,Y position of the outside arc.
    #
    radius_out = radius + float(img['chr_height'])
    (x, y) = get_x_y_coordinates(end_deg, radius_out)
    cr.line_to(x, y)
    cr.arc_negative(img['center_x'], img['center_y'], radius_out, math.radians(end_deg), math.radians(start_deg))
    cr.close_path()
    cr.set_source_rgb(0, 0, 0) # Black
    cr.stroke_preserve()
    cr.set_source_rgb(127.0/255.0, 127.0/255.0, 127.0/255.0) # Grey, RGB: #7f7f7f
    cr.fill()


def draw_chromosome_labels(cr, img, chr, radius):
    start_deg = img['chrs'][chr]['start_deg']
    end_deg   = img['chrs'][chr]['end_deg']

    #
    # Draw the chromosome name
    #
    regex    = re.compile('group')
    text     = regex.sub('', chr)

    if end_deg < start_deg:
        text_loc = (start_deg + float((360 - start_deg) + end_deg) / 2.0) % 360
    else:
        text_loc = float(end_deg + start_deg) / 2.0

    (x, y) = get_x_y_coordinates(text_loc, img['max_radius'] + 40)

    cr.select_font_face("Sans", cairo.FONT_SLANT_NORMAL, cairo.FONT_WEIGHT_NORMAL)
    cr.set_font_size(img['font_size'])
    cr.set_source_rgb(0, 0, 0) # Black

    #
    # Get the size of the text we want to write, returns a tuple:
    #   (x, y, width, height, dx, dy)
    #
    textents = cr.text_extents(text)
    text_width  = textents[2]
    text_height = textents[3]

    # Adjust it, depending on the quadrant of the circle we are in
    w = 0.0
    h = 0.0
    if text_loc <= 90:
        # Bottom right quadrant
        h += text_height / 2.0
    elif text_loc <= 180:
        # Bottom left quadrant
        w -= text_width 
        h += text_height
    elif text_loc <= 270:
        # Top left quadrant
        w -= text_width 

    cr.move_to(x+w, y+h);
    cr.show_text(text)

def hsv2rgb(H, s, v):
    #
    # Our value of $h is a fraction between 0 and 1, scale it to a number between 0 and 360.
    #
    h = H * 360.0

    #
    # Conversion routine from Perl Monks: http://www.perlmonks.org?node_id=139485 
    #
    if s == 0:
        return v, v, v

    h /= 60
    i = math.floor(h)
    f = h - i
    p = v * (1 - s)
    q = v * (1 - (s * f))
    t = v * (1 - (s * (1 - f)))

    if i == 0:
        return v, t, p
    elif i == 1:
        return q, v, p
    elif i == 2:
        return p, v, t
    elif i == 3:
        return p, q, v
    elif i == 4:
        return t, p, v
    else:
        return v, p, q

def scale_fst_color(val):
    min_range = 1.0
    max_range = 10.0
    min_color = 0.0
    max_color = 0.25

    if val == 0:
        return max_color

    val = (max_range - min_range) * val + 1.0

    color = (math.log(val)/ math.log(10)) * (max_color - min_color)

    return max_color - color

def scale_log_color(val):
    if val == 0:
        return 0 

    color = math.log(1.0 + (math.fabs(val) * 9.0)) / math.log(10)

    return color

def scale_sigmoid_color(val):
    val = math.fabs(val)

    #
    # Plot the values using a sigmoid function: S(t) = 1 / 1 + e^-t
    #
    # Scale val to be between -6 and 6, as that is the range on the X axis
    # for a response variable between 0 and 1.0.
    #
    #val = val * 12.0 - 6.0
    #val = 1.0 / (1.0 + math.exp(-1 * val))

    #
    # Plot the values using a Gompertz curve: y(t) = ae^(-be(^-ct))
    #   https://en.wikipedia.org/wiki/Gompertz_function
    #
    # The X-axis scales from 0 to 10, Y-axis from 0 to 1; we have shifted
    # the curve to the right by subtracting 2 from val.
    #
    a = 1
    b = 1
    c = 1
    val = val * 10
    val = a * math.exp(-1 * b * (math.exp(-1 * c * (val-2))))

    return val

def scale_sigmoid_color_mean(mean, val):
    #
    # Plot the values using a Gompertz curve: y(t) = ae^(-be(^-ct))
    #   https://en.wikipedia.org/wiki/Gompertz_function
    #
    # The X-axis scales from 0 to 10, Y-axis from 0 to 1; we have shifted
    # the curve to the right by subtracting 2 from val; we will shit it more
    # to line the sigmoid curve up with the mean value.
    #
    val  = math.fabs(val)
    val  = val * 10
    mean = mean * 10 / 2
    a    = 1
    b    = 1
    c    = 0.5
    val  = a * math.exp(-1 * b * (math.exp(-1 * c * (val - 2 - mean))))
    return val

def scale_hapdiv_color(val):
    min_range = 1.0
    max_range = 10.0

    min_color = 0.5
    max_color = 0.85

    #min_color = 0.666
    #max_color = 0.9
    
    if val == 0:
        return max_color

    #
    # Set maximum value for haplotype diversity to 0.6. Any haplotype diversity values
    # higher than that will be set to 0.6. Then scale the value to cover the full
    # min_color/max_color range.
    #
    if val > 0.7:
        val = 0.7
    val   = val * (10.0/7.0)
    color = val * (max_color - min_color)
        
    return min_color + color

def two_color_gradient((r1, g1, b1), (r2, g2, b2), alpha):
    r = ((1.0 - alpha) * r1) + (alpha * r2)
    g = ((1.0 - alpha) * g1) + (alpha * g2)
    b = ((1.0 - alpha) * b1) + (alpha * b2)
    return (r, g, b)

def three_color_gradient((r1, g1, b1), (r2, g2, b2), (r3, g3, b3), mean, alpha):
    alpha = scale_sigmoid_color_mean(mean, alpha)
    if (alpha < mean):
        scaled_alpha = float(alpha) / mean
        r = ((1.0 - scaled_alpha) * r1) + (scaled_alpha * r2)
        g = ((1.0 - scaled_alpha) * g1) + (scaled_alpha * g2)
        b = ((1.0 - scaled_alpha) * b1) + (scaled_alpha * b2)
    else:
        scaled_alpha = float(alpha - mean) / (1.0 - mean)
        r = (scaled_alpha * r3) + ((1.0 - scaled_alpha) * r2)
        g = (scaled_alpha * g3) + ((1.0 - scaled_alpha) * g2)
        b = (scaled_alpha * b3) + ((1.0 - scaled_alpha) * b2)
    return (r, g, b)

def draw_key(cr, img, key_type, index, mean):
    x1 = img['width'] - img['key_start_x'][index]
    x2 = x1 + img['key_width']
    y1 = img['key_start_y']
    y2 = y1 + img['key_height']
    y  = y2

    cr.set_line_width(img['chr_border'])
    cr.rectangle(x1, y1, img['key_width'], img['key_height'])
    cr.set_source_rgb(0, 0, 0) # Black
    cr.stroke_preserve()
    cr.set_source_rgb(127.0/255.0, 127.0/255.0, 127.0/255.0) # Grey, RGB: #7f7f7f
    cr.fill()

    cr.set_line_width(1)

    if key_type == "hapdiff" or key_type == "hapavg":
        scale_max  = 1.0
        color_step = 1.0 / (float(img['key_height']) / 2.0)
        text       = "Diff" if key_type == "hapdiff" else "Avg"
        i          = -1
        while (i <= 1.0):
            if i >= 0:
                if key_type == "hapdiff":
                    (r, g, b) = hsv2rgb(0.61, scale_sigmoid_color(i), 0.95)    
                else:
                    (r, g, b) = hsv2rgb(0.61, scale_sigmoid_color(i), 0.95)
            else:
                if key_type == "hapdiff":
                    (r, g, b) = hsv2rgb(0, scale_sigmoid_color(i), 1.0)
                else:
                    (r, g, b) = hsv2rgb(0.083, scale_sigmoid_color(i), 1.0)

            cr.set_source_rgb(r, g, b)
            cr.move_to(x1, y)
            cr.line_to(x2, y)
            cr.stroke()
            y -= 1
            i += color_step

    else:
        scale_max  = img['div_scale_max'] if key_type == "hapdiv" else img['fst_scale_max']
        color_step = 1.0 / float(img['key_height'])
        text       = "Div" if key_type == "hapdiv" else "Fst"
        i          = 0
        while (i <= 1.0):
            if key_type == "hapdiv":
                # (r, g, b) = hsv2rgb(scale_hapdiv_color(i), 1.0, 0.85)
                (r, g, b) = three_color_gradient(colors[0], colors[1], colors[2], mean, i)
            else:
                (r, g, b) = hsv2rgb(scale_fst_color(i), 1.0, 0.85)

            cr.set_source_rgb(r, g, b)
            cr.move_to(x1, y)
            cr.line_to(x2, y)
            cr.stroke()
            y -= 1
            i += color_step

    #
    # Draw the key title and scale bars
    #
    cr.set_font_size(img['font_size'] - 4)
    cr.set_source_rgb(0, 0, 0) # Black

    # Get the size of the text we want to write
    textents    = cr.text_extents(text)
    text_width  = textents[2]
    text_height = textents[3]

    cr.move_to((float(x1 + x2) / 2.0) - float(text_width / 2.0), y1 - 5)
    cr.show_text(text)

    cr.set_font_size(img['font_size'] - 6)

    if key_type == "hapdiff" or key_type == "hapavg":
        text = "-1.0"
    else:
        text = "0"
    textents    = cr.text_extents(text)
    text_width  = textents[2]
    text_height = textents[3]
    cr.move_to(x1 - text_width - 3, y2 + float(text_height / 2.0));
    cr.show_text(text)

    if key_type == "hapdiff" or key_type == "hapavg":
        text = "0"
    else:
        text = str(round(scale_max / 2.0, 2))
    textents    = cr.text_extents(text)
    text_width  = textents[2]
    text_height = textents[3]
    cr.move_to(x1 - text_width - 3, float((y1 + y2) / 2.0) + float(text_height / 2.0))
    cr.show_text(text)

    text        = str(round(scale_max, 2))
    textents    = cr.text_extents(text)
    text_width  = textents[2]
    text_height = textents[3]
    cr.move_to(x1 - text_width - 3, y1  + float(text_height / 2.0))
    cr.show_text(text)

def bucket_stat_values(img, chr, stats, radius, bcnt):

    buckets  = []
    c        = img['chrs'][chr]
    chr_len  = chrs[chr]
    stat_cnt = len(stats[chr]['bp'])

    #
    # How many basepairs can we fit into each bucket?
    #
    blimits = [];
    bps = math.ceil(chr_len / bcnt)

    # print >> sys.stderr, "Chromosome", chr, chrs[chr], "bps;", bcnt, "buckets;", bps, "basepairs per bucket."

    #
    # Create an array of buckets, divided into 1/8ths of a degree.
    #
    lim = 0;
    for i in range(0, int(bcnt)):
	lim += bps
	buckets.append([])
	blimits.append(lim)

    #
    # Sort the stat values.
    #
    sorted_stats = [];
    for i in range(0, stat_cnt):
	bp   = stats[chr]['bp'][i]
	stat = stats[chr]['stat'][i]
	sorted_stats.append((bp, stat))
    sorted_stats = sorted(sorted_stats, key=itemgetter(0))

    #
    # Iterate over the sorted statistical values and place them into the proper bucket.
    #
    lim = blimits.pop(0)
    i   = 0
    for (bp, stat) in sorted_stats:
	while (lim != None and bp > lim):
	    i += 1
	    lim = blimits.pop(0)
	buckets[i].append(stat)

    return buckets

def draw_stat_values(cr, img, stats, mean, chr, radius, file_type):
    #
    # Determine the number of degrees occupied by this chromosome.
    #
    buckets_per_degree = 8

    #
    # Handle the case where a chromosome crosses over the 0 degree boundary such that
    # the end degree is smaller than the start degree.
    #
    if img['chrs'][chr]['end_deg'] < img['chrs'][chr]['start_deg']:
        deg = (360 - img['chrs'][chr]['start_deg']) + img['chrs'][chr]['end_deg']
    else:
        deg = img['chrs'][chr]['end_deg'] - img['chrs'][chr]['start_deg']

    bcnt    = int(deg * buckets_per_degree)
    buckets = bucket_stat_values(img, chr, stats, radius, bcnt)

    cnt       = 0
    start_deg = float(img['chrs'][chr]['start_deg'])

    cr.set_line_width(img['stat_border'])
    
    for bucket in buckets:
        #print >> sys.stderr, "start deg:", start_deg, "bucket len:", len(bucket)
	end_deg = start_deg + 1.0/float(buckets_per_degree)

        if len(bucket) == 0:
            start_deg = end_deg
            cnt += 1
            continue

	mean_stat = 0
	if file_type == "hapavg":
	    #
	    # Convert the haplotype diversity numbers to the difference from the average value.
	    #
	    diffavg = []
	    for hdiv in bucket:
		diffavg.append(hdiv - mean)
            mean_stat = float(sum(diffavg)) / float(len(diffavg))
        else:
	    mean_stat = float(sum(bucket)) / float(len(bucket))

        if file_type == "hapdiff":
	    # Blue/Red, sigmoid scaled
            if mean_stat >= 0:
		(r, g, b) = hsv2rgb(0.61, scale_sigmoid_color(mean_stat), 0.95)
	    else:
		(r, g, b) = hsv2rgb(0, scale_sigmoid_color(mean_stat), 1.0)

        elif file_type == "hapavg":
	    # Blue/Brown
	    if mean_stat >= 0:
		(r, g, b) = hsv2rgb(0.61, scale_sigmoid_color(mean_stat), 0.95)
            else:
		(r, g, b) = hsv2rgb(0.083, scale_sigmoid_color(mean_stat), 1.0)

        elif file_type == "hapdiv":
            scale_factor = 1.0 / img['div_scale_max']
            if mean_stat > img['div_scale_max']:
                mean_stat = img['div_scale_max']
	    # Purple
	    # (r, g, b) = hsv2rgb(0.791, scale_hapdiv_color(mean_stat), 1.0)
	    # (r, g, b) = hsv2rgb(scale_hapdiv_color(mean_stat), 1.0, 0.85)
            (r, g, b) = three_color_gradient(colors[0], colors[1], colors[2], mean, mean_stat * scale_factor)
        else:
	    # Red/Green
            scale_factor = 1.0 / img['fst_scale_max']
            if mean_stat > img['fst_scale_max']:
                mean_stat = img['fst_scale_max']
	    (r, g, b) = hsv2rgb(scale_fst_color(mean_stat * scale_factor), 1.0, 0.85)

        x = 0.0
        y = 0.0
        #
        # Draw the inside arc.
        #
        (x, y) = get_x_y_coordinates(start_deg, radius)
        cr.move_to(x, y)
        cr.arc(img['center_x'], img['center_y'], radius, math.radians(start_deg), math.radians(end_deg))

        #
        # Determine the starting X,Y position of the outside arc.
        #
        radius_out = radius + float(img['chr_height'])
        (x, y) = get_x_y_coordinates(end_deg, radius_out)
        cr.line_to(x, y)
        cr.arc_negative(img['center_x'], img['center_y'], radius_out, math.radians(end_deg), math.radians(start_deg))
        cr.close_path()
        cr.set_source_rgb(r, g, b)
        cr.stroke_preserve()
        cr.fill()

	start_deg = end_deg
	cnt += 1

def generate_img(img, chrs, chrs_sorted, files, file_types):

    if outtype == "svg":
        ps = cairo.SVGSurface(outpath, img['height'], img['width'])
    else:
        ps = cairo.PDFSurface(outpath, img['height'], img['width'])
    cr = cairo.Context(ps)

    #
    # Determine the pixel boundaries for the image.
    #
    img['chr_margin']    = 2 # degrees
    img['center_x']      = int(img['width']  / 2.0)
    img['center_y']      = int(img['height'] / 2.0)
    img['circle_width']  = int(img['width']  * 0.8)
    img['circle_height'] = int(img['height'] * 0.8)
    img['min_radius']    = int(img['width']  / 2.0 * img['min_rad_pct'])
    img['max_radius']    = int(img['width']  / 2.0 * img['max_rad_pct'])
    img['chr_height']    = 24
    img['chr_border']    = 2 # px
    img['stat_border']   = 0.5 # px
    img['key_width']     = 25
    img['key_height']    = 125
    img['key_start_y']   = 20
    img['key_start_x']   = [40, 85, 130, 175]
    
    if draw_filenames:
        draw_filenames(cr, img, files, file_types)

    if draw_plotkey:
        keys = {}
        j    = 0
        #
        # Determine the keys that need to be drawn.
        #
        for i in range(len(file_types)):
            t = file_types[i]
            if t == "phist" or t == "fstp":
                t = "fst"
            if t not in keys:
                keys[t] = 1
                draw_key(cr, img, t, j, means[files[i]])
                j += 1

    #
    # Determine the length of each chromosome in degrees.
    #
    determine_chromosome_lengths(img, chrs, chrs_sorted)

    #
    # Draw scale bars across the circles.
    #
    draw_scale_bars(cr, img, chrs_sorted)

    #
    # Determine the number of concentric circles to draw, and radial spacing for each.
    #
    num_circles = len(files)
    radius_dist = 1
    if num_circles > 1:
        radius_dist = (img['max_radius'] - img['min_radius']) / (num_circles - 1)

    #
    # Draw the individual chromosomes as arcs along a circle
    #
    radius = float(img['min_radius'])
    labels = False

    i = 1
    for j in range(num_circles):
	print >> sys.stderr, "Processing sample", files[j], "; type:", file_types[j], "..."
	fpath = files[j]
	ftype = file_types[j]

	labels = True if i == num_circles else False

	for chr in chrs_sorted:
            # print >> sys.stderr, "Drawing chromosome", chr
	    draw_chromosome_arc(cr, img, chr, radius, labels)
            if labels:
                draw_chromosome_labels(cr, img, chr, radius)
            draw_stat_values(cr, img, stats[fpath], means[fpath], chr, radius, ftype)

	radius += radius_dist
	i      += 1

    cr.show_page()
    ps.finish()

def parse_chromosome_definitions(chrpath, chrs, chrs_default):

    if len(chrpath) > 0:
        fh = open(chrpath, "r")

        for line in fh:
            line = line.strip("\n")

            if len(line) == 0 or line[0] == "#":
                continue

            parts = line.split("\t")
            if len(parts) != 2:
                print >> sys.stderr, ("Cannot parse chromosome definition file, '" +
                                      chrpath +
                                      "'; should be two, tab-separated columns, found " + len(parts))
                exit(1)
            chrs[parts[0]] = int(parts[1])
            chrs_default.append(parts[0])

        fh.close()
    else:
        #
        # Define chromosome sizes in stickleback.
        #
        c = {'groupI'    : 28185914, 'groupII'   : 23295652, 'groupIII'   : 16798506, 'groupIV'  : 32632948,
             'groupIX'   : 20249479, 'groupV'    : 12251397, 'groupVI'    : 17083675, 'groupVII' : 27937443,
             'groupVIII' : 19368704, 'groupX'    : 15657440, 'groupXI'    : 16706052, 'groupXII' : 18401067,
             'groupXIII' : 20083130, 'groupXIV'  : 15246461, 'groupXIX'   : 20240660, 'groupXV'  : 16198764,
             'groupXVI'  : 18115788, 'groupXVII' : 14603141, 'groupXVIII' : 16282716, 'groupXX'  : 19732071,
             'groupXXI'  : 11717487}
        d = ['groupI', 'groupII', 'groupIII', 'groupIV', 'groupV', 'groupVI', 'groupVII',
             'groupVIII', 'groupIX', 'groupX', 'groupXI', 'groupXII', 'groupXIII', 'groupXIV', 
             'groupXV', 'groupXVI', 'groupXVII', 'groupXVIII', 'groupXIX', 'groupXX', 'groupXXI']
        for chr in d:
            chrs_default.append(chr)
            chrs[chr] = c[chr]

    print >> sys.stderr, "Parsed " + str(len(chrs_default)) + " chromosomes."


def populate_stats(files, file_types, stats, means):

    for i in range(len(files)):
	file  = files[i]
	ftype = file_types[i]
	path  = file[3:]

	stats[file] = {}
        mean = []

        fh = open(path, "r")

	for line in fh:
	    line = line.strip("\n")

	    if len(line) == 0 or line[0] == "#":
                continue

	    parts = line.split("\t")

	    if ftype == "hapdiv" or ftype == "hapavg":
                chr = parts[2]
                bp  = int(parts[3])
            else:
                chr = parts[4]
                bp  = int(parts[5])

            if chr not in chrs:
                continue

            if chr not in stats[file]:
                stats[file][chr] = {}
                stats[file][chr]['bp']   = []
                stats[file][chr]['stat'] = []

            stats[file][chr]['bp'].append(bp)
            
            if ftype == "hapdiff":
		if parts[0] == "oc":
		    stats[file][chr]['stat'].append(-1 * float(parts[cols[ftype]]))
		else:
		    stats[file][chr]['stat'].append(float(parts[cols[ftype]]))
	    else:
		stats[file][chr]['stat'].append(float(parts[cols[ftype]]))

	    mean.append(float(parts[cols[ftype]]))

	fh.close()

	means[file] = float(sum(mean)) / float(len(mean))

            
parse_command_line(img)

print >> sys.stderr, "Parsing chromosome definitions from '" + chrpath + "'..."
parse_chromosome_definitions(chrpath, chrs, chrs_default)
print >> sys.stderr, "Populating statistical measures for each input file..."
populate_stats(files, file_types, stats, means)
print >> sys.stderr,  "done."

generate_img(img, chrs, chrs_sorted, files, file_types)
