#!/usr/bin/perl
#
# Load a set of output files from the Stacks pipeline into a Stacks MySQL database.
#
# In order to differentiate parents from progeny, the script expects parents
# to have 'male' and/or 'female' as part of the filename.
#
# By Julian Catchen <jcatchen@uoregon.edu>
#

use strict;

my $debug      = 0;
my $db         = "";
my $in_path    = ".";
my $sample_id  = 0;
my $desc       = ""; #"Lepisosteus oculatus RAD-Tag Samples";
my $date       = ""; #"2009-05-31";
my $batch_id   = 0;
my $batch      = 0;
my $catalog    = 0;
my $white_list = "";

parse_command_line();

my (@files, @catalog);

build_file_list(\@files, \@catalog);

print STDERR scalar(@files), " files to process: ", join(", ", @files), "\n";
print STDERR "Catalog files to process: ", join(", ", @catalog), "\n";

if ($batch) {
    `mysql $db -e "INSERT INTO batches SET id=$batch_id, description='$desc', date='$date'"`;
}

my ($file, $f, $i, $type);

$i = 1;
foreach $file (@files) {
    $f = $in_path . "/$file" . ".tags.tsv";

    #
    # Pull out the sample ID and insert it into the database
    #
    $type = ($file =~ /f?e?male/) ? 'parent' : 'progeny';
    $sample_id = extract_sample_id($f);

    print STDERR "Sample: $file; Type: $type; Sample ID: $sample_id; i: $i\n";
    `mysql $db -e "INSERT INTO samples SET id=$sample_id, sample_id=$i, batch_id=$batch_id, type='$type', file='$file'"`;

    #
    # Import the unique tag files.
    #
    import_sql_file($f, "unique_tags");

    $f = $in_path . "/$file" . ".snps.tsv";
    import_sql_file($f, "snps");

    $f = $in_path . "/$file" . ".alleles.tsv";
    import_sql_file($f, "alleles");

    $f = $in_path . "/$file" . ".matches.tsv";
    import_sql_file($f, "matches");

    $i++;
}

#
# Import the catalog
#
if ($catalog) {
    foreach $file (@catalog) {
        $f = $in_path . "/$file" . ".catalog.tags.tsv";
        import_sql_file($f, "catalog_tags");

        $f = $in_path . "/$file" . ".catalog.snps.tsv";
        import_sql_file($f, "catalog_snps");

        $f = $in_path . "/$file" . ".catalog.alleles.tsv";
        import_sql_file($f, "catalog_alleles");
    }
}

sub extract_sample_id {
    my ($file) = @_;

    my ($line, @results, @parts);

    @results = `head -n 1 $file`;
    chomp $results[0];
    @parts = split(/\t/, $results[0]);

    #
    # Sample ID is expected to be the first column in the *.tags.tsv file.
    #
    return $parts[1];
}

sub import_sql_file {
    my ($file, $table) = @_;

    my (@results);

    @results = `mysql $db -e "LOAD DATA LOCAL INFILE '$file' INTO TABLE $table"`;
    print STDERR "mysql $db -e \"LOAD DATA LOCAL INFILE '$file' INTO TABLE $table\"\n", @results, "\n";
}

sub build_file_list {
    my ($files, $catalog_files) = @_;

    my (@wl, @ls, $line, $prefix);

    # Load a white list of files to process if it is supplied.
    if (length($white_list) > 0) {
        load_white_list(\@wl);
    }

    @ls = `ls -1 $in_path/*.tags.tsv 2> /dev/null`;

    if (scalar(@ls) == 0) {
	print STDERR "Unable to locate any input files to process within '$in_path'\n";
	usage();
    }

    foreach $line (@ls) {
	chomp $line;

	($prefix) = ($line =~ /$in_path\/(.+)\.tags\.tsv/);

        next if ($prefix =~ /catalog/);

        if (scalar(@wl) > 0) {
            next if (!grep(/^$prefix$/, @wl));
        }

	push(@{$files}, $prefix);
    }

    if ($catalog > 0) {
        @ls = `ls -1 $in_path/*.catalog.tags.tsv 2> /dev/null`;

        if (scalar(@ls) == 0) {
            print STDERR "Unable to locate any catalog input files to process within '$in_path'\n";
            usage();
        }

        foreach $line (@ls) {
            chomp $line;

            ($prefix) = ($line =~ /$in_path\/(.+)\.catalog.tags\.tsv/);

            if (scalar(@wl) > 0) {
                next if (!grep(/^$prefix$/, @wl));
            }

            push(@{$catalog_files}, $prefix);
        }
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
	if    ($_ =~ /^-p$/) { $in_path    = shift @ARGV; }
	elsif ($_ =~ /^-D$/) { $db         = shift @ARGV; }
	elsif ($_ =~ /^-c$/) { $catalog++; }
        elsif ($_ =~ /^-B$/) { $batch++; }
	elsif ($_ =~ /^-b$/) { $batch_id   = shift @ARGV; }
	elsif ($_ =~ /^-e$/) { $desc       = shift @ARGV; }
	elsif ($_ =~ /^-a$/) { $date       = shift @ARGV; }
        elsif ($_ =~ /^-W$/) { $white_list = shift @ARGV; }
	elsif ($_ =~ /^-d$/) { $debug++; }
	elsif ($_ =~ /^-h$/) { usage(); }
	else {
	    print STDERR "Unknown command line option: '$_'\n";
	    usage();
	}
    }

    $in_path  = substr($in_path, 0, -1)  if (substr($in_path, -1)  eq "/");

    if (length($db) == 0) {
	print STDERR "You must specify a database.\n";
	usage();
    }
}

sub usage {
    print STDERR <<EOQ; 
load_radtags.pl -D db -p path -b batch_id [-B -a date -e desc] [-c] [-W path] [-d] [-h]
    D: Database to load data into.
    p: Path to input files.
    b: Batch ID.
    c: Load the catalog into the database.
    B: Load information into batch table.
    e: batch dEscription.
    a: batch run dAte, yyyy-mm-dd.
    W: white list of files to load.
    h: display this help message.
    d: turn on debug output.

EOQ

exit(0);
}
