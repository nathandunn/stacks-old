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

#
# Load a set of output files from the Stacks pipeline into a Stacks MySQL database.
#
# In order to differentiate parents from progeny, the script expects parents
# to have 'male' and/or 'female' as part of the filename.
#
# By Julian Catchen <jcatchen@uoregon.edu>
#

use strict;

use constant stacks_version => "_VERSION_";

my $mysql_config  = "_PKGDATADIR_" . "sql/mysql.cnf";
my $dry_run       = 0;
my $db            = "";
my $in_path       = ".";
my $sample_id     = 0;
my $desc          = ""; #"Lepisosteus oculatus RAD-Tag Samples";
my $date          = ""; #"2009-05-31";
my $batch_id      = 0;
my $batch         = 0;
my $catalog       = 0;
my $white_list    = "";

parse_command_line();

my (@results, @files, @catalog, @parent_ids);

build_file_list(\@files, \@catalog);
extract_parental_ids(\@catalog, \@parent_ids);

print STDERR 
    scalar(@files), " files to process: ", join(", ", @files), "\n",
    "Catalog files to process: ", join(", ", @catalog), "\n",
    "Parent IDs identified: ", join(", ", @parent_ids), "\n";

if ($batch) {
    if (!$dry_run) {
	@results = `mysql --defaults-file=$mysql_config $db -e "INSERT INTO batches SET id=$batch_id, description='$desc', date='$date'"`;
    }
    print STDERR
	"mysql --defaults-file=$mysql_config $db ",
	"-e \"INSERT INTO batches SET id=$batch_id, description='$desc', date='$date'\"\n",
	@results, "\n";
}

my ($file, $f, $i, $type);

$i = 1;
foreach $file (@files) {
    $f = $in_path . "/$file" . ".tags.tsv";

    #
    # Pull out the sample ID and insert it into the database
    #
    $sample_id = extract_sample_id($f);
    $type = (grep(/^$sample_id$/, @parent_ids) > 0) ? 'parent' : 'progeny';

    if (!$dry_run) {
	@results = `mysql --defaults-file=$mysql_config $db -e "INSERT INTO samples SET id=$sample_id, sample_id=$i, batch_id=$batch_id, type='$type', file='$file'"`;
    }
    print STDERR 
	"mysql --defaults-file=$mysql_config $db ",
	"-e \"INSERT INTO samples SET id=$sample_id, sample_id=$i, batch_id=$batch_id, type='$type', file='$file'\"\n", 
	@results, "\n";

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

sub extract_parental_ids {
    my ($catalog, $parental_ids) = @_;

    my ($fh, $prefix, $path, $line, @parts, $tag_id, @tag_ids, $id, $tag, %ids);

    foreach $prefix (@catalog) {
        $path = $in_path . "/" . $prefix . ".catalog.tags.tsv";

	open($fh, "<$path") or die("Unable to open catalog file: '$path', $!\n");

	while ($line = <$fh>) {
	    chomp $line;
	    @parts = split(/\t/, $line);

	    @tag_ids = split(/,/, $parts[8]);

	    foreach $tag_id (@tag_ids) {
		($id, $tag) = split(/_/, $tag_id);
		$ids{$id}++;
	    }
	}

	close($fh);
    }

    @{$parental_ids} = keys %ids;
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

    if (!$dry_run) {
	@results = `mysql --defaults-file=$mysql_config $db -e "LOAD DATA LOCAL INFILE '$file' INTO TABLE $table"`;
    }
    print STDERR 
	"mysql --defaults-file=$mysql_config $db ",
	"-e \"LOAD DATA LOCAL INFILE '$file' INTO TABLE $table\"\n", 
	@results, "\n";
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
	elsif ($_ =~ /^-d$/) { $dry_run++; }
	elsif ($_ =~ /^-v$/) { version(); exit(); }
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

sub version {
    print STDERR "load_radtags.pl ", stacks_version, "\n";
}

sub usage {
    version();

    print STDERR <<EOQ; 
load_radtags.pl -D db -p path -b batch_id [-B -a date -e desc] [-c] [-d] [-W path] [-d] [-h]
    D: Database to load data into.
    p: Path to input files.
    b: Batch ID.
    c: Load the catalog into the database.
    B: Load information into batch table.
    e: batch dEscription.
    a: batch run dAte, yyyy-mm-dd.
    d: perform a dry run. Do not actually load any data, just print what would be executed.
    W: only load file found on this white list.
    h: display this help message.

EOQ

exit(0);
}
