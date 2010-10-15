#!/usr/bin/perl
#
# Execute the Stacks Pipeline exporter script, wait for it to finish, send
# a notification email to the person who submitted the export.
#
# Written by Julian Catchen <jcatchen@uoregon.edu>
#

use strict;
use File::Temp qw/tempfile/;
use NET::SMTP;

#
# Configuration:
#  exe_path:    Path to the export executable.
#  output_path: Path to web-accessible directory to output the export data
#  url:         URL to reach the directory specified by output_path
#  local_host:  Name of localhost to present to SMTP server
#  smtp_host:   Name of SMTP server through which to send mail
#  from:        email address to use in the 'From' field of the message
#
my $exe_path      = "/opt/local/bin/perl /research/acos/export_sql.pl";
my $output_path   = "/Library/WebServer/Documents/teleost/radtags/export/"; 
my $url           = "http://teleost.cs.uoregon.edu/radtags/export/";
my $local_host    = "teleost.cs.uoregon.edu";
my $smtp_host     = "teleost.cs.uoregon.edu";
my $from          = "jcatchen\@uoregon.edu";


my $debug         = 0;
my $db            = "";
my $batch_id      = 0;
my $out_file      = "";     # Generated file name for output data
my $out_type      = "tsv";  # Output format
my $email         = "";     # Email address to send message to
my $filter_str    = "";     # Comma separated list of filters

parse_command_line();

my (@filters, @cmd_opts, $filter);

#
# Generate file name 
#
my $template = "stacks_export_XXXXXXXX";
my $suffix   = $out_type eq "xls" ? ".xls" : ".tsv";
my (undef, $out_file) = tempfile($template, OPEN => 0, DIR => $output_path, SUFFIX => $suffix);

#
# Prepare the command line parameters
#
push(@cmd_opts, "-D $db", "-b $batch_id", "-f $out_file", "-o $out_type");

if (length($filter_str) > 0) {
    @filters = split(/,/, $filter_str);

    foreach $filter (@filters) {
        push(@cmd_opts, "-F $filter");
    }
}

my $cmd = join(" ", @cmd_opts);
$cmd = $exe_path . " " . $cmd;
print STDERR "CMD: $cmd\n" if ($debug);

#
# Execute the exporter program
#
my @results = `$cmd`;
#my @results = `echo Success`;

#
# Check the results, we expext a one line result: either 'Success' or 'Failure'
#
chomp $results[0];

if ($results[0] eq "Success") {
    send_email('success', $out_file);
} else {
    send_email('failure', $out_file);
}


sub send_email {
    my ($result) = @_;

    my $smtp = Net::SMTP->new($smtp_host,
                              'Hello' => $local_host,
                              'Timeout' => 60);
    $smtp->mail($from);
    $smtp->recipient($email);

    my $msg .= 
        "From: $from\r\n" .
        "To: $email\r\n" .
        "Subject: Stacks pipeline export complete\r\n" .
        "\r\n";

    if ($result eq "success") {
        #
        # Trim the path off the output file
        #
        my ($f) = ($out_file =~ /.*\/(stacks_export_\w{8}\.\w{3})$/);

        $msg .= "Your data has been exported and can be downloaded from: " . $url . $f . "\r\n";
    } else {
        
        $msg .= "There has been an error exporting your data, please contact the system administrator.\r\n";
    }

    $smtp->data($msg);

    $smtp->quit();
}

sub parse_command_line {
    while (@ARGV) {
	$_ = shift @ARGV;
	if    ($_ =~ /^-d$/) { $debug++; }
	elsif ($_ =~ /^-D$/) { $db         = shift @ARGV; }
	elsif ($_ =~ /^-b$/) { $batch_id   = shift @ARGV; }
	elsif ($_ =~ /^-e$/) { $email      = shift @ARGV; }
	elsif ($_ =~ /^-t$/) { $out_type   = shift @ARGV; }
	elsif ($_ =~ /^-F$/) { $filter_str = shift @ARGV; }
	elsif ($_ =~ /^-h$/) { usage(); }
	else {
	    print STDERR "Unknown command line options received: $_\n";
	    usage();
	}
    }

    if (length($db) == 0) {
        print STDERR "You must specify a database to index.\n";
        usage();
    }

    if (length($email) == 0) {
        print STDERR "You must specify an email for notification.\n";
        usage();
    }

    if ($out_type ne "tsv" && $out_type ne "xls") {
        print STDERR "Output type can only be 'tsv' or 'xls'.\n";
        usage();
    }

    if ($batch_id !~ /^\d{1,4}$/) {
        print STDERR "Batch ID must be a numeric value.\n";
        usage();
    }

}

sub usage {
	print << "EOQ";
stacks_export_notify.pl -e email -D db -b batch_id [-t type] [-F filters] [-d] [-h]
  e: email to use for notification.
  D: radtag database to examine.
  b: batch_id of data set to export.
  t: output type, either 'tsv' or 'xls'.
  F: comma separated list of filters to apply to the data.
  h: display this help message.
  d: turn on debug output.


EOQ

    exit(0);
}
