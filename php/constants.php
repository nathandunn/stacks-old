<?php
//
// Name to print in page header
//
$site_title  = "Stacks Analysis Pipeline";

//
// URL path to root of PHP/HTML files
//
$root_path   = "/stacks";
$img_path    = "/stacks/images";

//
// File system location of PHP/HTML files
//
$system_path = "/Library/WebServer/Documents/teleost/stacks";

//
// Credentials to access Stacks MySQL databases
//
$db_user     = "acos_user";
$db_pass     = "orthologs";
$db_host     = "localhost";

//
// Location of perl script to run the exporting jobs.
//
$export_cmd  = "/opt/local/bin/perl /research/acos/stacks_export_notify.pl";

//
// Colors for printing alleles/haplotypes
//
$colors = array("#008000",
		"#c00000",
		"#ffc600",
		"#29356c",
		"#860000",
		"#dc6200",
		"#4b398e",
		"#008f56",
		"#bf1e25",
		"#4cb8ff");
$color_size = count($colors);

?>
