<?php
//
// Copyright 2010, Julian Catchen <jcatchen@uoregon.edu>
//
// This file is part of Stacks.
//
// Stacks is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// Stacks is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with Stacks.  If not, see <http://www.gnu.org/licenses/>.
//

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
// Name to print in page header
//
$site_title  = "Stacks Analysis Pipeline";

//
// URL path to root of PHP/HTML files
//
$root_path   = "/stacks";
$img_path    = "/stacks/images";

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
