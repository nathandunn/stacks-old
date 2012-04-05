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
require_once("header.php");
require_once("CatalogClass.php");

$database  = isset($_GET['db'])     ? $_GET['db']     : "";
$batch_id  = isset($_GET['id'])     ? $_GET['id']     : 0;
$email     = isset($_GET['email'])  ? $_GET['email']  : "";
$data_type = isset($_GET['dtype'])  ? $_GET['dtype']  : "haplo";
$map_type  = isset($_GET['mtype'])  ? $_GET['mtype']  : "gen";
$depth_lim = isset($_GET['dlim'])   ? $_GET['dlim']   : "1";
$man_cor   = isset($_GET['mcor'])   ? $_GET['mcor']   : "0";
$ex_type   = isset($_GET['otype'])  ? $_GET['otype']  : "tsv";

// Connect to the database
$db = db_connect($database);

// Save these variables for automatic URL formation later on.
$display = array();
$display['id']          = $batch_id;
$display['db']          = $database;
$display['filter_type'] = array();

$filters = array();

process_filter($display, $filters);

//
// Create a new catalog object to hold filtered loci
//
$catalog = new Catalog($db, $batch_id, $display);

$catalog->determine_count();
$loci_cnt = $catalog->num_loci();

$mc = ($man_cor > 0) ? "-c" : "";

if ($data_type == "haplo") {
  $dt = "-a haplo";
  $dl = "-L $depth_lim";
} else if ($data_type == "geno") {
  $dt = "-a geno -m $map_type $mc";
  $dl = "";
}

$cmd = $export_cmd . " -D $database -b $batch_id $dt $dl -e $email -t $ex_type -F " . implode(",", $filters);

header("Content-type: text/xml");
$xml_output = 
    "<?xml version=\"1.0\"?>\n" .
    "<export>\n" . 
    "<loci>" . number_format($loci_cnt) . "</loci>\n" .
    "<email>$email</email>\n" .
    "<msg>$cmd</msg>\n" .
    "</export>\n";

echo $xml_output;

//
// Execute the email notification program, which will run the export and email
// the submitter upon its conclusion.
//
// export_cmd is defined in the constants file.
//
exec($cmd . " >/dev/null &");

function process_filter(&$display_params, &$filters) {

    if (!isset($_GET['filter_type'])) 
	return;

    foreach ($_GET['filter_type'] as $filter) {
	array_push($display_params['filter_type'], $filter);

	if ($filter == "alle") {
	    $display_params['filter_alle_l'] = $_GET['filter_alle_l'];
            array_push($filters, "alle_l=" . $_GET['filter_alle_l']);
	    $display_params['filter_alle_u'] = $_GET['filter_alle_u'];
            array_push($filters, "alle_u=" . $_GET['filter_alle_u']);

	} else if ($filter == "snps") {
	    $display_params['filter_snps_l'] = $_GET['filter_snps_l'];
            array_push($filters, "snps_l=" . $_GET['filter_snps_l']);
	    $display_params['filter_snps_u'] = $_GET['filter_snps_u'];
            array_push($filters, "snps_u=" . $_GET['filter_snps_u']);

	} else if ($filter == "pare") {
	    $display_params['filter_pare_l'] = $_GET['filter_pare_l'];
            array_push($filters, "pare_l=" . $_GET['filter_pare_l']);
	    $display_params['filter_pare_u'] = $_GET['filter_pare_u'];
            array_push($filters, "pare_u=" . $_GET['filter_pare_u']);

	} else if ($filter == "prog") {
	    $display_params['filter_prog'] = $_GET['filter_prog'];
            array_push($filters, "prog=" . $_GET['filter_prog']);

	} else if ($filter == "vprog") {
	    $display_params['filter_vprog'] = $_GET['filter_vprog'];
            array_push($filters, "vprog=" . $_GET['filter_vprog']);

	} else if ($filter == "cata") {
	    $display_params['filter_cata'] = $_GET['filter_cata'];
            array_push($filters, "cata=" . $_GET['filter_cata']);

	} else if ($filter == "mark") {
	    $display_params['filter_mark'] = $_GET['filter_mark'];
            array_push($filters, "mark=" . $_GET['filter_mark']);

	} else if ($filter == "gcnt") {
	    $display_params['filter_gcnt'] = $_GET['filter_gcnt'];
	    array_push($filters, "gcnt=" . $_GET['filter_gcnt']);

	} else if ($filter == "ref") {
	    $display_params['filter_ref'] = $_GET['filter_ref'];
	    array_push($filters, "ref=" . $_GET['filter_ref']);

	} else if ($filter == "loc") {
	    $display_params['filter_chr'] = $_GET['filter_chr'];
	    array_push($filters, "chr=" . $_GET['filter_chr']);
	    $display_params['filter_sbp'] = $_GET['filter_sbp'];
	    array_push($filters, "sbp=" . $_GET['filter_sbp']);
	    $display_params['filter_ebp'] = $_GET['filter_ebp'];
	    array_push($filters, "ebp=" . $_GET['filter_ebp']);
	}
    }
}

?>
