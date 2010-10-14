<?php
require_once("header.php");
require_once("CatalogClass.php");
require_once("radtag_functions.php");

$database  = isset($_GET['db'])     ? $_GET['db']     : "";
$tag_id    = isset($_GET['tag_id']) ? $_GET['tag_id'] : 0;
$batch_id  = isset($_GET['id'])     ? $_GET['id']     : 0;
$email     = isset($_GET['email'])  ? $_GET['email']  : "";
$ex_type   = isset($_GET['type'])   ? $_GET['type']   : "";

// Connect to the database
$db = db_connect($database);

// Save these variables for automatic URL formation later on.
$display = array();
$display['id']          = $batch_id;
$display['db']          = $database;
$display['p']           = $page;
$display['pp']          = $per_page;
$display['filter_type'] = array();

process_filter($display);

//
// Create a new catalog object to hold filtered loci
//
$catalog = new Catalog($db, $batch_id, $display);

$catalog->determine_count();
$loci_cnt = $catalog->num_loci();

header("Content-type: text/xml");
$xml_output = 
    "<?xml version=\"1.0\"?>\n" .
    "<export>\n" . 
    "<loci>" . number_format($loci_cnt) . "</loci>\n" .
    "<email>$email</email>\n" .
    "</export>\n";

echo $xml_output;

function process_filter(&$display_params) {

    if (!isset($_GET['filter_type'])) 
	return;

    foreach ($_GET['filter_type'] as $filter) {
	array_push($display_params['filter_type'], $filter);

	if ($filter == "alle") {
	    $display_params['filter_alle'] = $_GET['filter_alle'];

	} else if ($filter == "snps") {
	    $display_params['filter_snps'] = $_GET['filter_snps'];

	} else if ($filter == "pare") {
	    $display_params['filter_pare'] = $_GET['filter_pare'];

	} else if ($filter == "prog") {
	    $display_params['filter_prog'] = $_GET['filter_prog'];

	} else if ($filter == "vprog") {
	    $display_params['filter_vprog'] = $_GET['filter_vprog'];

	} else if ($filter == "cata") {
	    $display_params['filter_cata'] = $_GET['filter_cata'];

	} else if ($filter == "mark") {
	    $display_params['filter_mark'] = $_GET['filter_mark'];

	}
    }
}

?>
