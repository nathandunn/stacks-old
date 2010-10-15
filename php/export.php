<?php
require_once("header.php");
require_once("CatalogClass.php");

// $_GET['id']=2;
// $_GET['db']="trout_radtags";
// $_GET['p']=1;
// $_GET['filter_type']=array("snps");
// $_GET['filter_snps']=1;
// $_GET['pp']="all";

$batch_id  = isset($_GET['id']) ? $_GET['id'] : 0;
$database  = isset($_GET['db']) ? $_GET['db'] : "";
$page      = isset($_GET['p'])  ? $_GET['p']  : 1;
$per_page  = isset($_GET['pp']) ? $_GET['pp'] : 10;

// Connect to the database
$db = db_connect($database);

//
// Prepare some SQL queries
//
$query = 
    "SELECT id, file FROM samples " . 
    "WHERE batch_id=? ORDER BY id";
$db['samp_sth'] = $db['dbh']->prepare($query);
check_db_error($db['samp_sth'], __FILE__, __LINE__);

//
// Pull list of samples for this batch
//
$samples = array();
$result = $db['samp_sth']->execute($batch_id);
check_db_error($result, __FILE__, __LINE__);

while ($row = $result->fetchRow()) {
    $samples[$row['file']] = $row['id'];
}

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
$pagination_count = $catalog->num_loci();

$start_group = 0;
$end_group   = 0;

calculate_page_size($pagination_count, $start_group, $end_group, $num_pages, $per_page);

// Open the temporary file to write data to
$date = strftime("%m%d%y%H%M%S", time());

$tmp_file = $database . "_export_" . $date . ".xls";
$tmp_path = $system_path . "/export/" . $tmp_file;

// Prepare the command line to execute
$excel = "/opt/local/bin/perl /research/acos/export_excel.pl -f $tmp_path";

// Execute the perl program that will produce our excel file. The perl script
// will read in tab-seperated records from STDIN and print them to an excel
// file.
//$excel_handle = popen($excel, "w");
$excel_handle = fopen($tmp_path, "w");

//
// Populate our ParaGroups object, giving the object a hint
// as to which genes it needs to load.
//
$result = $catalog->populate($start_group - 1, $end_group - $start_group + 1);
$loci   = $catalog->loci();

// Print the heading out for the spreadsheet
$str = 
    "Catalog ID\t" . 
    "Annotation\t" . 
    "Chr\t" .
    "BP\t" .
    "Consensus Sequence\t" .
    "Num Parents\t" .
    "Num Progeny\t" . 
    "Num SNPs\t" .
    "SNPs\t" .
    "Num Alleles\t" .
    "Alleles\t";

foreach ($samples as $sample => $id) {
    $str .= $id . "\t";
}
$str  = substr($str, 0, -1);
$str .= "\n";

fwrite($excel_handle, $str);

foreach ($loci as $cat_id => $locus) {
    $str =
        $cat_id . "\t" .
        $locus->annotation . "\t" .
        $locus->chr . "\t" .
        $locus->bp . "\t" .
        $locus->seq . "\t" .
        $locus->num_parents . "\t" .
        $locus->num_progeny . "\t" .
        $locus->num_snps . "\t" .
        $locus->snps . "\t" . 
        $locus->num_alleles . "\t" .
        $locus->alleles . "\t";

    fwrite($excel_handle, $str);

    $gtypes = "";
    foreach ($samples as $id => $sample) {
        $types = $locus->genotype($id);

        if (!isset($types)) {
            $gtypes .= "\t";
            continue;
        }

        foreach ($types as $type)
            $gtypes .= $type['allele'] . "/";

        $gtypes = substr($gtypes, 0, -1);
        $gtypes .= "\t";
    }
    $gtypes = substr($gtypes, 0, -1);
    $gtypes .= "\n";

    fwrite($excel_handle, $gtypes);
}

$str = "\n";
foreach ($samples as $id => $sample) {
    $str .= $sample . "\t" . $id . "\n";
}
fwrite($excel_handle, $str);

//pclose($excel_handle);
fclose($excel_handle);

//
// Upload the file back to the client
//
header('Content-type: application/vnd.ms-excel');
header("Content-Disposition: attachment; filename=\"" . $tmp_file . "\";");
header("Content-Transfer-Encoding: binary");
header("Content-Length: " . filesize($tmp_path)); 
readfile($tmp_path);

// Delete the temp file
unlink($tmp_path);


function calculate_page_size($num_genes, &$start_locus, &$end_locus, &$num_pages, &$per_page) {
    global $display;

    $page     = $display['p'];
    $per_page = $display['pp'];

    if ($per_page == "all") 
        $per_page = $num_genes == 0 ? 1 : $num_genes;

    //
    // First figure out the total number of pages. If there are
    // additional genes left over, add an extra page.
    //
    $num_pages = floor($num_genes / $per_page);
    $num_pages += $num_genes % $per_page >= 1 ? 1 : 0;

    if ($page > $num_pages) {
        $page = $num_pages;
    }

    // Determine the start and end gene numbers
    $start_locus = 1 + (($page - 1) * $per_page);
    $end_locus   = 
        $start_locus + $per_page > $num_genes ? 
        $num_genes : 
        ($start_locus + $per_page - 1);

}

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
