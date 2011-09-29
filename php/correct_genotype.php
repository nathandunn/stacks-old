<?php
//
// Copyright 2011, Julian Catchen <jcatchen@uoregon.edu>
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

$database  = isset($_GET['db'])        ? $_GET['db']        : "";
$tag_id    = isset($_GET['tag_id'])    ? $_GET['tag_id']    : 0;
$batch_id  = isset($_GET['batch_id'])  ? $_GET['batch_id']  : 0;
$sample_id = isset($_GET['sample_id']) ? $_GET['sample_id'] : 0;
$gtype     = isset($_GET['gtype'])     ? $_GET['gtype']     : 0;

// Connect to the database
$db = db_connect($database);

//
// Prepare some SQL queries
//
$query = 
    "UPDATE genotype_corrections SET genotype=? WHERE id=?";
$db['upd_sth'] = $db['dbh']->prepare($query);
check_db_error($db['upd_sth'], __FILE__, __LINE__);

$query = 
    "INSERT INTO genotype_corrections SET batch_id=?, catalog_id=?, sample_id=?, genotype=?";
$db['ins_sth'] = $db['dbh']->prepare($query);
check_db_error($db['ins_sth'], __FILE__, __LINE__);

$query = 
    "DELETE FROM genotype_corrections WHERE id=?";
$db['del_sth'] = $db['dbh']->prepare($query);
check_db_error($db['del_sth'], __FILE__, __LINE__);

$query = 
    "SELECT catalog_genotypes.id as id, catalog_genotypes.sample_id, catalog_genotypes.genotype, " . 
    "genotype_corrections.genotype as correction, genotype_corrections.id as cid " . 
    "FROM catalog_genotypes " . 
    "LEFT JOIN genotype_corrections ON " . 
    "(genotype_corrections.catalog_id=catalog_genotypes.catalog_id AND " .
    "genotype_corrections.sample_id=catalog_genotypes.sample_id AND " .
    "genotype_corrections.batch_id=catalog_genotypes.batch_id) " .
    "WHERE catalog_genotypes.batch_id=? AND " . 
    "catalog_genotypes.catalog_id=? AND " .
    "catalog_genotypes.sample_id=?";
$db['geno_sth'] = $db['dbh']->prepare($query);
check_db_error($db['geno_sth'], __FILE__, __LINE__);

//
// Fetch the existing genotypes from the database
//
$result = $db['geno_sth']->execute(array($batch_id, $tag_id, $sample_id));
check_db_error($result, __FILE__, __LINE__);

if ($row = $result->fetchRow()) {
    $sample = array('id'           => $row['id'], 
		    'genotype'     => strtolower($row['genotype']),
		    'corrected'    => $row['correction'],
		    'corrected_id' => $row['cid']);
} else {
    return;
}

$corrected = "false";

if ($gtype == "clr") {
    $result = $db['del_sth']->execute($sample['corrected_id']);
    check_db_error($result, __FILE__, __LINE__);

    $corrected = "false";
    $gtype = $sample['genotype'];
//
// Is this genotype being reset to the original value? If so, delete the corrected record.
//
} else if ($gtype == $sample['genotype'] && 
    strlen($sample['corrected_id']) > 0) {
    $result = $db['del_sth']->execute($sample['corrected_id']);
    check_db_error($result, __FILE__, __LINE__);

    $corrected = "false";
//
// Is the corrected value for this genotype being changed? If so, update the corrected record.
//
} else if ($gtype != $sample['genotype'] && 
	   strlen($sample['corrected_id']) > 0) {
    $result = $db['upd_sth']->execute(array(strtoupper($gtype), $sample['corrected_id']));
    check_db_error($result, __FILE__, __LINE__);

    $corrected = "true";
//
// Otherwise, add a new correction.
//
} else if ($gtype != $sample['genotype']) {
    $result = $db['ins_sth']->execute(array($batch_id, $tag_id, $sample_id, strtoupper($gtype)));
    check_db_error($result, __FILE__, __LINE__);

    $corrected = "true";
}

$id = "gtype_" . $batch_id . "_" . $tag_id . "_" . $sample_id;
if ($corrected == "true")
    $gtype = strtoupper($gtype);

header("Content-type: text/xml");
$xml_output = 
    "<?xml version=\"1.0\"?>\n" .
    "<correction>\n" . 
    "<corrected>$corrected</corrected>\n" .
    "<gtype>$gtype</gtype>\n" .
    "<div_id>$id</div_id>\n" .
    "</correction>\n";

echo $xml_output;

?>
