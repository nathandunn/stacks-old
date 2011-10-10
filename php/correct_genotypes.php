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

$database  = isset($_POST['db'])       ? $_POST['db']       : "";
$tag_id    = isset($_POST['tag_id'])   ? $_POST['tag_id']   : 0;
$batch_id  = isset($_POST['batch_id']) ? $_POST['batch_id'] : 0;
$op        = isset($_POST['op'])       ? $_POST['op']       : "display";

// Connect to the database
if (!isset($db)) $db = db_connect($database);

// Save these variables for automatic URL formation later on.
$display = array();
$display['db']       = $database;
$display['tag_id']   = $tag_id;
$display['batch_id'] = $batch_id;
$diplsay['op']       = $op;

//
// Prepare some SQL queries
//
$query = 
    "SELECT catalog_genotypes.id as id, catalog_genotypes.sample_id, catalog_genotypes.genotype, " . 
    "genotype_corrections.genotype as correction, genotype_corrections.id as cid, file " . 
    "FROM catalog_genotypes " . 
    "LEFT JOIN genotype_corrections ON " . 
    "(genotype_corrections.catalog_id=catalog_genotypes.catalog_id AND " .
    "genotype_corrections.sample_id=catalog_genotypes.sample_id AND " .
    "genotype_corrections.batch_id=catalog_genotypes.batch_id) " .
    "JOIN samples ON (catalog_genotypes.sample_id=samples.id) " . 
    "WHERE catalog_genotypes.batch_id=? and catalog_genotypes.catalog_id=? " . 
    "ORDER BY catalog_genotypes.sample_id";
$db['geno_sth'] = $db['dbh']->prepare($query);
check_db_error($db['geno_sth'], __FILE__, __LINE__);

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
    "SELECT genotype_corrections.id FROM genotype_corrections " . 
    "JOIN catalog_genotypes ON " . 
    "(genotype_corrections.catalog_id=catalog_genotypes.catalog_id AND " .
    "genotype_corrections.sample_id=catalog_genotypes.sample_id AND " .
    "genotype_corrections.batch_id=catalog_genotypes.batch_id) " .
    "WHERE genotype_corrections.batch_id=? AND genotype_corrections.catalog_id=?";
$db['res_sth'] = $db['dbh']->prepare($query);
check_db_error($db['res_sth'], __FILE__, __LINE__);

if ($op == "reset") {
    reset_marker($display);

} else if ($op == "correct") {
    correct_marker($display);
}

include_once("catalog_genotypes.php");

function reset_marker($display) {
    global $db;

    $result = $db['res_sth']->execute(array($display['batch_id'], $display['tag_id']));
    check_db_error($result, __FILE__, __LINE__);

    while ($row = $result->fetchRow()) {
        $r = $db['del_sth']->execute($row['id']);
        check_db_error($r, __FILE__, __LINE__);
    }
}

function correct_marker($display) {
    global $db;

    $gtypes      = array();
    $form_gtypes = array();
    //
    // Fetch the existing genotypes from the database
    //
    $result = $db['geno_sth']->execute(array($display['batch_id'], $display['tag_id']));
    check_db_error($result, __FILE__, __LINE__);

    while ($row = $result->fetchRow()) {
        $gtypes[$row['sample_id']] = array('id'           => $row['id'],
                                           'file'         => $row['file'], 
                                           'genotype'     => strtolower($row['genotype']),
                                           'corrected'    => $row['correction'],
                                           'corrected_id' => $row['cid']);
    }

    //
    // Fetch the corrected genotypes from the submitted form
    //
    foreach ($_POST as $key => $value) {
        if (substr($key, 0, 5) != "gtype") continue;

        // ID should look like: 'gtype_batchid_catalogid_sampleid'
        $parts = explode("_", $key);
        $form_gtypes[$parts[3]] = strtolower($value);
	//print "Assigning $value to $parts[3]<br />\n";
    }

    foreach ($form_gtypes as $sample_id => $sample) {
      //print "LOOKING at sample ID: $sample_id: $sample, original value: " .  $gtypes[$sample_id]['genotype'] . "<br />\n";

        //
        // Is this genotype being reset to the original value? If so, delete the corrected record.
        //
        if ($sample == $gtypes[$sample_id]['genotype'] && 
            strlen($gtypes[$sample_id]['corrected_id']) > 0) {
            $result = $db['del_sth']->execute($gtypes[$sample_id]['corrected_id']);
            check_db_error($result, __FILE__, __LINE__);
        //
        // Is the corrected value for this genotype being changed? If so, update the corrected record.
        //
        } else if ($sample != $gtypes[$sample_id]['genotype'] && 
                   strlen($gtypes[$sample_id]['corrected_id']) > 0) {
            $result = $db['upd_sth']->execute(array(strtoupper($sample), $gtypes[$sample_id]['corrected_id']));
            check_db_error($result, __FILE__, __LINE__);
        //
        // Otherwise, add a new correction.
        //
        } else if ($sample != $gtypes[$sample_id]['genotype']) {
            $result = $db['ins_sth']->execute(array($display['batch_id'], $display['tag_id'], $sample_id, strtoupper($sample)));
            check_db_error($result, __FILE__, __LINE__);
        }
    }
}

?>
