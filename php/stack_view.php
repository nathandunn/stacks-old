<?php
//
// Copyright 2010-2015, Julian Catchen <jcatchen@illinois.edu>
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

$tag_str    = isset($_GET['tags'])     ? $_GET['tags']     : "";
$batch_id   = isset($_GET['batch_id']) ? $_GET['batch_id'] : 0;
$cat_id     = isset($_GET['tag_id'])   ? $_GET['tag_id']   : "";
$sample_str = isset($_GET['samples'])  ? $_GET['samples']  : "";
$database   = isset($_GET['db'])       ? $_GET['db']       : "";

$tags    = explode(",", $tag_str);
$samples = explode(",", $sample_str);

// Connect to the database
$db = db_connect($database);

//
// Prepare some SQL queries
//
$query = 
    "SELECT tag_id, sub_id, relationship, seq_id, seq, deleveraged, blacklisted, removed, file, pop_id " . 
    "FROM unique_tags " . 
    "JOIN samples ON (unique_tags.sample_id=samples.id) " . 
    "JOIN batches ON (samples.batch_id=batches.id) " .
    "WHERE batch_id=? AND unique_tags.sample_id=? AND tag_id=?";
$db['seq_sth'] = $db['dbh']->prepare($query);
check_db_error($db['seq_sth'], __FILE__, __LINE__);

$query = 
    "SELECT col, rank_1, rank_2 FROM snps " . 
    "JOIN samples ON (snps.sample_id=samples.id) " . 
    "JOIN batches ON (samples.batch_id=batches.id) " .
    "WHERE batch_id=? AND snps.sample_id=? AND tag_id=? AND snps.type='E' ORDER BY col";
$db['snp_sth'] = $db['dbh']->prepare($query);
check_db_error($db['snp_sth'], __FILE__, __LINE__);

$query = 
    "SELECT col FROM catalog_snps " . 
    "WHERE batch_id=? AND tag_id=?";
$db['cat_snp_sth'] = $db['dbh']->prepare($query);
check_db_error($db['cat_snp_sth'], __FILE__, __LINE__);

$json_str = 
    "{" .
    "\"id\": \"$cat_id\"," .
    "\"stacks\": [";

for ($i = 0; $i < count($tags); $i++) {
    $sample_id = $samples[$i];
    $tag_id    = $tags[$i];

    //
    // Fetch and store the SNPs.
    //
    $snps   = array();
    $result = $db['snp_sth']->execute(array($batch_id, $sample_id, $tag_id));
    check_db_error($result, __FILE__, __LINE__);
    while ($row = $result->fetchRow()) {
        $snps[$row['col']] = array('col' => $row['col'], 'rank_1' => $row['rank_1'], 'rank_2' => $row['rank_2']);
    }

    //
    // Fetch and store the catalog SNPs.
    //
    $cat_snps = array();
    $result   = $db['cat_snp_sth']->execute(array($batch_id, $cat_id));
    check_db_error($result, __FILE__, __LINE__);
    while ($row = $result->fetchRow()) {
        if (!isset($cat_snps[$row['col']]))
            $cat_snps[$row['col']] = 0;
        $cat_snps[$row['col']]++;
    }

    //
    // Fetch the sequences.
    //
    $deleveraged     = 0;
    $lumberjackstack = 0;
    $blacklisted     = 0;
    $consensus       = "";
    $model           = "";
    $file            = "";
    $stacks          = array();
    $secondary       = array();
    $result = $db['seq_sth']->execute(array($batch_id, $sample_id, $tag_id));
    check_db_error($result, __FILE__, __LINE__);

    while ($row = $result->fetchRow()) {
    
        if ($row['relationship'] == "consensus") {
            $deleveraged     = $row['deleveraged'];
            $lumberjackstack = $row['removed'];
            $blacklisted     = $row['blacklisted'];
            $consensus       = $row['seq'];
            $file            = $row['file'];

        } else if ($row['relationship'] == "model") {
            $model = $row['seq'];

        } else if ($row['relationship'] == "secondary") {
            array_push($secondary, array('s' => $row['seq'], 'id' => $row['seq_id']));

        } else {
            if (!isset($stacks[$row['sub_id']]))
                $stacks[$row['sub_id']] = array();
            array_push($stacks[$row['sub_id']], array('s' => $row['seq'], 'id' => $row['seq_id']));
        }
    }

    $con_len = strlen($consensus);

    $c = print_snps($tags[$i], $consensus, $consensus, $snps, false);
    $c = addslashes($c);

    $scale = print_scale($con_len);
    $scale = addslashes($scale);

    $json_str .= 
        "{" .
        "\"sample_id\": \"$sample_id\"," .
        "\"sample_name\": \"$file\"," .
        "\"tag_id\": \"$tag_id\"," .
        "\"consensus\": \"$c\"," .
        "\"model\": \"$model\"," .
        "\"scale\": \"$scale\"," .
        "\"primary\": [";

    foreach ($stacks as $sub_id => $stack) {
        $s = print_snps_errs($consensus, $stack[0]['s'], $snps, $cat_snps);
        $s = addslashes($s);

        $json_str .=
            "{" .
            "\"seq\": \"$s\"," .
            "\"ids\": [";

        foreach ($stack as $seq) {
            $json_str .=
                "\"$seq[id]\",";
        }

        if (count($stack) > 0)
            $json_str = substr($json_str, 0, -1);
        $json_str .= 
            "]" .
            "},";
    }

    if (count($stacks) > 0)
        $json_str = substr($json_str, 0, -1);
    $json_str .=
        "]," .
        "\"secondary\": [";

    foreach ($secondary as $seq) {
        $s = print_snps_errs($consensus, $seq['s'], $snps, $cat_snps);
        $s = addslashes($s);

        $json_str .=
            "{" .
            "\"id\": \"$seq[id]\"," .
            "\"seq\": \"$s\"" .
            "},";
    }

    if (count($secondary) > 0)
        $json_str = substr($json_str, 0, -1);
    $json_str .= 
        "]},";
}

$json_str  = substr($json_str, 0, -1);
$json_str .= "]}";

echo $json_str;

?>
