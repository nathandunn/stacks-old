<?php
//
// Copyright 2015-2016, Julian Catchen <jcatchen@illinois.edu>
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

$database   = isset($_GET['db'])       ? $_GET['db']       : "";
$tag_id     = isset($_GET['tag_id'])   ? $_GET['tag_id']   : 0;
$batch_id   = isset($_GET['batch_id']) ? $_GET['batch_id'] : 0;
$batch_type = isset($_GET['type'])     ? $_GET['type']     : "map";

// Connect to the database
$db = db_connect($database);

// Save these variables for automatic URL formation later on.
$display = array();
$display['db']       = $database;
$display['tag_id']   = $tag_id;
$display['batch_id'] = $batch_id;

//
// Prepare some SQL queries
//
$query = 
    "SELECT samples.id, samples.sample_id, samples.type, file, tag_id, allele, depth, lnl, pop_id " . 
    "FROM matches " . 
    "JOIN samples ON (matches.sample_id=samples.id) " . 
    "WHERE matches.batch_id=? AND catalog_id=? ORDER BY samples.id";
if (!($db['mat_sth'] = $db['dbh']->prepare($query)))
    write_db_error($db['dbh'], __FILE__, __LINE__);

$query = 
    "SELECT col, rank_1, rank_2, rank_3, rank_4 FROM catalog_snps " . 
    "WHERE batch_id=? AND tag_id=? ORDER BY col";
if (!($db['snp_sth'] = $db['dbh']->prepare($query)))
    write_db_error($db['dbh'], __FILE__, __LINE__);

$query = 
    "SELECT allele FROM catalog_alleles " . 
    "WHERE batch_id=? AND tag_id=? ";
if (!($db['all_sth'] = $db['dbh']->prepare($query)))
      write_db_error($db['dbh'], __FILE__, __LINE__);

$query = 
    "SELECT geno_map FROM markers " . 
    "WHERE batch_id=? AND catalog_id=? ";
if (!($db['map_sth'] = $db['dbh']->prepare($query)))
      write_db_error($db['dbh'], __FILE__, __LINE__);

$query = 
    "SELECT pop_id, pop_name FROM populations " . 
    "WHERE batch_id=?";
if (!($db['pop_sth'] = $db['dbh']->prepare($query)))
      write_db_error($db['dbh'], __FILE__, __LINE__);

$query = 
    "SELECT count(batch_id) as cnt FROM sumstats " . 
    "WHERE batch_id=? AND tag_id=?";
if (!($db['stats_sth'] = $db['dbh']->prepare($query)))
      write_db_error($db['dbh'], __FILE__, __LINE__);

$query = 
    "SELECT count(batch_id) as cnt FROM fst " . 
    "WHERE batch_id=? AND tag_id=?";
if (!($db['fst_sth'] = $db['dbh']->prepare($query)))
    write_db_error($db['dbh'], __FILE__, __LINE__);

$query = 
    "SELECT count(batch_id) as cnt FROM hapstats " . 
    "WHERE batch_id=? AND tag_id=?";
if (!($db['hapstats_sth'] = $db['dbh']->prepare($query)))
    write_db_error($db['dbh'], __FILE__, __LINE__);

$query = 
    "SELECT count(batch_id) as cnt FROM phist " . 
    "WHERE batch_id=? AND tag_id=?";
if (!($db['phist_sth'] = $db['dbh']->prepare($query)))
    write_db_error($db['dbh'], __FILE__, __LINE__);


$query = 
    "SELECT marker, catalog_genotypes.sample_id, file, " . 
    "catalog_genotypes.genotype, genotype_corrections.genotype as corrected " . 
    "FROM catalog_genotypes " . 
    "LEFT JOIN genotype_corrections ON " . 
    "(genotype_corrections.catalog_id=catalog_genotypes.catalog_id AND " .
    "genotype_corrections.sample_id=catalog_genotypes.sample_id AND " .
    "genotype_corrections.batch_id=catalog_genotypes.batch_id) " .
    "JOIN samples ON (catalog_genotypes.sample_id=samples.id) " . 
    "JOIN catalog_index ON (catalog_genotypes.catalog_id=catalog_index.tag_id AND " . 
    "catalog_genotypes.batch_id=catalog_index.batch_id) " .
    "WHERE catalog_genotypes.batch_id=? and catalog_genotypes.catalog_id=? " . 
    "ORDER BY catalog_genotypes.sample_id";
if (!($db['geno_sth'] = $db['dbh']->prepare($query)))
    write_db_error($db['dbh'], __FILE__, __LINE__);

//
// Check for the existence of SNP summary statistics or Fst data.
//
$snp_sumstats = 0;
$snp_fst_vals = 0;
$hap_sumstats = 0;
$hap_fst_vals = 0;
if ($batch_type == "population") {
    if (!$db['stats_sth']->bind_param("ii", $batch_id, $tag_id))
	write_db_error($db['stats_sth'], __FILE__, __LINE__);
    if (!$db['stats_sth']->execute())
	write_db_error($db['stats_sth'], __FILE__, __LINE__);
    $res = $db['stats_sth']->get_result();

    if ($row = $res->fetch_assoc()) {
	if ($row['cnt'] > 0)
	    $snp_sumstats = 1;
    }

    if (!$db['fst_sth']->bind_param("ii", $batch_id, $tag_id))
	write_db_error($db['fst_sth'], __FILE__, __LINE__);
    if (!$db['fst_sth']->execute())
	write_db_error($db['fst_sth'], __FILE__, __LINE__);
    $res = $db['fst_sth']->get_result();

    if ($row = $res->fetch_assoc()) {
	if ($row['cnt'] > 0)
	    $snp_fst_vals = 1;
    }

    if (!$db['hapstats_sth']->bind_param("ii", $batch_id, $tag_id))
	write_db_error($db['hapstats_sth'], __FILE__, __LINE__);
    if (!$db['hapstats_sth']->execute())
	write_db_error($db['hapstats_sth'], __FILE__, __LINE__);
    $res = $db['hapstats_sth']->get_result();

    if ($row = $res->fetch_assoc()) {
	if ($row['cnt'] > 0)
	    $hap_sumstats = 1;
    }

    if (!$db['phist_sth']->bind_param("ii", $batch_id, $tag_id))
	write_db_error($db['phist_sth'], __FILE__, __LINE__);
    if (!$db['phist_sth']->execute())
	write_db_error($db['phist_sth'], __FILE__, __LINE__);
    $res = $db['phist_sth']->get_result();

    if ($row = $res->fetch_assoc()) {
	if ($row['cnt'] > 0)
	    $hap_fst_vals = 1;
    }
}

//
// Fetch population names if available.
//
$pop_names = array();
if ($batch_type == "population") {
    if (!$db['pop_sth']->bind_param("i", $batch_id))
	write_db_error($db['pop_sth'], __FILE__, __LINE__);
    if (!$db['pop_sth']->execute())
	write_db_error($db['pop_sth'], __FILE__, __LINE__);
    $res = $db['pop_sth']->get_result();

    while ($row = $res->fetch_assoc())
        $pop_names[$row['pop_id']] = $row['pop_name'];
}

if (!$db['snp_sth']->bind_param("ii", $batch_id, $tag_id))
    write_db_error($db['snp_sth'], __FILE__, __LINE__);
if (!$db['snp_sth']->execute())
    write_db_error($db['snp_sth'], __FILE__, __LINE__);
$res = $db['snp_sth']->get_result();

$json_str = 
  "{" .
  "\"path\": \"$root_path\"," .
  "\"batch_id\": \"$batch_id\"," .
  "\"db\": \"$database\"," .
  "\"id\": \"$tag_id\"," .
  "\"type\": \"$batch_type\"," .
  "\"snpstat\": \"$snp_sumstats\"," .
  "\"snpfst\": \"$snp_fst_vals\"," .
  "\"hapstat\": \"$snp_sumstats\"," .
  "\"hapfst\": \"$snp_fst_vals\",";

$json_str .= "\"snps\": [";
$rows = 0;
while ($row = $res->fetch_assoc()) {
    $json_str .= 
      "{" .
      "\"col\": \"$row[col]\"," .
      "\"rank_1\": \"$row[rank_1]\"," .
      "\"rank_2\": \"$row[rank_2]\"," .
      "\"rank_3\": \"$row[rank_3]\"," .
      "\"rank_4\": \"$row[rank_4]\"" .
      "},";
    $rows++;
}
if ($rows > 0)
    $json_str  = substr($json_str, 0, -1);
$json_str .= 
  "]," .
  "\"alleles\": [";

if (!$db['map_sth']->bind_param("ii", $batch_id, $tag_id))
    write_db_error($db['map_sth'], __FILE__, __LINE__);
if (!$db['map_sth']->execute())
    write_db_error($db['map_sth'], __FILE__, __LINE__);
$res = $db['map_sth']->get_result();

if ($res->num_rows > 0)
    $row = $res->fetch_assoc();

$rows = 0;
if ($res->num_rows > 0 && strlen($row['geno_map']) > 0) {
    $map = array();

    $genos = explode(";", $row['geno_map']);
    $i     = 0;
    foreach ($genos as $g) {
	if (strlen($g) == 0) continue;
	$m = explode(":", $g);
	$map[$m[0]] = $m[1];
	$alleles[$m[0]] = $colors[$i % $color_size];
	$i++;
    }
    asort($map);
    foreach ($map as $hapl => $geno) {
	$json_str .= 
	"{" .
	"\"gtype\": \"$geno\"," .
	"\"htype\": \"$hapl\"" .
	"},";
    }
    $rows++;

} else {
    if (!$db['all_sth']->bind_param("ii", $batch_id, $tag_id))
	write_db_error($db['all_sth'], __FILE__, __LINE__);
    if (!$db['all_sth']->execute())
	write_db_error($db['all_sth'], __FILE__, __LINE__);
    $res = $db['all_sth']->get_result();

    $gtypes = "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ";

    $i = 0;
    while ($row = $res->fetch_assoc()) {
        $json_str .= 
	  "{" .
	  "\"gtype\": \"" . $gtypes[$i % 52]  . "\"," .
	  "\"htype\": \"$row[allele]\"" .
	  "},";
	$i++;
	$rows++;
    }
}

if ($rows > 0)
    $json_str  = substr($json_str, 0, -1);
$json_str .= 
  "]," .
  "\"popkey\": {";

$htypes = array();
$gtypes = array();
//
// Fetch and record Observed Haplotypes
//
if (!$db['mat_sth']->bind_param("ii", $batch_id, $tag_id))
    write_db_error($db['mat_sth'], __FILE__, __LINE__);
if (!$db['mat_sth']->execute())
    write_db_error($db['mat_sth'], __FILE__, __LINE__);
$res = $db['mat_sth']->get_result();

while ($row = $res->fetch_assoc()) {
    $a = array('id'     => $row['id'],
               'file'   => $row['file'], 
               'allele' => $row['allele'], 
               'tag_id' => $row['tag_id'],
               'depth'  => $row['depth'],
	       'lnl'    => $row['lnl'],
	       'pop_id' => $row['pop_id']);
    if (!isset($htypes[$row['pop_id']]))
        $htypes[$row['pop_id']] = array();

    if (!isset($htypes[$row['pop_id']][$row['file']]))
        $htypes[$row['pop_id']][$row['file']] = array();

    array_push($htypes[$row['pop_id']][$row['file']], $a);
}

//
// Fetch and record Genotypes
//
if (!$db['geno_sth']->bind_param("ii", $batch_id, $tag_id))
    write_db_error($db['geno_sth'], __FILE__, __LINE__);
if (!$db['geno_sth']->execute())
    write_db_error($db['geno_sth'], __FILE__, __LINE__);
$res = $db['geno_sth']->get_result();

while ($row = $res->fetch_assoc()) {
    $gtypes[$row['file']] = array('id'        => $row['sample_id'],
				  'file'      => $row['file'],
				  'genotype'  => $row['genotype'],
				  'corrected' => $row['corrected'],
				  'marker'    => $row['marker']);
}

ksort($htypes);

//
// Print the population key.
//
foreach ($htypes as $pop_id => $population) {

	if (isset($pop_names[$pop_id]))
	  $json_str .=
	    "\"$pop_id\": \"$pop_names[$pop_id]\",";
	else
	  $json_str .=
	    "\"$pop_id\": \"\",";
}

$json_str  = substr($json_str, 0, -1);
$json_str .= 
  "}," .
  "\"populations\": {";

//
// Print the observed haplotypes grouped by population.
//
foreach ($htypes as $pop_id => $population) {

    $json_str .= 
      "\"$pop_id\": [";

    foreach ($population as $sample => $match) {

        $genotype  = "";
        $corrected = 0;
	$marker    = "";

	if (count($gtypes) > 0 && isset($gtypes[$sample])) {
	    $marker = $gtypes[$sample]['marker'];

	    if (strlen($gtypes[$sample]['corrected']) > 0) {
	        $genotype  = $gtypes[$sample]['corrected'];
		$corrected = 1;
	    } else {
	        $genotype  = $gtypes[$sample]['genotype'];
	    }
	}

        $json_str .= 
	  "{" .
	  "\"sample\": \"$sample\"," .
	  "\"sample_id\": \"" . $match[0]['id'] . "\"," .
	  "\"lnl\": \"" . $match[0]['lnl'] . "\",";

	if ($batch_type == "map") {
	    $json_str .=
	        "\"marker\": \"$marker\"," .
	        "\"genotype\": \"$genotype\"," .
	        "\"corrected\": \"$corrected\",";
	}
	$json_str .=
	    "\"obshap\": [";

	foreach ($match as $m) {
	    $json_str .= 
	      "{" .
	      "\"tag_id\": \"$m[tag_id]\"," .
	      "\"allele\": \"$m[allele]\"," .
	      "\"depth\": \"$m[depth]\"" .
	      "},";
	}

	$json_str  = substr($json_str, 0, -1);
	$json_str .= 
	  "]},";
    }

    $json_str  = substr($json_str, 0, -1);
    $json_str .= 
      "],";
}

$json_str  = substr($json_str, 0, -1);
$json_str .= "}}";

echo $json_str;

?>
