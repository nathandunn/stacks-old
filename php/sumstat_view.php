<?php
//
// Copyright 2015, Julian Catchen <jcatchen@illinois.edu>
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
    "SELECT pop_id, pop_name FROM populations " . 
    "WHERE batch_id=?";
$db['pop_sth'] = $db['dbh']->prepare($query);
check_db_error($db['pop_sth'], __FILE__, __LINE__);

$query = 
   "SELECT pop_id, col, bp, p_nuc, q_nuc, n, p, obs_het, obs_hom, exp_het, exp_hom, pi, fis FROM sumstats " . 
   "WHERE batch_id=? AND tag_id=?";
$db['stats_sth'] = $db['dbh']->prepare($query);
check_db_error($db['stats_sth'], __FILE__, __LINE__);

//
// Fetch population names if available.
//
$pop_names = array();
if ($batch_type == "population") {
    $result = $db['pop_sth']->execute($batch_id);
    check_db_error($result, __FILE__, __LINE__);

    while ($row = $result->fetchRow())
        $pop_names[$row['pop_id']] = $row['pop_name'];
}

$result = $db['stats_sth']->execute(array($batch_id, $tag_id));
check_db_error($result, __FILE__, __LINE__);

$stats = array();
$pops  = array();

while ($row = $result->fetchRow()) {
  $a = array('col'     => $row['col'],
	     'bp'      => $row['bp'],
	     'p_nuc'   => $row['p_nuc'],
	     'q_nuc'   => $row['q_nuc'],
	     'n'       => $row['n'],
	     'p'       => $row['p'],
	     'obs_het' => $row['obs_het'],
	     'obs_hom' => $row['obs_hom'],
	     'exp_het' => $row['exp_het'],
	     'exp_hom' => $row['exp_hom'],
	     'pi'      => $row['pi'],
	     'fis'     => $row['fis'],
	     'pop_id'  => $row['pop_id']);

  if (!isset($stats[$row['col']]))
    $stats[$row['col']] = array();

  $stats[$row['col']][$row['pop_id']] = $a;

  $pops[$row['pop_id']] = $row['pop_id'];
}

ksort($stats);
ksort($pops);

$json_str = 
  "{" .
  "\"path\": \"$root_path\"," .
  "\"batch_id\": \"$batch_id\"," .
  "\"db\": \"$database\"," .
  "\"id\": \"$tag_id\"," .
  "\"type\": \"$batch_type\",";

$json_str .= "\"sumstats\": [";

foreach ($pops as $pop_id)
  if (!isset($pop_names[$pop_id])) 
    $pop_names[$pop_id] = $pop_id;

$index = count($pop_names) == 0 ? $pops : $pop_names;

$rows = 0;
foreach ($stats as $col => $stat) {    

  foreach ($index as $pop_id => $pop_name) {
    if (!isset($stat[$pop_id])) 
      continue;

    $s    = $stat[$pop_id];
    $p    = $s['p']       < 1 ? sprintf("%.5f", $s['p']) : $s['p'];
    $ohet = $s['obs_het'] > 0 ? sprintf("%.3f", $s['obs_het']) : $s['obs_het'];
    $ohom = $s['obs_hom'] < 1 ? sprintf("%.3f", $s['obs_hom']) : $s['obs_hom'];
    $ehet = $s['exp_het'] > 0 ? sprintf("%.3f", $s['exp_het']) : $s['exp_het'];
    $ehom = $s['exp_hom'] < 1 ? sprintf("%.3f", $s['exp_hom']) : $s['exp_hom'];
    $pi   = $s['pi']      > 0 ? sprintf("%.3f", $s['pi']) : $s['pi'];
    $fis  = $s['fis']    != 0 ? sprintf("%.3f", $s['fis']) : "0";

    $json_str .=
      "{" .
      "\"pop_id\": \"$pop_name\"," .
      "\"col\": \"$s[col]\"," .
      "\"bp\": \"$s[bp]\"," .
      "\"p_allele\": \"$s[p_nuc]\"," .
      "\"q_allele\": \"$s[q_nuc]\"," .
      "\"p\": \"$p\"," .
      "\"n\": \"$s[n]\"," .
      "\"obshet\": \"$ohet\"," .
      "\"obshom\": \"$ohom\"," .
      "\"exphet\": \"$ehet\"," .
      "\"exphom\": \"$ehom\"," .
      "\"pi\": \"$pi\"," .
      "\"fis\": \"$fis\"" .
      "},";
    $rows++;
  }
}

if ($rows > 0)
  $json_str  = substr($json_str, 0, -1);
$json_str .= 
  "]}";

echo $json_str;

?>
