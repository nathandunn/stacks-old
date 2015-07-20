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
    "SELECT col, pop_id_1, pop_id_2, pi_o, amova_fst_c as fst, fishers_p, lod FROM fst " . 
    "WHERE batch_id=? AND tag_id=?";
$db['fst_sth'] = $db['dbh']->prepare($query);
check_db_error($db['fst_sth'], __FILE__, __LINE__);

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

$result = $db['fst_sth']->execute(array($batch_id, $tag_id));
check_db_error($result, __FILE__, __LINE__);

$stats = array();
$pops  = array();

while ($row = $result->fetchRow()) {
  $a = array('col'       => $row['col'],
	     'pid_1'     => $row['pop_id_1'],
	     'pid_2'     => $row['pop_id_2'],
	     'pi_o'      => $row['pi_o'],
	     'fishers_p' => $row['fishers_p'],
	     'lod'       => $row['lod'],
	     'fst'       => $row['fst']);

  if (!isset($stats[$row['col']]))
    $stats[$row['col']] = array();

  array_push($stats[$row['col']], $a);

  $pops[$row['pop_id_1']] = $row['pop_id_1'];
  $pops[$row['pop_id_2']] = $row['pop_id_2'];
}

ksort($stats);
ksort($pops);

//
// Assign population IDs for any missing population names.
//
foreach ($pops as $pop_id)
  if (!isset($pop_names[$pop_id])) 
    $pop_names[$pop_id] = $pop_id;

ksort($pop_names);

$json_str = 
  "{" .
  "\"path\": \"$root_path\"," .
  "\"batch_id\": \"$batch_id\"," .
  "\"db\": \"$database\"," .
  "\"id\": \"$tag_id\"," .
  "\"type\": \"$batch_type\",";

$json_str .=  "\"popkey\": {";
//
// Print the population key.
//
foreach ($pop_names as $pop_id => $population) {
  $json_str .=
    "\"$pop_id\": \"$population\",";
}

$json_str  = substr($json_str, 0, -1);
$json_str .= 
  "}," .
  "\"columns\": {";

foreach ($stats as $col => $stat) {

  $json_str .= "\"$col\": [";

  foreach ($stat as $s) {

    $fst = $s['fst']  != 0 ? sprintf("%.3f", $s['fst'])  : "0";
    $lod = $s['lod']  != 0 ? sprintf("%.3f", $s['lod'])  : "0";
    $pio = $s['pi_o'] != 0 ? sprintf("%.3f", $s['pi_o']) : "0";
    $p   = $s['fishers_p'] != 0 ? sprintf("%.3f", $s['fishers_p']) : "0";

    $json_str .=
      "{" .
      "\"pid_1\": \"$s[pid_1]\"," .
      "\"pid_2\": \"$s[pid_2]\"," .
      "\"pi_o\": \"$pio\"," .
      "\"p\": \"$p\"," .
      "\"lod\": \"$lod\"," .
      "\"fst\": \"$fst\"" .
      "},";
  }
  
  if (count($stat) > 0)
    $json_str  = substr($json_str, 0, -1);
  $json_str .= "],";
}

if (count($stats) > 0) 
  $json_str  = substr($json_str, 0, -1);
$json_str .= 
  "}}";

echo $json_str;

?>
