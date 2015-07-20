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
   "SELECT pop_id, bp, n, hapcnt, gene_div, hap_div FROM hapstats " . 
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

while ($row = $result->fetchRow()) {
  $a = array('bp'       => $row['bp'],
	     'n'        => $row['n'],
	     'hapcnt'   => $row['hapcnt'],
	     'gene_div' => $row['gene_div'],
	     'hap_div'  => $row['hap_div'],
	     'pop_id'   => $row['pop_id']);

  $stats[$row['pop_id']] = $a;
}

ksort($stats);

$json_str = 
  "{" .
  "\"path\": \"$root_path\"," .
  "\"batch_id\": \"$batch_id\"," .
  "\"db\": \"$database\"," .
  "\"id\": \"$tag_id\"," .
  "\"type\": \"$batch_type\",";

$json_str .= "\"hapstats\": [";

foreach ($stats as $pop_id => $stat)
  if (!isset($pop_names[$pop_id])) 
    $pop_names[$pop_id] = $pop_id;

$rows = 0;
foreach ($stats as $pop_id => $s) {

    $gdiv = $s['gene_div'] > 0 ? sprintf("%.3f", $s['gene_div']) : $s['gene_div'];
    $hdiv = $s['hap_div']  > 0 ? sprintf("%.3f", $s['hap_div'])  : $s['hap_div'];

    $json_str .=
      "{" .
      "\"pop_id\": \"" . $pop_names[$pop_id] . "\"," .
      "\"bp\": \"$s[bp]\"," .
      "\"n\": \"$s[n]\"," .
      "\"hapcnt\": \"$s[hapcnt]\"," .
      "\"genediv\": \"$gdiv\"," .
      "\"hapdiv\": \"$hdiv\"" .
      "},";
    $rows++;
}

if ($rows > 0) 
  $json_str  = substr($json_str, 0, -1);
$json_str .= 
  "]}";

echo $json_str;

?>
