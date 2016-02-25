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
    "SELECT pop_id, pop_name FROM populations " . 
    "WHERE batch_id=?";
if (!($db['pop_sth'] = $db['dbh']->prepare($query)))
    write_db_error($db['dbh'], __FILE__, __LINE__);

$query = 
    "SELECT pop_id_1, pop_id_2, phist, fpst FROM phist " . 
    "WHERE batch_id=? AND tag_id=?";
if (!($db['fst_sth'] = $db['dbh']->prepare($query)))
    write_db_error($db['dbh'], __FILE__, __LINE__);

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

if (!$db['fst_sth']->bind_param("ii", $batch_id, $tag_id))
    write_db_error($db['fst_sth'], __FILE__, __LINE__);
if (!$db['fst_sth']->execute())
    write_db_error($db['fst_sth'], __FILE__, __LINE__);
$res = $db['fst_sth']->get_result();

$stats = array();
$pops  = array();

while ($row = $res->fetch_assoc()) {
    $a = array('pid_1'  => $row['pop_id_1'],
	       'pid_2'  => $row['pop_id_2'],
	       'phist'  => $row['phist'],
	       'fpst'   => $row['fpst']);

    array_push($stats, $a);

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
  "\"phist\": [";

foreach ($stats as $s) {

    $phist = $s['phist'] != 0 ? sprintf("%.3f", $s['phist']) : "0";
    $fpst  = $s['fpst']  != 0 ? sprintf("%.3f", $s['fpst'])  : "0";

    $json_str .=
      "{" .
      "\"pid_1\": \"$s[pid_1]\"," .
      "\"pid_2\": \"$s[pid_2]\"," .
      "\"phist\": \"$phist\"," .
      "\"fstp\": \"$fpst\"" .
      "},";
}

if (count($stats) > 0) 
    $json_str  = substr($json_str, 0, -1);
$json_str .= 
  "]}";

echo $json_str;

?>
