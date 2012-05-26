<?php
//
// Copyright 2012, Julian Catchen <jcatchen@uoregon.edu>
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

$database     = isset($_GET['db'])       ? $_GET['db']       : "";
$new_pop_name = isset($_GET['pop_name']) ? $_GET['pop_name'] : "";
$batch_id     = isset($_GET['batch_id']) ? $_GET['batch_id'] : 0;
$pop_id       = isset($_GET['pop_id'])   ? $_GET['pop_id']   : "";

// Connect to the database
if (!isset($db)) $db = db_connect($database);

// Save these variables for automatic URL formation later on.
$display = array();
$display['db']       = $database;
$display['tag_id']   = $tag_id;
$display['batch_id'] = $batch_id;
$display['pop_id']   = $pop_id;

//
// Prepare some SQL queries
//
$query = 
    "SELECT id, pop_name " . 
    "FROM populations " . 
    "WHERE batch_id=? and pop_id=?";
$db['sel_sth'] = $db['dbh']->prepare($query);
check_db_error($db['sel_sth'], __FILE__, __LINE__);

$query = 
    "UPDATE populations SET pop_name=? WHERE id=?";
$db['upd_sth'] = $db['dbh']->prepare($query);
check_db_error($db['upd_sth'], __FILE__, __LINE__);

$query = 
    "INSERT INTO populations SET batch_id=?, pop_id=?, pop_name=?";
$db['ins_sth'] = $db['dbh']->prepare($query);
check_db_error($db['ins_sth'], __FILE__, __LINE__);

$query = 
    "DELETE FROM populations WHERE id=?";
$db['del_sth'] = $db['dbh']->prepare($query);
check_db_error($db['del_sth'], __FILE__, __LINE__);

//
// Fetch any existing annotation for this marker
//
$result = $db['sel_sth']->execute(array($display['batch_id'], $display['pop_id']));
check_db_error($result, __FILE__, __LINE__);

$pop_name = "";
$sql_id   = 0;

if ($row = $result->fetchRow()) {
    $pop_name = $row['pop_name'];
    $sql_id   = $row['id'];
}

if ($new_pop_name != $pop_name) {
    //
    // Is this annotation being reset to the original value? If so, delete the corrected record.
    //
    if (strlen($pop_name) > 0 && strlen($new_pop_name) == 0) {
        $result = $db['del_sth']->execute($sql_id);
        check_db_error($result, __FILE__, __LINE__);
    //
    // Are we changing an existing annotation?
    //
    } else if (strlen($pop_name) > 0 && strlen($new_pop_name) > 0) {
        $result = $db['upd_sth']->execute(array($new_pop_name, $sql_id));
        check_db_error($result, __FILE__, __LINE__);
    //
    // Otherwise, add a new annotation.
    //
    } else if (strlen($new_pop_name) > 0) {
        $result = $db['ins_sth']->execute(array($display['batch_id'], $display['pop_id'], $new_pop_name));
        check_db_error($result, __FILE__, __LINE__);
    }
}

header("Content-type: text/xml");
$xml_output = 
    "<?xml version=\"1.0\"?>\n" .
    "<annotation>\n" . 
    "<text>$new_pop_name</text>\n" .
    "<pop_id>$pop_id</pop_id>\n" .
    "</annotation>\n";

echo $xml_output;

?>
