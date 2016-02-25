<?php
//
// Copyright 2010-2016, Julian Catchen <jcatchen@illinois.edu>
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

$database  = isset($_GET['db'])       ? $_GET['db']       : "";
$tag_id    = isset($_GET['tag_id'])   ? $_GET['tag_id']   : 0;
$batch_id  = isset($_GET['batch_id']) ? $_GET['batch_id'] : 0;
$ext_id    = isset($_GET['ext_id'])   ? $_GET['ext_id']   : "";

// Connect to the database
if (!isset($db)) $db = db_connect($database);

// Save these variables for automatic URL formation later on.
$display = array();
$display['db']       = $database;
$display['tag_id']   = $tag_id;
$display['batch_id'] = $batch_id;
$display['ext_id']   = $ext_id;

//
// Prepare some SQL queries
//
$query = 
    "SELECT id, external_id " . 
    "FROM catalog_annotations " . 
    "WHERE batch_id=? and catalog_id=?";
if (!($db['sel_sth'] = $db['dbh']->prepare($query)))
    write_db_error($db['dbh'], __FILE__, __LINE__);

$query = 
    "UPDATE catalog_annotations SET external_id=? WHERE id=?";
if (!($db['upd_sth'] = $db['dbh']->prepare($query)))
    write_db_error($db['dbh'], __FILE__, __LINE__);

$query = 
    "INSERT INTO catalog_annotations SET batch_id=?, catalog_id=?, external_id=?";
if (!($db['ins_sth'] = $db['dbh']->prepare($query)))
    write_db_error($db['dbh'], __FILE__, __LINE__);

$query = 
    "DELETE FROM catalog_annotations WHERE id=?";
if (!($db['del_sth'] = $db['dbh']->prepare($query)))
    write_db_error($db['dbh'], __FILE__, __LINE__);

//
// Fetch any existing annotation for this marker
//
if (!$db['sel_sth']->bind_param("ii", $display['batch_id'], $display['tag_id']))
    write_db_error($db['sel_sth'], __FILE__, __LINE__);
if (!$db['sel_sth']->execute())
    write_db_error($db['sel_sth'], __FILE__, __LINE__);
$res = $db['sel_sth']->get_result();

$external_id = "";
$sql_id      = 0;

if ($row = $res->fetch_assoc()) {
    $external_id = $row['external_id'];
    $sql_id      = $row['id'];
}

if ($external_id != $ext_id) {
    //
    // Is this annotation being reset to the original value? If so, delete the corrected record.
    //
    if (strlen($external_id) > 0 && strlen($ext_id) == 0) {
        if (!$db['del_sth']->bind_param("i", $sql_id))
            write_db_error($db['del_sth'], __FILE__, __LINE__);
        if (!$db['del_sth']->execute())
            write_db_error($db['del_sth'], __FILE__, __LINE__);
    //
    // Are we changing an existing annotation?
    //
    } else if (strlen($external_id) > 0 && strlen($ext_id) > 0) {
        if (!$db['upd_sth']->bind_param("si", $ext_id, $sql_id))
            write_db_error($db['upd_sth'], __FILE__, __LINE__);
        if (!$db['upd_sth']->execute())
            write_db_error($db['upd_sth'], __FILE__, __LINE__);
    //
    // Otherwise, add a new annotation.
    //
    } else if (strlen($ext_id) > 0) {
        if (!$db['ins_sth']->bind_param("iis", $display['batch_id'], $display['tag_id'], $ext_id))
            write_db_error($db['ins_sth'], __FILE__, __LINE__);
        if (!$db['ins_sth']->execute())
            write_db_error($db['ins_sth'], __FILE__, __LINE__);
    }
}

header("Content-type: text/xml");
$xml_output = 
    "<?xml version=\"1.0\"?>\n" .
    "<annotation>\n" . 
    "<text>$ext_id</text>\n" .
    "<marker_id>$tag_id</marker_id>\n" .
    "</annotation>\n";

echo $xml_output;

?>
