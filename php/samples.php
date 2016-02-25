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

$batch_id = isset($_GET['id']) ? $_GET['id'] : 0;
$database = isset($_GET['db']) ? $_GET['db'] : "";

// Connect to the database
$db = db_connect($database);

//
// Prepare some SQL queries
//
$query = 
    "SELECT id, date, description FROM batches WHERE id=?";
if (!($db['batch_sth'] = $db['dbh']->prepare($query)))
    write_db_error($db['dbh'], __FILE__, __LINE__);

$query = 
    "SELECT id, sample_id, type, file FROM samples " . 
    "WHERE batch_id=? ORDER BY id";
if (!($db['samp_sth'] = $db['dbh']->prepare($query)))
    write_db_error($db['dbh'], __FILE__, __LINE__);

$query = 
    "SELECT COUNT(tag_id) AS count FROM tag_index " . 
    "WHERE batch_id=? AND sample_id=?";
if (!($db['count_sth'] = $db['dbh']->prepare($query)))
    write_db_error($db['dbh'], __FILE__, __LINE__);

$query = 
    "SELECT snps FROM tag_index " .
    "WHERE batch_id=? AND sample_id=? AND snps>0";
if (!($db['snp_sth'] = $db['dbh']->prepare($query)))
    write_db_error($db['dbh'], __FILE__, __LINE__);

//
// Pull information about this batch
//
if (!$db['batch_sth']->bind_param("i", $batch_id))
    write_db_error($db['batch_sth'], __FILE__, __LINE__);
if (!$db['batch_sth']->execute())
    write_db_error($db['batch_sth'], __FILE__, __LINE__);
$res = $db['batch_sth']->get_result();
$row = $res->fetch_assoc();

$batch  = array();
$batch['id']   = $row['id'];
$batch['desc'] = $row['description'];
$batch['date'] = $row['date'];

$page_title = "RAD-Tag Samples";
write_header($page_title);

echo <<< EOQ
<form id="page_state" name="page_state"></form>

<h4 class="info_head" style="margin-left: 1em;">
  <img id="sources_img" src="$img_path/caret-d.png" />
  <a href="$root_path/index.php?db=$database">
  RAD-Tag Samples for batch #$batch[id] <span class="s">[$batch[date]; $batch[desc]]</span></a>
</h4>

<div id="sources" style="width: 100%; margin-top: 2em;">
<a name="results_top"></a>
<table class="db" style="width: 90%;">
<tr>
  <th style="width: 10%;">Id</th>
  <th style="width: 20%;">Type</th>
  <th style="width: 16%;">Unique Stacks</th>
  <th style="width: 12%;">Polymorphic Loci</th>
  <th style="width: 12%;">SNPs Found</th>
  <th style="width: 30%;">Source</th>
</tr>

EOQ;

if (!$db['samp_sth']->bind_param("i", $batch_id))
    write_db_error($db['samp_sth'], __FILE__, __LINE__);
if (!$db['samp_sth']->execute())
    write_db_error($db['samp_sth'], __FILE__, __LINE__);
$res = $db['samp_sth']->get_result();

$tag_id = 0;
if (!$db['count_sth']->bind_param("ii", $batch_id, $tag_id))
    write_db_error($db['count_sth'], __FILE__, __LINE__);
if (!$db['snp_sth']->bind_param("ii", $batch_id, $tag_id))
    write_db_error($db['count_sth'], __FILE__, __LINE__);

while ($row = $res->fetch_assoc()) {
    $snps   = 0;
    $poly   = 0;
    $tag_id = $row['id'];
    
    //
    // Query the database to determine how many tags belong to this sample.
    //
    if (!$db['count_sth']->execute())
        write_db_error($db['count_sth'], __FILE__, __LINE__);
    $count_res = $db['count_sth']->get_result();
    $count_row = $count_res->fetch_assoc();
    $count = $count_row['count'];

    //
    // Query the database to find how many SNPs were found in this sample.
    //
    if (!$db['snp_sth']->execute())
        write_db_error($db['snp_sth'], __FILE__, __LINE__);
    $count_res = $db['snp_sth']->get_result();

    while ($count_row = $count_res->fetch_assoc()) {
      $snps += $count_row['snps'];
      $poly += $count_row['snps'] > 0 ? 1 : 0;
    }

    print
	"<tr>\n" .
	"  <td><a href=\"$root_path/tags.php?db=$database&batch_id=$batch[id]&sample_id=$row[id]\">$row[id]</a></td>\n" .
	"  <td>" . ucfirst($row['type']) . "</td>\n" .
	"  <td>" . $count . "</td>\n" .
	"  <td>" . $poly . "</td>\n" .
	"  <td>" . $snps . "</td>\n" .
	"  <td>" . $row['file'] . "</td>\n" .
	"</tr>\n";
}

echo <<< EOQ
</table>
</div>

EOQ;

write_footer();

?>
