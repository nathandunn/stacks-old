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

$database  = isset($_GET['db'])       ? $_GET['db']       : "";
$tag_id    = isset($_GET['tag_id'])   ? $_GET['tag_id']   : 0;
$batch_id  = isset($_GET['batch_id']) ? $_GET['batch_id'] : 0;
$page      = isset($_GET['p'])        ? $_GET['p']        : 1;
$per_page  = isset($_GET['pp'])       ? $_GET['pp']       : 10;

// Connect to the database
$db = db_connect($database);

// Save these variables for automatic URL formation later on.
$display = array();
$display['db']       = $database;
$display['tag_id']   = $tag_id;
$display['batch_id'] = $batch_id;
$display['p']        = $page;
$display['pp']       = $per_page;

//
// Prepare some SQL queries
//
$query = 
    "SELECT id, catalog_id, seq_id, type FROM sequence " . 
    "WHERE batch_id=? AND catalog_id=?";
$db['seq_sth'] = $db['dbh']->prepare($query);
check_db_error($db['seq_sth'], __FILE__, __LINE__);

$query = 
    "SELECT id, catalog_id, seq_id, algorithm, query_id, query_len, hit_id, hit_len, score, e_value, percent_ident, hsp_rank, " . 
    "aln_len, query_aln_start, query_aln_end, hit_aln_start, hit_aln_end " . 
    "FROM sequence_blast " . 
    "WHERE batch_id=? AND catalog_id=?";
$db['blast_sth'] = $db['dbh']->prepare($query);
check_db_error($db['blast_sth'], __FILE__, __LINE__);

$page_title = "Catalog RAD-Tag Sequence/BLAST Hits Viewer";
write_compact_header($page_title);

$result = $db['seq_sth']->execute(array($batch_id, $tag_id));
check_db_error($result, __FILE__, __LINE__);

$seqs = array();

// Add the marker itself to the sequence array, it won't
// be directly listed in the sequence table.
$a = array('id'         => -1,
           'catalog_id' => $tag_id,
           'seq_id'     => "",
           'type'       => "se_radtag");
array_push($seqs, $a);

while ($row = $result->fetchRow()) {
    $a = array('id'         => $row['id'],
               'catalog_id' => $row['catalog_id'],
               'seq_id'     => $row['seq_id'],
               'type'       => $row['type']);

    array_push($seqs, $a);
}

$result = $db['blast_sth']->execute(array($batch_id, $tag_id));
check_db_error($result, __FILE__, __LINE__);

$hits = array();

while ($row = $result->fetchRow()) {
    $a = array('sql_id'    => $row['id'],
               'query_id'  => $row['query_id'],
               'hit_id'    => $row['hit_id'],
               'query_len' => $row['query_len'],
               'hit_len'   => $row['hit_len'],
               'score'     => $row['score'],
               'e_value'   => $row['e_value'],
               'pctid'     => $row['percent_ident'],
               'hsp_rank'  => $row['hsp_rank'],
               'aln_len'   => $row['aln_len']);
    if (!isset($hits[$row['query_id']]))
        $hits[$row['query_id']] = array();

    array_push($hits[$row['query_id']], $a);
}

foreach ($seqs as $seq) {
    if (strlen($seq['seq_id']) > 0)
        $query_id = $seq['catalog_id'] . "|" . $seq['seq_id'];
    else
        $query_id = $seq['catalog_id'];
    $hsps = $hits[$query_id];
    $i    = 0;

    if ($seq['type'] == "se_radtag" && !isset($hsps))
      continue;
    else if ($seq['type'] == "se_radtag")
      print "<h4><img src=\"$img_path/caret-u.png\" />$seq[seq_id] [$seq[type]]</h4>\n";
    else
      echo <<< EOQ
<h4><img src="$img_path/caret-u.png" />
    <a target="_blank" href="$root_path/view_sequence.php?db=$database&id=$seq[id]">$seq[seq_id] [$seq[type]]</a>
</h4>

EOQ;

    if (!isset($hsps)) {
        print "</table>\n";
	continue;
    }

    echo <<< EOQ
<table class="catalog">
<tr>
  <th>&nbsp;</th>
  <th>Hit ID</th>
  <th>Score</th>
  <th>E-Value</th>
  <th>Percent Ident</th>
  <th>Aln Len</th>
  <th>Query Len</th>
  <th>Hit Len</th>
  <th>HSP Rank</th>
</tr>

EOQ;

    foreach ($hsps as $hsp) {
        $i++;
        echo <<< EOQ
<tr>
<td>$i</td>
<td>$hsp[hit_id]</td>
<td>$hsp[score]</td>
<td>$hsp[e_value]</td>
<td>$hsp[pctid]</td>
<td>$hsp[aln_len]</td>
<td>$hsp[query_len]</td>
<td>$hsp[hit_len]</td>
<td>$hsp[hsp_rank]</td>
</tr>

EOQ;
    }

    print "</table>\n";
}

write_compact_footer();

?>
