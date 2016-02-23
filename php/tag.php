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

$tag_id    = isset($_GET['tag_id'])    ? $_GET['tag_id']    : 0;
$batch_id  = isset($_GET['batch_id'])  ? $_GET['batch_id']  : 0;
$sample_id = isset($_GET['sample_id']) ? $_GET['sample_id'] : 0;
$database  = isset($_GET['db'])        ? $_GET['db']        : "";
$page      = isset($_GET['p'])         ? $_GET['p']         : 1;
$per_page  = isset($_GET['pp'])        ? $_GET['pp']        : 10;

// Connect to the database
$db = db_connect($database);

// Save these variables for automatic URL formation later on.
$display = array();
$display['tag_id']    = $tag_id;
$display['batch_id']  = $batch_id;
$display['sample_id'] = $sample_id;
$display['db']        = $database;
$display['p']         = $page;
$display['pp']        = $per_page;

//
// Prepare some SQL queries
//
$query = 
    "SELECT batches.id as id, date, description, samples.id as sample_id, samples.sample_id as samp_id, samples.type, file " . 
    "FROM batches " . 
    "JOIN samples ON (batch_id=batches.id) " . 
    "WHERE batches.id=? AND samples.id=?";
if (!($db['batch_sth'] = $db['dbh']->prepare($query)))
    write_db_error($db['dbh'], __FILE__, __LINE__);

$query = 
    "SELECT count(id) as depth FROM unique_tags " . 
    "WHERE relationship!='consensus' AND relationship!='model' AND sample_id=? AND tag_id=?";
if (!($db['depth_sth'] = $db['dbh']->prepare($query)))
    write_db_error($db['dbh'], __FILE__, __LINE__);

$query = 
    "SELECT tag_id, sub_id, relationship, seq_id, seq, deleveraged, blacklisted, removed " . 
    "FROM unique_tags " . 
    "JOIN samples ON (unique_tags.sample_id=samples.id) " . 
    "JOIN batches ON (samples.batch_id=batches.id) " .
    "WHERE batch_id=? AND unique_tags.sample_id=? AND tag_id=?";
if (!($db['seq_sth'] = $db['dbh']->prepare($query)))
    write_db_error($db['dbh'], __FILE__, __LINE__);

$query = 
    "SELECT col, rank_1, rank_2 FROM snps " . 
    "JOIN samples ON (snps.sample_id=samples.id) " . 
    "JOIN batches ON (samples.batch_id=batches.id) " .
    "WHERE batch_id=? AND snps.sample_id=? AND tag_id=? AND snps.type='E' ORDER BY col";
if (!($db['snp_sth'] = $db['dbh']->prepare($query)))
    write_db_error($db['dbh'], __FILE__, __LINE__);

$query = 
    "SELECT allele, read_pct FROM alleles " . 
    "WHERE sample_id=? AND tag_id=? ";
if (!($db['all_sth'] = $db['dbh']->prepare($query)))
    write_db_error($db['dbh'], __FILE__, __LINE__);

$query = 
    "SELECT catalog_id FROM matches " . 
    "WHERE batch_id=? AND sample_id=? AND tag_id=? " .
    "GROUP BY catalog_id";
if (!($db['cat_sth'] = $db['dbh']->prepare($query)))
    write_db_error($db['dbh'], __FILE__, __LINE__);

$query = 
    "SELECT col FROM catalog_snps " . 
    "WHERE batch_id=? AND tag_id=?";
if (!($db['cat_snp_sth'] = $db['dbh']->prepare($query)))
    write_db_error($db['dbh'], __FILE__, __LINE__);

//
// Pull information about this batch
//
if (!$db['batch_sth']->bind_param("ii", $batch_id, $sample_id))
    write_db_error($db['batch_sth'], __FILE__, __LINE__);
if (!$db['batch_sth']->execute())
    write_db_error($db['batch_sth'], __FILE__, __LINE__);
$res = $db['batch_sth']->get_result();
$row = $res->fetch_assoc();

$batch                = array();
$batch['id']          = $row['id'];
$batch['desc']        = $row['description'];
$batch['date']        = $row['date'];
$batch['file']        = $row['file'];
$batch['sample_id']   = $row['sample_id'];
$batch['samp_id']     = $row['samp_id'];
$batch['sample_type'] = $row['type'];

$page_title = "RAD-Tag Viewer";
write_header($page_title, $batch);

echo <<< EOQ
<form id="page_state" name="page_state"></form>
<h3 style="margin-left: 1em;">
<a href="$root_path/samples.php?db=$database&id=$batch[id]">
Batch #$batch[id] <span class="s">[$batch[date]; $batch[desc]]</span></a></h3>

<h4 style="margin-left: 1em;">
<a href="$root_path/tags.php?db=$database&batch_id=$batch[id]&sample_id=$batch[sample_id]&p=$display[p]&pp=$display[pp]">
  RAD-Tag Sample #$batch[samp_id] [<span class="s">$batch[file]</span>]</a>
</h4>

<h4 class="info_head" style="margin-left: 1em;">
  <img id="sources_img" src="$img_path/caret-d.png" />
  <a onclick="toggle_div('sources', '$img_path', 'page_state');">
  Sequence #$display[tag_id]</a>
</h4>

<div id="sources" style="width: 100%;">
<table class="radtag" style="width: 90%; margin-left: 2.5%; margin-right: auto; margin-bottom: 1em;">
<tr>
  <th style="width: 10%;">Catalog ID</th>
  <th style="width: 5%;">Depth</th>
  <th style="width: 25%;">SNPs</th>
  <th style="width: 25%;">Alleles</th>
  <th style="width: 10%;">Deleveraged?</th>
  <th style="width: 15%;">Lumberjackstack?</th>
  <th style="width: 10%;">Blacklisted?</th>
</tr>
<tr>

EOQ;

//
// Fetch the seqeunces
//
$deleveraged     = 0;
$lumberjackstack = 0;
$blacklisted     = 0;
$seqs = array('consensus' => array(),
	      'primary'   => array(),
	      'secondary' => array(),
	      'model'     => array());

if (!$db['seq_sth']->bind_param("iii", $batch_id, $sample_id, $tag_id))
    write_db_error($db['seq_sth'], __FILE__, __LINE__);
if (!$db['seq_sth']->execute())
    write_db_error($db['seq_sth'], __FILE__, __LINE__);
$res = $db['seq_sth']->get_result();

while ($row = $res->fetch_assoc()) {
    array_push($seqs[$row['relationship']], array('s' => $row['seq'], 'id' => $row['seq_id'], 'sub_id' => $row['sub_id']));

    if ($row['relationship'] == "consensus") {
	$deleveraged     = $row['deleveraged'];
	$lumberjackstack = $row['removed'];
	$blacklisted     = $row['blacklisted'];
    }
}

if (!$db['cat_sth']->bind_param("iii", $batch_id, $sample_id, $tag_id))
    write_db_error($db['cat_sth'], __FILE__, __LINE__);
if (!$db['cat_sth']->execute())
    write_db_error($db['cat_sth'], __FILE__, __LINE__);
$res = $db['cat_sth']->get_result();
$row = $res->fetch_assoc();

if (isset($row['catalog_id'])) 
    $catalog_id = $row['catalog_id'];
else
    $catalog_id = -1;

print "<td style=\"text-align: center; vertical-align: top;\">\n";
if ($catalog_id >= 0)
  print "  <a href=\"$root_path/catalog.php?db=$database&id=$batch_id&filter_type[]=cata&filter_cata=$catalog_id\">#$catalog_id</a>\n";
print "</td>\n";

if (!$db['depth_sth']->bind_param("ii", $sample_id, $tag_id))
    write_db_error($db['depth_sth'], __FILE__, __LINE__);
if (!$db['depth_sth']->execute())
    write_db_error($db['depth_sth'], __FILE__, __LINE__);
$res = $db['depth_sth']->get_result();
$row = $res->fetch_assoc();

echo <<< EOQ
<td style="text-align: center; vertical-align: top;">
  $row[depth]x
</td>

EOQ;

$snps   = array();
if (!$db['snp_sth']->bind_param("iii", $batch_id, $sample_id, $tag_id))
    write_db_error($db['snp_sth'], __FILE__, __LINE__);
if (!$db['snp_sth']->execute())
    write_db_error($db['snp_sth'], __FILE__, __LINE__);
$res = $db['snp_sth']->get_result();

while ($row = $res->fetch_assoc()) {
  $snps[$row['col']] = array('col' => $row['col'], 'rank_1' => $row['rank_1'], 'rank_2' => $row['rank_2']);
}

$cat_snps = array();
if ($catalog_id >= 0) {

    if (!$db['cat_snp_sth']->bind_param("ii", $batch_id, $catalog_id))
	write_db_error($db['cat_snp_sth'], __FILE__, __LINE__);
    if (!$db['cat_snp_sth']->execute())
	write_db_error($db['cat_snp_sth'], __FILE__, __LINE__);
    $res = $db['cat_snp_sth']->get_result();

    while ($row = $res->fetch_assoc()) {
	if (!isset($cat_snps[$row['col']]))
  	    $cat_snps[$row['col']] = 0;
	$cat_snps[$row['col']]++;
    }
}

print 
"  <td style=\"vertical-align: top;\">\n".
"    <table style=\"margin-left: auto; margin-right: auto; width: 60%;\">\n";

foreach ($snps as $snp) {
    print 
      "<tr>\n" .
      "  <td>Column: $snp[col]</td><td style=\"text-align: right;\">$snp[rank_1]/$snp[rank_2]</td>\n" .
      "</tr>\n";
}

print 
"  </table>\n" .
"  </td>\n";

if (!$db['all_sth']->bind_param("ii", $sample_id, $tag_id))
    write_db_error($db['all_sth'], __FILE__, __LINE__);
if (!$db['all_sth']->execute())
    write_db_error($db['all_sth'], __FILE__, __LINE__);
$res = $db['all_sth']->get_result();

$alleles = array();
$i       = 0;

print
"  <td style=\"vertical-align: top;\">\n" .
"    <table style=\"margin-left: auto; margin-right: auto; width: 50%;\">\n";

while ($row = $res->fetch_assoc()) {
    $alleles[$row['allele']] = $colors[$i % $color_size];

    print 
      "    <tr>\n" .
      "      <td><span style=\"color: " . $alleles[$row['allele']] . ";\">$row[allele]</span></td>\n" .
      "      <td style=\"text-align: right;\">" . sprintf("%.2f%%", $row['read_pct']) . "</td>\n" .
      "    </tr>\n";

    $i++;
}
print 
"    </table>\n" .
"  </td>\n";

$deleveraged == 1 ? 
    print "  <td style=\"text-align: center; color: #870214;\">True</td>" : 
    print "  <td style=\"text-align: center\">False</td>";
$lumberjackstack == 1 ? 
    print "  <td style=\"text-align: center; color: #870214;\">True</td>" : 
    print "  <td style=\"text-align: center\">False</td>";
$blacklisted == 1 ? 
    print "  <td style=\"text-align: center; color: #870214;\">True</td>" : 
    print "  <td style=\"text-align: center\">False</td>";

echo <<< EOQ
</tr>
</table>

<a name="results_top"></a>
<div class="seq_frame_head">
<table class="radtag">
<tr>
  <th style="width: 5%;">&nbsp;</th>
  <th style="width: 15%;">Relationship</th>
  <th style="width: 20%;">Seq ID</th>
  <th style="width: 60%;">Sequence</th>
</tr>

EOQ;

$con_len = isset($seqs['consensus'][0]) ? strlen($seqs['consensus'][0]['s']) : 0;
$con_seq = isset($seqs['consensus'][0]) ? $seqs['consensus'][0]['s'] : "";

$s = print_scale($con_len);
print
    "<tr>\n" .
    "  <td class=\"num\">&nbsp;</td>\n" .
    "  <td class=\"con\">&nbsp;</td>\n" .
    "  <td class=\"id\">&nbsp;</td>\n" .
    "  <td class=\"tag\">" . $s . "</td>\n" .
    "</tr>\n";

$s = print_snps($display['tag_id'], $con_seq, $con_seq, $snps, false);
print
    "<tr>\n" .
    "  <td class=\"num\">&nbsp;</td>\n" .
    "  <td class=\"con\">consensus</td>\n" .
    "  <td class=\"id\">"  . (isset($seqs['consensus'][0]) ? $seqs['consensus'][0]['id'] : "") . "</td>\n" .
    "  <td class=\"tag\">" . $s . "</td>\n" .
    "</tr>\n";

$s = isset($seqs['model'][0]) ? $seqs['model'][0]['s'] : "";
echo <<< EOQ
</table>
</div>
<div class="seq_frame">
<table class="radtag">
<tr>
  <td class="num">&nbsp;</td>
  <td class="con">model</td>
  <td class="id"></td>
  <td class="tag">$s</td>
</tr>

EOQ;

$i = 1;
foreach (array('primary', 'secondary') as $key) {
    foreach ($seqs[$key] as $seq) {
        $s = print_snps_errs($con_seq, $seq['s'], $snps, $cat_snps);

	if ($key == "primary" && $seq['sub_id'] % 2 == 1) 
	    $bg = "style=\"background-color: #dddddd;\"";
	else
	    $bg = "";

	print
	    "<tr>\n" .
	    "  <td class=\"num\">$i</td>\n" .
	    "  <td class=\"$key\">$key</td>\n" .
	    "  <td class=\"id\">$seq[id]</td>\n" .
	    "  <td class=\"tag\"$bg>$s</td>\n" .
	    "</tr>\n";
	$i++;
    }
}

echo <<< EOQ
</table>
</div>
</div>

EOQ;

write_footer();

?>
