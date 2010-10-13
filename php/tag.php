<?php
require_once("header.php");
require_once("radtag_functions.php");

$tag_id    = isset($_GET['tag_id'])    ? $_GET['tag_id']    : 0;
$batch_id  = isset($_GET['batch_id'])  ? $_GET['batch_id']  : 0;
$sample_id = isset($_GET['sample_id']) ? $_GET['sample_id'] : 0;
$database  = isset($_GET['db'])        ? $_GET['db']        : "radtags";
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
    "SELECT batches.id as id, date, description, samples.id as sample_id, samples.sample_id as samp_id, type, file " . 
    "FROM batches " . 
    "JOIN samples ON (batch_id=batches.id) " . 
    "WHERE batches.id=? AND samples.id=?";
$db['batch_sth'] = $db['dbh']->prepare($query);
check_db_error($db['batch_sth'], __FILE__, __LINE__);

$query = 
    "SELECT count(id) as depth FROM unique_tags " . 
    "WHERE relationship!='consensus' AND sample_id=? AND tag_id=?";
$db['depth_sth'] = $db['dbh']->prepare($query);
check_db_error($db['depth_sth'], __FILE__, __LINE__);

$query = 
    "SELECT tag_id, sub_id, relationship, seq_id, seq, deleveraged, blacklisted, removed " . 
    "FROM unique_tags " . 
    "JOIN samples ON (unique_tags.sample_id=samples.id) " . 
    "JOIN batches ON (samples.batch_id=batches.id) " .
    "WHERE batch_id=? AND unique_tags.sample_id=? AND tag_id=?";
$db['seq_sth'] = $db['dbh']->prepare($query);
check_db_error($db['seq_sth'], __FILE__, __LINE__);

$query = 
    "SELECT col, rank_1, rank_2 FROM snps " . 
    "JOIN samples ON (snps.sample_id=samples.id) " . 
    "JOIN batches ON (samples.batch_id=batches.id) " .
    "WHERE batch_id=? AND snps.sample_id=? AND tag_id=? ORDER BY col";
$db['snp_sth'] = $db['dbh']->prepare($query);
check_db_error($db['snp_sth'], __FILE__, __LINE__);

$query = 
    "SELECT allele, read_pct FROM alleles " . 
    "WHERE sample_id=? AND tag_id=? ";
$db['all_sth'] = $db['dbh']->prepare($query);
check_db_error($db['all_sth'], __FILE__, __LINE__);

$query = 
    "SELECT catalog_id FROM matches " . 
    "WHERE batch_id=? AND sample_id=? AND tag_id=? " .
    "GROUP BY catalog_id";
$db['cat_sth'] = $db['dbh']->prepare($query);
check_db_error($db['par_sth'], __FILE__, __LINE__);

//
// Pull information about this batch
//
$result = $db['batch_sth']->execute(array($batch_id, $sample_id));
check_db_error($result, __FILE__, __LINE__);
$row    = $result->fetchRow();
$batch  = array();
$batch['id']   = $row['id'];
$batch['desc'] = $row['description'];
$batch['date'] = $row['date'];
$batch['file'] = $row['file'];
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
  <img id="sources_img" src="/acos/images/caret-d.png" />
  <a onclick="toggle_div('sources', '/acos/images', 'page_state');">
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
	      'tertiary'  => array());

$result = $db['seq_sth']->execute(array($batch_id, $sample_id, $tag_id));
check_db_error($result, __FILE__, __LINE__);

while ($row = $result->fetchRow()) {
  array_push($seqs[$row['relationship']], array('s' => $row['seq'], 'id' => $row['seq_id'], 'sub_id' => $row['sub_id']));

  if ($row['relationship'] == "consensus") {
    $deleveraged     = $row['deleveraged'];
    $lumberjackstack = $row['removed'];
    $blacklisted     = $row['blacklisted'];
  }
}

$result = $db['cat_sth']->execute(array($batch_id, $sample_id, $tag_id));
check_db_error($result, __FILE__, __LINE__);
$row = $result->fetchRow();

echo <<< EOQ
<td style="text-align: center; vertical-align: top;">
  <a href="$root_path/catalog.php?db=$database&id=$batch_id&filter_type[]=cata&filter_cata=$row[catalog_id]">#$row[catalog_id]</a>
</td>

EOQ;

$result = $db['depth_sth']->execute(array($sample_id, $tag_id));
check_db_error($result, __FILE__, __LINE__);
$row = $result->fetchRow();

echo <<< EOQ
<td style="text-align: center; vertical-align: top;">
  $row[depth]x
</td>

EOQ;

$snps   = array();
$result = $db['snp_sth']->execute(array($batch_id, $sample_id, $tag_id));
check_db_error($result, __FILE__, __LINE__);

while ($row = $result->fetchRow()) {
  $snps[$row['col']] = array('col' => $row['col'], 'rank_1' => $row['rank_1'], 'rank_2' => $row['rank_2']);
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

$result = $db['all_sth']->execute(array($sample_id, $tag_id));
check_db_error($result, __FILE__, __LINE__);

$alleles = array();
$i       = 0;

print
"  <td style=\"vertical-align: top;\">\n" .
"    <table style=\"margin-left: auto; margin-right: auto; width: 50%;\">\n";

while ($row = $result->fetchRow()) {
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
<table class="radtag" style="width: 95%;">
<tr>
  <th style="width: 5%;">&nbsp;</th>
  <th style="width: 15%;">Relationship</th>
  <th style="width: 20%;">Seq ID</th>
  <th style="width: 60%;">Sequence</th>
</tr>

EOQ;

$s = print_scale(strlen($seqs['consensus'][0]['s']));
print
    "<tr>\n" .
    "  <td class=\"num\">&nbsp;</td>\n" .
    "  <td class=\"con\">&nbsp;</td>\n" .
    "  <td class=\"id\">&nbsp;</td>\n" .
    "  <td class=\"tag\">" . $s . "</td>\n" .
    "</tr>\n";

$s = print_snps($seqs['consensus'][0]['s'], $seqs['consensus'][0]['s'], $snps);
print
    "<tr>\n" .
    "  <td class=\"num\">&nbsp;</td>\n" .
    "  <td class=\"con\">consensus</td>\n" .
    "  <td class=\"id\">"  . $seqs['consensus'][0]['id'] . "</td>\n" .
    "  <td class=\"tag\">" . $s . "</td>\n" .
    "</tr>\n";

$i = 1;
foreach (array('primary', 'secondary', 'tertiary') as $key) {
    foreach ($seqs[$key] as $seq) {
	$s = print_snps_errs($seqs['consensus'][0]['s'], $seq['s'], $snps);

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

EOQ;

write_footer();

?>
