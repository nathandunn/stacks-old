<?php
//
// Copyright 2010-2015, Julian Catchen <jcatchen@illinois.edu>
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

$batch_id  = isset($_GET['id']) ? $_GET['id'] : 0;
$database  = isset($_GET['db']) ? $_GET['db'] : "";
$page      = isset($_GET['p'])  ? $_GET['p']  : 1;
$per_page  = isset($_GET['pp']) ? $_GET['pp'] : 10;

// Connect to the database
$db = db_connect($database);

// Save these variables for automatic URL formation later on.
$display = array();
$display['id']          = $batch_id;
$display['db']          = $database;
$display['p']           = $page;
$display['pp']          = $per_page;
$display['filter_type'] = array();

//
// Process the filtering parameters
//
$param = array($batch_id);
process_filter($display);
prepare_filter_parameters($display, $param);

//
// Prepare some SQL queries
//
$query = 
    "SELECT batches.id as id, date, description, type FROM batches " . 
    "WHERE batches.id=?";
$db['batch_sth'] = $db['dbh']->prepare($query);
check_db_error($db['batch_sth'], __FILE__, __LINE__);

$query = 
    "SELECT COUNT(tag_id) as count FROM catalog_index " . 
    "WHERE batch_id=?";
$query .= apply_query_filters($display);
$db['count_sth'] = $db['dbh']->prepare($query);
check_db_error($db['count_sth'], __FILE__, __LINE__);

$query = 
    "SELECT alleles as count FROM catalog_index " . 
    "WHERE batch_id=? AND tag_id=?";
$db['allele_sth'] = $db['dbh']->prepare($query);
check_db_error($db['allele_sth'], __FILE__, __LINE__);

$query = 
    "SELECT col, rank_2 FROM catalog_snps " . 
    "JOIN batches ON (catalog_snps.batch_id=batches.id) " . 
    "WHERE batch_id=? AND tag_id=? ORDER BY col";
$db['snp_sth'] = $db['dbh']->prepare($query);
check_db_error($db['snp_sth'], __FILE__, __LINE__);

$query = 
  "SELECT chr, max_len FROM chr_index " . 
  "WHERE batch_id=?";
$db['chrs_sth'] = $db['dbh']->prepare($query);
check_db_error($db['chrs_sth'], __FILE__, __LINE__);

$query =
    "SELECT max(ests) as ests, max(pe_radtags) as pe_radtags, max(blast_hits) as blast_hits " . 
    "FROM catalog_index WHERE batch_id=?";
$db['seq_sth'] = $db['dbh']->prepare($query);
check_db_error($db['seq_sth'], __FILE__, __LINE__);

$result = $db['seq_sth']->execute($batch_id);
check_db_error($result, __FILE__, __LINE__);
$row = $result->fetchRow();


$cols = array();

if ($row['ests'] == 0 && 
    $row['pe_radtags'] == 0 && 
    $row['blast_hits'] == 0)
    $cols['seq'] = false;
else
    $cols['seq'] = true;

$query =
    "SELECT count(id) as cnt FROM catalog_genotypes WHERE batch_id=?";
$db['gcnt_sth'] = $db['dbh']->prepare($query);
check_db_error($db['gcnt_sth'], __FILE__, __LINE__);

$result = $db['gcnt_sth']->execute($batch_id);
check_db_error($result, __FILE__, __LINE__);
$row = $result->fetchRow();

if ($row['cnt'] > 0)
  $cols['gcnt'] = true;
else
  $cols['gcnt'] = false;

//
// Pull information about this batch
//
$result = $db['batch_sth']->execute($batch_id);
check_db_error($result, __FILE__, __LINE__);
$row    = $result->fetchRow();
$batch  = array();
$batch['id']   = $row['id'];
$batch['desc'] = $row['description'];
$batch['date'] = $row['date'];
$batch['type'] = $row['type'];

$page_title = "RAD-Tag Catalog Viewer";
write_header($page_title, $batch);

echo <<< EOQ
<script type="text/javascript"
        src="https://ajax.googleapis.com/ajax/libs/jquery/2.1.3/jquery.min.js"></script>
<script type="text/javascript"
        src="$root_path/population_view.js"></script>

<h3 class="info_head" style="margin-left: 1em;">
  <a href="$root_path/index.php?db=$database">
  Batch #$batch[id] <span class="s">[$batch[date]; $batch[desc]]</span></a>
</h3>

EOQ;

if ($batch['type'] == "population")
  write_pop_filter($cols);
else
  write_map_filter($cols);

//
// How many columns will we print
//
$num_cols = 9;
foreach ($cols as $col)
    if ($col == false)
        $num_cols--;

//
// Generate Excel export URL
//
$excel_export = generate_url("export_batch.php", false);

echo <<< EOQ
<div id="sources" style="width: 100%;">
<a name="results_top"></a>
<table class="db" style="width: 100%; border: none;">
<tr>
  <td colspan="$num_cols" style="border: none; padding-bottom: 0px;" class="export_icon">

<div id="export_popup" style="display: none;">
<h3>Export</h3>
<div id="export_popup_txt">
<p>This will take some time. Enter your email and you will be notified
when the results are ready.
</p>
<p>
  <form id="export_popup_frm">
  <input type="hidden" name="url" value="$excel_export" />
  <strong>E-mail:</strong> <input type="input" size=25 name="email" value="" />
</p>
<p>
  <strong>Data:</strong>
    <input type="radio" name="dtype" value="haplo" checked="checked" onchange="toggle_vis('export_popup_frm', 'dtype')" />Observed Haplotypes
    <input type="radio" name="dtype" value="geno" onchange="toggle_vis('export_popup_frm', 'dtype')" />Genotypes<br />
</p>
<p id="hopts" style="display: none;">
<strong>Min Stack Depth:</strong>
<select name="dlim">
  <option checked="checked">1</option>
  <option>3</option>
  <option>5</option>
  <option>10</option>
  <option>15</option>
</select>
</p>
<p id="gopts" style="display: none;">
  <strong><acronym title="Gen: Generic, F2: F2 cross, CP: F1 cross, BC1: Back cross, DH: Double Haploid cross">Map Type</acronym>:</strong>
<select name="mtype">
  <option value="gen">Gen</option>
  <option value="f2">F2</option>
  <option value="cp">CP</option>
  <option value="bc1">BC1</option>
  <option value="dh">DH</option>
</select>
<input type="checkbox" name="mcor" value="1" />
<strong><a onclick="toggle_cb('export_popup_frm', 'mcor')">Include manual corrections</a></strong>
</p>

<p>
  <strong>Output type:</strong>
    <input type="radio" name="otype" value="tsv" /><acronym title="Tab-separated Values Format">TSV</acronym>
    <input type="radio" name="otype" value="xls" /><acronym title="Microsoft Excel Format">XLS</acronym><br />
</p>
<p>
  <a onclick="export_data('export_popup')">submit</a> | <a onclick="close_export_popup('export_popup')">cancel</a>
</p>
  </form>
</div>
</div>

    <a onclick="toggle_export_popup('export_popup'); toggle_vis('export_popup_frm', 'dtype')">
    <img style="float: right;" src="$root_path/images/excel_icon.png" 
         title="Export data to Microsoft Excel format"/></a>
  </td>
</tr>
<tr>
  <td colspan="$num_cols" style="border: none; padding-bottom: 0px;">

EOQ;

//
// Figure out how many results there are (including filtering)
// and write out the proper pagination links
//
$result = $db['count_sth']->execute($param);
check_db_error($result, __FILE__, __LINE__);
$row = $result->fetchRow();
$pagination_count = $row['count'];
$start_group = 0;
$end_group   = 0;

write_pagination($pagination_count, $start_group, $end_group, "catalog.php");

if ($batch['type'] == "map") {
  echo <<< EOQ
  </td>
</tr>
<tr>
  <th style="width: 10%;">Id</th>
  <th style="width: 5%;">SNP</th>
  <th style="width: 45%;">Consensus</th>
  <th style="width: 12%;">Matches</th>
  <th style="width: 5%;">Marker</th>
  <th style="width: 10%;">Ratio</th>
  <th style="width: 3%; font-size: smaller;">ChiSq P-value</th>

EOQ;
} else {
  echo <<< EOQ
  </td>
</tr>
<tr>
  <th style="width: 10%;">Id</th>
  <th style="width: 5%;">SNP</th>
  <th style="width: 50%;">Consensus</th>
  <th style="width: 5%;">Matches</th>
  <th style="width: 15%;">Ratio</th>

EOQ;
}

if ($cols['seq'] == true)
  print "  <th style=\"width: 10%\">Sequence</th>\n";
print  "</tr>\n";

$db['dbh']->setLimit($display['pp'], $start_group - 1);
check_db_error($db['dbh'], __FILE__, __LINE__);

$query = 
    "SELECT catalog_index.tag_id as tag_id, alleles, parents, progeny, valid_progeny, " . 
    "seq, marker, uncor_marker, chisq_pval, lnl, ratio, ests, pe_radtags, blast_hits, external_id, geno_cnt, " .
    "catalog_index.chr, catalog_index.bp, catalog_tags.strand, catalog_index.type, gene, ext_id, ex_start, ex_end, ex_index " .
    "FROM catalog_index " .
    "JOIN catalog_tags ON (catalog_index.cat_id=catalog_tags.id) " . 
    "LEFT JOIN catalog_annotations ON (catalog_index.batch_id=catalog_annotations.batch_id AND catalog_index.tag_id=catalog_annotations.catalog_id) " .
    "LEFT JOIN ref_radome ON (catalog_index.ref_id=ref_radome.id) " .
    "WHERE catalog_index.batch_id=?";
$query .= apply_query_filters($display);

$db['tag_sth'] = $db['dbh']->prepare($query);
check_db_error($db['tag_sth'], __FILE__, __LINE__);

$result = $db['tag_sth']->execute($param);
check_db_error($result, __FILE__, __LINE__);

while ($row = $result->fetchRow()) {

    // Query the database to find how many SNPs were found in this sample.
    $snps = array();
    $snp_res = $db['snp_sth']->execute(array($batch_id, $row['tag_id']));
    check_db_error($snp_res, __FILE__, __LINE__);
    while ($snp_row = $snp_res->fetchRow()) {
	array_push($snps, array('col' => $snp_row['col'], 'rank' => $snp_row['rank_2']));
    }

    $url = "$root_path/pop_view.php?db=$database&batch_id=$batch_id&type=$batch[type]&tag_id=$row[tag_id]";
    $annotation = strlen($row['external_id']) > 0 ? $row['external_id'] : "annotate";

    echo <<< EOQ
<tr class="catrow">
  <td class="catlink">
<div onclick="ajax_locus_population_view('$row[tag_id]', '$img_path', '$url');">
<img id="{$row['tag_id']}_img" src="$img_path/caret-u.png" />
<a>$row[tag_id]</a>
</div>

<div class="ann_parent">
  <a class="annotation" id="{$row['tag_id']}_ann" onclick="toggle_annotation($row[tag_id])">$annotation</a>
  <div id="{$row['tag_id']}_div" style="display: none; font-size: small;">
    <form id="{$row['tag_id']}_frm">
      <input type="hidden" name="url" value="$root_path/annotate_marker.php?db=$database&batch_id=$batch_id&tag_id=$row[tag_id]" />
      <input type="input" size=15 name="ext_id" value="" />
      <a onclick="annotate_marker('$row[tag_id]')">save</a>|<a onclick="toggle_annotation('$row[tag_id]')">cancel</a>
    </form>
  </div>
</div>
  </td>

EOQ;

    if (count($snps) == 0)
	print "  <td>No</td>\n";
    else
	print "  <td>Yes <span class=\"s\">[" . count($snps) . "nuc]</span></td>\n";

    $s = print_snps($row['tag_id'], $row['seq'], $row['seq'], $snps, true);

    $ratio = explode(";", $row['ratio']);
    $ratio_parsed = "";
    $i            = 0;
    foreach ($ratio as $r) {
      if (strlen($r) == 0) continue;

      preg_match("/([a-z]+):(\d+)\((\d+\.?\d*%)\)/", $r, $matches);

      for ($j = 0; $j < strlen($matches[1]); $j++) {
	  $color = $color_map[$matches[1][$j]];
	  $ratio_parsed .= "<span style=\"color: $color; font-weight: bold\">" . $matches[1][$j] . "</span>";
      }

      $ratio_parsed .= ": $matches[2]";

      if ($matches[3] > 0)
	$ratio_parsed .= " ($matches[3])";

      $ratio_parsed .= "<br />";
      $i++;
    }

    $url = "$root_path/sequence_blast.php?db=$database&batch_id=$batch_id&tag_id=$row[tag_id]";
    if ($row['blast_hits'] > 0 || $row['pe_radtags'] > 0 || $row['ests'] > 0) {
        $blast_hits_str =
            "<div class=\"catlink\"><img id=\"{$row['tag_id']}_blast_img\" src=\"$img_path/caret-u.png\" />" .
            "<a onclick=\"toggle_aln_tr('{$row['tag_id']}_blast', '$img_path', '$url');\">" .
            "blast hits: $row[blast_hits]</a>";
    } else {
        $blast_hits_str = "blast hits: $row[blast_hits]";
    }

    if (strlen($row['chr']) > 0) {
      print 
	"<td class=\"seq\">\n" .
	"<div class=\"seq\">$s</div>\n" .
	"<div class=\"gloc\">Chr: $row[chr], <acronym title=\"" . number_format($row['bp']) . "bp\">" . print_bp($row['bp']) . "</acronym>, $row[strand]\n";

      if ($row['type'] == "exon") {
  	  if (strlen($row['ext_id']) == 0)
	    $gene = $row['gene'];
	  else
	    $gene = $row['ext_id'];
	  print 
	    ", Gene: <acronym title=\"Exon $row[ex_index]: " . number_format($row['ex_start']) . "-" . number_format($row['ex_end']) . "bp\">" . $gene . "</acronym>\n";
      }
      print
	", LnL: $row[lnl]</div></td>\n";

    } else {
      print 
	"<td class=\"seq\"><div class=\"seq\">$s</div>\n" .
	"<div class=\"gloc\">LnL: $row[lnl]</div>\n" .
	"</td>\n";
    }

    if ($batch['type'] == "map") {
      print 
	"  <td style=\"font-size: smaller;\"><acronym title=\"Matching Parents\"><span style=\"color: $colors[4]\">$row[parents]</span></acronym>" .
	"<strong> / </strong> " .  
	"<acronym title=\"Matching Progeny\"><span style=\"color: $colors[5]\">$row[progeny]</span></acronym>" . 
	"<strong> / </strong> " . 
	"<acronym title=\"Mappable Progeny\"><span style=\"color: $colors[6]\">$row[valid_progeny]</span></acronym>" . 
	"<strong> / </strong> " . 
	"<acronym title=\"Assigned Genotypes\"><span style=\"color: $colors[7]\">$row[geno_cnt]</span></acronym></td>\n";

      if (strlen($row['uncor_marker']) > 0 && $row['marker'] != $row['uncor_marker'])
	print "  <td><acronym title=\"Uncorrected marker: $row[uncor_marker]\">$row[marker]*</acronym></td>\n";
      else
	print "  <td>$row[marker]</td>\n";
    } else {
      print
	"  <td>$row[parents]</td>\n";
    }

    echo <<< EOQ
  <td style="text-align: left; font-size: smaller;">
    $ratio_parsed
  </td>

EOQ;

    if ($batch['type'] == "map")
      print "<td>$row[chisq_pval]</td>\n";

    if ($cols['seq'] == true)
        echo <<< EOQ
  <td>
    <table class="int" style="font-size: smaller;">
    <tr><td>ests: $row[ests]</td><td>pe: $row[pe_radtags]</td></tr>
    <tr>
      <td colspan="2">
      $blast_hits_str
      </td></tr>
    </table>
  </td>

EOQ;

echo <<< EOQ
</tr>
<tr id="{$row['tag_id']}" style="display: none">
  <td colspan="$num_cols">
    <div id="{$row['tag_id']}_popview_div" class="popview"></div>
  </td>
</tr>
<tr id="{$row['tag_id']}_blast" style="display: none">
  <td colspan="$num_cols">
    <iframe id="{$row['tag_id']}_blast_iframe" 
            frameborder="0" 
            scrolling="auto" 
            onload="this.style.height = (this.contentWindow.document.body.offsetHeight+25) + 'px';" 
            src=""></iframe>
  </td>
</tr>
EOQ;
}

print 
"<tr>\n" .
"  <td colspan=\"$num_cols\" style=\"border: none; padding-top: 0px;\">\n";

write_pagination($pagination_count, $start_group, $end_group, "catalog.php");

echo <<< EOQ
  </td>
</tr>
</table>
</div>

EOQ;

write_footer();

function generate_hidden_form_vars($var) {
    global $root_path, $display;

    $vars = "";

    foreach ($display as $key => $d) {
	if (strstr($key, $var))
	    continue;

	if (is_array($d)) {
	    foreach ($d as $e) {
		$vars .= "  <input type=\"hidden\" name=\"{$key}[]\" value=\"$e\" />\n";
	    }
	} else {
	    $vars .= "  <input type=\"hidden\" name=\"$key\" value=\"$d\" />\n";
	}
    }

    return $vars;
}

function generate_per_page_select($name, $per_page) {
    
    $pages = array("10", "50", "100", "all");

    $ctl = "  <select name=\"$name\" onchange=\"this.form.submit();\">\n";

    foreach ($pages as $p) {
	if ($p == $per_page) 
	    $ctl .= "  <option selected=\"selected\">$p</option>\n";
	else
	    $ctl .= "  <option>$p</option>\n";
    }

    $ctl .= "  </select>\n";

    return $ctl;
}

function generate_url($destination, $prefix) {
    global $root_path, $display;

    if ($prefix) 
        $url = "href=\"" . $root_path . "/" . $destination . "?";
    else 
        $url = $root_path . "/" . $destination . "?";

    foreach ($display as $key => $d) {
	if (is_array($d)) {
	    foreach ($d as $e) 
		$url .= "{$key}[]=$e&";
	} else {
	    $url .= "$key=$d&";
	}
    }

    // Remove the hanging '&'
    $url = substr($url, 0, -1);

    if ($prefix)
        $url .= "\"";

    return $url;
}

function generate_page_list($page, $num_pages, $destination) {
    global $display;

    $page_list = "";

    if ($page <= 4) {
	for ($i = 1; $i < $page; $i++) {
	    $display['p'] = $i;
	    $p            = generate_url($destination, true);
	    $page_list   .= "<a $p>$i</a>\n"; 
	} 
    } else {
	$display['p'] = 1;
	$p            = generate_url($destination, true);
	$page_list   .= "<a $p>1</a> ...\n"; 

	foreach (array($page - 3, $page - 2, $page - 1) as $i) {
	    $display['p'] = $i;
	    $p            = generate_url($destination, true);
	    $page_list   .= "<a $p>$i</a>\n"; 
	} 
    }

    $page_list .= " <strong>$page</strong>\n";

    if ($page <= $num_pages - 4) {
	for ($i = $page + 1; $i <= $page + 3; $i++) {
	    $display['p'] = $i;
	    $p            = generate_url($destination, true);
	    $page_list   .= "<a $p>$i</a>\n"; 
	} 

	$display['p'] = $num_pages;
	$p            = generate_url($destination, true);
	$page_list   .= "... <a $p>$num_pages</a>\n"; 

    } else {
	for ($i = $page + 1; $i <= $num_pages; $i++) {
	    $display['p'] = $i;
	    $p            = generate_url($destination, true);
	    $page_list   .= "<a $p>$i</a>\n"; 
	} 
    }

    $display['p'] = $page;

    $page_list =
	"<td class=\"s\">\n" .
	$page_list .
	"</td>\n";

    return $page_list;
}

function write_pagination($num_tags, &$start_gene, &$end_gene, $destination) {
    global $img_path, $root_path, $display;

    $cur_page = $display['p'];
    $page     = $display['p'];
    $per_page = $display['pp'];

    if ($per_page == "all") 
	$per_page = $num_tags;

    //
    // First figure out the total number of pages. If there are
    // additional genes left over, add an extra page.
    //
    $num_pages = floor($num_tags / $per_page);
    $num_pages += $num_tags % $per_page >= 1 ? 1 : 0;

    if ($page > $num_pages) {
	$page = $num_pages;
	$cur_page = $num_pages;
    }

    // Determine the start and end gene numbers
    $start_gene = 1 + (($page - 1) * $per_page);
    $end_gene   = 
	$start_gene + $per_page > $num_tags ? 
	$num_tags : 
	($start_gene + $per_page - 1);

    // Generate the URLs for our links
    $display['p'] -= 1;
    $prev_page     = generate_url($destination, true);
    $display['p'] += 2;
    $next_page     = generate_url($destination, true);
    $display['p']  = $cur_page;

    print 
	"<table class=\"int\">\n" .
	"<tr>\n" .
	"<td style=\"text-align: left; width: 20%;\">\n";

    if ($page == 1) {
	if ($num_pages == 1) {
	    echo <<< EOQ
<img style="vertical-align: bottom;" src="$img_path/l-arrow-disabled.png" alt="No Previous Page" />
 <strong>$page</strong>
<img style="vertical-align: bottom;" src="$img_path/r-arrow-disabled.png" alt="No Next Page" />

EOQ;
	} else {

	    echo <<< EOQ
<img style="vertical-align: bottom;" src="$img_path/l-arrow-disabled.png" alt="No Previous Page" />
 <strong>$page</strong>
<a $next_page title="View Next Page">
<img style="vertical-align: bottom;" src="$img_path/r-arrow.png" alt="View Next Page" /></a>

EOQ;
	}
    } else if ($page == $num_pages) {
	echo <<< EOQ
<a $prev_page title="View Previous Page">
<img style="vertical-align: bottom;" src="$img_path/l-arrow.png" alt="View Previous Page" /></a>
 <strong>$page</strong>
<img style="vertical-align: bottom;" src="$img_path/r-arrow-disabled.png" alt="No Next Page" />

EOQ;
    } else {
	echo <<< EOQ
<a $prev_page title="View Previous Page">
<img style="vertical-align: bottom;" src="$img_path/l-arrow.png" alt="View Previous Page" /></a>
 <strong>$page</strong>
<a $next_page title="View Next Page">
<img style="vertical-align: bottom;" src="$img_path/r-arrow.png" alt="View Next Page" /></a>

EOQ;
    }

    print
	" <span class=\"s\">($num_tags tags)</span>\n" .
	"</td>\n";

    $page_list = "";

    if ($num_pages > 1) 
	$page_list = generate_page_list($page, $num_pages, $destination);

    print $page_list;

    $hidden_vars  = generate_hidden_form_vars("pp");
    $per_page_ctl = generate_per_page_select("pp", $display['pp']);

    echo <<< EOQ
  <td class="s" style="width: 20%; text-align: right;">
  <form id="per_page" name="per_page" method="get" action="$root_path/$destination">
$hidden_vars
tags per page &nbsp;
$per_page_ctl
  </form>
  </td>
</tr>
</table>

EOQ;
}

function write_map_filter($cols) {
    global $img_path, $root_path, $display;

    $max_chr_len = 0;
    $hidden_vars = generate_hidden_form_vars("filter");
    $chrs        = fetch_chrs($max_chr_len);

    $filters = array("cata"  => array(),
		     "alle"  => array(),
		     "snps"  => array(),
		     "pare"  => array(),
		     "prog"  => array(),
		     "vprog" => array(),
		     "mark"  => array(),
		     "gcnt"  => array(),
		     "chisq" => array(),
		     "loc"   => array(),
		     "ref"   => array(),
                     "est"   => array(),
                     "pe"    => array(),
                     "blast" => array());

    $fall = isset($display['filter_alle_l'])  ? $display['filter_alle_l']  : "";
    $falu = isset($display['filter_alle_u'])  ? $display['filter_alle_u']  : "100";
    $fsnl = isset($display['filter_snps_l'])  ? $display['filter_snps_l']  : "";
    $fsnu = isset($display['filter_snps_u'])  ? $display['filter_snps_u']  : "100";
    $fpal = isset($display['filter_pare_l'])  ? $display['filter_pare_l']  : "";
    $fpau = isset($display['filter_pare_u'])  ? $display['filter_pare_u']  : "2";
    $fpr  = isset($display['filter_prog'])    ? $display['filter_prog']    : "";
    $fvp  = isset($display['filter_vprog'])   ? $display['filter_vprog']   : "";
    $fgc  = isset($display['filter_gcnt'])    ? $display['filter_gcnt']    : "";
    $csql = isset($display['filter_chisq_l']) ? $display['filter_chisq_l'] : "0.0";
    $csqu = isset($display['filter_chisq_u']) ? $display['filter_chisq_u'] : "1.0";
    $fma  = isset($display['filter_mark'])    ? $display['filter_mark']    : "";
    $ref  = isset($display['filter_ref'])     ? $display['filter_ref']     : "";
    $fch  = isset($display['filter_chr'])     ? $display['filter_chr']     : "";
    $fsb  = isset($display['filter_sbp'])     ? $display['filter_sbp']     : 0;
    $feb  = isset($display['filter_ebp'])     ? $display['filter_ebp']     : $max_chr_len;

    $r = range(1, 9);
    $r = array_merge($r, range(10, 100, 5));
    array_push($r, 1000);
    $alle_l_ctl  = generate_element_select("filter_alle_l",  $r,   $fall, "");
    $alle_u_ctl  = generate_element_select("filter_alle_u",  $r,   $falu, "");
    $snps_l_ctl  = generate_element_select("filter_snps_l",  $r,   $fsnl, "");
    $snps_u_ctl  = generate_element_select("filter_snps_u",  $r,   $fsnu, "");
    $r = range(1, 9);
    $r = array_merge($r, range(10, 500, 10));
    array_push($r, 1000, 2000, 10000);
    $pare_l_ctl  = generate_element_select("filter_pare_l",  $r, $fpal, "");
    $pare_u_ctl  = generate_element_select("filter_pare_u",  $r, $fpau, "");
    $prog_ctl    = generate_element_select("filter_prog",  $r, $fpr, "");
    $vprog_ctl   = generate_element_select("filter_vprog", $r, $fvp, "");
    $gcnt_ctl    = generate_element_select("filter_gcnt",  $r, $fgc, "");
    $csql_ctl    = generate_element_select("filter_chisq_l", array(0.0, 0.0005, 0.001, 0.005, 0.01, 0.05, 0.1, 0.15, 0.2, 0.25, 0.5, 1.0), $csql, "");
    $csqu_ctl    = generate_element_select("filter_chisq_u", array(0.0, 0.0005, 0.001, 0.005, 0.01, 0.05, 0.1, 0.15, 0.2, 0.25, 0.5, 1.0), $csqu, "");
    $chr_ctl     = generate_element_select("filter_chr",   $chrs, $fch, "");
    $sbp_ctl     = generate_element_select("filter_sbp",   range(0, $max_chr_len), $fsb, "");
    $ebp_ctl     = generate_element_select("filter_ebp",   range(0, $max_chr_len), $feb, "");
    $ref_ctl     = generate_element_select("filter_ref",  array("exon", "intron", "genomic"), $ref, "");
    $mark_ctl    = generate_element_select("filter_mark", 
					   array('Any', 'aa/bb', 'ab/--', '--/ab', 'aa/ab', 'ab/aa', 'ab/ab', 'ab/ac', 'ab/cd', 'ab/cc', 'cc/ab'), 
					   $fma, "");

    if (isset($display['filter_type'])) {

	foreach ($filters as $key => $f)
	    if (in_array($key, $display['filter_type'])) {
		$filters[$key]['sel'] = "checked=\"checked\"";
		$filters[$key]['tr']  = "class=\"active_filter\"";
	    } else {
		$filters[$key]['sel'] = "";
		$filters[$key]['tr']  = "";
	    }

    } else {
	$filters['none']['sel'] = "checked=\"checked\"";
    }

    echo <<< EOQ
<h4 class="info_head">
  <img id="stacks_filter_img" src="$img_path/caret-d.png" />
  <a onclick="toggle_div('stacks_filter', '$img_path', 'page_state');">Filter Results By</a>
</h4>
<div id="stacks_filter">
<form id="filter_results" name="filter_results" method="get" action="$root_path/catalog.php">
$hidden_vars

<table style="width: 100%; vertical-align: top;">
<tr>
EOQ;

    if (count($chrs) > 0) {
      echo <<< EOQ
<td style="width: 25%; vertical-align: top;">
<table class="loc_filter">
<tr>
  <td {$filters['loc']['tr']}>
      <input type="checkbox" name="filter_type[]" value="loc" onchange="rebuild_display_select()" {$filters['loc']['sel']} /> 
      <a onclick="toggle_cb('filter_results', 'loc')">
      <acronym title="Filter by genomic location">Location</acronym>:</a></td>
  <td {$filters['loc']['tr']}>
      $chr_ctl
  </td>
</tr>
<tr>
  <td {$filters['loc']['tr']}><span style="padding-left: 1.5em;">Start:</span></td>
  <td {$filters['loc']['tr']}>
      $sbp_ctl Mb
  </td>
</tr>
<tr>
  <td {$filters['loc']['tr']}><span style="padding-left: 1.5em;">End:</span></td>
  <td {$filters['loc']['tr']}>
      $ebp_ctl Mb
  </td>
</tr>
<tr>
  <td {$filters['ref']['tr']}>
      <input type="checkbox" name="filter_type[]" value="ref" onchange="rebuild_display_select()" {$filters['ref']['sel']} /> 
      <a onclick="toggle_cb('filter_results', 'ref')">
      <acronym title="Filter by type of RAD locus">Type</acronym>:</a></td>
  <td {$filters['ref']['tr']}>
      $ref_ctl
  </td>
</tr>
</table>
</td>

EOQ;

    }

    $cat_id_filter = isset($display['filter_cata']) ? $display['filter_cata'] : "";

    echo <<< EOQ
<td>
<table class="filter">
<tr>
  <td colspan="2" {$filters['cata']['tr']}>
      <input type="checkbox" name="filter_type[]" value="cata" onchange="rebuild_display_select()" {$filters['cata']['sel']} /> 
      <a onclick="toggle_cb('filter_results', 'cata')">
      <acronym title="Show a locus with a particular ID.">Catalog ID</acronym>:</a>
      <input name="filter_cata" value="$cat_id_filter" size="15" />
  </td>
  <td style="text-align: right; padding-right: 10px;">
      <input type="submit" value="filter" />
  </td>
</tr>
<tr>
  <td style="width: 27%; text-align: left;">
  <table style="text-align: left;">
  <tr>
  <td {$filters['alle']['tr']}>
      <input type="checkbox" name="filter_type[]" value="alle" onchange="rebuild_display_select()" {$filters['alle']['sel']} /> 
      <a onclick="toggle_cb('filter_results', 'alle')">
      <acronym title="Filter the catalog according to the number of alleles identified for a locus.">Alleles</acronym>:</a>
  </td>
  <td {$filters['alle']['tr']}>
$alle_l_ctl $alle_u_ctl
  </td>
  </tr>
  <tr>
  <td {$filters['snps']['tr']}>
      <input type="checkbox" name="filter_type[]" value="snps" onchange="rebuild_display_select()" {$filters['snps']['sel']} /> 
      <a onclick="toggle_cb('filter_results', 'snps')">
      <acronym title="Filter the catalog according to the number of SNPs found at a locus.">SNPs</acronym>:</a>
  </td>
  <td {$filters['snps']['tr']}>
$snps_l_ctl $snps_u_ctl
  </td>
  </tr>
  <tr>
    <td>&nbsp;</td>
    <td>&nbsp;</td>
  </tr>
  </table>
  </td>

  <td style="width: 39%; text-align: center;">
  <table style="text-align: left;">
  <tr>
  <td {$filters['pare']['tr']}>
      <input type="checkbox" name="filter_type[]" value="pare" onchange="rebuild_display_select()" {$filters['pare']['sel']} /> 
      <a onclick="toggle_cb('filter_results', 'pare')">
      <acronym title="Filter the catalog according to the number of parental samples that are matched to a locus.">Parental matches</acronym>:</a>
  </td>
  <td {$filters['pare']['tr']}>
$pare_l_ctl $pare_u_ctl
  </td>
  </tr>
  <tr>
  <td {$filters['prog']['tr']}>
      <input type="checkbox" name="filter_type[]" value="prog" onchange="rebuild_display_select()" {$filters['prog']['sel']} /> 
      <a onclick="toggle_cb('filter_results', 'prog')">
      <acronym title="Filter the catalog according to the number of progeny samples that are matched to a locus.">Progeny matches</acronym>:</a>
  </td>
  <td {$filters['prog']['tr']}>
$prog_ctl
  </td>
  </tr>
  <tr>
  <td {$filters['chisq']['tr']}>
      <input type="checkbox" name="filter_type[]" value="chisq" onchange="rebuild_display_select()" {$filters['chisq']['sel']} /> 
      <a onclick="toggle_cb('filter_results', 'chisq')">
      <acronym title="A chi-square test for segregation distortion. The smaller the p-value, the more likely the genotype ratios are distorted from the expectation for this marker.">Segregation distortion</acronym>:</a></td>
  <td {$filters['chisq']['tr']}>
$csql_ctl $csqu_ctl
  </td>
  </tr>
  </table>
  </td>

  <td style="width: 33%; text-align: right;">
  <table style="text-align: left;">
  <tr>
  <td {$filters['vprog']['tr']}>
      <input type="checkbox" name="filter_type[]" value="vprog" onchange="rebuild_display_select()" {$filters['vprog']['sel']} /> 
      <a onclick="toggle_cb('filter_results', 'vprog')">
      <acronym title="Filter the catalog according to the number of progeny for which a genotype could be inferred.">Mappable progeny</acronym>:</a></td>
  <td {$filters['vprog']['tr']}>
$vprog_ctl
  </td>
  </tr>
  <tr>
  <td {$filters['mark']['tr']}>
      <input type="checkbox" name="filter_type[]" value="mark" onchange="rebuild_display_select()" {$filters['mark']['sel']} /> 
      <a onclick="toggle_cb('filter_results', 'mark')">
      <acronym title="Filter the catalog to display loci for which a mappable marker type could be inferred.">Mappable markers</acronym>:</a>
  </td>
  <td {$filters['mark']['tr']}>
$mark_ctl
  </td>
  </tr>
  <tr>
  <td {$filters['gcnt']['tr']}>
      <input type="checkbox" name="filter_type[]" value="gcnt" onchange="rebuild_display_select()" {$filters['gcnt']['sel']} /> 
      <a onclick="toggle_cb('filter_results', 'gcnt')">
      <acronym title="Filter the catalog to show loci for which genotypes have been called.">Genotypes</acronym>:</a></td>
  <td {$filters['gcnt']['tr']}>
$gcnt_ctl
  </td>
  </tr>
  </table>
  </td>
</tr>

EOQ;

  if ($cols['seq'] == true)
    echo <<< EOQ
<tr>
  <td {$filters['est']['tr']}>
      <input type="checkbox" name="filter_type[]" value="est" onchange="rebuild_display_select()" {$filters['est']['sel']} /> 
      <a onclick="toggle_cb('filter_results', 'est')">
      <acronym title="Filter the catalog to show loci for which ESTs have been associated.">Contains ESTs</acronym></a>
  </td>
  <td {$filters['pe']['tr']}>
      <input type="checkbox" name="filter_type[]" value="pe" onchange="rebuild_display_select()" {$filters['pe']['sel']} /> 
      <a onclick="toggle_cb('filter_results', 'pe')">
      <acronym title="Filter the catalog to show loci for which paired-end RAD-Tags have been associated.">Contains Paired-end RAD-Tags</acronym></a>
  </td>
  <td {$filters['blast']['tr']}>
      <input type="checkbox" name="filter_type[]" value="blast" onchange="rebuild_display_select()" {$filters['blast']['sel']} /> 
      <a onclick="toggle_cb('filter_results', 'blast')">
      <acronym title="Filter the catalog to show loci for which BLAST hits  have been associated.">Contains BLAST Hits</acronym></a>
  </td>
</tr>

EOQ;

    echo <<< EOQ
</table>
</td>
</tr>
</table>
</form>
</div>

EOQ;

}

function write_pop_filter($cols) {
    global $img_path, $root_path, $display;

    $max_chr_len = 0;
    $hidden_vars = generate_hidden_form_vars("filter");
    $chrs        = fetch_chrs($max_chr_len);

    $filters = array("cata"  => array(),
		     "alle"  => array(),
		     "snps"  => array(),
		     "pare"  => array(),
		     "gcnt"  => array(),
		     "lnl"   => array(),
		     "loc"   => array(),
		     "ref"   => array(),
                     "est"   => array(),
                     "pe"    => array(),
                     "blast" => array());

    $fall = isset($display['filter_alle_l'])  ? $display['filter_alle_l'] : "";
    $falu = isset($display['filter_alle_u'])  ? $display['filter_alle_u'] : "100";
    $fsnl = isset($display['filter_snps_l'])  ? $display['filter_snps_l'] : "";
    $fsnu = isset($display['filter_snps_u'])  ? $display['filter_snps_u'] : "100";
    $fpal = isset($display['filter_pare_l'])  ? $display['filter_pare_l'] : "";
    $fpau = isset($display['filter_pare_u'])  ? $display['filter_pare_u'] : "1000";
    $ref  = isset($display['filter_ref'])     ? $display['filter_ref']    : "";
    $fch  = isset($display['filter_chr'])     ? $display['filter_chr']    : "";
    $fsb  = isset($display['filter_sbp'])     ? $display['filter_sbp']    : 0;
    $feb  = isset($display['filter_ebp'])     ? $display['filter_ebp']    : $max_chr_len;
    $flnl = isset($display['filter_lnl_l'])   ? $display['filter_lnl_l']  : 0;
    $flnu = isset($display['filter_lnl_u'])   ? $display['filter_lnl_u']  : -500;

    $r = range(1, 9);
    $r = array_merge($r, range(10, 100, 5));
    array_push($r, 1000);
    $alle_l_ctl  = generate_element_select("filter_alle_l",  $r,  $fall, "");
    $alle_u_ctl  = generate_element_select("filter_alle_u",  $r,  $falu, "");
    $snps_l_ctl  = generate_element_select("filter_snps_l",  $r,  $fsnl, "");
    $snps_u_ctl  = generate_element_select("filter_snps_u",  $r,  $fsnu, "");
    $r = range(1, 9);
    $r = array_merge($r, range(10, 500, 10));
    array_push($r, 1000, 2000, 10000);
    $pare_l_ctl  = generate_element_select("filter_pare_l",  $r, $fpal, "");
    $pare_u_ctl  = generate_element_select("filter_pare_u",  $r, $fpau, "");
    $chr_ctl   = generate_element_select("filter_chr",   $chrs, $fch, "");
    $sbp_ctl   = generate_element_select("filter_sbp",   range(0, $max_chr_len), $fsb, "");
    $ebp_ctl   = generate_element_select("filter_ebp",   range(0, $max_chr_len), $feb, "");
    $ref_ctl   = generate_element_select("filter_ref",   array("exon", "intron", "genomic"), $ref, "");
    $r = range(0, 9);
    $r = array_merge($r, range(10, 100, 5));
    $r = array_merge($r, range(200, 500, 100));
    for ($i = 0; $i < count($r); $i++)
      $r[$i] = $r[$i] * -1;
    $lnl_l_ctl = generate_element_select("filter_lnl_l",  $r,  $flnl, "");
    $lnl_u_ctl = generate_element_select("filter_lnl_u",  $r,  $flnu, "");

    if (isset($display['filter_type'])) {

	foreach ($filters as $key => $f)
	    if (in_array($key, $display['filter_type'])) {
		$filters[$key]['sel'] = "checked=\"checked\"";
		$filters[$key]['tr']  = "class=\"active_filter\"";
	    } else {
		$filters[$key]['sel'] = "";
		$filters[$key]['tr']  = "";
	    }

    } else {
	$filters['none']['sel'] = "checked=\"checked\"";
    }

    echo <<< EOQ
<h4 class="info_head">
  <img id="stacks_filter_img" src="$img_path/caret-d.png" />
  <a onclick="toggle_div('stacks_filter', '$img_path', 'page_state');">Filter Results By</a>
</h4>
<div id="stacks_filter">
<form id="filter_results" name="filter_results" method="get" action="$root_path/catalog.php">
$hidden_vars

<table style="width: 100%; vertical-align: top;">
<tr>
EOQ;

    if (count($chrs) > 0) {
      echo <<< EOQ
<td style="width: 25%;">
<table class="loc_filter">
<tr>
  <td {$filters['loc']['tr']}>
      <input type="checkbox" name="filter_type[]" value="loc" onchange="rebuild_display_select()" {$filters['loc']['sel']} /> 
      <a onclick="toggle_cb('filter_results', 'loc')">
      <acronym title="Filter by genomic location">Location</acronym>:</a></td>
  <td {$filters['loc']['tr']}>
      $chr_ctl
  </td>
</tr>
<tr>
  <td {$filters['loc']['tr']}><span style="padding-left: 1.5em;">Start:</span></td>
  <td {$filters['loc']['tr']}>
      $sbp_ctl Mb
  </td>
</tr>
<tr>
  <td {$filters['loc']['tr']}><span style="padding-left: 1.5em;">End:</span></td>
  <td {$filters['loc']['tr']}>
      $ebp_ctl Mb
  </td>
</tr>
<tr>
  <td {$filters['ref']['tr']}>
      <input type="checkbox" name="filter_type[]" value="ref" onchange="rebuild_display_select()" {$filters['ref']['sel']} /> 
      <a onclick="toggle_cb('filter_results', 'ref')">
      <acronym title="Filter by type of RAD locus">Type</acronym>:</a></td>
  <td {$filters['ref']['tr']}>
      $ref_ctl
  </td>
</tr>
</table>
</td>

EOQ;
    }

    $cat_id_filter = isset($display['filter_cata']) ? $display['filter_cata'] : "";

    echo <<< EOQ
<td>
<table class="filter">
<tr>
  <td colspan="2" {$filters['cata']['tr']}>
      <input type="checkbox" name="filter_type[]" value="cata" onchange="rebuild_display_select()" {$filters['cata']['sel']} /> 
      <a onclick="toggle_cb('filter_results', 'cata')">
      <acronym title="Show a locus with a particular ID.">Catalog ID</acronym>:</a>
      <input name="filter_cata" value="$cat_id_filter" size="15" />
  </td>
  <td style="text-align: right; padding-right: 10px;">
      <input type="submit" value="filter" />
  </td>
</tr>
<tr>
  <td style="width: 40%; text-align: left;">
  <table style="text-align: left;">
  <tr>
  <td {$filters['alle']['tr']}>
      <input type="checkbox" name="filter_type[]" value="alle" onchange="rebuild_display_select()" {$filters['alle']['sel']} /> 
      <a onclick="toggle_cb('filter_results', 'alle')">
      <acronym title="Filter the catalog according to the number of alleles identified for a locus.">Alleles</acronym>:</a>
  </td>
  <td {$filters['alle']['tr']}>
$alle_l_ctl $alle_u_ctl
  </td>
  </tr>
  <tr>
  <td {$filters['snps']['tr']}>
      <input type="checkbox" name="filter_type[]" value="snps" onchange="rebuild_display_select()" {$filters['snps']['sel']} /> 
      <a onclick="toggle_cb('filter_results', 'snps')">
      <acronym title="Filter the catalog according to the number of SNPs found at a locus.">SNPs</acronym>:</a>
  </td>
  <td {$filters['snps']['tr']}>
$snps_l_ctl $snps_u_ctl
  </td>
  </tr>
  <tr>
    <td>&nbsp;</td>
    <td>&nbsp;</td>
  </tr>
  </table>
  </td>

  <td style="width: 45%; text-align: center;">
  <table style="text-align: left;">
  <tr>
  <td {$filters['pare']['tr']}>
      <input type="checkbox" name="filter_type[]" value="pare" onchange="rebuild_display_select()" {$filters['pare']['sel']} /> 
      <a onclick="toggle_cb('filter_results', 'pare')">
      <acronym title="Filter the catalog according to the number of population samples that are matched to a locus.">Matching samples</acronym>:</a>
  </td>
  <td {$filters['pare']['tr']}>
$pare_l_ctl $pare_u_ctl
  </td>
  </tr>
  <tr>
<td {$filters['lnl']['tr']}>
      <input type="checkbox" name="filter_type[]" value="lnl" onchange="rebuild_display_select()" {$filters['lnl']['sel']} /> 
      <a onclick="toggle_cb('filter_results', 'lnl')">
      <acronym title="The mean log likelihood of each catalog locus.">LnL</acronym>:</a>
  </td>
  <td {$filters['lnl']['tr']}>
$lnl_l_ctl $lnl_u_ctl
  </td>
  </tr>
  <tr>
    <td>&nbsp;</td>
    <td>&nbsp;</td>
  </tr>
  </table>
  </td>

  <td style="width: 15%; text-align: right;">
  <table style="text-align: left;">
  <tr>
    <td>&nbsp;</td>
    <td>&nbsp;</td>
  </tr>
  <tr>
    <td>&nbsp;</td>
    <td>&nbsp;</td>
  </tr>

EOQ;

if ($cols['gcnt'] == true) {
  echo <<< EOQ
  <tr>
  <td {$filters['gcnt']['tr']}>
  <td {$filters['gcnt']['tr']}>
$gcnt_ctl
  </td>
  </tr>

EOQ;

} else {
  echo <<< EOQ
  <tr>
    <td>&nbsp;</td>
    <td>&nbsp;</td>
  </tr>

EOQ;
}

echo <<< EOQ
  </table>
  </td>
</tr>

EOQ;

  if ($cols['seq'] == true)
    echo <<< EOQ
<tr>
  <td {$filters['est']['tr']}>
      <input type="checkbox" name="filter_type[]" value="est" onchange="rebuild_display_select()" {$filters['est']['sel']} /> 
      <a onclick="toggle_cb('filter_results', 'est')">
      <acronym title="Filter the catalog to show loci for which ESTs have been associated.">Contains ESTs</acronym></a>
  </td>
  <td {$filters['pe']['tr']}>
      <input type="checkbox" name="filter_type[]" value="pe" onchange="rebuild_display_select()" {$filters['pe']['sel']} /> 
      <a onclick="toggle_cb('filter_results', 'pe')">
      <acronym title="Filter the catalog to show loci for which paired-end RAD-Tags have been associated.">Contains Paired-end RAD-Tags</acronym></a>
  </td>
  <td {$filters['blast']['tr']}>
      <input type="checkbox" name="filter_type[]" value="blast" onchange="rebuild_display_select()" {$filters['blast']['sel']} /> 
      <a onclick="toggle_cb('filter_results', 'blast')">
      <acronym title="Filter the catalog to show loci for which BLAST hits  have been associated.">Contains BLAST Hits</acronym></a>
  </td>
</tr>

EOQ;
    echo <<< EOQ
</table>
</td>
</tr>
</table>
</form>
</div>

EOQ;

}

function process_filter(&$display_params) {

    if (!isset($_GET['filter_type'])) 
	return;

    foreach ($_GET['filter_type'] as $filter) {
	array_push($display_params['filter_type'], $filter);

	if ($filter == "alle") {
	    $display_params['filter_alle_l'] = $_GET['filter_alle_l'];
	    $display_params['filter_alle_u'] = $_GET['filter_alle_u'];

	} else if ($filter == "snps") {
	    $display_params['filter_snps_l'] = $_GET['filter_snps_l'];
	    $display_params['filter_snps_u'] = $_GET['filter_snps_u'];

	} else if ($filter == "pare") {
	    $display_params['filter_pare_l'] = $_GET['filter_pare_l'];
	    $display_params['filter_pare_u'] = $_GET['filter_pare_u'];

	} else if ($filter == "lnl") {
	    $display_params['filter_lnl_l'] = $_GET['filter_lnl_l'];
	    $display_params['filter_lnl_u'] = $_GET['filter_lnl_u'];

	} else if ($filter == "prog") {
	    $display_params['filter_prog'] = $_GET['filter_prog'];

	} else if ($filter == "vprog") {
	    $display_params['filter_vprog'] = $_GET['filter_vprog'];

	} else if ($filter == "cata") {
	    $display_params['filter_cata'] = $_GET['filter_cata'];

	} else if ($filter == "mark") {
	    $display_params['filter_mark'] = $_GET['filter_mark'];

	} else if ($filter == "gcnt") {
	    $display_params['filter_gcnt'] = $_GET['filter_gcnt'];

	} else if ($filter == "chisq") {
	    $display_params['filter_chisq_l'] = $_GET['filter_chisq_l'];
	    $display_params['filter_chisq_u'] = $_GET['filter_chisq_u'];

	} else if ($filter == "ref") {
	    $display_params['filter_ref'] = $_GET['filter_ref'];

	} else if ($filter == "loc") {
	    $display_params['filter_chr'] = $_GET['filter_chr'];
	    $display_params['filter_sbp'] = $_GET['filter_sbp'];
	    $display_params['filter_ebp'] = $_GET['filter_ebp'];
	}
    }
}

function prepare_filter_parameters($display_params, &$param) {
    $filters = $display_params['filter_type'];

    if (!isset($filters))
	return;

    foreach ($filters as $filter) {

	if ($filter == "snps") {
	    array_push($param, $display_params['filter_snps_l']);
	    array_push($param, $display_params['filter_snps_u']);

	} else if ($filter == "alle") {
	    array_push($param, $display_params['filter_alle_l']);
	    array_push($param, $display_params['filter_alle_u']);

	} else if ($filter == "pare") {
	    array_push($param, $display_params['filter_pare_l']);
	    array_push($param, $display_params['filter_pare_u']);

	} else if ($filter == "lnl") {
	    array_push($param, $display_params['filter_lnl_l']);
	    array_push($param, $display_params['filter_lnl_u']);

	} else if ($filter == "prog") {
	    array_push($param, $display_params['filter_prog']);
	
	} else if ($filter == "vprog") {
	    array_push($param, $display_params['filter_vprog']);
	
	} else if ($filter == "cata") {
	    array_push($param, $display_params['filter_cata']);
	
	} else if ($filter == "est") {
	    array_push($param, 0);
	
	} else if ($filter == "pe") {
	    array_push($param, 0);
	
	} else if ($filter == "blast") {
	    array_push($param, 0);
	
	} else if ($filter == "gcnt") {
	    array_push($param, $display_params['filter_gcnt']);
	
	} else if ($filter == "chisq") {
	    array_push($param, $display_params['filter_chisq_l']);
	    array_push($param, $display_params['filter_chisq_u']);
	
	} else if ($filter == "ref") {
	    array_push($param, $display_params['filter_ref']);
	
	} else if ($filter == "loc") {
	    array_push($param, $display_params['filter_chr']);
	    array_push($param, $display_params['filter_sbp'] * 1000000);
	    array_push($param, $display_params['filter_ebp'] * 1000000);
	
	} else if ($filter == "mark") {
	  if ($display_params['filter_mark'] == "Any") 
	    array_push($param, "%/%");
	  else 
	    array_push($param, $display_params['filter_mark']);
	
	}
    }
}

function apply_query_filters($display_params) {
    $order = 0;
    $query = "";
    $sql_filters =
	array("cata"  => "(catalog_index.tag_id = ?)", 
	      "alle"  => "(alleles >= ? AND alleles <= ?)", 
	      "snps"  => "(snps >= ? AND snps <= ?)",
	      "pare"  => "(parents >= ? AND parents <= ?)",
	      "prog"  => "(progeny >= ?)",
	      "vprog" => "(valid_progeny >= ?)",
	      "lnl"   => "(lnl >= ? AND lnl <= ?)",
	      "mark"  => "(marker LIKE ?)", 
              "est"   => "(ests > ?)",
              "pe"    => "(pe_radtags > ?)",
              "blast" => "(blast_hits > ?)",
	      "gcnt"  => "(geno_cnt >= ?)",
	      "chisq" => "(chisq_pval >= ? AND chisq_pval <= ?)",
	      "ref"   => "(catalog_index.type = ?)",
	      "loc"   => "(catalog_index.chr = ? && catalog_index.bp >= ? && catalog_index.bp <= ?)");

    $filters = $display_params['filter_type'];

    if (count($filters) > 0) {
	$query = " AND ";

	while (count($filters) > 0) {
	    $filter = array_shift($filters);
	    $query .= $sql_filters[$filter];
	    $query .= count($filters) > 0 ? " AND " : "";

	    if ($filter == "loc") $order++;
	}

	if ($order)
	  $query .= " ORDER BY chr, bp";
    }

    return $query;
}

function fetch_chrs(&$max_len) {
    global $db, $batch_id;

    $max_len = 0;

    $chrs = array();
    $res = $db['chrs_sth']->execute($batch_id);
    check_db_error($res, __FILE__, __LINE__);

    while ($row = $res->fetchRow()) {
      if ($row['max_len'] > $max_len) 
	$max_len = $row['max_len'];

      array_push($chrs, $row['chr']);
    }

    return $chrs;
}

?>
