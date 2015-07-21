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

$database  = isset($_GET['db'])        ? $_GET['db']        : "";
$batch_id  = isset($_GET['batch_id'])  ? $_GET['batch_id']  : 0;
$sample_id = isset($_GET['sample_id']) ? $_GET['sample_id'] : 0;
$page      = isset($_GET['p'])         ? $_GET['p']         : 1;
$per_page  = isset($_GET['pp'])        ? $_GET['pp']        : 10;

// Connect to the database
$db = db_connect($database);

// Save these variables for automatic URL formation later on.
$display = array();
$display['batch_id']    = $batch_id;
$display['sample_id']   = $sample_id;
$display['db']          = $database;
$display['p']           = $page;
$display['pp']          = $per_page;
$display['filter_type'] = array();

//
// Process the filtering parameters
//
$param = array($batch_id, $sample_id);
process_filter($display);
prepare_filter_parameters($display, $param);

//
// Prepare some SQL queries
//
$query = 
    "SELECT batches.id as id, date, description, samples.id as sample_id, file FROM batches " . 
    "JOIN samples ON (batch_id=batches.id) " . 
    "WHERE batches.id=? AND samples.id=?";
$db['batch_sth'] = $db['dbh']->prepare($query);
check_db_error($db['batch_sth'], __FILE__, __LINE__);

$query = 
    "SELECT COUNT(tag_id) as count FROM tag_index " . 
    "WHERE batch_id=? AND sample_id=?";
$query .= apply_query_filters($display);
$db['count_sth'] = $db['dbh']->prepare($query);
check_db_error($db['count_sth'], __FILE__, __LINE__);

$query = 
    "SELECT depth as count FROM tag_index " . 
    "WHERE batch_id=? AND sample_id=? AND tag_id=?";
$db['depth_sth'] = $db['dbh']->prepare($query);
check_db_error($db['depth_sth'], __FILE__, __LINE__);

$query = 
    "SELECT col, rank_2 FROM snps " . 
    "JOIN samples ON (snps.sample_id=samples.id) " . 
    "JOIN batches ON (samples.batch_id=batches.id) " . 
    "WHERE snps.type='E' AND batch_id=? AND snps.sample_id=? AND tag_id=? ORDER BY col";
$db['snp_sth'] = $db['dbh']->prepare($query);
check_db_error($db['snp_sth'], __FILE__, __LINE__);

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
$batch['sample_id'] = $row['sample_id'];
$batch['file']      = $row['file'];

$page_title = "RAD-Tag Sample Viewer";
write_header($page_title, $batch);

echo <<< EOQ
<form id="page_state" name="page_state"></form>
<h3 style="margin-left: 1em;">
<a href="$root_path/samples.php?db=$database&id=$batch[id]">
Batch #$batch[id] <span class="s">[$batch[date]; $batch[desc]]</span></a></h3>

EOQ;

write_filter();

echo <<< EOQ
<h4 class="info_head" style="margin-left: 1em;">
  <img id="sources_img" src="$img_path/caret-d.png" />
  <a onclick="toggle_div('sources', '$img_path', 'page_state');">
  RAD-Tag Sample #$batch[sample_id] [<span class="s">$batch[file]</span>]</a>
</h4>

<div id="sources" style="width: 100%;">
<a name="results_top"></a>
<table class="db" style="width: 95%; border: none;">
<tr>
  <td colspan="5" style="border: none; padding-bottom: 0px;">

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

write_pagination($pagination_count, $start_group, $end_group, "tags.php");

echo <<< EOQ
  </td>
</tr>
<tr>
  <th style="width: 5%;">Id</th>
  <th style="width: 10%;">Depth</th>
  <th style="width: 10%;">SNP</th>
  <th style="width: 55%;">Consensus</th>
  <th style="width: 10%;">Catalog ID</th>
</tr>

EOQ;

$db['dbh']->setLimit($display['pp'], $start_group - 1);
check_db_error($db['dbh'], __FILE__, __LINE__);

$query = 
    "SELECT con_tag_id as id, tag_index.tag_id as tag_id, " . 
    "depth, seq, catalog_id, " . 
    "tag_index.deleveraged, tag_index.removed, tag_index.blacklisted " . 
    "FROM tag_index " . 
    "JOIN unique_tags ON (con_tag_id=unique_tags.id) " . 
    "WHERE tag_index.batch_id=? AND tag_index.sample_id=?";
$query .= apply_query_filters($display);

$db['tag_sth'] = $db['dbh']->prepare($query);
check_db_error($db['tag_sth'], __FILE__, __LINE__);
$result = $db['tag_sth']->execute($param);
check_db_error($result, __FILE__, __LINE__);

while ($row = $result->fetchRow()) {

    // Query the database to find how many SNPs were found in this sample.
    $snps = array();
    $snp_res = $db['snp_sth']->execute(array($batch_id, $sample_id, $row['tag_id']));
    check_db_error($snp_res, __FILE__, __LINE__);
    while ($snp_row = $snp_res->fetchRow()) {
      array_push($snps, array('col' => $snp_row['col'], 'rank' => $snp_row['rank_2']));
    }

    print
	"<tr>\n" .
	"  <td><a href=\"$root_path/tag.php?db=$database&batch_id=$batch_id&sample_id=$sample_id&tag_id=$row[tag_id]&p=$display[p]&pp=$display[pp]\">" . 
	"$row[tag_id]</a></td>\n" .
	"  <td>$row[depth]</td>\n";

    if (count($snps) == 0)
	print "  <td>No</td>\n";
    else
	print "  <td>Yes <span class=\"s\">[" . count($snps) . "nuc]</span></td>\n";

    $s = print_snps($row['tag_id'], $row['seq'], $row['seq'], $snps, true);

    print
      "  <td class=\"seq\"><div class=\"seq\">" . $s . "</div></td>\n";

    $row['catalog_id'] > 0 ?
        print 
	  "  <td><a href=\"$root_path/catalog.php?db=$database&id=$batch_id&filter_type[]=cata&filter_cata=$row[catalog_id]\">$row[catalog_id]</a></td>\n" :
        print 
          "  <td><span style=\"color: #888888; font-weight: bold;\">absent</span></td>\n";
    print "</tr>\n";
}

print 
"<tr>\n" .
"  <td colspan=\"5\" style=\"border: none; padding-top: 0px;\">\n";

write_pagination($pagination_count, $start_group, $end_group, "tags.php");

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

function generate_url($destination) {
    global $root_path, $display;

    $url = "href=\"" . $root_path . "/" . $destination . "?";

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

    $url .= "\"";

    return $url;
}

function generate_page_list($page, $num_pages, $destination) {
    global $display;

    $page_list = "";

    if ($page <= 4) {
	for ($i = 1; $i < $page; $i++) {
	    $display['p'] = $i;
	    $p            = generate_url($destination);
	    $page_list   .= "<a $p>$i</a>\n"; 
	} 
    } else {
	$display['p'] = 1;
	$p            = generate_url($destination);
	$page_list   .= "<a $p>1</a> ...\n"; 

	foreach (array($page - 3, $page - 2, $page - 1) as $i) {
	    $display['p'] = $i;
	    $p            = generate_url($destination);
	    $page_list   .= "<a $p>$i</a>\n"; 
	} 
    }

    $page_list .= " <strong>$page</strong>\n";

    if ($page <= $num_pages - 4) {
	for ($i = $page + 1; $i <= $page + 3; $i++) {
	    $display['p'] = $i;
	    $p            = generate_url($destination);
	    $page_list   .= "<a $p>$i</a>\n"; 
	} 

	$display['p'] = $num_pages;
	$p            = generate_url($destination);
	$page_list   .= "... <a $p>$num_pages</a>\n"; 

    } else {
	for ($i = $page + 1; $i <= $num_pages; $i++) {
	    $display['p'] = $i;
	    $p            = generate_url($destination);
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
    $prev_page     = generate_url($destination);
    $display['p'] += 2;
    $next_page     = generate_url($destination);
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

function write_filter() {
    global $img_path, $root_path, $display;

    $hidden_vars  = generate_hidden_form_vars("filter");

    $filters = array("depth" => array(),
		     "snps"  => array(),
		     "tagid" => array(),
		     "black" => array(),
		     "delv"  => array(),
		     "rem"   => array());

    $ele_name  = isset($display['filter_depth']) ? $display['filter_depth'] : "";
    $depth_ctl = generate_element_select("filter_depth", array(1, 5, 10, 20), $ele_name, "");

    $ele_name  = isset($display['filter_snps']) ? $display['filter_snps'] : "";
    $snps_ctl  = generate_element_select("filter_snps",  array(1, 2, 3, 4, 5), $ele_name, "");

    $ele_name  = isset($display['filter_delv']) ? $display['filter_delv'] : "";
    $delv_ctl  = generate_key_element_select("filter_delv",  array(1 => "True", 0 => "False"), $ele_name, "");

    $ele_name  = isset($display['filter_rem']) ? $display['filter_rem'] : "";
    $rem_ctl   = generate_key_element_select("filter_rem",   array(1 => "True", 0 => "False"), $ele_name, "");

    $ele_name  = isset($display['filter_black']) ? $display['filter_black'] : "";
    $black_ctl = generate_key_element_select("filter_black", array(1 => "True", 0 => "False"), $ele_name, "");

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

    $tagid = isset($display['filter_tagid']) ? $display['filter_tagid'] : "";

    echo <<< EOQ
<h4 class="info_head">
  <img id="filter_img" src="$img_path/caret-d.png" />
  <a onclick="toggle_div('filter', '$img_path', 'page_state');">Filter Results</a>
</h4>
<div class="filter">
<form id="filter_results" name="filter_results" method="get" action="$root_path/tags.php">
$hidden_vars
<table class="filter">
<tr {$filters['tagid']['tr']}>
  <td><input type="checkbox" name="filter_type[]" value="tagid" onchange="rebuild_display_select()" {$filters['tagid']['sel']} /> 
      <a onclick="toggle_cb('filter_results', 'tagid')">Filter by Tag ID:</a></td>
  <td>
    <input name="filter_tagid" value="$tagid" size="15" />
  </td>
</tr>
<tr {$filters['depth']['tr']}>
  <td><input type="checkbox" name="filter_type[]" value="depth" onchange="rebuild_display_select()" {$filters['depth']['sel']} /> 
      <a onclick="toggle_cb('filter_results', 'depth')">Filter by depth of coverage:</a></td>
  <td>
$depth_ctl
  </td>
</tr>
<tr {$filters['snps']['tr']}>
  <td><input type="checkbox" name="filter_type[]" value="snps" onchange="rebuild_display_select()" {$filters['snps']['sel']} /> 
      <a onclick="toggle_cb('filter_results', 'snps')">Filter by presence of SNPs:</a></td>
  <td>
$snps_ctl
  </td>
</tr>
<tr {$filters['delv']['tr']}>
  <td><input type="checkbox" name="filter_type[]" value="delv" onchange="rebuild_display_select()" {$filters['delv']['sel']} /> 
      <a onclick="toggle_cb('filter_results', 'delv')">View deleveraged tags:</a></td>
  <td>
$delv_ctl
  </td>
</tr>
<tr {$filters['black']['tr']}>
  <td><input type="checkbox" name="filter_type[]" value="black" onchange="rebuild_display_select()" {$filters['black']['sel']} /> 
      <a onclick="toggle_cb('filter_results', 'black')">View blacklisted tags:</a></td>
  <td>
$black_ctl
  </td>
</tr>
<tr {$filters['rem']['tr']}>
  <td><input type="checkbox" name="filter_type[]" value="rem" onchange="rebuild_display_select()" {$filters['rem']['sel']} /> 
      <a onclick="toggle_cb('filter_results', 'rem')">View lumberjack stacks:</a></td>
  <td>
$rem_ctl
  </td>
</tr>
<tr>
  <td colspan="2" style="text-align: right; padding-right: 10px;">
    <input type="submit" value="filter" onclick="update_page_state_form(this.form.id, 'page_state')" />
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

	if ($filter == "depth") {
	    $display_params['filter_depth'] = $_GET['filter_depth'];

	} else if ($filter == "snps") {
	    $display_params['filter_snps'] = $_GET['filter_snps'];

	} else if ($filter == "tagid") {
	    $display_params['filter_tagid'] = $_GET['filter_tagid'];

	} else if ($filter == "delv") {
	    $display_params['filter_delv'] = $_GET['filter_delv'];

	} else if ($filter == "black") {
	    $display_params['filter_black'] = $_GET['filter_black'];

	} else if ($filter == "rem") {
	    $display_params['filter_rem'] = $_GET['filter_rem'];

	}
    }
}

function prepare_filter_parameters($display_params, &$param) {
    $filters = $display_params['filter_type'];

    if (!isset($filters))
	return;

    foreach ($filters as $filter) {

	if ($filter == "snps") {
	    array_push($param, $display_params['filter_snps']);

	} else if ($filter == "depth") {
	    array_push($param, $display_params['filter_depth']);
	
	} else if ($filter == "tagid") {
	    array_push($param, $display_params['filter_tagid']);

	} else if ($filter == "delv") {
	    array_push($param, $display_params['filter_delv']);

	} else if ($filter == "black") {
	    array_push($param, $display_params['filter_black']);

	} else if ($filter == "rem") {
	    array_push($param, $display_params['filter_rem']);
	}
    }
}

function apply_query_filters($display_params) {
    $sql_filters =
	array("depth" => "(depth >= ?)", 
	      "snps"  => "(snps >= ?)",
	      "tagid" => "(tag_index.tag_id = ?)",
	      "delv"  => "(tag_index.deleveraged = ?)",
	      "black" => "(tag_index.blacklisted = ?)",
	      "rem"   => "(tag_index.removed = ?)");

    $filters = $display_params['filter_type'];
    $query   = "";

    if (count($filters) > 0) {
	$query = " AND ";

	while (count($filters) > 0) {
	    $filter = array_shift($filters);
	    $query .= $sql_filters[$filter];
	    $query .= count($filters) > 0 ? " AND " : "";
	}
    }

    return $query;
}

?>
