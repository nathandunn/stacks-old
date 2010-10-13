<?php
require_once("header.php");
require_once("radtag_functions.php");

$batch_id  = isset($_GET['id']) ? $_GET['id'] : 0;
$database  = isset($_GET['db']) ? $_GET['db'] : "radtags";
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
    "SELECT batches.id as id, date, description FROM batches " . 
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

$page_title = "RAD-Tag Catalog Viewer";
write_header($page_title, $batch);

echo <<< EOQ
<form id="page_state" name="page_state"></form>
<h3 style="margin-left: 1em;">
Catalog RAD-Tags</h3>

EOQ;

write_filter();

// Generate Excel export URL
$excel_export = generate_url("export.php");

echo <<< EOQ
<h4 class="info_head" style="margin-left: 1em;">
  <img id="sources_img" src="/acos/images/caret-d.png" />
  <a href="$root_path/index.php?db=$database">
  Batch #$batch[id] <span class="s">[$batch[date]; $batch[desc]]</span></a></h3>
</h4>

<div id="sources" style="width: 100%;">
<a name="results_top"></a>
<table class="db" style="width: 100%; border: none;">
<tr>
  <td colspan="8" style="border: none; padding-bottom: 0px;">
    <a $excel_export>
    <img style="float: right;" src="$root_path/images/excel_icon.png" title="Export data to Microsoft Excel format"/></a>
  </td>
</tr>
<tr>
  <td colspan="8" style="border: none; padding-bottom: 0px;">

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

echo <<< EOQ
  </td>
</tr>
<tr>
  <th style="width: 10%;">Id</th>
  <th style="width: 5%;">SNP</th>
  <th style="width: 50%;">Consensus</th>
  <th style="width: 5%;">Matching Parents</th>
  <th style="width: 5%;">Progeny</th>
  <th style="width: 5%;">Marker</th>
  <th style="width: 10%;">Ratio</th>
  <th style="width: 10%">Sequence</th>
</tr>

EOQ;

$db['dbh']->setLimit($display['pp'], $start_group - 1);
check_db_error($db['dbh'], __FILE__, __LINE__);

$query = 
    "SELECT catalog_index.tag_id as tag_id, alleles, parents, progeny, valid_progeny, " . 
    "seq, marker, max_pct, ratio, ests, pe_radtags, blast_hits, external_id " .
    "FROM catalog_index " .
    "JOIN catalog_tags ON (catalog_index.cat_id=catalog_tags.id) " . 
    "LEFT JOIN catalog_annotations ON (catalog_index.batch_id=catalog_annotations.batch_id AND catalog_index.tag_id=catalog_annotations.catalog_id) " .
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

    $url = "$root_path/catalog_tag.php?db=$database&batch_id=$batch_id&tag_id=$row[tag_id]";
    $annotation = strlen($row['external_id']) > 0 ? $row['external_id'] : "annotate";

    echo <<< EOQ
<tr>
  <td class="catlink">
<img id="{$row[tag_id]}_img" src="$img_path/caret-u.png" />
<a onclick="toggle_aln_tr('$row[tag_id]', '$img_path', '$url');">$row[tag_id]</a><br />

<a class="annotation" id="{$row[tag_id]}_ann" onclick="toggle_annotation($row[tag_id])">$annotation</a>
<div id="{$row[tag_id]}_div" style="display: none; font-size: small;">
<form id="{$row[tag_id]}_frm">
<input type="hidden" name="url" value="$root_path/annotate_marker.php?db=$database&batch_id=$batch_id&tag_id=$row[tag_id]" />
<input type="input" size=15 name="ext_id" value="" />
<a onclick="annotate_marker('$row[tag_id]')">save</a>|<a onclick="toggle_annotation('$row[tag_id]')">cancel</a>
</form>
</div>
  </td>

EOQ;

    if (count($snps) == 0)
	print "  <td>No</td>\n";
    else
	print "  <td>Yes <span class=\"s\">[" . count($snps) . "nuc]</span></td>\n";

    $s = print_snps($row['seq'], $row['seq'], $snps);

    $ratio = explode(";", $row['ratio']);
    $ratio_parsed = "";
    $i            = 0;
    foreach ($ratio as $r) {
      if (strlen($r) == 0) continue;

      preg_match("/([abcd]{1,2}):(\d+)\((\d+\.?\d*%)\)/", $r, $matches);

      $color = $colors[$i % $color_size];

      $ratio_parsed .= "<span style=\"color: $color;\">$matches[1]: $matches[2]";

      if ($matches[3] > 0)
	$ratio_parsed .= " ($matches[3])";

      $ratio_parsed .= "</span><br />";
      $i++;
    }

    if (strlen($ratio_parsed) > 0) {
        $url = "$root_path/catalog_genotypes.php?db=$database&batch_id=$batch_id&tag_id=$row[tag_id]";
        $ratio_parsed =
            $ratio_parsed .
            "<div class=\"catlink\"><img id=\"{$row[tag_id]}_gtypes_img\" src=\"$img_path/caret-u.png\" />" .
            "<a onclick=\"toggle_aln_tr('{$row[tag_id]}_gtypes', '$img_path', '$url');\">" .
            "view genotypes</a></div>";
    }

    $url = "$root_path/sequence_blast.php?db=$database&batch_id=$batch_id&tag_id=$row[tag_id]";
    if ($row[blast_hits] > 0) {
        $blast_hits_str =
            "<div class=\"catlink\"><img id=\"{$row[tag_id]}_blast_img\" src=\"$img_path/caret-u.png\" />" .
            "<a onclick=\"toggle_aln_tr('{$row[tag_id]}_blast', '$img_path', '$url');\">" .
            "blast hits: $row[blast_hits]</a>";
    } else {
        $blast_hits_str = "blast hits: $row[blast_hits]";
    }

    echo <<< EOQ
  <td class="seq">$s</td>
  <td>$row[parents]</td>
  <td><span title="Matching Progeny">$row[progeny]</span> <strong>/</strong> <span title="Mappable Progeny">$row[valid_progeny]</span></td>
  <td>$row[marker]</td>
  <td style="text-align: left; font-size: smaller;">
    $ratio_parsed
  </td>
  <td>
    <table class="int" style="font-size: smaller;">
    <tr><td>ests: $row[ests]</td><td>pe: $row[pe_radtags]</td></tr>
    <tr>
      <td colspan="2">
      $blast_hits_str
      </td></tr>
    </table>
  </td>
</tr>
<tr id="{$row[tag_id]}" style="display: none">
  <td colspan="8">
    <iframe id="{$row[tag_id]}_iframe" 
            frameborder="0" 
            scrolling="no" 
            onload="this.style.height = this.contentDocument.height + 'px';" 
            src=""></iframe>
  </td>
</tr>
<tr id="{$row[tag_id]}_gtypes" style="display: none">
  <td colspan="8">
    <iframe id="{$row[tag_id]}_gtypes_iframe" 
            frameborder="0" 
            scrolling="no" 
            onload="this.style.height = (this.contentDocument.height+25) + 'px';" 
            src=""></iframe>
  </td>
</tr>
<tr id="{$row[tag_id]}_blast" style="display: none">
  <td colspan="8">
    <iframe id="{$row[tag_id]}_blast_iframe" 
            frameborder="0" 
            scrolling="no" 
            onload="this.style.height = (this.contentDocument.height+25) + 'px';" 
            src=""></iframe>
  </td>
</tr>
EOQ;
}

print 
"<tr>\n" .
"  <td colspan=\"8\" style=\"border: none; padding-top: 0px;\">\n";

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

    $filters = array("cata"  => array(),
		     "alle"  => array(),
		     "snps"  => array(),
		     "pare"  => array(),
		     "prog"  => array(),
		     "vprog" => array(),
		     "mark"  => array(),
                     "est"   => array(),
                     "pe"    => array(),
                     "blast" => array());

    $alle_ctl  = generate_element_select("filter_alle", array(1, 2, 3, 4), $display['filter_alle'], "");
    $snps_ctl  = generate_element_select("filter_snps", array(1, 2, 3, 4, 8), $display['filter_snps'], "");
    $pare_ctl  = generate_element_select("filter_pare", array(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18), $display['filter_pare'], "");
    $prog_ctl  = generate_element_select("filter_prog", array(1, 2, 4, 8, 16, 32, 64, 70, 85, 90), $display['filter_prog'], "");
    $vprog_ctl = generate_element_select("filter_vprog", array(1, 2, 4, 8, 16, 32, 64, 70, 85, 90), $display['filter_vprog'], "");
    $mark_ctl  = generate_element_select("filter_mark", 
                                         array('Any', 'aa/bb', 'ab/--', '--/ab', 'aa/ab', 'ab/aa', 'ab/ab', 'ab/ac', 'ab/cd', 'ab/cc', 'cc/ab'), 
                                         $display['filter_mark'], "");

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
  <img id="filter_img" src="$img_path/caret-d.png" />
  <a onclick="toggle_div('filter', '$img_path', 'page_state');">Filter Results</a>
</h4>
<div id="filter" $filter_vis>
<form id="filter_results" name="filter_results" method="get" action="$root_path/catalog.php">
$hidden_vars
<table class="filter">
<tr {$filters['cata']['tr']}>
  <td><input type="checkbox" name="filter_type[]" value="cata" onchange="rebuild_display_select()" {$filters['cata']['sel']} /> 
      <a onclick="toggle_cb('filter_results', 'cata')">Filter by catalog ID:</a></td>
  <td>
    <input name="filter_cata" value="$display[filter_cata]" size="15" />
  </td>
</tr>
<tr {$filters['alle']['tr']}>
  <td><input type="checkbox" name="filter_type[]" value="alle" onchange="rebuild_display_select()" {$filters['alle']['sel']} /> 
      <a onclick="toggle_cb('filter_results', 'alle')">Filter by number of alleles:</a></td>
  <td>
$alle_ctl
  </td>
</tr>
<tr {$filters['snps']['tr']}>
  <td><input type="checkbox" name="filter_type[]" value="snps" onchange="rebuild_display_select()" {$filters['snps']['sel']} /> 
      <a onclick="toggle_cb('filter_results', 'snps')">Filter by presence of SNPs:</a></td>
  <td>
$snps_ctl
  </td>
</tr>
<tr {$filters['pare']['tr']}>
  <td><input type="checkbox" name="filter_type[]" value="pare" onchange="rebuild_display_select()" {$filters['pare']['sel']} /> 
      <a onclick="toggle_cb('filter_results', 'pare')">Filter by the number of parental matches:</a></td>
  <td>
$pare_ctl
  </td>
</tr>
<tr {$filters['prog']['tr']}>
  <td><input type="checkbox" name="filter_type[]" value="prog" onchange="rebuild_display_select()" {$filters['prog']['sel']} /> 
      <a onclick="toggle_cb('filter_results', 'prog')">Filter by the number of progeny matches:</a></td>
  <td>
$prog_ctl
  </td>
</tr>
<tr {$filters['vprog']['tr']}>
  <td><input type="checkbox" name="filter_type[]" value="vprog" onchange="rebuild_display_select()" {$filters['vprog']['sel']} /> 
      <a onclick="toggle_cb('filter_results', 'vprog')">Filter by the number of mappable progeny:</a></td>
  <td>
$vprog_ctl
  </td>
</tr>
<tr {$filters['mark']['tr']}>
  <td><input type="checkbox" name="filter_type[]" value="mark" onchange="rebuild_display_select()" {$filters['mark']['sel']} /> 
      <a onclick="toggle_cb('filter_results', 'mark')">Filter by mapable markers:</a></td>
  <td>
$mark_ctl
  </td>
</tr>
<tr {$filters['est']['tr']}>
  <td colspan="2"><input type="checkbox" name="filter_type[]" value="est" onchange="rebuild_display_select()" {$filters['est']['sel']} /> 
      <a onclick="toggle_cb('filter_results', 'est')">Contains ESTs</a></td>
</tr>
<tr {$filters['pe']['tr']}>
  <td colspan="2"><input type="checkbox" name="filter_type[]" value="pe" onchange="rebuild_display_select()" {$filters['pe']['sel']} /> 
      <a onclick="toggle_cb('filter_results', 'pe')">Contains Paired-end RAD-Tags</a></td>
</tr>
<tr {$filters['blast']['tr']}>
  <td colspan="2"><input type="checkbox" name="filter_type[]" value="blast" onchange="rebuild_display_select()" {$filters['blast']['sel']} /> 
      <a onclick="toggle_cb('filter_results', 'blast')">Contains BLAST Hits</a></td>
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

	if ($filter == "alle") {
	    $display_params['filter_alle'] = $_GET['filter_alle'];

	} else if ($filter == "snps") {
	    $display_params['filter_snps'] = $_GET['filter_snps'];

	} else if ($filter == "pare") {
	    $display_params['filter_pare'] = $_GET['filter_pare'];

	} else if ($filter == "prog") {
	    $display_params['filter_prog'] = $_GET['filter_prog'];

	} else if ($filter == "vprog") {
	    $display_params['filter_vprog'] = $_GET['filter_vprog'];

	} else if ($filter == "cata") {
	    $display_params['filter_cata'] = $_GET['filter_cata'];

	} else if ($filter == "mark") {
	    $display_params['filter_mark'] = $_GET['filter_mark'];

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

	} else if ($filter == "alle") {
	    array_push($param, $display_params['filter_alle']);
	
	} else if ($filter == "pare") {
	    array_push($param, $display_params['filter_pare']);
	
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
	
	} else if ($filter == "mark") {
	  if ($display_params['filter_mark'] == "Any") 
	    array_push($param, "%/%");
	  else 
	    array_push($param, $display_params['filter_mark']);
	
	}
    }
}

function apply_query_filters($display_params) {
    $sql_filters =
	array("cata"  => "(catalog_index.tag_id = ?)", 
	      "alle"  => "(alleles >= ?)", 
	      "snps"  => "(snps >= ?)",
	      "pare"  => "(parents = ?)",
	      "prog"  => "(progeny >= ?)",
	      "vprog" => "(valid_progeny >= ?)",
	      "mark"  => "(marker LIKE ?)", 
              "est"   => "(ests > ?)",
              "pe"    => "(pe_radtags > ?)",
              "blast" => "(blast_hits > ?)");

    $filters = $display_params['filter_type'];

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
