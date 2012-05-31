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

$database   = isset($_GET['db'])       ? $_GET['db']       : "";
$tag_id     = isset($_GET['tag_id'])   ? $_GET['tag_id']   : 0;
$batch_id   = isset($_GET['batch_id']) ? $_GET['batch_id'] : 0;
$batch_type = isset($_GET['type'])     ? $_GET['type']     : "map";
$page       = isset($_GET['p'])        ? $_GET['p']        : 1;
$per_page   = isset($_GET['pp'])       ? $_GET['pp']       : 10; 

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
// Prepare the possible select lists we will want to construct to manualy correct genotypes.
//
$marker_types = array('ab/--' => array('aa', 'bb', '-'),
		      '--/ab' => array('aa', 'bb', '-'),
		      'ab/aa' => array('aa', 'ab', '-', 'clr'),
                      'aa/ab' => array('aa', 'ab', '-'),
                      'ab/ab' => array('aa', 'ab', 'bb', '-'),
                      'ab/ac' => array('aa', 'ab', 'ac', 'bc', '-'),
                      'ab/cd' => array('aa', 'bb', 'cc', 'dd', 'ac', 'ad', 'bc', 'bd', '-'),
                      'aa/bb' => array('aa', 'bb', 'ab', '-'),
                      'ab/cc' => array('aa', 'bb', 'ab', 'ac', 'bc', 'cc', '-'),
                      'cc/ab' => array('aa', 'bb', 'ab', 'ac', 'bc', 'cc', '-'));

//
// Prepare some SQL queries
//
$query = 
    "SELECT samples.id, samples.sample_id, samples.type, file, tag_id, allele, depth, pop_id " . 
    "FROM matches " . 
    "JOIN samples ON (matches.sample_id=samples.id) " . 
    "WHERE matches.batch_id=? AND catalog_id=? ORDER BY samples.id";
$db['mat_sth'] = $db['dbh']->prepare($query);
check_db_error($db['mat_sth'], __FILE__, __LINE__);

$query = 
    "SELECT col, rank_1, rank_2, rank_3, rank_4 FROM catalog_snps " . 
    "WHERE batch_id=? AND tag_id=? ORDER BY col";
$db['snp_sth'] = $db['dbh']->prepare($query);
check_db_error($db['snp_sth'], __FILE__, __LINE__);

$query = 
    "SELECT allele FROM catalog_alleles " . 
    "WHERE batch_id=? AND tag_id=? ";
$db['all_sth'] = $db['dbh']->prepare($query);
check_db_error($db['all_sth'], __FILE__, __LINE__);

$query = 
    "SELECT geno_map FROM markers " . 
    "WHERE batch_id=? AND catalog_id=? ";
$db['map_sth'] = $db['dbh']->prepare($query);
check_db_error($db['map_sth'], __FILE__, __LINE__);

$query = 
    "SELECT pop_id, pop_name FROM populations " . 
    "WHERE batch_id=?";
$db['pop_sth'] = $db['dbh']->prepare($query);
check_db_error($db['pop_sth'], __FILE__, __LINE__);

$query = 
    "SELECT col, pop_id, n, p, obs_het, obs_hom, exp_het, exp_hom, pi, fis FROM sumstats " . 
    "WHERE batch_id=? AND tag_id=?";
$db['stats_sth'] = $db['dbh']->prepare($query);
check_db_error($db['stats_sth'], __FILE__, __LINE__);

$query = 
    "SELECT col, pop_id_1, pop_id_2, pi_o, fst, fst_s FROM fst " . 
    "WHERE batch_id=? AND tag_id=?";
$db['fst_sth'] = $db['dbh']->prepare($query);
check_db_error($db['stats_sth'], __FILE__, __LINE__);

$query = 
    "SELECT marker, catalog_genotypes.sample_id, file, " . 
    "catalog_genotypes.genotype, genotype_corrections.genotype as corrected " . 
    "FROM catalog_genotypes " . 
    "LEFT JOIN genotype_corrections ON " . 
    "(genotype_corrections.catalog_id=catalog_genotypes.catalog_id AND " .
    "genotype_corrections.sample_id=catalog_genotypes.sample_id AND " .
    "genotype_corrections.batch_id=catalog_genotypes.batch_id) " .
    "JOIN samples ON (catalog_genotypes.sample_id=samples.id) " . 
    "JOIN catalog_index ON (catalog_genotypes.catalog_id=catalog_index.tag_id AND " . 
    "catalog_genotypes.batch_id=catalog_index.batch_id) " .
    "WHERE catalog_genotypes.batch_id=? and catalog_genotypes.catalog_id=? " . 
    "ORDER BY catalog_genotypes.sample_id";
$db['geno_sth'] = $db['dbh']->prepare($query);
check_db_error($db['geno_sth'], __FILE__, __LINE__);

$page_title = "Catalog RAD-Tag Viewer";
write_compact_header($page_title);

//
// Fetch population names if available.
//
$pop_names = array();
if ($batch_type == "population") {
    $result = $db['pop_sth']->execute($batch_id);
    check_db_error($result, __FILE__, __LINE__);

    while ($row = $result->fetchRow())
        $pop_names[$row['pop_id']] = $row['pop_name'];
}

echo <<< EOQ
<table class="catalog">
<tr>

EOQ;

$result = $db['snp_sth']->execute(array($batch_id, $tag_id));
check_db_error($result, __FILE__, __LINE__);

echo <<< EOQ
  <td style="width: 30%;">
  <div id="snps" class="snps">
  <h3>SNPs</h3>

EOQ;

$i = 0;
while ($row = $result->fetchRow()) {
    $snp = "$row[rank_1]/$row[rank_2]";
    if (strlen($row['rank_3']) > 0) $snp .= "/$row[rank_3]";
    if (strlen($row['rank_4']) > 0) $snp .= "/$row[rank_4]";

    if ($i == 0) {
      print "    <input type=\"hidden\" id=\"snp_index\" value=\"$row[col]\" />\n";
    }
    $class = $i == 0 ? "class=\"snp_sel\"" : "class=\"snp\"";
    print "    <span $class id=\"{$row[col]}_snp\" onclick=\"toggle_sumstats(this, $row[col])\">Column: $row[col]; $snp</span><br />\n";
    $i++;
}

echo <<< EOQ
  </div>
  <div class="alleles">
  <h3>Alleles</h3>

EOQ;

$alleles = array();
$i       = 0;

$result = $db['map_sth']->execute(array($batch_id, $tag_id));
check_db_error($result, __FILE__, __LINE__);

if ($result->numRows() > 0) {
  $row = $result->fetchRow();
} else {
  $row = array();
}

if (isset($row['geno_map'])) {
  $map = array();

  $genos = explode(";", $row['geno_map']);
  foreach ($genos as $g) {
      if (strlen($g) == 0) continue;
      $m = explode(":", $g);
      $map[$m[0]] = $m[1];
      $alleles[$m[0]] = $colors[$i % $color_size];
      $i++;
  }
  asort($map);
  foreach ($map as $hapl => $geno) {
      print 
	"  <div class=\"haplotype_def\" style=\"color: " . $alleles[$hapl] . ";\">" . 
	"<acronym title=\"Genotype\">$geno</acronym> : " . 
	"<acronym title=\"Observed Haplotype\">$hapl</acronym></div>\n";
  }
} else {
    $result = $db['all_sth']->execute(array($batch_id, $tag_id));
    check_db_error($result, __FILE__, __LINE__);

    while ($row = $result->fetchRow()) {
        $alleles[$row['allele']] = $colors[$i % $color_size];

	print "  <span style=\"color: " . $alleles[$row['allele']] . ";\">$row[allele]</span><br />\n";

	$i++;
    }
}

echo <<< EOQ
  </div>
  <div class="sumstats">

EOQ;

if ($batch_type == "population") {
    print_sumstats($db, $pop_names);
    print_fst($db, $pop_names);
}

echo <<< EOQ
  </div>
</td>

EOQ;

$htypes = array();
$gtypes = array();
//
// Fetch and record Observed Haplotypes
//
$result = $db['mat_sth']->execute(array($batch_id, $tag_id));
check_db_error($result, __FILE__, __LINE__);

while ($row = $result->fetchRow()) {
    $a = array('id'     => $row['id'],
               'file'   => $row['file'], 
               'allele' => $row['allele'], 
               'tag_id' => $row['tag_id'],
	       'depth'  => $row['depth'],
	       'pop_id' => $row['pop_id']);
    if (!isset($htypes[$row['pop_id']]))
        $htypes[$row['pop_id']] = array();

    if (!isset($htypes[$row['pop_id']][$row['file']]))
        $htypes[$row['pop_id']][$row['file']] = array();

    array_push($htypes[$row['pop_id']][$row['file']], $a);
}

//
// Fetch and record Genotypes
//
$result = $db['geno_sth']->execute(array($batch_id, $tag_id));
check_db_error($result, __FILE__, __LINE__);

while ($row = $result->fetchRow()) {
    $gtypes[$row['file']] = array('id'        => $row['sample_id'],
				  'file'      => $row['file'],
				  'genotype'  => $row['genotype'],
				  'corrected' => $row['corrected'],
				  'marker'    => $row['marker']);
}

if (count($gtypes) > 0)
  $gtype_str = 
    "<input id=\"gen_cb\" type=\"checkbox\" onclick=\"toggle_genotypes($tag_id, 'locus_gtypes', 'gen')\" /> " .
    "<a onclick=\"document.getElementById('gen_cb').click()\">Genotypes</a>";
else
  $gtype_str = "";

echo <<< EOQ
<td class="matches" style="text-align: right;">

<table id="locus_gtypes" class="genotypes">
<tr>
  <td colspan="10" class="gtype_toggle">
    <strong>View:</strong>
    <input type="checkbox" id="hap_cb" checked="checked" onclick="toggle_genotypes($tag_id, 'locus_gtypes', 'hap')" />
    <a onclick="document.getElementById('hap_cb').click()">Haplotypes</a>
    <input type="checkbox" id="dep_cb" onclick="toggle_genotypes($tag_id, 'locus_gtypes', 'dep')" /> 
    <a onclick="document.getElementById('dep_cb').click()">Allele Depths</a>
    $gtype_str
  </td>
</tr>

EOQ;

$i        = 0;
$num_pops = count($htypes);
$num_cols = 0;

foreach ($htypes as $pop_id => $population)
  $num_cols = count($population) > $num_cols ? count($population) : $num_cols;
$num_cols = $num_cols < 10 ? $num_cols : 10;

ksort($htypes);

foreach ($htypes as $pop_id => $population) {
    print "<tr>\n";
    if ($num_pops > 1) {
      print 
	"  <td class=\"pop_id\" colspan=\"$num_cols\">";
      print_population_name($database, $batch_id, $pop_id, $pop_names);
      print 
	"</td>\n" .
	"</tr>\n";
    }

    foreach ($population as $sample => $match) {
        $i++;

	print 
	  "  <td>" .
	  "<span class=\"title\">" . ucfirst(str_replace("_", " ", $match[0]['file'])) . "</span><br />\n";

	$hap_strs = array();
	$dep_strs = array();

	foreach ($match as $m) {
	    $a = 
	      "<a target=\"blank\" href=\"$root_path/tag.php?db=$database&batch_id=$batch_id&sample_id=$m[id]&tag_id=$m[tag_id]\" " . 
	      "title=\"#$m[tag_id]\" style=\"color: " . $alleles[$m['allele']] . ";\">$m[allele]</a>";
	    array_push($hap_strs, $a);
	    $a = 
	      "<a target=\"blank\" href=\"$root_path/tag.php?db=$database&batch_id=$batch_id&sample_id=$m[id]&tag_id=$m[tag_id]\" " . 
	      "title=\"#$m[tag_id]\" style=\"color: " . $alleles[$m['allele']] . ";\">$m[depth]</a>";
	    array_push($dep_strs, $a);
	}

	$hap_str = implode(" / ", $hap_strs);
	$dep_str = implode(" / ", $dep_strs);

	if (count($gtypes) > 0 && isset($gtypes[$sample])) {

	  $id      = "gtype_" . $batch_id . "_" . $tag_id . "_" . $gtypes[$sample]['id'];
	  $url     = "$root_path/correct_genotype.php?db=$database&batch_id=$batch_id&tag_id=$tag_id&sample_id=" . $gtypes[$sample]['id'];
	  $jscript = "correct_genotype('$id', '$url')";
	  $blur_js = "cancel_correction('$id')";

	  if (strlen($gtypes[$sample]['corrected']) > 0) {
  	      $sel = generate_element_select($id, $marker_types[$gtypes[$sample]['marker']], strtolower($gtypes[$sample]['corrected']), $jscript, $blur_js);
	      $genotype = "<span class=\"corrected\">" . $gtypes[$sample]['corrected'] . "</span>";
	  } else {
	      $sel = generate_element_select($id, $marker_types[$gtypes[$sample]['marker']], strtolower($gtypes[$sample]['genotype']), $jscript, $blur_js);
	      $genotype = $gtypes[$sample]['genotype'];
	  }

	  $gen_str = 
	    "<div id=\"gen_{$i}\" style=\"display: none;\">" .
	    "<div id=\"{$id}_div\"><a onclick=\"toggle_correction('$id')\">" . $genotype . "</a></div>\n" . 
	    "  <div id=\"{$id}_sel\" style=\"display: none;\">\n" . 
	    $sel . 
	    "  </div>";
	} else {
	  $gen_str = "";
	}

	print
	  "<div class=\"haplotype\" id=\"hap_{$i}\">$hap_str</div><div id=\"dep_{$i}\" style=\"display: none;\">$dep_str</div>$gen_str" .
	  "</td>\n";

	if ($i % $num_cols == 0)
	  print 
            "</tr>\n" .
            "<tr>\n";
    }

    while ($i % $num_cols != 0) {
      print "  <td></td>\n";
      $i++;
    }

    print "</tr>\n";
}

echo <<< EOQ
</table>
</td>
</tr>
</table>

EOQ;

write_compact_footer();

function print_sumstats($db, $pop_names) {
    global $batch_id, $tag_id;

    $result = $db['stats_sth']->execute(array($batch_id, $tag_id));
    check_db_error($result, __FILE__, __LINE__);

    $stats = array();
    $pops  = array();

    while ($row = $result->fetchRow()) {
        $a = array('col'     => $row['col'],
		   'n'       => $row['n'],
		   'p'       => $row['p'],
		   'obs_het' => $row['obs_het'],
		   'obs_hom' => $row['obs_hom'],
		   'exp_het' => $row['exp_het'],
		   'exp_hom' => $row['exp_hom'],
		   'pi'      => $row['pi'],
		   'fis'     => $row['fis'],
		   'pop_id'  => $row['pop_id']);

	if (!isset($stats[$row['col']]))
  	    $stats[$row['col']] = array();

	$stats[$row['col']][$row['pop_id']] = $a;

	$pops[$row['pop_id']] = $row['pop_id'];
    }

    ksort($stats);
    ksort($pops);

    $index = count($pop_names) == 0 ? $pops : $pop_names;

    $n = 0;
    foreach ($stats as $col => $stat) {    
      $style = $n == 0 ? "" : "style=\"display: none;\"";

      print
	"    <table id=\"{$col}_sumstats\" class=\"sumstats\" $style>\n" .
	  "    <tr>\n" .
	  "      <td class=\"key\" style=\"background-color: #fff; border-left: 0px;\">&nbsp;</td>\n" .
	  "      <td class=\"key\">P</td>\n" .
	  "      <td class=\"key\">N</td>\n" .
	  "      <td class=\"key\">Obs Het</td>\n" .
	  "      <td class=\"key\">Obs Hom</td>\n" .
	  "      <td class=\"key\">Exp Het</td>\n" .
	  "      <td class=\"key\">Exp Hom</td>\n" .
	  "      <td class=\"key\">&pi;</td>\n" .
	  "      <td class=\"key\">F<sub>is</sub></td>\n" .
	  "    </tr>\n";

      foreach ($index as $pop_id => $pop_name) {
	  $s    = $stat[$pop_id];
	  $p    = $s['p']       < 1 ? sprintf("%.5f", $s['p']) : $s['p'];
	  $ohet = $s['obs_het'] > 0 ? sprintf("%.3f", $s['obs_het']) : $s['obs_het'];
	  $ohom = $s['obs_hom'] < 1 ? sprintf("%.3f", $s['obs_hom']) : $s['obs_hom'];
	  $ehet = $s['exp_het'] > 0 ? sprintf("%.3f", $s['exp_het']) : $s['exp_het'];
	  $ehom = $s['exp_hom'] < 1 ? sprintf("%.3f", $s['exp_hom']) : $s['exp_hom'];
	  $pi   = $s['pi']      > 0 ? sprintf("%.3f", $s['pi']) : $s['pi'];
	  $fis  = $s['fis']    != 0 ? sprintf("%.3f", $s['fis']) : "0";
	  print
	    "    <tr>\n" .
	    "      <td class=\"pop_id\">$pop_name</td>\n" .
	    "      <td>$p</td>\n" .
	    "      <td>" . $s['n'] . "</td>\n" .
	    "      <td>$ohet</td>\n" .
	    "      <td>$ohom</td>\n" .
	    "      <td>$ehet</td>\n" .
	    "      <td>$ehom</td>\n" .
	    "      <td>$pi</td>\n" .
	    "      <td>$fis</td>\n" .
	    "    </tr>\n";
      }
      print
	"    </table>\n";
      $n++;
    }
}

function print_fst($db, $pop_names) {
    global $batch_id, $tag_id;

    $result = $db['fst_sth']->execute(array($batch_id, $tag_id));
    check_db_error($result, __FILE__, __LINE__);

    $stats = array();
    $pops  = array();

    while ($row = $result->fetchRow()) {
        $a = array('col'   => $row['col'],
		   'pid_1' => $row['pop_id_1'],
		   'pid_2' => $row['pop_id_2'],
		   'pi_o'  => $row['pi_o'],
		   'fst'   => $row['fst']);

	if (!isset($stats[$row['col']]))
  	    $stats[$row['col']] = array();

	array_push($stats[$row['col']], $a);

	$pops[$row['pop_id_1']] = $row['pop_id_1'];
	$pops[$row['pop_id_2']] = $row['pop_id_2'];
    }

    ksort($stats);
    ksort($pops);

    print
	"    <table class=\"fst_key\">\n" .
	"    <tr>\n" .
	"      <td class=\"fst\">F<sub>st</sub></td>\n" .
	"      <td class=\"pi\">&pi;<sub>overall</sub></td>\n" .
	"    </tr>\n" .
	"    </table>\n";

    $index = count($pop_names) == 0 ? $pops : $pop_names;

    $n = 0;
    foreach ($stats as $col => $stat) {
      $style = $n == 0 ? "" : "style=\"display: none;\"";

      $matrix = array();
      foreach ($index as $pop_id => $pop_name)
	$matrix[$pop_id] = array();

      foreach ($stat as $s) {
	  $matrix[$s['pid_1']][$s['pid_2']] = $s['pi_o'];
	  $matrix[$s['pid_2']][$s['pid_1']] = $s['fst'];
      }

      $keys    = array_keys($pops);
      $pop_cnt = count($pops);

      print
	"    <table id=\"{$col}_fst\" class=\"fst\" $style>\n" .
	"    <tr>\n" .
	"      <td class=\"key\" style=\"background-color: #fff; border-left: 0px; border-right: 0px;\">&nbsp;</td>\n";

      for ($i = 0; $i < $pop_cnt; $i++)
	print "      <td class=\"key\">" . $index[$keys[$i]] . "</td>\n";

      print "    </tr>\n";

      for ($i = 0; $i < $pop_cnt; $i++) {
	print 
	  "    <tr>\n" .
	  "      <td class=\"pop_id\">" . $index[$keys[$i]] . "</td>\n";

	for ($j = 0; $j < $pop_cnt; $j++) {
	  $class = ($i < $j) ? "class=\"pi\"" : "class=\"fst\"";
	  
	  if ($i == $j)
	    print "      <td class=\"diagonal\">&nbsp;</td>\n";
	  else if (!isset($matrix[$keys[$i]][$keys[$j]]))
	    print "      <td $class>&nbsp;</td>\n";
	  else
	    print 
	      "      <td id=\"{$col}_{$i}_{$j}\" $class onmouseover=\"highlight_cells($col, $i, $j)\" onmouseout=\"unhighlight_cells($col, $i, $j)\">" . 
	      sprintf("%0.5f", $matrix[$keys[$i]][$keys[$j]]) . 
	      "</td>\n";
	}

	print "    </tr>\n";
      }

      print "    </table>\n";
      $n++;
    }
}

function print_population_name($db, $batch_id, $pop_id, $pop_names) {
    global $root_path;

    $pop_str = isset($pop_names[$pop_id]) ? $pop_names[$pop_id] : "Population $pop_id";

    echo <<< EOQ
<div class="pop_annotation" id="{$pop_id}_pop" onclick="toggle_population($pop_id)">$pop_str</div>
<div id="{$pop_id}_div" style="display: none; font-size: small;">
  <form id="{$pop_id}_frm">
    <input type="hidden" name="url" value="$root_path/annotate_population.php?db=$db&batch_id=$batch_id&pop_id=$pop_id" />
    <input type="input" size=15 name="pop_name" value="" />
    <a onclick="annotate_population('$pop_id')">save</a>|<a onclick="toggle_population('$pop_id')">cancel</a>
  </form>
</div>

EOQ;
}

?>
