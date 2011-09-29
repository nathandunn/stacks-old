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
// Prepare the possible select lists we will want to construct to manualy correct genotypes.
//
$marker_types = array('ab/--' => array('aa', 'bb', '-'),
		      '--/ab' => array('aa', 'bb', '-'),
		      'ab/aa' => array('aa', 'ab', '-'),
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
    "SELECT samples.id, samples.sample_id, samples.type, file, tag_id, allele, depth " . 
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

echo <<< EOQ
<table class="catalog">
<tr>
  <th style="width: 15%;">SNPs</th>
  <th style="width: 15%;">Alleles</th>
  <th style="width: 70%;">Matching Samples</th>
</tr>
<tr>

EOQ;

$result = $db['snp_sth']->execute(array($batch_id, $tag_id));
check_db_error($result, __FILE__, __LINE__);

print "  <td>\n";
while ($row = $result->fetchRow()) {
    $snp = "$row[rank_1]/$row[rank_2]";
    if (strlen($row['rank_3']) > 0) $snp .= "/$row[rank_3]";
    if (strlen($row['rank_4']) > 0) $snp .= "/$row[rank_4]";

    print "    Column: $row[col]; $snp<br />\n";
}
print "  </td>\n" .
      "  <td>\n";

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
	"  <span style=\"color: " . $alleles[$hapl] . ";\">" . 
	"<acronym title=\"Genotype\">$geno</acronym> : " . 
	"<acronym title=\"Observed Haplotype\">$hapl</acronym></span><br />\n";
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

print "  </td>\n";

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
	       'depth'  => $row['depth']);
    if (!isset($htypes[$row['file']]))
        $htypes[$row['file']] = array();

    array_push($htypes[$row['file']], $a);
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
  $gtype_str = "<input type=\"checkbox\" onclick=\"toggle_genotypes($tag_id, 'locus_gtypes', 'gen')\" /> Genotypes";
else
  $gtype_str = "";

echo <<< EOQ
<td style="text-align: right;">

<table id="locus_gtypes" class="genotypes">
<tr>
  <td colspan="10" class="gtype_toggle">
    <strong>View:</strong>
    <input type="checkbox" checked="checked" onclick="toggle_genotypes($tag_id, 'locus_gtypes', 'hap')" />Haplotypes
    <input type="checkbox" onclick="toggle_genotypes($tag_id, 'locus_gtypes', 'dep')" /> Allele Depths
    $gtype_str
  </td>
</tr>
<tr>

EOQ;

$i        = 0;
$num_cols = count($htypes) < 10 ? count($htypes) : 10;

foreach ($htypes as $sample => $match) {
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

	if (strlen($gtypes[$sample]['corrected']) > 0) {
	    $sel = generate_element_select($id, $marker_types[$gtypes[$sample]['marker']], strtolower($gtypes[$sample]['corrected']), $jscript);
	    $genotype = "<span class=\"corrected\">" . $gtypes[$sample]['corrected'] . "</span>";
	} else {
	    $sel = generate_element_select($id, $marker_types[$gtypes[$sample]['marker']], strtolower($gtypes[$sample]['genotype']), $jscript);
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
        "<div id=\"hap_{$i}\">$hap_str</div><div id=\"dep_{$i}\" style=\"display: none;\">$dep_str</div>$gen_str" .
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

echo <<< EOQ
</tr>
</table>
</td>
</tr>
</table>

EOQ;

write_compact_footer();

?>
