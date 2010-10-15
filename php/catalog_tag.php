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
    "SELECT samples.id, samples.sample_id, samples.type, file, tag_id, allele " . 
    "FROM matches " . 
    "JOIN samples ON (matches.sample_id=samples.id) " . 
    "WHERE matches.batch_id=? AND catalog_id=? ORDER BY samples.id";
$db['mat_sth'] = $db['dbh']->prepare($query);
check_db_error($db['mat_sth'], __FILE__, __LINE__);

$query = 
    "SELECT col, rank_1, rank_2 FROM catalog_snps " . 
    "WHERE batch_id=? AND tag_id=? ORDER BY col";
$db['snp_sth'] = $db['dbh']->prepare($query);
check_db_error($db['snp_sth'], __FILE__, __LINE__);

$query = 
    "SELECT allele FROM catalog_alleles " . 
    "WHERE batch_id=? AND tag_id=? ";
$db['all_sth'] = $db['dbh']->prepare($query);
check_db_error($db['all_sth'], __FILE__, __LINE__);

$page_title = "Catalog RAD-Tag Viewer";
write_compact_header($page_title, $batch);

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
    print "    Column: $row[col]; $row[rank_1]/$row[rank_2]<br />\n";
}
print "  </td>\n";

$result = $db['all_sth']->execute(array($batch_id, $tag_id));
check_db_error($result, __FILE__, __LINE__);

$alleles = array();
$i       = 0;

print "  <td>\n";
while ($row = $result->fetchRow()) {
    $alleles[$row['allele']] = $colors[$i % $color_size];

    print "  <span style=\"color: " . $alleles[$row['allele']] . ";\">$row[allele]</span><br />\n";

    $i++;
}
print "  </td>\n";

$result = $db['mat_sth']->execute(array($batch_id, $tag_id));
check_db_error($result, __FILE__, __LINE__);

while ($row = $result->fetchRow()) {
    $a = array('id'     => $row['id'],
               'file'   => $row['file'], 
               'allele' => $row['allele'], 
               'tag_id' => $row['tag_id']);
    if (!isset($gtypes[$row['file']]))
        $gtypes[$row['file']] = array();

    array_push($gtypes[$row['file']], $a);
}

print 
    "<td style=\"text-align: right;\">\n" .
    "<table class=\"genotypes\" style=\"padding: 1em; margin-top: 1em; margin-bottom: 1em;\">\n" .
    "<tr>\n";

$i        = 0;
$num_cols = count($gtypes) < 10 ? count($gtypes) : 10;

foreach ($gtypes as $sample => $match) {
    $i++;

    print 
        "  <td>" .
        "<span class=\"title\">" . ucfirst(str_replace("_", " ", $match[0]['file'])) . "</span><br />\n";

    $strs = array();

    foreach ($match as $m) {
        $a = 
            "<a target=\"blank\" href=\"$root_path/tag.php?db=$database&batch_id=$batch_id&sample_id=$m[id]&tag_id=$m[tag_id]\" " . 
            "title=\"#$m[tag_id]\" style=\"color: " . $alleles[$m['allele']] . ";\">$m[allele]</a>";
        array_push($strs, $a);
    }

    $str = implode(" / ", $strs);
    print
        $str .
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
