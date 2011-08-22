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

if (isset($_GET['db']))
    $database  = $_GET['db'];
else if (isset($_POST['db']))
    $database = $_POST['db'];

if (isset($_GET['tag_id']))
    $tag_id  = $_GET['tag_id'];
else if (isset($_POST['tag_id']))
    $tag_id = $_POST['tag_id'];
else 
    $tag_id = 0;

if (isset($_GET['batch_id']))
    $batch_id  = $_GET['batch_id'];
else if (isset($_POST['batch_id']))
    $batch_id = $_POST['batch_id'];
else 
    $batch_id = 0;

// Connect to the database
$db = db_connect($database);

// Save these variables for automatic URL formation later on.
$display = array();
$display['db']       = $database;
$display['tag_id']   = $tag_id;
$display['batch_id'] = $batch_id;

//
// Prepare the possible select lists we will want to construct
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
    "SELECT count(samples.id) as count FROM samples WHERE batch_id=?";
$db['samp_sth'] = $db['dbh']->prepare($query);
check_db_error($db['samp_sth'], __FILE__, __LINE__);

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

$page_title = "Catalog Genotype Viewer";
write_compact_header($page_title);

// 
// Get number of samples so we can determine how many rows to display 
// in the genotype table.
//
$result = $db['samp_sth']->execute($batch_id);
check_db_error($result, __FILE__, __LINE__);
$row = $result->fetchRow();
$num_samples = $row['count'];
$num_cols    = 10;
$num_rows    = ceil($num_samples / $num_cols);
$gtypes      = array();

$result = $db['geno_sth']->execute(array($batch_id, $tag_id));
check_db_error($result, __FILE__, __LINE__);

if ($result->numRows() == 0) {
    print 
        "<h4 style=\"margin-left: auto; margin-right: auto; text-align: center;\">" . 
        "This marker has no genotypes, probably because this tag does not have enough mappable progeny.</h4>\n";
    write_compact_footer();
    return;
}

while ($row = $result->fetchRow()) {
    $gtypes[$row['sample_id']] = array('file'      => $row['file'], 
                                       'genotype'  => $row['genotype'], 
                                       'corrected' => $row['corrected'],
				       'marker'    => $row['marker']);
}

print 
    "<form id=\"genotypes\" name=\"genotypes\" method=\"post\" action=\"$root_path/correct_genotypes.php\">\n" .
    "<input type=\"hidden\" name=\"op\" id=\"op\" value=\"\" />\n" .
    "<input type=\"hidden\" name=\"batch_id\" value=\"$batch_id\" />\n" .
    "<input type=\"hidden\" name=\"tag_id\" value=\"$tag_id\" />\n" .
    "<input type=\"hidden\" name=\"db\" value=\"$database\" />\n" .
    "<table class=\"genotypes\">\n" .
    "<tr>\n";
$i = 0;
foreach ($gtypes as $sample_id => $sample) {
    $i++;

    $id  = "gtype_" . $batch_id . "_" . $tag_id . "_" . $sample_id;

    if (strlen($sample['corrected']) > 0) {
        $sel = generate_element_select($id, $marker_types[$sample['marker']], strtolower($sample['corrected']), "");
        $genotype = "<span class=\"corrected\">$sample[corrected]</span>";

    } else {
        $sel = generate_element_select($id, $marker_types[$sample['marker']], strtolower($sample['genotype']), "");
        $genotype = $sample['genotype'];
    }

    print 
        "  <td>" .
        "<span class=\"title\">" . ucfirst(str_replace("_", " ", $sample['file'])) . "</span><br />" .
        "<a onclick=\"toggle_sel('{$id}_div')\">" . $genotype . "</a>\n" . 
        "  <div id=\"{$id}_div\" style=\"display: none;\">\n" . 
        $sel . 
        "  </div>" .
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
<tr>
  <td colspan="$num_cols" style="text-align: right;">
    <input type="button" value="Reset" onclick="set_operation('genotypes', 'reset')" />
    <input type="button" value="Submit" onclick="set_operation('genotypes', 'correct')" />
  </td>
</tr>
</table>
</form>

EOQ;

write_compact_footer();

?>
