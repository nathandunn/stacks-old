<?php
require_once("header.php");

$database = isset($_GET['db']) ? $_GET['db'] : "";

if (strlen($database) == 0)
    write_database_list($database);
else
    write_database($database);

function write_database_list($database) {
    global $db_user, $db_pass, $root_path;

    $databases = array();
    
    exec("mysql --user=$db_user --password=$db_pass -N -B -e \"SHOW DATABASES LIKE '%_radtags'\"", $databases);

    $page_title = "RAD-Tag Databases";
    write_header($page_title);

    echo <<< EOQ
<h4 class="info_head">
  <img id="sources_img" src="/acos/images/caret-d.png" />
  <a onclick="toggle_div('sources', '/acos/images', 'page_state');">RAD-Tag Databases</a>
</h4>

<a name="results_top"></a>
<table class="db" style="width: 35%;">
<tr>
  <th style="width: 100%;">Database</th>
</tr>

EOQ;

    foreach ($databases as $dbase) {
        print
            "<tr>\n" .
            "  <td><a href=\"$root_path/index.php?db=" . $dbase . "\">$dbase</a></td>\n" .
            "</tr>\n";
    }

    echo <<< EOQ
</table>

EOQ;

    write_footer();
}

function write_database($database) {
    global $root_path;

    //
    // Connect to the database
    //
    $db = db_connect($database);

    //
    // Prepare some SQL queries
    //
    $query = 
        "SELECT id, date, description FROM batches";
    $db['batch_sth'] = $db['dbh']->prepare($query);
    check_db_error($db['batch_sth'], __FILE__, __LINE__);

    $page_title = "RAD-Tag Analyses";
    write_header($page_title);

    echo <<< EOQ
<h4 class="info_head">
  <img id="sources_img" src="/acos/images/caret-d.png" />
  <a onclick="toggle_div('sources', '/acos/images', 'page_state');">RAD-Tag Samples</a>
</h4>

<a name="results_top"></a>
<table class="db" style="width: 75%;">
<tr>
  <th style="width: 10%;">&nbsp;</th>
  <th style="width: 10%;">&nbsp;</th>
  <th style="width: 10%;">Batch ID</th>
  <th style="width: 20%;">Date</th>
  <th style="width: 50%;">Description</th>
</tr>

EOQ;

    $result = $db['batch_sth']->execute();
    check_db_error($result, __FILE__, __LINE__);

    while ($row = $result->fetchRow()) {
        print
            "<tr>\n" .
            "  <td class=\"s\"><a href=\"$root_path/catalog.php?db=$database&id=$row[id]\">Catalog</a></td>\n" .
            "  <td class=\"s\"><a href=\"$root_path/samples.php?db=$database&id=$row[id]\">Samples</a></td>\n" .
            "  <td>" . $row['id'] . "</td>\n" .
            "  <td>" . $row['date'] . "</td>\n" .
            "  <td>" . $row['description'] . "</td>\n" .
            "</tr>\n";
    }

    echo <<< EOQ
</table>

EOQ;

    write_footer();
}
 
?>
