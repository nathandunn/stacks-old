<?php
require_once("header.php");
require_once("radtag_functions.php");

$database = isset($_GET['db']) ? $_GET['db'] : "radtags";
$seq_id   = isset($_GET['id']) ? $_GET['id'] : 0;

// Connect to the database
$db = db_connect($database);

// Save these variables for automatic URL formation later on.
$display = array();
$display['db']     = $database;
$display['seq_id'] = $seq_id;

//
// Prepare some SQL queries
//
$query = 
    "SELECT id, catalog_id, seq_id, type, seq FROM sequence " . 
    "WHERE id=?";
$db['seq_sth'] = $db['dbh']->prepare($query);
check_db_error($db['seq_sth'], __FILE__, __LINE__);


$page_title = "Catalog RAD-Tag Sequence Viewer";
write_compact_header($page_title, $batch);

$result = $db['seq_sth']->execute($seq_id);
check_db_error($result, __FILE__, __LINE__);

$row = $result->fetchRow();

//
// Line wrap the sequence to 80 characters
//
$full_seq = $row['seq'];
$seq = "";

do {
    $seq     .= substr($full_seq, 0, 80) . "\n";
    $full_seq = substr($full_seq, 80);
} while (strlen($full_seq) > 80);
        
if (strlen($seq) > 0) {
    $seq .= $full_seq . "\n";
}

    echo <<< EOQ
<pre>
>$row[catalog_id]|$row[seq_id]|$row[type]
$seq
</pre>

EOQ;

write_compact_footer();

?>
