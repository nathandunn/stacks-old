<?php

function db_connect($database) {
    $dsn = array(
		 'phptype'  => 'mysql',
		 'username' => 'acos_user',
		 'password' => 'orthologs',
		 'hostspec' => 'localhost',
		 'database' => $database,
		 );

    $dbh =& MDB2::connect($dsn, $options);

    if (PEAR::isError($dbh)) {
	die("File: " . __FILE__ . " (line " . __LINE__ . ") " . $dbh->getMessage());
    }

    // Set the database package to always return 
    // results as an associative array
    $dbh->setFetchMode(MDB2_FETCHMODE_ASSOC);

    // The $db array will hold the database handle and 
    // common, prepared SQL statements.
    $db = array();
    $db['dbh']  = $dbh;
    $db['name'] = $database;

    return $db;
}

function check_db_error($sth, $file, $line) {

    if (PEAR::isError($sth)) {

	$error_str = 
	    "File: $file (line $line)<br />\n " .
	    "<strong>" . $sth->getMessage() . "</strong><br />\n" .
	    $sth->getUserInfo() . "<br />\n";

	die($error_str);
    }
}

?>
