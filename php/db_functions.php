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

function db_connect($database) {
    global $db_user, $db_pass, $db_host;

    $dsn = array(
                 'phptype'  => 'mysql',
                 'username' => $db_user,
                 'password' => $db_pass,
                 'hostspec' => $db_host,
		 'port'     => 3306
                 );
    $options = array();

    if (strlen($database) > 0)
        $dsn['database'] = $database;
    else
        $dsn['database'] = false;

    $dbh = MDB2::connect($dsn, $options);

    if (MDB2::isError($dbh)) {
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

    if (MDB2::isError($sth)) {

	$error_str = 
	    "File: $file (line $line)<br />\n " .
	    "<strong>" . $sth->getMessage() . "</strong><br />\n" .
	    $sth->getUserInfo() . "<br />\n";

	die($error_str);
    }
}

?>
