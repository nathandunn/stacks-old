<?php
//
// Copyright 2010-2016, Julian Catchen <jcatchen@illinois.edu>
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

    $dbh = new mysqli($db_host, $db_user, $db_pass, $database);

    if ($dbh->connect_errno) {
        die("File: " . __FILE__ . " (line " . __LINE__ . ") Failed to connect to MySQL: (" . $dbh->connect_errno . ") " . $dbh->connect_error);
    }

    //
    // The $db array will hold the database handle and 
    // common, prepared SQL statements.
    //
    $db = array();
    $db['dbh']  = $dbh;
    $db['name'] = $database;

    return $db;
}

function write_db_error($sth, $file, $line) {
    $error_str = 
        "File: $file (line $line)<br />\n " .
        $sth->errno . ": <strong>" . $sth->error . "</strong><br />\n";
    die($error_str);
}

?>
