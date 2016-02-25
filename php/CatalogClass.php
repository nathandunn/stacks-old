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
require_once("Locus.php");

class Catalog {
    var $db;       // Hash of database statement handles.
    var $display;  // An array of form variables necessary to write URLs.
    var $params;   // An array of parameters to fetch our results from the database.
    var $queries;  // An array of SQL queries.
    var $loci;     // Array of groups objects.
    var $num_loci; // Total number of groups, after filtering

    function Catalog($db, $batch_id, $display_params) {
	$this->db       = $db;
	$this->batch    = $batch_id;
	$this->display  = $display_params;
	$this->queries  = array();
	$this->params   = array();
	$this->loci     = array();
	$this->num_loci = 0;

	$this->prepare_queries();
    }

    function prepare_queries() {

        $query = 
            "SELECT allele FROM catalog_alleles " . 
            "WHERE batch_id=? AND tag_id=?";
        if (!($this->db['allele_sth'] = $this->db['dbh']->prepare($query)))
            write_db_error($this->db['dbh'], __FILE__, __LINE__);

        $query = 
            "SELECT col, rank_1, rank_2 FROM catalog_snps " . 
            "WHERE batch_id=? AND tag_id=? ORDER BY col";
        if (!($this->db['snp_sth'] = $this->db['dbh']->prepare($query)))
            write_db_error($this->db['dbh'], __FILE__, __LINE__);

        $query = 
            "SELECT samples.id, samples.sample_id, samples.type, file, tag_id, allele " . 
            "FROM matches " . 
            "JOIN samples ON (matches.sample_id=samples.id) " . 
            "WHERE matches.batch_id=? AND catalog_id=? ORDER BY samples.id";
        if (!($this->db['mat_sth'] = $this->db['dbh']->prepare($query)))
            write_db_error($this->db['dbh'], __FILE__, __LINE__);

        $query = 
            "SELECT COUNT(tag_id) as count FROM catalog_index " . 
            "WHERE batch_id=?";
        $query .= $this->apply_query_filters(); 
        $this->queries['tag_count'] = $query;

        $query = 
            "SELECT catalog_index.tag_id as tag_id, chr, bp, snps, alleles, parents, progeny, valid_progeny, " . 
            "seq, marker, max_pct, ratio, ests, pe_radtags, blast_hits, external_id " .
            "FROM catalog_index " .
            "JOIN catalog_tags ON (catalog_index.cat_id=catalog_tags.id) " . 
            "LEFT JOIN catalog_annotations ON " . 
            "(" . 
            "catalog_index.batch_id=catalog_annotations.batch_id AND " . 
            "catalog_index.tag_id=catalog_annotations.catalog_id" . 
            ") " .
            "WHERE catalog_index.batch_id=?";
        $query .= $this->apply_query_filters();

        $this->queries['tag'] = $query;
    }

    function prepare_filter_parameters() {
        $this->params[] = &$this->batch;
	$typestr = "i";

	$filters = $this->display['filter_type'];

	if (!isset($filters)) {
	    array_unshift($param, $typestr);
	    return;
	}

	foreach ($filters as $filter) {

            if ($filter == "snps") {
                $this->params[] = &$this->display['filter_snps_l'];
                $this->params[] = &$this->display['filter_snps_u'];
		$typestr .= "ii";

            } else if ($filter == "alle") {
                $this->params[] = &$this->display['filter_alle_l'];
                $this->params[] = &$this->display['filter_alle_u'];
		$typestr .= "ii";

            } else if ($filter == "pare") {
                $this->params[] = &$this->display['filter_pare_l'];
                $this->params[] = &$this->display['filter_pare_u'];
		$typestr .= "ii";

            } else if ($filter == "prog") {
                $this->params[] = &$this->display['filter_prog'];
		$typestr .= "i";

            } else if ($filter == "vprog") {
                $this->params[] = &$this->display['filter_vprog'];
		$typestr .= "i";

            } else if ($filter == "cata") {
                $this->params[] = &$this->display['filter_cata'];
		$typestr .= "i";

            } else if ($filter == "gcnt") {
	        $this->params[] = &$this->display['filter_gcnt'];
		$typestr .= "i";

            } else if ($filter == "ref") {
		$this->params[] = &$this->display['filter_ref'];
		$typestr .= "i";

	    } else if ($filter == "loc") {
		$this->params[] = &$this->display['filter_chr'];
		$this->params[] = &$this->display['filter_sbp'];
		$this->params[] = &$this->display['filter_ebp'];
		$typestr .= "sii";
	
	    } else if ($filter == "mark") {
		$this->params[] = &$this->display['filter_mark'];
            }
	}

	array_unshift($this->params, $typestr);
    }

    function apply_query_filters() {
        $sql_filters =
            array("cata"  => "(catalog_index.tag_id = ?)", 
		  "alle"  => "(alleles >= ? AND alleles <= ?)", 
		  "snps"  => "(snps >= ? AND snps <= ?)",
		  "pare"  => "(parents >= ? AND parents <= ?)",
                  "prog"  => "(progeny >= ?)",
                  "vprog" => "(valid_progeny >= ?)",
                  "mark"  => "(marker LIKE ?)", 
                  "est"   => "(ests > 0)",
                  "pe"    => "(pe_radtags > 0)",
                  "blast" => "(blast_hits > 0)",
		  "gcnt"  => "(geno_cnt >= ?)",
		  "ref"   => "(catalog_index.type = ?)",
		  "loc"   => "(catalog_index.chr = ? && catalog_index.bp >= ? && catalog_index.bp <= ?)");

        $filters = $this->display['filter_type'];

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

    function &loci() {
	return $this->loci;
    }

    function &locus($id) {
	return $this->loci[$id];
    }

    function num_loci() {
	return $this->num_loci;
    }

    function determine_count() {

	if (!($this->db['tag_count_sth'] = $this->db['dbh']->prepare($this->queries['tag_count'])))
	    write_db_error($this->db['dbh'], __FILE__, __LINE__);

	$this->prepare_filter_parameters();

	array_unshift($this->params, $this->db['tag_count_sth']);
	call_user_func_array("mysqli_stmt_bind_param", $this->params);
	array_shift($this->params);
	if (!$this->db['tag_count_sth']->execute())
	    write_db_error($this->db['tag_count_sth'], __FILE__, __LINE__);
	$res = $this->db['tag_count_sth']->get_result();

	$row = $res->fetch_assoc();
	$this->num_loci = $row['count'];
    }

    function populate($start_group, $num_groups) {
	//
	// We only want to load genes between $start_gene and $end_gene. 
	//
	$this->queries['tag'] .= " LIMIT " . $start_group . ", " . $num_groups;
	
	if (!($this->db['tag_sth'] = $this->db['dbh']->prepare($this->queries['tag'])))
	    write_db_error($this->db['tag_sth'], __FILE__, __LINE__);

	$this->prepare_filter_parameters();

	//
	// Fetch the results and populate the array of groups.
	//
	array_unshift($this->params, $this->db['tag_sth']);
	call_user_func_array("mysqli_stmt_bind_param", $this->params);
	array_shift($this->params);
	if (!$this->db['tag_sth']->execute())
	    write_db_error($this->db['tag_sth'], __FILE__, __LINE__);
	$res = $this->db['tag_sth']->get_result();

	while ($row = $res->fetch_assoc()) {
            $locus = new Locus();
            $locus->id         = $row['tag_id'];
            $locus->annotation = $row['external_id'];
            $locus->chr        = $row['chr'];
            $locus->bp         = $row['bp'];
            $locus->marker     = $row['marker'];
            $locus->seq        = $row['seq'];

            $locus->num_alleles   = $row['alleles'];
            $locus->num_snps      = $row['snps'];
            $locus->num_parents   = $row['parents'];
            $locus->num_progeny   = $row['progeny'];
            $locus->valid_progeny = $row['valid_progeny'];
            $locus->num_ests      = $row['ests'];
            $locus->num_pe_tags   = $row['pe_radtags'];
            $locus->num_blast     = $row['blast_hits'];

            //
            // Fetch SNPs and Alleles
            //
	    if (!$this->db['snp_sth']->bind_param("ii", $this->batch, $locus->id))
		write_db_error($db['snp_sth'], __FILE__, __LINE__);
	    if (!$this->db['snp_sth']->execute())
		write_db_error($db['snp_sth'], __FILE__, __LINE__);
	    $snp_res = $this->db['snp_sth']->get_result();

            while ($snp_row = $snp_res->fetch_assoc()) {
                $locus->snps .= $snp_row['col'] . "," . $snp_row['rank_1'] . ">" . $snp_row['rank_2'] . ";";
            }
            $locus->snps = substr($locus->snps, 0, -1);

	    if (!$this->db['all_sth']->bind_param("ii", $this->batch, $locus->id))
		write_db_error($db['all_sth'], __FILE__, __LINE__);
	    if (!$this->db['all_sth']->execute())
		write_db_error($db['all_sth'], __FILE__, __LINE__);
	    $all_res = $this->db['all_sth']->get_result();

            while ($all_row = $all_res->fetch_assoc()) {
                $locus->alleles .= $all_row['allele'] . ";";
            }
            $locus->alleles = substr($locus->alleles, 0, -1);

            //
            // Add genotypes
            //
	    if (!$this->db['mat_sth']->bind_param("ii", $this->batch, $locus->id))
		write_db_error($db['mat_sth'], __FILE__, __LINE__);
	    if (!$this->db['mat_sth']->execute())
		write_db_error($db['mat_sth'], __FILE__, __LINE__);
	    $gen_res = $this->db['mat_sth']->get_result();

            while ($gen_row = $gen_res->fetch_assoc())
                $locus->add_genotype($gen_row['id'], $gen_row['file'], $gen_row['allele']);

	    $this->loci[$row['tag_id']] = $locus;
	}
    }
}
