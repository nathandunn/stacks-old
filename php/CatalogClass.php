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
        $this->db['allele_sth'] = $this->db['dbh']->prepare($query);
        check_db_error($this->db['allele_sth'], __FILE__, __LINE__);

        $query = 
            "SELECT col, rank_1, rank_2 FROM catalog_snps " . 
            "WHERE batch_id=? AND tag_id=? ORDER BY col";
        $this->db['snp_sth'] = $this->db['dbh']->prepare($query);
        check_db_error($this->db['snp_sth'], __FILE__, __LINE__);

        $query = 
            "SELECT samples.id, samples.sample_id, samples.type, file, tag_id, allele " . 
            "FROM matches " . 
            "JOIN samples ON (matches.sample_id=samples.id) " . 
            "WHERE matches.batch_id=? AND catalog_id=? ORDER BY samples.id";
        $this->db['mat_sth'] = $this->db['dbh']->prepare($query);
        check_db_error($this->db['mat_sth'], __FILE__, __LINE__);

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
        array_push($this->params, $this->batch);

	$filters = $this->display['filter_type'];

	if (!isset($filters))
	    return;

	foreach ($filters as $filter) {

            if ($filter == "snps") {
                array_push($this->params, $this->display['filter_snps_l']);
                array_push($this->params, $this->display['filter_snps_u']);

            } else if ($filter == "alle") {
                array_push($this->params, $this->display['filter_alle_l']);
                array_push($this->params, $this->display['filter_alle_u']);

            } else if ($filter == "pare") {
                array_push($this->params, $this->display['filter_pare_l']);
                array_push($this->params, $this->display['filter_pare_u']);
	
            } else if ($filter == "prog") {
                array_push($this->params, $this->display['filter_prog']);

            } else if ($filter == "vprog") {
                array_push($this->params, $this->display['filter_vprog']);

            } else if ($filter == "cata") {
                array_push($this->params, $this->display['filter_cata']);

            } else if ($filter == "gcnt") {
	        array_push($this->params, $this->display['filter_gcnt']);
	
	    } else if ($filter == "est") {
                array_push($this->params, 0);

            } else if ($filter == "pe") {
                array_push($this->params, 0);

            } else if ($filter == "blast") {
                array_push($this->params, 0);

            } else if ($filter == "ref") {
	      array_push($this->params, $this->display['filter_ref']);
	
	    } else if ($filter == "loc") {
	      array_push($this->params, $this->display['filter_chr']);
	      array_push($this->params, $this->display['filter_sbp'] * 1000000);
	      array_push($this->params, $this->display['filter_ebp'] * 1000000);
	
	    } else if ($filter == "mark") {
	        if ($this->display['filter_mark'] == "Any") 
		    array_push($this->params, "%/%");
                else 
                    array_push($this->params, $this->display['filter_mark']);
            }
	}
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
                  "est"   => "(ests > ?)",
                  "pe"    => "(pe_radtags > ?)",
                  "blast" => "(blast_hits > ?)",
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

	$this->db['tag_count_sth'] = $this->db['dbh']->prepare($this->queries['tag_count']);
	check_db_error($this->db['tag_count_sth'], __FILE__, __LINE__);

	$this->prepare_filter_parameters();

	$result =& $this->db['tag_count_sth']->execute($this->params);
	check_db_error($result, __FILE__, __LINE__);
	   
	$row = $result->fetchRow();
	$this->num_loci = $row['count'];
    }

    function populate($start_group, $num_groups) {
	//
	// We only want to load genes between $start_gene and $end_gene. 
	//
	$this->db['dbh']->setLimit($num_groups, $start_group);
	check_db_error($this->db['dbh'], __FILE__, __LINE__);

	$this->db['tag_sth'] = $this->db['dbh']->prepare($this->queries['tag']);
	check_db_error($this->db['tag_sth'], __FILE__, __LINE__);

	$this->prepare_filter_parameters();

	//
	// Fetch the results and populate the array of groups.
	//
	$result = $this->db['tag_sth']->execute($this->params);
	check_db_error($result, __FILE__, __LINE__);

	while ($row = $result->fetchRow()) {
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
            $snp_res = $this->db['snp_sth']->execute(array($this->batch, $locus->id));
            check_db_error($snp_res, __FILE__, __LINE__);
            while ($snp_row = $snp_res->fetchRow()) {
                $locus->snps .= $snp_row['col'] . "," . $snp_row['rank_1'] . ">" . $snp_row['rank_2'] . ";";
            }
            $locus->snps = substr($locus->snps, 0, -1);

            $all_res = $this->db['allele_sth']->execute(array($this->batch, $locus->id));
            check_db_error($all_res, __FILE__, __LINE__);
            while ($all_row = $all_res->fetchRow()) {
                $locus->alleles .= $all_row['allele'] . ";";
            }
            $locus->alleles = substr($locus->alleles, 0, -1);

            //
            // Add genotypes
            //
            $gen_res = $this->db['mat_sth']->execute(array($this->batch, $locus->id));
            check_db_error($genres, __FILE__, __LINE__);

            while ($gen_row = $gen_res->fetchRow())
                $locus->add_genotype($gen_row['id'], $gen_row['file'], $gen_row['allele']);

	    $this->loci[$row['tag_id']] = $locus;
	}
    }
}
