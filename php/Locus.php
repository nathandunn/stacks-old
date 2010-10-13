<?php

class Locus {
    var $id;
    var $seq;
    var $annotation;
    var $chr;
    var $bp;
    var $snps;
    var $alleles;

    var $genotypes;   // Array of samples that possess this locus
    var $num_parents;
    var $num_progeny;
    var $num_snps;
    var $num_alleles;
    var $num_ests;
    var $num_pe_tags;
    var $num_blast;
    var $valid_progeny;
    var $marker;

    function Locus() {
	$this->genotypes     = array();
	$this->annotation    = "";
        $this->chr           = "";
        $this->bp            = 0;
        $this->marker        = "";
	$this->snps          = "";
	$this->alleles       = "";
        $this->num_parents   = 0;
        $this->num_progeny   = 0;
        $this->valid_progeny = 0;
        $this->num_alleles   = 0;
        $this->num_snps      = 0;
        $this->num_ests      = 0;
        $this->num_pe_tags   = 0;
        $this->num_blast     = 0;
    }

    function &genotypes() {
	return $this->genotypes;
    }

    function &genotype($id) {
        if (isset($this->genotypes[$id]))
            return $this->genotypes[$id];
        else
            return NULL;
    }

    function add_genotype($id, $file, $allele) {
        $a = array('file'   => $file, 
                   'allele' => $allele, 
                   'tag_id' => $id);
        if (!isset($this->genotypes[$file]))
            $this->genotypes[$file] = array();

        array_push($this->genotypes[$file], $a);
    }
}

?>
