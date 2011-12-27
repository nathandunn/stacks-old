create table catalog_index (
       id	int unsigned not null primary key auto_increment,
       batch_id int unsigned not null,
       cat_id   int unsigned not null,
       tag_id   int unsigned not null,
       snps     int unsigned not null,
       parents  int unsigned not null,
       progeny  int unsigned not null,
       alleles  int unsigned not null,
       marker   enum('aa/bb', 'ab/--', '--/ab', 'aa/ab', 'ab/aa', 'ab/ab', 'ab/ac', 'ab/cd', 'ab/cc', 'cc/ab', ''),
       valid_progeny int unsigned not null default 0,
       max_pct    float not null default 0,
       ratio      varchar(256),
       ests       int unsigned not null,
       pe_radtags int unsigned not null,
       blast_hits int unsigned not null,
       geno_cnt   int unsigned not null,
       chr        varchar(32),
       bp         int unsigned default 0,
       type	  enum('genomic', 'exon', 'intron'),
       ref_id     int unsigned not null,
       INDEX batch_index (batch_id),
       INDEX tag_index (tag_id),
       INDEX snps_index (snps),
       INDEX parents_index (parents),
       INDEX progeny_index (progeny),
       INDEX allele_index (alleles),
       INDEX marker_index (marker),
       INDEX valid_index (valid_progeny),
       INDEX max_pct_index (max_pct),
       INDEX ests_index (ests),
       INDEX pe_rad_index (pe_radtags),
       INDEX hits_index (blast_hits),
       INDEX geno_index (geno_cnt),
       INDEX chr_index (chr),
       INDEX bp_index (bp),
       INDEX type_index (type)
);
