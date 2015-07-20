SET default_storage_engine=MYISAM;

create table batches (
       id           int unsigned not null primary key auto_increment,
       date         DATE not null,
       description  tinytext,
       type         enum('map', 'population')
);

create table samples (
       id        int unsigned not null primary key auto_increment,
       batch_id	 int unsigned not null,
       sample_id int unsigned not null,
       type      enum('parent', 'progeny', 'sample'),
       file      varchar(128),
       pop_id	 varchar(32),
       group_id	 varchar(32)
);

create table populations (
       id        int unsigned not null primary key auto_increment,
       batch_id	 int unsigned not null,
       pop_id    int unsigned not null,
       pop_name  varchar(128)
);

create table sumstats (
       batch_id	  int unsigned not null,
       tag_id	  int unsigned not null,
       chr	  varchar(32),
       bp	  int unsigned default 0,
       col	  int unsigned default 0,
       pop_id	  varchar(32),
       p_nuc      char,
       q_nuc      char,
       n	  int unsigned default 0,
       p	  double,
       obs_het	  double,
       obs_hom	  double,
       exp_het	  double,
       exp_hom	  double,
       pi	  double,
       pi_s       double,
       pi_s_pval  double,
       fis	  double,
       fis_s      double,
       fis_s_pval double,
       private    int unsigned default 0,
       INDEX      batch_id_index (batch_id),
       INDEX      tag_id_index (tag_id)
);

create table fst (
       batch_id	   int unsigned not null,
       tag_id	   int unsigned not null,
       pop_id_1	   varchar(32),
       pop_id_2	   varchar(32),
       chr	   varchar(32),
       bp	   int unsigned default 0,
       col	   int unsigned default 0,
       pi_o	   double,
       fst	   double,
       fishers_p   double,
       odds_ratio  double,
       ci_low	   double,
       ci_high	   double,
       lod	   double,
       fst_c	   double,
       fst_s	   double,
       amova_fst   double, 
       amova_fst_c double,
       amova_fst_s double,
       amova_fst_s_pval double,
       snp_cnt     int unsigned default 0,
       INDEX       batch_id_index (batch_id),
       INDEX       tag_id_index (tag_id)
);

create table hapstats (
       batch_id	     int unsigned not null,
       tag_id	     int unsigned not null,
       chr	     varchar(32),
       bp	     int unsigned default 0,
       pop_id	     varchar(32),
       n	     int unsigned default 0,
       hapcnt	     double,
       gene_div	     double,
       gene_div_s    double,
       gene_div_pval double,
       hap_div	     double,
       hap_div_s     double,
       hap_div_pval  double,
       haplotypes    tinytext,
       INDEX         batch_id_index (batch_id),
       INDEX         tag_id_index (tag_id)
);

create table phist (
       batch_id	   int unsigned not null,
       tag_id	   int unsigned not null,
       pop_id_1	   varchar(32),
       pop_id_2	   varchar(32),
       chr	   varchar(32),
       bp	   int unsigned default 0,
       phist	   double,
       phist_s     double,
       phist_pval  double,
       fpst	   double,
       fpst_s      double,
       fpst_pval   double,
       dest	   double,
       dest_s      double,
       dest_pval   double,
       INDEX       batch_id_index (batch_id),
       INDEX       tag_id_index (tag_id)
);

create table catalog_tags (
       id    	    int unsigned not null primary key auto_increment,
       batch_id     int unsigned not null,
       tag_id	    int unsigned not null,
       chr          varchar(32),
       bp           int unsigned default 0,
       strand	    enum('+', '-'),
       relationship enum('consensus', 'primary', 'secondary', 'tertiary'),
       sub_id	    int unsigned not null,
       merge_type   tinytext,
       seq 	    text,
       INDEX        batch_id_index (batch_id),
       INDEX        tag_id_index (tag_id)
);

create table catalog_snps (
       id    	   int unsigned not null primary key auto_increment,
       batch_id    int unsigned not null,
       tag_id      int unsigned not null,
       col         int unsigned not null,
       type        enum('E', 'O', 'U'),
       lratio      float,
       rank_1	   char(1),
       rank_2	   char(1),
       rank_3	   char(1),
       rank_4	   char(1),
       INDEX       batch_index (batch_id),
       INDEX       tag_index (tag_id)
);

create table catalog_alleles (
       id    	   int unsigned not null primary key auto_increment,
       batch_id    int unsigned not null,
       tag_id      int unsigned not null,
       allele      varchar(32),
       read_pct    float,
       read_cnt    int unsigned,
       INDEX       batch_index (batch_id),
       INDEX       tag_index (tag_id)
);

create table catalog_genotypes (
       id    	   int unsigned not null primary key auto_increment,
       batch_id    int unsigned not null,
       catalog_id  int unsigned not null,
       sample_id   int unsigned not null,
       genotype    char(2),
       INDEX       batch_index (batch_id),
       INDEX       catalog_index (catalog_id),
       INDEX       sample_index (sample_id)
);

create table genotype_corrections (
       id    	   int unsigned not null primary key auto_increment,
       batch_id    int unsigned not null,
       catalog_id  int unsigned not null,
       sample_id   int unsigned not null,
       genotype    char(2),
       INDEX       cat_index (catalog_id),
       INDEX       batch_index (batch_id),
       INDEX       sample_index (sample_id)
);

create table catalog_annotations (
       id    	   int unsigned not null primary key auto_increment,
       batch_id    int unsigned not null,
       catalog_id  int unsigned not null,
       external_id varchar(64),
       INDEX       batch_index (batch_id),
       INDEX       catalog_index (catalog_id),
       INDEX       external_index (external_id)
);

create table pileup (
       id    	    int unsigned not null primary key auto_increment,
       sample_id    int unsigned not null,
       tag_id	    int unsigned not null,
       chr          varchar(32),
       bp           int unsigned default 0,
       relationship enum('consensus', 'model', 'primary', 'secondary'),
       sub_id	    int unsigned not null,
       seq_id	    varchar(32),
       seq 	    text,
       INDEX        tag_id_index (tag_id),
       INDEX        sample_id_index (sample_id),
       INDEX        rel_index (relationship)
);

create table unique_tags (
       id    	    int unsigned not null primary key auto_increment,
       sample_id    int unsigned not null,
       tag_id	    int unsigned not null,
       chr          varchar(32),
       bp           int unsigned default 0,
       strand	    enum('+', '-'),
       relationship enum('consensus', 'model', 'primary', 'secondary'),
       sub_id	    int unsigned not null,
       seq_id	    varchar(32),
       seq 	    text,
       deleveraged  bool default false,
       blacklisted  bool default false,
       removed      bool default false,
       lnl	    float,
       INDEX        tag_id_index (tag_id),
       INDEX        sample_id_index (sample_id),
       INDEX        rel_index (relationship)
);

create table snps (
       id    	   int unsigned not null primary key auto_increment,
       sample_id   int unsigned not null,
       tag_id      int unsigned not null,
       col         int unsigned not null,
       type        enum('E', 'O', 'U'),
       lratio      float,
       rank_1	   char(1),
       rank_2	   char(1),
       rank_3	   char(1),
       rank_4	   char(1),
       INDEX       samp_index (sample_id),
       INDEX       tag_index (tag_id)
);

create table alleles (
       id    	   int unsigned not null primary key auto_increment,
       sample_id   int unsigned not null,
       tag_id      int unsigned not null,
       allele      varchar(32),
       read_pct    float,
       read_cnt    int unsigned,
       INDEX       samp_index (sample_id),
       INDEX       tag_index (tag_id)
);

create table matches (
       id    	   int unsigned not null primary key auto_increment,
       batch_id    int unsigned not null,
       catalog_id  int unsigned not null,
       sample_id   int unsigned not null,
       tag_id      int unsigned not null,
       allele	   varchar(256),
       depth       int unsigned not null,
       lnl	   float,
       INDEX	   batch_id_index (batch_id),
       INDEX	   catalog_id_index (catalog_id),
       INDEX	   sample_id_index (sample_id),
       INDEX	   tag_id_index (tag_id)
);

create table markers (
       id         int unsigned not null primary key auto_increment,
       batch_id   int unsigned not null,
       catalog_id int unsigned not null,
       type       enum('aa/bb', 'aa/b-', '-a/bb', 'ab/--', '--/ab', 'aa/ab', 'ab/aa', 'ab/a-', '-a/ab', 'ab/ab', 'ab/ac', 'ab/cd', 'ab/cc', 'cc/ab', 'ab/c-', '-c/ab'),
       progeny    int unsigned not null default 0,
       max_pct    float,
       ratio      varchar(512),
       chisq_pval float,
       lnl        float,
       geno_map   text,
       uncor_type enum('aa/bb', 'aa/b-', '-a/bb', 'ab/--', '--/ab', 'aa/ab', 'ab/aa', 'ab/a-', '-a/ab', 'ab/ab', 'ab/ac', 'ab/cd', 'ab/cc', 'cc/ab', 'ab/c-', '-c/ab')
);

create table sequence (
       id         int unsigned not null primary key auto_increment,
       batch_id   int unsigned not null,
       catalog_id int unsigned not null,
       type       enum('pe_radtag', 'est'),
       seq_id     varchar(64),
       seq        text,
       INDEX      catalog_id_index (catalog_id)
);

create table sequence_blast (
       id               int unsigned not null primary key auto_increment,
       batch_id         int unsigned not null default 0,
       catalog_id       int unsigned not null default 0,
       seq_id           int unsigned not null default 0,
       algorithm        enum('blastn', 'blastp', 'blastx', 'tblastn', 'tblastx'),
       query_id         varchar(64),
       query_len        int unsigned not null default 0,
       hit_id           varchar(128),
       hit_len          int unsigned not null default 0,
       score            double,
       e_value          double,
       percent_ident    double,
       hsp_rank         int unsigned not null default 0,
       aln_len          int unsigned not null default 0,
       aln_homology_str text,
       query_aln        text,
       query_aln_start  int unsigned not null default 0,
       query_aln_end    int unsigned not null default 0,
       hit_aln          text,
       hit_aln_start    int unsigned not null default 0,
       hit_aln_end      int unsigned not null default 0
);

create table chr_index (
       id    	int unsigned not null primary key auto_increment,
       batch_id	int unsigned not null,
       chr      varchar(32),
       max_len  int unsigned not null
);

create table ref_radome (
       id       int unsigned not null primary key,
       chr	varchar(32),
       bp	int unsigned not null,
       strand	enum('1', '-1'),
       type	enum('genomic', 'exon', 'intron'),
       gene	varchar(32),
       ext_id	varchar(32),
       exon	varchar(32),
       ex_start	int unsigned not null,
       ex_end	int unsigned not null,
       ex_index	int unsigned not null,
       tran_cnt int unsigned not null,
       dist	int unsigned not null,
       uniq	bool default false,
       tran	varchar(32)
);
