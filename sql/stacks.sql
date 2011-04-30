create table batches (
       id           int unsigned not null primary key auto_increment,
       date         DATE not null,
       description  tinytext
);

create table samples (
       id        int unsigned not null primary key auto_increment,
       batch_id	 int unsigned not null,
       sample_id int unsigned not null,
       type      enum('parent', 'progeny', 'sample'),
       file      varchar(128)
);

create table catalog_tags (
       id    	    int unsigned not null primary key auto_increment,
       batch_id     int unsigned not null,
       tag_id	    int unsigned not null,
       chr          varchar(32),
       bp           int unsigned default 0,
       relationship enum('consensus', 'primary', 'secondary', 'tertiary'),
       sub_id	    int unsigned not null,
       merge_type   tinytext,
       seq 	    tinytext,
       INDEX        batch_id_index (batch_id),
       INDEX        tag_id_index (tag_id)
);

create table catalog_snps (
       id    	   int unsigned not null primary key auto_increment,
       batch_id    int unsigned not null,
       tag_id      int unsigned not null,
       col         int unsigned not null,
       lratio      float,
       rank_1	   char(1),
       rank_2	   char(1),
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

create table unique_tags (
       id    	    int unsigned not null primary key auto_increment,
       sample_id    int unsigned not null,
       tag_id	    int unsigned not null,
       chr          varchar(32),
       bp           int unsigned default 0,
       relationship enum('consensus', 'primary', 'secondary', 'tertiary'),
       sub_id	    int unsigned not null,
       seq_id	    varchar(32),
       seq 	    tinytext,
       deleveraged  bool default false,
       blacklisted  bool default false,
       removed      bool default false,
       INDEX        tag_id_index (tag_id),
       INDEX        sample_id_index (sample_id),
       INDEX        rel_index (relationship)
);

create table snps (
       id    	   int unsigned not null primary key auto_increment,
       sample_id   int unsigned not null,
       tag_id      int unsigned not null,
       col         int unsigned not null,
       lratio      float,
       rank_1	   char(1),
       rank_2	   char(1),
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
       allele	   varchar(32),
       INDEX	   batch_id_index (batch_id),
       INDEX	   catalog_id_index (catalog_id),
       INDEX	   sample_id_index (sample_id),
       INDEX	   tag_id_index (tag_id)
);

create table markers (
       id         int unsigned not null primary key auto_increment,
       batch_id   int unsigned not null,
       catalog_id int unsigned not null,
       type       enum('aa/bb', 'ab/--', '--/ab', 'aa/ab', 'ab/aa', 'ab/ab', 'ab/ac', 'ab/cd', 'ab/cc', 'cc/ab'),
       progeny    int unsigned not null default 0,
       max_pct    float,
       ratio      varchar(128)
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
