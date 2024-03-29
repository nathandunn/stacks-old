SET default_storage_engine=MYISAM;

create table tag_index (
       id	   int unsigned not null primary key auto_increment,
       batch_id    int unsigned not null,
       sample_id   int unsigned not null,
       tag_id      int unsigned not null,
       con_tag_id  int unsigned not null,
       depth       int unsigned not null,
       snps        int unsigned not null,
       catalog_id  int unsigned not null,
       deleveraged bool default false,
       blacklisted bool default false,
       removed     bool default false,
       INDEX batch_index (batch_id),
       INDEX sample_index (sample_id),
       INDEX tag_index (tag_id),
       INDEX con_tag_index (con_tag_id),
       INDEX depth_index (depth),
       INDEX snps_index (snps),
       INDEX delv_index (deleveraged),
       INDEX black_index (blacklisted),
       INDEX rem_index (removed),
       INDEX catalog_index (catalog_id)
);
