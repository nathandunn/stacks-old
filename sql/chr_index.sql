create table chr_index (
       id    	int unsigned not null primary key auto_increment,
       batch_id	int unsigned not null,
       chr      varchar(32),
       max_len  int unsigned not null
);