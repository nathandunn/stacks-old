bin_PROGRAMS = ustacks pstacks cstacks sstacks hstacks

ustacks_SOURCES = src/ustacks.h src/ustacks.cc src/stacks.h src/stacks.cc src/constants.h \
	src/models.h src/models.cc \
	src/kmers.h src/kmers.cc \
	src/input.h src/input.cc src/Tsv.h src/Bowtie.h src/Fasta.h src/Fastq.h \
	src/cluster.c src/cluster.h

pstacks_SOURCES = src/pstacks.h src/pstacks.cc src/stacks.h src/stacks.cc src/constants.h \
	src/models.h src/models.cc \
	src/input.h src/input.cc src/Tsv.h src/Bowtie.h src/Fasta.h src/Fastq.h

cstacks_SOURCES = src/cstacks.h src/cstacks.cc src/stacks.h src/stacks.cc src/constants.h \
	src/kmers.h src/kmers.cc \
	src/input.h src/input.cc src/sql_utilities.h

hstacks_SOURCES = src/hstacks.h src/hstacks.cc src/stacks.h src/stacks.cc src/constants.h \
	src/kmers.h src/kmers.cc \
	src/input.h src/input.cc src/sql_utilities.h

sstacks_SOURCES = src/sstacks.h src/sstacks.cc src/stacks.h src/stacks.cc src/constants.h \
	src/input.h src/input.cc src/sql_utilities.h

ustacks_CXXFLAGS = $(OPENMP_CFLAGS)
pstacks_CXXFLAGS = $(OPENMP_CFLAGS)
cstacks_CXXFLAGS = $(OPENMP_CFLAGS)
hstacks_CXXFLAGS = $(OPENMP_CFLAGS)
sstacks_CXXFLAGS = $(OPENMP_CFLAGS)

dist_bin_SCRIPTS = scripts/denovo_map.pl scripts/export_sql.pl scripts/genotypes.pl \
	scripts/index_radtags.pl scripts/load_radtags.pl scripts/markers.pl 

dist_noinst_SCRIPTS = autogen.sh scripts/stacks_export_notify.pl php/CatalogClass.php php/annotate_marker.php php/constants.php \
	php/export.php php/index.php php/stacks.css php/tags.php php/Locus.php php/catalog.php php/correct_genotypes.php \
	php/export_batch.php php/last_modified.php php/stacks.js php/version.php php/ajax.js php/catalog_genotypes.php \
	php/db_functions.php php/header.php php/samples.php php/stacks_functions.php php/view_sequence.php \
	php/annotate.js php/catalog_tag.php php/export.js php/sequence_blast.php php/tag.php

EXTRA_DIST = php/images/caret-d.png php/images/caret-u.png php/images/excel_icon.png php/images/l-arrow-disabled.png \
	php/images/l-arrow.png php/images/r-arrow-disabled.png php/images/r-arrow.png php/images/stacks_bg.png \
	php/images/stacks_logo_rev_small.png sql/catalog_index.sql sql/radtags.sql sql/tag_index.sql