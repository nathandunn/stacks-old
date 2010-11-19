bin_PROGRAMS = ustacks pstacks cstacks sstacks hstacks

ustacks_SOURCES = src/ustacks.h src/ustacks.cc src/stacks.h src/stacks.cc src/constants.h \
	src/models.h src/models.cc \
	src/kmers.h src/kmers.cc \
	src/input.h src/input.cc src/Tsv.h src/Bowtie.h src/Fasta.h src/Fastq.h src/Sam.h \
	src/cluster.c src/cluster.h

pstacks_SOURCES = src/pstacks.h src/pstacks.cc src/stacks.h src/stacks.cc src/constants.h \
	src/models.h src/models.cc \
	src/input.h src/input.cc src/Tsv.h src/Bowtie.h src/Fasta.h src/Fastq.h src/Sam.h 

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

dist_noinst_SCRIPTS = autogen.sh scripts/stacks_export_notify.pl 

nobase_pkgdata_DATA = sql/mysql.cnf.dist sql/catalog_index.sql sql/stacks.sql sql/tag_index.sql \
	php/CatalogClass.php php/annotate_marker.php php/constants.php php/export.php php/index.php php/tags.php \
	php/Locus.php php/catalog.php php/correct_genotypes.php php/export_batch.php php/last_modified.php  \
	php/version.php php/catalog_genotypes.php php/db_functions.php php/header.php php/samples.php \
	php/stacks_functions.php php/view_sequence.php php/catalog_tag.php php/sequence_blast.php php/tag.php \
	php/ajax.js php/annotate.js php/stacks.js php/export.js php/stacks.css \
	php/images/caret-d.png php/images/caret-u.png php/images/excel_icon.png php/images/l-arrow-disabled.png \
	php/images/l-arrow.png php/images/r-arrow-disabled.png php/images/r-arrow.png php/images/stacks_bg.png \
	php/images/stacks_logo_rev_small.png 

EXTRA_DIST = $(nobase_pkgdata_DATA) LICENSE INSTALL

pkglocalstatedir = $(localstatedir)/$(PACKAGE)

install-data-hook:
	sed -e 's,_BINDIR_,$(bindir)/,g' $(DESTDIR)$(bindir)/denovo_map.pl > $(DESTDIR)$(bindir)/denovo_map.pl.subst
	mv $(DESTDIR)$(bindir)/denovo_map.pl.subst $(DESTDIR)$(bindir)/denovo_map.pl
	chmod +x $(DESTDIR)$(bindir)/denovo_map.pl
	sed -e 's,_PKGDATADIR_,$(pkgdatadir)/,g' $(DESTDIR)$(bindir)/export_sql.pl > $(DESTDIR)$(bindir)/export_sql.pl.subst
	mv $(DESTDIR)$(bindir)/export_sql.pl.subst $(DESTDIR)$(bindir)/export_sql.pl
	chmod +x $(DESTDIR)$(bindir)/export_sql.pl
	sed -e 's,_PKGDATADIR_,$(pkgdatadir)/,g' $(DESTDIR)$(bindir)/index_radtags.pl > $(DESTDIR)$(bindir)/index_radtags.pl.subst
	mv $(DESTDIR)$(bindir)/index_radtags.pl.subst $(DESTDIR)$(bindir)/index_radtags.pl
	chmod +x $(DESTDIR)$(bindir)/index_radtags.pl
	sed -e 's,_PKGDATADIR_,$(pkgdatadir)/,g' $(DESTDIR)$(bindir)/markers.pl > $(DESTDIR)$(bindir)/markers.pl.subst
	mv $(DESTDIR)$(bindir)/markers.pl.subst $(DESTDIR)$(bindir)/markers.pl
	chmod +x $(DESTDIR)$(bindir)/markers.pl
	sed -e 's,_PKGDATADIR_,$(pkgdatadir)/,g' -e 's,_LOCALSTATEDIR_,$(pkglocalstatedir)/,g' \
		$(DESTDIR)$(bindir)/genotypes.pl > $(DESTDIR)$(bindir)/genotypes.pl.subst
	mv $(DESTDIR)$(bindir)/genotypes.pl.subst $(DESTDIR)$(bindir)/genotypes.pl
	chmod +x $(DESTDIR)$(bindir)/genotypes.pl
	echo $(VERSION) > $(DESTDIR)$(pkgdatadir)/php/version.php
	$(install_sh) -d -m 755 $(DESTDIR)$(pkglocalstatedir)/caches

uninstall-hook:
	rm -rf $(DESTDIR)$(pkglocalstatedir)
