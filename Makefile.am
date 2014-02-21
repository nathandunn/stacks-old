bin_PROGRAMS = ustacks pstacks estacks cstacks sstacks hstacks process_radtags process_shortreads \
	kmer_filter clone_filter genotypes populations phasedstacks

ustacks_SOURCES = src/ustacks.h src/ustacks.cc src/stacks.h src/stacks.cc src/constants.h \
	src/mstack.h src/mstack.cc src/mst.h src/mst.cc src/cmb.h src/cmb.cc \
	src/DNASeq.h src/DNASeq.cc src/DNANSeq.h src/DNANSeq.cc \
	src/models.h src/models.cc src/utils.h src/utils.cc \
	src/kmers.h src/kmers.cc \
	src/input.h src/input.cc src/Tsv.h src/BowtieI.h src/FastaI.h src/FastqI.h src/SamI.h \
	src/gzFasta.h src/gzFastq.h

pstacks_SOURCES = src/pstacks.h src/pstacks.cc src/stacks.h src/stacks.cc src/constants.h \
	src/mstack.h src/mstack.cc \
	src/DNANSeq.h src/DNANSeq.cc src/DNASeq.h src/DNASeq.cc \
	src/models.h src/models.cc src/utils.h src/utils.cc \
	src/input.h src/input.cc src/Tsv.h src/BowtieI.h src/FastaI.h src/FastqI.h src/SamI.h \
	src/BamI.h

estacks_SOURCES = src/estacks.h src/estacks.cc src/stacks.h src/stacks.cc src/constants.h \
	src/mstack.h src/mstack.cc \
        src/models.h src/models.cc src/utils.h src/utils.cc \
	src/DNASeq.h src/DNASeq.cc src/DNANSeq.h src/DNANSeq.cc \
        src/input.h src/input.cc src/Tsv.h src/BowtieI.h src/FastaI.h src/FastqI.h src/SamI.h

cstacks_SOURCES = src/cstacks.h src/cstacks.cc src/constants.h \
	src/stacks.h src/stacks.cc src/locus.h src/locus.cc \
	src/kmers.h src/kmers.cc src/utils.h src/utils.cc \
	src/DNASeq.h src/DNASeq.cc src/DNANSeq.h src/DNANSeq.cc \
	src/input.h src/input.cc src/sql_utilities.h

hstacks_SOURCES = src/hstacks.h src/hstacks.cc src/constants.h \
	src/stacks.h src/stacks.cc src/locus.h src/locus.cc \
	src/kmers.h src/kmers.cc \
	src/DNASeq.h src/DNASeq.cc src/DNANSeq.h src/DNANSeq.cc \
	src/input.h src/input.cc src/sql_utilities.h

sstacks_SOURCES = src/sstacks.h src/sstacks.cc src/constants.h \
	src/stacks.h src/stacks.cc src/locus.h src/locus.cc \
	src/DNASeq.h src/DNASeq.cc src/DNANSeq.h src/DNANSeq.cc src/utils.h src/utils.cc \
	src/input.h src/input.cc src/sql_utilities.h

process_radtags_SOURCES = src/process_radtags.h src/process_radtags.cc src/constants.h \
	src/utils.h src/utils.cc src/log_utils.h src/log_utils.cc \
	src/write.h src/write.cc \
	src/clean.h src/clean.cc \
	src/file_io.h src/file_io.cc \
	src/input.h src/input.cc src/BustardI.h src/FastqI.h src/gzFastq.h src/renz.h

process_shortreads_SOURCES = src/process_shortreads.h src/process_shortreads.cc src/constants.h \
	src/clean.h src/clean.cc \
	src/file_io.h src/file_io.cc \
	src/utils.h src/utils.cc src/log_utils.h src/log_utils.cc \
	src/write.h src/write.cc \
	src/kmers.h src/kmers.cc \
	src/input.h src/input.cc src/BustardI.h src/FastqI.h src/gzFastq.h 

kmer_filter_SOURCES = src/kmer_filter.h src/kmer_filter.cc src/constants.h \
	src/utils.h src/utils.cc \
	src/write.h src/write.cc \
	src/kmers.h src/kmers.cc \
	src/input.h src/input.cc src/BustardI.h src/FastqI.h src/FastaI.h \
	src/gzFastq.h src/gzFasta.h

clone_filter_SOURCES = src/clone_filter.h src/clone_filter.cc src/constants.h \
	src/kmers.h src/kmers.cc \
	src/input.h src/input.cc src/BustardI.h src/FastqI.h src/FastaI.h \
	src/gzFastq.h src/gzFasta.h

genotypes_SOURCES = src/genotypes.h src/genotypes.cc src/constants.h \
	src/utils.h src/utils.cc src/catalog_utils.h src/catalog_utils.cc \
	src/log_utils.h src/log_utils.cc \
	src/stacks.h src/stacks.cc src/locus.h src/locus.cc \
	src/DNASeq.h src/DNASeq.cc src/DNANSeq.h src/DNANSeq.cc \
	src/PopMap.h src/genotype_dictionaries.h \
	src/input.h src/input.cc src/sql_utilities.h src/renz.h

populations_SOURCES = src/populations.h src/populations.cc src/constants.h \
	src/utils.h src/utils.cc src/catalog_utils.h src/catalog_utils.cc \
	src/log_utils.h src/log_utils.cc \
	src/stacks.h src/stacks.cc src/locus.h src/locus.cc \
	src/DNASeq.h src/DNASeq.cc src/DNANSeq.h src/DNANSeq.cc \
	src/PopMap.h src/PopSum.h src/genotype_dictionaries.h \
	src/input.h src/input.cc src/sql_utilities.h src/renz.h

phasedstacks_SOURCES = src/phasedstacks.h src/phasedstacks.cc src/constants.h \
	src/locus.h src/locus.cc \
	src/input.h src/input.cc src/sql_utilities.h \
	src/log_utils.h src/log_utils.cc \
	src/utils.h src/utils.cc src/catalog_utils.h src/catalog_utils.cc 

ustacks_CXXFLAGS = $(OPENMP_CFLAGS) $(SPARSEHASH_CFLAGS)
ustacks_LDFLAGS  = $(OPENMP_CFLAGS)
pstacks_CXXFLAGS = $(OPENMP_CFLAGS) $(BAM_CFLAGS)
pstacks_LDFLAGS  = $(OPENMP_CFLAGS)
pstacks_LDADD    = $(BAM_LIBS)
estacks_CXXFLAGS = $(OPENMP_CFLAGS)
estacks_LDFLAGS  = $(OPENMP_CFLAGS)
cstacks_CXXFLAGS = $(OPENMP_CFLAGS) $(SPARSEHASH_CFLAGS)
cstacks_LDFLAGS  = $(OPENMP_CFLAGS)
hstacks_CXXFLAGS = $(OPENMP_CFLAGS)
hstacks_LDFLAGS  = $(OPENMP_CFLAGS)
sstacks_CXXFLAGS = $(OPENMP_CFLAGS)
sstacks_LDFLAGS  = $(OPENMP_CFLAGS)
process_radtags_CXXFLAGS = $(OPENMP_CFLAGS)
process_radtags_LDFLAGS  = $(OPENMP_CFLAGS)
process_shortreads_CXXFLAGS = $(OPENMP_CFLAGS)
process_shortreads_LDFLAGS  = $(OPENMP_CFLAGS)
kmer_filter_CXXFLAGS  = $(OPENMP_CFLAGS)
kmer_filter_LDFLAGS   = $(OPENMP_CFLAGS)
genotypes_CXXFLAGS    = $(OPENMP_CFLAGS)
genotypes_LDFLAGS     = $(OPENMP_CFLAGS)
populations_CXXFLAGS  = $(OPENMP_CFLAGS)
populations_LDFLAGS   = $(OPENMP_CFLAGS)
phasedstacks_CXXFLAGS = $(OPENMP_CFLAGS)
phasedstacks_LDFLAGS  = $(OPENMP_CFLAGS)

dist_bin_SCRIPTS = scripts/denovo_map.pl scripts/ref_map.pl scripts/export_sql.pl \
	scripts/sort_read_pairs.pl scripts/exec_velvet.pl scripts/load_sequences.pl \
	scripts/index_radtags.pl scripts/load_radtags.pl scripts/stacks_export_notify.pl

dist_noinst_SCRIPTS = autogen.sh scripts/extract_interpop_chars.pl

nobase_pkgdata_DATA = sql/mysql.cnf.dist sql/catalog_index.sql sql/stacks.sql sql/tag_index.sql sql/chr_index.sql \
	php/CatalogClass.php php/annotate_marker.php php/constants.php.dist php/index.php php/tags.php \
	php/Locus.php php/catalog.php php/correct_genotypes.php php/correct_genotype.php php/export_batch.php php/last_modified.php  \
	php/version.php php/catalog_genotypes.php php/db_functions.php php/header.php php/samples.php \
	php/stacks_functions.php php/view_sequence.php php/catalog_tag.php php/sequence_blast.php php/tag.php \
	php/ajax.js php/annotate.js php/stacks.js php/export.js php/stacks.css \
	php/images/caret-d.png php/images/caret-u.png php/images/excel_icon.png php/images/l-arrow-disabled.png \
	php/images/l-arrow.png php/images/r-arrow-disabled.png php/images/r-arrow.png php/images/stacks_bg.png \
	php/images/stacks_logo_rev_small.png

LOG_DRIVER = env AM_TAP_AWK='$(AWK)' $(SHELL) $(abs_top_srcdir)/config/tap-driver.sh
TESTS = tests/process_radtags.t tests/kmer_filter.t

EXTRA_DIST = $(nobase_pkgdata_DATA) LICENSE INSTALL README ChangeLog $(TESTS)

pkglocalstatedir = $(localstatedir)/$(PACKAGE)

debug:
	$(MAKE) all "CXXFLAGS=-g -Wall -DDEBUG"

install-data-hook:
	sed -e 's,_VERSION_,$(VERSION),' -e 's,_BINDIR_,$(bindir)/,g' -e 's,_PKGDATADIR_,$(pkgdatadir)/,g' $(DESTDIR)$(bindir)/denovo_map.pl > $(DESTDIR)$(bindir)/denovo_map.pl.subst
	mv $(DESTDIR)$(bindir)/denovo_map.pl.subst $(DESTDIR)$(bindir)/denovo_map.pl
	chmod +x $(DESTDIR)$(bindir)/denovo_map.pl
	sed -e 's,_VERSION_,$(VERSION),' -e 's,_BINDIR_,$(bindir)/,g' -e 's,_PKGDATADIR_,$(pkgdatadir)/,g' $(DESTDIR)$(bindir)/ref_map.pl > $(DESTDIR)$(bindir)/ref_map.pl.subst
	mv $(DESTDIR)$(bindir)/ref_map.pl.subst $(DESTDIR)$(bindir)/ref_map.pl
	chmod +x $(DESTDIR)$(bindir)/ref_map.pl
	sed -e 's,_VERSION_,$(VERSION),' -e 's,_PKGDATADIR_,$(pkgdatadir)/,g' $(DESTDIR)$(bindir)/export_sql.pl > $(DESTDIR)$(bindir)/export_sql.pl.subst
	mv $(DESTDIR)$(bindir)/export_sql.pl.subst $(DESTDIR)$(bindir)/export_sql.pl
	chmod +x $(DESTDIR)$(bindir)/export_sql.pl
	sed -e 's,_VERSION_,$(VERSION),' -e 's,_PKGDATADIR_,$(pkgdatadir)/,g' $(DESTDIR)$(bindir)/index_radtags.pl > $(DESTDIR)$(bindir)/index_radtags.pl.subst
	mv $(DESTDIR)$(bindir)/index_radtags.pl.subst $(DESTDIR)$(bindir)/index_radtags.pl
	chmod +x $(DESTDIR)$(bindir)/index_radtags.pl
	sed -e 's,_VERSION_,$(VERSION),' -e 's,_PKGDATADIR_,$(pkgdatadir)/,g' $(DESTDIR)$(bindir)/load_radtags.pl > $(DESTDIR)$(bindir)/load_radtags.pl.subst
	mv $(DESTDIR)$(bindir)/load_radtags.pl.subst $(DESTDIR)$(bindir)/load_radtags.pl
	chmod +x $(DESTDIR)$(bindir)/load_radtags.pl
	sed -e 's,_VERSION_,$(VERSION),' $(DESTDIR)$(bindir)/sort_read_pairs.pl > $(DESTDIR)$(bindir)/sort_read_pairs.pl.subst
	mv $(DESTDIR)$(bindir)/sort_read_pairs.pl.subst $(DESTDIR)$(bindir)/sort_read_pairs.pl
	chmod +x $(DESTDIR)$(bindir)/sort_read_pairs.pl
	sed -e 's,_VERSION_,$(VERSION),' $(DESTDIR)$(bindir)/exec_velvet.pl > $(DESTDIR)$(bindir)/exec_velvet.pl.subst
	mv $(DESTDIR)$(bindir)/exec_velvet.pl.subst $(DESTDIR)$(bindir)/exec_velvet.pl
	chmod +x $(DESTDIR)$(bindir)/exec_velvet.pl
	sed -e 's,_VERSION_,$(VERSION),' -e 's,_PKGDATADIR_,$(pkgdatadir)/,g' $(DESTDIR)$(bindir)/load_sequences.pl > $(DESTDIR)$(bindir)/load_sequences.pl.subst
	mv $(DESTDIR)$(bindir)/load_sequences.pl.subst $(DESTDIR)$(bindir)/load_sequences.pl
	chmod +x $(DESTDIR)$(bindir)/load_sequences.pl
	sed -e 's,_VERSION_,$(VERSION),' -e 's,_PKGDATADIR_,$(pkgdatadir)/,g' -e 's,_BINDIR_,$(bindir)/,g' \
		$(DESTDIR)$(bindir)/stacks_export_notify.pl > $(DESTDIR)$(bindir)/stacks_export_notify.pl.subst
	mv $(DESTDIR)$(bindir)/stacks_export_notify.pl.subst $(DESTDIR)$(bindir)/stacks_export_notify.pl
	chmod +x $(DESTDIR)$(bindir)/stacks_export_notify.pl
	sed -e 's,_PKGDATADIR_,$(pkgdatadir)/,g' -e 's,_BINDIR_,$(bindir)/,g' \
		$(DESTDIR)$(pkgdatadir)/php/constants.php.dist > $(DESTDIR)$(pkgdatadir)/php/constants.php.dist.subst
	mv $(DESTDIR)$(pkgdatadir)/php/constants.php.dist.subst $(DESTDIR)$(pkgdatadir)/php/constants.php.dist
	echo $(VERSION) > $(DESTDIR)$(pkgdatadir)/php/version.php
	$(install_sh) -d -m 755 $(DESTDIR)$(pkgdatadir)/php/export

uninstall-hook:
	rm -rf $(DESTDIR)$(pkglocalstatedir)
