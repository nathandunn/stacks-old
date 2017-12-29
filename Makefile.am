bin_PROGRAMS = ustacks cstacks sstacks process_radtags process_shortreads \
	kmer_filter clone_filter populations phasedstacks \
	tsv2bam gstacks

noinst_LIBRARIES = libcore.a libclean.a libpop.a

libcore_a_SOURCES = \
	src/aln_utils.h src/aln_utils.cc \
	src/BamI.h src/BamI.cc \
	src/BamUnalignedI.h \
	src/BowtieI.h \
	src/BustardI.h \
	src/catalog_utils.h src/catalog_utils.cc \
	src/cmb.h src/cmb.cc \
	src/constants.h src/constants.cc \
	src/DNANSeq.h src/DNANSeq.cc \
	src/DNASeq4.h src/DNASeq4.cc \
	src/DNASeq.h src/DNASeq.cc \
	src/FastaI.h \
	src/FastqI.h \
	src/GappedAln.h \
	src/gzFasta.h src/gzFasta.cc \
	src/gzFastq.h \
	src/input.h src/input.cc \
	src/kmers.h src/kmers.cc \
	src/locus.h src/locus.cc \
	src/locus_readers.h \
	src/log_utils.h src/log_utils.cc \
	src/MetaPopInfo.h src/MetaPopInfo.cc \
	src/models.h src/models.cc \
	src/mstack.h src/mstack.cc \
	src/mst.h src/mst.cc \
	src/nucleotides.h src/nucleotides.cc \
	src/renz.h src/renz.cc \
	src/SamI.h \
	src/Seq.h src/Seq.cc \
	src/sql_utilities.h src/sql_utilities.cc \
	src/stacks.h src/stacks.cc \
	src/Tsv.h \
	src/utils.h src/utils.cc \
	src/Vcf.h src/Vcf.cc

libclean_a_SOURCES = \
	src/clean.h src/clean.cc \
	src/file_io.h src/file_io.cc \
	src/write.h src/write.cc

libpop_a_SOURCES = \
	src/bootstrap.h \
	src/export_formats.h src/export_formats.cc \
	src/genotype_dictionaries.h src/genotype_dictionaries.cc \
	src/ordered.h \
	src/PopMap.h \
	src/PopSum.h src/PopSum.cc \
	src/smoothing.h \
	src/smoothing_utils.h

libcore_a_CXXFLAGS = $(OPENMP_CFLAGS) $(BAM_CFLAGS)
libpop_a_CXXFLAGS  = $(OPENMP_CFLAGS) $(BAM_CFLAGS)

phasedstacks_SOURCES = libcore.a src/phasedstacks.h src/phasedstacks.cc

process_radtags_SOURCES    = libcore.a libclean.a src/process_radtags.h src/process_radtags.cc
process_shortreads_SOURCES = libcore.a libclean.a src/process_shortreads.h src/process_shortreads.cc
kmer_filter_SOURCES        = libcore.a libclean.a src/kmer_filter.h src/kmer_filter.cc
clone_filter_SOURCES       = libcore.a libclean.a src/clone_filter.h src/clone_filter.cc

ustacks_SOURCES = libcore.a src/ustacks.h src/ustacks.cc
cstacks_SOURCES = libcore.a src/cstacks.h src/cstacks.cc
sstacks_SOURCES = libcore.a src/sstacks.h src/sstacks.cc

populations_SOURCES = libcore.a libpop.a src/populations.h src/populations.cc

# tsv2bam & gstacks
tsv2bam_SOURCES = libcore.a src/tsv2bam.cc
gstacks_SOURCES = libcore.a src/gstacks.h src/gstacks.cc src/debruijn.h src/debruijn.cc src/Alignment.h src/SuffixTree.h src/SuffixTree.cc

phasedstacks_CXXFLAGS = $(OPENMP_CFLAGS)
phasedstacks_LDFLAGS  = $(OPENMP_CFLAGS)
phasedstacks_LDADD    = libcore.a $(BAM_LIBS)

ustacks_CXXFLAGS = $(OPENMP_CFLAGS)
ustacks_LDFLAGS  = $(OPENMP_CFLAGS)
ustacks_LDADD    = libcore.a $(BAM_LIBS)
cstacks_CXXFLAGS = $(OPENMP_CFLAGS)
cstacks_LDFLAGS  = $(OPENMP_CFLAGS)
cstacks_LDADD    = libcore.a $(BAM_LIBS)
sstacks_CXXFLAGS = $(OPENMP_CFLAGS) $(BAM_CFLAGS)
sstacks_LDFLAGS  = $(OPENMP_CFLAGS)
sstacks_LDADD    = libcore.a $(BAM_LIBS)

process_radtags_CXXFLAGS = $(OPENMP_CFLAGS) $(BAM_CFLAGS)
process_radtags_LDFLAGS  = $(OPENMP_CFLAGS)
process_radtags_LDADD    =  libcore.a libclean.a $(BAM_LIBS)
process_shortreads_CXXFLAGS = $(OPENMP_CFLAGS) $(BAM_CFLAGS)
process_shortreads_LDFLAGS  = $(OPENMP_CFLAGS)
process_shortreads_LDADD    =  libcore.a libclean.a $(BAM_LIBS)
clone_filter_CXXFLAGS = $(OPENMP_CFLAGS) $(BAM_CFLAGS)
clone_filter_LDFLAGS  = $(OPENMP_CFLAGS)
clone_filter_LDADD    = libcore.a libclean.a $(BAM_LIBS)
kmer_filter_CXXFLAGS = $(OPENMP_CFLAGS)
kmer_filter_LDFLAGS  = $(OPENMP_CFLAGS)
kmer_filter_LDADD    = libcore.a libclean.a $(BAM_LIBS)

populations_CXXFLAGS = $(OPENMP_CFLAGS) $(BAM_CFLAGS)
populations_LDFLAGS  = $(OPENMP_CFLAGS)
populations_LDADD    = libcore.a libpop.a $(BAM_LIBS)

# tsv2bam & gstacks
tsv2bam_CXXFLAGS = $(OPENMP_CFLAGS) $(BAM_CFLAGS)
tsv2bam_LDFLAGS  = $(OPENMP_CFLAGS)
tsv2bam_LDADD    = libcore.a $(BAM_LIBS)
gstacks_CXXFLAGS = $(OPENMP_CFLAGS) $(BAM_CFLAGS)
gstacks_LDFLAGS  = $(OPENMP_CFLAGS)
gstacks_LDADD    = libcore.a $(BAM_LIBS)

dist_bin_SCRIPTS = scripts/denovo_map.pl scripts/ref_map.pl	scripts/export_sql.pl \
	scripts/sort_read_pairs.pl scripts/exec_velvet.pl scripts/load_sequences.pl \
	scripts/index_radtags.pl scripts/load_radtags.pl scripts/stacks_export_notify.pl \
	scripts/integrate_alignments.py scripts/count_fixed_catalog_snps.py \
	scripts/stacks-integrate-alignments scripts/stacks-dist-extract scripts/stacks-gdb

dist_noinst_SCRIPTS = autogen.sh scripts/extract_interpop_chars.pl scripts/convert_stacks.pl

nobase_pkgdata_DATA = sql/mysql.cnf.dist sql/catalog_index.sql sql/stacks.sql sql/tag_index.sql sql/chr_index.sql \
	php/CatalogClass.php php/annotate_marker.php php/constants.php.dist php/index.php php/tags.php php/tag.php \
	php/Locus.php php/catalog.php php/correct_genotypes.php php/correct_genotype.php php/export_batch.php php/last_modified.php  \
	php/version.php php/catalog_genotypes.php php/db_functions.php php/header.php php/samples.php \
	php/stacks_functions.php php/view_sequence.php php/sequence_blast.php \
	php/pop_view.php php/sumstat_view.php php/hapstat_view.php php/fst_view.php php/phist_view.php php/stack_view.php \
	php/population_view.js php/ajax.js php/annotate.js php/stacks.js php/export.js php/stacks.css \
	php/images/caret-d.png php/images/caret-u.png php/images/excel_icon.png php/images/l-arrow-disabled.png \
	php/images/l-arrow.png php/images/r-arrow-disabled.png php/images/r-arrow.png php/images/stacks_bg.png \
	php/images/stacks_logo_rev_small.png

LOG_DRIVER = env AM_TAP_AWK='$(AWK)' $(SHELL) $(abs_top_srcdir)/tests/tap-driver.sh

TESTS = tests/process_radtags.t tests/kmer_filter.t tests/ustacks.t tests/pstacks.t

EXTRA_DIST = $(nobase_pkgdata_DATA) LICENSE INSTALL README ChangeLog $(TESTS)

pkglocalstatedir = $(localstatedir)/$(PACKAGE)

SUBDIRS = htslib .

debug:
	$(MAKE) all "CXXFLAGS=-g -Wall -DDEBUG -std=gnu++11"

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
	sed -i.bkp -e 's,_VERSION_,$(VERSION),' $(DESTDIR)$(bindir)/stacks-integrate-alignments
	rm -f $(DESTDIR)$(bindir)/stacks-integrate-alignments.bkp
	chmod +x $(DESTDIR)$(bindir)/stacks-integrate-alignments
	chmod +x $(DESTDIR)$(bindir)/stacks-dist-extract
	chmod +x $(DESTDIR)$(bindir)/stacks-gdb

uninstall-hook:
	rm -rf $(DESTDIR)$(pkglocalstatedir)
