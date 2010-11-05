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

ustacks_CPPFLAGS = -fopenmp
ustacks_LDFLAGS = -fopenmp

pstacks_CPPFLAGS = -fopenmp
pstacks_LDFLAGS = -fopenmp

cstacks_CPPFLAGS = -fopenmp
cstacks_LDFLAGS = -fopenmp

hstacks_CPPFLAGS = -fopenmp
hstacks_LDFLAGS = -fopenmp

sstacks_CPPFLAGS = -fopenmp
sstacks_LDFLAGS = -fopenmp

dist_noinst_SCRIPTS = autogen.sh

