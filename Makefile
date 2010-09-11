# Stacks Makefile
#
# Julian Catchen
#

#OMP = -fopenmp
OMP =
CXX = g++
CC  = gcc
CXXFLAGS += -g -Wall ${OMP}
LDFLAGS += 
TARGET = ustacks pstacks cstacks sstacks exstacks

all: ustacks cstacks pstacks sstacks exstacks

exstacks: exstacks.o input.o stacks.o
	${CXX} -fopenmp -o exstacks exstacks.o input.o stacks.o

exstacks.o:	exstacks.h exstacks.cc stacks.h constants.h sql_utilities.h
	${CXX} ${CXXFLAGS} -c exstacks.cc

cstacks: cstacks.o input.o stacks.o kmers.o
	${CXX} -fopenmp -o cstacks cstacks.o input.o stacks.o kmers.o

cstacks.o:	cstacks.h cstacks.cc stacks.h kmers.h constants.h sql_utilities.h
	${CXX} ${CXXFLAGS} -c cstacks.cc

sstacks: sstacks.o input.o stacks.o
	${CXX} -fopenmp -o sstacks sstacks.o input.o stacks.o

sstacks.o:	sstacks.h sstacks.cc stacks.h constants.h sql_utilities.h
	${CXX} ${CXXFLAGS} -c sstacks.cc

pstacks: pstacks.o input.o models.o stacks.o
	${CXX} -fopenmp -o pstacks pstacks.o input.o stacks.o models.o

pstacks.o:	pstacks.h pstacks.cc stacks.h constants.h Tsv.h Bowtie.h Sam.h Fasta.h Fastq.h
	${CXX} ${CXXFLAGS} -c pstacks.cc

ustacks: ustacks.o input.o models.o stacks.o kmers.o cluster.o
	${CXX} -fopenmp -o ustacks ustacks.o input.o stacks.o models.o kmers.o cluster.o

ustacks.o:	ustacks.h ustacks.cc stacks.h constants.h sql_utilities.h
	${CXX} ${CXXFLAGS} -c ustacks.cc

input.o:	input.h input.cc constants.h
	${CXX} ${CXXFLAGS} -c input.cc

cluster.o:	cluster.h cluster.c
	${CC} ${CXXFLAGS} -c cluster.c

stacks.o:	stacks.h stacks.cc constants.h
	${CXX} ${CXXFLAGS} -c stacks.cc

kmers.o:	kmers.h kmers.cc stacks.h constants.h
	${CXX} ${CXXFLAGS} -c kmers.cc

models.o:	models.h models.cc
	${CXX} ${CXXFLAGS} -c models.cc

.PHONY:	clean

clean:
	-/bin/rm -f ${TARGET} *.o *~ core

