# Stacks Makefile
#
# Julian Catchen
#

OMP = -fopenmp
#OMP =
CXX = g++
CC  = gcc
CXXFLAGS += -g -Wall ${OMP}
LDFLAGS += 
TARGET = pstacks cstacks exstacks

all: cstacks pstacks exstacks

exstacks: exstacks.o input.o stacks.o
	${CXX} -fopenmp -o exstacks exstacks.o input.o stacks.o

exstacks.o:	exstacks.h exstacks.cc stacks.h constants.h sql_utilities.h
	${CXX} ${CXXFLAGS} -c exstacks.cc

cstacks: cstacks.o input.o stacks.o
	${CXX} -fopenmp -o cstacks cstacks.o input.o stacks.o

cstacks.o:	cstacks.h cstacks.cc stacks.h constants.h sql_utilities.h
	${CXX} ${CXXFLAGS} -c cstacks.cc

pstacks: pstacks.o input.o models.o stacks.o
	${CXX} -fopenmp -o pstacks pstacks.o input.o stacks.o models.o

pstacks.o:	pstacks.h pstacks.cc stacks.h constants.h Tsv.h Bowtie.h Sam.h Fasta.h Fastq.h
	${CXX} ${CXXFLAGS} -c pstacks.cc

input.o:	input.h input.cc constants.h
	${CXX} ${CXXFLAGS} -c input.cc

stacks.o:	stacks.h stacks.cc constants.h
	${CXX} ${CXXFLAGS} -c stacks.cc

models.o:	models.h models.cc
	${CXX} ${CXXFLAGS} -c models.cc

.PHONY:	clean

clean:
	-/bin/rm -f ${TARGET} *.o *~ core

