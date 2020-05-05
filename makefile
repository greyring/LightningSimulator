DEBUG=0
CC=gcc
OMP=-fopenmp -DOMP
MPI=-DMPI
MPICC = mpicc
CFLAGS=-g -O3 -Wall -DDEBUG=$(DEBUG)
#CFLAGS=-g -O3 -Wall -DDEBUG=$(DEBUG) -DDYNAMIC
LDFLAGS= -lm
DDIR = ./data

SEQCFILES = light-seq.c graph.c sim-seq.c instrument.c cycletimer.c
OPENMPCFILES = light-openmp.c graph.c sim-openmp.c instrument.c cycletimer.c
MPICFILES = light-mpi.c graph.c sim-mpi.c instrument.c cycletimer.c mpiutil.c
HFILES = graph.h sim.h instrument.h cycletimer.h
MPIHFILES = graph.h sim-mpi.h instrument.h cycletimer.h mpiutil.h

all: light-seq light-openmp light-mpi

light-seq: $(SEQCFILES) $(HFILES)
	$(CC) $(CFLAGS) -o $@ $(SEQCFILES) $(LDFLAGS)

light-openmp: $(OPENMPCFILES) $(HFILES)
	$(CC) $(CFLAGS) $(OMP) -o $@ $(OPENMPCFILES) $(LDFLAGS)

light-mpi: $(MPICFILES) $(MPIHFILES)
	$(MPICC) $(CFLAGS) $(MPI) -o $@ $(MPICFILES) $(LDFLAGS)