DEBUG=0
CC=gcc
MPICC=mpicc
NVCC=nvcc

CFLAGS=-g -O3 -Wall -DDEBUG=$(DEBUG)
#CFLAGS=-g -O3 -Wall -DDEBUG=$(DEBUG) -DDYNAMIC
LDFLAGS= -lm -L/usr/local/depot/cuda-10.2/lib64/ -lcudart

OMP=-fopenmp -DOMP
MPI=-DMPI
NVCCFLAGS=-O3 -m64 --gpu-architecture compute_61

SEQCFILES=light-seq.c graph.c sim-seq.c instrument.c cycletimer.c
OPENMPCFILES=light-openmp.c graph.c sim-openmp.c instrument.c cycletimer.c
MPICFILES=light-mpi.c graph.c sim-mpi.c instrument.c cycletimer.c mpiutil.c
CUDACFILES=light-cuda.c graph.c instrument.c cycletimer.c
CUDAFILES=sim-cuda.cu

HFILES=graph.h sim.h instrument.h cycletimer.h
MPIHFILES=graph.h sim-mpi.h instrument.h cycletimer.h mpiutil.h

all: light-seq light-openmp light-mpi light-cuda

light-seq: $(SEQCFILES) $(HFILES)
	$(CC) $(CFLAGS) -o $@ $(SEQCFILES) $(LDFLAGS)

light-openmp: $(OPENMPCFILES) $(HFILES)
	$(CC) $(CFLAGS) $(OMP) -o $@ $(OPENMPCFILES) $(LDFLAGS)

light-mpi: $(MPICFILES) $(MPIHFILES)
	$(MPICC) $(CFLAGS) $(MPI) -o $@ $(MPICFILES) $(LDFLAGS)

light-cuda: $(CUDACFILES) $(HFILES) sim-cuda.o
	$(CC) $(CFLAGS) -o $@ $(CUDACFILES) sim-cuda.o $(LDFLAGS)

sim-cuda.o: $(CUDAFILES)
	$(NVCC) $(NVCCFLAGS) $(CUDAFILES) -c -o $@