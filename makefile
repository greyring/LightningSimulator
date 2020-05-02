DEBUG=0
CC=gcc
OMP=-fopenmp -DOMP
CFLAGS=-g -O3 -Wall -DDEBUG=$(DEBUG)
#CFLAGS=-g -O3 -Wall -DDEBUG=$(DEBUG) -DDYNAMIC
LDFLAGS= -lm
DDIR = ./data

SEQCFILES = light-seq.c graph.c sim-seq.c instrument.c cycletimer.c
SEQHFILES = graph.h sim.h instrument.h cycletimer.h

all: light-seq

light-seq: $(SEQCFILES) $(SEQHFILES)
	$(CC) $(CFLAGS) -o $@ $(SEQCFILES) $(LDFLAGS)