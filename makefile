DEBUG=0
CC=gcc
OMP=-fopenmp -DOMP
CFLAGS=-g -O3 -Wall -DDEBUG=$(DEBUG)
#CFLAGS=-g -O3 -Wall -DDEBUG=$(DEBUG) -DDYNAMIC
LDFLAGS= -lm
DDIR = ./data

CFILES = main.c graph.c sim.c instrument.c cycletimer.c
HFILES = graph.h instrument.h cycletimer.h

all: light-seq

light-seq: $(CFILES) $(HFILES)
	$(CC) $(CFLAGS) -o $@ $(CFILES) $(LDFLAGS)