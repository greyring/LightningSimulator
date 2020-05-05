#ifndef __SIM_H__
#define __SIM_H__
#include <stdio.h>
#include "graph.h"
#include "mpiutil.h"
void simulate(int process_count, bool mpi_master, graph_t *g, zonedef_t *zlist, zone_t *z, int count, FILE *ofile);
#endif