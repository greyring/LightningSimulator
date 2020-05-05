#ifndef __MPIUTIL_H__
#define __MPIUTIL_H__
#include <mpi.h>
#include "graph.h"

typedef struct {
    int this_zone; // never used

    int gheight; // used in init update charge
    int gwidth; // used in init update charge

    int start_row;
    int start_col;
    int width;
    int height;
    int eta; // shape of lightning
    int adj[4]; // zoneid of up, left, right, down used in exchange charges

    int power; // used in scatter power

    // electrical potential
    double *charge;
    double *charge_buffer;
    double *ghost_charge;

    // charge density (poisson equation)
    double *boundary;

    int *reset_bolt;
    int *bolt;
    int num_choice;
    int *choice_idxs; // choosed point, zoneidx

    // MPI buffer
    double *left_buf; // send buf used in exchange charge
    double *right_buf; // send buf used in exchange charge
    double *prob_buf; // send buf used in gatter_probs
    MPI_Request mpi_r;
}zone_t;

typedef struct {
    int start_row;
    int start_col;
    int width;
    int height;
    int eta;
    int adj[4]; // zoneid of up, left, right, down
    double *charge; // used for setup_zone, gather_charge
    int *bolt; // used for setup_zone, scatter_bolt
    double *boundary; // used for scatter_boundary

    int num_choice; // used for scatter_choice, gather_probs
    int *choice_idxs; // used for scatter_choice
    int *choice_idx_map; // used for gather_probs, map to g->choice_idxs's index
    double *probs; // used for gather_probs
    MPI_Request mpi_r;
}zonedef_t;

zonedef_t *generate_zones(graph_t *g, int process_count);
void send_zone(graph_t *g, zonedef_t *zonedef_list, int zone_id);
zone_t *setup_zone(int this_zone);

double get_charge(zone_t *z, int y, int x);
void exchange_charge(zone_t *z);
void scatter_power(int *power);
void scatter_bolt(int process_count, bool mpi_master, graph_t *g, zonedef_t *zlist, zone_t *z);
void scatter_choices(int process_count, bool mpi_master, graph_t *g, zonedef_t *zlist, zone_t *z);

void gather_probs(int process_count, bool mpi_master, graph_t *g, zonedef_t *zlist, zone_t *z);
void gather_charge(int process_count, bool mpi_master, graph_t *g, zonedef_t *zlist, zone_t *z);
void free_zonedef_list(zonedef_t *zonedef_list, int process_count);

#endif