#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>
#include "graph.h"
#include "mpiutil.h"
#include "instrument.h"

static zone_t *new_zone(int this_zone, int gheight, int gwidth, int start_row, int start_col, int height, int width, int eta) {
    zone_t *res = (zone_t*)calloc(1, sizeof(zone_t));
    res->this_zone = this_zone;
    res->gheight = gheight;
    res->gwidth = gwidth;
    res->start_row = start_row;
    res->start_col = start_col;
    res->height = height;
    res->width = width;
    res->eta = eta;

    res->charge = (double *)calloc(height * width, sizeof(double));
    res->charge_buffer = (double*)calloc(height * width, sizeof(double));
    res->ghost_charge = (double*)calloc(2 * (height + width), sizeof(double));
    res->boundary = (double*)calloc(height * width, sizeof(double));
    res->reset_bolt = (int*)calloc(height * width, sizeof(int));
    res->bolt = (int*)calloc(height * width, sizeof(int));
    
    res->left_buf = (double*)calloc(height, sizeof(double));
    res->right_buf = (double*)calloc(height, sizeof(double));
    res->choice_idxs = (int*)calloc(height * width, sizeof(int));
    res->prob_buf = (double*)calloc(height * width, sizeof(double));
    return res;
}

// only called by master
// divide graph into zones
zonedef_t *generate_zones(graph_t *g, int process_count) {
    zonedef_t *res;
    int i, j, idx;
    int num_row = (int)sqrt(process_count); // how many rows of zones
    int num_col = process_count / num_row; // how many columns of zones

    while (num_row * num_col != process_count) {
        num_row--;
        num_col = process_count / num_row;
    }

    res = (zonedef_t*) calloc(process_count, sizeof(zonedef_t));
    if (res == NULL) {
        return NULL;
    }

    // get boundary of each zone
    for (i = 0; i < num_row; i++) {
        for (j = 0; j < num_col; j++) {
            idx = i * num_col + j;
            res[idx].start_row = g->height / num_row * i;
            res[idx].start_col = g->width / num_col * j;
            if (i == num_row - 1) {
                res[idx].height = g->height - res[idx].start_row;
            } else {
                res[idx].height = g->height / num_row;
            }
            
            if (j == num_col - 1) {
                res[idx].width = g->width - res[idx].start_col;
            } else {
                res[idx].width = g->width / num_col;
            }

            res[idx].adj[0] = i > 0? (i - 1) * num_row + j : -1;
            res[idx].adj[1] = j > 0? i * num_row + j - 1 : -1;
            res[idx].adj[2] = j < num_col - 1? i * num_row + j + 1 : -1;
            res[idx].adj[3] = i < num_row - 1? (i + 1) * num_row + j : -1;
        }
    }

    // store data
    for (idx = 0; idx < process_count; idx++) {
        int start_row = res[idx].start_row;
        int start_col = res[idx].start_col;
        int height = res[idx].height;
        int width = res[idx].width;

        res[idx].eta = g->eta;
        res[idx].charge = calloc(width * height, sizeof(double));
        res[idx].boundary = calloc(width * height, sizeof(double));
        res[idx].bolt = calloc(width * height, sizeof(int));
        res[idx].choice_idxs = calloc(width * height, sizeof(int));
        res[idx].choice_idx_map = calloc(width * height, sizeof(int));
        res[idx].probs = calloc(width * height, sizeof(double));

        int b_idx = 0;
        int g_idx;
        for (i = start_row; i < start_row + height; i++) {
            for (j = start_col; j < start_col + width; j++) {
                g_idx = i * g->width + j;
                res[idx].charge[b_idx] = g->charge[g_idx];
                res[idx].bolt[b_idx] = g->reset_bolt[g_idx];
                b_idx++;
            }
        }
    }

    return res;
}

void free_zonedef_list(zonedef_t *zonedef_list, int process_count) {
    int i;
    for (i = 0; i < process_count; i++) {
        free(zonedef_list[i].charge);
        free(zonedef_list[i].bolt);
    }
    free(zonedef_list);
}

// only called by master
// send data to zones
void send_zone(graph_t *g, zonedef_t *zonedef_list, int zone_id) {
    int height = zonedef_list[zone_id].height;
    int width = zonedef_list[zone_id].width;
    MPI_Request dummy_request;
    
    MPI_Isend(&g->height, 1, MPI_INT, zone_id, 0, MPI_COMM_WORLD, &dummy_request);
    MPI_Isend(&g->width, 1, MPI_INT, zone_id, 0, MPI_COMM_WORLD, &dummy_request);
    MPI_Isend(&zonedef_list[zone_id].start_row, 1, MPI_INT, zone_id, 0, MPI_COMM_WORLD, &dummy_request);
    MPI_Isend(&zonedef_list[zone_id].start_col, 1, MPI_INT, zone_id, 0, MPI_COMM_WORLD, &dummy_request);
    MPI_Isend(&zonedef_list[zone_id].height, 1, MPI_INT, zone_id, 0, MPI_COMM_WORLD, &dummy_request);
    MPI_Isend(&zonedef_list[zone_id].width, 1, MPI_INT, zone_id, 0, MPI_COMM_WORLD, &dummy_request);
    MPI_Isend(&zonedef_list[zone_id].eta, 1, MPI_INT, zone_id, 0, MPI_COMM_WORLD, &dummy_request);
    MPI_Isend(zonedef_list[zone_id].adj, 4, MPI_INT, zone_id, 0, MPI_COMM_WORLD, &dummy_request);
    MPI_Isend(zonedef_list[zone_id].charge, width * height, MPI_DOUBLE, zone_id, 0, MPI_COMM_WORLD, &dummy_request);
    MPI_Isend(zonedef_list[zone_id].bolt, width * height, MPI_INT, zone_id, 0, MPI_COMM_WORLD, &dummy_request);
}

// called by all threads
// get data of the zone from master process
zone_t *setup_zone(int this_zone) {
    int gheight, gwidth;
    int start_row, start_col, height, width, eta;
    zone_t *zone = NULL;
    MPI_Recv(&gheight, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, NULL);
    MPI_Recv(&gwidth, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, NULL);
    MPI_Recv(&start_row, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, NULL);
    MPI_Recv(&start_col, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, NULL);
    MPI_Recv(&height, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, NULL);
    MPI_Recv(&width, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, NULL);
    MPI_Recv(&eta, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, NULL);
    zone = new_zone(this_zone, gheight, gwidth, start_row, start_col, height, width, eta);
    MPI_Recv(zone->adj, 4, MPI_INT, 0, 0, MPI_COMM_WORLD, NULL);
    MPI_Recv(zone->charge, height * width, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, NULL);
    MPI_Recv(zone->reset_bolt, height * width, MPI_INT, 0, 0, MPI_COMM_WORLD, NULL);

    // printf("%d: %d %d %d %d %d %d %d %d %d\n", this_zone, start_row, start_col, height, width, eta, zone->adj[0], zone->adj[1], zone->adj[2], zone->adj[3]);
    return zone;
}

// y and x are in the zone, [-1, height], [-1, width]
double get_charge(zone_t *z, int y, int x) {
    int height = z->height;
    int width = z->width;

    if (x >= 0 && x < width &&
        y >= 0 && y < height) {
        return z->charge[y * width + x];
    }

    if (y < 0) {
        if (z->adj[0] == -1)
            return 0;
        // ghost up row
        return z->ghost_charge[x];
    }

    if (x < 0) {
        if (z->adj[1] == -1)
            return 0;
        // ghost left col
        return z->ghost_charge[width + y];
    }

    if (x >= width) {
        if (z->adj[2] == -1)
            return 0;
        // ghost right col
        return z->ghost_charge[width + height + y];
    }

    if (y >= height) {
        if (z->adj[3] == -1)
            return 0;
        // ghost down row
        return z->ghost_charge[width + height + height + x];
    }

    // never reach here
    return 0.0;
}

// exchange charges between zones
void exchange_charge(zone_t *z) {
    int i;
    int width = z->width;
    int height = z->height;
    MPI_Request send_r[4];
    MPI_Request recv_r[4];

    START_ACTIVITY(ACTIVITY_COMM);
    // send charges
    if (z->adj[0] != -1) {
        // send to up
        MPI_Isend(z->charge, width, MPI_DOUBLE, z->adj[0], 1, MPI_COMM_WORLD, &send_r[0]);
    }
    if (z->adj[1] != -1) {
        // send to left
        for (i = 0; i < height; i++) {
            z->left_buf[i] = z->charge[i * width];
        }
        MPI_Isend(z->left_buf, height, MPI_DOUBLE, z->adj[1], 1, MPI_COMM_WORLD, &send_r[1]);
    }
    if (z->adj[2] != -1) {
        // send to right
        for (i = 0; i < height; i++) {
            z->right_buf[i] = z->charge[i * width + width - 1];
        }
        MPI_Isend(z->right_buf, height, MPI_DOUBLE, z->adj[2], 1, MPI_COMM_WORLD, &send_r[2]);
    }
    if (z->adj[3] != -1) {
        // send to down
        MPI_Isend(z->charge + (height - 1) * width, width, MPI_DOUBLE, z->adj[3], 1, MPI_COMM_WORLD, &send_r[3]);
    }

    // receive charges
    if (z->adj[0] != -1) {
        // receive up
        MPI_Irecv(z->ghost_charge, width, MPI_DOUBLE, z->adj[0], 1, MPI_COMM_WORLD, &recv_r[0]);
    }
    if (z->adj[1] != -1) {
        MPI_Irecv(z->ghost_charge + width, height, MPI_DOUBLE, z->adj[1], 1, MPI_COMM_WORLD, &recv_r[1]);
    }
    if (z->adj[2] != -1) {
        MPI_Irecv(z->ghost_charge + width + height, height, MPI_DOUBLE, z->adj[2], 1, MPI_COMM_WORLD, &recv_r[2]);
    }
    if (z->adj[3] != -1) {
        MPI_Irecv(z->ghost_charge + width + height + height, width, MPI_DOUBLE, z->adj[3], 1, MPI_COMM_WORLD, &recv_r[3]);
    }

    // wait send recv finish
    for (i = 0; i < 4; i++) {
        if (z->adj[i] != -1) {
            MPI_Wait(&send_r[i], MPI_STATUS_IGNORE);
            MPI_Wait(&recv_r[i], MPI_STATUS_IGNORE);
        }
    }
    FINISH_ACTIVITY(ACTIVITY_COMM);
}

void scatter_power(int *power) {
    START_ACTIVITY(ACTIVITY_COMM);
    MPI_Bcast(power, 1, MPI_INT, 0, MPI_COMM_WORLD);
    FINISH_ACTIVITY(ACTIVITY_COMM);
}

void scatter_bolt(int process_count, bool mpi_master, graph_t *g, zonedef_t *zlist, zone_t *z) {
    int i, j, idx, g_idx, b_idx;

    START_ACTIVITY(ACTIVITY_COMM);
    if (mpi_master) {
        for (idx = 0; idx < process_count; idx++) {
            int start_row = zlist[idx].start_row;
            int start_col = zlist[idx].start_col;
            int height = zlist[idx].height;
            int width = zlist[idx].width;
            b_idx = 0;
            for (i = start_row; i < start_row + height; i++) {
                for (j = start_col; j < start_col + width; j++) {
                    g_idx = i * g->width + j;
                    zlist[idx].bolt[b_idx] = g->bolt[g_idx];
                    b_idx++;
                }
            }
            MPI_Isend(zlist[idx].bolt, width * height, MPI_INT, idx, 3, MPI_COMM_WORLD, &zlist[idx].mpi_r);
        }
    }

    MPI_Recv(z->bolt, z->width * z->height, MPI_INT, 0, 3, MPI_COMM_WORLD, NULL);

    if (mpi_master) {
        for (idx = 0; idx < process_count; idx++) {
            MPI_Wait(&zlist[idx].mpi_r, MPI_STATUS_IGNORE);
        }
    }
    FINISH_ACTIVITY(ACTIVITY_COMM);
}

void scatter_choices(int process_count, bool mpi_master, graph_t *g, zonedef_t *zlist, zone_t *z) {
    int zid;
    MPI_Status status;

    START_ACTIVITY(ACTIVITY_COMM);
    if (mpi_master) {
        for (zid = 0; zid < process_count; zid++) {
            MPI_Isend(zlist[zid].choice_idxs, zlist[zid].num_choice, MPI_INT, zid, 3, MPI_COMM_WORLD, &zlist[zid].mpi_r);
        }
    }

    int width = z->width;
    int height = z->height;
    MPI_Recv(z->choice_idxs, width * height, MPI_INT, 0, 3, MPI_COMM_WORLD, &status);
    MPI_Get_count(&status, MPI_INT, &z->num_choice);

    if (mpi_master) {
        for (zid = 0; zid < process_count; zid++) {
            MPI_Wait(&zlist[zid].mpi_r, MPI_STATUS_IGNORE);
        }
    }
    FINISH_ACTIVITY(ACTIVITY_COMM);
}

void gather_probs(int process_count, bool mpi_master, graph_t *g, zonedef_t *zlist, zone_t *z) {
    int idx, b_idx, g_idx;

    START_ACTIVITY(ACTIVITY_COMM);
    MPI_Isend(z->prob_buf, z->num_choice, MPI_DOUBLE, 0, 3, MPI_COMM_WORLD, &z->mpi_r);

    if (mpi_master) {
        for (idx = 0; idx < process_count; idx++) {
            MPI_Irecv(zlist[idx].probs, zlist[idx].num_choice, MPI_DOUBLE, idx, 3, MPI_COMM_WORLD, &zlist[idx].mpi_r);
        }

        for (idx = 0; idx < process_count; idx++) {
            MPI_Wait(&zlist[idx].mpi_r, MPI_STATUS_IGNORE);
            for (b_idx = 0; b_idx < zlist[idx].num_choice; b_idx++) {
                g_idx = zlist[idx].choice_idx_map[b_idx];
                g->choice_probs[g_idx] = zlist[idx].probs[b_idx];
            }
        }
    }

    MPI_Wait(&z->mpi_r, MPI_STATUS_IGNORE);
    FINISH_ACTIVITY(ACTIVITY_COMM);
}

void gather_charge(int process_count, bool mpi_master, graph_t *g, zonedef_t *zlist, zone_t *z) {
    MPI_Request r;
    int b_idx, g_idx, idx;

    START_ACTIVITY(ACTIVITY_COMM);
    // send charge to master
    MPI_Isend(z->charge, z->width * z->height, MPI_DOUBLE, 0, 3, MPI_COMM_WORLD, &r);
    if (mpi_master) {
        // gather charge from all zones
        for (idx = 0; idx < process_count; idx++) {
            MPI_Irecv(zlist[idx].charge, zlist[idx].width * zlist[idx].height, MPI_DOUBLE, idx, 3, MPI_COMM_WORLD, &zlist[idx].mpi_r);
        }

        // store data in zlist's buffer and set charge in graph
        for (idx = 0; idx < process_count; idx++) {
            int start_row = zlist[idx].start_row;
            int start_col = zlist[idx].start_col;
            int width = zlist[idx].width;
            int height = zlist[idx].height;
            MPI_Wait(&zlist[idx].mpi_r, MPI_STATUS_IGNORE);
            for (b_idx = 0; b_idx < height * width; b_idx++) {
                g_idx = ((b_idx / width) + start_row) * g->width + (b_idx % width) + start_col;
                g->charge[g_idx] = zlist[idx].charge[b_idx];
            }
        }
    }

    MPI_Wait(&r, MPI_STATUS_IGNORE);
    FINISH_ACTIVITY(ACTIVITY_COMM);
}