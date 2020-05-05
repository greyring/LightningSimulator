#include <stdbool.h>
#include <math.h>
#include <stdlib.h>
#include <mpi.h>
#include "graph.h"
#include "mpiutil.h"
#include "sim-mpi.h"
#include "instrument.h"

/* What is the crossover between binary and linear search */
#define BINARY_THRESHOLD 4

/*
  Linear search
 */
static inline int locate_value_linear(double target, double *list, int len) {
    int i;
    for (i = 0; i < len; i++)
	    if (target < list[i])
	        return i;
    /* Shouldn't get here */
    return -1;
}

/*
  Binary search down to threshold, and then linear
 */
static inline int locate_value(double target, double *list, int len) {
    int left = 0;
    int right = len-1;
    while (left < right) {
	    if (right-left+1 < BINARY_THRESHOLD)
	        return left + locate_value_linear(target, list+left, right-left+1);
	    int mid = left + (right-left)/2;
	    if (target < list[mid])
	        right = mid;
	    else
	        left = mid+1;
    }
    return right;
}

static void reset_charge(zone_t *z) {
    int i;
    for (i = 0; i < z->height * z->width; i++) {
        z->charge[i] = z->charge_buffer[i] = 0;
    }
}

static void reset_boundary(zone_t *z) {
    int i;
    for (i = 0; i < z->height * z->width; i++) {
        z->boundary[i] = 0.0;
    }
}

static void reset_bolt(graph_t *g) {
    int i;
    for (i = 0; i < g->height * g->width; i++) {
        g->bolt[i] = g->reset_bolt[i];
    }
}

static void update_boundary(zone_t *z) {
    int idx;
    for (idx = 0; idx < z->height * z->width; idx++) {
        if (z->bolt[idx] > 1) {
            z->boundary[idx] = z->bolt[idx] * 0.0001;
        } else {
            z->boundary[idx] = 0;
        }
    }
}

static void reset_path(graph_t *g) {
    int i;
    for (i = 0; i < g->height * g->width; i++) {
        g->path[i] = -1;
    }
}

static void choose_helper(int process_count, graph_t *g, zonedef_t *zlist, int bolt_idx, int i, int j) {
    int idx = i * g->width + j;
    int zid;
    if (i >= 0 && i < g->height && j >= 0 && j < g->width &&
        g->choosed[idx] == 0 && g->bolt[idx] <= 0) {
        g->choosed[idx] = 1;
        g->choice_idxs[g->num_choice] = idx;
        for (zid = 0; zid < process_count; zid++) {
            if (i >= zlist[zid].start_row && i < zlist[zid].start_row + zlist[zid].height &&
            j >= zlist[zid].start_col && j < zlist[zid].start_col + zlist[zid].width) {
                zlist[zid].choice_idxs[zlist[zid].num_choice] = (i - zlist[zid].start_row) * zlist[zid].width + j - zlist[zid].start_col;
                zlist[zid].choice_idx_map[zlist[zid].num_choice] = g->num_choice;
                zlist[zid].num_choice++;
                break;
            }
        }
        g->num_choice++;
        g->path[idx] = bolt_idx;
    }
}

static void find_choice(int process_count, graph_t *g, zonedef_t *zlist, int idx) {
    int i, j;
    if (g->bolt[idx] > 0) {
        i = idx / g->width;
        j = idx % g->width;
        choose_helper(process_count, g, zlist, idx, i - 1, j);
        choose_helper(process_count, g, zlist, idx, i, j - 1);
        choose_helper(process_count, g, zlist, idx, i, j + 1);
        choose_helper(process_count, g, zlist, idx, i + 1, j);
    }
}

static void reset_choice(int process_count, graph_t *g, zonedef_t *zlist) {
    int i;
    g->num_choice = 0;
    for (i = 0; i < g->height * g->width; i++) {
        g->choosed[i] = 0;
    }
    for (i = 0; i < process_count; i++) {
        zlist[i].num_choice = 0;
    }
    // get choices idxs
    for (i = 0; i < g->width * g->height; i++) {
        find_choice(process_count, g, zlist, i);
    }
}
// get bolt at x, y
// if bolt < 0.0, charge = 1.0 // boundary
// if bolt > 0.0, charge = 0.0 // boundary
// else charge = (boundary + neighbor's charge) / 4
static void update_charge(zone_t *z) {
    int height = z->height;
    int width = z->width;
    int i, j, idx;
    double sum;

    exchange_charge(z);

    START_ACTIVITY(ACTIVITY_UPDATE);
    for (idx = 0; idx < height * width; idx++) {
        i = idx / width;
        j = idx % width;

        // boundary condition
        if (z->bolt[idx] < 0) {
            z->charge_buffer[idx] = 1;
        } else if (z->bolt[idx] > 0) {
            z->charge_buffer[idx] = 0;
        } else {
            sum = z->boundary[idx]; // poisson equation
            sum += get_charge(z, i - 1, j);
            sum += get_charge(z, i, j - 1);
            sum += get_charge(z, i, j + 1);
            sum += get_charge(z, i + 1, j);
            z->charge_buffer[idx] = sum / 4;
        }
    }

    // replace origin
    for (idx = 0; idx < height * width; idx++) {
        z->charge[idx] = z->charge_buffer[idx];
    }
    FINISH_ACTIVITY(ACTIVITY_UPDATE);
}

// add charge to bolt along the path
static void discharge(graph_t *g, int index, int charge) {
    int count = 500;
    while (index != -1 && count > 0) {
        count -= 1;
        g->bolt[index] += charge;
        index = g->path[index];
    }
}

static void calc_prob(zone_t *z) {
    int i, idx;

    START_ACTIVITY(ACTIVITY_NEXT);
    // calculate probability based on latest charge
    for (i = 0; i < z->num_choice; i++) {
        idx = z->choice_idxs[i];
        if (z->bolt[idx] > 0) {
            z->prob_buf[i] = 0;
        } else {
            z->prob_buf[i] = pow(z->charge[idx], z->eta);
        }
    }
    FINISH_ACTIVITY(ACTIVITY_NEXT);
}

static int find_next(graph_t *g) {
    int i, choice;
    double breach;
    for (i = 1; i < g->num_choice; i++) {
        g->choice_probs[i] += g->choice_probs[i - 1];
    }
    breach = (double)rand()/RAND_MAX * g->choice_probs[g->num_choice - 1];
    choice = locate_value(breach, g->choice_probs, g->num_choice);
    if (choice == -1)
        return -1;
    return g->choice_idxs[choice];
}

static void simulate_one(int process_count, bool mpi_master, graph_t *g, zonedef_t *zlist, zone_t *z) {
    int power;

    START_ACTIVITY(ACTIVITY_RECOVER);
    if (mpi_master) {
        power = g->power;
        reset_bolt(g);
        reset_path(g);
        reset_choice(process_count, g, zlist);
    }
    FINISH_ACTIVITY(ACTIVITY_RECOVER);

    scatter_power(&power);

    while (power > 0) {
        scatter_bolt(process_count, mpi_master, g, zlist, z);
        update_charge(z);

        scatter_choices(process_count, mpi_master, g, zlist, z);
        calc_prob(z);
        gather_probs(process_count, mpi_master, g, zlist, z);

        if (mpi_master) {
            START_ACTIVITY(ACTIVITY_NEXT);
            int next_bolt = -1;
            next_bolt = find_next(g);
            if (next_bolt != -1) {
                if (g->bolt[next_bolt] < 0) {
                    power += g->bolt[next_bolt];
                    discharge(g, next_bolt, -g->bolt[next_bolt]);
                }
                g->bolt[next_bolt] = 1;
                find_choice(process_count, g, zlist, next_bolt);
            }
            FINISH_ACTIVITY(ACTIVITY_NEXT);
        }
        scatter_power(&power);
    }

    // one lightning is generated
    scatter_bolt(process_count, mpi_master, g, zlist, z);

    START_ACTIVITY(ACTIVITY_RECOVER);
    update_boundary(z);
    FINISH_ACTIVITY(ACTIVITY_RECOVER);
}

void simulate(int process_count, bool mpi_master, graph_t *g, zonedef_t *zlist, zone_t *z, int count, FILE *ofile) {
    int i;

    // init graph
    if (mpi_master) {
        reset_bolt(g);
    }
    reset_charge(z);
    reset_boundary(z);
    scatter_bolt(process_count, mpi_master, g, zlist, z);

    for (i = 0; i < z->gheight + z->gwidth; i++) {
        update_charge(z);
    }

    // generate lightnings
    for (i = 0; i < count; i++) {
        simulate_one(process_count, mpi_master, g, zlist, z);

        START_ACTIVITY(ACTIVITY_PRINT);
        // print bolt
        if (mpi_master) {
            print_graph(g, ofile);
            fprintf(ofile, "\n");
        }
        FINISH_ACTIVITY(ACTIVITY_PRINT);
    }
}
