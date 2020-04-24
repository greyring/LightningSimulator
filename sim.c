#include <math.h>
#include <stdlib.h>
#include "graph.h"
#include "sim.h"
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

// get bolt at x, y
// if bolt < 0.0, charge = 1.0 // boundary
// if bolt > 0.0, charge = 0.0 // boundary
// else charge = (boundary + neighbor's charge) / 4
static void update_charge(graph_t *g) {
    int idx;
    double sum;
    for (int i = 0; i < g->height; i++) {
        for (int j = 0; j < g->width; j++) {
            idx = i * g->width + j;
            // boundary condition
            if (g->bolt[idx] < 0) {
                g->charge_buffer[idx] = 1.0;
            } else if (g->bolt[idx] > 0) {
                g->charge_buffer[idx] = 0.0;
            } else {
                sum = g->boundary[idx]; // poisson equation
                if (i > 0)
                    sum += g->charge[(i - 1) * g->width + j];
                if (i < g->height - 1)
                    sum += g->charge[(i + 1) * g->width + j];
                if (j > 0)
                    sum += g->charge[i * g->width + j - 1];
                if (j < g->width - 1)
                    sum += g->charge[i * g->width + j + 1];
                g->charge_buffer[idx] = sum / 4;
            }
        }
    }

    // replace origin
    for (int i = 0; i < g->height; i++) {
        for (int j = 0; j < g->width; j++) {
            idx = i * g->width + j;
            g->charge[idx] = g->charge_buffer[idx];
        }
    }
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

static void simulate_one(graph_t *g) {
    int power = g->power;
    double prob;
    int num_choice;
    int adj, idx;
    
    while (power > 0) {
        START_ACTIVITY(ACTIVITY_UPDATE);
        update_charge(g);
        FINISH_ACTIVITY(ACTIVITY_UPDATE);

        START_ACTIVITY(ACTIVITY_NEXT);
        num_choice = 0;
        for (int i = 0; i < g->height; i++) {
            for (int j = 0; j < g->width; j++) {
                idx = i * g->width + j;
                // if (g->charge[idx] == 0)
                //         continue;
                if ((adj = adjacent_pos(g, i, j)) != -1) {
                    prob = pow(g->charge[idx], g->eta);
                    if (num_choice == 0)
                        g->choice_probs[num_choice] = prob;
                    else
                        g->choice_probs[num_choice] = g->choice_probs[num_choice - 1] + prob;
                    g->choice_idxs[num_choice] = idx;
                    num_choice++;
                }
            }
        }

        if (num_choice != 0) {
            double breach = (double)rand()/RAND_MAX * g->choice_probs[num_choice - 1];
            int choice = locate_value(breach, g->choice_probs, num_choice);
            int choice_idx = g->choice_idxs[choice];
            int adj = adjacent_pos(g, choice_idx / g->width, choice_idx % g->width);
            g->path[choice_idx] = adj;
            if (g->bolt[choice_idx] < 0) {
                power += g->bolt[choice_idx];
                discharge(g, choice_idx, -g->bolt[choice_idx]);
            }
            g->bolt[choice_idx] = 1;
        }
        FINISH_ACTIVITY(ACTIVITY_NEXT);
    }
}

void simulate(graph_t *g, int count, FILE *ofile) {
    int idx;

    // init graph
    reset_charge(g);
    reset_boundary(g);
    reset_bolt(g);
    reset_path(g);

    // generate lightnings
    for (int i = 0; i < count; i++) {
        simulate_one(g);

        START_ACTIVITY(ACTIVITY_RECOVER);
        // one lightning is generated
        for (int y = 0; y < g->height; y++) {
            for (int x = 0; x < g->width; x++) {
                idx = y * g->width + x;
                if (g->bolt[idx] > 1) {
                    g->boundary[idx] = g->bolt[idx] * 0.001;
                } else {
                    g->boundary[idx] = 0;
                }
            }
        }
        FINISH_ACTIVITY(ACTIVITY_RECOVER);

        START_ACTIVITY(ACTIVITY_PRINT);
        // print bolt
        print_graph(g, ofile);
        fprintf(ofile, "\n");
        FINISH_ACTIVITY(ACTIVITY_PRINT);

        START_ACTIVITY(ACTIVITY_RECOVER);
        reset_bolt(g);
        reset_path(g);
        FINISH_ACTIVITY(ACTIVITY_RECOVER);
    }
}
