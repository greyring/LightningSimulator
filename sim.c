#include <math.h>
#include <stdlib.h>
#include "graph.h"
#include "sim.h"

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
    double total, prob;
    int num_choice;
    int adj, idx;
    
    while (power > 0) {
        update_charge(g);
        total = 0.0;
        num_choice = 0;
        for (int i = 0; i < g->height; i++) {
            for (int j = 0; j < g->width; j++) {
                if ((adj = adjacent_pos(g, i, j)) != -1) {
                    idx = i * g->width + j;
                    if (g->charge[idx] == 0)
                        continue;
                    g->path[idx] = adj;
                    prob = pow(g->charge[idx], g->eta);
                    g->choice_probs[num_choice] = prob;
                    g->choice_idxs[num_choice] = idx;
                    num_choice++;
                    total += prob;
                }
            }
        }

        if (num_choice != 0) {
            double breach = (double)rand()/RAND_MAX * total;
            int choice = locate_value(breach, g->choice_probs, num_choice);
            int choice_idx = g->choice_idxs[choice];
            if (g->bolt[choice_idx] < 0) {
                power += g->bolt[choice_idx];
                discharge(g, choice_idx, -g->bolt[choice_idx]);
            }
            g->bolt[choice_idx] = 1;
        }
    }
}

void simulate(graph_t *g, int count) {
    int idx;

    // init graph
    init_charge(g);
    init_boundary(g);
    srand(1);

    // generate lightnings
    for (int i = 0; i < count; i++) {
        init_bolt(g);
        init_path(g);
        simulate_one(g);
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

        // print bolt
        print_graph(g, stdout);
        fprintf(stdout, "\n");
    }
}
