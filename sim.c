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
    int i, j;
    int idx;
    double sum;
    START_ACTIVITY(ACTIVITY_UPDATE);
    for (i = 0; i < g->height; i++) {
        for (j = 0; j < g->width; j++) {
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
    for (i = 0; i < g->height * g->width; i++) {
        g->charge[i] = g->charge_buffer[i];
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

static int find_next(graph_t *g) {
    double prob, breach;
    int num_choice;
    int adj, idx, choice;
    int i, j;

    START_ACTIVITY(ACTIVITY_NEXT);
    num_choice = 0;
    for (i = 0; i < g->height; i++) {
        for (j = 0; j < g->width; j++) {
            idx = i * g->width + j;
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

    breach = (double)rand()/RAND_MAX * g->choice_probs[num_choice - 1];
    choice = locate_value(breach, g->choice_probs, num_choice);
    FINISH_ACTIVITY(ACTIVITY_NEXT);

    if (choice == -1)
        return -1;
    return g->choice_idxs[choice];
}

static void simulate_one(graph_t *g) {
    int power = g->power;
    int next_bolt = -1;

    while (power > 0) {
        update_charge(g);
        next_bolt = find_next(g);
        if (next_bolt != -1) {
            g->path[next_bolt] = adjacent_pos(g, next_bolt / g->width, next_bolt % g->width);
            if (g->bolt[next_bolt] < 0) {
                power += g->bolt[next_bolt];
                discharge(g, next_bolt, -g->bolt[next_bolt]);
            }
            g->bolt[next_bolt] = 1;
        }
    }
}

void simulate(graph_t *g, int count, FILE *ofile) {
    int i, idx;

    for (i = 0; i < g->width + g->height; i++) {
        update_charge(g);
    }

    // generate lightnings
    for (i = 0; i < count; i++) {
        simulate_one(g);

        START_ACTIVITY(ACTIVITY_RECOVER);
        // one lightning is generated
        for (idx = 0; idx < g->height * g->width; idx++) {
            if (g->bolt[idx] > 1) {
                g->boundary[idx] = g->bolt[idx] * 0.0001;
            } else {
                g->boundary[idx] = 0;
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
