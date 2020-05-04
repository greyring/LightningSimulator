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

static void reset_charge(graph_t *g) {
    int i;
    for (i = 0; i < g->height * g->width; i++) {
        g->charge[i] = g->charge_buffer[i] = 0;
    }
}

static void reset_boundary(graph_t *g) {
    int i;
    for (i = 0; i < g->height * g->width; i++) {
        g->boundary[i] = 0.0;
    }
}

static void reset_bolt(graph_t *g) {
    int i;
    for (i = 0; i < g->height * g->width; i++) {
        g->bolt[i] = g->reset_bolt[i];
    }
}

static void reset_path(graph_t *g) {
    int i;
    for (i = 0; i < g->height * g->width; i++) {
        g->path[i] = -1;
    }
}

static void choose_helper(graph_t *g, int bolt_idx, int i, int j) {
    int idx = i * g->width + j;
    if (i >= 0 && i < g->height && j >= 0 && j < g->width &&
        g->choosed[idx] == 0 && g->bolt[idx] <= 0) {
        g->choosed[idx] = 1;
        g->choice_idxs[g->num_choice] = idx;
        g->num_choice++;
        g->path[idx] = bolt_idx;
    }
}

static void find_choice(graph_t *g, int idx) {
    int i, j;
    if (g->bolt[idx] > 0) {
        i = idx / g->width;
        j = idx % g->width;
        choose_helper(g, idx, i - 1, j);
        choose_helper(g, idx, i, j - 1);
        choose_helper(g, idx, i, j + 1);
        choose_helper(g, idx, i + 1, j);
    }
}

static void reset_choice(graph_t *g) {
    int i;
    g->num_choice = 0;
    for (i = 0; i < g->height * g->width; i++) {
        g->choosed[i] = 0;
    }
    // get choices idxs
    for (i = 0; i < g->width * g->height; i++) {
        find_choice(g, i);
    }
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
    for (idx = 0; idx < g->height * g->width; idx++) {
        i = idx / g->width;
        j = idx % g->width;

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

    // replace origin
    for (idx = 0; idx < g->height * g->width; idx++) {
        g->charge[idx] = g->charge_buffer[idx];
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
    int i, idx, choice;

    // calculate probability based on latest charge
    for (i = 0; i < g->num_choice; i++) {
        idx = g->choice_idxs[i];

        if (g->bolt[idx] > 0) {
            prob = 0;
        } else {
            prob = pow(g->charge[idx], g->eta);
        }
    
        if (i == 0) {
            g->choice_probs[i] = prob;
        } else {
            g->choice_probs[i] = g->choice_probs[i - 1] + prob;
        }
    }

    // choose one as bolt
    breach = (double)rand()/RAND_MAX * g->choice_probs[g->num_choice - 1];
    choice = locate_value(breach, g->choice_probs, g->num_choice);

    if (choice == -1)
        return -1;
    return g->choice_idxs[choice];
}

static void simulate_one(graph_t *g) {
    int power = g->power;
    int next_bolt = -1;
    int idx;

    START_ACTIVITY(ACTIVITY_RECOVER);
    reset_bolt(g);
    reset_path(g);
    reset_choice(g);
    FINISH_ACTIVITY(ACTIVITY_RECOVER);

    while (power > 0) {
        update_charge(g);

        START_ACTIVITY(ACTIVITY_NEXT);
        next_bolt = find_next(g);
        if (next_bolt != -1) {
            if (g->bolt[next_bolt] < 0) {
                power += g->bolt[next_bolt];
                discharge(g, next_bolt, -g->bolt[next_bolt]);
            }
            g->bolt[next_bolt] = 1;
            find_choice(g, next_bolt);
        }
        FINISH_ACTIVITY(ACTIVITY_NEXT);
    }

    // one lightning is generated
    START_ACTIVITY(ACTIVITY_RECOVER);
    for (idx = 0; idx < g->height * g->width; idx++) {
        if (g->bolt[idx] > 1) {
            g->boundary[idx] = g->bolt[idx] * 0.0001;
        } else {
            g->boundary[idx] = 0;
        }
    }
    FINISH_ACTIVITY(ACTIVITY_RECOVER);
}

void simulate(graph_t *g, int count, FILE *ofile) {
    int i;

    // init graph
    reset_bolt(g);
    reset_charge(g);
    reset_boundary(g);

    for (i = 0; i < g->width + g->height; i++) {
        update_charge(g);
    }

    // generate lightnings
    for (i = 0; i < count; i++) {
        simulate_one(g);

        START_ACTIVITY(ACTIVITY_PRINT);
        // print bolt
        print_graph(g, ofile);
        fprintf(ofile, "\n");
        FINISH_ACTIVITY(ACTIVITY_PRINT);
    }
}
