#include <math.h>
#include <stdlib.h>
#include <omp.h>
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
    int g_width = g->width;
    int g_height = g->height;
    #pragma omp master
    {
        START_ACTIVITY(ACTIVITY_UPDATE);
    }
    #pragma omp for schedule(dynamic,32)
    for(idx = 0; idx < g_width*g_height; idx++){
        int i = idx / g_width;
        int j = idx % g_width;
        if (g->bolt[idx] < 0) {
            g->charge_buffer[idx] = 1.0;
        } else if (g->bolt[idx] > 0) {
            g->charge_buffer[idx] = 0.0;
        } else {
            double sum = g->boundary[idx]; // poisson equation
            if (i > 0)
                sum += g->charge[(i - 1) * g_width + j];
            if (i < g_height - 1)
                sum += g->charge[(i + 1) * g_width+ j];
            if (j > 0)
                sum += g->charge[i * g_width + j - 1];
            if (j < g_width - 1)
                sum += g->charge[i * g_width + j + 1];
            g->charge_buffer[idx] = sum / 4;
        }
    }

    // replace origin
    #pragma omp for
    for (idx = 0; idx < g->height * g->width; idx++) {
        g->charge[idx] = g->charge_buffer[idx];
    }
    #pragma omp master
    {
        FINISH_ACTIVITY(ACTIVITY_UPDATE);
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

static void find_next(graph_t *g, int* choice_point) {
    double breach, prob;
    int i, idx, choice;
    int g_eta = g->eta;

    #pragma omp for 
    for (i = 0; i < g->num_choice; i++) {
        idx = g->choice_idxs[i];

        if (g->bolt[idx] > 0) {
            prob = 0;
        } else {
            prob = pow(g->charge[idx], g_eta);
        }
        g->choice_probs[i] = prob;
    }
    
    #pragma omp master
    {
        
        for (i = 1; i < g->num_choice; i++) {
            g->choice_probs[i] += g->choice_probs[i - 1];
        }
        breach = (double)rand()/RAND_MAX * g->choice_probs[g->num_choice - 1];
        choice = locate_value(breach, g->choice_probs, g->num_choice); 
        if (choice == -1)
            *choice_point = -1;
        *choice_point = g->choice_idxs[choice];
    }
    #pragma omp barrier
}

static void simulate_one(graph_t *g, int *g_power, int *g_num_choice) {

    int idx;
    
    #pragma omp master
    {
        START_ACTIVITY(ACTIVITY_RECOVER);
        reset_bolt(g);
        reset_path(g);
        reset_choice(g);
        FINISH_ACTIVITY(ACTIVITY_RECOVER);
    }
    #pragma omp barrier
    while (*g_power > 0) {
        update_charge(g);
        int next_bolt = -1;
        #pragma omp master
        {
            START_ACTIVITY(ACTIVITY_NEXT);
        }
        find_next(g, &next_bolt);

        #pragma omp master
        {
            // printf("next_bolt: %d\n", next_bolt);   
            if (next_bolt != -1) {
                if (g->bolt[next_bolt] < 0) {
                    *g_power += g->bolt[next_bolt];
                    discharge(g, next_bolt, -g->bolt[next_bolt]);
                }
                g->bolt[next_bolt] = 1;
                find_choice(g, next_bolt);
            }
            FINISH_ACTIVITY(ACTIVITY_NEXT);

        }
        #pragma omp barrier
    }
    #pragma omp master
    {
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
    }
    #pragma omp barrier
}

void simulate(graph_t *g, int count, FILE *ofile) {
    int g_power;
    int g_num_choice;
    #pragma omp parallel
    {
        int i;
        for (i = 0; i < g->width + g->height; i++) {
            update_charge(g);
        }
        // generate lightnings
        for (i = 0; i < count; i++) {
            g_power = g->power;
            simulate_one(g, &g_power, &g_num_choice);
            #pragma omp master
            {
                // print bolt
                START_ACTIVITY(ACTIVITY_PRINT);
                print_graph(g, ofile);
                fprintf(ofile, "\n");
                FINISH_ACTIVITY(ACTIVITY_PRINT);

            }
            #pragma omp barrier
        }
    }
}
