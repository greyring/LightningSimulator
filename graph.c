#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "graph.h"

#define MAXLINE 1024

/**
 * store the whole graph and buffers
 */

// initialize buffer
static graph_t *new_graph(int width, int height, int power, int eta) {
    graph_t *g = malloc(sizeof(graph_t));
    int nnode = width * height;
    if (g == NULL)
        return NULL;
    g->width = width;
    g->height = height;
    g->power = power;
    g->eta = eta;
    g->charge = (double*)calloc(nnode, sizeof(double));
    g->charge_buffer = (double*)calloc(nnode, sizeof(double));
    g->boundary = (double*)calloc(nnode, sizeof(double));
    g->reset_bolt = (int*)calloc(nnode, sizeof(int));
    g->bolt = (int*)calloc(nnode, sizeof(int));
    g->num_choice = 0;
    g->choice_probs = (double*)calloc(nnode, sizeof(double));
    g->choice_idxs = (int*)calloc(nnode, sizeof(int));   
    g->choosed = (int*)calloc(nnode, sizeof(int));  
    g->path = (int*)calloc(nnode, sizeof(int));

    return g;
}

void free_graph(graph_t *g) {
    free(g->charge);
    free(g->charge_buffer);
    free(g->boundary);
    free(g->reset_bolt);
    free(g->bolt);
    free(g->choice_probs);
    free(g->choice_idxs);
    free(g->path);
    free(g);
}

/* Read in graph file and build graph data structure */
graph_t *read_graph(FILE *infile) {
    graph_t *g = NULL;
    char linebuf[MAXLINE];
    int width, height, power, eta;
    int num_negative, num_positive;
    int x, y;
    int i;

    // Read header information
    if (fgets(linebuf, MAXLINE, infile) == NULL) {
        return NULL;
    }
    if (sscanf(linebuf, "%d %d %d %d", &width, &height, &power, &eta) != 4) {
        fprintf(stderr, "Bad graph input Line 1\n");
        return NULL;
    }

    g = new_graph(width, height, power, eta);
    if (g == NULL) {
        fprintf(stderr, "Create graph failed\n");
	    return NULL;
    }

    // read positive bolts
    if (fgets(linebuf, MAXLINE, infile) == NULL) {
        return NULL;
    }
    if (sscanf(linebuf, "%d", &num_positive) != 1) {
        fprintf(stderr, "Bad graph input positive bolts\n");
        return NULL;
    }
    for (i = 0; i < num_positive; i++) {
        if (fgets(linebuf, MAXLINE, infile) == NULL) {
            return NULL;
        }
        if (sscanf(linebuf, "%d %d", &y, &x) != 2) {
            fprintf(stderr, "Bad graph input positive bolts\n");
            return NULL;
        }
        g->reset_bolt[y * g->width + x] = 1;
    }

    // read negative bolts
    if (fgets(linebuf, MAXLINE, infile) == NULL) {
        return NULL;
    }
    if (sscanf(linebuf, "%d", &num_negative) != 1) {
        fprintf(stderr, "Bad graph input negative bolts\n");
        return NULL;
    }
    for (i = 0; i < num_negative; i++) {
        if (fgets(linebuf, MAXLINE, infile) == NULL) {
            return NULL;
        }
        if (sscanf(linebuf, "%d %d", &y, &x) != 2) {
            fprintf(stderr, "Bad graph input negative bolts\n");
            return NULL;
        }
        g->reset_bolt[y * g->width + x] = -1;
    }

    return g;
}

/* print the bolt value to outfile */
void print_graph(graph_t *g, FILE *outfile) {
    int i, j;
    for (i = 0; i < g->height; i++) {
        for (j = 0; j < g->width; j++) {
            fprintf(outfile, "%d", g->bolt[i * g->width + j]);
            if (j != g->width) {
                fprintf(outfile, " ");
            }
        }
        fprintf(outfile, "\n");
    }
}

void print_charge(graph_t *g, FILE *outfile) {
    int i, j;
    for (i = 0; i < g->height; i++) {
        for (j = 0; j < g->width; j++) {
            fprintf(outfile, "%.2lf", g->charge[i * g->width + j]);
            if (j != g->width) {
                fprintf(outfile, " ");
            }
        }
        fprintf(outfile, "\n");
    }
}
