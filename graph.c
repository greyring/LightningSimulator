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
    g->init_bolt = (int*)calloc(nnode, sizeof(int));
    g->bolt = (int*)calloc(nnode, sizeof(int));
    g->choice_probs = (double*)calloc(nnode, sizeof(double));
    g->choice_idxs = (int*)calloc(nnode, sizeof(int));    
    g->path = (int*)calloc(nnode, sizeof(int));

    return g;
}

void free_graph(graph_t *g) {
    free(g->charge);
    free(g->charge_buffer);
    free(g->boundary);
    free(g->init_bolt);
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
    for (int i = 0; i < num_positive; i++) {
        if (fgets(linebuf, MAXLINE, infile) == NULL) {
            return NULL;
        }
        if (sscanf(linebuf, "%d %d", &y, &x) != 2) {
            fprintf(stderr, "Bad graph input positive bolts\n");
            return NULL;
        }
        g->init_bolt[y * g->width + x] = 1;
    }

    // read negative bolts
    if (fgets(linebuf, MAXLINE, infile) == NULL) {
        return NULL;
    }
    if (sscanf(linebuf, "%d", &num_negative) != 1) {
        fprintf(stderr, "Bad graph input negative bolts\n");
        return NULL;
    }
    for (int i = 0; i < num_negative; i++) {
        if (fgets(linebuf, MAXLINE, infile) == NULL) {
            return NULL;
        }
        if (sscanf(linebuf, "%d %d", &y, &x) != 2) {
            fprintf(stderr, "Bad graph input negative bolts\n");
            return NULL;
        }
        g->init_bolt[y * g->width + x] = -1;
    }

    return g;
}

void init_charge(graph_t *g) {
    for (int i = 0; i < g->height; i++) {
        for (int j = 0; j < g->width; j++) {
            g->charge[i * g->width + j] = 0;
            g->charge_buffer[i * g->width + j] = 0;
        }
    }
}

void init_boundary(graph_t *g) {
    for (int i = 0; i < g->height; i++) {
        for (int j = 0; j < g->width; j++) {
            g->boundary[i * g->width + j] = 0.0;
        }
    }
}

void init_bolt(graph_t *g) {
    for (int i = 0; i < g->height; i++) {
        for (int j = 0; j < g->width; j++) {
            g->bolt[i * g->width + j] = g->init_bolt[i * g->width + j];
        }
    }
}

void init_path(graph_t *g) {
    for (int i = 0; i < g->height; i++) {
        for (int j = 0; j < g->width; j++) {
            g->path[i * g->width + j] = -1;
        }
    }
}

// return index of adjacent node with bolt > 0, else return -1
// up, left, right, down
int adjacent_pos(graph_t *g, int y, int x) {
    if (g->bolt[y * g->width + x] > 0)
        return -1;
    if (y > 0 && g->bolt[(y - 1) * g->width + x] > 0)
        return (y - 1) * g->width + x;
    if (x > 0 && g->bolt[y * g->width + x - 1] > 0)
        return y * g->width + x - 1;
    if (x < g->width - 1 && g->bolt[y * g->width + x + 1] > 0)
        return y * g->width + x + 1;
    if (y < g->height - 1 && g->bolt[(y + 1) * g->width + x] > 0)
        return (y + 1) * g->width + x;
    return -1;
}

/* print the bolt value to outfile */
void print_graph(graph_t *g, FILE *outfile) {
    for (int i = 0; i < g->height; i++) {
        for (int j = 0; j < g->width; j++) {
            fprintf(outfile, "%d", g->bolt[i * g->width + j]);
            if (j != g->width) {
                fprintf(outfile, " ");
            }
        }
        fprintf(outfile, "\n");
    }
}

void print_charge(graph_t *g, FILE *outfile) {
    for (int i = 0; i < g->height; i++) {
        for (int j = 0; j < g->width; j++) {
            fprintf(outfile, "%.2lf", g->charge[i * g->width + j]);
            if (j != g->width) {
                fprintf(outfile, " ");
            }
        }
        fprintf(outfile, "\n");
    }
}
