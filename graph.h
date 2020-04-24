#ifndef __GRAPH_H__
#define __GRAPH_H__
#include <stdio.h>

typedef struct {
    int width;
    int height;
    int power; // #branchs of lightning
    int eta; // shape of lightning

    // electrical potential
    double *charge;
    double *charge_buffer;

    // charge density (poisson equation)
    double *boundary;

    int *reset_bolt;
    int *bolt;
    
    double *choice_probs;
    int *choice_idxs;
    int *path;
}graph_t;

graph_t *read_graph(FILE *infile);
void reset_charge(graph_t *g);
void reset_boundary(graph_t *g);
void reset_bolt(graph_t *g);
void reset_path(graph_t *g);
void free_graph(graph_t *g);
int adjacent_pos(graph_t *g, int y, int x);
void print_graph(graph_t *g, FILE *outfile);
void print_charge(graph_t *g, FILE *outfile) ;

#endif