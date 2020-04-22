#include <stdio.h>
#include "graph.h"
#include "sim.h"

int main() {
    FILE *infile;
    graph_t *g;

    infile = fopen("data/small_3.graph", "r");
    if (infile == NULL) {
        fprintf(stderr, "Cannot open graph\n");
        return -1;
    }

    g = read_graph(infile);
    if (g == NULL) {
        return -1;
    }
    fclose(infile);

    simulate(g, 30);

    free_graph(g);

    return 0;
}

