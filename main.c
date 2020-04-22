#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>
#include <stdbool.h>
#include "graph.h"
#include "sim.h"

static void usage(char *name) {
    char *use_string = "-g GFILE [-n STEPS] [-s SEED] [-u (r|b|s)] [-q] [-t THD] [-I]";
    fprintf(stdout, "Usage: %s %s\n", name, use_string);
    fprintf(stdout, "   -h        Print this message\n");
    fprintf(stdout, "   -g GFILE  Graph file\n");
    fprintf(stdout, "   -n STEPS  Number of simulation steps\n");
    fprintf(stdout, "   -s SEED   Initial RNG seed\n");
    fprintf(stdout, "   -t THD    Set number of threads\n");
    fprintf(stdout, "   -I        Instrument simulation activities\n");
    exit(0);
}

int main(int argc, char *argv[]) {
    FILE *gfile = NULL;
    graph_t *g = NULL;
    int count = 10;
    int thread_count = 1;
    unsigned long seed = 1;
    bool instrument = false;

    char c;
    char *optstring = "hg:n:s:t:I";
    while ((c = getopt(argc, argv, optstring)) != -1) {
        switch(c) {
        case 'h':
            usage(argv[0]);
            break;
        case 'g':
            gfile = fopen(optarg, "r");
            break;
        case 'n':
            count = atoi(optarg);
            break;
        case 's':
            seed = strtoul(optarg, NULL, 0);
            break;
        case 't':
            thread_count = atoi(optarg);
            break;
        case 'I':
            instrument = true;
            break;
        default:
            fprintf(stdout, "Unknown option '%c'\n", c);
            usage(argv[0]);
            exit(1);
        }
    }

    if (gfile == NULL) {
        fprintf(stdout, "Couldn't open graph file\n");
        exit(1);
    }
    g = read_graph(gfile);
    if (g == NULL) {
        exit(1);
    }
    fclose(gfile);

    srand(seed);

    simulate(g, count);

    free_graph(g);

    return 0;
}

