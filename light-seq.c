#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>
#include <stdbool.h>
#include "graph.h"
#include "sim.h"
#include "instrument.h"

static void usage(char *name) {
    char *use_string = "-g GFILE [-n STEPS] [-s SEED] [-u (r|b|s)] [-q] [-I]";
    fprintf(stdout, "Usage: %s %s\n", name, use_string);
    fprintf(stdout, "   -h        Print this message\n");
    fprintf(stdout, "   -g GFILE  Graph file\n");
    fprintf(stdout, "   -n STEPS  Number of simulation steps\n");
    fprintf(stdout, "   -s SEED   Initial RNG seed\n");
    fprintf(stdout, "   -I        Instrument simulation activities\n");
    exit(0);
}

int main(int argc, char *argv[]) {
    FILE *gfile = NULL;
    FILE *ofile = stdout;
    graph_t *g = NULL;
    int count = 10;
    unsigned long seed = 1;
    bool instrument = false;

    char c;
    char *optstring = "hg:o:n:s:I";
    while ((c = getopt(argc, argv, optstring)) != -1) {
        switch(c) {
        case 'h':
            usage(argv[0]);
            break;
        case 'g':
            gfile = fopen(optarg, "r");
            break;
        case 'o':
            ofile = fopen(optarg, "w");
            break;
        case 'n':
            count = atoi(optarg);
            break;
        case 's':
            seed = strtoul(optarg, NULL, 0);
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

    track_activity(instrument);
    START_ACTIVITY(ACTIVITY_STARTUP);
    if (gfile == NULL) {
        fprintf(stdout, "Couldn't open graph file\n");
        exit(1);
    }
    g = read_graph(gfile);
    if (g == NULL) {
        exit(1);
    }
    fclose(gfile);
    // init graph
    reset_charge(g);
    reset_boundary(g);
    srand(seed);

    fprintf(ofile, "%d %d %d\n", g->height, g->width, count);
    FINISH_ACTIVITY(ACTIVITY_STARTUP);

    simulate(g, count, ofile);

    SHOW_ACTIVITY(stderr, instrument);

    free_graph(g);
    fclose(ofile);

    return 0;
}

