#include "graphs.h"

Flow *edKarpSeq(Graph *g, int s, int t);

Flow *dinicSeq(Graph *g, int s, int t);

Flow *pushRelabelSeq(Graph *g, int s, int t);

int dinicSearch(Graph *g, int *flowMatrix, int *levels, int *curNeighbor, 
                int u, int t, int curFlow);