#include "graphs.h"

struct BfsResult;

int BFS(Graph *g, int *flowMatrix, int *parents, int s, int t);

bool bfsLayeredGraph(Graph *g, int *flowMatrix, int *levels, int s, int t);

int dinicSearch(Graph *g, int *flowMatrix, int *levels, int *curNeighbors, 
                int u, int t, int curFlow);

Flow *edKarpSeq(Graph *g, int s, int t);

Flow *dinicSeq(Graph *g, int s, int t);