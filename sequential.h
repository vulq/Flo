#include "graphs.h"

struct BfsResult;

BfsResult *BFS(Graph *g, int *flowMatrix, int s, int t);

Flow *edKarpSeq(Graph *g, int s, int t);

Flow *dinicSeq(Graph *g, int s, int t);