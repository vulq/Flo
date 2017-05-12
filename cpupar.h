#include "graphs.h"

Flow *edKarpCPU(Graph *g, int s, int t);

Flow *dinicCPU(Graph *g, int s, int t);

Flow *pushRelabelLockFree(Graph *g, int s, int t);