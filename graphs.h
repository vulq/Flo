#ifndef GRAPH_H
#define GRAPH_H

struct Graph {
    int n;
    int *capacities; // n x n matrix, entry (i, j) is c(i, j)
};

struct Flow {
    int maxFlow;
    int *finalEdgeFlows; // n x n matrix, entry (i, j) is flow pushed on (i, j)
};

#endif 