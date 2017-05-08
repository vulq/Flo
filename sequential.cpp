#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <array>
#include <queue>
#include <limits>
#include "sequential.h"

#define IDX(i, j, n) ((i) * (n) + (j))

struct BfsResult {
    int capacity;
    int *parents;
};

BfsResult *BFS(Graph *g, int *flowMatrix, int s, int t) {
    int *parents = (int *)malloc(g->n * sizeof(int));
    memset(parents, -1, (g->n * sizeof(int)));
    parents[s] = s;
    int *pathCapacities = (int *)calloc(g->n, sizeof(int));
    pathCapacities[s] = std::numeric_limits<int>::max();
    std::queue<int> bfsQueue;
    bfsQueue.push(s);
    while (!bfsQueue.empty()) {
        int u = bfsQueue.front();
        bfsQueue.pop();
        for (int v = 0; v < g->n; v++) {
            if (u == v) continue;
            int residual = g->capacities[IDX(u, v, g->n)] - flowMatrix[IDX(u, v, g->n)];
            if ((residual > 0) && (parents[v] == -1)) {
                parents[v] = u;
                pathCapacities[v] = std::min(pathCapacities[u], residual);
                if (v != t) {
                    bfsQueue.push(v);
                } else {
                    BfsResult *result = (BfsResult *)malloc(sizeof(BfsResult));
                    result->capacity = pathCapacities[t];
                    result->parents = parents;
                    free(pathCapacities);
                    return result;
                }
            }
        }
    }
    BfsResult *result = (BfsResult *)malloc(sizeof(BfsResult));
    result->capacity = 0;
    result->parents = parents;
    return result;
}

// Edmonds-Karp algorithm to find max s-t flow.
Flow *edKarpSeq(Graph *g, int s, int t) {
    int flow = 0;
    int *flowMatrix = (int *)calloc((g->n * g->n), sizeof(int));
    while (true) {
        BfsResult *tempResult = BFS(g, flowMatrix, s, t);
        if ((tempResult->capacity) == 0) {
            break;
        }
        flow += tempResult->capacity;
        int v = t;
        // backtrack
        while (v != s) {
            int u = tempResult->parents[v];
            flowMatrix[IDX(u, v, g->n)] += tempResult->capacity;
            flowMatrix[IDX(v, u, g->n)] -= tempResult->capacity;
            v = u;
        }
        free(tempResult->parents);
        free(tempResult);
    }
    Flow *result = (Flow *)malloc(sizeof(Flow));
    result->maxFlow = flow;
    result->finalEdgeFlows = flowMatrix;
    return result;
}

// Dinic's algorithm to find max s-t flow.
Flow *dinicSeq(Graph *g, int s, int t) {
    return;
}