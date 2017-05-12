#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits>
#include <algorithm>
#include <omp.h>
#include "cpupar.h"
#include "sequential.h"

#define NOT_VISITED_MARKER -1

#define IDX(i, j, n) ((i) * (n) + (j))

struct vertexSet {
  int count;
  int maxVertices;
  int *vertices;
};

void vertexSetClear(vertexSet *list) {
    list->count = 0;
}

void vertexSetInit(vertexSet *list, int count) {
    list->maxVertices = count;
    list->vertices = (int *)malloc(sizeof(int) * list->maxVertices);
    vertexSetClear(list);
}

// Take one step of "top-down" BFS. Modified from assignment 3.
bool topDownStep(Graph *g, vertexSet *frontier, vertexSet *threadFrontiers,
                 int *flowMatrix, int *parents, int *pathCapacities, int t) {
    int maxThreads = omp_get_max_threads();
    bool foundT = false;

    #pragma omp parallel for schedule(static)
    for (int i=0; i < frontier->count; i++) {
        // hack to get early loop termination
        if (foundT) {
            i = frontier->count;
            continue;
        }
        int threadNum = omp_get_thread_num();

        int u = frontier->vertices[i];

        // attempt to add all neighbors to the thread's new frontier
        for (int v = 0; ((v < g->n) && (!foundT)); v++) {
            if (u == v) continue;
            int residual = g->capacities[IDX(u, v, g->n)] - flowMatrix[IDX(u, v, g->n)];
            if ((residual > 0) && (parents[v] == NOT_VISITED_MARKER) &&
                __sync_bool_compare_and_swap((&parents[v]), NOT_VISITED_MARKER, u)) {
                pathCapacities[v] = std::min(pathCapacities[u], residual);
                if (v != t) {
                    vertexSet *curFrontier = &threadFrontiers[threadNum];
                    int index = curFrontier->count++;
                    curFrontier->vertices[index] = v;
                } else {
                    foundT = true;
                }                
            }
        }
    }

    vertexSetClear(frontier);

    if (foundT) {
        return true;
    }

    // copy individual new frontiers back to main, single frontier
    #pragma omp parallel for schedule(static)
    for (int i=0; i < maxThreads; i++) {
        vertexSet *curFrontier = &threadFrontiers[i];
        int count = curFrontier->count;
        if (count > 0) {
            int index = __sync_fetch_and_add(&(frontier->count), count);
            for (int j = index; j < count + index; j++) {
                frontier->vertices[j] = curFrontier->vertices[j - index];
            }
        }
        vertexSetClear(curFrontier);
    }
    return false;
}

// Implements top-down BFS. Modified from assignment 3.
int bfsTopDown(Graph *g, vertexSet *frontier, vertexSet *threadFrontiers,
               int *flowMatrix, int *parents, int *pathCapacities, int s, int t) {
    int maxThreads = omp_get_max_threads();

    vertexSetClear(frontier);

    #pragma omp parallel for schedule(static)
    for (int i=0; i < maxThreads; i++) {
        vertexSetClear(&threadFrontiers[i]);
    }

    #pragma omp parallel for schedule(static)
    for (int i=0; i < g->n; i++) {
        parents[i] = NOT_VISITED_MARKER;
        pathCapacities[i] = 0;
    }

    frontier->vertices[frontier->count++] = s;
    parents[s] = s;
    pathCapacities[s] = std::numeric_limits<int>::max();

    while (frontier->count != 0) {
        bool foundT = topDownStep(g, frontier, threadFrontiers, flowMatrix, 
                             parents, pathCapacities, t);
        if (foundT) {
            return pathCapacities[t];
        }
    }
    return 0;
}

// Edmonds-Karp algorithm to find max s-t flow
Flow *edKarpCPU(Graph *g, int s, int t) {
    int flow = 0;
    int *flowMatrix = (int *)malloc((g->n * g->n) * sizeof(int));
    int *parents = (int *)malloc(g->n * sizeof(int));
    int *pathCapacities = (int *)malloc(g->n * sizeof(int));
    int maxThreads = omp_get_max_threads();

    vertexSet list1;
    vertexSetInit(&list1, g->n);

    vertexSet *frontier = &list1;

    vertexSet *threadFrontiers = (vertexSet *)malloc(maxThreads *
                                    sizeof(vertexSet));

    // initialize all thread frontiers
    #pragma omp parallel for schedule(static)
    for (int i=0; i < maxThreads; i++) {
        vertexSetInit(&threadFrontiers[i], g->n);
    }

    #pragma omp parallel for schedule(static)
    for (int i=0; i < (g->n * g->n); i++) {
        flowMatrix[i] = 0;
    }

    while (true) {
        int tempCapacity = bfsTopDown(g, frontier, threadFrontiers, flowMatrix, 
                                      parents, pathCapacities, s, t);
        if (tempCapacity == 0) {
            break;
        }
        flow += tempCapacity;
        int v = t;
        // backtrack
        while (v != s) {
            int u = parents[v];
            flowMatrix[IDX(u, v, g->n)] += tempCapacity;
            flowMatrix[IDX(v, u, g->n)] -= tempCapacity;
            v = u;
        }
    }
    Flow *result = (Flow *)malloc(sizeof(Flow));
    result->maxFlow = flow;
    result->finalEdgeFlows = flowMatrix;

    #pragma omp parallel for schedule(static)
    for (int i=0; i < maxThreads; i++) {
        free(threadFrontiers[i].vertices);
    }

    free(frontier->vertices);
    free(threadFrontiers);
    free(parents);
    free(pathCapacities);
    return result;
}

bool topDownStepDinic(Graph *g, vertexSet *frontier, vertexSet *threadFrontiers,
                      int *flowMatrix, int *levels, int t) {
    int maxThreads = omp_get_max_threads();
    bool foundT = false;

    #pragma omp parallel for schedule(static)
    for (int i=0; i < frontier->count; i++) {
        // hack to get early loop termination
        if (foundT) {
            i = frontier->count;
            continue;
        }
        int threadNum = omp_get_thread_num();

        int u = frontier->vertices[i];

        // attempt to add all neighbors to the thread's new frontier
        for (int v = 0; ((v < g->n) && (!foundT)); v++) {
            if (u == v) continue;
            int residual = g->capacities[IDX(u, v, g->n)] - flowMatrix[IDX(u, v, g->n)];
            if ((residual > 0) && (levels[v] == NOT_VISITED_MARKER) &&
                __sync_bool_compare_and_swap((&levels[v]), NOT_VISITED_MARKER, levels[u] + 1)) {
                if (v != t) {
                    vertexSet *curFrontier = &threadFrontiers[threadNum];
                    int index = curFrontier->count++;
                    curFrontier->vertices[index] = v;
                } else {
                    foundT = true;
                }                
            }
        }
    }

    vertexSetClear(frontier);

    if (foundT) {
        return true;
    }

    // copy individual new frontiers back to main, single frontier
    #pragma omp parallel for schedule(static)
    for (int i=0; i < maxThreads; i++) {
        vertexSet *curFrontier = &threadFrontiers[i];
        int count = curFrontier->count;
        if (count > 0) {
            int index = __sync_fetch_and_add(&(frontier->count), count);
            for (int j = index; j < count + index; j++) {
                frontier->vertices[j] = curFrontier->vertices[j - index];
            }
        }
        vertexSetClear(curFrontier);
    }
    return false;
}

// Subroutine for Dinic's algorithm, used to update the levels of nodes in the
// layered graph with BFS. Returns true if a path from s to t with positive
// capacity still exists.
bool bfsTopDownDinicLayered(Graph *g, vertexSet *frontier, vertexSet *threadFrontiers,
                            int *flowMatrix, int *levels, int s, int t) {
    int maxThreads = omp_get_max_threads();
    bool foundT = false;

    vertexSetClear(frontier);

    #pragma omp parallel for schedule(static)
    for (int i=0; i < maxThreads; i++) {
        vertexSetClear(&threadFrontiers[i]);
    }

    #pragma omp parallel for schedule(static)
    for (int i=0; i < g->n; i++) {
        levels[i] = NOT_VISITED_MARKER;
    }

    frontier->vertices[frontier->count++] = s;
    levels[s] = 0;

    while (frontier->count != 0) {
        foundT = topDownStepDinic(g, frontier, threadFrontiers, flowMatrix, 
                                  levels, t);
        if (foundT) {
            break;
        }
    }
    return foundT;
}

// Dinic's algorithm to find max s-t flow
Flow *dinicCPU(Graph *g, int s, int t) {
    int flow = 0;
    int maxInt = std::numeric_limits<int>::max();
    int *curNeighbor = (int *)malloc(g->n * sizeof(int));
    int *levels = (int *)malloc(g->n * sizeof(int));
    int *flowMatrix = (int *)calloc((g->n * g->n), sizeof(int));
    int maxThreads = omp_get_max_threads();

    vertexSet list1;
    vertexSetInit(&list1, g->n);

    vertexSet *frontier = &list1;

    vertexSet *threadFrontiers = (vertexSet *)malloc(maxThreads *
                                    sizeof(vertexSet));

    // initialize all thread frontiers
    #pragma omp parallel for schedule(static)
    for (int i=0; i < maxThreads; i++) {
        vertexSetInit(&threadFrontiers[i], g->n);
    }
    while (true) {
        #pragma omp parallel for schedule(static)
        for (int i=0; i < g->n; i++) {
            curNeighbor[i] = 0;
        }
        bool existsPath = bfsTopDownDinicLayered(g, frontier, threadFrontiers, 
                                                 flowMatrix, levels, s, t);
        if (!existsPath) {
            break;
        }
        while (true) {
            int tempFlow = dinicSearch(g, flowMatrix, levels, curNeighbor, s, t, maxInt);
            if (tempFlow == 0) {
                break;
            }
            flow += tempFlow;
        }
    }
    Flow *result = (Flow *)malloc(sizeof(Flow));
    result->maxFlow = flow;
    result->finalEdgeFlows = flowMatrix;

    #pragma omp parallel for schedule(static)
    for (int i=0; i < maxThreads; i++) {
        free(threadFrontiers[i].vertices);
    }

    free(frontier->vertices);
    free(threadFrontiers);
    free(curNeighbor);
    free(levels);
    return result;
}

// Push-relabel algorithm to find max s-t flow. Based on lock-free implementation
// specified by Bo Hong.
Flow *pushRelabelLockFree(Graph *g, int s, int t) {
    int *residualFlow = (int *)malloc((g->n * g->n) * sizeof(int));
    memcpy(residualFlow, g->capacities, (g->n * g->n) * sizeof(int));
    int *height = (int *)calloc(g->n, sizeof(int));
    int *excessFlow = (int *)calloc(g->n, sizeof(int));

    // new termination condition: net flow out of s equals net flow into t
    int netFlowOutS = 0;
    int netFlowInT = 0;

    // first, initialize preflow
    height[s] = g->n;
    #pragma omp parallel for reduction(+:netFlowOutS)
    for (int v = 0; v < g->n; v++) {
        int cap = residualFlow[IDX(s, v, g->n)];
        if (cap > 0 && (s != v)) {
            residualFlow[IDX(s, v, g->n)] = 0;
            residualFlow[IDX(v, s, g->n)] += cap;
            netFlowOutS += cap;
            excessFlow[v] = cap;
            if (v == t) {
                netFlowInT += cap;
            }
        }
    }

    // now break up vertices in chunks
    #pragma omp parallel for schedule(static)
    for (int i = 0; i < std::min(g->n, omp_get_num_threads()); i++) {
        int startVx, endVx;
        if (g->n < omp_get_num_threads()) {
            startVx = i;
            endVx = i+1;
        } else {
            int gap = g->n / omp_get_num_threads();
            startVx = i * gap;
            endVx = startVx + gap; // exclusive
            if (i == (omp_get_num_threads() - 1)) {
                endVx = g->n;
            }
        }
        
        while (netFlowOutS != netFlowInT) {
            for (int u = startVx; u < endVx; u++) {
                if ((u == s) || (u == t)) continue;
                while (excessFlow[u] > 0) {
                    int curExcess = excessFlow[u];
                    int curLowestNeighbor = -1;
                    int neighborMinHeight = std::numeric_limits<int>::max();
                    for (int v = 0; v < g->n; v++) {
                        if (u == v) continue;
                        if (residualFlow[IDX(u, v, g->n)] > 0) {
                            int tempHeight = height[v];
                            if (tempHeight < neighborMinHeight) {
                                curLowestNeighbor = v;
                                neighborMinHeight = tempHeight;
                            }
                        }
                    }
                    if (height[u] > neighborMinHeight) {
                        int delta = std::min(curExcess, residualFlow[IDX(u, curLowestNeighbor, g->n)]);
                        __sync_fetch_and_add(&residualFlow[IDX(u, curLowestNeighbor, g->n)], -delta);
                        __sync_fetch_and_add(&residualFlow[IDX(curLowestNeighbor, u, g->n)], delta);
                        __sync_fetch_and_add(&excessFlow[u], -delta);
                        __sync_fetch_and_add(&excessFlow[curLowestNeighbor], delta);
                        if (curLowestNeighbor == s) {
                            __sync_fetch_and_add(&netFlowOutS, -delta);
                        } else if (curLowestNeighbor == t) {
                            __sync_fetch_and_add(&netFlowInT, delta);
                        }
                    } else {
                        // doesn't need lock protection, only thread writing 
                        // to height[u]
                        height[u] = neighborMinHeight + 1;
                    }
                }
            }
        }
    }

    // now update residualFlow to represent actual flow
    #pragma omp parallel for schedule(static)
    for (int i = 0; i < (g->n * g-> n); i++) {
        residualFlow[i] = g->capacities[i] - residualFlow[i];
    }

    Flow *result = (Flow *)malloc(sizeof(Flow));
    result->maxFlow = netFlowInT;
    result->finalEdgeFlows = residualFlow;
    free(height);
    free(excessFlow);
    return result;
}