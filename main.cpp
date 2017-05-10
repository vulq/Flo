#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include <time.h>
#include <random>
#include <vector>
#include <iostream>
#include <algorithm>
#include "CycleTimer.h"
#include "sequential.h"

#define IDX(i, j, n) ((i) * (n) + (j))

// TODO (maybe): test with one disconnect graph to ensure flow is 0?
// Generates a random directed graph that ensures there exists some s-t path, where
// s is node 0 and t is node numVxs - 1. 
Graph *generateGraph(int numVxs, int numEdges, int maxCap) {
    Graph *g = (Graph *)malloc(sizeof(Graph));
    g->n = numVxs;
    g->capacities = (int *)calloc((numVxs * numVxs), sizeof(int));
    std::vector<int> path;
    for (int i = 1; i < numVxs; i++) {
        path.push_back(i);
    }
    std::random_shuffle(path.begin(), path.end());
    int t = numVxs - 1;
    int first = path.at(0);
    int capacity = (rand() % (maxCap - 1)) + 1;
    g->capacities[IDX(0, first, g->n)] = capacity;
    int remaining = numEdges - 1;
    int last = first;
    for (auto it = (path.begin()++); ((it != path.end()) && (*it != t)); it++) {
        capacity = (rand() % (maxCap - 1)) + 1;
        g->capacities[IDX(last, (*it), g->n)] = capacity;
        last = *it;
        remaining--;
    }
    capacity = (rand() % (maxCap - 1)) + 1;
    g->capacities[IDX(last, t, g->n)] = capacity;
    remaining--;
    for (int i = 0; i < remaining; i++) {
        capacity = (rand() % (maxCap - 1)) + 1;
        int j = rand() % g->n;
        int k = rand() % g->n;
        g->capacities[IDX(j, k, g->n)] = capacity;
    }
    return g;
}

// for each non s or t node, check flow in == flow out, and check flow
// out of s equals flow into t equals total flow. Also check that flows
// are symmetric (i.e. F[u][v] == -F[v][u])
bool checkFlow(int totalFlow, int *flows, int n) {
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            if (flows[IDX(i, j, n)] != -flows[IDX(j, i, n)]) {
                return false;
            }
        }
    }
    for (int i = 1; i < (n-1); i++) {
        int flowIn = 0;
        int flowOut = 0;
        for (int j = 0; j < n; j++) {
            int edgeFlow = flows[IDX(i, j, n)];
            if (edgeFlow > 0) {
                flowIn += edgeFlow;
            } else {
                flowOut -= edgeFlow;
            }
        }
        if (flowIn != flowOut) {
            return false;
        }
    }
    int sFlow = 0;
    int tFlow = 0;
    for (int i = 0; i < n; i++) {
        sFlow += flows[IDX(0, i, n)];
        tFlow += flows[IDX(i, (n-1), n)];
    }
    return (sFlow == tFlow) && (sFlow == totalFlow);
}

int main() {
    int refFlow;
    Flow *result;
    bool check;
    srand(0); // TODO: change this to time(NULL) for randomness btwn trials
    int numGraphs = 5;
    // TODO: also maybe use different large graphs for benchmarking different algs? i.e. use
    // the 20k ones and maybe even more, larger (higher V?) ones for Dinic's and Push-relabel 
    // but not EdKarp
    int smallGraphNum = 1000; // TODO: shrink/remove this when have correctness, don't run
                              // all these tests when just trying to benchmark
    int totalGraphs = numGraphs + smallGraphNum;
    double start, finalTime;
    int numVxs[] = {500, 10000, 10000, 20000, 20000};
    int numEdges[] = {10000, 50000, 1000000, 12000000, 100000000};
    int maxCap[] = {500, 100, 50, 30, 20};
    Graph *graphs[totalGraphs];
    double edKarpSeqTimes[numGraphs];
    double dinicSeqTimes[numGraphs];
    double pushRelabelSeqTimes[numGraphs];
    for (int i = 0; i < numGraphs; i++) {
        graphs[i] = generateGraph(numVxs[i], numEdges[i], maxCap[i]);
    }
    // generate small graphs too
    for (int i = numGraphs; i < totalGraphs; i++) {
        int vxs = (rand() % 1000) + 100;
        int maxEdges = (vxs * (vxs - 1)) / 2;
        int edges = (rand() % 20000) + vxs;
        edges = std::min(edges, maxEdges);
        int cap = (rand() % 1000) + 20;
        graphs[i] = generateGraph(vxs, edges, cap);
    }
    
    for (int i = 0; i < totalGraphs; i++) {
        refFlow = -1;

        // first, test Edmonds-Karp
        if (i < 3 || i > 4) {
            start = CycleTimer::currentSeconds();
            result = edKarpSeq(graphs[i], 0, (graphs[i]->n)-1);
            finalTime = CycleTimer::currentSeconds() - start;
            if (i < numGraphs) {
                edKarpSeqTimes[i] = finalTime;
                std::cout << "edKarpSeq time: " << edKarpSeqTimes[i] << std::endl;
                std::cout << "edKarpSeq flow: " << result->maxFlow << std::endl;
            }
            check = checkFlow(result->maxFlow, result->finalEdgeFlows, graphs[i]->n);
            if (!check) {
                std::cout << "EdKarp flows don't agree with max flow on graph " << i << std::endl;
            }
            refFlow = result->maxFlow;

            free(result->finalEdgeFlows);
            free(result);
        }
        
        // now Dinics
        start = CycleTimer::currentSeconds();
        result = dinicSeq(graphs[i], 0, (graphs[i]->n)-1);
        finalTime = CycleTimer::currentSeconds() - start;
        if (i < numGraphs) {
            dinicSeqTimes[i] = finalTime;
            std::cout << "dinicSeq time: " << dinicSeqTimes[i] << std::endl;
            std::cout << "dinicSeq flow: " << result->maxFlow << std::endl;
        }
        check = checkFlow(result->maxFlow, result->finalEdgeFlows, graphs[i]->n);
        if (!check) {
            std::cout << "Dinic flows don't agree with max flow on graph " << i << std::endl;
        }
        if ((refFlow != -1) && (result->maxFlow != refFlow)) {
            std::cout << "Dinic flow doesn't agree with refFlow on graph " << i << std::endl;
        }

        free(result->finalEdgeFlows);
        free(result);

        // now Push-relabel
        if (i < 3 || i > 4) {
            start = CycleTimer::currentSeconds();
            result = pushRelabelSeq(graphs[i], 0, (graphs[i]->n)-1);
            finalTime = CycleTimer::currentSeconds() - start;
            if (i < numGraphs) {
                pushRelabelSeqTimes[i] = finalTime;
                std::cout << "pushRelabelSeq time: " << pushRelabelSeqTimes[i] << std::endl;
                std::cout << "pushRelabelSeq flow: " << result->maxFlow << std::endl;
            }
            check = checkFlow(result->maxFlow, result->finalEdgeFlows, graphs[i]->n);
            if (!check) {
                std::cout << "Push-relabel flows don't agree with max flow on graph " << i << std::endl;
            }
            if ((refFlow != -1) && (result->maxFlow != refFlow)) {
                std::cout << "Push-relabel flow doesn't agree with refFlow on graph " << i << std::endl;
            }

            free(result->finalEdgeFlows);
            free(result);
        }

        // TODO: repeat this for CPU and GPU, and then compare times, output speedup, etc.

        // TODO: when want to test speedup with different number of CPU/GPU threads, can do 
        // dynamically here, or need to recompile with different numthreads constants?
    }

    for (int i = 0; i < totalGraphs; i++) {
        free(graphs[i]->capacities);
        free(graphs[i]);
    }
    return 0;
}
