#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <cuda.h>
#include <omp.h>
#include <cuda_runtime.h>
#include <driver_functions.h>
#include <math_constants.h>

#include "pushRelabelGPU.h"

#ifdef DEBUG
#define cudaCheckError(ans) { cudaAssert((ans), __FILE__, __LINE__); }
inline void cudaAssert(cudaError_t code, const char *file, int line,
        bool abort=true)
{
   if (code != cudaSuccess)
   {
      fprintf(stderr, "CUDA Error: %s at %s:%d\n",
        cudaGetErrorString(code), file, line);
      if (abort) exit(code);
   }
}
#else
#define cudaCheckError(ans) ans
#endif

#define IDX(i, j, n) ((i) * (n) + (j))

#define UPDIV(n, d)   (((n) + (d) - 1) / (d))

static dim3 threadsPerBlock(1024, 1, 1);

__global__ void pushRelabelLockFreeKernel(int *residualFlow,
        int *height, int *excessFlow, int *netFlowOutS, int *netFlowInT, 
        int s, int t, int n) {
    int index = blockDim.x * blockIdx.x + threadIdx.x;

    int u = index;
    if (u >= s) {
        u++;
    }
    if (u >= t) {
        u++;
    }

    // one thread here for all vertices not s or t
    while (*netFlowOutS != *netFlowInT) {
        if (u < n && excessFlow[u] > 0) {
            int curExcess = excessFlow[u];
            int curLowestNeighbor = -1;
            int neighborMinHeight = (int) CUDART_INF;
            for (int v = 0; v < n; v++) {
                if (u == v) continue;
                if (residualFlow[IDX(u, v, n)] > 0) {
                    int tempHeight = height[v];
                    if (tempHeight < neighborMinHeight) {
                        curLowestNeighbor = v;
                        neighborMinHeight = tempHeight;
                    }
                }
            }
            if (height[u] > neighborMinHeight) {
                int delta = min(curExcess, residualFlow[IDX(u, curLowestNeighbor, n)]);
                atomicSub(&residualFlow[IDX(u, curLowestNeighbor, n)], delta);
                atomicAdd(&residualFlow[IDX(curLowestNeighbor, u, n)], delta);
                atomicSub(&excessFlow[u], delta);
                atomicAdd(&excessFlow[curLowestNeighbor], delta);
                if (curLowestNeighbor == s) {
                    atomicSub(netFlowOutS, delta);
                } else if (curLowestNeighbor == t) {
                    atomicAdd(netFlowInT, delta);
                }
            } else {
                height[u] = neighborMinHeight + 1;
            }
        }
    }
}

// Push-relabel algorithm to find max s-t flow. Based on lock-free implementation
// specified by Bo Hong. Uses one CUDA thread per vertex.
Flow *pushRelabelLockFreeGPU(Graph *g, int s, int t) {
    int *residualFlow;
    int *height;
    int *excessFlow;
    int *netFlowOutS;
    int *netFlowInT;
    int *tempHeights = (int *)calloc(g->n,  sizeof(int));
    int *tempExcessFlows = (int *)calloc(g->n,  sizeof(int));
    int *finalFlow = (int *)malloc((g->n * g->n) * sizeof(int));
    memcpy(finalFlow, g->capacities, (g->n * g->n) * sizeof(int));

    cudaCheckError(cudaMalloc((void **)&residualFlow, sizeof(int) * (g->n * g->n)));
    cudaCheckError(cudaMalloc((void **)&height, sizeof(int) * g->n));
    cudaCheckError(cudaMalloc((void **)&excessFlow, sizeof(int) *  g->n));
    cudaCheckError(cudaMalloc((void **)&netFlowOutS, sizeof(int)));
    cudaCheckError(cudaMalloc((void **)&netFlowInT, sizeof(int)));

    // initialize preflow
    int flowOutS = 0;
    int flowInT = 0;
    tempHeights[s] = g->n;
    #pragma omp parallel for reduction(+:flowOutS)
    for (int v = 0; v < g->n; v++) {
        int cap = g->capacities[IDX(s, v, g->n)];
        if (cap > 0 && (s != v)) {
            finalFlow[IDX(s, v, g->n)] = 0;
            finalFlow[IDX(v, s, g->n)] += cap;
            flowOutS += cap;
            tempExcessFlows[v] = cap;
            if (v == t) {
                flowInT += cap;
            }
        }
    }

    cudaCheckError(cudaMemcpy(residualFlow, finalFlow, sizeof(int) * (g->n * g->n),
                   cudaMemcpyHostToDevice));
    cudaCheckError(cudaMemcpy(height, tempHeights, sizeof(int) * g->n,
                   cudaMemcpyHostToDevice));
    cudaCheckError(cudaMemcpy(excessFlow, tempExcessFlows, sizeof(int) * g->n,
                   cudaMemcpyHostToDevice));
    cudaCheckError(cudaMemcpy(netFlowInT, &flowInT, sizeof(int),
                   cudaMemcpyHostToDevice));
    cudaCheckError(cudaMemcpy(netFlowOutS, &flowOutS, sizeof(int),
                   cudaMemcpyHostToDevice));

    int numBlocks = UPDIV((g->n - 2), threadsPerBlock.x);
    pushRelabelLockFreeKernel<<<numBlocks, threadsPerBlock>>>(residualFlow,
        height, excessFlow, netFlowOutS, netFlowInT, s, t, g->n);

    free(tempHeights);
    free(tempExcessFlows);

    cudaCheckError(cudaThreadSynchronize());

    cudaCheckError(cudaMemcpy(finalFlow, residualFlow, sizeof(int) * (g->n * g->n),
               cudaMemcpyDeviceToHost));

    cudaCheckError(cudaMemcpy(&flowInT, netFlowInT, sizeof(int),
               cudaMemcpyDeviceToHost));

    // now update flow to represent actual flow
    #pragma omp parallel for schedule(static)
    for (int i = 0; i < (g->n * g-> n); i++) {
        finalFlow[i] = g->capacities[i] - finalFlow[i];
    }

    Flow *result = (Flow *)malloc(sizeof(Flow));
    result->maxFlow = flowInT;
    result->finalEdgeFlows = finalFlow;
    cudaCheckError(cudaFree(residualFlow));
    cudaCheckError(cudaFree(height));
    cudaCheckError(cudaFree(excessFlow));
    cudaCheckError(cudaFree(netFlowOutS));
    cudaCheckError(cudaFree(netFlowInT));
    return result;
}