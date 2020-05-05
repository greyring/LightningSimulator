#include <math.h>
#include <stdlib.h>
#include "graph.h"
#include <cuda.h>
#include <cuda_runtime.h>
#include <driver_functions.h>
#include "sim.h"
#include "instrument.h"




/* What is the crossover between binary and linear search */
#define BINARY_THRESHOLD 4
#define BLOCK_HEIGHT 32
#define BLOCK_WIDTH 32
#define BLOCK_SIZE (BLOCK_HEIGHT*BLOCK_WIDTH)
#define PLUS_BLOCK_SIZE ((BLOCK_HEIGHT+2)*(BLOCK_WIDTH+2))
#define GET_INDEX(row, col, width) ((row)*(width) + (col))

// This stores the global constants
struct GlobalConstants {

    int width;
    int height;
    int eta;

    double *charge;
    double *charge_buffer;

    double *boundary;
    int *bolt;

    double* choice_probs;
    int* choice_inv_map;
    int* choosed;

};

__constant__ GlobalConstants cuConstGraph;
GlobalConstants params;

int *choice_map;

/*
  Linear search
 */
static inline int locate_value_linear(double target, double *list, int len) {
    int i;
    for (i = 0; i < len; i++)
	    if (target < list[i])
	        return i;
    /* Shouldn't get here */
    return -1;
}

/*
  Binary search down to threshold, and then linear
 */
static inline int locate_value(double target, double *list, int len) {
    int left = 0;
    int right = len-1;
    while (left < right) {
	    if (right-left+1 < BINARY_THRESHOLD)
	        return left + locate_value_linear(target, list+left, right-left+1);
	    int mid = left + (right-left)/2;
	    if (target < list[mid])
	        right = mid;
	    else
	        left = mid+1;
    }
    return right;
}
static void reset_charge(graph_t *g) {
    int i;
    for (i = 0; i < g->height * g->width; i++) {
        g->charge[i] = g->charge_buffer[i] = 0;
    }
}

static void reset_boundary(graph_t *g) {
    int i;
    for (i = 0; i < g->height * g->width; i++) {
        g->boundary[i] = 0.0;
    }
}

static void reset_bolt(graph_t *g) {
    int i;
    for (i = 0; i < g->height * g->width; i++) {
        g->bolt[i] = g->reset_bolt[i];
    }
}

static void reset_path(graph_t *g) {
    int i;
    for (i = 0; i < g->height * g->width; i++) {
        g->path[i] = -1;
    }
}

static void choose_helper(graph_t *g, int bolt_idx, int i, int j) {
    int idx = i * g->width + j;
    if (i >= 0 && i < g->height && j >= 0 && j < g->width &&
        g->choosed[idx] == 0 && g->bolt[idx] <= 0) {
        g->choosed[idx] = 1;
        g->choice_idxs[g->num_choice] = idx;
        choice_map[idx] = g->num_choice;
        g->num_choice++;
        g->path[idx] = bolt_idx;
    }
}

static void find_choice(graph_t *g, int idx) {
    int i, j;
    if (g->bolt[idx] > 0) {
        i = idx / g->width;
        j = idx % g->width;
        choose_helper(g, idx, i - 1, j);
        choose_helper(g, idx, i, j - 1);
        choose_helper(g, idx, i, j + 1);
        choose_helper(g, idx, i + 1, j);
    }
}

static void reset_choice(graph_t *g) {
    int i;
    g->num_choice = 0;
    for (i = 0; i < g->height * g->width; i++) {
        g->choosed[i] = 0;
    }
    // get choices idxs
    for (i = 0; i < g->width * g->height; i++) {
        find_choice(g, i);
    }
}


__global__ void kernel_update_value(){
    int imageX = blockIdx.x * blockDim.x + threadIdx.x;
    int imageY = blockIdx.y * blockDim.y + threadIdx.y;
    int globalIdx = GET_INDEX(imageY, imageX, cuConstGraph.width);
    int linearThreadIndex = GET_INDEX(threadIdx.y+1, threadIdx.x+1, blockDim.x+2);
    __shared__ double old_charge[PLUS_BLOCK_SIZE];
    __shared__ double new_charge[BLOCK_SIZE];
    
    if(imageX < cuConstGraph.width && imageY < cuConstGraph.height){
        
        old_charge[linearThreadIndex] = cuConstGraph.charge[globalIdx];
        if(threadIdx.x == 0){
            old_charge[GET_INDEX(threadIdx.y+1, threadIdx.x, blockDim.x+2)] = imageX > 0 ? cuConstGraph.charge[GET_INDEX(imageY, imageX-1, cuConstGraph.width)] : 0;
        }
        if(threadIdx.x == blockDim.x-1){
            old_charge[GET_INDEX(threadIdx.y+1, threadIdx.x+2, blockDim.x+2)] = imageX < cuConstGraph.width-1 ? cuConstGraph.charge[GET_INDEX(imageY, imageX+1, cuConstGraph.width)] : 0;
        } else if (threadIdx.x == cuConstGraph.width - 1) {
            old_charge[GET_INDEX(threadIdx.y+1, threadIdx.x+2, blockDim.x+2)] = 0;
        }
        if(threadIdx.y == 0){
            old_charge[GET_INDEX(threadIdx.y, threadIdx.x+1, blockDim.x+2)] = imageY > 0 ? cuConstGraph.charge[GET_INDEX(imageY-1, imageX, cuConstGraph.width)] : 0;
        }
        if(threadIdx.y == blockDim.y-1){
            old_charge[GET_INDEX(threadIdx.y+2, threadIdx.x+1, blockDim.x+2)] = imageY < cuConstGraph.height-1 ? cuConstGraph.charge[GET_INDEX(imageY+1, imageX, cuConstGraph.width)] : 0;
        } else if (threadIdx.y == cuConstGraph.height - 1) {
            old_charge[GET_INDEX(threadIdx.y+2, threadIdx.x+1, blockDim.x+2)] = 0;
        }
    }
    __syncthreads();
    if(imageX < cuConstGraph.width && imageY < cuConstGraph.height){
        linearThreadIndex = GET_INDEX(threadIdx.y, threadIdx.x, blockDim.x);
        if(cuConstGraph.bolt[globalIdx] < 0){
            new_charge[linearThreadIndex] = 1.0;
        }else if(cuConstGraph.bolt[globalIdx] > 0){
            new_charge[linearThreadIndex] = 0.0;
        }else{
            new_charge[linearThreadIndex] = cuConstGraph.boundary[globalIdx]+ old_charge[GET_INDEX(threadIdx.y+1, threadIdx.x, blockDim.x+2)] + old_charge[GET_INDEX(threadIdx.y+1, threadIdx.x+2, blockDim.x+2)] + old_charge[GET_INDEX(threadIdx.y, threadIdx.x+1, blockDim.x+2)] + old_charge[GET_INDEX(threadIdx.y+2, threadIdx.x+1, blockDim.x+2)];
            new_charge[linearThreadIndex] /= 4;
        }
        cuConstGraph.charge_buffer[globalIdx] = new_charge[linearThreadIndex];
        if(cuConstGraph.choosed[globalIdx] == 1){
            cuConstGraph.choice_probs[cuConstGraph.choice_inv_map[globalIdx]] = cuConstGraph.bolt[globalIdx] > 0 ? 0 : pow(new_charge[linearThreadIndex], cuConstGraph.eta);
        }
    }
    __syncthreads();
}

__global__ void kernel_replace_charge(){
    int imageX = blockIdx.x * blockDim.x + threadIdx.x;
    int imageY = blockIdx.y * blockDim.y + threadIdx.y;
    int globalIdx = GET_INDEX(imageY, imageX, cuConstGraph.width);
    if(imageX < cuConstGraph.width && imageY < cuConstGraph.height){
        cuConstGraph.charge[globalIdx] = cuConstGraph.charge_buffer[globalIdx];
    }
}
__global__ void kernel_update_boundary(){
    int imageX = blockIdx.x * blockDim.x + threadIdx.x;
    int imageY = blockIdx.y * blockDim.y + threadIdx.y;
    int globalIdx = GET_INDEX(imageY, imageX, cuConstGraph.width);
    if(imageX < cuConstGraph.width && imageY < cuConstGraph.height){
        if (cuConstGraph.bolt[globalIdx] > 1) {
            cuConstGraph.boundary[globalIdx] = cuConstGraph.bolt[globalIdx] * 0.0001;
        } else {
            cuConstGraph.boundary[globalIdx] = 0;
        }
    }

}
// get bolt at x, y
// if bolt < 0.0, charge = 1.0 // boundary
// if bolt > 0.0, charge = 0.0 // boundary
// else charge = (boundary + neighbor's charge) / 4
static void update_charge(graph_t *g) {
    dim3 blockDim(BLOCK_WIDTH, BLOCK_HEIGHT); // 16*16 = 256
    dim3 gridDim((g->width+blockDim.x-1)/blockDim.x, (g->height+blockDim.y-1)/blockDim.y);
    START_ACTIVITY(ACTIVITY_UPDATE);
    kernel_update_value<<<gridDim, blockDim>>>();
    cudaDeviceSynchronize();

    // replace origin
    kernel_replace_charge<<<gridDim, blockDim>>>();

    cudaDeviceSynchronize();
    FINISH_ACTIVITY(ACTIVITY_UPDATE);
}
static void update_boundary(graph_t *g){
    START_ACTIVITY(ACTIVITY_RECOVER);
    dim3 blockDim(BLOCK_WIDTH, BLOCK_HEIGHT); // 16*16 = 256
    dim3 gridDim((g->width+blockDim.x-1)/blockDim.x, (g->height+blockDim.y-1)/blockDim.y);
    kernel_update_boundary<<<gridDim, blockDim>>>();
    cudaDeviceSynchronize();
    FINISH_ACTIVITY(ACTIVITY_RECOVER);

}
// add charge to bolt along the path
static void discharge(graph_t *g, int index, int charge) {
    int count = 500;
    while (index != -1 && count > 0) {
        count -= 1;
        g->bolt[index] += charge;
        cudaMemcpy(&(params.bolt[index]), &(g->bolt[index]), sizeof(int), cudaMemcpyHostToDevice);
        index = g->path[index];
    }
}
static __inline__ void update_kernel_choosed(graph_t *g, int i, int j){
    int idx = i * g->width + j;
    if(i >= 0 && i < g->height && j >= 0 && j < g->height && g->bolt[idx] <= 0){
        cudaMemcpy(&(params.choosed[idx]), &(g->choosed[idx]), sizeof(int), cudaMemcpyHostToDevice);
        cudaMemcpy(&(params.choice_inv_map[idx]), &(choice_map[idx]), sizeof(int), cudaMemcpyHostToDevice);
    }

}
static __inline__ void update_kernel_state(graph_t *g, int next_bolt){
    cudaMemcpy(&(params.bolt[next_bolt]), &(g->bolt[next_bolt]), sizeof(int), cudaMemcpyHostToDevice);
    int i = next_bolt / g->width;
    int j = next_bolt % g->width;
    update_kernel_choosed(g, i-1, j);
    update_kernel_choosed(g, i+1, j);
    update_kernel_choosed(g, i, j-1);
    update_kernel_choosed(g, i, j+1);
}
static void find_next(graph_t *g, int* power) {
    int idx, choice, next_bolt;
    double breach;

    START_ACTIVITY(ACTIVITY_NEXT);
    cudaMemcpy(g->choice_probs, params.choice_probs, sizeof(double)*g->num_choice, cudaMemcpyDeviceToHost);
    // calculate probability based on latest charge
    for(idx = 1; idx < g->num_choice; idx++){
        g->choice_probs[idx] += g->choice_probs[idx-1];
    }
    breach = (double)rand()/RAND_MAX * g->choice_probs[g->num_choice - 1];
    choice = locate_value(breach, g->choice_probs, g->num_choice);
    // choose one as bolt
    if (choice != -1){
        next_bolt = g->choice_idxs[choice];
        if (g->bolt[next_bolt] < 0) {
            *power += g->bolt[next_bolt];
            discharge(g, next_bolt, -g->bolt[next_bolt]);
        }
        g->bolt[next_bolt] = 1;
        find_choice(g, next_bolt);
        //copy new bolt and the choose;
        update_kernel_state(g, next_bolt);
    }
    FINISH_ACTIVITY(ACTIVITY_NEXT);
}

static void simulate_one(graph_t *g) {
    int power = g->power;
    int graphSize = g->width*g->height;


    START_ACTIVITY(ACTIVITY_RECOVER);
    reset_bolt(g);
    reset_path(g);
    reset_choice(g);
    cudaMemcpy(params.bolt, g->bolt, sizeof(int)*graphSize, cudaMemcpyHostToDevice);
    cudaMemcpy(params.choosed, g->choosed, sizeof(int)*graphSize, cudaMemcpyHostToDevice);
    cudaMemcpy(params.choice_inv_map, choice_map, sizeof(int)*graphSize, cudaMemcpyHostToDevice);
    FINISH_ACTIVITY(ACTIVITY_RECOVER);

    while (power > 0) {
        update_charge(g);
        find_next(g, &power);
    }
    // one lightning is generated
    update_boundary(g);
    
}

void simulate(graph_t *g, int count, FILE *ofile) {
    int i;

    int graphSize = g->width*g->height;

    double *cuda_charge_buffer;
    double *cuda_charge;
    double *cuda_boundary;
    int *cuda_bolt;
    double* cuda_choice_probs;
    int* cuda_choosed;
    int* cuda_choice_map;


    reset_bolt(g);
    reset_charge(g);
    reset_boundary(g);
    
    params.width = g->width;
    params.height = g->height;
    params.eta = g->eta;
    
    choice_map = (int*)malloc(sizeof(int)*graphSize);

    start_activity(ACTIVITY_STARTUP);
    cudaMalloc(&cuda_charge_buffer, sizeof(double)*graphSize);
    cudaMalloc(&cuda_charge, sizeof(double)*graphSize);
    cudaMalloc(&cuda_boundary, sizeof(double)*graphSize);
    cudaMalloc(&cuda_bolt, sizeof(int)*graphSize);
    cudaMalloc(&cuda_choice_probs, sizeof(double)*graphSize);
    cudaMalloc(&cuda_choosed, sizeof(int)*graphSize);
    cudaMalloc(&cuda_choice_map, sizeof(int)*graphSize);

    cudaMemcpy(cuda_charge_buffer, g->charge_buffer, sizeof(double)*graphSize, cudaMemcpyHostToDevice);
    cudaMemcpy(cuda_charge, g->charge, sizeof(double)*graphSize, cudaMemcpyHostToDevice);
    cudaMemcpy(cuda_boundary, g->boundary, sizeof(double)*graphSize, cudaMemcpyHostToDevice);
    cudaMemcpy(cuda_bolt, g->bolt, sizeof(int)*graphSize, cudaMemcpyHostToDevice);

    params.charge = cuda_charge;
    params.charge_buffer = cuda_charge_buffer;
    params.boundary = cuda_boundary;
    params.bolt = cuda_bolt;
    params.choice_probs = cuda_choice_probs;
    params.choosed = cuda_choosed;
    params.choice_inv_map = cuda_choice_map;

    cudaMemcpyToSymbol(cuConstGraph, &params, sizeof(GlobalConstants));
    finish_activity(ACTIVITY_STARTUP);
   
    for (i = 0; i < g->width + g->height; i++) {
        update_charge(g);
    }
   
    // generate lightnings
    for (i = 0; i < count; i++) {
        simulate_one(g);

        START_ACTIVITY(ACTIVITY_PRINT);
        // print bolt
        print_graph(g, ofile);
        fprintf(ofile, "\n");
        FINISH_ACTIVITY(ACTIVITY_PRINT);
    }
}
