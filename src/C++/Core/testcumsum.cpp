// !nvcc -D_FORCE_INLINES -O3 -c --compiler-options '-fPIC' -w gpuCumsum.cu -o gpuCumsum.o
//!ar -rvs gpuCumsum.a gpuCumsum.o
// g++ -std=c++11 testcumsum.cpp -o testcumsum gpuCumsum.a -O3 -lcudart -lcublas
        
#include <string>
#include <fstream>
#include <sstream>
#include <iostream>
#include <math.h>
#include <string.h>

/*
#include <stdlib.h>
#include <stdio.h>
#include <fstream>
#include <sstream>
#include <iostream>
#include <time.h>
*/
using namespace std;
#include <cuda_runtime.h>
#include <cuda.h>

// #include "cuda_runtime.h"
// #include "device_launch_parameters.h"
// #include "device_functions.h"

void gpuCumsum(int *output, int *input, int length, bool bcao);
void gpuCumsum(int *output, int *input, int *d_sums, int *d_incr, int length, bool optimized);

//void gpuArrayCopy(int *output, int *input, int length);
template <typename T> void gpuArrayCopy(T *y, T *x, int n);
void gpuPrintArray(int * a, int m);

#define CHECK(call)                                                            \
{                                                                              \
    const cudaError_t error = call;                                            \
    if (error != cudaSuccess)                                                  \
    {                                                                          \
        fprintf(stderr, "Error: %s:%d, ", __FILE__, __LINE__);                 \
        fprintf(stderr, "code: %d, reason: %s\n", error,                       \
                cudaGetErrorString(error));                                    \
        exit(1);                                                               \
    }                                                                          \
}

template <typename T> static void cudaTemplateMalloc(T **d_data, int n)
{
    // allocate the memory on the GPU            
    CHECK( cudaMalloc( (void**)d_data, n * sizeof(T) ) );
}

long get_nanos() {
	struct timespec ts;
	timespec_get(&ts, TIME_UTC);
	return (long)ts.tv_sec * 1000000000L + ts.tv_nsec;
}

void printResult(const char* prefix, int result, long nanoseconds) {
	printf("  ");
	printf(prefix);
	printf(" : %i in %ld ms \n", result, nanoseconds / 1000);
}

void printResult(const char* prefix, int result, long nanoseconds, long cpuseconds) {
	printf("  ");
	printf(prefix);
	printf(" : %i in %ld ms with speedup factor %ld \n", result, nanoseconds / 1000, cpuseconds/nanoseconds);
}

void printResult(const char* prefix, int result, float milliseconds) {
	printf("  ");
	printf(prefix);
	printf(" : %i in %f ms \n", result, milliseconds);
}

void printiarray(int* a, int m)
{    
    for (int i=0; i<m; i++)
        printf("%i  ", a[i]);
    printf("\n");
}

long cpu_cumsum(int* output, int* input, int length) {
	long start_time = get_nanos();

	output[0] = 0; // since this is a precumsum, not a cumsum
	for (int j = 1; j < length; ++j)
	{
		output[j] = input[j - 1] + output[j - 1];
	}

	long end_time = get_nanos();
	return end_time - start_time;
}

long gpu_cumsum(int* output, int* input, int length) {
	
	int *d_out, *d_in, *d_sums, *d_incr;
	const int arraySize = length * sizeof(int);

	cudaMalloc((void **)&d_out, arraySize);
	cudaMalloc((void **)&d_in, arraySize);
	cudaMemcpy(d_out, output, arraySize, cudaMemcpyHostToDevice);
	cudaMemcpy(d_in, input, arraySize, cudaMemcpyHostToDevice);
	    
    cudaMalloc((void **)&d_sums, length * sizeof(int)/1024);
	cudaMalloc((void **)&d_incr, length * sizeof(int)/1024);
    
    long start_time = get_nanos();    
    gpuCumsum(d_out, d_in, d_sums, d_incr, length, true);    
	long end_time = get_nanos();
    
    cudaMemcpy(input, d_in, arraySize, cudaMemcpyDeviceToHost);
    cudaMemcpy(output, d_out, arraySize, cudaMemcpyDeviceToHost);
    
    if (length<100) {
        printiarray(input, length);
        printiarray(output, length);
    }
    
	cudaFree(d_out);
	cudaFree(d_in);
    cudaFree(d_sums);
	cudaFree(d_incr);
    
	return end_time - start_time;
}

long gpu_cumsum_memory(int* output, int* input, int length) {
	
	int *d_out, *d_in, *d_sums, *d_incr;
	const int arraySize = length * sizeof(int);

	cudaMalloc((void **)&d_out, arraySize);
	cudaMalloc((void **)&d_in, arraySize);
	cudaMemcpy(d_out, output, arraySize, cudaMemcpyHostToDevice);
	cudaMemcpy(d_in, input, arraySize, cudaMemcpyHostToDevice);
	        
    long start_time = get_nanos();    
    gpuCumsum(d_out, d_in, length, true);    
	long end_time = get_nanos();
    
    cudaMemcpy(input, d_in, arraySize, cudaMemcpyDeviceToHost);
    cudaMemcpy(output, d_out, arraySize, cudaMemcpyDeviceToHost);
        
	cudaFree(d_out);
	cudaFree(d_in);
    
	return end_time - start_time;
}

/*///////////////////////////////////*/
/*            Main.cpp               */
/*///////////////////////////////////*/
void test(int N) {
	bool canBeBlockscanned = N <= 1024;

	time_t t;
	srand((unsigned)time(&t));
	int *in = new int[N+1];
    in[N] = 0;
	for (int i = 0; i < N; i++) {
		in[i] = rand() % 10;
	}

	printf("%i Elements \n", N);

	// sequential cumsum on CPU
	int *outHost = new int[N+1]();
	long time_host = cpu_cumsum(outHost, in, N+1);
	printResult("CPU Performance", outHost[N], time_host);

	// gpu cumsum without allocating memory during the iteration
	int *outGPU = new int[N+1]();
	long time_gpu = gpu_cumsum(outGPU, in, N+1);
	printResult("Pre-allocated Memory GPU Performance", outGPU[N], time_gpu, time_host);
    
    // gpu cumsum with memory allocated during the iteration.
	long time_gpu_memory = gpu_cumsum_memory(outGPU, in, N+1);
	printResult("Allocated Memory GPU Performance", outGPU[N], time_gpu_memory, time_host);
    
	printf("\n");

	delete[] in;
	delete[] outHost;
	delete[] outGPU;
}

int main()
{
	int TEN_MILLION = 10000000;
	int ONE_MILLION = 1000000;
	int TEN_THOUSAND = 10000;

	int elements[] = {
		TEN_MILLION * 2,
		TEN_MILLION,
		ONE_MILLION,
		TEN_THOUSAND,
		5000,
		4096,
		2048,
		2000,
		1000,
		500,
		100,
		64,
		8,
		5
	};

	int numElements = sizeof(elements) / sizeof(elements[0]);

    int shmrank = 0;
    int device;    
    cudaSetDevice(shmrank); 
    cudaGetDevice( &device );
    size_t available, total;
    cudaMemGetInfo(&available, &total);
    cout<<"Available GPU Memory: "<<available<<" Total GPU Memory: "<<total<<endl;
    
	for (int i = 0; i < numElements; i++) {
		test(elements[i]);
	}

	return 0;
}
