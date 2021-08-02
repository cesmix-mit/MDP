/***************************************************************************                               
                    Molecular Dynamics Potentials (MDP)
                           CESMIX-MIT Project  
 
 Contributing authors: Ngoc-Cuong Nguyen (cuongng@mit.edu, exapde@gmail.com)
 ***************************************************************************/

#ifndef MDP_CUDAMEM_H
#define MDP_CUDAMEM_H

#define CUDA_CHECK(call)                                                            \
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

#define CHECK_CUBLAS(call)                                                     \
{                                                                              \
    cublasStatus_t err;                                                        \
    if ((err = (call)) != CUBLAS_STATUS_SUCCESS)                               \
    {                                                                          \
        fprintf(stderr, "Got CUBLAS error %d at %s:%d\n", err, __FILE__,       \
                __LINE__);                                                     \
        exit(1);                                                               \
    }                                                                          \
}

#define GPUFREE(x)                                                       \
{                                                                         \
    if (x != NULL) {                                                      \
        cudaTemplateFree(x);                                              \
        x = NULL;                                                         \
    }                                                                     \
}

#define CUDA_SYNC cudaDeviceSynchronize();  

template <typename T> static void cudaTemplateMalloc(T **d_data, int n)
{
    // allocate the memory on the GPU            
    CUDA_CHECK( cudaMalloc( (void**)d_data, n * sizeof(T) ) );
}

template <typename T> static void cudaTemplateMallocManaged(T **d_data, int n)
{
    // allocate unified memory 
    CUDA_CHECK( cudaMallocManaged( (void**)d_data, n * sizeof(T) ) );        
}

template <typename T> static void cudaTemplateHostAlloc(T **h_data, int n, unsigned int flags)
{
    // allocate zero-copy memory on host    
    CUDA_CHECK(cudaHostAlloc((void **)h_data, n * sizeof(T), flags));                
}

template <typename T> static void cudaTemplateHostAllocMappedMemory(T **h_data, int n)
{
    // allocate zero-copy memory on host    
    CUDA_CHECK(cudaHostAlloc((void **)h_data, n * sizeof(T), cudaHostAllocMapped));                
}

template <typename T> static void cudaTemplateHostAllocPinnedMemory(T **h_data, int n)
{
    // allocate pinned memory on host        
    CUDA_CHECK(cudaHostAlloc((void **)h_data, n * sizeof(T), cudaHostAllocDefault));                
}

template <typename T> static void cudaTemplateFree(T *d_data)
{
    // free the memory on the GPU            
    CUDA_CHECK( cudaFree( d_data ) );    
}

template <typename T> static void cudaCopytoDevice(T *d_data, T *h_data, int n)
{
    // copy data from CPU to GPU
    CUDA_CHECK( cudaMemcpy( d_data, h_data, n * sizeof(T), cudaMemcpyHostToDevice ) );    
}

template <typename T> static void cudaCopytoHost(T *h_data, T *d_data, int n)
{
    // copy data from GPU to CPU
    CUDA_CHECK( cudaMemcpy( h_data, d_data, n * sizeof(T), cudaMemcpyDeviceToHost ) );    
}

#endif

