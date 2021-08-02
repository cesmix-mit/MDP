/***************************************************************************                               
                    Molecular Dynamics Potentials (MDP)
                           CESMIX-MIT Project  
 
 Contributing authors: Ngoc-Cuong Nguyen (cuongng@mit.edu, exapde@gmail.com)
 ***************************************************************************/

#include "hip/hip_runtime.h"

#ifndef __HIPARRAYPERMUTE
#define __HIPARRAYPERMUTE

template <typename T> __global__ void hipTemplateKron(T *C, T *A, T *B, 
        int M1, int M2, int M)
{                    
    int idx = threadIdx.x + blockIdx.x * blockDim.x; // global thread index  
    while (idx<M) {    
        int ib = idx%M2;
        int ia = (idx-ib)/M2;        
        C[idx] += A[ia]*B[ib];      
        idx += blockDim.x * gridDim.x;
    }
}
template <typename T> void hipKron(T *C, T *A, T *B, int M1, int M2)
{                
    int M = M1*M2;    
    int BLOCKDIM = 256;
    int gridDim = (M + BLOCKDIM - 1) / BLOCKDIM;
    gridDim = (gridDim>1024)? 1024 : gridDim;
    hipLaunchKernelGGL(hipTemplateKron, gridDim, BLOCKDIM, 0, 0, C, A, B, M1, M2, M);
}

template <typename T> __global__ void hipTemplateKron(T *C, T *A, T *B, 
        int M1, int N1, int M2, int N2, int M, int N, int P)
{            
    int idx = threadIdx.x + blockIdx.x * blockDim.x; // global thread index  
    while (idx<P) {
        int i = idx%M;
        int j = (idx-i)/M;
        int ib = i%M2;
        int ia = (i-ib)/M2;
        int jb = j%N2;
        int ja = (j-jb)/N2;
        C[idx] = A[ia+M1*ja]*B[ib+M2*jb];
        idx += blockDim.x * gridDim.x;
    }
}
template <typename T> void hipKron(T *C, T *A, T *B, int M1, int N1, int M2, int N2)
{                
    int M = M1*M2;
    int N = N1*N2;
    int P = M*N;    
    int BLOCKDIM = 256;
    int gridDim = (P + BLOCKDIM - 1) / BLOCKDIM;
    gridDim = (gridDim>1024)? 1024 : gridDim;
    hipLaunchKernelGGL(hipTemplateKron, gridDim, BLOCKDIM, 0, 0, C, A, B, M1, N1, M2, N2, M, N, P);
}

__global__ void hipKernelIndexPermute12(int *index, int I1, int I2, int I3)
{                 
    int idx = threadIdx.x + blockIdx.x * blockDim.x; // global thread index  
    int M = I1*I2;
    int N = M*I3;
    while (idx<N) {
        int l = idx%M;
        int i = l%I1;
        int j = (l-i)/I1;
        int k = (idx-l)/M;        
        index[idx] = j+I2*i+M*k;
        idx += blockDim.x * gridDim.x;
    }
}

void hipIndexPermute12(int *index, int I1, int I2, int I3)
{                
    int BLOCKDIM = 256;
    int gridDim = (I1*I2*I3 + BLOCKDIM - 1) / BLOCKDIM;
    gridDim = (gridDim>1024)? 1024 : gridDim;
    hipLaunchKernelGGL(hipKernelIndexPermute12, gridDim, BLOCKDIM, 0, 0, index, I1, I2, I3);
}

__global__ void hipKernelIndexPermute13(int *index, int I1, int I2, int I3, int I4)
{            
    int idx = threadIdx.x + blockIdx.x * blockDim.x; // global thread index  
    int M = I1*I2;
    int N = M*I3;
    int P = N*I4;
    while (idx<P) {
        int n = idx%N;
        int i4 = (idx-n)/N;        
        int m = n%M;
        int i3 = (n-m)/M;      
        int i1 = m%I1;
        int i2 = (m-i1)/I1;      
        index[idx] = i3+I3*i2+I3*I2*i1+N*i4;
        idx += blockDim.x * gridDim.x;
    }
}

void hipIndexPermute13(int *index, int I1, int I2, int I3, int I4)
{                
    int BLOCKDIM = 256;
    int gridDim = (I1*I2*I3*I4 + BLOCKDIM - 1) / BLOCKDIM;
    gridDim = (gridDim>1024)? 1024 : gridDim;
    hipLaunchKernelGGL(hipKernelIndexPermute13, gridDim, BLOCKDIM, 0, 0, index, I1, I2, I3, I4);
}

__global__ void hipKernelIndexPermute23(int *index, int I1, int I2, int I3, int I4)
{            
    int idx = threadIdx.x + blockIdx.x * blockDim.x; // global thread index  
    int M = I1*I2;
    int N = M*I3;
    int P = N*I4;
    while (idx<P) {
        int n = idx%N;
        int i4 = (idx-n)/N;        
        int m = n%M;
        int i3 = (n-m)/M;      
        int i1 = m%I1;
        int i2 = (m-i1)/I1;      
        index[idx] = i1+I1*i3+I1*I3*i2+N*i4;
        idx += blockDim.x * gridDim.x;
    }
}

void hipIndexPermute23(int *index, int I1, int I2, int I3, int I4)
{                
    int BLOCKDIM = 256;
    int gridDim = (I1*I2*I3*I4 + BLOCKDIM - 1) / BLOCKDIM;
    gridDim = (gridDim>1024)? 1024 : gridDim;
    hipLaunchKernelGGL(hipKernelIndexPermute23, gridDim, BLOCKDIM, 0, 0, index, I1, I2, I3, I4);
}

template <typename T> __global__ void hipTemplatePermute(T *B, T *A, int *index, int N)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x; // global thread index      
    while (idx<N) {
        B[index[idx]] = A[idx];
        idx += blockDim.x * gridDim.x;
    }
}

template <typename T> void hipPermute(T *B, T *A, int *index, int N)
{                
    int BLOCKDIM = 256;
    int gridDim = (N + BLOCKDIM - 1) / BLOCKDIM;
    gridDim = (gridDim>1024)? 1024 : gridDim;
    hipLaunchKernelGGL(hipTemplatePermute, gridDim, BLOCKDIM, 0, 0, B, A, index, N);
}

template <typename T> __global__ void hipTemplatePermuteSharedMem(T *B, T *A, int *index, int N)
{
    int tx = threadIdx.x;
    int nx = blockDim.x * gridDim.x;    
    int idx = tx + blockIdx.x * blockDim.x; // global thread index      
    __shared__ T Ashared[1025];     
    while (idx<N) {
        Ashared[index[tx]] = A[idx]; 
        __syncthreads();
        B[idx] = Ashared[tx];
        idx += nx;
    }
}

template <typename T> void hipPermuteSharedMem(T *B, T *A, int *index, int N, int BLOCKDIM)
{                    
    int gridDim = (N + BLOCKDIM - 1) / BLOCKDIM;
    //gridDim = (gridDim>1024)? 1024 : gridDim;
    hipLaunchKernelGGL(hipTemplatePermuteSharedMem, gridDim, BLOCKDIM, 0, 0, B, A, index, N);
}

template <typename T> __global__ void hipTemplatePermute12(T *B, T *A, int I1, int I2, int I3)
{            
    //int tx = threadIdx.x;                   // block thread index    
    int idx = threadIdx.x + blockIdx.x * blockDim.x; // global thread index  
    int M = I1*I2;
    int N = M*I3;
    
    while (idx<N) {
        int l = idx%M;
        int i = l%I1;
        int j = (l-i)/I1;
        int k = (idx-l)/M;        
        B[j+I2*i+M*k] = A[idx];
    }
}

template <typename T> void hipPermute12(T *B, T *A, int I1, int I2, int I3)
{                
    int BLOCKDIM = 256;
    int gridDim = (I1*I2*I3 + BLOCKDIM - 1) / BLOCKDIM;
    gridDim = (gridDim>1024)? 1024 : gridDim;
    hipLaunchKernelGGL(hipTemplatePermute12, gridDim, BLOCKDIM, 0, 0, B, A, I1, I2, I3);
}


template <typename T> __global__ void hipTemplatePermute13(T *B, T *A, int I1, int I2, int I3, int I4)
{            
    int idx = threadIdx.x + blockIdx.x * blockDim.x; // global thread index  
    int M = I1*I2;
    int N = M*I3;
    int P = N*I4;
    while (idx<P) {
        int n = idx%N;
        int i4 = (idx-n)/N;        
        int m = n%M;
        int i3 = (n-m)/M;      
        int i1 = m%I1;
        int i2 = (m-i1)/I1;      
        B[i3+I3*i2+I3*I2*i1+N*i4] = A[idx];
        idx += blockDim.x * gridDim.x;
    }
}

template <typename T> void hipPermute13(T *B, T *A, int I1, int I2, int I3, int I4)
{                
    int BLOCKDIM = 256;
    int gridDim = (I1*I2*I3 + BLOCKDIM - 1) / BLOCKDIM;
    gridDim = (gridDim>1024)? 1024 : gridDim;
    hipLaunchKernelGGL(hipTemplatePermute13, gridDim, BLOCKDIM, 0, 0, B, A, I1, I2, I3, I4);
}

template <typename T> __global__ void hipTemplatePermute23(T *B, T *A, int I1, int I2, int I3, int I4)
{            
    int idx = threadIdx.x + blockIdx.x * blockDim.x; // global thread index  
    int M = I1*I2;
    int N = M*I3;
    int P = N*I4;
    while (idx<P) {
        int n = idx%N;
        int i4 = (idx-n)/N;        
        int m = n%M;
        int i3 = (n-m)/M;      
        int i1 = m%I1;
        int i2 = (m-i1)/I1;              
        B[i1+I1*i3+I1*I3*i2+N*i4] = A[idx];
        idx += blockDim.x * gridDim.x;
    }
}

template <typename T> void hipPermute23(T *B, T *A, int I1, int I2, int I3, int I4)
{                
    int BLOCKDIM = 256;
    int gridDim = (I1*I2*I3 + BLOCKDIM - 1) / BLOCKDIM;
    gridDim = (gridDim>1024)? 1024 : gridDim;
    hipLaunchKernelGGL(hipTemplatePermute23, gridDim, BLOCKDIM, 0, 0, B, A, I1, I2, I3, I4);
}

//void hipIndexPermute12(int*, int, int, int);
//void hipIndexPermute13(int*, int, int, int, int);
//void hipIndexPermute23(int*, int, int, int, int);

template void hipKron(double*, double*, double*, int, int);
template void hipKron(float*, float*, float*, int, int);
        
template void hipKron(double*, double*, double*, int, int, int, int);
template void hipKron(float*, float*, float*, int, int, int, int);

template void hipPermute(double*, double*, int*, int);
template void hipPermute(float*, float*, int*, int);
template void hipPermuteSharedMem(double*, double*, int*, int, int);
template void hipPermuteSharedMem(float*, float*, int*, int, int);

template void hipPermute12(double*, double*, int, int, int);
template void hipPermute12(float*, float*, int, int, int);
template void hipPermute13(double*, double*, int, int, int, int);
template void hipPermute13(float*, float*, int, int, int, int);
template void hipPermute23(double*, double*, int, int, int, int);
template void hipPermute23(float*, float*, int, int, int, int);

#endif
