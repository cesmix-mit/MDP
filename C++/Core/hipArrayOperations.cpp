/***************************************************************************                               
                    Molecular Dynamics Potentials (MDP)
                           CESMIX-MIT Project  
 
 Contributing authors: Ngoc-Cuong Nguyen (cuongng@mit.edu, exapde@gmail.com)
 ***************************************************************************/

#include "hip/hip_runtime.h"

#ifndef __HIPARRAYOPERATIONS
#define __HIPARRAYOPERATIONS

// template <typename T>
// __global__ void hipPrintArray2D(T* a, int m, int n)
// {        
//     for (int i=0; i<m; i++) {
//         for (int j=0; j<n; j++)
//             printf("%g   ", a[j*m+i]);                
//         printf("\n");
//     }
//     printf("\n");
// }

// template <typename T> void hipPrint2DArray(T* a, int m, int n)
// {
//     hipLaunchKernelGGL(hipPrintArray2D, 1, 1, 0, 0, a, m, n);
// }

// template <typename T>
// __global__ void hipPrintArray3D(T* a, int m, int n, int p)
// {    
//     for (int k=0; k<p; k++) {
//         for (int i=0; i<m; i++) {
//             for (int j=0; j<n; j++)
//                 printf("%g   ", a[k*n*m+j*m+i]);                
//             printf("\n");
//         }
//         printf("\n");
//     }
//     printf("\n");
// }

// template <typename T> void hipPrint3DArray(T* a, int m, int n, int p)
// {
//     hipLaunchKernelGGL(hipPrintArray3D, 1, 1, 0, 0, a, m, n, p);
// }

__global__ void hipKernelTripletnum(int *output, int *input, int n)
{
    int tid = threadIdx.x + blockIdx.x * blockDim.x;
    while (tid < n) {
        output[tid] = (input[tid]-1)*input[tid]/2;       
        tid += blockDim.x * gridDim.x;
    }
}
void hipTripletnum(int *output, int *input, int n)
{        
    int blockDim = 256;
    int gridDim = (n + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;
    hipLaunchKernelGGL(hipKernelTripletnum, gridDim, blockDim, 0, 0, output, input, n);
}

__global__ void hipKernelQuadrupletnum(int *output, int *input, int n)
{
    int tid = threadIdx.x + blockIdx.x * blockDim.x;
    while (tid < n) {
        output[tid] = (input[tid]-2)*(input[tid]-1)*input[tid]/6;       
        tid += blockDim.x * gridDim.x;
    }
}
void hipQuadrupletnum(int *output, int *input, int n)
{        
    int blockDim = 256;
    int gridDim = (n + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;
    hipLaunchKernelGGL(hipKernelQuadrupletnum, gridDim, blockDim, 0, 0, output, input, n);
}

__global__ void hipKernelIndexInit(int *ind, int n)
{
    int tid = threadIdx.x + blockIdx.x * blockDim.x;
    while (tid < n) {
        ind[tid] = tid;      
        tid += blockDim.x * gridDim.x;
    }
}
void hipIndexInit(int *ind, int n)
{        
    int blockDim = 256;
    int gridDim = (n + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;
    hipLaunchKernelGGL(hipKernelIndexInit, gridDim, blockDim, 0, 0, ind, n);
}

__global__ void hipKernelArrayFill(int* output, int start, int n) 
{	
    int tid = threadIdx.x + blockIdx.x * blockDim.x;
    while (tid < n) {
        output[tid] = start + tid;
        tid += blockDim.x * gridDim.x;
    }		
}
void hipArrayFill(int* output, int start, int n) 
{        
    int blockDim = 256;
    int gridDim = (n + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;
    hipLaunchKernelGGL(hipKernelArrayFill, gridDim, blockDim, 0, 0, output, start, n);
}
 
template <typename T>
__global__ void hipTemplateGetArrayAtIndex(T *y, T *x, int *ind, int n)
{
    int tid = threadIdx.x + blockIdx.x * blockDim.x;
    while (tid < n) {
        y[tid] = x[ind[tid]];      
        tid += blockDim.x * gridDim.x;
    }
}

template <typename T> void hipGetArrayAtIndex(T *y, T *x, int *ind, int n)
{        
    int blockDim = 256;
    int gridDim = (n + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;
    hipLaunchKernelGGL(hipTemplateGetArrayAtIndex, gridDim, blockDim, 0, 0, y, x, ind, n);
}

template <typename T>
__global__ void hipTemplatePutArrayAtIndex(T *y, T *x, int *ind, int n)
{
    int tid = threadIdx.x + blockIdx.x * blockDim.x;
    while (tid < n) {        
        y[ind[tid]] = x[tid];     
        tid += blockDim.x * gridDim.x;
    }
}

template <typename T> void hipPutArrayAtIndex(T *y, T *x, int *ind, int n)
{        
    int blockDim = 256;
    int gridDim = (n + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;
    hipLaunchKernelGGL(hipTemplatePutArrayAtIndex, gridDim, blockDim, 0, 0, y, x, ind, n);
}


template <typename T>
__global__ void hipTemplateArrayAXPYAtIndex(T *y, T *x, T a, int *ind, int n)
{
    int tid = threadIdx.x + blockIdx.x * blockDim.x;
    while (tid < n) {        
        y[ind[tid]] += a*x[tid];     
        tid += blockDim.x * gridDim.x;
    }
}

template <typename T> void hipArrayAXPYAtIndex(T *y, T *x, T a, int *ind, int n)
{        
    int blockDim = 256;
    int gridDim = (n + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;
    hipLaunchKernelGGL(hipTemplateArrayAXPYAtIndex, gridDim, blockDim, 0, 0, y, x, a, ind, n);
}

template <typename T>
__global__ void hipTemplateArrayPlusXAtIndex(T *y, T *x, int *ind, int n)
{
    int tid = threadIdx.x + blockIdx.x * blockDim.x;
    while (tid < n) {        
        y[ind[tid]] += x[tid];     
        tid += blockDim.x * gridDim.x;
    }
}

template <typename T> void hipArrayPlusXAtIndex(T *y, T *x, int *ind, int n)
{        
    int blockDim = 256;
    int gridDim = (n + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;
    hipLaunchKernelGGL(hipTemplateArrayPlusXAtIndex, gridDim, blockDim, 0, 0, y, x, ind, n);
}

template <typename T>
__global__ void hipTemplateArrayMinusXAtIndex(T *y, T *x, int *ind, int n)
{
    int tid = threadIdx.x + blockIdx.x * blockDim.x;
    while (tid < n) {        
        y[ind[tid]] -= x[tid];     
        tid += blockDim.x * gridDim.x;
    }
}

template <typename T> void hipArrayMinusXAtIndex(T *y, T *x, int *ind, int n)
{        
    int blockDim = 256;
    int gridDim = (n + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;
    hipLaunchKernelGGL(hipTemplateArrayMinusXAtIndex, gridDim, blockDim, 0, 0, y, x, ind, n);
}

template <typename T>
__global__ void hipTemplateArraySetValueAtIndex(T *y, T a, int n)
{    
    y[n] = a;          
}

template <typename T> void hipArraySetValueAtIndex(T *y, T a, int n)
{            
    hipLaunchKernelGGL(hipTemplateArraySetValueAtIndex, 1, 1, 0, 0, y, a, n);
}


template <typename T> T hipArrayGetValueAtIndex(T *y, int n)
{          
    T val; 
    hipMemcpy(&val, &y[n], sizeof(T), hipMemcpyDeviceToHost); 
    return val;
}

template <typename T>
__global__ void hipTemplateArraySetValue(T *y, T a, int n)
{
    int tid = threadIdx.x + blockIdx.x * blockDim.x;
    while (tid < n) {
        y[tid] = a;      
        tid += blockDim.x * gridDim.x;
    }
}

template <typename T> void hipArraySetValue(T *y, T a, int n)
{        
    int blockDim = 256;
    int gridDim = (n + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;
    hipLaunchKernelGGL(hipTemplateArraySetValue, gridDim, blockDim, 0, 0, y, a, n);
}

template <typename T>
__global__ void hipTemplateArrayMultiplyScalar(T *y, T a, int n)
{
    int tid = threadIdx.x + blockIdx.x * blockDim.x;
    while (tid < n) {
        y[tid] = a*y[tid];      
        tid += blockDim.x * gridDim.x;
    }
}

template <typename T> void hipArrayMultiplyScalar(T *y, T a, int n)
{        
    int blockDim = 256;
    int gridDim = (n + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;
    hipLaunchKernelGGL(hipTemplateArrayMultiplyScalar, gridDim, blockDim, 0, 0, y, a, n);
}

template <typename T>
__global__ void hipTemplateArrayAddScalar(T *y, T a, int n)
{
    int tid = threadIdx.x + blockIdx.x * blockDim.x;
    while (tid < n) {
        y[tid] += a;      
        tid += blockDim.x * gridDim.x;
    }
}

template <typename T> void hipArrayAddScalar(T *y, T a, int n)
{        
    int blockDim = 256;
    int gridDim = (n + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;
    hipLaunchKernelGGL(hipTemplateArrayAddScalar, gridDim, blockDim, 0, 0, y, a, n);
}

template <typename T>
__global__ void hipTemplateArrayCopy(T *y, T *x, int n)
{
    int tid = threadIdx.x + blockIdx.x * blockDim.x;
    while (tid < n) {
        y[tid] = x[tid];      
        tid += blockDim.x * gridDim.x;
    }
}

template <typename T> void hipArrayCopy(T *y, T *x, int n)
{        
    int blockDim = 256;
    int gridDim = (n + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;
    hipLaunchKernelGGL(hipTemplateArrayCopy, gridDim, blockDim, 0, 0, y, x, n);
}

template <typename T>
__global__ void hipTemplateArrayMinus(T *y, T *x, int n)
{
    int tid = threadIdx.x + blockIdx.x * blockDim.x;
    while (tid < n) {
        y[tid] = -x[tid];      
        tid += blockDim.x * gridDim.x;
    }
}

template <typename T> void hipArrayMinus(T *y, T *x, int n)
{        
    int blockDim = 256;
    int gridDim = (n + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;
    hipLaunchKernelGGL(hipTemplateArrayMinus, gridDim, blockDim, 0, 0, y, x, n);
}

template <typename T>
__global__ void hipTemplateArrayAbs(T *y, T *x, int n)
{
    int tid = threadIdx.x + blockIdx.x * blockDim.x;
    while (tid < n) {
        y[tid] = fabs(x[tid]);      
        tid += blockDim.x * gridDim.x;
    }
}

template <typename T> void hipArrayAbs(T *y, T *x, int n)
{        
    int blockDim = 256;
    int gridDim = (n + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;
    hipLaunchKernelGGL(hipTemplateArrayAbs, gridDim, blockDim, 0, 0, y, x, n);
}

template <typename T>
__global__ void hipTemplateArraySqrt(T *y, T *x, int n)
{
    int tid = threadIdx.x + blockIdx.x * blockDim.x;
    while (tid < n) {
        y[tid] = sqrt(x[tid]);      
        tid += blockDim.x * gridDim.x;
    }
}

template <typename T> void hipArraySqrt(T *y, T *x, int n)
{        
    int blockDim = 256;
    int gridDim = (n + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;
    hipLaunchKernelGGL(hipTemplateArraySqrt, gridDim, blockDim, 0, 0, y, x, n);
}

template <typename T>
__global__ void hipTemplateArraySin(T *y, T *x, int n)
{
    int tid = threadIdx.x + blockIdx.x * blockDim.x;
    while (tid < n) {
        y[tid] = sin(x[tid]);      
        tid += blockDim.x * gridDim.x;
    }
}

template <typename T> void hipArraySin(T *y, T *x, int n)
{        
    int blockDim = 256;
    int gridDim = (n + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;
    hipLaunchKernelGGL(hipTemplateArraySin, gridDim, blockDim, 0, 0, y, x, n);
}


template <typename T>
__global__ void hipTemplateArrayCos(T *y, T *x, int n)
{
    int tid = threadIdx.x + blockIdx.x * blockDim.x;
    while (tid < n) {
        y[tid] = cos(x[tid]);      
        tid += blockDim.x * gridDim.x;
    }
}

template <typename T> void hipArrayCos(T *y, T *x, int n)
{        
    int blockDim = 256;
    int gridDim = (n + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;
    hipLaunchKernelGGL(hipTemplateArrayCos, gridDim, blockDim, 0, 0, y, x, n);
}


template <typename T>
__global__ void hipTemplateArrayTan(T *y, T *x, int n)
{
    int tid = threadIdx.x + blockIdx.x * blockDim.x;
    while (tid < n) {
        y[tid] = tan(x[tid]);      
        tid += blockDim.x * gridDim.x;
    }
}

template <typename T> void hipArrayTan(T *y, T *x, int n)
{        
    int blockDim = 256;
    int gridDim = (n + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;
    hipLaunchKernelGGL(hipTemplateArrayTan, gridDim, blockDim, 0, 0, y, x, n);
}


template <typename T>
__global__ void hipTemplateArrayAsin(T *y, T *x, int n)
{
    int tid = threadIdx.x + blockIdx.x * blockDim.x;
    while (tid < n) {
        y[tid] = asin(x[tid]);      
        tid += blockDim.x * gridDim.x;
    }
}

template <typename T> void hipArrayAsin(T *y, T *x, int n)
{        
    int blockDim = 256;
    int gridDim = (n + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;
    hipLaunchKernelGGL(hipTemplateArrayAsin, gridDim, blockDim, 0, 0, y, x, n);
}


template <typename T>
__global__ void hipTemplateArrayAcos(T *y, T *x, int n)
{
    int tid = threadIdx.x + blockIdx.x * blockDim.x;
    while (tid < n) {
        y[tid] = acos(x[tid]);      
        tid += blockDim.x * gridDim.x;
    }
}

template <typename T> void hipArrayAcos(T *y, T *x, int n)
{        
    int blockDim = 256;
    int gridDim = (n + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;
    hipLaunchKernelGGL(hipTemplateArrayAcos, gridDim, blockDim, 0, 0, y, x, n);
}

template <typename T>
__global__ void hipTemplateArrayAtan(T *y, T *x, int n)
{
    int tid = threadIdx.x + blockIdx.x * blockDim.x;
    while (tid < n) {
        y[tid] = atan(x[tid]);      
        tid += blockDim.x * gridDim.x;
    }
}

template <typename T> void hipArrayAtan(T *y, T *x, int n)
{        
    int blockDim = 256;
    int gridDim = (n + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;
    hipLaunchKernelGGL(hipTemplateArrayAtan, gridDim, blockDim, 0, 0, y, x, n);
}

template <typename T>
__global__ void hipTemplateArraySinh(T *y, T *x, int n)
{
    int tid = threadIdx.x + blockIdx.x * blockDim.x;
    while (tid < n) {
        y[tid] = sinh(x[tid]);      
        tid += blockDim.x * gridDim.x;
    }
}

template <typename T> void hipArraySinh(T *y, T *x, int n)
{        
    int blockDim = 256;
    int gridDim = (n + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;
    hipLaunchKernelGGL(hipTemplateArraySinh, gridDim, blockDim, 0, 0, y, x, n);
}


template <typename T>
__global__ void hipTemplateArrayCosh(T *y, T *x, int n)
{
    int tid = threadIdx.x + blockIdx.x * blockDim.x;
    while (tid < n) {
        y[tid] = cosh(x[tid]);      
        tid += blockDim.x * gridDim.x;
    }
}

template <typename T> void hipArrayCosh(T *y, T *x, int n)
{        
    int blockDim = 256;
    int gridDim = (n + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;
    hipLaunchKernelGGL(hipTemplateArrayCosh, gridDim, blockDim, 0, 0, y, x, n);
}


template <typename T>
__global__ void hipTemplateArrayTanh(T *y, T *x, int n)
{
    int tid = threadIdx.x + blockIdx.x * blockDim.x;
    while (tid < n) {
        y[tid] = tanh(x[tid]);      
        tid += blockDim.x * gridDim.x;
    }
}

template <typename T> void hipArrayTanh(T *y, T *x, int n)
{        
    int blockDim = 256;
    int gridDim = (n + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;
    hipLaunchKernelGGL(hipTemplateArrayTanh, gridDim, blockDim, 0, 0, y, x, n);
}


template <typename T>
__global__ void hipTemplateArrayAsinh(T *y, T *x, int n)
{
    int tid = threadIdx.x + blockIdx.x * blockDim.x;
    while (tid < n) {
        y[tid] = asinh(x[tid]);      
        tid += blockDim.x * gridDim.x;
    }
}

template <typename T> void hipArrayAsinh(T *y, T *x, int n)
{        
    int blockDim = 256;
    int gridDim = (n + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;
    hipLaunchKernelGGL(hipTemplateArrayAsinh, gridDim, blockDim, 0, 0, y, x, n);
}


template <typename T>
__global__ void hipTemplateArrayAcosh(T *y, T *x, int n)
{
    int tid = threadIdx.x + blockIdx.x * blockDim.x;
    while (tid < n) {
        y[tid] = acosh(x[tid]);      
        tid += blockDim.x * gridDim.x;
    }
}

template <typename T> void hipArrayAcosh(T *y, T *x, int n)
{        
    int blockDim = 256;
    int gridDim = (n + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;
    hipLaunchKernelGGL(hipTemplateArrayAcosh, gridDim, blockDim, 0, 0, y, x, n);
}

template <typename T>
__global__ void hipTemplateArrayAtanh(T *y, T *x, int n)
{
    int tid = threadIdx.x + blockIdx.x * blockDim.x;
    while (tid < n) {
        y[tid] = atanh(x[tid]);      
        tid += blockDim.x * gridDim.x;
    }
}

template <typename T> void hipArrayAtanh(T *y, T *x, int n)
{        
    int blockDim = 256;
    int gridDim = (n + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;
    hipLaunchKernelGGL(hipTemplateArrayAtanh, gridDim, blockDim, 0, 0, y, x, n);
}

template <typename T>
__global__ void hipTemplateArrayExp(T *y, T *x, int n)
{
    int tid = threadIdx.x + blockIdx.x * blockDim.x;
    while (tid < n) {
        y[tid] = exp(x[tid]);      
        tid += blockDim.x * gridDim.x;
    }
}

template <typename T> void hipArrayExp(T *y, T *x, int n)
{        
    int blockDim = 256;
    int gridDim = (n + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;
    hipLaunchKernelGGL(hipTemplateArrayExp, gridDim, blockDim, 0, 0, y, x, n);
}

template <typename T>
__global__ void hipTemplateArrayLog(T *y, T *x, int n)
{
    int tid = threadIdx.x + blockIdx.x * blockDim.x;
    while (tid < n) {
        y[tid] = log(x[tid]);      
        tid += blockDim.x * gridDim.x;
    }
}

template <typename T> void hipArrayLog(T *y, T *x, int n)
{        
    int blockDim = 256;
    int gridDim = (n + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;
    hipLaunchKernelGGL(hipTemplateArrayLog, gridDim, blockDim, 0, 0, y, x, n);
}

template <typename T>
__global__ void hipTemplateArrayCeil(T *y, T *x, int n)
{
    int tid = threadIdx.x + blockIdx.x * blockDim.x;
    while (tid < n) {
        y[tid] = ceil(x[tid]);      
        tid += blockDim.x * gridDim.x;
    }
}

template <typename T> void hipArrayCeil(T *y, T *x, int n)
{        
    int blockDim = 256;
    int gridDim = (n + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;
    hipLaunchKernelGGL(hipTemplateArrayCeil, gridDim, blockDim, 0, 0, y, x, n);
}

template <typename T>
__global__ void hipTemplateArrayFloor(T *y, T *x, int n)
{
    int tid = threadIdx.x + blockIdx.x * blockDim.x;
    while (tid < n) {
        y[tid] = floor(x[tid]);      
        tid += blockDim.x * gridDim.x;
    }
}

template <typename T> void hipArrayFloor(T *y, T *x, int n)
{        
    int blockDim = 256;
    int gridDim = (n + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;
    hipLaunchKernelGGL(hipTemplateArrayFloor, gridDim, blockDim, 0, 0, y, x, n);
}

template <typename T>
__global__ void hipTemplateArrayErf(T *y, T *x, int n)
{
    int tid = threadIdx.x + blockIdx.x * blockDim.x;
    while (tid < n) {
        y[tid] = erf(x[tid]);      
        tid += blockDim.x * gridDim.x;
    }
}

template <typename T> void hipArrayErf(T *y, T *x, int n)
{        
    int blockDim = 256;
    int gridDim = (n + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;
    hipLaunchKernelGGL(hipTemplateArrayErf, gridDim, blockDim, 0, 0, y, x, n);
}

template <typename T>
__global__ void hipTemplateArrayErfc(T *y, T *x, int n)
{
    int tid = threadIdx.x + blockIdx.x * blockDim.x;
    while (tid < n) {
        y[tid] = erfc(x[tid]);      
        tid += blockDim.x * gridDim.x;
    }
}

template <typename T> void hipArrayErfc(T *y, T *x, int n)
{        
    int blockDim = 256;
    int gridDim = (n + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;
    hipLaunchKernelGGL(hipTemplateArrayErfc, gridDim, blockDim, 0, 0, y, x, n);
}

template <typename T>
__global__ void hipTemplateArraySquare(T *y, T *x, int n)
{
    int tid = threadIdx.x + blockIdx.x * blockDim.x;
    while (tid < n) {
        y[tid] = x[tid]*x[tid];      
        tid += blockDim.x * gridDim.x;
    }
}

template <typename T> void hipArraySquare(T *y, T *x, int n)
{        
    int blockDim = 256;
    int gridDim = (n + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;
    hipLaunchKernelGGL(hipTemplateArraySquare, gridDim, blockDim, 0, 0, y, x, n);
}

template <typename T>
__global__ void hipTemplateArrayPower(T *y, T *x, int p, int n)
{
    int tid = threadIdx.x + blockIdx.x * blockDim.x;
    while (tid < n) {
        y[tid] = x[tid];      
        for (int j=1; j<p; j++)
            y[tid] = y[tid]*x[tid];
        tid += blockDim.x * gridDim.x;
    }
}

template <typename T> void hipArrayPower(T *y, T *x, int p, int n)
{        
    int blockDim = 256;
    int gridDim = (n + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;
    hipLaunchKernelGGL(hipTemplateArrayPower, gridDim, blockDim, 0, 0, y, x, p, n);
}


template <typename T>
__global__ void hipTemplateArrayMultiplyScalarDiagonal(T *C, T a, int n)
{        
    int tid = threadIdx.x + blockIdx.x * blockDim.x;
    while (tid < n) {    
        C[tid+n*tid] = a*C[tid+n*tid];         
        tid += blockDim.x * gridDim.x;       
    }
}

template <typename T> void hipArrayMultiplyScalarDiagonal(T *C, T a, int n)
{        
    int blockDim = 256;
    int gridDim = (n + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;
    hipLaunchKernelGGL(hipTemplateArrayMultiplyScalarDiagonal, gridDim, blockDim, 0, 0, C, a, n);
}


template <typename T>
__global__ void hipTemplateArrayAddVectorToDiagonal(T *C, T *x, T a, int n)
{        
    int tid = threadIdx.x + blockIdx.x * blockDim.x;
    while (tid < n) {    
        C[tid+n*tid] += a*x[tid];         
        tid += blockDim.x * gridDim.x;       
    }
}

template <typename T> void hipArrayAddVectorToDiagonal(T *C, T *x, T a, int n)
{        
    int blockDim = 256;
    int gridDim = (n + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;
    hipLaunchKernelGGL(hipTemplateArrayAddVectorToDiagonal, gridDim, blockDim, 0, 0, C, x, a, n);
}


template <typename T>
__global__ void hipTemplateArrayRowAverage(T *y, T *x, int m, int n)
{    
    int j = threadIdx.x + blockIdx.x * blockDim.x;    
    while (j < n) {        
        T avg = 0;
        int i;
        for (i=0; i<m; i++)
            avg = avg + x[i + m*j];
        avg = avg/((T) m);
        for (i=0; i<m; i++)
            y[i + m*j] = avg;         
        j += blockDim.x * gridDim.x;
    }        
}

template <typename T> void hipArrayRowAverage(T *y, T *x, int m, int n)
{        
    int blockDim = 256;
    int gridDim = (n + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;
    hipLaunchKernelGGL(hipTemplateArrayRowAverage, gridDim, blockDim, 0, 0, y, x, m, n);
}


template <typename T>
__global__ void hipTemplateArrayAXPB(T *y, T *x, T a, T b, int n)
{
    int tid = threadIdx.x + blockIdx.x * blockDim.x;
    while (tid < n) {
        y[tid] = a*x[tid]+b;      
        tid += blockDim.x * gridDim.x;
    }
}

template <typename T> void hipArrayAXPB(T *y, T *x, T a, T b, int n)
{        
    int blockDim = 256;
    int gridDim = (n + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;
    hipLaunchKernelGGL(hipTemplateArrayAXPB, gridDim, blockDim, 0, 0, y, x, a, b, n);
}


template <typename T>
__global__ void hipTemplateArrayAXPBY(T *z, T *x, T *y, T a, T b, int n)
{
    int tid = threadIdx.x + blockIdx.x * blockDim.x;
    while (tid < n) {
        z[tid] = a*x[tid]+b*y[tid];      
        tid += blockDim.x * gridDim.x;
    }
}

template <typename T> void hipArrayAXPBY(T *z, T *x, T *y, T a, T b, int n)
{        
    int blockDim = 256;
    int gridDim = (n + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;
    hipLaunchKernelGGL(hipTemplateArrayAXPBY, gridDim, blockDim, 0, 0, z, x, y, a, b, n);
}

template <typename T>
__global__ void hipTemplateArrayAXY(T *s, T *x, T *y, T a, int n)
{
    int tid = threadIdx.x + blockIdx.x * blockDim.x;
    while (tid < n) {
        s[tid] = a*x[tid]*y[tid];      
        tid += blockDim.x * gridDim.x;
    }
}

template <typename T> void hipArrayAXY(T *s, T *x, T *y, T a, int n)
{        
    int blockDim = 256;
    int gridDim = (n + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;
    hipLaunchKernelGGL(hipTemplateArrayAXY, gridDim, blockDim, 0, 0, s, x, y, a, n);
}

template <typename T>
__global__ void hipTemplateArrayAXYZ(T *s, T *x, T *y, T *z, T a, int n)
{
    int tid = threadIdx.x + blockIdx.x * blockDim.x;
    while (tid < n) {
        s[tid] = a*x[tid]*y[tid]*z[tid];      
        tid += blockDim.x * gridDim.x;
    }
}

template <typename T> void hipArrayAXYZ(T *s, T *x, T *y, T *z, T a, int n)
{        
    int blockDim = 256;
    int gridDim = (n + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;
    hipLaunchKernelGGL(hipTemplateArrayAXYZ, gridDim, blockDim, 0, 0, s, x, y, z, a, n);
}

template <typename T>
__global__ void hipTemplateArrayAXYPBZ(T *s, T *x, T *y, T *z, T a, T b, int n)
{
    int tid = threadIdx.x + blockIdx.x * blockDim.x;
    while (tid < n) {
        s[tid] = a*x[tid]*y[tid] + b*z[tid];      
        tid += blockDim.x * gridDim.x;
    }
}

template <typename T> void hipArrayAXYPBZ(T *s, T *x, T *y, T *z, T a, T b, int n)
{        
    int blockDim = 256;
    int gridDim = (n + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;
    hipLaunchKernelGGL(hipTemplateArrayAXYPBZ, gridDim, blockDim, 0, 0, s, x, y, z, a, b, n);
}

template <typename T>
__global__ void hipTemplateArrayAdd3Vectors(T *s, T *x, T *y, T *z, T a, T b, T c, int n)
{
    int tid = threadIdx.x + blockIdx.x * blockDim.x;
    while (tid < n) {
        s[tid] = a*x[tid] + b*y[tid] + c*z[tid];      
        tid += blockDim.x * gridDim.x;
    }
}

template <typename T> void hipArrayAdd3Vectors(T *s, T *x, T *y, T *z, T a, T b, T c, int n)
{        
    int blockDim = 256;
    int gridDim = (n + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;
    hipLaunchKernelGGL(hipTemplateArrayAdd3Vectors, gridDim, blockDim, 0, 0, s, x, y, z, a, b, c, n);
}

template <typename T>
__global__ void  hipTemplateArrayAdd3Vector(T *a, T *b, T *c, T *d, int n)
{
    int tid = threadIdx.x + blockIdx.x * blockDim.x;
    while (tid < n) {        
        a[tid] = b[tid] + c[tid] + d[tid];
        tid += blockDim.x * gridDim.x;
    }
}

template <typename T> void hipArrayAdd3Vector(T *a, T *b, T *c, T *d, int n)
{        
    int blockDim = 256;
    int gridDim = (n + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;
    hipLaunchKernelGGL(hipTemplateArrayAdd3Vector, gridDim, blockDim, 0, 0, a, b, c, d, n);
}

template <typename T>
__global__ void hipTemplateArrayExtract(T *un, T *u, int I, int J, int M, int N,  
        int i1, int j1, int k1, int ni)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;
    while (idx < N) {
        int l = idx%M;        
        int i = l%ni+i1;
        int j = (l-i)/ni+j1;
        int k = (idx-l)/M+k1;
        un[idx] = u[i+I*j+I*J*k];      
        idx += blockDim.x * gridDim.x;
    }
}

template <typename T> void hipArrayExtract(T *un, T *u, int I, int J, int K, 
        int i1, int i2, int j1, int j2, int k1, int k2)
{        
    int ni = i2-i1;
    int nj = j2-j1;
    int nk = k2-k1;    
    int M = ni*nj;
    int N = M*nk;
    int blockDim = 256;
    int gridDim = (N + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;
    hipLaunchKernelGGL(hipTemplateArrayExtract, gridDim, blockDim, 0, 0, un, u, I, J, M, N, i1, j1, k1, ni);
}

template <typename T>
__global__ void hipTemplateArrayInsert(T *u, T *un, int I, int J, int M, int N,  
        int i1, int j1, int k1, int ni)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;
    while (idx < N) {
        int l = idx%M;        
        int i = l%ni+i1;
        int j = (l-i)/ni+j1;
        int k = (idx-l)/M+k1;
        u[i+I*j+I*J*k] = un[idx];      
        idx += blockDim.x * gridDim.x;
    }
}

template <typename T> void hipArrayInsert(T *u, T *un, int I, int J, int K, 
        int i1, int i2, int j1, int j2, int k1, int k2)
{        
    int ni = i2-i1;
    int nj = j2-j1;
    int nk = k2-k1;    
    int M = ni*nj;
    int N = M*nk;
    int blockDim = 256;
    int gridDim = (N + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;
    hipLaunchKernelGGL(hipTemplateArrayInsert, gridDim, blockDim, 0, 0, u, un, I, J, M, N, i1, j1, k1, ni);
}

template <typename T> 
__global__ void hipTemplateArrayGemmSharedMem(T *C, T *A, T *B, int I, int J, int K, int N, int Q)
{        
    // static shared memory
    __shared__ T Ashared[256];

    if (threadIdx.x<Q)
    {
      // load data from global memory to shared memory
      Ashared[threadIdx.x] = A[threadIdx.x];
    }

    // thread synchronization
    __syncthreads();

    int idx = threadIdx.x + blockIdx.x * blockDim.x;
    while (idx < N) {        
        int i = idx%I;
        int j = (idx-i)/I;                
        int m = K*j;
        C[idx] = 0.0;
        for (int k=0; k<K; k++)
            C[idx] += Ashared[i+I*k]*B[k+m];
        idx += blockDim.x * gridDim.x;
    }            
}

template <typename T> void hipArrayGemmSharedMem(T *C, T *A, T *B, int I, int J, int K)
{        
    // C[I*J] = A[I*K] x B[K*J]
    int N = I*J;    
    int Q = I*K;
    int blockDim = 256;
    int gridDim = (N + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;
    hipLaunchKernelGGL(hipTemplateArrayGemmSharedMem, gridDim, blockDim, 0, 0, C, A, B, I, J, K, N, Q);
}

template <typename T> 
__global__ void hipTemplateArrayGemmBatch(T *C, T *A, T *B, int I, int J, int K, int S, int M, int N, int P, int Q)
{        
    int idx = threadIdx.x + blockIdx.x * blockDim.x;
    while (idx < N) {
        int l = idx%M;        
        int i = l%I;
        int j = (l-i)/I;
        int s = (idx-l)/M;
        int a = i+Q*s;
        int b = K*j+P*s;
        C[idx] = 0.0;
        for (int k=0; k<K; k++)
            C[idx] += A[a+I*k]*B[k+b];
        idx += blockDim.x * gridDim.x;
    }            
}

template <typename T> void hipArrayGemmBatch(T *C, T *A, T *B, int I, int J, int K, int S)
{        
    // C[I*J*S] = A[I*K*S] x B[K*J*S]
    int M = I*J;
    int N = M*S;
    int Q = I*K;
    int P = K*J;
    int blockDim = 256;
    int gridDim = (N + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;
    hipLaunchKernelGGL(hipTemplateArrayGemmBatch, gridDim, blockDim, 0, 0, C, A, B, I, J, K, S, M, N, P, Q);
}

template <typename T> 
__global__ void hipTemplateArrayGemmBatch1(T *C, T *A, T *B, int I, int J, int K, int S, int M, int N, int P, int Q)
{        
    int idx = threadIdx.x + blockIdx.x * blockDim.x;
    while (idx < N) {
        int l = idx%M;        
        int i = l%I;
        int j = (l-i)/I;
        int s = (idx-l)/M;
        int a = i+Q*s;
        int b = K*j+P*s;
        for (int k=0; k<K; k++)
            C[idx] += A[a+I*k]*B[k+b];
        idx += blockDim.x * gridDim.x;
    }            
}

template <typename T> void hipArrayGemmBatch1(T *C, T *A, T *B, int I, int J, int K, int S)
{        
    // C[I*J*S] = A[I*K*S] x B[K*J*S] + C[I*J*S]
    int M = I*J;
    int N = M*S;
    int Q = I*K;
    int P = K*J;
    int blockDim = 256;
    int gridDim = (N + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;
    hipLaunchKernelGGL(hipTemplateArrayGemmBatch1, gridDim, blockDim, 0, 0, C, A, B, I, J, K, S, M, N, P, Q);
}

template <typename T> 
__global__ void hipTemplateArrayDG2CG(T *ucg, T *udg, int *cgent2dgent, int *rowent2elem, int nent)
{        
    int i = threadIdx.x + blockIdx.x * blockDim.x;
    while (i < nent) {
        ucg[i] = 0.0;
        int nelem = rowent2elem[i+1]-rowent2elem[i];
        for (int k=0; k<nelem; k++)
            ucg[i] += udg[cgent2dgent[rowent2elem[i]+k]]; 
        ucg[i] = ucg[i]/((T) nelem);        
        i += blockDim.x * gridDim.x;
    }            
}

template <typename T> void hipArrayDG2CG(T *ucg, T *udg, int *cgent2dgent, int *rowent2elem, int nent)
{        
    int blockDim = 256;
    int gridDim = (nent + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;
    hipLaunchKernelGGL(hipTemplateArrayDG2CG, gridDim, blockDim, 0, 0, ucg, udg, cgent2dgent, rowent2elem, nent);
}

template <typename T> 
__global__ void hipTemplateArrayDG2CG2(T *ucg, T *udg, int *colent2elem, int *rowent2elem, int nent, int npe)
{        
    int i = threadIdx.x + blockIdx.x * blockDim.x;
    while (i < nent) {
        ucg[i] = 0.0;
        int nelem = rowent2elem[i+1]-rowent2elem[i];
        for (int k=0; k<nelem; k++) {
            int e = colent2elem[rowent2elem[i]+k];
            for (int j=0; j<npe; j++)
                ucg[i] += udg[j+npe*e]; 
        }
        ucg[i] = ucg[i]/((T) (nelem*npe));
        i += blockDim.x * gridDim.x;
    }            
}
template <typename T> void hipArrayDG2CG2(T *ucg, T *udg, int *colent2elem, int *rowent2elem, int nent, int npe)
{        
    int blockDim = 256;
    int gridDim = (nent + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;
    hipLaunchKernelGGL(hipTemplateArrayDG2CG2, gridDim, blockDim, 0, 0, ucg, udg, colent2elem, rowent2elem, nent, npe);
}

template <typename T>
__global__ void hipTemplateArrayABPXYZ(T *fij, T *c3ij, T *eij, T *d3ij, T *g3ij, int dim, int n)
{
    int i = threadIdx.x + blockIdx.x * blockDim.x;
    while (i < n) {        
        for (int j=0; j<dim; j++)
            fij[j+3*i] = fij[j+3*i]*c3ij[i] + eij[i]*d3ij[i]*g3ij[j+3*i];  
        i += blockDim.x * gridDim.x;
    }
}
template <typename T> void hipArrayABPXYZ(T *fij, T *c3ij, T *eij, T *d3ij, T *g3ij, int dim, int n)
{        
    int blockDim = 256;
    int gridDim = (n + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;
    hipLaunchKernelGGL(hipTemplateArrayABPXYZ, gridDim, blockDim, 0, 0, fij, c3ij, eij, d3ij, g3ij, dim, n);
}
template void hipArrayABPXYZ(double*, double*, double*, double*, double*, int, int);
template void hipArrayABPXYZ(float*, float*, float*, float*, float*, int, int);

// template void hipPrint2DArray(double*, int, int);
// template void hipPrint3DArray(double*, int, int, int);
template void hipArraySetValue(double*, double, int);
template void hipArraySetValueAtIndex(double*, double, int);
template double hipArrayGetValueAtIndex(double*, int);
template void hipArrayAddScalar(double*, double, int);
template void hipArrayMultiplyScalar(double*, double, int);

template void hipGetArrayAtIndex(double*, double*, int*, int);
template void hipPutArrayAtIndex(double*, double*, int*, int);
template void hipArrayAXPYAtIndex(double*, double*, double, int*, int);
template void hipArrayPlusXAtIndex(double*, double*, int*, int);
template void hipArrayMinusXAtIndex(double*, double*, int*, int);

template void hipArrayCopy(double*, double*, int);
template void hipArrayMinus(double*, double*, int);
template void hipArrayAbs(double*, double*, int);
template void hipArraySqrt(double*, double*, int);
template void hipArraySin(double*, double*, int);
template void hipArrayCos(double*, double*, int);
template void hipArrayTan(double*, double*, int);
template void hipArrayAsin(double*, double*, int);
template void hipArrayAcos(double*, double*, int);
template void hipArrayAtan(double*, double*, int);
template void hipArraySinh(double*, double*, int);
template void hipArrayCosh(double*, double*, int);
template void hipArrayTanh(double*, double*, int);
template void hipArrayAsinh(double*, double*, int);
template void hipArrayAcosh(double*, double*, int);
template void hipArrayAtanh(double*, double*, int);
template void hipArrayExp(double*, double*, int);
template void hipArrayLog(double*, double*, int);
template void hipArrayCeil(double*, double*, int);
template void hipArrayFloor(double*, double*, int);
template void hipArrayErf(double*, double*, int);
template void hipArrayErfc(double*, double*, int);
template void hipArraySquare(double*, double*, int);
template void hipArrayPower(double*, double*, int, int);

template void hipArrayMultiplyScalarDiagonal(double*, double, int);
template void hipArrayAddVectorToDiagonal(double*, double*, double, int);
template void hipArrayRowAverage(double*, double*, int, int);
template void hipArrayAXPB(double*, double*, double, double, int);
template void hipArrayAXPBY(double*, double*, double*, double, double, int);
template void hipArrayAXY(double*, double*, double*, double, int);
template void hipArrayAXYZ(double*, double*, double*, double*, double, int);
template void hipArrayAXYPBZ(double*, double*, double*, double*, double, double, int);
template void hipArrayAdd3Vectors(double*, double*, double*, double*, double, double, double, int);
template void hipArrayAdd3Vector(double*, double*, double*, double*, int);
template void hipArrayExtract(double*, double*, int, int, int, int, int, int, int, int, int);
template void hipArrayInsert(double*, double*, int, int, int, int, int, int, int, int, int);
template void hipArrayGemmSharedMem(double*, double*, double*, int, int, int);
template void hipArrayGemmBatch(double*, double*, double*, int, int, int, int);
template void hipArrayGemmBatch1(double*, double*, double*, int, int, int, int);
template void hipArrayDG2CG(double*, double*, int*, int*, int);
template void hipArrayDG2CG2(double*, double*, int*, int*, int, int);

// template void hipPrint2DArray(float*, int, int);
// template void hipPrint3DArray(float*, int, int, int);
template void hipArraySetValue(float*, float, int);
template void hipArraySetValueAtIndex(float*, float, int);
template float hipArrayGetValueAtIndex(float*, int);
template void hipArrayAddScalar(float*, float, int);
template void hipArrayMultiplyScalar(float*, float, int);

template void hipGetArrayAtIndex(float*, float*, int*, int);
template void hipPutArrayAtIndex(float*, float*, int*, int);
template void hipArrayAXPYAtIndex(float*, float*, float, int*, int);
template void hipArrayPlusXAtIndex(float*, float*, int*, int);
template void hipArrayMinusXAtIndex(float*, float*, int*, int);

template void hipArrayCopy(float*, float*, int);
template void hipArrayMinus(float*, float*, int);
template void hipArrayAbs(float*, float*, int);
template void hipArraySqrt(float*, float*, int);
template void hipArraySin(float*, float*, int);
template void hipArrayCos(float*, float*, int);
template void hipArrayTan(float*, float*, int);
template void hipArrayAsin(float*, float*, int);
template void hipArrayAcos(float*, float*, int);
template void hipArrayAtan(float*, float*, int);
template void hipArraySinh(float*, float*, int);
template void hipArrayCosh(float*, float*, int);
template void hipArrayTanh(float*, float*, int);
template void hipArrayAsinh(float*, float*, int);
template void hipArrayAcosh(float*, float*, int);
template void hipArrayAtanh(float*, float*, int);
template void hipArrayExp(float*, float*, int);
template void hipArrayLog(float*, float*, int);
template void hipArrayCeil(float*, float*, int);
template void hipArrayFloor(float*, float*, int);
template void hipArrayErf(float*, float*, int);
template void hipArrayErfc(float*, float*, int);
template void hipArraySquare(float*, float*, int);
template void hipArrayPower(float*, float*, int, int);

template void hipArrayMultiplyScalarDiagonal(float*, float, int);
template void hipArrayAddVectorToDiagonal(float*, float*, float, int);
template void hipArrayRowAverage(float*, float*, int, int);
template void hipArrayAXPB(float*, float*, float, float, int);
template void hipArrayAXPBY(float*, float*, float*, float, float, int);
template void hipArrayAXY(float*, float*, float*, float, int);
template void hipArrayAXYZ(float*, float*, float*, float*, float, int);
template void hipArrayAXYPBZ(float*, float*, float*, float*, float, float, int);
template void hipArrayAdd3Vectors(float*, float*, float*, float*, float, float, float, int);
template void hipArrayAdd3Vector(float*, float*, float*, float*, int);
template void hipArrayExtract(float*, float*, int, int, int, int, int, int, int, int, int);
template void hipArrayInsert(float*, float*, int, int, int, int, int, int, int, int, int);
template void hipArrayGemmSharedMem(float*, float*, float*, int, int, int);
template void hipArrayGemmBatch(float*, float*, float*, int, int, int, int);
template void hipArrayGemmBatch1(float*, float*, float*, int, int, int, int);
template void hipArrayDG2CG(float*, float*, int*, int*, int);
template void hipArrayDG2CG2(float*, float*, int*, int*, int, int);

template int hipArrayGetValueAtIndex(int*, int);
template void hipArrayCopy(int*, int*, int);

#endif


