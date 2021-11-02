/***************************************************************************                               
                    Molecular Dynamics Potentials (MDP)
                           CESMIX-MIT Project  
 
 Contributing authors: Ngoc-Cuong Nguyen (cuongng@mit.edu, exapde@gmail.com)
 ***************************************************************************/

#ifndef __GPUARRAYOPERATIONS
#define __GPUARRAYOPERATIONS

// template <typename T>
// __global__ void gpuPrintArray2D(T* a, int m, int n)
// {        
//     for (int i=0; i<m; i++) {
//         for (int j=0; j<n; j++)
//             printf("%g   ", a[j*m+i]);                
//         printf("\n");
//     }
//     printf("\n");
// }

// template <typename T> void gpuPrint2DArray(T* a, int m, int n)
// {
//     gpuPrintArray2D<<<1, 1>>>(a, m, n);
// }

// template <typename T>
// __global__ void gpuPrintArray3D(T* a, int m, int n, int p)
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

// template <typename T> void gpuPrint3DArray(T* a, int m, int n, int p)
// {
//     gpuPrintArray3D<<<1, 1>>>(a, m, n, p);
// }

__global__ void gpuKernelTripletnum(int *output, int *input, int n)
{
    int tid = threadIdx.x + blockIdx.x * blockDim.x;
    while (tid < n) {
        output[tid] = (input[tid]-1)*input[tid]/2;       
        tid += blockDim.x * gridDim.x;
    }
}
void gpuTripletnum(int *output, int *input, int n)
{        
    int blockDim = 256;
    int gridDim = (n + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;
    gpuKernelTripletnum<<<gridDim, blockDim>>>(output, input, n);
}

__global__ void gpuKernelQuadrupletnum(int *output, int *input, int n)
{
    int tid = threadIdx.x + blockIdx.x * blockDim.x;
    while (tid < n) {
        output[tid] = (input[tid]-2)*(input[tid]-1)*input[tid]/6;       
        tid += blockDim.x * gridDim.x;
    }
}
void gpuQuadrupletnum(int *output, int *input, int n)
{        
    int blockDim = 256;
    int gridDim = (n + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;
    gpuKernelQuadrupletnum<<<gridDim, blockDim>>>(output, input, n);
}

__global__ void gpuKernelIndexInit(int *ind, int n)
{
    int tid = threadIdx.x + blockIdx.x * blockDim.x;
    while (tid < n) {
        ind[tid] = tid;      
        tid += blockDim.x * gridDim.x;
    }
}
void gpuIndexInit(int *ind, int n)
{        
    int blockDim = 256;
    int gridDim = (n + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;
    gpuKernelIndexInit<<<gridDim, blockDim>>>(ind, n);
}

__global__ void gpuKernelArrayFill(int* output, int start, int n) 
{	
    int tid = threadIdx.x + blockIdx.x * blockDim.x;
    while (tid < n) {
        output[tid] = start + tid;
        tid += blockDim.x * gridDim.x;
    }		
}
void gpuArrayFill(int* output, int start, int n) 
{        
    int blockDim = 256;
    int gridDim = (n + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;
    gpuKernelArrayFill<<<gridDim, blockDim>>>(output, start, n);
}
 
template <typename T>
__global__ void gpuKernelGetArraySpacing(T *y, T *x, int m, int n)
{
    int tid = threadIdx.x + blockIdx.x * blockDim.x;
    while (tid < n) {
        y[tid] = x[m*tid];      
        tid += blockDim.x * gridDim.x;
    }
}
template <typename T> void gpuGetArraySpacing(T *y, T *x, int m, int n)
{        
    int blockDim = 256;
    int gridDim = (n + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;
    gpuKernelGetArraySpacing<<<gridDim, blockDim>>>(y, x, m, n);
}

template <typename T> 
__global__ void gpuKernelArrayTranspose(T *A, T *B, int m, int n, int k)
{  
    int ii = threadIdx.x + blockIdx.x * blockDim.x;  
    while (ii < k) {    
        int i = ii%m;        
        int j = (ii-i)/m;                
        A[j+n*i] = B[ii]; // ii = i + m*j
        ii += blockDim.x * gridDim.x;
    }
}
template <typename T> void gpuArrayTranspose(T *A, T *B, int m, int n)
{
    int blockDim = 256;
    int k = m*n;
    int gridDim = (k + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;
    gpuKernelArrayTranspose<<<gridDim, blockDim>>>(A, B, m, n, k);    
}

template <typename T> 
__global__ void gpuKernelArrayPlusAtColumnIndex(T *A, T *B, int *colind, int m, int n)
{  
    int ii = threadIdx.x + blockIdx.x * blockDim.x;  
    while (ii < n) {    
        int i = colind[ii];   
        for (int j=0; j<m; j++)
            A[ii*m+j] += B[i*m+j];
        ii += blockDim.x * gridDim.x;
    }
}
template <typename T> void gpuArrayPlusAtColumnIndex(T *A, T *B, int *colind, int m, int n)
{
    int blockDim = 256;
    int gridDim = (n + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;
    gpuKernelArrayPlusAtColumnIndex<<<gridDim, blockDim>>>(A, B, colind, m, n);    
}

template <typename T> 
__global__ void gpuKernelArrayMinusAtColumnIndex(T *A, T *B, int *colind, int m, int n)
{  
    int ii = threadIdx.x + blockIdx.x * blockDim.x;  
    while (ii < n) {    
        int i = colind[ii];   
        for (int j=0; j<m; j++)
            A[ii*m+j] = B[i*m+j] - A[ii*m+j];
        ii += blockDim.x * gridDim.x;
    }
}
template <typename T> void gpuArrayMinusAtColumnIndex(T *A, T *B, int *colind, int m, int n)
{
    int blockDim = 256;
    int gridDim = (n + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;
    gpuKernelArrayMinusAtColumnIndex<<<gridDim, blockDim>>>(A, B, colind, m, n);    
}

template <typename T> 
__global__ void gpuKernelGetArrayAtColumnIndex(T *A, T *B, int *colind, int m, int n)
{  
    int ii = threadIdx.x + blockIdx.x * blockDim.x;  
    while (ii < n) {    
        int i = colind[ii];   
        for (int j=0; j<m; j++)
            A[ii*m+j] = B[i*m+j];
        ii += blockDim.x * gridDim.x;
    }
}
template <typename T> void gpuGetArrayAtColumnIndex(T *A, T *B, int *colind, int m, int n)
{
    int blockDim = 256;
    int gridDim = (n + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;
    gpuKernelGetArrayAtColumnIndex<<<gridDim, blockDim>>>(A, B, colind, m, n);    
}

template <typename T> 
__global__ void gpuKernelArrayTransposeAtColumnIndex(T *A, T *B, int *colind, int m, int n)
{  
    int ii = threadIdx.x + blockIdx.x * blockDim.x;  
    while (ii < n) {    
        int i = colind[ii];   
        for (int j=0; j<m; j++)
            A[ii+j*n] = B[i*m+j];
        ii += blockDim.x * gridDim.x;
    }
}
template <typename T> void gpuArrayTransposeAtColumnIndex(T *A, T *B, int *colind, int m, int n)
{
    int blockDim = 256;
    int gridDim = (n + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;
    gpuKernelArrayTransposeAtColumnIndex<<<gridDim, blockDim>>>(A, B, colind, m, n);    
}

template <typename T> 
__global__ void gpuKernelGetArrayAtRowIndex(T *A, T *B, int *rowind, int m, int n)
{  
    int ii = threadIdx.x + blockIdx.x * blockDim.x;  
    while (ii < m) {    
        int i = rowind[ii];   
        for (int j=0; j<n; j++)
            A[j*m+ii] = B[j*m+i];
        ii += blockDim.x * gridDim.x;
    }
}
template <typename T> void gpuGetArrayAtRowIndex(T *A, T *B, int *rowind, int m, int n)
{
    int blockDim = 256;
    int gridDim = (m + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;
    gpuKernelGetArrayAtRowIndex<<<gridDim, blockDim>>>(A, B, rowind, m, n);    
}

template <typename T> 
__global__ void gpuKernelArrayTransposeAtRowIndex(T *A, T *B, int *rowind, int m, int n)
{  
    int ii = threadIdx.x + blockIdx.x * blockDim.x;  
    while (ii < m) {    
        int i = rowind[ii];   
        for (int j=0; j<n; j++)
            A[ii*n+j] = B[j*m+i];
        ii += blockDim.x * gridDim.x;
    }
}
template <typename T> void gpuArrayTransposeAtRowIndex(T *A, T *B, int *rowind, int m, int n)
{
    int blockDim = 256;
    int gridDim = (m + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;
    gpuKernelArrayTransposeAtRowIndex<<<gridDim, blockDim>>>(A, B, rowind, m, n);    
}

template <typename T>
__global__ void gpuKernelGetArrayAtIndex(T *y, T *x, int *ind, int n)
{
    int tid = threadIdx.x + blockIdx.x * blockDim.x;
    while (tid < n) {
        y[tid] = x[ind[tid]];      
        tid += blockDim.x * gridDim.x;
    }
}

template <typename T> void gpuGetArrayAtIndex(T *y, T *x, int *ind, int n)
{        
    int blockDim = 256;
    int gridDim = (n + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;
    gpuKernelGetArrayAtIndex<<<gridDim, blockDim>>>(y, x, ind, n);
}

template <typename T>
__global__ void gpuKernelPutArrayAtIndex(T *y, T *x, int *ind, int n)
{
    int tid = threadIdx.x + blockIdx.x * blockDim.x;
    while (tid < n) {        
        y[ind[tid]] = x[tid];     
        tid += blockDim.x * gridDim.x;
    }
}

template <typename T> void gpuPutArrayAtIndex(T *y, T *x, int *ind, int n)
{        
    int blockDim = 256;
    int gridDim = (n + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;
    gpuKernelPutArrayAtIndex<<<gridDim, blockDim>>>(y, x, ind, n);
}


template <typename T>
__global__ void gpuKernelArrayAXPYAtIndex(T *y, T *x, T a, int *ind, int n)
{
    int tid = threadIdx.x + blockIdx.x * blockDim.x;
    while (tid < n) {        
        y[ind[tid]] += a*x[tid];     
        tid += blockDim.x * gridDim.x;
    }
}

template <typename T> void gpuArrayAXPYAtIndex(T *y, T *x, T a, int *ind, int n)
{        
    int blockDim = 256;
    int gridDim = (n + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;
    gpuKernelArrayAXPYAtIndex<<<gridDim, blockDim>>>(y, x, a, ind, n);
}

template <typename T>
__global__ void gpuKernelArrayPlusXAtIndex(T *y, T *x, int *ind, int n)
{
    int tid = threadIdx.x + blockIdx.x * blockDim.x;
    while (tid < n) {        
        y[ind[tid]] += x[tid];     
        tid += blockDim.x * gridDim.x;
    }
}

template <typename T> void gpuArrayPlusXAtIndex(T *y, T *x, int *ind, int n)
{        
    int blockDim = 256;
    int gridDim = (n + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;
    gpuKernelArrayPlusXAtIndex<<<gridDim, blockDim>>>(y, x, ind, n);
}

template <typename T>
__global__ void gpuKernelArrayMinusXAtIndex(T *y, T *x, int *ind, int n)
{
    int tid = threadIdx.x + blockIdx.x * blockDim.x;
    while (tid < n) {        
        y[ind[tid]] -= x[tid];     
        tid += blockDim.x * gridDim.x;
    }
}

template <typename T> void gpuArrayMinusXAtIndex(T *y, T *x, int *ind, int n)
{        
    int blockDim = 256;
    int gridDim = (n + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;
    gpuKernelArrayMinusXAtIndex<<<gridDim, blockDim>>>(y, x, ind, n);
}

template <typename T>
__global__ void gpuKernelArraySetValueAtIndex(T *y, T a, int n)
{    
    y[n] = a;          
}

template <typename T> void gpuArraySetValueAtIndex(T *y, T a, int n)
{            
    gpuKernelArraySetValueAtIndex<<<1, 1>>>(y, a, n);
}


template <typename T> T gpuArrayGetValueAtIndex(T *y, int n)
{          
    T val; 
    cudaMemcpy(&val, &y[n], sizeof(T), cudaMemcpyDeviceToHost); 
    return val;
}

template <typename T>
__global__ void gpuKernelArraySetValue(T *y, T a, int n)
{
    int tid = threadIdx.x + blockIdx.x * blockDim.x;
    while (tid < n) {
        y[tid] = a;      
        tid += blockDim.x * gridDim.x;
    }
}

template <typename T> void gpuArraySetValue(T *y, T a, int n)
{        
    int blockDim = 256;
    int gridDim = (n + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;
    gpuKernelArraySetValue<<<gridDim, blockDim>>>(y, a, n);
}

template <typename T>
__global__ void gpuKernelArrayMultiplyScalar(T *y, T a, int n)
{
    int tid = threadIdx.x + blockIdx.x * blockDim.x;
    while (tid < n) {
        y[tid] = a*y[tid];      
        tid += blockDim.x * gridDim.x;
    }
}

template <typename T> void gpuArrayMultiplyScalar(T *y, T a, int n)
{        
    int blockDim = 256;
    int gridDim = (n + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;
    gpuKernelArrayMultiplyScalar<<<gridDim, blockDim>>>(y, a, n);
}

template <typename T>
__global__ void gpuKernelArrayAddScalar(T *y, T a, int n)
{
    int tid = threadIdx.x + blockIdx.x * blockDim.x;
    while (tid < n) {
        y[tid] += a;      
        tid += blockDim.x * gridDim.x;
    }
}

template <typename T> void gpuArrayAddScalar(T *y, T a, int n)
{        
    int blockDim = 256;
    int gridDim = (n + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;
    gpuKernelArrayAddScalar<<<gridDim, blockDim>>>(y, a, n);
}

template <typename T>
__global__ void gpuKernelArrayCopy(T *y, T *x, int n)
{
    int tid = threadIdx.x + blockIdx.x * blockDim.x;
    while (tid < n) {
        y[tid] = x[tid];      
        tid += blockDim.x * gridDim.x;
    }
}

template <typename T> void gpuArrayCopy(T *y, T *x, int n)
{        
    int blockDim = 256;
    int gridDim = (n + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;
    gpuKernelArrayCopy<<<gridDim, blockDim>>>(y, x, n);
}

template <typename T>
__global__ void gpuKernelArrayMinus(T *y, T *x, int n)
{
    int tid = threadIdx.x + blockIdx.x * blockDim.x;
    while (tid < n) {
        y[tid] = -x[tid];      
        tid += blockDim.x * gridDim.x;
    }
}

template <typename T> void gpuArrayMinus(T *y, T *x, int n)
{        
    int blockDim = 256;
    int gridDim = (n + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;
    gpuKernelArrayMinus<<<gridDim, blockDim>>>(y, x, n);
}

template <typename T>
__global__ void gpuKernelArrayAbs(T *y, T *x, int n)
{
    int tid = threadIdx.x + blockIdx.x * blockDim.x;
    while (tid < n) {
        y[tid] = fabs(x[tid]);      
        tid += blockDim.x * gridDim.x;
    }
}

template <typename T> void gpuArrayAbs(T *y, T *x, int n)
{        
    int blockDim = 256;
    int gridDim = (n + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;
    gpuKernelArrayAbs<<<gridDim, blockDim>>>(y, x, n);
}

template <typename T>
__global__ void gpuKernelArraySqrt(T *y, T *x, int n)
{
    int tid = threadIdx.x + blockIdx.x * blockDim.x;
    while (tid < n) {
        y[tid] = sqrt(x[tid]);      
        tid += blockDim.x * gridDim.x;
    }
}

template <typename T> void gpuArraySqrt(T *y, T *x, int n)
{        
    int blockDim = 256;
    int gridDim = (n + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;
    gpuKernelArraySqrt<<<gridDim, blockDim>>>(y, x, n);
}

template <typename T>
__global__ void gpuKernelArraySin(T *y, T *x, int n)
{
    int tid = threadIdx.x + blockIdx.x * blockDim.x;
    while (tid < n) {
        y[tid] = sin(x[tid]);      
        tid += blockDim.x * gridDim.x;
    }
}

template <typename T> void gpuArraySin(T *y, T *x, int n)
{        
    int blockDim = 256;
    int gridDim = (n + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;
    gpuKernelArraySin<<<gridDim, blockDim>>>(y, x, n);
}


template <typename T>
__global__ void gpuKernelArrayCos(T *y, T *x, int n)
{
    int tid = threadIdx.x + blockIdx.x * blockDim.x;
    while (tid < n) {
        y[tid] = cos(x[tid]);      
        tid += blockDim.x * gridDim.x;
    }
}

template <typename T> void gpuArrayCos(T *y, T *x, int n)
{        
    int blockDim = 256;
    int gridDim = (n + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;
    gpuKernelArrayCos<<<gridDim, blockDim>>>(y, x, n);
}


template <typename T>
__global__ void gpuKernelArrayTan(T *y, T *x, int n)
{
    int tid = threadIdx.x + blockIdx.x * blockDim.x;
    while (tid < n) {
        y[tid] = tan(x[tid]);      
        tid += blockDim.x * gridDim.x;
    }
}

template <typename T> void gpuArrayTan(T *y, T *x, int n)
{        
    int blockDim = 256;
    int gridDim = (n + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;
    gpuKernelArrayTan<<<gridDim, blockDim>>>(y, x, n);
}


template <typename T>
__global__ void gpuKernelArrayAsin(T *y, T *x, int n)
{
    int tid = threadIdx.x + blockIdx.x * blockDim.x;
    while (tid < n) {
        y[tid] = asin(x[tid]);      
        tid += blockDim.x * gridDim.x;
    }
}

template <typename T> void gpuArrayAsin(T *y, T *x, int n)
{        
    int blockDim = 256;
    int gridDim = (n + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;
    gpuKernelArrayAsin<<<gridDim, blockDim>>>(y, x, n);
}


template <typename T>
__global__ void gpuKernelArrayAcos(T *y, T *x, int n)
{
    int tid = threadIdx.x + blockIdx.x * blockDim.x;
    while (tid < n) {
        y[tid] = acos(x[tid]);      
        tid += blockDim.x * gridDim.x;
    }
}

template <typename T> void gpuArrayAcos(T *y, T *x, int n)
{        
    int blockDim = 256;
    int gridDim = (n + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;
    gpuKernelArrayAcos<<<gridDim, blockDim>>>(y, x, n);
}

template <typename T>
__global__ void gpuKernelArrayAtan(T *y, T *x, int n)
{
    int tid = threadIdx.x + blockIdx.x * blockDim.x;
    while (tid < n) {
        y[tid] = atan(x[tid]);      
        tid += blockDim.x * gridDim.x;
    }
}

template <typename T> void gpuArrayAtan(T *y, T *x, int n)
{        
    int blockDim = 256;
    int gridDim = (n + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;
    gpuKernelArrayAtan<<<gridDim, blockDim>>>(y, x, n);
}

template <typename T>
__global__ void gpuKernelArraySinh(T *y, T *x, int n)
{
    int tid = threadIdx.x + blockIdx.x * blockDim.x;
    while (tid < n) {
        y[tid] = sinh(x[tid]);      
        tid += blockDim.x * gridDim.x;
    }
}

template <typename T> void gpuArraySinh(T *y, T *x, int n)
{        
    int blockDim = 256;
    int gridDim = (n + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;
    gpuKernelArraySinh<<<gridDim, blockDim>>>(y, x, n);
}


template <typename T>
__global__ void gpuKernelArrayCosh(T *y, T *x, int n)
{
    int tid = threadIdx.x + blockIdx.x * blockDim.x;
    while (tid < n) {
        y[tid] = cosh(x[tid]);      
        tid += blockDim.x * gridDim.x;
    }
}

template <typename T> void gpuArrayCosh(T *y, T *x, int n)
{        
    int blockDim = 256;
    int gridDim = (n + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;
    gpuKernelArrayCosh<<<gridDim, blockDim>>>(y, x, n);
}


template <typename T>
__global__ void gpuKernelArrayTanh(T *y, T *x, int n)
{
    int tid = threadIdx.x + blockIdx.x * blockDim.x;
    while (tid < n) {
        y[tid] = tanh(x[tid]);      
        tid += blockDim.x * gridDim.x;
    }
}

template <typename T> void gpuArrayTanh(T *y, T *x, int n)
{        
    int blockDim = 256;
    int gridDim = (n + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;
    gpuKernelArrayTanh<<<gridDim, blockDim>>>(y, x, n);
}


template <typename T>
__global__ void gpuKernelArrayAsinh(T *y, T *x, int n)
{
    int tid = threadIdx.x + blockIdx.x * blockDim.x;
    while (tid < n) {
        y[tid] = asinh(x[tid]);      
        tid += blockDim.x * gridDim.x;
    }
}

template <typename T> void gpuArrayAsinh(T *y, T *x, int n)
{        
    int blockDim = 256;
    int gridDim = (n + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;
    gpuKernelArrayAsinh<<<gridDim, blockDim>>>(y, x, n);
}


template <typename T>
__global__ void gpuKernelArrayAcosh(T *y, T *x, int n)
{
    int tid = threadIdx.x + blockIdx.x * blockDim.x;
    while (tid < n) {
        y[tid] = acosh(x[tid]);      
        tid += blockDim.x * gridDim.x;
    }
}

template <typename T> void gpuArrayAcosh(T *y, T *x, int n)
{        
    int blockDim = 256;
    int gridDim = (n + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;
    gpuKernelArrayAcosh<<<gridDim, blockDim>>>(y, x, n);
}

template <typename T>
__global__ void gpuKernelArrayAtanh(T *y, T *x, int n)
{
    int tid = threadIdx.x + blockIdx.x * blockDim.x;
    while (tid < n) {
        y[tid] = atanh(x[tid]);      
        tid += blockDim.x * gridDim.x;
    }
}

template <typename T> void gpuArrayAtanh(T *y, T *x, int n)
{        
    int blockDim = 256;
    int gridDim = (n + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;
    gpuKernelArrayAtanh<<<gridDim, blockDim>>>(y, x, n);
}

template <typename T>
__global__ void gpuKernelArrayExp(T *y, T *x, int n)
{
    int tid = threadIdx.x + blockIdx.x * blockDim.x;
    while (tid < n) {
        y[tid] = exp(x[tid]);      
        tid += blockDim.x * gridDim.x;
    }
}

template <typename T> void gpuArrayExp(T *y, T *x, int n)
{        
    int blockDim = 256;
    int gridDim = (n + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;
    gpuKernelArrayExp<<<gridDim, blockDim>>>(y, x, n);
}

template <typename T>
__global__ void gpuKernelArrayLog(T *y, T *x, int n)
{
    int tid = threadIdx.x + blockIdx.x * blockDim.x;
    while (tid < n) {
        y[tid] = log(x[tid]);      
        tid += blockDim.x * gridDim.x;
    }
}

template <typename T> void gpuArrayLog(T *y, T *x, int n)
{        
    int blockDim = 256;
    int gridDim = (n + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;
    gpuKernelArrayLog<<<gridDim, blockDim>>>(y, x, n);
}

template <typename T>
__global__ void gpuKernelArrayCeil(T *y, T *x, int n)
{
    int tid = threadIdx.x + blockIdx.x * blockDim.x;
    while (tid < n) {
        y[tid] = ceil(x[tid]);      
        tid += blockDim.x * gridDim.x;
    }
}

template <typename T> void gpuArrayCeil(T *y, T *x, int n)
{        
    int blockDim = 256;
    int gridDim = (n + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;
    gpuKernelArrayCeil<<<gridDim, blockDim>>>(y, x, n);
}

template <typename T>
__global__ void gpuKernelArrayFloor(T *y, T *x, int n)
{
    int tid = threadIdx.x + blockIdx.x * blockDim.x;
    while (tid < n) {
        y[tid] = floor(x[tid]);      
        tid += blockDim.x * gridDim.x;
    }
}

template <typename T> void gpuArrayFloor(T *y, T *x, int n)
{        
    int blockDim = 256;
    int gridDim = (n + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;
    gpuKernelArrayFloor<<<gridDim, blockDim>>>(y, x, n);
}

template <typename T>
__global__ void gpuKernelArrayErf(T *y, T *x, int n)
{
    int tid = threadIdx.x + blockIdx.x * blockDim.x;
    while (tid < n) {
        y[tid] = erf(x[tid]);      
        tid += blockDim.x * gridDim.x;
    }
}

template <typename T> void gpuArrayErf(T *y, T *x, int n)
{        
    int blockDim = 256;
    int gridDim = (n + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;
    gpuKernelArrayErf<<<gridDim, blockDim>>>(y, x, n);
}

template <typename T>
__global__ void gpuKernelArrayErfc(T *y, T *x, int n)
{
    int tid = threadIdx.x + blockIdx.x * blockDim.x;
    while (tid < n) {
        y[tid] = erfc(x[tid]);      
        tid += blockDim.x * gridDim.x;
    }
}

template <typename T> void gpuArrayErfc(T *y, T *x, int n)
{        
    int blockDim = 256;
    int gridDim = (n + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;
    gpuKernelArrayErfc<<<gridDim, blockDim>>>(y, x, n);
}

template <typename T>
__global__ void gpuKernelArraySquare(T *y, T *x, int n)
{
    int tid = threadIdx.x + blockIdx.x * blockDim.x;
    while (tid < n) {
        y[tid] = x[tid]*x[tid];      
        tid += blockDim.x * gridDim.x;
    }
}

template <typename T> void gpuArraySquare(T *y, T *x, int n)
{        
    int blockDim = 256;
    int gridDim = (n + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;
    gpuKernelArraySquare<<<gridDim, blockDim>>>(y, x, n);
}

template <typename T>
__global__ void gpuKernelArrayPower(T *y, T *x, int p, int n)
{
    int tid = threadIdx.x + blockIdx.x * blockDim.x;
    while (tid < n) {
        y[tid] = x[tid];      
        for (int j=1; j<p; j++)
            y[tid] = y[tid]*x[tid];
        tid += blockDim.x * gridDim.x;
    }
}

template <typename T> void gpuArrayPower(T *y, T *x, int p, int n)
{        
    int blockDim = 256;
    int gridDim = (n + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;
    gpuKernelArrayPower<<<gridDim, blockDim>>>(y, x, p, n);
}


template <typename T>
__global__ void gpuKernelArrayMultiplyScalarDiagonal(T *C, T a, int n)
{        
    int tid = threadIdx.x + blockIdx.x * blockDim.x;
    while (tid < n) {    
        C[tid+n*tid] = a*C[tid+n*tid];         
        tid += blockDim.x * gridDim.x;       
    }
}

template <typename T> void gpuArrayMultiplyScalarDiagonal(T *C, T a, int n)
{        
    int blockDim = 256;
    int gridDim = (n + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;
    gpuKernelArrayMultiplyScalarDiagonal<<<gridDim, blockDim>>>(C, a, n);
}


template <typename T>
__global__ void gpuKernelArrayAddVectorToDiagonal(T *C, T *x, T a, int n)
{        
    int tid = threadIdx.x + blockIdx.x * blockDim.x;
    while (tid < n) {    
        C[tid+n*tid] += a*x[tid];         
        tid += blockDim.x * gridDim.x;       
    }
}

template <typename T> void gpuArrayAddVectorToDiagonal(T *C, T *x, T a, int n)
{        
    int blockDim = 256;
    int gridDim = (n + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;
    gpuKernelArrayAddVectorToDiagonal<<<gridDim, blockDim>>>(C, x, a, n);
}


template <typename T>
__global__ void gpuKernelArrayRowAverage(T *y, T *x, int m, int n)
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

template <typename T> void gpuArrayRowAverage(T *y, T *x, int m, int n)
{        
    int blockDim = 256;
    int gridDim = (n + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;
    gpuKernelArrayRowAverage<<<gridDim, blockDim>>>(y, x, m, n);
}

template <typename T>
__global__ void gpuKernelArrayRowSum(T *y, T *x, int m, int n)
{    
    int j = threadIdx.x + blockIdx.x * blockDim.x;    
    while (j < n) {        
        y[j] = 0;        
        for (int i=0; i<m; i++)
            y[j] += x[i + m*j];
        j += blockDim.x * gridDim.x;
    }        
}

template <typename T> void gpuArrayRowSum(T *y, T *x, int m, int n)
{        
    int blockDim = 256;
    int gridDim = (n + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;
    gpuKernelArrayRowSum<<<gridDim, blockDim>>>(y, x, m, n);
}

template <typename T>
__global__ void gpuKernelArrayRowSquareSum(T *y, T *x, int m, int n)
{    
    int j = threadIdx.x + blockIdx.x * blockDim.x;    
    while (j < n) {        
        y[j] = 0;        
        for (int i=0; i<m; i++)
            y[j] += x[i + m*j]*x[i + m*j];
        j += blockDim.x * gridDim.x;
    }        
}

template <typename T> void gpuArrayRowSquareSum(T *y, T *x, int m, int n)
{        
    int blockDim = 256;
    int gridDim = (n + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;
    gpuKernelArrayRowSquareSum<<<gridDim, blockDim>>>(y, x, m, n);
}

template <typename T>
__global__ void gpuKernelArrayDistSquareSum(T *y, T *x1, T *x2, int m, int n)
{    
    int j = threadIdx.x + blockIdx.x * blockDim.x;    
    while (j < n) {        
        y[j] = 0;        
        for (int i=0; i<m; i++)
            y[j] += (x2[i + m*j]-x1[i + m*j])*(x2[i + m*j]-x1[i + m*j]);
        j += blockDim.x * gridDim.x;
    }        
}

template <typename T> void gpuArrayDistSquareSum(T *y, T *x1, T *x2, int m, int n)
{        
    int blockDim = 256;
    int gridDim = (n + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;
    gpuKernelArrayDistSquareSum<<<gridDim, blockDim>>>(y, x1, x2, m, n);
}

template <typename T>
__global__ void gpuKernelArrayAXPB(T *y, T *x, T a, T b, int n)
{
    int tid = threadIdx.x + blockIdx.x * blockDim.x;
    while (tid < n) {
        y[tid] = a*x[tid]+b;      
        tid += blockDim.x * gridDim.x;
    }
}

template <typename T> void gpuArrayAXPB(T *y, T *x, T a, T b, int n)
{        
    int blockDim = 256;
    int gridDim = (n + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;
    gpuKernelArrayAXPB<<<gridDim, blockDim>>>(y, x, a, b, n);
}


template <typename T>
__global__ void gpuKernelArrayAXPBY(T *z, T *x, T *y, T a, T b, int n)
{
    int tid = threadIdx.x + blockIdx.x * blockDim.x;
    while (tid < n) {
        z[tid] = a*x[tid]+b*y[tid];      
        tid += blockDim.x * gridDim.x;
    }
}

template <typename T> void gpuArrayAXPBY(T *z, T *x, T *y, T a, T b, int n)
{        
    int blockDim = 256;
    int gridDim = (n + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;
    gpuKernelArrayAXPBY<<<gridDim, blockDim>>>(z, x, y, a, b, n);
}

template <typename T>
__global__ void gpuKernelArrayAXY(T *s, T *x, T *y, T a, int n)
{
    int tid = threadIdx.x + blockIdx.x * blockDim.x;
    while (tid < n) {
        s[tid] = a*x[tid]*y[tid];      
        tid += blockDim.x * gridDim.x;
    }
}

template <typename T> void gpuArrayAXY(T *s, T *x, T *y, T a, int n)
{        
    int blockDim = 256;
    int gridDim = (n + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;
    gpuKernelArrayAXY<<<gridDim, blockDim>>>(s, x, y, a, n);
}

template <typename T>
__global__ void gpuKernelArrayAXYZ(T *s, T *x, T *y, T *z, T a, int n)
{
    int tid = threadIdx.x + blockIdx.x * blockDim.x;
    while (tid < n) {
        s[tid] = a*x[tid]*y[tid]*z[tid];      
        tid += blockDim.x * gridDim.x;
    }
}

template <typename T> void gpuArrayAXYZ(T *s, T *x, T *y, T *z, T a, int n)
{        
    int blockDim = 256;
    int gridDim = (n + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;
    gpuKernelArrayAXYZ<<<gridDim, blockDim>>>(s, x, y, z, a, n);
}

template <typename T>
__global__ void gpuKernelArrayAXYPBZ(T *s, T *x, T *y, T *z, T a, T b, int n)
{
    int tid = threadIdx.x + blockIdx.x * blockDim.x;
    while (tid < n) {
        s[tid] = a*x[tid]*y[tid] + b*z[tid];      
        tid += blockDim.x * gridDim.x;
    }
}

template <typename T> void gpuArrayAXYPBZ(T *s, T *x, T *y, T *z, T a, T b, int n)
{        
    int blockDim = 256;
    int gridDim = (n + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;
    gpuKernelArrayAXYPBZ<<<gridDim, blockDim>>>(s, x, y, z, a, b, n);
}

template <typename T>
__global__ void gpuKernelArrayAdd3Vectors(T *s, T *x, T *y, T *z, T a, T b, T c, int n)
{
    int tid = threadIdx.x + blockIdx.x * blockDim.x;
    while (tid < n) {
        s[tid] = a*x[tid] + b*y[tid] + c*z[tid];      
        tid += blockDim.x * gridDim.x;
    }
}

template <typename T> void gpuArrayAdd3Vectors(T *s, T *x, T *y, T *z, T a, T b, T c, int n)
{        
    int blockDim = 256;
    int gridDim = (n + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;
    gpuKernelArrayAdd3Vectors<<<gridDim, blockDim>>>(s, x, y, z, a, b, c, n);
}

template <typename T>
__global__ void  gpuKernelArrayAdd3Vector(T *a, T *b, T *c, T *d, int n)
{
    int tid = threadIdx.x + blockIdx.x * blockDim.x;
    while (tid < n) {        
        a[tid] = b[tid] + c[tid] + d[tid];
        tid += blockDim.x * gridDim.x;
    }
}

template <typename T> void gpuArrayAdd3Vector(T *a, T *b, T *c, T *d, int n)
{        
    int blockDim = 256;
    int gridDim = (n + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;
    gpuKernelArrayAdd3Vector<<<gridDim, blockDim>>>(a, b, c, d, n);
}

template <typename T>
__global__ void gpuKernelArrayExtract(T *un, T *u, int I, int J, int M, int N,  
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

template <typename T> void gpuArrayExtract(T *un, T *u, int I, int J, int K, 
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
    gpuKernelArrayExtract<<<gridDim, blockDim>>>(un, u, I, J, M, N, i1, j1, k1, ni);
}

template <typename T>
__global__ void gpuKernelArrayInsert(T *u, T *un, int I, int J, int M, int N,  
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

template <typename T> void gpuArrayInsert(T *u, T *un, int I, int J, int K, 
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
    gpuKernelArrayInsert<<<gridDim, blockDim>>>(u, un, I, J, M, N, i1, j1, k1, ni);
}

template <typename T> 
__global__ void gpuKernelArrayGemmSharedMem(T *C, T *A, T *B, int I, int J, int K, int N, int Q)
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

template <typename T> void gpuArrayGemmSharedMem(T *C, T *A, T *B, int I, int J, int K)
{        
    // C[I*J] = A[I*K] x B[K*J]
    int N = I*J;    
    int Q = I*K;
    int blockDim = 256;
    int gridDim = (N + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;
    gpuKernelArrayGemmSharedMem<<<gridDim, blockDim>>>(C, A, B, I, J, K, N, Q);
}

template <typename T> 
__global__ void gpuKernelArrayGemmBatch(T *C, T *A, T *B, int I, int J, int K, int S, int M, int N, int P, int Q)
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

template <typename T> void gpuArrayGemmBatch(T *C, T *A, T *B, int I, int J, int K, int S)
{        
    // C[I*J*S] = A[I*K*S] x B[K*J*S]
    int M = I*J;
    int N = M*S;
    int Q = I*K;
    int P = K*J;
    int blockDim = 256;
    int gridDim = (N + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;
    gpuKernelArrayGemmBatch<<<gridDim, blockDim>>>(C, A, B, I, J, K, S, M, N, P, Q);
}

template <typename T> 
__global__ void gpuKernelArrayGemmBatch1(T *C, T *A, T *B, int I, int J, int K, int S, int M, int N, int P, int Q)
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

template <typename T> void gpuArrayGemmBatch1(T *C, T *A, T *B, int I, int J, int K, int S)
{        
    // C[I*J*S] = A[I*K*S] x B[K*J*S] + C[I*J*S]
    int M = I*J;
    int N = M*S;
    int Q = I*K;
    int P = K*J;
    int blockDim = 256;
    int gridDim = (N + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;
    gpuKernelArrayGemmBatch1<<<gridDim, blockDim>>>(C, A, B, I, J, K, S, M, N, P, Q);
}

template <typename T> 
__global__ void gpuKernelArrayDG2CG(T *ucg, T *udg, int *cgent2dgent, int *rowent2elem, int nent)
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

template <typename T> void gpuArrayDG2CG(T *ucg, T *udg, int *cgent2dgent, int *rowent2elem, int nent)
{        
    int blockDim = 256;
    int gridDim = (nent + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;
    gpuKernelArrayDG2CG<<<gridDim, blockDim>>>(ucg, udg, cgent2dgent, rowent2elem, nent);
}

template <typename T> 
__global__ void gpuKernelArrayDG2CG2(T *ucg, T *udg, int *colent2elem, int *rowent2elem, int nent, int npe)
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
template <typename T> void gpuArrayDG2CG2(T *ucg, T *udg, int *colent2elem, int *rowent2elem, int nent, int npe)
{        
    int blockDim = 256;
    int gridDim = (nent + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;
    gpuKernelArrayDG2CG2<<<gridDim, blockDim>>>(ucg, udg, colent2elem, rowent2elem, nent, npe);
}

template <typename T>
__global__ void gpuKernelArrayABPXYZ(T *fij, T *c3ij, T *eij, T *d3ij, T *g3ij, int dim, int n)
{
    int i = threadIdx.x + blockIdx.x * blockDim.x;
    while (i < n) {        
        for (int j=0; j<dim; j++)
            fij[j+3*i] = fij[j+3*i]*c3ij[i] + eij[i]*d3ij[i]*g3ij[j+3*i];  
        i += blockDim.x * gridDim.x;
    }
}
template <typename T> void gpuArrayABPXYZ(T *fij, T *c3ij, T *eij, T *d3ij, T *g3ij, int dim, int n)
{        
    int blockDim = 256;
    int gridDim = (n + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;
    gpuKernelArrayABPXYZ<<<gridDim, blockDim>>>(fij, c3ij, eij, d3ij, g3ij, dim, n);
}
template void gpuArrayABPXYZ(double*, double*, double*, double*, double*, int, int);
template void gpuArrayABPXYZ(float*, float*, float*, float*, float*, int, int);

// template void gpuPrint2DArray(double*, int, int);
// template void gpuPrint3DArray(double*, int, int, int);
template void gpuArraySetValue(double*, double, int);
template void gpuArraySetValueAtIndex(double*, double, int);
template double gpuArrayGetValueAtIndex(double*, int);
template void gpuArrayAddScalar(double*, double, int);
template void gpuArrayMultiplyScalar(double*, double, int);

template void gpuGetArraySpacing(double*, double*, int, int);
template void gpuGetArrayAtIndex(double*, double*, int*, int);
template void gpuPutArrayAtIndex(double*, double*, int*, int);
template void gpuArrayAXPYAtIndex(double*, double*, double, int*, int);
template void gpuArrayPlusXAtIndex(double*, double*, int*, int);
template void gpuArrayMinusXAtIndex(double*, double*, int*, int);

template void gpuArrayCopy(double*, double*, int);
template void gpuArrayMinus(double*, double*, int);
template void gpuArrayAbs(double*, double*, int);
template void gpuArraySqrt(double*, double*, int);
template void gpuArraySin(double*, double*, int);
template void gpuArrayCos(double*, double*, int);
template void gpuArrayTan(double*, double*, int);
template void gpuArrayAsin(double*, double*, int);
template void gpuArrayAcos(double*, double*, int);
template void gpuArrayAtan(double*, double*, int);
template void gpuArraySinh(double*, double*, int);
template void gpuArrayCosh(double*, double*, int);
template void gpuArrayTanh(double*, double*, int);
template void gpuArrayAsinh(double*, double*, int);
template void gpuArrayAcosh(double*, double*, int);
template void gpuArrayAtanh(double*, double*, int);
template void gpuArrayExp(double*, double*, int);
template void gpuArrayLog(double*, double*, int);
template void gpuArrayCeil(double*, double*, int);
template void gpuArrayFloor(double*, double*, int);
template void gpuArrayErf(double*, double*, int);
template void gpuArrayErfc(double*, double*, int);
template void gpuArraySquare(double*, double*, int);
template void gpuArrayPower(double*, double*, int, int);

template void gpuArrayMultiplyScalarDiagonal(double*, double, int);
template void gpuArrayAddVectorToDiagonal(double*, double*, double, int);
template void gpuArrayRowAverage(double*, double*, int, int);
template void gpuArrayRowSum(double*, double*, int, int);
template void gpuArrayRowSquareSum(double*, double*, int, int);
template void gpuArrayDistSquareSum(double*, double*, double *, int, int);
template void gpuArrayAXPB(double*, double*, double, double, int);
template void gpuArrayAXPBY(double*, double*, double*, double, double, int);
template void gpuArrayAXY(double*, double*, double*, double, int);
template void gpuArrayAXYZ(double*, double*, double*, double*, double, int);
template void gpuArrayAXYPBZ(double*, double*, double*, double*, double, double, int);
template void gpuArrayAdd3Vectors(double*, double*, double*, double*, double, double, double, int);
template void gpuArrayAdd3Vector(double*, double*, double*, double*, int);
template void gpuArrayExtract(double*, double*, int, int, int, int, int, int, int, int, int);
template void gpuArrayInsert(double*, double*, int, int, int, int, int, int, int, int, int);
template void gpuArrayGemmSharedMem(double*, double*, double*, int, int, int);
template void gpuArrayGemmBatch(double*, double*, double*, int, int, int, int);
template void gpuArrayGemmBatch1(double*, double*, double*, int, int, int, int);
template void gpuArrayDG2CG(double*, double*, int*, int*, int);
template void gpuArrayDG2CG2(double*, double*, int*, int*, int, int);

// template void gpuPrint2DArray(float*, int, int);
// template void gpuPrint3DArray(float*, int, int, int);
template void gpuArraySetValue(float*, float, int);
template void gpuArraySetValueAtIndex(float*, float, int);
template float gpuArrayGetValueAtIndex(float*, int);
template void gpuArrayAddScalar(float*, float, int);
template void gpuArrayMultiplyScalar(float*, float, int);

template void gpuGetArraySpacing(float*, float*, int, int);
template void gpuGetArrayAtIndex(float*, float*, int*, int);
template void gpuPutArrayAtIndex(float*, float*, int*, int);
template void gpuArrayAXPYAtIndex(float*, float*, float, int*, int);
template void gpuArrayPlusXAtIndex(float*, float*, int*, int);
template void gpuArrayMinusXAtIndex(float*, float*, int*, int);

template void gpuArrayCopy(float*, float*, int);
template void gpuArrayMinus(float*, float*, int);
template void gpuArrayAbs(float*, float*, int);
template void gpuArraySqrt(float*, float*, int);
template void gpuArraySin(float*, float*, int);
template void gpuArrayCos(float*, float*, int);
template void gpuArrayTan(float*, float*, int);
template void gpuArrayAsin(float*, float*, int);
template void gpuArrayAcos(float*, float*, int);
template void gpuArrayAtan(float*, float*, int);
template void gpuArraySinh(float*, float*, int);
template void gpuArrayCosh(float*, float*, int);
template void gpuArrayTanh(float*, float*, int);
template void gpuArrayAsinh(float*, float*, int);
template void gpuArrayAcosh(float*, float*, int);
template void gpuArrayAtanh(float*, float*, int);
template void gpuArrayExp(float*, float*, int);
template void gpuArrayLog(float*, float*, int);
template void gpuArrayCeil(float*, float*, int);
template void gpuArrayFloor(float*, float*, int);
template void gpuArrayErf(float*, float*, int);
template void gpuArrayErfc(float*, float*, int);
template void gpuArraySquare(float*, float*, int);
template void gpuArrayPower(float*, float*, int, int);

template void gpuArrayMultiplyScalarDiagonal(float*, float, int);
template void gpuArrayAddVectorToDiagonal(float*, float*, float, int);
template void gpuArrayRowAverage(float*, float*, int, int);
template void gpuArrayRowSum(float*, float*, int, int);
template void gpuArrayRowSquareSum(float*, float*, int, int);
template void gpuArrayDistSquareSum(float*, float*, float*, int, int);
template void gpuArrayAXPB(float*, float*, float, float, int);
template void gpuArrayAXPBY(float*, float*, float*, float, float, int);
template void gpuArrayAXY(float*, float*, float*, float, int);
template void gpuArrayAXYZ(float*, float*, float*, float*, float, int);
template void gpuArrayAXYPBZ(float*, float*, float*, float*, float, float, int);
template void gpuArrayAdd3Vectors(float*, float*, float*, float*, float, float, float, int);
template void gpuArrayAdd3Vector(float*, float*, float*, float*, int);
template void gpuArrayExtract(float*, float*, int, int, int, int, int, int, int, int, int);
template void gpuArrayInsert(float*, float*, int, int, int, int, int, int, int, int, int);
template void gpuArrayGemmSharedMem(float*, float*, float*, int, int, int);
template void gpuArrayGemmBatch(float*, float*, float*, int, int, int, int);
template void gpuArrayGemmBatch1(float*, float*, float*, int, int, int, int);
template void gpuArrayDG2CG(float*, float*, int*, int*, int);
template void gpuArrayDG2CG2(float*, float*, int*, int*, int, int);

template int gpuArrayGetValueAtIndex(int*, int);
template void gpuArrayCopy(int*, int*, int);
template void gpuArraySetValue(int*, int, int);

template void gpuArrayTranspose(double *A, double *B, int m, int n);
template void gpuArrayTranspose(float *A, float *B, int m, int n);
template void gpuArrayTranspose(int *A, int *B, int m, int n);

template void gpuArrayPlusAtColumnIndex(double *A, double *B, int *colind, int m, int n);
template void gpuArrayPlusAtColumnIndex(float *A, float *B, int *colind, int m, int n);
template void gpuArrayPlusAtColumnIndex(int *A, int *B, int *colind, int m, int n);

template void gpuArrayMinusAtColumnIndex(double *A, double *B, int *colind, int m, int n);
template void gpuArrayMinusAtColumnIndex(float *A, float *B, int *colind, int m, int n);
template void gpuArrayMinusAtColumnIndex(int *A, int *B, int *colind, int m, int n);

template void gpuGetArrayAtColumnIndex(double *A, double *B, int *colind, int m, int n);
template void gpuGetArrayAtColumnIndex(float *A, float *B, int *colind, int m, int n);
template void gpuGetArrayAtColumnIndex(int *A, int *B, int *colind, int m, int n);

template void gpuArrayTransposeAtColumnIndex(double *A, double *B, int *colind, int m, int n);
template void gpuArrayTransposeAtColumnIndex(float *A, float *B, int *colind, int m, int n);
template void gpuArrayTransposeAtColumnIndex(int *A, int *B, int *colind, int m, int n);

template void gpuGetArrayAtRowIndex(double *A, double *B, int *colind, int m, int n);
template void gpuGetArrayAtRowIndex(float *A, float *B, int *colind, int m, int n);
template void gpuGetArrayAtRowIndex(int *A, int *B, int *colind, int m, int n);

template void gpuArrayTransposeAtRowIndex(double *A, double *B, int *colind, int m, int n);
template void gpuArrayTransposeAtRowIndex(float *A, float *B, int *colind, int m, int n);
template void gpuArrayTransposeAtRowIndex(int *A, int *B, int *colind, int m, int n);


#endif


