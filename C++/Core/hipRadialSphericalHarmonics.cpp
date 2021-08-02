/***************************************************************************                               
                    Molecular Dynamics Potentials (MDP)
                           CESMIX-MIT Project  
 
 Contributing authors: Ngoc-Cuong Nguyen (cuongng@mit.edu, exapde@gmail.com)
 ***************************************************************************/

#include "hip/hip_runtime.h"

#ifndef __HIPRADIALSPHERICALHARMONICS
#define __HIPRADIALSPHERICALHARMONICS

#include "hipcub/device/device_radix_sort.hpp"
#include "hipcub/device/device_run_length_encode.hpp"
#include "hipcub/device/device_scan.hpp"
#include "hipcub/device/device_segmented_reduce.hpp"

template <typename T> __global__ void hipKernelRadialSphericalHarmonics(T *Sr, T *Si, T *Ylmr, T *Ylmi, 
        T *g, int K, int N, int N3)
{                        
    // Ylmr  [N*(L+1)*(L+1)] : real spherical harmonics functions
    // Ylmi  [N*(L+1)*(L+1)] : imag spherical harmonics functions    
    // g  [N*K*(L+1)]:         spherical Bessel functions
    // Sr  [N*K*(L+1)*(L+1)] : real radial spherical harmonics functions
    // Si  [N*K*(L+1)*(L+1)] : imag radial spherical harmonics functions
        
    int tid = threadIdx.x + blockIdx.x * blockDim.x;    
    while (tid < N3) { //N3 = K*(L+1)^2*N,  N2 = K*(L+1)^2
        int i = tid%N;     // [1,N]           
        int q = (tid-i)/N; // [1, K*(L+1)^2] 
        int k = q%K;        // [1, K] 
        int j = (q-k)/K;    // [1, (L+1)^2]                                          
//         int q = tid%N2;     // [1, K*(L+1)^2]           
//         int i = (tid-q)/N2; // [1, N] 
//         int k = q%K;        // [1, K] 
//         int j = (q-k)/K;    // [1, (L+1)^2]                              
        // j = 0 -> l = 0, m = 0
        // j = 1 -> l = 1, m = -1
        // j = 2 -> l = 1, m = 0          
        // j = 3 -> l = 1, m = 1    
        // j = 4 -> l = 2, m = -2
        // j = 5 -> l = 2, m = -1          
        // j = 6 -> l = 2, m = 0                          
        // j = 7 -> l = 2, m = 1          
        // j = 8 -> l = 2, m = 2                          
        // j = (l*l + l + m)    
        // determine l
        // l*l <= j 
        int l=0;
        for (l=0; l<1000; l++)                 
            if (l*l > j)
                break;        
        l = l - 1;
        int m = j - (l*l + l);                
        int Nl = (l*l + l + m)*N + i;
        int Nk = (l*K + k)*N + i;
        int Nj = ((l*l + l + m)*K + k)*N + i;                
        Sr[Nj] = g[Nk]*Ylmr[Nl];
        Si[Nj] = g[Nk]*Ylmi[Nl];

        tid += blockDim.x * gridDim.x;
    }
}
template <typename T> void hipRadialSphericalHarmonics(T *Sr, T *Si, T *Ylmr, T *Ylmi, 
        T *g, int L, int K, int N)
{        
    //N3 = K*(L+1)^2*N,  N2 = K*(L+1)^2
    int N2 = K*(L+1)*(L+1);
    int N3 = N2*N;                        
    int blockDim = 256;
    int gridDim = (N3 + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;
    hipLaunchKernelGGL(hipKernelRadialSphericalHarmonics, gridDim, blockDim, 0, 0, Sr, Si, Ylmr, Ylmi, 
        g, K, N, N3);
}
template void hipRadialSphericalHarmonics(double*, double*, double*, double*, double*, int, int, int);
template void hipRadialSphericalHarmonics(float*, float*, float*, float*, float*, int, int, int);

template <typename T> __global__ void hipKernelRadialSphericalHarmonicsWithDeriv(T *Sr, T *Si, 
        T *Srx, T *Six, T *Sry, T *Siy, T *Srz, T *Siz, T *Ylmr, T *Ylmi, T *Ylmrx, T *Ylmix, T *Ylmry, T *Ylmiy, 
        T *Ylmrz, T *Ylmiz,  T *g, T *gx, T *gy, T *gz, int K, int N, int N3)
{                        
    // Ylmr  [N*(L+1)*(L+1)] : real spherical harmonics functions
    // Ylmi  [N*(L+1)*(L+1)] : imag spherical harmonics functions    
    // g  [N*K*(L+1)]:         spherical Bessel functions
    // Sr  [N*K*(L+1)*(L+1)] : real radial spherical harmonics functions
    // Si  [N*K*(L+1)*(L+1)] : imag radial spherical harmonics functions
    
    int tid = threadIdx.x + blockIdx.x * blockDim.x;    
    while (tid < N3) { //N3 = K*(L+1)^2*N,  N2 = K*(L+1)^2
        int i = tid%N;     // [1,N]           
        int q = (tid-i)/N; // [1, K*(L+1)^2] 
        int k = q%K;        // [1, K] 
        int j = (q-k)/K;    // [1, (L+1)^2]                                          
        int l=0;
        for (l=0; l<1000; l++)                 
            if (l*l > j)
                break;        
        l = l - 1;
        int m = j - (l*l + l);                
        int Nl = (l*l + l + m)*N + i;
        int Nk = (l*K + k)*N + i;
        int Nj = ((l*l + l + m)*K + k)*N + i;       
        Sr[Nj] = g[Nk]*Ylmr[Nl];
        Si[Nj] = g[Nk]*Ylmi[Nl];    
        Srx[Nj] = gx[Nk]*Ylmr[Nl] + g[Nk]*Ylmrx[Nl];
        Sry[Nj] = gy[Nk]*Ylmr[Nl] + g[Nk]*Ylmry[Nl];
        Srz[Nj] = gz[Nk]*Ylmr[Nl] + g[Nk]*Ylmrz[Nl];
        Six[Nj] = gx[Nk]*Ylmi[Nl] + g[Nk]*Ylmix[Nl];
        Siy[Nj] = gy[Nk]*Ylmi[Nl] + g[Nk]*Ylmiy[Nl];
        Siz[Nj] = gz[Nk]*Ylmi[Nl] + g[Nk]*Ylmiz[Nl];                                   
        tid += blockDim.x * gridDim.x;
    }
}
template <typename T> void hipRadialSphericalHarmonicsWithDeriv(T *Sr, T *Si, 
        T *Srx, T *Six, T *Sry, T *Siy, T *Srz, T *Siz, T *Ylmr, T *Ylmi, T *Ylmrx, T *Ylmix, T *Ylmry, T *Ylmiy, 
        T *Ylmrz, T *Ylmiz,  T *g, T *gx, T *gy, T *gz, int L, int K, int N)
{        
    //N3 = K*(L+1)^2*N,  N2 = K*(L+1)^2
    int N2 = K*(L+1)*(L+1);
    int N3 = N2*N;                        
    int blockDim = 256;
    int gridDim = (N3 + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;
    hipLaunchKernelGGL(hipKernelRadialSphericalHarmonicsWithDeriv, gridDim, blockDim, 0, 0, Sr, Si, 
            Srx, Six, Sry, Siy, Srz, Siz, Ylmr, Ylmi, Ylmrx, Ylmix, Ylmry, Ylmiy, Ylmrz, Ylmiz, 
            g, gx, gy, gz, K, N, N3);
}
template void hipRadialSphericalHarmonicsWithDeriv(double*, double*, double*, double*, double*, double*, double*, double*, 
        double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, int, int, int);
template void hipRadialSphericalHarmonicsWithDeriv(float*, float*, float*, float*, float*, float*, float*, float*, 
        float*, float*, float*, float*, float*, float*, float*, float*, float*, float*, float*, float*, int, int, int);

template <typename T> __global__ void hipKernelRadialSphericalHarmonicsNeighborSum(T *ar, T *ai, T *Sr, T *Si, 
        int *Nnb, int Na, int L, int K, int N3)
{                            
    // L                         : the maximum degree of spherical harmonics
    // K                         : the maximum degree of spherical Bessel functions
    // Na                        : number of atoms in the simulation domain 
    // Nnb [Na+1]                : a list containing the number of neighbors for each global atom    
    // N = Nnb[Na]-Nnb[0]        : total number of neighbors
    // Sr [N*K*(L+1)*(L+1)]      : real spherical harmonics Bessel functions
    // Si [N*K*(L+1)*(L+1)]      : imag spherical harmonics Bessel functions        
    // ar [Na*K*(L+1)*(L+1)]     : sum of real spherical harmonics Bessel functions
    // ai [Na*K*(L+1)*(L+1)]     : sum of imag spherical harmonics Bessel functions    
        
    // Sum the spherical harmonics Bessel functions over neighbors            
    int tid = threadIdx.x + blockIdx.x * blockDim.x;    
    int N = Nnb[Na]-Nnb[0]; // total number of neighbors    
    while (tid < N3) { //N3 = K*(L+1)^2*N,  N2 = K*(L+1)^2
        int n = tid%Na;     // [1,Na]           
        int q = (tid-n)/Na; // [1, K*(L+1)^2] 
        int k = q%K;        // [1, K] 
        int l = (q-k)/K;    // [1, (L+1)^2]                                          
        int j = (l*K + k)*Na + n; // index of atom n
        int m = (l*K + k)*N + (Nnb[n]-Nnb[0]); //  starting index of neighbors
        ar[j] = 0.0; 
        ai[j] = 0.0;
        int numnb = Nnb[n+1]-Nnb[n];  // number of neighbors for atom n
        for (int i=0; i<numnb; i++) { // loop over neighbors 
            ar[j] += Sr[m + i]; // sum over neighbors of atom n
            ai[j] += Si[m + i];                                                    
        }
        tid += blockDim.x * gridDim.x;
    }                
}
template <typename T> void hipRadialSphericalHarmonicsNeighborSum(T *ar, T *ai, T *Sr, T *Si, 
        int *Nnb, int Na, int L, int K)
{
    int N2 = K*(L+1)*(L+1);
    int N3 = N2*Na;                        
    int blockDim = 256;
    int gridDim = (N3 + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;
    hipLaunchKernelGGL(hipKernelRadialSphericalHarmonicsNeighborSum, gridDim, blockDim, 0, 0, ar, ai, Sr, Si, Nnb, Na, L, K, N3);            
}
template void hipRadialSphericalHarmonicsNeighborSum(double*, double*, double*, double*, int*, int, int, int);
template void hipRadialSphericalHarmonicsNeighborSum(float*, float*, float*, float*, int *, int, int, int);

template <typename T> __global__ void hipKernelRadialSphericalHarmonicsPowerSpectrum(T *p, T *ar, T *ai, int *indk, 
        int Na, int L, int K, int K2, int N3)
{   
    // L                      : the maximum degree of spherical harmonics
    // K                      : the maximum degree of spherical Bessel functions
    // Na                     : number of atoms in the simulation domain 
    // indk [K*(K+1)/2]       : store indices for symmetric tensor products of spherical Bessel functions
    // ar [Na*K*(L+1)*(L+1)]  : sum of real spherical harmonics Bessel functions
    // ai [Na*K*(L+1)*(L+1)]  : sum of imag spherical harmonics Bessel functions    
    // p [Na*(L+1)*K*(K+1)/2] : power spectrum components
    
    // Compute the power spectrum components for radial spherical harmonics    
    int tid = threadIdx.x + blockIdx.x * blockDim.x;    
    while (tid < N3) { //N3 = K*(L+1)^2*Na,  N2 = K*(L+1)^2
        int n = tid%Na;     // [1,Na]           
        int q = (tid-n)/Na; // [1, K2*(L+1)] 
        int l = q%(L+1);    // [1, (L+1)] 
        int k = (q-l)/(L+1); // [1, K*(K+1)/2]                                          
        int k2 = indk[k];
        int k1 = indk[K2+k];    
        int l1 = (k*(L+1) + l)*Na + n; // global index of atom n
        p[l1] = (T) 0.0;
        for (int m=0; m<(2*l+1); m++) { // loop over magnetic quantum number
            int i1 = ((l*l + m)*K + k1)*Na + n;
            int j1 = ((l*l + m)*K + k2)*Na + n;
            p[l1] += ar[i1]*ar[j1] + ai[i1]*ai[j1]; // power spectrum           
        }
        tid += blockDim.x * gridDim.x;
    }
}
template <typename T> void hipRadialSphericalHarmonicsPowerSpectrum(T *p, T *ar, T *ai, int *indk, int Na, int L, int K)
{
    int K2 = K*(K+1)/2;
    int N2 = K2*(L+1);
    int N3 = N2*Na;                        
    int blockDim = 256;
    int gridDim = (N3 + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;
    hipLaunchKernelGGL(hipKernelRadialSphericalHarmonicsPowerSpectrum, gridDim, blockDim, 0, 0, p, ar, ai, indk, Na, L, K, K2, N3);            
}
template void hipRadialSphericalHarmonicsPowerSpectrum(double*, double*, double*, int*, int, int, int);
template void hipRadialSphericalHarmonicsPowerSpectrum(float*, float*, float*, int*, int, int, int);

template <typename T> __global__ void hipKernelRadialSphericalHarmonicsPowerSpectrumDeriv(T *px, T *py, T *pz, T *ar, T *ai, 
        T *arx, T *aix, T *ary, T *aiy, T *arz, T *aiz, int *indk, int *Nnb, int Na, int L, int K, int K2, int N3)
{   
    // L                      : the maximum degree of spherical harmonics
    // K                      : the maximum degree of spherical Bessel functions    
    // N                      : total number of neighbors 
    // Na                     : number of global atoms in the simulation domain 
    // Nnb [Na+1]             : a list containing the number of neighbors for each global atom    
    // indk [K*(K+1)/2]       : store indices for symmetric tensor products of spherical Bessel functions
    // ar   [Na*K*(L+1)*(L+1)]: sum of real spherical harmonics Bessel functions
    // ai   [Na*K*(L+1)*(L+1)]: sum of imag spherical harmonics Bessel functions    
    // arx  [N*K*(L+1)*(L+1)] : x-derivative of real spherical harmonics Bessel functions
    // aix  [N*K*(L+1)*(L+1)] : x-derivative of imag spherical harmonics Bessel functions
    // ary  [N*K*(L+1)*(L+1)] : y-derivative of real spherical harmonics Bessel functions
    // aiy  [N*K*(L+1)*(L+1)] : y-derivative of imag spherical harmonics Bessel functions
    // arz  [N*K*(L+1)*(L+1)] : z-derivative of real spherical harmonics Bessel functions
    // aiz  [N*K*(L+1)*(L+1)] : z-derivative of imag spherical harmonics Bessel functions
    // p [(L+1)*K*(K+1)/2]    : power spectrum components
    // px [N*(L+1)*K*(K+1)/2] : x-derivative of power spectrum components
    // py [N*(L+1)*K*(K+1)/2] : y-derivative of power spectrum components
    // pz [N*(L+1)*K*(K+1)/2] : z-derivative of power spectrum components
    
    // Compute partial derivatives of the power spectrum components    
    int tid = threadIdx.x + blockIdx.x * blockDim.x;    
    int N = Nnb[Na]-Nnb[0]; // total number of neighbors    
    while (tid < N3) { //N3 = K*(L+1)^2*Na,  N2 = K*(L+1)^2
        int n = tid%Na;     // [1,Na]           
        int q = (tid-n)/Na; // [1, K2*(L+1)] 
        int l = q%(L+1);    // [1, (L+1)] 
        int k = (q-l)/(L+1); // [1, K*(K+1)/2]                           
        int k2 = indk[k];
        int k1 = indk[K2+k];    
        int numnb = Nnb[n+1]-Nnb[n];  // number of neighbors for atom n
        for (int i=0; i<numnb; i++) {// loop over each neighbor of atom n                         
            int j = (k*(L+1) + l)*N + Nnb[n] + i; // index of px, py, pz        
            px[j] = (T) 0.0;
            py[j] = (T) 0.0;
            pz[j] = (T) 0.0;
            for (int m=0; m<(2*l+1); m++) {
                int i1 = ((l*l + m)*K + k1)*Na + n; // index of ar and ai 
                int j1 = ((l*l + m)*K + k2)*Na + n; // index of ar and ai                                         
                int i2 = ((l*l + m)*K + k1)*N + Nnb[n] + i;  // index of arx and aix     
                int j2 = ((l*l + m)*K + k2)*N + Nnb[n] + i;  // index of arx and aix                                                    
                px[j] += ar[i1]*arx[j2] + arx[i2]*ar[j1] + ai[i1]*aix[j2] + aix[i2]*ai[j1];
                py[j] += ar[i1]*ary[j2] + ary[i2]*ar[j1] + ai[i1]*aiy[j2] + aiy[i2]*ai[j1];
                pz[j] += ar[i1]*arz[j2] + arz[i2]*ar[j1] + ai[i1]*aiz[j2] + aiz[i2]*ai[j1];                    
            }                    
        }
        tid += blockDim.x * gridDim.x;
    }
}
template <typename T> void hipRadialSphericalHarmonicsPowerSpectrumDeriv(T *px, T *py, T *pz, T *ar, T *ai, 
        T *arx, T *aix, T *ary, T *aiy, T *arz, T *aiz, int *indk, int *Nnb, int Na, int L, int K)
{    
    int K2 = K*(K+1)/2;
    int N2 = K2*(L+1);
    int N3 = N2*Na;                        
    int blockDim = 256;
    int gridDim = (N3 + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;
    hipLaunchKernelGGL(hipKernelRadialSphericalHarmonicsPowerSpectrumDeriv, gridDim, blockDim, 0, 0, px, py, pz, ar, ai, 
            arx, aix, ary, aiy, arz, aiz, indk, Nnb, Na, L, K, K2, N3);            
}
template void hipRadialSphericalHarmonicsPowerSpectrumDeriv(double*, double*, double*, double*, double*, double*, 
        double*, double*, double*, double*, double*, int*, int*, int, int, int);
template void hipRadialSphericalHarmonicsPowerSpectrumDeriv(float*, float*, float*, float*, float*, float*,  
        float*, float*, float*, float*, float*, int*, int*, int, int, int);

template <typename T> __global__ void hipKernelRadialSphericalHarmonicsPowerSpectrumDeriv2(T *pd, T *ar, T *ai, 
        T *arx, T *aix, T *ary, T *aiy, T *arz, T *aiz, int *indk, int *Nnb, int Na, int L, int K, int K2, int N3)
{   
    // L                      : the maximum degree of spherical harmonics
    // K                      : the maximum degree of spherical Bessel functions    
    // N                      : total number of neighbors 
    // Na                     : number of global atoms in the simulation domain 
    // Nnb [Na+1]             : a list containing the number of neighbors for each global atom    
    // indk [K*(K+1)/2]       : store indices for symmetric tensor products of spherical Bessel functions
    // ar   [Na*K*(L+1)*(L+1)]: sum of real spherical harmonics Bessel functions
    // ai   [Na*K*(L+1)*(L+1)]: sum of imag spherical harmonics Bessel functions    
    // arx  [N*K*(L+1)*(L+1)] : x-derivative of real spherical harmonics Bessel functions
    // aix  [N*K*(L+1)*(L+1)] : x-derivative of imag spherical harmonics Bessel functions
    // ary  [N*K*(L+1)*(L+1)] : y-derivative of real spherical harmonics Bessel functions
    // aiy  [N*K*(L+1)*(L+1)] : y-derivative of imag spherical harmonics Bessel functions
    // arz  [N*K*(L+1)*(L+1)] : z-derivative of real spherical harmonics Bessel functions
    // aiz  [N*K*(L+1)*(L+1)] : z-derivative of imag spherical harmonics Bessel functions
    // p [(L+1)*K*(K+1)/2]    : power spectrum components
    // pd [3*N*(L+1)*K*(K+1)/2] : derivatives of power spectrum components
    
    // Compute partial derivatives of the power spectrum components
    int tid = threadIdx.x + blockIdx.x * blockDim.x;   
    int N = Nnb[Na]-Nnb[0]; // total number of neighbors    
    while (tid < N3) { //N3 = K*(L+1)^2*Na,  N2 = K*(L+1)^2
        int n = tid%Na;     // [1,Na]           
        int q = (tid-n)/Na; // [1, K2*(L+1)] 
        int l = q%(L+1);    // [1, (L+1)] 
        int k = (q-l)/(L+1); // [1, K*(K+1)/2]                           
        int k2 = indk[k];
        int k1 = indk[K2+k];    
        int numnb = Nnb[n+1]-Nnb[n];  // number of neighbors for atom n
        for (int i=0; i<numnb; i++) {// loop over each neighbor of atom n        
            int j = (k*(L+1) + l)*N + Nnb[n] + i; // index of px, py, pz
            pd[3*j+0] = (T) 0.0;
            pd[3*j+1] = (T) 0.0;
            pd[3*j+2] = (T) 0.0;
            for (int m=0; m<(2*l+1); m++) {
                int i1 = ((l*l + m)*K + k1)*Na + n; // index of ar and ai 
                int j1 = ((l*l + m)*K + k2)*Na + n; // index of ar and ai                                         
                int i2 = ((l*l + m)*K + k1)*N + Nnb[n] + i;  // index of arx and aix     
                int j2 = ((l*l + m)*K + k2)*N + Nnb[n] + i;  // index of arx and aix                                                    
                pd[3*j+0] += ar[i1]*arx[j2] + arx[i2]*ar[j1] + ai[i1]*aix[j2] + aix[i2]*ai[j1];
                pd[3*j+1] += ar[i1]*ary[j2] + ary[i2]*ar[j1] + ai[i1]*aiy[j2] + aiy[i2]*ai[j1];
                pd[3*j+2] += ar[i1]*arz[j2] + arz[i2]*ar[j1] + ai[i1]*aiz[j2] + aiz[i2]*ai[j1];                    
            }                
        }
        tid += blockDim.x * gridDim.x;
    }
}
template <typename T> void hipRadialSphericalHarmonicsPowerSpectrumDeriv2(T *pd, T *ar, T *ai, 
        T *arx, T *aix, T *ary, T *aiy, T *arz, T *aiz, int *indk, int *Nnb, int Na, int L, int K)
{
    int K2 = K*(K+1)/2;
    int N2 = K2*(L+1);
    int N3 = N2*Na;                        
    int blockDim = 256;
    int gridDim = (N3 + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;
    hipLaunchKernelGGL(hipKernelRadialSphericalHarmonicsPowerSpectrumDeriv2, gridDim, blockDim, 0, 0, pd, ar, ai, 
            arx, aix, ary, aiy, arz, aiz, indk, Nnb, Na, L, K, K2, N3);            
}
template void hipRadialSphericalHarmonicsPowerSpectrumDeriv2(double*, double*, double*, double*, 
        double*, double*, double*, double*, double*, int*, int*, int, int, int);
template void hipRadialSphericalHarmonicsPowerSpectrumDeriv2(float*, float*, float*, float*,  
        float*, float*, float*, float*, float*, int*, int*, int, int, int);

template <typename T> __global__ void hipKernelRadialSphericalHarmonicsBispectrum(T *b, T *ar, T *ai, T *cg, int *indk, 
        int *indl, int *indm, int *rowm, int Nub, int Ncg, int Na, int L, int K, int K2, int N3)
{   
    // L                      : the maximum degree of spherical harmonics
    // K                      : the maximum degree of spherical Bessel functions
    // Na                     : number of atoms in the simulation domain    
    // Ncg                    : the total number of non-zero Clebsch-Gordan coefficients 
    // Nub                    : the number of non-zero unique spherical harmonics bispectrum components
    // cg   [Ncg*1]           : store non-zero Clebsch-Gordan coefficients
    // indl [Nub*3]           : store the indices of tuple (l2,l1,l)        
    // indk [K*(K+1)/2]       : store indices for symmetric tensor products of spherical Bessel functions
    // indm [Ncg*3]           : store the indices of tuple (m2, m1, m) for non-zero Clebsch-Gordan coefficients
    // rowm [Nub+1]           : store the number of indices of tuple (m2, m1, m) for each tuple (l2, l1, l)    
    // ar  [Na*K*(L+1)*(L+1)] : sum of real spherical harmonics Bessel functions
    // ai  [Na*K*(L+1)*(L+1)] : sum of imag spherical harmonics Bessel functions
    // b [Na*Nub*K*(K+1)/2]   : bispectrum components    
    
    // Compute the bispectrum components for radial spherical harmonics
    
    int tid = threadIdx.x + blockIdx.x * blockDim.x;    
    while (tid < N3) { //N3 = K*(L+1)^2*Na,  N2 = K*(L+1)^2
        int n = tid%Na;     // [1,Na]           
        int q = (tid-n)/Na; // [1, K2*Nub] 
        int i = q%Nub;      // [1, Nub] 
        int k = (q-i)/Nub;  // [1, K*(K+1)/2]                                          
        int k2 = indk[k];
        int k1 = indk[K2+k];    
        int l2 = indl[i];
        int l1 = indl[Nub+i];
        int l = indl[2*Nub+i];   
        int nm = rowm[i+1]-rowm[i];
        int indn = (k*Nub + i)*Na + n; // global index of atom n
        b[indn] = (T) 0.0;
        for (int j = 0; j<nm; j++) {
            int m2 = indm[rowm[i]+j];
            int m1 = indm[Ncg+rowm[i]+j];
            int m = indm[2*Ncg + rowm[i]+j];                  
            T a1, b1, a2, b2, a3, b3;                                
            a1 = ar[((l*l + l + m)*K + k1)*Na + n];
            b1 = ai[((l*l + l + m)*K + k1)*Na + n];
            a2 = ar[((l1*l1 + l1 + m1)*K + k2)*Na + n];
            b2 = ai[((l1*l1 + l1 + m1)*K + k2)*Na + n];
            a3 = ar[((l2*l2 + l2 + m2)*K + k2)*Na + n];
            b3 = ai[((l2*l2 + l2 + m2)*K + k2)*Na + n];
            b[indn] += cg[rowm[i]+j]*(a1*a2*a3 + a2*b1*b3 + a3*b1*b2 - a1*b2*b3);                                              
        }                                
        tid += blockDim.x * gridDim.x;        
    }            
}
template <typename T> void hipRadialSphericalHarmonicsBispectrum(T *b, T *ar, T *ai, T *cg, int *indk, 
        int *indl, int *indm, int *rowm, int Nub, int Ncg, int Na, int L, int K)
{
    int K2 = K*(K+1)/2;
    int N2 = K2*Nub;
    int N3 = N2*Na;                        
    int blockDim = 256;
    int gridDim = (N3 + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;
    hipLaunchKernelGGL(hipKernelRadialSphericalHarmonicsBispectrum, gridDim, blockDim, 0, 0, b, ar, ai, cg, indk, indl, indm, rowm, 
            Nub, Ncg, Na, L, K, K2, N3);            
}
template void hipRadialSphericalHarmonicsBispectrum(double*, double*, double*, double*, int*, int*, int*, int*, int, int, int, int, int);
template void hipRadialSphericalHarmonicsBispectrum(float*, float*, float*, float*, int*, int*, int*, int*, int, int, int, int, int);

template <typename T> __global__ void hipKernelRadialSphericalHarmonicsBispectrumDeriv(T *bx, T *by, T *bz, T *ar, T *ai, 
        T *arx, T *aix, T *ary, T *aiy, T *arz, T *aiz, T*cg, int *indk, int *indl,
        int *indm, int *rowm, int *Nnb, int Nub, int Ncg, int Na, int K, int K2, int N3)
{   
    // K                      : the maximum degree of spherical Bessel functions
    // N                      : total number of neighbors 
    // Na                     : number of global atoms in the simulation domain     
    // Ncg                    : the total number of non-zero Clebsch-Gordan coefficients 
    // Nub                    : the number of non-zero unique spherical harmonics bispectrum components
    // Nnb [Na+1]             : a list containing the number of neighbors for each global atom        
    // indl [Nub*3]           : store the indices of tuple (l2,l1,l)        
    // indk [K*(K+1)/2]       : store indices for symmetric tensor products of spherical Bessel functions
    // indm [Ncg*3]           : store the indices of tuple (m2, m1, m) for non-zero Clebsch-Gordan coefficients
    // rowm [Nub+1]           : store the number of indices of tuple (m2, m1, m) for each tuple (l2, l1, l)    
    // cg   [Ncg*1]           : store non-zero Clebsch-Gordan coefficients
    // ar   [Na*K*(L+1)*(L+1)]: real part of spherical harmonics Bessel functions
    // ai   [Na*K*(L+1)*(L+1)]: imag part of spherical harmonics Bessel functions    
    // arx  [N*K*(L+1)*(L+1)] : x-derivative of real spherical harmonics Bessel functions
    // aix  [N*K*(L+1)*(L+1)] : x-derivative of imag spherical harmonics Bessel functions
    // ary  [N*K*(L+1)*(L+1)] : y-derivative of real spherical harmonics Bessel functions
    // aiy  [N*K*(L+1)*(L+1)] : y-derivative of imag spherical harmonics Bessel functions
    // arz  [N*K*(L+1)*(L+1)] : z-derivative of real spherical harmonics Bessel functions
    // aiz  [N*K*(L+1)*(L+1)] : z-derivative of imag spherical harmonics Bessel functions
    // b    [Na*Nub*K*(K+1)/2]: bispectrum components        
    // bx   [N*Nub*K*(K+1)/2] : x-derivative of power spectrum components
    // by   [N*Nub*K*(K+1)/2] : y-derivative of power spectrum components
    // bz   [N*Nub*K*(K+1)/2] : z-derivative of power spectrum components
        
    // Compute partial derivatives of the bispectrum components
    int tid = threadIdx.x + blockIdx.x * blockDim.x;    
    int N = Nnb[Na]-Nnb[0]; // total number of neighbors    
    while (tid < N3) { //N3 = K*(L+1)^2*Na,  N2 = K*(L+1)^2
        int n = tid%Na;     // [1,Na]           
        int q = (tid-n)/Na; // [1, K2*Nub] 
        int i = q%Nub;      // [1, Nub] 
        int k = (q-i)/Nub;  // [1, K*(K+1)/2]                                          
        int k2 = indk[k];
        int k1 = indk[K2+k];    
        int l2 = indl[i];
        int l1 = indl[Nub+i];
        int l = indl[2*Nub+i];   
        int nm = rowm[i+1]-rowm[i];
        int numnb = Nnb[n+1]-Nnb[n];  // number of neighbors for atom n
        for (int m=0; m<numnb; m++) {// loop over each neighbor of atom n        
            int ii = (k*Nub + i)*N + Nnb[n] + m; // index of bx, by, bz                                        
            bx[ii] = (T) 0.0;
            by[ii] = (T) 0.0;
            bz[ii] = (T) 0.0;                                    
            for (int j = 0; j<nm; j++) {
                int m2 = indm[rowm[i]+j];
                int m1 = indm[Ncg+rowm[i]+j];
                int m = indm[2*Ncg + rowm[i]+j];                            
                int n1 = (l*l + l + m)*K + k1;    
                int n2 = (l1*l1 + l1 + m1)*K + k2;
                int n3 = (l2*l2 + l2 + m2)*K + k2;                        
                int lm1 = n1*Na + n; // index of ar and ai 
                int lm2 = n2*Na + n; // index of ar and ai                                         
                int lm3 = n3*Na + n; // index of ar and ai                                                                 
                int mlk1 = n1*N + Nnb[n] + q; // index of arx and aix      
                int mlk2 = n2*N + Nnb[n] + q; // index of arx and aix            
                int mlk3 = n3*N + Nnb[n] + q; // index of arx and aix                                                  

                T c = cg[rowm[i]+j];                    
                T a1, b1, a2, b2, a3, b3;                                
                T a1x, b1x, a2x, b2x, a3x, b3x;
                T a1y, b1y, a2y, b2y, a3y, b3y;
                T a1z, b1z, a2z, b2z, a3z, b3z;
                a1 = ar[lm1];
                b1 = ai[lm1];
                a2 = ar[lm2];
                b2 = ai[lm2];
                a3 = ar[lm3];
                b3 = ai[lm3];
                a1x = arx[mlk1];
                a1y = ary[mlk1];
                a1z = arz[mlk1];
                b1x = aix[mlk1];
                b1y = aiy[mlk1];
                b1z = aiz[mlk1];
                a2x = arx[mlk2];
                a2y = ary[mlk2];
                a2z = arz[mlk2];
                b2x = aix[mlk2];
                b2y = aiy[mlk2];
                b2z = aiz[mlk2];
                a3x = arx[mlk3];
                a3y = ary[mlk3];
                a3z = arz[mlk3];
                b3x = aix[mlk3];
                b3y = aiy[mlk3];
                b3z = aiz[mlk3];

                T t1 = a1x*a2*a3 + a1*a2x*a3 + a1*a2*a3x;                
                T t2 = a2x*b1*b3 + a2*b1x*b3 + a2*b1*b3x;
                T t3 = a3x*b1*b2 + a3*b1x*b2 + a3*b1*b2x;
                T t4 = a1x*b2*b3 + a1*b2x*b3 + a1*b2*b3x;
                bx[ii] += c*(t1 + t2 + t3 - t4);

                t1 = a1y*a2*a3 + a1*a2y*a3 + a1*a2*a3y;                
                t2 = a2y*b1*b3 + a2*b1y*b3 + a2*b1*b3y;
                t3 = a3y*b1*b2 + a3*b1y*b2 + a3*b1*b2y;
                t4 = a1y*b2*b3 + a1*b2y*b3 + a1*b2*b3y;
                by[ii] += c*(t1 + t2 + t3 - t4);

                t1 = a1z*a2*a3 + a1*a2z*a3 + a1*a2*a3z;                
                t2 = a2z*b1*b3 + a2*b1z*b3 + a2*b1*b3z;
                t3 = a3z*b1*b2 + a3*b1z*b2 + a3*b1*b2z;
                t4 = a1z*b2*b3 + a1*b2z*b3 + a1*b2*b3z;
                bz[ii] += c*(t1 + t2 + t3 - t4);                    
            }                    
        }
        tid += blockDim.x * gridDim.x;   
    }
}
template <typename T> void hipRadialSphericalHarmonicsBispectrumDeriv(T *bx, T *by, T *bz, T *ar, T *ai, 
        T *arx, T *aix, T *ary, T *aiy, T *arz, T *aiz, T*cg, int *indk, int *indl,
        int *indm, int *rowm, int *Nnb, int Nub, int Ncg, int Na, int K)
{    
    int K2 = K*(K+1)/2;
    int N2 = K2*Nub;
    int N3 = N2*Na;                        
    int blockDim = 256;
    int gridDim = (N3 + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;
    hipLaunchKernelGGL(hipKernelRadialSphericalHarmonicsBispectrumDeriv, gridDim, blockDim, 0, 0, bx, by, bz, ar, ai, arx, aix, 
           ary, aiy, arz, aiz, cg, indk, indl, indm, rowm, Nnb, Nub, Ncg, Na, K, K2, N3);            
}
template void hipRadialSphericalHarmonicsBispectrumDeriv(double*, double*, double*, double*, double*, double*, 
        double*, double*, double*, double*, double*, double*, int*, int*, int*, int*, int*, int, int, int, int);
template void hipRadialSphericalHarmonicsBispectrumDeriv(float*, float*, float*, float*, float*, float*,  
        float*, float*, float*, float*, float*, float*, int*, int*, int*, int*, int*, int, int, int, int);

template <typename T> __global__ void hipKernelRadialSphericalHarmonicsBispectrumDeriv2(T *bd, T *ar, T *ai, 
        T *arx, T *aix, T *ary, T *aiy, T *arz, T *aiz, T*cg, int *indk, int *indl,
        int *indm, int *rowm, int *Nnb, int Nub, int Ncg, int Na, int K, int K2, int N3, int npower)
{   
    // K                      : the maximum degree of spherical Bessel functions
    // N                      : total number of neighbors 
    // Na                     : number of global atoms in the simulation domain     
    // Ncg                    : the total number of non-zero Clebsch-Gordan coefficients 
    // Nub                    : the number of non-zero unique spherical harmonics bispectrum components
    // Nnb [Na+1]             : a list containing the number of neighbors for each global atom        
    // indl [Nub*3]           : store the indices of tuple (l2,l1,l)        
    // indk [K*(K+1)/2]       : store indices for symmetric tensor products of spherical Bessel functions
    // indm [Ncg*3]           : store the indices of tuple (m2, m1, m) for non-zero Clebsch-Gordan coefficients
    // rowm [Nub+1]           : store the number of indices of tuple (m2, m1, m) for each tuple (l2, l1, l)    
    // cg   [Ncg*1]           : store non-zero Clebsch-Gordan coefficients
    // ar   [Na*K*(L+1)*(L+1)]: real part of spherical harmonics Bessel functions
    // ai   [Na*K*(L+1)*(L+1)]: imag part of spherical harmonics Bessel functions    
    // arx  [N*K*(L+1)*(L+1)] : x-derivative of real spherical harmonics Bessel functions
    // aix  [N*K*(L+1)*(L+1)] : x-derivative of imag spherical harmonics Bessel functions
    // ary  [N*K*(L+1)*(L+1)] : y-derivative of real spherical harmonics Bessel functions
    // aiy  [N*K*(L+1)*(L+1)] : y-derivative of imag spherical harmonics Bessel functions
    // arz  [N*K*(L+1)*(L+1)] : z-derivative of real spherical harmonics Bessel functions
    // aiz  [N*K*(L+1)*(L+1)] : z-derivative of imag spherical harmonics Bessel functions
    // b    [Na*Nub*K*(K+1)/2]: bispectrum components        
    // bd   [3*N*Nub*K*(K+1)/2] : derivatives of power spectrum components
        
    // Compute partial derivatives of the bispectrum components
    
    int tid = threadIdx.x + blockIdx.x * blockDim.x;    
    int N = Nnb[Na]-Nnb[0]; // total number of neighbors    
    while (tid < N3) { //N3 = K*(L+1)^2*Na,  N2 = K*(L+1)^2
        int n = tid%Na;     // [1,Na]           
        int q = (tid-n)/Na; // [1, K2*Nub] 
        int i = q%Nub;      // [1, Nub] 
        int k = (q-i)/Nub;  // [1, K*(K+1)/2]                                          
        int k2 = indk[k];
        int k1 = indk[K2+k];    
        int l2 = indl[i];
        int l1 = indl[Nub+i];
        int l = indl[2*Nub+i];   
        int nm = rowm[i+1]-rowm[i];
        int numnb = Nnb[n+1]-Nnb[n];  // number of neighbors for atom n
        for (int m=0; m<numnb; m++) {// loop over each neighbor of atom n        
            int ii = (npower + k*Nub + i)*N + Nnb[n] + m; // index of bx, by, bz     
            bd[3*ii+0] = (T) 0.0;
            bd[3*ii+1] = (T) 0.0;
            bd[3*ii+2] = (T) 0.0;                                    
            for (int j = 0; j<nm; j++) {
                int m2 = indm[rowm[i]+j];
                int m1 = indm[Ncg+rowm[i]+j];
                int m = indm[2*Ncg + rowm[i]+j];                            
                int n1 = (l*l + l + m)*K + k1;    
                int n2 = (l1*l1 + l1 + m1)*K + k2;
                int n3 = (l2*l2 + l2 + m2)*K + k2;                        
                int lm1 = n1*Na + n; // index of ar and ai 
                int lm2 = n2*Na + n; // index of ar and ai                                         
                int lm3 = n3*Na + n; // index of ar and ai                                                                 
                int mlk1 = n1*N + Nnb[n] + q; // index of arx and aix      
                int mlk2 = n2*N + Nnb[n] + q; // index of arx and aix            
                int mlk3 = n3*N + Nnb[n] + q; // index of arx and aix                                                  

                T c = cg[rowm[i]+j];                    
                T a1, b1, a2, b2, a3, b3;                                
                T a1x, b1x, a2x, b2x, a3x, b3x;
                T a1y, b1y, a2y, b2y, a3y, b3y;
                T a1z, b1z, a2z, b2z, a3z, b3z;
                a1 = ar[lm1];
                b1 = ai[lm1];
                a2 = ar[lm2];
                b2 = ai[lm2];
                a3 = ar[lm3];
                b3 = ai[lm3];
                a1x = arx[mlk1];
                a1y = ary[mlk1];
                a1z = arz[mlk1];
                b1x = aix[mlk1];
                b1y = aiy[mlk1];
                b1z = aiz[mlk1];
                a2x = arx[mlk2];
                a2y = ary[mlk2];
                a2z = arz[mlk2];
                b2x = aix[mlk2];
                b2y = aiy[mlk2];
                b2z = aiz[mlk2];
                a3x = arx[mlk3];
                a3y = ary[mlk3];
                a3z = arz[mlk3];
                b3x = aix[mlk3];
                b3y = aiy[mlk3];
                b3z = aiz[mlk3];

                T t1 = a1x*a2*a3 + a1*a2x*a3 + a1*a2*a3x;                
                T t2 = a2x*b1*b3 + a2*b1x*b3 + a2*b1*b3x;
                T t3 = a3x*b1*b2 + a3*b1x*b2 + a3*b1*b2x;
                T t4 = a1x*b2*b3 + a1*b2x*b3 + a1*b2*b3x;
                bd[3*ii+0] += c*(t1 + t2 + t3 - t4);

                t1 = a1y*a2*a3 + a1*a2y*a3 + a1*a2*a3y;                
                t2 = a2y*b1*b3 + a2*b1y*b3 + a2*b1*b3y;
                t3 = a3y*b1*b2 + a3*b1y*b2 + a3*b1*b2y;
                t4 = a1y*b2*b3 + a1*b2y*b3 + a1*b2*b3y;
                bd[3*ii+1] += c*(t1 + t2 + t3 - t4);

                t1 = a1z*a2*a3 + a1*a2z*a3 + a1*a2*a3z;                
                t2 = a2z*b1*b3 + a2*b1z*b3 + a2*b1*b3z;
                t3 = a3z*b1*b2 + a3*b1z*b2 + a3*b1*b2z;
                t4 = a1z*b2*b3 + a1*b2z*b3 + a1*b2*b3z;
                bd[3*ii+2] += c*(t1 + t2 + t3 - t4);                    
            }
        }
        tid += blockDim.x * gridDim.x;   
    }
}
template <typename T> void hipRadialSphericalHarmonicsBispectrumDeriv2(T *bd, T *ar, T *ai, 
        T *arx, T *aix, T *ary, T *aiy, T *arz, T *aiz, T*cg, int *indk, int *indl,
        int *indm, int *rowm, int *Nnb, int Nub, int Ncg, int Na, int K, int npower)
{    
    int K2 = K*(K+1)/2;
    int N2 = K2*Nub;
    int N3 = N2*Na;                        
    int blockDim = 256;
    int gridDim = (N3 + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;
    hipLaunchKernelGGL(hipKernelRadialSphericalHarmonicsBispectrumDeriv2, gridDim, blockDim, 0, 0, bd, ar, ai, arx, aix, 
           ary, aiy, arz, aiz, cg, indk, indl, indm, rowm, Nnb, Nub, Ncg, Na, K, K2, N3, npower);            
}
template void hipRadialSphericalHarmonicsBispectrumDeriv2(double*, double*, double*, double*, 
        double*, double*, double*, double*, double*, double*, int*, int*, int*, int*, int*, int, int, int, int, int);
template void hipRadialSphericalHarmonicsBispectrumDeriv2(float*, float*, float*, float*,  
        float*, float*, float*, float*, float*, float*, int*, int*, int*, int*, int*, int, int, int, int, int);

template <typename T> void hipRadialSphericalHarmonicsSpectrum(T *c, T *ar, T *ai, T *sr, T *si, T *cg, int *indk, 
        int *indl, int *indm, int *rowm, int *Nnb, int Nub, int Ncg, int Na, int L, int K, int spectrum)
{
    //hipSphericalHarmonicsBessel(sr, si, xij, x0, P, tmp, f, fac, M_PI, L, K, N);
    hipRadialSphericalHarmonicsNeighborSum(ar, ai, sr, si, Nnb, Na, L, K);
    
    if (spectrum==0) {  // power spectrum          
        hipRadialSphericalHarmonicsPowerSpectrum(c, ar, ai, indk, Na, L, K);
    }
    else if (spectrum==1) { // bispectrum            
        hipRadialSphericalHarmonicsBispectrum(c, ar, ai, cg, indk, indl, indm, rowm, Nub, Ncg, Na, L, K);
    }
    else if (spectrum==2) { // power spectrum and bispectrum         
        int npower = (L+1)*K*(K+1)/2;      // number of power components
        hipRadialSphericalHarmonicsPowerSpectrum(c, ar, ai, indk, Na, L, K);
        hipRadialSphericalHarmonicsBispectrum(&c[Na*npower], ar, ai, cg, indk, indl, indm, rowm, Nub, Ncg, Na, L, K);            
    }
}
template void hipRadialSphericalHarmonicsSpectrum(double*, double*, double*, double*, double*, double*, int*, int*, 
        int*, int*, int*, int, int, int, int, int, int);
template void hipRadialSphericalHarmonicsSpectrum(float*, float*, float*, float*, float*, float*, int*, int*, 
        int*, int*, int*, int, int, int, int, int, int);
 
template <typename T> void hipRadialSphericalHarmonicsSpectrumDeriv(T *cd, T *ar, T *ai, T *sr, T *si, 
        T *Srx, T *Six, T *Sry, T *Siy, T *Srz, T *Siz, T *cg, int *indk, int *indl, int *indm, 
        int *rowm, int *Nnb, int Nub, int Ncg, int Na, int L, int K, int spectrum)
{    
//     hipSphericalHarmonicsBesselWithDeriv(sr, si, Srx, Six, Sry, Siy, Srz, Siz, xij,
//         x0, P, tmp, f, dP, dtmp, df, fac, M_PI, L, K, N)
    hipRadialSphericalHarmonicsNeighborSum(ar, ai, sr, si, Nnb, Na, L, K);
    
    if (spectrum==0) {  // power spectrum          
        hipRadialSphericalHarmonicsPowerSpectrumDeriv2(cd, ar, ai, Srx, Six, Sry, Siy, Srz, Siz, 
                indk, Nnb, Na, L, K);        
    }
    else if (spectrum==1) { // bispectrum                    
        hipRadialSphericalHarmonicsBispectrumDeriv2(cd, ar, ai, Srx, Six, Sry, Siy, Srz, Siz, 
               cg, indk, indl, indm, rowm, Nnb, Nub, Ncg,  Na, K, (int) 0);                
    }
    else if (spectrum==2) { // power spectrum and bispectrum         
        int npower = (L+1)*K*(K+1)/2;      // number of power components
        hipRadialSphericalHarmonicsPowerSpectrumDeriv2(cd, ar, ai, Srx, Six, Sry, Siy, Srz, Siz, 
                indk, Nnb, Na, L, K);        
        hipRadialSphericalHarmonicsBispectrumDeriv2(cd, ar, ai, Srx, Six, Sry, Siy, Srz, Siz, 
               cg, indk, indl, indm, rowm, Nnb, Nub, Ncg, Na, K, npower);                        
    }
}
template void hipRadialSphericalHarmonicsSpectrumDeriv(double*, double*, double*, double*, double*, double*, double*, 
        double*, double*, double*, double*, double*, int*, int*, int*, int*, int*, int, int, int, int, int, int);
template void hipRadialSphericalHarmonicsSpectrumDeriv(float*, float*, float*, float*, float*, float*, float*,
        float*, float*, float*, float*, float*, int*, int*, int*, int*, int*, int, int, int, int, int, int);

template <typename T> void hipRadialSphericalHarmonicsSpectrumDeriv(T *cd, T *ar, T *ai, 
        T *Srx, T *Six, T *Sry, T *Siy, T *Srz, T *Siz, T *cg, int *indk, int *indl, int *indm, 
        int *rowm, int *Nnb, int Nub, int Ncg, int Na, int L, int K, int spectrum)
{    
//     hipSphericalHarmonicsBesselWithDeriv(sr, si, Srx, Six, Sry, Siy, Srz, Siz, xij,
//         x0, P, tmp, f, dP, dtmp, df, fac, M_PI, L, K, N)    
    if (spectrum==0) {  // power spectrum          
        hipRadialSphericalHarmonicsPowerSpectrumDeriv2(cd, ar, ai, Srx, Six, Sry, Siy, Srz, Siz, 
                indk, Nnb, Na, L, K);        
    }
    else if (spectrum==1) { // bispectrum                    
        hipRadialSphericalHarmonicsBispectrumDeriv2(cd, ar, ai, Srx, Six, Sry, Siy, Srz, Siz, 
               cg, indk, indl, indm, rowm, Nnb, Nub, Ncg,  Na, K, (int) 0);                
    }
    else if (spectrum==2) { // power spectrum and bispectrum         
        int npower = (L+1)*K*(K+1)/2;      // number of power components
        hipRadialSphericalHarmonicsPowerSpectrumDeriv2(cd, ar, ai, Srx, Six, Sry, Siy, Srz, Siz, 
                indk, Nnb, Na, L, K);        
        hipRadialSphericalHarmonicsBispectrumDeriv2(cd, ar, ai, Srx, Six, Sry, Siy, Srz, Siz, 
               cg, indk, indl, indm, rowm, Nnb, Nub, Ncg,  Na, K, npower);                        
    }
}
template void hipRadialSphericalHarmonicsSpectrumDeriv(double*, double*, double*, double*, double*, double*, double*, 
        double*, double*, double*, int*, int*, int*, int*, int*, int, int, int, int, int, int);
template void hipRadialSphericalHarmonicsSpectrumDeriv(float*, float*, float*, float*,  float*, float*, float*,
        float*, float*, float*, int*, int*, int*, int*, int*, int, int, int, int, int, int);

template <typename T> __global__ void hipKernelForceDecomposition(T *f, T *fij, int *ai, int *aj, 
        int inum, int ijnum, int Nbf, int N)
{    
    // fij : 3*ijnum*Nbf 
    // f   : 3*inum*Nbf
    int tid = threadIdx.x + blockIdx.x * blockDim.x;    
    while (tid < N) { 
        int ii = tid%ijnum;     // [1,ijnum]           
        int m = (tid-ii)/ijnum; // [1, Nbf] 
        int n = 3*inum*m;
        int p = 3*ijnum*m;
        int i = 3*ai[ii] + n;       // atom i
        int j = 3*aj[ii] + n;       // atom j
        int l = 3*ii + p;        
        // use atomicAdd on hip
        atomicAdd(&f[0+i], -fij[0+l]);
        atomicAdd(&f[1+i], -fij[1+l]);
        atomicAdd(&f[2+i], -fij[2+l]);                
        atomicAdd(&f[0+j],  fij[0+l]);
        atomicAdd(&f[1+j],  fij[1+l]);
        atomicAdd(&f[2+j],  fij[2+l]);   
        tid += blockDim.x * gridDim.x;           
    }
}
template <typename T> void hipForceDecomposition(T *f, T *fij, int *ai, int *aj, 
        int inum, int ijnum, int Nbf)
{    
    int N = ijnum*Nbf;
    int blockDim = 256;
    int gridDim = (N + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;
    hipLaunchKernelGGL(hipKernelForceDecomposition, gridDim, blockDim, 0, 0, f, fij, ai, aj, inum, ijnum, Nbf, N);  
}
template void hipForceDecomposition(double*, double*, int*, int*, int, int, int);
template void hipForceDecomposition(float*, float*, int*, int*, int, int, int);

template <typename T> __global__ void hipKernelForceDecomposition(T *f, T *fij, int *ai, int *aj, int ijnum)
{    
    // fij : 3*ijnum 
    // f   : 3*inum
    int ii = threadIdx.x + blockIdx.x * blockDim.x;    
    while (ii < ijnum) { // for each atom pair ij in the simulation box     
        int i = 3*ai[ii];       // atom i
        int j = 3*aj[ii];       // atom j
        int l = 3*ii; 
        // use atomicAdd on hip
        atomicAdd(&f[0+i], -fij[0+l]);
        atomicAdd(&f[1+i], -fij[1+l]);
        atomicAdd(&f[2+i], -fij[2+l]);                
        atomicAdd(&f[0+j],  fij[0+l]);
        atomicAdd(&f[1+j],  fij[1+l]);
        atomicAdd(&f[2+j],  fij[2+l]);   
        ii += blockDim.x * gridDim.x;   
    }
}
template <typename T> void hipForceDecomposition(T *f, T *fij, int *ai, int *aj, int ijnum)
{    
    int blockDim = 256;
    int gridDim = (ijnum + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;
    hipLaunchKernelGGL(hipKernelForceDecomposition, gridDim, blockDim, 0, 0, f, fij, ai, aj, ijnum);  
}
template void hipForceDecomposition(double*, double*, int*, int*, int);
template void hipForceDecomposition(float*, float*, int*, int*, int);

template <typename T> __global__ void hipKernelCenterAtomDecomposition(T *e, T *ei, int *ilist, int Na)
{    
    int ii = threadIdx.x + blockIdx.x * blockDim.x;    
    while (ii < Na) {  // for each atom i in the simulation box     
        int i = ilist[ii];       // atom i        
        e[i] += ei[ii];             
        ii += blockDim.x * gridDim.x;   
    }    
}
template <typename T> void hipCenterAtomDecomposition(T *e, T *ei, int *ilist, int Na)
{    
    int blockDim = 256;
    int gridDim = (Na + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;
     hipLaunchKernelGGL(hipKernelCenterAtomDecomposition, gridDim, blockDim, 0, 0, e, ei, ilist, Na);  
}
template void hipCenterAtomDecomposition(double*, double*, int*, int);
template void hipCenterAtomDecomposition(float*, float*, int*, int);

template <typename T> __global__ void hipKernelCenterAtomDecomposition(T *f, T *fij, int *ilist, int *anumsum, 
        int inum, int ijnum, int Na, int Nbf, int N)
{    
    int tid = threadIdx.x + blockIdx.x * blockDim.x;    
    while (tid < N) {          
        int ii = tid%Na;     // [1, Na]           
        int m = (tid-ii)/Na; // [1, Nbf] 
        int n = 3*inum*m;
        int p = 3*ijnum*m;
        int i = ilist[ii];       // atom i        
        int start = anumsum[ii];   
        int k = anumsum[ii+1]-start;    // number of neighbors around i             
        int idim = 3*i + n;        
        for (int l=0; l<k ; l++) {   // loop over each atom j around atom i
            int ndim = 3*(start + l) + p;                
            f[0+idim] += -fij[0+ndim]; 
            f[1+idim] += -fij[1+ndim]; 
            f[2+idim] += -fij[2+ndim]; 
        }                
        tid += blockDim.x * gridDim.x;    
    }
}
template <typename T> void hipCenterAtomDecomposition(T *f, T *fij, int *ilist, int *anumsum, 
        int inum, int ijnum, int Na, int Nbf)
{    
    int N = Na*Nbf;
    int blockDim = 256;
    int gridDim = (N + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;
     hipLaunchKernelGGL(hipKernelCenterAtomDecomposition, gridDim, blockDim, 0, 0, f, fij, ilist, anumsum, 
             inum, ijnum, Na, Nbf, N);  
}
template void hipCenterAtomDecomposition(double*, double*, int*, int*, int, int, int, int);
template void hipCenterAtomDecomposition(float*, float*, int*, int*, int, int, int, int);

template <typename T> __global__ void hipKernelNeighborAtomDecomposition(T *f, T *fij, int *jlist, int *bnumsum, int *index, 
        int inum, int ijnum, int jnum, int Nbf, int N)
{        
    int tid = threadIdx.x + blockIdx.x * blockDim.x;    
    while (tid < N) {          
        int ii = tid%jnum;     // [1, jnum]           
        int m = (tid-ii)/jnum; // [1, Nbf]
        int n = 3*inum*m;
        int p = 3*ijnum*m;    
        int j = jlist[ii];       // atom j        
        int start = bnumsum[ii];   
        int k = bnumsum[ii+1]-start;        // number of neighbors around j             
        int jdim = 3*j + n;        
        for (int l=0; l<k ; l++) {   // loop over each atom l around atom j                
            int ndim = 3*index[start + l] + p;     
            f[0+jdim] += fij[0+ndim]; 
            f[1+jdim] += fij[1+ndim]; 
            f[2+jdim] += fij[2+ndim];                 
        }                        
        tid += blockDim.x * gridDim.x;                    
    }
}
template <typename T> void hipNeighborAtomDecomposition(T *f, T *fij, int *jlist, int *bnumsum, 
        int *index, int inum, int ijnum, int jnum, int Nbf)
{    
    int N = jnum*Nbf;
    int blockDim = 256;
    int gridDim = (N + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;
     hipLaunchKernelGGL(hipKernelNeighborAtomDecomposition, gridDim, blockDim, 0, 0, f, fij, jlist, bnumsum, 
             index, inum, ijnum, jnum, Nbf, N);  
}
template void hipNeighborAtomDecomposition(double*, double*, int*, int*, int*, int, int, int, int);
template void hipNeighborAtomDecomposition(float*, float*, int*, int*, int*, int, int, int, int);

template <typename T> __global__ void hipKernelCenterAtomDecomposition(T *f, T *fij, int *ilist, int *anumsum, int inum)
{    
    int ii = threadIdx.x + blockIdx.x * blockDim.x;    
    while (ii < inum) { // for each atom pair ij in the simulation box                 
        int i = ilist[ii];       // atom i        
        int start = anumsum[ii];   
        int k = anumsum[ii+1]-start;    // number of neighbors around i             
        int idim = 3*i;        
        for (int l=0; l<k ; l++) {   // loop over each atom j around atom i
            int ndim = 3*(start + l);                
            f[0+idim] += -fij[0+ndim]; 
            f[1+idim] += -fij[1+ndim]; 
            f[2+idim] += -fij[2+ndim]; 
        }        
        ii += blockDim.x * gridDim.x;   
    }
}
template <typename T> void hipCenterAtomDecomposition(T *f, T *fij, int *ilist, int *anumsum, int inum)
{    
    int blockDim = 256;
    int gridDim = (inum + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;
     hipLaunchKernelGGL(hipKernelCenterAtomDecomposition, gridDim, blockDim, 0, 0, f, fij, ilist, anumsum, inum);  
}
template void hipCenterAtomDecomposition(double*, double*, int*, int*, int);
template void hipCenterAtomDecomposition(float*, float*, int*, int*, int);

template <typename T> __global__ void hipKernelNeighborAtomDecomposition(T *f, T *fij, int *jlist, int *bnumsum, int *index, int jnum)
{        
    int ii = threadIdx.x + blockIdx.x * blockDim.x;    
    while (ii < jnum) { // for each atom pair ij in the simulation box                 
        int j = jlist[ii];       // atom j        
        int start = bnumsum[ii];   
        int k = bnumsum[ii+1]-start;        // number of neighbors around j             
        int jdim = 3*j;        
        for (int l=0; l<k ; l++) {   // loop over each atom l around atom j                
            int ndim = 3*index[start + l];     
            f[0+jdim] += fij[0+ndim]; 
            f[1+jdim] += fij[1+ndim]; 
            f[2+jdim] += fij[2+ndim];                 
        }        
        ii += blockDim.x * gridDim.x;   
    }    
}
template <typename T> void hipNeighborAtomDecomposition(T *f, T *fij, int *jlist, int *bnumsum, int *index, int jnum)
{    
    int blockDim = 256;
    int gridDim = (jnum + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;
     hipLaunchKernelGGL(hipKernelNeighborAtomDecomposition, gridDim, blockDim, 0, 0, f, fij, jlist, bnumsum, index, jnum);  
}
template void hipNeighborAtomDecomposition(double*, double*, int*, int*, int*, int);
template void hipNeighborAtomDecomposition(float*, float*, int*, int*, int*, int);


template <typename T> void hipNeighborAtomDecomposition(T *c, T *b, T *b_sorted, int *a, int *a_sorted,
        int *ia, int *index, int *d_temp, int n)
{           
    // Inputs: a is an integer array of size n
    //         b is a float array of size n
    // output: c is a float array of size m such that c[a[i]] = c[a[i]] + b[i] for i = 1,2,...,n          
            
    // ia = [0,1,2,...,n]        
    hipIndexInit(ia, n);

    // sort array a
    size_t  temp_storage_bytes  = 0;
    void  *d_temp_storage = NULL;
    hipcub::DeviceRadixSort::SortPairs(d_temp_storage, temp_storage_bytes, a, a_sorted, ia, index, n);                
    hipcub::DeviceRadixSort::SortPairs((void*) d_temp, temp_storage_bytes, a, a_sorted, ia, index, n);
    
    //b_sorted[index[i]] = b[i] for i = 1,2,...,n
    hipPutArrayAtIndex(b_sorted, b, index, n);

    // a contains unique elements of a_sorted
    // index[0] is the number of unique elements        
    // ia contains the counts for unique elements    
    temp_storage_bytes  = 0;
    hipcub::DeviceRunLengthEncode::Encode(d_temp_storage, temp_storage_bytes, a_sorted, a, ia, index, n);                
    hipcub::DeviceRunLengthEncode::Encode((void*) d_temp, temp_storage_bytes, a_sorted, a, ia, index, n);

    // number of unique elements
    int m = hipArrayGetValueAtIndex(index, 0);

    // index is inclusive scan of ia
    temp_storage_bytes  = 0;
    hipcub::DeviceScan::InclusiveSum(d_temp_storage, temp_storage_bytes, ia, index, m);                
    hipcub::DeviceScan::InclusiveSum((void*) d_temp, temp_storage_bytes, ia, index, m);
    
    //  b is segmented reduced sum of b_sorted 
    temp_storage_bytes  = 0;
    hipcub::DeviceSegmentedReduce::Sum(d_temp_storage, temp_storage_bytes, b_sorted, b, m, index, index+1);                
    hipcub::DeviceSegmentedReduce::Sum((void*) d_temp, temp_storage_bytes, b_sorted, b, m, index, index+1);

    // c[a[i]] = c[a[i]] + b[i] for i = 1,2,...,m
    hipPutArrayAtIndex(c, b, a, m);
    

//     template <
//         typename            InputIteratorT,
//         typename            OutputIteratorT>
//     CUB_RUNTIME_FUNCTION
//     static hipError_t InclusiveSum(
//         void*               d_temp_storage,                 ///< [in] %Device-accessible allocation of temporary storage.  When NULL, the required allocation size is written to \p temp_storage_bytes and no work is done.
//         size_t&             temp_storage_bytes,             ///< [in,out] Reference to size in bytes of \p d_temp_storage allocation
//         InputIteratorT      d_in,                           ///< [in] Pointer to the input sequence of data items
//         OutputIteratorT     d_out,                          ///< [out] Pointer to the output sequence of data items
//         int                 num_items,                      ///< [in] Total number of input items (i.e., the length of \p d_in)
//         hipStream_t        stream             = 0,         ///< [in] <b>[optional]</b> CUDA stream to launch kernels within.  Default is stream<sub>0</sub>.
//         bool                debug_synchronous  = false)     ///< [in] <b>[optional]</b> Whether or not to synchronize the stream after every kernel launch to check for errors.  May cause significant slowdown.  Default is \p false.
    
//     static hipError_t Sum(
//         void                *d_temp_storage,                    ///< [in] %Device-accessible allocation of temporary storage.  When NULL, the required allocation size is written to \p temp_storage_bytes and no work is done.
//         size_t              &temp_storage_bytes,                ///< [in,out] Reference to size in bytes of \p d_temp_storage allocation
//         InputIteratorT      d_in,                               ///< [in] Pointer to the input sequence of data items
//         OutputIteratorT     d_out,                              ///< [out] Pointer to the output aggregate
//         int                 num_segments,                       ///< [in] The number of segments that comprise the sorting data
//         OffsetIteratorT     d_begin_offsets,                    ///< [in] Pointer to the sequence of beginning offsets of length \p num_segments, such that <tt>d_begin_offsets[i]</tt> is the first element of the <em>i</em><sup>th</sup> data segment in <tt>d_keys_*</tt> and <tt>d_values_*</tt>
//         OffsetIteratorT     d_end_offsets,                      ///< [in] Pointer to the sequence of ending offsets of length \p num_segments, such that <tt>d_end_offsets[i]-1</tt> is the last element of the <em>i</em><sup>th</sup> data segment in <tt>d_keys_*</tt> and <tt>d_values_*</tt>.  If <tt>d_end_offsets[i]-1</tt> <= <tt>d_begin_offsets[i]</tt>, the <em>i</em><sup>th</sup> is considered empty.
//         hipStream_t        stream              = 0,            ///< [in] <b>[optional]</b> CUDA stream to launch kernels within.  Default is stream<sub>0</sub>.
//         bool                debug_synchronous   = false)        ///< [in] <b>[optional]</b> Whether or not to synchronize the stream after every kernel launch to check for errors.  Also causes launch configurations to be printed to the console.  Default is \p false.
        
//     static hipError_t Encode(
//         void*                       d_temp_storage,                ///< [in] %Device-accessible allocation of temporary storage.  When NULL, the required allocation size is written to \p temp_storage_bytes and no work is done.
//         size_t                      &temp_storage_bytes,            ///< [in,out] Reference to size in bytes of \p d_temp_storage allocation
//         InputIteratorT              d_in,                           ///< [in] Pointer to the input sequence of keys
//         UniqueOutputIteratorT       d_unique_out,                   ///< [out] Pointer to the output sequence of unique keys (one key per run)
//         LengthsOutputIteratorT      d_counts_out,                   ///< [out] Pointer to the output sequence of run-lengths (one count per run)
//         NumRunsOutputIteratorT      d_num_runs_out,                     ///< [out] Pointer to total number of runs
//         int                         num_items,                      ///< [in] Total number of associated key+value pairs (i.e., the length of \p d_in_keys and \p d_in_values)
//         hipStream_t                stream             = 0,         ///< [in] <b>[optional]</b> CUDA stream to launch kernels within.  Default is stream<sub>0</sub>.
//         bool                        debug_synchronous  = false)     ///< [in] <b>[optional]</b> Whether or not to synchronize the stream after every kernel launch to check for errors.  May cause significant slowdown.  Default is \p false.    
}
template void hipNeighborAtomDecomposition(double*, double*, double*, int*, int*, int*, int*, int*, int);
template void hipNeighborAtomDecomposition(float*, float*, float*, int*, int*, int*, int*, int*, int);


#endif

