/***************************************************************************                               
                    Molecular Dynamics Potentials (MDP)
                           CESMIX-MIT Project  
 
 Contributing authors: Ngoc-Cuong Nguyen (cuongng@mit.edu, exapde@gmail.com)
 ***************************************************************************/

#include "hip/hip_runtime.h"

#ifndef __HIPCOORDINATETRANSFORMATIONS
#define __HIPCOORDINATETRANSFORMATIONS

template <typename T> __global__ void hipKernelCart2Sphere(T *the, T *phi, T *r, T *x, T *y, T *z, int N)
{    
    //#pragma omp parallel for
    for (int i=0; i<N; i++) {
        r[i] = sqrt(x[i]*x[i] + y[i]*y[i] + z[i]*z[i]);
        the[i] = acos(z[i]/r[i]);
        phi[i] = atan2(y[i],x[i]);
    }
}
template <typename T> void hipCart2Sphere(T *the, T *phi, T *r, T *x, T *y, T *z, int N)
{        
    int blockDim = 256;
    int gridDim = (N + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;
    hipLaunchKernelGGL(hipKernelCart2Sphere, gridDim, blockDim, 0, 0, the, phi, r, x, y, z, N);
}
template void hipCart2Sphere(double*, double*, double*, double*, double*, double*, int);
template void hipCart2Sphere(float*, float*, float*, float*, float*, float*, int);

template <typename T> __global__ void hipKernelCart2SphereDeriv(T *the, T *phi, T *r, T *thex, T *they, T *thez, T *phix, T *phiy, T *phiz, T *rx, T *ry, T *rz, T *x, T *y, T *z, int N)
{    
    //#pragma omp parallel for
    for (int i=0; i<N; i++) {
        r[i] = sqrt(x[i]*x[i] + y[i]*y[i] + z[i]*z[i]);
        the[i] = acos(z[i]/r[i]);
        phi[i] = atan2(y[i],x[i]);
        
        T r2 = r[i]*r[i];
        T rxy = sqrt(x[i]*x[i] + y[i]*y[i]);
        T rxy2 = rxy*rxy;
        T rr2 = rxy*r2;

        rx[i] = x[i]/r[i];
        ry[i] = y[i]/r[i];
        rz[i] = z[i]/r[i];
        thex[i] = x[i]*z[i]/rr2;
        they[i] = y[i]*z[i]/rr2;
        thez[i] = -rxy/r2;
        phix[i] = -y[i]/rxy2;
        phiy[i] = x[i]/rxy2;
        phiz[i] = 0.0;        
    }
}
template <typename T> void hipCart2SphereDeriv(T *the, T *phi, T *r, T *thex, T *they, T *thez, 
        T *phix, T *phiy, T *phiz, T *rx, T *ry, T *rz, T *x, T *y, T *z, int N)
{        
    int blockDim = 256;
    int gridDim = (N + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;
    hipLaunchKernelGGL(hipKernelCart2SphereDeriv, gridDim, blockDim, 0, 0, the, phi, r, thex, they, thez, phix, phiy, phiz,
            rx, ry, rz, x, y, z, N);
}
template void hipCart2SphereDeriv(double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, int);
template void hipCart2SphereDeriv(float*, float*, float*, float*, float*, float*, float*, float*, float*, float*, float*, float*, float*, float*, float*, int);

template <typename T> __global__ void hipKernelSphere2Cart(T *x, T *y, T *z, T *the, T *phi, T *r, int N)
{    
    //#pragma omp parallel for
    for (int i=0; i<N; i++) {
        x[i] = r[i]*sin(the[i])*cos(phi[i]);
        y[i] = r[i]*sin(the[i])*sin(phi[i]);
        z[i] = r[i]*cos(the[i]);
    }
}
template <typename T> void hipSphere2Cart(T *x, T *y, T *z, T *the, T *phi, T *r, int N)
{        
    int blockDim = 256;
    int gridDim = (N + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;
    hipLaunchKernelGGL(hipKernelSphere2Cart, gridDim, blockDim, 0, 0, x, y, z, the, phi, r, N);
}
template void hipSphere2Cart(double*, double*, double*, double*, double*, double*, int);
template void hipSphere2Cart(float*, float*, float*, float*, float*, float*, int);

template <typename T> __global__ void hipKernelEuler2Rotm(T *R11, T *R12, T *R13, T *R21, 
                T *R22, T *R23, T *R31, T *R32, T *R33, T *alpha, T *beta, T *gamma, int N)
{            
    //#pragma omp parallel for
    for (int i=0; i<N; i++) {
        T ca = cos(alpha[i]);
        T cb = cos(beta[i]);
        T cg = cos(gamma[i]);
        T sa = sin(alpha[i]);
        T sb = sin(beta[i]);
        T sg = sin(gamma[i]);
        
        R11[i] = ca*cg*cb - sa*sg;
        R12[i] = -ca*cb*sg - sa*cg;
        R13[i] = ca*sb;
        R21[i] = sa*cg*cb + ca*sg;
        R22[i] = -sa*cb*sg + ca*cg;
        R23[i] = sa*sb;
        R31[i] = -sb*cg;
        R32[i] = sb*sg;
        R33[i] = cb;
    }        
}
template <typename T> void hipEuler2Rotm(T *R11, T *R12, T *R13, T *R21, 
                T *R22, T *R23, T *R31, T *R32, T *R33, T *alpha, T *beta, T *gamma, int N)
{        
    int blockDim = 256;
    int gridDim = (N + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;
    hipLaunchKernelGGL(hipKernelEuler2Rotm, gridDim, blockDim, 0, 0, R11, R12, R13, R21, R22, R23, R31, R32, R33, alpha, beta, gamma, N);
}
template void hipEuler2Rotm(double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, int);
template void hipEuler2Rotm(float*, float*, float*, float*, float*, float*, float*, float*, float*, float*, float*, float*, int);

template <typename T> __global__ void hipKernelRotc(T *X, T *Y, T *Z, T *R, T *x, T *y, T *z, int N)
{            
    //#pragma omp parallel for
    for (int i=0; i<N; i++) {        
        X[i] = R[0]*x[i] + R[3]*y[i] + R[6]*z[i];
        Y[i] = R[1]*x[i] + R[4]*y[i] + R[7]*z[i];
        Z[i] = R[2]*x[i] + R[5]*y[i] + R[8]*z[i];
    }        
}
template <typename T> void hipRotc(T *X, T *Y, T *Z, T *R, T *x, T *y, T *z, int N)
{        
    int blockDim = 256;
    int gridDim = (N + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;
    hipLaunchKernelGGL(hipKernelRotc, gridDim, blockDim, 0, 0, X, Y, Z, R, x, y, z, N);
}
template void hipRotc(double*, double*, double*, double*, double*, double*, double*, int);
template void hipRotc(float*, float*, float*, float*, float*, float*, float*, int);

#endif

