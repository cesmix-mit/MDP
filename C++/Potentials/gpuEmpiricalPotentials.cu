/***************************************************************************                               
                    Molecular Dynamics Potentials (MDP)
                           CESMIX-MIT Project  
 
 Contributing authors: Ngoc-Cuong Nguyen (cuongng@mit.edu, exapde@gmail.com)
 ***************************************************************************/

#ifndef __GPUEMPIRICALPOTENTIALS
#define __GPUEMPIRICALPOTENTIALS

template <typename... Args>
__device__ void __enzyme_autodiff(void*, Args... args);
int __device__ enzyme_const;
int __device__ enzyme_dup;

#include "gpuSinglea.cu"
#include "gpuSingleb.cu"
#include "gpuPaira.cu"
#include "gpuPairb.cu"
#include "gpuPairc.cu"
#include "gpuTripleta.cu"
#include "gpuTripletb.cu"
#include "gpuTripletc.cu"
#include "gpuQuadrupleta.cu"
#include "gpuQuadrupletb.cu"

template <typename T> __global__ void gpuKernelElectronDensity(T *rhoi, T *rhoij, int *pairnum, 
        int *pairnumsum, int inum) 
{
    int i = threadIdx.x + blockIdx.x * blockDim.x;
	while (i<inum) {    
        int jnum = pairnum[i];
        int start = pairnumsum[i];
        rhoi[i] = 0.0;
        for (int j=0; j<jnum; j++) 
            rhoi[i] += rhoij[start+j];        
        i += blockDim.x * gridDim.x;
    }
}
template <typename T> void gpuElectronDensity(T *rhoi, T *rhoij, int *pairnum, 
        int *pairnumsum, int inum) 
{
	int blockDim = 256;
	int gridDim = (inum + blockDim - 1) / blockDim;
	gridDim = (gridDim>1024)? 1024 : gridDim;
	gpuKernelElectronDensity<<<gridDim, blockDim>>>(rhoi, rhoij, pairnum, pairnumsum, inum);
}
template void gpuElectronDensity(double *, double *, int *, int *, int);
template void gpuElectronDensity(float *, float *, int *, int *, int);

template <typename T> __global__ void gpuKernelEmbedingForce(T *fij, T *d_rhoi, int *pairnum, 
        int *pairnumsum, int inum)
{    
    int i = threadIdx.x + blockIdx.x * blockDim.x;
	while (i<inum) {    
        int jnum = pairnum[i];
        int start = pairnumsum[i];
        for (int j=0; j<jnum; j++)             
            fij[start+j] = d_rhoi[i]*fij[start+j];      
        i += blockDim.x * gridDim.x;
    }        
}
template <typename T> void gpuEmbedingForce(T *fij, T *d_rhoi, int *pairnum, 
        int *pairnumsum, int inum) 
{
	int blockDim = 256;
	int gridDim = (inum + blockDim - 1) / blockDim;
	gridDim = (gridDim>1024)? 1024 : gridDim;
	gpuKernelEmbedingForce<<<gridDim, blockDim>>>(fij, d_rhoi, pairnum, pairnumsum, inum);
}
template void gpuEmbedingForce(double *, double *, int *, int *, int);
template void gpuEmbedingForce(float *, float *, int *, int *, int);
 
//template <typename T> __device__ void gpuKernelLJ_device(T * u, T * xij, T *  qi, T * qj, 
//   int * ti, int * tj, int * ai, int * aj, T * mu, T * eta, int * kappa, int dim, int ncq, 
//   int nmu, int neta, int nkappa, int ng)
template <typename T> __device__ void gpuKernelLJ_device(T * __restrict__ u, T *__restrict__ xij,
 T* __restrict__ qi, T *__restrict__ qj, 
        int *__restrict__ ti, int *__restrict__ tj, int *__restrict__ ai, int *__restrict__ aj,
 T *__restrict__ mu, T *__restrict__ eta, int *__restrict__ kappa, int dim, int ncq, 
                int nmu, int neta, int nkappa, int ng)
{
	int i = threadIdx.x + blockIdx.x * blockDim.x;

	while (i<ng) {
		T mu1 = mu[0];
		T mu2 = mu[1];
		T xij1 = xij[0 + i*3];
		T xij2 = xij[1 + i*3];
		T xij3 = xij[2 + i*3];
		T t2 = xij1*xij1;
		T t3 = xij2*xij2;
		T t4 = xij3*xij3;
		T t5 = t2+t3+t4; // rij^2
        T t6 = t5*t5*t5; // rij^6
        T t8 = t6*t5;    // rij^8
        T t12 = t6*t6;   // rij^12
        T t14 = t12*t5;  // rij^14
		u[i] = mu1*1.0/(t12) - mu2*1.0/(t6);
        i += blockDim.x * gridDim.x;
	}	
}

template <typename T> __global__ void gradient_gpuKernelLJ(T *u, T *d_u, T *xij, T *u_x, T *qi, T *qj, 
        int *ti, int *tj, int *ai, int *aj, T *mu, T *eta, int *kappa, int dim, int ncq, 
                int nmu, int neta, int nkappa, int ng)
{
    __enzyme_autodiff((void*)gpuKernelLJ_device<T>, 
                        enzyme_dup, u, d_u,
                        enzyme_dup, xij, u_x,
                        enzyme_const, qi,
                        enzyme_const, qj,
                        enzyme_const, ti,
                        enzyme_const, tj,
                        enzyme_const, ai,
                        enzyme_const, aj,
                        // TODO can make this non const to get du/dmu
                        enzyme_const, mu,
                        enzyme_const, eta,
                        enzyme_const, kappa,
                        dim, ncq, nmu, neta, nkappa, ng);  
  //if (threadIdx.x == 0 && blockIdx.x == 0)
  //  printf("Called gpuGradient u_x[0]=%f d_u[0]=%f\n", u_x[0], d_u[0]);

}
template <typename T> void gpuGradientLJ(T *u, T *du, T *xij, T *u_x, T *qi, T *qj, int *ti, int *tj, 
        int *ai, int *aj, T *mu, T *eta, int *kappa, int dim, int ncq, 
                int nmu, int neta, int nkappa, int ng)
{
	int blockDim = 256;
	int gridDim = (ng + blockDim - 1) / blockDim;
	gridDim = (gridDim>1024)? 1024 : gridDim;
	gradient_gpuKernelLJ<<<gridDim, blockDim>>>(u, du, xij, u_x, qi, qj, ti, tj, ai, aj, mu, eta, kappa, 
            dim, ncq, nmu, neta, nkappa, ng);
}
template void gpuGradientLJ(double *, double *, double *, double *, double *, double *, int *, int *, 
int *, int *, double *, double *, int*, int, int, int, int, int, int);
template void gpuGradientLJ(float *, float *, float *, float *, float *, float *, int *, int *, 
int *, int *, float *, float *, int *, int, int, int, int, int, int);


template <typename T> __global__ void gpuKernelLJ(T *u, T *f, T *xij, T *qi, T *qj, 
        int *ti, int *tj, int *ai, int *aj, T *mu, T *eta, int *kappa, int dim, int ncq, 
                int nmu, int neta, int nkappa, int ng)
{
	int i = threadIdx.x + blockIdx.x * blockDim.x;
	while (i<ng) {
		T mu1 = mu[0];
		T mu2 = mu[1];
		T xij1 = xij[0 + i*3];
		T xij2 = xij[1 + i*3];
		T xij3 = xij[2 + i*3];
		T t2 = xij1*xij1;
		T t3 = xij2*xij2;
		T t4 = xij3*xij3;
		T t5 = t2+t3+t4; // rij^2
        T t6 = t5*t5*t5; // rij^6
        T t8 = t6*t5;    // rij^8
        T t12 = t6*t6;   // rij^12
        T t14 = t12*t5;  // rij^14
		u[i] = mu1*1.0/(t12) - mu2*1.0/(t6);
        f[0 + 3*i] = (12*mu1*xij1)/t14 - (6*mu2*xij1)/t8;
        f[1 + 3*i] = (12*mu1*xij2)/t14 - (6*mu2*xij2)/t8;
        f[2 + 3*i] = (12*mu1*xij3)/t14 - (6*mu2*xij3)/t8;
        i += blockDim.x * gridDim.x;
	}	
}

template <typename T> void gpuLJ(T *u, T *f, T *xij, T *qi, T *qj, int *ti, int *tj, 
        int *ai, int *aj, T *mu, T *eta, int *kappa, int dim, int ncq, 
                int nmu, int neta, int nkappa, int ng)
{
	int blockDim = 256;
	int gridDim = (ng + blockDim - 1) / blockDim;
	gridDim = (gridDim>1024)? 1024 : gridDim;
	gpuKernelLJ<<<gridDim, blockDim>>>(u, f, xij, qi, qj, ti, tj, ai, aj, mu, eta, kappa, 
            dim, ncq, nmu, neta, nkappa, ng);
}
template void gpuLJ(double *, double *, double *, double *, double *, int *, int *, 
int *, int *, double *, double *, int*, int, int, int, int, int, int);
template void gpuLJ(float *, float *, float *, float *, float *, int *, int *, 
int *, int *, float *, float *, int *, int, int, int, int, int, int);


#endif
