/***************************************************************************                               
                    Molecular Dynamics Potentials (MDP)
                           CESMIX-MIT Project  
 
 Contributing authors: Ngoc-Cuong Nguyen (cuongng@mit.edu, exapde@gmail.com)
 ***************************************************************************/

#include "hip/hip_runtime.h"

#ifndef __HIPFORCEDECOMPOSITION
#define __HIPFORCEDECOMPOSITION

template <typename T> __global__ void hipKernelSingleDecomposition(T *e, T *ei, int *ai, int inum, int dim)
{    
    int ii = threadIdx.x + blockIdx.x * blockDim.x;         
    while (ii < inum) { 
        int i = ai[ii];       // atom i        
        int l = dim*ii;                
        e[i] += ei[ii];
        ii += blockDim.x * gridDim.x;   
    }
}
template <typename T> void hipSingleDecomposition(T *e, T *ei, int *ai, int inum, int dim)
{                        
    int blockDim = 256;
    int gridDim = (inum + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;
      
    hipLaunchKernelGGL(hipKernelSingleDecomposition, gridDim, blockDim, 0, 0, e, ei, ai, inum, dim);
}
template void hipSingleDecomposition(double*, double*, int*, int, int);
template void hipSingleDecomposition(float*, float*, int*, int, int);

template <typename T> __global__ void hipKernelFullNeighPairDecomposition(T *e, T *eij, int *ai, int ijnum, int dim)
{    
    int ii = threadIdx.x + blockIdx.x * blockDim.x;         
    while (ii < ijnum)  {  // for each atom pair ij in the simulation box     
        int i = ai[ii];       // atom i
        int l = dim*ii;
        // use atomicAdd on hip
        //e[i] += 0.5*eij[ii];
        atomicAdd(&e[i], 0.5*eij[ii]);
        ii += blockDim.x * gridDim.x;   
    }
}
template <typename T> void hipFullNeighPairDecomposition(T *e, T *eij, int *ai, int ijnum, int dim)
{                        
    int blockDim = 256;
    int gridDim = (ijnum + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;
      
    hipLaunchKernelGGL(hipKernelFullNeighPairDecomposition, gridDim, blockDim, 0, 0, e, eij, ai, ijnum, dim);
}
template void hipFullNeighPairDecomposition(double*, double*, int*, int, int);
template void hipFullNeighPairDecomposition(float*, float*, int*, int, int);

template <typename T> __global__ void hipKernelCenterAtomPairDecomposition(T *e, T *eij, int *ilist, 
        int *anumsum, int inum, int dim)
{    
    int ii = threadIdx.x + blockIdx.x * blockDim.x;         
    while (ii < inum)  {  // for each atom i in the simulation box     
        int i = ilist[ii];       // atom i        
        int start = anumsum[ii];   
        int m = anumsum[ii+1]-start;        // number of neighbors around i             
        int idim = dim*i;        
        for (int l=0; l<m ; l++)    // loop over each atom around atom i            
            e[i] += 0.5*eij[start + l];             
        ii += blockDim.x * gridDim.x;   
    }
}
template <typename T> void hipCenterAtomPairDecomposition(T *e, T *eij, int *ilist, 
        int *anumsum, int inum, int dim)
{                        
    int blockDim = 256;
    int gridDim = (inum + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;
      
    hipLaunchKernelGGL(hipKernelCenterAtomPairDecomposition, gridDim, blockDim, 0, 0, e, eij, ilist, anumsum, inum, dim);
}
template void hipCenterAtomPairDecomposition(double*, double*, int*, int*, int, int);
template void hipCenterAtomPairDecomposition(float*, float*, int*, int*, int, int);

template <typename T> __global__ void hipKernelHalfNeighPairDecomposition(T *e, T *eij, int *ai, 
        int *aj, int ijnum, int dim)
{    
    int ii = threadIdx.x + blockIdx.x * blockDim.x;         
    while (ii < ijnum)  {  // for each atom pair ij in the simulation box     
        int i = ai[ii];       // atom i
        int j = aj[ii];       // atom j
        int l = dim*ii;        
        // use atomicAdd on hip        
        //e[i] += 0.5*eij[ii];
        //e[j] += 0.5*eij[ii];
        atomicAdd(&e[i], 0.5*eij[ii]);                
        atomicAdd(&e[j], 0.5*eij[ii]);                
        ii += blockDim.x * gridDim.x;   
    }
}
template <typename T> void hipHalfNeighPairDecomposition(T *e, T *eij, int *ai, 
        int *aj, int ijnum, int dim)
{                        
    int blockDim = 256;
    int gridDim = (ijnum + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;
      
    hipLaunchKernelGGL(hipKernelHalfNeighPairDecomposition, gridDim, blockDim, 0, 0, e, eij, ai, aj, ijnum, dim);
}
template void hipHalfNeighPairDecomposition(double*, double*, int*, int*, int, int);
template void hipHalfNeighPairDecomposition(float*, float*, int*, int*, int, int);

template <typename T> __global__ void hipKernelNeighborAtomPairDecomposition(T *e, T *eij, int *jlist, 
        int *bnumsum, int *index, int jnum, int dim)
{        
    int ii = threadIdx.x + blockIdx.x * blockDim.x;         
    while (ii < jnum)  {  // for each atom j in the simulation box     
        int j = jlist[ii];       // atom j        
        int start = bnumsum[ii];   
        int m = bnumsum[ii+1]-start;        // number of neighbors around i             
        for (int l=0; l<m ; l++) {   // loop over each atom around atom i
            int n = index[start + l];
            e[j] += 0.5*eij[n];                                
        }        
        ii += blockDim.x * gridDim.x;   
    }    
}
template <typename T> void hipNeighborAtomPairDecomposition(T *e, T *eij, int *jlist, 
        int *bnumsum, int *index, int jnum, int dim)
{                        
    int blockDim = 256;
    int gridDim = (jnum + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;
      
    hipLaunchKernelGGL(hipKernelNeighborAtomPairDecomposition, gridDim, blockDim, 0, 0, e, eij, jlist, bnumsum, index, jnum, dim);
}
template void hipNeighborAtomPairDecomposition(double*, double*, int*, int*, int*, int, int);
template void hipNeighborAtomPairDecomposition(float*, float*, int*, int*, int*, int, int);

template <typename T> __global__ void hipKernelTripletDecomposition(T *e, T *eijk, int *ai, int *aj, 
        int *ak, int ijknum, int dim)
{    
    int ii = threadIdx.x + blockIdx.x * blockDim.x;         
    while (ii < ijknum)  {  // for each atom pair ij in the simulation box     
        int i = ai[ii];       // atom i
        int j = aj[ii];       // atom j
        int k = ak[ii];       // atom k
        // use atomicAdd on hip
        T ee = eijk[ii]/3.0;
//         e[i] += ee;
//         e[j] += ee;
//         e[k] += ee;
        atomicAdd(&e[i], ee);                
        atomicAdd(&e[j], ee);                
        atomicAdd(&e[k], ee);                
        ii += blockDim.x * gridDim.x;   
    }
}
template <typename T> void hipTripletDecomposition(T *e, T *eijk, int *ai, int *aj, 
        int *ak, int ijknum, int dim)
{                        
    int blockDim = 256;
    int gridDim = (ijknum + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;
      
    hipLaunchKernelGGL(hipKernelTripletDecomposition, gridDim, blockDim, 0, 0, e, eijk, ai, aj, ak, ijknum, dim);
}
template void hipTripletDecomposition(double*, double*, int*, int*, int*, int, int);
template void hipTripletDecomposition(float*, float*, int*, int*, int*, int, int);

template <typename T> __global__ void hipKernelCenterAtomTripletDecomposition(T *e, T *eijk, int *ilist, 
        int *anumsum, int inum, int dim)
{    
    int ii = threadIdx.x + blockIdx.x * blockDim.x;         
    while (ii < inum)  {  // for each atom i in the simulation box     
        int i = ilist[ii];       // atom i        
        int start = anumsum[ii];   
        int m = anumsum[ii+1]-start;        // number of neighbors around i             
        for (int l=0; l<m ; l++)    // loop over each atom around atom i            
            e[i] += eijk[start + l]/3.0;                                
        ii += blockDim.x * gridDim.x;   
    }
}
template <typename T> void hipCenterAtomTripletDecomposition(T *e, T *eijk, int *ilist, 
        int *anumsum, int inum, int dim)
{                        
    int blockDim = 256;
    int gridDim = (inum + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;
      
    hipLaunchKernelGGL(hipKernelCenterAtomTripletDecomposition, gridDim, blockDim, 0, 0, e, eijk, ilist, anumsum, inum, dim);
}
template void hipCenterAtomTripletDecomposition(double*, double*, int*, int*, int, int);
template void hipCenterAtomTripletDecomposition(float*, float*,  int*, int*, int, int);

template <typename T> __global__ void hipKernelNeighborAtomTripletDecomposition(T *e, T *eij, int *jlist, 
        int *bnumsum, int *index, int jnum, int dim)
{        
    int ii = threadIdx.x + blockIdx.x * blockDim.x;         
    while (ii < jnum)  {  // for each atom j in the simulation box     
        int j = jlist[ii];       // atom j        
        int start = bnumsum[ii];   
        int m = bnumsum[ii+1]-start;        // number of neighbors around i             
        for (int l=0; l<m ; l++) {   // loop over each atom around atom i
            int n = index[start + l];
            e[j] += eij[n]/3.0;                        
        }        
        ii += blockDim.x * gridDim.x;   
    }    
}
template <typename T> void hipNeighborAtomTripletDecomposition(T *e, T *eij, int *jlist, 
        int *bnumsum, int *index, int jnum, int dim)
{                        
    int blockDim = 256;
    int gridDim = (jnum + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;
      
    hipLaunchKernelGGL(hipKernelNeighborAtomTripletDecomposition, gridDim, blockDim, 0, 0, e, eij, jlist, bnumsum, index, jnum, dim);
}
template void hipNeighborAtomTripletDecomposition(double*, double*, int*, int*, int*, int, int);
template void hipNeighborAtomTripletDecomposition(float*, float*, int*, int*, int*, int, int);

template <typename T> __global__ void hipKernelQuadrupletDecomposition(T *e, T *eijkl, 
        int *ai, int *aj, int *ak, int *al, int ijklnum, int dim)
{    
    int ii = threadIdx.x + blockIdx.x * blockDim.x;         
    while (ii < ijklnum)  {  // for each atom quadruplet ijkl in the simulation box     
        int i = ai[ii];       // atom i
        int j = aj[ii];       // atom j
        int k = ak[ii];       // atom k
        int l = al[ii];       // atom l
        // use atomicAdd on hip
        T ee = eijkl[ii]/4.0;
//         e[i] += ee;
//         e[j] += ee;
//         e[k] += ee;
//         e[l] += ee;
        atomicAdd(&e[i], ee);                
        atomicAdd(&e[j], ee);                
        atomicAdd(&e[k], ee);                
        atomicAdd(&e[l], ee);                
        ii += blockDim.x * gridDim.x;   
    }
}
template <typename T> void hipQuadrupletDecomposition(T *e, T *eijkl, 
        int *ai, int *aj, int *ak, int *al, int ijklnum, int dim)
{                        
    int blockDim = 256;
    int gridDim = (ijklnum + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;
      
    hipLaunchKernelGGL(hipKernelQuadrupletDecomposition, gridDim, blockDim, 0, 0, e, eijkl, ai, aj, ak, al, ijklnum, dim);
}
template void hipQuadrupletDecomposition(double*, double*, int*, int*, int*, int*, int, int);
template void hipQuadrupletDecomposition(float*, float*, int*, int*, int*, int*, int, int);

template <typename T> __global__ void hipKernelCenterAtomQuadrupletDecomposition(T *e, T *eijkl, int *ilist, 
        int *anumsum, int inum, int dim)
{    
    int ii = threadIdx.x + blockIdx.x * blockDim.x;         
    while (ii < inum)  {  // for each atom i in the simulation box     
        int i = ilist[ii];       // atom i        
        int start = anumsum[ii];   
        int m = anumsum[ii+1]-start;        // number of neighbors around i             
        for (int l=0; l<m ; l++)    // loop over each atom around atom i            
            e[i] += eijkl[start + l]/4.0;                                
        ii += blockDim.x * gridDim.x;   
    }
}
template <typename T> void hipCenterAtomQuadrupletDecomposition(T *e, T *eijkl, int *ilist, 
        int *anumsum, int inum, int dim)
{                        
    int blockDim = 256;
    int gridDim = (inum + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;
      
    hipLaunchKernelGGL(hipKernelCenterAtomQuadrupletDecomposition, gridDim, blockDim, 0, 0, e, eijkl, ilist, anumsum, inum, dim);
}
template void hipCenterAtomQuadrupletDecomposition(double*, double*, int*, int*, int, int);
template void hipCenterAtomQuadrupletDecomposition(float*, float*, int*, int*, int, int);

template <typename T> __global__ void hipKernelNeighborAtomQuadrupletDecomposition(T *e, T *eij, int *jlist, 
        int *bnumsum, int *index, int jnum, int dim)
{        
    int ii = threadIdx.x + blockIdx.x * blockDim.x;         
    while (ii < jnum)  {  // for each atom j in the simulation box     
        int j = jlist[ii];       // atom j        
        int start = bnumsum[ii];   
        int m = bnumsum[ii+1]-start;        // number of neighbors around i                     
        for (int l=0; l<m ; l++) {   // loop over each atom around atom i
            int n = index[start + l];
            e[j] += eij[n]/4.0;                        
        }        
        ii += blockDim.x * gridDim.x;   
    }    
}
template <typename T> void hipNeighborAtomQuadrupletDecomposition(T *e, T *eij, int *jlist, 
        int *bnumsum, int *index, int jnum, int dim)
{                        
    int blockDim = 256;
    int gridDim = (jnum + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;
      
    hipLaunchKernelGGL(hipKernelNeighborAtomQuadrupletDecomposition, gridDim, blockDim, 0, 0, e, eij, jlist, bnumsum, index, jnum, dim);
}
template void hipNeighborAtomQuadrupletDecomposition(double*, double*, int*, int*, int*, int, int);
template void hipNeighborAtomQuadrupletDecomposition(float*, float*, int*, int*, int*, int, int);

template <typename T> __global__ void hipKernelSingleDecomposition(T *e, T *f, T *ei, T *fi, int *ai, int inum, int dim)
{    
    int ii = threadIdx.x + blockIdx.x * blockDim.x;         
    while (ii < inum)  {  // for each atom pair ij in the simulation box     
        int i = ai[ii];       // atom i        
        int l = dim*ii;                
        e[i] += ei[ii];
        i = dim*i;
        for (int d=0; d<dim; d++) 
            f[i+d] +=  fi[l+d];        
        ii += blockDim.x * gridDim.x;   
    }
}
template <typename T> void hipSingleDecomposition(T *e, T *f, T *ei, T *fi, int *ai, int inum, int dim)
{                        
    int blockDim = 256;
    int gridDim = (inum + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;
      
    hipLaunchKernelGGL(hipKernelSingleDecomposition, gridDim, blockDim, 0, 0, e, f, ei, fi, ai, inum, dim);
}
template void hipSingleDecomposition(double*, double*, double*, double*, int*, int, int);
template void hipSingleDecomposition(float*, float*, float*, float*, int*, int, int);

template <typename T> __global__ void hipKernelCenterAtomPairDecomposition(T *e, T *f, T *eij, 
        T *fij, int *ilist, int *anumsum, int inum, int dim)
{    
    int ii = threadIdx.x + blockIdx.x * blockDim.x;         
    while (ii < inum)  {  // for each atom i in the simulation box     
        int i = ilist[ii];       // atom i        
        int start = anumsum[ii];   
        int m = anumsum[ii+1]-start;        // number of neighbors around i             
        int idim = dim*i;        
        for (int l=0; l<m ; l++) {   // loop over each atom around atom i
            int n = start + l;
            int ndim = dim*n;
            e[i] += 0.5*eij[n];                        
            for (int d=0; d<dim; d++) 
                f[idim+d] +=  -fij[ndim+d];        
        }        
        ii += blockDim.x * gridDim.x;   
    }
}
template <typename T> void hipCenterAtomPairDecomposition(T *e, T *f, T *eij, 
        T *fij, int *ilist, int *anumsum, int inum, int dim)
{                        
    int blockDim = 256;
    int gridDim = (inum + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;
      
    hipLaunchKernelGGL(hipKernelCenterAtomPairDecomposition, gridDim, blockDim, 0, 0, e, f, eij, fij, ilist, anumsum, inum, dim);
}
template void hipCenterAtomPairDecomposition(double*, double*, double*, double*, int*, int*, int, int);
template void hipCenterAtomPairDecomposition(float*, float*, float*, float*, int*, int*, int, int);

template <typename T> __global__ void hipKernelNeighborAtomPairDecomposition(T *e, T *f, T *eij, T *fij, 
        int *jlist, int *bnumsum, int *index, int jnum, int dim)
{        
    int ii = threadIdx.x + blockIdx.x * blockDim.x;         
    while (ii < jnum)  {  // for each atom j in the simulation box     
        int j = jlist[ii];       // atom j        
        int start = bnumsum[ii];   
        int m = bnumsum[ii+1]-start;        // number of neighbors around i             
        int jdim = dim*j;        
        for (int l=0; l<m ; l++) {   // loop over each atom around atom i
            int n = index[start + l];
            int ndim = dim*n;
            e[j] += 0.5*eij[n];                        
            for (int d=0; d<dim; d++) 
                f[jdim+d] +=  fij[ndim+d];        
        }        
        ii += blockDim.x * gridDim.x;   
    }    
}
template <typename T> void hipNeighborAtomPairDecomposition(T *e, T *f, T *eij, T *fij, 
        int *jlist, int *bnumsum, int *index, int jnum, int dim)
{                        
    int blockDim = 256;
    int gridDim = (jnum + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;
      
    hipLaunchKernelGGL(hipKernelNeighborAtomPairDecomposition, gridDim, blockDim, 0, 0, e, f, eij, fij, jlist, bnumsum, index, jnum, dim);
}
template void hipNeighborAtomPairDecomposition(double*, double*, double*, double*, int*, int*, int*, int, int);
template void hipNeighborAtomPairDecomposition(float*, float*, float*, float*, int*, int*, int*, int, int);

template <typename T> __global__ void hipKernelFullNeighPairDecomposition(T *e, T *f, T *eij, T *fij, 
        int *ai, int ijnum, int dim)
{    
    int ii = threadIdx.x + blockIdx.x * blockDim.x;         
    while (ii < ijnum)  {  // for each atom pair ij in the simulation box     
        int i = ai[ii];       // atom i
        int l = dim*ii;
        // use atomicAdd on hip
        //e[i] += 0.5*eij[ii];
        atomicAdd(&e[i], 0.5*eij[ii]);                
        i = dim*i;
        for (int d=0; d<dim; d++) 
            atomicAdd(&f[i+d], -fij[l+d]);                
            //f[i+d] +=  -fij[l+d];        
        ii += blockDim.x * gridDim.x;   
    }
}
template <typename T> void hipFullNeighPairDecomposition(T *e, T *f, T *eij, T *fij, 
        int *ai, int ijnum, int dim)
{                        
    int blockDim = 256;
    int gridDim = (ijnum + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;
      
    hipLaunchKernelGGL(hipKernelFullNeighPairDecomposition, gridDim, blockDim, 0, 0, e, f, eij, fij, ai, ijnum, dim);
}
template void hipFullNeighPairDecomposition(double*, double*, double*, double*, int*, int, int);
template void hipFullNeighPairDecomposition(float*, float*, float*, float*, int*, int, int);

template <typename T> __global__ void hipKernelHalfNeighPairDecomposition(T *e, T *f, T *eij, T *fij, 
        int *ai, int *aj, int ijnum, int dim)
{    
    int ii = threadIdx.x + blockIdx.x * blockDim.x;         
    while (ii < ijnum)  {  // for each atom pair ij in the simulation box     
        int i = ai[ii];       // atom i
        int j = aj[ii];       // atom j
        int l = dim*ii;        
        // use atomicAdd on hip
        //e[i] += 0.5*eij[ii];
        //e[j] += 0.5*eij[ii];
        atomicAdd(&e[i], 0.5*eij[ii]);     
        atomicAdd(&e[j], 0.5*eij[ii]);     
        i = dim*i;
        j = dim*j;
        for (int d=0; d<dim; d++) {
            //f[i+d] += -fij[l+d]; 
            //f[j+d] +=  fij[l+d]; 
            atomicAdd(&f[i+d], -fij[l+d]);    
            atomicAdd(&f[j+d],  fij[l+d]);    
        }        
        ii += blockDim.x * gridDim.x;   
    }
}
template <typename T> void hipHalfNeighPairDecomposition(T *e, T *f, T *eij, T *fij, 
        int *ai, int *aj, int ijnum, int dim)
{                        
    int blockDim = 256;
    int gridDim = (ijnum + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;
      
    hipLaunchKernelGGL(hipKernelHalfNeighPairDecomposition, gridDim, blockDim, 0, 0, e, f, eij, fij, ai, aj, ijnum, dim);
}
template void hipHalfNeighPairDecomposition(double*, double*, double*, double*, int*, int*, int, int);
template void hipHalfNeighPairDecomposition(float*, float*, float*, float*, int*, int*, int, int);

template <typename T> __global__ void hipKernelTripletDecomposition(T *e, T *f, T *eijk, T *fij, 
        T *fik, int *ai, int *aj, int *ak, int ijknum, int dim)
{    
    int ii = threadIdx.x + blockIdx.x * blockDim.x;         
    while (ii < ijknum)  {  // for each atom pair ij in the simulation box     
        int i = ai[ii];       // atom i
        int j = aj[ii];       // atom j
        int k = ak[ii];       // atom k
        int n = dim*ii;        
        // use atomicAdd on hip
        T ee = eijk[ii]/3.0;
        //e[i] += ee;
        //e[j] += ee;
        //e[k] += ee;
        atomicAdd(&e[i], ee);             
        atomicAdd(&e[j], ee);             
        atomicAdd(&e[k], ee);             
        i = dim*i;
        j = dim*j;
        k = dim*k;
        for (int d=0; d<dim; d++) {
            //f[i+d] += -(fij[n+d] + fik[n+d]); 
            //f[j+d] +=  fij[n+d]; 
            //f[k+d] +=  fik[n+d]; 
            atomicAdd(&f[i+d], -(fij[n+d] + fik[n+d])); 
            atomicAdd(&f[j+d],  fij[n+d]);    
            atomicAdd(&f[k+d],  fik[n+d]);    
        }        
        ii += blockDim.x * gridDim.x;   
    }
}
template <typename T> void hipTripletDecomposition(T *e, T *f, T *eijk, T *fij, 
        T *fik, int *ai, int *aj, int *ak, int ijknum, int dim)
{                        
    int blockDim = 256;
    int gridDim = (ijknum + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;
      
    hipLaunchKernelGGL(hipKernelTripletDecomposition, gridDim, blockDim, 0, 0, e, f, eijk, fij, fik, ai, aj, ak, ijknum, dim);
}
template void hipTripletDecomposition(double*, double*, double*, double*, double*, int*, int*, int*, int, int);
template void hipTripletDecomposition(float*, float*, float*, float*, float*, int*, int*, int*, int, int);

template <typename T> __global__ void hipKernelCenterAtomTripletDecomposition(T *e, T *f, T *eijk, T *fij, 
        T *fik, int *ilist, int *anumsum, int inum, int dim)
{    
    int ii = threadIdx.x + blockIdx.x * blockDim.x;         
    while (ii < inum)  {  // for each atom i in the simulation box     
        int i = ilist[ii];       // atom i        
        int start = anumsum[ii];   
        int m = anumsum[ii+1]-start;        // number of neighbors around i             
        int idim = dim*i;        
        for (int l=0; l<m ; l++) {   // loop over each atom around atom i
            int n = start + l;
            int ndim = dim*n;
            e[i] += eijk[n]/3.0;                        
            for (int d=0; d<dim; d++) 
                f[idim+d] +=  -(fij[ndim+d]+fik[ndim+d]);        
        }        
        ii += blockDim.x * gridDim.x;   
    }
}
template <typename T> void hipCenterAtomTripletDecomposition(T *e, T *f, T *eijk, T *fij, 
        T *fik, int *ilist, int *anumsum, int inum, int dim)
{                        
    int blockDim = 256;
    int gridDim = (inum + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;
      
    hipLaunchKernelGGL(hipKernelCenterAtomTripletDecomposition, gridDim, blockDim, 0, 0, e, f, eijk, fij, fik, ilist, anumsum, inum, dim);
}
template void hipCenterAtomTripletDecomposition(double*, double*, double*, double*, double*, int*, int*, int, int);
template void hipCenterAtomTripletDecomposition(float*, float*, float*, float*, float*, int*, int*, int, int);

template <typename T> __global__ void hipKernelNeighborAtomTripletDecomposition(T *e, T *f, T *eij, T *fij, 
        int *jlist, int *bnumsum, int *index, int jnum, int dim)
{        
    int ii = threadIdx.x + blockIdx.x * blockDim.x;         
    while (ii < jnum)  {  // for each atom j in the simulation box     
        int j = jlist[ii];       // atom j        
        int start = bnumsum[ii];   
        int m = bnumsum[ii+1]-start;        // number of neighbors around i             
        int jdim = dim*j;        
        for (int l=0; l<m ; l++) {   // loop over each atom around atom i
            int n = index[start + l];
            int ndim = dim*n;
            e[j] += eij[n]/3.0;                        
            for (int d=0; d<dim; d++) 
                f[jdim+d] +=  fij[ndim+d];        
        }        
        ii += blockDim.x * gridDim.x;   
    }    
}
template <typename T> void hipNeighborAtomTripletDecomposition(T *e, T *f, T *eij, T *fij, 
        int *jlist, int *bnumsum, int *index, int jnum, int dim)
{                        
    int blockDim = 256;
    int gridDim = (jnum + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;
      
    hipLaunchKernelGGL(hipKernelNeighborAtomTripletDecomposition, gridDim, blockDim, 0, 0, e, f, eij, fij, jlist, bnumsum, index, jnum, dim);
}
template void hipNeighborAtomTripletDecomposition(double*, double*, double*, double*, int*, int*, int*, int, int);
template void hipNeighborAtomTripletDecomposition(float*, float*, float*, float*, int*, int*, int*, int, int);

template <typename T> __global__ void hipKernelQuadrupletDecomposition(T *e, T *f, T *eijkl, T *fij, T *fik, T *fil, 
        int *ai, int *aj, int *ak, int *al, int ijklnum, int dim)
{    
    int ii = threadIdx.x + blockIdx.x * blockDim.x;         
    while (ii < ijklnum)  {  // for each atom quadruplet ijkl in the simulation box     
        int i = ai[ii];       // atom i
        int j = aj[ii];       // atom j
        int k = ak[ii];       // atom k
        int l = al[ii];       // atom l
        int n = dim*ii;        
        // use atomicAdd on hip
        T ee = eijkl[ii]/4.0;
//         e[i] += ee;
//         e[j] += ee;
//         e[k] += ee;
//         e[l] += ee;
        atomicAdd(&e[i], ee);             
        atomicAdd(&e[j], ee);             
        atomicAdd(&e[k], ee);             
        atomicAdd(&e[l], ee);            
        i = dim*i;
        j = dim*j;
        k = dim*k;
        l = dim*l;
        for (int d=0; d<dim; d++) {
//             f[i+d] += -(fij[n+d] + fik[n+d] + fil[n+d]); 
//             f[j+d] +=  fij[n+d]; 
//             f[k+d] +=  fik[n+d]; 
//             f[l+d] +=  fil[n+d]; 
            atomicAdd(&f[i+d], -(fij[n+d] + fik[n+d] + fil[n+d])); 
            atomicAdd(&f[j+d],  fij[n+d]);    
            atomicAdd(&f[k+d],  fik[n+d]);   
            atomicAdd(&f[l+d],  fil[n+d]);   
        }        
        ii += blockDim.x * gridDim.x;   
    }
}
template <typename T> void hipQuadrupletDecomposition(T *e, T *f, T *eijkl, T *fij, T *fik, T *fil, 
        int *ai, int *aj, int *ak, int *al, int ijklnum, int dim)
{                        
    int blockDim = 256;
    int gridDim = (ijklnum + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;
      
    hipLaunchKernelGGL(hipKernelQuadrupletDecomposition, gridDim, blockDim, 0, 0, e, f, eijkl, fij, fik, fil, ai, aj, ak, al, ijklnum, dim);
}
template void hipQuadrupletDecomposition(double*, double*, double*, double*, double*, double*, int*, int*, int*, int*, int, int);
template void hipQuadrupletDecomposition(float*, float*, float*, float*, float*, float*, int*, int*, int*, int*, int, int);

template <typename T> __global__ void hipKernelCenterAtomQuadrupletDecomposition(T *e, T *f, T *eijkl, T *fij, 
        T *fik, T *fil, int *ilist, int *anumsum, int inum, int dim)
{    
    int ii = threadIdx.x + blockIdx.x * blockDim.x;         
    while (ii < inum)  {  // for each atom i in the simulation box     
        int i = ilist[ii];       // atom i        
        int start = anumsum[ii];   
        int m = anumsum[ii+1]-start;        // number of neighbors around i             
        int idim = dim*i;        
        for (int l=0; l<m ; l++) {   // loop over each atom around atom i
            int n = start + l;
            int ndim = dim*n;
            e[i] += eijkl[n]/4.0;                        
            for (int d=0; d<dim; d++) 
                f[idim+d] +=  -(fij[ndim+d]+fik[ndim+d]+fil[ndim+d]);        
        }        
        ii += blockDim.x * gridDim.x;   
    }
}
template <typename T> void hipCenterAtomQuadrupletDecomposition(T *e, T *f, T *eijkl, T *fij, 
        T *fik, T *fil, int *ilist, int *anumsum, int inum, int dim)
{                        
    int blockDim = 256;
    int gridDim = (inum + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;
      
    hipLaunchKernelGGL(hipKernelCenterAtomQuadrupletDecomposition, gridDim, blockDim, 0, 0, e, f, eijkl, fij, fik, fil, ilist, anumsum, inum, dim);
}
template void hipCenterAtomQuadrupletDecomposition(double*, double*, double*, double*, double*, double*, int*, int*, int, int);
template void hipCenterAtomQuadrupletDecomposition(float*, float*, float*, float*, float*, float*, int*, int*, int, int);

template <typename T> __global__ void hipKernelNeighborAtomQuadrupletDecomposition(T *e, T *f, T *eij, T *fij, 
        int *jlist, int *bnumsum, int *index, int jnum, int dim)
{        
    int ii = threadIdx.x + blockIdx.x * blockDim.x;         
    while (ii < jnum)  {  // for each atom j in the simulation box     
        int j = jlist[ii];       // atom j        
        int start = bnumsum[ii];   
        int m = bnumsum[ii+1]-start;        // number of neighbors around i             
        int jdim = dim*j;        
        for (int l=0; l<m ; l++) {   // loop over each atom around atom i
            int n = index[start + l];
            int ndim = dim*n;
            e[j] += eij[n]/4.0;                        
            for (int d=0; d<dim; d++) 
                f[jdim+d] +=  fij[ndim+d];        
        }        
        ii += blockDim.x * gridDim.x;   
    }    
}
template <typename T> void hipNeighborAtomQuadrupletDecomposition(T *e, T *f, T *eij, T *fij, 
        int *jlist, int *bnumsum, int *index, int jnum, int dim)
{                        
    int blockDim = 256;
    int gridDim = (jnum + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;
      
    hipLaunchKernelGGL(hipKernelNeighborAtomQuadrupletDecomposition, gridDim, blockDim, 0, 0, e, f, eij, fij, jlist, bnumsum, index, jnum, dim);
}
template void hipNeighborAtomQuadrupletDecomposition(double*, double*, double*, double*, int*, int*, int*, int, int);
template void hipNeighborAtomQuadrupletDecomposition(float*, float*, float*, float*, int*, int*, int*, int, int);


template <typename T> __global__ void hipKernelTripletForceDecomposition(T *f3ik, T *eij, T *d3ij, 
        int *tripletnum, int *tripletnumsum, int npairs, int dim)
{        
    int i = threadIdx.x + blockIdx.x * blockDim.x;   
    while (i < npairs)  {  
        int m = tripletnum[i];
        int s = tripletnumsum[i];
        for (int j=0; j<m; j++) 
            for (int d=0; d<dim; d++)
                f3ik[d+3*(s+j)] = eij[i]*d3ij[i]*f3ik[d+3*(s+j)];    

        i += blockDim.x * gridDim.x;  
    }                            
}
template <typename T> void hipTripletForceDecomposition(T *f3ik, T *eij, T *d3ij, 
        int *tripletnum, int *tripletnumsum, int npairs, int dim)
{                        
    int blockDim = 256;
    int gridDim = (npairs + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;
      
    hipLaunchKernelGGL(hipKernelTripletForceDecomposition, gridDim, blockDim, 0, 0, f3ik, eij, d3ij, tripletnum, tripletnumsum, npairs, dim);
}
template void hipTripletForceDecomposition(double*, double*, double*, int*, int*, int, int);
template void hipTripletForceDecomposition(float*, float*, float*, int*, int*, int, int);


template <typename T> __global__ void hipKernelHalfForceDecomposition2D(T *f, T *fij, int *ai, int *aj, int ijnum)
{    
    int ii = threadIdx.x + blockIdx.x * blockDim.x;     
    while (ii < ijnum)  {  // for each atom pair ij in the simulation box     
        int i = 2*ai[ii];       // atom i
        int j = 2*aj[ii];       // atom j
        int l = 2*ii;
        // use atomicAdd on hip
//         f[i+0] +=  -fij[l+0];
//         f[i+1] +=  -fij[l+1];
//         f[j+0] +=  fij[l+0];
//         f[j+1] +=  fij[l+1];
        atomicAdd(&f[0+i], -fij[0+l]);
        atomicAdd(&f[1+i], -fij[1+l]);
        atomicAdd(&f[0+j],  fij[0+l]);
        atomicAdd(&f[1+j],  fij[1+l]);
        ii += blockDim.x * gridDim.x;   
    }
}
template <typename T> void hipHalfForceDecomposition2D(T *f, T *fij, int *ai, int *aj, int ijnum)
{                        
    int blockDim = 256;
    int gridDim = (ijnum + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;
      
    hipLaunchKernelGGL(hipKernelHalfForceDecomposition2D, gridDim, blockDim, 0, 0, f, fij, ai, aj, ijnum);
}
template void hipHalfForceDecomposition2D(double*, double*, int*, int*, int);
template void hipHalfForceDecomposition2D(float*, float*, int*, int*, int);

template <typename T> __global__ void hipKernelFullForceDecomposition2D(T *f, T *fij, int *ai, int *aj, int ijnum)
{    
    int ii = threadIdx.x + blockIdx.x * blockDim.x;     
    while (ii < ijnum)  {  // for each atom pair ij in the simulation box     
        int i = 2*ai[ii];       // atom i
        int l = 2*ii;
        // use atomicAdd on hip
//         f[i+0] +=  -fij[l+0];
//         f[i+1] +=  -fij[l+1];
        atomicAdd(&f[0+i], -fij[0+l]);
        atomicAdd(&f[1+i], -fij[1+l]);
        ii += blockDim.x * gridDim.x;   
    }
}
template <typename T> void hipFullForceDecomposition2D(T *f, T *fij, int *ai, int *aj, int ijnum)
{                        
    int blockDim = 256;
    int gridDim = (ijnum + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;
      
    hipLaunchKernelGGL(hipKernelFullForceDecomposition2D, gridDim, blockDim, 0, 0, f, fij, ai, aj, ijnum);
}
template void hipFullForceDecomposition2D(double*, double*, int*, int*, int);
template void hipFullForceDecomposition2D(float*, float*, int*, int*, int);

template <typename T> __global__ void hipKernelIAtomDecomposition2D(T *f, T *fij, int *ilist, int *anumsum, int inum)
{    
    int ii = threadIdx.x + blockIdx.x * blockDim.x;     
    while (ii < inum)  {  // for each atom i in the simulation box     
        int i = 2*ilist[ii];       // atom i        
        int start = anumsum[ii];   
        int m = anumsum[ii+1]-start;        // number of neighbors around i             
        T f1 = 0.0;
        T f2 = 0.0;
        for (int l=0; l<m ; l++) {   // loop over each atom around atom i
            int k = 2*(start + l);                     
            f1 +=  -fij[k+0];
            f2 +=  -fij[k+1];
        }
        f[i+0] = f1;
        f[i+1] = f2;
        ii += blockDim.x * gridDim.x;   
    }
}
template <typename T> void hipIAtomDecomposition2D(T *f, T *fij, int *ilist, int *anumsum, int inum)
{                        
    int blockDim = 256;
    int gridDim = (inum + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;
      
    hipLaunchKernelGGL(hipKernelIAtomDecomposition2D, gridDim, blockDim, 0, 0, f, fij, ilist, anumsum, inum);
}
template void hipIAtomDecomposition2D(double*, double*, int*, int*, int);
template void hipIAtomDecomposition2D(float*, float*, int*, int*, int);

template <typename T> __global__ void hipKernelJAtomDecomposition2D(T *f, T *fij, 
        int *jlist, int *bnumsum, int *index, int jnum)
{        
    int ii = threadIdx.x + blockIdx.x * blockDim.x;     
    while (ii < jnum)  {  // for each atom j in the simulation box     
        int j = 2*jlist[ii];       // atom j
        T f1 = 0.0;
        T f2 = 0.0;                
        // need to determine bnumsum and index 
        int start = bnumsum[ii];   
        int m = bnumsum[ii+1]-start;   // number of neighbors i around j  (j < i)           
        for (int l=0; l<m ; l++) { // loop over each atom around i atom j -> pair ij 
            int k = 2*index[start + l];                     
            f1 +=  fij[k+0];
            f2 +=  fij[k+1];
        }
        f[j+0] = f1;
        f[j+1] = f2;
        ii += blockDim.x * gridDim.x;   
    }
}
template <typename T> void hipJAtomDecomposition2D(T *f, T *fij, 
        int *jlist, int *bnumsum, int *index, int jnum)
{                        
    int blockDim = 256;
    int gridDim = (jnum + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;
      
    hipLaunchKernelGGL(hipKernelJAtomDecomposition2D, gridDim, blockDim, 0, 0, f, fij, jlist, bnumsum, index, jnum);
}
template void hipJAtomDecomposition2D(double*, double*, int*, int*, int*, int);
template void hipJAtomDecomposition2D(float*, float*, int*, int*, int*, int);

template <typename T> __global__ void hipKernelForceDecompositionTriplet2D(T *f, T *fij, T *fik, 
        int *ai, int *aj, int *ak, int ijknum)
{    
    int ii = threadIdx.x + blockIdx.x * blockDim.x;     
    while (ii < ijknum)  {  // for each atom pair ij in the simulation box     
        int i = 2*ai[ii];       // atom i
        int j = 2*aj[ii];       // atom j
        int k = 2*ak[ii];       // atom k
        int l = 2*ii;
        // use atomicAdd on hip
//         f[i+0] +=  -(fij[l+0]+fik[l+0]);
//         f[i+1] +=  -(fij[l+1]+fik[l+1]);
//         f[j+0] +=  fij[l+0];
//         f[j+1] +=  fij[l+1];
//         f[k+0] +=  fik[l+0];
//         f[k+1] +=  fik[l+1];
        atomicAdd(&f[0+i], -(fij[0+l]+fik[l+0]));
        atomicAdd(&f[1+i], -(fij[1+l]+fik[l+1]));
        atomicAdd(&f[0+j], fij[0+l]);
        atomicAdd(&f[1+j], fij[1+l]);
        atomicAdd(&f[0+k], fik[0+l]);
        atomicAdd(&f[1+k], fik[1+l]);
        ii += blockDim.x * gridDim.x;   
    }
}
template <typename T> void hipForceDecompositionTriplet2D(T *f, T *fij, T *fik, 
        int *ai, int *aj, int *ak, int ijknum)
{                        
    int blockDim = 256;
    int gridDim = (ijknum + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;
      
    hipLaunchKernelGGL(hipKernelForceDecompositionTriplet2D, gridDim, blockDim, 0, 0, f, fij, fik, ai, aj, ak, ijknum);
}
template void hipForceDecompositionTriplet2D(double*, double*, double*, int*, int*, int*, int);
template void hipForceDecompositionTriplet2D(float*, float*, float*, int*, int*, int*, int);

template <typename T> __global__ void hipKernelAtomDecompositionTriplet2D(T *f, T *fij, T *fik, 
        int *ilist, int *anumsum, int inum)
{    
    int ii = threadIdx.x + blockIdx.x * blockDim.x;     
    while (ii < inum)  {  // for each atom i in the simulation box     
        int i = 2*ilist[ii];       // atom i        
        int start = anumsum[ii];   
        int m = anumsum[ii+1]-start;        // number of neighbors around i             
        T f1 = 0.0;
        T f2 = 0.0;
        for (int l=0; l<m ; l++) {   // loop over each atom pair jk around atom i
            int s = start + l;
            int n = 2*s;                     
            f1 +=  -(fij[n+0] + fik[n+0]);
            f2 +=  -(fij[n+1] + fik[n+1]);
//             // use atomicAdd on hip
//             int j = 2*aj[s];
//             int k = 2*ak[s];            
//             f[j+0] +=  fij[n+0];
//             f[j+1] +=  fij[n+1];
//             f[k+0] +=  fik[n+0];
//             f[k+1] +=  fik[n+1];            
        }
        f[i+0] = f1;
        f[i+1] = f2;                      
        ii += blockDim.x * gridDim.x;   
    }
}
template <typename T> void hipAtomDecompositionTriplet2D(T *f, T *fij, T *fik, 
        int *ilist, int *anumsum, int inum)
{                        
    int blockDim = 256;
    int gridDim = (inum + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;
      
    hipLaunchKernelGGL(hipKernelAtomDecompositionTriplet2D, gridDim, blockDim, 0, 0, f, fij, fik, ilist, anumsum, inum);
}
template void hipAtomDecompositionTriplet2D(double*, double*, double*, int*, int*, int);
template void hipAtomDecompositionTriplet2D(float*, float*, float*, int*, int*, int);

template <typename T> __global__ void hipKernelForceDecompositionQuadruplet2D(T *f, T *fij, T *fik,  T *fil, 
        int *ai, int *aj, int *ak, int *al, int ijklnum)
{    
    for (int ii=0; ii<ijklnum; ii++) {  // for each atom pair ij in the simulation box     
        int i = 2*ai[ii];       // atom i
        int j = 2*aj[ii];       // atom j
        int k = 2*ak[ii];       // atom k
        int l = 2*al[ii];       // atom l
        int n = 2*ii;
        // use atomicAdd on hip
//         f[i+0] +=  -(fij[n+0]+fik[n+0]+fil[n+0]);
//         f[i+1] +=  -(fij[n+1]+fik[n+1]+fil[n+1]);
//         f[j+0] +=  fij[n+0];
//         f[j+1] +=  fij[n+1];
//         f[k+0] +=  fik[n+0];
//         f[k+1] +=  fik[n+1];
//         f[l+0] +=  fil[n+0];
//         f[l+1] +=  fil[n+1];
        atomicAdd(&f[0+i], -(fij[0+n]+fik[n+0]+fil[n+0]));
        atomicAdd(&f[1+i], -(fij[1+n]+fik[n+1]+fil[n+1]));
        atomicAdd(&f[0+j], fij[0+n]);
        atomicAdd(&f[1+j], fij[1+n]);
        atomicAdd(&f[0+k], fik[0+n]);
        atomicAdd(&f[1+k], fik[1+n]);
        atomicAdd(&f[0+l], fil[0+n]);
        atomicAdd(&f[1+l], fil[1+n]);
        ii += blockDim.x * gridDim.x;   
    }
}
template <typename T> void hipForceDecompositionQuadruplet2D(T *f, T *fij, T *fik,  T *fil, 
        int *ai, int *aj, int *ak, int *al, int ijklnum)
{                        
    int blockDim = 256;
    int gridDim = (ijklnum + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;
      
    hipLaunchKernelGGL(hipKernelForceDecompositionQuadruplet2D, gridDim, blockDim, 0, 0, f, fij, fik, fil, ai, aj, ak, al, ijklnum);
}
template void hipForceDecompositionQuadruplet2D(double*, double*, double*, double*, int*, int*,  int*, int*, int);
template void hipForceDecompositionQuadruplet2D(float*, float*, float*, float*, int*, int*, int*, int*, int);

template <typename T> __global__ void hipKernelAtomDecompositionQuadruplet2D(T *f, T *fij, T *fik, T *fil,
        int *ilist, int *anumsum, int inum)
{    
    int ii = threadIdx.x + blockIdx.x * blockDim.x;     
    while (ii < inum)  {  // for each atom i in the simulation box     
        int i = 2*ilist[ii];       // atom i        
        int start = anumsum[ii];   
        int m = anumsum[ii+1]-start;        // number of neighbors around i             
        T f1 = 0.0;
        T f2 = 0.0;
        for (int l=0; l<m ; l++) {   // loop over each atom pair jk around atom i
            int s = start + l;
            int n = 2*s;                     
            f1 +=  -(fij[n+0] + fik[n+0] + fil[n+0]);
            f2 +=  -(fij[n+1] + fik[n+1] + fil[n+1]);
        }
        f[i+0] = f1;
        f[i+1] = f2;                        
        ii += blockDim.x * gridDim.x;   
    }
}
template <typename T> void hipAtomDecompositionQuadruplet2D(T *f, T *fij, T *fik, T *fil,
        int *ilist, int *anumsum, int inum)
{                        
    int blockDim = 256;
    int gridDim = (inum + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;
      
    hipLaunchKernelGGL(hipKernelAtomDecompositionQuadruplet2D, gridDim, blockDim, 0, 0, f, fij, fik, fil, ilist, anumsum, inum);
}
template void hipAtomDecompositionQuadruplet2D(double*, double*, double*, double*, int*, int*, int);
template void hipAtomDecompositionQuadruplet2D(float*, float*, float*, float*, int*, int*, int);

template <typename T> __global__ void hipKernelHalfForceDecomposition3D(T *f, T *fij, int *ai, int *aj, int ijnum)
{    
    int ii = threadIdx.x + blockIdx.x * blockDim.x;     
    while (ii < ijnum)  {  // for each atom pair ij in the simulation box     
        int i = 3*ai[ii];       // atom i
        int j = 3*aj[ii];       // atom j
        int l = 3*ii;
        // use atomicAdd on hip
//         f[i+0] +=  -fij[l+0];
//         f[i+1] +=  -fij[l+1];
//         f[i+2] +=  -fij[l+2];                 
//         f[j+0] +=  fij[l+0];
//         f[j+1] +=  fij[l+1];
//         f[j+2] +=  fij[l+2];              
        atomicAdd(&f[0+i], -fij[0+l]);
        atomicAdd(&f[1+i], -fij[1+l]);
        atomicAdd(&f[2+i], -fij[2+l]);                
        atomicAdd(&f[0+j],  fij[0+l]);
        atomicAdd(&f[1+j],  fij[1+l]);
        atomicAdd(&f[2+j],  fij[2+l]);   
        ii += blockDim.x * gridDim.x;   
    }
}
template <typename T> void hipHalfForceDecomposition3D(T *f, T *fij, int *ai, int *aj, int ijnum)
{                        
    int blockDim = 256;
    int gridDim = (ijnum + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;
      
    hipLaunchKernelGGL(hipKernelHalfForceDecomposition3D, gridDim, blockDim, 0, 0, f, fij, ai, aj, ijnum);
}
template void hipHalfForceDecomposition3D(double*, double*, int*, int*, int);
template void hipHalfForceDecomposition3D(float*, float*, int*, int*, int);

template <typename T> __global__ void hipKernelFullForceDecomposition3D(T *f, T *fij, int *ai, int *aj, int ijnum)
{    
    int ii = threadIdx.x + blockIdx.x * blockDim.x;     
    while (ii < ijnum)  {  // for each atom pair ij in the simulation box     
        int i = 3*ai[ii];             // atom i
        int l = 3*ii;
        // use atomicAdd on hip
//         f[i+0] +=  -fij[l+0];
//         f[i+1] +=  -fij[l+1];
//         f[i+2] +=  -fij[l+2];                  
        atomicAdd(&f[0+i], -fij[0+l]);
        atomicAdd(&f[1+i], -fij[1+l]);
        atomicAdd(&f[2+i], -fij[2+l]);                
        ii += blockDim.x * gridDim.x;   
    }
}
template <typename T> void hipFullForceDecomposition3D(T *f, T *fij, int *ai, int *aj, int ijnum)
{                        
    int blockDim = 256;
    int gridDim = (ijnum + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;
      
    hipLaunchKernelGGL(hipKernelFullForceDecomposition3D, gridDim, blockDim, 0, 0, f, fij, ai, aj, ijnum);
}
template void hipFullForceDecomposition3D(double*, double*, int*, int*, int);
template void hipFullForceDecomposition3D(float*, float*, int*, int*, int);

template <typename T> __global__ void hipKernelIAtomDecomposition3D(T *f, T *fij, int *ilist, int *anumsum, int inum)
{    
    int ii = threadIdx.x + blockIdx.x * blockDim.x;     
    while (ii < inum)  {  // for each atom i in the simulation box     
        int i = 3*ilist[ii];       // atom i        
        T f1 = 0.0;
        T f2 = 0.0;        
        T f3 = 0.0;        
        int start = anumsum[ii];   
        int m = anumsum[ii+1]-start;        // number of neighbors around i             
        for (int l=0; l<m ; l++) {   // loop over each atom around atom i
            int k = 3*(start + l);                     
            f1 +=  -fij[k+0];
            f2 +=  -fij[k+1];
            f3 +=  -fij[k+2];             
        }
        f[i+0] = f1;
        f[i+1] = f2;
        f[i+2] = f3;
        ii += blockDim.x * gridDim.x;   
    }
}
template <typename T> void hipIAtomDecomposition3D(T *f, T *fij, int *ilist, int *anumsum, int inum)
{                        
    int blockDim = 256;
    int gridDim = (inum + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;
      
    hipLaunchKernelGGL(hipKernelIAtomDecomposition3D, gridDim, blockDim, 0, 0, f, fij, ilist, anumsum, inum);
}
template void hipIAtomDecomposition3D(double*, double*, int*, int*, int);
template void hipIAtomDecomposition3D(float*, float*, int*, int*, int);

template <typename T> __global__ void hipKernelJAtomDecomposition3D(T *f, T *fij, int *jlist, 
        int *bnumsum, int *index, int jnum)
{        
    int ii = threadIdx.x + blockIdx.x * blockDim.x;     
    while (ii < jnum)  {  // for each atom j in the simulation box     
        int j = 3*jlist[ii];       // atom j
        T f1 = 0.0;
        T f2 = 0.0;                
        T f3 = 0.0;                
        // need to determine bnumsum and index 
        int start = bnumsum[ii];   
        int m = bnumsum[ii+1]-start;   // number of neighbors i around j  (j < i)           
        for (int l=0; l<m ; l++) { // loop over each atom around i atom j -> pair ij 
            int k = 3*index[start + l];                     
            f1 +=  fij[k+0];
            f2 +=  fij[k+1];
            f3 +=  fij[k+2];
        }
        f[j+0] = f1;
        f[j+1] = f2;
        f[j+2] = f3;
        ii += blockDim.x * gridDim.x;   
    }
}
template <typename T> void hipJAtomDecomposition3D(T *f, T *fij, 
        int *jlist, int *bnumsum, int *index, int jnum)
{                        
    int blockDim = 256;
    int gridDim = (jnum + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;
      
    hipLaunchKernelGGL(hipKernelJAtomDecomposition3D, gridDim, blockDim, 0, 0, f, fij, jlist, bnumsum, index, jnum);
}
template void hipJAtomDecomposition3D(double*, double*, int*, int*, int*, int);
template void hipJAtomDecomposition3D(float*, float*, int*, int*, int*, int);

template <typename T> __global__ void hipKernelForceDecompositionTriplet3D(T *f, T *fij, T *fik, 
        int *ai, int *aj, int *ak, int ijknum)
{    
    int ii = threadIdx.x + blockIdx.x * blockDim.x;     
    while (ii < ijknum)  {  // for each atom pair ij in the simulation box     
        int i = 3*ai[ii];       // atom i
        int j = 3*aj[ii];       // atom j
        int k = 3*ak[ii];       // atom k
        int l = 3*ii;
        // use atomicAdd on hip
//         f[i+0] +=  -(fij[l+0]+fik[l+0]);
//         f[i+1] +=  -(fij[l+1]+fik[l+1]);
//         f[i+2] +=  -(fij[l+2]+fik[l+2]);
//         f[j+0] +=  fij[l+0];
//         f[j+1] +=  fij[l+1];
//         f[j+2] +=  fij[l+2];
//         f[k+0] +=  fik[l+0];
//         f[k+1] +=  fik[l+1];
//         f[k+2] +=  fik[l+2];
        atomicAdd(&f[0+i], -(fij[0+l]+fik[l+0]));
        atomicAdd(&f[1+i], -(fij[1+l]+fik[l+1]));
        atomicAdd(&f[2+i], -(fij[2+l]+fik[l+2]));
        atomicAdd(&f[0+j], fij[0+l]);
        atomicAdd(&f[1+j], fij[1+l]);
        atomicAdd(&f[2+j], fij[2+l]);
        atomicAdd(&f[0+k], fik[0+l]);
        atomicAdd(&f[1+k], fik[1+l]);
        atomicAdd(&f[2+k], fik[2+l]);
        ii += blockDim.x * gridDim.x;   
    }
}
template <typename T> void hipForceDecompositionTriplet3D(T *f, T *fij, T *fik, 
        int *ai, int *aj, int *ak, int ijknum)
{                        
    int blockDim = 256;
    int gridDim = (ijknum + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;
      
    hipLaunchKernelGGL(hipKernelForceDecompositionTriplet3D, gridDim, blockDim, 0, 0, f, fij, fik, ai, aj, ak, ijknum);
}
template void hipForceDecompositionTriplet3D(double*, double*, double*, int*, int*, int*, int);
template void hipForceDecompositionTriplet3D(float*, float*, float*, int*, int*, int*, int);

template <typename T> __global__ void hipKernelAtomDecompositionTriplet3D(T *f, T *fij, T *fik, 
        int *ilist, int *anumsum, int inum)
{    
    int ii = threadIdx.x + blockIdx.x * blockDim.x;     
    while (ii < inum)  {  // for each atom i in the simulation box     
        int i = 3*ilist[ii];       // atom i        
        int start = anumsum[ii];   
        int m = anumsum[ii+1]-start;        // number of neighbors around i             
        T f1 = 0.0;
        T f2 = 0.0;
        T f3 = 0.0;
        for (int l=0; l<m ; l++) {   // loop over each atom pair jk around atom i
            int s = start + l;
            int n = 3*s;                     
            f1 +=  -(fij[n+0] + fik[n+0]);
            f2 +=  -(fij[n+1] + fik[n+1]);
            f3 +=  -(fij[n+2] + fik[n+2]);
        }
        f[i+0] = f1;
        f[i+1] = f2;                        
        f[i+2] = f3;                        
        ii += blockDim.x * gridDim.x;   
    }
}
template <typename T> void hipAtomDecompositionTriplet3D(T *f, T *fij, T *fik, 
        int *ilist, int *anumsum, int inum)
{                        
    int blockDim = 256;
    int gridDim = (inum + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;
      
    hipLaunchKernelGGL(hipKernelAtomDecompositionTriplet3D, gridDim, blockDim, 0, 0, f, fij, fik, ilist, anumsum, inum);
}
template void hipAtomDecompositionTriplet3D(double*, double*, double*, int*, int*, int);
template void hipAtomDecompositionTriplet3D(float*, float*, float*, int*, int*, int);

template <typename T> __global__ void hipKernelForceDecompositionQuadruplet3D(T *f, T *fij, T *fik,  
        T *fil, int *ai, int *aj, int *ak, int *al, int ijklnum)
{    
    int ii = threadIdx.x + blockIdx.x * blockDim.x;     
    while (ii < ijklnum)  {  // for each atom pair ij in the simulation box     
        int i = 3*ai[ii];       // atom i
        int j = 3*aj[ii];       // atom j
        int k = 3*ak[ii];       // atom k
        int l = 3*al[ii];       // atom l
        int n = 3*ii;
        // use atomicAdd on hip
//         f[i+0] +=  -(fij[n+0]+fik[n+0]+fil[n+0]);
//         f[i+1] +=  -(fij[n+1]+fik[n+1]+fil[n+1]);
//         f[i+2] +=  -(fij[n+2]+fik[n+2]+fil[n+2]);
//         f[j+0] +=  fij[n+0];
//         f[j+1] +=  fij[n+1];
//         f[j+2] +=  fij[n+2];
//         f[k+0] +=  fik[n+0];
//         f[k+1] +=  fik[n+1];
//         f[k+2] +=  fik[n+2];
//         f[l+0] +=  fil[n+0];
//         f[l+1] +=  fil[n+1];
//         f[l+2] +=  fil[n+2];
        atomicAdd(&f[0+i], -(fij[0+n]+fik[n+0]+fil[n+0]));
        atomicAdd(&f[1+i], -(fij[1+n]+fik[n+1]+fil[n+1]));
        atomicAdd(&f[2+i], -(fij[2+n]+fik[n+2]+fil[n+2]));
        atomicAdd(&f[0+j], fij[0+n]);
        atomicAdd(&f[1+j], fij[1+n]);
        atomicAdd(&f[2+j], fij[2+n]);
        atomicAdd(&f[0+k], fik[0+n]);
        atomicAdd(&f[1+k], fik[1+n]);
        atomicAdd(&f[2+k], fik[2+n]);
        atomicAdd(&f[0+l], fil[0+n]);
        atomicAdd(&f[1+l], fil[1+n]);
        atomicAdd(&f[2+l], fil[2+n]);
        ii += blockDim.x * gridDim.x;   
    }
}
template <typename T> void hipForceDecompositionQuadruplet3D(T *f, T *fij, T *fik,  T *fil, 
        int *ai, int *aj, int *ak, int *al, int ijklnum)
{                        
    int blockDim = 256;
    int gridDim = (ijklnum + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;
      
    hipLaunchKernelGGL(hipKernelForceDecompositionQuadruplet3D, gridDim, blockDim, 0, 0, f, fij, fik, fil, ai, aj, ak, al, ijklnum);
}
template void hipForceDecompositionQuadruplet3D(double*, double*, double*, double*, int*, int*,  int*, int*, int);
template void hipForceDecompositionQuadruplet3D(float*, float*, float*, float*, int*, int*, int*, int*, int);

template <typename T> __global__ void hipKernelAtomDecompositionQuadruplet3D(T *f, T *fij, T *fik, T *fil, 
        int *ilist, int *anumsum, int inum)
{    
    int ii = threadIdx.x + blockIdx.x * blockDim.x;     
    while (ii < inum)  {  // for each atom i in the simulation box     
        int i = 3*ilist[ii];       // atom i        
        int start = anumsum[ii];   
        int m = anumsum[ii+1]-start;        // number of neighbors around i             
        T f1 = 0.0;
        T f2 = 0.0;
        T f3 = 0.0;
        for (int l=0; l<m ; l++) {   // loop over each atom pair jk around atom i
            int s = start + l;
            int n = 3*s;                     
            f1 +=  -(fij[n+0] + fik[n+0] + fil[n+0]);
            f2 +=  -(fij[n+1] + fik[n+1] + fil[n+1]);
            f3 +=  -(fij[n+2] + fik[n+2] + fil[n+2]);
        }
        f[i+0] = f1;
        f[i+1] = f2;                        
        f[i+2] = f3;                        
        ii += blockDim.x * gridDim.x;   
    }
}
template <typename T> void hipAtomDecompositionQuadruplet3D(T *f, T *fij, T *fik, T *fil,
        int *ilist, int *anumsum, int inum)
{                        
    int blockDim = 256;
    int gridDim = (inum + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;
      
    hipLaunchKernelGGL(hipKernelAtomDecompositionQuadruplet3D, gridDim, blockDim, 0, 0, f, fij, fik, fil, ilist, anumsum, inum);
}
template void hipAtomDecompositionQuadruplet3D(double*, double*, double*, double*, int*, int*, int);
template void hipAtomDecompositionQuadruplet3D(float*, float*, float*, float*, int*, int*, int);

#endif


