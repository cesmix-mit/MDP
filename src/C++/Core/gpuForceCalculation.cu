// assign threads to same-colored cells
// each thread gets its cell number
// each thread loops over atoms in its cell
//     compute xij, ti, tj
//     compute forces
//     compute basis function derivatives

// assign threads to each atom
// each thread loops over atoms 
//     compute basis function derivatives
//     assemble forces using atomicAdd to avoid 

// analyze aj

#ifndef __GPUFORCECALCULATION
#define __GPUFORCECALCULATION

template <typename T> __global__ void gpuKernelSingleDecomposition(T *e, T *ei, int *ai, int inum, int dim)
{    
    int ii = threadIdx.x + blockIdx.x * blockDim.x;         
    while (ii < inum) { 
        int i = ai[ii];       // atom i        
        int l = dim*ii;                
        e[i] += ei[ii];
        ii += blockDim.x * gridDim.x;   
    }
}
template <typename T> void gpuSingleDecomposition(T *e, T *ei, int *ai, int inum, int dim)
{                        
    int blockDim = 256;
    int gridDim = (inum + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;
      
    gpuKernelSingleDecomposition<<<gridDim, blockDim>>>(e, ei, ai, inum, dim);
}
template void gpuSingleDecomposition(double*, double*, int*, int, int);
template void gpuSingleDecomposition(float*, float*, int*, int, int);

template <typename T> __global__ void gpuKernelFullNeighPairDecomposition(T *e, T *eij, int *ai, int ijnum, int dim)
{    
    int ii = threadIdx.x + blockIdx.x * blockDim.x;         
    while (ii < ijnum)  {  // for each atom pair ij in the simulation box     
        int i = ai[ii];       // atom i
        int l = dim*ii;
        // use atomicAdd on gpu
        //e[i] += 0.5*eij[ii];
        atomicAdd(&e[i], 0.5*eij[ii]);
        ii += blockDim.x * gridDim.x;   
    }
}
template <typename T> void gpuFullNeighPairDecomposition(T *e, T *eij, int *ai, int ijnum, int dim)
{                        
    int blockDim = 256;
    int gridDim = (ijnum + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;
      
    gpuKernelFullNeighPairDecomposition<<<gridDim, blockDim>>>(e, eij, ai, ijnum, dim);
}
template void gpuFullNeighPairDecomposition(double*, double*, int*, int, int);
template void gpuFullNeighPairDecomposition(float*, float*, int*, int, int);

template <typename T> __global__ void gpuKernelCenterAtomPairDecomposition(T *e, T *eij, int *ilist, 
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
template <typename T> void gpuCenterAtomPairDecomposition(T *e, T *eij, int *ilist, 
        int *anumsum, int inum, int dim)
{                        
    int blockDim = 256;
    int gridDim = (inum + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;
      
    gpuKernelCenterAtomPairDecomposition<<<gridDim, blockDim>>>(e, eij, ilist, anumsum, inum, dim);
}
template void gpuCenterAtomPairDecomposition(double*, double*, int*, int*, int, int);
template void gpuCenterAtomPairDecomposition(float*, float*, int*, int*, int, int);

template <typename T> __global__ void gpuKernelHalfNeighPairDecomposition(T *e, T *eij, int *ai, 
        int *aj, int ijnum, int dim)
{    
    int ii = threadIdx.x + blockIdx.x * blockDim.x;         
    while (ii < ijnum)  {  // for each atom pair ij in the simulation box     
        int i = ai[ii];       // atom i
        int j = aj[ii];       // atom j
        int l = dim*ii;        
        // use atomicAdd on gpu        
        //e[i] += 0.5*eij[ii];
        //e[j] += 0.5*eij[ii];
        atomicAdd(&e[i], 0.5*eij[ii]);                
        atomicAdd(&e[j], 0.5*eij[ii]);                
        ii += blockDim.x * gridDim.x;   
    }
}
template <typename T> void gpuHalfNeighPairDecomposition(T *e, T *eij, int *ai, 
        int *aj, int ijnum, int dim)
{                        
    int blockDim = 256;
    int gridDim = (ijnum + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;
      
    gpuKernelHalfNeighPairDecomposition<<<gridDim, blockDim>>>(e, eij, ai, aj, ijnum, dim);
}
template void gpuHalfNeighPairDecomposition(double*, double*, int*, int*, int, int);
template void gpuHalfNeighPairDecomposition(float*, float*, int*, int*, int, int);

template <typename T> __global__ void gpuKernelNeighborAtomPairDecomposition(T *e, T *eij, int *jlist, 
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
template <typename T> void gpuNeighborAtomPairDecomposition(T *e, T *eij, int *jlist, 
        int *bnumsum, int *index, int jnum, int dim)
{                        
    int blockDim = 256;
    int gridDim = (jnum + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;
      
    gpuKernelNeighborAtomPairDecomposition<<<gridDim, blockDim>>>(e, eij, jlist, bnumsum, index, jnum, dim);
}
template void gpuNeighborAtomPairDecomposition(double*, double*, int*, int*, int*, int, int);
template void gpuNeighborAtomPairDecomposition(float*, float*, int*, int*, int*, int, int);

template <typename T> __global__ void gpuKernelTripletDecomposition(T *e, T *eijk, int *ai, int *aj, 
        int *ak, int ijknum, int dim)
{    
    int ii = threadIdx.x + blockIdx.x * blockDim.x;         
    while (ii < ijknum)  {  // for each atom pair ij in the simulation box     
        int i = ai[ii];       // atom i
        int j = aj[ii];       // atom j
        int k = ak[ii];       // atom k
        // use atomicAdd on gpu
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
template <typename T> void gpuTripletDecomposition(T *e, T *eijk, int *ai, int *aj, 
        int *ak, int ijknum, int dim)
{                        
    int blockDim = 256;
    int gridDim = (ijknum + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;
      
    gpuKernelTripletDecomposition<<<gridDim, blockDim>>>(e, eijk, ai, aj, ak, ijknum, dim);
}
template void gpuTripletDecomposition(double*, double*, int*, int*, int*, int, int);
template void gpuTripletDecomposition(float*, float*, int*, int*, int*, int, int);

template <typename T> __global__ void gpuKernelCenterAtomTripletDecomposition(T *e, T *eijk, int *ilist, 
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
template <typename T> void gpuCenterAtomTripletDecomposition(T *e, T *eijk, int *ilist, 
        int *anumsum, int inum, int dim)
{                        
    int blockDim = 256;
    int gridDim = (inum + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;
      
    gpuKernelCenterAtomTripletDecomposition<<<gridDim, blockDim>>>(e, eijk, ilist, anumsum, inum, dim);
}
template void gpuCenterAtomTripletDecomposition(double*, double*, int*, int*, int, int);
template void gpuCenterAtomTripletDecomposition(float*, float*,  int*, int*, int, int);

template <typename T> __global__ void gpuKernelNeighborAtomTripletDecomposition(T *e, T *eij, int *jlist, 
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
template <typename T> void gpuNeighborAtomTripletDecomposition(T *e, T *eij, int *jlist, 
        int *bnumsum, int *index, int jnum, int dim)
{                        
    int blockDim = 256;
    int gridDim = (jnum + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;
      
    gpuKernelNeighborAtomTripletDecomposition<<<gridDim, blockDim>>>(e, eij, jlist, bnumsum, index, jnum, dim);
}
template void gpuNeighborAtomTripletDecomposition(double*, double*, int*, int*, int*, int, int);
template void gpuNeighborAtomTripletDecomposition(float*, float*, int*, int*, int*, int, int);

template <typename T> __global__ void gpuKernelQuadrupletDecomposition(T *e, T *eijkl, 
        int *ai, int *aj, int *ak, int *al, int ijklnum, int dim)
{    
    int ii = threadIdx.x + blockIdx.x * blockDim.x;         
    while (ii < ijklnum)  {  // for each atom quadruplet ijkl in the simulation box     
        int i = ai[ii];       // atom i
        int j = aj[ii];       // atom j
        int k = ak[ii];       // atom k
        int l = al[ii];       // atom l
        // use atomicAdd on gpu
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
template <typename T> void gpuQuadrupletDecomposition(T *e, T *eijkl, 
        int *ai, int *aj, int *ak, int *al, int ijklnum, int dim)
{                        
    int blockDim = 256;
    int gridDim = (ijklnum + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;
      
    gpuKernelQuadrupletDecomposition<<<gridDim, blockDim>>>(e, eijkl, ai, aj, ak, al, ijklnum, dim);
}
template void gpuQuadrupletDecomposition(double*, double*, int*, int*, int*, int*, int, int);
template void gpuQuadrupletDecomposition(float*, float*, int*, int*, int*, int*, int, int);

template <typename T> __global__ void gpuKernelCenterAtomQuadrupletDecomposition(T *e, T *eijkl, int *ilist, 
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
template <typename T> void gpuCenterAtomQuadrupletDecomposition(T *e, T *eijkl, int *ilist, 
        int *anumsum, int inum, int dim)
{                        
    int blockDim = 256;
    int gridDim = (inum + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;
      
    gpuKernelCenterAtomQuadrupletDecomposition<<<gridDim, blockDim>>>(e, eijkl, ilist, anumsum, inum, dim);
}
template void gpuCenterAtomQuadrupletDecomposition(double*, double*, int*, int*, int, int);
template void gpuCenterAtomQuadrupletDecomposition(float*, float*, int*, int*, int, int);

template <typename T> __global__ void gpuKernelNeighborAtomQuadrupletDecomposition(T *e, T *eij, int *jlist, 
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
template <typename T> void gpuNeighborAtomQuadrupletDecomposition(T *e, T *eij, int *jlist, 
        int *bnumsum, int *index, int jnum, int dim)
{                        
    int blockDim = 256;
    int gridDim = (jnum + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;
      
    gpuKernelNeighborAtomQuadrupletDecomposition<<<gridDim, blockDim>>>(e, eij, jlist, bnumsum, index, jnum, dim);
}
template void gpuNeighborAtomQuadrupletDecomposition(double*, double*, int*, int*, int*, int, int);
template void gpuNeighborAtomQuadrupletDecomposition(float*, float*, int*, int*, int*, int, int);

template <typename T> __global__ void gpuKernelSingleDecomposition(T *e, T *f, T *ei, T *fi, int *ai, int inum, int dim)
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
template <typename T> void gpuSingleDecomposition(T *e, T *f, T *ei, T *fi, int *ai, int inum, int dim)
{                        
    int blockDim = 256;
    int gridDim = (inum + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;
      
    gpuKernelSingleDecomposition<<<gridDim, blockDim>>>(e, f, ei, fi, ai, inum, dim);
}
template void gpuSingleDecomposition(double*, double*, double*, double*, int*, int, int);
template void gpuSingleDecomposition(float*, float*, float*, float*, int*, int, int);

template <typename T> __global__ void gpuKernelCenterAtomPairDecomposition(T *e, T *f, T *eij, 
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
template <typename T> void gpuCenterAtomPairDecomposition(T *e, T *f, T *eij, 
        T *fij, int *ilist, int *anumsum, int inum, int dim)
{                        
    int blockDim = 256;
    int gridDim = (inum + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;
      
    gpuKernelCenterAtomPairDecomposition<<<gridDim, blockDim>>>(e, f, eij, fij, ilist, anumsum, inum, dim);
}
template void gpuCenterAtomPairDecomposition(double*, double*, double*, double*, int*, int*, int, int);
template void gpuCenterAtomPairDecomposition(float*, float*, float*, float*, int*, int*, int, int);

template <typename T> __global__ void gpuKernelNeighborAtomPairDecomposition(T *e, T *f, T *eij, T *fij, 
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
template <typename T> void gpuNeighborAtomPairDecomposition(T *e, T *f, T *eij, T *fij, 
        int *jlist, int *bnumsum, int *index, int jnum, int dim)
{                        
    int blockDim = 256;
    int gridDim = (jnum + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;
      
    gpuKernelNeighborAtomPairDecomposition<<<gridDim, blockDim>>>(e, f, eij, fij, jlist, bnumsum, index, jnum, dim);
}
template void gpuNeighborAtomPairDecomposition(double*, double*, double*, double*, int*, int*, int*, int, int);
template void gpuNeighborAtomPairDecomposition(float*, float*, float*, float*, int*, int*, int*, int, int);

template <typename T> __global__ void gpuKernelFullNeighPairDecomposition(T *e, T *f, T *eij, T *fij, 
        int *ai, int ijnum, int dim)
{    
    int ii = threadIdx.x + blockIdx.x * blockDim.x;         
    while (ii < ijnum)  {  // for each atom pair ij in the simulation box     
        int i = ai[ii];       // atom i
        int l = dim*ii;
        // use atomicAdd on gpu
        //e[i] += 0.5*eij[ii];
        atomicAdd(&e[i], 0.5*eij[ii]);                
        i = dim*i;
        for (int d=0; d<dim; d++) 
            atomicAdd(&f[i+d], -fij[l+d]);                
            //f[i+d] +=  -fij[l+d];        
        ii += blockDim.x * gridDim.x;   
    }
}
template <typename T> void gpuFullNeighPairDecomposition(T *e, T *f, T *eij, T *fij, 
        int *ai, int ijnum, int dim)
{                        
    int blockDim = 256;
    int gridDim = (ijnum + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;
      
    gpuKernelFullNeighPairDecomposition<<<gridDim, blockDim>>>(e, f, eij, fij, ai, ijnum, dim);
}
template void gpuFullNeighPairDecomposition(double*, double*, double*, double*, int*, int, int);
template void gpuFullNeighPairDecomposition(float*, float*, float*, float*, int*, int, int);

template <typename T> __global__ void gpuKernelHalfNeighPairDecomposition(T *e, T *f, T *eij, T *fij, 
        int *ai, int *aj, int ijnum, int dim)
{    
    int ii = threadIdx.x + blockIdx.x * blockDim.x;         
    while (ii < ijnum)  {  // for each atom pair ij in the simulation box     
        int i = ai[ii];       // atom i
        int j = aj[ii];       // atom j
        int l = dim*ii;        
        // use atomicAdd on gpu
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
template <typename T> void gpuHalfNeighPairDecomposition(T *e, T *f, T *eij, T *fij, 
        int *ai, int *aj, int ijnum, int dim)
{                        
    int blockDim = 256;
    int gridDim = (ijnum + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;
      
    gpuKernelHalfNeighPairDecomposition<<<gridDim, blockDim>>>(e, f, eij, fij, ai, aj, ijnum, dim);
}
template void gpuHalfNeighPairDecomposition(double*, double*, double*, double*, int*, int*, int, int);
template void gpuHalfNeighPairDecomposition(float*, float*, float*, float*, int*, int*, int, int);

template <typename T> __global__ void gpuKernelTripletDecomposition(T *e, T *f, T *eijk, T *fij, 
        T *fik, int *ai, int *aj, int *ak, int ijknum, int dim)
{    
    int ii = threadIdx.x + blockIdx.x * blockDim.x;         
    while (ii < ijknum)  {  // for each atom pair ij in the simulation box     
        int i = ai[ii];       // atom i
        int j = aj[ii];       // atom j
        int k = ak[ii];       // atom k
        int n = dim*ii;        
        // use atomicAdd on gpu
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
template <typename T> void gpuTripletDecomposition(T *e, T *f, T *eijk, T *fij, 
        T *fik, int *ai, int *aj, int *ak, int ijknum, int dim)
{                        
    int blockDim = 256;
    int gridDim = (ijknum + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;
      
    gpuKernelTripletDecomposition<<<gridDim, blockDim>>>(e, f, eijk, fij, fik, ai, aj, ak, ijknum, dim);
}
template void gpuTripletDecomposition(double*, double*, double*, double*, double*, int*, int*, int*, int, int);
template void gpuTripletDecomposition(float*, float*, float*, float*, float*, int*, int*, int*, int, int);

template <typename T> __global__ void gpuKernelCenterAtomTripletDecomposition(T *e, T *f, T *eijk, T *fij, 
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
template <typename T> void gpuCenterAtomTripletDecomposition(T *e, T *f, T *eijk, T *fij, 
        T *fik, int *ilist, int *anumsum, int inum, int dim)
{                        
    int blockDim = 256;
    int gridDim = (inum + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;
      
    gpuKernelCenterAtomTripletDecomposition<<<gridDim, blockDim>>>(e, f, eijk, fij, fik, ilist, anumsum, inum, dim);
}
template void gpuCenterAtomTripletDecomposition(double*, double*, double*, double*, double*, int*, int*, int, int);
template void gpuCenterAtomTripletDecomposition(float*, float*, float*, float*, float*, int*, int*, int, int);

template <typename T> __global__ void gpuKernelNeighborAtomTripletDecomposition(T *e, T *f, T *eij, T *fij, 
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
template <typename T> void gpuNeighborAtomTripletDecomposition(T *e, T *f, T *eij, T *fij, 
        int *jlist, int *bnumsum, int *index, int jnum, int dim)
{                        
    int blockDim = 256;
    int gridDim = (jnum + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;
      
    gpuKernelNeighborAtomTripletDecomposition<<<gridDim, blockDim>>>(e, f, eij, fij, jlist, bnumsum, index, jnum, dim);
}
template void gpuNeighborAtomTripletDecomposition(double*, double*, double*, double*, int*, int*, int*, int, int);
template void gpuNeighborAtomTripletDecomposition(float*, float*, float*, float*, int*, int*, int*, int, int);

template <typename T> __global__ void gpuKernelQuadrupletDecomposition(T *e, T *f, T *eijkl, T *fij, T *fik, T *fil, 
        int *ai, int *aj, int *ak, int *al, int ijklnum, int dim)
{    
    int ii = threadIdx.x + blockIdx.x * blockDim.x;         
    while (ii < ijklnum)  {  // for each atom quadruplet ijkl in the simulation box     
        int i = ai[ii];       // atom i
        int j = aj[ii];       // atom j
        int k = ak[ii];       // atom k
        int l = al[ii];       // atom l
        int n = dim*ii;        
        // use atomicAdd on gpu
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
template <typename T> void gpuQuadrupletDecomposition(T *e, T *f, T *eijkl, T *fij, T *fik, T *fil, 
        int *ai, int *aj, int *ak, int *al, int ijklnum, int dim)
{                        
    int blockDim = 256;
    int gridDim = (ijklnum + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;
      
    gpuKernelQuadrupletDecomposition<<<gridDim, blockDim>>>(e, f, eijkl, fij, fik, fil, ai, aj, ak, al, ijklnum, dim);
}
template void gpuQuadrupletDecomposition(double*, double*, double*, double*, double*, double*, int*, int*, int*, int*, int, int);
template void gpuQuadrupletDecomposition(float*, float*, float*, float*, float*, float*, int*, int*, int*, int*, int, int);

template <typename T> __global__ void gpuKernelCenterAtomQuadrupletDecomposition(T *e, T *f, T *eijkl, T *fij, 
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
template <typename T> void gpuCenterAtomQuadrupletDecomposition(T *e, T *f, T *eijkl, T *fij, 
        T *fik, T *fil, int *ilist, int *anumsum, int inum, int dim)
{                        
    int blockDim = 256;
    int gridDim = (inum + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;
      
    gpuKernelCenterAtomQuadrupletDecomposition<<<gridDim, blockDim>>>(e, f, eijkl, fij, fik, fil, ilist, anumsum, inum, dim);
}
template void gpuCenterAtomQuadrupletDecomposition(double*, double*, double*, double*, double*, double*, int*, int*, int, int);
template void gpuCenterAtomQuadrupletDecomposition(float*, float*, float*, float*, float*, float*, int*, int*, int, int);

template <typename T> __global__ void gpuKernelNeighborAtomQuadrupletDecomposition(T *e, T *f, T *eij, T *fij, 
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
template <typename T> void gpuNeighborAtomQuadrupletDecomposition(T *e, T *f, T *eij, T *fij, 
        int *jlist, int *bnumsum, int *index, int jnum, int dim)
{                        
    int blockDim = 256;
    int gridDim = (jnum + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;
      
    gpuKernelNeighborAtomQuadrupletDecomposition<<<gridDim, blockDim>>>(e, f, eij, fij, jlist, bnumsum, index, jnum, dim);
}
template void gpuNeighborAtomQuadrupletDecomposition(double*, double*, double*, double*, int*, int*, int*, int, int);
template void gpuNeighborAtomQuadrupletDecomposition(float*, float*, float*, float*, int*, int*, int*, int, int);


template <typename T> __global__ void gpuKernelTripletForceDecomposition(T *f3ik, T *eij, T *d3ij, 
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
template <typename T> void gpuTripletForceDecomposition(T *f3ik, T *eij, T *d3ij, 
        int *tripletnum, int *tripletnumsum, int npairs, int dim)
{                        
    int blockDim = 256;
    int gridDim = (npairs + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;
      
    gpuKernelTripletForceDecomposition<<<gridDim, blockDim>>>(f3ik, eij, d3ij, tripletnum, tripletnumsum, npairs, dim);
}
template void gpuTripletForceDecomposition(double*, double*, double*, int*, int*, int, int);
template void gpuTripletForceDecomposition(float*, float*, float*, int*, int*, int, int);


template <typename T> __global__ void gpuKernelHalfForceDecomposition2D(T *f, T *fij, int *ai, int *aj, int ijnum)
{    
    int ii = threadIdx.x + blockIdx.x * blockDim.x;     
    while (ii < ijnum)  {  // for each atom pair ij in the simulation box     
        int i = 2*ai[ii];       // atom i
        int j = 2*aj[ii];       // atom j
        int l = 2*ii;
        // use atomicAdd on gpu
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
template <typename T> void gpuHalfForceDecomposition2D(T *f, T *fij, int *ai, int *aj, int ijnum)
{                        
    int blockDim = 256;
    int gridDim = (ijnum + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;
      
    gpuKernelHalfForceDecomposition2D<<<gridDim, blockDim>>>(f, fij, ai, aj, ijnum);
}
template void gpuHalfForceDecomposition2D(double*, double*, int*, int*, int);
template void gpuHalfForceDecomposition2D(float*, float*, int*, int*, int);

template <typename T> __global__ void gpuKernelFullForceDecomposition2D(T *f, T *fij, int *ai, int *aj, int ijnum)
{    
    int ii = threadIdx.x + blockIdx.x * blockDim.x;     
    while (ii < ijnum)  {  // for each atom pair ij in the simulation box     
        int i = 2*ai[ii];       // atom i
        int l = 2*ii;
        // use atomicAdd on gpu
//         f[i+0] +=  -fij[l+0];
//         f[i+1] +=  -fij[l+1];
        atomicAdd(&f[0+i], -fij[0+l]);
        atomicAdd(&f[1+i], -fij[1+l]);
        ii += blockDim.x * gridDim.x;   
    }
}
template <typename T> void gpuFullForceDecomposition2D(T *f, T *fij, int *ai, int *aj, int ijnum)
{                        
    int blockDim = 256;
    int gridDim = (ijnum + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;
      
    gpuKernelFullForceDecomposition2D<<<gridDim, blockDim>>>(f, fij, ai, aj, ijnum);
}
template void gpuFullForceDecomposition2D(double*, double*, int*, int*, int);
template void gpuFullForceDecomposition2D(float*, float*, int*, int*, int);

template <typename T> __global__ void gpuKernelIAtomDecomposition2D(T *f, T *fij, int *ilist, int *anumsum, int inum)
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
template <typename T> void gpuIAtomDecomposition2D(T *f, T *fij, int *ilist, int *anumsum, int inum)
{                        
    int blockDim = 256;
    int gridDim = (inum + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;
      
    gpuKernelIAtomDecomposition2D<<<gridDim, blockDim>>>(f, fij, ilist, anumsum, inum);
}
template void gpuIAtomDecomposition2D(double*, double*, int*, int*, int);
template void gpuIAtomDecomposition2D(float*, float*, int*, int*, int);

template <typename T> __global__ void gpuKernelJAtomDecomposition2D(T *f, T *fij, 
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
template <typename T> void gpuJAtomDecomposition2D(T *f, T *fij, 
        int *jlist, int *bnumsum, int *index, int jnum)
{                        
    int blockDim = 256;
    int gridDim = (jnum + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;
      
    gpuKernelJAtomDecomposition2D<<<gridDim, blockDim>>>(f, fij, jlist, bnumsum, index, jnum);
}
template void gpuJAtomDecomposition2D(double*, double*, int*, int*, int*, int);
template void gpuJAtomDecomposition2D(float*, float*, int*, int*, int*, int);

template <typename T> __global__ void gpuKernelForceDecompositionTriplet2D(T *f, T *fij, T *fik, 
        int *ai, int *aj, int *ak, int ijknum)
{    
    int ii = threadIdx.x + blockIdx.x * blockDim.x;     
    while (ii < ijknum)  {  // for each atom pair ij in the simulation box     
        int i = 2*ai[ii];       // atom i
        int j = 2*aj[ii];       // atom j
        int k = 2*ak[ii];       // atom k
        int l = 2*ii;
        // use atomicAdd on gpu
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
template <typename T> void gpuForceDecompositionTriplet2D(T *f, T *fij, T *fik, 
        int *ai, int *aj, int *ak, int ijknum)
{                        
    int blockDim = 256;
    int gridDim = (ijknum + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;
      
    gpuKernelForceDecompositionTriplet2D<<<gridDim, blockDim>>>(f, fij, fik, ai, aj, ak, ijknum);
}
template void gpuForceDecompositionTriplet2D(double*, double*, double*, int*, int*, int*, int);
template void gpuForceDecompositionTriplet2D(float*, float*, float*, int*, int*, int*, int);

template <typename T> __global__ void gpuKernelAtomDecompositionTriplet2D(T *f, T *fij, T *fik, 
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
//             // use atomicAdd on gpu
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
template <typename T> void gpuAtomDecompositionTriplet2D(T *f, T *fij, T *fik, 
        int *ilist, int *anumsum, int inum)
{                        
    int blockDim = 256;
    int gridDim = (inum + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;
      
    gpuKernelAtomDecompositionTriplet2D<<<gridDim, blockDim>>>(f, fij, fik, ilist, anumsum, inum);
}
template void gpuAtomDecompositionTriplet2D(double*, double*, double*, int*, int*, int);
template void gpuAtomDecompositionTriplet2D(float*, float*, float*, int*, int*, int);

template <typename T> __global__ void gpuKernelForceDecompositionQuadruplet2D(T *f, T *fij, T *fik,  T *fil, 
        int *ai, int *aj, int *ak, int *al, int ijklnum)
{    
    for (int ii=0; ii<ijklnum; ii++) {  // for each atom pair ij in the simulation box     
        int i = 2*ai[ii];       // atom i
        int j = 2*aj[ii];       // atom j
        int k = 2*ak[ii];       // atom k
        int l = 2*al[ii];       // atom l
        int n = 2*ii;
        // use atomicAdd on gpu
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
template <typename T> void gpuForceDecompositionQuadruplet2D(T *f, T *fij, T *fik,  T *fil, 
        int *ai, int *aj, int *ak, int *al, int ijklnum)
{                        
    int blockDim = 256;
    int gridDim = (ijklnum + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;
      
    gpuKernelForceDecompositionQuadruplet2D<<<gridDim, blockDim>>>(f, fij, fik, fil, ai, aj, ak, al, ijklnum);
}
template void gpuForceDecompositionQuadruplet2D(double*, double*, double*, double*, int*, int*,  int*, int*, int);
template void gpuForceDecompositionQuadruplet2D(float*, float*, float*, float*, int*, int*, int*, int*, int);

template <typename T> __global__ void gpuKernelAtomDecompositionQuadruplet2D(T *f, T *fij, T *fik, T *fil,
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
template <typename T> void gpuAtomDecompositionQuadruplet2D(T *f, T *fij, T *fik, T *fil,
        int *ilist, int *anumsum, int inum)
{                        
    int blockDim = 256;
    int gridDim = (inum + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;
      
    gpuKernelAtomDecompositionQuadruplet2D<<<gridDim, blockDim>>>(f, fij, fik, fil, ilist, anumsum, inum);
}
template void gpuAtomDecompositionQuadruplet2D(double*, double*, double*, double*, int*, int*, int);
template void gpuAtomDecompositionQuadruplet2D(float*, float*, float*, float*, int*, int*, int);

template <typename T> __global__ void gpuKernelHalfForceDecomposition3D(T *f, T *fij, int *ai, int *aj, int ijnum)
{    
    int ii = threadIdx.x + blockIdx.x * blockDim.x;     
    while (ii < ijnum)  {  // for each atom pair ij in the simulation box     
        int i = 3*ai[ii];       // atom i
        int j = 3*aj[ii];       // atom j
        int l = 3*ii;
        // use atomicAdd on gpu
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
template <typename T> void gpuHalfForceDecomposition3D(T *f, T *fij, int *ai, int *aj, int ijnum)
{                        
    int blockDim = 256;
    int gridDim = (ijnum + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;
      
    gpuKernelHalfForceDecomposition3D<<<gridDim, blockDim>>>(f, fij, ai, aj, ijnum);
}
template void gpuHalfForceDecomposition3D(double*, double*, int*, int*, int);
template void gpuHalfForceDecomposition3D(float*, float*, int*, int*, int);

template <typename T> __global__ void gpuKernelFullForceDecomposition3D(T *f, T *fij, int *ai, int *aj, int ijnum)
{    
    int ii = threadIdx.x + blockIdx.x * blockDim.x;     
    while (ii < ijnum)  {  // for each atom pair ij in the simulation box     
        int i = 3*ai[ii];             // atom i
        int l = 3*ii;
        // use atomicAdd on gpu
//         f[i+0] +=  -fij[l+0];
//         f[i+1] +=  -fij[l+1];
//         f[i+2] +=  -fij[l+2];                  
        atomicAdd(&f[0+i], -fij[0+l]);
        atomicAdd(&f[1+i], -fij[1+l]);
        atomicAdd(&f[2+i], -fij[2+l]);                
        ii += blockDim.x * gridDim.x;   
    }
}
template <typename T> void gpuFullForceDecomposition3D(T *f, T *fij, int *ai, int *aj, int ijnum)
{                        
    int blockDim = 256;
    int gridDim = (ijnum + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;
      
    gpuKernelFullForceDecomposition3D<<<gridDim, blockDim>>>(f, fij, ai, aj, ijnum);
}
template void gpuFullForceDecomposition3D(double*, double*, int*, int*, int);
template void gpuFullForceDecomposition3D(float*, float*, int*, int*, int);

template <typename T> __global__ void gpuKernelIAtomDecomposition3D(T *f, T *fij, int *ilist, int *anumsum, int inum)
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
template <typename T> void gpuIAtomDecomposition3D(T *f, T *fij, int *ilist, int *anumsum, int inum)
{                        
    int blockDim = 256;
    int gridDim = (inum + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;
      
    gpuKernelIAtomDecomposition3D<<<gridDim, blockDim>>>(f, fij, ilist, anumsum, inum);
}
template void gpuIAtomDecomposition3D(double*, double*, int*, int*, int);
template void gpuIAtomDecomposition3D(float*, float*, int*, int*, int);

template <typename T> __global__ void gpuKernelJAtomDecomposition3D(T *f, T *fij, int *jlist, 
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
template <typename T> void gpuJAtomDecomposition3D(T *f, T *fij, 
        int *jlist, int *bnumsum, int *index, int jnum)
{                        
    int blockDim = 256;
    int gridDim = (jnum + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;
      
    gpuKernelJAtomDecomposition3D<<<gridDim, blockDim>>>(f, fij, jlist, bnumsum, index, jnum);
}
template void gpuJAtomDecomposition3D(double*, double*, int*, int*, int*, int);
template void gpuJAtomDecomposition3D(float*, float*, int*, int*, int*, int);

template <typename T> __global__ void gpuKernelForceDecompositionTriplet3D(T *f, T *fij, T *fik, 
        int *ai, int *aj, int *ak, int ijknum)
{    
    int ii = threadIdx.x + blockIdx.x * blockDim.x;     
    while (ii < ijknum)  {  // for each atom pair ij in the simulation box     
        int i = 3*ai[ii];       // atom i
        int j = 3*aj[ii];       // atom j
        int k = 3*ak[ii];       // atom k
        int l = 3*ii;
        // use atomicAdd on gpu
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
template <typename T> void gpuForceDecompositionTriplet3D(T *f, T *fij, T *fik, 
        int *ai, int *aj, int *ak, int ijknum)
{                        
    int blockDim = 256;
    int gridDim = (ijknum + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;
      
    gpuKernelForceDecompositionTriplet3D<<<gridDim, blockDim>>>(f, fij, fik, ai, aj, ak, ijknum);
}
template void gpuForceDecompositionTriplet3D(double*, double*, double*, int*, int*, int*, int);
template void gpuForceDecompositionTriplet3D(float*, float*, float*, int*, int*, int*, int);

template <typename T> __global__ void gpuKernelAtomDecompositionTriplet3D(T *f, T *fij, T *fik, 
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
template <typename T> void gpuAtomDecompositionTriplet3D(T *f, T *fij, T *fik, 
        int *ilist, int *anumsum, int inum)
{                        
    int blockDim = 256;
    int gridDim = (inum + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;
      
    gpuKernelAtomDecompositionTriplet3D<<<gridDim, blockDim>>>(f, fij, fik, ilist, anumsum, inum);
}
template void gpuAtomDecompositionTriplet3D(double*, double*, double*, int*, int*, int);
template void gpuAtomDecompositionTriplet3D(float*, float*, float*, int*, int*, int);

template <typename T> __global__ void gpuKernelForceDecompositionQuadruplet3D(T *f, T *fij, T *fik,  
        T *fil, int *ai, int *aj, int *ak, int *al, int ijklnum)
{    
    int ii = threadIdx.x + blockIdx.x * blockDim.x;     
    while (ii < ijklnum)  {  // for each atom pair ij in the simulation box     
        int i = 3*ai[ii];       // atom i
        int j = 3*aj[ii];       // atom j
        int k = 3*ak[ii];       // atom k
        int l = 3*al[ii];       // atom l
        int n = 3*ii;
        // use atomicAdd on gpu
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
template <typename T> void gpuForceDecompositionQuadruplet3D(T *f, T *fij, T *fik,  T *fil, 
        int *ai, int *aj, int *ak, int *al, int ijklnum)
{                        
    int blockDim = 256;
    int gridDim = (ijklnum + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;
      
    gpuKernelForceDecompositionQuadruplet3D<<<gridDim, blockDim>>>(f, fij, fik, fil, ai, aj, ak, al, ijklnum);
}
template void gpuForceDecompositionQuadruplet3D(double*, double*, double*, double*, int*, int*,  int*, int*, int);
template void gpuForceDecompositionQuadruplet3D(float*, float*, float*, float*, int*, int*, int*, int*, int);

template <typename T> __global__ void gpuKernelAtomDecompositionQuadruplet3D(T *f, T *fij, T *fik, T *fil, 
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
template <typename T> void gpuAtomDecompositionQuadruplet3D(T *f, T *fij, T *fik, T *fil,
        int *ilist, int *anumsum, int inum)
{                        
    int blockDim = 256;
    int gridDim = (inum + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;
      
    gpuKernelAtomDecompositionQuadruplet3D<<<gridDim, blockDim>>>(f, fij, fik, fil, ilist, anumsum, inum);
}
template void gpuAtomDecompositionQuadruplet3D(double*, double*, double*, double*, int*, int*, int);
template void gpuAtomDecompositionQuadruplet3D(float*, float*, float*, float*, int*, int*, int);

#endif


