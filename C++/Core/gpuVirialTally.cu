/***************************************************************************                               
                    Molecular Dynamics Potentials (MDP)
                           CESMIX-MIT Project  
 
 Contributing authors: Ngoc-Cuong Nguyen (cuongng@mit.edu, exapde@gmail.com)
 ***************************************************************************/

/* ----------------------------------------------------------------------
   tally virial into per-atom accumulators for pair interactions
   pair virial = riFi = ri*fi (pairwise forces have been accumulated in fi)
   pair virial = riFi + rjFj = (ri-rj) Fi = rij*fij
 ------------------------------------------------------------------------- */
template <typename T> __global__ void gpuKernelVirialPairTally(T *vatom, T *fij, T *rij, T factor, 
        int *ai, int dim, int ijnum)
{ 
    int ii = threadIdx.x + blockIdx.x * blockDim.x;
    while (ii < ijnum) {
        int i = ai[ii];        
        T dx = rij[0+dim*ii];
        T dy = rij[1+dim*ii];
        T dz = (dim==3) ? rij[2+dim*ii] : 0.0;    
        T fx = fij[0+dim*ii];
        T fy = fij[1+dim*ii];
        T fz = (dim==3) ? fij[2+dim*ii] : 0.0;    
        vatom[0+6*i] += factor*dx*fx;
        vatom[1+6*i] += factor*dy*fy;
        vatom[2+6*i] += factor*dz*fz;
        vatom[3+6*i] += factor*dx*fy;
        vatom[4+6*i] += factor*dx*fz;
        vatom[5+6*i] += factor*dy*fz;        

        ii += blockDim.x * gridDim.x;
    }
}
template <typename T> void gpuVirialPairTally(T *vatom, T *fij, T *rij, T factor, 
        int *ai, int dim, int ijnum)
{
    int blockDim = 256;
    int gridDim = (ijnum + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;
    gpuKernelVirialPairTally<<<gridDim, blockDim>>>(vatom, fij, rij, factor, 
        ai, dim, ijnum);
}
template void gpuVirialPairTally(double*, double*, double*, double, int*, int, int);
template void gpuVirialPairTally(float*, float*, float*, float, int*, int, int);

/* ----------------------------------------------------------------------
   tally virial into per-atom accumulators for pair interactions
   virial = riFi + rjFj = (ri-rj) Fi 
 ------------------------------------------------------------------------- */
template <typename T> __global__ void gpuKernelVirialPairTally(T *vatom, T *fij, T *rij, T factor, 
        int *ai, int *aj, int dim, int ijnum)
{ 
    int ii = threadIdx.x + blockIdx.x * blockDim.x;
    while (ii < ijnum) {
        int i = ai[ii];        
        int j = aj[ii];        
        T dx = rij[0+dim*ii];
        T dy = rij[1+dim*ii];
        T dz = (dim==3) ? rij[2+dim*ii] : 0.0;    
        T fx = fij[0+dim*ii];
        T fy = fij[1+dim*ii];
        T fz = (dim==3) ? fij[2+dim*ii] : 0.0;    
        T v0 = factor*dx*fx;
        T v1 = factor*dy*fy;
        T v2 = factor*dz*fz;
        T v3 = factor*dx*fy;
        T v4 = factor*dx*fz;
        T v5 = factor*dy*fz;      
        
        vatom[0+6*i] += v0; vatom[1+6*i] += v1;
        vatom[2+6*i] += v2; vatom[3+6*i] += v3;
        vatom[4+6*i] += v4; vatom[5+6*i] += v5;        
        
        vatom[0+6*j] += v0; vatom[1+6*j] += v1;
        vatom[2+6*j] += v2; vatom[3+6*j] += v3;
        vatom[4+6*j] += v4; vatom[5+6*j] += v5;        

        ii += blockDim.x * gridDim.x;
    }
}
template <typename T> void gpuVirialPairTally(T *vatom, T *fij, T *rij, T factor, 
        int *ai, int *aj, int dim, int ijnum)
{
    int blockDim = 256;
    int gridDim = (ijnum + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;
    gpuKernelVirialPairTally<<<gridDim, blockDim>>>(vatom, fij, rij, factor, 
        ai, aj, dim, ijnum);
}
template void gpuVirialPairTally(double*, double*, double*, double, int*, int*, int, int);
template void gpuVirialPairTally(float*, float*, float*, float, int*, int*, int, int);

/* ----------------------------------------------------------------------
   tally virial into per-atom accumulators for triplet interactions
   virial = riFi + rjFj + rkFk = (rj-ri) Fj + (rk-ri) Fk = rij*fij + rik*fik
 ------------------------------------------------------------------------- */
template <typename T> __global__ void gpuKernelVirialTripletTally(T *vatom, T *fij, T *fik, T *rij, 
        T *rik, T factor, int *ai, int *aj, int *ak, int dim, int ijnum)
{ 
    int ii = threadIdx.x + blockIdx.x * blockDim.x;
    while (ii < ijnum) {
        int i = ai[ii];        
        int j = aj[ii];        
        int k = ak[ii];        
        T v0,v1,v3;
        T v2 = 0.0;
        T v4 = 0.0;
        T v5 = 0.0;
        v0 = factor*(rij[0]*fij[0] + rik[0]*fik[0]);
        v1 = factor*(rij[1]*fij[1] + rik[1]*fik[1]);        
        v3 = factor*(rij[0]*fij[1] + rik[0]*fik[1]);        
        if (dim==3) {
            v2 = factor*(rij[2]*fij[2] + rik[2]*fik[2]);
            v4 = factor*(rij[0]*fij[2] + rik[0]*fik[2]);
            v5 = factor*(rij[1]*fij[2] + rik[1]*fik[2]);
        }                
        vatom[0+6*i] += v0; vatom[1+6*i] += v1;
        vatom[2+6*i] += v2; vatom[3+6*i] += v3;
        vatom[4+6*i] += v4; vatom[5+6*i] += v5;        
        
        vatom[0+6*j] += v0; vatom[1+6*j] += v1;
        vatom[2+6*j] += v2; vatom[3+6*j] += v3;
        vatom[4+6*j] += v4; vatom[5+6*j] += v5;        
        
        vatom[0+6*k] += v0; vatom[1+6*k] += v1;
        vatom[2+6*k] += v2; vatom[3+6*k] += v3;
        vatom[4+6*k] += v4; vatom[5+6*k] += v5;        

        ii += blockDim.x * gridDim.x;
    }
}
template <typename T> void gpuVirialTripletTally(T *vatom, T *fij, T *fik, T *rij, 
        T *rik, T factor, int *ai, int *aj, int *ak, int dim, int ijnum)
{
    int blockDim = 256;
    int gridDim = (ijnum + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;
    gpuKernelVirialTripletTally<<<gridDim, blockDim>>>(vatom, fij, fik, rij, rik, factor, 
        ai, aj, ak, dim, ijnum);
}
template void gpuVirialTripletTally(double*, double*, double*, double*, double*,
        double, int*, int*, int*, int, int);
template void gpuVirialTripletTally(float*, float*, float*,  float*, float*, 
        float, int*, int*, int*, int, int);

/* ----------------------------------------------------------------------
   tally virial into per-atom accumulators for quadruplet interactions
   virial = riFi + rjFj + rkFk + rmFm = (rj-ri) Fj + (rk-ri) Fk + (rm-ri) Fm = rij*fij + rik*fik + rim*fim
 ------------------------------------------------------------------------- */
template <typename T> __global__ void gpuKernelVirialQuadrupletTally(T *vatom, T *fij, T *fik, T *fim, 
        T *rij, T *rik, T *rim, T factor, int *ai, int *aj, int *ak, int *am, int dim, int ijnum)
{ 
    int ii = threadIdx.x + blockIdx.x * blockDim.x;
    while (ii < ijnum) {
        int i = ai[ii];        
        int j = aj[ii];        
        int k = ak[ii];        
        int m = ak[ii];        
        T v0,v1,v3;
        T v2 = 0.0;
        T v4 = 0.0;
        T v5 = 0.0;
        v0 = factor*(rij[0]*fij[0] + rik[0]*fik[0] + rim[0]*fim[0]);
        v1 = factor*(rij[1]*fij[1] + rik[1]*fik[1] + rim[1]*fim[1]);        
        v3 = factor*(rij[0]*fij[1] + rik[0]*fik[1] + rim[0]*fim[1]);        
        if (dim==3) {
            v2 = factor*(rij[2]*fij[2] + rik[2]*fik[2] + rim[2]*fim[2]);
            v4 = factor*(rij[0]*fij[2] + rik[0]*fik[2] + rim[0]*fim[2]);
            v5 = factor*(rij[1]*fij[2] + rik[1]*fik[2] + rim[1]*fim[2]);
        }                
        vatom[0+6*i] += v0; vatom[1+6*i] += v1;
        vatom[2+6*i] += v2; vatom[3+6*i] += v3;
        vatom[4+6*i] += v4; vatom[5+6*i] += v5;        
        
        vatom[0+6*j] += v0; vatom[1+6*j] += v1;
        vatom[2+6*j] += v2; vatom[3+6*j] += v3;
        vatom[4+6*j] += v4; vatom[5+6*j] += v5;        
        
        vatom[0+6*k] += v0; vatom[1+6*k] += v1;
        vatom[2+6*k] += v2; vatom[3+6*k] += v3;
        vatom[4+6*k] += v4; vatom[5+6*k] += v5;        
        
        vatom[0+6*m] += v0; vatom[1+6*m] += v1;
        vatom[2+6*m] += v2; vatom[3+6*m] += v3;
        vatom[4+6*m] += v4; vatom[5+6*m] += v5;        

        ii += blockDim.x * gridDim.x;
    }
}
template <typename T> void gpuVirialQuadrupletTally(T *vatom, T *fij, T *fik, T *fim, 
        T *rij, T *rik, T *rim, T factor, int *ai, int *aj, int *ak, int *am, int dim, int ijnum)
{
    int blockDim = 256;
    int gridDim = (ijnum + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;
    gpuKernelVirialQuadrupletTally<<<gridDim, blockDim>>>(vatom, fij, fik, fim, rij, rik, rim, factor, 
        ai, aj, ak, am, dim, ijnum);
}
template void gpuVirialQuadrupletTally(double*, double*, double*, double*, double*, double*, double*,
        double, int*, int*, int*, int*, int, int);
template void gpuVirialQuadrupletTally(float*, float*, float*, float*, float*, float*, float*,  
        float, int*, int*, int*, int*, int, int);

