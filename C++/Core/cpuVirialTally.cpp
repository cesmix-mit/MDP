/***************************************************************************                               
                    Molecular Dynamics Potentials (MDP)
                           CESMIX-MIT Project  
 
 Contributing authors: Ngoc-Cuong Nguyen (cuongng@mit.edu, exapde@gmail.com)
 ***************************************************************************/

// template <typename T> void cpuVirialFdotR(T *vatom, T *f, T *x, int dim, int inum)
// { 
//     for(int ii = 0; ii < inum; ii++) {
//         T dx = x[0+dim*ii];
//         T dy = x[1+dim*ii];
//         T dz = (dim==3) ? x[2+dim*ii] : 0.0;    
//         T fx = f[0+dim*ii];
//         T fy = f[1+dim*ii];
//         T fz = (dim==3) ? f[2+dim*ii] : 0.0;    
//         vatom[0+6*ii] += dx*fx;
//         vatom[1+6*ii] += dy*fy;
//         vatom[2+6*ii] += dz*fz;
//         vatom[3+6*ii] += dx*fy;
//         vatom[4+6*ii] += dx*fz;
//         vatom[5+6*ii] += dy*fz;        
//     }
// }
// template void cpuVirialFdotR(double*, double*, double*, int, int);
// template void cpuVirialFdotR(float*, float*, float*, int, int);

/* ----------------------------------------------------------------------
   tally virial into per-atom accumulators for pair interactions
   pair virial = riFi = ri*fi (pairwise forces have been accumulated in fi)
   pair virial = riFi + rjFj = (ri-rj) Fi = rij*fij
 *cpuVirialPairTally(v, fij, xij, -0.5, ai, dim, common.inum, ntuples);        
 ------------------------------------------------------------------------- */
template <typename T> void cpuVirialPairTally(T *vatom, T *fij, T *rij, T factor, 
        int *ai, int dim, int inum, int ijnum)
{ 
    for(int ii = 0; ii < ijnum; ii++) {
       int i = ai[ii];        
        T dx = rij[0+dim*ii];
        T dy = rij[1+dim*ii];
        T dz = (dim==3) ? rij[2+dim*ii] : 0.0;    
        T fx = fij[0+dim*ii];
        T fy = fij[1+dim*ii];
        T fz = (dim==3) ? fij[2+dim*ii] : 0.0;    
//         vatom[0+6*i] += factor*dx*fx;
//         vatom[1+6*i] += factor*dy*fy;
//         vatom[2+6*i] += factor*dz*fz;
//         vatom[3+6*i] += factor*dx*fy;
//         vatom[4+6*i] += factor*dx*fz;
//         vatom[5+6*i] += factor*dy*fz;        
        vatom[0*inum+i] += factor*dx*fx;
        vatom[1*inum+i] += factor*dy*fy;
        vatom[2*inum+i] += factor*dz*fz;
        vatom[3*inum+i] += factor*dx*fy;
        vatom[4*inum+i] += factor*dx*fz;
        vatom[5*inum+i] += factor*dy*fz;                
    }
}
template void cpuVirialPairTally(double*, double*, double*, double, int*, int, int, int);
template void cpuVirialPairTally(float*, float*, float*, float, int*, int, int, int);

/* ----------------------------------------------------------------------
   tally virial into per-atom accumulators for pair interactions
   virial = riFi + rjFj = (ri-rj) Fi 
 ------------------------------------------------------------------------- */
template <typename T> void cpuVirialPairTally(T *vatom, T *fij, T *rij, T factor, 
        int *ai, int *aj, int dim, int inum, int ijnum)
{ 
    for(int ii = 0; ii < ijnum; ii++) {
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
        
//         vatom[0+6*i] += v0; vatom[1+6*i] += v1;
//         vatom[2+6*i] += v2; vatom[3+6*i] += v3;
//         vatom[4+6*i] += v4; vatom[5+6*i] += v5;        
//         
//         vatom[0+6*j] += v0; vatom[1+6*j] += v1;
//         vatom[2+6*j] += v2; vatom[3+6*j] += v3;
//         vatom[4+6*j] += v4; vatom[5+6*j] += v5;       
        vatom[0*inum+i] += v0; vatom[1*inum+i] += v1;
        vatom[2*inum+i] += v2; vatom[3*inum+i] += v3;
        vatom[4*inum+i] += v4; vatom[5*inum+i] += v5;        
        
        vatom[0*inum+j] += v0; vatom[1*inum+j] += v1;
        vatom[2*inum+j] += v2; vatom[3*inum+j] += v3;
        vatom[4*inum+j] += v4; vatom[5*inum+j] += v5;                
    }
}
template void cpuVirialPairTally(double*, double*, double*, double, int*, int*, int, int, int);
template void cpuVirialPairTally(float*, float*, float*, float, int*, int*, int, int, int);

/* ----------------------------------------------------------------------
   tally virial into per-atom accumulators for triplet interactions
   virial = riFi + rjFj + rkFk = (rj-ri) Fj + (rk-ri) Fk = rij*fij + rik*fik
 ------------------------------------------------------------------------- */
template <typename T> void cpuVirialTripletTally(T *vatom, T *fij, T *fik, T *rij, T *rik, T factor, 
        int *ai, int *aj, int *ak, int dim, int inum, int ijknum)
{ 
    for(int ii = 0; ii < ijknum; ii++) {
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
//         vatom[0+6*i] += v0; vatom[1+6*i] += v1;
//         vatom[2+6*i] += v2; vatom[3+6*i] += v3;
//         vatom[4+6*i] += v4; vatom[5+6*i] += v5;        
//         
//         vatom[0+6*j] += v0; vatom[1+6*j] += v1;
//         vatom[2+6*j] += v2; vatom[3+6*j] += v3;
//         vatom[4+6*j] += v4; vatom[5+6*j] += v5;        
//         
//         vatom[0+6*k] += v0; vatom[1+6*k] += v1;
//         vatom[2+6*k] += v2; vatom[3+6*k] += v3;
//         vatom[4+6*k] += v4; vatom[5+6*k] += v5;        
        vatom[0*inum+i] += v0; vatom[1*inum+i] += v1;
        vatom[2*inum+i] += v2; vatom[3*inum+i] += v3;
        vatom[4*inum+i] += v4; vatom[5*inum+i] += v5;        
        
        vatom[0*inum+j] += v0; vatom[1*inum+j] += v1;
        vatom[2*inum+j] += v2; vatom[3*inum+j] += v3;
        vatom[4*inum+j] += v4; vatom[5*inum+j] += v5;        
        
        vatom[0*inum+k] += v0; vatom[1*inum+k] += v1;
        vatom[2*inum+k] += v2; vatom[3*inum+k] += v3;
        vatom[4*inum+k] += v4; vatom[5*inum+k] += v5;                
    }
}
template void cpuVirialTripletTally(double*, double*, double*, double*, double*,
        double, int*, int*, int*, int, int, int);
template void cpuVirialTripletTally(float*, float*, float*,  float*, float*, 
        float, int*, int*, int*, int, int, int);

/* ----------------------------------------------------------------------
   tally virial into per-atom accumulators for quadruplet interactions
   virial = riFi + rjFj + rkFk + rmFm = (rj-ri) Fj + (rk-ri) Fk + (rm-ri) Fm = rij*fij + rik*fik + rim*fim
 ------------------------------------------------------------------------- */
template <typename T> void cpuVirialQuadrupletTally(T *vatom, T *fij, T *fik, T *fim, 
        T *rij, T *rik, T *rim, T factor, int *ai, int *aj, int *ak, int *am, int dim, int inum, int ijnum)
{ 
    for(int ii = 0; ii < ijnum; ii++) {
        int i = ai[ii];        
        int j = aj[ii];        
        int k = ak[ii];        
        int m = am[ii];        
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
//         vatom[0+6*i] += v0; vatom[1+6*i] += v1;
//         vatom[2+6*i] += v2; vatom[3+6*i] += v3;
//         vatom[4+6*i] += v4; vatom[5+6*i] += v5;        
//         
//         vatom[0+6*j] += v0; vatom[1+6*j] += v1;
//         vatom[2+6*j] += v2; vatom[3+6*j] += v3;
//         vatom[4+6*j] += v4; vatom[5+6*j] += v5;        
//         
//         vatom[0+6*k] += v0; vatom[1+6*k] += v1;
//         vatom[2+6*k] += v2; vatom[3+6*k] += v3;
//         vatom[4+6*k] += v4; vatom[5+6*k] += v5;        
//         
//         vatom[0+6*m] += v0; vatom[1+6*m] += v1;
//         vatom[2+6*m] += v2; vatom[3+6*m] += v3;
//         vatom[4+6*m] += v4; vatom[5+6*m] += v5;        
        vatom[0*inum+i] += v0; vatom[1*inum+i] += v1;
        vatom[2*inum+i] += v2; vatom[3*inum+i] += v3;
        vatom[4*inum+i] += v4; vatom[5*inum+i] += v5;        
        
        vatom[0*inum+j] += v0; vatom[1*inum+j] += v1;
        vatom[2*inum+j] += v2; vatom[3*inum+j] += v3;
        vatom[4*inum+j] += v4; vatom[5*inum+j] += v5;        
        
        vatom[0*inum+k] += v0; vatom[1*inum+k] += v1;
        vatom[2*inum+k] += v2; vatom[3*inum+k] += v3;
        vatom[4*inum+k] += v4; vatom[5*inum+k] += v5;        
        
        vatom[0*inum+m] += v0; vatom[1*inum+m] += v1;
        vatom[2*inum+m] += v2; vatom[3*inum+m] += v3;
        vatom[4*inum+m] += v4; vatom[5*inum+m] += v5;                
    }
}
template void cpuVirialQuadrupletTally(double*, double*, double*, double*, double*, double*, double*,
        double, int*, int*, int*, int*, int, int, int);
template void cpuVirialQuadrupletTally(float*, float*, float*, float*, float*, float*, float*,  
        float, int*, int*, int*, int*, int, int, int);

