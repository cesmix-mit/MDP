#ifndef __OPUAPP
#define __OPUAPP

#include <math.h>

#include "opuSinglea.cpp"
#include "opuSingleb.cpp"
#include "opuPaira.cpp"
#include "opuPairb.cpp"
#include "opuPairc.cpp"
#include "opuTripleta.cpp"
#include "opuTripletb.cpp"
#include "opuTripletc.cpp"
#include "opuQuadrupleta.cpp"
#include "opuQuadrupletb.cpp"

// TO DO: Use Enzyme to calculate the potential derivatives

void cpuComputeSingleEnergyForce(dstype *ei, dstype *fi, dstype *xi, dstype *qi, int *ti, int *ai, 
        dstype *mu, dstype *eta, int *kappa, int dim, int ncq, int nmu, int neta, int nkappa, 
        int inum, int potnum, int bondtype)        
{
    if (bondtype==0)  {
        opuSinglea(ei, xi, qi, ti, ai, mu, eta, kappa, dim, ncq, nmu, neta, nkappa, inum, potnum);    
    }
    
    if (bondtype==1) {
        opuSingleb(ei, xi, qi, ti, ai, mu, eta, kappa, dim, ncq, nmu, neta, nkappa, inum, potnum);    
    }
}

template <typename T> void opuLJ1pot(T *u, T *f, T *xij, T *qi, T *qj, int *ti, int *tj, int *ai, int *aj, T *mu, T *eta, int *kappa, int dim, int ncq, int nmu, int neta, int nkappa, int ng, int potnum)
{
	for (int i = 0; i <ng; i++) {
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
	}	
}
template void opuLJ1pot(double *, double *, double *, double *, double *, int *, int *, int *, int *, double *, double *, int*, int, int, int, int, int, int, int);
template void opuLJ1pot(float *, float *, float *, float *, float *, int *, int *, int *, int *, float *, float *, int *, int, int, int, int, int, int, int);

template <typename T> void opuLJ2pot(T *u, T *f, T *xij, T *qi, T *qj, int *ti, int *tj, int *ai, int *aj, T *mu, T *eta, int *kappa, int dim, int ncq, int nmu, int neta, int nkappa, int ng, int potnum)
{
	for (int i = 0; i <ng; i++) {                
        int ntype = kappa[0]; // number of atom types
        int a = ti[i]-1 + ntype*(tj[i]-1);  // ti-tj interaction      
		T mu1 = mu[0 + 2*a];
		T mu2 = mu[1 + 2*a];
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
	}	
}
template void opuLJ2pot(double *, double *, double *, double *, double *, int *, int *, int *, int *, double *, double *, int*, int, int, int, int, int, int, int);
template void opuLJ2pot(float *, float *, float *, float *, float *, int *, int *, int *, int *, float *, float *, int *, int, int, int, int, int, int, int);

void cpuComputePairEnergyForce(dstype *eij, dstype *fij, dstype *xij, dstype *qi, dstype *qj, 
        int *ti, int *tj, int *ai, int *aj, dstype *mu, dstype *eta, int *kappa, int dim, 
        int ncq, int nmu, int neta, int nkappa, int ijnum, int potnum, int bondtype)
{
    if (bondtype==0) {
        //opuPaira(eij, xij, qi, qj, ti, tj, ai, aj, mu, eta, kappa, dim, ncq, nmu, neta, nkappa, ijnum, potnum);    
        opuLJ1pot(eij, fij, xij, qi, qj, ti, tj, ai, aj, mu, eta, kappa, dim, ncq, nmu, neta, nkappa, ijnum, potnum);    
    }
    
    if (bondtype==1) {
        opuPairb(eij, xij, qi, qj, ti, tj, ai, aj, mu, eta, kappa, dim, ncq, nmu, neta, nkappa, ijnum, potnum);    
    }
    
    if (bondtype==2) {
        opuPairc(eij, xij, qi, qj, ti, tj, ai, aj, mu, eta, kappa, dim, ncq, nmu, neta, nkappa, ijnum, potnum);    
    }
    
    if (bondtype==3) {
        opuTripletcPair(eij, xij, qi, qj, ti, tj, ai, aj, mu, eta, kappa, dim, ncq, nmu, neta, nkappa, ijnum, potnum);    
    }
}

void cpuElectronDensity(dstype *rhoi, dstype *rhoij, int *pairnum, int *pairnumsum, int inum) 
{
    for (int i=0; i<inum; i++) {
        int jnum = pairnum[i];
        int start = pairnumsum[i];
        rhoi[i] = 0.0;
        for (int j=0; j<jnum; j++) 
            rhoi[i] += rhoij[start+j];        
    }
}
void cpuPaircDensity(dstype *u, dstype *urho, dstype *rho, dstype *mu, dstype *eta, int *kappa, 
        int nrho, int nmu, int neta, int nkappa, int ng, int potnum)
{
    opuPaircDensity(u, rho, mu, eta, kappa, 1, nmu, neta, nkappa, ng, potnum);
}
void cpuEmbedingForce(dstype *fij, dstype *drhoi, int *pairnum, int *pairnumsum, int inum)
{    
    for (int i=0; i<inum; i++) {
        int jnum = pairnum[i];
        int start = pairnumsum[i];
        for (int j=0; j<jnum; j++)             
            fij[start+j] = drhoi[i]*fij[start+j];                
    }        
}
        
void cpuComputeTripletEnergyForce(dstype *eijk, dstype *fij, dstype *fik, dstype *xij, dstype *xik, dstype *qi, dstype *qj, dstype *qk, 
     int *ti, int *tj, int *tk, int *ai, int *aj, int *ak, dstype *mu, dstype *eta, int *kappa, int dim, int ncq, 
     int nmu, int neta, int nkappa, int ijknum, int potnum, int bondtype)        
{
    if (bondtype==0) {
        opuTripleta(eijk, xij, xik, qi, qj, qk, ti, tj, tk, ai, aj, ak, mu, eta, kappa, dim, ncq, nmu, neta, nkappa, ijknum, potnum);        
    }
    
    if (bondtype==1) {
        opuTripletb(eijk, xij, xik, qi, qj, qk, ti, tj, tk, ai, aj, ak, mu, eta, kappa, dim, ncq, nmu, neta, nkappa, ijknum, potnum);
    }
    
    if (bondtype==2) {
        opuTripletc(eijk, xij, xik, qi, qj, qk, ti, tj, tk, ai, aj, ak, mu, eta, kappa, dim, ncq, nmu, neta, nkappa, ijknum, potnum);
    }
}
void cpuTripletcDensity(dstype *u, dstype *urho, dstype *rho, dstype *mu, dstype *eta, int *kappa, 
        int nrho, int nmu, int neta, int nkappa, int ng, int potnum)
{
    opuTripletcDensity(u, rho, mu, eta, kappa, 1, nmu, neta, nkappa, ng, potnum);
}

void cpuComputeQuadrupletEnergyForce(dstype *eijkl, dstype *fij, dstype *fik, dstype *fil, dstype *xij, dstype *xik, dstype *xil, 
      dstype *qi, dstype *qj, dstype *qk, dstype *ql, int *ti, int *tj, int *tk, int *tl, int *ai, int *aj, int *ak, int *al,
      dstype *mu, dstype *eta, int *kappa, int dim, int ncq, int nmu, int neta, int nkappa, int ijklnum, int potnum, int bondtype)        
{
    if (bondtype==0) {
        opuQuadrupleta(eijkl, xij, xik, xil, qi, qj, qk, ql, ti, tj, tk, tl, ai, aj, ak, al,
                mu, eta, kappa, dim, ncq, nmu, neta, nkappa, ijklnum, potnum);        
    }
    
    if (bondtype==1) {
        opuQuadrupletb(eijkl, xij, xik, xil, qi, qj, qk, ql, ti, tj, tk, tl, ai, aj, ak, al, 
                mu, eta, kappa, dim, ncq, nmu, neta, nkappa, ijklnum, potnum);
    }        
}

        
#endif


