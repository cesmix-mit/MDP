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

#ifdef HAVE_ENZYME                
template <typename... Args>
void __enzyme_autodiff(void*, Args... args);
int enzyme_const, enzyme_dup;
#endif        

 
// TO DO: Use Enzyme to calculate the potential derivatives
template <typename T> void cpuSingle(T *__restrict__ u, 
            T *__restrict__ xi, 
            T *__restrict__ qi, 
            int *__restrict__ ti, 
            int *__restrict__ ai, 
            T *__restrict__ mu, 
            T *__restrict__ eta, 
            int *__restrict__ kappa, 
            int dim, int ncq, int nmu, int neta, int nkappa, int inum, int potnum, int bondtype)
{
    if (bondtype==0)  {
        opuSinglea(u, xi, qi, ti, ai, mu, eta, kappa, dim, ncq, nmu, neta, nkappa, inum, potnum);    
    }    
    else if (bondtype==1) {
        opuSingleb(u, xi, qi, ti, ai, mu, eta, kappa, dim, ncq, nmu, neta, nkappa, inum, potnum);    
    }    
}
template void cpuSingle(double *, double *, double *, int *, int *, double *, double *,
        int*, int, int, int, int, int, int, int, int);
template void cpuSingle(float *, float *, float *, int *, int *, float *, float *, 
        int *, int, int, int, int, int, int, int, int);

template <typename T> void cpuComputeSingleEnergyForce(T *__restrict__ u,
            T *__restrict__ d_u,
            T *__restrict__ u_x,
            T *__restrict__ xi, 
            T *__restrict__ qi, 
            int *__restrict__ ti, 
            int *__restrict__ ai, 
            T *__restrict__ mu, 
            T *__restrict__ eta, 
            int *__restrict__ kappa, 
            int dim, int ncq, int nmu, int neta, int nkappa, int inum, int potnum, int bondtype)
{
    //T *d_u = new T[inum];
    cpuArraySetValue(d_u, (T) 1.0, inum);
    cpuArraySetValue(u_x, (T) 0.0, dim*inum);
 
#ifdef HAVE_ENZYME    
    __enzyme_autodiff((void*)cpuSingle<T>, 
                        enzyme_dup, u, d_u,
                        enzyme_dup, xi, u_x,
                        enzyme_const, qi,
                        enzyme_const, ti,
                        enzyme_const, ai,
                        enzyme_const, mu,
                        enzyme_const, eta,
                        enzyme_const, kappa,
                        dim, ncq, nmu, neta, nkappa, inum, potnum, bondtype);
#endif
}
template void cpuComputeSingleEnergyForce(double *, double *, double *, double *, double *, 
        int *, int *, double *, double *, int*, int, int, int, int, int, int, int, int);
template void cpuComputeSingleEnergyForce(float *, float *, float *, float *, float *,
        int *, int *, float *, float *, int *, int, int, int, int, int, int, int, int);

template <typename T> void cpuPair(T *__restrict__ u,
                                    T *__restrict__ xij,
                                    T *__restrict__ qi, 
                                    T *__restrict__ qj, 
                                    int *__restrict__ ti, 
                                    int *__restrict__ tj, 
                                    int *__restrict__ ai, 
                                    int *__restrict__ aj, 
                                    T *__restrict__ mu, 
                                    T *__restrict__ eta, 
                                    int *__restrict__ kappa, 
                                    int dim, int ncq, int nmu, int neta, 
                                    int nkappa, int ijnum, int potnum, int bondtype)
{
    if (bondtype==0) {
        opuPaira(u, xij, qi, qj, ti, tj, ai, aj, mu, eta, kappa, 
                dim, ncq, nmu, neta, nkappa, ijnum, potnum);    
    }    
    else if (bondtype==1) {
        opuPairb(u, xij, qi, qj, ti, tj, ai, aj, mu, eta, kappa, 
                dim, ncq, nmu, neta, nkappa, ijnum, potnum);    
    }    
    else if (bondtype==2) {
        opuPairc(u, xij, qi, qj, ti, tj, ai, aj, mu, eta, kappa, 
                dim, ncq, nmu, neta, nkappa, ijnum, potnum);    
    }    
    else if (bondtype==3) {
        opuTripletcPair(u, xij, qi, qj, ti, tj, ai, aj, mu, eta, kappa, 
                dim, ncq, nmu, neta, nkappa, ijnum, potnum);    
    }    
}
template void cpuPair(double *, double *, double *, double *, int *, int *, int *, int *, 
        double *, double *, int*, int, int, int, int, int, int, int, int);
template void cpuPair(float *, float *, float *, float *, int *, int *, int *, int *, 
        float *, float *, int *, int, int, int, int, int, int, int, int);

template <typename T> void cpuComputePairEnergyForce(T *__restrict__ u,
                                    T *__restrict__ d_u,
                                    T *__restrict__ u_x,
                                    T *__restrict__ xij,
                                    T *__restrict__ qi, 
                                    T *__restrict__ qj, 
                                    int *__restrict__ ti, 
                                    int *__restrict__ tj, 
                                    int *__restrict__ ai, 
                                    int *__restrict__ aj, 
                                    T *__restrict__ mu, 
                                    T *__restrict__ eta, 
                                    int *__restrict__ kappa, 
                                    int dim, int ncq, int nmu, int neta, 
                                    int nkappa, int ijnum, int potnum, int bondtype)
{
    //T *d_u = new T[ijnum];
    cpuArraySetValue(d_u, (T) 1.0, ijnum);
    cpuArraySetValue(u_x, (T) 0.0, dim*ijnum);

#ifdef HAVE_ENZYME    
    __enzyme_autodiff((void*)cpuPair<T>, 
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
                        dim, ncq, nmu, neta, nkappa, ijnum, potnum, bondtype);  
#endif        
}
template void cpuComputePairEnergyForce(double *, double *, double *, double *, double *, double *, 
        int *, int *, int *, int *, double *, double *, int*, int, int, int, int, int, int, int, int);
template void cpuComputePairEnergyForce(float *, float *, float *, float *, float *, float *,
        int *, int *, int *, int *, float *, float *, int *, int, int, int, int, int, int, int, int);

template <typename T> void cpuTriplet(T *__restrict__ u,
                                    T *__restrict__ xij,
                                    T *__restrict__ xik,
                                    T *__restrict__ qi, 
                                    T *__restrict__ qj, 
                                    T *__restrict__ qk,         
                                    int *__restrict__ ti, 
                                    int *__restrict__ tj, 
                                    int *__restrict__ tk,         
                                    int *__restrict__ ai, 
                                    int *__restrict__ aj, 
                                    int *__restrict__ ak,         
                                    T *__restrict__ mu, 
                                    T *__restrict__ eta, 
                                    int *__restrict__ kappa, 
                                    int dim, int ncq, int nmu, int neta, 
                                    int nkappa, int ijknum, int potnum, int bondtype)
{
    if (bondtype==0) {
        opuTripleta(u, xij, xik, qi, qj, qk, ti, tj, tk, ai, aj, ak, mu, eta, 
                kappa, dim, ncq, nmu, neta, nkappa, ijknum, potnum);        
    }    
    else if (bondtype==1) {
        opuTripletb(u, xij, xik, qi, qj, qk, ti, tj, tk, ai, aj, ak, mu, eta, 
                kappa, dim, ncq, nmu, neta, nkappa, ijknum, potnum);
    }    
    else if (bondtype==2) {
        opuTripletc(u, xij, xik, qi, qj, qk, ti, tj, tk, ai, aj, ak, mu, eta, 
                kappa, dim, ncq, nmu, neta, nkappa, ijknum, potnum);
    }    
}
template void cpuTriplet(double *, double *, double *, double *, double *, double *, int *, int *, 
        int *, int *, int *, int *, double *, double *, int*, int, int, int, int, int, int, int, int);
template void cpuTriplet(float *, float *, float *, float *, float *, float *, int *, int *, 
        int *, int *, int *, int *, float *, float *, int *, int, int, int, int, int, int, int, int);

template <typename T> void cpuComputeTripletEnergyForce(T *__restrict__ u,
                                    T *__restrict__ d_u,        
                                    T *__restrict__ u_xij,   
                                    T *__restrict__ u_xik,   
                                    T *__restrict__ xij,
                                    T *__restrict__ xik,
                                    T *__restrict__ qi, 
                                    T *__restrict__ qj, 
                                    T *__restrict__ qk,         
                                    int *__restrict__ ti, 
                                    int *__restrict__ tj, 
                                    int *__restrict__ tk,         
                                    int *__restrict__ ai, 
                                    int *__restrict__ aj, 
                                    int *__restrict__ ak,         
                                    T *__restrict__ mu, 
                                    T *__restrict__ eta, 
                                    int *__restrict__ kappa, 
                                    int dim, int ncq, int nmu, int neta, 
                                    int nkappa, int ijknum, int potnum, int bondtype)
{    
    //T *d_u = new T[ijknum];
    cpuArraySetValue(d_u, (T) 1.0, ijknum);
    cpuArraySetValue(u_xij, (T) 0.0, dim*ijknum);
    cpuArraySetValue(u_xik, (T) 0.0, dim*ijknum);

#ifdef HAVE_ENZYME    
    __enzyme_autodiff((void*)cpuTriplet<T>, 
                        enzyme_dup, u, d_u,
                        enzyme_dup, xij, u_xij,
                        enzyme_dup, xik, u_xik,
                        enzyme_const, qi,
                        enzyme_const, qj,
                        enzyme_const, qk,            
                        enzyme_const, ti,
                        enzyme_const, tj,
                        enzyme_const, tk,            
                        enzyme_const, ai,
                        enzyme_const, aj,
                        enzyme_const, ak,            
                        // TODO can make this non const to get du/dmu
                        enzyme_const, mu,
                        enzyme_const, eta,
                        enzyme_const, kappa,
                        dim, ncq, nmu, neta, nkappa, ijknum, potnum, bondtype);  
#endif    
}
template void cpuComputeTripletEnergyForce(double *, double *, double *, double *, double *, double *, double *, 
        double *, double *, int *, int *, int *, int *, int *, int *, double *, double *, int*, int, int, 
        int, int, int, int, int, int);
template void cpuComputeTripletEnergyForce(float *, float *, float *, float *, float *, float *, float *, 
        float *, float *, int *, int *, int *, int *, int *, int *, float *, float *, int *, int, int, 
        int, int, int, int, int, int);

template <typename T> void cpuQuadruplet(T *__restrict__ u,
                                    T *__restrict__ xij,
                                    T *__restrict__ xik,
                                    T *__restrict__ xil,        
                                    T *__restrict__ qi, 
                                    T *__restrict__ qj, 
                                    T *__restrict__ qk,
                                    T *__restrict__ ql,        
                                    int *__restrict__ ti, 
                                    int *__restrict__ tj, 
                                    int *__restrict__ tk,    
                                    int *__restrict__ tl,            
                                    int *__restrict__ ai, 
                                    int *__restrict__ aj, 
                                    int *__restrict__ ak,    
                                    int *__restrict__ al,            
                                    T *__restrict__ mu, 
                                    T *__restrict__ eta, 
                                    int *__restrict__ kappa, 
                                    int dim, int ncq, int nmu, int neta, 
                                    int nkappa, int ijklnum, int potnum, int bondtype)
{
    if (bondtype==0) {
        opuQuadrupleta(u, xij, xik, xil, qi, qj, qk, ql, ti, tj, tk, tl, ai, aj, ak, al,
                mu, eta, kappa, dim, ncq, nmu, neta, nkappa, ijklnum, potnum);        
    }    
    else if (bondtype==1) {
        opuQuadrupletb(u, xij, xik, xil, qi, qj, qk, ql, ti, tj, tk, tl, ai, aj, ak, al, 
                mu, eta, kappa, dim, ncq, nmu, neta, nkappa, ijklnum, potnum);
    }            
}
template void cpuQuadruplet(double *, double *, double *, double *, double *, double *, double *, 
        double *, int *, int *, int *, int *, int *, int *, int *, int *, double *, double *, 
        int*, int, int, int, int, int, int, int, int);
template void cpuQuadruplet(float *, float *, float *, float *, float *, float *, float *, 
        float *, int *, int *, int *, int *, int *, int *, int *, int *, float *, float *, 
        int *, int, int, int, int, int, int, int, int);

template <typename T> void cpuComputeQuadrupletEnergyForce(T *__restrict__ u,
                                    T *__restrict__ d_u,             
                                    T *__restrict__ u_xij,   
                                    T *__restrict__ u_xik,   
                                    T *__restrict__ u_xil,   
                                    T *__restrict__ xij,
                                    T *__restrict__ xik,
                                    T *__restrict__ xil,        
                                    T *__restrict__ qi, 
                                    T *__restrict__ qj, 
                                    T *__restrict__ qk,
                                    T *__restrict__ ql,        
                                    int *__restrict__ ti, 
                                    int *__restrict__ tj, 
                                    int *__restrict__ tk,    
                                    int *__restrict__ tl,            
                                    int *__restrict__ ai, 
                                    int *__restrict__ aj, 
                                    int *__restrict__ ak,    
                                    int *__restrict__ al,            
                                    T *__restrict__ mu, 
                                    T *__restrict__ eta, 
                                    int *__restrict__ kappa, 
                                    int dim, int ncq, int nmu, int neta, 
                                    int nkappa, int ijklnum, int potnum, int bondtype)
{
    //T *d_u = new T[ijklnum];
    cpuArraySetValue(d_u, (T) 1.0, ijklnum);
    cpuArraySetValue(u_xij, (T) 0.0, dim*ijklnum);
    cpuArraySetValue(u_xik, (T) 0.0, dim*ijklnum);
    cpuArraySetValue(u_xil, (T) 0.0, dim*ijklnum);
    
#ifdef HAVE_ENZYME    
    __enzyme_autodiff((void*)cpuQuadruplet<T>, 
                        enzyme_dup, u, d_u,
                        enzyme_dup, xij, u_xij,
                        enzyme_dup, xik, u_xik,
                        enzyme_dup, xil, u_xil,
                        enzyme_const, qi,
                        enzyme_const, qj,
                        enzyme_const, qk,
                        enzyme_const, ql,            
                        enzyme_const, ti,
                        enzyme_const, tj,
                        enzyme_const, tk,      
                        enzyme_const, tl,                  
                        enzyme_const, ai,
                        enzyme_const, aj,
                        enzyme_const, ak,        
                        enzyme_const, al,                    
                        // TODO can make this non const to get du/dmu
                        enzyme_const, mu,
                        enzyme_const, eta,
                        enzyme_const, kappa,
                        dim, ncq, nmu, neta, nkappa, ijklnum, potnum, bondtype);  
#endif    
}
template void cpuComputeQuadrupletEnergyForce(double *, double *, double *, double *, double *, double *, 
        double *, double *, double *, double *, double *, double *, int *, int *, int *, int *, int *, 
        int *, int *, int *, double *, double *, int*, int, int, int, int, int, int, int, int);
template void cpuComputeQuadrupletEnergyForce(float *, float *, float *, float *, float *, float *, 
        float *, float *, float *, float *, float *, float *, int *, int *, int *, int *, int *, 
        int *, int *, int *, float *, float *, int *, int, int, int, int, int, int, int, int);


template <typename T> void cpuPaircDensityGradient(T *__restrict__ u,
                                    T *__restrict__ d_u,                     
                                    T *__restrict__ u_rho,   
                                    T *__restrict__ rho,   
                                    T *__restrict__ mu, 
                                    T *__restrict__ eta, 
                                    int *__restrict__ kappa, 
                                    int nrho, int nmu, int neta, 
                                    int nkappa, int inum, int potnum)
{
    //T *d_u = new T[inum];
    cpuArraySetValue(d_u, (T) 1.0, inum);
    cpuArraySetValue(u_rho, (T) 0.0, inum);    
    
#ifdef HAVE_ENZYME    
    __enzyme_autodiff((void*)opuPaircDensity<T>, 
                        enzyme_dup, u, d_u,
                        enzyme_dup, rho, u_rho,                        
                        enzyme_const, mu,
                        enzyme_const, eta,
                        enzyme_const, kappa,
                        1, nmu, neta, nkappa, inum, potnum);  
#endif
}
template void cpuPaircDensityGradient(double *, double *, double *, double *, double *, double *,
        int*, int, int, int, int, int, int);
template void cpuPaircDensityGradient(float *, float *, float *, float *, float *, float *,
        int *, int, int, int, int, int, int);
        
template <typename T> void cpuTripletcDensityGradient(T *__restrict__ u,
                                    T *__restrict__ d_u,                             
                                    T *__restrict__ u_rho,   
                                    T *__restrict__ rho,   
                                    T *__restrict__ mu, 
                                    T *__restrict__ eta, 
                                    int *__restrict__ kappa, 
                                    int nrho, int nmu, int neta, 
                                    int nkappa, int inum, int potnum)
{
    //T *d_u = new T[inum];
    cpuArraySetValue(d_u, (T) 1.0, inum);
    cpuArraySetValue(u_rho, (T) 0.0, inum);    
    
#ifdef HAVE_ENZYME    
    __enzyme_autodiff((void*)opuTripletcDensity<T>, 
                        enzyme_dup, u, d_u,
                        enzyme_dup, rho, u_rho,                        
                        enzyme_const, mu,
                        enzyme_const, eta,
                        enzyme_const, kappa,
                        1, nmu, neta, nkappa, inum, potnum);  
#endif    
}
template void cpuTripletcDensityGradient(double *, double *, double *, double *, double *, double *,
        int*, int, int, int, int, int, int);
template void cpuTripletcDensityGradient(float *, float *, float *, float *, float *, float *,
        int *, int, int, int, int, int, int);

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

void cpuEmbedingForce(dstype *fij, dstype *d_rhoi, int *pairnum, int *pairnumsum, int inum)
{    
    for (int i=0; i<inum; i++) {
        int jnum = pairnum[i];
        int start = pairnumsum[i];
        for (int j=0; j<jnum; j++)             
            fij[start+j] = d_rhoi[i]*fij[start+j];                
    }        
}

// template <typename T> void primal_opuLJ1pot(T *__restrict__ u,
//                                     T *__restrict__ xij,
//                                     T *__restrict__ qi, 
//                                     T *__restrict__ qj, 
//                                     int *__restrict__ ti, 
//                                     int *__restrict__ tj, 
//                                     int *__restrict__ ai, 
//                                     int *__restrict__ aj, 
//                                     T *__restrict__ mu, 
//                                     T *__restrict__ eta, 
//                                     int *__restrict__ kappa, 
//                                     int dim, int ncq, int nmu, int neta, int nkappa, int ng, int potnum)
// {
// 	for (int i = 0; i <ng; i++) {
// 		T mu1 = mu[0]; // mu_1
// 		T mu2 = mu[1]; // mu_2
// 		T xij1 = xij[0 + i*3]; // x
// 		T xij2 = xij[1 + i*3]; // y
// 		T xij3 = xij[2 + i*3]; // z
// 		T t2 = xij1*xij1;
// 		T t3 = xij2*xij2;
// 		T t4 = xij3*xij3;
// 		T t5 = t2+t3+t4; // rij^2
//         T t6 = t5*t5*t5; // rij^6
//         //T t8 = t6*t5;    // rij^8
//         T t12 = t6*t6;   // rij^12
//         //T t14 = t12*t5;  // rij^14
// 		u[i] = mu1*1.0/(t12) - mu2*1.0/(t6);        
// 	}	
// }
// 
// #if 1
// 
// template <typename T> void opuLJ1pot(T *__restrict__ u,
//                                     T *__restrict__ d_xij,
//                                     T *__restrict__ xij,
//                                     T *__restrict__ qi, 
//                                     T *__restrict__ qj, 
//                                     int *__restrict__ ti, 
//                                     int *__restrict__ tj, 
//                                     int *__restrict__ ai, 
//                                     int *__restrict__ aj, 
//                                     T *__restrict__ mu, 
//                                     T *__restrict__ eta, 
//                                     int *__restrict__ kappa, 
//                                     int dim, int ncq, int nmu, int neta, int nkappa, int ng, int potnum) {
//     T *d_u = new T[ng];
//     // TODO double check
//     for (int i=0; i<ng; ++i) {
//         d_u[i] = 1.0;
//     }
// 
//     // We assume that the mem is already allocated and zero'd
//     //d_xij = new T[3*ng]
//     for (int i=0; i<3*ng; ++i) {
//         d_xij[i] = 0.0;
//     }
// 
//     // The following does
//     // d_xij[..] += d(u[..])/d(xij[...]) * d_u[...]
//     // d_u[...] = 0
// 
// 
//     __enzyme_autodiff((void*)primal_opuLJ1pot<T>, 
//                         enzyme_dup, u, d_u,
//                         enzyme_dup, xij, d_xij,
//                         enzyme_const, qi,
//                         enzyme_const, qj,
//                         enzyme_const, ti,
//                         enzyme_const, tj,
//                         enzyme_const, ai,
//                         enzyme_const, aj,
//                         // TODO can make this non const to get du/dmu
//                         enzyme_const, mu,
//                         enzyme_const, eta,
//                         enzyme_const, kappa,
//                         dim, ncq, nmu, neta, nkappa, ng, potnum);
//     delete[] d_u;
// }
// #else
// template <typename T> void opuLJ1pot(T *__restrict__ u, T *__restrict__ f, T *__restrict__ xij, T *__restrict__ qi, T *__restrict__ qj, int *__restrict__ ti, int *__restrict__ tj, int *__restrict__ ai, int *__restrict__ aj, T *__restrict__ mu, T *__restrict__ eta, int *__restrict__ kappa, int dim, int ncq, int nmu, int neta, int nkappa, int ng, int potnum)
// {
// 	for (int i = 0; i <ng; i++) {
// 		T mu1 = mu[0]; // mu_1
// 		T mu2 = mu[1]; // mu_2
// 		T xij1 = xij[0 + i*3]; // x
// 		T xij2 = xij[1 + i*3]; // y
// 		T xij3 = xij[2 + i*3]; // z
// 		T t2 = xij1*xij1;
// 		T t3 = xij2*xij2;
// 		T t4 = xij3*xij3;
// 		T t5 = t2+t3+t4; // rij^2
//         T t6 = t5*t5*t5; // rij^6
//         T t8 = t6*t5;    // rij^8
//         T t12 = t6*t6;   // rij^12
//         T t14 = t12*t5;  // rij^14
// 		u[i] = mu1*1.0/(t12) - mu2*1.0/(t6);
//         f[0 + 3*i] = (12*mu1*xij1)/t14 - (6*mu2*xij1)/t8;
//         f[1 + 3*i] = (12*mu1*xij2)/t14 - (6*mu2*xij2)/t8;
//         f[2 + 3*i] = (12*mu1*xij3)/t14 - (6*mu2*xij3)/t8;
// 	}	
// }
// #endif
// template void opuLJ1pot(double *, double *, double *, double *, double *, int *, int *, int *, int *, double *, double *, int*, int, int, int, int, int, int, int);
// template void opuLJ1pot(float *, float *, float *, float *, float *, int *, int *, int *, int *, float *, float *, int *, int, int, int, int, int, int, int);
// 
// template <typename T> void opuLJ2pot(T *u, T *f, T *xij, T *qi, T *qj, int *ti, int *tj, int *ai, int *aj, T *mu, T *eta, int *kappa, int dim, int ncq, int nmu, int neta, int nkappa, int ng, int potnum)
// {
// 	for (int i = 0; i <ng; i++) {                
//         int ntype = kappa[0]; // number of atom types
//         int a = ti[i]-1 + ntype*(tj[i]-1);  // ti-tj interaction      
// 		T mu1 = mu[0 + 2*a];
// 		T mu2 = mu[1 + 2*a];
// 		T xij1 = xij[0 + i*3];
// 		T xij2 = xij[1 + i*3];
// 		T xij3 = xij[2 + i*3];
// 		T t2 = xij1*xij1;
// 		T t3 = xij2*xij2;
// 		T t4 = xij3*xij3;
// 		T t5 = t2+t3+t4; // rij^2
//         T t6 = t5*t5*t5; // rij^6
//         T t8 = t6*t5;    // rij^8
//         T t12 = t6*t6;   // rij^12
//         T t14 = t12*t5;  // rij^14
// 		u[i] = mu1*1.0/(t12) - mu2*1.0/(t6);
//         f[0 + 3*i] = (12*mu1*xij1)/t14 - (6*mu2*xij1)/t8;
//         f[1 + 3*i] = (12*mu1*xij2)/t14 - (6*mu2*xij2)/t8;
//         f[2 + 3*i] = (12*mu1*xij3)/t14 - (6*mu2*xij3)/t8;                
// 	}	
// }
// template void opuLJ2pot(double *, double *, double *, double *, double *, int *, int *, int *, int *, double *, double *, int*, int, int, int, int, int, int, int);
// template void opuLJ2pot(float *, float *, float *, float *, float *, int *, int *, int *, int *, float *, float *, int *, int, int, int, int, int, int, int);
// 
// void cpuComputePair1EnergyForce(dstype *eij, dstype *fij, dstype *xij, dstype *qi, dstype *qj, 
//         int *ti, int *tj, int *ai, int *aj, dstype *mu, dstype *eta, int *kappa, int dim, 
//         int ncq, int nmu, int neta, int nkappa, int ijnum, int potnum, int bondtype)
// {
//     if (bondtype==0) {
//         //opuPaira(eij, xij, qi, qj, ti, tj, ai, aj, mu, eta, kappa, dim, ncq, nmu, neta, nkappa, ijnum, potnum);    
//         opuLJ1pot(eij, fij, xij, qi, qj, ti, tj, ai, aj, mu, eta, kappa, dim, ncq, nmu, neta, nkappa, ijnum, potnum);    
//     }
//     
//     if (bondtype==1) {
//         opuPairb(eij, xij, qi, qj, ti, tj, ai, aj, mu, eta, kappa, dim, ncq, nmu, neta, nkappa, ijnum, potnum);    
//     }
//     
//     if (bondtype==2) {
//         opuPairc(eij, xij, qi, qj, ti, tj, ai, aj, mu, eta, kappa, dim, ncq, nmu, neta, nkappa, ijnum, potnum);    
//     }
//     
//     if (bondtype==3) {
//         opuTripletcPair(eij, xij, qi, qj, ti, tj, ai, aj, mu, eta, kappa, dim, ncq, nmu, neta, nkappa, ijnum, potnum);    
//     }
// }

// void cpuElectronDensity(dstype *rhoi, dstype *rhoij, int *pairnum, int *pairnumsum, int inum) 
// {
//     for (int i=0; i<inum; i++) {
//         int jnum = pairnum[i];
//         int start = pairnumsum[i];
//         rhoi[i] = 0.0;
//         for (int j=0; j<jnum; j++) 
//             rhoi[i] += rhoij[start+j];        
//     }
// }
// void cpuPaircDensity(dstype *u, dstype *urho, dstype *rho, dstype *mu, dstype *eta, int *kappa, 
//         int nrho, int nmu, int neta, int nkappa, int ng, int potnum)
// {
//     opuPaircDensity(u, rho, mu, eta, kappa, 1, nmu, neta, nkappa, ng, potnum);
// }
// void cpuEmbedingForce(dstype *fij, dstype *drhoi, int *pairnum, int *pairnumsum, int inum)
// {    
//     for (int i=0; i<inum; i++) {
//         int jnum = pairnum[i];
//         int start = pairnumsum[i];
//         for (int j=0; j<jnum; j++)             
//             fij[start+j] = drhoi[i]*fij[start+j];                
//     }        
// }
        
// void cpuComputeTripletEnergyForce(dstype *eijk, dstype *fij, dstype *fik, dstype *xij, dstype *xik, dstype *qi, dstype *qj, dstype *qk, 
//      int *ti, int *tj, int *tk, int *ai, int *aj, int *ak, dstype *mu, dstype *eta, int *kappa, int dim, int ncq, 
//      int nmu, int neta, int nkappa, int ijknum, int potnum, int bondtype)        
// {
//     if (bondtype==0) {
//         opuTripleta(eijk, xij, xik, qi, qj, qk, ti, tj, tk, ai, aj, ak, mu, eta, kappa, dim, ncq, nmu, neta, nkappa, ijknum, potnum);        
//     }
//     
//     if (bondtype==1) {
//         opuTripletb(eijk, xij, xik, qi, qj, qk, ti, tj, tk, ai, aj, ak, mu, eta, kappa, dim, ncq, nmu, neta, nkappa, ijknum, potnum);
//     }
//     
//     if (bondtype==2) {
//         opuTripletc(eijk, xij, xik, qi, qj, qk, ti, tj, tk, ai, aj, ak, mu, eta, kappa, dim, ncq, nmu, neta, nkappa, ijknum, potnum);
//     }
// }
// 
// // void cpuTripletcDensity(dstype *u, dstype *urho, dstype *rho, dstype *mu, dstype *eta, int *kappa, 
// //         int nrho, int nmu, int neta, int nkappa, int ng, int potnum)
// // {
// //     opuTripletcDensity(u, rho, mu, eta, kappa, 1, nmu, neta, nkappa, ng, potnum);
// // }
// 
// void cpuComputeQuadrupletEnergyForce(dstype *eijkl, dstype *fij, dstype *fik, dstype *fil, dstype *xij, dstype *xik, dstype *xil, 
//       dstype *qi, dstype *qj, dstype *qk, dstype *ql, int *ti, int *tj, int *tk, int *tl, int *ai, int *aj, int *ak, int *al,
//       dstype *mu, dstype *eta, int *kappa, int dim, int ncq, int nmu, int neta, int nkappa, int ijklnum, int potnum, int bondtype)        
// {
//     if (bondtype==0) {
//         opuQuadrupleta(eijkl, xij, xik, xil, qi, qj, qk, ql, ti, tj, tk, tl, ai, aj, ak, al,
//                 mu, eta, kappa, dim, ncq, nmu, neta, nkappa, ijklnum, potnum);        
//     }
//     
//     if (bondtype==1) {
//         opuQuadrupletb(eijkl, xij, xik, xil, qi, qj, qk, ql, ti, tj, tk, tl, ai, aj, ak, al, 
//                 mu, eta, kappa, dim, ncq, nmu, neta, nkappa, ijklnum, potnum);
//     }        
// }
// 

#endif


