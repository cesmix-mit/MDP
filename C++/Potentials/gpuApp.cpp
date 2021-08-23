/***************************************************************************                               
                    Molecular Dynamics Potentials (MDP)
                           CESMIX-MIT Project  
 
 Contributing authors: Ngoc-Cuong Nguyen (cuongng@mit.edu, exapde@gmail.com)
 ***************************************************************************/

#ifndef __GPUAPP
#define __GPUAPP

//#include <math.h>


// #ifdef HAVE_ENZYME                
// template <typename... Args>
// void __enzyme_autodiff(void*, Args... args);
// int enzyme_const, enzyme_dup;
// #endif        

// TO DO: Use Enzyme to calculate the potential derivatives
template <typename T> void gpuSingle(T *__restrict__ u, 
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
        gpuSinglea(u, xi, qi, ti, ai, mu, eta, kappa, dim, ncq, nmu, neta, nkappa, inum, potnum);    
    }    
    else if (bondtype==1) {
        gpuSingleb(u, xi, qi, ti, ai, mu, eta, kappa, dim, ncq, nmu, neta, nkappa, inum, potnum);    
    }    
}
template void gpuSingle(double *, double *, double *, int *, int *, double *, double *,
        int*, int, int, int, int, int, int, int, int);
template void gpuSingle(float *, float *, float *, int *, int *, float *, float *, 
        int *, int, int, int, int, int, int, int, int);

template <typename T> void gpuComputeSingleEnergyForce(T *__restrict__ u,
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
    gpuArraySetValue(u_x, (T) 0.0, dim*inum);
 
#ifdef HAVE_ENZYME        
    gpuArraySetValue(d_u, (T) 1.0, inum);
    if (bondtype==0)  {
        gpuSingleaGradient(u, d_u, u_x, xi, qi, ti, ai, mu, eta, kappa, dim, ncq, nmu, neta, nkappa, inum, potnum);    
    }    
    else if (bondtype==1) {
        gpuSinglebGradient(u, d_u, u_x, xi, qi, ti, ai, mu, eta, kappa, dim, ncq, nmu, neta, nkappa, inum, potnum);    
    }    
#else
    if (bondtype==0)  {
        gpuSingleaGradient(u, u_x, xi, qi, ti, ai, mu, eta, kappa, dim, ncq, nmu, neta, nkappa, inum, potnum);    
    }    
    else if (bondtype==1) {
        gpuSinglebGradient(u, u_x, xi, qi, ti, ai, mu, eta, kappa, dim, ncq, nmu, neta, nkappa, inum, potnum);    
    }        
#endif
    
//     //T *d_u = new T[inum];
//     gpuArraySetValue(d_u, (T) 1.0, inum);
//     gpuArraySetValue(u_x, (T) 0.0, dim*inum);
//  
// #ifdef HAVE_ENZYME    
//     __enzyme_autodiff((void*)gpuSingle<T>, 
//                         enzyme_dup, u, d_u,
//                         enzyme_dup, xi, u_x,
//                         enzyme_const, qi,
//                         enzyme_const, ti,
//                         enzyme_const, ai,
//                         enzyme_const, mu,
//                         enzyme_const, eta,
//                         enzyme_const, kappa,
//                         dim, ncq, nmu, neta, nkappa, inum, potnum, bondtype);
// #endif
}
template void gpuComputeSingleEnergyForce(double *, double *, double *, double *, double *, 
        int *, int *, double *, double *, int*, int, int, int, int, int, int, int, int);
template void gpuComputeSingleEnergyForce(float *, float *, float *, float *, float *,
        int *, int *, float *, float *, int *, int, int, int, int, int, int, int, int);

template <typename T> void gpuPair(T *__restrict__ u,
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
        gpuPaira(u, xij, qi, qj, ti, tj, ai, aj, mu, eta, kappa, 
                dim, ncq, nmu, neta, nkappa, ijnum, potnum);    
    }    
    else if (bondtype==1) {
        gpuPairb(u, xij, qi, qj, ti, tj, ai, aj, mu, eta, kappa, 
                dim, ncq, nmu, neta, nkappa, ijnum, potnum);    
    }    
    else if (bondtype==2) {
        gpuPairc(u, xij, qi, qj, ti, tj, ai, aj, mu, eta, kappa, 
                dim, ncq, nmu, neta, nkappa, ijnum, potnum);    
    }    
    else if (bondtype==3) {
        gpuTripletcPair(u, xij, qi, qj, ti, tj, ai, aj, mu, eta, kappa, 
                dim, ncq, nmu, neta, nkappa, ijnum, potnum);    
    }    
}
template void gpuPair(double *, double *, double *, double *, int *, int *, int *, int *, 
        double *, double *, int*, int, int, int, int, int, int, int, int);
template void gpuPair(float *, float *, float *, float *, int *, int *, int *, int *, 
        float *, float *, int *, int, int, int, int, int, int, int, int);

template <typename T> void gpuComputePairEnergyForce(T *__restrict__ u,
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
    //printf("called gpuComputePairEnergyForce\n");    
    gpuArraySetValue(u_x, (T) 0.0, dim*ijnum);

    //gpuLJ(u, u_x, xij, qi, qj, ti, tj, ai, aj, mu, eta, kappa, 
    //        dim, ncq, nmu, neta, nkappa, ijnum);  
//     gpuGradientLJ(u, d_u, xij, u_x, qi, qj, ti, tj, ai, aj, mu, eta, kappa, 
//             dim, ncq, nmu, neta, nkappa, ijnum);
    //printf("Called gpuGradient u_x[0]=%f d_u[0]=%f\n", u_x[0], d_u[0]);
    
#ifdef HAVE_ENZYME        
    gpuArraySetValue(d_u, (T) 1.0, ijnum);
    if (bondtype==0) {
        gpuPairaGradient(u, d_u, u_x, xij, qi, qj, ti, tj, ai, aj, mu, eta, kappa, 
                dim, ncq, nmu, neta, nkappa, ijnum, potnum);    
    }    
    else if (bondtype==1) {
        gpuPairbGradient(u, d_u, u_x, xij, qi, qj, ti, tj, ai, aj, mu, eta, kappa, 
                dim, ncq, nmu, neta, nkappa, ijnum, potnum);    
    }    
    else if (bondtype==2) {
        gpuPaircGradient(u, d_u, u_x, xij, qi, qj, ti, tj, ai, aj, mu, eta, kappa, 
                dim, ncq, nmu, neta, nkappa, ijnum, potnum);     
    }    
    else if (bondtype==3) {
        gpuTripletcPairGradient(u, d_u, u_x, xij, qi, qj, ti, tj, ai, aj, mu, eta, kappa, 
                dim, ncq, nmu, neta, nkappa, ijnum, potnum);     
    }    
#else
    if (bondtype==0) {
        gpuPairaGradient(u, u_x, xij, qi, qj, ti, tj, ai, aj, mu, eta, kappa, 
                dim, ncq, nmu, neta, nkappa, ijnum, potnum);    
    }    
    else if (bondtype==1) {
        gpuPairbGradient(u, u_x, xij, qi, qj, ti, tj, ai, aj, mu, eta, kappa, 
                dim, ncq, nmu, neta, nkappa, ijnum, potnum);    
    }    
    else if (bondtype==2) {
        gpuPaircGradient(u, u_x, xij, qi, qj, ti, tj, ai, aj, mu, eta, kappa, 
                dim, ncq, nmu, neta, nkappa, ijnum, potnum);     
    }    
    else if (bondtype==3) {
        gpuTripletcPairGradient(u, u_x, xij, qi, qj, ti, tj, ai, aj, mu, eta, kappa, 
                dim, ncq, nmu, neta, nkappa, ijnum, potnum);     
    }        
#endif            
// #ifdef HAVE_ENZYME    
//     __enzyme_autodiff((void*)gpuPair<T>, 
//                         enzyme_dup, u, d_u,
//                         enzyme_dup, xij, u_x,
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
//                         dim, ncq, nmu, neta, nkappa, ijnum, potnum, bondtype);  
// #endif        
}
template void gpuComputePairEnergyForce(double *, double *, double *, double *, double *, double *, 
        int *, int *, int *, int *, double *, double *, int*, int, int, int, int, int, int, int, int);
template void gpuComputePairEnergyForce(float *, float *, float *, float *, float *, float *,
        int *, int *, int *, int *, float *, float *, int *, int, int, int, int, int, int, int, int);

template <typename T> void gpuTriplet(T *__restrict__ u,
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
        gpuTripleta(u, xij, xik, qi, qj, qk, ti, tj, tk, ai, aj, ak, mu, eta, 
                kappa, dim, ncq, nmu, neta, nkappa, ijknum, potnum);        
    }    
    else if (bondtype==1) {
        gpuTripletb(u, xij, xik, qi, qj, qk, ti, tj, tk, ai, aj, ak, mu, eta, 
                kappa, dim, ncq, nmu, neta, nkappa, ijknum, potnum);
    }    
    else if (bondtype==2) {
        gpuTripletc(u, xij, xik, qi, qj, qk, ti, tj, tk, ai, aj, ak, mu, eta, 
                kappa, dim, ncq, nmu, neta, nkappa, ijknum, potnum);
    }    
}
template void gpuTriplet(double *, double *, double *, double *, double *, double *, int *, int *, 
        int *, int *, int *, int *, double *, double *, int*, int, int, int, int, int, int, int, int);
template void gpuTriplet(float *, float *, float *, float *, float *, float *, int *, int *, 
        int *, int *, int *, int *, float *, float *, int *, int, int, int, int, int, int, int, int);

template <typename T> void gpuComputeTripletEnergyForce(T *__restrict__ u,
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
    gpuArraySetValue(u_xij, (T) 0.0, dim*ijknum);
    gpuArraySetValue(u_xik, (T) 0.0, dim*ijknum);
    
#ifdef HAVE_ENZYME        
    gpuArraySetValue(d_u, (T) 1.0, ijknum);
    if (bondtype==0) {
        gpuTripletaGradient(u, d_u, u_xij, u_xik, xij, xik, qi, qj, qk, ti, tj, tk, ai, aj, ak, mu, eta, 
                kappa, dim, ncq, nmu, neta, nkappa, ijknum, potnum);        
    }    
    else if (bondtype==1) {
        gpuTripletbGradient(u, d_u, u_xij, u_xik, xij, xik, qi, qj, qk, ti, tj, tk, ai, aj, ak, mu, eta, 
                kappa, dim, ncq, nmu, neta, nkappa, ijknum, potnum);
    }    
    else if (bondtype==2) {
        gpuTripletcGradient(u, d_u, u_xij, u_xik, xij, xik, qi, qj, qk, ti, tj, tk, ai, aj, ak, mu, eta, 
                kappa, dim, ncq, nmu, neta, nkappa, ijknum, potnum);
    }    
#else
    if (bondtype==0) {
        gpuTripletaGradient(u, u_xij, u_xik, xij, xik, qi, qj, qk, ti, tj, tk, ai, aj, ak, mu, eta, 
                kappa, dim, ncq, nmu, neta, nkappa, ijknum, potnum);        
    }    
    else if (bondtype==1) {
        gpuTripletbGradient(u, u_xij, u_xik, xij, xik, qi, qj, qk, ti, tj, tk, ai, aj, ak, mu, eta, 
                kappa, dim, ncq, nmu, neta, nkappa, ijknum, potnum);
    }    
    else if (bondtype==2) {
        gpuTripletcGradient(u, u_xij, u_xik, xij, xik, qi, qj, qk, ti, tj, tk, ai, aj, ak, mu, eta, 
                kappa, dim, ncq, nmu, neta, nkappa, ijknum, potnum);
    }        
#endif        
//     //T *d_u = new T[ijknum];
//     gpuArraySetValue(d_u, (T) 1.0, ijknum);
//     gpuArraySetValue(u_xij, (T) 0.0, dim*ijknum);
//     gpuArraySetValue(u_xik, (T) 0.0, dim*ijknum);
// 
// #ifdef HAVE_ENZYME    
//     __enzyme_autodiff((void*)gpuTriplet<T>, 
//                         enzyme_dup, u, d_u,
//                         enzyme_dup, xij, u_xij,
//                         enzyme_dup, xik, u_xik,
//                         enzyme_const, qi,
//                         enzyme_const, qj,
//                         enzyme_const, qk,            
//                         enzyme_const, ti,
//                         enzyme_const, tj,
//                         enzyme_const, tk,            
//                         enzyme_const, ai,
//                         enzyme_const, aj,
//                         enzyme_const, ak,            
//                         // TODO can make this non const to get du/dmu
//                         enzyme_const, mu,
//                         enzyme_const, eta,
//                         enzyme_const, kappa,
//                         dim, ncq, nmu, neta, nkappa, ijknum, potnum, bondtype);  
// #endif    
}
template void gpuComputeTripletEnergyForce(double *, double *, double *, double *, double *, double *, double *, 
        double *, double *, int *, int *, int *, int *, int *, int *, double *, double *, int*, int, int, 
        int, int, int, int, int, int);
template void gpuComputeTripletEnergyForce(float *, float *, float *, float *, float *, float *, float *, 
        float *, float *, int *, int *, int *, int *, int *, int *, float *, float *, int *, int, int, 
        int, int, int, int, int, int);

template <typename T> void gpuQuadruplet(T *__restrict__ u,
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
        gpuQuadrupleta(u, xij, xik, xil, qi, qj, qk, ql, ti, tj, tk, tl, ai, aj, ak, al,
                mu, eta, kappa, dim, ncq, nmu, neta, nkappa, ijklnum, potnum);        
    }    
    else if (bondtype==1) {
        gpuQuadrupletb(u, xij, xik, xil, qi, qj, qk, ql, ti, tj, tk, tl, ai, aj, ak, al, 
                mu, eta, kappa, dim, ncq, nmu, neta, nkappa, ijklnum, potnum);
    }            
}
template void gpuQuadruplet(double *, double *, double *, double *, double *, double *, double *, 
        double *, int *, int *, int *, int *, int *, int *, int *, int *, double *, double *, 
        int*, int, int, int, int, int, int, int, int);
template void gpuQuadruplet(float *, float *, float *, float *, float *, float *, float *, 
        float *, int *, int *, int *, int *, int *, int *, int *, int *, float *, float *, 
        int *, int, int, int, int, int, int, int, int);

template <typename T> void gpuComputeQuadrupletEnergyForce(T *__restrict__ u,
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
    gpuArraySetValue(u_xij, (T) 0.0, dim*ijklnum);
    gpuArraySetValue(u_xik, (T) 0.0, dim*ijklnum);
    gpuArraySetValue(u_xil, (T) 0.0, dim*ijklnum);
    
#ifdef HAVE_ENZYME        
    gpuArraySetValue(d_u, (T) 1.0, ijklnum);
    if (bondtype==0) {
        gpuQuadrupletaGradient(u, d_u, u_xij, u_xik, u_xil, xij, xik, xil, qi, qj, qk, ql, ti, tj, tk, tl, ai, aj, ak, al,
                mu, eta, kappa, dim, ncq, nmu, neta, nkappa, ijklnum, potnum);        
    }    
    else if (bondtype==1) {
        gpuQuadrupletbGradient(u, d_u, u_xij, u_xik, u_xil, xij, xik, xil, qi, qj, qk, ql, ti, tj, tk, tl, ai, aj, ak, al, 
                mu, eta, kappa, dim, ncq, nmu, neta, nkappa, ijklnum, potnum);
    } 
#else
    if (bondtype==0) {
        gpuQuadrupletaGradient(u, u_xij, u_xik, u_xil, xij, xik, xil, qi, qj, qk, ql, ti, tj, tk, tl, ai, aj, ak, al,
                mu, eta, kappa, dim, ncq, nmu, neta, nkappa, ijklnum, potnum);        
    }    
    else if (bondtype==1) {
        gpuQuadrupletbGradient(u, u_xij, u_xik, u_xil, xij, xik, xil, qi, qj, qk, ql, ti, tj, tk, tl, ai, aj, ak, al, 
                mu, eta, kappa, dim, ncq, nmu, neta, nkappa, ijklnum, potnum);
    }     
#endif    
    //T *d_u = new T[ijklnum];
//     gpuArraySetValue(d_u, (T) 1.0, ijklnum);
//     gpuArraySetValue(u_xij, (T) 0.0, dim*ijklnum);
//     gpuArraySetValue(u_xik, (T) 0.0, dim*ijklnum);
//     gpuArraySetValue(u_xil, (T) 0.0, dim*ijklnum);
//     
// #ifdef HAVE_ENZYME    
//     __enzyme_autodiff((void*)gpuQuadruplet<T>, 
//                         enzyme_dup, u, d_u,
//                         enzyme_dup, xij, u_xij,
//                         enzyme_dup, xik, u_xik,
//                         enzyme_dup, xil, u_xil,
//                         enzyme_const, qi,
//                         enzyme_const, qj,
//                         enzyme_const, qk,
//                         enzyme_const, ql,            
//                         enzyme_const, ti,
//                         enzyme_const, tj,
//                         enzyme_const, tk,      
//                         enzyme_const, tl,                  
//                         enzyme_const, ai,
//                         enzyme_const, aj,
//                         enzyme_const, ak,        
//                         enzyme_const, al,                    
//                         // TODO can make this non const to get du/dmu
//                         enzyme_const, mu,
//                         enzyme_const, eta,
//                         enzyme_const, kappa,
//                         dim, ncq, nmu, neta, nkappa, ijklnum, potnum, bondtype);  
// #endif    
}
template void gpuComputeQuadrupletEnergyForce(double *, double *, double *, double *, double *, double *, 
        double *, double *, double *, double *, double *, double *, int *, int *, int *, int *, int *, 
        int *, int *, int *, double *, double *, int*, int, int, int, int, int, int, int, int);
template void gpuComputeQuadrupletEnergyForce(float *, float *, float *, float *, float *, float *, 
        float *, float *, float *, float *, float *, float *, int *, int *, int *, int *, int *, 
        int *, int *, int *, float *, float *, int *, int, int, int, int, int, int, int, int);

template <typename T> void gpuPairDensityGradient(T *__restrict__ u,
                                    T *__restrict__ d_u,                     
                                    T *__restrict__ u_rho,   
                                    T *__restrict__ rho,   
                                    T *__restrict__ mu, 
                                    T *__restrict__ eta, 
                                    int *__restrict__ kappa, 
                                    int nrho, int nmu, int neta, 
                                    int nkappa, int inum, int potnum)
{
    gpuArraySetValue(u_rho, (T) 0.0, inum);    
    
#ifdef HAVE_ENZYME            
    gpuArraySetValue(d_u, (T) 1.0, inum);
    gpuPaircDensityGradient(u, du, u_rho, rho, mu, eta, kappa, nrho, nmu, neta, nkappa, inum, potnum);
#else
    gpuPaircDensityGradient(u, u_rho, rho, mu, eta, kappa, nrho, nmu, neta, nkappa, inum, potnum);
#endif    
}
template void gpuPairDensityGradient(double *, double *, double *, double *, double *, double *,
        int*, int, int, int, int, int, int);
template void gpuPairDensityGradient(float *, float *, float *, float *, float *, float *,
        int *, int, int, int, int, int, int);

template <typename T> void gpuTripletDensityGradient(T *__restrict__ u,
                                    T *__restrict__ d_u,                     
                                    T *__restrict__ u_rho,   
                                    T *__restrict__ rho,   
                                    T *__restrict__ mu, 
                                    T *__restrict__ eta, 
                                    int *__restrict__ kappa, 
                                    int nrho, int nmu, int neta, 
                                    int nkappa, int inum, int potnum)
{
    gpuArraySetValue(u_rho, (T) 0.0, inum);    
    
#ifdef HAVE_ENZYME            
    gpuArraySetValue(d_u, (T) 1.0, inum);
    gpuTripletcDensityGradient(u, du, u_rho, rho, mu, eta, kappa, nrho, nmu, neta, nkappa, inum, potnum);
#else
    gpuTripletcDensityGradient(u, u_rho, rho, mu, eta, kappa, nrho, nmu, neta, nkappa, inum, potnum);
#endif    
}
template void gpuTripletDensityGradient(double *, double *, double *, double *, double *, double *,
        int*, int, int, int, int, int, int);
template void gpuTripletDensityGradient(float *, float *, float *, float *, float *, float *,
        int *, int, int, int, int, int, int);

// void gpuElectronDensity(dstype *rhoi, dstype *rhoij, int *pairnum, int *pairnumsum, int inum) 
// {
//     for (int i=0; i<inum; i++) {
//         int jnum = pairnum[i];
//         int start = pairnumsum[i];
//         rhoi[i] = 0.0;
//         for (int j=0; j<jnum; j++) 
//             rhoi[i] += rhoij[start+j];        
//     }
// }
// 
// 
// void gpuEmbedingForce(dstype *fij, dstype *d_rhoi, int *pairnum, int *pairnumsum, int inum)
// {    
//     for (int i=0; i<inum; i++) {
//         int jnum = pairnum[i];
//         int start = pairnumsum[i];
//         for (int j=0; j<jnum; j++)             
//             fij[start+j] = d_rhoi[i]*fij[start+j];                
//     }        
// }
//  
// template <typename T> void gpuPairDensity(T *__restrict__ u,
//                                     T *__restrict__ d_u,                     
//                                     T *__restrict__ u_rho,   
//                                     T *__restrict__ rho,   
//                                     T *__restrict__ mu, 
//                                     T *__restrict__ eta, 
//                                     int *__restrict__ kappa, 
//                                     int nrho, int nmu, int neta, 
//                                     int nkappa, int inum, int potnum)
// {
//     //T *d_u = new T[inum];
//     gpuArraySetValue(d_u, (T) 1.0, inum);
//     gpuArraySetValue(u_rho, (T) 0.0, inum);    
//     
// #ifdef HAVE_ENZYME    
//     __enzyme_autodiff((void*)gpuPaircDensity<T>, 
//                         enzyme_dup, u, d_u,
//                         enzyme_dup, rho, u_rho,                        
//                         // TODO can make this non const to get du/dmu
//                         enzyme_const, mu,
//                         enzyme_const, eta,
//                         enzyme_const, kappa,
//                         1, nmu, neta, nkappa, inum, potnum);  
// #endif
// }
// template void gpuPairDensity(double *, double *, double *, double *, double *, double *,
//         int*, int, int, int, int, int, int);
// template void gpuPairDensity(float *, float *, float *, float *, float *, float *,
//         int *, int, int, int, int, int, int);
       
// template <typename T> void gpuTripletDensity(T *__restrict__ u,
//                                     T *__restrict__ d_u,                             
//                                     T *__restrict__ u_rho,   
//                                     T *__restrict__ rho,   
//                                     T *__restrict__ mu, 
//                                     T *__restrict__ eta, 
//                                     int *__restrict__ kappa, 
//                                     int nrho, int nmu, int neta, 
//                                     int nkappa, int inum, int potnum)
// {
//     //T *d_u = new T[inum];
//     gpuArraySetValue(d_u, (T) 1.0, inum);
//     gpuArraySetValue(u_rho, (T) 0.0, inum);    
//     
// #ifdef HAVE_ENZYME    
//     __enzyme_autodiff((void*)gpuPaircDensity<T>, 
//                         enzyme_dup, u, d_u,
//                         enzyme_dup, rho, u_rho,                        
//                         // TODO can make this non const to get du/dmu
//                         enzyme_const, mu,
//                         enzyme_const, eta,
//                         enzyme_const, kappa,
//                         1, nmu, neta, nkappa, inum, potnum);  
// #endif    
// }
// template void gpuTripletDensity(double *, double *, double *, double *, double *, double *,
//         int*, int, int, int, int, int, int);
// template void gpuTripletDensity(float *, float *, float *, float *, float *, float *,
//         int *, int, int, int, int, int, int);


// template <typename T> void primal_gpuLJ1pot(T *__restrict__ u,
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
// template <typename T> void gpuLJ1pot(T *__restrict__ u,
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
//     __enzyme_autodiff((void*)primal_gpuLJ1pot<T>, 
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
// template <typename T> void gpuLJ1pot(T *__restrict__ u, T *__restrict__ f, T *__restrict__ xij, T *__restrict__ qi, T *__restrict__ qj, int *__restrict__ ti, int *__restrict__ tj, int *__restrict__ ai, int *__restrict__ aj, T *__restrict__ mu, T *__restrict__ eta, int *__restrict__ kappa, int dim, int ncq, int nmu, int neta, int nkappa, int ng, int potnum)
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
// template void gpuLJ1pot(double *, double *, double *, double *, double *, int *, int *, int *, int *, double *, double *, int*, int, int, int, int, int, int, int);
// template void gpuLJ1pot(float *, float *, float *, float *, float *, int *, int *, int *, int *, float *, float *, int *, int, int, int, int, int, int, int);
// 
// template <typename T> void gpuLJ2pot(T *u, T *f, T *xij, T *qi, T *qj, int *ti, int *tj, int *ai, int *aj, T *mu, T *eta, int *kappa, int dim, int ncq, int nmu, int neta, int nkappa, int ng, int potnum)
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
// template void gpuLJ2pot(double *, double *, double *, double *, double *, int *, int *, int *, int *, double *, double *, int*, int, int, int, int, int, int, int);
// template void gpuLJ2pot(float *, float *, float *, float *, float *, int *, int *, int *, int *, float *, float *, int *, int, int, int, int, int, int, int);
// 
// void gpuComputePair1EnergyForce(dstype *eij, dstype *fij, dstype *xij, dstype *qi, dstype *qj, 
//         int *ti, int *tj, int *ai, int *aj, dstype *mu, dstype *eta, int *kappa, int dim, 
//         int ncq, int nmu, int neta, int nkappa, int ijnum, int potnum, int bondtype)
// {
//     if (bondtype==0) {
//         //gpuPaira(eij, xij, qi, qj, ti, tj, ai, aj, mu, eta, kappa, dim, ncq, nmu, neta, nkappa, ijnum, potnum);    
//         gpuLJ1pot(eij, fij, xij, qi, qj, ti, tj, ai, aj, mu, eta, kappa, dim, ncq, nmu, neta, nkappa, ijnum, potnum);    
//     }
//     
//     if (bondtype==1) {
//         gpuPairb(eij, xij, qi, qj, ti, tj, ai, aj, mu, eta, kappa, dim, ncq, nmu, neta, nkappa, ijnum, potnum);    
//     }
//     
//     if (bondtype==2) {
//         gpuPairc(eij, xij, qi, qj, ti, tj, ai, aj, mu, eta, kappa, dim, ncq, nmu, neta, nkappa, ijnum, potnum);    
//     }
//     
//     if (bondtype==3) {
//         gpuTripletcPair(eij, xij, qi, qj, ti, tj, ai, aj, mu, eta, kappa, dim, ncq, nmu, neta, nkappa, ijnum, potnum);    
//     }
// }

// void gpuElectronDensity(dstype *rhoi, dstype *rhoij, int *pairnum, int *pairnumsum, int inum) 
// {
//     for (int i=0; i<inum; i++) {
//         int jnum = pairnum[i];
//         int start = pairnumsum[i];
//         rhoi[i] = 0.0;
//         for (int j=0; j<jnum; j++) 
//             rhoi[i] += rhoij[start+j];        
//     }
// }
// void gpuPaircDensity(dstype *u, dstype *urho, dstype *rho, dstype *mu, dstype *eta, int *kappa, 
//         int nrho, int nmu, int neta, int nkappa, int ng, int potnum)
// {
//     gpuPaircDensity(u, rho, mu, eta, kappa, 1, nmu, neta, nkappa, ng, potnum);
// }
// void gpuEmbedingForce(dstype *fij, dstype *drhoi, int *pairnum, int *pairnumsum, int inum)
// {    
//     for (int i=0; i<inum; i++) {
//         int jnum = pairnum[i];
//         int start = pairnumsum[i];
//         for (int j=0; j<jnum; j++)             
//             fij[start+j] = drhoi[i]*fij[start+j];                
//     }        
// }
        
// void gpuComputeTripletEnergyForce(dstype *eijk, dstype *fij, dstype *fik, dstype *xij, dstype *xik, dstype *qi, dstype *qj, dstype *qk, 
//      int *ti, int *tj, int *tk, int *ai, int *aj, int *ak, dstype *mu, dstype *eta, int *kappa, int dim, int ncq, 
//      int nmu, int neta, int nkappa, int ijknum, int potnum, int bondtype)        
// {
//     if (bondtype==0) {
//         gpuTripleta(eijk, xij, xik, qi, qj, qk, ti, tj, tk, ai, aj, ak, mu, eta, kappa, dim, ncq, nmu, neta, nkappa, ijknum, potnum);        
//     }
//     
//     if (bondtype==1) {
//         gpuTripletb(eijk, xij, xik, qi, qj, qk, ti, tj, tk, ai, aj, ak, mu, eta, kappa, dim, ncq, nmu, neta, nkappa, ijknum, potnum);
//     }
//     
//     if (bondtype==2) {
//         gpuTripletc(eijk, xij, xik, qi, qj, qk, ti, tj, tk, ai, aj, ak, mu, eta, kappa, dim, ncq, nmu, neta, nkappa, ijknum, potnum);
//     }
// }
// 
// // void gpuTripletcDensity(dstype *u, dstype *urho, dstype *rho, dstype *mu, dstype *eta, int *kappa, 
// //         int nrho, int nmu, int neta, int nkappa, int ng, int potnum)
// // {
// //     gpuTripletcDensity(u, rho, mu, eta, kappa, 1, nmu, neta, nkappa, ng, potnum);
// // }
// 
// void gpuComputeQuadrupletEnergyForce(dstype *eijkl, dstype *fij, dstype *fik, dstype *fil, dstype *xij, dstype *xik, dstype *xil, 
//       dstype *qi, dstype *qj, dstype *qk, dstype *ql, int *ti, int *tj, int *tk, int *tl, int *ai, int *aj, int *ak, int *al,
//       dstype *mu, dstype *eta, int *kappa, int dim, int ncq, int nmu, int neta, int nkappa, int ijklnum, int potnum, int bondtype)        
// {
//     if (bondtype==0) {
//         gpuQuadrupleta(eijkl, xij, xik, xil, qi, qj, qk, ql, ti, tj, tk, tl, ai, aj, ak, al,
//                 mu, eta, kappa, dim, ncq, nmu, neta, nkappa, ijklnum, potnum);        
//     }
//     
//     if (bondtype==1) {
//         gpuQuadrupletb(eijkl, xij, xik, xil, qi, qj, qk, ql, ti, tj, tk, tl, ai, aj, ak, al, 
//                 mu, eta, kappa, dim, ncq, nmu, neta, nkappa, ijklnum, potnum);
//     }        
// }
// 

#endif


