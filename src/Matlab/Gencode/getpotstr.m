%***************************************************************************                               
%                     Molecular Dynamics Potentials (MDP)
%                            CESMIX-MIT Project  
%  
% Contributing authors: Ngoc-Cuong Nguyen (cuongng@mit.edu, exapde@gmail.com)
%***************************************************************************

function [sp0, sp1, sp2, sp3] = getpotstr(pot)

if pot==0
    sp0 = "(T *u, T *rho, T *mu, T *eta, int *kappa, int nrho, int nmu, int neta, int nkappa, int ng)\n";
    sp1 = "(double *, double *, double *, double *, int*, int, int, int, int, int);\n";
    sp2 = "(float *, float *, float *, float *, int *, int, int, int, int, int);\n";
    sp3 = "(u, rho, mu, eta, kappa, nrho, nmu, neta, nkappa, ng);\n";
elseif pot==1
    sp0 = "(T *u, T *xi, T *qi, int *ti, int *ai, T *mu, T *eta, int *kappa, int dim, int ncq, int nmu, int neta, int nkappa, int ng)\n";
    sp1 = "(double *, double *, double *, int *, int *, double *, double *, int*, int, int, int, int, int, int);\n";
    sp2 = "(float *, float *, float *, int *, int *, float *, float *, int *, int, int, int, int, int, int);\n";
    sp3 = "(u, xi, qi, ti, ai, mu, eta, kappa, dim, ncq, nmu, neta, nkappa, ng);\n";
elseif pot==2
    sp0 = "(T *u, T *xij, T *qi, T *qj, int *ti, int *tj, int *ai, int *aj, T *mu, T *eta, int *kappa, int dim, int ncq, int nmu, int neta, int nkappa, int ng)\n";
    sp1 = "(double *, double *, double *, double *, int *, int *, int *, int *, double *, double *, int*, int, int, int, int, int, int);\n";
    sp2 = "(float *, float *, float *, float *, int *, int *, int *, int *, float *, float *, int *, int, int, int, int, int, int);\n";
    sp3 = "(u, xij, qi, qj, ti, tj, ai, aj, mu, eta, kappa, dim, ncq, nmu, neta, nkappa, ng);\n";
elseif pot==3
    sp0 = "(T *u, T *xij, T *xik, T *qi, T *qj, T *qk, int *ti, int *tj, int *tk, int *ai, int *aj, int *ak, T *mu, T *eta, int *kappa, int dim, int ncq, int nmu, int neta, int nkappa, int ng)\n";
    sp1 = "(double *, double *, double *, double *, double *, double *, int *, int *, int *, int *, int *, int *, double *, double *, int*, int, int, int, int, int, int);\n";
    sp2 = "(float *, float *, float *, float *, float *, float *, int *, int *, int *, int *, int *, int *, float *, float *, int *, int, int, int, int, int, int);\n";
    sp3 = "(u, xij, xik, qi, qj, qk, ti, tj, tk, ai, aj, ak, mu, eta, kappa, dim, ncq, nmu, neta, nkappa, ng);\n";
elseif pot==4
    sp0 = "(T *u, T *xij, T *xik, T *xil, T *qi, T *qj, T *qk, T *ql, int *ti, int *tj, int *tk, int *tl, int *ai, int *aj, int *ak, int *al, T *mu, T *eta, int *kappa, int dim, int ncq, int nmu, int neta, int nkappa, int ng)\n";
    sp1 = "(double *, double *, double *, double *, double *, double *, double *, double *, int *, int *, int *, int *, int *, int *, int *, int *, double *, double *, int*, int, int, int, int, int, int);\n";
    sp2 = "(float *, float *, float *, float *, float *, float *, float *, float *, int *, int *, int *, int *, int *, int *, int *, int *, float *, float *, int *, int, int, int, int, int, int);\n";
    sp3 = "(u, xij, xik, xil, qi, qj, qk, ql, ti, tj, tk, tl, ai, aj, ak, al, mu, eta, kappa, dim, ncq, nmu, neta, nkappa, ng);\n";    
end

