template <typename T> void gpuQuadrupleta(T *u, T *xij, T *xik, T *xil, T *qi, T *qj, T *qk, T *ql, int *ti, int *tj, int *tk, int *tl, int *ai, int *aj, int *ak, int *al, T *mu, T *eta, int *kappa, int dim, int ncq, int nmu, int neta, int nkappa, int ng, int potnum)
{
}
template void gpuQuadrupleta(double *, double *, double *, double *, double *, double *, double *, double *, int *, int *, int *, int *, int *, int *, int *, int *, double *, double *, int*, int, int, int, int, int, int, int);
template void gpuQuadrupleta(float *, float *, float *, float *, float *, float *, float *, float *, int *, int *, int *, int *, int *, int *, int *, int *, float *, float *, int *, int, int, int, int, int, int, int);

template <typename T> void gpuQuadrupletaGradient(T *u, T *du, T *u_xij, T *u_xik, T *u_xil, T *xij, T *xik, T *xil, T *qi, T *qj, T *qk, T *ql, int *ti, int *tj, int *tk, int *tl, int *ai, int *aj, int *ak, int *al, T *mu, T *eta, int *kappa, int dim, int ncq, int nmu, int neta, int nkappa, int ng, int potnum)
{
}
template void gpuQuadrupletaGradient(double *, double *, double *, double *, double *, double *, double *, double *, double *, double *, double *, double *, int *, int *, int *, int *, int *, int *, int *, int *, double *, double *, int*, int, int, int, int, int, int, int);
template void gpuQuadrupletaGradient(float *, float *, float *, float *, float *, float *, float *, float *, float *, float *, float *, float *, int *, int *, int *, int *, int *, int *, int *, int *, float *, float *, int *, int, int, int, int, int, int, int);
