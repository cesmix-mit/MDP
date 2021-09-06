template <typename T> void gpuSinglea(T *u, T *xi, T *qi, int *ti, int *ai, T *mu, T *eta, int *kappa, int dim, int ncq, int nmu, int neta, int nkappa, int ng, int potnum)
{
}
template void gpuSinglea(double *, double *, double *, int *, int *, double *, double *, int*, int, int, int, int, int, int, int);
template void gpuSinglea(float *, float *, float *, int *, int *, float *, float *, int *, int, int, int, int, int, int, int);

template <typename T> void gpuSingleaGradient(T *u, T *du, T *u_xi, T *xi, T *qi, int *ti, int *ai, T *mu, T *eta, int *kappa, int dim, int ncq, int nmu, int neta, int nkappa, int ng, int potnum)
{
}
template void gpuSingleaGradient(double *, double *, double *, double *, double *, int *, int *, double *, double *, int*, int, int, int, int, int, int, int);
template void gpuSingleaGradient(float *, float *, float *, float *, float *, int *, int *, float *, float *, int *, int, int, int, int, int, int, int);
