template <typename T> void gpuTripleta(T *u, T *xij, T *xik, T *qi, T *qj, T *qk, int *ti, int *tj, int *tk, int *ai, int *aj, int *ak, T *mu, T *eta, int *kappa, int dim, int ncq, int nmu, int neta, int nkappa, int ng, int potnum)
{
}
template void gpuTripleta(double *, double *, double *, double *, double *, double *, int *, int *, int *, int *, int *, int *, double *, double *, int*, int, int, int, int, int, int, int);
template void gpuTripleta(float *, float *, float *, float *, float *, float *, int *, int *, int *, int *, int *, int *, float *, float *, int *, int, int, int, int, int, int, int);
