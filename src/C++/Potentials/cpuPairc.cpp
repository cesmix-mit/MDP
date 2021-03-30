template <typename T> void cpuPairc(T *u, T *xij, T *qi, T *qj, int *ti, int *tj, int *ai, int *aj, T *mu, T *eta, int *kappa, int dim, int ncq, int nmu, int neta, int nkappa, int ng, int potnum)
{
}
template void cpuPairc(double *, double *, double *, double *, int *, int *, int *, int *, double *, double *, int*, int, int, int, int, int, int, int);
template void cpuPairc(float *, float *, float *, float *, int *, int *, int *, int *, float *, float *, int *, int, int, int, int, int, int, int);

template <typename T> void cpuPaircDensity(T *u, T *rho, T *mu, T *eta, int *kappa, int nrho, int nmu, int neta, int nkappa, int ng, int potnum)
{
}
template void cpuPaircDensity(double *, double *, double *, double *, int*, int, int, int, int, int, int);
template void cpuPaircDensity(float *, float *, float *, float *, int *, int, int, int, int, int, int);

