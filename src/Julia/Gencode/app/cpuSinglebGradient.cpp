template <typename T> void cpuSinglebGradient(T *u, T *u_xi, T *xi, T *qi, int *ti, int *ai, T *mu, T *eta, int *kappa, int dim, int ncq, int nmu, int neta, int nkappa, int ng, int potnum)
{
}
template void cpuSinglebGradient(double *, double *, double *, double *, int *, int *, double *, double *, int*, int, int, int, int, int, int, int);
template void cpuSinglebGradient(float *, float *, float *, float *, int *, int *, float *, float *, int *, int, int, int, int, int, int, int);
