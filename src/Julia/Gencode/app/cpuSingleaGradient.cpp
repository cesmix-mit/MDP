template <typename T> void cpuSingleaGradient(T *u, T *u_xi, T *xi, T *qi, int *ti, int *ai, T *mu, T *eta, int *kappa, int dim, int ncq, int nmu, int neta, int nkappa, int ng, int potnum)
{
}
template void cpuSingleaGradient(double *, double *, double *, double *, int *, int *, double *, double *, int*, int, int, int, int, int, int, int);
template void cpuSingleaGradient(float *, float *, float *, float *, int *, int *, float *, float *, int *, int, int, int, int, int, int, int);
