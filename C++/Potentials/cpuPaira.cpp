template <typename T> void cpuPaira1(T *u, T *xij, T *qi, T *qj, int *ti, int *tj, int *ai, int *aj, T *mu, T *eta, int *kappa, int dim, int ncq, int nmu, int neta, int nkappa, int ng)
{
	#pragma omp parallel for
	for (int i = 0; i <ng; i++) {
		T xij1 = xij[0 + i*3];
		T xij2 = xij[1 + i*3];
		T xij3 = xij[2 + i*3];
		T x0 = sqrt(pow(xij1, 2) + pow(xij2, 2) + pow(xij3, 2));
		T x1 = tanh(100.0*x0 - 400.0);
		u[i] = 0.0081259386830030578*x1 + (0.5*x1 + 0.5)*((0.18277826364023042 - 0.034282917879905832*x0)*pow(x0 - 4, 3) - 0.016251877366006116) - 0.0081259386830030578 + (13946.714966258751*exp(-36.644389638721982*x0) + 39124.468185401303*exp(-10.791187546931475*x0) + 21502.880153205104*exp(-4.6140460608291409*x0) + 2161.6449001348501*exp(-2.3089698852925573*x0))/x0;
	}
}

template <typename T> void cpuPaira(T *u, T *xij, T *qi, T *qj, int *ti, int *tj, int *ai, int *aj, T *mu, T *eta, int *kappa, int dim, int ncq, int nmu, int neta, int nkappa, int ng, int potnum)
{
	if (potnum == 1)
		cpuPaira1(u, xij, qi, qj, ti, tj, ai, aj, mu, eta, kappa, dim, ncq, nmu, neta, nkappa, ng);
}
template void cpuPaira(double *, double *, double *, double *, int *, int *, int *, int *, double *, double *, int*, int, int, int, int, int, int, int);
template void cpuPaira(float *, float *, float *, float *, int *, int *, int *, int *, float *, float *, int *, int, int, int, int, int, int, int);
