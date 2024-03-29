template <typename T> void cpuQuadrupletb1(T *u, T *xij, T *xik, T *xil, T *qi, T *qj, T *qk, T *ql, int *ti, int *tj, int *tk, int *tl, int *ai, int *aj, int *ak, int *al, T *mu, T *eta, int *kappa, int dim, int ncq, int nmu, int neta, int nkappa, int ng)
{
	#pragma omp parallel for
	for (int i = 0; i <ng; i++) {
		T mu1 = mu[0];
		T mu2 = mu[1];
		T xij1 = xij[0 + i*3];
		T xij2 = xij[1 + i*3];
		T xij3 = xij[2 + i*3];
		T xik1 = xik[0 + i*3];
		T xik2 = xik[1 + i*3];
		T xik3 = xik[2 + i*3];
		T xil1 = xil[0 + i*3];
		T xil2 = xil[1 + i*3];
		T xil3 = xil[2 + i*3];
		T x0 = -xij2;
		T x1 = x0 + xik2;
		T x2 = -xij3;
		T x3 = x2 + xil3;
		T x4 = x0 + xil2;
		T x5 = x2 + xik3;
		T x6 = x1*x3 - x4*x5;
		T x7 = -xij1;
		T x8 = x7 + xil1;
		T x9 = x7 + xik1;
		T x10 = -x3*x9 + x5*x8;
		T x11 = -x1*x8 + x4*x9;
		T x12 = x10*xij2 + x11*xij3 + x6*xij1;
		T x13 = pow(x10, 2) + pow(x11, 2) + pow(x6, 2);
		u[i] = mu1*pow(x12, 2)/x13 + mu2*pow(x12, 4)/pow(x13, 2);
	}
}


template <typename T> void cpuQuadrupletb(T *u, T *xij, T *xik, T *xil, T *qi, T *qj, T *qk, T *ql, int *ti, int *tj, int *tk, int *tl, int *ai, int *aj, int *ak, int *al, T *mu, T *eta, int *kappa, int dim, int ncq, int nmu, int neta, int nkappa, int ng, int potnum)
{
	if (potnum == 1)
		cpuQuadrupletb1(u, xij, xik, xil, qi, qj, qk, ql, ti, tj, tk, tl, ai, aj, ak, al, mu, eta, kappa, dim, ncq, nmu, neta, nkappa, ng);
}
template void cpuQuadrupletb(double *, double *, double *, double *, double *, double *, double *, double *, int *, int *, int *, int *, int *, int *, int *, int *, double *, double *, int*, int, int, int, int, int, int, int);
template void cpuQuadrupletb(float *, float *, float *, float *, float *, float *, float *, float *, int *, int *, int *, int *, int *, int *, int *, int *, float *, float *, int *, int, int, int, int, int, int, int);
