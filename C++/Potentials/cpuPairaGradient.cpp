template <typename T> void cpuPairaGradient1(T *u, T *u_xij, T *xij, T *qi, T *qj, int *ti, int *tj, int *ai, int *aj, T *mu, T *eta, int *kappa, int dim, int ncq, int nmu, int neta, int nkappa, int ng)
{
	#pragma omp parallel for
	for (int i = 0; i <ng; i++) {
		T xij1 = xij[0 + i*3];
		T xij2 = xij[1 + i*3];
		T xij3 = xij[2 + i*3];
		T x0 = pow(xij1, 2) + pow(xij2, 2) + pow(xij3, 2);
		T x1 = sqrt(x0);
		T x2 = tanh(100.0*x1 - 400.0);
		T x3 = 0.5*x2 + 0.5;
		T x4 = 0.18277826364023042 - 0.034282917879905832*x1;
		T x5 = x1 - 4;
		T x6 = pow(x5, 3);
		T x7 = x4*x6 - 0.016251877366006116;
		T x8 = 1.0/x1;
		T x9 = exp(-36.644389638721982*x1);
		T x10 = exp(-10.791187546931475*x1);
		T x11 = exp(-4.6140460608291409*x1);
		T x12 = exp(-2.3089698852925573*x1);
		T x13 = 39124.468185401303*x10 + 21502.880153205104*x11 + 2161.6449001348501*x12 + 13946.714966258751*x9;
		T x14 = x8*xij1;
		T x15 = 1 - pow(x2, 2);
		T x16 = 0.81259386830030578*x15;
		T x17 = 50.0*x15*x7;
		T x18 = x13/pow(x0, 3.0/2.0);
		T x19 = 0.034282917879905832*x6;
		T x20 = 3*x4*pow(x5, 2);
		T x21 = 511068.85740378097*x9;
		T x22 = 422199.47386261926*x10;
		T x23 = 99215.279467377128*x11;
		T x24 = 4991.1729771076061*x12;
		T x25 = x8*xij2;
		T x26 = x8*xij3;
		u[i] = x13*x8 + 0.0081259386830030578*x2 + x3*x7 - 0.0081259386830030578;
		u_xij[0 + i*3] = x14*x16 + x14*x17 - x18*xij1 + x3*(-x14*x19 + x14*x20) + x8*(-x14*x21 - x14*x22 - x14*x23 - x14*x24);
		u_xij[1 + i*3] = x16*x25 + x17*x25 - x18*xij2 + x3*(-x19*x25 + x20*x25) + x8*(-x21*x25 - x22*x25 - x23*x25 - x24*x25);
		u_xij[2 + i*3] = x16*x26 + x17*x26 - x18*xij3 + x3*(-x19*x26 + x20*x26) + x8*(-x21*x26 - x22*x26 - x23*x26 - x24*x26);
	}
}


template <typename T> void cpuPairaGradient(T *u, T *u_xij, T *xij, T *qi, T *qj, int *ti, int *tj, int *ai, int *aj, T *mu, T *eta, int *kappa, int dim, int ncq, int nmu, int neta, int nkappa, int ng, int potnum)
{
	if (potnum == 1)
		cpuPairaGradient1(u, u_xij, xij, qi, qj, ti, tj, ai, aj, mu, eta, kappa, dim, ncq, nmu, neta, nkappa, ng);
}
template void cpuPairaGradient(double *, double *, double *, double *, double *, int *, int *, int *, int *, double *, double *, int*, int, int, int, int, int, int, int);
template void cpuPairaGradient(float *, float *, float *, float *, float *, int *, int *, int *, int *, float *, float *, int *, int, int, int, int, int, int, int);
