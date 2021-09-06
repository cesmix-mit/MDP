template <typename T> void cpuTripletcGradient1(T *u, T *u_xij, T *u_xik, T *xij, T *xik, T *qi, T *qj, T *qk, int *ti, int *tj, int *tk, int *ai, int *aj, int *ak, T *mu, T *eta, int *kappa, int dim, int ncq, int nmu, int neta, int nkappa, int ng)
{
	#pragma omp parallel for
	for (int i = 0; i <ng; i++) {
		T mu3 = mu[2];
		T mu4 = mu[3];
		T mu5 = mu[4];
		T mu6 = mu[5];
		T mu7 = mu[6];
		T mu8 = mu[7];
		T xij1 = xij[0 + i*3];
		T xij2 = xij[1 + i*3];
		T xij3 = xij[2 + i*3];
		T xik1 = xik[0 + i*3];
		T xik2 = xik[1 + i*3];
		T xik3 = xik[2 + i*3];
		T x0 = pow(mu5, 2);
		T x1 = pow(mu6, 2);
		T x2 = xij1*xik1 + xij2*xik2 + xij3*xik3;
		T x3 = pow(xij1, 2) + pow(xij2, 2) + pow(xij3, 2);
		T x4 = sqrt(x3);
		T x5 = 1.0/x4;
		T x6 = pow(xik1, 2) + pow(xik2, 2) + pow(xik3, 2);
		T x7 = sqrt(x6);
		T x8 = 1.0/x7;
		T x9 = x5*x8;
		T x10 = -mu7 + x2*x9;
		T x11 = x1 + pow(x10, 2);
		T x12 = x4 - x7;
		T x13 = pow(mu3, mu8)*pow(x12, mu8);
		T x14 = mu4*exp(x13);
		T x15 = x14*(-x0/x11 + x0/x1 + 1);
		T x16 = mu8*x13*x15/x12;
		T x17 = x16*x5;
		T x18 = 2*x9;
		T x19 = 2*x2;
		T x20 = x19*x8/pow(x3, 3.0/2.0);
		T x21 = x0*x10*x14/pow(x11, 2);
		T x22 = x16*x8;
		T x23 = x19*x5/pow(x6, 3.0/2.0);
		u[i] = x15;
		u_xij[0 + i*3] = x17*xij1 + x21*(x18*xik1 - x20*xij1);
		u_xij[1 + i*3] = x17*xij2 + x21*(x18*xik2 - x20*xij2);
		u_xij[2 + i*3] = x17*xij3 + x21*(x18*xik3 - x20*xij3);
		u_xik[0 + i*3] = x21*(x18*xij1 - x23*xik1) - x22*xik1;
		u_xik[1 + i*3] = x21*(x18*xij2 - x23*xik2) - x22*xik2;
		u_xik[2 + i*3] = x21*(x18*xij3 - x23*xik3) - x22*xik3;
	}
}

template <typename T> void cpuTripletcPairGradient1(T *u, T *u_xij, T *xij, T *qi, T *qj, int *ti, int *tj, int *ai, int *aj, T *mu, T *eta, int *kappa, int dim, int ncq, int nmu, int neta, int nkappa, int ng)
{
	#pragma omp parallel for
	for (int i = 0; i <ng; i++) {
		T mu1 = mu[0];
		T mu2 = mu[1];
		T xij1 = xij[0 + i*3];
		T xij2 = xij[1 + i*3];
		T xij3 = xij[2 + i*3];
		T x0 = sqrt(pow(xij1, 2) + pow(xij2, 2) + pow(xij3, 2));
		T x1 = mu1*exp(-mu2*x0);
		T x2 = mu2*x1/x0;
		u[0 + i*1] = -x1;
		u_xij[0 + i*3] = x2*xij1;
		u_xij[1 + i*3] = x2*xij2;
		u_xij[2 + i*3] = x2*xij3;
	}
}

template <typename T> void cpuTripletcDensityGradient1(T *u, T *u_rho, T *rho, T *mu, T *eta, int *kappa, int nrho, int nmu, int neta, int nkappa, int ng)
{
	#pragma omp parallel for
	for (int i = 0; i <ng; i++) {
		T mu1 = mu[0];
		T mu9 = mu[8];
		T mu10 = mu[9];
		T rho1 = rho[0 + i*1];
		T x0 = pow(mu10, mu9)*pow(rho1, mu9);
		T x1 = x0 + 1;
		T x2 = pow(x1, -0.5/mu9);
		u[i] = x2;
		u_rho[i] = -0.5*x0*x2/(rho1*x1);
	}
}


template <typename T> void cpuTripletcGradient(T *u, T *u_xij, T *u_xik, T *xij, T *xik, T *qi, T *qj, T *qk, int *ti, int *tj, int *tk, int *ai, int *aj, int *ak, T *mu, T *eta, int *kappa, int dim, int ncq, int nmu, int neta, int nkappa, int ng, int potnum)
{
	if (potnum == 1)
		cpuTripletcGradient1(u, u_xij, u_xik, xij, xik, qi, qj, qk, ti, tj, tk, ai, aj, ak, mu, eta, kappa, dim, ncq, nmu, neta, nkappa, ng);
}
template void cpuTripletcGradient(double *, double *, double *, double *, double *, double *, double *, double *, int *, int *, int *, int *, int *, int *, double *, double *, int*, int, int, int, int, int, int, int);
template void cpuTripletcGradient(float *, float *, float *, float *, float *, float *, float *, float *, int *, int *, int *, int *, int *, int *, float *, float *, int *, int, int, int, int, int, int, int);

template <typename T> void cpuTripletcPairGradient(T *u, T *u_xij, T *xij, T *qi, T *qj, int *ti, int *tj, int *ai, int *aj, T *mu, T *eta, int *kappa, int dim, int ncq, int nmu, int neta, int nkappa, int ng, int potnum)
{
	if (potnum == 1)
		cpuTripletcPairGradient1(u, u_xij, xij, qi, qj, ti, tj, ai, aj, mu, eta, kappa, dim, ncq, nmu, neta, nkappa, ng);
}
template void cpuTripletcPairGradient(double *, double *, double *, double *, double *, int *, int *, int *, int *, double *, double *, int*, int, int, int, int, int, int, int);
template void cpuTripletcPairGradient(float *, float *, float *, float *, float *, int *, int *, int *, int *, float *, float *, int *, int, int, int, int, int, int, int);

template <typename T> void cpuTripletcDensityGradient(T *u, T *u_rho, T *rho, T *mu, T *eta, int *kappa, int nrho, int nmu, int neta, int nkappa, int ng, int potnum)
{
	if (potnum == 1)
		cpuTripletcDensityGradient1(u, u_rho, rho, mu, eta, kappa, nrho, nmu, neta, nkappa, ng);
}
template void cpuTripletcDensityGradient(double *, double *, double *, double *, double *, int*, int, int, int, int, int, int);
template void cpuTripletcDensityGradient(float *, float *, float *, float *, float *, int *, int, int, int, int, int, int);
