template <typename T> void cpuTripletbGradient1(T *u, T *u_xij, T *u_xik, T *xij, T *xik, T *qi, T *qj, T *qk, int *ti, int *tj, int *tk, int *ai, int *aj, int *ak, T *mu, T *eta, int *kappa, int dim, int ncq, int nmu, int neta, int nkappa, int ng)
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
		T x0 = xij1*xik1 + xij2*xik2 + xij3*xik3;
		T x1 = pow(xij1, 2) + pow(xij2, 2) + pow(xij3, 2);
		T x2 = pow(x1, -1.0/2.0);
		T x3 = pow(xik1, 2) + pow(xik2, 2) + pow(xik3, 2);
		T x4 = pow(x3, -1.0/2.0);
		T x5 = x2*x4;
		T x6 = -mu2 + x0*x5;
		T x7 = 2*x5;
		T x8 = 2*x0;
		T x9 = x4*x8/pow(x1, 3.0/2.0);
		T x10 = mu1*x6;
		T x11 = x2*x8/pow(x3, 3.0/2.0);
		u[i] = mu1*pow(x6, 2);
		u_xij[0 + i*3] = x10*(x7*xik1 - x9*xij1);
		u_xij[1 + i*3] = x10*(x7*xik2 - x9*xij2);
		u_xij[2 + i*3] = x10*(x7*xik3 - x9*xij3);
		u_xik[0 + i*3] = x10*(-x11*xik1 + x7*xij1);
		u_xik[1 + i*3] = x10*(-x11*xik2 + x7*xij2);
		u_xik[2 + i*3] = x10*(-x11*xik3 + x7*xij3);
	}
}

template <typename T> void cpuTripletbGradient2(T *u, T *u_xij, T *u_xik, T *xij, T *xik, T *qi, T *qj, T *qk, int *ti, int *tj, int *tk, int *ai, int *aj, int *ak, T *mu, T *eta, int *kappa, int dim, int ncq, int nmu, int neta, int nkappa, int ng)
{
	#pragma omp parallel for
	for (int i = 0; i <ng; i++) {
		T mu1 = mu[0];
		T mu2 = mu[1];
		T mu3 = mu[2];
		T mu4 = mu[3];
		T xij1 = xij[0 + i*3];
		T xij2 = xij[1 + i*3];
		T xij3 = xij[2 + i*3];
		T xik1 = xik[0 + i*3];
		T xik2 = xik[1 + i*3];
		T xik3 = xik[2 + i*3];
		T x0 = xij1*xik1;
		T x1 = xij2*xik2;
		T x2 = xij3*xik3;
		T x3 = pow(xij1, 2) + pow(xij2, 2) + pow(xij3, 2);
		T x4 = pow(xik1, 2) + pow(xik2, 2) + pow(xik3, 2);
		T x5 = sqrt(-2*x0 - 2*x1 - 2*x2 + x3 + x4);
		T x6 = mu4 - x5;
		T x7 = x0 + x1 + x2;
		T x8 = pow(x3, -1.0/2.0);
		T x9 = pow(x4, -1.0/2.0);
		T x10 = x8*x9;
		T x11 = -mu2 + x10*x7;
		T x12 = 2*mu3*x6/x5;
		T x13 = 2*x10;
		T x14 = 2*x7;
		T x15 = x14*x9/pow(x3, 3.0/2.0);
		T x16 = mu1*x11;
		T x17 = x14*x8/pow(x4, 3.0/2.0);
		u[i] = mu1*pow(x11, 2) + mu3*pow(x6, 2);
		u_xij[0 + i*3] = -x12*(xij1 - xik1) + x16*(x13*xik1 - x15*xij1);
		u_xij[1 + i*3] = -x12*(xij2 - xik2) + x16*(x13*xik2 - x15*xij2);
		u_xij[2 + i*3] = -x12*(xij3 - xik3) + x16*(x13*xik3 - x15*xij3);
		u_xik[0 + i*3] = -x12*(-xij1 + xik1) + x16*(x13*xij1 - x17*xik1);
		u_xik[1 + i*3] = -x12*(-xij2 + xik2) + x16*(x13*xij2 - x17*xik2);
		u_xik[2 + i*3] = -x12*(-xij3 + xik3) + x16*(x13*xij3 - x17*xik3);
	}
}


template <typename T> void cpuTripletbGradient(T *u, T *u_xij, T *u_xik, T *xij, T *xik, T *qi, T *qj, T *qk, int *ti, int *tj, int *tk, int *ai, int *aj, int *ak, T *mu, T *eta, int *kappa, int dim, int ncq, int nmu, int neta, int nkappa, int ng, int potnum)
{
	if (potnum == 1)
		cpuTripletbGradient1(u, u_xij, u_xik, xij, xik, qi, qj, qk, ti, tj, tk, ai, aj, ak, mu, eta, kappa, dim, ncq, nmu, neta, nkappa, ng);
	else if (potnum == 2)
		cpuTripletbGradient2(u, u_xij, u_xik, xij, xik, qi, qj, qk, ti, tj, tk, ai, aj, ak, mu, eta, kappa, dim, ncq, nmu, neta, nkappa, ng);
}
template void cpuTripletbGradient(double *, double *, double *, double *, double *, double *, double *, double *, int *, int *, int *, int *, int *, int *, double *, double *, int*, int, int, int, int, int, int, int);
template void cpuTripletbGradient(float *, float *, float *, float *, float *, float *, float *, float *, int *, int *, int *, int *, int *, int *, float *, float *, int *, int, int, int, int, int, int, int);
