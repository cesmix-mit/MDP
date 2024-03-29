template <typename T> void cpuPairb1(T *u, T *xij, T *qi, T *qj, int *ti, int *tj, int *ai, int *aj, T *mu, T *eta, int *kappa, int dim, int ncq, int nmu, int neta, int nkappa, int ng)
{
	#pragma omp parallel for
	for (int i = 0; i <ng; i++) {
		T mu1 = mu[0];
		T mu2 = mu[1];
		T mu3 = mu[2];
		T xij1 = xij[0 + i*3];
		T xij2 = xij[1 + i*3];
		T xij3 = xij[2 + i*3];
		u[i] = -mu1*exp(2*mu2*(mu3 - sqrt(pow(xij1, 2) + pow(xij2, 2) + pow(xij3, 2))));
	}
}

template <typename T> void cpuPairb2(T *u, T *xij, T *qi, T *qj, int *ti, int *tj, int *ai, int *aj, T *mu, T *eta, int *kappa, int dim, int ncq, int nmu, int neta, int nkappa, int ng)
{
	#pragma omp parallel for
	for (int i = 0; i <ng; i++) {
		T mu4 = mu[3];
		T mu5 = mu[4];
		T mu6 = mu[5];
		T xij1 = xij[0 + i*3];
		T xij2 = xij[1 + i*3];
		T xij3 = xij[2 + i*3];
		T x0 = -mu5;
		u[i] = mu4*(-pow(mu6 + x0, 2) + pow(x0 + sqrt(pow(xij1, 2) + pow(xij2, 2) + pow(xij3, 2)), 2));
	}
}

template <typename T> void cpuPairb(T *u, T *xij, T *qi, T *qj, int *ti, int *tj, int *ai, int *aj, T *mu, T *eta, int *kappa, int dim, int ncq, int nmu, int neta, int nkappa, int ng, int potnum)
{
	if (potnum == 1)
		cpuPairb1(u, xij, qi, qj, ti, tj, ai, aj, mu, eta, kappa, dim, ncq, nmu, neta, nkappa, ng);
	else if (potnum == 2)
		cpuPairb2(u, xij, qi, qj, ti, tj, ai, aj, mu, eta, kappa, dim, ncq, nmu, neta, nkappa, ng);
}
template void cpuPairb(double *, double *, double *, double *, int *, int *, int *, int *, double *, double *, int*, int, int, int, int, int, int, int);
template void cpuPairb(float *, float *, float *, float *, int *, int *, int *, int *, float *, float *, int *, int, int, int, int, int, int, int);
