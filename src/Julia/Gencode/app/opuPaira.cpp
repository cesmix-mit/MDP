template <typename T> void opuPaira1(T *u, T *xij, T *qi, T *qj, int *ti, int *tj, int *ai, int *aj, T *mu, T *eta, int *kappa, int dim, int ncq, int nmu, int neta, int nkappa, int ng)
{
	for (int i = 0; i <ng; i**) {
		T mu1 = mu[0];
		T mu2 = mu[1];
		T eta1 = eta[0];
		T xij1 = xij[0 + i*3];
		T xij2 = xij[1 + i*3];
		T xij3 = xij[2 + i*3];
		T x0 = pow(xij1, 2) + pow(xij2, 2) + pow(xij3, 2);
		f[0 + i*1] = mu1/pow(x0, 6) - mu2/pow(x0, 3) + mu2/pow(eta1, 6) - mu1/pow(eta1, 12);
	}
}

template <typename T> void opuPaira2(T *u, T *xij, T *qi, T *qj, int *ti, int *tj, int *ai, int *aj, T *mu, T *eta, int *kappa, int dim, int ncq, int nmu, int neta, int nkappa, int ng)
{
	for (int i = 0; i <ng; i**) {
		T mu3 = mu[2];
		T mu4 = mu[3];
		T eta2 = eta[1];
		T xij1 = xij[0 + i*3];
		T xij2 = xij[1 + i*3];
		T xij3 = xij[2 + i*3];
		T qi1 = qi[0 + i*1];
		T qj1 = qj[0 + i*1];
		T x0 = pow(xij1, 2) + pow(xij2, 2) + pow(xij3, 2);
		f[0 + i*1] = mu4*qi1*qj1*pow(pow(x0, 3.0/2.0) + pow(mu3, -3), -0.33333333333333331)*(1 - 35*pow(x0, 2)/pow(eta2, 4) + 84*pow(x0, 5.0/2.0)/pow(eta2, 5) - 70*pow(x0, 3)/pow(eta2, 6) + 20*pow(x0, 7.0/2.0)/pow(eta2, 7));
	}
}

template <typename T> void opuPaira3(T *u, T *xij, T *qi, T *qj, int *ti, int *tj, int *ai, int *aj, T *mu, T *eta, int *kappa, int dim, int ncq, int nmu, int neta, int nkappa, int ng)
{
	for (int i = 0; i <ng; i**) {
		T mu5 = mu[4];
		T mu6 = mu[5];
		T mu7 = mu[6];
		T mu8 = mu[7];
		T mu9 = mu[8];
		T eta3 = eta[2];
		T xij1 = xij[0 + i*3];
		T xij2 = xij[1 + i*3];
		T xij3 = xij[2 + i*3];
		T x0 = pow(xij1, 2) + pow(xij2, 2) + pow(xij3, 2);
		T x1 = mu7*(1 - pow(pow(x0, (1.0/2.0)*mu9) + pow(1.0/mu5, mu9), 1.0/mu9)/mu8);
		f[0 + i*1] = mu6*(-2*exp(0.5*x1) + exp(x1))*(1 - 35*pow(x0, 2)/pow(eta3, 4) + 84*pow(x0, 5.0/2.0)/pow(eta3, 5) - 70*pow(x0, 3)/pow(eta3, 6) + 20*pow(x0, 7.0/2.0)/pow(eta3, 7));
	}
}

template <typename T> void opuPaira(T *u, T *xij, T *qi, T *qj, int *ti, int *tj, int *ai, int *aj, T *mu, T *eta, int *kappa, int dim, int ncq, int nmu, int neta, int nkappa, int ng, int potnum)
{
	if (potnum == 1)
		opuPaira1(u, xij, qi, qj, ti, tj, ai, aj, mu, eta, kappa, dim, ncq, nmu, neta, nkappa, ng);
	else if (potnum == 2)
		opuPaira2(u, xij, qi, qj, ti, tj, ai, aj, mu, eta, kappa, dim, ncq, nmu, neta, nkappa, ng);
	else if (potnum == 3)
		opuPaira3(u, xij, qi, qj, ti, tj, ai, aj, mu, eta, kappa, dim, ncq, nmu, neta, nkappa, ng);
}
template void opuPaira(double *, double *, double *, double *, int *, int *, int *, int *, double *, double *, int*, int, int, int, int, int, int, int);
template void opuPaira(float *, float *, float *, float *, int *, int *, int *, int *, float *, float *, int *, int, int, int, int, int, int, int);
