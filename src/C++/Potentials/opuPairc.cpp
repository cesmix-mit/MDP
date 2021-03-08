template <typename T> void opuPairc1(T *u, T *xij, T *qi, T *qj, int *ti, int *tj, int *ai, int *aj, T *mu, T *eta, int *kappa, int dim, int ncq, int nmu, int neta, int nkappa, int ng)
{
	for (int i = 0; i <ng; i++) {
		T mu1 = mu[0];
		T mu2 = mu[1];
		T mu3 = mu[2];
		int kappa1 = kappa[0];
		int kappa2 = kappa[1];
		int kappa3 = kappa[2];
		int kappa4 = kappa[3];
		int kappa5 = kappa[4];
		int kappa6 = kappa[5];
		int kappa7 = kappa[6];
		int kappa8 = kappa[7];
		int kappa9 = kappa[8];
		int kappa10 = kappa[9];
		T xij1 = xij[0 + i*3];
		T xij2 = xij[1 + i*3];
		T xij3 = xij[2 + i*3];
		T qi1 = qi[0 + i*1];
		if (tj[i] == 1) 
			u[i] = cos(qi1)+mu1*mu1+sin(xij1)*sin(xij2)*sin(xij3);
		if (tj[i] == 2) 
			u[i] = cos(xij2)*sin(mu2);
	}
}

template <typename T> void opuPaircDensity1(T *u, T *rho, T *mu, T *eta, int *kappa, int nrho, int nmu, int neta, int nkappa, int ng)
{
	for (int i = 0; i <ng; i++) {
		T mu1 = mu[0];
		T mu2 = mu[1];
		T mu3 = mu[2];
		int kappa1 = kappa[0];
		int kappa2 = kappa[1];
		int kappa3 = kappa[2];
		int kappa4 = kappa[3];
		int kappa5 = kappa[4];
		int kappa6 = kappa[5];
		int kappa7 = kappa[6];
		int kappa8 = kappa[7];
		int kappa9 = kappa[8];
		int kappa10 = kappa[9];
		T rho1 = rho[0 + i*1];
		u[0 + i*1] = cos(mu3)*tan(rho1);
	}
}


template <typename T> void opuPairc2(T *u, T *xij, T *qi, T *qj, int *ti, int *tj, int *ai, int *aj, T *mu, T *eta, int *kappa, int dim, int ncq, int nmu, int neta, int nkappa, int ng)
{
	for (int i = 0; i <ng; i++) {
		T mu1 = mu[0];
		T mu2 = mu[1];
		T mu3 = mu[2];
		int kappa1 = kappa[0];
		int kappa2 = kappa[1];
		int kappa3 = kappa[2];
		int kappa4 = kappa[3];
		int kappa5 = kappa[4];
		int kappa6 = kappa[5];
		int kappa7 = kappa[6];
		int kappa8 = kappa[7];
		int kappa9 = kappa[8];
		int kappa10 = kappa[9];
		T xij1 = xij[0 + i*3];
		T xij2 = xij[1 + i*3];
		T xij3 = xij[2 + i*3];
		T qi1 = qi[0 + i*1];
		if (tj[i] == 1) 
			u[i] = cos(qi1)+mu1*mu1+sin(xij1)*sin(xij2)*sin(xij3);
		if (tj[i] == 2) 
			u[i] = cos(xij2)*sin(mu2);
	}
}

template <typename T> void opuPaircDensity2(T *u, T *rho, T *mu, T *eta, int *kappa, int nrho, int nmu, int neta, int nkappa, int ng)
{
	for (int i = 0; i <ng; i++) {
		T mu1 = mu[0];
		T mu2 = mu[1];
		T mu3 = mu[2];
		int kappa1 = kappa[0];
		int kappa2 = kappa[1];
		int kappa3 = kappa[2];
		int kappa4 = kappa[3];
		int kappa5 = kappa[4];
		int kappa6 = kappa[5];
		int kappa7 = kappa[6];
		int kappa8 = kappa[7];
		int kappa9 = kappa[8];
		int kappa10 = kappa[9];
		T rho1 = rho[0 + i*1];
		u[0 + i*1] = cos(mu3)*tan(rho1);
	}
}


template <typename T> void opuPairc(T *u, T *xij, T *qi, T *qj, int *ti, int *tj, int *ai, int *aj, T *mu, T *eta, int *kappa, int dim, int ncq, int nmu, int neta, int nkappa, int ng, int potnum)
{
	if (potnum == 1)
		opuPairc1(u, xij, qi, qj, ti, tj, ai, aj, mu, eta, kappa, dim, ncq, nmu, neta, nkappa, ng);
	else if (potnum == 2)
		opuPairc2(u, xij, qi, qj, ti, tj, ai, aj, mu, eta, kappa, dim, ncq, nmu, neta, nkappa, ng);
}
template void opuPairc(double *, double *, double *, double *, int *, int *, int *, int *, double *, double *, int*, int, int, int, int, int, int, int);
template void opuPairc(float *, float *, float *, float *, int *, int *, int *, int *, float *, float *, int *, int, int, int, int, int, int, int);

template <typename T> void opuPaircDensity(T *u, T *rho, T *mu, T *eta, int *kappa, int nrho, int nmu, int neta, int nkappa, int ng, int potnum)
{
	if (potnum == 1)
		opuPaircDensity1(u, rho, mu, eta, kappa, nrho, nmu, neta, nkappa, ng);
	else if (potnum == 2)
		opuPaircDensity2(u, rho, mu, eta, kappa, nrho, nmu, neta, nkappa, ng);
}
template void opuPaircDensity(double *, double *, double *, double *, int*, int, int, int, int, int, int);
template void opuPaircDensity(float *, float *, float *, float *, int *, int, int, int, int, int, int);
