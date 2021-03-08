template <typename T> void opuTripletc1(T *u, T *xij, T *xik, T *qi, T *qj, T *qk, int *ti, int *tj, int *tk, int *ai, int *aj, int *ak, T *mu, T *eta, int *kappa, int dim, int ncq, int nmu, int neta, int nkappa, int ng)
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
		if (tk[i] == 1) 
			u[i] = sin(xij1)*sin(xij2)*sin(xij3);
		if (tk[i] == 2) 
			u[i] = cos(qi1)+mu1*mu1;
	}
}


template <typename T> void opuTripletcPair1(T *u, T *xij, T *qi, T *qj, int *ti, int *tj, int *ai, int *aj, T *mu, T *eta, int *kappa, int dim, int ncq, int nmu, int neta, int nkappa, int ng)
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
		u[0 + i*1] = cos(xij2)*sin(mu2);
	}
}

template <typename T> void opuTripletcDensity1(T *u, T *rho, T *mu, T *eta, int *kappa, int nrho, int nmu, int neta, int nkappa, int ng)
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


template <typename T> void opuTripletc2(T *u, T *xij, T *xik, T *qi, T *qj, T *qk, int *ti, int *tj, int *tk, int *ai, int *aj, int *ak, T *mu, T *eta, int *kappa, int dim, int ncq, int nmu, int neta, int nkappa, int ng)
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
		if (tk[i] == 1) 
			u[i] = sin(xij1)*sin(xij2)*sin(xij3);
		if (tk[i] == 2) 
			u[i] = cos(qi1)+mu1*mu1;
	}
}


template <typename T> void opuTripletcPair2(T *u, T *xij, T *qi, T *qj, int *ti, int *tj, int *ai, int *aj, T *mu, T *eta, int *kappa, int dim, int ncq, int nmu, int neta, int nkappa, int ng)
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
		u[0 + i*1] = cos(xij2)*sin(mu2);
	}
}

template <typename T> void opuTripletcDensity2(T *u, T *rho, T *mu, T *eta, int *kappa, int nrho, int nmu, int neta, int nkappa, int ng)
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


template <typename T> void opuTripletc3(T *u, T *xij, T *xik, T *qi, T *qj, T *qk, int *ti, int *tj, int *tk, int *ai, int *aj, int *ak, T *mu, T *eta, int *kappa, int dim, int ncq, int nmu, int neta, int nkappa, int ng)
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
		if (tk[i] == 1) 
			u[i] = sin(xij1)*sin(xij2)*sin(xij3);
		if (tk[i] == 2) 
			u[i] = cos(qi1)+mu1*mu1;
	}
}


template <typename T> void opuTripletcPair3(T *u, T *xij, T *qi, T *qj, int *ti, int *tj, int *ai, int *aj, T *mu, T *eta, int *kappa, int dim, int ncq, int nmu, int neta, int nkappa, int ng)
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
		u[0 + i*1] = cos(xij2)*sin(mu2);
	}
}

template <typename T> void opuTripletcDensity3(T *u, T *rho, T *mu, T *eta, int *kappa, int nrho, int nmu, int neta, int nkappa, int ng)
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


template <typename T> void opuTripletc(T *u, T *xij, T *xik, T *qi, T *qj, T *qk, int *ti, int *tj, int *tk, int *ai, int *aj, int *ak, T *mu, T *eta, int *kappa, int dim, int ncq, int nmu, int neta, int nkappa, int ng, int potnum)
{
	if (potnum == 1)
		opuTripletc1(u, xij, xik, qi, qj, qk, ti, tj, tk, ai, aj, ak, mu, eta, kappa, dim, ncq, nmu, neta, nkappa, ng);
	else if (potnum == 2)
		opuTripletc2(u, xij, xik, qi, qj, qk, ti, tj, tk, ai, aj, ak, mu, eta, kappa, dim, ncq, nmu, neta, nkappa, ng);
	else if (potnum == 3)
		opuTripletc3(u, xij, xik, qi, qj, qk, ti, tj, tk, ai, aj, ak, mu, eta, kappa, dim, ncq, nmu, neta, nkappa, ng);
}
template void opuTripletc(double *, double *, double *, double *, double *, double *, int *, int *, int *, int *, int *, int *, double *, double *, int*, int, int, int, int, int, int, int);
template void opuTripletc(float *, float *, float *, float *, float *, float *, int *, int *, int *, int *, int *, int *, float *, float *, int *, int, int, int, int, int, int, int);

template <typename T> void opuTripletcPair(T *u, T *xij, T *qi, T *qj, int *ti, int *tj, int *ai, int *aj, T *mu, T *eta, int *kappa, int dim, int ncq, int nmu, int neta, int nkappa, int ng, int potnum)
{
	if (potnum == 1)
		opuTripletcPair1(u, xij, qi, qj, ti, tj, ai, aj, mu, eta, kappa, dim, ncq, nmu, neta, nkappa, ng);
	else if (potnum == 2)
		opuTripletcPair2(u, xij, qi, qj, ti, tj, ai, aj, mu, eta, kappa, dim, ncq, nmu, neta, nkappa, ng);
	else if (potnum == 3)
		opuTripletcPair3(u, xij, qi, qj, ti, tj, ai, aj, mu, eta, kappa, dim, ncq, nmu, neta, nkappa, ng);
}
template void opuTripletcPair(double *, double *, double *, double *, int *, int *, int *, int *, double *, double *, int*, int, int, int, int, int, int, int);
template void opuTripletcPair(float *, float *, float *, float *, int *, int *, int *, int *, float *, float *, int *, int, int, int, int, int, int, int);

template <typename T> void opuTripletcDensity(T *u, T *rho, T *mu, T *eta, int *kappa, int nrho, int nmu, int neta, int nkappa, int ng, int potnum)
{
	if (potnum == 1)
		opuTripletcDensity1(u, rho, mu, eta, kappa, nrho, nmu, neta, nkappa, ng);
	else if (potnum == 2)
		opuTripletcDensity2(u, rho, mu, eta, kappa, nrho, nmu, neta, nkappa, ng);
	else if (potnum == 3)
		opuTripletcDensity3(u, rho, mu, eta, kappa, nrho, nmu, neta, nkappa, ng);
}
template void opuTripletcDensity(double *, double *, double *, double *, int*, int, int, int, int, int, int);
template void opuTripletcDensity(float *, float *, float *, float *, int *, int, int, int, int, int, int);
