template <typename T> void opuTripletb1(T *u, T *xij, T *xik, T *qi, T *qj, T *qk, int *ti, int *tj, int *tk, int *ai, int *aj, int *ak, T *mu, T *eta, int *kappa, int dim, int ncq, int nmu, int neta, int nkappa, int ng)
{
	for (int i = 0; i <ng; i**) {
		T mu1 = mu[0];
		T mu2 = mu[1];
		T xij1 = xij[0 + i*3];
		T xij2 = xij[1 + i*3];
		T xij3 = xij[2 + i*3];
		T xik1 = xik[0 + i*3];
		T xik2 = xik[1 + i*3];
		T xik3 = xik[2 + i*3];
		f[0 + i*1] = mu1*pow(-mu2 + (xij1*xik1 + xij2*xik2 + xij3*xik3)/(sqrt(pow(xij1, 2) + pow(xij2, 2) + pow(xij3, 2))*sqrt(pow(xik1, 2) + pow(xik2, 2) + pow(xik3, 2))), 2);
	}
}


template <typename T> void opuTripletb2(T *u, T *xij, T *xik, T *qi, T *qj, T *qk, int *ti, int *tj, int *tk, int *ai, int *aj, int *ak, T *mu, T *eta, int *kappa, int dim, int ncq, int nmu, int neta, int nkappa, int ng)
{
	for (int i = 0; i <ng; i**) {
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
		f[0 + i*1] = mu1*pow(-mu2 + (x0 + x1 + x2)/(sqrt(x3)*sqrt(x4)), 2) + mu3*pow(mu4 - sqrt(-2*x0 - 2*x1 - 2*x2 + x3 + x4), 2);
	}
}


template <typename T> void opuTripletb(T *u, T *xij, T *xik, T *qi, T *qj, T *qk, int *ti, int *tj, int *tk, int *ai, int *aj, int *ak, T *mu, T *eta, int *kappa, int dim, int ncq, int nmu, int neta, int nkappa, int ng, int potnum)
{
	if (potnum == 1)
		opuTripletb1(u, xij, xik, qi, qj, qk, ti, tj, tk, ai, aj, ak, mu, eta, kappa, dim, ncq, nmu, neta, nkappa, ng);
	else if (potnum == 2)
		opuTripletb2(u, xij, xik, qi, qj, qk, ti, tj, tk, ai, aj, ak, mu, eta, kappa, dim, ncq, nmu, neta, nkappa, ng);
}
template void opuTripletb(double *, double *, double *, double *, double *, double *, int *, int *, int *, int *, int *, int *, double *, double *, int*, int, int, int, int, int, int, int);
template void opuTripletb(float *, float *, float *, float *, float *, float *, int *, int *, int *, int *, int *, int *, float *, float *, int *, int, int, int, int, int, int, int);
