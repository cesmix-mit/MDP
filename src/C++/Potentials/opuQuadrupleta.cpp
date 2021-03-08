template <typename T> void opuQuadrupleta1(T *u, T *xij, T *xik, T *xil, T *qi, T *qj, T *qk, T *ql, int *ti, int *tj, int *tk, int *tl, int *ai, int *aj, int *ak, int *al, T *mu, T *eta, int *kappa, int dim, int ncq, int nmu, int neta, int nkappa, int ng)
{
	for (int i = 0; i <ng; i++) {
		T mu1 = mu[0];
		T mu2 = mu[1];
		T mu3 = mu[2];
		T mu4 = mu[3];
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
		u[i] = cos(qi1)+mu1*mu1+sin(xij1)*sin(xij2)*sin(xij3);
	}
}


template <typename T> void opuQuadrupleta2(T *u, T *xij, T *xik, T *xil, T *qi, T *qj, T *qk, T *ql, int *ti, int *tj, int *tk, int *tl, int *ai, int *aj, int *ak, int *al, T *mu, T *eta, int *kappa, int dim, int ncq, int nmu, int neta, int nkappa, int ng)
{
	for (int i = 0; i <ng; i++) {
		T mu1 = mu[0];
		T mu2 = mu[1];
		T mu3 = mu[2];
		T mu4 = mu[3];
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
		T xik1 = xik[0 + i*3];
		T xik2 = xik[1 + i*3];
		T xik3 = xik[2 + i*3];
		u[i] = cos(xik2)*sin(mu2);
	}
}


template <typename T> void opuQuadrupleta3(T *u, T *xij, T *xik, T *xil, T *qi, T *qj, T *qk, T *ql, int *ti, int *tj, int *tk, int *tl, int *ai, int *aj, int *ak, int *al, T *mu, T *eta, int *kappa, int dim, int ncq, int nmu, int neta, int nkappa, int ng)
{
	for (int i = 0; i <ng; i++) {
		T mu1 = mu[0];
		T mu2 = mu[1];
		T mu3 = mu[2];
		T mu4 = mu[3];
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
		T xil1 = xil[0 + i*3];
		T xil2 = xil[1 + i*3];
		T xil3 = xil[2 + i*3];
		u[i] = cos(mu3)*tan(xil2);
	}
}


template <typename T> void opuQuadrupleta(T *u, T *xij, T *xik, T *xil, T *qi, T *qj, T *qk, T *ql, int *ti, int *tj, int *tk, int *tl, int *ai, int *aj, int *ak, int *al, T *mu, T *eta, int *kappa, int dim, int ncq, int nmu, int neta, int nkappa, int ng, int potnum)
{
	if (potnum == 1)
		opuQuadrupleta1(u, xij, xik, xil, qi, qj, qk, ql, ti, tj, tk, tl, ai, aj, ak, al, mu, eta, kappa, dim, ncq, nmu, neta, nkappa, ng);
	else if (potnum == 2)
		opuQuadrupleta2(u, xij, xik, xil, qi, qj, qk, ql, ti, tj, tk, tl, ai, aj, ak, al, mu, eta, kappa, dim, ncq, nmu, neta, nkappa, ng);
	else if (potnum == 3)
		opuQuadrupleta3(u, xij, xik, xil, qi, qj, qk, ql, ti, tj, tk, tl, ai, aj, ak, al, mu, eta, kappa, dim, ncq, nmu, neta, nkappa, ng);
}
template void opuQuadrupleta(double *, double *, double *, double *, double *, double *, double *, double *, int *, int *, int *, int *, int *, int *, int *, int *, double *, double *, int*, int, int, int, int, int, int, int);
template void opuQuadrupleta(float *, float *, float *, float *, float *, float *, float *, float *, int *, int *, int *, int *, int *, int *, int *, int *, float *, float *, int *, int, int, int, int, int, int, int);
