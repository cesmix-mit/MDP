template <typename T> void opuQuadrupletb1(T *u, T *xij, T *xik, T *xil, T *qi, T *qj, T *qk, T *ql, int *ti, int *tj, int *tk, int *tl, int *ai, int *aj, int *ak, int *al, T *mu, T *eta, int *kappa, int dim, int ncq, int nmu, int neta, int nkappa, int ng)
{
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
		T t2 = -xik1;
		T t3 = -xik2;
		T t4 = -xik3;
		T t5 = -xil1;
		T t6 = -xil2;
		T t7 = -xil3;
		T t8 = t2+xij1;
		T t9 = t3+xij2;
		T t10 = t4+xij3;
		T t11 = t5+xij1;
		T t12 = t6+xij2;
		T t13 = t7+xij3;
		T t14 = t8*t12;
		T t15 = t9*t11;
		T t16 = t8*t13;
		T t17 = t10*t11;
		T t18 = t9*t13;
		T t19 = t10*t12;
		T t20 = -t15;
		T t21 = -t17;
		T t22 = -t19;
		T t23 = t14+t20;
		T t24 = t16+t21;
		T t25 = t18+t22;
		T t26 = t23*t23;
		T t27 = t24*t24;
		T t28 = t25*t25;
		T t29 = t23*xij3;
		T t30 = t24*xij2;
		T t31 = t25*xij1;
		T t32 = -t30;
		T t33 = t26+t27+t28;
		T t34 = t29+t31+t32;
		u[i] = (mu1*(t34*t34))/t33+mu2*1.0/(t33*t33)*(t34*t34*t34*t34);
	}
}


template <typename T> void opuQuadrupletb(T *u, T *xij, T *xik, T *xil, T *qi, T *qj, T *qk, T *ql, int *ti, int *tj, int *tk, int *tl, int *ai, int *aj, int *ak, int *al, T *mu, T *eta, int *kappa, int dim, int ncq, int nmu, int neta, int nkappa, int ng, int potnum)
{
	if (potnum == 1)
		opuQuadrupletb1(u, xij, xik, xil, qi, qj, qk, ql, ti, tj, tk, tl, ai, aj, ak, al, mu, eta, kappa, dim, ncq, nmu, neta, nkappa, ng);
}
template void opuQuadrupletb(double *, double *, double *, double *, double *, double *, double *, double *, int *, int *, int *, int *, int *, int *, int *, int *, double *, double *, int*, int, int, int, int, int, int, int);
template void opuQuadrupletb(float *, float *, float *, float *, float *, float *, float *, float *, int *, int *, int *, int *, int *, int *, int *, int *, float *, float *, int *, int, int, int, int, int, int, int);
