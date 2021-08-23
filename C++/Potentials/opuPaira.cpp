template <typename T> void opuPaira1(T *u, T *xij, T *qi, T *qj, int *ti, int *tj, int *ai, int *aj, T *mu, T *eta, int *kappa, int dim, int ncq, int nmu, int neta, int nkappa, int ng)
{
	for (int i = 0; i <ng; i++) {
		T xij1 = xij[0 + i*3];
		T xij2 = xij[1 + i*3];
		T xij3 = xij[2 + i*3];
		T t2 = xij1*xij1;
		T t3 = xij2*xij2;
		T t4 = xij3*xij3;
		T t5 = t2+t3+t4;
		T t6 = sqrt(t5);
		T t7 = t6*1.0E2;
		T t8 = t7-4.0E2;
		T t9 = tanh(t8);
		T t10 = t6-4.0;
		T t11 = t10*t10;
		u[i] = t9*8.125938683003058E-3-(t9*(1.0/2.0)+1.0/2.0)*(t10*t11*(t6*3.428291787990583E-2-1.827782636402304E-1)+1.625187736600612E-2)+1.0/sqrt(t5)*(exp(t6*(-2.308969885292557))*2.16164490013485E3+exp(t6*(-4.614046060829141))*2.15028801532051E4+exp(t6*(-3.664438963872198E1))*1.394671496625875E4+exp(t6*(-1.079118754693148E1))*3.91244681854013E4)-8.125938683003058E-3;
	}
}

template <typename T> void opuPaira(T *u, T *xij, T *qi, T *qj, int *ti, int *tj, int *ai, int *aj, T *mu, T *eta, int *kappa, int dim, int ncq, int nmu, int neta, int nkappa, int ng, int potnum)
{
	if (potnum == 1)
		opuPaira1(u, xij, qi, qj, ti, tj, ai, aj, mu, eta, kappa, dim, ncq, nmu, neta, nkappa, ng);
}
template void opuPaira(double *, double *, double *, double *, int *, int *, int *, int *, double *, double *, int*, int, int, int, int, int, int, int);
template void opuPaira(float *, float *, float *, float *, int *, int *, int *, int *, float *, float *, int *, int, int, int, int, int, int, int);
