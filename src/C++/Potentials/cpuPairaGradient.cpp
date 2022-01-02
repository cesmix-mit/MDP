template <typename T> void cpuPairaGradient1(T *u, T *u_xij, T *xij, T *qi, T *qj, int *ti, int *tj, int *ai, int *aj, T *mu, T *eta, int *kappa, int dim, int ncq, int nmu, int neta, int nkappa, int ng)
{
	#pragma omp parallel for
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
		T t12 = 1.0/sqrt(t5);
		T t13 = t6*3.428291787990583E-2;
		T t14 = t13-1.827782636402304E-1;
		T t15 = t9*(1.0/2.0);
		T t16 = t15+1.0/2.0;
		T t21 = t6*2.308969885292557;
		T t17 = exp(-t21);
		T t23 = t6*4.614046060829141;
		T t18 = exp(-t23);
		T t25 = t6*3.664438963872198E1;
		T t19 = exp(-t25);
		T t27 = t6*1.079118754693148E1;
		T t20 = exp(-t27);
		T t22 = t17*2.16164490013485E3;
		T t24 = t18*2.15028801532051E4;
		T t26 = t19*1.394671496625875E4;
		T t28 = t20*3.91244681854013E4;
		T t29 = t22+t24+t26+t28;
		T t30 = t9*t9;
		T t31 = t30-1.0;
		T t32 = t10*t11*t14;
		T t33 = t32+1.625187736600612E-2;
		T t34 = 1.0/pow(t5,3.0/2.0);
		u[i] = t9*8.125938683003058E-3+t12*t29-t16*t33-8.125938683003058E-3;
		u_xij[0 + i*3] = -t12*(t12*t17*xij1*4.991172977107606E3+t12*t18*xij1*9.921527946737711E4+t12*t19*xij1*5.11068857403781E5+t12*t20*xij1*4.221994738626192E5)-t16*(t10*t11*t12*xij1*3.428291787990583E-2+t11*t12*t14*xij1*3.0)-t12*t31*xij1*8.125938683003058E-1-t29*t34*xij1+t12*t31*t33*xij1*5.0E1;
		u_xij[1 + i*3] = -t12*(t12*t17*xij2*4.991172977107606E3+t12*t18*xij2*9.921527946737711E4+t12*t19*xij2*5.11068857403781E5+t12*t20*xij2*4.221994738626192E5)-t16*(t10*t11*t12*xij2*3.428291787990583E-2+t11*t12*t14*xij2*3.0)-t12*t31*xij2*8.125938683003058E-1-t29*t34*xij2+t12*t31*t33*xij2*5.0E1;
		u_xij[2 + i*3] = -t12*(t12*t17*xij3*4.991172977107606E3+t12*t18*xij3*9.921527946737711E4+t12*t19*xij3*5.11068857403781E5+t12*t20*xij3*4.221994738626192E5)-t16*(t10*t11*t12*xij3*3.428291787990583E-2+t11*t12*t14*xij3*3.0)-t12*t31*xij3*8.125938683003058E-1-t29*t34*xij3+t12*t31*t33*xij3*5.0E1;
	}
}


template <typename T> void cpuPairaGradient(T *u, T *u_xij, T *xij, T *qi, T *qj, int *ti, int *tj, int *ai, int *aj, T *mu, T *eta, int *kappa, int dim, int ncq, int nmu, int neta, int nkappa, int ng, int potnum)
{
	if (potnum == 1)
		cpuPairaGradient1(u, u_xij, xij, qi, qj, ti, tj, ai, aj, mu, eta, kappa, dim, ncq, nmu, neta, nkappa, ng);
}
template void cpuPairaGradient(double *, double *, double *, double *, double *, int *, int *, int *, int *, double *, double *, int*, int, int, int, int, int, int, int);
template void cpuPairaGradient(float *, float *, float *, float *, float *, int *, int *, int *, int *, float *, float *, int *, int, int, int, int, int, int, int);
