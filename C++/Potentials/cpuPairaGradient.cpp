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
		T t7 = 1.0/t6;
		T t9 = t6-4.0;
		T t10 = t6*1.0E+2;
		T t19 = t6*3.664438963872198E+1;
		T t20 = t6*1.079118754693148E+1;
		T t23 = t6*2.308969885292557;
		T t24 = t6*4.614046060829141;
		T t31 = t6*3.428291787990583E-2;
		T t8 = t7*t7*t7;
		T t11 = t9*t9;
		T t12 = t9*t9*t9;
		T t13 = t9*1.0E+2;
		T t21 = -t19;
		T t22 = -t20;
		T t27 = -t23;
		T t28 = -t24;
		T t35 = t31-1.827782636402304E-1;
		T t14 = tanh(t13);
		T t25 = exp(t21);
		T t26 = exp(t22);
		T t29 = exp(t27);
		T t30 = exp(t28);
		T t37 = t12*t35;
		T t15 = t14*t14;
		T t16 = t14/2.0;
		T t32 = t25*1.394671496625875E+4;
		T t33 = t26*3.91244681854013E+4;
		T t34 = t29*2.16164490013485E+3;
		T t36 = t30*2.15028801532051E+4;
		T t38 = t37+1.625187736600612E-2;
		T t17 = t15-1.0;
		T t18 = t16+1.0/2.0;
		T t39 = t32+t33+t34+t36;
		u[i] = t14*8.125938683003058E-3+t7*t39-t18*t38-8.125938683003058E-3;
		u_xij[0 + i*3] = -t7*(t7*t25*xij1*5.11068857403781E+5+t7*t26*xij1*4.221994738626192E+5+t7*t29*xij1*4.991172977107606E+3+t7*t30*xij1*9.921527946737711E+4)-t18*(t7*t12*xij1*3.428291787990583E-2+t7*t11*t35*xij1*3.0)-t7*t17*xij1*8.125938683003058E-1-t8*t39*xij1+t7*t17*t38*xij1*5.0E+1;
		u_xij[1 + i*3] = -t7*(t7*t25*xij2*5.11068857403781E+5+t7*t26*xij2*4.221994738626192E+5+t7*t29*xij2*4.991172977107606E+3+t7*t30*xij2*9.921527946737711E+4)-t18*(t7*t12*xij2*3.428291787990583E-2+t7*t11*t35*xij2*3.0)-t7*t17*xij2*8.125938683003058E-1-t8*t39*xij2+t7*t17*t38*xij2*5.0E+1;
		u_xij[2 + i*3] = -t7*(t7*t25*xij3*5.11068857403781E+5+t7*t26*xij3*4.221994738626192E+5+t7*t29*xij3*4.991172977107606E+3+t7*t30*xij3*9.921527946737711E+4)-t18*(t7*t12*xij3*3.428291787990583E-2+t7*t11*t35*xij3*3.0)-t7*t17*xij3*8.125938683003058E-1-t8*t39*xij3+t7*t17*t38*xij3*5.0E+1;
	}
}


template <typename T> void cpuPairaGradient(T *u, T *u_xij, T *xij, T *qi, T *qj, int *ti, int *tj, int *ai, int *aj, T *mu, T *eta, int *kappa, int dim, int ncq, int nmu, int neta, int nkappa, int ng, int potnum)
{
	if (potnum == 1)
		cpuPairaGradient1(u, u_xij, xij, qi, qj, ti, tj, ai, aj, mu, eta, kappa, dim, ncq, nmu, neta, nkappa, ng);
}
template void cpuPairaGradient(double *, double *, double *, double *, double *, int *, int *, int *, int *, double *, double *, int*, int, int, int, int, int, int, int);
template void cpuPairaGradient(float *, float *, float *, float *, float *, int *, int *, int *, int *, float *, float *, int *, int, int, int, int, int, int, int);
