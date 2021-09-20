template <typename T>  __global__  void kernelgpuPaira1(T *u, T *xij, T *qi, T *qj, int *ti, int *tj, int *ai, int *aj, T *mu, T *eta, int *kappa, int dim, int ncq, int nmu, int neta, int nkappa, int ng)
{
	int i = threadIdx.x + blockIdx.x * blockDim.x;
	while (i<ng) {
		T xij1 = xij[0 + i*3];
		T xij2 = xij[1 + i*3];
		T xij3 = xij[2 + i*3];
		int ti1 = ti[0 + i*1];
		int tj1 = tj[0 + i*1];
		T t2 = xij1*xij1;
		T t3 = xij2*xij2;
		T t4 = xij3*xij3;
		T t5 = t2+t3+t4;
		T t6 = sqrt(t5);
		T t7 = t6-4.0;
		T t8 = t7*t7;
		T t9 = ti1*3.4E1;
		T t10 = tj1*3.4E1;
		T t11 = -t9+8.3E1;
		T t12 = pow(t11,2.3E1/1.0E2);
		T t13 = -t10+8.3E1;
		T t14 = pow(t13,2.3E1/1.0E2);
		T t29 = t12*2.868550693703308E1;
		T t30 = t14*2.868550693703308E1;
		T t15 = -t29-t30;
		T t16 = exp(t15);
		T t26 = t12*1.807479188900747;
		T t27 = t14*1.807479188900747;
		T t17 = -t26-t27;
		T t18 = exp(t17);
		T t23 = t12*3.611910352187834;
		T t24 = t14*3.611910352187834;
		T t19 = -t23-t24;
		T t20 = exp(t19);
		T t32 = t12*8.447423692636073;
		T t33 = t14*8.447423692636073;
		T t21 = -t32-t33;
		T t22 = exp(t21);
		T t57 = t12*8.599786552828175E-1;
		T t58 = t14*8.599786552828175E-1;
		T t25 = t57+t58;
		T t61 = t12*4.303521878335112E-1;
		T t62 = t14*4.303521878335112E-1;
		T t28 = t61+t62;
		T t65 = t12*6.829882604055496;
		T t66 = t14*6.829882604055496;
		T t31 = t65+t66;
		T t69 = t12*2.011291355389541;
		T t70 = t14*2.011291355389541;
		T t34 = t69+t70;
		T t35 = t9-8.3E1;
		T t36 = t10-8.3E1;
		T t37 = t12*1.241331163287086;
		T t38 = t14*1.241331163287086;
		T t39 = t37+t38;
		T t40 = t12*1.212302113127001E-2;
		T t41 = t14*1.212302113127001E-2;
		T t42 = t40+t41;
		T t43 = t12*2.409832187833511E-1;
		T t44 = t14*2.409832187833511E-1;
		T t45 = t43+t44;
		T t46 = t12*1.025477010458911;
		T t47 = t14*1.025477010458911;
		T t48 = t46+t47;
		T t49 = t20*3.177097505668934E-2;
		T t50 = t18*3.193877551020408E-3;
		T t51 = t16*2.060657596371882E-2;
		T t52 = t22*5.780725623582766E-2;
		T t53 = t16*t39*(1.0E1/2.1E1);
		T t54 = t18*t42*(1.0E1/2.1E1);
		T t55 = t20*t45*(1.0E1/2.1E1);
		T t56 = t22*t48*(1.0E1/2.1E1);
		T t59 = t25*t25;
		T t60 = t20*t59*2.8022E-1;
		T t63 = t28*t28;
		T t64 = t18*t63*2.817E-2;
		T t67 = t31*t31;
		T t68 = t16*t67*1.8175E-1;
		T t71 = t34*t34;
		T t72 = t22*t71*5.0986E-1;
		T t73 = t49+t50+t51+t52+t53+t54+t55+t56+t60+t64+t68+t72;
		T t74 = t20*6.671904761904762E-2;
		T t75 = t18*6.707142857142857E-3;
		T t76 = t16*4.327380952380952E-2;
		T t77 = t22*1.213952380952381E-1;
		T t78 = t16*t39;
		T t79 = t18*t42;
		T t80 = t20*t45;
		T t81 = t22*t48;
		T t82 = t74+t75+t76+t77+t78+t79+t80+t81;
		T t83 = t6*1.0E2;
		T t84 = t83-4.0E2;
		T t85 = tanh(t84);
		T t86 = t85*(1.0/2.0);
		T t87 = t20*2.8022E-1;
		T t88 = t18*2.817E-2;
		T t89 = t16*1.8175E-1;
		T t90 = t22*5.0986E-1;
		T t91 = t87+t88+t89+t90;
		T t92 = t35*t36*t91*3.428486904761905;
		T t93 = t35*t36*t73*1.142828968253968E-2;
		T t94 = t35*t36*t82*3.428486904761905E-1;
		u[i] = -(t86+1.0/2.0)*(t92+t93+t94-t7*t8*(-t7*(t35*t36*t73*2.14280431547619E1+t35*t36*t82*2.14280431547619E2)+t35*t36*t73*5.714144841269841+t35*t36*t82*8.571217261904762E1))+(t86-1.0/2.0)*(t92+t93+t94)+1.0/sqrt(t5)*t35*t36*(exp(-t6*t25)*2.8022E-1+exp(-t6*t28)*2.817E-2+exp(-t6*t31)*1.8175E-1+exp(-t6*t34)*5.0986E-1)*1.4399645E1;
		i += blockDim.x * gridDim.x;
	}
}

template <typename T> void gpuPaira1(T *u, T *xij, T *qi, T *qj, int *ti, int *tj, int *ai, int *aj, T *mu, T *eta, int *kappa, int dim, int ncq, int nmu, int neta, int nkappa, int ng)
{
	int blockDim = 256;
	int gridDim = (ng + blockDim - 1) / blockDim;
	gridDim = (gridDim>1024)? 1024 : gridDim;
	kernelgpuPaira1<<<gridDim, blockDim>>>(u, xij, qi, qj, ti, tj, ai, aj, mu, eta, kappa, dim, ncq, nmu, neta, nkappa, ng);
}

#ifdef _ENZYME
template <typename T>  __device__  void devicegpuPaira1(T *__restrict__ u, T *__restrict__ xij, T *__restrict__ qi, T *__restrict__ qj, int *__restrict__ ti, int *__restrict__ tj, int *__restrict__ ai, int *__restrict__ aj, T *__restrict__ mu, T *__restrict__ eta, int *__restrict__ kappa, int dim, int ncq, int nmu, int neta, int nkappa, int ng)
{
	int i = threadIdx.x + blockIdx.x * blockDim.x;
	while (i<ng) {
		T xij1 = xij[0 + i*3];
		T xij2 = xij[1 + i*3];
		T xij3 = xij[2 + i*3];
		int ti1 = ti[0 + i*1];
		int tj1 = tj[0 + i*1];
		T t2 = xij1*xij1;
		T t3 = xij2*xij2;
		T t4 = xij3*xij3;
		T t5 = t2+t3+t4;
		T t6 = sqrt(t5);
		T t7 = t6-4.0;
		T t8 = t7*t7;
		T t9 = ti1*3.4E1;
		T t10 = tj1*3.4E1;
		T t11 = -t9+8.3E1;
		T t12 = pow(t11,2.3E1/1.0E2);
		T t13 = -t10+8.3E1;
		T t14 = pow(t13,2.3E1/1.0E2);
		T t29 = t12*2.868550693703308E1;
		T t30 = t14*2.868550693703308E1;
		T t15 = -t29-t30;
		T t16 = exp(t15);
		T t26 = t12*1.807479188900747;
		T t27 = t14*1.807479188900747;
		T t17 = -t26-t27;
		T t18 = exp(t17);
		T t23 = t12*3.611910352187834;
		T t24 = t14*3.611910352187834;
		T t19 = -t23-t24;
		T t20 = exp(t19);
		T t32 = t12*8.447423692636073;
		T t33 = t14*8.447423692636073;
		T t21 = -t32-t33;
		T t22 = exp(t21);
		T t57 = t12*8.599786552828175E-1;
		T t58 = t14*8.599786552828175E-1;
		T t25 = t57+t58;
		T t61 = t12*4.303521878335112E-1;
		T t62 = t14*4.303521878335112E-1;
		T t28 = t61+t62;
		T t65 = t12*6.829882604055496;
		T t66 = t14*6.829882604055496;
		T t31 = t65+t66;
		T t69 = t12*2.011291355389541;
		T t70 = t14*2.011291355389541;
		T t34 = t69+t70;
		T t35 = t9-8.3E1;
		T t36 = t10-8.3E1;
		T t37 = t12*1.241331163287086;
		T t38 = t14*1.241331163287086;
		T t39 = t37+t38;
		T t40 = t12*1.212302113127001E-2;
		T t41 = t14*1.212302113127001E-2;
		T t42 = t40+t41;
		T t43 = t12*2.409832187833511E-1;
		T t44 = t14*2.409832187833511E-1;
		T t45 = t43+t44;
		T t46 = t12*1.025477010458911;
		T t47 = t14*1.025477010458911;
		T t48 = t46+t47;
		T t49 = t20*3.177097505668934E-2;
		T t50 = t18*3.193877551020408E-3;
		T t51 = t16*2.060657596371882E-2;
		T t52 = t22*5.780725623582766E-2;
		T t53 = t16*t39*(1.0E1/2.1E1);
		T t54 = t18*t42*(1.0E1/2.1E1);
		T t55 = t20*t45*(1.0E1/2.1E1);
		T t56 = t22*t48*(1.0E1/2.1E1);
		T t59 = t25*t25;
		T t60 = t20*t59*2.8022E-1;
		T t63 = t28*t28;
		T t64 = t18*t63*2.817E-2;
		T t67 = t31*t31;
		T t68 = t16*t67*1.8175E-1;
		T t71 = t34*t34;
		T t72 = t22*t71*5.0986E-1;
		T t73 = t49+t50+t51+t52+t53+t54+t55+t56+t60+t64+t68+t72;
		T t74 = t20*6.671904761904762E-2;
		T t75 = t18*6.707142857142857E-3;
		T t76 = t16*4.327380952380952E-2;
		T t77 = t22*1.213952380952381E-1;
		T t78 = t16*t39;
		T t79 = t18*t42;
		T t80 = t20*t45;
		T t81 = t22*t48;
		T t82 = t74+t75+t76+t77+t78+t79+t80+t81;
		T t83 = t6*1.0E2;
		T t84 = t83-4.0E2;
		T t85 = tanh(t84);
		T t86 = t85*(1.0/2.0);
		T t87 = t20*2.8022E-1;
		T t88 = t18*2.817E-2;
		T t89 = t16*1.8175E-1;
		T t90 = t22*5.0986E-1;
		T t91 = t87+t88+t89+t90;
		T t92 = t35*t36*t91*3.428486904761905;
		T t93 = t35*t36*t73*1.142828968253968E-2;
		T t94 = t35*t36*t82*3.428486904761905E-1;
		u[i] = -(t86+1.0/2.0)*(t92+t93+t94-t7*t8*(-t7*(t35*t36*t73*2.14280431547619E1+t35*t36*t82*2.14280431547619E2)+t35*t36*t73*5.714144841269841+t35*t36*t82*8.571217261904762E1))+(t86-1.0/2.0)*(t92+t93+t94)+1.0/sqrt(t5)*t35*t36*(exp(-t6*t25)*2.8022E-1+exp(-t6*t28)*2.817E-2+exp(-t6*t31)*1.8175E-1+exp(-t6*t34)*5.0986E-1)*1.4399645E1;
		i += blockDim.x * gridDim.x;
	}
}


template <typename T>  __global__  void kernelgpuPaira1Gradient(T *__restrict__ u, T *__restrict__ du, T *__restrict__ u_xij, T *__restrict__ xij, T *__restrict__ qi, T *__restrict__ qj, int *__restrict__ ti, int *__restrict__ tj, int *__restrict__ ai, int *__restrict__ aj, T *__restrict__ mu, T *__restrict__ eta, int *__restrict__ kappa, int dim, int ncq, int nmu, int neta, int nkappa, int ng)
{
	__enzyme_autodiff((void*)devicegpuPaira1<T>, 
		enzyme_dup, u, du, 
		enzyme_dup, xij, u_xij, 
		enzyme_const, qi, 
		enzyme_const, qj, 
		enzyme_const, ti, 
		enzyme_const, tj, 
		enzyme_const, ai, 
		enzyme_const, aj, 
		enzyme_const, mu, 
		enzyme_const, eta, 
		enzyme_const, kappa, 
		dim, ncq, nmu, neta, nkappa, ng); 
}

template <typename T> void gpuPaira1Gradient(T *u, T *du, T *u_xij, T *xij, T *qi, T *qj, int *ti, int *tj, int *ai, int *aj, T *mu, T *eta, int *kappa, int dim, int ncq, int nmu, int neta, int nkappa, int ng)
{
	int blockDim = 256;
	int gridDim = (ng + blockDim - 1) / blockDim;
	gridDim = (gridDim>1024)? 1024 : gridDim;
	kernelgpuPaira1Gradient<<<gridDim, blockDim>>>(u, du, u_xij, xij, qi, qj, ti, tj, ai, aj, mu, eta, kappa, dim, ncq, nmu, neta, nkappa, ng);
}
#endif

template <typename T> void gpuPaira(T *u, T *xij, T *qi, T *qj, int *ti, int *tj, int *ai, int *aj, T *mu, T *eta, int *kappa, int dim, int ncq, int nmu, int neta, int nkappa, int ng, int potnum)
{
	if (potnum == 1)
		gpuPaira1(u, xij, qi, qj, ti, tj, ai, aj, mu, eta, kappa, dim, ncq, nmu, neta, nkappa, ng);
}
template void gpuPaira(double *, double *, double *, double *, int *, int *, int *, int *, double *, double *, int*, int, int, int, int, int, int, int);
template void gpuPaira(float *, float *, float *, float *, int *, int *, int *, int *, float *, float *, int *, int, int, int, int, int, int, int);

#ifdef _ENZYME
template <typename T> void gpuPairaGradient(T *u, T *du, T *u_xij, T *xij, T *qi, T *qj, int *ti, int *tj, int *ai, int *aj, T *mu, T *eta, int *kappa, int dim, int ncq, int nmu, int neta, int nkappa, int ng, int potnum)
{
	if (potnum == 1)
		gpuPaira1Gradient(u, du, u_xij, xij, qi, qj, ti, tj, ai, aj, mu, eta, kappa, dim, ncq, nmu, neta, nkappa, ng);
}
template void gpuPairaGradient(double *, double *, double*, double *, double *, double *, int *, int *, int *, int *, double *, double *, int*, int, int, int, int, int, int, int);
template void gpuPairaGradient(float *, float *, float*, float *, float *, float *, int *, int *, int *, int *, float *, float *, int *, int, int, int, int, int, int, int);
#endif
