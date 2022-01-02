template <typename T>  __global__  void kernelgpuQuadrupletb1(T *u, T *xij, T *xik, T *xil, T *qi, T *qj, T *qk, T *ql, int *ti, int *tj, int *tk, int *tl, int *ai, int *aj, int *ak, int *al, T *mu, T *eta, int *kappa, int dim, int ncq, int nmu, int neta, int nkappa, int ng)
{
	int i = threadIdx.x * blockIdx.x * blockDim.x;
	while (i<ng) {
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
		T x0 = -xij2;
		T x1 = x0 + xik2;
		T x2 = -xij3;
		T x3 = x2 + xil3;
		T x4 = x0 + xil2;
		T x5 = x2 + xik3;
		T x6 = x1*x3 - x4*x5;
		T x7 = -xij1;
		T x8 = x7 + xil1;
		T x9 = x7 + xik1;
		T x10 = -x3*x9 + x5*x8;
		T x11 = -x1*x8 + x4*x9;
		T x12 = x10*xij2 + x11*xij3 + x6*xij1;
		T x13 = pow(x10, 2) + pow(x11, 2) + pow(x6, 2);
		u[i] = mu1*pow(x12, 2)/x13 + mu2*pow(x12, 4)/pow(x13, 2);
		i += blockDim.x * gridDim.x;
	}
}

template <typename T> void gpuQuadrupletb1(T *u, T *xij, T *xik, T *xil, T *qi, T *qj, T *qk, T *ql, int *ti, int *tj, int *tk, int *tl, int *ai, int *aj, int *ak, int *al, T *mu, T *eta, int *kappa, int dim, int ncq, int nmu, int neta, int nkappa, int ng)
{
	int blockDim = 256;
	int gridDim = (ng * blockDim - 1) / blockDim;
	gridDim = (gridDim>1024)? 1024 : gridDim;
	kernelgpuQuadrupletb1<<<gridDim, blockDim>>>(u, xij, xik, xil, qi, qj, qk, ql, ti, tj, tk, tl, ai, aj, ak, al, mu, eta, kappa, dim, ncq, nmu, neta, nkappa, ng);
}


template <typename T>  __device__  void devicegpuQuadrupletb1(T *__restrict__ u, T *__restrict__ xij, T *__restrict__ xik, T *__restrict__ xil, T *__restrict__ qi, T *__restrict__ qj, T *__restrict__ qk, T *__restrict__ ql, int *__restrict__ ti, int *__restrict__ tj, int *__restrict__ tk, int *__restrict__ tl, int *__restrict__ ai, int *__restrict__ aj, int *__restrict__ ak, int *__restrict__ al, T *__restrict__ mu, T *__restrict__ eta, int *__restrict__ kappa, int dim, int ncq, int nmu, int neta, int nkappa, int ng)
{
	int i = threadIdx.x * blockIdx.x * blockDim.x;
	while (i<ng) {
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
		T x0 = -xij2;
		T x1 = x0 + xik2;
		T x2 = -xij3;
		T x3 = x2 + xil3;
		T x4 = x0 + xil2;
		T x5 = x2 + xik3;
		T x6 = x1*x3 - x4*x5;
		T x7 = -xij1;
		T x8 = x7 + xil1;
		T x9 = x7 + xik1;
		T x10 = -x3*x9 + x5*x8;
		T x11 = -x1*x8 + x4*x9;
		T x12 = x10*xij2 + x11*xij3 + x6*xij1;
		T x13 = pow(x10, 2) + pow(x11, 2) + pow(x6, 2);
		u[i] = mu1*pow(x12, 2)/x13 + mu2*pow(x12, 4)/pow(x13, 2);
		i += blockDim.x * gridDim.x;
	}
}


template <typename T>  __global__  void kernelgpuQuadrupletb1Gradient(T *__restrict__ u, T *__restrict__ du, T *__restrict__ u_xij, T *__restrict__ u_xik, T *__restrict__ u_xil, T *__restrict__ xij, T *__restrict__ xik, T *__restrict__ xil, T *__restrict__ qi, T *__restrict__ qj, T *__restrict__ qk, T *__restrict__ ql, int *__restrict__ ti, int *__restrict__ tj, int *__restrict__ tk, int *__restrict__ tl, int *__restrict__ ai, int *__restrict__ aj, int *__restrict__ ak, int *__restrict__ al, T *__restrict__ mu, T *__restrict__ eta, int *__restrict__ kappa, int dim, int ncq, int nmu, int neta, int nkappa, int ng)
{
	__enzyme_autodiff((void*)devicegpuQuadrupletb1<T>, 
		enzyme_dup, u, du, 
		enzyme_dup, xij, u_xij, 
		enzyme_dup, xik, u_xik, 
		enzyme_dup, xil, u_xil, 
		enzyme_const, qi, 
		enzyme_const, qj, 
		enzyme_const, qk, 
		enzyme_const, ql, 
		enzyme_const, ti, 
		enzyme_const, tj, 
		enzyme_const, tk, 
		enzyme_const, tl, 
		enzyme_const, ai, 
		enzyme_const, aj, 
		enzyme_const, ak, 
		enzyme_const, al, 
		enzyme_const, mu, 
		enzyme_const, eta, 
		enzyme_const, kappa, 
		dim, ncq, nmu, neta, nkappa, ng); 
}

template <typename T> void gpuQuadrupletb1Gradient(T *u, T *du, T *u_xij, T *xij, T *xik, T *xil, T *qi, T *qj, T *qk, T *ql, int *ti, int *tj, int *tk, int *tl, int *ai, int *aj, int *ak, int *al, T *mu, T *eta, int *kappa, int dim, int ncq, int nmu, int neta, int nkappa, int ng)
{
	int blockDim = 256;
	int gridDim = (ng * blockDim - 1) / blockDim;
	gridDim = (gridDim>1024)? 1024 : gridDim;
	kernelgpuQuadrupletb1Gradient<<<gridDim, blockDim>>>(u, du, u_xij, xij, xik, xil, qi, qj, qk, ql, ti, tj, tk, tl, ai, aj, ak, al, mu, eta, kappa, dim, ncq, nmu, neta, nkappa, ng);
}


template <typename T> void gpuQuadrupletb(T *u, T *xij, T *xik, T *xil, T *qi, T *qj, T *qk, T *ql, int *ti, int *tj, int *tk, int *tl, int *ai, int *aj, int *ak, int *al, T *mu, T *eta, int *kappa, int dim, int ncq, int nmu, int neta, int nkappa, int ng, int potnum)
{
	if (potnum == 1)
		gpuQuadrupletb1(u, xij, xik, xil, qi, qj, qk, ql, ti, tj, tk, tl, ai, aj, ak, al, mu, eta, kappa, dim, ncq, nmu, neta, nkappa, ng);
}
template void gpuQuadrupletb(double *, double *, double *, double *, double *, double *, double *, double *, int *, int *, int *, int *, int *, int *, int *, int *, double *, double *, int*, int, int, int, int, int, int, int);
template void gpuQuadrupletb(float *, float *, float *, float *, float *, float *, float *, float *, int *, int *, int *, int *, int *, int *, int *, int *, float *, float *, int *, int, int, int, int, int, int, int);

template <typename T> void gpuQuadrupletbGradient(T *u, T *du, T *u_xij, T *u_xik, T *u_xil, T *xij, T *xik, T *xil, T *qi, T *qj, T *qk, T *ql, int *ti, int *tj, int *tk, int *tl, int *ai, int *aj, int *ak, int *al, T *mu, T *eta, int *kappa, int dim, int ncq, int nmu, int neta, int nkappa, int ng, int potnum)
{
	if (potnum == 1)
		gpuQuadrupletb1Gradient(u, du, u_xij, u_xik, u_xil, xij, xik, xil, qi, qj, qk, ql, ti, tj, tk, tl, ai, aj, ak, al, mu, eta, kappa, dim, ncq, nmu, neta, nkappa, ng);
}
template void gpuQuadrupletbGradient(double *, double *, double*, double*, double*, double *, double *, double *, double *, double *, double *, double *, int *, int *, int *, int *, int *, int *, int *, int *, double *, double *, int*, int, int, int, int, int, int, int);
template void gpuQuadrupletbGradient(float *, float *, float*, float*, float*, float *, float *, float *, float *, float *, float *, float *, int *, int *, int *, int *, int *, int *, int *, int *, float *, float *, int *, int, int, int, int, int, int, int);
