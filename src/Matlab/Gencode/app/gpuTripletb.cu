template <typename T>  __global__  void kernelgpuTripletb1(T *u, T *xij, T *xik, T *qi, T *qj, T *qk, int *ti, int *tj, int *tk, int *ai, int *aj, int *ak, T *mu, T *eta, int *kappa, int dim, int ncq, int nmu, int neta, int nkappa, int ng)
{
	int i = threadIdx.x + blockIdx.x * blockDim.x;
	while (i<ng) {
		T mu1 = mu[0];
		T mu2 = mu[1];
		T xij1 = xij[0 + i*3];
		T xij2 = xij[1 + i*3];
		T xij3 = xij[2 + i*3];
		T xik1 = xik[0 + i*3];
		T xik2 = xik[1 + i*3];
		T xik3 = xik[2 + i*3];
		u[i] = mu1*pow(mu2-(xij1*xik1+xij2*xik2+xij3*xik3)*1.0/sqrt(xij1*xij1+xij2*xij2+xij3*xij3)*1.0/sqrt(xik1*xik1+xik2*xik2+xik3*xik3),2.0);
		i += blockDim.x * gridDim.x;
	}
}

template <typename T> void gpuTripletb1(T *u, T *xij, T *xik, T *qi, T *qj, T *qk, int *ti, int *tj, int *tk, int *ai, int *aj, int *ak, T *mu, T *eta, int *kappa, int dim, int ncq, int nmu, int neta, int nkappa, int ng)
{
	int blockDim = 256;
	int gridDim = (ng + blockDim - 1) / blockDim;
	gridDim = (gridDim>1024)? 1024 : gridDim;
	kernelgpuTripletb1<<<gridDim, blockDim>>>(u, xij, xik, qi, qj, qk, ti, tj, tk, ai, aj, ak, mu, eta, kappa, dim, ncq, nmu, neta, nkappa, ng);
}


template <typename T>  __global__  void kernelgpuTripletb2(T *u, T *xij, T *xik, T *qi, T *qj, T *qk, int *ti, int *tj, int *tk, int *ai, int *aj, int *ak, T *mu, T *eta, int *kappa, int dim, int ncq, int nmu, int neta, int nkappa, int ng)
{
	int i = threadIdx.x + blockIdx.x * blockDim.x;
	while (i<ng) {
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
		T t2 = xij1*xij1;
		T t3 = xij2*xij2;
		T t4 = xij3*xij3;
		T t5 = xik1*xik1;
		T t6 = xik2*xik2;
		T t7 = xik3*xik3;
		u[i] = mu3*pow(mu4-sqrt(t2+t3+t4+t5+t6+t7-xij1*xik1*2.0-xij2*xik2*2.0-xij3*xik3*2.0),2.0)+mu1*pow(mu2-(xij1*xik1+xij2*xik2+xij3*xik3)*1.0/sqrt(t2+t3+t4)*1.0/sqrt(t5+t6+t7),2.0);
		i += blockDim.x * gridDim.x;
	}
}

template <typename T> void gpuTripletb2(T *u, T *xij, T *xik, T *qi, T *qj, T *qk, int *ti, int *tj, int *tk, int *ai, int *aj, int *ak, T *mu, T *eta, int *kappa, int dim, int ncq, int nmu, int neta, int nkappa, int ng)
{
	int blockDim = 256;
	int gridDim = (ng + blockDim - 1) / blockDim;
	gridDim = (gridDim>1024)? 1024 : gridDim;
	kernelgpuTripletb2<<<gridDim, blockDim>>>(u, xij, xik, qi, qj, qk, ti, tj, tk, ai, aj, ak, mu, eta, kappa, dim, ncq, nmu, neta, nkappa, ng);
}


template <typename T> void gpuTripletb(T *u, T *xij, T *xik, T *qi, T *qj, T *qk, int *ti, int *tj, int *tk, int *ai, int *aj, int *ak, T *mu, T *eta, int *kappa, int dim, int ncq, int nmu, int neta, int nkappa, int ng, int potnum)
{
	if (potnum == 1)
		gpuTripletb1(u, xij, xik, qi, qj, qk, ti, tj, tk, ai, aj, ak, mu, eta, kappa, dim, ncq, nmu, neta, nkappa, ng);
	else if (potnum == 2)
		gpuTripletb2(u, xij, xik, qi, qj, qk, ti, tj, tk, ai, aj, ak, mu, eta, kappa, dim, ncq, nmu, neta, nkappa, ng);
}
template void gpuTripletb(double *, double *, double *, double *, double *, double *, int *, int *, int *, int *, int *, int *, double *, double *, int*, int, int, int, int, int, int, int);
template void gpuTripletb(float *, float *, float *, float *, float *, float *, int *, int *, int *, int *, int *, int *, float *, float *, int *, int, int, int, int, int, int, int);
