template <typename T>  __global__  void kernelgpuPairb1(T *u, T *xij, T *qi, T *qj, int *ti, int *tj, int *ai, int *aj, T *mu, T *eta, int *kappa, int dim, int ncq, int nmu, int neta, int nkappa, int ng)
{
	int i = threadIdx.x + blockIdx.x * blockDim.x;
	while (i<ng) {
		T mu1 = mu[0];
		T mu2 = mu[1];
		T mu3 = mu[2];
		T xij1 = xij[0 + i*3];
		T xij2 = xij[1 + i*3];
		T xij3 = xij[2 + i*3];
		u[i] = -mu1*exp(mu2*(mu3-sqrt(xij1*xij1+xij2*xij2+xij3*xij3))*2.0);
		i += blockDim.x * gridDim.x;
	}
}

template <typename T> void gpuPairb1(T *u, T *xij, T *qi, T *qj, int *ti, int *tj, int *ai, int *aj, T *mu, T *eta, int *kappa, int dim, int ncq, int nmu, int neta, int nkappa, int ng)
{
	int blockDim = 256;
	int gridDim = (ng + blockDim - 1) / blockDim;
	gridDim = (gridDim>1024)? 1024 : gridDim;
	kernelgpuPairb1<<<gridDim, blockDim>>>(u, xij, qi, qj, ti, tj, ai, aj, mu, eta, kappa, dim, ncq, nmu, neta, nkappa, ng);
}

template <typename T>  __global__  void kernelgpuPairb2(T *u, T *xij, T *qi, T *qj, int *ti, int *tj, int *ai, int *aj, T *mu, T *eta, int *kappa, int dim, int ncq, int nmu, int neta, int nkappa, int ng)
{
	int i = threadIdx.x + blockIdx.x * blockDim.x;
	while (i<ng) {
		T mu4 = mu[3];
		T mu5 = mu[4];
		T mu6 = mu[5];
		T xij1 = xij[0 + i*3];
		T xij2 = xij[1 + i*3];
		T xij3 = xij[2 + i*3];
		u[i] = mu4*(pow(mu5-sqrt(xij1*xij1+xij2*xij2+xij3*xij3),2.0)-pow(mu5-mu6,2.0));
		i += blockDim.x * gridDim.x;
	}
}

template <typename T> void gpuPairb2(T *u, T *xij, T *qi, T *qj, int *ti, int *tj, int *ai, int *aj, T *mu, T *eta, int *kappa, int dim, int ncq, int nmu, int neta, int nkappa, int ng)
{
	int blockDim = 256;
	int gridDim = (ng + blockDim - 1) / blockDim;
	gridDim = (gridDim>1024)? 1024 : gridDim;
	kernelgpuPairb2<<<gridDim, blockDim>>>(u, xij, qi, qj, ti, tj, ai, aj, mu, eta, kappa, dim, ncq, nmu, neta, nkappa, ng);
}

template <typename T> void gpuPairb(T *u, T *xij, T *qi, T *qj, int *ti, int *tj, int *ai, int *aj, T *mu, T *eta, int *kappa, int dim, int ncq, int nmu, int neta, int nkappa, int ng, int potnum)
{
	if (potnum == 1)
		gpuPairb1(u, xij, qi, qj, ti, tj, ai, aj, mu, eta, kappa, dim, ncq, nmu, neta, nkappa, ng);
	else if (potnum == 2)
		gpuPairb2(u, xij, qi, qj, ti, tj, ai, aj, mu, eta, kappa, dim, ncq, nmu, neta, nkappa, ng);
}
template void gpuPairb(double *, double *, double *, double *, int *, int *, int *, int *, double *, double *, int*, int, int, int, int, int, int, int);
template void gpuPairb(float *, float *, float *, float *, int *, int *, int *, int *, float *, float *, int *, int, int, int, int, int, int, int);
