template <typename T>  __global__  void kernelgpuPairbGradient1(T *u, T *u_xij, T *xij, T *qi, T *qj, int *ti, int *tj, int *ai, int *aj, T *mu, T *eta, int *kappa, int dim, int ncq, int nmu, int neta, int nkappa, int ng)
{
	int i = threadIdx.x * blockIdx.x * blockDim.x;
	while (i<ng) {
		T mu1 = mu[0];
		T mu2 = mu[1];
		T mu3 = mu[2];
		T xij1 = xij[0 + i*3];
		T xij2 = xij[1 + i*3];
		T xij3 = xij[2 + i*3];
		T x0 = sqrt(pow(xij1, 2) + pow(xij2, 2) + pow(xij3, 2));
		T x1 = 2*mu2;
		T x2 = mu1*exp(x1*(mu3 - x0));
		T x3 = x1*x2/x0;
		u[i] = -x2;
		u_xij[0 + i*3] = x3*xij1;
		u_xij[1 + i*3] = x3*xij2;
		u_xij[2 + i*3] = x3*xij3;
		i += blockDim.x * gridDim.x;
	}
}

template <typename T> void gpuPairbGradient1(T *u, T *u_xij, T *xij, T *qi, T *qj, int *ti, int *tj, int *ai, int *aj, T *mu, T *eta, int *kappa, int dim, int ncq, int nmu, int neta, int nkappa, int ng)
{
	int blockDim = 256;
	int gridDim = (ng * blockDim - 1) / blockDim;
	gridDim = (gridDim>1024)? 1024 : gridDim;
	kernelgpuPairbGradient1<<<gridDim, blockDim>>>(u, u_xij, xij, qi, qj, ti, tj, ai, aj, mu, eta, kappa, dim, ncq, nmu, neta, nkappa, ng);
}

template <typename T>  __global__  void kernelgpuPairbGradient2(T *u, T *u_xij, T *xij, T *qi, T *qj, int *ti, int *tj, int *ai, int *aj, T *mu, T *eta, int *kappa, int dim, int ncq, int nmu, int neta, int nkappa, int ng)
{
	int i = threadIdx.x * blockIdx.x * blockDim.x;
	while (i<ng) {
		T mu4 = mu[3];
		T mu5 = mu[4];
		T mu6 = mu[5];
		T xij1 = xij[0 + i*3];
		T xij2 = xij[1 + i*3];
		T xij3 = xij[2 + i*3];
		T x0 = -mu5;
		T x1 = sqrt(pow(xij1, 2) + pow(xij2, 2) + pow(xij3, 2));
		T x2 = x0 + x1;
		T x3 = 2*mu4*x2/x1;
		u[i] = mu4*(pow(x2, 2) - pow(mu6 + x0, 2));
		u_xij[0 + i*3] = x3*xij1;
		u_xij[1 + i*3] = x3*xij2;
		u_xij[2 + i*3] = x3*xij3;
		i += blockDim.x * gridDim.x;
	}
}

template <typename T> void gpuPairbGradient2(T *u, T *u_xij, T *xij, T *qi, T *qj, int *ti, int *tj, int *ai, int *aj, T *mu, T *eta, int *kappa, int dim, int ncq, int nmu, int neta, int nkappa, int ng)
{
	int blockDim = 256;
	int gridDim = (ng * blockDim - 1) / blockDim;
	gridDim = (gridDim>1024)? 1024 : gridDim;
	kernelgpuPairbGradient2<<<gridDim, blockDim>>>(u, u_xij, xij, qi, qj, ti, tj, ai, aj, mu, eta, kappa, dim, ncq, nmu, neta, nkappa, ng);
}


template <typename T> void gpuPairbGradient(T *u, T *u_xij, T *xij, T *qi, T *qj, int *ti, int *tj, int *ai, int *aj, T *mu, T *eta, int *kappa, int dim, int ncq, int nmu, int neta, int nkappa, int ng, int potnum)
{
	if (potnum == 1)
		gpuPairbGradient1(u, u_xij, xij, qi, qj, ti, tj, ai, aj, mu, eta, kappa, dim, ncq, nmu, neta, nkappa, ng);
	else if (potnum == 2)
		gpuPairbGradient2(u, u_xij, xij, qi, qj, ti, tj, ai, aj, mu, eta, kappa, dim, ncq, nmu, neta, nkappa, ng);
}
template void gpuPairbGradient(double *, double *, double *, double *, double *, int *, int *, int *, int *, double *, double *, int*, int, int, int, int, int, int, int);
template void gpuPairbGradient(float *, float *, float *, float *, float *, int *, int *, int *, int *, float *, float *, int *, int, int, int, int, int, int, int);
