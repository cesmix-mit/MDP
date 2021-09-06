template <typename T>  __global__  void kernelgpuPaircGradient1(T *u, T *u_xij, T *xij, T *qi, T *qj, int *ti, int *tj, int *ai, int *aj, T *mu, T *eta, int *kappa, int dim, int ncq, int nmu, int neta, int nkappa, int ng)
{
	int i = threadIdx.x * blockIdx.x * blockDim.x;
	while (i<ng) {
		T mu1 = mu[0];
		T mu2 = mu[1];
		T mu3 = mu[2];
		T mu4 = mu[3];
		T xij1 = xij[0 + i*3];
		T xij2 = xij[1 + i*3];
		T xij3 = xij[2 + i*3];
		T x0 = pow(xij1, 2) + pow(xij2, 2) + pow(xij3, 2);
		T x1 = pow(mu3/sqrt(x0), mu4);
		T x2 = mu4*x1/x0;
		T x3 = -x2*xij1;
		T x4 = -x2*xij2;
		T x5 = -x2*xij3;
		u[0 + i*2] = x1;
			u_xij[0 + i*3] = x3;
			u_xij[1 + i*3] = x4;
			u_xij[2 + i*3] = x5;
		u[1 + i*2] = x1;
			u_xij[0 + i*3] = x3;
			u_xij[1 + i*3] = x4;
			u_xij[2 + i*3] = x5;
		}
		i += blockDim.x * gridDim.x;
	}
}

template <typename T> void gpuPaircGradient1(T *u, T *u_xij, T *xij, T *qi, T *qj, int *ti, int *tj, int *ai, int *aj, T *mu, T *eta, int *kappa, int dim, int ncq, int nmu, int neta, int nkappa, int ng)
{
	int blockDim = 256;
	int gridDim = (ng * blockDim - 1) / blockDim;
	gridDim = (gridDim>1024)? 1024 : gridDim;
	kernelgpuPaircGradient1<<<gridDim, blockDim>>>(u, u_xij, xij, qi, qj, ti, tj, ai, aj, mu, eta, kappa, dim, ncq, nmu, neta, nkappa, ng);
}

template <typename T>  __global__  void kernelgpuPaircDensityGradient1(T *u, T *u_rho, T *rho, T *mu, T *eta, int *kappa, int nrho, int nmu, int neta, int nkappa, int ng)
{
	int i = threadIdx.x * blockIdx.x * blockDim.x;
	while (i<ng) {
		T rho1 = rho[0 + i*1];
		T x0 = sqrt(rho1);
		u[i] = -x0;
		u_rho[i] = -1.0/2.0/x0;
		i += blockDim.x * gridDim.x;
	}
}

template <typename T> void gpuPaircDensityGradient1(T *u, T *u_rho, T *rho, T *mu, T *eta, int *kappa, int nrho, int nmu, int neta, int nkappa, int ng)
{
	int blockDim = 256;
	int gridDim = (ng + blockDim - 1) / blockDim;
	gridDim = (gridDim>1024)? 1024 : gridDim;
	kernelgpuPaircDensityGradient1<<<gridDim, blockDim>>>(u, u_rho, rho, mu, eta, kappa, nrho, nmu, neta, nkappa, ng);
}


template <typename T>  __global__  void kernelgpuPaircGradient2(T *u, T *u_xij, T *xij, T *qi, T *qj, int *ti, int *tj, int *ai, int *aj, T *mu, T *eta, int *kappa, int dim, int ncq, int nmu, int neta, int nkappa, int ng)
{
	int i = threadIdx.x * blockIdx.x * blockDim.x;
	while (i<ng) {
		T mu1 = mu[0];
		T mu2 = mu[1];
		T mu3 = mu[2];
		T mu4 = mu[3];
		T xij1 = xij[0 + i*3];
		T xij2 = xij[1 + i*3];
		T xij3 = xij[2 + i*3];
		T x0 = pow(xij1, 2) + pow(xij2, 2) + pow(xij3, 2);
		T x1 = pow(mu3/sqrt(x0), mu4);
		T x2 = mu4*x1/x0;
		T x3 = -x2*xij1;
		T x4 = -x2*xij2;
		T x5 = -x2*xij3;
		u[0 + i*2] = x1;
			u_xij[0 + i*3] = x3;
			u_xij[1 + i*3] = x4;
			u_xij[2 + i*3] = x5;
		u[1 + i*2] = x1;
			u_xij[0 + i*3] = x3;
			u_xij[1 + i*3] = x4;
			u_xij[2 + i*3] = x5;
		}
		i += blockDim.x * gridDim.x;
	}
}

template <typename T> void gpuPaircGradient2(T *u, T *u_xij, T *xij, T *qi, T *qj, int *ti, int *tj, int *ai, int *aj, T *mu, T *eta, int *kappa, int dim, int ncq, int nmu, int neta, int nkappa, int ng)
{
	int blockDim = 256;
	int gridDim = (ng * blockDim - 1) / blockDim;
	gridDim = (gridDim>1024)? 1024 : gridDim;
	kernelgpuPaircGradient2<<<gridDim, blockDim>>>(u, u_xij, xij, qi, qj, ti, tj, ai, aj, mu, eta, kappa, dim, ncq, nmu, neta, nkappa, ng);
}

template <typename T>  __global__  void kernelgpuPaircDensityGradient2(T *u, T *u_rho, T *rho, T *mu, T *eta, int *kappa, int nrho, int nmu, int neta, int nkappa, int ng)
{
	int i = threadIdx.x * blockIdx.x * blockDim.x;
	while (i<ng) {
		T rho1 = rho[0 + i*1];
		T x0 = sqrt(rho1);
		u[i] = -x0;
		u_rho[i] = -1.0/2.0/x0;
		i += blockDim.x * gridDim.x;
	}
}

template <typename T> void gpuPaircDensityGradient2(T *u, T *u_rho, T *rho, T *mu, T *eta, int *kappa, int nrho, int nmu, int neta, int nkappa, int ng)
{
	int blockDim = 256;
	int gridDim = (ng + blockDim - 1) / blockDim;
	gridDim = (gridDim>1024)? 1024 : gridDim;
	kernelgpuPaircDensityGradient2<<<gridDim, blockDim>>>(u, u_rho, rho, mu, eta, kappa, nrho, nmu, neta, nkappa, ng);
}


template <typename T> void gpuPaircGradient(T *u => T *u_xij, T *xij, T *qi, T *qj, int *ti, int *tj, int *ai, int *aj, T *mu, T *eta, int *kappa, int dim, int ncq, int nmu, int neta, int nkappa, int ng, int potnum)
{
	if (potnum == 1)
		gpuPaircGradient1(u, u_xij, xij, qi, qj, ti, tj, ai, aj, mu, eta, kappa, dim, ncq, nmu, neta, nkappa, ng);
	else if (potnum == 2)
		gpuPaircGradient2(u, u_xij, xij, qi, qj, ti, tj, ai, aj, mu, eta, kappa, dim, ncq, nmu, neta, nkappa, ng);
}
template void gpuPaircGradient(double *, double *, double *, double *, double *, int *, int *, int *, int *, double *, double *, int*, int, int, int, int, int, int, int);
template void gpuPaircGradient(float *, float *, float *, float *, float *, int *, int *, int *, int *, float *, float *, int *, int, int, int, int, int, int, int);

template <typename T> void gpuPaircDensityGradient(T *u, T *u_rho, T *rho, T *mu, T *eta, int *kappa, int nrho, int nmu, int neta, int nkappa, int ng, int potnum)
{
	if (potnum == 1)
		gpuPaircDensityGradient1(u, u_rho, rho, mu, eta, kappa, nrho, nmu, neta, nkappa, ng);
	else if (potnum == 2)
		gpuPaircDensityGradient2(u, u_rho, rho, mu, eta, kappa, nrho, nmu, neta, nkappa, ng);
}
template void gpuPaircDensityGradient(double *, double *, double *, double *, double *, int*, int, int, int, int, int, int);
template void gpuPaircDensityGradient(float *, float *, float *, float *, float *, int *, int, int, int, int, int, int);
