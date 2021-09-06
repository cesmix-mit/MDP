template <typename T>  __global__  void kernelgpuPairc1(T *u, T *xij, T *qi, T *qj, int *ti, int *tj, int *ai, int *aj, T *mu, T *eta, int *kappa, int dim, int ncq, int nmu, int neta, int nkappa, int ng)
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
		T x0 = pow(pow(xij1, 2) + pow(xij2, 2) + pow(xij3, 2), -1.0/2.0);
		u[0 + i*2] = pow(mu1*x0, mu2);
		u[1 + i*2] = pow(mu3*x0, mu4);
		i += blockDim.x * gridDim.x;
	}
}

template <typename T> void gpuPairc1(T *u, T *xij, T *qi, T *qj, int *ti, int *tj, int *ai, int *aj, T *mu, T *eta, int *kappa, int dim, int ncq, int nmu, int neta, int nkappa, int ng)
{
	int blockDim = 256;
	int gridDim = (ng * blockDim - 1) / blockDim;
	gridDim = (gridDim>1024)? 1024 : gridDim;
	kernelgpuPairc1<<<gridDim, blockDim>>>(u, xij, qi, qj, ti, tj, ai, aj, mu, eta, kappa, dim, ncq, nmu, neta, nkappa, ng);
}

template <typename T>  __device__  void devicegpuPairc1(T *__restrict__ u, T *__restrict__ xij, T *__restrict__ qi, T *__restrict__ qj, int *__restrict__ ti, int *__restrict__ tj, int *__restrict__ ai, int *__restrict__ aj, T *__restrict__ mu, T *__restrict__ eta, int *__restrict__ kappa, int dim, int ncq, int nmu, int neta, int nkappa, int ng)
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
		T x0 = pow(pow(xij1, 2) + pow(xij2, 2) + pow(xij3, 2), -1.0/2.0);
		u[0 + i*2] = pow(mu1*x0, mu2);
		u[1 + i*2] = pow(mu3*x0, mu4);
		i += blockDim.x * gridDim.x;
	}
}


template <typename T>  __global__  void kernelgpuPairc1Gradient(T *__restrict__ u, T *__restrict__ du, T *__restrict__ u_xij, T *__restrict__ xij, T *__restrict__ qi, T *__restrict__ qj, int *__restrict__ ti, int *__restrict__ tj, int *__restrict__ ai, int *__restrict__ aj, T *__restrict__ mu, T *__restrict__ eta, int *__restrict__ kappa, int dim, int ncq, int nmu, int neta, int nkappa, int ng)
{
	__enzyme_autodiff((void*)devicegpuPairc1<T>, 
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

template <typename T> void gpuPairc1Gradient(T *u, T *du, T *u_xij, T *xij, T *qi, T *qj, int *ti, int *tj, int *ai, int *aj, T *mu, T *eta, int *kappa, int dim, int ncq, int nmu, int neta, int nkappa, int ng)
{
	int blockDim = 256;
	int gridDim = (ng * blockDim - 1) / blockDim;
	gridDim = (gridDim>1024)? 1024 : gridDim;
	kernelgpuPairc1Gradient<<<gridDim, blockDim>>>(u, du, u_xij, xij, qi, qj, ti, tj, ai, aj, mu, eta, kappa, dim, ncq, nmu, neta, nkappa, ng);
}

template <typename T>  __global__  void kernelgpuPaircDensity1(T *u, T *rho, T *mu, T *eta, int *kappa, int nrho, int nmu, int neta, int nkappa, int ng)
{
	int i = threadIdx.x * blockIdx.x * blockDim.x;
	while (i<ng) {
		T rho1 = rho[0 + i*1];
		u[0 + i*1] = -sqrt(rho1);
		i += blockDim.x * gridDim.x;
	}
}

template <typename T> void gpuPaircDensity1(T *u, T *rho, T *mu, T *eta, int *kappa, int nrho, int nmu, int neta, int nkappa, int ng)
{
	int blockDim = 256;
	int gridDim = (ng * blockDim - 1) / blockDim;
	gridDim = (gridDim>1024)? 1024 : gridDim;
	kernelgpuPaircDensity1<<<gridDim, blockDim>>>(u, rho, mu, eta, kappa, nrho, nmu, neta, nkappa, ng);
}


template <typename T>  __device__  void devicegpuPaircDensity1(T *__restrict__ u, T *__restrict__ rho, T *__restrict__ mu, T *__restrict__ eta, int *__restrict__ kappa, int nrho, int nmu, int neta, int nkappa, int ng)
{
	int i = threadIdx.x * blockIdx.x * blockDim.x;
	while (i<ng) {
		T rho1 = rho[0 + i*1];
		u[0 + i*1] = -sqrt(rho1);
		i += blockDim.x * gridDim.x;
	}
}


template <typename T>  __global__  void kernelgpuPaircDensity1Gradient(T *__restrict__ u, T *__restrict__ du, T *__restrict__ u_rho, T *__restrict__ rho, T *__restrict__ mu, T *__restrict__ eta, int *__restrict__ kappa, int nrho, int nmu, int neta, int nkappa, int ng)
{
	__enzyme_autodiff((void*)devicegpuPaircDensity1<T>, 
		enzyme_dup, u, du, 
		enzyme_dup, rho, u_rho, 
		enzyme_const, mu, 
		enzyme_const, eta, 
		enzyme_const, kappa, 
		dim, nrho, nmu, neta, nkappa, ng); 
}

template <typename T> void gpuPaircDensity1Gradient(T *u, T *du, T *u_rho, T *rho, T *mu, T *eta, int *kappa, int nrho, int nmu, int neta, int nkappa, int ng)
{
	int blockDim = 256;
	int gridDim = (ng * blockDim - 1) / blockDim;
	gridDim = (gridDim>1024)? 1024 : gridDim;
	kernelgpuPaircDensity1Gradient<<<gridDim, blockDim>>>(u, du, u_rho, rho, mu, eta, kappa, nrho, nmu, neta, nkappa, ng);
}


template <typename T>  __global__  void kernelgpuPairc2(T *u, T *xij, T *qi, T *qj, int *ti, int *tj, int *ai, int *aj, T *mu, T *eta, int *kappa, int dim, int ncq, int nmu, int neta, int nkappa, int ng)
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
		T x0 = pow(pow(xij1, 2) + pow(xij2, 2) + pow(xij3, 2), -1.0/2.0);
		u[0 + i*2] = pow(mu1*x0, mu2);
		u[1 + i*2] = pow(mu3*x0, mu4);
		i += blockDim.x * gridDim.x;
	}
}

template <typename T> void gpuPairc2(T *u, T *xij, T *qi, T *qj, int *ti, int *tj, int *ai, int *aj, T *mu, T *eta, int *kappa, int dim, int ncq, int nmu, int neta, int nkappa, int ng)
{
	int blockDim = 256;
	int gridDim = (ng * blockDim - 1) / blockDim;
	gridDim = (gridDim>1024)? 1024 : gridDim;
	kernelgpuPairc2<<<gridDim, blockDim>>>(u, xij, qi, qj, ti, tj, ai, aj, mu, eta, kappa, dim, ncq, nmu, neta, nkappa, ng);
}

template <typename T>  __device__  void devicegpuPairc2(T *__restrict__ u, T *__restrict__ xij, T *__restrict__ qi, T *__restrict__ qj, int *__restrict__ ti, int *__restrict__ tj, int *__restrict__ ai, int *__restrict__ aj, T *__restrict__ mu, T *__restrict__ eta, int *__restrict__ kappa, int dim, int ncq, int nmu, int neta, int nkappa, int ng)
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
		T x0 = pow(pow(xij1, 2) + pow(xij2, 2) + pow(xij3, 2), -1.0/2.0);
		u[0 + i*2] = pow(mu1*x0, mu2);
		u[1 + i*2] = pow(mu3*x0, mu4);
		i += blockDim.x * gridDim.x;
	}
}


template <typename T>  __global__  void kernelgpuPairc2Gradient(T *__restrict__ u, T *__restrict__ du, T *__restrict__ u_xij, T *__restrict__ xij, T *__restrict__ qi, T *__restrict__ qj, int *__restrict__ ti, int *__restrict__ tj, int *__restrict__ ai, int *__restrict__ aj, T *__restrict__ mu, T *__restrict__ eta, int *__restrict__ kappa, int dim, int ncq, int nmu, int neta, int nkappa, int ng)
{
	__enzyme_autodiff((void*)devicegpuPairc2<T>, 
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

template <typename T> void gpuPairc2Gradient(T *u, T *du, T *u_xij, T *xij, T *qi, T *qj, int *ti, int *tj, int *ai, int *aj, T *mu, T *eta, int *kappa, int dim, int ncq, int nmu, int neta, int nkappa, int ng)
{
	int blockDim = 256;
	int gridDim = (ng * blockDim - 1) / blockDim;
	gridDim = (gridDim>1024)? 1024 : gridDim;
	kernelgpuPairc2Gradient<<<gridDim, blockDim>>>(u, du, u_xij, xij, qi, qj, ti, tj, ai, aj, mu, eta, kappa, dim, ncq, nmu, neta, nkappa, ng);
}

template <typename T>  __global__  void kernelgpuPaircDensity2(T *u, T *rho, T *mu, T *eta, int *kappa, int nrho, int nmu, int neta, int nkappa, int ng)
{
	int i = threadIdx.x * blockIdx.x * blockDim.x;
	while (i<ng) {
		T rho1 = rho[0 + i*1];
		u[0 + i*1] = -sqrt(rho1);
		i += blockDim.x * gridDim.x;
	}
}

template <typename T> void gpuPaircDensity2(T *u, T *rho, T *mu, T *eta, int *kappa, int nrho, int nmu, int neta, int nkappa, int ng)
{
	int blockDim = 256;
	int gridDim = (ng * blockDim - 1) / blockDim;
	gridDim = (gridDim>1024)? 1024 : gridDim;
	kernelgpuPaircDensity2<<<gridDim, blockDim>>>(u, rho, mu, eta, kappa, nrho, nmu, neta, nkappa, ng);
}


template <typename T>  __device__  void devicegpuPaircDensity2(T *__restrict__ u, T *__restrict__ rho, T *__restrict__ mu, T *__restrict__ eta, int *__restrict__ kappa, int nrho, int nmu, int neta, int nkappa, int ng)
{
	int i = threadIdx.x * blockIdx.x * blockDim.x;
	while (i<ng) {
		T rho1 = rho[0 + i*1];
		u[0 + i*1] = -sqrt(rho1);
		i += blockDim.x * gridDim.x;
	}
}


template <typename T>  __global__  void kernelgpuPaircDensity2Gradient(T *__restrict__ u, T *__restrict__ du, T *__restrict__ u_rho, T *__restrict__ rho, T *__restrict__ mu, T *__restrict__ eta, int *__restrict__ kappa, int nrho, int nmu, int neta, int nkappa, int ng)
{
	__enzyme_autodiff((void*)devicegpuPaircDensity2<T>, 
		enzyme_dup, u, du, 
		enzyme_dup, rho, u_rho, 
		enzyme_const, mu, 
		enzyme_const, eta, 
		enzyme_const, kappa, 
		dim, nrho, nmu, neta, nkappa, ng); 
}

template <typename T> void gpuPaircDensity2Gradient(T *u, T *du, T *u_rho, T *rho, T *mu, T *eta, int *kappa, int nrho, int nmu, int neta, int nkappa, int ng)
{
	int blockDim = 256;
	int gridDim = (ng * blockDim - 1) / blockDim;
	gridDim = (gridDim>1024)? 1024 : gridDim;
	kernelgpuPaircDensity2Gradient<<<gridDim, blockDim>>>(u, du, u_rho, rho, mu, eta, kappa, nrho, nmu, neta, nkappa, ng);
}


template <typename T> void gpuPairc(T *u, T *xij, T *qi, T *qj, int *ti, int *tj, int *ai, int *aj, T *mu, T *eta, int *kappa, int dim, int ncq, int nmu, int neta, int nkappa, int ng, int potnum)
{
	if (potnum == 1)
		gpuPairc1(u, xij, qi, qj, ti, tj, ai, aj, mu, eta, kappa, dim, ncq, nmu, neta, nkappa, ng);
	else if (potnum == 2)
		gpuPairc2(u, xij, qi, qj, ti, tj, ai, aj, mu, eta, kappa, dim, ncq, nmu, neta, nkappa, ng);
}
template void gpuPairc(double *, double *, double *, double *, int *, int *, int *, int *, double *, double *, int*, int, int, int, int, int, int, int);
template void gpuPairc(float *, float *, float *, float *, int *, int *, int *, int *, float *, float *, int *, int, int, int, int, int, int, int);


template <typename T> void gpuPaircGradient(T *u, T *du, T *u_xij, T *xij, T *qi, T *qj, int *ti, int *tj, int *ai, int *aj, T *mu, T *eta, int *kappa, int dim, int ncq, int nmu, int neta, int nkappa, int ng, int potnum)
{
	if (potnum == 1)
		gpuPairc1Gradient(u, du, u_xij, xij, qi, qj, ti, tj, ai, aj, mu, eta, kappa, dim, ncq, nmu, neta, nkappa, ng);
	else if (potnum == 2)
		gpuPairc2Gradient(u, du, u_xij, xij, qi, qj, ti, tj, ai, aj, mu, eta, kappa, dim, ncq, nmu, neta, nkappa, ng);
}
template void gpuPaircGradient(double *, double *, double*, double *, double *, double *, int *, int *, int *, int *, double *, double *, int*, int, int, int, int, int, int, int);
template void gpuPaircGradient(float *, float *, float*, float *, float *, float *, int *, int *, int *, int *, float *, float *, int *, int, int, int, int, int, int, int);
template <typename T> void gpuPaircDensity(T *u, T *rho, T *mu, T *eta, int *kappa, int nrho, int nmu, int neta, int nkappa, int ng, int potnum)
{
	if (potnum == 1)
		gpuPaircDensity1(u, rho, mu, eta, kappa, nrho, nmu, neta, nkappa, ng);
	else if (potnum == 2)
		gpuPaircDensity2(u, rho, mu, eta, kappa, nrho, nmu, neta, nkappa, ng);
}
template void gpuPaircDensity(double *, double *, double *, double *, int*, int, int, int, int, int, int);
template void gpuPaircDensity(float *, float *, float *, float *, int *, int, int, int, int, int, int);

template <typename T> void gpuPaircDensityGradient(T *u, T *du, T *u_rho, T *rho, T *mu, T *eta, int *kappa, int nrho, int nmu, int neta, int nkappa, int ng, int potnum)
{
	if (potnum == 1)
		gpuPaircDensity1Gradient(u, du, u_rho, rho, mu, eta, kappa, nrho, nmu, neta, nkappa, ng);
	else if (potnum == 2)
		gpuPaircDensity2Gradient(u, du, u_rho, rho, mu, eta, kappa, nrho, nmu, neta, nkappa, ng);
}
template void gpuPaircDensityGradient(double *, double *, double*, double *, double *, double *, int*, int, int, int, int, int, int);
template void gpuPaircDensityGradient(float *, float *, float*, float *, float *, float *, int *, int, int, int, int, int, int);
