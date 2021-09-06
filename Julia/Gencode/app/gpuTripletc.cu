template <typename T>  __global__  void kernelgpuTripletc1(T *u, T *xij, T *xik, T *qi, T *qj, T *qk, int *ti, int *tj, int *tk, int *ai, int *aj, int *ak, T *mu, T *eta, int *kappa, int dim, int ncq, int nmu, int neta, int nkappa, int ng)
{
	int i = threadIdx.x * blockIdx.x * blockDim.x;
	while (i<ng) {
		T mu3 = mu[2];
		T mu4 = mu[3];
		T mu5 = mu[4];
		T mu6 = mu[5];
		T mu7 = mu[6];
		T mu8 = mu[7];
		T xij1 = xij[0 + i*3];
		T xij2 = xij[1 + i*3];
		T xij3 = xij[2 + i*3];
		T xik1 = xik[0 + i*3];
		T xik2 = xik[1 + i*3];
		T xik3 = xik[2 + i*3];
		T x0 = sqrt(pow(xij1, 2) + pow(xij2, 2) + pow(xij3, 2));
		T x1 = sqrt(pow(xik1, 2) + pow(xik2, 2) + pow(xik3, 2));
		T x2 = pow(mu5, 2);
		T x3 = pow(mu6, 2);
		u[0 + i*1] = mu4*(-x2/(x3 + pow(-mu7 + (xij1*xik1 + xij2*xik2 + xij3*xik3)/(x0*x1), 2)) + x2/x3 + 1)*exp(pow(mu3, mu8)*pow(x0 - x1, mu8));
		i += blockDim.x * gridDim.x;
	}
}

template <typename T> void gpuTripletc1(T *u, T *xij, T *xik, T *qi, T *qj, T *qk, int *ti, int *tj, int *tk, int *ai, int *aj, int *ak, T *mu, T *eta, int *kappa, int dim, int ncq, int nmu, int neta, int nkappa, int ng)
{
	int blockDim = 256;
	int gridDim = (ng * blockDim - 1) / blockDim;
	gridDim = (gridDim>1024)? 1024 : gridDim;
	kernelgpuTripletc1<<<gridDim, blockDim>>>(u, xij, xik, qi, qj, qk, ti, tj, tk, ai, aj, ak, mu, eta, kappa, dim, ncq, nmu, neta, nkappa, ng);
}


template <typename T>  __device__  void devicegpuTripletc1(T *__restrict__ u, T *__restrict__ xij, T *__restrict__ xik, T *__restrict__ qi, T *__restrict__ qj, T *__restrict__ qk, int *__restrict__ ti, int *__restrict__ tj, int *__restrict__ tk, int *__restrict__ ai, int *__restrict__ aj, int *__restrict__ ak, T *__restrict__ mu, T *__restrict__ eta, int *__restrict__ kappa, int dim, int ncq, int nmu, int neta, int nkappa, int ng)
{
	int i = threadIdx.x * blockIdx.x * blockDim.x;
	while (i<ng) {
		T mu3 = mu[2];
		T mu4 = mu[3];
		T mu5 = mu[4];
		T mu6 = mu[5];
		T mu7 = mu[6];
		T mu8 = mu[7];
		T xij1 = xij[0 + i*3];
		T xij2 = xij[1 + i*3];
		T xij3 = xij[2 + i*3];
		T xik1 = xik[0 + i*3];
		T xik2 = xik[1 + i*3];
		T xik3 = xik[2 + i*3];
		T x0 = sqrt(pow(xij1, 2) + pow(xij2, 2) + pow(xij3, 2));
		T x1 = sqrt(pow(xik1, 2) + pow(xik2, 2) + pow(xik3, 2));
		T x2 = pow(mu5, 2);
		T x3 = pow(mu6, 2);
		u[0 + i*1] = mu4*(-x2/(x3 + pow(-mu7 + (xij1*xik1 + xij2*xik2 + xij3*xik3)/(x0*x1), 2)) + x2/x3 + 1)*exp(pow(mu3, mu8)*pow(x0 - x1, mu8));
		i += blockDim.x * gridDim.x;
	}
}


template <typename T>  __global__  void kernelgpuTripletc1Gradient(T *__restrict__ u, T *__restrict__ du, T *__restrict__ u_xij, T *__restrict__ u_xik, T *__restrict__ xij, T *__restrict__ xik, T *__restrict__ qi, T *__restrict__ qj, T *__restrict__ qk, int *__restrict__ ti, int *__restrict__ tj, int *__restrict__ tk, int *__restrict__ ai, int *__restrict__ aj, int *__restrict__ ak, T *__restrict__ mu, T *__restrict__ eta, int *__restrict__ kappa, int dim, int ncq, int nmu, int neta, int nkappa, int ng)
{
	__enzyme_autodiff((void*)devicegpuTripletc1<T>, 
		enzyme_dup, u, du, 
		enzyme_dup, xij, u_xij, 
		enzyme_dup, xik, u_xik, 
		enzyme_const, qi, 
		enzyme_const, qj, 
		enzyme_const, qk, 
		enzyme_const, ti, 
		enzyme_const, tj, 
		enzyme_const, tk, 
		enzyme_const, ai, 
		enzyme_const, aj, 
		enzyme_const, ak, 
		enzyme_const, mu, 
		enzyme_const, eta, 
		enzyme_const, kappa, 
		dim, ncq, nmu, neta, nkappa, ng); 
}

template <typename T> void gpuTripletc1Gradient(T *u, T *du, T *u_xij, T *xij, T *xik, T *qi, T *qj, T *qk, int *ti, int *tj, int *tk, int *ai, int *aj, int *ak, T *mu, T *eta, int *kappa, int dim, int ncq, int nmu, int neta, int nkappa, int ng)
{
	int blockDim = 256;
	int gridDim = (ng * blockDim - 1) / blockDim;
	gridDim = (gridDim>1024)? 1024 : gridDim;
	kernelgpuTripletc1Gradient<<<gridDim, blockDim>>>(u, du, u_xij, xij, xik, qi, qj, qk, ti, tj, tk, ai, aj, ak, mu, eta, kappa, dim, ncq, nmu, neta, nkappa, ng);
}


template <typename T>  __global__  void kernelgpuTripletcPair1(T *u, T *xij, T *qi, T *qj, int *ti, int *tj, int *ai, int *aj, T *mu, T *eta, int *kappa, int dim, int ncq, int nmu, int neta, int nkappa, int ng)
{
	int i = threadIdx.x * blockIdx.x * blockDim.x;
	while (i<ng) {
		T mu1 = mu[0];
		T mu2 = mu[1];
		T xij1 = xij[0 + i*3];
		T xij2 = xij[1 + i*3];
		T xij3 = xij[2 + i*3];
		u[0 + i*1] = -mu1*exp(-mu2*sqrt(pow(xij1, 2) + pow(xij2, 2) + pow(xij3, 2)));
		i += blockDim.x * gridDim.x;
	}
}

template <typename T> void gpuTripletcPair1(T *u, T *xij, T *qi, T *qj, int *ti, int *tj, int *ai, int *aj, T *mu, T *eta, int *kappa, int dim, int ncq, int nmu, int neta, int nkappa, int ng)
{
	int blockDim = 256;
	int gridDim = (ng * blockDim - 1) / blockDim;
	gridDim = (gridDim>1024)? 1024 : gridDim;
	kernelgpuTripletcPair1<<<gridDim, blockDim>>>(u, xij, qi, qj, ti, tj, ai, aj, mu, eta, kappa, dim, ncq, nmu, neta, nkappa, ng);
}

template <typename T>  __device__  void devicegpuTripletcPair1(T *__restrict__ u, T *__restrict__ xij, T *__restrict__ qi, T *__restrict__ qj, int *__restrict__ ti, int *__restrict__ tj, int *__restrict__ ai, int *__restrict__ aj, T *__restrict__ mu, T *__restrict__ eta, int *__restrict__ kappa, int dim, int ncq, int nmu, int neta, int nkappa, int ng)
{
	int i = threadIdx.x * blockIdx.x * blockDim.x;
	while (i<ng) {
		T mu1 = mu[0];
		T mu2 = mu[1];
		T xij1 = xij[0 + i*3];
		T xij2 = xij[1 + i*3];
		T xij3 = xij[2 + i*3];
		u[0 + i*1] = -mu1*exp(-mu2*sqrt(pow(xij1, 2) + pow(xij2, 2) + pow(xij3, 2)));
		i += blockDim.x * gridDim.x;
	}
}


template <typename T>  __global__  void kernelgpuTripletcPair1Gradient(T *__restrict__ u, T *__restrict__ du, T *__restrict__ u_xij, T *__restrict__ xij, T *__restrict__ qi, T *__restrict__ qj, int *__restrict__ ti, int *__restrict__ tj, int *__restrict__ ai, int *__restrict__ aj, T *__restrict__ mu, T *__restrict__ eta, int *__restrict__ kappa, int dim, int ncq, int nmu, int neta, int nkappa, int ng)
{
	__enzyme_autodiff((void*)devicegpuTripletcPair1<T>, 
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

template <typename T> void gpuTripletcPair1Gradient(T *u, T *du, T *u_xij, T *xij, T *qi, T *qj, int *ti, int *tj, int *ai, int *aj, T *mu, T *eta, int *kappa, int dim, int ncq, int nmu, int neta, int nkappa, int ng)
{
	int blockDim = 256;
	int gridDim = (ng * blockDim - 1) / blockDim;
	gridDim = (gridDim>1024)? 1024 : gridDim;
	kernelgpuTripletcPair1Gradient<<<gridDim, blockDim>>>(u, du, u_xij, xij, qi, qj, ti, tj, ai, aj, mu, eta, kappa, dim, ncq, nmu, neta, nkappa, ng);
}

template <typename T>  __global__  void kernelgpuTripletcDensity1(T *u, T *rho, T *mu, T *eta, int *kappa, int nrho, int nmu, int neta, int nkappa, int ng)
{
	int i = threadIdx.x * blockIdx.x * blockDim.x;
	while (i<ng) {
		T mu1 = mu[0];
		T mu9 = mu[8];
		T mu10 = mu[9];
		T rho1 = rho[0 + i*1];
		u[0 + i*1] = pow(pow(mu10, mu9)*pow(rho1, mu9) + 1, -0.5/mu9);
		i += blockDim.x * gridDim.x;
	}
}

template <typename T> void gpuTripletcDensity1(T *u, T *rho, T *mu, T *eta, int *kappa, int nrho, int nmu, int neta, int nkappa, int ng)
{
	int blockDim = 256;
	int gridDim = (ng * blockDim - 1) / blockDim;
	gridDim = (gridDim>1024)? 1024 : gridDim;
	kernelgpuTripletcDensity1<<<gridDim, blockDim>>>(u, rho, mu, eta, kappa, nrho, nmu, neta, nkappa, ng);
}


template <typename T>  __device__  void devicegpuTripletcDensity1(T *__restrict__ u, T *__restrict__ rho, T *__restrict__ mu, T *__restrict__ eta, int *__restrict__ kappa, int nrho, int nmu, int neta, int nkappa, int ng)
{
	int i = threadIdx.x * blockIdx.x * blockDim.x;
	while (i<ng) {
		T mu1 = mu[0];
		T mu9 = mu[8];
		T mu10 = mu[9];
		T rho1 = rho[0 + i*1];
		u[0 + i*1] = pow(pow(mu10, mu9)*pow(rho1, mu9) + 1, -0.5/mu9);
		i += blockDim.x * gridDim.x;
	}
}


template <typename T>  __global__  void kernelgpuTripletcDensity1Gradient(T *__restrict__ u, T *__restrict__ du, T *__restrict__ u_rho, T *__restrict__ rho, T *__restrict__ mu, T *__restrict__ eta, int *__restrict__ kappa, int nrho, int nmu, int neta, int nkappa, int ng)
{
	__enzyme_autodiff((void*)devicegpuTripletcDensity1<T>, 
		enzyme_dup, u, du, 
		enzyme_dup, rho, u_rho, 
		enzyme_const, mu, 
		enzyme_const, eta, 
		enzyme_const, kappa, 
		dim, nrho, nmu, neta, nkappa, ng); 
}

template <typename T> void gpuTripletcDensity1Gradient(T *u, T *du, T *u_rho, T *rho, T *mu, T *eta, int *kappa, int nrho, int nmu, int neta, int nkappa, int ng)
{
	int blockDim = 256;
	int gridDim = (ng * blockDim - 1) / blockDim;
	gridDim = (gridDim>1024)? 1024 : gridDim;
	kernelgpuTripletcDensity1Gradient<<<gridDim, blockDim>>>(u, du, u_rho, rho, mu, eta, kappa, nrho, nmu, neta, nkappa, ng);
}


template <typename T> void gpuTripletc(T *u, T *xij, T *xik, T *qi, T *qj, T *qk, int *ti, int *tj, int *tk, int *ai, int *aj, int *ak, T *mu, T *eta, int *kappa, int dim, int ncq, int nmu, int neta, int nkappa, int ng, int potnum)
{
	if (potnum == 1)
		gpuTripletc1(u, xij, xik, qi, qj, qk, ti, tj, tk, ai, aj, ak, mu, eta, kappa, dim, ncq, nmu, neta, nkappa, ng);
}
template void gpuTripletc(double *, double *, double *, double *, double *, double *, int *, int *, int *, int *, int *, int *, double *, double *, int*, int, int, int, int, int, int, int);
template void gpuTripletc(float *, float *, float *, float *, float *, float *, int *, int *, int *, int *, int *, int *, float *, float *, int *, int, int, int, int, int, int, int);


template <typename T> void gpuTripletcGradient(T *u, T *du, T *u_xij, T *u_xik, T *xij, T *xik, T *qi, T *qj, T *qk, int *ti, int *tj, int *tk, int *ai, int *aj, int *ak, T *mu, T *eta, int *kappa, int dim, int ncq, int nmu, int neta, int nkappa, int ng, int potnum)
{
	if (potnum == 1)
		gpuTripletc1Gradient(u, du, u_xij, u_xik, xij, xik, qi, qj, qk, ti, tj, tk, ai, aj, ak, mu, eta, kappa, dim, ncq, nmu, neta, nkappa, ng);
}
template void gpuTripletcGradient(double *, double *, double*, double*, double *, double *, double *, double *, double *, int *, int *, int *, int *, int *, int *, double *, double *, int*, int, int, int, int, int, int, int);
template void gpuTripletcGradient(float *, float *, float*, float*, float *, float *, float *, float *, float *, int *, int *, int *, int *, int *, int *, float *, float *, int *, int, int, int, int, int, int, int);
template <typename T> void gpuTripletcPair(T *u, T *xij, T *qi, T *qj, int *ti, int *tj, int *ai, int *aj, T *mu, T *eta, int *kappa, int dim, int ncq, int nmu, int neta, int nkappa, int ng, int potnum)
{
	if (potnum == 1)
		gpuTripletcPair1(u, xij, qi, qj, ti, tj, ai, aj, mu, eta, kappa, dim, ncq, nmu, neta, nkappa, ng);
}
template void gpuTripletcPair(double *, double *, double *, double *, int *, int *, int *, int *, double *, double *, int*, int, int, int, int, int, int, int);
template void gpuTripletcPair(float *, float *, float *, float *, int *, int *, int *, int *, float *, float *, int *, int, int, int, int, int, int, int);


template <typename T> void gpuTripletcPairGradient(T *u, T *du, T *u_xij, T *xij, T *qi, T *qj, int *ti, int *tj, int *ai, int *aj, T *mu, T *eta, int *kappa, int dim, int ncq, int nmu, int neta, int nkappa, int ng, int potnum)
{
	if (potnum == 1)
		gpuTripletcPair1Gradient(u, du, u_xij, xij, qi, qj, ti, tj, ai, aj, mu, eta, kappa, dim, ncq, nmu, neta, nkappa, ng);
}
template void gpuTripletcPairGradient(double *, double *, double*, double *, double *, double *, int *, int *, int *, int *, double *, double *, int*, int, int, int, int, int, int, int);
template void gpuTripletcPairGradient(float *, float *, float*, float *, float *, float *, int *, int *, int *, int *, float *, float *, int *, int, int, int, int, int, int, int);
template <typename T> void gpuTripletcDensity(T *u, T *rho, T *mu, T *eta, int *kappa, int nrho, int nmu, int neta, int nkappa, int ng, int potnum)
{
	if (potnum == 1)
		gpuTripletcDensity1(u, rho, mu, eta, kappa, nrho, nmu, neta, nkappa, ng);
}
template void gpuTripletcDensity(double *, double *, double *, double *, int*, int, int, int, int, int, int);
template void gpuTripletcDensity(float *, float *, float *, float *, int *, int, int, int, int, int, int);

template <typename T> void gpuTripletcDensityGradient(T *u, T *du, T *u_rho, T *rho, T *mu, T *eta, int *kappa, int nrho, int nmu, int neta, int nkappa, int ng, int potnum)
{
	if (potnum == 1)
		gpuTripletcDensity1Gradient(u, du, u_rho, rho, mu, eta, kappa, nrho, nmu, neta, nkappa, ng);
}
template void gpuTripletcDensityGradient(double *, double *, double*, double *, double *, double *, int*, int, int, int, int, int, int);
template void gpuTripletcDensityGradient(float *, float *, float*, float *, float *, float *, int *, int, int, int, int, int, int);
