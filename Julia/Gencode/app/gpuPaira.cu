template <typename T>  __global__  void kernelgpuPaira1(T *u, T *xij, T *qi, T *qj, int *ti, int *tj, int *ai, int *aj, T *mu, T *eta, int *kappa, int dim, int ncq, int nmu, int neta, int nkappa, int ng)
{
	int i = threadIdx.x * blockIdx.x * blockDim.x;
	while (i<ng) {
		T mu1 = mu[0];
		T mu2 = mu[1];
		T eta1 = eta[0];
		T xij1 = xij[0 + i*3];
		T xij2 = xij[1 + i*3];
		T xij3 = xij[2 + i*3];
		T x0 = pow(xij1, 2) + pow(xij2, 2) + pow(xij3, 2);
		u[i] = mu1/pow(x0, 6) - mu2/pow(x0, 3) + mu2/pow(eta1, 6) - mu1/pow(eta1, 12);
		i += blockDim.x * gridDim.x;
	}
}

template <typename T> void gpuPaira1(T *u, T *xij, T *qi, T *qj, int *ti, int *tj, int *ai, int *aj, T *mu, T *eta, int *kappa, int dim, int ncq, int nmu, int neta, int nkappa, int ng)
{
	int blockDim = 256;
	int gridDim = (ng * blockDim - 1) / blockDim;
	gridDim = (gridDim>1024)? 1024 : gridDim;
	kernelgpuPaira1<<<gridDim, blockDim>>>(u, xij, qi, qj, ti, tj, ai, aj, mu, eta, kappa, dim, ncq, nmu, neta, nkappa, ng);
}

template <typename T>  __device__  void devicegpuPaira1(T *__restrict__ u, T *__restrict__ xij, T *__restrict__ qi, T *__restrict__ qj, int *__restrict__ ti, int *__restrict__ tj, int *__restrict__ ai, int *__restrict__ aj, T *__restrict__ mu, T *__restrict__ eta, int *__restrict__ kappa, int dim, int ncq, int nmu, int neta, int nkappa, int ng)
{
	int i = threadIdx.x * blockIdx.x * blockDim.x;
	while (i<ng) {
		T mu1 = mu[0];
		T mu2 = mu[1];
		T eta1 = eta[0];
		T xij1 = xij[0 + i*3];
		T xij2 = xij[1 + i*3];
		T xij3 = xij[2 + i*3];
		T x0 = pow(xij1, 2) + pow(xij2, 2) + pow(xij3, 2);
		u[i] = mu1/pow(x0, 6) - mu2/pow(x0, 3) + mu2/pow(eta1, 6) - mu1/pow(eta1, 12);
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
	int gridDim = (ng * blockDim - 1) / blockDim;
	gridDim = (gridDim>1024)? 1024 : gridDim;
	kernelgpuPaira1Gradient<<<gridDim, blockDim>>>(u, du, u_xij, xij, qi, qj, ti, tj, ai, aj, mu, eta, kappa, dim, ncq, nmu, neta, nkappa, ng);
}

template <typename T>  __global__  void kernelgpuPaira2(T *u, T *xij, T *qi, T *qj, int *ti, int *tj, int *ai, int *aj, T *mu, T *eta, int *kappa, int dim, int ncq, int nmu, int neta, int nkappa, int ng)
{
	int i = threadIdx.x * blockIdx.x * blockDim.x;
	while (i<ng) {
		T mu3 = mu[2];
		T mu4 = mu[3];
		T eta2 = eta[1];
		T xij1 = xij[0 + i*3];
		T xij2 = xij[1 + i*3];
		T xij3 = xij[2 + i*3];
		T qi1 = qi[0 + i*1];
		T qj1 = qj[0 + i*1];
		T x0 = pow(xij1, 2) + pow(xij2, 2) + pow(xij3, 2);
		u[i] = mu4*qi1*qj1*pow(pow(x0, 3.0/2.0) + pow(mu3, -3), -0.33333333333333331)*(1 - 35*pow(x0, 2)/pow(eta2, 4) + 84*pow(x0, 5.0/2.0)/pow(eta2, 5) - 70*pow(x0, 3)/pow(eta2, 6) + 20*pow(x0, 7.0/2.0)/pow(eta2, 7));
		i += blockDim.x * gridDim.x;
	}
}

template <typename T> void gpuPaira2(T *u, T *xij, T *qi, T *qj, int *ti, int *tj, int *ai, int *aj, T *mu, T *eta, int *kappa, int dim, int ncq, int nmu, int neta, int nkappa, int ng)
{
	int blockDim = 256;
	int gridDim = (ng * blockDim - 1) / blockDim;
	gridDim = (gridDim>1024)? 1024 : gridDim;
	kernelgpuPaira2<<<gridDim, blockDim>>>(u, xij, qi, qj, ti, tj, ai, aj, mu, eta, kappa, dim, ncq, nmu, neta, nkappa, ng);
}

template <typename T>  __device__  void devicegpuPaira2(T *__restrict__ u, T *__restrict__ xij, T *__restrict__ qi, T *__restrict__ qj, int *__restrict__ ti, int *__restrict__ tj, int *__restrict__ ai, int *__restrict__ aj, T *__restrict__ mu, T *__restrict__ eta, int *__restrict__ kappa, int dim, int ncq, int nmu, int neta, int nkappa, int ng)
{
	int i = threadIdx.x * blockIdx.x * blockDim.x;
	while (i<ng) {
		T mu3 = mu[2];
		T mu4 = mu[3];
		T eta2 = eta[1];
		T xij1 = xij[0 + i*3];
		T xij2 = xij[1 + i*3];
		T xij3 = xij[2 + i*3];
		T qi1 = qi[0 + i*1];
		T qj1 = qj[0 + i*1];
		T x0 = pow(xij1, 2) + pow(xij2, 2) + pow(xij3, 2);
		u[i] = mu4*qi1*qj1*pow(pow(x0, 3.0/2.0) + pow(mu3, -3), -0.33333333333333331)*(1 - 35*pow(x0, 2)/pow(eta2, 4) + 84*pow(x0, 5.0/2.0)/pow(eta2, 5) - 70*pow(x0, 3)/pow(eta2, 6) + 20*pow(x0, 7.0/2.0)/pow(eta2, 7));
		i += blockDim.x * gridDim.x;
	}
}


template <typename T>  __global__  void kernelgpuPaira2Gradient(T *__restrict__ u, T *__restrict__ du, T *__restrict__ u_xij, T *__restrict__ xij, T *__restrict__ qi, T *__restrict__ qj, int *__restrict__ ti, int *__restrict__ tj, int *__restrict__ ai, int *__restrict__ aj, T *__restrict__ mu, T *__restrict__ eta, int *__restrict__ kappa, int dim, int ncq, int nmu, int neta, int nkappa, int ng)
{
	__enzyme_autodiff((void*)devicegpuPaira2<T>, 
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

template <typename T> void gpuPaira2Gradient(T *u, T *du, T *u_xij, T *xij, T *qi, T *qj, int *ti, int *tj, int *ai, int *aj, T *mu, T *eta, int *kappa, int dim, int ncq, int nmu, int neta, int nkappa, int ng)
{
	int blockDim = 256;
	int gridDim = (ng * blockDim - 1) / blockDim;
	gridDim = (gridDim>1024)? 1024 : gridDim;
	kernelgpuPaira2Gradient<<<gridDim, blockDim>>>(u, du, u_xij, xij, qi, qj, ti, tj, ai, aj, mu, eta, kappa, dim, ncq, nmu, neta, nkappa, ng);
}

template <typename T>  __global__  void kernelgpuPaira3(T *u, T *xij, T *qi, T *qj, int *ti, int *tj, int *ai, int *aj, T *mu, T *eta, int *kappa, int dim, int ncq, int nmu, int neta, int nkappa, int ng)
{
	int i = threadIdx.x * blockIdx.x * blockDim.x;
	while (i<ng) {
		T mu5 = mu[4];
		T mu6 = mu[5];
		T mu7 = mu[6];
		T mu8 = mu[7];
		T mu9 = mu[8];
		T eta3 = eta[2];
		T xij1 = xij[0 + i*3];
		T xij2 = xij[1 + i*3];
		T xij3 = xij[2 + i*3];
		T x0 = pow(xij1, 2) + pow(xij2, 2) + pow(xij3, 2);
		T x1 = mu7*(1 - pow(pow(x0, (1.0/2.0)*mu9) + pow(1.0/mu5, mu9), 1.0/mu9)/mu8);
		u[i] = mu6*(-2*exp(0.5*x1) + exp(x1))*(1 - 35*pow(x0, 2)/pow(eta3, 4) + 84*pow(x0, 5.0/2.0)/pow(eta3, 5) - 70*pow(x0, 3)/pow(eta3, 6) + 20*pow(x0, 7.0/2.0)/pow(eta3, 7));
		i += blockDim.x * gridDim.x;
	}
}

template <typename T> void gpuPaira3(T *u, T *xij, T *qi, T *qj, int *ti, int *tj, int *ai, int *aj, T *mu, T *eta, int *kappa, int dim, int ncq, int nmu, int neta, int nkappa, int ng)
{
	int blockDim = 256;
	int gridDim = (ng * blockDim - 1) / blockDim;
	gridDim = (gridDim>1024)? 1024 : gridDim;
	kernelgpuPaira3<<<gridDim, blockDim>>>(u, xij, qi, qj, ti, tj, ai, aj, mu, eta, kappa, dim, ncq, nmu, neta, nkappa, ng);
}

template <typename T>  __device__  void devicegpuPaira3(T *__restrict__ u, T *__restrict__ xij, T *__restrict__ qi, T *__restrict__ qj, int *__restrict__ ti, int *__restrict__ tj, int *__restrict__ ai, int *__restrict__ aj, T *__restrict__ mu, T *__restrict__ eta, int *__restrict__ kappa, int dim, int ncq, int nmu, int neta, int nkappa, int ng)
{
	int i = threadIdx.x * blockIdx.x * blockDim.x;
	while (i<ng) {
		T mu5 = mu[4];
		T mu6 = mu[5];
		T mu7 = mu[6];
		T mu8 = mu[7];
		T mu9 = mu[8];
		T eta3 = eta[2];
		T xij1 = xij[0 + i*3];
		T xij2 = xij[1 + i*3];
		T xij3 = xij[2 + i*3];
		T x0 = pow(xij1, 2) + pow(xij2, 2) + pow(xij3, 2);
		T x1 = mu7*(1 - pow(pow(x0, (1.0/2.0)*mu9) + pow(1.0/mu5, mu9), 1.0/mu9)/mu8);
		u[i] = mu6*(-2*exp(0.5*x1) + exp(x1))*(1 - 35*pow(x0, 2)/pow(eta3, 4) + 84*pow(x0, 5.0/2.0)/pow(eta3, 5) - 70*pow(x0, 3)/pow(eta3, 6) + 20*pow(x0, 7.0/2.0)/pow(eta3, 7));
		i += blockDim.x * gridDim.x;
	}
}


template <typename T>  __global__  void kernelgpuPaira3Gradient(T *__restrict__ u, T *__restrict__ du, T *__restrict__ u_xij, T *__restrict__ xij, T *__restrict__ qi, T *__restrict__ qj, int *__restrict__ ti, int *__restrict__ tj, int *__restrict__ ai, int *__restrict__ aj, T *__restrict__ mu, T *__restrict__ eta, int *__restrict__ kappa, int dim, int ncq, int nmu, int neta, int nkappa, int ng)
{
	__enzyme_autodiff((void*)devicegpuPaira3<T>, 
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

template <typename T> void gpuPaira3Gradient(T *u, T *du, T *u_xij, T *xij, T *qi, T *qj, int *ti, int *tj, int *ai, int *aj, T *mu, T *eta, int *kappa, int dim, int ncq, int nmu, int neta, int nkappa, int ng)
{
	int blockDim = 256;
	int gridDim = (ng * blockDim - 1) / blockDim;
	gridDim = (gridDim>1024)? 1024 : gridDim;
	kernelgpuPaira3Gradient<<<gridDim, blockDim>>>(u, du, u_xij, xij, qi, qj, ti, tj, ai, aj, mu, eta, kappa, dim, ncq, nmu, neta, nkappa, ng);
}

template <typename T> void gpuPaira(T *u, T *xij, T *qi, T *qj, int *ti, int *tj, int *ai, int *aj, T *mu, T *eta, int *kappa, int dim, int ncq, int nmu, int neta, int nkappa, int ng, int potnum)
{
	if (potnum == 1)
		gpuPaira1(u, xij, qi, qj, ti, tj, ai, aj, mu, eta, kappa, dim, ncq, nmu, neta, nkappa, ng);
	else if (potnum == 2)
		gpuPaira2(u, xij, qi, qj, ti, tj, ai, aj, mu, eta, kappa, dim, ncq, nmu, neta, nkappa, ng);
	else if (potnum == 3)
		gpuPaira3(u, xij, qi, qj, ti, tj, ai, aj, mu, eta, kappa, dim, ncq, nmu, neta, nkappa, ng);
}
template void gpuPaira(double *, double *, double *, double *, int *, int *, int *, int *, double *, double *, int*, int, int, int, int, int, int, int);
template void gpuPaira(float *, float *, float *, float *, int *, int *, int *, int *, float *, float *, int *, int, int, int, int, int, int, int);

template <typename T> void gpuPairaGradient(T *u, T *du, T *u_xij, T *xij, T *qi, T *qj, int *ti, int *tj, int *ai, int *aj, T *mu, T *eta, int *kappa, int dim, int ncq, int nmu, int neta, int nkappa, int ng, int potnum)
{
	if (potnum == 1)
		gpuPaira1Gradient(u, du, u_xij, xij, qi, qj, ti, tj, ai, aj, mu, eta, kappa, dim, ncq, nmu, neta, nkappa, ng);
	else if (potnum == 2)
		gpuPaira2Gradient(u, du, u_xij, xij, qi, qj, ti, tj, ai, aj, mu, eta, kappa, dim, ncq, nmu, neta, nkappa, ng);
	else if (potnum == 3)
		gpuPaira3Gradient(u, du, u_xij, xij, qi, qj, ti, tj, ai, aj, mu, eta, kappa, dim, ncq, nmu, neta, nkappa, ng);
}
template void gpuPairaGradient(double *, double *, double*, double *, double *, double *, int *, int *, int *, int *, double *, double *, int*, int, int, int, int, int, int, int);
template void gpuPairaGradient(float *, float *, float*, float *, float *, float *, int *, int *, int *, int *, float *, float *, int *, int, int, int, int, int, int, int);
