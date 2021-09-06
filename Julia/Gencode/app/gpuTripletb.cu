template <typename T>  __global__  void kernelgpuTripletb1(T *u, T *xij, T *xik, T *qi, T *qj, T *qk, int *ti, int *tj, int *tk, int *ai, int *aj, int *ak, T *mu, T *eta, int *kappa, int dim, int ncq, int nmu, int neta, int nkappa, int ng)
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
		u[i] = mu1*pow(-mu2 + (xij1*xik1 + xij2*xik2 + xij3*xik3)/(sqrt(pow(xij1, 2) + pow(xij2, 2) + pow(xij3, 2))*sqrt(pow(xik1, 2) + pow(xik2, 2) + pow(xik3, 2))), 2);
		i += blockDim.x * gridDim.x;
	}
}

template <typename T> void gpuTripletb1(T *u, T *xij, T *xik, T *qi, T *qj, T *qk, int *ti, int *tj, int *tk, int *ai, int *aj, int *ak, T *mu, T *eta, int *kappa, int dim, int ncq, int nmu, int neta, int nkappa, int ng)
{
	int blockDim = 256;
	int gridDim = (ng * blockDim - 1) / blockDim;
	gridDim = (gridDim>1024)? 1024 : gridDim;
	kernelgpuTripletb1<<<gridDim, blockDim>>>(u, xij, xik, qi, qj, qk, ti, tj, tk, ai, aj, ak, mu, eta, kappa, dim, ncq, nmu, neta, nkappa, ng);
}


template <typename T>  __device__  void devicegpuTripletb1(T *__restrict__ u, T *__restrict__ xij, T *__restrict__ xik, T *__restrict__ qi, T *__restrict__ qj, T *__restrict__ qk, int *__restrict__ ti, int *__restrict__ tj, int *__restrict__ tk, int *__restrict__ ai, int *__restrict__ aj, int *__restrict__ ak, T *__restrict__ mu, T *__restrict__ eta, int *__restrict__ kappa, int dim, int ncq, int nmu, int neta, int nkappa, int ng)
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
		u[i] = mu1*pow(-mu2 + (xij1*xik1 + xij2*xik2 + xij3*xik3)/(sqrt(pow(xij1, 2) + pow(xij2, 2) + pow(xij3, 2))*sqrt(pow(xik1, 2) + pow(xik2, 2) + pow(xik3, 2))), 2);
		i += blockDim.x * gridDim.x;
	}
}


template <typename T>  __global__  void kernelgpuTripletb1Gradient(T *__restrict__ u, T *__restrict__ du, T *__restrict__ u_xij, T *__restrict__ u_xik, T *__restrict__ xij, T *__restrict__ xik, T *__restrict__ qi, T *__restrict__ qj, T *__restrict__ qk, int *__restrict__ ti, int *__restrict__ tj, int *__restrict__ tk, int *__restrict__ ai, int *__restrict__ aj, int *__restrict__ ak, T *__restrict__ mu, T *__restrict__ eta, int *__restrict__ kappa, int dim, int ncq, int nmu, int neta, int nkappa, int ng)
{
	__enzyme_autodiff((void*)devicegpuTripletb1<T>, 
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

template <typename T> void gpuTripletb1Gradient(T *u, T *du, T *u_xij, T *xij, T *xik, T *qi, T *qj, T *qk, int *ti, int *tj, int *tk, int *ai, int *aj, int *ak, T *mu, T *eta, int *kappa, int dim, int ncq, int nmu, int neta, int nkappa, int ng)
{
	int blockDim = 256;
	int gridDim = (ng * blockDim - 1) / blockDim;
	gridDim = (gridDim>1024)? 1024 : gridDim;
	kernelgpuTripletb1Gradient<<<gridDim, blockDim>>>(u, du, u_xij, xij, xik, qi, qj, qk, ti, tj, tk, ai, aj, ak, mu, eta, kappa, dim, ncq, nmu, neta, nkappa, ng);
}


template <typename T>  __global__  void kernelgpuTripletb2(T *u, T *xij, T *xik, T *qi, T *qj, T *qk, int *ti, int *tj, int *tk, int *ai, int *aj, int *ak, T *mu, T *eta, int *kappa, int dim, int ncq, int nmu, int neta, int nkappa, int ng)
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
		T xik1 = xik[0 + i*3];
		T xik2 = xik[1 + i*3];
		T xik3 = xik[2 + i*3];
		T x0 = xij1*xik1;
		T x1 = xij2*xik2;
		T x2 = xij3*xik3;
		T x3 = pow(xij1, 2) + pow(xij2, 2) + pow(xij3, 2);
		T x4 = pow(xik1, 2) + pow(xik2, 2) + pow(xik3, 2);
		u[i] = mu1*pow(-mu2 + (x0 + x1 + x2)/(sqrt(x3)*sqrt(x4)), 2) + mu3*pow(mu4 - sqrt(-2*x0 - 2*x1 - 2*x2 + x3 + x4), 2);
		i += blockDim.x * gridDim.x;
	}
}

template <typename T> void gpuTripletb2(T *u, T *xij, T *xik, T *qi, T *qj, T *qk, int *ti, int *tj, int *tk, int *ai, int *aj, int *ak, T *mu, T *eta, int *kappa, int dim, int ncq, int nmu, int neta, int nkappa, int ng)
{
	int blockDim = 256;
	int gridDim = (ng * blockDim - 1) / blockDim;
	gridDim = (gridDim>1024)? 1024 : gridDim;
	kernelgpuTripletb2<<<gridDim, blockDim>>>(u, xij, xik, qi, qj, qk, ti, tj, tk, ai, aj, ak, mu, eta, kappa, dim, ncq, nmu, neta, nkappa, ng);
}


template <typename T>  __device__  void devicegpuTripletb2(T *__restrict__ u, T *__restrict__ xij, T *__restrict__ xik, T *__restrict__ qi, T *__restrict__ qj, T *__restrict__ qk, int *__restrict__ ti, int *__restrict__ tj, int *__restrict__ tk, int *__restrict__ ai, int *__restrict__ aj, int *__restrict__ ak, T *__restrict__ mu, T *__restrict__ eta, int *__restrict__ kappa, int dim, int ncq, int nmu, int neta, int nkappa, int ng)
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
		T xik1 = xik[0 + i*3];
		T xik2 = xik[1 + i*3];
		T xik3 = xik[2 + i*3];
		T x0 = xij1*xik1;
		T x1 = xij2*xik2;
		T x2 = xij3*xik3;
		T x3 = pow(xij1, 2) + pow(xij2, 2) + pow(xij3, 2);
		T x4 = pow(xik1, 2) + pow(xik2, 2) + pow(xik3, 2);
		u[i] = mu1*pow(-mu2 + (x0 + x1 + x2)/(sqrt(x3)*sqrt(x4)), 2) + mu3*pow(mu4 - sqrt(-2*x0 - 2*x1 - 2*x2 + x3 + x4), 2);
		i += blockDim.x * gridDim.x;
	}
}


template <typename T>  __global__  void kernelgpuTripletb2Gradient(T *__restrict__ u, T *__restrict__ du, T *__restrict__ u_xij, T *__restrict__ u_xik, T *__restrict__ xij, T *__restrict__ xik, T *__restrict__ qi, T *__restrict__ qj, T *__restrict__ qk, int *__restrict__ ti, int *__restrict__ tj, int *__restrict__ tk, int *__restrict__ ai, int *__restrict__ aj, int *__restrict__ ak, T *__restrict__ mu, T *__restrict__ eta, int *__restrict__ kappa, int dim, int ncq, int nmu, int neta, int nkappa, int ng)
{
	__enzyme_autodiff((void*)devicegpuTripletb2<T>, 
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

template <typename T> void gpuTripletb2Gradient(T *u, T *du, T *u_xij, T *xij, T *xik, T *qi, T *qj, T *qk, int *ti, int *tj, int *tk, int *ai, int *aj, int *ak, T *mu, T *eta, int *kappa, int dim, int ncq, int nmu, int neta, int nkappa, int ng)
{
	int blockDim = 256;
	int gridDim = (ng * blockDim - 1) / blockDim;
	gridDim = (gridDim>1024)? 1024 : gridDim;
	kernelgpuTripletb2Gradient<<<gridDim, blockDim>>>(u, du, u_xij, xij, xik, qi, qj, qk, ti, tj, tk, ai, aj, ak, mu, eta, kappa, dim, ncq, nmu, neta, nkappa, ng);
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

template <typename T> void gpuTripletbGradient(T *u, T *du, T *u_xij, T *u_xik, T *xij, T *xik, T *qi, T *qj, T *qk, int *ti, int *tj, int *tk, int *ai, int *aj, int *ak, T *mu, T *eta, int *kappa, int dim, int ncq, int nmu, int neta, int nkappa, int ng, int potnum)
{
	if (potnum == 1)
		gpuTripletb1Gradient(u, du, u_xij, u_xik, xij, xik, qi, qj, qk, ti, tj, tk, ai, aj, ak, mu, eta, kappa, dim, ncq, nmu, neta, nkappa, ng);
	else if (potnum == 2)
		gpuTripletb2Gradient(u, du, u_xij, u_xik, xij, xik, qi, qj, qk, ti, tj, tk, ai, aj, ak, mu, eta, kappa, dim, ncq, nmu, neta, nkappa, ng);
}
template void gpuTripletbGradient(double *, double *, double*, double*, double *, double *, double *, double *, double *, int *, int *, int *, int *, int *, int *, double *, double *, int*, int, int, int, int, int, int, int);
template void gpuTripletbGradient(float *, float *, float*, float*, float *, float *, float *, float *, float *, int *, int *, int *, int *, int *, int *, float *, float *, int *, int, int, int, int, int, int, int);
