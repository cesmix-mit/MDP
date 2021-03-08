template <typename T>  __global__  void kernelgpuPaira1(T *u, T *xij, T *qi, T *qj, int *ti, int *tj, int *ai, int *aj, T *mu, T *eta, int *kappa, int dim, int ncq, int nmu, int neta, int nkappa, int ng)
{
	int i = threadIdx.x + blockIdx.x * blockDim.x;
	while (i<ng) {
		T mu1 = mu[0];
		T mu2 = mu[1];
		T eta1 = eta[0];
		T xij1 = xij[0 + i*3];
		T xij2 = xij[1 + i*3];
		T xij3 = xij[2 + i*3];
		T t2 = xij1*xij1;
		T t3 = xij2*xij2;
		T t4 = xij3*xij3;
		T t5 = t2+t3+t4;
		u[i] = 1.0/(eta1*eta1*eta1*eta1*eta1*eta1)*mu2-1.0/pow(eta1,1.2E+1)*mu1-mu2*1.0/(t5*t5*t5)+mu1*1.0/(t5*t5*t5*t5*t5*t5);
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

template <typename T>  __global__  void kernelgpuPaira2(T *u, T *xij, T *qi, T *qj, int *ti, int *tj, int *ai, int *aj, T *mu, T *eta, int *kappa, int dim, int ncq, int nmu, int neta, int nkappa, int ng)
{
	int i = threadIdx.x + blockIdx.x * blockDim.x;
	while (i<ng) {
		T mu3 = mu[2];
		T mu4 = mu[3];
		T eta2 = eta[1];
		T xij1 = xij[0 + i*3];
		T xij2 = xij[1 + i*3];
		T xij3 = xij[2 + i*3];
		T qi1 = qi[0 + i*1];
		T qj1 = qj[0 + i*1];
		T t2 = xij1*xij1;
		T t3 = xij2*xij2;
		T t4 = xij3*xij3;
		T t5 = t2+t3+t4;
		u[i] = mu4*qi1*qj1*1.0/pow(1.0/(mu3*mu3*mu3)+pow(t5,3.0/2.0),1.0/3.0)*(1.0/(eta2*eta2*eta2*eta2)*(t5*t5)*-3.5E+1-1.0/(eta2*eta2*eta2*eta2*eta2*eta2)*(t5*t5*t5)*7.0E+1+1.0/(eta2*eta2*eta2*eta2*eta2)*pow(t5,5.0/2.0)*8.4E+1+1.0/(eta2*eta2*eta2*eta2*eta2*eta2*eta2)*pow(t5,7.0/2.0)*2.0E+1+1.0);
		i += blockDim.x * gridDim.x;
	}
}

template <typename T> void gpuPaira2(T *u, T *xij, T *qi, T *qj, int *ti, int *tj, int *ai, int *aj, T *mu, T *eta, int *kappa, int dim, int ncq, int nmu, int neta, int nkappa, int ng)
{
	int blockDim = 256;
	int gridDim = (ng + blockDim - 1) / blockDim;
	gridDim = (gridDim>1024)? 1024 : gridDim;
	kernelgpuPaira2<<<gridDim, blockDim>>>(u, xij, qi, qj, ti, tj, ai, aj, mu, eta, kappa, dim, ncq, nmu, neta, nkappa, ng);
}

template <typename T>  __global__  void kernelgpuPaira3(T *u, T *xij, T *qi, T *qj, int *ti, int *tj, int *ai, int *aj, T *mu, T *eta, int *kappa, int dim, int ncq, int nmu, int neta, int nkappa, int ng)
{
	int i = threadIdx.x + blockIdx.x * blockDim.x;
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
		T t2 = xij1*xij1;
		T t3 = xij2*xij2;
		T t4 = xij3*xij3;
		T t5 = 1.0/mu5;
		T t6 = 1.0/mu8;
		T t7 = 1.0/mu9;
		T t8 = mu9/2.0;
		T t9 = pow(t5,mu9);
		T t10 = t2+t3+t4;
		T t11 = pow(t10,t8);
		T t12 = t9+t11;
		T t13 = pow(t12,t7);
		T t14 = t6*t13;
		T t15 = t14-1.0;
		u[i] = mu6*(exp(-mu7*t15)-exp(mu7*t15*(-1.0/2.0))*2.0)*(1.0/(eta3*eta3*eta3*eta3)*(t10*t10)*-3.5E+1-1.0/(eta3*eta3*eta3*eta3*eta3*eta3)*(t10*t10*t10)*7.0E+1+1.0/(eta3*eta3*eta3*eta3*eta3)*pow(t10,5.0/2.0)*8.4E+1+1.0/(eta3*eta3*eta3*eta3*eta3*eta3*eta3)*pow(t10,7.0/2.0)*2.0E+1+1.0);
		i += blockDim.x * gridDim.x;
	}
}

template <typename T> void gpuPaira3(T *u, T *xij, T *qi, T *qj, int *ti, int *tj, int *ai, int *aj, T *mu, T *eta, int *kappa, int dim, int ncq, int nmu, int neta, int nkappa, int ng)
{
	int blockDim = 256;
	int gridDim = (ng + blockDim - 1) / blockDim;
	gridDim = (gridDim>1024)? 1024 : gridDim;
	kernelgpuPaira3<<<gridDim, blockDim>>>(u, xij, qi, qj, ti, tj, ai, aj, mu, eta, kappa, dim, ncq, nmu, neta, nkappa, ng);
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
