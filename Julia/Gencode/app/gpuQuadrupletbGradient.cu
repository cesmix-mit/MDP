template <typename T>  __global__  void kernelgpuQuadrupletbGradient1(T *u, T *u_xij, T *u_xik, T *u_xil, T *xij, T *xik, T *xil, T *qi, T *qj, T *qk, T *ql, int *ti, int *tj, int *tk, int *tl, int *ai, int *aj, int *ak, int *al, T *mu, T *eta, int *kappa, int dim, int ncq, int nmu, int neta, int nkappa, int ng)
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
		T x4 = x1*x3;
		T x5 = x0 + xil2;
		T x6 = x2 + xik3;
		T x7 = x5*x6;
		T x8 = x4 - x7;
		T x9 = -xij1;
		T x10 = x9 + xil1;
		T x11 = x10*x6;
		T x12 = x9 + xik1;
		T x13 = x12*x3;
		T x14 = x11 - x13;
		T x15 = x12*x5;
		T x16 = x1*x10;
		T x17 = x15 - x16;
		T x18 = x14*xij2 + x17*xij3 + x8*xij1;
		T x19 = pow(x18, 2);
		T x20 = pow(x14, 2) + pow(x17, 2) + pow(x8, 2);
		T x21 = mu1/x20;
		T x22 = pow(x20, -2);
		T x23 = mu2*pow(x18, 4);
		T x24 = -xik3;
		T x25 = x24 + xil3;
		T x26 = 2*xij2;
		T x27 = -xil2;
		T x28 = x27 + xik2;
		T x29 = 2*xij3;
		T x30 = x18*x21;
		T x31 = 4*xij2;
		T x32 = 4*xij3;
		T x33 = mu2*pow(x18, 3)*x22;
		T x34 = 2*xik2;
		T x35 = 2*xil2;
		T x36 = -x35;
		T x37 = x17*(x34 + x36);
		T x38 = 2*xik3;
		T x39 = -x38;
		T x40 = 2*xil3;
		T x41 = x14*(x39 + x40);
		T x42 = mu1*x19*x22;
		T x43 = x23/pow(x20, 3);
		T x44 = -xil3;
		T x45 = x44 + xik3;
		T x46 = 2*xij1;
		T x47 = -xik1;
		T x48 = x47 + xil1;
		T x49 = 4*xij1;
		T x50 = 2*xik1;
		T x51 = -x50;
		T x52 = 2*xil1;
		T x53 = x17*(x51 + x52);
		T x54 = -x40;
		T x55 = x8*(x38 + x54);
		T x56 = -xik2;
		T x57 = x56 + xil2;
		T x58 = -xil1;
		T x59 = x58 + xik1;
		T x60 = -x52;
		T x61 = x14*(x50 + x60);
		T x62 = -x34;
		T x63 = x8*(x35 + x62);
		T x64 = x44 + xij3;
		T x65 = -x26;
		T x66 = x17*(x35 + x65);
		T x67 = x14*(x29 + x54);
		T x68 = x58 + xij1;
		T x69 = x17*(x46 + x60);
		T x70 = -x29;
		T x71 = x8*(x40 + x70);
		T x72 = x27 + xij2;
		T x73 = -x46;
		T x74 = x14*(x52 + x73);
		T x75 = x8*(x26 + x36);
		T x76 = x56 + xij2;
		T x77 = x17*(x26 + x62);
		T x78 = x14*(x38 + x70);
		T x79 = x24 + xij3;
		T x80 = x17*(x50 + x73);
		T x81 = x8*(x29 + x39);
		T x82 = x47 + xij1;
		T x83 = x14*(x46 + x51);
		T x84 = x8*(x34 + x65);
		u[i] = x19*x21 + x22*x23;
		u_xij[0 + i*3] = x30*(x25*x26 + x28*x29 + 2*x4 - 2*x7) + x33*(x25*x31 + x28*x32 + 4*x4 - 4*x7) + x42*(-x37 - x41) + x43*(-2*x37 - 2*x41);
		u_xij[1 + i*3] = x30*(2*x11 - 2*x13 + x29*x48 + x45*x46) + x33*(4*x11 - 4*x13 + x32*x48 + x45*x49) + x42*(-x53 - x55) + x43*(-2*x53 - 2*x55);
		u_xij[2 + i*3] = x30*(2*x15 - 2*x16 + x26*x59 + x46*x57) + x33*(4*x15 - 4*x16 + x31*x59 + x49*x57) + x42*(-x61 - x63) + x43*(-2*x61 - 2*x63);
		u_xik[0 + i*3] = x30*(x26*x64 + x29*x5) + x33*(x31*x64 + x32*x5) + x42*(-x66 - x67) + x43*(-2*x66 - 2*x67);
		u_xik[1 + i*3] = x30*(x29*x68 + x3*x46) + x33*(x3*x49 + x32*x68) + x42*(-x69 - x71) + x43*(-2*x69 - 2*x71);
		u_xik[2 + i*3] = x30*(x10*x26 + x46*x72) + x33*(x10*x31 + x49*x72) + x42*(-x74 - x75) + x43*(-2*x74 - 2*x75);
		u_xil[0 + i*3] = x30*(x26*x6 + x29*x76) + x33*(x31*x6 + x32*x76) + x42*(-x77 - x78) + x43*(-2*x77 - 2*x78);
		u_xil[1 + i*3] = x30*(x12*x29 + x46*x79) + x33*(x12*x32 + x49*x79) + x42*(-x80 - x81) + x43*(-2*x80 - 2*x81);
		u_xil[2 + i*3] = x30*(x1*x46 + x26*x82) + x33*(x1*x49 + x31*x82) + x42*(-x83 - x84) + x43*(-2*x83 - 2*x84);
		i += blockDim.x * gridDim.x;
	}
}

template <typename T> void gpuQuadrupletbGradient1(T *u, T *u_xij, T *u_xik, T *u_xil, T *xij, T *xik, T *xil, T *qi, T *qj, T *qk, T *ql, int *ti, int *tj, int *tk, int *tl, int *ai, int *aj, int *ak, int *al, T *mu, T *eta, int *kappa, int dim, int ncq, int nmu, int neta, int nkappa, int ng)
{
	int blockDim = 256;
	int gridDim = (ng * blockDim - 1) / blockDim;
	gridDim = (gridDim>1024)? 1024 : gridDim;
	kernelgpuQuadrupletbGradient1<<<gridDim, blockDim>>>(u, u_xij, u_xik, u_xil, xij, xik, xil, qi, qj, qk, ql, ti, tj, tk, tl, ai, aj, ak, al, mu, eta, kappa, dim, ncq, nmu, neta, nkappa, ng);
}


template <typename T> void gpuQuadrupletbGradient(T *u, T *u_xij, T *u_xik, T *u_xil, T *xij, T *xik, T *xil, T *qi, T *qj, T *qk, T *ql, int *ti, int *tj, int *tk, int *tl, int *ai, int *aj, int *ak, int *al, T *mu, T *eta, int *kappa, int dim, int ncq, int nmu, int neta, int nkappa, int ng, int potnum)
{
	if (potnum == 1)
		gpuQuadrupletbGradient1(u, u_xij, u_xik, u_xil, xij, xik, xil, qi, qj, qk, ql, ti, tj, tk, tl, ai, aj, ak, al, mu, eta, kappa, dim, ncq, nmu, neta, nkappa, ng);
}
template void gpuQuadrupletbGradient(double *, double *, double *, double *, double *, double *, double *, double *, double *, double *, double *, int *, int *, int *, int *, int *, int *, int *, int *, double *, double *, int*, int, int, int, int, int, int, int);
template void gpuQuadrupletbGradient(float *, float *, float *, float *, float *, float *, float *, float *, float *, float *, float *, int *, int *, int *, int *, int *, int *, int *, int *, float *, float *, int *, int, int, int, int, int, int, int);
