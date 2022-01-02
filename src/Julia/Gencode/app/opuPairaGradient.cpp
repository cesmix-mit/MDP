template <typename T> void opuPairaGradient1(T *u, T *u_xij, T *xij, T *qi, T *qj, int *ti, int *tj, int *ai, int *aj, T *mu, T *eta, int *kappa, int dim, int ncq, int nmu, int neta, int nkappa, int ng)
{
	for (int i = 0; i <ng; i++) {
		T mu1 = mu[0];
		T mu2 = mu[1];
		T eta1 = eta[0];
		T xij1 = xij[0 + i*3];
		T xij2 = xij[1 + i*3];
		T xij3 = xij[2 + i*3];
		T x0 = pow(xij1, 2) + pow(xij2, 2) + pow(xij3, 2);
		T x1 = 12*mu1/pow(x0, 7);
		T x2 = 6*mu2/pow(x0, 4);
		u[i] = mu1/pow(x0, 6) - mu2/pow(x0, 3) + mu2/pow(eta1, 6) - mu1/pow(eta1, 12);
		u_xij[0 + i*3] = -x1*xij1 + x2*xij1;
		u_xij[1 + i*3] = -x1*xij2 + x2*xij2;
		u_xij[2 + i*3] = -x1*xij3 + x2*xij3;
	}
}

template <typename T> void opuPairaGradient2(T *u, T *u_xij, T *xij, T *qi, T *qj, int *ti, int *tj, int *ai, int *aj, T *mu, T *eta, int *kappa, int dim, int ncq, int nmu, int neta, int nkappa, int ng)
{
	for (int i = 0; i <ng; i++) {
		T mu3 = mu[2];
		T mu4 = mu[3];
		T eta2 = eta[1];
		T xij1 = xij[0 + i*3];
		T xij2 = xij[1 + i*3];
		T xij3 = xij[2 + i*3];
		T qi1 = qi[0 + i*1];
		T qj1 = qj[0 + i*1];
		T x0 = pow(eta2, -7);
		T x1 = pow(xij1, 2) + pow(xij2, 2) + pow(xij3, 2);
		T x2 = pow(eta2, -6);
		T x3 = pow(eta2, -5);
		T x4 = pow(x1, 5.0/2.0);
		T x5 = pow(eta2, -4);
		T x6 = pow(x1, 2);
		T x7 = 20*x0*pow(x1, 7.0/2.0) - 70*pow(x1, 3)*x2 + 84*x3*x4 - 35*x5*x6 + 1;
		T x8 = pow(x1, 3.0/2.0);
		T x9 = x8 + pow(mu3, -3);
		T x10 = mu4*qi1*qj1*pow(x9, -0.33333333333333331);
		T x11 = 140*xij1;
		T x12 = x1*x5;
		T x13 = x0*x4;
		T x14 = 420*xij1;
		T x15 = x2*x6;
		T x16 = x3*x8;
		T x17 = 1.0*mu4*qi1*qj1*sqrt(x1)*x7*pow(x9, -1.3333333333333333);
		T x18 = 140*xij2;
		T x19 = 420*xij2;
		T x20 = 140*xij3;
		T x21 = 420*xij3;
		u[i] = x10*x7;
		u_xij[0 + i*3] = x10*(-x11*x12 + x11*x13 - x14*x15 + x14*x16) - x17*xij1;
		u_xij[1 + i*3] = x10*(-x12*x18 + x13*x18 - x15*x19 + x16*x19) - x17*xij2;
		u_xij[2 + i*3] = x10*(-x12*x20 + x13*x20 - x15*x21 + x16*x21) - x17*xij3;
	}
}

template <typename T> void opuPairaGradient3(T *u, T *u_xij, T *xij, T *qi, T *qj, int *ti, int *tj, int *ai, int *aj, T *mu, T *eta, int *kappa, int dim, int ncq, int nmu, int neta, int nkappa, int ng)
{
	for (int i = 0; i <ng; i++) {
		T mu5 = mu[4];
		T mu6 = mu[5];
		T mu7 = mu[6];
		T mu8 = mu[7];
		T mu9 = mu[8];
		T eta3 = eta[2];
		T xij1 = xij[0 + i*3];
		T xij2 = xij[1 + i*3];
		T xij3 = xij[2 + i*3];
		T x0 = pow(eta3, -7);
		T x1 = pow(xij1, 2) + pow(xij2, 2) + pow(xij3, 2);
		T x2 = pow(eta3, -6);
		T x3 = pow(eta3, -5);
		T x4 = pow(x1, 5.0/2.0);
		T x5 = pow(eta3, -4);
		T x6 = pow(x1, 2);
		T x7 = 20*x0*pow(x1, 7.0/2.0) - 70*pow(x1, 3)*x2 + 84*x3*x4 - 35*x5*x6 + 1;
		T x8 = pow(x1, (1.0/2.0)*mu9);
		T x9 = x8 + pow(1.0/mu5, mu9);
		T x10 = pow(x9, 1.0/mu9)/mu8;
		T x11 = mu7*(1 - x10);
		T x12 = exp(x11);
		T x13 = exp(0.5*x11);
		T x14 = mu6*(x12 - 2*x13);
		T x15 = 140*xij1;
		T x16 = x1*x5;
		T x17 = x0*x4;
		T x18 = 420*xij1;
		T x19 = x2*x6;
		T x20 = pow(x1, 3.0/2.0)*x3;
		T x21 = 1.0*mu7*x10*x8/(x1*x9);
		T x22 = x21*xij1;
		T x23 = mu6*x7;
		T x24 = 140*xij2;
		T x25 = 420*xij2;
		T x26 = x21*xij2;
		T x27 = 140*xij3;
		T x28 = 420*xij3;
		T x29 = x21*xij3;
		u[i] = x14*x7;
		u_xij[0 + i*3] = x14*(-x15*x16 + x15*x17 - x18*x19 + x18*x20) + x23*(-x12*x22 + x13*x22);
		u_xij[1 + i*3] = x14*(-x16*x24 + x17*x24 - x19*x25 + x20*x25) + x23*(-x12*x26 + x13*x26);
		u_xij[2 + i*3] = x14*(-x16*x27 + x17*x27 - x19*x28 + x20*x28) + x23*(-x12*x29 + x13*x29);
	}
}


template <typename T> void opuPairaGradient(T *u, T *u_xij, T *xij, T *qi, T *qj, int *ti, int *tj, int *ai, int *aj, T *mu, T *eta, int *kappa, int dim, int ncq, int nmu, int neta, int nkappa, int ng, int potnum)
{
	if (potnum == 1)
		opuPairaGradient1(u, u_xij, xij, qi, qj, ti, tj, ai, aj, mu, eta, kappa, dim, ncq, nmu, neta, nkappa, ng);
	else if (potnum == 2)
		opuPairaGradient2(u, u_xij, xij, qi, qj, ti, tj, ai, aj, mu, eta, kappa, dim, ncq, nmu, neta, nkappa, ng);
	else if (potnum == 3)
		opuPairaGradient3(u, u_xij, xij, qi, qj, ti, tj, ai, aj, mu, eta, kappa, dim, ncq, nmu, neta, nkappa, ng);
}
template void opuPairaGradient(double *, double *, double *, double *, double *, int *, int *, int *, int *, double *, double *, int*, int, int, int, int, int, int, int);
template void opuPairaGradient(float *, float *, float *, float *, float *, int *, int *, int *, int *, float *, float *, int *, int, int, int, int, int, int, int);
