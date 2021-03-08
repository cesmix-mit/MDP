template <typename T> void opuPaira1(T *u, T *xij, T *qi, T *qj, int *ti, int *tj, int *ai, int *aj, T *mu, T *eta, int *kappa, int dim, int ncq, int nmu, int neta, int nkappa, int ng)
{
	for (int i = 0; i <ng; i++) {
		T mu1 = mu[0];
		T mu2 = mu[1];
		T mu3 = mu[2];
		T mu4 = mu[3];
		T mu5 = mu[4];
		T mu6 = mu[5];
		T mu7 = mu[6];
		T mu8 = mu[7];
		T mu9 = mu[8];
		T mu10 = mu[9];
		int kappa1 = kappa[0];
		int kappa2 = kappa[1];
		int kappa3 = kappa[2];
		int kappa4 = kappa[3];
		int kappa5 = kappa[4];
		int kappa6 = kappa[5];
		int kappa7 = kappa[6];
		int kappa8 = kappa[7];
		int kappa9 = kappa[8];
		int kappa10 = kappa[9];
		T xij1 = xij[0 + i*3];
		T xij2 = xij[1 + i*3];
		T xij3 = xij[2 + i*3];
		T t2 = xij1*xij1;
		T t3 = xij2*xij2;
		T t4 = xij3*xij3;
		T t5 = t2+t3+t4;
		u[i] = 1.0/(mu1*mu1*mu1*mu1*mu1*mu1)*mu3-1.0/pow(mu1,1.2E+1)*mu2-mu3*1.0/(t5*t5*t5)+mu2*1.0/(t5*t5*t5*t5*t5*t5);
	}
}

template <typename T> void opuPaira2(T *u, T *xij, T *qi, T *qj, int *ti, int *tj, int *ai, int *aj, T *mu, T *eta, int *kappa, int dim, int ncq, int nmu, int neta, int nkappa, int ng)
{
	for (int i = 0; i <ng; i++) {
		T mu1 = mu[0];
		T mu2 = mu[1];
		T mu3 = mu[2];
		T mu4 = mu[3];
		T mu5 = mu[4];
		T mu6 = mu[5];
		T mu7 = mu[6];
		T mu8 = mu[7];
		T mu9 = mu[8];
		T mu10 = mu[9];
		int kappa1 = kappa[0];
		int kappa2 = kappa[1];
		int kappa3 = kappa[2];
		int kappa4 = kappa[3];
		int kappa5 = kappa[4];
		int kappa6 = kappa[5];
		int kappa7 = kappa[6];
		int kappa8 = kappa[7];
		int kappa9 = kappa[8];
		int kappa10 = kappa[9];
		T xij1 = xij[0 + i*3];
		T xij2 = xij[1 + i*3];
		T xij3 = xij[2 + i*3];
		T qi1 = qi[0 + i*1];
		T qj1 = qj[0 + i*1];
		T t2 = xij1*xij1;
		T t3 = xij2*xij2;
		T t4 = xij3*xij3;
		T t5 = t2+t3+t4;
		u[i] = mu6*qi1*qj1*1.0/pow(1.0/(mu5*mu5*mu5)+pow(t5,3.0/2.0),1.0/3.0)*(1.0/(mu4*mu4*mu4*mu4)*(t5*t5)*-3.5E+1-1.0/(mu4*mu4*mu4*mu4*mu4*mu4)*(t5*t5*t5)*7.0E+1+1.0/(mu4*mu4*mu4*mu4*mu4)*pow(t5,5.0/2.0)*8.4E+1+1.0/(mu4*mu4*mu4*mu4*mu4*mu4*mu4)*pow(t5,7.0/2.0)*2.0E+1+1.0);
	}
}

template <typename T> void opuPaira(T *u, T *xij, T *qi, T *qj, int *ti, int *tj, int *ai, int *aj, T *mu, T *eta, int *kappa, int dim, int ncq, int nmu, int neta, int nkappa, int ng, int potnum)
{
	if (potnum == 1)
		opuPaira1(u, xij, qi, qj, ti, tj, ai, aj, mu, eta, kappa, dim, ncq, nmu, neta, nkappa, ng);
	else if (potnum == 2)
		opuPaira2(u, xij, qi, qj, ti, tj, ai, aj, mu, eta, kappa, dim, ncq, nmu, neta, nkappa, ng);
}
template void opuPaira(double *, double *, double *, double *, int *, int *, int *, int *, double *, double *, int*, int, int, int, int, int, int, int);
template void opuPaira(float *, float *, float *, float *, int *, int *, int *, int *, float *, float *, int *, int, int, int, int, int, int, int);
