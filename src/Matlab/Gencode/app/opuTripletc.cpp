template <typename T> void opuTripletc1(T *u, T *xij, T *xik, T *qi, T *qj, T *qk, int *ti, int *tj, int *tk, int *ai, int *aj, int *ak, T *mu, T *eta, int *kappa, int dim, int ncq, int nmu, int neta, int nkappa, int ng)
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
		T xij1 = xij[0 + i*3];
		T xij2 = xij[1 + i*3];
		T xij3 = xij[2 + i*3];
		T xik1 = xik[0 + i*3];
		T xik2 = xik[1 + i*3];
		T xik3 = xik[2 + i*3];
		T t2 = mu5*mu5;
		T t3 = xij1*xij1;
		T t4 = xij2*xij2;
		T t5 = xij3*xij3;
		T t6 = xik1*xik1;
		T t7 = xik2*xik2;
		T t8 = xik3*xik3;
		T t9 = t3+t4+t5;
		T t10 = t6+t7+t8;
		T t11 = sqrt(t9);
		if (tk[i] == 1) 
			u[i] = -mu1*exp(-mu2*t11);
		if (tk[i] == 2) 
			u[i] = mu4*exp(pow(mu3,mu8)*pow(t11-sqrt(t10),mu8))*(1.0/(mu6*mu6)*t2-t2/(pow(mu7-(1.0/sqrt(t10)*(xij1*xik1+xij2*xik2+xij3*xik3))/t11,2.0)+mu6*mu6)+1.0);
	}
}


template <typename T> void opuTripletcPair1(T *u, T *xij, T *qi, T *qj, int *ti, int *tj, int *ai, int *aj, T *mu, T *eta, int *kappa, int dim, int ncq, int nmu, int neta, int nkappa, int ng)
{
	for (int i = 0; i <ng; i++) {
		T mu3 = mu[2];
		T mu4 = mu[3];
		T mu5 = mu[4];
		T mu6 = mu[5];
		T mu7 = mu[6];
		T mu8 = mu[7];
		T xij1 = xij[0 + i*3];
		T xij2 = xij[1 + i*3];
		T xij3 = xij[2 + i*3];
		T t2 = mu5*mu5;
		T t3 = xij1*xij1;
		T t4 = xij2*xij2;
		T t5 = xij3*xij3;
		T t6 = xik1*xik1;
		T t7 = xik2*xik2;
		T t8 = xik3*xik3;
		T t9 = t3+t4+t5;
		T t10 = t6+t7+t8;
		u[0 + i*1] = mu4*exp(pow(mu3,mu8)*pow(sqrt(t9)-sqrt(t10),mu8))*(1.0/(mu6*mu6)*t2-t2/(pow(mu7-1.0/sqrt(t9)*1.0/sqrt(t10)*(xij1*xik1+xij2*xik2+xij3*xik3),2.0)+mu6*mu6)+1.0);
	}
}

template <typename T> void opuTripletcDensity1(T *u, T *rho, T *mu, T *eta, int *kappa, int nrho, int nmu, int neta, int nkappa, int ng)
{
	for (int i = 0; i <ng; i++) {
		T mu1 = mu[0];
		T mu9 = mu[8];
		T mu10 = mu[9];
		T rho1 = rho[0 + i*1];
		u[0 + i*1] = pow(pow(mu10,mu9)*pow(rho1,mu9)+1.0,(-1.0/2.0)/mu9);
	}
}


template <typename T> void opuTripletc(T *u, T *xij, T *xik, T *qi, T *qj, T *qk, int *ti, int *tj, int *tk, int *ai, int *aj, int *ak, T *mu, T *eta, int *kappa, int dim, int ncq, int nmu, int neta, int nkappa, int ng, int potnum)
{
	if (potnum == 1)
		opuTripletc1(u, xij, xik, qi, qj, qk, ti, tj, tk, ai, aj, ak, mu, eta, kappa, dim, ncq, nmu, neta, nkappa, ng);
}
template void opuTripletc(double *, double *, double *, double *, double *, double *, int *, int *, int *, int *, int *, int *, double *, double *, int*, int, int, int, int, int, int, int);
template void opuTripletc(float *, float *, float *, float *, float *, float *, int *, int *, int *, int *, int *, int *, float *, float *, int *, int, int, int, int, int, int, int);

template <typename T> void opuTripletcPair(T *u, T *xij, T *qi, T *qj, int *ti, int *tj, int *ai, int *aj, T *mu, T *eta, int *kappa, int dim, int ncq, int nmu, int neta, int nkappa, int ng, int potnum)
{
	if (potnum == 1)
		opuTripletcPair1(u, xij, qi, qj, ti, tj, ai, aj, mu, eta, kappa, dim, ncq, nmu, neta, nkappa, ng);
}
template void opuTripletcPair(double *, double *, double *, double *, int *, int *, int *, int *, double *, double *, int*, int, int, int, int, int, int, int);
template void opuTripletcPair(float *, float *, float *, float *, int *, int *, int *, int *, float *, float *, int *, int, int, int, int, int, int, int);

template <typename T> void opuTripletcDensity(T *u, T *rho, T *mu, T *eta, int *kappa, int nrho, int nmu, int neta, int nkappa, int ng, int potnum)
{
	if (potnum == 1)
		opuTripletcDensity1(u, rho, mu, eta, kappa, nrho, nmu, neta, nkappa, ng);
}
template void opuTripletcDensity(double *, double *, double *, double *, int*, int, int, int, int, int, int);
template void opuTripletcDensity(float *, float *, float *, float *, int *, int, int, int, int, int, int);
