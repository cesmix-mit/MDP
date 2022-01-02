void cpuCart2Sphere(double *the, double *phi, double *r, double *x, double *y, double *z, int N)
{    
    //#pragma omp parallel for
    for (int i=0; i<N; i++) {
        r[i] = sqrt(x[i]*x[i] + y[i]*y[i] + z[i]*z[i]);
        the[i] = acos(z[i]/r[i]);
        phi[i] = atan2(y[i],x[i]);
    }
}

void cpuCart2SphereDeriv(double *the, double *phi, double *r, double *thex, double *they, double *thez, double *phix, double *phiy, double *phiz, double *rx, double *ry, double *rz, double *x, double *y, double *z, int N)
{    
    //#pragma omp parallel for
    for (int i=0; i<N; i++) {
        r[i] = sqrt(x[i]*x[i] + y[i]*y[i] + z[i]*z[i]);
        the[i] = acos(z[i]/r[i]);
        phi[i] = atan2(y[i],x[i]);
        
        double r2 = r[i]*r[i];
        double rxy = sqrt(x[i]*x[i] + y[i]*y[i]);
        double rxy2 = rxy*rxy;
        double rr2 = rxy*r2;

        rx[i] = x[i]/r[i];
        ry[i] = y[i]/r[i];
        rz[i] = z[i]/r[i];
        thex[i] = x[i]*z[i]/rr2;
        they[i] = y[i]*z[i]/rr2;
        thez[i] = -rxy/r2;
        phix[i] = -y[i]/rxy2;
        phiy[i] = x[i]/rxy2;
        phiz[i] = 0.0;        
    }
}

void cpuSphere2Cart(double *x, double *y, double *z, double *the, double *phi, double *r, int N)
{    
    //#pragma omp parallel for
    for (int i=0; i<N; i++) {
        x[i] = r[i]*sin(the[i])*cos(phi[i]);
        y[i] = r[i]*sin(the[i])*sin(phi[i]);
        z[i] = r[i]*cos(the[i]);
    }
}

void cpuEuler2Rotm(double *R11, double *R12, double *R13, double *R21, 
                double *R22, double *R23, double *R31, double *R32, double *R33, double *alpha, double *beta, double *gamma, int N)
{            
    //#pragma omp parallel for
    for (int i=0; i<N; i++) {
        double ca = cos(alpha[i]);
        double cb = cos(beta[i]);
        double cg = cos(gamma[i]);
        double sa = sin(alpha[i]);
        double sb = sin(beta[i]);
        double sg = sin(gamma[i]);
        
        R11[i] = ca*cg*cb - sa*sg;
        R12[i] = -ca*cb*sg - sa*cg;
        R13[i] = ca*sb;
        R21[i] = sa*cg*cb + ca*sg;
        R22[i] = -sa*cb*sg + ca*cg;
        R23[i] = sa*sb;
        R31[i] = -sb*cg;
        R32[i] = sb*sg;
        R33[i] = cb;
    }        
}

void cpuRotc(double *X, double *Y, double *Z, double *R, double *x, double *y, double *z, int N)
{            
    //#pragma omp parallel for
    for (int i=0; i<N; i++) {        
        X[i] = R[0]*x[i] + R[3]*y[i] + R[6]*z[i];
        Y[i] = R[1]*x[i] + R[4]*y[i] + R[7]*z[i];
        Z[i] = R[2]*x[i] + R[5]*y[i] + R[8]*z[i];
    }        
}

