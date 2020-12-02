#ifndef __CPUSPHERICALHARMONICS
#define __CPUSPHERICALHARMONICS

template <typename T> void cpuSphericalHarmonics(T *Ylmr, T *Ylmi, T *the, T *phi, 
                T *P, T *tmp, T *fac, T *C, T pi, int L, int N)
{                        
    //#pragma omp parallel for
    for (int i=0; i<N; i++) {                
        // l = 0
        Ylmr[i] = sqrt(1/(4*pi));
        Ylmi[i] = 0.0;
        
        T costhe = cos(the[i]);
        T a = -sin(the[i]);
        
        int l = 1;
        int k = 2*l+1;
        int m = 0;    
        P[0] = costhe;
        P[1] = a;
        
        m = 0;
        C[0] = sqrt(k*fac[l-m]/(4*pi*fac[l+m]));
        Ylmr[N+i] = C[0]*P[0];
        Ylmi[N+i] = 0.0;
        m = 1;
        C[1] = sqrt(k*fac[l-m]/(4*pi*fac[l+m]));
        Ylmr[2*N+i] = C[1]*cos(phi[i])*P[1];
        Ylmi[2*N+i] = C[1]*sin(phi[i])*P[1];                
        
        for (l=2; l<=L; l++) {                            
            
            T Pll = P[l-1];
            tmp[(l-1)] = Pll;
            P[(l-1)] = (2*(l-1)+1)*(costhe*Pll);
            P[l] = (2*(l-1)+1)*(a*Pll); 
            for (m=0; m<l-1; m++) {
                Pll = P[m];
                P[m] = ((2*(l-1)+1)*(costhe*Pll) - (l+m-1)*tmp[m])/(T (l-m));
                tmp[m] = Pll;
            }
            
            int k = 2*l+1;                
            for (m=0; m<=l; m++) {
                C[m] = sqrt(k*fac[l-m]/(4*pi*fac[l+m])); 
                Ylmr[(l*(l+1)/2+m)*N + i] = C[m]*(cos(m*phi[i])*P[m]);
                Ylmi[(l*(l+1)/2+m)*N + i] = C[m]*(sin(m*phi[i])*P[m]);
            }
        }
    }                        
}
template void cpuSphericalHarmonics(double*, double*, double*, double*, double*, double*, double*, double*, double, int, int);
template void cpuSphericalHarmonics(float*, float*, float*, float*, float*, float*, float*, float*, float, int, int);


template <typename T> void cpuSphericalHarmonicsDeriv(T *Ylmr, T *Ylmi, T *YlmrThe, T *YlmiThe, T *YlmrPhi, 
        T *YlmiPhi, T *the, T *phi, T *P, T *tmp, T *dP, T *dtmp, T *fac, T *C, T pi, int L, int N)
{                        
    //#pragma omp parallel for
    for (int i=0; i<N; i++) {                
        // l = 0
        Ylmr[i] = sqrt(1/(4*pi));
        Ylmi[i] = 0.0;
        YlmrThe[i] = 0.0;
        YlmrPhi[i] = 0.0;
        YlmiThe[i] = 0.0;        
        YlmiPhi[i] = 0.0;
        
        T costhe = cos(the[i]);
        T a = -sin(the[i]);
        T dcosthe = a;
        T da = -costhe;

        int l = 1;
        int k = 2*l+1;
        int m = 0;    
        P[0] = costhe;
        P[1] = a;
        dP[0] = dcosthe;    
        dP[1] = da;

        m = 0;
        C[0] = sqrt(k*fac[l-m]/(4*pi*fac[l+m]));
        Ylmr[N+i] = C[0]*P[0];
        Ylmi[N+i] = 0.0;
        YlmrThe[N+i] = C[0]*a;    
        YlmiThe[N+i] = 0.0;    
        YlmrPhi[N+i] = 0.0;          
        YlmiPhi[N+i] = 0.0;          

        m = 1;
        C[1] = sqrt(k*fac[l-m]/(4*pi*fac[l+m]));
        Ylmr[2*N+i] = C[1]*cos(phi[i])*P[1];
        Ylmi[2*N+i] = C[1]*sin(phi[i])*P[1];                
        YlmrThe[2*N+i] = C[1]*(cos(phi[i])*da);      
        YlmiThe[2*N+i] = C[1]*(sin(phi[i])*da);      
        YlmrPhi[2*N+i] = -C[1]*(sin(phi[i])*P[1]); 
        YlmiPhi[2*N+i] = C[1]*(cos(phi[i])*P[1]);  
        
        for (l=2; l<=L; l++) {                            
            
            T Pll = P[l-1];
            tmp[l-1] = Pll;
            P[l-1] = (2*(l-1)+1)*(costhe*Pll);
            P[l] = (2*(l-1)+1)*(a*Pll); 
            T dPll = dP[l-1];
            dtmp[l-1] = dPll;
            dP[l-1] = (2*(l-1)+1)*(costhe*dPll + dcosthe*Pll);   
            dP[l] = (2*(l-1)+1)*(a*dPll + da*Pll);             
            for (m=0; m<l-1; m++) {    
                Pll = P[m];
                P[m] = ((2*(l-1)+1)*(costhe*Pll) - (l+m-1)*tmp[(m)])/(T (l-m));
                tmp[m] = Pll;
                dPll = dP[m];
                dP[m] = ((2*(l-1)+1)*(costhe*dPll+dcosthe*Pll) - (l+m-1)*dtmp[m])/(T (l-m));
                dtmp[m] = dPll;
            }
                        
            int k = 2*l+1;  
            for (m=0; m<=l; m++) {
                int i1 = (l*(l+1)/2+m)*N + i;
                
                C[m] = sqrt(k*fac[l-m]/(4*pi*fac[l+m])); 
                Ylmr[i1] = C[m]*(cos(m*phi[i])*P[m]);
                Ylmi[i1] = C[m]*(sin(m*phi[i])*P[m]);
                YlmrThe[i1] = C[m]*(cos(m*phi[i])*dP[m]); 
                YlmiThe[i1] = C[m]*(sin(m*phi[i])*dP[m]); 
                YlmrPhi[i1] = -(m*C[m])*(sin(m*phi[i])*P[m]); 
                YlmiPhi[i1] = (m*C[m])*(cos(m*phi[i])*P[m]); 
            }
        }
    }                         
}
template void cpuSphericalHarmonicsDeriv(double*, double*, double*, double*, double*, double*, double*, double*, 
        double*, double*, double*, double*, double*, double*, double, int, int);
template void cpuSphericalHarmonicsDeriv(float*, float*, float*, float*, float*, float*, float*, float*, 
        float*, float*, float*, float*, float*, float*, float, int, int);

template <typename T> void cpuSphericalBessel(T *g, T *r, T *x0, T *f, int L, int K, int N)
{                        
    //#pragma omp parallel for
    for (int i=0; i<N; i++) {                
        int l, k;
        T x;        
        for (k=0; k<K; k++) {
            l = 0;
            x = x0[l*K+k]*r[i];                        
            g[l*N*K + k*N + i] = cos(x)/x;
            
            l = 1;
            x = x0[l*K+k]*r[i];                        
            g[l*N*K + k*N + i] = cos(x)/(x*x) + sin(x)/x;
            
            for (l=2; l<(L+1); l++) {
                x = x0[l*K+k]*r[i];                     
                f[0] = -cos(x)/x;
                f[1] = -cos(x)/(x*x) - sin(x)/x;
                for (int j=2; j<(l+1); j++) 
                    f[j] = ((2*(j+1)-3)/x)*f[(j-1)] - f[j-2];                 
                g[l*N*K + k*N + i] = -f[l];
            }            
        }        
    }
}
template void cpuSphericalBessel(double*, double*, double*, double*, int, int, int);
template void cpuSphericalBessel(float*, float*, float*, float*, int, int, int);

template <typename T> void cpuSphericalBesselDeriv(T *g, T *dg, T *r, T *x0, T *f, T *df, int L, int K, int N)
{                        
    //#pragma omp parallel for
    for (int i=0; i<N; i++) {                
        int l, k;
        T x;        
        for (k=0; k<K; k++) {
            l = 0;
            x = x0[l*K+k]*r[i];                        
            g[l*N*K + k*N + i] = cos(x)/x;
            dg[l*N*K + k*N + i] = x0[l*K+k]*(-cos(x)/(x*x) - sin(x)/x);
            
            l = 1;
            x = x0[l*K+k]*r[i];                        
            g[l*N*K + k*N + i] = cos(x)/(x*x) + sin(x)/x;
            dg[l*N*K + k*N + i] = x0[l*K+k]*(cos(x)/x - (2*cos(x))/(x*x*x) - (2*sin(x))/(x*x));
            
            for (l=2; l<(L+1); l++) {
                x = x0[l*K+k]*r[i];                     
                f[0] = -cos(x)/x;
                f[1] = -cos(x)/(x*x) - sin(x)/x;
                df[0] = x0[l*K+k]*(cos(x)/(x*x) + sin(x)/x);
                df[1] = x0[l*K+k]*((2*cos(x))/(x*x*x) - cos(x)/x  + (2*sin(x))/(x*x));
        
                for (int j=2; j<(l+1); j++) {
                    f[j] = ((2*(j+1)-3)/x)*f[j-1] - f[j-2];                 
                    df[j] = ((2*(j+1)-3)/x)*df[j-1] - x0[l*K+k]*((2*(j+1)-3)/(x*x))*f[j-1] - df[j-2];        
                }
                g[l*N*K + k*N + i] = -f[l];
                dg[l*N*K + k*N + i] = -df[l];
            }            
        }        
    }
}
template void cpuSphericalBesselDeriv(double*, double*, double*, double*, double*, double*, int, int, int);
template void cpuSphericalBesselDeriv(float*, float*, float*, float*, float*, float*, int, int, int);

template <typename T> void cpuRadialSphericalHarmonics(T *Sr, T *Si, T *g, T *Ylmr, T *Ylmi, int L, int K, int N)
{                        
    //#pragma omp parallel for
    for (int i=0; i<N; i++)                                
        for (int k=0; k<K; k++) 
            for (int l=0; l<(L+1); l++) 
                for (int m=0; m<(l+1); m++) {                  
                    Sr[(l*(l+1)/2+m)*N*K + k*N + i] = g[l*N*K + k*N + i]*Ylmr[(l*(l+1)/2+m)*N + i];
                    Si[(l*(l+1)/2+m)*N*K + k*N + i] = g[l*N*K + k*N + i]*Ylmi[(l*(l+1)/2+m)*N + i];                                
                }
}
template void cpuRadialSphericalHarmonics(double*, double*, double*, double*, double*, int, int, int);
template void cpuRadialSphericalHarmonics(float*, float*, float*, float*, float*, int, int, int);

template <typename T> void cpuWeightedRadialSphericalHarmonics(T *Sr, T *Si, T *g, T *Ylmr, T *Ylmi, T *w, int L, int K, int N)
{                        
    //#pragma omp parallel for
    for (int i=0; i<N; i++)                                
        for (int k=0; k<K; k++) 
            for (int l=0; l<(L+1); l++) 
                for (int m=0; m<(l+1); m++) {                  
                    Sr[(l*(l+1)/2+m)*N*K + k*N + i] = w[i]*g[l*N*K + k*N + i]*Ylmr[(l*(l+1)/2+m)*N + i];
                    Si[(l*(l+1)/2+m)*N*K + k*N + i] = w[i]*g[l*N*K + k*N + i]*Ylmi[(l*(l+1)/2+m)*N + i];                                
                }
}
template void cpuWeightedRadialSphericalHarmonics(double*, double*, double*, double*, double*, double*, int, int, int);
template void cpuWeightedRadialSphericalHarmonics(float*, float*, float*, float*, float*, float*, int, int, int);

template <typename T> void cpuRadialSphericalHarmonicsSum(T *Sr, T *Si, T *g, T *Ylmr, T *Ylmi, int L, int K, int N)
{                        
    //#pragma omp parallel for                
    for (int k=0; k<K; k++) 
        for (int l=0; l<(L+1); l++) 
            for (int m=0; m<(l+1); m++) {                  
                Sr[(l*(l+1)/2+m)*K + k] = 0.0;
                Si[(l*(l+1)/2+m)*K + k] = 0.0;
                for (int i=0; i<N; i++) {                   
                    Sr[(l*(l+1)/2+m)*K + k] += g[l*N*K + k*N + i]*Ylmr[(l*(l+1)/2+m)*N + i];
                    Si[(l*(l+1)/2+m)*K + k] += g[l*N*K + k*N + i]*Ylmi[(l*(l+1)/2+m)*N + i];                                
                }
            }
}
template void cpuRadialSphericalHarmonicsSum(double*, double*, double*, double*, double*, int, int, int);
template void cpuRadialSphericalHarmonicsSum(float*, float*, float*, float*, float*, int, int, int);

template <typename T> void cpuWeightedRadialSphericalHarmonicsSum(T *Sr, T *Si, T *g, T *Ylmr, T *Ylmi, T *w, int L, int K, int N)
{                        
    //#pragma omp parallel for                
    for (int k=0; k<K; k++) 
        for (int l=0; l<(L+1); l++) 
            for (int m=0; m<(l+1); m++) {                  
                Sr[(l*(l+1)/2+m)*K + k] = 0.0;
                Si[(l*(l+1)/2+m)*K + k] = 0.0;
                for (int i=0; i<N; i++) {                   
                    Sr[(l*(l+1)/2+m)*K + k] += w[i]*g[l*N*K + k*N + i]*Ylmr[(l*(l+1)/2+m)*N + i];
                    Si[(l*(l+1)/2+m)*K + k] += w[i]*g[l*N*K + k*N + i]*Ylmi[(l*(l+1)/2+m)*N + i];                                
                }
            }
}
template void cpuWeightedRadialSphericalHarmonicsSum(double*, double*, double*, double*, double*, double*, int, int, int);
template void cpuWeightedRadialSphericalHarmonicsSum(float*, float*, float*, float*, float*, float*, int, int, int);

template <typename T> void cpuRadialSphericalHarmonicsDeriv(T *Srx, T *Six, T *Sry, T *Siy, T *Srz, T *Siz, 
        T *x, T *y, T *z, T *g, T *Ylmr, T *Ylmi, T *gR, T *YlmrThe, T *YlmiThe, T *YlmrPhi, T *YlmiPhi, 
        int L, int K, int N)
{                        
    //#pragma omp parallel for
    for (int i=0; i<N; i++) {                
        T r = sqrt(x[i]*x[i] + y[i]*y[i] + z[i]*z[i]);
        
        T r2 = r*r;
        T rxy = sqrt(x[i]*x[i] + y[i]*y[i]);
        T rxy2 = rxy*rxy;
        T rr2 = rxy*r2;

        T Rx = x[i]/r;
        T Ry = y[i]/r;
        T Rz = z[i]/r;
        T Thex = x[i]*z[i]/rr2;
        T They = y[i]*z[i]/rr2;
        T Thez = -rxy/r2;
        T Phix = -y[i]/rxy2;
        T Phiy = x[i]/rxy2;
        T Phiz = 0.0;        
        
        for (int k=0; k<K; k++) {
            for (int l=0; l<(L+1); l++) {
                int i0 = l*N*K + k*N + i;
                T gx = gR[i0]*Rx;
                T gy = gR[i0]*Ry;
                T gz = gR[i0]*Rz;
                for (int m=0; m<(l+1); m++) {                  
                    int i1 = (l*(l+1)/2+m)*N*K + k*N + i;
                    int i2 = (l*(l+1)/2+m)*N + i;
                    
                    T Ylmrx = YlmrThe[i2]*Thex + YlmrPhi[i2]*Phix;
                    T Ylmry = YlmrThe[i2]*They + YlmrPhi[i2]*Phiy;
                    T Ylmrz = YlmrThe[i2]*Thez + YlmrPhi[i2]*Phiz;
                    T Ylmix = YlmiThe[i2]*Thex + YlmiPhi[i2]*Phix;
                    T Ylmiy = YlmiThe[i2]*They + YlmiPhi[i2]*Phiy;
                    T Ylmiz = YlmiThe[i2]*Thez + YlmiPhi[i2]*Phiz;                    
                    
                    Srx[i1] = gx*Ylmr[i2] + g[i0]*Ylmrx;
                    Sry[i1] = gy*Ylmr[i2] + g[i0]*Ylmry;
                    Srz[i1] = gz*Ylmr[i2] + g[i0]*Ylmrz;
                    Six[i1] = gx*Ylmi[i2] + g[i0]*Ylmix;
                    Siy[i1] = gy*Ylmi[i2] + g[i0]*Ylmiy;
                    Siz[i1] = gz*Ylmi[i2] + g[i0]*Ylmiz;
                }
            }
        }
    }
}
template void cpuRadialSphericalHarmonicsDeriv(double*, double*, double*, double*, double*, double*, double*, 
        double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, int, int, int);
template void cpuRadialSphericalHarmonicsDeriv(float*, float*, float*, float*, float*, float*, float*, 
        float*, float*, float*, float*, float*, float*, float*, float*, float*, float*, int, int, int);

template <typename T> void cpuWeightedRadialSphericalHarmonicsDeriv(T *Srx, T *Six, T *Sry, T *Siy, T *Srz, T *Siz, 
        T *x, T *y, T *z, T *g, T *Ylmr, T *Ylmi, T *gR, T *YlmrThe, T *YlmiThe, T *YlmrPhi, T *YlmiPhi, T *w, 
        int L, int K, int N)
{                        
    //#pragma omp parallel for
    for (int i=0; i<N; i++) {                
        T r = sqrt(x[i]*x[i] + y[i]*y[i] + z[i]*z[i]);
        T r2 = r*r;
        T rxy = sqrt(x[i]*x[i] + y[i]*y[i]);
        T rxy2 = rxy*rxy;
        T rr2 = rxy*r2;

        T Rx = x[i]/r;
        T Ry = y[i]/r;
        T Rz = z[i]/r;
        T Thex = x[i]*z[i]/rr2;
        T They = y[i]*z[i]/rr2;
        T Thez = -rxy/r2;
        T Phix = -y[i]/rxy2;
        T Phiy = x[i]/rxy2;
        T Phiz = 0.0;        
        
        for (int k=0; k<K; k++) {
            for (int l=0; l<(L+1); l++) {
                int i0 = l*N*K + k*N + i;
                T gx = w[i]*gR[i0]*Rx;
                T gy = w[i]*gR[i0]*Ry;
                T gz = w[i]*gR[i0]*Rz;
                T wg = w[i]*g[i0];
                for (int m=0; m<(l+1); m++) {                  
                    int i1 = (l*(l+1)/2+m)*N*K + k*N + i;
                    int i2 = (l*(l+1)/2+m)*N + i;
                                        
                    T Ylmrx = YlmrThe[i2]*Thex + YlmrPhi[i2]*Phix;
                    T Ylmry = YlmrThe[i2]*They + YlmrPhi[i2]*Phiy;
                    T Ylmrz = YlmrThe[i2]*Thez + YlmrPhi[i2]*Phiz;
                    T Ylmix = YlmiThe[i2]*Thex + YlmiPhi[i2]*Phix;
                    T Ylmiy = YlmiThe[i2]*They + YlmiPhi[i2]*Phiy;
                    T Ylmiz = YlmiThe[i2]*Thez + YlmiPhi[i2]*Phiz;                    
                    
                    Srx[i1] = gx*Ylmr[i2] + wg*Ylmrx;
                    Sry[i1] = gy*Ylmr[i2] + wg*Ylmry;
                    Srz[i1] = gz*Ylmr[i2] + wg*Ylmrz;
                    Six[i1] = gx*Ylmi[i2] + wg*Ylmix;
                    Siy[i1] = gy*Ylmi[i2] + wg*Ylmiy;
                    Siz[i1] = gz*Ylmi[i2] + wg*Ylmiz;
                }
            }
        }
    }
}
template void cpuWeightedRadialSphericalHarmonicsDeriv(double*, double*, double*, double*, double*, 
        double*, double*, double*, double*, double*, double*, double*, double*, double*, 
        double*, double*, double*, double*, int, int, int);
template void cpuWeightedRadialSphericalHarmonicsDeriv(float*, float*, float*, float*, float*, 
        float*, float*, float*, float*, float*, float*, float*, float*, float*, float*, 
        float*, float*, float*, int, int, int);


#endif


