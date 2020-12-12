// Example usage of Enzyme to differentiate cpuSphericalHarmonics function
// Example ways to run (setting path as appropriate:
//  1)
//     clang++ -ffast-math -O3 enzymeSphericalHarmoics.cpp -Xclang -load -Xclang /path/to/ClangEnzyme-7.so -lm
//     ./a.out
//  2)
//     clang++ -ffast-math -O3 enzymeSphericalHarmoics. -S -emit-llvm | opt - -load=/path/to/LLVMEnzyme-7.so -enzyme -O3 -S | clang -x ir - -lm
//     ./a.out
// The first method only runs optimization prior to Enzyme, whereas the second
// method runs before and after Enzyme

#include <stdio.h>
#include <math.h>

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

template <typename T> void cpuRadialSphericalHarmonicsPower(T *p, T *Sr, T *Si, int *indk, int L, int K)
{
    int K2 = K*(K+1)/2;

    //#pragma omp parallel for
    for (int k=0; k<K2; k++) {
        int j = indk[k];
        int i = indk[K2+k];
        for (int l=0; l<(L+1); l++) {
            int m = 0;
            int i1 = (l*(l+1)/2+m)*K + i;
            int j1 = (l*(l+1)/2+m)*K + j;
            int l1 = k*(L+1) + l;
            p[l1] = Sr[i1]*Sr[j1] + Si[i1]*Si[j1];
            for (m=1; m<(l+1); m++) {
                i1 = (l*(l+1)/2+m)*K + i;
                j1 = (l*(l+1)/2+m)*K + j;
                p[l1] += 2.0*(Sr[i1]*Sr[j1] + Si[i1]*Si[j1]);
            }
        }
    }
}
template void cpuRadialSphericalHarmonicsPower(double*, double*, double*, int*, int, int);
template void cpuRadialSphericalHarmonicsPower(float*, float*, float*, int*, int, int);

template <typename T> void cpuRadialSphericalHarmonicsPowerDeriv(T *px, T *py, T *pz, T *Sr, T *Si,
        T *Srx, T *Six, T *Sry, T *Siy, T *Srz, T *Siz, int *indk, int L, int K, int N)
{
    int K2 = K*(K+1)/2;

    //#pragma omp parallel for
    for (int n=0; n<N; n++)
        for (int k=0; k<K2; k++) {
            int j = indk[k];
            int i = indk[K2+k];
            for (int l=0; l<(L+1); l++) {
                int m = 0;
                int i1 = (l*(l+1)/2+m)*K + i;
                int j1 = (l*(l+1)/2+m)*K + j;
                int l1 = k*(L+1) + l;
                int i2 = (l*(l+1)/2+m)*N*K + i*N + n;
                int j2 = (l*(l+1)/2+m)*N*K + j*N + n;
                //p[l1] = Sr[i1]*Sr[j1] + Si[i1]*Si[j1];
                px[l1*N + n] = Sr[i1]*Srx[j2] + Srx[i2]*Sr[j1] + Si[i1]*Six[j2] + Six[i2]*Si[j1];
                py[l1*N + n] = Sr[i1]*Sry[j2] + Sry[i2]*Sr[j1] + Si[i1]*Siy[j2] + Siy[i2]*Si[j1];
                pz[l1*N + n] = Sr[i1]*Srz[j2] + Srz[i2]*Sr[j1] + Si[i1]*Siz[j2] + Siz[i2]*Si[j1];
                for (m=1; m<(l+1); m++) {
                    i1 = (l*(l+1)/2+m)*K + i;
                    j1 = (l*(l+1)/2+m)*K + j;
                    i2 = (l*(l+1)/2+m)*N*K + i*N + n;
                    j2 = (l*(l+1)/2+m)*N*K + j*N + n;
                    //p[l1] += 2.0*(Sr[i1]*Sr[j1] + Si[i1]*Si[j1]);
                    px[l1*N + n] += 2*(Sr[i1]*Srx[j2] + Srx[i2]*Sr[j1] + Si[i1]*Six[j2] + Six[i2]*Si[j1]);
                    py[l1*N + n] += 2*(Sr[i1]*Sry[j2] + Sry[i2]*Sr[j1] + Si[i1]*Siy[j2] + Siy[i2]*Si[j1]);
                    pz[l1*N + n] += 2*(Sr[i1]*Srz[j2] + Srz[i2]*Sr[j1] + Si[i1]*Siz[j2] + Siz[i2]*Si[j1]);
                }
            }
        }
}
template void cpuRadialSphericalHarmonicsPowerDeriv(double*, double*, double*, double*, double*, double*,
        double*, double*, double*, double*, double*, int*, int, int, int);
template void cpuRadialSphericalHarmonicsPowerDeriv(float*, float*, float*, float*, float*, float*,
        float*, float*, float*, float*, float*, int*, int, int, int);

template <typename T> void cpuRadialSphericalHarmonicsBispectrum(T *b, T *Sr, T *Si, T*cg, int *indk, int *indl,
        int *indm, int *rowm, int K, int L, int M)
{
    int K2 = K*(K+1)/2;

    //#pragma omp parallel for
    for (int k=0; k<K2; k++) {
        int k2 = indk[k];
        int k1 = indk[K2+k];
        for (int i=0; i<L; i++) {
            int l2 = indl[i];
            int l1 = indl[L+i];
            int l = indl[2*L+i];
            T tmp = 0;
            int nm = rowm[i+1]-rowm[i];
            for (int j = 0; j<nm; j++) {
                int m2 = indm[rowm[i]+j];
                int m1 = indm[M+rowm[i]+j];
                int m = indm[2*M + rowm[i]+j];

                int mm = (m>=0) ? m : -m;
                int mm1 = (m1>=0) ? m1 : -m1;
                int mm2 = (m2>=0) ? m2 : -m2;
                int i1 = (l*(l+1)/2+mm)*K + k1;
                int i2 = (l1*(l1+1)/2+mm1)*K + k2;
                int i3 = (l2*(l2+1)/2+mm2)*K + k2;

                T a1, b1, c1, a2, b2, c2, a3, b3, c3;
                c1 = ( mm % 2 == 0) ?  1.0 : -1.0;
                c2 = ( mm1 % 2 == 0) ?  1.0 : -1.0;
                c3 = ( mm2 % 2 == 0) ?  1.0 : -1.0;

                if (m>=0) {
                    a1 = Sr[i1];
                    b1 = Si[i1];
                }
                else {
                    a1 =-c1*Sr[i1];
                    b1 = c1*Si[i1];
                }
                if (m1>=0) {
                    a2 = Sr[i2];
                    b2 = Si[i2];
                }
                else {
                    a2 = -c2*Sr[i2];
                    b2 =  c2*Si[i2];
                }
                if (m2>=0) {
                    a3 = Sr[i3];
                    b3 = Si[i3];
                }
                else {
                    a3 = -c3*Sr[i3];
                    b3 =  c3*Si[i3];
                }
                tmp = tmp + cg[rowm[i]+j]*(a1*a2*a3 + a2*b1*b3 + a3*b1*b2 - a1*b2*b3);
            }
            b[i+k*L] = tmp;
        }
    }
}
template void cpuRadialSphericalHarmonicsBispectrum(double*, double*, double*, double*, int*, int*, int*, int*, int, int, int);
template void cpuRadialSphericalHarmonicsBispectrum(float*, float*, float*, float*, int*, int*, int*, int*, int, int, int);

template <typename T> void cpuRadialSphericalHarmonicsBispectrumDeriv(T *bx, T *by, T *bz, T *Sr, T *Si,
        T *Srx, T *Six, T *Sry, T *Siy, T *Srz, T *Siz, T*cg, int *indk, int *indl,
        int *indm, int *rowm, int L, int K, int M, int N)
{
    int K2 = K*(K+1)/2;

    //#pragma omp parallel for
    for (int n=0; n<N; n++)
        for (int k=0; k<K2; k++) {
            int k2 = indk[k];
            int k1 = indk[K2+k];
            for (int i=0; i<L; i++) {
                int l2 = indl[i];
                int l1 = indl[L+i];
                int l = indl[2*L+i];
                int ii = (i+k*L)*N + n;
                bx[ii] = 0.0;
                by[ii] = 0.0;
                bz[ii] = 0.0;
                int nm = rowm[i+1]-rowm[i];
                for (int j = 0; j<nm; j++) {
                    int m2 = indm[rowm[i]+j];
                    int m1 = indm[M+rowm[i]+j];
                    int m = indm[2*M + rowm[i]+j];
                    int mm = (m>=0) ? m : -m;
                    int mm1 = (m1>=0) ? m1 : -m1;
                    int mm2 = (m2>=0) ? m2 : -m2;
                    int i1 = (l*(l+1)/2+mm)*K + k1;
                    int i2 = (l1*(l1+1)/2+mm1)*K + k2;
                    int i3 = (l2*(l2+1)/2+mm2)*K + k2;
                    int mlk1 = (l*(l+1)/2+mm)*N*K + k1*N + n;
                    int mlk2 = (l1*(l1+1)/2+mm1)*N*K + k2*N + n;
                    int mlk3 = (l2*(l2+1)/2+mm2)*N*K + k2*N + n;

                    T a1, b1, c1, a2, b2, c2, a3, b3, c3;
                    T a1x, b1x, a2x, b2x, a3x, b3x;
                    T a1y, b1y, a2y, b2y, a3y, b3y;
                    T a1z, b1z, a2z, b2z, a3z, b3z;
                    c1 = ( mm % 2 == 0) ?  1.0 : -1.0;
                    c2 = ( mm1 % 2 == 0) ?  1.0 : -1.0;
                    c3 = ( mm2 % 2 == 0) ?  1.0 : -1.0;

                    if (m>=0) {
                        a1 = Sr[i1];
                        b1 = Si[i1];
                        a1x = Srx[mlk1];
                        a1y = Sry[mlk1];
                        a1z = Srz[mlk1];
                        b1x = Six[mlk1];
                        b1y = Siy[mlk1];
                        b1z = Siz[mlk1];
                    }
                    else {
                        a1 = -c1*Sr[i1];
                        b1 =  c1*Si[i1];
                        a1x = -c1*Srx[mlk1];
                        a1y = -c1*Sry[mlk1];
                        a1z = -c1*Srz[mlk1];
                        b1x = c1*Six[mlk1];
                        b1y = c1*Siy[mlk1];
                        b1z = c1*Siz[mlk1];
                    }
                    if (m1>=0) {
                        a2 = Sr[i2];
                        b2 = Si[i2];
                        a2x = Srx[mlk2];
                        a2y = Sry[mlk2];
                        a2z = Srz[mlk2];
                        b2x = Six[mlk2];
                        b2y = Siy[mlk2];
                        b2z = Siz[mlk2];
                    }
                    else {
                        a2 = -c2*Sr[i2];
                        b2 =  c2*Si[i2];
                        a2x = -c2*Srx[mlk2];
                        a2y = -c2*Sry[mlk2];
                        a2z = -c2*Srz[mlk2];
                        b2x = c2*Six[mlk2];
                        b2y = c2*Siy[mlk2];
                        b2z = c2*Siz[mlk2];
                    }
                    if (m2>=0) {
                        a3 = Sr[i3];
                        b3 = Si[i3];
                        a3x = Srx[mlk3];
                        a3y = Sry[mlk3];
                        a3z = Srz[mlk3];
                        b3x = Six[mlk3];
                        b3y = Siy[mlk3];
                        b3z = Siz[mlk3];
                    }
                    else {
                        a3 = -c3*Sr[i3];
                        b3 =  c3*Si[i3];
                        a3x = -c3*Srx[mlk3];
                        a3y = -c3*Sry[mlk3];
                        a3z = -c3*Srz[mlk3];
                        b3x = c3*Six[mlk3];
                        b3y = c3*Siy[mlk3];
                        b3z = c3*Siz[mlk3];
                    }
                    T c = cg[rowm[i]+j];

                    //tmp = tmp + cg[rowm[i]+j]*(a1*a2*a3 + a2*b1*b3 + a3*b1*b2 - a1*b2*b3);
                    T t1 = a1x*a2*a3 + a1*a2x*a3 + a1*a2*a3x;
                    T t2 = a2x*b1*b3 + a2*b1x*b3 + a2*b1*b3x;
                    T t3 = a3x*b1*b2 + a3*b1x*b2 + a3*b1*b2x;
                    T t4 = a1x*b2*b3 + a1*b2x*b3 + a1*b2*b3x;
                    bx[ii] = bx[ii] + c*(t1 + t2 + t3 - t4);

                    t1 = a1y*a2*a3 + a1*a2y*a3 + a1*a2*a3y;
                    t2 = a2y*b1*b3 + a2*b1y*b3 + a2*b1*b3y;
                    t3 = a3y*b1*b2 + a3*b1y*b2 + a3*b1*b2y;
                    t4 = a1y*b2*b3 + a1*b2y*b3 + a1*b2*b3y;
                    by[ii] = by[ii] + c*(t1 + t2 + t3 - t4);

                    t1 = a1z*a2*a3 + a1*a2z*a3 + a1*a2*a3z;
                    t2 = a2z*b1*b3 + a2*b1z*b3 + a2*b1*b3z;
                    t3 = a3z*b1*b2 + a3*b1z*b2 + a3*b1*b2z;
                    t4 = a1z*b2*b3 + a1*b2z*b3 + a1*b2*b3z;
                    bz[ii] = bz[ii] + c*(t1 + t2 + t3 - t4);
                }
            }
        }
}
template void cpuRadialSphericalHarmonicsBispectrumDeriv(double*, double*, double*, double*, double*, double*,
        double*, double*, double*, double*, double*, double*, int*, int*, int*, int*, int, int, int, int);
template void cpuRadialSphericalHarmonicsBispectrumDeriv(float*, float*, float*, float*, float*, float*,
        float*, float*, float*, float*, float*, float*, int*, int*, int*, int*, int, int, int, int);

template<typename...T>
void __enzyme_autodiff(void*, T...);
extern "C" {
extern int enzyme_const;
}

#define N 10
int main() {
    double Ylmr[3*N];
    double Ylmi[3*N];

    double the[3*N] = { 0.0 };
    double phi[3*N] = { 0.0 };
    double P[3*N];
    double tmp[3*N];
    double fac[3*N];
    double C[3*N];
    double pi = M_PI;
    int L = 1;
    cpuSphericalHarmonics(Ylmr, Ylmi, the, phi, P, tmp, fac, C, pi, L, N);
    printf("Ylmr[0]=%f\n", Ylmr[0]);

    double d_Ylmr[3*N] = { 0.0 };
    double d_Ylmi[3*N] = { 0.0 };
    double d_the[3*N] = { 1.0 };
    double d_phi[3*N] = { 0.0 };

    __enzyme_autodiff((void*)cpuSphericalHarmonics<double>, Ylmr, d_Ylmr, Ylmi, d_Ylmi, the, d_the, phi, d_phi,
                                            enzyme_const, P,
                                            enzyme_const, tmp,
                                            enzyme_const, fac,
                                            enzyme_const, C,
                                            enzyme_const, pi,
                                            enzyme_const, L,
                                            enzyme_const, N);
    printf("d_Ylmr[0]=%f\n", d_Ylmr[0]);
}
