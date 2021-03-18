#ifndef __CPUSPHERICALBESSEL
#define __CPUSPHERICALBESSEL

// https://en.wikipedia.org/wiki/Bessel_function    
template <typename T> void cpuSphericalBessel(T *g, T *x, T *y, T *z, T *x0, T *f, int L, int K, int N)
{                     
    // L             : the maximum degree of spherical harmonics
    // K             : the maximum degree of spherical Bessel functions
    // N             : length of x, y, z
    // x0 [(L+1)*K]  : zeros of pherical Bessel functions
    // f  [L+1]      : temporary storage for the recurrence formula
    // g  [N*K*(L+1)]: spherical Bessel functions
    
    // Compute spherical Bessel functions, g_{lk}(x,y,z), for l = 0,1,...,L and k = 1,2,...,K
    //  l = 0:    0    1     2    ....   K-1
    //  l = 1:    K   K+1   K+2   ....  2K-1
    //  l = 2:   2K  2K+1  2K+2   ....  3K-1        
    //  l = 3:   3K  3K+1  3K+2   ....  4K-1        
    //  l = 4:   4K  4K+1  4K+2   ....  5K-1        
    //  ....
    //  l = L:   LK  LK+1  LK+2   .... (L+1)K-1        
    // The total number of spherical Bessel functions is K*(L+1).
    
    //#pragma omp parallel for
    for (int i=0; i<N; i++) { // loop over each neighbor atom                               
        int l, k;        
        T r = sqrt(x[i]*x[i] + y[i]*y[i] + z[i]*z[i]);
        T xr;        
        for (k=0; k<K; k++) {
            // Compute spherical Bessel functions, g_{lk}(x,y,z), for l = 0
            l = 0;
            xr = x0[k*(L+1)+l]*r;                        
            g[l*N*K + k*N + i] = cos(xr)/xr;
            
            // Compute spherical Bessel functions, g_{lk}(x,y,z), for l = 1
            l = 1;
            xr = x0[k*(L+1)+l]*r;                        
            g[l*N*K + k*N + i] = cos(xr)/(xr*xr) + sin(xr)/xr;
            
            for (l=2; l<(L+1); l++) {
                // Compute spherical Bessel functions, g_{lk}(x,y,z), for l > 1 using recurrence formula
                xr = x0[k*(L+1)+l]*r;                     
                f[0] = -cos(xr)/xr;
                f[1] = -cos(xr)/(xr*xr) - sin(xr)/xr;
                for (int j=2; j<(l+1); j++) 
                    f[j] = ((2*(j+1)-3)/xr)*f[(j-1)] - f[j-2];                 
                g[l*N*K + k*N + i] = -f[l];
            }            
        }        
    }
}
template void cpuSphericalBessel(double*, double*, double*, double*, double*, double*, int, int, int);
template void cpuSphericalBessel(float*, float*, float*, float*, float*, float*, int, int, int);

template <typename T> void cpuSphericalBessel(T *g, T *xij, T *x0, T *f, int L, int K, int N)
{                         
    //#pragma omp parallel for
    for (int i=0; i<N; i++) { // loop over each neighbor atom                               
        int l, k;        
        T x = xij[i*3];
        T y = xij[i*3+1];
        T z = xij[i*3+2];                
        T r = sqrt(x*x + y*y + z*z);
        T xr;        
        for (k=0; k<K; k++) {
            // Compute spherical Bessel functions, g_{lk}(x,y,z), for l = 0
            l = 0;
            xr = x0[k*(L+1)+l]*r;                        
            g[l*N*K + k*N + i] = cos(xr)/xr;
            
            // Compute spherical Bessel functions, g_{lk}(x,y,z), for l = 1
            l = 1;
            xr = x0[k*(L+1)+l]*r;                        
            g[l*N*K + k*N + i] = cos(xr)/(xr*xr) + sin(xr)/xr;
            
            for (l=2; l<(L+1); l++) {
                // Compute spherical Bessel functions, g_{lk}(x,y,z), for l > 1 using recurrence formula
                xr = x0[k*(L+1)+l]*r;                     
                f[0] = -cos(xr)/xr;
                f[1] = -cos(xr)/(xr*xr) - sin(xr)/xr;
                for (int j=2; j<(l+1); j++) 
                    f[j] = ((2*(j+1)-3)/xr)*f[(j-1)] - f[j-2];                 
                g[l*N*K + k*N + i] = -f[l];
            }            
        }        
    }
}
template void cpuSphericalBessel(double*, double*, double*, double*, int, int, int);
template void cpuSphericalBessel(float*, float*,  float*, float*, int, int, int);

template <typename T> void cpuSphericalBesselWithDeriv(T *gr, T *grx, T *gry, T *grz, 
        T *x, T *y, T *z, T *x0, T *f, T *df, int L, int K, int N)
{                        
    // L                      : the maximum degree of spherical harmonics
    // K                      : the maximum degree of spherical Bessel functions
    // N                      : length of x, y, z
    // x0  [(L+1)*K]          : zeros of pherical Bessel functions
    // f   [L+1]              : temporary storage for the recurrence formula
    // df   [L+1]             : temporary storage for the recurrence formula    
    // grx  [N*K*(L+1)] : x-derivative of real spherical harmonics Bessel functions
    // gry  [N*K*(L+1)] : y-derivative of real spherical harmonics Bessel functions
    // grz  [N*K*(L+1)] : z-derivative of real spherical harmonics Bessel functions
        
    // https://en.wikipedia.org/wiki/Bessel_function
    // Spherical Bessel functions, g_{lk}(x,y,z), for l = 0,1,...,L and k = 1,2,...,K
    //  l = 0:    0    1     2    ....   K-1
    //  l = 1:    K   K+1   K+2   ....  2K-1
    //  l = 2:   2K  2K+1  2K+2   ....  3K-1        
    //  l = 3:   3K  3K+1  3K+2   ....  4K-1        
    //  l = 4:   4K  4K+1  4K+2   ....  5K-1        
    //  ....
    //  l = L:   LK  LK+1  LK+2   .... (L+1)K-1        
    // The total number of spherical Bessel functions is K*(L+1).
                
    //#pragma omp parallel for
    for (int i=0; i<N; i++) { // loop over each neighbor atom                               
        int k, j;
        
        T r = sqrt(x[i]*x[i] + y[i]*y[i] + z[i]*z[i]);             
        T Rx = x[i]/r;
        T Ry = y[i]/r;
        T Rz = z[i]/r;                
     
        int l = 0;               
        T xr, g, gR, gx, gy, gz;
        for (k=0; k<K; k++) {
            // https://en.wikipedia.org/wiki/Bessel_function#Spherical_Bessel_functions:_jn,_yn             
            // Spherical Bessel g_{lk}(x,y,z) and their derivatives for l = 0             
            xr = x0[k*(L+1)+l]*r;                        
            g = cos(xr)/xr;
            gR = x0[k*(L+1)+l]*(-cos(xr)/(xr*xr) - sin(xr)/xr);
            gx = gR*Rx;
            gy = gR*Ry;
            gz = gR*Rz;            
            
            // derivatives of g_{lk}(x,y,z) for l = 0
            j = k*N + i;                
            gr[j] = g;  
            grx[j] = gx;
            gry[j] = gy;
            grz[j] = gz;            
        }
                        
        l = 1;
        for (k=0; k<K; k++) {
            // Spherical Bessel g_{lk}(x,y,z) and their derivatives for l = 1
            xr = x0[k*(L+1)+l]*r;                        
            g = cos(xr)/(xr*xr) + sin(xr)/xr;
            gR = x0[k*(L+1)+l]*(cos(xr)/xr - (2*cos(xr))/(xr*xr*xr) - (2*sin(xr))/(xr*xr));
            gx = gR*Rx;
            gy = gR*Ry;
            gz = gR*Rz;            
                        
            // derivatives of g_{lk}(x,y,z) for l = 1 
            j = i + N*k + N*K*l;      
            gr[j] = g; 
            grx[j] = gx;
            gry[j] = gy;
            grz[j] = gz;
        }        
                        
        for (l=2; l<=L; l++) {        
            for (k=0; k<K; k++) {
                // Spherical Bessel g_{lk}(x,y,z) and their derivatives for l >=2
                xr = x0[k*(L+1)+l]*r;                                            
                f[0] = -cos(xr)/xr;
                f[1] = -cos(xr)/(xr*xr) - sin(xr)/xr;
                df[0] = x0[k*(L+1)+l]*(cos(xr)/(xr*xr) + sin(xr)/xr);
                df[1] = x0[k*(L+1)+l]*((2*cos(xr))/(xr*xr*xr) - cos(xr)/xr  + (2*sin(xr))/(xr*xr));
                for (int jk=2; jk<(l+1); jk++) {
                    f[jk] = ((2*(jk+1)-3)/xr)*f[jk-1] - f[jk-2];                 
                    df[jk] = ((2*(jk+1)-3)/xr)*df[jk-1] - x0[k*(L+1)+l]*((2*(jk+1)-3)/(xr*xr))*f[jk-1] - df[jk-2];        
                }
                g = -f[l];
                gR = -df[l];                    
                gx = gR*Rx;
                gy = gR*Ry;
                gz = gR*Rz; 

                // derivatives of g_{lk}(x,y,z) for l >= 2          
                j = i + N*k + N*K*l;      
                gr[j] = g; 
                grx[j] = gx;
                gry[j] = gy;
                grz[j] = gz;
            }                       
        }
    }                        
}
template void cpuSphericalBesselWithDeriv(double*, double*, double*, double*, double*, double*, double*, 
        double*, double*, double*, int, int, int);
template void cpuSphericalBesselWithDeriv(float*, float*, float*, float*, float*, float*, float*, 
        float*, float*, float*, int, int, int);

template <typename T> void cpuSphericalBesselWithDeriv(T *gr, T *grx, T *gry, T *grz, 
        T *xij, T *x0, T *f, T *df, int L, int K, int N)
{                                        
    //#pragma omp parallel for
    for (int i=0; i<N; i++) { // loop over each neighbor atom                               
        int k, j;
        
        T x = xij[i*3];
        T y = xij[i*3+1];
        T z = xij[i*3+2];                
        T r = sqrt(x*x + y*y + z*z);                     
        T Rx = x/r;
        T Ry = y/r;
        T Rz = z/r;                
     
        int l = 0;               
        T xr, g, gR, gx, gy, gz;
        for (k=0; k<K; k++) {
            // https://en.wikipedia.org/wiki/Bessel_function#Spherical_Bessel_functions:_jn,_yn             
            // Spherical Bessel g_{lk}(x,y,z) and their derivatives for l = 0             
            xr = x0[k*(L+1)+l]*r;                        
            g = cos(xr)/xr;
            gR = x0[k*(L+1)+l]*(-cos(xr)/(xr*xr) - sin(xr)/xr);
            gx = gR*Rx;
            gy = gR*Ry;
            gz = gR*Rz;            
            
            // derivatives of g_{lk}(x,y,z) for l = 0
            j = k*N + i;                
            gr[j] = g;  
            grx[j] = gx;
            gry[j] = gy;
            grz[j] = gz;            
        }
                        
        l = 1;
        for (k=0; k<K; k++) {
            // Spherical Bessel g_{lk}(x,y,z) and their derivatives for l = 1
            xr = x0[k*(L+1)+l]*r;                        
            g = cos(xr)/(xr*xr) + sin(xr)/xr;
            gR = x0[k*(L+1)+l]*(cos(xr)/xr - (2*cos(xr))/(xr*xr*xr) - (2*sin(xr))/(xr*xr));
            gx = gR*Rx;
            gy = gR*Ry;
            gz = gR*Rz;            
                        
            // derivatives of g_{lk}(x,y,z) for l = 1 
            j = i + N*k + N*K*l;      
            gr[j] = g; 
            grx[j] = gx;
            gry[j] = gy;
            grz[j] = gz;
        }        
                        
        for (l=2; l<=L; l++) {        
            for (k=0; k<K; k++) {
                // Spherical Bessel g_{lk}(x,y,z) and their derivatives for l >=2
                xr = x0[k*(L+1)+l]*r;                                            
                f[0] = -cos(xr)/xr;
                f[1] = -cos(xr)/(xr*xr) - sin(xr)/xr;
                df[0] = x0[k*(L+1)+l]*(cos(xr)/(xr*xr) + sin(xr)/xr);
                df[1] = x0[k*(L+1)+l]*((2*cos(xr))/(xr*xr*xr) - cos(xr)/xr  + (2*sin(xr))/(xr*xr));
                for (int jk=2; jk<(l+1); jk++) {
                    f[jk] = ((2*(jk+1)-3)/xr)*f[jk-1] - f[jk-2];                 
                    df[jk] = ((2*(jk+1)-3)/xr)*df[jk-1] - x0[k*(L+1)+l]*((2*(jk+1)-3)/(xr*xr))*f[jk-1] - df[jk-2];        
                }
                g = -f[l];
                gR = -df[l];                    
                gx = gR*Rx;
                gy = gR*Ry;
                gz = gR*Rz; 

                // derivatives of g_{lk}(x,y,z) for l >= 2          
                j = i + N*k + N*K*l;      
                gr[j] = g; 
                grx[j] = gx;
                gry[j] = gy;
                grz[j] = gz;
            }                       
        }
    }                        
}
template void cpuSphericalBesselWithDeriv(double*, double*, double*, double*, double*, 
        double*, double*, double*, int, int, int);
template void cpuSphericalBesselWithDeriv(float*, float*, float*, float*, float*, 
        float*, float*, float*, int, int, int);

#endif


