#ifndef __GPUSPHERICALHARMONICSBESSEL
#define __GPUSPHERICALHARMONICSBESSEL

#include <stdio.h>

#define Lmax 10

// core function
template <typename T> __global__ void gpuKernelSphericalHarmonicsBessel(T *Sr, T *Si, T *x, T *y, T *z, 
                T *x0, T *fac, T pi, int L, int K, int N)
{                        
    // L                     : the maximum degree of spherical harmonics
    // K                     : the maximum degree of spherical Bessel functions
    // N                     : length of x, y, z
    // x0  [(L+1)*K]         : zeros of pherical Bessel functions
    // P   [(L+1)*(L+2)/2]   : temporary storage for the recurrence formula
    // tmp [(L+1)*(L+2)/2]   : temporary storage for the recurrence formula
    // f   [L+1]             : temporary storage for the recurrence formula
    // fac                   : factorial look-up table
    // Sr  [N*K*(L+1)*(L+1)] : real part of spherical harmonics Bessel functions
    // Si  [N*K*(L+1)*(L+1)] : imag part of spherical harmonics Bessel functions
        
    // Compute spherical harmonics Bessel functions: S_{klm}(x,y,z) = g_{lk}(x,y,z)*Y_{lm}(x,y,z)
    
    // https://en.wikipedia.org/wiki/Spherical_harmonics
    // Spherical harmonics functions, Y_{lm}(x,y,z), for l = 0,1,...,L and m = -l,...,0,...,l        
    //  l = 0:                0
    //  l = 1:              1 2 3
    //  l = 2:            4 5 6 7 8
    //  l = 3:        9 10 11 12 13 14 15
    //  l = 4:    16 17 18 19 20 21 22 23 24
    //  ....
    //  l = L:  L*L L*L+1 ....   ....  (L+1)*(L+1)-1 
    // The total number of spherical harmonics is (L+1)*(L+1).
    
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
        
    // Hence, the total number of spherical harmonics Bessel functions is K*(L+1)*(L+1).            
            
    //#pragma omp parallel for    
    int i = threadIdx.x + blockIdx.x * blockDim.x;     
    while (i < N) { // loop over each neighbor atom               
        int k, j, jm;
        T Pa[Lmax];
        T tmpa[Lmax];
        T fa[Lmax]; 
        
        tmpa[0] = 1.0;
        
        // Cartersian -> Spherical
        T r = sqrt(x[i]*x[i] + y[i]*y[i] + z[i]*z[i]);
        T the = acos(z[i]/r);
        T phi = atan2(y[i],x[i]);
        
        // Spherical harmonics for l = 0        
        int l = 0;
        T Ylmr = sqrt(1/(4*pi));
        T Ylmi = 0.0;                
        T xr, g;
        // Spherical Bessel for l = 0, see https://en.wikipedia.org/wiki/Bessel_function#Spherical_Bessel_functions:_jn,_yn 
        for (k=0; k<K; k++) {
            xr = x0[k*(L+1)+l]*r;                        
            g = cos(xr)/xr;            
            j = k*N + i;                      
            Sr[j] = g*Ylmr;  // real part                   
            Si[j] = g*Ylmi;  // imag part       
        }
                
        // Spherical harmonics for l = 1;
        l = 1;
        T costhe = cos(the);
        T a = -sin(the);        
        int m = 0;    
        Pa[0] = costhe;
        Pa[1] = a;
        
        T C = sqrt((2*l+1)*fac[l-m]/(4*pi*fac[l+m]));
        Ylmr = C*Pa[0];
        Ylmi = 0.0;
        // Spherical Bessel for l = 1, see https://en.wikipedia.org/wiki/Bessel_function#Spherical_Bessel_functions:_jn,_yn
        for (k=0; k<K; k++) {
            xr = x0[k*(L+1)+l]*r;                        
            g = cos(xr)/(xr*xr) + sin(xr)/xr;
            j = ((l*l + l + m)*K + k)*N + i;                
            Sr[j] = g*Ylmr; // real part           
            Si[j] = g*Ylmi; // imag part                                    
        }        
        
        m = 1;
        C = sqrt((2*l+1)*fac[l-m]/(4*pi*fac[l+m]));
        Ylmr = C*cos(phi)*Pa[1];
        Ylmi = C*sin(phi)*Pa[1];                
        // Spherical Bessel for l = 1, see https://en.wikipedia.org/wiki/Bessel_function#Spherical_Bessel_functions:_jn,_yn
        for (k=0; k<K; k++) {
            xr = x0[k*(L+1)+l]*r;                        
            g = cos(xr)/(xr*xr) + sin(xr)/xr; 
            j = ((l*l + l + m)*K + k)*N + i; // spherical harmonics Bessel functions for m > 0                
            Sr[j] = g*Ylmr; // real part           
            Si[j] = g*Ylmi; // imag part                                                            
            jm = ((l*l + l - m)*K + k)*N + i; // spherical harmonics Bessel functions for m < 0                
            Sr[jm] = -Sr[j]; // real part                      
            Si[jm] =  Si[j]; // imag part                                                                       
        }        
                
        for (l=2; l<=L; l++) {                                        
            // Compute associated Legendre polynomial using recurrence formula
            // https://en.wikipedia.org/wiki/Associated_Legendre_polynomials            
            T Pll = Pa[l-1];
            tmpa[(l-1)] = Pll;
            Pa[(l-1)] = (2*(l-1)+1)*(costhe*Pll);
            Pa[l] = (2*(l-1)+1)*(a*Pll);             
            for (m=0; m<l-1; m++) {
                Pll = Pa[m];
                Pa[m] = ((2*(l-1)+1)*(costhe*Pll) - (l+m-1)*tmpa[m])/(T (l-m));
                tmpa[m] = Pll;                
            }
//             for (k=0; k<K; k++) {
//                 // Compute the spherical Bessel functions using recurrence formula
//                 xr = x0[k*(L+1)+l]*r;                                            
//                 fa[0] = -cos(xr)/xr;
//                 fa[1] = -cos(xr)/(xr*xr) - sin(xr)/xr;                    
//                 for (int jk=2; jk<(l+1); jk++) 
//                     fa[jk] = ((2*(jk+1)-3)/xr)*fa[jk-1] - fa[jk-2];   
//                 tmpa[k] = -fa[l];
//             }

            // Compute spherical harmonics Bessel function at level l for m = 0, 1, ..., l
            for (m=0; m<=l; m++) {
                C = sqrt((2*l+1)*fac[l-m]/(4*pi*fac[l+m])); 
                Ylmr = C*(cos(m*phi)*Pa[m]);
                Ylmi = C*(sin(m*phi)*Pa[m]);                
                for (k=0; k<K; k++) {
                    // Compute the spherical Bessel functions using recurrence formula
//                     xr = x0[k*(L+1)+l]*r;                                            
//                     f[0] = -cos(xr)/xr;
//                     f[1] = -cos(xr)/(xr*xr) - sin(xr)/xr;                    
//                     for (int jk=2; jk<(l+1); jk++) 
//                         f[jk] = ((2*(jk+1)-3)/xr)*f[jk-1] - f[jk-2];                                     
//                     g = -f[l];              
                    xr = x0[k*(L+1)+l]*r;                                            
                    fa[0] = -cos(xr)/xr;
                    fa[1] = -cos(xr)/(xr*xr) - sin(xr)/xr;                    
                    for (int jk=2; jk<(l+1); jk++) 
                        fa[jk] = ((2*(jk+1)-3)/xr)*fa[jk-1] - fa[jk-2];   
                    g = -fa[l];              
                    //g = tmpa[k];
                    j = ((l*l + l + m)*K + k)*N + i;   
                    Sr[j] = g*Ylmr; // real part            
                    Si[j] = g*Ylmi; // imag part                                                            
                }                
            }        
            
            // Compute the spherical harmonics Bessel functions for m < 0 using symmetry properties
            C = -1.0;
            for (int m=1; m<=l; m++)  {                 
                for (k=0; k<K; k++) {
                    j =  ((l*l + l + m)*K + k)*N + i;  // spherical harmonics Bessel functions for m > 0 
                    jm = ((l*l + l - m)*K + k)*N + i;  // spherical harmonics Bessel functions for m < 0                 
                    Sr[jm] = C*Sr[j]; // real part
                    Si[jm] =-C*Si[j]; // imag part                                                                             
                }                     
                C = C*(-1.0);        
            }
        }
        i += blockDim.x * gridDim.x;
    }          
}
template <typename T> void gpuSphericalHarmonicsBessel(T *Sr, T *Si, T *x, T *y, T *z, 
                T *x0, T *P, T *tmp, T *f, T *fac, T pi, int L, int K, int N)
{        
    int blockDim = 256;
    int gridDim = (N + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;
    gpuKernelSphericalHarmonicsBessel<<<gridDim, blockDim>>>(Sr, Si, x, y, z, 
            x0, fac, pi, L, K, N);
}
template void gpuSphericalHarmonicsBessel(double*, double*, double*, double*, double*, double*, 
        double*, double*, double*, double*, double, int, int, int);
template void gpuSphericalHarmonicsBessel(float*, float*, float*, float*, float*, float*, 
        float*, float*, float*, float*, float, int, int, int);

// core function
template <typename T> __global__ void gpuKernelSphericalHarmonicsBessel(T *Sr, T *Si, T *xij, 
                T *x0, T *fac, T pi, int L, int K, int N)
{                        
    // L                     : the maximum degree of spherical harmonics
    // K                     : the maximum degree of spherical Bessel functions
    // N                     : length of x, y, z
    // x0  [(L+1)*K]         : zeros of pherical Bessel functions
    // P   [(L+1)*(L+2)/2]   : temporary storage for the recurrence formula
    // tmp [(L+1)*(L+2)/2]   : temporary storage for the recurrence formula
    // f   [L+1]             : temporary storage for the recurrence formula
    // fac                   : factorial look-up table
    // Sr  [N*K*(L+1)*(L+1)] : real part of spherical harmonics Bessel functions
    // Si  [N*K*(L+1)*(L+1)] : imag part of spherical harmonics Bessel functions
        
    // Compute spherical harmonics Bessel functions: S_{klm}(x,y,z) = g_{lk}(x,y,z)*Y_{lm}(x,y,z)
    
    // https://en.wikipedia.org/wiki/Spherical_harmonics
    // Spherical harmonics functions, Y_{lm}(x,y,z), for l = 0,1,...,L and m = -l,...,0,...,l        
    //  l = 0:                0
    //  l = 1:              1 2 3
    //  l = 2:            4 5 6 7 8
    //  l = 3:        9 10 11 12 13 14 15
    //  l = 4:    16 17 18 19 20 21 22 23 24
    //  ....
    //  l = L:  L*L L*L+1 ....   ....  (L+1)*(L+1)-1 
    // The total number of spherical harmonics is (L+1)*(L+1).
    
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
        
    // Hence, the total number of spherical harmonics Bessel functions is K*(L+1)*(L+1).            
            
    //#pragma omp parallel for    
    int i = threadIdx.x + blockIdx.x * blockDim.x;     
    while (i < N) { // loop over each neighbor atom               
        int k, j, jm;
        T Pa[Lmax];
        T tmpa[Lmax];
        T fa[Lmax]; 
        T tmpb[Lmax];

        tmpa[0] = 1.0;
        
        T x = xij[i*3];
        T y = xij[i*3+1];
        T z = xij[i*3+2];        
        
        if (sqrt(x*x + y*y) < 2e-12) {
            T signx = (x >= 0.0) ? 1.0 : -1.0;
            T signy = (y >= 0.0) ? -1.0 : 1.0;
            x = fabs(x) > 1e-12 ? x : signx*1e-12;
            y = fabs(y) > 1e-12 ? y : signy*1e-12;
        }
                
        // Cartersian -> Spherical
        T r = sqrt(x*x + y*y + z*z);
        //T the = acos(z/r);
        T phi = atan2(y,x);
        
        // Spherical harmonics for l = 0        
        int l = 0;
        T Ylmr = sqrt(1/(4*pi));
        T Ylmi = 0.0;                
        T xr, g;
        // Spherical Bessel for l = 0, see https://en.wikipedia.org/wiki/Bessel_function#Spherical_Bessel_functions:_jn,_yn 
        for (k=0; k<K; k++) {
            xr = x0[k*(L+1)+l]*r;                        
            g = cos(xr)/xr;            
            j = k*N + i;                      
            Sr[j] = g*Ylmr;  // real part                   
            Si[j] = g*Ylmi;  // imag part       
        }
                
        // Spherical harmonics for l = 1;
        l = 1;
//         T costhe = cos(the);
//         T a = -sin(the);        
        T costhe = z/r; 
        T a = -sqrt((r-z)*(r+z))/r;             
        int m = 0;    
        Pa[0] = costhe;
        Pa[1] = a;
        
        T C = sqrt((2*l+1)*fac[l-m]/(4*pi*fac[l+m]));
        Ylmr = C*Pa[0];
        Ylmi = 0.0;
        // Spherical Bessel for l = 1, see https://en.wikipedia.org/wiki/Bessel_function#Spherical_Bessel_functions:_jn,_yn
        for (k=0; k<K; k++) {
            xr = x0[k*(L+1)+l]*r;                        
            g = cos(xr)/(xr*xr) + sin(xr)/xr;
            j = ((l*l + l + m)*K + k)*N + i;                
            Sr[j] = g*Ylmr; // real part           
            Si[j] = g*Ylmi; // imag part                                    
        }        
        
        m = 1;
        C = sqrt((2*l+1)*fac[l-m]/(4*pi*fac[l+m]));
        Ylmr = C*cos(phi)*Pa[1];
        Ylmi = C*sin(phi)*Pa[1];                
        //cout<<m<<"  "<<C<<"  "<<cos(phi)<<"  "<<P[1]<<endl;
        // Spherical Bessel for l = 1, see https://en.wikipedia.org/wiki/Bessel_function#Spherical_Bessel_functions:_jn,_yn
        for (k=0; k<K; k++) {
            xr = x0[k*(L+1)+l]*r;                        
            g = cos(xr)/(xr*xr) + sin(xr)/xr; 
            j = ((l*l + l + m)*K + k)*N + i; // spherical harmonics Bessel functions for m > 0                
            Sr[j] = g*Ylmr; // real part           
            Si[j] = g*Ylmi; // imag part              
            //cout<<j<<"  "<<g<<"  "<<Ylmr<<"  "<<Sr[j]<<endl;
            jm = ((l*l + l - m)*K + k)*N + i; // spherical harmonics Bessel functions for m < 0                
            Sr[jm] = -Sr[j]; // real part                      
            Si[jm] =  Si[j]; // imag part                                                                       
        }        
                
        for (l=2; l<=L; l++) {                                        
            // Compute associated Legendre polynomial using recurrence formula
            // https://en.wikipedia.org/wiki/Associated_Legendre_polynomials            
            T Pll = Pa[l-1];
            tmpa[(l-1)] = Pll;
            Pa[(l-1)] = (2*(l-1)+1)*(costhe*Pll);
            Pa[l] = (2*(l-1)+1)*(a*Pll);             
            for (m=0; m<l-1; m++) {
                Pll = Pa[m];
                Pa[m] = ((2*(l-1)+1)*(costhe*Pll) - (l+m-1)*tmpa[m])/(T (l-m));
                tmpa[m] = Pll;                
            }
            for (k=0; k<K; k++) {
                // Compute the spherical Bessel functions using recurrence formula
                xr = x0[k*(L+1)+l]*r;                                            
                fa[0] = -cos(xr)/xr;
                fa[1] = -cos(xr)/(xr*xr) - sin(xr)/xr;                    
                for (int jk=2; jk<(l+1); jk++) 
                    fa[jk] = ((2*(jk+1)-3)/xr)*fa[jk-1] - fa[jk-2];   
                tmpb[k] = -fa[l];
            }
                            
            // Compute spherical harmonics Bessel function at level l for m = 0, 1, ..., l
            for (m=0; m<=l; m++) {
                C = sqrt((2*l+1)*fac[l-m]/(4*pi*fac[l+m])); 
                Ylmr = C*(cos(m*phi)*Pa[m]);
                Ylmi = C*(sin(m*phi)*Pa[m]);                
                for (k=0; k<K; k++) {
                    // Compute the spherical Bessel functions using recurrence formula
//                     xr = x0[k*(L+1)+l]*r;                                            
//                     fa[0] = -cos(xr)/xr;
//                     fa[1] = -cos(xr)/(xr*xr) - sin(xr)/xr;                    
//                     for (int jk=2; jk<(l+1); jk++) 
//                         fa[jk] = ((2*(jk+1)-3)/xr)*fa[jk-1] - fa[jk-2];                               
//                     g = -fa[l];     
                    g = tmpb[k];
                    j = ((l*l + l + m)*K + k)*N + i;   
                    Sr[j] = g*Ylmr; // real part            
                    Si[j] = g*Ylmi; // imag part                                                            
                }                
            }        
            
            // Compute the spherical harmonics Bessel functions for m < 0 using symmetry properties
            C = -1.0;
            for (int m=1; m<=l; m++)  {                 
                for (k=0; k<K; k++) {
                    j =  ((l*l + l + m)*K + k)*N + i;  // spherical harmonics Bessel functions for m > 0 
                    jm = ((l*l + l - m)*K + k)*N + i;  // spherical harmonics Bessel functions for m < 0                 
                    Sr[jm] = C*Sr[j]; // real part
                    Si[jm] =-C*Si[j]; // imag part                                                                             
                }                     
                C = C*(-1.0);        
            }
        }
        i += blockDim.x * gridDim.x;
    }          
}
template <typename T> void gpuSphericalHarmonicsBessel(T *Sr, T *Si, T *xij, 
                T *x0, T *P, T *tmp, T *f, T *fac, T pi, int L, int K, int N)
{        
    int blockDim = 256;
    int gridDim = (N + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;
    gpuKernelSphericalHarmonicsBessel<<<gridDim, blockDim>>>(Sr, Si, xij, 
            x0, fac, pi, L, K, N);
}
template void gpuSphericalHarmonicsBessel(double*, double*, double*, double*, 
        double*, double*, double*, double*, double, int, int, int);
template void gpuSphericalHarmonicsBessel(float*, float*, float*, float*, 
        float*, float*, float*, float*, float, int, int, int);

// core function
template <typename T> __global__ void gpuKernelSphericalHarmonicsBesselWithDeriv(T *Sr, T *Si,
        T *Srx, T *Six, T *Sry, T *Siy, T *Srz, T *Siz, T *x, T *y, T *z, 
        T *x0, T *fac, T pi, int L, int K, int N)
{                        
    // L                      : the maximum degree of spherical harmonics
    // K                      : the maximum degree of spherical Bessel functions
    // N                      : length of x, y, z
    // x0  [(L+1)*K]          : zeros of pherical Bessel functions
    // P   [(L+1)*(L+2)/2]    : temporary storage for the recurrence formula
    // tmp [(L+1)*(L+2)/2]    : temporary storage for the recurrence formula
    // f   [L+1]              : temporary storage for the recurrence formula
    // dP   [(L+1)*(L+2)/2]   : temporary storage for the recurrence formula
    // dtmp [(L+1)*(L+2)/2]   : temporary storage for the recurrence formula
    // df   [L+1]             : temporary storage for the recurrence formula    
    // fac                    : factorial look-up table
    // Srx  [N*K*(L+1)*(L+1)] : x-derivative of real spherical harmonics Bessel functions
    // Six  [N*K*(L+1)*(L+1)] : x-derivative of imag spherical harmonics Bessel functions
    // Sry  [N*K*(L+1)*(L+1)] : y-derivative of real spherical harmonics Bessel functions
    // Siy  [N*K*(L+1)*(L+1)] : y-derivative of imag spherical harmonics Bessel functions
    // Srz  [N*K*(L+1)*(L+1)] : z-derivative of real spherical harmonics Bessel functions
    // Siz  [N*K*(L+1)*(L+1)] : z-derivative of imag spherical harmonics Bessel functions
    
    // Compute partial derivatives of spherical harmonics Bessel functions: S_{klm}(x,y,z) = g_{lk}(x,y,z)*Y_{lm}(x,y,z)
    
    // https://en.wikipedia.org/wiki/Spherical_harmonics
    // Spherical harmonics functions, Y_{lm}(x,y,z), for l = 0,1,...,L and m = -l,...,0,...,l        
    //  l = 0:                0
    //  l = 1:              1 2 3
    //  l = 2:            4 5 6 7 8
    //  l = 3:        9 10 11 12 13 14 15
    //  l = 4:    16 17 18 19 20 21 22 23 24
    //  ....
    //  l = L:  L*L L*L+1 ....   ....  (L+1)*(L+1)-1 
    // The total number of spherical harmonics is (L+1)*(L+1).
    
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
            
    // Hence, the total number of spherical harmonics Bessel functions, S_{klm}(x,y,z), is K*(L+1)*(L+1).            
    
    //#pragma omp parallel for
    int i = threadIdx.x + blockIdx.x * blockDim.x;     
    while (i < N) { // loop over each neighbor atom                               
        int k, j, jm;
        T Pa[Lmax];
        T tmpa[Lmax];
        T fa[Lmax]; 
        T dPa[Lmax];
        T dtmpa[Lmax];
        T dfa[Lmax];
        
        // Cartersian -> Spherical
        T r = sqrt(x[i]*x[i] + y[i]*y[i] + z[i]*z[i]);
        //T the = acos(z[i]/r);
        T phi = atan2(y[i],x[i]);
                
//         cout<<z[i]<<endl;
//         cout<<r<<endl;
//         cout<<the<<endl;
//         cout<<sin(the)<<endl;
//         cout<<cos(the)<<endl;
        
        // partial derivatives of spherical coordinates
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
        
        // Spherical harmonics Y_{lm}(x,y,z) and their derivatives for l = 0        
        int l = 0;
        T Ylmr = sqrt(1/(4*pi));
        T Ylmi = 0.0;
        T YlmrThe = 0.0;
        T YlmrPhi = 0.0;
        T YlmiThe = 0.0;        
        T YlmiPhi = 0.0;        
        T Ylmrx = YlmrThe*Thex + YlmrPhi*Phix;
        T Ylmry = YlmrThe*They + YlmrPhi*Phiy;
        T Ylmrz = YlmrThe*Thez + YlmrPhi*Phiz;
        T Ylmix = YlmiThe*Thex + YlmiPhi*Phix;
        T Ylmiy = YlmiThe*They + YlmiPhi*Phiy;
        T Ylmiz = YlmiThe*Thez + YlmiPhi*Phiz;
        
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
            
            // derivatives of S_{klm}(x,y,z) = g_{lk}(x,y,z)*Y_{lm}(x,y,z) for l = 0
            j = k*N + i;                
            Sr[j] = g*Ylmr;  // real part                   
            Si[j] = g*Ylmi;  // imag part                   
            Srx[j] = gx*Ylmr + g*Ylmrx;
            Sry[j] = gy*Ylmr + g*Ylmry;
            Srz[j] = gz*Ylmr + g*Ylmrz;
            Six[j] = gx*Ylmi + g*Ylmix;
            Siy[j] = gy*Ylmi + g*Ylmiy;
            Siz[j] = gz*Ylmi + g*Ylmiz;            
        }
                        
        l = 1;
        T costhe = z[i]/r; //cos(the);
        T a = -sqrt(1-costhe*costhe);        
        //T a = -sin(the);        
        T dcosthe = a;
        T da = -costhe;        
        int m = 0;    
        Pa[0] = costhe;
        Pa[1] = a;
        dPa[0] = dcosthe;    
        dPa[1] = da;
        tmpa[0] = 1.0;
        dtmpa[0] = 0.0;
        
        // Spherical harmonics Y_{lm}(x,y,z) and their derivatives for l = 1 and m = 0  
        m = 0;
        T C = sqrt((2*l+1)*fac[l-m]/(4*pi*fac[l+m]));
        Ylmr = C*Pa[0];
        Ylmi = 0.0;
        YlmrThe = C*a;    
        YlmiThe = 0.0;    
        YlmrPhi = 0.0;          
        YlmiPhi = 0.0;                  
        Ylmrx = YlmrThe*Thex + YlmrPhi*Phix;
        Ylmry = YlmrThe*They + YlmrPhi*Phiy;
        Ylmrz = YlmrThe*Thez + YlmrPhi*Phiz;
        Ylmix = YlmiThe*Thex + YlmiPhi*Phix;
        Ylmiy = YlmiThe*They + YlmiPhi*Phiy;
        Ylmiz = YlmiThe*Thez + YlmiPhi*Phiz;        
        for (k=0; k<K; k++) {
            // Spherical Bessel g_{lk}(x,y,z) and their derivatives for l = 1
            xr = x0[k*(L+1)+l]*r;                        
            g = cos(xr)/(xr*xr) + sin(xr)/xr;
            gR = x0[k*(L+1)+l]*(cos(xr)/xr - (2*cos(xr))/(xr*xr*xr) - (2*sin(xr))/(xr*xr));
            gx = gR*Rx;
            gy = gR*Ry;
            gz = gR*Rz;            
                        
            // derivatives of S_{klm}(x,y,z) = g_{lk}(x,y,z)*Y_{lm}(x,y,z) for l = 1 and m = 0
            j = ((l*l + l + m)*K + k)*N + i;      
            Sr[j] = g*Ylmr; // real part           
            Si[j] = g*Ylmi; // imag part                                                
            Srx[j] = gx*Ylmr + g*Ylmrx;
            Sry[j] = gy*Ylmr + g*Ylmry;
            Srz[j] = gz*Ylmr + g*Ylmrz;
            Six[j] = gx*Ylmi + g*Ylmix;
            Siy[j] = gy*Ylmi + g*Ylmiy;
            Siz[j] = gz*Ylmi + g*Ylmiz;            
        }        
        
        // Spherical harmonics Y_{lm}(x,y,z) and their derivatives for l = 1 and m = 1
        m = 1;
        C = sqrt((2*l+1)*fac[l-m]/(4*pi*fac[l+m]));
        Ylmr = C*cos(phi)*Pa[1];
        Ylmi = C*sin(phi)*Pa[1];                
        YlmrThe = C*(cos(phi)*da);      
        YlmiThe = C*(sin(phi)*da);      
        YlmrPhi = -C*(sin(phi)*Pa[1]); 
        YlmiPhi = C*(cos(phi)*Pa[1]);  
        Ylmrx = YlmrThe*Thex + YlmrPhi*Phix;
        Ylmry = YlmrThe*They + YlmrPhi*Phiy;
        Ylmrz = YlmrThe*Thez + YlmrPhi*Phiz;
        Ylmix = YlmiThe*Thex + YlmiPhi*Phix;
        Ylmiy = YlmiThe*They + YlmiPhi*Phiy;
        Ylmiz = YlmiThe*Thez + YlmiPhi*Phiz;
        //cout<<m<<"  "<<C<<"  "<<cos(phi)<<"  "<<Pa[1]<<"  "<<the<<"  "<<sin(the)<<endl;
        for (k=0; k<K; k++) {
            // Spherical Bessel g_{lk}(x,y,z) and their derivatives for l = 1
            xr = x0[k*(L+1)+l]*r;                        
            g = cos(xr)/(xr*xr) + sin(xr)/xr;
            gR = x0[k*(L+1)+l]*(cos(xr)/xr - (2*cos(xr))/(xr*xr*xr) - (2*sin(xr))/(xr*xr));
            gx = gR*Rx;
            gy = gR*Ry;
            gz = gR*Rz;            
            
            // derivatives of S_{klm}(x,y,z) = g_{lk}(x,y,z)*Y_{lm}(x,y,z) for l = 1 and m = 1
            j = ((l*l + l + m)*K + k)*N + i;          
            Sr[j] = g*Ylmr; // real part           
            Si[j] = g*Ylmi; // imag part        
            //cout<<j<<"  "<<g<<"  "<<Ylmr<<"  "<<Sr[j]<<endl;
            Srx[j] = gx*Ylmr + g*Ylmrx;
            Sry[j] = gy*Ylmr + g*Ylmry;
            Srz[j] = gz*Ylmr + g*Ylmrz;
            Six[j] = gx*Ylmi + g*Ylmix;
            Siy[j] = gy*Ylmi + g*Ylmiy;
            Siz[j] = gz*Ylmi + g*Ylmiz;        
            
            // derivatives of S_{klm}(x,y,z) = g_{lk}(x,y,z)*Y_{lm}(x,y,z) for l = 1 and m = -1
            jm = ((l*l + l - m)*K + k)*N + i;       
            Sr[jm] = -Sr[j]; // real part                      
            Si[jm] =  Si[j]; // imag part                                                                                   
            Srx[jm] = -Srx[j];
            Sry[jm] = -Sry[j];
            Srz[jm] = -Srz[j];
            Six[jm] = Six[j];
            Siy[jm] = Siy[j];
            Siz[jm] = Siz[j];            
        }        
                
        for (l=2; l<=L; l++) {        
            // Compute associated Legendre polynomials P_m and their derivatives using recurrence formula
            // https://en.wikipedia.org/wiki/Associated_Legendre_polynomials                        
            T Pll = Pa[l-1];
            tmpa[(l-1)] = Pll;
            Pa[(l-1)] = (2*(l-1)+1)*(costhe*Pll);
            Pa[l] = (2*(l-1)+1)*(a*Pll); 
            T dPll = dPa[l-1];
            dtmpa[l-1] = dPll;
            dPa[l-1] = (2*(l-1)+1)*(costhe*dPll + dcosthe*Pll);   
            dPa[l] = (2*(l-1)+1)*(a*dPll + da*Pll);             
            for (m=0; m<l-1; m++) {
                Pll = Pa[m];
                Pa[m] = ((2*(l-1)+1)*(costhe*Pll) - (l+m-1)*tmpa[m])/(T (l-m));
                tmpa[m] = Pll;
                dPll = dPa[m];
                dPa[m] = ((2*(l-1)+1)*(costhe*dPll+dcosthe*Pll) - (l+m-1)*dtmpa[m])/(T (l-m));
                dtmpa[m] = dPll;
            }
             
//             for (k=0; k<K; k++) {
//                 // Spherical Bessel g_{lk}(x,y,z) and their derivatives for l >=2
//                 xr = x0[k*(L+1)+l]*r;                                            
//                 fa[0] = -cos(xr)/xr;
//                 fa[1] = -cos(xr)/(xr*xr) - sin(xr)/xr;
//                 dfa[0] = x0[k*(L+1)+l]*(cos(xr)/(xr*xr) + sin(xr)/xr);
//                 dfa[1] = x0[k*(L+1)+l]*((2*cos(xr))/(xr*xr*xr) - cos(xr)/xr  + (2*sin(xr))/(xr*xr));
//                 for (int jk=2; jk<(l+1); jk++) {
//                     fa[jk] = ((2*(jk+1)-3)/xr)*fa[jk-1] - fa[jk-2];                 
//                     dfa[jk] = ((2*(jk+1)-3)/xr)*dfa[jk-1] - x0[k*(L+1)+l]*((2*(jk+1)-3)/(xr*xr))*fa[jk-1] - dfa[jk-2];        
//                 }
//                 tmpa[k] = -fa[l];
//                 dtmpa[k] = -dfa[l];                                    
//             }

            // Compute spherical harmonics Bessel functions and their derivatives at level l for m = 0, 1, ..., l
            for (m=0; m<=l; m++) {
                // Spherical harmonics Y_{lm}(x,y,z) and their derivatives for l and m = 0,1,..,l  
                C = sqrt((2*l+1)*fac[l-m]/(4*pi*fac[l+m])); 
                Ylmr = C*(cos(m*phi)*Pa[m]);
                Ylmi = C*(sin(m*phi)*Pa[m]);
                YlmrThe = C*(cos(m*phi)*dPa[m]); 
                YlmiThe = C*(sin(m*phi)*dPa[m]); 
                YlmrPhi = -(m*C)*(sin(m*phi)*Pa[m]); 
                YlmiPhi = (m*C)*(cos(m*phi)*Pa[m]);          
                Ylmrx = YlmrThe*Thex + YlmrPhi*Phix;
                Ylmry = YlmrThe*They + YlmrPhi*Phiy;
                Ylmrz = YlmrThe*Thez + YlmrPhi*Phiz;
                Ylmix = YlmiThe*Thex + YlmiPhi*Phix;
                Ylmiy = YlmiThe*They + YlmiPhi*Phiy;
                Ylmiz = YlmiThe*Thez + YlmiPhi*Phiz; 

                for (k=0; k<K; k++) {    
                    xr = x0[k*(L+1)+l]*r;                                            
                    fa[0] = -cos(xr)/xr;
                    fa[1] = -cos(xr)/(xr*xr) - sin(xr)/xr;
                    dfa[0] = x0[k*(L+1)+l]*(cos(xr)/(xr*xr) + sin(xr)/xr);
                    dfa[1] = x0[k*(L+1)+l]*((2*cos(xr))/(xr*xr*xr) - cos(xr)/xr  + (2*sin(xr))/(xr*xr));
                    for (int jk=2; jk<(l+1); jk++) {
                        fa[jk] = ((2*(jk+1)-3)/xr)*fa[jk-1] - fa[jk-2];                 
                        dfa[jk] = ((2*(jk+1)-3)/xr)*dfa[jk-1] - x0[k*(L+1)+l]*((2*(jk+1)-3)/(xr*xr))*fa[jk-1] - dfa[jk-2];        
                    }                    
                    g = -fa[l];
                    gR = -dfa[l];                                        
                    //g  = tmpa[k];
                    //gR = dtmpa[k];        
                    gx = gR*Rx;
                    gy = gR*Ry;
                    gz = gR*Rz; 
                    
                    // derivatives of S_{klm}(x,y,z) = g_{lk}(x,y,z)*Y_{lm}(x,y,z) for l >= 2 and m = 0,1,...,l
                    j = ((l*l + l + m)*K + k)*N + i;      
                    Sr[j] = g*Ylmr; // real part            
                    Si[j] = g*Ylmi; // imag part                                                                                
                    Srx[j] = gx*Ylmr + g*Ylmrx;
                    Sry[j] = gy*Ylmr + g*Ylmry;
                    Srz[j] = gz*Ylmr + g*Ylmrz;
                    Six[j] = gx*Ylmi + g*Ylmix;
                    Siy[j] = gy*Ylmi + g*Ylmiy;
                    Siz[j] = gz*Ylmi + g*Ylmiz;        
                }                
            }
            
            // derivatives of S_{klm}(x,y,z) = g_{lk}(x,y,z)*Y_{lm}(x,y,z) for l >= 2 and m = -1,-2,...,-l            
            C = -1.0;
            for (int m=0; m<l; m++) {  
                for (k=0; k<K; k++) {
                    j  = ((l*l + l + m + 1)*K + k)*N + i;  // spherical harmonics Bessel functions for m > 0 
                    jm = ((l*l + l - m - 1)*K + k)*N + i;  // spherical harmonics Bessel functions for m < 0                                     
                    Sr[jm] = C*Sr[j]; // real part
                    Si[jm] =-C*Si[j]; // imag part                                                                                                                             
                    Srx[jm] = C*Srx[j]; 
                    Six[jm] =-C*Six[j];
                    Sry[jm] = C*Sry[j]; 
                    Siy[jm] =-C*Siy[j];
                    Srz[jm] = C*Srz[j]; 
                    Siz[jm] =-C*Siz[j];                    
                }                      
                C = C*(-1.0);        
            }
        }
        i += blockDim.x * gridDim.x;
    }                        
}
template <typename T> void gpuSphericalHarmonicsBesselWithDeriv(T *Sr, T *Si,
        T *Srx, T *Six, T *Sry, T *Siy, T *Srz, T *Siz, T *x, T *y, T *z, 
        T *x0, T *P, T *tmp, T *f, T *dP, T *dtmp, T *df, T *fac, T pi, int L, int K, int N)
{        
    int blockDim = 256;
    int gridDim = (N + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;
    gpuKernelSphericalHarmonicsBesselWithDeriv<<<gridDim, blockDim>>>(Sr, Si, Srx, Six, Sry, Siy, Srz, Siz, x, y, z,
                 x0, fac, pi, L, K, N);
}
template void gpuSphericalHarmonicsBesselWithDeriv(double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, 
        double*, double*, double*, double*, double*, double*, double*, double*, double, int, int, int);
template void gpuSphericalHarmonicsBesselWithDeriv(float*, float*, float*, float*, float*, float*, float*, float*, float*, float*, float*, 
        float*, float*, float*, float*, float*, float*, float*, float*, float, int, int, int);

// core function
template <typename T> __global__ void gpuKernelSphericalHarmonicsBesselWithDeriv(T *Sr, T *Si,
        T *Srx, T *Six, T *Sry, T *Siy, T *Srz, T *Siz, T *xij, T *x0, T *fac, T pi, int L, int K, int N)
{                        
    // L                      : the maximum degree of spherical harmonics
    // K                      : the maximum degree of spherical Bessel functions
    // N                      : length of x, y, z
    // x0  [(L+1)*K]          : zeros of pherical Bessel functions
    // P   [(L+1)*(L+2)/2]    : temporary storage for the recurrence formula
    // tmp [(L+1)*(L+2)/2]    : temporary storage for the recurrence formula
    // f   [L+1]              : temporary storage for the recurrence formula
    // dP   [(L+1)*(L+2)/2]   : temporary storage for the recurrence formula
    // dtmp [(L+1)*(L+2)/2]   : temporary storage for the recurrence formula
    // df   [L+1]             : temporary storage for the recurrence formula    
    // fac                    : factorial look-up table
    // Srx  [N*K*(L+1)*(L+1)] : x-derivative of real spherical harmonics Bessel functions
    // Six  [N*K*(L+1)*(L+1)] : x-derivative of imag spherical harmonics Bessel functions
    // Sry  [N*K*(L+1)*(L+1)] : y-derivative of real spherical harmonics Bessel functions
    // Siy  [N*K*(L+1)*(L+1)] : y-derivative of imag spherical harmonics Bessel functions
    // Srz  [N*K*(L+1)*(L+1)] : z-derivative of real spherical harmonics Bessel functions
    // Siz  [N*K*(L+1)*(L+1)] : z-derivative of imag spherical harmonics Bessel functions
    
    // Compute partial derivatives of spherical harmonics Bessel functions: S_{klm}(x,y,z) = g_{lk}(x,y,z)*Y_{lm}(x,y,z)
    
    // https://en.wikipedia.org/wiki/Spherical_harmonics
    // Spherical harmonics functions, Y_{lm}(x,y,z), for l = 0,1,...,L and m = -l,...,0,...,l        
    //  l = 0:                0
    //  l = 1:              1 2 3
    //  l = 2:            4 5 6 7 8
    //  l = 3:        9 10 11 12 13 14 15
    //  l = 4:    16 17 18 19 20 21 22 23 24
    //  ....
    //  l = L:  L*L L*L+1 ....   ....  (L+1)*(L+1)-1 
    // The total number of spherical harmonics is (L+1)*(L+1).
    
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
            
    // Hence, the total number of spherical harmonics Bessel functions, S_{klm}(x,y,z), is K*(L+1)*(L+1).            
    
    //#pragma omp parallel for
    int i = threadIdx.x + blockIdx.x * blockDim.x;     
    while (i < N) { // loop over each neighbor atom                               
        int k, j, jm;
        T Pa[Lmax];
        T tmpa[Lmax];
        T fa[Lmax]; 
        T dPa[Lmax];
        T dtmpa[Lmax];
        T dfa[Lmax];
        T dtmpb[Lmax];
        T tmpb[Lmax];

        T x = xij[i*3];
        T y = xij[i*3+1];
        T z = xij[i*3+2];        
                
        if (sqrt(x*x + y*y) < 2e-12) {
            T signx = (x >= 0.0) ? 1.0 : -1.0;
            T signy = (y >= 0.0) ? -1.0 : 1.0;
            x = fabs(x) > 1e-12 ? x : signx*1e-12;
            y = fabs(y) > 1e-12 ? y : signy*1e-12;
        }
                
        // Cartersian -> Spherical
        T r = sqrt(x*x + y*y + z*z);
        //T the = acos(z/r);
        T phi = atan2(y,x);
                                
        // partial derivatives of spherical coordinates
        T r2 = r*r;
        T rxy = sqrt(x*x + y*y);
        T rxy2 = rxy*rxy;
        T rr2 = rxy*r2;
        T Rx = x/r;
        T Ry = y/r;
        T Rz = z/r;
        T Thex = x*z/rr2;
        T They = y*z/rr2;
        T Thez = -rxy/r2;
        T Phix = -y/rxy2;
        T Phiy = x/rxy2;
        T Phiz = 0.0;        
        
        // Spherical harmonics Y_{lm}(x,y,z) and their derivatives for l = 0        
        int l = 0;
        T Ylmr = sqrt(1/(4*pi));
        T Ylmi = 0.0;
        T YlmrThe = 0.0;
        T YlmrPhi = 0.0;
        T YlmiThe = 0.0;        
        T YlmiPhi = 0.0;        
        T Ylmrx = YlmrThe*Thex + YlmrPhi*Phix;
        T Ylmry = YlmrThe*They + YlmrPhi*Phiy;
        T Ylmrz = YlmrThe*Thez + YlmrPhi*Phiz;
        T Ylmix = YlmiThe*Thex + YlmiPhi*Phix;
        T Ylmiy = YlmiThe*They + YlmiPhi*Phiy;
        T Ylmiz = YlmiThe*Thez + YlmiPhi*Phiz;
        
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
            
            // derivatives of S_{klm}(x,y,z) = g_{lk}(x,y,z)*Y_{lm}(x,y,z) for l = 0
            j = k*N + i;                
            Sr[j] = g*Ylmr;  // real part                   
            Si[j] = g*Ylmi;  // imag part                   
            Srx[j] = gx*Ylmr + g*Ylmrx;
            Sry[j] = gy*Ylmr + g*Ylmry;
            Srz[j] = gz*Ylmr + g*Ylmrz;
            Six[j] = gx*Ylmi + g*Ylmix;
            Siy[j] = gy*Ylmi + g*Ylmiy;
            Siz[j] = gz*Ylmi + g*Ylmiz;            
            //cout<<i<<" "<<j<<" "<<Sr[j]<<"  "<<Srx[j]<<" "<<gx<<"  "<<Ylmrx<<"  "<<Thex<<"  "<<Phix<<endl;
        }
                        
        l = 1;
//         T costhe = cos(the);
//         T a = -sin(the);        
        T costhe = z/r; 
        T a = -sqrt((r-z)*(r+z))/r;     
        T dcosthe = a;
        T da = -costhe;        
        int m = 0;    
        Pa[0] = costhe;
        Pa[1] = a;
        dPa[0] = dcosthe;    
        dPa[1] = da;
        tmpa[0] = 1.0;
        dtmpa[0] = 0.0;
        // Spherical harmonics Y_{lm}(x,y,z) and their derivatives for l = 1 and m = 0  
        m = 0;
        T C = sqrt((2*l+1)*fac[l-m]/(4*pi*fac[l+m]));
        Ylmr = C*Pa[0];
        Ylmi = 0.0;
        YlmrThe = C*a;    
        YlmiThe = 0.0;    
        YlmrPhi = 0.0;          
        YlmiPhi = 0.0;                  
        Ylmrx = YlmrThe*Thex + YlmrPhi*Phix;
        Ylmry = YlmrThe*They + YlmrPhi*Phiy;
        Ylmrz = YlmrThe*Thez + YlmrPhi*Phiz;
        Ylmix = YlmiThe*Thex + YlmiPhi*Phix;
        Ylmiy = YlmiThe*They + YlmiPhi*Phiy;
        Ylmiz = YlmiThe*Thez + YlmiPhi*Phiz;        
        //printf("%g \n", Ylmr);
        //printf("%d %g %g %g %g %g %g %g\n", i, C, x, y, z, r, costhe, Pa[0], Ylmr);
        for (k=0; k<K; k++) {
            // Spherical Bessel g_{lk}(x,y,z) and their derivatives for l = 1
            xr = x0[k*(L+1)+l]*r;                        
            g = cos(xr)/(xr*xr) + sin(xr)/xr;
            gR = x0[k*(L+1)+l]*(cos(xr)/xr - (2*cos(xr))/(xr*xr*xr) - (2*sin(xr))/(xr*xr));
            gx = gR*Rx;
            gy = gR*Ry;
            gz = gR*Rz;            
                        
            // derivatives of S_{klm}(x,y,z) = g_{lk}(x,y,z)*Y_{lm}(x,y,z) for l = 1 and m = 0
            j = ((l*l + l + m)*K + k)*N + i;      
            Sr[j] = g*Ylmr; // real part           
            Si[j] = g*Ylmi; // imag part                                                
            Srx[j] = gx*Ylmr + g*Ylmrx;
            Sry[j] = gy*Ylmr + g*Ylmry;
            Srz[j] = gz*Ylmr + g*Ylmrz;
            Six[j] = gx*Ylmi + g*Ylmix;
            Siy[j] = gy*Ylmi + g*Ylmiy;
            Siz[j] = gz*Ylmi + g*Ylmiz;            
        }        
        
        // Spherical harmonics Y_{lm}(x,y,z) and their derivatives for l = 1 and m = 1
        m = 1;
        C = sqrt((2*l+1)*fac[l-m]/(4*pi*fac[l+m]));
        Ylmr = C*cos(phi)*Pa[1];
        Ylmi = C*sin(phi)*Pa[1];                
        YlmrThe = C*(cos(phi)*da);      
        YlmiThe = C*(sin(phi)*da);      
        YlmrPhi = -C*(sin(phi)*Pa[1]); 
        YlmiPhi = C*(cos(phi)*Pa[1]);  
        Ylmrx = YlmrThe*Thex + YlmrPhi*Phix;
        Ylmry = YlmrThe*They + YlmrPhi*Phiy;
        Ylmrz = YlmrThe*Thez + YlmrPhi*Phiz;
        Ylmix = YlmiThe*Thex + YlmiPhi*Phix;
        Ylmiy = YlmiThe*They + YlmiPhi*Phiy;
        Ylmiz = YlmiThe*Thez + YlmiPhi*Phiz;
        //printf("%d %g %g %g %g %g %g %g %g\n", i, C, x, y, z, r, cos(phi), Pa[1], Ylmr);
        for (k=0; k<K; k++) {
            // Spherical Bessel g_{lk}(x,y,z) and their derivatives for l = 1
            xr = x0[k*(L+1)+l]*r;                        
            g = cos(xr)/(xr*xr) + sin(xr)/xr;
            gR = x0[k*(L+1)+l]*(cos(xr)/xr - (2*cos(xr))/(xr*xr*xr) - (2*sin(xr))/(xr*xr));
            gx = gR*Rx;
            gy = gR*Ry;
            gz = gR*Rz;            
            
            // derivatives of S_{klm}(x,y,z) = g_{lk}(x,y,z)*Y_{lm}(x,y,z) for l = 1 and m = 1
            j = ((l*l + l + m)*K + k)*N + i;          
            Sr[j] = g*Ylmr; // real part           
            Si[j] = g*Ylmi; // imag part                                                                        
            Srx[j] = gx*Ylmr + g*Ylmrx;
            Sry[j] = gy*Ylmr + g*Ylmry;
            Srz[j] = gz*Ylmr + g*Ylmrz;
            Six[j] = gx*Ylmi + g*Ylmix;
            Siy[j] = gy*Ylmi + g*Ylmiy;
            Siz[j] = gz*Ylmi + g*Ylmiz;        
            
            // derivatives of S_{klm}(x,y,z) = g_{lk}(x,y,z)*Y_{lm}(x,y,z) for l = 1 and m = -1
            jm = ((l*l + l - m)*K + k)*N + i;       
            Sr[jm] = -Sr[j]; // real part                      
            Si[jm] =  Si[j]; // imag part                                                                                   
            Srx[jm] = -Srx[j];
            Sry[jm] = -Sry[j];
            Srz[jm] = -Srz[j];
            Six[jm] = Six[j];
            Siy[jm] = Siy[j];
            Siz[jm] = Siz[j];            
        }        
                
        for (l=2; l<=L; l++) {        
            // Compute associated Legendre polynomials P_m and their derivatives using recurrence formula
            // https://en.wikipedia.org/wiki/Associated_Legendre_polynomials                        
            T Pll = Pa[l-1];
            tmpa[(l-1)] = Pll;
            Pa[(l-1)] = (2*(l-1)+1)*(costhe*Pll);
            Pa[l] = (2*(l-1)+1)*(a*Pll); 
            T dPll = dPa[l-1];
            dtmpa[l-1] = dPll;
            dPa[l-1] = (2*(l-1)+1)*(costhe*dPll + dcosthe*Pll);   
            dPa[l] = (2*(l-1)+1)*(a*dPll + da*Pll);             
            for (m=0; m<l-1; m++) {
                Pll = Pa[m];
                Pa[m] = ((2*(l-1)+1)*(costhe*Pll) - (l+m-1)*tmpa[m])/(T (l-m));
                tmpa[m] = Pll;
                dPll = dPa[m];
                dPa[m] = ((2*(l-1)+1)*(costhe*dPll+dcosthe*Pll) - (l+m-1)*dtmpa[m])/(T (l-m));
                dtmpa[m] = dPll;
            }
                            
            for (k=0; k<K; k++) {
                // Spherical Bessel g_{lk}(x,y,z) and their derivatives for l >=2
                xr = x0[k*(L+1)+l]*r;                                            
                fa[0] = -cos(xr)/xr;
                fa[1] = -cos(xr)/(xr*xr) - sin(xr)/xr;
                dfa[0] = x0[k*(L+1)+l]*(cos(xr)/(xr*xr) + sin(xr)/xr);
                dfa[1] = x0[k*(L+1)+l]*((2*cos(xr))/(xr*xr*xr) - cos(xr)/xr  + (2*sin(xr))/(xr*xr));
                for (int jk=2; jk<(l+1); jk++) {
                    fa[jk] = ((2*(jk+1)-3)/xr)*fa[jk-1] - fa[jk-2];                 
                    dfa[jk] = ((2*(jk+1)-3)/xr)*dfa[jk-1] - x0[k*(L+1)+l]*((2*(jk+1)-3)/(xr*xr))*fa[jk-1] - dfa[jk-2];        
                }
                tmpb[k] = -fa[l];
                dtmpb[k] = -dfa[l];                                    
            }

            // Compute spherical harmonics Bessel functions and their derivatives at level l for m = 0, 1, ..., l
            for (m=0; m<=l; m++) {
                // Spherical harmonics Y_{lm}(x,y,z) and their derivatives for l and m = 0,1,..,l  
                C = sqrt((2*l+1)*fac[l-m]/(4*pi*fac[l+m])); 
                Ylmr = C*(cos(m*phi)*Pa[m]);
                Ylmi = C*(sin(m*phi)*Pa[m]);
                YlmrThe = C*(cos(m*phi)*dPa[m]); 
                YlmiThe = C*(sin(m*phi)*dPa[m]); 
                YlmrPhi = -(m*C)*(sin(m*phi)*Pa[m]); 
                YlmiPhi = (m*C)*(cos(m*phi)*Pa[m]);          
                Ylmrx = YlmrThe*Thex + YlmrPhi*Phix;
                Ylmry = YlmrThe*They + YlmrPhi*Phiy;
                Ylmrz = YlmrThe*Thez + YlmrPhi*Phiz;
                Ylmix = YlmiThe*Thex + YlmiPhi*Phix;
                Ylmiy = YlmiThe*They + YlmiPhi*Phiy;
                Ylmiz = YlmiThe*Thez + YlmiPhi*Phiz; 

                for (k=0; k<K; k++) {
                    // Spherical Bessel g_{lk}(x,y,z) and their derivatives for l >=2
//                     xr = x0[k*(L+1)+l]*r;                                            
//                     fa[0] = -cos(xr)/xr;
//                     fa[1] = -cos(xr)/(xr*xr) - sin(xr)/xr;
//                     dfa[0] = x0[k*(L+1)+l]*(cos(xr)/(xr*xr) + sin(xr)/xr);
//                     dfa[1] = x0[k*(L+1)+l]*((2*cos(xr))/(xr*xr*xr) - cos(xr)/xr  + (2*sin(xr))/(xr*xr));
//                     for (int jk=2; jk<(l+1); jk++) {
//                         fa[jk] = ((2*(jk+1)-3)/xr)*fa[jk-1] - fa[jk-2];                 
//                         dfa[jk] = ((2*(jk+1)-3)/xr)*dfa[jk-1] - x0[k*(L+1)+l]*((2*(jk+1)-3)/(xr*xr))*fa[jk-1] - dfa[jk-2];        
//                     }
//                     g = -fa[l];
//                     gR = -dfa[l];   
                    g  = tmpb[k];
                    gR = dtmpb[k];
                    gx = gR*Rx;
                    gy = gR*Ry;
                    gz = gR*Rz; 
                    
                    // derivatives of S_{klm}(x,y,z) = g_{lk}(x,y,z)*Y_{lm}(x,y,z) for l >= 2 and m = 0,1,...,l
                    j = ((l*l + l + m)*K + k)*N + i;      
                    Sr[j] = g*Ylmr; // real part            
                    Si[j] = g*Ylmi; // imag part                                                                                
                    Srx[j] = gx*Ylmr + g*Ylmrx;
                    Sry[j] = gy*Ylmr + g*Ylmry;
                    Srz[j] = gz*Ylmr + g*Ylmrz;
                    Six[j] = gx*Ylmi + g*Ylmix;
                    Siy[j] = gy*Ylmi + g*Ylmiy;
                    Siz[j] = gz*Ylmi + g*Ylmiz;        
                }                
            }
            
            // derivatives of S_{klm}(x,y,z) = g_{lk}(x,y,z)*Y_{lm}(x,y,z) for l >= 2 and m = -1,-2,...,-l            
            C = -1.0;
            for (int m=0; m<l; m++) {  
                for (k=0; k<K; k++) {
                    j  = ((l*l + l + m + 1)*K + k)*N + i;  // spherical harmonics Bessel functions for m > 0 
                    jm = ((l*l + l - m - 1)*K + k)*N + i;  // spherical harmonics Bessel functions for m < 0                                     
                    Sr[jm] = C*Sr[j]; // real part
                    Si[jm] =-C*Si[j]; // imag part                                                                                                                             
                    Srx[jm] = C*Srx[j]; 
                    Six[jm] =-C*Six[j];
                    Sry[jm] = C*Sry[j]; 
                    Siy[jm] =-C*Siy[j];
                    Srz[jm] = C*Srz[j]; 
                    Siz[jm] =-C*Siz[j];                    
                }                      
                C = C*(-1.0);        
            }
        }
        i += blockDim.x * gridDim.x;
    }                        
}
template <typename T> void gpuSphericalHarmonicsBesselWithDeriv(T *Sr, T *Si,
        T *Srx, T *Six, T *Sry, T *Siy, T *Srz, T *Siz, T *xij, 
        T *x0, T *P, T *tmp, T *f, T *dP, T *dtmp, T *df, T *fac, T pi, int L, int K, int N)
{        
    int blockDim = 256;
    int gridDim = (N + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;
    gpuKernelSphericalHarmonicsBesselWithDeriv<<<gridDim, blockDim>>>(Sr, Si, Srx, Six, Sry, Siy, Srz, Siz, xij,
                 x0, fac, pi, L, K, N);
}
template void gpuSphericalHarmonicsBesselWithDeriv(double*, double*, double*, double*, double*, 
        double*, double*, double*, double*, double*, double*, double*, double*, double*, 
        double*, double*, double*, double, int, int, int);
template void gpuSphericalHarmonicsBesselWithDeriv(float*, float*, float*, float*, float*, 
        float*, float*, float*, float*, float*, float*, float*, float*, 
        float*, float*, float*, float*, float, int, int, int);


// core function
template <typename T> __global__ void gpuKernelSphericalHarmonicsBesselJ(T *Sr, T *Si, T *xij, 
                T *x0, T *fac, T pi, int L, int K, int N)
{                        
    // L                     : the maximum degree of spherical harmonics
    // K                     : the maximum degree of spherical Bessel functions
    // N                     : length of x, y, z
    // x0  [(L+1)*K]         : zeros of pherical Bessel functions
    // P   [(L+1)*(L+2)/2]   : temporary storage for the recurrence formula
    // tmp [(L+1)*(L+2)/2]   : temporary storage for the recurrence formula
    // f   [L+1]             : temporary storage for the recurrence formula
    // fac                   : factorial look-up table
    // Sr  [N*K*(L+1)*(L+1)] : real part of spherical harmonics Bessel functions
    // Si  [N*K*(L+1)*(L+1)] : imag part of spherical harmonics Bessel functions
        
    // Compute spherical harmonics Bessel functions: S_{klm}(x,y,z) = g_{lk}(x,y,z)*Y_{lm}(x,y,z)
    
    // https://en.wikipedia.org/wiki/Spherical_harmonics
    // Spherical harmonics functions, Y_{lm}(x,y,z), for l = 0,1,...,L and m = -l,...,0,...,l        
    //  l = 0:                0
    //  l = 1:              1 2 3
    //  l = 2:            4 5 6 7 8
    //  l = 3:        9 10 11 12 13 14 15
    //  l = 4:    16 17 18 19 20 21 22 23 24
    //  ....
    //  l = L:  L*L L*L+1 ....   ....  (L+1)*(L+1)-1 
    // The total number of spherical harmonics is (L+1)*(L+1).
    
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
        
    // Hence, the total number of spherical harmonics Bessel functions is K*(L+1)*(L+1).            
            
    //#pragma omp parallel for    
    int i = threadIdx.x + blockIdx.x * blockDim.x;     
    while (i < N) { // loop over each neighbor atom               
        int k, j, jm;
        T Pa[Lmax];
        T tmpa[Lmax];
        T fa[Lmax]; 

        tmpa[0] = 1.0;
        
        T x = xij[i*3];
        T y = xij[i*3+1];
        T z = xij[i*3+2];        
        
        if (sqrt(x*x + y*y) < 2e-12) {
            T signx = (x >= 0.0) ? 1.0 : -1.0;
            T signy = (y >= 0.0) ? -1.0 : 1.0;
            x = fabs(x) > 1e-12 ? x : signx*1e-12;
            y = fabs(y) > 1e-12 ? y : signy*1e-12;
        }
                
        // Cartersian -> Spherical
        T r = sqrt(x*x + y*y + z*z);
        //T the = acos(z/r);
        T phi = atan2(y,x);
        
        // Spherical harmonics for l = 0        
        int l = 0;
        T Ylmr = sqrt(1/(4*pi));
        T Ylmi = 0.0;                
        T xr, g;
        // Spherical Bessel for l = 0, see https://en.wikipedia.org/wiki/Bessel_function#Spherical_Bessel_functions:_jn,_yn 
        for (k=0; k<K; k++) {
            xr = x0[k*(L+1)+l]*r;                        
            g = sin(xr)/xr;            
            j = k*N + i;                      
            Sr[j] = g*Ylmr;  // real part                   
            Si[j] = g*Ylmi;  // imag part       
        }
                
        // Spherical harmonics for l = 1;
        l = 1;
//         T costhe = cos(the);
//         T a = -sin(the);        
        T costhe = z/r; 
        T a = -sqrt((r-z)*(r+z))/r;             
        int m = 0;    
        Pa[0] = costhe;
        Pa[1] = a;
        
        T C = sqrt((2*l+1)*fac[l-m]/(4*pi*fac[l+m]));
        Ylmr = C*Pa[0];
        Ylmi = 0.0;
        // Spherical Bessel for l = 1, see https://en.wikipedia.org/wiki/Bessel_function#Spherical_Bessel_functions:_jn,_yn
        for (k=0; k<K; k++) {
            xr = x0[k*(L+1)+l]*r;                        
            g = sin(xr)/(xr*xr) - cos(xr)/xr;
            j = ((l*l + l + m)*K + k)*N + i;                
            Sr[j] = g*Ylmr; // real part           
            Si[j] = g*Ylmi; // imag part                                    
        }        
        
        m = 1;
        C = sqrt((2*l+1)*fac[l-m]/(4*pi*fac[l+m]));
        Ylmr = C*cos(phi)*Pa[1];
        Ylmi = C*sin(phi)*Pa[1];                
        //cout<<m<<"  "<<C<<"  "<<cos(phi)<<"  "<<Pa[1]<<endl;
        // Spherical Bessel for l = 1, see https://en.wikipedia.org/wiki/Bessel_function#Spherical_Bessel_functions:_jn,_yn
        for (k=0; k<K; k++) {
            xr = x0[k*(L+1)+l]*r;                        
            g = sin(xr)/(xr*xr) - cos(xr)/xr; 
            j = ((l*l + l + m)*K + k)*N + i; // spherical harmonics Bessel functions for m > 0                
            Sr[j] = g*Ylmr; // real part           
            Si[j] = g*Ylmi; // imag part              
            //cout<<j<<"  "<<g<<"  "<<Ylmr<<"  "<<Sr[j]<<endl;
            jm = ((l*l + l - m)*K + k)*N + i; // spherical harmonics Bessel functions for m < 0                
            Sr[jm] = -Sr[j]; // real part                      
            Si[jm] =  Si[j]; // imag part                                                                       
        }        
                
        for (l=2; l<=L; l++) {                                        
            // Compute associated Legendre polynomial using recurrence formula
            // https://en.wikipedia.org/wiki/Associated_Legendre_polynomials            
            T Pll = Pa[l-1];
            tmpa[(l-1)] = Pll;
            Pa[(l-1)] = (2*(l-1)+1)*(costhe*Pll);
            Pa[l] = (2*(l-1)+1)*(a*Pll);             
            for (m=0; m<l-1; m++) {
                Pll = Pa[m];
                Pa[m] = ((2*(l-1)+1)*(costhe*Pll) - (l+m-1)*tmpa[m])/(T (l-m));
                tmpa[m] = Pll;                
            }
                            
//             for (k=0; k<K; k++) {
//                 // Compute the spherical Bessel functions using recurrence formula
//                 xr = x0[k*(L+1)+l]*r;                                            
//                 fa[0] = sin(xr)/xr;
//                 fa[1] = sin(xr)/(xr*xr) - cos(xr)/xr;                    
//                 for (int jk=2; jk<(l+1); jk++) 
//                     fa[jk] = ((2*(jk+1)-3)/xr)*fa[jk-1] - fa[jk-2];                                     
//                 tmpa[k] = fa[l];                                                                              
//             }                

            // Compute spherical harmonics Bessel function at level l for m = 0, 1, ..., l
            for (m=0; m<=l; m++) {
                C = sqrt((2*l+1)*fac[l-m]/(4*pi*fac[l+m])); 
                Ylmr = C*(cos(m*phi)*Pa[m]);
                Ylmi = C*(sin(m*phi)*Pa[m]);                
                for (k=0; k<K; k++) {
                    // Compute the spherical Bessel functions using recurrence formula
                    xr = x0[k*(L+1)+l]*r;                                            
                    fa[0] = sin(xr)/xr;
                    fa[1] = sin(xr)/(xr*xr) - cos(xr)/xr;                    
                    for (int jk=2; jk<(l+1); jk++) 
                        fa[jk] = ((2*(jk+1)-3)/xr)*fa[jk-1] - fa[jk-2];                                     
                    g = fa[l];                                                                                  
                    //g = tmpa[k];                                                              
                    j = ((l*l + l + m)*K + k)*N + i;   
                    Sr[j] = g*Ylmr; // real part            
                    Si[j] = g*Ylmi; // imag part                                                            
                }                
            }        
            
            // Compute the spherical harmonics Bessel functions for m < 0 using symmetry properties
            C = -1.0;
            for (int m=1; m<=l; m++)  {                 
                for (k=0; k<K; k++) {
                    j =  ((l*l + l + m)*K + k)*N + i;  // spherical harmonics Bessel functions for m > 0 
                    jm = ((l*l + l - m)*K + k)*N + i;  // spherical harmonics Bessel functions for m < 0                 
                    Sr[jm] = C*Sr[j]; // real part
                    Si[jm] =-C*Si[j]; // imag part                                                                             
                }                     
                C = C*(-1.0);        
            }
        }
        i += blockDim.x * gridDim.x;
    }          
}
template <typename T> void gpuSphericalHarmonicsBesselJ(T *Sr, T *Si, T *xij, 
                T *x0, T *P, T *tmp, T *f, T *fac, T pi, int L, int K, int N)
{        
    int blockDim = 256;
    int gridDim = (N + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;
    gpuKernelSphericalHarmonicsBesselJ<<<gridDim, blockDim>>>(Sr, Si, xij, 
            x0, fac, pi, L, K, N);
}
template void gpuSphericalHarmonicsBesselJ(double*, double*, double*, double*, 
        double*, double*, double*, double*, double, int, int, int);
template void gpuSphericalHarmonicsBesselJ(float*, float*, float*, float*, 
        float*, float*, float*, float*, float, int, int, int);

// core function
template <typename T> __global__ void gpuKernelSphericalHarmonicsBesselJWithDeriv(T *Sr, T *Si,
        T *Srx, T *Six, T *Sry, T *Siy, T *Srz, T *Siz, T *xij, T *x0, 
                T *fac, T pi, int L, int K, int N)
{                        
    // L                      : the maximum degree of spherical harmonics
    // K                      : the maximum degree of spherical Bessel functions
    // N                      : length of x, y, z
    // x0  [(L+1)*K]          : zeros of pherical Bessel functions
    // P   [(L+1)*(L+2)/2]    : temporary storage for the recurrence formula
    // tmp [(L+1)*(L+2)/2]    : temporary storage for the recurrence formula
    // f   [L+1]              : temporary storage for the recurrence formula
    // dP   [(L+1)*(L+2)/2]   : temporary storage for the recurrence formula
    // dtmp [(L+1)*(L+2)/2]   : temporary storage for the recurrence formula
    // df   [L+1]             : temporary storage for the recurrence formula    
    // fac                    : factorial look-up table
    // Srx  [N*K*(L+1)*(L+1)] : x-derivative of real spherical harmonics Bessel functions
    // Six  [N*K*(L+1)*(L+1)] : x-derivative of imag spherical harmonics Bessel functions
    // Sry  [N*K*(L+1)*(L+1)] : y-derivative of real spherical harmonics Bessel functions
    // Siy  [N*K*(L+1)*(L+1)] : y-derivative of imag spherical harmonics Bessel functions
    // Srz  [N*K*(L+1)*(L+1)] : z-derivative of real spherical harmonics Bessel functions
    // Siz  [N*K*(L+1)*(L+1)] : z-derivative of imag spherical harmonics Bessel functions
    
    // Compute partial derivatives of spherical harmonics Bessel functions: S_{klm}(x,y,z) = g_{lk}(x,y,z)*Y_{lm}(x,y,z)
    
    // https://en.wikipedia.org/wiki/Spherical_harmonics
    // Spherical harmonics functions, Y_{lm}(x,y,z), for l = 0,1,...,L and m = -l,...,0,...,l        
    //  l = 0:                0
    //  l = 1:              1 2 3
    //  l = 2:            4 5 6 7 8
    //  l = 3:        9 10 11 12 13 14 15
    //  l = 4:    16 17 18 19 20 21 22 23 24
    //  ....
    //  l = L:  L*L L*L+1 ....   ....  (L+1)*(L+1)-1 
    // The total number of spherical harmonics is (L+1)*(L+1).
    
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
            
    // Hence, the total number of spherical harmonics Bessel functions, S_{klm}(x,y,z), is K*(L+1)*(L+1).            
    
    //#pragma omp parallel for
    int i = threadIdx.x + blockIdx.x * blockDim.x;     
    while (i < N) { // loop over each neighbor atom                               
        int k, j, jm;
        T Pa[Lmax];
        T tmpa[Lmax];
        T fa[Lmax]; 
        T dPa[Lmax];
        T dtmpa[Lmax];
        T dfa[Lmax];
        
        T x = xij[i*3];
        T y = xij[i*3+1];
        T z = xij[i*3+2];        
                
        if (sqrt(x*x + y*y) < 2e-12) {
            T signx = (x >= 0.0) ? 1.0 : -1.0;
            T signy = (y >= 0.0) ? -1.0 : 1.0;
            x = fabs(x) > 1e-12 ? x : signx*1e-12;
            y = fabs(y) > 1e-12 ? y : signy*1e-12;
        }
        
        // Cartersian -> Spherical
        T r = sqrt(x*x + y*y + z*z);
        //T the = acos(z/r);
        T phi = atan2(y,x);
                                
        // partial derivatives of spherical coordinates
        T r2 = r*r;
        T rxy = sqrt(x*x + y*y);
        T rxy2 = rxy*rxy;
        T rr2 = rxy*r2;
        T Rx = x/r;
        T Ry = y/r;
        T Rz = z/r;
        T Thex = x*z/rr2;
        T They = y*z/rr2;
        T Thez = -rxy/r2;
        T Phix = -y/rxy2;
        T Phiy = x/rxy2;
        T Phiz = 0.0;        
        
        // Spherical harmonics Y_{lm}(x,y,z) and their derivatives for l = 0        
        int l = 0;
        T Ylmr = sqrt(1/(4*pi));
        T Ylmi = 0.0;
        T YlmrThe = 0.0;
        T YlmrPhi = 0.0;
        T YlmiThe = 0.0;        
        T YlmiPhi = 0.0;        
        T Ylmrx = YlmrThe*Thex + YlmrPhi*Phix;
        T Ylmry = YlmrThe*They + YlmrPhi*Phiy;
        T Ylmrz = YlmrThe*Thez + YlmrPhi*Phiz;
        T Ylmix = YlmiThe*Thex + YlmiPhi*Phix;
        T Ylmiy = YlmiThe*They + YlmiPhi*Phiy;
        T Ylmiz = YlmiThe*Thez + YlmiPhi*Phiz;
        
        T xr, g, gR, gx, gy, gz;
        for (k=0; k<K; k++) {
            // https://en.wikipedia.org/wiki/Bessel_function#Spherical_Bessel_functions:_jn,_yn             
            // Spherical Bessel g_{lk}(x,y,z) and their derivatives for l = 0             
            xr = x0[k*(L+1)+l]*r;                        
            g = sin(xr)/xr;
            gR = x0[k*(L+1)+l]*(cos(xr)/xr - sin(xr)/(xr*xr));
            gx = gR*Rx;
            gy = gR*Ry;
            gz = gR*Rz;            
            
            // derivatives of S_{klm}(x,y,z) = g_{lk}(x,y,z)*Y_{lm}(x,y,z) for l = 0
            j = k*N + i;                
            Sr[j] = g*Ylmr;  // real part                   
            Si[j] = g*Ylmi;  // imag part                   
            Srx[j] = gx*Ylmr + g*Ylmrx;
            Sry[j] = gy*Ylmr + g*Ylmry;
            Srz[j] = gz*Ylmr + g*Ylmrz;
            Six[j] = gx*Ylmi + g*Ylmix;
            Siy[j] = gy*Ylmi + g*Ylmiy;
            Siz[j] = gz*Ylmi + g*Ylmiz;            
            //cout<<i<<" "<<j<<" "<<Sr[j]<<"  "<<Srx[j]<<" "<<gx<<"  "<<Ylmrx<<"  "<<Thex<<"  "<<Phix<<endl;
        }
                        
        l = 1;
//         T costhe = cos(the);
//         T a = -sin(the);        
        T costhe = z/r; 
        T a = -sqrt((r-z)*(r+z))/r;     
        T dcosthe = a;
        T da = -costhe;        
        int m = 0;    
        Pa[0] = costhe;
        Pa[1] = a;
        dPa[0] = dcosthe;    
        dPa[1] = da;
        tmpa[0] = 1.0;
        dtmpa[0] = 0.0;
        // Spherical harmonics Y_{lm}(x,y,z) and their derivatives for l = 1 and m = 0  
        m = 0;
        T C = sqrt((2*l+1)*fac[l-m]/(4*pi*fac[l+m]));
        Ylmr = C*Pa[0];
        Ylmi = 0.0;
        YlmrThe = C*a;    
        YlmiThe = 0.0;    
        YlmrPhi = 0.0;          
        YlmiPhi = 0.0;                  
        Ylmrx = YlmrThe*Thex + YlmrPhi*Phix;
        Ylmry = YlmrThe*They + YlmrPhi*Phiy;
        Ylmrz = YlmrThe*Thez + YlmrPhi*Phiz;
        Ylmix = YlmiThe*Thex + YlmiPhi*Phix;
        Ylmiy = YlmiThe*They + YlmiPhi*Phiy;
        Ylmiz = YlmiThe*Thez + YlmiPhi*Phiz;        
        for (k=0; k<K; k++) {
            // Spherical Bessel g_{lk}(x,y,z) and their derivatives for l = 1
            xr = x0[k*(L+1)+l]*r;                        
            g = sin(xr)/(xr*xr) - cos(xr)/xr;
            gR = x0[k*(L+1)+l]*(sin(xr)/xr - (2*sin(xr))/(xr*xr*xr) + (2*cos(xr))/(xr*xr));
            gx = gR*Rx;
            gy = gR*Ry;
            gz = gR*Rz;            
                        
            // derivatives of S_{klm}(x,y,z) = g_{lk}(x,y,z)*Y_{lm}(x,y,z) for l = 1 and m = 0
            j = ((l*l + l + m)*K + k)*N + i;      
            Sr[j] = g*Ylmr; // real part           
            Si[j] = g*Ylmi; // imag part                                                
            Srx[j] = gx*Ylmr + g*Ylmrx;
            Sry[j] = gy*Ylmr + g*Ylmry;
            Srz[j] = gz*Ylmr + g*Ylmrz;
            Six[j] = gx*Ylmi + g*Ylmix;
            Siy[j] = gy*Ylmi + g*Ylmiy;
            Siz[j] = gz*Ylmi + g*Ylmiz;            
        }        
        
        // Spherical harmonics Y_{lm}(x,y,z) and their derivatives for l = 1 and m = 1
        m = 1;
        C = sqrt((2*l+1)*fac[l-m]/(4*pi*fac[l+m]));
        Ylmr = C*cos(phi)*Pa[1];
        Ylmi = C*sin(phi)*Pa[1];                
        YlmrThe = C*(cos(phi)*da);      
        YlmiThe = C*(sin(phi)*da);      
        YlmrPhi = -C*(sin(phi)*Pa[1]); 
        YlmiPhi = C*(cos(phi)*Pa[1]);  
        Ylmrx = YlmrThe*Thex + YlmrPhi*Phix;
        Ylmry = YlmrThe*They + YlmrPhi*Phiy;
        Ylmrz = YlmrThe*Thez + YlmrPhi*Phiz;
        Ylmix = YlmiThe*Thex + YlmiPhi*Phix;
        Ylmiy = YlmiThe*They + YlmiPhi*Phiy;
        Ylmiz = YlmiThe*Thez + YlmiPhi*Phiz;
        for (k=0; k<K; k++) {
            // Spherical Bessel g_{lk}(x,y,z) and their derivatives for l = 1
            xr = x0[k*(L+1)+l]*r;                        
            g = sin(xr)/(xr*xr) - cos(xr)/xr;
            gR = x0[k*(L+1)+l]*(sin(xr)/xr - (2*sin(xr))/(xr*xr*xr) + (2*cos(xr))/(xr*xr));
            gx = gR*Rx;
            gy = gR*Ry;
            gz = gR*Rz;            
            
            // derivatives of S_{klm}(x,y,z) = g_{lk}(x,y,z)*Y_{lm}(x,y,z) for l = 1 and m = 1
            j = ((l*l + l + m)*K + k)*N + i;          
            Sr[j] = g*Ylmr; // real part           
            Si[j] = g*Ylmi; // imag part                                                                        
            Srx[j] = gx*Ylmr + g*Ylmrx;
            Sry[j] = gy*Ylmr + g*Ylmry;
            Srz[j] = gz*Ylmr + g*Ylmrz;
            Six[j] = gx*Ylmi + g*Ylmix;
            Siy[j] = gy*Ylmi + g*Ylmiy;
            Siz[j] = gz*Ylmi + g*Ylmiz;        
            
            // derivatives of S_{klm}(x,y,z) = g_{lk}(x,y,z)*Y_{lm}(x,y,z) for l = 1 and m = -1
            jm = ((l*l + l - m)*K + k)*N + i;       
            Sr[jm] = -Sr[j]; // real part                      
            Si[jm] =  Si[j]; // imag part                                                                                   
            Srx[jm] = -Srx[j];
            Sry[jm] = -Sry[j];
            Srz[jm] = -Srz[j];
            Six[jm] = Six[j];
            Siy[jm] = Siy[j];
            Siz[jm] = Siz[j];            
        }        
                
        for (l=2; l<=L; l++) {        
            // Compute associated Legendre polynomials P_m and their derivatives using recurrence formula
            // https://en.wikipedia.org/wiki/Associated_Legendre_polynomials                        
            T Pll = Pa[l-1];
            tmpa[(l-1)] = Pll;
            Pa[(l-1)] = (2*(l-1)+1)*(costhe*Pll);
            Pa[l] = (2*(l-1)+1)*(a*Pll); 
            T dPll = dPa[l-1];
            dtmpa[l-1] = dPll;
            dPa[l-1] = (2*(l-1)+1)*(costhe*dPll + dcosthe*Pll);   
            dPa[l] = (2*(l-1)+1)*(a*dPll + da*Pll);             
            for (m=0; m<l-1; m++) {
                Pll = Pa[m];
                Pa[m] = ((2*(l-1)+1)*(costhe*Pll) - (l+m-1)*tmpa[m])/(T (l-m));
                tmpa[m] = Pll;
                dPll = dPa[m];
                dPa[m] = ((2*(l-1)+1)*(costhe*dPll+dcosthe*Pll) - (l+m-1)*dtmpa[m])/(T (l-m));
                dtmpa[m] = dPll;
            }
                            
//             for (k=0; k<K; k++) {
//                 // Spherical Bessel g_{lk}(x,y,z) and their derivatives for l >=2
//                 xr = x0[k*(L+1)+l]*r;                                            
//                 fa[0] = sin(xr)/xr;
//                 fa[1] = sin(xr)/(xr*xr) - cos(xr)/xr;                    
//                 dfa[0] = x0[k*(L+1)+l]*(cos(xr)/xr - sin(xr)/(xr*xr));
//                 dfa[1] = x0[k*(L+1)+l]*(sin(xr)/xr - (2*sin(xr))/(xr*xr*xr) + (2*cos(xr))/(xr*xr));
//                 for (int jk=2; jk<(l+1); jk++) {
//                     fa[jk] = ((2*(jk+1)-3)/xr)*fa[jk-1] - fa[jk-2];                 
//                     dfa[jk] = ((2*(jk+1)-3)/xr)*dfa[jk-1] - x0[k*(L+1)+l]*((2*(jk+1)-3)/(xr*xr))*fa[jk-1] - dfa[jk-2];        
//                 }
//                 tmpa[k] = fa[l];
//                 dtmpa[k] = dfa[l];                    
//             }

            // Compute spherical harmonics Bessel functions and their derivatives at level l for m = 0, 1, ..., l
            for (m=0; m<=l; m++) {
                // Spherical harmonics Y_{lm}(x,y,z) and their derivatives for l and m = 0,1,..,l  
                C = sqrt((2*l+1)*fac[l-m]/(4*pi*fac[l+m])); 
                Ylmr = C*(cos(m*phi)*Pa[m]);
                Ylmi = C*(sin(m*phi)*Pa[m]);
                YlmrThe = C*(cos(m*phi)*dPa[m]); 
                YlmiThe = C*(sin(m*phi)*dPa[m]); 
                YlmrPhi = -(m*C)*(sin(m*phi)*Pa[m]); 
                YlmiPhi = (m*C)*(cos(m*phi)*Pa[m]);          
                Ylmrx = YlmrThe*Thex + YlmrPhi*Phix;
                Ylmry = YlmrThe*They + YlmrPhi*Phiy;
                Ylmrz = YlmrThe*Thez + YlmrPhi*Phiz;
                Ylmix = YlmiThe*Thex + YlmiPhi*Phix;
                Ylmiy = YlmiThe*They + YlmiPhi*Phiy;
                Ylmiz = YlmiThe*Thez + YlmiPhi*Phiz; 

                for (k=0; k<K; k++) {
                    // Spherical Bessel g_{lk}(x,y,z) and their derivatives for l >=2
                    xr = x0[k*(L+1)+l]*r;                                            
                    fa[0] = sin(xr)/xr;
                    fa[1] = sin(xr)/(xr*xr) - cos(xr)/xr;                    
                    dfa[0] = x0[k*(L+1)+l]*(cos(xr)/xr - sin(xr)/(xr*xr));
                    dfa[1] = x0[k*(L+1)+l]*(sin(xr)/xr - (2*sin(xr))/(xr*xr*xr) + (2*cos(xr))/(xr*xr));
                    for (int jk=2; jk<(l+1); jk++) {
                        fa[jk] = ((2*(jk+1)-3)/xr)*fa[jk-1] - fa[jk-2];                 
                        dfa[jk] = ((2*(jk+1)-3)/xr)*dfa[jk-1] - x0[k*(L+1)+l]*((2*(jk+1)-3)/(xr*xr))*fa[jk-1] - dfa[jk-2];        
                    }
                    g = fa[l];
                    gR = dfa[l];                                        
                    //g = tmpa[k];
                    //gR = dtmpa[k];                    
                    gx = gR*Rx;
                    gy = gR*Ry;
                    gz = gR*Rz; 
                    
                    // derivatives of S_{klm}(x,y,z) = g_{lk}(x,y,z)*Y_{lm}(x,y,z) for l >= 2 and m = 0,1,...,l
                    j = ((l*l + l + m)*K + k)*N + i;      
                    Sr[j] = g*Ylmr; // real part            
                    Si[j] = g*Ylmi; // imag part                                                                                
                    Srx[j] = gx*Ylmr + g*Ylmrx;
                    Sry[j] = gy*Ylmr + g*Ylmry;
                    Srz[j] = gz*Ylmr + g*Ylmrz;
                    Six[j] = gx*Ylmi + g*Ylmix;
                    Siy[j] = gy*Ylmi + g*Ylmiy;
                    Siz[j] = gz*Ylmi + g*Ylmiz;        
                }                
            }
            
            // derivatives of S_{klm}(x,y,z) = g_{lk}(x,y,z)*Y_{lm}(x,y,z) for l >= 2 and m = -1,-2,...,-l            
            C = -1.0;
            for (int m=0; m<l; m++) {  
                for (k=0; k<K; k++) {
                    j  = ((l*l + l + m + 1)*K + k)*N + i;  // spherical harmonics Bessel functions for m > 0 
                    jm = ((l*l + l - m - 1)*K + k)*N + i;  // spherical harmonics Bessel functions for m < 0                                     
                    Sr[jm] = C*Sr[j]; // real part
                    Si[jm] =-C*Si[j]; // imag part                                                                                                                             
                    Srx[jm] = C*Srx[j]; 
                    Six[jm] =-C*Six[j];
                    Sry[jm] = C*Sry[j]; 
                    Siy[jm] =-C*Siy[j];
                    Srz[jm] = C*Srz[j]; 
                    Siz[jm] =-C*Siz[j];                    
                }                      
                C = C*(-1.0);        
            }
        }
        i += blockDim.x * gridDim.x;
    }                        
}
template <typename T> void gpuSphericalHarmonicsBesselJWithDeriv(T *Sr, T *Si,
        T *Srx, T *Six, T *Sry, T *Siy, T *Srz, T *Siz, T *xij, 
        T *x0, T *P, T *tmp, T *f, T *dP, T *dtmp, T *df, T *fac, T pi, int L, int K, int N)
{        
    int blockDim = 256;
    int gridDim = (N + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;
    gpuKernelSphericalHarmonicsBesselJWithDeriv<<<gridDim, blockDim>>>(Sr, Si, Srx, Six, Sry, Siy, Srz, Siz, xij,
                 x0, fac, pi, L, K, N);
}
template void gpuSphericalHarmonicsBesselJWithDeriv(double*, double*, double*, double*, double*, 
        double*, double*, double*, double*, double*, double*, double*, double*, double*, 
        double*, double*, double*, double, int, int, int);
template void gpuSphericalHarmonicsBesselJWithDeriv(float*, float*, float*, float*, float*, 
        float*, float*, float*, float*, float*, float*, float*, float*, 
        float*, float*, float*, float*, float, int, int, int);

#endif


