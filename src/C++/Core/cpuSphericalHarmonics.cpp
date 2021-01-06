#ifndef __CPUSPHERICALHARMONICS
#define __CPUSPHERICALHARMONICS

// template <typename T>void printarr(T* a, Int m, Int n)
// {    
//     for (Int i=0; i<m; i++) {
//         for (Int j=0; j<n; j++)
//             cout << scientific << a[j*m+i] << "   ";
//         cout << endl;
//     }
//     cout << endl;
// }
// template void printarr(double*, int, int);
// template void printarr(float*, int, int);

// core function
void cpuGetIndk(int *indk, int K)
{               
    // K    : the number (degree) of spherical Bessel functions
    // indk : store indices for symmetric tensor products of spherical Bessel functions
    // compute indices for symmetric tensor products of spherical Bessel functions
    
    int N = K*(K+1)/2;
    int count = 0;
    for (int k1=0; k1<K; k1++)
        for (int k2=0; k2<=k1; k2++) { // k2 <= k1
            indk[count] = k2;
            indk[N+count] = k1;
            count += 1;
        }
}

// core function
template <typename T> T clebschgordan(int j1, int m1, int j2, int m2, int j, int m, T *fac)
{               
    // compute the Clebsch-Gordan coefficient for two tuple pairs (j2, j1, j) and (m2, m1, m)
    
    T C;
    
    int jm = max(j1-j2,j2-j1);
    
    if ((m1+m2-m != 0) || (j < jm)  || (j > j1+j2)) { // zero Clebsch-Gordan coefficient 
        C = (T) 0.0;
    }
    else { // non-zero Clebsch-Gordan coefficient 
        int k1 = max(max(0,j2-j-m1),j1-j+m2);
        int k2 = min(min(j1+j2-j,j1-m1),j2+m2);
        int n = k2-k1+1;        
        C = 0.0;
        for (int i=0; i<n; i++) {
            int k = k1 + i;        
            T a = ( k % 2 == 0) ? (T) 1.0 : (T) -1.0;
            C += a/(fac[k]*fac[j1+j2-j-k]*fac[j1-m1-k]*fac[j2+m2-k]*fac[j-j2+m1+k]*fac[j-j1-m2+k]);
        }
        C = C*sqrt((2*j+1)*fac[j+j1-j2]*fac[j+j2-j1]*fac[j1+j2-j]/fac[j1+j2+j+1])*
              sqrt(fac[j+m]*fac[j-m]*fac[j1+m1]*fac[j1-m1]*fac[j2+m2]*fac[j2-m2]);        
    }    
    return C;
}
template double clebschgordan(int, int, int, int, int, int, double*);
template float clebschgordan(int, int, int, int, int, int, float*);

// core function
int cgcoefficients(int *indl, int N)
{                   
    // N          : the number of tuples (l2, l1, l)
    // indl [N*3] : store the indices of tuple (l2,l1,l)
    // compute the total number of non-zero Clebsch-Gordan coefficients for N tuples (l2, l1, l)
    
    int M = 0;
    for (int n=0; n<N; n++) { // loop over each tuple (l2, l1, l)
        int l2 = indl[n];
        int l1 = indl[N+n];
        int l = indl[2*N+n];
        int lm = max(l1-l2,l2-l1);
        // loop over all possible tuples (m2, m1, m)
        for (int m=-l; m<=l; m++) 
            for (int m1=-l1; m1<=l1; m1++) 
                for (int m2=-l2; m2<=l2; m2++) 
                    if (!((m1+m2-m != 0) || (l < lm) || (l > l1+l2)))                        
                        M += 1; // increment M since the Clebsch-Gordan coefficient is non-zero                     
    }    
    return M;
}

// core function
template <typename T> void cgcoefficients(T *cg, int *indm, int *rowm, int *indl, T *fac, int M, int N)
{                   
    // M             : the total number of non-zero Clebsch-Gordan coefficients computed in the previous functions
    // N             : the number of tuples (l2, l1, l)
    // cg   [M*1]    : store non-zero Clebsch-Gordan coefficients
    // indm [M*3]    : store the indices of tuple (m2, m1, m) for non-zero Clebsch-Gordan coefficients
    // rowm [(N+1)*1]: store the number of indices of tuple (m2, m1, m) for each tuple (l2, l1, l)
    // indl [N*3]    : store the indices of tuple (l2,l1,l)
    // fac           : factorial look-up table
    // compute the total number of non-zero Clebsch-Gordan coefficients for N tuples (l2, l1, l)
    
    int nc = 0;
    rowm[0] = 0;
    for (int n=0; n<N; n++) { // loop over each tuple (l2, l1, l)
        int l2 = indl[n];
        int l1 = indl[N+n];
        int l = indl[2*N+n];
        int lm = max(l1-l2,l2-l1);
        int count = 0;
        // loop over all possible tuples (m2, m1, m)
        for (int m=-l; m<=l; m++) 
            for (int m1=-l1; m1<=l1; m1++) 
                for (int m2=-l2; m2<=l2; m2++) 
                    if (!((m1+m2-m != 0) || (l < lm) || (l > l1+l2))) {                    
                        T C = clebschgordan(l1, m1, l2, m2, l, m, fac);   
                        cg[nc] = C;
                        indm[nc] = m2;   
                        indm[M+nc] = m1;   
                        indm[2*M+nc] = m;   
                        nc += 1;    // increment nc since the Clebsch-Gordan coefficient is non-zero                     
                        count += 1; // increment count since the Clebsch-Gordan coefficient is non-zero                     
                    }
        rowm[n+1] = count;
    }
    
    for (int n=2; n<=N; n++)
        rowm[n] = rowm[n] + rowm[n-1];    
}
template void cgcoefficients(double*, int*, int*, int*, double*, int, int);
template void cgcoefficients(float*, int*, int*, int*, float*, int, int);

// core function
template <typename T> void cpuSphericalHarmonicsBispectrum(T *b, T *Sr, T *Si, T *fac, int L)
{                   
    // L                    : the maximum degree of spherical harmonics
    // Sr      [(L+1)*(L+1)]: real part of spherical harmonics
    // Si      [(L+1)*(L+1)]: imag part of spherical harmonics
    // fac                  : factorial look-up table
    // b [(L+1)*(L+1)*(L+1)]: store the bispectrum coeffcients for spherical harmonic functions
    // compute the bispectrum coeffcients for spherical harmonic functions
    
    int N = L + 1;
    for (int l=0; l<N; l++)     
        for (int l1=0; l1<N; l1++)     
            for (int l2=0; l2<N; l2++)  {                         
                T tmp = (T) 0.0;
                for (int m=-l; m<=l; m++) 
                    for (int m1=-l1; m1<=l1; m1++) 
                        for (int m2=-l2; m2<=l2; m2++) {
                            T cg = clebschgordan(l2,m2,l1,m1,l,m,fac);
                            if (fabs(cg)>0) {
                                T a1, b1, a2, b2, a3, b3;                                
                                a1 = Sr[l*l + l + m];
                                b1 = Si[l*l + l + m];
                                a2 = Sr[l1*l1 + l1 + m1];
                                b2 = Si[l1*l1 + l1 + m1];
                                a3 = Sr[l2*l2 + l2 + m2];
                                b3 = Si[l2*l2 + l2 + m2];
                                tmp = tmp + cg*(a1*a2*a3 + a2*b1*b3 + a3*b1*b2 - a1*b2*b3);         
                            }
                        }
                b[l*N*N+l1*N+l2] = tmp;
            }
}
template void cpuSphericalHarmonicsBispectrum(double*, double*, double*, double*, int);
template void cpuSphericalHarmonicsBispectrum(float*, float*, float*, float*, int);

// core function
template <typename T> int cpuSphericalHarmonicsBispectrumIndex(T *b, int L)
{               
    // L : the maximum degree of spherical harmonics    
    // b : the bispectrum coeffcients of spherical harmonic functions computed by the above function
    // compute the number of non-zero unqiue bispectrum coeffcients for spherical harmonic functions
    
    // Symmetry conditions for the bispectrum components
    // [l2 l1 l]
    // [l1 l2 l]
    // [l l2 l1]
    // [l2 l l1]
    // [l1 l l2]
    // [l l1 l2]
    // Use the symmetry conditions to remove zero and duplicate bispectrum compoments
    
    int N = L+1;
    int inc = 0;
    for (int l=0; l<N; l++)     
        for (int l1=0; l1<N; l1++)     
            for (int l2=0; l2<N; l2++)  {                         
                T tm = b[l*N*N+l1*N+l2];            
                if ((fabs(tm)>1e-10) && (l2<=l1)) {
                    if (l1==l2) {
                        inc = inc + 1;
                    }
                    else {
                        if ((l2<l1) && (l2<l) && (l1<l))
                            inc = inc + 1;                                                
                    }
                }
            }
    return inc;
}
template int cpuSphericalHarmonicsBispectrumIndex(double *, int);
template int cpuSphericalHarmonicsBispectrumIndex(float *, int);

// core function
template <typename T> void cpuSphericalHarmonicsBispectrumIndex(int *indl, T *b, int M, int L)
{    
    // L    : the maximum degree of spherical harmonics    
    // M    : the number of non-zero unqiue bispectrum coeffcients for spherical harmonic functions
    // b    : the bispectrum coeffcients of spherical harmonic functions computed by the above function
    // indl : store indices (l2,l1,l) of non-zero unqiue bispectrum coeffcients for spherical harmonic functions
    // compute the indices of non-zero unqiue bispectrum coeffcients for spherical harmonic functions
    
    int N = L+1;
    int inc = 0;
    for (int l=0; l<N; l++)     
        for (int l1=0; l1<N; l1++)     
            for (int l2=0; l2<N; l2++)  {                         
                T tm = b[l*N*N+l1*N+l2];            
                if ((fabs(tm)>1e-10) && (l2<=l1)) {
                    if (l1==l2) {                        
                        indl[inc] = l2;
                        indl[M+inc] = l1;
                        indl[2*M+inc] = l;
                        inc = inc + 1;
                    }
                    else {
                        if ((l2<l1) && (l2<l) && (l1<l)) {
                            indl[inc] = l2;
                            indl[M+inc] = l1;
                            indl[2*M+inc] = l;
                            inc = inc + 1;                
                        }
                    }
                }
            }            
}
template void cpuSphericalHarmonicsBispectrumIndex(int*, double *, int, int);
template void cpuSphericalHarmonicsBispectrumIndex(int*, float *, int, int);

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

// https://en.wikipedia.org/wiki/Spherical_harmonics
template <typename T> void cpuSphericalHarmonics(T *Ylmr, T *Ylmi, T *x, T *y, T *z, 
                T *P, T *tmp, T *fac, T pi, int L, int N)
{                           
    // L             : the maximum degree of spherical harmonics
    // N             : length of x, y, z
    // P    [(L+1)*(L+2)/2]  : temporary storage for the recurrence formula
    // tmp  [(L+1)*(L+2)/2]  : temporary storage for the recurrence formula
    // fac                   : factorial look-up table
    // Ylmr [N*(L+1)*(L+2)/2]: real part of spherical harmonics
    // Ylmi [N*(L+1)*(L+2)/2]: imag part of spherical harmonics
    
    // Compute spherical harmonics functions, Y_{lm}(x,y,z), for l = 0,1,...,L and m = -l,...,0,...,l        
    //  l = 0:                0
    //  l = 1:              1 2 3
    //  l = 2:            4 5 6 7 8
    //  l = 3:        9 10 11 12 13 14 15
    //  l = 4:    16 17 18 19 20 21 22 23 24
    //  ....
    //  l = L:  L*L L*L+1 ....   ....  (L+1)*(L+1)-1 
    // The total number of spherical harmonics is (L+1)*(L+1).    
    
    // Due to symmetry properties of the spherical harmonics, we actually compute the following spherical harmonics for m >=0 
    //  l = 0:                0
    //  l = 1:                2 3
    //  l = 2:                6 7 8
    //  l = 3:                12 13 14 15
    //  l = 4:                20 21 22 23 24
    //  ....
    //  l = L:                L*L+L ....  (L+1)*(L+1)-1     
    //  Hence, the total number of computed spherical harmonics is (L+1)*(L+2)/2.
    
    //#pragma omp parallel for
    for (int i=0; i<N; i++) {  // loop over each neighbor atom                             
        tmp[0] = 1.0;
        
        // Cartersian -> Spherical
        T r = sqrt(x[i]*x[i] + y[i]*y[i] + z[i]*z[i]);
        T the = acos(z[i]/r);
        T phi = atan2(y[i],x[i]);
                
        // Spherical harmonics for l = 0        
        Ylmr[i] = sqrt(1/(4*pi)); // real part
        Ylmi[i] = 0.0;            // imag part 
                
        // Spherical harmonics for l = 1
        int l = 1;
        int m = 0;    
        T costhe = cos(the);
        T a = -sin(the);        
        P[0] = costhe;
        P[1] = a;
        
        m = 0;
        T C = sqrt((2*l+1)*fac[l-m]/(4*pi*fac[l+m]));
        Ylmr[N+i] = C*P[0]; // real part
        Ylmi[N+i] = 0.0;    // imag part 
        
        m = 1;
        C = sqrt((2*l+1)*fac[l-m]/(4*pi*fac[l+m]));
        Ylmr[2*N+i] = C*cos(phi)*P[1]; // real part
        Ylmi[2*N+i] = C*sin(phi)*P[1]; // imag part                
        
        // Spherical harmonics for l > 1
        for (l=2; l<=L; l++) {                                        
            // Compute associated Legendre polynomial using recurrence formula
            // https://en.wikipedia.org/wiki/Associated_Legendre_polynomials
            T Pll = P[l-1];
            tmp[(l-1)] = Pll;
            P[(l-1)] = (2*(l-1)+1)*(costhe*Pll);
            P[l] = (2*(l-1)+1)*(a*Pll); 
            for (m=0; m<l-1; m++) {
                Pll = P[m];
                P[m] = ((2*(l-1)+1)*(costhe*Pll) - (l+m-1)*tmp[m])/(T (l-m));
                tmp[m] = Pll;
            }
                                         
            // compute spherical harmonics at level l for m = 0, 1, ..., l
            for (m=0; m<=l; m++) {
                C = sqrt((2*l+1)*fac[l-m]/(4*pi*fac[l+m])); 
                Ylmr[(l*(l+1)/2+m)*N + i] = C*(cos(m*phi)*P[m]); // real part
                Ylmi[(l*(l+1)/2+m)*N + i] = C*(sin(m*phi)*P[m]); // imag part 
            }
        }
    }                        
}
template void cpuSphericalHarmonics(double*, double*, double*, double*, double*, double*, double*, double*, double, int, int);
template void cpuSphericalHarmonics(float*, float*, float*, float*, float*, float*, float*, float*, float, int, int);

// core function
template <typename T> void cpuSphericalHarmonicsBessel(T *Sr, T *Si, T *x, T *y, T *z, 
                T *x0, T *P, T *tmp, T *f, T *fac, T pi, int L, int K, int N)
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
    for (int i=0; i<N; i++) { // loop over each neighbor atom               
        int k, j, jm;
        tmp[0] = 1.0;
        
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
        P[0] = costhe;
        P[1] = a;
        
        T C = sqrt((2*l+1)*fac[l-m]/(4*pi*fac[l+m]));
        Ylmr = C*P[0];
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
        Ylmr = C*cos(phi)*P[1];
        Ylmi = C*sin(phi)*P[1];                
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
            T Pll = P[l-1];
            tmp[(l-1)] = Pll;
            P[(l-1)] = (2*(l-1)+1)*(costhe*Pll);
            P[l] = (2*(l-1)+1)*(a*Pll);             
            for (m=0; m<l-1; m++) {
                Pll = P[m];
                P[m] = ((2*(l-1)+1)*(costhe*Pll) - (l+m-1)*tmp[m])/(T (l-m));
                tmp[m] = Pll;                
            }
                            
            // Compute spherical harmonics Bessel function at level l for m = 0, 1, ..., l
            for (m=0; m<=l; m++) {
                C = sqrt((2*l+1)*fac[l-m]/(4*pi*fac[l+m])); 
                Ylmr = C*(cos(m*phi)*P[m]);
                Ylmi = C*(sin(m*phi)*P[m]);                
                for (k=0; k<K; k++) {
                    // Compute the spherical Bessel functions using recurrence formula
                    xr = x0[k*(L+1)+l]*r;                                            
                    f[0] = -cos(xr)/xr;
                    f[1] = -cos(xr)/(xr*xr) - sin(xr)/xr;                    
                    for (int jk=2; jk<(l+1); jk++) 
                        f[jk] = ((2*(jk+1)-3)/xr)*f[jk-1] - f[jk-2];                                     
                    g = -f[l];                                                              
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
    }          
}
template void cpuSphericalHarmonicsBessel(double*, double*, double*, double*, double*, double*, 
        double*, double*, double*, double*, double, int, int, int);
template void cpuSphericalHarmonicsBessel(float*, float*, float*, float*, float*, float*, 
        float*, float*, float*, float*, float, int, int, int);

// core function
template <typename T> void cpuSphericalHarmonicsSum(T *Ylmr, T *Ylmi, T *x, T *y, T *z, 
                T *P, T *tmp, T *fac, T pi, int L, int N)
{                            
    // L                   : the maximum degree of spherical harmonics
    // N                   : length of x, y, z
    // P   [(L+1)*(L+2)/2] : temporary storage for the recurrence formula
    // tmp [(L+1)*(L+2)/2] : temporary storage for the recurrence formula
    // fac                 : factorial look-up table
    // Ylmr  [(L+1)*(L+1)] : real part of spherical harmonics functions
    // Ylmi  [(L+1)*(L+1)] : imag part of spherical harmonics functions
    
    // Compute spherical harmonics functions, sum_i Y_{lm}(x_i,y_i,z_i), for l = 0,1,...,L and m = -l,...,0,...,l        
    //  l = 0:                0
    //  l = 1:              1 2 3
    //  l = 2:            4 5 6 7 8
    //  l = 3:        9 10 11 12 13 14 15
    //  l = 4:    16 17 18 19 20 21 22 23 24
    //  ....
    //  l = L:  L*L L*L+1 ....   ....  (L+1)*(L+1)-1 
    // The total number of spherical harmonics is (L+1)*(L+1).        
    for (int l=0; l<(L+1)*(L+1); l++)  {
        Ylmr[l] = 0.0;
        Ylmi[l] = 0.0;
    }

    //#pragma omp parallel for
    for (int i=0; i<N; i++) {  // loop over each neighbor atom                             
        tmp[0] = 1.0;
        
        // Cartersian -> Spherical
        T r = sqrt(x[i]*x[i] + y[i]*y[i] + z[i]*z[i]);
        T the = acos(z[i]/r);
        T phi = atan2(y[i],x[i]);
                
        // Spherical harmonics for l = 0        
        Ylmr[0] += sqrt(1/(4*pi)); // real part
        Ylmi[0] += 0.0;            // imag part 
                
        // Spherical harmonics for l = 1
        int l = 1;
        int m = 0;    
        T costhe = cos(the);
        T a = -sin(the);        
        P[0] = costhe;
        P[1] = a;
        
        m = 0;
        T C = sqrt((2*l+1)*fac[l-m]/(4*pi*fac[l+m]));
        int j = (l*l + l + m);
        Ylmr[j] += C*P[0]; // real part
        Ylmi[j] += 0.0;    // imag part 
        
        m = 1;
        j = (l*l + l + m);
        C = sqrt((2*l+1)*fac[l-m]/(4*pi*fac[l+m]));
        Ylmr[j] += C*cos(phi)*P[1]; // real part
        Ylmi[j] += C*sin(phi)*P[1]; // imag part                
        
        // Spherical harmonics for l > 1
        for (l=2; l<=L; l++) {                                        
            // Compute associated Legendre polynomial using recurrence formula
            // https://en.wikipedia.org/wiki/Associated_Legendre_polynomials
            T Pll = P[l-1];
            tmp[(l-1)] = Pll;
            P[(l-1)] = (2*(l-1)+1)*(costhe*Pll);
            P[l] = (2*(l-1)+1)*(a*Pll); 
            for (m=0; m<l-1; m++) {
                Pll = P[m];
                P[m] = ((2*(l-1)+1)*(costhe*Pll) - (l+m-1)*tmp[m])/(T (l-m));
                tmp[m] = Pll;
            }
                                         
            // compute spherical harmonics at level l for m = 0, 1, ..., l
            for (m=0; m<=l; m++) {
                j = (l*l + l + m);
                C = sqrt((2*l+1)*fac[l-m]/(4*pi*fac[l+m])); 
                Ylmr[j] += C*(cos(m*phi)*P[m]); // real part
                Ylmi[j] += C*(sin(m*phi)*P[m]); // imag part 
            }
        }
    }
    
    // Compute the spherical harmonics for m < 0 using symmetry properties
    for (int l=1; l<=L; l++) {
        T C = -1.0;
        for (int m=1; m<=l; m++)  {                 
            int j =  (l*l + l + m); // spherical harmonics for m > 0 
            int jm = (l*l + l - m); // spherical harmonics for m < 0 
            Ylmr[jm] = C*Ylmr[j]; // real part
            Ylmi[jm] =-C*Ylmi[j]; // imag part                                       
            C = C*(-1.0);        
        }
    }        
}
template void cpuSphericalHarmonicsSum(double*, double*, double*, double*, double*, double*, double*, double*, double, int, int);
template void cpuSphericalHarmonicsSum(float*, float*, float*, float*, float*, float*, float*, float*, float, int, int);

template <typename T> void cpuSphericalHarmonicsSum(T *Sr, T *Si, T *Ylmr, T *Ylmi, int L, int N)
{                        
    for (int l=0; l<(L+1); l++) {
        // Sum the spherical harmonics for m >= 0 
        for (int m=0; m<(l+1); m++) {
            int Nl = (l*(l+1)/2 + m)*N;
            int j = l*l + l + m;
            Sr[j] = 0.0; 
            Si[j] = 0.0;
            for (int i=0; i<N; i++) {                   
                Sr[j] += Ylmr[Nl + i];
                Si[j] += Ylmi[Nl + i];                                
            }
        }
        // Sum the spherical harmonics for m < 0 using symmetry properties
        T C = -1.0;
        for (int m=0; m<l; m++) {            
            int j = l*l + l + m + 1;
            int k = l*l + l - m - 1;
            Sr[k] = C*Sr[j]; 
            Si[k] =-C*Si[j];
            C = C*(-1.0);        
        }
    }    
}
template void cpuSphericalHarmonicsSum(double*, double*, double*, double*, int, int);
template void cpuSphericalHarmonicsSum(float*, float*, float*, float*, int, int);

template <typename T> void cpuRadialSphericalHarmonicsSum(T *Sr, T *Si, T *Ylmr, T *Ylmi, T *g, int L, int K, int N)
{                        
    for (int k=0; k<K; k++) 
        for (int l=0; l<(L+1); l++) {
            // Sum the radial spherical harmonics for m >= 0 
            for (int m=0; m<(l+1); m++) {
                int Nl = (l*(l+1)/2 + m)*N;
                int Nk = l*N*K + k*N;
                int j = (l*l + l + m)*K + k;                
                Sr[j] = 0.0; 
                Si[j] = 0.0;
                for (int i=0; i<N; i++) {                   
                    Sr[j] += g[Nk + i]*Ylmr[Nl + i];
                    Si[j] += g[Nk + i]*Ylmi[Nl + i];                                
                }
            }
            // Sum the radial spherical harmonics for m < 0 using symmetry properties
            T C = -1.0;
            for (int m=0; m<l; m++) {            
                int j = (l*l + l + m + 1)*K + k;
                int n = (l*l + l - m - 1)*K + k;
                Sr[n] = C*Sr[j]; 
                Si[n] =-C*Si[j];
                C = C*(-1.0);        
            }
        }    
}
template void cpuRadialSphericalHarmonicsSum(double*, double*, double*, double*, double*, int, int, int);
template void cpuRadialSphericalHarmonicsSum(float*, float*, float*, float*, float*, int, int, int);

template <typename T> void cpuSphericalHarmonicsBesselSum(T *Sr, T *Si, T *x, T *y, T *z, 
                T *x0, T *P, T *tmp, T *f, T *fac, T pi, int L, int K, int N)
{                        
    // L                   : the maximum degree of spherical harmonics
    // K                   : the maximum degree of spherical Bessel functions
    // N                   : length of x, y, z
    // x0  [(L+1)*K]       : zeros of pherical Bessel functions
    // P   [(L+1)*(L+2)/2] : temporary storage for the recurrence formula
    // tmp [(L+1)*(L+2)/2] : temporary storage for the recurrence formula
    // f   [L+1]           : temporary storage for the recurrence formula
    // fac                 : factorial look-up table
    // Sr  [K*(L+1)*(L+1)] : real part of spherical harmonics Bessel functions
    // Si  [K*(L+1)*(L+1)] : imag part of spherical harmonics Bessel functions
    
    // Compute spherical harmonics Bessel functions 
    // sum_i S_{klm}(x_i,y_i,z_i) = sum_i g_{lk}(x_i,y_i,z_i)*Y_{lm}(x_i,y_i,z_i)
    
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
    for (int i=0; i<K*(L+1)*(L+1); i++)  {
        Sr[i] = 0.0;
        Si[i] = 0.0;
    }
        
    //printarr(x0, L+1, K);    
    
    //#pragma omp parallel for    
    for (int i=0; i<N; i++) { // loop over each neighbor atom               
        int k, j, jm;
        tmp[0] = 1.0;
        
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
            j = k;                
            Sr[j] += g*Ylmr;  // real part          
            Si[j] += g*Ylmi;  // imag part                 
        }
                
        // Spherical harmonics for l = 1;
        l = 1;
        T costhe = cos(the);
        T a = -sin(the);        
        int m = 0;    
        P[0] = costhe;
        P[1] = a;
        
        T C = sqrt((2*l+1)*fac[l-m]/(4*pi*fac[l+m]));
        Ylmr = C*P[0];
        Ylmi = 0.0;
        // Spherical Bessel for l = 1, see https://en.wikipedia.org/wiki/Bessel_function#Spherical_Bessel_functions:_jn,_yn
        for (k=0; k<K; k++) {
            xr = x0[k*(L+1)+l]*r;                        
            g = cos(xr)/(xr*xr) + sin(xr)/xr;
            j = (l*l + l + m)*K + k;                
            Sr[j] += g*Ylmr; // real part           
            Si[j] += g*Ylmi; // imag part                                    
        }        
        
        m = 1;
        C = sqrt((2*l+1)*fac[l-m]/(4*pi*fac[l+m]));
        Ylmr = C*cos(phi)*P[1];
        Ylmi = C*sin(phi)*P[1];                
        // Spherical Bessel for l = 1, see https://en.wikipedia.org/wiki/Bessel_function#Spherical_Bessel_functions:_jn,_yn
        for (k=0; k<K; k++) {
            xr = x0[k*(L+1)+l]*r;                        
            g = cos(xr)/(xr*xr) + sin(xr)/xr; 
            j = (l*l + l + m)*K + k;                
            Sr[j] += g*Ylmr; // real part           
            Si[j] += g*Ylmi; // imag part                                                            
//             jm = (l*l + l - m)*K + k;                
//             Sr[jm] = -Sr[j];            
//             Si[jm] =  Si[j];            
        }        
                
        for (l=2; l<=L; l++) {                                        
            // Compute associated Legendre polynomial using recurrence formula
            // https://en.wikipedia.org/wiki/Associated_Legendre_polynomials            
            T Pll = P[l-1];
            tmp[(l-1)] = Pll;
            P[(l-1)] = (2*(l-1)+1)*(costhe*Pll);
            P[l] = (2*(l-1)+1)*(a*Pll);             
            for (m=0; m<l-1; m++) {
                Pll = P[m];
                P[m] = ((2*(l-1)+1)*(costhe*Pll) - (l+m-1)*tmp[m])/(T (l-m));
                tmp[m] = Pll;                
            }
                            
            // Compute spherical harmonics Bessel function at level l for m = 0, 1, ..., l
            for (m=0; m<=l; m++) {
                C = sqrt((2*l+1)*fac[l-m]/(4*pi*fac[l+m])); 
                Ylmr = C*(cos(m*phi)*P[m]);
                Ylmi = C*(sin(m*phi)*P[m]);                
                for (k=0; k<K; k++) {
                    // Compute the spherical Bessel functions using recurrence formula
                    xr = x0[k*(L+1)+l]*r;                                            
                    f[0] = -cos(xr)/xr;
                    f[1] = -cos(xr)/(xr*xr) - sin(xr)/xr;                    
                    for (int jk=2; jk<(l+1); jk++) 
                        f[jk] = ((2*(jk+1)-3)/xr)*f[jk-1] - f[jk-2];                                     
                    g = -f[l];                                                              
                    j = (l*l + l + m)*K + k;   
                    Sr[j] += g*Ylmr; // real part            
                    Si[j] += g*Ylmi; // imag part                                                            
                }                
            }                        
//             C = -1.0;
//             for (m=1; m<=l; m++)  {                 
//                 for (k=0; k<K; k++) {
//                     j =  (l*l + l + m)*K + k;  
//                     jm = (l*l + l - m)*K + k;                  
//                     Sr[jm] = C*Sr[j]; 
//                     Si[jm] =-C*Si[j];                                        
//                 }                     
//                 C = C*(-1.0);        
//             }
        }
    }
    
    // Compute the spherical harmonics Bessel functions for m < 0 using symmetry properties
    for (int l=1; l<=L; l++) {
        T C = -1.0;
        for (int m=1; m<=l; m++)  {                 
            for (int k=0; k<K; k++) {
                int j =  (l*l + l + m)*K + k; // spherical harmonics Bessel functions for m > 0 
                int jm = (l*l + l - m)*K + k; // spherical harmonics Bessel functions for m < 0 
                Sr[jm] = C*Sr[j]; // real part
                Si[jm] =-C*Si[j]; // imag part                                       
            }                     
            C = C*(-1.0);        
        }
    }    
}
template void cpuSphericalHarmonicsBesselSum(double*, double*, double*, double*, double*, double*, 
        double*, double*, double*, double*, double, int, int, int);
template void cpuSphericalHarmonicsBesselSum(float*, float*, float*, float*, float*, float*, 
        float*, float*, float*, float*, float, int, int, int);

// template <typename T> void cpuSphericalHarmonicsBesselSum(T *Sr, T *Si, T *Sklmr, T *Sklmi, int L, int K, int N)
// {                            
//     for (int l=0; l<(L+1); l++) {
//         // Sum the spherical harmonics Bessel functions
//         for (int m=0; m<(2*l+1); m++) {
//             for (int k=0; k<K; k++) {                
//                 int j = (l*l + m)*K + k;
//                 int Nj = j*N;                
//                 Sr[j] = 0.0; 
//                 Si[j] = 0.0;
//                 for (int i=0; i<N; i++) {                   
//                     Sr[j] += Sklmr[Nj + i];
//                     Si[j] += Sklmi[Nj + i];                                
//                 }
//             }
//         }
//     }    
// }
// template void cpuSphericalHarmonicsBesselSum(double*, double*, double*, double*, int, int, int);
// template void cpuSphericalHarmonicsBesselSum(float*, float*, float*, float*, int, int, int);

// core function
template <typename T> void cpuSphericalHarmonicsBesselDeriv(T *Srx, T *Six, T *Sry, T *Siy, T *Srz, T *Siz, T *x, T *y, T *z, 
                T *x0, T *P, T *tmp, T *f, T *dP, T *dtmp, T *df, T *fac, T pi, int L, int K, int N)
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
    for (int i=0; i<N; i++) { // loop over each neighbor atom                               
        int k, j, jm;
        
        // Cartersian -> Spherical
        T r = sqrt(x[i]*x[i] + y[i]*y[i] + z[i]*z[i]);
        T the = acos(z[i]/r);
        T phi = atan2(y[i],x[i]);
                
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
            Srx[j] = gx*Ylmr + g*Ylmrx;
            Sry[j] = gy*Ylmr + g*Ylmry;
            Srz[j] = gz*Ylmr + g*Ylmrz;
            Six[j] = gx*Ylmi + g*Ylmix;
            Siy[j] = gy*Ylmi + g*Ylmiy;
            Siz[j] = gz*Ylmi + g*Ylmiz;            
        }
                        
        l = 1;
        T costhe = cos(the);
        T a = -sin(the);        
        T dcosthe = a;
        T da = -costhe;        
        int m = 0;    
        P[0] = costhe;
        P[1] = a;
        dP[0] = dcosthe;    
        dP[1] = da;
        tmp[0] = 1.0;
        dtmp[0] = 0.0;
        
        // Spherical harmonics Y_{lm}(x,y,z) and their derivatives for l = 1 and m = 0  
        m = 0;
        T C = sqrt((2*l+1)*fac[l-m]/(4*pi*fac[l+m]));
        Ylmr = C*P[0];
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
        Ylmr = C*cos(phi)*P[1];
        Ylmi = C*sin(phi)*P[1];                
        YlmrThe = C*(cos(phi)*da);      
        YlmiThe = C*(sin(phi)*da);      
        YlmrPhi = -C*(sin(phi)*P[1]); 
        YlmiPhi = C*(cos(phi)*P[1]);  
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
            
            // derivatives of S_{klm}(x,y,z) = g_{lk}(x,y,z)*Y_{lm}(x,y,z) for l = 1 and m = 1
            j = ((l*l + l + m)*K + k)*N + i;                
            Srx[j] = gx*Ylmr + g*Ylmrx;
            Sry[j] = gy*Ylmr + g*Ylmry;
            Srz[j] = gz*Ylmr + g*Ylmrz;
            Six[j] = gx*Ylmi + g*Ylmix;
            Siy[j] = gy*Ylmi + g*Ylmiy;
            Siz[j] = gz*Ylmi + g*Ylmiz;        
            
            // derivatives of S_{klm}(x,y,z) = g_{lk}(x,y,z)*Y_{lm}(x,y,z) for l = 1 and m = -1
            jm = ((l*l + l - m)*K + k)*N + i;                
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
            T Pll = P[l-1];
            tmp[(l-1)] = Pll;
            P[(l-1)] = (2*(l-1)+1)*(costhe*Pll);
            P[l] = (2*(l-1)+1)*(a*Pll); 
            T dPll = dP[l-1];
            dtmp[l-1] = dPll;
            dP[l-1] = (2*(l-1)+1)*(costhe*dPll + dcosthe*Pll);   
            dP[l] = (2*(l-1)+1)*(a*dPll + da*Pll);             
            for (m=0; m<l-1; m++) {
                Pll = P[m];
                P[m] = ((2*(l-1)+1)*(costhe*Pll) - (l+m-1)*tmp[m])/(T (l-m));
                tmp[m] = Pll;
                dPll = dP[m];
                dP[m] = ((2*(l-1)+1)*(costhe*dPll+dcosthe*Pll) - (l+m-1)*dtmp[m])/(T (l-m));
                dtmp[m] = dPll;
            }
                            
            // Compute spherical harmonics Bessel functions and their derivatives at level l for m = 0, 1, ..., l
            for (m=0; m<=l; m++) {
                // Spherical harmonics Y_{lm}(x,y,z) and their derivatives for l and m = 0,1,..,l  
                C = sqrt((2*l+1)*fac[l-m]/(4*pi*fac[l+m])); 
                Ylmr = C*(cos(m*phi)*P[m]);
                Ylmi = C*(sin(m*phi)*P[m]);
                YlmrThe = C*(cos(m*phi)*dP[m]); 
                YlmiThe = C*(sin(m*phi)*dP[m]); 
                YlmrPhi = -(m*C)*(sin(m*phi)*P[m]); 
                YlmiPhi = (m*C)*(cos(m*phi)*P[m]);          
                Ylmrx = YlmrThe*Thex + YlmrPhi*Phix;
                Ylmry = YlmrThe*They + YlmrPhi*Phiy;
                Ylmrz = YlmrThe*Thez + YlmrPhi*Phiz;
                Ylmix = YlmiThe*Thex + YlmiPhi*Phix;
                Ylmiy = YlmiThe*They + YlmiPhi*Phiy;
                Ylmiz = YlmiThe*Thez + YlmiPhi*Phiz; 

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
                    
                    // derivatives of S_{klm}(x,y,z) = g_{lk}(x,y,z)*Y_{lm}(x,y,z) for l >= 2 and m = 0,1,...,l
                    j = ((l*l + l + m)*K + k)*N + i;                
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
    }                        
}
template void cpuSphericalHarmonicsBesselDeriv(double*, double*, double*, double*, double*, double*, double*, double*, double*, 
        double*, double*, double*, double*, double*, double*, double*, double*, double, int, int, int);
template void cpuSphericalHarmonicsBesselDeriv(float*, float*, float*, float*, float*, float*, float*, float*, float*, 
        float*, float*, float*, float*, float*, float*, float*, float*, float, int, int, int);


template <typename T> void cpuRadialSphericalHarmonicsPower(T *p, T *ar, T *ai, int *indk, int L, int K)
{   
    // L                   : the maximum degree of spherical harmonics
    // K                   : the maximum degree of spherical Bessel functions
    // indk [K*(K+1)/2]    : store indices for symmetric tensor products of spherical Bessel functions
    // ar  [K*(L+1)*(L+1)] : sum of real spherical harmonics Bessel functions
    // ai  [K*(L+1)*(L+1)] : sum of imag spherical harmonics Bessel functions
    // p [(L+1)*K*(K+1)/2] : power spectrum components
    
    // Compute the power spectrum components for radial spherical harmonics
    
    int K2 = K*(K+1)/2;    
    //#pragma omp parallel for    
    for (int k=0; k<K2; k++) {
        int j = indk[k];
        int i = indk[K2+k];
        for (int l=0; l<(L+1); l++) {        
            int l1 = k*(L+1) + l;
            p[l1] = (T) 0.0;
            for (int m=0; m<(2*l+1); m++) {
                int i1 = (l*l + m)*K + i;
                int j1 = (l*l + m)*K + j;
                p[l1] += ar[i1]*ar[j1] + ai[i1]*ai[j1];            
            }
        }    
    }                
}
template void cpuRadialSphericalHarmonicsPower(double*, double*, double*, int*, int, int);
template void cpuRadialSphericalHarmonicsPower(float*, float*, float*, int*, int, int);

template <typename T> void cpuRadialSphericalHarmonicsPowerDeriv(T *px, T *py, T *pz, T *ar, T *ai, 
        T *arx, T *aix, T *ary, T *aiy, T *arz, T *aiz, int *indk, int L, int K, int N)
{   
    // L                      : the maximum degree of spherical harmonics
    // K                      : the maximum degree of spherical Bessel functions    
    // N                      : length of x, y, z
    // indk [K*(K+1)/2]       : store indices for symmetric tensor products of spherical Bessel functions
    // ar   [K*(L+1)*(L+1)]   : sum of real spherical harmonics Bessel functions
    // ai   [K*(L+1)*(L+1)]   : sum of imag spherical harmonics Bessel functions    
    // arx  [N*K*(L+1)*(L+1)] : x-derivative of real spherical harmonics Bessel functions
    // aix  [N*K*(L+1)*(L+1)] : x-derivative of imag spherical harmonics Bessel functions
    // ary  [N*K*(L+1)*(L+1)] : y-derivative of real spherical harmonics Bessel functions
    // aiy  [N*K*(L+1)*(L+1)] : y-derivative of imag spherical harmonics Bessel functions
    // arz  [N*K*(L+1)*(L+1)] : z-derivative of real spherical harmonics Bessel functions
    // aiz  [N*K*(L+1)*(L+1)] : z-derivative of imag spherical harmonics Bessel functions
    // p [(L+1)*K*(K+1)/2]    : power spectrum components
    // px [N*(L+1)*K*(K+1)/2] : x-derivative of power spectrum components
    // py [N*(L+1)*K*(K+1)/2] : y-derivative of power spectrum components
    // pz [N*(L+1)*K*(K+1)/2] : z-derivative of power spectrum components
    
    // Compute partial derivatives of the power spectrum components
    
    int K2 = K*(K+1)/2;    
    //#pragma omp parallel for    
    for (int n=0; n<N; n++)       
        for (int k=0; k<K2; k++) {
            int j = indk[k];
            int i = indk[K2+k];
            for (int l=0; l<(L+1); l++) {        
                int l1 = k*(L+1) + l;
                px[l1*N + n] = (T) 0.0;
                py[l1*N + n] = (T) 0.0;
                pz[l1*N + n] = (T) 0.0;
                for (int m=0; m<(2*l+1); m++) {
                    int i1 = (l*l + m)*K + i;
                    int j1 = (l*l + m)*K + j;
                    int i2 = (l*l + m)*N*K + i*N + n;      
                    int j2 = (l*l + m)*N*K + j*N + n;                                                 
                    px[l1*N + n] += ar[i1]*arx[j2] + arx[i2]*ar[j1] + ai[i1]*aix[j2] + aix[i2]*ai[j1];
                    py[l1*N + n] += ar[i1]*ary[j2] + ary[i2]*ar[j1] + ai[i1]*aiy[j2] + aiy[i2]*ai[j1];
                    pz[l1*N + n] += ar[i1]*arz[j2] + arz[i2]*ar[j1] + ai[i1]*aiz[j2] + aiz[i2]*ai[j1];                    
                }
            }    
        }                 
}
template void cpuRadialSphericalHarmonicsPowerDeriv(double*, double*, double*, double*, double*, double*, 
        double*, double*, double*, double*, double*, int*, int, int, int);
template void cpuRadialSphericalHarmonicsPowerDeriv(float*, float*, float*, float*, float*, float*,  
        float*, float*, float*, float*, float*, int*, int, int, int);

template <typename T> void cpuRadialSphericalHarmonicsBispectrum(T *b, T *ar, T *ai, T *cg, int *indk, int *indl,
        int *indm, int *rowm, int Nub, int Ncg, int K)
{   
    // K                   : the maximum degree of spherical Bessel functions
    // Ncg                 : the total number of non-zero Clebsch-Gordan coefficients 
    // Nub                 : the number of non-zero unique spherical harmonics bispectrum components
    // cg   [Ncg*1]        : store non-zero Clebsch-Gordan coefficients
    // indl [Nub*3]        : store the indices of tuple (l2,l1,l)        
    // indk [K*(K+1)/2]    : store indices for symmetric tensor products of spherical Bessel functions
    // indm [Ncg*3]        : store the indices of tuple (m2, m1, m) for non-zero Clebsch-Gordan coefficients
    // rowm [Nub+1]        : store the number of indices of tuple (m2, m1, m) for each tuple (l2, l1, l)    
    // ar  [K*(L+1)*(L+1)] : sum of real spherical harmonics Bessel functions
    // ai  [K*(L+1)*(L+1)] : sum of imag spherical harmonics Bessel functions
    // b [Nub*K*(K+1)/2]   : bispectrum components    
    
    // Compute the bispectrum components for radial spherical harmonics
    
    int K2 = K*(K+1)/2;    
    //#pragma omp parallel for                
    for (int k=0; k<K2; k++) {
        int k2 = indk[k];
        int k1 = indk[K2+k];
        for (int i=0; i<Nub; i++) {        
            int l2 = indl[i];
            int l1 = indl[Nub+i];
            int l = indl[2*Nub+i];     
            T tmp = (T) 0.0;
            int nm = rowm[i+1]-rowm[i];
            for (int j = 0; j<nm; j++) {
                int m2 = indm[rowm[i]+j];
                int m1 = indm[Ncg+rowm[i]+j];
                int m = indm[2*Ncg + rowm[i]+j];                                
                T a1, b1, a2, b2, a3, b3;                                
                a1 = ar[(l*l + l + m)*K + k1];
                b1 = ai[(l*l + l + m)*K + k1];
                a2 = ar[(l1*l1 + l1 + m1)*K + k2];
                b2 = ai[(l1*l1 + l1 + m1)*K + k2];
                a3 = ar[(l2*l2 + l2 + m2)*K + k2];
                b3 = ai[(l2*l2 + l2 + m2)*K + k2];
                tmp += cg[rowm[i]+j]*(a1*a2*a3 + a2*b1*b3 + a3*b1*b2 - a1*b2*b3);                                              
            }               
            b[i+k*Nub] = tmp;   
        }
    }        
}
template void cpuRadialSphericalHarmonicsBispectrum(double*, double*, double*, double*, int*, int*, int*, int*, int, int, int);
template void cpuRadialSphericalHarmonicsBispectrum(float*, float*, float*, float*, int*, int*, int*, int*, int, int, int);

template <typename T> void cpuRadialSphericalHarmonicsBispectrumDeriv(T *bx, T *by, T *bz, T *ar, T *ai, 
        T *arx, T *aix, T *ary, T *aiy, T *arz, T *aiz, T*cg, int *indk, int *indl,
        int *indm, int *rowm, int Nub, int Ncg, int K, int N)
{   
    // K                      : the maximum degree of spherical Bessel functions
    // N                      : length of x, y, z
    // Ncg                    : the total number of non-zero Clebsch-Gordan coefficients 
    // Nub                    : the number of non-zero unique spherical harmonics bispectrum components
    // cg   [Ncg*1]           : store non-zero Clebsch-Gordan coefficients
    // indl [Nub*3]           : store the indices of tuple (l2,l1,l)        
    // indk [K*(K+1)/2]       : store indices for symmetric tensor products of spherical Bessel functions
    // indm [Ncg*3]           : store the indices of tuple (m2, m1, m) for non-zero Clebsch-Gordan coefficients
    // rowm [Nub+1]           : store the number of indices of tuple (m2, m1, m) for each tuple (l2, l1, l)    
    // ar   [K*(L+1)*(L+1)]   : real part of spherical harmonics Bessel functions
    // ai   [K*(L+1)*(L+1)]   : imag part of spherical harmonics Bessel functions    
    // arx  [N*K*(L+1)*(L+1)] : x-derivative of real spherical harmonics Bessel functions
    // aix  [N*K*(L+1)*(L+1)] : x-derivative of imag spherical harmonics Bessel functions
    // ary  [N*K*(L+1)*(L+1)] : y-derivative of real spherical harmonics Bessel functions
    // aiy  [N*K*(L+1)*(L+1)] : y-derivative of imag spherical harmonics Bessel functions
    // arz  [N*K*(L+1)*(L+1)] : z-derivative of real spherical harmonics Bessel functions
    // aiz  [N*K*(L+1)*(L+1)] : z-derivative of imag spherical harmonics Bessel functions
    // b    [Nub*K*(K+1)/2]   : bispectrum components        
    // bx   [N*Nub*K*(K+1)/2] : x-derivative of power spectrum components
    // by   [N*Nub*K*(K+1)/2] : y-derivative of power spectrum components
    // bz   [N*Nub*K*(K+1)/2] : z-derivative of power spectrum components
        
    // Compute partial derivatives of the bispectrum components
    
    int K2 = K*(K+1)/2;    
    for (int n=0; n<N; n++)       
        for (int k=0; k<K2; k++) {
            int k2 = indk[k];
            int k1 = indk[K2+k];
            for (int i=0; i<Nub; i++) {        
                int l2 = indl[i];
                int l1 = indl[Nub+i];
                int l = indl[2*Nub+i];     
                int ii = (i+k*Nub)*N + n;
                bx[ii] = (T) 0.0;
                by[ii] = (T) 0.0;
                bz[ii] = (T) 0.0;                
                int nm = rowm[i+1]-rowm[i];
                for (int j = 0; j<nm; j++) {
                    int m2 = indm[rowm[i]+j];
                    int m1 = indm[Ncg+rowm[i]+j];
                    int m = indm[2*Ncg + rowm[i]+j];    
                    int lm1 = (l*l + l + m)*K + k1;
                    int lm2 = (l1*l1 + l1 + m1)*K + k2;
                    int lm3 = (l2*l2 + l2 + m2)*K + k2;
                    int mlk1 = lm1*N + n;      
                    int mlk2 = lm2*N + n;       
                    int mlk3 = lm3*N + n;                                             
                    
                    T c = cg[rowm[i]+j];                    
                    T a1, b1, a2, b2, a3, b3;                                
                    T a1x, b1x, a2x, b2x, a3x, b3x;
                    T a1y, b1y, a2y, b2y, a3y, b3y;
                    T a1z, b1z, a2z, b2z, a3z, b3z;
                    a1 = ar[lm1];
                    b1 = ai[lm1];
                    a2 = ar[lm2];
                    b2 = ai[lm2];
                    a3 = ar[lm3];
                    b3 = ai[lm3];
                    a1x = arx[mlk1];
                    a1y = ary[mlk1];
                    a1z = arz[mlk1];
                    b1x = aix[mlk1];
                    b1y = aiy[mlk1];
                    b1z = aiz[mlk1];
                    a2x = arx[mlk2];
                    a2y = ary[mlk2];
                    a2z = arz[mlk2];
                    b2x = aix[mlk2];
                    b2y = aiy[mlk2];
                    b2z = aiz[mlk2];
                    a3x = arx[mlk3];
                    a3y = ary[mlk3];
                    a3z = arz[mlk3];
                    b3x = aix[mlk3];
                    b3y = aiy[mlk3];
                    b3z = aiz[mlk3];

                    T t1 = a1x*a2*a3 + a1*a2x*a3 + a1*a2*a3x;                
                    T t2 = a2x*b1*b3 + a2*b1x*b3 + a2*b1*b3x;
                    T t3 = a3x*b1*b2 + a3*b1x*b2 + a3*b1*b2x;
                    T t4 = a1x*b2*b3 + a1*b2x*b3 + a1*b2*b3x;
                    bx[ii] += c*(t1 + t2 + t3 - t4);
                    
                    t1 = a1y*a2*a3 + a1*a2y*a3 + a1*a2*a3y;                
                    t2 = a2y*b1*b3 + a2*b1y*b3 + a2*b1*b3y;
                    t3 = a3y*b1*b2 + a3*b1y*b2 + a3*b1*b2y;
                    t4 = a1y*b2*b3 + a1*b2y*b3 + a1*b2*b3y;
                    by[ii] += c*(t1 + t2 + t3 - t4);
                    
                    t1 = a1z*a2*a3 + a1*a2z*a3 + a1*a2*a3z;                
                    t2 = a2z*b1*b3 + a2*b1z*b3 + a2*b1*b3z;
                    t3 = a3z*b1*b2 + a3*b1z*b2 + a3*b1*b2z;
                    t4 = a1z*b2*b3 + a1*b2z*b3 + a1*b2*b3z;
                    bz[ii] += c*(t1 + t2 + t3 - t4);                    
                }               
            }
        }                
}
template void cpuRadialSphericalHarmonicsBispectrumDeriv(double*, double*, double*, double*, double*, double*, 
        double*, double*, double*, double*, double*, double*, int*, int*, int*, int*, int, int, int, int);
template void cpuRadialSphericalHarmonicsBispectrumDeriv(float*, float*, float*, float*, float*, float*,  
        float*, float*, float*, float*, float*, float*, int*, int*, int*, int*, int, int, int, int);

// core function
template <typename T> void cpuRadialSphericalHarmonicsSum(T *ar, T *ai, T *Sr, T *Si, 
        int *Nnb, int Na, int L, int K)
{                            
    // L                         : the maximum degree of spherical harmonics
    // K                         : the maximum degree of spherical Bessel functions
    // Na                        : number of atoms in the simulation domain 
    // Nnb [Na+1]                : a list containing the number of neighbors for each global atom    
    // N = Nnb[Na]-Nnb[0]      : total number of neighbors
    // Sr [N*K*(L+1)*(L+1)]      : real spherical harmonics Bessel functions
    // Si [N*K*(L+1)*(L+1)]      : imag spherical harmonics Bessel functions        
    // ar [Na*K*(L+1)*(L+1)]     : sum of real spherical harmonics Bessel functions
    // ai [Na*K*(L+1)*(L+1)]     : sum of imag spherical harmonics Bessel functions    
        
    // Sum the spherical harmonics Bessel functions over neighbors
    int L2 = (L+1)*(L+1);
    int N = Nnb[Na]-Nnb[0]; // total number of neighbors 
    for (int l=0; l<L2; l++)
        for (int k=0; k<K; k++) 
            for (int n=0; n<Na; n++) {    // loop over each atom
                int j = (l*K + k)*Na + n; // index of atom n
                int m = (l*K + k)*N + (Nnb[n]-Nnb[0]); //  starting index of neighbors
                ar[j] = 0.0; 
                ai[j] = 0.0;
                int numnb = Nnb[n+1]-Nnb[n];  // number of neighbors for atom n
                for (int i=0; i<numnb; i++) { // loop over neighbors 
                    ar[j] += Sr[m + i]; // sum over neighbors of atom n
                    ai[j] += Si[m + i];                                                    
                }                    
            }                    
}
template void cpuRadialSphericalHarmonicsSum(double*, double*, double*, double*, int*, int, int, int);
template void cpuRadialSphericalHarmonicsSum(float*, float*, float*, float*, int *, int, int, int);

// core function
template <typename T> void cpuRadialSphericalHarmonicsPower(T *p, T *ar, T *ai, int *indk, int Na, int L, int K)
{   
    // L                      : the maximum degree of spherical harmonics
    // K                      : the maximum degree of spherical Bessel functions
    // Na                     : number of atoms in the simulation domain 
    // indk [K*(K+1)/2]       : store indices for symmetric tensor products of spherical Bessel functions
    // ar [Na*K*(L+1)*(L+1)]  : sum of real spherical harmonics Bessel functions
    // ai [Na*K*(L+1)*(L+1)]  : sum of imag spherical harmonics Bessel functions    
    // p [Na*(L+1)*K*(K+1)/2] : power spectrum components
    
    // Compute the power spectrum components for radial spherical harmonics
    
    int K2 = K*(K+1)/2;    
    //#pragma omp parallel for    
    for (int k=0; k<K2; k++) {
        int k2 = indk[k];
        int k1 = indk[K2+k];
        for (int l=0; l<(L+1); l++) { // loop over orbital quantum number
            for (int n=0; n<Na; n++) {  // loop over each atom          
                int l1 = (k*(L+1) + l)*Na + n; // global index of atom n
                p[l1] = (T) 0.0;
                for (int m=0; m<(2*l+1); m++) { // loop over magnetic quantum number
                    int i1 = ((l*l + m)*K + k1)*Na + n;
                    int j1 = ((l*l + m)*K + k2)*Na + n;
                    p[l1] += ar[i1]*ar[j1] + ai[i1]*ai[j1]; // power spectrum           
                }
            }    
        }
    }                 
}
template void cpuRadialSphericalHarmonicsPower(double*, double*, double*, int*, int, int, int);
template void cpuRadialSphericalHarmonicsPower(float*, float*, float*, int*, int, int, int);

// core function
template <typename T> void cpuRadialSphericalHarmonicsPowerDeriv(T *px, T *py, T *pz, T *ar, T *ai, 
        T *arx, T *aix, T *ary, T *aiy, T *arz, T *aiz, int *indk, int *Nnb, int Na, int L, int K)
{   
    // L                      : the maximum degree of spherical harmonics
    // K                      : the maximum degree of spherical Bessel functions    
    // N                      : total number of neighbors 
    // Na                     : number of global atoms in the simulation domain 
    // Nnb [Na+1]             : a list containing the number of neighbors for each global atom    
    // indk [K*(K+1)/2]       : store indices for symmetric tensor products of spherical Bessel functions
    // ar   [Na*K*(L+1)*(L+1)]: sum of real spherical harmonics Bessel functions
    // ai   [Na*K*(L+1)*(L+1)]: sum of imag spherical harmonics Bessel functions    
    // arx  [N*K*(L+1)*(L+1)] : x-derivative of real spherical harmonics Bessel functions
    // aix  [N*K*(L+1)*(L+1)] : x-derivative of imag spherical harmonics Bessel functions
    // ary  [N*K*(L+1)*(L+1)] : y-derivative of real spherical harmonics Bessel functions
    // aiy  [N*K*(L+1)*(L+1)] : y-derivative of imag spherical harmonics Bessel functions
    // arz  [N*K*(L+1)*(L+1)] : z-derivative of real spherical harmonics Bessel functions
    // aiz  [N*K*(L+1)*(L+1)] : z-derivative of imag spherical harmonics Bessel functions
    // p [(L+1)*K*(K+1)/2]    : power spectrum components
    // px [N*(L+1)*K*(K+1)/2] : x-derivative of power spectrum components
    // py [N*(L+1)*K*(K+1)/2] : y-derivative of power spectrum components
    // pz [N*(L+1)*K*(K+1)/2] : z-derivative of power spectrum components
    
    // Compute partial derivatives of the power spectrum components
    
    int N = Nnb[Na]-Nnb[0]; // total number of neighbors 
    int K2 = K*(K+1)/2;    
    //#pragma omp parallel for         
    for (int k=0; k<K2; k++) {
        int k2 = indk[k];
        int k1 = indk[K2+k];
        for (int l=0; l<(L+1); l++) {     
            for (int n=0; n<Na; n++)  {  // loop over each atom
                int numnb = Nnb[n+1]-Nnb[n];  // number of neighbors for atom n
                for (int i=0; i<numnb; i++) {// loop over each neighbor of atom n                         
                    int j = (k*(L+1) + l)*N + Nnb[n] + i; // index of px, py, pz
                    px[j] = (T) 0.0;
                    py[j] = (T) 0.0;
                    pz[j] = (T) 0.0;
                    for (int m=0; m<(2*l+1); m++) {
                        int i1 = ((l*l + m)*K + k1)*Na + n; // index of ar and ai 
                        int j1 = ((l*l + m)*K + k2)*Na + n; // index of ar and ai                                         
                        int i2 = ((l*l + m)*K + k1)*N + Nnb[n] + i;  // index of arx and aix     
                        int j2 = ((l*l + m)*K + k2)*N + Nnb[n] + i;  // index of arx and aix                                                    
                        px[j] += ar[i1]*arx[j2] + arx[i2]*ar[j1] + ai[i1]*aix[j2] + aix[i2]*ai[j1];
                        py[j] += ar[i1]*ary[j2] + ary[i2]*ar[j1] + ai[i1]*aiy[j2] + aiy[i2]*ai[j1];
                        pz[j] += ar[i1]*arz[j2] + arz[i2]*ar[j1] + ai[i1]*aiz[j2] + aiz[i2]*ai[j1];                    
                    }
                }
            }    
        }
    }
}
template void cpuRadialSphericalHarmonicsPowerDeriv(double*, double*, double*, double*, double*, double*, 
        double*, double*, double*, double*, double*, int*, int*, int, int, int);
template void cpuRadialSphericalHarmonicsPowerDeriv(float*, float*, float*, float*, float*, float*,  
        float*, float*, float*, float*, float*, int*, int*, int, int, int);

// core function
template <typename T> void cpuRadialSphericalHarmonicsBispectrum(T *b, T *ar, T *ai, T *cg, int *indk, 
        int *indl, int *indm, int *rowm, int Nub, int Ncg, int Na, int L, int K)
{   
    // L                      : the maximum degree of spherical harmonics
    // K                      : the maximum degree of spherical Bessel functions
    // Na                     : number of atoms in the simulation domain    
    // Ncg                    : the total number of non-zero Clebsch-Gordan coefficients 
    // Nub                    : the number of non-zero unique spherical harmonics bispectrum components
    // cg   [Ncg*1]           : store non-zero Clebsch-Gordan coefficients
    // indl [Nub*3]           : store the indices of tuple (l2,l1,l)        
    // indk [K*(K+1)/2]       : store indices for symmetric tensor products of spherical Bessel functions
    // indm [Ncg*3]           : store the indices of tuple (m2, m1, m) for non-zero Clebsch-Gordan coefficients
    // rowm [Nub+1]           : store the number of indices of tuple (m2, m1, m) for each tuple (l2, l1, l)    
    // ar  [Na*K*(L+1)*(L+1)] : sum of real spherical harmonics Bessel functions
    // ai  [Na*K*(L+1)*(L+1)] : sum of imag spherical harmonics Bessel functions
    // b [Na*Nub*K*(K+1)/2]   : bispectrum components    
    
    // Compute the bispectrum components for radial spherical harmonics
    
    int K2 = K*(K+1)/2;    
    //#pragma omp parallel for      

    for (int k=0; k<K2; k++) {
        int k2 = indk[k];
        int k1 = indk[K2+k];
        for (int i=0; i<Nub; i++) {        
            int l2 = indl[i];
            int l1 = indl[Nub+i];
            int l = indl[2*Nub+i];                 
            int nm = rowm[i+1]-rowm[i];
            for (int n=0; n<Na; n++) {  // loop over each atom          
                int indn = (k*Nub + i)*Na + n; // global index of atom n
                b[indn] = (T) 0.0;
                for (int j = 0; j<nm; j++) {
                    int m2 = indm[rowm[i]+j];
                    int m1 = indm[Ncg+rowm[i]+j];
                    int m = indm[2*Ncg + rowm[i]+j];                  
                    T a1, b1, a2, b2, a3, b3;                                
                    a1 = ar[((l*l + l + m)*K + k1)*Na + n];
                    b1 = ai[((l*l + l + m)*K + k1)*Na + n];
                    a2 = ar[((l1*l1 + l1 + m1)*K + k2)*Na + n];
                    b2 = ai[((l1*l1 + l1 + m1)*K + k2)*Na + n];
                    a3 = ar[((l2*l2 + l2 + m2)*K + k2)*Na + n];
                    b3 = ai[((l2*l2 + l2 + m2)*K + k2)*Na + n];
                    b[indn] += cg[rowm[i]+j]*(a1*a2*a3 + a2*b1*b3 + a3*b1*b2 - a1*b2*b3);                                              
                }                                
            }
        }
    }                
}
template void cpuRadialSphericalHarmonicsBispectrum(double*, double*, double*, double*, int*, int*, int*, int*, int, int, int, int, int);
template void cpuRadialSphericalHarmonicsBispectrum(float*, float*, float*, float*, int*, int*, int*, int*, int, int, int, int, int);

// core function
template <typename T> void cpuRadialSphericalHarmonicsBispectrumDeriv(T *bx, T *by, T *bz, T *ar, T *ai, 
        T *arx, T *aix, T *ary, T *aiy, T *arz, T *aiz, T*cg, int *indk, int *indl,
        int *indm, int *rowm, int *Nnb, int Na, int Nub, int Ncg, int K)
{   
    // K                      : the maximum degree of spherical Bessel functions
    // N                      : total number of neighbors 
    // Na                     : number of global atoms in the simulation domain     
    // Ncg                    : the total number of non-zero Clebsch-Gordan coefficients 
    // Nub                    : the number of non-zero unique spherical harmonics bispectrum components
    // Nnb [Na+1]             : a list containing the number of neighbors for each global atom        
    // indl [Nub*3]           : store the indices of tuple (l2,l1,l)        
    // indk [K*(K+1)/2]       : store indices for symmetric tensor products of spherical Bessel functions
    // indm [Ncg*3]           : store the indices of tuple (m2, m1, m) for non-zero Clebsch-Gordan coefficients
    // rowm [Nub+1]           : store the number of indices of tuple (m2, m1, m) for each tuple (l2, l1, l)    
    // cg   [Ncg*1]           : store non-zero Clebsch-Gordan coefficients
    // ar   [Na*K*(L+1)*(L+1)]: real part of spherical harmonics Bessel functions
    // ai   [Na*K*(L+1)*(L+1)]: imag part of spherical harmonics Bessel functions    
    // arx  [N*K*(L+1)*(L+1)] : x-derivative of real spherical harmonics Bessel functions
    // aix  [N*K*(L+1)*(L+1)] : x-derivative of imag spherical harmonics Bessel functions
    // ary  [N*K*(L+1)*(L+1)] : y-derivative of real spherical harmonics Bessel functions
    // aiy  [N*K*(L+1)*(L+1)] : y-derivative of imag spherical harmonics Bessel functions
    // arz  [N*K*(L+1)*(L+1)] : z-derivative of real spherical harmonics Bessel functions
    // aiz  [N*K*(L+1)*(L+1)] : z-derivative of imag spherical harmonics Bessel functions
    // b    [Na*Nub*K*(K+1)/2]: bispectrum components        
    // bx   [N*Nub*K*(K+1)/2] : x-derivative of power spectrum components
    // by   [N*Nub*K*(K+1)/2] : y-derivative of power spectrum components
    // bz   [N*Nub*K*(K+1)/2] : z-derivative of power spectrum components
        
    // Compute partial derivatives of the bispectrum components
    
    int N = Nnb[Na]-Nnb[0]; // total number of neighbors 
    int K2 = K*(K+1)/2;    
    for (int k=0; k<K2; k++) {
        int k2 = indk[k];
        int k1 = indk[K2+k];
        for (int i=0; i<Nub; i++) {        
            int l2 = indl[i];
            int l1 = indl[Nub+i];
            int l = indl[2*Nub+i];     
            int nm = rowm[i+1]-rowm[i];
            for (int n=0; n<Na; n++)  {  // loop over each atom
                int numnb = Nnb[n+1]-Nnb[n];  // number of neighbors for atom n
                for (int q=0; q<numnb; q++) {// loop over each neighbor of atom n                         
                    int ii = (k*Nub + i)*N + Nnb[n] + q; // index of bx, by, bz                                        
                    bx[ii] = (T) 0.0;
                    by[ii] = (T) 0.0;
                    bz[ii] = (T) 0.0;                                    
                    for (int j = 0; j<nm; j++) {
                        int m2 = indm[rowm[i]+j];
                        int m1 = indm[Ncg+rowm[i]+j];
                        int m = indm[2*Ncg + rowm[i]+j];                            
                        int n1 = (l*l + l + m)*K + k1;    
                        int n2 = (l1*l1 + l1 + m1)*K + k2;
                        int n3 = (l2*l2 + l2 + m2)*K + k2;                        
                        int lm1 = n1*Na + n; // index of ar and ai 
                        int lm2 = n2*Na + n; // index of ar and ai                                         
                        int lm3 = n3*Na + n; // index of ar and ai                                                                 
                        int mlk1 = n1*N + Nnb[n] + q; // index of arx and aix      
                        int mlk2 = n2*N + Nnb[n] + q; // index of arx and aix            
                        int mlk3 = n3*N + Nnb[n] + q; // index of arx and aix                                                  
                        
                        T c = cg[rowm[i]+j];                    
                        T a1, b1, a2, b2, a3, b3;                                
                        T a1x, b1x, a2x, b2x, a3x, b3x;
                        T a1y, b1y, a2y, b2y, a3y, b3y;
                        T a1z, b1z, a2z, b2z, a3z, b3z;
                        a1 = ar[lm1];
                        b1 = ai[lm1];
                        a2 = ar[lm2];
                        b2 = ai[lm2];
                        a3 = ar[lm3];
                        b3 = ai[lm3];
                        a1x = arx[mlk1];
                        a1y = ary[mlk1];
                        a1z = arz[mlk1];
                        b1x = aix[mlk1];
                        b1y = aiy[mlk1];
                        b1z = aiz[mlk1];
                        a2x = arx[mlk2];
                        a2y = ary[mlk2];
                        a2z = arz[mlk2];
                        b2x = aix[mlk2];
                        b2y = aiy[mlk2];
                        b2z = aiz[mlk2];
                        a3x = arx[mlk3];
                        a3y = ary[mlk3];
                        a3z = arz[mlk3];
                        b3x = aix[mlk3];
                        b3y = aiy[mlk3];
                        b3z = aiz[mlk3];

                        T t1 = a1x*a2*a3 + a1*a2x*a3 + a1*a2*a3x;                
                        T t2 = a2x*b1*b3 + a2*b1x*b3 + a2*b1*b3x;
                        T t3 = a3x*b1*b2 + a3*b1x*b2 + a3*b1*b2x;
                        T t4 = a1x*b2*b3 + a1*b2x*b3 + a1*b2*b3x;
                        bx[ii] += c*(t1 + t2 + t3 - t4);

                        t1 = a1y*a2*a3 + a1*a2y*a3 + a1*a2*a3y;                
                        t2 = a2y*b1*b3 + a2*b1y*b3 + a2*b1*b3y;
                        t3 = a3y*b1*b2 + a3*b1y*b2 + a3*b1*b2y;
                        t4 = a1y*b2*b3 + a1*b2y*b3 + a1*b2*b3y;
                        by[ii] += c*(t1 + t2 + t3 - t4);

                        t1 = a1z*a2*a3 + a1*a2z*a3 + a1*a2*a3z;                
                        t2 = a2z*b1*b3 + a2*b1z*b3 + a2*b1*b3z;
                        t3 = a3z*b1*b2 + a3*b1z*b2 + a3*b1*b2z;
                        t4 = a1z*b2*b3 + a1*b2z*b3 + a1*b2*b3z;
                        bz[ii] += c*(t1 + t2 + t3 - t4);                    
                    }
                }               
            }
        }
    }
}
template void cpuRadialSphericalHarmonicsBispectrumDeriv(double*, double*, double*, double*, double*, double*, 
        double*, double*, double*, double*, double*, double*, int*, int*, int*, int*, int*, int, int, int, int);
template void cpuRadialSphericalHarmonicsBispectrumDeriv(float*, float*, float*, float*, float*, float*,  
        float*, float*, float*, float*, float*, float*, int*, int*, int*, int*, int*, int, int, int, int);

// core function
template <typename T> void cpuRadialSphericalHarmonicsBasis(T *d, T *c, int *atomtype, 
        int Ntype, int Na, int Nbf)
{       
    // Ntype                : number of atom types 
    // Na                   : number of atoms in the simulation domain 
    // Nbf                  : number of basis functions per atom type    
    // atomtype [Ntype]     : a list containing the types of atoms
    // c        [Na*Nbf]    : spectrum components for all atoms
    // d        [Nbf*Ntype] : basis functions based on atom types
        
    for (int t=0; t< Ntype; t++) // for each atom type
        for (int m=0; m<Nbf; m++) // for each basis function
        {
            int k = Nbf*t + m;       // index of the basis function            
            d[k] = (T) 0.0;
            for (int n=0; n<Na; n++) // for each atom n        
            {
                int tn = atomtype[n]; // type of atom n      
                if (tn==t) { // type of atom n is equal to t
                    d[k] += c[m*Nbf + n];
                }
            }                
        }    
}
template void cpuRadialSphericalHarmonicsBasis(double*, double*, int*, int, int, int);
template void cpuRadialSphericalHarmonicsBasis(float*, float*, int*, int, int, int);

// core function
template <typename T> void cpuRadialSphericalHarmonicsBasisDeriv(T *dx, T *dy, T *dz, T *cx, T *cy, T *cz,
        int *atomtype, int *neighlist, int *Nnb, int Ntype, int Na, int Nbf)
{       
    // Ntype                : number of atom types 
    // Na                   : number of atoms in the simulation domain 
    // Nbf                  : number of basis functions per atom type
    // N = Nnb[Na]-Nnb[0] : total number of neighbors     
    // atomtype [Na]        : a list containing the types of atoms
    // neighborlist [N]     : a list containing the indices of neighbors
    // Nnb [Na+1]           : a list containing the number of neighbors for each global atom    
    // c [Na*Nbf]           : spectrum components for all atoms
    // d [Nbf*Ntype]        : basis functions based on atom types
    // cx, cy, cz [N*Nbf]   : derivatives of spectrum components  
    // dx, dy, dz [Na*Nbf*Ntype]: derivatives of the basis functions
    
    int N = Nnb[Na]-Nnb[0]; // total number of neighbors     
    for (int t=0; t< Ntype; t++) // for each atom type
        for (int m=0; m<Nbf; m++) // for each basis function
            for (int n=0; n<Na; n++) // for each atom n            
            {
                int nmt = Na*Nbf*t + Na*m + n; //index of the derivatives of the basis function                
                dx[nmt] = (T) 0.0;
                dy[nmt] = (T) 0.0;
                dz[nmt] = (T) 0.0;                
                int tn = atomtype[n]; // type of atom n                
                int numnb = Nnb[n+1]-Nnb[n];  // number of neighbors for atom n                                                   
                
                for (int q=0; q<numnb; q++) { // loop over each neighbor of atom n                    
                    if (tn==t) {// type of atom n is equal to t
                        // atom n is self, atom i is neighbor
                        int qnm = m*N + (Nnb[n] + q);  // index of the spectrum components of the atom pair (i, n)                    
                        dx[nmt] +=  -cx[qnm];
                        dy[nmt] +=  -cy[qnm];
                        dz[nmt] +=  -cz[qnm]; 
                    }                    
                    
                    int i = neighlist[Nnb[n] + q]; // atom i 
                    int ti = atomtype[i];             // type of atom i                                        
                    if (ti==t) { // type of atom i is equal to t                    
                        int Nnbi = Nnb[i+1]-Nnb[i];  // number of neighbors for atom i
                        int r;
                        for (r=0; r<Nnbi; r++) // loop over each neighbor of atom i
                            if (neighlist[Nnb[i] + r] == n) // if a neighbor of atom i matchs atom n
                                break;

                        int rim = m*N + (Nnb[i] + r); // index of the spectrum components of the atom pair (n, i)                                        
                        // atom n is neighbor, atom i is self
                        dx[nmt] += cx[rim];
                        dy[nmt] += cy[rim];
                        dz[nmt] += cz[rim];
                    }
                }
            }                
}
template void cpuRadialSphericalHarmonicsBasisDeriv(double*, double*, double*, double*, double*, double*, int*, int*, int*, int, int, int);
template void cpuRadialSphericalHarmonicsBasisDeriv(float*, float*, float*, float*, float*, float*, int*, int*, int*, int, int, int);



// template <typename T> void cpuRadialSphericalHarmonicsPowerForce(T *p, T *ar, T *ai, int *indk, int Na, 
//         int Ng, int L, int K)
// {   
//     // L                      : the maximum degree of spherical harmonics
//     // K                      : the maximum degree of spherical Bessel functions
//     // Ng                     : number of global atoms in the simulation domain 
//     // Na                     : number of atoms in the spherical harmonics Bessel functions
//     // indk [K*(K+1)/2]       : store indices for symmetric tensor products of spherical Bessel functions
//     // ar [Na*K*(L+1)*(L+1)]  : sum of real spherical harmonics Bessel functions
//     // ai [Na*K*(L+1)*(L+1)]  : sum of imag spherical harmonics Bessel functions    
//     // p [Ng*(L+1)*K*(K+1)/2] : power spectrum components
//     
//     // px [N*(L+1)*K*(K+1)/2] : x-derivative of power spectrum components
//     // py [N*(L+1)*K*(K+1)/2] : y-derivative of power spectrum components
//     // pz [N*(L+1)*K*(K+1)/2] : z-derivative of power spectrum components
// 
//     // Fpx [Ng*(L+1)*K*(K+1)/2] : x-force of power spectrum components
//     // Fpy [Ng*(L+1)*K*(K+1)/2] : y-force of power spectrum components
//     // Fpz [Ng*(L+1)*K*(K+1)/2] : z-force of power spectrum components    
//     
//     // Compute the power spectrum components for radial spherical harmonics
//     
//     // Sum the spherical harmonics Bessel functions over neighbors
//     int L2 = (L+1)*(L+1);
//     int N = Nnb[Na]-Nnb[0]; // total number of neighbors 
//     for (int l=0; l<L2; l++)
//         for (int k=0; k<K; k++) 
//             for (int n=0; n<Na; n++) {    // loop over each atom
//                 int j = (l*K + k)*Na + n; // index of atom n
//                 int m = (l*K + k)*N + (Nnb[n]-Nnb[0]); //  starting index of neighbors
//                 ar[j] = 0.0; 
//                 ai[j] = 0.0;
//                 int numnb = Nnb[n+1]-Nnb[n];  // number of neighbors for atom n
//                 for (int i=0; i<numnb; i++) { // loop over neighbors 
//                     ar[j] += Sr[m + i]; // sum over neighbors of atom n
//                     ai[j] += Si[m + i];                                                    
//                 }                    
//             }                    
//         
//     //#pragma omp parallel for    
//     for (int k=0; k<Nbf; k++) { // loop over each basis function
//         for (int n=0; n<Na; n++) {  // loop over each atom      
//             int numnb = Nnb[n+1]-Nnb[n];  // number of neighbors for atom n     
//             int m = (l*K2 + k)*N + (Nnb[n]-Nnb[0]); //  starting index of neighbors
//             int l1 = (k*(L+1) + l)*Ng + n; // global index of atom n
//             Fpx[l1] = (T) 0.0;
//             Fpy[l1] = (T) 0.0;
//             Fpz[l1] = (T) 0.0;
//             for (int q=0; q<numnb; q++) { // loop over neighbors                     
//                 Fpx[l1] += px[m + q];
//                 Fpy[l1] += py[m + q];
//                 Fpz[l1] += pz[m + q];
//             }
//         }            
//     }                 
// }
// template void cpuRadialSphericalHarmonicsPower(double*, double*, double*, int*, int, int, int, int);
// template void cpuRadialSphericalHarmonicsPower(float*, float*, float*, int*, int, int, int, int);


// 
// template <typename T> void cpuSphericalHarmonicsBesselSum(T *Sr, T *Si, T *x, T *y, T *z, T *x0, 
//         T *P, T *tmp, T *f, T *fac, T pi, int *Nnb, int Na, int L, int K, int N1, int N2)
// {                        
//     // L                     : the maximum degree of spherical harmonics
//     // K                     : the maximum degree of spherical Bessel functions
//     // N1                    : start index for x, y, z
//     // N2                    : end index for x, y, z
//     // Na                    : number of global atoms
//     // Nnb [Na+1]            : a list containing the number of neighbors for each global atom
//     // x0  [(L+1)*K]         : zeros of pherical Bessel functions
//     // P   [(L+1)*(L+2)/2]   : temporary storage for the recurrence formula
//     // tmp [(L+1)*(L+2)/2]   : temporary storage for the recurrence formula
//     // f   [L+1]             : temporary storage for the recurrence formula
//     // fac                   : factorial look-up table
//     // Sr [K*(L+1)*(L+1)*Na] : real part of spherical harmonics Bessel functions
//     // Si [K*(L+1)*(L+1)*Na] : imag part of spherical harmonics Bessel functions
//     
//     // Compute spherical harmonics Bessel functions: S_{klm}(x,y,z) = g_{lk}(x,y,z)*Y_{lm}(x,y,z)
//     
//     // https://en.wikipedia.org/wiki/Spherical_harmonics
//     // Spherical harmonics functions, Y_{lm}(x,y,z), for l = 0,1,...,L and m = -l,...,0,...,l        
//     //  l = 0:                0
//     //  l = 1:              1 2 3
//     //  l = 2:            4 5 6 7 8
//     //  l = 3:        9 10 11 12 13 14 15
//     //  l = 4:    16 17 18 19 20 21 22 23 24
//     //  ....
//     //  l = L:  L*L L*L+1 ....   ....  (L+1)*(L+1)-1 
//     // The total number of spherical harmonics is (L+1)*(L+1).
//     
//     // https://en.wikipedia.org/wiki/Bessel_function
//     // Spherical Bessel functions, g_{lk}(x,y,z), for l = 0,1,...,L and k = 1,2,...,K
//     //  l = 0:    0    1     2    ....   K-1
//     //  l = 1:    K   K+1   K+2   ....  2K-1
//     //  l = 2:   2K  2K+1  2K+2   ....  3K-1        
//     //  l = 3:   3K  3K+1  3K+2   ....  4K-1        
//     //  l = 4:   4K  4K+1  4K+2   ....  5K-1        
//     //  ....
//     //  l = L:   LK  LK+1  LK+2   .... (L+1)K-1        
//     // The total number of spherical Bessel functions is K*(L+1).
//         
//     // Hence, the total number of spherical harmonics Bessel functions is K*(L+1)*(L+1).    
//     
//     int L1 = L+1;
//     int Ns = K*L1*L1;    // number of spherical harmonics Bessel functions per atom
// //     int Nt = Na*K*L1*L1; // total number of spherical harmonics Bessel functions for Na atoms
// //     for (int i=0; i<Nt; i++)  {
// //         Sr[i] = 0.0;
// //         Si[i] = 0.0;
// //     }
//         
//     //#pragma omp parallel for    
//     for (int i=N1; i<N2; i++) { // loop over neighbors
//                 
//         int ia; // index of global atom connected to neighbor i 
//         for (ia=0; ia<Na; ia++)
//             if ((Nnb[ia]<=i) && (i<Nnb[ia+1]))
//                 break;                 
//         int Nsi = ia*Ns;
//         
//         int k, j, jm;
//         tmp[0] = 1.0;
//         
//         // Cartersian -> Spherical
//         T r = sqrt(x[i]*x[i] + y[i]*y[i] + z[i]*z[i]);
//         T the = acos(z[i]/r);
//         T phi = atan2(y[i],x[i]);
//         
//         // Spherical harmonics for l = 0        
//         int l = 0;
//         T Ylmr = sqrt(1/(4*pi));
//         T Ylmi = 0.0;                
//         T xr, g;
//         // Spherical Bessel for l = 0, see https://en.wikipedia.org/wiki/Bessel_function#Spherical_Bessel_functions:_jn,_yn 
//         for (k=0; k<K; k++) {
//             xr = x0[k*(L+1)+l]*r;                        
//             g = cos(xr)/xr;            
//             j = Nsi + k;                
//             Sr[j] += g*Ylmr;  // real part          
//             Si[j] += g*Ylmi;  // imag part           
//         }
//                 
//         // Spherical harmonics for l = 1;
//         l = 1;
//         T costhe = cos(the);
//         T a = -sin(the);        
//         int m = 0;    
//         P[0] = costhe;
//         P[1] = a;
//         
//         T C = sqrt((2*l+1)*fac[l-m]/(4*pi*fac[l+m]));
//         Ylmr = C*P[0];
//         Ylmi = 0.0;
//         // Spherical Bessel for l = 1, see https://en.wikipedia.org/wiki/Bessel_function#Spherical_Bessel_functions:_jn,_yn
//         for (k=0; k<K; k++) {
//             xr = x0[k*(L+1)+l]*r;                        
//             g = cos(xr)/(xr*xr) + sin(xr)/xr; // Spherical Bessel for l = 1
//             j = Nsi + (l*l + l + m)*K + k;                
//             Sr[j] += g*Ylmr; // real part           
//             Si[j] += g*Ylmi; // imag part                                    
//         }        
//         
//         m = 1;
//         C = sqrt((2*l+1)*fac[l-m]/(4*pi*fac[l+m]));
//         Ylmr = C*cos(phi)*P[1];
//         Ylmi = C*sin(phi)*P[1];                
//         // Spherical Bessel for l = 1, see https://en.wikipedia.org/wiki/Bessel_function#Spherical_Bessel_functions:_jn,_yn
//         for (k=0; k<K; k++) {
//             xr = x0[k*(L+1)+l]*r;                        
//             g = cos(xr)/(xr*xr) + sin(xr)/xr; // Spherical Bessel for l = 1
//             j = Nsi + (l*l + l + m)*K + k; // spherical harmonics Bessel functions for m > 0                
//             Sr[j] += g*Ylmr; // real part           
//             Si[j] += g*Ylmi; // imag part                                                            
//             jm = Nsi + (l*l + l - m)*K + k; // spherical harmonics Bessel functions for m < 0                                 
//             Sr[jm] = -Sr[j]; // real part                      
//             Si[jm] =  Si[j]; // imag part                                                                       
//         }        
//                 
//         for (l=2; l<=L; l++) {                                        
//             // Compute associated Legendre polynomial using recurrence formula
//             // https://en.wikipedia.org/wiki/Associated_Legendre_polynomials            
//             T Pll = P[l-1];
//             tmp[(l-1)] = Pll;
//             P[(l-1)] = (2*(l-1)+1)*(costhe*Pll);
//             P[l] = (2*(l-1)+1)*(a*Pll);             
//             for (m=0; m<l-1; m++) {
//                 Pll = P[m];
//                 P[m] = ((2*(l-1)+1)*(costhe*Pll) - (l+m-1)*tmp[m])/(T (l-m));
//                 tmp[m] = Pll;                
//             }
//                             
//             // Compute spherical harmonics Bessel function at level l for m = 0, 1, ..., l
//             for (m=0; m<=l; m++) {
//                 C = sqrt((2*l+1)*fac[l-m]/(4*pi*fac[l+m])); 
//                 Ylmr = C*(cos(m*phi)*P[m]);
//                 Ylmi = C*(sin(m*phi)*P[m]);                
//                 for (k=0; k<K; k++) {
//                     // Compute the spherical Bessel functions using recurrence formula
//                     xr = x0[k*(L+1)+l]*r;                                            
//                     f[0] = -cos(xr)/xr;
//                     f[1] = -cos(xr)/(xr*xr) - sin(xr)/xr;                    
//                     for (int jk=2; jk<(l+1); jk++) 
//                         f[jk] = ((2*(jk+1)-3)/xr)*f[jk-1] - f[jk-2];                                     
//                     g = -f[l];                                                              
//                     j = Nsi + (l*l + l + m)*K + k;   
//                     Sr[j] += g*Ylmr; // real part            
//                     Si[j] += g*Ylmi; // imag part                                                            
//                 }                
//             }        
//             // Compute the spherical harmonics Bessel functions for m < 0 using symmetry properties
//             C = -1.0;
//             for (m=1; m<=l; m++)  {                 
//                 for (k=0; k<K; k++) {
//                     j =  Nsi + (l*l + l + m)*K + k;  // spherical harmonics Bessel functions for m > 0 
//                     jm = Nsi + (l*l + l - m)*K + k;  // spherical harmonics Bessel functions for m < 0                 
//                     Sr[jm] = C*Sr[j]; // real part
//                     Si[jm] =-C*Si[j]; // imag part                                                                              
//                 }                     
//                 C = C*(-1.0);        
//             }
//         }
//     }    
// }
// template void cpuSphericalHarmonicsBesselSum(double*, double*, double*, double*, double*, double*, 
//         double*, double*, double*, double*, double, int *, int, int, int, int, int);
// template void cpuSphericalHarmonicsBesselSum(float*, float*, float*, float*, float*, float*, 
//         float*, float*, float*, float*, float, int*, int, int, int, int, int);

// template <typename T> void cpuSphericalHarmonics(T *Ylmr, T *Ylmi, T *the, T *phi, 
//                 T *P, T *tmp, T *fac, T pi, int L, int N)
// {                        
//     //#pragma omp parallel for
//     for (int i=0; i<N; i++) {                
//         tmp[0] = 1.0;
//         
//         // l = 0        
//         Ylmr[i] = sqrt(1/(4*pi));
//         Ylmi[i] = 0.0;
//         
//         T costhe = cos(the[i]);
//         T a = -sin(the[i]);
//         
//         int l = 1;
//         int k = 2*l+1;
//         int m = 0;    
//         P[0] = costhe;
//         P[1] = a;
//         
//         m = 0;
//         T C = sqrt(k*fac[l-m]/(4*pi*fac[l+m]));
//         Ylmr[N+i] = C*P[0];
//         Ylmi[N+i] = 0.0;
//         m = 1;
//         C = sqrt(k*fac[l-m]/(4*pi*fac[l+m]));
//         Ylmr[2*N+i] = C*cos(phi[i])*P[1];
//         Ylmi[2*N+i] = C*sin(phi[i])*P[1];                
//         
//         for (l=2; l<=L; l++) {                            
//             
//             T Pll = P[l-1];
//             tmp[(l-1)] = Pll;
//             P[(l-1)] = (2*(l-1)+1)*(costhe*Pll);
//             P[l] = (2*(l-1)+1)*(a*Pll); 
//             for (m=0; m<l-1; m++) {
//                 Pll = P[m];
//                 P[m] = ((2*(l-1)+1)*(costhe*Pll) - (l+m-1)*tmp[m])/(T (l-m));
//                 tmp[m] = Pll;
//             }
//                         
//             int k = 2*l+1;                
//             for (m=0; m<=l; m++) {
//                 C = sqrt(k*fac[l-m]/(4*pi*fac[l+m])); 
//                 Ylmr[(l*(l+1)/2+m)*N + i] = C*(cos(m*phi[i])*P[m]);
//                 Ylmi[(l*(l+1)/2+m)*N + i] = C*(sin(m*phi[i])*P[m]);
//             }
//         }
//     }                        
// }
// template void cpuSphericalHarmonics(double*, double*, double*, double*, double*, double*, double*, double, int, int);
// template void cpuSphericalHarmonics(float*, float*, float*, float*, float*, float*, float*, float, int, int);
// 
// template <typename T> void cpuSphericalHarmonicsDeriv(T *Ylmr, T *Ylmi, T *YlmrThe, T *YlmiThe, T *YlmrPhi, 
//         T *YlmiPhi, T *the, T *phi, T *P, T *tmp, T *dP, T *dtmp, T *fac, T *C, T pi, int L, int N)
// {                        
//     //#pragma omp parallel for
//     for (int i=0; i<N; i++) {    
//         tmp[0] = 1.0;
//         
//         // l = 0
//         Ylmr[i] = sqrt(1/(4*pi));
//         Ylmi[i] = 0.0;
//         YlmrThe[i] = 0.0;
//         YlmrPhi[i] = 0.0;
//         YlmiThe[i] = 0.0;        
//         YlmiPhi[i] = 0.0;
//         
//         T costhe = cos(the[i]);
//         T a = -sin(the[i]);
//         T dcosthe = a;
//         T da = -costhe;
// 
//         int l = 1;
//         int k = 2*l+1;
//         int m = 0;    
//         P[0] = costhe;
//         P[1] = a;
//         dP[0] = dcosthe;    
//         dP[1] = da;
// 
//         m = 0;
//         C[0] = sqrt(k*fac[l-m]/(4*pi*fac[l+m]));
//         Ylmr[N+i] = C[0]*P[0];
//         Ylmi[N+i] = 0.0;
//         YlmrThe[N+i] = C[0]*a;    
//         YlmiThe[N+i] = 0.0;    
//         YlmrPhi[N+i] = 0.0;          
//         YlmiPhi[N+i] = 0.0;          
// 
//         m = 1;
//         C[1] = sqrt(k*fac[l-m]/(4*pi*fac[l+m]));
//         Ylmr[2*N+i] = C[1]*cos(phi[i])*P[1];
//         Ylmi[2*N+i] = C[1]*sin(phi[i])*P[1];                
//         YlmrThe[2*N+i] = C[1]*(cos(phi[i])*da);      
//         YlmiThe[2*N+i] = C[1]*(sin(phi[i])*da);      
//         YlmrPhi[2*N+i] = -C[1]*(sin(phi[i])*P[1]); 
//         YlmiPhi[2*N+i] = C[1]*(cos(phi[i])*P[1]);  
//         
//         for (l=2; l<=L; l++) {                            
//             
//             T Pll = P[l-1];
//             tmp[l-1] = Pll;
//             P[l-1] = (2*(l-1)+1)*(costhe*Pll);
//             P[l] = (2*(l-1)+1)*(a*Pll); 
//             T dPll = dP[l-1];
//             dtmp[l-1] = dPll;
//             dP[l-1] = (2*(l-1)+1)*(costhe*dPll + dcosthe*Pll);   
//             dP[l] = (2*(l-1)+1)*(a*dPll + da*Pll);             
//             for (m=0; m<l-1; m++) {    
//                 Pll = P[m];
//                 P[m] = ((2*(l-1)+1)*(costhe*Pll) - (l+m-1)*tmp[(m)])/(T (l-m));
//                 tmp[m] = Pll;
//                 dPll = dP[m];
//                 dP[m] = ((2*(l-1)+1)*(costhe*dPll+dcosthe*Pll) - (l+m-1)*dtmp[m])/(T (l-m));
//                 dtmp[m] = dPll;
//             }
//                         
//             int k = 2*l+1;  
//             for (m=0; m<=l; m++) {
//                 int i1 = (l*(l+1)/2+m)*N + i;
//                 
//                 C[m] = sqrt(k*fac[l-m]/(4*pi*fac[l+m])); 
//                 Ylmr[i1] = C[m]*(cos(m*phi[i])*P[m]);
//                 Ylmi[i1] = C[m]*(sin(m*phi[i])*P[m]);
//                 YlmrThe[i1] = C[m]*(cos(m*phi[i])*dP[m]); 
//                 YlmiThe[i1] = C[m]*(sin(m*phi[i])*dP[m]); 
//                 YlmrPhi[i1] = -(m*C[m])*(sin(m*phi[i])*P[m]); 
//                 YlmiPhi[i1] = (m*C[m])*(cos(m*phi[i])*P[m]); 
//             }
//         }
//     }                         
// }
// template void cpuSphericalHarmonicsDeriv(double*, double*, double*, double*, double*, double*, double*, double*, 
//         double*, double*, double*, double*, double*, double*, double, int, int);
// template void cpuSphericalHarmonicsDeriv(float*, float*, float*, float*, float*, float*, float*, float*, 
//         float*, float*, float*, float*, float*, float*, float, int, int);

// template <typename T> void cpuSphericalBessel(T *g, T *r, T *x0, T *f, int L, int K, int N)
// {                        
//     //#pragma omp parallel for
//     for (int i=0; i<N; i++) {                
//         int l, k;
//         T x;        
//         for (k=0; k<K; k++) {
//             l = 0;
//             x = x0[k*(L+1)+l]*r[i];                        
//             g[l*N*K + k*N + i] = cos(x)/x;
//             
//             l = 1;
//             x = x0[k*(L+1)+l]*r[i];                        
//             g[l*N*K + k*N + i] = cos(x)/(x*x) + sin(x)/x;
//             
//             for (l=2; l<(L+1); l++) {
//                 x = x0[k*(L+1)+l]*r[i];                     
//                 f[0] = -cos(x)/x;
//                 f[1] = -cos(x)/(x*x) - sin(x)/x;
//                 for (int j=2; j<(l+1); j++) 
//                     f[j] = ((2*(j+1)-3)/x)*f[(j-1)] - f[j-2];                 
//                 g[l*N*K + k*N + i] = -f[l];
//             }            
//         }        
//     }
// }
// template void cpuSphericalBessel(double*, double*, double*, double*, int, int, int);
// template void cpuSphericalBessel(float*, float*, float*, float*, int, int, int);
// 
// template <typename T> void cpuSphericalBesselDeriv(T *g, T *dg, T *r, T *x0, T *f, T *df, int L, int K, int N)
// {                        
//     //#pragma omp parallel for
//     for (int i=0; i<N; i++) {                
//         int l, k;
//         T x;        
//         for (k=0; k<K; k++) {
//             l = 0;
//             x = x0[k*(L+1)+l]*r[i];                        
//             g[l*N*K + k*N + i] = cos(x)/x;
//             dg[l*N*K + k*N + i] = x0[l*K+k]*(-cos(x)/(x*x) - sin(x)/x);
//             
//             l = 1;
//             x = x0[k*(L+1)+l]*r[i];                        
//             g[l*N*K + k*N + i] = cos(x)/(x*x) + sin(x)/x;
//             dg[l*N*K + k*N + i] = x0[l*K+k]*(cos(x)/x - (2*cos(x))/(x*x*x) - (2*sin(x))/(x*x));
//             
//             for (l=2; l<(L+1); l++) {
//                 x = x0[k*(L+1)+l]*r[i];                     
//                 f[0] = -cos(x)/x;
//                 f[1] = -cos(x)/(x*x) - sin(x)/x;
//                 df[0] = x0[k*(L+1)+l]*(cos(x)/(x*x) + sin(x)/x);
//                 df[1] = x0[k*(L+1)+l]*((2*cos(x))/(x*x*x) - cos(x)/x  + (2*sin(x))/(x*x));
//         
//                 for (int j=2; j<(l+1); j++) {
//                     f[j] = ((2*(j+1)-3)/x)*f[j-1] - f[j-2];                 
//                     df[j] = ((2*(j+1)-3)/x)*df[j-1] - x0[k*(L+1)+l]*((2*(j+1)-3)/(x*x))*f[j-1] - df[j-2];        
//                 }
//                 g[l*N*K + k*N + i] = -f[l];
//                 dg[l*N*K + k*N + i] = -df[l];
//             }            
//         }        
//     }
// }
// template void cpuSphericalBesselDeriv(double*, double*, double*, double*, double*, double*, int, int, int);
// template void cpuSphericalBesselDeriv(float*, float*, float*, float*, float*, float*, int, int, int);

// template <typename T> void cpuRadialSphericalHarmonics(T *Sr, T *Si, T *g, T *Ylmr, T *Ylmi, int L, int K, int N)
// {                        
//     //#pragma omp parallel for
//     for (int i=0; i<N; i++)                                
//         for (int k=0; k<K; k++) 
//             for (int l=0; l<(L+1); l++) 
//                 for (int m=0; m<(l+1); m++) {                  
//                     Sr[(l*(l+1)/2+m)*N*K + k*N + i] = g[l*N*K + k*N + i]*Ylmr[(l*(l+1)/2+m)*N + i];
//                     Si[(l*(l+1)/2+m)*N*K + k*N + i] = g[l*N*K + k*N + i]*Ylmi[(l*(l+1)/2+m)*N + i];                                
//                 }
// }
// template void cpuRadialSphericalHarmonics(double*, double*, double*, double*, double*, int, int, int);
// template void cpuRadialSphericalHarmonics(float*, float*, float*, float*, float*, int, int, int);

// template <typename T> void cpuWeightedRadialSphericalHarmonics(T *Sr, T *Si, T *g, T *Ylmr, T *Ylmi, T *w, int L, int K, int N)
// {                        
//     //#pragma omp parallel for
//     for (int i=0; i<N; i++)                                
//         for (int k=0; k<K; k++) 
//             for (int l=0; l<(L+1); l++) 
//                 for (int m=0; m<(l+1); m++) {                  
//                     Sr[(l*(l+1)/2+m)*N*K + k*N + i] = w[i]*g[l*N*K + k*N + i]*Ylmr[(l*(l+1)/2+m)*N + i];
//                     Si[(l*(l+1)/2+m)*N*K + k*N + i] = w[i]*g[l*N*K + k*N + i]*Ylmi[(l*(l+1)/2+m)*N + i];                                
//                 }
// }
// template void cpuWeightedRadialSphericalHarmonics(double*, double*, double*, double*, double*, double*, int, int, int);
// template void cpuWeightedRadialSphericalHarmonics(float*, float*, float*, float*, float*, float*, int, int, int);

// template <typename T> void cpuRadialSphericalHarmonicsSum(T *Sr, T *Si, T *g, T *Ylmr, T *Ylmi, int L, int K, int N)
// {                        
//     //#pragma omp parallel for                
//     for (int k=0; k<K; k++) 
//         for (int l=0; l<(L+1); l++) 
//             for (int m=0; m<(l+1); m++) {                  
//                 Sr[(l*(l+1)/2+m)*K + k] = 0.0;
//                 Si[(l*(l+1)/2+m)*K + k] = 0.0;
//                 for (int i=0; i<N; i++) {                   
//                     Sr[(l*(l+1)/2+m)*K + k] += g[l*N*K + k*N + i]*Ylmr[(l*(l+1)/2+m)*N + i];
//                     Si[(l*(l+1)/2+m)*K + k] += g[l*N*K + k*N + i]*Ylmi[(l*(l+1)/2+m)*N + i];                                
//                 }
//             }
// }
// template void cpuRadialSphericalHarmonicsSum(double*, double*, double*, double*, double*, int, int, int);
// template void cpuRadialSphericalHarmonicsSum(float*, float*, float*, float*, float*, int, int, int);
// 
// template <typename T> void cpuWeightedRadialSphericalHarmonicsSum(T *Sr, T *Si, T *g, T *Ylmr, T *Ylmi, T *w, int L, int K, int N)
// {                        
//     //#pragma omp parallel for                
//     for (int k=0; k<K; k++) 
//         for (int l=0; l<(L+1); l++) 
//             for (int m=0; m<(l+1); m++) {                  
//                 Sr[(l*(l+1)/2+m)*K + k] = 0.0;
//                 Si[(l*(l+1)/2+m)*K + k] = 0.0;
//                 for (int i=0; i<N; i++) {                   
//                     Sr[(l*(l+1)/2+m)*K + k] += w[i]*g[l*N*K + k*N + i]*Ylmr[(l*(l+1)/2+m)*N + i];
//                     Si[(l*(l+1)/2+m)*K + k] += w[i]*g[l*N*K + k*N + i]*Ylmi[(l*(l+1)/2+m)*N + i];                                
//                 }
//             }
// }
// template void cpuWeightedRadialSphericalHarmonicsSum(double*, double*, double*, double*, double*, double*, int, int, int);
// template void cpuWeightedRadialSphericalHarmonicsSum(float*, float*, float*, float*, float*, float*, int, int, int);

// template <typename T> void cpuRadialSphericalHarmonicsDeriv(T *Srx, T *Six, T *Sry, T *Siy, T *Srz, T *Siz, 
//         T *x, T *y, T *z, T *g, T *Ylmr, T *Ylmi, T *gR, T *YlmrThe, T *YlmiThe, T *YlmrPhi, T *YlmiPhi, 
//         int L, int K, int N)
// {                        
//     //#pragma omp parallel for
//     for (int i=0; i<N; i++) {                
//         T r = sqrt(x[i]*x[i] + y[i]*y[i] + z[i]*z[i]);
//         
//         T r2 = r*r;
//         T rxy = sqrt(x[i]*x[i] + y[i]*y[i]);
//         T rxy2 = rxy*rxy;
//         T rr2 = rxy*r2;
// 
//         T Rx = x[i]/r;
//         T Ry = y[i]/r;
//         T Rz = z[i]/r;
//         T Thex = x[i]*z[i]/rr2;
//         T They = y[i]*z[i]/rr2;
//         T Thez = -rxy/r2;
//         T Phix = -y[i]/rxy2;
//         T Phiy = x[i]/rxy2;
//         T Phiz = 0.0;        
//         
//         for (int k=0; k<K; k++) {
//             for (int l=0; l<(L+1); l++) {
//                 int i0 = l*N*K + k*N + i;
//                 T gx = gR[i0]*Rx;
//                 T gy = gR[i0]*Ry;
//                 T gz = gR[i0]*Rz;
//                 for (int m=0; m<(l+1); m++) {                  
//                     int i1 = (l*(l+1)/2+m)*N*K + k*N + i;
//                     int i2 = (l*(l+1)/2+m)*N + i;
//                     
//                     T Ylmrx = YlmrThe[i2]*Thex + YlmrPhi[i2]*Phix;
//                     T Ylmry = YlmrThe[i2]*They + YlmrPhi[i2]*Phiy;
//                     T Ylmrz = YlmrThe[i2]*Thez + YlmrPhi[i2]*Phiz;
//                     T Ylmix = YlmiThe[i2]*Thex + YlmiPhi[i2]*Phix;
//                     T Ylmiy = YlmiThe[i2]*They + YlmiPhi[i2]*Phiy;
//                     T Ylmiz = YlmiThe[i2]*Thez + YlmiPhi[i2]*Phiz;                    
//                     
//                     Srx[i1] = gx*Ylmr[i2] + g[i0]*Ylmrx;
//                     Sry[i1] = gy*Ylmr[i2] + g[i0]*Ylmry;
//                     Srz[i1] = gz*Ylmr[i2] + g[i0]*Ylmrz;
//                     Six[i1] = gx*Ylmi[i2] + g[i0]*Ylmix;
//                     Siy[i1] = gy*Ylmi[i2] + g[i0]*Ylmiy;
//                     Siz[i1] = gz*Ylmi[i2] + g[i0]*Ylmiz;
//                 }
//             }
//         }
//     }
// }
// template void cpuRadialSphericalHarmonicsDeriv(double*, double*, double*, double*, double*, double*, double*, 
//         double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, int, int, int);
// template void cpuRadialSphericalHarmonicsDeriv(float*, float*, float*, float*, float*, float*, float*, 
//         float*, float*, float*, float*, float*, float*, float*, float*, float*, float*, int, int, int);
// 
// template <typename T> void cpuWeightedRadialSphericalHarmonicsDeriv(T *Srx, T *Six, T *Sry, T *Siy, T *Srz, T *Siz, 
//         T *x, T *y, T *z, T *g, T *Ylmr, T *Ylmi, T *gR, T *YlmrThe, T *YlmiThe, T *YlmrPhi, T *YlmiPhi, T *w, 
//         int L, int K, int N)
// {                        
//     //#pragma omp parallel for
//     for (int i=0; i<N; i++) {                
//         T r = sqrt(x[i]*x[i] + y[i]*y[i] + z[i]*z[i]);
//         T r2 = r*r;
//         T rxy = sqrt(x[i]*x[i] + y[i]*y[i]);
//         T rxy2 = rxy*rxy;
//         T rr2 = rxy*r2;
// 
//         T Rx = x[i]/r;
//         T Ry = y[i]/r;
//         T Rz = z[i]/r;
//         T Thex = x[i]*z[i]/rr2;
//         T They = y[i]*z[i]/rr2;
//         T Thez = -rxy/r2;
//         T Phix = -y[i]/rxy2;
//         T Phiy = x[i]/rxy2;
//         T Phiz = 0.0;        
//         
//         for (int k=0; k<K; k++) {
//             for (int l=0; l<(L+1); l++) {
//                 int i0 = l*N*K + k*N + i;
//                 T gx = w[i]*gR[i0]*Rx;
//                 T gy = w[i]*gR[i0]*Ry;
//                 T gz = w[i]*gR[i0]*Rz;
//                 T wg = w[i]*g[i0];
//                 for (int m=0; m<(l+1); m++) {                  
//                     int i1 = (l*(l+1)/2+m)*N*K + k*N + i;
//                     int i2 = (l*(l+1)/2+m)*N + i;
//                                         
//                     T Ylmrx = YlmrThe[i2]*Thex + YlmrPhi[i2]*Phix;
//                     T Ylmry = YlmrThe[i2]*They + YlmrPhi[i2]*Phiy;
//                     T Ylmrz = YlmrThe[i2]*Thez + YlmrPhi[i2]*Phiz;
//                     T Ylmix = YlmiThe[i2]*Thex + YlmiPhi[i2]*Phix;
//                     T Ylmiy = YlmiThe[i2]*They + YlmiPhi[i2]*Phiy;
//                     T Ylmiz = YlmiThe[i2]*Thez + YlmiPhi[i2]*Phiz;                    
//                     
//                     Srx[i1] = gx*Ylmr[i2] + wg*Ylmrx;
//                     Sry[i1] = gy*Ylmr[i2] + wg*Ylmry;
//                     Srz[i1] = gz*Ylmr[i2] + wg*Ylmrz;
//                     Six[i1] = gx*Ylmi[i2] + wg*Ylmix;
//                     Siy[i1] = gy*Ylmi[i2] + wg*Ylmiy;
//                     Siz[i1] = gz*Ylmi[i2] + wg*Ylmiz;
//                 }
//             }
//         }
//     }
// }
// template void cpuWeightedRadialSphericalHarmonicsDeriv(double*, double*, double*, double*, double*, 
//         double*, double*, double*, double*, double*, double*, double*, double*, double*, 
//         double*, double*, double*, double*, int, int, int);
// template void cpuWeightedRadialSphericalHarmonicsDeriv(float*, float*, float*, float*, float*, 
//         float*, float*, float*, float*, float*, float*, float*, float*, float*, float*, 
//         float*, float*, float*, int, int, int);

#endif


