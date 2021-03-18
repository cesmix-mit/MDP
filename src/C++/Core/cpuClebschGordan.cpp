#ifndef __CPUCLEBSCHGORDAN
#define __CPUCLEBSCHGORDAN 

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
template <typename T> void cpuSphericalHarmonicsAtomSum(T *Ylmr, T *Ylmi, T *x, T *y, T *z, 
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
template void cpuSphericalHarmonicsAtomSum(double*, double*, double*, double*, double*, double*, double*, double*, double, int, int);
template void cpuSphericalHarmonicsAtomSum(float*, float*, float*, float*, float*, float*, float*, float*, float, int, int);

// core function
template <typename T> void cpuSphericalHarmonicsBispectrum(T *b, T *Ylmr, T *Ylmi, T *fac, int L)
{                   
    // L                    : the maximum degree of spherical harmonics
    // Ylmr      [(L+1)*(L+1)]: real part of spherical harmonics
    // Ylmi      [(L+1)*(L+1)]: imag part of spherical harmonics
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
                                a1 = Ylmr[l*l + l + m];
                                b1 = Ylmi[l*l + l + m];
                                a2 = Ylmr[l1*l1 + l1 + m1];
                                b2 = Ylmi[l1*l1 + l1 + m1];
                                a3 = Ylmr[l2*l2 + l2 + m2];
                                b3 = Ylmi[l2*l2 + l2 + m2];
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

#endif


