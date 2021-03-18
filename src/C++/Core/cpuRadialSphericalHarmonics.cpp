#ifndef __CPURADIALSPHERICALHARMONICS
#define __CPURADIALSPHERICALHARMONICS
                       
template <typename T> void cpuRadialSphericalHarmonics(T *Sr, T *Si, T *Ylmr, T *Ylmi, T *g, int L, int K, int N)
{                        
    // Ylmr  [N*(L+1)*(L+1)] : real spherical harmonics functions
    // Ylmi  [N*(L+1)*(L+1)] : imag spherical harmonics functions    
    // g  [N*K*(L+1)]:         spherical Bessel functions
    // Sr  [N*K*(L+1)*(L+1)] : real radial spherical harmonics functions
    // Si  [N*K*(L+1)*(L+1)] : imag radial spherical harmonics functions
        
    for (int k=0; k<K; k++) 
        for (int l=0; l<(L+1); l++) 
            for (int m=-l; m<(l+1); m++) 
                for (int i=0; i<N; i++) {
                    int Nl = (l*l + l + m)*N + i;
                    int Nk = (l*K + k)*N + i;
                    int Nj = ((l*l + l + m)*K + k)*N + i;                
                    Sr[Nj] = g[Nk]*Ylmr[Nl];
                    Si[Nj] = g[Nk]*Ylmi[Nl];                                                
                }        
}
template void cpuRadialSphericalHarmonics(double*, double*, double*, double*, double*, int, int, int);
template void cpuRadialSphericalHarmonics(float*, float*, float*, float*, float*, int, int, int);

template <typename T> void cpuRadialSphericalHarmonicsWithDeriv(T *Sr, T *Si, T *Srx, T *Six, T *Sry, T *Siy, 
        T *Srz, T *Siz, T *Ylmr, T *Ylmi, T *Ylmrx, T *Ylmix, T *Ylmry, T *Ylmiy, 
        T *Ylmrz, T *Ylmiz,  T *g, T *gx, T *gy, T *gz, int L, int K, int N)
{                        
    // Ylmr  [N*(L+1)*(L+1)] : real spherical harmonics functions
    // Ylmi  [N*(L+1)*(L+1)] : imag spherical harmonics functions    
    // g  [N*K*(L+1)]:         spherical Bessel functions
    // Sr  [N*K*(L+1)*(L+1)] : real radial spherical harmonics functions
    // Si  [N*K*(L+1)*(L+1)] : imag radial spherical harmonics functions
    
    for (int k=0; k<K; k++) 
        for (int l=0; l<(L+1); l++) 
            for (int m=-l; m<(l+1); m++) 
                for (int i=0; i<N; i++) {
                    int Nl = (l*l + l + m)*N + i;
                    int Nk = (l*K + k)*N + i;
                    int Nj = ((l*l + l + m)*K + k)*N + i;                
                    Sr[Nj] = g[Nk]*Ylmr[Nl];
                    Si[Nj] = g[Nk]*Ylmi[Nl];    
                    Srx[Nj] = gx[Nk]*Ylmr[Nl] + g[Nk]*Ylmrx[Nl];
                    Sry[Nj] = gy[Nk]*Ylmr[Nl] + g[Nk]*Ylmry[Nl];
                    Srz[Nj] = gz[Nk]*Ylmr[Nl] + g[Nk]*Ylmrz[Nl];
                    Six[Nj] = gx[Nk]*Ylmi[Nl] + g[Nk]*Ylmix[Nl];
                    Siy[Nj] = gy[Nk]*Ylmi[Nl] + g[Nk]*Ylmiy[Nl];
                    Siz[Nj] = gz[Nk]*Ylmi[Nl] + g[Nk]*Ylmiz[Nl];                                   
                }                
}
template void cpuRadialSphericalHarmonicsWithDeriv(double*, double*, double*, double*, double*, double*, double*, double*, 
        double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, int, int, int);
template void cpuRadialSphericalHarmonicsWithDeriv(float*, float*, float*, float*, float*, float*, float*, float*, 
        float*, float*, float*, float*, float*, float*, float*, float*, float*, float*, float*, float*, int, int, int);

template <typename T> void cpuRadialSphericalHarmonicsNeighborSum(T *ar, T *ai, T *Sr, T *Si, 
        int *Nnb, int Na, int L, int K)
{                            
    // L                         : the maximum degree of spherical harmonics
    // K                         : the maximum degree of spherical Bessel functions
    // Na                        : number of atoms in the simulation domain 
    // Nnb [Na+1]                : a list containing the number of neighbors for each global atom    
    // N = Nnb[Na]-Nnb[0]        : total number of neighbors
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
template void cpuRadialSphericalHarmonicsNeighborSum(double*, double*, double*, double*, int*, int, int, int);
template void cpuRadialSphericalHarmonicsNeighborSum(float*, float*, float*, float*, int *, int, int, int);

template <typename T> void cpuRadialSphericalHarmonicsPowerSpectrum(T *p, T *ar, T *ai, int *indk, int Na, int L, int K)
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
template void cpuRadialSphericalHarmonicsPowerSpectrum(double*, double*, double*, int*, int, int, int);
template void cpuRadialSphericalHarmonicsPowerSpectrum(float*, float*, float*, int*, int, int, int);

template <typename T> void cpuRadialSphericalHarmonicsPowerSpectrumDeriv(T *px, T *py, T *pz, T *ar, T *ai, 
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
template void cpuRadialSphericalHarmonicsPowerSpectrumDeriv(double*, double*, double*, double*, double*, double*, 
        double*, double*, double*, double*, double*, int*, int*, int, int, int);
template void cpuRadialSphericalHarmonicsPowerSpectrumDeriv(float*, float*, float*, float*, float*, float*,  
        float*, float*, float*, float*, float*, int*, int*, int, int, int);

template <typename T> void cpuRadialSphericalHarmonicsPowerSpectrumDeriv2(T *pd, T *ar, T *ai, 
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
    // pd [3*N*(L+1)*K*(K+1)/2] : derivatives of power spectrum components
    
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
                    pd[3*j+0] = (T) 0.0;
                    pd[3*j+1] = (T) 0.0;
                    pd[3*j+2] = (T) 0.0;
                    for (int m=0; m<(2*l+1); m++) {
                        int i1 = ((l*l + m)*K + k1)*Na + n; // index of ar and ai 
                        int j1 = ((l*l + m)*K + k2)*Na + n; // index of ar and ai                                         
                        int i2 = ((l*l + m)*K + k1)*N + Nnb[n] + i;  // index of arx and aix     
                        int j2 = ((l*l + m)*K + k2)*N + Nnb[n] + i;  // index of arx and aix                                                    
                        pd[3*j+0] += ar[i1]*arx[j2] + arx[i2]*ar[j1] + ai[i1]*aix[j2] + aix[i2]*ai[j1];
                        pd[3*j+1] += ar[i1]*ary[j2] + ary[i2]*ar[j1] + ai[i1]*aiy[j2] + aiy[i2]*ai[j1];
                        pd[3*j+2] += ar[i1]*arz[j2] + arz[i2]*ar[j1] + ai[i1]*aiz[j2] + aiz[i2]*ai[j1];                    
                    }
                }
            }    
        }
    }
}
template void cpuRadialSphericalHarmonicsPowerSpectrumDeriv2(double*, double*, double*, double*, 
        double*, double*, double*, double*, double*, int*, int*, int, int, int);
template void cpuRadialSphericalHarmonicsPowerSpectrumDeriv2(float*, float*, float*, float*,  
        float*, float*, float*, float*, float*, int*, int*, int, int, int);

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

template <typename T> void cpuRadialSphericalHarmonicsBispectrumDeriv(T *bx, T *by, T *bz, T *ar, T *ai, 
        T *arx, T *aix, T *ary, T *aiy, T *arz, T *aiz, T*cg, int *indk, int *indl,
        int *indm, int *rowm, int *Nnb, int Nub, int Ncg, int Na, int K)
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

template <typename T> void cpuRadialSphericalHarmonicsBispectrumDeriv2(T *bd, T *ar, T *ai, 
        T *arx, T *aix, T *ary, T *aiy, T *arz, T *aiz, T*cg, int *indk, int *indl,
        int *indm, int *rowm, int *Nnb, int Nub, int Ncg, int Na, int K)
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
    // bd   [3*N*Nub*K*(K+1)/2] : derivatives of power spectrum components
        
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
                    bd[3*ii+0] = (T) 0.0;
                    bd[3*ii+1] = (T) 0.0;
                    bd[3*ii+2] = (T) 0.0;                                    
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
                        bd[3*ii+0] += c*(t1 + t2 + t3 - t4);

                        t1 = a1y*a2*a3 + a1*a2y*a3 + a1*a2*a3y;                
                        t2 = a2y*b1*b3 + a2*b1y*b3 + a2*b1*b3y;
                        t3 = a3y*b1*b2 + a3*b1y*b2 + a3*b1*b2y;
                        t4 = a1y*b2*b3 + a1*b2y*b3 + a1*b2*b3y;
                        bd[3*ii+1] += c*(t1 + t2 + t3 - t4);

                        t1 = a1z*a2*a3 + a1*a2z*a3 + a1*a2*a3z;                
                        t2 = a2z*b1*b3 + a2*b1z*b3 + a2*b1*b3z;
                        t3 = a3z*b1*b2 + a3*b1z*b2 + a3*b1*b2z;
                        t4 = a1z*b2*b3 + a1*b2z*b3 + a1*b2*b3z;
                        bd[3*ii+2] += c*(t1 + t2 + t3 - t4);                    
                    }
                }               
            }
        }
    }
}
template void cpuRadialSphericalHarmonicsBispectrumDeriv2(double*, double*, double*, double*, 
        double*, double*, double*, double*, double*, double*, int*, int*, int*, int*, int*, int, int, int, int);
template void cpuRadialSphericalHarmonicsBispectrumDeriv2(float*, float*, float*, float*,  
        float*, float*, float*, float*, float*, float*, int*, int*, int*, int*, int*, int, int, int, int);

template <typename T> void cpuRadialSphericalHarmonicsSpectrum(T *c, T *ar, T *ai, T *sr, T *si, T *cg, int *indk, 
        int *indl, int *indm, int *rowm, int *Nnb, int Nub, int Ncg, int Na, int L, int K, int spectrum)
{
    //cpuSphericalHarmonicsBessel(sr, si, xij, x0, P, tmp, f, fac, M_PI, L, K, N);
    cpuRadialSphericalHarmonicsNeighborSum(ar, ai, sr, si, Nnb, Na, L, K);

    if (spectrum==0) {  // power spectrum          
        cpuRadialSphericalHarmonicsPowerSpectrum(c, ar, ai, indk, Na, L, K);
    }
    else if (spectrum==1) { // bispectrum            
        cpuRadialSphericalHarmonicsBispectrum(c, ar, ai, cg, indk, indl, indm, rowm, Nub, Ncg, Na, L, K);
    }
    else if (spectrum==2) { // power spectrum and bispectrum         
        int npower = (L+1)*K*(K+1)/2;      // number of power components
        cpuRadialSphericalHarmonicsPowerSpectrum(c, ar, ai, indk, Na, L, K);
        cpuRadialSphericalHarmonicsBispectrum(&c[Na*npower], ar, ai, cg, indk, indl, indm, rowm, Nub, Ncg, Na, L, K);            
    }
}
template void cpuRadialSphericalHarmonicsSpectrum(double*, double*, double*, double*, double*, double*, int*, int*, 
        int*, int*, int*, int, int, int, int, int, int);
template void cpuRadialSphericalHarmonicsSpectrum(float*, float*, float*, float*, float*, float*, int*, int*, 
        int*, int*, int*, int, int, int, int, int, int);
 
template <typename T> void cpuRadialSphericalHarmonicsSpectrumDeriv(T *cd, T *ar, T *ai, T *sr, T *si, 
        T *srx, T *six, T *sry, T *siy, T *srz, T *siz, T *cg, int *indk, int *indl, int *indm, 
        int *rowm, int *Nnb, int Nub, int Ncg, int Na, int L, int K, int spectrum)
{    
//     cpuSphericalHarmonicsBesselWithDeriv(sr, si, srx, six, sry, siy, srz, siz, xij,
//         x0, P, tmp, f, dP, dtmp, df, fac, M_PI, L, K, N)
    cpuRadialSphericalHarmonicsNeighborSum(ar, ai, sr, si, Nnb, Na, L, K);
    
    if (spectrum==0) {  // power spectrum          
        cpuRadialSphericalHarmonicsPowerSpectrumDeriv2(cd, ar, ai, srx, six, sry, siy, srz, siz, 
                indk, Nnb, Na, L, K);        
    }
    else if (spectrum==1) { // bispectrum                    
        cpuRadialSphericalHarmonicsBispectrumDeriv2(cd, ar, ai, srx, six, sry, siy, srz, siz, 
               cg, indk, indl, indm, rowm, Nnb, Nub, Ncg,  Na, K);                
    }
    else if (spectrum==2) { // power spectrum and bispectrum         
        int npower = (L+1)*K*(K+1)/2;      // number of power components
        int N = Nnb[Na]-Nnb[0]; // total number of neighbors 
        cpuRadialSphericalHarmonicsPowerSpectrumDeriv2(cd, ar, ai, srx, six, sry, siy, srz, siz, 
                indk, Nnb, Na, L, K);        
        cpuRadialSphericalHarmonicsBispectrumDeriv2(&cd[3*N*npower], ar, ai, srx, six, sry, siy, srz, siz, 
               cg, indk, indl, indm, rowm, Nnb, Nub, Ncg,  Na, K);                        
    }
}
template void cpuRadialSphericalHarmonicsSpectrumDeriv(double*, double*, double*, double*, double*, double*, double*, 
        double*, double*, double*, double*, double*, int*, int*, int*, int*, int*, int, int, int, int, int, int);
template void cpuRadialSphericalHarmonicsSpectrumDeriv(float*, float*, float*, float*, float*, float*, float*,
        float*, float*, float*, float*, float*, int*, int*, int*, int*, int*, int, int, int, int, int, int);

template <typename T> void cpuRadialSphericalHarmonicsSpectrumDeriv(T *cd, T *ar, T *ai, 
        T *srx, T *six, T *sry, T *siy, T *srz, T *siz, T *cg, int *indk, int *indl, int *indm, 
        int *rowm, int *Nnb, int Nub, int Ncg, int Na, int L, int K, int spectrum)
{    
//     cpuSphericalHarmonicsBesselWithDeriv(sr, si, srx, six, sry, siy, srz, siz, xij,
//         x0, P, tmp, f, dP, dtmp, df, fac, M_PI, L, K, N)    
    if (spectrum==0) {  // power spectrum          
        cpuRadialSphericalHarmonicsPowerSpectrumDeriv2(cd, ar, ai, srx, six, sry, siy, srz, siz, 
                indk, Nnb, Na, L, K);        
    }
    else if (spectrum==1) { // bispectrum                    
        cpuRadialSphericalHarmonicsBispectrumDeriv2(cd, ar, ai, srx, six, sry, siy, srz, siz, 
               cg, indk, indl, indm, rowm, Nnb, Nub, Ncg,  Na, K);                
    }
    else if (spectrum==2) { // power spectrum and bispectrum         
        int npower = (L+1)*K*(K+1)/2;      // number of power components
        int N = Nnb[Na]-Nnb[0]; // total number of neighbors 
        cpuRadialSphericalHarmonicsPowerSpectrumDeriv2(cd, ar, ai, srx, six, sry, siy, srz, siz, 
                indk, Nnb, Na, L, K);        
        cpuRadialSphericalHarmonicsBispectrumDeriv2(&cd[3*N*npower], ar, ai, srx, six, sry, siy, srz, siz, 
               cg, indk, indl, indm, rowm, Nnb, Nub, Ncg,  Na, K);                        
    }
}
template void cpuRadialSphericalHarmonicsSpectrumDeriv(double*, double*, double*, double*, double*, double*, double*, 
        double*, double*, double*, int*, int*, int*, int*, int*, int, int, int, int, int, int);
template void cpuRadialSphericalHarmonicsSpectrumDeriv(float*, float*, float*, float*,  float*, float*, float*,
        float*, float*, float*, int*, int*, int*, int*, int*, int, int, int, int, int, int);

template <typename T> void cpuForceDecomposition(T *f, T *fij, int *ai, int *aj, 
        int inum, int ijnum, int Nbf)
{    
    // fij : 3*ijnum*Nbf 
    // f   : 3*inum*Nbf
    for (int m=0; m<Nbf; m++) {// for each basis function
        int n = 3*inum*m;
        int p = 3*ijnum*m;
        for (int ii=0; ii<ijnum; ii++) {  // for each atom pair ij in the simulation box     
            int i = 3*ai[ii] + n;       // atom i
            int j = 3*aj[ii] + n;       // atom j
            int l = 3*ii + p;        
            // use atomicAdd on gpu
            f[0+i] += -fij[0+l]; 
            f[1+i] += -fij[1+l]; 
            f[2+i] += -fij[2+l]; 
            f[0+j] +=  fij[0+l];         
            f[1+j] +=  fij[1+l];         
            f[2+j] +=  fij[2+l];         
        }
    }
}
template void cpuForceDecomposition(double*, double*, int*, int*, int, int, int);
template void cpuForceDecomposition(float*, float*, int*, int*, int, int, int);

template <typename T> void cpuForceDecomposition(T *f, T *fij, int *ai, int *aj, int ijnum)
{    
    // fij : 3*ijnum 
    // f   : 3*inum
    for (int ii=0; ii<ijnum; ii++) {  // for each atom pair ij in the simulation box     
        int i = 3*ai[ii];       // atom i
        int j = 3*aj[ii];       // atom j
        int l = 3*ii; 
        // use atomicAdd on gpu
        f[0+i] += -fij[0+l]; 
        f[1+i] += -fij[1+l]; 
        f[2+i] += -fij[2+l]; 
        f[0+j] +=  fij[0+l];         
        f[1+j] +=  fij[1+l];         
        f[2+j] +=  fij[2+l];         
    }
}
template void cpuForceDecomposition(double*, double*, int*, int*, int);
template void cpuForceDecomposition(float*, float*, int*, int*, int);

template <typename T> void cpuCenterAtomDecomposition(T *e, T *ei, int *ilist, int Na)
{    
    for (int ii=0; ii<Na; ii++) {  // for each atom i in the simulation box     
        int i = ilist[ii];       // atom i        
        e[i] += ei[ii];             
    }    
//     for (int m=0; m<Nbf; m++) {// for each basis function
//         int n = inum*m;
//         int k = Na*m;
//         for (int ii=0; ii<Na; ii++) {  // for each atom i in the simulation box     
//             int i = ilist[ii];       // atom i        
//             e[i + n] = ei[ii + k];             
//         }
//     }
}
template void cpuCenterAtomDecomposition(double*, double*, int*, int);
template void cpuCenterAtomDecomposition(float*, float*, int*, int);

template <typename T> void cpuCenterAtomDecomposition(T *f, T *fij, int *ilist, int *anumsum, 
        int inum, int ijnum, int Na, int Nbf)
{    
    for (int m=0; m<Nbf; m++) {// for each basis function
        int n = 3*inum*m;
        int p = 3*ijnum*m;
        for (int ii=0; ii<Na; ii++) {  // for each atom i in the simulation box     
            int i = ilist[ii];       // atom i        
            int start = anumsum[ii];   
            int k = anumsum[ii+1]-start;    // number of neighbors around i             
            int idim = 3*i + n;        
            for (int l=0; l<k ; l++) {   // loop over each atom j around atom i
                int ndim = 3*(start + l) + p;                
                f[0+idim] += -fij[0+ndim]; 
                f[1+idim] += -fij[1+ndim]; 
                f[2+idim] += -fij[2+ndim]; 
            }        
        }
    }
}
template void cpuCenterAtomDecomposition(double*, double*, int*, int*, int, int, int, int);
template void cpuCenterAtomDecomposition(float*, float*, int*, int*, int, int, int, int);

template <typename T> void cpuNeighborAtomDecomposition(T *f, T *fij, int *jlist, int *bnumsum, int *index, 
        int inum, int ijnum, int jnum, int Nbf)
{        
    for (int m=0; m<Nbf; m++) {// for each basis function
        int n = 3*inum*m;
        int p = 3*ijnum*m;    
        for (int ii=0; ii<jnum; ii++) {  // for each atom j in the simulation box     
            int j = jlist[ii];       // atom j        
            int start = bnumsum[ii];   
            int k = bnumsum[ii+1]-start;        // number of neighbors around j             
            int jdim = 3*j + n;        
            for (int l=0; l<k ; l++) {   // loop over each atom l around atom j                
                int ndim = 3*index[start + l] + p;     
                f[0+jdim] += fij[0+ndim]; 
                f[1+jdim] += fij[1+ndim]; 
                f[2+jdim] += fij[2+ndim];                 
            }        
        }    
    }
}
template void cpuNeighborAtomDecomposition(double*, double*, int*, int*, int*, int, int, int, int);
template void cpuNeighborAtomDecomposition(float*, float*, int*, int*, int*, int, int, int, int);

template <typename T> void cpuCenterAtomDecomposition(T *f, T *fij, int *ilist, int *anumsum, int inum)
{    
    for (int ii=0; ii<inum; ii++) {  // for each atom i in the simulation box     
        int i = ilist[ii];       // atom i        
        int start = anumsum[ii];   
        int k = anumsum[ii+1]-start;    // number of neighbors around i             
        int idim = 3*i;        
        for (int l=0; l<k ; l++) {   // loop over each atom j around atom i
            int ndim = 3*(start + l);                
            f[0+idim] += -fij[0+ndim]; 
            f[1+idim] += -fij[1+ndim]; 
            f[2+idim] += -fij[2+ndim]; 
        }        
    }
}
template void cpuCenterAtomDecomposition(double*, double*, int*, int*, int);
template void cpuCenterAtomDecomposition(float*, float*, int*, int*, int);


template <typename T> void cpuNeighborAtomDecomposition(T *f, T *fij, int *jlist, int *bnumsum, int *index, int jnum)
{        
    for (int ii=0; ii<jnum; ii++) {  // for each atom j in the simulation box     
        int j = jlist[ii];       // atom j        
        int start = bnumsum[ii];   
        int k = bnumsum[ii+1]-start;        // number of neighbors around j             
        int jdim = 3*j;        
        for (int l=0; l<k ; l++) {   // loop over each atom l around atom j                
            int ndim = 3*index[start + l];     
            f[0+jdim] += fij[0+ndim]; 
            f[1+jdim] += fij[1+ndim]; 
            f[2+jdim] += fij[2+ndim];                 
        }        
    }    
}
template void cpuNeighborAtomDecomposition(double*, double*, int*, int*, int*, int);
template void cpuNeighborAtomDecomposition(float*, float*, int*, int*, int*, int);

#endif

