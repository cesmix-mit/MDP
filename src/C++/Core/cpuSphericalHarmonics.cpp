#ifndef __CPUSPHERICALHARMONICS
#define __CPUSPHERICALHARMONICS

template <typename T> void cpuSphericalHarmonics(T *Sr, T *Si, T *x, T *y, T *z, 
                T *P, T *tmp, T *fac, T pi, int L, int N)
{                        
    // L                     : the maximum degree of spherical harmonics
    // N                     : length of x, y, z
    // P   [(L+1)*(L+2)/2]   : temporary storage for the recurrence formula
    // tmp [(L+1)*(L+2)/2]   : temporary storage for the recurrence formula
    // fac                   : factorial look-up table
    // Sr  [N*(L+1)*(L+1)] : real part of spherical harmonics  functions
    // Si  [N*(L+1)*(L+1)] : imag part of spherical harmonics  functions
        
    // Compute spherical harmonics  functions: S_{klm}(x,y,z) = g_{lk}(x,y,z)*Y_{lm}(x,y,z)
    
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
                
    //#pragma omp parallel for    
    for (int i=0; i<N; i++) { // loop over each neighbor atom               
        int j, jm;
        tmp[0] = 1.0;
        
        // Cartersian -> Spherical
        T r = sqrt(x[i]*x[i] + y[i]*y[i] + z[i]*z[i]);
        T the = acos(z[i]/r);
        T phi = atan2(y[i],x[i]);
        
        // Spherical harmonics for l = 0        
        int l = 0;
        T Ylmr = sqrt(1/(4*pi));
        T Ylmi = 0.0;                
        Sr[i] = Ylmr;  // real part                   
        Si[i] = Ylmi;  // imag part       
                
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
        j = (l*l + l + m)*N + i;                
        Sr[j] = Ylmr; // real part           
        Si[j] = Ylmi; // imag part                                    
                
        m = 1;
        C = sqrt((2*l+1)*fac[l-m]/(4*pi*fac[l+m]));
        Ylmr = C*cos(phi)*P[1];
        Ylmi = C*sin(phi)*P[1];                
        j = (l*l + l + m)*N + i; // spherical harmonics functions for m > 0                
        Sr[j] = Ylmr; // real part           
        Si[j] = Ylmi; // imag part     
        
        jm = (l*l + l - m)*N + i; // spherical harmonics  functions for m < 0                
        Sr[jm] = -Sr[j]; // real part                      
        Si[jm] =  Si[j]; // imag part                                                                       
                
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
                            
            // Compute spherical harmonics  function at level l for m = 0, 1, ..., l
            for (m=0; m<=l; m++) {
                C = sqrt((2*l+1)*fac[l-m]/(4*pi*fac[l+m])); 
                Ylmr = C*(cos(m*phi)*P[m]);
                Ylmi = C*(sin(m*phi)*P[m]);                                
                j = (l*l + l + m)*N + i;   
                Sr[j] = Ylmr; // real part            
                Si[j] = Ylmi; // imag part                                                            
            }        
            
            // Compute the spherical harmonics  functions for m < 0 using symmetry properties
            C = -1.0;
            for (int m=1; m<=l; m++)  {                 
                j =  (l*l + l + m)*N + i;  // spherical harmonics  functions for m > 0 
                jm = (l*l + l - m)*N + i;  // spherical harmonics  functions for m < 0                 
                Sr[jm] = C*Sr[j]; // real part
                Si[jm] =-C*Si[j]; // imag part                                                                             
                C = C*(-1.0);        
            }
        }
    }          
}
template void cpuSphericalHarmonics(double*, double*, double*, double*, double*, double*, 
         double*, double*, double, int, int);
template void cpuSphericalHarmonics(float*, float*, float*, float*, float*, float*, 
         float*, float*, float, int, int);

template <typename T> void cpuSphericalHarmonics(T *Sr, T *Si, T *xij, 
                T *P, T *tmp, T *fac, T pi, int L, int N)
{                                        
    //#pragma omp parallel for    
    for (int i=0; i<N; i++) { // loop over each neighbor atom               
        int j, jm;
        tmp[0] = 1.0;
        
        T x = xij[i*3];
        T y = xij[i*3+1];
        T z = xij[i*3+2];        
        // Cartersian -> Spherical
        T r = sqrt(x*x + y*y + z*z);
        T the = acos(z/r);
        T phi = atan2(y,x);
        
        // Spherical harmonics for l = 0        
        int l = 0;
        T Ylmr = sqrt(1/(4*pi));
        T Ylmi = 0.0;                
        Sr[i] = Ylmr;  // real part                   
        Si[i] = Ylmi;  // imag part       
                
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
        j = (l*l + l + m)*N + i;                
        Sr[j] = Ylmr; // real part           
        Si[j] = Ylmi; // imag part                                    
                
        m = 1;
        C = sqrt((2*l+1)*fac[l-m]/(4*pi*fac[l+m]));
        Ylmr = C*cos(phi)*P[1];
        Ylmi = C*sin(phi)*P[1];                
        j = (l*l + l + m)*N + i; // spherical harmonics  functions for m > 0                
        Sr[j] = Ylmr; // real part           
        Si[j] = Ylmi; // imag part     
        
        jm = (l*l + l - m)*N + i; // spherical harmonics  functions for m < 0                
        Sr[jm] = -Sr[j]; // real part                      
        Si[jm] =  Si[j]; // imag part                                                                       
                
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
                            
            // Compute spherical harmonics  function at level l for m = 0, 1, ..., l
            for (m=0; m<=l; m++) {
                C = sqrt((2*l+1)*fac[l-m]/(4*pi*fac[l+m])); 
                Ylmr = C*(cos(m*phi)*P[m]);
                Ylmi = C*(sin(m*phi)*P[m]);                                
                j = (l*l + l + m)*N + i;   
                Sr[j] = Ylmr; // real part            
                Si[j] = Ylmi; // imag part                                                            
            }        
            
            // Compute the spherical harmonics  functions for m < 0 using symmetry properties
            C = -1.0;
            for (int m=1; m<=l; m++)  {                 
                j =  (l*l + l + m)*N + i;  // spherical harmonics  functions for m > 0 
                jm = (l*l + l - m)*N + i;  // spherical harmonics  functions for m < 0                 
                Sr[jm] = C*Sr[j]; // real part
                Si[jm] =-C*Si[j]; // imag part                                                                             
                C = C*(-1.0);        
            }
        }
    }          
}
template void cpuSphericalHarmonics(double*, double*, double*, double*, 
         double*, double*, double, int, int);
template void cpuSphericalHarmonics(float*, float*, float*, float*, 
         float*, float*, float, int, int);

template <typename T> void cpuSphericalHarmonicsWithDeriv(T *Sr, T *Si,
        T *Srx, T *Six, T *Sry, T *Siy, T *Srz, T *Siz, T *x, T *y, T *z, 
        T *P, T *tmp, T *dP, T *dtmp, T *fac, T pi, int L, int N)
{                        
    // L                      : the maximum degree of spherical harmonics
    // N                      : length of x, y, z
    // P   [(L+1)*(L+2)/2]    : temporary storage for the recurrence formula
    // tmp [(L+1)*(L+2)/2]    : temporary storage for the recurrence formula
    // dP   [(L+1)*(L+2)/2]   : temporary storage for the recurrence formula
    // dtmp [(L+1)*(L+2)/2]   : temporary storage for the recurrence formula
    // fac                    : factorial look-up table
    // Srx  [N*K*(L+1)*(L+1)] : x-derivative of real spherical harmonics  functions
    // Six  [N*K*(L+1)*(L+1)] : x-derivative of imag spherical harmonics  functions
    // Sry  [N*K*(L+1)*(L+1)] : y-derivative of real spherical harmonics  functions
    // Siy  [N*K*(L+1)*(L+1)] : y-derivative of imag spherical harmonics  functions
    // Srz  [N*K*(L+1)*(L+1)] : z-derivative of real spherical harmonics  functions
    // Siz  [N*K*(L+1)*(L+1)] : z-derivative of imag spherical harmonics  functions
        
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
        
    //#pragma omp parallel for
    for (int i=0; i<N; i++) { // loop over each neighbor atom                               
        int j, jm;
        
        // Cartersian -> Spherical
        T r = sqrt(x[i]*x[i] + y[i]*y[i] + z[i]*z[i]);
        T the = acos(z[i]/r);
        T phi = atan2(y[i],x[i]);
                
        // partial derivatives of spherical coordinates
        T r2 = r*r;
        T rxy = sqrt(x[i]*x[i] + y[i]*y[i]);
        T rxy2 = rxy*rxy;
        T rr2 = rxy*r2;
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
            
        // derivatives of S_{klm}(x,y,z) = g_{lk}(x,y,z)*Y_{lm}(x,y,z) for l = 0
        j = i;                
        Sr[j] = Ylmr;  // real part                   
        Si[j] = Ylmi;  // imag part                   
        Srx[j] = Ylmrx;
        Sry[j] = Ylmry;
        Srz[j] = Ylmrz;
        Six[j] = Ylmix;
        Siy[j] = Ylmiy;
        Siz[j] = Ylmiz;                    
                        
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
        j = (l*l + l + m)*N + i;      
        Sr[j] = Ylmr; // real part           
        Si[j] = Ylmi; // imag part                                                
        Srx[j] = Ylmrx;
        Sry[j] = Ylmry;
        Srz[j] = Ylmrz;
        Six[j] = Ylmix;
        Siy[j] = Ylmiy;
        Siz[j] = Ylmiz;            
        
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
        j = (l*l + l + m)*N + i;          
        Sr[j] = Ylmr; // real part           
        Si[j] = Ylmi; // imag part     
        Srx[j] = Ylmrx;
        Sry[j] = Ylmry;
        Srz[j] = Ylmrz;
        Six[j] = Ylmix;
        Siy[j] = Ylmiy;
        Siz[j] = Ylmiz;            

        // derivatives of Y_{lm}(x,y,z) for l = 1 and m = -1
        jm = (l*l + l - m)*N + i;       
        Sr[jm] = -Sr[j]; // real part                      
        Si[jm] =  Si[j]; // imag part                                                                                   
        Srx[jm] = -Srx[j];
        Sry[jm] = -Sry[j];
        Srz[jm] = -Srz[j];
        Six[jm] = Six[j];
        Siy[jm] = Siy[j];
        Siz[jm] = Siz[j];                    
                
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
                            
            // Compute spherical harmonics  functions and their derivatives at level l for m = 0, 1, ..., l
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

                // derivatives of Y_{lm}(x,y,z) for l >= 2 and m = 0,1,...,l
                j = (l*l + l + m)*N + i;      
                Sr[j] = Ylmr; // real part           
                Si[j] = Ylmi; // imag part     
                Srx[j] = Ylmrx;
                Sry[j] = Ylmry;
                Srz[j] = Ylmrz;
                Six[j] = Ylmix;
                Siy[j] = Ylmiy;
                Siz[j] = Ylmiz;                            
            }
            
            // derivatives of S_{klm}(x,y,z) = g_{lk}(x,y,z)*Y_{lm}(x,y,z) for l >= 2 and m = -1,-2,...,-l            
            C = -1.0;
            for (int m=0; m<l; m++) {  
                j  = (l*l + l + m + 1)*N + i;  // spherical harmonics  functions for m > 0 
                jm = (l*l + l - m - 1)*N + i;  // spherical harmonics  functions for m < 0                                     
                Sr[jm] = C*Sr[j]; // real part
                Si[jm] =-C*Si[j]; // imag part                                                                                                                             
                Srx[jm] = C*Srx[j]; 
                Six[jm] =-C*Six[j];
                Sry[jm] = C*Sry[j]; 
                Siy[jm] =-C*Siy[j];
                Srz[jm] = C*Srz[j]; 
                Siz[jm] =-C*Siz[j];                                    
                C = C*(-1.0);        
            }
        }
    }                        
}
template void cpuSphericalHarmonicsWithDeriv(double*, double*, double*, double*, double*, double*, double*, double*, 
        double*, double*, double*, double*, double*, double*, double*, double*, double, int, int);
template void cpuSphericalHarmonicsWithDeriv(float*, float*, float*, float*, float*, float*, float*, float*, 
        float*, float*, float*, float*, float*, float*, float*, float*, float, int, int);


template <typename T> void cpuSphericalHarmonicsWithDeriv(T *Sr, T *Si, T *Srx, T *Six, T *Sry, T *Siy, 
        T *Srz, T *Siz, T *xij, T *P, T *tmp, T *dP, T *dtmp, T *fac, T pi, int L, int N)
{                                
    //#pragma omp parallel for
    for (int i=0; i<N; i++) { // loop over each neighbor atom                               
        int j, jm;
        
        T x = xij[i*3];
        T y = xij[i*3+1];
        T z = xij[i*3+2];        
        // Cartersian -> Spherical
        T r = sqrt(x*x + y*y + z*z);
        T the = acos(z/r);
        T phi = atan2(y,x);
                
        // partial derivatives of spherical coordinates
        T r2 = r*r;
        T rxy = sqrt(x*x + y*y);
        T rxy2 = rxy*rxy;
        T rr2 = rxy*r2;
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
            
        // derivatives of Y_{lm}(x,y,z) for l = 0
        j = i;                
        Sr[j] = Ylmr;  // real part                   
        Si[j] = Ylmi;  // imag part                   
        Srx[j] = Ylmrx;
        Sry[j] = Ylmry;
        Srz[j] = Ylmrz;
        Six[j] = Ylmix;
        Siy[j] = Ylmiy;
        Siz[j] = Ylmiz;                    
                        
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
        j = (l*l + l + m)*N + i;      
        Sr[j] = Ylmr; // real part           
        Si[j] = Ylmi; // imag part                                                
        Srx[j] = Ylmrx;
        Sry[j] = Ylmry;
        Srz[j] = Ylmrz;
        Six[j] = Ylmix;
        Siy[j] = Ylmiy;
        Siz[j] = Ylmiz;            
        
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
        j = (l*l + l + m)*N + i;          
        Sr[j] = Ylmr; // real part           
        Si[j] = Ylmi; // imag part     
        Srx[j] = Ylmrx;
        Sry[j] = Ylmry;
        Srz[j] = Ylmrz;
        Six[j] = Ylmix;
        Siy[j] = Ylmiy;
        Siz[j] = Ylmiz;            

        // derivatives of Y_{lm}(x,y,z) for l = 1 and m = -1
        jm = (l*l + l - m)*N + i;       
        Sr[jm] = -Sr[j]; // real part                      
        Si[jm] =  Si[j]; // imag part                                                                                   
        Srx[jm] = -Srx[j];
        Sry[jm] = -Sry[j];
        Srz[jm] = -Srz[j];
        Six[jm] = Six[j];
        Siy[jm] = Siy[j];
        Siz[jm] = Siz[j];                    
                
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
                            
            // Compute spherical harmonics  functions and their derivatives at level l for m = 0, 1, ..., l
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

                // derivatives of Y_{lm}(x,y,z) for l >= 2 and m = 0,1,...,l
                j = (l*l + l + m)*N + i;      
                Sr[j] = Ylmr; // real part           
                Si[j] = Ylmi; // imag part     
                Srx[j] = Ylmrx;
                Sry[j] = Ylmry;
                Srz[j] = Ylmrz;
                Six[j] = Ylmix;
                Siy[j] = Ylmiy;
                Siz[j] = Ylmiz;                            
            }
            
            // derivatives of Y_{lm}(x,y,z) for l >= 2 and m = -1,-2,...,-l            
            C = -1.0;
            for (int m=0; m<l; m++) {  
                j  = (l*l + l + m + 1)*N + i;  // spherical harmonics  functions for m > 0 
                jm = (l*l + l - m - 1)*N + i;  // spherical harmonics  functions for m < 0                                     
                Sr[jm] = C*Sr[j]; // real part
                Si[jm] =-C*Si[j]; // imag part                                                                                                                             
                Srx[jm] = C*Srx[j]; 
                Six[jm] =-C*Six[j];
                Sry[jm] = C*Sry[j]; 
                Siy[jm] =-C*Siy[j];
                Srz[jm] = C*Srz[j]; 
                Siz[jm] =-C*Siz[j];                                    
                C = C*(-1.0);        
            }
        }
    }                        
}
template void cpuSphericalHarmonicsWithDeriv(double*, double*, double*, double*, double*, double*, 
        double*, double*, double*, double*, double*, double*, double*, double*, double, int, int);
template void cpuSphericalHarmonicsWithDeriv(float*, float*, float*, float*, float*, float*, 
        float*, float*, float*, float*, float*, float*, float*, float*, float, int, int);


#endif


