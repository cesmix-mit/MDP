#ifndef __CPUPRECOMPUTE
#define __CPUPRECOMPUTE

void cpuGetIndk(int *indk, int K)
{               
    int N = K*(K+1)/2;
    int count = 0;
    for (int k1=0; k1<K; k1++)
        for (int k2=0; k2<K; k2++) {
            indk[count] = k2;
            indk[N+count] = k1;
            count += 1;
        }
}

// int cpuGetM(int *indl, int N)
// {               
//     int M = 0;
//     for (int n=0; n<N; n++) {
//         int l2 = indl[n];
//         int l1 = indl[N+n];
//         int l = indl[2*N+n];
//         int lm = max(l1-l2,l2-l1);
//         for (int m=-l; m<=l; m++) 
//             for (int m1=-l1; m1<=l1; m1++) 
//                 for (int m2=-l2; m2<=l2; m2++) 
//                     if !((m1+m2-m != 0) || (l < lm) || (l > l1+l2))                        
//                         M = M + 1;                    
//     }
// 
//     return M;
// }

double clebschgordan(int j1, int m1, int j2, int m2, int j, int m, double *fac)
{               
    double C;
    
    int jm = max(j1-j2,j2-j1);
    
    if ((m1+m2-m != 0) || (j < jm)  || (j > j1+j2)) {
        C = 0.0;
    }
    else {
        int k1 = max(max(0,j2-j-m1),j1-j+m2);
        int k2 = min(min(j1+j2-j,j1-m1),j2+m2);
        int n = k2-k1+1;        
        C = 0.0;
        for (int i=0; i<n; i++) {
            int k = k1 + i;        
            double a = ( k % 2 == 0) ?  1.0 : -1.0;
            C += a/(fac[k]*fac[j1+j2-j-k]*fac[j1-m1-k]*fac[j2+m2-k]*fac[j-j2+m1+k]*fac[j-j1-m2+k]);
        }
        C = C*sqrt((2*j+1)*fac[j+j1-j2]*fac[j+j2-j1]*fac[j1+j2-j]/fac[j1+j2+j+1])*
              sqrt(fac[j+m]*fac[j-m]*fac[j1+m1]*fac[j1-m1]*fac[j2+m2]*fac[j2-m2]);        
    }    
    return C;
}


int cgcoefficients(double **cg, int **indm, int *rowm, int *indl, double *fac, int N)
{                   
    int M = 0;
    for (int n=0; n<N; n++) {
        int l2 = indl[n];
        int l1 = indl[N+n];
        int l = indl[2*N+n];
        int lm = max(l1-l2,l2-l1);
        for (int m=-l; m<=l; m++) 
            for (int m1=-l1; m1<=l1; m1++) 
                for (int m2=-l2; m2<=l2; m2++) 
                    if (!((m1+m2-m != 0) || (l < lm) || (l > l1+l2)))                        
                        M = M + 1;                    
    }
    *indm = (int *) malloc(3*M*sizeof(int));  
    *cg = (double *) malloc(M*sizeof(double));  
    
    int nc = 0;
    rowm[0] = 0;
    for (int n=0; n<N; n++) {
        int l2 = indl[n];
        int l1 = indl[N+n];
        int l = indl[2*N+n];
        int lm = max(l1-l2,l2-l1);
        int count = 0;
        for (int m=-l; m<=l; m++) 
            for (int m1=-l1; m1<=l1; m1++) 
                for (int m2=-l2; m2<=l2; m2++) 
                    if (!((m1+m2-m != 0) || (l < lm) || (l > l1+l2))) {                    
                        double C = clebschgordan(l1, m1, l2, m2, l, m, fac);   
                        *cg[nc] = C;
                        *indm[nc] = m2;   
                        *indm[M+nc] = m1;   
                        *indm[2*M+nc] = m;   
                        nc += 1; 
                        count = count + 1;
                    }
        rowm[n+1] = count;
    }
    
    for (int n=2; n<=N; n++)
        rowm[n] = rowm[n-2] + rowm[n-1];
    
    return M;
}

void cpuSphericalHarmonicsBispectrum(double *b, double *Sr, double *Si, double *fac, int L)
{               
    int N = L + 1;
    for (int l=0; l<N; l++)     
        for (int l1=0; l1<N; l1++)     
            for (int l2=0; l2<N; l2++)  {                         
                double tmp = 0;
                for (int m=-l; m<=l; m++) 
                    for (int m1=-l1; m1<=l1; m1++) 
                        for (int m2=-l2; m2<=l2; m2++) {
                            double cg = clebschgordan(l2,m2,l1,m1,l,m,fac);
                            if (fabs(cg)>0) {
                                int mm = -m;
                                double a = ( mm % 2 == 0) ?  1.0 : -1.0;
                                double a1, b1, a2, b2, a3, b3;
                                if (m>=0) {
                                    a1 = Sr[l*(l+1)/2+m];
                                    b1 = Si[l*(l+1)/2+m];
                                }
                                else {                                                                        
                                    a1 =-a*Sr[l*(l+1)/2+mm];
                                    b1 = a*Si[l*(l+1)/2+mm];
                                }
                                mm = -m1;
                                a = ( mm % 2 == 0) ?  1.0 : -1.0;
                                if (m1>=0) {
                                    a2 = Sr[l1*(l1+1)/2+m1];
                                    b2 = Si[l1*(l1+1)/2+m1];
                                }
                                else {
                                    a2 = -a*Sr[l1*(l1+1)/2+m1];
                                    b2 =  a*Si[l1*(l1+1)/2+m1];
                                }
                                mm = -m2;
                                a = ( mm % 2 == 0) ?  1.0 : -1.0;
                                if (m2>=0) {
                                    a3 = Sr[l2*(l2+1)/2+m2];
                                    b3 = Si[l2*(l2+1)/2+m2];      
                                }
                                else {
                                    a3 = -a*Sr[l2*(l2+1)/2+m2];
                                    b3 = -a*Si[l2*(l2+1)/2+m2];      
                                }                                
                                tmp = tmp + cg*(a1*a2*a3 + a2*b1*b3 + a3*b1*b2 - a1*b2*b3);                                                 
                            }
                        }
                b[l*N*N+l1*N+l2] = tmp;
            }
}

void cpuSphericalHarmonicsBispectrumIndex(int **indl, double *b, int L)
{               
    int N = L+1;
    int inc = 0;
    for (int l=0; l<N; l++)     
        for (int l1=0; l1<N; l1++)     
            for (int l2=0; l2<N; l2++)  {                         
                double tm = b[l*N*N+l1*N+l2];            
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
    
    int M = inc;
    *indl = (int *) malloc(3*M*sizeof(int));      
        
    inc = 0;
    for (int l=0; l<N; l++)     
        for (int l1=0; l1<N; l1++)     
            for (int l2=0; l2<N; l2++)  {                         
                double tm = b[l*N*N+l1*N+l2];            
                if ((fabs(tm)>1e-10) && (l2<=l1)) {
                    if (l1==l2) {                        
                        *indl[inc] = l2;
                        *indl[M+inc] = l1;
                        *indl[2*M+inc] = l;
                        inc = inc + 1;
                    }
                    else {
                        if ((l2<l1) && (l2<l) && (l1<l)) {
                            *indl[inc] = l2;
                            *indl[M+inc] = l1;
                            *indl[2*M+inc] = l;
                            inc = inc + 1;                
                        }
                    }
                }
            }
        
}


#endif


