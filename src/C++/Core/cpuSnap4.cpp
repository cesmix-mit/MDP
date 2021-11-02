/***************************************************************************                               
                    Molecular Dynamics Potentials (MDP)
                           CESMIX-MIT Project  
 
 Contributing authors: Ngoc-Cuong Nguyen (cuongng@mit.edu, exapde@gmail.com)
 ***************************************************************************/

#ifndef MDP_CPUSNAP4
#define MDP_CPUSNAP4

template <typename T> void cpuSnapComputeUi(T *Utotr, T *Utoti, T *rootpqarray, T *rij, T *wjelem, T *radelem, 
        T rmin0, T rfac0, T rcutfac, int *idxu_block, int *map, int *aii, int *ti, int *tj, 
        int twojmax, int idxu_max, int inum, int ijnum, int switchflag, int chemflag)                
{    
  for(int ij=0; ij<ijnum; ij++) {        
    T x = rij[ij*3+0];
    T y = rij[ij*3+1];
    T z = rij[ij*3+2];    
    T r = sqrt(x * x + y * y + z * z);

    T rcutij = (radelem[ti[ij]]+radelem[tj[ij]])*rcutfac; //(radelem[type[ii]]+radelem[type[jj]])*rcutfac;
    //T rscale0 = rfac0 * M_PI / (rcutij - rmin0);
    //T theta0 = (r - rmin0) * rscale0;
    T z0 = r / tan((r - rmin0) * rfac0 * M_PI / (rcutij - rmin0));                
            
    T sfac = 0.0;
    if (switchflag == 0) 
        sfac = 1.0;    
    else if (switchflag == 1) {
        if (r <= rmin0) {
            sfac = 1.0;
        }
        else if(r > rcutij) {
            sfac = 1.0;
        }
        else {
            T rcutfac = M_PI / (rcutij - rmin0);
            sfac =  0.5 * (cos((r - rmin0) * rcutfac) + 1.0);   
        }
    } 
    sfac *= wjelem[tj[ij]];
    
    T a_r, a_i, b_r, b_i;
    T rootpq;

    //T r0inv;    
    rcutij = 1.0 / sqrt(r * r + z0 * z0);
    a_r = rcutij * z0;
    a_i = -rcutij * z;
    b_r = rcutij * y;
    b_i = -rcutij * x;

    // 2Jmax = 10
    T Pr[11], Pi[11], Qr[9], Qi[9];
    Pr[0] = 1.0;
    Pi[0] = 0.0;    
    
    int jdim = twojmax + 1;
    int njelem = (chemflag==0) ? 0 : (inum*idxu_max*map[tj[ij]]);    
    int i = aii[ij] + njelem;                        
    Utotr[i] += sfac; // atomic add   
    
    int mb = 0;    
    for (int j = 1; j <= twojmax; j++) {        
        // fill in left side of matrix layer from previous layer
        int ma = 0;
        // x y z z0 
        // T p_r, p_i; // -> x, y
        // T u_r = Pr[ma]; // -> z
        // T u_i = Pi[ma]; // -> z0
        z = Pr[ma];
        z0 = Pi[ma];
        rootpq = rootpqarray[(j - ma)*jdim + (j - mb)];
        Pr[ma] = rootpq * (a_r * z + a_i * z0);
        Pi[ma] = rootpq * (a_r * z0 - a_i * z);            
        for (ma = 1; ma < j; ma++) {
            x = Pr[ma];
            y = Pi[ma];
            rootpq = rootpqarray[(j - ma)*jdim + (j - mb)];
            rcutij = rootpqarray[ma*jdim + (j - mb)];
            Pr[ma] = rootpq * (a_r * x + a_i * y) -rcutij * (b_r * z + b_i * z0);
            Pi[ma] = rootpq * (a_r * y - a_i * x) -rcutij * (b_r * z0 - b_i * z);
            z = x;
            z0 = y;
        }
        ma = j;
        rcutij = rootpqarray[ma*jdim + (j - mb)];
        Pr[ma] = -rcutij * (b_r * z + b_i * z0);
        Pi[ma] = -rcutij * (b_r * z0 - b_i * z);                        
                                
        if (j==(2*mb+1)) { // store Qr, Qi, for the next mb level
            int mapar = (((j+1)/2)%2 == 0) ? -1 : 1;
            for (int ma = 0; ma <= j; ma++) {     
                if (mapar == 1) {                    
                    Qr[j-ma] = Pr[ma];
                    Qi[j-ma] = -Pi[ma];
                } else {
                    Qr[j-ma] = -Pr[ma];
                    Qi[j-ma] =  Pi[ma];
                }
                mapar = -mapar;
            }                                                
        }
        
        int k =  1 + (j+1)*mb;
        for (int ma = 2; ma <= j; ma++)
            k += ma*ma;                    
        for (int ma = 0; ma <= j; ma++) {
            int in = i + inum*k;                
            Utotr[in] += sfac*Pr[ma]; // atomic add   
            Utoti[in] += sfac*Pi[ma]; // atomic add                       
            k += 1;
        }                   
    }
    
    for (mb = 1; 2*mb <= twojmax; mb++) {     
        for (int ma = 0; ma < 2*mb; ma++) {                      
            Pr[ma] = Qr[ma];
            Pi[ma] = Qi[ma];
        }                
        for (int j = 2*mb; j <= twojmax; j++) { 
            int ma = 0;
            //T p_r, p_i;
            //T u_r = Pr[ma];
            //T u_i = Pi[ma];
            z = Pr[ma];
            z0 = Pi[ma];
            rootpq = rootpqarray[(j - ma)*jdim + (j - mb)];
            Pr[ma] = rootpq * (a_r * z + a_i * z0);
            Pi[ma] = rootpq * (a_r * z0 - a_i * z);            
            for (ma = 1; ma < j; ma++) {
                x = Pr[ma];
                y = Pi[ma];
                rootpq = rootpqarray[(j - ma)*jdim + (j - mb)];
                rcutij = rootpqarray[ma*jdim + (j - mb)];
                Pr[ma] = rootpq * (a_r * x + a_i * y) -rcutij * (b_r * z + b_i * z0);
                Pi[ma] = rootpq * (a_r * y - a_i * x) -rcutij * (b_r * z0 - b_i * z);
                z = x;
                z0 = y;
            }
            ma = j;
            rcutij = rootpqarray[ma*jdim + (j - mb)];
            Pr[ma] = -rcutij * (b_r * z + b_i * z0);
            Pi[ma] = -rcutij * (b_r * z0 - b_i * z);       
            
            if (j==(2*mb)) {
                int mapar = 1;
                for (int ma = 0; ma <= j/2; ma++) {
                    if (mapar == 1) {                    
                        Pr[j/2+ma] = Pr[j/2-ma];
                        Pi[j/2+ma] = -Pi[j/2-ma];
                    } else {
                        Pr[j/2+ma] = -Pr[j/2-ma];
                        Pi[j/2+ma] = Pi[j/2-ma];
                    }
                    mapar = -mapar;        
                }                                                        
            }
            
            // store Qr, Qi, for the next mb level
            if (j==(2*mb+1)) {
                int mapar = (((j+1)/2)%2 == 0) ? -1 : 1;
                for (int ma = 0; ma <= j; ma++) {     
                    if (mapar == 1) {                    
                        Qr[j-ma] = Pr[ma];
                        Qi[j-ma] = -Pi[ma];
                    } else {
                        Qr[j-ma] = -Pr[ma];
                        Qi[j-ma] =  Pi[ma];
                    }
                    mapar = -mapar;
                }                                                
            }
            
            int k =  1 + (j+1)*mb;
            for (int ma = 2; ma <= j; ma++)
                k += ma*ma;            
            for (int ma = 0; ma <= j; ma++) {
                int in = i + inum*k;                
                Utotr[in] += sfac*Pr[ma]; // atomic add   
                Utoti[in] += sfac*Pi[ma]; // atomic add                       
                k += 1; 
            }                                                           
        }
    }        
  }
};
template void cpuSnapComputeUi(double*, double*, double*, double*, double*, double*, double, double, double, 
        int*, int*, int*, int*, int*, int, int, int, int, int, int);
template void cpuSnapComputeUi(float*, float*, float*, float*, float*, float*, float, float, float, 
        int*, int*, int*, int*, int*, int, int, int, int, int, int);

template <typename T> void cpuAddWself2Ui(T *Utotr, T *Utoti, T wself, int *idxu_block, 
        int *type, int *map, int *ai, int wselfall_flag, int chemflag, int idxu_max, int nelements, 
         int twojmax, int inum)
{
    int N1 = inum;
    int N2 = N1*(twojmax+1);
    int N3 = N2*nelements;                                
    
    for (int idx=0; idx < N3; idx++) {
        int l = idx%N2;  // inum*(twojmax+1)
        int ii = l%N1;    // inum
        int j = (l-ii)/N1; // (twojmax+1)
        int jelem = (idx-l)/N2; // nelements   
        int ielem = (chemflag) ? map[type[ii]]: 0;                
        int nmax = ii + inum*idxu_max*jelem;
        
        int jju = idxu_block[j];                
        for (int mb = 0; mb <= j; mb++) {
            for (int ma = 0; ma <= j; ma++) {                
                if (jelem == ielem || wselfall_flag)
                    if (ma==mb)                        
                        Utotr[inum*jju + nmax] += wself;                                     
                jju++;                                
            }
        }
        
        // copy left side to right side with inversion symmetry VMK 4.4(2)
        // u[ma-j][mb-j] = (-1)^(ma-mb)*Conj([u[ma][mb])
        
        jju = idxu_block[j];        
        int jjup = jju+(j+1)*(j+1)-1;
        int mbpar = 1;
        for (int mb = 0; 2*mb < j; mb++) {
            int mapar = mbpar;
            for (int ma = 0; ma <= j; ma++) {
                int njju =  inum*jju + nmax;
                int njjup = inum*jjup + nmax;
                if (mapar == 1) {
                    Utotr[njjup] = Utotr[njju];
                    Utoti[njjup] = -Utoti[njju];
                } else {
                    Utotr[njjup] = -Utotr[njju];
                    Utoti[njjup] =  Utoti[njju];
                }
                mapar = -mapar;
                jju++;
                jjup--;
            }
            mbpar = -mbpar;
        }        
    }                    
};
template void cpuAddWself2Ui(double*, double*, double, int*, int*, 
        int*, int*, int, int, int, int, int, int);
template void cpuAddWself2Ui(float*, float*, float, int*, int*, 
        int*, int*, int, int, int, int, int, int);

template <typename T> void cpuSnapComputeEi(T *eatom, T *Utotr, T *Utoti, T *cglist, 
        T *bzero, T *coeffelem, int *ilist, int *map, int *type, int *idxb, int *idxcg_block, int *idxu_block, 
        int twojmax, int idxb_max, int idxu_max, int nelements, int ncoeffall, int bnorm_flag, 
        int bzero_flag, int wselfall_flag, int quadraticflag, int inum)
{    
    int nelemsq = nelements*nelements;    
    int nu_max = idxu_max*inum;
    int nb_max = idxb_max*inum;    
    int jdim = twojmax+1;
    int N2 = inum*idxb_max;
    int N3 = N2*nelements*nelemsq;    
    for (int idx=0; idx < N3; idx++) {
        int l = idx%N2;        
        int ii = l%inum;
        int jjb = (l-ii)/inum;
        int jelem = (idx-l)/N2;                    
        int k = jelem%nelemsq;      
        int elem3 = k%nelements;
        int elem2 = (k-elem3)/nelements;
        int elem1 = (jelem-k)/nelemsq;    
        //int itriple = elem3 + nelements*elem2 + nelemsq*elem1;    
        int idouble = elem2 + nelements*elem1;  
        const int j1 = idxb[jjb*3 + 0];
        const int j2 = idxb[jjb*3 + 1];
        const int j = idxb[jjb*3 + 2];  
        const T *cgblock = cglist + idxcg_block[j + j2*jdim + j1*jdim*jdim];

        int itype = type[ii]; // element type of atom i
        int ielem = map[itype];  // index of that element type        
        int icoeff = jjb + idxb_max*(jelem);
        int itriple = 1+icoeff+ielem*ncoeffall;
        
        int jju = idxu_block[j];
        int ii1 = ii + inum*idxu_max*elem1;
        int ii2 = ii + inum*idxu_max*elem2;        
        int idu, jju1, jju2, icga, icgb;        
        int ia, ib, ma, mb, ma1min, ma2max, na, mb1min, mb2max, nb, ma1, ma2;
        T betaj, ztmp_r, ztmp_i, sumzu = 0.0;
        for (mb = 0; 2 * mb < j; mb++)
            for (ma = 0; ma <= j; ma++) {
                idu = ii + inum*jju + nu_max*elem3;        
                // calculate Z_{j1,j2,j}^{elem1,elem2}(ma,mb)   
                ma1min = max(0, (2 * ma - j - j2 + j1) / 2);
                ma2max = (2 * ma - j - (2 * ma1min - j1) + j2) / 2;
                na = min(j1, (2 * ma - j + j2 + j1) / 2) - ma1min + 1;
                mb1min = max(0, (2 * mb - j - j2 + j1) / 2);
                mb2max = (2 * mb - j - (2 * mb1min - j1) + j2) / 2;
                nb = min(j1, (2 * mb - j + j2 + j1) / 2) - mb1min + 1;                
                jju1 = idxu_block[j1] + (j1 + 1) * mb1min;
                jju2 = idxu_block[j2] + (j2 + 1) * mb2max;
                icgb = mb1min * (j2 + 1) + mb2max;
                ztmp_r = 0.0;
                ztmp_i = 0.0;
                for (ib = 0; ib < nb; ib++) {
                    ma1 = ma1min;
                    ma2 = ma2max;
                    icga = ma1min * (j2 + 1) + ma2max;
                    for (ia = 0; ia < na; ia++) {
                        ztmp_r += cgblock[icgb]*cgblock[icga] * (Utotr[ii1+inum*(jju1+ma1)] * Utotr[ii2+inum*(jju2+ma2)] - Utoti[ii1+inum*(jju1+ma1)] * Utoti[ii2+inum*(jju2+ma2)]);
                        ztmp_i += cgblock[icgb]*cgblock[icga] * (Utotr[ii1+inum*(jju1+ma1)] * Utoti[ii2+inum*(jju2+ma2)] + Utoti[ii1+inum*(jju1+ma1)] * Utotr[ii2+inum*(jju2+ma2)]);                    
                      ma1++;
                      ma2--;
                      icga += j2;
                    } // end loop over ia
                    jju1 += j1 + 1;
                    jju2 -= j2 + 1;
                    icgb += j2;
                } // end loop over ib
                if (bnorm_flag) {
                  ztmp_i /= j+1;
                  ztmp_r /= j+1;
                }            
                sumzu += Utotr[idu] * ztmp_r + Utoti[idu] * ztmp_i;
                //sumzu += Utotr[idu] * zlist_r[idz] + Utoti[idu] * zlist_i[idz];                                              
                jju++;
            } // end loop over ma, mb

        // For j even, handle middle column
        if (j % 2 == 0) {
            mb = j / 2;
            for (ma = 0; ma < mb; ma++) {
                idu = ii + inum*jju + nu_max*elem3;        
                ma1min = max(0, (2 * ma - j - j2 + j1) / 2);
                ma2max = (2 * ma - j - (2 * ma1min - j1) + j2) / 2;
                na = min(j1, (2 * ma - j + j2 + j1) / 2) - ma1min + 1;
                mb1min = max(0, (2 * mb - j - j2 + j1) / 2);
                mb2max = (2 * mb - j - (2 * mb1min - j1) + j2) / 2;
                nb = min(j1, (2 * mb - j + j2 + j1) / 2) - mb1min + 1;                
                jju1 = idxu_block[j1] + (j1 + 1) * mb1min;
                jju2 = idxu_block[j2] + (j2 + 1) * mb2max;
                icgb = mb1min * (j2 + 1) + mb2max;
                ztmp_r = 0.0;
                ztmp_i = 0.0;
                for (ib = 0; ib < nb; ib++) {
                    ma1 = ma1min;
                    ma2 = ma2max;
                    icga = ma1min * (j2 + 1) + ma2max;
                    for (ia = 0; ia < na; ia++) {
                        ztmp_r += cgblock[icgb]*cgblock[icga] * (Utotr[ii1+inum*(jju1+ma1)] * Utotr[ii2+inum*(jju2+ma2)] - Utoti[ii1+inum*(jju1+ma1)] * Utoti[ii2+inum*(jju2+ma2)]);
                        ztmp_i += cgblock[icgb]*cgblock[icga] * (Utotr[ii1+inum*(jju1+ma1)] * Utoti[ii2+inum*(jju2+ma2)] + Utoti[ii1+inum*(jju1+ma1)] * Utotr[ii2+inum*(jju2+ma2)]);                    
                      ma1++;
                      ma2--;
                      icga += j2;
                    } // end loop over ia
                    jju1 += j1 + 1;
                    jju2 -= j2 + 1;
                    icgb += j2;
                } // end loop over ib
                if (bnorm_flag) {
                  ztmp_i /= j+1;
                  ztmp_r /= j+1;
                }            
                sumzu += Utotr[idu] * ztmp_r + Utoti[idu] * ztmp_i;                
                //sumzu += Utotr[idu] * zlist_r[idz] + Utoti[idu] * zlist_i[idz];                                
                jju++;
            }
            ma = mb;
            idu = ii + inum*jju + nu_max*elem3;        
            ma1min = max(0, (2 * ma - j - j2 + j1) / 2);
            ma2max = (2 * ma - j - (2 * ma1min - j1) + j2) / 2;
            na = min(j1, (2 * ma - j + j2 + j1) / 2) - ma1min + 1;
            mb1min = max(0, (2 * mb - j - j2 + j1) / 2);
            mb2max = (2 * mb - j - (2 * mb1min - j1) + j2) / 2;
            nb = min(j1, (2 * mb - j + j2 + j1) / 2) - mb1min + 1;                
            jju1 = idxu_block[j1] + (j1 + 1) * mb1min;
            jju2 = idxu_block[j2] + (j2 + 1) * mb2max;
            icgb = mb1min * (j2 + 1) + mb2max;
            ztmp_r = 0.0;
            ztmp_i = 0.0;
            for (ib = 0; ib < nb; ib++) {
                ma1 = ma1min;
                ma2 = ma2max;
                icga = ma1min * (j2 + 1) + ma2max;
                for (ia = 0; ia < na; ia++) {
                    ztmp_r += cgblock[icgb]*cgblock[icga] * (Utotr[ii1+inum*(jju1+ma1)] * Utotr[ii2+inum*(jju2+ma2)] - Utoti[ii1+inum*(jju1+ma1)] * Utoti[ii2+inum*(jju2+ma2)]);
                    ztmp_i += cgblock[icgb]*cgblock[icga] * (Utotr[ii1+inum*(jju1+ma1)] * Utoti[ii2+inum*(jju2+ma2)] + Utoti[ii1+inum*(jju1+ma1)] * Utotr[ii2+inum*(jju2+ma2)]);                    
                  ma1++;
                  ma2--;
                  icga += j2;
                } // end loop over ia
                jju1 += j1 + 1;
                jju2 -= j2 + 1;
                icgb += j2;
            } // end loop over ib
            if (bnorm_flag) {
              ztmp_i /= j+1;
              ztmp_r /= j+1;
            }            
            sumzu += 0.5*(Utotr[idu] * ztmp_r + Utoti[idu] * ztmp_i);            
            //sumzu += 0.5 * (Utotr[idu] * zlist_r[idz] + Utoti[idu] * zlist_i[idz]);            
        } // end if jeven
        
        sumzu *= 2.0;                
        
        if (bzero_flag) {
          if (!wselfall_flag) {
            if (elem1 == elem2 && elem1 == elem3) {
              sumzu -= bzero[j];
            }
          } 
          else {
            sumzu -= bzero[j];
          }
        }
        
        //blist[idx] = sumzu;                                
        
        if (icoeff==0)
            eatom[ilist[ii]] += coeffelem[ielem*ncoeffall];
        eatom[ilist[ii]] += coeffelem[itriple]*sumzu;                           
    }
}
template void cpuSnapComputeEi(double*, double*, double*, double*, double*, double*,
        int*, int*, int*, int*, int*, int*, int, int, int, int, int, int, int, int, int, int);
template void cpuSnapComputeEi(float*, float*, float*, float*, float*, float*, 
        int*, int*, int*, int*, int*, int*, int, int, int, int, int, int, int, int, int, int);

template <typename T> void cpuSnapComputeBi(T *blist, T *Utotr, T *Utoti, T *cglist, T *bzero, 
        int *idxb, int *idxcg_block, int *idxu_block, int *idxz_block, int twojmax, int idxb_max, int idxu_max, 
        int idxz_max, int nelements, int bnorm_flag, int bzero_flag, int wselfall_flag, int inum)
{    
    int nelemsq = nelements*nelements;
    int nz_max = idxz_max*inum;
    int nu_max = idxu_max*inum;
    int nb_max = idxb_max*inum;    
    int jdim = twojmax+1;
    int N2 = inum*idxb_max;
    int N3 = N2*nelements*nelemsq;    
    for (int idx=0; idx < N3; idx++) {
        int l = idx%N2;        
        int ii = l%inum;
        int jjb = (l-ii)/inum;
        int jelem = (idx-l)/N2;                    
        int k = jelem%nelemsq;      
        int elem3 = k%nelements;
        int elem2 = (k-elem3)/nelements;
        int elem1 = (jelem-k)/nelemsq;    
        //int itriple = elem3 + nelements*elem2 + nelemsq*elem1;    
        int idouble = elem2 + nelements*elem1;  
        const int j1 = idxb[jjb*3 + 0];
        const int j2 = idxb[jjb*3 + 1];
        const int j = idxb[jjb*3 + 2];  
        const T *cgblock = cglist + idxcg_block[j + j2*jdim + j1*jdim*jdim];

        int jjz = idxz_block[j + j2*jdim + j1*jdim*jdim];
        int jju = idxu_block[j];
        int ii1 = ii + inum*idxu_max*elem1;
        int ii2 = ii + inum*idxu_max*elem2;        
        int idu, idz, jju1, jju2, icga, icgb;        
        int ia, ib, ma, mb, ma1min, ma2max, na, mb1min, mb2max, nb, ma1, ma2;
        T ztmp_r, ztmp_i, sumzu = 0.0;
        for (mb = 0; 2 * mb < j; mb++)
            for (ma = 0; ma <= j; ma++) {
                idu = ii + inum*jju + nu_max*elem3;        
                idz = ii + inum*jjz + nz_max*idouble;        
                // calculate Z_{j1,j2,j}^{elem1,elem2}(ma,mb)   
                ma1min = max(0, (2 * ma - j - j2 + j1) / 2);
                ma2max = (2 * ma - j - (2 * ma1min - j1) + j2) / 2;
                na = min(j1, (2 * ma - j + j2 + j1) / 2) - ma1min + 1;
                mb1min = max(0, (2 * mb - j - j2 + j1) / 2);
                mb2max = (2 * mb - j - (2 * mb1min - j1) + j2) / 2;
                nb = min(j1, (2 * mb - j + j2 + j1) / 2) - mb1min + 1;                
                jju1 = idxu_block[j1] + (j1 + 1) * mb1min;
                jju2 = idxu_block[j2] + (j2 + 1) * mb2max;
                icgb = mb1min * (j2 + 1) + mb2max;
                ztmp_r = 0.0;
                ztmp_i = 0.0;
                for (ib = 0; ib < nb; ib++) {
                    ma1 = ma1min;
                    ma2 = ma2max;
                    icga = ma1min * (j2 + 1) + ma2max;
                    for (ia = 0; ia < na; ia++) {
                        ztmp_r += cgblock[icgb]*cgblock[icga] * (Utotr[ii1+inum*(jju1+ma1)] * Utotr[ii2+inum*(jju2+ma2)] - Utoti[ii1+inum*(jju1+ma1)] * Utoti[ii2+inum*(jju2+ma2)]);
                        ztmp_i += cgblock[icgb]*cgblock[icga] * (Utotr[ii1+inum*(jju1+ma1)] * Utoti[ii2+inum*(jju2+ma2)] + Utoti[ii1+inum*(jju1+ma1)] * Utotr[ii2+inum*(jju2+ma2)]);                    
                      ma1++;
                      ma2--;
                      icga += j2;
                    } // end loop over ia
                    jju1 += j1 + 1;
                    jju2 -= j2 + 1;
                    icgb += j2;
                } // end loop over ib
                if (bnorm_flag) {
                  ztmp_i /= j+1;
                  ztmp_r /= j+1;
                }            
                sumzu += Utotr[idu] * ztmp_r + Utoti[idu] * ztmp_i;
                //sumzu += Utotr[idu] * zlist_r[idz] + Utoti[idu] * zlist_i[idz];
                jjz++;
                jju++;
            } // end loop over ma, mb

        // For j even, handle middle column
        if (j % 2 == 0) {
            mb = j / 2;
            for (ma = 0; ma < mb; ma++) {
                idu = ii + inum*jju + nu_max*elem3;        
                idz = ii + inum*jjz + nz_max*idouble;        
                ma1min = max(0, (2 * ma - j - j2 + j1) / 2);
                ma2max = (2 * ma - j - (2 * ma1min - j1) + j2) / 2;
                na = min(j1, (2 * ma - j + j2 + j1) / 2) - ma1min + 1;
                mb1min = max(0, (2 * mb - j - j2 + j1) / 2);
                mb2max = (2 * mb - j - (2 * mb1min - j1) + j2) / 2;
                nb = min(j1, (2 * mb - j + j2 + j1) / 2) - mb1min + 1;                
                jju1 = idxu_block[j1] + (j1 + 1) * mb1min;
                jju2 = idxu_block[j2] + (j2 + 1) * mb2max;
                icgb = mb1min * (j2 + 1) + mb2max;
                ztmp_r = 0.0;
                ztmp_i = 0.0;
                for (ib = 0; ib < nb; ib++) {
                    ma1 = ma1min;
                    ma2 = ma2max;
                    icga = ma1min * (j2 + 1) + ma2max;
                    for (ia = 0; ia < na; ia++) {
                        ztmp_r += cgblock[icgb]*cgblock[icga] * (Utotr[ii1+inum*(jju1+ma1)] * Utotr[ii2+inum*(jju2+ma2)] - Utoti[ii1+inum*(jju1+ma1)] * Utoti[ii2+inum*(jju2+ma2)]);
                        ztmp_i += cgblock[icgb]*cgblock[icga] * (Utotr[ii1+inum*(jju1+ma1)] * Utoti[ii2+inum*(jju2+ma2)] + Utoti[ii1+inum*(jju1+ma1)] * Utotr[ii2+inum*(jju2+ma2)]);                    
                      ma1++;
                      ma2--;
                      icga += j2;
                    } // end loop over ia
                    jju1 += j1 + 1;
                    jju2 -= j2 + 1;
                    icgb += j2;
                } // end loop over ib
                if (bnorm_flag) {
                  ztmp_i /= j+1;
                  ztmp_r /= j+1;
                }            
                sumzu += Utotr[idu] * ztmp_r + Utoti[idu] * ztmp_i;                
                //sumzu += Utotr[idu] * zlist_r[idz] + Utoti[idu] * zlist_i[idz];
                jjz++;
                jju++;
            }
            ma = mb;
            idu = ii + inum*jju + nu_max*elem3;        
            idz = ii + inum*jjz + nz_max*idouble;        
            ma1min = max(0, (2 * ma - j - j2 + j1) / 2);
            ma2max = (2 * ma - j - (2 * ma1min - j1) + j2) / 2;
            na = min(j1, (2 * ma - j + j2 + j1) / 2) - ma1min + 1;
            mb1min = max(0, (2 * mb - j - j2 + j1) / 2);
            mb2max = (2 * mb - j - (2 * mb1min - j1) + j2) / 2;
            nb = min(j1, (2 * mb - j + j2 + j1) / 2) - mb1min + 1;                
            jju1 = idxu_block[j1] + (j1 + 1) * mb1min;
            jju2 = idxu_block[j2] + (j2 + 1) * mb2max;
            icgb = mb1min * (j2 + 1) + mb2max;
            ztmp_r = 0.0;
            ztmp_i = 0.0;
            for (ib = 0; ib < nb; ib++) {
                ma1 = ma1min;
                ma2 = ma2max;
                icga = ma1min * (j2 + 1) + ma2max;
                for (ia = 0; ia < na; ia++) {
                    ztmp_r += cgblock[icgb]*cgblock[icga] * (Utotr[ii1+inum*(jju1+ma1)] * Utotr[ii2+inum*(jju2+ma2)] - Utoti[ii1+inum*(jju1+ma1)] * Utoti[ii2+inum*(jju2+ma2)]);
                    ztmp_i += cgblock[icgb]*cgblock[icga] * (Utotr[ii1+inum*(jju1+ma1)] * Utoti[ii2+inum*(jju2+ma2)] + Utoti[ii1+inum*(jju1+ma1)] * Utotr[ii2+inum*(jju2+ma2)]);                    
                  ma1++;
                  ma2--;
                  icga += j2;
                } // end loop over ia
                jju1 += j1 + 1;
                jju2 -= j2 + 1;
                icgb += j2;
            } // end loop over ib
            if (bnorm_flag) {
              ztmp_i /= j+1;
              ztmp_r /= j+1;
            }            
            sumzu += 0.5*(Utotr[idu] * ztmp_r + Utoti[idu] * ztmp_i);            
            //sumzu += 0.5 * (Utotr[idu] * zlist_r[idz] + Utoti[idu] * zlist_i[idz]);
        } // end if jeven
        
        sumzu *= 2.0;                
        
        if (bzero_flag) {
          if (!wselfall_flag) {
            if (elem1 == elem2 && elem1 == elem3) {
              sumzu -= bzero[j];
            }
          } 
          else {
            sumzu -= bzero[j];
          }
        }
        
        blist[idx] = sumzu;                                        
    }
}
template void cpuSnapComputeBi(double*, double*, double*, double*, double*, 
        int*, int*, int*, int*, int, int, int, int, int, int, int, int, int);
template void cpuSnapComputeBi(float*, float*, float*, float*, float*, 
        int*, int*, int*, int*, int, int, int, int, int, int, int, int, int);

template <typename T> void cpuSnapComputeYi(T *ylist_r, T *ylist_i, T *Utotr, T *Utoti, 
        T *cglist, T* beta, int *idxz, int *idxb_block, int *idxu_block, int *idxcg_block, int twojmax, 
        int idxb_max, int idxu_max, int idxz_max, int nelements, int bnorm_flag, int inum)
{    
    int N1 = idxu_max*nelements*inum;
    cpuArraySetValue(ylist_r, (T) 0.0, N1);
    cpuArraySetValue(ylist_i, (T) 0.0, N1);
    
    int jdim = twojmax + 1;         
    int N2 = idxz_max*inum;                          
    for (int idx=0; idx < N2; idx++) {
      int ii = idx%inum;              
      int jjz = (idx-ii)/inum;         
      int jjz10 = jjz*10;
      const int j1 = idxz[jjz10+0];
      const int j2 = idxz[jjz10+1];
      const int j = idxz[jjz10+2];
      const int ma1min = idxz[jjz10+3];
      const int ma2max = idxz[jjz10+4];
      const int na = idxz[jjz10+5];
      const int mb1min = idxz[jjz10+6];
      const int mb2max = idxz[jjz10+7];
      const int nb = idxz[jjz10+8];      
      const int jju = idxz[jjz10+9];
      const T *cgblock = cglist + idxcg_block[j + j2*jdim + j1*jdim*jdim];
      
      for(int elem1 = 0; elem1 < nelements; elem1++)
        for (int elem2 = 0; elem2 < nelements; elem2++) {        
          
            T ztmp_r = 0.0;
            T ztmp_i = 0.0;
            
            int ii1 = ii + inum*idxu_max*elem1;
            int ii2 = ii + inum*idxu_max*elem2;
            int jju1 = idxu_block[j1] + (j1 + 1) * mb1min;
            int jju2 = idxu_block[j2] + (j2 + 1) * mb2max;
            int icgb = mb1min * (j2 + 1) + mb2max;
            for (int ib = 0; ib < nb; ib++) {

                int ma1 = ma1min;
                int ma2 = ma2max;
                int icga = ma1min * (j2 + 1) + ma2max;
                for (int ia = 0; ia < na; ia++) {
                    ztmp_r += cgblock[icgb]*cgblock[icga] * (Utotr[ii1+inum*(jju1+ma1)] * Utotr[ii2+inum*(jju2+ma2)] - Utoti[ii1+inum*(jju1+ma1)] * Utoti[ii2+inum*(jju2+ma2)]);
                    ztmp_i += cgblock[icgb]*cgblock[icga] * (Utotr[ii1+inum*(jju1+ma1)] * Utoti[ii2+inum*(jju2+ma2)] + Utoti[ii1+inum*(jju1+ma1)] * Utotr[ii2+inum*(jju2+ma2)]);                    
                  ma1++;
                  ma2--;
                  icga += j2;
                } // end loop over ia

                jju1 += j1 + 1;
                jju2 -= j2 + 1;
                icgb += j2;
            } // end loop over ib


            if (bnorm_flag) {
              ztmp_i /= j+1;
              ztmp_r /= j+1;
            }            
                   
            for(int elem3 = 0; elem3 < nelements; elem3++) {
              int itriple;  
              T betaj;
              if (j >= j1) {
                const int jjb = idxb_block[j + j2*jdim + j1*jdim*jdim]; 
                itriple = ((elem1 * nelements + elem2) * nelements + elem3) * idxb_max*inum + jjb*inum + ii;
                if (j1 == j) {
                  if (j2 == j) betaj = 3*beta[itriple];
                  else betaj = 2*beta[itriple];
                } else betaj = beta[itriple];          
              } else if (j >= j2) {
                const int jjb = idxb_block[j1 + j2*jdim + j*jdim*jdim];
                itriple = ((elem3 * nelements + elem2) * nelements + elem1) * idxb_max*inum + jjb*inum + ii;
                if (j2 == j) betaj = 2*beta[itriple];
                else betaj = beta[itriple];
              } else {
                const int jjb = idxb_block[j1 + j*jdim + j2*jdim*jdim];
                itriple = ((elem2 * nelements + elem3) * nelements + elem1) * idxb_max*inum + jjb*inum + ii;
                betaj = beta[itriple];
              }
              
              if (!bnorm_flag && j1 > j)
                betaj *= (j1 + 1) / (j + 1.0);
                         
              ylist_r[ii + inum*jju + inum*idxu_max*elem3] += betaj * ztmp_r;
              ylist_i[ii + inum*jju + inum*idxu_max*elem3] += betaj * ztmp_i;        
           }
        }         
    }  
}
template void cpuSnapComputeYi(double*, double*, double*, double*, double*, double*, 
        int*, int*, int*, int*, int, int, int, int, int, int, int);
template void cpuSnapComputeYi(float*, float*, float*, float*, float*, float*,
        int*, int*, int*, int*, int, int, int, int, int, int, int);

template <typename T> void cpuSnapComputeYi(T *ylist_r, T *ylist_i, T *Utotr, T *Utoti, T *cglist, T* beta, 
        int *map, int *type, int *idxz, int *idxb_block, int *idxu_block, int *idxcg_block, int twojmax, 
        int idxb_max, int idxu_max, int idxz_max, int nelements, int ncoeffall, int bnorm_flag, int inum)
{    
    
    int N1 = idxu_max*nelements*inum;
    //printf("%i %i %i %i %i %i\n", inum, idxu_max, nelements, N1, idxz_max, idxb_max);
    cpuArraySetValue(ylist_r, (T) 0.0, N1);
    cpuArraySetValue(ylist_i, (T) 0.0, N1);
    //printf("%i %i %i %i %i %i\n", inum, idxu_max, nelements, N1, idxz_max, idxb_max);
    
    int jdim = twojmax + 1;         
    int N2 = idxz_max*inum;                          
    for (int idx=0; idx < N2; idx++) {
      int ii = idx%inum;              
      int jjz = (idx-ii)/inum;         
      int jjz10 = jjz*10;
      const int j1 = idxz[jjz10+0];
      const int j2 = idxz[jjz10+1];
      const int j = idxz[jjz10+2];
      const int ma1min = idxz[jjz10+3];
      const int ma2max = idxz[jjz10+4];
      const int na = idxz[jjz10+5];
      const int mb1min = idxz[jjz10+6];
      const int mb2max = idxz[jjz10+7];
      const int nb = idxz[jjz10+8];          
      const int jju = idxz[jjz10+9];
      const T *cgblock = cglist + idxcg_block[j + j2*jdim + j1*jdim*jdim];
      
      int itype = type[ii]; // element type of atom i
      int ielem = map[itype];  // index of that element type                        

      for(int elem1 = 0; elem1 < nelements; elem1++)
        for (int elem2 = 0; elem2 < nelements; elem2++) {        
          
            T ztmp_r = 0.0;
            T ztmp_i = 0.0;
            
            int ii1 = ii + inum*idxu_max*elem1;
            int ii2 = ii + inum*idxu_max*elem2;
            int jju1 = idxu_block[j1] + (j1 + 1) * mb1min;
            int jju2 = idxu_block[j2] + (j2 + 1) * mb2max;
            int icgb = mb1min * (j2 + 1) + mb2max;
            for (int ib = 0; ib < nb; ib++) {
                int ma1 = ma1min;
                int ma2 = ma2max;
                int icga = ma1min * (j2 + 1) + ma2max;
                for (int ia = 0; ia < na; ia++) {
                    ztmp_r += cgblock[icgb]*cgblock[icga] * (Utotr[ii1+inum*(jju1+ma1)] * Utotr[ii2+inum*(jju2+ma2)] - Utoti[ii1+inum*(jju1+ma1)] * Utoti[ii2+inum*(jju2+ma2)]);
                    ztmp_i += cgblock[icgb]*cgblock[icga] * (Utotr[ii1+inum*(jju1+ma1)] * Utoti[ii2+inum*(jju2+ma2)] + Utoti[ii1+inum*(jju1+ma1)] * Utotr[ii2+inum*(jju2+ma2)]);                    
                  ma1++;
                  ma2--;
                  icga += j2;
                } // end loop over ia

                jju1 += j1 + 1;
                jju2 -= j2 + 1;
                icgb += j2;
            } // end loop over ib

            if (bnorm_flag) {
              ztmp_i /= j+1;
              ztmp_r /= j+1;
            }            
                               
            for(int elem3 = 0; elem3 < nelements; elem3++) {
              int itriple;  
              T betaj;
              if (j >= j1) {
                const int jjb = idxb_block[j + j2*jdim + j1*jdim*jdim]; 
                itriple = ((elem1 * nelements + elem2) * nelements + elem3) * idxb_max + jjb + 1 + ielem*ncoeffall;
                if (j1 == j) {
                  if (j2 == j) betaj = 3*beta[itriple];
                  else betaj = 2*beta[itriple];
                } else betaj = beta[itriple];          
              } else if (j >= j2) {
                const int jjb = idxb_block[j1 + j2*jdim + j*jdim*jdim];
                itriple = ((elem3 * nelements + elem2) * nelements + elem1) * idxb_max + jjb + 1 + ielem*ncoeffall;
                if (j2 == j) betaj = 2*beta[itriple];
                else betaj = beta[itriple];
              } else {
                const int jjb = idxb_block[j1 + j*jdim + j2*jdim*jdim];
                itriple = ((elem2 * nelements + elem3) * nelements + elem1) * idxb_max + jjb + 1 + ielem*ncoeffall;
                betaj = beta[itriple];
              }
              
              if (!bnorm_flag && j1 > j)
                betaj *= (j1 + 1) / (j + 1.0);
                         
              //printf("%i %i %i %i %i %i\n", inum, ii, jju, elem1, elem2, elem3);
              ylist_r[ii + inum*jju + inum*idxu_max*elem3] += betaj * ztmp_r;
              ylist_i[ii + inum*jju + inum*idxu_max*elem3] += betaj * ztmp_i;        
           }
        }         
    }  
}
template void cpuSnapComputeYi(double*, double*, double*, double*, double*, double*, 
        int*, int*, int*, int*, int*, int*, int, int, int, int, int, int, int, int);
template void cpuSnapComputeYi(float*, float*, float*, float*, float*, float*,
        int*, int*, int*, int*, int*, int*, int, int, int, int, int, int, int, int);

template <typename T> void cpuSnapComputeFi(T *fatom, T *vatom, T *ylist_r, T *ylist_i, T *rootpqarray, T *rij, 
        T *wjelem, T *radelem,  T rmin0, T rfac0, T rcutfac, int *idxu_block, int *map, int *aii, int *ai, int *aj, 
        int *ti, int *tj, int twojmax, int idxu_max, int inum, int anum, int ijnum, int switchflag, int chemflag) 
{                 
  for(int ij=0; ij<ijnum; ij++) {        
    T x = rij[ij*3+0];
    T y = rij[ij*3+1];
    T z = rij[ij*3+2];    
    T rsq = x * x + y * y + z * z;
    T r = sqrt(rsq);
    T rinv = 1.0 / r;
    T ux = x * rinv;
    T uy = y * rinv;
    T uz = z * rinv;

    T rcutij = (radelem[ti[ij]]+radelem[tj[ij]])*rcutfac; //(radelem[type[ii]]+radelem[type[jj]])*rcutfac;
    T z0 = r / tan((r - rmin0) * rfac0 * M_PI / (rcutij - rmin0));                
    T dz0dr = z0 / r - (r*rfac0 * M_PI / (rcutij - rmin0)) * (rsq + z0 * z0) / rsq;

    T sfac = 0.0, dsfac = 0.0;        
    if (switchflag == 0) {
        sfac = 1.0;
        dsfac = 0.0;
    }
    else if (switchflag == 1) {
        if (r <= rmin0) {
            sfac = 1.0;
            dsfac = 0.0;
        }
        else if(r > rcutij) {
            sfac = 1.0;
            dsfac = 0.0;
        }
        else {
            T rcutfac = M_PI / (rcutij - rmin0);
            sfac =  0.5 * (cos((r - rmin0) * rcutfac) + 1.0);   
            dsfac = -0.5 * sin((r - rmin0) * rcutfac) * rcutfac;
        }
    } 
    sfac *= wjelem[tj[ij]];
    dsfac *= wjelem[tj[ij]];
    
    T a_r, a_i, b_r, b_i, rootpq;        
    rcutij = 1.0 / sqrt(r * r + z0 * z0);
    a_r = rcutij * z0;
    a_i = -rcutij * z;
    b_r = rcutij * y;
    b_i = -rcutij * x;

    T u_r, u_i, ux_r, ux_i, uy_r, uy_i, uz_r, uz_i;
    T w_r, w_i, wx_r, wx_i, wy_r, wy_i, wz_r, wz_i;
    u_r = -pow(rcutij, 3.0) * (r + z0 * dz0dr);
    wx_r = u_r * ux;
    wy_r = u_r * uy;
    wz_r = u_r * uz;
    ux_r = dz0dr * ux;
    uy_r = dz0dr * uy;
    uz_r = dz0dr * uz;

    T dardx, daidx, dardy, daidy, dardz, daidz;
    dardx = ux_r * rcutij + z0 * wx_r;
    daidx = -z * wx_r;
    dardy = uy_r * rcutij + z0 * wy_r;
    daidy = -z * wy_r;
    dardz = uz_r * rcutij + z0 * wz_r;
    daidz = -z * wz_r;    
    daidz += -rcutij;

    T dbrdx, dbidx, dbrdy, dbidy, dbrdz, dbidz;
    dbrdx = y * wx_r;
    dbidx = -x * wx_r;    
    dbrdy = y * wy_r;
    dbidy = -x * wy_r;    
    dbrdz = y * wz_r;
    dbidz = -x * wz_r;        
    dbidx += -rcutij;
    dbrdy += rcutij;
    
    // 2Jmax = 10    
    T Pr[11], Pi[11], Qr[9], Qi[9];
    T Prx[11], Pix[11], Qrx[9], Qix[9];
    T Pry[11], Piy[11], Qry[9], Qiy[9];        
    T Prz[11], Piz[11], Qrz[9], Qiz[9];
    Pr[0] = 1.0; Pi[0] = 0.0;    
    Prx[0] = 0.0; Pix[0] = 0.0;    
    Pry[0] = 0.0; Piy[0] = 0.0;    
    Prz[0] = 0.0; Piz[0] = 0.0;        
    
    int jdim = twojmax + 1;
    int njelem = (chemflag==0) ? 0 : (inum*idxu_max*map[tj[ij]]);    
    int i = aii[ij] + njelem;                        
                    
    T e, dedx, dedy, dedz;       
    u_r = 0.5*ylist_r[i];  
    //e    =  sfac*u_r;
    dedx = (dsfac * ux) * u_r;
    dedy = (dsfac * uy) * u_r;
    dedz = (dsfac * uz) * u_r; 
            
    int j, k, ma, mb, mapar;    
    mb = 0;
    for (j = 1; j <= twojmax; j++) {        
        // fill in left side of matrix layer from previous layer
        ma = 0;
        u_r = Pr[ma];
        u_i = Pi[ma];
        ux_r = Prx[ma];
        ux_i = Pix[ma];            
        uy_r = Pry[ma];
        uy_i = Piy[ma];            
        uz_r = Prz[ma];
        uz_i = Piz[ma];                    
        rootpq = rootpqarray[(j - ma)*jdim + (j - mb)];
        Pr[ma] = rootpq * (a_r * u_r + a_i * u_i);
        Pi[ma] = rootpq * (a_r * u_i - a_i * u_r);        
        Prx[ma] = rootpq * (dardx * u_r + daidx * u_i + a_r * ux_r + a_i * ux_i);
        Pix[ma] = rootpq * (dardx * u_i - daidx * u_r + a_r * ux_i - a_i * ux_r);
        Pry[ma] = rootpq * (dardy * u_r + daidy * u_i + a_r * uy_r + a_i * uy_i);
        Piy[ma] = rootpq * (dardy * u_i - daidy * u_r + a_r * uy_i - a_i * uy_r);
        Prz[ma] = rootpq * (dardz * u_r + daidz * u_i + a_r * uz_r + a_i * uz_i);
        Piz[ma] = rootpq * (dardz * u_i - daidz * u_r + a_r * uz_i - a_i * uz_r);                    
        for (ma = 1; ma < j; ma++) {
            w_r = Pr[ma];
            w_i = Pi[ma];
            wx_r = Prx[ma];
            wx_i = Pix[ma];            
            wy_r = Pry[ma];
            wy_i = Piy[ma];            
            wz_r = Prz[ma];
            wz_i = Piz[ma];                                        
            rootpq = rootpqarray[(j - ma)*jdim + (j - mb)];
            rcutij = rootpqarray[ma*jdim + (j - mb)];
            Pr[ma] = rootpq * (a_r * w_r + a_i * w_i) -rcutij * (b_r * u_r + b_i * u_i);
            Pi[ma] = rootpq * (a_r * w_i - a_i * w_r) -rcutij * (b_r * u_i - b_i * u_r);
            Prx[ma] = rootpq * (dardx * w_r + daidx * w_i + a_r * wx_r + a_i * wx_i) -rcutij * (dbrdx * u_r + dbidx * u_i + b_r * ux_r + b_i * ux_i);
            Pix[ma] = rootpq * (dardx * w_i - daidx * w_r + a_r * wx_i - a_i * wx_r) -rcutij * (dbrdx * u_i - dbidx * u_r + b_r * ux_i - b_i * ux_r);
            Pry[ma] = rootpq * (dardy * w_r + daidy * w_i + a_r * wy_r + a_i * wy_i) -rcutij * (dbrdy * u_r + dbidy * u_i + b_r * uy_r + b_i * uy_i);
            Piy[ma] = rootpq * (dardy * w_i - daidy * w_r + a_r * wy_i - a_i * wy_r) -rcutij * (dbrdy * u_i - dbidy * u_r + b_r * uy_i - b_i * uy_r);
            Prz[ma] = rootpq * (dardz * w_r + daidz * w_i + a_r * wz_r + a_i * wz_i) -rcutij * (dbrdz * u_r + dbidz * u_i + b_r * uz_r + b_i * uz_i);
            Piz[ma] = rootpq * (dardz * w_i - daidz * w_r + a_r * wz_i - a_i * wz_r) -rcutij * (dbrdz * u_i - dbidz * u_r + b_r * uz_i - b_i * uz_r);            
            u_r = w_r;
            u_i = w_i;
            ux_r = wx_r;
            ux_i = wx_i;
            uy_r = wy_r;
            uy_i = wy_i;
            uz_r = wz_r;
            uz_i = wz_i;            
        }
        ma = j;
        rcutij = rootpqarray[ma*jdim + (j - mb)];
        Pr[ma] = -rcutij * (b_r * u_r + b_i * u_i);
        Pi[ma] = -rcutij * (b_r * u_i - b_i * u_r);                        
        Prx[ma] =-rcutij * (dbrdx * u_r + dbidx * u_i + b_r * ux_r + b_i * ux_i);
        Pix[ma] =-rcutij * (dbrdx * u_i - dbidx * u_r + b_r * ux_i - b_i * ux_r);
        Pry[ma] =-rcutij * (dbrdy * u_r + dbidy * u_i + b_r * uy_r + b_i * uy_i);
        Piy[ma] =-rcutij * (dbrdy * u_i - dbidy * u_r + b_r * uy_i - b_i * uy_r);
        Prz[ma] =-rcutij * (dbrdz * u_r + dbidz * u_i + b_r * uz_r + b_i * uz_i);
        Piz[ma] =-rcutij * (dbrdz * u_i - dbidz * u_r + b_r * uz_i - b_i * uz_r);
                                
        if (j==(2*mb+1)) { // store Qr, Qi, for the next mb level
            mapar = (((j+1)/2)%2 == 0) ? -1 : 1;
            for (ma = 0; ma <= j; ma++) {     
                if (mapar == 1) {                    
                    Qr[j-ma] = Pr[ma];
                    Qi[j-ma] = -Pi[ma];
                    Qrx[j-ma] =  Prx[ma];
                    Qix[j-ma] = -Pix[ma];
                    Qry[j-ma] =  Pry[ma];
                    Qiy[j-ma] = -Piy[ma];
                    Qrz[j-ma] =  Prz[ma];
                    Qiz[j-ma] = -Piz[ma];                    
                } else {
                    Qr[j-ma] = -Pr[ma];
                    Qi[j-ma] =  Pi[ma];
                    Qrx[j-ma] = -Prx[ma];
                    Qix[j-ma] =  Pix[ma];
                    Qry[j-ma] = -Pry[ma];
                    Qiy[j-ma] =  Piy[ma];
                    Qrz[j-ma] = -Prz[ma];
                    Qiz[j-ma] =  Piz[ma];                    
                }
                mapar = -mapar;
            }                              
        }
        
        k =  1 + (j+1)*mb;
        for (ma = 2; ma <= j; ma++)
            k += ma*ma;                    
        for (ma = 0; ma <= j; ma++) {                            
            rsq = ylist_r[i + inum*k]; 
            rinv = ylist_i[i + inum*k];
            //e += sfac*(Pr[ma] * rsq + Pi[ma] * rinv);
            dedx += (dsfac * Pr[ma] * ux + sfac * Prx[ma]) * rsq + (dsfac * Pi[ma] * ux + sfac * Pix[ma]) * rinv;
            dedy += (dsfac * Pr[ma] * uy + sfac * Pry[ma]) * rsq + (dsfac * Pi[ma] * uy + sfac * Piy[ma]) * rinv;
            dedz += (dsfac * Pr[ma] * uz + sfac * Prz[ma]) * rsq + (dsfac * Pi[ma] * uz + sfac * Piz[ma]) * rinv;            
            k += 1;
        }                   
    }
    
    for (mb = 1; 2*mb <= twojmax; mb++) {     
        for (ma = 0; ma < 2*mb; ma++) {                      
            Pr[ma] = Qr[ma];
            Pi[ma] = Qi[ma];
            Prx[ma] = Qrx[ma];
            Pix[ma] = Qix[ma];
            Pry[ma] = Qry[ma];
            Piy[ma] = Qiy[ma];
            Prz[ma] = Qrz[ma];
            Piz[ma] = Qiz[ma];            
        }                
        for (j = 2*mb; j <= twojmax; j++) { 
            ma = 0;
            u_r = Pr[ma];
            u_i = Pi[ma];
            ux_r = Prx[ma];
            ux_i = Pix[ma];            
            uy_r = Pry[ma];
            uy_i = Piy[ma];            
            uz_r = Prz[ma];
            uz_i = Piz[ma];                                
            rootpq = rootpqarray[(j - ma)*jdim + (j - mb)];
            Pr[ma] = rootpq * (a_r * u_r + a_i * u_i);
            Pi[ma] = rootpq * (a_r * u_i - a_i * u_r);            
            Prx[ma] = rootpq * (dardx * u_r + daidx * u_i + a_r * ux_r + a_i * ux_i);
            Pix[ma] = rootpq * (dardx * u_i - daidx * u_r + a_r * ux_i - a_i * ux_r);
            Pry[ma] = rootpq * (dardy * u_r + daidy * u_i + a_r * uy_r + a_i * uy_i);
            Piy[ma] = rootpq * (dardy * u_i - daidy * u_r + a_r * uy_i - a_i * uy_r);
            Prz[ma] = rootpq * (dardz * u_r + daidz * u_i + a_r * uz_r + a_i * uz_i);
            Piz[ma] = rootpq * (dardz * u_i - daidz * u_r + a_r * uz_i - a_i * uz_r);                                
            for (ma = 1; ma < j; ma++) {
                w_r = Pr[ma];
                w_i = Pi[ma];
                wx_r = Prx[ma];
                wx_i = Pix[ma];            
                wy_r = Pry[ma];
                wy_i = Piy[ma];            
                wz_r = Prz[ma];
                wz_i = Piz[ma];                                            
                rootpq = rootpqarray[(j - ma)*jdim + (j - mb)];
                rcutij = rootpqarray[ma*jdim + (j - mb)];
                Pr[ma] = rootpq * (a_r * w_r + a_i * w_i) -rcutij * (b_r * u_r + b_i * u_i);
                Pi[ma] = rootpq * (a_r * w_i - a_i * w_r) -rcutij * (b_r * u_i - b_i * u_r);
                Prx[ma] = rootpq * (dardx * w_r + daidx * w_i + a_r * wx_r + a_i * wx_i) -rcutij * (dbrdx * u_r + dbidx * u_i + b_r * ux_r + b_i * ux_i);
                Pix[ma] = rootpq * (dardx * w_i - daidx * w_r + a_r * wx_i - a_i * wx_r) -rcutij * (dbrdx * u_i - dbidx * u_r + b_r * ux_i - b_i * ux_r);
                Pry[ma] = rootpq * (dardy * w_r + daidy * w_i + a_r * wy_r + a_i * wy_i) -rcutij * (dbrdy * u_r + dbidy * u_i + b_r * uy_r + b_i * uy_i);
                Piy[ma] = rootpq * (dardy * w_i - daidy * w_r + a_r * wy_i - a_i * wy_r) -rcutij * (dbrdy * u_i - dbidy * u_r + b_r * uy_i - b_i * uy_r);
                Prz[ma] = rootpq * (dardz * w_r + daidz * w_i + a_r * wz_r + a_i * wz_i) -rcutij * (dbrdz * u_r + dbidz * u_i + b_r * uz_r + b_i * uz_i);
                Piz[ma] = rootpq * (dardz * w_i - daidz * w_r + a_r * wz_i - a_i * wz_r) -rcutij * (dbrdz * u_i - dbidz * u_r + b_r * uz_i - b_i * uz_r);                            
                u_r = w_r;
                u_i = w_i;
                ux_r = wx_r;
                ux_i = wx_i;
                uy_r = wy_r;
                uy_i = wy_i;
                uz_r = wz_r;
                uz_i = wz_i;                
            }
            ma = j;
            rcutij = rootpqarray[ma*jdim + (j - mb)];
            Pr[ma] = -rcutij * (b_r * u_r + b_i * u_i);
            Pi[ma] = -rcutij * (b_r * u_i - b_i * u_r);       
            Prx[ma] =-rcutij * (dbrdx * u_r + dbidx * u_i + b_r * ux_r + b_i * ux_i);
            Pix[ma] =-rcutij * (dbrdx * u_i - dbidx * u_r + b_r * ux_i - b_i * ux_r);
            Pry[ma] =-rcutij * (dbrdy * u_r + dbidy * u_i + b_r * uy_r + b_i * uy_i);
            Piy[ma] =-rcutij * (dbrdy * u_i - dbidy * u_r + b_r * uy_i - b_i * uy_r);
            Prz[ma] =-rcutij * (dbrdz * u_r + dbidz * u_i + b_r * uz_r + b_i * uz_i);
            Piz[ma] =-rcutij * (dbrdz * u_i - dbidz * u_r + b_r * uz_i - b_i * uz_r);
            
            if (j==(2*mb)) {
                mapar = 1;
                for (ma = 0; ma <= j/2; ma++) {
                    if (mapar == 1) {                    
                        Pr[j/2+ma] = Pr[j/2-ma];
                        Pi[j/2+ma] = -Pi[j/2-ma];
                        Prx[j/2+ma] = Prx[j/2-ma];
                        Pix[j/2+ma] = -Pix[j/2-ma];
                        Pry[j/2+ma] = Pry[j/2-ma];
                        Piy[j/2+ma] = -Piy[j/2-ma];
                        Prz[j/2+ma] = Prz[j/2-ma];
                        Piz[j/2+ma] = -Piz[j/2-ma];                        
                    } else {
                        Pr[j/2+ma] = -Pr[j/2-ma];
                        Pi[j/2+ma] = Pi[j/2-ma];
                        Prx[j/2+ma] = -Prx[j/2-ma];
                        Pix[j/2+ma] =  Pix[j/2-ma];
                        Pry[j/2+ma] = -Pry[j/2-ma];
                        Piy[j/2+ma] =  Piy[j/2-ma];
                        Prz[j/2+ma] = -Prz[j/2-ma];
                        Piz[j/2+ma] =  Piz[j/2-ma];                        
                    }
                    mapar = -mapar;        
                }                                                        
            }
            
            if (j==(2*mb)) {
                mapar = 1;
                for (ma = 0; ma <= j; ma++) {
                    if (mapar == 1) {                    
                        Prx[j/2+ma] = Prx[j/2-ma];
                        Pix[j/2+ma] = -Pix[j/2-ma];
                        Pry[j/2+ma] = Pry[j/2-ma];
                        Piy[j/2+ma] = -Piy[j/2-ma];
                        Prz[j/2+ma] = Prz[j/2-ma];
                        Piz[j/2+ma] = -Piz[j/2-ma];                        
                    } else {
                        Prx[j/2+ma] = -Prx[j/2-ma];
                        Pix[j/2+ma] =  Pix[j/2-ma];
                        Pry[j/2+ma] = -Pry[j/2-ma];
                        Piy[j/2+ma] =  Piy[j/2-ma];
                        Prz[j/2+ma] = -Prz[j/2-ma];
                        Piz[j/2+ma] =  Piz[j/2-ma];                        
                    }
                    mapar = -mapar;        
                }                                                        
            }
                        
            // store Qr, Qi, for the next mb level
            if (j==(2*mb+1)) {
                mapar = (((j+1)/2)%2 == 0) ? -1 : 1;
                for (ma = 0; ma <= j; ma++) {     
                    if (mapar == 1) {                    
                        Qr[j-ma] = Pr[ma];
                        Qi[j-ma] = -Pi[ma];  
                        Qrx[j-ma] =  Prx[ma];
                        Qix[j-ma] = -Pix[ma];
                        Qry[j-ma] =  Pry[ma];
                        Qiy[j-ma] = -Piy[ma];
                        Qrz[j-ma] =  Prz[ma];
                        Qiz[j-ma] = -Piz[ma];                                            
                    } else {
                        Qr[j-ma] = -Pr[ma];
                        Qi[j-ma] =  Pi[ma];
                        Qrx[j-ma] = -Prx[ma];
                        Qix[j-ma] =  Pix[ma];
                        Qry[j-ma] = -Pry[ma];
                        Qiy[j-ma] =  Piy[ma];
                        Qrz[j-ma] = -Prz[ma];
                        Qiz[j-ma] =  Piz[ma];                                            
                    }
                    mapar = -mapar;
                }                                                
            }
            
            k =  1 + (j+1)*mb;
            for (ma = 2; ma <= j; ma++)
                k += ma*ma;                            
            if (j==(2*mb)) {
                for (ma = 0; ma < mb; ma++) {
                    rsq = ylist_r[i + inum*k]; 
                    rinv = ylist_i[i + inum*k];         
                    //e += sfac*(Pr[ma] * rsq + Pi[ma] * rinv);
                    dedx += (dsfac * Pr[ma] * ux + sfac * Prx[ma]) * rsq + (dsfac * Pi[ma] * ux + sfac * Pix[ma]) * rinv;
                    dedy += (dsfac * Pr[ma] * uy + sfac * Pry[ma]) * rsq + (dsfac * Pi[ma] * uy + sfac * Piy[ma]) * rinv;
                    dedz += (dsfac * Pr[ma] * uz + sfac * Prz[ma]) * rsq + (dsfac * Pi[ma] * uz + sfac * Piz[ma]) * rinv;
                    k += 1; 
                }     
                ma = mb;
                rsq = 0.5*ylist_r[i + inum*k]; 
                rinv = 0.5*ylist_i[i + inum*k];        
                //e += sfac*(Pr[ma] * rsq + Pi[ma] * rinv);
                dedx += (dsfac * Pr[ma] * ux + sfac * Prx[ma]) * rsq + (dsfac * Pi[ma] * ux + sfac * Pix[ma]) * rinv;
                dedy += (dsfac * Pr[ma] * uy + sfac * Pry[ma]) * rsq + (dsfac * Pi[ma] * uy + sfac * Piy[ma]) * rinv;
                dedz += (dsfac * Pr[ma] * uz + sfac * Prz[ma]) * rsq + (dsfac * Pi[ma] * uz + sfac * Piz[ma]) * rinv;
            }
            else {
                for (ma = 0; ma <= j; ma++) {
                    rsq = ylist_r[i + inum*k]; 
                    rinv = ylist_i[i + inum*k];             
                    //e += sfac*(Pr[ma] * rsq + Pi[ma] * rinv);
                    dedx += (dsfac * Pr[ma] * ux + sfac * Prx[ma]) * rsq + (dsfac * Pi[ma] * ux + sfac * Pix[ma]) * rinv;
                    dedy += (dsfac * Pr[ma] * uy + sfac * Pry[ma]) * rsq + (dsfac * Pi[ma] * uy + sfac * Piy[ma]) * rinv;
                    dedz += (dsfac * Pr[ma] * uz + sfac * Prz[ma]) * rsq + (dsfac * Pi[ma] * uz + sfac * Piz[ma]) * rinv;
                    k += 1; 
                }                
            }            
        }
    }      

    ma = ai[ij];        
    mb = aj[ij];                    
    //eatom[ma] += 2.0*e;
    
    dedx = 2.0*dedx;      
    dedy = 2.0*dedy;      
    dedz = 2.0*dedz;             
    fatom[0+3*ma] += dedx;
    fatom[1+3*ma] += dedy;
    fatom[2+3*ma] += dedz;
    fatom[0+3*mb] -= dedx;
    fatom[1+3*mb] -= dedy;
    fatom[2+3*mb] -= dedz;    
    
    rootpq = -0.5;
    dbrdx = rootpq*x*dedx;
    dbidx = rootpq*y*dedy;
    dbrdy = rootpq*z*dedz;
    dbidy = rootpq*x*dedy;
    dbrdz = rootpq*x*dedz;
    dbidz = rootpq*y*dedz;        
    vatom[0*anum+ma] += dbrdx;
    vatom[1*anum+ma] += dbidx;
    vatom[2*anum+ma] += dbrdy;
    vatom[3*anum+ma] += dbidy;
    vatom[4*anum+ma] += dbrdz;
    vatom[5*anum+ma] += dbidz;        
    vatom[0*anum+mb] += dbrdx;
    vatom[1*anum+mb] += dbidx;
    vatom[2*anum+mb] += dbrdy;
    vatom[3*anum+mb] += dbidy;
    vatom[4*anum+mb] += dbrdz;
    vatom[5*anum+mb] += dbidz;            
  }   
}
template void cpuSnapComputeFi(double*, double*, double*, double*, double*, double*,  double*, double*, 
         double, double, double, int*, int*, int*, int*, int*, int*, int*, int, int, int, int, int, int, int);
template void cpuSnapComputeFi(float*, float*, float*, float*, float*, float*, float*, float*, 
        float, float, float, int*, int*, int*, int*, int*, int*, int*, int, int, int, int, int, int, int);

#endif
        