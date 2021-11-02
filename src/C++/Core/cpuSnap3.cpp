/***************************************************************************                               
                    Molecular Dynamics Potentials (MDP)
                           CESMIX-MIT Project  
 
 Contributing authors: Ngoc-Cuong Nguyen (cuongng@mit.edu, exapde@gmail.com)
 ***************************************************************************/

#ifndef MDP_CPUSNAP3
#define MDP_CPUSNAP3

template <typename T> void cpuComputeUi(T *Stotr, T *Stoti, T *rootpqarray, T *rij, T *wjelem, T *radelem, 
        T rmin0, T rfac0, T rcutfac, int *idxu_block, int *map, int *ai, int *ti, int *tj, 
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
    int i = ai[ij] + njelem;                        
    Stotr[i] += sfac; // atomic add   
    
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
            Stotr[in] += sfac*Pr[ma]; // atomic add   
            Stoti[in] += sfac*Pi[ma]; // atomic add                       
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
                Stotr[in] += sfac*Pr[ma]; // atomic add   
                Stoti[in] += sfac*Pi[ma]; // atomic add                       
                k += 1; 
            }                                                           
        }
    }        
  }
};
template void cpuComputeUi(double*, double*, double*, double*, double*, double*, double, double, double, 
        int*, int*, int*, int*, int*, int, int, int, int, int, int);
template void cpuComputeUi(float*, float*, float*, float*, float*, float*, float, float, float, 
        int*, int*, int*, int*, int*, int, int, int, int, int, int);

template <typename T> void cpuAddWself2Ui(T *Stotr, T *Stoti, T wself, int *idxu_block, 
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
        int ielem = (chemflag) ? map[type[ai[ii]]]: 0;                
        int nmax = ii + inum*idxu_max*jelem;
        
        int jju = idxu_block[j];                
        for (int mb = 0; mb <= j; mb++) {
            for (int ma = 0; ma <= j; ma++) {                
                if (jelem == ielem || wselfall_flag)
                    if (ma==mb)                        
                        Stotr[inum*jju + nmax] += wself;                                     
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
                    Stotr[njjup] = Stotr[njju];
                    Stoti[njjup] = -Stoti[njju];
                } else {
                    Stotr[njjup] = -Stotr[njju];
                    Stoti[njjup] =  Stoti[njju];
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

template <typename T> void cpuSnapComputeEi(T *eatom, T *Stotr, T *Stoti, T *cglist, 
        T *bzero, T *coeffelem, int *map, int *type, int *idxb, int *idxcg_block, int *idxu_block, 
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
                        ztmp_r += cgblock[icgb]*cgblock[icga] * (Stotr[ii1+inum*(jju1+ma1)] * Stotr[ii2+inum*(jju2+ma2)] - Stoti[ii1+inum*(jju1+ma1)] * Stoti[ii2+inum*(jju2+ma2)]);
                        ztmp_i += cgblock[icgb]*cgblock[icga] * (Stotr[ii1+inum*(jju1+ma1)] * Stoti[ii2+inum*(jju2+ma2)] + Stoti[ii1+inum*(jju1+ma1)] * Stotr[ii2+inum*(jju2+ma2)]);                    
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
                sumzu += Stotr[idu] * ztmp_r + Stoti[idu] * ztmp_i;
                //sumzu += Stotr[idu] * zlist_r[idz] + Stoti[idu] * zlist_i[idz];                                              
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
                        ztmp_r += cgblock[icgb]*cgblock[icga] * (Stotr[ii1+inum*(jju1+ma1)] * Stotr[ii2+inum*(jju2+ma2)] - Stoti[ii1+inum*(jju1+ma1)] * Stoti[ii2+inum*(jju2+ma2)]);
                        ztmp_i += cgblock[icgb]*cgblock[icga] * (Stotr[ii1+inum*(jju1+ma1)] * Stoti[ii2+inum*(jju2+ma2)] + Stoti[ii1+inum*(jju1+ma1)] * Stotr[ii2+inum*(jju2+ma2)]);                    
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
                sumzu += Stotr[idu] * ztmp_r + Stoti[idu] * ztmp_i;                
                //sumzu += Stotr[idu] * zlist_r[idz] + Stoti[idu] * zlist_i[idz];                                
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
                    ztmp_r += cgblock[icgb]*cgblock[icga] * (Stotr[ii1+inum*(jju1+ma1)] * Stotr[ii2+inum*(jju2+ma2)] - Stoti[ii1+inum*(jju1+ma1)] * Stoti[ii2+inum*(jju2+ma2)]);
                    ztmp_i += cgblock[icgb]*cgblock[icga] * (Stotr[ii1+inum*(jju1+ma1)] * Stoti[ii2+inum*(jju2+ma2)] + Stoti[ii1+inum*(jju1+ma1)] * Stotr[ii2+inum*(jju2+ma2)]);                    
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
            sumzu += 0.5*(Stotr[idu] * ztmp_r + Stoti[idu] * ztmp_i);            
            //sumzu += 0.5 * (Stotr[idu] * zlist_r[idz] + Stoti[idu] * zlist_i[idz]);            
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
            eatom[ii] += coeffelem[ielem*ncoeffall];
        eatom[ii] += coeffelem[itriple]*sumzu;                           
    }
}
template void cpuSnapComputeEi(double*, double*, double*, double*, double*, double*,
        int*, int*, int*, int*, int*, int, int, int, int, int, int, int, int, int, int);
template void cpuSnapComputeEi(float*, float*, float*, float*, float*, float*, 
        int*, int*, int*, int*, int*, int, int, int, int, int, int, int, int, int, int);

template <typename T> void cpuComputeBi(T *blist, T *zlist_r, T *zlist_i, T *Stotr, T *Stoti, T *cglist, T *bzero, 
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
                        ztmp_r += cgblock[icgb]*cgblock[icga] * (Stotr[ii1+inum*(jju1+ma1)] * Stotr[ii2+inum*(jju2+ma2)] - Stoti[ii1+inum*(jju1+ma1)] * Stoti[ii2+inum*(jju2+ma2)]);
                        ztmp_i += cgblock[icgb]*cgblock[icga] * (Stotr[ii1+inum*(jju1+ma1)] * Stoti[ii2+inum*(jju2+ma2)] + Stoti[ii1+inum*(jju1+ma1)] * Stotr[ii2+inum*(jju2+ma2)]);                    
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
                sumzu += Stotr[idu] * ztmp_r + Stoti[idu] * ztmp_i;
                //sumzu += Stotr[idu] * zlist_r[idz] + Stoti[idu] * zlist_i[idz];
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
                        ztmp_r += cgblock[icgb]*cgblock[icga] * (Stotr[ii1+inum*(jju1+ma1)] * Stotr[ii2+inum*(jju2+ma2)] - Stoti[ii1+inum*(jju1+ma1)] * Stoti[ii2+inum*(jju2+ma2)]);
                        ztmp_i += cgblock[icgb]*cgblock[icga] * (Stotr[ii1+inum*(jju1+ma1)] * Stoti[ii2+inum*(jju2+ma2)] + Stoti[ii1+inum*(jju1+ma1)] * Stotr[ii2+inum*(jju2+ma2)]);                    
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
                sumzu += Stotr[idu] * ztmp_r + Stoti[idu] * ztmp_i;                
                //sumzu += Stotr[idu] * zlist_r[idz] + Stoti[idu] * zlist_i[idz];
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
                    ztmp_r += cgblock[icgb]*cgblock[icga] * (Stotr[ii1+inum*(jju1+ma1)] * Stotr[ii2+inum*(jju2+ma2)] - Stoti[ii1+inum*(jju1+ma1)] * Stoti[ii2+inum*(jju2+ma2)]);
                    ztmp_i += cgblock[icgb]*cgblock[icga] * (Stotr[ii1+inum*(jju1+ma1)] * Stoti[ii2+inum*(jju2+ma2)] + Stoti[ii1+inum*(jju1+ma1)] * Stotr[ii2+inum*(jju2+ma2)]);                    
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
            sumzu += 0.5*(Stotr[idu] * ztmp_r + Stoti[idu] * ztmp_i);            
            //sumzu += 0.5 * (Stotr[idu] * zlist_r[idz] + Stoti[idu] * zlist_i[idz]);
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
template void cpuComputeBi(double*, double*, double*, double*, double*, double*, double*, 
        int*, int*, int*, int*, int, int, int, int, int, int, int, int, int);
template void cpuComputeBi(float*, float*, float*, float*, float*, float*, float*,
        int*, int*, int*, int*, int, int, int, int, int, int, int, int, int);

template <typename T> void cpuComputeBi(T *blist, T *zlist_r, T *zlist_i, T *Stotr, T *Stoti, T *bzero, 
        int *idxb, int *idxu_block, int *idxz_block, int twojmax, int idxb_max, int idxu_max, 
        int idxz_max, int nelements, int bzero_flag, int wselfall_flag, int inum)
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

        int jjz = idxz_block[j + j2*jdim + j1*jdim*jdim];
        int jju = idxu_block[j];
        int idu;
        int idz;
        T sumzu = 0.0;
        for (int mb = 0; 2 * mb < j; mb++)
            for (int ma = 0; ma <= j; ma++) {
                idu = ii + inum*jju + nu_max*elem3;        
                idz = ii + inum*jjz + nz_max*idouble;        
                sumzu += Stotr[idu] * zlist_r[idz] + Stoti[idu] * zlist_i[idz];
                jjz++;
                jju++;
            } // end loop over ma, mb

        // For j even, handle middle column
        if (j % 2 == 0) {
            int mb = j / 2;
            for (int ma = 0; ma < mb; ma++) {
                idu = ii + inum*jju + nu_max*elem3;        
                idz = ii + inum*jjz + nz_max*idouble;        
                sumzu += Stotr[idu] * zlist_r[idz] + Stoti[idu] * zlist_i[idz];
                jjz++;
                jju++;
            }
            idu = ii + inum*jju + nu_max*elem3;        
            idz = ii + inum*jjz + nz_max*idouble;        
            sumzu += 0.5 * (Stotr[idu] * zlist_r[idz] + Stoti[idu] * zlist_i[idz]);
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
template void cpuComputeBi(double*, double*, double*, double*, double*, double*, 
        int*, int*, int*, int, int, int, int, int, int, int, int);
template void cpuComputeBi(float*, float*, float*, float*, float*, float*,
        int*, int*, int*, int, int, int, int, int, int, int, int);

template <typename T> void cpuComputeYi(T *ylist_r, T *ylist_i, T *Stotr, T *Stoti, 
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
                    ztmp_r += cgblock[icgb]*cgblock[icga] * (Stotr[ii1+inum*(jju1+ma1)] * Stotr[ii2+inum*(jju2+ma2)] - Stoti[ii1+inum*(jju1+ma1)] * Stoti[ii2+inum*(jju2+ma2)]);
                    ztmp_i += cgblock[icgb]*cgblock[icga] * (Stotr[ii1+inum*(jju1+ma1)] * Stoti[ii2+inum*(jju2+ma2)] + Stoti[ii1+inum*(jju1+ma1)] * Stotr[ii2+inum*(jju2+ma2)]);                    
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
            // store Z_{j1,j2,j}^{elem1,elem2}(ma,mb)                
            //zlist_r[idx + N2*elem2 + N2*nelements*elem1] = ztmp_r;
            //zlist_i[idx + N2*elem2 + N2*nelements*elem1] = ztmp_i;          
       
            int jju = idxz[jjz10+9];
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
template void cpuComputeYi(double*, double*, double*, double*, double*, double*, 
        int*, int*, int*, int*, int, int, int, int, int, int, int);
template void cpuComputeYi(float*, float*, float*, float*, float*, float*,
        int*, int*, int*, int*, int, int, int, int, int, int, int);

template <typename T> void cpuSnapComputeFi(T *fatom, T *ylist_r, T *ylist_i, T *rootpqarray, T *rij, 
        T *wjelem, T *radelem,  T rmin0, T rfac0, T rcutfac, int *idxu_block, int *map, int *ai, int *aj, 
        int *ti, int *tj, int twojmax, int idxu_max, int inum, int ijnum, int switchflag, int chemflag) 
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
    int i = ai[ij] + njelem;                        
                    
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
  }   
}
template void cpuSnapComputeFi(double*, double*, double*, double*, double*,  double*, double*, 
         double, double, double, int*, int*, int*, int*, int*, int*, int, int, int, int, int, int);
template void cpuSnapComputeFi(float*, float*, float*, float*, float*, float*, float*, 
        float, float, float, int*, int*, int*, int*, int*, int*, int, int, int, int, int, int);

template <typename T> void cpuComputeEi(T *ei, T *ylist_r, T *ylist_i, 
         T *Stot_r, T *Stot_i, int twojmax, int idxu_max, int nelem, int inum)
{       
    //printf("%i %i %i %i\n", twojmax,  idxu_max,  nelem,  inum);                    
    int inumax = inum*idxu_max;
    for (int i=0; i < inum; i++) {
        
        T e = 0.0;
        for (int j = 0; j <= twojmax; j++) {
            int m, k =  0;
            for (int ma = 1; ma <= j; ma++)
                k += ma*ma;                                                                
            
            for (int mb = 0; 2 * mb < j; mb++)
                for (int ma = 0; ma <= j; ma++) {
                    for (int n = 0; n< nelem; n++) {
                        m = i + inum*k + inumax*n;        
                        e += (ylist_r[m]*Stot_r[m] + ylist_i[m]*Stot_i[m]);                          
                    }
                    k += 1;
                }
        
            if (j % 2 == 0) {
                int mb = j / 2;
                for (int ma = 0; ma < mb; ma++) {
                    for (int n = 0; n< nelem; n++) {
                        m = i + inum*k + inumax*n;        
                        e += (ylist_r[m]*Stot_r[m] + ylist_i[m]*Stot_i[m]);                          
                    }                    
                    k += 1;
                }
                for (int n = 0; n< nelem; n++) {
                    m = i + inum*k + inumax*n;        
                    e += 0.5*(ylist_r[m]*Stot_r[m] + ylist_i[m]*Stot_i[m]);                      
                }                
            } 
        }
        ei[i] += 2.0*e;
        
//         T e = 0.0;
//         int j, k, m, n;
//         k = 0;
//         for (n = 0; n< nelem; n++) {
//             m = i + inum*k + inumax*n;
//             e += 0.5*(ylist_r[m]*Stot_r[m] + ylist_i[m]*Stot_i[m]);  
//         }
//         
//         int ma, mb = 0;
//         for (j = 1; j <= twojmax; j++) {        
//             k =  1 + (j+1)*mb;
//             for (ma = 2; ma <= j; ma++)
//                 k += ma*ma;                    
//             for (ma = 0; ma <= j; ma++) {                                            
//                 for (n = 0; n< nelem; n++) {
//                     m = i + inum*k + inumax*n;
//                     e += (ylist_r[m]*Stot_r[m] + ylist_i[m]*Stot_i[m]);  
//                 }                
//                 k += 1;
//             }                               
//         }        
//         
//         for (mb = 1; 2*mb <= twojmax; mb++) {     
//             for (j = 2*mb; j <= twojmax; j++) {            
//                 k =  1 + (j+1)*mb;
//                 for (ma = 2; ma <= j; ma++)
//                     k += ma*ma;                            
//                 if (j==(2*mb)) {
//                     for (ma = 0; ma < mb; ma++) {                        
//                         for (n = 0; n< nelem; n++) {
//                             m = i + inum*k + inumax*n;
//                             e += (ylist_r[m]*Stot_r[m] + ylist_i[m]*Stot_i[m]);  
//                         }                                        
//                         k += 1; 
//                     }     
//                     ma = mb;
//                     for (n = 0; n< nelem; n++) {
//                         m = i + inum*k + inumax*n;
//                         e += 0.5*(ylist_r[m]*Stot_r[m] + ylist_i[m]*Stot_i[m]);  
//                     }                                    
//                 }
//                 else {
//                     for (ma = 0; ma <= j; ma++) {
//                         for (n = 0; n< nelem; n++) {
//                             m = i + inum*k + inumax*n;
//                             e += (ylist_r[m]*Stot_r[m] + ylist_i[m]*Stot_i[m]);  
//                         }                                        
//                         k += 1; 
//                     }                
//                 }                    
//             }
//         }
//         ei[i] += 2.0*e;
     }
}
template void cpuComputeEi(double*, double*, double*, double*, double*, int, int, int, int);
template void cpuComputeEi(float*, float*, float*, float*, float*, int, int, int, int);

template <typename T> void cpuComputeDeidrj(T *fatom, T *ylist_r, T *ylist_i, T *rootpqarray, T *rij, 
        T *wjelem, T *radelem,  T rmin0, T rfac0, T rcutfac, int *idxu_block, int *map, int *ai, int *aj, 
        int *ti, int *tj, int twojmax, int idxu_max, int inum, int ijnum, int switchflag, int chemflag) 
{               
  //cpuArraySetValue(dedr, (T) 0.0, ijnum*3);
  
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
    int i = ai[ij] + njelem;                        
                    
    T dedx, dedy, dedz;       
    dedx = 0.0; dedy = 0.0; dedz = 0.0;                
    u_r = 0.5*ylist_r[i];    
    dedx += (dsfac * ux) * u_r;
    dedy += (dsfac * uy) * u_r;
    dedz += (dsfac * uz) * u_r; 
            
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
            //T y_r = ylist_r[in]; 
            //T y_i = ylist_i[in]; 
            rsq = ylist_r[i + inum*k]; 
            rinv = ylist_i[i + inum*k];
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
                    //T y_r = ylist_r[in];
                    //T y_i = ylist_i[in];
                    rsq = ylist_r[i + inum*k]; 
                    rinv = ylist_i[i + inum*k];                    
                    dedx += (dsfac * Pr[ma] * ux + sfac * Prx[ma]) * rsq + (dsfac * Pi[ma] * ux + sfac * Pix[ma]) * rinv;
                    dedy += (dsfac * Pr[ma] * uy + sfac * Pry[ma]) * rsq + (dsfac * Pi[ma] * uy + sfac * Piy[ma]) * rinv;
                    dedz += (dsfac * Pr[ma] * uz + sfac * Prz[ma]) * rsq + (dsfac * Pi[ma] * uz + sfac * Piz[ma]) * rinv;
                    k += 1; 
                }     
                ma = mb;
                //T y_r = 0.5*ylist_r[in];
                //T y_i = 0.5*ylist_i[in];
                rsq = 0.5*ylist_r[i + inum*k]; 
                rinv = 0.5*ylist_i[i + inum*k];                
                dedx += (dsfac * Pr[ma] * ux + sfac * Prx[ma]) * rsq + (dsfac * Pi[ma] * ux + sfac * Pix[ma]) * rinv;
                dedy += (dsfac * Pr[ma] * uy + sfac * Pry[ma]) * rsq + (dsfac * Pi[ma] * uy + sfac * Piy[ma]) * rinv;
                dedz += (dsfac * Pr[ma] * uz + sfac * Prz[ma]) * rsq + (dsfac * Pi[ma] * uz + sfac * Piz[ma]) * rinv;
            }
            else {
                for (ma = 0; ma <= j; ma++) {
                    //T y_r = ylist_r[in];
                    //T y_i = ylist_i[in];
                    rsq = ylist_r[i + inum*k]; 
                    rinv = ylist_i[i + inum*k];                                        
                    dedx += (dsfac * Pr[ma] * ux + sfac * Prx[ma]) * rsq + (dsfac * Pi[ma] * ux + sfac * Pix[ma]) * rinv;
                    dedy += (dsfac * Pr[ma] * uy + sfac * Pry[ma]) * rsq + (dsfac * Pi[ma] * uy + sfac * Piy[ma]) * rinv;
                    dedz += (dsfac * Pr[ma] * uz + sfac * Prz[ma]) * rsq + (dsfac * Pi[ma] * uz + sfac * Piz[ma]) * rinv;
                    k += 1; 
                }                
            }            
        }
    }      
//     dedr[0 + 3*ij] = 2.0*dedx;      
//     dedr[1 + 3*ij] = 2.0*dedy;      
//     dedr[2 + 3*ij] = 2.0*dedz;         

    dedx = 2.0*dedx;      
    dedy = 2.0*dedy;      
    dedz = 2.0*dedz;             
    ma = ai[ij];        
    mb = aj[ij];        
    fatom[0+3*ma] += dedx;
    fatom[1+3*ma] += dedy;
    fatom[2+3*ma] += dedz;
    fatom[0+3*mb] -= dedx;
    fatom[1+3*mb] -= dedy;
    fatom[2+3*mb] -= dedz;    
  }   
}
template void cpuComputeDeidrj(double*, double*, double*, double*, double*,  double*, double*, 
         double, double, double, int*, int*, int*, int*, int*, int*, int, int, int, int, int, int);
template void cpuComputeDeidrj(float*, float*, float*, float*, float*, float*, float*, 
        float, float, float, int*, int*, int*, int*, int*, int*, int, int, int, int, int, int);

// template <typename T> void cpuComputeDeidrj(T *dedr, T *Stotr, T *Stoti, T *dulist_r, T *dulist_i, T *ylist_r, T *ylist_i, T *rootpqarray, T *rij, 
//         T *wjelem, T *radelem,  T rmin0, T rfac0, T rcutfac, int *idxu_block, int *map, int *ai, int *aj, 
//         int *ti, int *tj, int twojmax, int idxu_max, int inum, int ijnum, int switchflag, int chemflag) 
// {               
//   cpuArraySetValue(dedr, (T) 0.0, ijnum*3);
//   
//   for(int ij=0; ij<ijnum; ij++) {        
//     T x = rij[ij*3+0];
//     T y = rij[ij*3+1];
//     T z = rij[ij*3+2];    
//     T rsq = x * x + y * y + z * z;
//     T r = sqrt(rsq);
//     T rinv = 1.0 / r;
//     T ux = x * rinv;
//     T uy = y * rinv;
//     T uz = z * rinv;
// 
//     T rcutij = (radelem[ti[ij]]+radelem[tj[ij]])*rcutfac; //(radelem[type[ii]]+radelem[type[jj]])*rcutfac;
//     T z0 = r / tan((r - rmin0) * rfac0 * M_PI / (rcutij - rmin0));                
//     T dz0dr = z0 / r - (r*rfac0 * M_PI / (rcutij - rmin0)) * (rsq + z0 * z0) / rsq;
// 
//     T sfac = 0.0, dsfac = 0.0;        
//     if (switchflag == 0) {
//         sfac = 1.0;
//         dsfac = 0.0;
//     }
//     else if (switchflag == 1) {
//         if (r <= rmin0) {
//             sfac = 1.0;
//             dsfac = 0.0;
//         }
//         else if(r > rcutij) {
//             sfac = 1.0;
//             dsfac = 0.0;
//         }
//         else {
//             T rcutfac = M_PI / (rcutij - rmin0);
//             sfac =  0.5 * (cos((r - rmin0) * rcutfac) + 1.0);   
//             dsfac = -0.5 * sin((r - rmin0) * rcutfac) * rcutfac;
//         }
//     } 
//     sfac *= wjelem[tj[ij]];
//     dsfac *= wjelem[tj[ij]];
//     
//     T a_r, a_i, b_r, b_i, rootpq;        
//     rcutij = 1.0 / sqrt(r * r + z0 * z0);
//     a_r = rcutij * z0;
//     a_i = -rcutij * z;
//     b_r = rcutij * y;
//     b_i = -rcutij * x;
// 
//     T u_r, u_i, ux_r, ux_i, uy_r, uy_i, uz_r, uz_i;
//     T w_r, w_i, wx_r, wx_i, wy_r, wy_i, wz_r, wz_i;
//     u_r = -pow(rcutij, 3.0) * (r + z0 * dz0dr);
//     wx_r = u_r * ux;
//     wy_r = u_r * uy;
//     wz_r = u_r * uz;
//     ux_r = dz0dr * ux;
//     uy_r = dz0dr * uy;
//     uz_r = dz0dr * uz;
// 
//     T dardx, daidx, dardy, daidy, dardz, daidz;
//     dardx = ux_r * rcutij + z0 * wx_r;
//     daidx = -z * wx_r;
//     dardy = uy_r * rcutij + z0 * wy_r;
//     daidy = -z * wy_r;
//     dardz = uz_r * rcutij + z0 * wz_r;
//     daidz = -z * wz_r;    
//     daidz += -rcutij;
// 
//     T dbrdx, dbidx, dbrdy, dbidy, dbrdz, dbidz;
//     dbrdx = y * wx_r;
//     dbidx = -x * wx_r;    
//     dbrdy = y * wy_r;
//     dbidy = -x * wy_r;    
//     dbrdz = y * wz_r;
//     dbidz = -x * wz_r;        
//     dbidx += -rcutij;
//     dbrdy += rcutij;
//     
//     // 2Jmax = 10    
//     T Pr[11], Pi[11], Qr[9], Qi[9];
//     T Prx[11], Pix[11], Qrx[9], Qix[9];
//     T Pry[11], Piy[11], Qry[9], Qiy[9];        
//     T Prz[11], Piz[11], Qrz[9], Qiz[9];
//     Pr[0] = 1.0; Pi[0] = 0.0;    
//     Prx[0] = 0.0; Pix[0] = 0.0;    
//     Pry[0] = 0.0; Piy[0] = 0.0;    
//     Prz[0] = 0.0; Piz[0] = 0.0;        
//     
//     int jdim = twojmax + 1;
//     int njelem = (chemflag==0) ? 0 : (inum*idxu_max*map[tj[ij]]);    
//     int i = ai[ij] + njelem;                        
//                     
//     //T dedx, dedy, dedz;       
//     x = 0.0; y = 0.0; z = 0.0;                
//     u_r = 0.5*ylist_r[i];    
//     x += (dsfac * ux) * u_r;
//     y += (dsfac * uy) * u_r;
//     z += (dsfac * uz) * u_r; 
//             
//     int j, k, ma, mb, mapar;    
//     mb = 0;
//     for (j = 1; j <= twojmax; j++) {        
//         // fill in left side of matrix layer from previous layer
//         ma = 0;
//         u_r = Pr[ma];
//         u_i = Pi[ma];
//         ux_r = Prx[ma];
//         ux_i = Pix[ma];            
//         uy_r = Pry[ma];
//         uy_i = Piy[ma];            
//         uz_r = Prz[ma];
//         uz_i = Piz[ma];                    
//         rootpq = rootpqarray[(j - ma)*jdim + (j - mb)];
//         Pr[ma] = rootpq * (a_r * u_r + a_i * u_i);
//         Pi[ma] = rootpq * (a_r * u_i - a_i * u_r);        
//         Prx[ma] = rootpq * (dardx * u_r + daidx * u_i + a_r * ux_r + a_i * ux_i);
//         Pix[ma] = rootpq * (dardx * u_i - daidx * u_r + a_r * ux_i - a_i * ux_r);
//         Pry[ma] = rootpq * (dardy * u_r + daidy * u_i + a_r * uy_r + a_i * uy_i);
//         Piy[ma] = rootpq * (dardy * u_i - daidy * u_r + a_r * uy_i - a_i * uy_r);
//         Prz[ma] = rootpq * (dardz * u_r + daidz * u_i + a_r * uz_r + a_i * uz_i);
//         Piz[ma] = rootpq * (dardz * u_i - daidz * u_r + a_r * uz_i - a_i * uz_r);                    
//         for (ma = 1; ma < j; ma++) {
//             w_r = Pr[ma];
//             w_i = Pi[ma];
//             wx_r = Prx[ma];
//             wx_i = Pix[ma];            
//             wy_r = Pry[ma];
//             wy_i = Piy[ma];            
//             wz_r = Prz[ma];
//             wz_i = Piz[ma];                                        
//             rootpq = rootpqarray[(j - ma)*jdim + (j - mb)];
//             rcutij = rootpqarray[ma*jdim + (j - mb)];
//             Pr[ma] = rootpq * (a_r * w_r + a_i * w_i) -rcutij * (b_r * u_r + b_i * u_i);
//             Pi[ma] = rootpq * (a_r * w_i - a_i * w_r) -rcutij * (b_r * u_i - b_i * u_r);
//             Prx[ma] = rootpq * (dardx * w_r + daidx * w_i + a_r * wx_r + a_i * wx_i) -rcutij * (dbrdx * u_r + dbidx * u_i + b_r * ux_r + b_i * ux_i);
//             Pix[ma] = rootpq * (dardx * w_i - daidx * w_r + a_r * wx_i - a_i * wx_r) -rcutij * (dbrdx * u_i - dbidx * u_r + b_r * ux_i - b_i * ux_r);
//             Pry[ma] = rootpq * (dardy * w_r + daidy * w_i + a_r * wy_r + a_i * wy_i) -rcutij * (dbrdy * u_r + dbidy * u_i + b_r * uy_r + b_i * uy_i);
//             Piy[ma] = rootpq * (dardy * w_i - daidy * w_r + a_r * wy_i - a_i * wy_r) -rcutij * (dbrdy * u_i - dbidy * u_r + b_r * uy_i - b_i * uy_r);
//             Prz[ma] = rootpq * (dardz * w_r + daidz * w_i + a_r * wz_r + a_i * wz_i) -rcutij * (dbrdz * u_r + dbidz * u_i + b_r * uz_r + b_i * uz_i);
//             Piz[ma] = rootpq * (dardz * w_i - daidz * w_r + a_r * wz_i - a_i * wz_r) -rcutij * (dbrdz * u_i - dbidz * u_r + b_r * uz_i - b_i * uz_r);            
//             u_r = w_r;
//             u_i = w_i;
//             ux_r = wx_r;
//             ux_i = wx_i;
//             uy_r = wy_r;
//             uy_i = wy_i;
//             uz_r = wz_r;
//             uz_i = wz_i;            
//         }
//         ma = j;
//         rcutij = rootpqarray[ma*jdim + (j - mb)];
//         Pr[ma] = -rcutij * (b_r * u_r + b_i * u_i);
//         Pi[ma] = -rcutij * (b_r * u_i - b_i * u_r);                        
//         Prx[ma] =-rcutij * (dbrdx * u_r + dbidx * u_i + b_r * ux_r + b_i * ux_i);
//         Pix[ma] =-rcutij * (dbrdx * u_i - dbidx * u_r + b_r * ux_i - b_i * ux_r);
//         Pry[ma] =-rcutij * (dbrdy * u_r + dbidy * u_i + b_r * uy_r + b_i * uy_i);
//         Piy[ma] =-rcutij * (dbrdy * u_i - dbidy * u_r + b_r * uy_i - b_i * uy_r);
//         Prz[ma] =-rcutij * (dbrdz * u_r + dbidz * u_i + b_r * uz_r + b_i * uz_i);
//         Piz[ma] =-rcutij * (dbrdz * u_i - dbidz * u_r + b_r * uz_i - b_i * uz_r);
//                                 
//         if (j==(2*mb+1)) { // store Qr, Qi, for the next mb level
//             mapar = (((j+1)/2)%2 == 0) ? -1 : 1;
//             for (ma = 0; ma <= j; ma++) {     
//                 if (mapar == 1) {                    
//                     Qr[j-ma] = Pr[ma];
//                     Qi[j-ma] = -Pi[ma];
//                     Qrx[j-ma] =  Prx[ma];
//                     Qix[j-ma] = -Pix[ma];
//                     Qry[j-ma] =  Pry[ma];
//                     Qiy[j-ma] = -Piy[ma];
//                     Qrz[j-ma] =  Prz[ma];
//                     Qiz[j-ma] = -Piz[ma];                    
//                 } else {
//                     Qr[j-ma] = -Pr[ma];
//                     Qi[j-ma] =  Pi[ma];
//                     Qrx[j-ma] = -Prx[ma];
//                     Qix[j-ma] =  Pix[ma];
//                     Qry[j-ma] = -Pry[ma];
//                     Qiy[j-ma] =  Piy[ma];
//                     Qrz[j-ma] = -Prz[ma];
//                     Qiz[j-ma] =  Piz[ma];                    
//                 }
//                 mapar = -mapar;
//             }                              
//         }
//         
//         k =  1 + (j+1)*mb;
//         for (ma = 2; ma <= j; ma++)
//             k += ma*ma;                    
//         for (ma = 0; ma <= j; ma++) {                            
//             //T y_r = ylist_r[in]; 
//             //T y_i = ylist_i[in]; 
//             rsq = ylist_r[i + inum*k]; 
//             rinv = ylist_i[i + inum*k];
//             x += (dsfac * Pr[ma] * ux + sfac * Prx[ma]) * rsq + (dsfac * Pi[ma] * ux + sfac * Pix[ma]) * rinv;
//             y += (dsfac * Pr[ma] * uy + sfac * Pry[ma]) * rsq + (dsfac * Pi[ma] * uy + sfac * Piy[ma]) * rinv;
//             z += (dsfac * Pr[ma] * uz + sfac * Prz[ma]) * rsq + (dsfac * Pi[ma] * uz + sfac * Piz[ma]) * rinv;            
//             k += 1;
//         }                   
//     }
//     
//     for (mb = 1; 2*mb <= twojmax; mb++) {     
//         for (ma = 0; ma < 2*mb; ma++) {                      
//             Pr[ma] = Qr[ma];
//             Pi[ma] = Qi[ma];
//             Prx[ma] = Qrx[ma];
//             Pix[ma] = Qix[ma];
//             Pry[ma] = Qry[ma];
//             Piy[ma] = Qiy[ma];
//             Prz[ma] = Qrz[ma];
//             Piz[ma] = Qiz[ma];            
//         }                
//         for (j = 2*mb; j <= twojmax; j++) { 
//             ma = 0;
//             u_r = Pr[ma];
//             u_i = Pi[ma];
//             ux_r = Prx[ma];
//             ux_i = Pix[ma];            
//             uy_r = Pry[ma];
//             uy_i = Piy[ma];            
//             uz_r = Prz[ma];
//             uz_i = Piz[ma];                                
//             rootpq = rootpqarray[(j - ma)*jdim + (j - mb)];
//             Pr[ma] = rootpq * (a_r * u_r + a_i * u_i);
//             Pi[ma] = rootpq * (a_r * u_i - a_i * u_r);            
//             Prx[ma] = rootpq * (dardx * u_r + daidx * u_i + a_r * ux_r + a_i * ux_i);
//             Pix[ma] = rootpq * (dardx * u_i - daidx * u_r + a_r * ux_i - a_i * ux_r);
//             Pry[ma] = rootpq * (dardy * u_r + daidy * u_i + a_r * uy_r + a_i * uy_i);
//             Piy[ma] = rootpq * (dardy * u_i - daidy * u_r + a_r * uy_i - a_i * uy_r);
//             Prz[ma] = rootpq * (dardz * u_r + daidz * u_i + a_r * uz_r + a_i * uz_i);
//             Piz[ma] = rootpq * (dardz * u_i - daidz * u_r + a_r * uz_i - a_i * uz_r);                                
//             for (ma = 1; ma < j; ma++) {
//                 w_r = Pr[ma];
//                 w_i = Pi[ma];
//                 wx_r = Prx[ma];
//                 wx_i = Pix[ma];            
//                 wy_r = Pry[ma];
//                 wy_i = Piy[ma];            
//                 wz_r = Prz[ma];
//                 wz_i = Piz[ma];                                            
//                 rootpq = rootpqarray[(j - ma)*jdim + (j - mb)];
//                 rcutij = rootpqarray[ma*jdim + (j - mb)];
//                 Pr[ma] = rootpq * (a_r * w_r + a_i * w_i) -rcutij * (b_r * u_r + b_i * u_i);
//                 Pi[ma] = rootpq * (a_r * w_i - a_i * w_r) -rcutij * (b_r * u_i - b_i * u_r);
//                 Prx[ma] = rootpq * (dardx * w_r + daidx * w_i + a_r * wx_r + a_i * wx_i) -rcutij * (dbrdx * u_r + dbidx * u_i + b_r * ux_r + b_i * ux_i);
//                 Pix[ma] = rootpq * (dardx * w_i - daidx * w_r + a_r * wx_i - a_i * wx_r) -rcutij * (dbrdx * u_i - dbidx * u_r + b_r * ux_i - b_i * ux_r);
//                 Pry[ma] = rootpq * (dardy * w_r + daidy * w_i + a_r * wy_r + a_i * wy_i) -rcutij * (dbrdy * u_r + dbidy * u_i + b_r * uy_r + b_i * uy_i);
//                 Piy[ma] = rootpq * (dardy * w_i - daidy * w_r + a_r * wy_i - a_i * wy_r) -rcutij * (dbrdy * u_i - dbidy * u_r + b_r * uy_i - b_i * uy_r);
//                 Prz[ma] = rootpq * (dardz * w_r + daidz * w_i + a_r * wz_r + a_i * wz_i) -rcutij * (dbrdz * u_r + dbidz * u_i + b_r * uz_r + b_i * uz_i);
//                 Piz[ma] = rootpq * (dardz * w_i - daidz * w_r + a_r * wz_i - a_i * wz_r) -rcutij * (dbrdz * u_i - dbidz * u_r + b_r * uz_i - b_i * uz_r);                            
//                 u_r = w_r;
//                 u_i = w_i;
//                 ux_r = wx_r;
//                 ux_i = wx_i;
//                 uy_r = wy_r;
//                 uy_i = wy_i;
//                 uz_r = wz_r;
//                 uz_i = wz_i;                
//             }
//             ma = j;
//             rcutij = rootpqarray[ma*jdim + (j - mb)];
//             Pr[ma] = -rcutij * (b_r * u_r + b_i * u_i);
//             Pi[ma] = -rcutij * (b_r * u_i - b_i * u_r);       
//             Prx[ma] =-rcutij * (dbrdx * u_r + dbidx * u_i + b_r * ux_r + b_i * ux_i);
//             Pix[ma] =-rcutij * (dbrdx * u_i - dbidx * u_r + b_r * ux_i - b_i * ux_r);
//             Pry[ma] =-rcutij * (dbrdy * u_r + dbidy * u_i + b_r * uy_r + b_i * uy_i);
//             Piy[ma] =-rcutij * (dbrdy * u_i - dbidy * u_r + b_r * uy_i - b_i * uy_r);
//             Prz[ma] =-rcutij * (dbrdz * u_r + dbidz * u_i + b_r * uz_r + b_i * uz_i);
//             Piz[ma] =-rcutij * (dbrdz * u_i - dbidz * u_r + b_r * uz_i - b_i * uz_r);
//             
//             if (j==(2*mb)) {
//                 mapar = 1;
//                 for (ma = 0; ma <= j/2; ma++) {
//                     if (mapar == 1) {                    
//                         Pr[j/2+ma] = Pr[j/2-ma];
//                         Pi[j/2+ma] = -Pi[j/2-ma];
//                         Prx[j/2+ma] = Prx[j/2-ma];
//                         Pix[j/2+ma] = -Pix[j/2-ma];
//                         Pry[j/2+ma] = Pry[j/2-ma];
//                         Piy[j/2+ma] = -Piy[j/2-ma];
//                         Prz[j/2+ma] = Prz[j/2-ma];
//                         Piz[j/2+ma] = -Piz[j/2-ma];                        
//                     } else {
//                         Pr[j/2+ma] = -Pr[j/2-ma];
//                         Pi[j/2+ma] = Pi[j/2-ma];
//                         Prx[j/2+ma] = -Prx[j/2-ma];
//                         Pix[j/2+ma] =  Pix[j/2-ma];
//                         Pry[j/2+ma] = -Pry[j/2-ma];
//                         Piy[j/2+ma] =  Piy[j/2-ma];
//                         Prz[j/2+ma] = -Prz[j/2-ma];
//                         Piz[j/2+ma] =  Piz[j/2-ma];                        
//                     }
//                     mapar = -mapar;        
//                 }                                                        
//             }
//             
//             if (j==(2*mb)) {
//                 mapar = 1;
//                 for (ma = 0; ma <= j; ma++) {
//                     if (mapar == 1) {                    
//                         Prx[j/2+ma] = Prx[j/2-ma];
//                         Pix[j/2+ma] = -Pix[j/2-ma];
//                         Pry[j/2+ma] = Pry[j/2-ma];
//                         Piy[j/2+ma] = -Piy[j/2-ma];
//                         Prz[j/2+ma] = Prz[j/2-ma];
//                         Piz[j/2+ma] = -Piz[j/2-ma];                        
//                     } else {
//                         Prx[j/2+ma] = -Prx[j/2-ma];
//                         Pix[j/2+ma] =  Pix[j/2-ma];
//                         Pry[j/2+ma] = -Pry[j/2-ma];
//                         Piy[j/2+ma] =  Piy[j/2-ma];
//                         Prz[j/2+ma] = -Prz[j/2-ma];
//                         Piz[j/2+ma] =  Piz[j/2-ma];                        
//                     }
//                     mapar = -mapar;        
//                 }                                                        
//             }
//                         
//             // store Qr, Qi, for the next mb level
//             if (j==(2*mb+1)) {
//                 mapar = (((j+1)/2)%2 == 0) ? -1 : 1;
//                 for (ma = 0; ma <= j; ma++) {     
//                     if (mapar == 1) {                    
//                         Qr[j-ma] = Pr[ma];
//                         Qi[j-ma] = -Pi[ma];  
//                         Qrx[j-ma] =  Prx[ma];
//                         Qix[j-ma] = -Pix[ma];
//                         Qry[j-ma] =  Pry[ma];
//                         Qiy[j-ma] = -Piy[ma];
//                         Qrz[j-ma] =  Prz[ma];
//                         Qiz[j-ma] = -Piz[ma];                                            
//                     } else {
//                         Qr[j-ma] = -Pr[ma];
//                         Qi[j-ma] =  Pi[ma];
//                         Qrx[j-ma] = -Prx[ma];
//                         Qix[j-ma] =  Pix[ma];
//                         Qry[j-ma] = -Pry[ma];
//                         Qiy[j-ma] =  Piy[ma];
//                         Qrz[j-ma] = -Prz[ma];
//                         Qiz[j-ma] =  Piz[ma];                                            
//                     }
//                     mapar = -mapar;
//                 }                                                
//             }
//             
//             k =  1 + (j+1)*mb;
//             for (ma = 2; ma <= j; ma++)
//                 k += ma*ma;                            
//             if (j==(2*mb)) {
//                 for (ma = 0; ma < mb; ma++) {
//                     //T y_r = ylist_r[in];
//                     //T y_i = ylist_i[in];
//                     rsq = ylist_r[i + inum*k]; 
//                     rinv = ylist_i[i + inum*k];                    
//                     x += (dsfac * Pr[ma] * ux + sfac * Prx[ma]) * rsq + (dsfac * Pi[ma] * ux + sfac * Pix[ma]) * rinv;
//                     y += (dsfac * Pr[ma] * uy + sfac * Pry[ma]) * rsq + (dsfac * Pi[ma] * uy + sfac * Piy[ma]) * rinv;
//                     z += (dsfac * Pr[ma] * uz + sfac * Prz[ma]) * rsq + (dsfac * Pi[ma] * uz + sfac * Piz[ma]) * rinv;
//                     k += 1; 
//                 }     
//                 ma = mb;
//                 //T y_r = 0.5*ylist_r[in];
//                 //T y_i = 0.5*ylist_i[in];
//                 rsq = 0.5*ylist_r[i + inum*k]; 
//                 rinv = 0.5*ylist_i[i + inum*k];                
//                 x += (dsfac * Pr[ma] * ux + sfac * Prx[ma]) * rsq + (dsfac * Pi[ma] * ux + sfac * Pix[ma]) * rinv;
//                 y += (dsfac * Pr[ma] * uy + sfac * Pry[ma]) * rsq + (dsfac * Pi[ma] * uy + sfac * Piy[ma]) * rinv;
//                 z += (dsfac * Pr[ma] * uz + sfac * Prz[ma]) * rsq + (dsfac * Pi[ma] * uz + sfac * Piz[ma]) * rinv;
//             }
//             else {
//                 for (ma = 0; ma <= j; ma++) {
//                     //T y_r = ylist_r[in];
//                     //T y_i = ylist_i[in];
//                     rsq = ylist_r[i + inum*k]; 
//                     rinv = ylist_i[i + inum*k];                                        
//                     x += (dsfac * Pr[ma] * ux + sfac * Prx[ma]) * rsq + (dsfac * Pi[ma] * ux + sfac * Pix[ma]) * rinv;
//                     y += (dsfac * Pr[ma] * uy + sfac * Pry[ma]) * rsq + (dsfac * Pi[ma] * uy + sfac * Piy[ma]) * rinv;
//                     z += (dsfac * Pr[ma] * uz + sfac * Prz[ma]) * rsq + (dsfac * Pi[ma] * uz + sfac * Piz[ma]) * rinv;
//                     k += 1; 
//                 }                
//             }            
//         }
//     }
//       
//     dedr[0 + 3*ij] = 2.0*x;      
//     dedr[1 + 3*ij] = 2.0*y;      
//     dedr[2 + 3*ij] = 2.0*z;         
//   }   
// }
// // template void cpuComputeDeidrj(double*, double*, double*, double*, double*,  double*, double*, 
// //          double, double, double, int*, int*, int*, int*, int*, int*, int, int, int, int, int, int);
// // template void cpuComputeDeidrj(float*, float*, float*, float*, float*, float*, float*, 
// //         float, float, float, int*, int*, int*, int*, int*, int*, int, int, int, int, int, int);
// template void cpuComputeDeidrj(double*, double*, double*, double*, double*, double*, double*,  double*, double*, 
//          double*, double*, double, double, double, int*, int*, int*, int*, int*, int*, int, int, int, int, int, int);
// template void cpuComputeDeidrj(float*, float*, float*, float*, float*, float*, float*, float*, float*, 
//         float*, float*, float, float, float, int*, int*, int*, int*, int*, int*, int, int, int, int, int, int);


// template <typename T> void cpuComputeDeidrj(T *dedr, T *Stotr, T *Stoti, T *dulist_r, T *dulist_i, T *ylist_r, T *ylist_i, T *rootpqarray, T *rij, 
//         T *wjelem, T *radelem,  T rmin0, T rfac0, T rcutfac, int *idxu_block, int *map, int *ai, int *aj, 
//         int *ti, int *tj, int twojmax, int idxu_max, int inum, int ijnum, int switchflag, int chemflag) 
// {               
//   cpuArraySetValue(dedr, (T) 0.0, ijnum*3);
//   
//   for(int ij=0; ij<ijnum; ij++) {        
//     T x = rij[ij*3+0];
//     T y = rij[ij*3+1];
//     T z = rij[ij*3+2];    
//     T rsq = x * x + y * y + z * z;
//     T r = sqrt(rsq);
//     T rinv = 1.0 / r;
//     T ux = x * rinv;
//     T uy = y * rinv;
//     T uz = z * rinv;
// 
//     T rcutij = (radelem[ti[ij]]+radelem[tj[ij]])*rcutfac; //(radelem[type[ii]]+radelem[type[jj]])*rcutfac;
//     //T rscale0 = rfac0 * M_PI / (rcutij - rmin0);
//     //T theta0 = (r - rmin0) * rscale0;
//     T z0 = r / tan((r - rmin0) * rfac0 * M_PI / (rcutij - rmin0));                
//     T dz0dr = z0 / r - (r*rfac0 * M_PI / (rcutij - rmin0)) * (rsq + z0 * z0) / rsq;
// 
//     T sfac = 0.0, dsfac = 0.0;        
//     if (switchflag == 0) {
//         sfac = 1.0;
//         dsfac = 0.0;
//     }
//     else if (switchflag == 1) {
//         if (r <= rmin0) {
//             sfac = 1.0;
//             dsfac = 0.0;
//         }
//         else if(r > rcutij) {
//             sfac = 1.0;
//             dsfac = 0.0;
//         }
//         else {
//             T rcutfac = M_PI / (rcutij - rmin0);
//             sfac =  0.5 * (cos((r - rmin0) * rcutfac) + 1.0);   
//             dsfac = -0.5 * sin((r - rmin0) * rcutfac) * rcutfac;
//         }
//     } 
//     sfac *= wjelem[tj[ij]];
//     dsfac *= wjelem[tj[ij]];
//                 
//     //sfac = 1.0; 
//     //dsfac = 0.0;
//     
//     T r0inv, dr0invdr;
//     T a_r, a_i, b_r, b_i;
//     T da_r[3], da_i[3], db_r[3], db_i[3];
//     T dz0[3], dr0inv[3];
//     T rootpq;       
// 
//     r0inv = 1.0 / sqrt(r * r + z0 * z0);
//     a_r = r0inv * z0;
//     a_i = -r0inv * z;
//     b_r = r0inv * y;
//     b_i = -r0inv * x;
// 
//     dr0invdr = -pow(r0inv, 3.0) * (r + z0 * dz0dr);
// 
//     dr0inv[0] = dr0invdr * ux;
//     dr0inv[1] = dr0invdr * uy;
//     dr0inv[2] = dr0invdr * uz;
// 
//     dz0[0] = dz0dr * ux;
//     dz0[1] = dz0dr * uy;
//     dz0[2] = dz0dr * uz;
// 
//     for (int k = 0; k < 3; k++) {
//         da_r[k] = dz0[k] * r0inv + z0 * dr0inv[k];
//         da_i[k] = -z * dr0inv[k];
//     }
//     da_i[2] += -r0inv;
// 
//     for (int k = 0; k < 3; k++) {
//         db_r[k] = y * dr0inv[k];
//         db_i[k] = -x * dr0inv[k];
//     }
//     db_i[0] += -r0inv;
//     db_r[1] += r0inv;
//     
//     // 2Jmax = 10    
//     T u_r, u_i, ux_r, ux_i, uy_r, uy_i, uz_r, uz_i;
//     T w_r, w_i, wx_r, wx_i, wy_r, wy_i, wz_r, wz_i;
//     T Pr[11], Pi[11], Qr[9], Qi[9];
//     T Prx[11], Pix[11], Qrx[9], Qix[9];
//     T Pry[11], Piy[11], Qry[9], Qiy[9];        
//     T Prz[11], Piz[11], Qrz[9], Qiz[9];
//     Pr[0] = 1.0; Pi[0] = 0.0;    
//     Prx[0] = 0.0; Pix[0] = 0.0;    
//     Pry[0] = 0.0; Piy[0] = 0.0;    
//     Prz[0] = 0.0; Piz[0] = 0.0;    
//     T de[3];       
//     de[0] = 0.0; de[1] = 0.0; de[2] = 0.0;                
//     
//     int jdim = twojmax + 1;
//     int njelem = (chemflag==0) ? 0 : (inum*idxu_max*map[tj[ij]]);    
//     int i = ai[ij] + njelem;                        
//     Stotr[i] += sfac; // atomic add               
//     
//     dulist_r[ij+ijnum*0+ijnum*idxu_max*0] = dsfac * 1.0 * ux;
//     dulist_i[ij+ijnum*0+ijnum*idxu_max*0] = dsfac * 0.0 * ux;
//     dulist_r[ij+ijnum*0+ijnum*idxu_max*1] = dsfac * 1.0 * uy;
//     dulist_i[ij+ijnum*0+ijnum*idxu_max*1] = dsfac * 0.0 * uy;
//     dulist_r[ij+ijnum*0+ijnum*idxu_max*2] = dsfac * 1.0 * uz;
//     dulist_i[ij+ijnum*0+ijnum*idxu_max*2] = dsfac * 0.0 * uz;                        
//                 
//     T y_r = 0.5*ylist_r[i];    
//     de[0] += (dsfac * 1.0 * ux) * y_r;
//     de[1] += (dsfac * 1.0 * uy) * y_r;
//     de[2] += (dsfac * 1.0 * uz) * y_r; 
//             
//     int mb = 0;    
//     for (int j = 1; j <= twojmax; j++) {        
//         // fill in left side of matrix layer from previous layer
//         int ma = 0;
//         u_r = Pr[ma];
//         u_i = Pi[ma];
//         ux_r = Prx[ma];
//         ux_i = Pix[ma];            
//         uy_r = Pry[ma];
//         uy_i = Piy[ma];            
//         uz_r = Prz[ma];
//         uz_i = Piz[ma];                    
//         rootpq = rootpqarray[(j - ma)*jdim + (j - mb)];
//         Pr[ma] = rootpq * (a_r * u_r + a_i * u_i);
//         Pi[ma] = rootpq * (a_r * u_i - a_i * u_r);        
//         Prx[ma] = rootpq * (da_r[0] * u_r + da_i[0] * u_i + a_r * ux_r + a_i * ux_i);
//         Pix[ma] = rootpq * (da_r[0] * u_i - da_i[0] * u_r + a_r * ux_i - a_i * ux_r);
//         Pry[ma] = rootpq * (da_r[1] * u_r + da_i[1] * u_i + a_r * uy_r + a_i * uy_i);
//         Piy[ma] = rootpq * (da_r[1] * u_i - da_i[1] * u_r + a_r * uy_i - a_i * uy_r);
//         Prz[ma] = rootpq * (da_r[2] * u_r + da_i[2] * u_i + a_r * uz_r + a_i * uz_i);
//         Piz[ma] = rootpq * (da_r[2] * u_i - da_i[2] * u_r + a_r * uz_i - a_i * uz_r);                    
//         for (ma = 1; ma < j; ma++) {
//             w_r = Pr[ma];
//             w_i = Pi[ma];
//             wx_r = Prx[ma];
//             wx_i = Pix[ma];            
//             wy_r = Pry[ma];
//             wy_i = Piy[ma];            
//             wz_r = Prz[ma];
//             wz_i = Piz[ma];                                        
//             rootpq = rootpqarray[(j - ma)*jdim + (j - mb)];
//             rcutij = rootpqarray[ma*jdim + (j - mb)];
//             Pr[ma] = rootpq * (a_r * w_r + a_i * w_i) -rcutij * (b_r * u_r + b_i * u_i);
//             Pi[ma] = rootpq * (a_r * w_i - a_i * w_r) -rcutij * (b_r * u_i - b_i * u_r);
//             Prx[ma] = rootpq * (da_r[0] * w_r + da_i[0] * w_i + a_r * wx_r + a_i * wx_i) -rcutij * (db_r[0] * u_r + db_i[0] * u_i + b_r * ux_r + b_i * ux_i);
//             Pix[ma] = rootpq * (da_r[0] * w_i - da_i[0] * w_r + a_r * wx_i - a_i * wx_r) -rcutij * (db_r[0] * u_i - db_i[0] * u_r + b_r * ux_i - b_i * ux_r);
//             Pry[ma] = rootpq * (da_r[1] * w_r + da_i[1] * w_i + a_r * wy_r + a_i * wy_i) -rcutij * (db_r[1] * u_r + db_i[1] * u_i + b_r * uy_r + b_i * uy_i);
//             Piy[ma] = rootpq * (da_r[1] * w_i - da_i[1] * w_r + a_r * wy_i - a_i * wy_r) -rcutij * (db_r[1] * u_i - db_i[1] * u_r + b_r * uy_i - b_i * uy_r);
//             Prz[ma] = rootpq * (da_r[2] * w_r + da_i[2] * w_i + a_r * wz_r + a_i * wz_i) -rcutij * (db_r[2] * u_r + db_i[2] * u_i + b_r * uz_r + b_i * uz_i);
//             Piz[ma] = rootpq * (da_r[2] * w_i - da_i[2] * w_r + a_r * wz_i - a_i * wz_r) -rcutij * (db_r[2] * u_i - db_i[2] * u_r + b_r * uz_i - b_i * uz_r);            
//             u_r = w_r;
//             u_i = w_i;
//             ux_r = wx_r;
//             ux_i = wx_i;
//             uy_r = wy_r;
//             uy_i = wy_i;
//             uz_r = wz_r;
//             uz_i = wz_i;            
//         }
//         ma = j;
//         rcutij = rootpqarray[ma*jdim + (j - mb)];
//         Pr[ma] = -rcutij * (b_r * u_r + b_i * u_i);
//         Pi[ma] = -rcutij * (b_r * u_i - b_i * u_r);                        
//         Prx[ma] =-rcutij * (db_r[0] * u_r + db_i[0] * u_i + b_r * ux_r + b_i * ux_i);
//         Pix[ma] =-rcutij * (db_r[0] * u_i - db_i[0] * u_r + b_r * ux_i - b_i * ux_r);
//         Pry[ma] =-rcutij * (db_r[1] * u_r + db_i[1] * u_i + b_r * uy_r + b_i * uy_i);
//         Piy[ma] =-rcutij * (db_r[1] * u_i - db_i[1] * u_r + b_r * uy_i - b_i * uy_r);
//         Prz[ma] =-rcutij * (db_r[2] * u_r + db_i[2] * u_i + b_r * uz_r + b_i * uz_i);
//         Piz[ma] =-rcutij * (db_r[2] * u_i - db_i[2] * u_r + b_r * uz_i - b_i * uz_r);
//                                 
//         if (j==(2*mb+1)) { // store Qr, Qi, for the next mb level
//             int mapar = (((j+1)/2)%2 == 0) ? -1 : 1;
//             for (int ma = 0; ma <= j; ma++) {     
//                 if (mapar == 1) {                    
//                     Qr[j-ma] = Pr[ma];
//                     Qi[j-ma] = -Pi[ma];
//                     Qrx[j-ma] =  Prx[ma];
//                     Qix[j-ma] = -Pix[ma];
//                     Qry[j-ma] =  Pry[ma];
//                     Qiy[j-ma] = -Piy[ma];
//                     Qrz[j-ma] =  Prz[ma];
//                     Qiz[j-ma] = -Piz[ma];                    
//                 } else {
//                     Qr[j-ma] = -Pr[ma];
//                     Qi[j-ma] =  Pi[ma];
//                     Qrx[j-ma] = -Prx[ma];
//                     Qix[j-ma] =  Pix[ma];
//                     Qry[j-ma] = -Pry[ma];
//                     Qiy[j-ma] =  Piy[ma];
//                     Qrz[j-ma] = -Prz[ma];
//                     Qiz[j-ma] =  Piz[ma];                    
//                 }
//                 mapar = -mapar;
//             }                              
//         }
//         
//         int k =  1 + (j+1)*mb;
//         for (int ma = 2; ma <= j; ma++)
//             k += ma*ma;                    
//         for (int ma = 0; ma <= j; ma++) {
//             int in = i + inum*k;                
//             Stotr[in] += sfac*Pr[ma]; // atomic add   
//             Stoti[in] += sfac*Pi[ma]; // atomic add       
//             T y_r = ylist_r[in];
//             T y_i = ylist_i[in];
//             dulist_r[ij+ijnum*k+ijnum*idxu_max*0] = dsfac * Pr[ma] * ux + sfac * Prx[ma];
//             dulist_i[ij+ijnum*k+ijnum*idxu_max*0] = dsfac * Pi[ma] * ux + sfac * Pix[ma];
//             dulist_r[ij+ijnum*k+ijnum*idxu_max*1] = dsfac * Pr[ma] * uy + sfac * Pry[ma];
//             dulist_i[ij+ijnum*k+ijnum*idxu_max*1] = dsfac * Pi[ma] * uy + sfac * Piy[ma];
//             dulist_r[ij+ijnum*k+ijnum*idxu_max*2] = dsfac * Pr[ma] * uz + sfac * Prz[ma];
//             dulist_i[ij+ijnum*k+ijnum*idxu_max*2] = dsfac * Pi[ma] * uz + sfac * Piz[ma];                        
//             de[0] += (dsfac * Pr[ma] * ux + sfac * Prx[ma]) * y_r + (dsfac * Pi[ma] * ux + sfac * Pix[ma]) * y_i;
//             de[1] += (dsfac * Pr[ma] * uy + sfac * Pry[ma]) * y_r + (dsfac * Pi[ma] * uy + sfac * Piy[ma]) * y_i;
//             de[2] += (dsfac * Pr[ma] * uz + sfac * Prz[ma]) * y_r + (dsfac * Pi[ma] * uz + sfac * Piz[ma]) * y_i;            
//             k += 1;
//         }                   
//     }
//     
//     for (mb = 1; 2*mb <= twojmax; mb++) {     
//         for (int ma = 0; ma < 2*mb; ma++) {                      
//             Pr[ma] = Qr[ma];
//             Pi[ma] = Qi[ma];
//             Prx[ma] = Qrx[ma];
//             Pix[ma] = Qix[ma];
//             Pry[ma] = Qry[ma];
//             Piy[ma] = Qiy[ma];
//             Prz[ma] = Qrz[ma];
//             Piz[ma] = Qiz[ma];            
//         }                
//         for (int j = 2*mb; j <= twojmax; j++) { 
//             int ma = 0;
//             u_r = Pr[ma];
//             u_i = Pi[ma];
//             ux_r = Prx[ma];
//             ux_i = Pix[ma];            
//             uy_r = Pry[ma];
//             uy_i = Piy[ma];            
//             uz_r = Prz[ma];
//             uz_i = Piz[ma];                                
//             rootpq = rootpqarray[(j - ma)*jdim + (j - mb)];
//             Pr[ma] = rootpq * (a_r * u_r + a_i * u_i);
//             Pi[ma] = rootpq * (a_r * u_i - a_i * u_r);            
//             Prx[ma] = rootpq * (da_r[0] * u_r + da_i[0] * u_i + a_r * ux_r + a_i * ux_i);
//             Pix[ma] = rootpq * (da_r[0] * u_i - da_i[0] * u_r + a_r * ux_i - a_i * ux_r);
//             Pry[ma] = rootpq * (da_r[1] * u_r + da_i[1] * u_i + a_r * uy_r + a_i * uy_i);
//             Piy[ma] = rootpq * (da_r[1] * u_i - da_i[1] * u_r + a_r * uy_i - a_i * uy_r);
//             Prz[ma] = rootpq * (da_r[2] * u_r + da_i[2] * u_i + a_r * uz_r + a_i * uz_i);
//             Piz[ma] = rootpq * (da_r[2] * u_i - da_i[2] * u_r + a_r * uz_i - a_i * uz_r);                                
//             for (ma = 1; ma < j; ma++) {
//                 w_r = Pr[ma];
//                 w_i = Pi[ma];
//                 wx_r = Prx[ma];
//                 wx_i = Pix[ma];            
//                 wy_r = Pry[ma];
//                 wy_i = Piy[ma];            
//                 wz_r = Prz[ma];
//                 wz_i = Piz[ma];                                            
//                 rootpq = rootpqarray[(j - ma)*jdim + (j - mb)];
//                 rcutij = rootpqarray[ma*jdim + (j - mb)];
//                 Pr[ma] = rootpq * (a_r * w_r + a_i * w_i) -rcutij * (b_r * u_r + b_i * u_i);
//                 Pi[ma] = rootpq * (a_r * w_i - a_i * w_r) -rcutij * (b_r * u_i - b_i * u_r);
//                 Prx[ma] = rootpq * (da_r[0] * w_r + da_i[0] * w_i + a_r * wx_r + a_i * wx_i) -rcutij * (db_r[0] * u_r + db_i[0] * u_i + b_r * ux_r + b_i * ux_i);
//                 Pix[ma] = rootpq * (da_r[0] * w_i - da_i[0] * w_r + a_r * wx_i - a_i * wx_r) -rcutij * (db_r[0] * u_i - db_i[0] * u_r + b_r * ux_i - b_i * ux_r);
//                 Pry[ma] = rootpq * (da_r[1] * w_r + da_i[1] * w_i + a_r * wy_r + a_i * wy_i) -rcutij * (db_r[1] * u_r + db_i[1] * u_i + b_r * uy_r + b_i * uy_i);
//                 Piy[ma] = rootpq * (da_r[1] * w_i - da_i[1] * w_r + a_r * wy_i - a_i * wy_r) -rcutij * (db_r[1] * u_i - db_i[1] * u_r + b_r * uy_i - b_i * uy_r);
//                 Prz[ma] = rootpq * (da_r[2] * w_r + da_i[2] * w_i + a_r * wz_r + a_i * wz_i) -rcutij * (db_r[2] * u_r + db_i[2] * u_i + b_r * uz_r + b_i * uz_i);
//                 Piz[ma] = rootpq * (da_r[2] * w_i - da_i[2] * w_r + a_r * wz_i - a_i * wz_r) -rcutij * (db_r[2] * u_i - db_i[2] * u_r + b_r * uz_i - b_i * uz_r);                            
//                 u_r = w_r;
//                 u_i = w_i;
//                 ux_r = wx_r;
//                 ux_i = wx_i;
//                 uy_r = wy_r;
//                 uy_i = wy_i;
//                 uz_r = wz_r;
//                 uz_i = wz_i;                
//             }
//             ma = j;
//             rcutij = rootpqarray[ma*jdim + (j - mb)];
//             Pr[ma] = -rcutij * (b_r * u_r + b_i * u_i);
//             Pi[ma] = -rcutij * (b_r * u_i - b_i * u_r);       
//             Prx[ma] =-rcutij * (db_r[0] * u_r + db_i[0] * u_i + b_r * ux_r + b_i * ux_i);
//             Pix[ma] =-rcutij * (db_r[0] * u_i - db_i[0] * u_r + b_r * ux_i - b_i * ux_r);
//             Pry[ma] =-rcutij * (db_r[1] * u_r + db_i[1] * u_i + b_r * uy_r + b_i * uy_i);
//             Piy[ma] =-rcutij * (db_r[1] * u_i - db_i[1] * u_r + b_r * uy_i - b_i * uy_r);
//             Prz[ma] =-rcutij * (db_r[2] * u_r + db_i[2] * u_i + b_r * uz_r + b_i * uz_i);
//             Piz[ma] =-rcutij * (db_r[2] * u_i - db_i[2] * u_r + b_r * uz_i - b_i * uz_r);
//             
//             if (j==(2*mb)) {
//                 int mapar = 1;
//                 for (int ma = 0; ma <= j/2; ma++) {
//                     if (mapar == 1) {                    
//                         Pr[j/2+ma] = Pr[j/2-ma];
//                         Pi[j/2+ma] = -Pi[j/2-ma];
//                         Prx[j/2+ma] = Prx[j/2-ma];
//                         Pix[j/2+ma] = -Pix[j/2-ma];
//                         Pry[j/2+ma] = Pry[j/2-ma];
//                         Piy[j/2+ma] = -Piy[j/2-ma];
//                         Prz[j/2+ma] = Prz[j/2-ma];
//                         Piz[j/2+ma] = -Piz[j/2-ma];                        
//                     } else {
//                         Pr[j/2+ma] = -Pr[j/2-ma];
//                         Pi[j/2+ma] = Pi[j/2-ma];
//                         Prx[j/2+ma] = -Prx[j/2-ma];
//                         Pix[j/2+ma] =  Pix[j/2-ma];
//                         Pry[j/2+ma] = -Pry[j/2-ma];
//                         Piy[j/2+ma] =  Piy[j/2-ma];
//                         Prz[j/2+ma] = -Prz[j/2-ma];
//                         Piz[j/2+ma] =  Piz[j/2-ma];                        
//                     }
//                     mapar = -mapar;        
//                 }                                                        
//             }
//             
//             if (j==(2*mb)) {
//                 int mapar = 1;
//                 for (int ma = 0; ma <= j; ma++) {
//                     if (mapar == 1) {                    
//                         Prx[j/2+ma] = Prx[j/2-ma];
//                         Pix[j/2+ma] = -Pix[j/2-ma];
//                         Pry[j/2+ma] = Pry[j/2-ma];
//                         Piy[j/2+ma] = -Piy[j/2-ma];
//                         Prz[j/2+ma] = Prz[j/2-ma];
//                         Piz[j/2+ma] = -Piz[j/2-ma];                        
//                     } else {
//                         Prx[j/2+ma] = -Prx[j/2-ma];
//                         Pix[j/2+ma] =  Pix[j/2-ma];
//                         Pry[j/2+ma] = -Pry[j/2-ma];
//                         Piy[j/2+ma] =  Piy[j/2-ma];
//                         Prz[j/2+ma] = -Prz[j/2-ma];
//                         Piz[j/2+ma] =  Piz[j/2-ma];                        
//                     }
//                     mapar = -mapar;        
//                 }                                                        
//             }
//                         
//             // store Qr, Qi, for the next mb level
//             if (j==(2*mb+1)) {
//                 int mapar = (((j+1)/2)%2 == 0) ? -1 : 1;
//                 for (int ma = 0; ma <= j; ma++) {     
//                     if (mapar == 1) {                    
//                         Qr[j-ma] = Pr[ma];
//                         Qi[j-ma] = -Pi[ma];  
//                         Qrx[j-ma] =  Prx[ma];
//                         Qix[j-ma] = -Pix[ma];
//                         Qry[j-ma] =  Pry[ma];
//                         Qiy[j-ma] = -Piy[ma];
//                         Qrz[j-ma] =  Prz[ma];
//                         Qiz[j-ma] = -Piz[ma];                                            
//                     } else {
//                         Qr[j-ma] = -Pr[ma];
//                         Qi[j-ma] =  Pi[ma];
//                         Qrx[j-ma] = -Prx[ma];
//                         Qix[j-ma] =  Pix[ma];
//                         Qry[j-ma] = -Pry[ma];
//                         Qiy[j-ma] =  Piy[ma];
//                         Qrz[j-ma] = -Prz[ma];
//                         Qiz[j-ma] =  Piz[ma];                                            
//                     }
//                     mapar = -mapar;
//                 }                                                
//             }
//             
//             int k =  1 + (j+1)*mb;
//             for (int ma = 2; ma <= j; ma++)
//                 k += ma*ma;            
//             for (int ma = 0; ma <= j; ma++) {
//                 int in = i + inum*k;                
//                 Stotr[in] += sfac*Pr[ma]; // atomic add   
//                 Stoti[in] += sfac*Pi[ma]; // atomic add     
//                 dulist_r[ij+ijnum*k+ijnum*idxu_max*0] = dsfac * Pr[ma] * ux + sfac * Prx[ma];
//                 dulist_i[ij+ijnum*k+ijnum*idxu_max*0] = dsfac * Pi[ma] * ux + sfac * Pix[ma];
//                 dulist_r[ij+ijnum*k+ijnum*idxu_max*1] = dsfac * Pr[ma] * uy + sfac * Pry[ma];
//                 dulist_i[ij+ijnum*k+ijnum*idxu_max*1] = dsfac * Pi[ma] * uy + sfac * Piy[ma];
//                 dulist_r[ij+ijnum*k+ijnum*idxu_max*2] = dsfac * Pr[ma] * uz + sfac * Prz[ma];
//                 dulist_i[ij+ijnum*k+ijnum*idxu_max*2] = dsfac * Pi[ma] * uz + sfac * Piz[ma];                            
//                 k += 1; 
//             }                                  
//             
//             k =  1 + (j+1)*mb;
//             for (int ma = 2; ma <= j; ma++)
//                 k += ma*ma;                            
//             if (j==(2*mb)) {
//                 for (int ma = 0; ma < mb; ma++) {
//                     int in = i + inum*k;                
//                     T y_r = ylist_r[in];
//                     T y_i = ylist_i[in];
//                     de[0] += (dsfac * Pr[ma] * ux + sfac * Prx[ma]) * y_r + (dsfac * Pi[ma] * ux + sfac * Pix[ma]) * y_i;
//                     de[1] += (dsfac * Pr[ma] * uy + sfac * Pry[ma]) * y_r + (dsfac * Pi[ma] * uy + sfac * Piy[ma]) * y_i;
//                     de[2] += (dsfac * Pr[ma] * uz + sfac * Prz[ma]) * y_r + (dsfac * Pi[ma] * uz + sfac * Piz[ma]) * y_i;
//                     k += 1; 
//                 }     
//                 ma = mb;
//                 int in = i + inum*k;                
//                 T y_r = 0.5*ylist_r[in];
//                 T y_i = 0.5*ylist_i[in];
//                 de[0] += (dsfac * Pr[ma] * ux + sfac * Prx[ma]) * y_r + (dsfac * Pi[ma] * ux + sfac * Pix[ma]) * y_i;
//                 de[1] += (dsfac * Pr[ma] * uy + sfac * Pry[ma]) * y_r + (dsfac * Pi[ma] * uy + sfac * Piy[ma]) * y_i;
//                 de[2] += (dsfac * Pr[ma] * uz + sfac * Prz[ma]) * y_r + (dsfac * Pi[ma] * uz + sfac * Piz[ma]) * y_i;
//             }
//             else {
//                 for (int ma = 0; ma <= j; ma++) {
//                     int in = i + inum*k;                
//                     T y_r = ylist_r[in];
//                     T y_i = ylist_i[in];
//                     de[0] += (dsfac * Pr[ma] * ux + sfac * Prx[ma]) * y_r + (dsfac * Pi[ma] * ux + sfac * Pix[ma]) * y_i;
//                     de[1] += (dsfac * Pr[ma] * uy + sfac * Pry[ma]) * y_r + (dsfac * Pi[ma] * uy + sfac * Piy[ma]) * y_i;
//                     de[2] += (dsfac * Pr[ma] * uz + sfac * Prz[ma]) * y_r + (dsfac * Pi[ma] * uz + sfac * Piz[ma]) * y_i;
//                     k += 1; 
//                 }                
//             }            
//         }
//     }
//   
// //     if (ij<20) 
// //         printf("%g ", dulist_r[ij+ijnum*0+ijnum*idxu_max*0]);            
// //     if (ij==20) printf("\n");
//     
//     dedr[0 + 3*ij] = 2.0*de[0];      
//     dedr[1 + 3*ij] = 2.0*de[1];      
//     dedr[2 + 3*ij] = 2.0*de[2];         
//   }   
// }
// // template void cpuComputeDeidrj(double*, double*, double*, double*, double*,  double*, double*, 
// //          double, double, double, int*, int*, int*, int*, int*, int*, int, int, int, int, int, int);
// // template void cpuComputeDeidrj(float*, float*, float*, float*, float*, float*, float*, 
// //         float, float, float, int*, int*, int*, int*, int*, int*, int, int, int, int, int, int);
// template void cpuComputeDeidrj(double*, double*, double*, double*, double*, double*, double*,  double*, double*, 
//          double*, double*, double, double, double, int*, int*, int*, int*, int*, int*, int, int, int, int, int, int);
// template void cpuComputeDeidrj(float*, float*, float*, float*, float*, float*, float*, float*, float*, 
//         float*, float*, float, float, float, int*, int*, int*, int*, int*, int*, int, int, int, int, int, int);

template <typename T> void cpuComputeEi(T *ei, T *ylist_r, T *ylist_i, 
         T *Stotr, T *Stoti, int idxu_max, int nelem, int inum)
{    
    int m = idxu_max*nelem;
    for (int i=0; i<inum; i++) {           
        ei[i] = 0.0;
        for (int j=0; j<m; j++) {
            int k = i + inum*j;
            ei[i] += ylist_r[k]*Stotr[k] + ylist_i[k]*Stoti[k]; 
        }
    }
}
template void cpuComputeEi(double*, double*, double*, double*, double*, int, int, int);
template void cpuComputeEi(float*, float*, float*, float*, float*, int, int, int);


// template <typename T> void cpuComputeEi(T *ei, T *ylist_r, T *ylist_i, 
//          T *Stotr, T *Stoti, int idxu_max, int nelem, int inum)
// {    
//     int N1 = inum;
//     int N2 = N1*(twojmax+1);
//     int N3 = N2*nelements;                                
//     
//     for (int idx=0; idx < N3; idx++) {
//         int l = idx%N2;  // inum*(twojmax+1)
//         int ii = l%N1;    // inum
//         int j = (l-ii)/N1; // (twojmax+1)
//         int jelem = (idx-l)/N2; // nelements   
//         int nmax = ii + inum*idxu_max*jelem;
// 
//         
// //         int jju = idxu_block[j];                
// //         for (int mb = 0; mb <= j; mb++) {
// //             for (int ma = 0; ma <= j; ma++) {                
// //                 if (jelem == ielem || wselfall_flag)
// //                     if (ma==mb)                        
// //                         Stotr[inum*jju + nmax] += wself;                                     
// //                 jju++;                                
// //             }
// //         }
// //         
// //         // copy left side to right side with inversion symmetry VMK 4.4(2)
// //         // u[ma-j][mb-j] = (-1)^(ma-mb)*Conj([u[ma][mb])
// //         
// //         jju = idxu_block[j];        
// //         int jjup = jju+(j+1)*(j+1)-1;
// //         int mbpar = 1;
// //         for (int mb = 0; 2*mb < j; mb++) {
// //             int mapar = mbpar;
// //             for (int ma = 0; ma <= j; ma++) {
// //                 int njju =  inum*jju + nmax;
// //                 int njjup = inum*jjup + nmax;
// //                 if (mapar == 1) {
// //                     Stotr[njjup] = Stotr[njju];
// //                     Stoti[njjup] = -Stoti[njju];
// //                 } else {
// //                     Stotr[njjup] = -Stotr[njju];
// //                     Stoti[njjup] =  Stoti[njju];
// //                 }
// //                 mapar = -mapar;
// //                 jju++;
// //                 jjup--;
// //             }
// //             mbpar = -mbpar;
// //         }        
//     }                    
// }
// 
// template <typename T> void cpuComputeEi(T *ei, T *ylist_r, T *ylist_i, 
//         T *Stotr, T *Stoti, int *idxu_block, int *map, int *ai, int *tj,
//         int twojmax, int idxu_max, int chemflag, int inum, int ijnum) 
// {                
//     cpuArraySetValue(dedr, (T) 0.0, ijnum*3);
//     
//     int niimax = inum*idxu_max;
//     int nijmax = ijnum*idxu_max;  
//     for(int ij=0; ij<ijnum; ij++) {                        
//       int jelem = (chemflag) ? map[tj[ij]] : 0; //(chemflag) ? map[type[alist[aj[ij]]]] : 0;
//       int nimax = niimax*jelem;              
//       int i = ai[ij]; // atom i        
//       T de[3];
//       de[0] = 0.0; 
//       de[1] = 0.0; 
//       de[2] = 0.0;
//       for(int j = 0; j <= twojmax; j++) {
//         int jju = idxu_block[j];        
//         
//         for(int mb = 0; 2*mb < j; mb++)
//           for(int ma = 0; ma <= j; ma++) {
//             int n1 = i + inum*jju + nimax;
//             int n2 = ij+ ijnum*jju;
//             T y_r = ylist_r[n1];
//             T y_i = ylist_i[n1];
//             for(int k = 0; k < 3; k++)
//               de[k] += dulist_r[n2+nijmax*k] * y_r + dulist_i[n2+nijmax*k] * y_i;
//             jju++;
//           } //end loop over ma mb
// 
//         // For j even, handle middle column
// 
//         if (j%2 == 0) {
//           int mb = j/2;
//           for(int ma = 0; ma < mb; ma++) {
//             int n1 = i + inum*jju + nimax;
//             int n2 = ij+ ijnum*jju;
//             T y_r = ylist_r[n1];
//             T y_i = ylist_i[n1];
//             for(int k = 0; k < 3; k++)
//                 de[k] += dulist_r[n2+nijmax*k] * y_r + dulist_i[n2+nijmax*k] * y_i;
//             jju++;
//           }
// 
//           int n1 = i + inum*jju + nimax;
//           int n2 = ij+ ijnum*jju;
//           T y_r = ylist_r[n1];
//           T y_i = ylist_i[n1];
//           for(int k = 0; k < 3; k++)
//             de[k] += (dulist_r[n2+nijmax*k] * y_r + dulist_i[n2+nijmax*k] * y_i)*0.5;
//           // jju++;
//         } // end if jeven
//       } // end loop over j
// 
//       for(int k = 0; k < 3; k++)
//         dedr[k + 3*ij] = 2.0*de[k];      
//     }  
// }
// template void cpuComputeEi(double*, double*, double*, double*, double*, 
//         int*, int*, int*, int*, int, int, int, int, int);
// template void cpuComputeEi(float*, float*, float*, float*, float*, 
//         int*, int*, int*, int*, int, int, int, int, int);



// template <typename T> void cpuComputeUi(T *Stotr, T *Stoti, T *rootpqarray, T *rij, T *wjelem, T *radelem, 
//         T rmin0, T rfac0, T rcutfac, int *idxu_block, int *map, int *ai, int *ti, int *tj, 
//         int twojmax, int idxu_max, int inum, int ijnum, int switchflag, int chemflag)                
// {    
//   for(int ij=0; ij<ijnum; ij++) {        
//     T x = rij[ij*3+0];
//     T y = rij[ij*3+1];
//     T z = rij[ij*3+2];
//     T rsq = x * x + y * y + z * z;
//     T r = sqrt(rsq);
// 
//     T rcutij = (radelem[ti[ij]]+radelem[tj[ij]])*rcutfac; //(radelem[type[ii]]+radelem[type[jj]])*rcutfac;
//     T rscale0 = rfac0 * M_PI / (rcutij - rmin0);
//     T theta0 = (r - rmin0) * rscale0;
//     T z0 = r / tan(theta0);                
//             
//     T sfac = 0.0;
//     if (switchflag == 0) 
//         sfac = 1.0;    
//     else if (switchflag == 1) {
//         if (r <= rmin0) {
//             sfac = 1.0;
//         }
//         else if(r > rcutij) {
//             sfac = 1.0;
//         }
//         else {
//             T rcutfac = M_PI / (rcutij - rmin0);
//             sfac =  0.5 * (cos((r - rmin0) * rcutfac) + 1.0);   
//         }
//     } 
//     sfac *= wjelem[tj[ij]];
// 
//     T r0inv;
//     T a_r, a_i, b_r, b_i;
//     T rootpq, rootpq2;
//     int jdim = twojmax + 1;
//   
//     r0inv = 1.0 / sqrt(r * r + z0 * z0);
//     a_r = r0inv * z0;
//     a_i = -r0inv * z;
//     b_r = r0inv * y;
//     b_i = -r0inv * x;
//         
//     int i, j, ni, njelem, ma, mb;        
//     T up_r[20];
//     T up_i[20];
//         
//     njelem = (chemflag==0) ? 0 : (inum*idxu_max*map[tj[ij]]);    
//     i = ai[ij] + njelem;             
//     
//     //j = 0;    
//     Stotr[i] += sfac; // atomic add   
//     
//     j = 1;    
//     mb = 0;
//     ma = 0;
//     rootpq = rootpqarray[(j - ma)*jdim + (j - mb)];        
//     up_r[0] = rootpq * (a_r * sfac + a_i * 0.0);
//     up_i[0] = rootpq * (a_r * 0.0 - a_i * sfac);
//     ni = i + inum*1;          
//     Stotr[ni] += up_r[0]; // atomic add   
//     Stoti[ni] += up_i[0]; // atomic add   
//     
//     rootpq = rootpqarray[(ma + 1)*jdim + (j - mb)];    
//     up_r[1] = -rootpq * (b_r * sfac + b_i * 0.0);
//     up_i[1] = -rootpq * (b_r * 0.0 - b_i * sfac);
//     ni = i + inum*2;     
//     Stotr[ni] += up_r[1]; // atomic add   
//     Stoti[ni] += up_i[1]; // atomic add      
//                         
//     j = 2;                
//     T t1 = a_r * up_r[0] + a_i * up_i[0];
//     T t2 = a_r * up_i[0] - a_i * up_r[0];
//     T t3 = a_r * up_r[1] + a_i * up_i[1];
//     T t4 = a_r * up_i[1] - a_i * up_r[1];        
//     T t5 = b_r * up_r[0] + b_i * up_i[0];
//     T t6 = b_r * up_i[0] - b_i * up_r[0];
//     T t7 = b_r * up_r[1] + b_i * up_i[1];
//     T t8 = b_r * up_i[1] - b_i * up_r[1];
// 
//     mb = 0;        
//     ma = 0;
//     rootpq = rootpqarray[(j - ma)*jdim + (j - mb)];            
//     up_r[0] = rootpq * t1; // (a_r * up_r[0] + a_i * up_i[0])
//     up_i[0] = rootpq * t2; // (a_r * up_i[0] - a_i * up_r[0])       
//     ni = i + inum*5;          
//     Stotr[ni] += up_r[0]; // atomic add   
//     Stoti[ni] += up_i[0]; // atomic add   
// 
//     ma = 1;        
//     rootpq = rootpqarray[(j - ma)*jdim + (j - mb)];               
//     rootpq2 = rootpqarray[ma*jdim + (j - mb)];      
//     // (b_r * up_r[0] + b_i * up_i[0]) - (a_r * up_r[1] + a_i * up_i[1])
//     up_r[1] = rootpq2 * t5 - rootpq * t3;
//     // (b_r * up_i[0] - b_i * up_i[0]) - (a_r * up_r[1] - a_i * up_i[1])
//     up_i[1] = rootpq2 * t6 - rootpq * t4;
//     ni = i + inum*6;       
//     Stotr[ni] += up_r[1]; // atomic add   
//     Stoti[ni] += up_i[1]; // atomic add   
// 
//     ma = 2;                   
//     rootpq2 = rootpqarray[ma*jdim + (j - mb)];    
//     up_r[2] = -rootpq2 * t7; // (b_r * up_r[1] + b_i * up_i[1])
//     up_i[2] = -rootpq2 * t8; // (b_r * up_i[1] - b_i * up_i[1])
//     ni = i + inum*7;       
//     Stotr[ni] += up_r[2]; // atomic add   
//     Stoti[ni] += up_i[2]; // atomic add        
// 
//     mb = 1;
//     ma = 0;
//     rootpq = rootpqarray[(j - ma)*jdim + (j - mb)];            
//     up_r[4] = -rootpq * t3; // -(a_r * up_r[1] + a_i * up_i[1])
//     up_i[4] = -rootpq * t4; // -(a_r * up_i[1] - a_i * up_r[1])        
//     ni = i + inum*8;          
//     Stotr[ni] += up_r[4]; // atomic add   
//     Stoti[ni] += up_i[4]; // atomic add   
// 
//     ma = 1;        
//     rootpq = rootpqarray[(j - ma)*jdim + (j - mb)];               
//     rootpq2 = rootpqarray[ma*jdim + (j - mb)];      
//     // -(b_r * up_r[1] + b_i * up_i[1]) + (a_r * up_r[0] + a_i * up_i[0])
//     up_r[5] = rootpq * t1 - rootpq2 * t7;
//     // -(b_r * up_i[1] - b_i * up_i[1]) + (a_r * up_r[0] - a_i * up_i[0])
//     up_i[5] = rootpq * t2 - rootpq2 * t8;
//     ni = i + inum*9;       
//     Stotr[ni] += up_r[5]; // atomic add   
//     Stoti[ni] += up_i[5]; // atomic add   
// 
//     ma = 2;                   
//     rootpq2 = rootpqarray[ma*jdim + (j - mb)];    
//     up_r[6] = rootpq2 * t5; //-rootpq2 * (b_r * u2_r + b_i * u2_i);
//     up_i[6] = rootpq2 * t6; //-rootpq2 * (b_r * u2_i - b_i * u2_r);
//     ni = i + inum*10;       
//     Stotr[ni] += up_r[6]; // atomic add   
//     Stoti[ni] += up_i[6]; // atomic add        
//         
//     if (twojmax >= 4) {
//         j = 3;
//     
//         T t9, t10, t11, t12;
//         
//         mb = 0;    
//         t1  = a_r * up_r[0] + a_i * up_i[0];
//         t2  = a_r * up_i[0] - a_i * up_r[0];
//         t3  = a_r * up_r[1] + a_i * up_i[1];
//         t4  = a_r * up_i[1] - a_i * up_r[1];  
//         t5  = a_r * up_r[2] + a_i * up_i[2];
//         t6  = a_r * up_i[2] - a_i * up_r[2];        
//         t7  = b_r * up_r[0] + b_i * up_i[0];
//         t8  = b_r * up_i[0] - b_i * up_r[0];
//         t9  = b_r * up_r[1] + b_i * up_i[1];
//         t10 = b_r * up_i[1] - b_i * up_r[1];
//         t11 = b_r * up_r[2] + b_i * up_i[2];
//         t12 = b_r * up_i[2] - b_i * up_r[2];
//                     
//         ma = 0;
//         rootpq = rootpqarray[(j - ma)*jdim + (j - mb)];            
//         up_r[0] = rootpq * t1; // (a_r * up_r[0] + a_i * up_i[0])
//         up_i[0] = rootpq * t2; // (a_r * up_i[0] - a_i * up_r[0])       
//         ni = i + inum*14;          
//         Stotr[ni] += up_r[0]; // atomic add   
//         Stoti[ni] += up_i[0]; // atomic add   
//         
//         ma = 1;        
//         rootpq = rootpqarray[(j - ma)*jdim + (j - mb)];               
//         rootpq2 = rootpqarray[ma*jdim + (j - mb)];      
//         up_r[1] = rootpq2 * t7 - rootpq * t3;        
//         up_i[1] = rootpq2 * t8 - rootpq * t4;
//         ni = i + inum*15;       
//         Stotr[ni] += up_r[1]; // atomic add   
//         Stoti[ni] += up_i[1]; // atomic add   
//                 
//         ma = 2;        
//         rootpq = rootpqarray[(j - ma)*jdim + (j - mb)];               
//         rootpq2 = rootpqarray[ma*jdim + (j - mb)];      
//         up_r[2] = rootpq2 * t9 - rootpq * t5;        
//         up_i[2] = rootpq2 * t10 - rootpq * t6;
//         ni = i + inum*16;       
//         Stotr[ni] += up_r[2]; // atomic add   
//         Stoti[ni] += up_i[2]; // atomic add   
//         
//         ma = 3;                   
//         rootpq2 = rootpqarray[ma*jdim + (j - mb)];    
//         up_r[3] = -rootpq2 * t11; //-rootpq2 * (b_r * u2_r + b_i * u2_i);
//         up_i[3] = -rootpq2 * t12; //-rootpq2 * (b_r * u2_i - b_i * u2_r);
//         ni = i + inum*17;       
//         Stotr[ni] += up_r[3]; // atomic add   
//         Stoti[ni] += up_i[3]; // atomic add        
//         
//         mb = 1;        
//         t1  = a_r * up_r[6] + a_i * up_i[6];
//         t2  = a_r * up_i[6] - a_i * up_r[6];
//         t3  = a_r * up_r[5] + a_i * up_i[5];
//         t4  = a_r * up_i[5] - a_i * up_r[5];  
//         t5  = a_r * up_r[4] + a_i * up_i[4];
//         t6  = a_r * up_i[4] - a_i * up_r[4];        
//         t7  = b_r * up_r[6] + b_i * up_i[6];
//         t8  = b_r * up_i[6] - b_i * up_r[6];
//         t9  = b_r * up_r[5] + b_i * up_i[5];
//         t10 = b_r * up_i[5] - b_i * up_r[5];
//         t11 = b_r * up_r[4] + b_i * up_i[4];
//         t12 = b_r * up_i[4] - b_i * up_r[4];
//                 
//         ma = 0;
//         rootpq = rootpqarray[(j - ma)*jdim + (j - mb)];            
//         up_r[4] = rootpq * t1; // (a_r * up_r[0] + a_i * up_i[0])
//         up_i[4] = rootpq * t2; // (a_r * up_i[0] - a_i * up_r[0])       
//         ni = i + inum*18;          
//         Stotr[ni] += up_r[4]; // atomic add   
//         Stoti[ni] += up_i[4]; // atomic add   
//         
//         ma = 1;        
//         rootpq = rootpqarray[(j - ma)*jdim + (j - mb)];               
//         rootpq2 = rootpqarray[ma*jdim + (j - mb)];      
//         up_r[5] = rootpq2 * t7 - rootpq * t3;        
//         up_i[5] = rootpq2 * t8 - rootpq * t4;
//         ni = i + inum*19;       
//         Stotr[ni] += up_r[5]; // atomic add   
//         Stoti[ni] += up_i[5]; // atomic add   
//                 
//         ma = 2;        
//         rootpq = rootpqarray[(j - ma)*jdim + (j - mb)];               
//         rootpq2 = rootpqarray[ma*jdim + (j - mb)];      
//         up_r[6] = rootpq2 * t9 - rootpq * t5;        
//         up_i[6] = rootpq2 * t10 - rootpq * t6;
//         ni = i + inum*20;       
//         Stotr[ni] += up_r[6]; // atomic add   
//         Stoti[ni] += up_i[6]; // atomic add   
//         
//         ma = 3;                   
//         rootpq2 = rootpqarray[ma*jdim + (j - mb)];    
//         up_r[7] = -rootpq2 * t11; //-rootpq2 * (b_r * u2_r + b_i * u2_i);
//         up_i[7] = -rootpq2 * t12; //-rootpq2 * (b_r * u2_i - b_i * u2_r);
//         ni = i + inum*21;       
//         Stotr[ni] += up_r[7]; // atomic add   
//         Stoti[ni] += up_i[7]; // atomic add        
//         
//     }
//     
// //     T Sr[300];
// //     T Si[300];
// //         
// //     Sr[0] = 1.0;
// //     Si[0] = 0.0;
// //     
// //     int njelem = (chemflag==0) ? 0 : (inum*idxu_max*map[tj[ij]]);    
// //     int i = ai[ij] + njelem;                
// //     Stotr[i] += sfac; // atomic add   
// //     
// //     int k = 1;
// //     for (int j = 1; j <= twojmax; j++) {
// //         int jju = idxu_block[j];
// //         int jjup = idxu_block[j-1];
// //         
// //         // fill in left side of matrix layer from previous layer
// //         for (int mb = 0; 2*mb <= j; mb++) {
// //             Sr[jju] = 0.0;
// //             Si[jju] = 0.0;
// //             for (int ma = 0; ma < j; ma++) {
// //                 rootpq = rootpqarray[(j - ma)*jdim + (j - mb)];
// //                 int njju = jju;
// //                 int njju1 = (jju+1);
// //                 int njjup = jjup;
// //                 T u_r = Sr[njjup];
// //                 T u_i = Si[njjup];
// // 
// //                 Sr[njju] += rootpq * (a_r * u_r + a_i * u_i);
// //                 Si[njju] += rootpq * (a_r * u_i - a_i * u_r);
// // 
// //                 rootpq = rootpqarray[(ma + 1)*jdim + (j - mb)];
// //                 Sr[njju1] = -rootpq * (b_r * u_r + b_i * u_i);
// //                 Si[njju1] = -rootpq * (b_r * u_i - b_i * u_r);
// //                 jju++;
// //                 jjup++;
// //             }
// //             jju++;
// //         }
// // 
// //     // copy left side to right side with inversion symmetry VMK 4.4(2)
// //     // u[ma-j][mb-j] = (-1)^(ma-mb)*Conj([u[ma][mb])
// //         
// //         jju = idxu_block[j];
// //         jjup = jju+(j+1)*(j+1)-1;
// //         int mbpar = 1;
// //         for (int mb = 0; 2*mb <= j; mb++) {
// //             int mapar = mbpar;
// //             for (int ma = 0; ma <= j; ma++) {
// //                 int njju = jju;
// //                 int njjup = jjup;
// //                 if (mapar == 1) {
// //                     Sr[njjup] = Sr[njju];
// //                     Si[njjup] = -Si[njju];
// //                 } else {
// //                     Sr[njjup] = -Sr[njju];
// //                     Si[njjup] =  Si[njju];
// //                 }
// //                 mapar = -mapar;
// //                 jju++;
// //                 jjup--;
// //             }
// //             mbpar = -mbpar;
// //         }
// //         
// //         for (int mb = 0; mb <= j; mb++) 
// //             for (int ma = 0; ma <= j; ma++) {
// //                 i = ai[ij] + inum*k + njelem;                
// //                 Stotr[i] += sfac*Sr[k];
// //                 Stoti[i] += sfac*Si[k];
// //                 k += 1; 
// //             }
// //     }        
//        
// //     for (int k=0; k<idxu_max; k++) {                             
// //         int i = ai[ij] + inum*k + njelem;                
// //         Stotr[i] += sfac*Sr[k];
// //         Stoti[i] += sfac*Si[k];
// //     }                
// //     if (chemflag==0) {
// //         for (int k=0; k<idxu_max; k++) {        
// //             int i = ai[ij] + inum*k;  
// //             Stotr[i] += sfac*Sr[k];
// //             Stoti[i] += sfac*Si[k];
// //         }    
// //     }
// //     else {
// //         int jelem = map[tj[ij]];  
// //         int njelem = inum*idxu_max*jelem;
// //         for (int k=0; k<idxu_max; k++) {                                   
// //             int i = ai[ij] + inum*k + njelem;                
// //             Stotr[i] += sfac*Sr[k];
// //             Stoti[i] += sfac*Si[k];
// //         }            
// //     }        
//   }
// };
// template void cpuComputeUi(double*, double*, double*, double*, double*, double*, double, double, double, 
//         int*, int*, int*, int*, int*, int, int, int, int, int, int);
// template void cpuComputeUi(float*, float*, float*, float*, float*, float*, float, float, float, 
//         int*, int*, int*, int*, int*, int, int, int, int, int, int);


// template <typename T> void cpuComputeUi(T *Stotr, T *Stoti, T *rootpqarray, T *rij, T *wjelem, T *radelem, 
//         T rmin0, T rfac0, T rcutfac, int *idxu_block, int *map, int *ai, int *ti, int *tj, 
//         int twojmax, int idxu_max, int inum, int ijnum, int switchflag, int chemflag)                
// {    
//   for(int ij=0; ij<ijnum; ij++) {        
//     T x = rij[ij*3+0];
//     T y = rij[ij*3+1];
//     T z = rij[ij*3+2];
//     T rsq = x * x + y * y + z * z;
//     T r = sqrt(rsq);
// 
//     T rcutij = (radelem[ti[ij]]+radelem[tj[ij]])*rcutfac; //(radelem[type[ii]]+radelem[type[jj]])*rcutfac;
//     T rscale0 = rfac0 * M_PI / (rcutij - rmin0);
//     T theta0 = (r - rmin0) * rscale0;
//     T z0 = r / tan(theta0);                
//             
//     T sfac = 0.0;
//     if (switchflag == 0) 
//         sfac = 1.0;    
//     else if (switchflag == 1) {
//         if (r <= rmin0) {
//             sfac = 1.0;
//         }
//         else if(r > rcutij) {
//             sfac = 1.0;
//         }
//         else {
//             T rcutfac = M_PI / (rcutij - rmin0);
//             sfac =  0.5 * (cos((r - rmin0) * rcutfac) + 1.0);   
//         }
//     } 
//     sfac *= wjelem[tj[ij]];
// 
//     T r0inv;
//     T a_r, a_i, b_r, b_i;
//     T rootpq;
//     int jdim = twojmax + 1;
//   
//     r0inv = 1.0 / sqrt(r * r + z0 * z0);
//     a_r = r0inv * z0;
//     a_i = -r0inv * z;
//     b_r = r0inv * y;
//     b_i = -r0inv * x;
//         
//     T Sr[300];
//     T Si[300];
//         
//     Sr[0] = 1.0;
//     Si[0] = 0.0;
//     
//     int njelem = (chemflag==0) ? 0 : (inum*idxu_max*map[tj[ij]]);    
//     int i = ai[ij] + njelem;                
//     Stotr[i] += sfac; // atomic add   
//     
//     int k = 1;
//     for (int j = 1; j <= twojmax; j++) {
//         int jju = idxu_block[j];
//         int jjup = idxu_block[j-1];
//         
//         // fill in left side of matrix layer from previous layer
//         for (int mb = 0; 2*mb <= j; mb++) {
//             Sr[jju] = 0.0;
//             Si[jju] = 0.0;
//             for (int ma = 0; ma < j; ma++) {
//                 rootpq = rootpqarray[(j - ma)*jdim + (j - mb)];
//                 int njju = jju;
//                 int njju1 = (jju+1);
//                 int njjup = jjup;
//                 T u_r = Sr[njjup];
//                 T u_i = Si[njjup];
// 
//                 Sr[njju] += rootpq * (a_r * u_r + a_i * u_i);
//                 Si[njju] += rootpq * (a_r * u_i - a_i * u_r);
// 
//                 rootpq = rootpqarray[(ma + 1)*jdim + (j - mb)];
//                 Sr[njju1] = -rootpq * (b_r * u_r + b_i * u_i);
//                 Si[njju1] = -rootpq * (b_r * u_i - b_i * u_r);
//                 jju++;
//                 jjup++;
//             }
//             jju++;
//         }
// 
//     // copy left side to right side with inversion symmetry VMK 4.4(2)
//     // u[ma-j][mb-j] = (-1)^(ma-mb)*Conj([u[ma][mb])
//         
//         jju = idxu_block[j];
//         jjup = jju+(j+1)*(j+1)-1;
//         int mbpar = 1;
//         for (int mb = 0; 2*mb <= j; mb++) {
//             int mapar = mbpar;
//             for (int ma = 0; ma <= j; ma++) {
//                 int njju = jju;
//                 int njjup = jjup;
//                 if (mapar == 1) {
//                     Sr[njjup] = Sr[njju];
//                     Si[njjup] = -Si[njju];
//                 } else {
//                     Sr[njjup] = -Sr[njju];
//                     Si[njjup] =  Si[njju];
//                 }
//                 mapar = -mapar;
//                 jju++;
//                 jjup--;
//             }
//             mbpar = -mbpar;
//         }
//         
//         for (int mb = 0; mb <= j; mb++) 
//             for (int ma = 0; ma <= j; ma++) {
//                 i = ai[ij] + inum*k + njelem;                
//                 Stotr[i] += sfac*Sr[k];
//                 Stoti[i] += sfac*Si[k];
//                 k += 1; 
//             }
//     }        
//        
// //     for (int k=0; k<idxu_max; k++) {                             
// //         int i = ai[ij] + inum*k + njelem;                
// //         Stotr[i] += sfac*Sr[k];
// //         Stoti[i] += sfac*Si[k];
// //     }                
// //     if (chemflag==0) {
// //         for (int k=0; k<idxu_max; k++) {        
// //             int i = ai[ij] + inum*k;  
// //             Stotr[i] += sfac*Sr[k];
// //             Stoti[i] += sfac*Si[k];
// //         }    
// //     }
// //     else {
// //         int jelem = map[tj[ij]];  
// //         int njelem = inum*idxu_max*jelem;
// //         for (int k=0; k<idxu_max; k++) {                                   
// //             int i = ai[ij] + inum*k + njelem;                
// //             Stotr[i] += sfac*Sr[k];
// //             Stoti[i] += sfac*Si[k];
// //         }            
// //     }        
//   }
// };
// template void cpuComputeUi(double*, double*, double*, double*, double*, double*, double, double, double, 
//         int*, int*, int*, int*, int*, int, int, int, int, int, int);
// template void cpuComputeUi(float*, float*, float*, float*, float*, float*, float, float, float, 
//         int*, int*, int*, int*, int*, int, int, int, int, int, int);

// template <typename T> void cpuComputeUi(T *Stotr, T *Stoti, T *rootpqarray, T *rij, T *wjelem, T *radelem, 
//         T rmin0, T rfac0, T rcutfac, int *idxu_block, int *map, int *ai, int *ti, int *tj, 
//         int twojmax, int idxu_max, int inum, int ijnum, int switchflag, int chemflag)                
// {    
//   for(int ij=0; ij<ijnum; ij++) {        
//     T x = rij[ij*3+0];
//     T y = rij[ij*3+1];
//     T z = rij[ij*3+2];
//     T rsq = x * x + y * y + z * z;
//     T r = sqrt(rsq);
// 
//     T rcutij = (radelem[ti[ij]]+radelem[tj[ij]])*rcutfac; //(radelem[type[ii]]+radelem[type[jj]])*rcutfac;
//     T rscale0 = rfac0 * M_PI / (rcutij - rmin0);
//     T theta0 = (r - rmin0) * rscale0;
//     T z0 = r / tan(theta0);                
//             
//     T sfac = 0.0;
//     if (switchflag == 0) 
//         sfac = 1.0;    
//     else if (switchflag == 1) {
//         if (r <= rmin0) {
//             sfac = 1.0;
//         }
//         else if(r > rcutij) {
//             sfac = 1.0;
//         }
//         else {
//             T rcutfac = M_PI / (rcutij - rmin0);
//             sfac =  0.5 * (cos((r - rmin0) * rcutfac) + 1.0);   
//         }
//     } 
//     sfac *= wjelem[tj[ij]];
// 
//     T r0inv;
//     T a_r, a_i, b_r, b_i;
//     T rootpq;
// 
//     r0inv = 1.0 / sqrt(r * r + z0 * z0);
//     a_r = r0inv * z0;
//     a_i = -r0inv * z;
//     b_r = r0inv * y;
//     b_i = -r0inv * x;
//         
//     T Pr[100];
//     T Pi[100];
//     T Sr[120];
//     T Si[120];        
//     Pr[0] = 1.0;
//     Pi[0] = 0.0;    
//     
//     int jdim = twojmax + 1;
//     int njelem = (chemflag==0) ? 0 : (inum*idxu_max*map[tj[ij]]);    
//     int i = ai[ij] + njelem;                        
//     Stotr[i] += sfac; // atomic add   
//     
//     int k = 1;
//     for (int j = 1; j <= twojmax; j++) {
//         int jju = 0;
//         int jjup = 0;
//         
//         // fill in left side of matrix layer from previous layer
//         for (int mb = 0; 2*mb <= j; mb++) {
//             int ma = 0;
//             T p_r, p_i, rootpq2;
//             T u_r = Pr[ma+jdim*mb];
//             T u_i = Pi[ma+jdim*mb];
//             rootpq = rootpqarray[(j - ma)*jdim + (j - mb)];
//             Pr[ma+jdim*mb] = rootpq * (a_r * u_r + a_i * u_i);
//             Pi[ma+jdim*mb] = rootpq * (a_r * u_i - a_i * u_r);            
//             for (ma = 1; ma < j; ma++) {
//                 p_r = Pr[ma+jdim*mb];
//                 p_i = Pi[ma+jdim*mb];
//                 rootpq = rootpqarray[(j - ma)*jdim + (j - mb)];
//                 rootpq2 = rootpqarray[ma*jdim + (j - mb)];
//                 Pr[ma+jdim*mb] = rootpq * (a_r * p_r + a_i * p_i) -rootpq2 * (b_r * u_r + b_i * u_i);
//                 Pi[ma+jdim*mb] = rootpq * (a_r * p_i - a_i * p_r) -rootpq2 * (b_r * u_i - b_i * u_r);
//                 u_r = p_r;
//                 u_i = p_i;
//             }
//             ma = j;
//             rootpq2 = rootpqarray[ma*jdim + (j - mb)];
//             Pr[ma+jdim*mb] = -rootpq2 * (b_r * u_r + b_i * u_i);
//             Pi[ma+jdim*mb] = -rootpq2 * (b_r * u_i - b_i * u_r);
//             
//             
// //             Sr[jju] = 0.0;
// //             Si[jju] = 0.0;
// //             for (int ma = 0; ma < j; ma++) {
// //                 rootpq = rootpqarray[(j - ma)*jdim + (j - mb)];
// //                 int jju1 = (jju+1);
// //                 T u_r = Pr[ma+jdim*mb];
// //                 T u_i = Pi[ma+jdim*mb];
// // 
// //                 Sr[jju] += rootpq * (a_r * u_r + a_i * u_i);
// //                 Si[jju] += rootpq * (a_r * u_i - a_i * u_r);
// // 
// //                 rootpq = rootpqarray[(ma + 1)*jdim + (j - mb)];
// //                 Sr[jju1] = -rootpq * (b_r * u_r + b_i * u_i);
// //                 Si[jju1] = -rootpq * (b_r * u_i - b_i * u_r);
// //                 
// //                 Pr[ma+jdim*mb] = Sr[jju];
// //                 Pi[ma+jdim*mb] = Si[jju];
// //                 
// //                 jju++;
// //                 jjup++;                                
// //             }
// //             Pr[j+jdim*mb] = Sr[jju];
// //             Pi[j+jdim*mb] = Si[jju];            
// //             jju++;
//         }
//             
//         if (j%2 == 0) {
//             jju = (j+1)*j/2 + j/2;
//             jjup = jju;
//             int mapar = 1;
//             for (int ma = 0; ma <= j/2; ma++) {
//                 if (mapar == 1) {                    
//                     //Sr[jju] = Sr[jjup];
//                     //Si[jju] = -Si[jjup];
//                     Pr[j/2+ma+jdim*j/2] = Pr[j/2-ma+jdim*j/2];
//                     Pi[j/2+ma+jdim*j/2] = -Pi[j/2-ma+jdim*j/2];
//                 } else {
//                     //Sr[jju] = -Sr[jjup];
//                     //Si[jju] =  Si[jjup];
//                     Pr[j/2+ma+jdim*j/2] = -Pr[j/2-ma+jdim*j/2];
//                     Pi[j/2+ma+jdim*j/2] = Pi[j/2-ma+jdim*j/2];
//                 }
//                 mapar = -mapar;        
//                 jju++;
//                 jjup--;
//             }                
//         }
//         else {
//             jju = (j+1)*(j-1)/2;
//             jjup = jju+2*(j+1)-1;
//             int mapar = (((j+1)/2)%2 == 0) ? -1 : 1;
//             for (int ma = 0; ma <= j; ma++) {     
//                 if (mapar == 1) {                    
//                     //Pr[j-ma+jdim*(j+1)/2] = Sr[jju];
//                     //Pi[j-ma+jdim*(j+1)/2] = -Si[jju];
//                     Pr[j-ma+jdim*(j+1)/2] = Pr[ma+jdim*(j-1)/2];
//                     Pi[j-ma+jdim*(j+1)/2] = -Pi[ma+jdim*(j-1)/2];
//                 } else {
//                     //Pr[j-ma+jdim*(j+1)/2] = -Sr[jju];
//                     //Pi[j-ma+jdim*(j+1)/2] =  Si[jju];
//                     Pr[j-ma+jdim*(j+1)/2] = -Pr[ma+jdim*(j-1)/2];
//                     Pi[j-ma+jdim*(j+1)/2] =  Pi[ma+jdim*(j-1)/2];
//                 }
//                 mapar = -mapar;
//                 jju++;
//                 jjup--;
//             }            
//         }
//         
//         jju = 0;
//         for (int mb = 0; mb <= j; mb++) 
//             for (int ma = 0; ma <= j; ma++) {
//                 if (2*mb <= j) { // j = 3 -> mb = 0, 1
//                     //printf("%i %i %i %i %i\n", j, ma, mb, jju, k);
//                     int in = i + inum*k;                
//                     //Stotr[in] += sfac*Sr[jju]; // atomic add   
//                     //Stoti[in] += sfac*Si[jju]; // atomic add                       
//                     //Pr[ma+jdim*mb] = Sr[jju];
//                     //Pi[ma+jdim*mb] = Si[jju];                    
//                     Stotr[in] += sfac*Pr[ma+jdim*mb]; // atomic add   
//                     Stoti[in] += sfac*Pi[ma+jdim*mb]; // atomic add                       
//                     jju += 1;
//                 }
//                 k += 1; 
//             }
//     
// //     int k = 1;
// //     for (int j = 1; j <= twojmax; j++) {
// //         int jju = 0;
// //         int jjup = 0;
// //         
// //         // fill in left side of matrix layer from previous layer
// //         for (int mb = 0; 2*mb <= j; mb++) {
// //             Sr[jju] = 0.0;
// //             Si[jju] = 0.0;
// //             for (int ma = 0; ma < j; ma++) {
// //                 rootpq = rootpqarray[(j - ma)*jdim + (j - mb)];
// //                 int jju1 = (jju+1);
// //                 T u_r = Pr[jjup];
// //                 T u_i = Pi[jjup];
// // 
// //                 Sr[jju] += rootpq * (a_r * u_r + a_i * u_i);
// //                 Si[jju] += rootpq * (a_r * u_i - a_i * u_r);
// // 
// //                 rootpq = rootpqarray[(ma + 1)*jdim + (j - mb)];
// //                 Sr[jju1] = -rootpq * (b_r * u_r + b_i * u_i);
// //                 Si[jju1] = -rootpq * (b_r * u_i - b_i * u_r);
// //                 jju++;
// //                 jjup++;                                
// //             }
// //             jju++;
// //         }
// //             
// //         if (j%2 == 0) {
// //             jju = (j+1)*j/2 + j/2;
// //             jjup = jju;
// //             int mapar = 1;
// //             for (int ma = 0; ma <= j/2; ma++) {
// //                 if (mapar == 1) {                    
// //                     Sr[jjup] = Sr[jju];
// //                     Si[jjup] = -Si[jju];
// //                 } else {
// //                     Sr[jjup] = -Sr[jju];
// //                     Si[jjup] =  Si[jju];
// //                 }
// //                 mapar = -mapar;        
// //                 jju++;
// //                 jjup--;
// //             }                
// //         }
// //         else {
// //             jju = (j+1)*(j-1)/2;
// //             jjup = jju+2*(j+1)-1;
// //             int mapar = (((j+1)/2)%2 == 0) ? -1 : 1;
// //             for (int ma = 0; ma <= j; ma++) {     
// //                 if (mapar == 1) {                    
// //                     Pr[jjup] = Sr[jju];
// //                     Pi[jjup] = -Si[jju];
// //                 } else {
// //                     Pr[jjup] = -Sr[jju];
// //                     Pi[jjup] =  Si[jju];
// //                 }
// //                 mapar = -mapar;
// //                 jju++;
// //                 jjup--;
// //             }            
// //         }
// //         
// //         jju = 0;
// //         for (int mb = 0; mb <= j; mb++) 
// //             for (int ma = 0; ma <= j; ma++) {
// //                 if (2*mb <= j) { // j = 3 -> mb = 0, 1
// //                     //printf("%i %i %i %i %i\n", j, ma, mb, jju, k);
// //                     int in = i + inum*k;                
// //                     Stotr[in] += sfac*Sr[jju]; // atomic add   
// //                     Stoti[in] += sfac*Si[jju]; // atomic add                       
// //                     Pr[jju] = Sr[jju];
// //                     Pi[jju] = Si[jju];                    
// //                     jju += 1;
// //                 }
// //                 k += 1; 
// //             }
//         
//     // copy left side to right side with inversion symmetry VMK 4.4(2)
//     // u[ma-j][mb-j] = (-1)^(ma-mb)*Conj([u[ma][mb])
//         
// //         if (j%2 == 0) {
// //             jju = (j+1)*j/2 + j/2;
// //             jjup = jju;
// //             int mapar = 1;
// //             for (int ma = 0; ma <= j/2; ma++) {
// //                 if (mapar == 1) {                    
// //                     Sr[jjup] = Sr[jju];
// //                     Si[jjup] = -Si[jju];
// //                 } else {
// //                     Sr[jjup] = -Sr[jju];
// //                     Si[jjup] =  Si[jju];
// //                 }
// //                 mapar = -mapar;        
// //                 jju++;
// //                 jjup--;
// //             }                
// //         }
// //         
// //         jju = 0;
// //         jjup = jju+(j+1)*(j+1)-1;
// //         int mbpar = 1;
// //         for (int mb = 0; 2*mb < j; mb++) {
// //             int mapar = mbpar;
// //             for (int ma = 0; ma <= j; ma++) {                
// //                 //printf("%i %i %i %i %i %i %i\n", j, ma, mb, mapar, mbpar, jju, jjup);
// //                 if (mapar == 1) {                    
// //                     Sr[jjup] = Sr[jju];
// //                     Si[jjup] = -Si[jju];
// //                 } else {
// //                     Sr[jjup] = -Sr[jju];
// //                     Si[jjup] =  Si[jju];
// //                 }
// //                 mapar = -mapar;
// //                 jju++;
// //                 jjup--;
// //             }
// //             mbpar = -mbpar;
// //         }
// //                 
// //         jju = 0;
// //         for (int mb = 0; mb <= j; mb++) 
// //             for (int ma = 0; ma <= j; ma++) {
// //                 i = ai[ij] + inum*k + njelem;                
// //                 Stotr[i] += sfac*Sr[jju]; // atomic add   
// //                 Stoti[i] += sfac*Si[jju]; // atomic add   
// //                 Pr[jju] = Sr[jju];
// //                 Pi[jju] = Si[jju];
// //                 jju += 1;
// //                 k += 1; 
// //             }
//         
//     }               
//   }
// };
// template void cpuComputeUi(double*, double*, double*, double*, double*, double*, double, double, double, 
//         int*, int*, int*, int*, int*, int, int, int, int, int, int);
// template void cpuComputeUi(float*, float*, float*, float*, float*, float*, float, float, float, 
//         int*, int*, int*, int*, int*, int, int, int, int, int, int);


// template <typename T> void cpuComputeUi(T *Stotr, T *Stoti, T *rootpqarray, T *rij, T *wjelem, T *radelem, 
//         T rmin0, T rfac0, T rcutfac, int *idxu_block, int *map, int *ai, int *ti, int *tj, 
//         int twojmax, int idxu_max, int inum, int ijnum, int switchflag, int chemflag)                
// {    
//   for(int ij=0; ij<ijnum; ij++) {        
//     T x = rij[ij*3+0];
//     T y = rij[ij*3+1];
//     T z = rij[ij*3+2];
//     T rsq = x * x + y * y + z * z;
//     T r = sqrt(rsq);
// 
//     T rcutij = (radelem[ti[ij]]+radelem[tj[ij]])*rcutfac; //(radelem[type[ii]]+radelem[type[jj]])*rcutfac;
//     T rscale0 = rfac0 * M_PI / (rcutij - rmin0);
//     T theta0 = (r - rmin0) * rscale0;
//     T z0 = r / tan(theta0);                
//             
//     T sfac = 0.0;
//     if (switchflag == 0) 
//         sfac = 1.0;    
//     else if (switchflag == 1) {
//         if (r <= rmin0) {
//             sfac = 1.0;
//         }
//         else if(r > rcutij) {
//             sfac = 1.0;
//         }
//         else {
//             T rcutfac = M_PI / (rcutij - rmin0);
//             sfac =  0.5 * (cos((r - rmin0) * rcutfac) + 1.0);   
//         }
//     } 
//     sfac *= wjelem[tj[ij]];
// 
//     T r0inv;
//     T a_r, a_i, b_r, b_i;
//     T rootpq;
// 
//     r0inv = 1.0 / sqrt(r * r + z0 * z0);
//     a_r = r0inv * z0;
//     a_i = -r0inv * z;
//     b_r = r0inv * y;
//     b_i = -r0inv * x;
//         
//     T Pr[300];
//     T Pi[300];
//     Pr[0] = 1.0;
//     Pi[0] = 0.0;    
//     
//     int jdim = twojmax + 1;
//     int njelem = (chemflag==0) ? 0 : (inum*idxu_max*map[tj[ij]]);    
//     int i = ai[ij] + njelem;                        
//     Stotr[i] += sfac; // atomic add   
//     
//     int k = 1;
//     for (int j = 1; j <= twojmax; j++) {        
//         // fill in left side of matrix layer from previous layer
//         for (int mb = 0; 2*mb <= j; mb++) {
//             int ma = 0;
//             T p_r, p_i, rootpq2;
//             T u_r = Pr[ma+jdim*mb];
//             T u_i = Pi[ma+jdim*mb];
//             rootpq = rootpqarray[(j - ma)*jdim + (j - mb)];
//             Pr[ma+jdim*mb] = rootpq * (a_r * u_r + a_i * u_i);
//             Pi[ma+jdim*mb] = rootpq * (a_r * u_i - a_i * u_r);            
//             for (ma = 1; ma < j; ma++) {
//                 p_r = Pr[ma+jdim*mb];
//                 p_i = Pi[ma+jdim*mb];
//                 rootpq = rootpqarray[(j - ma)*jdim + (j - mb)];
//                 rootpq2 = rootpqarray[ma*jdim + (j - mb)];
//                 Pr[ma+jdim*mb] = rootpq * (a_r * p_r + a_i * p_i) -rootpq2 * (b_r * u_r + b_i * u_i);
//                 Pi[ma+jdim*mb] = rootpq * (a_r * p_i - a_i * p_r) -rootpq2 * (b_r * u_i - b_i * u_r);
//                 u_r = p_r;
//                 u_i = p_i;
//             }
//             ma = j;
//             rootpq2 = rootpqarray[ma*jdim + (j - mb)];
//             Pr[ma+jdim*mb] = -rootpq2 * (b_r * u_r + b_i * u_i);
//             Pi[ma+jdim*mb] = -rootpq2 * (b_r * u_i - b_i * u_r);                        
//         }
//             
//         if (j%2 == 0) { // handle middle column when j is even
//             int mapar = 1;
//             for (int ma = 0; ma <= j/2; ma++) {
//                 if (mapar == 1) {                    
//                     Pr[j/2+ma+jdim*j/2] = Pr[j/2-ma+jdim*j/2];
//                     Pi[j/2+ma+jdim*j/2] = -Pi[j/2-ma+jdim*j/2];
//                 } else {
//                     Pr[j/2+ma+jdim*j/2] = -Pr[j/2-ma+jdim*j/2];
//                     Pi[j/2+ma+jdim*j/2] = Pi[j/2-ma+jdim*j/2];
//                 }
//                 mapar = -mapar;        
//             }                                        
//         }
//         else { // use symmetry to copy an extra column when j is odd
//             int mapar = (((j+1)/2)%2 == 0) ? -1 : 1;
//             for (int ma = 0; ma <= j; ma++) {     
//                 if (mapar == 1) {                    
//                     Pr[j-ma+jdim*(j+1)/2] = Pr[ma+jdim*(j-1)/2];
//                     Pi[j-ma+jdim*(j+1)/2] = -Pi[ma+jdim*(j-1)/2];
//                 } else {
//                     Pr[j-ma+jdim*(j+1)/2] = -Pr[ma+jdim*(j-1)/2];
//                     Pi[j-ma+jdim*(j+1)/2] =  Pi[ma+jdim*(j-1)/2];
//                 }
//                 mapar = -mapar;
//             }                                    
//         }
//         
//         int mbpar = 1;
//         for (int mb = 0; 2*mb < j; mb++) {
//             int mapar = mbpar;
//             for (int ma = 0; ma <= j; ma++) {                
//                 if (mapar == 1) {                    
//                     Pr[j-ma+jdim*(j-mb)] = Pr[ma+jdim*mb];
//                     Pi[j-ma+jdim*(j-mb)] = -Pi[ma+jdim*mb];
//                 } else {
//                     Pr[j-ma+jdim*(j-mb)] = -Pr[ma+jdim*mb];
//                     Pi[j-ma+jdim*(j-mb)] = Pi[ma+jdim*mb];
//                 }
//                 mapar = -mapar;
//             }
//             mbpar = -mbpar;
//         }
//         
//         for (int mb = 0; mb <= j; mb++) 
//             for (int ma = 0; ma <= j; ma++) {
//                 //if (2*mb <= j) { // j = 3 -> mb = 0, 1
//                     //printf("%i %i %i %i %i\n", j, ma, mb, jju, k);
//                     int in = i + inum*k;                
//                     Stotr[in] += sfac*Pr[ma+jdim*mb]; // atomic add   
//                     Stoti[in] += sfac*Pi[ma+jdim*mb]; // atomic add                       
//                 //}
//                 k += 1; 
//             }           
//     }               
//   }
// };
// template void cpuComputeUi(double*, double*, double*, double*, double*, double*, double, double, double, 
//         int*, int*, int*, int*, int*, int, int, int, int, int, int);
// template void cpuComputeUi(float*, float*, float*, float*, float*, float*, float, float, float, 
//         int*, int*, int*, int*, int*, int, int, int, int, int, int);

#endif
        