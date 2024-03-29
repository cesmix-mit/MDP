/***************************************************************************                               
                    Molecular Dynamics Potentials (MDP)
                           CESMIX-MIT Project  
 
 Contributing authors: Ngoc-Cuong Nguyen (cuongng@mit.edu, exapde@gmail.com)
 ***************************************************************************/

#ifndef MDP_CPUSNAP2
#define MDP_CPUSNAP2
        
// Sr, Si, Srx, Sry, Srz, Six, Siy, Siz [ijnum*idxu_max]
// Stotr, Stoti [inum*idxu_max]       
// Ztotr, Ztoti [inum*idxz_max*nelements*nelements]        
// blist [inum*idxb_max*nelements*nelements*nelements]        
// ylist [inum*idxu_max*nelements]           
// beta [inum*ncoeff]         
// dedr [ijnum*3]        
//   if (chem_flag)
//     nelements = nelements_in;
//   else
//     nelements = 1;        
        
template <typename T> void cpuComputeSij(T *Sr, T *Si, T *Srx, T *Six, T *Sry, T *Siy, T *Srz, T *Siz, 
        T *rootpqarray, T *rij, T *wjelem, T *radelem, T rmin0, T rfac0, T rcutfac, int *idxu_block,  
        int *ti, int *tj, int twojmax, int idxu_max, int ijnum, int switch_flag)                
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
    T rscale0 = rfac0 * M_PI / (rcutij - rmin0);
    T theta0 = (r - rmin0) * rscale0;
    T z0 = r / tan(theta0);                
    T dz0dr = z0 / r - (r*rscale0) * (rsq + z0 * z0) / rsq;
            
    T sfac = 0.0, dsfac = 0.0;        
    if (switch_flag == 0) {
        sfac = 1.0;
        dsfac = 0.0;
    }
    else if (switch_flag == 1) {
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

    //sfac = 1.0; 
    //dsfac = 0.0;
    
    T r0inv, dr0invdr;
    T a_r, a_i, b_r, b_i;
    T da_r[3], da_i[3], db_r[3], db_i[3];
    T dz0[3], dr0inv[3];
    T rootpq;
    int jdim = twojmax + 1;
  
    r0inv = 1.0 / sqrt(r * r + z0 * z0);
    a_r = r0inv * z0;
    a_i = -r0inv * z;
    b_r = r0inv * y;
    b_i = -r0inv * x;

    dr0invdr = -pow(r0inv, 3.0) * (r + z0 * dz0dr);

    dr0inv[0] = dr0invdr * ux;
    dr0inv[1] = dr0invdr * uy;
    dr0inv[2] = dr0invdr * uz;

    dz0[0] = dz0dr * ux;
    dz0[1] = dz0dr * uy;
    dz0[2] = dz0dr * uz;

    for (int k = 0; k < 3; k++) {
        da_r[k] = dz0[k] * r0inv + z0 * dr0inv[k];
        da_i[k] = -z * dr0inv[k];
    }
    da_i[2] += -r0inv;

    for (int k = 0; k < 3; k++) {
        db_r[k] = y * dr0inv[k];
        db_i[k] = -x * dr0inv[k];
    }
    db_i[0] += -r0inv;
    db_r[1] += r0inv;
    
    Sr[ij+0*ijnum] = 1.0;
    Si[ij+0*ijnum] = 0.0;
    Srx[ij+0*ijnum] = 0.0;
    Six[ij+0*ijnum] = 0.0;
    Sry[ij+0*ijnum] = 0.0;
    Siy[ij+0*ijnum] = 0.0;
    Srz[ij+0*ijnum] = 0.0;
    Siz[ij+0*ijnum] = 0.0;
    for (int j = 1; j <= twojmax; j++) {
        int jju = idxu_block[j];
        int jjup = idxu_block[j-1];
        
        // fill in left side of matrix layer from previous layer
        for (int mb = 0; 2*mb <= j; mb++) {
            Sr[ij+jju*ijnum] = 0.0;
            Si[ij+jju*ijnum] = 0.0;
            Srx[ij+jju*ijnum] = 0.0;
            Six[ij+jju*ijnum] = 0.0;
            Sry[ij+jju*ijnum] = 0.0;
            Siy[ij+jju*ijnum] = 0.0;
            Srz[ij+jju*ijnum] = 0.0;
            Siz[ij+jju*ijnum] = 0.0;
            for (int ma = 0; ma < j; ma++) {
                rootpq = rootpqarray[(j - ma)*jdim + (j - mb)];
                int njju = ij+jju*ijnum;
                int njju1 = ij+(jju+1)*ijnum;
                int njjup = ij+jjup*ijnum;
                T u_r = Sr[njjup];
                T u_i = Si[njjup];
                T ux_r = Srx[njjup];
                T ux_i = Six[njjup];
                T uy_r = Sry[njjup];
                T uy_i = Siy[njjup];
                T uz_r = Srz[njjup];
                T uz_i = Siz[njjup];

                Sr[njju] += rootpq * (a_r * u_r + a_i * u_i);
                Si[njju] += rootpq * (a_r * u_i - a_i * u_r);
                Srx[njju] += rootpq * (da_r[0] * u_r + da_i[0] * u_i + a_r * ux_r + a_i * ux_i);
                Six[njju] += rootpq * (da_r[0] * u_i - da_i[0] * u_r + a_r * ux_i - a_i * ux_r);
                Sry[njju] += rootpq * (da_r[1] * u_r + da_i[1] * u_i + a_r * uy_r + a_i * uy_i);
                Siy[njju] += rootpq * (da_r[1] * u_i - da_i[1] * u_r + a_r * uy_i - a_i * uy_r);
                Srz[njju] += rootpq * (da_r[2] * u_r + da_i[2] * u_i + a_r * uz_r + a_i * uz_i);
                Siz[njju] += rootpq * (da_r[2] * u_i - da_i[2] * u_r + a_r * uz_i - a_i * uz_r);

                rootpq = rootpqarray[(ma + 1)*jdim + (j - mb)];
                Sr[njju1] = -rootpq * (b_r * u_r + b_i * u_i);
                Si[njju1] = -rootpq * (b_r * u_i - b_i * u_r);
                Srx[njju1] = -rootpq * (db_r[0] * u_r + db_i[0] * u_i + b_r * ux_r + b_i * ux_i);
                Six[njju1] = -rootpq * (db_r[0] * u_i - db_i[0] * u_r + b_r * ux_i - b_i * ux_r);
                Sry[njju1] = -rootpq * (db_r[1] * u_r + db_i[1] * u_i + b_r * uy_r + b_i * uy_i);
                Siy[njju1] = -rootpq * (db_r[1] * u_i - db_i[1] * u_r + b_r * uy_i - b_i * uy_r);
                Srz[njju1] = -rootpq * (db_r[2] * u_r + db_i[2] * u_i + b_r * uz_r + b_i * uz_i);
                Siz[njju1] = -rootpq * (db_r[2] * u_i - db_i[2] * u_r + b_r * uz_i - b_i * uz_r);
                jju++;
                jjup++;
            }
            jju++;
        }

    // copy left side to right side with inversion symmetry VMK 4.4(2)
    // u[ma-j][mb-j] = (-1)^(ma-mb)*Conj([u[ma][mb])
                   
        jju = idxu_block[j];
        jjup = jju+(j+1)*(j+1)-1;
        int mbpar = 1;
        for (int mb = 0; 2*mb <= j; mb++) {
            int mapar = mbpar;
            for (int ma = 0; ma <= j; ma++) {
                int njju = ij+jju*ijnum;
                int njjup = ij+jjup*ijnum;
                if (mapar == 1) {
                    Sr[njjup] = Sr[njju];
                    Si[njjup] = -Si[njju];
                    if (j%2==1 && mb==(j/2)) {                    
                    Srx[njjup] =  Srx[njju];
                    Six[njjup] = -Six[njju];
                    Sry[njjup] =  Sry[njju];
                    Siy[njjup] = -Siy[njju];
                    Srz[njjup] =  Srz[njju];
                    Siz[njjup] = -Siz[njju];
                    }
                } else {
                    Sr[njjup] = -Sr[njju];
                    Si[njjup] =  Si[njju];
                    if (j%2==1 && mb==(j/2)) {
                    Srx[njjup] = -Srx[njju];
                    Six[njjup] =  Six[njju];
                    Sry[njjup] = -Sry[njju];
                    Siy[njjup] =  Siy[njju];
                    Srz[njjup] = -Srz[njju];
                    Siz[njjup] =  Siz[njju];                    
                    }
                }
                mapar = -mapar;
                jju++;
                jjup--;
            }
            mbpar = -mbpar;
        }        
    }        

    for (int j = 0; j <= twojmax; j++) {
        int jju = idxu_block[j];
        for (int mb = 0; 2*mb <= j; mb++)
          for (int ma = 0; ma <= j; ma++) {
            int ijk = ij+jju*ijnum;               
            Srx[ijk] = dsfac * Sr[ijk] * ux + sfac * Srx[ijk]; 
            Six[ijk] = dsfac * Si[ijk] * ux + sfac * Six[ijk]; 
            Sry[ijk] = dsfac * Sr[ijk] * uy + sfac * Sry[ijk]; 
            Siy[ijk] = dsfac * Si[ijk] * uy + sfac * Siy[ijk]; 
            Srz[ijk] = dsfac * Sr[ijk] * uz + sfac * Srz[ijk]; 
            Siz[ijk] = dsfac * Si[ijk] * uz + sfac * Siz[ijk];                  
            jju++;
          }
    }
    
    for (int k=0; k<idxu_max; k++) {
        int ijk = ij + ijnum*k;
        Sr[ijk] = sfac*Sr[ijk];
        Si[ijk] = sfac*Si[ijk];
    }            
  }
};
template void cpuComputeSij(double*, double*, double*, double*, double*, double*, double*, double*, 
        double*, double*, double*, double*, double, double, double, int*, int*, int*, int, int, int, int);
template void cpuComputeSij(float*, float*, float*, float*, float*, float*, float*, float*, 
        float*, float*, float*, float*, float, float, float, int*, int*, int*, int, int, int, int);

template <typename T> void cpuZeroUarraytot2(T *Stotr, T *Stoti, T wself, int *idxu_block, 
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
        int jju = idxu_block[j];        
        for (int mb = 0; mb <= j; mb++) {
            for (int ma = 0; ma <= j; ma++) {
                int n = ii + inum*jju + inum*idxu_max*jelem;        
                Stotr[n] = 0.0;
                Stoti[n] = 0.0;
                if (jelem == ielem || wselfall_flag)
                    if (ma==mb)
                        Stotr[n] = wself; ///// T check this
                jju++;
            }
        }
        
    }                    
};
template void cpuZeroUarraytot2(double*, double*, double, int*, int*, 
        int*, int*, int, int, int, int, int, int);
template void cpuZeroUarraytot2(float*, float*, float, int*, int*, 
        int*, int*, int, int, int, int, int, int);

template <typename T> void cpuKernelAddUarraytot(T *Stotr, T *Stoti, T *Sr, T *Si, 
        int *pairnum, int *pairnumsum, int *map, int *tj, 
        int idxu_max, int inum, int ijnum, int N2, int chemflag)
{    
    for (int idx=0; idx < N2; idx++) {
        int ii = idx%inum;  // inum
        int jju = (idx-ii)/inum;    // idxu_max
        int m = pairnum[ii];        // number of neighbors around i             
        int start = pairnumsum[ii];         
        int n = inum*idxu_max;
        for (int l=0; l<m; l++) {   // loop over each atom around atom i    
            int k = start + l;                
            int jelem = (chemflag) ? map[tj[k]] : 0;                                        
            int kl = idx + n*jelem;        
            int kr = k + ijnum*jju;        
            Stotr[kl] += Sr[kr];
            Stoti[kl] += Si[kr];                
        }        
    }
};
template <typename T>  void cpuKernelAddUarraytot(T *Stotr, T *Stoti, T *Sr, T *Si, 
        int *pairnum, int *pairnumsum, int *map, int *tj, int idxu_max, int inum, int ijnum, int N2)
{    
    for (int idx=0; idx < N2; idx++) {
        int ii = idx%inum;  // inum
        int jju = (idx-ii)/inum;    // idxu_max
        int m = pairnum[ii];        // number of neighbors around i             
        int start = pairnumsum[ii];         
        int n = start + ijnum*jju;
        T qr = 0.0, qi = 0.0;    
        for (int l=0; l<m; l++) {   // loop over each atom around atom i                                       
            qr += Sr[l + n];
            qi += Si[l + n];
//             qr += Sr[idx+N2*l];
//             qi += Si[idx+N2*l];        
        }
        Stotr[idx] += qr;
        Stoti[idx] += qi;                        
    }
};
template <typename T> void cpuAddUarraytot(T *Stotr, T *Stoti, T *Sr, 
        T *Si, int *pairnum, int *pairnumsum, int *map, int *tj, 
        int idxu_max, int nelements, int inum, int ijnum, int chemflag)
{    
    int N2 = inum*idxu_max;    
    if (chemflag==0) {
        cpuKernelAddUarraytot(Stotr, Stoti, Sr, Si, 
            pairnum, pairnumsum, map, tj, idxu_max, inum, ijnum, N2);          
    } else
        cpuKernelAddUarraytot(Stotr, Stoti, Sr, Si, 
                pairnum, pairnumsum, map, tj, idxu_max, inum, ijnum, N2, chemflag);  
}
template void cpuAddUarraytot(double*, double*, double*, double*, 
        int*, int*, int*, int*, int, int, int, int, int);
template void cpuAddUarraytot(float*, float*, float*, float*, 
        int*, int*, int*, int*, int, int, int, int, int);

template <typename T>  void cpuKernelAddUarraytot(T *Stotr, T *Stoti, T *Sr, T *Si, 
        int *map, int *ai, int *tj, int inum, int ijnum, int N1, int N2, int chemflag)
{    
    for (int idx=0; idx < N2; idx++) {
        int ij = idx%ijnum;  // ijnum
        int jju = (idx-ij)/ijnum;    // idxu_max   
        int jelem = (chemflag) ? map[tj[ij]] : 0;     
        int i = ai[ij] + inum*jju + N1*jelem;                
        Stotr[i] += Sr[idx];
        Stoti[i] += Si[idx];                
    }
};
template <typename T>  void cpuKernelAddUarraytot(T *Stotr, T *Stoti, T *Sr, T *Si, 
        int *map, int *ai, int *tj, int inum, int ijnum, int N1, int N2)
{    
    for (int idx=0; idx < N2; idx++) {
        int ij = idx%ijnum;  // ijnum
        int jju = (idx-ij)/ijnum;    // idxu_max        
        int i = ai[ij] + inum*jju;                
        Stotr[i] += Sr[idx];
        Stoti[i] += Si[idx];                
    }
};
template <typename T> void cpuAddUarraytot(T *Stotr, T *Stoti, T *Sr, 
        T *Si, int *map, int *ai, int *tj, int idxu_max, int inum, int ijnum, int chemflag)
{   
    int N1 = inum*idxu_max;    
    int N2 = ijnum*idxu_max;    
    if (chemflag==0) {
        cpuKernelAddUarraytot(Stotr, Stoti, Sr, Si, map, ai, tj, inum, ijnum, N1, N2);          
    } else
        cpuKernelAddUarraytot(Stotr, Stoti, Sr, Si, map, ai, tj, inum, ijnum, N1, N2, chemflag);  
}
template void cpuAddUarraytot(double*, double*, double*, double*, 
        int*, int*, int*, int, int, int, int);
template void cpuAddUarraytot(float*, float*, float*, float*, 
        int*, int*, int*, int, int, int, int);

template <typename T> void cpuComputeZi2(T *zlist_r, T *zlist_i, T *Stotr, T *Stoti, 
        T *cglist, int *idxz, int *idxu_block, int *idxcg_block, int twojmax, int idxu_max, 
        int idxz_max, int nelements, int bnorm_flag, int inum)
{
    int jdim = twojmax + 1;    
    int N1 = inum;    
    int N2 = N1*idxz_max;
    int N3 = N2*nelements*nelements;                                
    for (int idx=0; idx < N3; idx++) {
        int l = idx%N2;   //  inum*idxz_max
        int ii = l%inum;    // inum
        int jjz = (l-ii)/inum; // idxz_max
        int ielem = (idx-l)/N2;  // nelements*nelements  
        int elem2 = ielem%nelements; // nelements
        int elem1 = (ielem-elem2)/nelements; // nelements
              
        const int j1 = idxz[jjz*10+0];
        const int j2 = idxz[jjz*10+1];
        const int j = idxz[jjz*10+2];
        const int ma1min = idxz[jjz*10+3];
        const int ma2max = idxz[jjz*10+4];
        const int na = idxz[jjz*10+5];
        const int mb1min = idxz[jjz*10+6];
        const int mb2max = idxz[jjz*10+7];
        const int nb = idxz[jjz*10+8];
        int jju1 = idxu_block[j1] + (j1 + 1) * mb1min;
        int jju2 = idxu_block[j2] + (j2 + 1) * mb2max;
        int icgb = mb1min * (j2 + 1) + mb2max;

        const T *cgblock = cglist + idxcg_block[j + j2*jdim + j1*jdim*jdim];
        T qr = 0.0;
        T qi = 0.0;          
        for (int ib = 0; ib < nb; ib++) {
            T suma1_r = 0.0;
            T suma1_i = 0.0;

            // Stotr: inum*idxu_max*nelements  
            const T *u1_r = &Stotr[ii + inum*jju1 + inum*idxu_max*elem1];
            const T *u1_i = &Stoti[ii + inum*jju1 + inum*idxu_max*elem1];
            const T *u2_r = &Stotr[ii + inum*jju2 + inum*idxu_max*elem2];
            const T *u2_i = &Stoti[ii + inum*jju2 + inum*idxu_max*elem2];

            int ma1 = ma1min;
            int ma2 = ma2max;
            int icga = ma1min * (j2 + 1) + ma2max;

            for (int ia = 0; ia < na; ia++) {
                suma1_r += cgblock[icga] * (u1_r[inum*ma1] * u2_r[inum*ma2] - u1_i[inum*ma1] * u2_i[inum*ma2]);
                suma1_i += cgblock[icga] * (u1_r[inum*ma1] * u2_i[inum*ma2] + u1_i[inum*ma1] * u2_r[inum*ma2]);
                ma1++;
                ma2--;
                icga += j2;
            } // end loop over ia

            qr += cgblock[icgb] * suma1_r;
            qi += cgblock[icgb] * suma1_i;

            jju1 += j1 + 1;
            jju2 -= j2 + 1;
            icgb += j2;
        } // end loop over ib
        
        if (bnorm_flag) {
            qr /= (j+1);
            qi /= (j+1);
        }        
        
        zlist_r[idx] = qr;
        zlist_i[idx] = qi;          
    }
};
template void cpuComputeZi2(double*, double*, double*, double*, double*, 
        int*, int*, int*, int, int, int, int, int, int);
template void cpuComputeZi2(float*, float*, float*, float*, float*, 
        int*, int*, int*, int, int, int, int, int, int);

template <typename T> void cpuKernelComputeBi1(T *blist, T *zlist_r, T *zlist_i, 
        T *Stotr, T *Stoti, int *idxb, int *idxu_block, int *idxz_block, int jdim,         
        int nelements, int nelemsq, int nz_max, int nu_max, int nb_max, int inum, int N2, int N3)
{    
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

        blist[idx] = 2.0 * sumzu;                
    }
}
template <typename T> void cpuKernelComputeBi2(T *blist, T *bzero,int *ilist, int *type,
       int *map, int *idxb, int nelements, int nb_max, int inum, int N2, int chemflag)
{        
    for (int idx=0; idx < N2; idx++) {        
        int ii = idx%inum;        
        int jjb = (idx-ii)/inum;    
        
        int ielem = (chemflag) ? map[type[ilist[ii]]]: 0;                
        int itriple = (ielem*nelements+ielem)*nelements+ielem;

        const int j = idxb[jjb*3 + 2];  
        blist[ii + inum*jjb + nb_max*itriple] -= bzero[j];                
    }
}
template <typename T> void cpuKernelComputeBi4(T *blist, T *bzero,
       int *idxb, int inum, int N2, int N3)
{        
    for (int idx=0; idx < N3; idx++) {
        int l = idx%N2;
        int ii = l%inum;        
        int jjb = (l-ii)/inum;    
        int j = idxb[jjb*3 + 2];  
        blist[idx] -= bzero[j];        
    }
}
template <typename T> void cpuComputeBi2(T *blist, T *zlist_r, T *zlist_i, T *Stotr, T *Stoti, 
        T *bzero, int *ilist, int *type, int *map, int *idxb, int *idxu_block, int *idxz_block, 
        int twojmax, int idxb_max, int idxu_max, int idxz_max, int nelements, int bzero_flag, 
        int wselfall_flag, int chemflag, int inum)
{                
    int nelemsq = nelements*nelements;
    int nz_max = idxz_max*inum;
    int nu_max = idxu_max*inum;
    int nb_max = idxb_max*inum;
    int N2 = inum*idxb_max;
    int N3 = N2*nelements*nelemsq;
    int jdim = twojmax+1;

    cpuKernelComputeBi1(blist, zlist_r, zlist_i, Stotr, Stoti, 
            idxb, idxu_block, idxz_block, jdim, nelements, nelemsq, nz_max, nu_max, nb_max, inum, N2, N3);

    if (bzero_flag) {
        if (!wselfall_flag) {
            cpuKernelComputeBi2(blist, bzero, ilist, type, map, 
                    idxb, nelements, nb_max, inum, N2, chemflag);
        }
        else {
            cpuKernelComputeBi4(blist, bzero, idxb, inum, N2, N3);            
        }
    }
}
template void cpuComputeBi2(double*, double*, double*, double*, double*, double*, 
        int*, int*, int*, int*, int*, int*, int, int, int, int, int, int, int, int, int);
template void cpuComputeBi2(float*, float*, float*, float*, float*, float*,
        int*, int*, int*, int*, int*, int*, int, int, int, int, int, int, int, int, int);

template <typename T> void cpuKernelComputeBeta1(T *beta, T *coeffelem, int *ilist, int *map, int *type, 
        int inum, int ncoeff, int ncoeffall, int N2)
{    
    for (int idx=0; idx < N2; idx++) {
        int ii = idx%inum;              
        int icoeff = (idx-ii)/inum;         
        int i = ilist[ii]; // index of atom i
        const int itype = type[i]; // element type of atom i
        const int ielem = map[itype];  // index of that element type
        beta[icoeff*inum + ii] = coeffelem[icoeff+1+ielem*ncoeffall];
    }
}
template <typename T> void cpuKernelComputeBeta2(T *beta, T *bispectrum, T *coeffelem, int *ilist, int *map, int *type, 
        int inum, int ncoeff, int ncoeffall, int N2)
{    
    for (int idx=0; idx < N2; idx++) {
        int ii = idx%inum;              
        int icoeff = (idx-ii)/inum;         
        int i = ilist[ii]; // index of atom i
        const int itype = type[i]; // element type of atom i
        const int ielem = map[itype];  // index of that element type
        T bveci = bispectrum[icoeff*inum + ii];        
        int k = ncoeff + 1 + icoeff*ncoeff - (icoeff-1)*icoeff/2;
        beta[icoeff*inum + ii] = coeffelem[k+ielem*ncoeffall];
        for (int jcoeff = icoeff+1; jcoeff < ncoeff; jcoeff++) {
          T bvecj = bispectrum[jcoeff*inum + ii];
          beta[icoeff*inum + ii] += coeffelem[k+ielem*ncoeffall]*bvecj;
          beta[jcoeff*inum + ii] += coeffelem[k+ielem*ncoeffall]*bveci;
          k++;
        }
    }
}
template <typename T> void cpuComputeBeta(T *beta, T *bispectrum, T *coeffelem, int *ilist, int *map, int *type, 
        int inum, int ncoeff, int ncoeffall, int quadraticflag)
{
    int N1 = inum;
    int N2 = ncoeff*inum;
    cpuKernelComputeBeta1(beta, coeffelem, ilist,
            map, type, inum, ncoeff, ncoeffall, N2);

    if (quadraticflag) {
        cpuKernelComputeBeta2(beta, bispectrum, coeffelem, ilist,
            map, type, inum, ncoeff, ncoeffall, N2);
    }
}
template void cpuComputeBeta(double*, double*, double*, int*, int*, int*,
        int, int, int, int);
template void cpuComputeBeta(float*, float*, float*, int*, int*, int*,
        int, int, int, int);

// template <typename T> void cpuComputeYi(T *ylist_r, T *ylist_i, T *Stotr, T *Stoti, T *cglist, 
//         T* beta, int *idxz, int *idxb_block, int *idxu_block, int *idxcg_block, int twojmax, 
//         int idxb_max, int idxu_max, int idxz_max, int nelements, int bnorm_flag, int inum)
// {    
//     int N1 = idxu_max*nelements*inum;
//     cpuArraySetValue(ylist_r, (T) 0.0, N1);
//     cpuArraySetValue(ylist_i, (T) 0.0, N1);
//     
//     int jdim = twojmax + 1;         
//     int N2 = idxz_max*inum;                        
//   
//   for (int idx=0; idx < N2; idx++) {
//       int ii = idx%inum;              
//       int jjz = (idx-ii)/inum;                    
//       for(int elem1 = 0; elem1 < nelements; elem1++)
//         for (int elem2 = 0; elem2 < nelements; elem2++) {        
//             const int j1 = idxz[jjz*10+0];
//             const int j2 = idxz[jjz*10+1];
//             const int j = idxz[jjz*10+2];
//             const int ma1min = idxz[jjz*10+3];
//             const int ma2max = idxz[jjz*10+4];
//             const int na = idxz[jjz*10+5];
//             const int mb1min = idxz[jjz*10+6];
//             const int mb2max = idxz[jjz*10+7];
//             const int nb = idxz[jjz*10+8];
//           
//             const T *cgblock = cglist + idxcg_block[j + j2*jdim + j1*jdim*jdim];
//           
//             T ztmp_r = 0.0;
//             T ztmp_i = 0.0;
// 
//             int jju1 = idxu_block[j1] + (j1 + 1) * mb1min;
//             int jju2 = idxu_block[j2] + (j2 + 1) * mb2max;
//             int icgb = mb1min * (j2 + 1) + mb2max;
//             for (int ib = 0; ib < nb; ib++) {
// 
//                 T suma1_r = 0.0;
//                 T suma1_i = 0.0;
// 
//                 const T *u1_r = &Stotr[ii + inum*jju1 + inum*idxu_max*elem1];
//                 const T *u1_i = &Stoti[ii + inum*jju1 + inum*idxu_max*elem1];
//                 const T *u2_r = &Stotr[ii + inum*jju2 + inum*idxu_max*elem2];
//                 const T *u2_i = &Stoti[ii + inum*jju2 + inum*idxu_max*elem2];
// 
//                 int ma1 = ma1min;
//                 int ma2 = ma2max;
//                 int icga = ma1min * (j2 + 1) + ma2max;
// 
//                 for (int ia = 0; ia < na; ia++) {
//                     suma1_r += cgblock[icga] * (u1_r[inum*ma1] * u2_r[inum*ma2] - u1_i[inum*ma1] * u2_i[inum*ma2]);
//                     suma1_i += cgblock[icga] * (u1_r[inum*ma1] * u2_i[inum*ma2] + u1_i[inum*ma1] * u2_r[inum*ma2]);                    
//                   ma1++;
//                   ma2--;
//                   icga += j2;
//                 } // end loop over ia
// 
//                 ztmp_r += cgblock[icgb] * suma1_r;
//                 ztmp_i += cgblock[icgb] * suma1_i;
// 
//                 jju1 += j1 + 1;
//                 jju2 -= j2 + 1;
//                 icgb += j2;
//             } // end loop over ib
// 
// 
//             if (bnorm_flag) {
//               ztmp_i /= j+1;
//               ztmp_r /= j+1;
//             }
//        
//             int jju = idxz[jjz*10+9];
//             for(int elem3 = 0; elem3 < nelements; elem3++) {
//               int itriple;  
//               T betaj;
//               if (j >= j1) {
//                 const int jjb = idxb_block[j + j2*jdim + j1*jdim*jdim];
//                 itriple = ((elem1 * nelements + elem2) * nelements + elem3) * idxb_max*inum + jjb*inum + ii;
//                 if (j1 == j) {
//                   if (j2 == j) betaj = 3*beta[itriple];
//                   else betaj = 2*beta[itriple];
//                 } else betaj = beta[itriple];          
//               } else if (j >= j2) {
//                 const int jjb = idxb_block[j1 + j2*jdim + j*jdim*jdim];
//                 itriple = ((elem3 * nelements + elem2) * nelements + elem1) * idxb_max*inum + jjb*inum + ii;
//                 if (j2 == j) betaj = 2*beta[itriple];
//                 else betaj = beta[itriple];
//               } else {
//                 const int jjb = idxb_block[j1 + j*jdim + j2*jdim*jdim];
//                 itriple = ((elem2 * nelements + elem3) * nelements + elem1) * idxb_max*inum + jjb*inum + ii;
//                 betaj = beta[itriple];
//               }
//               
//               if (!bnorm_flag && j1 > j)
//                 betaj *= (j1 + 1) / (j + 1.0);
//            
//               ylist_r[ii + inum*jju + inum*idxu_max*elem3] += betaj * ztmp_r;
//               ylist_i[ii + inum*jju + inum*idxu_max*elem3] += betaj * ztmp_i;        
//            }
//         }         
//     }  
// }
// template void cpuComputeYi(double*, double*, double*, double*, double*, double*, 
//         int*, int*, int*, int*, int, int, int, int, int, int, int);
// template void cpuComputeYi(float*, float*, float*, float*, float*, float*,
//         int*, int*, int*, int*, int, int, int, int, int, int, int);

template <typename T> void cpuComputeYi(T *ylist_r, T *ylist_i, T *zlist_r, T *zlist_i, 
        T* beta, int *idxz, int *idxb_block, int twojmax, int idxb_max, int idxu_max, int idxz_max,
        int nelements, int bnorm_flag, int inum)
{        
    cpuArraySetValue(ylist_r, (T) 0.0, idxu_max*nelements*inum);
    cpuArraySetValue(ylist_i, (T) 0.0, idxu_max*nelements*inum);
            
    int jdim = twojmax + 1;    
    int N1 = inum;    
    int N2 = N1*idxz_max;
    int N3 = N2*nelements*nelements;                                    
    for (int idx=0; idx < N3; idx++) {
        int l = idx%N2;   //  inum*idxz_max
        int ii = l%inum;    // inum
        int jjz = (l-ii)/inum; // idxz_max
        int ielem = (idx-l)/N2;  // nelements*nelements  
        int elem2 = ielem%nelements; // nelements
        int elem1 = (ielem-elem2)/nelements; // nelements

        const int j1 = idxz[jjz*10+0];
        const int j2 = idxz[jjz*10+1];
        const int j = idxz[jjz*10+2];
        const int jju = idxz[jjz*10+9];
        
        T ztmp_r = zlist_r[idx];
        T ztmp_i = zlist_i[idx];                  
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
template void cpuComputeYi(double*, double*, double*, double*, double*, 
        int*, int*, int, int, int, int, int, int, int);
template void cpuComputeYi(float*, float*, float*, float*, float*, 
        int*, int*, int, int, int, int, int, int, int);

template <typename T> void cpuComputeDeidrj(T *dedr, T *ylist_r, T *ylist_i, 
        T *dulist_r, T *dulist_i, int *idxu_block, int *map, int *ai, int *tj,
        int twojmax, int idxu_max, int chemflag, int inum, int ijnum) 
{                
    cpuArraySetValue(dedr, (T) 0.0, ijnum*3);
    
    int niimax = inum*idxu_max;
    int nijmax = ijnum*idxu_max;  
  for(int ij=0; ij<ijnum; ij++) {                        
      int jelem = (chemflag) ? map[tj[ij]] : 0; //(chemflag) ? map[type[alist[aj[ij]]]] : 0;
      int nimax = niimax*jelem;              
      int i = ai[ij]; // atom i        
      T de[3];
      de[0] = 0.0; 
      de[1] = 0.0; 
      de[2] = 0.0;
      for(int j = 0; j <= twojmax; j++) {
        int jju = idxu_block[j];        
        
        for(int mb = 0; 2*mb < j; mb++)
          for(int ma = 0; ma <= j; ma++) {
            int n1 = i + inum*jju + nimax;
            int n2 = ij+ ijnum*jju;
            T y_r = ylist_r[n1];
            T y_i = ylist_i[n1];
            for(int k = 0; k < 3; k++)
              de[k] += dulist_r[n2+nijmax*k] * y_r + dulist_i[n2+nijmax*k] * y_i;
            jju++;
          } //end loop over ma mb

        // For j even, handle middle column

        if (j%2 == 0) {
          int mb = j/2;
          for(int ma = 0; ma < mb; ma++) {
            int n1 = i + inum*jju + nimax;
            int n2 = ij+ ijnum*jju;
            T y_r = ylist_r[n1];
            T y_i = ylist_i[n1];
            for(int k = 0; k < 3; k++)
                de[k] += dulist_r[n2+nijmax*k] * y_r + dulist_i[n2+nijmax*k] * y_i;
            jju++;
          }

          int n1 = i + inum*jju + nimax;
          int n2 = ij+ ijnum*jju;
          T y_r = ylist_r[n1];
          T y_i = ylist_i[n1];
          for(int k = 0; k < 3; k++)
            de[k] += (dulist_r[n2+nijmax*k] * y_r + dulist_i[n2+nijmax*k] * y_i)*0.5;
          // jju++;
        } // end if jeven
      } // end loop over j
      
      for(int k = 0; k < 3; k++)
        dedr[k + 3*ij] = 2.0*de[k];      
  }  
}
template void cpuComputeDeidrj(double*, double*, double*, double*, double*, 
        int*, int*, int*, int*, int, int, int, int, int);
template void cpuComputeDeidrj(float*, float*, float*, float*, float*, 
        int*, int*, int*, int*, int, int, int, int, int);

template <typename T> void cpuComputeDbidrj(T *dblist, T *zlist_r, T *zlist_i, 
        T *dulist_r, T *dulist_i, int *idxb, int *idxu_block, int *idxz_block, 
        int *map, int *ai, int *tj, int twojmax, int idxb_max, int idxu_max, int idxz_max, 
        int nelements, int bnorm_flag, int chemflag, int inum, int ijnum)
{                    
    int nz_max = idxz_max*inum;
    int nb_max = idxb_max*ijnum;    
    int nu_max = idxu_max*ijnum;    
    int N2 = ijnum*idxb_max;
    int jdim = twojmax+1;

    cpuArraySetValue(dblist, (T) 0.0, nb_max*3*nelements*nelements*nelements);
        
    for (int idx=0; idx < N2; idx++) {
        int ij = idx%ijnum;              
        int jjb = (idx-ij)/ijnum;                              
        int elem3 = (chemflag) ? map[tj[ij]] : 0;//(chemflag) ? map[type[alist[aj[ij]]]] : 0;
        int i = ai[ij]; // atom i
        const int j1 = idxb[jjb*3 + 0];
        const int j2 = idxb[jjb*3 + 1];
        const int j = idxb[jjb*3 + 2];  

       // Sum terms Conj(dudr(j,ma,mb))*z(j1,j2,j,ma,mb)
        for(int elem1 = 0; elem1 < nelements; elem1++)
            for(int elem2 = 0; elem2 < nelements; elem2++) {

            int jjz = idxz_block[j + j2*jdim + j1*jdim*jdim];
            int jju = idxu_block[j];
            int idouble = elem1*nelements+elem2;
            int itriple = (elem1*nelements+elem2)*nelements+elem3;
            int nimax = nz_max*idouble;                      

            T *dbdr = &dblist[nb_max*3*itriple];
            T sumzdu_r[3];
            for (int k = 0; k < 3; k++)
              sumzdu_r[k] = 0.0;

            for (int mb = 0; 2 * mb < j; mb++)
              for (int ma = 0; ma <= j; ma++) {
                int n1 = i + inum*jjz + nimax;
                int n2 = ij+ ijnum*jju;
                T z_r = zlist_r[n1];
                T z_i = zlist_i[n1];
                for (int k = 0; k < 3; k++)
                  sumzdu_r[k] += dulist_r[n2 + nu_max*k] * z_r + dulist_i[n2 + nu_max*k] * z_i;
                jjz++;
                jju++;
              } //end loop over ma mb

            // For j even, handle middle column

            if (j % 2 == 0) {
              int mb = j / 2;
              for (int ma = 0; ma < mb; ma++) {
                int n1 = i + inum*jjz + nimax;
                int n2 = ij+ ijnum*jju;
                T z_r = zlist_r[n1];
                T z_i = zlist_i[n1];
                for (int k = 0; k < 3; k++)
                  sumzdu_r[k] += dulist_r[n2 + nu_max*k] * z_r + dulist_i[n2 + nu_max*k] * z_i;
                jjz++;
                jju++;
              }

              int n1 = i + inum*jjz + nimax;
              int n2 = ij+ ijnum*jju;
              T z_r = zlist_r[n1];
              T z_i = zlist_i[n1];
              for (int k = 0; k < 3; k++)
                sumzdu_r[k] += (dulist_r[n2 + nu_max*k] * z_r + dulist_i[n2 + nu_max*k] * z_i) * 0.5;
              // jjz++;
              // jju++;
            } // end if jeven
                                
            for (int k = 0; k < 3; k++)
              dbdr[ij + ijnum*k + ijnum*3*jjb] += 2.0 * sumzdu_r[k];
            
            // Sum over Conj(dudr(j1,ma1,mb1))*z(j,j2,j1,ma1,mb1)
            T j1fac = (j + 1) / (j1 + 1.0);
            idouble = elem1*nelements+elem2;
            itriple = (elem3*nelements+elem2)*nelements+elem1;            
            //jjz = idxz_block[j][j2][j1];
            jjz = idxz_block[j1 + j2*jdim + j*jdim*jdim];
            jju = idxu_block[j1];

            //dbdr = &dblist[nb_max*3*itriple];
            for (int k = 0; k < 3; k++)
              sumzdu_r[k] = 0.0;

            for (int mb = 0; 2 * mb < j1; mb++)
              for (int ma = 0; ma <= j1; ma++) {
                int n1 = i + inum*jjz + nimax;
                int n2 = ij+ ijnum*jju;
                T z_r = zlist_r[n1];
                T z_i = zlist_i[n1];
                for (int k = 0; k < 3; k++)
                  sumzdu_r[k] += dulist_r[n2 + nu_max*k] * z_r + dulist_i[n2 + nu_max*k] * z_i;                       
                jjz++;
                jju++;
              } //end loop over ma mb

            // For j1 even, handle middle column

            if (j1 % 2 == 0) {
              int mb = j1 / 2;
              for (int ma = 0; ma < mb; ma++) {
                int n1 = i + inum*jjz + nimax;
                int n2 = ij+ ijnum*jju;
                T z_r = zlist_r[n1];
                T z_i = zlist_i[n1];
                for (int k = 0; k < 3; k++)
                  sumzdu_r[k] += dulist_r[n2 + nu_max*k] * z_r + dulist_i[n2 + nu_max*k] * z_i;                       
                jjz++;
                jju++;
              }
              int n1 = i + inum*jjz + nimax;
              int n2 = ij+ ijnum*jju;
              T z_r = zlist_r[n1];
              T z_i = zlist_i[n1];
              for (int k = 0; k < 3; k++)
                sumzdu_r[k] += (dulist_r[n2 + nu_max*k] * z_r + dulist_i[n2 + nu_max*k] * z_i) * 0.5;
              // jjz++;
              // jju++;
            } // end if j1even

            for (int k = 0; k < 3; k++)
              if (bnorm_flag)
                dbdr[ij + ijnum*k + ijnum*3*jjb] += 2.0 * sumzdu_r[k];
              else
                dbdr[ij + ijnum*k + ijnum*3*jjb] += 2.0 * sumzdu_r[k] * j1fac;

            // Sum over Conj(dudr(j2,ma2,mb2))*z(j,j1,j2,ma2,mb2)
            T j2fac = (j + 1) / (j2 + 1.0);
            idouble = elem2*nelements+elem1;
            itriple = (elem1*nelements+elem3)*nelements+elem2;
            //jjz = idxz_block[j][j1][j2];
            jjz = idxz_block[j2 + j1*jdim + j*jdim*jdim];        
            jju = idxu_block[j2];
            
            //dbdr = &dblist[nb_max*3*itriple];
            for (int k = 0; k < 3; k++)
              sumzdu_r[k] = 0.0;

            for (int mb = 0; 2 * mb < j2; mb++)
              for (int ma = 0; ma <= j2; ma++) {
                int n1 = i + inum*jjz + nimax;
                int n2 = ij+ ijnum*jju;
                T z_r = zlist_r[n1];
                T z_i = zlist_i[n1];
                for (int k = 0; k < 3; k++)
                  sumzdu_r[k] += dulist_r[n2 + nu_max*k] * z_r + dulist_i[n2 + nu_max*k] * z_i;                                            
                jjz++;
                jju++;
              } //end loop over ma mb

            // For j2 even, handle middle column

            if (j2 % 2 == 0) {
              int mb = j2 / 2;
              for (int ma = 0; ma < mb; ma++) {
                int n1 = i + inum*jjz + nimax;
                int n2 = ij+ ijnum*jju;
                T z_r = zlist_r[n1];
                T z_i = zlist_i[n1];
                for (int k = 0; k < 3; k++)
                  sumzdu_r[k] += dulist_r[n2 + nu_max*k] * z_r + dulist_i[n2 + nu_max*k] * z_i;                                                                 
                jjz++;
                jju++;
              }
              
              int n1 = i + inum*jjz + nimax;
              int n2 = ij+ ijnum*jju;
              T z_r = zlist_r[n1];
              T z_i = zlist_i[n1];
              for (int k = 0; k < 3; k++)
                sumzdu_r[k] += (dulist_r[n2 + nu_max*k] * z_r + dulist_i[n2 + nu_max*k] * z_i) * 0.5;
              // jjz++;
              // jju++;
            } // end if j2even

            for (int k = 0; k < 3; k++)
              if (bnorm_flag)
                dbdr[ij + ijnum*k + ijnum*3*jjb] += 2.0 * sumzdu_r[k];
              else
                dbdr[ij + ijnum*k + ijnum*3*jjb] += 2.0 * sumzdu_r[k] * j2fac;
          }        
    }
}
template void cpuComputeDbidrj(double*, double*, double*, double*, double*, 
        int*, int*, int*, int*, int*, int*, int, int, int, int, int, int, int, int, int);
template void cpuComputeDbidrj(float*, float*, float*, float*, float*, 
        int*, int*, int*, int*, int*, int*, int, int, int, int, int, int, int, int, int);

template <typename T>  void cpuKernelComputeSna1(T *sna, T *blist, int *ilist, 
        int *type, int ncoeff, int anum, int nperdim, int inum, int N2)
{        
    for (int idx=0; idx < N2; idx++) {
        int ii = idx%inum;              
        int icoeff = (idx-ii)/inum;                 
        int i = ilist[ii];
        int itype = type[i]-1;        
        sna[i + anum*icoeff + anum*nperdim*itype] = blist[ii + inum*icoeff];               
    }
}
template <typename T>  void cpuKernelComputeSna2(T *sna, T *blist, int *ilist, 
        int *type, int ncoeff, int anum, int nperdim, int inum, int N2)
{        
    for (int idx=0; idx < N2; idx++) {
        int ii = idx%inum;              
        int icoeff = (idx-ii)/inum;                 
        int i = ilist[ii];
        int itype = type[i]-1;        
        T bi = blist[ii + inum*icoeff];

        // diagonal element of quadratic matrix
        int ncount = ncoeff + icoeff*ncoeff - icoeff*(icoeff-1)/2;
        sna[i + anum*ncount + anum*nperdim*itype] = 0.5*bi*bi;
        ncount++;

        // upper-triangular elements of quadratic matrix
        for (int jcoeff = icoeff+1; jcoeff < ncoeff; jcoeff++) {
            sna[i + anum*ncount + anum*nperdim*itype] = bi*blist[jcoeff + ncoeff*ii];                        
            ncount++;
        }        
    }
}
template <typename T> void cpuComputeSna(T *sna, T *blist, int *ilist, int *type,
        int ncoeff, int anum, int nperdim, int inum, int quadraticflag)
{
    int N2 = ncoeff*inum;
    cpuKernelComputeSna1(sna, blist, ilist, type, ncoeff, anum, 
            nperdim, inum, N2);

    if (quadraticflag)
        cpuKernelComputeSna2(sna, blist, ilist, type, ncoeff, 
            anum, nperdim, inum, N2);
}
template void cpuComputeSna(double*, double*, int*, int*, int, int, int, int, int);
template void cpuComputeSna(float*, float*, int*, int*, int, int, int, int, int);

template <typename T>  void cpuKernelComputeSnad1(T *snad, T *dblist, 
        int *ai, int *aj, int *ti, int anum, int nperdim, int ijnum, int N2)
{                   
    for (int idx=0; idx < N2; idx++) {
        int k = idx%ijnum;              
        int icoeff = (idx-k)/ijnum;                 
        int i = ai[k];
        int j = aj[k];
        int itype = ti[k];
        int nb = k + ijnum*3*icoeff;
        int ni = 3*i + 3*anum*icoeff + 3*anum*nperdim*(itype-1);
        int nj = 3*j + 3*anum*icoeff + 3*anum*nperdim*(itype-1);
        snad[0+ni] += dblist[0*ijnum+nb];
        snad[1+ni] += dblist[1*ijnum+nb];
        snad[2+ni] += dblist[2*ijnum+nb];
        snad[0+nj] += -dblist[0*ijnum+nb];
        snad[1+nj] += -dblist[1*ijnum+nb];
        snad[2+nj] += -dblist[2*ijnum+nb];        
    }
}
template <typename T>  void cpuKernelComputeSnad2(T *snad, T *dblist, T *blist, int *aii, 
        int *ai, int *aj, int *ti, int anum, int nperdim, int ncoeff, int inum, int ijnum, int N2)
{        
    for (int idx=0; idx < N2; idx++) {
        int k = idx%ijnum;              
        int icoeff = (idx-k)/ijnum;                 
        int i = ai[k];
        int j = aj[k];
        int itype = ti[k];

        T bi = blist[aii[k] + inum*icoeff];
        T bix = dblist[k + ijnum*0 + ijnum*3*icoeff];
        T biy = dblist[k + ijnum*1 + ijnum*3*icoeff];
        T biz = dblist[k + ijnum*2 + ijnum*3*icoeff];

        // diagonal elements of quadratic matrix
        T dbxtmp = bi*bix;
        T dbytmp = bi*biy;
        T dbztmp = bi*biz;

        int ncount = ncoeff + icoeff*ncoeff - (icoeff-1)*icoeff/2;            
        int ni = 3*i + 3*anum*ncount + 3*anum*nperdim*(itype-1);
        int nj = 3*j + 3*anum*ncount + 3*anum*nperdim*(itype-1);
        snad[0+ni] += dbxtmp;
        snad[1+ni] += dbytmp;
        snad[2+ni] += dbztmp;
        snad[0+nj] += -dbxtmp;
        snad[1+nj] += -dbytmp;
        snad[2+nj] += -dbztmp;
        ncount++;
        // upper-triangular elements of quadratic matrix

        for (int jcoeff = icoeff+1; jcoeff < ncoeff; jcoeff++) {
            dbxtmp = bi*dblist[k + ijnum*0 + ijnum*3*jcoeff] + bix*blist[aii[k] + inum*jcoeff]; 
            dbytmp = bi*dblist[k + ijnum*1 + ijnum*3*jcoeff] + biy*blist[aii[k] + inum*jcoeff];
            dbztmp = bi*dblist[k + ijnum*2 + ijnum*3*jcoeff] + biz*blist[aii[k] + inum*jcoeff];
            int ni = 3*i + 3*anum*ncount + 3*anum*nperdim*(itype-1);
            int nj = 3*j + 3*anum*ncount + 3*anum*nperdim*(itype-1);
            snad[0+ni] += dbxtmp;
            snad[1+ni] += dbytmp;
            snad[2+ni] += dbztmp;
            snad[0+nj] += -dbxtmp;
            snad[1+nj] += -dbytmp;
            snad[2+nj] += -dbztmp;
            ncount++;              
        }                                    
    }
}
template <typename T> void cpuComputeSnad(T *snad, T *dblist, T *blist, int *aii, int *ai, int *aj, int *ti, 
       int anum, int nperdim, int ncoeff, int inum, int ijnum, int quadraticflag)
{        
    int N2 = ncoeff*ijnum;
    cpuKernelComputeSnad1(snad, dblist, ai, aj, 
          ti, anum, nperdim, ijnum, N2);

    if (quadraticflag)
        cpuKernelComputeSnad2(snad, dblist, blist, aii, ai, aj, 
          ti, anum, nperdim, ncoeff, inum, ijnum, N2);
}
template void cpuComputeSnad(double*, double*, double*, int*, int*, int*, int*,
        int, int, int, int, int, int);
template void cpuComputeSnad(float*, float*, float*, int*, int*, int*, int*, 
        int, int, int, int, int, int);

template <typename T>  void cpuKernelComputeSnav1(T *snav, T *dblist, T *x, 
        int *ai, int *aj, int *ti, int anum, int nperdim, int ijnum, int N2)
{        
    for (int idx=0; idx < N2; idx++) {
        int k = idx%ijnum;              
        int icoeff = (idx-k)/ijnum;                 
        int i = ai[k];
        int j = aj[k];
        int itype = ti[k];
        int nb = k + ijnum*3*icoeff;
        int ni = i + 6*anum*icoeff + 6*anum*nperdim*(itype-1);
        int nj = j + 6*anum*icoeff + 6*anum*nperdim*(itype-1);
        T xi = x[0+3*i]; 
        T yi = x[1+3*i]; 
        T zi = x[2+3*i]; 
        T xj = x[0+3*j]; 
        T yj = x[1+3*j]; 
        T zj = x[2+3*j];             
        snav[0*anum+ni] +=  dblist[0*ijnum+nb]*xi;
        snav[1*anum+ni] +=  dblist[1*ijnum+nb]*yi;
        snav[2*anum+ni] +=  dblist[2*ijnum+nb]*zi;
        snav[3*anum+ni] +=  dblist[1*ijnum+nb]*zi;
        snav[4*anum+ni] +=  dblist[0*ijnum+nb]*zi;
        snav[5*anum+ni] +=  dblist[0*ijnum+nb]*yi;
        snav[0*anum+nj] += -dblist[0*ijnum+nb]*xj;
        snav[1*anum+nj] += -dblist[1*ijnum+nb]*yj;
        snav[2*anum+nj] += -dblist[2*ijnum+nb]*zj;
        snav[3*anum+nj] += -dblist[1*ijnum+nb]*zj;
        snav[4*anum+nj] += -dblist[0*ijnum+nb]*zj;
        snav[5*anum+nj] += -dblist[0*ijnum+nb]*yj;        
    }
}
template <typename T>  void cpuKernelComputeSnav2(T *snav, T *dblist, T *blist, T *x, 
        int *aii, int *ai, int *aj,int *ti, int anum, int nperdim, int ncoeff, int inum, int ijnum, int N2)
{        
    for (int idx=0; idx < N2; idx++) {
        int k = idx%ijnum;              
        int icoeff = (idx-k)/ijnum;                 
        int i = ai[k];        
        int j = aj[k];
        int itype = ti[k];

        T xi = x[0+3*i]; 
        T yi = x[1+3*i]; 
        T zi = x[2+3*i]; 
        T xj = x[0+3*j]; 
        T yj = x[1+3*j]; 
        T zj = x[2+3*j];             

        T bi = blist[aii[k] + inum*icoeff];
        T bix = dblist[k + ijnum*0 + ijnum*3*icoeff];
        T biy = dblist[k + ijnum*1 + ijnum*3*icoeff];
        T biz = dblist[k + ijnum*2 + ijnum*3*icoeff];

        // diagonal elements of quadratic matrix
        T dbxtmp = bi*bix;
        T dbytmp = bi*biy;
        T dbztmp = bi*biz;

        int ncount = ncoeff + icoeff*ncoeff - (icoeff-1)*icoeff/2;            
        int ni = i + 6*anum*ncount + 6*anum*nperdim*(itype-1);
        int nj = j + 6*anum*ncount + 6*anum*nperdim*(itype-1);
        snav[0*anum+ni] +=  dbxtmp*xi;
        snav[1*anum+ni] +=  dbytmp*yi;
        snav[2*anum+ni] +=  dbztmp*zi;
        snav[3*anum+ni] +=  dbytmp*zi;
        snav[4*anum+ni] +=  dbxtmp*zi;
        snav[5*anum+ni] +=  dbxtmp*yi;
        snav[0*anum+nj] += -dbxtmp*xj;
        snav[1*anum+nj] += -dbytmp*yj;
        snav[2*anum+nj] += -dbztmp*zj;
        snav[3*anum+nj] += -dbytmp*zj;
        snav[4*anum+nj] += -dbxtmp*zj;
        snav[5*anum+nj] += -dbxtmp*yj;
        ncount++;      

        for (int jcoeff = icoeff+1; jcoeff < ncoeff; jcoeff++) {
            dbxtmp = bi*dblist[k + ijnum*0 + ijnum*3*jcoeff] + bix*blist[aii[k] + inum*jcoeff]; 
            dbytmp = bi*dblist[k + ijnum*1 + ijnum*3*jcoeff] + biy*blist[aii[k] + inum*jcoeff];
            dbztmp = bi*dblist[k + ijnum*2 + ijnum*3*jcoeff] + biz*blist[aii[k] + inum*jcoeff];
            int ni = i + 6*anum*ncount + 6*anum*nperdim*(itype-1);
            int nj = j + 6*anum*ncount + 6*anum*nperdim*(itype-1);
            snav[0*anum+ni] +=  dbxtmp*xi;
            snav[1*anum+ni] +=  dbytmp*yi;
            snav[2*anum+ni] +=  dbztmp*zi;
            snav[3*anum+ni] +=  dbytmp*zi;
            snav[4*anum+ni] +=  dbxtmp*zi;
            snav[5*anum+ni] +=  dbxtmp*yi;
            snav[0*anum+nj] += -dbxtmp*xj;
            snav[1*anum+nj] += -dbytmp*yj;
            snav[2*anum+nj] += -dbztmp*zj;
            snav[3*anum+nj] += -dbytmp*zj;
            snav[4*anum+nj] += -dbxtmp*zj;
            snav[5*anum+nj] += -dbxtmp*yj;
            ncount++;      
        }                                            
    }
}

template <typename T> void cpuComputeSnav(T *snav, T *dblist, T *blist, T *x, int *aii, int *ai, int *aj,
        int *ti, int anum, int nperdim, int ncoeff, int inum, int ijnum, int quadraticflag)
{        
    int N2 = ncoeff*ijnum;
    cpuKernelComputeSnav1(snav, dblist, x, ai, aj, 
          ti, anum, nperdim, ijnum, N2);

    if (quadraticflag)
        cpuKernelComputeSnav2(snav, dblist, blist, x, aii, ai, aj, 
          ti, anum, nperdim, ncoeff, inum, ijnum, N2);
}
template void cpuComputeSnav(double*, double*, double*, double*, int*, int*, int*, int*, 
        int, int, int, int, int, int);
template void cpuComputeSnav(float*, float*, float*, float*, int*, int*, int*, int*, 
        int, int, int, int, int, int);

template <typename T> void cpuSnapTallyEnergyFull2(T *eatom, T *bispectrum, T *coeffelem, int *ilist, 
        int *map, int *type, int inum, int ncoeff, int ncoeffall, int quadraticflag)
{  
    for (int ii=0; ii<inum; ii++) {
        int i = ilist[ii]; // index of atom i
        const int itype = type[i]; // element type of atom i
        const int ielem = map[itype];  // index of that element type
        T* coeffi = &coeffelem[ielem*ncoeffall]; // coefficient for that particular element
        
        T evdwl = coeffi[0];
        for (int icoeff = 0; icoeff < ncoeff; icoeff++)
            evdwl += coeffi[icoeff+1]*bispectrum[ii + inum*icoeff];
                  
        if (quadraticflag) {
            int k = ncoeff+1;
            for (int icoeff = 0; icoeff < ncoeff; icoeff++) {
                T bveci = bispectrum[ii + inum*icoeff];
                evdwl += 0.5*coeffi[k++]*bveci*bveci;
                for (int jcoeff = icoeff+1; jcoeff < ncoeff; jcoeff++) {
                    T bvecj = bispectrum[ii + inum*jcoeff];
                    evdwl += coeffi[k++]*bveci*bvecj;
                }
            }
        }                
        eatom[ii] += evdwl;
    }
}
template void cpuSnapTallyEnergyFull2(double*, double*, double*, int*, int*, int*,
        int, int, int, int);
template void cpuSnapTallyEnergyFull2(float*, float*, float*, int*, int*, int*,
        int, int, int, int);

template <typename T> void cpuSnapTallyBispectrum(T *bi, T *bispectrum, int *ilist, 
        int *type, int inum, int ncoeff, int nperdim, int ntype, int quadraticflag)
{      
    cpuArraySetValue(bi, (T) 0.0, nperdim*ntype);
    
    int N2 = inum*ncoeff;
    for (int idx=0; idx<N2; idx++) {
        int ii = idx%inum;
        int icoeff = (idx-ii)/inum;
        int i = ilist[ii]; // index of atom i
        int itype = type[i]; // element type of atom i        
        int n = nperdim*(itype-1);
        bi[icoeff+n] += bispectrum[ii + inum*icoeff];                  
        if (quadraticflag) {
            int k = n+ncoeff + ncoeff*(ncoeff+1)/2 - (ncoeff-icoeff)*(ncoeff-icoeff+1)/2;
            T bveci = bispectrum[ii + inum*icoeff];
            bi[k] += 0.5*bveci*bveci;
            for (int jcoeff = icoeff+1; jcoeff < ncoeff; jcoeff++) {
                T bvecj = bispectrum[ii + inum*jcoeff];
                bi[k++] += bveci*bvecj;
            }
        }                
    }
}
template void cpuSnapTallyBispectrum(double*, double*, int*, int*, int, int, int, int, int);
template void cpuSnapTallyBispectrum(float*, float*, int*, int*, int, int, int, int, int);

template <typename T> void cpuSnapTallyBispectrumDeriv(T *db, T *bispectrum, T *dbdr, int *aii, 
        int *ai, int *aj, int *ti, int inum, int ijnum, int ncoeff, int nperdim, int ntype, int quadraticflag)
{   
    cpuArraySetValue(db, (T) 0.0, inum*3*nperdim*ntype);
        
    int N2 = ijnum*ncoeff;
    for (int idx=0; idx<N2; idx++) {
        int ij = idx%ijnum;
        int icoeff = (idx-ij)/ijnum;        
        int ii = aii[ij]; // index of atom i
        int i = ai[ij]; // index of atom i
        int j = aj[ij]; // index of atom i
        int itype = ti[ij]; // element type of atom i       
        int n = nperdim*(itype-1);        
        int nii = inum*3*(icoeff + n);  
        int nij = ijnum*3*icoeff;
        
        T bix = dbdr[ij + ijnum*0 + nij];
        T biy = dbdr[ij + ijnum*1 + nij];
        T biz = dbdr[ij + ijnum*2 + nij];        
        db[0 + 3*i + nii] += bix; 
        db[1 + 3*i + nii] += biy;
        db[2 + 3*i + nii] += biz;
        db[0 + 3*j + nii] -= bix;
        db[1 + 3*j + nii] -= biy;
        db[2 + 3*j + nii] -= biz;
        
        if (quadraticflag) {
            T bi = bispectrum[ii + inum*icoeff];
            T dbxtmp = bi*bix;
            T dbytmp = bi*biy;
            T dbztmp = bi*biz;
            int k = ncoeff + ncoeff*(ncoeff+1)/2 - (ncoeff-icoeff)*(ncoeff-icoeff+1)/2;
            nii = inum*3*(k + n);  
            db[0 + 3*i + nii] += dbxtmp;
            db[1 + 3*i + nii] += dbytmp;
            db[2 + 3*i + nii] += dbztmp;
            db[0 + 3*j + nii] -= dbxtmp;
            db[1 + 3*j + nii] -= dbytmp;
            db[2 + 3*j + nii] -= dbztmp;          
         
            for (int jcoeff = icoeff+1; jcoeff < ncoeff; jcoeff++) {
                int nj = ijnum*3*jcoeff;
                T bjx = dbdr[ij + ijnum*0 + nj];
                T bjy = dbdr[ij + ijnum*1 + nj];
                T bjz = dbdr[ij + ijnum*2 + nj];        
                T bj = bispectrum[ii + inum*jcoeff];                
                dbxtmp = bi*bjx + bix*bj;
                dbytmp = bi*bjy + biy*bj;
                dbztmp = bi*bjz + biz*bj;

                k += 1;
                nii = inum*3*(k + n);  
                db[0 + 3*i + nii] += dbxtmp;
                db[1 + 3*i + nii] += dbytmp;
                db[2 + 3*i + nii] += dbztmp;
                db[0 + 3*j + nii] -= dbxtmp;
                db[1 + 3*j + nii] -= dbytmp;
                db[2 + 3*j + nii] -= dbztmp;                          
//                 db[i + inum*0 + nii] += dbxtmp;
//                 db[i + inum*1 + nii] += dbytmp;
//                 db[i + inum*2 + nii] += dbztmp;
//                 db[j + inum*0 + nii] -= dbxtmp;
//                 db[j + inum*1 + nii] -= dbytmp;
//                 db[j + inum*2 + nii] -= dbztmp;                                          
            }            
        }        
    }
}
template void cpuSnapTallyBispectrumDeriv(double*, double*, double*, int*, int*, int*, int*, 
        int, int, int, int, int, int);
template void cpuSnapTallyBispectrumDeriv(float*, float*, float*, int*, int*, int*, int*, 
        int, int, int, int, int, int);

template <typename T> void cpuSnapTallyBispectrumVirial(T *bv, T *bispectrum, T *dbdr, T *rij, int *aii, 
        int *ai, int *aj, int *ti, int inum, int ijnum, int ncoeff, int nperdim, int ntype, int quadraticflag)
{      
    cpuArraySetValue(bv, (T) 0.0, 6*nperdim*ntype);
    
    int N2 = ijnum*ncoeff;
    for (int idx=0; idx<N2; idx++) {
        int ij = idx%ijnum;
        int icoeff = (idx-ij)/ijnum;        
        int ii = aii[ij]; // index of atom i
        int i = ai[ij]; // index of atom i
        int j = aj[ij]; // index of atom i
        int itype = ti[ij]; // element type of atom i       
        int n = nperdim*(itype-1);        
        int nii = 6*(icoeff + n);  
        int nij = ijnum*3*icoeff;
        
        T factor = 1.0;
        T dx = -rij[0+3*ij];
        T dy = -rij[1+3*ij];
        T dz = -rij[2+3*ij];                    
        T bix = dbdr[ij + ijnum*0 + nij];
        T biy = dbdr[ij + ijnum*1 + nij];
        T biz = dbdr[ij + ijnum*2 + nij];                
        T v0 = factor*dx*bix;
        T v1 = factor*dy*biy;
        T v2 = factor*dz*biz;
        T v3 = factor*dx*biy;
        T v4 = factor*dx*biz;
        T v5 = factor*dy*biz;        
        
        bv[0 + nii] += v0;
        bv[1 + nii] += v1;
        bv[2 + nii] += v2;
        bv[3 + nii] += v3;
        bv[4 + nii] += v4;
        bv[5 + nii] += v5;        
//         bv[0 + nii] += v0;
//         bv[1 + nii] += v1;
//         bv[2 + nii] += v2;
//         bv[3 + nii] += v3;
//         bv[4 + nii] += v4;
//         bv[5 + nii] += v5;                       
        if (quadraticflag) {
            T bi = bispectrum[ii + inum*icoeff];
            T dbxtmp = bi*bix;
            T dbytmp = bi*biy;
            T dbztmp = bi*biz;
            int k = ncoeff + ncoeff*(ncoeff+1)/2 - (ncoeff-icoeff)*(ncoeff-icoeff+1)/2;
            nii = 6*(k + n);  
            T v0 = factor*dx*dbxtmp;
            T v1 = factor*dy*dbytmp;
            T v2 = factor*dz*dbztmp;
            T v3 = factor*dx*dbytmp;
            T v4 = factor*dx*dbztmp;
            T v5 = factor*dy*dbztmp;        
            bv[0 + nii] += v0;
            bv[1 + nii] += v1;
            bv[2 + nii] += v2;
            bv[3 + nii] += v3;
            bv[4 + nii] += v4;
            bv[5 + nii] += v5;        
//             bv[0 + nii] += v0;
//             bv[1 + nii] += v1;
//             bv[2 + nii] += v2;
//             bv[3 + nii] += v3;
//             bv[4 + nii] += v4;
//             bv[5 + nii] += v5;                                   
            for (int jcoeff = icoeff+1; jcoeff < ncoeff; jcoeff++) {
                int nj = ijnum*3*jcoeff;
                T bjx = dbdr[ij + ijnum*0 + nj];
                T bjy = dbdr[ij + ijnum*1 + nj];
                T bjz = dbdr[ij + ijnum*2 + nj];        
                T bj = bispectrum[ii + inum*jcoeff];                
                dbxtmp = bi*bjx + bix*bj;
                dbytmp = bi*bjy + biy*bj;
                dbztmp = bi*bjz + biz*bj;

                k += 1;                
                nii = 6*(k + n);  
                T v0 = factor*dx*dbxtmp;
                T v1 = factor*dy*dbytmp;
                T v2 = factor*dz*dbztmp;
                T v3 = factor*dx*dbytmp;
                T v4 = factor*dx*dbztmp;
                T v5 = factor*dy*dbztmp;        
                bv[0 + nii] += v0;
                bv[1 + nii] += v1;
                bv[2 + nii] += v2;
                bv[3 + nii] += v3;
                bv[4 + nii] += v4;
                bv[5 + nii] += v5;        
//                 bv[0 + nii] += v0;
//                 bv[1 + nii] += v1;
//                 bv[2 + nii] += v2;
//                 bv[3 + nii] += v3;
//                 bv[4 + nii] += v4;
//                 bv[5 + nii] += v5;                                                   
            }            
        }        
    }
}
template void cpuSnapTallyBispectrumVirial(double*, double*, double*, double*, int*, int*, int*, int*, 
        int, int, int, int, int, int);
template void cpuSnapTallyBispectrumVirial(float*, float*, float*, float*, int*, int*, int*, int*, 
        int, int, int, int, int, int);

template <typename T> void cpuSnapTallyForceFull2(T *fatom, T *fij, int *ai, int *aj, int *alist, int ijnum)
{ 
    for(int k = 0; k < ijnum; k++) {
        int i = ai[k];        
        int j = aj[k];        
        T fx = fij[k];
        T fy = fij[k+ijnum];
        T fz = fij[k+ijnum*2];    
        fatom[0+3*i] += fx;
        fatom[1+3*i] += fy;
        fatom[2+3*i] += fz;
        fatom[0+3*j] -= fx;
        fatom[1+3*j] -= fy;
        fatom[2+3*j] -= fz;
    }
};
template void cpuSnapTallyForceFull2(double*, double*, int*, int*, int*, int);
template void cpuSnapTallyForceFull2(float*, float*, int*, int*, int*, int);

template <typename T> void cpuSnapTallyVirialFull2(T *vatom, T *fij, T *rij, int *ai, int *aj,
        int inum, int ijnum)
{ 
    for(int k = 0; k < ijnum; k++) {
        int i = ai[k];        
        int j = aj[k];        
        T factor = 0.5;        
        T dx = -rij[0+3*k];
        T dy = -rij[1+3*k];
        T dz = -rij[2+3*k];    
        T fx = fij[k];
        T fy = fij[k+ijnum];
        T fz = fij[k+ijnum*2];    
        T v0 = factor*dx*fx;
        T v1 = factor*dy*fy;
        T v2 = factor*dz*fz;
        T v3 = factor*dx*fy;
        T v4 = factor*dx*fz;
        T v5 = factor*dy*fz;        
        vatom[0*inum+i] += v0;
        vatom[1*inum+i] += v1;
        vatom[2*inum+i] += v2;
        vatom[3*inum+i] += v3;
        vatom[4*inum+i] += v4;
        vatom[5*inum+i] += v5;        
        vatom[0*inum+j] += v0;
        vatom[1*inum+j] += v1;
        vatom[2*inum+j] += v2;
        vatom[3*inum+j] += v3;
        vatom[4*inum+j] += v4;
        vatom[5*inum+j] += v5;        
    }
};
template void cpuSnapTallyVirialFull2(double*, double*, double*, int*, int*, int, int);
template void cpuSnapTallyVirialFull2(float*, float*, float*, int*, int*, int, int);


// template <typename T>  void cpuKernelComputeUij(T *slist, T *Sr, T *Si, T *rootpqarray, 
//         T *rij, T *wjelem, T *radelem, T rmin0, T rfac0, T rcutfac, int *idxu_block,  int *ai, int *aj, 
//         int *ti, int *tj, int twojmax, int idxu_max, int ijnum, int switch_flag)
// { 
//    
//   for(int ij=0; ij<ijnum; ij++) {        
//     T x = rij[ij*3+0];
//     T y = rij[ij*3+1];
//     T z = rij[ij*3+2];
//     T rsq = x * x + y * y + z * z;
//     T r = sqrt(rsq);
// 
//     T rcutij = (radelem[ti[ij]]+radelem[tj[ij]])*rcutfac; //(radelem[type[ii]]+radelem[type[jj]])*rcutfac;
//     T theta0 = (r - rmin0) * rfac0 * M_PI / (rcutij - rmin0);
//     //    theta0 = (r - rmin0) * rscale0;
//     T z0 = r / tan(theta0);    
//             
//     T sfac = 0.0, dsfac = 0.0;        
//     if (switch_flag == 0) {
//         sfac = 1.0;
//         dsfac = 0.0;
//     }
//     else if (switch_flag == 1) {
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
//     slist[ij] = sfac;
//     slist[ij+ijnum] = dsfac;
// 
//     T r0inv;
//     T a_r, b_r, a_i, b_i;
//     T rootpq;
//     int jdim = twojmax + 1;
//   
//     r0inv = 1.0 / sqrt(r * r + z0 * z0);
//     a_r = r0inv * z0;
//     a_i = -r0inv * z;
//     b_r = r0inv * y;
//     b_i = -r0inv * x;
// 
//     Sr[ij+0*ijnum] = 1.0;
//     Si[ij+0*ijnum] = 0.0;
// 
//     for (int j = 1; j <= twojmax; j++) {
//         int jju = idxu_block[j];
//         int jjup = idxu_block[j-1];
//         
//         // fill in left side of matrix layer from previous layer
//         for (int mb = 0; 2*mb <= j; mb++) {
//             Sr[ij+jju*ijnum] = 0.0;
//             Si[ij+jju*ijnum] = 0.0;
//             for (int ma = 0; ma < j; ma++) {
//                 rootpq = rootpqarray[(j - ma)*jdim + (j - mb)];
//                 int njju = ij+jju*ijnum;
//                 int njjup = ij+jjup*ijnum;
//                 T u_r = Sr[njjup];
//                 T u_i = Si[njjup];
//                 Sr[njju] += rootpq * (a_r * u_r + a_i * u_i);
//                 Si[njju] += rootpq * (a_r * u_i - a_i * u_r);
// 
//                 rootpq = rootpqarray[(ma + 1)*jdim + (j - mb)];
//                 Sr[njju+1] = -rootpq * (b_r * u_r + b_i * u_i);
//                 Si[njju+1] = -rootpq * (b_r * u_i - b_i * u_r);
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
//                 int njju = ij+jju*ijnum;
//                 int njjup = ij+jjup*ijnum;
//                 if (mapar == 1) {
//                     Sr[njjup] = Sr[njju];
//                     Si[njjup] = -Si[njju];
//                 } else {
//                     Sr[njjup] = -Sr[njju];
//                     Si[njjup] = Si[njju];
//                 }
//                 mapar = -mapar;
//                 jju++;
//                 jjup--;
//             }
//             mbpar = -mbpar;
//         }
//     }        
//     
//   }
// };
// template <typename T> void cpuComputeUij(T *slist, T *Sr, T *Si, T *rootpqarray, 
//         T *rij, T *wjelem, T *radelem, T rmin0, T rfac0, T rcutfac, int *idxu_block,  int *ai, int *aj, 
//         int *ti, int *tj, int twojmax, int idxu_max, int ijnum, int switch_flag)                
// {
//     int blockDim = 256;
//     int gridDim = (ijnum + blockDim - 1) / blockDim;
//     gridDim = (gridDim>1024)? 1024 : gridDim;
//     cpuKernelComputeUij(slist, Sr, Si, rootpqarray, rij, wjelem, radelem, 
//             rmin0, rfac0, rcutfac, idxu_block, ai, aj, ti, tj, twojmax, idxu_max, ijnum, switch_flag);
// }
// template void cpuComputeUij(double*, double*, double*, double*, double*, double*, double*, double, 
//         double, double, int*, int*, int*, int*, int*, int, int, int, int);
// template void cpuComputeUij(float*, float*, float*, float*, float*, float*, float*, float,  float,
//         float, int*, int*, int*, int*, int*, int, int, int, int);
// 
// template <typename T>  void cpuKernelComputeDuijdrj(T *dSr, T *dSi, T *Sr, 
//      T *Si, T *sfac, T *dsfac, T *rootpqarray, T* rij, T *radelem, T rmin0, T rfac0, T rcutfac, 
//      int *idxu_block, int *ti, int *tj, int twojmax, int idxu_max, int ijnum, int switch_flag)
// {
//   int ij = threadIdx.x + blockIdx.x * blockDim.x;
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
//     T dz0dr = z0 / r - (r*rscale0) * (rsq + z0 * z0) / rsq;
// 
//     // dulist [ijnum*3*idxu_max], ulist [ijnum*idxu_max]
//     T r0inv;
//     T a_r, a_i, b_r, b_i;
//     T da_r[3], da_i[3], db_r[3], db_i[3];
//     T dz0[3], dr0inv[3], dr0invdr;
//     T rootpq;
// 
//     T rinv = 1.0 / r;
//     T ux = x * rinv;
//     T uy = y * rinv;
//     T uz = z * rinv;
// 
//     int jdim = twojmax + 1;
// 
//     r0inv = 1.0 / sqrt(r * r + z0 * z0);
//     a_r = z0 * r0inv;
//     a_i = -z * r0inv;
//     b_r = y * r0inv;
//     b_i = -x * r0inv;
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
// 
//     da_i[2] += -r0inv;
// 
//     for (int k = 0; k < 3; k++) {
//         db_r[k] = y * dr0inv[k];
//         db_i[k] = -x * dr0inv[k];
//     }
// 
//     db_i[0] += -r0inv;
//     db_r[1] += r0inv;
//     
//     int n3 = ijnum*3; 
//     int j0 = ij + ijnum*0;
//     int j1 = ij + ijnum*1;
//     int j2 = ij + ijnum*2;
//     dSr[j0] = 0.0;
//     dSr[j1] = 0.0;
//     dSr[j2] = 0.0;
//     dSi[j0] = 0.0;
//     dSi[j1] = 0.0;
//     dSi[j2] = 0.0;
// 
//     for (int j = 1; j <= twojmax; j++) {
//         int jju = idxu_block[j];
//         int jjup = idxu_block[j-1];
//         for (int mb = 0; 2*mb <= j; mb++) {
//           dSr[j0+n3*jju] = 0.0;
//           dSr[j1+n3*jju] = 0.0;
//           dSr[j2+n3*jju] = 0.0;
//           dSi[j0+n3*jju] = 0.0;
//           dSi[j1+n3*jju] = 0.0;
//           dSi[j2+n3*jju] = 0.0;
// 
//           for (int ma = 0; ma < j; ma++) {
//             rootpq = rootpqarray[(j - ma)*jdim + (j - mb)];
//             int njju = ij+jju*ijnum;
//             int njjup = ij+jjup*ijnum;
//             T u_r = Sr[njjup];
//             T u_i = Si[njjup];
//             for (int k = 0; k < 3; k++) {
//               dSr[ij + ijnum*k + n3*jju] +=
//                 rootpq * (da_r[k] * u_r + da_i[k] * u_i +
//                           a_r * dSr[ij + ijnum*k + n3*jjup] + a_i * dSi[ij + ijnum*k + n3*jjup]);
//               dSi[ij + ijnum*k + n3*jju] +=
//                 rootpq * (da_r[k] * u_i - da_i[k] * u_r +
//                           a_r * dSi[ij + ijnum*k + n3*jjup] - a_i * dSr[ij + ijnum*k + n3*jjup]);
//             }
// 
//             rootpq = rootpqarray[(ma + 1)*jdim + (j - mb)];
//             for (int k = 0; k < 3; k++) {
//               dSr[ij + ijnum*k + n3*(jju+1)] =
//                 -rootpq * (db_r[k] * u_r + db_i[k] * u_i +
//                            b_r * dSr[ij + ijnum*k + n3*jjup] + b_i * dSi[ij + ijnum*k + n3*jjup]);
//               dSi[ij + ijnum*k + n3*(jju+1)] =
//                 -rootpq * (db_r[k] * u_i - db_i[k] * u_r +
//                            b_r * dSi[ij + ijnum*k + n3*jjup] - b_i * dSr[ij + ijnum*k + n3*jjup]);
//             }
//             jju++;
//             jjup++;
//           }
//           jju++;
//         }
// 
//         // copy left side to right side with inversion symmetry VMK 4.4(2)
//         // u[ma-j][mb-j] = (-1)^(ma-mb)*Conj([u[ma][mb])
// 
//         jju = idxu_block[j];
//         jjup = jju+(j+1)*(j+1)-1;
//         int mbpar = 1;
//         for (int mb = 0; 2*mb <= j; mb++) {
//           int mapar = mbpar;
//           for (int ma = 0; ma <= j; ma++) {
//             if (mapar == 1) {
//               for (int k = 0; k < 3; k++) {
//                 dSr[ij + ijnum*k + n3*jjup] =  dSr[ij + ijnum*k + n3*jju];
//                 dSi[ij + ijnum*k + n3*jjup] = -dSi[ij + ijnum*k + n3*jju];
//               }
//             } else {
//               for (int k = 0; k < 3; k++) {
//                 dSr[ij + ijnum*k + n3*jjup] = -dSr[ij + ijnum*k + n3*jju];
//                 dSi[ij + ijnum*k + n3*jjup] =  dSi[ij + ijnum*k + n3*jju];
//               }
//             }
//             mapar = -mapar;
//             jju++;
//             jjup--;
//           }
//           mbpar = -mbpar;
//         }
//     }
// 
//     for (int j = 0; j <= twojmax; j++) {
//     int jju = idxu_block[j];
//     int njju = ij+jju*ijnum;
//     T u_r = Sr[njju];
//     T u_i = Si[njju];
//     for (int mb = 0; 2*mb <= j; mb++)
//       for (int ma = 0; ma <= j; ma++) {
//         dSr[j0 + n3*jju] = dsfac[ij] * u_r * ux +
//                                   sfac[ij] * dSr[j0 + n3*jju];
//         dSi[j0 + n3*jju] = dsfac[ij] * u_i * ux +
//                                   sfac[ij] * dSi[j0 + n3*jju];
//         dSr[j1 + n3*jju] = dsfac[ij] * u_r * uy +
//                                   sfac[ij] * dSr[j1 + n3*jju];
//         dSi[j1 + n3*jju] = dsfac[ij] * u_i * uy +
//                                   sfac[ij] * dSi[j1 + n3*jju];
//         dSr[j2 + n3*jju] = dsfac[ij] * u_r * uz +
//                                   sfac[ij] * dSr[j2 + n3*jju];
//         dSi[j2 + n3*jju] = dsfac[ij] * u_i * uz +
//                                   sfac[ij] * dSi[j2 + n3*jju];
//         jju++;
//       }
//     }        
//     
//   }  
// }
// template <typename T> void cpuComputeDuijdrj(T *dSr, T *dSi, T *Sr, T *Si, 
//     T *sfac, T *dsfac, T *rootpqarray, T* rij, T *radelem, T rmin0, T rfac0, T rcutfac, 
//      int *idxu_block, int *ti, int *tj, int twojmax, int idxu_max, int ijnum, int switch_flag)
// {                
//     int blockDim = 256;
//     int gridDim = (ijnum + blockDim - 1) / blockDim;
//     gridDim = (gridDim>1024)? 1024 : gridDim;    
//     cpuKernelComputeDuijdrj(dSr, dSi, Sr, Si, 
//           sfac, dsfac, rootpqarray, rij, radelem, rmin0, rfac0, rcutfac, idxu_block, 
//           ti, tj, twojmax, idxu_max, ijnum, switch_flag);
// }
// template void cpuComputeDuijdrj(double*, double*, double*, double*, double*, double*, 
//         double*, double*, double*, double, double, double, int*, int*, int*, int, int, int, int);
// template void cpuComputeDuijdrj(float*, float*, float*, float*, float*, float*, float*, 
//         float*, float*, float, float, float, int*, int*, int*, int, int, int, int);
// 

#endif

