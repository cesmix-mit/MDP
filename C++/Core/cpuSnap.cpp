/***************************************************************************                               
                    Molecular Dynamics Potentials (MDP)
                           CESMIX-MIT Project  
 
 Contributing authors: Ngoc-Cuong Nguyen (cuongng@mit.edu, exapde@gmail.com)
 ***************************************************************************/

#ifndef MDP_CPUSNAP
#define MDP_CPUSNAP

#define MIN(a,b) ((a) < (b) ? (a) : (b))
#define MAX(a,b) ((a) > (b) ? (a) : (b))

void cpuBuildIndexList(int *idx_max, int *idxz, int *idxz_block, int *idxb, int *idxb_block, int *idxu_block, 
        int *idxcg_block, int twojmax)
{
  // index list for cglist

  int jdim = twojmax + 1;

  int idxcg_count = 0;
  for(int j1 = 0; j1 <= twojmax; j1++)
    for(int j2 = 0; j2 <= j1; j2++)
      for(int j = j1 - j2; j <= MIN(twojmax, j1 + j2); j += 2) {
        idxcg_block[j + j2*jdim + j1*jdim*jdim] = idxcg_count;  
        for (int m1 = 0; m1 <= j1; m1++)
          for (int m2 = 0; m2 <= j2; m2++)
            idxcg_count++;
      }
  idx_max[0] = idxcg_count;
          
  int idxu_count = 0;

  for(int j = 0; j <= twojmax; j++) {
    idxu_block[j] = idxu_count;
    for(int mb = 0; mb <= j; mb++)
      for(int ma = 0; ma <= j; ma++)
        idxu_count++;
  }
  //idxu_max = idxu_count;
  idx_max[1] = idxu_count;
  
  // index list for beta and B

  int idxb_count = 0;
  for(int j1 = 0; j1 <= twojmax; j1++)
    for(int j2 = 0; j2 <= j1; j2++)
      for(int j = j1 - j2; j <= MIN(twojmax, j1 + j2); j += 2)
        if (j >= j1) idxb_count++;

  int idxb_max = idxb_count;
  idx_max[2] = idxb_max;

  idxb_count = 0;
  for(int j1 = 0; j1 <= twojmax; j1++)
    for(int j2 = 0; j2 <= j1; j2++)
      for(int j = j1 - j2; j <= MIN(twojmax, j1 + j2); j += 2)
        if (j >= j1) {
          idxb[idxb_count*3 + 0] = j1;
          idxb[idxb_count*3 + 1] = j2;
          idxb[idxb_count*3 + 2] = j;  
          idxb_count++;
        }
  
  idxb_count = 0;
  for(int j1 = 0; j1 <= twojmax; j1++)
    for(int j2 = 0; j2 <= j1; j2++)
      for(int j = j1 - j2; j <= MIN(twojmax, j1 + j2); j += 2) {
        if (j >= j1) {
          idxb_block[j + j2*jdim + j1*jdim*jdim] = idxb_count;    
          idxb_count++;
        }
      }

  // index list for zlist

  int idxz_count = 0;

  for(int j1 = 0; j1 <= twojmax; j1++)
    for(int j2 = 0; j2 <= j1; j2++)
      for(int j = j1 - j2; j <= MIN(twojmax, j1 + j2); j += 2)
        for (int mb = 0; 2*mb <= j; mb++)
          for (int ma = 0; ma <= j; ma++)
            idxz_count++;

  int idxz_max = idxz_count;
  idx_max[3] = idxz_max;

  idxz_count = 0;
  for(int j1 = 0; j1 <= twojmax; j1++)
    for(int j2 = 0; j2 <= j1; j2++)
      for(int j = j1 - j2; j <= MIN(twojmax, j1 + j2); j += 2) {
        idxz_block[j + j2*jdim + j1*jdim*jdim] = idxz_count;    

        for (int mb = 0; 2*mb <= j; mb++)
          for (int ma = 0; ma <= j; ma++) {

            idxz[idxz_count*10 + 0] = j1;
            idxz[idxz_count*10 + 1] = j2;
            idxz[idxz_count*10 + 2] = j;
            idxz[idxz_count*10 + 3] = MAX(0, (2 * ma - j - j2 + j1) / 2);
            idxz[idxz_count*10 + 4] = (2 * ma - j - (2 * idxz[idxz_count*10 + 3] - j1) + j2) / 2;
            idxz[idxz_count*10 + 5] = MIN(j1, (2 * ma - j + j2 + j1) / 2) - idxz[idxz_count*10 + 3] + 1;
            idxz[idxz_count*10 + 6] = MAX(0, (2 * mb - j - j2 + j1) / 2);
            idxz[idxz_count*10 + 7] = (2 * mb - j - (2 * idxz[idxz_count*10 + 6] - j1) + j2) / 2;
            idxz[idxz_count*10 + 8] = MIN(j1, (2 * mb - j + j2 + j1) / 2) - idxz[idxz_count*10 + 6] + 1;
            // apply to z(j1,j2,j,ma,mb) to unique element of y(j)

            const int jju = idxu_block[j] + (j+1)*mb + ma;
            idxz[idxz_count*10 + 9] = jju;
              
            idxz_count++;
          }
      }
};

template <typename T> void cpuInitRootpqArray(T *rootpqarray, int twojmax)
{
  int jdim = twojmax + 1;
  for (int p = 1; p <= twojmax; p++)
    for (int q = 1; q <= twojmax; q++)
      rootpqarray[p*jdim + q] = sqrt(static_cast<T>(p)/q);
};
template void cpuInitRootpqArray(double*, int);
template void cpuInitRootpqArray(float*, int);

template <typename T> T cpuDeltacg(T *factorial, int j1, int j2, int j)
{
  T sfaccg = factorial[(j1 + j2 + j) / 2 + 1];
  return sqrt(factorial[(j1 + j2 - j) / 2] *
              factorial[(j1 - j2 + j) / 2] *
              factorial[(-j1 + j2 + j) / 2] / sfaccg);
};
template double cpuDeltacg(double*, int, int, int);
template float cpuDeltacg(float*, int, int, int);

template <typename T> void cpuInitClebschGordan(T *cglist, T *factorial, int twojmax)
{
  T sum,dcg,sfaccg;
  int m, aa2, bb2, cc2;
  int ifac;

  int idxcg_count = 0;
  for(int j1 = 0; j1 <= twojmax; j1++)
    for(int j2 = 0; j2 <= j1; j2++)
      for(int j = j1 - j2; j <= MIN(twojmax, j1 + j2); j += 2) {
        for (int m1 = 0; m1 <= j1; m1++) {
          aa2 = 2 * m1 - j1;

          for (int m2 = 0; m2 <= j2; m2++) {

            bb2 = 2 * m2 - j2;
            m = (aa2 + bb2 + j) / 2;

            if(m < 0 || m > j) {
              cglist[idxcg_count] = 0.0;
              idxcg_count++;
              continue;
            }

            sum = 0.0;

            for (int z = MAX(0, MAX(-(j - j2 + aa2)
                                    / 2, -(j - j1 - bb2) / 2));
                 z <= MIN((j1 + j2 - j) / 2,
                          MIN((j1 - aa2) / 2, (j2 + bb2) / 2));
                 z++) {
              ifac = z % 2 ? -1 : 1;
              sum += ifac /
                (factorial[z] *
                 factorial[(j1 + j2 - j) / 2 - z] *
                 factorial[(j1 - aa2) / 2 - z] *
                 factorial[(j2 + bb2) / 2 - z] *
                 factorial[(j - j2 + aa2) / 2 + z] *
                 factorial[(j - j1 - bb2) / 2 + z]);
            }

            cc2 = 2 * m - j;
            dcg = cpuDeltacg(factorial, j1, j2, j);
            sfaccg = sqrt(factorial[(j1 + aa2) / 2] *
                          factorial[(j1 - aa2) / 2] *
                          factorial[(j2 + bb2) / 2] *
                          factorial[(j2 - bb2) / 2] *
                          factorial[(j  + cc2) / 2] *
                          factorial[(j  - cc2) / 2] *
                          (j + 1));

            cglist[idxcg_count] = sum * dcg * sfaccg;
            idxcg_count++;
          }
        }
      }
}
template void cpuInitClebschGordan(double*, double*, int);
template void cpuInitClebschGordan(float*, float*, int);

template <typename T> void cpuInitSna(T *rootpqarray, T *cglist, T *factorial, int *idx_max, int *idxz, 
      int *idxz_block, int *idxb, int *idxb_block, int *idxu_block, int *idxcg_block, int twojmax)
{
    cpuBuildIndexList(idx_max, idxz, idxz_block, idxb, 
            idxb_block, idxu_block, idxcg_block, twojmax);
    
    cpuInitRootpqArray(rootpqarray, twojmax);
    cpuInitClebschGordan(cglist, factorial, twojmax);        
}
template void cpuInitSna(double*, double*, double *, int*, int*, int*, int*, int*, int*, int*, int);
template void cpuInitSna(float*, float*, float*, int*, int*, int*, int*, int*, int*, int*, int);

template <typename T> T cpuComputeSfac(T r, T rcut, T rmin0, int switch_flag)
{
  if (switch_flag == 0) return 1.0;
  if (switch_flag == 1) {
    if(r <= rmin0) return 1.0;
    else if(r > rcut) return 0.0;
    else {
      T rcutfac = M_PI / (rcut - rmin0);
      return 0.5 * (cos((r - rmin0) * rcutfac) + 1.0);
    }
  }
  return 0.0;
};
template double cpuComputeSfac(double, double, double, int);
template float cpuComputeSfac(float, float, float, int);

template <typename T> T cpuComputeDsfac(T r, T rcut, T rmin0, int switch_flag)
{
  if (switch_flag == 0) return 0.0;
  if (switch_flag == 1) {
    if(r <= rmin0) return 0.0;
    else if(r > rcut) return 0.0;
    else {
      T rcutfac = M_PI / (rcut - rmin0);
      return -0.5 * sin((r - rmin0) * rcutfac) * rcutfac;
    }
  }
  return 0.0;
}
template double cpuComputeDsfac(double, double, double, int);
template float cpuComputeDsfac(float, float, float, int);

template <typename T> void cpuZeroUarraytot(T *ulisttot_r, T *ulisttot_i, T wself, int *idxu_block, int *type,
        int *map, int *ai, int wselfall_flag, int chemflag, int idxu_max, int nelements, int twojmax, int inum)
{
  for(int ii = 0; ii < inum; ii++)
  for(int jelem = 0; jelem < nelements; jelem++)
  for (int j = 0; j <= twojmax; j++) {
    int jju = idxu_block[j];
    for (int mb = 0; mb <= j; mb++) {
      for (int ma = 0; ma <= j; ma++) {
        ulisttot_r[ii*nelements*idxu_max+jelem*idxu_max+jju] = 0.0;
        ulisttot_i[ii*nelements*idxu_max+jelem*idxu_max+jju] = 0.0;
                
        int ielem = (chemflag) ? map[type[ai[ii]]]: 0;                
        // utot(j,ma,ma) = wself, sometimes
        if (jelem == ielem || wselfall_flag)
          if (ma==mb)
          ulisttot_r[ii*nelements*idxu_max+jelem*idxu_max+jju] = wself; ///// T check this
        jju++;
      }
    }
  }
};
template void cpuZeroUarraytot(double*, double*, double, int*, int*,
        int*, int*, int, int, int, int, int, int);
template void cpuZeroUarraytot(float*, float*, float, int*, int*,
        int*, int*, int, int, int, int, int, int);

template <typename T> void cpuComputeUarray(T *ulist_r, T *ulist_i, T *rootpqarray, T x, T y, T z,
                         T z0, T r, int *idxu_block, int twojmax)
{
  T r0inv;
  T a_r, b_r, a_i, b_i;
  T rootpq;
  int jdim = twojmax + 1;
  
  // compute Cayley-Klein parameters for unit quaternion

  r0inv = 1.0 / sqrt(r * r + z0 * z0);
  a_r = r0inv * z0;
  a_i = -r0inv * z;
  b_r = r0inv * y;
  b_i = -r0inv * x;

  // VMK Section 4.8.2


  //T* ulist_r = ulist_r_ij[jj];
  //T* ulist_i = ulist_i_ij[jj];

  ulist_r[0] = 1.0;
  ulist_i[0] = 0.0;

  for (int j = 1; j <= twojmax; j++) {
    int jju = idxu_block[j];
    int jjup = idxu_block[j-1];

    // fill in left side of matrix layer from previous layer

    for (int mb = 0; 2*mb <= j; mb++) {
      ulist_r[jju] = 0.0;
      ulist_i[jju] = 0.0;

      for (int ma = 0; ma < j; ma++) {
        rootpq = rootpqarray[(j - ma)*jdim + (j - mb)];
        ulist_r[jju] +=
          rootpq *
          (a_r * ulist_r[jjup] +
           a_i * ulist_i[jjup]);
        ulist_i[jju] +=
          rootpq *
          (a_r * ulist_i[jjup] -
           a_i * ulist_r[jjup]);

        rootpq = rootpqarray[(ma + 1)*jdim + (j - mb)];
        ulist_r[jju+1] =
          -rootpq *
          (b_r * ulist_r[jjup] +
           b_i * ulist_i[jjup]);
        ulist_i[jju+1] =
          -rootpq *
          (b_r * ulist_i[jjup] -
           b_i * ulist_r[jjup]);
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
        if (mapar == 1) {
          ulist_r[jjup] = ulist_r[jju];
          ulist_i[jjup] = -ulist_i[jju];
        } else {
          ulist_r[jjup] = -ulist_r[jju];
          ulist_i[jjup] = ulist_i[jju];
        }
        mapar = -mapar;
        jju++;
        jjup--;
      }
      mbpar = -mbpar;
    }
  };
};
template void cpuComputeUarray(double*, double*, double*, double, double, double, double, double,
        int*, int);
template void cpuComputeUarray(float*, float*, float*, float, float, float, float, float,
        int*, int);

template <typename T> void cpuAddUarraytot(T *ulisttot_r, T *ulisttot_i, T *ulist_r, T *ulist_i, T *rij, 
        T *wjelem, T *radelem, T rmin0, T rcutfac, int *idxu_block, int *ilist, 
        int *type, int *pairnum, int *pairnumsum, int *map, int *tj, int twojmax, int idxu_max, 
        int nelements, int inum, int switch_flag, int chemflag)
{
    T rsq, r, x, y, z, rcut, sfac;

    for (int ii=0; ii<inum; ii++) {  // for each atom i in the simulation box     
        int i = ilist[ii];       // atom i
        int itype = type[i];
        int m = pairnum[ii];        // number of neighbors around i             
        int start = pairnumsum[ii];   
        for (int l=0; l<m ; l++) {   // loop over each atom around atom i            
            int k = start + l;    
            int jtype = tj[k];
            int jelem = (chemflag) ? map[jtype] : 0;                
            x = rij[k*3+0];
            y = rij[k*3+1];
            z = rij[k*3+2];
            rsq = x * x + y * y + z * z;
            r = sqrt(rsq);
            rcut = (radelem[itype]+radelem[jtype])*rcutfac;
            sfac = cpuComputeSfac(r, rcut, rmin0, switch_flag);
            sfac *= wjelem[jtype];
                                    
            for (int j = 0; j <= twojmax; j++) {
              int jju = idxu_block[j];
              for (int mb = 0; mb <= j; mb++)
                for (int ma = 0; ma <= j; ma++) {
                  int kl = ii*nelements*idxu_max + jelem*idxu_max+jju;
                  int kr = k*idxu_max + jju;
                  ulisttot_r[kl] += sfac * ulist_r[kr];
                  ulisttot_i[kl] += sfac * ulist_i[kr];
                  jju++;
                }
            }            
        }
    }      
};
template void cpuAddUarraytot(double*, double*, double*, double*, double*, double*, double*, 
        double, double, int*, int*, int*, int*, int*, int*, int*, int, int, int, int, int, int);
template void cpuAddUarraytot(float*, float*, float*, float*, float*, float*, float*, float,
        float, int*, int*, int*, int*, int*, int*, int*, int, int, int, int, int, int);
        
template <typename T> void cpuComputeUij(T *ulist_r, T *ulist_i, T *rootpqarray, T *rij, 
        T *radelem, T rmin0, T rfac0, T rcutfac, int *idxu_block, 
        int *type, int *ai, int *aj, int twojmax, int idxu_max, int ijnum)
{
  T rsq, r, x, y, z, z0, theta0, rcutij;
 
  for(int j = 0; j < ijnum; j++) {
    x = rij[j*3+0];
    y = rij[j*3+1];
    z = rij[j*3+2];
    rsq = x * x + y * y + z * z;
    r = sqrt(rsq);

    int ii = ai[j];
    int jj = aj[j];
    rcutij = (radelem[type[ii]]+radelem[type[jj]])*rcutfac;
    theta0 = (r - rmin0) * rfac0 * M_PI / (rcutij - rmin0);
    //    theta0 = (r - rmin0) * rscale0;
    z0 = r / tan(theta0);    
            
    //compute_uarray(x, y, z, z0, r, j);
    cpuComputeUarray(&ulist_r[idxu_max*j], &ulist_i[idxu_max*j], rootpqarray, 
            x, y, z, z0, r, idxu_block, twojmax);
    
  }
};
template void cpuComputeUij(double*, double*, double*, double*, double*, double, 
        double, double, int*, int*, int*, int*, int, int, int);
template void cpuComputeUij(float*, float*, float*, float*, float*, float,  float,
        float, int*, int*, int*, int*, int, int, int);

template <typename T> void cpuComputeZi(T *zlist_r, T *zlist_i, T *ulisttot_r, T *ulisttot_i, T *cglist,
        int *idxz, int *idxu_block, int *idxcg_block, int twojmax, int idxu_max, int idxz_max, int nelements, 
        int bnorm_flag, int inum)
{

  int jdim = twojmax + 1;  
  
  T * zptr_r;
  T * zptr_i;
  
  for (int ii=0; ii<inum; ii++)       
  for(int elem1 = 0; elem1 < nelements; elem1++)
    for(int elem2 = 0; elem2 < nelements; elem2++) {

      int iT = elem2 + nelements*elem1;  
      zptr_r = &zlist_r[iT*idxz_max + idxz_max*nelements*nelements*ii];
      zptr_i = &zlist_i[iT*idxz_max + idxz_max*nelements*nelements*ii];
      
      for (int jjz = 0; jjz < idxz_max; jjz++) {
        const int j1 = idxz[jjz*10+0];
        const int j2 = idxz[jjz*10+1];
        const int j = idxz[jjz*10+2];
        const int ma1min = idxz[jjz*10+3];
        const int ma2max = idxz[jjz*10+4];
        const int na = idxz[jjz*10+5];
        const int mb1min = idxz[jjz*10+6];
        const int mb2max = idxz[jjz*10+7];
        const int nb = idxz[jjz*10+8];

        //const T *cgblock = cglist + idxcg_block[j1][j2][j];
        const T *cgblock = cglist + idxcg_block[j + j2*jdim + j1*jdim*jdim];

        zptr_r[jjz] = 0.0;
        zptr_i[jjz] = 0.0;      

        int jju1 = idxu_block[j1] + (j1 + 1) * mb1min;
        int jju2 = idxu_block[j2] + (j2 + 1) * mb2max;
        int icgb = mb1min * (j2 + 1) + mb2max;
        for (int ib = 0; ib < nb; ib++) {

          T suma1_r = 0.0;
          T suma1_i = 0.0;

          const T *u1_r = &ulisttot_r[elem1*idxu_max+jju1 + idxu_max*nelements*ii];
          const T *u1_i = &ulisttot_i[elem1*idxu_max+jju1 + idxu_max*nelements*ii];
          const T *u2_r = &ulisttot_r[elem2*idxu_max+jju2 + idxu_max*nelements*ii];
          const T *u2_i = &ulisttot_i[elem2*idxu_max+jju2 + idxu_max*nelements*ii];

          int ma1 = ma1min;
          int ma2 = ma2max;
          int icga = ma1min * (j2 + 1) + ma2max;

          for (int ia = 0; ia < na; ia++) {
            suma1_r += cgblock[icga] * (u1_r[ma1] * u2_r[ma2] - u1_i[ma1] * u2_i[ma2]);
            suma1_i += cgblock[icga] * (u1_r[ma1] * u2_i[ma2] + u1_i[ma1] * u2_r[ma2]);
            ma1++;
            ma2--;
            icga += j2;
          } // end loop over ia

          zptr_r[jjz] += cgblock[icgb] * suma1_r;
          zptr_i[jjz] += cgblock[icgb] * suma1_i;

          jju1 += j1 + 1;
          jju2 -= j2 + 1;
          icgb += j2;
        } // end loop over ib
        if (bnorm_flag) {
          zptr_r[jjz] /= (j+1);
          zptr_i[jjz] /= (j+1);
        }
      } // end loop over jjz
    }
};
template void cpuComputeZi(double*, double*, double*, double*, double*, 
        int*, int*, int*, int, int, int, int, int, int);
template void cpuComputeZi(float*, float*, float*, float*, float*, 
        int*, int*, int*, int, int, int, int, int, int);

template <typename T> void cpuComputeYi(T *ylist_r, T *ylist_i, T *ulisttot_r, T *ulisttot_i, T *cglist, 
        T* betalist, int *idxz, int *idxb_block, int *idxu_block, int *idxcg_block, int twojmax, 
        int idxb_max, int idxu_max, int idxz_max, int nelements, int bnorm_flag, int ncoeff, int inum)
{
  int jdim = twojmax + 1;  
  int jju;
  T betaj;
  int itriple;

//   for (int ii=0; ii<inum; ii++)       
//   for(int ielem1 = 0; ielem1 < nelements; ielem1++)
//     for(int j = 0; j <= twojmax; j++) {
//       jju = idxu_block[j];
//       for(int mb = 0; 2*mb <= j; mb++)
//         for(int ma = 0; ma <= j; ma++) {
//           ylist_r[ielem1*idxu_max+jju + idxu_max*nelements*ii] = 0.0;
//           ylist_i[ielem1*idxu_max+jju + idxu_max*nelements*ii] = 0.0;
//           jju++;
//         } // end loop over ma, mb
//     } // end loop over j

  for (int ii=0; ii<inum; ii++)       
  for(int ielem1 = 0; ielem1 < nelements; ielem1++)
    for(int j = 0; j <= idxu_max; j++) {
          ylist_r[ielem1*idxu_max+j + idxu_max*nelements*ii] = 0.0;
          ylist_i[ielem1*idxu_max+j + idxu_max*nelements*ii] = 0.0;
          jju++;
    } // end loop over j

  for (int ii=0; ii<inum; ii++) {   
  T *beta = &betalist[ncoeff*ii];    
  for(int elem1 = 0; elem1 < nelements; elem1++)
    for (int elem2 = 0; elem2 < nelements; elem2++) {
        for (int jjz = 0; jjz < idxz_max; jjz++) {
        const int j1 = idxz[jjz*10+0];
        const int j2 = idxz[jjz*10+1];
        const int j = idxz[jjz*10+2];
        const int ma1min = idxz[jjz*10+3];
        const int ma2max = idxz[jjz*10+4];
        const int na = idxz[jjz*10+5];
        const int mb1min = idxz[jjz*10+6];
        const int mb2max = idxz[jjz*10+7];
        const int nb = idxz[jjz*10+8];
          
          //const T *cgblock = cglist + idxcg_block[j1][j2][j];
          const T *cgblock = cglist + idxcg_block[j + j2*jdim + j1*jdim*jdim];
          
          T ztmp_r = 0.0;
          T ztmp_i = 0.0;

          int jju1 = idxu_block[j1] + (j1 + 1) * mb1min;
          int jju2 = idxu_block[j2] + (j2 + 1) * mb2max;
          int icgb = mb1min * (j2 + 1) + mb2max;
          for (int ib = 0; ib < nb; ib++) {

            T suma1_r = 0.0;
            T suma1_i = 0.0;

            const T *u1_r = &ulisttot_r[elem1*idxu_max+jju1 + idxu_max*nelements*ii];
            const T *u1_i = &ulisttot_i[elem1*idxu_max+jju1 + idxu_max*nelements*ii];
            const T *u2_r = &ulisttot_r[elem2*idxu_max+jju2 + idxu_max*nelements*ii];
            const T *u2_i = &ulisttot_i[elem2*idxu_max+jju2 + idxu_max*nelements*ii];

            int ma1 = ma1min;
            int ma2 = ma2max;
            int icga = ma1min * (j2 + 1) + ma2max;

            for (int ia = 0; ia < na; ia++) {
              suma1_r += cgblock[icga] * (u1_r[ma1] * u2_r[ma2] - u1_i[ma1] * u2_i[ma2]);
              suma1_i += cgblock[icga] * (u1_r[ma1] * u2_i[ma2] + u1_i[ma1] * u2_r[ma2]);
              ma1++;
              ma2--;
              icga += j2;
            } // end loop over ia


            ztmp_r += cgblock[icgb] * suma1_r;
            ztmp_i += cgblock[icgb] * suma1_i;

            jju1 += j1 + 1;
            jju2 -= j2 + 1;
            icgb += j2;
          } // end loop over ib

          // apply to z(j1,j2,j,ma,mb) to unique element of y(j)
          // find right y_list[jju] and beta[jjb] entries
          // multiply and divide by j+1 factors
          // account for multiplicity of 1, 2, or 3

        if (bnorm_flag) {
          ztmp_i /= j+1;
          ztmp_r /= j+1;
        }

        //jju = idxz[jjz].jju;
        jju = idxz[jjz*10+9];
        for(int elem3 = 0; elem3 < nelements; elem3++) {
        // pick out right beta value
          if (j >= j1) {
            //const int jjb = idxb_block[j1][j2][j];
            const int jjb = idxb_block[j + j2*jdim + j1*jdim*jdim];
            itriple = ((elem1 * nelements + elem2) * nelements + elem3) * idxb_max + jjb;
            if (j1 == j) {
              if (j2 == j) betaj = 3*beta[itriple];
              else betaj = 2*beta[itriple];
            } else betaj = beta[itriple];
          } else if (j >= j2) {
            //const int jjb = idxb_block[j][j2][j1];
              const int jjb = idxb_block[j1 + j2*jdim + j*jdim*jdim];
            itriple = ((elem3 * nelements + elem2) * nelements + elem1) * idxb_max + jjb;
            if (j2 == j) betaj = 2*beta[itriple];
            else betaj = beta[itriple];
          } else {
            //const int jjb = idxb_block[j2][j][j1];
            const int jjb = idxb_block[j1 + j*jdim + j2*jdim*jdim];
            itriple = ((elem2 * nelements + elem3) * nelements + elem1) * idxb_max + jjb;
            betaj = beta[itriple];
          }

          if (!bnorm_flag && j1 > j)
            betaj *= (j1 + 1) / (j + 1.0);

          ylist_r[elem3 * idxu_max + jju + idxu_max*nelements*ii] += betaj * ztmp_r;
          ylist_i[elem3 * idxu_max + jju + idxu_max*nelements*ii] += betaj * ztmp_i;
        }
      } // end loop over jjz
    }
  }
}
template void cpuComputeYi(double*, double*, double*, double*, double*, double*, 
        int*, int*, int*, int*, int, int, int, int, int, int, int, int);
template void cpuComputeYi(float*, float*, float*, float*, float*, float*,
        int*, int*, int*, int*, int, int, int, int, int, int, int, int);

template <typename T> void cpuComputeBi(T *blist, T *zlist_r, T *zlist_i, T *ulisttot_r, T *ulisttot_i, 
        T *bzero,int *ilist, int *type, int *map, int *idxb, int *idxu_block, int *idxz_block, 
        int twojmax, int idxb_max, int idxu_max, int idxz_max, int nelements, int bzero_flag, 
        int wselfall_flag, int chemflag, int inum)
{
  // for j1 = 0,...,twojmax
  //   for j2 = 0,twojmax
  //     for j = |j1-j2|,Min(twojmax,j1+j2),2
  //        b(j1,j2,j) = 0
  //        for mb = 0,...,jmid
  //          for ma = 0,...,j
  //            b(j1,j2,j) +=
  //              2*Conj(u(j,ma,mb))*z(j1,j2,j,ma,mb)

  int jdim = twojmax + 1;  
  
  for (int ii=0; ii<inum; ii++)       
  for(int elem1 = 0; elem1 < nelements; elem1++)
    for(int elem2 = 0; elem2 < nelements; elem2++) {

      int iT = elem2 + nelements*elem1;  
      T *zptr_r = &zlist_r[iT*idxz_max + idxz_max*nelements*nelements*ii];
      T *zptr_i = &zlist_i[iT*idxz_max + idxz_max*nelements*nelements*ii];
      
      for (int elem3 = 0; elem3 < nelements; elem3++) {
        int itriple = elem3 + nelements*elem2 + nelements*nelements*elem1;    
        for (int jjb = 0; jjb < idxb_max; jjb++) {
          //const int j1 = idxb[jjb].j1;
          //const int j2 = idxb[jjb].j2;
          //const int j = idxb[jjb].j;
          const int j1 = idxb[jjb*3 + 0];
          const int j2 = idxb[jjb*3 + 1];
          const int j = idxb[jjb*3 + 2];  
          
          //int jjz = idxz_block[j1][j2][j];
          int jjz = idxz_block[j + j2*jdim + j1*jdim*jdim];
          int jju = idxu_block[j];
          T sumzu = 0.0;
          for (int mb = 0; 2 * mb < j; mb++)
            for (int ma = 0; ma <= j; ma++) {
              sumzu += ulisttot_r[elem3*idxu_max+jju + idxu_max*nelements*ii] * zptr_r[jjz] +
                       ulisttot_i[elem3*idxu_max+jju + idxu_max*nelements*ii] * zptr_i[jjz];
              jjz++;
              jju++;
            } // end loop over ma, mb

          // For j even, handle middle column

          if (j % 2 == 0) {
            int mb = j / 2;
            for (int ma = 0; ma < mb; ma++) {
              sumzu += ulisttot_r[elem3*idxu_max+jju + idxu_max*nelements*ii] * zptr_r[jjz] +
                       ulisttot_i[elem3*idxu_max+jju + idxu_max*nelements*ii] * zptr_i[jjz];
              jjz++;
              jju++;
            }

            sumzu += 0.5 * (ulisttot_r[elem3*idxu_max+jju + idxu_max*nelements*ii] * zptr_r[jjz] +
                            ulisttot_i[elem3*idxu_max+jju + idxu_max*nelements*ii] * zptr_i[jjz]);
          } // end if jeven

          blist[itriple*idxb_max+jjb + idxb_max*nelements*nelements*nelements*ii] = 2.0 * sumzu;

        }        
      }      
    }

  // apply bzero shift

  if (bzero_flag) {
    if (!wselfall_flag) {
      for (int ii=0; ii<inum; ii++) {        
          int ielem = (chemflag) ? map[type[ilist[ii]]]: 0;                
          int itriple = (ielem*nelements+ielem)*nelements+ielem;
          for (int jjb = 0; jjb < idxb_max; jjb++) {
            //const int j = idxb[jjb].j;
            const int j = idxb[jjb*3 + 2];  
            blist[itriple*idxb_max+jjb + idxb_max*nelements*nelements*nelements*ii] -= bzero[j];
          } // end loop over JJ
      }
    } else {
      for (int ii=0; ii<inum; ii++) {                  
          for(int elem1 = 0; elem1 < nelements; elem1++)
            for(int elem2 = 0; elem2 < nelements; elem2++) {
              for(int elem3 = 0; elem3 < nelements; elem3++) {
                int itriple = elem3 + nelements*elem2 + nelements*nelements*elem1;      
                for (int jjb = 0; jjb < idxb_max; jjb++) {
                  //const int j = idxb[jjb].j;
                  const int j = idxb[jjb*3 + 2];  
                  blist[itriple*idxb_max+jjb+ idxb_max*nelements*nelements*nelements*ii] -= bzero[j];
                } // end loop over JJ
              } // end loop over elem3
            } // end loop over elem1,elem2
      }
    }
  }
};
template void cpuComputeBi(double*, double*, double*, double*, double*, double*, 
        int*, int*, int*, int*, int*, int*, int, int, int, int, int, int, int, int, int);
template void cpuComputeBi(float*, float*, float*, float*, float*, float*,
        int*, int*, int*, int*, int*, int*, int, int, int, int, int, int, int, int, int);

template <typename T> void cpuComputeDuarray(T *dulist_r, T *dulist_i, T *ulist_r, T *ulist_i, T *rootpqarray,
                     T x, T y, T z, T z0, T r, T dz0dr, T wj, 
                     T rcut, T rmin0, int *idxu_block, int twojmax, int switch_flag)
{
  T r0inv;
  T a_r, a_i, b_r, b_i;
  T da_r[3], da_i[3], db_r[3], db_i[3];
  T dz0[3], dr0inv[3], dr0invdr;
  T rootpq;

  T rinv = 1.0 / r;
  T ux = x * rinv;
  T uy = y * rinv;
  T uz = z * rinv;

  int jdim = twojmax + 1;
  
  r0inv = 1.0 / sqrt(r * r + z0 * z0);
  a_r = z0 * r0inv;
  a_i = -z * r0inv;
  b_r = y * r0inv;
  b_i = -x * r0inv;

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

  //T* ulist_r = ulist_r_ij[jj];
  //T* ulist_i = ulist_i_ij[jj];

  dulist_r[0*3+0] = 0.0;
  dulist_r[0*3+1] = 0.0;
  dulist_r[0*3+2] = 0.0;
  dulist_i[0*3+0] = 0.0;
  dulist_i[0*3+1] = 0.0;
  dulist_i[0*3+2] = 0.0;

  for (int j = 1; j <= twojmax; j++) {
    int jju = idxu_block[j];
    int jjup = idxu_block[j-1];
    for (int mb = 0; 2*mb <= j; mb++) {
      dulist_r[jju*3+0] = 0.0;
      dulist_r[jju*3+1] = 0.0;
      dulist_r[jju*3+2] = 0.0;
      dulist_i[jju*3+0] = 0.0;
      dulist_i[jju*3+1] = 0.0;
      dulist_i[jju*3+2] = 0.0;

      for (int ma = 0; ma < j; ma++) {
        rootpq = rootpqarray[(j - ma)*jdim + (j - mb)];
        for (int k = 0; k < 3; k++) {
          dulist_r[jju*3+k] +=
            rootpq * (da_r[k] * ulist_r[jjup] +
                      da_i[k] * ulist_i[jjup] +
                      a_r * dulist_r[jjup*3+k] +
                      a_i * dulist_i[jjup*3+k]);
          dulist_i[jju*3+k] +=
            rootpq * (da_r[k] * ulist_i[jjup] -
                      da_i[k] * ulist_r[jjup] +
                      a_r * dulist_i[jjup*3+k] -
                      a_i * dulist_r[jjup*3+k]);
        }

        rootpq = rootpqarray[(ma + 1)*jdim + (j - mb)];
        for (int k = 0; k < 3; k++) {
          dulist_r[(jju+1)*3+k] =
            -rootpq * (db_r[k] * ulist_r[jjup] +
                       db_i[k] * ulist_i[jjup] +
                       b_r * dulist_r[jjup*3+k] +
                       b_i * dulist_i[jjup*3+k]);
          dulist_i[(jju+1)*3+k] =
            -rootpq * (db_r[k] * ulist_i[jjup] -
                       db_i[k] * ulist_r[jjup] +
                       b_r * dulist_i[jjup*3+k] -
                       b_i * dulist_r[jjup*3+k]);
        }
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
        if (mapar == 1) {
          for (int k = 0; k < 3; k++) {
            dulist_r[jjup*3+k] = dulist_r[jju*3+k];
            dulist_i[jjup*3+k] = -dulist_i[jju*3+k];
          }
        } else {
          for (int k = 0; k < 3; k++) {
            dulist_r[jjup*3+k] = -dulist_r[jju*3+k];
            dulist_i[jjup*3+k] = dulist_i[jju*3+k];
          }
        }
        mapar = -mapar;
        jju++;
        jjup--;
      }
      mbpar = -mbpar;
    }
  }

  //T sfac = compute_sfac(r, rcut);
  //T dsfac = compute_dsfac(r, rcut);
  T sfac = cpuComputeSfac(r, rcut, rmin0, switch_flag);  
  T dsfac = cpuComputeDsfac(r, rcut, rmin0, switch_flag);  
    
  sfac *= wj;
  dsfac *= wj;
  for (int j = 0; j <= twojmax; j++) {
    int jju = idxu_block[j];
    for (int mb = 0; 2*mb <= j; mb++)
      for (int ma = 0; ma <= j; ma++) {
        dulist_r[jju*3+0] = dsfac * ulist_r[jju] * ux +
                                  sfac * dulist_r[jju*3+0];
        dulist_i[jju*3+0] = dsfac * ulist_i[jju] * ux +
                                  sfac * dulist_i[jju*3+0];
        dulist_r[jju*3+1] = dsfac * ulist_r[jju] * uy +
                                  sfac * dulist_r[jju*3+1];
        dulist_i[jju*3+1] = dsfac * ulist_i[jju] * uy +
                                  sfac * dulist_i[jju*3+1];
        dulist_r[jju*3+2] = dsfac * ulist_r[jju] * uz +
                                  sfac * dulist_r[jju*3+2];
        dulist_i[jju*3+2] = dsfac * ulist_i[jju] * uz +
                                  sfac * dulist_i[jju*3+2];
        jju++;
      }
  }
}
template void cpuComputeDuarray(double*, double*, double*, double*, double*, double, double, double, 
        double, double, double, double, double, double, int*, int, int);
template void cpuComputeDuarray(float*, float*, float*, float*, float*, float, float, float,
        float, float, float, float, float, float, int*, int, int);

template <typename T> int cpuComputeDuidrj(T *dulist_r, T *dulist_i, T *ulist_r, T *ulist_i, T *rootpqarray, 
        T* rij, T wj, T rcut, T rmin0, T rfac0, int *idxu_block, 
        int jelem, int twojmax, int switch_flag)
{
  T rsq, r, x, y, z, z0, theta0, cs, sn;
  T dz0dr;

  x = rij[0];
  y = rij[1];
  z = rij[2];
  rsq = x * x + y * y + z * z;
  r = sqrt(rsq);
  T rscale0 = rfac0 * M_PI / (rcut - rmin0);
  theta0 = (r - rmin0) * rscale0;
  cs = cos(theta0);
  sn = sin(theta0);
  z0 = r * cs / sn;
  dz0dr = z0 / r - (r*rscale0) * (rsq + z0 * z0) / rsq;

  //elem_duarray = jelem;
  cpuComputeDuarray(dulist_r, dulist_i, ulist_r, ulist_i, rootpqarray, x, y, z, z0, r, dz0dr, wj, rcut, rmin0,
          idxu_block, twojmax, switch_flag);
  
  return jelem;
}
template int cpuComputeDuidrj(double*, double*, double*, double*, double*, double*, double, double, 
        double, double, int*, int, int, int);
template int cpuComputeDuidrj(float*, float*, float*, float*, float*, float*, float, float, 
        float, float, int*, int, int, int);

template <typename T> void cpuComputeDuijdrj(T *dulist_r, T *dulist_i, T *ulist_r, T *ulist_i, T *rootpqarray, 
        T* rij, T *wjelem, T *radelem, T rmin0, T rfac0, T rcutfac, int *idxu_block, 
        int *type, int *ai, int *aj, int twojmax, int idxu_max, int ijnum, int switch_flag)
{
  T rsq, r, x, y, z, z0, theta0, cs, sn, rcutij, rscale0;
  T dz0dr;

  for(int j = 0; j < ijnum; j++) {
    x = rij[j*3+0];
    y = rij[j*3+1];
    z = rij[j*3+2];
    rsq = x * x + y * y + z * z;
    r = sqrt(rsq);

    int ii = ai[j];
    int jj = aj[j];
    rcutij = (radelem[type[ii]]+radelem[type[jj]])*rcutfac;
    rscale0 = rfac0 * M_PI / (rcutij - rmin0);
    theta0 = (r - rmin0) * rscale0;
    z0 = r / tan(theta0);    
    cs = cos(theta0);
    sn = sin(theta0);
    z0 = r * cs / sn;
    dz0dr = z0 / r - (r*rscale0) * (rsq + z0 * z0) / rsq;
    
    cpuComputeDuarray(&dulist_r[3*idxu_max*j], &dulist_i[3*idxu_max*j], &ulist_r[idxu_max*j], 
            &ulist_i[idxu_max*j], rootpqarray, x, y, z, z0, r, dz0dr, wjelem[type[jj]], rcutij, rmin0,
            idxu_block, twojmax, switch_flag);    
  }  
}
template void cpuComputeDuijdrj(double*, double*, double*, double*, double*, double*, 
        double*, double*, double, double, double, int*, int*, int*, int*, int, int, int, int);
template void cpuComputeDuijdrj(float*, float*, float*, float*, float*, float*, float*, 
        float*, float, float, float, int*, int*, int*, int*, int, int, int, int);

template <typename T> void cpuComputeDeidrj(T *dedr, T *ylist_r, T *ylist_i, T *dulist_r, T *dulist_i,         
        int *idxu_block, int *type, int *map, int *ai, int *aj, int nelements, int twojmax, int idxu_max, 
        int chemflag, int ijnum) 
{

  for(int ij = 0; ij < ijnum; ij++)  
    for(int k = 0; k < 3; k++)
        dedr[k + 3*ij] = 0.0;

  for(int ij = 0; ij < ijnum; ij++) {
  int jelem = (chemflag) ? map[type[aj[ij]]] : 0;
  int i = ai[ij]; // atom i  
  for(int j = 0; j <= twojmax; j++) {
    int jju = idxu_block[j];

    for(int mb = 0; 2*mb < j; mb++)
      for(int ma = 0; ma <= j; ma++) {
          
        T* dudr_r = &dulist_r[jju*3 + 3*idxu_max*ij];
        T* dudr_i = &dulist_i[jju*3 + 3*idxu_max*ij];
        T jjjmambyarray_r = ylist_r[jelem*idxu_max+jju+idxu_max*nelements*i];
        T jjjmambyarray_i = ylist_i[jelem*idxu_max+jju+idxu_max*nelements*i];

        for(int k = 0; k < 3; k++)
          dedr[k+ 3*ij] +=
            dudr_r[k] * jjjmambyarray_r +
            dudr_i[k] * jjjmambyarray_i;
        jju++;
      } //end loop over ma mb

    // For j even, handle middle column

    if (j%2 == 0) {

      int mb = j/2;
      for(int ma = 0; ma < mb; ma++) {
        T* dudr_r = &dulist_r[jju*3 + 3*idxu_max*ij];
        T* dudr_i = &dulist_i[jju*3 + 3*idxu_max*ij];
        T jjjmambyarray_r = ylist_r[jelem*idxu_max+jju+idxu_max*nelements*i];
        T jjjmambyarray_i = ylist_i[jelem*idxu_max+jju+idxu_max*nelements*i];

        for(int k = 0; k < 3; k++)
          dedr[k + 3*ij] +=
            dudr_r[k] * jjjmambyarray_r +
            dudr_i[k] * jjjmambyarray_i;
        jju++;
      }

      T* dudr_r = &dulist_r[jju*3 + 3*idxu_max*ij];
      T* dudr_i = &dulist_i[jju*3 + 3*idxu_max*ij];
      T jjjmambyarray_r = ylist_r[jelem*idxu_max+jju+idxu_max*nelements*i];
      T jjjmambyarray_i = ylist_i[jelem*idxu_max+jju+idxu_max*nelements*i];

      for(int k = 0; k < 3; k++)
        dedr[k + 3*ij] +=
          (dudr_r[k] * jjjmambyarray_r +
           dudr_i[k] * jjjmambyarray_i)*0.5;
      // jju++;

    } // end if jeven

  } // end loop over j

  for(int k = 0; k < 3; k++)
    dedr[k + 3*ij] *= 2.0;
  }
  
}
template void cpuComputeDeidrj(double*, double*, double*, double*, double*, 
        int*, int*, int*, int*, int*, int, int, int, int, int);
template void cpuComputeDeidrj(float*, float*, float*, float*, float*, 
        int*, int*, int*, int*, int*, int, int, int, int, int);

template <typename T> void cpuComputeDbidrj(T *dblist, T *zlist_r, T *zlist_i, T *dulist_r, T *dulist_i, 
        int *idxb, int *idxu_block, int *idxz_block, int *type, int *map, int *ai, int *aj, int twojmax, 
        int idxb_max, int idxu_max, int idxz_max, int nelements, int bnorm_flag, int chemflag, int ijnum)
{  
  int jdim = twojmax + 1;  
  T* dbdr;
  T* dudr_r, *dudr_i;
  T sumzdu_r[3];
  T* zptr_r;
  T* zptr_i;
  int jjz, jju;

  int iT;
  int itriple;

  // set all the derivatives to zero once

  for(int ij = 0; ij < ijnum; ij++)
  for(int jjb = 0; jjb < idxb_max; jjb++) {

    for(int elem1 = 0; elem1 < nelements; elem1++)
      for(int elem2 = 0; elem2 < nelements; elem2++)
        for(int elem3 = 0; elem3 < nelements; elem3++) {

          itriple = (elem1 * nelements + elem2) * nelements + elem3;

          dbdr = &dblist[(itriple*idxb_max+jjb)*3 + 3*idxb_max*nelements*nelements*nelements*ij];
          dbdr[0] = 0.0;
          dbdr[1] = 0.0;
          dbdr[2] = 0.0;
        }

  }

  for(int ij = 0; ij < ijnum; ij++) {
  int elem3 = (chemflag) ? map[type[aj[ij]]] : 0;
  int i = ai[ij]; // atom i
  for(int jjb = 0; jjb < idxb_max; jjb++) {
    const int j1 = idxb[jjb*3 + 0];
    const int j2 = idxb[jjb*3 + 1];
    const int j = idxb[jjb*3 + 2];  

    // Sum terms Conj(dudr(j,ma,mb))*z(j1,j2,j,ma,mb)

    for(int elem1 = 0; elem1 < nelements; elem1++)
      for(int elem2 = 0; elem2 < nelements; elem2++) {

        //jjz = idxz_block[j1][j2][j];
        jjz = idxz_block[j + j2*jdim + j1*jdim*jdim];
        jju = idxu_block[j];
        iT = elem1*nelements+elem2;
        itriple = (elem1*nelements+elem2)*nelements+elem3;
        dbdr = &dblist[(itriple*idxb_max+jjb)*3 + 3*idxb_max*nelements*nelements*nelements*ij];
        zptr_r = &zlist_r[iT*idxz_max + idxz_max*nelements*nelements*i];
        zptr_i = &zlist_i[iT*idxz_max + idxz_max*nelements*nelements*i];

        for (int k = 0; k < 3; k++)
          sumzdu_r[k] = 0.0;

        for (int mb = 0; 2 * mb < j; mb++)
          for (int ma = 0; ma <= j; ma++) {
            dudr_r = &dulist_r[jju*3 + 3*idxu_max*ij];
            dudr_i = &dulist_i[jju*3 + 3*idxu_max*ij];
            for (int k = 0; k < 3; k++)
              sumzdu_r[k] +=
                  dudr_r[k] * zptr_r[jjz] +
                  dudr_i[k] * zptr_i[jjz];
            jjz++;
            jju++;
          } //end loop over ma mb

        // For j even, handle middle column

        if (j % 2 == 0) {
          int mb = j / 2;
          for (int ma = 0; ma < mb; ma++) {
            dudr_r = &dulist_r[jju*3 + 3*idxu_max*ij];
            dudr_i = &dulist_i[jju*3 + 3*idxu_max*ij];
            for (int k = 0; k < 3; k++)
              sumzdu_r[k] +=
                  dudr_r[k] * zptr_r[jjz] +
                  dudr_i[k] * zptr_i[jjz];
            jjz++;
            jju++;
          }
          dudr_r = &dulist_r[jju*3 + 3*idxu_max*ij];
          dudr_i = &dulist_i[jju*3 + 3*idxu_max*ij];
          for (int k = 0; k < 3; k++)
            sumzdu_r[k] +=
                (dudr_r[k] * zptr_r[jjz] +
                 dudr_i[k] * zptr_i[jjz]) * 0.5;
          // jjz++;
          // jju++;
        } // end if jeven

        for (int k = 0; k < 3; k++)
          dbdr[k] += 2.0 * sumzdu_r[k];
        // Sum over Conj(dudr(j1,ma1,mb1))*z(j,j2,j1,ma1,mb1)

        T j1fac = (j + 1) / (j1 + 1.0);

        iT = elem1*nelements+elem2;
        itriple = (elem3*nelements+elem2)*nelements+elem1;
        dbdr = &dblist[(itriple*idxb_max+jjb)*3 + 3*idxb_max*nelements*nelements*nelements*ij];
        //jjz = idxz_block[j][j2][j1];
        jjz = idxz_block[j1 + j2*jdim + j*jdim*jdim];
        jju = idxu_block[j1];
        zptr_r = &zlist_r[iT*idxz_max  + idxz_max*nelements*nelements*i];
        zptr_i = &zlist_i[iT*idxz_max  + idxz_max*nelements*nelements*i];

        for (int k = 0; k < 3; k++)
          sumzdu_r[k] = 0.0;

        for (int mb = 0; 2 * mb < j1; mb++)
          for (int ma = 0; ma <= j1; ma++) {
            dudr_r = &dulist_r[jju*3 + 3*idxu_max*ij];
            dudr_i = &dulist_i[jju*3 + 3*idxu_max*ij];
            for (int k = 0; k < 3; k++)
              sumzdu_r[k] +=
                  dudr_r[k] * zptr_r[jjz] +
                  dudr_i[k] * zptr_i[jjz];
            jjz++;
            jju++;
          } //end loop over ma mb

        // For j1 even, handle middle column

        if (j1 % 2 == 0) {
          int mb = j1 / 2;
          for (int ma = 0; ma < mb; ma++) {
            dudr_r = &dulist_r[jju*3 + 3*idxu_max*ij];
            dudr_i = &dulist_i[jju*3 + 3*idxu_max*ij];
            for (int k = 0; k < 3; k++)
              sumzdu_r[k] +=
                  dudr_r[k] * zptr_r[jjz] +
                  dudr_i[k] * zptr_i[jjz];
            jjz++;
            jju++;
          }
          dudr_r = &dulist_r[jju*3 + 3*idxu_max*ij];
          dudr_i = &dulist_i[jju*3 + 3*idxu_max*ij];
          for (int k = 0; k < 3; k++)
            sumzdu_r[k] +=
                (dudr_r[k] * zptr_r[jjz] +
                 dudr_i[k] * zptr_i[jjz]) * 0.5;
          // jjz++;
          // jju++;
        } // end if j1even

        for (int k = 0; k < 3; k++)
          if (bnorm_flag)
            dbdr[k] += 2.0 * sumzdu_r[k];
          else
            dbdr[k] += 2.0 * sumzdu_r[k] * j1fac;

        // Sum over Conj(dudr(j2,ma2,mb2))*z(j,j1,j2,ma2,mb2)

        T j2fac = (j + 1) / (j2 + 1.0);

        iT = elem2*nelements+elem1;
        itriple = (elem1*nelements+elem3)*nelements+elem2;
        dbdr = &dblist[(itriple*idxb_max+jjb)*3 + 3*idxb_max*nelements*nelements*nelements*ij];
        //jjz = idxz_block[j][j1][j2];
        int jjz = idxz_block[j2 + j1*jdim + j*jdim*jdim];        
        jju = idxu_block[j2];
        zptr_r = &zlist_r[iT*idxz_max  + idxz_max*nelements*nelements*i];
        zptr_i = &zlist_i[iT*idxz_max  + idxz_max*nelements*nelements*i];

        for (int k = 0; k < 3; k++)
          sumzdu_r[k] = 0.0;

        for (int mb = 0; 2 * mb < j2; mb++)
          for (int ma = 0; ma <= j2; ma++) {
            dudr_r = &dulist_r[jju*3 + 3*idxu_max*ij];
            dudr_i = &dulist_i[jju*3 + 3*idxu_max*ij];
            for (int k = 0; k < 3; k++)
              sumzdu_r[k] +=
                  dudr_r[k] * zptr_r[jjz] +
                  dudr_i[k] * zptr_i[jjz];
            jjz++;
            jju++;
          } //end loop over ma mb

        // For j2 even, handle middle column

        if (j2 % 2 == 0) {
          int mb = j2 / 2;
          for (int ma = 0; ma < mb; ma++) {
            dudr_r = &dulist_r[jju*3 + 3*idxu_max*ij];
            dudr_i = &dulist_i[jju*3 + 3*idxu_max*ij];
            for (int k = 0; k < 3; k++)
              sumzdu_r[k] +=
                  dudr_r[k] * zptr_r[jjz] +
                  dudr_i[k] * zptr_i[jjz];
            jjz++;
            jju++;
          }
          dudr_r = &dulist_r[jju*3 + 3*idxu_max*ij];
          dudr_i = &dulist_i[jju*3 + 3*idxu_max*ij];
          for (int k = 0; k < 3; k++)
            sumzdu_r[k] +=
                (dudr_r[k] * zptr_r[jjz] +
                 dudr_i[k] * zptr_i[jjz]) * 0.5;
          // jjz++;
          // jju++;
        } // end if j2even

        for (int k = 0; k < 3; k++)
          if (bnorm_flag)
            dbdr[k] += 2.0 * sumzdu_r[k];
          else
            dbdr[k] += 2.0 * sumzdu_r[k] * j2fac;

      }
  } //end loop over j1 j2 j
  }
}
template void cpuComputeDbidrj(double*, double*, double*, double*, double*, 
        int*, int*, int*, int*, int*, int*, int*, int, int, int, int, int, int, int, int);
template void cpuComputeDbidrj(float*, float*, float*, float*, float*, 
        int*, int*, int*, int*, int*, int*, int*, int, int, int, int, int, int, int, int);

template <typename T> void cpuComputeSna(T *sna, T *blist, int *ilist, int *mask, 
        int ncoeff, int nrows, int inum, int quadraticflag)
{
   for (int ii = 0; ii<inum; ii++) {
      int i = ilist[ii];
      for (int icoeff = 0; icoeff < ncoeff; icoeff++) {          
          T bi = blist[icoeff + ncoeff*ii];
          sna[icoeff + nrows*i] = (mask[i]) ? bi : 0.0;
      }
   }
   
   if (quadraticflag) {
        for (int ii = 0; ii<inum; ii++) {
            int i = ilist[ii];       
            //int ncount = ncoeff;
            for (int icoeff = 0; icoeff < ncoeff; icoeff++) {             
                T bi = blist[icoeff + ncoeff*ii];
                
                // diagonal element of quadratic matrix
                int ncount = ncoeff + icoeff*ncoeff - icoeff*(icoeff-1)/2;
                sna[ncount + nrows*i] = (mask[i]) ? 0.5*bi*bi : 0.0;
                ncount++;
                
                // upper-triangular elements of quadratic matrix
                for (int jcoeff = icoeff+1; jcoeff < ncoeff; jcoeff++) {
                    sna[ncount + nrows*i] = (mask[i]) ? bi*blist[jcoeff + ncoeff*ii] : 0.0;                        
                    ncount++;
                }
            }       
        }
    }   
}
template void cpuComputeSna(double*, double*, int*, int*, int, int, int, int);
template void cpuComputeSna(float*, float*, int*, int*, int, int, int, int);

template <typename T> void cpuComputeSna(T *sna, T *blist, int *ilist, int *mask, int *type,
        int ncoeff, int ntype, int nperdim, int inum, int quadraticflag)
{
   for (int ii = 0; ii<inum; ii++) {
      int i = ilist[ii];
      int itype = type[i]-1;
      int imsk = mask[i];
      for (int icoeff = 0; icoeff < ncoeff; icoeff++)                     
          sna[icoeff + nperdim*itype + nperdim*ntype*i] = (imsk) ? blist[icoeff + ncoeff*ii] : 0.0;      
   }
   
   if (quadraticflag) {
        for (int ii = 0; ii<inum; ii++) {
            int i = ilist[ii];       
            int itype = type[i]-1;
            int imsk = mask[i];
            for (int icoeff = 0; icoeff < ncoeff; icoeff++) {             
                T bi = blist[icoeff + ncoeff*ii];
                
                // diagonal element of quadratic matrix
                int ncount = ncoeff + icoeff*ncoeff - icoeff*(icoeff-1)/2;
                sna[ncount + nperdim*itype + nperdim*ntype*i] = (imsk) ? 0.5*bi*bi : 0.0;
                ncount++;
                
                // upper-triangular elements of quadratic matrix
                for (int jcoeff = icoeff+1; jcoeff < ncoeff; jcoeff++) {
                    sna[ncount + nperdim*itype + nperdim*ntype*i] = (imsk) ? bi*blist[jcoeff + ncoeff*ii] : 0.0;                        
                    ncount++;
                }
            }       
        }
    }   
}
template void cpuComputeSna(double*, double*, int*, int*, int*, int, int, int, int, int);
template void cpuComputeSna(float*, float*, int*, int*, int*, int, int, int, int, int);

template <typename T> void cpuComputeSnad(T *snad, T *dblist, T *blist, int *aii, int *ai, int *aj, int *ti, 
        int *mask, int ncoeff, int ntypes, int nperdim, int ijnum, int quadraticflag)
{
    for (int k=0; k<ijnum ; k++) {  
        int i = ai[k];
        if (mask[i]) {
            int j = aj[k];
            int itype = ti[k];
            T *db = &dblist[3*ncoeff*k];                        
            T *snadi = &snad[3*nperdim*ntypes*i + 3*nperdim*(itype-1)];
            T *snadj = &snad[3*nperdim*ntypes*j + 3*nperdim*(itype-1)];
            for (int icoeff = 0; icoeff < ncoeff; icoeff++) {
              snadi[0+3*icoeff] += db[0+3*icoeff];
              snadi[1+3*icoeff] += db[1+3*icoeff];
              snadi[2+3*icoeff] += db[2+3*icoeff];
              snadj[0+3*icoeff] -= db[0+3*icoeff];
              snadj[1+3*icoeff] -= db[1+3*icoeff];
              snadj[2+3*icoeff] -= db[2+3*icoeff];
            }
        }
    }
    
    if (quadraticflag) {
        for (int k=0; k<ijnum ; k++) {  
            int i = ai[k];
            if (mask[i]) {
                int j = aj[k];
                int itype = ti[k];
                T *b = &blist[ncoeff*aii[k]];         
                T *db = &dblist[3*ncoeff*k];                        
                T *snadi = &snad[3*nperdim*ntypes*i + 3*nperdim*(itype-1)];
                T *snadj = &snad[3*nperdim*ntypes*j + 3*nperdim*(itype-1)];
                int ncount = ncoeff;
                for (int icoeff = 0; icoeff < ncoeff; icoeff++) {
                    T bi = b[icoeff];
                    T bix = db[0+3*icoeff];
                    T biy = db[1+3*icoeff];
                    T biz = db[2+3*icoeff];

                    // diagonal elements of quadratic matrix
                    T dbxtmp = bi*bix;
                    T dbytmp = bi*biy;
                    T dbztmp = bi*biz;

                    snadi[0+3*ncount] += dbxtmp;
                    snadi[1+3*ncount] += dbytmp;
                    snadi[2+3*ncount] += dbztmp;
                    snadj[0+3*ncount] -= dbxtmp;
                    snadj[1+3*ncount] -= dbytmp;
                    snadj[2+3*ncount] -= dbztmp;
                    ncount++;

                    // upper-triangular elements of quadratic matrix

                    for (int jcoeff = icoeff+1; jcoeff < ncoeff; jcoeff++) {
                        dbxtmp = bi*db[0+3*jcoeff] + bix*b[jcoeff];
                        dbytmp = bi*db[1+3*jcoeff] + biy*b[jcoeff];
                        dbztmp = bi*db[2+3*jcoeff] + biz*b[jcoeff];

                        snadi[0+3*ncount] += dbxtmp;
                        snadi[1+3*ncount] += dbytmp;
                        snadi[2+3*ncount] += dbztmp;
                        snadj[0+3*ncount] -= dbxtmp;
                        snadj[1+3*ncount] -= dbytmp;
                        snadj[2+3*ncount] -= dbztmp;
                        ncount++;              
                    }                
                }
            }    
        }
    }    
}
template void cpuComputeSnad(double*, double*, double*, int*, int*, int*, int*, int*,
        int, int, int, int, int);
template void cpuComputeSnad(float*, float*, float*, int*, int*, int*, int*, int*,
        int, int, int, int, int);

template <typename T> void cpuComputeSnad(T *snad, T *dblist, T *blist, int *aii, int *ai, int *aj, int *ti, 
        int *mask, int *tag, int ncoeff, int ntypes, int nperdim, int ijnum, int quadraticflag)
{
    for (int k=0; k<ijnum ; k++) {  
        int i = ai[k];
        int ig = tag[i]-1; // global index of atom i       
        if (mask[i]) {
            int j = aj[k];
            int jg = tag[j]-1; // global index of atom j
            int itype = ti[k];            
            T *db = &dblist[3*ncoeff*k];                        
            T *snadi = &snad[nperdim*ntypes*3*ig + nperdim*(itype-1)];
            T *snadj = &snad[nperdim*ntypes*3*jg + nperdim*(itype-1)];
            for (int icoeff = 0; icoeff < ncoeff; icoeff++) {
              snadi[icoeff] += db[0+3*icoeff];
              snadi[icoeff + nperdim*ntypes] += db[1+3*icoeff];
              snadi[icoeff + 2*nperdim*ntypes] += db[2+3*icoeff];
              snadj[icoeff] -= db[0+3*icoeff];
              snadj[icoeff + nperdim*ntypes] -= db[1+3*icoeff];
              snadj[icoeff + 2*nperdim*ntypes] -= db[2+3*icoeff];
            }
        }
    }
    
    if (quadraticflag) {
        for (int k=0; k<ijnum ; k++) {  
            int i = ai[k];
            int ig = tag[i]-1; // global index of atom i       
            if (mask[i]) {
                int j = aj[k];
                int jg = tag[j]-1; // global index of atom j
                int itype = ti[k];
                T *b = &blist[ncoeff*aii[k]];         
                T *db = &dblist[3*ncoeff*k];                        
                T *snadi = &snad[nperdim*ntypes*3*ig + nperdim*(itype-1)];
                T *snadj = &snad[nperdim*ntypes*3*jg + nperdim*(itype-1)];
                int ncount = ncoeff;
                for (int icoeff = 0; icoeff < ncoeff; icoeff++) {
                    T bi = b[icoeff];
                    T bix = db[0+3*icoeff];
                    T biy = db[1+3*icoeff];
                    T biz = db[2+3*icoeff];

                    // diagonal elements of quadratic matrix
                    T dbxtmp = bi*bix;
                    T dbytmp = bi*biy;
                    T dbztmp = bi*biz;

                    snadi[ncount]                    += dbxtmp;
                    snadi[ncount +   nperdim*ntypes] += dbytmp;
                    snadi[ncount + 2*nperdim*ntypes] += dbztmp;
                    snadj[ncount]                    -= dbxtmp;
                    snadj[ncount +   nperdim*ntypes] -= dbytmp;
                    snadj[ncount + 2*nperdim*ntypes] -= dbztmp;
                    ncount++;

                    // upper-triangular elements of quadratic matrix

                    for (int jcoeff = icoeff+1; jcoeff < ncoeff; jcoeff++) {
                        dbxtmp = bi*db[0+3*jcoeff] + bix*b[jcoeff];
                        dbytmp = bi*db[1+3*jcoeff] + biy*b[jcoeff];
                        dbztmp = bi*db[2+3*jcoeff] + biz*b[jcoeff];

                        snadi[ncount]                    += dbxtmp;
                        snadi[ncount +   nperdim*ntypes] += dbytmp;
                        snadi[ncount + 2*nperdim*ntypes] += dbztmp;
                        snadj[ncount]                    -= dbxtmp;
                        snadj[ncount +   nperdim*ntypes] -= dbytmp;
                        snadj[ncount + 2*nperdim*ntypes] -= dbztmp;                        
                        ncount++;              
                    }                
                }
            }    
        }
    }    
}
template void cpuComputeSnad(double*, double*, double*, int*, int*, int*, int*, int*, int*,
        int, int, int, int, int);
template void cpuComputeSnad(float*, float*, float*, int*, int*, int*, int*, int*, int*,
        int, int, int, int, int);

template <typename T> void cpuComputeSnav(T *snav, T *dblist, T *blist, T *x, int *aii, int *ai, int *aj,
        int *ti, int *mask, int ncoeff, int ntypes, int nperdim, int ijnum, int quadraticflag)
{
    for (int k=0; k<ijnum ; k++) {  
        int i = ai[k];        
        if (mask[i]) {
            int j = aj[k];
            int itype = ti[k];
            T *db = &dblist[3*ncoeff*k];                        
            T *snavi = &snav[6*nperdim*ntypes*i + 6*nperdim*(itype-1)];
            T *snavj = &snav[6*nperdim*ntypes*j + 6*nperdim*(itype-1)];
            T xi = x[0+3*i]; 
            T yi = x[1+3*i]; 
            T zi = x[2+3*i]; 
            T xj = x[0+3*j]; 
            T yj = x[1+3*j]; 
            T zj = x[2+3*j];             
            for (int icoeff = 0; icoeff < ncoeff; icoeff++) {
              snavi[0+6*icoeff] += db[0+3*icoeff]*xi;
              snavi[1+6*icoeff] += db[1+3*icoeff]*yi;
              snavi[2+6*icoeff] += db[2+3*icoeff]*zi;
              snavi[3+6*icoeff] += db[1+3*icoeff]*zi;
              snavi[4+6*icoeff] += db[0+3*icoeff]*zi;
              snavi[5+6*icoeff] += db[0+3*icoeff]*yi;              
              snavj[0+6*icoeff] -= db[0+3*icoeff]*xj;
              snavj[1+6*icoeff] -= db[1+3*icoeff]*yj;
              snavj[2+6*icoeff] -= db[2+3*icoeff]*zj;
              snavj[3+6*icoeff] -= db[1+3*icoeff]*zj;
              snavj[4+6*icoeff] -= db[0+3*icoeff]*zj;
              snavj[5+6*icoeff] -= db[0+3*icoeff]*yj;
            }
        }
    }
        
    if (quadraticflag) {
        for (int k=0; k<ijnum ; k++) {  
            int i = ai[k];
            if (mask[i]) {
                int j = aj[k];
                int itype = ti[k];
                T *b = &blist[ncoeff*aii[k]];         
                T *db = &dblist[3*ncoeff*k];                        
                T *snavi = &snav[6*nperdim*ntypes*i + 6*nperdim*(itype-1)];
                T *snavj = &snav[6*nperdim*ntypes*j + 6*nperdim*(itype-1)];
                T xi = x[0+3*i]; 
                T yi = x[1+3*i]; 
                T zi = x[2+3*i]; 
                T xj = x[0+3*j]; 
                T yj = x[1+3*j]; 
                T zj = x[2+3*j];                             
                int ncount = ncoeff;
                for (int icoeff = 0; icoeff < ncoeff; icoeff++) {
                    T bi = b[icoeff];
                    T bix = db[0+3*icoeff];
                    T biy = db[1+3*icoeff];
                    T biz = db[2+3*icoeff];

                    // diagonal elements of quadratic matrix
                    T dbxtmp = bi*bix;
                    T dbytmp = bi*biy;
                    T dbztmp = bi*biz;

                    snavi[0+6*ncount] += dbxtmp*xi;
                    snavi[1+6*ncount] += dbytmp*yi;
                    snavi[2+6*ncount] += dbztmp*zi;
                    snavi[3+6*ncount] += dbytmp*zi;
                    snavi[4+6*ncount] += dbxtmp*zi;
                    snavi[5+6*ncount] += dbxtmp*yi;                    
                    snavj[0+6*ncount] -= dbxtmp*xj;
                    snavj[1+6*ncount] -= dbytmp*yj;
                    snavj[2+6*ncount] -= dbztmp*zj;
                    snavj[3+6*ncount] -= dbytmp*zj;
                    snavj[4+6*ncount] -= dbxtmp*zj;
                    snavj[5+6*ncount] -= dbxtmp*yj;
                    ncount++;

                    // upper-triangular elements of quadratic matrix

                    for (int jcoeff = icoeff+1; jcoeff < ncoeff; jcoeff++) {
                        dbxtmp = bi*db[0+3*jcoeff] + bix*b[jcoeff];
                        dbytmp = bi*db[1+3*jcoeff] + biy*b[jcoeff];
                        dbztmp = bi*db[2+3*jcoeff] + biz*b[jcoeff];
                        snavi[0+6*ncount] += dbxtmp*xi;
                        snavi[1+6*ncount] += dbytmp*yi;
                        snavi[2+6*ncount] += dbztmp*zi;
                        snavi[3+6*ncount] += dbytmp*zi;
                        snavi[4+6*ncount] += dbxtmp*zi;
                        snavi[5+6*ncount] += dbxtmp*yi;                    
                        snavj[0+6*ncount] -= dbxtmp*xj;
                        snavj[1+6*ncount] -= dbytmp*yj;
                        snavj[2+6*ncount] -= dbztmp*zj;
                        snavj[3+6*ncount] -= dbytmp*zj;
                        snavj[4+6*ncount] -= dbxtmp*zj;
                        snavj[5+6*ncount] -= dbxtmp*yj;
                        ncount++;              
                    }                
                }
            }    
        }
    }    
}
template void cpuComputeSnav(double*, double*, double*, double*, int*, int*, int*, int*, int*,
        int, int, int, int, int);
template void cpuComputeSnav(float*, float*, float*, float*, int*, int*, int*, int*, int*,
        int, int, int, int, int);

template <typename T> void cpuComputeSnav2(T *snav, T *dblist, T *blist, T *x, int *aii, int *ai, int *aj,
        int *ti, int *mask, int ncoeff, int ntypes, int nperdim, int ijnum, int quadraticflag)
{
    for (int k=0; k<ijnum ; k++) {  
        int i = ai[k];        
        if (mask[i]) {
            int j = aj[k];
            int itype = ti[k];
            T *db = &dblist[3*ncoeff*k];                        
            T *snavi = &snav[6*nperdim*ntypes*i + nperdim*(itype-1)];
            T *snavj = &snav[6*nperdim*ntypes*j + nperdim*(itype-1)];
            T xi = x[0+3*i]; 
            T yi = x[1+3*i]; 
            T zi = x[2+3*i]; 
            T xj = x[0+3*j]; 
            T yj = x[1+3*j]; 
            T zj = x[2+3*j];             
            for (int icoeff = 0; icoeff < ncoeff; icoeff++) {
              snavi[icoeff] += db[0+3*icoeff]*xi;
              snavi[icoeff + nperdim*ntypes] += db[1+3*icoeff]*yi;
              snavi[icoeff + 2*nperdim*ntypes] += db[2+3*icoeff]*zi;
              snavi[icoeff + 3*nperdim*ntypes] += db[1+3*icoeff]*zi;
              snavi[icoeff + 4*nperdim*ntypes] += db[0+3*icoeff]*zi;
              snavi[icoeff + 5*nperdim*ntypes] += db[0+3*icoeff]*yi;              
              snavj[icoeff] -= db[0+3*icoeff]*xj;
              snavj[icoeff + nperdim*ntypes] -= db[1+3*icoeff]*yj;
              snavj[icoeff + 2*nperdim*ntypes] -= db[2+3*icoeff]*zj;
              snavj[icoeff + 3*nperdim*ntypes] -= db[1+3*icoeff]*zj;
              snavj[icoeff + 4*nperdim*ntypes] -= db[0+3*icoeff]*zj;
              snavj[icoeff + 5*nperdim*ntypes] -= db[0+3*icoeff]*yj;
            }
        }
    }
        
    if (quadraticflag) {
        for (int k=0; k<ijnum ; k++) {  
            int i = ai[k];
            if (mask[i]) {
                int j = aj[k];
                int itype = ti[k];
                T *b = &blist[ncoeff*aii[k]];         
                T *db = &dblist[3*ncoeff*k];                        
                T *snavi = &snav[6*nperdim*ntypes*i + nperdim*(itype-1)];
                T *snavj = &snav[6*nperdim*ntypes*j + nperdim*(itype-1)];
                T xi = x[0+3*i]; 
                T yi = x[1+3*i]; 
                T zi = x[2+3*i]; 
                T xj = x[0+3*j]; 
                T yj = x[1+3*j]; 
                T zj = x[2+3*j];                             
                int ncount = ncoeff;
                for (int icoeff = 0; icoeff < ncoeff; icoeff++) {
                    T bi = b[icoeff];
                    T bix = db[0+3*icoeff];
                    T biy = db[1+3*icoeff];
                    T biz = db[2+3*icoeff];

                    // diagonal elements of quadratic matrix
                    T dbxtmp = bi*bix;
                    T dbytmp = bi*biy;
                    T dbztmp = bi*biz;

                    snavi[ncount] += dbxtmp*xi;
                    snavi[ncount + nperdim*ntypes] += dbytmp*yi;
                    snavi[ncount + 2*nperdim*ntypes] += dbztmp*zi;
                    snavi[ncount + 3*nperdim*ntypes] += dbytmp*zi;
                    snavi[ncount + 4*nperdim*ntypes] += dbxtmp*zi;
                    snavi[ncount + 5*nperdim*ntypes] += dbxtmp*yi;                    
                    snavj[ncount] -= dbxtmp*xj;
                    snavj[ncount + nperdim*ntypes] -= dbytmp*yj;
                    snavj[ncount + 2*nperdim*ntypes] -= dbztmp*zj;
                    snavj[ncount + 3*nperdim*ntypes] -= dbytmp*zj;
                    snavj[ncount + 4*nperdim*ntypes] -= dbxtmp*zj;
                    snavj[ncount + 5*nperdim*ntypes] -= dbxtmp*yj;
                    ncount++;

                    // upper-triangular elements of quadratic matrix

                    for (int jcoeff = icoeff+1; jcoeff < ncoeff; jcoeff++) {
                        dbxtmp = bi*db[0+3*jcoeff] + bix*b[jcoeff];
                        dbytmp = bi*db[1+3*jcoeff] + biy*b[jcoeff];
                        dbztmp = bi*db[2+3*jcoeff] + biz*b[jcoeff];
                        snavi[ncount] += dbxtmp*xi;
                        snavi[ncount + nperdim*ntypes] += dbytmp*yi;
                        snavi[ncount + 2*nperdim*ntypes] += dbztmp*zi;
                        snavi[ncount + 3*nperdim*ntypes] += dbytmp*zi;
                        snavi[ncount + 4*nperdim*ntypes] += dbxtmp*zi;
                        snavi[ncount + 5*nperdim*ntypes] += dbxtmp*yi;                    
                        snavj[ncount] -= dbxtmp*xj;
                        snavj[ncount + nperdim*ntypes] -= dbytmp*yj;
                        snavj[ncount + 2*nperdim*ntypes] -= dbztmp*zj;
                        snavj[ncount + 3*nperdim*ntypes] -= dbytmp*zj;
                        snavj[ncount + 4*nperdim*ntypes] -= dbxtmp*zj;
                        snavj[ncount + 5*nperdim*ntypes] -= dbxtmp*yj;
                        ncount++;              
                    }                
                }
            }    
        }
    }    
}
template void cpuComputeSnav2(double*, double*, double*, double*, int*, int*, int*, int*, int*,
        int, int, int, int, int);
template void cpuComputeSnav2(float*, float*, float*, float*, int*, int*, int*, int*, int*,
        int, int, int, int, int);

template <typename T> void cpuComputeBeta2(T *beta, T *bispectrum, T *coeffelem, int *ilist, int *map, int *type, 
        int inum, int ncoeff, int ncoeffall, int quadraticflag)
{
  for (int ii = 0; ii < inum; ii++) { // for each atom in the list
    int i = ilist[ii]; // index of atom i
    const int itype = type[i]; // element type of atom i
    const int ielem = map[itype];  // index of that element type
    T* coeffi = &coeffelem[ielem*ncoeffall]; // coefficient for that particular element
    
    for (int icoeff = 0; icoeff < ncoeff; icoeff++) // for each coefficient 
      beta[ii*ncoeff + icoeff] = coeffi[icoeff+1]; // beta is the same for all atoms of the same element type

    if (quadraticflag) {
      int k = ncoeff+1;
      for (int icoeff = 0; icoeff < ncoeff; icoeff++) {
        T bveci = bispectrum[ii*ncoeff + icoeff];
        beta[ii*ncoeff + icoeff] += coeffi[k]*bveci;
        k++;
        for (int jcoeff = icoeff+1; jcoeff < ncoeff; jcoeff++) {
          T bvecj = bispectrum[ii*ncoeff + jcoeff];
          beta[ii*ncoeff + icoeff] += coeffi[k]*bvecj;
          beta[ii*ncoeff + jcoeff] += coeffi[k]*bveci;
          k++;
        }
      }
    }
  }
}
template void cpuComputeBeta2(double*, double*, double*, int*, int*, int*,
        int, int, int, int);
template void cpuComputeBeta2(float*, float*, float*, int*, int*, int*,
        int, int, int, int);

template <typename T> void cpuSnapTallyEnergyFull(T *eatom, T *bispectrum, T *coeffelem, int *ilist, int *map, int *type, 
        int inum, int ncoeff, int ncoeffall, int quadraticflag)
{  
    for (int ii=0; ii<inum; ii++) {
        int i = ilist[ii]; // index of atom i
        const int itype = type[i]; // element type of atom i
        const int ielem = map[itype];  // index of that element type
        T* coeffi = &coeffelem[ielem*ncoeffall]; // coefficient for that particular element
        
        T evdwl = coeffi[0];
        for (int icoeff = 0; icoeff < ncoeff; icoeff++)
            evdwl += coeffi[icoeff+1]*bispectrum[ii*ncoeff + icoeff];
     
        if (quadraticflag) {
            int k = ncoeff+1;
            for (int icoeff = 0; icoeff < ncoeff; icoeff++) {
                T bveci = bispectrum[ii*ncoeff + icoeff];
                evdwl += 0.5*coeffi[k++]*bveci*bveci;
                for (int jcoeff = icoeff+1; jcoeff < ncoeff; jcoeff++) {
                    T bvecj = bispectrum[ii*ncoeff + jcoeff];
                    evdwl += coeffi[k++]*bveci*bvecj;
                }
            }
        }                
        eatom[ii] = evdwl;
    }
}
template void cpuSnapTallyEnergyFull(double*, double*, double*, int*, int*, int*,
        int, int, int, int);
template void cpuSnapTallyEnergyFull(float*, float*, float*, int*, int*, int*,
        int, int, int, int);

template <typename T> void cpuSnapTallyForceFull(T *fatom, T *fij, int *ai, int *aj, int ijnum)
{ 
    for(int k = 0; k < ijnum; k++) {
        int i = ai[k];        
        int j = aj[k];        
        T fx = fij[0+3*k];
        T fy = fij[1+3*k];
        T fz = fij[2+3*k];    
        fatom[0+3*i] += fx;
        fatom[1+3*i] += fy;
        fatom[2+3*i] += fz;
        fatom[0+3*j] -= fx;
        fatom[1+3*j] -= fy;
        fatom[2+3*j] -= fz;
    }
};
template void cpuSnapTallyForceFull(double*, double*, int*, int*, int);
template void cpuSnapTallyForceFull(float*, float*, int*, int*, int);

template <typename T> void cpuSnapTallyVirialFull(T *vatom, T *fij, T *rij, int *ai, int *aj, int ijnum)
{ 
    for(int k = 0; k < ijnum; k++) {
        int i = ai[k];        
        int j = aj[k];        
        T dx = -rij[0+3*k];
        T dy = -rij[1+3*k];
        T dz = -rij[2+3*k];    
        T fx = fij[0+3*k];
        T fy = fij[1+3*k];
        T fz = fij[2+3*k];    
        T v0 = dx*fx;
        T v1 = dy*fy;
        T v2 = dz*fz;
        T v3 = dx*fy;
        T v4 = dx*fz;
        T v5 = dy*fz;        
        vatom[0+6*i] += v0;
        vatom[1+6*i] += v1;
        vatom[2+6*i] += v2;
        vatom[3+6*i] += v3;
        vatom[4+6*i] += v4;
        vatom[5+6*i] += v5;        
        vatom[0+6*j] += v0;
        vatom[1+6*j] += v1;
        vatom[2+6*j] += v2;
        vatom[3+6*j] += v3;
        vatom[4+6*j] += v4;
        vatom[5+6*j] += v5;        
    }
};
template void cpuSnapTallyVirialFull(double*, double*, double*, int*, int*, int);
template void cpuSnapTallyVirialFull(float*, float*, float*, int*, int*, int);

template <typename T> void cpuNeighPairList(int *pairnum, int *pairlist, T *x, T rcutsq, int *ilist, int *neighlist, 
        int *neighnum, int inum, int jnum, int dim)
{    
    for (int ii=0; ii<inum; ii++) {
        int i = ilist[ii];       // atom i
        int m = neighnum[i];     // number of neighbors around i                       
        int count = 0;              
        for (int l=0; l<m ; l++) {   // loop over each atom around atom i {
            int j = neighlist[l + jnum*i];            
            // distance between atom i and atom j                                    
            T xij0 = x[j*dim] - x[i*dim];  // xj - xi
            T xij1 = x[j*dim+1] - x[i*dim+1]; // xj - xi               
            T xij2 = x[j*dim+2] - x[i*dim+2]; // xj - xi               
            T dij = xij0*xij0 + xij1*xij1 + xij2*xij2;                        
            if (dij <= rcutsq) {
                pairlist[count + jnum*ii] = j;  // atom j     
                count += 1;
            }
        }        
        pairnum[ii] = count;       
    }    
};
template void cpuNeighPairList(int*, int*, double*, double, int*, int*, int*, int, int, int);
template void cpuNeighPairList(int*, int*, float*, float, int*, int*, int*, int, int, int);

template <typename T> void cpuNeighPairList(int *pairnum, int *pairlist, T *x, T *rcutsq, int *ilist, int *neighlist, 
        int *neighnum, int inum, int jnum, int dim)
{    
    for (int ii=0; ii<inum; ii++) {
        int i = ilist[ii];       // atom i
        int m = neighnum[i];     // number of neighbors around i                       
        int count = 0;              
        for (int l=0; l<m ; l++) {   // loop over each atom around atom i {
            int j = neighlist[l + jnum*i];            
            // distance between atom i and atom j                                    
            T xij0 = x[j*dim] - x[i*dim];  // xj - xi
            T xij1 = x[j*dim+1] - x[i*dim+1]; // xj - xi               
            T xij2 = x[j*dim+2] - x[i*dim+2]; // xj - xi               
            T dij = xij0*xij0 + xij1*xij1 + xij2*xij2;                        
            if (dij <= rcutsq[0]) {
                pairlist[count + jnum*ii] = j;  // atom j     
                count += 1;
            }
        }        
        pairnum[ii] = count;       
    }    
};
template void cpuNeighPairList(int*, int*, double*, double*, int*, int*, int*, int, int, int);
template void cpuNeighPairList(int*, int*, float*, float*, int*, int*, int*, int, int, int);

template <typename T> void cpuNeighPairList(int *pairnum, int *pairlist, T *x, T *rcutsq, int *ilist, int *neighlist, 
        int *neighnum, int *atomtype, int inum, int jnum, int dim, int ntypes)
{    
    for (int ii=0; ii<inum; ii++) {
        int i = ilist[ii];       // atom i
        int itype = atomtype[i];
        int m = neighnum[i];     // number of neighbors around i                       
        int count = 0;              
        for (int l=0; l<m ; l++) {   // loop over each atom around atom i {
            int j = neighlist[l + jnum*i];         
            int jtype  = atomtype[j];        
            // distance between atom i and atom j                                    
            T xij0 = x[j*dim] - x[i*dim];  // xj - xi
            T xij1 = x[j*dim+1] - x[i*dim+1]; // xj - xi               
            T xij2 = x[j*dim+2] - x[i*dim+2]; // xj - xi               
            T dij = xij0*xij0 + xij1*xij1 + xij2*xij2;                        
            if (dij <= rcutsq[(jtype-1) + (itype-1)*ntypes]) {
                pairlist[count + jnum*ii] = j;  // atom j     
                count += 1;
            }
        }        
        pairnum[ii] = count;       
    }    
};
template void cpuNeighPairList(int*, int*, double*, double*, int*, int*, int*, int*, 
        int, int, int, int);
template void cpuNeighPairList(int*, int*, float*, float*, int*, int*, int*, int*, 
        int, int, int, int);

template <typename T> void cpuNeighPairs(T *xij, T *x, int *aii, int *ai, int *aj,  
      int *ti, int *tj, int *pairnum, int *pairlist, int *pairnumsum, int *ilist, 
      int *atomtype, int inum, int jnum, int dim)
{        
    for (int ii=0; ii<inum; ii++) {  // for each atom i in the simulation box     
        int i = ilist[ii];       // atom i
        int itype = atomtype[i];
        int m = pairnum[ii];        // number of neighbors around i             
        int start = pairnumsum[ii];   
        for (int l=0; l<m ; l++) {   // loop over each atom around atom i
            int j = pairlist[l + jnum*ii];  // atom j              
            int k = start + l;                                     
            aii[k]       = ii;
            ai[k]        = i;
            aj[k]        = j;          
            ti[k]        = itype;       
            tj[k]        = atomtype[j];        
            for (int d=0; d<dim; d++) 
                xij[k*dim+d]   = x[j*dim+d] -  x[i*dim+d];  // xj - xi            
        }
    }    
};
template void cpuNeighPairs(double*, double*, int*, int*, int*, int*, int*, int*, int*, int*, 
        int*, int*, int, int, int);
template void cpuNeighPairs(float*, float*, int*, int*, int*, int*, int*, int*, int*, int*, 
        int*, int*, int, int, int);

#endif

