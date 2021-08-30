/***************************************************************************                               
                    Molecular Dynamics Potentials (MDP)
                           CESMIX-MIT Project  
 
 Contributing authors: Ngoc-Cuong Nguyen (cuongng@mit.edu, exapde@gmail.com)
 ***************************************************************************/

#ifndef MDP_GPUSNAP
#define MDP_GPUSNAP

#define MIN(a,b) ((a) < (b) ? (a) : (b))
#define MAX(a,b) ((a) > (b) ? (a) : (b))

void gpuBuildIndexList(int *idx_max, int *idxz, int *idxz_block, int *idxb, int *idxb_block, int *idxu_block, 
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

template <typename T> void gpuInitRootpqArray(T *rootpqarray, int twojmax)
{
  int jdim = twojmax + 1;
  for (int p = 1; p <= twojmax; p++)
    for (int q = 1; q <= twojmax; q++)
      rootpqarray[p*jdim + q] = sqrt(static_cast<T>(p)/q);
};
template void gpuInitRootpqArray(double*, int);
template void gpuInitRootpqArray(float*, int);

template <typename T> T gpuDeltacg(T *factorial, int j1, int j2, int j)
{
  T sfaccg = factorial[(j1 + j2 + j) / 2 + 1];
  return sqrt(factorial[(j1 + j2 - j) / 2] *
              factorial[(j1 - j2 + j) / 2] *
              factorial[(-j1 + j2 + j) / 2] / sfaccg);
};
template double gpuDeltacg(double*, int, int, int);
template float gpuDeltacg(float*, int, int, int);

template <typename T> void gpuInitClebschGordan(T *cglist, T *factorial, int twojmax)
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
            dcg = gpuDeltacg(factorial, j1, j2, j);
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
template void gpuInitClebschGordan(double*, double*, int);
template void gpuInitClebschGordan(float*, float*, int);

template <typename T> void gpuInitSna(T *rootpqarray, T *cglist, T *factorial, int *idx_max, int *idxz, 
      int *idxz_block, int *idxb, int *idxb_block, int *idxu_block, int *idxcg_block, int twojmax)
{
    gpuBuildIndexList(idx_max, idxz, idxz_block, idxb, 
            idxb_block, idxu_block, idxcg_block, twojmax);
    
    gpuInitRootpqArray(rootpqarray, twojmax);
    gpuInitClebschGordan(cglist, factorial, twojmax);        
}
template void gpuInitSna(double*, double*, double *, int*, int*, int*, int*, int*, int*, int*, int);
template void gpuInitSna(float*, float*, float*, int*, int*, int*, int*, int*, int*, int*, int);

template <typename T> __device__ T gpuComputeSfac(T r, T rcut, T rmin0, int switch_flag)
{
  if (switch_flag == 0) return (T) 1.0;
  if (switch_flag == 1) {
    if(r <= rmin0) return (T) 1.0;
    else if(r > rcut) return (T) 0.0;
    else {
      T rcutfac = M_PI / (rcut - rmin0);
      return 0.5 * (cos((r - rmin0) * rcutfac) + 1.0);
    }
  }
  return (T) 0.0;
};
// template __device__ double gpuComputeSfac(double, double, double, int);
// template __device__ float gpuComputeSfac(float, float, float, int);

template <typename T> __device__ T gpuComputeDsfac(T r, T rcut, T rmin0, int switch_flag)
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
// template __device__ double gpuComputeDsfac(double, double, double, int);
// template __device__ float gpuComputeDsfac(float, float, float, int);

template <typename T> __global__ void gpuKernelZeroUarraytot(T *ulisttot_r, T *ulisttot_i, T wself, int *idxu_block, 
        int *type, int *map, int *ai, int wselfall_flag, int chemflag, int idxu_max, int nelements,
        int N1, int N2, int N3)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;
    while (idx < N3) {
        int l = idx%N2;        
        int j = l%N1;
        int jelem = (l-j)/N1;
        int ii = (idx-l)/N2;    
        int ielem = (chemflag) ? map[type[ai[ii]]]: 0;                
        int jju = idxu_block[j];        
        for (int mb = 0; mb <= j; mb++) {
            for (int ma = 0; ma <= j; ma++) {
                ulisttot_r[ii*nelements*idxu_max+jelem*idxu_max+jju] = 0.0;
                ulisttot_i[ii*nelements*idxu_max+jelem*idxu_max+jju] = 0.0;                                
                
                if (jelem == ielem || wselfall_flag)
                    if (ma==mb)
                        ulisttot_r[ii*nelements*idxu_max+jelem*idxu_max+jju] = wself; ///// T check this
                jju++;
            }
        }
        idx += blockDim.x * gridDim.x;
    }                    
};
template <typename T> void gpuZeroUarraytot(T *ulisttot_r, T *ulisttot_i, T wself, int *idxu_block, int *type,
        int *map, int *ai, int wselfall_flag, int chemflag, int idxu_max, int nelements, 
                int twojmax, int inum)
{
    int N1 = twojmax+1;
    int N2 = N1*nelements;
    int N3 = N2*inum;                        
    int blockDim = 256;
    int gridDim = (N3 + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;
    gpuKernelZeroUarraytot<<<gridDim, blockDim>>>(ulisttot_r, ulisttot_i, wself, idxu_block, 
            type, map, ai, wselfall_flag, chemflag, idxu_max, nelements, N1, N2, N3);            
}
template void gpuZeroUarraytot(double*, double*, double, int*, int*, 
        int*, int*, int, int, int, int, int, int);
template void gpuZeroUarraytot(float*, float*, float, int*, int*, 
        int*, int*, int, int, int, int, int, int);

// template <typename T> __global__ void gpuKernelZeroUarraytot2(T *ulisttot_r, T *ulisttot_i, T wself, int *idxu_block, 
//         int *type, int *map, int *ai, int wselfall_flag, int chemflag, int inum, int idxu_max, int nelements,
//         int N1, int N2, int N3)
// {
//     int idx = threadIdx.x + blockIdx.x * blockDim.x;
//     while (idx < N3) {
//         int l = idx%N2;  // inum*(twojmax+1)
//         int ii = l%N1;    // inum
//         int j = (l-ii)/N1; // (twojmax+1)
//         int jelem = (idx-l)/N2; // nelements   
//         int ielem = (chemflag) ? map[type[ai[ii]]]: 0;                
//         int jju = idxu_block[j];        
//         for (int mb = 0; mb <= j; mb++) {
//             for (int ma = 0; ma <= j; ma++) {
//                 int n = ii + inum*jju + inum*idxu_max*jelem;        
//                 ulisttot_r[n] = 0.0;
//                 if (jelem == ielem || wselfall_flag)
//                     if (ma==mb)
//                         ulisttot_r[n] = wself; ///// T check this
//                 jju++;
//             }
//         }
//         idx += blockDim.x * gridDim.x;
//     }                    
// };
// template <typename T> void gpuZeroUarraytot2(T *ulisttot_r, T *ulisttot_i, T wself, int *idxu_block, int *type,
//         int *map, int *ai, int wselfall_flag, int chemflag, int idxu_max, int nelements, 
//                 int twojmax, int inum)
// {
//     int N1 = inum;
//     int N2 = N1*(twojmax+1);
//     int N3 = N2*nelements;                                
//     int blockDim = 256;
//     int gridDim = (N3 + blockDim - 1) / blockDim;
//     gridDim = (gridDim>1024)? 1024 : gridDim;
//     gpuKernelZeroUarraytot2<<<gridDim, blockDim>>>(ulisttot_r, ulisttot_i, wself, idxu_block, 
//             type, map, ai, wselfall_flag, chemflag, inum, idxu_max, nelements, N1, N2, N3);            
// }
// template void gpuZeroUarraytot2(double*, double*, double, int*, int*, 
//         int*, int*, int, int, int, int, int, int);
// template void gpuZeroUarraytot2(float*, float*, float, int*, int*, 
//         int*, int*, int, int, int, int, int, int);
// 
template <typename T> __global__ void gpuKernelAddUarraytot(T *ulisttot_r, T *ulisttot_i, T *ulist_r, T *ulist_i, 
        T *rij, T *wjelem, T *radelem, T rmin0, T rcutfac, int *idxu_block, int *ilist, 
        int *type, int *pairnum, int *pairnumsum, int *map, int *tj, int twojmax, int idxu_max, 
        int nelements, int inum, int switch_flag, int chemflag)
{
    int ii = threadIdx.x + blockIdx.x * blockDim.x;    
    while (ii < inum) {     
    //for (int ii=0; ii<inum; ii++) {  // for each atom i in the simulation box     
        int i = ilist[ii];       // atom i
        int itype = type[i];
        int m = pairnum[ii];        // number of neighbors around i             
        int start = pairnumsum[ii];   
        for (int l=0; l<m ; l++) {   // loop over each atom around atom i            
            int k = start + l;    
            int jtype = tj[k];
            int jelem = (chemflag) ? map[jtype] : 0;                
            T x = rij[k*3+0];
            T y = rij[k*3+1];
            T z = rij[k*3+2];
            T rsq = x * x + y * y + z * z;
            T r = sqrt(rsq);
            T rcut = (radelem[itype]+radelem[jtype])*rcutfac;
            T sfac = gpuComputeSfac(r, rcut, rmin0, switch_flag);
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
        ii += blockDim.x * gridDim.x;
    }      
};
template <typename T> void gpuAddUarraytot(T *ulisttot_r, T *ulisttot_i, T *ulist_r, T *ulist_i, 
        T *rij, T *wjelem, T *radelem, T rmin0, T rcutfac, int *idxu_block, int *ilist, 
        int *type, int *pairnum, int *pairnumsum, int *map, int *tj, int twojmax, int idxu_max, 
        int nelements, int inum, int switch_flag, int chemflag)
{
    int blockDim = 256;
    int gridDim = (inum + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;
    gpuKernelAddUarraytot<<<gridDim, blockDim>>>(ulisttot_r, ulisttot_i, ulist_r, ulist_i, rij, wjelem, radelem, rmin0, 
                rcutfac, idxu_block, ilist, type, pairnum, pairnumsum, map, tj, twojmax, 
                idxu_max, nelements, inum, switch_flag, chemflag);  
}
template void gpuAddUarraytot(double*, double*, double*, double*, double*, double*, double*, 
        double, double, int*, int*, int*, int*, int*, int*, int*, int, int, int, int, int, int);
template void gpuAddUarraytot(float*, float*, float*, float*, float*, float*, float*, float,
        float, int*, int*, int*, int*, int*, int*, int*, int, int, int, int, int, int);

template <typename T> __global__ void gpuKernelAddUarraytot2(T *ulisttot_r, T *ulisttot_i, T *ulist_r, T *ulist_i, 
        T *rij, T *wjelem, T *radelem, T rmin0, T rcutfac, int *idxu_block, int *ilist, 
        int *type, int *pairnum, int *pairnumsum, int *map, int *tj, int twojmax, int idxu_max, 
        int nelements, int inum, int ijnum, int N2, int switch_flag, int chemflag)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;
    while (idx < N2) {
        int ii = idx%inum;  // inum
        int j = (idx-ii)/inum;    // (twojmax+1)        
        int i = ilist[ii];       // atom i
        int itype = type[i];
        int m = pairnum[ii];        // number of neighbors around i             
        int start = pairnumsum[ii];           
        for (int l=0; l<m; l++) {   // loop over each atom around atom i    
            int k = start + l;    
            int jtype = tj[k];
            int jelem = (chemflag) ? map[jtype] : 0;                            
            T x = rij[k*3+0];
            T y = rij[k*3+1];
            T z = rij[k*3+2];
            T rsq = x * x + y * y + z * z;
            T r = sqrt(rsq);
            T rcut = (radelem[itype]+radelem[jtype])*rcutfac;
            //T sfac = gpuComputeSfac(r, rcut, rmin0, switch_flag);
            T sfac = 0.0;        
            if (switch_flag == 0) 
                sfac = 1.0;
            else if (switch_flag == 1) {
                if (r <= rmin0) sfac = 1.0;
                else if(r > rcut) sfac = 1.0;
                else sfac = 0.5 * (cos((r - rmin0) * M_PI / (rcut - rmin0)) + 1.0);                
            } 
            sfac *= wjelem[jtype];
            
            int jju = idxu_block[j];
            for (int mb = 0; mb <= j; mb++)
                for (int ma = 0; ma <= j; ma++) {
                      int kl = ii + inum*jju + inum*idxu_max*jelem;        
                      int kr = k + ijnum*jju;        
                      ulisttot_r[kl] += sfac * ulist_r[kr];
                      ulisttot_i[kl] += sfac * ulist_i[kr];
                      jju++;
                }
        }
        idx += blockDim.x * gridDim.x;
    }
};
template <typename T> void gpuAddUarraytot2(T *ulisttot_r, T *ulisttot_i, T *ulist_r, T *ulist_i, 
        T *rij, T *wjelem, T *radelem, T rmin0, T rcutfac, int *idxu_block, int *ilist, 
        int *type, int *pairnum, int *pairnumsum, int *map, int *tj, int twojmax, int idxu_max, 
        int nelements, int inum, int ijnum, int switch_flag, int chemflag)
{    
    int N2 = inum*(twojmax+1);
    int blockDim = 256;
    int gridDim = (N2 + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;
    gpuKernelAddUarraytot2<<<gridDim, blockDim>>>(ulisttot_r, ulisttot_i, ulist_r, ulist_i, rij, wjelem, radelem, rmin0, 
                rcutfac, idxu_block, ilist, type, pairnum, pairnumsum, map, tj, twojmax, 
                idxu_max, nelements, inum, ijnum, N2, switch_flag, chemflag);  
}
template void gpuAddUarraytot2(double*, double*, double*, double*, double*, double*, double*, 
        double, double, int*, int*, int*, int*, int*, int*, int*, int, int, int, int, int, int, int);
template void gpuAddUarraytot2(float*, float*, float*, float*, float*, float*, float*, float,
        float, int*, int*, int*, int*, int*, int*, int*, int, int, int, int, int, int, int);

template <typename T> __global__ void gpuKernelAddUarraytot3(T *ulisttot_r, T *ulisttot_i, T *ulist_r, T *ulist_i, 
        T *rij, T *wjelem, T *radelem, T rmin0, T rcutfac, int *idxu_block, int *ilist, 
        int *type, int *pairnum, int *pairnumsum, int *map, int *tj, int twojmax, int idxu_max, 
        int nelements, int inum, int ijnum, int N2, int switch_flag, int chemflag)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;
    while (idx < N2) {
        int ii = idx%inum;  // inum
        int jju = (idx-ii)/inum;    // idxu_max
        int i = ilist[ii];       // atom i
        int itype = type[i];
        int m = pairnum[ii];        // number of neighbors around i             
        int start = pairnumsum[ii];           
        for (int l=0; l<m; l++) {   // loop over each atom around atom i    
            int k = start + l;    
            int jtype = tj[k];
            int jelem = (chemflag) ? map[jtype] : 0;                            
            T x = rij[k*3+0];
            T y = rij[k*3+1];
            T z = rij[k*3+2];
            T rsq = x * x + y * y + z * z;
            T r = sqrt(rsq);
            T rcut = (radelem[itype]+radelem[jtype])*rcutfac;
            T sfac = gpuComputeSfac(r, rcut, rmin0, switch_flag);
            sfac *= wjelem[jtype];
            
            int kl = ii + inum*jju + inum*idxu_max*jelem;        
            int kr = k + ijnum*jju;        
            ulisttot_r[kl] += sfac * ulist_r[kr];
            ulisttot_i[kl] += sfac * ulist_i[kr];                
        }
        idx += blockDim.x * gridDim.x;
    }
};
template <typename T> void gpuAddUarraytot3(T *ulisttot_r, T *ulisttot_i, T *ulist_r, T *ulist_i, 
        T *rij, T *wjelem, T *radelem, T rmin0, T rcutfac, int *idxu_block, int *ilist, 
        int *type, int *pairnum, int *pairnumsum, int *map, int *tj, int twojmax, int idxu_max, 
        int nelements, int inum, int ijnum, int switch_flag, int chemflag)
{    
    int N2 = inum*idxu_max;
    int blockDim = 256;
    int gridDim = (N2 + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;
    gpuKernelAddUarraytot3<<<gridDim, blockDim>>>(ulisttot_r, ulisttot_i, ulist_r, ulist_i, rij, wjelem, radelem, rmin0, 
                rcutfac, idxu_block, ilist, type, pairnum, pairnumsum, map, tj, twojmax, 
                idxu_max, nelements, inum, ijnum, N2, switch_flag, chemflag);  
}
template void gpuAddUarraytot3(double*, double*, double*, double*, double*, double*, double*, 
        double, double, int*, int*, int*, int*, int*, int*, int*, int, int, int, int, int, int, int);
template void gpuAddUarraytot3(float*, float*, float*, float*, float*, float*, float*, float,
        float, int*, int*, int*, int*, int*, int*, int*, int, int, int, int, int, int, int);

template <typename T> __device__ void gpuComputeUarray(T *ulist_r, T *ulist_i, T *rootpqarray, T x, T y, T z,
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
  }
};
template void gpuComputeUarray(double*, double*, double*, double, double, double, double, double,
        int*, int);
template void gpuComputeUarray(float*, float*, float*, float, float, float, float, float,
        int*, int);
        
template <typename T> __global__ void gpuKernelComputeUij(T *ulist_r, T *ulist_i, T *rootpqarray, 
        T *rij, T *radelem, T rmin0, T rfac0, T rcutfac, int *idxu_block, 
        int *ai, int *aj, int *ti, int *tj, int twojmax, int idxu_max, int ijnum)
{ 
  int ij = threadIdx.x + blockIdx.x * blockDim.x;     
  while (ij < ijnum) {        
  //for(int ij = 0; ij < ijnum; ij++) {
    T x = rij[ij*3+0];
    T y = rij[ij*3+1];
    T z = rij[ij*3+2];
    T rsq = x * x + y * y + z * z;
    T r = sqrt(rsq);

    //int ii = ai[ij];
    //int jj = alist[aj[ij]];
    int nij = idxu_max*ij;    
    T rcutij = (radelem[ti[ij]]+radelem[tj[ij]])*rcutfac; //(radelem[type[ii]]+radelem[type[jj]])*rcutfac;
    T theta0 = (r - rmin0) * rfac0 * M_PI / (rcutij - rmin0);
    //    theta0 = (r - rmin0) * rscale0;
    T z0 = r / tan(theta0);    
            
//     gpuComputeUarray(&ulist_r[idxu_max*ij], &ulist_i[idxu_max*ij], rootpqarray, 
//             x, y, z, z0, r, idxu_block, twojmax);    

    T r0inv;
    T a_r, b_r, a_i, b_i;
    T rootpq;
    int jdim = twojmax + 1;
  
    r0inv = 1.0 / sqrt(r * r + z0 * z0);
    a_r = r0inv * z0;
    a_i = -r0inv * z;
    b_r = r0inv * y;
    b_i = -r0inv * x;

    ulist_r[nij+0] = 1.0;
    ulist_i[nij+0] = 0.0;

    for (int j = 1; j <= twojmax; j++) {
        int jju = idxu_block[j];
        int jjup = idxu_block[j-1];
        
        // fill in left side of matrix layer from previous layer
        for (int mb = 0; 2*mb <= j; mb++) {
            ulist_r[nij+jju] = 0.0;
            ulist_i[nij+jju] = 0.0;
            for (int ma = 0; ma < j; ma++) {
                rootpq = rootpqarray[(j - ma)*jdim + (j - mb)];
                ulist_r[nij+jju] += rootpq * (a_r * ulist_r[nij+jjup] + a_i * ulist_i[nij+jjup]);
                ulist_i[nij+jju] += rootpq * (a_r * ulist_i[nij+jjup] - a_i * ulist_r[nij+jjup]);

                rootpq = rootpqarray[(ma + 1)*jdim + (j - mb)];
                ulist_r[nij+jju+1] = -rootpq * (b_r * ulist_r[nij+jjup] + b_i * ulist_i[nij+jjup]);
                ulist_i[nij+jju+1] = -rootpq * (b_r * ulist_i[nij+jjup] - b_i * ulist_r[nij+jjup]);
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
                    ulist_r[nij+jjup] = ulist_r[nij+jju];
                    ulist_i[nij+jjup] = -ulist_i[nij+jju];
                } else {
                    ulist_r[nij+jjup] = -ulist_r[nij+jju];
                    ulist_i[nij+jjup] = ulist_i[nij+jju];
                }
                mapar = -mapar;
                jju++;
                jjup--;
            }
            mbpar = -mbpar;
        }
    }        
    ij += blockDim.x * gridDim.x;
  }
};
template <typename T> void gpuComputeUij(T *ulist_r, T *ulist_i, T *rootpqarray, T *rij, 
        T *radelem, T rmin0, T rfac0, T rcutfac, int *idxu_block, 
        int *ai, int *aj, int *ti, int *tj, int twojmax, int idxu_max, int ijnum)
{
    int blockDim = 256;
    int gridDim = (ijnum + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;
    gpuKernelComputeUij<<<gridDim, blockDim>>>(ulist_r, ulist_i, rootpqarray, rij, radelem, 
            rmin0, rfac0, rcutfac, idxu_block, ai, aj, ti, tj, twojmax, idxu_max, ijnum);
}
template void gpuComputeUij(double*, double*, double*, double*, double*, double, 
        double, double, int*, int*, int*, int*, int*, int, int, int);
template void gpuComputeUij(float*, float*, float*, float*, float*, float,  float,
        float, int*, int*, int*, int*, int*, int, int, int);
        
template <typename T> __global__ void gpuKernelComputeUij2(T *ulist_r, T *ulist_i, T *rootpqarray, 
        T *rij, T *radelem, T rmin0, T rfac0, T rcutfac, int *idxu_block, 
        int *ai, int *aj, int *ti, int *tj, int twojmax, int idxu_max, int ijnum)
{ 
  int ij = threadIdx.x + blockIdx.x * blockDim.x;     
  while (ij < ijnum) {        
  //for(int ij = 0; ij < ijnum; ij++) {
    T x = rij[ij*3+0];
    T y = rij[ij*3+1];
    T z = rij[ij*3+2];
    T rsq = x * x + y * y + z * z;
    T r = sqrt(rsq);

    //int ii = ai[ij];
    //int jj = alist[aj[ij]];
    T rcutij = (radelem[ti[ij]]+radelem[tj[ij]])*rcutfac; //(radelem[type[ii]]+radelem[type[jj]])*rcutfac;
    T theta0 = (r - rmin0) * rfac0 * M_PI / (rcutij - rmin0);
    //    theta0 = (r - rmin0) * rscale0;
    T z0 = r / tan(theta0);    
            
    T r0inv;
    T a_r, b_r, a_i, b_i;
    T rootpq;
    int jdim = twojmax + 1;
  
    r0inv = 1.0 / sqrt(r * r + z0 * z0);
    a_r = r0inv * z0;
    a_i = -r0inv * z;
    b_r = r0inv * y;
    b_i = -r0inv * x;

    //int nij = idxu_max*ij;    
    ulist_r[ij+0*ijnum] = 1.0;
    ulist_i[ij+0*ijnum] = 0.0;

    for (int j = 1; j <= twojmax; j++) {
        int jju = idxu_block[j];
        int jjup = idxu_block[j-1];
        
        // fill in left side of matrix layer from previous layer
        for (int mb = 0; 2*mb <= j; mb++) {
            ulist_r[ij+jju*ijnum] = 0.0;
            ulist_i[ij+jju*ijnum] = 0.0;
            for (int ma = 0; ma < j; ma++) {
                rootpq = rootpqarray[(j - ma)*jdim + (j - mb)];
                int njju = ij+jju*ijnum;
                int njjup = ij+jjup*ijnum;
                T u_r = ulist_r[njjup];
                T u_i = ulist_i[njjup];
                ulist_r[njju] += rootpq * (a_r * u_r + a_i * u_i);
                ulist_i[njju] += rootpq * (a_r * u_i - a_i * u_r);

                rootpq = rootpqarray[(ma + 1)*jdim + (j - mb)];
                ulist_r[njju+1] = -rootpq * (b_r * u_r + b_i * u_i);
                ulist_i[njju+1] = -rootpq * (b_r * u_i - b_i * u_r);
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
                    ulist_r[njjup] = ulist_r[njju];
                    ulist_i[njjup] = -ulist_i[njju];
                } else {
                    ulist_r[njjup] = -ulist_r[njju];
                    ulist_i[njjup] = ulist_i[njju];
                }
                mapar = -mapar;
                jju++;
                jjup--;
            }
            mbpar = -mbpar;
        }
    }        
    ij += blockDim.x * gridDim.x;
  }
};
template <typename T> void gpuComputeUij2(T *ulist_r, T *ulist_i, T *rootpqarray, T *rij, 
        T *radelem, T rmin0, T rfac0, T rcutfac, int *idxu_block, 
        int *ai, int *aj, int *ti, int *tj, int twojmax, int idxu_max, int ijnum)
{
    int blockDim = 256;
    int gridDim = (ijnum + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;
    gpuKernelComputeUij2<<<gridDim, blockDim>>>(ulist_r, ulist_i, rootpqarray, rij, radelem, 
            rmin0, rfac0, rcutfac, idxu_block, ai, aj, ti, tj, twojmax, idxu_max, ijnum);
}
template void gpuComputeUij2(double*, double*, double*, double*, double*, double, 
        double, double, int*, int*, int*, int*, int*, int, int, int);
template void gpuComputeUij2(float*, float*, float*, float*, float*, float,  float,
        float, int*, int*, int*, int*, int*, int, int, int);

template <typename T> __global__ void gpuKernelComputeZi(T *zlist_r, T *zlist_i, T *ulisttot_r, 
        T *ulisttot_i, T *cglist, int *idxz, int *idxu_block, int *idxcg_block, int jdim, 
        int idxu_max, int idxz_max, int nelements, int bnorm_flag, int N1, int N2, int N3)
{  
  int idx = threadIdx.x + blockIdx.x * blockDim.x;
  while (idx < N3) {
      int l = idx%N2;        
      int elem2 = l%N1;
      int elem1 = (l-elem2)/N1;
      int ii = (idx-l)/N2;    

      int iT = elem2 + nelements*elem1;  
      T * zptr_r = &zlist_r[iT*idxz_max + idxz_max*nelements*nelements*ii];
      T * zptr_i = &zlist_i[iT*idxz_max + idxz_max*nelements*nelements*ii];
      
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
      idx += blockDim.x * gridDim.x;
    }
};
template <typename T> void gpuComputeZi(T *zlist_r, T *zlist_i, T *ulisttot_r, T *ulisttot_i, 
        T *cglist, int *idxz, int *idxu_block, int *idxcg_block, int twojmax, int idxu_max, 
        int idxz_max, int nelements, int bnorm_flag, int inum)
{
    int jdim = twojmax + 1;    
    int N1 = nelements;
    int N2 = N1*nelements;
    int N3 = N2*inum;                        
    int blockDim = 256;
    int gridDim = (N3 + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;
    gpuKernelComputeZi<<<gridDim, blockDim>>>(zlist_r, zlist_i, ulisttot_r, ulisttot_i, cglist, 
       idxz, idxu_block, idxcg_block, jdim, idxu_max, idxz_max, nelements, bnorm_flag, N1, N2, N3);
}
template void gpuComputeZi(double*, double*, double*, double*, double*, 
        int*, int*, int*, int, int, int, int, int, int);
template void gpuComputeZi(float*, float*, float*, float*, float*, 
        int*, int*, int*, int, int, int, int, int, int);

template <typename T> __global__ void gpuKernelComputeYi(T *ylist_r, T *ylist_i, T *ulisttot_r, 
        T *ulisttot_i, T *cglist, T* betalist, int *idxz, int *idxb_block, int *idxu_block, 
        int *idxcg_block, int jdim, int idxb_max, int idxu_max, int idxz_max, int nelements, 
        int bnorm_flag, int ncoeff, int inum)
{
    
  int ii = threadIdx.x + blockIdx.x * blockDim.x;
  while (ii < inum) {  
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
        int jju = idxz[jjz*10+9];
        for(int elem3 = 0; elem3 < nelements; elem3++) {
        // pick out right beta value
          int itriple;
          T betaj;
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
    ii += blockDim.x * gridDim.x;
  }
}

template <typename T> __global__ void gpuKernelComputeYi(T *ylist_r, T *ylist_i, T *ulisttot_r, 
        T *ulisttot_i, T *cglist, T* betalist, int *idxz, int *idxb_block, int *idxu_block, 
        int *idxcg_block, int jdim, int idxb_max, int idxu_max, int idxz_max, int nelements, 
        int bnorm_flag, int ncoeff, int N1, int N2)
{
  int idx = threadIdx.x + blockIdx.x * blockDim.x;
  while (idx < N2) {
      int jjz = idx%N1;              
      int ii = (idx-jjz)/N1;              
      T *beta = &betalist[ncoeff*ii];    
      for(int elem1 = 0; elem1 < nelements; elem1++)
        for (int elem2 = 0; elem2 < nelements; elem2++) {        
            const int j1 = idxz[jjz*10+0];
            const int j2 = idxz[jjz*10+1];
            const int j = idxz[jjz*10+2];
            const int ma1min = idxz[jjz*10+3];
            const int ma2max = idxz[jjz*10+4];
            const int na = idxz[jjz*10+5];
            const int mb1min = idxz[jjz*10+6];
            const int mb2max = idxz[jjz*10+7];
            const int nb = idxz[jjz*10+8];
          
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


            if (bnorm_flag) {
              ztmp_i /= j+1;
              ztmp_r /= j+1;
            }
       
            int jju = idxz[jjz*10+9];
            for(int elem3 = 0; elem3 < nelements; elem3++) {
              int itriple;  
              T betaj;
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

              //ylist_r[elem3 * idxu_max + jju + idxu_max*nelements*ii] += betaj * ztmp_r;
              //ylist_i[elem3 * idxu_max + jju + idxu_max*nelements*ii] += betaj * ztmp_i;
           
              atomicAdd(&ylist_r[elem3 * idxu_max + jju + idxu_max*nelements*ii], betaj * ztmp_r);
              atomicAdd(&ylist_i[elem3 * idxu_max + jju + idxu_max*nelements*ii], betaj * ztmp_i);
           }
        } 
        idx += blockDim.x * gridDim.x;
    }  
}
template <typename T> void gpuComputeYi(T *ylist_r, T *ylist_i, T *ulisttot_r, T *ulisttot_i, T *cglist, 
        T* betalist, int *idxz, int *idxb_block, int *idxu_block, int *idxcg_block, int twojmax, 
        int idxb_max, int idxu_max, int idxz_max, int nelements, int bnorm_flag, int ncoeff, int inum)
{    
    int N1 = idxu_max*nelements*inum;
    gpuArraySetValue(ylist_r, (T) 0.0, N1);
    gpuArraySetValue(ylist_i, (T) 0.0, N1);

     int jdim = twojmax + 1;    
     N1 = idxz_max;
     int N2 = N1*inum;                        
     int blockDim = 256;
     int gridDim = (N2 + blockDim - 1) / blockDim;
     gridDim = (gridDim>1024)? 1024 : gridDim;
     gpuKernelComputeYi<<<gridDim, blockDim>>>(ylist_r, ylist_i, ulisttot_r, ulisttot_i, cglist, 
             betalist, idxz, idxb_block, idxu_block, idxcg_block, jdim, idxb_max, idxu_max, 
                     idxz_max, nelements, bnorm_flag, ncoeff, N1, N2);

// template <typename T> __global__ void gpuKernelComputeYi(T *ylist_r, T *ylist_i, T *ulisttot_r, 
//         T *ulisttot_i, T *cglist, T* betalist, int *idxz, int *idxb_block, int *idxu_block, 
//         int *idxcg_block, int jdim, int idxb_max, int idxu_max, int idxz_max, int nelements, 
//         int bnorm_flag, int ncoeff, int inum)

// template <typename T> __global__ void gpuKernelComputeYi(T *ylist_r, T *ylist_i, T *ulisttot_r, 
//         T *ulisttot_i, T *cglist, T* betalist, int *idxz, int *idxb_block, int *idxu_block, 
//         int *idxcg_block, int jdim, int idxb_max, int idxu_max, int idxz_max, int nelements, 
//         int bnorm_flag, int ncoeff, int N1, int N2)

  //  int jdim = twojmax + 1;                        
  //  int blockDim = 256;
  //  int gridDim = (inum + blockDim - 1) / blockDim;
  //  gridDim = (gridDim>1024)? 1024 : gridDim;
  //  gpuKernelComputeYi<<<gridDim, blockDim>>>(ylist_r, ylist_i, ulisttot_r, ulisttot_i, cglist, 
  //          betalist, idxz, idxb_block, idxu_block, idxcg_block, jdim, idxb_max, idxu_max, 
  //                  idxz_max, nelements, bnorm_flag, ncoeff, inum);
}
template void gpuComputeYi(double*, double*, double*, double*, double*, double*, 
        int*, int*, int*, int*, int, int, int, int, int, int, int, int);
template void gpuComputeYi(float*, float*, float*, float*, float*, float*,
        int*, int*, int*, int*, int, int, int, int, int, int, int, int);

template <typename T> __global__ void gpuKernelComputeBi1(T *blist, T *zlist_r, T *zlist_i, 
        T *ulisttot_r, T *ulisttot_i, int *ilist, int *type, int *map, int *idxb, 
        int *idxu_block, int *idxz_block, int jdim, int idxb_max, int idxu_max, int idxz_max, 
        int nelements, int nelemsq, int nz_max, int nu_max, int nb_max, int N1, int N2, int N3)
{
    // nz_max = idxz_max*nelemsq
    // nu_max = idxu_max*nelements
    // nb_max = idxb_max*nelements*nelemsq        
    int idx = threadIdx.x + blockIdx.x * blockDim.x;
    while (idx < N3) {
        int l = idx%N2;        
        int jjb = l%N1;
        int jelem = (l-jjb)/N1;
        int ii = (idx-l)/N2;    
        
        int k = jelem%nelemsq;      
        int elem3 = k%nelements;
        int elem2 = (k-elem3)/nelements;
        int elem1 = (jelem-k)/nelemsq;    
        int itriple = elem3 + nelements*elem2 + nelemsq*elem1;    

        int iT = elem2 + nelements*elem1;  
        T *zptr_r = &zlist_r[iT*idxz_max + nz_max*ii];
        T *zptr_i = &zlist_i[iT*idxz_max + nz_max*ii];
        
        const int j1 = idxb[jjb*3 + 0];
        const int j2 = idxb[jjb*3 + 1];
        const int j = idxb[jjb*3 + 2];  

        int jjz = idxz_block[j + j2*jdim + j1*jdim*jdim];
        int jju = idxu_block[j];
        T sumzu = 0.0;
        for (int mb = 0; 2 * mb < j; mb++)
            for (int ma = 0; ma <= j; ma++) {
                sumzu += ulisttot_r[elem3*idxu_max+jju + nu_max*ii] * zptr_r[jjz] +
                   ulisttot_i[elem3*idxu_max+jju + nu_max*ii] * zptr_i[jjz];
                jjz++;
                jju++;
            } // end loop over ma, mb

        // For j even, handle middle column
        if (j % 2 == 0) {
            int mb = j / 2;
            for (int ma = 0; ma < mb; ma++) {
                sumzu += ulisttot_r[elem3*idxu_max+jju + nu_max*ii] * zptr_r[jjz] +
                   ulisttot_i[elem3*idxu_max+jju + nu_max*ii] * zptr_i[jjz];
                jjz++;
                jju++;
            }
            sumzu += 0.5 * (ulisttot_r[elem3*idxu_max+jju + nu_max*ii] * zptr_r[jjz] +
                        ulisttot_i[elem3*idxu_max+jju + nu_max*ii] * zptr_i[jjz]);
        } // end if jeven

        blist[itriple*idxb_max + jjb + nb_max*ii] = 2.0 * sumzu;
        
        idx += blockDim.x * gridDim.x;
    }
}

template <typename T> __global__ void gpuKernelComputeBi2(T *blist, T *bzero,int *ilist, int *type,
       int *map, int *idxb, int idxb_max, int nelements, int nb_max, int N1, int N2, int chemflag)
{    
    // nb_max = idxb_max*nelements*nelements*nelements
    int idx = threadIdx.x + blockIdx.x * blockDim.x;
    while (idx < N2) {        
        int jjb = idx%N1;        
        int ii = (idx-jjb)/N1;    
        
        int ielem = (chemflag) ? map[type[ilist[ii]]]: 0;                
        int itriple = (ielem*nelements+ielem)*nelements+ielem;

        const int j = idxb[jjb*3 + 2];  
        blist[itriple*idxb_max + jjb + nb_max*ii] -= bzero[j];
        
        idx += blockDim.x * gridDim.x;
    }
}

template <typename T> __global__ void gpuKernelComputeBi3(T *blist, T *bzero,
       int *idxb, int N1, int N2, int N3)
{    
    // nb_max = idxb_max*nelements*nelements*nelements
    int idx = threadIdx.x + blockIdx.x * blockDim.x;
    while (idx < N3) {
        int l = idx%N2;        
        int jjb = l%N1;
//         int jelem = (l-jjb)/N1;
//         int ii = (idx-l)/N2;            
//         int k = jelem%nelemsq;      
//         int elem3 = k%nelements;
//         int elem2 = (k-elem3)/nelements;
//         int elem1 = (jelem-k)/nelemsq;    
//         int itriple = elem3 + nelements*elem2 + nelemsq*elem1;    
        
        const int j = idxb[jjb*3 + 2];  
        // idx = itriple*idxb_max + jjb + nb_max*ii
        blist[idx] -= bzero[j];
        idx += blockDim.x * gridDim.x;
    }
}

template <typename T> void gpuComputeBi(T *blist, T *zlist_r, T *zlist_i, T *ulisttot_r, T *ulisttot_i, 
        T *bzero,int *ilist, int *type, int *map, int *idxb, int *idxu_block, int *idxz_block, 
        int twojmax, int idxb_max, int idxu_max, int idxz_max, int nelements, int bzero_flag, 
        int wselfall_flag, int chemflag, int inum)
{                
    int nelemsq = nelements*nelements;
    int nz_max = idxz_max*nelemsq;
    int nu_max = idxu_max*nelements;
    int nb_max = idxb_max*nelements*nelemsq;
    int N1 = idxb_max;
    int N2 = N1*nelements*nelemsq;
    int N3 = N2*inum;
    int jdim = twojmax+1;

    int blockDim = 256;
    int gridDim = (N3 + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;    
    gpuKernelComputeBi1<<<gridDim, blockDim>>>(blist, zlist_r, zlist_i, ulisttot_r, ulisttot_i, 
            ilist, type, map, idxb, idxu_block, idxz_block, jdim, idxb_max, idxu_max, idxz_max,
             nelements, nelemsq, nz_max, nu_max, nb_max, N1, N2, N3);

    if (bzero_flag) {
        if (!wselfall_flag) {
            N2 = N1*inum;
            gridDim = (N2 + blockDim - 1) / blockDim;
            gridDim = (gridDim>1024)? 1024 : gridDim;    
            gpuKernelComputeBi2<<<gridDim, blockDim>>>(blist, bzero, ilist, type, map, 
                    idxb, idxb_max, nelements, nb_max, N1, N2, chemflag);
        }
        else {
            gpuKernelComputeBi3<<<gridDim, blockDim>>>(blist, bzero, idxb, N1, N2, N3);            
        }
    }
}
template void gpuComputeBi(double*, double*, double*, double*, double*, double*, 
        int*, int*, int*, int*, int*, int*, int, int, int, int, int, int, int, int, int);
template void gpuComputeBi(float*, float*, float*, float*, float*, float*,
        int*, int*, int*, int*, int*, int*, int, int, int, int, int, int, int, int, int);

template <typename T> __device__ void gpuComputeDuarray(T *dulist_r, T *dulist_i, T *ulist_r, T *ulist_i, 
                   T *rootpqarray, T x, T y, T z, T z0, T r, T dz0dr, T wj, T rcut, T rmin0, 
                   int *idxu_block, int twojmax, int switch_flag)
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
  T sfac = gpuComputeSfac(r, rcut, rmin0, switch_flag);  
  T dsfac = gpuComputeDsfac(r, rcut, rmin0, switch_flag);  
    
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

template <typename T> __global__ void gpuKernelComputeDuijdrj(T *dulist_r, T *dulist_i, T *ulist_r, 
     T *ulist_i, T *rootpqarray, T* rij, T *wjelem, T *radelem, T rmin0, T rfac0, T rcutfac, 
     int *idxu_block, int *ai, int *aj, int *ti, int *tj, int twojmax, int idxu_max, int ijnum, int switch_flag)
{
  int j = threadIdx.x + blockIdx.x * blockDim.x;
  while (j < ijnum) {              
    T x = rij[j*3+0];
    T y = rij[j*3+1];
    T z = rij[j*3+2];
    T rsq = x * x + y * y + z * z;
    T r = sqrt(rsq);

    //int ii = ai[j];
    //int jj = alist[aj[j]];
    T rcutij = (radelem[ti[j]]+radelem[tj[j]])*rcutfac; //(radelem[type[ii]]+radelem[type[jj]])*rcutfac;
    T rscale0 = rfac0 * M_PI / (rcutij - rmin0);
    T theta0 = (r - rmin0) * rscale0;
    T z0 = r / tan(theta0);    
    //T cs = cos(theta0);
    //T sn = sin(theta0);
    //T z0 = r * cs / sn;
    T dz0dr = z0 / r - (r*rscale0) * (rsq + z0 * z0) / rsq;
    
    gpuComputeDuarray(&dulist_r[3*idxu_max*j], &dulist_i[3*idxu_max*j], &ulist_r[idxu_max*j], 
            &ulist_i[idxu_max*j], rootpqarray, x, y, z, z0, r, dz0dr, wjelem[tj[j]], rcutij, rmin0,
            idxu_block, twojmax, switch_flag);    

    j += blockDim.x * gridDim.x;
  }  
}
template <typename T> void gpuComputeDuijdrj(T *dulist_r, T *dulist_i, T *ulist_r, T *ulist_i, 
    T *rootpqarray, T* rij, T *wjelem, T *radelem, T rmin0, T rfac0, T rcutfac, int *idxu_block,
    int *ai, int *aj, int *ti, int *tj, int twojmax, int idxu_max, int ijnum, int switch_flag)
{                
    int blockDim = 256;
    int gridDim = (ijnum + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;    
    gpuKernelComputeDuijdrj<<<gridDim, blockDim>>>(dulist_r, dulist_i, ulist_r, ulist_i, 
          rootpqarray, rij, wjelem, radelem, rmin0, rfac0, rcutfac, idxu_block, 
          ai, aj, ti, tj, twojmax, idxu_max, ijnum, switch_flag);
}
template void gpuComputeDuijdrj(double*, double*, double*, double*, double*, double*, 
        double*, double*, double, double, double, int*, int*, int*, int*, int*, int, int, int, int);
template void gpuComputeDuijdrj(float*, float*, float*, float*, float*, float*, float*, 
        float*, float, float, float, int*, int*, int*, int*, int*, int, int, int, int);

template <typename T> __global__ void gpuKernelComputeDeidrj(T *dedr, T *ylist_r, T *ylist_i, 
        T *dulist_r, T *dulist_i, int *idxu_block, int *map, int *ai, int *aj, int *ti, int *tj,
        int nelements, int twojmax, int idxu_max, int chemflag, int ijnum) 
{

//   for(int ij = 0; ij < ijnum; ij++)  
//     for(int k = 0; k < 3; k++)
//         dedr[k + 3*ij] = 0.0;

  //for(int ij = 0; ij < ijnum; ij++) {
  int ij = threadIdx.x + blockIdx.x * blockDim.x;
  while (ij < ijnum) {                        
      int jelem = (chemflag) ? map[tj[ij]] : 0; //(chemflag) ? map[type[alist[aj[ij]]]] : 0;
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

      ij += blockDim.x * gridDim.x;
  }  
}
template <typename T> void gpuComputeDeidrj(T *dedr, T *ylist_r, T *ylist_i, 
        T *dulist_r, T *dulist_i, int *idxu_block, int *map, int *ai, int *aj, int *ti, int *tj,
        int nelements, int twojmax, int idxu_max, int chemflag, int ijnum) 
{                
    gpuArraySetValue(dedr, (T) 0.0, ijnum*3);

    int blockDim = 256;
    int gridDim = (ijnum + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;    
    gpuKernelComputeDeidrj<<<gridDim, blockDim>>>(dedr, ylist_r, ylist_i, dulist_r, dulist_i, 
            idxu_block, map, ai, aj, ti, tj, nelements, twojmax, idxu_max, chemflag, ijnum); 
}
template void gpuComputeDeidrj(double*, double*, double*, double*, double*, 
        int*, int*, int*, int*, int*, int*, int, int, int, int, int);
template void gpuComputeDeidrj(float*, float*, float*, float*, float*, 
        int*, int*, int*, int*, int*, int*, int, int, int, int, int);

template <typename T> __global__ void gpuKernelComputeDbidrj(T *dblist, T *zlist_r, T *zlist_i, 
        T *dulist_r, T *dulist_i, int *idxb, int *idxu_block, int *idxz_block, 
        int *map, int *ai, int *aj, int *ti, int *tj, int jdim, int idxb_max, int idxu_max, int idxz_max, 
        int nelements, int bnorm_flag, int chemflag, int nb_max, int nz_max, int N1, int N2)
{  
  // set all the derivatives to zero once
//   for(int ij = 0; ij < ijnum; ij++)
//   for(int jjb = 0; jjb < idxb_max; jjb++) {
//     for(int elem1 = 0; elem1 < nelements; elem1++)
//       for(int elem2 = 0; elem2 < nelements; elem2++)
//         for(int elem3 = 0; elem3 < nelements; elem3++) {
// 
//           itriple = (elem1 * nelements + elem2) * nelements + elem3;
// 
//           dbdr = &dblist[(itriple*idxb_max+jjb)*3 + 3*idxb_max*nelements*nelements*nelements*ij];
//           dbdr[0] = 0.0;
//           dbdr[1] = 0.0;
//           dbdr[2] = 0.0;
//         }
//   }

    int idx = threadIdx.x + blockIdx.x * blockDim.x;
    while (idx < N2) {
        int jjb = idx%N1;              
        int ij = (idx-jjb)/N1;                      
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
            int iT = elem1*nelements+elem2;
            int itriple = (elem1*nelements+elem2)*nelements+elem3;
            T *dbdr = &dblist[(itriple*idxb_max+jjb)*3 + 3*nb_max*ij];
            T *zptr_r = &zlist_r[iT*idxz_max + nz_max*i];
            T *zptr_i = &zlist_i[iT*idxz_max + nz_max*i];

            T *dudr_r, *dudr_i;
            T sumzdu_r[3];
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
            dbdr = &dblist[(itriple*idxb_max+jjb)*3 + 3*nb_max*ij];
            //jjz = idxz_block[j][j2][j1];
            jjz = idxz_block[j1 + j2*jdim + j*jdim*jdim];
            jju = idxu_block[j1];
            zptr_r = &zlist_r[iT*idxz_max  + nz_max*i];
            zptr_i = &zlist_i[iT*idxz_max  + nz_max*i];

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
            dbdr = &dblist[(itriple*idxb_max+jjb)*3 + 3*nb_max*ij];
            //jjz = idxz_block[j][j1][j2];
            jjz = idxz_block[j2 + j1*jdim + j*jdim*jdim];        
            jju = idxu_block[j2];
            zptr_r = &zlist_r[iT*idxz_max  + nz_max*i];
            zptr_i = &zlist_i[iT*idxz_max  + nz_max*i];

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
        idx += blockDim.x * gridDim.x;
    }
}
template <typename T> void gpuComputeDbidrj(T *dblist, T *zlist_r, T *zlist_i, 
        T *dulist_r, T *dulist_i, int *idxb, int *idxu_block, int *idxz_block, 
        int *map, int *ai, int *aj, int *ti, int *tj, int twojmax, int idxb_max, int idxu_max, int idxz_max, 
        int nelements, int bnorm_flag, int chemflag, int ijnum)
{                
    
    int nelemsq = nelements*nelements;
    int nz_max = idxz_max*nelemsq;
    int nb_max = idxb_max*nelements*nelemsq;
    int N1 = idxb_max;
    int N2 = N1*ijnum;
    int jdim = twojmax+1;

    gpuArraySetValue(dblist, (T) 0.0, nb_max*ijnum*3);
    
    int blockDim = 256;
    int gridDim = (N2 + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;        
    gpuKernelComputeDbidrj<<<gridDim, blockDim>>>(dblist, zlist_r, zlist_i, dulist_r, dulist_i, 
            idxb, idxu_block, idxz_block, map, ai, aj, ti, tj, jdim, idxb_max, idxu_max, 
            idxz_max, nelements, bnorm_flag, chemflag, nb_max, nz_max, N1, N2);
}
template void gpuComputeDbidrj(double*, double*, double*, double*, double*, 
        int*, int*, int*, int*, int*, int*, int*, int*, int, int, int, int, int, int, int, int);
template void gpuComputeDbidrj(float*, float*, float*, float*, float*, 
        int*, int*, int*, int*, int*, int*, int*, int*, int, int, int, int, int, int, int, int);

template <typename T> __global__ void gpuKernelComputeSna1(T *sna, T *blist, int *ilist, 
        int *mask, int ncoeff, int nrows, int N1, int N2)
{    
    //N1 = ncoeff, N2 = ncoeff*inum            
    int idx = threadIdx.x + blockIdx.x * blockDim.x;
    while (idx < N2) {
        int icoeff = idx%N1;              
        int ii = (idx-icoeff)/N1;         
        int i = ilist[ii];
        sna[icoeff + nrows*i] = (mask[i]) ? blist[icoeff + ncoeff*ii] : 0.0;      
        idx += blockDim.x * gridDim.x;
    }
}

template <typename T> __global__ void gpuKernelComputeSna1(T *sna, T *blist, int *ilist, 
        int *mask, int *type, int ncoeff, int ntype, int nperdim, int N1, int N2)
{    
    //N1 = ncoeff, N2 = ncoeff*inum            
    int idx = threadIdx.x + blockIdx.x * blockDim.x;
    while (idx < N2) {
        int icoeff = idx%N1;              
        int ii = (idx-icoeff)/N1;         
        int i = ilist[ii];
        int itype = type[i]-1;        
        sna[icoeff + nperdim*itype + nperdim*ntype*i] = (mask[i]) ? blist[icoeff + ncoeff*ii] : 0.0;      
        idx += blockDim.x * gridDim.x;
    }
}

template <typename T> __global__ void gpuKernelComputeSna2(T *sna, T *blist, int *ilist, 
        int *mask, int ncoeff, int nrows, int N1, int N2)
{    
    //N1 = ncoeff, N2 = ncoeff*inum            
    int idx = threadIdx.x + blockIdx.x * blockDim.x;
    while (idx < N2) {
        int icoeff = idx%N1;              
        int ii = (idx-icoeff)/N1;         
        int i = ilist[ii];
        int imsk = mask[i];
        T bi = blist[icoeff + ncoeff*ii];

        // diagonal element of quadratic matrix
        int ncount = ncoeff + icoeff*ncoeff - icoeff*(icoeff-1)/2;
        sna[ncount + nrows*i] = (imsk) ? 0.5*bi*bi : 0.0;
        ncount++;

        // upper-triangular elements of quadratic matrix
        for (int jcoeff = icoeff+1; jcoeff < ncoeff; jcoeff++) {
            sna[ncount + nrows*i] = (imsk) ? bi*blist[jcoeff + ncoeff*ii] : 0.0;                        
            ncount++;
        }
        
        idx += blockDim.x * gridDim.x;
    }
}

template <typename T> __global__ void gpuKernelComputeSna2(T *sna, T *blist, int *ilist, 
        int *mask, int *type, int ncoeff, int ntype, int nperdim, int N1, int N2)
{    
    //N1 = ncoeff, N2 = ncoeff*inum            
    int idx = threadIdx.x + blockIdx.x * blockDim.x;
    while (idx < N2) {
        int icoeff = idx%N1;              
        int ii = (idx-icoeff)/N1;         
        int i = ilist[ii];
        int itype = type[i]-1;        
        int imsk = mask[i];
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

        idx += blockDim.x * gridDim.x;
    }
}
template <typename T> void gpuComputeSna(T *sna, T *blist, int *ilist, int *mask, 
        int ncoeff, int nrows, int inum, int quadraticflag)
{
    int N1 = ncoeff;
    int N2 = ncoeff*inum;
    int blockDim = 256;
    int gridDim = (N2 + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;        
    gpuKernelComputeSna1<<<gridDim, blockDim>>>(sna, blist, ilist, mask, ncoeff, nrows, N1, N2);

    if (quadraticflag)
        gpuKernelComputeSna2<<<gridDim, blockDim>>>(sna, blist, ilist, mask, ncoeff, nrows, N1, N2);            
}
template void gpuComputeSna(double*, double*, int*, int*, int, int, int, int);
template void gpuComputeSna(float*, float*, int*, int*, int, int, int, int);

template <typename T> void gpuComputeSna(T *sna, T *blist, int *ilist, int *mask, int *type,
        int ncoeff, int ntype, int nperdim, int inum, int quadraticflag)
{
    int N1 = ncoeff;
    int N2 = ncoeff*inum;
    int blockDim = 256;
    int gridDim = (N2 + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;        
    gpuKernelComputeSna1<<<gridDim, blockDim>>>(sna, blist, ilist, mask, type, ncoeff, 
            ntype, nperdim, N1, N2);

    if (quadraticflag)
        gpuKernelComputeSna2<<<gridDim, blockDim>>>(sna, blist, ilist, mask, type, ncoeff, 
            ntype, nperdim, N1, N2);
}
template void gpuComputeSna(double*, double*, int*, int*, int*, int, int, int, int, int);
template void gpuComputeSna(float*, float*, int*, int*, int*, int, int, int, int, int);

template <typename T> __global__ void gpuKernelComputeSnad1(T *snad, T *dblist, T *blist, 
        int *aii, int *ai, int *aj, int *ti, int *mask, int ncoeff, int ntypes, 
        int nperdim, int N1, int N2)
{    
    //N1 = ncoeff, N2 = ncoeff*ijnum            
    int idx = threadIdx.x + blockIdx.x * blockDim.x;
    while (idx < N2) {
        int icoeff = idx%N1;              
        int k = (idx-icoeff)/N1;         
        int i = ai[k];
        if (mask[i]) {
            int j = aj[k];
            int itype = ti[k];
            int nb = 3*icoeff + 3*ncoeff*k;
            int ni = 3*icoeff + 3*nperdim*(itype-1) + 3*nperdim*ntypes*i;
            int nj = 3*icoeff + 3*nperdim*(itype-1) + 3*nperdim*ntypes*j;

            atomicAdd(&snad[0+ni], dblist[0+nb]);
            atomicAdd(&snad[1+ni], dblist[1+nb]);
            atomicAdd(&snad[2+ni], dblist[2+nb]);
            atomicAdd(&snad[0+nj],-dblist[0+nb]);
            atomicAdd(&snad[1+nj],-dblist[1+nb]);
            atomicAdd(&snad[2+nj],-dblist[2+nb]);
        }
        idx += blockDim.x * gridDim.x;
    }
}

template <typename T> __global__ void gpuKernelComputeSnad1(T *snad, T *dblist, T *blist, 
        int *aii, int *ai, int *aj, int *ti, int *mask, int *tag, int ncoeff, int ntypes, 
        int nperdim, int N1, int N2)
{    
    //N1 = ncoeff, N2 = ncoeff*ijnum            
    int idx = threadIdx.x + blockIdx.x * blockDim.x;
    while (idx < N2) {
        int icoeff = idx%N1;              
        int k = (idx-icoeff)/N1;         
        int i = ai[k];
        if (mask[i]) {            
            int j = aj[k];
            int ig = tag[i]-1; // global index of atom i                   
            int jg = tag[j]-1; // global index of atom j
            int itype = ti[k];
            int nb = 3*icoeff + 3*ncoeff*k;
            int ni = nperdim*(itype-1) + nperdim*ntypes*3*ig;
            int nj = nperdim*(itype-1) + nperdim*ntypes*3*jg;

            atomicAdd(&snad[icoeff+ni], dblist[0+nb]);
            atomicAdd(&snad[icoeff+nperdim*ntypes+ni], dblist[1+nb]);
            atomicAdd(&snad[icoeff+2*nperdim*ntypes+ni], dblist[2+nb]);
            atomicAdd(&snad[icoeff+nj],-dblist[0+nb]);
            atomicAdd(&snad[icoeff+nperdim*ntypes+nj],-dblist[1+nb]);
            atomicAdd(&snad[icoeff+2*nperdim*ntypes+nj],-dblist[2+nb]);
        }
        idx += blockDim.x * gridDim.x;
    }
}

template <typename T> __global__ void gpuKernelComputeSnad2(T *snad, T *dblist, T *blist, 
        int *aii, int *ai, int *aj, int *ti, int *mask, int ncoeff, int ntypes, 
        int nperdim, int N1, int N2)
{    
    //N1 = ncoeff, N2 = ncoeff*ijnum            
    int idx = threadIdx.x + blockIdx.x * blockDim.x;
    while (idx < N2) {
        int icoeff = idx%N1;              
        int k = (idx-icoeff)/N1;         
        int i = ai[k];
        if (mask[i]) {
            int j = aj[k];
            int itype = ti[k];
            int nb = ncoeff*aii[k];
            int ndb = 3*ncoeff*k;            

            T bi = blist[icoeff+nb];
            T bix = dblist[0+3*icoeff+ndb];
            T biy = dblist[1+3*icoeff+ndb];
            T biz = dblist[2+3*icoeff+ndb];

            // diagonal elements of quadratic matrix
            T dbxtmp = bi*bix;
            T dbytmp = bi*biy;
            T dbztmp = bi*biz;
            
            int ni = 3*nperdim*ntypes*i + 3*nperdim*(itype-1);
            int nj = 3*nperdim*ntypes*j + 3*nperdim*(itype-1);
            int ncount = ncoeff + icoeff*ncoeff - (icoeff-1)*icoeff/2;            
            atomicAdd(&snad[0+3*ncount+ni], dbxtmp);
            atomicAdd(&snad[1+3*ncount+ni], dbytmp);
            atomicAdd(&snad[2+3*ncount+ni], dbztmp);
            atomicAdd(&snad[0+3*ncount+nj],-dbxtmp);
            atomicAdd(&snad[1+3*ncount+nj],-dbytmp);
            atomicAdd(&snad[2+3*ncount+nj],-dbztmp);
            ncount++;
            // upper-triangular elements of quadratic matrix

            for (int jcoeff = icoeff+1; jcoeff < ncoeff; jcoeff++) {
                dbxtmp = bi*dblist[0+3*jcoeff+ndb] + bix*blist[jcoeff+nb];
                dbytmp = bi*dblist[1+3*jcoeff+ndb] + biy*blist[jcoeff+nb];
                dbztmp = bi*dblist[2+3*jcoeff+ndb] + biz*blist[jcoeff+nb];
                atomicAdd(&snad[0+3*ncount+ni], dbxtmp);
                atomicAdd(&snad[1+3*ncount+ni], dbytmp);
                atomicAdd(&snad[2+3*ncount+ni], dbztmp);
                atomicAdd(&snad[0+3*ncount+nj],-dbxtmp);
                atomicAdd(&snad[1+3*ncount+nj],-dbytmp);
                atomicAdd(&snad[2+3*ncount+nj],-dbztmp);
                ncount++;              
            }                            
        }
        idx += blockDim.x * gridDim.x;
    }
}

template <typename T> __global__ void gpuKernelComputeSnad2(T *snad, T *dblist, T *blist, 
        int *aii, int *ai, int *aj, int *ti, int *mask, int *tag, int ncoeff, int ntypes, 
        int nperdim, int N1, int N2)
{    
    //N1 = ncoeff, N2 = ncoeff*ijnum            
    int idx = threadIdx.x + blockIdx.x * blockDim.x;
    while (idx < N2) {
        int icoeff = idx%N1;              
        int k = (idx-icoeff)/N1;         
        int i = ai[k];
        if (mask[i]) {            
            int j = aj[k];
            int ig = tag[i]-1; // global index of atom i                   
            int jg = tag[j]-1; // global index of atom j
            int itype = ti[k];
            int nb = ncoeff*aii[k];
            int ndb = 3*ncoeff*k;            

            T bi = blist[icoeff+nb];
            T bix = dblist[0+3*icoeff+ndb];
            T biy = dblist[1+3*icoeff+ndb];
            T biz = dblist[2+3*icoeff+ndb];

            // diagonal elements of quadratic matrix
            T dbxtmp = bi*bix;
            T dbytmp = bi*biy;
            T dbztmp = bi*biz;

            int ni = nperdim*(itype-1) + nperdim*ntypes*3*ig;
            int nj = nperdim*(itype-1) + nperdim*ntypes*3*jg;
            int ncount = ncoeff + icoeff*ncoeff - (icoeff-1)*icoeff/2;            

            atomicAdd(&snad[ncount+ni], dbxtmp);
            atomicAdd(&snad[ncount+nperdim*ntypes+ni], dbytmp);
            atomicAdd(&snad[ncount+2*nperdim*ntypes+ni], dbztmp);
            atomicAdd(&snad[ncount+nj],-dbxtmp);
            atomicAdd(&snad[ncount+nperdim*ntypes+nj],-dbytmp);
            atomicAdd(&snad[ncount+2*nperdim*ntypes+nj],-dbztmp);
            ncount++;

            for (int jcoeff = icoeff+1; jcoeff < ncoeff; jcoeff++) {
                dbxtmp = bi*dblist[0+3*jcoeff+ndb] + bix*blist[jcoeff+nb];
                dbytmp = bi*dblist[1+3*jcoeff+ndb] + biy*blist[jcoeff+nb];
                dbztmp = bi*dblist[2+3*jcoeff+ndb] + biz*blist[jcoeff+nb];
                atomicAdd(&snad[ncount+ni], dbxtmp);
                atomicAdd(&snad[ncount+nperdim*ntypes+ni], dbytmp);
                atomicAdd(&snad[ncount+2*nperdim*ntypes+ni], dbztmp);
                atomicAdd(&snad[ncount+nj],-dbxtmp);
                atomicAdd(&snad[ncount+nperdim*ntypes+nj],-dbytmp);
                atomicAdd(&snad[ncount+2*nperdim*ntypes+nj],-dbztmp);
                ncount++;              
            }                            
        }
        idx += blockDim.x * gridDim.x;
    }
}

template <typename T> void gpuComputeSnad(T *snad, T *dblist, T *blist, int *aii, int *ai, int *aj, int *ti, 
        int *mask, int ncoeff, int ntypes, int nperdim, int ijnum, int quadraticflag)
{        
    int N1 = ncoeff;
    int N2 = ncoeff*ijnum;
    int blockDim = 256;
    int gridDim = (N2 + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;        
    gpuKernelComputeSnad1<<<gridDim, blockDim>>>(snad, dblist, blist, aii, ai, aj, 
          ti, mask, ncoeff, ntypes, nperdim, N1, N2);

    if (quadraticflag)
        gpuKernelComputeSnad2<<<gridDim, blockDim>>>(snad, dblist, blist, aii, ai, aj, 
          ti, mask, ncoeff, ntypes, nperdim, N1, N2);
}
template void gpuComputeSnad(double*, double*, double*, int*, int*, int*, int*, int*,
        int, int, int, int, int);
template void gpuComputeSnad(float*, float*, float*, int*, int*, int*, int*, int*,
        int, int, int, int, int);

template <typename T> void gpuComputeSnad(T *snad, T *dblist, T *blist, int *aii, int *ai, int *aj, int *ti, 
        int *mask, int *tag, int ncoeff, int ntypes, int nperdim, int ijnum, int quadraticflag)
{        
    int N1 = ncoeff;
    int N2 = ncoeff*ijnum;
    int blockDim = 256;
    int gridDim = (N2 + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;        
    gpuKernelComputeSnad1<<<gridDim, blockDim>>>(snad, dblist, blist, aii, ai, aj, 
          ti, mask, tag, ncoeff, ntypes, nperdim, N1, N2);

    if (quadraticflag)
        gpuKernelComputeSnad2<<<gridDim, blockDim>>>(snad, dblist, blist, aii, ai, aj, 
          ti, mask, tag, ncoeff, ntypes, nperdim, N1, N2);
}
template void gpuComputeSnad(double*, double*, double*, int*, int*, int*, int*, int*, int*,
        int, int, int, int, int);
template void gpuComputeSnad(float*, float*, float*, int*, int*, int*, int*, int*, int*,
        int, int, int, int, int);

template <typename T> __global__ void gpuKernelComputeSnav11(T *snav, T *dblist, T *blist, T *x, 
        int *aii, int *ai, int *aj,int *ti, int *mask, int ncoeff, int ntypes, int nperdim, 
        int N1, int N2)
{    
    int idx = threadIdx.x + blockIdx.x * blockDim.x;
    while (idx < N2) {
        int icoeff = idx%N1;              
        int k = (idx-icoeff)/N1;         
        int i = ai[k];
        if (mask[i]) {
            int j = aj[k];
            int itype = ti[k];
            int ndb = 3*icoeff + 3*ncoeff*k;
            int ni = 6*icoeff + 6*nperdim*(itype-1) + 6*nperdim*ntypes*i;
            int nj = 6*icoeff + 6*nperdim*(itype-1) + 6*nperdim*ntypes*j;
            T xi = x[0+3*i]; 
            T yi = x[1+3*i]; 
            T zi = x[2+3*i]; 
            T xj = x[0+3*j]; 
            T yj = x[1+3*j]; 
            T zj = x[2+3*j];             

            atomicAdd(&snav[0+ni], dblist[0+ndb]*xi);
            atomicAdd(&snav[1+ni], dblist[1+ndb]*yi);
            atomicAdd(&snav[2+ni], dblist[2+ndb]*zi);
            atomicAdd(&snav[3+ni], dblist[1+ndb]*zi);
            atomicAdd(&snav[4+ni], dblist[0+ndb]*zi);
            atomicAdd(&snav[5+ni], dblist[0+ndb]*yi);
            atomicAdd(&snav[0+nj],-dblist[0+ndb]*xj);
            atomicAdd(&snav[1+nj],-dblist[1+ndb]*yj);
            atomicAdd(&snav[2+nj],-dblist[2+ndb]*zj);
            atomicAdd(&snav[3+nj],-dblist[1+ndb]*zj);
            atomicAdd(&snav[4+nj],-dblist[0+ndb]*zj);
            atomicAdd(&snav[5+nj],-dblist[0+ndb]*yj);
        }
        idx += blockDim.x * gridDim.x;
    }
}

template <typename T> __global__ void gpuKernelComputeSnav21(T *snav, T *dblist, T *blist, T *x, 
        int *aii, int *ai, int *aj, int *ti, int *mask, int ncoeff, int ntypes, int nperdim, 
        int N1, int N2)
{    
    int idx = threadIdx.x + blockIdx.x * blockDim.x;
    while (idx < N2) {
        int icoeff = idx%N1;              
        int k = (idx-icoeff)/N1;         
        int i = ai[k];
        if (mask[i]) {
            int j = aj[k];
            int itype = ti[k];
            int ndb = 3*icoeff + 3*ncoeff*k;
            int ni = nperdim*(itype-1) + 6*nperdim*ntypes*i;
            int nj = nperdim*(itype-1) + 6*nperdim*ntypes*j;
            T xi = x[0+3*i]; 
            T yi = x[1+3*i]; 
            T zi = x[2+3*i]; 
            T xj = x[0+3*j]; 
            T yj = x[1+3*j]; 
            T zj = x[2+3*j];             

            atomicAdd(&snav[icoeff+ni], dblist[0+ndb]*xi);
            atomicAdd(&snav[icoeff+nperdim*ntypes+ni], dblist[1+ndb]*yi);
            atomicAdd(&snav[icoeff+2*nperdim*ntypes+ni], dblist[2+ndb]*zi);
            atomicAdd(&snav[icoeff+3*nperdim*ntypes+ni], dblist[1+ndb]*zi);
            atomicAdd(&snav[icoeff+4*nperdim*ntypes+ni], dblist[0+ndb]*zi);
            atomicAdd(&snav[icoeff+5*nperdim*ntypes+ni], dblist[0+ndb]*yi);
            atomicAdd(&snav[icoeff+nj],-dblist[0+ndb]*xj);
            atomicAdd(&snav[icoeff+nperdim*ntypes+nj],-dblist[1+ndb]*yj);
            atomicAdd(&snav[icoeff+2*nperdim*ntypes+nj],-dblist[2+ndb]*zj);
            atomicAdd(&snav[icoeff+3*nperdim*ntypes+nj],-dblist[1+ndb]*zj);
            atomicAdd(&snav[icoeff+4*nperdim*ntypes+nj],-dblist[0+ndb]*zj);
            atomicAdd(&snav[icoeff+5*nperdim*ntypes+nj],-dblist[0+ndb]*yj);
        }
        idx += blockDim.x * gridDim.x;
    }
}

template <typename T> __global__ void gpuKernelComputeSnav12(T *snav, T *dblist, T *blist, T *x, 
        int *aii, int *ai, int *aj,int *ti, int *mask, int ncoeff, int ntypes, int nperdim, 
        int N1, int N2)
{    
    int idx = threadIdx.x + blockIdx.x * blockDim.x;
    while (idx < N2) {
        int icoeff = idx%N1;              
        int k = (idx-icoeff)/N1;         
        int i = ai[k];
        if (mask[i]) {
            int j = aj[k];
            int itype = ti[k];
            int nb = ncoeff*aii[k];
            int ndb = 3*icoeff + 3*ncoeff*k;
            T xi = x[0+3*i]; 
            T yi = x[1+3*i]; 
            T zi = x[2+3*i]; 
            T xj = x[0+3*j]; 
            T yj = x[1+3*j]; 
            T zj = x[2+3*j];             

            T bi = blist[icoeff+nb];
            T bix = dblist[0+ndb];
            T biy = dblist[1+ndb];
            T biz = dblist[2+ndb];

            // diagonal elements of quadratic matrix
            T dbxtmp = bi*bix;
            T dbytmp = bi*biy;
            T dbztmp = bi*biz;

            int ncount = ncoeff + icoeff*ncoeff - (icoeff-1)*icoeff/2;         
            int mi = 6*nperdim*(itype-1) + 6*nperdim*ntypes*i;
            int mj = 6*nperdim*(itype-1) + 6*nperdim*ntypes*j;
            int ni = 6*ncount + mi;
            int nj = 6*ncount + mj;

            atomicAdd(&snav[0+ni], dbxtmp*xi);
            atomicAdd(&snav[1+ni], dbytmp*yi);
            atomicAdd(&snav[2+ni], dbztmp*zi);
            atomicAdd(&snav[3+ni], dbytmp*zi);
            atomicAdd(&snav[4+ni], dbxtmp*zi);
            atomicAdd(&snav[5+ni], dbxtmp*yi);
            atomicAdd(&snav[0+nj],-dbxtmp*xj);
            atomicAdd(&snav[1+nj],-dbytmp*yj);
            atomicAdd(&snav[2+nj],-dbztmp*zj);
            atomicAdd(&snav[3+nj],-dbytmp*zj);
            atomicAdd(&snav[4+nj],-dbxtmp*zj);
            atomicAdd(&snav[5+nj],-dbxtmp*yj);
            ncount++;      
            ni = 6*ncount + mi;
            nj = 6*ncount + mj;

            for (int jcoeff = icoeff+1; jcoeff < ncoeff; jcoeff++) {
                dbxtmp = bi*dblist[0+3*jcoeff+ndb] + bix*blist[jcoeff+nb];
                dbytmp = bi*dblist[1+3*jcoeff+ndb] + biy*blist[jcoeff+nb];
                dbztmp = bi*dblist[2+3*jcoeff+ndb] + biz*blist[jcoeff+nb];

                atomicAdd(&snav[0+ni], dbxtmp*xi);
                atomicAdd(&snav[1+ni], dbytmp*yi);
                atomicAdd(&snav[2+ni], dbztmp*zi);
                atomicAdd(&snav[3+ni], dbytmp*zi);
                atomicAdd(&snav[4+ni], dbxtmp*zi);
                atomicAdd(&snav[5+ni], dbxtmp*yi);
                atomicAdd(&snav[0+nj],-dbxtmp*xj);
                atomicAdd(&snav[1+nj],-dbytmp*yj);
                atomicAdd(&snav[2+nj],-dbztmp*zj);
                atomicAdd(&snav[3+nj],-dbytmp*zj);
                atomicAdd(&snav[4+nj],-dbxtmp*zj);
                atomicAdd(&snav[5+nj],-dbxtmp*yj);
                ncount++;      
                ni = 6*ncount + mi;
                nj = 6*ncount + mj;                
            }                            
        }
        idx += blockDim.x * gridDim.x;
    }
}

template <typename T> __global__ void gpuKernelComputeSnav22(T *snav, T *dblist, T *blist, T *x, 
        int *aii, int *ai, int *aj,int *ti, int *mask, int ncoeff, int ntypes, int nperdim, 
        int N1, int N2)
{    
    int idx = threadIdx.x + blockIdx.x * blockDim.x;
    while (idx < N2) {
        int icoeff = idx%N1;              
        int k = (idx-icoeff)/N1;         
        int i = ai[k];
        if (mask[i]) {
            int j = aj[k];
            int itype = ti[k];
            int nb = ncoeff*aii[k];
            int ndb = 3*icoeff + 3*ncoeff*k;
            T xi = x[0+3*i]; 
            T yi = x[1+3*i]; 
            T zi = x[2+3*i]; 
            T xj = x[0+3*j]; 
            T yj = x[1+3*j]; 
            T zj = x[2+3*j];             

            T bi = blist[icoeff+nb];
            T bix = dblist[0+ndb];
            T biy = dblist[1+ndb];
            T biz = dblist[2+ndb];

            // diagonal elements of quadratic matrix
            T dbxtmp = bi*bix;
            T dbytmp = bi*biy;
            T dbztmp = bi*biz;

            int ncount = ncoeff + icoeff*ncoeff - (icoeff-1)*icoeff/2;         
            int ni = nperdim*(itype-1) + 6*nperdim*ntypes*i;
            int nj = nperdim*(itype-1) + 6*nperdim*ntypes*j;
            
            atomicAdd(&snav[ncount+ni], dblist[0+ndb]*xi);
            atomicAdd(&snav[ncount+nperdim*ntypes+ni], dblist[1+ndb]*yi);
            atomicAdd(&snav[ncount+2*nperdim*ntypes+ni], dblist[2+ndb]*zi);
            atomicAdd(&snav[ncount+3*nperdim*ntypes+ni], dblist[1+ndb]*zi);
            atomicAdd(&snav[ncount+4*nperdim*ntypes+ni], dblist[0+ndb]*zi);
            atomicAdd(&snav[ncount+5*nperdim*ntypes+ni], dblist[0+ndb]*yi);
            atomicAdd(&snav[ncount+nj],-dblist[0+ndb]*xj);
            atomicAdd(&snav[ncount+nperdim*ntypes+nj],-dblist[1+ndb]*yj);
            atomicAdd(&snav[ncount+2*nperdim*ntypes+nj],-dblist[2+ndb]*zj);
            atomicAdd(&snav[ncount+3*nperdim*ntypes+nj],-dblist[1+ndb]*zj);
            atomicAdd(&snav[ncount+4*nperdim*ntypes+nj],-dblist[0+ndb]*zj);
            atomicAdd(&snav[ncount+5*nperdim*ntypes+nj],-dblist[0+ndb]*yj);
            ncount++;      

            for (int jcoeff = icoeff+1; jcoeff < ncoeff; jcoeff++) {
                dbxtmp = bi*dblist[0+3*jcoeff+ndb] + bix*blist[jcoeff+nb];
                dbytmp = bi*dblist[1+3*jcoeff+ndb] + biy*blist[jcoeff+nb];
                dbztmp = bi*dblist[2+3*jcoeff+ndb] + biz*blist[jcoeff+nb];

                atomicAdd(&snav[ncount+ni], dblist[0+ndb]*xi);
                atomicAdd(&snav[ncount+nperdim*ntypes+ni], dblist[1+ndb]*yi);
                atomicAdd(&snav[ncount+2*nperdim*ntypes+ni], dblist[2+ndb]*zi);
                atomicAdd(&snav[ncount+3*nperdim*ntypes+ni], dblist[1+ndb]*zi);
                atomicAdd(&snav[ncount+4*nperdim*ntypes+ni], dblist[0+ndb]*zi);
                atomicAdd(&snav[ncount+5*nperdim*ntypes+ni], dblist[0+ndb]*yi);
                atomicAdd(&snav[ncount+nj],-dblist[0+ndb]*xj);
                atomicAdd(&snav[ncount+nperdim*ntypes+nj],-dblist[1+ndb]*yj);
                atomicAdd(&snav[ncount+2*nperdim*ntypes+nj],-dblist[2+ndb]*zj);
                atomicAdd(&snav[ncount+3*nperdim*ntypes+nj],-dblist[1+ndb]*zj);
                atomicAdd(&snav[ncount+4*nperdim*ntypes+nj],-dblist[0+ndb]*zj);
                atomicAdd(&snav[ncount+5*nperdim*ntypes+nj],-dblist[0+ndb]*yj);
                ncount++;      
            }                            
        }
        idx += blockDim.x * gridDim.x;
    }
}

template <typename T> void gpuComputeSnav(T *snav, T *dblist, T *blist, T *x, int *aii, int *ai, int *aj,
        int *ti, int *mask, int ncoeff, int ntypes, int nperdim, int ijnum, int quadraticflag)
{        
    int N1 = ncoeff;
    int N2 = ncoeff*ijnum;
    int blockDim = 256;
    int gridDim = (N2 + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;        
    gpuKernelComputeSnav11<<<gridDim, blockDim>>>(snav, dblist, blist, x, aii, ai, aj, 
          ti, mask, ncoeff, ntypes, nperdim, N1, N2);

    if (quadraticflag)
        gpuKernelComputeSnav12<<<gridDim, blockDim>>>(snav, dblist, blist, x, aii, ai, aj, 
          ti, mask, ncoeff, ntypes, nperdim, N1, N2);
}
template void gpuComputeSnav(double*, double*, double*, double*, int*, int*, int*, int*, int*,
        int, int, int, int, int);
template void gpuComputeSnav(float*, float*, float*, float*, int*, int*, int*, int*, int*,
        int, int, int, int, int);

template <typename T> void gpuComputeSnav2(T *snav, T *dblist, T *blist, T *x, int *aii, int *ai, int *aj,
        int *ti, int *mask, int ncoeff, int ntypes, int nperdim, int ijnum, int quadraticflag)
{        
    int N1 = ncoeff;
    int N2 = ncoeff*ijnum;
    int blockDim = 256;
    int gridDim = (N2 + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;        
    gpuKernelComputeSnav21<<<gridDim, blockDim>>>(snav, dblist, blist, x, aii, ai, aj, 
          ti, mask, ncoeff, ntypes, nperdim, N1, N2);

    if (quadraticflag)
        gpuKernelComputeSnav22<<<gridDim, blockDim>>>(snav, dblist, blist, x, aii, ai, aj, 
          ti, mask, ncoeff, ntypes, nperdim, N1, N2);
}
template void gpuComputeSnav2(double*, double*, double*, double*, int*, int*, int*, int*, int*,
        int, int, int, int, int);
template void gpuComputeSnav2(float*, float*, float*, float*, int*, int*, int*, int*, int*,
        int, int, int, int, int);

template <typename T> __global__ void gpuKernelComputeBeta1(T *beta, T *bispectrum, T *coeffelem, int *ilist, int *map, int *type, 
        int inum, int ncoeff, int ncoeffall, int N1, int N2)
{    
    int idx = threadIdx.x + blockIdx.x * blockDim.x;
    while (idx < N2) {
        int icoeff = idx%N1;              
        int ii = (idx-icoeff)/N1;         
        int i = ilist[ii]; // index of atom i
        const int itype = type[i]; // element type of atom i
        const int ielem = map[itype];  // index of that element type
        beta[ii*ncoeff + icoeff] = coeffelem[icoeff+1+ielem*ncoeffall];
        idx += blockDim.x * gridDim.x;
    }
}

template <typename T> __global__ void gpuKernelComputeBeta2(T *beta, T *bispectrum, T *coeffelem, int *ilist, int *map, int *type, 
        int inum, int ncoeff, int ncoeffall, int N1, int N2)
{    
    int idx = threadIdx.x + blockIdx.x * blockDim.x;
    while (idx < N2) {
        int icoeff = idx%N1;              
        int ii = (idx-icoeff)/N1;         
        int i = ilist[ii]; // index of atom i
        const int itype = type[i]; // element type of atom i
        const int ielem = map[itype];  // index of that element type
        T bveci = bispectrum[ii*ncoeff + icoeff];        
        int k = ncoeff + 1 + icoeff*ncoeff - (icoeff-1)*icoeff/2;
        beta[ii*ncoeff + icoeff] = coeffelem[k+ielem*ncoeffall];
        for (int jcoeff = icoeff+1; jcoeff < ncoeff; jcoeff++) {
          T bvecj = bispectrum[ii*ncoeff + jcoeff];
          beta[ii*ncoeff + icoeff] += coeffelem[k+ielem*ncoeffall]*bvecj;
          beta[ii*ncoeff + jcoeff] += coeffelem[k+ielem*ncoeffall]*bveci;
          k++;
        }
        idx += blockDim.x * gridDim.x;
    }
}

template <typename T> void gpuComputeBeta2(T *beta, T *bispectrum, T *coeffelem, int *ilist, int *map, int *type, 
        int inum, int ncoeff, int ncoeffall, int quadraticflag)
{
    int N1 = ncoeff;
    int N2 = ncoeff*inum;
    int blockDim = 256;
    int gridDim = (N2 + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;        
    gpuKernelComputeBeta1<<<gridDim, blockDim>>>(beta, bispectrum, coeffelem, ilist,
            map, type, inum, ncoeff, ncoeffall, N1, N2);

    if (quadraticflag)
        gpuKernelComputeBeta2<<<gridDim, blockDim>>>(beta, bispectrum, coeffelem, ilist,
            map, type, inum, ncoeff, ncoeffall, N1, N2);
    }
template void gpuComputeBeta2(double*, double*, double*, int*, int*, int*,
        int, int, int, int);
template void gpuComputeBeta2(float*, float*, float*, int*, int*, int*,
        int, int, int, int);

template <typename T> __global__ void gpuKernelSnapTallyEnergyFull(T *eatom, T *bispectrum, 
        T *coeffelem, int *ilist, int *map, int *type, 
        int inum, int ncoeff, int ncoeffall, int quadraticflag)
{      
    int ii = threadIdx.x + blockIdx.x * blockDim.x;
    while (ii < inum) {        
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
        eatom[ii] += evdwl;
        ii += blockDim.x * gridDim.x;
    }
}
template <typename T> void gpuSnapTallyEnergyFull(T *eatom, T *bispectrum, T *coeffelem, int *ilist, int *map, int *type, 
        int inum, int ncoeff, int ncoeffall, int quadraticflag)
{
    int blockDim = 256;
    int gridDim = (inum + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;        
    gpuKernelSnapTallyEnergyFull<<<gridDim, blockDim>>>(eatom, bispectrum, coeffelem, ilist, map, type, 
                inum, ncoeff, ncoeffall, quadraticflag);
    
}
template void gpuSnapTallyEnergyFull(double*, double*, double*, int*, int*, int*,
        int, int, int, int);
template void gpuSnapTallyEnergyFull(float*, float*, float*, int*, int*, int*,
        int, int, int, int);

template <typename T> __global__ void gpuKernelSnapTallyForceFull(T *fatom, T *fij, int *ai, 
        int *aj, int ijnum)
{ 
    int k = threadIdx.x + blockIdx.x * blockDim.x;
    while (k < ijnum) {            
        int i = ai[k];        
        int j = aj[k];        
        T fx = fij[0+3*k];
        T fy = fij[1+3*k];
        T fz = fij[2+3*k];    
        atomicAdd(&fatom[0+3*i], fx);
        atomicAdd(&fatom[1+3*i], fy);
        atomicAdd(&fatom[2+3*i], fz);
        atomicAdd(&fatom[0+3*j], -fx);
        atomicAdd(&fatom[1+3*j], -fy);
        atomicAdd(&fatom[2+3*j], -fz);
        k += blockDim.x * gridDim.x;
    }
};
template <typename T> __global__ void gpuKernelSnapTallyForceI(T *fatom, T *fij, int *ai, int ijnum)
{ 
    int k = threadIdx.x + blockIdx.x * blockDim.x;
    while (k < ijnum) {            
        int i = ai[k];                    
        atomicAdd(&fatom[0+3*i], fij[0+3*k]);
        atomicAdd(&fatom[1+3*i], fij[1+3*k]);
        atomicAdd(&fatom[2+3*i], fij[2+3*k]);
        k += blockDim.x * gridDim.x;
    }
};
template <typename T> __global__ void gpuKernelSnapTallyForceJ(T *fatom, T *fij, int *aj, int ijnum)
{ 
    int k = threadIdx.x + blockIdx.x * blockDim.x;
    while (k < ijnum) {            
        int j = aj[k];        
        atomicAdd(&fatom[0+3*j], -fij[0+3*k]);
        atomicAdd(&fatom[1+3*j], -fij[1+3*k]);
        atomicAdd(&fatom[2+3*j], -fij[2+3*k]);
        k += blockDim.x * gridDim.x;
    }
};
template <typename T> void gpuSnapTallyForceFull(T *fatom, T *fij, int *ai, 
        int *aj, int *alist, int ijnum)
{
    int blockDim = 256;
    int gridDim = (ijnum + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;        
    gpuKernelSnapTallyForceFull<<<gridDim, blockDim>>>(fatom, fij, ai, aj, ijnum);
    //gpuKernelSnapTallyForceI<<<gridDim, blockDim>>>(fatom, fij, ai, ijnum);
    //gpuKernelSnapTallyForceJ<<<gridDim, blockDim>>>(fatom, fij, aj, ijnum);
}
template void gpuSnapTallyForceFull(double*, double*, int*, int*, int*, int);
template void gpuSnapTallyForceFull(float*, float*, int*, int*, int*, int);

template <typename T> __global__ void gpuKernelSnapTallyVirialFull(T *vatom, T *fij, T *rij, 
        int *ai, int *aj, int inum, int ijnum)
{ 
    int k = threadIdx.x + blockIdx.x * blockDim.x;
    while (k < ijnum) {                
        int i = ai[k];        
        int j = aj[k];        
        T factor = 0.5;
        T dx = -rij[0+3*k];
        T dy = -rij[1+3*k];
        T dz = -rij[2+3*k];    
        T fx = fij[0+3*k];
        T fy = fij[1+3*k];
        T fz = fij[2+3*k];    
        T v0 = factor*dx*fx;
        T v1 = factor*dy*fy;
        T v2 = factor*dz*fz;
        T v3 = factor*dx*fy;
        T v4 = factor*dx*fz;
        T v5 = factor*dy*fz;        
        atomicAdd(&vatom[0*inum+i], v0);
        atomicAdd(&vatom[1*inum+i], v1);
        atomicAdd(&vatom[2*inum+i], v2);
        atomicAdd(&vatom[3*inum+i], v3);
        atomicAdd(&vatom[4*inum+i], v4);
        atomicAdd(&vatom[5*inum+i], v5);
        atomicAdd(&vatom[0*inum+j], v0);
        atomicAdd(&vatom[1*inum+j], v1);
        atomicAdd(&vatom[2*inum+j], v2);
        atomicAdd(&vatom[3*inum+j], v3);
        atomicAdd(&vatom[4*inum+j], v4);
        atomicAdd(&vatom[5*inum+j], v5);
        k += blockDim.x * gridDim.x;
    }
};
template <typename T> void gpuSnapTallyVirialFull(T *vatom, T *fij, T *rij, int *ai, int *aj, int inum, int ijnum)
{ 
    int blockDim = 256;
    int gridDim = (ijnum + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;        
    gpuKernelSnapTallyVirialFull<<<gridDim, blockDim>>>(vatom, fij, rij, ai, aj, inum, ijnum);
}
template void gpuSnapTallyVirialFull(double*, double*, double*, int*, int*, int, int);
template void gpuSnapTallyVirialFull(float*, float*, float*, int*, int*, int, int);

template <typename T> __global__ void gpuKernelNeighPairList(int *pairnum, int *pairlist, T *x, T rcutsq, int *ilist, int *neighlist, 
        int *neighnum, int inum, int jnum, int dim)
{    
    int ii = threadIdx.x + blockIdx.x * blockDim.x;
    while (ii < inum) {        
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
        ii += blockDim.x * gridDim.x;
    }    
};
template <typename T> void gpuNeighPairList(int *pairnum, int *pairlist, T *x, T rcutsq, int *ilist, int *neighlist, 
        int *neighnum, int inum, int jnum, int dim)
{ 
    int blockDim = 256;
    int gridDim = (inum + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;        
    gpuKernelNeighPairList<<<gridDim, blockDim>>>(pairnum, pairlist, x, rcutsq, ilist, neighlist, 
                            neighnum, inum, jnum, dim);
}
template void gpuNeighPairList(int*, int*, double*, double, int*, int*, int*, int, int, int);
template void gpuNeighPairList(int*, int*, float*, float, int*, int*, int*, int, int, int);

template <typename T> __global__ void gpuKernelNeighPairList(int *pairnum, int *pairlist, T *x, 
        T *rcutsq, int *ilist, int *neighlist, int *neighnum, int inum, int jnum, int dim)
{    
    int ii = threadIdx.x + blockIdx.x * blockDim.x;
    while (ii < inum) {        
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
        ii += blockDim.x * gridDim.x;
    }    
};
template <typename T> void gpuNeighPairList(int *pairnum, int *pairlist, T *x, T *rcutsq, 
        int *ilist, int *neighlist, int *neighnum, int inum, int jnum, int dim)
{ 
    int blockDim = 256;
    int gridDim = (inum + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;        
    gpuKernelNeighPairList<<<gridDim, blockDim>>>(pairnum, pairlist, x, rcutsq, ilist, neighlist, 
                            neighnum, inum, jnum, dim);
}
template void gpuNeighPairList(int*, int*, double*, double*, int*, int*, int*, int, int, int);
template void gpuNeighPairList(int*, int*, float*, float*, int*, int*, int*, int, int, int);

template <typename T> __global__ void gpuKernelNeighPairList(int *pairnum, int *pairlist, T *x, T *rcutsq, int *ilist, int *neighlist, 
        int *neighnum, int *atomtype, int *alist, int inum, int jnum, int dim, int ntypes)
{    
    int ii = threadIdx.x + blockIdx.x * blockDim.x;
    while (ii < inum) {        
        int i = ilist[ii];       // atom i
        int itype = atomtype[i];
        int m = neighnum[i];     // number of neighbors around i                       
        int count = 0;              
        for (int l=0; l<m ; l++) {   // loop over each atom around atom i {
            int j = neighlist[l + jnum*i];         
            int jtype  = atomtype[alist[j]];        
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
        ii += blockDim.x * gridDim.x;
    }    
};
template <typename T> void gpuNeighPairList(int *pairnum, int *pairlist, T *x, T *rcutsq, int *ilist, 
        int *neighlist, int *neighnum, int *atomtype, int *alist, int inum, int jnum, int dim, int ntypes)
{ 
    int blockDim = 256;
    int gridDim = (inum + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;        
    gpuKernelNeighPairList<<<gridDim, blockDim>>>(pairnum, pairlist, x, rcutsq, ilist, neighlist, 
                            neighnum, atomtype, alist, inum, jnum, dim, ntypes);
}
template void gpuNeighPairList(int*, int*, double*, double*, int*, int*, int*, int*, int*, 
        int, int, int, int);
template void gpuNeighPairList(int*, int*, float*, float*, int*, int*, int*, int*, int*, 
        int, int, int, int);

template <typename T> __global__ void gpuKernelNeighPairs(T *xij, T *x, int *aii, int *ai, 
      int *aj, int *ti, int *tj, int *pairnum, int *pairlist, int *pairnumsum, int *ilist, 
      int *atomtype, int *alist, int inum, int jnum, int dim)
{        
    int ii = threadIdx.x + blockIdx.x * blockDim.x;
    while (ii < inum) {        
        int i = ilist[ii];       // atom i
        int itype = atomtype[i];
        int m = pairnum[ii];        // number of neighbors around i             
        int start = pairnumsum[ii];   
        for (int l=0; l<m ; l++) {   // loop over each atom around atom i
            int j = pairlist[l + jnum*ii];  // atom j              
            int k = start + l;                                     
            aii[k]       = ii;
            ai[k]        = i;
            aj[k]        = alist[j];          
            ti[k]        = itype;       
            tj[k]        = atomtype[alist[j]];        
            for (int d=0; d<dim; d++) 
                xij[k*dim+d]   = x[j*dim+d] -  x[i*dim+d];  // xj - xi            
        }
        ii += blockDim.x * gridDim.x;
    }    
};
template <typename T> void gpuNeighPairs(T *xij, T *x, int *aii, int *ai, int *aj,  
      int *ti, int *tj, int *pairnum, int *pairlist, int *pairnumsum, int *ilist, 
      int *atomtype, int *alist, int inum, int jnum, int dim)
{ 
    int blockDim = 256;
    int gridDim = (inum + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;        
    gpuKernelNeighPairs<<<gridDim, blockDim>>>(xij, x, aii, ai, aj,  ti, tj, pairnum, pairlist, pairnumsum, ilist, 
                        atomtype, alist, inum, jnum, dim);
}
template void gpuNeighPairs(double*, double*, int*, int*, int*, int*, int*, int*, int*, int*, int*, 
        int*, int*, int, int, int);
template void gpuNeighPairs(float*, float*, int*, int*, int*, int*, int*, int*, int*, int*, int*, 
        int*, int*, int, int, int);

#endif

