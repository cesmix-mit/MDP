/***************************************************************************                               
                    Molecular Dynamics Potentials (MDP)
                           CESMIX-MIT Project  
 
 Contributing authors: Ngoc-Cuong Nguyen (cuongng@mit.edu, exapde@gmail.com)
 ***************************************************************************/
  
#ifndef MDP_CPUFIXES
#define MDP_CPUFIXES

#include <math.h>
#include <algorithm>

//#include "cpuRandom.cpp"

using std::min;
using std::max;
        
// cpuFixSetForce, cpuFixLineForce, cpuFixPlaneForce, cpuFixAddForce, cpuFixAveForce, cpuFixDrag, 
// cpuFixWallReflect, cpuFixWallHarmonic, cpuFixWallLJ93, cpuFixWallLJ126, cpuFixWallLJ1043, cpuFixWallMorse 
        
template <typename T> void cpuFixSetForce(T *x, T *v, T *f, T *eatom, T *vatom, T *fparam, 
        int *iparam, int *ilist, int eflag_atom, int vflag_atom, int dim, int inum)        
{
    if (dim==3) {
        for (int ii = 0; ii < inum; ii++) {
            int i = ilist[ii]; 
            if (iparam[0]) f[i*dim+0] = fparam[0];
            if (iparam[1]) f[i*dim+1] = fparam[1];
            if (iparam[2]) f[i*dim+2] = fparam[2];            
        }
    } else {
        for (int ii = 0; ii < inum; ii++) {
            int i = ilist[ii]; 
            if (iparam[0]) f[i*dim+0] = fparam[0];
            if (iparam[1]) f[i*dim+1] = fparam[1];
        }        
    }
}
template void cpuFixSetForce(double *x, double *v, double *f, double *eatom, double *vatom, 
        double *fparam, int *iparam, int *ilist, int eflag_atom, int vflag_atom, int dim, int inum);
template void cpuFixSetForce(float *x, float *v, float *f, float *eatom, float *vatom, 
        float *fparam, int *iparam, int *ilist, int eflag_atom, int vflag_atom, int dim, int inum);
         
template <typename T> void cpuFixLineForce(T *x, T *v, T *f, T *eatom, T *vatom, T *fparam, 
        int *iparam, int *ilist, int eflag_atom, int vflag_atom, int dim, int inum)        
{
    if (dim==3) {
        for (int ii = 0; ii < inum; ii++) {
            int i = ilist[ii]; 
            T dot = f[i*dim+0]*fparam[0] + f[i*dim+1]*fparam[1] + f[i*dim+2]*fparam[2];
            f[i*dim+0] = dot * fparam[0];
            f[i*dim+1] = dot * fparam[1];
            f[i*dim+2] = dot * fparam[1];            
        }
    } else {
        for (int ii = 0; ii < inum; ii++) {
            int i = ilist[ii]; 
            T dot = f[i*dim+0]*fparam[0] + f[i*dim+1]*fparam[1];
            f[i*dim+0] = dot * fparam[0];
            f[i*dim+1] = dot * fparam[1];
        }        
    }
}
template void cpuFixLineForce(double *x, double *v, double *f, double *eatom, double *vatom, 
        double *fparam, int *iparam, int *ilist, int eflag_atom, int vflag_atom, int dim, int inum);
template void cpuFixLineForce(float *x, float *v, float *f, float *eatom, float *vatom, 
        float *fparam, int *iparam, int *ilist, int eflag_atom, int vflag_atom, int dim, int inum);

template <typename T> void cpuFixPlaneForce(T *x, T *v, T *f, T *eatom, T *vatom, T *fparam, 
        int *iparam, int *ilist, int eflag_atom, int vflag_atom, int dim, int inum)        
{
    if (dim==3) {
        for (int ii = 0; ii < inum; ii++) {
            int i = ilist[ii]; 
            T dot = f[i*dim+0]*fparam[0] + f[i*dim+1]*fparam[1] + f[i*dim+2]*fparam[2];
            f[i*dim+0] -= dot * fparam[0];
            f[i*dim+1] -= dot * fparam[1];
            f[i*dim+2] -= dot * fparam[1];            
        }
    } else {
        for (int ii = 0; ii < inum; ii++) {
            int i = ilist[ii]; 
            T dot = f[i*dim+0]*fparam[0] + f[i*dim+1]*fparam[1];
            f[i*dim+0] -= dot * fparam[0];
            f[i*dim+1] -= dot * fparam[1];
        }        
    }
}
template void cpuFixPlaneForce(double *x, double *v, double *f, double *eatom, double *vatom, 
        double *fparam, int *iparam, int *ilist, int eflag_atom, int vflag_atom, int dim, int inum);
template void cpuFixPlaneForce(float *x, float *v, float *f, float *eatom, float *vatom, 
        float *fparam, int *iparam, int *ilist, int eflag_atom, int vflag_atom, int dim, int inum);

template <typename T> void cpuFixAddForce(T *x, T *v, T *f, T *eatom, T *vatom, T *fparam, 
        int *iparam, int *ilist, int eflag_atom, int vflag_atom, int dim, int inum)
{    
    if (dim==3) {
        for (int ii = 0; ii < inum; ii++) {
            int i = ilist[ii]; 
            T y0 = x[i*dim+0];
            T y1 = x[i*dim+1];
            T y2 = x[i*dim+2];  
            f[i*dim+0] += fparam[0];
            f[i*dim+1] += fparam[1];
            f[i*dim+2] += fparam[2];            
            if (eflag_atom)
                eatom[i] = -(fparam[0]*y0 + fparam[1]*y1 + fparam[2]*y2);
            if (vflag_atom) {
              vatom[i*6+0] += fparam[0] * y0;
              vatom[i*6+1] += fparam[1] * y1;
              vatom[i*6+2] += fparam[2] * y2;
              vatom[i*6+3] += fparam[0] * y1;
              vatom[i*6+4] += fparam[0] * y2;
              vatom[i*6+5] += fparam[1] * y2;
            }                    
        }
    } else {
        for (int ii = 0; ii < inum; ii++) {
            int i = ilist[ii]; 
            T y0 = x[i*dim+0];
            T y1 = x[i*dim+1];
            f[i*dim+0] += fparam[0];
            f[i*dim+1] += fparam[1];
            if (eflag_atom)
                eatom[i] = -(fparam[0]*y0 + fparam[1]*y1);
            if (vflag_atom) {
              vatom[i*6+0] += fparam[0] * y0;
              vatom[i*6+1] += fparam[1] * y1;
              vatom[i*6+2] += 0.0;
              vatom[i*6+3] += fparam[0] * y1;
              vatom[i*6+4] += 0.0;
              vatom[i*6+5] += 0.0;
            }                    
        }        
    }    
}
template void cpuFixAddForce(double *x, double *v, double *f, double *eatom, double *vatom, 
        double *fparam, int *iparam, int *ilist, int eflag_atom, int vflag_atom, int dim, int inum);
template void cpuFixAddForce(float *x, float *v, float *f, float *eatom, float *vatom, 
        float *fparam, int *iparam, int *ilist, int eflag_atom, int vflag_atom, int dim, int inum);

template <typename T> void cpuFixAddForce(T *x, T *v, T *f, T *eatom, T *vatom, T *fparam, T *box,
        int *iparam, int *ilist, int *image, int triclinic, int eflag_atom, int vflag_atom, int dim, int inum)
{    
    if (dim==3) {
        for (int ii = 0; ii < inum; ii++) {
            int i = ilist[ii]; 
            T y0, y1, y2;  
            f[i*dim+0] += fparam[0];
            f[i*dim+1] += fparam[1];
            f[i*dim+2] += fparam[2];            
//             int xbox = (image[i] & IMGMASK) - IMGMAX;
//             int ybox = (image[i] >> IMGBITS & IMGMASK) - IMGMAX;
//             int zbox = (image[i] >> IMG2BITS) - IMGMAX;            
            if (triclinic == 0) {
                y0 = x[0+i*dim] + image[i*dim+0]*box[0];
                y1 = x[1+i*dim] + image[i*dim+1]*box[1];
                y2 = x[2+i*dim] + image[i*dim+2]*box[2];
            } else {
                y0 = x[0+i*dim] + box[0]*image[i*dim+0] + box[5]*image[i*dim+1] + box[4]*image[i*dim+2];
                y1 = x[1+i*dim] + box[1]*image[i*dim+1] + box[3]*image[i*dim+2];
                y2 = x[2+i*dim] + box[2]*image[i*dim+2];
            }                                   
            if (eflag_atom)
                eatom[i] = -(fparam[0]*y0 + fparam[1]*y1 + fparam[2]*y2);
            if (vflag_atom) {
              vatom[i*6+0] += fparam[0] * y0;
              vatom[i*6+1] += fparam[1] * y1;
              vatom[i*6+2] += fparam[2] * y2;
              vatom[i*6+3] += fparam[0] * y1;
              vatom[i*6+4] += fparam[0] * y2;
              vatom[i*6+5] += fparam[1] * y2;
            }                    
        }
    } else {
        for (int ii = 0; ii < inum; ii++) {
            int i = ilist[ii]; 
            T y0, y1;  
            f[i*dim+0] += fparam[0];
            f[i*dim+1] += fparam[1];
//             int xbox = (image[i] & IMGMASK) - IMGMAX;
//             int ybox = (image[i] >> IMGBITS & IMGMASK) - IMGMAX;
            if (triclinic == 0) {
                y0 = x[0+i*dim] + image[i*dim+0]*box[0];
                y1 = x[1+i*dim] + image[i*dim+1]*box[1];
            } else {
                y0 = x[0+i*dim] + box[0]*image[i*dim+0] + box[5]*image[i*dim+1];
                y1 = x[1+i*dim] + box[1]*image[i*dim+1];
            }                                                         
            if (eflag_atom)
                eatom[i] = -(fparam[0]*y0 + fparam[1]*y1);
            if (vflag_atom) {
              vatom[i*6+0] += fparam[0] * y0;
              vatom[i*6+1] += fparam[1] * y1;
              vatom[i*6+2] += 0.0;
              vatom[i*6+3] += fparam[0] * y1;
              vatom[i*6+4] += 0.0;
              vatom[i*6+5] += 0.0;
            }                    
        }        
    }    
}
template void cpuFixAddForce(double *x, double *v, double *f, double *eatom, double *vatom, 
        double *fparam, double *h, int *iparam, int *ilist, int *image, int triclinic, 
        int eflag_atom, int vflag_atom, int dim, int inum);
template void cpuFixAddForce(float *x, float *v, float *f, float *eatom, float *vatom, 
        float *fparam, float *h, int *iparam, int *ilist, int *image, int triclinic, int eflag_atom, 
        int vflag_atom, int dim, int inum);

// template <typename T> T cpuFixAddForce(T *x, T *f, T *eatom, T *virial, T *vatom, T *h, T *fvalue, 
//         int *ilist, int *image, int *fdim, int evflag, int triclinic, 
//         int vflag_global, int vflag_atom, int dim, int inum)
// {    
//     T unwrap[3];  
//     T v[6];    
//     for (int ii = 0; ii < 6; ii++)
//         v[ii] = 0.0;
//     T energy = 0.0;    
//     if (dim==3) {
//         for (int ii = 0; ii < inum; ii++) {
//             int i = ilist[ii]; 
//             if (fdim[0]) f[i*dim+0] += fvalue[0];
//             if (fdim[1]) f[i*dim+1] += fvalue[1];
//             if (fdim[2]) f[i*dim+2] += fvalue[2];            
//             cpuUnmap(unwrap, &x[i*dim], h, image[i], triclinic, dim);
//             eatom[i] = -(fvalue[0]*unwrap[0] + fvalue[1]*unwrap[1] + fvalue[2]*unwrap[2]);
//             energy += eatom[i];            
//             if (evflag) {
//               v[0] = fvalue[0] * unwrap[0];
//               v[1] = fvalue[1] * unwrap[1];
//               v[2] = fvalue[2] * unwrap[2];
//               v[3] = fvalue[0] * unwrap[1];
//               v[4] = fvalue[0] * unwrap[2];
//               v[5] = fvalue[1] * unwrap[2];
//               cpuVirialTally(virial, vatom, v, i, vflag_global, vflag_atom);
//             }                    
//         }
//     } else {
//         for (int ii = 0; ii < inum; ii++) {
//             int i = ilist[ii]; 
//             if (fdim[0]) f[i*dim+0] += fvalue[0];
//             if (fdim[1]) f[i*dim+1] += fvalue[1];            
//             cpuUnmap(unwrap, &x[i*dim], h, image[i], triclinic, dim);
//             //energy -= fvalue[0]*unwrap[0] + fvalue[1]*unwrap[1];
//             eatom[i] = -(fvalue[0]*unwrap[0] + fvalue[1]*unwrap[1]);
//             energy += eatom[i];            
//             if (evflag) {
//               v[0] = fvalue[0] * unwrap[0];
//               v[1] = fvalue[1] * unwrap[1];
//               v[3] = fvalue[0] * unwrap[1];                            
//               cpuVirialTally(virial, vatom, v, i, vflag_global, vflag_atom);
//             }                    
//         }        
//     }
//     
//     return energy;
// }
// template double cpuFixAddForce(double *x, double *f, double *eatom, double *virial, double *vatom, 
//         double *h, double *fvalue, int *ilist, int *image, int *fdim, int evflag, int triclinic, 
//         int vflag_global, int vflag_atom, int dim, int inum);
// template float cpuFixAddForce(float *x, float *f, float *eatom, float *virial, float *vatom, 
//         float *h, float *fvalue, int *ilist, int *image, int *fdim, int evflag, int triclinic, 
//         int vflag_global, int vflag_atom, int dim, int inum);

template <typename T> void cpuFixAveForce(T *x, T *v, T *f, T *eatom, T *vatom, T *fparam, 
        int *iparam, int *ilist, int eflag_atom, int vflag_atom, int dim, int inum)
{
  T foriginal[4];
  foriginal[0] = foriginal[1] = foriginal[2] = foriginal[3] = 0.0;

  for (int ii = 0; ii < inum; ii++) {
      int i = ilist[ii]; 
      foriginal[0] += f[i*dim+0];
      foriginal[1] += f[i*dim+1];
      if (dim==3) foriginal[2] += f[i*dim+2];
      foriginal[3] += 1.0;
  }
    
  //MPI_Allreduce(foriginal,foriginal_all,4,MPI_DOUBLE,MPI_SUM,world);
  
  int ncount = static_cast<int> (foriginal[3]);
  if (ncount == 0) return;
    
  foriginal[0] = foriginal[0]/ncount + fparam[0];
  foriginal[1] = foriginal[1]/ncount + fparam[1];
  if (dim==3) 
    foriginal[2] = foriginal[2]/ncount + fparam[2];
  
  for (int ii = 0; ii < inum; ii++) {
        int i = ilist[ii]; 
        if (iparam[0]) f[i*dim+0] = foriginal[0];
        if (iparam[1]) f[i*dim+1] = foriginal[1];
        if (iparam[2]) f[i*dim+2] = foriginal[2];    
   }
}
template void cpuFixAveForce(double *x, double *v, double *f, double *eatom, double *vatom, 
        double *fparam, int *iparam, int *ilist, int eflag_atom, int vflag_atom, int dim, int inum);
template void cpuFixAveForce(float *x, float *v, float *f, float *eatom, float *vatom, 
        float *fparam, int *iparam, int *ilist, int eflag_atom, int vflag_atom, int dim, int inum);

// template <typename T> void cpuFixAveForce(T *x, T *f, T *foriginal, T *fvalue, 
//         int *ilist, int *fdim, int dim, int inum)
// {
//   //T foriginal[4];
//   foriginal[0] = foriginal[1] = foriginal[2] = foriginal[3] = 0.0;
// 
//   for (int ii = 0; ii < inum; ii++) {
//       int i = ilist[ii]; 
//       foriginal[0] += f[i*dim+0];
//       foriginal[1] += f[i*dim+1];
//       if (dim==3) foriginal[2] += f[i*dim+2];
//       foriginal[3] += 1.0;
//   }
//     
//   //MPI_Allreduce(foriginal,foriginal_all,4,MPI_DOUBLE,MPI_SUM,world);
//   
//   int ncount = static_cast<int> (foriginal[3]);
//   if (ncount == 0) return;
//     
//   foriginal[0] = foriginal[0]/ncount + fvalue[0];
//   foriginal[1] = foriginal[1]/ncount + fvalue[1];
//   if (dim==3) 
//     foriginal[2] = foriginal[2]/ncount + fvalue[2];
//   
//   for (int ii = 0; ii < inum; ii++) {
//         int i = ilist[ii]; 
//         if (fdim[0]) f[i*dim+0] = foriginal[0];
//         if (fdim[1]) f[i*dim+1] = foriginal[1];
//         if (fdim[2]) f[i*dim+2] = foriginal[2];    
//    }
// }
// template void cpuFixAveForce(double *x, double *f, double *foriginal, double *fvalue,
//         int *ilist, int *fdim, int dim, int inum);
// template void cpuFixAveForce(float *x, float *f, float *foriginal, float *fvalue, 
//         int *ilist, int *fdim, int dim, int inum);

template <typename T> void cpuFixDragForce(T *x, T *v, T *f, T *eatom, T *vatom, T *fparam, T *box,
        int *iparam, int *ilist, int *pbc, int triclinic, int eflag_atom, int vflag_atom, int dim, int inum)
{
  if (dim==3) {
      for (int ii = 0; ii < inum; ii++) {
          int i = ilist[ii]; 
          T f_mag = fparam[0];
          T delta = fparam[1];
          T dp[3]; 
          dp[0] = x[i*dim+0] - fparam[2];
          dp[1] = x[i*dim+1] - fparam[3];
          dp[2] = x[i*dim+2] - fparam[4];
          if (!iparam[0]) dp[0] = 0.0;
          if (!iparam[1]) dp[1] = 0.0;
          if (!iparam[2]) dp[2] = 0.0;      
          cpuMinimumImage(dp, box, pbc, triclinic, dim);
          T r = sqrt(dp[0]*dp[0] + dp[1]*dp[1] + dp[2]*dp[2]);
          if (r > delta) {
            T prefactor = f_mag/r;
            T fx = prefactor*dp[0];
            T fy = prefactor*dp[1];
            T fz = prefactor*dp[2];
            f[i*dim+0] -= fx;
            f[i*dim+1] -= fy;
            f[i*dim+2] -= fz;
          }
      }
  } else {
      for (int ii = 0; ii < inum; ii++) {
          int i = ilist[ii]; 
          T f_mag = fparam[0];
          T delta = fparam[1];          
          T dp[3]; 
          dp[0] = x[i*dim+0] - fparam[2];
          dp[1] = x[i*dim+1] - fparam[3];
          if (!iparam[0]) dp[0] = 0.0;
          if (!iparam[1]) dp[1] = 0.0;
          cpuMinimumImage(dp, box, pbc, triclinic, dim);
          T r = sqrt(dp[0]*dp[0] + dp[1]*dp[1]);
          if (r > delta) {
            T prefactor = f_mag/r;
            T fx = prefactor*dp[0];
            T fy = prefactor*dp[1];
            f[i*dim+0] -= fx;
            f[i*dim+1] -= fy;
          }
      }      
  }  
}
template void cpuFixDragForce(double *x, double *v, double *f, double *eatom, double *vatom, 
        double *fparam, double *box, int *iparam, int *ilist, int *pbc, int triclinic, 
        int eflag_atom, int vflag_atom, int dim, int inum);
template void cpuFixDragForce(float *x, float *v, float *f, float *eatom, float *vatom, 
        float *fparam, float *box, int *iparam, int *ilist, int *pbc, int triclinic, int eflag_atom, 
        int vflag_atom, int dim, int inum);


// template <typename T> void cpuFixDragForce(T *x, T *f, T *box, T *p, T *ftotal, T *param,
//         int *ilist, int *pbc, int *dimflag, int triclinic, int dim, int inum)
// {
//   ftotal[0] = ftotal[1] = ftotal[2] = 0.0;
// 
//   if (dim==3) {
//       for (int ii = 0; ii < inum; ii++) {
//           int i = ilist[ii]; 
//           T f_mag = param[0];
//           T delta = param[1];
//           T dp[3]; 
//           dp[0] = x[i*dim+0] - p[0];
//           dp[1] = x[i*dim+1] - p[1];
//           dp[2] = x[i*dim+2] - p[2];
//           if (!dimflag[0]) dp[0] = 0.0;
//           if (!dimflag[1]) dp[1] = 0.0;
//           if (!dimflag[2]) dp[2] = 0.0;      
//           cpuMinimumImage(dp, box, pbc, triclinic, dim);
//           T r = sqrt(dp[0]*dp[0] + p[1]*dp[1] + p[2]*dp[2]);
//           if (r > delta) {
//             T prefactor = f_mag/r;
//             T fx = prefactor*dp[0];
//             T fy = prefactor*dp[1];
//             T fz = prefactor*dp[2];
//             f[i*dim+0] -= fx;
//             f[i*dim+1] -= fy;
//             f[i*dim+2] -= fz;
//             ftotal[0] -= fx;
//             ftotal[1] -= fy;
//             ftotal[2] -= fz;
//           }
//       }
//   } else {
//       for (int ii = 0; ii < inum; ii++) {
//           int i = ilist[ii]; 
//           T f_mag = param[0];
//           T delta = param[1];          
//           T dp[3]; 
//           dp[0] = x[i*dim+0] - p[0];
//           dp[1] = x[i*dim+1] - p[1];
//           if (!dimflag[0]) dp[0] = 0.0;
//           if (!dimflag[1]) dp[1] = 0.0;
//           cpuMinimumImage(dp, box, pbc, triclinic, dim);
//           T r = sqrt(dp[0]*dp[0] + p[1]*dp[1]);
//           if (r > delta) {
//             T prefactor = f_mag/r;
//             T fx = prefactor*dp[0];
//             T fy = prefactor*dp[1];
//             f[i*dim+0] -= fx;
//             f[i*dim+1] -= fy;
//             ftotal[0] -= fx;
//             ftotal[1] -= fy;
//           }
//       }      
//   }  
// }
// template void cpuFixDragForce(double *x, double *f, double *box, double *p, double *ftotal, double *param, 
//         int *ilist, int *pbc, int *dimflag, int triclinic, int dim, int inum);
// template void cpuFixDragForce(float *x, float *f, float *box, float *p, float *ftotal, float *param,
//         int *ilist, int *pbc, int *dimflag, int triclinic, int dim, int inum);
        
template <typename T> void cpuFixWallReflect(T *x, T *v, T *f, T *eatom, T *vatom, T *fparam, 
        int *iparam, int *ilist, int eflag_atom, int vflag_atom, int dim, int inum)
{  
  for (int ii = 0; ii < inum; ii++) {
      int i = ilist[ii]; 
      int direction = iparam[0] / 2;
      int side = iparam[0] % 2;      
      T coord = fparam[0];
      if (side == 0) {
        if (x[i*dim+direction] < coord) {
          x[i*dim+direction] = coord + (coord - x[i*dim+direction]);
          v[i*dim+direction] = -v[i*dim+direction];
        }
      } else {
        if (x[i*dim+direction] > coord) {
          x[i*dim+direction] = coord - (x[i*dim+direction] - coord);
          v[i*dim+direction] = -v[i*dim+direction];
        }
      }    
  }  
}
template void cpuFixWallReflect(double *x, double *v, double *f, double *eatom, double *vatom, 
        double *fparam, int *iparam, int *ilist, int eflag_atom, int vflag_atom, int dim, int inum);
template void cpuFixWallReflect(float *x, float *v, float *f, float *eatom, float *vatom, 
        float *fparam, int *iparam, int *ilist, int eflag_atom, int vflag_atom, int dim, int inum);
 
template <typename T> void cpuFixWallHarmonic(T *x, T *v, T *f, T *eatom, T *vatom, T *fparam, 
        int *iparam, int *ilist, int eflag_atom, int vflag_atom, int dim, int inum)
{
  for (int ii = 0; ii < inum; ii++) {
      int i = ilist[ii]; 
      int direction = iparam[0] / 2;
      int side = iparam[0] % 2;     
      if (side == 0) side = -1;
      T delta,dr,fwall,vn;
      T coord = fparam[0];
      T cutoff = fparam[1];
      T epsilon = fparam[2];
      if (side < 0) delta = x[i*dim+direction] - coord;
      else delta = coord - x[i*dim+direction];
      if (delta >= cutoff) continue;
      if (delta <= 0.0) continue;      
      dr = cutoff-delta;
      fwall = side * 2.0*epsilon*dr;
      f[i*dim+direction] -= fwall;
      
      if (eflag_atom) 
          eatom[i] += epsilon*dr*dr;      
      if (vflag_atom) {
        if (side < 0) vn = -fwall*delta;
        else vn = fwall*delta;
        vatom[i*6+direction] += vn; 
     }
  }  
}
template void cpuFixWallHarmonic(double *x, double *v, double *f, double *eatom, double *vatom, 
        double *fparam, int *iparam, int *ilist, int eflag_atom, int vflag_atom, int dim, int inum);
template void cpuFixWallHarmonic(float *x, float *v, float *f, float *eatom, float *vatom, 
        float *fparam, int *iparam, int *ilist, int eflag_atom, int vflag_atom, int dim, int inum);

template <typename T> void cpuFixWallLJ93(T *x, T *v, T *f, T *eatom, T *vatom, T *fparam, 
        int *iparam, int *ilist, int eflag_atom, int vflag_atom, int dim, int inum)
{
  for (int ii = 0; ii < inum; ii++) {
      int i = ilist[ii]; 
      int direction = iparam[0] / 2;
      int side = iparam[0] % 2;     
      if (side == 0) side = -1;
      T delta,fwall,vn;    
      T coord = fparam[0];
      T cutoff = fparam[1];
      T offset = fparam[2];
      T coeff1 = fparam[3];
      T coeff2 = fparam[4];
      T coeff3 = fparam[5];
      T coeff4 = fparam[6];            
      if (side < 0) delta = x[i*dim+direction] - coord;
      else delta = coord - x[i*dim+direction];
      if (delta >= cutoff) continue;
      if (delta <= 0.0) continue;      
      T rinv = 1.0/delta;
      T r2inv = rinv*rinv;
      T r4inv = r2inv*r2inv;
      T r10inv = r4inv*r4inv*r2inv;
      fwall = side * (coeff1*r10inv - coeff2*r4inv);
      f[i*dim+direction] -= fwall;
      
      if (eflag_atom) 
        eatom[i] += coeff3*r4inv*r4inv*rinv - coeff4*r2inv*rinv - offset;
            
      if (vflag_atom) {
        if (side < 0) vn = -fwall*delta;
        else vn = fwall*delta;
        vatom[i*6+direction] += vn; 
      }    
  }  
}
template void cpuFixWallLJ93(double *x, double *v, double *f, double *eatom, double *vatom, 
        double *fparam, int *iparam, int *ilist, int eflag_atom, int vflag_atom, int dim, int inum);
template void cpuFixWallLJ93(float *x, float *v, float *f, float *eatom, float *vatom, 
        float *fparam, int *iparam, int *ilist, int eflag_atom, int vflag_atom, int dim, int inum);

template <typename T> void cpuFixWallLJ126(T *x, T *v, T *f, T *eatom, T *vatom, T *fparam, 
        int *iparam, int *ilist, int eflag_atom, int vflag_atom, int dim, int inum)
{
  for (int ii = 0; ii < inum; ii++) {
      int i = ilist[ii];     
      int direction = iparam[0] / 2;
      int side = iparam[0] % 2;     
      if (side == 0) side = -1;
      T delta,fwall,vn;      
      T coord = fparam[0];
      T cutoff = fparam[1];
      T offset = fparam[2];
      T coeff1 = fparam[3];
      T coeff2 = fparam[4];
      T coeff3 = fparam[5];
      T coeff4 = fparam[6];            
      if (side < 0) delta = x[i*dim+direction] - coord;
      else delta = coord - x[i*dim+direction];
      if (delta >= cutoff) continue;
      if (delta <= 0.0) continue;      
      T rinv = 1.0/delta;
      T r2inv = rinv*rinv;
      T r6inv = r2inv*r2inv*r2inv;
      fwall = side * r6inv*(coeff1*r6inv - coeff2) * rinv;
      f[i*dim+direction] -= fwall;
      
      if (eflag_atom) 
        eatom[i] += r6inv*(coeff3*r6inv - coeff4) - offset;
            
      if (vflag_atom) {
        if (side < 0) vn = -fwall*delta;
        else vn = fwall*delta;
        vatom[i*6+direction] += vn;         
      }    
  }    
}
template void cpuFixWallLJ126(double *x, double *v, double *f, double *eatom, double *vatom, 
        double *fparam, int *iparam, int *ilist, int eflag_atom, int vflag_atom, int dim, int inum);
template void cpuFixWallLJ126(float *x, float *v, float *f, float *eatom, float *vatom, 
        float *fparam, int *iparam, int *ilist, int eflag_atom, int vflag_atom, int dim, int inum);

template <typename T> void cpuFixWallLJ1043(T *x, T *v, T *f, T *eatom, T *vatom, T *fparam, 
        int *iparam, int *ilist, int eflag_atom, int vflag_atom, int dim, int inum)
{
  for (int ii = 0; ii < inum; ii++) {
      int i = ilist[ii]; 
      int direction = iparam[0] / 2;
      int side = iparam[0] % 2;     
      if (side == 0) side = -1;
      T delta,fwall,vn;      
      T coord = fparam[0];
      T cutoff = fparam[1];
      T offset = fparam[2];
      T coeff1 = fparam[3];
      T coeff2 = fparam[4];
      T coeff3 = fparam[5];
      T coeff4 = fparam[6];            
      T coeff5 = fparam[7];            
      T coeff6 = fparam[8];            
      T coeff7 = fparam[9];            
      if (side < 0) delta = x[i*dim+direction] - coord;
      else delta = coord - x[i*dim+direction];
      if (delta >= cutoff) continue;
      if (delta <= 0.0) continue;      
      T rinv = 1.0/delta;
      T r2inv = rinv*rinv;
      T r4inv = r2inv*r2inv;
      T r10inv = r4inv*r4inv*r2inv;
      fwall = side * (coeff5*r10inv*rinv - coeff6*r4inv*rinv -
        coeff7*pow(delta+coeff4,-4.0));
      f[i*dim+direction] -= fwall;
      
      if (eflag_atom)
        eatom[i] += coeff1*r10inv - coeff2*r4inv - coeff3*pow(delta+coeff4,-3.0) - offset;
            
      if (vflag_atom) {
        if (side < 0) vn = -fwall*delta;
        else vn = fwall*delta;
        vatom[i*6+direction] += vn;         
      }    
  }  
}
template void cpuFixWallLJ1043(double *x, double *v, double *f, double *eatom, double *vatom, 
        double *fparam, int *iparam, int *ilist, int eflag_atom, int vflag_atom, int dim, int inum);
template void cpuFixWallLJ1043(float *x, float *v, float *f, float *eatom, float *vatom, 
        float *fparam, int *iparam, int *ilist, int eflag_atom, int vflag_atom, int dim, int inum);

template <typename T> void cpuFixWallMorse(T *x, T *v, T *f, T *eatom, T *vatom, T *fparam, 
        int *iparam, int *ilist, int eflag_atom, int vflag_atom, int dim, int inum)
{
  for (int ii = 0; ii < inum; ii++) {
      int i = ilist[ii]; 
      int direction = iparam[0] / 2;
      int side = iparam[0] % 2;     
      if (side == 0) side = -1;
      T delta,dr,fwall,vn, dexp;     
      T coord = fparam[0];
      T cutoff = fparam[1];
      T offset = fparam[2];
      T epsilon = fparam[3];
      T alpha = fparam[4];
      T coeff1 = fparam[5];
      if (side < 0) delta = x[i*dim+direction] - coord;
      else delta = coord - x[i*dim+direction];
      if (delta >= cutoff) continue;
      if (delta <= 0.0) continue;      
      dr = cutoff-delta;
      dexp = exp(-alpha * dr);
      fwall = side * coeff1 * (dexp*dexp - dexp) / delta;
      f[i*dim+direction] -= fwall;
      
      if (vflag_atom)
        eatom[i] += epsilon * (dexp*dexp - 2.0*dexp) - offset;      
      
      if (vflag_atom) {
        if (side < 0) vn = -fwall*delta;
        else vn = fwall*delta;
        vatom[i*6+direction] += vn;         
      }    
  }  
}
template void cpuFixWallMorse(double *x, double *v, double *f, double *eatom, double *vatom, 
        double *fparam, int *iparam, int *ilist, int eflag_atom, int vflag_atom, int dim, int inum);
template void cpuFixWallMorse(float *x, float *v, float *f, float *eatom, float *vatom, 
        float *fparam, int *iparam, int *ilist, int eflag_atom, int vflag_atom, int dim, int inum);


// template <typename T> int cpuFixWallHarmonic(T *x, T *v, T *f, T *eatom, T *virial, T *vatom, T *ewall, 
//         T *param, int *ilist, int evflag, int which, int vflag_global, int vflag_atom, int m, int dim, int inum)
// {
//   int direction = which / 2;
//   int side = which % 2;
//   if (side == 0) side = -1;
//   int onflag = 0;
// 
//   for (int ii = 0; ii < inum; ii++) {
//       int i = ilist[ii]; 
//       T delta,dr,fwall,vn;
//       T coord = param[0];
//       T cutoff = param[1];
//       T epsilon = param[2];
//       if (side < 0) delta = x[i*dim+direction] - coord;
//       else delta = coord - x[i*dim+direction];
//       if (delta >= cutoff) continue;
//       if (delta <= 0.0) {
//         onflag = 1;
//         continue;
//       }
//       dr = cutoff-delta;
//       fwall = side * 2.0*epsilon*dr;
//       f[i*dim+direction] -= fwall;
//       ewall[0] += epsilon*dr*dr;
//       ewall[m+1] += fwall;
// 
//       if (evflag) {
//         if (side < 0) vn = -fwall*delta;
//         else vn = fwall*delta;
//         //v_tally(direction, i, vn);
//         if (vflag_global) 
//             virial[direction] += vn;
//         if (vflag_atom) 
//             vatom[i*6+direction] += vn; 
//      }
//   }
//   
//   return onflag;    
// }
// template int cpuFixWallHarmonic(double *x, double *v, double *f, double *eatom, double *virial, 
//         double *vatom, double *ewall, double *param,  int *ilist, int evflag, int which, int vflag_global,
//         int vflag_atom, int m, int dim, int inum);
// template int cpuFixWallHarmonic(float *x, float *v, float *f, float *eatom, float *virial, 
//         float *vatom, float *ewall, float *param, int *ilist, int evflag, int which, int vflag_global,
//         int vflag_atom, int m, int dim, int inum);

// template <typename T> int cpuFixWallLJ93(T *x, T *v, T *f, T *eatom, T *virial, T *vatom, T *ewall, 
//         T *param, int *ilist, int evflag, int vflag_global, int vflag_atom, int which, int m, int dim, int inum)
// {
//   int direction = which / 2;
//   int side = which % 2;
//   if (side == 0) side = -1;
//   int onflag = 0;
// 
//   for (int ii = 0; ii < inum; ii++) {
//       int i = ilist[ii]; 
//       T delta,fwall,vn;    
//       T coord = param[0];
//       T cutoff = param[1];
//       T offset = param[2];
//       T coeff1 = param[3];
//       T coeff2 = param[4];
//       T coeff3 = param[5];
//       T coeff4 = param[6];            
//       if (side < 0) delta = x[i*dim+direction] - coord;
//       else delta = coord - x[i*dim+direction];
//       if (delta >= cutoff) continue;
//       if (delta <= 0.0) {
//         onflag = 1;
//         continue;
//       }
//       T rinv = 1.0/delta;
//       T r2inv = rinv*rinv;
//       T r4inv = r2inv*r2inv;
//       T r10inv = r4inv*r4inv*r2inv;
//       fwall = side * (coeff1*r10inv - coeff2*r4inv);
//       f[i*dim+direction] -= fwall;
//       eatom[i] = coeff3*r4inv*r4inv*rinv - coeff4*r2inv*rinv - offset;
//       ewall[0] += eatom[i];
//       ewall[m+1] += fwall;
//       
//       if (evflag) {
//         if (side < 0) vn = -fwall*delta;
//         else vn = fwall*delta;
//         //v_tally(direction, i, vn);
//         if (vflag_global) 
//             virial[direction] += vn;
//         if (vflag_atom) 
//             vatom[i*6+direction] += vn;         
//       }    
//   }
//   
//   return onflag;
// }
// template int cpuFixWallLJ93(double *x, double *v, double *f, double *eatom, double *virial, 
//         double *vatom, double *ewall, double *param, int *ilist, int evflag, int vflag_global, 
//         int vflag_atom, int which, int m, int dim, int inum);
// template int cpuFixWallLJ93(float *x, float *v, float *f, float *eatom, float *virial, 
//         float *vatom, float *ewall, float *param, int *ilist, int evflag, int vflag_global, 
//         int vflag_atom, int which, int m, int dim, int inum);
// 
// template <typename T> int cpuFixWallLJ126(T *x, T *v, T *f, T *eatom, T *virial, T *vatom, T *ewall, 
//         T *param, int *ilist, int evflag, int vflag_global, int vflag_atom, int which, int m, int dim, int inum)        
// {
//   int direction = which / 2;
//   int side = which % 2;
//   if (side == 0) side = -1;
//   int onflag = 0;
// 
//   for (int ii = 0; ii < inum; ii++) {
//       int i = ilist[ii];     
//       T delta,fwall,vn;      
//       T coord = param[0];
//       T cutoff = param[1];
//       T offset = param[2];
//       T coeff1 = param[3];
//       T coeff2 = param[4];
//       T coeff3 = param[5];
//       T coeff4 = param[6];            
//       if (side < 0) delta = x[i*dim+direction] - coord;
//       else delta = coord - x[i*dim+direction];
//       if (delta >= cutoff) continue;
//       if (delta <= 0.0) {
//         onflag = 1;
//         continue;
//       }
//       T rinv = 1.0/delta;
//       T r2inv = rinv*rinv;
//       T r6inv = r2inv*r2inv*r2inv;
//       fwall = side * r6inv*(coeff1*r6inv - coeff2) * rinv;
//       f[i*dim+direction] -= fwall;
//       eatom[i] = r6inv*(coeff3*r6inv - coeff4) - offset;
//       ewall[0] += eatom[i];
//       ewall[m+1] += fwall;
//             
//       if (evflag) {
//         if (side < 0) vn = -fwall*delta;
//         else vn = fwall*delta;
//         if (vflag_global) 
//             virial[direction] += vn;
//         if (vflag_atom) 
//             vatom[i*6+direction] += vn;         
//       }    
//   }
//   
//   return onflag;
// }
// template int cpuFixWallLJ126(double *x, double *v, double *f, double *eatom, double *virial, 
//         double *vatom, double *ewall, double *param, int *ilist, int evflag, int vflag_global, 
//         int vflag_atom, int which, int m, int dim, int inum);
// template int cpuFixWallLJ126(float *x, float *v, float *f, float *eatom, float *virial, 
//         float *vatom, float *ewall, float *param, int *ilist, int evflag, int vflag_global, 
//         int vflag_atom, int which, int m, int dim, int inum);
// 
// template <typename T> int cpuFixWallLJ1043(T *x, T *v, T *f, T *eatom, T *virial, T *vatom, T *ewall, 
//         T *param, int *ilist, int evflag, int vflag_global, int vflag_atom, int which, int m, int dim, int inum)        
// {
//   int direction = which / 2;
//   int side = which % 2;
//   if (side == 0) side = -1;
//   int onflag = 0;
// 
//   for (int ii = 0; ii < inum; ii++) {
//       int i = ilist[ii]; 
//       T delta,fwall,vn;      
//       T coord = param[0];
//       T cutoff = param[1];
//       T offset = param[2];
//       T coeff1 = param[3];
//       T coeff2 = param[4];
//       T coeff3 = param[5];
//       T coeff4 = param[6];            
//       T coeff5 = param[7];            
//       T coeff6 = param[8];            
//       T coeff7 = param[9];            
//       if (side < 0) delta = x[i*dim+direction] - coord;
//       else delta = coord - x[i*dim+direction];
//       if (delta >= cutoff) continue;
//       if (delta <= 0.0) {
//         onflag = 1;
//         continue;
//       }
//       T rinv = 1.0/delta;
//       T r2inv = rinv*rinv;
//       T r4inv = r2inv*r2inv;
//       T r10inv = r4inv*r4inv*r2inv;
//       fwall = side * (coeff5*r10inv*rinv - coeff6*r4inv*rinv -
//         coeff7*pow(delta+coeff4,-4.0));
//       f[i*dim+direction] -= fwall;
//       eatom[i] = coeff1*r10inv - coeff2*r4inv - coeff3*pow(delta+coeff4,-3.0) - offset;
//       ewall[0] += eatom[i];
//       ewall[m+1] += fwall;
//             
//       if (evflag) {
//         if (side < 0) vn = -fwall*delta;
//         else vn = fwall*delta;
//         if (vflag_global) 
//             virial[direction] += vn;
//         if (vflag_atom) 
//             vatom[i*6+direction] += vn;         
//       }    
//   }
//   
//   return onflag;
// }
// template int cpuFixWallLJ1043(double *x, double *v, double *f, double *eatom, double *virial, 
//         double *vatom, double *ewall, double *param, int *ilist, int evflag, int vflag_global, 
//         int vflag_atom, int which, int m, int dim, int inum);
// template int cpuFixWallLJ1043(float *x, float *v, float *f, float *eatom, float *virial, 
//         float *vatom, float *ewall, float *param, int *ilist, int evflag, int vflag_global, 
//         int vflag_atom, int which, int m, int dim, int inum);
// 
// template <typename T> int cpuFixWallMorse(T *x, T *v, T *f, T *eatom, T *virial, T *vatom, T *ewall, 
//         T *param, int *ilist, int evflag, int vflag_global, int vflag_atom, int which, int m, int dim, int inum)        
// {
//   int direction = which / 2;
//   int side = which % 2;
//   if (side == 0) side = -1;
//   int onflag = 0;
// 
//   for (int ii = 0; ii < inum; ii++) {
//       int i = ilist[ii]; 
//       T delta,dr,fwall,vn, dexp;     
//       T coord = param[0];
//       T cutoff = param[1];
//       T offset = param[2];
//       T epsilon = param[3];
//       T alpha = param[4];
//       T coeff1 = param[5];
//       if (side < 0) delta = x[i*dim+direction] - coord;
//       else delta = coord - x[i*dim+direction];
//       if (delta >= cutoff) continue;
//       if (delta <= 0.0) {
//         onflag = 1;
//         continue;
//       }      
//       dr = cutoff-delta;
//       dexp = exp(-alpha * dr);
//       fwall = side * coeff1 * (dexp*dexp - dexp) / delta;
//       eatom[i] = epsilon * (dexp*dexp - 2.0*dexp) - offset;
//       ewall[0] += eatom[i];
//       f[i*dim+direction] -= fwall;
//       ewall[m+1] += fwall;
//       
//       if (evflag) {
//         if (side < 0) vn = -fwall*delta;
//         else vn = fwall*delta;
//         if (vflag_global) 
//             virial[direction] += vn;
//         if (vflag_atom) 
//             vatom[i*6+direction] += vn;         
//       }    
//   }
//   
//   return onflag;     
// }
// template int cpuFixWallMorse(double *x, double *v, double *f, double *eatom, double *virial, 
//         double *vatom, double *ewall, double *param, int *ilist, int evflag, int vflag_global, 
//         int vflag_atom, int which, int m, int dim, int inum);
// template int cpuFixWallMorse(float *x, float *v, float *f, float *eatom, float *virial, 
//         float *vatom, float *ewall, float *param, int *ilist, int evflag, int vflag_global, 
//         int vflag_atom, int which, int m, int dim, int inum);

// template <typename T> int cpuFixMomentum(T *x, T *v, T *box, T *p, T *ftotal, T f_mag, T delta,
//         int *ilist, int *pbc, int xflag, int yflag, int zflag,
//         int triclinic, int dim, int inum)
// {
//   T ekin_old,ekin_new;
//   ekin_old = ekin_new = 0.0;
// 
//   // compute kinetic energy before momentum removal, if needed
//   if (rescale) {
//     T ke=0.0;    
//     for (int i = 0; i < inum; i++)
//        if (ilist[i] & groupbit)
//           ke +=  mass[type[i]] *
//             (v[i][0]*v[i][0] + v[i][1]*v[i][1] + v[i][2]*v[i][2]);    
//     //MPI_Allreduce(&ke,&ekin_old,1,MPI_DOUBLE,MPI_SUM,world);
//      ekin_old = ke;
//   }
// 
//   if (linear) {
//     // adjust velocities by vcm to zero linear momentum
//     // only adjust a component if flag is set
//     for (int i = 0; i < inum; i++)
//       if (ilist[i] & groupbit) {
//         if (xflag) v[i][0] -= vcm[0];
//         if (yflag) v[i][1] -= vcm[1];
//         if (zflag) v[i][2] -= vcm[2];
//       }
//   }
// 
//   if (angular) {
//     // adjust velocities to zero omega
//     // vnew_i = v_i - w x r_i
//     // must use unwrapped coords to compute r_i correctly
//     double dx,dy,dz;
//     double unwrap[3];
//     for (int i = 0; i < inum; i++)
//       if (ilist[i] & groupbit) {
//         domain->unmap(x[i],image[i],unwrap);
//         dx = unwrap[0] - xcm[0];
//         dy = unwrap[1] - xcm[1];
//         dz = unwrap[2] - xcm[2];
//         v[i][0] -= omega[1]*dz - omega[2]*dy;
//         v[i][1] -= omega[2]*dx - omega[0]*dz;
//         v[i][2] -= omega[0]*dy - omega[1]*dx;
//       }
//   }
// 
//   // compute kinetic energy after momentum removal, if needed
// 
//   if (rescale) {
//     double ke=0.0, factor=1.0;
//     for (int i = 0; i < inum; i++)
//       if (ilist[i] & groupbit)
//         ke +=  mass[type[i]] *
//           (v[i][0]*v[i][0] + v[i][1]*v[i][1] + v[i][2]*v[i][2]);
// 
//     ekin_new = ke;
//     //MPI_Allreduce(&ke,&ekin_new,1,MPI_DOUBLE,MPI_SUM,world);
// 
//     if (ekin_new != 0.0) factor = sqrt(ekin_old/ekin_new);
//     for (int i = 0; i < inum; i++) {
//       if (ilist[i] & groupbit) {
//         v[i][0] *= factor;
//         v[i][1] *= factor;
//         v[i][2] *= factor;
//       }
//     }
//   }
// }
// 
// void FixGravity::post_force(int /*vflag*/)
// {
//   // update gravity due to variables
// 
//   if (varflag != CONSTANT) {
//     modify->clearstep_compute();
//     if (mstyle == EQUAL) magnitude = input->variable->compute_equal(mvar);
//     if (vstyle == EQUAL) vert = input->variable->compute_equal(vvar);
//     if (pstyle == EQUAL) phi = input->variable->compute_equal(pvar);
//     if (tstyle == EQUAL) theta = input->variable->compute_equal(tvar);
//     if (xstyle == EQUAL) xdir = input->variable->compute_equal(xvar);
//     if (ystyle == EQUAL) ydir = input->variable->compute_equal(yvar);
//     if (zstyle == EQUAL) zdir = input->variable->compute_equal(zvar);
//     modify->addstep_compute(update->ntimestep + 1);
// 
//     set_acceleration();
//   }
// 
//   // just exit if application of force is disabled
// 
//   if (disable) return;
// 
//   // apply gravity force to each particle
// 
//   double **x = atom->x;
//   double **f = atom->f;
//   double *rmass = atom->rmass;
//   double *mass = atom->mass;
//   int *ilist = atom->ilist;
//   int *type = atom->type;
//   int inum = atom->inum;
//   double massone;
// 
//   eflag = 0;
//   egrav = 0.0;
// 
//   if (rmass) {
//     for (int i = 0; i < inum; i++)
//       if (ilist[i] & groupbit) {
//         massone = rmass[i];
//         f[i][0] += massone*xacc;
//         f[i][1] += massone*yacc;
//         f[i][2] += massone*zacc;
//         egrav -= massone * (xacc*x[i][0] + yacc*x[i][1] + zacc*x[i][2]);
//       }
//   } else {
//     for (int i = 0; i < inum; i++)
//       if (ilist[i] & groupbit) {
//         massone = mass[type[i]];
//         f[i][0] += massone*xacc;
//         f[i][1] += massone*yacc;
//         f[i][2] += massone*zacc;
//         egrav -= massone * (xacc*x[i][0] + yacc*x[i][1] + zacc*x[i][2]);
//       }
//   }
// }
// 
// /* ---------------------------------------------------------------------- */
// 
// void FixGravity::post_force_respa(int vflag, int ilevel, int /*iloop*/)
// {
//   if (ilevel == ilevel_respa) post_force(vflag);
// }
// 
// /* ---------------------------------------------------------------------- */
// 
// void FixGravity::set_acceleration()
// {
//   if (style == CHUTE || style == SPHERICAL) {
//     if (style == CHUTE) {
//       phi = 0.0;
//       theta = 180.0 - vert;
//     }
//     if (domain->dimension == 3) {
//       xgrav = sin(degree2rad * theta) * cos(degree2rad * phi);
//       ygrav = sin(degree2rad * theta) * sin(degree2rad * phi);
//       zgrav = cos(degree2rad * theta);
//     } else {
//       xgrav = sin(degree2rad * theta);
//       ygrav = cos(degree2rad * theta);
//       zgrav = 0.0;
//     }
//   } else if (style == VECTOR) {
//     if (domain->dimension == 3) {
//       double length = sqrt(xdir*xdir + ydir*ydir + zdir*zdir);
//       xgrav = xdir/length;
//       ygrav = ydir/length;
//       zgrav = zdir/length;
//     } else {
//       double length = sqrt(xdir*xdir + ydir*ydir);
//       xgrav = xdir/length;
//       ygrav = ydir/length;
//       zgrav = 0.0;
//     }
//   }
// 
//   gvec[0] = xacc = magnitude*xgrav;
//   gvec[1] = yacc = magnitude*ygrav;
//   gvec[2] = zacc = magnitude*zgrav;
// }
// 
// /* ----------------------------------------------------------------------
//    potential energy in gravity field
// ------------------------------------------------------------------------- */
// 
// double FixGravity::compute_scalar()
// {
//   // only sum across procs one time
// 
//   if (eflag == 0) {
//     MPI_Allreduce(&egrav,&egrav_all,1,MPI_DOUBLE,MPI_SUM,world);
//     eflag = 1;
//   }
//   return egrav_all;
// }
// 
// /* ----------------------------------------------------------------------
//    extract current gravity direction vector
// ------------------------------------------------------------------------- */
// 
// void *FixGravity::extract(const char *name, int &dim)
// {
//   if (strcmp(name,"gvec") == 0) {
//     dim = 1;
//     return (void *) gvec;
//   }
//   return nullptr;
// }
// 
// void FixHeat::end_of_step()
// {
//   int i;
//   double heat,ke,massone;
//   double vsub[3],vcm[3];
// 
//   double **x = atom->x;
//   double **v = atom->v;
//   int *ilist = atom->ilist;
//   int inum = atom->inum;
//   int *type = atom->type;
//   double *mass = atom->mass;
//   double *rmass = atom->rmass;
// 
//   // reallocate per-atom arrays if necessary
// 
//   if (hstyle == ATOM && atom->nmax > maxatom) {
//     maxatom = atom->nmax;
//     memory->destroy(vheat);
//     memory->destroy(vscale);
//     memory->create(vheat,maxatom,"heat:vheat");
//     memory->create(vscale,maxatom,"heat:vscale");
//   }
// 
//   // evaluate variable
// 
//   if (hstyle != CONSTANT) {
//     modify->clearstep_compute();
//     if (hstyle == EQUAL) heat_input = input->variable->compute_equal(hvar);
//     else input->variable->compute_atom(hvar,igroup,vheat,1,0);
//     modify->addstep_compute(update->ntimestep + nevery);
//   }
// 
//   // vcm = center-of-mass velocity of scaled atoms
// 
//   if (iregion < 0) {
//     ke = group->ke(igroup)*force->ftm2v;
//     group->vcm(igroup,masstotal,vcm);
//   } else {
//     masstotal = group->mass(igroup,iregion);
//     if (masstotal == 0.0) error->all(FLERR,"Fix heat group has no atoms");
//     ke = group->ke(igroup,iregion)*force->ftm2v;
//     group->vcm(igroup,masstotal,vcm,iregion);
//   }
//   double vcmsq = vcm[0]*vcm[0] + vcm[1]*vcm[1] + vcm[2]*vcm[2];
// 
//   // add heat via scale factor on velocities for CONSTANT and EQUAL cases
//   // scale = velocity scale factor to accomplish eflux change in energy
//   // vsub = velocity subtracted from each atom to preserve momentum
//   // overall KE cannot go negative
// 
//   Region *region = nullptr;
//   if (iregion >= 0) {
//     region = domain->regions[iregion];
//     region->prematch();
//   }
// 
//   if (hstyle != ATOM) {
//     heat = heat_input*nevery*update->dt*force->ftm2v;
//     double escale =
//       (ke + heat - 0.5*vcmsq*masstotal)/(ke - 0.5*vcmsq*masstotal);
//     if (escale < 0.0) error->all(FLERR,"Fix heat kinetic energy went negative");
//     scale = sqrt(escale);
//     vsub[0] = (scale-1.0) * vcm[0];
//     vsub[1] = (scale-1.0) * vcm[1];
//     vsub[2] = (scale-1.0) * vcm[2];
// 
//     if (iregion < 0) {
//       for (i = 0; i < inum; i++)
//         if (ilist[i] & groupbit) {
//           v[i][0] = scale*v[i][0] - vsub[0];
//           v[i][1] = scale*v[i][1] - vsub[1];
//           v[i][2] = scale*v[i][2] - vsub[2];
//         }
//     } else {
//       for (int i = 0; i < inum; i++)
//         if (ilist[i] & groupbit && region->match(x[i][0],x[i][1],x[i][2])) {
//           v[i][0] = scale*v[i][0] - vsub[0];
//           v[i][1] = scale*v[i][1] - vsub[1];
//           v[i][2] = scale*v[i][2] - vsub[2];
//         }
//     }
// 
//   // add heat via per-atom scale factor on velocities for ATOM case
//   // vscale = velocity scale factor to accomplish eflux change in energy
//   // vsub = velocity subtracted from each atom to preserve momentum
//   // KE of an atom cannot go negative
// 
//   } else {
//     vsub[0] = vsub[1] = vsub[2] = 0.0;
//     if (iregion < 0) {
//       for (i = 0; i < inum; i++) {
//         if (ilist[i] & groupbit) {
//           heat = vheat[i]*nevery*update->dt*force->ftm2v;
//           vscale[i] =
//             (ke + heat - 0.5*vcmsq*masstotal)/(ke - 0.5*vcmsq*masstotal);
//           if (vscale[i] < 0.0)
//             error->all(FLERR,
//                        "Fix heat kinetic energy of an atom went negative");
//           scale = sqrt(vscale[i]);
//           if (rmass) massone = rmass[i];
//           else massone = mass[type[i]];
//           vsub[0] += (scale-1.0) * v[i][0]*massone;
//           vsub[1] += (scale-1.0) * v[i][1]*massone;
//           vsub[2] += (scale-1.0) * v[i][2]*massone;
//         }
//       }
// 
//       vsub[0] /= masstotal;
//       vsub[1] /= masstotal;
//       vsub[2] /= masstotal;
// 
//       for (i = 0; i < inum; i++)
//         if (ilist[i] & groupbit) {
//           scale = sqrt(vscale[i]);
//           v[i][0] = scale*v[i][0] - vsub[0];
//           v[i][1] = scale*v[i][1] - vsub[1];
//           v[i][2] = scale*v[i][2] - vsub[2];
//         }
// 
//     } else {
//       for (i = 0; i < inum; i++) {
//         if (ilist[i] & groupbit && region->match(x[i][0],x[i][1],x[i][2])) {
//           heat = vheat[i]*nevery*update->dt*force->ftm2v;
//           vscale[i] =
//             (ke + heat - 0.5*vcmsq*masstotal)/(ke - 0.5*vcmsq*masstotal);
//           if (vscale[i] < 0.0)
//             error->all(FLERR,
//                        "Fix heat kinetic energy of an atom went negative");
//           scale = sqrt(vscale[i]);
//           if (rmass) massone = rmass[i];
//           else massone = mass[type[i]];
//           vsub[0] += (scale-1.0) * v[i][0]*massone;
//           vsub[1] += (scale-1.0) * v[i][1]*massone;
//           vsub[2] += (scale-1.0) * v[i][2]*massone;
//         }
//       }
// 
//       vsub[0] /= masstotal;
//       vsub[1] /= masstotal;
//       vsub[2] /= masstotal;
// 
//       for (i = 0; i < inum; i++)
//         if (ilist[i] & groupbit && region->match(x[i][0],x[i][1],x[i][2])) {
//           scale = sqrt(vscale[i]);
//           v[i][0] = scale*v[i][0] - vsub[0];
//           v[i][1] = scale*v[i][1] - vsub[1];
//           v[i][2] = scale*v[i][2] - vsub[2];
//         }
//     }
//   }
// }
// 
// /* ---------------------------------------------------------------------- */
// 
// double FixHeat::compute_scalar()
// {
//   double average_scale = scale;
//   if (hstyle == ATOM) {
//     if (!vscale) return 1.0;
//     double scale_sum = 0.0;
//     int ncount = 0;
//     int *ilist = atom->ilist;
//     double **x = atom->x;
//     int inum = atom->inum;
//     if (iregion < 0) {
//       for (int i = 0; i < inum; i++) {
//         if (ilist[i] & groupbit) {
//           scale_sum += sqrt(vscale[i]);
//           ncount++;
//         }
//       }
//     } else {
//       Region *region = domain->regions[iregion];
//       region->prematch();
//       for (int i = 0; i < inum; i++) {
//         if (ilist[i] & groupbit && region->match(x[i][0],x[i][1],x[i][2])) {
//           scale_sum += sqrt(vscale[i]);
//           ncount++;
//         }
//       }
//     }
//     double scale_sum_all = 0.0;
//     int ncount_all = 0;
//     MPI_Allreduce(&scale_sum,&scale_sum_all,1,MPI_DOUBLE,MPI_SUM,world);
//     MPI_Allreduce(&ncount,&ncount_all,1,MPI_INT,MPI_SUM,world);
//     if (ncount_all == 0) average_scale = 0.0;
//     else average_scale = scale_sum_all/static_cast<double>(ncount_all);
//   }
//   return average_scale;
// }
// 
// /* ----------------------------------------------------------------------
//    memory usage of local atom-based arrays
// ------------------------------------------------------------------------- */
// 
// double FixHeat::memory_usage()
// {
//   double bytes = 0.0;
//   if (hstyle == ATOM) bytes = atom->nmax*2 * sizeof(double);
//   return bytes;
// }
// 
// void FixExternal::post_force(int vflag)
// {
//   bigint ntimestep = update->ntimestep;
// 
//   int eflag = eflag_caller;
//   ev_init(eflag,vflag);
// 
//   // invoke the callback in driver program
//   // it will fill fexternal with forces
// 
//   if (mode == PF_CALLBACK && ntimestep % ncall == 0)
//     (this->callback)(ptr_caller,update->ntimestep,
//                      atom->inum,atom->tag,atom->x,fexternal);
// 
//   // add forces from fexternal to atoms in group
// 
//   if (ntimestep % napply == 0) {
//     double **f = atom->f;
//     int *ilist = atom->ilist;
//     int inum = atom->inum;
// 
//     for (int i = 0; i < inum; i++)
//       if (ilist[i] & groupbit) {
//         f[i][0] += fexternal[i][0];
//         f[i][1] += fexternal[i][1];
//         f[i][2] += fexternal[i][2];
//       }
//   }
// }

#endif

