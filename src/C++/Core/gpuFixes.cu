/***************************************************************************                               
                    Molecular Dynamics Potentials (MDP)
                           CESMIX-MIT Project  
 
 Contributing authors: Ngoc-Cuong Nguyen (cuongng@mit.edu, exapde@gmail.com)
 ***************************************************************************/
               
#ifndef MDP_GPUFIXES
#define MDP_GPUFIXES

#include <math.h>
#include <algorithm>

//#include "gpuRandom.cpp"

using std::min;
using std::max;
        
// gpuFixSetForce, gpuFixLineForce, gpuFixPlaneForce, gpuFixAddForce, gpuFixAveForce, gpuFixDrag, 
// gpuFixWallReflect, gpuFixWallHarmonic, gpuFixWallLJ93, gpuFixWallLJ126, gpuFixWallLJ1043, gpuFixWallMorse 
        
template <typename T> __global__ void gpuKernelFixSetForce(T *x, T *v, T *f, T *eatom, T *vatom, T *fparam, 
        int *iparam, int *ilist, int eflag_atom, int vflag_atom, int dim, int inum)        
{
    int ii = threadIdx.x + blockIdx.x * blockDim.x;
    while (ii < inum) {
        int i = ilist[ii]; 
        if (iparam[0]) f[i*dim+0] = fparam[0];
        if (iparam[1]) f[i*dim+1] = fparam[1];
        if (dim==3) {
            if (iparam[2]) f[i*dim+2] = fparam[2];            
        }
        ii += blockDim.x * gridDim.x;
    }
}
template <typename T> void gpuFixSetForce(T *x, T *v, T *f, T *eatom, T *vatom, T *fparam, 
        int *iparam, int *ilist, int eflag_atom, int vflag_atom, int dim, int inum)        
{        
    int blockDim = 256;
    int gridDim = (inum + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;
    gpuKernelFixSetForce<<<gridDim, blockDim>>>(x, v, f, eatom, vatom, fparam, 
            iparam, ilist, eflag_atom, vflag_atom, dim, inum);
}
template void gpuFixSetForce(double *x, double *v, double *f, double *eatom, double *vatom, 
        double *fparam, int *iparam, int *ilist, int eflag_atom, int vflag_atom, int dim, int inum);
template void gpuFixSetForce(float *x, float *v, float *f, float *eatom, float *vatom, 
        float *fparam, int *iparam, int *ilist, int eflag_atom, int vflag_atom, int dim, int inum);
         
template <typename T> __global__ void gpuKernelFixLineForce(T *x, T *v, T *f, T *eatom, T *vatom, T *fparam, 
        int *iparam, int *ilist, int eflag_atom, int vflag_atom, int dim, int inum)        
{    
    int ii = threadIdx.x + blockIdx.x * blockDim.x;
    while (ii < inum) {
        int i = ilist[ii];         
        T dot = f[i*dim+0]*fparam[0] + f[i*dim+1]*fparam[1];
        if (dim==3) dot += f[i*dim+2]*fparam[2];
        f[i*dim+0] = dot * fparam[0];
        f[i*dim+1] = dot * fparam[1];
        if (dim==3) f[i*dim+2] = dot * fparam[2];
        ii += blockDim.x * gridDim.x;
    }            
}
template <typename T> void gpuFixLineForce(T *x, T *v, T *f, T *eatom, T *vatom, T *fparam, 
        int *iparam, int *ilist, int eflag_atom, int vflag_atom, int dim, int inum)        
{        
    int blockDim = 256;
    int gridDim = (inum + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;
    gpuKernelFixLineForce<<<gridDim, blockDim>>>(x, v, f, eatom, vatom, fparam, 
            iparam, ilist, eflag_atom, vflag_atom, dim, inum);
}
template void gpuFixLineForce(double *x, double *v, double *f, double *eatom, double *vatom, 
        double *fparam, int *iparam, int *ilist, int eflag_atom, int vflag_atom, int dim, int inum);
template void gpuFixLineForce(float *x, float *v, float *f, float *eatom, float *vatom, 
        float *fparam, int *iparam, int *ilist, int eflag_atom, int vflag_atom, int dim, int inum);

template <typename T> __global__ void gpuKernelFixPlaneForce(T *x, T *v, T *f, T *eatom, T *vatom, T *fparam, 
        int *iparam, int *ilist, int eflag_atom, int vflag_atom, int dim, int inum)        
{    
    int ii = threadIdx.x + blockIdx.x * blockDim.x;
    while (ii < inum) {
        int i = ilist[ii];         
        T dot = f[i*dim+0]*fparam[0] + f[i*dim+1]*fparam[1];
        if (dim==3) dot += f[i*dim+2]*fparam[2];
        f[i*dim+0] -= dot * fparam[0];
        f[i*dim+1] -= dot * fparam[1];
        if (dim==3) f[i*dim+2] -= dot * fparam[2];
        ii += blockDim.x * gridDim.x;
    }            
}
template <typename T> void gpuFixPlaneForce(T *x, T *v, T *f, T *eatom, T *vatom, T *fparam, 
        int *iparam, int *ilist, int eflag_atom, int vflag_atom, int dim, int inum)        
{        
    int blockDim = 256;
    int gridDim = (inum + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;
    gpuKernelFixPlaneForce<<<gridDim, blockDim>>>(x, v, f, eatom, vatom, fparam, 
            iparam, ilist, eflag_atom, vflag_atom, dim, inum);
}
template void gpuFixPlaneForce(double *x, double *v, double *f, double *eatom, double *vatom, 
        double *fparam, int *iparam, int *ilist, int eflag_atom, int vflag_atom, int dim, int inum);
template void gpuFixPlaneForce(float *x, float *v, float *f, float *eatom, float *vatom, 
        float *fparam, int *iparam, int *ilist, int eflag_atom, int vflag_atom, int dim, int inum);

template <typename T> __global__ void gpuKernelFixAddForce3D(T *x, T *v, T *f, T *eatom, T *vatom, T *fparam, 
        int *iparam, int *ilist, int eflag_atom, int vflag_atom, int dim, int inum)
{    
    int ii = threadIdx.x + blockIdx.x * blockDim.x;
    while (ii < inum) {
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
        ii += blockDim.x * gridDim.x;
    }
}

template <typename T> __global__ void gpuKernelFixAddForce2D(T *x, T *v, T *f, T *eatom, T *vatom, T *fparam, 
        int *iparam, int *ilist, int eflag_atom, int vflag_atom, int dim, int inum)
{    
    int ii = threadIdx.x + blockIdx.x * blockDim.x;
    while (ii < inum) {
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
        ii += blockDim.x * gridDim.x;
    }            
}
template <typename T> void gpuFixAddForce(T *x, T *v, T *f, T *eatom, T *vatom, T *fparam, 
        int *iparam, int *ilist, int eflag_atom, int vflag_atom, int dim, int inum)        
{        
    int blockDim = 256;
    int gridDim = (inum + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;
    
    if (dim==3)
        gpuKernelFixAddForce3D<<<gridDim, blockDim>>>(x, v, f, eatom, vatom, fparam, 
                iparam, ilist, eflag_atom, vflag_atom, dim, inum);
    else
        gpuKernelFixAddForce2D<<<gridDim, blockDim>>>(x, v, f, eatom, vatom, fparam, 
                iparam, ilist, eflag_atom, vflag_atom, dim, inum);
}
template void gpuFixAddForce(double *x, double *v, double *f, double *eatom, double *vatom, 
        double *fparam, int *iparam, int *ilist, int eflag_atom, int vflag_atom, int dim, int inum);
template void gpuFixAddForce(float *x, float *v, float *f, float *eatom, float *vatom, 
        float *fparam, int *iparam, int *ilist, int eflag_atom, int vflag_atom, int dim, int inum);

template <typename T> __global__ void gpuKernelFixAddForce3D(T *x, T *v, T *f, T *eatom, T *vatom, T *fparam, T *box,
        int *iparam, int *ilist, int *image, int triclinic, int eflag_atom, int vflag_atom, int dim, int inum)
{        
    int ii = threadIdx.x + blockIdx.x * blockDim.x;
    while (ii < inum) {
        int i = ilist[ii]; 
        T y0, y1, y2;  
        f[i*dim+0] += fparam[0];
        f[i*dim+1] += fparam[1];
        f[i*dim+2] += fparam[2];            
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
        ii += blockDim.x * gridDim.x;
    }
}

template <typename T> __global__ void gpuKernelFixAddForce2D(T *x, T *v, T *f, T *eatom, T *vatom, T *fparam, T *box,
        int *iparam, int *ilist, int *image, int triclinic, int eflag_atom, int vflag_atom, int dim, int inum)
{    
    int ii = threadIdx.x + blockIdx.x * blockDim.x;
    while (ii < inum) {
        int i = ilist[ii]; 
        T y0, y1;  
        f[i*dim+0] += fparam[0];
        f[i*dim+1] += fparam[1];
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
        ii += blockDim.x * gridDim.x;
    }        
}
template <typename T> void gpuFixAddForce(T *x, T *v, T *f, T *eatom, T *vatom, T *fparam, T *box,
        int *iparam, int *ilist, int *image, int triclinic, int eflag_atom, int vflag_atom, int dim, int inum)
{        
    int blockDim = 256;
    int gridDim = (inum + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;
    
    if (dim==3)
        gpuKernelFixAddForce3D<<<gridDim, blockDim>>>(x, v, f, eatom, vatom, fparam, box,
                iparam, ilist, image, triclinic, eflag_atom, vflag_atom, dim, inum);
    else
        gpuKernelFixAddForce2D<<<gridDim, blockDim>>>(x, v, f, eatom, vatom, fparam, box,
                iparam, ilist, image, triclinic, eflag_atom, vflag_atom, dim, inum);
}
template void gpuFixAddForce(double *x, double *v, double *f, double *eatom, double *vatom, 
        double *fparam, double *box, int *iparam, int *ilist, int *image, int triclinic, 
        int eflag_atom, int vflag_atom, int dim, int inum);
template void gpuFixAddForce(float *x, float *v, float *f, float *eatom, float *vatom, 
        float *fparam, float *box, int *iparam, int *ilist, int *image, int triclinic, int eflag_atom, 
        int vflag_atom, int dim, int inum);

// template <typename T> void gpuFixAveForce(T *x, T *v, T *f, T *eatom, T *vatom, T *fparam, 
//         int *iparam, int *ilist, int eflag_atom, int vflag_atom, int dim, int inum)
// {
//   T foriginal[3];
//   foriginal[0] = foriginal[1] = foriginal[2] = 0.0;
// 
//   for (int ii = 0; ii < inum; ii++) {
//       int i = ilist[ii]; 
//       foriginal[0] += f[i*dim+0];
//       foriginal[1] += f[i*dim+1];
//       if (dim==3) foriginal[2] += f[i*dim+2];
//   }
//     
//   //MPI_Allreduce(foriginal,foriginal_all,4,MPI_DOUBLE,MPI_SUM,world);
//       
//   foriginal[0] = foriginal[0]/inum + fparam[0];
//   foriginal[1] = foriginal[1]/inum + fparam[1];
//   if (dim==3) 
//     foriginal[2] = foriginal[2]/inum + fparam[2];
//   
//   for (int ii = 0; ii < inum; ii++) {
//         int i = ilist[ii]; 
//         if (iparam[0]) f[i*dim+0] = foriginal[0];
//         if (iparam[1]) f[i*dim+1] = foriginal[1];
//         if (iparam[2]) f[i*dim+2] = foriginal[2];    
//    }
// }
// template void gpuFixAveForce(double *x, double *v, double *f, double *eatom, double *vatom, 
//         double *fparam, int *iparam, int *ilist, int eflag_atom, int vflag_atom, int dim, int inum);
// template void gpuFixAveForce(float *x, float *v, float *f, float *eatom, float *vatom, 
//         float *fparam, int *iparam, int *ilist, int eflag_atom, int vflag_atom, int dim, int inum);

template <typename T> __global__ void gpuKernelFixDragForce3D(T *x, T *v, T *f, T *eatom, T *vatom, T *fparam, T *box,
        int *iparam, int *ilist, int *pbc, int triclinic, int eflag_atom, int vflag_atom, int dim, int inum)
{
    int ii = threadIdx.x + blockIdx.x * blockDim.x;
    while (ii < inum) {
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
      gpuMinimumImage(dp, box, pbc, triclinic, dim);
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
      ii += blockDim.x * gridDim.x;
   }
}
template <typename T> __global__ void gpuKernelFixDragForce2D(T *x, T *v, T *f, T *eatom, T *vatom, T *fparam, T *box,
        int *iparam, int *ilist, int *pbc, int triclinic, int eflag_atom, int vflag_atom, int dim, int inum)
{
  int ii = threadIdx.x + blockIdx.x * blockDim.x;
  while (ii < inum) {
      int i = ilist[ii]; 
      T f_mag = fparam[0];
      T delta = fparam[1];          
      T dp[3]; 
      dp[0] = x[i*dim+0] - fparam[2];
      dp[1] = x[i*dim+1] - fparam[3];
      if (!iparam[0]) dp[0] = 0.0;
      if (!iparam[1]) dp[1] = 0.0;
      gpuMinimumImage(dp, box, pbc, triclinic, dim);
      T r = sqrt(dp[0]*dp[0] + dp[1]*dp[1]);
      if (r > delta) {
        T prefactor = f_mag/r;
        T fx = prefactor*dp[0];
        T fy = prefactor*dp[1];
        f[i*dim+0] -= fx;
        f[i*dim+1] -= fy;
      }
      ii += blockDim.x * gridDim.x;
  }        
}
template <typename T> void gpuFixDragForce(T *x, T *v, T *f, T *eatom, T *vatom, T *fparam, T *box,
        int *iparam, int *ilist, int *pbc, int triclinic, int eflag_atom, int vflag_atom, int dim, int inum)
{        
    int blockDim = 256;
    int gridDim = (inum + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;
    
    if (dim==3)
        gpuKernelFixDragForce3D<<<gridDim, blockDim>>>(x, v, f, eatom, vatom, fparam, box,
                iparam, ilist, pbc, triclinic, eflag_atom, vflag_atom, dim, inum);
    else
        gpuKernelFixDragForce2D<<<gridDim, blockDim>>>(x, v, f, eatom, vatom, fparam, box,
                iparam, ilist, pbc, triclinic, eflag_atom, vflag_atom, dim, inum);
}
template void gpuFixDragForce(double *x, double *v, double *f, double *eatom, double *vatom, 
        double *fparam, double *box, int *iparam, int *ilist, int *pbc, int triclinic, 
        int eflag_atom, int vflag_atom, int dim, int inum);
template void gpuFixDragForce(float *x, float *v, float *f, float *eatom, float *vatom, 
        float *fparam, float *box, int *iparam, int *ilist, int *pbc, int triclinic, int eflag_atom, 
        int vflag_atom, int dim, int inum);
        
template <typename T> __global__ void gpuKernelFixWallReflect(T *x, T *v, T *f, T *eatom, T *vatom, T *fparam, 
        int *iparam, int *ilist, int eflag_atom, int vflag_atom, int dim, int inum)
{  
  int ii = threadIdx.x + blockIdx.x * blockDim.x;
  while (ii < inum) {
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
      ii += blockDim.x * gridDim.x;
  }  
}
template <typename T> void gpuFixWallReflect(T *x, T *v, T *f, T *eatom, T *vatom, T *fparam, 
        int *iparam, int *ilist, int eflag_atom, int vflag_atom, int dim, int inum)        
{        
    int blockDim = 256;
    int gridDim = (inum + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;
    gpuKernelFixWallReflect<<<gridDim, blockDim>>>(x, v, f, eatom, vatom, fparam, 
            iparam, ilist, eflag_atom, vflag_atom, dim, inum);
}
template void gpuFixWallReflect(double *x, double *v, double *f, double *eatom, double *vatom, 
        double *fparam, int *iparam, int *ilist, int eflag_atom, int vflag_atom, int dim, int inum);
template void gpuFixWallReflect(float *x, float *v, float *f, float *eatom, float *vatom, 
        float *fparam, int *iparam, int *ilist, int eflag_atom, int vflag_atom, int dim, int inum);
 
template <typename T> __global__ void gpuKernelFixWallHarmonic(T *x, T *v, T *f, T *eatom, T *vatom, T *fparam, 
        int *iparam, int *ilist, int eflag_atom, int vflag_atom, int dim, int inum)
{
  int ii = threadIdx.x + blockIdx.x * blockDim.x;
  while (ii < inum) {
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
     ii += blockDim.x * gridDim.x;
  }  
}
template <typename T> void gpuFixWallHarmonic(T *x, T *v, T *f, T *eatom, T *vatom, T *fparam, 
        int *iparam, int *ilist, int eflag_atom, int vflag_atom, int dim, int inum)        
{        
    int blockDim = 256;
    int gridDim = (inum + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;
    gpuKernelFixWallHarmonic<<<gridDim, blockDim>>>(x, v, f, eatom, vatom, fparam, 
            iparam, ilist, eflag_atom, vflag_atom, dim, inum);
}
template void gpuFixWallHarmonic(double *x, double *v, double *f, double *eatom, double *vatom, 
        double *fparam, int *iparam, int *ilist, int eflag_atom, int vflag_atom, int dim, int inum);
template void gpuFixWallHarmonic(float *x, float *v, float *f, float *eatom, float *vatom, 
        float *fparam, int *iparam, int *ilist, int eflag_atom, int vflag_atom, int dim, int inum);

template <typename T> __global__ void gpuKernelFixWallLJ93(T *x, T *v, T *f, T *eatom, T *vatom, T *fparam, 
        int *iparam, int *ilist, int eflag_atom, int vflag_atom, int dim, int inum)
{
  int ii = threadIdx.x + blockIdx.x * blockDim.x;
  while (ii < inum) {
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
      ii += blockDim.x * gridDim.x;
  }  
}
template <typename T> void gpuFixWallLJ93(T *x, T *v, T *f, T *eatom, T *vatom, T *fparam, 
        int *iparam, int *ilist, int eflag_atom, int vflag_atom, int dim, int inum)        
{        
    int blockDim = 256;
    int gridDim = (inum + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;
    gpuKernelFixWallLJ93<<<gridDim, blockDim>>>(x, v, f, eatom, vatom, fparam, 
            iparam, ilist, eflag_atom, vflag_atom, dim, inum);
}
template void gpuFixWallLJ93(double *x, double *v, double *f, double *eatom, double *vatom, 
        double *fparam, int *iparam, int *ilist, int eflag_atom, int vflag_atom, int dim, int inum);
template void gpuFixWallLJ93(float *x, float *v, float *f, float *eatom, float *vatom, 
        float *fparam, int *iparam, int *ilist, int eflag_atom, int vflag_atom, int dim, int inum);

template <typename T> __global__ void gpuKernelFixWallLJ126(T *x, T *v, T *f, T *eatom, T *vatom, T *fparam, 
        int *iparam, int *ilist, int eflag_atom, int vflag_atom, int dim, int inum)
{
  int ii = threadIdx.x + blockIdx.x * blockDim.x;
  while (ii < inum) {
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
      ii += blockDim.x * gridDim.x;
  }    
}
template <typename T> void gpuFixWallLJ126(T *x, T *v, T *f, T *eatom, T *vatom, T *fparam, 
        int *iparam, int *ilist, int eflag_atom, int vflag_atom, int dim, int inum)        
{        
    int blockDim = 256;
    int gridDim = (inum + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;
    gpuKernelFixWallLJ126<<<gridDim, blockDim>>>(x, v, f, eatom, vatom, fparam, 
            iparam, ilist, eflag_atom, vflag_atom, dim, inum);
}
template void gpuFixWallLJ126(double *x, double *v, double *f, double *eatom, double *vatom, 
        double *fparam, int *iparam, int *ilist, int eflag_atom, int vflag_atom, int dim, int inum);
template void gpuFixWallLJ126(float *x, float *v, float *f, float *eatom, float *vatom, 
        float *fparam, int *iparam, int *ilist, int eflag_atom, int vflag_atom, int dim, int inum);

template <typename T> __global__ void gpuKernelFixWallLJ1043(T *x, T *v, T *f, T *eatom, T *vatom, T *fparam, 
        int *iparam, int *ilist, int eflag_atom, int vflag_atom, int dim, int inum)
{
  int ii = threadIdx.x + blockIdx.x * blockDim.x;
  while (ii < inum) {
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
      ii += blockDim.x * gridDim.x;
  }  
}
template <typename T> void gpuFixWallLJ1043(T *x, T *v, T *f, T *eatom, T *vatom, T *fparam, 
        int *iparam, int *ilist, int eflag_atom, int vflag_atom, int dim, int inum)        
{        
    int blockDim = 256;
    int gridDim = (inum + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;
    gpuKernelFixWallLJ1043<<<gridDim, blockDim>>>(x, v, f, eatom, vatom, fparam, 
            iparam, ilist, eflag_atom, vflag_atom, dim, inum);
}
template void gpuFixWallLJ1043(double *x, double *v, double *f, double *eatom, double *vatom, 
        double *fparam, int *iparam, int *ilist, int eflag_atom, int vflag_atom, int dim, int inum);
template void gpuFixWallLJ1043(float *x, float *v, float *f, float *eatom, float *vatom, 
        float *fparam, int *iparam, int *ilist, int eflag_atom, int vflag_atom, int dim, int inum);

template <typename T> __global__ void gpuKernelFixWallMorse(T *x, T *v, T *f, T *eatom, T *vatom, T *fparam, 
        int *iparam, int *ilist, int eflag_atom, int vflag_atom, int dim, int inum)
{
  int ii = threadIdx.x + blockIdx.x * blockDim.x;
  while (ii < inum) {
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
      ii += blockDim.x * gridDim.x;
  }  
}
template <typename T> void gpuFixWallMorse(T *x, T *v, T *f, T *eatom, T *vatom, T *fparam, 
        int *iparam, int *ilist, int eflag_atom, int vflag_atom, int dim, int inum)        
{        
    int blockDim = 256;
    int gridDim = (inum + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;
    gpuKernelFixWallMorse<<<gridDim, blockDim>>>(x, v, f, eatom, vatom, fparam, 
            iparam, ilist, eflag_atom, vflag_atom, dim, inum);
}
template void gpuFixWallMorse(double *x, double *v, double *f, double *eatom, double *vatom, 
        double *fparam, int *iparam, int *ilist, int eflag_atom, int vflag_atom, int dim, int inum);
template void gpuFixWallMorse(float *x, float *v, float *f, float *eatom, float *vatom, 
        float *fparam, int *iparam, int *ilist, int eflag_atom, int vflag_atom, int dim, int inum);

#endif

