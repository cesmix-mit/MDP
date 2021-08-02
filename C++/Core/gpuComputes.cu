/***************************************************************************                               
                    Molecular Dynamics Potentials (MDP)
                           CESMIX-MIT Project  
 
 Contributing authors: Ngoc-Cuong Nguyen (cuongng@mit.edu, exapde@gmail.com)
 ***************************************************************************/

#ifndef MDP_GPUCOMPUTES
#define MDP_GPUCOMPUTES

#include <math.h>
#include <algorithm>

using std::min;
using std::max;

#define SWAP(a,b) do {       \
    tmp = a; a = b; b = tmp; \
  } while(0)

#define ISWAP(a,b) do {        \
    itmp = a; a = b; b = itmp; \
  } while(0)

// template <typename T> __device__ void gpuMinimumImage(T *dp, T *box, int *pbc, int triclinic, int dim)
// {
//   if (dim==3) {
//       T xprd = box[0];
//       T yprd = box[1];
//       T zprd = box[2];  
// 
//       if (triclinic == 0) {
//         if (pbc[0]) {
//           while (fabs(dp[0]) > 0.5*xprd) {
//             if (dp[0] < 0.0) dp[0] += xprd;
//             else dp[0] -= xprd;
//           }
//         }
//         if (pbc[1]) {
//           while (fabs(dp[1]) > 0.5*yprd) {
//             if (dp[1] < 0.0) dp[1] += yprd;
//             else dp[1] -= yprd;
//           }
//         }
//         if (pbc[2]) {
//           while (fabs(dp[2]) > 0.5*zprd) {
//             if (dp[2] < 0.0) dp[2] += zprd;
//             else dp[2] -= zprd;
//           }
//         }
// 
//       } else {
//         T xy = box[3];
//         T xz = box[4];
//         T yz = box[5];        
//         if (pbc[2]) {
//           while (fabs(dp[2]) > 0.5*zprd) {
//             if (dp[2] < 0.0) {
//               dp[2] += zprd;
//               dp[1] += yz;
//               dp[0] += xz;
//             } else {
//               dp[2] -= zprd;
//               dp[1] -= yz;
//               dp[0] -= xz;
//             }
//           }
//         }
//         if (pbc[1]) {
//           while (fabs(dp[1]) > 0.5*yprd) {
//             if (dp[1] < 0.0) {
//               dp[1] += yprd;
//               dp[0] += xy;
//             } else {
//               dp[1] -= yprd;
//               dp[0] -= xy;
//             }
//           }
//         }
//         if (pbc[0]) {
//           while (fabs(dp[0]) > 0.5*xprd) {
//             if (dp[0] < 0.0) dp[0] += xprd;
//             else dp[0] -= xprd;
//           }
//         }
//       }
//   } else {
//       T xprd = box[0];
//       T yprd = box[1];
// 
//       if (triclinic == 0) {
//         if (pbc[0]) {
//           while (fabs(dp[0]) > 0.5*xprd) {
//             if (dp[0] < 0.0) dp[0] += xprd;
//             else dp[0] -= xprd;
//           }
//         }
//         if (pbc[1]) {
//           while (fabs(dp[1]) > 0.5*yprd) {
//             if (dp[1] < 0.0) dp[1] += yprd;
//             else dp[1] -= yprd;
//           }
//         }    
//       } else {
//         T xy = box[3];
//         T xz = box[4];
//         if (pbc[1]) {
//           while (fabs(dp[1]) > 0.5*yprd) {
//             if (dp[1] < 0.0) {
//               dp[1] += yprd;
//               dp[0] += xy;
//             } else {
//               dp[1] -= yprd;
//               dp[0] -= xy;
//             }
//           }
//         }
//         if (pbc[0]) {
//           while (fabs(dp[0]) > 0.5*xprd) {
//             if (dp[0] < 0.0) dp[0] += xprd;
//             else dp[0] -= xprd;
//           }
//         }
//       }      
//   }
// }
        
template <typename T> __global__ void gpuKernelPackIntProperty(T *buf, int *prop, int *ilist, 
         int m, int mvalues, int n, int nvalues, int inum)
{
  int ii = threadIdx.x + blockIdx.x * blockDim.x;
  while (ii < inum) {
    int i = ilist[ii];
    buf[n + nvalues*i] = (T) prop[m + mvalues*i];
    ii += blockDim.x * gridDim.x;
  }
}
template <typename T> void gpuPackIntProperty(T *buf, int *prop, int *ilist, 
         int m, int mvalues, int n, int nvalues, int inum)
{        
    int blockDim = 256;
    int gridDim = (inum + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;
    gpuKernelPackIntProperty<<<gridDim, blockDim>>>(buf, prop, ilist, m, mvalues, n, nvalues, inum);
}
template void gpuPackIntProperty(double *buf, int *prop, int *ilist, int m, int mvalues, int n, int nvalues, int inum);
template void gpuPackIntProperty(float *buf, int *prop, int *ilist, int m, int mvalues, int n, int nvalues, int inum);

template <typename T> __global__ void gpuKernelPackIntProperty(T *buf, int *prop, int *type, int *ilist, 
         int n, int nvalues, int inum)
{
  int ii = threadIdx.x + blockIdx.x * blockDim.x;  
  while (ii < inum) {
    int i = ilist[ii]; 
    buf[n + nvalues*i] = (T) prop[type[i]];
    ii += blockDim.x * gridDim.x;
  }
}
template <typename T> void gpuPackIntProperty(T *buf, int *prop, int *type, int *ilist, 
         int n, int nvalues, int inum)
{        
    int blockDim = 256;
    int gridDim = (inum + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;
    gpuKernelPackIntProperty<<<gridDim, blockDim>>>(buf, prop, type, ilist, n, nvalues, inum);
}
template void gpuPackIntProperty(double *buf, int *prop, int *type, int *ilist, int n, int nvalues, int inum);
template void gpuPackIntProperty(float *buf, int *prop, int *type, int *ilist, int n, int nvalues, int inum);

template <typename T> __global__ void gpuKernelPackFloatProperty(T *buf, T *prop, int *ilist, 
         int m, int mvalues, int n, int nvalues, int inum)
{
  int ii = threadIdx.x + blockIdx.x * blockDim.x;  
  while (ii < inum) {
    int i = ilist[ii]; 
    buf[n + nvalues*i] = prop[m + mvalues*i];
    ii += blockDim.x * gridDim.x;
  }
}
template <typename T> void gpuPackFloatProperty(T *buf, T *prop, int *ilist, 
         int m, int mvalues, int n, int nvalues, int inum)
{        
    int blockDim = 256;
    int gridDim = (inum + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;
    gpuKernelPackFloatProperty<<<gridDim, blockDim>>>(buf, prop, ilist, m, mvalues, n, nvalues, inum);
}
template void gpuPackFloatProperty(double *buf, double *prop, int *ilist, int m, int mvalues, int n, int nvalues, int inum);
template void gpuPackFloatProperty(float *buf, float *prop, int *ilist, int m, int mvalues, int n, int nvalues, int inum);

template <typename T> __global__ void gpuKernelPackFloatProperty(T *buf, T *prop, T a, T b, int *ilist, 
         int m, int mvalues, int n, int nvalues, int inum)
{
  int ii = threadIdx.x + blockIdx.x * blockDim.x;  
  while (ii < inum) {
    int i = ilist[ii]; 
    buf[n + nvalues*i] = a*prop[m + mvalues*i] + b;
    ii += blockDim.x * gridDim.x;
  }
}
template <typename T> void gpuPackFloatProperty(T *buf, T *prop, T a, T b, int *ilist, 
         int m, int mvalues, int n, int nvalues, int inum)
{        
    int blockDim = 256;
    int gridDim = (inum + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;
    gpuKernelPackFloatProperty<<<gridDim, blockDim>>>(buf, prop, a, b, ilist, m, mvalues, n, nvalues, inum);
}
template void gpuPackFloatProperty(double *buf, double *prop, double a, double b, int *ilist, 
         int m, int mvalues, int n, int nvalues, int inum);
template void gpuPackFloatProperty(float *buf, float *prop, float a, float b, int *ilist, 
         int m, int mvalues, int n, int nvalues, int inum);

template <typename T> __global__ void gpuKernelPackFloatProperty(T *buf, T *prop, int *type, int *ilist, 
         int n, int nvalues, int inum)
{
  int ii = threadIdx.x + blockIdx.x * blockDim.x;   
  while (ii < inum) {
    int i = ilist[ii]; 
    buf[n + nvalues*i] = prop[type[i]];
    ii += blockDim.x * gridDim.x;
  }
}
template <typename T> void gpuPackFloatProperty(T *buf, T *prop, int *type, int *ilist, 
         int n, int nvalues, int inum)
{        
    int blockDim = 256;
    int gridDim = (inum + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;
    gpuKernelPackFloatProperty<<<gridDim, blockDim>>>(buf, prop, type, ilist, n, nvalues, inum);
}
template void gpuPackFloatProperty(double *buf, double *prop, int *type, int *ilist, 
         int n, int nvalues, int inum);
template void gpuPackFloatProperty(float *buf, float *prop, int *type, int *ilist, 
         int n, int nvalues, int inum);

template <typename T> T gpuComputeMass(T *amass, T *mass, T *tmp, int *type, int *ilist, int inum)
{
    T masstotal;
    gpuPackFloatProperty(amass, mass, type, ilist, 0, 1, inum);
    masstotal = gpuArraySum(amass, tmp, inum);
    return masstotal;
}         
template double gpuComputeMass(double *amass, double *mass, double *tmp, int *type, int *ilist, int inum);
template float gpuComputeMass(float *amass, float *mass, float *tmp, int *type, int *ilist, int inum);

template <typename T> __global__ void gpuKernelComputeXCM(T *axcm, T *x, T *mass, T *box, 
        int *ilist, int *type, int *image, int triclinic, int dim, int inum)
{
  int ii = threadIdx.x + blockIdx.x * blockDim.x;
  while (ii < inum) {
        int i = ilist[ii]; 
        T massone = mass[type[i]];
        T y0, y1, y2;  
        if (triclinic == 0) {
            y0 = x[0+i*dim] + image[i*dim+0]*box[0];
            y1 = x[1+i*dim] + image[i*dim+1]*box[1];
            y2 = (dim==3) ? (x[2+i*dim] + image[i*dim+2]*box[2]) : 0.0;
        } else {
            y0 = x[0+i*dim] + box[0]*image[i*dim+0] + box[5]*image[i*dim+1] + box[4]*image[i*dim+2];
            y1 = x[1+i*dim] + box[1]*image[i*dim+1] + box[3]*image[i*dim+2];            
            y2 = (dim==3) ? (x[2+i*dim] + image[i*dim+2]*box[2]) : 0.0;
        }                           
        axcm[ii] = y0 * massone;
        axcm[ii+inum] = y1 * massone;
        if (dim==3) axcm[ii+2*inum] = y2 * massone;      
        ii += blockDim.x * gridDim.x;
  }
}
template <typename T> void gpuComputeXCM(T *xcm, T *axcm, T *x, T *tmp, T *mass, T *box, 
        T masstotal, int *ilist, int *type, int *image, int triclinic, int dim, int inum)
{        
    int blockDim = 256;
    int gridDim = (inum + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;
    gpuKernelComputeXCM<<<gridDim, blockDim>>>(axcm, x, mass, box, ilist, type, image, triclinic, dim, inum);

    if (masstotal > 0.0) 
        for (int i=0; i<dim; i++)
            xcm[i] = gpuArraySum(&axcm[inum*i], tmp, inum)/masstotal;    
}
template void  gpuComputeXCM(double *xcm, double *axcm, double *x, double *tmp, double *mass, double *box, 
        double masstotal, int *ilist, int *type, int *image, int triclinic, int dim, int inum);
template void  gpuComputeXCM(float *xcm, float *axcm, float *x, float *tmp, float *mass, float *box, 
        float masstotal, int *ilist, int *type, int *image, int triclinic, int dim, int inum);

template <typename T> __global__ void gpuKernelComputeVCM(T *avcm, T *v, T *mass, 
        int *ilist, int *type, int dim, int inum)
{  
    int ii = threadIdx.x + blockIdx.x * blockDim.x;
    while (ii < inum) {
        int i = ilist[ii]; 
        T massone = mass[type[i]];
        avcm[ii] = v[i*dim+0]*massone;
        avcm[ii+inum] = v[i*dim+1]*massone;
        if (dim==3) avcm[ii+2*inum] = v[i*dim+2]*massone;
        ii += blockDim.x * gridDim.x;
    }  
}
template <typename T> void gpuComputeVCM(T *vcm, T *avcm, T *v, T *tmp, T *mass,
        T masstotal, int *ilist, int *type, int dim, int inum)
{        
    int blockDim = 256;
    int gridDim = (inum + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;
    gpuKernelComputeVCM<<<gridDim, blockDim>>>(avcm, v, mass, ilist, type, dim, inum);

    if (masstotal > 0.0) 
        for (int i=0; i<dim; i++)
            vcm[i] = gpuArraySum(&avcm[inum*i], tmp, inum)/masstotal;    
}
template void  gpuComputeVCM(double *vcm, double *avcm, double *v, double *tmp, double *mass, double masstotal, 
        int *ilist, int *type, int dim, int inum);
template void  gpuComputeVCM(float *vcm, float *avcm, float *v, float *tmp, float *mass, float masstotal, 
        int *ilist, int *type, int dim, int inum);

template <typename T> __global__ void gpuKernelComputeGyration(T * ag, T *xcm, T *x, T *mass, T *box, 
        int *ilist, int *type, int *image, int triclinic, int dim, int inum)
{
  int ii = threadIdx.x + blockIdx.x * blockDim.x;
  while (ii < inum) {
        int i = ilist[ii]; 
        T massone = mass[type[i]];
        T y0, y1, y2;  
        if (triclinic == 0) {
            y0 = x[0+i*dim] + image[i*dim+0]*box[0];
            y1 = x[1+i*dim] + image[i*dim+1]*box[1];
            y2 = (dim==3) ? (x[2+i*dim] + image[i*dim+2]*box[2]) : 0.0;
        } else {
            y0 = x[0+i*dim] + box[0]*image[i*dim+0] + box[5]*image[i*dim+1] + box[4]*image[i*dim+2];
            y1 = x[1+i*dim] + box[1]*image[i*dim+1] + box[3]*image[i*dim+2];            
            y2 = (dim==3) ? (x[2+i*dim] + image[i*dim+2]*box[2]) : 0.0;
        }                           
        y0 = y0 - xcm[0];
        y1 = y1 - xcm[1];        
        y2 = (dim==3) ? (y2 - xcm[2]) : 0.0;
        ag[ii] = (y0*y0 + y1*y1 + y2*y2) * massone;
        ii += blockDim.x * gridDim.x;
  }    
}
template <typename T> T gpuComputeGyration(T * ag, T *xcm, T *x, T *tmp, T *mass, T *box, T masstotal,
        int *ilist, int *type, int *image, int triclinic, int dim, int inum)
{        
    int blockDim = 256;
    int gridDim = (inum + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;
    gpuKernelComputeGyration<<<gridDim, blockDim>>>(ag, xcm, x, mass, box, ilist, 
            type, image, triclinic, dim, inum);

    if (masstotal > 0.0)         
        return sqrt(gpuArraySum(ag, tmp, inum)/masstotal);    
    
    return 0.0;        
}
template double gpuComputeGyration(double *ag, double *xcm, double *x, double *tmp, double *mass, double *box, 
        double masstotal, int *ilist, int *type, int *image, int triclinic, int dim, int inum);
template float  gpuComputeGyration(float *ag, float *xcm, float *x, float *tmp, float *mass, float *box, 
        float masstotal, int *ilist, int *type, int *image, int triclinic, int dim, int inum);

template <typename T> __global__ void gpuKernelComputeAngmom(T *p, T *xcm, T *x, T *v, T *mass, T *box, 
        int *ilist, int *type, int *image, int triclinic, int dim, int inum)
{
  int ii = threadIdx.x + blockIdx.x * blockDim.x;          
  while (ii < inum) {
        int i = ilist[ii]; 
        T massone = mass[type[i]];
        T y0, y1, y2;  
        if (triclinic == 0) {
            y0 = x[0+i*dim] + image[i*dim+0]*box[0];
            y1 = x[1+i*dim] + image[i*dim+1]*box[1];
            y2 = (dim==3) ? (x[2+i*dim] + image[i*dim+2]*box[2]) : 0.0;
        } else {
            y0 = x[0+i*dim] + box[0]*image[i*dim+0] + box[5]*image[i*dim+1] + box[4]*image[i*dim+2];
            y1 = x[1+i*dim] + box[1]*image[i*dim+1] + box[3]*image[i*dim+2];            
            y2 = (dim==3) ? (x[2+i*dim] + image[i*dim+2]*box[2]) : 0.0;
        }                           
        y0 = y0 - xcm[0];
        y1 = y1 - xcm[1];        
        p[ii+2*inum] = massone * (y0*v[i*dim+1] - y1*v[i*dim+0]);
        if (dim==3) {
            y2 = y2 - xcm[2];
            p[ii] = massone * (y1*v[i*dim+2] - y2*v[i*dim+1]);
            p[ii+inum] = massone * (y2*v[i*dim+0] - y0*v[i*dim+2]);
        }
        ii += blockDim.x * gridDim.x;
  }   
}
template <typename T> void gpuComputeAngmom(T *lmom, T *p, T *xcm, T *x, T *v, T *tmp, T *mass, T *box, 
        int *ilist, int *type, int *image, int triclinic, int dim, int inum)
{        
    int blockDim = 256;
    int gridDim = (inum + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;
    gpuKernelComputeAngmom<<<gridDim, blockDim>>>(p, xcm, x, v, mass, box, ilist, 
            type, image, triclinic, dim, inum);
    
    for (int i=0; i<dim; i++)
        lmom[i] = gpuArraySum(&p[inum*i], tmp, inum);        
}
template void gpuComputeAngmom(double *lmom, double *p, double *xcm, double *x, double *v, double *tmp, 
        double *mass, double *box, int *ilist, int *type, int *image, int triclinic, int dim, int inum);
template void gpuComputeAngmom(float *lmom, float *p, float *xcm, float *x, float *v, float *tmp, 
        float *mass, float *box, int *ilist, int *type, int *image, int triclinic, int dim, int inum);

template <typename T> __global__ void gpuKernelComputeTorque(T *tlocal, T *xcm, T *x, T *f, T *box, 
        int *ilist, int *image, int triclinic, int dim, int inum)
{
  T unwrap[3];
  int ii = threadIdx.x + blockIdx.x * blockDim.x;          
  while (ii < inum) {
        int i = ilist[ii]; 
        T y0, y1, y2;  
        if (triclinic == 0) {
            y0 = x[0+i*dim] + image[i*dim+0]*box[0];
            y1 = x[1+i*dim] + image[i*dim+1]*box[1];
            y2 = (dim==3) ? (x[2+i*dim] + image[i*dim+2]*box[2]) : 0.0;
        } else {
            y0 = x[0+i*dim] + box[0]*image[i*dim+0] + box[5]*image[i*dim+1] + box[4]*image[i*dim+2];
            y1 = x[1+i*dim] + box[1]*image[i*dim+1] + box[3]*image[i*dim+2];            
            y2 = (dim==3) ? (x[2+i*dim] + image[i*dim+2]*box[2]) : 0.0;
        }                           
        y0 = y0 - xcm[0];
        y1 = y1 - xcm[1];        
        tlocal[ii+2*inum] = y0*f[i*dim+1] - y1*f[i*dim+0];
        if (dim==3) {
            y2 = y2 - xcm[2];
            tlocal[ii] = y1*f[i*dim+2] - y2*f[i*dim+1];
            tlocal[ii+inum] = y2*f[i*dim+0] - y0*f[i*dim+2];
        }
        ii += blockDim.x * gridDim.x;
  }     
}
template <typename T> void gpuComputeTorque(T *tq, T *q, T *xcm, T *x, T *f, T *tmp, T *box, 
        int *ilist, int *image, int triclinic, int dim, int inum)
{        
    int blockDim = 256;
    int gridDim = (inum + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;
    gpuKernelComputeTorque<<<gridDim, blockDim>>>(q, xcm, x, f, box, ilist, 
            image, triclinic, dim, inum);
            
    for (int i=0; i<3; i++)
        tq[i] = gpuArraySum(&q[inum*i], tmp, inum);        
}
template void gpuComputeTorque(double *tq, double *q, double *xcm, double *x, double *f, double *tmp, 
         double *box, int *ilist, int *image, int triclinic, int dim, int inum);
template void gpuComputeTorque(float *tq, float *q, float *xcm, float *x, float *f, float *tmp, 
        float *box, int *ilist, int *image, int triclinic, int dim, int inum);

template <typename T> __global__ void gpuKernelComputeInertia(T *ione, T *xcm, T *x, T *mass, T *box, 
        int *ilist, int *type, int *image, int triclinic, int dim, int inum)
{
  int ii = threadIdx.x + blockIdx.x * blockDim.x;        
  while (ii < inum) {
        int i = ilist[ii]; 
        T massone = mass[type[i]];
        T y0, y1, y2;  
        if (triclinic == 0) {
            y0 = x[0+i*dim] + image[i*dim+0]*box[0];
            y1 = x[1+i*dim] + image[i*dim+1]*box[1];
            y2 = (dim==3) ? (x[2+i*dim] + image[i*dim+2]*box[2]) : 0.0;
        } else {
            y0 = x[0+i*dim] + box[0]*image[i*dim+0] + box[5]*image[i*dim+1] + box[4]*image[i*dim+2];
            y1 = x[1+i*dim] + box[1]*image[i*dim+1] + box[3]*image[i*dim+2];            
            y2 = (dim==3) ? (x[2+i*dim] + image[i*dim+2]*box[2]) : 0.0;
        }                           
        y0 = y0 - xcm[0];
        y1 = y1 - xcm[1];                        
        ione[ii] = massone * (y1*y1);
        ione[ii+inum] = massone * (y0*y0);
        ione[ii+2*inum] = massone * (y0*y0 + y1*y1);
        ione[ii+3*inum] = -massone * y0*y1;
        if (dim==3) {
            y2 = y2 - xcm[2];
            ione[ii] += massone * (y2*y2);
            ione[ii+inum] += massone * (y2*y2);
            ione[ii+4*inum] = -massone * y0*y2;
            ione[ii+5*inum] = -massone * y1*y2;            
        }
  }    
  //MPI_Allreduce(&ione[0*dim+0],&itensor[0*dim+0],9,MPI_DOUBLE,MPI_SUM,world);
}
template <typename T> void gpuComputeInertia(T *inertia, T *ione, T *xcm, T *x, T *tmp, T *mass, T *box, 
        int *ilist, int *type, int *image, int triclinic, int dim, int inum)
{        
    int blockDim = 256;
    int gridDim = (inum + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;
    gpuKernelComputeInertia<<<gridDim, blockDim>>>(ione, xcm, x, mass, box, ilist, 
            type, image, triclinic, dim, inum);
    
    if (dim==2) {
        for (int i=0; i<4; i++)
            inertia[i] = gpuArraySum(&ione[inum*i], tmp, inum);        
    } else {
        for (int i=0; i<6; i++)
            inertia[i] = gpuArraySum(&ione[inum*i], tmp, inum);        
    }
}
template void gpuComputeInertia(double *inertia, double *ione, double *xcm, double *x, double *tmp, 
        double *mass, double *box, int *ilist, int *type, int *image, int triclinic, int dim, int inum);
template void gpuComputeInertia(float *inertia, float *ione, float *xcm, float *x, float *tmp,
        float *mass, float *box, int *ilist, int *type, int *image, int triclinic, int dim, int inum);

template <typename T> __global__ void gpuKernelComputeKEAtom(T *ke, T *mass,  T *v, 
        T mvv2e, int *type, int *ilist,  int dim, int inum)
{    
    int ii = threadIdx.x + blockIdx.x * blockDim.x;        
    while (ii < inum) {
        int i = ilist[ii]; 
        ke[ii] = v[i*dim+0]*v[i*dim+0] + v[i*dim+1]*v[i*dim+1];        
        if (dim==3)
            ke[ii] += v[i*dim+2]*v[i*dim+2];
        ke[ii] *=  0.5 * mvv2e * mass[type[i]];
        ii += blockDim.x * gridDim.x;
    }    
}
template <typename T> void gpuComputeKEAtom(T *ke, T *mass, T *v, 
        T mvv2e, int *type, int *ilist, int dim, int inum)
{        
    int blockDim = 256;
    int gridDim = (inum + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;
    gpuKernelComputeKEAtom<<<gridDim, blockDim>>>(ke, mass, v, mvv2e, 
            type, ilist, dim, inum);
    
}
template void gpuComputeKEAtom(double *ke, double *mass,  double *v, 
        double mvv2e, int *type, int *ilist,  int dim, int inum);
template void gpuComputeKEAtom(float *ke, float *mass,  float *v, 
        float mvv2e, int *type, int *ilist,  int dim, int inum);

template <typename T> T gpuComputeTempScalar(T *ke, T *v, T *tmp, T *mass, T tfactor, 
        int *type, int *ilist, int dim, int inum)
{
    gpuComputeKEAtom(ke, mass, v, (T) 2.0*tfactor, type, ilist, dim, inum);

    return gpuArraySum(ke, tmp, inum);    
}
template double gpuComputeTempScalar(double *ke, double *v, double *tmp, double *mass, double tfactor, 
        int *type, int *ilist, int dim, int inum);
template float gpuComputeTempScalar(float *ke, float *v, float *tmp, float *mass, float tfactor, 
        int *type, int *ilist, int dim, int inum);

template <typename T> __global__ void gpuKernelComputeKETensorAtom(T *stress, T *mass, T *v, 
        T mvv2e, int *type, int *ilist, int dim, int inum)
{  
  // mvv2e converts mv^2 to energy
    int ii = threadIdx.x + blockIdx.x * blockDim.x;  
    while (ii < inum) {    
        int i = ilist[ii];   
        T onemass = mvv2e * mass[type[i]];
        stress[ii] += onemass*v[i*dim+0]*v[i*dim+0];
        stress[ii+inum] += onemass*v[i*dim+1]*v[i*dim+1];
        stress[ii+3*inum] += onemass*v[i*dim+0]*v[i*dim+1];
        if (dim==3) {
        stress[ii+2*inum] += onemass*v[i*dim+2]*v[i*dim+2];
        stress[ii+4*inum] += onemass*v[i*dim+0]*v[i*dim+2];
        stress[ii+5*inum] += onemass*v[i*dim+1]*v[i*dim+2];
        }
        ii += blockDim.x * gridDim.x;
    }      
}

template <typename T> __global__ void gpuKernelComputeCentroidTensorAtom(T *stress, T *mass, T *v, 
        T mvv2e, int *type, int *ilist, int dim, int inum)
{  
  // mvv2e converts mv^2 to energy
    int ii = threadIdx.x + blockIdx.x * blockDim.x;  
    while (ii < inum) {    
        int i = ilist[ii];   
        T onemass = mvv2e * mass[type[i]];
        stress[ii] += onemass*v[i*dim+0]*v[i*dim+0];
        stress[ii+inum] += onemass*v[i*dim+1]*v[i*dim+1];
        stress[ii+3*inum] += onemass*v[i*dim+0]*v[i*dim+1];
        stress[ii+6*inum] += onemass*v[i*dim+1]*v[i*dim+0];
        if (dim==3) {
        stress[ii+2*inum] += onemass*v[i*dim+2]*v[i*dim+2];
        stress[ii+4*inum] += onemass*v[i*dim+0]*v[i*dim+2];
        stress[ii+5*inum] += onemass*v[i*dim+1]*v[i*dim+2];
        stress[ii+7*inum] += onemass*v[i*dim+2]*v[i*dim+0];
        stress[ii+8*inum] += onemass*v[i*dim+2]*v[i*dim+1];
        }
        ii += blockDim.x * gridDim.x;
    }      
}

template <typename T> void gpuComputeStressAtom(T *stress, T *mass, T *vatom, T *v, T mvv2e, 
        T nktv2p, int *type, int *ilist,  int vflag, int keflag, int dim, int inum)
{
    if (vflag)
        gpuArrayTransposeAtColumnIndex(stress, vatom, ilist, 6, inum);
    else
        gpuArraySetValue(stress, (T) 0.0, inum*6);

    if (keflag) {
        int blockDim = 256;
        int gridDim = (inum + blockDim - 1) / blockDim;
        gridDim = (gridDim>1024)? 1024 : gridDim;
        gpuKernelComputeKETensorAtom<<<gridDim, blockDim>>>(stress, mass, v, mvv2e, 
                type, ilist, dim, inum);
    }

    gpuArrayMultiplyScalar(stress, nktv2p, inum*6);
}
template void gpuComputeStressAtom(double *stress, double *mass, double *vatom, double *v, double mvv2e, 
        double nktv2p, int *type, int *ilist, int vflag, int keflag, int dim, int inum);
template void gpuComputeStressAtom(float *stress, float *mass, float *vatom, float *v, float mvv2e, 
        float nktv2p, int *type, int *ilist, int vflag, int keflag, int dim, int inum);

template <typename T> void gpuComputeCentroidStressAtom(T *stress, T *mass, T *cvatom, T *v, T mvv2e, 
        T nktv2p, int *type, int *ilist,  int vflag, int keflag, int dim, int inum)
{
    if (vflag)        
        gpuArrayTransposeAtColumnIndex(stress, cvatom, ilist, 9, inum);
    else
        gpuArraySetValue(stress, (T) 0.0, inum*9);

    if (keflag) {
        int blockDim = 256;
        int gridDim = (inum + blockDim - 1) / blockDim;
        gridDim = (gridDim>1024)? 1024 : gridDim;
        gpuKernelComputeCentroidTensorAtom<<<gridDim, blockDim>>>(stress, mass, v, mvv2e, 
                type, ilist, dim, inum);
    }

    gpuArrayMultiplyScalar(stress, nktv2p, inum*9);
}
template void gpuComputeCentroidStressAtom(double *stress, double *mass, double *cvatom, double *v, double mvv2e, 
        double nktv2p, int *type, int *ilist, int vflag, int keflag, int dim, int inum);
template void gpuComputeCentroidStressAtom(float *stress, float *mass, float *cvatom, float *v, float mvv2e, 
        float nktv2p, int *type, int *ilist, int vflag, int keflag, int dim, int inum);

template <typename T> void gpuComputeTempSymTensor(T *ke_tensor, T *stress, T *v, T *tmp, T *mass, T tfactor, 
        int *type, int *ilist, int dim, int inum)
{
    for (int i = 0; i < 6; i++) ke_tensor[i] = 0.0;

    gpuArraySetValue(stress, (T) 0.0, inum*6);
    
    int blockDim = 256;
    int gridDim = (inum + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;
    gpuKernelComputeKETensorAtom<<<gridDim, blockDim>>>(stress, mass, v, tfactor, 
            type, ilist, dim, inum);

    if (dim==2) {
        for (int i=0; i<4; i++)
            ke_tensor[i] = gpuArraySum(&stress[inum*i], tmp, inum);        
    } else {
        for (int i=0; i<6; i++)
            ke_tensor[i] = gpuArraySum(&stress[inum*i], tmp, inum);        
    }        
}
template void gpuComputeTempSymTensor(double *ke_tensor, double *stress, double *v, double *tmp, double *mass, 
        double tfactor, int *type, int *ilist, int dim, int inum);
template void gpuComputeTempSymTensor(float *ke_tensor, float *stress, float *v, float *tmp, float *mass, 
        float tfactor, int *type, int *ilist, int dim, int inum);

template <typename T> T gpuComputePressureScalar(T *virial, T volume, T temp, T tempdof, 
        T boltz, T nktv2p, int dim)
{
  T factor = nktv2p /(dim*volume);
  if (dim == 3) {
      return (tempdof * boltz * temp + virial[0] + virial[1] + virial[2]) * factor;
  } else {
      return (tempdof * boltz * temp + virial[0] + virial[1]) *factor;
  }  
}
template double gpuComputePressureScalar(double *virial, double volume, double temp, double tempdof, 
        double boltz, double nktv2p, int dim);
template float gpuComputePressureScalar(float *virial, float volume, float temp, float tempdof, 
        float boltz, float nktv2p, int dim);

template <typename T> void gpuComputePressureSymTensor(T *vector, T *virial, T *ke_tensor, 
        T volume, T nktv2p, int dim)
{
  T factor = nktv2p /volume;
  
  if (dim == 3) {
      for (int i = 0; i < 6; i++)
        vector[i] = (ke_tensor[i] + virial[i]) * factor;
  } else {
      vector[0] = (ke_tensor[0] + virial[0]) * factor;
      vector[1] = (ke_tensor[1] + virial[1]) * factor;
      vector[3] = (ke_tensor[3] + virial[3]) * factor;
      vector[2] = vector[4] = vector[5] = 0.0;
  }
}
template void gpuComputePressureSymTensor(double *vector, double *virial, double *ke_tensor, 
        double volume, double nktv2p, int dim);
template void gpuComputePressureSymTensor(float *vector, float *virial, float *ke_tensor, 
        float volume, float nktv2p, int dim);

template <typename T> __global__ void gpuKernelComputeHeatCentroidAtom(T *jc, T *ke, T *pe, T *stress, T *v, 
        int *ilist, int dim, int inum)
{   
    int ii = threadIdx.x + blockIdx.x * blockDim.x;
    while (ii < inum) {
        int i = ilist[ii]; 
        T eng = pe[i] + ke[i];
        jc[ii+0*inum] = eng*v[i*dim+0];
        jc[ii+1*inum] = eng*v[i*dim+1];
        jc[ii+2*inum] = eng*v[i*dim+2];
        jc[ii+3*inum] = -(stress[ii+0*inum]*v[i*dim+0] + stress[ii+3*inum]*v[i*dim+1] +
          stress[ii+4*inum]*v[i*dim+2]);
        jc[ii+4*inum] = -(stress[ii+6*inum]*v[i*dim+0] + stress[ii+1*inum]*v[i*dim+1] +
          stress[ii+5*inum]*v[i*dim+2]);        
        jc[ii+5*inum] = -(stress[ii+7*inum]*v[i*dim+0] + stress[ii+8*inum]*v[i*dim+1] +
          stress[ii+2*inum]*v[i*dim+2]);      
        ii += blockDim.x * gridDim.x;
    }    
}

template <typename T> __global__ void gpuKernelComputeHeatTensorAtom(T *jc, T *ke, T *pe, T *stress, T *v, 
        int *ilist, int dim, int inum)
{   
    int ii = threadIdx.x + blockIdx.x * blockDim.x;
    while (ii < inum) {
        int i = ilist[ii]; 
        T eng = pe[i] + ke[i];
        jc[ii+0*inum] = eng*v[i*dim+0];
        jc[ii+1*inum] = eng*v[i*dim+1];
        jc[ii+2*inum] = eng*v[i*dim+2];
        jc[ii+3*inum] = -(stress[ii+0*inum]*v[i*dim+0] + stress[ii+3*inum]*v[i*dim+1] +
          stress[ii+4*inum]*v[i*dim+2]);
        jc[ii+4*inum] = -(stress[ii+3*inum]*v[i*dim+0] + stress[ii+1*inum]*v[i*dim+1] +
          stress[ii+5*inum]*v[i*dim+2]);
        jc[ii+5*inum] = -(stress[ii+4*inum]*v[i*dim+0] + stress[ii+5*inum]*v[i*dim+1] +
          stress[ii+2*inum]*v[i*dim+2]);      
        ii += blockDim.x * gridDim.x;
    }    
}

template <typename T> void gpuComputeHeatFlux(T *vector, T *jc, T *ke, T *pe, T *stress, T *v, 
        T *tmp, T nktv2p, int *ilist,  int pressatomflag, int dim, int inum)
{   
    int blockDim = 256;
    int gridDim = (inum + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;
    
    if (pressatomflag == 2) {
        gpuKernelComputeHeatCentroidAtom<<<gridDim, blockDim>>>(jc, ke, pe, stress, v, ilist, dim, inum);
    } else {
        gpuKernelComputeHeatTensorAtom<<<gridDim, blockDim>>>(jc, ke, pe, stress, v, ilist, dim, inum);
    }

    T jcsum[6];
    for (int i=0; i<6; i++)
        jcsum[i] = gpuArraySum(&jc[inum*i], tmp, inum);        
    
    jcsum[3] /= nktv2p;
    jcsum[4] /= nktv2p;
    jcsum[5] /= nktv2p;
    
    for (int i=0; i<3; i++) {
        vector[i] = jcsum[i]+jcsum[i+3]; // 1st 3 terms are total heat flux
        vector[i+3] = jcsum[i]; // 2nd 3 terms are just conductive portion
    }  
}
template void gpuComputeHeatFlux(double *vector, double *jc, double *ke, double *pe, double *stress, 
        double *v, double *tmp, double nktv2p, int *ilist,  int pressatomflag, int dim, int inum);
template void gpuComputeHeatFlux(float *vector, float *jc, float *ke, float *pe, float *stress, float *v, 
        float *tmp, float nktv2p, int *ilist,  int pressatomflag, int dim, int inum);

template <typename T> __global__ void gpuKernelComputeDisplaceAtom(T *displace, T *x, T *xoriginal, T *box, 
        int *pbc, int *ilist,  int triclinic, int dim, int inum)
{    
    int ii = threadIdx.x + blockIdx.x * blockDim.x;
    while (ii < inum) {    
        int i = ilist[ii];   
        T dp[3];
        dp[0] = x[i*dim+0] - xoriginal[i*dim+0];
        dp[1] = x[i*dim+1] - xoriginal[i*dim+1];        
        dp[2] = (dim==3) ? (x[i*dim+2] - xoriginal[i*dim+2]) : 0.0;
        gpuMinimumImage(dp, box, pbc, triclinic, dim);
        displace[ii+0*inum] = dp[0];
        displace[ii+1*inum] = dp[1];
        displace[ii+2*inum] = dp[2];        
        displace[ii+3*inum] = sqrt(dp[0]*dp[0] +dp[1]*dp[1] + dp[2]*dp[2]);
        ii += blockDim.x * gridDim.x;
    }
}
template <typename T> void gpuComputeDisplaceAtom(T *displace, T *x, T *xoriginal, T *box, 
        int *pbc, int *ilist,  int triclinic, int dim, int inum)
{
    int blockDim = 256;
    int gridDim = (inum + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;
    gpuKernelComputeDisplaceAtom<<<gridDim, blockDim>>>(displace, x, xoriginal, box, pbc, 
            ilist, triclinic, dim, inum);
}
template void gpuComputeDisplaceAtom(double *displace, double *x, double *xoriginal, double *box, 
        int *pbc, int *ilist,  int triclinic, int dim, int inum);
template void gpuComputeDisplaceAtom(float *displace, float *x, float *xoriginal, float *box, 
        int *pbc, int *ilist,  int triclinic, int dim, int inum);

template <typename T> __device__ void gpuSelect3(int k, int n, T *arr, int *iarr, T *arr3)
{
  int dim=3;
  int i,ir,j,l,mid,ia,itmp;
  T a,tmp,a3[3];

  arr--;
  iarr--;
  for (i=0; i<3; i++)
    arr3--;
  l = 1;
  ir = n;
  for (;;) {
    if (ir <= l+1) {
      if (ir == l+1 && arr[ir] < arr[l]) {
        SWAP(arr[l],arr[ir]);
        ISWAP(iarr[l],iarr[ir]);
        SWAP(arr3[dim*l+0],arr3[dim*ir+0]);
        SWAP(arr3[dim*l+1],arr3[dim*ir+1]);
        SWAP(arr3[dim*l+2],arr3[dim*ir+2]);
        //SWAP3(arr3[l],arr3[ir]);        
      }
      return;
    } else {
      mid=(l+ir) >> 1;
      SWAP(arr[mid],arr[l+1]);
      ISWAP(iarr[mid],iarr[l+1]);
      SWAP(arr3[dim*mid+0],arr3[dim*(l+1)+0]);
      SWAP(arr3[dim*mid+1],arr3[dim*(l+1)+1]);
      SWAP(arr3[dim*mid+2],arr3[dim*(l+1)+2]);
      //SWAP3(arr3[mid],arr3[l+1]);      
      if (arr[l] > arr[ir]) {
        SWAP(arr[l],arr[ir]);
        ISWAP(iarr[l],iarr[ir]);
        SWAP(arr3[dim*l+0],arr3[dim*ir+0]);
        SWAP(arr3[dim*l+1],arr3[dim*ir+1]);
        SWAP(arr3[dim*l+2],arr3[dim*ir+2]);
        //SWAP3(arr3[l],arr3[ir]);
      }
      if (arr[l+1] > arr[ir]) {
        SWAP(arr[l+1],arr[ir]);
        ISWAP(iarr[l+1],iarr[ir]);
        //SWAP3(arr3[l+1],arr3[ir]);
        SWAP(arr3[dim*(l+1)+0],arr3[dim*ir+0]);
        SWAP(arr3[dim*(l+1)+1],arr3[dim*ir+1]);
        SWAP(arr3[dim*(l+1)+2],arr3[dim*ir+2]);
      }
      if (arr[l] > arr[l+1]) {
        SWAP(arr[l],arr[l+1]);
        ISWAP(iarr[l],iarr[l+1]);
        SWAP(arr3[dim*l+0],arr3[dim*(l+1)+0]);
        SWAP(arr3[dim*l+1],arr3[dim*(l+1)+1]);
        SWAP(arr3[dim*l+2],arr3[dim*(l+1)+2]);
        //SWAP3(arr3[l],arr3[l+1]);
      }
      i = l+1;
      j = ir;
      a = arr[l+1];
      ia = iarr[l+1];
      a3[0] = arr3[dim*(l+1)+0];
      a3[1] = arr3[dim*(l+1)+1];
      a3[2] = arr3[dim*(l+1)+2];
      for (;;) {
        do i++; while (arr[i] < a);
        do j--; while (arr[j] > a);
        if (j < i) break;
        SWAP(arr[i],arr[j]);
        ISWAP(iarr[i],iarr[j]);
        SWAP(arr3[dim*i+0],arr3[dim*j+0]);
        SWAP(arr3[dim*i+1],arr3[dim*j+1]);
        SWAP(arr3[dim*i+2],arr3[dim*j+2]);
        //SWAP3(arr3[i],arr3[j]);
      }
      arr[l+1] = arr[j];
      arr[j] = a;
      iarr[l+1] = iarr[j];
      iarr[j] = ia;
      arr3[dim*(l+1)+0] = arr3[j*dim+0];
      arr3[dim*(l+1)+1] = arr3[j*dim+1];
      arr3[dim*(l+1)+2] = arr3[j*dim+2];
      arr3[j*dim+0] = a3[0];
      arr3[j*dim+1] = a3[1];
      arr3[j*dim+2] = a3[2];      
      if (j >= k) ir = j-1;
      if (j <= k) l = i;
    }
  }
}

template <typename T> __device__ T gpuAssociatedLegendre(int l, int m, T x)
{
  if (l < m) return 0.0;

  T p = 1.0;
  T pm1 = 0.0;
  T pm2 = 0.0;

  if (m != 0) {
    T msqx = -sqrt(1.0-x*x);
    for (int i=1; i < m+1; ++i)
      p *= (2*i-1) * msqx;
  }

  for (int i=m+1; i < l+1; ++i) {
    pm2 = pm1;
    pm1 = p;
    p = ((2*i-1)*x*pm1 - (i+m-1)*pm2) / ((T) (i-m));
  }

  return p;
}

template <typename T> __device__ T gpuPolarPrefactor(int l, int m, T MY_4PI, T costheta)
{
  const int mabs = abs(m);

  T prefactor = 1.0;
  for (int i=l-mabs+1; i < l+mabs+1; ++i)
    prefactor *= (T) i;

  prefactor = sqrt(static_cast<T>(2*l+1)/(MY_4PI*prefactor))
    * gpuAssociatedLegendre(l,mabs,costheta);

  if ((m < 0) && (m % 2)) prefactor = -prefactor;

  return prefactor;
}

template <typename T> __device__ void gpuCalcBoop(T *qn, T *rlist, T *cglist, T *qnm_r, T *qnm_i,
        T MY_EPSILON, T QEPSILON, T MY_4PI, int *qlist, int nqlist, int qmax,  int wlflag, 
        int wlhatflag, int qlcompflag, int iqlcomp, int qlcomp, int ncount) 
{

  int mmax = 2*qmax+1;
  for (int il = 0; il < nqlist; il++) {
    int l = qlist[il];
    for(int m = 0; m < 2*l+1; m++) {
      qnm_r[m+mmax*il] = 0.0;
      qnm_i[m+mmax*il] = 0.0;
    }
  }

  for(int ineigh = 0; ineigh < ncount; ineigh++) {
    const T * const r = &rlist[3*ineigh];
    T rmag = sqrt(r[0]*r[0]+r[1]*r[1]+r[2]*r[2]);
    if(rmag <= MY_EPSILON) {
      return;
    }

    T costheta = r[2] / rmag;
    T expphi_r = r[0];
    T expphi_i = r[1];
    T rxymag = sqrt(expphi_r*expphi_r+expphi_i*expphi_i);
    if(rxymag <= MY_EPSILON) {
      expphi_r = 1.0;
      expphi_i = 0.0;
    } else {
      T rxymaginv = 1.0/rxymag;
      expphi_r *= rxymaginv;
      expphi_i *= rxymaginv;
    }

    for (int il = 0; il < nqlist; il++) {
      int l = qlist[il];

      // calculate spherical harmonics
      // Ylm, -l <= m <= l
      // sign convention: sign(Yll(0,0)) = (-1)^l

      qnm_r[il*mmax+l] += gpuPolarPrefactor(l, 0, MY_4PI, costheta);
      T expphim_r = expphi_r;
      T expphim_i = expphi_i;
      for(int m = 1; m <= +l; m++) {

        T prefactor = gpuPolarPrefactor(l, m, MY_4PI, costheta);
        T ylm_r = prefactor * expphim_r;
        T ylm_i = prefactor * expphim_i;
        qnm_r[il*mmax+m+l] += ylm_r;
        qnm_i[il*mmax+m+l] += ylm_i;
        if(m & 1) {
          qnm_r[il*mmax-m+l] -= ylm_r;
          qnm_i[il*mmax-m+l] += ylm_i;
        } else {
          qnm_r[il*mmax-m+l] += ylm_r;
          qnm_i[il*mmax-m+l] -= ylm_i;
        }
        T tmp_r = expphim_r*expphi_r - expphim_i*expphi_i;
        T tmp_i = expphim_r*expphi_i + expphim_i*expphi_r;
        expphim_r = tmp_r;
        expphim_i = tmp_i;
      }

    }
  }

  // convert sums to averages
  T facn = 1.0 / ncount;
  for (int il = 0; il < nqlist; il++) {
    int l = qlist[il];
    for(int m = 0; m < 2*l+1; m++) {
      qnm_r[il*mmax+m] *= facn;
      qnm_i[il*mmax+m] *= facn;
    }
  }

  // calculate Q_l
  // NOTE: optional W_l_hat and components of Q_qlcomp use these stored Q_l values

  int jj = 0;
  for (int il = 0; il < nqlist; il++) {
    int l = qlist[il];
    T qnormfac = sqrt(MY_4PI/(2*l+1));
    T qm_sum = 0.0;
    for(int m = 0; m < 2*l+1; m++)
      qm_sum += qnm_r[il*mmax+m]*qnm_r[il*mmax+m] + qnm_i[il*mmax+m]*qnm_i[il*mmax+m];
    qn[jj++] = qnormfac * sqrt(qm_sum);
  }

  // calculate W_l

  if (wlflag) {
    int idxcg_count = 0;
    for (int il = 0; il < nqlist; il++) {
      int l = qlist[il];
      T wlsum = 0.0;
      for(int m1 = 0; m1 < 2*l+1; m1++) {
        for(int m2 = max(0,l-m1); m2 < min(2*l+1,3*l-m1+1); m2++) {
          int m = m1 + m2 - l;
          T qm1qm2_r = qnm_r[il*mmax+m1]*qnm_r[il*mmax+m2] - qnm_i[il*mmax+m1]*qnm_i[il*mmax+m2];
          T qm1qm2_i = qnm_r[il*mmax+m1]*qnm_i[il*mmax+m2] + qnm_i[il*mmax+m1]*qnm_r[il*mmax+m2];
          wlsum += (qm1qm2_r*qnm_r[il*mmax+m] + qm1qm2_i*qnm_i[il*mmax+m])*cglist[idxcg_count];
          idxcg_count++;
        }
      }
      qn[jj++] = wlsum/sqrt(2*l+1);
    }
  }

  // calculate W_l_hat

  if (wlhatflag) {
    int idxcg_count = 0;
    for (int il = 0; il < nqlist; il++) {
      int l = qlist[il];
      T wlsum = 0.0;
      for(int m1 = 0; m1 < 2*l+1; m1++) {
        for(int m2 = max(0,l-m1); m2 < min(2*l+1,3*l-m1+1); m2++) {
          int m = m1 + m2 - l;
          T qm1qm2_r = qnm_r[il*mmax+m1]*qnm_r[il*mmax+m2] - qnm_i[il*mmax+m1]*qnm_i[il*mmax+m2];
          T qm1qm2_i = qnm_r[il*mmax+m1]*qnm_i[il*mmax+m2] + qnm_i[il*mmax+m1]*qnm_r[il*mmax+m2];
          wlsum += (qm1qm2_r*qnm_r[il*mmax+m] + qm1qm2_i*qnm_i[il*mmax+m])*cglist[idxcg_count];
          idxcg_count++;
        }
      }
      if (qn[il] < QEPSILON)
        qn[jj++] = 0.0;
      else {
        T qnormfac = sqrt(MY_4PI/(2*l+1));
        T qnfac = qnormfac/qn[il];
        qn[jj++] = wlsum/sqrt(2*l+1)*(qnfac*qnfac*qnfac);
      }
    }
  }

  // Calculate components of Q_l/|Q_l|, for l=qlcomp

  if (qlcompflag) {
    int il = iqlcomp;
    int l = qlcomp;
    if (qn[il] < QEPSILON)
      for(int m = 0; m < 2*l+1; m++) {
        qn[jj++] = 0.0;
        qn[jj++] = 0.0;
      }
    else {
      T qnormfac = sqrt(MY_4PI/(2*l+1));
      T qnfac = qnormfac/qn[il];
      for(int m = 0; m < 2*l+1; m++) {
        qn[jj++] = qnm_r[il*mmax+m] * qnfac;
        qn[jj++] = qnm_i[il*mmax+m] * qnfac;
      }
    }
  }
}

template <typename T> __global__ void gpuKernelComputeOrientOrderAtom(T *qnarray, T *x, T *rlist, T *cglist, 
        T *fac, T *qnm_r, T *qnm_i, T *distsq, T cutsq, T MY_EPSILON, T QEPSILON, T MY_4PI, int *neighlist, 
        int *neighnum,  int *ilist, int *qlist, int *nearest,  int nqlist, int qmax, int wlflag,  
        int wlhatflag,  int qlcompflag, int iqlcomp, int qlcomp, int nnn, int ncol, int jnum, int dim, int inum) 
{
  int ii = threadIdx.x + blockIdx.x * blockDim.x;       
  while (ii < inum) {
      int i = ilist[ii];
      T *qn = &qnarray[ncol*i];    
      T xtmp = x[i*dim+0];
      T ytmp = x[i*dim+1];
      T ztmp = x[i*dim+2];
      int *jlist = &neighlist[jnum*i];
      int m = neighnum[i];
            
      int ncount = 0;
      for (int jj = 0; jj < m; jj++) {
        int j = jlist[jj];
        //j &= NEIGHMASK;

        T delx = xtmp - x[j*dim+0];
        T dely = ytmp - x[j*dim+1];
        T delz = ztmp - x[j*dim+2];
        T rsq = delx*delx + dely*dely + delz*delz;
        if (rsq < cutsq) {
          distsq[ncount] = rsq;
          rlist[ncount*dim+0] = delx;
          rlist[ncount*dim+1] = dely;
          rlist[ncount*dim+2] = delz;
          nearest[ncount++] = j;
        }
      }

      // if not nnn neighbors, order parameter = 0;
      if ((ncount == 0) || (ncount < nnn)) {
        for (int jj = 0; jj < ncol; jj++)
          qn[jj] = 0.0;
        continue;
      }

      // if nnn > 0, use only nearest nnn neighbors
      if (nnn > 0) {
        gpuSelect3(nnn,ncount,distsq,nearest,rlist);
        ncount = nnn;
      }
      
      gpuCalcBoop(qn, rlist, cglist, qnm_r, qnm_i, MY_EPSILON, QEPSILON, MY_4PI,
              qlist, nqlist, qmax, wlflag, wlhatflag, qlcompflag, iqlcomp, qlcomp, ncount);          

     ii += blockDim.x * gridDim.x;
  }
}
template <typename T> void gpuComputeOrientOrderAtom(T *qnarray, T *x, T *rlist, T *cglist, T *fac, 
        T *qnm_r, T *qnm_i, T *distsq, T cutsq, T MY_EPSILON, T QEPSILON, T MY_4PI, int *neighlist, 
        int *neighnum,  int *ilist, int *qlist, int *nearest,  int nqlist, int qmax, int wlflag,  
        int wlhatflag,  int qlcompflag, int iqlcomp, int qlcomp, int nnn, int jnum, int dim, int inum) 
{
    int ncol = nqlist;
    if (wlflag) ncol += nqlist;
    if (wlhatflag) ncol += nqlist;
    if (qlcompflag) ncol += 2*(2*qlcomp+1);    
    
    int blockDim = 256;
    int gridDim = (inum + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;
    gpuKernelComputeOrientOrderAtom<<<gridDim, blockDim>>>(qnarray, x, rlist, cglist, fac, qnm_r, qnm_i, 
        distsq, cutsq, MY_EPSILON, QEPSILON, MY_4PI, neighlist, neighnum, ilist, qlist, nearest,  
        nqlist, qmax, wlflag, wlhatflag, qlcompflag, iqlcomp, qlcomp, nnn, ncol, jnum, dim, inum); 
        
}
template void gpuComputeOrientOrderAtom(double *qnarray, double *x, double *rlist, double *cglist, double *fac, 
        double *qnm_r, double *qnm_i, double *distsq, double cutsq, double MY_EPSILON, double QEPSILON, 
        double MY_4PI, int *neighlist, int *neighnum,  int *ilist, int *qlist, int *nearest,  int nqlist, 
        int qmax, int wlflag, int wlhatflag,  int qlcompflag, int iqlcomp, int qlcomp, int nnn, int jnum, int dim, int inum); 
template void gpuComputeOrientOrderAtom(float *qnarray, float *x, float *rlist, float *cglist, float *fac, 
        float *qnm_r, float *qnm_i, float *distsq, float cutsq, float MY_EPSILON, float QEPSILON, float MY_4PI, 
        int *neighlist, int *neighnum,  int *ilist, int *qlist, int *nearest,  int nqlist, int qmax, int wlflag,  
        int wlhatflag,  int qlcompflag, int iqlcomp, int qlcomp, int nnn, int jnum, int dim, int inum); 

template <typename T> __global__ void gpuKernelComputeCoordAtomCutoff(int *cvec, T *x, T *rcutsq, 
        int *type, int *ilist, int *neighlist, int *neighnum, int *typelo, int *typehi, 
         int *jgroupbit, int dim, int ntypes, int jnum, int inum)
{
   int ii = threadIdx.x + blockIdx.x * blockDim.x;    
   while (ii < inum) {
      int i = ilist[ii];
      T xtmp = x[i*dim];
      T ytmp = x[i*dim+1];
      T ztmp = x[i*dim+2];
      int itype = type[i];
      int kk = neighnum[i];
      int n = 0;
      for (int jj = 0; jj < kk; jj++) {
        int j = neighlist[jj + jnum*i];     
        //j &= NEIGHMASK;
        if (jgroupbit[j]) {
          int jtype = type[j];
          T delx = xtmp - x[j*dim+0];
          T dely = ytmp - x[j*dim+1];
          T delz = ztmp - x[j*dim+2];
          T rsq = delx*delx + dely*dely + delz*delz;
          if (rsq < rcutsq[(jtype-1) + (itype-1)*ntypes] && jtype >= typelo[0] && jtype <= typehi[0])
            n++;
        }
      }
      cvec[i] = n;    
      ii += blockDim.x * gridDim.x;
   }
}
template <typename T> void gpuComputeCoordAtomCutoff(int *cvec, T *x, T *rcutsq, int *type, 
        int *ilist, int *neighlist, int *neighnum, int *typelo, int *typehi, 
         int *jgroupbit, int dim, int ntypes, int jnum, int inum)
{
    int blockDim = 256;
    int gridDim = (inum + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;
    gpuKernelComputeCoordAtomCutoff<<<gridDim, blockDim>>>(cvec, x, rcutsq, type, ilist, 
       neighlist, neighnum, typelo, typehi, jgroupbit, dim, ntypes, jnum, inum);     
}
template void gpuComputeCoordAtomCutoff(int *cvec, double *x, double *rcutsq, int *type, 
        int *ilist, int *neighlist, int *neighnum, int *typelo, int *typehi, 
         int *jgroupbit, int dim, int ntypes, int jnum, int inum);
template void gpuComputeCoordAtomCutoff(int *cvec, float *x, float *rcutsq, int *type, 
        int *ilist, int *neighlist, int *neighnum, int *typelo, int *typehi, 
         int *jgroupbit, int dim, int ntypes, int jnum, int inum);

template <typename T> __global__ void gpuKernelComputeCoordAtomCutoff(int *carray, T *x, T *rcutsq, 
        int *type, int *ilist, int *neighlist, int *neighnum, int *typelo, int *typehi, 
         int *jgroupbit, int ncol, int dim, int ntypes, int jnum, int inum)
{
   int ii = threadIdx.x + blockIdx.x * blockDim.x;     
   while (ii < inum) {
      int i = ilist[ii];    
      for (int m = 0; m < ncol; m++) carray[m+ncol*i] = 0;    
      T xtmp = x[i*dim];
      T ytmp = x[i*dim+1];
      T ztmp = x[i*dim+2];
      int itype = type[i];
      int kk = neighnum[i];
      int n = 0;
      for (int jj = 0; jj < kk; jj++) {
        int j = neighlist[jj + jnum*i];     
        //j &= NEIGHMASK;
        if (jgroupbit[j]) {
          int jtype = type[j];
          T delx = xtmp - x[j*dim+0];
          T dely = ytmp - x[j*dim+1];
          T delz = ztmp - x[j*dim+2];
          T rsq = delx*delx + dely*dely + delz*delz;
          if (rsq < rcutsq[(jtype-1) + (itype-1)*ntypes]) {
            for (int m = 0; m < ncol; m++)
                if (jtype >= typelo[m] && jtype <= typehi[m])
                    carray[m+ncol*i] += 1;                                  
          }          
        }
      }    
      ii += blockDim.x * gridDim.x;
   }
}
template <typename T> void gpuComputeCoordAtomCutoff(int *carray, T *x, T *rcutsq, int *type, 
        int *ilist, int *neighlist, int *neighnum, int *typelo, int *typehi, 
         int *jgroupbit, int ncol, int dim, int ntypes, int jnum, int inum)
{
    int blockDim = 256;
    int gridDim = (inum + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;
    gpuKernelComputeCoordAtomCutoff<<<gridDim, blockDim>>>(carray, x, rcutsq, type, ilist, 
       neighlist, neighnum, typelo, typehi, jgroupbit, ncol, dim, ntypes, jnum, inum);        
}
template void gpuComputeCoordAtomCutoff(int *carray, double *x, double *rcutsq, int *type, 
        int *ilist, int *neighlist, int *neighnum, int *typelo, int *typehi, 
         int *jgroupbit, int ncol, int dim, int ntypes, int jnum, int inum);
template void gpuComputeCoordAtomCutoff(int *carray, float *x, float *rcutsq, int *type, 
        int *ilist, int *neighlist, int *neighnum, int *typelo, int *typehi, 
         int *jgroupbit, int ncol, int dim, int ntypes, int jnum, int inum);

template <typename T> __global__ void gpuKernelComputeCoordAtomOrient(int *cvec, T *x, T *rcutsq, 
        T *normv, T threshold, int *type, int *ilist, int *neighlist, int *neighnum, int *typelo, 
        int *typehi, int *jgroupbit, int dim, int ntypes, int nqlist, int ncol, int l, int jnum, int inum)
{
   int ii = threadIdx.x + blockIdx.x * blockDim.x;      
   while (ii < inum) {
      int i = ilist[ii];
      T xtmp = x[i*dim];
      T ytmp = x[i*dim+1];
      T ztmp = x[i*dim+2];
      int itype = type[i];
      int kk = neighnum[i];
      int n = 0;
      for (int jj = 0; jj < kk; jj++) {
        int j = neighlist[jj + jnum*i];     
        //j &= NEIGHMASK;
        if (jgroupbit[j]) {
          int jtype = type[j];
          T delx = xtmp - x[j*dim+0];
          T dely = ytmp - x[j*dim+1];
          T delz = ztmp - x[j*dim+2];
          T rsq = delx*delx + dely*dely + delz*delz;
          if (rsq < rcutsq[(jtype-1) + (itype-1)*ntypes]) {
            T dot_product = 0.0;
            for (int m=0; m < 2*(2*l+1); m++) {
              dot_product += normv[(nqlist+m) + ncol*i]*normv[(nqlist+m) +  ncol*j];
            }
            if (dot_product > threshold) n++;
          }          
        }
      }
      cvec[i] = n;   
      ii += blockDim.x * gridDim.x;
   }
}
template <typename T> void gpuComputeCoordAtomOrient(int *cvec, T *x, T *rcutsq, T *normv, T threshold, 
        int *type, int *ilist, int *neighlist, int *neighnum, int *typelo, int *typehi, 
         int *jgroupbit, int dim, int ntypes, int nqlist, int ncol, int l, int jnum, int inum)
{
    int blockDim = 256;
    int gridDim = (inum + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;
    gpuKernelComputeCoordAtomOrient<<<gridDim, blockDim>>>(cvec, x, rcutsq, normv, threshold, type, ilist, 
       neighlist, neighnum, typelo, typehi, jgroupbit, dim, ntypes, nqlist, ncol, l, jnum, inum);            
}
template void gpuComputeCoordAtomOrient(int *cvec, double *x, double *rcutsq, double *normv, double threshold, 
        int *type, int *ilist, int *neighlist, int *neighnum, int *typelo, int *typehi, 
        int *jgroupbit, int dim, int ntypes, int nqlist, int ncol, int l, int jnum, int inum);
template void gpuComputeCoordAtomOrient(int *cvec, float *x, float *rcutsq, float *normv, float threshold, 
        int *type, int *ilist, int *neighlist, int *neighnum, int *typelo, int *typehi, 
        int *jgroupbit, int dim, int ntypes, int nqlist, int ncol, int l, int jnum, int inum);

template <typename T> __global__ void gpuKernelComputeMSD(T *vec, T *x, T *xoriginal, T *box, T *xcm,
         T navfac, int *ilist, int *image, int naverage, int avflag, int triclinic, int nmsd, int dim,  int inum)
{
    int ii = threadIdx.x + blockIdx.x * blockDim.x;      
    while (ii < inum) {
        int i = ilist[ii];
        if (dim==3) {
            T xtmp, ytmp, ztmp;
            if (triclinic == 0) {
                xtmp = x[i*dim+0] + image[i*dim+0]*box[0] - xcm[0];
                ytmp = x[i*dim+1] + image[i*dim+1]*box[1] - xcm[1];
                ztmp = x[i*dim+2] + image[i*dim+2]*box[2] - xcm[2];
            } else {
                xtmp = x[i*dim+0] + box[0]*image[i*dim+0] + box[5]*image[i*dim+1] + box[4]*image[i*dim+2] - xcm[0];
                ytmp = x[i*dim+1] + box[1]*image[i*dim+1] + box[3]*image[i*dim+2] - xcm[1];
                ztmp = x[i*dim+2] + box[2]*image[i*dim+2] - xcm[2];            
            }
            if (avflag) {
              xoriginal[i*dim+0] = (xoriginal[i*dim+0]*naverage + xtmp)*navfac;
              xoriginal[i*dim+1] = (xoriginal[i*dim+1]*naverage + ytmp)*navfac;
              xoriginal[i*dim+2] = (xoriginal[i*dim+2]*naverage + ztmp)*navfac;
            }
            T dx = xtmp - xoriginal[i*dim+0];
            T dy = ytmp - xoriginal[i*dim+1];
            T dz = ztmp - xoriginal[i*dim+2];
            vec[ii+0*inum] = dx*dx;
            vec[ii+1*inum] = dy*dy;
            vec[ii+2*inum] = dz*dz;
        } else {
            T xtmp, ytmp;
            if (triclinic == 0) {
                xtmp = x[i*dim+0] + image[i*dim+0]*box[0] - xcm[0];
                ytmp = x[i*dim+1] + image[i*dim+1]*box[1] - xcm[1];
            } else {
                xtmp = x[i*dim+0] + box[0]*image[i*dim+0] + box[5]*image[i*dim+1] - xcm[0];
                ytmp = x[i*dim+1] + box[1]*image[i*dim+1] - xcm[1];
            }
            if (avflag) {
              xoriginal[i*dim+0] = (xoriginal[i*dim+0]*naverage + xtmp)*navfac;
              xoriginal[i*dim+1] = (xoriginal[i*dim+1]*naverage + ytmp)*navfac;
            }
            T dx = xtmp - xoriginal[i*dim+0];
            T dy = ytmp - xoriginal[i*dim+1];
            vec[ii+0*inum] = dx*dx;
            vec[ii+1*inum] = dy*dy;
        }
        ii += blockDim.x * gridDim.x;
    }    
  //MPI_Allreduce(msd,vector,4,MPI_DOUBLE,MPI_SUM,world);  
}
template <typename T> void gpuComputeMSD(T *msd, T *vec, T *x, T *xoriginal, T *box, T *xcm, T *tmp,
         int *ilist, int *image, int naverage, int avflag, int triclinic, int nmsd, int dim,  int inum)
{
    T navfac;  
    if (avflag) {
        naverage += 1;
        navfac = 1.0/(naverage+1);
    }
    
    int blockDim = 256;
    int gridDim = (inum + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;
    gpuKernelComputeMSD<<<gridDim, blockDim>>>(vec, x, xoriginal, box, xcm, navfac,
         ilist, image, naverage, avflag, triclinic, nmsd,  dim,  inum);
    
    msd[2] = 0.0;
    for (int i=0; i<dim; i++)
        msd[i] = gpuArraySum(&vec[inum*i], tmp, inum);        
        
    msd[0] /= nmsd;
    msd[1] /= nmsd;
    msd[2] /= nmsd;
    msd[3] = msd[0] + msd[1] + msd[2];          
}
template void gpuComputeMSD(double *msd, double *vec, double *x, double *xoriginal, double *box, double *xcm,
       double *tmp, int *ilist, int *image, int naverage, int avflag, int triclinic, int nmsd, int dim,  int inum);
template void gpuComputeMSD(float *msd, float *vec, float *x, float *xoriginal, float *box, float *xcm,
       float *tmp,  int *ilist, int *image, int naverage, int avflag, int triclinic, int nmsd, int dim,  int inum);

template <typename T> __global__ void gpuKernelComputeVACF(T *vec, T *v, T *voriginal, 
         int *ilist, int dim,  int inum)
{  
  int ii = threadIdx.x + blockIdx.x * blockDim.x;        
  while (ii < inum) {
      int i = ilist[ii];
      vec[ii+0*inum] = v[i*dim+0] * voriginal[i*dim+0];
      vec[ii+1*inum] = v[i*dim+1] * voriginal[i*dim+1];
      if (dim==3)      
        vec[ii+2*inum] = v[i*dim+2] * voriginal[i*dim+2];
      ii += blockDim.x * gridDim.x;
  }
  //MPI_Allreduce(vacf,vector,4,MPI_DOUBLE,MPI_SUM,world);
}
template <typename T> void gpuComputeVACF(T *vacf, T *vec, T *v, T *voriginal, T *tmp,
         int *ilist, int nvacf, int dim,  int inum)
{
    
    int blockDim = 256;
    int gridDim = (inum + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;
    gpuKernelComputeVACF<<<gridDim, blockDim>>>(vec, v, voriginal, 
         ilist, dim, inum);    

    vacf[2] = 0.0;
    for (int i=0; i<dim; i++)
        vacf[i] = gpuArraySum(&vec[inum*i], tmp, inum);        
        
    vacf[0] /= nvacf;
    vacf[1] /= nvacf;
    vacf[2] /= nvacf;
    vacf[3] = vacf[0] + vacf[1] + vacf[2];          
}
template void gpuComputeVACF(double *vacf, double *vec, double *v, double *voriginal, double *tmp,
         int *ilist, int nvacf, int dim,  int inum);
template void gpuComputeVACF(float *vacf, float *vec, float *v, float *voriginal, float *tmp,
         int *ilist, int nvacf, int dim,  int inum);


#endif

