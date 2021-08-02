/***************************************************************************                               
                    Molecular Dynamics Potentials (MDP)
                           CESMIX-MIT Project  
 
 Contributing authors: Ngoc-Cuong Nguyen (cuongng@mit.edu, exapde@gmail.com)
 ***************************************************************************/

template <typename T> __global__ void gpuKernelLamda2Box(T *x, T *lambda, T *h, T *boxlo, int dim, int n)
{    
    int i = threadIdx.x + blockIdx.x * blockDim.x;
    while (i < n) {
        int k = i*dim;
        x[k+0] = h[0]*lambda[k+0] + h[5]*lambda[k+1] + boxlo[0];
        x[k+1] = h[1]*lambda[k+1] + boxlo[1];    
        if (dim==3) {
            x[k+0] += h[4]*lambda[k+2];
            x[k+1] += h[3]*lambda[k+2];    
            x[k+2] = h[2]*lambda[k+2] + boxlo[2];
        }
        i += blockDim.x * gridDim.x;
    }        
}
template <typename T> void gpuLamda2Box(T *x, T *lambda, T *h, T *boxlo, int dim, int n)
{        
    int blockDim = 256;
    int gridDim = (n + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;
    gpuKernelLamda2Box<<<gridDim, blockDim>>>(x, lambda, h, boxlo, dim, n);
}
template void gpuLamda2Box(double *x, double *lambda, double *h, double *boxlo, int dim, int n);
template void gpuLamda2Box(float *x, float *lambda, float *h, float *boxlo, int dim, int n);

template <typename T> __global__ void gpuKernelBox2Lamda(T *lambda, T *x, T *h_inv, T *boxlo, int dim, int n)
{
  int i = threadIdx.x + blockIdx.x * blockDim.x;
    while (i < n) {
    int k = i*dim;
    T deltax = x[k+0] - boxlo[0];
    T deltay = x[k+1] - boxlo[1];
    lambda[k+0] = h_inv[0]*deltax + h_inv[5]*deltay;
    lambda[k+1] = h_inv[1]*deltay;    
    if (dim==3) {
        T deltaz = x[k+2] - boxlo[2];
        lambda[k+0] += h_inv[4]*deltaz;
        lambda[k+1] += h_inv[3]*deltaz;
        lambda[k+2] = h_inv[2]*deltaz;
    }
    i += blockDim.x * gridDim.x;
  }
}
template <typename T> void gpuBox2Lamda(T *lambda, T *x, T *h_inv, T *boxlo, int dim, int n)
{        
    int blockDim = 256;
    int gridDim = (n + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;
    gpuKernelBox2Lamda<<<gridDim, blockDim>>>(lambda, x, h_inv, boxlo, dim, n);
}
template void gpuBox2Lamda(double *lambda, double *x, double *h_inv, double *boxlo, int dim, int n);
template void gpuBox2Lamda(float *lambda, float *x, float *h_inv, float *boxlo, int dim, int n);

template <typename T> __global__ void gpuKernelPBCOrthBox(T *x, T *v, int *image, T *hi, T *lo, 
        T *h, T *h_rate, int *pbc, int vdeform, int dim, int nlocal)
{  
  for (int i = 0; i < nlocal; i++) {
    int idim,otherdims;
    int k = i*dim;
    if (pbc[0]) {
      if (x[k+0] < lo[0]) {
        x[k+0] += h[0];
        if (vdeform) v[k+0] += h_rate[0];
        image[k+0] -= 1;
      }
      if (x[k+0] >= hi[0]) {
        x[k+0] -= h[0];
        x[k+0] = MAX(x[k+0],lo[0]);
        if (vdeform) v[k+0] -= h_rate[0];        
        image[k+0] += 1;
      }
    }

    if (pbc[1]) {
      if (x[k+1] < lo[1]) {
        x[k+1] += h[1];
        if (vdeform) {
          v[k+0] += h_rate[5];
          v[k+1] += h_rate[1];
        }        
        image[k+1] -= 1;
      }
      if (x[k+1] >= hi[1]) {
        x[k+1] -= h[1];
        x[k+1] = MAX(x[k+1],lo[1]);
        if (vdeform) {
          v[k+0] -= h_rate[5];
          v[k+1] -= h_rate[1];
        }        
        image[k+1] += 1;
      }
    }

    if (pbc[2]) {
      if (x[k+2] < lo[2]) {
        x[k+2] += h[2];
        if (vdeform) {
          v[k+0] += h_rate[4];
          v[k+1] += h_rate[3];
          v[k+2] += h_rate[2];
        }        
        image[k+2] -= 1;
      }
      if (x[k+2] >= hi[2]) {
        x[k+2] -= h[2];
        x[k+2] = MAX(x[k+2],lo[2]);
        if (vdeform) {
          v[k+0] -= h_rate[4];
          v[k+1] -= h_rate[3];
          v[k+2] -= h_rate[2];
        }
        image[k+2] += 1;
      }
    }
  }
}

template <typename T> __global__ void gpuKernelPBCTiltBox(T *x, T *v, int *image, T *boxlo, T *h, T *h_inv, T *h_rate, 
        T *hi_lambda, T *lo_lambda, int *pbc, int vdeform, int dim, int nlocal)
{  
  for (int i = 0; i < nlocal; i++) {
    // map x to lambda  
    int k = i*dim;  
    T deltax = x[k+0] - boxlo[0];
    T deltay = x[k+1] - boxlo[1];
    int idim,otherdims;
    T lambda[3];
    lambda[0] = h_inv[0]*deltax + h_inv[5]*deltay;
    lambda[1] = h_inv[1]*deltay;        
    if (dim==3) {
        T deltaz = x[k+2] - boxlo[2];
        lambda[0] += h_inv[4]*deltaz;
        lambda[1] += h_inv[3]*deltaz;
        lambda[2] = h_inv[2]*deltaz;
    }
    
    if (pbc[0]) {
      if (lambda[0] < lo_lambda[0]) {
        lambda[0] += h[0];
        if (vdeform) v[k+0] += h_rate[0];
        image[k+0] -= 1;
      }
      if (lambda[0] >= hi_lambda[0]) {
        lambda[0] -= h[0];
        lambda[0] = MAX(lambda[0],lo_lambda[0]);
        if (vdeform) v[k+0] -= h_rate[0];
        image[k+0] += 1;
      }
    }

    if (pbc[1]) {
      if (lambda[1] < lo_lambda[1]) {
        lambda[1] += h[1];
        if (vdeform) {
          v[k+0] += h_rate[5];
          v[k+1] += h_rate[1];
        }
        image[k+1] -= 1;
      }
      if (lambda[1] >= hi_lambda[1]) {
        lambda[1] -= h[1];
        lambda[1] = MAX(lambda[1],lo_lambda[1]);
        if (vdeform) {
          v[k+0] -= h_rate[5];
          v[k+1] -= h_rate[1];
        }
        image[k+1] += 1;
      }
    }

    if (pbc[2]) {
      if (lambda[2] < lo_lambda[2]) {
        lambda[2] += h[2];
        if (vdeform) {
          v[k+0] += h_rate[4];
          v[k+1] += h_rate[3];
          v[k+2] += h_rate[2];
        }
        image[k+2] -= 1;
      }
      if (lambda[2] >= hi_lambda[2]) {
        lambda[2] -= h[2];
        lambda[2] = MAX(lambda[2],lo_lambda[2]);
        if (vdeform) {
          v[k+0] -= h_rate[4];
          v[k+1] -= h_rate[3];
          v[k+2] -= h_rate[2];
        }
        image[k+2] += 1;
      }
    }
    
    // map lambda  to x
    x[k+0] = h[0]*lambda[0] + h[5]*lambda[1] + boxlo[0];
    x[k+1] = h[1]*lambda[1] + boxlo[1];    
    if (dim==3) {
        x[k+0] += h[4]*lambda[2];
        x[k+1] += h[3]*lambda[2];    
        x[k+2] = h[2]*lambda[2] + boxlo[2];
    }    
  }
}

template <typename T> void gpuPBC(T *x, T *v, int *image, T *boxhi, T *boxlo, T *hi_lambda, T *lo_lambda,  
        T *h, T *h_inv, T *h_rate, int *pbc, int vdeform, int triclinic, int dim, int nlocal)
{
    int blockDim = 256;
    int gridDim = (nlocal + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;        
    if (triclinic==0) 
        gpuKernelPBCOrthBox<<<gridDim, blockDim>>>(x, v, image, boxhi, boxlo, h, h_rate, pbc, vdeform, dim, nlocal);
    else
        gpuKernelPBCTiltBox<<<gridDim, blockDim>>>(x, v, image, boxlo, h, h_inv, h_rate, hi_lambda, lo_lambda, pbc, vdeform, dim, nlocal);
}
template void gpuPBC(double *x, double *v, int *image, double *boxhi, double *boxlo, double *hi_lambda, double *lo_lambda,  
        double *h, double *h_inv, double *h_rate, int *pbc, int vdeform, int triclinic, int dim, int nlocal);
template void gpuPBC(float *x, float *v, int *image, float *boxhi, float *boxlo, float *hi_lambda, float *lo_lambda,  
        float *h, float *h_inv, float *h_rate, int *pbc, int vdeform, int triclinic, int dim, int nlocal); 

template <typename T> __global__ void gpuKernelInsideOrthBox(T *inside, T *x, T *lo, T *hi, int dim, int n)
{  
    for (int i=0; i<n; i++) {
        int k = i*dim;
        if (dim==3) {
            if (x[0+k] < lo[0] || x[0+k] >= hi[0] ||
                x[1+k] < lo[1] || x[1+k] >= hi[1] ||
                x[2+k] < lo[2] || x[2+k] >= hi[2]) inside[i] = 0;
            else inside[i] = 1;            
        }
        else {
            if (x[0+k] < lo[0] || x[0+k] >= hi[0] ||
                x[1+k] < lo[1] || x[1+k] >= hi[1]) inside[i] = 0;
            else inside[i] = 1;            
        }
    }
}

template <typename T> __global__ void gpuKernelInsideTiltBox(T *inside, T *x, T *lo_lamda, T *hi_lamda, 
        T *boxlo, T *h_inv, int dim, int n)
{
  for (int i = 0; i < n; i++) {
    int k = i*dim;
    T deltax = x[k+0] - boxlo[0];
    T deltay = x[k+1] - boxlo[1];
    T lamda0 = h_inv[0]*deltax + h_inv[5]*deltay;
    T lamda1 = h_inv[1]*deltay;        
    if (dim==3) {
        T deltaz = x[k+2] - boxlo[2];
        lamda0 += h_inv[4]*deltaz;
        lamda1 += h_inv[3]*deltaz;
        T lamda2 = h_inv[2]*deltaz;
        if (lamda0 < lo_lamda[0] || lamda0 >= hi_lamda[0] ||
            lamda1 < lo_lamda[1] || lamda1 >= hi_lamda[1] ||
            lamda2 < lo_lamda[2] || lamda2 >= hi_lamda[2]) inside[i] = 0;
        else inside[i] = 1;            
    }
    else {
        if (lamda0 < lo_lamda[0] || lamda0 >= hi_lamda[0] ||
            lamda1 < lo_lamda[1] || lamda1 >= hi_lamda[1]) inside[i] = 0;
        else inside[i] = 1;                    
    }
  }    
}
template <typename T> void gpuInsideBox(T *inside, T *x, T *boxlo, T *boxhi, T *lo_lamda, T *hi_lamda, 
        T *h_inv, int triclinic, int dim, int n)
{
    int blockDim = 256;
    int gridDim = (n + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;        
    
    if (triclinic==0) 
        gpuKernelInsideOrthBox<<<gridDim, blockDim>>>(inside, x, boxlo, boxhi, dim, n);
    else
        gpuKernelInsideTiltBox<<<gridDim, blockDim>>>(inside, x, boxlo, lo_lamda, hi_lamda, h_inv, dim, n);
}
template void gpuInsideBox(double *inside, double *x, double *boxlo, double *boxhi, double *lo_lamda, double *hi_lamda, double *h_inv, int triclinic, int dim, int n);
template void gpuInsideBox(float *inside, float *x, float *boxlo, float *boxhi, float *lo_lamda, float *hi_lamda, float *h_inv, int triclinic, int dim, int n);

template <typename T> __device__ void gpuMinimumImage(T *dp, T *h, int *pbc, int triclinic, int dim)
{
    if (triclinic == 0) {
        if (pbc[0]) {
          while (fabs(dp[0]) > 0.5*h[0]) {
            if (dp[0] < 0.0) dp[0] += h[0];
            else dp[0] -= h[0];
          }
        }
        if (pbc[1]) {
          while (fabs(dp[1]) > 0.5*h[1]) {
            if (dp[1] < 0.0) dp[1] += h[1];
            else dp[1] -= h[1];
          }
        }
        if ((dim==3) && pbc[2]) {
          while (fabs(dp[2]) > 0.5*h[2]) {
            if (dp[2] < 0.0) dp[2] += h[2];
            else dp[2] -= h[2];
          }
        }        
    }
    else {
        if ((dim==3) && pbc[2]) {
          while (fabs(dp[2]) > 0.5*h[2]) {
            if (dp[2] < 0.0) {
              dp[2] += h[2];
              dp[1] += h[5];
              dp[0] += h[4];
            } else {
              dp[2] -= h[2];
              dp[1] -= h[5];
              dp[0] -= h[4];
            }
          }
        }
        if (pbc[1]) {
          while (fabs(dp[1]) > 0.5*h[1]) {
            if (dp[1] < 0.0) {
              dp[1] += h[1];
              dp[0] += h[3];
            } else {
              dp[1] -= h[1];
              dp[0] -= h[3];
            }
          }
        }
        if (pbc[0]) {
          while (fabs(dp[0]) > 0.5*h[0]) {
            if (dp[0] < 0.0) dp[0] += h[0];
            else dp[0] -= h[0];
          }
        }          
    }
}
template __device__ void gpuMinimumImage(double *dp, double *h, int *pbc, int triclinic, int dim);
template __device__ void gpuMinimumImage(float *dp, float *h, int *pbc, int triclinic, int dim);

template <typename T> __global__ void gpuKernelMinimumImageOrthBox(T *dp, T *h, int *pbc, int dim, int n)
{
    for (int i=0; i<n; i++) {
        int k = i*dim;
        if (pbc[0]) {
          while (fabs(dp[0+k]) > 0.5*h[0]) {
            if (dp[0+k] < 0.0) dp[0+k] += h[0];
            else dp[0+k] -= h[0];
          }
        }
        if (pbc[1]) {
          while (fabs(dp[1+k]) > 0.5*h[1]) {
            if (dp[1+k] < 0.0) dp[1+k] += h[1];
            else dp[1+k] -= h[1];
          }
        }
        if (pbc[2]) {
          while (fabs(dp[2+k]) > 0.5*h[2]) {
            if (dp[2+k] < 0.0) dp[2+k] += h[2];
            else dp[2+k] -= h[2];
          }
        }
    }
}

template <typename T> __global__ void gpuKernelMinimumImageTiltBox(T *dp, T *h, int *pbc, int dim, int n)
{
    for (int i=0; i<n; i++) {
        int k = i*dim;
        if (pbc[2]) {
          while (fabs(dp[2+k]) > 0.5*h[2]) {
            if (dp[2+k] < 0.0) {
              dp[2+k] += h[2];
              dp[1+k] += h[5];
              dp[0+k] += h[4];
            } else {
              dp[2+k] -= h[2];
              dp[1+k] -= h[5];
              dp[0+k] -= h[4];
            }
          }
        }
        if (pbc[1]) {
          while (fabs(dp[1+k]) > 0.5*h[1]) {
            if (dp[1+k] < 0.0) {
              dp[1+k] += h[1];
              dp[0+k] += h[3];
            } else {
              dp[1+k] -= h[1];
              dp[0+k] -= h[3];
            }
          }
        }
        if (pbc[0]) {
          while (fabs(dp[0+k]) > 0.5*h[0]) {
            if (dp[0+k] < 0.0) dp[0+k] += h[0];
            else dp[0+k] -= h[0];
          }
        }  
    }
}

template <typename T> void gpuMinimumImage(T *dp, T *h, int *pbc, int triclinic, int dim, int n)
{
    int blockDim = 256;
    int gridDim = (n + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;        
        
    if (triclinic == 0)
        gpuKernelMinimumImageOrthBox<<<gridDim, blockDim>>>(dp, h, pbc, dim, n);
    else
        gpuKernelMinimumImageTiltBox<<<gridDim, blockDim>>>(dp, h, pbc, dim, n);
}
template void gpuMinimumImage(double *dp, double *h, int *pbc, int triclinic, int dim, int n);
template void gpuMinimumImage(float *dp, float *h, int *pbc, int triclinic, int dim, int n);

template <typename T> __device__ void gpuUnmap(T *y, T *x, T *h, int *image, int triclinic, int dim)
{
  if (triclinic == 0) {
    y[0] = x[0] + image[0]*h[0];
    y[1] = x[1] + image[1]*h[1];
    if (dim==3)
        y[2] = x[2] + image[2]*h[2];
  } else {
    y[0] = x[0] + h[0]*image[0] + h[5]*image[1];
    y[1] = x[1] + h[1]*image[1];
    if (dim==3) {
        y[0] += h[4]*image[2];
        y[1] += h[3]*image[2];
        y[2] = x[2] + h[2]*image[2];
    }    
  }
}
template __device__ void gpuUnmap(double *y, double *x, double *h, int *image, int triclinic, int dim);
template __device__ void gpuUnmap(float *y, float *x, float *h, int *image, int triclinic, int dim);

template <typename T> __global__ void gpuKernelUnmapOrthBox(T *y, T *x, T *h, int *image, int dim, int n)
{
     for (int i=0; i<n; i++) {
        int k = i*dim; 
        y[k+0] = x[k+0] + image[k+0]*h[0];
        y[k+1] = x[k+1] + image[k+1]*h[1];
        if (dim==3)
            y[k+2] = x[k+2] + image[k+2]*h[2];
     }
}

template <typename T> __global__ void gpuKernelUnmapTiltBox(T *y, T *x, T *h, int *image, int dim, int n)
{
     for (int i=0; i<n; i++) {
        int k = i*dim; 
        y[k+0] = x[k+0] + h[0]*image[k+0] + h[5]*image[k+1];
        y[k+1] = x[k+1] + h[1]*image[k+1];
        if (dim==3) {
            y[k+0] += h[4]*image[k+2];
            y[k+1] += h[3]*image[k+2];
            y[k+2] = x[k+2] + h[2]*image[k+2];
        }            
     }
}

template <typename T> void gpuUnmap(T *y, T *x, T *h, int *image, int triclinic, int dim, int n)
{
    int blockDim = 256;
    int gridDim = (n + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;        
            
    if (triclinic == 0)
        gpuKernelUnmapOrthBox<<<gridDim, blockDim>>>(y, x, h, image, dim, n);
    else
        gpuKernelUnmapTiltBox<<<gridDim, blockDim>>>(y, x, h, image, dim, n);
}
template void gpuUnmap(double *y, double *x, double *h, int *image, int triclinic, int dim, int n);
template void gpuUnmap(float *y, float *x, float *h, int *image, int triclinic, int dim, int n);

