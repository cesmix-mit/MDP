/***************************************************************************                               
                    Molecular Dynamics Potentials (MDP)
                           CESMIX-MIT Project  
 
 Contributing authors: Ngoc-Cuong Nguyen (cuongng@mit.edu, exapde@gmail.com)
 ***************************************************************************/

template <typename T> void cpuSetGlobalBox(T *h, T *h_inv, T *boxlo_bound, 
        T *boxhi_bound, T *boxhi, T *boxlo, T *boxtilt, int triclinic)
{
  T xprd = boxhi[0] - boxlo[0];
  T yprd = boxhi[1] - boxlo[1];
  T zprd = boxhi[2] - boxlo[2];

  h[0] = xprd;
  h[1] = yprd;
  h[2] = zprd;
  h_inv[0] = 1.0/h[0];
  h_inv[1] = 1.0/h[1];
  h_inv[2] = 1.0/h[2];

  if (triclinic) {
    T xy = boxtilt[0];
    T xz = boxtilt[1];
    T yz = boxtilt[2];
    h[3] = yz;
    h[4] = xz;
    h[5] = xy;
    h_inv[3] = -h[3] / (h[1]*h[2]);
    h_inv[4] = (h[3]*h[5] - h[1]*h[4]) / (h[0]*h[1]*h[2]);
    h_inv[5] = -h[5] / (h[0]*h[1]);

    boxlo_bound[0] = MIN(boxlo[0],boxlo[0]+xy);
    boxlo_bound[0] = MIN(boxlo_bound[0],boxlo_bound[0]+xz);
    boxlo_bound[1] = MIN(boxlo[1],boxlo[1]+yz);
    boxlo_bound[2] = boxlo[2];

    boxhi_bound[0] = MAX(boxhi[0],boxhi[0]+xy);
    boxhi_bound[0] = MAX(boxhi_bound[0],boxhi_bound[0]+xz);
    boxhi_bound[1] = MAX(boxhi[1],boxhi[1]+yz);
    boxhi_bound[2] = boxhi[2];
  }
  else {
    h[3] = 0;
    h[4] = 0;
    h[5] = 0;
    h_inv[3] = 0;
    h_inv[4] = 0;
    h_inv[5] = 0;      
    for (int i = 0; i<3; i ++) {
        boxhi_bound[i] = boxhi[i];
        boxlo_bound[i] = boxlo[i];
    }
  }
}
template void cpuSetGlobalBox(double *h, double *h_inv, double *boxlo_bound, double *boxhi_bound, double *boxhi, double *boxlo, double *boxtilt, int triclinic);
template void cpuSetGlobalBox(float *h, float *h_inv, float *boxlo_bound, float *boxhi_bound, float *boxhi, float *boxlo, float *boxtilt, int triclinic);

template <typename T> void cpuSetLocalOrthBox(T *subhi, T *sublo, T *boxhi, T *boxlo, T *subhi_lamda, T *sublo_lamda, int dim)
{
    for (int i=0; i<dim; i++) {
        sublo[i] = boxlo[i] + (boxhi[i] - boxlo[i])*sublo_lamda[i];
        subhi[i] = boxlo[i] + (boxhi[i] - boxlo[i])*subhi_lamda[i];
    }    
}
template void cpuSetLocalOrthBox(double *subhi, double *sublo, double *boxhi, double *boxlo, 
        double *subhi_lamda, double *sublo_lamda, int dim);
template void cpuSetLocalOrthBox(float *subhi, float *sublo, float *boxhi, float *boxlo, 
        float *subhi_lamda, float *sublo_lamda, int dim);

template <typename T> void cpuShiftLocalOrthBox(T *subhi, T *sublo, T *boxhi, T *boxlo, T *epsilon, int *pbc, int dim)
{
    for (int i=0; i<dim; i++)
      if (pbc[i]) {
        if (fabs(sublo[i] - boxlo[i]) < epsilon[i]) sublo[i] -= epsilon[i];
        if (fabs(subhi[i] - boxhi[i]) < epsilon[i]) subhi[i] -= epsilon[i];
      }
}

template <typename T> void cpuLamda2Box(T *x, T *lambda, T *h, T *boxlo, int dim, int n)
{    
    for (int i = 0; i < n; i++) {
        int k = i*dim;
        x[k+0] = h[0]*lambda[k+0] + h[5]*lambda[k+1] + boxlo[0];
        x[k+1] = h[1]*lambda[k+1] + boxlo[1];    
        if (dim==3) {
            x[k+0] += h[4]*lambda[k+2];
            x[k+1] += h[3]*lambda[k+2];    
            x[k+2] = h[2]*lambda[k+2] + boxlo[2];
        }
    }        
}
template void cpuLamda2Box(double *x, double *lambda, double *h, double *boxlo, int dim, int n);
template void cpuLamda2Box(float *x, float *lambda, float *h, float *boxlo, int dim, int n);

template <typename T> void cpuBox2Lamda(T *lambda, T *x, T *h_inv, T *boxlo, int dim, int n)
{
  for (int i = 0; i < n; i++) {
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
  }
}
template void cpuBox2Lamda(double *lambda, double *x, double *h_inv, double *boxlo, int dim, int n);
template void cpuBox2Lamda(float *lambda, float *x, float *h_inv, float *boxlo, int dim, int n);

/* ----------------------------------------------------------------------
   convert 8 lamda corner pts of lo/hi box to box coords
   return bboxlo/hi = bounding box around 8 corner pts in box coords
------------------------------------------------------------------------- */
template <typename T> void cpuDomainBbox(T *bboxlo, T *bboxhi, T *lo_lamda, T *hi_lamda, T *boxlo, T *h, int dim)
{
  T x[3];

  bboxlo[0] = bboxlo[1] = bboxlo[2] = BIG;
  bboxhi[0] = bboxhi[1] = bboxhi[2] = -BIG;

  x[0] = lo_lamda[0]; x[1] = lo_lamda[1]; 
  if (dim==3) x[2] = lo_lamda[2]; 
  cpuLamda2Box(x, x, h, boxlo, dim, 1);
  bboxlo[0] = MIN(bboxlo[0],x[0]); bboxhi[0] = MAX(bboxhi[0],x[0]);
  bboxlo[1] = MIN(bboxlo[1],x[1]); bboxhi[1] = MAX(bboxhi[1],x[1]);
  if (dim==3)
    bboxlo[2] = MIN(bboxlo[2],x[2]); bboxhi[2] = MAX(bboxhi[2],x[2]);

  x[0] = hi_lamda[0]; x[1] = lo_lamda[1]; 
  if (dim==3) x[2] = lo_lamda[2];
  cpuLamda2Box(x, x, h, boxlo, dim, 1);
  bboxlo[0] = MIN(bboxlo[0],x[0]); bboxhi[0] = MAX(bboxhi[0],x[0]);
  bboxlo[1] = MIN(bboxlo[1],x[1]); bboxhi[1] = MAX(bboxhi[1],x[1]);
  if (dim==3)
    bboxlo[2] = MIN(bboxlo[2],x[2]); bboxhi[2] = MAX(bboxhi[2],x[2]);

  x[0] = lo_lamda[0]; x[1] = hi_lamda[1]; 
  if (dim==3) x[2] = lo_lamda[2];
  cpuLamda2Box(x, x, h, boxlo, dim, 1);
  bboxlo[0] = MIN(bboxlo[0],x[0]); bboxhi[0] = MAX(bboxhi[0],x[0]);
  bboxlo[1] = MIN(bboxlo[1],x[1]); bboxhi[1] = MAX(bboxhi[1],x[1]);
  if (dim==3)
    bboxlo[2] = MIN(bboxlo[2],x[2]); bboxhi[2] = MAX(bboxhi[2],x[2]);

  x[0] = hi_lamda[0]; x[1] = hi_lamda[1]; 
  if (dim==3) x[2] = lo_lamda[2];
  cpuLamda2Box(x, x, h, boxlo, dim, 1);
  bboxlo[0] = MIN(bboxlo[0],x[0]); bboxhi[0] = MAX(bboxhi[0],x[0]);
  bboxlo[1] = MIN(bboxlo[1],x[1]); bboxhi[1] = MAX(bboxhi[1],x[1]);
  if (dim==3)
    bboxlo[2] = MIN(bboxlo[2],x[2]); bboxhi[2] = MAX(bboxhi[2],x[2]);

  x[0] = lo_lamda[0]; x[1] = lo_lamda[1]; 
  if (dim==3) x[2] = hi_lamda[2];
  cpuLamda2Box(x, x, h, boxlo, dim, 1);
  bboxlo[0] = MIN(bboxlo[0],x[0]); bboxhi[0] = MAX(bboxhi[0],x[0]);
  bboxlo[1] = MIN(bboxlo[1],x[1]); bboxhi[1] = MAX(bboxhi[1],x[1]);
  if (dim==3)
    bboxlo[2] = MIN(bboxlo[2],x[2]); bboxhi[2] = MAX(bboxhi[2],x[2]);

  x[0] = hi_lamda[0]; x[1] = lo_lamda[1]; 
  if (dim==3) x[2] = hi_lamda[2];
  cpuLamda2Box(x, x, h, boxlo, dim, 1);
  bboxlo[0] = MIN(bboxlo[0],x[0]); bboxhi[0] = MAX(bboxhi[0],x[0]);
  bboxlo[1] = MIN(bboxlo[1],x[1]); bboxhi[1] = MAX(bboxhi[1],x[1]);
  if (dim==3)
    bboxlo[2] = MIN(bboxlo[2],x[2]); bboxhi[2] = MAX(bboxhi[2],x[2]);

  x[0] = lo_lamda[0]; x[1] = hi_lamda[1]; 
  if (dim==3) x[2] = hi_lamda[2];
  cpuLamda2Box(x, x, h, boxlo, dim, 1);
  bboxlo[0] = MIN(bboxlo[0],x[0]); bboxhi[0] = MAX(bboxhi[0],x[0]);
  bboxlo[1] = MIN(bboxlo[1],x[1]); bboxhi[1] = MAX(bboxhi[1],x[1]);
  if (dim==3)
    bboxlo[2] = MIN(bboxlo[2],x[2]); bboxhi[2] = MAX(bboxhi[2],x[2]);

  x[0] = hi_lamda[0]; x[1] = hi_lamda[1]; 
  if (dim==3) x[2] = hi_lamda[2];
  cpuLamda2Box(x, x, h, boxlo, dim, 1);
  bboxlo[0] = MIN(bboxlo[0],x[0]); bboxhi[0] = MAX(bboxhi[0],x[0]);
  bboxlo[1] = MIN(bboxlo[1],x[1]); bboxhi[1] = MAX(bboxhi[1],x[1]);
  if (dim==3)
    bboxlo[2] = MIN(bboxlo[2],x[2]); bboxhi[2] = MAX(bboxhi[2],x[2]);
}

// return shifted sub box in box coords if triclinic == 0
//                        in lambda coordinates if triclinic == 1
template <typename T> void cpuShiftedSubbox(T *ssublo, T *ssubhi, T *boxlo, T *boxhi, 
        T *boxlo_lamda, T *boxhi_lamda, T *sublo, T *subhi, 
        T *sublo_lamda, T *subhi_lamda, T *epsilon, int *pbc, int triclinic)
{
  if (triclinic == 0) {
    // shited subbox in atom coordinates
    ssublo[0] = sublo[0]; ssubhi[0] = subhi[0];
    ssublo[1] = sublo[1]; ssubhi[1] = subhi[1];
    ssublo[2] = sublo[2]; ssubhi[2] = subhi[2];
    cpuShiftLocalOrthBox(ssubhi, ssublo, boxhi, boxlo, epsilon, pbc, 3);
  } else {
    // shited subbox in lambda coordinates
    ssublo[0] = sublo_lamda[0]; ssubhi[0] = subhi_lamda[0];
    ssublo[1] = sublo_lamda[1]; ssubhi[1] = subhi_lamda[1];
    ssublo[2] = sublo_lamda[2]; ssubhi[2] = subhi_lamda[2];
    cpuShiftLocalOrthBox(ssubhi, ssublo, boxhi_lamda, boxlo_lamda, epsilon, pbc, 3);
  }
}
template void cpuShiftedSubbox(double *ssublo, double *ssubhi, double *boxlo, double *boxhi, 
        double *boxlo_lamda, double *boxhi_lamda, double *sublo, double *subhi, 
        double *sublo_lamda, double *subhi_lamda, double *epsilon, int *pbc, int triclinic);
template void cpuShiftedSubbox(float *ssublo, float *ssubhi, float *boxlo, float *boxhi, 
        float *boxlo_lamda, float *boxhi_lamda, float *sublo, float *subhi, 
        float *sublo_lamda, float *subhi_lamda, float *epsilon, int *pbc, int triclinic);

// return bounding sub box in box coords 
template <typename T> void cpuBoundingSubbox(T *bsublo, T *bsubhi, T *sublo, T *subhi, 
        T *sublo_lamda, T *subhi_lamda, T *boxlo, T *h, int triclinic)
{        
  // bounding box for a subbox in atom coordinates
  if (triclinic == 0) {
    bsublo[0] = sublo[0]; bsubhi[0] = subhi[0];
    bsublo[1] = sublo[1]; bsubhi[1] = subhi[1];
    bsublo[2] = sublo[2]; bsubhi[2] = subhi[2];
  } else // use subbox in lambda coordinates to compute the bounding box in atom coordinates
      cpuDomainBbox(bsublo, bsubhi, sublo_lamda, subhi_lamda, boxlo, h, 3);
}
template void cpuBoundingSubbox(double *bsublo, double *bsubhi, double *sublo, double *subhi, 
        double *sublo_lamda, double *subhi_lamda, double *boxlo, double *h, int triclinic);
template void cpuBoundingSubbox(float *bsublo, float *bsubhi, float *sublo, float *subhi, 
        float *sublo_lamda, float *subhi_lamda, float *boxlo, float *h, int triclinic);

template <typename T> void cpuPBCOrthBox(T *x, T *v, int *image, T *hi, T *lo, 
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
        printf("%i %g %g\n",i,h[1],x[k+1]);
        x[k+1] += h[1];
        if (vdeform) {
          v[k+0] += h_rate[5];
          v[k+1] += h_rate[1];
        }        
        image[k+1] -= 1;
        printf("%i %g %g\n",i,h[1],x[k+1]);
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

template <typename T> void cpuPBCTiltBox(T *x, T *v, int *image, T *boxlo, T *h, T *h_inv, T *h_rate, 
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
        lambda[0] += hi_lambda[0];
        if (vdeform) v[k+0] += h_rate[0];
        image[k+0] -= 1;
      }
      if (lambda[0] >= hi_lambda[0]) {
        lambda[0] -= hi_lambda[0];
        lambda[0] = MAX(lambda[0],lo_lambda[0]);
        if (vdeform) v[k+0] -= h_rate[0];
        image[k+0] += 1;
      }
    }

    if (pbc[1]) {
      if (lambda[1] < lo_lambda[1]) {
        lambda[1] += hi_lambda[1];
        if (vdeform) {
          v[k+0] += h_rate[5];
          v[k+1] += h_rate[1];
        }
        image[k+1] -= 1;
      }
      if (lambda[1] >= hi_lambda[1]) {
        lambda[1] -= hi_lambda[1];
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
        lambda[2] += hi_lambda[2];
        if (vdeform) {
          v[k+0] += h_rate[4];
          v[k+1] += h_rate[3];
          v[k+2] += h_rate[2];
        }
        image[k+2] -= 1;
      }
      if (lambda[2] >= hi_lambda[2]) {
        lambda[2] -= hi_lambda[2];
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

template <typename T> void cpuPBC(T *x, T *v, int *image, T *boxhi, T *boxlo, T *hi_lambda, T *lo_lambda,  
        T *h, T *h_inv, T *h_rate, int *pbc, int vdeform, int triclinic, int dim, int nlocal)
{
    if (triclinic==0) 
        cpuPBCOrthBox(x, v, image, boxhi, boxlo, h, h_rate, pbc, vdeform, dim, nlocal);
    else
        cpuPBCTiltBox(x, v, image, boxlo, h, h_inv, h_rate, hi_lambda, lo_lambda, pbc, vdeform, dim, nlocal);
}
template void cpuPBC(double *x, double *v, int *image, double *boxhi, double *boxlo, double *hi_lambda, double *lo_lambda,  
        double *h, double *h_inv, double *h_rate, int *pbc, int vdeform, int triclinic, int dim, int nlocal);
template void cpuPBC(float *x, float *v, int *image, float *boxhi, float *boxlo, float *hi_lambda, float *lo_lambda,  
        float *h, float *h_inv, float *h_rate, int *pbc, int vdeform, int triclinic, int dim, int nlocal); 

template <typename T> void cpuInsideOrthBox(T *inside, T *x, T *lo, T *hi, int dim, int n)
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

template <typename T> void cpuInsideTiltBox(T *inside, T *x, T *lo_lamda, T *hi_lamda, 
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

template <typename T> void cpuInsideBox(T *inside, T *x, T *boxlo, T *boxhi, T *lo_lamda, T *hi_lamda, 
        T *h_inv, int triclinic, int dim, int n)
{
    if (triclinic==0) 
        cpuInsideOrthBox(inside, x, boxlo, boxhi, dim, n);
    else
        cpuInsideTiltBox(inside, x, boxlo, lo_lamda, hi_lamda, h_inv, dim, n);
}
template void cpuInsideBox(double *inside, double *x, double *boxlo, double *boxhi, double *lo_lamda, double *hi_lamda, double *h_inv, int triclinic, int dim, int n);
template void cpuInsideBox(float *inside, float *x, float *boxlo, float *boxhi, float *lo_lamda, float *hi_lamda, float *h_inv, int triclinic, int dim, int n);

template <typename T> void cpuMinimumImage(T *dp, T *h, int *pbc, int triclinic, int dim)
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
template void cpuMinimumImage(double *dp, double *h, int *pbc, int triclinic, int dim);
template void cpuMinimumImage(float *dp, float *h, int *pbc, int triclinic, int dim);

template <typename T> void cpuMinimumImageOrthBox(T *dp, T *h, int *pbc, int dim, int n)
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

template <typename T> void cpuMinimumImageTiltBox(T *dp, T *h, int *pbc, int dim, int n)
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

template <typename T> void cpuMinimumImage(T *dp, T *h, int *pbc, int triclinic, int dim, int n)
{
    if (triclinic == 0)
        cpuMinimumImageOrthBox(dp, h, pbc, dim, n);
    else
        cpuMinimumImageTiltBox(dp, h, pbc, dim, n);
}
template void cpuMinimumImage(double *dp, double *h, int *pbc, int triclinic, int dim, int n);
template void cpuMinimumImage(float *dp, float *h, int *pbc, int triclinic, int dim, int n);

template <typename T> void cpuUnmap(T *y, T *x, T *h, int *image, int triclinic, int dim)
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
template void cpuUnmap(double *y, double *x, double *h, int *image, int triclinic, int dim);
template void cpuUnmap(float *y, float *x, float *h, int *image, int triclinic, int dim);

template <typename T> void cpuUnmapOrthBox(T *y, T *x, T *h, int *image, int dim, int n)
{
     for (int i=0; i<n; i++) {
        int k = i*dim; 
        y[k+0] = x[k+0] + image[k+0]*h[0];
        y[k+1] = x[k+1] + image[k+1]*h[1];
        if (dim==3)
            y[k+2] = x[k+2] + image[k+2]*h[2];
     }
}

template <typename T> void cpuUnmapTiltBox(T *y, T *x, T *h, int *image, int dim, int n)
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

template <typename T> void cpuUnmap(T *y, T *x, T *h, int *image, int triclinic, int dim, int n)
{
    if (triclinic == 0)
        cpuUnmapOrthBox(y, x, h, image, dim, n);
    else
        cpuUnmapTiltBox(y, x, h, image, dim, n);
}
template void cpuUnmap(double *y, double *x, double *h, int *image, int triclinic, int dim, int n);
template void cpuUnmap(float *y, float *x, float *h, int *image, int triclinic, int dim, int n);


// void cpuCreateBox()
// {    
//   for (int i=0; i<3; i++) {
//     boxlo[i] = lo[i];
//     boxhi[i] = hi[i];    
//   }
//   
//   if (triclinic==1)
//     for (int i=0; i<3; i++)  
//         boxtilt[i] = tilt[i];
//   
//   cpuSetGlobalBox(h, h_inv, boxlo_bound, boxhi_bound, boxhi, boxlo, boxtilt, triclinic);        
// }

