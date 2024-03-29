/***************************************************************************                               
                    Molecular Dynamics Potentials (MDP)
                           CESMIX-MIT Project  
 
 Contributing authors: Ngoc-Cuong Nguyen (cuongng@mit.edu, exapde@gmail.com)
 ***************************************************************************/
               
#ifndef MDP_GPUVELOCITYINTEGRATION
#define MDP_GPUVELOCITYINTEGRATION

#define IA 16807
#define IM 2147483647
#define AM (1.0/IM)
#define IQ 127773
#define IR 2836
        
double gpuRandomUniform(int *seed)
{
  int k = seed[0]/IQ;
  seed[0] = IA*(seed[0]-k*IQ) - IR*k;
  if (seed[0] < 0) seed[0] += IM;
  double ans = AM*seed[0];
  return ans;
}

/* ----------------------------------------------------------------------
   gaussian RN
------------------------------------------------------------------------- */

template <typename T> T gpuRandomGaussian(int *seed, int *save, T *second)
{
  T first,v1,v2,rsq,fac;

  if (!save[0]) {
    do {
      v1 = 2.0*gpuRandomUniform(seed)-1.0;
      v2 = 2.0*gpuRandomUniform(seed)-1.0;
      rsq = v1*v1 + v2*v2;
    } while ((rsq >= 1.0) || (rsq == 0.0));
    fac = sqrt(-2.0*log(rsq)/rsq);
    second[0] = v1*fac;
    first = v2*fac;
    save[0] = 1;
  } else {
    first = second[0];
    save[0] = 0;
  }
  return first;
}

double gpuGamdev(int ia, int *seed)
{
  int j;
  double am,e,s,v1,v2,x,y;

  if (ia < 1) return 0.0;
  if (ia < 6) {
    x=1.0;
    for (j=1; j<=ia; j++)
      x *= gpuRandomUniform(seed);

    // make certain, that -log() doesn't overflow.
    if (x < 2.2250759805e-308)
      x = 708.4;
    else
      x = -log(x);
  } else {
  restart:
    do {
      do {
        do {
          v1 = gpuRandomUniform(seed);
          v2 = 2.0*gpuRandomUniform(seed) - 1.0;
        } while (v1*v1 + v2*v2 > 1.0);

        y=v2/v1;
        am=ia-1;
        s=sqrt(2.0*am+1.0);
        x=s*y+am;
      } while (x <= 0.0);

      if (am*log(x/am)-s*y < -700 || v1<0.00001) {
        goto restart;
      }

      e=(1.0+y*y)*exp(am*log(x/am)-s*y);
    } while (gpuRandomUniform(seed) > e);
  }
  return x;
}

template <typename T> __global__ void gpuKernelNVEInitialIntegrate(T *x, T *v, T *f, T *mass, T dtf, T dtv,
        int *type, int *ilist, int dim, int inum)
{
    int ii = threadIdx.x + blockIdx.x * blockDim.x;
    while (ii < inum) {
        int i = ilist[ii]; 
        T dtfm = dtf / mass[type[i]];
        v[i*dim+0] += dtfm * f[i*dim+0];
        v[i*dim+1] += dtfm * f[i*dim+1];
        x[i*dim+0] += dtv * v[i*dim+0];
        x[i*dim+1] += dtv * v[i*dim+1];        
        if (dim==3) {
            v[i*dim+2] += dtfm * f[i*dim+2];
            x[i*dim+2] += dtv * v[i*dim+2];
        }                
        ii += blockDim.x * gridDim.x;
    }  
}
template <typename T> void gpuNVEInitialIntegrate(T *x, T *v, T *f, T *mass, 
        T dtf, T dtv, int *type, int *ilist, int dim, int inum)
{        
    int blockDim = 256;
    int gridDim = (inum + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;
    gpuKernelNVEInitialIntegrate<<<gridDim, blockDim>>>(x, v, f, mass, 
            dtf, dtv, type, ilist,  dim, inum);
}
// template void gpuNVEInitialIntegrate(double *x, double *v, double *f, double *mass,
//         double dtf, double dtv, int *type, int *ilist, int dim, int inum);
// template void gpuNVEInitialIntegrate(float *x, float *v, float *f, float *mass, 
//         float dtf, float dtv, int *type, int *ilist, int dim, int inum);

template <typename T> __global__ void gpuKernelNVEFinalIntegrate(T *x, T *v, T *f, T *mass, T dtf,
        int *type, int *ilist, int dim, int inum)
{
  // update v of atoms in group  
    int ii = threadIdx.x + blockIdx.x * blockDim.x;          
    while (ii < inum) {
        int i = ilist[ii]; 
        //T dtf = dtarray[1];            
        T dtfm = dtf / mass[type[i]];
        v[i*dim+0] += dtfm * f[i*dim+0];
        v[i*dim+1] += dtfm * f[i*dim+1];
        if (dim==3) v[i*dim+2] += dtfm * f[i*dim+2];
        ii += blockDim.x * gridDim.x;
    }  
}
template <typename T> void gpuNVEFinalIntegrate(T *x, T *v, T *f, T *mass, 
        T dtf, int *type, int *ilist, int dim, int inum)
{        
    int blockDim = 256;
    int gridDim = (inum + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;
    gpuKernelNVEFinalIntegrate<<<gridDim, blockDim>>>(x, v, f, mass, 
            dtf,  type, ilist,  dim, inum);
}

template <typename T> __global__ void gpuKernelNVELimitInitialIntegrate(T *x, T *v, T *f, T *mass, T dtf, T dtv,
        T vlimitsq, int *type, int *ilist, int dim, int inum)
{
    int ii = threadIdx.x + blockIdx.x * blockDim.x;   
    while (ii < inum) {
        int i = ilist[ii]; 
        int j;
        T dtfm = dtf / mass[type[i]];
        T vsq = 0; 
        for (j=0; j<dim; j++) {
            v[i*dim+j] += dtfm * f[i*dim+j];        
            vsq += v[i*dim+j]*v[i*dim+j];
        }

        if (vsq > vlimitsq) {
          T scale = sqrt(vlimitsq/vsq);
          for (j=0; j<dim; j++) 
            v[i*dim+j] *= scale;        
        }

        for (j=0; j<dim; j++) 
            x[i*dim+j] += dtv * v[i*dim+j];        
        ii += blockDim.x * gridDim.x;
    }      
}
template <typename T> void gpuNVELimitInitialIntegrate(T *x, T *v, T *f, T *mass, 
        T dtf, T dtv, T vlimitsq, int *type, int *ilist, int dim, int inum)
{        
    int blockDim = 256;
    int gridDim = (inum + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;
    gpuKernelNVELimitInitialIntegrate<<<gridDim, blockDim>>>(x, v, f, mass, 
            dtf, dtv, vlimitsq, type, ilist,  dim, inum);
}

template <typename T> __global__ void gpuKernelNVELimitFinalIntegrate(T *x, T *v, T *f, T *mass, T dtf,
        T vlimitsq, int *type, int *ilist, int dim, int inum)
{ 
    int ii = threadIdx.x + blockIdx.x * blockDim.x;   
    while (ii < inum) {
        int i = ilist[ii]; 
        int j;
        T dtfm = dtf / mass[type[i]];
        T vsq = 0; 
        for (j=0; j<dim; j++) {
            v[i*dim+j] += dtfm * f[i*dim+j];        
            vsq += v[i*dim+j]*v[i*dim+j];
        }

        if (vsq > vlimitsq) {          
          T scale = sqrt(vlimitsq/vsq);
          for (j=0; j<dim; j++) 
            v[i*dim+j] *= scale;        
        }      
        ii += blockDim.x * gridDim.x;
    }  
}
template <typename T> void gpuNVELimitFinalIntegrate(T *x, T *v, T *f, T *mass, 
        T dtf, T vlimitsq, int *type, int *ilist, int dim, int inum)
{        
    int blockDim = 256;
    int gridDim = (inum + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;
    gpuKernelNVELimitFinalIntegrate<<<gridDim, blockDim>>>(x, v, f, mass, 
            dtf, vlimitsq, type, ilist,  dim, inum);
}

template <typename T> __global__ void gpuKernelSetVelocityInitialIntegrate(T *x, T *v, T *f, T *mass, T *fparam,
        T dtf, T dtv, int *type, int *ilist, int *iparam, int dim, int inum)
{    
    int ii = threadIdx.x + blockIdx.x * blockDim.x;   
    while (ii < inum) {
        int i = ilist[ii]; 
        T dtfm = dtf / mass[type[i]];        
        for (int j=0; j<dim; j++) {
            v[i*dim+j] = (iparam[j]) ? fparam[j] : (v[i*dim+j] + dtfm * f[i*dim+j]);                
            x[i*dim+j] += dtv * v[i*dim+j];
        }
        ii += blockDim.x * gridDim.x;
    }  
}
template <typename T> void gpuSetVelocityInitialIntegrate(T *x, T *v, T *f, T *mass, T *fparam,
        T dtf, T dtv, int *type, int *ilist, int *iparam, int dim, int inum)
{        
    int blockDim = 256;
    int gridDim = (inum + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;
    gpuKernelSetVelocityInitialIntegrate<<<gridDim, blockDim>>>(x, v, f, mass, fparam, 
            dtf, dtv, type, ilist, iparam, dim, inum);
}
template void gpuSetVelocityInitialIntegrate(double *x, double *v, double *f, double *mass, double *fparam,
        double dtf, double dtv, int *type, int *ilist, int *iparam, int dim, int inum);
template void gpuSetVelocityInitialIntegrate(float *x, float *v, float *f, float *mass, float *fparam,
        float dtf, float dtv, int *type, int *ilist, int *iparam, int dim, int inum);

template <typename T> __global__ void gpuKernelSetVelocityFinalIntegrate(T *x, T *v, T *f, T *mass, 
        T dtf, int *type, int *ilist, int *iparam, int dim, int inum)
{
  // update v of atoms in group  
    int ii = threadIdx.x + blockIdx.x * blockDim.x;            
    while (ii < inum) {
        int i = ilist[ii]; 
        T dtfm = dtf / mass[type[i]];
        for (int j=0; j<dim; j++) 
            v[i*dim+j] = (iparam[j]) ? v[i*dim+j] : (v[i*dim+j] + dtfm * f[i*dim+j]);                         
        ii += blockDim.x * gridDim.x;
    }  
}
template <typename T> void gpuSetVelocityFinalIntegrate(T *x, T *v, T *f, T *mass, 
        T dtf, int *type, int *ilist, int *iparam, int dim, int inum)
{        
    int blockDim = 256;
    int gridDim = (inum + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;
    gpuKernelSetVelocityFinalIntegrate<<<gridDim, blockDim>>>(x, v, f, mass, 
            dtf, type, ilist, iparam, dim, inum);
}
template void gpuSetVelocityFinalIntegrate(double *x, double *v, double *f, double *mass, 
        double dtf, int *type, int *ilist, int *iparam, int dim, int inum);
template void gpuSetVelocityFinalIntegrate(float *x, float *v, float *f, float *mass, 
        float dtf, int *type, int *ilist, int *iparam, int dim, int inum);

template <typename T> __global__ void gpuKernelVelocityFactor(T *v, T factor_eta, int *ilist, 
        int biasflag, int dim, int inum)
{
    int ii = threadIdx.x + blockIdx.x * blockDim.x;     
    while (ii < inum) {
        int i = ilist[ii]; 
        for (int j=0; j<dim; j++) 
            v[i*dim+j] *= factor_eta;        
        ii += blockDim.x * gridDim.x;
    }    
}
template <typename T> void gpuVelocityFactor(T *v, T factor, int *ilist, 
         int biasflag, int dim, int inum)
{        
    int blockDim = 256;
    int gridDim = (inum + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;
    gpuKernelVelocityFactor<<<gridDim, blockDim>>>(v, factor, ilist, biasflag, dim, inum);
}
template void gpuVelocityFactor(double *v, double factor, int *ilist, int biasflag, int dim, int inum);
template void gpuVelocityFactor(float *v, float factor, int *ilist, int biasflag, int dim, int inum);

template <typename T> void gpuNoseHooverThermostat(T *v, T *dtarray, T *tarray, T *eta_mass, T *eta, 
        T *eta_dot, T *eta_dotdot, int *ilist, int eta_mass_flag, int biasflag, int mtchain, 
        int nc_tchain, int dim, int inum)
{
  int ich;
  T dt = dtarray[0];
  T dtf = dtarray[1];
  T dtv = dtarray[2];
  T dthalf = 0.5 * dt;
  T dt4 = 0.25 * dt;
  T dt8 = 0.125 * dt;
  
  //T t_start = tarray[0];  
  //T t_stop = tarray[1];     
  //T t_freq = tarray[2];
  //T t_target = tarray[3]; 
  T t_current = tarray[4];   
  //T tdof = tarray[5];
  //T boltz = tarray[6];
  //T drag = tarray[7];
  //T mvv2e = tarray[8];    
  //T tfactor = mvv2e / (tdof * boltz);
    
  T t_freq = tarray[2];
  T t_target = tarray[3]; 
  //T t_current = tarray[4];   
  T tdof = tarray[5];
  T boltz = tarray[6];
  T drag = tarray[7];
  T tdrag_factor = 1.0 - (dt * t_freq * drag / nc_tchain);
  
  T expfac;
  T kecurrent = tdof * boltz * tarray[4];
  T ke_target = tdof * boltz * t_target;
  
  // Update masses, to preserve initial freq, if flag set

  if (eta_mass_flag) {
    eta_mass[0] = tdof * boltz * t_target / (t_freq*t_freq);
    for (int ich = 1; ich < mtchain; ich++)
      eta_mass[ich] = boltz * t_target / (t_freq*t_freq);
  }

  if (eta_mass[0] > 0.0)
    eta_dotdot[0] = (kecurrent - ke_target)/eta_mass[0];
  else eta_dotdot[0] = 0.0;

  T ncfac = 1.0/nc_tchain;
  for (int iloop = 0; iloop < nc_tchain; iloop++) {

    for (ich = mtchain-1; ich > 0; ich--) {
      expfac = exp(-ncfac*dt8*eta_dot[ich+1]);
      eta_dot[ich] *= expfac;
      eta_dot[ich] += eta_dotdot[ich] * ncfac*dt4;
      eta_dot[ich] *= tdrag_factor;
      eta_dot[ich] *= expfac;
    }

    expfac = exp(-ncfac*dt8*eta_dot[1]);
    eta_dot[0] *= expfac;
    eta_dot[0] += eta_dotdot[0] * ncfac*dt4;
    eta_dot[0] *= tdrag_factor;
    eta_dot[0] *= expfac;

    T factor_eta = exp(-ncfac*dthalf*eta_dot[0]);
    gpuVelocityFactor(v, factor_eta, ilist, biasflag, dim, inum);

    // rescale temperature due to velocity scaling
    // should not be necessary to explicitly recompute the temperature

    tarray[4] *= factor_eta*factor_eta;
    kecurrent = tdof * boltz * tarray[4];

    if (eta_mass[0] > 0.0)
      eta_dotdot[0] = (kecurrent - ke_target)/eta_mass[0];
    else eta_dotdot[0] = 0.0;

    for (ich = 0; ich < mtchain; ich++)
      eta[ich] += ncfac*dthalf*eta_dot[ich];

    eta_dot[0] *= expfac;
    eta_dot[0] += eta_dotdot[0] * ncfac*dt4;
    eta_dot[0] *= expfac;

    for (ich = 1; ich < mtchain; ich++) {
      expfac = exp(-ncfac*dt8*eta_dot[ich+1]);
      eta_dot[ich] *= expfac;
      eta_dotdot[ich] = (eta_mass[ich-1]*eta_dot[ich-1]*eta_dot[ich-1]
                         - boltz * t_target)/eta_mass[ich];
      eta_dot[ich] += eta_dotdot[ich] * ncfac*dt4;
      eta_dot[ich] *= expfac;
    }
  }
}

template <typename T> void gpuNVTInitialIntegrate(T *x, T *v, T *f, T *mass, T *dtarray, T *tarray,
        T *eta_mass, T *eta, T *eta_dot, T *eta_dotdot, T *ke, T *tmp, int *type, int *ilist, 
        int eta_mass_flag, int biasflag, int mtchain, int nc_tchain, int dim, int inum)
{  
  T dtf = dtarray[1];
  T dtv = dtarray[2];
  T beginstep = dtarray[3];
  T endstep = dtarray[4];
  T ntimestep = dtarray[5];      
  
  T t_start = tarray[0];  
  T t_stop = tarray[1];   
  T delta = (ntimestep - beginstep)/(endstep - beginstep);  
  tarray[3] = t_start + delta * (t_stop - t_start); // t_target
    
  T tdof = tarray[5];
  T boltz = tarray[6];
  T mvv2e = tarray[8];    
  T tfactor = mvv2e / (tdof * boltz);    
  T t_current = gpuComputeTempScalar(ke, v, tmp, mass, tfactor, type, ilist, dim, inum);
  tarray[4] = t_current;   
        
  gpuNoseHooverThermostat(v, dtarray, tarray, eta_mass, eta, eta_dot, eta_dotdot, 
           ilist, eta_mass_flag, biasflag, mtchain, nc_tchain, dim, inum);

  gpuNVEInitialIntegrate(x, v, f, mass, dtf, dtv, type, ilist, dim, inum);  
}


template <typename T> void gpuNVTFinalIntegrate(T *x, T *v, T *f, T *mass, T *dtarray, T *tarray,
        T *eta_mass, T *eta, T *eta_dot, T *eta_dotdot, T *ke, T *tmp, int *type, int *ilist, 
        int eta_mass_flag, int biasflag, int mtchain, int nc_tchain, int dim, int inum)
{
  T dtf = dtarray[1];
  T dtv = dtarray[2];    
  gpuNVEFinalIntegrate(x, v, f, mass, dtf, type, ilist, dim, inum);    

  T tdof = tarray[5];
  T boltz = tarray[6];
  T mvv2e = tarray[8];    
  T tfactor = mvv2e / (tdof * boltz);    
  T t_current = gpuComputeTempScalar(ke, v, tmp, mass, tfactor, type, ilist, dim, inum);        
  tarray[4] = t_current;   
    
  gpuNoseHooverThermostat(v, dtarray, tarray, eta_mass, eta, eta_dot, eta_dotdot, 
          ilist, eta_mass_flag, biasflag, mtchain, nc_tchain, dim, inum);
}

template <typename T> void gpuInitialIntegrate(T *x, T *v, T *f, T *mass, T *dtarray, T *tarray,
        T *eta_mass, T *eta, T *eta_dot, T *eta_dotdot, T *ke, T *tmp, T vlimitsq, int *type, int *ilist,
        int eta_mass_flag, int biasflag, int mtchain, int nc_tchain, int mode, int dim, int inum)
{
    T dtf = dtarray[1];
    T dtv = dtarray[2];        
    if (mode==0) {
        gpuNVEInitialIntegrate(x, v, f, mass, dtf, dtv,
            type, ilist, dim, inum);
    }
    else if (mode==1) {
        gpuNVELimitInitialIntegrate(x, v, f, mass, dtf, dtv,
            vlimitsq, type, ilist, dim, inum);
    }
    else if (mode==2) {
        gpuNVTInitialIntegrate(x, v, f, mass, dtarray, tarray,
            eta_mass, eta, eta_dot, eta_dotdot, ke, tmp, type, ilist, 
            eta_mass_flag, biasflag, mtchain, nc_tchain, dim, inum);        
//         gpuNVTInitialIntegrate(x, v, f, mass, dtarray, tarray,
//             eta_mass, eta, eta_dot, eta_dotdot, ke, tmp, type, ilist, 
//             eta_mass_flag, biasflag, mtchain, nc_tchain, dim, inum);
    }
//     else if (mode==3) {
//         gpuMoveInitialIntegrate(x, v, f, mass, fparam, dtf, dtv, type, ilist, iparam, dim, inum);
//     }
}
template void gpuInitialIntegrate(double *x, double *v, double *f, double *mass, 
        double *dtarray, double *tarray, double *eta_mass, double *eta, double *eta_dot, 
        double *eta_dotdot, double *ke, double *tmp, double vlimitsq, int *type, int *ilist, 
        int eta_mass_flag, int biasflag,  int mtchain, int nc_tchain, int mode, int dim, int inum);
template void gpuInitialIntegrate(float *x, float *v, float *f, float *mass, 
        float *dtarray, float *tarray, float *eta_mass, float *eta, float *eta_dot, 
        float *eta_dotdot, float *ke, float *tmp, float vlimitsq, int *type, int *ilist, 
        int eta_mass_flag, int biasflag, int mtchain, int nc_tchain, int mode, int dim, int inum);

template <typename T> void gpuFinalIntegrate(T *x, T *v, T *f, T *mass, T *dtarray, T *tarray,
        T *eta_mass, T *eta, T *eta_dot, T *eta_dotdot, T *ke, T *tmp, T vlimitsq, int *type, int *ilist, 
        int eta_mass_flag, int biasflag, int mtchain, int nc_tchain, int mode, int dim, int inum)
{    
    T dtf = dtarray[1];
    if (mode==0) {
        gpuNVEFinalIntegrate(x, v, f, mass, dtf, type, ilist, dim, inum);
    }
    else if (mode==1) {
        gpuNVELimitFinalIntegrate(x, v, f, mass, dtf, vlimitsq, type, ilist, dim, inum);        
    }
    else if (mode==2) {
        gpuNVTFinalIntegrate(x, v, f, mass, dtarray, tarray, eta_mass, eta, eta_dot, eta_dotdot, 
                ke, tmp, type, ilist, eta_mass_flag, biasflag, mtchain, nc_tchain, dim, inum);
    }
//     else if (mode==3) {
//         gpuMoveFinalIntegrate(x, v, f, mass, fparam, dtf, type, ilist, iparam, dim, inum);
//     }    
}
template void gpuFinalIntegrate(double *x, double *v, double *f, double *mass, 
        double *dtarray, double *tarray, double *eta_mass, double *eta, double *eta_dot, 
        double *eta_dotdot, double *ke, double *tmp, double vlimitsq, int *type, int *ilist, 
        int eta_mass_flag, int biasflag, int mtchain, int nc_tchain, int mode, int dim, int inum);
template void gpuFinalIntegrate(float *x, float *v, float *f, float *mass, 
        float *dtarray, float *tarray, float *eta_mass, float *eta, float *eta_dot, 
        float *eta_dotdot, float *ke, float *tmp, float vlimitsq, int *type, int *ilist, 
        int eta_mass_flag, int biasflag, int mtchain, int nc_tchain, int mode, int dim, int inum);

template <typename T> T gpuBerendsenThermostat(T *v, T *mass, T *dtarray, T *tarray, T *ke, T *tmp, 
        T energy, int *type, int *ilist, int biasflag, int dim, int inum)
{
  T tdof = tarray[5];
  T boltz = tarray[6];
  T mvv2e = tarray[8];      
  T tfactor = mvv2e / (tdof * boltz);        
  T t_current = gpuComputeTempScalar(ke, v, tmp, mass, tfactor, type, ilist, dim, inum);
  tarray[4] = t_current;         
  
  if (tdof < 1) return energy;
  //if (t_current == 0.0)
  //  error("Computed temperature for fix temp/berendsen cannot be 0.0");

  T dt = dtarray[0];
  T beginstep = dtarray[3];
  T endstep = dtarray[4];
  T ntimestep = dtarray[5];        
  T t_start = tarray[0];  
  T t_stop = tarray[1];   
  T t_freq = tarray[2];
  T delta = (ntimestep - beginstep)/(endstep - beginstep);  
  T t_target  = t_start + delta * (t_stop - t_start); // 
  tarray[3] = t_target;
    
  T lamda = sqrt(1.0 + dt*t_freq*(t_target/t_current - 1.0));
  T efactor = 0.5 * boltz * tdof;
  energy += t_current * (1.0-lamda*lamda) * efactor;

  gpuVelocityFactor(v, lamda, ilist, biasflag, dim, inum);
    
  return energy;
}

template <typename T> T gpuRescaleThermostat(T *v, T *mass, T *dtarray, T *tarray, T *ke, T *tmp, 
        T energy, int *type, int *ilist, int biasflag, int dim, int inum)
{
  T tdof = tarray[5];
  T boltz = tarray[6];
  T mvv2e = tarray[8];      
  T t_window = tarray[9];    
  T fraction = tarray[10];    
  T tfactor = mvv2e / (tdof * boltz);        
  T t_current = gpuComputeTempScalar(ke, v, tmp, mass, tfactor, type, ilist, dim, inum);
  tarray[4] = t_current;         
  
  if (tdof < 1) return energy;
  //if (t_current == 0.0)
  //  error("Computed temperature for fix temp/berendsen cannot be 0.0");

  T dt = dtarray[0];
  T beginstep = dtarray[3];
  T endstep = dtarray[4];
  T ntimestep = dtarray[5];        
  T t_start = tarray[0];  
  T t_stop = tarray[1];   
  T t_freq = tarray[2];
  T delta = (ntimestep - beginstep)/(endstep - beginstep);  
  T t_target  = t_start + delta * (t_stop - t_start); // 
  tarray[3] = t_target;
    
  T lamda = sqrt(1.0 + dt*t_freq*(t_target/t_current - 1.0));
  T efactor = 0.5 * boltz * tdof;
  energy += t_current * (1.0-lamda*lamda) * efactor;

  if (fabs(t_current-t_target) > t_window) {
    t_target = t_current - fraction*(t_current-t_target);
    T factor = sqrt(t_target/t_current);
    T efactor = 0.5 * boltz * tdof;
    energy += (t_current-t_target) * efactor;      
    gpuVelocityFactor(v, factor, ilist, biasflag, dim, inum);    
  }  
  
  return energy;
}

template <typename T> T gpuCsvrThermostat(T *v, T *mass, T *dtarray, T *tarray, T *ke, T *tmp, 
        T *second, T energy, int *type, int *ilist, int *seed, int *save, int biasflag, int dim, int inum)
{
  T tdof = tarray[5];
  T boltz = tarray[6];
  T mvv2e = tarray[8];      
  T tfactor = mvv2e / (tdof * boltz);        
  T t_current = gpuComputeTempScalar(ke, v, tmp, mass, tfactor, type, ilist, dim, inum);
  tarray[4] = t_current;         
  
  if (tdof < 1) return energy;
  //if (t_current == 0.0)
  //  error("Computed temperature for fix temp/berendsen cannot be 0.0");

  T dt = dtarray[0];
  T beginstep = dtarray[3];
  T endstep = dtarray[4];
  T ntimestep = dtarray[5];        
  T t_start = tarray[0];  
  T t_stop = tarray[1];   
  T t_freq = tarray[2];
  T delta = (ntimestep - beginstep)/(endstep - beginstep);  
  T t_target  = t_start + delta * (t_stop - t_start); // 
  tarray[3] = t_target;
    
  //T lamda = sqrt(1.0 + dt*t_freq*(t_target/t_current - 1.0));
  T efactor = 0.5 * boltz * tdof;
  T ekin_old = t_current * efactor;
  T ekin_new = t_target * efactor;
  
  T c1 = exp(-dt*t_freq);
  T c2 = (1.0-c1)*ekin_new/ekin_old/tdof;
  T r1 = gpuRandomGaussian(seed, save, second);
  //T r2 = sumnoises(tdof - 1);
  T r2;
  int nn = (int) (tdof-1);
  if (nn == 0) {
    r2 = 0.0;
  } else if (nn == 1) {
    T rr = gpuRandomGaussian(seed, save, second);
    r2 = rr*rr;
  } else if (nn % 2 == 0) {
    r2 = 2.0 * gpuGamdev(nn / 2, seed);
  } else {
    T rr = gpuRandomGaussian(seed, save, second);
    r2 =  2.0 * gpuGamdev((nn-1) / 2, seed) + rr*rr;
  }
  
  T scale = c1 + c2*(r1*r1+r2) + 2.0*r1*sqrt(c1*c2);
  T lamda = sqrt(scale);
  
  energy += t_current * (1.0-lamda*lamda) * efactor;

  gpuVelocityFactor(v, lamda, ilist, biasflag, dim, inum);
  
  return energy;
}

template <typename T> T gpuVelocityRescalingThermostat(T *v, T *mass, T *dtarray, T *tarray, 
        T *ke, T *tmp, T *second, T energy, int *type, int *ilist, int *seed, int *save, 
        int biasflag, int mode, int dim, int inum)
{
    if (mode==0) {
        return gpuBerendsenThermostat(v, mass, dtarray, tarray, ke, tmp, energy, type, ilist, biasflag, dim, inum);
    }
    else if (mode==1) {
        return gpuRescaleThermostat(v, mass, dtarray, tarray, ke, tmp, energy, type, ilist, biasflag, dim, inum);
    }
    else if (mode==2) {
        return gpuCsvrThermostat(v, mass, dtarray, tarray, ke, tmp, second, energy, type, ilist, 
                seed, save, biasflag, dim, inum);
    }
    return 0.0;
}
template double gpuVelocityRescalingThermostat(double *v, double *mass, double *dtarray, double *tarray,
        double *ke, double *tmp, double *second, double energy, int *type, int *ilist, int *seed, 
        int *save, int biasflag, int mode, int dim, int inum);
template float gpuVelocityRescalingThermostat(float *v, float *mass, float *dtarray, float *tarray, 
        float *ke, float *tmp, float *second, float energy, int *type, int *ilist, int *seed, 
        int *save, int biasflag, int mode, int dim, int inum);


#endif

