/***************************************************************************                               
                    Molecular Dynamics Potentials (MDP)
                           CESMIX-MIT Project  
 
 Contributing authors: Ngoc-Cuong Nguyen (cuongng@mit.edu, exapde@gmail.com)
 ***************************************************************************/
               
#ifndef MDP_CPUVELOCITYINTEGRATION
#define MDP_CPUVELOCITYINTEGRATION

#define WARMUP 100

template <typename T> void cpuVelocityZeroMomentum(T *v, T *vcm, int dim, int nlocal)
{  
  for (int i = 0; i < nlocal; i++) {
      v[i*dim+0] -= vcm[0];
      v[i*dim+1] -= vcm[1];
      if (dim==3) v[i*dim+2] -= vcm[2];
  }
}
template void cpuVelocityZeroMomentum(double *v, double *vcm, int dim, int nlocal);
template void cpuVelocityZeroMomentum(float *v, float *vcm, int dim, int nlocal);

template <typename T> void cpuVelocityZeroRotation(T *x, T *v, T *box, T *xcm, T *omega, 
        int *image, int triclinic, int dim, int nlocal)
{  
  T unwrap[3];  
  if (dim==3) {
      for (int i = 0; i < nlocal; i++) {        
          cpuUnmap(unwrap, &x[i*dim], box, &image[i*dim], triclinic, dim);        
          T dx = unwrap[0] - xcm[0];
          T dy = unwrap[1] - xcm[1];
          T dz = unwrap[2] - xcm[2];
          v[i*dim+0] -= omega[1]*dz - omega[2]*dy;
          v[i*dim+1] -= omega[2]*dx - omega[0]*dz;
          v[i*dim+2] -= omega[0]*dy - omega[1]*dx;
     }
  }
  else {
      for (int i = 0; i < nlocal; i++) {  
          cpuUnmap(unwrap, &x[i*dim], box, &image[i*dim], triclinic, dim);        
          T dx = unwrap[0] - xcm[0];
          T dy = unwrap[1] - xcm[1];
          v[i*dim+0] -= -omega[2]*dy;
          v[i*dim+1] -=  omega[2]*dx;
      }      
  }
}
template void cpuVelocityZeroRotation(double *x, double *v, double *box, double *xcm, double *omega, 
        int *image, int triclinic, int dim, int nlocal);
template void cpuVelocityZeroRotation(float *x, float *v, float *box, float *xcm, float *omega, 
        int *image, int triclinic, int dim, int nlocal);

template <typename T> void cpuVelocityCreate(T *x, T *v, T *mass, T *second, 
         int *seed, int *save, int *map, int *type, int seed0, int sum_flag, int dist_flag, int loop_flag, 
        int dim, int mpiRank, int nlocal, int natoms)
{
  int i;
  
  if (sum_flag==0) {
      for (i = 0; i < dim*nlocal; i++)
          v[i] = 0.0;
  }

  int m;
  T vx,vy,vz,factor;

  if (loop_flag == 0) { // ALL

    // loop over all atoms in system
    // generate RNGs for all atoms, only assign to ones I own
    // use either per-type mass or per-atom rmass
    save[0] = 0;
    for (i = 1; i <= natoms; i++) {
      if (dist_flag == 0) {
        vx = cpuRandomUniform(seed) - 0.5;
        vy = cpuRandomUniform(seed) - 0.5;
        vz = cpuRandomUniform(seed) - 0.5;
      } else {
        vx = cpuRandomGaussian(seed, save, second);
        vy = cpuRandomGaussian(seed, save, second);
        vz = cpuRandomGaussian(seed, save, second);
      }
      m = map[i]; // map global ID to local ID
      if (m >= 0 && m < nlocal) {        
          factor = 1.0/sqrt(mass[type[m]]);
          v[m*dim+0] += vx * factor;
          v[m*dim+1] += vy * factor;
          if (dim == 3) 
              v[m*dim+2] += vz * factor;        
      }
    }
  } else if (loop_flag == 1) { // LOCAL
    save[0] = 0;
    seed[0] = seed[0] + mpiRank;
    for (i = 0; i < WARMUP; i++) cpuRandomUniform(seed);

    for (i = 0; i < nlocal; i++) {
        if (dist_flag == 0) {
          vx = cpuRandomUniform(seed) - 0.5;
          vy = cpuRandomUniform(seed) - 0.5;
          vz = cpuRandomUniform(seed) - 0.5;
        } else {
          vx = cpuRandomGaussian(seed, save, second);
          vy = cpuRandomGaussian(seed, save, second);
          vz = cpuRandomGaussian(seed, save, second);
        }
        factor = 1.0/sqrt(mass[type[i]]);
        v[i*dim+0] += vx * factor;
        v[i*dim+1] += vy * factor;
        if (dim == 3) v[i*dim+2] += vz * factor;      
    }

  } else if (loop_flag == 2) {      // GEOM
    save[0] = 0;
    seed[0] = 1;    
    for (i = 0; i < nlocal; i++) {
        cpuRandomResetSeed(seed,save,seed0,&x[i*dim]);        
        if (dist_flag == 0) {
          vx = cpuRandomUniform(seed) - 0.5;
          vy = cpuRandomUniform(seed) - 0.5;
          vz = cpuRandomUniform(seed) - 0.5;
        } else {
          vx = cpuRandomGaussian(seed, save, second);
          vy = cpuRandomGaussian(seed, save, second);
          vz = cpuRandomGaussian(seed, save, second);
        }

        factor = 1.0/sqrt(mass[type[i]]);
        v[i*dim+0] += vx * factor;
        v[i*dim+1] += vy * factor;
        if (dim == 3) v[i*dim+2] += vz * factor;      
    }
  }
}
template void cpuVelocityCreate(double *x, double *v, double *mass, double *second, 
         int *seed, int *save, int *map, int *type, int seed0, int sum_flag, int dist_flag, int loop_flag, 
        int dim, int mpiRank, int nlocal, int natoms);
template void cpuVelocityCreate(float *x, float *v, float *mass, float *second, 
         int *seed, int *save, int *map, int *type, int seed0, int sum_flag, int dist_flag, int loop_flag, 
        int dim, int mpiRank, int nlocal, int natoms);

template <typename T> void cpuVelocityCreate(T *x, T *v, T *mass, T *second, T *omega, 
        T *box, T *xcm, T *vcm, T t_desired, T t_current, int *seed, int *save, int *map, int *image, 
        int *type, int sum_flag, int dist_flag, int loop_flag, int rotation_flag, 
        int momentum_flag, int triclinic, int dim, int mpiRank, int nlocal, int natoms)
{
  int i;
  
  if (sum_flag==0) {
      for (i = 0; i < dim*nlocal; i++)
          v[i] = 0.0;
  }

  int m;
  T vx,vy,vz,factor;

  if (loop_flag == 0) { // ALL

    // loop over all atoms in system
    // generate RNGs for all atoms, only assign to ones I own
    // use either per-type mass or per-atom rmass
    save[0] = 0;
    for (i = 1; i <= natoms; i++) {
      if (dist_flag == 0) {
        vx = cpuRandomUniform(seed) - 0.5;
        vy = cpuRandomUniform(seed) - 0.5;
        vz = cpuRandomUniform(seed) - 0.5;
      } else {
        vx = cpuRandomGaussian(seed, save, second);
        vy = cpuRandomGaussian(seed, save, second);
        vz = cpuRandomGaussian(seed, save, second);
      }
      m = map[i]; // map global ID to local ID
      if (m >= 0 && m < nlocal) {        
          factor = 1.0/sqrt(mass[type[m]]);
          v[m*dim+0] += vx * factor;
          v[m*dim+1] += vy * factor;
          if (dim == 3) 
              v[m*dim+2] += vz * factor;        
      }
    }
  } else if (loop_flag == 1) { // LOCAL
    save[0] = 0;
    seed[0] = seed[0] + mpiRank;
    for (i = 0; i < WARMUP; i++) cpuRandomUniform(seed);

    for (i = 0; i < nlocal; i++) {
        if (dist_flag == 0) {
          vx = cpuRandomUniform(seed) - 0.5;
          vy = cpuRandomUniform(seed) - 0.5;
          vz = cpuRandomUniform(seed) - 0.5;
        } else {
          vx = cpuRandomGaussian(seed, save, second);
          vy = cpuRandomGaussian(seed, save, second);
          vz = cpuRandomGaussian(seed, save, second);
        }
        factor = 1.0/sqrt(mass[type[i]]);
        v[i*dim+0] += vx * factor;
        v[i*dim+1] += vy * factor;
        if (dim == 3) v[i*dim+2] += vz * factor;      
    }

  } else if (loop_flag == 2) {      // GEOM
    save[0] = 0;
    seed[0] = 1;    
    for (i = 0; i < nlocal; i++) {
        cpuRandomResetSeed(seed,save,seed[0],&x[i*dim]);
        if (dist_flag == 0) {
          vx = cpuRandomUniform(seed) - 0.5;
          vy = cpuRandomUniform(seed) - 0.5;
          vz = cpuRandomUniform(seed) - 0.5;
        } else {
          vx = cpuRandomGaussian(seed, save, second);
          vy = cpuRandomGaussian(seed, save, second);
          vz = cpuRandomGaussian(seed, save, second);
        }

        factor = 1.0/sqrt(mass[type[i]]);
        v[i*dim+0] += vx * factor;
        v[i*dim+1] += vy * factor;
        if (dim == 3) v[i*dim+2] += vz * factor;      
    }
  }

  // apply momentum and rotation zeroing

  if (momentum_flag) cpuVelocityZeroMomentum(v, vcm, dim, nlocal);
  if (rotation_flag) cpuVelocityZeroRotation(x, v, box, xcm, omega, image, 
          triclinic, dim, nlocal);
  
  // scale temp to desired value
  //cpuVelocityScale(v, t_current, t_desired, mask, dim, nlocal);          
  cpuArrayMultiplyScalar(v, sqrt(t_desired/t_current), nlocal*dim);
}

template <typename T> void cpuVelocitySet(T *v, T *vext, int *vdim,
        int sum_flag, int dim, int nlocal)
{
    for (int i = 0; i < nlocal; i++) {
        if (sum_flag == 0) {
          if (vdim[0]) v[i*dim+0] = vext[0];
          if (vdim[1]) v[i*dim+1] = vext[1];
          if (vdim[2]) v[i*dim+2] = vext[2];
        } else {
          if (vdim[0]) v[i*dim+0] += vext[0];
          if (vdim[1]) v[i*dim+1] += vext[1];
          if (vdim[2]) v[i*dim+2] += vext[2];
        }      
    }
}
template void cpuVelocitySet(double *v, double *vext, int *vdim, int sum_flag, int dim, int nlocal);
template void cpuVelocitySet(float *v, float *vext, int *vdim, int sum_flag, int dim, int nlocal);

template <typename T> void cpuVelocityRamp(T *x, T *v, T *v_lo, T *v_hi, T *coord_lo, T *coord_hi,
        int *coord_dim, int *v_dim, int sum_flag, int dim, int nlocal)
{
  for (int i = 0; i < nlocal; i++)
      for (int j = 0; j < dim; j++) {  
          T fraction = (x[i*dim+coord_dim[j]] - coord_lo[j]) / (coord_hi[j] - coord_lo[j]);
          fraction = (fraction > 0.0) ? fraction : 0.0;
          fraction = (fraction < 1.0) ? fraction : 1.0;
          T vramp = v_lo[j] + fraction*(v_hi[j] - v_lo[j]);
          if (sum_flag) v[i*dim+v_dim[j]] += vramp;
          else v[i*dim+v_dim[j]] = vramp;
      }    
}
template void cpuVelocityRamp(double *x, double *v, double *v_lo, double *v_hi, double *coord_lo, double *coord_hi,
        int *coord_dim, int *v_dim, int sum_flag, int dim, int nlocal);
template void cpuVelocityRamp(float *x, float *v, float *v_lo, float *v_hi, float *coord_lo, float *coord_hi,
        int *coord_dim, int *v_dim, int sum_flag, int dim, int nlocal);

// template <typename T> void cpuVelocity(T *x, T *v, T *box, T *xcm, T *vcm, 
//         T *mass, T *second, T *omega, T *vext, T *v_lo, T *v_hi, T *coord_lo, T *coord_hi, 
//         T t_desired, T t_current, int *seed, int *save, int *map, int *image, int *type, 
//         int *coord_dim, int *vdim, int sum_flag, int dist_flag, int loop_flag, 
//         int rotation_flag, int momentum_flag, int triclinic, int dim, int mpiRank, 
//         int vmode, int nlocal, int natoms)
// {
//     if (vmode==0) {
//         cpuVelocityCreate(x, v, mass, second, omega,  box, xcm, vcm, t_desired, 
//                 t_current, seed, save, map, image, type, sum_flag, dist_flag, loop_flag, 
//                 rotation_flag, momentum_flag, triclinic, dim, mpiRank, nlocal, natoms);
//     } else if (vmode==1) {
//         cpuVelocitySet(v, vext,  vdim, sum_flag, dim, nlocal);        
//     } else if (vmode==2) {
//         cpuVelocityRamp(x, v, v_lo, v_hi, coord_lo, coord_hi, 
//                 coord_dim, vdim, sum_flag, dim, nlocal);        
//     } else if (vmode==3) {
//         //cpuVelocityScale(v, t_current, t_desired, mask, dim, nlocal);
//         cpuArrayMultiplyScalar(v, sqrt(t_desired/t_current), nlocal*dim);
//     } else if (vmode==4) {
//         cpuVelocityZeroMomentum(v, vcm, dim, nlocal);
//     } else if (vmode==5) {
//         cpuVelocityZeroRotation(x, v, box, xcm, omega, image, triclinic, dim, nlocal);
//     }    
// }
// template void cpuVelocity(double *, double *, double *, double *, 
//         double *, double *, double *, double *, double *, double *, double *, double *, double *, 
//         double, double, int *, int *, int *, int *, int *, int *, int *, int, int, int, 
//         int, int, int, int, int, int, int, int);
// template void cpuVelocity(float *, float *, float *, float *, 
//         float *, float *, float *, float *, float *, float *, float *, float *, float *, 
//         float, float, int *, int *, int *, int *, int *, int *, int *, int, int, int, 
//         int, int, int, int, int, int, int, int);

template <typename T> void cpuNVEInitialIntegrate(T *x, T *v, T *f, T *mass, T dtf, T dtv,
        int *type, int *ilist, int dim, int inum)
{
    for (int ii = 0; ii < inum; ii++) {
        int i = ilist[ii]; 
//         T dtf = dtarray[1];
//         T dtv = dtarray[2];
        T dtfm = dtf / mass[type[i]];
        v[i*dim+0] += dtfm * f[i*dim+0];
        v[i*dim+1] += dtfm * f[i*dim+1];
        x[i*dim+0] += dtv * v[i*dim+0];
        x[i*dim+1] += dtv * v[i*dim+1];        
        if (dim==3) {
            v[i*dim+2] += dtfm * f[i*dim+2];
            x[i*dim+2] += dtv * v[i*dim+2];
        }                
        //printf("%i %g %g %g",i,x[i*dim+0],x[i*dim+1],x[i*dim+2]);
    }      
}

template <typename T> void cpuNVEFinalIntegrate(T *x, T *v, T *f, T *mass, T dtf,
        int *type, int *ilist, int dim, int inum)
{
  // update v of atoms in group  
    for (int ii = 0; ii < inum; ii++) {
        int i = ilist[ii]; 
        //T dtf = dtarray[1];            
        T dtfm = dtf / mass[type[i]];
        v[i*dim+0] += dtfm * f[i*dim+0];
        v[i*dim+1] += dtfm * f[i*dim+1];
        if (dim==3) v[i*dim+2] += dtfm * f[i*dim+2];        
    }  
}

template <typename T> void cpuNVELimitInitialIntegrate(T *x, T *v, T *f, T *mass, T dtf, T dtv,
        T vlimitsq, int *type, int *ilist, int dim, int inum)
{
    if (dim==2) {
        for (int ii = 0; ii < inum; ii++) {
            int i = ilist[ii]; 
//             T dtf = dtarray[1];
//             T dtv = dtarray[2];                
            T dtfm = dtf / mass[type[i]];
            v[i*dim+0] += dtfm * f[i*dim+0];
            v[i*dim+1] += dtfm * f[i*dim+1];

            T vsq = v[i*dim+0]*v[i*dim+0] + v[i*dim+1]*v[i*dim+1];
            if (vsq > vlimitsq) {
              T scale = sqrt(vlimitsq/vsq);
              v[i*dim+0] *= scale;
              v[i*dim+1] *= scale;
            }

            x[i*dim+0] += dtv * v[i*dim+0];
            x[i*dim+1] += dtv * v[i*dim+1];
        }  
    }
    else {
        for (int ii = 0; ii < inum; ii++) {
            int i = ilist[ii]; 
//             T dtf = dtarray[1];
//             T dtv = dtarray[2];
            T dtfm = dtf / mass[type[i]];
            v[i*dim+0] += dtfm * f[i*dim+0];
            v[i*dim+1] += dtfm * f[i*dim+1];
            v[i*dim+2] += dtfm * f[i*dim+2];

            T vsq = v[i*dim+0]*v[i*dim+0] + v[i*dim+1]*v[i*dim+1] + v[i*dim+2]*v[i*dim+2];
            if (vsq > vlimitsq) {
              T scale = sqrt(vlimitsq/vsq);
              v[i*dim+0] *= scale;
              v[i*dim+1] *= scale;
              v[i*dim+2] *= scale;
            }

            x[i*dim+0] += dtv * v[i*dim+0];
            x[i*dim+1] += dtv * v[i*dim+1];
            x[i*dim+2] += dtv * v[i*dim+2];      
        }  
    }
}

template <typename T> void cpuNVELimitFinalIntegrate(T *x, T *v, T *f, T *mass, T dtf,
        T vlimitsq, int *type, int *ilist, int dim, int inum)
{ 
    if (dim==3) {
        for (int ii = 0; ii < inum; ii++) {
            int i = ilist[ii]; 
            //T dtf = dtarray[1];
            T dtfm = dtf / mass[type[i]];
            v[i*dim+0] += dtfm * f[i*dim+0];
            v[i*dim+1] += dtfm * f[i*dim+1];
            v[i*dim+2] += dtfm * f[i*dim+2];

            T vsq = v[i*dim+0]*v[i*dim+0] + v[i*dim+1]*v[i*dim+1] + v[i*dim+2]*v[i*dim+2];
            if (vsq > vlimitsq) {
              //ncount[0]++;
              T scale = sqrt(vlimitsq/vsq);
              v[i*dim+0] *= scale;
              v[i*dim+1] *= scale;
              v[i*dim+2] *= scale;
            }      
        }  
    }
    else {
        for (int ii = 0; ii < inum; ii++) {
            int i = ilist[ii]; 
            //T dtf = dtarray[1];            
            T dtfm = dtf / mass[type[i]];
            v[i*dim+0] += dtfm * f[i*dim+0];
            v[i*dim+1] += dtfm * f[i*dim+1];

            T vsq = v[i*dim+0]*v[i*dim+0] + v[i*dim+1]*v[i*dim+1];
            if (vsq > vlimitsq) {
              //ncount[0]++;
              T scale = sqrt(vlimitsq/vsq);
              v[i*dim+0] *= scale;
              v[i*dim+1] *= scale;
            }      
        }          
    }
}

template <typename T> void cpuSetVelocityInitialIntegrate(T *x, T *v, T *f, T *mass, T *fparam,
        T dtf, T dtv, int *type, int *ilist, int *iparam, int dim, int inum)
{    
    for (int ii = 0; ii < inum; ii++) {
        int i = ilist[ii]; 
//         T dtf = dtarray[1];
//         T dtv = dtarray[2];
        T dtfm = dtf / mass[type[i]];        
        for (int j=0; j<dim; j++) {
            v[i*dim+j] = (iparam[j]) ? fparam[j] : (v[i*dim+j] + dtfm * f[i*dim+j]);                
            x[i*dim+j] += dtv * v[i*dim+j];
        }
    }  
}
template void cpuSetVelocityInitialIntegrate(double *x, double *v, double *f, double *mass, double *fparam,
        double dtf, double dtv, int *type, int *ilist, int *iparam, int dim, int inum);
template void cpuSetVelocityInitialIntegrate(float *x, float *v, float *f, float *mass, float *fparam,
        float dtf, float dtv, int *type, int *ilist, int *iparam, int dim, int inum);

template <typename T> void cpuSetVelocityFinalIntegrate(T *x, T *v, T *f, T *mass, 
        T dtf, int *type, int *ilist, int *iparam, int dim, int inum)
{
  // update v of atoms in group  
    for (int ii = 0; ii < inum; ii++) {
        int i = ilist[ii]; 
        T dtfm = dtf / mass[type[i]];
        for (int j=0; j<dim; j++) 
            v[i*dim+j] = (iparam[j]) ? v[i*dim+j] : (v[i*dim+j] + dtfm * f[i*dim+j]);                         
    }  
}
template void cpuSetVelocityFinalIntegrate(double *x, double *v, double *f, double *mass, 
        double dtf, int *type, int *ilist, int *iparam, int dim, int inum);
template void cpuSetVelocityFinalIntegrate(float *x, float *v, float *f, float *mass, 
        float dtf, int *type, int *ilist, int *iparam, int dim, int inum);

template <typename T> void cpuVelocityFactor(T *v, T factor_eta, int *ilist, 
        int biasflag, int dim, int inum)
{
  if (dim==3) {
      if (biasflag == 0) {
        for (int ii = 0; ii < inum; ii++) {
            int i = ilist[ii]; 
            v[i*dim+0] *= factor_eta;
            v[i*dim+1] *= factor_eta;
            v[i*dim+2] *= factor_eta;      
        }
      } else if (biasflag == 1)  {
        for (int ii = 0; ii < inum; ii++) {
            int i = ilist[ii]; 
            //remove_bias(i,&v[i*dim]);
            v[i*dim+0] *= factor_eta;
            v[i*dim+1] *= factor_eta;
            v[i*dim+2] *= factor_eta;
            //restore_bias(i,&v[i*dim]);      
        }
      }
  } else {
      if (biasflag == 0) {
        for (int ii = 0; ii < inum; ii++) {
            int i = ilist[ii]; 
            v[i*dim+0] *= factor_eta;
            v[i*dim+1] *= factor_eta;
        }
      } else if (biasflag == 1)  {
        for (int ii = 0; ii < inum; ii++) {
            int i = ilist[ii]; 
            //remove_bias(i,&v[i*dim]);
            v[i*dim+0] *= factor_eta;
            v[i*dim+1] *= factor_eta;
            //restore_bias(i,&v[i*dim]);      
        }
      }      
  }
}

template <typename T> void cpuNoseHooverThermostat(T *v, T *dtarray, T *tarray, T *eta_mass, T *eta, 
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

  //printf("%i %i %i %i %g %g %g %g %g \n", biasflag, eta_mass_flag, mtchain, nc_tchain, tdof, boltz, t_freq, t_target, t_current);
  
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
    cpuVelocityFactor(v, factor_eta, ilist, biasflag, dim, inum);

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

template <typename T> void cpuNVTInitialIntegrate(T *x, T *v, T *f, T *mass, T *dtarray, T *tarray,
        T *eta_mass, T *eta, T *eta_dot, T *eta_dotdot, int *type, int *ilist, 
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
  T t_current = cpuComputeTempScalar(v, mass, tfactor, type, ilist, dim, inum);
  tarray[4] = t_current; 
  
  cpuNoseHooverThermostat(v, dtarray, tarray, eta_mass, eta, eta_dot, eta_dotdot, 
           ilist, eta_mass_flag, biasflag, mtchain, nc_tchain, dim, inum);

  cpuNVEInitialIntegrate(x, v, f, mass, dtf, dtv, type, ilist, dim, inum);  
}


template <typename T> void cpuNVTFinalIntegrate(T *x, T *v, T *f, T *mass, T *dtarray, T *tarray,
        T *eta_mass, T *eta, T *eta_dot, T *eta_dotdot, int *type, int *ilist, 
        int eta_mass_flag, int biasflag, int mtchain, int nc_tchain, int dim, int inum)
{
  T dtf = dtarray[1];
  T dtv = dtarray[2];    
  cpuNVEFinalIntegrate(x, v, f, mass, dtf, type, ilist, dim, inum);  

  T tdof = tarray[5];
  T boltz = tarray[6];
  T mvv2e = tarray[8];    
  T tfactor = mvv2e / (tdof * boltz);    
  T t_current = cpuComputeTempScalar(v, mass, tfactor, type, ilist, dim, inum);
  tarray[4] = t_current; 
      
  cpuNoseHooverThermostat(v, dtarray, tarray, eta_mass, eta, eta_dot, eta_dotdot, 
          ilist, eta_mass_flag, biasflag, mtchain, nc_tchain, dim, inum);
}

template <typename T> void cpuInitialIntegrate(T *x, T *v, T *f, T *mass, T *dtarray, T *tarray,
        T *eta_mass, T *eta, T *eta_dot, T *eta_dotdot, T vlimitsq, int *type, int *ilist,
        int eta_mass_flag, int biasflag, int mtchain, int nc_tchain, int mode, int dim, int inum)
{
    T dtf = dtarray[1];
    T dtv = dtarray[2];        
    if (mode==0) {
        cpuNVEInitialIntegrate(x, v, f, mass, dtf, dtv,
            type, ilist, dim, inum);
    }
    else if (mode==1) {
        cpuNVELimitInitialIntegrate(x, v, f, mass, dtf, dtv,
            vlimitsq, type, ilist, dim, inum);
    }
    else if (mode==2) {
        cpuNVTInitialIntegrate(x, v, f, mass, dtarray, tarray,
            eta_mass, eta, eta_dot, eta_dotdot, type, ilist, 
            eta_mass_flag, biasflag, mtchain, nc_tchain, dim, inum);
    }
//     else if (mode==3) {
//         cpuMoveInitialIntegrate(x, v, f, mass, fparam, dtf, dtv, type, ilist, iparam, dim, inum);
//     }
}
template void cpuInitialIntegrate(double *x, double *v, double *f, double *mass, 
        double *dtarray, double *tarray, double *eta_mass, double *eta, double *eta_dot, 
        double *eta_dotdot, double vlimitsq, int *type, int *ilist, int eta_mass_flag, 
        int biasflag,  int mtchain, int nc_tchain, int mode, int dim, int inum);
template void cpuInitialIntegrate(float *x, float *v, float *f, float *mass, 
        float *dtarray, float *tarray, float *eta_mass, float *eta, float *eta_dot, 
        float *eta_dotdot, float vlimitsq, int *type, int *ilist, int eta_mass_flag, 
        int biasflag, int mtchain, int nc_tchain, int mode, int dim, int inum);


template <typename T> void cpuFinalIntegrate(T *x, T *v, T *f, T *mass, T *dtarray, T *tarray,
        T *eta_mass, T *eta, T *eta_dot, T *eta_dotdot, T vlimitsq, int *type, int *ilist, 
        int eta_mass_flag, int biasflag, int mtchain, int nc_tchain, int mode, int dim, int inum)
{    
    T dtf = dtarray[1];
    if (mode==0) {
        cpuNVEFinalIntegrate(x, v, f, mass, dtf, type, ilist, dim, inum);
    }
    else if (mode==1) {
        cpuNVELimitFinalIntegrate(x, v, f, mass, dtf, vlimitsq, type, ilist, dim, inum);        
    }
    else if (mode==2) {
        cpuNVTFinalIntegrate(x, v, f, mass, dtarray, tarray, eta_mass, eta, eta_dot, eta_dotdot, 
                type, ilist, eta_mass_flag, biasflag, mtchain, nc_tchain, dim, inum);        
    }
}
template void cpuFinalIntegrate(double *x, double *v, double *f, double *mass, 
        double *dtarray, double *tarray, double *eta_mass, double *eta, double *eta_dot, 
        double *eta_dotdot, double vlimitsq, int *type, int *ilist, int eta_mass_flag, 
        int biasflag, int mtchain, int nc_tchain, int mode, int dim, int inum);
template void cpuFinalIntegrate(float *x, float *v, float *f, float *mass, 
        float *dtarray, float *tarray, float *eta_mass, float *eta, float *eta_dot, 
        float *eta_dotdot, float vlimitsq, int *type, int *ilist, int eta_mass_flag, 
        int biasflag, int mtchain, int nc_tchain, int mode, int dim, int inum);

template <typename T> T cpuBerendsenThermostat(T *v, T *mass, T *dtarray, T *tarray, T energy, 
        int *type, int *ilist, int biasflag, int dim, int inum)
{
  T tdof = tarray[5];
  T boltz = tarray[6];
  T mvv2e = tarray[8];      
  T tfactor = mvv2e / (tdof * boltz);        
  T t_current = cpuComputeTempScalar(v, mass, tfactor, type, ilist, dim, inum);
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

  cpuVelocityFactor(v, lamda, ilist, biasflag, dim, inum);
    
  return energy;
}

template <typename T> T cpuRescaleThermostat(T *v, T *mass, T *dtarray, T *tarray, T energy, 
        int *type, int *ilist, int biasflag, int dim, int inum)
{
  T tdof = tarray[5];
  T boltz = tarray[6];
  T mvv2e = tarray[8];      
  T t_window = tarray[9];    
  T fraction = tarray[10];    
  T tfactor = mvv2e / (tdof * boltz);        
  T t_current = cpuComputeTempScalar(v, mass, tfactor, type, ilist, dim, inum);
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
    cpuVelocityFactor(v, factor, ilist, biasflag, dim, inum);    
  }  
  
  return energy;
}

template <typename T> T cpuCsvrThermostat(T *v, T *mass, T *dtarray, T *tarray, T *second, 
        T energy, int *type, int *ilist, int *seed, int *save, int biasflag, int dim, int inum)
{
  T tdof = tarray[5];
  T boltz = tarray[6];
  T mvv2e = tarray[8];      
  T tfactor = mvv2e / (tdof * boltz);        
  T t_current = cpuComputeTempScalar(v, mass, tfactor, type, ilist, dim, inum);
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
  T r1 = cpuRandomGaussian(seed, save, second);
  //T r2 = sumnoises(tdof - 1);
  T r2;
  int nn = (int) (tdof-1);
  if (nn == 0) {
    r2 = 0.0;
  } else if (nn == 1) {
    T rr = cpuRandomGaussian(seed, save, second);
    r2 = rr*rr;
  } else if (nn % 2 == 0) {
    r2 = 2.0 * cpuGamdev(nn / 2, seed);
  } else {
    T rr = cpuRandomGaussian(seed, save, second);
    r2 =  2.0 * cpuGamdev((nn-1) / 2, seed) + rr*rr;
  }
  
  T scale = c1 + c2*(r1*r1+r2) + 2.0*r1*sqrt(c1*c2);
  T lamda = sqrt(scale);
  
  energy += t_current * (1.0-lamda*lamda) * efactor;

  cpuVelocityFactor(v, lamda, ilist, biasflag, dim, inum);
  
  return energy;
}

template <typename T> T cpuVelocityRescalingThermostat(T *v, T *mass, T *dtarray, T *tarray, T *second, 
        T energy, int *type, int *ilist, int *seed, int *save, int biasflag, int mode, int dim, int inum)
{    
    if (mode==0) {
        return cpuBerendsenThermostat(v, mass, dtarray, tarray, energy, type, ilist, biasflag, dim, inum);
    }
    else if (mode==1) {
        return cpuRescaleThermostat(v, mass, dtarray, tarray,  energy, type, ilist, biasflag, dim, inum);
    }
    else if (mode==2) {
        return cpuCsvrThermostat(v, mass, dtarray, tarray, second, energy, type, ilist, 
                seed, save, biasflag, dim, inum);
    }
    return 0.0;
}
template double cpuVelocityRescalingThermostat(double *v, double *mass, double *dtarray, double *tarray, double *second, 
        double energy, int *type, int *ilist, int *seed, int *save, int biasflag, int mode, int dim, int inum);

template float cpuVelocityRescalingThermostat(float *v, float *mass, float *dtarray, float *tarray, float *second, 
        float energy, int *type, int *ilist, int *seed, int *save, int biasflag, int mode, int dim, int inum);


#endif

