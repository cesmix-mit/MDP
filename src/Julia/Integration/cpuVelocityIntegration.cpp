/***************************************************************************                               
                    Molecular Dynamics Potentials (MDP)
                           CESMIX-MIT Project  
 
 Contributing authors: Ngoc-Cuong Nguyen (cuongng@mit.edu, exapde@gmail.com)
 ***************************************************************************/
               
#define WARMUP 100

void cpuVelocityZeroMomentum(double *v, double *vcm, int dim, int nlocal)
{  
  for (int i = 0; i < nlocal; i++) {
      v[i*dim+0] -= vcm[0];
      v[i*dim+1] -= vcm[1];
      if (dim==3) v[i*dim+2] -= vcm[2];
  }
}

void cpuVelocityZeroRotation(double *x, double *v, double *box, double *xcm, double *omega, 
        int *image, int triclinic, int dim, int nlocal)
{  
  double unwrap[3];  
  if (dim==3) {
      for (int i = 0; i < nlocal; i++) {        
          cpuUnmap(unwrap, &x[i*dim], box, &image[i*dim], triclinic, dim);        
          double dx = unwrap[0] - xcm[0];
          double dy = unwrap[1] - xcm[1];
          double dz = unwrap[2] - xcm[2];
          v[i*dim+0] -= omega[1]*dz - omega[2]*dy;
          v[i*dim+1] -= omega[2]*dx - omega[0]*dz;
          v[i*dim+2] -= omega[0]*dy - omega[1]*dx;
     }
  }
  else {
      for (int i = 0; i < nlocal; i++) {  
          cpuUnmap(unwrap, &x[i*dim], box, &image[i*dim], triclinic, dim);        
          double dx = unwrap[0] - xcm[0];
          double dy = unwrap[1] - xcm[1];
          v[i*dim+0] -= -omega[2]*dy;
          v[i*dim+1] -=  omega[2]*dx;
      }      
  }
}

void cpuVelocityCreate(double *x, double *v, double *mass, double *second, 
         int *seed, int *save, int *map, int *type, int seed0, int sum_flag, int dist_flag, int loop_flag, 
        int dim, int mpiRank, int nlocal, int natoms)
{
  int i;
  
  if (sum_flag==0) {
      for (i = 0; i < dim*nlocal; i++)
          v[i] = 0.0;
  }

  int m;
  double vx,vy,vz,factor;

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

void cpuVelocityCreate(double *x, double *v, double *mass, double *second, double *omega, 
        double *box, double *xcm, double *vcm, double t_desired, double t_current, int *seed, int *save, int *map, int *image, 
        int *type, int sum_flag, int dist_flag, int loop_flag, int rotation_flag, 
        int momentum_flag, int triclinic, int dim, int mpiRank, int nlocal, int natoms)
{
  int i;
  
  if (sum_flag==0) {
      for (i = 0; i < dim*nlocal; i++)
          v[i] = 0.0;
  }

  int m;
  double vx,vy,vz,factor;

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

void cpuVelocitySet(double *v, double *vext, int *vdim,
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

void cpuVelocityRamp(double *x, double *v, double *v_lo, double *v_hi, double *coord_lo, double *coord_hi,
        int *coord_dim, int *v_dim, int sum_flag, int dim, int nlocal)
{
  for (int i = 0; i < nlocal; i++)
      for (int j = 0; j < dim; j++) {  
          double fraction = (x[i*dim+coord_dim[j]] - coord_lo[j]) / (coord_hi[j] - coord_lo[j]);
          fraction = (fraction > 0.0) ? fraction : 0.0;
          fraction = (fraction < 1.0) ? fraction : 1.0;
          double vramp = v_lo[j] + fraction*(v_hi[j] - v_lo[j]);
          if (sum_flag) v[i*dim+v_dim[j]] += vramp;
          else v[i*dim+v_dim[j]] = vramp;
      }    
}

void cpuNVEInitialIntegrate(double *x, double *v, double *f, double *mass, double dtf, double dtv,
        int *type, int *ilist, int dim, int inum)
{
    for (int ii = 0; ii < inum; ii++) {
        int i = ilist[ii]; 
//         double dtf = dtarray[1];
//         double dtv = dtarray[2];
        double dtfm = dtf / mass[type[i]];
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

void cpuNVEFinalIntegrate(double *x, double *v, double *f, double *mass, double dtf,
        int *type, int *ilist, int dim, int inum)
{
  // update v of atoms in group  
    for (int ii = 0; ii < inum; ii++) {
        int i = ilist[ii]; 
        //double dtf = dtarray[1];            
        double dtfm = dtf / mass[type[i]];
        v[i*dim+0] += dtfm * f[i*dim+0];
        v[i*dim+1] += dtfm * f[i*dim+1];
        if (dim==3) v[i*dim+2] += dtfm * f[i*dim+2];        
    }  
}

void cpuNVELimitInitialIntegrate(double *x, double *v, double *f, double *mass, double dtf, double dtv,
        double vlimitsq, int *type, int *ilist, int dim, int inum)
{
    if (dim==2) {
        for (int ii = 0; ii < inum; ii++) {
            int i = ilist[ii]; 
//             double dtf = dtarray[1];
//             double dtv = dtarray[2];                
            double dtfm = dtf / mass[type[i]];
            v[i*dim+0] += dtfm * f[i*dim+0];
            v[i*dim+1] += dtfm * f[i*dim+1];

            double vsq = v[i*dim+0]*v[i*dim+0] + v[i*dim+1]*v[i*dim+1];
            if (vsq > vlimitsq) {
              double scale = sqrt(vlimitsq/vsq);
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
//             double dtf = dtarray[1];
//             double dtv = dtarray[2];
            double dtfm = dtf / mass[type[i]];
            v[i*dim+0] += dtfm * f[i*dim+0];
            v[i*dim+1] += dtfm * f[i*dim+1];
            v[i*dim+2] += dtfm * f[i*dim+2];

            double vsq = v[i*dim+0]*v[i*dim+0] + v[i*dim+1]*v[i*dim+1] + v[i*dim+2]*v[i*dim+2];
            if (vsq > vlimitsq) {
              double scale = sqrt(vlimitsq/vsq);
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

void cpuNVELimitFinalIntegrate(double *x, double *v, double *f, double *mass, double dtf,
        double vlimitsq, int *type, int *ilist, int dim, int inum)
{ 
    if (dim==3) {
        for (int ii = 0; ii < inum; ii++) {
            int i = ilist[ii]; 
            //double dtf = dtarray[1];
            double dtfm = dtf / mass[type[i]];
            v[i*dim+0] += dtfm * f[i*dim+0];
            v[i*dim+1] += dtfm * f[i*dim+1];
            v[i*dim+2] += dtfm * f[i*dim+2];

            double vsq = v[i*dim+0]*v[i*dim+0] + v[i*dim+1]*v[i*dim+1] + v[i*dim+2]*v[i*dim+2];
            if (vsq > vlimitsq) {
              //ncount[0]++;
              double scale = sqrt(vlimitsq/vsq);
              v[i*dim+0] *= scale;
              v[i*dim+1] *= scale;
              v[i*dim+2] *= scale;
            }      
        }  
    }
    else {
        for (int ii = 0; ii < inum; ii++) {
            int i = ilist[ii]; 
            //double dtf = dtarray[1];            
            double dtfm = dtf / mass[type[i]];
            v[i*dim+0] += dtfm * f[i*dim+0];
            v[i*dim+1] += dtfm * f[i*dim+1];

            double vsq = v[i*dim+0]*v[i*dim+0] + v[i*dim+1]*v[i*dim+1];
            if (vsq > vlimitsq) {
              //ncount[0]++;
              double scale = sqrt(vlimitsq/vsq);
              v[i*dim+0] *= scale;
              v[i*dim+1] *= scale;
            }      
        }          
    }
}

void cpuSetVelocityInitialIntegrate(double *x, double *v, double *f, double *mass, double *fparam,
        double dtf, double dtv, int *type, int *ilist, int *iparam, int dim, int inum)
{    
    for (int ii = 0; ii < inum; ii++) {
        int i = ilist[ii]; 
//         double dtf = dtarray[1];
//         double dtv = dtarray[2];
        double dtfm = dtf / mass[type[i]];        
        for (int j=0; j<dim; j++) {
            v[i*dim+j] = (iparam[j]) ? fparam[j] : (v[i*dim+j] + dtfm * f[i*dim+j]);                
            x[i*dim+j] += dtv * v[i*dim+j];
        }
    }  
}

void cpuSetVelocityFinalIntegrate(double *x, double *v, double *f, double *mass, 
        double dtf, int *type, int *ilist, int *iparam, int dim, int inum)
{
  // update v of atoms in group  
    for (int ii = 0; ii < inum; ii++) {
        int i = ilist[ii]; 
        double dtfm = dtf / mass[type[i]];
        for (int j=0; j<dim; j++) 
            v[i*dim+j] = (iparam[j]) ? v[i*dim+j] : (v[i*dim+j] + dtfm * f[i*dim+j]);                         
    }  
}

void cpuVelocityFactor(double *v, double factor_eta, int *ilist, 
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

void cpuNoseHooverThermostat(double *v, double *dtarray, double *tarray, double *eta_mass, double *eta, 
        double *eta_dot, double *eta_dotdot, int *ilist, int eta_mass_flag, int biasflag, int mtchain, 
        int nc_tchain, int dim, int inum)
{
  int ich;
  double dt = dtarray[0];
  //double dtf = dtarray[1];
  //double dtv = dtarray[2];
  double dthalf = 0.5 * dt;
  double dt4 = 0.25 * dt;
  double dt8 = 0.125 * dt;
  
  //double t_start = tarray[0];  
  //double t_stop = tarray[1];     
  //double t_freq = tarray[2];
  //double t_target = tarray[3]; 
  //double t_current = tarray[4];   
  //double tdof = tarray[5];
  //double boltz = tarray[6];
  //double drag = tarray[7];
  //double mvv2e = tarray[8];    
  //double tfactor = mvv2e / (tdof * boltz);
      
  double t_freq = tarray[2];
  double t_target = tarray[3]; 
  //double t_current = tarray[4];   
  double tdof = tarray[5];
  double boltz = tarray[6];
  double drag = tarray[7];
  double tdrag_factor = 1.0 - (dt * t_freq * drag / nc_tchain);
  
  double expfac;
  double kecurrent = tdof * boltz * tarray[4];
  double ke_target = tdof * boltz * t_target;
  
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
  
  double ncfac = 1.0/nc_tchain;
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

    double factor_eta = exp(-ncfac*dthalf*eta_dot[0]);
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

void cpuNVTInitialIntegrate(double *x, double *v, double *f, double *mass, double *dtarray, double *tarray,
        double *eta_mass, double *eta, double *eta_dot, double *eta_dotdot, int *type, int *ilist, 
        int eta_mass_flag, int biasflag, int mtchain, int nc_tchain, int dim, int inum)
{  
  double dtf = dtarray[1];
  double dtv = dtarray[2];
  double beginstep = dtarray[3];
  double endstep = dtarray[4];
  double ntimestep = dtarray[5];      
  
  double t_start = tarray[0];  
  double t_stop = tarray[1];   
  double delta = (ntimestep - beginstep)/(endstep - beginstep);  
  tarray[3] = t_start + delta * (t_stop - t_start); // t_target
    
  double tdof = tarray[5];
  double boltz = tarray[6];
  double mvv2e = tarray[8];    
  double tfactor = mvv2e / (tdof * boltz);    
  double t_current = cpuComputeTempScalar(v, mass, tfactor, type, ilist, dim, inum);
  tarray[4] = t_current; 
  
  cpuNoseHooverThermostat(v, dtarray, tarray, eta_mass, eta, eta_dot, eta_dotdot, 
           ilist, eta_mass_flag, biasflag, mtchain, nc_tchain, dim, inum);

  cpuNVEInitialIntegrate(x, v, f, mass, dtf, dtv, type, ilist, dim, inum);  
}


void cpuNVTFinalIntegrate(double *x, double *v, double *f, double *mass, double *dtarray, double *tarray,
        double *eta_mass, double *eta, double *eta_dot, double *eta_dotdot, int *type, int *ilist, 
        int eta_mass_flag, int biasflag, int mtchain, int nc_tchain, int dim, int inum)
{
  double dtf = dtarray[1];
  //double dtv = dtarray[2];    
  cpuNVEFinalIntegrate(x, v, f, mass, dtf, type, ilist, dim, inum);  

  double tdof = tarray[5];
  double boltz = tarray[6];
  double mvv2e = tarray[8];    
  double tfactor = mvv2e / (tdof * boltz);    
  double t_current = cpuComputeTempScalar(v, mass, tfactor, type, ilist, dim, inum);
  tarray[4] = t_current; 
      
  cpuNoseHooverThermostat(v, dtarray, tarray, eta_mass, eta, eta_dot, eta_dotdot, 
          ilist, eta_mass_flag, biasflag, mtchain, nc_tchain, dim, inum);
}

void cpuInitialIntegrate(double *x, double *v, double *f, double *mass, double *dtarray, double *tarray,
        double *eta_mass, double *eta, double *eta_dot, double *eta_dotdot, double vlimitsq, int *type, int *ilist,
        int eta_mass_flag, int biasflag, int mtchain, int nc_tchain, int mode, int dim, int inum)
{
    double dtf = dtarray[1];
    double dtv = dtarray[2];        
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

void cpuFinalIntegrate(double *x, double *v, double *f, double *mass, double *dtarray, double *tarray,
        double *eta_mass, double *eta, double *eta_dot, double *eta_dotdot, double vlimitsq, int *type, int *ilist, 
        int eta_mass_flag, int biasflag, int mtchain, int nc_tchain, int mode, int dim, int inum)
{    
    double dtf = dtarray[1];
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

double cpuBerendsenThermostat(double *v, double *mass, double *dtarray, double *tarray, double energy, 
        int *type, int *ilist, int biasflag, int dim, int inum)
{
  double tdof = tarray[5];
  double boltz = tarray[6];
  double mvv2e = tarray[8];      
  double tfactor = mvv2e / (tdof * boltz);        
  double t_current = cpuComputeTempScalar(v, mass, tfactor, type, ilist, dim, inum);
  tarray[4] = t_current;         
  
  if (tdof < 1) return energy;
  //if (t_current == 0.0)
  //  error("Computed temperature for fix temp/berendsen cannot be 0.0");

  double dt = dtarray[0];
  double beginstep = dtarray[3];
  double endstep = dtarray[4];
  double ntimestep = dtarray[5];        
  double t_start = tarray[0];  
  double t_stop = tarray[1];   
  double t_freq = tarray[2];
  double delta = (ntimestep - beginstep)/(endstep - beginstep);  
  double t_target  = t_start + delta * (t_stop - t_start); // 
  tarray[3] = t_target;
    
  double lamda = sqrt(1.0 + dt*t_freq*(t_target/t_current - 1.0));
  double efactor = 0.5 * boltz * tdof;
  energy += t_current * (1.0-lamda*lamda) * efactor;

  cpuVelocityFactor(v, lamda, ilist, biasflag, dim, inum);
    
  return energy;
}

double cpuRescaleThermostat(double *v, double *mass, double *dtarray, double *tarray, double energy, 
        int *type, int *ilist, int biasflag, int dim, int inum)
{
  double tdof = tarray[5];
  double boltz = tarray[6];
  double mvv2e = tarray[8];      
  double t_window = tarray[9];    
  double fraction = tarray[10];    
  double tfactor = mvv2e / (tdof * boltz);        
  double t_current = cpuComputeTempScalar(v, mass, tfactor, type, ilist, dim, inum);
  tarray[4] = t_current;         
  
  if (tdof < 1) return energy;
  //if (t_current == 0.0)
  //  error("Computed temperature for fix temp/berendsen cannot be 0.0");

  double dt = dtarray[0];
  double beginstep = dtarray[3];
  double endstep = dtarray[4];
  double ntimestep = dtarray[5];        
  double t_start = tarray[0];  
  double t_stop = tarray[1];   
  double t_freq = tarray[2];
  double delta = (ntimestep - beginstep)/(endstep - beginstep);  
  double t_target  = t_start + delta * (t_stop - t_start); // 
  tarray[3] = t_target;
    
  double lamda = sqrt(1.0 + dt*t_freq*(t_target/t_current - 1.0));
  double efactor = 0.5 * boltz * tdof;
  energy += t_current * (1.0-lamda*lamda) * efactor;

  if (fabs(t_current-t_target) > t_window) {
    t_target = t_current - fraction*(t_current-t_target);
    double factor = sqrt(t_target/t_current);
    double efactor = 0.5 * boltz * tdof;
    energy += (t_current-t_target) * efactor;      
    cpuVelocityFactor(v, factor, ilist, biasflag, dim, inum);    
  }  
  
  return energy;
}

double cpuCsvrThermostat(double *v, double *mass, double *dtarray, double *tarray, double *second, 
        double energy, int *type, int *ilist, int *seed, int *save, int biasflag, int dim, int inum)
{
  double tdof = tarray[5];
  double boltz = tarray[6];
  double mvv2e = tarray[8];      
  double tfactor = mvv2e / (tdof * boltz);        
  double t_current = cpuComputeTempScalar(v, mass, tfactor, type, ilist, dim, inum);
  tarray[4] = t_current;         
  
  if (tdof < 1) return energy;
  //if (t_current == 0.0)
  //  error("Computed temperature for fix temp/berendsen cannot be 0.0");

  double dt = dtarray[0];
  double beginstep = dtarray[3];
  double endstep = dtarray[4];
  double ntimestep = dtarray[5];        
  double t_start = tarray[0];  
  double t_stop = tarray[1];   
  double t_freq = tarray[2];
  double delta = (ntimestep - beginstep)/(endstep - beginstep);  
  double t_target  = t_start + delta * (t_stop - t_start); // 
  tarray[3] = t_target;
    
  //double lamda = sqrt(1.0 + dt*t_freq*(t_target/t_current - 1.0));
  double efactor = 0.5 * boltz * tdof;
  double ekin_old = t_current * efactor;
  double ekin_new = t_target * efactor;
  
  double c1 = exp(-dt*t_freq);
  double c2 = (1.0-c1)*ekin_new/ekin_old/tdof;
  double r1 = cpuRandomGaussian(seed, save, second);
  //double r2 = sumnoises(tdof - 1);
  double r2;
  int nn = (int) (tdof-1);
  if (nn == 0) {
    r2 = 0.0;
  } else if (nn == 1) {
    double rr = cpuRandomGaussian(seed, save, second);
    r2 = rr*rr;
  } else if (nn % 2 == 0) {
    r2 = 2.0 * cpuGamdev(nn / 2, seed);
  } else {
    double rr = cpuRandomGaussian(seed, save, second);
    r2 =  2.0 * cpuGamdev((nn-1) / 2, seed) + rr*rr;
  }
  
  double scale = c1 + c2*(r1*r1+r2) + 2.0*r1*sqrt(c1*c2);
  double lamda = sqrt(scale);
  
  energy += t_current * (1.0-lamda*lamda) * efactor;

  cpuVelocityFactor(v, lamda, ilist, biasflag, dim, inum);
  
  return energy;
}

double cpuVelocityRescalingThermostat(double *v, double *mass, double *dtarray, double *tarray, double *second, 
        double energy, int *type, int *ilist, int *seed, int *save, int biasflag, int mode, int dim, int inum)
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


