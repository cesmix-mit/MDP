/***************************************************************************                               
                    Molecular Dynamics Potentials (MDP)
                           CESMIX-MIT Project  
 
 Contributing authors: Ngoc-Cuong Nguyen (cuongng@mit.edu, exapde@gmail.com)
 ***************************************************************************/
          
// cpuFixSetForce, cpuFixLineForce, cpuFixPlaneForce, cpuFixAddForce, cpuFixAveForce, cpuFixDrag, 
// cpuFixWallReflect, cpuFixWallHarmonic, cpuFixWallLJ93, cpuFixWallLJ126, cpuFixWallLJ1043, cpuFixWallMorse 
        
void cpuFixSetForce(double *x, double *v, double *f, double *eatom, double *vatom, double *fparam, 
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
         
void cpuFixLineForce(double *x, double *v, double *f, double *eatom, double *vatom, double *fparam, 
        int *iparam, int *ilist, int eflag_atom, int vflag_atom, int dim, int inum)        
{
    if (dim==3) {
        for (int ii = 0; ii < inum; ii++) {
            int i = ilist[ii]; 
            double dot = f[i*dim+0]*fparam[0] + f[i*dim+1]*fparam[1] + f[i*dim+2]*fparam[2];
            f[i*dim+0] = dot * fparam[0];
            f[i*dim+1] = dot * fparam[1];
            f[i*dim+2] = dot * fparam[1];            
        }
    } else {
        for (int ii = 0; ii < inum; ii++) {
            int i = ilist[ii]; 
            double dot = f[i*dim+0]*fparam[0] + f[i*dim+1]*fparam[1];
            f[i*dim+0] = dot * fparam[0];
            f[i*dim+1] = dot * fparam[1];
        }        
    }
}

void cpuFixPlaneForce(double *x, double *v, double *f, double *eatom, double *vatom, double *fparam, 
        int *iparam, int *ilist, int eflag_atom, int vflag_atom, int dim, int inum)        
{
    if (dim==3) {
        for (int ii = 0; ii < inum; ii++) {
            int i = ilist[ii]; 
            double dot = f[i*dim+0]*fparam[0] + f[i*dim+1]*fparam[1] + f[i*dim+2]*fparam[2];
            f[i*dim+0] -= dot * fparam[0];
            f[i*dim+1] -= dot * fparam[1];
            f[i*dim+2] -= dot * fparam[1];            
        }
    } else {
        for (int ii = 0; ii < inum; ii++) {
            int i = ilist[ii]; 
            double dot = f[i*dim+0]*fparam[0] + f[i*dim+1]*fparam[1];
            f[i*dim+0] -= dot * fparam[0];
            f[i*dim+1] -= dot * fparam[1];
        }        
    }
}

void cpuFixAddForce(double *x, double *v, double *f, double *eatom, double *vatom, double *fparam, 
        int *iparam, int *ilist, int eflag_atom, int vflag_atom, int dim, int inum)
{    
    if (dim==3) {
        for (int ii = 0; ii < inum; ii++) {
            int i = ilist[ii]; 
            double y0 = x[i*dim+0];
            double y1 = x[i*dim+1];
            double y2 = x[i*dim+2];  
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
            double y0 = x[i*dim+0];
            double y1 = x[i*dim+1];
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

void cpuFixAddForce(double *x, double *v, double *f, double *eatom, double *vatom, double *fparam, double *box,
        int *iparam, int *ilist, int *image, int triclinic, int eflag_atom, int vflag_atom, int dim, int inum)
{    
    if (dim==3) {
        for (int ii = 0; ii < inum; ii++) {
            int i = ilist[ii]; 
            double y0, y1, y2;  
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
            double y0, y1;  
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

void cpuFixAveForce(double *x, double *v, double *f, double *eatom, double *vatom, double *fparam, 
        int *iparam, int *ilist, int eflag_atom, int vflag_atom, int dim, int inum)
{
  double foriginal[4];
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

void cpuFixDragForce(double *x, double *v, double *f, double *eatom, double *vatom, double *fparam, double *box,
        int *iparam, int *ilist, int *pbc, int triclinic, int eflag_atom, int vflag_atom, int dim, int inum)
{
  if (dim==3) {
      for (int ii = 0; ii < inum; ii++) {
          int i = ilist[ii]; 
          double f_mag = fparam[0];
          double delta = fparam[1];
          double dp[3]; 
          dp[0] = x[i*dim+0] - fparam[2];
          dp[1] = x[i*dim+1] - fparam[3];
          dp[2] = x[i*dim+2] - fparam[4];
          if (!iparam[0]) dp[0] = 0.0;
          if (!iparam[1]) dp[1] = 0.0;
          if (!iparam[2]) dp[2] = 0.0;      
          cpuMinimumImage(dp, box, pbc, triclinic, dim);
          double r = sqrt(dp[0]*dp[0] + dp[1]*dp[1] + dp[2]*dp[2]);
          if (r > delta) {
            double prefactor = f_mag/r;
            double fx = prefactor*dp[0];
            double fy = prefactor*dp[1];
            double fz = prefactor*dp[2];
            f[i*dim+0] -= fx;
            f[i*dim+1] -= fy;
            f[i*dim+2] -= fz;
          }
      }
  } else {
      for (int ii = 0; ii < inum; ii++) {
          int i = ilist[ii]; 
          double f_mag = fparam[0];
          double delta = fparam[1];          
          double dp[3]; 
          dp[0] = x[i*dim+0] - fparam[2];
          dp[1] = x[i*dim+1] - fparam[3];
          if (!iparam[0]) dp[0] = 0.0;
          if (!iparam[1]) dp[1] = 0.0;
          cpuMinimumImage(dp, box, pbc, triclinic, dim);
          double r = sqrt(dp[0]*dp[0] + dp[1]*dp[1]);
          if (r > delta) {
            double prefactor = f_mag/r;
            double fx = prefactor*dp[0];
            double fy = prefactor*dp[1];
            f[i*dim+0] -= fx;
            f[i*dim+1] -= fy;
          }
      }      
  }  
}
        
void cpuFixWallReflect(double *x, double *v, double *f, double *eatom, double *vatom, double *fparam, 
        int *iparam, int *ilist, int eflag_atom, int vflag_atom, int dim, int inum)
{  
  for (int ii = 0; ii < inum; ii++) {
      int i = ilist[ii]; 
      int direction = iparam[0] / 2;
      int side = iparam[0] % 2;      
      double coord = fparam[0];
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
 
void cpuFixWallHarmonic(double *x, double *v, double *f, double *eatom, double *vatom, double *fparam, 
        int *iparam, int *ilist, int eflag_atom, int vflag_atom, int dim, int inum)
{
  for (int ii = 0; ii < inum; ii++) {
      int i = ilist[ii]; 
      int direction = iparam[0] / 2;
      int side = iparam[0] % 2;     
      if (side == 0) side = -1;
      double delta,dr,fwall,vn;
      double coord = fparam[0];
      double cutoff = fparam[1];
      double epsilon = fparam[2];
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

void cpuFixWallLJ93(double *x, double *v, double *f, double *eatom, double *vatom, double *fparam, 
        int *iparam, int *ilist, int eflag_atom, int vflag_atom, int dim, int inum)
{
  for (int ii = 0; ii < inum; ii++) {
      int i = ilist[ii]; 
      int direction = iparam[0] / 2;
      int side = iparam[0] % 2;     
      if (side == 0) side = -1;
      double delta,fwall,vn;    
      double coord = fparam[0];
      double cutoff = fparam[1];
      double offset = fparam[2];
      double coeff1 = fparam[3];
      double coeff2 = fparam[4];
      double coeff3 = fparam[5];
      double coeff4 = fparam[6];            
      if (side < 0) delta = x[i*dim+direction] - coord;
      else delta = coord - x[i*dim+direction];
      if (delta >= cutoff) continue;
      if (delta <= 0.0) continue;      
      double rinv = 1.0/delta;
      double r2inv = rinv*rinv;
      double r4inv = r2inv*r2inv;
      double r10inv = r4inv*r4inv*r2inv;
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

void cpuFixWallLJ126(double *x, double *v, double *f, double *eatom, double *vatom, double *fparam, 
        int *iparam, int *ilist, int eflag_atom, int vflag_atom, int dim, int inum)
{
  for (int ii = 0; ii < inum; ii++) {
      int i = ilist[ii];     
      int direction = iparam[0] / 2;
      int side = iparam[0] % 2;     
      if (side == 0) side = -1;
      double delta,fwall,vn;      
      double coord = fparam[0];
      double cutoff = fparam[1];
      double offset = fparam[2];
      double coeff1 = fparam[3];
      double coeff2 = fparam[4];
      double coeff3 = fparam[5];
      double coeff4 = fparam[6];            
      if (side < 0) delta = x[i*dim+direction] - coord;
      else delta = coord - x[i*dim+direction];
      if (delta >= cutoff) continue;
      if (delta <= 0.0) continue;      
      double rinv = 1.0/delta;
      double r2inv = rinv*rinv;
      double r6inv = r2inv*r2inv*r2inv;
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

void cpuFixWallLJ1043(double *x, double *v, double *f, double *eatom, double *vatom, double *fparam, 
        int *iparam, int *ilist, int eflag_atom, int vflag_atom, int dim, int inum)
{
  for (int ii = 0; ii < inum; ii++) {
      int i = ilist[ii]; 
      int direction = iparam[0] / 2;
      int side = iparam[0] % 2;     
      if (side == 0) side = -1;
      double delta,fwall,vn;      
      double coord = fparam[0];
      double cutoff = fparam[1];
      double offset = fparam[2];
      double coeff1 = fparam[3];
      double coeff2 = fparam[4];
      double coeff3 = fparam[5];
      double coeff4 = fparam[6];            
      double coeff5 = fparam[7];            
      double coeff6 = fparam[8];            
      double coeff7 = fparam[9];            
      if (side < 0) delta = x[i*dim+direction] - coord;
      else delta = coord - x[i*dim+direction];
      if (delta >= cutoff) continue;
      if (delta <= 0.0) continue;      
      double rinv = 1.0/delta;
      double r2inv = rinv*rinv;
      double r4inv = r2inv*r2inv;
      double r10inv = r4inv*r4inv*r2inv;
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

void cpuFixWallMorse(double *x, double *v, double *f, double *eatom, double *vatom, double *fparam, 
        int *iparam, int *ilist, int eflag_atom, int vflag_atom, int dim, int inum)
{
  for (int ii = 0; ii < inum; ii++) {
      int i = ilist[ii]; 
      int direction = iparam[0] / 2;
      int side = iparam[0] % 2;     
      if (side == 0) side = -1;
      double delta,dr,fwall,vn, dexp;     
      double coord = fparam[0];
      double cutoff = fparam[1];
      double offset = fparam[2];
      double epsilon = fparam[3];
      double alpha = fparam[4];
      double coeff1 = fparam[5];
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

