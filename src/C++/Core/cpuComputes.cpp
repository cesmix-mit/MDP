/***************************************************************************                               
                    Molecular Dynamics Potentials (MDP)
                           CESMIX-MIT Project  
 
 Contributing authors: Ngoc-Cuong Nguyen (cuongng@mit.edu, exapde@gmail.com)
 ***************************************************************************/

#ifndef MDP_CPUCOMPUTES
#define MDP_CPUCOMPUTES

#include <math.h>
#include <algorithm>

// #define IMGMASK 1023
// #define IMGMAX 512
// #define IMGBITS 10
// #define IMG2BITS 20

// using std::min;
// using std::max;

#define SWAP(a,b) do {       \
    tmp = a; a = b; b = tmp; \
  } while(0)

#define ISWAP(a,b) do {        \
    itmp = a; a = b; b = itmp; \
  } while(0)


template <typename T> void cpuPackIntProperty(T *buf, int *prop, int *ilist, 
         int m, int mvalues, int n, int nvalues, int inum)
{
  for (int ii = 0; ii < inum; ii++) {
    int i = ilist[ii];
    buf[n + nvalues*i] = (T) prop[m + mvalues*i];
  }
}
template void cpuPackIntProperty(double *buf, int *prop, int *ilist, int m, int mvalues, int n, int nvalues, int inum);
template void cpuPackIntProperty(float *buf, int *prop, int *ilist, int m, int mvalues, int n, int nvalues, int inum);

template <typename T> void cpuPackIntProperty(T *buf, int *prop, int *type, int *ilist, 
         int n, int nvalues, int inum)
{
  for (int ii = 0; ii < inum; ii++) {
    int i = ilist[ii]; 
    buf[n + nvalues*i] = (T) prop[type[i]];
  }
}
template void cpuPackIntProperty(double *buf, int *prop, int *type, int *ilist, int n, int nvalues, int inum);
template void cpuPackIntProperty(float *buf, int *prop, int *type, int *ilist, int n, int nvalues, int inum);

template <typename T> void cpuPackFloatProperty(T *buf, T *prop, int *ilist, 
         int m, int mvalues, int n, int nvalues, int inum)
{
  for (int ii = 0; ii < inum; ii++) {
    int i = ilist[ii]; 
    buf[n + nvalues*i] = prop[m + mvalues*i];
  }
}
template void cpuPackFloatProperty(double *buf, double *prop, int *ilist, int m, int mvalues, int n, int nvalues, int inum);
template void cpuPackFloatProperty(float *buf, float *prop, int *ilist, int m, int mvalues, int n, int nvalues, int inum);

template <typename T> void cpuPackFloatProperty(T *buf, T *prop, T a, T b, int *ilist, 
         int m, int mvalues, int n, int nvalues, int inum)
{
  for (int ii = 0; ii < inum; ii++) {
    int i = ilist[ii]; 
    buf[n + nvalues*i] = a*prop[m + mvalues*i] + b;
  }
}
template void cpuPackFloatProperty(double *buf, double *prop, double a, double b, int *ilist, 
         int m, int mvalues, int n, int nvalues, int inum);
template void cpuPackFloatProperty(float *buf, float *prop, float a, float b, int *ilist, 
         int m, int mvalues, int n, int nvalues, int inum);

template <typename T> void cpuPackFloatProperty(T *buf, T *prop, int *type, int *ilist, 
         int n, int nvalues, int inum)
{
  for (int ii = 0; ii < inum; ii++) {
    int i = ilist[ii]; 
    buf[n + nvalues*i] = prop[type[i]];
  }
}
template void cpuPackFloatProperty(double *buf, double *prop, int *type, int *ilist, 
         int n, int nvalues, int inum);
template void cpuPackFloatProperty(float *buf, float *prop, int *type, int *ilist, 
         int n, int nvalues, int inum);

template <typename T> T cpuComputeMass(T *amass, T *mass, int *type, int *ilist, int inum)
{
    T masstotal;
    cpuPackFloatProperty(amass, mass, type, ilist, 0, 1, inum);
    masstotal = cpuArraySum(amass, inum);
    return masstotal;
}         
template double cpuComputeMass(double *amass, double *mass, int *type, int *ilist, int inum);
template float cpuComputeMass(float *amass, float *mass, int *type, int *ilist, int inum);

template <typename T> void  cpuComputeXCM(T *xcm, T *x, T *mass, T *box, T masstotal, 
        int *ilist, int *type, int *image, int triclinic, int dim, int inum)
{
  xcm[0] = xcm[1] = xcm[2] = 0.0;
  T unwrap[3];
  for (int ii = 0; ii < inum; ii++) {
        int i = ilist[ii]; 
        T massone = mass[type[i]];
        cpuUnmap(unwrap, &x[i*dim], box, &image[i*dim], triclinic, dim);        
        xcm[0] += unwrap[0] * massone;
        xcm[1] += unwrap[1] * massone;
        if (dim==3) xcm[2] += unwrap[2] * massone;      
  }

  //MPI_Allreduce(cmone,xcm,3,MPI_DOUBLE,MPI_SUM,world);
  if (masstotal > 0.0) {
    xcm[0] /= masstotal;
    xcm[1] /= masstotal;
    if (dim==3) xcm[2] /= masstotal;
  }
}
template void  cpuComputeXCM(double *xcm, double *x, double *mass, double *box, double masstotal, 
        int *ilist, int *type, int *image, int triclinic, int dim, int inum);
template void  cpuComputeXCM(float *xcm, float *x, float *mass, float *box, float masstotal, 
        int *ilist, int *type, int *image, int triclinic, int dim, int inum);

template <typename T> void  cpuComputeVCM(T *vcm, T *v, T *mass, T masstotal, 
        int *ilist, int *type, int dim, int inum)
{  
    vcm[0] = vcm[1] = vcm[2] = 0.0;

    for (int ii = 0; ii < inum; ii++) {
        int i = ilist[ii]; 
        T massone = mass[type[i]];
        vcm[0] += v[i*dim+0]*massone;
        vcm[1] += v[i*dim+1]*massone;
        if (dim==3) vcm[2] += v[i*dim+2]*massone;
    }
  
  //MPI_Allreduce(p,xcm,3,MPI_DOUBLE,MPI_SUM,world);
  if (masstotal > 0.0) {
    vcm[0] /= masstotal;
    vcm[1] /= masstotal;
    if (dim==3) vcm[2] /= masstotal;
  }
}
template void  cpuComputeVCM(double *vcm, double *v, double *mass, double masstotal, 
        int *ilist, int *type, int dim, int inum);
template void  cpuComputeVCM(float *vcm, float *v, float *mass, float masstotal, 
        int *ilist, int *type, int dim, int inum);

template <typename T> T cpuComputeGyration(T *xcm, T *x, T *mass, T *box, T masstotal, 
        int *ilist, int *type, int *image, int triclinic, int dim, int inum)
{
  T unwrap[3];
  T rg = 0.0;

  for (int ii = 0; ii < inum; ii++) {
        int i = ilist[ii]; 
        T massone = mass[type[i]];
        cpuUnmap(unwrap, &x[i*dim], box, &image[i*dim], triclinic, dim);        
        T dx = unwrap[0] - xcm[0];
        T dy = unwrap[1] - xcm[1];        
        T dz = (dim==3) ? (unwrap[2] - xcm[2]) : 0.0;
        rg += (dx*dx + dy*dy + dz*dz) * massone;
  }
  
  if (masstotal > 0.0) return sqrt(rg/masstotal);
  return 0.0;
}
template double cpuComputeGyration(double *xcm, double *x, double *mass, double *box, double masstotal, 
        int *ilist, int *type, int *image, int triclinic, int dim, int inum);
template float  cpuComputeGyration(float *xcm, float *x, float *mass, float *box, float masstotal, 
        int *ilist, int *type, int *image, int triclinic, int dim, int inum);

template <typename T> void cpuComputeAngmom(T *p, T *xcm, T *x, T *v, T *mass, T *box, 
        int *ilist, int *type, int *image, int triclinic, int dim, int inum)
{
  T unwrap[3];  
  p[0] = p[1] = p[2] = 0.0;
  for (int ii = 0; ii < inum; ii++) {
        int i = ilist[ii]; 
        T massone = mass[type[i]];
        cpuUnmap(unwrap, &x[i*dim], box, &image[i*dim], triclinic, dim);        
        T dx = unwrap[0] - xcm[0];
        T dy = unwrap[1] - xcm[1];        
        p[2] += massone * (dx*v[i*dim+1] - dy*v[i*dim+0]);
        if (dim==3) {
            T dz = unwrap[2] - xcm[2];
            p[0] += massone * (dy*v[i*dim+2] - dz*v[i*dim+1]);
            p[1] += massone * (dz*v[i*dim+0] - dx*v[i*dim+2]);
        }
  }   
  //MPI_Allreduce(p,lmom,3,MPI_DOUBLE,MPI_SUM,world);
}
template void cpuComputeAngmom(double *p, double *xcm, double *x, double *v, double *mass, double *box, 
        int *ilist, int *type, int *image, int triclinic, int dim, int inum);
template void cpuComputeAngmom(float *p, float *xcm, float *x, float *v, float *mass, float *box, 
        int *ilist, int *type, int *image, int triclinic, int dim, int inum);

template <typename T> void cpuComputeTorque(T *tlocal, T *xcm, T *x, T *f, T *box, 
        int *ilist, int *image, int triclinic, int dim, int inum)
{
  T unwrap[3];
  tlocal[0] = tlocal[1] = tlocal[2] = 0.0;
  for (int ii = 0; ii < inum; ii++) {
        int i = ilist[ii]; 
        cpuUnmap(unwrap, &x[i*dim], box, &image[i*dim], triclinic, dim);        
        T dx = unwrap[0] - xcm[0];
        T dy = unwrap[1] - xcm[1];        
        tlocal[2] += dx*f[i*dim+1] - dy*f[i*dim+0];
        if (dim==3) {
            T dz = unwrap[2] - xcm[2];
            tlocal[0] += dy*f[i*dim+2] - dz*f[i*dim+1];
            tlocal[1] += dz*f[i*dim+0] - dx*f[i*dim+2];
        }
  }     
  //MPI_Allreduce(tlocal,tq,3,MPI_DOUBLE,MPI_SUM,world);
}
template void cpuComputeTorque(double *tlocal, double *xcm, double *x, double *f, double *box, 
        int *ilist, int *image, int triclinic, int dim, int inum);
template void cpuComputeTorque(float *tlocal, float *xcm, float *x, float *f, float *box, 
        int *ilist, int *image, int triclinic, int dim, int inum);

template <typename T> void cpuComputeInertia(T *ione, T *xcm, T *x, T *mass, T *box, 
        int *ilist, int *type, int *image, int triclinic, int dim, int inum)
{
  
  T unwrap[3];
  for (int i = 0; i < dim; i++)
    for (int j = 0; j < dim; j++)
      ione[i*dim+j] = 0.0;

  for (int ii = 0; ii < inum; ii++) {
        int i = ilist[ii]; 
        T massone = mass[type[i]];
        cpuUnmap(unwrap, &x[i*dim], box, &image[i*dim], triclinic, dim);        
        T dx = unwrap[0] - xcm[0];
        T dy = unwrap[1] - xcm[1];                        
        ione[0] += massone * (dy*dy);
        ione[1] += massone * (dx*dx);        
        ione[2] += massone * (dx*dx + dy*dy);
        ione[3] -= massone * dx*dy;
        if (dim==3) {
            T dz = unwrap[2] - xcm[2];
            ione[0] += massone * (dz*dz);
            ione[1] += massone * (dz*dz);
            ione[4] -= massone * dy*dz;
            ione[5] -= massone * dx*dz;
        }
  }    
//   ione[1*dim+0] = ione[0*dim+1];
//   ione[2*dim+1] = ione[1*dim+2];
//   ione[2*dim+0] = ione[0*dim+2];
  //MPI_Allreduce(&ione[0],&itensor[0],6,MPI_DOUBLE,MPI_SUM,world);
}
template void cpuComputeInertia(double *ione, double *xcm, double *x, double *mass, double *box, 
        int *ilist, int *type, int *image, int triclinic, int dim, int inum);
template void cpuComputeInertia(float *ione, float *xcm, float *x, float *mass, float *box, 
        int *ilist, int *type, int *image, int triclinic, int dim, int inum);

template <typename T> T cpuComputeNVTEnergy(T *tarray, T *eta, T *eta_mass, T *eta_dot, int mtchain)
{
    T t_target = tarray[3]; 
    T tdof = tarray[5];
    T boltz = tarray[6];  
    T kt = boltz * t_target;
    T ke_target = tdof * kt;        
    T energy = ke_target * eta[0] + 0.5*eta_mass[0]*eta_dot[0]*eta_dot[0];
    for (int ich = 1; ich < mtchain; ich++)
        energy += kt * eta[ich] + 0.5*eta_mass[ich]*eta_dot[ich]*eta_dot[ich];        
    return energy;
}
template double cpuComputeNVTEnergy(double *tarray, double *eta, double *eta_mass, double *eta_dot, int mtchain);
template float cpuComputeNVTEnergy(float *tarray, float *eta, float *eta_mass, float *eta_dot, int mtchain);

template <typename T> T cpuComputeTempScalar(T *v, T *mass, T tfactor, int *type, int *ilist, 
         int dim, int inum)
{
    T t = 0.0;
    if (dim==3) {    
        for (int ii = 0; ii < inum; ii++) {
            int i = ilist[ii];
            t += (v[0+dim*i]*v[0+dim*i] + v[1+dim*i]*v[1+dim*i] + v[2+dim*i]*v[2+dim*i]) * mass[type[i]];
        }
    }
    else if (dim==2) {
        for (int ii = 0; ii < inum; ii++) {
            int i = ilist[ii];
            t += (v[0+dim*i]*v[0+dim*i] + v[1+dim*i]*v[1+dim*i]) * mass[type[i]];        
        }
    }
    
  t = tfactor*t;   
  
  return t;    
}
template double cpuComputeTempScalar(double *v, double *mass, double tfactor, int *type, int *ilist, 
         int dim, int inum);
template float cpuComputeTempScalar(float *v, float *mass, float tfactor, int *type, int *ilist, 
         int dim, int inum);

template <typename T> void cpuComputeTempSymTensor(T *t, T *v, T *mass, T tfactor, int *type, int *ilist, 
         int dim, int inum)
{
    for (int i = 0; i < 6; i++) t[i] = 0.0;
    
    if (dim==3) {    
        for (int ii = 0; ii < inum; ii++) {
          int i = ilist[ii]; 
          T massone = mass[type[i]];
          t[0] += massone * v[0+dim*i]*v[0+dim*i];
          t[1] += massone * v[1+dim*i]*v[1+dim*i];
          t[2] += massone * v[2+dim*i]*v[2+dim*i];
          t[3] += massone * v[0+dim*i]*v[1+dim*i];
          t[4] += massone * v[0+dim*i]*v[2+dim*i];
          t[5] += massone * v[1+dim*i]*v[2+dim*i];           
        }            
    }
    else if (dim==2) {
        for (int ii = 0; ii < inum; ii++) {
          int i = ilist[ii]; 
          T massone = mass[type[i]];
          t[0] += massone * v[0+dim*i]*v[0+dim*i];
          t[1] += massone * v[1+dim*i]*v[1+dim*i];
          t[3] += massone * v[0+dim*i]*v[1+dim*i];
        }
    }
  
  for (int i = 0; i < 6; i++) t[i] = tfactor*t[i];      
}
template void cpuComputeTempSymTensor(double *t, double *v, double *mass, double tfactor, int *type, int *ilist, 
         int dim, int inum);
template void cpuComputeTempSymTensor(float *t, float *v, float *mass, float tfactor, int *type, int *ilist, 
         int dim, int inum);

template <typename T> T cpuComputePressureScalar(T *virial, T volume, T temp, T tempdof, 
        T boltz, T nktv2p, int dim)
{
  T factor = nktv2p /(dim*volume);
  if (dim == 3) {
      return (tempdof * boltz * temp + virial[0] + virial[1] + virial[2]) * factor;
  } else {
      return (tempdof * boltz * temp + virial[0] + virial[1]) *factor;
  }  
}
template double cpuComputePressureScalar(double *virial, double volume, double temp, double tempdof, 
        double boltz, double nktv2p, int dim);
template float cpuComputePressureScalar(float *virial, float volume, float temp, float tempdof, 
        float boltz, float nktv2p, int dim);

template <typename T> void cpuComputePressureSymTensor(T *vector, T *virial, T *ke_tensor, 
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
template void cpuComputePressureSymTensor(double *vector, double *virial, double *ke_tensor, 
        double volume, double nktv2p, int dim);
template void cpuComputePressureSymTensor(float *vector, float *virial, float *ke_tensor, 
        float volume, float nktv2p, int dim);

template <typename T> void cpuComputeHeatFlux(T *vector, T *ke, T *pe, T *stress, T *v, 
        T nktv2p, int *ilist,  int pressatomflag, int dim, int inum)
{   
  T jc[3] = {0.0,0.0,0.0};
  T jv[3] = {0.0,0.0,0.0};
  T eng;

  if (pressatomflag == 2) {
    for (int ii = 0; ii < inum; ii++) {
        int i = ilist[ii]; 
        eng = pe[i] + ke[i];
        jc[0] += eng*v[i*dim+0];
        jc[1] += eng*v[i*dim+1];
        jc[2] += eng*v[i*dim+2];
        // stress[0]: rijx*fijx
        // stress[1]: rijy*fijy
        // stress[2]: rijz*fijz
        // stress[3]: rijx*fijy
        // stress[4]: rijx*fijz
        // stress[5]: rijy*fijz
        // stress[6]: rijy*fijx
        // stress[7]: rijz*fijx
        // stress[8]: rijz*fijy
        // jv[0]  = rijx fijx vjx + rijx fijy vjy + rijx fijz vjz
        jv[0] -= stress[i*9+0]*v[i*dim+0] + stress[i*9+3]*v[i*dim+1] +
          stress[i*9+4]*v[i*dim+2];
        // jv[1]  = rijy fijx vjx + rijy fijy vjy + rijy fijz vjz
        jv[1] -= stress[i*9+6]*v[i*dim+0] + stress[i*9+1]*v[i*dim+1] +
          stress[i*9+5]*v[i*dim+2];
        // jv[2]  = rijz fijx vjx + rijz fijy vjy + rijz fijz vjz
        jv[2] -= stress[i*9+7]*v[i*dim+0] + stress[i*9+8]*v[i*dim+1] +
          stress[i*9+2]*v[i*dim+2];      
    }
  } else {
    for (int ii = 0; ii < inum; ii++) {
        int i = ilist[ii]; 
        eng = pe[i] + ke[i];
        jc[0] += eng*v[i*dim+0];
        jc[1] += eng*v[i*dim+1];
        jc[2] += eng*v[i*dim+2];
        jv[0] -= stress[i*6+0]*v[i*dim+0] + stress[i*6+3]*v[i*dim+1] +
          stress[i*6+4]*v[i*dim+2];
        jv[1] -= stress[i*6+3]*v[i*dim+0] + stress[i*6+1]*v[i*dim+1] +
          stress[i*6+5]*v[i*dim+2];
        jv[2] -= stress[i*6+4]*v[i*dim+0] + stress[i*6+5]*v[i*dim+1] +
          stress[i*6+2]*v[i*dim+2];      
    }
  }

  // convert jv from stress*volume to energy units via nktv2p factor
  jv[0] /= nktv2p;
  jv[1] /= nktv2p;
  jv[2] /= nktv2p;
    
  for (int i=0; i<3; i++) {
    vector[i] = jc[i]+jv[i]; // 1st 3 terms are total heat flux
    vector[i+3] = jc[i]; // 2nd 3 terms are just conductive portion
  }  
}
template void cpuComputeHeatFlux(double *vector, double *ke, double *pe, double *stress, double *v, 
        double nktv2p, int *ilist,  int pressatomflag, int dim, int inum);
template void cpuComputeHeatFlux(float *vector, float *ke, float *pe, float *stress, float *v, 
        float nktv2p, int *ilist,  int pressatomflag, int dim, int inum);

template <typename T> void cpuComputeKEAtom(T *ke, T *mass,  T *v, 
        T mvv2e, int *type, int *ilist,  int dim, int inum)
{    
    if (dim==3) {
        for (int ii = 0; ii < inum; ii++) {
            int i = ilist[ii]; 
            ke[ii] = 0.5 * mvv2e * mass[type[i]] *
              (v[i*dim+0]*v[i*dim+0] + v[i*dim+1]*v[i*dim+1] + v[i*dim+2]*v[i*dim+2]);        
        }
    } else {
        for (int ii = 0; ii < inum; ii++) {
            int i = ilist[ii]; 
            ke[ii] = 0.5 * mvv2e * mass[type[i]] *
              (v[i*dim+0]*v[i*dim+0] + v[i*dim+1]*v[i*dim+1]);        
        }        
    }       
}
template void cpuComputeKEAtom(double *ke, double *mass,  double *v, 
        double mvv2e, int *type, int *ilist,  int dim, int inum);
template void cpuComputeKEAtom(float *ke, float *mass,  float *v, 
        float mvv2e, int *type, int *ilist,  int dim, int inum);

template <typename T> void cpuComputeStressAtom(T *stress, T *mass, T *vatom, T *v, 
        T mvv2e, T nktv2p, int *type, int *ilist,  int vflag, int keflag, int dim, int inum)
{
  if (vflag) {  
      for (int ii = 0; ii < inum; ii++)
      { 
        int i = ilist[ii];   
        for (int j = 0; j < 6; j++)
          stress[ii+j*inum] = vatom[i*6+j];
      }
  }
  else {
      cpuArraySetValue(stress, (T) 0.0, inum*6);
  }
  
  // zero virial of atoms not in group
  // only do this after comm since ghost contributions must be included
//   for (int ii = 0; ii < inum; ii++)
//   {
//       int i = ilist[ii];   
//       stress[i*6+0] = 0.0;
//       stress[i*6+1] = 0.0;
//       stress[i*6+2] = 0.0;
//       stress[i*6+3] = 0.0;
//       stress[i*6+4] = 0.0;
//       stress[i*6+5] = 0.0;
//   }

  // include kinetic energy term for each atom in group
  // mvv2e converts mv^2 to energy
  if (keflag) {
    for (int ii = 0; ii < inum; ii++)
    {
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
//         stress[i*6+0] += onemass*v[i*dim+0]*v[i*dim+0];
//         stress[i*6+1] += onemass*v[i*dim+1]*v[i*dim+1];
//         stress[i*6+3] += onemass*v[i*dim+0]*v[i*dim+1];
//         if (dim==3) {
//         stress[i*6+2] += onemass*v[i*dim+2]*v[i*dim+2];
//         stress[i*6+4] += onemass*v[i*dim+0]*v[i*dim+2];
//         stress[i*6+5] += onemass*v[i*dim+1]*v[i*dim+2];
//         }
    }      
  }

  cpuArrayMultiplyScalar(stress, nktv2p, inum*6);    
}
template void cpuComputeStressAtom(double *stress, double *mass, double *vatom, double *v, 
        double mvv2e, double nktv2p, int *type, int *ilist, int vflag, int keflag, int dim, int inum);
template void cpuComputeStressAtom(float *stress, float *mass, float *vatom, float *v, 
        float mvv2e, float nktv2p, int *type, int *ilist, int vflag, int keflag, int dim, int inum);

template <typename T> void cpuComputeCentroidStressAtom(T *stress, T *mass, T *cvatom, T *v, 
        T mvv2e, T nktv2p, int *type, int *ilist, int vflag, int keflag, int dim, int inum)
{
  if (vflag) {
      for (int ii = 0; ii < inum; ii++)
      {
        int i = ilist[ii];   
        for (int j = 0; j < 9; j++)
          stress[inum*j+ii] = cvatom[i*9+j];
      }
  } 
  else {
      cpuArraySetValue(stress, (T) 0.0, inum*9);
  }
  
  // zero virial of atoms not in group
  // only do this after comm since ghost contributions must be included
//   for (int ii = 0; ii < inum; ii++)
//   {
//       int i = ilist[ii];   
//       stress[i*9+0] = 0.0;
//       stress[i*9+1] = 0.0;
//       stress[i*9+2] = 0.0;
//       stress[i*9+3] = 0.0;
//       stress[i*9+4] = 0.0;
//       stress[i*9+5] = 0.0;
//       stress[i*9+6] = 0.0;
//       stress[i*9+7] = 0.0;
//       stress[i*9+8] = 0.0;
//   }    

  // include kinetic energy term for each atom in group
  // apply temperature bias is applicable
  // mvv2e converts mv^2 to energy
  if (keflag) {
    for (int ii = 0; ii < inum; ii++)
    {
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
//         stress[i*9+0] += onemass*v[i*dim+0]*v[i*dim+0];
//         stress[i*9+1] += onemass*v[i*dim+1]*v[i*dim+1];        
//         stress[i*9+3] += onemass*v[i*dim+0]*v[i*dim+1];        
//         stress[i*9+6] += onemass*v[i*dim+1]*v[i*dim+0];
//         if (dim==3) {
//         stress[i*9+2] += onemass*v[i*dim+2]*v[i*dim+2];
//         stress[i*9+4] += onemass*v[i*dim+0]*v[i*dim+2];
//         stress[i*9+5] += onemass*v[i*dim+1]*v[i*dim+2];
//         stress[i*9+7] += onemass*v[i*dim+2]*v[i*dim+0];
//         stress[i*9+8] += onemass*v[i*dim+2]*v[i*dim+1];
//         }
    }          
  }

  // convert to stress*volume units = -pressure*volume
  cpuArrayMultiplyScalar(stress, nktv2p, inum*9);
//   for (int ii = 0; ii < inum; ii++)
//   {
//       int i = ilist[ii];   
//       stress[i*9+0] *= nktv2p;
//       stress[i*9+1] *= nktv2p;
//       stress[i*9+2] *= nktv2p;
//       stress[i*9+3] *= nktv2p;
//       stress[i*9+4] *= nktv2p;
//       stress[i*9+5] *= nktv2p;
//       stress[i*9+6] *= nktv2p;
//       stress[i*9+7] *= nktv2p;
//       stress[i*9+8] *= nktv2p;
//     }
}
template void cpuComputeCentroidStressAtom(double *stress, double *mass, double *cvatom, double *v, 
        double mvv2e, double nktv2p, int *type, int *ilist, int vflag, int keflag, int dim, int inum);
template void cpuComputeCentroidStressAtom(float *stress, float *mass, float *cvatom, float *v, 
        float mvv2e, float nktv2p, int *type, int *ilist, int vflag, int keflag, int dim, int inum);

template <typename T> void cpuComputeDisplaceAtom(T *displace, T *x, T *xoriginal, T *box, 
        int *pbc, int *ilist,  int triclinic, int dim, int inum)
{
    for (int ii = 0; ii < inum; ii++)
    {
        int i = ilist[ii];   
        T dp[3];
        dp[0] = x[i*dim+0] - xoriginal[i*dim+0];
        dp[1] = x[i*dim+1] - xoriginal[i*dim+1];        
        dp[2] = (dim==3) ? (x[i*dim+2] - xoriginal[i*dim+2]) : 0.0;
        cpuMinimumImage(dp, box, pbc, triclinic, dim);
        displace[ii+0*inum] = dp[0];
        displace[ii+1*inum] = dp[1];
        displace[ii+2*inum] = dp[2];        
        displace[ii+3*inum] = sqrt(dp[0]*dp[0] +dp[1]*dp[1] + dp[2]*dp[2]);        
//         displace[i*4+0] = x[i*dim+0] - xoriginal[i*dim+0];
//         displace[i*4+1] = x[i*dim+1] - xoriginal[i*dim+1];
//         displace[i*4+2] = x[i*dim+2] - xoriginal[i*dim+2];
//         cpuMinimumImage(&displace[i*4], box, pbc, triclinic, dim);
//         displace[i*4+3] = sqrt(displace[i*4+0]*displace[i*4+0] +
//                                displace[i*4+1]*displace[i*4+1] + 
//                                displace[i*4+2]*displace[i*4+2]);
//       } else {
//           displace[i*4+0] = 0.0;
//           displace[i*4+1] = 0.0;
//           displace[i*4+2] = 0.0;
//           displace[i*4+3] = 0.0;  
    }
}
template void cpuComputeDisplaceAtom(double *displace, double *x, double *xoriginal, double *box, 
        int *pbc, int *ilist,  int triclinic, int dim, int inum);
template void cpuComputeDisplaceAtom(float *displace, float *x, float *xoriginal, float *box, 
        int *pbc, int *ilist,  int triclinic, int dim, int inum);

template <typename T> void cpuSelect3(int k, int n, T *arr, int *iarr, T *arr3)
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
//       arr3[l+1][0] = arr3[j][0];
//       arr3[l+1][1] = arr3[j][1];
//       arr3[l+1][2] = arr3[j][2];
//       arr3[j][0] = a3[0];
//       arr3[j][1] = a3[1];
//       arr3[j][2] = a3[2];
      if (j >= k) ir = j-1;
      if (j <= k) l = i;
    }
  }
}

template <typename T> T cpuAssociatedLegendre(int l, int m, T x)
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

template <typename T> T cpuPolarPrefactor(int l, int m, T MY_4PI, T costheta)
{
  const int mabs = abs(m);

  T prefactor = 1.0;
  for (int i=l-mabs+1; i < l+mabs+1; ++i)
    prefactor *= (T) i;

  prefactor = sqrt(static_cast<T>(2*l+1)/(MY_4PI*prefactor))
    * cpuAssociatedLegendre(l,mabs,costheta);

  if ((m < 0) && (m % 2)) prefactor = -prefactor;

  return prefactor;
}

/* ----------------------------------------------------------------------
   associated legendre polynomial
   sign convention: P(l,l) = (2l-1)!!(-sqrt(1-x^2))^l
------------------------------------------------------------------------- */

template <typename T> void cpuCalcBoop(T *qn, T *rlist, T *cglist, T *qnm_r, T *qnm_i,
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

      qnm_r[il*mmax+l] += cpuPolarPrefactor(l, 0, MY_4PI, costheta);
      T expphim_r = expphi_r;
      T expphim_i = expphi_i;
      for(int m = 1; m <= +l; m++) {

        T prefactor = cpuPolarPrefactor(l, m, MY_4PI, costheta);
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


template <typename T> int cpuInitClebschGordan(T *cglist, T *factorial, int *qlist, int nqlist)
{
  T sum,dcg,sfaccg, sfac1, sfac2;
  int m, aa2, bb2, cc2;
  int ifac, idxcg_count;

  idxcg_count = 0;
  for (int il = 0; il < nqlist; il++) {
    int l = qlist[il];
    for(int m1 = 0; m1 < 2*l+1; m1++)
      for(int m2 = max(0,l-m1); m2 < min(2*l+1,3*l-m1+1); m2++)
        idxcg_count++;
  }
  int idxcg_max = idxcg_count;

  idxcg_count = 0;
  for (int il = 0; il < nqlist; il++) {
    int l = qlist[il];
    for(int m1 = 0; m1 < 2*l+1; m1++) {
        aa2 = m1 - l;
        for(int m2 = max(0,l-m1); m2 < min(2*l+1,3*l-m1+1); m2++) {
          bb2 = m2 - l;
          m = aa2 + bb2 + l;

          sum = 0.0;
          for (int z = max(0, max(-aa2, bb2));
               z <= min(l, min(l - aa2, l + bb2)); z++) {
            ifac = z % 2 ? -1 : 1;
            sum += ifac /
              (factorial[z] *
               factorial[l - z] *
               factorial[l - aa2 - z] *
               factorial[l + bb2 - z] *
               factorial[aa2 + z] *
               factorial[-bb2 + z]);
          }

          cc2 = m - l;
          sfaccg = sqrt(factorial[l + aa2] *
                        factorial[l - aa2] *
                        factorial[l + bb2] *
                        factorial[l - bb2] *
                        factorial[l + cc2] *
                        factorial[l - cc2] *
                        (2*l + 1));

          sfac1 = factorial[3*l + 1];
          sfac2 = factorial[l];
          dcg = sqrt(sfac2*sfac2*sfac2 / sfac1);

          cglist[idxcg_count] = sum * dcg * sfaccg;
          idxcg_count++;
        }
      }
  }
  
  return idxcg_max;
}

template <typename T> void cpuComputeOrientOrderAtom(T* qnarray, T *x, T *rlist, T *cglist, T *fac, 
        T *qnm_r, T *qnm_i, T *distsq, T cutsq, T MY_EPSILON, T QEPSILON, T MY_4PI, int *neighlist, 
        int *neighnum,  int *ilist, int *qlist, int *nearest,  int nqlist, int qmax, int wlflag,  
        int wlhatflag,  int qlcompflag, int iqlcomp, int qlcomp, int nnn, int jnum, int dim, int inum) 
{
    
  int ncol = nqlist;
  if (wlflag) ncol += nqlist;
  if (wlhatflag) ncol += nqlist;
  if (qlcompflag) ncol += 2*(2*qlcomp+1);    
  
  cpuInitClebschGordan(cglist, fac, qlist, nqlist);
  
  for (int ii = 0; ii < inum; ii++) {
      int i = ilist[ii];
      T* qn = &qnarray[ncol*i];    
      T xtmp = x[i*dim+0];
      T ytmp = x[i*dim+1];
      T ztmp = x[i*dim+2];
      int *jlist = &neighlist[jnum*i];
      int m = neighnum[i];
      
//         memory->create(distsq,maxneigh,"orientorder/atom:distsq");
//         memory->create(rlist,maxneigh,3,"orientorder/atom:rlist");
//         memory->create(nearest,maxneigh,"orientorder/atom:nearest");
      
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
        cpuSelect3(nnn,ncount,distsq,nearest,rlist);
        ncount = nnn;
      }
      
      cpuCalcBoop(qn, rlist, cglist, qnm_r, qnm_i, MY_EPSILON, QEPSILON, MY_4PI,
              qlist, nqlist, qmax, wlflag, wlhatflag, qlcompflag, iqlcomp, qlcomp, ncount);
          
  }
}
template void cpuComputeOrientOrderAtom(double *qnarray, double *x, double *rlist, double *cglist, double *fac, 
        double *qnm_r, double *qnm_i, double *distsq, double cutsq, double MY_EPSILON, double QEPSILON, double MY_4PI, int *neighlist, 
        int *neighnum,  int *ilist, int *qlist, int *nearest,  int nqlist, int qmax, int wlflag,  
        int wlhatflag,  int qlcompflag, int iqlcomp, int qlcomp, int nnn, int jnum, int dim, int inum); 
template void cpuComputeOrientOrderAtom(float *qnarray, float *x, float *rlist, float *cglist, float *fac, 
        float *qnm_r, float *qnm_i, float *distsq, float cutsq, float MY_EPSILON, float QEPSILON, float MY_4PI, int *neighlist, 
        int *neighnum,  int *ilist, int *qlist, int *nearest,  int nqlist, int qmax, int wlflag,  
        int wlhatflag,  int qlcompflag, int iqlcomp, int qlcomp, int nnn, int jnum, int dim, int inum); 

template <typename T> void cpuComputeCoordAtomCutoff(int *cvec, T *x, T *rcutsq, int *type, 
        int *ilist, int *neighlist, int *neighnum, int *typelo, int *typehi, 
         int *jgroupbit, int dim, int ntypes, int jnum, int inum)
{
  for (int ii = 0; ii < inum; ii++) {
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
  }
}
template void cpuComputeCoordAtomCutoff(int *cvec, double *x, double *rcutsq, int *type, 
        int *ilist, int *neighlist, int *neighnum, int *typelo, int *typehi, 
         int *jgroupbit, int dim, int ntypes, int jnum, int inum);
template void cpuComputeCoordAtomCutoff(int *cvec, float *x, float *rcutsq, int *type, 
        int *ilist, int *neighlist, int *neighnum, int *typelo, int *typehi, 
         int *jgroupbit, int dim, int ntypes, int jnum, int inum);

template <typename T> void cpuComputeCoordAtomCutoff(int *carray, T *x, T *rcutsq, int *type, 
        int *ilist, int *neighlist, int *neighnum, int *typelo, int *typehi, 
         int *jgroupbit, int ncol, int dim, int ntypes, int jnum, int inum)
{
  for (int ii = 0; ii < inum; ii++) {
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
  }
}
template void cpuComputeCoordAtomCutoff(int *carray, double *x, double *rcutsq, int *type, 
        int *ilist, int *neighlist, int *neighnum, int *typelo, int *typehi, 
         int *jgroupbit, int ncol, int dim, int ntypes, int jnum, int inum);
template void cpuComputeCoordAtomCutoff(int *carray, float *x, float *rcutsq, int *type, 
        int *ilist, int *neighlist, int *neighnum, int *typelo, int *typehi, 
         int *jgroupbit, int ncol, int dim, int ntypes, int jnum, int inum);

template <typename T> void cpuComputeCoordAtomOrient(int *cvec, T *x, T *rcutsq, T *normv, T threshold, 
        int *type, int *ilist, int *neighlist, int *neighnum, int *typelo, int *typehi, 
         int *jgroupbit, int dim, int ntypes, int nqlist, int ncol, int l, int jnum, int inum)
{
  for (int ii = 0; ii < inum; ii++) {
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
  }
}
template void cpuComputeCoordAtomOrient(int *cvec, double *x, double *rcutsq, double *normv, double threshold, 
        int *type, int *ilist, int *neighlist, int *neighnum, int *typelo, int *typehi, 
         int *jgroupbit, int dim, int ntypes, int nqlist, int ncol, int l, int jnum, int inum);
template void cpuComputeCoordAtomOrient(int *cvec, float *x, float *rcutsq, float *normv, float threshold, 
        int *type, int *ilist, int *neighlist, int *neighnum, int *typelo, int *typehi, 
         int *jgroupbit, int dim, int ntypes, int nqlist, int ncol, int l, int jnum, int inum);

template <typename T> void cpuComputeMSD(T *vector, T *x, T *xoriginal, T *h, T *xcm,
         int *ilist, int *image, int naverage, int avflag, int triclinic, int nmsd, int dim,  int inum)
{
  T dx,dy,dz;
  int xbox,ybox,zbox;

  T msd[4];
  msd[0] = msd[1] = msd[2] = msd[3] = 0.0;

  T xtmp, ytmp, ztmp;

  // update number of averages if requested
  T navfac;  
  if (avflag) {
    naverage += 1;
    navfac = 1.0/(naverage+1);
  }

  if (dim==3) {
      if (triclinic == 0) {
        for (int ii = 0; ii < inum; ii++) {
            int i = ilist[ii];
//             xbox = (image[i] & IMGMASK) - IMGMAX;
//             ybox = (image[i] >> IMGBITS & IMGMASK) - IMGMAX;
//             zbox = (image[i] >> IMG2BITS) - IMGMAX;
            xtmp = x[i*dim+0] + image[i*dim+0]*h[0] - xcm[0];
            ytmp = x[i*dim+1] + image[i*dim+1]*h[1] - xcm[1];
            ztmp = x[i*dim+2] + image[i*dim+2]*h[2] - xcm[2];
            if (avflag) {
              xoriginal[i*dim+0] = (xoriginal[i*dim+0]*naverage + xtmp)*navfac;
              xoriginal[i*dim+1] = (xoriginal[i*dim+1]*naverage + ytmp)*navfac;
              xoriginal[i*dim+2] = (xoriginal[i*dim+2]*naverage + ztmp)*navfac;
            }

            dx = xtmp - xoriginal[i*dim+0];
            dy = ytmp - xoriginal[i*dim+1];
            dz = ztmp - xoriginal[i*dim+2];
            msd[0] += dx*dx;
            msd[1] += dy*dy;
            msd[2] += dz*dz;
            msd[3] += dx*dx + dy*dy + dz*dz;
          }
      } else {
        for (int ii = 0; ii < inum; ii++) {
            int i = ilist[ii];
//             xbox = (image[i] & IMGMASK) - IMGMAX;
//             ybox = (image[i] >> IMGBITS & IMGMASK) - IMGMAX;
//             zbox = (image[i] >> IMG2BITS) - IMGMAX;
            xtmp = x[i*dim+0] + h[0]*image[i*dim+0] + h[5]*image[i*dim+1] + h[4]*image[i*dim+2] - xcm[0];
            ytmp = x[i*dim+1] + h[1]*image[i*dim+1] + h[3]*zbox - xcm[1];
            ztmp = x[i*dim+2] + h[2]*image[i*dim+2] - xcm[2];

            // use running average position for reference if requested

            if (avflag) {
              xoriginal[i*dim+0] = (xoriginal[i*dim+0]*naverage + xtmp)*navfac;
              xoriginal[i*dim+1] = (xoriginal[i*dim+0]*naverage + xtmp)*navfac;
              xoriginal[i*dim+2] = (xoriginal[i*dim+0]*naverage + xtmp)*navfac;
            }

            dx = xtmp - xoriginal[i*dim+0];
            dy = ytmp - xoriginal[i*dim+1];
            dz = ztmp - xoriginal[i*dim+2];
            msd[0] += dx*dx;
            msd[1] += dy*dy;
            msd[2] += dz*dz;
            msd[3] += dx*dx + dy*dy + dz*dz;
          }
      }      
  } else {
      if (triclinic == 0) {
        for (int ii = 0; ii < inum; ii++) {
            int i = ilist[ii];
//             xbox = (image[i] & IMGMASK) - IMGMAX;
//             ybox = (image[i] >> IMGBITS & IMGMASK) - IMGMAX;
            xtmp = x[i*dim+0] + image[i*dim+0]*h[0] - xcm[0];
            ytmp = x[i*dim+1] + image[i*dim+1]*h[1] - xcm[1];
            if (avflag) {
              xoriginal[i*dim+0] = (xoriginal[i*dim+0]*naverage + xtmp)*navfac;
              xoriginal[i*dim+1] = (xoriginal[i*dim+1]*naverage + ytmp)*navfac;
            }
            dx = xtmp - xoriginal[i*dim+0];
            dy = ytmp - xoriginal[i*dim+1];
            msd[0] += dx*dx;
            msd[1] += dy*dy;
            msd[3] += dx*dx + dy*dy;
          }
      } else {
        for (int ii = 0; ii < inum; ii++) {
            int i = ilist[ii];
//             xbox = (image[i] & IMGMASK) - IMGMAX;
//             ybox = (image[i] >> IMGBITS & IMGMASK) - IMGMAX;
            xtmp = x[i*dim+0] + h[0]*image[i*dim+0] + h[5]*image[i*dim+1] - xcm[0];
            ytmp = x[i*dim+1] + h[1]*image[i*dim+1] - xcm[1];

            // use running average position for reference if requested
            if (avflag) {
              xoriginal[i*dim+0] = (xoriginal[i*dim+0]*naverage + xtmp)*navfac;
              xoriginal[i*dim+1] = (xoriginal[i*dim+0]*naverage + xtmp)*navfac;
            }

            dx = xtmp - xoriginal[i*dim+0];
            dy = ytmp - xoriginal[i*dim+1];
            msd[0] += dx*dx;
            msd[1] += dy*dy;
            msd[3] += dx*dx + dy*dy;
          }
      }      
  }
  
  //MPI_Allreduce(msd,vector,4,MPI_DOUBLE,MPI_SUM,world);  
    vector[0] /= nmsd;
    vector[1] /= nmsd;
    vector[2] /= nmsd;
    vector[3] /= nmsd;  
}
template void cpuComputeMSD(double *vector, double *x, double *xoriginal, double *h, double *xcm,
         int *ilist, int *image, int naverage, int avflag, int triclinic, int nmsd, int dim,  int inum);
template void cpuComputeMSD(float *vector, float *x, float *xoriginal, float *h, float *xcm,
         int *ilist, int *image, int naverage, int avflag, int triclinic, int nmsd, int dim,  int inum);

template <typename T> void cpuComputeVACF(T *vacf, T *v, T *voriginal, 
         int *ilist, int nvacf, int dim,  int inum)
{
  T vxsq,vysq,vzsq;
  vacf[0] = vacf[1] = vacf[2] = vacf[3] = 0.0;

  for (int ii = 0; ii < inum; ii++) {
      int i = ilist[ii];
      vxsq = v[i*dim+0] * voriginal[i*dim+0];
      vysq = v[i*dim+1] * voriginal[i*dim+1];
      vzsq = (dim==3) ? (v[i*dim+2] * voriginal[i*dim+2]) : 0.0;
      vacf[0] += vxsq;
      vacf[1] += vysq;
      vacf[2] += vzsq;
      vacf[3] += vxsq + vysq + vzsq;
  }

  //MPI_Allreduce(vacf,vector,4,MPI_DOUBLE,MPI_SUM,world);
    vacf[0] /= nvacf;
    vacf[1] /= nvacf;
    vacf[2] /= nvacf;
    vacf[3] /= nvacf;
}
template void cpuComputeVACF(double *vacf, double *v, double *voriginal, 
         int *ilist, int nvacf, int dim,  int inum);
template void cpuComputeVACF(float *vacf, float *v, float *voriginal, 
         int *ilist, int nvacf, int dim,  int inum);

#endif

