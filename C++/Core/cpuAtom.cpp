template <typename T> void cpuAtomInside(int *inside, T *x, T *h_inv, T *boxlo, T *lo_lamda,
        T *hi_lamda, int dim, int n)
{
    for (int i = 0; i < n; i++) {
        T lamda[3];
        int k = i*dim;
        T deltax = x[k+0] - boxlo[0];
        T deltay = x[k+1] - boxlo[1];
        lamda[0] = h_inv[0]*deltax + h_inv[5]*deltay;
        lamda[1] = h_inv[1]*deltay;    
        int insidesub = 0;
        if (dim==3) {
            T deltaz = x[k+2] - boxlo[2];
            lamda[0] += h_inv[4]*deltaz;
            lamda[1] += h_inv[3]*deltaz;
            lamda[2] = h_inv[2]*deltaz;
            
            if (lamda[0] >= lo_lamda[0] && lamda[0] < hi_lamda[0] &&
                lamda[1] >= lo_lamda[1] && lamda[1] < hi_lamda[1] &&
                lamda[2] >= lo_lamda[2] && lamda[2] < hi_lamda[2]) insidesub = 1;            
        }
        else {
            if (lamda[0] >= lo_lamda[0] && lamda[0] < hi_lamda[0] &&
                lamda[1] >= lo_lamda[1] && lamda[1] < hi_lamda[1]) insidesub = 1;                        
        }
            
        inside[i] = (insidesub) ? 1 : 0;
    }
}
template void cpuAtomInside(int *inside, double *x, double *h_inv, double *boxlo, double *lo_lamda,
        double *hi_lamda, int dim, int n);
template void cpuAtomInside(int *inside, float *x, float *h_inv, float *boxlo, float *lo_lamda,
        float *hi_lamda, int dim, int n);

template <typename T> int cpuAtomAdd(T *y, int *atomtype, T *x, T *h_inv, T *boxlo, T *lo_lamda,
        T *hi_lamda, int *type, int dim, int nlocal, int n)
{
    int ninside = 0;
    for (int i = 0; i < n; i++) {
        T lamda[dim];
        int k = i*dim;
        T deltax = x[k+0] - boxlo[0];
        T deltay = x[k+1] - boxlo[1];
        lamda[0] = h_inv[0]*deltax + h_inv[5]*deltay;
        lamda[1] = h_inv[1]*deltay;    
        int insidesub = 0;
        if (dim==3) {
            T deltaz = x[k+2] - boxlo[2];
            lamda[0] += h_inv[4]*deltaz;
            lamda[1] += h_inv[3]*deltaz;
            lamda[2] = h_inv[2]*deltaz;
            
            if (lamda[0] >= lo_lamda[0] && lamda[0] < hi_lamda[0] &&
                lamda[1] >= lo_lamda[1] && lamda[1] < hi_lamda[1] &&
                lamda[2] >= lo_lamda[2] && lamda[2] < hi_lamda[2]) insidesub = 1;            
        }
        else {
            if (lamda[0] >= lo_lamda[0] && lamda[0] < hi_lamda[0] &&
                lamda[1] >= lo_lamda[1] && lamda[1] < hi_lamda[1]) insidesub = 1;                        
        }
            
        if (insidesub) {
            atomtype[ninside+nlocal] = type[i];   
            for (int j=0; j<dim; j++) 
                y[(ninside+nlocal)*dim+j] = x[k+j];                                 
            ninside += 1;
        }
    }
    
    return ninside;
}
template int cpuAtomAdd(double *y, int *atomtype, double *x, double *h_inv, double *boxlo, double *lo_lamda,
        double *hi_lamda, int *type, int dim, int nlocal, int n);
template int cpuAtomAdd(float *y, int *atomtype, float *x, float *h_inv, float *boxlo, float *lo_lamda,
        float *hi_lamda, int *type, int dim, int nlocal, int n);

// template <typename T> int cpuAtomAdd(T *y, int *atomtype, T *x, T *h_inv, T *boxlo, T *lo_lamda,
//         T *hi_lamda, int *inside, int *type, int *tmp, int dim, int nlocal, int n)
// {    
//     cpuArrayFill(tmp, 0, n) ;
//     cpuAtomInside(&tmp[n], x, h_inv, boxlo, lo_lamda, hi_lamda, dim, n);    
//     int ninside = cpuBoolFlagged(inside, tmp, &tmp[n], &tmp[2*n], n);    
//     cpuGetArrayAtColumnIndex(&y[nlocal*dim], x, inside, dim, ninside);    
//     cpuGetArrayAtIndex(atomtype, type, inside, ninside);    
//     return ninside;
// }

template <typename T> void cpuAtomLattice(T *y, int *atomtype, T *basis, T *primitive, T *rotaterow, T *origin, 
        T *latticespacing, T scale, int *basistype, int nbasis, int nlocal, 
        int ilo, int ihi, int jlo, int jhi, int klo, int khi, int dim)
{
  int i,j,k,m;

  int natom = 0;
  for (k = klo; k <= khi; k++) {
    for (j = jlo; j <= jhi; j++) {
      for (i = ilo; i <= ihi; i++) {
        for (m = 0; m < nbasis; m++) {
          T lamda[3];
          T x[3];

          x[0] = i + basis[m*3+0];
          x[1] = j + basis[m*3+1];
          x[2] = k + basis[m*3+2];

          // convert from lattice coords to box coords
          //cpuLattice2Box(x, primitive, rotaterow, origin, latticespacing, scale, dim);          
          if (dim==3) {  
              T x1 = primitive[0*3+0]*x[0] + primitive[0*3+1]*x[1] + primitive[0*3+2]*x[2];
              T y1 = primitive[1*3+0]*x[0] + primitive[1*3+1]*x[1] + primitive[1*3+2]*x[2];
              T z1 = primitive[2*3+0]*x[0] + primitive[2*3+1]*x[1] + primitive[2*3+2]*x[2];

              x1 *= scale;
              y1 *= scale;
              z1 *= scale;

              T xnew = rotaterow[0*3+0]*x1 + rotaterow[0*3+1]*y1 + rotaterow[0*3+2]*z1;
              T ynew = rotaterow[1*3+0]*x1 + rotaterow[1*3+1]*y1 + rotaterow[1*3+2]*z1;
              T znew = rotaterow[2*3+0]*x1 + rotaterow[2*3+1]*y1 + rotaterow[2*3+2]*z1;

              x[0] = xnew + latticespacing[0]*origin[0];
              x[1] = ynew + latticespacing[1]*origin[1];
              x[2] = znew + latticespacing[2]*origin[2];
          }
          else {
              T x1 = primitive[0*3+0]*x[0] + primitive[0*3+1]*x[1];
              T y1 = primitive[1*3+0]*x[0] + primitive[1*3+1]*x[1];      

              x1 *= scale;
              y1 *= scale;      

              T xnew = rotaterow[0*3+0]*x1 + rotaterow[0*3+1]*y1;
              T ynew = rotaterow[1*3+0]*x1 + rotaterow[1*3+1]*y1;      

              x[0] = xnew + latticespacing[0]*origin[0];
              x[1] = ynew + latticespacing[1]*origin[1];      
          }
          
          atomtype[nlocal+natom] = basistype[m];                  
          for (int d=0; d<dim; d++) 
              y[(nlocal+natom)*dim+d] = x[d];      

          natom++;         
        }
      }
    }
  }  
}
template void cpuAtomLattice(double *y, int *atomtype, double *basis, double *primitive, double *rotaterow, 
        double *origin, double *latticespacing, double scale, int *basistype, int nbasis, int nlocal, 
        int ilo, int ihi, int jlo, int jhi, int klo, int khi, int dim);
template void cpuAtomLattice(float *y, int *atomtype, float *basis, float *primitive, float *rotaterow, 
        float *origin, float *latticespacing, float scale, int *basistype, int nbasis, int nlocal, 
        int ilo, int ihi, int jlo, int jhi, int klo, int khi, int dim);


// template <typename T> int cpuAtomAddLattice(T *y, T *basis, T *boxlo, T *h_inv, T *lo_lamda, T *hi_lamda, T *primitive, 
//         T *rotaterow, T *origin, T *latticespacing, T scale, int *atomtype, int *basistype, int nbasis, 
//         int nlocal, int ilo, int ihi, int jlo, int jhi, int klo, int khi, int dim)
// {
//   int i,j,k,m;
// 
//   int nlatt = 0;
//   for (k = klo; k <= khi; k++) {
//     for (j = jlo; j <= jhi; j++) {
//       for (i = ilo; i <= ihi; i++) {
//         for (m = 0; m < nbasis; m++) {
//           T lamda[3];
//           T x[3];
// 
//           x[0] = i + basis[m*3+0];
//           x[1] = j + basis[m*3+1];
//           x[2] = k + basis[m*3+2];
// 
//           // convert from lattice coords to box coords
//           cpuLattice2Box(x, primitive, rotaterow, origin, latticespacing, scale, dim);
//           
//             T deltax = x[k+0] - boxlo[0];
//             T deltay = x[k+1] - boxlo[1];
//             lamda[0] = h_inv[0]*deltax + h_inv[5]*deltay;
//             lamda[1] = h_inv[1]*deltay;    
//             int insidesub = 0;
//             if (dim==3) {
//                 T deltaz = x[k+2] - boxlo[2];
//                 lamda[0] += h_inv[4]*deltaz;
//                 lamda[1] += h_inv[3]*deltaz;
//                 lamda[2] = h_inv[2]*deltaz;
// 
//                 if (lamda[0] >= lo_lamda[0] && lamda[0] < hi_lamda[0] &&
//                     lamda[1] >= lo_lamda[1] && lamda[1] < hi_lamda[1] &&
//                     lamda[2] >= lo_lamda[2] && lamda[2] < hi_lamda[2]) insidesub = 1;            
//             }
//             else {
//                 if (lamda[0] >= lo_lamda[0] && lamda[0] < hi_lamda[0] &&
//                     lamda[1] >= lo_lamda[1] && lamda[1] < hi_lamda[1]) insidesub = 1;                        
//             }
//                                   
//           if (insidesub) {
//               atomtype[nlocal+nlatt] = basistype[m];                  
//               for (int d=0; d<dim; d++) 
//                 y[(nlocal+nlatt)*dim+d] = x[d];      
// 
//               nlatt++;
//           }
//         }
//       }
//     }
//   }
//   
//   return nlatt;
// }

// template <typename T> int cpuAtomAddLattice(T *y, T *basis, T *boxlo, T *h_inv, T *slo, T *shi, T *primitive, 
//         T *rotaterow, T *origin, T *latticespacing, T scale, int *atomtype, int *basistype, int nbasis, 
//         int nlocal, int ilo, int ihi, int jlo, int jhi, int klo, int khi, int dim, int triclinic)
// {
//   int i,j,k,m;
// 
//   int nlatt = 0;
//   for (k = klo; k <= khi; k++) {
//     for (j = jlo; j <= jhi; j++) {
//       for (i = ilo; i <= ihi; i++) {
//         for (m = 0; m < nbasis; m++) {
//           T coord[3];
//           T x[3];
// 
//           x[0] = i + basis[m*3+0];
//           x[1] = j + basis[m*3+1];
//           x[2] = k + basis[m*3+2];
// 
//           // convert from lattice coords to box coords
//           cpuLattice2Box(x, primitive, rotaterow, origin, latticespacing, scale, dim);
//           
//           // test if atom/molecule position is in my subbox
//           if (triclinic) {            
//             cpuBox2Lamda(coord, x, h_inv, boxlo, dim, 1);
//           } else {
//               for (int d=0; d<dim; d++)
//                 coord[d] = x[d];
//           }
//           
//           if (coord[0] < slo[0] || coord[0] >= shi[0] ||
//               coord[1] < slo[1] || coord[1] >= shi[1] ||
//               coord[2] < slo[2] || coord[2] >= shi[2]) continue;
// 
//           atomtype[nlocal+nlatt] = basistype[m];                  
//           for (int d=0; d<dim; d++) 
//             y[(nlocal+nlatt)*dim+d] = x[d];      
//           
//           nlatt++;
//         }
//       }
//     }
//   }
//   
//   return nlatt;
// }
// 
// 
