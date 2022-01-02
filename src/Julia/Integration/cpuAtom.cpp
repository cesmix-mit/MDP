int cpuSetAtomType(double *x, double fraction, int *atomtype,
         int *seed, int *save, int seed0, int newtype, int dim, int nlocal)
{
    int count = 0;
    for (int i = 0; i < nlocal; i++) {
        cpuRandomResetSeed(seed,save,seed0,&x[i*dim]);                
        if (cpuRandomUniform(seed) > fraction) continue;
        atomtype[i] = newtype;
        count++;   
    }
    return count;
}

void cpuAtomInside(int *inside, double *x, double *h_inv, double *boxlo, double *lo_lamda,
        double *hi_lamda, int dim, int n)
{
    for (int i = 0; i < n; i++) {
        double lamda[3];
        int k = i*dim;
        double deltax = x[k+0] - boxlo[0];
        double deltay = x[k+1] - boxlo[1];
        lamda[0] = h_inv[0]*deltax + h_inv[5]*deltay;
        lamda[1] = h_inv[1]*deltay;    
        int insidesub = 0;
        if (dim==3) {
            double deltaz = x[k+2] - boxlo[2];
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

int cpuAtomAdd(double *y, int *atomtype, double *x, double *h_inv, double *boxlo, double *lo_lamda,
        double *hi_lamda, int *type, int dim, int nlocal, int n)
{
    int ninside = 0;
    for (int i = 0; i < n; i++) {
        double lamda[3];
        int k = i*dim;
        double deltax = x[k+0] - boxlo[0];
        double deltay = x[k+1] - boxlo[1];
        lamda[0] = h_inv[0]*deltax + h_inv[5]*deltay;
        lamda[1] = h_inv[1]*deltay;    
        int insidesub = 0;
        if (dim==3) {
            double deltaz = x[k+2] - boxlo[2];
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

void cpuAtomLattice(double *y, int *atomtype, double *basis, double *primitive, double *rotaterow, double *origin, 
        double *latticespacing, double scale, int *basistype, int nbasis, int nlocal, 
        int ilo, int ihi, int jlo, int jhi, int klo, int khi, int dim)
{
  int i,j,k,m;

  int natom = 0;
  for (k = klo; k <= khi; k++) {
    for (j = jlo; j <= jhi; j++) {
      for (i = ilo; i <= ihi; i++) {
        for (m = 0; m < nbasis; m++) {
          //double lamda[3];
          double x[3];

          x[0] = i + basis[m*3+0];
          x[1] = j + basis[m*3+1];
          x[2] = k + basis[m*3+2];

          // convert from lattice coords to box coords
          //cpuLattice2Box(x, primitive, rotaterow, origin, latticespacing, scale, dim);          
          if (dim==3) {  
              double x1 = primitive[0*3+0]*x[0] + primitive[0*3+1]*x[1] + primitive[0*3+2]*x[2];
              double y1 = primitive[1*3+0]*x[0] + primitive[1*3+1]*x[1] + primitive[1*3+2]*x[2];
              double z1 = primitive[2*3+0]*x[0] + primitive[2*3+1]*x[1] + primitive[2*3+2]*x[2];

              x1 *= scale;
              y1 *= scale;
              z1 *= scale;

              double xnew = rotaterow[0*3+0]*x1 + rotaterow[0*3+1]*y1 + rotaterow[0*3+2]*z1;
              double ynew = rotaterow[1*3+0]*x1 + rotaterow[1*3+1]*y1 + rotaterow[1*3+2]*z1;
              double znew = rotaterow[2*3+0]*x1 + rotaterow[2*3+1]*y1 + rotaterow[2*3+2]*z1;

              x[0] = xnew + latticespacing[0]*origin[0];
              x[1] = ynew + latticespacing[1]*origin[1];
              x[2] = znew + latticespacing[2]*origin[2];
          }
          else {
              double x1 = primitive[0*3+0]*x[0] + primitive[0*3+1]*x[1];
              double y1 = primitive[1*3+0]*x[0] + primitive[1*3+1]*x[1];      

              x1 *= scale;
              y1 *= scale;      

              double xnew = rotaterow[0*3+0]*x1 + rotaterow[0*3+1]*y1;
              double ynew = rotaterow[1*3+0]*x1 + rotaterow[1*3+1]*y1;      

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
