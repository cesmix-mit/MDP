/***************************************************************************                               
                    Molecular Dynamics Potentials (MDP)
                           CESMIX-MIT Project  
 
 Contributing authors: Ngoc-Cuong Nguyen (cuongng@mit.edu, exapde@gmail.com)
 ***************************************************************************/
//#define BIG 1.0e30

/* ----------------------------------------------------------------------
   check if 3 orientation vectors are mutually orthogonal
------------------------------------------------------------------------- */

int cpuLatticeOrthogonal(int *orientx, int *orienty, int *orientz)
{
  if (orientx[0]*orienty[0] + orientx[1]*orienty[1] +
      orientx[2]*orienty[2]) return 0;
  if (orienty[0]*orientz[0] + orienty[1]*orientz[1] +
      orienty[2]*orientz[2]) return 0;
  if (orientx[0]*orientz[0] + orientx[1]*orientz[1] +
      orientx[2]*orientz[2]) return 0;
  return 1;
}

/* ----------------------------------------------------------------------
   check righthandedness of orientation vectors
   x cross y must be in same direction as z
------------------------------------------------------------------------- */

int cpuLatticeRightHanded(int *orientx, int *orienty, int *orientz)
{
  int xy0 = orientx[1]*orienty[2] - orientx[2]*orienty[1];
  int xy1 = orientx[2]*orienty[0] - orientx[0]*orienty[2];
  int xy2 = orientx[0]*orienty[1] - orientx[1]*orienty[0];
  if (xy0*orientz[0] + xy1*orientz[1] + xy2*orientz[2] <= 0) return 0;
  return 1;
}

/* ----------------------------------------------------------------------
   check collinearity of each pair of primitive vectors
------------------------------------------------------------------------- */
template <typename T> int cpuLatticeCollinear(T *a1, T *a2, T *a3)
{
  T vec[3];
  cross(a1,a2,vec);
  if (dot(vec,vec) == 0.0) return 1;
  cross(a2,a3,vec);
  if (dot(vec,vec) == 0.0) return 1;
  cross(a1,a3,vec);
  if (dot(vec,vec) == 0.0) return 1;
  return 0;
}

template <typename T> void cpuLatticeTransform(T *primitive, T *priminv, T *rotaterow, T *rotatecol, 
        T *a1, T *a2, T *a3, int *orientx, int *orienty, int *orientz)
{
  T length;

  // primitive = 3x3 matrix with primitive vectors as columns

  primitive[0*3+0] = a1[0];
  primitive[1*3+0] = a1[1];
  primitive[2*3+0] = a1[2];
  primitive[0*3+1] = a2[0];
  primitive[1*3+1] = a2[1];
  primitive[2*3+1] = a2[2];
  primitive[0*3+2] = a3[0];
  primitive[1*3+2] = a3[1];
  primitive[2*3+2] = a3[2];

  // priminv = inverse of primitive

  T determinant = primitive[0*3+0]*primitive[1*3+1]*primitive[2*3+2] +
    primitive[0*3+1]*primitive[1*3+2]*primitive[2*3+0] +
    primitive[0*3+2]*primitive[1*3+0]*primitive[2*3+1] -
    primitive[0*3+0]*primitive[1*3+2]*primitive[2*3+1] -
    primitive[0*3+1]*primitive[1*3+0]*primitive[2*3+2] -
    primitive[0*3+2]*primitive[1*3+1]*primitive[2*3+0];

  priminv[0*3+0] = (primitive[1*3+1]*primitive[2*3+2] -
                   primitive[1*3+2]*primitive[2*3+1]) / determinant;
  priminv[1*3+0] = (primitive[1*3+2]*primitive[2*3+0] -
                   primitive[1*3+0]*primitive[2*3+2]) / determinant;
  priminv[2*3+0] = (primitive[1*3+0]*primitive[2*3+1] -
                   primitive[1*3+1]*primitive[2*3+0]) / determinant;

  priminv[0*3+1] = (primitive[0*3+2]*primitive[2*3+1] -
                   primitive[0*3+1]*primitive[2*3+2]) / determinant;
  priminv[1*3+1] = (primitive[0*3+0]*primitive[2*3+2] -
                   primitive[0*3+2]*primitive[2*3+0]) / determinant;
  priminv[2*3+1] = (primitive[0*3+1]*primitive[2*3+0] -
                   primitive[0*3+0]*primitive[2*3+1]) / determinant;

  priminv[0*3+2] = (primitive[0*3+1]*primitive[1*3+2] -
                   primitive[0*3+2]*primitive[1*3+1]) / determinant;
  priminv[1*3+2] = (primitive[0*3+2]*primitive[1*3+0] -
                   primitive[0*3+0]*primitive[1*3+2]) / determinant;
  priminv[2*3+2] = (primitive[0*3+0]*primitive[1*3+1] -
                   primitive[0*3+1]*primitive[1*3+0]) / determinant;

  // rotaterow = 3x3 matrix with normalized orient vectors as rows

  int lensq = orientx[0]*orientx[0] + orientx[1]*orientx[1] +
    orientx[2]*orientx[2];
  length = sqrt((T) lensq);

  rotaterow[0*3+0] = orientx[0] / length;
  rotaterow[0*3+1] = orientx[1] / length;
  rotaterow[0*3+2] = orientx[2] / length;

  lensq = orienty[0]*orienty[0] + orienty[1]*orienty[1] +
    orienty[2]*orienty[2];
  length = sqrt((T) lensq);

  rotaterow[1*3+0] = orienty[0] / length;
  rotaterow[1*3+1] = orienty[1] / length;
  rotaterow[1*3+2] = orienty[2] / length;

  lensq = orientz[0]*orientz[0] + orientz[1]*orientz[1] +
    orientz[2]*orientz[2];
  length = sqrt((T) lensq);

  rotaterow[2*3+0] = orientz[0] / length;
  rotaterow[2*3+1] = orientz[1] / length;
  rotaterow[2*3+2] = orientz[2] / length;

  // rotatecol = 3x3 matrix with normalized orient vectors as columns

  rotatecol[0*3+0] = rotaterow[0*3+0];
  rotatecol[1*3+0] = rotaterow[0*3+1];
  rotatecol[2*3+0] = rotaterow[0*3+2];

  rotatecol[0*3+1] = rotaterow[1*3+0];
  rotatecol[1*3+1] = rotaterow[1*3+1];
  rotatecol[2*3+1] = rotaterow[1*3+2];

  rotatecol[0*3+2] = rotaterow[2*3+0];
  rotatecol[1*3+2] = rotaterow[2*3+1];
  rotatecol[2*3+2] = rotaterow[2*3+2];
}

/* ----------------------------------------------------------------------
   convert lattice coords to box coords
   input x,y,z = point in lattice coords
   output x,y,z = point in box coords
   transformation: xyz_box = Rotate_row * scale * P * xyz_lattice + offset
     xyz_box = 3-vector of output box coords
     Rotate_row = 3x3 matrix = normalized orient vectors as rows
     scale = scale factor
     P = 3x3 matrix = primitive vectors as columns
     xyz_lattice = 3-vector of input lattice coords
     offset = 3-vector = (xlatt*origin[0], ylatt*origin[1], zlatt*origin[2])
------------------------------------------------------------------------- */

template <typename T> void cpuLattice2Box(T *x, T *primitive, T *rotaterow, T *origin, 
        T *spacing, T scale, int dim)
{
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

      x[0] = xnew + spacing[0]*origin[0];
      x[1] = ynew + spacing[1]*origin[1];
      x[2] = znew + spacing[2]*origin[2];
  }
  else {
      T x1 = primitive[0*3+0]*x[0] + primitive[0*3+1]*x[1];
      T y1 = primitive[1*3+0]*x[0] + primitive[1*3+1]*x[1];      

      x1 *= scale;
      y1 *= scale;      

      T xnew = rotaterow[0*3+0]*x1 + rotaterow[0*3+1]*y1;
      T ynew = rotaterow[1*3+0]*x1 + rotaterow[1*3+1]*y1;      

      x[0] = xnew + spacing[0]*origin[0];
      x[1] = ynew + spacing[1]*origin[1];      
  }
}

/* ----------------------------------------------------------------------
   convert box coords to lattice coords
   input x,y,z = point in box coords
   output x,y,z = point in lattice coords
   transformation: xyz_latt = P_inv * 1/scale * Rotate_col * (xyz_box - offset)
     xyz_lattice = 3-vector of output lattice coords
     P_inv = 3x3 matrix = inverse of primitive vectors as columns
     scale = scale factor
     Rotate_col = 3x3 matrix = normalized orient vectors as columns
     xyz_box = 3-vector of input box coords
     offset = 3-vector = (xlatt*origin[0], ylatt*origin[1], zlatt*origin[2])
------------------------------------------------------------------------- */

template <typename T> void cpuBox2Lattice(T *x, T *priminv, T *rotatecol, T *origin, 
        T *spacing, T scale, int dim)
{
  if (dim==3) {
      x[0] -= spacing[0]*origin[0];
      x[1] -= spacing[1]*origin[1];
      x[2] -= spacing[2]*origin[2];

      T x1 = rotatecol[0*3+0]*x[0] + rotatecol[0*3+1]*x[1] + rotatecol[0*3+2]*x[2];
      T y1 = rotatecol[1*3+0]*x[0] + rotatecol[1*3+1]*x[1] + rotatecol[1*3+2]*x[2];
      T z1 = rotatecol[2*3+0]*x[0] + rotatecol[2*3+1]*x[1] + rotatecol[2*3+2]*x[2];

      x1 /= scale;
      y1 /= scale;
      z1 /= scale;

      x[0] = priminv[0*3+0]*x1 + priminv[0*3+1]*y1 + priminv[0*3+2]*z1;
      x[1] = priminv[1*3+0]*x1 + priminv[1*3+1]*y1 + priminv[1*3+2]*z1;
      x[2] = priminv[2*3+0]*x1 + priminv[2*3+1]*y1 + priminv[2*3+2]*z1;
  } else {
      x[0] -= spacing[0]*origin[0];
      x[1] -= spacing[1]*origin[1];

      T x1 = rotatecol[0*3+0]*x[0] + rotatecol[0*3+1]*x[1];
      T y1 = rotatecol[1*3+0]*x[0] + rotatecol[1*3+1]*x[1];

      x1 /= scale;
      y1 /= scale;

      x[0] = priminv[0*3+0]*x1 + priminv[0*3+1]*y1;
      x[1] = priminv[1*3+0]*x1 + priminv[1*3+1]*y1;
  }
}

/* ----------------------------------------------------------------------
   convert x,y,z from lattice coords to box coords (flag = 0) or vice versa
   use new point to expand bounding box (min to max)
------------------------------------------------------------------------- */
template <typename T> void cpuLatticeBbox(T *pmin, T *pmax, T *x, T *primitive, T *rotaterow, T *priminv, T *rotatecol,
        T *origin, T *spacing, T scale, int flag, int dim)
{    
  if (flag == 0) cpuLattice2Box(x, primitive, rotaterow, origin, spacing, scale, dim);
  else cpuBox2Lattice(x, priminv, rotatecol, origin, spacing, scale, dim);
        
  for (int i =0; i<dim; i++) {
    pmin[i] = MIN(x[i],pmin[i]); 
    pmax[i] = MAX(x[i],pmax[i]);
  }
}

template <typename T> int cpuLattice(T *basis, T *primitive, T *rotaterow, T *priminv, T *rotatecol, T *origin, 
        T *spacing, T *a1, T *a2, T *a3, T &scale, int *orientx, int *orienty, int *orientz, 
        int style, int unit_style, int spaceflag, int dimension)
{
  int nbasis = 0;  
  int NONE = 0, SC = 1, BCC = 2, FCC = 3, HCP = 4, DIAMOND = 5, SQ = 6, SQ2 = 7, HEX = 8, CUSTOM = 9;
  int LJ = 0;
  
  if (style == NONE) {    
    spacing[0] = spacing[1] = spacing[2] = 1.0;        
    return nbasis;
  }
  
    // set defaults 
//   origin[0] = origin[1] = origin[2] = 0.0;
//   orientx[0] = 1;  orientx[1] = 0;  orientx[2] = 0;
//   orienty[0] = 0;  orienty[1] = 1;  orienty[2] = 0;
//   orientz[0] = 0;  orientz[1] = 0;  orientz[2] = 1;
  
  if (style != CUSTOM) {
      a1[0] = 1.0;  a1[1] = 0.0;  a1[2] = 0.0;
      a2[0] = 0.0;  a2[1] = 1.0;  a2[2] = 0.0;
      a3[0] = 0.0;  a3[1] = 0.0;  a3[2] = 1.0;

      if (style == HEX) a2[1] = sqrt(3.0);
      if (style == HCP) {
        a2[1] = sqrt(3.0);
        a3[2] = sqrt(8.0/3.0);
      }

      if (style == SC) {        
        basis[nbasis*3+0] = 0.0; basis[nbasis*3+1] = 0.0; basis[nbasis*3+2] = 0.0; nbasis++;  
      } else if (style == BCC) {        
        basis[nbasis*3+0] = 0.0; basis[nbasis*3+1] = 0.0; basis[nbasis*3+2] = 0.0; nbasis++;        
        basis[nbasis*3+0] = 0.5; basis[nbasis*3+1] = 0.5; basis[nbasis*3+2] = 0.5; nbasis++;
      } else if (style == FCC) {
        basis[nbasis*3+0] = 0.0; basis[nbasis*3+1] = 0.0; basis[nbasis*3+2] = 0.0; nbasis++;        
        basis[nbasis*3+0] = 0.5; basis[nbasis*3+1] = 0.5; basis[nbasis*3+2] = 0.0; nbasis++;        
        basis[nbasis*3+0] = 0.5; basis[nbasis*3+1] = 0.0; basis[nbasis*3+2] = 0.5; nbasis++;        
        basis[nbasis*3+0] = 0.0; basis[nbasis*3+1] = 0.5; basis[nbasis*3+2] = 0.5; nbasis++;        
      } else if (style == HCP) {
        basis[nbasis*3+0] = 0.0; basis[nbasis*3+1] = 0.0; basis[nbasis*3+2] = 0.0; nbasis++;        
        basis[nbasis*3+0] = 0.5; basis[nbasis*3+1] = 0.5; basis[nbasis*3+2] = 0.0; nbasis++;        
        basis[nbasis*3+0] = 0.5; basis[nbasis*3+1] = 5.0/6.0; basis[nbasis*3+2] = 0.5; nbasis++;        
        basis[nbasis*3+0] = 0.0; basis[nbasis*3+1] = 1.0/3.0; basis[nbasis*3+2] = 0.5; nbasis++;                
      } else if (style == SQ) {
        basis[nbasis*3+0] = 0.0; basis[nbasis*3+1] = 0.0; basis[nbasis*3+2] = 0.0; nbasis++;        
      } else if (style == SQ2) {
        basis[nbasis*3+0] = 0.0; basis[nbasis*3+1] = 0.0; basis[nbasis*3+2] = 0.0; nbasis++;        
        basis[nbasis*3+0] = 0.5; basis[nbasis*3+1] = 0.5; basis[nbasis*3+2] = 0.0; nbasis++;        
      } else if (style == HEX) {
        basis[nbasis*3+0] = 0.0; basis[nbasis*3+1] = 0.0; basis[nbasis*3+2] = 0.0; nbasis++;        
        basis[nbasis*3+0] = 0.5; basis[nbasis*3+1] = 0.5; basis[nbasis*3+2] = 0.0; nbasis++;                  
      } else if (style == DIAMOND) {
        basis[nbasis*3+0] = 0.0; basis[nbasis*3+1] = 0.0; basis[nbasis*3+2] = 0.0; nbasis++;        
        basis[nbasis*3+0] = 0.5; basis[nbasis*3+1] = 0.5; basis[nbasis*3+2] = 0.0; nbasis++;        
        basis[nbasis*3+0] = 0.5; basis[nbasis*3+1] = 0.0; basis[nbasis*3+2] = 0.5; nbasis++;        
        basis[nbasis*3+0] = 0.0; basis[nbasis*3+1] = 0.5; basis[nbasis*3+2] = 0.5; nbasis++;                
        basis[nbasis*3+0] = 0.25; basis[nbasis*3+1] = 0.25; basis[nbasis*3+2] = 0.25; nbasis++;        
        basis[nbasis*3+0] = 0.25; basis[nbasis*3+1] = 0.75; basis[nbasis*3+2] = 0.75; nbasis++;        
        basis[nbasis*3+0] = 0.75; basis[nbasis*3+1] = 0.25; basis[nbasis*3+2] = 0.75; nbasis++;        
        basis[nbasis*3+0] = 0.75; basis[nbasis*3+1] = 0.75; basis[nbasis*3+2] = 0.25; nbasis++;                
      }
  }
  
  if (unit_style == LJ) {
    T vec[3];
    cross3(a2,a3,vec);
    T volume = dot3(a1,vec);
    scale = pow(nbasis/volume/scale,1.0/dimension);
  }

  // initialize lattice <-> box transformation matrices  
  cpuLatticeTransform(primitive, priminv, rotaterow, rotatecol, 
        a1, a2, a3, orientx, orienty, orientz);
        
  
// template <typename T> void cpuLatticeTransform(T **primitive, T **priminv, T **rotaterow, T **rotatecol, 
//         T *a1, T *a2, T *a3, int *orientx, int *orienty, int *orientz)
  
  // convert 8 corners of primitive unit cell from lattice coords to box coords
  // min to max = bounding box around the pts in box space
  // xlattice,ylattice,zlattice = extent of bbox in box space
  // set xlattice,ylattice,zlattice to 0.0 initially
  //   since bbox uses them to shift origin (irrelevant for this computation)

  if (spaceflag == 0) {
    T pmin[3], pmax[3], x[3];
    pmin[0] = pmin[1] = pmin[2] = BIG;
    pmax[0] = pmax[1] = pmax[2] = -BIG;
    spacing[0] = spacing[1] = spacing[2] = 0.0;

// template <typename T> void cpuLatticeBbox(T *pmin, T *pmax, T *x, T *primitive, T *rotaterow, T *priminv, T *rotatecol,
//         T *origin, T *spacing, T scale)
    
    //printf("%g \n", scale);
    x[0] = 0.0; x[1] = 0.0; x[2] = 0.0;
    cpuLatticeBbox(pmin, pmax, x, primitive, rotaterow, priminv, rotatecol,
        origin, spacing, scale, 0, 3);
    //printf("%g %g %g %g %g %g\n", pmin[0], pmin[1], pmin[2], pmax[0], pmax[1], pmax[2]);
    x[0] = 1.0; x[1] = 0.0; x[2] = 0.0;
    cpuLatticeBbox(pmin, pmax, x, primitive, rotaterow, priminv, rotatecol,
        origin, spacing, scale, 0, 3);
    //printf("%g %g %g %g %g %g\n", pmin[0], pmin[1], pmin[2], pmax[0], pmax[1], pmax[2]);
    x[0] = 0.0; x[1] = 1.0; x[2] = 0.0;
    cpuLatticeBbox(pmin, pmax, x, primitive, rotaterow, priminv, rotatecol,
        origin, spacing, scale, 0, 3);
    //printf("%g %g %g %g %g %g\n", pmin[0], pmin[1], pmin[2], pmax[0], pmax[1], pmax[2]);
    x[0] = 1.0; x[1] = 1.0; x[2] = 0.0;
    cpuLatticeBbox(pmin, pmax, x, primitive, rotaterow, priminv, rotatecol,
        origin, spacing, scale, 0, 3);
    x[0] = 0.0; x[1] = 0.0; x[2] = 1.0;
    cpuLatticeBbox(pmin, pmax, x, primitive, rotaterow, priminv, rotatecol,
        origin, spacing, scale, 0, 3);
    x[0] = 1.0; x[1] = 0.0; x[2] = 1.0;
    cpuLatticeBbox(pmin, pmax, x, primitive, rotaterow, priminv, rotatecol,
        origin, spacing, scale, 0, 3);
    x[0] = 0.0; x[1] = 1.0; x[2] = 1.0;
    cpuLatticeBbox(pmin, pmax, x, primitive, rotaterow, priminv, rotatecol,
        origin, spacing, scale, 0, 3);
    x[0] = 1.0; x[1] = 1.0; x[2] = 1.0;
    cpuLatticeBbox(pmin, pmax, x, primitive, rotaterow, priminv, rotatecol,
        origin, spacing, scale, 0, 3);

    for (int i=0; i<3; i++)
        spacing[i] = pmax[i] - pmin[i];

  } else {
      for (int i=0; i<3; i++)
        spacing[i] *= scale;
  }
  
  return nbasis;
}
template int cpuLattice(double *basis, double *primitive, double *rotaterow, double *priminv, double *rotatecol, 
        double *origin, double *spacing, double *a1, double *a2, double *a3, double &scale, int *orientx,
        int *orienty, int *orientz, int style, int unit_style, int spaceflag, int dimension);
template int cpuLattice(float *basis, float *primitive, float *rotaterow, float *priminv, float *rotatecol, 
        float *origin, float *spacing, float *a1, float *a2, float *a3, float &scale, int *orientx,
        int *orienty, int *orientz, int style, int unit_style, int spaceflag, int dimension);        
/* ----------------------------------------------------------------------
   return bounding box around 8 corner pts in lattice coords 
------------------------------------------------------------------------- */

template <typename T> void cpuLatticeBoundingBox(T *lmin, T *lmax, T *bsublo, T *bsubhi, T *primitive, 
        T *rotaterow, T *priminv, T *rotatecol, T *origin, T *spacing, T scale)
{
  int dim = 3;
  // bounding box around the box corners in lattice space
  for (int i=0; i<dim; i++) {
      lmin[i] = BIG;
      lmax[i] = -BIG;
  }            
  
  T x[3];  
  // convert to lattice coordinates and set bounding box
  x[0] = bsublo[0]; x[1] = bsublo[1]; x[2] = bsublo[2];
  cpuLatticeBbox(lmin, lmax, x, primitive, rotaterow, priminv, 
          rotatecol, origin, spacing, scale, 1, dim);          
  x[0] = bsubhi[0]; x[1] = bsublo[1]; x[2] = bsublo[2];
  cpuLatticeBbox(lmin, lmax, x, primitive, rotaterow, priminv, 
          rotatecol, origin, spacing, scale, 1, dim);
  x[0] = bsublo[0]; x[1] = bsubhi[1]; x[2] = bsublo[2];
  cpuLatticeBbox(lmin, lmax, x, primitive, rotaterow, priminv, 
          rotatecol, origin, spacing, scale, 1, dim);          
  x[0] = bsubhi[0]; x[1] = bsubhi[1]; x[2] = bsublo[2];
  cpuLatticeBbox(lmin, lmax, x, primitive, rotaterow, priminv, 
          rotatecol, origin, spacing, scale, 1, dim);          
  x[0] = bsublo[0]; x[1] = bsublo[1]; x[2] = bsubhi[2];
  cpuLatticeBbox(lmin, lmax, x, primitive, rotaterow, priminv, 
          rotatecol, origin, spacing, scale, 1, dim);          
  x[0] = bsubhi[0]; x[1] = bsublo[1]; x[2] = bsubhi[2];
  cpuLatticeBbox(lmin, lmax, x, primitive, rotaterow, priminv, 
          rotatecol, origin, spacing, scale, 1, dim);          
  x[0] = bsublo[0]; x[1] = bsubhi[1]; x[2] = bsubhi[2];
  cpuLatticeBbox(lmin, lmax, x, primitive, rotaterow, priminv, 
          rotatecol, origin, spacing, scale, 1, dim);          
  x[0] = bsubhi[0]; x[1] = bsubhi[1]; x[2] = bsubhi[2];
  cpuLatticeBbox(lmin, lmax, x, primitive, rotaterow, priminv, 
          rotatecol, origin, spacing, scale, 1, dim);      
}
template void cpuLatticeBoundingBox(double *lmin, double *lmax, double *bsublo, double *bsubhi, double *primitive, 
        double *rotaterow, double *priminv, double *rotatecol, double *origin, double *spacing, double scale);
template void cpuLatticeBoundingBox(float *lmin, float *lmax, float *bsublo, float *bsubhi, float *primitive, 
        float *rotaterow, float *priminv, float *rotatecol, float *origin, float *spacing, float scale);

int cpuLatticeCount(int nbasis, int ilo, int ihi, int jlo, int jhi, int klo, int khi)
{
//   ilo = static_cast<int> (lmin[0]) - 1;
//   jlo = static_cast<int> (lmin[1]) - 1;
//   klo = static_cast<int> (lmin[2]) - 1;
//   ihi = static_cast<int> (lmax[0]) + 1;
//   jhi = static_cast<int> (lmax[1]) + 1;
//   khi = static_cast<int> (lmax[2]) + 1;
// 
//   if (lmin[0] < 0.0) ilo--;
//   if (lmin[1] < 0.0) jlo--;
//   if (lmin[2] < 0.0) klo--;
    
  int i,j,k,m;
  int natom = 0;
  for (k = klo; k <= khi; k++) 
    for (j = jlo; j <= jhi; j++) 
      for (i = ilo; i <= ihi; i++) 
        for (m = 0; m < nbasis; m++) 
            natom += 1;
      
  return natom;  
}

          
      