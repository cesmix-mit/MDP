//#define BIG 1.0e20


/* ----------------------------------------------------------------------
   test if x is inside triangle with vertices v1,v2,v3
   norm = normal to triangle, defined by right-hand rule for v1,v2,v3 ordering
   edge = edge vector of triangle
   pvec = vector from vertex to x
   xproduct = cross product of edge with pvec
   if xproduct dot norm < 0.0 for any of 3 edges, then x is outside triangle
------------------------------------------------------------------------- */

int cpuInsideTri(double *x, double *v1, double *v2, double *v3, double *norm)
{
  double edge[3],pvec[3],xproduct[3];

  sub3(v2,v1,edge);
  sub3(x,v1,pvec);
  cross3(edge,pvec,xproduct);
  if (dot3(xproduct,norm) < 0.0) return 0;

  sub3(v3,v2,edge);
  sub3(x,v2,pvec);
  cross3(edge,pvec,xproduct);
  if (dot3(xproduct,norm) < 0.0) return 0;

  sub3(v1,v3,edge);
  sub3(x,v3,pvec);
  cross3(edge,pvec,xproduct);
  if (dot3(xproduct,norm) < 0.0) return 0;

  return 1;
}

/* ---------------------------------------------------------------------- */

double cpuClosest(double *x, double *near, double *nearest, double dsq)
{
  double delx = x[0] - near[0];
  double dely = x[1] - near[1];
  double delz = x[2] - near[2];
  double rsq = delx*delx + dely*dely + delz*delz;
  if (rsq >= dsq) return dsq;

  nearest[0] = near[0];
  nearest[1] = near[1];
  nearest[2] = near[2];
  return rsq;
}

/* ----------------------------------------------------------------------
   add a single contact at Nth location in contact array
   x = particle position
   xp,yp,zp = region surface point
------------------------------------------------------------------------- */

void cpuAddContact(double *r, double *radius, double *delx, double *dely, double *delz, double *x, double *p, int n)
{
  double dx = x[0] - p[0];
  double dy = x[1] - p[1];
  double dz = x[2] - p[2];
  r[n] = sqrt(dx*dx + dy*dy + dz*dz);
  radius[n] = 0;
  delx[n] = dx;
  dely[n] = dy;
  delz[n] = dz;
}

/* ----------------------------------------------------------------------
   rotate x,y,z by angle via right-hand rule around point and runit normal
   sign of angle determines whether rotating forward/backward in time
   return updated x,y,z
   R = vector axis of rotation
   P = point = point to rotate around
   R0 = runit = unit vector for R
   X0 = x,y,z = initial coord of atom
   D = X0 - P = vector from P to X0
   C = (D dot R0) R0 = projection of D onto R, i.e. Dparallel
   A = D - C = vector from R line to X0, i.e. Dperp
   B = R0 cross A = vector perp to A in plane of rotation, same len as A
   A,B define plane of circular rotation around R line
   new x,y,z = P + C + A cos(angle) + B sin(angle)
------------------------------------------------------------------------- */

void cpuRotate(double *x, double *point, double *runit, double angle)
{
  double a[3],b[3],c[3],d[3],disp[3];

  double sine = sin(angle);
  double cosine = cos(angle);
  d[0] = x[0] - point[0];
  d[1] = x[1] - point[1];
  d[2] = x[2] - point[2];
  double x0dotr = d[0]*runit[0] + d[1]*runit[1] + d[2]*runit[2];
  c[0] = x0dotr * runit[0];
  c[1] = x0dotr * runit[1];
  c[2] = x0dotr * runit[2];
  a[0] = d[0] - c[0];
  a[1] = d[1] - c[1];
  a[2] = d[2] - c[2];
  b[0] = runit[1]*a[2] - runit[2]*a[1];
  b[1] = runit[2]*a[0] - runit[0]*a[2];
  b[2] = runit[0]*a[1] - runit[1]*a[0];
  disp[0] = a[0]*cosine  + b[0]*sine;
  disp[1] = a[1]*cosine  + b[1]*sine;
  disp[2] = a[2]*cosine  + b[2]*sine;
  x[0] = point[0] + c[0] + disp[0];
  x[1] = point[1] + c[1] + disp[1];
  x[2] = point[2] + c[2] + disp[2];
}

/* ----------------------------------------------------------------------
   transform a point x,y,z in region space to moved space
   rotate first (around original P), then displace
------------------------------------------------------------------------- */

void cpuForwardTransform(double *x, double *dx, double *point, double *runit, double theta, int rotateflag, int moveflag)
{
  //(double &x, double &y, double &z, double *point, double *runit, double angle)
  if (rotateflag) cpuRotate(x,point,runit,theta);
  if (moveflag) {
    x[0] += dx[0];
    x[1] += dx[1];
    x[2] += dx[2];
  }
}

/* ----------------------------------------------------------------------
   transform a point x,y,z in moved space back to region space
   undisplace first, then unrotate (around original P)
------------------------------------------------------------------------- */

void cpuInverseTransform(double *x, double *dx, double *point, double *runit, double theta, int rotateflag, int moveflag)
{
  if (moveflag) {
    x[0] -= dx[0];
    x[1] -= dx[1];
    x[2] -= dx[2];
  }
  if (rotateflag) cpuRotate(x,point,runit,-theta);
}


/* ----------------------------------------------------------------------
   parse optional parameters at end of region input line
------------------------------------------------------------------------- */

int cpuRegionSetup(double *point, double *runit, double *scale, double *axis, double *lattice, 
        int scaleflag, int rotateflag, int moveflag)
{  
  for (int i=0; i<3; i++)
    scale[i] = (scaleflag) ? lattice[i] : 1.0;
  
  if (rotateflag) {
    // runit = unit vector along rotation axis
    double len = sqrt(axis[0]*axis[0] + axis[1]*axis[1] + axis[2]*axis[2]);
    for (int i=0; i<3; i++) {
        point[i] *= scale[i];
        runit[i] = axis[i]/len;
    }
  }

  int dynamic;
  if (moveflag || rotateflag) dynamic = 1;
  else dynamic = 0;
    
  return dynamic;  
}

/* ----------------------------------------------------------------------
   find nearest point to C on line segment A,B and return it as D
   project (C-A) onto (B-A)
   t = length of that projection, normalized by length of (B-A)
   t <= 0, C is closest to A
   t >= 1, C is closest to B
   else closest point is between A and B
------------------------------------------------------------------------- */

void cpuPointOnLineSegment(double *a, double *b, double *c, double *d)
{
  double ba[3],ca[3];
  for (int i=0; i<3; i++) {
      ba[i] = b[i] - a[i];
      ca[i] = c[i] - a[i];
  }  
  double t = (ca[0] * ba[0] + ca[1] * ba[1] + ca[2] * ba[2])/(ba[0] * ba[0] + ba[1] * ba[1] + ba[2] * ba[2]);  
  if (t <= 0.0) {
    d[0] = a[0];
    d[1] = a[1];
    d[2] = a[2];
  } else if (t >= 1.0) {
    d[0] = b[0];
    d[1] = b[1];
    d[2] = b[2];
  } else {
    d[0] = a[0] + t*ba[0];
    d[1] = a[1] + t*ba[1];
    d[2] = a[2] + t*ba[2];
  }
}


/* ----------------------------------------------------------------------
   generate list of contact points for interior or exterior regions
   if region has variable shape, invoke shape_update() once per timestep
   if region is dynamic:
     before: inverse transform x,y,z (unmove, then unrotate)
     after: forward transform contact point xs,yx,zs (rotate, then move),
            then reset contact delx,dely,delz based on new contact point
            no need to do this if no rotation since delxyz doesn't change
   caller is responsible for wrapping this call with
     modify->clearstep_compute() and modify->addstep_compute() if needed
------------------------------------------------------------------------- */

//(double *x, double *dx, double *point, double *runit, double theta, int rotateflag, int moveflag)
// int cpuRegionSurface(double *x,  double *dx, double *delx, double *dely, double *delz, 
//         double *point, double *runit, double theta, double cutoff, int dynamic, int moveflag, int rotateflag, 
//         int openflag, int interior)
// {
//   int ncontact;  
//   double s[3],xnear[3],xorig[3];
// 
//   if (dynamic) {
//     xorig[0] = x[0];
//     xorig[1] = x[1];
//     xorig[2] = x[2];
//     //cpuInverseTransform(x,y,z);
//     cpuInverseTransform(x, dx, point, runit, theta, rotateflag, moveflag);
//   }
// 
//   xnear[0] = x[0];
//   xnear[1] = x[1];
//   xnear[2] = x[2];
// 
//   if (!openflag) {
//     if (interior) ncontact = cpuSurfaceInterior(xnear,cutoff);
//     else ncontact = cpuSurfaceExterior(xnear,cutoff);
//   } else {
//     // one of surface_int/ext() will return 0
//     // so no need to worry about offset of contact indices
//     ncontact = cpuSurfaceExterior(xnear,cutoff) + cpuSurfaceInterior(xnear,cutoff);
//   }
// 
//   if (rotateflag && ncontact) {
//     for (int i = 0; i < ncontact; i++) {
//       s[0] = xnear[0] - delx[i];
//       s[1] = xnear[1] - dely[i];
//       s[2] = xnear[2] - delz[i];
//       //cpuForwardTransform(xs,ys,zs);
//       cpuForwardTransform(s, dx, point, runit, theta, rotateflag, moveflag);
//       delx[i] = xorig[0] - s[0];
//       dely[i] = xorig[1] - s[1];
//       delz[i] = xorig[2] - s[2];
//     }
//   }
// 
//   return ncontact;
// }

/* ----------------------------------------------------------------------
   infer translational and angular velocity of region
   necessary b/c motion variables are for displacement & theta
     there is no analytic formula for v & omega
   prev[4] contains values of dx,dy,dz,theta at previous step
     used for difference, then updated to current step values
   dt is time elapsed since previous step
   rpoint = point updated by current displacement
   called by fix wall/gran/region every timestep
------------------------------------------------------------------------- */

void cpuRegionSetVelocity(double *v, double *prev, double *omega, double *dx, double *point, double *rpoint, 
        double *axis, double dt, double theta, int moveflag, int rotateflag, int ntimestep)
{
  if (moveflag) {
    if (ntimestep > 0) {
      v[0] = (dx[0] - prev[0])/dt;
      v[1] = (dx[1] - prev[1])/dt;
      v[2] = (dx[2] - prev[2])/dt;
    }
    else v[0] = v[1] = v[2] = 0.0;
    prev[0] = dx[0];
    prev[1] = dx[1];
    prev[2] = dx[2];
  }

  if (rotateflag) {
    rpoint[0] = point[0] + dx[0];
    rpoint[1] = point[1] + dx[1];
    rpoint[2] = point[2] + dx[2];
    if (ntimestep > 0) {
      double angvel = (theta-prev[3]) / dt;
      omega[0] = angvel*axis[0];
      omega[1] = angvel*axis[1];
      omega[2] = angvel*axis[2];
    }
    else omega[0] = omega[1] = omega[2] = 0.0;
    prev[3] = theta;
  }
}

/* ----------------------------------------------------------------------
   compute velocity of wall for given contact
   since contacts only store delx/y/z, need to pass particle coords
     to compute contact point
   called by fix/wall/gran/region every contact every timestep
------------------------------------------------------------------------- */

void cpuRegionVelocityContact(double *vwall, double *x, double *v, double *omega, double *delx, double *dely, double *delz, 
        double *rpoint,int moveflag, int rotateflag, int ic)
{
  double xc[3];

  vwall[0] = vwall[1] = vwall[2] = 0.0;

  if (moveflag) {
    vwall[0] = v[0];
    vwall[1] = v[1];
    vwall[2] = v[2];
  }
  if (rotateflag) {
    xc[0] = x[0] - delx[ic];
    xc[1] = x[1] - dely[ic];
    xc[2] = x[2] - delz[ic];
    vwall[0] += omega[1]*(xc[2] - rpoint[2]) - omega[2]*(xc[1] - rpoint[1]);
    vwall[1] += omega[2]*(xc[0] - rpoint[0]) - omega[0]*(xc[2] - rpoint[2]);
    vwall[2] += omega[0]*(xc[1] - rpoint[1]) - omega[1]*(xc[0] - rpoint[0]);
  }
}

/* ---------------------------------------------------------------------- */

void cpuBlockSetup(double *extent_lo, double *extent_hi, double **face, double ***corners, double *lo, double *hi, double *scale)
{
    for (int i=0; i<3; i++) {
        lo[i] = scale[i]*lo[i];
        hi[i] = scale[i]*hi[i];
    }
        
    extent_lo[0] = lo[0];
    extent_hi[0] = hi[0];
    extent_lo[1] = lo[1];
    extent_hi[1] = hi[1];
    extent_lo[2] = lo[2];
    extent_hi[2] = hi[2];
    
  // open face data structs
  face[0][0] = -1.0;
  face[0][1] = 0.0;
  face[0][2] = 0.0;
  face[1][0] = 1.0;
  face[1][1] = 0.0;
  face[1][2] = 0.0;
  face[2][0] = 0.0;
  face[2][1] = -1.0;
  face[2][2] = 0.0;
  face[3][0] = 0.0;
  face[3][1] = 1.0;
  face[3][2] = 0.0;
  face[4][0] = 0.0;
  face[4][1] = 0.0;
  face[4][2] = -1.0;
  face[5][0] = 0.0;
  face[5][1] = 0.0;
  face[5][2] = 1.0;

  // face[0]

  corners[0][0][0] = lo[0];
  corners[0][0][1] = lo[1];
  corners[0][0][2] = lo[2];
  corners[0][1][0] = lo[0];
  corners[0][1][1] = lo[1];
  corners[0][1][2] = hi[2];
  corners[0][2][0] = lo[0];
  corners[0][2][1] = hi[1];
  corners[0][2][2] = hi[2];
  corners[0][3][0] = lo[0];
  corners[0][3][1] = hi[1];
  corners[0][3][2] = lo[2];

  // face[1]

  corners[1][0][0] = hi[0];
  corners[1][0][1] = lo[1];
  corners[1][0][2] = lo[2];
  corners[1][1][0] = hi[0];
  corners[1][1][1] = lo[1];
  corners[1][1][2] = hi[2];
  corners[1][2][0] = hi[0];
  corners[1][2][1] = hi[1];
  corners[1][2][2] = hi[2];
  corners[1][3][0] = hi[0];
  corners[1][3][1] = hi[1];
  corners[1][3][2] = lo[2];

  for (int i=0; i<3; i++) {
      // face[2]  
      corners[2][0][i] = corners[0][0][i];
      corners[2][1][i] = corners[1][0][i];
      corners[2][2][i] = corners[1][1][i];
      corners[2][3][i] = corners[0][1][i];

      // face[3]
      corners[3][0][i] = corners[0][3][i];
      corners[3][1][i] = corners[0][2][i];
      corners[3][2][i] = corners[1][2][i];
      corners[3][3][i] = corners[1][3][i];

      // face[4]
      corners[4][0][i] = corners[0][0][i];
      corners[4][1][i] = corners[0][3][i];
      corners[4][2][i] = corners[1][3][i];
      corners[4][3][i] = corners[1][0][i];

      // face[5]
      corners[5][0][i] = corners[0][1][i];
      corners[5][1][i] = corners[1][1][i];
      corners[5][2][i] = corners[1][2][i];
      corners[5][3][i] = corners[0][2][i];
  }
}


/* ----------------------------------------------------------------------
   inside = 1 if x,y,z is inside or on surface
   inside = 0 if x,y,z is outside and not on surface
------------------------------------------------------------------------- */
int cpuBlockInside(double *x, double *lo, double *hi)
{
  if (x[0] >= lo[0] && x[0] <= hi[0] && x[1] >= lo[1] && x[1] <= hi[1] && x[2] >= lo[2] && x[2] <= hi[2])
    return 1;
  return 0;
}

/* ----------------------------------------------------------------------
   contact if 0 <= x < cutoff from one or more inner surfaces of block
   can be one contact for each of 6 faces
   no contact if outside (possible if called from union/intersect)
   delxyz = vector from nearest point on block to x
------------------------------------------------------------------------- */

int cpuBlockSurfaceInterior(double *r, double *delx, double *dely, double *delz, double *radius, 
        double *x, double *lo, double *hi, double cutoff, int *iwall, int *open_faces)
{
  double delta;

  // x is exterior to block

  if (x[0] < lo[0] || x[0] > hi[0] || x[1] < lo[1] || x[1] > hi[1] ||
      x[2] < lo[2] || x[2] > hi[2]) return 0;

  // x is interior to block or on its surface

  int n = 0;

  delta = x[0] - lo[0];
  if (delta < cutoff && !open_faces[0]) {
    r[n] = delta;
    delx[n] = delta;
    dely[n] = delz[n] = 0.0;
    radius[n] = 0;
    iwall[n] = 0;
    n++;
  }
  delta = hi[0] - x[0];
  if (delta < cutoff && !open_faces[1]) {
    r[n] = delta;
    delx[n] = -delta;
    dely[n] = delz[n] = 0.0;
    radius[n] = 0;
    iwall[n] = 1;
    n++;
  }

  delta = x[1] - lo[1];
  if (delta < cutoff && !open_faces[2]) {
    r[n] = delta;
    dely[n] = delta;
    delx[n] = delz[n] = 0.0;
    radius[n] = 0;
    iwall[n] = 2;
    n++;
  }
  delta = hi[1] - x[1];
  if (delta < cutoff && !open_faces[3]) {
    r[n] = delta;
    dely[n] = -delta;
    delx[n] = delz[n] = 0.0;
    radius[n] = 0;
    iwall[n] = 3;
    n++;
  }

  delta = x[2] - lo[2];
  if (delta < cutoff && !open_faces[4]) {
    r[n] = delta;
    delz[n] = delta;
    delx[n] = dely[n] = 0.0;
    radius[n] = 0;
    iwall[n] = 4;
    n++;
  }
  delta = hi[2] - x[2];
  if (delta < cutoff && !open_faces[5]) {
    r[n] = delta;
    delz[n] = -delta;
    delx[n] = dely[n] = 0.0;
    radius[n] = 0;
    iwall[n] = 5;
    n++;
  }

  return n;
}

/* ----------------------------------------------------------------------
   one contact if 0 <= x < cutoff from outer surface of block
   no contact if inside (possible if called from union/intersect)
   delxyz = vector from nearest point on block to x
------------------------------------------------------------------------- */

// int cpuBlockSurfaceExterior(double *r, double *delx, double *dely, double *delz, double *radius, 
//         double *x, double *lo, double *hi, double cutoff, int *iwall, int *open_faces, int openflag)
// {
//   double p[3];
//   double c[3],dist,mindist;
// 
//   // x is far enough from block that there is no contact
//   // x is interior to block
// 
//   if (x[0] <= lo[0]-cutoff || x[0] >= hi[0]+cutoff ||
//       x[1] <= lo[1]-cutoff || x[1] >= hi[1]+cutoff ||
//       x[2] <= lo[2]-cutoff || x[2] >= hi[2]+cutoff) return 0;
//   if (x[0] > lo[0] && x[0] < hi[0] && x[1] > lo[1] && x[1] < hi[1] &&
//       x[2] > lo[2] && x[2] < hi[2]) return 0;
// 
//   // x is exterior to block or on its surface
//   // p[0],p[1],p[2] = point on surface of block that x is closest to
//   //            could be edge or corner pt of block
//   // do not add contact point if r >= cutoff
// 
//   if (!openflag) {
//     if (x[0] < lo[0]) p[0] = lo[0];
//     else if (x[0] > hi[0]) p[0] = hi[0];
//     else p[0] = x[0];
//     if (x[1] < lo[1]) p[1] = lo[1];
//     else if (x[1] > hi[1]) p[1] = hi[1];
//     else p[1] = x[1];
//     if (x[2] < lo[2]) p[2] = lo[2];
//     else if (x[2] > hi[2]) p[2] = hi[2];
//     else p[2] = x[2];
//   } else {
//     mindist = BIG;
//     for (int i = 0; i < 6; i++) {
//       if (open_faces[i]) continue;
//       dist = cpuBlockFindClosestPoint(i,x,c);
//       if (dist < mindist) {
//         p[0] = c[0];
//         p[1] = c[1];
//         p[2] = c[2];
//         mindist = dist;
//       }
//     }
//   }
// 
//   //add_contact(0,x,p[0],p[1],p[2]);
//   cpuAddContact(r, radius, delx, dely, delz, x, p, 0);
//   iwall[0] = 0;
//   if (r[0] < cutoff) return 1;
//   return 0;
// }

/*------------------------------------------------------------------------
  return distance to closest point on surface I of block region
  store closest point in xc,yc,zc
--------------------------------------------------------------------------*/

// double cpuBlockFindClosestPoint(double *x, double **face, double ***corners, double *lo, double *hi, double *c, int i)
// {
//   double dot,d2,d2min;
//   double xr[3],xproj[3],p[3];
// 
//   xr[0] = x[0] - corners[i][0][0];
//   xr[1] = x[1] - corners[i][0][1];
//   xr[2] = x[2] - corners[i][0][2];
//   dot = face[i][0]*xr[0] + face[i][1]*xr[1] + face[i][2]*xr[2];
//   xproj[0] = xr[0] - dot*face[i][0];
//   xproj[1] = xr[1] - dot*face[i][1];
//   xproj[2] = xr[2] - dot*face[i][2];
// 
//   d2min = BIG;
// 
//   // check if point projects inside of face
// 
//   if (cpuBlockInsideFace(xproj,lo, hi, i)) {
//     d2 = d2min = dot*dot;
//     c[0] = xproj[0] + corners[i][0][0];
//     c[1] = xproj[1] + corners[i][0][1];
//     c[2] = xproj[2] + corners[i][0][2];
// 
//  // check each edge
//   } else {
//     cpuPoinOnLineSegment(corners[i][0],corners[i][1],x,p);
//     d2 = (p[0]-x[0])*(p[0]-x[0]) + (p[1]-x[1])*(p[1]-x[1]) +
//       (p[2]-x[2])*(p[2]-x[2]);
//     if (d2 < d2min) {
//       d2min = d2;
//       c[0] = p[0];
//       c[1] = p[1];
//       c[2] = p[2];
//     }
// 
//     cpuPoinOnLineSegment(corners[i][1],corners[i][2],x,p);
//     d2 = (p[0]-x[0])*(p[0]-x[0]) + (p[1]-x[1])*(p[1]-x[1]) +
//       (p[2]-x[2])*(p[2]-x[2]);
//     if (d2 < d2min) {
//       d2min = d2;
//       c[0] = p[0];
//       c[1] = p[1];
//       c[2] = p[2];      
//     }
// 
//     cpuPoinOnLineSegment(corners[i][2],corners[i][3],x,p);
//     d2 = (p[0]-x[0])*(p[0]-x[0]) + (p[1]-x[1])*(p[1]-x[1]) +
//       (p[2]-x[2])*(p[2]-x[2]);
//     if (d2 < d2min) {
//       d2min = d2;
//       c[0] = p[0];
//       c[1] = p[1];
//       c[2] = p[2];      
//     }
// 
//     cpuPoinOnLineSegment(corners[i][3],corners[i][0],x,p);
//     d2 = (p[0]-x[0])*(p[0]-x[0]) + (p[1]-x[1])*(p[1]-x[1]) +
//       (p[2]-x[2])*(p[2]-x[2]);
//     if (d2 < d2min) {
//       d2min = d2;
//       c[0] = p[0];
//       c[1] = p[1];
//       c[2] = p[2];      
//     }
//   }
// 
//   return d2min;
// }

/*------------------------------------------------------------------------
  determine if projected point is inside given face of the block
--------------------------------------------------------------------------*/

int cpuBlockInsideFace(double *xproj, double *lo, double *hi, int iface)
{
  if (iface < 2) {
    if (xproj[1] > 0 && (xproj[1] < hi[1]-lo[1]) &&
        xproj[2] > 0 && (xproj[2] < hi[2]-lo[2])) return 1;
  } else if (iface < 4) {
    if (xproj[0] > 0 && (xproj[0] < (hi[0]-lo[0])) &&
        xproj[2] > 0 && (xproj[2] < (hi[2]-lo[2]))) return 1;
  } else {
    if (xproj[0] > 0 && xproj[0] < (hi[0]-lo[0]) &&
        xproj[1] > 0 && xproj[1] < (hi[1]-lo[1])) return 1;
  }

  return 0;
}


// int cpuPrismSetup(double **h, double **hinv, double **face, double *corners, double *a, double *b, double*c, 
//         double *clo, double *chi, double *extent_lo, double *extent_hi, double *lo, double *hi, double *tilt, double *scale,
//         int **tri, int *open_faces, int interior, int openflag)
// {
//     for (int i=0; i<3; i++) {
//         lo[i] = scale[i]*lo[i];
//         hi[i] = scale[i]*hi[i];
//     }
//     
//   tilt[0] = scale[0]*tilt[0];
//   tilt[1] = scale[0]*tilt[1];
//   tilt[2] = scale[1]*tilt[2];
//     
//   // extent of prism
//   int bboxflag;
//   if (interior) {
//     bboxflag = 1;
//     extent_lo[0] = MIN(lo[0],lo[0]+tilt[0]);
//     extent_lo[0] = MIN(extent_lo[0],extent_lo[0]+tilt[1]);
//     extent_lo[1] = MIN(lo[1],lo[1]+tilt[2]);
//     extent_lo[2] = lo[2];
// 
//     extent_hi[0] = MAX(hi[0],hi[0]+tilt[0]);
//     extent_hi[0] = MAX(extent_hi[0],extent_hi[0]+tilt[1]);
//     extent_hi[1] = MAX(hi[1],hi[1]+tilt[2]);
//     extent_hi[2] = hi[2];
//   } else bboxflag = 0;
// 
//   // h = transformation matrix from tilt coords (0-1) to box coords (xyz)
//   // columns of h are edge vectors of tilted box
//   // hinv = transformation matrix from box coords to tilt coords
//   // both h and hinv are upper triangular
//   //   since 1st edge of prism is along x-axis
//   //   and bottom face of prism is in xy plane
// 
//   h[0][0] = hi[0] - lo[0];
//   h[0][1] = tilt[0];
//   h[0][2] = tilt[1];
//   h[1][1] = hi[1] - lo[1];
//   h[1][2] = tilt[2];
//   h[2][2] = hi[2] - lo[2];
// 
//   hinv[0][0] = 1.0/h[0][0];
//   hinv[0][1] = -h[0][1] / (h[0][0]*h[1][1]);
//   hinv[0][2] = (h[0][1]*h[1][2] - h[0][2]*h[1][1]) / (h[0][0]*h[1][1]*h[2][2]);
//   hinv[1][1] = 1.0/h[1][1];
//   hinv[1][2] = -h[1][2] / (h[1][1]*h[2][2]);
//   hinv[2][2] = 1.0/h[2][2];
// 
//   // corners = 8 corner points of prism
//   // order = x varies fastest, then y, finally z
//   // clo/chi = lo and hi corner pts of prism
// 
//   a[0] = hi[0]-lo[0];
//   a[1] = 0.0;
//   a[2] = 0.0;
//   b[0] = tilt[0];
//   b[1] = hi[1]-lo[1];
//   b[2] = 0.0;
//   c[0] = tilt[1];
//   c[1] = tilt[2];
//   c[2] = hi[2]-lo[2];
// 
//   clo[0] = corners[0][0] = lo[0];
//   clo[1] = corners[0][1] = lo[1];
//   clo[2] = corners[0][2] = lo[2];
// 
//   corners[1][0] = lo[0] + a[0];
//   corners[1][1] = lo[1] + a[1];
//   corners[1][2] = lo[2] + a[2];
// 
//   corners[2][0] = lo[0] + b[0];
//   corners[2][1] = lo[1] + b[1];
//   corners[2][2] = lo[2] + b[2];
// 
//   corners[3][0] = lo[0] + a[0] + b[0];
//   corners[3][1] = lo[1] + a[1] + b[1];
//   corners[3][2] = lo[2] + a[2] + b[2];
// 
//   corners[4][0] = lo[0] + c[0];
//   corners[4][1] = lo[1] + c[1];
//   corners[4][2] = lo[2] + c[2];
// 
//   corners[5][0] = lo[0] + a[0] + c[0];
//   corners[5][1] = lo[1] + a[1] + c[1];
//   corners[5][2] = lo[2] + a[2] + c[2];
// 
//   corners[6][0] = lo[0] + b[0] + c[0];
//   corners[6][1] = lo[1] + b[1] + c[1];
//   corners[6][2] = lo[2] + b[2] + c[2];
// 
//   chi[0] = corners[7][0] = lo[0] + a[0] + b[0] + c[0];
//   chi[1] = corners[7][1] = lo[1] + a[1] + b[1] + c[1];
//   chi[2] = corners[7][2] = lo[2] + a[2] + b[2] + c[2];
// 
//   // face = 6 inward-facing unit normals to prism faces
//   // order = xy plane, xz plane, yz plane
// 
//   cross3(a,b,face[0]);
//   cross3(b,a,face[1]);
//   cross3(c,a,face[2]);
//   cross3(a,c,face[3]);
//   cross3(b,c,face[4]);
//   cross3(c,b,face[5]);
// 
//   // remap open face indices to be consistent
// 
//   if (openflag) {
//     int temp[6];
//     for (int i = 0; i < 6; i++)
//       temp[i] = open_faces[i];
//     open_faces[0] = temp[4];
//     open_faces[1] = temp[5];
//     open_faces[2] = temp[2];
//     open_faces[3] = temp[3];
//     open_faces[4] = temp[0];
//     open_faces[5] = temp[1];
//   }
// 
//   for (int i = 0; i < 6; i++) norm3(face[i]);
// 
//   // tri = 3 vertices (0-7) in each of 12 triangles on 6 faces
//   // verts in each tri are ordered so that right-hand rule gives inward norm
//   // order = xy plane, xz plane, yz plane
// 
//   tri[0][0] =  0; tri[0][1] =  1; tri[0][2] =  3;
//   tri[1][0] =  0; tri[1][1] =  3; tri[1][2] =  2;
//   tri[2][0] =  4; tri[2][1] =  7; tri[2][2] =  5;
//   tri[3][0] =  4; tri[3][1] =  6; tri[3][2] =  7;
// 
//   tri[4][0] =  0; tri[4][1] =  4; tri[4][2] =  5;
//   tri[5][0] =  0; tri[5][1] =  5; tri[5][2] =  1;
//   tri[6][0] =  2; tri[6][1] =  7; tri[6][2] =  6;
//   tri[7][0] =  2; tri[7][1] =  3; tri[7][2] =  7;
// 
//   tri[8][0] =  2; tri[8][1] =  6; tri[8][2] =  4;
//   tri[9][0] =  2; tri[9][1] =  4; tri[9][2] =  0;
//   tri[10][0] = 1; tri[10][1] = 5; tri[10][2] = 7;
//   tri[11][0] = 1; tri[11][1] = 7; tri[11][2] = 3;
//   
//   return bboxflag;
// }


/* ----------------------------------------------------------------------
   inside = 1 if x,y,z is inside or on surface
   inside = 0 if x,y,z is outside and not on surface
   abc = Hinv * (xyz - xyz/lo)
   abc = tilt coords (0-1)
   Hinv = transformation matrix from box coords to tilt coords
   xyz = box coords
   xyz/lo = lower-left corner of prism
------------------------------------------------------------------------- */

// int cpuPrismInside(double *x, double **hinv, double *lo)
// {
//   double a = hinv[0][0]*(x[0]-lo[0]) + hinv[0][1]*(x[1]-lo[1]) + hinv[0][2]*(x[2]-lo[2]);
//   double b = hinv[1][1]*(x[1]-lo[1]) + hinv[1][2]*(x[2]-lo[2]);
//   double c = hinv[2][2]*(x[2]-lo[2]);
// 
//   if (a >= 0.0 && a <= 1.0 && b >= 0.0 && b <= 1.0 && c >= 0.0 && c <= 1.0)
//     return 1;
//   return 0;
// }
// 
// /* ----------------------------------------------------------------------
//    contact if 0 <= x < cutoff from one or more inner surfaces of prism
//    can be one contact for each of 6 faces
//    no contact if outside (possible if called from union/intersect)
//    delxyz = vector from nearest point on prism to x
// ------------------------------------------------------------------------- */
// 
// int cpuPrismSurfaceInterior(double *r, double *delx, double *dely, double *delz, double *radius, 
//          double **face, double *x, double *clo, double *chi, double cutoff, int *iwall, int *open_faces)
// {
//   int i;
//   double dot;
//   double *corner;
// 
//   // x is exterior to prism
// 
//   for (i = 0; i < 6; i++) {
//     if (i % 2) corner = chi;
//     else corner = clo;
//     dot = (x[0]-corner[0])*face[i][0] + (x[1]-corner[1])*face[i][1] +
//       (x[2]-corner[2])*face[i][2];
//     if (dot < 0.0) return 0;
//   }
// 
//   // x is interior to prism or on its surface
// 
//   int n = 0;
// 
//   for (i = 0; i < 6; i++) {
//     if (open_faces[i]) continue;
//     if (i % 2) corner = chi;
//     else corner = clo;
//     dot = (x[0]-corner[0])*face[i][0] + (x[1]-corner[1])*face[i][1] +
//       (x[2]-corner[2])*face[i][2];
//     if (dot < cutoff) {
//       r[n] = dot;
//       delx[n] = dot*face[i][0];
//       dely[n] = dot*face[i][1];
//       delz[n] = dot*face[i][2];
//       radius[n] = 0;
//       iwall[n] = i;
//       n++;
//     }
//   }
// 
//   return n;
// }
// 
// /* ----------------------------------------------------------------------
//    x is exterior to prism or on its surface
//    return (xp,yp,zp) = nearest pt to x that is on surface of prism
// ------------------------------------------------------------------------- */
// 
// void cpuPrismFindNearest(double *nearest, double *x, double **face, double **corners, double **tri, int *open_faces)
// {
//   int i,j,k,iface;
//   double xproj[3],xline[3];
//   double dot;
// 
//   double distsq = BIG;
// 
//   for (int itri = 0; itri < 12; itri++) {
//     iface = itri/2;
//     if (open_faces[iface]) continue;
//     i = tri[itri][0];
//     j = tri[itri][1];
//     k = tri[itri][2];
//     dot = (x[0]-corners[i][0])*face[iface][0] +
//       (x[1]-corners[i][1])*face[iface][1] +
//       (x[2]-corners[i][2])*face[iface][2];
//     xproj[0] = x[0] - dot*face[iface][0];
//     xproj[1] = x[1] - dot*face[iface][1];
//     xproj[2] = x[2] - dot*face[iface][2];
//     if (cpuInsideTri(xproj,corners[i],corners[j],corners[k],face[iface])) {
//       distsq = cpuClosest(x,xproj,nearest,distsq);
//     }
//     else {
//       cpuPointOnLineSegment(corners[i],corners[j],xproj,xline);
//       distsq = cpuClosest(x,xline,nearest,distsq);
//       cpuPointOnLineSegment(corners[j],corners[k],xproj,xline);
//       distsq = cpuClosest(x,xline,nearest,distsq);
//       cpuPointOnLineSegment(corners[i],corners[k],xproj,xline);
//       distsq = cpuClosest(x,xline,nearest,distsq);
//     }
//   }
// 
// }
// 
// /* ----------------------------------------------------------------------
//    one contact if 0 <= x < cutoff from outer surface of prism
//    no contact if inside (possible if called from union/intersect)
//    delxyz = vector from nearest point on prism to x
// ------------------------------------------------------------------------- */
// 
// int cpuPrismSurfaceExterior(double *r, double *delx, double *dely, double *delz, double *radius, 
//           double **face, double *corners, double **tri, double *x, double *clo, double *chi, double cutoff, 
//           int *iwall, int *open_faces)
// {
//   int i;
//   double dot;
//   double *corner;
//   double p[3];
// 
//   // x is far enough from prism that there is no contact
// 
//   for (i = 0; i < 6; i++) {
//     if (i % 2) corner = chi;
//     else corner = clo;
//     dot = (x[0]-corner[0])*face[i][0] + (x[1]-corner[1])*face[i][1] +
//       (x[2]-corner[2])*face[i][2];
//     if (dot <= -cutoff) return 0;
//   }
// 
//   // x is interior to prism
// 
//   for (i = 0; i < 6; i++) {
//     if (i % 2) corner = chi;
//     else corner = clo;
//     dot = (x[0]-corner[0])*face[i][0] + (x[1]-corner[1])*face[i][1] +
//       (x[2]-corner[2])*face[i][2];
//     if (dot <= 0.0) break;
//   }
// 
//   if (i == 6) return 0;
// 
//   // x is exterior to prism or on its surface
//   // xp,yp,zp = point on surface of prism that x is closest to
//   //            could be edge or corner pt of prism
//   // do not add contact point if r >= cutoff
//   cpuPrismFindNearest(p, x, face, corners, tri, open_faces);
//   cpuAddContact(r, radius, delx, dely, delz, x, p, 0);
//   radius[0] = 0;
//   iwall[0] = 0;
//   if (r[0] < cutoff) return 1;
//   return 0;
// }
// 
