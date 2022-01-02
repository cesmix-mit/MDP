#pragma once

#ifdef __cplusplus
extern "C" {
#endif

// double cpuRandomUniform(int *seed);
// double cpuRandomGaussian(int *seed, int *save, double *second);
// void cpuRandomResetSeed(int *seed, int *save, int ibase, double *coord);
// double cpuGamdev(const int ia, int *seed);

// int cpuLatticeOrthogonal(int *orientx, int *orienty, int *orientz);
// int cpuLatticeRightHanded(int *orientx, int *orienty, int *orientz);
// int cpuLatticeCollinear(double *a1, double *a2, double *a3);
// void cpuLatticeTransform(double *primitive, double *priminv, double *rotaterow, double *rotatecol, 
//         double *a1, double *a2, double *a3, int *orientx, int *orienty, int *orientz);
// void cpuLattice2Box(double *x, double *primitive, double *rotaterow, double *origin, 
//         double *spacing, double scale, int dim);
// void cpuBox2Lattice(double *x, double *priminv, double *rotatecol, double *origin, 
//         double *spacing, double scale, int dim);
// void cpuLatticeBbox(double *pmin, double *pmax, double *x, double *primitive, double *rotaterow, double *priminv, double *rotatecol,
//         double *origin, double *spacing, double scale, int flag, int dim);
    
// int cpuInsideTri(double *x, double *v1, double *v2, double *v3, double *norm);
// double cpuClosest(double *x, double *near, double *nearest, double dsq);
// void cpuAddContact(double *r, double *radius, double *delx, double *dely, double *delz, double *x, double *p, int n);
// void cpuRotate(double *x, double *point, double *runit, double angle);
// void cpuForwardTransform(double *x, double *dx, double *point, double *runit, double theta, int rotateflag, int moveflag);
// void cpuInverseTransform(double *x, double *dx, double *point, double *runit, double theta, int rotateflag, int moveflag);
    
// int cpuRegionSetup(double *point, double *runit, double *scale, double *axis, double *lattice, 
//         int scaleflag, int rotateflag, int moveflag);
// void cpuPointOnLineSegment(double *a, double *b, double *c, double *d);
// // int cpuRegionSurface(double *x,  double *dx, double *delx, double *dely, double *delz, 
// //         double *point, double *runit, double theta, double cutoff, int dynamic, int moveflag, int rotateflag, 
// //         int openflag, int interior);
// void cpuRegionSetVelocity(double *v, double *prev, double *omega, double *dx, double *point, double *rpoint, 
//         double *axis, double dt, double theta, int moveflag, int rotateflag, int ntimestep);
// void cpuRegionVelocityContact(double *vwall, double *x, double *v, double *omega, double *delx, double *dely, double *delz, 
//         double *rpoint,int moveflag, int rotateflag, int ic);
// void cpuBlockSetup(double *extent_lo, double *extent_hi, double **face, double ***corners, double *lo, double *hi, double *scale);
// int cpuBlockInside(double *x, double *lo, double *hi);
// int cpuBlockSurfaceInterior(double *r, double *delx, double *dely, double *delz, double *radius, 
//         double *x, double *lo, double *hi, double cutoff, int *iwall, int *open_faces);
// int cpuBlockSurfaceExterior(double *r, double *delx, double *dely, double *delz, double *radius, 
//         double *x, double *lo, double *hi, double cutoff, int *iwall, int *open_faces, int openflag);
// //double cpuBlockFindClosestPoint(double *x, double **face, double ***corners, double *lo, double *hi, double *c, int i);
// int cpuBlockInsideFace(double *xproj, double *lo, double *hi, int iface);
// int cpuPrismSetup(double **h, double **hinv, double **face, double *corners, double *a, double *b, double*c, 
//         double *clo, double *chi, double *extent_lo, double *extent_hi, double *lo, double *hi, double *tilt, double *scale,
//         int **tri, int *open_faces, int interior, int openflag);
// int cpuPrismInside(double *x, double **hinv, double *lo);
// int cpuPrismSurfaceInterior(double *r, double *delx, double *dely, double *delz, double *radius, 
//          double **face, double *x, double *clo, double *chi, double cutoff, int *iwall, int *open_faces);
// void cpuPrismFindNearest(double *nearest, double *x, double **face, double **corners, double **tri, int *open_faces);
// int cpuPrismSurfaceExterior(double *r, double *delx, double *dely, double *delz, double *radius, 
//           double **face, double *corners, double **tri, double *x, double *clo, double *chi, double cutoff, 
//           int *iwall, int *open_faces);

int cpuLattice(double *basis, double *primitive, double *rotaterow, double *priminv, double *rotatecol, double *origin, 
        double *spacing, double *a1, double *a2, double *a3, double &scale, int *orientx, int *orienty, int *orientz, 
        int style, int unit_style, int spaceflag, int dimension);
void cpuLatticeBoundingBox(double *lmin, double *lmax, double *bsublo, double *bsubhi, double *primitive, 
        double *rotaterow, double *priminv, double *rotatecol, double *origin, double *spacing, double scale);
int cpuLatticeCount(int nbasis, int ilo, int ihi, int jlo, int jhi, int klo, int khi);

void cpuSetGlobalBox(double *h, double *h_inv, double *boxlo_bound, 
        double *boxhi_bound, double *boxhi, double *boxlo, double *boxtilt, int triclinic);
void cpuSetLocalOrthBox(double *subhi, double *sublo, double *boxhi, double *boxlo, double *subhi_lamda, double *sublo_lamda, int dim);
void cpuShiftLocalOrthBox(double *subhi, double *sublo, double *boxhi, double *boxlo, double *epsilon, int *pbc, int dim);
void cpuLamda2Box(double *x, double *lambda, double *h, double *boxlo, int dim, int n);
void cpuBox2Lamda(double *lambda, double *x, double *h_inv, double *boxlo, int dim, int n);
//void cpuDomainBbox(double *bboxlo, double *bboxhi, double *lo_lamda, double *hi_lamda, double *boxlo, double *h, int dim);
void cpuShiftedSubbox(double *ssublo, double *ssubhi, double *boxlo, double *boxhi, 
        double *boxlo_lamda, double *boxhi_lamda, double *sublo, double *subhi, 
        double *sublo_lamda, double *subhi_lamda, double *epsilon, int *pbc, int triclinic);
void cpuBoundingSubbox(double *bsublo, double *bsubhi, double *sublo, double *subhi, 
        double *sublo_lamda, double *subhi_lamda, double *boxlo, double *h, int triclinic);
// void cpuPBCOrthBox(double *x, double *v, int *image, double *hi, double *lo, 
//         double *h, double *h_rate, int *pbc, int vdeform, int dim, int nlocal);
// void cpuPBCTiltBox(double *x, double *v, int *image, double *boxlo, double *h, double *h_inv, double *h_rate, 
//         double *hi_lambda, double *lo_lambda, int *pbc, int vdeform, int dim, int nlocal);
void cpuPBC(double *x, double *v, int *image, double *boxhi, double *boxlo, double *hi_lambda, double *lo_lambda,  
        double *h, double *h_inv, double *h_rate, int *pbc, int vdeform, int triclinic, int dim, int nlocal);
// void cpuInsideOrthBox(double *inside, double *x, double *lo, double *hi, int dim, int n);
// void cpuInsideTiltBox(double *inside, double *x, double *lo_lamda, double *hi_lamda, 
//         double *boxlo, double *h_inv, int dim, int n);
void cpuInsideBox(double *inside, double *x, double *boxlo, double *boxhi, double *lo_lamda, double *hi_lamda, 
        double *h_inv, int triclinic, int dim, int n);
//void cpuMinimumImage(double *dp, double *h, int *pbc, int triclinic, int dim);
// void cpuMinimumImageOrthBox(double *dp, double *h, int *pbc, int dim, int n);
// void cpuMinimumImageTiltBox(double *dp, double *h, int *pbc, int dim, int n);
void cpuMinimumImage(double *dp, double *h, int *pbc, int triclinic, int dim, int n);
//void cpuUnmap(double *y, double *x, double *h, int *image, int triclinic, int dim);
// void cpuUnmapOrthBox(double *y, double *x, double *h, int *image, int dim, int n);
// void cpuUnmapTiltBox(double *y, double *x, double *h, int *image, int dim, int n);
void cpuUnmap(double *y, double *x, double *h, int *image, int triclinic, int dim, int n);

// void cpuPackIntProperty(double *buf, int *prop, int *ilist, 
// int m, int mvalues, int n, int nvalues, int inum);
// void cpuPackIntProperty(double *buf, int *prop, int *type, int *ilist, 
//          int n, int nvalues, int inum);
// void cpuPackFloatProperty(double *buf, double *prop, int *ilist, 
//          int m, int mvalues, int n, int nvalues, int inum);
// void cpuPackFloatProperty(double *buf, double *prop, double a, double b, int *ilist, 
//          int m, int mvalues, int n, int nvalues, int inum);
// void cpuPackFloatProperty(double *buf, double *prop, int *type, int *ilist, 
//          int n, int nvalues, int inum);
double cpuComputeMass(double *amass, double *mass, int *type, int *ilist, int inum);
void cpuComputeXCM(double *xcm, double *x, double *mass, double *box, double masstotal, 
        int *ilist, int *type, int *image, int triclinic, int dim, int inum);
void cpuComputeVCM(double *vcm, double *v, double *mass, double masstotal, 
        int *ilist, int *type, int dim, int inum);
double cpuComputeGyration(double *xcm, double *x, double *mass, double *box, double masstotal, 
        int *ilist, int *type, int *image, int triclinic, int dim, int inum);
void cpuComputeAngmom(double *p, double *xcm, double *x, double *v, double *mass, double *box, 
        int *ilist, int *type, int *image, int triclinic, int dim, int inum);
void cpuComputeTorque(double *tlocal, double *xcm, double *x, double *f, double *box, 
        int *ilist, int *image, int triclinic, int dim, int inum);
void cpuComputeInertia(double *ione, double *xcm, double *x, double *mass, double *box, 
        int *ilist, int *type, int *image, int triclinic, int dim, int inum);
double cpuComputeNVTEnergy(double *tarray, double *eta, double *eta_mass, double *eta_dot, int mtchain);
double cpuComputeTempScalar(double *v, double *mass, double tfactor, int *type, int *ilist, 
         int dim, int inum);
void cpuComputeTempSymTensor(double *t, double *v, double *mass, double tfactor, int *type, int *ilist, 
         int dim, int inum);
double cpuComputePressureScalar(double *virial, double volume, double temp, double tempdof, 
        double boltz, double nktv2p, int dim);
void cpuComputePressureSymTensor(double *vector, double *virial, double *ke_tensor, 
        double volume, double nktv2p, int dim);
void cpuComputeHeatFlux(double *vector, double *ke, double *pe, double *stress, double *v, 
        double nktv2p, int *ilist,  int pressatomflag, int dim, int inum);
void cpuComputeKEAtom(double *ke, double *mass,  double *v, 
        double mvv2e, int *type, int *ilist,  int dim, int inum);
void cpuComputeStressAtom(double *stress, double *mass, double *vatom, double *v, 
        double mvv2e, double nktv2p, int *type, int *ilist,  int vflag, int keflag, int dim, int inum);
void cpuComputeCentroidStressAtom(double *stress, double *mass, double *cvatom, double *v, 
        double mvv2e, double nktv2p, int *type, int *ilist, int vflag, int keflag, int dim, int inum);
void cpuComputeOrientOrderAtom(double* qnarray, double *x, double *rlist, double *cglist, double *fac, 
        double *qnm_r, double *qnm_i, double *distsq, double cutsq, double MY_EPSILON, double QEPSILON, double MY_4PI, int *neighlist, 
        int *neighnum,  int *ilist, int *qlist, int *nearest,  int nqlist, int qmax, int wlflag,  
        int wlhatflag,  int qlcompflag, int iqlcomp, int qlcomp, int nnn, int jnum, int dim, int inum); 
// void cpuComputeCoordAtomCutoff(int *cvec, double *x, double *rcutsq, int *type, 
//         int *ilist, int *neighlist, int *neighnum, int *typelo, int *typehi, 
//          int *jgroupbit, int dim, int ntypes, int jnum, int inum);
void cpuComputeCoordAtomCutoff(int *carray, double *x, double *rcutsq, int *type, 
        int *ilist, int *neighlist, int *neighnum, int *typelo, int *typehi, 
         int *jgroupbit, int ncol, int dim, int ntypes, int jnum, int inum);
void cpuComputeCoordAtomOrient(int *cvec, double *x, double *rcutsq, double *normv, double threshold, 
        int *type, int *ilist, int *neighlist, int *neighnum, 
         int *jgroupbit, int dim, int ntypes, int nqlist, int ncol, int l, int jnum, int inum);
void cpuComputeMSD(double *vector, double *x, double *xoriginal, double *h, double *xcm,
         int *ilist, int *image, int naverage, int avflag, int triclinic, int nmsd, int dim,  int inum);
void cpuComputeVACF(double *vacf, double *v, double *voriginal, 
         int *ilist, int nvacf, int dim,  int inum);

int cpuSetAtomType(double *x, double fraction, int *atomtype,
         int *seed, int *save, int seed0, int newtype, int dim, int nlocal);
void cpuAtomInside(int *inside, double *x, double *h_inv, double *boxlo, double *lo_lamda,
        double *hi_lamda, int dim, int n);
int cpuAtomAdd(double *y, int *atomtype, double *x, double *h_inv, double *boxlo, double *lo_lamda,
        double *hi_lamda, int *type, int dim, int nlocal, int n);
void cpuAtomLattice(double *y, int *atomtype, double *basis, double *primitive, double *rotaterow, double *origin, 
        double *latticespacing, double scale, int *basistype, int nbasis, int nlocal, 
        int ilo, int ihi, int jlo, int jhi, int klo, int khi, int dim);

void cpuFixSetForce(double *x, double *v, double *f, double *eatom, double *vatom, double *fparam, 
        int *iparam, int *ilist, int eflag_atom, int vflag_atom, int dim, int inum);        
void cpuFixLineForce(double *x, double *v, double *f, double *eatom, double *vatom, double *fparam, 
        int *iparam, int *ilist, int eflag_atom, int vflag_atom, int dim, int inum);
void cpuFixPlaneForce(double *x, double *v, double *f, double *eatom, double *vatom, double *fparam, 
        int *iparam, int *ilist, int eflag_atom, int vflag_atom, int dim, int inum);        
// void cpuFixAddForce(double *x, double *v, double *f, double *eatom, double *vatom, double *fparam, 
//         int *iparam, int *ilist, int eflag_atom, int vflag_atom, int dim, int inum);
void cpuFixAddForce(double *x, double *v, double *f, double *eatom, double *vatom, double *fparam, double *box,
        int *iparam, int *ilist, int *image, int triclinic, int eflag_atom, int vflag_atom, int dim, int inum);
void cpuFixAveForce(double *x, double *v, double *f, double *eatom, double *vatom, double *fparam, 
        int *iparam, int *ilist, int eflag_atom, int vflag_atom, int dim, int inum);
void cpuFixDragForce(double *x, double *v, double *f, double *eatom, double *vatom, double *fparam, double *box,
        int *iparam, int *ilist, int *pbc, int triclinic, int eflag_atom, int vflag_atom, int dim, int inum);
void cpuFixWallReflect(double *x, double *v, double *f, double *eatom, double *vatom, double *fparam, 
        int *iparam, int *ilist, int eflag_atom, int vflag_atom, int dim, int inum);
void cpuFixWallHarmonic(double *x, double *v, double *f, double *eatom, double *vatom, double *fparam, 
        int *iparam, int *ilist, int eflag_atom, int vflag_atom, int dim, int inum);
void cpuFixWallLJ93(double *x, double *v, double *f, double *eatom, double *vatom, double *fparam, 
        int *iparam, int *ilist, int eflag_atom, int vflag_atom, int dim, int inum);
void cpuFixWallLJ126(double *x, double *v, double *f, double *eatom, double *vatom, double *fparam, 
        int *iparam, int *ilist, int eflag_atom, int vflag_atom, int dim, int inum);
void cpuFixWallLJ1043(double *x, double *v, double *f, double *eatom, double *vatom, double *fparam, 
        int *iparam, int *ilist, int eflag_atom, int vflag_atom, int dim, int inum);
void cpuFixWallMorse(double *x, double *v, double *f, double *eatom, double *vatom, double *fparam, 
        int *iparam, int *ilist, int eflag_atom, int vflag_atom, int dim, int inum);

void cpuVelocityZeroMomentum(double *v, double *vcm, int dim, int nlocal);
void cpuVelocityZeroRotation(double *x, double *v, double *box, double *xcm, double *omega, 
        int *image, int triclinic, int dim, int nlocal);
void cpuVelocityCreate(double *x, double *v, double *mass, double *second, 
         int *seed, int *save, int *map, int *type, int seed0, int sum_flag, int dist_flag, int loop_flag, 
        int dim, int mpiRank, int nlocal, int natoms);
// void cpuVelocityCreate(double *x, double *v, double *mass, double *second, double *omega, 
//         double *box, double *xcm, double *vcm, double t_desired, double t_current, int *seed, int *save, int *map, int *image, 
//         int *type, int sum_flag, int dist_flag, int loop_flag, int rotation_flag, 
//         int momentum_flag, int triclinic, int dim, int mpiRank, int nlocal, int natoms);
void cpuVelocitySet(double *v, double *vext, int *vdim,
        int sum_flag, int dim, int nlocal);
void cpuVelocityRamp(double *x, double *v, double *v_lo, double *v_hi, double *coord_lo, double *coord_hi,
        int *coord_dim, int *v_dim, int sum_flag, int dim, int nlocal);
// void cpuNVEInitialIntegrate(double *x, double *v, double *f, double *mass, double dtf, double dtv,
//         int *type, int *ilist, int dim, int inum);
// void cpuNVEFinalIntegrate(double *x, double *v, double *f, double *mass, double dtf,
//         int *type, int *ilist, int dim, int inum);
// void cpuNVELimitInitialIntegrate(double *x, double *v, double *f, double *mass, double dtf, double dtv,
//         double vlimitsq, int *type, int *ilist, int dim, int inum);
// void cpuNVELimitFinalIntegrate(double *x, double *v, double *f, double *mass, double dtf,
//         double vlimitsq, int *type, int *ilist, int dim, int inum);
void cpuSetVelocityInitialIntegrate(double *x, double *v, double *f, double *mass, double *fparam,
        double dtf, double dtv, int *type, int *ilist, int *iparam, int dim, int inum);
void cpuSetVelocityFinalIntegrate(double *x, double *v, double *f, double *mass, 
        double dtf, int *type, int *ilist, int *iparam, int dim, int inum);
// void cpuVelocityFactor(double *v, double factor_eta, int *ilist, 
//         int biasflag, int dim, int inum);
// void cpuNoseHooverThermostat(double *v, double *dtarray, double *tarray, double *eta_mass, double *eta, 
//         double *eta_dot, double *eta_dotdot, int *ilist, int eta_mass_flag, int biasflag, int mtchain, 
//         int nc_tchain, int dim, int inum);
// void cpuNVTInitialIntegrate(double *x, double *v, double *f, double *mass, double *dtarray, double *tarray,
//         double *eta_mass, double *eta, double *eta_dot, double *eta_dotdot, int *type, int *ilist, 
//         int eta_mass_flag, int biasflag, int mtchain, int nc_tchain, int dim, int inum);
// void cpuNVTFinalIntegrate(double *x, double *v, double *f, double *mass, double *dtarray, double *tarray,
//         double *eta_mass, double *eta, double *eta_dot, double *eta_dotdot, int *type, int *ilist, 
//         int eta_mass_flag, int biasflag, int mtchain, int nc_tchain, int dim, int inum);
void cpuInitialIntegrate(double *x, double *v, double *f, double *mass, double *dtarray, double *tarray,
        double *eta_mass, double *eta, double *eta_dot, double *eta_dotdot, double vlimitsq, int *type, int *ilist,
        int eta_mass_flag, int biasflag, int mtchain, int nc_tchain, int mode, int dim, int inum);
void cpuFinalIntegrate(double *x, double *v, double *f, double *mass, double *dtarray, double *tarray,
        double *eta_mass, double *eta, double *eta_dot, double *eta_dotdot, double vlimitsq, int *type, int *ilist, 
        int eta_mass_flag, int biasflag, int mtchain, int nc_tchain, int mode, int dim, int inum);
// double cpuBerendsenThermostat(double *v, double *mass, double *dtarray, double *tarray, double energy, 
//         int *type, int *ilist, int biasflag, int dim, int inum);
// double cpuRescaleThermostat(double *v, double *mass, double *dtarray, double *tarray, double energy, 
//         int *type, int *ilist, int biasflag, int dim, int inum);
// double cpuCsvrThermostat(double *v, double *mass, double *dtarray, double *tarray, double *second, 
//         double energy, int *type, int *ilist, int *seed, int *save, int biasflag, int dim, int inum);
double cpuVelocityRescalingThermostat(double *v, double *mass, double *dtarray, double *tarray, double *second, 
        double energy, int *type, int *ilist, int *seed, int *save, int biasflag, int mode, int dim, int inum);

#ifdef __cplusplus
}
#endif


