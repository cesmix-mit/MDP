/***************************************************************************                               
                    Molecular Dynamics Potentials (MDP)
                           CESMIX-MIT Project  
 
 Contributing authors: Ngoc-Cuong Nguyen (cuongng@mit.edu, exapde@gmail.com)
 ***************************************************************************/

#ifndef __GPUCORE_H__
#define __GPUCORE_H__

template <typename T> void gpuKron(T *C, T *A, T *B, int M1, int M2);
template <typename T> void gpuKron(T *C, T *A, T *B, int M1, int N1, int M2, int N2);
void gpuIndexPermute12(int *index, int I1, int I2, int I3);
void gpuIndexPermute13(int *index, int I1, int I2, int I3, int I4);
void gpuIndexPermute23(int *index, int I1, int I2, int I3, int I4);
template <typename T> void gpuPermute(T *B, T *A, int *index, int N);
template <typename T> void gpuPermuteSharedMem(T *B, T *A, int *index, int N, int BLOCKDIM);
template <typename T> void gpuPermute12(T *B, T *A, int I1, int I2, int I3);
template <typename T> void gpuPermute13(T *B, T *A, int I1, int I2, int I3, int I4);
template <typename T> void gpuPermute23(T *B, T *A, int I1, int I2, int I3, int I4);
template <typename T> void gpuGEMM(T *C, T *A, T *B, int I, int J, int K);
template <typename T> void gpuTensorGEMM(T *C, T *A1, T *A2, T *A3, T *B, T* Ctmp, int *index, 
        int M1, int N1, int M2, int N2, int M3, int N3, int K);

template <typename T> void gpuGetArrayAtIndex(T *y, T *x, int *ind, int n);
template <typename T> void gpuPutArrayAtIndex(T *y, T *x, int *ind, int n);
template <typename T> void gpuArrayPlusXAtIndex(T *y, T *x, int *ind, int n);
template <typename T> void gpuArrayMinusXAtIndex(T *y, T *x, int *ind, int n);
template <typename T> void gpuArrayAXPYAtIndex(T *y, T *x, T a, int *ind, int n);

template <typename T> void gpuArrayTranspose(T *A, T *B, int m, int n);
template <typename T> void gpuArrayPlusAtColumnIndex(T *A, T *B, int *colind, int m, int n);
template <typename T> void gpuArrayMinusAtColumnIndex(T *A, T *B, int *colind, int m, int n);
template <typename T> void gpuGetArrayAtColumnIndex(T *A, T *B, int *colind, int m, int n);
template <typename T> void gpuArrayTransposeAtColumnIndex(T *A, T *B, int *colind, int m, int n);
template <typename T> void gpuGetArrayAtRowIndex(T *A, T *B, int *rowind, int m, int n);
template <typename T> void gpuArrayTransposeAtRowIndex(T *A, T *B, int *rowind, int m, int n);
template <typename T> void gpuArrayRowSum(T *y, T *x, int m, int n);
template <typename T> void gpuArrayRowSquareSum(T *y, T *x, int m, int n);
template <typename T> void gpuArrayDistSquareSum(T *y, T *x1, T *x2, int m, int n);

void gpuTripletnum(int *a, int *b, int n);
void gpuQuadrupletnum(int *a, int *b, int n);
void gpuArrayFill(int *a, int start, int n);
template <typename T> void gpuPrint2DArray(T* a, int m, int n);
template <typename T> void gpuPrint3DArray(T* a, int m, int n, int p);
template <typename T> void gpuArraySetValue(T *y, T a, int n);
template <typename T> void gpuArraySetValueAtIndex(T *y, T a, int n);
template <typename T> T gpuArrayGetValueAtIndex(T *y, int n);
template <typename T> void gpuArrayAddScalar(T *y, T a, int n);
template <typename T> void gpuArrayMultiplyScalar(T *y, T a, int n);

template <typename T> T gpuArraySum(T *a, T *b, int n);
template <typename T> T gpuArrayMax(T *a, T *b, int n);
template <typename T> T gpuArrayMin(T *a, T *b, int n);
template <typename T> void gpuArraySumEveryColumn(T *b, T *a, int m, int n);
template <typename T> void gpuArrayMaxEveryColumn(T *b, T *a, int m, int n);
template <typename T> void gpuArrayMinEveryColumn(T *b, T *a, int m, int n);
template <typename T> void gpuArraySumEveryRow(T *b, T *a, T *c, int m, int n);
template <typename T> void gpuArrayMaxEveryRow(T *b, T *a, T *c, int m, int n);
template <typename T> void gpuArrayMinEveryRow(T *b, T *a, T *c, int m, int n);

template <typename T> void gpuArrayCopy(T *y, T *x, int n);
template <typename T> void gpuArrayMinus(T *y, T *x, int n);
template <typename T> void gpuArrayAbs(T *y, T *x, int n);
template <typename T> void gpuArraySqrt(T *y, T *x, int n);
template <typename T> void gpuArraySin(T *y, T *x, int n);
template <typename T> void gpuArrayCos(T *y, T *x, int n);
template <typename T> void gpuArrayTan(T *y, T *x, int n);
template <typename T> void gpuArrayAsin(T *y, T *x, int n);
template <typename T> void gpuArrayAcos(T *y, T *x, int n);
template <typename T> void gpuArrayAtan(T *y, T *x, int n);
template <typename T> void gpuArraySinh(T *y, T *x, int n);
template <typename T> void gpuArrayCosh(T *y, T *x, int n);
template <typename T> void gpuArrayTanh(T *y, T *x, int n);
template <typename T> void gpuArrayAsinh(T *y, T *x, int n);
template <typename T> void gpuArrayAcosh(T *y, T *x, int n);
template <typename T> void gpuArrayAtanh(T *y, T *x, int n);
template <typename T> void gpuArrayExp(T *y, T *x, int n);
template <typename T> void gpuArrayLog(T *y, T *x, int n);
template <typename T> void gpuArrayCeil(T *y, T *x, int n);
template <typename T> void gpuArrayFloor(T *y, T *x, int n);
template <typename T> void gpuArrayErf(T *y, T *x, int n);
template <typename T> void gpuArrayErfc(T *y, T *x, int n);
template <typename T> void gpuArraySquare(T *y, T *x, int n);
template <typename T> void gpuArrayPower(T *y, T *x, int p, int n);

template <typename T> void gpuArrayMultiplyScalarDiagonal(T *C, T a, int n);
template <typename T> void gpuArrayAddVectorToDiagonal(T *C, T *x, T a, int n);
template <typename T> void gpuArrayRowAverage(T *a, T *b, int M, int N);
template <typename T> void gpuArrayAXPB(T *y, T *x, T a, T b, int n);
template <typename T> void gpuArrayAXPBY(T *z, T *x, T *y, T a, T b, int n);
template <typename T> void gpuArrayAXY(T *s, T *x, T *y, T a, int n);
template <typename T> void gpuArrayAXYZ(T *s, T *x, T *y, T *z, T a, int n);
template <typename T> void gpuArrayAXYPBZ(T *s, T *x, T *y, T *z, T a, T b, int n);
template <typename T> void gpuArrayAdd3Vectors(T *s, T *x, T *y, T *z, T a, T b, T c, int n);
template <typename T> void gpuArrayAdd3Vector(T *a, T *b, T *c, T *d, Int n);
template <typename T> void gpuArrayExtract(T *un, T *u, int I, int J, int K, 
        int i1, int i2, int j1, int j2, int k1, int k2);
template <typename T> void gpuArrayInsert(T *u, T *un, int I, int J, int K, 
        int i1, int i2, int j1, int j2, int k1, int k2);
template <typename T> void gpuArrayGemmBatch(T *C, T *A, T *B, int I, int J, int K, int S);
template <typename T> void gpuArrayGemmBatch1(T *C, T *A, T *B, int I, int J, int K, int S);
template <typename T> void gpuArrayDG2CG(T *ucg, T *udg, int *cgent2dgent, int *rowent2elem, int nent);
template <typename T> void gpuArrayDG2CG2(T *ucg, T *udg, int *colent2elem, int *rowent2elem, int nent, int npe);

template <typename T> void gpuArrayABPXYZ(T *fij, T *c3ij, T *eij, T *d3ij, T *g3ij, int dim, int n);
void gpuCumsum(int *d_out, int *d_in, int *d_sums, int *d_incr, int length);

template <typename T> void gpuCart2Sphere(T *the, T *phi, T *r, T *x, T *y, T *z, int N);
template <typename T> void gpuCart2SphereDeriv(T *the, T *phi, T *r, T *thex, T *they, T *thez, T *phix, T *phiy, T *phiz, T *rx, T *ry, T *rz, T *x, T *y, T *z, int N);
template <typename T> void gpuSphere2Cart(T *x, T *y, T *z, T *the, T *phi, T *r, int N);
template <typename T> void gpuEuler2Rotm(T *R11, T *R12, T *R13, T *R21, 
                T *R22, T *R23, T *R31, T *R32, T *R33, T *alpha, T *beta, T *gamma, int N);
template <typename T> void gpuRotc(T *X, T *Y, T *Z, T *R, T *x, T *y, T *z, int N);


int gpuUniqueSort(int *b, int *c, int *d, int *e, int *a, int *p, int *t, int *q, int n);
int gpuFindAtomType(int *tlist, int* ilist, int *atomtype, int *p, int *q, int typei, int inum);

//********************************* simulation Box ****************************************//
template <typename T> void gpuLamda2Box(T *x, T *lambda, T *h, T *boxlo, int dim, int n);
template <typename T> void gpuBox2Lamda(T *lambda, T *x, T *h_inv, T *boxlo, int dim, int n);
template <typename T> void gpuInsideBox(T *inside, T *x, T *boxlo, T *boxhi, T *lo_lamda, T *hi_lamda, 
        T *h_inv, int triclinic, int dim, int n);
template <typename T> void gpuMinimumImage(T *dp, T *h, int *pbc, int triclinic, int dim, int n);
template <typename T> void gpuUnmap(T *y, T *x, T *h, int *image, int triclinic, int dim, int n);

//********************************* Neighbor Lists ****************************************//
template <typename T> void gpuAtomList2D(int *alist,  int *inside, int *glistnumsum, int *glistnum, 
        int *d_sums, int *d_incr, T *x, T *pimages, T *wc, T *s2rmap, int inum, int pnum, int dim);
template <typename T> void gpuAtomList3D(int *alist,  int *inside, int *glistnumsum, int *glistnum, 
        int *d_sums, int *d_incr, T *x, T *pimages, T *wc, T *s2rmap, int inum, int pnum, int dim);
template <typename T> void gpuCellList2D(int *clist, T *xi, T *eta1, T *eta2, T *eta3, T *s2rmap, 
        int *nc, int inum, int natom, int dim);
template <typename T> void gpuCellList3D(int *clist, T *xi, T *eta1, T *eta2, T *eta3, T *s2rmap, 
        int *nc, int inum, int natom, int dim);
void gpuCell2AtomList(int *c2alist, int *c2anumsum, int *c2anum, int *clist, int *d_temp, int natom, int ncell);
template <typename T> void gpuFullNeighborList2D(int *neighborlist, int *neighbornum, T *x, T *rcutsq, 
        int anum, int inum, int jnum, int dim);
template <typename T> void gpuFullNeighborList3D(int *neighborlist, int *neighbornum, T *x, T *rcutsq, 
        int anum, int inum, int jnum, int dim);
template <typename T> void gpuFullNeighborList2D(int *neighborlist, int *neighbornum, T *x, T *rcutsq,
        int *alist, int *clist, int *c2alist, int *c2anumsum, int *nc, int inum, int jnum, int dim);
template <typename T> void gpuFullNeighborList3D(int *neighborlist, int *neighbornum, T *x, T *rcutsq,
        int *alist, int *clist, int *c2alist, int *c2anumsum, int *nc, int inum, int jnum, int dim);

template <typename T> void gpuNeighSingles(T *xi, T *qi, T *x, T *q, int *ai, 
      int *ti, int *ilist, int *atomtype, int inum, int ncq, int dim);
template <typename T> void gpuFullNeighPairList(int *pairnum, int *pairlist, T *x, T *rcutsq, int *ilist, int *neighlist, 
        int *neighnum, int inum, int jnum, int dim);
template <typename T> void gpuFullNeighPairList(int *pairnum, int *pairlist, T *x, T *rcutsq,  int *atomtype, int *ilist, int *alist, 
        int *neighlist, int *neighnum, int inum, int jnum, int typej, int dim);
template <typename T> void gpuFullNeighPairList(int *pairnum, int *pairlist, T *x, T *rcutsq,  int *atomtype, int *ilist, int *alist, 
        int *neighlist, int *neighnum, int inum, int jnum, int typej, int typek, int dim);
template <typename T> void gpuFullNeighPairList(int *pairnum, int *pairlist, T *x, T *rcutsq,  int *atomtype, int *ilist, int *alist, 
        int *neighlist, int *neighnum, int inum, int jnum, int typej, int typek, int typel, int dim);
template <typename T> void gpuHalfNeighPairList(int *pairnum, int *pairlist, T *x, T *rcutsq, int *ilist, int *alist, 
        int *neighlist, int *neighnum, int inum, int jnum, int dim);
template <typename T> void gpuHalfNeighPairList(int *pairnum, int *pairlist, T *x, T *rcutsq, int *atomtype, int *ilist, 
        int *alist, int *neighlist, int *neighnum, int inum, int jnum, int typej, int dim);
template <typename T> void gpuNeighPairs(T *xij, T *qi, T *qj, T *x, T *q, int *ai, int *aj,  
      int *ti, int *tj, int *pairnum, int *pairlist, int *pairnumsum, int *ilist, int *alist,  
      int *atomtype, int inum, int jnum, int ncq, int dim);

void gpuNeighTripletList(int *tripletlist, int *tripletnumsum, int *pairnum, 
        int *pairlist, int *ilist, int *alist, int inum, int jnum);
template <typename T> void gpuNeighTriplets(T *xij, T *xik, T *qi, T *qj, T *qk, T *x, T *q, int *ai, int *aj, int *ak,  
      int *ti, int *tj, int *tk, int *tripletnum, int *tripletlist, int *tripletnumsum, int *ilist, int *alist,  
      int *atomtype, int inum, int ncq, int dim);
template <typename T> void gpuNeighTripletList(int *tripletnum, int *tripletlist, T *x, T *rcutsq, int *pairnum, int *pairnumsum, 
        int *pairlist, int *atomtype, int *ilist, int *alist, int *neighlist, int *neighnum, int inum, 
        int jnum, int typek, int dim);
template <typename T> void gpuNeighTriplets(T *xij, T *xik, T *qi, T *qj, T *qk, T *x, T *q, int *ai, int *aj, int *ak,  
      int *ti, int *tj, int *tk, int *tripletnum, int *tripletlist, int *tripletnumsum, int *alist,  
      int *atomtype, int jnum, int ijnum, int ncq, int dim);
void gpuNeighQuadrupletList(int *quadrupletlist, int *quadrupletnumsum, int *pairnum, 
        int *pairlist, int *ilist, int *alist, int inum, int jnum);
template <typename T> void gpuNeighQuadruplets(T *xij, T *xik, T *xil, T *qi, T *qj, T *qk, T *ql, T *x, T *q, 
      int *ai, int *aj, int *ak, int *al, int *ti, int *tj, int *tk, int *tl, int *quadrupletnum, int *quadrupletlist, 
      int *quadrupletnumsum, int *ilist, int *alist,  int *atomtype, int inum, int ncq, int dim);

//********************************* Per-Atom Enegry, Force, Virial Tally ****************************************//
template <typename T> void gpuSingleDecomposition(T *e, T *ei, int *ai, int inum, int dim);
template <typename T> void gpuFullNeighPairDecomposition(T *e, T *eij, int *ai, int ijnum, int dim);
template <typename T> void gpuCenterAtomPairDecomposition(T *e, T *eij, int *ilist, int *anumsum, int inum, int dim);
template <typename T> void gpuHalfNeighPairDecomposition(T *e, T *eij, int *ai, int *aj, int ijnum, int dim);
template <typename T> void gpuNeighborAtomPairDecomposition(T *e, T *eij, int *jlist, int *bnumsum, int *index, int jnum, int dim);
template <typename T> void gpuTripletDecomposition(T *e, T *eijk, int *ai, int *aj, int *ak, int ijknum, int dim);
template <typename T> void gpuCenterAtomTripletDecomposition(T *e, T *eijk, int *ilist, int *anumsum, int inum, int dim);
template <typename T> void gpuNeighborAtomTripletDecomposition(T *e, T *eij, int *jlist, int *bnumsum, int *index, int jnum, int dim);
template <typename T> void gpuQuadrupletDecomposition(T *e, T *eijkl, int *ai, int *aj, int *ak, int *al, int ijklnum, int dim);
template <typename T> void gpuCenterAtomQuadrupletDecomposition(T *e, T *eijkl, int *ilist, int *anumsum, int inum, int dim);
template <typename T> void gpuNeighborAtomQuadrupletDecomposition(T *e, T *eij, int *jlist, int *bnumsum, int *index, int jnum, int dim);

template <typename T> void gpuSingleDecomposition(T *e, T *f, T *ei, T *fi, int *ai, int inum, int dim);
template <typename T> void gpuCenterAtomPairDecomposition(T *e, T *f, T *eij, T *fij, int *ilist, int *anumsum, int inum, int dim);
template <typename T> void gpuNeighborAtomPairDecomposition(T *e, T *f, T *eij, T *fij, int *jlist, int *bnumsum, int *index, int jnum, int dim);
template <typename T> void gpuFullNeighPairDecomposition(T *e, T *f, T *eij, T *fij, int *ai, int ijnum, int dim);
template <typename T> void gpuHalfNeighPairDecomposition(T *e, T *f, T *eij, T *fij, int *ai, int *aj, int ijnum, int dim);
template <typename T> void gpuTripletDecomposition(T *e, T *f, T *eijk, T *fij, T *fik, int *ai, int *aj, int *ak, int ijknum, int dim);
template <typename T> void gpuCenterAtomTripletDecomposition(T *e, T *f, T *eijk, T *fij, T *fik, int *ilist, int *anumsum, int inum, int dim);
template <typename T> void gpuNeighborAtomTripletDecomposition(T *e, T *f, T *eij, T *fij, int *jlist, int *bnumsum, int *index, int jnum, int dim);
template <typename T> void gpuQuadrupletDecomposition(T *e, T *f, T *eijkl, T *fij, T *fik, T *fil, 
        int *ai, int *aj, int *ak, int *al, int ijklnum, int dim);
template <typename T> void gpuCenterAtomQuadrupletDecomposition(T *e, T *f, T *eijkl, T *fij, T *fik, T *fil, 
        int *ilist, int *anumsum, int inum, int dim);
template <typename T> void gpuNeighborAtomQuadrupletDecomposition(T *e, T *f, T *eij, T *fij, int *jlist, int *bnumsum, int *index, int jnum, int dim);

template <typename T> void gpuHalfForceDecomposition(T *f, T *fij, int *ai, int *aj, int dim, int ijnum);
template <typename T> void gpuIAtomDecomposition(T *f, T *fij, int *ilist, int *anumsum, int dim, int inum);
template <typename T> void gpuJAtomDecomposition(T *f, T *fij, int *jlist, int *bnumsum, int *index, int dim, int jnum);
template <typename T> void gpuTripletForceDecomposition(T *f3ik, T *eij, T *d3ij, 
        int *tripletnum, int *tripletnumsum, int npairs, int dim);

// template <typename T> void gpuVirialPairTally(T *vatom, T *fij, T *rij, T factor, 
//         int *ai, int dim, int ijnum);
// template <typename T> void gpuVirialPairTally(T *vatom, T *fij, T *rij, T factor, 
//         int *ai, int *aj, int dim, int ijnum);
// template <typename T> void gpuVirialTripletTally(T *vatom, T *fij, T *fik, T *rij, T *rik, T factor, 
//         int *ai, int *aj, int *ak, int dim, int ijnum);
// template <typename T> void gpuVirialQuadrupletTally(T *vatom, T *fij, T *fik, T *fim, 
//         T *rij, T *rik, T *rim, T factor, int *ai, int *aj, int *ak, int *am, int dim, int ijnum);

template <typename T> void gpuVirialPairTally(T *vatom, T *fij, T *rij, T factor, 
        int *ai, int dim, int inum, int ijnum);
template <typename T> void gpuVirialPairTally(T *vatom, T *fij, T *rij, T factor, 
        int *ai, int *aj, int dim, int inum, int ijnum);
template <typename T> void gpuVirialTripletTally(T *vatom, T *fij, T *fik, T *rij, 
        T *rik, T factor, int *ai, int *aj, int *ak, int dim, int inum, int ijnum);
template <typename T> void gpuVirialQuadrupletTally(T *vatom, T *fij, T *fik, T *fim, 
        T *rij, T *rik, T *rim, T factor, int *ai, int *aj, int *ak, int *am, int dim, int inum, int ijnum);

//***************************  Compute Outputs ***************************************//
template <typename T> void gpuPackIntProperty(T *buf, int *prop, int *ilist, 
         int m, int mvalues, int n, int nvalues, int inum);
template <typename T> void gpuPackIntProperty(T *buf, int *prop, int *type, int *ilist, 
         int n, int nvalues, int inum);
template <typename T> void gpuPackFloatProperty(T *buf, T *prop, int *ilist, 
         int m, int mvalues, int n, int nvalues, int inum);
template <typename T> void gpuPackFloatProperty(T *buf, T *prop, T a, T b, int *ilist, 
         int m, int mvalues, int n, int nvalues, int inum);
template <typename T> void gpuPackFloatProperty(T *buf, T *prop, int *type, int *ilist, 
         int n, int nvalues, int inum);

template <typename T> T gpuComputeMass(T *amass, T *mass, T *tmp, int *type, int *ilist, int inum);
template <typename T> void gpuComputeXCM(T *xcm, T *axcm, T *x, T *tmp, T *mass, T *box, 
        T masstotal, int *ilist, int *type, int *image, int triclinic, int dim, int inum);
template <typename T> void gpuComputeVCM(T *vcm, T *avcm, T *v, T *tmp, T *mass,
        T masstotal, int *ilist, int *type, int dim, int inum);
template <typename T> T gpuComputeGyration(T * ag, T *xcm, T *x, T *tmp, T *mass, T *box, T masstotal,
        int *ilist, int *type, int *image, int triclinic, int dim, int inum);
template <typename T> void gpuComputeAngmom(T *lmom, T *p, T *xcm, T *x, T *v, T *tmp, T *mass, T *box, 
        int *ilist, int *type, int *image, int triclinic, int dim, int inum);
template <typename T> void gpuComputeTorque(T *tq, T *q, T *xcm, T *x, T *f, T *tmp, T *box, 
        int *ilist, int *image, int triclinic, int dim, int inum);
template <typename T> void gpuComputeInertia(T *inertia, T *ione, T *xcm, T *x, T *tmp, T *mass, T *box, 
        int *ilist, int *type, int *image, int triclinic, int dim, int inum);
template <typename T> void gpuComputeKEAtom(T *ke, T *mass, T *v, 
        T mvv2e, int *type, int *ilist, int dim, int inum);
template <typename T> void gpuComputeStressAtom(T *stress, T *mass, T *vatom, T *v, T mvv2e, 
        T nktv2p, int *type, int *ilist,  int vflag, int keflag, int dim, int inum);
template <typename T> void gpuComputeCentroidStressAtom(T *stress, T *mass, T *cvatom, T *v, T mvv2e, 
        T nktv2p, int *type, int *ilist,  int vflag, int keflag, int dim, int inum);
template <typename T> void gpuComputeDisplaceAtom(T *displace, T *x, T *xoriginal, T *box, 
        int *pbc, int *ilist,  int triclinic, int dim, int inum);
template <typename T> void gpuComputeTempSymTensor(T *ke_tensor, T *stress, T *v, T *tmp, T *mass, T tfactor, 
        int *type, int *ilist, int dim, int inum);
template <typename T> T gpuComputeTempScalar(T *ke, T *v, T *tmp, T *mass, T tfactor, 
        int *type, int *ilist, int dim, int inum);
template <typename T> T gpuComputePressureScalar(T *virial, T volume, T temp, T tempdof, 
        T boltz, T nktv2p, int dim);
template <typename T> void gpuComputePressureSymTensor(T *vector, T *virial, T *ke_tensor, 
        T volume, T nktv2p, int dim);
template <typename T> void gpuComputeHeatFlux(T *vector, T *jc, T *ke, T *pe, T *stress, T *v, 
        T *tmp, T nktv2p, int *ilist,  int pressatomflag, int dim, int inum);
template <typename T> void gpuComputeOrientOrderAtom(T *qnarray, T *x, T *rlist, T *cglist, T *fac, 
        T *qnm_r, T *qnm_i, T *distsq, T cutsq, T MY_EPSILON, T QEPSILON, T MY_4PI, int *neighlist, 
        int *neighnum,  int *ilist, int *qlist, int *nearest,  int nqlist, int qmax, int wlflag,  
        int wlhatflag,  int qlcompflag, int iqlcomp, int qlcomp, int nnn, int jnum, int dim, int inum) ;
template <typename T> void gpuComputeCoordAtomCutoff(int *cvec, T *x, T *rcutsq, int *type, 
        int *ilist, int *neighlist, int *neighnum, int *typelo, int *typehi, 
         int *jgroupbit, int dim, int ntypes, int jnum, int inum);
template <typename T> void gpuComputeCoordAtomCutoff(int *carray, T *x, T *rcutsq, int *type, 
        int *ilist, int *neighlist, int *neighnum, int *typelo, int *typehi, 
         int *jgroupbit, int ncol, int dim, int ntypes, int jnum, int inum);
template <typename T> void gpuComputeCoordAtomOrient(int *cvec, T *x, T *rcutsq, T *normv, T threshold, 
        int *type, int *ilist, int *neighlist, int *neighnum, int *typelo, int *typehi, 
         int *jgroupbit, int dim, int ntypes, int nqlist, int ncol, int l, int jnum, int inum);
template <typename T> void gpuComputeMSD(T *msd, T *vec, T *x, T *xoriginal, T *box, T *xcm, T *tmp,
         int *ilist, int *image, int naverage, int avflag, int triclinic, int nmsd, int dim,  int inum);
template <typename T> void gpuComputeVACF(T *vacf, T *vec, T *v, T *voriginal, T *tmp,
         int *ilist, int nvacf, int dim,  int inum);

//********************************* Velocity integration ****************************************//
template <typename T> void gpuVelocity(T *x, T *v, T *f, T *box, T *xcm, T *vcm, 
        T *mass, T *second, T *omega, T *vext, T *v_lo, T *v_hi, T *coord_lo, T *coord_hi, 
        T t_desired, T t_current, int *seed, int *save, int *map, int *image, int *type, 
        int *coord_dim, int *vdim, int sum_flag, int dist_flag, int loop_flag, 
        int rotation_flag, int momentum_flag, int triclinic, int dim, int mpiRank, 
        int vmode, int nlocal, int natoms);

template <typename T> void gpuSetVelocityInitialIntegrate(T *x, T *v, T *f, T *mass, T *fparam,
        T dtf, T dtv, int *type, int *ilist, int *iparam, int dim, int inum);

template <typename T> void gpuSetVelocityFinalIntegrate(T *x, T *v, T *f, T *mass, 
        T dtf, int *type, int *ilist, int *iparam, int dim, int inum);

template <typename T> void gpuInitialIntegrate(T *x, T *v, T *f, T *mass, T *dtarray, T *tarray,
        T *eta_mass, T *eta, T *eta_dot, T *eta_dotdot, T vlimitsq, int *type, int *ilist,
        int eta_mass_flag, int biasflag, int mtchain, int nc_tchain, int mode, int dim, int inum);

template <typename T> void gpuFinalIntegrate(T *x, T *v, T *f, T *mass, T *dtarray, T *tarray,
        T *eta_mass, T *eta, T *eta_dot, T *eta_dotdot, T *ke, T *tmp, T vlimitsq, int *type, int *ilist, 
        int eta_mass_flag, int biasflag, int mtchain, int nc_tchain, int mode, int dim, int inum);

template <typename T> T gpuVelocityRescalingThermostat(T *v, T *mass, T *dtarray, T *tarray, 
        T *ke, T *tmp, T *second, T energy, int *type, int *ilist, int *seed, int *save, 
        int biasflag, int mode, int dim, int inum);

template <typename T> void gpuPBC(T *x, T *v, int *image, T *boxhi, T *boxlo, T *hi_lambda, T *lo_lambda,  
        T *h, T *h_inv, T *h_rate, int *pbc, int vdeform, int triclinic, int dim, int nlocal);

//********************************* Force Constraints ****************************************//
template <typename T> void gpuFixSetForce(T *x, T *v, T *f, T *eatom, T *vatom, T *fparam, 
        int *iparam, int *ilist, int eflag_atom, int vflag_atom, int dim, int inum);
template <typename T> void gpuFixLineForce(T *x, T *v, T *f, T *eatom, T *vatom, T *fparam, 
        int *iparam, int *ilist, int eflag_atom, int vflag_atom, int dim, int inum);
template <typename T> void gpuFixPlaneForce(T *x, T *v, T *f, T *eatom, T *vatom, T *fparam, 
        int *iparam, int *ilist, int eflag_atom, int vflag_atom, int dim, int inum);
template <typename T> void gpuFixAveForce(T *x, T *v, T *f, T *eatom, T *vatom, T *fparam, 
        int *iparam, int *ilist, int eflag_atom, int vflag_atom, int dim, int inum);

template <typename T> void gpuFixAddForce(T *x, T *v, T *f, T *eatom, T *vatom, T *fparam, T *box,
        int *iparam, int *ilist, int *image, int triclinic, int eflag_atom, int vflag_atom, int dim, int inum);
template <typename T> void gpuFixDragForce(T *x, T *v, T *f, T *eatom, T *vatom, T *fparam, T *box,
        int *iparam, int *ilist, int *pbc, int triclinic, int eflag_atom, int vflag_atom, int dim, int inum);

template <typename T> void gpuFixWallReflect(T *x, T *v, T *f, T *eatom, T *vatom, T *fparam, 
        int *iparam, int *ilist, int eflag_atom, int vflag_atom, int dim, int inum);
template <typename T> void gpuFixWallHarmonic(T *x, T *v, T *f, T *eatom, T *vatom, T *fparam, 
        int *iparam, int *ilist, int eflag_atom, int vflag_atom, int dim, int inum);
template <typename T> void gpuFixWallLJ93(T *x, T *v, T *f, T *eatom, T *vatom, T *fparam, 
        int *iparam, int *ilist, int eflag_atom, int vflag_atom, int dim, int inum);
template <typename T> void gpuFixWallLJ126(T *x, T *v, T *f, T *eatom, T *vatom, T *fparam, 
        int *iparam, int *ilist, int eflag_atom, int vflag_atom, int dim, int inum);
template <typename T> void gpuFixWallLJ1043(T *x, T *v, T *f, T *eatom, T *vatom, T *fparam, 
        int *iparam, int *ilist, int eflag_atom, int vflag_atom, int dim, int inum);
template <typename T> void gpuFixWallMorse(T *x, T *v, T *f, T *eatom, T *vatom, T *fparam, 
        int *iparam, int *ilist, int eflag_atom, int vflag_atom, int dim, int inum);

//***************************  Empirical Potentials ***************************************//
template <typename T> void gpuSinglea(T *u, T *xi, T *qi, int *ti, int *ai, T *mu, T *eta, int *kappa, int dim, int ncq, int nmu, int neta, int nkappa, int ng, int potnum);
template <typename T> void gpuSingleb(T *u, T *xi, T *qi, int *ti, int *ai, T *mu, T *eta, int *kappa, int dim, int ncq, int nmu, int neta, int nkappa, int ng, int potnum);
template <typename T> void gpuPaira(T *u, T *xij, T *qi, T *qj, int *ti, int *tj, int *ai, int *aj, T *mu, T *eta, int *kappa, int dim, int ncq, int nmu, int neta, int nkappa, int ng, int potnum);
template <typename T> void gpuPairb(T *u, T *xij, T *qi, T *qj, int *ti, int *tj, int *ai, int *aj, T *mu, T *eta, int *kappa, int dim, int ncq, int nmu, int neta, int nkappa, int ng, int potnum);
template <typename T> void gpuPairc(T *u, T *xij, T *qi, T *qj, int *ti, int *tj, int *ai, int *aj, T *mu, T *eta, int *kappa, int dim, int ncq, int nmu, int neta, int nkappa, int ng, int potnum);
template <typename T> void gpuPaircDensity(T *u, T *rho, T *mu, T *eta, int *kappa, int nrho, int nmu, int neta, int nkappa, int ng, int potnum);
template <typename T> void gpuTripleta(T *u, T *xij, T *xik, T *qi, T *qj, T *qk, int *ti, int *tj, int *tk, int *ai, int *aj, int *ak, T *mu, T *eta, int *kappa, int dim, int ncq, int nmu, int neta, int nkappa, int ng, int potnum);
template <typename T> void gpuTripletb(T *u, T *xij, T *xik, T *qi, T *qj, T *qk, int *ti, int *tj, int *tk, int *ai, int *aj, int *ak, T *mu, T *eta, int *kappa, int dim, int ncq, int nmu, int neta, int nkappa, int ng, int potnum);
template <typename T> void gpuTripletc(T *u, T *xij, T *xik, T *qi, T *qj, T *qk, int *ti, int *tj, int *tk, int *ai, int *aj, int *ak, T *mu, T *eta, int *kappa, int dim, int ncq, int nmu, int neta, int nkappa, int ng, int potnum);
template <typename T> void gpuTripletcPair(T *u, T *xij, T *qi, T *qj, int *ti, int *tj, int *ai, int *aj, T *mu, T *eta, int *kappa, int dim, int ncq, int nmu, int neta, int nkappa, int ng, int potnum);
template <typename T> void gpuTripletcDensity(T *u, T *rho, T *mu, T *eta, int *kappa, int nrho, int nmu, int neta, int nkappa, int ng, int potnum);
template <typename T> void gpuQuadrupleta(T *u, T *xij, T *xik, T *xil, T *qi, T *qj, T *qk, T *ql, int *ti, int *tj, int *tk, int *tl, int *ai, int *aj, int *ak, int *al, T *mu, T *eta, int *kappa, int dim, int ncq, int nmu, int neta, int nkappa, int ng, int potnum);
template <typename T> void gpuQuadrupletb(T *u, T *xij, T *xik, T *xil, T *qi, T *qj, T *qk, T *ql, int *ti, int *tj, int *tk, int *tl, int *ai, int *aj, int *ak, int *al, T *mu, T *eta, int *kappa, int dim, int ncq, int nmu, int neta, int nkappa, int ng, int potnum);
// template <typename T> void gpuLJ(T *u, T *f, T *xij, T *qi, T *qj, int *ti, int *tj, int *ai, int *aj, T *mu, T *eta, int *kappa, int dim, int ncq, int nmu, int neta, int nkappa, int ng);
// template <typename T> void gpuGradientLJ(T *u, T *du, T *xij, T *u_x, T *qi, T *qj, int *ti, int *tj, 
//         int *ai, int *aj, T *mu, T *eta, int *kappa, int dim, int ncq, 
//                 int nmu, int neta, int nkappa, int ng);

template <typename T> void gpuSingleaGradient(T *u, T *u_xi, T *xi, T *qi, int *ti, int *ai, T *mu, T *eta, int *kappa, int dim, int ncq, int nmu, int neta, int nkappa, int ng, int potnum);
template <typename T> void gpuSinglebGradient(T *u, T *u_xi, T *xi, T *qi, int *ti, int *ai, T *mu, T *eta, int *kappa, int dim, int ncq, int nmu, int neta, int nkappa, int ng, int potnum);
template <typename T> void gpuPairaGradient(T *u, T *u_xij, T *xij, T *qi, T *qj, int *ti, int *tj, int *ai, int *aj, T *mu, T *eta, int *kappa, int dim, int ncq, int nmu, int neta, int nkappa, int ng, int potnum);
template <typename T> void gpuPairbGradient(T *u, T *u_xij, T *xij, T *qi, T *qj, int *ti, int *tj, int *ai, int *aj, T *mu, T *eta, int *kappa, int dim, int ncq, int nmu, int neta, int nkappa, int ng, int potnum);
template <typename T> void gpuPaircGradient(T *u, T *u_xij, T *xij, T *qi, T *qj, int *ti, int *tj, int *ai, int *aj, T *mu, T *eta, int *kappa, int dim, int ncq, int nmu, int neta, int nkappa, int ng, int potnum);
template <typename T> void gpuPaircDensityGradient(T *u, T *u_rho, T *rho, T *mu, T *eta, int *kappa, int nrho, int nmu, int neta, int nkappa, int ng, int potnum);
template <typename T> void gpuTripletaGradient(T *u, T *u_xij, T *u_xik, T *xij, T *xik, T *qi, T *qj, T *qk, int *ti, int *tj, int *tk, int *ai, int *aj, int *ak, T *mu, T *eta, int *kappa, int dim, int ncq, int nmu, int neta, int nkappa, int ng, int potnum);
template <typename T> void gpuTripletbGradient(T *u, T *u_xij, T *u_xik, T *xij, T *xik, T *qi, T *qj, T *qk, int *ti, int *tj, int *tk, int *ai, int *aj, int *ak, T *mu, T *eta, int *kappa, int dim, int ncq, int nmu, int neta, int nkappa, int ng, int potnum);
template <typename T> void gpuTripletcGradient(T *u, T *u_xij, T *u_xik, T *xij, T *xik, T *qi, T *qj, T *qk, int *ti, int *tj, int *tk, int *ai, int *aj, int *ak, T *mu, T *eta, int *kappa, int dim, int ncq, int nmu, int neta, int nkappa, int ng, int potnum);
template <typename T> void gpuTripletcPairGradient(T *u, T *u_xij, T *xij, T *qi, T *qj, int *ti, int *tj, int *ai, int *aj, T *mu, T *eta, int *kappa, int dim, int ncq, int nmu, int neta, int nkappa, int ng, int potnum);
template <typename T> void gpuTripletcDensityGradient(T *u, T *u_rho, T *rho, T *mu, T *eta, int *kappa, int nrho, int nmu, int neta, int nkappa, int ng, int potnum);
template <typename T> void gpuQuadrupletaGradient(T *u, T *u_xij, T *u_xik, T *u_xil, T *xij, T *xik, T *xil, T *qi, T *qj, T *qk, T *ql, int *ti, int *tj, int *tk, int *tl, int *ai, int *aj, int *ak, int *al, T *mu, T *eta, int *kappa, int dim, int ncq, int nmu, int neta, int nkappa, int ng, int potnum);
template <typename T> void gpuQuadrupletbGradient(T *u, T *u_xij, T *u_xik, T *u_xil, T *xij, T *xik, T *xil, T *qi, T *qj, T *qk, T *ql, int *ti, int *tj, int *tk, int *tl, int *ai, int *aj, int *ak, int *al, T *mu, T *eta, int *kappa, int dim, int ncq, int nmu, int neta, int nkappa, int ng, int potnum);

template <typename T> void gpuSingleaGradient(T *u, T *du, T *u_xi, T *xi, T *qi, int *ti, int *ai, T *mu, T *eta, int *kappa, int dim, int ncq, int nmu, int neta, int nkappa, int ng, int potnum);
template <typename T> void gpuSinglebGradient(T *u, T *du, T *u_xi, T *xi, T *qi, int *ti, int *ai, T *mu, T *eta, int *kappa, int dim, int ncq, int nmu, int neta, int nkappa, int ng, int potnum);
template <typename T> void gpuPairaGradient(T *u, T *du, T *u_xij, T *xij, T *qi, T *qj, int *ti, int *tj, int *ai, int *aj, T *mu, T *eta, int *kappa, int dim, int ncq, int nmu, int neta, int nkappa, int ng, int potnum);
template <typename T> void gpuPairbGradient(T *u, T *du, T *u_xij, T *xij, T *qi, T *qj, int *ti, int *tj, int *ai, int *aj, T *mu, T *eta, int *kappa, int dim, int ncq, int nmu, int neta, int nkappa, int ng, int potnum);
template <typename T> void gpuPaircGradient(T *u, T *du, T *u_xij, T *xij, T *qi, T *qj, int *ti, int *tj, int *ai, int *aj, T *mu, T *eta, int *kappa, int dim, int ncq, int nmu, int neta, int nkappa, int ng, int potnum);
template <typename T> void gpuPaircDensityGradient(T *u, T *du, T *u_rho, T *rho, T *mu, T *eta, int *kappa, int nrho, int nmu, int neta, int nkappa, int ng, int potnum);
template <typename T> void gpuTripletaGradient(T *u, T *du, T *u_xij, T *u_xik, T *xij, T *xik, T *qi, T *qj, T *qk, int *ti, int *tj, int *tk, int *ai, int *aj, int *ak, T *mu, T *eta, int *kappa, int dim, int ncq, int nmu, int neta, int nkappa, int ng, int potnum);
template <typename T> void gpuTripletbGradient(T *u, T *du, T *u_xij, T *u_xik, T *xij, T *xik, T *qi, T *qj, T *qk, int *ti, int *tj, int *tk, int *ai, int *aj, int *ak, T *mu, T *eta, int *kappa, int dim, int ncq, int nmu, int neta, int nkappa, int ng, int potnum);
template <typename T> void gpuTripletcGradient(T *u, T *du, T *u_xij, T *u_xik, T *xij, T *xik, T *qi, T *qj, T *qk, int *ti, int *tj, int *tk, int *ai, int *aj, int *ak, T *mu, T *eta, int *kappa, int dim, int ncq, int nmu, int neta, int nkappa, int ng, int potnum);
template <typename T> void gpuTripletcPairGradient(T *u, T *du, T *u_xij, T *xij, T *qi, T *qj, int *ti, int *tj, int *ai, int *aj, T *mu, T *eta, int *kappa, int dim, int ncq, int nmu, int neta, int nkappa, int ng, int potnum);
template <typename T> void gpuTripletcDensityGradient(T *u, T *du, T *u_rho, T *rho, T *mu, T *eta, int *kappa, int nrho, int nmu, int neta, int nkappa, int ng, int potnum);
template <typename T> void gpuQuadrupletaGradient(T *u, T *du, T *u_xij, T *u_xik, T *u_xil, T *xij, T *xik, T *xil, T *qi, T *qj, T *qk, T *ql, int *ti, int *tj, int *tk, int *tl, int *ai, int *aj, int *ak, int *al, T *mu, T *eta, int *kappa, int dim, int ncq, int nmu, int neta, int nkappa, int ng, int potnum);
template <typename T> void gpuQuadrupletbGradient(T *u, T *du, T *u_xij, T *u_xik, T *u_xil, T *xij, T *xik, T *xil, T *qi, T *qj, T *qk, T *ql, int *ti, int *tj, int *tk, int *tl, int *ai, int *aj, int *ak, int *al, T *mu, T *eta, int *kappa, int dim, int ncq, int nmu, int neta, int nkappa, int ng, int potnum);
template <typename T> void gpuElectronDensity(T *rhoi, T *rhoij, int *pairnum, int *pairnumsum, int inum);
template <typename T> void gpuEmbedingForce(T *fij, T *d_rhoi, int *pairnum, int *pairnumsum, int inum);

//***************************  spherical harmonic potentials ***************************************//
template <typename T> void gpuSphericalHarmonicsBessel(T *Sr, T *Si, T *x, T *y, T *z, 
                T *x0, T *P, T *tmp, T *f, T *fac, T pi, int L, int K, int N);
template <typename T> void gpuSphericalHarmonicsBessel(T *Sr, T *Si, T *xij,
                T *x0, T *P, T *tmp, T *f, T *fac, T pi, int L, int K, int N);
template <typename T> void gpuSphericalHarmonicsBesselWithDeriv(T *Sr, T *Si, T *Srx, T *Six, T *Sry, T *Siy, T *Srz, T *Siz, T *x, T *y, T *z, 
                T *x0, T *P, T *tmp, T *f, T *dP, T *dtmp, T *df, T *fac, T pi, int L, int K, int N);
template <typename T> void gpuSphericalHarmonicsBesselWithDeriv(T *Sr, T *Si, T *Srx, T *Six, T *Sry, T *Siy, T *Srz, T *Siz, T *xij, 
                T *x0, T *P, T *tmp, T *f, T *dP, T *dtmp, T *df, T *fac, T pi, int L, int K, int N);
template <typename T> void gpuRadialSphericalHarmonicsSpectrum(T *c, T *ar, T *ai, T *sr, T *si, T *cg, int *indk, 
        int *indl, int *indm, int *rowm, int *Nnb, int Nub, int Ncg, int Na, int L, int K, int spectrum);
template <typename T> void gpuRadialSphericalHarmonicsSpectrumDeriv(T *cd, T *ar, T *ai, 
        T *srx, T *six, T *sry, T *siy, T *srz, T *siz, T *cg, int *indk, int *indl, int *indm, 
        int *rowm, int *Nnb, int Nub, int Ncg, int Na, int L, int K, int spectrum);
template <typename T> void gpuRadialSphericalHarmonicsSpectrumDeriv(T *cd, T *ar, T *ai, 
        T *Srx, T *Six, T *Sry, T *Siy, T *Srz, T *Siz, T *cg, int *indk, int *indl, int *indm, 
        int *rowm, int *idxi, int *Nnb, int Nub, int Ncg, int Na, int Nij, int L, int K, int spectrum);

template <typename T> void gpuForceDecomposition(T *f, T *fij, int *ai, int *aj, 
        int inum, int ijnum, int Nbf);
template <typename T> void gpuForceDecomposition(T *f, T *fij, int *ai, int *aj, int ijnum);
template <typename T> void gpuCenterAtomDecomposition(T *f, T *fij, int *ilist, int *anumsum, 
        int inum, int ijnum, int Na, int Nbf);
template <typename T> void gpuNeighborAtomDecomposition(T *f, T *fij, int *jlist, int *bnumsum, int *index, 
        int inum, int ijnum, int jnum, int Nbf);
template <typename T> void gpuCenterAtomDecomposition(T *f, T *fij, int *ilist, int *anumsum, int inum);
template <typename T> void gpuNeighborAtomDecomposition(T *f, T *fij, int *jlist, int *bnumsum, int *index, int jnum);
template <typename T> void gpuCenterAtomDecomposition(T *e, T *ei, int *ilist, int Na);

//********************************** SNAP POTENTIAL ***********************************************//
void gpuBuildIndexList(int *idx_max, int *idxz, int *idxz_block, int *idxb, int *idxb_block, int *idxu_block, 
        int *idxcg_block, int twojmax);

template <typename T> void gpuInitRootpqArray(T *rootpqarray, int twojmax);
template <typename T> void gpuInitClebschGordan(T *cglist, T *factorial, int twojmax);

template <typename T> void gpuInitSna(T *rootpqarray, T *cglist, T *factorial, int *idx_max, int *idxz, 
     int *idxz_block, int *idxb, int *idxb_block, int *idxu_block, int *idxcg_block, int twojmax);

template <typename T> void gpuZeroUarraytot(T *ulisttot_r, T *ulisttot_i, T wself, int *idxu_block, int *type,
        int *map, int *ai, int wselfall_flag, int chemflag, int idxu_max, int nelements, int twojmax, int inum);

template <typename T> void gpuAddUarraytot(T *ulisttot_r, T *ulisttot_i, T *ulist_r, T *ulist_i, T *rij, 
        T *wjelem, T *radelem, T rmin0, T rcutfac, int *idxu_block, int *ilist, 
        int *type, int *pairnum, int *pairnumsum, int *map, int *tj, int twojmax, int idxu_max, 
        int nelements, int inum, int switch_flag, int chemflag);

template <typename T> void gpuComputeBeta2(T *beta, T *bispectrum, T *coeffelem, int *ilist, int *map, int *type, 
        int inum, int ncoeff, int ncoeffall, int quadraticflag);

template <typename T> void gpuComputeUij(T *ulist_r, T *ulist_i, T *rootpqarray, T *rij, 
        T *radelem, T rmin0, T rfac0, T rcutfac, int *idxu_block, 
        int *ai, int *aj, int *ti, int *tj, int twojmax, int idxu_max, int ijnum);

template <typename T> void gpuComputeZi(T *zlist_r, T *zlist_i, T *ulisttot_r, T *ulisttot_i, T *cglist,
        int *idxz, int *idxu_block, int *idxcg_block, int twojmax, int idxu_max, int idxz_max, int nelements, 
        int bnorm_flag, int inum);

template <typename T> void gpuComputeYi(T *ylist_r, T *ylist_i, T *ulisttot_r, T *ulisttot_i, T *cglist, 
        T* betalist, int *idxz, int *idxb_block, int *idxu_block, int *idxcg_block, int twojmax, 
        int idxb_max, int idxu_max, int idxz_max, int nelements, int bnorm_flag, int ncoeff, int inum);

template <typename T> void gpuComputeBi(T *blist, T *zlist_r, T *zlist_i, T *ulisttot_r, T *ulisttot_i, 
        T *bzero,int *ilist, int *type, int *map, int *idxb, int *idxu_block, int *idxz_block, 
        int twojmax, int idxb_max, int idxu_max, int idxz_max, int nelements, int bzero_flag, 
        int wselfall_flag, int chemflag, int inum);

template <typename T> void gpuComputeDbidrj(T *dblist, T *zlist_r, T *zlist_i, 
        T *dulist_r, T *dulist_i, int *idxb, int *idxu_block, int *idxz_block, 
        int *map, int *ai, int *aj, int *ti, int *tj, int twojmax, int idxb_max, int idxu_max, int idxz_max, 
        int nelements, int bnorm_flag, int chemflag, int ijnum);

template <typename T> void gpuComputeDbidrj(T *dblist, T *zlist_r, T *zlist_i, 
        T *dulist_r, T *dulist_i, int *idxb, int *idxu_block, int *idxz_block, 
        int *map, int *ai, int *tj, int twojmax, int idxb_max, int idxu_max, int idxz_max, 
        int nelements, int bnorm_flag, int chemflag, int inum, int ijnum);

template <typename T> void gpuComputeDuijdrj(T *dulist_r, T *dulist_i, T *ulist_r, T *ulist_i, 
    T *rootpqarray, T* rij, T *wjelem, T *radelem, T rmin0, T rfac0, T rcutfac, int *idxu_block,
    int *ai, int *aj, int *ti, int *tj, int twojmax, int idxu_max, int ijnum, int switch_flag);

template <typename T> void gpuComputeDeidrj(T *dedr, T *ylist_r, T *ylist_i, 
        T *dulist_r, T *dulist_i, int *idxu_block, int *map, int *ai, int *aj, int *ti, int *tj,
        int nelements, int twojmax, int idxu_max, int chemflag, int ijnum); 

// template <typename T> void gpuComputeDbidrj(T *dblist, T *zlist_r, T *zlist_i, T *dulist_r, T *dulist_i, 
//         int *idxb, int *idxu_block, int *idxz_block, int *type, int *map, int *ai, int *aj, int twojmax, 
//         int idxb_max, int idxu_max, int idxz_max, int nelements, int bnorm_flag, int chemflag, int ijnum);

template <typename T> void gpuComputeSna(T *sna, T *blist, int *ilist, int *mask, 
        int ncoeff, int nrows, int inum, int quadraticflag);

template <typename T> void gpuComputeSna(T *sna, T *blist, int *ilist, int *mask, int *type,
        int ncoeff, int ntype, int nperdim, int inum, int quadraticflag);

template <typename T> void gpuComputeSnad(T *snad, T *dblist, T *blist, int *aii, int *ai, int *aj, int *ti, 
        int *mask, int ncoeff, int ntypes, int nperdim, int ijnum, int quadraticflag);

template <typename T> void gpuComputeSnad(T *snad, T *dblist, T *blist, int *aii, int *ai, int *aj, int *ti, 
        int *mask, int *tag, int ncoeff, int ntypes, int nperdim, int ijnum, int quadraticflag);

template <typename T> void gpuComputeSnav(T *snav, T *dblist, T *blist, T *x, int *aii, int *ai, int *aj,
        int *ti, int *mask, int ncoeff, int ntypes, int nperdim, int ijnum, int quadraticflag);

template <typename T> void gpuComputeSnav2(T *snav, T *dblist, T *blist, T *x, int *aii, int *ai, int *aj,
        int *ti, int *mask, int ncoeff, int ntypes, int nperdim, int ijnum, int quadraticflag);

template <typename T> void gpuSnapTallyEnergyFull(T *eatom, T *bispectrum, T *coeffelem, int *ilist, int *map, int *type, 
        int inum, int ncoeff, int ncoeffall, int quadraticflag);

template <typename T> void gpuSnapTallyForceFull(T *fatom, T *fij, int *ai, 
        int *aj, int *alist, int ijnum);

template <typename T> void gpuSnapTallyVirialFull(T *vatom, T *fij, T *rij, int *ai, int *aj, int inum, int ijnum);

template <typename T> void gpuNeighPairList(int *pairnum, int *pairlist, T *x, T rcutsq, int *ilist, int *neighlist, 
        int *neighnum, int inum, int jnum, int dim);

template <typename T> void gpuNeighPairList(int *pairnum, int *pairlist, T *x, T *rcutsq, int *ilist, int *neighlist, 
        int *neighnum, int *atomtype, int inum, int jnum, int dim, int ntypes);

// template <typename T> void gpuNeighPairs(T *xij, T *x, int *aii, int *ai, int *aj,  
//       int *ti, int *tj, int *pairnum, int *pairlist, int *pairnumsum, int *ilist, 
//       int *atomtype, int inum, int jnum, int dim);

template <typename T> void gpuNeighPairs(T *xij, T *x, int *aii, int *ai, int *aj,  
      int *ti, int *tj, int *pairnum, int *pairlist, int *pairnumsum, int *ilist, 
      int *atomtype, int *alist, int inum, int jnum, int dim);

// GPU-optimized snap implementation
template <typename T> void gpuComputeSij(T *Sr, T *Si, T *Srx, T *Six, T *Sry, T *Siy, T *Srz, T *Siz, 
        T *rootpqarray, T *rij, T *wjelem, T *radelem, T rmin0, T rfac0, T rcutfac, int *idxu_block,  
        int *ti, int *tj, int twojmax, int idxu_max, int ijnum, int switch_flag);                

template <typename T> void gpuZeroUarraytot2(T *Stotr, T *Stoti, T wself, int *idxu_block, 
        int *type, int *map, int *ai, int wselfall_flag, int chemflag, int idxu_max, int nelements, 
         int twojmax, int inum);

template <typename T> void gpuAddUarraytot(T *Stotr, T *Stoti, T *Sr, 
        T *Si, int *pairnum, int *pairnumsum, int *map, int *tj, 
        int idxu_max, int nelements, int inum, int ijnum, int chemflag);

template <typename T> void gpuAddUarraytot(T *Stotr, T *Stoti, T *Sr, 
        T *Si, int *map, int *ai, int *tj, int idxu_max, int inum, int ijnum, int chemflag);

template <typename T> void gpuComputeZi2(T *zlist_r, T *zlist_i, T *Stotr, T *Stoti, 
        T *cglist, int *idxz, int *idxu_block, int *idxcg_block, int twojmax, int idxu_max, 
        int idxz_max, int nelements, int bnorm_flag, int inum);

template <typename T> void gpuComputeBi2(T *blist, T *zlist_r, T *zlist_i, T *Stotr, T *Stoti, 
        T *bzero, int *ilist, int *type, int *map, int *idxb, int *idxu_block, int *idxz_block, 
        int twojmax, int idxb_max, int idxu_max, int idxz_max, int nelements, int bzero_flag, 
        int wselfall_flag, int chemflag, int inum);

template <typename T> void gpuComputeBeta(T *beta, T *bispectrum, T *coeffelem, int *ilist, int *map, int *type, 
        int inum, int ncoeff, int ncoeffall, int quadraticflag);

template <typename T> void gpuComputeYi(T *ylist_r, T *ylist_i, T *Stotr, T *Stoti, T *cglist, 
        T* beta, int *idxz, int *idxb_block, int *idxu_block, int *idxcg_block, int twojmax, 
        int idxb_max, int idxu_max, int idxz_max, int nelements, int bnorm_flag, int inum);

template <typename T> void gpuComputeYi(T *ylist_r, T *ylist_i, T *zlist_r, T *zlist_i, 
        T* beta, int *idxz, int *idxb_block, int twojmax, int idxb_max, int idxu_max, int idxz_max,
        int nelements, int bnorm_flag, int inum);

template <typename T> void gpuComputeDeidrj(T *dedr, T *ylist_r, T *ylist_i, 
        T *dulist_r, T *dulist_i, int *idxu_block, int *map, int *ai, int *tj,
        int twojmax, int idxu_max, int chemflag, int inum, int ijnum); 

template <typename T> void gpuSnapTallyEnergyFull2(T *eatom, T *bispectrum, T *coeffelem, int *ilist, 
        int *map, int *type, int inum, int ncoeff, int ncoeffall, int quadraticflag);

template <typename T> void gpuSnapTallyBispectrum(T *bi, T *bispectrum, int *ilist, 
        int *type, int inum, int ncoeff, int nperdim, int ntype, int quadraticflag);

template <typename T> void gpuSnapTallyBispectrumDeriv(T *db, T *bispectrum, T *dbdr, int *aii, 
        int *ai, int *aj, int *ti, int inum, int ijnum, int ncoeff, int nperdim, int ntype, int quadraticflag);

template <typename T> void gpuSnapTallyBispectrumVirial(T *bv, T *bispectrum, T *dbdr, T *rij, int *aii, 
        int *ai, int *aj, int *ti, int inum, int ijnum, int ncoeff, int nperdim, int ntype, int quadraticflag);

#endif  

