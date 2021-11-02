/***************************************************************************                               
                    Molecular Dynamics Potentials (MDP)
                           CESMIX-MIT Project  
 
 Contributing authors: Ngoc-Cuong Nguyen (cuongng@mit.edu, exapde@gmail.com)
 ***************************************************************************/

#ifndef __MDPCPUCORE_H__
#define __MDPCPUCORE_H__

template <typename T> void cpuKron(T *C, T *A, T *B, int M1, int M2);
template <typename T> void cpuKron(T *C, T *A, T *B, int M1, int N1, int M2, int N2);
void cpuIndexPermute12(int *index, int I1, int I2, int I3);
void cpuIndexPermute13(int *index, int I1, int I2, int I3, int I4);
void cpuIndexPermute23(int *index, int I1, int I2, int I3, int I4);
template <typename T> void cpuPermute(T *B, T *A, int *index, int N);
template <typename T> void cpuPermuteSharedMem(T *B, T *A, int *index, int N, int BLOCKDIM);
template <typename T> void cpuPermute12(T *B, T *A, int I1, int I2, int I3);
template <typename T> void cpuPermute13(T *B, T *A, int I1, int I2, int I3, int I4);
template <typename T> void cpuPermute23(T *B, T *A, int I1, int I2, int I3, int I4);

template <typename T> void cpuGetArrayAtIndex(T *y, T *x, int *ind, int n);
template <typename T> void cpuPutArrayAtIndex(T *y, T *x, int *ind, int n);
template <typename T> void cpuArrayPlusXAtIndex(T *y, T *x, int *ind, int n);
template <typename T> void cpuArrayMinusXAtIndex(T *y, T *x, int *ind, int n);
template <typename T> void cpuArrayAXPYAtIndex(T *y, T *x, T a, int *ind, int n);

template <typename T> void cpuArrayTranspose(T *A, T *B, int m, int n);
template <typename T> void cpuArrayPlusAtColumnIndex(T *A, T *B, int *colind, int m, int n);
template <typename T> void cpuArrayMinusAtColumnIndex(T *A, T *B, int *colind, int m, int n);
template <typename T> void cpuGetArrayAtColumnIndex(T *A, T *B, int *colind, int m, int n);
template <typename T> void cpuArrayTransposeAtColumnIndex(T *A, T *B, int *colind, int m, int n);
template <typename T> void cpuGetArrayAtRowIndex(T *A, T *B, int *rowind, int m, int n);
template <typename T> void cpuArrayTransposeAtRowIndex(T *A, T *B, int *rowind, int m, int n);
template <typename T> void cpuArrayRowSum(T *y, T *x, int m, int n);
template <typename T> void cpuArrayRowSquareSum(T *y, T *x, int m, int n);
template <typename T> void cpuArrayDistSquareSum(T *y, T *x1, T *x2, int m, int n);

void cpuQuadrupletnum(int *a, int *b, int n);
void cpuTripletnum(int *a, int *b, int n);
void cpuArrayFill(int *a, int start, int n);
template <typename T> T cpuArrayMaxError(T *a, T *b, int n);
template <typename T> T cpuArrayMin(T *a, int n);
template <typename T> T cpuArrayMax(T *a, int n);
template <typename T> T cpuArraySum(T *a, int n);
template <typename T> void cpuArraySumEveryColumn(T *b, T *a, int m, int n);
template <typename T> void cpuArrayMaxEveryColumn(T *b, T *a, int m, int n);
template <typename T> void cpuArrayMinEveryColumn(T *b, T *a, int m, int n);
template <typename T> void cpuArraySumEveryRow(T *b, T *a, int m, int n);
template <typename T> void cpuArrayMaxEveryRow(T *b, T *a, int m, int n);
template <typename T> void cpuArrayMinEveryRow(T *b, T *a, int m, int n);
template <typename T> T cpuArraySquareSum(T *a, int n);
template <typename T> void cpuArrayDistSquareSum(T *y, T *x1, T *x2, int m, int n);
template <typename T> T cpuArrayMean(T *a, int n);
template <typename T> void cpuArraySetValue(T *y, T a, int n);
template <typename T> void cpuArrayAddScalar(T *y, T a, int n);
template <typename T> void cpuArrayMultiplyScalar(T *y, T a, int n);

template <typename T> void cpuArrayCopy(T *y, T *x, int n);
template <typename T> void cpuArrayMinus(T *y, T *x, int n);
template <typename T> void cpuArrayAbs(T *y, T *x, int n);
template <typename T> void cpuArraySqrt(T *y, T *x, int n);
template <typename T> void cpuArraySin(T *y, T *x, int n);
template <typename T> void cpuArrayCos(T *y, T *x, int n);
template <typename T> void cpuArrayTan(T *y, T *x, int n);
template <typename T> void cpuArrayAsin(T *y, T *x, int n);
template <typename T> void cpuArrayAcos(T *y, T *x, int n);
template <typename T> void cpuArrayAtan(T *y, T *x, int n);
template <typename T> void cpuArraySinh(T *y, T *x, int n);
template <typename T> void cpuArrayCosh(T *y, T *x, int n);
template <typename T> void cpuArrayTanh(T *y, T *x, int n);
template <typename T> void cpuArrayAsinh(T *y, T *x, int n);
template <typename T> void cpuArrayAcosh(T *y, T *x, int n);
template <typename T> void cpuArrayAtanh(T *y, T *x, int n);
template <typename T> void cpuArrayExp(T *y, T *x, int n);
template <typename T> void cpuArrayLog(T *y, T *x, int n);
template <typename T> void cpuArrayCeil(T *y, T *x, int n);
template <typename T> void cpuArrayFloor(T *y, T *x, int n);
template <typename T> void cpuArrayErf(T *y, T *x, int n);
template <typename T> void cpuArrayErfc(T *y, T *x, int n);
template <typename T> void cpuArraySquare(T *y, T *x, int n);
template <typename T> void cpuArrayPower(T *y, T *x, int p, int n);

template <typename T> void cpuArrayMultiplyScalarDiagonal(T *C, T a, int n);
template <typename T> void cpuArrayAddVectorToDiagonal(T *C, T *x, T a, int n);
template <typename T> void cpuArrayRowAverage(T *a, T *b, int m, int n);
template <typename T> void cpuArrayAXPB(T *y, T *x, T a, T b, int n);
template <typename T> void cpuArrayRowkAXPB(T *y, T *x, T a, T b, int m, int n, int k);
template <typename T> void cpuArrayAXPBY(T *z, T *x, T *y, T a, T b, int n);
template <typename T> void cpuArrayAXY(T *s, T *x, T *y, T a, int n);
template <typename T> void cpuArrayAXYZ(T *s, T *x, T *y, T *z, T a, int n);
template <typename T> void cpuArrayAXYPBZ(T *s, T *x, T *y, T *z, T a, T b, int n);
template <typename T> void cpuArrayAdd3Vectors(T *s, T *x, T *y, T *z, T a, T b, T c, int n);
template <typename T> void cpuArrayExtract(T *un, T *u, int I, int J, int K, 
        int i1, int i2, int j1, int j2, int k1, int k2);
template <typename T> void cpuArrayInsert(T *u, T *un, int I, int J, int K, 
        int i1, int i2, int j1, int j2, int k1, int k2);
template <typename T> void cpuArrayGemmBatch(T *C, T *A, T *B, int I, int J, int K, int S);
template <typename T> void cpuArrayGemmBatch1(T *C, T *A, T *B, int I, int J, int K, int S);
template <typename T> void cpuArrayDG2CG(T *ucg, T *udg, int *cgent2dgent, int *rowent2elem, int nent);
template <typename T> void cpuArrayDG2CG2(T *ucg, T *udg, int *colent2elem, int *rowent2elem, int nent, int npe);

template <typename T> void cpuCart2Sphere(T *the, T *phi, T *r, T *x, T *y, T *z, int N);
template <typename T> void cpuCart2SphereDeriv(T *the, T *phi, T *r, T *thex, T *they, T *thez, T *phix, T *phiy, T *phiz, T *rx, T *ry, T *rz, T *x, T *y, T *z, int N);
template <typename T> void cpuSphere2Cart(T *x, T *y, T *z, T *the, T *phi, T *r, int N);
template <typename T> void cpuEuler2Rotm(T *R11, T *R12, T *R13, T *R21, 
                T *R22, T *R23, T *R31, T *R32, T *R33, T *alpha, T *beta, T *gamma, int N);
template <typename T> void cpuRotc(T *X, T *Y, T *Z, T *R, T *x, T *y, T *z, int N);

void cpuCumsum(int* output, int* input, int length);
template <typename T> void cpuSmallMatrixInverse(T *invA, T *A, int dim);
int cpuUniqueSort(int *b, int *c, int *d, int *e, int *a, int *p, int *t, int *q, int n);
int cpuFindAtomType(int *tlist, int* ilist, int *atomtype, int *p, int *q, int typei, int inum);
void cpuMergeSort(int *output, int *index, int *input, int length);

//********************************* simulation Box ****************************************//
template <typename T> void cpuSetGlobalBox(T *h, T *h_inv, T *boxlo_bound, 
        T *boxhi_bound, T *boxhi, T *boxlo, T *boxtilt, int triclinic);
template <typename T> void cpuLamda2Box(T *x, T *lambda, T *h, T *boxlo, int dim, int n);
template <typename T> void cpuBox2Lamda(T *lambda, T *x, T *h_inv, T *boxlo, int dim, int n);
template <typename T> void cpuInsideBox(T *inside, T *x, T *boxlo, T *boxhi, T *lo_lamda, T *hi_lamda, 
        T *h_inv, int triclinic, int dim, int n);
template <typename T> void cpuMinimumImage(T *dp, T *h, int *pbc, int triclinic, int dim, int n);
template <typename T> void cpuUnmap(T *y, T *x, T *h, int *image, int triclinic, int dim, int n);

//********************************* Neighbor Lists ****************************************//
template <typename T> void cpuBoundingBox2D(T *vc, T*wc, T *v, T *w, T *a, T *b, T *r, int *pbc);
template <typename T> int cpuPeriodicImages2D(T *pimages, T *a, T *b, int *pbc);
template <typename T> void cpuBoundingBox3D(T *vc, T *wc, T *v, T *w, T *a, T *b, T *c, T *r, int *pbc);
template <typename T> int cpuPeriodicImages3D(T *pimages, T *a, T *b, T *c, int *pbc);
template <typename T> void cpuMakeReferenceGrid(T *eta, T smin, T smax, T ds, int nc);

void cpuCell2AtomList(int *c2alist, int *c2anumsum, int *c2anum, int *clist, int natom, int ncell);
template <typename T> void cpuAtomList2D(int *alist, int *inside, int *glistnumsum, int *glistnum, 
        T *x, T *pimages, T *wc, T *s2rmap, int inum, int pnum, int dim);
template <typename T> void cpuCellList2D(int *clist, T *x, T *eta1, T *eta2, T *eta3, T *s2rmap, 
        int *nc, int inum, int natom, int dim);
template <typename T> void cpuFullNeighborList2D(int *neighborlist, int *neighbornum, T *x, T *rcutsq, int *alist, 
        int *clist, int *c2alist, int *c2anumsum, int *nc, int inum, int jnum, int dim);
// template <typename T> void cpuFullNeighborList2D(int *neighborlist, int *neighbornum, T *x, T *rcutsq, int *alist, 
//         int *clist, int *c2alist, int *c2anumsum, int *atomtype, int *nc, int ntype, int inum, int jnum, int dim);

template <typename T> void cpuAtomList3D(int *alist, int *inside, int *glistnumsum, int *glistnum, 
        T *x, T *pimages, T *wc, T *s2rmap, int inum, int pnum, int dim);
template <typename T> void cpuCellList3D(int *clist, T *x, T *eta1, T *eta2, T *eta3, T *s2rmap, 
        int *nc, int inum, int natom, int dim);
template <typename T> void cpuFullNeighborList3D(int *neighborlist, int *neighbornum, T *x, T *rcutsq, int *alist, 
        int *clist, int *c2alist, int *c2anumsum, int *nc, int inum, int jnum, int dim);
// template <typename T> void cpuFullNeighborList3D(int *neighborlist, int *neighbornum, T *x, T *rcutsq, int *alist, 
//         int *clist, int *c2alist, int *c2anumsum, int *atomtype, int *nc, int ntype, int inum, int jnum, int dim);

template <typename T> void cpuFullNeighborList2D(int *neighborlist, int *neighbornum, T *x, T *rcutsq, 
        int anum, int inum, int jnum, int dim);
template <typename T> void cpuFullNeighborList3D(int *neighborlist, int *neighbornum, T *x, T *rcutsq, 
        int anum, int inum, int jnum, int dim);

template <typename T> void cpuFullNeighPairList(int *pairnum, int *pairlist, T *x, T *rcutsq, int *ilist, int *neighlist, 
        int *neighnum, int inum, int jnum, int dim);
template <typename T> void cpuFullNeighPairList(int *pairnum, int *pairlist, T *x, T *rcutsq,  int *atomtype, int *ilist, int *alist, 
        int *neighlist, int *neighnum, int inum, int jnum, int typej, int dim);
template <typename T> void cpuFullNeighPairList(int *pairnum, int *pairlist, T *x, T *rcutsq,  int *atomtype, int *ilist, int *alist, 
        int *neighlist, int *neighnum, int inum, int jnum, int typej, int typek, int dim);
template <typename T> void cpuFullNeighPairList(int *pairnum, int *pairlist, T *x, T *rcutsq,  int *atomtype, int *ilist, int *alist, 
        int *neighlist, int *neighnum, int inum, int jnum, int typej, int typek, int typel, int dim);
template <typename T> void cpuHalfNeighPairList(int *pairnum, int *pairlist, T *x, T *rcutsq, int *ilist, int *alist, 
        int *neighlist, int *neighnum, int inum, int jnum, int dim);
template <typename T> void cpuHalfNeighPairList(int *pairnum, int *pairlist, T *x, T *rcutsq, int *atomtype, int *ilist, 
        int *alist, int *neighlist, int *neighnum, int inum, int jnum, int typej, int dim);
void cpuNeighTripletList(int *tripletlist, int *tripletnumsum, int *pairnum, 
        int *pairlist, int *ilist, int *alist, int inum, int jnum);
void cpuNeighQuadrupletList(int *quadrupletlist, int *quadrupletnumsum, int *pairnum, 
        int *pairlist, int *ilist, int *alist, int inum, int jnum);

template <typename T> void cpuNeighTripletList(int *tripletnum, int *tripletlist, T *x, T *rcutsq, int *pairnum, int *pairnumsum, 
        int *pairlist, int *atomtype, int *ilist, int *alist, int *neighlist, int *neighnum, int inum, 
        int jnum, int typek, int dim);
template <typename T> void cpuNeighTriplets(T *xij, T *xik, T *qi, T *qj, T *qk, T *x, T *q, int *ai, int *aj, int *ak,  
      int *ti, int *tj, int *tk, int *tripletnum, int *tripletlist, int *tripletnumsum, int *alist,  
      int *atomtype, int jnum, int ijnum, int ncq, int dim);

template <typename T> void cpuNeighSingles(T *xi, T *qi, T *x, T *q, int *ai, 
      int *ti, int *ilist, int *atomtype, int inum, int ncq, int dim);
template <typename T> void cpuNeighPairs(T *xij, T *qi, T *qj, T *x, T *q, int *ai, int *aj,  
      int *ti, int *tj, int *pairnum, int *pairlist, int *pairnumsum, int *ilist, int *alist,  
      int *atomtype, int inum, int jnum, int ncq, int dim);
template <typename T> void cpuNeighTriplets(T *xij, T *xik, T *qi, T *qj, T *qk, T *x, T *q, int *ai, int *aj, int *ak,  
      int *ti, int *tj, int *tk, int *tripletnum, int *tripletlist, int *tripletnumsum, int *ilist, int *alist,  
      int *atomtype, int inum, int ncq, int dim);
template <typename T> void cpuNeighQuadruplets(T *xij, T *xik, T *xil, T *qi, T *qj, T *qk, T *ql, T *x, T *q, 
      int *ai, int *aj, int *ak, int *al, int *ti, int *tj, int *tk, int *tl, int *quadrupletnum, int *quadrupletlist, 
      int *quadrupletnumsum, int *ilist, int *alist,  int *atomtype, int inum, int ncq, int dim);

//********************************* Per-Atom Enegry, Force, Virial Tally ****************************************//
template <typename T> void cpuSingleDecomposition(T *e, T *ei, int *ai, int inum, int dim);
template <typename T> void cpuFullNeighPairDecomposition(T *e, T *eij, int *ai, int ijnum, int dim);
template <typename T> void cpuCenterAtomPairDecomposition(T *e, T *eij, int *ilist, int *anumsum, int inum, int dim);
template <typename T> void cpuHalfNeighPairDecomposition(T *e, T *eij, int *ai, int *aj, int ijnum, int dim);
template <typename T> void cpuNeighborAtomPairDecomposition(T *e, T *eij, int *jlist, int *bnumsum, int *index, int jnum, int dim);
template <typename T> void cpuTripletDecomposition(T *e, T *eijk, int *ai, int *aj, int *ak, int ijknum, int dim);
template <typename T> void cpuCenterAtomTripletDecomposition(T *e, T *eijk, int *ilist, int *anumsum, int inum, int dim);
template <typename T> void cpuNeighborAtomTripletDecomposition(T *e, T *eij, int *jlist, int *bnumsum, int *index, int jnum, int dim);
template <typename T> void cpuQuadrupletDecomposition(T *e, T *eijkl, int *ai, int *aj, int *ak, int *al, int ijklnum, int dim);
template <typename T> void cpuCenterAtomQuadrupletDecomposition(T *e, T *eijkl, int *ilist, int *anumsum, int inum, int dim);
template <typename T> void cpuNeighborAtomQuadrupletDecomposition(T *e, T *eij, int *jlist, int *bnumsum, int *index, int jnum, int dim);

template <typename T> void cpuSingleDecomposition(T *e, T *f, T *ei, T *fi, int *ai, int inum, int dim);
template <typename T> void cpuCenterAtomPairDecomposition(T *e, T *f, T *eij, T *fij, int *ilist, int *anumsum, int inum, int dim);
template <typename T> void cpuNeighborAtomPairDecomposition(T *e, T *f, T *eij, T *fij, int *jlist, int *bnumsum, int *index, int jnum, int dim);
template <typename T> void cpuFullNeighPairDecomposition(T *e, T *f, T *eij, T *fij, int *ai, int ijnum, int dim);
template <typename T> void cpuHalfNeighPairDecomposition(T *e, T *f, T *eij, T *fij, int *ai, int *aj, int ijnum, int dim);
template <typename T> void cpuTripletDecomposition(T *e, T *f, T *eijk, T *fij, T *fik, int *ai, int *aj, int *ak, int ijknum, int dim);
template <typename T> void cpuCenterAtomTripletDecomposition(T *e, T *f, T *eijk, T *fij, T *fik, int *ilist, int *anumsum, int inum, int dim);
template <typename T> void cpuNeighborAtomTripletDecomposition(T *e, T *f, T *eij, T *fij, int *jlist, int *bnumsum, int *index, int jnum, int dim);
template <typename T> void cpuQuadrupletDecomposition(T *e, T *f, T *eijkl, T *fij, T *fik, T *fil, 
        int *ai, int *aj, int *ak, int *al, int ijklnum, int dim);
template <typename T> void cpuCenterAtomQuadrupletDecomposition(T *e, T *f, T *eijkl, T *fij, T *fik, T *fil, 
        int *ilist, int *anumsum, int inum, int dim);
template <typename T> void cpuNeighborAtomQuadrupletDecomposition(T *e, T *f, T *eij, T *fij, int *jlist, int *bnumsum, int *index, int jnum, int dim);

template <typename T> void cpuHalfForceDecomposition(T *f, T *fij, int *ai, int *aj, int dim, int ijnum);
template <typename T> void cpuIAtomDecomposition(T *f, T *fij, int *ilist, int *anumsum, int dim, int inum);
template <typename T> void cpuJAtomDecomposition(T *f, T *fij, int *jlist, int *bnumsum, int *index, int dim, int jnum);

template <typename T> void cpuVirialPairTally(T *vatom, T *fij, T *rij, T factor, 
        int *ai, int dim, int inum, int ijnum);
template <typename T> void cpuVirialPairTally(T *vatom, T *fij, T *rij, T factor, 
        int *ai, int *aj, int dim, int inum, int ijnum);
template <typename T> void cpuVirialTripletTally(T *vatom, T *fij, T *fik, T *rij, T *rik, T factor, 
        int *ai, int *aj, int *ak, int dim, int inum, int ijnum);
template <typename T> void cpuVirialQuadrupletTally(T *vatom, T *fij, T *fik, T *fim, 
        T *rij, T *rik, T *rim, T factor, int *ai, int *aj, int *ak, int *am, int dim, int inum, int ijnum);

// template <typename T> void cpuHalfForceDecomposition2D(T *fi, T *fij, int *ai, int *aj, int ijnum);
// template <typename T> void cpuIAtomDecomposition2D(T *fi, T *fij, int *ilist, int *anumsum, int inum);
// template <typename T> void cpuJAtomDecomposition2D(T *fi, T *fij, int *jlist, int *bnumsum, int *index, int jnum);
// 
// template <typename T> void cpuHalfForceDecomposition3D(T *fi, T *fij, int *ai, int *aj, int ijnum);
// template <typename T> void cpuIAtomDecomposition3D(T *fi, T *fij, int *ilist, int *anumsum, int inum);
// template <typename T> void cpuJAtomDecomposition3D(T *fi, T *fij, int *jlist, int *bnumsum, int *index, int jnum);

//***************************  Compute Outputs ***************************************//
template <typename T> void cpuPackIntProperty(T *buf, int *prop, int *ilist, 
         int m, int mvalues, int n, int nvalues, int inum);
template <typename T> void cpuPackIntProperty(T *buf, int *prop, int *type, int *ilist, 
         int n, int nvalues, int inum);
template <typename T> void cpuPackFloatProperty(T *buf, T *prop, int *ilist, 
         int m, int mvalues, int n, int nvalues, int inum);
template <typename T> void cpuPackFloatProperty(T *buf, T *prop, T a, T b, int *ilist, 
         int m, int mvalues, int n, int nvalues, int inum);
template <typename T> void cpuPackFloatProperty(T *buf, T *prop, int *type, int *ilist, 
         int n, int nvalues, int inum);
template <typename T> T cpuComputeMass(T *amass, T *mass, int *type, int *ilist, int inum);
template <typename T> void  cpuComputeXCM(T *xcm, T *x, T *mass, T *box, T masstotal, 
        int *ilist, int *type, int *image, int triclinic, int dim, int inum);
template <typename T> void  cpuComputeVCM(T *vcm, T *v, T *mass, T masstotal, 
        int *ilist, int *type, int dim, int inum);
template <typename T> T cpuComputeGyration(T *xcm, T *x, T *mass, T *box, T masstotal, 
        int *ilist, int *type, int *image, int triclinic, int dim, int inum);
template <typename T> void cpuComputeAngmom(T *p, T *xcm, T *x, T *v, T *mass, T *box, 
        int *ilist, int *type, int *image, int triclinic, int dim, int inum);
template <typename T> void cpuComputeTorque(T *tlocal, T *xcm, T *x, T *f, T *box, 
        int *ilist, int *image, int triclinic, int dim, int inum);
template <typename T> void cpuComputeInertia(T *ione, T *xcm, T *x, T *mass, T *box, 
        int *ilist, int *type, int *image, int triclinic, int dim, int inum);
template <typename T> T cpuComputeNVTEnergy(T *tarray, T *eta, T *eta_mass, T *eta_dot, int mtchain);
template <typename T> T cpuComputeTempScalar(T *v, T *mass, T tfactor, int *type, int *ilist, 
         int dim, int inum);
template <typename T> void cpuComputeTempSymTensor(T *t, T *v, T *mass, T tfactor, int *type, int *ilist, 
         int dim, int inum);
template <typename T> T cpuComputePressureScalar(T *virial, T volume, T temp, T tempdof, 
        T boltz, T nktv2p, int dim);
template <typename T> void cpuComputePressureSymTensor(T *vector, T *virial, T *ke_tensor, 
        T volume, T nktv2p, int dim);
template <typename T> void cpuComputeHeatFlux(T *vector, T *ke, T *pe, T *stress, T *v, 
        T nktv2p, int *ilist,  int pressatomflag, int dim, int inum);
template <typename T> void cpuComputeKEAtom(T *ke, T *mass,  T *v, 
        T mvv2e, int *type, int *ilist,  int dim, int inum);
template <typename T> void cpuComputeStressAtom(T *stress, T *mass, T *vatom, T *v, 
        T mvv2e, T nktv2p, int *type, int *ilist,  int vflag, int keflag, int dim, int inum);
template <typename T> void cpuComputeCentroidStressAtom(T *stress, T *mass, T *cvatom, T *v, 
        T mvv2e, T nktv2p, int *type, int *ilist, int vflag, int keflag, int dim, int inum);
template <typename T> void cpuComputeDisplaceAtom(T *displace, T *x, T *xoriginal, T *box, 
        int *pbc, int *ilist,  int triclinic, int dim, int inum);
template <typename T> void cpuComputeOrientOrderAtom(T* qnarray, T *x, T *rlist, T *cglist, T *fac, 
        T *qnm_r, T *qnm_i, T *distsq, T cutsq, T MY_EPSILON, T QEPSILON, T MY_4PI, int *neighlist, 
        int *neighnum,  int *ilist, int *qlist, int *nearest,  int nqlist, int qmax, int wlflag,  
        int wlhatflag,  int qlcompflag, int iqlcomp, int qlcomp, int nnn, int jnum, int dim, int inum); 
template <typename T> void cpuComputeCoordAtomCutoff(int *cvec, T *x, T *rcutsq, int *type, 
        int *ilist, int *neighlist, int *neighnum, int *typelo, int *typehi, 
         int *jgroupbit, int dim, int ntypes, int jnum, int inum);
template <typename T> void cpuComputeCoordAtomCutoff(int *carray, T *x, T *rcutsq, int *type, 
        int *ilist, int *neighlist, int *neighnum, int *typelo, int *typehi, 
         int *jgroupbit, int ncol, int dim, int ntypes, int jnum, int inum);
template <typename T> void cpuComputeCoordAtomOrient(int *cvec, T *x, T *rcutsq, T *normv, T threshold, 
        int *type, int *ilist, int *neighlist, int *neighnum, int *typelo, int *typehi, 
         int *jgroupbit, int dim, int ntypes, int nqlist, int ncol, int l, int jnum, int inum);
template <typename T> void cpuComputeMSD(T *vector, T *x, T *xoriginal, T *h, T *xcm,
         int *ilist, int *image, int naverage, int avflag, int triclinic, int nmsd, int dim,  int inum);
template <typename T> void cpuComputeVACF(T *vacf, T *v, T *voriginal, 
         int *ilist, int nvacf, int dim,  int inum);

//********************************* Velocity integration ****************************************//
template <typename T> void cpuVelocityZeroMomentum(T *v, T *vcm, int dim, int nlocal);
template <typename T> void cpuVelocityZeroRotation(T *x, T *v, T *box, T *xcm, T *omega, 
        int *image, int triclinic, int dim, int nlocal);
template <typename T> void cpuVelocitySet(T *v, T *vext, int *vdim,
        int sum_flag, int dim, int nlocal);
template <typename T> void cpuVelocityRamp(T *x, T *v, T *v_lo, T *v_hi, T *coord_lo, T *coord_hi,
        int *coord_dim, int *v_dim, int sum_flag, int dim, int nlocal);
template <typename T> void cpuVelocityCreate(T *x, T *v, T *mass, T *second, 
         int *seed, int *save, int *map, int *type, int seed0, int sum_flag, int dist_flag, int loop_flag, 
        int dim, int mpiRank, int nlocal, int natoms);

template <typename T> void cpuVelocity(T *x, T *v, T *f, T *box, T *xcm, T *vcm, 
        T *mass, T *second, T *omega, T *vext, T *v_lo, T *v_hi, T *coord_lo, T *coord_hi, 
        T t_desired, T t_current, int *seed, int *save, int *map, int *image, int *type, 
        int *coord_dim, int *vdim, int sum_flag, int dist_flag, int loop_flag, 
        int rotation_flag, int momentum_flag, int triclinic, int dim, int mpiRank, 
        int vmode, int nlocal, int natoms);

template <typename T> void cpuSetVelocityInitialIntegrate(T *x, T *v, T *f, T *mass, T *fparam,
        T dtf, T dtv, int *type, int *ilist, int *iparam, int dim, int inum);

template <typename T> void cpuSetVelocityFinalIntegrate(T *x, T *v, T *f, T *mass, 
        T dtf, int *type, int *ilist, int *iparam, int dim, int inum);

template <typename T> void cpuInitialIntegrate(T *x, T *v, T *f, T *mass, T *dtarray, T *tarray,
        T *eta_mass, T *eta, T *eta_dot, T *eta_dotdot, T vlimitsq, int *type, int *ilist,
        int eta_mass_flag, int biasflag, int mtchain, int nc_tchain, int mode, int dim, int inum);

template <typename T> void cpuFinalIntegrate(T *x, T *v, T *f, T *mass, T *dtarray, T *tarray,
        T *eta_mass, T *eta, T *eta_dot, T *eta_dotdot, T vlimitsq, int *type, int *ilist, 
        int eta_mass_flag, int biasflag, int mtchain, int nc_tchain, int mode, int dim, int inum);

template <typename T> T cpuVelocityRescalingThermostat(T *v, T *mass, T *dtarray, T *tarray, T *second, 
        T energy, int *type, int *ilist, int *seed, int *save, int biasflag, int mode, int dim, int inum);

template <typename T> void cpuPBC(T *x, T *v, int *image, T *boxhi, T *boxlo, T *hi_lambda, T *lo_lambda,  
        T *h, T *h_inv, T *h_rate, int *pbc, int vdeform, int triclinic, int dim, int nlocal);

//********************************* Force Constraints ****************************************//
template <typename T> void cpuFixSetForce(T *x, T *v, T *f, T *eatom, T *vatom, T *fparam, 
        int *iparam, int *ilist, int eflag_atom, int vflag_atom, int dim, int inum);
template <typename T> void cpuFixLineForce(T *x, T *v, T *f, T *eatom, T *vatom, T *fparam, 
        int *iparam, int *ilist, int eflag_atom, int vflag_atom, int dim, int inum);
template <typename T> void cpuFixPlaneForce(T *x, T *v, T *f, T *eatom, T *vatom, T *fparam, 
        int *iparam, int *ilist, int eflag_atom, int vflag_atom, int dim, int inum);
template <typename T> void cpuFixAveForce(T *x, T *v, T *f, T *eatom, T *vatom, T *fparam, 
        int *iparam, int *ilist, int eflag_atom, int vflag_atom, int dim, int inum);

template <typename T> void cpuFixAddForce(T *x, T *v, T *f, T *eatom, T *vatom, T *fparam, T *box,
        int *iparam, int *ilist, int *image, int triclinic, int eflag_atom, int vflag_atom, int dim, int inum);
template <typename T> void cpuFixDragForce(T *x, T *v, T *f, T *eatom, T *vatom, T *fparam, T *box,
        int *iparam, int *ilist, int *pbc, int triclinic, int eflag_atom, int vflag_atom, int dim, int inum);

template <typename T> void cpuFixWallReflect(T *x, T *v, T *f, T *eatom, T *vatom, T *fparam, 
        int *iparam, int *ilist, int eflag_atom, int vflag_atom, int dim, int inum);
template <typename T> void cpuFixWallHarmonic(T *x, T *v, T *f, T *eatom, T *vatom, T *fparam, 
        int *iparam, int *ilist, int eflag_atom, int vflag_atom, int dim, int inum);
template <typename T> void cpuFixWallLJ93(T *x, T *v, T *f, T *eatom, T *vatom, T *fparam, 
        int *iparam, int *ilist, int eflag_atom, int vflag_atom, int dim, int inum);
template <typename T> void cpuFixWallLJ126(T *x, T *v, T *f, T *eatom, T *vatom, T *fparam, 
        int *iparam, int *ilist, int eflag_atom, int vflag_atom, int dim, int inum);
template <typename T> void cpuFixWallLJ1043(T *x, T *v, T *f, T *eatom, T *vatom, T *fparam, 
        int *iparam, int *ilist, int eflag_atom, int vflag_atom, int dim, int inum);
template <typename T> void cpuFixWallMorse(T *x, T *v, T *f, T *eatom, T *vatom, T *fparam, 
        int *iparam, int *ilist, int eflag_atom, int vflag_atom, int dim, int inum);

//***************************  Empirical Potentials ***************************************//
template <typename T> void opuSinglea(T *u, T *xi, T *qi, int *ti, int *ai, T *mu, T *eta, int *kappa, int dim, int ncq, int nmu, int neta, int nkappa, int ng, int potnum);
template <typename T> void opuSingleb(T *u, T *xi, T *qi, int *ti, int *ai, T *mu, T *eta, int *kappa, int dim, int ncq, int nmu, int neta, int nkappa, int ng, int potnum);
template <typename T> void opuPaira(T *u, T *xij, T *qi, T *qj, int *ti, int *tj, int *ai, int *aj, T *mu, T *eta, int *kappa, int dim, int ncq, int nmu, int neta, int nkappa, int ng, int potnum);
template <typename T> void opuPairb(T *u, T *xij, T *qi, T *qj, int *ti, int *tj, int *ai, int *aj, T *mu, T *eta, int *kappa, int dim, int ncq, int nmu, int neta, int nkappa, int ng, int potnum);
template <typename T> void opuPairc(T *u, T *xij, T *qi, T *qj, int *ti, int *tj, int *ai, int *aj, T *mu, T *eta, int *kappa, int dim, int ncq, int nmu, int neta, int nkappa, int ng, int potnum);
template <typename T> void opuPaircDensity(T *u, T *rho, T *mu, T *eta, int *kappa, int nrho, int nmu, int neta, int nkappa, int ng, int potnum);
template <typename T> void opuTripleta(T *u, T *xij, T *xik, T *qi, T *qj, T *qk, int *ti, int *tj, int *tk, int *ai, int *aj, int *ak, T *mu, T *eta, int *kappa, int dim, int ncq, int nmu, int neta, int nkappa, int ng, int potnum);
template <typename T> void opuTripletb(T *u, T *xij, T *xik, T *qi, T *qj, T *qk, int *ti, int *tj, int *tk, int *ai, int *aj, int *ak, T *mu, T *eta, int *kappa, int dim, int ncq, int nmu, int neta, int nkappa, int ng, int potnum);
template <typename T> void opuTripletc(T *u, T *xij, T *xik, T *qi, T *qj, T *qk, int *ti, int *tj, int *tk, int *ai, int *aj, int *ak, T *mu, T *eta, int *kappa, int dim, int ncq, int nmu, int neta, int nkappa, int ng, int potnum);
template <typename T> void opuTripletcPair(T *u, T *xij, T *qi, T *qj, int *ti, int *tj, int *ai, int *aj, T *mu, T *eta, int *kappa, int dim, int ncq, int nmu, int neta, int nkappa, int ng, int potnum);
template <typename T> void opuTripletcDensity(T *u, T *rho, T *mu, T *eta, int *kappa, int nrho, int nmu, int neta, int nkappa, int ng, int potnum);
template <typename T> void opuQuadrupleta(T *u, T *xij, T *xik, T *xil, T *qi, T *qj, T *qk, T *ql, int *ti, int *tj, int *tk, int *tl, int *ai, int *aj, int *ak, int *al, T *mu, T *eta, int *kappa, int dim, int ncq, int nmu, int neta, int nkappa, int ng, int potnum);
template <typename T> void opuQuadrupletb(T *u, T *xij, T *xik, T *xil, T *qi, T *qj, T *qk, T *ql, int *ti, int *tj, int *tk, int *tl, int *ai, int *aj, int *ak, int *al, T *mu, T *eta, int *kappa, int dim, int ncq, int nmu, int neta, int nkappa, int ng, int potnum);

template <typename T> void opuSingleaGradient(T *u, T *u_xi, T *xi, T *qi, int *ti, int *ai, T *mu, T *eta, int *kappa, int dim, int ncq, int nmu, int neta, int nkappa, int ng, int potnum);
template <typename T> void opuSinglebGradient(T *u, T *u_xi, T *xi, T *qi, int *ti, int *ai, T *mu, T *eta, int *kappa, int dim, int ncq, int nmu, int neta, int nkappa, int ng, int potnum);
template <typename T> void opuPairaGradient(T *u, T *u_xij, T *xij, T *qi, T *qj, int *ti, int *tj, int *ai, int *aj, T *mu, T *eta, int *kappa, int dim, int ncq, int nmu, int neta, int nkappa, int ng, int potnum);
template <typename T> void opuPairbGradient(T *u, T *u_xij, T *xij, T *qi, T *qj, int *ti, int *tj, int *ai, int *aj, T *mu, T *eta, int *kappa, int dim, int ncq, int nmu, int neta, int nkappa, int ng, int potnum);
template <typename T> void opuPaircGradient(T *u, T *u_xij, T *xij, T *qi, T *qj, int *ti, int *tj, int *ai, int *aj, T *mu, T *eta, int *kappa, int dim, int ncq, int nmu, int neta, int nkappa, int ng, int potnum);
template <typename T> void opuPaircDensityGradient(T *u, T *u_rho, T *rho, T *mu, T *eta, int *kappa, int nrho, int nmu, int neta, int nkappa, int ng, int potnum);
template <typename T> void opuTripletaGradient(T *u, T *u_xij, T *u_xik, T *xij, T *xik, T *qi, T *qj, T *qk, int *ti, int *tj, int *tk, int *ai, int *aj, int *ak, T *mu, T *eta, int *kappa, int dim, int ncq, int nmu, int neta, int nkappa, int ng, int potnum);
template <typename T> void opuTripletbGradient(T *u, T *u_xij, T *u_xik, T *xij, T *xik, T *qi, T *qj, T *qk, int *ti, int *tj, int *tk, int *ai, int *aj, int *ak, T *mu, T *eta, int *kappa, int dim, int ncq, int nmu, int neta, int nkappa, int ng, int potnum);
template <typename T> void opuTripletcGradient(T *u, T *u_xij, T *u_xik, T *xij, T *xik, T *qi, T *qj, T *qk, int *ti, int *tj, int *tk, int *ai, int *aj, int *ak, T *mu, T *eta, int *kappa, int dim, int ncq, int nmu, int neta, int nkappa, int ng, int potnum);
template <typename T> void opuTripletcPairGradient(T *u, T *u_xij, T *xij, T *qi, T *qj, int *ti, int *tj, int *ai, int *aj, T *mu, T *eta, int *kappa, int dim, int ncq, int nmu, int neta, int nkappa, int ng, int potnum);
template <typename T> void opuTripletcDensityGradient(T *u, T *u_rho, T *rho, T *mu, T *eta, int *kappa, int nrho, int nmu, int neta, int nkappa, int ng, int potnum);
template <typename T> void opuQuadrupletaGradient(T *u, T *u_xij, T *u_xik, T *u_xil, T *xij, T *xik, T *xil, T *qi, T *qj, T *qk, T *ql, int *ti, int *tj, int *tk, int *tl, int *ai, int *aj, int *ak, int *al, T *mu, T *eta, int *kappa, int dim, int ncq, int nmu, int neta, int nkappa, int ng, int potnum);
template <typename T> void opuQuadrupletbGradient(T *u, T *u_xij, T *u_xik, T *u_xil, T *xij, T *xik, T *xil, T *qi, T *qj, T *qk, T *ql, int *ti, int *tj, int *tk, int *tl, int *ai, int *aj, int *ak, int *al, T *mu, T *eta, int *kappa, int dim, int ncq, int nmu, int neta, int nkappa, int ng, int potnum);

//***************************  spherical harmonic potentials ***************************************//
void cpuGetIndk(int *indk, int K);
template <typename T> T clebschgordan(int j1, int m1, int j2, int m2, int j, int m, T *fac);
int cgcoefficients(int *indl, int N);
template <typename T> void cgcoefficients(T *cg, int *indm, int *rowm, int *indl, T *fac, int M, int N);
template <typename T> void cpuSphericalHarmonicsBispectrum(T *b, T *Sr, T *Si, T *fac, int L);
template <typename T> void cpuSphericalHarmonicsAtomSum(T *Ylmr, T *Ylmi, T *x, T *y, T *z, 
                T *P, T *tmp, T *fac, T pi, int L, int N);
template <typename T> int cpuSphericalHarmonicsBispectrumIndex(T *b, int L);
template <typename T> void cpuSphericalHarmonicsBispectrumIndex(int *indl, T *b, int M, int L);

template <typename T> void cpuSphericalHarmonicsBessel(T *Sr, T *Si, T *xij,
                T *x0, T *P, T *tmp, T *f, T *fac, T pi, int L, int K, int N);
template <typename T> void cpuSphericalHarmonicsBesselWithDeriv(T *Sr, T *Si, T *Srx, T *Six, T *Sry, T *Siy, T *Srz, T *Siz, T *xij, 
                T *x0, T *P, T *tmp, T *f, T *dP, T *dtmp, T *df, T *fac, T pi, int L, int K, int N);
template <typename T> void cpuRadialSphericalHarmonicsSpectrum(T *c, T *ar, T *ai, T *sr, T *si, T *cg, int *indk, 
        int *indl, int *indm, int *rowm, int *Nnb, int Nub, int Ncg, int Na, int L, int K, int spectrum);
template <typename T> void cpuRadialSphericalHarmonicsSpectrumDeriv(T *cd, T *ar, T *ai, 
        T *srx, T *six, T *sry, T *siy, T *srz, T *siz, T *cg, int *indk, int *indl, int *indm, 
        int *rowm, int *Nnb, int Nub, int Ncg, int Na, int L, int K, int spectrum);
template <typename T> void cpuForceDecomposition(T *f, T *fij, int *ai, int *aj, 
        int inum, int ijnum, int Nbf);
template <typename T> void cpuForceDecomposition(T *f, T *fij, int *ai, int *aj, int ijnum);
template <typename T> void cpuCenterAtomDecomposition(T *f, T *fij, int *ilist, int *anumsum, 
        int inum, int ijnum, int Na, int Nbf);
template <typename T> void cpuNeighborAtomDecomposition(T *f, T *fij, int *jlist, int *bnumsum, int *index, 
        int inum, int ijnum, int jnum, int Nbf);
template <typename T> void cpuCenterAtomDecomposition(T *f, T *fij, int *ilist, int *anumsum, int inum);
template <typename T> void cpuNeighborAtomDecomposition(T *f, T *fij, int *jlist, int *bnumsum, int *index, int jnum);
template <typename T> void cpuCenterAtomDecomposition(T *e, T *ei, int *ilist, int Na);

//*************************** SNAP POTENTIAL ***********************************************//
void cpuBuildIndexList(int *idx_max, int *idxz, int *idxz_block, int *idxb, int *idxb_block, int *idxu_block, 
        int *idxcg_block, int twojmax);

template <typename T> void cpuInitRootpqArray(T *rootpqarray, int twojmax);
template <typename T> void cpuInitClebschGordan(T *cglist, T *factorial, int twojmax);

template <typename T> void cpuInitSna(T *rootpqarray, T *cglist, T *factorial, int *idx_max, int *idxz, 
      int *idxz_block, int *idxb, int *idxb_block, int *idxu_block, int *idxcg_block, int twojmax);

template <typename T> void cpuZeroUarraytot(T *ulisttot_r, T *ulisttot_i, T wself, int *idxu_block, int *type,
        int *map, int *ai, int wselfall_flag, int chemflag, int idxu_max, int nelements, int twojmax, int inum);

template <typename T> void cpuAddUarraytot(T *ulisttot_r, T *ulisttot_i, T *ulist_r, T *ulist_i, T *rij, 
        T *wjelem, T *radelem, T rmin0, T rcutfac, int *idxu_block, int *ilist, 
        int *type, int *pairnum, int *pairnumsum, int *map, int *tj, int twojmax, int idxu_max, 
        int nelements, int inum, int switch_flag, int chemflag);

template <typename T> void cpuComputeBeta2(T *beta, T *bispectrum, T *coeffelem, int *ilist, int *map, int *type, 
        int inum, int ncoeff, int ncoeffall, int quadraticflag);

// template <typename T> void cpuComputeUij(T *ulist_r, T *ulist_i, T *rootpqarray, T *rij, 
//         T *radelem, T rmin0, T rfac0, T rcutfac, int *idxu_block, 
//         int *type, int *alist, int *ai, int *aj, int twojmax, int idxu_max, int ijnum);
 
template <typename T> void cpuComputeUij(T *ulist_r, T *ulist_i, T *rootpqarray, T *rij, 
        T *radelem, T rmin0, T rfac0, T rcutfac, int *idxu_block, 
        int *ai, int *aj, int *ti, int *tj, int twojmax, int idxu_max, int ijnum);

template <typename T> void cpuComputeZi(T *zlist_r, T *zlist_i, T *ulisttot_r, T *ulisttot_i, T *cglist,
        int *idxz, int *idxu_block, int *idxcg_block, int twojmax, int idxu_max, int idxz_max, int nelements, 
        int bnorm_flag, int inum);

template <typename T> void cpuComputeYi(T *ylist_r, T *ylist_i, T *ulisttot_r, T *ulisttot_i, T *cglist, 
        T* betalist, int *idxz, int *idxb_block, int *idxu_block, int *idxcg_block, int twojmax, 
        int idxb_max, int idxu_max, int idxz_max, int nelements, int bnorm_flag, int ncoeff, int inum);

template <typename T> void cpuComputeBi(T *blist, T *zlist_r, T *zlist_i, T *ulisttot_r, T *ulisttot_i, 
        T *bzero,int *ilist, int *type, int *map, int *idxb, int *idxu_block, int *idxz_block, 
        int twojmax, int idxb_max, int idxu_max, int idxz_max, int nelements, int bzero_flag, 
        int wselfall_flag, int chemflag, int inum);

// template <typename T> void cpuComputeDbidrj(T *dblist, T *zlist_r, T *zlist_i, T *dulist_r, T *dulist_i, 
//         int *idxb, int *idxu_block, int *idxz_block, int *type, int *map, int *ai, int *aj, int twojmax, 
//         int idxb_max, int idxu_max, int idxz_max, int nelements, int bnorm_flag, int chemflag, int ijnum);

template <typename T> void cpuComputeDbidrj(T *dblist, T *zlist_r, T *zlist_i, T *dulist_r, T *dulist_i, 
        int *idxb, int *idxu_block, int *idxz_block, int *map, int *ai, int *aj, int *ti, int *tj, int twojmax, 
        int idxb_max, int idxu_max, int idxz_max, int nelements, int bnorm_flag, int chemflag, int ijnum);

template <typename T> void cpuComputeDbidrj(T *dblist, T *zlist_r, T *zlist_i, 
        T *dulist_r, T *dulist_i, int *idxb, int *idxu_block, int *idxz_block, 
        int *map, int *ai, int *tj, int twojmax, int idxb_max, int idxu_max, int idxz_max, 
        int nelements, int bnorm_flag, int chemflag, int inum, int ijnum);

template <typename T> void cpuComputeDuijdrj(T *dulist_r, T *dulist_i, T *ulist_r, T *ulist_i, T *rootpqarray, 
        T* rij, T *wjelem, T *radelem, T rmin0, T rfac0, T rcutfac, int *idxu_block, 
        int *ai, int *aj, int *ti, int *tj, int twojmax, int idxu_max, int ijnum, int switch_flag);

template <typename T> void cpuComputeDeidrj(T *dedr, T *ylist_r, T *ylist_i, T *dulist_r, T *dulist_i,         
        int *idxu_block, int *map, int *ai, int *aj, int *ti, int *tj, int nelements, int twojmax, int idxu_max, 
        int chemflag, int ijnum) ;

template <typename T> void cpuComputeSna(T *sna, T *blist, int *ilist, int *mask, 
        int ncoeff, int nrows, int inum, int quadraticflag);

template <typename T> void cpuComputeSna(T *sna, T *blist, int *ilist, int *mask, int *type,
        int ncoeff, int ntype, int nperdim, int inum, int quadraticflag);

template <typename T> void cpuComputeSnad(T *snad, T *dblist, T *blist, int *aii, int *ai, int *aj, int *ti, 
        int *mask, int ncoeff, int ntypes, int nperdim, int ijnum, int quadraticflag);

template <typename T> void cpuComputeSnad(T *snad, T *dblist, T *blist, int *aii, int *ai, int *aj, int *ti, 
        int *mask, int *tag, int ncoeff, int ntypes, int nperdim, int ijnum, int quadraticflag);

template <typename T> void cpuComputeSnav(T *snav, T *dblist, T *blist, T *x, int *aii, int *ai, int *aj,
        int *ti, int *mask, int ncoeff, int ntypes, int nperdim, int ijnum, int quadraticflag);

template <typename T> void cpuComputeSnav2(T *snav, T *dblist, T *blist, T *x, int *aii, int *ai, int *aj,
        int *ti, int *mask, int ncoeff, int ntypes, int nperdim, int ijnum, int quadraticflag);

template <typename T> void cpuSnapTallyEnergyFull(T *eatom, T *bispectrum, T *coeffelem, int *ilist, int *map, int *type, 
        int inum, int ncoeff, int ncoeffall, int quadraticflag);

template <typename T> void cpuSnapTallyForceFull(T *fatom, T *fij, int *ai, int *aj, int *alist, int ijnum);

//template <typename T> void cpuSnapTallyVirialFull(T *vatom, T *fij, T *rij, int *ai, int *aj, int ijnum);
template <typename T> void cpuSnapTallyVirialFull(T *vatom, T *fij, T *rij, int *ai, int *aj, int inum, int ijnum);

template <typename T> void cpuNeighPairList(int *pairnum, int *pairlist, T *x, T rcutsq, int *ilist, int *neighlist, 
        int *neighnum, int inum, int jnum, int dim);

// template <typename T> void cpuNeighPairList(int *pairnum, int *pairlist, T *x, T *rcutsq, int *ilist, int *neighlist, 
//         int *neighnum, int *atomtype, int *alist, int inum, int jnum, int dim, int ntypes)

template <typename T> void cpuNeighPairList(int *pairnum, int *pairlist, T *x, T *rcutsq, int *ilist, int *neighlist, 
        int *neighnum, int *atomtype, int *alist, int inum, int jnum, int dim, int ntypes);

// template <typename T> void cpuNeighPairs(T *xij, T *x, int *aii, int *ai, int *aj,  
//       int *ti, int *tj, int *pairnum, int *pairlist, int *pairnumsum, int *ilist, 
//       int *atomtype, int inum, int jnum, int dim);

template <typename T> void cpuNeighPairs(T *xij, T *x, int *aii, int *ai, int *aj,  
      int *ti, int *tj, int *pairnum, int *pairlist, int *pairnumsum, int *ilist, 
      int *atomtype, int *alist, int inum, int jnum, int dim);

template <typename T> void cpuComputeSij(T *Sr, T *Si, T *Srx, T *Six, T *Sry, T *Siy, T *Srz, T *Siz, 
        T *rootpqarray, T *rij, T *wjelem, T *radelem, T rmin0, T rfac0, T rcutfac, int *idxu_block,  
        int *ti, int *tj, int twojmax, int idxu_max, int ijnum, int switch_flag);                

template <typename T> void cpuZeroUarraytot2(T *Stotr, T *Stoti, T wself, int *idxu_block, 
        int *type, int *map, int *ai, int wselfall_flag, int chemflag, int idxu_max, int nelements, 
         int twojmax, int inum);

template <typename T> void cpuAddUarraytot(T *Stotr, T *Stoti, T *Sr, 
        T *Si, int *pairnum, int *pairnumsum, int *map, int *tj, 
        int idxu_max, int nelements, int inum, int ijnum, int chemflag);

template <typename T> void cpuAddUarraytot(T *Stotr, T *Stoti, T *Sr, 
        T *Si, int *map, int *ai, int *tj, int idxu_max, int inum, int ijnum, int chemflag);

template <typename T> void cpuComputeZi2(T *zlist_r, T *zlist_i, T *Stotr, T *Stoti, 
        T *cglist, int *idxz, int *idxu_block, int *idxcg_block, int twojmax, int idxu_max, 
        int idxz_max, int nelements, int bnorm_flag, int inum);

template <typename T> void cpuComputeBi2(T *blist, T *zlist_r, T *zlist_i, T *Stotr, T *Stoti, 
        T *bzero, int *ilist, int *type, int *map, int *idxb, int *idxu_block, int *idxz_block, 
        int twojmax, int idxb_max, int idxu_max, int idxz_max, int nelements, int bzero_flag, 
        int wselfall_flag, int chemflag, int inum);

template <typename T> void cpuComputeBeta(T *beta, T *bispectrum, T *coeffelem, int *ilist, int *map, int *type, 
        int inum, int ncoeff, int ncoeffall, int quadraticflag);

template <typename T> void cpuComputeYi(T *ylist_r, T *ylist_i, T *zlist_r, T *zlist_i, 
        T* beta, int *idxz, int *idxb_block, int twojmax, int idxb_max, int idxu_max, int idxz_max,
        int nelements, int bnorm_flag, int inum);

template <typename T> void cpuComputeDeidrj(T *dedr, T *ylist_r, T *ylist_i, 
        T *dulist_r, T *dulist_i, int *idxu_block, int *map, int *ai, int *tj,
        int twojmax, int idxu_max, int chemflag, int inum, int ijnum); 

template <typename T> void cpuSnapTallyEnergyFull2(T *eatom, T *bispectrum, T *coeffelem, int *ilist, 
        int *map, int *type, int inum, int ncoeff, int ncoeffall, int quadraticflag);

template <typename T> void cpuSnapTallyBispectrum(T *bi, T *bispectrum, int *ilist, 
        int *type, int inum, int ncoeff, int nperdim, int ntype, int quadraticflag);

template <typename T> void cpuSnapTallyBispectrumDeriv(T *db, T *bispectrum, T *dbdr, int *aii, 
        int *ai, int *aj, int *ti, int inum, int ijnum, int ncoeff, int nperdim, int ntype, int quadraticflag);

template <typename T> void cpuSnapTallyBispectrumVirial(T *bv, T *bispectrum, T *dbdr, T *rij, int *aii, 
        int *ai, int *aj, int *ti, int inum, int ijnum, int ncoeff, int nperdim, int ntype, int quadraticflag);

template <typename T> void cpuSnapComputeUi(T *Stotr, T *Stoti, T *rootpqarray, T *rij, T *wjelem, T *radelem, 
        T rmin0, T rfac0, T rcutfac, int *idxu_block, int *map, int *ai, int *ti, int *tj, 
        int twojmax, int idxu_max, int inum, int ijnum, int switchflag, int chemflag);                

template <typename T> void cpuAddWself2Ui(T *Stotr, T *Stoti, T wself, int *idxu_block, 
        int *type, int *map, int *ai, int wselfall_flag, int chemflag, int idxu_max, int nelements, 
         int twojmax, int inum);

template <typename T> void cpuSnapComputeBi(T *blist, T *Stotr, T *Stoti, T *cglist, T *bzero, 
        int *idxb, int *idxcg_block, int *idxu_block, int *idxz_block, int twojmax, int idxb_max, int idxu_max, 
        int idxz_max, int nelements, int bnorm_flag, int bzero_flag, int wselfall_flag, int inum);

template <typename T> void cpuSnapComputeEi(T *eatom, T *Stotr, T *Stoti, T *cglist, 
        T *bzero, T *coeffelem, int *ilist, int *map, int *type, int *idxb, int *idxcg_block, int *idxu_block, 
        int twojmax, int idxb_max, int idxu_max, int nelements, int ncoeffall, int bnorm_flag, 
        int bzero_flag, int wselfall_flag, int quadraticflag, int inum);

template <typename T> void cpuSnapComputeYi(T *ylist_r, T *ylist_i, T *Stotr, T *Stoti, T *cglist, 
        T* beta, int *idxz, int *idxb_block, int *idxu_block, int *idxcg_block, int twojmax, 
        int idxb_max, int idxu_max, int idxz_max, int nelements, int bnorm_flag, int inum);

template <typename T> void cpuSnapComputeYi(T *ylist_r, T *ylist_i, T *Stotr, T *Stoti, T *cglist, T* coeffelem, 
        int *map, int *type, int *idxz, int *idxb_block, int *idxu_block, int *idxcg_block, int twojmax, 
        int idxb_max, int idxu_max, int idxz_max, int nelements, int ncoeffall, int bnorm_flag, int inum);

template <typename T> void cpuSnapComputeFi(T *fatom, T *vatom, T *ylist_r, T *ylist_i, T *rootpqarray, T *rij, 
        T *wjelem, T *radelem,  T rmin0, T rfac0, T rcutfac, int *idxu_block, int *map, int *aii, int *ai, int *aj, 
        int *ti, int *tj, int twojmax, int idxu_max, int inum, int anum, int ijnum, int switchflag, int chemflag); 

//*************************** Lattice, Region, Domain ***********************************************//
template <typename T> int cpuLattice(T *basis, T *primitive, T *rotaterow, T *priminv, T *rotatecol, T *origin, 
        T *spacing, T *a1, T *a2, T *a3, T &scale, int *orientx, int *orienty, int *orientz, 
        int style, int unit_style, int spaceflag, int dimension);
template <typename T> void cpuLatticeBoundingBox(T *lmin, T *lmax, T *bsublo, T *bsubhi, T *primitive, 
        T *rotaterow, T *priminv, T *rotatecol, T *origin, T *spacing, T scale);
int cpuLatticeCount(int nbasis, int ilo, int ihi, int jlo, int jhi, int klo, int khi);

template <typename T> void cpuSetGlobalBox(T *h, T *h_inv, T *boxlo_bound, 
        T *boxhi_bound, T *boxhi, T *boxlo, T *boxtilt, int triclinic);
template <typename T> void cpuSetLocalOrthBox(T *subhi, T *sublo, T *boxhi, 
        T *boxlo, T *subhi_lamda, T *sublo_lamda, int dim);
template <typename T> void cpuShiftedSubbox(T *ssublo, T *ssubhi, T *boxlo, T *boxhi, T *boxlo_lamda, 
        T *boxhi_lamda, T *sublo, T *subhi, T *sublo_lamda, T *subhi_lamda, T *epsilon, int *pbc, int triclinic);
template <typename T> void cpuBoundingSubbox(T *bsublo, T *bsubhi, T *sublo, T *subhi, 
        T *sublo_lamda, T *subhi_lamda, T *boxlo, T *h, int triclinic);

template <typename T> int cpuSetAtomType(T *x, T fraction, int *atomtype,
         int *seed, int *save, int seed0, int newtype, int dim, int nlocal);
template <typename T> void cpuAtomInside(int *inside, T *x, T *h_inv, T *boxlo, T *lo_lamda,
        T *hi_lamda, int dim, int n);
template <typename T> int cpuAtomAdd(T *y, int *atomtype, T *x, T *h_inv, T *boxlo, T *lo_lamda,
        T *hi_lamda, int *type, int dim, int nlocal, int n);
template <typename T> void cpuAtomLattice(T *y, int *atomtype, T *basis, T *primitive, T *rotaterow, T *origin, 
        T *latticespacing, T scale, int *basistype, int nbasis, int nlocal, 
        int ilo, int ihi, int jlo, int jhi, int klo, int khi, int dim);

#endif  


