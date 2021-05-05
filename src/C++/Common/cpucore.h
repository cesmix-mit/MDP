#ifndef __CPUCORE_H__
#define __CPUCORE_H__

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

void cpuQuadrupletnum(int *a, int *b, int n);
void cpuTripletnum(int *a, int *b, int n);
void cpuArrayFill(int *a, int start, int n);
template <typename T> T cpuArrayMaxError(T *a, T *b, int n);
template <typename T> T cpuArrayMin(T *a, int n);
template <typename T> T cpuArrayMax(T *a, int n);
template <typename T> T cpuArraySum(T *a, int n);
template <typename T> T cpuArraySquareSum(T *a, int n);
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

// neighbor list
template <typename T> void cpuMakeReferenceGrid(T *eta, T smin, T smax, T ds, int nc);
void cpuCell2AtomList(int *c2alist, int *c2anumsum, int *c2anum, int *clist, int natom, int ncell);

void cpuGridColoring2D(int *cellcolor, int *nc, int *bsize, int dim);
template <typename T> void cpuBoundingBox2D(T *vc, T*wc, T *v, T *w, T *a, T *b, T *r, int *pbc);
template <typename T> int cpuPeriodicImages2D(T *pimages, T *a, T *b, int *pbc);
template <typename T> void cpuAtomList2D(int *alist, int *inside, int *glistnumsum, int *glistnum, 
        T *x, T *pimages, T *wc, T *s2rmap, int inum, int pnum, int dim);
template <typename T> void cpuCellList2D(int *clist, T *x, T *eta1, T *eta2, T *eta3, T *s2rmap, 
        int *nc, int inum, int natom, int dim);
template <typename T> void cpuFullNeighborList2D(int *neighborlist, int *neighbornum, T *x, T *ellipsoid, int *alist, 
        int *clist, int *c2alist, int *c2anumsum, int *nc, int inum, int jnum, int dim);
template <typename T> void cpuFullNeighborList2D(int *neighborlist, int *neighbornum, T *x, T *rcutsq, int *alist, 
        int *clist, int *c2alist, int *c2anumsum, int *atomtype, int *nc, int ntype, int inum, int jnum, int dim);
template <typename T> void cpuHalfNeighborList2D(int *neighborlist, int *neighbornum, T *x, T *ellipsoid, int *alist, 
        int *clist, int *c2alist, int *c2anumsum, int *nc, int inum, int jnum, int dim);
template <typename T> void cpuHalfNeighborList2D(int *neighborlist, int *neighbornum, T *x, T *rcutsq, int *alist, 
        int *clist, int *c2alist, int *c2anumsum, int *atomtype, int *nc, int ntype, int inum, int jnum, int dim);        

void cpuGridColoring3D(int *cellcolor, int *nc, int *bsize, int dim);
template <typename T> void cpuBoundingBox3D(T *vc, T *wc, T *v, T *w, T *a, T *b, T *c, T *r, int *pbc);
template <typename T> int cpuPeriodicImages3D(T *pimages, T *a, T *b, T *c, int *pbc);
template <typename T> void cpuAtomList3D(int *alist, int *inside, int *glistnumsum, int *glistnum, 
        T *x, T *pimages, T *wc, T *s2rmap, int inum, int pnum, int dim);
template <typename T> void cpuCellList3D(int *clist, T *x, T *eta1, T *eta2, T *eta3, T *s2rmap, 
        int *nc, int inum, int natom, int dim);
template <typename T> void cpuFullNeighborList3D(int *neighborlist, int *neighbornum, T *x, T *ellipsoid, int *alist, 
        int *clist, int *c2alist, int *c2anumsum, int *nc, int inum, int jnum, int dim);
template <typename T> void cpuFullNeighborList3D(int *neighborlist, int *neighbornum, T *x, T *rcutsq, int *alist, 
        int *clist, int *c2alist, int *c2anumsum, int *atomtype, int *nc, int ntype, int inum, int jnum, int dim);
template <typename T> void cpuHalfNeighborList3D(int *neighborlist, int *neighbornum, T *x, T *ellipsoid, int *alist, 
        int *clist, int *c2alist, int *c2anumsum, int *nc, int inum, int jnum, int dim);
template <typename T> void cpuHalfNeighborList3D(int *neighborlist, int *neighbornum, T *x, T *rcutsq, int *alist, 
        int *clist, int *c2alist, int *c2anumsum, int *atomtype, int *nc, int ntype, int inum, int jnum, int dim);        

template <typename T> void cpuFullNeighborList2D(int *neighborlist, int *neighbornum, T *x, T *rcutsq, 
        int anum, int inum, int jnum, int dim);
template <typename T> void cpuFullNeighborList3D(int *neighborlist, int *neighbornum, T *x, T *rcutsq, 
        int anum, int inum, int jnum, int dim);

// template <typename T> void cpuGetNeighPairs(T *xij, T *qi, T *qj, T *x, T *q, int *ai, int *aj,  
//       int *ti, int *tj, int *anum, int *anumsum, int *ilist, int *alist, int *neighlist, int *neighnum, 
//       int *atomtype, int ncq, int inum, int jnum, int dim);
// template <typename T> void cpuGetHalfNeighPairs(T *xij, T *qi, T *qj, T *x, T *q, int *ai, int *aj,  
//       int *ti, int *tj, int *anum, int *anumsum, int *ilist, int *alist, int *neighlist, int *neighnum, 
//       int *atomtype, int ncq, int inum, int jnum, int dim);
// template <typename T> void cpuGetNeighPairs(T *xij, T *qi, T *qj, T *x, T *q, int *ai, int *aj,  
//       int *ti, int *tj, int *anum, int *anumsum, int *ilist, int *alist, int *neighlist, int *neighnum, 
//       int *atomtype, int typej, int ncq, int inum, int jnum, int dim);
// template <typename T> void cpuGetNeighPairs(T *xij, T *qi, T *qj, T *x, T *q, int *ai, int *aj, 
//       int *ti, int *tj, int *anum, int *anumsum, int *ilist, int *tlist, int *alist, int *neighlist, int *neighnum, 
//       int *atomtype, int *slist, int *start, int *index, int typei, int inum, int jnum, int ncq, int dim);
// template <typename T> void cpuGetNeighPairs(T *xij, T *qi, T *qj, T *x, T *q, int *ai, int *aj, 
//       int *ti, int *tj, int *anum, int *anumsum, int *ilist, int *tlist, int *alist, int *neighlist, int *neighnum, 
//       int *atomtype, int *slist, int *start, int *index, int typei, int typej, int inum, int jnum, int ncq, int dim);
// template <typename T> void cpuGetNeighTriplets(T *xij, T *xik, T *qi, T *qj, T *qk, T *x, T *q, int *ai, int *aj,  
//       int *ak, int *ti, int *tj, int *tk, int *anum, int *anumsum, int *ilist,  int *tlist, int *alist, int *neighlist, int *neighnum, 
//       int *atomtype, int *slist, int *start, int *index, int typei, int typej, int typek, int inum, int jnum, int ncq, int dim);
// template <typename T> void cpuGetNeighTriplets(T *xij, T *xik, T *qi, T *qj, T *qk, T *x, T *q, int *ai, int *aj,  
//       int *ak, int *ti, int *tj, int *tk, int *anum, int *anumsum, int *ilist, int *alist, int *neighlist, int *neighnum, 
//       int *atomtype, int inum, int jnum, int ncq, int dim);
// void cpuFullNeighPairList(int *pairnum, int *pairlist, int *ilist, int *neighlist, 
//         int *neighnum, int inum, int jnum);
// void cpuFullNeighPairList(int *pairnum, int *pairlist, int *atomtype, int *ilist, int *alist, 
//         int *neighlist, int *neighnum, int inum, int jnum, int typej);
// void cpuHalfNeighPairList(int *pairnum, int *pairlist, int *ilist, int *alist, 
//         int *neighlist, int *neighnum, int inum, int jnum);
// void cpuHalfNeighPairList(int *pairnum, int *pairlist, int *atomtype, int *ilist, int *alist, 
//         int *neighlist, int *neighnum, int inum, int jnum, int typej);
// void cpuHalfNeighTripletList(int *tripletnum, int *tripletlist, int *atomtype, int *ilist, int *alist, 
//         int *neighlist, int *neighnum, int inum, int jnum, int jknum, int typej, int typek);

template <typename T> void cpuNeighSingles(T *xi, T *qi, T *x, T *q, int *ai, 
      int *ti, int *ilist, int *atomtype, int inum, int ncq, int dim);
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
template <typename T> void cpuNeighPairs(T *xij, T *qi, T *qj, T *x, T *q, int *ai, int *aj,  
      int *ti, int *tj, int *pairnum, int *pairlist, int *pairnumsum, int *ilist, int *alist,  
      int *atomtype, int inum, int jnum, int ncq, int dim);
// template <typename T> void cpuNeighTripletList(int *tripletnum, int *tripletlist, T *x, T rcutsq, int *atomtype, int *ilist, int *alist, 
//         int *neighlist, int *neighnum, int inum, int jnum, int jknum, int typej, int typek, int dim);
// template <typename T> void cpuNeighTriplets(T *xij, T *xik, T *qi, T *qj, T *qk, T *x, T *q, int *ai, int *aj, int *ak,  
//       int *ti, int *tj, int *tk, int *tripletnum, int *tripletlist, int *tripletnumsum, int *ilist, int *alist,  
//       int *atomtype, int inum, int jnum, int jknum, int ncq, int dim);
// template <typename T> void cpuNeighQuadrupletList(int *quadrupletnum, int *quadrupletlist, T *x, T rcutsq, int *atomtype, int *ilist,
//         int *alist, int *neighlist, int *neighnum, int inum, int jnum, int jklnum, int typej, int typek, int typel, int dim);
// template <typename T> void cpuNeighQuadruplets(T *xij, T *xik, T *xil, T *qi, T *qj, T *qk, T *ql, T *x, T *q, 
//       int *ai, int *aj, int *ak, int *al, int *ti, int *tj, int *tk, int *tl, int *quadrupletnum, int *quadrupletlist, 
//       int *quadrupletnumsum, int *ilist, int *alist,  int *atomtype, int inum, int jnum, int jklnum, int ncq, int dim);
void cpuNeighTripletList(int *tripletlist, int *tripletnumsum, int *pairnum, 
        int *pairlist, int *ilist, int *alist, int inum, int jnum);
template <typename T> void cpuNeighTriplets(T *xij, T *xik, T *qi, T *qj, T *qk, T *x, T *q, int *ai, int *aj, int *ak,  
      int *ti, int *tj, int *tk, int *tripletnum, int *tripletlist, int *tripletnumsum, int *ilist, int *alist,  
      int *atomtype, int inum, int ncq, int dim);
void cpuNeighQuadrupletList(int *quadrupletlist, int *quadrupletnumsum, int *pairnum, 
        int *pairlist, int *ilist, int *alist, int inum, int jnum);
template <typename T> void cpuNeighQuadruplets(T *xij, T *xik, T *xil, T *qi, T *qj, T *qk, T *ql, T *x, T *q, 
      int *ai, int *aj, int *ak, int *al, int *ti, int *tj, int *tk, int *tl, int *quadrupletnum, int *quadrupletlist, 
      int *quadrupletnumsum, int *ilist, int *alist,  int *atomtype, int inum, int ncq, int dim);

template <typename T> void cpuNeighTripletList(int *tripletnum, int *tripletlist, T *x, T *rcutsq, int *pairnum, int *pairnumsum, 
        int *pairlist, int *atomtype, int *ilist, int *alist, int *neighlist, int *neighnum, int inum, 
        int jnum, int typek, int dim);
template <typename T> void cpuNeighTriplets(T *xij, T *xik, T *qi, T *qj, T *qk, T *x, T *q, int *ai, int *aj, int *ak,  
      int *ti, int *tj, int *tk, int *tripletnum, int *tripletlist, int *tripletnumsum, int *alist,  
      int *atomtype, int jnum, int ijnum, int ncq, int dim);

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

template <typename T> void cpuSingleDecomposition2D(T *f, T *fi, int *ai, int inum);
template <typename T> void cpuHalfForceDecomposition2D(T *fi, T *fij, int *ai, int *aj, int ijnum);
template <typename T> void cpuFullForceDecomposition2D(T *fi, T *fij, int *ai, int *aj, int ijnum);
template <typename T> void cpuIAtomDecomposition2D(T *fi, T *fij, int *ilist, int *anumsum, int inum);
template <typename T> void cpuJAtomDecomposition2D(T *fi, T *fij, int *jlist, int *bnumsum, int *index, int jnum);
template <typename T> void cpuForceDecompositionTriplet2D(T *fi, T *fij, T *fik, int *ai, int *aj, int *ak, int ijknum);
template <typename T> void cpuAtomDecompositionTriplet2D(T *fi, T *fij, T *fik, int *ilist, int *anumsum, int inum);
template <typename T> void cpuForceDecompositionQuadruplet2D(T *f, T *fij, T *fik,  T *fil, int *ai, int *aj, int *ak, int *al, int ijklnum);
template <typename T> void cpuAtomDecompositionQuadruplet2D(T *f, T *fij, T *fik, T *fil, int *ilist, int *anumsum, int inum);

template <typename T> void cpuSingleDecomposition3D(T *f, T *fi, int *ai, int inum);
template <typename T> void cpuHalfForceDecomposition3D(T *fi, T *fij, int *ai, int *aj, int ijnum);
template <typename T> void cpuFullForceDecomposition3D(T *fi, T *fij, int *ai, int *aj, int ijnum);
template <typename T> void cpuIAtomDecomposition3D(T *fi, T *fij, int *ilist, int *anumsum, int inum);
template <typename T> void cpuJAtomDecomposition3D(T *fi, T *fij, int *jlist, int *bnumsum, int *index, int jnum);
template <typename T> void cpuForceDecompositionTriplet3D(T *fi, T *fij, T *fik, int *ai, int *aj, int *ak, int ijknum);
template <typename T> void cpuAtomDecompositionTriplet3D(T *fi, T *fij, T *fik, int *ilist, int *anumsum, int inum);
template <typename T> void cpuForceDecompositionQuadruplet3D(T *f, T *fij, T *fik,  T *fil, int *ai, int *aj, int *ak, int *al, int ijklnum);
template <typename T> void cpuAtomDecompositionQuadruplet3D(T *f, T *fij, T *fik, T *fil, int *ilist, int *anumsum, int inum);

template <typename T> void cpuMakeReferenceGrid(T *eta, T smin, T smax, int nc, int pbc);
template <typename T> void cpuGhostAtoms(int *glistnum, T *x, T *pimages, T *wc, T *B2C, int n, int m, int dim);
template <typename T> void cpuCreateAtomList(int *ilist, int *glistnumsum, int *glistnum, int *atomtype, T *x, T *pimages, T *wc, T *B2C, int n, int m, int dim);
template <typename T> void cpuCellList(int *clist, int *c2anum, T *x, T *eta1, T *eta2, T *eta3, T *B2C, int *nc, int inum, int gnum, int dim);
void cpuCell2AtomList(int *c2alist, int *c2anumsum, int *c2anum, int *clist, int *nc, int inum, int gnum, int dim);
template <typename T> void cpuVerletAtoms(int *verletnum, T *x, T *ellipsoid, int *ilist, 
        int *clist, int *c2alist, int *c2anum, int *c2anumsum, int *nc, int inum, int dim);
template <typename T> void cpuCreateVerletList(int *verletlist, T *x, T *ellipsoid, int *verletnum, int *verletnumsum, 
     int *ilist, int *clist, int *c2alist, int *c2anum, int *c2anumsum, int *nc, int inum, int dim);
template <typename T> void cpuFullNeighNum(int *neighnum, T *x, T* ellipsoid, int *ilist, 
        int *verletlist, int *verletnum, int *verletnumsum, int inum, int dim);
template <typename T> void cpuFullNeighList(int *neighlist, T *x, T* ellipsoid, int *neighnum, int *neighnumsum,
        int *ilist, int *verletlist, int *verletnum, int *verletnumsum, int inum, int dim);
template <typename T> void cpuHalfNeighNum(int *neighnum, T *x, T* ellipsoid, 
        int *ilist, int *verletlist, int *verletnum, int *verletnumsum, int inum, int dim);
template <typename T> void cpuHalfNeighList(int *neighlist, T *x, T* ellipsoid, int *neighnum, int *neighnumsum,
        int *ilist, int *verletlist, int *verletnum, int *verletnumsum, int inum, int dim);
template <typename T> void cpuGetNeighPairs(T *xij, T *x, int *ti, int *tj, int *ilist, int *neighlist,  
        int *neighnum, int *neighnumsum, int *atomtype, int ntype, int inum, int dim);
template <typename T> void cpuGetNeighPairs(T *xij, T *x, int *anum, int *anumsum, int *ai, int *aj, int *ti, 
      int *tj, int *ilist, int *neighlist, int *neighnum, int *neighnumsum, int *atomtype, int inum, int dim);
template <typename T> int cpuGetNeighPairs(T *xij, T *x, int *anum, int *anumsum, int *ai, int *aj, int *ti, 
      int *tj, int *ilist, int *tlist, int *neighlist, int *neighnum, int *neighnumsum, int *atomtype, 
       int typei, int inum, int dim);
template <typename T> int cpuGetNeighPairs(T *xij, T *x, int *anum, int *anumsum, int *ai, int *aj, int *ti, 
      int *tj, int *ilist, int *tlist, int *neighlist, int *neighnum, int *neighnumsum, int *atomtype, 
       int typei, int typej, int inum, int dim);
template <typename T> void cpuGetNeighPairs(T *xi, T *xj, T *x, int *ti, int *tj, int *ilist, int *neighlist,  
        int *neighnum, int *neighnumsum, int *atomtype, int ntype, int inum, int dim);
template <typename T> void cpuVerletAtoms(int *verletnum, T *x, T *ellipsoid, int *atomtype, int *ilist, 
        int *clist, int *c2alist, int *c2anum, int *c2anumsum, int *nc, int ntype, int inum, int dim);
template <typename T> void cpuCreateVerletList(int *verletlist, T *x, T *ellipsoid, int *verletnum,  int *verletnumsum, 
     int *atomtype, int *ilist, int *clist, int *c2alist, int *c2anum, int *c2anumsum, int *nc, int ntype, int inum, int dim);
template <typename T> void cpuFullNeighNum(int *neighnum, T *x, T* ellipsoid, int *atomtype, 
        int *ilist, int *verletlist, int *verletnum, int *verletnumsum, int ntype, int inum, int dim);
template <typename T> void cpuFullNeighList(int *neighlist, T *x, T* ellipsoid, int *neighnum, int *neighnumsum,
        int *atomtype, int *ilist, int *verletlist, int *verletnum, int *verletnumsum, int ntype, int inum, int dim);
template <typename T> void cpuHalfNeighNum(int *neighnum, T *x, T* ellipsoid, int *atomtype, 
        int *ilist, int *verletlist, int *verletnum, int *verletnumsum, int ntype, int inum, int dim);
template <typename T> void cpuHalfNeighList(int *neighlist, T *x, T* ellipsoid, int *neighnum, int *neighnumsum,
        int *atomtype, int *ilist, int *verletlist, int *verletnum, int *verletnumsum, int ntype, int inum, int dim);

// spherical harmonic potential
void cpuGetIndk(int *indk, int K);
template <typename T> T clebschgordan(int j1, int m1, int j2, int m2, int j, int m, T *fac);
int cgcoefficients(int *indl, int N);
template <typename T> void cgcoefficients(T *cg, int *indm, int *rowm, int *indl, T *fac, int M, int N);
template <typename T> void cpuSphericalHarmonicsBispectrum(T *b, T *Sr, T *Si, T *fac, int L);
template <typename T> int cpuSphericalHarmonicsBispectrumIndex(T *b, int L);
template <typename T> void cpuSphericalHarmonicsBispectrumIndex(int *indl, T *b, int M, int L);
template <typename T> void cpuSphericalBessel(T *g, T *x, T *y, T *z, T *x0, T *f, int L, int K, int N);
template <typename T> void cpuSphericalHarmonics(T *Ylmr, T *Ylmi, T *x, T *y, T *z, 
                T *P, T *tmp, T *fac, T pi, int L, int N);
template <typename T> void cpuSphericalHarmonicsBessel(T *Sr, T *Si, T *x, T *y, T *z, 
                T *x0, T *P, T *tmp, T *f, T *fac, T pi, int L, int K, int N);
template <typename T> void cpuSphericalHarmonicsBessel(T *Sr, T *Si, T *xij,
                T *x0, T *P, T *tmp, T *f, T *fac, T pi, int L, int K, int N);
template <typename T> void cpuSphericalHarmonicsBesselJ(T *Sr, T *Si, T *xij,
                T *x0, T *P, T *tmp, T *f, T *fac, T pi, int L, int K, int N);
template <typename T> void cpuSphericalHarmonicsAtomSum(T *Ylmr, T *Ylmi, T *x, T *y, T *z, 
                T *P, T *tmp, T *fac, T pi, int L, int N);
template <typename T> void cpuSphericalHarmonicsSum(T *Sr, T *Si, T *Ylmr, T *Ylmi, int L, int N);
template <typename T> void cpuRadialSphericalHarmonicsSum(T *Sr, T *Si, T *Ylmr, T *Ylmi, T *g, int L, int K, int N);
template <typename T> void cpuSphericalHarmonicsBesselSum(T *Sr, T *Si, T *x, T *y, T *z, 
                T *x0, T *P, T *tmp, T *f, T *fac, T pi, int L, int K, int N);
template <typename T> void cpuSphericalHarmonicsBesselDeriv(T *Srx, T *Six, T *Sry, T *Siy, T *Srz, T *Siz, T *x, T *y, T *z, 
                T *x0, T *P, T *tmp, T *f, T *dP, T *dtmp, T *df, T *fac, T pi, int L, int K, int N);
template <typename T> void cpuSphericalHarmonicsBesselWithDeriv(T *Sr, T *Si, T *Srx, T *Six, T *Sry, T *Siy, T *Srz, T *Siz, T *x, T *y, T *z, 
                T *x0, T *P, T *tmp, T *f, T *dP, T *dtmp, T *df, T *fac, T pi, int L, int K, int N);
template <typename T> void cpuSphericalHarmonicsBesselWithDeriv(T *Sr, T *Si, T *Srx, T *Six, T *Sry, T *Siy, T *Srz, T *Siz, T *xij, 
                T *x0, T *P, T *tmp, T *f, T *dP, T *dtmp, T *df, T *fac, T pi, int L, int K, int N);
template <typename T> void cpuSphericalHarmonicsBesselJWithDeriv(T *Sr, T *Si, T *Srx, T *Six, T *Sry, T *Siy, T *Srz, T *Siz, T *xij, 
                T *x0, T *P, T *tmp, T *f, T *dP, T *dtmp, T *df, T *fac, T pi, int L, int K, int N);
template <typename T> void cpuRadialSphericalHarmonicsPower(T *p, T *ar, T *ai, int *indk, int L, int K);
template <typename T> void cpuRadialSphericalHarmonicsPowerDeriv(T *px, T *py, T *pz, T *ar, T *ai, 
        T *arx, T *aix, T *ary, T *aiy, T *arz, T *aiz, int *indk, int L, int K, int N);
template <typename T> void cpuRadialSphericalHarmonicsBispectrum(T *b, T *ar, T *ai, T *cg, int *indk, int *indl,
        int *indm, int *rowm, int Nub, int Ncg, int K);
template <typename T> void cpuRadialSphericalHarmonicsBispectrumDeriv(T *bx, T *by, T *bz, T *ar, T *ai, 
        T *arx, T *aix, T *ary, T *aiy, T *arz, T *aiz, T*cg, int *indk, int *indl,
        int *indm, int *rowm, int Nub, int Ncg, int K, int N);
template <typename T> void cpuRadialSphericalHarmonicsSum(T *ar, T *ai, T *Sr, T *Si, 
        int *Nnb, int Na, int L, int K);
template <typename T> void cpuRadialSphericalHarmonicsPower(T *p, T *ar, T *ai, int *indk, int Na, int L, int K);
template <typename T> void cpuRadialSphericalHarmonicsPowerDeriv(T *px, T *py, T *pz, T *ar, T *ai, 
        T *arx, T *aix, T *ary, T *aiy, T *arz, T *aiz, int *indk, int *Nnb, int Na, int L, int K);
template <typename T> void cpuRadialSphericalHarmonicsPowerDeriv2(T *pd, T *ar, T *ai, 
        T *arx, T *aix, T *ary, T *aiy, T *arz, T *aiz, int *indk, int *Nnb, int Na, int L, int K);
template <typename T> void cpuRadialSphericalHarmonicsBispectrum(T *b, T *ar, T *ai, T *cg, int *indk, 
        int *indl, int *indm, int *rowm, int Nub, int Ncg, int Na, int L, int K);
template <typename T> void cpuRadialSphericalHarmonicsBispectrumDeriv(T *bx, T *by, T *bz, T *ar, T *ai, 
        T *arx, T *aix, T *ary, T *aiy, T *arz, T *aiz, T*cg, int *indk, int *indl,
        int *indm, int *rowm, int *Nnb, int Na, int Nub, int Ncg, int K);
template <typename T> void cpuRadialSphericalHarmonicsBispectrumDeriv2(T *bd, T *ar, T *ai, 
        T *arx, T *aix, T *ary, T *aiy, T *arz, T *aiz, T*cg, int *indk, int *indl,
        int *indm, int *rowm, int *Nnb, int Na, int Nub, int Ncg, int K);
template <typename T> void cpuRadialSphericalHarmonicsBasis(T *d, T *c, int *atomtype, 
        int Ntype, int Na, int Nbf);
template <typename T> void cpuRadialSphericalHarmonicsBasisDeriv(T *dx, T *dy, T *dz, T *cx, T *cy, T *cz,
        int *atomtype, int *neighlist, int *Nnb, int Ntype, int Na, int Nbf);
template <typename T> void cpuRadialSphericalHarmonicsBasisDeriv2(T *dd, T *cd,
        int *atomtype, int *neighlist, int *Nnb, int Ntype, int Na, int Nbf);

template <typename T> void cpuRadialSphericalHarmonicsSpectrum(T *c, T *ar, T *ai, T *sr, T *si, T *cg, int *indk, 
        int *indl, int *indm, int *rowm, int *Nnb, int Nub, int Ncg, int Na, int L, int K, int spectrum);
template <typename T> void cpuRadialSphericalHarmonicsSpectrumDeriv(T *cd, T *ar, T *ai, T *sr, T *si, 
        T *srx, T *six, T *sry, T *siy, T *srz, T *siz, T *cg, int *indk, int *indl, int *indm, 
        int *rowm, int *Nnb, int Nub, int Ncg, int Na, int L, int K, int spectrum);
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



#endif  

