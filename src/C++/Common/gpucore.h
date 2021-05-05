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

// neighbor list
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

int gpuFindAtomType(int *tlist, int* ilist, int *atomtype, int *p, int *q, int typei, int inum);
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


// spherical harmonic potential
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

int gpuUniqueSort(int *b, int *c, int *d, int *e, int *a, int *p, int *t, int *q, int n);

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

//gpuSingleDecomposition(e, f, ei, fi, ai, na, dim); 

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

template <typename T> void gpuSingleDecomposition2D(T *f, T *fi, int *ai, int inum);
template <typename T> void gpuHalfForceDecomposition2D(T *fi, T *fij, int *ai, int *aj, int ijnum);
template <typename T> void gpuFullForceDecomposition2D(T *fi, T *fij, int *ai, int *aj, int ijnum);
template <typename T> void gpuIAtomDecomposition2D(T *fi, T *fij, int *ilist, int *anumsum, int inum);
template <typename T> void gpuJAtomDecomposition2D(T *fi, T *fij, int *jlist, int *bnumsum, int *index, int jnum);
template <typename T> void gpuForceDecompositionTriplet2D(T *fi, T *fij, T *fik, int *ai, int *aj, int *ak, int ijknum);
template <typename T> void gpuAtomDecompositionTriplet2D(T *fi, T *fij, T *fik, int *ilist, int *anumsum, int inum);
template <typename T> void gpuForceDecompositionQuadruplet2D(T *f, T *fij, T *fik,  T *fil, int *ai, int *aj, int *ak, int *al, int ijklnum);
template <typename T> void gpuAtomDecompositionQuadruplet2D(T *f, T *fij, T *fik, T *fil, int *ilist, int *anumsum, int inum);

template <typename T> void gpuSingleDecomposition3D(T *f, T *fi, int *ai, int inum);
template <typename T> void gpuHalfForceDecomposition3D(T *fi, T *fij, int *ai, int *aj, int ijnum);
template <typename T> void gpuFullForceDecomposition3D(T *fi, T *fij, int *ai, int *aj, int ijnum);
template <typename T> void gpuIAtomDecomposition3D(T *fi, T *fij, int *ilist, int *anumsum, int inum);
template <typename T> void gpuJAtomDecomposition3D(T *fi, T *fij, int *jlist, int *bnumsum, int *index, int jnum);
template <typename T> void gpuForceDecompositionTriplet3D(T *fi, T *fij, T *fik, int *ai, int *aj, int *ak, int ijknum);
template <typename T> void gpuAtomDecompositionTriplet3D(T *fi, T *fij, T *fik, int *ilist, int *anumsum, int inum);
template <typename T> void gpuForceDecompositionQuadruplet3D(T *f, T *fij, T *fik,  T *fil, int *ai, int *aj, int *ak, int *al, int ijklnum);
template <typename T> void gpuAtomDecompositionQuadruplet3D(T *f, T *fij, T *fik, T *fil, int *ilist, int *anumsum, int inum);

template <typename T> void gpuTripletForceDecomposition(T *f3ik, T *eij, T *d3ij, 
        int *tripletnum, int *tripletnumsum, int npairs, int dim);


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
template <typename T> void gpuLJ(T *u, T *f, T *xij, T *qi, T *qj, int *ti, int *tj, int *ai, int *aj, T *mu, T *eta, int *kappa, int dim, int ncq, int nmu, int neta, int nkappa, int ng);
template <typename T> void gpuGradientLJ(T *u, T *du, T *xij, T *u_x, T *qi, T *qj, int *ti, int *tj, 
        int *ai, int *aj, T *mu, T *eta, int *kappa, int dim, int ncq, 
                int nmu, int neta, int nkappa, int ng);

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

// template <typename T> void gpuGhostAtoms(int *glistnum, T *x, T *pimages, T *wc, T *B2C, int n, int m, int dim);
// template <typename T> void gpuCreateAtomList(int *ilist, int *glistnumsum, int *glistnum, int *atomtype, T *x, T *pimages, T *wc, T *B2C, int n, int m, int dim);
// template <typename T> void gpuCellList(int *clist, int *c2anum, T *x, T *eta1, T *eta2, T *eta3, T *B2C, int *nc, int inum, int gnum, int dim);
// void gpuCell2AtomList(int *c2alist, int *c2anumsum, int *c2anum, int *clist, int *nc, int inum, int gnum, int dim);
// template <typename T> void gpuVerletAtoms(int *verletnum, T *x, T *ellipsoid, int *ilist, 
//         int *clist, int *c2alist, int *c2anum, int *c2anumsum, int *nc, int inum, int dim);
// template <typename T> void gpuCreateVerletList(int *verletlist, T *x, T *ellipsoid, int *verletnum, int *verletnumsum, 
//      int *ilist, int *clist, int *c2alist, int *c2anum, int *c2anumsum, int *nc, int inum, int dim);
// template <typename T> void gpuFullNeighNum(int *neighnum, T *x, T* ellipsoid, int *ilist, 
//         int *verletlist, int *verletnum, int *verletnumsum, int inum, int dim);
// template <typename T> void gpuFullNeighList(int *neighlist, T *x, T* ellipsoid, int *neighnum, int *neighnumsum,
//         int *ilist, int *verletlist, int *verletnum, int *verletnumsum, int inum, int dim);
// template <typename T> void gpuHalfNeighNum(int *neighnum, T *x, T* ellipsoid, 
//         int *ilist, int *verletlist, int *verletnum, int *verletnumsum, int inum, int dim);
// template <typename T> void gpuHalfNeighList(int *neighlist, T *x, T* ellipsoid, int *neighnum, int *neighnumsum,
//         int *ilist, int *verletlist, int *verletnum, int *verletnumsum, int inum, int dim);
// template <typename T> void gpuGetNeighPairs(T *xij, T *x, int *ti, int *tj, int *ilist, int *neighlist,  
//         int *neighnum, int *neighnumsum, int *atomtype, int ntype, int inum, int dim);
// template <typename T> void gpuGetNeighPairs(T *xi, T *xj, T *x, int *ti, int *tj, int *ilist, int *neighlist,  
//         int *neighnum, int *neighnumsum, int *atomtype, int ntype, int inum, int dim);
// template <typename T> void gpuVerletAtoms(int *verletnum, T *x, T *ellipsoid, int *atomtype, int *ilist, 
//         int *clist, int *c2alist, int *c2anum, int *c2anumsum, int *nc, int ntype, int inum, int dim);
// template <typename T> void gpuCreateVerletList(int *verletlist, T *x, T *ellipsoid, int *verletnum,  int *verletnumsum, 
//      int *atomtype, int *ilist, int *clist, int *c2alist, int *c2anum, int *c2anumsum, int *nc, int ntype, int inum, int dim);
// template <typename T> void gpuFullNeighNum(int *neighnum, T *x, T* ellipsoid, int *atomtype, 
//         int *ilist, int *verletlist, int *verletnum, int *verletnumsum, int ntype, int inum, int dim);
// template <typename T> void gpuFullNeighList(int *neighlist, T *x, T* ellipsoid, int *neighnum, int *neighnumsum,
//         int *atomtype, int *ilist, int *verletlist, int *verletnum, int *verletnumsum, int ntype, int inum, int dim);
// template <typename T> void gpuHalfNeighNum(int *neighnum, T *x, T* ellipsoid, int *atomtype, 
//         int *ilist, int *verletlist, int *verletnum, int *verletnumsum, int ntype, int inum, int dim);
// template <typename T> void gpuHalfNeighList(int *neighlist, T *x, T* ellipsoid, int *neighnum, int *neighnumsum,
//         int *atomtype, int *ilist, int *verletlist, int *verletnum, int *verletnumsum, int ntype, int inum, int dim);
// 
// template <typename T> void gpuSphericalHarmonicsBesselDeriv(T *Srx, T *Six, T *Sry, T *Siy, T *Srz, T *Siz, 
//       T *x, T *y, T *z, T *x0, T *P, T *tmp, T *f, T *dP, T *dtmp, T *df, T *fac, T pi, int L, int K, int N);
// template <typename T> void gpuSphericalHarmonicsBesselWithDeriv(T *Sr, T *Si, T *Srx, T *Six, T *Sry, T *Siy, T *Srz, T *Siz, T *x, T *y, T *z, 
//                 T *x0, T *P, T *tmp, T *f, T *dP, T *dtmp, T *df, T *fac, T pi, int L, int K, int N);
// template <typename T> void gpuSphericalHarmonicsBesselWithDeriv(T *Sr, T *Si, T *Srx, T *Six, T *Sry, T *Siy, T *Srz, T *Siz, T *xij, 
//                 T *x0, T *P, T *tmp, T *f, T *dP, T *dtmp, T *df, T *fac, T pi, int L, int K, int N);
// template <typename T> void gpuRadialSphericalHarmonicsSum(T *ar, T *ai, T *Sr, T *Si, 
//         int *Nnb, int Na, int L, int K);
// template <typename T> void gpuRadialSphericalHarmonicsPower(T *p, T *ar, T *ai, int *indk, int Na, int L, int K);
// template <typename T> void gpuRadialSphericalHarmonicsPowerDeriv(T *px, T *py, T *pz, T *ar, T *ai, 
//         T *arx, T *aix, T *ary, T *aiy, T *arz, T *aiz, int *indk, int *Nnb, int Na, int L, int K);
// template <typename T> void gpuRadialSphericalHarmonicsPowerDeriv2(T *pd, T *ar, T *ai, 
//         T *arx, T *aix, T *ary, T *aiy, T *arz, T *aiz, int *indk, int *Nnb, int Na, int L, int K);
// template <typename T> void gpuRadialSphericalHarmonicsBispectrum(T *b, T *ar, T *ai, T *cg, int *indk, 
//         int *indl, int *indm, int *rowm, int Nub, int Ncg, int Na, int L, int K);
// template <typename T> void gpuRadialSphericalHarmonicsBispectrumDeriv(T *bx, T *by, T *bz, 
//         T *ar, T *ai, T *arx, T *aix, T *ary, T *aiy, T *arz, T *aiz, T*cg, int *indk, int *indl,
//         int *indm, int *rowm, int *Nnb, int Na, int Nub, int Ncg, int K);
// template <typename T> void gpuRadialSphericalHarmonicsBispectrumDeriv2(T *bd, T *ar, T *ai, 
//         T *arx, T *aix, T *ary, T *aiy, T *arz, T *aiz, T*cg, int *indk, int *indl,
//         int *indm, int *rowm, int *Nnb, int Na, int Nub, int Ncg, int K);
// template <typename T> void gpuRadialSphericalHarmonicsBasis(T *d, T *c, int *atomtype, 
//         int Ntype, int Na, int Nbf);
// template <typename T> void gpuRadialSphericalHarmonicsBasisDeriv(T *dx, T *dy, T *dz, 
//         T *cx, T *cy, T *cz, int *atomtype, int *neighlist, int *Nnb, int Ntype, int Na, int Nbf);
// template <typename T> void gpuRadialSphericalHarmonicsBasisDeriv2(T *dd, T *cd,
//         int *atomtype, int *neighlist, int *Nnb, int Ntype, int Na, int Nbf);
                
#endif  

