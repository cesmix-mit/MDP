#ifndef __GPUCORE_H__
#define __GPUCORE_H__

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

template <typename T> void gpuPrint2DArray(T* a, int m, int n);
template <typename T> void gpuPrint3DArray(T* a, int m, int n, int p);
template <typename T> void gpuArraySetValue(T *y, T a, int n);
template <typename T> void gpuArraySetValueAtIndex(T *y, T a, int n);
template <typename T> T gpuArrayGetValueAtIndex(T *y, int n);
template <typename T> void gpuArrayAddScalar(T *y, T a, int n);
template <typename T> void gpuArrayMultiplyScalar(T *y, T a, int n);

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

template <typename T> void gpuCart2Sphere(T *the, T *phi, T *r, T *x, T *y, T *z, int N);
template <typename T> void gpuCart2SphereDeriv(T *the, T *phi, T *r, T *thex, T *they, T *thez, T *phix, T *phiy, T *phiz, T *rx, T *ry, T *rz, T *x, T *y, T *z, int N);
template <typename T> void gpuSphere2Cart(T *x, T *y, T *z, T *the, T *phi, T *r, int N);
template <typename T> void gpuEuler2Rotm(T *R11, T *R12, T *R13, T *R21, 
                T *R22, T *R23, T *R31, T *R32, T *R33, T *alpha, T *beta, T *gamma, int N);
template <typename T> void gpuRotc(T *X, T *Y, T *Z, T *R, T *x, T *y, T *z, int N);

// neighbor list
template <typename T> void gpuGhostAtoms(int *glistnum, T *x, T *pimages, T *wc, T *B2C, int n, int m, int dim);
template <typename T> void gpuCreateAtomList(int *ilist, int *glistnumsum, int *glistnum, int *atomtype, T *x, T *pimages, T *wc, T *B2C, int n, int m, int dim);
template <typename T> void gpuCellList(int *clist, int *c2anum, T *x, T *eta1, T *eta2, T *eta3, T *B2C, int *nc, int inum, int gnum, int dim);
void gpuCell2AtomList(int *c2alist, int *c2anumsum, int *c2anum, int *clist, int *nc, int inum, int gnum, int dim);
template <typename T> void gpuVerletAtoms(int *verletnum, T *x, T *ellipsoid, int *ilist, 
        int *clist, int *c2alist, int *c2anum, int *c2anumsum, int *nc, int inum, int dim);
template <typename T> void gpuCreateVerletList(int *verletlist, T *x, T *ellipsoid, int *verletnum, int *verletnumsum, 
     int *ilist, int *clist, int *c2alist, int *c2anum, int *c2anumsum, int *nc, int inum, int dim);
template <typename T> void gpuFullNeighNum(int *neighnum, T *x, T* ellipsoid, int *ilist, 
        int *verletlist, int *verletnum, int *verletnumsum, int inum, int dim);
template <typename T> void gpuFullNeighList(int *neighlist, T *x, T* ellipsoid, int *neighnum, int *neighnumsum,
        int *ilist, int *verletlist, int *verletnum, int *verletnumsum, int inum, int dim);
template <typename T> void gpuHalfNeighNum(int *neighnum, T *x, T* ellipsoid, 
        int *ilist, int *verletlist, int *verletnum, int *verletnumsum, int inum, int dim);
template <typename T> void gpuHalfNeighList(int *neighlist, T *x, T* ellipsoid, int *neighnum, int *neighnumsum,
        int *ilist, int *verletlist, int *verletnum, int *verletnumsum, int inum, int dim);
template <typename T> void gpuGetNeighPairs(T *xij, T *x, int *ti, int *tj, int *ilist, int *neighlist,  
        int *neighnum, int *neighnumsum, int *atomtype, int ntype, int inum, int dim);
template <typename T> void gpuGetNeighPairs(T *xi, T *xj, T *x, int *ti, int *tj, int *ilist, int *neighlist,  
        int *neighnum, int *neighnumsum, int *atomtype, int ntype, int inum, int dim);
template <typename T> void gpuVerletAtoms(int *verletnum, T *x, T *ellipsoid, int *atomtype, int *ilist, 
        int *clist, int *c2alist, int *c2anum, int *c2anumsum, int *nc, int ntype, int inum, int dim);
template <typename T> void gpuCreateVerletList(int *verletlist, T *x, T *ellipsoid, int *verletnum,  int *verletnumsum, 
     int *atomtype, int *ilist, int *clist, int *c2alist, int *c2anum, int *c2anumsum, int *nc, int ntype, int inum, int dim);
template <typename T> void gpuFullNeighNum(int *neighnum, T *x, T* ellipsoid, int *atomtype, 
        int *ilist, int *verletlist, int *verletnum, int *verletnumsum, int ntype, int inum, int dim);
template <typename T> void gpuFullNeighList(int *neighlist, T *x, T* ellipsoid, int *neighnum, int *neighnumsum,
        int *atomtype, int *ilist, int *verletlist, int *verletnum, int *verletnumsum, int ntype, int inum, int dim);
template <typename T> void gpuHalfNeighNum(int *neighnum, T *x, T* ellipsoid, int *atomtype, 
        int *ilist, int *verletlist, int *verletnum, int *verletnumsum, int ntype, int inum, int dim);
template <typename T> void gpuHalfNeighList(int *neighlist, T *x, T* ellipsoid, int *neighnum, int *neighnumsum,
        int *atomtype, int *ilist, int *verletlist, int *verletnum, int *verletnumsum, int ntype, int inum, int dim);

// spherical harmonic potential
template <typename T> void gpuSphericalHarmonicsBessel(T *Sr, T *Si, T *x, T *y, T *z, 
                T *x0, T *P, T *tmp, T *f, T *fac, T pi, int L, int K, int N);
template <typename T> void gpuSphericalHarmonicsBessel(T *Sr, T *Si, T *xij,
                T *x0, T *P, T *tmp, T *f, T *fac, T pi, int L, int K, int N);
template <typename T> void gpuSphericalHarmonicsBesselDeriv(T *Srx, T *Six, T *Sry, T *Siy, T *Srz, T *Siz, 
      T *x, T *y, T *z, T *x0, T *P, T *tmp, T *f, T *dP, T *dtmp, T *df, T *fac, T pi, int L, int K, int N);
template <typename T> void gpuSphericalHarmonicsBesselWithDeriv(T *Sr, T *Si, T *Srx, T *Six, T *Sry, T *Siy, T *Srz, T *Siz, T *x, T *y, T *z, 
                T *x0, T *P, T *tmp, T *f, T *dP, T *dtmp, T *df, T *fac, T pi, int L, int K, int N);
template <typename T> void cpuSphericalHarmonicsBesselWithDeriv(T *Sr, T *Si, T *Srx, T *Six, T *Sry, T *Siy, T *Srz, T *Siz, T *xij, 
                T *x0, T *P, T *tmp, T *f, T *dP, T *dtmp, T *df, T *fac, T pi, int L, int K, int N);
template <typename T> void gpuRadialSphericalHarmonicsSum(T *ar, T *ai, T *Sr, T *Si, 
        int *Nnb, int Na, int L, int K);
template <typename T> void gpuRadialSphericalHarmonicsPower(T *p, T *ar, T *ai, int *indk, int Na, int L, int K);
template <typename T> void gpuRadialSphericalHarmonicsPowerDeriv(T *px, T *py, T *pz, T *ar, T *ai, 
        T *arx, T *aix, T *ary, T *aiy, T *arz, T *aiz, int *indk, int *Nnb, int Na, int L, int K);
template <typename T> void gpuRadialSphericalHarmonicsPowerDeriv2(T *pd, T *ar, T *ai, 
        T *arx, T *aix, T *ary, T *aiy, T *arz, T *aiz, int *indk, int *Nnb, int Na, int L, int K);
template <typename T> void gpuRadialSphericalHarmonicsBispectrum(T *b, T *ar, T *ai, T *cg, int *indk, 
        int *indl, int *indm, int *rowm, int Nub, int Ncg, int Na, int L, int K);
template <typename T> void gpuRadialSphericalHarmonicsBispectrumDeriv(T *bx, T *by, T *bz, 
        T *ar, T *ai, T *arx, T *aix, T *ary, T *aiy, T *arz, T *aiz, T*cg, int *indk, int *indl,
        int *indm, int *rowm, int *Nnb, int Na, int Nub, int Ncg, int K);
template <typename T> void gpuRadialSphericalHarmonicsBispectrumDeriv2(T *bd, T *ar, T *ai, 
        T *arx, T *aix, T *ary, T *aiy, T *arz, T *aiz, T*cg, int *indk, int *indl,
        int *indm, int *rowm, int *Nnb, int Na, int Nub, int Ncg, int K);
template <typename T> void gpuRadialSphericalHarmonicsBasis(T *d, T *c, int *atomtype, 
        int Ntype, int Na, int Nbf);
template <typename T> void gpuRadialSphericalHarmonicsBasisDeriv(T *dx, T *dy, T *dz, 
        T *cx, T *cy, T *cz, int *atomtype, int *neighlist, int *Nnb, int Ntype, int Na, int Nbf);
template <typename T> void gpuRadialSphericalHarmonicsBasisDeriv2(T *dd, T *cd,
        int *atomtype, int *neighlist, int *Nnb, int Ntype, int Na, int Nbf);
                
#endif  

