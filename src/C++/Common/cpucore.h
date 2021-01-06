#ifndef __CPUCORE_H__
#define __CPUCORE_H__

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
template <typename T> void cpuSphericalHarmonicsSum(T *Ylmr, T *Ylmi, T *x, T *y, T *z, 
                T *P, T *tmp, T *fac, T pi, int L, int N);
template <typename T> void cpuSphericalHarmonicsSum(T *Sr, T *Si, T *Ylmr, T *Ylmi, int L, int N);
template <typename T> void cpuRadialSphericalHarmonicsSum(T *Sr, T *Si, T *Ylmr, T *Ylmi, T *g, int L, int K, int N);
template <typename T> void cpuSphericalHarmonicsBesselSum(T *Sr, T *Si, T *x, T *y, T *z, 
                T *x0, T *P, T *tmp, T *f, T *fac, T pi, int L, int K, int N);
template <typename T> void cpuSphericalHarmonicsBesselDeriv(T *Srx, T *Six, T *Sry, T *Siy, T *Srz, T *Siz, T *x, T *y, T *z, 
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
template <typename T> void cpuRadialSphericalHarmonicsBispectrum(T *b, T *ar, T *ai, T *cg, int *indk, 
        int *indl, int *indm, int *rowm, int Nub, int Ncg, int Na, int L, int K);
template <typename T> void cpuRadialSphericalHarmonicsBispectrumDeriv(T *bx, T *by, T *bz, T *ar, T *ai, 
        T *arx, T *aix, T *ary, T *aiy, T *arz, T *aiz, T*cg, int *indk, int *indl,
        int *indm, int *rowm, int *Nnb, int Na, int Nub, int Ncg, int K);
template <typename T> void cpuRadialSphericalHarmonicsBasis(T *d, T *c, int *atomtype, 
        int Ntype, int Na, int Nbf);
template <typename T> void cpuRadialSphericalHarmonicsBasisDeriv(T *dx, T *dy, T *dz, T *cx, T *cy, T *cz,
        int *atomtype, int *neighlist, int *Nnb, int Ntype, int Na, int Nbf);


#endif  

