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

template <typename T> void cpuGetElemNodes(T *un, T *u, int np, int nc, int nc1, int nc2, int e1, int e2);
template <typename T> void cpuPutElemNodes(T *u, T *un, int np, int nc, int nc1, int nc2, int e1, int e2);
template <typename T> void cpuGetElemNodes2(T *un, T *u, int np, int nc, int nc1, int nc2, int e1, int e2);
template <typename T> void cpuPutElemNodes2(T *u, T *un, int np, int nc, int nc1, int nc2, int e1, int e2);
template <typename T> void cpuGetFaceNodes(T *uh, T *udg, int *facecon, int npf, int ncu, int npe, int nc, int f1, int f2, int opts);
template <typename T> void cpuPutFaceNodes(T *udg, T *uh, int *facecon, int npf, int ncu, int npe, int nc, int f1, int f2, int opts);
template <typename T> void cpuPutFaceNodes(T *udg, T *uh, int *rowe2f1, int *cole2f1, int *ent2ind1,
        int *rowe2f2, int *cole2f2, int *ent2ind2, int npf, int npe, int nc, int e1, int e2, int opts);

#endif  

