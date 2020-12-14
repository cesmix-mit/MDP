#ifndef __OPUCORE_H__
#define __OPUCORE_H__

template <typename T> void opuKron(T *C, T *A, T *B, int M1, int N1, int M2, int N2);
void opuIndexPermute12(int *index, int I1, int I2, int I3);
void opuIndexPermute13(int *index, int I1, int I2, int I3, int I4);
void opuIndexPermute23(int *index, int I1, int I2, int I3, int I4);
template <typename T> void opuPermute(T *B, T *A, int *index, int N);
template <typename T> void opuPermuteSharedMem(T *B, T *A, int *index, int N, int BLOCKDIM);
template <typename T> void opuPermute12(T *B, T *A, int I1, int I2, int I3);
template <typename T> void opuPermute13(T *B, T *A, int I1, int I2, int I3, int I4);
template <typename T> void opuPermute23(T *B, T *A, int I1, int I2, int I3, int I4);

template <typename T> void opuGetArrayAtIndex(T *y, T *x, int *ind, int n);
template <typename T> void opuPutArrayAtIndex(T *y, T *x, int *ind, int n);
template <typename T> void opuArrayPlusXAtIndex(T *y, T *x, int *ind, int n);
template <typename T> void opuArrayMinusXAtIndex(T *y, T *x, int *ind, int n);
template <typename T> void opuArrayAXPYAtIndex(T *y, T *x, T a, int *ind, int n);
template <typename T> void opuArrayAverageAtIndex(T *y, T *x, int *ind1, int *ind2, int n);

template <typename T> T opuArrayMin(T *a, int n);
template <typename T> T opuArrayMax(T *a, int n);
template <typename T> T opuArraySum(T *a, int n);
template <typename T> T opuArraySquareSum(T *a, int n);
template <typename T> T opuArrayMean(T *a, int n);
template <typename T> void opuArraySetValue(T *y, T a, int n);
template <typename T> void opuArrayAddScalar(T *y, T a, int n);
template <typename T> void opuArrayMultiplyScalar(T *y, T a, int n);

template <typename T> void opuArrayCopy(T *y, T *x, int n);
template <typename T> void opuArrayMinus(T *y, T *x, int n);
template <typename T> void opuArrayAbs(T *y, T *x, int n);
template <typename T> void opuArraySqrt(T *y, T *x, int n);
template <typename T> void opuArraySin(T *y, T *x, int n);
template <typename T> void opuArrayCos(T *y, T *x, int n);
template <typename T> void opuArrayTan(T *y, T *x, int n);
template <typename T> void opuArrayAsin(T *y, T *x, int n);
template <typename T> void opuArrayAcos(T *y, T *x, int n);
template <typename T> void opuArrayAtan(T *y, T *x, int n);
template <typename T> void opuArraySinh(T *y, T *x, int n);
template <typename T> void opuArrayCosh(T *y, T *x, int n);
template <typename T> void opuArrayTanh(T *y, T *x, int n);
template <typename T> void opuArrayAsinh(T *y, T *x, int n);
template <typename T> void opuArrayAcosh(T *y, T *x, int n);
template <typename T> void opuArrayAtanh(T *y, T *x, int n);
template <typename T> void opuArrayExp(T *y, T *x, int n);
template <typename T> void opuArrayLog(T *y, T *x, int n);
template <typename T> void opuArrayCeil(T *y, T *x, int n);
template <typename T> void opuArrayFloor(T *y, T *x, int n);
template <typename T> void opuArrayErf(T *y, T *x, int n);
template <typename T> void opuArrayErfc(T *y, T *x, int n);
template <typename T> void opuArraySquare(T *y, T *x, int n);
template <typename T> void opuArrayPower(T *y, T *x, int p, int n);

template <typename T> void opuArrayMultiplyScalarDiagonal(T *C, T a, int n);
template <typename T> void opuArrayAddVectorToDiagonal(T *C, T *x, T a, int n);
template <typename T> void opuArrayRowAverage(T *a, T *b, int M, int N);
template <typename T> void opuArrayAXPB(T *y, T *x, T a, T b, int n);
template <typename T> void opuArrayAXPBY(T *z, T *x, T *y, T a, T b, int n);
template <typename T> void opuArrayAXY(T *s, T *x, T *y, T a, int n);
template <typename T> void opuArrayAXYZ(T *s, T *x, T *y, T *z, T a, int n);
template <typename T> void opuArrayAXYPBZ(T *s, T *x, T *y, T *z, T a, T b, int n);
template <typename T> void opuArrayAdd3Vectors(T *s, T *x, T *y, T *z, T a, T b, T c, int n);
template <typename T> void opuArrayExtract(T *un, T *u, int I, int J, int K, 
        int i1, int i2, int j1, int j2, int k1, int k2);
template <typename T> void opuArrayInsert(T *u, T *un, int I, int J, int K, 
        int i1, int i2, int j1, int j2, int k1, int k2);
template <typename T> void opuArrayGemmBatch(T *C, T *A, T *B, int I, int J, int K, int S);
template <typename T> void opuArrayGemmBatch1(T *C, T *A, T *B, int I, int J, int K, int S);
template <typename T> void opuArrayDG2CG(T *ucg, T *udg, int *cgent2dgent, int *rowent2elem, int nent);
template <typename T> void opuArrayDG2CG2(T *ucg, T *udg, int *colent2elem, int *rowent2elem, int nent, int npe);

template <typename T> void opuGetElemNodes(T *un, T *u, int np, int nc, int nc1, int nc2, int e1, int e2);
template <typename T> void opuPutElemNodes(T *u, T *un, int np, int nc, int nc1, int nc2, int e1, int e2);
template <typename T> void opuGetElemNodes2(T *un, T *u, int np, int nc, int nc1, int nc2, int e1, int e2);
template <typename T> void opuPutElemNodes2(T *u, T *un, int np, int nc, int nc1, int nc2, int e1, int e2);
template <typename T> void opuGetFaceNodes(T *uh, T *udg, int *facecon, int npf, int ncu, int npe, int nc, int f1, int f2, int opts);
template <typename T> void opuPutFaceNodes(T *udg, T *uh, int *facecon, int npf, int ncu, int npe, int nc, int f1, int f2, int opts);
template <typename T> void opuPutFaceNodes(T *udg, T *uh, int *rowe2f1, int *cole2f1, int *ent2ind1,
        int *rowe2f2, int *cole2f2, int *ent2ind2, int npf, int npe, int nc, int e1, int e2, int opts);


#endif  

