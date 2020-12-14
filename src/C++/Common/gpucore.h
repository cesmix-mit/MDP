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

template <typename T> void gpuGetElemNodes(T *un, T *u, int N, int nn, int np, int e1, int nc1, int nc);
template <typename T> void gpuPutElemNodes(T *un, T *u, int N, int nn, int np, int e1, int nc1, int nc);
template <typename T> void gpuGetFaceNodes0(T *uh, T *udg, int *facecon, int N, int ndf, int npf, int npe, int nc, int f1);
template <typename T> void gpuGetFaceNodes1(T *uh, T *udg, int *facecon, int N, int ndf, int npf, int npe, int nc, int f1);
template <typename T> void gpuGetFaceNodes2(T *uh, T *udg, int *facecon, int N, int ndf, int npf, int npe, int nc, int f1);
template <typename T> void gpuPutFaceNodes0(T *udg, T *uh, int *facecon, int N, int ndf, int npf, int npe, int nc, int f1);
template <typename T> void gpuPutFaceNodes1(T *udg, T *uh, int *facecon, int N, int ndf, int npf, int npe, int nc, int f1);
template <typename T> void gpuPutFaceNodes(T *udg, T *uh, int *rowe2f1, int *cole2f1, int *ent2ind1,
        int *rowe2f2, int *cole2f2, int *ent2ind2, int npf, int npe, int nc, int e1, int e2, int opts);

#endif  

