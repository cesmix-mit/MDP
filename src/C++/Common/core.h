#ifndef __CORE_H__
#define __CORE_H__

#include "cpucore.h"

#ifdef HAVE_CUDA      
#include "gpucore.h"
#endif                

static void GetArrayAtIndex(dstype *y, dstype *x, Int *ind, Int n, Int backend)
{
    if (backend == 1)  
        cpuGetArrayAtIndex(y, x, ind, n);
#ifdef HAVE_CUDA            
    if (backend == 2)  
        gpuGetArrayAtIndex(y, x, ind, n);
#endif                  
}

static void PutArrayAtIndex(dstype *y, dstype *x, Int *ind, Int n, Int backend)
{
    if (backend == 1)  
        cpuPutArrayAtIndex(y, x, ind, n);
#ifdef HAVE_CUDA            
    if (backend == 2)  
        gpuPutArrayAtIndex(y, x, ind, n);
#endif                  
}

static void ArrayPlusXAtIndex(dstype *y, dstype *x, Int *ind, Int n, Int backend)
{
    if (backend == 1)  
        cpuArrayPlusXAtIndex(y, x, ind, n);
#ifdef HAVE_CUDA            
    if (backend == 2)  
        gpuArrayPlusXAtIndex(y, x, ind, n);
#endif                  
}

static void ArrayMinusXAtIndex(dstype *y, dstype *x, Int *ind, Int n, Int backend)
{
    if (backend == 1)  
        cpuArrayMinusXAtIndex(y, x, ind, n);
#ifdef HAVE_CUDA            
    if (backend == 2)  
        gpuArrayMinusXAtIndex(y, x, ind, n);
#endif                  
}

static void ArraySquare(dstype *y, dstype *x, Int n, Int backend)
{
    if (backend == 1)  
        cpuArraySquare(y, x, n);   
#ifdef HAVE_CUDA            
    if (backend == 2)  
        gpuArraySquare(y, x, n);
#endif                  
}

static void ArraySetValue(dstype *y, dstype a, Int n, Int backend)
{
    if (backend == 1)  
        cpuArraySetValue(y, a, n);
#ifdef HAVE_CUDA            
    if (backend == 2)  
        gpuArraySetValue(y, a, n);
#endif                  
}

static void ArrayAddScalar(dstype *y, dstype a, Int n, Int backend)
{
    if (backend == 1)  
        cpuArrayAddScalar(y, a, n);
#ifdef HAVE_CUDA            
    if (backend == 2)  
        gpuArrayAddScalar(y, a, n);
#endif                  
}

static void ArrayMultiplyScalar(dstype *y, dstype a, Int n, Int backend)
{
    if (backend == 1)  
        cpuArrayMultiplyScalar(y, a, n);
#ifdef HAVE_CUDA            
    if (backend == 2)  
        gpuArrayMultiplyScalar(y, a, n);
#endif                  
}

static void ArrayCopy(dstype *y, dstype *x, Int n, Int backend)
{
    if (backend == 1)  
        cpuArrayCopy(y, x, n);   
#ifdef HAVE_CUDA            
    if (backend == 2)  
        gpuArrayCopy(y, x, n);
#endif                  
}

static void ArrayMinus(dstype *y, dstype *x, Int n, Int backend)
{
    if (backend == 1)  
        cpuArrayMinus(y, x, n);   
#ifdef HAVE_CUDA            
    if (backend == 2)  
        gpuArrayMinus(y, x, n);
#endif                  
}

static void ArrayAbs(dstype *y, dstype *x, Int n, Int backend)
{
    if (backend == 1)  
        cpuArrayAbs(y, x, n);   
#ifdef HAVE_CUDA            
    if (backend == 2)  
        gpuArrayAbs(y, x, n);
#endif                  
}

static void ArraySqrt(dstype *y, dstype *x, Int n, Int backend)
{
    if (backend == 1)  
        cpuArraySqrt(y, x, n);   
#ifdef HAVE_CUDA            
    if (backend == 2)  
        gpuArraySqrt(y, x, n);
#endif                  
}

static void ArrayMultiplyScalarDiagonal(dstype *C, dstype a, Int n, Int backend)
{
    if (backend == 1)  
        cpuArrayMultiplyScalarDiagonal(C, a, n);   
#ifdef HAVE_CUDA            
    if (backend == 2)  
        gpuArrayMultiplyScalarDiagonal(C, a, n);
#endif                  
}

static void ArrayAddVectorToDiagonal(dstype *C, dstype *x, dstype a, Int n, Int backend)
{
    if (backend == 1)  
        cpuArrayAddVectorToDiagonal(C, x, a, n);   
#ifdef HAVE_CUDA            
    if (backend == 2)  
        gpuArrayAddVectorToDiagonal(C, x, a, n);
#endif                  
}

static void ArrayAXPB(dstype *y, dstype *x, dstype a, dstype b, Int n, Int backend)
{
    if (backend == 1)  
        cpuArrayAXPB(y, x, a, b, n);
#ifdef HAVE_CUDA            
    if (backend == 2)  
        gpuArrayAXPB(y, x, a, b, n);
#endif                  
}

static void ArrayAXPBY(dstype *z, dstype *x, dstype *y, dstype a, dstype b, Int n, Int backend)
{
    if (backend == 1)  
        cpuArrayAXPBY(z, x, y, a, b, n);
#ifdef HAVE_CUDA            
    if (backend == 2)  
        gpuArrayAXPBY(z, x, y, a, b, n);
#endif                  
}

static void ArrayAXY(dstype *z, dstype *x, dstype *y, dstype a, Int n, Int backend)
{
    if (backend == 1)  
        cpuArrayAXY(z, x, y, a, n);
#ifdef HAVE_CUDA            
    if (backend == 2)  
        gpuArrayAXY(z, x, y, a, n);
#endif                  
}

static void ArrayAXYZ(dstype *s, dstype *x, dstype *y, dstype *z, dstype a, Int n, Int backend)
{
    if (backend == 1)  
        cpuArrayAXYZ(s, x, y, z, a, n);
#ifdef HAVE_CUDA            
    if (backend == 2)  
        gpuArrayAXYZ(s, x, y, z, a, n);
#endif                  
}

static void ArrayAXYPBZ(dstype *s, dstype *x, dstype *y, dstype *z, dstype a, dstype b, Int n, Int backend)
{
    if (backend == 1)  
        cpuArrayAXYPBZ(s, x, y, z, a, b, n);
#ifdef HAVE_CUDA            
    if (backend == 2)  
        gpuArrayAXYPBZ(s, x, y, z, a, b, n);
#endif                  
}

static void ArrayAdd3Vectors(dstype *s, dstype *x, dstype *y, dstype *z, dstype a, dstype b, dstype c, Int n, Int backend)
{
    if (backend == 1)  
        cpuArrayAdd3Vectors(s, x, y, z, a, b, c, n);
#ifdef HAVE_CUDA            
    if (backend == 2)  
        gpuArrayAdd3Vectors(s, x, y, z, a, b, c, n);
#endif                  
}

static void ArrayExtract(dstype *un, dstype *u, Int I, Int J, Int K, 
        Int i1, Int i2, Int j1, Int j2, Int k1, Int k2, Int backend)
{
    if (backend == 1)  
        cpuArrayExtract(un, u, I, J, K, i1, i2, j1, j2, k1, k2);
#ifdef HAVE_CUDA            
    if (backend == 2)  
        gpuArrayExtract(un, u, I, J, K, i1, i2, j1, j2, k1, k2);
#endif                  
}

static void ArrayInsert(dstype *u, dstype *un, Int I, Int J, Int K, 
        Int i1, Int i2, Int j1, Int j2, Int k1, Int k2, Int backend)
{
    if (backend == 1)  
        cpuArrayInsert(u, un, I, J, K, i1, i2, j1, j2, k1, k2);
#ifdef HAVE_CUDA            
    if (backend == 2)  
        gpuArrayInsert(u, un, I, J, K, i1, i2, j1, j2, k1, k2);    
#endif                  
}

static void ArrayGemmBatch(dstype *C, dstype *A, dstype *B, Int I, Int J, Int K, Int S, Int backend)
{
    if (backend == 1)  
        cpuArrayGemmBatch(C, A, B, I, J, K, S);    
#ifdef HAVE_CUDA            
    if (backend == 2)  
        gpuArrayGemmBatch(C, A, B, I, J, K, S);    
#endif                  
}

static void ArrayGemmBatch1(dstype *C, dstype *A, dstype *B, Int I, Int J, Int K, Int S, Int backend)
{
    if (backend == 1)  
        cpuArrayGemmBatch1(C, A, B, I, J, K, S);    
#ifdef HAVE_CUDA            
    if (backend == 2)  
        gpuArrayGemmBatch1(C, A, B, I, J, K, S);    
#endif                  
}

static void ArrayDG2CG(dstype *ucg, dstype *udg, Int *cgent2dgent, Int *rowent2elem, Int nent, Int backend)
{
    if (backend == 1)  
        cpuArrayDG2CG(ucg, udg, cgent2dgent, rowent2elem, nent);    
#ifdef HAVE_CUDA            
    if (backend == 2)  
        gpuArrayDG2CG(ucg, udg, cgent2dgent, rowent2elem, nent);    
#endif                  
}

static void ArrayDG2CG2(dstype *ucg, dstype *udg, Int *colent2elem, Int *rowent2elem, Int nent, Int npe, Int backend)
{
    if (backend == 1)  
        cpuArrayDG2CG2(ucg, udg, colent2elem, rowent2elem, nent, npe);    
#ifdef HAVE_CUDA            
    if (backend == 2)  
        gpuArrayDG2CG2(ucg, udg, colent2elem, rowent2elem, nent, npe);    
#endif                  
}

static void Cart2Sphere(dstype *the, dstype *phi, dstype *r, dstype *x, dstype *y, dstype *z, Int N, Int backend)
{
    if (backend == 1)  
        cpuCart2Sphere(the, phi, r, x, y, z, N);    
#ifdef HAVE_CUDA            
    if (backend == 2)  
        gpuCart2Sphere(the, phi, r, x, y, z, N);    
#endif                  
}

static void Cart2SphereDeriv(dstype *the, dstype *phi, dstype *r, dstype *thex, dstype *they, dstype *thez, dstype *phix, 
        dstype *phiy, dstype *phiz, dstype *rx, dstype *ry, dstype *rz, dstype *x, dstype *y, dstype *z, Int N, Int backend)
{
    if (backend == 1)  
        cpuCart2SphereDeriv(the, phi, r, thex, they, thez, phix, phiy, phiz,
            rx, ry, rz, x, y, z, N);
#ifdef HAVE_CUDA            
    if (backend == 2)  
        gpuCart2SphereDeriv(the, phi, r, thex, they, thez, phix, phiy, phiz,
            rx, ry, rz, x, y, z, N);
#endif                  
}

static void Sphere2Cart(dstype *x, dstype *y, dstype *z, dstype *the, dstype *phi, dstype *r, Int N, Int backend)
{
    if (backend == 1)  
        cpuSphere2Cart(x, y, z, the, phi, r, N);    
#ifdef HAVE_CUDA            
    if (backend == 2)  
        gpuSphere2Cart(x, y, z, the, phi, r, N);    
#endif                  
}

static void Euler2Rotm(dstype *R11, dstype *R12, dstype *R13, dstype *R21, 
                dstype *R22, dstype *R23, dstype *R31, dstype *R32, dstype *R33, dstype *alpha, dstype *beta, dstype *gamma, Int N, Int backend)
{
    if (backend == 1)  
        cpuEuler2Rotm(R11, R12, R13, R21, R22, R23, R31, R32, R33, alpha, beta, gamma, N);    
#ifdef HAVE_CUDA            
    if (backend == 2)  
        gpuEuler2Rotm(R11, R12, R13, R21, R22, R23, R31, R32, R33, alpha, beta, gamma, N);    
#endif                  
}

static void Rotc(dstype *X, dstype *Y, dstype *Z, dstype *R, dstype *x, dstype *y, dstype *z, Int N, Int backend)
{
    if (backend == 1)  
        cpuRotc(X, Y, Z, R, x, y, z, N);    
#ifdef HAVE_CUDA            
    if (backend == 2)  
        gpuRotc(X, Y, Z, R, x, y, z, N);    
#endif                  
}

static void coreSphericalHarmonicsBessel(dstype *Sr, dstype *Si, dstype *x, dstype *y, dstype *z, 
                dstype *x0, dstype *P, dstype *tmp, dstype *f, dstype *fac, dstype pi, Int L, Int K, Int N, Int backend)
{
    if (backend == 1)  
        cpuSphericalHarmonicsBessel(Sr, Si, x, y, z, x0, P, tmp, f, fac, pi, L, K, N);    
#ifdef HAVE_CUDA            
    if (backend == 2)  
        gpuSphericalHarmonicsBessel(Sr, Si, x, y, z, x0, P, tmp, f, fac, pi, L, K, N);    
#endif                  
}

static void coreSphericalHarmonicsBesselDeriv(dstype *Srx, dstype *Six, dstype *Sry, dstype *Siy, dstype *Srz, dstype *Siz, 
      dstype *x, dstype *y, dstype *z, dstype *x0, dstype *P, dstype *tmp, dstype *f, dstype *dP, dstype *dtmp, dstype *df, dstype *fac, dstype pi, Int L, Int K, Int N, Int backend)
{
    if (backend == 1)  
        cpuSphericalHarmonicsBesselDeriv(Srx, Six, Sry, Siy, Srz, Siz, x, y, z, 
            x0, P, tmp, f, dP, dtmp, df, fac, pi, L, K, N);    
#ifdef HAVE_CUDA            
    if (backend == 2)  
        gpuSphericalHarmonicsBesselDeriv(Srx, Six, Sry, Siy, Srz, Siz, x, y, z, 
            x0, P, tmp, f, dP, dtmp, df, fac, pi, L, K, N);    
#endif                  
}

static void coreSphericalHarmonicsBesselWithDeriv(dstype *Sr, dstype *Si, dstype *Srx, dstype *Six, dstype *Sry, dstype *Siy, dstype *Srz, dstype *Siz, 
      dstype *x, dstype *y, dstype *z, dstype *x0, dstype *P, dstype *tmp, dstype *f, dstype *dP, dstype *dtmp, dstype *df, dstype *fac, dstype pi, Int L, Int K, Int N, Int backend)
{
    if (backend == 1)  
        cpuSphericalHarmonicsBesselWithDeriv(Sr, Si, Srx, Six, Sry, Siy, Srz, Siz, x, y, z, 
            x0, P, tmp, f, dP, dtmp, df, fac, pi, L, K, N);    
#ifdef HAVE_CUDA            
    if (backend == 2)  
        gpuSphericalHarmonicsBesselWithDeriv(Sr, Si, Srx, Six, Sry, Siy, Srz, Siz, x, y, z, 
            x0, P, tmp, f, dP, dtmp, df, fac, pi, L, K, N);    
#endif                  
}

static void coreRadialSphericalHarmonicsSum(dstype *ar, dstype *ai, dstype *Sr, dstype *Si, Int *Nnb, Int Na, Int L, Int K, Int backend)
{
    if (backend == 1)  
        cpuRadialSphericalHarmonicsSum(ar, ai, Sr, Si, Nnb, Na, L, K);    
#ifdef HAVE_CUDA            
    if (backend == 2)  
        gpuRadialSphericalHarmonicsSum(ar, ai, Sr, Si, Nnb, Na, L, K);    
#endif                  
}

static void coreRadialSphericalHarmonicsPower(dstype *p, dstype *ar, dstype *ai, Int *indk, Int Na, Int L, Int K, Int backend)
{
    if (backend == 1)  
        cpuRadialSphericalHarmonicsPower(p, ar, ai, indk, Na, L, K);    
#ifdef HAVE_CUDA            
    if (backend == 2)  
        gpuRadialSphericalHarmonicsPower(p, ar, ai, indk, Na, L, K);    
#endif                  
}

static void coreRadialSphericalHarmonicsPowerDeriv(dstype *px, dstype *py, dstype *pz, dstype *ar, dstype *ai, 
        dstype *arx, dstype *aix, dstype *ary, dstype *aiy, dstype *arz, dstype *aiz, Int *indk, Int *Nnb, Int Na, Int L, Int K, Int backend)
{
    if (backend == 1)  
        cpuRadialSphericalHarmonicsPowerDeriv(px, py, pz, ar, ai, arx, aix, 
            ary, aiy, arz, aiz, indk, Nnb, Na, L, K);    
#ifdef HAVE_CUDA            
    if (backend == 2)  
        gpuRadialSphericalHarmonicsPowerDeriv(px, py, pz, ar, ai, arx, aix, 
            ary, aiy, arz, aiz, indk, Nnb, Na, L, K);
#endif                  
}

static void coreRadialSphericalHarmonicsPowerDeriv2(dstype *pd, dstype *ar, dstype *ai, dstype *arx, dstype *aix,
        dstype *ary, dstype *aiy, dstype *arz, dstype *aiz, Int *indk, Int *Nnb, Int Na, Int L, Int K, Int backend)
{
    if (backend == 1)  
        cpuRadialSphericalHarmonicsPowerDeriv2(pd, ar, ai, arx, aix, 
            ary, aiy, arz, aiz, indk, Nnb, Na, L, K);    
#ifdef HAVE_CUDA            
    if (backend == 2)  
        gpuRadialSphericalHarmonicsPowerDeriv2(pd, ar, ai, arx, aix, 
            ary, aiy, arz, aiz, indk, Nnb, Na, L, K);
#endif                  
}

static void coreRadialSphericalHarmonicsBispectrum(dstype *b, dstype *ar, dstype *ai, dstype *cg, Int *indk, 
        Int *indl, Int *indm, Int *rowm, Int Nub, Int Ncg, Int Na, Int L, Int K, Int backend)
{
    if (backend == 1)  
        cpuRadialSphericalHarmonicsBispectrum(b, ar, ai, cg, indk, indl, indm, rowm,
            Nub, Ncg, Na, L, K);    
#ifdef HAVE_CUDA            
    if (backend == 2)  
        gpuRadialSphericalHarmonicsBispectrum(b, ar, ai, cg, indk, indl, indm, rowm,
            Nub, Ncg, Na, L, K);
#endif                  
}

static void coreRadialSphericalHarmonicsBispectrumDeriv(dstype *bx, dstype *by, dstype *bz, 
        dstype *ar, dstype *ai, dstype *arx, dstype *aix, dstype *ary, dstype *aiy, dstype *arz, dstype *aiz, dstype*cg, Int *indk, Int *indl,
        Int *indm, Int *rowm, Int *Nnb, Int Na, Int Nub, Int Ncg, Int K, Int backend)
{
    if (backend == 1)  
        cpuRadialSphericalHarmonicsBispectrumDeriv(bx, by, bz, ar, ai, arx, aix, 
            ary, aiy, arz, aiz, cg, indk, indl, indm, rowm, Nnb, Na, Nub, Ncg, K);    
#ifdef HAVE_CUDA            
    if (backend == 2)  
        gpuRadialSphericalHarmonicsBispectrumDeriv(bx, by, bz, ar, ai, arx, aix, 
            ary, aiy, arz, aiz, cg, indk, indl, indm, rowm, Nnb, Na, Nub, Ncg, K);
#endif                  
}

static void coreRadialSphericalHarmonicsBispectrumDeriv2(dstype *bd, dstype *ar, dstype *ai, 
        dstype *arx, dstype *aix, dstype *ary, dstype *aiy, dstype *arz, dstype *aiz, dstype*cg, 
        Int *indk, Int *indl, Int *indm, Int *rowm, Int *Nnb, Int Na, Int Nub, Int Ncg, Int K, Int backend)
{
    if (backend == 1)  
        cpuRadialSphericalHarmonicsBispectrumDeriv2(bd, ar, ai, arx, aix, 
            ary, aiy, arz, aiz, cg, indk, indl, indm, rowm, Nnb, Na, Nub, Ncg, K);    
#ifdef HAVE_CUDA            
    if (backend == 2)  
        gpuRadialSphericalHarmonicsBispectrumDeriv2(bd, ar, ai, arx, aix, 
            ary, aiy, arz, aiz, cg, indk, indl, indm, rowm, Nnb, Na, Nub, Ncg, K);
#endif                  
}

static void coreRadialSphericalHarmonicsBasis(dstype *d, dstype *c, Int *atomtype, 
        Int Ntype, Int Na, Int Nbf, Int backend)
{
    if (backend == 1)  
        cpuRadialSphericalHarmonicsBasis(d, c, atomtype, Ntype, Na, Nbf);    
#ifdef HAVE_CUDA            
    if (backend == 2)  
        gpuRadialSphericalHarmonicsBasis(d, c, atomtype, Ntype, Na, Nbf);
#endif                  
}

static void coreRadialSphericalHarmonicsBasisDeriv(dstype *dx, dstype *dy, dstype *dz, 
        dstype *cx, dstype *cy, dstype *cz, Int *atomtype, Int *neighlist, Int *Nnb, Int Ntype, Int Na, Int Nbf, Int backend)
{
    if (backend == 1)  
        cpuRadialSphericalHarmonicsBasisDeriv(dx, dy, dz, cx, cy, cz, atomtype, 
            neighlist, Nnb, Ntype, Na, Nbf);    
#ifdef HAVE_CUDA            
    if (backend == 2)  
        gpuRadialSphericalHarmonicsBasisDeriv(dx, dy, dz, cx, cy, cz, atomtype, 
            neighlist, Nnb, Ntype, Na, Nbf);
#endif                  
}

static void coreRadialSphericalHarmonicsBasisDeriv2(dstype *dd, dstype *cd, 
        Int *atomtype, Int *neighlist, Int *Nnb, Int Ntype, Int Na, Int Nbf, Int backend)
{
    if (backend == 1)  
        cpuRadialSphericalHarmonicsBasisDeriv2(dd, cd, atomtype, neighlist, Nnb, Ntype, Na, Nbf);    
#ifdef HAVE_CUDA            
    if (backend == 2)  
        gpuRadialSphericalHarmonicsBasisDeriv2(dd, cd, atomtype, neighlist, Nnb, Ntype, Na, Nbf);
#endif                  
}

#endif  

