#ifndef __SPHERICALHARMONICS
#define __SPHERICALHARMONICS

#include "SphericalHarmonics.h"
#include "SphericalHarmonics_impl.cpp"

void CSphericalHarmonics::Init(Int Kin, Int Lin, Int backendin) 
{
    L = Lin;
    K = Kin;
    backend = backendin;
    InitSphericalHarmonic(sh, Kin, Lin, backendin); 
    Nub = sh.Nub;
    Ncg = sh.Ncg;    
}

// CPU constructor
CSphericalHarmonics::CSphericalHarmonics(Int Kin, Int Lin, Int backend) 
{
    this->Init(Kin, Lin, backend);
}

// destructor 
CSphericalHarmonics::~CSphericalHarmonics()
{    
     sh.freememory(backend);
}

void CSphericalHarmonics::SphericalHarmonicsBessel(dstype *Sr, dstype *Si, dstype *xij, Int N)
{
    implSphericalHarmonicsBessel(Sr, Si, xij, N, backend, sh);
}

void CSphericalHarmonics::SphericalHarmonicsBesselWithDeriv(dstype *Sr, dstype *Si, dstype *Srx, dstype *Six, dstype *Sry, dstype *Siy, dstype *Srz, dstype *Siz, dstype *xij, Int N)
{
    implSphericalHarmonicsBesselWithDeriv(Sr, Si, Srx, Six, Sry, Siy, Srz, Siz, xij, sh.x0, sh.P, N, backend, sh);    
}

void CSphericalHarmonics::SphericalHarmonicsBessel(dstype *Sr, dstype *Si, dstype *x, dstype *y, dstype *z, Int N)
{
    implSphericalHarmonicsBessel(Sr, Si, x, y, z, N, backend, sh);
}

void CSphericalHarmonics::SphericalHarmonicsBesselDeriv(dstype *Srx, dstype *Six, dstype *Sry, dstype *Siy, dstype *Srz, dstype *Siz, dstype *x, dstype *y, dstype *z, Int N)
{
    implSphericalHarmonicsBesselDeriv(Srx, Six, Sry, Siy, Srz, Siz, x, y, z, N, backend, sh);    
}

void CSphericalHarmonics::SphericalHarmonicsBesselWithDeriv(dstype *Sr, dstype *Si, dstype *Srx, dstype *Six, dstype *Sry, dstype *Siy, dstype *Srz, dstype *Siz, dstype *x, dstype *y, dstype *z, Int N)
{
    implSphericalHarmonicsBesselWithDeriv(Sr, Si, Srx, Six, Sry, Siy, Srz, Siz, x, y, z, N, backend, sh);    
}

void CSphericalHarmonics::RadialSphericalHarmonicsSum(dstype *ar, dstype *ai, dstype *Sr, dstype *Si, Int *numneigh, Int Na)
{
    implRadialSphericalHarmonicsSum(ar, ai, Sr, Si, numneigh, Na, backend, sh);
}

void CSphericalHarmonics::RadialSphericalHarmonicsPower(dstype *p, dstype *ar, dstype *ai, Int Na)
{
    implRadialSphericalHarmonicsPower(p, ar, ai, Na, backend, sh);
}

void CSphericalHarmonics::RadialSphericalHarmonicsPowerDeriv(dstype *px, dstype *py, dstype *pz, dstype *ar, dstype *ai, 
        dstype *arx, dstype *aix, dstype *ary, dstype *aiy, dstype *arz, dstype *aiz, Int *numneigh, Int Na)
{
    implRadialSphericalHarmonicsPowerDeriv(px, py, pz, ar, ai, arx, aix, ary, aiy, arz, aiz, numneigh, Na, backend, sh);
}

void CSphericalHarmonics::RadialSphericalHarmonicsPowerDeriv2(dstype *pd, dstype *ar, dstype *ai, 
        dstype *arx, dstype *aix, dstype *ary, dstype *aiy, dstype *arz, dstype *aiz, Int *numneigh, Int Na)
{
    implRadialSphericalHarmonicsPowerDeriv2(pd, ar, ai, arx, aix, ary, aiy, arz, aiz, numneigh, Na, backend, sh);
}

void CSphericalHarmonics::RadialSphericalHarmonicsBispectrum(dstype *b, dstype *ar, dstype *ai, Int Na)
{    
    implRadialSphericalHarmonicsBispectrum(b, ar, ai, Na, backend, sh);
}

void CSphericalHarmonics::RadialSphericalHarmonicsBispectrumDeriv(dstype *bx, dstype *by, dstype *bz, dstype *ar, dstype *ai, 
        dstype *arx, dstype *aix, dstype *ary, dstype *aiy, dstype *arz, dstype *aiz, Int *numneigh, Int Na)
{
    implRadialSphericalHarmonicsBispectrumDeriv(bx, by, bz, ar, ai, arx, aix, ary, aiy, arz, aiz, 
             numneigh, Na, backend, sh);
}

void CSphericalHarmonics::RadialSphericalHarmonicsBispectrumDeriv2(dstype *bd, dstype *ar, dstype *ai, 
        dstype *arx, dstype *aix, dstype *ary, dstype *aiy, dstype *arz, dstype *aiz, Int *numneigh, Int Na)
{
    implRadialSphericalHarmonicsBispectrumDeriv2(bd, ar, ai, arx, aix, ary, aiy, arz, aiz, 
             numneigh, Na, backend, sh);
}

void CSphericalHarmonics::RadialSphericalHarmonicsBasis(dstype *d, dstype *c, Int *atomtype, Int Ntype, Int Na, Int Nbf)
{
    implRadialSphericalHarmonicsBasis(d, c, atomtype, Ntype, Na, Nbf, backend);
}

void CSphericalHarmonics::RadialSphericalHarmonicsBasisDeriv(dstype *dx, dstype *dy, dstype *dz, dstype *cx, dstype *cy, dstype *cz,
        int *atomtype, int *neighlist, int *numneigh, int Ntype, int Na, int Nbf)
{
    implRadialSphericalHarmonicsBasisDeriv(dx, dy, dz, cx, cy, cx, atomtype, neighlist, numneigh, Ntype, Na, Nbf, backend);
}

void CSphericalHarmonics::RadialSphericalHarmonicsBasisDeriv2(dstype *dd, dstype *cd,
        int *atomtype, int *neighlist, int *numneigh, int Ntype, int Na, int Nbf)
{
    implRadialSphericalHarmonicsBasisDeriv2(dd, cd, atomtype, neighlist, numneigh, Ntype, Na, Nbf, backend);
}

#endif        
