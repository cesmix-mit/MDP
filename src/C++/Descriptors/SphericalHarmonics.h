#ifndef __SPHERICALHARMONICS_H__
#define __SPHERICALHARMONICS_H__

class CSphericalHarmonics {
private:
public:
    shstruct sh;        
    Int L;  // the maximum degree of spherical harmonics 
    Int K;  // the number of zeros of spherical Bessel functions
    Int Nub;// number of non-zero unqiue bispectrum compoments of spherical harmonics 
    Int Ncg;// the total number of non-zero Clebsch-Gordan coefficients 
    Int backend; // computing platform
        
    // constructor for both CPU and GPU
    CSphericalHarmonics(Int Kin, Int Lin, Int backendin); 
    
    // destructor        
    ~CSphericalHarmonics(); 
    
    // constructor for both CPU and GPU
    void Init(Int Kin, Int Lin, Int backendin); 
    
    void SphericalHarmonicsBessel(dstype *Sr, dstype *Si, dstype *xij, Int N);    
    
    // compute derivatives of spherical harmonics Bessel functions
    void SphericalHarmonicsBesselWithDeriv(dstype *Sr, dstype *Si, dstype *Srx, dstype *Six, dstype *Sry, dstype *Siy, dstype *Srz, dstype *Siz, dstype *xij, Int N);
    
    // compute spherical harmonics Bessel functions for N atoms
    void SphericalHarmonicsBessel(dstype *Sr, dstype *Si, dstype *x, dstype *y, dstype *z, Int N);    
    
    // compute derivatives of spherical harmonics Bessel functions
    void SphericalHarmonicsBesselDeriv(dstype *Srx, dstype *Six, dstype *Sry, dstype *Siy, dstype *Srz, dstype *Siz, dstype *x, dstype *y, dstype *z, Int N);
    
    // compute derivatives of spherical harmonics Bessel functions
    void SphericalHarmonicsBesselWithDeriv(dstype *Sr, dstype *Si, dstype *Srx, dstype *Six, dstype *Sry, dstype *Siy, dstype *Srz, dstype *Siz, dstype *x, dstype *y, dstype *z, Int N);
    
    // compute radial spherical harmonics functions and sum them over neighbors
    void RadialSphericalHarmonicsSum(dstype *ar, dstype *ai, dstype *Sr, dstype *Si, Int *numneigh, Int Na);        
        
    // compute power spectrum components of radial spherical harmonics functions for Na atoms
    void RadialSphericalHarmonicsPower(dstype *p, dstype *ar, dstype *ai, Int Na);        
    
    // compute derivatives of power spectrum components of radial spherical harmonics functions for Na atoms
    void RadialSphericalHarmonicsPowerDeriv(dstype *px, dstype *py, dstype *pz, dstype *ar, dstype *ai, 
        dstype *arx, dstype *aix, dstype *ary, dstype *aiy, dstype *arz, dstype *aiz, Int *numneigh, Int Na);

    // compute derivatives of power spectrum components of radial spherical harmonics functions for Na atoms
    void RadialSphericalHarmonicsPowerDeriv2(dstype *pd, dstype *ar, dstype *ai, 
        dstype *arx, dstype *aix, dstype *ary, dstype *aiy, dstype *arz, dstype *aiz, Int *numneigh, Int Na);
    
    // compute bispectrum components of radial spherical harmonics functions for Na atoms
    void RadialSphericalHarmonicsBispectrum(dstype *b, dstype *ar, dstype *ai, Int Na);        
    
    // compute derivatives of bispectrum components of radial spherical harmonics functions for Na atoms
    void RadialSphericalHarmonicsBispectrumDeriv(dstype *bx, dstype *by, dstype *bz, dstype *ar, dstype *ai, 
        dstype *arx, dstype *aix, dstype *ary, dstype *aiy, dstype *arz, dstype *aiz, Int *numneigh, Int Na);

    // compute derivatives of bispectrum components of radial spherical harmonics functions for Na atoms
    void RadialSphericalHarmonicsBispectrumDeriv2(dstype *bd, dstype *ar, dstype *ai, 
        dstype *arx, dstype *aix, dstype *ary, dstype *aiy, dstype *arz, dstype *aiz, Int *numneigh, Int Na);
    
    // compute invariant basis functions for Na atoms
    void RadialSphericalHarmonicsBasis(dstype *d, dstype *c, Int *atomtype, Int Ntype, Int Na, Int Nbf);                    
    
    // compute derivative of invariant basis functions for Na atoms
    void RadialSphericalHarmonicsBasisDeriv(dstype *dx, dstype *dy, dstype *dz, dstype *cx, dstype *cy, dstype *cz,
        int *atomtype, int *neighlist, int *numneigh, int Ntype, int Na, int Nbf);     
    
    // compute derivative of invariant basis functions for Na atoms
    void RadialSphericalHarmonicsBasisDeriv2(dstype *dd, dstype *cd,
        int *atomtype, int *neighlist, int *numneigh, int Ntype, int Na, int Nbf);            
};


#endif        

