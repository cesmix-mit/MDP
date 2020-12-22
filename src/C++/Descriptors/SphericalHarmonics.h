#ifndef __SPHERICALHARMONICS_H__
#define __SPHERICALHARMONICS_H__

class CSphericalHarmonics {
private:
public:
    shstruct sh;        
    Int L; 
    Int K;
    Int Nub;
    Int Ncg;
    Int backend;
        
    // constructor for both CPU and GPU
    CSphericalHarmonics(Int Kin, Int Lin, Int backend); 
    
    // destructor        
    ~CSphericalHarmonics(); 
    
    // constructor for both CPU and GPU
    void Init(Int Kin, Int Lin, Int backend); 
    
    // compute
    void SphericalHarmonicsBessel(dstype *Sr, dstype *Si, dstype *x, dstype *y, dstype *z, Int N);    
    
    void SphericalHarmonicsBesselDeriv(dstype *Srx, dstype *Six, dstype *Sry, dstype *Siy, dstype *Srz, dstype *Siz, dstype *x, dstype *y, dstype *z, Int N);
        
    void RadialSphericalHarmonicsSum(dstype *ar, dstype *ai, dstype *Sr, dstype *Si, Int *numneigh, Int Na);        
        
    void RadialSphericalHarmonicsPower(dstype *p, dstype *ar, dstype *ai, Int Na);        
    
    void RadialSphericalHarmonicsPowerDeriv(dstype *px, dstype *py, dstype *pz, dstype *ar, dstype *ai, 
        dstype *arx, dstype *aix, dstype *ary, dstype *aiy, dstype *arz, dstype *aiz, Int *numneigh, Int Na);
        
    void RadialSphericalHarmonicsBispectrum(dstype *b, dstype *ar, dstype *ai, Int Na);        
    
    void RadialSphericalHarmonicsBispectrumDeriv(dstype *bx, dstype *by, dstype *bz, dstype *ar, dstype *ai, 
        dstype *arx, dstype *aix, dstype *ary, dstype *aiy, dstype *arz, dstype *aiz, Int *numneigh, Int Na);
    
    void RadialSphericalHarmonicsBasis(dstype *d, dstype *c, Int *atomtype, Int Ntype, Int Na, Int Nbf);                    
    
    void RadialSphericalHarmonicsBasisDeriv(dstype *dx, dstype *dy, dstype *dz, dstype *cx, dstype *cy, dstype *cz,
        int *atomtype, int *neighlist, int *Nnb, int Ntype, int Na, int Nbf);        
};


#endif        

