#ifndef __DESCRIPTORS
#define __DESCRIPTORS

#include "SphericalHarmonics.cpp"
#include "Descriptors.h"

// destructor 
CDescriptors::~CDescriptors()
{    
     tmp.freememory(backend);
     nb.freememory(backend);
     sys.freememory(backend);     
}

void CDescriptors::SetNeighborStruct(CConfiguration& cconf)
{
    
}

void CDescriptors::SetTempStruct(CConfiguration& cconf)
{
    
}

void CDescriptors::SetSysStruct(CConfiguration& cconf)
{
    
}

void CDescriptors::Init(CConfiguration& cconf)
{
    backend = cconf.common.backend;
    descriptor = cconf.common.descriptor;
    spectrum = cconf.common.spectrum;             
    natomtypes = cconf.config.natomtypes;   // number of atom types;    
    
    if (descriptor==0) { 
        Int K = csh.K;
        Int L = csh.L;
        Int Nub = csh.Nub;        
        if (spectrum==0) {  // power spectrum          
            nbasis = (L+1)*K*(K+1)/2;      // number of power components              
        }
        else if (spectrum==1) { // bispectrum
            nbasis = Nub*K*(K+1)/2;        // number of bispectrum components
        }
        else if (spectrum==2) { // power spectrum and bispectrum         
            Int npower = (L+1)*K*(K+1)/2;      // number of power components
            Int nbispectrum = Nub*K*(K+1)/2;   // number of bispectrum components
            nbasis = npower + nbispectrum;
        }
        else {
        }        
    }
        
    this->SetNeighborStruct(cconf);
    this->SetSysStruct(cconf);    
    this->SetTempStruct(cconf);    
} 

void CDescriptors::NeighborList(dstype* x, Int* atomtype, Int numatoms)
{
    niatoms = numatoms;
}
    
void CDescriptors::BasisFunctions(dstype *d, dstype *x, Int* atomtype, Int numatoms)
{       
    this->NeighborList(x, atomtype, numatoms);    
    
    if (descriptor==0) { 
        Int K = csh.K;
        Int L = csh.L;
        Int Nsh = K*(L+1)*(L+1);
        
        dstype *sr = &tmp.tempmem[0];
        dstype *si = &tmp.tempmem[nijatoms*Nsh];
        dstype *ar = &tmp.tempmem[(2*nijatoms)*Nsh];
        dstype *ai = &tmp.tempmem[(2*nijatoms+niatoms)*Nsh];        
        dstype *c = &tmp.tempmem[(2*nijatoms+2*niatoms)*Nsh];
        
        // xij, numneigh, atomtype, 
        csh.SphericalHarmonicsBessel(sr, si, &nb.xij[0], &nb.xij[nijatoms], &nb.xij[2*nijatoms], nijatoms);               
        csh.RadialSphericalHarmonicsSum(ar, ai, sr, si, nb.numneigh, niatoms);    
        
        if (spectrum==0) {  // power spectrum          
            csh.RadialSphericalHarmonicsPower(c, ar, ai, niatoms);   
            csh.RadialSphericalHarmonicsBasis(d, c, nb.atomtype, natomtypes, niatoms, nbasis);                    
        }
        else if (spectrum==1) { // bispectrum
            csh.RadialSphericalHarmonicsBispectrum(c, ar, ai, niatoms);   
            csh.RadialSphericalHarmonicsBasis(d, c, nb.atomtype, natomtypes, niatoms, nbasis);                    
        }
        else if (spectrum==2) { // power spectrum and bispectrum         
            Int npower = (L+1)*K*(K+1)/2;      // number of power components
            csh.RadialSphericalHarmonicsPower(c, ar, ai, niatoms);   
            csh.RadialSphericalHarmonicsBispectrum(&c[niatoms*npower], ar, ai, niatoms);   
            csh.RadialSphericalHarmonicsBasis(d, c, nb.atomtype, natomtypes, niatoms, nbasis);                    
        }
    }
}

void CDescriptors::BasisFunctionsDeriv(dstype* dd, dstype* x, Int* atomtype, Int numatoms)
{
    this->NeighborList(x, atomtype, numatoms);  
    
    if (descriptor==0) { 
        Int K = csh.K;
        Int L = csh.L;
        Int Nsh = K*(L+1)*(L+1);
        dstype *sr = &tmp.tempmem[0];
        dstype *si = &tmp.tempmem[nijatoms*Nsh];
        dstype *srx = &tmp.tempmem[2*nijatoms*Nsh];
        dstype *six = &tmp.tempmem[3*nijatoms*Nsh];
        dstype *sry = &tmp.tempmem[4*nijatoms*Nsh];
        dstype *siy = &tmp.tempmem[5*nijatoms*Nsh];
        dstype *srz = &tmp.tempmem[6*nijatoms*Nsh];
        dstype *siz = &tmp.tempmem[7*nijatoms*Nsh];        
        dstype *ar = &tmp.tempmem[(8*nijatoms)*Nsh];
        dstype *ai = &tmp.tempmem[(8*nijatoms + niatoms)*Nsh];         
        dstype *cd = &tmp.tempmem[(8*nijatoms + 2*niatoms)*Nsh];
        
        csh.SphericalHarmonicsBesselWithDeriv(sr, si, srx, six, sry, siy, srz, siz,
                &nb.xij[0], &nb.xij[nijatoms], &nb.xij[2*nijatoms], nijatoms);               
        csh.RadialSphericalHarmonicsSum(ar, ai, sr, si, nb.numneigh, niatoms);    
        
        if (spectrum==0) {  // power spectrum    
            csh.RadialSphericalHarmonicsPowerDeriv2(cd, ar, ai, srx, six, sry, siy, srz, siz, nb.numneigh, niatoms);
            csh.RadialSphericalHarmonicsBasisDeriv2(dd, cd, nb.atomtype, nb.neighlist, nb.numneigh, 
                    natomtypes, niatoms, nbasis);            
        }
        else if (spectrum==1) { // bispectrum
            csh.RadialSphericalHarmonicsBispectrumDeriv2(cd, ar, ai, srx, six, sry, siy, srz, siz, nb.numneigh, niatoms);
            csh.RadialSphericalHarmonicsBasisDeriv2(dd, cd, nb.atomtype, nb.neighlist, nb.numneigh, 
                    natomtypes, niatoms, nbasis);                        
        }
        else if (spectrum==2) { // power spectrum and bispectrum         
            Int npower = (L+1)*K*(K+1)/2;      // number of power components
            csh.RadialSphericalHarmonicsPowerDeriv2(cd, ar, ai, srx, six, sry, siy, srz, siz, nb.numneigh, niatoms);
            csh.RadialSphericalHarmonicsBispectrumDeriv2(&cd[3*nijatoms*npower], ar, ai, srx, six, sry, siy, srz, siz, nb.numneigh, niatoms);
            csh.RadialSphericalHarmonicsBasisDeriv2(dd, cd, nb.atomtype, nb.neighlist, nb.numneigh, 
                    natomtypes, niatoms, nbasis);                        
        }
    }    
}

void CDescriptors::BasisFunctionsWithDeriv(dstype* d, dstype* dd, dstype* x, Int* atomtype, Int numatoms)
{
    this->NeighborList(x, atomtype, numatoms);  
    
    if (descriptor==0) { 
        Int K = csh.K;
        Int L = csh.L;
        Int Nsh = K*(L+1)*(L+1);
        dstype *sr = &tmp.tempmem[0];
        dstype *si = &tmp.tempmem[nijatoms*Nsh];
        dstype *srx = &tmp.tempmem[2*nijatoms*Nsh];
        dstype *six = &tmp.tempmem[3*nijatoms*Nsh];
        dstype *sry = &tmp.tempmem[4*nijatoms*Nsh];
        dstype *siy = &tmp.tempmem[5*nijatoms*Nsh];
        dstype *srz = &tmp.tempmem[6*nijatoms*Nsh];
        dstype *siz = &tmp.tempmem[7*nijatoms*Nsh];        
        dstype *ar = &tmp.tempmem[(8*nijatoms)*Nsh];
        dstype *ai = &tmp.tempmem[(8*nijatoms + niatoms)*Nsh];         
        dstype *c = &tmp.tempmem[(8*nijatoms + 2*niatoms)*Nsh];
        dstype *cd = &tmp.tempmem[(8*nijatoms + 2*niatoms)*Nsh + niatoms*nbasis];
        
        csh.SphericalHarmonicsBesselWithDeriv(sr, si, srx, six, sry, siy, srz, siz,
                &nb.xij[0], &nb.xij[nijatoms], &nb.xij[2*nijatoms], nijatoms);               
        csh.RadialSphericalHarmonicsSum(ar, ai, sr, si, nb.numneigh, niatoms);    
        
        if (spectrum==0) {  // power spectrum    
            csh.RadialSphericalHarmonicsPower(c, ar, ai, niatoms);               
            csh.RadialSphericalHarmonicsPowerDeriv2(cd, ar, ai, srx, six, sry, siy, srz, siz, nb.numneigh, niatoms);
            csh.RadialSphericalHarmonicsBasis(d, c, nb.atomtype, natomtypes, niatoms, nbasis);                                
            csh.RadialSphericalHarmonicsBasisDeriv2(dd, cd, nb.atomtype, nb.neighlist, nb.numneigh, 
                    natomtypes, niatoms, nbasis);            
        }
        else if (spectrum==1) { // bispectrum
            csh.RadialSphericalHarmonicsBispectrum(c, ar, ai, niatoms);               
            csh.RadialSphericalHarmonicsBispectrumDeriv2(cd, ar, ai, srx, six, sry, siy, srz, siz, nb.numneigh, niatoms);
            csh.RadialSphericalHarmonicsBasis(d, c, nb.atomtype, natomtypes, niatoms, nbasis);                   
            csh.RadialSphericalHarmonicsBasisDeriv2(dd, cd, nb.atomtype, nb.neighlist, nb.numneigh, 
                    natomtypes, niatoms, nbasis);                        
        }
        else if (spectrum==2) { // power spectrum and bispectrum         
            Int npower = (L+1)*K*(K+1)/2;      // number of power components
            csh.RadialSphericalHarmonicsPower(c, ar, ai, niatoms);   
            csh.RadialSphericalHarmonicsBispectrum(&c[niatoms*npower], ar, ai, niatoms);               
            csh.RadialSphericalHarmonicsPowerDeriv2(cd, ar, ai, srx, six, sry, siy, srz, siz, nb.numneigh, niatoms);
            csh.RadialSphericalHarmonicsBispectrumDeriv2(&cd[3*nijatoms*npower], ar, ai, srx, six, sry, siy, srz, siz, nb.numneigh, niatoms);
            csh.RadialSphericalHarmonicsBasis(d, c, nb.atomtype, natomtypes, niatoms, nbasis);    
            csh.RadialSphericalHarmonicsBasisDeriv2(dd, cd, nb.atomtype, nb.neighlist, nb.numneigh, 
                    natomtypes, niatoms, nbasis);                        
        }
    }    
}

void CDescriptors::Energy(dstype e, dstype* x, dstype* coeff, Int* atomtype, Int numatoms, Int ncoeff)
{    
    dstype *d = &tmp.tempmem[0];
    this->BasisFunctions(d, x, atomtype, numatoms);
    DOT(cublasHandle, ncoeff, d, inc1, coeff, inc1, &e, backend);
}

void CDescriptors::Forces(dstype* f, dstype* x, dstype* coeff, Int* atomtype, Int numatoms, Int ncoeff)
{
    dstype *dd = &tmp.tempmem[0];
    this->BasisFunctionsDeriv(dd, x, atomtype, numatoms);     
    PGEMNV(cublasHandle, 3*numatoms, ncoeff, &one, dd, 3*numatoms, coeff, inc1, &zero, f, inc1, backend);        
}

void CDescriptors::EnergyForces(dstype e, dstype* f, dstype* x, dstype* coeff, Int* atomtype, Int numatoms, Int ncoeff)
{
    dstype *d = &tmp.tempmem[0];
    dstype *dd = &tmp.tempmem[ncoeff];
    this->BasisFunctionsWithDeriv(d, dd, x, atomtype, numatoms);     
    DOT(cublasHandle, ncoeff, d, inc1, coeff, inc1, &e, backend);
    PGEMNV(cublasHandle, 3*numatoms, ncoeff, &one, dd, 3*numatoms, coeff, inc1, &zero, f, inc1, backend);            
}

void CDescriptors::Stresses(dstype* s, dstype* x, dstype* coeff, Int* atomtype, Int numatoms, Int ncoeff)
{
}
        
void CDescriptors::EnergyForcesStresses(dstype e, dstype* f, dstype* s, dstype* x, dstype* coeff, Int* atomtype, Int numatoms, Int ncoeff)
{
}

#endif        

