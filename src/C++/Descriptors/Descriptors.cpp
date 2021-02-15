#ifndef __DESCRIPTORS
#define __DESCRIPTORS

#include "SphericalHarmonics.cpp"
#include "Descriptors.h"
#include "Descriptors_impl.cpp"

// // destructor 
// CDescriptors::~CDescriptors()
// {    
//      tmp.freememory(backend);
//      nb.freememory(backend);
//      sys.freememory(backend);     
// }

// void CDescriptors::SetNeighborStruct(CConfiguration& cconf)
// {
// //     nimax = cconf.common.nimax;
// //     njmax = cconf.common.njmax;
// //     
// //     nb.freememory(cconf.common.backend);
// //     
// //     TemplateMalloc(&nb.bcs, 6, cconf.common.backend); 
// //     TemplateMalloc(&nb.simbox, 9, cconf.common.backend); 
// //     TemplateMalloc(&nb.neighnum, nimax, cconf.common.backend); 
// //     //TemplateMalloc(&nb.atomtype, nimax, cconf.common.backend); 
// //     TemplateMalloc(&nb.neighlist, nimax*njmax, cconf.common.backend); 
// //     TemplateMalloc(&nb.xij, nimax*njmax*3, cconf.common.backend); 
// //     
// //     if (backend==2) { // GPU
// // #ifdef HAVE_CUDA        
// //         CHECK( cudaMemcpy(nb.bcs, cconf.config.bcs, 6*sizeof(Int), cudaMemcpyHostToDevice ) );              
// //         CHECK( cudaMemcpy(nb.simbox, cconf.config.simbox, 9*sizeof(dstype), cudaMemcpyHostToDevice ) );              
// // #endif        
// //     }    
// //     else { // CPU
// //         for (int i=0; i<6; i++)
// //             nb.bcs[i] = cconf.config.bcs[i];                
// //         for (int i=0; i<9; i++)
// //             nb.simbox[i] = cconf.config.simbox[i];                
// //     }        
// }
// 
// void CDescriptors::SetTempStruct(CConfiguration& cconf)
// {
//     
// }
// 
// void CDescriptors::SetSysStruct(CConfiguration& cconf)
// {
//     
// }

// void CDescriptors::Init(CConfiguration& cconf)
// {
//     backend = cconf.common.backend;
//     descriptor = cconf.common.descriptor;
//     spectrum = cconf.common.spectrum;             
//     natomtypes = cconf.common.natomtypes;   // number of atom types;    
//     
//     if (descriptor==0) { 
//         Int K = csh.K;
//         Int L = csh.L;
//         Int Nub = csh.Nub;        
//         if (spectrum==0) {  // power spectrum          
//             nbasis = (L+1)*K*(K+1)/2;      // number of power components              
//         }
//         else if (spectrum==1) { // bispectrum
//             nbasis = Nub*K*(K+1)/2;        // number of bispectrum components
//         }
//         else if (spectrum==2) { // power spectrum and bispectrum         
//             Int npower = (L+1)*K*(K+1)/2;      // number of power components
//             Int nbispectrum = Nub*K*(K+1)/2;   // number of bispectrum components
//             nbasis = npower + nbispectrum;
//         }
//         else {
//         }        
//     }
//         
//     this->SetNeighborStruct(cconf);
//     this->SetSysStruct(cconf);    
//     this->SetTempStruct(cconf);    
// } 

// void CDescriptors::BuildVerletList(dstype* x, Int* atomtype, Int numatoms)
// { 
//    implBuildVerletList(x, atomtype, numatoms, backend, nb, tmp);
    
//     nb.inum = numatoms; // number of atoms in the simulation box
//        
//     GhostAtoms(nb.glistnum, x, nb.pimages, nb.rbvertices, nb.s2rmap, 
//             nb.inum, nb.pnum, nb.dim, backend);   
//     
//     CreateAtomList(nb.ilist, nb.glistnumsum, nb.glistnum, atomtype, x, nb.pimages, nb.rbvertices, nb.s2rmap, 
//             nb.inum, nb.pnum, nb.dim, backend);   
//     
//     // number of ghost atoms
//     nb.gnum = IntArrayGetValueAtIndex(nb.glistnumsum, nb.inum, backend);
//     
//     CellList(nb.clist, nb.c2anum, x, nb.eta1, nb.eta2, nb.eta3, nb.rbvertices, nb.cellnum, 
//             nb.inum, nb.gnum, nb.dim, backend);   
//     
//     Cell2AtomList(nb.c2alist, nb.c2anumsum, nb.c2anum, nb.clist, nb.cellnum, 
//             nb.inum, nb.gnum, nb.dim, backend);   
//     
//     if (nb.cutofftype==1) { // pairwise cut-off geometries between pairs of atom types
//         VerletAtoms(nb.verletnum, x, nb.ellipsoid, atomtype, nb.ilist, nb.clist, nb.c2alist, 
//                 nb.c2anum, nb.c2anumsum, nb.cellnum, nb.ntype, nb.inum, nb.dim, backend);   
// 
//         CreateVerletList(nb.verletlist, x, nb.ellipsoid, nb.verletnum, nb.verletnumsum, atomtype, 
//                 nb.ilist, nb.clist, nb.c2alist, nb.c2anum, nb.c2anumsum, nb.cellnum, nb.ntype, 
//                 nb.inum, nb.dim, backend);               
//     }
//     else if (nb.cutofftype==0) { // one cut-off geometry for all atoms
//         VerletAtoms(nb.verletnum, x, nb.ellipsoid, nb.ilist, nb.clist, nb.c2alist, nb.c2anum, nb.c2anumsum, 
//                     nb.cellnum, nb.inum, nb.dim, backend);   
// 
//         CreateVerletList(nb.verletlist, x, nb.ellipsoid, nb.verletnum, nb.verletnumsum, nb.ilist, nb.clist,
//                 nb.c2alist, nb.c2anum, nb.c2anumsum, nb.cellnum, nb.inum, nb.dim, backend);                       
//     }        
// }

// void CDescriptors::NeighborList(dstype* x, Int* atomtype, Int numatoms)
// {    
//     if (nb.cutofftype==1) { // pairwise cut-off geometries between pairs of atom types
//         if (nb.neightype==0) {
//             FullNeighNum(nb.neighnum, x, nb.ellipsoid, atomtype, nb.ilist, nb.verletlist, nb.verletnum, 
//                     nb.verletnumsum, nb.ntype, nb.inum, nb.dim, backend);   
//             FullNeighList(nb.neighlist, x, nb.ellipsoid, nb.neighnum, nb.neighnumsum, atomtype, nb.ilist, 
//                     nb.verletlist, nb.verletnum, nb.verletnumsum, nb.ntype, nb.inum, nb.dim, backend);  
//         }
//         else {
//             HalfNeighNum(nb.neighnum, x, nb.ellipsoid, atomtype, nb.ilist, nb.verletlist, nb.verletnum, 
//                     nb.verletnumsum, nb.ntype, nb.inum, nb.dim, backend);   
//             HalfNeighList(nb.neighlist, x, nb.ellipsoid, nb.neighnum, nb.neighnumsum, atomtype, nb.ilist, 
//                     nb.verletlist, nb.verletnum, nb.verletnumsum, nb.ntype, nb.inum, nb.dim, backend);              
//         }                
//     }
//     else if (nb.cutofftype==0) { // one cut-off geometry for all atoms
//         if (nb.neightype==0) {
//             FullNeighNum(nb.neighnum, x, nb.ellipsoid, nb.ilist, nb.verletlist, nb.verletnum, 
//                     nb.verletnumsum, nb.inum, nb.dim, backend);          
//             FullNeighList(nb.neighlist, x, nb.ellipsoid, nb.neighnum, nb.neighnumsum, nb.ilist, 
//                     nb.verletlist, nb.verletnum, nb.verletnumsum, nb.inum, nb.dim, backend); 
//         }
//         else {
//             HalfNeighNum(nb.neighnum, x, nb.ellipsoid, nb.ilist, nb.verletlist, nb.verletnum, 
//                     nb.verletnumsum, nb.inum, nb.dim, backend);          
//             HalfNeighList(nb.neighlist, x, nb.ellipsoid, nb.neighnum, nb.neighnumsum, nb.ilist, 
//                     nb.verletlist, nb.verletnum, nb.verletnumsum, nb.inum, nb.dim, backend);             
//         }                
//     }        
// }
// 
// int CDescriptors::NeighborPairs(dstype* x, Int* atomtype, Int *ilist, Int ilistsize, Int typei, Int typej)
// {        
//     int tnum;
//     if (nb.pairtype==0) { //
//         GetNeighPairs(nb.xij, x, nb.anum, nb.anumsum, nb.ai, nb.aj, nb.ti, nb.tj, ilist, nb.neighlist, 
//                 nb.neighnum, nb.neighnumsum, atomtype, ilistsize, nb.dim, backend);            
//         tnum = ilistsize;
//     }
//     else if (nb.pairtype==1) { // group self atoms according to their atom types
//         tnum = GetNeighPairs(nb.xij, x, nb.anum, nb.anumsum, nb.ai, nb.aj, nb.ti, nb.tj, ilist, nb.tlist,
//                 nb.neighlist, nb.neighnum, nb.neighnumsum, atomtype, typei, ilistsize, nb.dim, backend);                           
//     }                       
//     else if (nb.pairtype==2) { // group both self atoms and neighbors according to their atom types
//         tnum = GetNeighPairs(nb.xij, x, nb.anum, nb.anumsum, nb.ai, nb.aj, nb.ti, nb.tj, ilist, nb.tlist,
//                 nb.neighlist, nb.neighnum, nb.neighnumsum, atomtype, typei, typej, ilistsize, nb.dim, backend);                           
//     }                       
//     return tnum;
// }

// void CDescriptors::BasisFunctions(dstype *d, dstype *x, Int *atomtype, Int *ilist, Int listsize)
// {           
//         
//     if (descriptor==0) { 
//         this->NeighborPairs(x, atomtype, ilist, listsize, typei, typej);    
//         
//         Int K = csh.K;
//         Int L = csh.L;
//         Int Nsh = K*(L+1)*(L+1);
//         
//         dstype *sr = &tmp.tmpmem[0];
//         dstype *si = &tmp.tmpmem[nijatoms*Nsh];
//         dstype *ar = &tmp.tmpmem[(2*nijatoms)*Nsh];
//         dstype *ai = &tmp.tmpmem[(2*nijatoms+niatoms)*Nsh];        
//         dstype *c = &tmp.tmpmem[(2*nijatoms+2*niatoms)*Nsh];
//         
//         // xij, neighnum, atomtype, 
//         csh.SphericalHarmonicsBessel(sr, si, nb.xij, nijatoms);               
//         csh.RadialSphericalHarmonicsSum(ar, ai, sr, si, nb.neighnum, niatoms);    
//         
//         if (spectrum==0) {  // power spectrum          
//             csh.RadialSphericalHarmonicsPower(c, ar, ai, niatoms);   
//             csh.RadialSphericalHarmonicsBasis(d, c, atomtype, natomtypes, niatoms, nbasis);                    
//         }
//         else if (spectrum==1) { // bispectrum
//             csh.RadialSphericalHarmonicsBispectrum(c, ar, ai, niatoms);   
//             csh.RadialSphericalHarmonicsBasis(d, c, atomtype, natomtypes, niatoms, nbasis);                    
//         }
//         else if (spectrum==2) { // power spectrum and bispectrum         
//             Int npower = (L+1)*K*(K+1)/2;      // number of power components
//             csh.RadialSphericalHarmonicsPower(c, ar, ai, niatoms);   
//             csh.RadialSphericalHarmonicsBispectrum(&c[niatoms*npower], ar, ai, niatoms);   
//             csh.RadialSphericalHarmonicsBasis(d, c, atomtype, natomtypes, niatoms, nbasis);                    
//         }
//     }
// }
// 
// void CDescriptors::BasisFunctions(dstype *d, dstype *x, Int* atomtype, Int numatoms)
// {       
//     //this->NeighborPairs(x, atomtype, ilist, listsize, typei, typej);    
//         
//     if (descriptor==0) { 
//         Int K = csh.K;
//         Int L = csh.L;
//         Int Nsh = K*(L+1)*(L+1);
//         
//         dstype *sr = &tmp.tmpmem[0];
//         dstype *si = &tmp.tmpmem[nijatoms*Nsh];
//         dstype *ar = &tmp.tmpmem[(2*nijatoms)*Nsh];
//         dstype *ai = &tmp.tmpmem[(2*nijatoms+niatoms)*Nsh];        
//         dstype *c = &tmp.tmpmem[(2*nijatoms+2*niatoms)*Nsh];
//         
//         // xij, neighnum, atomtype, 
//         csh.SphericalHarmonicsBessel(sr, si, nb.xij, nijatoms);               
//         csh.RadialSphericalHarmonicsSum(ar, ai, sr, si, nb.neighnum, niatoms);    
//         
//         if (spectrum==0) {  // power spectrum          
//             csh.RadialSphericalHarmonicsPower(c, ar, ai, niatoms);   
//             csh.RadialSphericalHarmonicsBasis(d, c, atomtype, natomtypes, niatoms, nbasis);                    
//         }
//         else if (spectrum==1) { // bispectrum
//             csh.RadialSphericalHarmonicsBispectrum(c, ar, ai, niatoms);   
//             csh.RadialSphericalHarmonicsBasis(d, c, atomtype, natomtypes, niatoms, nbasis);                    
//         }
//         else if (spectrum==2) { // power spectrum and bispectrum         
//             Int npower = (L+1)*K*(K+1)/2;      // number of power components
//             csh.RadialSphericalHarmonicsPower(c, ar, ai, niatoms);   
//             csh.RadialSphericalHarmonicsBispectrum(&c[niatoms*npower], ar, ai, niatoms);   
//             csh.RadialSphericalHarmonicsBasis(d, c, atomtype, natomtypes, niatoms, nbasis);                    
//         }
//     }
// }
// 
// void CDescriptors::BasisFunctionsDeriv(dstype* dd, dstype* x, Int* atomtype, Int numatoms)
// {
//     //this->NeighborList(x, atomtype, numatoms);  
//     
//     if (descriptor==0) { 
//         Int K = csh.K;
//         Int L = csh.L;
//         Int Nsh = K*(L+1)*(L+1);
//         dstype *sr = &tmp.tmpmem[0];
//         dstype *si = &tmp.tmpmem[nijatoms*Nsh];
//         dstype *srx = &tmp.tmpmem[2*nijatoms*Nsh];
//         dstype *six = &tmp.tmpmem[3*nijatoms*Nsh];
//         dstype *sry = &tmp.tmpmem[4*nijatoms*Nsh];
//         dstype *siy = &tmp.tmpmem[5*nijatoms*Nsh];
//         dstype *srz = &tmp.tmpmem[6*nijatoms*Nsh];
//         dstype *siz = &tmp.tmpmem[7*nijatoms*Nsh];        
//         dstype *ar = &tmp.tmpmem[(8*nijatoms)*Nsh];
//         dstype *ai = &tmp.tmpmem[(8*nijatoms + niatoms)*Nsh];         
//         dstype *cd = &tmp.tmpmem[(8*nijatoms + 2*niatoms)*Nsh];
//         
//         csh.SphericalHarmonicsBesselWithDeriv(sr, si, srx, six, sry, siy, srz, siz, nb.xij, nijatoms);               
//         csh.RadialSphericalHarmonicsSum(ar, ai, sr, si, nb.neighnum, niatoms);    
//         
//         if (spectrum==0) {  // power spectrum    
//             csh.RadialSphericalHarmonicsPowerDeriv2(cd, ar, ai, srx, six, sry, siy, srz, siz, nb.neighnum, niatoms);
//             csh.RadialSphericalHarmonicsBasisDeriv2(dd, cd, atomtype, nb.neighlist, nb.neighnum, 
//                     natomtypes, niatoms, nbasis);            
//         }
//         else if (spectrum==1) { // bispectrum
//             csh.RadialSphericalHarmonicsBispectrumDeriv2(cd, ar, ai, srx, six, sry, siy, srz, siz, nb.neighnum, niatoms);
//             csh.RadialSphericalHarmonicsBasisDeriv2(dd, cd, atomtype, nb.neighlist, nb.neighnum, 
//                     natomtypes, niatoms, nbasis);                        
//         }
//         else if (spectrum==2) { // power spectrum and bispectrum         
//             Int npower = (L+1)*K*(K+1)/2;      // number of power components
//             csh.RadialSphericalHarmonicsPowerDeriv2(cd, ar, ai, srx, six, sry, siy, srz, siz, nb.neighnum, niatoms);
//             csh.RadialSphericalHarmonicsBispectrumDeriv2(&cd[3*nijatoms*npower], ar, ai, srx, six, sry, siy, srz, siz, nb.neighnum, niatoms);
//             csh.RadialSphericalHarmonicsBasisDeriv2(dd, cd, atomtype, nb.neighlist, nb.neighnum, 
//                     natomtypes, niatoms, nbasis);                        
//         }
//     }    
// }
// 
// void CDescriptors::BasisFunctionsWithDeriv(dstype* d, dstype* dd, dstype* x, Int* atomtype, Int numatoms)
// {
//     //this->NeighborList(x, atomtype, numatoms);  
//     
//     if (descriptor==0) { 
//         Int K = csh.K;
//         Int L = csh.L;
//         Int Nsh = K*(L+1)*(L+1);
//         dstype *sr = &tmp.tmpmem[0];
//         dstype *si = &tmp.tmpmem[nijatoms*Nsh];
//         dstype *srx = &tmp.tmpmem[2*nijatoms*Nsh];
//         dstype *six = &tmp.tmpmem[3*nijatoms*Nsh];
//         dstype *sry = &tmp.tmpmem[4*nijatoms*Nsh];
//         dstype *siy = &tmp.tmpmem[5*nijatoms*Nsh];
//         dstype *srz = &tmp.tmpmem[6*nijatoms*Nsh];
//         dstype *siz = &tmp.tmpmem[7*nijatoms*Nsh];        
//         dstype *ar = &tmp.tmpmem[(8*nijatoms)*Nsh];
//         dstype *ai = &tmp.tmpmem[(8*nijatoms + niatoms)*Nsh];         
//         dstype *c = &tmp.tmpmem[(8*nijatoms + 2*niatoms)*Nsh];
//         dstype *cd = &tmp.tmpmem[(8*nijatoms + 2*niatoms)*Nsh + niatoms*nbasis];
//         
//         csh.SphericalHarmonicsBesselWithDeriv(sr, si, srx, six, sry, siy, srz, siz,
//                 &nb.xij[0], &nb.xij[nijatoms], &nb.xij[2*nijatoms], nijatoms);               
//         csh.RadialSphericalHarmonicsSum(ar, ai, sr, si, nb.neighnum, niatoms);    
//         
//         if (spectrum==0) {  // power spectrum    
//             csh.RadialSphericalHarmonicsPower(c, ar, ai, niatoms);               
//             csh.RadialSphericalHarmonicsPowerDeriv2(cd, ar, ai, srx, six, sry, siy, srz, siz, nb.neighnum, niatoms);
//             csh.RadialSphericalHarmonicsBasis(d, c, atomtype, natomtypes, niatoms, nbasis);                                
//             csh.RadialSphericalHarmonicsBasisDeriv2(dd, cd, atomtype, nb.neighlist, nb.neighnum, 
//                     natomtypes, niatoms, nbasis);            
//         }
//         else if (spectrum==1) { // bispectrum
//             csh.RadialSphericalHarmonicsBispectrum(c, ar, ai, niatoms);               
//             csh.RadialSphericalHarmonicsBispectrumDeriv2(cd, ar, ai, srx, six, sry, siy, srz, siz, nb.neighnum, niatoms);
//             csh.RadialSphericalHarmonicsBasis(d, c, atomtype, natomtypes, niatoms, nbasis);                   
//             csh.RadialSphericalHarmonicsBasisDeriv2(dd, cd, atomtype, nb.neighlist, nb.neighnum, 
//                     natomtypes, niatoms, nbasis);                        
//         }
//         else if (spectrum==2) { // power spectrum and bispectrum         
//             Int npower = (L+1)*K*(K+1)/2;      // number of power components
//             csh.RadialSphericalHarmonicsPower(c, ar, ai, niatoms);   
//             csh.RadialSphericalHarmonicsBispectrum(&c[niatoms*npower], ar, ai, niatoms);               
//             csh.RadialSphericalHarmonicsPowerDeriv2(cd, ar, ai, srx, six, sry, siy, srz, siz, nb.neighnum, niatoms);
//             csh.RadialSphericalHarmonicsBispectrumDeriv2(&cd[3*nijatoms*npower], ar, ai, srx, six, sry, siy, srz, siz, nb.neighnum, niatoms);
//             csh.RadialSphericalHarmonicsBasis(d, c, atomtype, natomtypes, niatoms, nbasis);    
//             csh.RadialSphericalHarmonicsBasisDeriv2(dd, cd, atomtype, nb.neighlist, nb.neighnum, 
//                     natomtypes, niatoms, nbasis);                        
//         }
//     }    
// }
// 
// void CDescriptors::Energy(dstype e, dstype* x, dstype* coeff, Int* atomtype, Int numatoms, Int ncoeff)
// {    
//     dstype *d = &tmp.tmpmem[0];
//     this->BasisFunctions(d, x, atomtype, numatoms);
//     DOT(cublasHandle, ncoeff, d, inc1, coeff, inc1, &e, backend);
// }
// 
// void CDescriptors::Forces(dstype* f, dstype* x, dstype* coeff, Int* atomtype, Int numatoms, Int ncoeff)
// {
//     dstype *dd = &tmp.tmpmem[0];
//     this->BasisFunctionsDeriv(dd, x, atomtype, numatoms);     
//     PGEMNV(cublasHandle, 3*numatoms, ncoeff, &one, dd, 3*numatoms, coeff, inc1, &zero, f, inc1, backend);        
// }
// 
// void CDescriptors::EnergyForces(dstype e, dstype* f, dstype* x, dstype* coeff, Int* atomtype, Int numatoms, Int ncoeff)
// {
//     dstype *d = &tmp.tmpmem[0];
//     dstype *dd = &tmp.tmpmem[ncoeff];
//     this->BasisFunctionsWithDeriv(d, dd, x, atomtype, numatoms);     
//     DOT(cublasHandle, ncoeff, d, inc1, coeff, inc1, &e, backend);
//     PGEMNV(cublasHandle, 3*numatoms, ncoeff, &one, dd, 3*numatoms, coeff, inc1, &zero, f, inc1, backend);            
// }
// 
// void CDescriptors::Stresses(dstype* s, dstype* x, dstype* coeff, Int* atomtype, Int numatoms, Int ncoeff)
// {
// }
//         
// void CDescriptors::EnergyForcesStresses(dstype e, dstype* f, dstype* s, dstype* x, dstype* coeff, Int* atomtype, Int numatoms, Int ncoeff)
// {
// }

#endif        

