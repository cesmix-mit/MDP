#ifndef __DESCRIPTORS_IMPL
#define __DESCRIPTORS_IMPL

// void implSetNeighborStruct(appstruct &app, commonstruct &common, configstruct &config, neighborstruct &nb)
// {
// //     Int backend = common.backend;
// //     Int nimax = common.nimax;
// //     Int njmax = common.njmax;
// //     
// //     nb.freememory(common.backend);
// //     
// //     TemplateMalloc(&nb.bcs, 6, common.backend); 
// //     TemplateMalloc(&nb.simbox, 9, common.backend); 
// //     TemplateMalloc(&nb.neighnum, nimax, common.backend); 
// //     //TemplateMalloc(&nb.atomtype, nimax, common.backend); 
// //     TemplateMalloc(&nb.neighlist, nimax*njmax, common.backend); 
// //     TemplateMalloc(&nb.xij, nimax*njmax*3, common.backend); 
// //     
// //     if (backend==2) { // GPU
// // #ifdef HAVE_CUDA        
// //         CHECK( cudaMemcpy(nb.bcs, config.bcs, 6*sizeof(Int), cudaMemcpyHostToDevice ) );              
// //         CHECK( cudaMemcpy(nb.simbox, config.simbox, 9*sizeof(dstype), cudaMemcpyHostToDevice ) );              
// // #endif        
// //     }    
// //     else { // CPU
// //         for (int i=0; i<6; i++)
// //             nb.bcs[i] = config.bcs[i];                
// //         for (int i=0; i<9; i++)
// //             nb.simbox[i] = config.simbox[i];                
// //     }        
// }

// void implSetTempStruct(appstruct &app, commonstruct &common, configstruct &config, tempstruct &tmp)
// {
//     
// }
// 
// void implSetSysStruct(appstruct &app, commonstruct &common, configstruct &config, sysstruct &sys)
// {
//     
// }

void implInit(appstruct &app, commonstruct &common, configstruct &config, neighborstruct &nb, tempstruct &tmp, sysstruct &sys, shstruct &sh)
{
    Int backend = common.backend;
    Int descriptor = common.descriptor;
    Int spectrum = common.spectrum;             
    Int natomtypes = common.natomtypes;   // number of atom types;    
    
    if (descriptor==0) { 
        Int K = sh.K;
        Int L = sh.L;
        Int Nub = sh.Nub;        
        if (spectrum==0) {  // power spectrum          
            sh.nbasis = (L+1)*K*(K+1)/2;      // number of power components              
            sh.npower = (L+1)*K*(K+1)/2;      // number of power components
        }
        else if (spectrum==1) { // bispectrum
            sh.nbasis = Nub*K*(K+1)/2;        // number of bispectrum components
            sh.nbispectrum = Nub*K*(K+1)/2;   // number of bispectrum components
        }
        else if (spectrum==2) { // power spectrum and bispectrum         
            sh.npower = (L+1)*K*(K+1)/2;      // number of power components
            sh.nbispectrum = Nub*K*(K+1)/2;   // number of bispectrum components
            sh.nbasis = sh.npower + sh.nbispectrum;
        }
        else {
        }        
    }
        
//     implSetNeighborStruct(app, common, config, nb);
//     implSetSysStruct(app, common, config, sys);    
//     implSetTempStruct(app, common, config, tmp);    
} 

// void implBuildVerletList(dstype* x, Int* atomtype, Int numatoms, Int backend, neighborstruct &nb, tempstruct &tmp)
// { 
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
// 
// void implNeighborList(dstype* x, Int* atomtype, Int numatoms, Int backend, neighborstruct &nb, tempstruct &tmp)
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
// int implNeighborPairs(dstype* x, Int* atomtype, Int *ilist, Int ilistsize, Int typei, Int typej, Int backend, neighborstruct &nb, tempstruct &tmp)
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
// 

// void implBasisFunctions(dstype *d, dstype *x, Int* atomtype, Int numatoms, Int backend, appstruct &app, 
//         commonstruct &common, configstruct &config, neighborstruct &nb, tempstruct &tmp, shstruct &sh)
// {       
//     //this->NeighborPairs(x, atomtype, ilist, listsize, typei, typej);    
//         
//     if (common.descriptor==0) { 
//         Int K = sh.K;
//         Int L = sh.L;
//         Int Nsh = K*(L+1)*(L+1);
//         Int nijatoms = nb.nxij;
//         Int niatoms = numatoms;
//         Int spectrum = common.spectrum;
//         Int nbasis = sh.nbasis;
//         Int natomtypes = common.natomtypes;
//         
//         dstype *sr = &tmp.tmpmem[0];
//         dstype *si = &tmp.tmpmem[nijatoms*Nsh];
//         dstype *ar = &tmp.tmpmem[(2*nijatoms)*Nsh];
//         dstype *ai = &tmp.tmpmem[(2*nijatoms+niatoms)*Nsh];        
//         dstype *c = &tmp.tmpmem[(2*nijatoms+2*niatoms)*Nsh];
//         
//         // xij, neighnum, atomtype, 
//         implSphericalHarmonicsBessel(sr, si, nb.xij, nijatoms, backend, sh);               
//         implRadialSphericalHarmonicsSum(ar, ai, sr, si, nb.neighnum, niatoms, backend, sh);    
//         
//         if (spectrum==0) {  // power spectrum          
//             implRadialSphericalHarmonicsPower(c, ar, ai, niatoms, backend, sh);   
//             implRadialSphericalHarmonicsBasis(d, c, atomtype, natomtypes, niatoms, nbasis, backend);                    
//         }
//         else if (spectrum==1) { // bispectrum
//             implRadialSphericalHarmonicsBispectrum(c, ar, ai, niatoms, backend, sh);   
//             implRadialSphericalHarmonicsBasis(d, c, atomtype, natomtypes, niatoms, nbasis, backend);                    
//         }
//         else if (spectrum==2) { // power spectrum and bispectrum         
//             Int npower = (L+1)*K*(K+1)/2;      // number of power components
//             implRadialSphericalHarmonicsPower(c, ar, ai, niatoms, backend, sh);   
//             implRadialSphericalHarmonicsBispectrum(&c[niatoms*npower], ar, ai, niatoms, backend, sh);   
//             implRadialSphericalHarmonicsBasis(d, c, atomtype, natomtypes, niatoms, nbasis, backend);                    
//         }
//     }
// }
// 
// void implBasisFunctionsDeriv(dstype* dd, dstype* x, Int* atomtype, Int numatoms, Int backend, appstruct &app, 
//         commonstruct &common, configstruct &config, neighborstruct &nb, tempstruct &tmp, shstruct &sh)
// {
//     //this->NeighborList(x, atomtype, numatoms);  
//     
//     if (common.descriptor==0) { 
//         Int K = sh.K;
//         Int L = sh.L;
//         Int Nsh = K*(L+1)*(L+1);
//         Int nijatoms = nb.nxij;
//         Int niatoms = numatoms;
//         Int spectrum = common.spectrum;
//         Int nbasis = sh.nbasis;
//         Int natomtypes = common.natomtypes;
// 
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
//         implSphericalHarmonicsBesselWithDeriv(sr, si, srx, six, sry, siy, srz, siz, nb.xij, nijatoms, backend, sh);               
//         implRadialSphericalHarmonicsSum(ar, ai, sr, si, nb.neighnum, niatoms, backend, sh);    
//         
//         if (spectrum==0) {  // power spectrum    
//             implRadialSphericalHarmonicsPowerDeriv2(cd, ar, ai, srx, six, sry, siy, srz, siz, 
//                     nb.neighnum, niatoms, backend, sh);
//             implRadialSphericalHarmonicsBasisDeriv2(dd, cd, atomtype, nb.neighlist, nb.neighnum, 
//                     natomtypes, niatoms, nbasis, backend);            
//         }
//         else if (spectrum==1) { // bispectrum
//             implRadialSphericalHarmonicsBispectrumDeriv2(cd, ar, ai, srx, six, sry, siy, srz, siz, 
//                     nb.neighnum, niatoms, backend, sh);
//             implRadialSphericalHarmonicsBasisDeriv2(dd, cd, atomtype, nb.neighlist, nb.neighnum, 
//                     natomtypes, niatoms, nbasis, backend);                        
//         }
//         else if (spectrum==2) { // power spectrum and bispectrum         
//             Int npower = (L+1)*K*(K+1)/2;      // number of power components
//             implRadialSphericalHarmonicsPowerDeriv2(cd, ar, ai, srx, six, sry, siy, srz, siz, 
//                     nb.neighnum, niatoms, backend, sh);
//             implRadialSphericalHarmonicsBispectrumDeriv2(&cd[3*nijatoms*npower], ar, ai, srx, six, 
//                     sry, siy, srz, siz, nb.neighnum, niatoms, backend, sh);
//             implRadialSphericalHarmonicsBasisDeriv2(dd, cd, atomtype, nb.neighlist, nb.neighnum, 
//                     natomtypes, niatoms, nbasis, backend);                        
//         }
//     }    
// }
// 
// void implBasisFunctionsWithDeriv(dstype* d, dstype* dd, dstype* x, Int* atomtype, Int numatoms, Int backend,
//    appstruct &app, commonstruct &common, configstruct &config, neighborstruct &nb, tempstruct &tmp, shstruct &sh)
// {
//     //NeighborList(x, atomtype, numatoms);  
//     
//     if (common.descriptor==0) { 
//         Int K = sh.K;
//         Int L = sh.L;
//         Int Nsh = K*(L+1)*(L+1);
//         Int nijatoms = nb.nxij;
//         Int niatoms = numatoms;
//         Int spectrum = common.spectrum;
//         Int nbasis = sh.nbasis;
//         Int natomtypes = common.natomtypes;        
//         
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
//         implSphericalHarmonicsBesselWithDeriv(sr, si, srx, six, sry, siy, srz, siz,
//                 &nb.xij[0], &nb.xij[nijatoms], &nb.xij[2*nijatoms], nijatoms, backend, sh);               
//         implRadialSphericalHarmonicsSum(ar, ai, sr, si, nb.neighnum, niatoms, backend, sh);    
//         
//         if (spectrum==0) {  // power spectrum    
//             implRadialSphericalHarmonicsPower(c, ar, ai, niatoms, backend, sh);               
//             implRadialSphericalHarmonicsPowerDeriv2(cd, ar, ai, srx, six, sry, siy, srz, siz, nb.neighnum, niatoms, backend, sh);
//             implRadialSphericalHarmonicsBasis(d, c, atomtype, natomtypes, niatoms, nbasis, backend);                                
//             implRadialSphericalHarmonicsBasisDeriv2(dd, cd, atomtype, nb.neighlist, nb.neighnum, 
//                     natomtypes, niatoms, nbasis, backend);            
//         }
//         else if (spectrum==1) { // bispectrum
//             implRadialSphericalHarmonicsBispectrum(c, ar, ai, niatoms, backend, sh);               
//             implRadialSphericalHarmonicsBispectrumDeriv2(cd, ar, ai, srx, six, sry, siy, srz, siz, nb.neighnum, niatoms, backend, sh);
//             implRadialSphericalHarmonicsBasis(d, c, atomtype, natomtypes, niatoms, nbasis, backend);                   
//             implRadialSphericalHarmonicsBasisDeriv2(dd, cd, atomtype, nb.neighlist, nb.neighnum, 
//                     natomtypes, niatoms, nbasis, backend);                        
//         }
//         else if (spectrum==2) { // power spectrum and bispectrum         
//             Int npower = (L+1)*K*(K+1)/2;      // number of power components
//             implRadialSphericalHarmonicsPower(c, ar, ai, niatoms, backend, sh);   
//             implRadialSphericalHarmonicsBispectrum(&c[niatoms*npower], ar, ai, niatoms, backend, sh);               
//             implRadialSphericalHarmonicsPowerDeriv2(cd, ar, ai, srx, six, sry, siy, srz, siz, nb.neighnum, niatoms, backend, sh);
//             implRadialSphericalHarmonicsBispectrumDeriv2(&cd[3*nijatoms*npower], ar, ai, srx, six, sry, siy, srz, siz, nb.neighnum, niatoms, backend, sh);
//             implRadialSphericalHarmonicsBasis(d, c, atomtype, natomtypes, niatoms, nbasis, backend);    
//             implRadialSphericalHarmonicsBasisDeriv2(dd, cd, atomtype, nb.neighlist, nb.neighnum, 
//                     natomtypes, niatoms, nbasis, backend);                        
//         }
//     }    
// }
// 

// void implEnergy(dstype e, dstype* x, dstype* coeff, Int* atomtype, Int numatoms, Int ncoeff, Int backend, neighborstruct &nb, tempstruct &tmp)
// {    
//     dstype *d = &tmp.tmpmem[0];
//     this->BasisFunctions(d, x, atomtype, numatoms);
//     DOT(cublasHandle, ncoeff, d, inc1, coeff, inc1, &e, backend);
// }
// 
// void implForces(dstype* f, dstype* x, dstype* coeff, Int* atomtype, Int numatoms, Int ncoeff, Int backend, neighborstruct &nb, tempstruct &tmp)
// {
//     dstype *dd = &tmp.tmpmem[0];
//     this->BasisFunctionsDeriv(dd, x, atomtype, numatoms);     
//     PGEMNV(cublasHandle, 3*numatoms, ncoeff, &one, dd, 3*numatoms, coeff, inc1, &zero, f, inc1, backend);        
// }
// 
// void implEnergyForces(dstype e, dstype* f, dstype* x, dstype* coeff, Int* atomtype, Int numatoms, Int ncoeff, Int backend, neighborstruct &nb, tempstruct &tmp)
// {
//     dstype *d = &tmp.tmpmem[0];
//     dstype *dd = &tmp.tmpmem[ncoeff];
//     this->BasisFunctionsWithDeriv(d, dd, x, atomtype, numatoms);     
//     DOT(cublasHandle, ncoeff, d, inc1, coeff, inc1, &e, backend);
//     PGEMNV(cublasHandle, 3*numatoms, ncoeff, &one, dd, 3*numatoms, coeff, inc1, &zero, f, inc1, backend);            
// }
// 
// void implStresses(dstype* s, dstype* x, dstype* coeff, Int* atomtype, Int numatoms, Int ncoeff, Int backend, neighborstruct &nb, tempstruct &tmp)
// {
// }
//         
// void implEnergyForcesStresses(dstype e, dstype* f, dstype* s, dstype* x, dstype* coeff, Int* atomtype, Int numatoms, Int ncoeff)
// {
// }

#endif        

