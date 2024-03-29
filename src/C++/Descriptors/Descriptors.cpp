/***************************************************************************                               
                    Molecular Dynamics Potentials (MDP)
                           CESMIX-MIT Project  
 
 Contributing authors: Ngoc-Cuong Nguyen (cuongng@mit.edu, exapde@gmail.com)
 ***************************************************************************/

#ifndef __DESCRIPTORS
#define __DESCRIPTORS

#include "snappotential.cpp"

void InitSphericalHarmonics(shstruct &sh, commonstruct &common) 
{
//     L       the maximum degree of spherical harmonics 
//     K       the number of zeros of spherical Bessel functions
//     Nub     number of non-zero unqiue bispectrum compoments of spherical harmonics 
//     Ncg     the total number of non-zero Clebsch-Gordan coefficients 
//     backend computing platform
    
    Int L = common.L;
    Int K = common.K;
    Int backend = common.backend;
    Int Nub;// number of non-zero unqiue bispectrum compoments of spherical harmonics 
    Int Ncg;// the total number of non-zero Clebsch-Gordan coefficients 
        
    Int M = (L+1)*(L+2)/2;
    Int L2 = (L+1)*(L+1);     
    Int K2 = K*(K+1);
    
    // factorial table
    dstype fac[168];    
    for (Int i=0; i<168; i++)
        fac[i] = (dstype) factable[i];
    
    // get index pairs (k2, k1) with k2 <= k1 
    Int indk[K2]; 
    cpuGetIndk(indk, K);
                
    dstype f[L+1];        
    dstype P[M];
    dstype tmp[M];
    dstype df[L+1];
    dstype dP[M];
    dstype dtmp[M];    
    
    // x, y, z
    Int N = 100;
    string filein = "spherecoords.bin";
    ifstream in(filein.c_str(), ios::in | ios::binary);    
    if (!in) 
        error("Unable to open file " + filein);           
    dstype *sph;        
    readarray(in, &sph, N*3);
    in.close();    
    dstype *the = &sph[0];
    dstype *phi = &sph[N];
    dstype *r = &sph[2*N];    
    
    dstype x[N];
    dstype y[N];
    dstype z[N];        
    cpuSphere2Cart(x, y, z, the, phi, r, N);
        
    dstype Ylmr[L2];
    dstype Ylmi[L2];            
    dstype b[L2*(L+1)];    
    cpuSphericalHarmonicsAtomSum(Ylmr, Ylmi, x, y, z, P, tmp, fac, M_PI, L, N);
    cpuSphericalHarmonicsBispectrum(b, Ylmr, Ylmi, fac, L);    
    
    // compute non-zero Clebsch-Gordan coefficients 
    Nub = cpuSphericalHarmonicsBispectrumIndex(b, L);
    Int indl[Nub*3];
    cpuSphericalHarmonicsBispectrumIndex(indl, b, Nub, L);    
    Ncg = cgcoefficients(indl, Nub);        
    Int indm[Ncg*3];
    Int rowm[Nub+1];    
    dstype cg[Ncg];                        
    cgcoefficients(cg, indm, rowm, indl, fac, Ncg, Nub);        
    
    // roots of spherical Bessel functions
    string filename = "besselzeros.bin";
    ifstream fin(filename.c_str(), ios::in | ios::binary);    
    if (!fin) 
        error("Unable to open file " + filename);
    dstype *bzeros;
    readarray(fin, &bzeros, 25*20);
    fin.close();    
    dstype x0[(L+1)*K];
    for (int k=0; k<K; k++)
        for (int l=0; l<(L+1); l++)
            x0[k*(L+1) + l] = bzeros[k*25 + l]/common.rcutml;
    
    
    //printArray2D(x0, (1+L), K, backend);    
    //cout<<common.rcutml<<endl;
    //sh.freememory(backend);
         
    // allocate memory for sh struct
    TemplateMalloc(&sh.indk, K2, backend); 
    TemplateMalloc(&sh.indl, Nub*3, backend); 
    TemplateMalloc(&sh.indm, Ncg*3, backend); 
    TemplateMalloc(&sh.rowm, Nub+1, backend); 
    TemplateMalloc(&sh.cg, Ncg, backend); 
    
    TemplateMalloc(&sh.fac, 168, backend); 
    TemplateMalloc(&sh.x0, (L+1)*K, backend); 
    TemplateMalloc(&sh.f, (L+1), backend); 
    TemplateMalloc(&sh.P, M, backend); 
    TemplateMalloc(&sh.tmp, M, backend); 
    TemplateMalloc(&sh.df, (L+1), backend); 
    TemplateMalloc(&sh.dP, M, backend); 
    TemplateMalloc(&sh.dtmp, M, backend); 
            
    // copy data to sh struct 
    if (backend==2) { // GPU
#ifdef HAVE_CUDA        
        CUDA_CHECK( cudaMemcpy(sh.indk, indk, K2*sizeof(Int), cudaMemcpyHostToDevice ) );      
        CUDA_CHECK( cudaMemcpy(sh.indl, indl, Nub*3*sizeof(Int), cudaMemcpyHostToDevice ) );      
        CUDA_CHECK( cudaMemcpy(sh.indm, indm, Ncg*3*sizeof(Int), cudaMemcpyHostToDevice ) );      
        CUDA_CHECK( cudaMemcpy(sh.rowm, rowm, (Nub+1)*sizeof(Int), cudaMemcpyHostToDevice ) );   
        CUDA_CHECK( cudaMemcpy(sh.x0, x0, ((L+1)*K)*sizeof(dstype), cudaMemcpyHostToDevice ) );      
        CUDA_CHECK( cudaMemcpy(sh.cg, cg, Ncg*sizeof(dstype), cudaMemcpyHostToDevice ) );      
        CUDA_CHECK( cudaMemcpy(sh.fac, fac, 168*sizeof(dstype), cudaMemcpyHostToDevice ) );              
#endif        
    }    
    else { // CPU
        for (int i=0; i<K2; i++)
            sh.indk[i] = indk[i];
        for (int i=0; i<Nub*3; i++)
            sh.indl[i] = indl[i];
        for (int i=0; i<Ncg*3; i++)
            sh.indm[i] = indm[i];
        for (int i=0; i<=Nub; i++)
            sh.rowm[i] = rowm[i];
        for (int i=0; i<(L+1)*K; i++)
            sh.x0[i] = x0[i];
        for (int i=0; i<Ncg; i++)
            sh.cg[i] = cg[i];
        for (int i=0; i<168; i++)
            sh.fac[i] = fac[i];                
    }    
    
    common.Nub = Nub;// number of non-zero unqiue bispectrum compoments of spherical harmonics 
    common.Ncg = Ncg;// the total number of non-zero Clebsch-Gordan coefficients     
    
    if (common.descriptor==0) { 
        Int Nub = common.Nub;        
        if (common.spectrum==0) {  // power spectrum          
            common.Nbf = (L+1)*K*(K+1)/2;      // number of power components              
            common.Npower = (L+1)*K*(K+1)/2;      // number of power components
            common.Nbispectrum = 0;
        }
        else if (common.spectrum==1) { // bispectrum
            common.Nbf = Nub*K*(K+1)/2;        // number of bispectrum components
            common.Npower = 0;
            common.Nbispectrum = Nub*K*(K+1)/2;   // number of bispectrum components
        }
        else if (common.spectrum==2) { // power spectrum and bispectrum         
            common.Npower = (L+1)*K*(K+1)/2;      // number of power components
            common.Nbispectrum = Nub*K*(K+1)/2;   // number of bispectrum components
            common.Nbf = common.Npower + common.Nbispectrum;
        }
        else {
        }        
    }            
    if (common.chemtype == 0)
        common.Ncoeff = common.Nbf*common.natomtypes;
    else
        common.Ncoeff = common.Nbf*common.natomtypes*common.natomtypes;    
    
    sh.L= common.L; 
    sh.K= common.K; 
    sh.Nub= common.Nub;
    sh.Ncg= common.Ncg;
    sh.npower= common.Npower;
    sh.nbispectrum= common.Nbispectrum;
    sh.nbasis= common.Nbf;    
    
    //cout<<common.L<<" "<<common.K<<" "<<common.Nub<<" "<<common.Ncg<<" "<<common.Npower<<" "<<common.Nbispectrum<<" "<<common.Nbf<<" "<<common.Ncoeff<<endl;
}

void SphericalHarmonicBesselDescriptors(dstype *e, neighborstruct &nb, commonstruct &common, appstruct &app, tempstruct &tmp, 
       shstruct &sh, dstype* x, dstype *q, dstype* param, dstype *rcutsq, Int *atomtype, Int nparam, Int typei, Int typej, Int decomp)
{        
    Int Nbf = common.Nbf;
    Int ncq = 0;
    for (Int b=0; b<common.nba; b++) {
        Int e1 = common.ablks[b];
        Int e2 = common.ablks[b+1];            
        Int na = e2 - e1; // number of atoms in this block
        Int neighmax = common.neighmax;
        Int dim = common.dim;
        Int backend = common.backend;
                
        Int *ilist = &tmp.intmem[0]; //na     
        if (typei>0) {               
            Int *olist = &tmp.intmem[na]; //na        
            ArrayFill(olist, e1, na, backend);        

            Int *t0 = &tmp.intmem[2*na]; //na        
            Int *t1 = &tmp.intmem[3*na]; //na        
            na = FindAtomType(ilist, olist, atomtype, t0, t1, typei, na, backend);
        }
        else {
            ArrayFill(ilist, e1, na, backend);        
        }                
        
        Int *pairnum = &tmp.intmem[na]; // na
        Int *pairlist = &tmp.intmem[2*na]; // na*neighmax
        if (typej>0)
            FullNeighPairList(pairnum, pairlist, x, rcutsq, atomtype, ilist, nb.alist, nb.neighlist, nb.neighnum, na, neighmax, typej, dim, backend);
        else
            FullNeighPairList(pairnum, pairlist, x, rcutsq, ilist, nb.neighlist, nb.neighnum, na, neighmax, dim, backend);        
                
        //a list contains the starting positions of the first neighbor 
        Int *pairnumsum = &tmp.intmem[2*na+na*neighmax]; // na+1        
        //cpuCumsum(pairnumsum, pairnum, na+1);                                         
        Cumsum(pairnumsum, pairnum, &tmp.intmem[3*na+na*neighmax+1], &tmp.intmem[4*na+na*neighmax+2], na+1, backend);                                         
        int ntuples = IntArrayGetValueAtIndex(pairnumsum, na, backend);                             
        
        Int *ai = &tmp.intmem[1+3*na+na*neighmax]; // ntuples        
        Int *aj = &tmp.intmem[1+3*na+ntuples+na*neighmax]; // ntuples        
        Int *ti = &tmp.intmem[1+3*na+2*ntuples+na*neighmax]; // ntuples        
        Int *tj = &tmp.intmem[1+3*na+3*ntuples+na*neighmax]; // ntuples        
        //dstype *xij = &tmp.tmpmem[0]; // ntuples*dim
        Int nbasis = common.K*(common.L+1)*(common.L+1);
        dstype *xij = &tmp.tmpmem[2*na*nbasis+ntuples*(2*nbasis)]; // ntuples*dim
        dstype *qi;
        dstype *qj;
        NeighPairs(xij, qi, qj, x, q, ai, aj, ti, tj, pairnum, pairlist, pairnumsum, ilist, nb.alist, 
                atomtype, na, neighmax, ncq, dim, backend);       
                                                
        dstype *cr =  &tmp.tmpmem[0];         // na*nbasis
        dstype *ci =  &tmp.tmpmem[na*nbasis]; // na*nbasis              
        dstype *sr =  &tmp.tmpmem[2*na*nbasis]; // ntuples*nbasis
        dstype *si =  &tmp.tmpmem[2*na*nbasis+ntuples*(1*nbasis)]; // ntuples*nbasis                       
        SphericalHarmonicsBessel(sr, si, xij, sh.x0, sh.P, sh.tmp, sh.f, 
                sh.fac, M_PI, common.L, common.K, ntuples, backend);
                                
        dstype *ei = &tmp.tmpmem[2*na*nbasis+ntuples*(2*nbasis)]; // na * Nbf
        RadialSphericalHarmonicsSpectrum(ei, cr, ci, sr, si, sh.cg, sh.indk, sh.indl, sh.indm, 
                sh.rowm, pairnumsum, common.Nub, common.Ncg, na, common.L, common.K, common.spectrum, backend);
        
        // ei : na x Nbf,  onevec : na x 1, e : Nbf x 1
        // e = e + ei^T*onevec        
        dstype *onevec =  &tmp.tmpmem[0];  // na
        ArraySetValue(onevec, 1.0, na, common.backend);
        PGEMTV(common.cublasHandle, na, Nbf, &one, ei, na, onevec, inc1, &one, e, inc1, common.backend);                        
    }        
}

void implSphericalHarmonicBesselDescriptors(dstype *e, neighborstruct &nb, commonstruct &common, 
        appstruct &app, tempstruct &tmp, shstruct &sh, dstype* x, dstype *q, dstype* param, Int nparam)
{    
    //ArraySetValue(e, 0.0, common.Ncoeff, common.backend);
   // ArraySetValue(f, 0.0, common.dim*common.inum*common.Ncoeff, common.backend);
    Int natomtypes = common.natomtypes;
    if (natomtypes==1) {
        SphericalHarmonicBesselDescriptors(e, nb, common, app, tmp, sh, x, q, param, &app.rcutsqml[0], 
            nb.atomtype, nparam, 0, 0, common.decomposition);                    
    }                 
    else {
        Int Nbf = common.Nbf;
        Int inum = common.inum;
        if (common.chemtype == 0) {
            for (int i = 0; i < natomtypes; i++)
                SphericalHarmonicBesselDescriptors(&e[i*Nbf], nb, common, app, tmp, sh, 
                        x, q, param, &app.rcutsqml[0], nb.atomtype, nparam, i+1, 0, common.decomposition);                                    
        }
        else {
            for (int i = 0; i < natomtypes; i++)
                for (int j = 0; j < natomtypes; j++)
                    SphericalHarmonicBesselDescriptors(&e[(j+i*natomtypes)*Nbf], nb, common, app, tmp, sh, x, q, 
                        param, &app.rcutsqml[0], nb.atomtype, nparam, i+1, j+1, common.decomposition);                                    
        }
    }
}

void SphericalHarmonicBesselDescriptors(dstype *e, dstype *f, neighborstruct &nb, commonstruct &common, appstruct &app, tempstruct &tmp, 
       shstruct &sh, dstype* x, dstype *q, dstype* param, dstype *rcutsq, Int *atomtype, Int nparam, Int typei, Int typej, Int decomp)
{        
    Int Nbf = common.Nbf;
    Int ncq = 0;    
    INIT_TIMING;
    for (Int b=0; b<common.nba; b++) {
        Int e1 = common.ablks[b];
        Int e2 = common.ablks[b+1];            
        Int na = e2 - e1; // number of atoms in this block
        Int neighmax = common.neighmax;
        Int dim = common.dim;
        Int backend = common.backend;
        
        START_TIMING;
        Int *ilist = &tmp.intmem[0]; //na     
        if (typei>0) {               
            Int *olist = &tmp.intmem[na]; //na        
            ArrayFill(olist, e1, na, backend);        

            Int *t0 = &tmp.intmem[2*na]; //na        
            Int *t1 = &tmp.intmem[3*na]; //na        
            na = FindAtomType(ilist, olist, atomtype, t0, t1, typei, na, backend);
        }
        else {
            ArrayFill(ilist, e1, na, backend);        
        }                
                
        Int *pairnum = &tmp.intmem[na]; // na
        Int *pairlist = &tmp.intmem[2*na]; // na*neighmax
        if (typej>0)
            FullNeighPairList(pairnum, pairlist, x, rcutsq, atomtype, ilist, nb.alist, nb.neighlist, nb.neighnum, na, neighmax, typej, dim, backend);
        else
            FullNeighPairList(pairnum, pairlist, x, rcutsq, ilist, nb.neighlist, nb.neighnum, na, neighmax, dim, backend);        
                
                        //printArray2D(nb.neighnum, 1, na, common.backend);        
        //a list contains the starting positions of the first neighbor 
        Int *pairnumsum = &tmp.intmem[2*na+na*neighmax]; // na+1        
        //cpuCumsum(pairnumsum, pairnum, na+1);         
        Cumsum(pairnumsum, pairnum, &tmp.intmem[3*na+na*neighmax+1], &tmp.intmem[4*na+na*neighmax+2], na+1, backend);                                         
        int ntuples = IntArrayGetValueAtIndex(pairnumsum, na, common.backend);     
        END_TIMING(50);
        
        START_TIMING;
        Int *ai = &tmp.intmem[1+3*na+na*neighmax]; // ntuples        
        Int *aj = &tmp.intmem[1+3*na+ntuples+na*neighmax]; // ntuples        
        Int *ti = &tmp.intmem[1+3*na+2*ntuples+na*neighmax]; // ntuples        
        Int *tj = &tmp.intmem[1+3*na+3*ntuples+na*neighmax]; // ntuples        
        //dstype *xij = &tmp.tmpmem[0]; // ntuples*dim
        Int nbasis = common.K*(common.L+1)*(common.L+1);
        dstype *xij = &tmp.tmpmem[2*na*nbasis+ntuples*(8*nbasis)]; // ntuples*dim
        dstype *qi;
        dstype *qj;
        NeighPairs(xij, qi, qj, x, q, ai, aj, ti, tj, pairnum, pairlist, pairnumsum, ilist, nb.alist, 
                atomtype, na, neighmax, ncq, dim, backend);                                                               
        END_TIMING(51);
        
        START_TIMING;                
        dstype *cr =  &tmp.tmpmem[0];         // na*nbasis
        dstype *ci =  &tmp.tmpmem[na*nbasis]; // na*nbasis                        
        dstype *srx = &tmp.tmpmem[2*na*nbasis]; // ntuples*nbasis
        dstype *sry = &tmp.tmpmem[2*na*nbasis+ntuples*(1*nbasis)]; // ntuples*nbasis
        dstype *srz = &tmp.tmpmem[2*na*nbasis+ntuples*(2*nbasis)]; // ntuples*nbasis
        dstype *six = &tmp.tmpmem[2*na*nbasis+ntuples*(3*nbasis)]; // ntuples*nbasis
        dstype *siy = &tmp.tmpmem[2*na*nbasis+ntuples*(4*nbasis)]; // ntuples*nbasis
        dstype *siz = &tmp.tmpmem[2*na*nbasis+ntuples*(5*nbasis)]; // ntuples*nbasis         
        dstype *sr =  &tmp.tmpmem[2*na*nbasis+ntuples*(6*nbasis)]; // ntuples*nbasis
        dstype *si =  &tmp.tmpmem[2*na*nbasis+ntuples*(7*nbasis)]; // ntuples*nbasis               
        SphericalHarmonicsBesselWithDeriv(sr, si, srx, six, sry, siy, srz, siz, xij,
                 sh.x0, sh.P, sh.tmp, sh.f, sh.dP, sh.dtmp, sh.df, sh.fac, M_PI, common.L, common.K, ntuples, backend);
        END_TIMING(52);
        
//         string fn0 = (common.backend == 2) ? "srgpu.bin" : "srcpu.bin";
//         writearray2file(fn0, sr, ntuples*nbasis, common.backend); 
//         fn0 = (common.backend == 2) ? "sigpu.bin" : "sicpu.bin";
//         writearray2file(fn0, si, ntuples*nbasis, common.backend); 
//         fn0 = (common.backend == 2) ? "xijgpu.bin" : "xijcpu.bin";
//         writearray2file(fn0, xij, ntuples*dim, common.backend); 
        
// #ifdef HAVE_CHECK  // check              
//         dstype epsil = 1e-6;
//         dstype *sr1 =  new dstype[ntuples*nbasis];
//         dstype *si1 =  new dstype[ntuples*nbasis]; 
//         dstype *xi1 =  new dstype[ntuples*dim]; 
//         dstype *xii =  new dstype[ntuples*dim]; 
//         dstype *ei1 =  new dstype[na*Nbf]; 
//         dstype *onev =  new dstype[na]; 
//         dstype *et =  new dstype[Nbf]; 
//         cpuArrayCopy(xii, xij, dim*ntuples);
//         cpuArrayCopy(xi1, xij, dim*ntuples);
//         ArraySetValue(onev, 1.0, na, common.backend);
//         
//         cpuArrayRowkAXPB(xi1, xij, 1.0, epsil, dim, ntuples, 0);
//         cpuSphericalHarmonicsBessel(sr1, si1, xi1, sh.x0, sh.P, sh.tmp, sh.f, 
//                  sh.fac, M_PI, common.L, common.K, ntuples);
//         cpuArrayAXPBY(sr1, sr1, sr, 1.0/epsil, -1.0/epsil, ntuples*nbasis);
//         cpuArrayAXPBY(sr1, sr1, srx, 1.0, -1.0, ntuples*nbasis);
//         cpuArrayAbs(sr1, sr1, ntuples*nbasis);        
//         cout<<"Maximum absolute error: "<<cpuArrayMax(sr1, ntuples*nbasis)<<endl;        
//         
//         cpuArrayCopy(xi1, xij, dim*ntuples);
//         cpuArrayRowkAXPB(xi1, xij, 1.0, epsil, dim, ntuples, 1);
//         cpuSphericalHarmonicsBessel(sr1, si1, xi1, sh.x0, sh.P, sh.tmp, sh.f, 
//                  sh.fac, M_PI, common.L, common.K, ntuples);
//         cpuArrayAXPBY(sr1, sr1, sr, 1.0/epsil, -1.0/epsil, ntuples*nbasis);
//         cpuArrayAXPBY(sr1, sr1, sry, 1.0, -1.0, ntuples*nbasis);
//         cpuArrayAbs(sr1, sr1, ntuples*nbasis);        
//         cout<<"Maximum absolute error: "<<cpuArrayMax(sr1, ntuples*nbasis)<<endl;         
//         
//         cpuArrayCopy(xi1, xij, dim*ntuples);
//         cpuArrayRowkAXPB(xi1, xij, 1.0, epsil, dim, ntuples, 2);
//         cpuSphericalHarmonicsBessel(sr1, si1, xi1, sh.x0, sh.P, sh.tmp, sh.f, 
//                  sh.fac, M_PI, common.L, common.K, ntuples);
//         cpuArrayAXPBY(sr1, sr1, sr, 1.0/epsil, -1.0/epsil, ntuples*nbasis);
//         cpuArrayAXPBY(sr1, sr1, srz, 1.0, -1.0, ntuples*nbasis);
//         cpuArrayAbs(sr1, sr1, ntuples*nbasis);        
//         cout<<"Maximum absolute error: "<<cpuArrayMax(sr1, ntuples*nbasis)<<endl;    
// #endif                                      
        
        START_TIMING;
        dstype *ei = &tmp.tmpmem[2*na*nbasis+ntuples*(8*nbasis)]; // na * Nbf
        RadialSphericalHarmonicsSpectrum(ei, cr, ci, sr, si, sh.cg, sh.indk, sh.indl, sh.indm, 
                sh.rowm, pairnumsum, common.Nub, common.Ncg, na, common.L, common.K, common.spectrum, backend);
        
        // ei : na x Nbf,  onevec : na x 1, e : Nbf x 1
        // e = e + ei^T*onevec        
        dstype *onevec =  &tmp.tmpmem[2*na*nbasis+ntuples*(6*nbasis)];  // na
        ArraySetValue(onevec, one, na, common.backend);
        PGEMTV2(common.cublasHandle, na, Nbf, &one, ei, na, onevec, inc1, &one, e, inc1, common.backend);        
        END_TIMING(53);
             
//         string fn1 = (common.backend == 2) ? "eigpu.bin" : "eicpu.bin";
//         writearray2file(fn1, ei, na*Nbf, common.backend);              
//         fn1 = (common.backend == 2) ? "egpu.bin" : "ecpu.bin";
//         writearray2file(fn1, e, Nbf, common.backend); 
        //fn1 = (common.backend == 2) ? "ogpu.bin" : "ocpu.bin";
        //writearray2file(fn1, onevec, na, common.backend); 
        
        START_TIMING;
        dstype *fij = &tmp.tmpmem[2*na*nbasis+ntuples*(6*nbasis)]; // dim * ntuples * Nbf
//         RadialSphericalHarmonicsSpectrumDeriv(fij, cr, ci, srx, six, sry, siy, srz, siz, sh.cg, sh.indk, 
//             sh.indl, sh.indm, sh.rowm, pairnumsum, common.Nub, common.Ncg, na, common.L, common.K, common.spectrum, backend);                        
        RadialSphericalHarmonicsSpectrumDeriv(fij, cr, ci, srx, six, sry, siy, srz, siz, sh.cg, sh.indk, 
            sh.indl, sh.indm, sh.rowm, ai, pairnumsum, common.Nub, common.Ncg, na, ntuples, common.L, common.K, common.spectrum, backend);                        
        END_TIMING(54);
                
//         fn1 = (common.backend == 2) ? "fijgpu.bin" : "fijcpu.bin";
//         writearray2file(fn1, fij, dim*ntuples*Nbf, common.backend); 
//         error("here");

// #ifdef HAVE_CHECK    // check        
//         for (int i = 0; i<dim*ntuples; i++) {
//             cpuArrayCopy(xi1, xii, dim*ntuples);  
//             xi1[i] = xi1[i] + epsil;
//             cpuSphericalHarmonicsBessel(sr1, si1, xi1, sh.x0, sh.P, sh.tmp, sh.f, 
//                      sh.fac, M_PI, common.L, common.K, ntuples);
//                                  
//             cpuRadialSphericalHarmonicsSpectrum(ei1, cr, ci, sr1, si1, sh.cg, sh.indk, sh.indl, sh.indm, 
//                     sh.rowm, pairnumsum, common.Nub, common.Ncg, na, common.L, common.K, common.spectrum);
//             
//             PGEMTV(common.cublasHandle, na, Nbf, &one, ei1, na, onev, inc1, &zero, et, inc1, common.backend);                    
//             printArray2D(e, 1, Nbf, common.backend);
//             printArray2D(et, 1, Nbf, common.backend);
//             for (int j=0; j<Nbf; j++) 
//                 et[j] = fabs((et[j]-e[j])/epsil - fij[i + dim*ntuples*j]);
//             cout<<"Maximum absolute error: "<<cpuArrayMax(et, Nbf)<<endl;      
//         }
//         delete[] sr1;
//         delete[] si1;
//         delete[] xi1;
//         delete[] xii;
//         delete[] ei1;
//         delete[] onev;
//         delete[] et;        
//         error("here");            
// #endif
                
        START_TIMING;
        if (decomp==0)
            ForceDecomposition(f, fij, ai, aj, common.inum, ntuples, Nbf, backend);
        else {
            CenterAtomDecomposition(f, fij, ilist, pairnumsum, common.inum, ntuples, na, Nbf, backend);
            ArrayCopy(tmp.intmem, aj, ntuples, backend);
            Int *jlist = &tmp.intmem[ntuples];   
            Int *bnumsum = &tmp.intmem[2*ntuples]; 
            Int *index = &tmp.intmem[3*ntuples]; // ntuples       
            Int *p0 = &tmp.intmem[4*ntuples]; // ntuples       
            Int *p1 = &tmp.intmem[5*ntuples]; // ntuples       
            Int *p2 = &tmp.intmem[6*ntuples]; // ntuples       
            Int *p3 = &tmp.intmem[7*ntuples]; // ntuples       
            Int naj = UniqueSort(jlist, bnumsum, index, p0, tmp.intmem, p1, p2, p3, ntuples, backend);                        
            NeighborAtomDecomposition(f, fij, jlist, bnumsum, index, common.inum, ntuples, naj, Nbf, backend);
        }              
        END_TIMING(55);
    }     
    
#ifdef HAVE_DEBUG                      
    string fn = (common.backend == 2) ? "dgpu.bin" : "dcpu.bin";
    writearray2file(fn, e, Nbf, common.backend); 
    fn = (common.backend == 2) ? "ddgpu.bin" : "ddcpu.bin";
    writearray2file(fn, f, common.inum*Nbf, common.backend);
    //error("here");    
#endif                                

}

void implSphericalHarmonicBesselDescriptors(dstype *e, dstype *f, neighborstruct &nb, commonstruct &common, 
        appstruct &app, tempstruct &tmp, shstruct &sh, dstype* x, dstype *q, dstype* param, Int nparam)
{    
    //ArraySetValue(e, 0.0, common.Ncoeff, common.backend);
   // ArraySetValue(f, 0.0, common.dim*common.inum*common.Ncoeff, common.backend);
    Int natomtypes = common.natomtypes;
    if (natomtypes==1) {
        SphericalHarmonicBesselDescriptors(e, f, nb, common, app, tmp, sh, x, q, param, &app.rcutsqml[0], 
            nb.atomtype, nparam, 0, 0, common.decomposition);                    
    }                 
    else {
        Int Nbf = common.Nbf;
        Int inum = common.inum;
        if (common.chemtype == 0) {
            for (int i = 0; i < natomtypes; i++)
                SphericalHarmonicBesselDescriptors(&e[i*Nbf], &f[i*3*inum*Nbf], nb, common, app, tmp, sh, 
                        x, q, param, &app.rcutsqml[0], nb.atomtype, nparam, i+1, 0, common.decomposition);                                    
        }
        else {
            for (int i = 0; i < natomtypes; i++)
                for (int j = 0; j < natomtypes; j++)
                    SphericalHarmonicBesselDescriptors(&e[(j+i*natomtypes)*Nbf], &f[(j+i*natomtypes)*3*inum*Nbf], nb, common, app, tmp, sh, x, q, 
                        param, &app.rcutsqml[0], nb.atomtype, nparam, i+1, j+1, common.decomposition);                                    
        }
    }
}
 
void SphericalHarmonicBesselEnergyForce(dstype *e, dstype *f, neighborstruct &nb, commonstruct &common, appstruct &app, tempstruct &tmp, 
       shstruct &sh, dstype* x, dstype *coeff, dstype *q, dstype* param, dstype *rcutsq, Int *atomtype, Int nparam, Int typei, Int typej, Int decomp)
{        
    Int Nbf = common.Nbf;
    Int ncq = 0;
    for (Int b=0; b<common.nba; b++) {
        Int e1 = common.ablks[b];
        Int e2 = common.ablks[b+1];            
        Int na = e2 - e1; // number of atoms in this block
        Int neighmax = common.neighmax;
        Int dim = common.dim;
        Int backend = common.backend;
                
        Int *ilist = &tmp.intmem[0]; //na     
        if (typei>0) {               
            Int *olist = &tmp.intmem[na]; //na        
            ArrayFill(olist, e1, na, backend);        

            Int *t0 = &tmp.intmem[2*na]; //na        
            Int *t1 = &tmp.intmem[3*na]; //na        
            na = FindAtomType(ilist, olist, atomtype, t0, t1, typei, na, backend);
        }
        else {
            ArrayFill(ilist, e1, na, backend);        
        }                
                        
        Int *pairnum = &tmp.intmem[na]; // na
        Int *pairlist = &tmp.intmem[2*na]; // na*neighmax
        if (typej>0)
            FullNeighPairList(pairnum, pairlist, x, rcutsq, atomtype, ilist, nb.alist, nb.neighlist, nb.neighnum, na, neighmax, typej, dim, backend);
        else
            FullNeighPairList(pairnum, pairlist, x, rcutsq, ilist, nb.neighlist, nb.neighnum, na, neighmax, dim, backend);        
                                                
        //a list contains the starting positions of the first neighbor 
        Int *pairnumsum = &tmp.intmem[2*na+na*neighmax]; // na+1        
        Cumsum(pairnumsum, pairnum, &tmp.intmem[3*na+na*neighmax+1], &tmp.intmem[4*na+na*neighmax+2], na+1, backend);                                         
        int ntuples = IntArrayGetValueAtIndex(pairnumsum, na, common.backend);     
        
        Int *ai = &tmp.intmem[1+3*na+na*neighmax]; // ntuples        
        Int *aj = &tmp.intmem[1+3*na+ntuples+na*neighmax]; // ntuples        
        Int *ti = &tmp.intmem[1+3*na+2*ntuples+na*neighmax]; // ntuples        
        Int *tj = &tmp.intmem[1+3*na+3*ntuples+na*neighmax]; // ntuples        
        //dstype *xij = &tmp.tmpmem[0]; // ntuples*dim
        Int nbasis = common.K*(common.L+1)*(common.L+1);
        dstype *xij = &tmp.tmpmem[2*na*nbasis+ntuples*(8*nbasis)]; // ntuples*dim
        dstype *qi;
        dstype *qj;
        NeighPairs(xij, qi, qj, x, q, ai, aj, ti, tj, pairnum, pairlist, pairnumsum, ilist, nb.alist, 
                atomtype, na, neighmax, ncq, dim, backend);       
        
// #ifdef HAVE_CHECK         
//         dstype *xi1 =  new dstype[ntuples*dim]; 
//         dstype *xii =  new dstype[ntuples*dim];         
//         cpuArrayCopy(xii, xij, dim*ntuples);
//         cpuArrayCopy(xi1, xij, dim*ntuples);
// #endif
        
        dstype *cr =  &tmp.tmpmem[0];         // na*nbasis
        dstype *ci =  &tmp.tmpmem[na*nbasis]; // na*nbasis                        
        dstype *srx = &tmp.tmpmem[2*na*nbasis]; // ntuples*nbasis
        dstype *sry = &tmp.tmpmem[2*na*nbasis+ntuples*(1*nbasis)]; // ntuples*nbasis
        dstype *srz = &tmp.tmpmem[2*na*nbasis+ntuples*(2*nbasis)]; // ntuples*nbasis
        dstype *six = &tmp.tmpmem[2*na*nbasis+ntuples*(3*nbasis)]; // ntuples*nbasis
        dstype *siy = &tmp.tmpmem[2*na*nbasis+ntuples*(4*nbasis)]; // ntuples*nbasis
        dstype *siz = &tmp.tmpmem[2*na*nbasis+ntuples*(5*nbasis)]; // ntuples*nbasis         
        dstype *sr =  &tmp.tmpmem[2*na*nbasis+ntuples*(6*nbasis)]; // ntuples*nbasis
        dstype *si =  &tmp.tmpmem[2*na*nbasis+ntuples*(7*nbasis)]; // ntuples*nbasis                       
        SphericalHarmonicsBesselWithDeriv(sr, si, srx, six, sry, siy, srz, siz, xij,
                 sh.x0, sh.P, sh.tmp, sh.f, sh.dP, sh.dtmp, sh.df, sh.fac, M_PI, common.L, common.K, ntuples, backend);
                        
        dstype *di = &tmp.tmpmem[2*na*nbasis+ntuples*(8*nbasis)]; // na * Nbf                
        RadialSphericalHarmonicsSpectrum(di, cr, ci, sr, si, sh.cg, sh.indk, 
                sh.indl, sh.indm, sh.rowm, pairnumsum, common.Nub, common.Ncg, na, common.L, common.K, common.spectrum, backend);

        // di : na x Nbf,  coeff : Nbf x 1, ei : na x 1
        // ei = di*coeff                        
        dstype *ei = &tmp.tmpmem[2*na*nbasis+ntuples*(6*nbasis)]; // na         
        PGEMNV(common.cublasHandle, na, Nbf, &one, di, na, coeff, inc1, &zero, ei, inc1, common.backend);                
        CenterAtomDecomposition(e, ei, ilist, na, backend);
        
        dstype *dd = &tmp.tmpmem[2*na*nbasis+ntuples*(6*nbasis)]; // dim * ntuples * Nbf
        RadialSphericalHarmonicsSpectrumDeriv(dd, cr, ci, srx, six, sry, siy, srz, siz, sh.cg, sh.indk, 
                sh.indl, sh.indm, sh.rowm, pairnumsum, common.Nub, common.Ncg, na, common.L, common.K, common.spectrum, backend);                
        
        dstype *fij = &tmp.tmpmem[0]; // dim * ntuples 
        PGEMNV(common.cublasHandle, dim*ntuples, Nbf, &minusone, dd, dim*ntuples, coeff, inc1, &zero, fij, inc1, common.backend);    

// #ifdef HAVE_CHECK    
//         dstype en = cpuArraySum(e, na);
//         dstype epsil = 1e-6;        
//         dstype *sr1 =  new dstype[ntuples*nbasis];
//         dstype *si1 =  new dstype[ntuples*nbasis]; 
//         dstype *cr1 =  new dstype[na*nbasis]; 
//         dstype *ci1 =  new dstype[na*nbasis]; 
//         dstype *ei1 =  new dstype[na*Nbf]; 
//         dstype *onev =  new dstype[na]; 
//         dstype *et =  new dstype[Nbf]; 
//         ArraySetValue(onev, 1.0, na, common.backend);                
//         for (int i = 0; i<dim*ntuples; i++) {
//             cpuArrayCopy(xi1, xii, dim*ntuples);  
//             xi1[i] = xi1[i] + epsil;
//             
//             cpuSphericalHarmonicsBessel(sr1, si1, xi1, sh.x0, sh.P, sh.tmp, sh.f, 
//                      sh.fac, M_PI, common.L, common.K, ntuples);
//                                                          
//             cpuRadialSphericalHarmonicsSpectrum(ei1, cr1, ci1, sr1, si1, sh.cg, sh.indk, sh.indl, sh.indm, 
//                     sh.rowm, pairnumsum, common.Nub, common.Ncg, na, common.L, common.K, common.spectrum);
//                         
//             PGEMTV(common.cublasHandle, na, Nbf, &one, ei1, na, onev, inc1, &zero, et, inc1, common.backend);                    
// 
//             dstype em = 0.0;
//             for (int j=0; j<Nbf; j++)
//                 em += et[j]*coeff[j];                                    
//             dstype emax = fabs((en-em)/epsil - fij[i]);   
//                         
//             cout<<"Maximum absolute error: "<<emax<<endl;      
//         }
//         delete[] sr1;
//         delete[] si1;
//         delete[] cr1;
//         delete[] ci1;
//         delete[] xi1;
//         delete[] xii;
//         delete[] ei1;
//         delete[] onev;
//         delete[] et;        
//         error("here");            
// #endif
                        
        if (decomp==0)
            ForceDecomposition(f, fij, ai, aj, ntuples, backend);        
        else {
            CenterAtomDecomposition(f, fij, ilist, pairnumsum, na, backend);            
            ArrayCopy(tmp.intmem, aj, ntuples, backend);
            Int *jlist = &tmp.intmem[ntuples];   
            Int *bnumsum = &tmp.intmem[2*ntuples]; 
            Int *index = &tmp.intmem[3*ntuples]; // ntuples       
            Int *p0 = &tmp.intmem[4*ntuples]; // ntuples       
            Int *p1 = &tmp.intmem[5*ntuples]; // ntuples       
            Int *p2 = &tmp.intmem[6*ntuples]; // ntuples       
            Int *p3 = &tmp.intmem[7*ntuples]; // ntuples       
            Int naj = UniqueSort(jlist, bnumsum, index, p0, tmp.intmem, p1, p2, p3, ntuples, backend);                        
            NeighborAtomDecomposition(f, fij, jlist, bnumsum, index, naj, backend);            
        }
    }        
}

void implSphericalHarmonicBesselEnergyForce(dstype *e, dstype *f, neighborstruct &nb, commonstruct &common, 
        appstruct &app, tempstruct &tmp, shstruct &sh, dstype* x, dstype *coeff, dstype *q, dstype* param, Int nparam)
{    
    //ArraySetValue(e, 0.0, common.inum, common.backend);
    //ArraySetValue(f, 0.0, common.dim*common.inum, common.backend);
    Int natomtypes = common.natomtypes;
    if (natomtypes==1) {
        SphericalHarmonicBesselEnergyForce(e, f, nb, common, app, tmp, sh, x, coeff, q, param, &app.rcutsqml[0], 
            nb.atomtype, nparam, 0, 0, common.decomposition);                    
    }                 
    else {
        Int Nbf = common.Nbf;
        Int inum = common.inum;
        if (common.chemtype == 0) {
            for (int i = 0; i < natomtypes; i++)
                SphericalHarmonicBesselEnergyForce(e, f, nb, common, app, tmp, sh, 
                        x, &coeff[i*Nbf], q, param, &app.rcutsqml[0], nb.atomtype, nparam, i+1, 0, common.decomposition);                                    
        }
        else {
            for (int i = 0; i < natomtypes; i++)
                for (int j = 0; j < natomtypes; j++)
                    SphericalHarmonicBesselEnergyForce(e, f, nb, common, app, tmp, sh, 
                       x, &coeff[(j+i*natomtypes)*Nbf], q, param, &app.rcutsqml[0], nb.atomtype, nparam, i+1, j+1, common.decomposition);                                    
        }
    }
}

void SphericalHarmonicBesselEnergyForceVirial(dstype *e, dstype *f, dstype *v, neighborstruct &nb, commonstruct &common, appstruct &app, tempstruct &tmp, 
       shstruct &sh, dstype* x, dstype *coeff, dstype *q, dstype* param, dstype *rcutsq, Int *atomtype, Int nparam, Int typei, Int typej, Int decomp)
{        
    Int Nbf = common.Nbf;
    Int ncq = 0;
    for (Int b=0; b<common.nba; b++) {
        Int e1 = common.ablks[b];
        Int e2 = common.ablks[b+1];            
        Int na = e2 - e1; // number of atoms in this block
        Int neighmax = common.neighmax;
        Int dim = common.dim;
        Int backend = common.backend;
                
        Int *ilist = &tmp.intmem[0]; //na     
        if (typei>0) {               
            Int *olist = &tmp.intmem[na]; //na        
            ArrayFill(olist, e1, na, backend);        

            Int *t0 = &tmp.intmem[2*na]; //na        
            Int *t1 = &tmp.intmem[3*na]; //na        
            na = FindAtomType(ilist, olist, atomtype, t0, t1, typei, na, backend);
        }
        else {
            ArrayFill(ilist, e1, na, backend);        
        }                
                        
        Int *pairnum = &tmp.intmem[na]; // na
        Int *pairlist = &tmp.intmem[2*na]; // na*neighmax
        if (typej>0)
            FullNeighPairList(pairnum, pairlist, x, rcutsq, atomtype, ilist, nb.alist, nb.neighlist, nb.neighnum, na, neighmax, typej, dim, backend);
        else
            FullNeighPairList(pairnum, pairlist, x, rcutsq, ilist, nb.neighlist, nb.neighnum, na, neighmax, dim, backend);        
                                                
        //a list contains the starting positions of the first neighbor 
        Int *pairnumsum = &tmp.intmem[2*na+na*neighmax]; // na+1        
        Cumsum(pairnumsum, pairnum, &tmp.intmem[3*na+na*neighmax+1], &tmp.intmem[4*na+na*neighmax+2], na+1, backend);                                         
        int ntuples = IntArrayGetValueAtIndex(pairnumsum, na, common.backend);     
        
        Int *ai = &tmp.intmem[1+3*na+na*neighmax]; // ntuples        
        Int *aj = &tmp.intmem[1+3*na+ntuples+na*neighmax]; // ntuples        
        Int *ti = &tmp.intmem[1+3*na+2*ntuples+na*neighmax]; // ntuples        
        Int *tj = &tmp.intmem[1+3*na+3*ntuples+na*neighmax]; // ntuples        
        //dstype *xij = &tmp.tmpmem[0]; // ntuples*dim
        Int nbasis = common.K*(common.L+1)*(common.L+1);
        dstype *xij = &tmp.tmpmem[2*na*nbasis+ntuples*(8*nbasis)]; // ntuples*dim
        dstype *qi;
        dstype *qj;
        NeighPairs(xij, qi, qj, x, q, ai, aj, ti, tj, pairnum, pairlist, pairnumsum, ilist, nb.alist, 
                atomtype, na, neighmax, ncq, dim, backend);       
        
// #ifdef HAVE_CHECK         
//         dstype *xi1 =  new dstype[ntuples*dim]; 
//         dstype *xii =  new dstype[ntuples*dim];         
//         cpuArrayCopy(xii, xij, dim*ntuples);
//         cpuArrayCopy(xi1, xij, dim*ntuples);
// #endif
        
        dstype *cr =  &tmp.tmpmem[0];         // na*nbasis
        dstype *ci =  &tmp.tmpmem[na*nbasis]; // na*nbasis                        
        dstype *srx = &tmp.tmpmem[2*na*nbasis]; // ntuples*nbasis
        dstype *sry = &tmp.tmpmem[2*na*nbasis+ntuples*(1*nbasis)]; // ntuples*nbasis
        dstype *srz = &tmp.tmpmem[2*na*nbasis+ntuples*(2*nbasis)]; // ntuples*nbasis
        dstype *six = &tmp.tmpmem[2*na*nbasis+ntuples*(3*nbasis)]; // ntuples*nbasis
        dstype *siy = &tmp.tmpmem[2*na*nbasis+ntuples*(4*nbasis)]; // ntuples*nbasis
        dstype *siz = &tmp.tmpmem[2*na*nbasis+ntuples*(5*nbasis)]; // ntuples*nbasis         
        dstype *sr =  &tmp.tmpmem[2*na*nbasis+ntuples*(6*nbasis)]; // ntuples*nbasis
        dstype *si =  &tmp.tmpmem[2*na*nbasis+ntuples*(7*nbasis)]; // ntuples*nbasis                       
        SphericalHarmonicsBesselWithDeriv(sr, si, srx, six, sry, siy, srz, siz, xij,
                 sh.x0, sh.P, sh.tmp, sh.f, sh.dP, sh.dtmp, sh.df, sh.fac, M_PI, common.L, common.K, ntuples, backend);
                        
        dstype *di = &tmp.tmpmem[2*na*nbasis+ntuples*(8*nbasis)]; // na * Nbf                
        RadialSphericalHarmonicsSpectrum(di, cr, ci, sr, si, sh.cg, sh.indk, 
                sh.indl, sh.indm, sh.rowm, pairnumsum, common.Nub, common.Ncg, na, common.L, common.K, common.spectrum, backend);

        // di : na x Nbf,  coeff : Nbf x 1, ei : na x 1
        // ei = di*coeff                        
        dstype *ei = &tmp.tmpmem[2*na*nbasis+ntuples*(6*nbasis)]; // na         
        PGEMNV(common.cublasHandle, na, Nbf, &one, di, na, coeff, inc1, &zero, ei, inc1, common.backend);                
        CenterAtomDecomposition(e, ei, ilist, na, backend);
        
        dstype *dd = &tmp.tmpmem[2*na*nbasis+ntuples*(6*nbasis)]; // dim * ntuples * Nbf
        RadialSphericalHarmonicsSpectrumDeriv(dd, cr, ci, srx, six, sry, siy, srz, siz, sh.cg, sh.indk, 
                sh.indl, sh.indm, sh.rowm, pairnumsum, common.Nub, common.Ncg, na, common.L, common.K, common.spectrum, backend);                
        
        dstype *fij = &tmp.tmpmem[0]; // dim * ntuples 
        PGEMNV(common.cublasHandle, dim*ntuples, Nbf, &minusone, dd, dim*ntuples, coeff, inc1, &zero, fij, inc1, common.backend);    

// #ifdef HAVE_CHECK    
//         dstype en = cpuArraySum(e, na);
//         dstype epsil = 1e-6;        
//         dstype *sr1 =  new dstype[ntuples*nbasis];
//         dstype *si1 =  new dstype[ntuples*nbasis]; 
//         dstype *cr1 =  new dstype[na*nbasis]; 
//         dstype *ci1 =  new dstype[na*nbasis]; 
//         dstype *ei1 =  new dstype[na*Nbf]; 
//         dstype *onev =  new dstype[na]; 
//         dstype *et =  new dstype[Nbf]; 
//         ArraySetValue(onev, 1.0, na, common.backend);                
//         for (int i = 0; i<dim*ntuples; i++) {
//             cpuArrayCopy(xi1, xii, dim*ntuples);  
//             xi1[i] = xi1[i] + epsil;
//             
//             cpuSphericalHarmonicsBessel(sr1, si1, xi1, sh.x0, sh.P, sh.tmp, sh.f, 
//                      sh.fac, M_PI, common.L, common.K, ntuples);
//                                                          
//             cpuRadialSphericalHarmonicsSpectrum(ei1, cr1, ci1, sr1, si1, sh.cg, sh.indk, sh.indl, sh.indm, 
//                     sh.rowm, pairnumsum, common.Nub, common.Ncg, na, common.L, common.K, common.spectrum);
//                         
//             PGEMTV(common.cublasHandle, na, Nbf, &one, ei1, na, onev, inc1, &zero, et, inc1, common.backend);                    
// 
//             dstype em = 0.0;
//             for (int j=0; j<Nbf; j++)
//                 em += et[j]*coeff[j];                                    
//             dstype emax = fabs((en-em)/epsil - fij[i]);   
//                         
//             cout<<"Maximum absolute error: "<<emax<<endl;      
//         }
//         delete[] sr1;
//         delete[] si1;
//         delete[] cr1;
//         delete[] ci1;
//         delete[] xi1;
//         delete[] xii;
//         delete[] ei1;
//         delete[] onev;
//         delete[] et;        
//         error("here");            
// #endif
                        
        if (decomp==0)
            ForceDecomposition(f, fij, ai, aj, ntuples, backend);        
        else {
            CenterAtomDecomposition(f, fij, ilist, pairnumsum, na, backend);            
            ArrayCopy(tmp.intmem, aj, ntuples, backend);
            Int *jlist = &tmp.intmem[ntuples];   
            Int *bnumsum = &tmp.intmem[2*ntuples]; 
            Int *index = &tmp.intmem[3*ntuples]; // ntuples       
            Int *p0 = &tmp.intmem[4*ntuples]; // ntuples       
            Int *p1 = &tmp.intmem[5*ntuples]; // ntuples       
            Int *p2 = &tmp.intmem[6*ntuples]; // ntuples       
            Int *p3 = &tmp.intmem[7*ntuples]; // ntuples       
            Int naj = UniqueSort(jlist, bnumsum, index, p0, tmp.intmem, p1, p2, p3, ntuples, backend);                        
            NeighborAtomDecomposition(f, fij, jlist, bnumsum, index, naj, backend);            
        }
        SnapTallyVirialFull(v, fij, xij, ai, aj, common.inum, ntuples, backend);  
    }                
}

void implSphericalHarmonicBesselEnergyForceVirial(dstype *e, dstype *f, dstype *v, neighborstruct &nb, commonstruct &common, 
        appstruct &app, tempstruct &tmp, shstruct &sh, dstype* x, dstype *coeff, dstype *q, dstype* param, Int nparam)
{    
    //ArraySetValue(e, 0.0, common.inum, common.backend);
    //ArraySetValue(f, 0.0, common.dim*common.inum, common.backend);
    Int natomtypes = common.natomtypes;
    if (natomtypes==1) {
        SphericalHarmonicBesselEnergyForceVirial(e, f, v, nb, common, app, tmp, sh, x, coeff, q, param, &app.rcutsqml[0], 
            nb.atomtype, nparam, 0, 0, common.decomposition);                    
    }                 
    else {
        Int Nbf = common.Nbf;
        Int inum = common.inum;
        if (common.chemtype == 0) {
            for (int i = 0; i < natomtypes; i++)
                SphericalHarmonicBesselEnergyForceVirial(e, f, v, nb, common, app, tmp, sh, 
                        x, &coeff[i*Nbf], q, param, &app.rcutsqml[0], nb.atomtype, nparam, i+1, 0, common.decomposition);                                    
        }
        else {
            for (int i = 0; i < natomtypes; i++)
                for (int j = 0; j < natomtypes; j++)
                    SphericalHarmonicBesselEnergyForceVirial(e, f, v, nb, common, app, tmp, sh, 
                       x, &coeff[(j+i*natomtypes)*Nbf], q, param, &app.rcutsqml[0], nb.atomtype, nparam, i+1, j+1, common.decomposition);                                    
        }
    }
}

#endif        

