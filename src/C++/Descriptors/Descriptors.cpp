#ifndef __DESCRIPTORS
#define __DESCRIPTORS

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
    
    //printArray2D(x0, 1, K*(1+L), backend);    
    // free memory
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
        CHECK( cudaMemcpy(sh.indk, indk, K2*sizeof(Int), cudaMemcpyHostToDevice ) );      
        CHECK( cudaMemcpy(sh.indl, indl, Nub*3*sizeof(Int), cudaMemcpyHostToDevice ) );      
        CHECK( cudaMemcpy(sh.indm, indm, Ncg*3*sizeof(Int), cudaMemcpyHostToDevice ) );      
        CHECK( cudaMemcpy(sh.rowm, rowm, (Nub+1)*sizeof(Int), cudaMemcpyHostToDevice ) );   
        CHECK( cudaMemcpy(sh.x0, x0, ((L+1)*K)*sizeof(dstype), cudaMemcpyHostToDevice ) );      
        CHECK( cudaMemcpy(sh.cg, cg, Ncg*sizeof(dstype), cudaMemcpyHostToDevice ) );      
        CHECK( cudaMemcpy(sh.fac, fac, 168*sizeof(dstype), cudaMemcpyHostToDevice ) );              
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
        Int jnum = common.jnum;
        Int dim = common.dim;
                
        Int *ilist = &tmp.intmem[0]; //na     
        if (typei>0) {               
            Int *olist = &tmp.intmem[na]; //na        
            cpuArrayFill(olist, e1, na);        

            Int *t0 = &tmp.intmem[2*na]; //na        
            Int *t1 = &tmp.intmem[3*na]; //na        
            na = cpuFindAtomType(ilist, olist, atomtype, t0, t1, typei, na);
        }
        else {
            cpuArrayFill(ilist, e1, na);        
        }
        
        Int *pairnum = &tmp.intmem[na]; // na
        Int *pairlist = &tmp.intmem[2*na]; // na*jnum
        if (typej>0)
            cpuFullNeighPairList(pairnum, pairlist, x, rcutsq, atomtype, ilist, nb.alist, nb.neighlist, nb.neighnum, na, jnum, typej, dim);
        else
            cpuFullNeighPairList(pairnum, pairlist, x, rcutsq, ilist, nb.neighlist, nb.neighnum, na, jnum, dim);        
                
        //a list contains the starting positions of the first neighbor 
        Int *pairnumsum = &tmp.intmem[2*na+na*jnum]; // na+1        
        cpuCumsum(pairnumsum, pairnum, na+1);                                         
        int ntuples = IntArrayGetValueAtIndex(pairnumsum, na, common.backend);     
                
        Int *ai = &tmp.intmem[1+3*na+na*jnum]; // ntuples        
        Int *aj = &tmp.intmem[1+3*na+ntuples+na*jnum]; // ntuples        
        Int *ti = &tmp.intmem[1+3*na+2*ntuples+na*jnum]; // ntuples        
        Int *tj = &tmp.intmem[1+3*na+3*ntuples+na*jnum]; // ntuples        
        //dstype *xij = &tmp.tmpmem[0]; // ntuples*dim
        Int nbasis = common.K*(common.L+1)*(common.L+1);
        dstype *xij = &tmp.tmpmem[2*na*nbasis+ntuples*(2*nbasis)]; // ntuples*dim
        dstype *qi;
        dstype *qj;
        cpuNeighPairs(xij, qi, qj, x, q, ai, aj, ti, tj, pairnum, pairlist, pairnumsum, ilist, nb.alist, 
                atomtype, na, jnum, ncq, dim);       
                                                
        dstype *cr =  &tmp.tmpmem[0];         // na*nbasis
        dstype *ci =  &tmp.tmpmem[na*nbasis]; // na*nbasis              
        dstype *sr =  &tmp.tmpmem[2*na*nbasis]; // ntuples*nbasis
        dstype *si =  &tmp.tmpmem[2*na*nbasis+ntuples*(1*nbasis)]; // ntuples*nbasis                       
        cpuSphericalHarmonicsBessel(sr, si, xij, sh.x0, sh.P, sh.tmp, sh.f, 
                sh.fac, M_PI, common.L, common.K, ntuples);
                                
        dstype *ei = &tmp.tmpmem[2*na*nbasis+ntuples*(2*nbasis)]; // na * Nbf
        cpuRadialSphericalHarmonicsSpectrum(ei, cr, ci, sr, si, sh.cg, sh.indk, sh.indl, sh.indm, 
                sh.rowm, pairnumsum, common.Nub, common.Ncg, na, common.L, common.K, common.spectrum);
        
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
    for (Int b=0; b<common.nba; b++) {
        Int e1 = common.ablks[b];
        Int e2 = common.ablks[b+1];            
        Int na = e2 - e1; // number of atoms in this block
        Int jnum = common.jnum;
        Int dim = common.dim;
                
        Int *ilist = &tmp.intmem[0]; //na     
        if (typei>0) {               
            Int *olist = &tmp.intmem[na]; //na        
            cpuArrayFill(olist, e1, na);        

            Int *t0 = &tmp.intmem[2*na]; //na        
            Int *t1 = &tmp.intmem[3*na]; //na        
            na = cpuFindAtomType(ilist, olist, atomtype, t0, t1, typei, na);
        }
        else {
            cpuArrayFill(ilist, e1, na);        
        }
        
        Int *pairnum = &tmp.intmem[na]; // na
        Int *pairlist = &tmp.intmem[2*na]; // na*jnum
        if (typej>0)
            cpuFullNeighPairList(pairnum, pairlist, x, rcutsq, atomtype, ilist, nb.alist, nb.neighlist, nb.neighnum, na, jnum, typej, dim);
        else
            cpuFullNeighPairList(pairnum, pairlist, x, rcutsq, ilist, nb.neighlist, nb.neighnum, na, jnum, dim);        
                
        //a list contains the starting positions of the first neighbor 
        Int *pairnumsum = &tmp.intmem[2*na+na*jnum]; // na+1        
        cpuCumsum(pairnumsum, pairnum, na+1);                                         
        int ntuples = IntArrayGetValueAtIndex(pairnumsum, na, common.backend);     
                
        Int *ai = &tmp.intmem[1+3*na+na*jnum]; // ntuples        
        Int *aj = &tmp.intmem[1+3*na+ntuples+na*jnum]; // ntuples        
        Int *ti = &tmp.intmem[1+3*na+2*ntuples+na*jnum]; // ntuples        
        Int *tj = &tmp.intmem[1+3*na+3*ntuples+na*jnum]; // ntuples        
        //dstype *xij = &tmp.tmpmem[0]; // ntuples*dim
        Int nbasis = common.K*(common.L+1)*(common.L+1);
        dstype *xij = &tmp.tmpmem[2*na*nbasis+ntuples*(8*nbasis)]; // ntuples*dim
        dstype *qi;
        dstype *qj;
        cpuNeighPairs(xij, qi, qj, x, q, ai, aj, ti, tj, pairnum, pairlist, pairnumsum, ilist, nb.alist, 
                atomtype, na, jnum, ncq, dim);       
                                                
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
        cpuSphericalHarmonicsBesselWithDeriv(sr, si, srx, six, sry, siy, srz, siz, xij,
                 sh.x0, sh.P, sh.tmp, sh.f, sh.dP, sh.dtmp, sh.df, sh.fac, M_PI, common.L, common.K, ntuples);
                
#ifdef CHECK  // check              
        dstype epsil = 1e-6;
        dstype *sr1 =  new dstype[ntuples*nbasis];
        dstype *si1 =  new dstype[ntuples*nbasis]; 
        dstype *xi1 =  new dstype[ntuples*dim]; 
        dstype *xii =  new dstype[ntuples*dim]; 
        dstype *ei1 =  new dstype[na*Nbf]; 
        dstype *onev =  new dstype[na]; 
        dstype *et =  new dstype[Nbf]; 
        cpuArrayCopy(xii, xij, dim*ntuples);
        cpuArrayCopy(xi1, xij, dim*ntuples);
        ArraySetValue(onev, 1.0, na, common.backend);
        
        cpuArrayRowkAXPB(xi1, xij, 1.0, epsil, dim, ntuples, 0);
        cpuSphericalHarmonicsBessel(sr1, si1, xi1, sh.x0, sh.P, sh.tmp, sh.f, 
                 sh.fac, M_PI, common.L, common.K, ntuples);
        cpuArrayAXPBY(sr1, sr1, sr, 1.0/epsil, -1.0/epsil, ntuples*nbasis);
        cpuArrayAXPBY(sr1, sr1, srx, 1.0, -1.0, ntuples*nbasis);
        cpuArrayAbs(sr1, sr1, ntuples*nbasis);        
        cout<<"Maximum absolute error: "<<cpuArrayMax(sr1, ntuples*nbasis)<<endl;        
        
        cpuArrayCopy(xi1, xij, dim*ntuples);
        cpuArrayRowkAXPB(xi1, xij, 1.0, epsil, dim, ntuples, 1);
        cpuSphericalHarmonicsBessel(sr1, si1, xi1, sh.x0, sh.P, sh.tmp, sh.f, 
                 sh.fac, M_PI, common.L, common.K, ntuples);
        cpuArrayAXPBY(sr1, sr1, sr, 1.0/epsil, -1.0/epsil, ntuples*nbasis);
        cpuArrayAXPBY(sr1, sr1, sry, 1.0, -1.0, ntuples*nbasis);
        cpuArrayAbs(sr1, sr1, ntuples*nbasis);        
        cout<<"Maximum absolute error: "<<cpuArrayMax(sr1, ntuples*nbasis)<<endl;         
        
        cpuArrayCopy(xi1, xij, dim*ntuples);
        cpuArrayRowkAXPB(xi1, xij, 1.0, epsil, dim, ntuples, 2);
        cpuSphericalHarmonicsBessel(sr1, si1, xi1, sh.x0, sh.P, sh.tmp, sh.f, 
                 sh.fac, M_PI, common.L, common.K, ntuples);
        cpuArrayAXPBY(sr1, sr1, sr, 1.0/epsil, -1.0/epsil, ntuples*nbasis);
        cpuArrayAXPBY(sr1, sr1, srz, 1.0, -1.0, ntuples*nbasis);
        cpuArrayAbs(sr1, sr1, ntuples*nbasis);        
        cout<<"Maximum absolute error: "<<cpuArrayMax(sr1, ntuples*nbasis)<<endl;    
#endif        
        
        dstype *ei = &tmp.tmpmem[2*na*nbasis+ntuples*(8*nbasis)]; // na * Nbf
        cpuRadialSphericalHarmonicsSpectrum(ei, cr, ci, sr, si, sh.cg, sh.indk, sh.indl, sh.indm, 
                sh.rowm, pairnumsum, common.Nub, common.Ncg, na, common.L, common.K, common.spectrum);
        
        // ei : na x Nbf,  onevec : na x 1, e : Nbf x 1
        // e = e + ei^T*onevec        
        dstype *onevec =  &tmp.tmpmem[2*na*nbasis+ntuples*(6*nbasis)];  // na
        ArraySetValue(onevec, 1.0, na, common.backend);
        PGEMTV(common.cublasHandle, na, Nbf, &one, ei, na, onevec, inc1, &one, e, inc1, common.backend);        
                
        dstype *fij = &tmp.tmpmem[2*na*nbasis+ntuples*(6*nbasis)]; // dim * ntuples * Nbf
        cpuRadialSphericalHarmonicsSpectrumDeriv(fij, cr, ci, srx, six, sry, siy, srz, siz, sh.cg, sh.indk, 
                sh.indl, sh.indm, sh.rowm, pairnumsum, common.Nub, common.Ncg, na, common.L, common.K, common.spectrum);                        
        
#ifdef CHECK    // check        
        for (int i = 0; i<dim*ntuples; i++) {
            cpuArrayCopy(xi1, xii, dim*ntuples);  
            xi1[i] = xi1[i] + epsil;
            cpuSphericalHarmonicsBessel(sr1, si1, xi1, sh.x0, sh.P, sh.tmp, sh.f, 
                     sh.fac, M_PI, common.L, common.K, ntuples);
                                 
            cpuRadialSphericalHarmonicsSpectrum(ei1, cr, ci, sr1, si1, sh.cg, sh.indk, sh.indl, sh.indm, 
                    sh.rowm, pairnumsum, common.Nub, common.Ncg, na, common.L, common.K, common.spectrum);
            
            PGEMTV(common.cublasHandle, na, Nbf, &one, ei1, na, onev, inc1, &zero, et, inc1, common.backend);                    
            printArray2D(e, 1, Nbf, common.backend);
            printArray2D(et, 1, Nbf, common.backend);
            for (int j=0; j<Nbf; j++) 
                et[j] = fabs((et[j]-e[j])/epsil - fij[i + dim*ntuples*j]);
            cout<<"Maximum absolute error: "<<cpuArrayMax(et, Nbf)<<endl;      
        }
        delete[] sr1;
        delete[] si1;
        delete[] xi1;
        delete[] xii;
        delete[] ei1;
        delete[] onev;
        delete[] et;        
        error("here");            
#endif
                
        if (decomp==0)
            cpuForceDecomposition(f, fij, ai, aj, common.inum, ntuples, Nbf);
        else {
            cpuCenterAtomDecomposition(f, fij, ilist, pairnumsum, common.inum, ntuples, na, Nbf);
            cpuArrayCopy(tmp.intmem, aj, ntuples);
            Int *jlist = &tmp.intmem[ntuples];   
            Int *bnumsum = &tmp.intmem[2*ntuples]; 
            Int *index = &tmp.intmem[3*ntuples]; // ntuples       
            Int *p0 = &tmp.intmem[4*ntuples]; // ntuples       
            Int *p1 = &tmp.intmem[5*ntuples]; // ntuples       
            Int *p2 = &tmp.intmem[6*ntuples]; // ntuples       
            Int *p3 = &tmp.intmem[7*ntuples]; // ntuples       
            Int naj = cpuUniqueSort(jlist, bnumsum, index, p0, tmp.intmem, p1, p2, p3, ntuples);                        
            cpuNeighborAtomDecomposition(f, fij, jlist, bnumsum, index, common.inum, ntuples, naj, Nbf);
        }        
    }        
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
        Int jnum = common.jnum;
        Int dim = common.dim;
                
        Int *ilist = &tmp.intmem[0]; //na     
        if (typei>0) {               
            Int *olist = &tmp.intmem[na]; //na        
            cpuArrayFill(olist, e1, na);        

            Int *t0 = &tmp.intmem[2*na]; //na        
            Int *t1 = &tmp.intmem[3*na]; //na        
            na = cpuFindAtomType(ilist, olist, atomtype, t0, t1, typei, na);
        }
        else {
            cpuArrayFill(ilist, e1, na);        
        }
        
        Int *pairnum = &tmp.intmem[na]; // na
        Int *pairlist = &tmp.intmem[2*na]; // na*jnum
        if (typej>0)
            cpuFullNeighPairList(pairnum, pairlist, x, rcutsq, atomtype, ilist, nb.alist, nb.neighlist, nb.neighnum, na, jnum, typej, dim);
        else
            cpuFullNeighPairList(pairnum, pairlist, x, rcutsq, ilist, nb.neighlist, nb.neighnum, na, jnum, dim);        
                        
        //a list contains the starting positions of the first neighbor 
        Int *pairnumsum = &tmp.intmem[2*na+na*jnum]; // na+1        
        cpuCumsum(pairnumsum, pairnum, na+1);                                         
        int ntuples = IntArrayGetValueAtIndex(pairnumsum, na, common.backend);     
                
        Int *ai = &tmp.intmem[1+3*na+na*jnum]; // ntuples        
        Int *aj = &tmp.intmem[1+3*na+ntuples+na*jnum]; // ntuples        
        Int *ti = &tmp.intmem[1+3*na+2*ntuples+na*jnum]; // ntuples        
        Int *tj = &tmp.intmem[1+3*na+3*ntuples+na*jnum]; // ntuples        
        //dstype *xij = &tmp.tmpmem[0]; // ntuples*dim
        Int nbasis = common.K*(common.L+1)*(common.L+1);
        dstype *xij = &tmp.tmpmem[2*na*nbasis+ntuples*(8*nbasis)]; // ntuples*dim
        dstype *qi;
        dstype *qj;
        cpuNeighPairs(xij, qi, qj, x, q, ai, aj, ti, tj, pairnum, pairlist, pairnumsum, ilist, nb.alist, 
                atomtype, na, jnum, ncq, dim);       
                                
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
        cpuSphericalHarmonicsBesselWithDeriv(sr, si, srx, six, sry, siy, srz, siz, xij,
                 sh.x0, sh.P, sh.tmp, sh.f, sh.dP, sh.dtmp, sh.df, sh.fac, M_PI, common.L, common.K, ntuples);
                
        dstype *di = &tmp.tmpmem[2*na*nbasis+ntuples*(8*nbasis)]; // na * Nbf                
        cpuRadialSphericalHarmonicsSpectrum(di, cr, ci, sr, si, sh.cg, sh.indk, 
                sh.indl, sh.indm, sh.rowm, pairnumsum, common.Nub, common.Ncg, na, common.L, common.K, common.spectrum);

        // di : na x Nbf,  coeff : Nbf x 1, ei : na x 1
        // ei = di*coeff                        
        dstype *ei = &tmp.tmpmem[2*na*nbasis+ntuples*(6*nbasis)]; // na         
        PGEMNV(common.cublasHandle, na, Nbf, &one, di, na, coeff, inc1, &zero, ei, inc1, common.backend);                
        cpuCenterAtomDecomposition(e, ei, ilist, na);
        
        dstype *dd = &tmp.tmpmem[2*na*nbasis+ntuples*(6*nbasis)]; // dim * ntuples * Nbf
        cpuRadialSphericalHarmonicsSpectrumDeriv(dd, cr, ci, srx, six, sry, siy, srz, siz, sh.cg, sh.indk, 
                sh.indl, sh.indm, sh.rowm, pairnumsum, common.Nub, common.Ncg, na, common.L, common.K, common.spectrum);                
        
        dstype *fij = &tmp.tmpmem[0]; // dim * ntuples 
        PGEMNV(common.cublasHandle, dim*ntuples, Nbf, &minusone, dd, dim*ntuples, coeff, inc1, &zero, fij, inc1, common.backend);    
        
        if (decomp==0)
            cpuForceDecomposition(f, fij, ai, aj, ntuples);        
        else {
            cpuCenterAtomDecomposition(f, fij, ilist, pairnumsum, na);            
            cpuArrayCopy(tmp.intmem, aj, ntuples);
            Int *jlist = &tmp.intmem[ntuples];   
            Int *bnumsum = &tmp.intmem[2*ntuples]; 
            Int *index = &tmp.intmem[3*ntuples]; // ntuples       
            Int *p0 = &tmp.intmem[4*ntuples]; // ntuples       
            Int *p1 = &tmp.intmem[5*ntuples]; // ntuples       
            Int *p2 = &tmp.intmem[6*ntuples]; // ntuples       
            Int *p3 = &tmp.intmem[7*ntuples]; // ntuples       
            Int naj = cpuUniqueSort(jlist, bnumsum, index, p0, tmp.intmem, p1, p2, p3, ntuples);                        
            cpuNeighborAtomDecomposition(f, fij, jlist, bnumsum, index, naj);            
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

