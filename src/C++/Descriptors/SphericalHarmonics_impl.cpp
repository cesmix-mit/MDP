#ifndef __SPHERICALHARMONICS_IMPL
#define __SPHERICALHARMONICS_IMPL

void InitSphericalHarmonic(shstruct &sh, Int Kin, Int Lin, Int backendin) 
{
//     L       the maximum degree of spherical harmonics 
//     K       the number of zeros of spherical Bessel functions
//     Nub     number of non-zero unqiue bispectrum compoments of spherical harmonics 
//     Ncg     the total number of non-zero Clebsch-Gordan coefficients 
//     backend computing platform
    
    Int L = Lin;
    Int K = Kin;
    Int backend = backendin;
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
//     dstype the[N];
//     dstype phi[N];    
//     dstype r[N];    
//     std::random_device rd;
//     std::mt19937 gen(rd());
//     std::uniform_real_distribution<> dis(0, 1);//uniform distribution between 0 and 1
//     for (Int n = 0; n < N; ++n) {
//         the[n] = M_PI*dis(gen);
//         phi[n] = 2*M_PI*dis(gen) - M_PI;
//         r[n] = 0.1 + dis(gen);
//     }
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
    cpuSphericalHarmonicsSum(Ylmr, Ylmi, x, y, z, P, tmp, fac, M_PI, L, N);
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
            x0[k*(L+1) + l] = bzeros[k*25 + l];
    
    // free memory
    sh.freememory(backend);
         
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
        
//     Int ia[10] = {Nub, Ncg, L, K, M, L2, K2, N};
//     print1iarray(ia, 10);
    
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
    
    sh.L = L;    // the maximum degree of spherical harmonics 
    sh.K = K;    // the number of zeros of spherical Bessel functions
    sh.Nub = Nub;// number of non-zero unqiue bispectrum compoments of spherical harmonics 
    sh.Ncg = Ncg;// the total number of non-zero Clebsch-Gordan coefficients     
}


void implSphericalHarmonicsBessel(dstype *Sr, dstype *Si, dstype *xij, Int N, Int backend, shstruct &sh)
{
    coreSphericalHarmonicsBessel(Sr, Si, xij, sh.x0, sh.P, sh.tmp, sh.f, sh.fac, M_PI, sh.L, sh.K, N, backend);
}

void implSphericalHarmonicsBesselWithDeriv(dstype *Sr, dstype *Si, dstype *Srx, dstype *Six, dstype *Sry, dstype *Siy, dstype *Srz, dstype *Siz, dstype *xij, Int N, Int backend, shstruct &sh)
{
    coreSphericalHarmonicsBesselWithDeriv(Sr, Si, Srx, Six, Sry, Siy, Srz, Siz, xij, sh.x0, sh.P, sh.tmp, sh.f, sh.dP, sh.dtmp, sh.df, sh.fac, M_PI, sh.L, sh.K, N, backend);    
}

void implSphericalHarmonicsBessel(dstype *Sr, dstype *Si, dstype *x, dstype *y, dstype *z, Int N, Int backend, shstruct &sh)
{
    coreSphericalHarmonicsBessel(Sr, Si, x, y, z, sh.x0, sh.P, sh.tmp, sh.f, sh.fac, M_PI, sh.L, sh.K, N, backend);
}

void implSphericalHarmonicsBesselDeriv(dstype *Srx, dstype *Six, dstype *Sry, dstype *Siy, dstype *Srz, dstype *Siz, dstype *x, dstype *y, dstype *z, Int N, Int backend, shstruct &sh)
{
    coreSphericalHarmonicsBesselDeriv(Srx, Six, Sry, Siy, Srz, Siz, x, y, z, sh.x0, sh.P, sh.tmp, sh.f, sh.dP, sh.dtmp, sh.df, sh.fac, M_PI, sh.L, sh.K, N, backend);    
}

void implSphericalHarmonicsBesselWithDeriv(dstype *Sr, dstype *Si, dstype *Srx, dstype *Six, dstype *Sry, dstype *Siy, dstype *Srz, dstype *Siz, dstype *x, dstype *y, dstype *z, Int N, Int backend, shstruct &sh)
{
    coreSphericalHarmonicsBesselWithDeriv(Sr, Si, Srx, Six, Sry, Siy, Srz, Siz, x, y, z, sh.x0, sh.P, sh.tmp, sh.f, sh.dP, sh.dtmp, sh.df, sh.fac, M_PI, sh.L, sh.K, N, backend);    
}

void implRadialSphericalHarmonicsSum(dstype *ar, dstype *ai, dstype *Sr, dstype *Si, Int *numneigh, Int Na, Int backend, shstruct &sh)
{
    coreRadialSphericalHarmonicsSum(ar, ai, Sr, Si, numneigh, Na, sh.L, sh.K, backend);
}

void implRadialSphericalHarmonicsPower(dstype *p, dstype *ar, dstype *ai, Int Na, Int backend, shstruct &sh)
{
    coreRadialSphericalHarmonicsPower(p, ar, ai, sh.indk, Na, sh.L, sh.K, backend);
}

void implRadialSphericalHarmonicsPowerDeriv(dstype *px, dstype *py, dstype *pz, dstype *ar, dstype *ai, 
        dstype *arx, dstype *aix, dstype *ary, dstype *aiy, dstype *arz, dstype *aiz, Int *numneigh, Int Na, Int backend, shstruct &sh)
{
    coreRadialSphericalHarmonicsPowerDeriv(px, py, pz, ar, ai, arx, aix, ary, aiy, arz, aiz, sh.indk, numneigh, Na, sh.L, sh.K, backend);
}

void implRadialSphericalHarmonicsPowerDeriv2(dstype *pd, dstype *ar, dstype *ai, 
        dstype *arx, dstype *aix, dstype *ary, dstype *aiy, dstype *arz, dstype *aiz, Int *numneigh, Int Na, Int backend, shstruct &sh)
{
    coreRadialSphericalHarmonicsPowerDeriv2(pd, ar, ai, arx, aix, ary, aiy, arz, aiz, sh.indk, numneigh, Na, sh.L, sh.K, backend);
}

void implRadialSphericalHarmonicsBispectrum(dstype *b, dstype *ar, dstype *ai, Int Na, Int backend, shstruct &sh)
{    
    coreRadialSphericalHarmonicsBispectrum(b, ar, ai, sh.cg, sh.indk, sh.indl, sh.indm, sh.rowm, sh.Nub, sh.Ncg, Na, sh.L, sh.K, backend);
}

void implRadialSphericalHarmonicsBispectrumDeriv(dstype *bx, dstype *by, dstype *bz, dstype *ar, dstype *ai, 
        dstype *arx, dstype *aix, dstype *ary, dstype *aiy, dstype *arz, dstype *aiz, Int *numneigh, Int Na, Int backend, shstruct &sh)
{
    coreRadialSphericalHarmonicsBispectrumDeriv(bx, by, bz, ar, ai, arx, aix, ary, aiy, arz, aiz, 
             sh.cg, sh.indk, sh.indl, sh.indm, sh.rowm, numneigh, Na, sh.Nub, sh.Ncg, sh.K, backend);
}

void implRadialSphericalHarmonicsBispectrumDeriv2(dstype *bd, dstype *ar, dstype *ai, 
        dstype *arx, dstype *aix, dstype *ary, dstype *aiy, dstype *arz, dstype *aiz, Int *numneigh, Int Na, Int backend, shstruct &sh)
{
    coreRadialSphericalHarmonicsBispectrumDeriv2(bd, ar, ai, arx, aix, ary, aiy, arz, aiz, 
             sh.cg, sh.indk, sh.indl, sh.indm, sh.rowm, numneigh, Na, sh.Nub, sh.Ncg, sh.K, backend);
}

void implRadialSphericalHarmonicsBasis(dstype *d, dstype *c, Int *atomtype, Int Ntype, Int Na, Int Nbf, Int backend)
{
    coreRadialSphericalHarmonicsBasis(d, c, atomtype, Ntype, Na, Nbf, backend);
}

void implRadialSphericalHarmonicsBasisDeriv(dstype *dx, dstype *dy, dstype *dz, dstype *cx, dstype *cy, dstype *cz,
        Int *atomtype, Int *neighlist, Int *numneigh, Int Ntype, Int Na, Int Nbf, Int backend)
{
    coreRadialSphericalHarmonicsBasisDeriv(dx, dy, dz, cx, cy, cx, atomtype, neighlist, numneigh, Ntype, Na, Nbf, backend);
}

void implRadialSphericalHarmonicsBasisDeriv2(dstype *dd, dstype *cd,
        Int *atomtype, Int *neighlist, Int *numneigh, Int Ntype, Int Na, Int Nbf, Int backend)
{
    coreRadialSphericalHarmonicsBasisDeriv2(dd, cd, atomtype, neighlist, numneigh, Ntype, Na, Nbf, backend);
}

#endif        

