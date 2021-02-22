#ifndef __CONFIGURATION_IMPL
#define __CONFIGURATION_IMPL

void implReadAppStruct(appstruct &app, string filein, Int mpiprocs, Int mpirank, Int backend)
{
    string filename = filein + "app.bin";                    
    
    // Open file to read
    ifstream in(filename.c_str(), ios::in | ios::binary);
    
    if (!in) 
        error("Unable to open file " + filename);
       
    if (mpirank==0)
        printf("Read app struct from a binary file...\n");   
    
    /* Read data to app structure */            
    app.lsize = readiarrayfromdouble(in, 1);
    app.nsize = readiarrayfromdouble(in, app.lsize[0]);
    app.ndims = readiarrayfromdouble(in, app.nsize[0]);
    app.flags = readiarrayfromdouble(in, app.nsize[1]);
    app.bcs = readiarrayfromdouble(in, app.nsize[2]);
    app.pbc = readiarrayfromdouble(in, app.nsize[3]);
    app.atomnumber = readiarrayfromdouble(in, app.nsize[4]);    
    
    readarray(in, &app.atommass, app.nsize[5]);    
    readarray(in, &app.atomcharge, app.nsize[6]);    
    readarray(in, &app.physicsparam, app.nsize[7]);       
    readarray(in, &app.solversparam, app.nsize[8]);   
    readarray(in, &app.simulaparam, app.nsize[9]);   
    readarray(in, &app.singleparam, app.nsize[10]);       
    readarray(in, &app.pairparam, app.nsize[11]);   
    readarray(in, &app.tripletparam, app.nsize[12]);   
    readarray(in, &app.coeff, app.nsize[13]);   
    readarray(in, &app.cutoffellipsoid, app.nsize[14]);   
    readarray(in, &app.atomellipsoid, app.nsize[15]);   
    readarray(in, &app.rcutsq, app.nsize[16]);  
    readarray(in, &app.boxoffset, app.nsize[17]);  

//     dim = app.ndims[0];
//     TemplateMalloc(&app.a, dim, 0); 
//     TemplateMalloc(&app.b, dim, 0); 
//     TemplateMalloc(&app.c, dim, 0);    
//     TemplateMalloc(&app.s2rmap, dim*dim, 0);      
//     
//     Int n = (dim==2) ? dim*4 : dim*8;         
//     TemplateMalloc(&app.refvertices, n, 0); 
//     TemplateMalloc(&app.rbvertices,  n, 0); 
//     TemplateMalloc(&app.boxvertices, n, 0); 
//     TemplateMalloc(&app.bbvertices, n, 0);    
//     
//     n = (dim==2) ? dim*9 : dim*27;         
//     TemplateMalloc(&app.pimages, n, 0);    
//     
//     TemplateMalloc(&app.cellnum, dim, 0);     
//     for (Int i=0; i<dim; i++) {
//         app.cellnum[i] = (Int) floor(1.0/app.cellsize[i]);   
//         if (app.pbc[i]==1)
//             app.cellnum[i] += 2;        
//     }
    
    // Close file:
    in.close();        
}

void implSetAppStruct(appstruct &app, appstruct &happ, Int backend)
{    
    /* Allocate memory on GPU */            
    TemplateMalloc(&app.nsize, happ.lsize[0], backend); 
    TemplateMalloc(&app.ndims, happ.nsize[0], backend); 
    TemplateMalloc(&app.flags, happ.nsize[1], backend); 
    TemplateMalloc(&app.bcs, happ.nsize[2], backend); 
    TemplateMalloc(&app.pbc, happ.nsize[3], backend); 
    TemplateMalloc(&app.atomnumber, happ.nsize[4], backend); 
    TemplateMalloc(&app.atommass, happ.nsize[5], backend); 
    TemplateMalloc(&app.atomcharge, happ.nsize[6], backend); 
    TemplateMalloc(&app.physicsparam, happ.nsize[7], backend); 
    TemplateMalloc(&app.solversparam, happ.nsize[8], backend); 
    TemplateMalloc(&app.simulaparam, happ.nsize[9], backend); 
    TemplateMalloc(&app.singleparam, happ.nsize[10], backend); 
    TemplateMalloc(&app.pairparam, happ.nsize[11], backend); 
    TemplateMalloc(&app.tripletparam, happ.nsize[12], backend); 
    TemplateMalloc(&app.coeff, happ.nsize[13], backend); 
    TemplateMalloc(&app.cutoffellipsoid, happ.nsize[14], backend); 
    TemplateMalloc(&app.atomellipsoid, happ.nsize[15], backend); 
    TemplateMalloc(&app.rcutsq, happ.nsize[16], backend);     
    TemplateMalloc(&app.boxoffset, happ.nsize[17], backend);     
    
    if (backend==2) { // GPU
#ifdef HAVE_CUDA        
        CHECK( cudaMemcpy(app.nsize, happ.nsize, happ.lsize[0]*sizeof(Int), cudaMemcpyHostToDevice ) );              
        CHECK( cudaMemcpy(app.ndims, happ.ndims, happ.nsize[0]*sizeof(Int), cudaMemcpyHostToDevice ) );              
        CHECK( cudaMemcpy(app.flags, happ.flags, happ.nsize[1]*sizeof(Int), cudaMemcpyHostToDevice ) );              
        CHECK( cudaMemcpy(app.bcs, happ.bcs, happ.nsize[2]*sizeof(Int), cudaMemcpyHostToDevice ) );              
        CHECK( cudaMemcpy(app.pbc, happ.pbc, happ.nsize[3]*sizeof(Int), cudaMemcpyHostToDevice ) );              
        CHECK( cudaMemcpy(app.atomnumber, happ.atomnumber, happ.nsize[4]*sizeof(Int), cudaMemcpyHostToDevice ) );              
           
        CHECK( cudaMemcpy(app.atommass, happ.atommass, happ.nsize[5]*sizeof(dstype), cudaMemcpyHostToDevice ) );              
        CHECK( cudaMemcpy(app.atomcharge, happ.atomcharge, happ.nsize[6]*sizeof(dstype), cudaMemcpyHostToDevice ) );              
        CHECK( cudaMemcpy(app.physicsparam, happ.physicsparam, happ.nsize[7]*sizeof(dstype), cudaMemcpyHostToDevice ) );              
        CHECK( cudaMemcpy(app.silversparam, happ.solversparam, happ.nsize[8]*sizeof(dstype), cudaMemcpyHostToDevice ) );              
        CHECK( cudaMemcpy(app.simulaparam, happ.simulaparam, happ.nsize[9]*sizeof(dstype), cudaMemcpyHostToDevice ) );              
        CHECK( cudaMemcpy(app.singleparam, happ.singleparam, happ.nsize[10]*sizeof(dstype), cudaMemcpyHostToDevice ) );              
        CHECK( cudaMemcpy(app.pairparam, happ.pairparam, happ.nsize[11]*sizeof(dstype), cudaMemcpyHostToDevice ) );              
        CHECK( cudaMemcpy(app.tripletparam, happ.tripletparam, happ.nsize[12]*sizeof(dstype), cudaMemcpyHostToDevice ) );              
        CHECK( cudaMemcpy(app.coeff, happ.coeff, happ.nsize[13]*sizeof(dstype), cudaMemcpyHostToDevice ) );              
        CHECK( cudaMemcpy(app.cutoffellipsoid, happ.cutoffellipsoid, happ.nsize[14]*sizeof(dstype), cudaMemcpyHostToDevice ) );              
        CHECK( cudaMemcpy(app.atomellipsoid, happ.atomellipsoid, happ.nsize[15]*sizeof(dstype), cudaMemcpyHostToDevice ) );              
        CHECK( cudaMemcpy(app.rcutsq, happ.rcutsq, happ.nsize[16]*sizeof(dstype), cudaMemcpyHostToDevice ) );                              
        CHECK( cudaMemcpy(app.boxoffset, happ.boxoffset, happ.nsize[17]*sizeof(dstype), cudaMemcpyHostToDevice ) );                              
#endif        
    }        
}

void implReadConfigStruct(configstruct &config, string filein, Int mpiprocs, Int mpirank, Int backend)
{
    Int filenumber = mpirank+1; //file number     
    string filename = filein + "config" + NumberToString(filenumber) + ".bin";                    
    
    // Open file to read
    ifstream in(filename.c_str(), ios::in | ios::binary);
    
    if (!in) 
        error("Unable to open file " + filename);
       
    if (mpirank==0)
        printf("Read app struct from a binary file...\n");   
    
    /* Read data to config structure */            
    config.lsize = readiarrayfromdouble(in, 1);
    config.nsize = readiarrayfromdouble(in, config.lsize[0]);
    config.ndims = readiarrayfromdouble(in, config.nsize[0]);
    config.natoms = readiarrayfromdouble(in, config.nsize[1]);
    readarray(in, &config.simbox, config.nsize[2]);   
    
    config.t = readiarrayfromdouble(in, config.nsize[3]);
    readarray(in, &config.x, config.nsize[4]);   
    readarray(in, &config.v, config.nsize[5]);   
    readarray(in, &config.f, config.nsize[6]);   
    readarray(in, &config.q, config.nsize[7]);       
            
    TemplateMalloc(&config.natomssum, config.nsize[1]+1, 0); 
    cpuCumsum(config.natomssum, config.natoms, config.nsize[1]+1); 
    
    // Close file:
    in.close();            
}

void implSetCommonStruct(commonstruct &common, appstruct &app, configstruct &config, string filein, string fileout, Int mpiprocs, Int mpirank, Int backend)
{
    common.filein = filein;
    common.fileout = fileout;
    common.mpiProcs = mpiprocs;
    common.mpiRank = mpirank;
    common.backend = backend;
    
    common.dim = app.ndims[0];
    common.K = app.ndims[1];  // order of radial basis functions
    common.L = app.ndims[2];  // order of spherical harmonics          
    common.ntimesteps = app.ndims[3]; // number of time steps
    common.nab = app.ndims[4]; // number of atoms per block

    common.descriptor = app.flags[0];   // descriptor flag: 0 -> Spherical Harmonics Bessel
    common.spectrum = app.flags[1];     // spectrum flag: 0-> power spectrum, 1-> bispectrum, 2-> both power and bispectrum 
    common.training = app.flags[2];     // 0 -> no training, 1 -> Linear regression, 2 -> Gaussian process, 3 -> Neural net
    common.runMD = app.flags[3];        // 0 no MD simulation, 1 -> run MD simulation
    common.potential = app.flags[4];    // 0 -> empirical potential, 1 -> ML potential, 2 -> combined potential    
    common.cutofftype = app.flags[5]; // 0 -> single cut-off raidus for all atoms, 1 -> multiple cut-off radii for atom pairs  
    common.neightype = app.flags[6];  // 0 -> full neighbor list, 1 -> half neighbor list
    common.pairtype = app.flags[7];   // 0 -> neighbor pairs, 1 -> neighbor pairs based on type of atom i, 
                                    // 2 -> neighbor pairs based on type of atom i and type of atom j
    common.energycal = app.flags[8];    // turns energy calculation on or off
    common.forcecal = app.flags[9];     // turns force calculation on or off
    common.stresscal = app.flags[10];   // turns stress calculation on or off
    
    common.nflags = app.nsize[1];
    common.nbcs = app.nsize[2];
    common.npbc = app.nsize[3];    
    common.natomtypes = app.nsize[4]; // number of atomtypes    
    common.nphysicsparam = app.nsize[7];       
    common.nsolversparam = app.nsize[8];   
    common.nsimulaparam = app.nsize[9];   
    common.nsingleparam = app.nsize[10];       
    common.npairparam = app.nsize[11];   
    common.ntripletparam = app.nsize[12];   
    common.ncoeff = app.nsize[13];   
    common.ncutoffellipsoid = app.nsize[14];   
    common.natomellipsoid = app.nsize[15];   
    common.nrcutsq = app.nsize[16];           
            
    common.time = app.simulaparam[0];
    common.dt = app.simulaparam[1];
           
    common.nconfigs = config.nsize[1]; // number of configurations
    common.ns = config.nsize[2];
    common.nt = config.nsize[3];
    common.nx = config.nsize[4];
    common.nv = config.nsize[5];
    common.nf = config.nsize[6];
    common.nq = config.nsize[7];
    common.ncq = common.nq/(common.nt);        
    
    common.inum = config.natoms[0];    
    common.inummax = config.natoms[0];    
    for (int i=1; i<common.nconfigs; i++)         
        if (config.natoms[i] > common.inummax)
            common.inummax = config.natoms[i];   
    
    for (int i=1; i<common.dim; i++) {       
        common.pbc[i] = app.pbc[i];
        common.boxoffset[i] = app.boxoffset[i];
    }
}

void implSetNeighborStruct(neighborstruct &nb, commonstruct &common, configstruct &config, Int ci)
{               
    Int dim = common.dim;    
    Int m = (dim==2) ? dim*4 : dim*8;         
    Int n = (dim==2) ? dim*9 : dim*27;         
        
    Int ne[dim], cellnum[3];
    dstype a[dim], b[dim], c[dim], cellsize[dim], s2rmap[dim*dim];        
    dstype refvertices[m], rbvertices[m], boxvertices[m], bbvertices[m];
    dstype pimages[n];
    
    dstype *simbox = &config.simbox[dim*dim*ci];    
    cpuSmallMatrixInverse(s2rmap, simbox, dim);       
    for (Int i=0; i<dim; i++) {
        a[i] = simbox[i];
        b[i] = simbox[dim+i];        
        c[i] = (dim > 2) ? simbox[2*dim+i] : 0.0;        
    }                    
    
    if (dim==2) {
        cpuBoundingBox2D(refvertices, rbvertices, boxvertices, bbvertices, 
                a, b, common.boxoffset, common.pbc);
        
        common.pnum = cpuPeriodicImages2D(pimages, a, b, common.pbc);                
    }
    else if (dim==3){
        cpuBoundingBox3D(refvertices, rbvertices, boxvertices, bbvertices, 
                a, b, c, common.boxoffset, common.pbc);
        
        common.pnum = cpuPeriodicImages3D(pimages, a, b, c, common.pbc);                
    }    
                
    dstype *smin = &rbvertices[0];
    dstype *smax = (dim == 2) ? &rbvertices[dim*2] :  &rbvertices[dim*6];        
    common.cnum = 1;
    for (Int i=0; i<dim; i++) {
        ne[i] = (Int) floor(1.0/(smax[i]-1.0));   
        ne[i] = (ne[i] < 1) ? 1 : ne[i];            
        cellnum[i] = ne[i] + 2;                
        cellsize[i] = 1.0/((dstype) ne[i]);
        common.cnum = common.cnum*cellnum[i];
    }
    if (dim==3) 
        cellnum[2] = 0;    
    dstype eta1[cellnum[0]+1], eta2[cellnum[1]+1], eta3[cellnum[2]+1];
    cpuMakeReferenceGrid(eta1, smin[0], smax[0], cellsize[0], cellnum[0]);
    cpuMakeReferenceGrid(eta2, smin[1], smax[1], cellsize[1], cellnum[1]);
    if (dim==3)
        cpuMakeReferenceGrid(eta3, smin[2], smax[2], cellsize[2], cellnum[2]);    
    
    nb.freememory(common.backend);        
    TemplateMalloc(&nb.cellnum, dim, common.backend);         
    TemplateMalloc(&nb.cellsize, dim, common.backend);         
    TemplateMalloc(&nb.a, dim, common.backend); 
    TemplateMalloc(&nb.b, dim, common.backend); 
    TemplateMalloc(&nb.c, dim, common.backend);    
    TemplateMalloc(&nb.s2rmap, dim*dim, common.backend);              
    TemplateMalloc(&nb.refvertices, m, common.backend); 
    TemplateMalloc(&nb.rbvertices,  m, common.backend); 
    TemplateMalloc(&nb.boxvertices, m, common.backend); 
    TemplateMalloc(&nb.bbvertices, m, common.backend);    
    TemplateMalloc(&nb.pimages, n, common.backend);    
    TemplateMalloc(&nb.eta1, cellnum[0]+1, common.backend);   
    TemplateMalloc(&nb.eta2, cellnum[1]+1, common.backend);   
    TemplateMalloc(&nb.eta3, cellnum[2]+1, common.backend);  
        
    if (common.backend <= 1) {
        for (Int i=0; i<dim; i++)
            nb.cellnum[i] = cellnum[i];
        ArrayCopy(nb.cellsize, cellsize, dim, common.backend);
        ArrayCopy(nb.a, a, dim, common.backend);
        ArrayCopy(nb.b, b, dim, common.backend);
        ArrayCopy(nb.c, c, dim, common.backend);
        ArrayCopy(nb.s2rmap, s2rmap, dim*dim, common.backend);
        ArrayCopy(nb.refvertices, refvertices, m, common.backend);
        ArrayCopy(nb.rbvertices, rbvertices, m, common.backend);
        ArrayCopy(nb.boxvertices, boxvertices, m, common.backend);
        ArrayCopy(nb.bbvertices, bbvertices, m, common.backend);
        ArrayCopy(nb.pimages, pimages, n, common.backend);
        ArrayCopy(nb.eta1, eta1, cellnum[0]+1, common.backend);
        ArrayCopy(nb.eta2, eta2, cellnum[1]+1, common.backend);
        ArrayCopy(nb.eta3, eta3, cellnum[2]+1, common.backend);
    }
    else if (common.backend==2) { // GPU
#ifdef HAVE_CUDA        
        CHECK( cudaMemcpy(nb.cellnum, cellnum, dim*sizeof(Int), cudaMemcpyHostToDevice ) );              
        CHECK( cudaMemcpy(nb.cellsize, cellsize, dim*sizeof(dstype), cudaMemcpyHostToDevice ) );              
        CHECK( cudaMemcpy(nb.a, a, dim*sizeof(dstype), cudaMemcpyHostToDevice ) );              
        CHECK( cudaMemcpy(nb.b, b, dim*sizeof(dstype), cudaMemcpyHostToDevice ) );              
        CHECK( cudaMemcpy(nb.c, c, dim*sizeof(dstype), cudaMemcpyHostToDevice ) );              
        CHECK( cudaMemcpy(nb.s2rmap, s2rmap, dim*dim*sizeof(dstype), cudaMemcpyHostToDevice ) );              
        CHECK( cudaMemcpy(nb.refvertices, refvertices, m*sizeof(dstype), cudaMemcpyHostToDevice ) );              
        CHECK( cudaMemcpy(nb.rbvertices, rbvertices, m*sizeof(dstype), cudaMemcpyHostToDevice ) );              
        CHECK( cudaMemcpy(nb.boxvertices, boxvertices, m*sizeof(dstype), cudaMemcpyHostToDevice ) );              
        CHECK( cudaMemcpy(nb.bbvertices, bbvertices, m*sizeof(dstype), cudaMemcpyHostToDevice ) );              
        CHECK( cudaMemcpy(nb.pimages, pimages, n*sizeof(dstype), cudaMemcpyHostToDevice ) );              
        CHECK( cudaMemcpy(nb.eta1, eta1, (cellnum[0]+1)*sizeof(dstype), cudaMemcpyHostToDevice ) );              
        CHECK( cudaMemcpy(nb.eta2, eta2, (cellnum[1]+1)*sizeof(dstype), cudaMemcpyHostToDevice ) );              
        CHECK( cudaMemcpy(nb.eta3, eta3, (cellnum[2]+1)*sizeof(dstype), cudaMemcpyHostToDevice ) );              
#endif        
    }           
}

    
void implNeighborList(neighborstruct &nb, commonstruct &common, appstruct &app, tempstruct &tmp, 
        dstype* x, Int *atomtype, Int inum)
{
    common.inum = inum;
    Int ns = min(inum, common.nab);
    common.nba = max((Int) floor(inum/ns), 16); // number of blocks
    Int na = round(inum/common.nba); // number of atoms per block    
    for (int i=0; i<=common.nba; i++)
        common.ablks[i] = i*na;
    common.ablks[common.nba] = min(common.ablks[common.nba], inum);
            
    Int dim = common.dim;
    Int cnum = common.cnum;
    Int pnum = common.pnum;
    Int jnum = common.jnum;    
    Int ntype = common.natomtypes;
    Int gnum, anum;

    Int *inside = &tmp.intmem[0]; // inum*pnum
    Int *glistnum = &tmp.intmem[inum*common.pnum]; // inum
    Int *glistnumsum = &tmp.intmem[inum*common.pnum + inum]; // inum+1
    
    if (dim==2) {
        cpuAtomList2D(nb.alist, inside, glistnumsum, glistnum,  x, nb.pimages, nb.rbvertices, 
            nb.s2rmap, inum, pnum, dim);
        
        gnum = IntArrayGetValueAtIndex(glistnumsum, inum, common.backend);        
        anum = inum + common.gnum;
        common.gnum = gnum;
        common.anum = anum;
                
        Int *clist = &tmp.intmem[0]; //anum
        Int *c2alist = &tmp.intmem[anum]; // anum
        Int *c2anum = &tmp.intmem[2*anum]; // cnum
        Int *c2anumsum = &tmp.intmem[2*anum+cnum]; // cnum+1
        cpuCellList2D(clist, x, nb.eta1, nb.eta2, nb.eta3, nb.s2rmap, nb.cellnum, anum, dim);                
        cpuCell2AtomList(c2alist, c2anumsum, c2anum, clist, anum, cnum);
                
        if (common.cutofftype) {
            if (common.neightype)
                cpuHalfNeighborList2D(nb.neighlist, nb.neighnum, x, app.cutoffellipsoid, nb.alist, clist, 
                        c2alist, c2anumsum, nb.cellnum, inum, jnum, dim);                
            else
                cpuFullNeighborList2D(nb.neighlist, nb.neighnum, x, app.cutoffellipsoid, nb.alist, clist, 
                        c2alist, c2anumsum, nb.cellnum, inum, jnum, dim);                                                                            
        }
        else {
            if (common.neightype)
                cpuHalfNeighborList2D(nb.neighlist, nb.neighnum, x, app.rcutsq, nb.alist, clist, 
                        c2alist, c2anumsum, atomtype, nb.cellnum, ntype, inum, jnum, dim);                
            else
                cpuFullNeighborList2D(nb.neighlist, nb.neighnum, x, app.rcutsq, nb.alist, clist, 
                        c2alist, c2anumsum, atomtype, nb.cellnum, ntype, inum, jnum, dim);                            
        }               
    }
    else {
        cpuAtomList3D(nb.alist, inside, glistnumsum, glistnum,  x, nb.pimages, nb.rbvertices, 
            nb.s2rmap, inum, pnum, dim);
        
        gnum = IntArrayGetValueAtIndex(glistnumsum, inum, common.backend);        
        anum = inum + common.gnum;
        common.gnum = gnum;
        common.anum = anum;
                
        Int *clist = &tmp.intmem[0]; //anum
        Int *c2alist = &tmp.intmem[anum]; // anum
        Int *c2anum = &tmp.intmem[2*anum]; // cnum
        Int *c2anumsum = &tmp.intmem[2*anum+cnum]; // cnum+1
        cpuCellList3D(clist, x, nb.eta1, nb.eta2, nb.eta3, nb.s2rmap, nb.cellnum, anum, dim);                
        cpuCell2AtomList(c2alist, c2anumsum, c2anum, clist, anum, cnum);
        
        if (common.cutofftype) {
            if (common.neightype)
                cpuHalfNeighborList3D(nb.neighlist, nb.neighnum, x, app.cutoffellipsoid, nb.alist, clist, 
                        c2alist, c2anumsum, nb.cellnum, inum, jnum, dim);                
            else
                cpuFullNeighborList3D(nb.neighlist, nb.neighnum, x, app.cutoffellipsoid, nb.alist, clist, 
                        c2alist, c2anumsum, nb.cellnum, inum, jnum, dim);                                                                            
        }
        else {
            if (common.neightype)
                cpuHalfNeighborList3D(nb.neighlist, nb.neighnum, x, app.rcutsq, nb.alist, clist, 
                        c2alist, c2anumsum, atomtype, nb.cellnum, ntype, inum, jnum, dim);                
            else
                cpuFullNeighborList3D(nb.neighlist, nb.neighnum, x, app.rcutsq, nb.alist, clist, 
                        c2alist, c2anumsum, atomtype, nb.cellnum, ntype, inum, jnum, dim);                            
        }        
    }
    
    cpuCumsum(nb.neighnumsum, nb.neighnum, inum+1);    
}


// void implNeighborList(neighborstruct &nb, commonstruct &common, appstruct &app, tempstruct &tmp, dstype* x, Int inum)
// {
//     common.inum = inum;
//     Int dim = common.dim;
//     Int cnum = common.cnum;
//     Int pnum = common.pnum;
//     Int jnum = common.jnum;    
//     Int gnum, anum;
// 
//     Int *inside = &tmp.intmem[0]; // inum*pnum
//     Int *glistnum = &tmp.intmem[inum*common.pnum]; // inum
//     Int *glistnumsum = &tmp.intmem[inum*common.pnum + inum]; // inum+1
//     
//     if (dim==2) {
//         cpuAtomList2D(nb.alist, inside, glistnumsum, glistnum,  x, nb.pimages, nb.rbvertices, 
//             nb.s2rmap, inum, pnum, dim);
//         
//         gnum = IntArrayGetValueAtIndex(glistnumsum, inum, common.backend);        
//         anum = inum + common.gnum;
//         common.gnum = gnum;
//         common.anum = anum;
//                 
//         Int *clist = &tmp.intmem[0]; //anum
//         Int *c2alist = &tmp.intmem[anum]; // anum
//         Int *c2anum = &tmp.intmem[2*anum]; // cnum
//         Int *c2anumsum = &tmp.intmem[2*anum+cnum]; // cnum+1
//         cpuCellList2D(clist, x, nb.eta1, nb.eta2, nb.eta3, nb.s2rmap, nb.cellnum, anum, dim);                
//         cpuCell2AtomList(c2alist, c2anumsum, c2anum, clist, anum, cnum);
//         
//         if (common.neightype)
//             cpuHalfNeighborList2D(nb.neighlist, nb.neighnum, x, app.cutoffellipsoid, nb.alist, clist, 
//                     c2alist, c2anumsum, nb.cellnum, inum, jnum, dim);                
//         else
//             cpuFullNeighborList2D(nb.neighlist, nb.neighnum, x, app.cutoffellipsoid, nb.alist, clist, 
//                     c2alist, c2anumsum, nb.cellnum, inum, jnum, dim);                                            
//     }
//     else {
//         cpuAtomList3D(nb.alist, inside, glistnumsum, glistnum,  x, nb.pimages, nb.rbvertices, 
//             nb.s2rmap, inum, pnum, dim);
//         
//         gnum = IntArrayGetValueAtIndex(glistnumsum, inum, common.backend);        
//         anum = inum + common.gnum;
//         common.gnum = gnum;
//         common.anum = anum;
//                 
//         Int *clist = &tmp.intmem[0]; //anum
//         Int *c2alist = &tmp.intmem[anum]; // anum
//         Int *c2anum = &tmp.intmem[2*anum]; // cnum
//         Int *c2anumsum = &tmp.intmem[2*anum+cnum]; // cnum+1
//         cpuCellList3D(clist, x, nb.eta1, nb.eta2, nb.eta3, nb.s2rmap, nb.cellnum, anum, dim);                
//         cpuCell2AtomList(c2alist, c2anumsum, c2anum, clist, anum, cnum);
//         
//         if (common.neightype)
//             cpuHalfNeighborList3D(nb.neighlist, nb.neighnum, x, app.cutoffellipsoid, nb.alist, clist, 
//                     c2alist, c2anumsum, nb.cellnum, inum, jnum, dim);                
//         else
//             cpuFullNeighborList3D(nb.neighlist, nb.neighnum, x, app.cutoffellipsoid, nb.alist, clist, 
//                     c2alist, c2anumsum, nb.cellnum, inum, jnum, dim);                                                    
//     }
// }        
    

#endif

