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
    readarray(in, &app.boxoffset, app.nsize[4]);  
    app.atomnumber = readiarrayfromdouble(in, app.nsize[5]);            
    readarray(in, &app.atommass, app.nsize[6]);    
    readarray(in, &app.atomcharge, app.nsize[7]);    
    readarray(in, &app.simulaparam, app.nsize[8]);   
    readarray(in, &app.solversparam, app.nsize[9]);       
    readarray(in, &app.eta, app.nsize[10]);       
    app.kappa = readiarrayfromdouble(in, app.nsize[11]);    
    readarray(in, &app.muml, app.nsize[12]);       
    readarray(in, &app.mu1a, app.nsize[13]);       
    readarray(in, &app.mu1b, app.nsize[14]);       
    readarray(in, &app.mu2a, app.nsize[15]);       
    readarray(in, &app.mu2b, app.nsize[16]);       
    readarray(in, &app.mu2c, app.nsize[17]);       
    readarray(in, &app.mu3a, app.nsize[18]);       
    readarray(in, &app.mu3b, app.nsize[19]);       
    readarray(in, &app.mu3c, app.nsize[20]);       
    readarray(in, &app.mu4a, app.nsize[21]);       
    readarray(in, &app.mu4b, app.nsize[22]);       
    app.pot1a = readiarrayfromdouble(in, app.nsize[23]);    
    app.pot1b = readiarrayfromdouble(in, app.nsize[24]);    
    app.pot2a = readiarrayfromdouble(in, app.nsize[25]);    
    app.pot2b = readiarrayfromdouble(in, app.nsize[26]);    
    app.pot2c = readiarrayfromdouble(in, app.nsize[27]);    
    app.pot3a = readiarrayfromdouble(in, app.nsize[28]);    
    app.pot3b = readiarrayfromdouble(in, app.nsize[29]);    
    app.pot3c = readiarrayfromdouble(in, app.nsize[30]);    
    app.pot4a = readiarrayfromdouble(in, app.nsize[31]);    
    app.pot4b = readiarrayfromdouble(in, app.nsize[32]);        
    readarray(in, &app.rcutsqml, app.nsize[33]);       
    readarray(in, &app.rcutsq2a, app.nsize[34]);       
    readarray(in, &app.rcutsq2b, app.nsize[35]);       
    readarray(in, &app.rcutsq2c, app.nsize[36]);       
    readarray(in, &app.rcutsq3a, app.nsize[37]);       
    readarray(in, &app.rcutsq3b, app.nsize[38]);       
    readarray(in, &app.rcutsq3c, app.nsize[39]);       
    readarray(in, &app.rcutsq4a, app.nsize[40]);       
    readarray(in, &app.rcutsq4b, app.nsize[41]);       
    readarray(in, &app.rcutsq, app.nsize[42]);       
    app.atom1b = readiarrayfromdouble(in, app.nsize[43]);    
    app.atom2b = readiarrayfromdouble(in, app.nsize[44]);    
    app.atom2c = readiarrayfromdouble(in, app.nsize[45]);    
    app.atom3b = readiarrayfromdouble(in, app.nsize[46]);    
    app.atom3c = readiarrayfromdouble(in, app.nsize[47]);    
    app.atom4b = readiarrayfromdouble(in, app.nsize[48]);        
    
//     printArray2D(app.nsize, 1, app.lsize[0], 0);       
//     printArray2D(app.atomnumber, 1, app.nsize[5], 0);   
//     printArray2D(app.atommass, 1, app.nsize[6], 0);   
//     printArray2D(app.atomcharge, 1, app.nsize[7], 0);   
//     printArray2D(app.simulaparam, 1, app.nsize[8], 0);       
//     cout<<app.rcutsqml[0]<<endl;
//     cout<<app.rcutsq[0]<<endl;
//     error("here");
    
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
    TemplateMalloc(&app.boxoffset, happ.nsize[4], backend);  
    TemplateMalloc(&app.atomnumber, happ.nsize[5], backend); 
    TemplateMalloc(&app.atommass, happ.nsize[6], backend); 
    TemplateMalloc(&app.atomcharge, happ.nsize[7], backend); 
    TemplateMalloc(&app.simulaparam, happ.nsize[8], backend);     
    TemplateMalloc(&app.solversparam, happ.nsize[9], backend);
    TemplateMalloc(&app.eta, happ.nsize[10], backend); 
    TemplateMalloc(&app.kappa, happ.nsize[11], backend); 
    TemplateMalloc(&app.muml, happ.nsize[12], backend); 
    TemplateMalloc(&app.mu1a, happ.nsize[13], backend); 
    TemplateMalloc(&app.mu1b, happ.nsize[14], backend); 
    TemplateMalloc(&app.mu2a, happ.nsize[15], backend); 
    TemplateMalloc(&app.mu2a, happ.nsize[16], backend); 
    TemplateMalloc(&app.mu2c, happ.nsize[17], backend); 
    TemplateMalloc(&app.mu3a, happ.nsize[18], backend); 
    TemplateMalloc(&app.mu3b, happ.nsize[19], backend); 
    TemplateMalloc(&app.mu3c, happ.nsize[20], backend); 
    TemplateMalloc(&app.mu4a, happ.nsize[21], backend); 
    TemplateMalloc(&app.mu4b, happ.nsize[22], backend); 
    TemplateMalloc(&app.pot1a, happ.nsize[23], backend); 
    TemplateMalloc(&app.pot1b, happ.nsize[24], backend); 
    TemplateMalloc(&app.pot2a, happ.nsize[25], backend); 
    TemplateMalloc(&app.pot2b, happ.nsize[26], backend); 
    TemplateMalloc(&app.pot2c, happ.nsize[27], backend); 
    TemplateMalloc(&app.pot3a, happ.nsize[28], backend); 
    TemplateMalloc(&app.pot3b, happ.nsize[29], backend); 
    TemplateMalloc(&app.pot3c, happ.nsize[30], backend); 
    TemplateMalloc(&app.pot4a, happ.nsize[31], backend); 
    TemplateMalloc(&app.pot4b, happ.nsize[32], backend);     
    TemplateMalloc(&app.rcutsqml, happ.nsize[33], backend); 
    TemplateMalloc(&app.rcutsq2a, happ.nsize[34], backend); 
    TemplateMalloc(&app.rcutsq2b, happ.nsize[35], backend); 
    TemplateMalloc(&app.rcutsq2c, happ.nsize[36], backend); 
    TemplateMalloc(&app.rcutsq3a, happ.nsize[37], backend); 
    TemplateMalloc(&app.rcutsq3b, happ.nsize[38], backend); 
    TemplateMalloc(&app.rcutsq3c, happ.nsize[39], backend); 
    TemplateMalloc(&app.rcutsq4a, happ.nsize[40], backend); 
    TemplateMalloc(&app.rcutsq4b, happ.nsize[41], backend); 
    TemplateMalloc(&app.rcutsq, happ.nsize[42], backend);         
    TemplateMalloc(&app.atom1b, happ.nsize[43], backend); 
    TemplateMalloc(&app.atom2b, happ.nsize[44], backend); 
    TemplateMalloc(&app.atom2c, happ.nsize[45], backend); 
    TemplateMalloc(&app.atom3b, happ.nsize[46], backend); 
    TemplateMalloc(&app.atom3c, happ.nsize[47], backend); 
    TemplateMalloc(&app.atom4b, happ.nsize[48], backend);                
    
    if (backend==2) { // GPU
#ifdef HAVE_CUDA        
        CHECK( cudaMemcpy(app.nsize, happ.nsize, happ.lsize[0]*sizeof(Int), cudaMemcpyHostToDevice ) );              
        CHECK( cudaMemcpy(app.ndims, happ.ndims, happ.nsize[0]*sizeof(Int), cudaMemcpyHostToDevice ) );              
        CHECK( cudaMemcpy(app.flags, happ.flags, happ.nsize[1]*sizeof(Int), cudaMemcpyHostToDevice ) );              
        CHECK( cudaMemcpy(app.bcs, happ.bcs, happ.nsize[2]*sizeof(Int), cudaMemcpyHostToDevice ) );              
        CHECK( cudaMemcpy(app.pbc, happ.pbc, happ.nsize[3]*sizeof(Int), cudaMemcpyHostToDevice ) );              
        CHECK( cudaMemcpy(app.boxoffset, happ.boxoffset, happ.nsize[4]*sizeof(dstype), cudaMemcpyHostToDevice ) );                                      
        CHECK( cudaMemcpy(app.atomnumber, happ.atomnumber, happ.nsize[5]*sizeof(Int), cudaMemcpyHostToDevice ) );                         
        CHECK( cudaMemcpy(app.atommass, happ.atommass, happ.nsize[6]*sizeof(dstype), cudaMemcpyHostToDevice ) );              
        CHECK( cudaMemcpy(app.atomcharge, happ.atomcharge, happ.nsize[7]*sizeof(dstype), cudaMemcpyHostToDevice ) );                      
        CHECK( cudaMemcpy(app.simulaparam, happ.simulaparam, happ.nsize[8]*sizeof(dstype), cudaMemcpyHostToDevice ) );              
        CHECK( cudaMemcpy(app.solversparam, happ.solversparam, happ.nsize[9]*sizeof(dstype), cudaMemcpyHostToDevice ) );              
        CHECK( cudaMemcpy(app.eta, happ.eta, happ.nsize[10]*sizeof(dstype), cudaMemcpyHostToDevice ) );              
        CHECK( cudaMemcpy(app.kappa, happ.kappa, happ.nsize[11]*sizeof(Int), cudaMemcpyHostToDevice ) );   
        CHECK( cudaMemcpy(app.muml, happ.muml, happ.nsize[12]*sizeof(dstype), cudaMemcpyHostToDevice ) );              
        CHECK( cudaMemcpy(app.mu1a, happ.mu1a, happ.nsize[13]*sizeof(dstype), cudaMemcpyHostToDevice ) );              
        CHECK( cudaMemcpy(app.mu1b, happ.mu1b, happ.nsize[14]*sizeof(dstype), cudaMemcpyHostToDevice ) );              
        CHECK( cudaMemcpy(app.mu2a, happ.mu2a, happ.nsize[15]*sizeof(dstype), cudaMemcpyHostToDevice ) );              
        CHECK( cudaMemcpy(app.mu2b, happ.mu2b, happ.nsize[16]*sizeof(dstype), cudaMemcpyHostToDevice ) );              
        CHECK( cudaMemcpy(app.mu2c, happ.mu2c, happ.nsize[17]*sizeof(dstype), cudaMemcpyHostToDevice ) );              
        CHECK( cudaMemcpy(app.mu3a, happ.mu3a, happ.nsize[18]*sizeof(dstype), cudaMemcpyHostToDevice ) );              
        CHECK( cudaMemcpy(app.mu3b, happ.mu3b, happ.nsize[19]*sizeof(dstype), cudaMemcpyHostToDevice ) );              
        CHECK( cudaMemcpy(app.mu3c, happ.mu3c, happ.nsize[20]*sizeof(dstype), cudaMemcpyHostToDevice ) );              
        CHECK( cudaMemcpy(app.mu4a, happ.mu4a, happ.nsize[21]*sizeof(dstype), cudaMemcpyHostToDevice ) );              
        CHECK( cudaMemcpy(app.mu4b, happ.mu4b, happ.nsize[22]*sizeof(dstype), cudaMemcpyHostToDevice ) );                      
        CHECK( cudaMemcpy(app.pot1a, happ.pot1a, happ.nsize[23]*sizeof(Int), cudaMemcpyHostToDevice ) );              
        CHECK( cudaMemcpy(app.pot1b, happ.pot1b, happ.nsize[24]*sizeof(Int), cudaMemcpyHostToDevice ) );              
        CHECK( cudaMemcpy(app.pot2a, happ.pot2a, happ.nsize[25]*sizeof(Int), cudaMemcpyHostToDevice ) );              
        CHECK( cudaMemcpy(app.pot2b, happ.pot2b, happ.nsize[26]*sizeof(Int), cudaMemcpyHostToDevice ) );              
        CHECK( cudaMemcpy(app.pot2c, happ.pot2c, happ.nsize[27]*sizeof(Int), cudaMemcpyHostToDevice ) );              
        CHECK( cudaMemcpy(app.pot3a, happ.pot3a, happ.nsize[28]*sizeof(Int), cudaMemcpyHostToDevice ) );              
        CHECK( cudaMemcpy(app.pot3b, happ.pot3b, happ.nsize[29]*sizeof(Int), cudaMemcpyHostToDevice ) );              
        CHECK( cudaMemcpy(app.pot3c, happ.pot3c, happ.nsize[30]*sizeof(Int), cudaMemcpyHostToDevice ) );              
        CHECK( cudaMemcpy(app.pot4a, happ.pot4a, happ.nsize[31]*sizeof(Int), cudaMemcpyHostToDevice ) );              
        CHECK( cudaMemcpy(app.pot4b, happ.pot4b, happ.nsize[32]*sizeof(Int), cudaMemcpyHostToDevice ) );                      
        CHECK( cudaMemcpy(app.rcutsqml, happ.rcutsqml, happ.nsize[33]*sizeof(dstype), cudaMemcpyHostToDevice ) );              
        CHECK( cudaMemcpy(app.rcutsq2a, happ.rcutsq2a, happ.nsize[34]*sizeof(dstype), cudaMemcpyHostToDevice ) );              
        CHECK( cudaMemcpy(app.rcutsq2b, happ.rcutsq2b, happ.nsize[35]*sizeof(dstype), cudaMemcpyHostToDevice ) );              
        CHECK( cudaMemcpy(app.rcutsq2c, happ.rcutsq2c, happ.nsize[36]*sizeof(dstype), cudaMemcpyHostToDevice ) );              
        CHECK( cudaMemcpy(app.rcutsq3a, happ.rcutsq3a, happ.nsize[37]*sizeof(dstype), cudaMemcpyHostToDevice ) );              
        CHECK( cudaMemcpy(app.rcutsq3b, happ.rcutsq3b, happ.nsize[38]*sizeof(dstype), cudaMemcpyHostToDevice ) );              
        CHECK( cudaMemcpy(app.rcutsq3c, happ.rcutsq3c, happ.nsize[39]*sizeof(dstype), cudaMemcpyHostToDevice ) );              
        CHECK( cudaMemcpy(app.rcutsq4a, happ.rcutsq4a, happ.nsize[40]*sizeof(dstype), cudaMemcpyHostToDevice ) );              
        CHECK( cudaMemcpy(app.rcutsq4b, happ.rcutsq4b, happ.nsize[41]*sizeof(dstype), cudaMemcpyHostToDevice ) );                      
        CHECK( cudaMemcpy(app.rcutsq, happ.rcutsq, happ.nsize[42]*sizeof(dstype), cudaMemcpyHostToDevice ) );                          
        CHECK( cudaMemcpy(app.atom1b, happ.atom1b, happ.nsize[43]*sizeof(Int), cudaMemcpyHostToDevice ) );                      
        CHECK( cudaMemcpy(app.atom2b, happ.atom2b, happ.nsize[44]*sizeof(Int), cudaMemcpyHostToDevice ) );              
        CHECK( cudaMemcpy(app.atom2c, happ.atom2c, happ.nsize[45]*sizeof(Int), cudaMemcpyHostToDevice ) );                      
        CHECK( cudaMemcpy(app.atom3b, happ.atom3b, happ.nsize[46]*sizeof(Int), cudaMemcpyHostToDevice ) );              
        CHECK( cudaMemcpy(app.atom3c, happ.atom3c, happ.nsize[47]*sizeof(Int), cudaMemcpyHostToDevice ) );                      
        CHECK( cudaMemcpy(app.atom4b, happ.atom4b, happ.nsize[48]*sizeof(Int), cudaMemcpyHostToDevice ) );                                      
#endif        
    }        
}

void implReadConfigStruct(configstruct &config, string filein, Int mpiprocs, Int mpirank, Int backend)
{
    string filename;
    if (mpiprocs > 1) {
        Int filenumber = mpirank+1; //file number     
        filename = filein + "config" + NumberToString(filenumber) + ".bin";                    
    }
    else
        filename = filein + "config" + ".bin";                    
    
    // Open file to read
    ifstream in(filename.c_str(), ios::in | ios::binary);
    
    if (!in) 
        error("Unable to open file " + filename);
       
    if (mpirank==0)
        printf("Read config struct from a binary file...\n");   
    
    /* Read data to config structure */            
    config.lsize = readiarrayfromdouble(in, 1);
    config.nsize = readiarrayfromdouble(in, config.lsize[0]);
    config.ndims = readiarrayfromdouble(in, config.nsize[0]);
    config.natoms = readiarrayfromdouble(in, config.nsize[1]);        
    readarray(in, &config.a, config.nsize[2]);   
    readarray(in, &config.b, config.nsize[3]);   
    readarray(in, &config.c, config.nsize[4]);   
    readarray(in, &config.e, config.nsize[5]);   
    config.t = readiarrayfromdouble(in, config.nsize[6]);
    readarray(in, &config.x, config.nsize[7]);   
    readarray(in, &config.v, config.nsize[8]);   
    readarray(in, &config.f, config.nsize[9]);   
    readarray(in, &config.q, config.nsize[10]);       
            
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
    common.natomtypes = app.ndims[8]; // number of atom types
    common.nmoletypes = app.ndims[9]; // number of molecule types
    common.jnum = app.ndims[10]; // % maximum number of neighbors allowed
    
    common.descriptor = app.flags[0];   // descriptor flag: 0 -> Spherical Harmonics Bessel
    common.spectrum = app.flags[1];     // spectrum flag: 0-> power spectrum, 1-> bispectrum, 2-> both power and bispectrum 
    common.training = app.flags[2];     // 0 -> no training, 1 -> Linear regression, 2 -> Gaussian process, 3 -> Neural net
    common.runMD = app.flags[3];        // 0 no MD simulation, 1 -> run MD simulation
    common.potential = app.flags[4];    // 0 -> empirical potential, 1 -> ML potential, 2 -> combined potential        
    common.neighpair = app.flags[5];  // 0 -> full neighbor list, 1 -> half neighbor list for pair potentials
    common.energycal = app.flags[6];    // turns energy calculation on or off
    common.forcecal = app.flags[7];     // turns force calculation on or off
    common.stresscal = app.flags[8];   // turns stress calculation on or off
    common.neighcell = app.flags[9];   // 0 -> O(N^2) algorithm, 1 -> Cell-linked list algorithm to form neighbor list
    common.decomposition = app.flags[10]; // 0 -> force decomposition, 1 -> atom decomposition
    
    common.time = app.simulaparam[0];
    common.dt = app.simulaparam[1];
    
    common.nflags = app.nsize[1];
    common.nsimulaparam = app.nsize[8];   
    common.nsolversparam = app.nsize[9];       
    common.neta = app.nsize[10];  
    common.nkappa = app.nsize[11];  
    common.nmuml = app.nsize[12];  
    common.nmu1a = app.nsize[13];  
    common.nmu1b = app.nsize[14];  
    common.nmu2a = app.nsize[15];  
    common.nmu2b = app.nsize[16];  
    common.nmu2c = app.nsize[17];  
    common.nmu3a = app.nsize[18];  
    common.nmu3b = app.nsize[19];  
    common.nmu3c = app.nsize[20];  
    common.nmu4a = app.nsize[21];  
    common.nmu4b = app.nsize[22];      
    common.npot1a = app.nsize[23];  
    common.npot1b = app.nsize[24];  
    common.npot2a = app.nsize[25];  
    common.npot2b = app.nsize[26];  
    common.npot2c = app.nsize[27];  
    common.npot3a = app.nsize[28];  
    common.npot3b = app.nsize[29];  
    common.npot3c = app.nsize[30];  
    common.npot4a = app.nsize[31];  
    common.npot4b = app.nsize[32];                  
    common.natom1b = app.nsize[43];  
    common.natom2b = app.nsize[44];  
    common.natom2c = app.nsize[45];  
    common.natom3b = app.nsize[46];  
    common.natom3c = app.nsize[47];  
    common.natom4b = app.nsize[48];  
        
    common.nconfigs = config.nsize[1]; // number of configurations
    common.ne = config.nsize[5];
    common.nt = config.nsize[6];
    common.nx = config.nsize[7];
    common.nq = config.nsize[8];
    common.nv = config.nsize[9];
    common.nf = config.nsize[10];
    common.ncq = common.nq/(common.nt);        
            
    common.inum = config.natoms[0];    
    common.inummax = config.natoms[0];    
    for (int i=1; i<common.nconfigs; i++)         
        if (config.natoms[i] > common.inummax)
            common.inummax = config.natoms[i];   
    
    for (int i=0; i<common.dim; i++) {       
        common.pbc[i] = app.pbc[i];
        common.boxoffset[i] = app.boxoffset[i];
    }
    
    TemplateMalloc(&common.pot1a, app.nsize[23], 0); 
    TemplateMalloc(&common.pot1b, app.nsize[24], 0); 
    TemplateMalloc(&common.pot2a, app.nsize[25], 0); 
    TemplateMalloc(&common.pot2b, app.nsize[26], 0); 
    TemplateMalloc(&common.pot2c, app.nsize[27], 0); 
    TemplateMalloc(&common.pot3a, app.nsize[28], 0); 
    TemplateMalloc(&common.pot3b, app.nsize[29], 0); 
    TemplateMalloc(&common.pot3c, app.nsize[30], 0); 
    TemplateMalloc(&common.pot4a, app.nsize[31], 0); 
    TemplateMalloc(&common.pot4b, app.nsize[32], 0);         
    TemplateMalloc(&common.atom1b, app.nsize[43], 0); 
    TemplateMalloc(&common.atom2b, app.nsize[44], 0); 
    TemplateMalloc(&common.atom2c, app.nsize[45], 0); 
    TemplateMalloc(&common.atom3b, app.nsize[46], 0); 
    TemplateMalloc(&common.atom3c, app.nsize[47], 0); 
    TemplateMalloc(&common.atom4b, app.nsize[48], 0);        
    
    ArrayCopy(common.pot1a, app.pot1a, app.nsize[23], 0);
    ArrayCopy(common.pot1b, app.pot1b, app.nsize[24], 0);
    ArrayCopy(common.pot2a, app.pot2a, app.nsize[25], 0);
    ArrayCopy(common.pot2b, app.pot2b, app.nsize[26], 0);
    ArrayCopy(common.pot2c, app.pot2c, app.nsize[27], 0);
    ArrayCopy(common.pot3a, app.pot3a, app.nsize[28], 0);
    ArrayCopy(common.pot3b, app.pot3b, app.nsize[29], 0);
    ArrayCopy(common.pot3c, app.pot3c, app.nsize[30], 0);
    ArrayCopy(common.pot4a, app.pot4a, app.nsize[31], 0);
    ArrayCopy(common.pot4b, app.pot4b, app.nsize[32], 0);    
    ArrayCopy(common.atom1b, app.atom1b, app.nsize[43], 0);
    ArrayCopy(common.atom2b, app.atom2b, app.nsize[44], 0);
    ArrayCopy(common.atom2c, app.atom2c, app.nsize[45], 0);
    ArrayCopy(common.atom3b, app.atom3b, app.nsize[46], 0);
    ArrayCopy(common.atom3c, app.atom3c, app.nsize[47], 0);
    ArrayCopy(common.atom4b, app.atom4b, app.nsize[48], 0);
}

void implSetAtomBlocks(commonstruct &common) 
{
    Int inum = common.inum;
    Int ns = min(inum, common.nab);  // number of atoms per block
    common.nba = min((Int) round(inum/ns), 16);     // number of blocks    
    Int na = (Int) round(inum/common.nba); // number of atoms per block    
        
    for (int i=0; i<=common.nba; i++)
        common.ablks[i] = i*na;
    common.ablks[common.nba] = inum;
}

void implGetConfiguration(Int *atomtype, dstype *x, commonstruct &common, configstruct &config, Int ci)
{
    Int inum = config.natoms[ci];    
    Int start = config.natomssum[ci];    
    common.inum = inum;    
    
    if (common.backend <= 1) {               
        ArrayCopy(x, &config.x[common.dim*start], common.dim*inum, common.backend);
        ArrayCopy(atomtype, &config.t[start], inum, common.backend);
    }
    else {
#ifdef HAVE_CUDA                
        CHECK( cudaMemcpy(atomtype, &config.t[start], inum*sizeof(Int), cudaMemcpyHostToDevice ) );                      
        CHECK( cudaMemcpy(x, &config.x[common.dim*start], common.dim*inum*sizeof(dstype), cudaMemcpyHostToDevice ) );              
#endif        
    }        
}

void implGetAtomtypes(Int *atomtype, commonstruct &common, configstruct &config, Int ci)
{
    common.inum = config.natoms[ci];    
    Int start = config.natomssum[ci];        
    if (common.backend <= 1) {               
        ArrayCopy(atomtype, &config.t[start], common.inum, common.backend);
    }
    else {
#ifdef HAVE_CUDA                
        CHECK( cudaMemcpy(atomtype, &config.t[start], common.inum*sizeof(Int), cudaMemcpyHostToDevice ) );                      
#endif        
    }        
}

void implGetPositions(dstype *x, commonstruct &common, configstruct &config, Int ci)
{
    common.inum = config.natoms[ci];    
    Int start = config.natomssum[ci];      
    if (common.backend <= 1) {
        ArrayCopy(x, &config.x[common.dim*start], common.dim*common.inum, common.backend);    
    }
    else {
#ifdef HAVE_CUDA                        
        CHECK( cudaMemcpy(x, &config.x[common.dim*start], common.dim*inum*sizeof(dstype), cudaMemcpyHostToDevice ) );              
#endif    
    }
}

void implGetVelocities(dstype *v, commonstruct &common, configstruct &config, Int ci)
{
    Int start = config.natomssum[ci];      
    if (common.backend <= 1) {
        ArrayCopy(v, &config.v[common.dim*start], common.dim*common.inum, common.backend);    
    }
    else {
#ifdef HAVE_CUDA                        
        CHECK( cudaMemcpy(v, &config.v[common.dim*start], common.dim*inum*sizeof(dstype), cudaMemcpyHostToDevice ) );              
#endif    
    }
}

void implGetForces(dstype *f, commonstruct &common, configstruct &config, Int ci)
{
    Int start = config.natomssum[ci];        
    if (common.backend <= 1) {
        ArrayCopy(f, &config.f[common.dim*start], common.dim*common.inum, common.backend);    
    }
    else {
#ifdef HAVE_CUDA                                
        CHECK( cudaMemcpy(f, &config.f[common.dim*start], common.dim*inum*sizeof(dstype), cudaMemcpyHostToDevice ) );              
#endif    
    }
}

void implGetCharges(dstype *q, commonstruct &common, configstruct &config, Int ci)
{
    Int start = config.natomssum[ci]; 
    if (common.backend <= 1) {
        ArrayCopy(q, &config.q[common.dim*start], common.ncq*common.inum, common.backend);    
    }
    else {
#ifdef HAVE_CUDA                                
        CHECK( cudaMemcpy(q, &config.q[common.dim*start], common.ncq*inum*sizeof(dstype), cudaMemcpyHostToDevice ) );              
#endif            
    }        
}

void implGetEnergy(dstype *e, commonstruct &common, configstruct &config, Int ci)
{    
    if (common.backend <= 1) {
        ArrayCopy(e, &config.e[ci], 1, common.backend);    
    }
    else {
#ifdef HAVE_CUDA                                
        CHECK( cudaMemcpy(e, &config.e[ci], sizeof(dstype), cudaMemcpyHostToDevice ) );              
#endif                    
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
    
    dstype simbox[dim*dim];    
    for (Int i=0; i<dim; i++) {
        simbox[i] = config.a[i + dim*dim*ci];
        simbox[dim+i] = config.b[i + dim*dim*ci];
        if (dim > 2)
            simbox[2*dim+i] = config.c[i + dim*dim*ci];
    }                    
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
        ne[i] = (Int) round(1.0/(smax[i]-1.0));   
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
            
    common.inum = config.natoms[ci];    
    TemplateMalloc(&nb.atomtype, common.inum, common.backend);  
    implGetAtomtypes(nb.atomtype, common, config, ci);                
    
    // need to check size of nb.alist
    TemplateMalloc(&nb.alist, 7*common.inum, common.backend);  
    TemplateMalloc(&nb.neighlist, common.inum*common.jnum, common.backend);  
    TemplateMalloc(&nb.neighnum, common.inum, common.backend);      
    //TemplateMalloc(&nb.neighnumsum, common.inum+1, common.backend);      
    for (Int i=0; i<common.inum*common.jnum; i++)
        nb.neighlist[i] = 0;
        
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

void implSetTempStruct(tempstruct & tmp, commonstruct &common)
{
    TemplateMalloc(&tmp.intmem, 40*(2+common.pnum)*common.inum+1, common.backend);  
    TemplateMalloc(&tmp.tmpmem, 100*common.inum, common.backend);  
}

void implSetSysStruct(sysstruct & sys, commonstruct &common, configstruct &config, Int ci)
{    
    TemplateMalloc(&sys.e, common.inum, common.backend);  
    TemplateMalloc(&sys.f, common.dim*common.inum, common.backend);  
    
    if (common.nx>0) {
        // need to check size of sys.x
        TemplateMalloc(&sys.x, 7*common.dim*common.inum, common.backend);  
        implGetPositions(sys.x, common, config, ci);
    }
    
    if (common.nv > 0) {
        TemplateMalloc(&sys.v, common.dim*common.inum, common.backend);  
        implGetVelocities(sys.v, common, config, ci);
    }
    
    if (common.nq > 0) {
        TemplateMalloc(&sys.q, common.ncq*common.inum, common.backend);  
        implGetCharges(sys.q, common, config, ci);           
    }
}

void implNeighborList(neighborstruct &nb, commonstruct &common, appstruct &app, tempstruct &tmp, 
        dstype* x, Int inum)
{
            
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
        anum = inum + gnum;
        common.inum = inum;
        common.gnum = gnum;
        common.anum = anum;
 
        if (common.neighcell == 0) { // O(N^2) algorithm to form the neighbor list        
            cpuFullNeighborList2D(nb.neighlist, nb.neighnum, x, app.rcutsq, anum, inum, jnum, dim);    
//         printArray2D(nb.neighnum, 1, inum, common.backend);  
//         printArray2D(nb.neighlist, jnum, inum, common.backend);  
        }
        else { // form neighbor list using cell list
            Int *a2clist = &tmp.intmem[0]; //anum
            Int *c2alist = &tmp.intmem[anum]; // anum
            Int *c2anum = &tmp.intmem[2*anum]; // anum
            Int *c2anumsum = &tmp.intmem[3*anum]; // cnum+1

            cpuCellList2D(a2clist, x, nb.eta1, nb.eta2, nb.eta3, nb.s2rmap, nb.cellnum, inum, anum, dim);                
            cpuCell2AtomList(c2alist, c2anumsum, c2anum, a2clist, anum, cnum);                              
            cpuFullNeighborList2D(nb.neighlist, nb.neighnum, x, app.rcutsq, nb.alist, a2clist, 
                    c2alist, c2anumsum, nb.cellnum, inum, jnum, dim);              
//         printArray2D(a2clist, 1, anum, common.backend);  
//         printArray2D(c2anum, 1, anum, common.backend);  
//         printArray2D(c2alist, 1, anum, common.backend);          
//         printArray2D(c2anumsum, 1, cnum+1, common.backend);      
//         printArray2D(nb.neighnum, 1, inum, common.backend);  
//         printArray2D(nb.neighlist, jnum, inum, common.backend);  
        }
    }
    else {
        cpuAtomList3D(nb.alist, inside, glistnumsum, glistnum,  x, nb.pimages, nb.rbvertices, 
            nb.s2rmap, inum, pnum, dim);
        
        gnum = IntArrayGetValueAtIndex(glistnumsum, inum, common.backend);        
        anum = inum + gnum;
        common.inum = inum;
        common.gnum = gnum;
        common.anum = anum;

        if (common.neighcell == 0) { // O(N^2) algorithm to form the neighbor list
            cpuFullNeighborList3D(nb.neighlist, nb.neighnum, x, app.rcutsq, anum, inum, jnum, dim);    
//             printArray2D(nb.neighnum, 1, inum, common.backend);  
//             printArray2D(nb.neighlist, jnum, inum, common.backend);                  
        }
        else { // form neighbor list using cell list
            Int *clist = &tmp.intmem[0]; //anum
            Int *c2alist = &tmp.intmem[anum]; // anum
            Int *c2anum = &tmp.intmem[2*anum]; // cnum
            Int *c2anumsum = &tmp.intmem[2*anum+cnum]; // cnum+1

            cpuCellList3D(clist, x, nb.eta1, nb.eta2, nb.eta3, nb.s2rmap, nb.cellnum, inum, anum, dim);                
            cpuCell2AtomList(c2alist, c2anumsum, c2anum, clist, anum, cnum);        
            cpuFullNeighborList3D(nb.neighlist, nb.neighnum, x, app.rcutsq, nb.alist, clist, 
                    c2alist, c2anumsum, nb.cellnum, inum, jnum, dim);           
//             printArray2D(nb.neighnum, 1, inum, common.backend);  
//             printArray2D(nb.neighlist, jnum, inum, common.backend);          
        }
    }    
    //cpuCumsum(nb.neighnumsum, nb.neighnum, inum+1);    
}

void implReadInputFiles(appstruct &app, configstruct &config, commonstruct &common, 
        string filein, string fileout, Int mpiprocs, Int mpirank, Int backend)
{
    implReadAppStruct(app, filein, mpiprocs, mpirank, backend);    
    implReadConfigStruct(config, filein, mpiprocs, mpirank, backend);            
    implSetCommonStruct(common, app, config, filein,  fileout, mpiprocs, mpirank, backend);        
}

void implSetConfiguration(neighborstruct &nb, tempstruct &tmp, sysstruct &sys, appstruct &app, configstruct &config, commonstruct &common, Int ci)
{
    implSetNeighborStruct(nb, common, config, ci);            
    implSetSysStruct(sys, common, config, ci);       
    implSetAtomBlocks(common);
    implSetTempStruct(tmp, common);  
}







//         if (common.cutofftype) {
//             if (common.neightype)
//                 cpuHalfNeighborList2D(nb.neighlist, nb.neighnum, x, app.cutoffellipsoid, nb.alist, clist, 
//                         c2alist, c2anumsum, nb.cellnum, inum, jnum, dim);                
//             else
//                 cpuFullNeighborList2D(nb.neighlist, nb.neighnum, x, app.cutoffellipsoid, nb.alist, clist, 
//                         c2alist, c2anumsum, nb.cellnum, inum, jnum, dim);                                                                            
//         }
//         else {
//             if (common.neightype)
//                 cpuHalfNeighborList2D(nb.neighlist, nb.neighnum, x, app.rcutsq, nb.alist, clist, 
//                         c2alist, c2anumsum, atomtype, nb.cellnum, ntype, inum, jnum, dim);                
//             else
//                 cpuFullNeighborList2D(nb.neighlist, nb.neighnum, x, app.rcutsq, nb.alist, clist, 
//                         c2alist, c2anumsum, atomtype, nb.cellnum, ntype, inum, jnum, dim);                            
//         }               

//         if (common.cutofftype) {
//             if (common.neightype)
//                 cpuHalfNeighborList3D(nb.neighlist, nb.neighnum, x, app.cutoffellipsoid, nb.alist, clist, 
//                         c2alist, c2anumsum, nb.cellnum, inum, jnum, dim);                
//             else
//                 cpuFullNeighborList3D(nb.neighlist, nb.neighnum, x, app.cutoffellipsoid, nb.alist, clist, 
//                         c2alist, c2anumsum, nb.cellnum, inum, jnum, dim);                                                                            
//         }
//         else {
//             if (common.neightype)
//                 cpuHalfNeighborList3D(nb.neighlist, nb.neighnum, x, app.rcutsq, nb.alist, clist, 
//                         c2alist, c2anumsum, atomtype, nb.cellnum, ntype, inum, jnum, dim);                
//             else
//                 cpuFullNeighborList3D(nb.neighlist, nb.neighnum, x, app.rcutsq, nb.alist, clist, 
//                         c2alist, c2anumsum, atomtype, nb.cellnum, ntype, inum, jnum, dim);                            
//         }        


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

