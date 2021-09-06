/***************************************************************************                               
                    Molecular Dynamics Potentials (MDP)
                           CESMIX-MIT Project  
 
 Contributing authors: Ngoc-Cuong Nguyen (cuongng@mit.edu, exapde@gmail.com)
 ***************************************************************************/

#ifndef __READDATA
#define __READDATA

void implReadAppStruct(appstruct &app, commonstruct &common, string filein, Int mpiprocs, Int mpirank, Int backend)
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
    common.traininglist = readiarrayfromdouble(in, app.nsize[50]);        
    common.validatelist = readiarrayfromdouble(in, app.nsize[51]);        
    common.trainingnum = app.nsize[50];
    common.validatenum = app.nsize[51];
    readarray(in, &common.nveparam, app.nsize[52]);    
    readarray(in, &common.nvtparam, app.nsize[53]);    
    readarray(in, &common.snaparam, app.nsize[54]);    
    readarray(in, &common.snapelemradius, app.nsize[55]);    
    readarray(in, &common.snapelemweight, app.nsize[56]);    
    readarray(in, &common.snapcoeff, app.nsize[57]);    
    readarray(in, &common.createvelocity, app.nsize[58]);    
    
    int n = 0;
    for (int i=13; i<=22; i++) 
        n += app.nsize[i];
    TemplateMalloc(&app.muep, n, 0); 
    
    for (int i=0; i < app.nsize[13]; i++) 
        app.muep[i] = app.mu1a[i];
    n = app.nsize[13];    
    for (int i=0; i < app.nsize[14]; i++) 
        app.muep[i+n] = app.mu1b[i];
    n += app.nsize[14];    
    for (int i=0; i < app.nsize[15]; i++) 
        app.muep[i+n] = app.mu2a[i];
    n += app.nsize[15];    
    for (int i=0; i < app.nsize[16]; i++) 
        app.muep[i+n] = app.mu2b[i];
    n += app.nsize[16];    
    for (int i=0; i < app.nsize[17]; i++) 
        app.muep[i+n] = app.mu2c[i];
    n += app.nsize[17];        
    for (int i=0; i < app.nsize[18]; i++) 
        app.muep[i+n] = app.mu3a[i];
    n += app.nsize[18];    
    for (int i=0; i < app.nsize[19]; i++) 
        app.muep[i+n] = app.mu3b[i];
    n += app.nsize[19];    
    for (int i=0; i < app.nsize[20]; i++) 
        app.muep[i+n] = app.mu3c[i];
    n += app.nsize[20];    
    for (int i=0; i < app.nsize[21]; i++) 
        app.muep[i+n] = app.mu4a[i];
    n += app.nsize[21];    
    for (int i=0; i < app.nsize[22]; i++) 
        app.muep[i+n] = app.mu4b[i];
    n += app.nsize[22];    
    
    // Close file:
    in.close();        
}

int implReadConfigStruct(configstruct &config, string filein, Int mpiprocs, Int mpirank, Int backend)
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
        return 0;
        //error("Unable to open file " + filename);
       
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
    readarray(in, &config.q, config.nsize[8]);   
    readarray(in, &config.v, config.nsize[9]);   
    readarray(in, &config.f, config.nsize[10]);      
    readarray(in, &config.we, config.nsize[11]);      
    readarray(in, &config.wf, config.nsize[12]);      
                
    TemplateMalloc(&config.natomssum, config.nsize[1]+1, 0); 
    cpuCumsum(config.natomssum, config.natoms, config.nsize[1]+1); 
        
    // Close file:
    in.close();            
    
    return 1;    
}

int implReadLatticeStruct(latticestruct &lat, string filein)
{
    string filename = filein + "lattice.bin";                    
    
    // Open file to read
    ifstream in(filename.c_str(), ios::in | ios::binary);
    
    if (!in) {
        return 0;
    }
        
    // allocate memory
    lat.allocatememory(1);
    
    /* Read data */            
    dstype *tmp; 
    int *lsize = readiarrayfromdouble(in, 1);
    int *nsize = readiarrayfromdouble(in, lsize[0]);
    readarray(in, &tmp, nsize[0]);  
    readarraynomalloc(in, &lat.origin, nsize[1]);  
    int *orientx = readiarrayfromdouble(in, nsize[2]);
    int *orienty = readiarrayfromdouble(in, nsize[3]);
    int *orientz = readiarrayfromdouble(in, nsize[4]);
    readarraynomalloc(in, &lat.spacing, nsize[5]);  
    readarraynomalloc(in, &lat.a1, nsize[6]);
    readarraynomalloc(in, &lat.a2, nsize[7]);
    readarraynomalloc(in, &lat.a3, nsize[8]);
    readarraynomalloc(in, &lat.atombasis, nsize[9]);
    int *type = readiarrayfromdouble(in, nsize[10]);
    
    for (int i=0; i<nsize[10]; i ++)
        lat.atomtype[i] = type[i];
    for (int i=0; i<3; i ++) {
        lat.orientx[i] = orientx[i];
        lat.orienty[i] = orienty[i];
        lat.orientz[i] = orientz[i];
    }
    lat.style = (int) tmp[0];
    lat.nbasis = (int) tmp[1];
    lat.spaceflag = (int) tmp[2];
    lat.scale = tmp[3];
    free(tmp); free(lsize); free(nsize); free(type);  
    free(orientx); free(orienty); free(orientz);  
    
    // Close file:
    in.close();        
    
    return 1;    
}

int implReadRegionStruct(regionstruct &reg, string filein)
{
    string filename = filein + "region.bin";                    
    
    // Open file to read
    ifstream in(filename.c_str(), ios::in | ios::binary);
    
    if (!in) {
        return 0;
    }
    
    // allocate memory
    reg.allocatememory(1);
        
    /* Read data */            
    dstype *tmp; 
    int *lsize = readiarrayfromdouble(in, 1);
    int *nsize = readiarrayfromdouble(in, lsize[0]);
    readarray(in, &tmp, nsize[0]);  
    readarraynomalloc(in, &reg.boxlo, nsize[1]);  
    readarraynomalloc(in, &reg.boxhi, nsize[2]);  
    readarraynomalloc(in, &reg.boxtilt, nsize[3]);  
    
    reg.triclinic = (int) tmp[0];
    free(tmp); free(lsize); free(nsize);   
    
    // Close file:
    in.close();        
    
    return 1;    
}

int implReadDomainStruct(domainstruct &dom, string filein)
{
    string filename = filein + "domain.bin";                    
    
    // Open file to read
    ifstream in(filename.c_str(), ios::in | ios::binary);
    
    dom.allocatememory(1);    
    for (int i = 0; i<3; i++) {
        dom.boxlo_lamda[i] = 0.0;
        dom.boxhi_lamda[i] = 1.0;                   
        dom.sublo_lamda[i] = 0.0;
        dom.subhi_lamda[i] = 1.0;                            
    }                                
    
    if (!in) {
        return 0;
    }
            
    /* Read data */            
    dstype *tmp; 
    int *lsize = readiarrayfromdouble(in, 1);
    int *nsize = readiarrayfromdouble(in, lsize[0]);
    readarray(in, &tmp, nsize[0]);  
    readarraynomalloc(in, &dom.boxlo, nsize[1]);  
    readarraynomalloc(in, &dom.boxhi, nsize[2]);  
    readarraynomalloc(in, &dom.boxtilt, nsize[3]);  
    
    dom.triclinic = (int) tmp[0];
    free(tmp); free(lsize); free(nsize);   
    
    // Close file:
    in.close();        
    
    return 1;    
}

void implSetCommonStruct(commonstruct &common, appstruct &app, configstruct &config, string filein, string fileout, Int mpiprocs, Int mpirank, Int backend)
{
    common.filein = filein;
    common.fileout = fileout;
    common.mpiProcs = mpiprocs;
    common.mpiRank = mpirank;
    common.backend = backend;
    
    common.dim = app.ndims[0];
    common.L = app.ndims[1];  // order of radial basis functions
    common.K = app.ndims[2];  // order of spherical harmonics          
    common.ntimesteps = app.ndims[3]; // number of time steps
    common.nab = app.ndims[4]; // number of atoms per block
    common.natomtypes = app.ndims[8]; // number of atom types
    common.nmoletypes = app.ndims[9]; // number of molecule types
    common.neighmax = app.ndims[10]; // % maximum number of neighbors allowed
    common.neighevery = app.ndims[11]; // rebuild neighborlist every iterations
    common.neighdelay = app.ndims[12]; 
    common.globalfreq = app.ndims[13]; 
    common.peratomfreq = app.ndims[14]; 
        
    common.descriptor = app.flags[0];   // descriptor flag: 0 -> Spherical Harmonics Bessel, 1 -> SNAP
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
    common.chemtype = app.flags[11];      // 0 -> single atom-type basis functions, 1 -> double atom-type basis functions 
    common.dftdata = app.flags[12];       // 0 -> no data, 1 -> energies only, 2 -> forces only, 3 -> energies and forces
    common.unitstyle = app.flags[13];       
    common.ensemblemode = app.flags[14];   
    common.neighcheck = app.flags[15]; 
            
    // set unit parameters
    setunits(common, common.unitstyle);
        
    common.rcutmax = app.solversparam[0]; 
    common.neighskin = app.solversparam[1];        
    common.time = app.simulaparam[0]; // current time    
    common.dt = app.simulaparam[1];    
    common.currentstep = (int) common.time/common.dt;
    common.rcutml = sqrt(app.rcutsqml[0]);
    
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
    common.nmu[0] = 0;
    common.nmu[1] = common.nmu1a;
    common.nmu[2] = common.nmu[1] + common.nmu1b;
    common.nmu[3] = common.nmu[2] + common.nmu2a;
    common.nmu[4] = common.nmu[3] + common.nmu2b;
    common.nmu[5] = common.nmu[4] + common.nmu2c;
    common.nmu[6] = common.nmu[5] + common.nmu3a;
    common.nmu[7] = common.nmu[6] + common.nmu3b;
    common.nmu[8] = common.nmu[7] + common.nmu3c;
    common.nmu[9] = common.nmu[8] + common.nmu4a;
    common.nmu[10] = common.nmu[9] + common.nmu4b;
    common.Nempot = common.npot1a + common.npot1b;
    common.Nempot += common.npot2a + common.npot2b + common.npot2c;
    common.Nempot += common.npot3a + common.npot3b + common.npot3c;
    common.Nempot += common.npot4a + common.npot4b;
    
    if (common.readconfig) {
        common.nconfigs = config.nsize[1]; // number of configurations
        common.ne = config.nsize[5];
        common.nt = config.nsize[6];
        common.nx = config.nsize[7];
        common.nq = config.nsize[8];
        common.nv = config.nsize[9];
        common.nf = config.nsize[10];
        common.ncq = common.nq/(common.nt);        

        common.inum = config.natoms[0];    
        common.nlocal = common.inum;    
        common.inummax = config.natoms[0];        
        for (int i=1; i<common.nconfigs; i++)         
            if (config.natoms[i] > common.inummax)
                common.inummax = config.natoms[i];   
    }
                    
    for (int i=0; i<common.dim; i++) {       
        common.pbc[i] = app.pbc[i];
        common.boxoffset[i] = app.boxoffset[i];
    }
    
    common.dtarray[0] = common.dt;     // dt
    common.dtarray[1] = 0.5*common.dt*common.ftm2v; // dtf
    common.dtarray[2] = common.dt;     // dtv
    common.dtarray[3] = 0;             // beginstep 
    common.dtarray[4] = common.ntimesteps; // endstep
    common.dtarray[5] = common.currentstep;      
    
    if (common.ensemblemode == 1) { // NVE limit ensemble                
      common.vlimitsq = (common.nveparam[0]/common.dt) * (common.nveparam[0]/common.dt);
    }
    
    if (common.ensemblemode == 2) { // NVT ensemble                
        common.tarray[0] = common.nvtparam[0]; // temp start
        common.tarray[1] = common.nvtparam[1]; // temp stop
        common.tarray[2] = 1/common.nvtparam[2]; // temp frequency    
        common.tarray[5] = common.tdof; 
        common.tarray[6] = common.boltz; 
        common.tarray[7] = (app.nsize[53]>3) ? common.nvtparam[3] : 0.0; // drag factor added to barostat/thermostat         
        common.tarray[8] = common.mvv2e; 
        common.mtchain   = (app.nsize[53]>4) ? (int) common.nvtparam[4] : 3;         
        common.nc_tchain = (app.nsize[53]>5) ? (int) common.nvtparam[5] : 1;                 
    }
    
    if (common.readconfig==0) {
        common.nconfigs = 1; // number of configurations
        common.ne = 0;
        common.nt = 1;
        common.nx = 1;
        common.nq = 0;
        common.nv = 0;
        common.nf = 0;
        common.ncq = 0;                
        common.inum = 0;    
        common.nlocal = 0;    
        common.inummax = 0;        
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
    
    cpuArrayCopy(common.pot1a, app.pot1a, app.nsize[23]);
    cpuArrayCopy(common.pot1b, app.pot1b, app.nsize[24]);
    cpuArrayCopy(common.pot2a, app.pot2a, app.nsize[25]);
    cpuArrayCopy(common.pot2b, app.pot2b, app.nsize[26]);
    cpuArrayCopy(common.pot2c, app.pot2c, app.nsize[27]);
    cpuArrayCopy(common.pot3a, app.pot3a, app.nsize[28]);
    cpuArrayCopy(common.pot3b, app.pot3b, app.nsize[29]);
    cpuArrayCopy(common.pot3c, app.pot3c, app.nsize[30]);
    cpuArrayCopy(common.pot4a, app.pot4a, app.nsize[31]);
    cpuArrayCopy(common.pot4b, app.pot4b, app.nsize[32]);    
    cpuArrayCopy(common.atom1b, app.atom1b, app.nsize[43]);
    cpuArrayCopy(common.atom2b, app.atom2b, app.nsize[44]);
    cpuArrayCopy(common.atom2c, app.atom2c, app.nsize[45]);
    cpuArrayCopy(common.atom3b, app.atom3b, app.nsize[46]);
    cpuArrayCopy(common.atom3c, app.atom3c, app.nsize[47]);
    cpuArrayCopy(common.atom4b, app.atom4b, app.nsize[48]);    
}

void implSetDomainStruct(domainstruct &dom, commonstruct &common, regionstruct &reg, latticestruct &lat)
{
    if (common.readregion == 0) 
        error("Domain and region are not defined.");
    
    if (common.readlattice==0) {
        common.dom.triclinic = common.reg.triclinic;
        for (int i = 0; i< 3; i++) {
            common.dom.boxhi[i] = common.reg.boxhi[i]; 
            common.dom.boxlo[i] = common.reg.boxlo[i]; 
            common.dom.boxtilt[i] = common.reg.boxtilt[i]; 
            common.dom.boxlo_lamda[i] = 0.0;
            common.dom.boxhi_lamda[i] = 1.0;                    
        }
        cpuSetGlobalBox(common.dom.h, common.dom.h_inv, common.dom.boxlo_bound, common.dom.boxhi_bound, 
                common.dom.boxhi, common.dom.boxlo, common.dom.boxtilt, common.dom.triclinic);                
        cpuSetLocalOrthBox(common.dom.subhi, common.dom.sublo, common.dom.boxhi, common.dom.boxlo, 
                common.dom.subhi_lamda, common.dom.sublo_lamda, 3);
                
        dstype epsilon[3]; 
        for (int i=0; i<3; i++) epsilon[i] = 1e-6*(common.dom.boxhi[i]-common.dom.boxlo[i]);
        cpuShiftedSubbox(common.dom.ssublo, common.dom.ssubhi, common.dom.boxlo, common.dom.boxhi, 
                common.dom.boxlo_lamda, common.dom.boxhi_lamda, common.dom.sublo, common.dom.subhi, 
                common.dom.sublo_lamda, common.dom.subhi_lamda, epsilon, common.pbc, common.dom.triclinic);
        cpuBoundingSubbox(common.dom.bsublo, common.dom.bsubhi, common.dom.sublo, common.dom.subhi,   
                common.dom.sublo_lamda, common.dom.subhi_lamda, common.dom.boxlo, common.dom.h, common.dom.triclinic);                                                
    }
    else {
        cpuLattice(common.lat.atombasis, common.lat.primitive, common.lat.rotaterow, common.lat.primitinv, 
                common.lat.rotatecol, common.lat.origin, common.lat.spacing, common.lat.a1, common.lat.a2, 
                common.lat.a3, common.lat.scale, common.lat.orientx, common.lat.orienty, common.lat.orientz, 
                common.lat.style, common.unitstyle, common.lat.spaceflag, common.dim);

        common.dom.triclinic = common.reg.triclinic;
        common.dom.boxtilt[0] = common.lat.spacing[0]*common.reg.boxtilt[0];                     
        common.dom.boxtilt[1] = common.lat.spacing[0]*common.reg.boxtilt[1];                     
        common.dom.boxtilt[2] = common.lat.spacing[1]*common.reg.boxtilt[2];                                     
        for (int i = 0; i< 3; i++) {
            common.dom.boxhi[i] = common.lat.spacing[i]*common.reg.boxhi[i]; 
            common.dom.boxlo[i] = common.lat.spacing[i]*common.reg.boxlo[i];        
            common.dom.boxlo_lamda[i] = 0.0;
            common.dom.boxhi_lamda[i] = 1.0;                    
        }                          
        
        cpuSetGlobalBox(common.dom.h, common.dom.h_inv, common.dom.boxlo_bound, common.dom.boxhi_bound, 
                common.dom.boxhi, common.dom.boxlo, common.dom.boxtilt, common.dom.triclinic);                
        cpuSetLocalOrthBox(common.dom.subhi, common.dom.sublo, common.dom.boxhi, common.dom.boxlo, 
                common.dom.subhi_lamda, common.dom.sublo_lamda, 3);
                
        dstype epsilon[3]; 
        for (int i=0; i<3; i++) epsilon[i] = 1e-6*(common.dom.boxhi[i]-common.dom.boxlo[i]);
        cpuShiftedSubbox(common.dom.ssublo, common.dom.ssubhi, common.dom.boxlo, common.dom.boxhi, 
                common.dom.boxlo_lamda, common.dom.boxhi_lamda, common.dom.sublo, common.dom.subhi, 
                common.dom.sublo_lamda, common.dom.subhi_lamda, epsilon, common.pbc, common.dom.triclinic);
        cpuBoundingSubbox(common.dom.bsublo, common.dom.bsubhi, common.dom.sublo, common.dom.subhi,   
                common.dom.sublo_lamda, common.dom.subhi_lamda, common.dom.boxlo, common.dom.h, common.dom.triclinic);                                        
        cpuLatticeBoundingBox(common.lat.sublo, common.lat.subhi, common.dom.bsublo, common.dom.bsubhi, 
                common.lat.primitive, common.lat.rotaterow, common.lat.primitinv, common.lat.rotatecol, 
                common.lat.origin, common.lat.spacing, common.lat.scale);
        common.lat.latticebounds();        
    }     
}

void implSetConfigStruct(configstruct &config, domainstruct &dom, appstruct &app, commonstruct &common, regionstruct &reg, latticestruct &lat)
{
    int dim = common.dim;
    int mpiRank = common.mpiRank;
    int backend = 1;    
    int *ilist;
    dstype *y, *amass;        
        
    TemplateMalloc(&y, common.dim*common.lat.natom, backend);
    TemplateMalloc(&ilist, common.lat.natom, backend);           
    TemplateMalloc(&config.x, common.dim*common.lat.natom, backend);
    TemplateMalloc(&config.t, common.lat.natom, backend);                
    TemplateMalloc(&config.natoms, 1, backend);     
    TemplateMalloc(&config.natomssum, 2, backend);         
    TemplateMalloc(&config.a, 3, backend);     
    TemplateMalloc(&config.b, 3, backend);     
    TemplateMalloc(&config.c, 3, backend);     
    
    config.a[0] = dom.boxhi[0]-dom.boxlo[0]; config.a[1] = 0.0;                       config.a[2] = 0.0;
    config.b[0] = dom.boxtilt[0];            config.b[1] = dom.boxhi[1]-dom.boxlo[1]; config.b[2] = 0.0;
    config.c[0] = dom.boxtilt[1];            config.c[1] = dom.boxtilt[2];            config.c[2] = dom.boxhi[2] - dom.boxlo[2];
    
    cpuAtomLattice(y, ilist, common.lat.atombasis, common.lat.primitive, common.lat.rotaterow, common.lat.origin, 
            common.lat.spacing, common.lat.scale, common.lat.atomtype, common.lat.nbasis, 0, common.lat.ilo, 
            common.lat.ihi, common.lat.jlo, common.lat.jhi, common.lat.klo, common.lat.khi, dim);
        
    int nlocal = cpuAtomAdd(config.x, config.t, y, common.dom.h_inv, common.dom.boxlo, common.dom.sublo_lamda,
        common.dom.subhi_lamda, ilist, common.dim, 0, common.lat.natom);

    config.natoms[0] = nlocal;
    config.natomssum[0] = 0;
    config.natomssum[1] = nlocal;
    common.nlocal = nlocal;
    common.inum = nlocal;
    common.inummax = nlocal;
    common.tdof = (nlocal-1)*dim;    
    common.tfactor = common.mvv2e / (common.tdof * common.boltz);    
            
    TemplateMalloc(&config.v, dim*nlocal, backend);        
    TemplateMalloc(&amass, nlocal, backend);        
    for (int i=0; i<nlocal; i++) ilist[i] = i;
        
    if (app.nsize[58]>0) {
        dstype t_desired = common.createvelocity[0];
        int seed0 = (int) common.createvelocity[1];
        int dist_flag = (int) common.createvelocity[2];
        int sum_flag = (int) common.createvelocity[3];
        int loop_flag = (int) common.createvelocity[4];
        int momentum_flag = (int) common.createvelocity[5];         
        int rotation_flag = (int) common.createvelocity[6];     

        int natom=nlocal;
        cpuVelocityCreate(config.x, config.v, app.atommass, common.second, common.seed, common.save, ilist, 
                config.t, seed0, sum_flag, dist_flag, loop_flag, dim, mpiRank, nlocal, natom);
        
        common.masstotal = cpuComputeMass(amass, app.atommass, config.t, ilist, nlocal);                                
        cpuComputeVCM(common.vcm, config.v, app.atommass, common.masstotal, ilist, config.t, dim, nlocal);                
        cpuVelocityZeroMomentum(config.v, common.vcm, dim, nlocal);

        common.temp = cpuComputeTempScalar(config.v, app.atommass, common.tfactor, config.t, ilist, dim, nlocal);    
        cpuArrayMultiplyScalar(config.v, sqrt(t_desired/common.temp), nlocal*dim);                                              
        common.nv = 1;
    }
    
    free(y); free(ilist); free(amass);           
}

  
void implReadInputFiles(appstruct &app, configstruct &config, commonstruct &common, 
        string filein, string fileout, Int mpiprocs, Int mpirank, Int backend)
{        
    implReadAppStruct(app, common, filein, mpiprocs, mpirank, backend);        
    common.readconfig = implReadConfigStruct(config, filein, mpiprocs, mpirank, backend);    
    common.readlattice = implReadLatticeStruct(common.lat, filein);    
    common.readregion = implReadRegionStruct(common.reg, filein);    
    common.readdomain = implReadDomainStruct(common.dom, filein);   
    implSetCommonStruct(common, app, config, filein,  fileout, mpiprocs, mpirank, backend);      
    if (common.readdomain==0) {
        implSetDomainStruct(common.dom, common, common.reg, common.lat);        
    }    
    if (common.readconfig==0) {
        implSetConfigStruct(config, common.dom, app, common, common.reg, common.lat);            
    }
    common.triclinic = common.dom.triclinic;
    app.dom.allocatememory(backend);
    copydomain(app.dom, common.dom, backend);       
    
#ifdef HAVE_MPI        
    MPI_Allreduce(&common.nlocal, &common.natoms, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
#else    
    common.natoms = common.nlocal;
#endif    
    common.tdof = (common.natoms-1)*common.dim;        
}

#endif

