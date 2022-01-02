/***************************************************************************                               
                    Molecular Dynamics Potentials (MDP)
                           CESMIX-MIT Project  
 
 Contributing authors: Ngoc-Cuong Nguyen (cuongng@mit.edu, exapde@gmail.com)
 ***************************************************************************/

#ifndef __INTEGRATION
#define __INTEGRATION

#include "Integration.h"

CIntegration::CIntegration(CCalculation &CCal)
{
}

void CIntegration::VelocityVerlet(CCalculation &CCal)
{    
  INIT_TIMING;
  
  // store atom positions at the beginning of the time loop  
  dstype *xstart = CCal.tmp.tmpmem;
  
  this->IntegrationSetup(CCal);  
  
  auto starttime = std::chrono::high_resolution_clock::now();   
  
  int ntimesteps = CCal.common.ntimesteps; // # timesteps in this simulation 
  for (int istep = 0; istep < ntimesteps; istep++) {    
    CCal.common.currentstep += 1; // current time step    
    CCal.common.time += CCal.common.dt; // current simulation time          
    //printf("\nTimestep :  %d,   Time : %g\n",CCal.common.currentstep, CCal.common.time);        
    
    int dim = CCal.common.dim;
    int inum = CCal.common.inum;    
    int backend = CCal.common.backend;
    
    // copy atoms I to xstart  
    ArrayCopy(xstart, CCal.sys.x, dim*inum, backend);
    
//     if (backend<=1) {
//         writearray2file("x1a.bin", CCal.sys.x, dim*inum, backend);    
//         writearray2file("v1a.bin", CCal.sys.v, dim*inum, backend);    
//         writearray2file("f1a.bin", CCal.sys.f, dim*inum, backend);    
//     }
//     else {
//         writearray2file("x2a.bin", CCal.sys.x, dim*inum, backend);    
//         writearray2file("v2a.bin", CCal.sys.v, dim*inum, backend);    
//         writearray2file("f2a.bin", CCal.sys.f, dim*inum, backend);    
//     }
    
    // initial integration of the velocity vertlet method: update velocity and position        
    START_TIMING;
    this->InitialIntegration(CCal);
    END_TIMING_CCAL(0);    
    
    // impose position/velocity constraints and periodic boundary conditions
    START_TIMING;
    this->PostInitialIntegration(CCal);
    END_TIMING_CCAL(1);    
    
    // redistribute atoms among processors to balance computation load
    START_TIMING;
    this->DynamicLoadBalancing(CCal);
    END_TIMING_CCAL(2);    
    
    // exchange data between processors 
    START_TIMING;
    this->ProcessorsCommunication(CCal);
    END_TIMING_CCAL(3);    
    
    // rebuild neighbor list if necessary
    START_TIMING;
    this->NeighborListRebuild(CCal, xstart, istep);
    END_TIMING_CCAL(4);    
    
//     if (backend<=1) {
//         writearray2file("x1b.bin", CCal.sys.x, dim*inum, backend);    
//         writearray2file("v1b.bin", CCal.sys.v, dim*inum, backend);    
//     }
//     else {
//         writearray2file("x2b.bin", CCal.sys.x, dim*inum, backend);    
//         writearray2file("v2b.bin", CCal.sys.v, dim*inum, backend);    
//     }
    
    // calculate atomic forces
    START_TIMING;
    CCal.PotentialEnergyForceVirial(CCal.sys.e, CCal.sys.f, CCal.sys.vatom, CCal.sys.x, CCal.app.muep, CCal.common.nmu);                         
    END_TIMING_CCAL(5);    
    
//     if (backend<=1) {
//         writearray2file("f1.bin", CCal.sys.f, dim*inum, backend);    
//     }
//     else {
//         writearray2file("f2.bin", CCal.sys.f, dim*inum, backend);    
//     }
    
    // // impose force constraints if any
    START_TIMING;
    this->PostForceComputation(CCal);
    END_TIMING_CCAL(6);    
        
//     if (backend<=1) {
//         writearray2file("x1c.bin", CCal.sys.x, dim*inum, backend);    
//         writearray2file("v1c.bin", CCal.sys.v, dim*inum, backend);    
//         writearray2file("f1c.bin", CCal.sys.f, dim*inum, backend);    
//     }
//     else {
//         writearray2file("x2c.bin", CCal.sys.x, dim*inum, backend);    
//         writearray2file("v2c.bin", CCal.sys.v, dim*inum, backend);    
//         writearray2file("f2c.bin", CCal.sys.f, dim*inum, backend);    
//     }
    
    // final integration of the velocity vertlet method: update velocity after computing force    
    START_TIMING;
    this->FinalIntegration(CCal);    
    END_TIMING_CCAL(7);    
                
//     if (backend<=1) {
//         writearray2file("x1d.bin", CCal.sys.x, dim*inum, backend);    
//         writearray2file("v1d.bin", CCal.sys.v, dim*inum, backend);    
//         writearray2file("f1d.bin", CCal.sys.f, dim*inum, backend);    
//     }
//     else {
//         writearray2file("x2d.bin", CCal.sys.x, dim*inum, backend);    
//         writearray2file("v2d.bin", CCal.sys.v, dim*inum, backend);    
//         writearray2file("f2d.bin", CCal.sys.f, dim*inum, backend);    
//     }
    
    // fix velocity and rescale velocity if necessary to control the temperature
    START_TIMING;
    this->PostFinalIntegration(CCal);        
    END_TIMING_CCAL(8);    
    
    // compute output, print on screen, and save it into binary files
    START_TIMING;
    this->TimestepCompletion(CCal);                   
    END_TIMING_CCAL(9);            
  }    
  
  //CUDA_SYNC;
  auto endtime = std::chrono::high_resolution_clock::now();
  dstype time = std::chrono::duration_cast<std::chrono::nanoseconds>(endtime-starttime).count()/1e9;        
  if (CCal.common.mpiRank==0) {
    printf("\n\tTotal number of atoms = %i\n", CCal.common.natoms);
    printf("\tTotal simulation time = %g seconds\n", time);
    printf("\tSimulation time per timestep = %g seconds\n", time/ntimesteps);
  }  
  
#ifdef _TIMING  
  if (CCal.common.mpiRank==0) {
	printf("Timing MDP for peformance improvement: \n");
    printArray2D(CCal.common.timing, 10, 10, 1);    
  }
#endif        
}

void CIntegration::IntegrationSetup(CCalculation &CCal)
{
    writearray2file("x1.bin", CCal.sys.x, CCal.common.dim*CCal.common.inum, CCal.common.backend);    
    writearray2file("v1.bin", CCal.sys.v, CCal.common.dim*CCal.common.inum, CCal.common.backend);    
    
    ArraySetValue(CCal.common.timing, 0.0, 100, 1);
    CCal.NeighborList(CCal.sys.x);
    //CCal.nb.printout(CCal.common.inum, CCal.common.gnum, CCal.common.neighmax, CCal.common.backend);    
    CCal.PotentialEnergyForceVirial(CCal.sys.e, CCal.sys.f, CCal.sys.vatom, CCal.sys.x, CCal.app.muep, CCal.common.nmu);                         
    //CCal.sys.printout(CCal.common.dim, CCal.common.inum, CCal.common.backend);
    CCal.ThermoOutput(0);         
    if (CCal.common.mpiRank==0) {
        printf("Total number of atoms = %i\n", CCal.common.natoms);
        printf("Time step     Temperature     Potential energy     Total energy     Pressure \n"); 
        printf("\t%i  \t %g   \t\t  %g   \t     %g   \t  %g\n", 0, CCal.common.temp, CCal.common.pe, CCal.common.ke+CCal.common.pe, CCal.common.pres);                    
    }
    //error("here");
}

void CIntegration::InitialIntegration(CCalculation &CCal)
{
    dstype *x = CCal.sys.x;
    dstype *e = CCal.sys.e;
    dstype *f = CCal.sys.f;
    dstype *v = CCal.sys.v;
    dstype *mass = CCal.app.atommass;        
    dstype *dtarray = CCal.common.dtarray;
    dstype *tarray = CCal.common.tarray;
    dstype vlimitsq = CCal.common.vlimitsq;
        
    int backend = CCal.common.backend;
    int dim = CCal.common.dim;
    int ensemblemode = CCal.common.ensemblemode; 
    int biasflag = CCal.common.biasflag;
    int eta_mass_flag = CCal.common.eta_mass_flag;
    int mtchain = CCal.common.mtchain;
    int nc_tchain = CCal.common.nc_tchain;
    int *atomtype = CCal.nb.atomtype;    
    
    //int ensemblegroup = CCal.common.ensemblegroup;     
    //int inum = CCal.common.inumgroup[ensemblegroup];        
    //int *ilist = CCal.nb.atomgroups[ensemblegroup];
    int inum = CCal.common.inum;
    int *ilist = CCal.nb.alist;    
    
//     printf("%i \n", ensemblemode);
//     printArray2D(dtarray, 1, 6, backend);
//     printArray2D(tarray, 1, 9, backend);
    
    InitialIntegrate(x, v, f, mass, dtarray, tarray, CCal.common.eta_mass, CCal.common.eta, 
           CCal.common.eta_dot, CCal.common.eta_dotdot, &CCal.tmp.tmpmem[dim*inum], &CCal.tmp.tmpmem[dim*inum+inum], 
            vlimitsq, atomtype, ilist, eta_mass_flag, biasflag, mtchain, nc_tchain, ensemblemode, 
            dim, inum, backend);        
}

void CIntegration::PostInitialIntegration(CCalculation &CCal)
{
    dstype *x = CCal.sys.x;
    dstype *f = CCal.sys.f;
    dstype *v = CCal.sys.v;
    dstype *eatom = CCal.sys.eatom;
    dstype *vatom = CCal.sys.vatom;    
    int backend = CCal.common.backend;    
    int dim = CCal.common.dim;    
    int eflag_atom = 0;
    int vflag_atom = 0;            
    
    int nsetvelocity = CCal.common.nsetvelocity;    
    for (int i=0; i<nsetvelocity; i++) {       
        int g = CCal.common.gsetvelocity[i];
        int inum = CCal.common.inumgroup[g];        
        int *ilist = CCal.nb.atomgroups[g];                    
        int *iparam = &CCal.app.isetvelocity[i];
        dstype *fparam = &CCal.app.fsetvelocity[i];       
        dstype dtf = CCal.common.dtarray[1];
        dstype dtv = CCal.common.dtarray[2];            
        SetVelocityInitialIntegrate(x, v, f, CCal.app.atommass, fparam, dtf, 
                dtv, CCal.nb.atomtype, ilist, iparam, dim, inum, backend);        
    }            
    
    int nwallreflect = CCal.common.nwallreflect;    
    for (int i=0; i<nwallreflect ; i++) {       
        int g = CCal.common.gwallreflect[i];
        int inum = CCal.common.inumgroup[g];        
        int *ilist = CCal.nb.atomgroups[g];                    
        int *iparam = &CCal.app.iwallreflect[i];
        dstype *fparam = &CCal.app.fwallreflect[i];        
        FixWallReflect(x, v, f, eatom, vatom, fparam, 
            iparam, ilist, eflag_atom, vflag_atom, dim, inum, backend);    
    }        
}

void CIntegration::DynamicLoadBalancing(CCalculation &CCal)
{
    
}

void CIntegration::ProcessorsCommunication(CCalculation &CCal)
{
    
}

void CIntegration::NeighborListRebuild(CCalculation &CCal, dstype *xstart, int istep)
{    
    int build = 0;    
    int ago = istep + 1;
    int delay = CCal.common.neighdelay;
    int every = CCal.common.neighevery;
    int distcheck = CCal.common.neighcheck;
    int dim = CCal.common.dim;
    int backend = CCal.common.backend;
    int inum = CCal.common.inum;
            
    if ((ago >= delay) && (ago % every == 0)) {
        if (distcheck == 0) {
            build = 1;
        } else {
            dstype *y = &CCal.tmp.tmpmem[dim*inum];
            dstype *z = &CCal.tmp.tmpmem[dim*inum+inum];
            dstype skin = CCal.common.neighskin;
            
            //Unmap(z, CCal.sys.x, CCal.app.dom.h, CCal.sys.image, CCal.common.triclinic, dim, inum, backend);
            ArrayDistSquareSum(y, CCal.sys.x, CCal.sys.xhold, dim, inum, backend);
            dstype maxdistsquare = ArrayMax(y, z, inum, backend);
            //printf("Maximum Distance Square %g...\n",maxdistsquare);
            if (maxdistsquare > 0.25*skin*skin)
                build = 1;            
        }        
    } 

    if (build == 1) {
        if (CCal.common.mpiRank==0)
            printf("Rebuild neighbor list at timestep %i...\n",istep);
        
        PBC(CCal.sys.x, CCal.sys.v, CCal.sys.image, CCal.app.dom.boxhi, CCal.app.dom.boxlo, CCal.app.dom.boxhi_lamda, 
                CCal.app.dom.boxlo_lamda, CCal.app.dom.h, CCal.app.dom.h_inv, CCal.app.dom.h_rate, 
                CCal.app.pbc, CCal.common.vdeform, CCal.common.triclinic, dim, inum, backend);  
        CCal.NeighborList(CCal.sys.x);       
    } else {
        int gnum = CCal.common.gnum;
        // calculate the distance between x and xstart: xstart = x - xstart
        ArrayAXPBY(xstart, CCal.sys.x, xstart, one, minusone, dim*inum, backend);
        // shift periodic images of atoms I by the above distance
        ArrayPlusAtColumnIndex(&CCal.sys.x[dim*inum], xstart, &CCal.nb.alist[inum], dim, gnum, backend);                
    }    
}

void CIntegration::PostForceComputation(CCalculation &CCal)
{    
    // apply external forces on particular groups of atoms to 
    // impose wall/boundary conditions    
    dstype *x = CCal.sys.x;
    dstype *f = CCal.sys.f;
    dstype *v = CCal.sys.v;
    dstype *eatom = CCal.sys.eatom;
    dstype *vatom = CCal.sys.vatom;    
    dstype *box = CCal.app.dom.h;
    int *pbc = CCal.app.pbc;
    int *image = CCal.sys.image;    
    
    int dim = CCal.common.dim;    
    int backend = CCal.common.backend;    
    int triclinic = CCal.common.triclinic;
    int eflag_atom = 0;
    int vflag_atom = 0;            
    
    int nsetforce = CCal.common.nsetforce;    
    for (int i=0; i<nsetforce ; i++) {       
        int g = CCal.common.gsetforce[i];
        int inum = CCal.common.inumgroup[g];        
        int *ilist = CCal.nb.atomgroups[g];                    
        int *iparam = &CCal.app.isetforce[3*i];
        dstype *fparam = &CCal.app.fsetforce[3*i];
        FixSetForce(x, v, f, eatom, vatom, fparam, 
            iparam, ilist, eflag_atom, vflag_atom, dim, inum, backend);    
    }

    int nlineforce = CCal.common.nlineforce;    
    for (int i=0; i<nlineforce ; i++) {       
        int g = CCal.common.glineforce[i];
        int inum = CCal.common.inumgroup[g];        
        int *ilist = CCal.nb.atomgroups[g];                    
        int *iparam = CCal.app.ilineforce;
        dstype *fparam = &CCal.app.flineforce[3*i];
        FixLineForce(x, v, f, eatom, vatom, fparam, 
            iparam, ilist, eflag_atom, vflag_atom, dim, inum, backend);    
    }
    
    int nplaneforce = CCal.common.nplaneforce;    
    for (int i=0; i<nplaneforce ; i++) {       
        int g = CCal.common.gplaneforce[i];
        int inum = CCal.common.inumgroup[g];        
        int *ilist = CCal.nb.atomgroups[g];                    
        int *iparam = CCal.app.iplaneforce;
        dstype *fparam = &CCal.app.fplaneforce[3*i];
        FixPlaneForce(x, v, f, eatom, vatom, fparam, 
            iparam, ilist, eflag_atom, vflag_atom, dim, inum, backend);    
    }
    
//     int naveforce = CCal.common.naveforce;    
//     for (int i=0; i<naveforce ; i++) {       
//         int g = CCal.common.gaveforce[i];
//         int inum = CCal.common.inumgroup[g];        
//         int *ilist = CCal.nb.atomgroups[g];                    
//         int *iparam = &CCal.app.iaveforce[3*i];
//         dstype *fparam = &CCal.app.faveforce[3*i];
//         cpuFixAveForce(x, v, f, eatom, vatom, fparam, 
//             iparam, ilist, eflag_atom, vflag_atom, dim, inum);    
//     }

    int naddforce = CCal.common.naddforce;    
    for (int i=0; i<naddforce ; i++) {       
        int eflag_atom = CCal.common.iaddforce[0];
        int vflag_atom = CCal.common.iaddforce[1];                    
        int g = CCal.common.gaddforce[i];
        int inum = CCal.common.inumgroup[g];        
        int *ilist = CCal.nb.atomgroups[g];                    
        int *iparam = CCal.app.iaddforce;
        dstype *fparam = &CCal.app.faddforce[3*i];
        FixAddForce(x, v, f, eatom, vatom, fparam, box,
            iparam, ilist, image, triclinic, eflag_atom, vflag_atom, dim, inum, backend);    
    }

    int ndragforce = CCal.common.ndragforce;    
    for (int i=0; i<ndragforce ; i++) {       
        int eflag_atom = CCal.common.idragforce[0];
        int vflag_atom = CCal.common.idragforce[1];                    
        int g = CCal.common.gdragforce[i];
        int inum = CCal.common.inumgroup[g];        
        int *ilist = CCal.nb.atomgroups[g];                    
        int *iparam = &CCal.app.idragforce[3*i];
        dstype *fparam = &CCal.app.fdragforce[3*i];
        FixDragForce(x, v, f, eatom, vatom, fparam, box,
            iparam, ilist, pbc, triclinic, eflag_atom, vflag_atom, dim, inum, backend);    
    }
    
    int nwallharmonic = CCal.common.nwallharmonic;    
    for (int i=0; i<nwallharmonic ; i++) {       
        int eflag_atom = CCal.common.iwallharmonic[0];
        int vflag_atom = CCal.common.iwallharmonic[1];                    
        int g = CCal.common.gwallharmonic[i];
        int inum = CCal.common.inumgroup[g];        
        int *ilist = CCal.nb.atomgroups[g];                    
        int *iparam = &CCal.app.iwallharmonic[i];
        dstype *fparam = &CCal.app.fwallharmonic[3*i];
        FixWallHarmonic(x, v, f, eatom, vatom, fparam, 
            iparam, ilist, eflag_atom, vflag_atom, dim, inum, backend);    
    }

    int nwalllj93 = CCal.common.nwalllj93;    
    for (int i=0; i<nwalllj93 ; i++) {       
        int eflag_atom = CCal.common.iwalllj93[0];
        int vflag_atom = CCal.common.iwalllj93[1];                    
        int g = CCal.common.gwalllj93[i];
        int inum = CCal.common.inumgroup[g];        
        int *ilist = CCal.nb.atomgroups[g];                    
        int *iparam = &CCal.app.iwalllj93[i];
        dstype *fparam = &CCal.app.fwalllj93[7*i];
        FixWallLJ93(x, v, f, eatom, vatom, fparam, 
            iparam, ilist, eflag_atom, vflag_atom, dim, inum, backend);    
    }

    int nwalllj126 = CCal.common.nwalllj126;    
    for (int i=0; i<nwalllj126 ; i++) {       
        int eflag_atom = CCal.common.iwalllj126[0];
        int vflag_atom = CCal.common.iwalllj126[1];                    
        int g = CCal.common.gwalllj126[i];
        int inum = CCal.common.inumgroup[g];        
        int *ilist = CCal.nb.atomgroups[g];                    
        int *iparam = &CCal.app.iwalllj126[i];
        dstype *fparam = &CCal.app.fwalllj126[7*i];
        FixWallLJ126(x, v, f, eatom, vatom, fparam, 
            iparam, ilist, eflag_atom, vflag_atom, dim, inum, backend);    
    }

    int nwalllj1043 = CCal.common.nwalllj1043;    
    for (int i=0; i<nwalllj1043 ; i++) {       
        int eflag_atom = CCal.common.iwalllj1043[0];
        int vflag_atom = CCal.common.iwalllj1043[1];                    
        int g = CCal.common.gwalllj1043[i];
        int inum = CCal.common.inumgroup[g];        
        int *ilist = CCal.nb.atomgroups[g];                    
        int *iparam = &CCal.app.iwalllj1043[i];
        dstype *fparam = &CCal.app.fwalllj1043[10*i];
        FixWallLJ1043(x, v, f, eatom, vatom, fparam, 
            iparam, ilist, eflag_atom, vflag_atom, dim, inum, backend);    
    }
    
    int nwallmorse = CCal.common.nwallmorse;    
    for (int i=0; i<nwallmorse ; i++) {       
        int eflag_atom = CCal.common.iwallmorse[0];
        int vflag_atom = CCal.common.iwallmorse[1];                    
        int g = CCal.common.gwallmorse[i];
        int inum = CCal.common.inumgroup[g];        
        int *ilist = CCal.nb.atomgroups[g];                    
        int *iparam = &CCal.app.iwallmorse[i];
        dstype *fparam = &CCal.app.fwallmorse[6*i];
        FixWallMorse(x, v, f, eatom, vatom, fparam, 
            iparam, ilist, eflag_atom, vflag_atom, dim, inum, backend);    
    }    
}

void CIntegration::FinalIntegration(CCalculation &CCal)
{
    dstype *x = CCal.sys.x;
    dstype *e = CCal.sys.e;
    dstype *f = CCal.sys.f;
    dstype *v = CCal.sys.v;
    dstype *mass = CCal.app.atommass;        
    dstype *dtarray = CCal.common.dtarray;
    dstype *tarray = CCal.common.tarray;
    dstype vlimitsq = CCal.common.vlimitsq;
        
    int backend = CCal.common.backend;
    int dim = CCal.common.dim;
    int ensemblemode = CCal.common.ensemblemode;     
    int biasflag = CCal.common.biasflag;
    int eta_mass_flag = CCal.common.eta_mass_flag;
    int mtchain = CCal.common.mtchain;
    int nc_tchain = CCal.common.nc_tchain;
    int *atomtype = CCal.nb.atomtype;                
    
//     int ensemblegroup = CCal.common.ensemblegroup;     
//     int inum = CCal.common.inumgroup[ensemblegroup];        
//     int *ilist = CCal.nb.atomgroups[ensemblegroup];
    int inum = CCal.common.inum;
    int *ilist = CCal.nb.alist;    

//     printArray2D(v, dim, inum, backend);
//     printArray2D(x, dim, inum, backend);    
//     printArray2D(f, dim, inum, backend);    
    
    FinalIntegrate(x, v, f, mass, dtarray, tarray, CCal.common.eta_mass, CCal.common.eta, 
            CCal.common.eta_dot, CCal.common.eta_dotdot, &CCal.tmp.tmpmem[dim*inum], &CCal.tmp.tmpmem[dim*inum+inum], 
            vlimitsq, atomtype, ilist, eta_mass_flag, biasflag, mtchain, nc_tchain, ensemblemode, dim, inum, backend);               
}

void CIntegration::PostFinalIntegration(CCalculation &CCal)
{   
    dstype *x = CCal.sys.x;
    dstype *v = CCal.sys.v;
    dstype *f = CCal.sys.f;
    dstype *mass = CCal.app.atommass;       
    int *atomtype = CCal.nb.atomtype;      
    int backend = CCal.common.backend;
    int dim = CCal.common.dim;
        
    int nsetvelocity = CCal.common.nsetvelocity;    
    for (int i=0; i<nsetvelocity; i++) {       
        int g = CCal.common.gsetvelocity[i];
        int inum = CCal.common.inumgroup[g];        
        int *ilist = CCal.nb.atomgroups[g];                    
        int *iparam = &CCal.app.isetvelocity[i];               
        dstype dtf = CCal.common.dtarray[1];                
        SetVelocityFinalIntegrate(x, v, f, mass, dtf, 
                atomtype, ilist, iparam, dim, inum, backend);              
    }            
    
    int mode = CCal.common.vrtmode;
    dstype energy = 0.0;
    if (mode >= 0) {
        dstype *second = CCal.common.second;
        int *seed = CCal.common.seed;
        int *save = CCal.common.save;    
        int ensemblegroup = CCal.common.ensemblegroup;     
        int inum = CCal.common.inumgroup[ensemblegroup];        
        int biasflag = CCal.common.biasflag;          
        int *ilist = CCal.nb.atomgroups[ensemblegroup];            
        dstype *dtarray = CCal.common.dtarray;
        dstype *tarray = CCal.common.tarray;
        
        energy = VelocityRescalingThermostat(v, mass, dtarray, tarray, CCal.tmp.tmpmem, &CCal.tmp.tmpmem[inum],
                second, energy, atomtype, ilist, seed, save, biasflag, mode, dim, inum, backend);            
    }        
    CCal.common.vrtenergy = energy;
}

void CIntegration::TimestepCompletion(CCalculation &CCal)
{           
    if (CCal.common.currentstep % CCal.common.globalfreq == 0) {
        CCal.ThermoOutput(0);            
        if (CCal.common.mpiRank==0)
            printf("\t%i  \t %g   \t  %g   \t   %g   \t    %g\n", CCal.common.currentstep, CCal.common.temp, CCal.common.pe, CCal.common.ke+CCal.common.pe, CCal.common.pres);                        
    }
}
        
#endif        




