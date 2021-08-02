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
//     int M = CCal.common.M;
//     int backend = CCal.common.backend;    
//     TemplateMalloc(&CCal.sys.A, M*M, backend);   
}

void CIntegration::VelocityVerlet(CCalculation &CCal)
{    
    
  this->IntegrationSetup(CCal);
                  
  int ntimesteps = CCal.common.ntimesteps; // # timesteps in this simulation 
  for (int istep = 0; istep < ntimesteps; istep++) {    
    CCal.common.currentstep += 1; // current time step    
    CCal.common.time += CCal.common.dt; // current simulation time          
    printf("\nTimestep :  %d,   Time : %g\n",CCal.common.currentstep, CCal.common.time);        
    
    // initial integration of the velocity vertlet method: update velocity and position    
    this->InitialIntegration(CCal);
        
    // impose position/velocity constraints and periodic boundary conditions
    this->PostInitialIntegration(CCal);
        
    // redistribute atoms among processors to balance computation load
    this->DynamicLoadBalancing(CCal);
    
    // exchange data between processors 
    this->ProcessorsCommunication(CCal);
    
    // rebuild neighbor list if necessary
    this->NeighborListRebuild(CCal, istep);
                
    // calculate atomic forces
    CCal.PotentialEnergyForce(CCal.sys.e, CCal.sys.f, CCal.sys.x, CCal.sys.c, CCal.sys.q, CCal.app.muep, CCal.common.nmu);                     
    
    // // impose force constraints if any
    this->PostForceComputation(CCal);
        
    // final integration of the velocity vertlet method: update velocity after computing force    
    this->FinalIntegration(CCal);    
            
    // fix velocity and rescale velocity if necessary to control the temperature
    this->PostFinalIntegration(CCal);        
    
    // compute output, print on screen, and save it into binary files
    this->TimestepCompletion(CCal);           
  }    
}

void CIntegration::IntegrationSetup(CCalculation &CCal)
{
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
    int ensemblegroup = CCal.common.ensemblegroup;     
    int inum = CCal.common.inumgroup[ensemblegroup];        
    int biasflag = CCal.common.biasflag;
    int eta_mass_flag = CCal.common.eta_mass_flag;
    int mtchain = CCal.common.mtchain;
    int nc_tchain = CCal.common.nc_tchain;
    int *atomtype = CCal.nb.atomtype;            
    int *ilist = CCal.nb.atomgroups[ensemblegroup];
    
    InitialIntegrate(x, v, f, mass, dtarray, tarray, CCal.common.eta_mass, CCal.common.eta, 
           CCal.common.eta_dot, CCal.common.eta_dotdot, CCal.tmp.tmpmem, &CCal.tmp.tmpmem[inum],
            vlimitsq, atomtype, ilist, eta_mass_flag, biasflag, mtchain, nc_tchain, ensemblemode, dim, inum, backend);           
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
        
    int inum = CCal.common.nlocal;        
    // move out-of-box atoms back to the simulation box to impose periodic boundary conditions
    PBC(x, v, CCal.sys.image, CCal.app.boxhi, CCal.app.boxlo, CCal.app.boxhi_lamda, CCal.app.boxlo_lamda, 
            CCal.app.h, CCal.app.h_inv, CCal.app.h_rate, CCal.app.pbc, CCal.common.vdeform, 
            CCal.common.triclinic, dim, inum, backend);    
}

void CIntegration::DynamicLoadBalancing(CCalculation &CCal)
{
    
}

void CIntegration::ProcessorsCommunication(CCalculation &CCal)
{
    
}

void CIntegration::NeighborListRebuild(CCalculation &CCal, int istep)
{
    dstype *x = CCal.sys.x;
    dstype *xhold = CCal.sys.xhold;
    dstype *y = CCal.tmp.tmpmem;
    dstype skin = CCal.common.skin;
    
    int build = 0;    
    int ago = istep + 1;
    int delay = CCal.common.delay;
    int every = CCal.common.every;
    int distcheck = CCal.common.distcheck;
    int dim = CCal.common.dim;
    int backend = CCal.common.backend;
    int nlocal = CCal.common.nlocal;
            
    if ((ago >= delay) && (ago % every == 0)) {
        if (distcheck == 0) {
            build = 1;
        } else {
            ArrayDistSquareSum(y, x, xhold, dim, nlocal, backend);
            dstype maxdistsquare = ArrayMax(y, CCal.tmp.tmpmem, nlocal, backend);
            if (maxdistsquare > 0.25*skin*skin)
                build = 1;
        }        
    } 

    if (build == 1) {   
        CCal.NeighborList(x);
        ArrayCopy(xhold, x, nlocal*dim, backend);
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
    dstype *box = CCal.app.box;
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
    int ensemblegroup = CCal.common.ensemblegroup;     
    int inum = CCal.common.inumgroup[ensemblegroup];        
    int biasflag = CCal.common.biasflag;
    int eta_mass_flag = CCal.common.eta_mass_flag;
    int mtchain = CCal.common.mtchain;
    int nc_tchain = CCal.common.nc_tchain;
    int *atomtype = CCal.nb.atomtype;            
    int *ilist = CCal.nb.atomgroups[ensemblegroup];
    
    FinalIntegrate(x, v, f, mass, dtarray, tarray, CCal.common.eta_mass, CCal.common.eta, 
            CCal.common.eta_dot, CCal.common.eta_dotdot, CCal.tmp.tmpmem, &CCal.tmp.tmpmem[inum], 
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
            
    dstype *x = CCal.sys.x;
    dstype *e = CCal.sys.e;
    dstype *f = CCal.sys.f;
    dstype *v = CCal.sys.v;
    dstype *vatom = CCal.sys.vatom;
    dstype *mass = CCal.app.atommass;        
    dstype *tarray = CCal.common.tarray;
    dstype *scalars = CCal.common.scalars;
    dstype *tmpmem = CCal.tmp.tmpmem;
    
    int backend = CCal.common.backend;
    int dim = CCal.common.dim;
    int ensemblemode = CCal.common.ensemblemode; 
    int mtchain = CCal.common.mtchain;
    int *atomtype = CCal.nb.atomtype;            
    int *ilist = CCal.nb.atomgroups[0];
    int inum = CCal.common.inumgroup[0];        
            
    dstype mvv2e = CCal.common.mvv2e;            
    dstype boltz = CCal.common.boltz;            
    dstype nktv2p = CCal.common.nktv2p;    
    dstype volume = CCal.common.volume;    
    dstype masstotal = CCal.common.masstotal;    
            
    int ncol = 14;    
    // potential energy
    ArrayCopy(tmpmem, e, inum, backend);    
    // kinetic energy
    ComputeKEAtom(&tmpmem[inum], mass, v, mvv2e, atomtype, ilist, dim, inum, backend);
    // kinetic energy symmetric tensor
    ComputeStressAtom(&tmpmem[2*inum], mass, vatom, v, mvv2e, one, 
            atomtype, ilist, 0, 1, dim, inum, backend);   
    // per-atom virial    
    ArrayTranspose(&tmpmem[8*inum], vatom, 6, inum, backend);        
    
    // tmpmem : inum x ncol,  onevec : inum x 1
    dstype *onevec =  &tmpmem[ncol*inum];  // inum
    ArraySetValue(onevec, one, inum, backend);
    PGEMTV(CCal.common.cublasHandle, inum, ncol, &one, tmpmem, inum, onevec, 
            inc1, &one, &tmpmem[(ncol+1)*inum], inc1, &tmpmem[(ncol+1)*inum+ncol], backend);        
    
    // copy to host memory
    TemplateCopytoHost(scalars, &tmpmem[(ncol+1)*inum], ncol, backend);
        
    CCal.common.pe = scalars[0];
    CCal.common.ke = scalars[1];
    CCal.common.tdof = (CCal.common.natoms-1)*dim;
    CCal.common.temp = 2.0*CCal.common.ke/ (CCal.common.tdof * boltz);    
    for(int i=0; i<6; i++) 
        CCal.common.ke_tensor[i] = scalars[2+i];    
    
    for(int i=0; i<6; i++) 
        CCal.common.virial[i] = scalars[8+i];
    cpuComputePressureSymTensor(CCal.common.pres_tensor, CCal.common.virial, 
            CCal.common.ke_tensor, volume, nktv2p, dim);
    CCal.common.pres = cpuComputePressureScalar(CCal.common.virial, volume, CCal.common.temp, 
            CCal.common.tdof, boltz, nktv2p, dim);
    CCal.common.enthalpy = CCal.common.pe + CCal.common.ke + CCal.common.pres*volume;    

    // couple energy
    if (ensemblemode==2) 
        CCal.common.nvtenergy = cpuComputeNVTEnergy(tarray, CCal.common.eta, CCal.common.eta_mass,
                CCal.common.eta_dot, mtchain);        
    else 
        CCal.common.nvtenergy = 0.0;           
    dstype localce = CCal.common.nvtenergy + CCal.common.vrtenergy;    
#ifdef  HAVE_MPI          
#ifdef USE_FLOAT         
    MPI_Allreduce(&localce, &CCal.common.ce, 1, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);
#else            
    MPI_Allreduce(&localce, &CCal.common.ce, n, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#endif    
#endif
    
    printf("# atoms:  %d,   total mass: %g,   volume: %g\n",CCal.common.natoms, masstotal, volume);        
    printf("potential energy:  %g,   kinetic energy: %g,   couple energy: %g\n",CCal.common.pe, CCal.common.ke, CCal.common.ce);        
    printf("temperature:  %g,   total pressure: %g,   enthalpy: %g\n",CCal.common.temp, CCal.common.pres, CCal.common.enthalpy);            
    printf("temperature tensor: "); print1darray(CCal.common.ke_tensor, 6);
    printf("virial: "); print1darray(CCal.common.virial, 6);
    printf("pressure tensor: "); print1darray(CCal.common.pres_tensor, 6);
    
//     ComputeXCM(&scalars[10], tmpmem, x, &tmpmem[3*inum], mass, CCal.app.box, masstotal, 
//             ilist, atomtype, CCal.sys.image, CCal.common.triclinic, dim, inum, backend);
//     
//     ComputeVCM(&scalars[13], tmpmem, v, &tmpmem[3*inum], mass, 
//             masstotal, ilist, atomtype, dim, inum, backend);
    
//     printf("\nPotential energy :  %g,   kinetic energy : %g\n",CCal.common.currentstep, CCal.common.time);        
//     
//     ComputeXCM(xcm, axcm, x, tmp, mass, box, masstotal, ilist, type, image, triclinic, dim, inum, backend);
//     ComputeVCM(vcm, avcm, v, tmp, mass, masstotal, ilist, type, dim, inum);
//     gyr = ComputeGyration(ag, xcm, x, tmp, mass, box, masstotal, ilist, type, image, triclinic, dim, inum);
//     ComputeAngmom(lmom, p, xcm, x, v, tmp, mass, box, ilist, type, image, triclinic, dim, inum);
//     ComputeTorque(tq, q, xcm, x, f, tmp, box, ilist, image, triclinic, dim, inum);
//     ComputeInertia(inertia, ione, xcm, x, tmp, mass, box, ilist, type, image, triclinic, dim, inum);
//     
//     ComputeKEAtom(keatom, mass, v, mvv2e, type, ilist, dim, inum, backend);
//     ComputeStressAtom(stress, mass, vatom, v, mvv2e, nktv2p, type, ilist,  vflag, keflag, dim, inum, backend);
//     ComputeDisplaceAtom(displace, x, xoriginal, box, pbc, ilist,  triclinic, dim, inum);
//     
//     temp = ComputeTempScalar(keatom, v, tmp, mass, tfactor, type, ilist, dim, inum, backend);
//     press = ComputePressureScalar(virial, volume, temp, tempdof, boltz, nktv2p, dim, backend);
//     ComputeTempSymTensor(ke_tensor, stress, v, tmp, mass, tfactor, type, ilist, dim, inum);
//     ComputePressureSymTensor(press_tensor, virial, ke_tensor, volume, nktv2p, dim);
//     ComputeHeatFlux(vector, jc, ke, pe, stress, v, tmp, nktv2p, ilist,  pressatomflag, dim, inum);
//     ComputeMSD(msd, vec, x, xoriginal, box, xcm, tmp, ilist, image, naverage, avflag, triclinic, nmsd, dim,  inum);
//     ComputeVACF(vacf, vec, v, voriginal, tmp, ilist, nvacf, dim,  inum);
//     
//     ComputeCoordAtomCutoff(cvec, x, rcutsq, type, ilist, neighlist, neighnum, typelo, typehi, jgroupbit, dim, ntypes, jnum, inum);
//     ComputeCoordAtomCutoff(carray, x, rcutsq, type, ilist, neighlist, neighnum, typelo, typehi, jgroupbit, ncol, dim, ntypes, jnum, inum);
    
}

        
#endif        




