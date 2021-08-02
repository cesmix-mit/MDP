/***************************************************************************                               
                    Molecular Dynamics Potentials (MDP)
                           CESMIX-MIT Project  
 
 Contributing authors: Ngoc-Cuong Nguyen (cuongng@mit.edu, exapde@gmail.com)
 ***************************************************************************/

#ifndef MDP_COMMON_H
#define MDP_COMMON_H

using std::string;

#ifdef USE_CUDA
#include "cudamem.h"
#endif

#ifndef CUDA_SYNC
#define CUDA_SYNC
#endif

#ifndef USE_CUDA    
#define cublasHandle_t int
#define cudaEvent_t int
#endif
// #ifndef cublasHandle_t
// #define cublasHandle_t int
// #endif
// #ifndef cudaEvent_t
// #define cudaEvent_t int
// #endif

#define SDOT sdot_
#define SGEMV sgemv_
#define SGEMM sgemm_
#define SGETRF sgetrf_
#define SGETRI sgetri_

#define DDOT ddot_
#define DGEMV dgemv_
#define DGEMM dgemm_
#define DGETRF dgetrf_
#define DGETRI dgetri_

typedef int Int; 
#ifdef PREC_FLOAT
typedef float dstype;
#else
typedef double dstype; //  double is default precision 
#endif

#define NEIGHMASK 0x3FFFFFFF
#define MAX_ATOM_TYPES      30
#define MAX_MOLECULE_TYPES  20
#define MAX_MOLECULE_SIZE   10
#define MDPMIN(a,b) ((a) < (b) ? (a) : (b))
#define MDPMAX(a,b) ((a) > (b) ? (a) : (b))

#define CPUFREE(x)                                                           \
{                                                                         \
    if (x != NULL) {                                                      \
        free(x);                                                          \
        x = NULL;                                                         \
    }                                                                     \
}

extern "C" {
    double DNRM2(int*,double*,int*);
    double DDOT(int*,double*,int*,double*,int*);
    void DAXPY(int*,double*,double*,int*,double*,int*);
    void DGEMV(char*,int*,int*,double*,double*,int*,double*,int*,double*,double*,int*);  
    void DGEMM(char*,char*,int*,int*,int*,double*,double*,int*,
             double*,int*,double*,double*,int*);        
    void DGETRF(int*,int*,double*,int*,int*,int*);
    void DGETRI(int*,double*,int*,int*,double*,int*,int*);
    void DTRSM(char *, char*, char*, char *, int *, int *, double*, double*, int*,
             double*, int*);

    float SNRM2(int*,float*,int*);  
    float SDOT(int*,float*,int*,float*,int*);
    void SAXPY(int*,float*,float*,int*,float*,int*);
    void SGEMM(char*,char*,int*,int*,int*,float*,float*,int*,
             float*,int*,float*,float*,int*);  
    void SGEMV(char*,int*,int*,float*,float*,int*,float*,int*,float*,float*,int*);      
    void SGETRF(int*,int*,float*,int*,int*,int*);    
    void SGETRI(int*,float*,int*,int*,float*,int*,int*);
    void STRSM(char *, char*, char*, char *, int *, int*, float*, float*, int*,
             float*, int*);    
}

// // global variables for BLAS  
dstype one = 1.0;
dstype minusone = -1.0;
dstype zero = 0.0;
char chn = 'N';
char cht = 'T';
char chl = 'L';
char chu = 'U';
char chr = 'R';
int inc1 = 1;

// global variables for CUBLAS  
dstype cublasOne[1] = {one};
dstype cublasMinusone[1] = {minusone};
dstype cublasZero[1] = {zero};

#ifdef TIMING    
    #define INIT_TIMING auto begin = chrono::high_resolution_clock::now(); auto end = chrono::high_resolution_clock::now();
#else
    #define INIT_TIMING
#endif

#ifdef TIMING
   #define TIMING_START  begin = chrono::high_resolution_clock::now();   
#else 
   #define TIMING_START     
#endif       

#ifdef TIMING
   #define TIMING_END    end = chrono::high_resolution_clock::now();   
#else 
   #define TIMING_END     
#endif       

#ifdef TIMING       
   #define TIMING_GET(num) common.timing[num] += chrono::duration_cast<chrono::nanoseconds>(end-begin).count()/1e6;        
#else 
   #define TIMING_GET(num)  
#endif                      

#ifdef TIMING
   #define START_TIMING {CUDA_SYNC; TIMING_START;}       
#else 
   #define START_TIMING
#endif       

#ifdef TIMING
   #define END_TIMING(num) {CUDA_SYNC; TIMING_END; TIMING_GET(num)}   
#else 
   #define END_TIMING(num)
#endif       

#ifdef TIMING       
   #define TIMING_GET1(num) disc.common.timing[num] += chrono::duration_cast<chrono::nanoseconds>(end-begin).count()/1e6;        
#else 
   #define TIMING_GET1(num)  
#endif                      

#ifdef TIMING
   #define END_TIMING_DISC(num) {CUDA_SYNC; TIMING_END; TIMING_GET1(num)}   
#else 
   #define END_TIMING_DISC(num)
#endif       
                
template <typename T> static void TemplateMalloc(T **data, int n, int backend)
{
    if (backend == 0)       // One thread CPU
        // allocate the memory on the CPU        
        *data = (T *) malloc(n*sizeof(T));     
    if (backend == 1)  // Open MP
        // allocate the memory on the CPU        
        *data = (T *) malloc(n*sizeof(T));    
#ifdef USE_CUDA            
    if (backend == 2)  // CUDA C                
        // allocate the memory on the GPU            
        CUDA_CHECK( cudaMalloc( (void**)data, n * sizeof(T) ) );
#endif                  
}

template <typename T> static void TemplateCopytoDevice(T *d_data, T *h_data, int n, int backend)
{
    if (backend == 0)       
        //cpuArrayCopy(d_data, h_data, n);
        for (int i=0; i<n; i++) d_data[i] = h_data[i];
    if (backend == 1)         
        for (int i=0; i<n; i++) d_data[i] = h_data[i];
#ifdef USE_CUDA            
    if (backend == 2)  // CUDA          
        cudaCopytoDevice(d_data, h_data, n);
#endif                  
#ifdef USE_HIP            
    if (backend == 3)  // HIP
        hipCopytoDevice(d_data, h_data, n);
#endif                      
}

template <typename T> static void TemplateCopytoHost(T *h_data, T *d_data, int n, int backend)
{
    if (backend == 0)       
        //cpuArrayCopy(h_data, d_data, n);
        for (int i=0; i<n; i++) h_data[i] = d_data[i];
    if (backend == 1)         
        for (int i=0; i<n; i++) h_data[i] = d_data[i];
#ifdef USE_CUDA            
    if (backend == 2)  // CUDA          
        cudaCopytoHost(h_data, d_data, n);
#endif                  
#ifdef USE_HIP            
    if (backend == 3)  // HIP
        hipCopytoHost(h_data, d_data, n);
#endif                      
}

template <typename T> static void TemplateFree(T *data, int backend)
{
    if (backend == 0)       
        CPUFREE(data);
    if (backend == 1)         
        CPUFREE(data);
#ifdef USE_CUDA            
    if (backend == 2)  // CUDA          
        GPUFREE(data);
#endif                  
#ifdef USE_HIP            
    if (backend == 3)  // HIP          
        HIPFREE(data);
#endif                      
}

struct commonstruct {     
    string filein;       // Name of binary file with input data
    string fileout;      // Name of binary file to write the solution    
            
    int backend=1;   // 1: CPU; 2: CUDA GPU  
    int mpiRank=0;  // MPI rank      
    int mpiProcs=1; // number of MPI ranks
    
    int dim;     // physical dimensions
    int nconfigs=1; // number of configurations        
    int natomtypes=1;     // number of atom types     
    int nmoletypes=1;     // number of molecule types         
    int K=0;       // order of radial basis functions
    int L=0;       // order of spherical harmonics      
    int M=0;       // total number of potentials 
    int Nub;     // number of unique bispectrum components
    int Npower; // number of power spectrum components
    int Nbispectrum; // number of bispectrum components
    int Nbf=0;    // Nbf = Npower+Nbispectrum
    int Ncoeff=0; // number of machine learning descriptors coefficients
    int Ncg=0;// the total number of non-zero Clebsch-Gordan coefficients 
    int descriptor=0;   // descriptor flag: 0 -> Spherical Harmonics Bessel
    int spectrum=1;     // spectrum flag: 0-> power spectrum, 1-> bispectrum, 2-> both power and bispectrum 
    int training=0;     // 0 -> no training, 1 -> Linear regression, 2 -> Gaussian process, 3 -> Neural net
    int dftdata;      // 0 -> no data, 1 -> energies only, 2 -> forces only, 3 -> energies and forces
    int runMD;        // 0 no MD simulation, 1 -> run MD simulation
    int potential;    // 0 -> empirical potential, 1 -> empirical + LR, 2 -> empirical + GP, 3 -> empirical + NN
    int cutofftype=0; // 0 -> single cut-off raidus for all atoms, 1 -> multiple cut-off radii for atom pairs  
    int neighpair=0;  // 0 -> full neighbor list, 1 -> half neighbor list
    int neighcell=0;  // 0 -> O(N^2) algorithm, 1 -> Cell-linked list algorithm to form neighbor list
    int decomposition=0;// 0 -> force decomposition, 1 -> atom decomposition
    int bondtype=0;   // 0 -> non-bonded interaction, 1 -> bonded interaction
    int chemtype = 0;   // 0 -> single atom-type basis functions, 1 -> double atom-type basis functions 
    int pairsymmetry;  // 1 -> V(r_ij) equal V(r_ji), 0 -> V(r_ij) not equal V(r_ji) 
    int tripletsymmetry; // 1 -> V(r_ij, r_ik) equal V(r_ik, r_ij)
                         // 0 -> V(r_ij, r_ik) not equal V(r_ik, r_ij)
    int energycal;    // turns energy calculation on or off
    int forcecal;     // turns force calculation on or off
    int stresscal;     // turns stress calculation on or off
    int triclinic;     // types of the simulation box  
    int unitstyle;
    int vflag=0;
    int vdeform=0;
    
    int nflags=0;
    int nsolversparam=0;   
    int nsimulaparam=0;   
    int neta=0;   
    int nkappa=0;   
    int nmuep=0;
    int nmuml=0;
    int nmu1a=0;
    int nmu1b=0;
    int nmu2a=0;
    int nmu2b=0;
    int nmu2c=0;
    int nmu3a=0;
    int nmu3b=0;
    int nmu3c=0;
    int nmu4a=0;
    int nmu4b=0;    
    int Nempot=0; // number of empirical potentials
    int npot1a=0;
    int npot1b=0;
    int npot2a=0;
    int npot2b=0;
    int npot2c=0;
    int npot3a=0;
    int npot3b=0;
    int npot3c=0;
    int npot4a=0;
    int npot4b=0;    
    int natom1b=0;
    int natom2b=0;
    int natom2c=0;
    int natom3b=0;
    int natom3c=0;
    int natom4b=0;    
    int nmu[12];
    
    // cpuFixSetForce, cpuFixLineForce, cpuFixPlaneForce, cpuFixAddForce, cpuFixAveForce, cpuFixDragForce, 
    // cpuFixWallReflect, cpuFixWallHarmonic, cpuFixWallLJ93, cpuFixWallLJ126, cpuFixWallLJ1043, cpuFixWallMorse 

    // number of force constraints 
    int nsetforce=0; 
    int nlineforce=0;
    int nplaneforce=0;
    int naddforce=0; 
    int naveforce=0; 
    int ndragforce=0;
    int ngravityforce=0;
    int nwallreflect=0;
    int nwallharmonic=0;
    int nwalllj93=0;
    int nwalllj126=0;
    int nwalllj1043=0;
    int nwallmorse=0;
    int nfixforce=0;
     
    // list of atom groups for each force constraint
    int *gsetforce=NULL; 
    int *glineforce=NULL;
    int *gplaneforce=NULL;
    int *gaddforce=NULL; 
    int *gaveforce=NULL; 
    int *gdragforce=NULL;
    int *ggravityforce=NULL;
    int *gwallreflect=NULL;
    int *gwallharmonic=NULL;
    int *gwalllj93=NULL;
    int *gwalllj126=NULL;
    int *gwalllj1043=NULL;
    int *gwallmorse=NULL;
    int *gfixforce=NULL;
    
    // list of integer parameters for every force constraint
    int *isetforce=NULL; 
    int *ilineforce=NULL;
    int *iplaneforce=NULL;
    int *iaddforce=NULL; 
    int *iaveforce=NULL; 
    int *idragforce=NULL;
    int *igravityforce=NULL;
    int *iwallreflect=NULL;
    int *iwallharmonic=NULL;
    int *iwalllj93=NULL;
    int *iwalllj126=NULL;
    int *iwalllj1043=NULL;
    int *iwallmorse=NULL;
    int *ifixforce=NULL;

    // list of float parameters for every force constraint
    dstype *fsetforce=NULL; 
    dstype *flineforce=NULL;
    dstype *fplaneforce=NULL;
    dstype *faddforce=NULL; 
    dstype *faveforce=NULL; 
    dstype *fdragforce=NULL;
    dstype *fgravityforce=NULL;
    dstype *fwallreflect=NULL;
    dstype *fwallharmonic=NULL;
    dstype *fwalllj93=NULL;
    dstype *fwalllj126=NULL;
    dstype *fwalllj1043=NULL;
    dstype *fwallmorse=NULL;
    dstype *ffixforce=NULL;
    
    // number of velocity constraints 
    int nsetvelocity=0;
    int nfixvelocity=0;
    
    // list of atom groups for each velocity constraint
    int *gsetvelocity=NULL;
    int *gfixvelocity=NULL;
    
    // list of integer parameters for every velocity constraint
    int *isetvelocity=NULL;
    int *ifixvelocity=NULL;
    
    // list of float parameters for every velocity constraint
    dstype *fsetvelocity=NULL;
    dstype *ffixvelocity=NULL;
        
    // time integration parameters
    dstype dtarray[10];
    dstype tarray[10];
    dstype eta[10];
    dstype eta_dot[10];
    dstype eta_dotdot[10];
    dstype eta_mass[10];
    dstype nvtenergy;
    dstype vlimitsq;
    int eta_mass_flag;
    int biasflag;
    int mtchain=1;
    int nc_tchain=1;
    int ensemblemode=0;
    int ensemblegroup=0;
    int delay=10;
    int every=1;
    int distcheck=1;
                    
    int ngroups=1; // number of atom groups
    int *inumgroup=NULL; // a list of numbers of atoms in every group 
    
    int natoms; // number of atoms in the simulation box
    int nlocal; // number of atoms in the current processor    
    int inum;  // number of atoms in the group    
    int inummax;  // maximum number of atoms in the group
    int jnum;  // maximum number of neighbors    
    int gnum;  // number of ghost atoms
    int anum;  // (nlocal+gnum) total number of atoms 
    int anummax; // maximum  number of atoms allowed     
    int pnum;  // number of periodic images
    int cnum;  // number of cells    
    int nintmem;
    int ntmpmem;
    int ne; 
    int nt; 
    int nx; 
    int nv;
    int nf;
    int nq;
    int ncq;
    int nba; // number of blocks of atoms
    int nab; // number of atoms per block
    int nabmax; // maxmimum number of atoms per block
    int ablks[17];
        
    int ntimesteps; // number of time steps
    int currentstep=0; // current timestep
    dstype dt = 0.005;   // timestep size     
    dstype time = 0.0; // current simulation time    
    dstype skin = 0.3;
    dstype rcutml; // cutoff-radius for machine learning potential;
    dstype rbbvol; // volume of reference bounding box
    dstype volume; // volume of the simulation box
    
    dstype pe = 0.0; // potential energy
    dstype ke = 0.0; // kinetic energy
    dstype ce = 0.0; // couple energy
    dstype masstotal = 0.0; // total mass
    dstype temp = 0.0;
    dstype pres = 0.0;
    dstype enthalpy = 0.0;
    dstype tdof;
    dstype extra_dof;

    // default LJ unit system
    dstype boltz = 1.0;
    dstype hplanck = 1.0;
    dstype mvv2e = 1.0;
    dstype ftm2v = 1.0;
    dstype mv2d = 1.0;
    dstype nktv2p = 1.0;
    dstype qqr2e = 1.0;
    dstype qe2f = 1.0;
    dstype vxmu2f = 1.0;
    dstype xxt2kmu = 1.0;
    dstype e_mass = 0.0;    // not yet set
    dstype hhmrr2e = 0.0;
    dstype mvh2r = 0.0;
    dstype angstrom = 1.0;
    dstype femtosecond = 1.0;
    dstype qelectron = 1.0;    
    
    dstype vrtenergy=0.0;
    dstype second[1];
    int seed[1];
    int save[1];
    int vrtmode=-1;    
            
    dstype boxoffset[3];
    dstype scalars[32];
    dstype virial[6];
    dstype ke_tensor[6];
    dstype pres_tensor[6];
    int pbc[3];
    
    int *pot1a=NULL; 
    int *pot1b=NULL;
    int *pot2a=NULL; 
    int *pot2b=NULL;
    int *pot2c=NULL; 
    int *pot3a=NULL;
    int *pot3b=NULL; 
    int *pot3c=NULL;
    int *pot4a=NULL;
    int *pot4b=NULL; 
    int *atom1b=NULL;
    int *atom2b=NULL;
    int *atom2c=NULL; 
    int *atom3b=NULL; 
    int *atom3c=NULL;
    int *atom4b=NULL; 
    
    int *traininglist;
    int *validatelist;
    int trainingnum=0;
    int validatenum=0;
            
//     dstype pi = M_PI; // Pi  number
     
    cudaEvent_t eventHandle;
    cublasHandle_t cublasHandle;
    
#ifdef  HAVE_MPI
    MPI_Request * requests;
    MPI_Status * statuses;
#endif
    
    void freememory()
    {
        CPUFREE(pot1a);
        CPUFREE(pot1b);
        CPUFREE(pot2a);
        CPUFREE(pot2b);
        CPUFREE(pot2c);
        CPUFREE(pot3a);
        CPUFREE(pot3b);
        CPUFREE(pot3c);
        CPUFREE(pot4a);
        CPUFREE(pot4b);            
        CPUFREE(atom1b);
        CPUFREE(atom2b);
        CPUFREE(atom2c);
        CPUFREE(atom3b);
        CPUFREE(atom3c);
        CPUFREE(atom4b); 
        CPUFREE(traininglist); 
        CPUFREE(validatelist); 
    }                         
};

struct appstruct {              
    int *lsize=NULL;
    int *nsize=NULL;  // input data size    
    int *ndims=NULL;  // dimensions
    int *flags=NULL;        // flag parameters            
    int *bcs=NULL;           // boundary conditions
    int *pbc=NULL;           // periodic boundary conditions        
        
    dstype *boxhi=NULL;
    dstype *boxlo=NULL;
    dstype *boxhi_lamda=NULL;
    dstype *boxlo_lamda=NULL;
    dstype *h=NULL;
    dstype *h_inv=NULL;
    dstype *h_rate=NULL;
    dstype *box; // six parameters for the simulation box
    dstype *boxoffset=NULL;                
    dstype *atommass=NULL; //  a list of atomic mass for every atom type
    dstype *atomcharge=NULL; //  a list of atomic charge for every atom type
    dstype *simulaparam=NULL; // simulation parameters   
    dstype *solversparam=NULL; // solvers parameters              
    dstype *physicsparam=NULL; // general physical parameters
    dstype *nvtparam=NULL; // NVT parameters         
    dstype *eta=NULL; // hyperparameters
    int *atomnumber=NULL;   // a list of atomic numbers for every atom type
    int *kappa=NULL; // integer parameters                  
              
    dstype *muep=NULL;
    dstype *muml=NULL; 
    dstype *mu1a=NULL; 
    dstype *mu1b=NULL;
    dstype *mu2a=NULL; 
    dstype *mu2b=NULL;
    dstype *mu2c=NULL; 
    dstype *mu3a=NULL;
    dstype *mu3b=NULL; 
    dstype *mu3c=NULL;
    dstype *mu4a=NULL;
    dstype *mu4b=NULL;     
    int *pot1a=NULL; 
    int *pot1b=NULL;
    int *pot2a=NULL; 
    int *pot2b=NULL;
    int *pot2c=NULL; 
    int *pot3a=NULL;
    int *pot3b=NULL; 
    int *pot3c=NULL;
    int *pot4a=NULL;
    int *pot4b=NULL; 
    dstype *rcutsqml=NULL;
    dstype *rcutsq2a=NULL; 
    dstype *rcutsq2b=NULL;
    dstype *rcutsq2c=NULL; 
    dstype *rcutsq3a=NULL;
    dstype *rcutsq3b=NULL; 
    dstype *rcutsq3c=NULL;
    dstype *rcutsq4a=NULL;
    dstype *rcutsq4b=NULL;     
    dstype *rcutsq=NULL;  // cut-off radius for neighbor list formation
    int *atom1b=NULL;
    int *atom2b=NULL;
    int *atom2c=NULL; 
    int *atom3b=NULL; 
    int *atom3c=NULL;
    int *atom4b=NULL; 
    
    // list of integer parameters for every force constraint
    int *isetforce=NULL; 
    int *ilineforce=NULL;
    int *iplaneforce=NULL;
    int *iaddforce=NULL; 
    int *iaveforce=NULL; 
    int *idragforce=NULL;
    int *igravityforce=NULL;
    int *iwallreflect=NULL;
    int *iwallharmonic=NULL;
    int *iwalllj93=NULL;
    int *iwalllj126=NULL;
    int *iwalllj1043=NULL;
    int *iwallmorse=NULL;
    int *ifixforce=NULL;

    // list of float parameters for every force constraint
    dstype *fsetforce=NULL; 
    dstype *flineforce=NULL;
    dstype *fplaneforce=NULL;
    dstype *faddforce=NULL; 
    dstype *faveforce=NULL; 
    dstype *fdragforce=NULL;
    dstype *fgravityforce=NULL;
    dstype *fwallreflect=NULL;
    dstype *fwallharmonic=NULL;
    dstype *fwalllj93=NULL;
    dstype *fwalllj126=NULL;
    dstype *fwalllj1043=NULL;
    dstype *fwallmorse=NULL;
    dstype *ffixforce=NULL;
    
    // list of integer parameters for every velocity constraint
    int *isetvelocity=NULL;
    int *ifixvelocity=NULL;
    
    // list of float parameters for every velocity constraint
    dstype *fsetvelocity=NULL;
    dstype *ffixvelocity=NULL;
    
    // custom destructor
    void freememory(int backend)
    {
        TemplateFree(lsize, backend);
        TemplateFree(nsize, backend);
        TemplateFree(ndims, backend);   
        TemplateFree(flags, backend);    
        TemplateFree(bcs, backend);    
        TemplateFree(pbc, backend);    
        TemplateFree(boxoffset, backend);
        TemplateFree(atomnumber, backend);    
        TemplateFree(atommass, backend);    
        TemplateFree(atomcharge, backend);
        TemplateFree(simulaparam, backend);
        TemplateFree(solversparam, backend);
        TemplateFree(physicsparam, backend);
        TemplateFree(eta, backend);
        TemplateFree(kappa, backend);
        TemplateFree(muml, backend);
        TemplateFree(muep, backend);
        TemplateFree(mu1a, backend);
        TemplateFree(mu1b, backend);
        TemplateFree(mu2a, backend);
        TemplateFree(mu2b, backend);
        TemplateFree(mu2c, backend);
        TemplateFree(mu3a, backend);
        TemplateFree(mu3b, backend);
        TemplateFree(mu3c, backend);
        TemplateFree(mu4a, backend);
        TemplateFree(mu4b, backend);            
        TemplateFree(pot1a, backend);
        TemplateFree(pot1b, backend);
        TemplateFree(pot2a, backend);
        TemplateFree(pot2b, backend);
        TemplateFree(pot2c, backend);
        TemplateFree(pot3a, backend);
        TemplateFree(pot3b, backend);
        TemplateFree(pot3c, backend);
        TemplateFree(pot4a, backend);
        TemplateFree(pot4b, backend);            
        TemplateFree(rcutsqml, backend);
        TemplateFree(rcutsq2a, backend);
        TemplateFree(rcutsq2b, backend);
        TemplateFree(rcutsq2c, backend);
        TemplateFree(rcutsq3a, backend);
        TemplateFree(rcutsq3b, backend);
        TemplateFree(rcutsq3c, backend);
        TemplateFree(rcutsq4a, backend);
        TemplateFree(rcutsq4b, backend);                        
        TemplateFree(rcutsq, backend);            
        TemplateFree(atom1b, backend);
        TemplateFree(atom2b, backend);
        TemplateFree(atom2c, backend);
        TemplateFree(atom3b, backend);
        TemplateFree(atom3c, backend);
        TemplateFree(atom4b, backend);                        
    }
};

struct configstruct {    
    int *lsize=NULL;
    int *nsize=NULL;  // data size
    int *ndims=NULL;  // dimensions    
    int *natoms=NULL; // a list containing the number of atoms in each configuration 
    int *natomssum=NULL; // a list containing the number of atoms in each configuration     
    int *t=NULL;      // atomic types of atoms for all configurations 
    dstype *x=NULL;   // positions of atoms for all configurations        
    dstype *v=NULL;   // velocities of atoms for all configurations        
    dstype *e=NULL;   // energies acting on atoms for all configurations    
    dstype *f=NULL;   // forces acting on atoms for all configurations    
    dstype *q=NULL;   // charges acting on atoms for all configurations        
    dstype *a=NULL;   // principal vectors of the simulation box for all configurations       
    dstype *b=NULL;   // principal vectors of the simulation box for all configurations       
    dstype *c=NULL;   // principal vectors of the simulation box for all configurations       
    dstype *we=NULL;   // weights on energies
    dstype *wf=NULL;   // weights on forces 
    
    void freememory()
    {        
        CPUFREE(lsize);
        CPUFREE(nsize);
        CPUFREE(ndims);   
        CPUFREE(natoms);
        CPUFREE(natomssum);
        CPUFREE(t);
        CPUFREE(x);
        CPUFREE(v);
        CPUFREE(e);
        CPUFREE(f);
        CPUFREE(q);           
        CPUFREE(a);
        CPUFREE(b);
        CPUFREE(c);
        CPUFREE(we);
        CPUFREE(wf);
    }                         
};

struct neighborstruct {      
    dstype *a=NULL;
    dstype *b=NULL;
    dstype *c=NULL;        
    int *cellnum=NULL;
    dstype *cellsize=NULL;
    dstype *eta1=NULL;
    dstype *eta2=NULL;
    dstype *eta3=NULL;        
    dstype *boxvertices=NULL;// vertices of the simulation box    
    dstype *bbvertices=NULL; // vertices of the bounding box    
    dstype *refvertices=NULL;// vertices of the reference box    
    dstype *rbvertices=NULL; // vertices of the reference bounding box  
    dstype *s2rmap=NULL;     // map simulation domain to reference domain 
    dstype *pimages=NULL;    // coordinates of periodic images     
    
    int **atomgroups;     // groups of atoms       
    int *atomtype=NULL;  // type of each atom i      
    int *alist=NULL;     // list of atoms i in the processor (inclduding own and ghost atoms)
    int *neighnum=NULL;  // numbers of neighbors for each atom i 
    int *neighlist=NULL; // list of neighbors for each atom i    
    int *mask=NULL;   
    int *tag=NULL;       
    
    void freememory(int backend)
    {
        TemplateFree(a, backend); 
        TemplateFree(b, backend); 
        TemplateFree(c, backend); 
        TemplateFree(cellnum, backend); 
        TemplateFree(cellsize, backend); 
        TemplateFree(eta1, backend); 
        TemplateFree(eta2, backend); 
        TemplateFree(eta3, backend);
        TemplateFree(boxvertices, backend); 
        TemplateFree(bbvertices, backend); 
        TemplateFree(refvertices, backend); 
        TemplateFree(rbvertices, backend); 
        TemplateFree(s2rmap, backend); 
        TemplateFree(pimages, backend); 
        TemplateFree(neighnum, backend); 
        TemplateFree(neighlist, backend);
        TemplateFree(alist, backend);
        TemplateFree(atomtype, backend); 
        TemplateFree(mask, backend); 
        TemplateFree(tag, backend); 
    }
};

struct sysstruct {        
    dstype *x=NULL;    // cartesian coordinates of local atoms in processor mpiRank        
    dstype *v=NULL;    // velocities of local atoms in processor mpiRank        
    dstype *e=NULL;    // per-atom energies 
    dstype *ee=NULL;   // per-atom energies for empirical descriptors
    dstype *f=NULL;    // atomic forces of local atoms in processor mpiRank        
    dstype *q=NULL;    // atomic charges of local atoms in processor mpiRank               
    dstype *s=NULL;    // stresses    
    dstype *d=NULL;    // descriptors
    dstype *dd=NULL;   // derivatives of descriptors
    dstype *c=NULL;    // vector of coeffcients associated with the basis functions            
    dstype *A=NULL;    // Regression matrix A 
    dstype *b=NULL;    // Regression vector b
    dstype *eatom=NULL;
    dstype *vatom=NULL;
    dstype *xhold=NULL;
    int *image=NULL;   
    
    void freememory(int backend)
    {
        TemplateFree(x, backend); 
        TemplateFree(v, backend); 
        TemplateFree(e, backend); 
        TemplateFree(ee, backend); 
        TemplateFree(f, backend); 
        TemplateFree(q, backend); 
        TemplateFree(s, backend);             
        TemplateFree(d, backend);             
        TemplateFree(dd, backend);             
        TemplateFree(c, backend);       
        TemplateFree(A, backend);       
        TemplateFree(b, backend);              
    }                     
};
  
struct tempstruct {
    int *intmem=NULL;     
    dstype *tmpmem=NULL;     
    dstype *buffrecv=NULL;
    dstype *buffsend=NULL;
    
    void freememory(int backend)
    {
        TemplateFree(intmem, backend); 
        TemplateFree(tmpmem, backend); 
        TemplateFree(buffrecv, backend); 
        TemplateFree(buffsend, backend); 
    }                     
};

struct shstruct {  
    int L;  // the maximum degree of spherical harmonics 
    int K;  // the number of zeros of spherical Bessel functions
    int Nub;// number of non-zero unqiue bispectrum compoments of spherical harmonics 
    int Ncg;// the total number of non-zero Clebsch-Gordan coefficients 
    int npower;
    int nbispectrum;
    int nbasis;
    
    int *indk=NULL;
    int *indl=NULL;        
    int *indm=NULL;       
    int *rowm=NULL;    
  
    dstype *fac=NULL;
    dstype *cg=NULL;
    dstype *x0=NULL;
    dstype *f=NULL;
    dstype *P=NULL;
    dstype *tmp=NULL;
    dstype *df=NULL;
    dstype *dP=NULL;
    dstype *dtmp=NULL;   
    
    void freememory(int backend)
    {
        TemplateFree(indk,backend); 
        TemplateFree(indl,backend); 
        TemplateFree(indm,backend); 
        TemplateFree(rowm,backend); 
        TemplateFree(fac,backend); 
        TemplateFree(cg,backend); 
        TemplateFree(x0,backend);             
        TemplateFree(f,backend);             
        TemplateFree(P,backend);             
        TemplateFree(tmp,backend);             
        TemplateFree(df,backend);             
        TemplateFree(dP,backend);             
        TemplateFree(dtmp,backend);                     
    }                         
};

struct snastruct {        
    int twojmax;
    int ncoeff;
    int ncoeffall;
    int nperdim;
    int idxb_max;
    int idxu_max;
    int idxz_max;
    int idxcg_max;
    int ntypes;
    int nelements;    
    int ndoubles;   // number of multi-element pairs
    int ntriples;   // number of multi-element triplets      
    int beta_max;                 // length of beta
    int bnormflag;
    int chemflag;    
    int quadraticflag;
    int switchflag;
    int bzeroflag;
    int wselfallflag;
    int rcutfacflag;
    int twojmaxflag; // flags for required parameters
    
    dstype wself;
    dstype rmin0;
    dstype rfac0;
    dstype rcutfac;
    dstype rcutmax;    
        
    int *map=NULL;  // map types to [0,nelements)    
    int *element=NULL;  // index on [0,nelements)
    int *idx_max=NULL; 
    int *idxz=NULL;
    int *idxz_block=NULL;
    int *idxb=NULL;
    int *idxb_block=NULL;
    int *idxu_block=NULL;
    int *idxcg_block=NULL;
    
    dstype *rcutsq=NULL;    
    dstype *radelem=NULL;
    dstype *wjelem=NULL; 
    dstype *bzero=NULL;
    dstype *coeffelem=NULL;           // element bispectrum coefficients
    dstype *fac=NULL;
    dstype *rootpqarray=NULL; 
    dstype *cglist=NULL;
    
    void freememory(int backend)
    {   
        TemplateFree(map, backend);
        TemplateFree(element, backend);
        TemplateFree(idx_max, backend);
        TemplateFree(idxz, backend);
        TemplateFree(idxb, backend);
        TemplateFree(idxb_block, backend);
        TemplateFree(idxu_block, backend);
        TemplateFree(idxz_block, backend);
        TemplateFree(idxcg_block, backend);
        
        TemplateFree(rootpqarray, backend);
        TemplateFree(cglist, backend);
        TemplateFree(fac, backend);
        TemplateFree(bzero, backend);
        TemplateFree(coeffelem, backend);
        TemplateFree(wjelem, backend);
        TemplateFree(radelem, backend);
        TemplateFree(rcutsq, backend);
    }                         
};

const double factable[] = {
  1,
  1,
  2,
  6,
  24,
  120,
  720,
  5040,
  40320,
  362880,
  3628800,
  39916800,
  479001600,
  6227020800,
  87178291200,
  1307674368000,
  20922789888000,
  355687428096000,
  6.402373705728e+15,
  1.21645100408832e+17,
  2.43290200817664e+18,
  5.10909421717094e+19,
  1.12400072777761e+21,
  2.5852016738885e+22,
  6.20448401733239e+23,
  1.5511210043331e+25,
  4.03291461126606e+26,
  1.08888694504184e+28,
  3.04888344611714e+29,
  8.8417619937397e+30,
  2.65252859812191e+32,
  8.22283865417792e+33,
  2.63130836933694e+35,
  8.68331761881189e+36,
  2.95232799039604e+38,
  1.03331479663861e+40,
  3.71993326789901e+41,
  1.37637530912263e+43,
  5.23022617466601e+44,
  2.03978820811974e+46,
  8.15915283247898e+47,
  3.34525266131638e+49,
  1.40500611775288e+51,
  6.04152630633738e+52,
  2.65827157478845e+54,
  1.1962222086548e+56,
  5.50262215981209e+57,
  2.58623241511168e+59,
  1.24139155925361e+61,
  6.08281864034268e+62,
  3.04140932017134e+64,
  1.55111875328738e+66,
  8.06581751709439e+67,
  4.27488328406003e+69,
  2.30843697339241e+71,
  1.26964033536583e+73,
  7.10998587804863e+74,
  4.05269195048772e+76,
  2.35056133128288e+78,
  1.3868311854569e+80,
  8.32098711274139e+81,
  5.07580213877225e+83,
  3.14699732603879e+85,
  1.98260831540444e+87,
  1.26886932185884e+89,
  8.24765059208247e+90,
  5.44344939077443e+92,
  3.64711109181887e+94,
  2.48003554243683e+96,
  1.71122452428141e+98,
  1.19785716699699e+100,
  8.50478588567862e+101,
  6.12344583768861e+103,
  4.47011546151268e+105,
  3.30788544151939e+107,
  2.48091408113954e+109,
  1.88549470166605e+111,
  1.45183092028286e+113,
  1.13242811782063e+115,
  8.94618213078297e+116,
  7.15694570462638e+118,
  5.79712602074737e+120,
  4.75364333701284e+122,
  3.94552396972066e+124,
  3.31424013456535e+126,
  2.81710411438055e+128,
  2.42270953836727e+130,
  2.10775729837953e+132,
  1.85482642257398e+134,
  1.65079551609085e+136,
  1.48571596448176e+138,
  1.3520015276784e+140,
  1.24384140546413e+142,
  1.15677250708164e+144,
  1.08736615665674e+146,
  1.03299784882391e+148,
  9.91677934870949e+149,
  9.61927596824821e+151,
  9.42689044888324e+153,
  9.33262154439441e+155,
  9.33262154439441e+157,
  9.42594775983835e+159,
  9.61446671503512e+161,
  9.90290071648618e+163,
  1.02990167451456e+166,
  1.08139675824029e+168,
  1.14628056373471e+170,
  1.22652020319614e+172,
  1.32464181945183e+174,
  1.44385958320249e+176,
  1.58824554152274e+178,
  1.76295255109024e+180,
  1.97450685722107e+182,
  2.23119274865981e+184,
  2.54355973347219e+186,
  2.92509369349301e+188,
  3.3931086844519e+190,
  3.96993716080872e+192,
  4.68452584975429e+194,
  5.5745857612076e+196,
  6.68950291344912e+198,
  8.09429852527344e+200,
  9.8750442008336e+202,
  1.21463043670253e+205,
  1.50614174151114e+207,
  1.88267717688893e+209,
  2.37217324288005e+211,
  3.01266001845766e+213,
  3.8562048236258e+215,
  4.97450422247729e+217,
  6.46685548922047e+219,
  8.47158069087882e+221,
  1.118248651196e+224,
  1.48727070609069e+226,
  1.99294274616152e+228,
  2.69047270731805e+230,
  3.65904288195255e+232,
  5.01288874827499e+234,
  6.91778647261949e+236,
  9.61572319694109e+238,
  1.34620124757175e+241,
  1.89814375907617e+243,
  2.69536413788816e+245,
  3.85437071718007e+247,
  5.5502938327393e+249,
  8.04792605747199e+251,
  1.17499720439091e+254,
  1.72724589045464e+256,
  2.55632391787286e+258,
  3.80892263763057e+260,
  5.71338395644585e+262,
  8.62720977423323e+264,
  1.31133588568345e+267,
  2.00634390509568e+269,
  3.08976961384735e+271,
  4.78914290146339e+273,
  7.47106292628289e+275,
  1.17295687942641e+278,
  1.85327186949373e+280,
  2.94670227249504e+282,
  4.71472363599206e+284,
  7.59070505394721e+286,
  1.22969421873945e+289,
  2.0044015765453e+291,
  3.28721858553429e+293,
  5.42391066613159e+295,
  9.00369170577843e+297,
  1.503616514865e+300,
};

#endif  
