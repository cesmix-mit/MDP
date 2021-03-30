#ifndef __COMMON_H__
#define __COMMON_H__

#define SDOT sdot_
#define SGEMV sgemv_
#define SGEMM sgemm_
#define SGETRF sgetrf_
#define SGETRI sgetri_

#define DDOT ddot_
#define DGEMV dgemv_
#define DGEMM dgemm_

#ifndef HAVE_CUDA    
#define DGETRF dgetrf_
#define DGETRI dgetri_
#endif

#ifdef USE_FLOAT
typedef float dstype;
#define cpuNRM2 SNRM2
#define cpuDOT SDOT
#define cpuAXPY SAXPY
#define cpuGEMV SGEMV
#define cpuGEMM SGEMM
#define cpuGETRF SGETRF
#define cpuGETRI SGETRI
#define cpuTRSM STRSM
#else
typedef double dstype; //  double is default precision 
#define cpuNRM2 DNRM2
#define cpuDOT DDOT
#define cpuAXPY DAXPY
#define cpuGEMV DGEMV
#define cpuGEMM DGEMM
#define cpuGETRF DGETRF
#define cpuGETRI DGETRI
#define cpuTRSM DTRSM
#endif

#ifdef USE_LONG
typedef long Int;
#else
typedef int Int; 
#endif

#ifndef HAVE_CUDA    
#define cublasHandle_t int
#define cudaEvent_t int
#endif

#define MKL_INT int

#define NEIGHMASK 0x3FFFFFFF
#define MAX_ATOM_TYPES      30
#define MAX_MOLECULE_TYPES  20
#define MAX_MOLECULE_SIZE   10

#define CPUFREE(x)                                                           \
{                                                                         \
    if (x != NULL) {                                                      \
        free(x);                                                          \
        x = NULL;                                                         \
    }                                                                     \
}

extern "C" {
    double DNRM2(Int*,double*,Int*);
    double DDOT(Int*,double*,Int*,double*,Int*);
    void DAXPY(Int*,double*,double*,Int*,double*,Int*);
    void DGEMV(char*,Int*,Int*,double*,double*,Int*,double*,Int*,double*,double*,Int*);  
    void DGEMM(char*,char*,Int*,Int*,Int*,double*,double*,Int*,
             double*,Int*,double*,double*,Int*);        
    void DGETRF(Int*,Int*,double*,Int*,Int*,Int*);
    void DGETRI(Int*,double*,Int*,Int*,double*,Int*,Int*);
    void DTRSM(char *, char*, char*, char *, Int *, Int *, double*, double*, Int*,
             double*, Int*);

    float SNRM2(Int*,float*,Int*);  
    float SDOT(Int*,float*,Int*,float*,Int*);
    void SAXPY(Int*,float*,float*,Int*,float*,Int*);
    void SGEMM(char*,char*,Int*,Int*,Int*,float*,float*,Int*,
             float*,Int*,float*,float*,Int*);  
    void SGEMV(char*,Int*,Int*,float*,float*,Int*,float*,Int*,float*,float*,Int*);      
    void SGETRF(Int*,Int*,float*,Int*,Int*,Int*);    
    void SGETRI(Int*,float*,Int*,Int*,float*,Int*,Int*);
    void STRSM(char *, char*, char*, char *, Int *, Int*, float*, float*, Int*,
             float*, Int*);    
}

// global variables for BLAS  
dstype one = 1.0;
dstype minusone = -1.0;
dstype zero = 0.0;
char chn = 'N';
char cht = 'T';
char chl = 'L';
char chu = 'U';
char chr = 'R';
Int inc1 = 1;

// global variables for CUBLAS  
// dstype *cublasOne;
// dstype *cublasMinusone;
// dstype *cublasZero;
dstype cublasOne[1] = {one};
dstype cublasMinusone[1] = {minusone};
dstype cublasZero[1] = {zero};

#ifdef HAVE_CUDA       
   #define CUDA_SYNC cudaDeviceSynchronize();  
#else 
   #define CUDA_SYNC
#endif                      

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
                
#ifdef HAVE_CUDA     

#ifdef USE_FLOAT
#define cublasNRM2 cublasSnorm2
#define cublasDOT cublasSdot
#define cublasAXPY cublasSaxpy
#define cublasGEMV cublasSgemv
#define cublasGEMM cublasSgemm
#define cublasGEMVBatched cublasSgemvBatched
#define cublasGEMMBatched cublasSgemmBatched
#define cublasGEMVStridedBatched cublasSgemvStridedBatched
#define cublasGEMMStridedBatched cublasSgemmStridedBatched
#define cublasGETRF cublasSgetrf
#define cublasGETRI cublasSgetri
#define cublasGETRFBatched cublasSgetrfBatched
#define cublasGETRIBatched cublasSgetriBatched
#define cublasTRSM cublasStrsm 
#else
#define cublasNRM2 cublasDnorm2
#define cublasDOT cublasDdot
#define cublasAXPY cublasDaxpy
#define cublasGEMV cublasDgemv
#define cublasGEMM cublasDgemm
#define cublasGEMVBatched cublasDgemvBatched
#define cublasGEMMBatched cublasDgemmBatched
#define cublasGEMVStridedBatched cublasDgemvStridedBatched
#define cublasGEMMStridedBatched cublasDgemmStridedBatched
#define cublasGETRF cublasDgetrf
#define cublasGETRI cublasDgetri
#define cublasGETRFBatched cublasDgetrfBatched
#define cublasGETRIBatched cublasDgetriBatched
#define cublasTRSM cublasDtrsm 
#endif

#define CHECK(call)                                                            \
{                                                                              \
    const cudaError_t error = call;                                            \
    if (error != cudaSuccess)                                                  \
    {                                                                          \
        fprintf(stderr, "Error: %s:%d, ", __FILE__, __LINE__);                 \
        fprintf(stderr, "code: %d, reason: %s\n", error,                       \
                cudaGetErrorString(error));                                    \
        exit(1);                                                               \
    }                                                                          \
}

#define CHECK_CUBLAS(call)                                                     \
{                                                                              \
    cublasStatus_t err;                                                        \
    if ((err = (call)) != CUBLAS_STATUS_SUCCESS)                               \
    {                                                                          \
        fprintf(stderr, "Got CUBLAS error %d at %s:%d\n", err, __FILE__,       \
                __LINE__);                                                     \
        exit(1);                                                               \
    }                                                                          \
}

#define GPUFREE(x)                                                       \
{                                                                         \
    if (x != NULL) {                                                      \
        cudaTemplateFree(x);                                              \
        x = NULL;                                                         \
    }                                                                     \
}

template <typename T> static void cudaTemplateMalloc(T **d_data, Int n)
{
    // allocate the memory on the GPU            
    CHECK( cudaMalloc( (void**)d_data, n * sizeof(T) ) );
}

template <typename T> static void cudaTemplateMallocManaged(T **d_data, Int n)
{
    // allocate unified memory 
    CHECK( cudaMallocManaged( (void**)d_data, n * sizeof(T) ) );        
}

template <typename T> static void cudaTemplateHostAlloc(T **h_data, Int n, unsigned int flags)
{
    // allocate zero-copy memory on host    
    CHECK(cudaHostAlloc((void **)h_data, n * sizeof(T), flags));                
}

template <typename T> static void cudaTemplateHostAllocMappedMemory(T **h_data, Int n)
{
    // allocate zero-copy memory on host    
    CHECK(cudaHostAlloc((void **)h_data, n * sizeof(T), cudaHostAllocMapped));                
}

template <typename T> static void cudaTemplateHostAllocPinnedMemory(T **h_data, Int n)
{
    // allocate pinned memory on host        
    CHECK(cudaHostAlloc((void **)h_data, n * sizeof(T), cudaHostAllocDefault));                
}

template <typename T> static void cudaTemplateFree(T *d_data)
{
    // free the memory on the GPU            
    CHECK( cudaFree( d_data ) );    
}

template <typename T> static void cudaCopytoDevice(T *d_data, T *h_data, Int n)
{
    // copy data from CPU to GPU
    CHECK( cudaMemcpy( d_data, h_data, n * sizeof(T), cudaMemcpyHostToDevice ) );    
}

template <typename T> static void cudaCopytoHost(T *h_data, T *d_data, Int n)
{
    // copy data from GPU to CPU
    CHECK( cudaMemcpy( h_data, d_data, n * sizeof(T), cudaMemcpyDeviceToHost ) );    
}

#endif

template <typename T> static void TemplateMalloc(T **data, Int n, Int backend)
{
    if (backend == 0)       // One thread CPU
        // allocate the memory on the CPU        
        *data = (T *) malloc(n*sizeof(T));      
    if (backend == 1)  // Open MP
        // allocate the memory on the CPU        
        *data = (T *) malloc(n*sizeof(T));      
#ifdef HAVE_CUDA            
    if (backend == 2)  // CUDA C                
        // allocate the memory on the GPU            
        CHECK( cudaMalloc( (void**)data, n * sizeof(T) ) );
#endif                  
}

struct commonstruct {     
    string filein;       // Name of binary file with input data
    string fileout;      // Name of binary file to write the solution    
            
    Int backend=1;   // 1: CPU; 2: CUDA GPU  
    Int mpiRank;  // MPI rank      
    Int mpiProcs; // number of MPI ranks
    
    Int dim;     // physical dimensions
    Int nconfigs; // number of configurations        
    Int natomtypes;     // number of atom types     
    Int nmoletypes;     // number of molecule types     
    Int K;       // order of radial basis functions
    Int L;       // order of spherical harmonics             
    Int Nub;     // number of unique bispectrum components
    Int Npower; // number of power spectrum components
    Int Nbispectrum; // number of bispectrum components
    Int Nbf=0;    // Nbf = Npower+Nbispectrum
    Int Ncoeff; // number of descriptors coefficients
    Int Ncg;// the total number of non-zero Clebsch-Gordan coefficients 
    Int descriptor;   // descriptor flag: 0 -> Spherical Harmonics Bessel
    Int spectrum;     // spectrum flag: 0-> power spectrum, 1-> bispectrum, 2-> both power and bispectrum 
    Int training;     // 0 -> no training, 1 -> Linear regression, 2 -> Gaussian process, 3 -> Neural net
    Int dftdata;      // 0 -> no data, 1 -> energies only, 2 -> forces only, 3 -> energies and forces
    Int runMD;        // 0 no MD simulation, 1 -> run MD simulation
    Int potential;    // 0 -> empirical potential, 1 -> empirical + LR, 2 -> empirical + GP, 3 -> empirical + NN
    Int cutofftype=0; // 0 -> single cut-off raidus for all atoms, 1 -> multiple cut-off radii for atom pairs  
    Int neighpair=0;  // 0 -> full neighbor list, 1 -> half neighbor list
    Int neighcell=0;  // 0 -> O(N^2) algorithm, 1 -> Cell-linked list algorithm to form neighbor list
    Int decomposition=0;// 0 -> force decomposition, 1 -> atom decomposition
    Int bondtype=0;   // 0 -> non-bonded interaction, 1 -> bonded interaction
    Int chemtype = 0;   // 0 -> single atom-type basis functions, 1 -> double atom-type basis functions 
    Int pairsymmetry;  // 1 -> V(r_ij) equal V(r_ji), 0 -> V(r_ij) not equal V(r_ji) 
    Int tripletsymmetry; // 1 -> V(r_ij, r_ik) equal V(r_ik, r_ij)
                         // 0 -> V(r_ij, r_ik) not equal V(r_ik, r_ij)
    Int energycal;    // turns energy calculation on or off
    Int forcecal;     // turns force calculation on or off
    Int stresscal;     // turns stress calculation on or off
        
    Int nflags;
    Int nsolversparam;   
    Int nsimulaparam;   
    Int neta;   
    Int nkappa;   
    Int nmuep;
    Int nmuml;
    Int nmu1a;
    Int nmu1b;
    Int nmu2a;
    Int nmu2b;
    Int nmu2c;
    Int nmu3a;
    Int nmu3b;
    Int nmu3c;
    Int nmu4a;
    Int nmu4b;    
    Int npot1a;
    Int npot1b;
    Int npot2a;
    Int npot2b;
    Int npot2c;
    Int npot3a;
    Int npot3b;
    Int npot3c;
    Int npot4a;
    Int npot4b;    
    Int natom1b;
    Int natom2b;
    Int natom2c;
    Int natom3b;
    Int natom3c;
    Int natom4b;    
    Int nmu[12];
            
    Int pnum;  // number of periodic images
    Int inum;  // number of atoms in the simulation box
    Int jnum;  // maximum number of neighbors
    Int inummax;  // maximum number of atoms in the simulation box
    Int gnum;  // number of ghost atoms
    Int anum;  // (inum+gnum) total number of atoms 
    Int anummax; // maximum  number of atoms allowed 
    Int cnum;  // number of cells
    Int nintmem;
    Int ntmpmem;
    Int nxij; //  number of xij = xj - xi  
    Int ne; 
    Int nt; 
    Int nx; 
    Int nv;
    Int nf;
    Int nq;
    Int ncq;
    Int nba; // number of blocks of atoms
    Int nab; // number of atoms per block
    Int nabmax; // maxmimum number of atoms per block
    Int ablks[17];
        
    Int ntimesteps; // number of time steps
    dstype dt;   // timestep size     
    dstype time; // current simulation time    
    dstype rcutml; // cutoff-radius for machine learning potential;
    dstype rbbvol; // volume of reference bounding box
    
    dstype boxoffset[3];
    Int pbc[3];
    
    Int *pot1a=NULL; 
    Int *pot1b=NULL;
    Int *pot2a=NULL; 
    Int *pot2b=NULL;
    Int *pot2c=NULL; 
    Int *pot3a=NULL;
    Int *pot3b=NULL; 
    Int *pot3c=NULL;
    Int *pot4a=NULL;
    Int *pot4b=NULL; 
    Int *atom1b=NULL;
    Int *atom2b=NULL;
    Int *atom2c=NULL; 
    Int *atom3b=NULL; 
    Int *atom3c=NULL;
    Int *atom4b=NULL; 
    
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
    }                         
};

struct appstruct {              
    Int *lsize=NULL;
    Int *nsize=NULL;  // input data size    
    Int *ndims=NULL;  // dimensions
    Int *flags=NULL;        // flag parameters            
    Int *bcs=NULL;           // boundary conditions
    Int *pbc=NULL;           // periodic boundary conditions        
        
    dstype *boxoffset=NULL;            
    Int *atomnumber=NULL;   // a list of atomic numbers for every atom type
    dstype *atommass=NULL; //  a list of atomic mass for every atom type
    dstype *atomcharge=NULL; //  a list of atomic charge for every atom type
    dstype *simulaparam=NULL; // simulation parameters   
    dstype *solversparam=NULL; // solvers parameters              
    dstype *physicsparam=NULL; // general physical parameters
    dstype *eta=NULL; // hyperparameters      
    Int *kappa=NULL; // integer parameters                  
                        
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
    Int *pot1a=NULL; 
    Int *pot1b=NULL;
    Int *pot2a=NULL; 
    Int *pot2b=NULL;
    Int *pot2c=NULL; 
    Int *pot3a=NULL;
    Int *pot3b=NULL; 
    Int *pot3c=NULL;
    Int *pot4a=NULL;
    Int *pot4b=NULL; 
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
    Int *atom1b=NULL;
    Int *atom2b=NULL;
    Int *atom2c=NULL; 
    Int *atom3b=NULL; 
    Int *atom3c=NULL;
    Int *atom4b=NULL; 
    
    // custom destructor
    void freememory(Int backend)
    {
       if (backend<=1) {
            CPUFREE(lsize);
            CPUFREE(nsize);
            CPUFREE(ndims);   
            CPUFREE(flags);    
            CPUFREE(bcs);    
            CPUFREE(pbc);    
            CPUFREE(boxoffset);
            CPUFREE(atomnumber);    
            CPUFREE(atommass);    
            CPUFREE(atomcharge);
            CPUFREE(simulaparam);
            CPUFREE(solversparam);
            CPUFREE(physicsparam);
            CPUFREE(eta);
            CPUFREE(kappa);
            CPUFREE(muml);
            CPUFREE(muep);
            CPUFREE(mu1a);
            CPUFREE(mu1b);
            CPUFREE(mu2a);
            CPUFREE(mu2b);
            CPUFREE(mu2c);
            CPUFREE(mu3a);
            CPUFREE(mu3b);
            CPUFREE(mu3c);
            CPUFREE(mu4a);
            CPUFREE(mu4b);            
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
            CPUFREE(rcutsqml);
            CPUFREE(rcutsq2a);
            CPUFREE(rcutsq2b);
            CPUFREE(rcutsq2c);
            CPUFREE(rcutsq3a);
            CPUFREE(rcutsq3b);
            CPUFREE(rcutsq3c);
            CPUFREE(rcutsq4a);
            CPUFREE(rcutsq4b);                        
            CPUFREE(rcutsq);            
            CPUFREE(atom1b);
            CPUFREE(atom2b);
            CPUFREE(atom2c);
            CPUFREE(atom3b);
            CPUFREE(atom3c);
            CPUFREE(atom4b);                        
        }            
#ifdef HAVE_CUDA      
       else {
            GPUFREE(lsize);
            GPUFREE(nsize);
            GPUFREE(ndims);   
            GPUFREE(flags);    
            GPUFREE(bcs);    
            GPUFREE(pbc);    
            GPUFREE(boxoffset);
            GPUFREE(atomnumber);    
            GPUFREE(atommass);    
            GPUFREE(atomcharge);
            GPUFREE(simulaparam);
            GPUFREE(solversparam);
            GPUFREE(physicsparam);
            GPUFREE(eta);
            GPUFREE(kappa);
            GPUFREE(muml);
            GPUFREE(muep);
            GPUFREE(mu1a);
            GPUFREE(mu1b);
            GPUFREE(mu2a);
            GPUFREE(mu2b);
            GPUFREE(mu2c);
            GPUFREE(mu3a);
            GPUFREE(mu3b);
            GPUFREE(mu3c);
            GPUFREE(mu4a);
            GPUFREE(mu4b);            
            GPUFREE(pot1a);
            GPUFREE(pot1b);
            GPUFREE(pot2a);
            GPUFREE(pot2b);
            GPUFREE(pot2c);
            GPUFREE(pot3a);
            GPUFREE(pot3b);
            GPUFREE(pot3c);
            GPUFREE(pot4a);
            GPUFREE(pot4b);            
            GPUFREE(rcutsqml);
            GPUFREE(rcutsq2a);
            GPUFREE(rcutsq2b);
            GPUFREE(rcutsq2c);
            GPUFREE(rcutsq3a);
            GPUFREE(rcutsq3b);
            GPUFREE(rcutsq3c);
            GPUFREE(rcutsq4a);
            GPUFREE(rcutsq4b);                        
            GPUFREE(rcutsq);            
            GPUFREE(atom1b);
            GPUFREE(atom2b);
            GPUFREE(atom2c);
            GPUFREE(atom3b);
            GPUFREE(atom3c);
            GPUFREE(atom4b);                       
       }
#endif       
    }
};

struct configstruct {    
    Int *lsize=NULL;
    Int *nsize=NULL;  // data size
    Int *ndims=NULL;  // dimensions    
    Int *natoms=NULL; // a list containing the number of atoms in each configuration 
    Int *natomssum=NULL; // a list containing the number of atoms in each configuration     
    Int *t=NULL;      // atomic types of atoms for all configurations 
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
    Int *cellnum=NULL;
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
            
    Int *atomtype=NULL;       // type of each atom i      
    Int *alist=NULL;
    Int *neighnum=NULL;  // numbers of neighbors for each atom i 
    Int *neighlist=NULL; // list of neighbors for each atom i    
            
    void freememory(Int backend)
    {
        if (backend<=1) {            
            CPUFREE(a); 
            CPUFREE(b); 
            CPUFREE(c); 
            CPUFREE(cellnum); 
            CPUFREE(cellsize); 
            CPUFREE(eta1); 
            CPUFREE(eta2); 
            CPUFREE(eta3);
            CPUFREE(boxvertices); 
            CPUFREE(bbvertices); 
            CPUFREE(refvertices); 
            CPUFREE(rbvertices); 
            CPUFREE(s2rmap); 
            CPUFREE(pimages); 
            CPUFREE(neighnum); 
            CPUFREE(neighlist);
            CPUFREE(alist);
            CPUFREE(atomtype); 
        }
#ifdef HAVE_CUDA                 
        else {         
            GPUFREE(a); 
            GPUFREE(b); 
            GPUFREE(c); 
            GPUFREE(cellnum); 
            GPUFREE(cellsize); 
            GPUFREE(eta1); 
            GPUFREE(eta2); 
            GPUFREE(eta3);
            GPUFREE(boxvertices); 
            GPUFREE(bbvertices); 
            GPUFREE(refvertices); 
            GPUFREE(rbvertices); 
            GPUFREE(s2rmap); 
            GPUFREE(pimages);             
            GPUFREE(neighnum); 
            GPUFREE(neighlist); 
            GPUFREE(alist); 
            GPUFREE(atomtype); 
        }
#endif               
    }
};

struct sysstruct {        
    dstype *x=NULL;    // cartesian coordinates of local atoms in processor mpiRank        
    dstype *v=NULL;    // velocities of local atoms in processor mpiRank        
    dstype *e=NULL;    // atomic energies of local atoms in processor mpiRank        
    dstype *f=NULL;    // atomic forces of local atoms in processor mpiRank        
    dstype *q=NULL;    // atomic charges of local atoms in processor mpiRank               
    dstype *s=NULL;    // stresses
    dstype *d=NULL;    // basis functions
    dstype *dd=NULL;   // derivatives of basis functions
    dstype *c=NULL;    // vector of coeffcients associated with the basis functions            
    dstype *A=NULL;    // Regression matrix A 
    dstype *b=NULL;    // Regression vector b
    
    void freememory(Int backend)
    {
       if (backend<=1) {
            CPUFREE(x); 
            CPUFREE(v); 
            CPUFREE(e); 
            CPUFREE(f); 
            CPUFREE(q); 
            CPUFREE(s);             
            CPUFREE(d);             
            CPUFREE(dd);             
            CPUFREE(c);       
            CPUFREE(A);       
            CPUFREE(b);       
        }            
#ifdef HAVE_CUDA      
       else {
            GPUFREE(x);
            GPUFREE(v);   
            GPUFREE(e); 
            GPUFREE(f);
            GPUFREE(q);
            GPUFREE(s);
            GPUFREE(d);
            GPUFREE(dd);             
            GPUFREE(c);             
            GPUFREE(A);       
            GPUFREE(b);       
       }
#endif       
    }                     
};

  
struct tempstruct {
    Int *intmem=NULL;     
    dstype *tmpmem=NULL;     
    dstype *buffrecv=NULL;
    dstype *buffsend=NULL;
    
    void freememory(Int backend)
    {
       if (backend<=1) {
            CPUFREE(intmem); 
            CPUFREE(tmpmem); 
            CPUFREE(buffrecv); 
            CPUFREE(buffsend); 
        }            
#ifdef HAVE_CUDA      
       else {
            GPUFREE(intmem); 
            GPUFREE(tmpmem);
            GPUFREE(buffrecv); 
            GPUFREE(buffsend); 
       }
#endif       
    }                     
};

struct shstruct {  
    Int L;  // the maximum degree of spherical harmonics 
    Int K;  // the number of zeros of spherical Bessel functions
    Int Nub;// number of non-zero unqiue bispectrum compoments of spherical harmonics 
    Int Ncg;// the total number of non-zero Clebsch-Gordan coefficients 
    Int npower;
    Int nbispectrum;
    Int nbasis;
    
    Int *indk=NULL;
    Int *indl=NULL;        
    Int *indm=NULL;       
    Int *rowm=NULL;    
  
    dstype *fac=NULL;
    dstype *cg=NULL;
    dstype *x0=NULL;
    dstype *f=NULL;
    dstype *P=NULL;
    dstype *tmp=NULL;
    dstype *df=NULL;
    dstype *dP=NULL;
    dstype *dtmp=NULL;   
    
    void freememory(Int backend)
    {
       if (backend<=1) {
            CPUFREE(indk); 
            CPUFREE(indl); 
            CPUFREE(indm); 
            CPUFREE(rowm); 
            CPUFREE(fac); 
            CPUFREE(cg); 
            CPUFREE(x0);             
            CPUFREE(f);             
            CPUFREE(P);             
            CPUFREE(tmp);             
            CPUFREE(df);             
            CPUFREE(dP);             
            CPUFREE(dtmp);             
        }            
#ifdef HAVE_CUDA      
       else {
            GPUFREE(indk); 
            GPUFREE(indl); 
            GPUFREE(indm); 
            GPUFREE(rowm); 
            GPUFREE(fac); 
            GPUFREE(cg); 
            GPUFREE(x0);             
            GPUFREE(f);             
            GPUFREE(P);             
            GPUFREE(tmp);             
            GPUFREE(df);             
            GPUFREE(dP);             
            GPUFREE(dtmp);             
       }
#endif       
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
