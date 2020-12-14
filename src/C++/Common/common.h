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
#ifdef HAVE_ONETHREAD         
    if (backend == 0)       // One thread CPU
        // allocate the memory on the CPU        
        *data = (T *) malloc(n*sizeof(T));      
#endif                   
#ifdef HAVE_OPENMP        
    if (backend == 1)  // Open MP
        // allocate the memory on the CPU        
        *data = (T *) malloc(n*sizeof(T));      
#endif                 
#ifdef HAVE_CUDA            
    if (backend == 2)  // CUDA C                
        // allocate the memory on the GPU            
        CHECK( cudaMalloc( (void**)data, n * sizeof(T) ) );
#endif                  
}

struct trainingstruct {     
    string trainingfile;    // name of a mother file for training data
    Int Ntrainconfig;   // number of training configurations
    Int *Ntrain;        // a list containing number of global atoms in each training configuration          
    void freememory()
    {
        CPUFREE(Ntrain);
    }                         
};

struct validationstruct {     
    string validationfile;  // name of a mother file for validation data
    Int Nvalidconfig;   // number of validation configurations
    Int *Nvalid;        // a list containing number of global atoms in each validation configuration    
    void freememory()
    {
        CPUFREE(Nvalid);
    }                         
};

struct commonstruct {     
    string filein;       // Name of binary file with input data
    string fileout;      // Name of binary file to write the solution    
            
    Int backend;   // 0: Serial; 1: OpenMP; 2: CUDA  
    Int cpuMemory; // 1: data in CPU memory; 0: data in GPU memory
    
    Int mpiRank;  // MPI rank      
    Int mpiProcs; // number of MPI ranks

    Int Natype;  // number of atom types
    Int Nmtype;  // number of molecule types 
    Int Za[MAX_ATOM_TYPES];  // a list of atomic numbers for all atom types
    Int Zm[MAX_MOLECULE_TYPES][MAX_MOLECULE_SIZE];  // a list of atomic numbers for all molecule types        
        
    Int Natomglobal; // number of global atoms in the configuration
    Int Nmoleglobal; // number of global molecules in the global configuration
    
    Int Natomlocal;  // number of local atoms in processor mpiRank
    Int Nmolelocal;  // number of local molecules in processor mpiRank
    Int Natomghost;  // number of ghost atoms in processor mpiRank 
    Int Nmoleghost;  // number of ghost molecules in processor mpiRank 
    Int inum;        // number of local atoms in the neighbor list
    Int gnum;        // number of ghost atoms in the neighbor list
        
    Int dim;     // physical dimensions
    Int Nbf;     // number of basis functions per atom type
    Int K;       // order of radial basis functions
    Int L;       // order of spherical harmonics                      
    Int nflags;  // number of flags
    Int nphysicsparam; // number of physical parameters
    Int nsolverparam;  // number of solver parameters
                
    Int newton_pair; // turns Newton’s third law on or off for pairwise interactions.    
    Int eflag;       // turns energy calculation on or off
    Int vflag;       // turns virial calculation on or off
    Int pairflag;    // 0 -> Nbf basis functions per atom type -> Nbf*Natype basis functions
                     // 1 -> Nbf basis functions per pair of two atom types -> Nbf*Natype*Natype basis functions 
    
    dstype lattice[6]; // lattice structure parameters
    dstype boxmin[3];  // minimum coordinates of the simulation box  
    dstype boxmax[3];  // maximum coordinates of the simulation box    
    dstype boxang[3];  // three angles of the simulation box
    dstype boxtensor[3*3];   // principal vectors of the simulation box    
    dstype boxvertices[3*8]; // vertices of the simulation box
    Int    boxfaces[4*6];    // faces of the simulation box: bottom, top, left, right, front, end 
    dstype boxvolume; // volume of the simulation box
    bool periodic[3]; // 1 -> periodic faces, 0-> non periodic faces         
    
    dstype localenergy;  // potential energy from local atoms in processor mpiRank
    dstype globalenergy; // potential energy from global atoms in the configuration
    dstype time; // current simulation time    
    dstype dt;   // timestep size     
    dstype Trange[2];  // temperature range
    dstype Prange[2];  // pressure range
    
    dstype rmincut;  // min cutoff for all atom types
    dstype rmaxcut;  // max cutoff for all atom types
    
    cudaEvent_t eventHandle;
    cublasHandle_t cublasHandle;
    
#ifdef  HAVE_MPI
    MPI_Request * requests;
    MPI_Status * statuses;
#endif
    
    void freememory()
    {
    }                         
};

struct appstruct {              
    Int *flags=NULL;   // flag parameters    
    dstype *fac;       // factorial    
    dstype *physicsparam=NULL; // physical parameters
    dstype *solversparam=NULL; // solvers parameters
    
    dstype *rminsq;  // square of the minimum cutoff radii between two different types of atom
    dstype *rmaxsq;  // square of the maximum cutoff radii between two different types of atom
    dstype *atomweights=NULL; // weight per atom type
    dstype *interatomweights=NULL; // weight for the interaction between two different types of atom
    dstype *moleculeweights=NULL; // weight per molecule type
    dstype *intermoleculeweights=NULL; // weight for the interaction between two different types of molecules
            
    // custom destructor
    void freememory(Int hostmemory)
    {
       if (hostmemory==1) {
            CPUFREE(flags);    
            CPUFREE(fac);            
            CPUFREE(physicsparam);
            CPUFREE(solversparam);
            CPUFREE(rminsq);
            CPUFREE(rmaxsq);
            CPUFREE(atomweights);
            CPUFREE(interatomweights);
            CPUFREE(moleculeweights);
            CPUFREE(intermoleculeweights);
        }            
#ifdef HAVE_CUDA      
       else {
            GPUFREE(flags);    
            GPUFREE(fac);            
            GPUFREE(physicsparam);
            GPUFREE(solversparam);
            GPUFREE(rminsq);
            GPUFREE(rmaxsq);
            GPUFREE(atomweights);
            GPUFREE(interatomweights);
            GPUFREE(moleculeweights);
            GPUFREE(intermoleculeweights);
       }
#endif       
    }
};

struct neighborstruct {  
  Int *ilist;          // local indices of I atoms
  Int *numneigh;       // number of J neighbors for each I atom  
  Int *ptrneigh;       // local indices of J neighbors    
};

struct sysstruct {        
    Int *type; // a vector of atom types for every local atom    
    dstype *x=NULL;    // cartesian coordinates of local atoms in processor mpiRank        
    dstype *v=NULL;    // velocities of local atoms in processor mpiRank        
    dstype *f=NULL;    // atomic forces of local atoms in processor mpiRank        
    dstype *q=NULL;    // atomic charges of local atoms in processor mpiRank           
    dstype *s=NULL;    // stresses
    dstype *d=NULL;    // basis functions
    dstype *c=NULL;    // vector of coeffcients associated with the basis functions            
    
    void freememory(Int hostmemory)
    {
       if (hostmemory==1) {
            CPUFREE(type); 
            CPUFREE(x); 
            CPUFREE(v); 
            CPUFREE(f); 
            CPUFREE(q); 
            CPUFREE(s);             
            CPUFREE(d);             
            CPUFREE(c);             
        }            
#ifdef HAVE_CUDA      
       else {
            GPUFREE(type); 
            GPUFREE(x);
            GPUFREE(v);   
            GPUFREE(f);
            GPUFREE(q);
            GPUFREE(s);
            GPUFREE(d);
            GPUFREE(c);                   
       }
#endif       
    }                     
};

  
struct tempstruct {
    dstype *tempn=NULL; 
    dstype *tempg=NULL;
    dstype *buffrecv=NULL;
    dstype *buffsend=NULL;
    
    void freememory(Int hostmemory)
    {
       if (hostmemory==1) {
            CPUFREE(tempn); 
            CPUFREE(tempg); 
            CPUFREE(buffrecv); 
            CPUFREE(buffsend); 
        }            
#ifdef HAVE_CUDA      
       else {
            GPUFREE(tempn);
            GPUFREE(tempg);
            GPUFREE(buffrecv); 
            GPUFREE(buffsend); 
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
