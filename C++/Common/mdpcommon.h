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

void print1iarray(Int* a, Int m)
{    
    for (Int i=0; i<m; i++)
        cout << a[i] << "   ";
    cout << endl;
}

void print2iarray(Int* a, Int m, Int n)
{
    for (Int i=0; i<m; i++) {
        for (Int j=0; j<n; j++)
            cout << a[j*m+i] << "   ";
        cout << endl;
    }
    cout << endl;
}

void print3iarray(Int* a, Int m, Int n, Int p)
{    
    for (Int k=0; k<p; k++) {
        for (Int i=0; i<m; i++) {
            for (Int j=0; j<n; j++)
                cout << a[k*n*m+j*m+i] << "   ";
            cout << endl;
        }
        cout << endl;
    }
    cout << endl;
}


void print1darray(dstype* a, Int m)
{
    //cout.precision(4);
    for (Int i=0; i<m; i++)
        cout << scientific << a[i] << "   ";
    cout << endl;
}

void printArray2D(Int* a, Int m, Int n, Int backend)
{
    if (backend==2) {
#ifdef HAVE_CUDA
        Int N = m*n;
        Int *b = (Int*) malloc (sizeof (Int)*N);
        cudaMemcpy(b, a, N*sizeof(Int), cudaMemcpyDeviceToHost);
        print2iarray(b, m, n);
        free(b);
#endif
    }
    else
        print2iarray(a, m, n);
}

void printArray3D(Int* a, Int m, Int n, Int p, Int backend)
{
    if (backend==2) {
#ifdef HAVE_CUDA
        Int N = m*n*p;
        Int *b = (Int*) malloc (sizeof (Int)*N);
        cudaMemcpy(b, a, N*sizeof(Int), cudaMemcpyDeviceToHost);
        print3iarray(b, m, n, p);
        free(b);
#endif
    }
    else
        print3iarray(a, m, n, p);
}

void print2darray(dstype* a, Int m, Int n)
{
    //cout.precision(4);
    for (Int i=0; i<m; i++) {
        for (Int j=0; j<n; j++)
            cout << scientific << a[j*m+i] << "   ";
        cout << endl;
    }
    cout << endl;
}

void print3darray(dstype* a, Int m, Int n, Int p)
{
    //cout.precision(8);
    for (Int k=0; k<p; k++) {
        for (Int i=0; i<m; i++) {
            for (Int j=0; j<n; j++)
                cout << scientific << a[k*n*m+j*m+i] << "   ";
            cout << endl;
        }
        cout << endl;
    }
    cout << endl;
}

void printArray2D(dstype* a, Int m, Int n, Int backend)
{
    if (backend==2) {
#ifdef  HAVE_CUDA        
        Int N = m*n;
        dstype *b = (dstype*) malloc (sizeof (dstype)*N);
        cudaMemcpy(b, a, N*sizeof(dstype), cudaMemcpyDeviceToHost);    
        print2darray(b, m, n);
        free(b);
#endif        
    }
    else
        print2darray(a, m, n);
}

void printArray3D(dstype* a, Int m, Int n, Int p, Int backend)
{
    if (backend==2) {
#ifdef  HAVE_CUDA        
        Int N = m*n*p;
        dstype *b = (dstype*) malloc (sizeof (dstype)*N);
        cudaMemcpy(b, a, N*sizeof(dstype), cudaMemcpyDeviceToHost);    
        print3darray(b, m, n, p);
        free(b);
#endif        
    }
    else
        print3darray(a, m, n, p);
}

struct latticestruct {     
    dstype *atombasis=NULL;//[12*3]; // fractional coords of each basis atom within unit cell (0 <= coord < 1)
    dstype *primitive=NULL;//[9]; // lattice <-> box transform matrices
    dstype *primitinv=NULL;//[9];
    dstype *rotaterow=NULL;//[9];
    dstype *rotatecol=NULL;//[9];
    dstype *spacing=NULL;//[3]; // lattice scale factors in 3 dims
    dstype *origin=NULL;//[3]; // lattice origin
    dstype *sublo=NULL;//[3];  // sub-box bounds in lattice space on this proc
    dstype *subhi=NULL;//[3];  // sub-box bounds in lattice space on this proc
    dstype *a1=NULL;//[3]; // edge vectors of unit cell  
    dstype *a2=NULL;//[3]; // edge vectors of unit cell  
    dstype *a3=NULL;//[3]; // edge vectors of unit cell  
    dstype scale;
    int *atomtype;  // type of basis atoms
    int *orientx=NULL;//[3]; // lattice orientation vecs
    int *orienty=NULL;//[3]; // orientx = what lattice dir lies
    int *orientz=NULL;//[3]; //           along x dim in box
    int style;  // NONE,SC,FCC,etc
    int spaceflag;    
    int nbasis;                             // # of basis atoms in unit cell    
    int ilo, ihi, jlo, jhi, klo, khi;       // lattice bounds for sub-box on this proc
    int natom;                              // # of atoms in the sub-box on this proc

    void latticebounds()
    {
      ilo = static_cast<int> (sublo[0]) - 1;
      jlo = static_cast<int> (sublo[1]) - 1;
      klo = static_cast<int> (sublo[2]) - 1;
      ihi = static_cast<int> (subhi[0]) + 1;
      jhi = static_cast<int> (subhi[1]) + 1;
      khi = static_cast<int> (subhi[2]) + 1;

      if (sublo[0] < 0.0) ilo--;
      if (sublo[1] < 0.0) jlo--;
      if (sublo[2] < 0.0) klo--;
      
      int i,j,k,m;
      natom = 0;
      for (k = klo; k <= khi; k++) 
        for (j = jlo; j <= jhi; j++) 
          for (i = ilo; i <= ihi; i++) 
            for (m = 0; m < nbasis; m++) 
                natom += 1;            
    }
    
    void printout(int backend)
    {
        printf("style, spaceflag, nbasis, natom, ilo, ihi, jlo, jhi, klo, khi, scale\n");
        printf("%i %i %i %i %i %i %i %i %i %i %g\n", style, spaceflag, nbasis, natom, ilo, ihi, jlo, jhi, klo, khi, scale);        
        printf("origin: "); printArray2D(origin, 1, 3, backend);
        printf("spacing: "); printArray2D(spacing, 1, 3, backend);
        printf("orientx: "); printArray2D(orientx, 1, 3, backend);
        printf("orienty: "); printArray2D(orienty, 1, 3, backend);
        printf("orientz: "); printArray2D(orientz, 1, 3, backend);
        printf("a1: "); printArray2D(a1, 1, 3, backend);
        printf("a2: "); printArray2D(a2, 1, 3, backend);
        printf("a3: "); printArray2D(a3, 1, 3, backend);
        printf("sublo: "); printArray2D(sublo, 1, 3, backend);
        printf("subhi: "); printArray2D(subhi, 1, 3, backend);                
        printf("type: "); printArray2D(atomtype, 1, nbasis, backend);
        printf("basis: \n"); printArray2D(atombasis, 3, nbasis, backend);
        printf("primitive: \n"); printArray2D(primitive, 3, 3, backend);
        printf("primitinv: \n"); printArray2D(primitinv, 3, 3, backend);
        printf("rotaterow: \n"); printArray2D(rotaterow, 3, 3, backend);
        printf("rotatecol: \n"); printArray2D(rotatecol, 3, 3, backend);
    }    
        
    void allocatememory(int backend)
    {
        TemplateMalloc(&atombasis, 36, backend);
        TemplateMalloc(&atomtype, 12, backend);
        TemplateMalloc(&primitive, 9, backend);
        TemplateMalloc(&primitinv, 9, backend);   
        TemplateMalloc(&rotaterow, 9, backend);    
        TemplateMalloc(&rotatecol, 9, backend);    
        TemplateMalloc(&spacing, 3, backend);    
        TemplateMalloc(&origin, 3, backend);
        TemplateMalloc(&sublo, 3, backend);    
        TemplateMalloc(&subhi, 3, backend);    
        TemplateMalloc(&a1, 3, backend);
        TemplateMalloc(&a2, 3, backend);
        TemplateMalloc(&a3, 3, backend);
        TemplateMalloc(&orientx, 3, backend);
        TemplateMalloc(&orienty, 3, backend);
        TemplateMalloc(&orientz, 3, backend);
    }    
    
    void freememory(int backend)
    {
        TemplateFree(atombasis, backend);
        TemplateFree(atomtype, backend);
        TemplateFree(primitive, backend);
        TemplateFree(primitinv, backend);   
        TemplateFree(rotaterow, backend);    
        TemplateFree(rotatecol, backend);    
        TemplateFree(spacing, backend);    
        TemplateFree(origin, backend);
        TemplateFree(sublo, backend);    
        TemplateFree(subhi, backend);    
        TemplateFree(a1, backend);
        TemplateFree(a2, backend);
        TemplateFree(a3, backend);
        TemplateFree(orientx, backend);
        TemplateFree(orienty, backend);
        TemplateFree(orientz, backend);
    }    
};

void copylattice(latticestruct &dlat, latticestruct &hlat, int backend)
{
    dlat.style = hlat.style;  
    dlat.spaceflag = hlat.spaceflag;    
    dlat.nbasis = hlat.nbasis;                  
    dlat.ilo = hlat.ilo;                  
    dlat.ihi = hlat.ihi;                  
    dlat.jlo = hlat.jlo;                  
    dlat.jhi = hlat.jhi;                  
    dlat.klo = hlat.klo;                  
    dlat.khi = hlat.khi;                  
    TemplateCopytoDevice(dlat.atombasis, hlat.atombasis, 36, backend);
    TemplateCopytoDevice(dlat.primitive, hlat.primitive, 9, backend);
    TemplateCopytoDevice(dlat.primitive, hlat.primitinv, 9, backend);
    TemplateCopytoDevice(dlat.rotaterow, hlat.rotaterow, 9, backend);
    TemplateCopytoDevice(dlat.rotatecol, hlat.rotatecol, 9, backend);
    TemplateCopytoDevice(dlat.spacing, hlat.spacing, 3, backend);
    TemplateCopytoDevice(dlat.origin, hlat.origin, 3, backend);
    TemplateCopytoDevice(dlat.sublo, hlat.sublo, 3, backend);
    TemplateCopytoDevice(dlat.subhi, hlat.subhi, 3, backend);
    TemplateCopytoDevice(dlat.a1, hlat.a1, 3, backend);
    TemplateCopytoDevice(dlat.a2, hlat.a2, 3, backend);
    TemplateCopytoDevice(dlat.a3, hlat.a3, 3, backend);
    TemplateCopytoDevice(dlat.orientx, hlat.orientx, 3, backend);
    TemplateCopytoDevice(dlat.orienty, hlat.orienty, 3, backend);
    TemplateCopytoDevice(dlat.orientz, hlat.orientz, 3, backend);
}

struct domainstruct {     
    int box_exist;                          // 0 = not yet created, 1 = exists
    int dimension;                          // 2 = 2d, 3 = 3d
    int nonperiodic;                        // 0 = periodic in all 3 dims
                                            // 1 = periodic or fixed in all 6
                                            // 2 = shrink-wrap in any of 6
    int triclinic;    // 0 = orthog box, 1 = triclinic
    int tiltsmall;    // 1 if limit tilt, else 0  
    int box_change;           // 1 if any of next 3 flags are set, else 0
    int box_change_size;      // 1 if box size changes, 0 if not
    int box_change_shape;     // 1 if box shape changes, 0 if not
    int box_change_domain;    // 1 if proc sub-domains change, 0 if not
    int deform_flag;        // 1 if fix deform exist, else 0
    int deform_vremap;      // 1 if fix deform remaps v, else 0
    int deform_groupbit;    // atom group to perform v remap for    
    int *boundary=NULL;//[3*2];    // settings for 6 boundaries
                         // 0 = periodic
                         // 1 = fixed non-periodic
                         // 2 = shrink-wrap non-periodic
                         // 3 = shrink-wrap non-per w/ min
    int *pbc=NULL;//[3];          // xyz periodicity as array
    
    dstype *face=NULL;//[6*3];   // unit normals of 6 box faces
    dstype *corners=NULL;//[8*3];                     // 8 corner points
    dstype *faces=NULL;//[6*4*3];   // 4 corner pts of 6 prism faces   
    dstype *boxhi=NULL;//[3];   // orthogonal box global bounds
    dstype *boxlo=NULL;//[3];   // orthogonal box global bounds
    dstype *boxtilt=NULL;//[3]; // 3 tilt factors
    dstype *boxhi_lamda=NULL;//[3]; // lamda box = (0,1)
    dstype *boxlo_lamda=NULL;//[3]; // lamda box = (0,1)
    dstype *boxhi_bound=NULL;//[3];  // bounding box of tilted domain
    dstype *boxlo_bound=NULL;//[3];  // bounding box of tilted domain
    dstype *subhi=NULL;//[3]; // sub-box bounds on this proc
    dstype *sublo=NULL;//[3]; // sub-box bounds on this proc
    dstype *ssubhi=NULL;//[3]; // shifted sub-box on this proc
    dstype *ssublo=NULL;//[3]; // shifted sub-box on this proc
    dstype *bsubhi=NULL;//[3]; // bounding for sub-box on this proc
    dstype *bsublo=NULL;//[3]; // bounding for sub-box on this proc
    dstype *subhi_lamda=NULL;//[3];  // bounds of subbox in lamda
    dstype *sublo_lamda=NULL;//[3];  // bounds of subbox in lamda
    dstype *h=NULL;//[6]               // shape matrix in Voigt ordering
                                      // Voigt = xx,yy,zz,yz,xz,xy
    dstype *h_inv=NULL;//[6];  // inverse of h
    dstype *h_rate=NULL;//[6]; // rate of box size/shape change

    void printout(int backend)
    {
        printf("boxlo: "); printArray2D(boxlo, 1, 3, backend);
        printf("boxhi: "); printArray2D(boxhi, 1, 3, backend);
        printf("boxtilt: "); printArray2D(boxtilt, 1, 3, backend);
        printf("sublo: "); printArray2D(sublo, 1, 3, backend);
        printf("subhi: "); printArray2D(subhi, 1, 3, backend);
        printf("sublo_lamda: "); printArray2D(sublo_lamda, 1, 3, backend);
        printf("subhi_lamda: "); printArray2D(subhi_lamda, 1, 3, backend);
        printf("bsublo: "); printArray2D(bsublo, 1, 3, backend);
        printf("bsubhi: "); printArray2D(bsubhi, 1, 3, backend);
        printf("ssublo: "); printArray2D(ssublo, 1, 3, backend);
        printf("ssubhi: "); printArray2D(ssubhi, 1, 3, backend);
        printf("h: "); printArray2D(h, 1, 6, backend);
        printf("h_inv: "); printArray2D(h_inv, 1, 6, backend);
    }    
    
    void allocatememory(int backend)
    {
        TemplateMalloc(&boundary, 6, backend);
        TemplateMalloc(&pbc, 3, backend);
        TemplateMalloc(&face, 18, backend);
        TemplateMalloc(&corners, 24, backend);   
        TemplateMalloc(&faces, 72, backend);   
        TemplateMalloc(&boxlo, 3, backend);    
        TemplateMalloc(&boxhi, 3, backend);    
        TemplateMalloc(&boxtilt, 3, backend);    
        TemplateMalloc(&boxhi_lamda, 3, backend);
        TemplateMalloc(&boxlo_lamda, 3, backend);    
        TemplateMalloc(&boxhi_bound, 3, backend);    
        TemplateMalloc(&boxlo_bound, 3, backend);
        TemplateMalloc(&subhi, 3, backend);
        TemplateMalloc(&sublo, 3, backend);
        TemplateMalloc(&ssubhi, 3, backend);
        TemplateMalloc(&ssublo, 3, backend);
        TemplateMalloc(&bsubhi, 3, backend);
        TemplateMalloc(&bsublo, 3, backend);
        TemplateMalloc(&subhi_lamda, 3, backend);
        TemplateMalloc(&sublo_lamda, 3, backend);
        TemplateMalloc(&h, 6, backend);
        TemplateMalloc(&h_inv, 6, backend);
        TemplateMalloc(&h_rate, 6, backend);
    }    
    
    void freememory(int backend)
    {
        TemplateFree(boundary, backend);
        TemplateFree(pbc, backend);
        TemplateFree(face, backend);   
        TemplateFree(corners, backend);   
        TemplateFree(faces, backend);   
        TemplateFree(boxlo, backend);    
        TemplateFree(boxhi, backend);    
        TemplateFree(boxtilt, backend);    
        TemplateFree(boxhi_lamda, backend);
        TemplateFree(boxlo_lamda, backend);    
        TemplateFree(boxhi_bound, backend);    
        TemplateFree(boxlo_bound, backend);
        TemplateFree(subhi, backend);
        TemplateFree(sublo, backend);
        TemplateFree(ssubhi, backend);
        TemplateFree(ssublo, backend);
        TemplateFree(bsubhi, backend);
        TemplateFree(bsublo, backend);
        TemplateFree(subhi_lamda, backend);
        TemplateFree(sublo_lamda, backend);
        TemplateFree(h, backend);
        TemplateFree(h_inv, backend);
        TemplateFree(h_rate, backend);
    }        
};

void copydomain(domainstruct &ddom, domainstruct &hdom, int backend)
{
    ddom.box_exist = hdom.box_exist;  
    ddom.dimension = hdom.dimension;    
    ddom.nonperiodic = hdom.nonperiodic;                  
    ddom.triclinic = hdom.triclinic;   
    ddom.tiltsmall = hdom.tiltsmall;   
    ddom.box_change = hdom.box_change;   
    ddom.box_change_size = hdom.box_change_size;   
    ddom.box_change_shape = hdom.box_change_shape;   
    ddom.box_change_domain = hdom.box_change_domain;   
    ddom.deform_flag = hdom.deform_flag;   
    ddom.deform_vremap = hdom.deform_vremap;   
    ddom.deform_groupbit = hdom.deform_groupbit;   
    
    TemplateCopytoDevice(ddom.boundary, hdom.boundary, 6, backend);
    TemplateCopytoDevice(ddom.pbc, hdom.pbc, 3, backend);
    TemplateCopytoDevice(ddom.face, hdom.face, 18, backend);
    TemplateCopytoDevice(ddom.corners, hdom.corners, 24, backend);
    TemplateCopytoDevice(ddom.faces, hdom.faces, 72, backend);
    TemplateCopytoDevice(ddom.boxlo, hdom.boxlo, 3, backend);
    TemplateCopytoDevice(ddom.boxhi, hdom.boxhi, 3, backend);
    TemplateCopytoDevice(ddom.boxtilt, hdom.boxtilt, 3, backend);
    TemplateCopytoDevice(ddom.boxlo_lamda, hdom.boxlo_lamda, 3, backend);
    TemplateCopytoDevice(ddom.boxhi_lamda, hdom.boxhi_lamda, 3, backend);
    TemplateCopytoDevice(ddom.boxlo_bound, hdom.boxlo_bound, 3, backend);
    TemplateCopytoDevice(ddom.boxhi_bound, hdom.boxhi_bound, 3, backend);        
    TemplateCopytoDevice(ddom.sublo, hdom.sublo, 3, backend);
    TemplateCopytoDevice(ddom.subhi, hdom.subhi, 3, backend);
    TemplateCopytoDevice(ddom.ssublo, hdom.ssublo, 3, backend);
    TemplateCopytoDevice(ddom.ssubhi, hdom.ssubhi, 3, backend);
    TemplateCopytoDevice(ddom.bsublo, hdom.bsublo, 3, backend);
    TemplateCopytoDevice(ddom.bsubhi, hdom.bsubhi, 3, backend);
    TemplateCopytoDevice(ddom.sublo_lamda, hdom.sublo_lamda, 3, backend);
    TemplateCopytoDevice(ddom.subhi_lamda, hdom.subhi_lamda, 3, backend);    
    TemplateCopytoDevice(ddom.h, hdom.h, 6, backend);
    TemplateCopytoDevice(ddom.h_inv, hdom.h_inv, 6, backend);
    TemplateCopytoDevice(ddom.h_rate, hdom.h_rate, 6, backend);
}

struct regionstruct {       
    int style;
    int triclinic; 
    int interior;                     // 1 for interior, 0 for exterior
    int scaleflag;                    // 1 for lattice, 0 for box  
    int bboxflag;                // 1 if bounding box is computable
    int varshape;                // 1 if region shape changes over time
    int dynamic;                 // 1 if position/orient changes over time
    int moveflag, rotateflag;    // 1 if position/orientation changes
    int openflag;                // 1 if any face is open
    
    dstype *r=NULL;                   // distance between particle & surf, r > 0.0
    dstype *delx=NULL;
    dstype *dely=NULL;
    dstype *delz=NULL;    // vector from surface pt to particle
    dstype *radius=NULL;              // curvature of region at contact point
    int *iwall=NULL;                  // unique id of wall for storing shear history    
    int *openfaces=NULL;//[6];           // flags for which faces are open
    int *tri=NULL;//[12*3];    // 3 corner pts of 12 triangles (2 per face)
    dstype theta; // orientation
    dstype rprev; // speed of time-dependent radius, if applicable
    
    dstype *face=NULL;//[6*3];   // unit normals of 6 prism faces
    dstype *corners=NULL;//[8*3];  // 8 corner pts of prism     
    dstype *faces=NULL;//[6*4*3];   // 4 corner pts of 6 prism faces   
    dstype *boxlo=NULL;//[3];
    dstype *boxhi=NULL;//[3];
    dstype *boxtilt=NULL;//[3];
    dstype *clo=NULL;//[3]; // opposite corners of prism
    dstype *chi=NULL;//[3]; // opposite corners of prism   
    dstype *scale=NULL;//[3];
    dstype *extent_lo=NULL;//[3];
    dstype *extent_hi=NULL;//[3];
    dstype *a=NULL;//[3];  // edge vectors of region
    dstype *b=NULL;//[3];  // edge vectors of region
    dstype *c=NULL;//[3];  // edge vectors of region
    dstype *h=NULL;//[6];
    dstype *h_inv=NULL;//[6];    
    dstype *dx=NULL;//[3];    // current displacement 
    dstype *v=NULL;//[3];                 // translational velocity
    dstype *rpoint=NULL;//[3];            // current origin of rotation axis
    dstype *omega=NULL;//[3];             // angular velocity    
    dstype *xcenter=NULL;//[3];    // translated/rotated center of cylinder/sphere (only used if varshape)
    dstype *prev=NULL;//[5];       // stores displacement (X3), angle and if
    
    void allocatememory(int backend)
    {
        TemplateMalloc(&face, 18, backend);
        TemplateMalloc(&corners, 24, backend);
        TemplateMalloc(&faces, 72, backend);   
        TemplateMalloc(&boxlo, 3, backend);    
        TemplateMalloc(&boxhi, 3, backend);    
        TemplateMalloc(&boxtilt, 3, backend);    
        TemplateMalloc(&clo, 3, backend);
        TemplateMalloc(&chi, 3, backend);    
        TemplateMalloc(&scale, 3, backend);    
        TemplateMalloc(&extent_lo, 3, backend);
        TemplateMalloc(&extent_hi, 3, backend);
        TemplateMalloc(&a, 3, backend);
        TemplateMalloc(&b, 3, backend);
        TemplateMalloc(&c, 3, backend);
        TemplateMalloc(&h, 6, backend);
        TemplateMalloc(&h_inv, 6, backend);
        TemplateMalloc(&dx, 3, backend);
        TemplateMalloc(&v, 3, backend);
        TemplateMalloc(&rpoint, 3, backend);
        TemplateMalloc(&omega, 3, backend);
        TemplateMalloc(&xcenter, 3, backend);
        TemplateMalloc(&prev, 5, backend);
        TemplateMalloc(&openfaces, 6, backend);
        TemplateMalloc(&tri, 36, backend);
    }    
                
    void freememory(int backend)
    {
        TemplateFree(face, backend);
        TemplateFree(corners, backend);
        TemplateFree(faces, backend);   
        TemplateFree(boxlo, backend);    
        TemplateFree(boxhi, backend);    
        TemplateFree(boxtilt, backend);    
        TemplateFree(clo, backend);
        TemplateFree(chi, backend);    
        TemplateFree(scale, backend);    
        TemplateFree(extent_lo, backend);
        TemplateFree(extent_hi, backend);
        TemplateFree(a, backend);
        TemplateFree(b, backend);
        TemplateFree(c, backend);
        TemplateFree(h, backend);
        TemplateFree(h_inv, backend);
        TemplateFree(dx, backend);
        TemplateFree(v, backend);
        TemplateFree(rpoint, backend);
        TemplateFree(omega, backend);
        TemplateFree(xcenter, backend);
        TemplateFree(prev, backend);
        TemplateFree(openfaces, backend);
        TemplateFree(tri, backend);
    }            
};

void copyregion(regionstruct &dreg, regionstruct &hreg, int backend)
{
    dreg.style = hreg.style;  
    dreg.interior = hreg.interior;    
    dreg.scaleflag = hreg.scaleflag;                  
    dreg.bboxflag = hreg.bboxflag;   
    dreg.varshape = hreg.varshape;   
    dreg.dynamic = hreg.dynamic;   
    dreg.moveflag = hreg.moveflag;   
    dreg.rotateflag = hreg.rotateflag;   
    dreg.openflag = hreg.openflag;   
        
    TemplateCopytoDevice(dreg.face, hreg.face, 18, backend);
    TemplateCopytoDevice(dreg.corners, hreg.corners, 24, backend);
    TemplateCopytoDevice(dreg.faces, hreg.faces, 72, backend);
    TemplateCopytoDevice(dreg.boxlo, hreg.boxlo, 3, backend);
    TemplateCopytoDevice(dreg.boxhi, hreg.boxhi, 3, backend);
    TemplateCopytoDevice(dreg.boxtilt, hreg.boxtilt, 3, backend);
    TemplateCopytoDevice(dreg.clo, hreg.clo, 3, backend);
    TemplateCopytoDevice(dreg.chi, hreg.chi, 3, backend);
    TemplateCopytoDevice(dreg.scale, hreg.scale, 3, backend);
    TemplateCopytoDevice(dreg.extent_lo, hreg.extent_lo, 3, backend);        
    TemplateCopytoDevice(dreg.extent_hi, hreg.extent_hi, 3, backend);        
    TemplateCopytoDevice(dreg.a, hreg.a, 3, backend);
    TemplateCopytoDevice(dreg.b, hreg.b, 3, backend);
    TemplateCopytoDevice(dreg.c, hreg.c, 3, backend);
    TemplateCopytoDevice(dreg.h, hreg.h, 6, backend);
    TemplateCopytoDevice(dreg.h_inv, hreg.h_inv, 6, backend);    
    TemplateCopytoDevice(dreg.dx, hreg.dx, 3, backend);
    TemplateCopytoDevice(dreg.v, hreg.v, 3, backend);    
    TemplateCopytoDevice(dreg.rpoint, hreg.rpoint, 3, backend);
    TemplateCopytoDevice(dreg.omega, hreg.omega, 3, backend);
    TemplateCopytoDevice(dreg.xcenter, hreg.xcenter, 3, backend);        
    TemplateCopytoDevice(dreg.prev, hreg.prev, 5, backend);
    TemplateCopytoDevice(dreg.openfaces, hreg.openfaces, 6, backend);
    TemplateCopytoDevice(dreg.tri, hreg.tri, 36, backend);    
}

struct commonstruct {     
    string filein;       // Name of binary file with input data
    string fileout;      // Name of binary file to write the solution    
            
    int backend=1;   // 1: CPU; 2: CUDA GPU  
    int mpiRank=0;  // MPI rank      
    int mpiProcs=1; // number of MPI ranks
    int readlattice=0;
    int readregion=0;
    int readdomain=0;
    
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
    int chemtype = 0;   // 0 -> single atom-type basis functions, 1 -> dstype atom-type basis functions 
    int pairsymmetry;  // 1 -> V(r_ij) equal V(r_ji), 0 -> V(r_ij) not equal V(r_ji) 
    int tripletsymmetry; // 1 -> V(r_ij, r_ik) equal V(r_ik, r_ij)
                         // 0 -> V(r_ij, r_ik) not equal V(r_ik, r_ij)
    int energycal;    // turns energy calculation on or off
    int forcecal;     // turns force calculation on or off
    int stresscal;     // turns stress calculation on or off
    int triclinic;     // types of the simulation box  
    int unitstyle=0;
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
    int eta_mass_flag=1;
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
    
    int *traininglist;
    int *validatelist;
    int trainingnum=0;
    int validatenum=0;
            
//     dstype pi = M_PI; // Pi  number
     
    latticestruct lat;
    regionstruct reg;
    domainstruct dom;
        
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
        lat.freememory(0);
        reg.freememory(0);
        dom.freememory(0);        
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
    dstype *boxtilt=NULL;
    dstype *boxhi_lamda=NULL;
    dstype *boxlo_lamda=NULL;
//     dstype *boxhi_bound=NULL;
//     dstype *boxlo_bound=NULL;
//     dstype *subhi=NULL;
//     dstype *sublo=NULL;
//     dstype *subhi_lamda=NULL;
//     dstype *sublo_lamda=NULL;
    dstype *h=NULL;
    dstype *h_inv=NULL;
    dstype *h_rate=NULL;
    dstype *box; // 9 parameters for the simulation box
    dstype *boxoffset=NULL;  
    
    dstype *atommass=NULL; //  a list of atomic mass for every atom type
    dstype *atomcharge=NULL; //  a list of atomic charge for every atom type
    dstype *simulaparam=NULL; // simulation parameters   
    dstype *solversparam=NULL; // solvers parameters              
    dstype *physicsparam=NULL; // general physical parameters
    dstype *nveparam=NULL; // NVE parameters     
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
    
    latticestruct lat;
    regionstruct reg;
    domainstruct dom;

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
        lat.freememory(backend);
        reg.freememory(backend);
        dom.freememory(backend);        
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

const dstype factable[] = {
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
