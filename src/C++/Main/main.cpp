#include <string>
#include <fstream>
#include <sstream>
#include <iostream>
#include <math.h>
#include <string.h>
//#include <stdlib.h>
//#include <chrono>
#include <sys/unistd.h>


#ifdef _OPENMP
#define HAVE_OPENMP
#else
#define HAVE_ONETHREAD
#endif

#ifdef _CUDA
#define HAVE_CUDA
#endif

#ifdef _MPI
#define HAVE_MPI
#endif

#ifdef HAVE_OPENMP
#include <omp.h>
#endif

#ifdef HAVE_CUDA    
#include <cuda_runtime.h>
#include <cuda.h>
#include "cublas_v2.h"
#endif

#ifdef HAVE_MPI
#include <mpi.h>
#endif

#ifdef TIMING
#include <chrono>
#endif

using namespace std;

#include "../Common/common.h"     // declaration of variables and structs
#include "../Common/core.h"       // interface to core library
#include "../Common/pblas.h"      // wrappers for blas libaries and MPI     
// #include "../Common/application.h"// interface to application library

#include "../Core/commonCore.cpp" 
#include "../Core/opuArrayPermute.cpp" 
#include "../Core/opuArrayOperations.cpp" 
#include "../Core/opuElemFaceNode.cpp" 

#include "../Core/cpuCoordinateTransformations.cpp" 
#include "../Core/cpuPrecompute.cpp" 
#include "../Core/cpuSphericalHarmonics.cpp" 

#include "../Preprocessing/errormsg.cpp"
#include "../Preprocessing/ioutilities.cpp"
//#include "../preprocessing/readbinaryfiles.cpp"

// #include "../Core/cpuArrayPermute.cpp" 
// #include "../Core/cpuArrayOperations.cpp" 

// #include "../Discretization/discretization.cpp" // discretization class
// #include "../Preconditioning/preconditioner.cpp" // preconditioner class
// #include "../Solver/solver.cpp"                 // solver class
// #include "../Solution/solution.cpp"             // solution class

int main(int argc, char** argv) 
{   
    if( argc >= 3 ) {
    }
    else {
      printf("Usage: ./cppfile InputFile OutputFile\n");
      return 1;
    }                
    
    string filein  = string(argv[1]); // input files
    string fileout  = string(argv[2]); // output files           
    Int restart, mpiprocs, mpirank, shmrank, ncores, nthreads, backend;    
    //int compnodes = 50;
 
    restart = 0;
    if (argc>=4) {
        string mystr = string(argv[3]);
        restart = stoi(mystr);
    }             
    
#ifdef HAVE_MPI    
    // Initialize the MPI environment
    MPI_Init(NULL, NULL);

    // Get the number of processes    
    MPI_Comm_size(MPI_COMM_WORLD, &mpiprocs);

    // Get the rank of the process
    MPI_Comm_rank(MPI_COMM_WORLD, &mpirank);
   
    MPI_Comm shmcomm;
    MPI_Comm_split_type(MPI_COMM_WORLD, MPI_COMM_TYPE_SHARED, 0,
                    MPI_INFO_NULL, &shmcomm);
    
    MPI_Comm_rank(shmcomm, &shmrank);
#else        
    // get number of MPI processes and MPI ranks
    mpiprocs = 1;
    mpirank = 0;
    shmrank = 0;
#endif                
        
#ifdef HAVE_OPENMP    
    // set OpenMP threads
    ncores = omp_get_num_procs();
    nthreads = omp_get_max_threads();
    omp_set_num_threads(nthreads);            
#else
    ncores = 1; 
    nthreads = 1;    
#endif
    
#ifdef HAVE_OPENMP
    backend=1; // multi-thread
#else
    backend=0; // single thread
#endif
#ifdef HAVE_CUDA  // CUDA          
    backend=2;
#endif
          
    if (backend==2) {
        if (mpirank==0) 
            printf("Using %d processors to solve the problem on GPU platform...\n", mpiprocs);
    }
    else {
        if (mpirank==0) 
            printf("Using %d processors to solve the problem on CPU platform...\n", mpiprocs);
    }
    
#ifdef HAVE_CUDA            
    int device;    
    cudaSetDevice(shmrank); 
    //gpuDeviceInfo(shmrank);
    cudaGetDevice( &device );
    size_t available, total;
    cudaMemGetInfo(&available, &total);
    //cout<<"Available GPU Memory: "<<available<<" Total GPU Memory: "<<total<<endl;
#endif                           
    
    Int ngpus = 0;
    Int gpuid = 0;        
    if (backend==2)
        ngpus = 1;
       
    // Open file to read
    ifstream in(filein.c_str(), ios::in | ios::binary);    
    if (!in) 
        error("Unable to open file " + filein);
           
    double *sph;
    int n = 100;
    int dim = 3;
    readarray(in, &sph, n*dim);
    //printArray2D(sph, n, dim, 0);
    
    double *the = &sph[0];
    double *phi = &sph[n];
    double *r = &sph[2*n];
    
    double xcoord[n*dim];
    double *x = &xcoord[0];
    double *y = &xcoord[n];
    double *z = &xcoord[2*n];
        
    cpuArraySphere2Cart(x, y, z, the, phi, r, n);
    //printArray2D(xcoord, n, dim, 0);
    
    double fac[168];    
    for (int i=0; i<168; i++)
        fac[i] = factable[i];
    //print1darray(fac,30);
    
    int L = 3;
    int M = (L+1)*(L+2)/2;
    double Ylmr[n*M];
    double Ylmi[n*M];
    double P[M];
    double tmp[M];
    double C[M];
    
    cpuSphericalHarmonics(Ylmr, Ylmi, the, phi, P, tmp, fac, C, M_PI, L, n);
    printArray2D(Ylmr, n, M, 0);    
    printArray2D(Ylmi, n, M, 0);
    

#ifdef HAVE_MPI
  MPI_Finalize();
#endif 
    
    return 0;             
}
