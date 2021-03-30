#include <string>
#include <fstream>
#include <sstream>
#include <iostream>
#include <math.h>
//#include <string.h>
//#include <stdlib.h>
//#include <sys/unistd.h>
//#include <random>

#ifdef _CUDA
#define HAVE_CUDA
#endif
        
#ifdef _MPI
#define HAVE_MPI
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

//using namespace std;
using std::cout;
using std::endl;
using std::string;
using std::ios;
using std::ofstream;
using std::ifstream;
using std::ostringstream;
using std::min;
using std::max;
using std::scientific;

#include "../Common/common.h"     // declaration of variables and structs
#include "../Common/core.h"       // interface to core libraries
#include "../Common/pblas.h"      // wrappers for blas libaries and MPI     

#include "../Configuration/Configuration.cpp" // read and preprocess input files 
#include "../Potentials/Potentials.cpp" // empirical potentials
#include "../Descriptors/Descriptors.cpp"   // nonparameteric potentials 
#include "../Calculation/Calculation.cpp"  // calculate energy and forces
#include "../Regression/Regression.cpp"  // perform linear/gaussian/dnn regressions 
#include "../Integration/Integration.cpp"  // perform MD simulations 

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
    
    backend=1; 
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
    cout<<"Available GPU Memory: "<<available<<" Total GPU Memory: "<<total<<endl;
#endif                           
    
    Int ngpus = 0;
    Int gpuid = 0;        
    if (backend==2)
        ngpus = 1;
    
    Int ci = 0;            
    dstype *x, *e, *f, *q, *param;
    
    // Read input files
    CCalculation CCal(filein, fileout, mpiprocs, mpirank, backend);       
    
    // set up configuration and allocate memory
    CCal.SetConfiguration(ci);
        
    x = &CCal.sys.x[0];
    e = &CCal.sys.e[0];
    f = &CCal.sys.f[0];
    q = &CCal.sys.q[0];
    param = &CCal.app.mu2a[0];
    Int nparam = CCal.common.nmu2a;
    
    // LJ potential
    for (ci=0; ci<CCal.common.nconfigs; ci++) {  // loop over each configuration                              
        // get atom positions for configuration ci       
        CCal.GetPositions(x, ci);   
        
        // get atom types for configuration ci
        CCal.GetAtomtypes(CCal.nb.atomtype, ci);           
        
        // initialize per-atom energy and forces
        ArraySetValue(e, 0.0, CCal.common.inum, CCal.common.backend);  
        ArraySetValue(f, 0.0, CCal.common.dim*CCal.common.inum, CCal.common.backend);  
        
        // form neighbor list
        CCal.NeighborList(x);
        
        // compute energy and forces for a nonbonded pair potential
        //CCal.NonbondedPairEnergyForce(e, f, x, q, param, nparam);
        CCal.EmpiricalPotentialEnergyForce(e, f, x, q, CCal.app.muep, CCal.common.nmu); 
        
        dstype energy = cpuArraySum(e, CCal.common.inum);
        cout<<"Configuration # "<<ci+1<<": "<<endl;
        cout<<"Potential energy: "<<energy<<endl;
        cout<<"Per-atom energies: "<<endl;
        printArray2D(e, 1, CCal.common.inum, backend);
        cout<<"Atom forces: "<<endl;
        printArray2D(f, CCal.common.dim, CCal.common.inum, backend);
    }
            
    // construct regression object
    CRegression CReg(CCal);
    
    // train machine learning potential using linear regression
    CReg.LinearRegression(CCal);
        
    // check linear regression errors
    for (int ci=0; ci<CCal.common.nconfigs; ci++) { // loop over each configuration             
        // get atom positions for configuration ci   
        CCal.GetPositions(x, ci);   

        // get atom types for configuration ci
        CCal.GetAtomtypes(CCal.nb.atomtype, ci);           

        // form neighbor list
        CCal.NeighborList(x);

        // Calculate energies and forces using ML potential
        ArraySetValue(e, 0.0, CCal.common.inum, CCal.common.backend);  
        ArraySetValue(f, 0.0, CCal.common.dim*CCal.common.inum, CCal.common.backend);  
        CCal.RadialSphericalHarmonicEnergyForce(e, f, x, CCal.sys.c, q, param, 0); 
        
        // check errors
        CCal.GetForces(x, ci);
        cout<<"Configuration # "<<ci+1<<": "<<cpuArraySum(e, CCal.common.inum)<<"  "<<CCal.config.e[ci]<<endl;
        printArray2D(f, CCal.common.dim, 10, CCal.common.backend);
        printArray2D(x, CCal.common.dim, 10, CCal.common.backend);                
    }            
    
    
    
//     if (CCal.common.K > 0) {        
//         dstype *d, *dd;        
//         TemplateMalloc(&d, CCal.common.Ncoeff, backend);         
//         TemplateMalloc(&dd, CCal.common.dim*CCal.common.inummax*CCal.common.Ncoeff, backend);         
//         for (ci=0; ci<CCal.common.nconfigs; ci++) {                  
//             // get atom positions for configuration ci   
//             CCal.GetPositions(x, ci);   
//             
//             // get atom types for configuration ci
//             CCal.GetAtomtypes(CCal.nb.atomtype, ci);           
//             
//             // form neighbor list
//             CCal.NeighborList(x);
// 
//             // spherical harmonic descriptors
//             ArraySetValue(d, 0.0, CCal.common.Ncoeff, CCal.common.backend);  
//             // derivatives of spherical harmonic descriptors
//             ArraySetValue(dd, 0.0, CCal.common.dim*CCal.common.inum*CCal.common.Ncoeff, CCal.common.backend);          
//             
//             // compute descriptors and their derivatives
//             CCal.RadialSphericalHarmonicDescriptors(d, dd, x, q, param, 0);
//             
//             cout<<"Descriptors for configuration # "<<ci+1<<": "<<endl;
//             printArray2D(d, 1, CCal.common.Ncoeff, CCal.common.backend);
//             //printArray2D(dd, CCal.common.dim*CCal.common.inum, CCal.common.Ncoeff, CCal.common.backend);            
//         }
//     }
    
//     int M = CCal.common.Ncoeff;
//     int dim = CCal.common.dim;
//     int K = CCal.common.nconfigs;
//     //int backend = CCal.common.backend;
//     
//     dstype *d, *dd, *Q;        
//     TemplateMalloc(&d, M, backend);         
//     TemplateMalloc(&dd, dim*CCal.common.inummax*M, backend);             
    
//     TemplateMalloc(&Q, K*M, backend);          
//     for (int ci=0; ci<CCal.common.nconfigs; ci++) {      
//         // get atom positions for configuration ci   
//         CCal.GetPositions(x, ci);   
// 
//         // get atom types for configuration ci
//         CCal.GetAtomtypes(CCal.nb.atomtype, ci);           
// 
//         // form neighbor list
//         CCal.NeighborList(x);
// 
//         // spherical harmonic descriptors
//         ArraySetValue(d, 0.0, M, backend);  
//         // derivatives of spherical harmonic descriptors
//         ArraySetValue(dd, 0.0, dim*CCal.common.inum*M, backend);          
//         
//         // compute descriptors and their derivatives
//         CCal.RadialSphericalHarmonicDescriptors(d, dd, x, q, param, 0);
//         
//         // insert d into Q at row ci
//         ArrayInsert(Q, d, K, M, 1, ci, ci+1, 0, M, 0, 1, backend);
//     }
    
    // total number of atoms for all configurations
//     int ki = 0, N = 0;  
//     for (int ci=0; ci<CCal.common.nconfigs; ci++)
//         N += CCal.config.natoms[ci];
//         
//     TemplateMalloc(&Q, (dim*N)*M, backend);          
//     for (int ci=0; ci<CCal.common.nconfigs; ci++) {      
//         // get atom positions for configuration ci   
//         CCal.GetPositions(x, ci);   
// 
//         // get atom types for configuration ci
//         CCal.GetAtomtypes(CCal.nb.atomtype, ci);           
// 
//         // form neighbor list
//         CCal.NeighborList(x);
// 
//         // number of atoms in configuration ci        
//         int inum = CCal.common.inum; 
//         // spherical harmonic descriptors
//         ArraySetValue(d, 0.0, M, backend);  
//         // derivatives of spherical harmonic descriptors
//         ArraySetValue(dd, 0.0, dim*inum*M, backend);          
//         
//         // compute descriptors and their derivatives
//         CCal.RadialSphericalHarmonicDescriptors(d, dd, x, q, param, 0);
//                 
//         // insert dd into Q at row ki
//         ArrayInsert(Q, dd, dim*N, M, 1, ki, ki+dim*inum, 0, M, 0, 1, backend);
//         ki += dim*inum;
//     }
    
//     TemplateMalloc(&Q, (K + dim*N)*M, backend);          
//     for (int ci=0; ci<CCal.common.nconfigs; ci++) {      
//         // get atom positions for configuration ci   
//         CCal.GetPositions(x, ci);   
// 
//         // get atom types for configuration ci
//         CCal.GetAtomtypes(CCal.nb.atomtype, ci);           
// 
//         // form neighbor list
//         CCal.NeighborList(x);
// 
//         // number of atoms in configuration ci        
//         int inum = CCal.common.inum; 
//         // spherical harmonic descriptors
//         ArraySetValue(d, 0.0, M, backend);  
//         // derivatives of spherical harmonic descriptors
//         ArraySetValue(dd, 0.0, dim*inum*M, backend);          
//         
//         // compute descriptors and their derivatives
//         CCal.RadialSphericalHarmonicDescriptors(d, dd, x, q, param, 0);
//                 
//         // insert d and dd into Q at row ki
//         ArrayInsert(Q, d, K, M, 1, ki, ki+1, 0, M, 0, 1, backend);
//         ArrayInsert(Q, dd, dim*N, M, 1, ki+1, ki+1+dim*inum, 0, M, 0, 1, backend);
//         ki += 1 + dim*inum;
//     }
    
#ifdef HAVE_MPI
  MPI_Finalize();
#endif 
    
    return 0;             
}
