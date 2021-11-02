/***************************************************************************                               
                    Molecular Dynamics Potentials (MDP)
                           CESMIX-MIT Project  
 
 Contributing authors: Ngoc-Cuong Nguyen (cuongng@mit.edu, exapde@gmail.com)
 ***************************************************************************/

#include <string>
#include <fstream>
#include <sstream>
#include <iostream>
#include <math.h>
#include <sys/time.h>
//#include <string.h>
//#include <stdlib.h>
//#include <sys/unistd.h>
//#include <random>

#ifdef _DEBUG
#define HAVE_DEBUG
#endif

#ifdef USE_CUDA
#define HAVE_CUDA
#endif

#ifdef _CUDA
#define HAVE_CUDA
#define USE_CUDA
#endif

#ifdef _ENZYME
#define HAVE_ENZYME
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

//#ifdef _TIMING
#include <chrono>
//#endif

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

#include "../Common/mdpcommon.h"     // declaration of variables and structs
#include "../Common/mdpcore.h"       // interface to core libraries
#include "../Common/pblas.h"      // wrappers for blas libaries and MPI     

#include "../Configuration/Configuration.cpp" // read and preprocess input files 
#include "../Potentials/cpuPotentials.cpp" // empirical potentials

#ifdef HAVE_CUDA    
#include "../Configuration/gpuDeviceInfo.cpp" // read and preprocess input files 
#include "../Potentials/gpuPotentials.cpp" // empirical potentials
#endif

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
            printf("MDP runs with %d MPI processes on GPU platform...\n", mpiprocs);
    }
    else {
        if (mpirank==0) 
            printf("MDP runs with %d MPI processes on CPU platform...\n", mpiprocs);
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
             
    // Read input files and set up Calculation class
    CCalculation CCal(filein, fileout, mpiprocs, mpirank, backend);       
    
    // set up configuration and allocate memory
    CCal.SetConfiguration(0);
    
    if (CCal.common.runMD) { 
        // construct integration object
        CIntegration CInt(CCal);

        // perform MD simulation using velocity verlet algorithm
        CInt.VelocityVerlet(CCal);
    } else if (CCal.common.training>0) {
        // construct regression object
        CRegression CReg(CCal);
        
        CReg.PotentialDescriptors(CCal);
//         // fit the potential using linear regression 
//         CReg.LinearRegression(CCal);    
// 
//         do {
//            cout <<"Potential fitting is done! Press RETURN key to continue potential validation.";
//         } while (std::cin.get() != '\n');
// 
//         // Validate linear regression potential
//         CReg.ValidateLinearRegression(CCal);
    }        
            
#ifdef HAVE_MPI
  MPI_Finalize();
#endif 
    
    return 0;             
}
























//     // check linear regression errors
//     for (int i=0; i<CCal.common.validatenum; i++) { // loop over each configuration             
//         int ci = CCal.common.validatelist[i]; // configuration ci
//         int N = CCal.common.dim*CCal.common.inum;
//         
//         // get atom positions for configuration ci   
//         CCal.GetPositions(x, ci);   
// 
//         // get atom types for configuration ci
//         CCal.GetAtomtypes(CCal.nb.atomtype, ci);           
// 
//         // form neighbor list
//         CCal.NeighborList(x);
// 
//         // Calculate energies and forces using ML potential
//         ArraySetValue(e, 0.0, CCal.common.inum, CCal.common.backend);  
//         ArraySetValue(f, 0.0, N, CCal.common.backend);  
//         CCal.PotentialEnergyForce(e, f, x, CCal.sys.c, q, CCal.app.muep, CCal.common.nmu); 
//                 
//         // check errors
//         CCal.GetForces(x, ci);
//         cout<<"Configuration # "<<ci+1<<": "<<cpuArraySum(e, CCal.common.inum)<<"  "<<CCal.config.e[ci]<<endl;
//         //printArray2D(f, CCal.common.dim, 10, CCal.common.backend);
//         //printArray2D(x, CCal.common.dim, 10, CCal.common.backend);                
//         dstype energyerror = fabs((cpuArraySum(e, CCal.common.inum)-CCal.config.e[ci])/CCal.config.e[ci]);
//         cout<<"Relative error in energy : "<<energyerror<<endl;        
//         dstype normf = PNORM(CCal.common.cublasHandle, N, x, CCal.common.backend);
//         ArrayAXPBY(f, f, x, one, minusone, N, CCal.common.backend);    
//         dstype norme = PNORM(CCal.common.cublasHandle, N, f, CCal.common.backend);
//         dstype forceerror = norme/normf;
//         cout<<"Relative error in forces : "<<forceerror<<endl;        
//     }            
        
    
//     // LJ potential
//     for (ci=0; ci<CCal.common.nconfigs; ci++) {  // loop over each configuration                              
//         // get atom positions for configuration ci       
//         CCal.GetPositions(x, ci);   
//         
//         // get atom types for configuration ci
//         CCal.GetAtomtypes(CCal.nb.atomtype, ci);           
//         
//         // initialize per-atom energy and forces
//         ArraySetValue(e, 0.0, CCal.common.inum, CCal.common.backend);  
//         ArraySetValue(f, 0.0, CCal.common.dim*CCal.common.inum, CCal.common.backend);  
//         
//         // form neighbor list
//         CCal.NeighborList(x);
//         
//         // compute energy and forces for a nonbonded pair potential
//         //CCal.NonbondedPairEnergyForce(e, f, x, q, param, nparam);
//         CCal.EmpiricalPotentialEnergyForce(e, f, x, q, CCal.app.muep, CCal.common.nmu); 
//         
//         dstype energy = cpuArraySum(e, CCal.common.inum);
//         cout<<"Configuration # "<<ci+1<<": "<<endl;
//         cout<<"Potential energy: "<<energy<<endl;
//         cout<<"Per-atom energies: "<<endl;
//         printArray2D(e, 1, CCal.common.inum, backend);
//         cout<<"Atom forces: "<<endl;
//         printArray2D(f, CCal.common.dim, CCal.common.inum, backend);
//     }
