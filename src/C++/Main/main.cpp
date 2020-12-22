#include <string>
#include <fstream>
#include <sstream>
#include <iostream>
#include <math.h>
#include <string.h>
//#include <stdlib.h>
//#include <chrono>
#include <sys/unistd.h>
#include <random>


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

#include "../Preprocessing/errormsg.cpp"
#include "../Preprocessing/ioutilities.cpp"

#include "../Core/commonCore.cpp" 
#include "../Core/opuArrayPermute.cpp" 
#include "../Core/opuArrayOperations.cpp" 

#include "../Core/cpuCoordinateTransformations.cpp" 
#include "../Core/cpuPrecompute.cpp" 
#include "../Core/cpuSphericalHarmonics.cpp" 

//#include "../preprocessing/readbinaryfiles.cpp"

#include "../Descriptors/SphericalHarmonics.cpp"

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
           
    dstype *sph;
    int n = 100;
    int dim = 3;
    readarray(in, &sph, n*dim);
    in.close();
        
    dstype *the = &sph[0];
    dstype *phi = &sph[n];
    dstype *r = &sph[2*n];
    
    dstype xcoord[n*dim];
    dstype *x = &xcoord[0];
    dstype *y = &xcoord[n];
    dstype *z = &xcoord[2*n];
        
    cpuArraySphere2Cart(x, y, z, the, phi, r, n);
    
    int K = 5, L = 3;
    int L2 = (L+1)*(L+1); 
    int numneigh[2];
    numneigh[1] = n;
    numneigh[0] = 0;
    
    CSphericalHarmonics shobj(K, L, backend);
    
    dstype sr[n*K*L2], si[n*K*L2];
    dstype ar[K*L2], ai[K*L2];    
    dstype arx[n*K*L2], ary[n*K*L2], arz[n*K*L2];
    dstype aix[n*K*L2], aiy[n*K*L2], aiz[n*K*L2];    
    shobj.SphericalHarmonicsBessel(sr, si, x, y, z, n);       
    shobj.SphericalHarmonicsBesselDeriv(arx, aix, ary, aiy, arz, aiz, x, y, z, n);                    
    shobj.RadialSphericalHarmonicsSum(ar, ai, sr, si, numneigh, 1);      
    
    dstype p[(L+1)*K*(K+1)/2];
    shobj.RadialSphericalHarmonicsPower(p, ar, ai, 1);    
    
    dstype px[n*(L+1)*K*(K+1)/2], py[n*(L+1)*K*(K+1)/2], pz[n*(L+1)*K*(K+1)/2];
    shobj.RadialSphericalHarmonicsPowerDeriv(px, py, pz, ar, ai, 
        arx, aix, ary, aiy, arz, aiz, numneigh, 1);        
    
    int Nub = shobj.Nub;
    dstype b[Nub*K*(K+1)/2];
    shobj.RadialSphericalHarmonicsBispectrum(b, ar, ai, 1);        
            
    dstype bx[n*Nub*K*(K+1)/2], by[n*Nub*K*(K+1)/2], bz[n*Nub*K*(K+1)/2];
    shobj.RadialSphericalHarmonicsBispectrumDeriv(bx, by, bz, ar, ai, 
        arx, aix, ary, aiy, arz, aiz, numneigh, 1);        
    
    
    
    
    
//     string filename = "besselzeros.bin";
//     ifstream fin(filename.c_str(), ios::in | ios::binary);    
//     if (!fin) 
//         error("Unable to open file " + filename);
//     dstype *bzeros;
//     readarray(fin, &bzeros, 25*20);
//     fin.close();
//    
//     dstype fac[168];    
//     for (int i=0; i<168; i++)
//         fac[i] = factable[i];
//     //print1darray(fac,30);
//     
//     int L = 3;
//     int M = (L+1)*(L+2)/2;
//     int L2 = (L+1)*(L+1); 
//     dstype Ylmr[n*M];
//     dstype Ylmi[n*M];
//     dstype P[M];
//     dstype tmp[M];
//     dstype C[M];
//     dstype Sr[L2];
//     dstype Si[L2];
//         
//     cpuSphericalHarmonics(Ylmr, Ylmi, x, y, z, P, tmp, fac, M_PI, L, n);
// //     printArray2D(Ylmr, n, M, 0);    
// //     printArray2D(Ylmi, n, M, 0);
//     
//     cpuSphericalHarmonicsSum(Sr, Si, Ylmr, Ylmi, L, n);
//     print1darray(Sr,L2);
//     print1darray(Si,L2);    
//         
//     dstype Qr[L2];
//     dstype Qi[L2];
//     cpuSphericalHarmonicsSum(Qr, Qi, x, y, z, P, tmp, fac, M_PI, L, n);
//     print1darray(Qr,L2);
//     print1darray(Qi,L2);    
// 
//     dstype b[L2*(L+1)];    
//     cpuSphericalHarmonicsBispectrum(b, Sr, Si, fac, L);    
//     int Nub = cpuSphericalHarmonicsBispectrumIndex(b, L);
//     int indl[Nub*3];
//     cpuSphericalHarmonicsBispectrumIndex(indl, b, Nub, L);
// //     printArray3D(b, L+1, L+1, L+1, 0);    
// //     printArray2D(indl, Nub, 3, 0);    
//     
//     int Ncg = cgcoefficients(indl, Nub);
//     int indm[Ncg*3];
//     int rowm[Nub+1];
//     dstype cg[Ncg];
//     cgcoefficients(cg, indm, rowm, indl, fac, Ncg, Nub);    
// //     cout<<Ncg<<endl;
// //     printArray2D(rowm, 1, Nub+1, 0);    
// //     printArray2D(cg, 1, Ncg, 0);    
// //     printArray2D(indm, Ncg, 3, 0);    
//     
//     int K = 5;
//     dstype g[n*K*(L+1)];
//     dstype f[L+1];
//     dstype x0[(L+1)*K];
//     for (int k=0; k<K; k++)
//         for (int l=0; l<(L+1); l++)
//             x0[k*(L+1) + l] = bzeros[k*25 + l];
// //     //printArray2D(x0, L+1, K, 0);    
// //     cpuSphericalBessel(g, r, x0, f, L, K, n);
// //     printArray2D(g, n, (L+1)*K, 0);
//     cpuSphericalBessel(g, x, y, z, x0, f, L, K, n);
// //     printArray2D(g, n, (L+1)*K, 0);
//  
//     dstype ar[K*L2];
//     dstype ai[K*L2];
//     cpuRadialSphericalHarmonicsSum(ar, ai, Ylmr, Ylmi, g, L, K, n);
// //     printArray2D(ar, K, L2, 0);    
// //     printArray2D(ai, K, L2, 0);    
//         
//     dstype cr[K*L2];
//     dstype ci[K*L2];
//     cpuSphericalHarmonicsBesselSum(cr, ci, x, y, z, x0, P, tmp, f, fac, M_PI, L, K, n);
// //     printArray2D(cr, K, L2, 0);    
// //     printArray2D(ci, K, L2, 0);    
//     
//     dstype Srx[n*K*L2], Sry[n*K*L2], Srz[n*K*L2];
//     dstype Six[n*K*L2], Siy[n*K*L2], Siz[n*K*L2];;
//     dstype df[L+1];
//     dstype dP[M];
//     dstype dtmp[M];
//     cpuSphericalHarmonicsBesselDeriv(Srx, Six, Sry, Siy, Srz, Siz, x, y, z, 
//                 x0, P, tmp, f, dP, dtmp, df, fac, M_PI, L, K, n);
//     //printArray2D(Srx, n, K*L2, 0);    
//     
//     int indk[K*(K+1)];
//     cpuGetIndk(indk, K);
//     //printArray2D(indk, K*(K+1)/2, 2, 0);    
//     
//     dstype p[(L+1)*K*(K+1)/2];
//     cpuRadialSphericalHarmonicsPower(p, ar, ai, indk, L, K);
//     //printArray2D(p, (L+1), K*(K+1)/2, 0);    
//     
//     dstype px[n*(L+1)*K*(K+1)/2], py[n*(L+1)*K*(K+1)/2], pz[n*(L+1)*K*(K+1)/2];
//     cpuRadialSphericalHarmonicsPowerDeriv(px, py, pz, ar, ai, 
//         Srx, Six, Sry, Siy, Srz, Siz, indk, L, K, n);        
//     //printArray2D(px, n, (L+1)*K*(K+1)/2, 0);    
//     
//     dstype bi[Nub*K*(K+1)/2];
//     cpuRadialSphericalHarmonicsBispectrum(bi, ar, ai, cg, indk, indl,
//         indm, rowm, Nub, Ncg, K);
//     //printArray2D(bi, Nub, K*(K+1)/2, 0);        
//     
//     dstype bx[n*Nub*K*(K+1)/2], by[n*Nub*K*(K+1)/2], bz[n*Nub*K*(K+1)/2];
//     cpuRadialSphericalHarmonicsBispectrumDeriv(bx, by, bz, ar, ai, 
//         Srx, Six, Sry, Siy, Srz, Siz, cg, indk, indl, indm, rowm, Nub, Ncg, K, n);
//     //printArray2D(bx, n, Nub*K*(K+1)/2, 0);        
//     
//     CSphericalHarmonics shobj(K, L, 0);
// //     cout<<Nub<<" "<<Ncg<<endl;
// //     cout<<shobj.Nub<<" "<<shobj.Ncg<<endl;
// 
//     //printArray2D(shobj.sh.indk, 1, K*(K+1)/2, 0);    
//     //printArray2D(shobj.sh.indm, Ncg, 3, 0);    
//     
//     dstype arx[n*K*L2], ary[n*K*L2], arz[n*K*L2];
//     dstype aix[n*K*L2], aiy[n*K*L2], aiz[n*K*L2];;    
//     shobj.SphericalHarmonicsBesselDeriv(arx, aix, ary, aiy, arz, aiz, x, y, z, n);        
//     opuArrayAXPBY(Srx, Srx, arx, 1.0, -1.0, n*K*L2);
//     opuArrayAbs(Srx, Srx, n*K*L2);    
//     cout<<opuArrayMax(Srx, n*K*L2)<<endl;
//     
//     int numneigh[2];
//     numneigh[1] = n;
//     numneigh[0] = 0;
//     dstype sr[n*K*L2], si[n*K*L2];
//     shobj.SphericalHarmonicsBessel(sr, si, x, y, z, n);       
//     shobj.RadialSphericalHarmonicsSum(ar, ai, sr, si, numneigh, 1);      
//     dstype ptm[(L+1)*K*(K+1)/2];
//     shobj.RadialSphericalHarmonicsPower(ptm, ar, ai, 1);    
//     opuArrayAXPBY(ptm, ptm, p, 1.0, -1.0, (L+1)*K*(K+1)/2);
//     opuArrayAbs(ptm, ptm, (L+1)*K*(K+1)/2);    
//     cout<<opuArrayMax(ptm, (L+1)*K*(K+1)/2)<<endl;
//     
//     dstype ptx[n*(L+1)*K*(K+1)/2], pty[n*(L+1)*K*(K+1)/2], ptz[n*(L+1)*K*(K+1)/2];
//     shobj.RadialSphericalHarmonicsPowerDeriv(ptx, pty, ptz, ar, ai, 
//         arx, aix, ary, aiy, arz, aiz, numneigh, 1);        
//     opuArrayAXPBY(pty, pty, py, 1.0, -1.0, n*(L+1)*K*(K+1)/2);
//     opuArrayAbs(pty, pty, n*(L+1)*K*(K+1)/2);    
//     cout<<opuArrayMax(pty, n*(L+1)*K*(K+1)/2)<<endl;
//     
//     dstype bt[Nub*K*(K+1)/2];
//     shobj.RadialSphericalHarmonicsBispectrum(bt, ar, ai, 1);        
//     opuArrayAXPBY(bt, bt, bi, 1.0, -1.0, Nub*K*(K+1)/2);
//     opuArrayAbs(bt, bt, Nub*K*(K+1)/2);    
//     cout<<opuArrayMax(bt, Nub*K*(K+1)/2)<<endl;    
//         
//     dstype btx[n*Nub*K*(K+1)/2], bty[n*Nub*K*(K+1)/2], btz[n*Nub*K*(K+1)/2];
//     shobj.RadialSphericalHarmonicsBispectrumDeriv(btx, bty, btz, ar, ai, 
//         arx, aix, ary, aiy, arz, aiz, numneigh, 1);        
//     opuArrayAXPBY(btz, btz, bz, 1.0, -1.0, n*Nub*K*(K+1)/2);
//     opuArrayAbs(btz, btz, n*Nub*K*(K+1)/2);    
//     cout<<opuArrayMax(btz, n*Nub*K*(K+1)/2)<<endl;
        
#ifdef HAVE_MPI
  MPI_Finalize();
#endif 
    
    return 0;             
}