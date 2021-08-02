/***************************************************************************                               
                    Molecular Dynamics Potentials (MDP)
                           CESMIX-MIT Project  
 
 Contributing authors: Ngoc-Cuong Nguyen (cuongng@mit.edu, exapde@gmail.com)
 ***************************************************************************/

#ifndef __REGRESSION
#define __REGRESSION

#include "Regression.h"

CRegression::CRegression(CCalculation &CCal)
{
    int M = CCal.common.M;
    int backend = CCal.common.backend;    
    TemplateMalloc(&CCal.sys.A, M*M, backend);   
    TemplateMalloc(&CCal.sys.b, M, backend);  
}

void CRegression::LinearRegression(CCalculation &CCal)
{
    int M = CCal.common.M;
    int dim = CCal.common.dim;
    int nc = CCal.common.trainingnum;
    int backend = CCal.common.backend;
    int dftdata =  CCal.common.dftdata;
    int nparam = CCal.common.nmu[0];    
    struct timeval tv1, tv2;
    
    dstype *A, *b, *c, *x, *f, *d, *dd, *q, *param;
    //param = &CCal.app.muml[0];    
    param = &CCal.app.muep[0];
    A = &CCal.sys.A[0];
    b = &CCal.sys.b[0];
    c = &CCal.sys.c[0];
    x = &CCal.sys.x[0];
    //e = &CCal.sys.e[0];
    f = &CCal.sys.f[0];
    q = &CCal.sys.q[0];
    d = &CCal.sys.d[0];
    dd = &CCal.sys.dd[0];
    
    // intialize regression matrix and vector
    ArraySetValue(A, 0.0, M*M, backend);  
    ArraySetValue(b, 0.0, M, backend);  
            
    if (dftdata == 1) { // energies only
        for (int i=0; i<nc; i++) { // loop over each configuration     
            int ci = CCal.common.traininglist[i]; // configuration ci
            cout<<"Configuration (energies) # "<<ci+1<<": "<<endl;
            
            // get atom positions for configuration ci   
            CCal.GetPositions(x, ci); 
            
            // get atom types for configuration ci
            CCal.GetAtomtypes(CCal.nb.atomtype, ci);           

            gettimeofday(&tv1, NULL); 
            // form neighbor list
            CCal.NeighborList(x);

            // spherical harmonic descriptors
            ArraySetValue(d, 0.0, M, backend);  
            gettimeofday(&tv2, NULL);            
            printf("\nExecution time (in millisec) for constructing the neighbor list:  %g\n", 
                (double)(tv2.tv_usec-tv1.tv_usec)/1000 + 
                (double)(tv2.tv_sec -tv1.tv_sec )*1000);

            gettimeofday(&tv1, NULL); 
            // compute descriptors 
            //CCal.RadialSphericalHarmonicDescriptors(d, x, q, param, 0);
            CCal.PotentialDescriptors(d, x, q, param, CCal.common.nmu);
            gettimeofday(&tv2, NULL);          
            printf("\nExecution time (in millisec) for computing the potential descrtiptors:  %g\n", 
                (double)(tv2.tv_usec-tv1.tv_usec)/1000 + 
                (double)(tv2.tv_sec -tv1.tv_sec )*1000);
                        
            gettimeofday(&tv1, NULL); 
            // apply a weight to the descriptor vector             
            cpuArrayMultiplyScalar(d, CCal.config.we[ci], M);                                
            
            // form the regression vector b = b + we[ci]*e[ci]*d           
            cpuArrayAXPBY(b, b, d, 1.0, CCal.config.we[ci]*CCal.config.e[ci], M);                                
            
            // form the regression matrix A = A + d * d^T
            cpuKron(A, d, d, M, M);            
            
            gettimeofday(&tv2, NULL);            
            printf("\nExecution time (in millisec) for forming the linear system:  %g\n", 
                (double)(tv2.tv_usec-tv1.tv_usec)/1000 + 
                (double)(tv2.tv_sec -tv1.tv_sec )*1000);            
        }        
    }
    else if (dftdata == 2) { // forces only
        for (int i=0; i<nc; i++) { // loop over each configuration     
            int ci = CCal.common.traininglist[i]; // configuration ci      
            cout<<"Configuration (forces) # "<<ci+1<<": "<<endl;
            
            Int N = dim*CCal.common.inum;
           
            // get atom positions for configuration ci   
            CCal.GetPositions(x, ci);   

            // get atom types for configuration ci
            CCal.GetAtomtypes(CCal.nb.atomtype, ci);           

             // get DFT forces
            CCal.GetForces(f, ci);                       
            
            gettimeofday(&tv1, NULL); 
            // form neighbor list
            CCal.NeighborList(x);

            // spherical harmonic descriptors
            ArraySetValue(d, 0.0, M, backend);  
            ArraySetValue(dd, 0.0, N*M, backend);  
            gettimeofday(&tv2, NULL);   
            printf("\nExecution time (in millisec) for constructing the neighbor list:  %g\n", 
                (double)(tv2.tv_usec-tv1.tv_usec)/1000 + 
                (double)(tv2.tv_sec -tv1.tv_sec )*1000);
            
            gettimeofday(&tv1, NULL); 
            // compute descriptors and their derivatives
            //CCal.RadialSphericalHarmonicDescriptors(&d[0], &dd[0], x, q, param, 0);
            CCal.PotentialDescriptors(d, dd, x, q, param, CCal.common.nmu);
            gettimeofday(&tv2, NULL);          
            printf("\nExecution time (in millisec) for computing the potential descrtiptors:  %g\n", 
                (double)(tv2.tv_usec-tv1.tv_usec)/1000 + 
                (double)(tv2.tv_sec -tv1.tv_sec )*1000);
            
            gettimeofday(&tv1, NULL); 
            // apply a weight to the descriptors' derivatives             
            cpuArrayMultiplyScalar(dd, -CCal.config.wf[ci], N*M);
            
            // form the regression vector b = b + wf[ci]*dd^T*f                                
            PGEMTV(CCal.common.cublasHandle, N, M, &CCal.config.wf[ci], dd, N, f, inc1, &one, b, inc1, backend);    
            
            // form the regression matrix A = A + dd^T * dd
            PGEMTM(CCal.common.cublasHandle, M, M, N, &one, dd, N, dd, N, &one, A, M, backend);                             
            gettimeofday(&tv2, NULL);            
            printf("\nExecution time (in millisec) for forming the linear system:  %g\n", 
                (double)(tv2.tv_usec-tv1.tv_usec)/1000 + 
                (double)(tv2.tv_sec -tv1.tv_sec )*1000);                        
        }                
    }
    else if (dftdata == 3) { // enegies and forces
        for (int i=0; i<nc; i++) { // loop over each configuration     
            int ci = CCal.common.traininglist[i]; // configuration ci       
            cout<<"Configuration (energies + forces) # "<<ci+1<<": "<<endl;
            
            Int N = dim*CCal.common.inum;
            
            // get atom positions for configuration ci   
            CCal.GetPositions(x, ci);   

             // get DFT forces 
            CCal.GetForces(f, ci);
                                   
            // get atom types for configuration ci
            CCal.GetAtomtypes(CCal.nb.atomtype, ci);           
 
            gettimeofday(&tv1, NULL); 
            // form neighbor list
            CCal.NeighborList(x);
            
            // spherical harmonic descriptors
            ArraySetValue(d, 0.0, M, backend);  
            ArraySetValue(dd, 0.0, N*M, backend);  
            gettimeofday(&tv2, NULL);            
            printf("\nExecution time (in millisec) for constructing the neighbor list:  %g\n", 
                (double)(tv2.tv_usec-tv1.tv_usec)/1000 + 
                (double)(tv2.tv_sec -tv1.tv_sec )*1000);
            
            gettimeofday(&tv1, NULL); 
            // compute descriptors             
            //CCal.RadialSphericalHarmonicDescriptors(d, dd, x, q, param, 0);
            CCal.PotentialDescriptors(d, dd, x, q, param, CCal.common.nmu);
            gettimeofday(&tv2, NULL);            
            printf("\nExecution time (in millisec) for computing the potential descrtiptors:  %g\n", 
                (double)(tv2.tv_usec-tv1.tv_usec)/1000 + 
                (double)(tv2.tv_sec -tv1.tv_sec )*1000);
                        
//             string fn = (backend == 2) ? "dgpu.bin" : "dcpu.bin";
//             writearray2file(fn, d, M, backend); 
//             fn = (backend == 2) ? "ddgpu.bin" : "ddcpu.bin";
//             writearray2file(fn, dd, N*M, backend);
//             error("here");
            
            gettimeofday(&tv1, NULL); 
            // apply a weight to the descriptor vectors      
            ArrayMultiplyScalar(d, CCal.config.we[ci], M, backend);
            
            // form the regression vector b = b + we[ci]*e[ci]*d            
            ArrayAXPBY(b, b, d, 1.0, CCal.config.we[ci]*CCal.config.e[ci], M, backend);
            
            // form the regression matrix A = A + d * d^T
            Kron(A, d, d, M, M, backend);                        
                        
            // apply a weight to the descriptors' derivatives       
            ArrayMultiplyScalar(dd, -CCal.config.wf[ci], N*M, backend);            
                        
            // form the regression vector b = b + wf[ci]*dd^T*f                                    
            PGEMTV2(CCal.common.cublasHandle, N, M, &CCal.config.wf[ci], dd, N, f, inc1, &one, b, inc1, backend);    
            
            // form the regression matrix A = A + dd^T * dd
            PGEMTM(CCal.common.cublasHandle, M, M, N, &one, dd, N, dd, N, &one, A, M, backend);                                     
            gettimeofday(&tv2, NULL);            
            printf("\nExecution time (in millisec) for forming the linear system:  %g\n", 
                (double)(tv2.tv_usec-tv1.tv_usec)/1000 + 
                (double)(tv2.tv_sec -tv1.tv_sec )*1000);
           
//             printArray2D(b, 1, M, backend);
//             printArray2D(A, M, M, backend);
//             error("here");            
        }                
    }

    // solve the linear system A*c = b         
    dstype *work = &CCal.tmp.tmpmem[0];  
    Int *ipiv = &CCal.tmp.intmem[0];
    Inverse(CCal.common.cublasHandle, A, work, ipiv, M, 1, backend);     
    PGEMNV(CCal.common.cublasHandle, M, M, &one, A, M, b, inc1, &zero, c, inc1, backend);    
    
    // write coefficients to a binary file
    writearray2file("coefficients.bin", c, M, backend);        
    
#ifdef HAVE_DEBUG                      
    string fn = (backend == 2) ? "cgpu.bin" : "ccpu.bin";
    writearray2file(fn, c, M, backend); 
    fn = (backend == 2) ? "bgpu.bin" : "bcpu.bin";
    writearray2file(fn, b, M, backend); 
#endif                                
    
}

void CRegression::ValidateLinearRegression(CCalculation &CCal, dstype *x, dstype *e, dstype *f, dstype *q)
{
    string filename = "validation.bin";
    ofstream out(filename.c_str(), ios::out | ios::binary);    
    if (!out) {
        error("Unable to open file " + filename);
    }    
    writedouble(out, (double) CCal.common.validatenum);        
            
    // check linear regression errors
    for (int i=0; i<CCal.common.validatenum; i++) { // loop over each configuration             
        int ci = CCal.common.validatelist[i]; // configuration ci
        int N = CCal.common.dim*CCal.common.inum;
                        
        // get atom positions for configuration ci   
        CCal.GetPositions(x, ci);   
        
        // get atom types for configuration ci
        CCal.GetAtomtypes(CCal.nb.atomtype, ci);           

        // form neighbor list
        CCal.NeighborList(x);                
        
        // Calculate energies and forces using ML potential
        ArraySetValue(e, 0.0, CCal.common.inum, CCal.common.backend);  
        ArraySetValue(f, 0.0, N, CCal.common.backend);  
                
        CCal.PotentialEnergyForce(e, f, x, CCal.sys.c, q, CCal.app.muep, CCal.common.nmu);                 
        dstype energy = ArraySum(e, CCal.tmp.tmpmem, CCal.common.inum, CCal.common.backend);
        
        writedouble(out, (double) ci);    
        writedouble(out, (double) N);                   
        writearray(out, x, N, CCal.common.backend);        
        writearray(out, f, N, CCal.common.backend);
        writedouble(out, energy); 
        
        // check errors
        CCal.GetForces(x, ci);
        cout<<"Configuration # "<<ci+1<<": "<<energy<<"  "<<CCal.config.e[ci]<<endl;
        printArray2D(f, CCal.common.dim, 10, CCal.common.backend);
        printArray2D(x, CCal.common.dim, 10, CCal.common.backend);                
        dstype energyerror = fabs((energy-CCal.config.e[ci])/CCal.config.e[ci]);
        cout<<"Relative error in energy : "<<energyerror<<endl;        
        dstype normf = PNORM(CCal.common.cublasHandle, N, x, CCal.common.backend);
        ArrayAXPBY(f, f, x, one, minusone, N, CCal.common.backend);    
        dstype norme = PNORM(CCal.common.cublasHandle, N, f, CCal.common.backend);
        dstype forceerror = norme/normf;
        cout<<"Relative error in forces : "<<forceerror<<endl;        
    }                
    
    out.close();
}

void CRegression::GaussianRegression(CCalculation &CCal)
{
}

void CRegression::NeuralNetRegression(CCalculation &CCal)
{
}
        
#endif        



        // Create CCal object
        // Pass the CCal object to Julia    

        // In C++, write a function that returns the object of a class
        // In Julia, I would call that function to get the object

        // create a pointer in Julia that points to x
        // pass that pointer
        // in C++, write a dummy interface function
        // double* calculateforce(CCal, double *x)
        // {
               // get  the pointer to x from Julia
               // CCal.PotentialEnergyForce(e, f, x, CCal.sys.c, q, CCal.app.muep, CCal.common.nmu);                 
               // return f;
        // }
