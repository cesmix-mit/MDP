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

            // form neighbor list
            CCal.NeighborList(x);

            // spherical harmonic descriptors
            ArraySetValue(d, 0.0, M, backend);  

            // compute descriptors 
            //CCal.RadialSphericalHarmonicDescriptors(d, x, q, param, 0);
            CCal.PotentialDescriptors(d, x, q, param, CCal.common.nmu);
            
            // apply a weight to the descriptor vector             
            cpuArrayMultiplyScalar(d, CCal.config.we[ci], M);                                
            
            // form the regression vector b = b + we[ci]*e[ci]*d           
            cpuArrayAXPBY(b, b, d, 1.0, CCal.config.we[ci]*CCal.config.e[ci], M);                                
            
            // form the regression matrix A = A + d * d^T
            cpuKron(A, d, d, M, M);                   
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

            // form neighbor list
            CCal.NeighborList(x);

            // spherical harmonic descriptors
            ArraySetValue(d, 0.0, M, backend);  
            ArraySetValue(dd, 0.0, N*M, backend);  
            
            // compute descriptors and their derivatives
            //CCal.RadialSphericalHarmonicDescriptors(&d[0], &dd[0], x, q, param, 0);
            CCal.PotentialDescriptors(d, dd, x, q, param, CCal.common.nmu);
            
            // apply a weight to the descriptors' derivatives             
            cpuArrayMultiplyScalar(dd, -CCal.config.wf[ci], N*M);
            
            // get DFT forces
            CCal.GetForces(f, ci);                       
            
            // form the regression vector b = b + wf[ci]*dd^T*f                                
            PGEMTV(CCal.common.cublasHandle, N, M, &CCal.config.wf[ci], dd, N, f, inc1, &one, b, inc1, backend);    
            
            // form the regression matrix A = A + dd^T * dd
            PGEMTM(CCal.common.cublasHandle, M, M, N, &one, dd, N, dd, N, &one, A, M, backend);                             
        }                
    }
    else if (dftdata == 3) { // enegies and forces
        for (int i=0; i<nc; i++) { // loop over each configuration     
            int ci = CCal.common.traininglist[i]; // configuration ci       
            cout<<"Configuration (energies + forces) # "<<ci+1<<": "<<endl;
            
            Int N = dim*CCal.common.inum;
            
            // get atom positions for configuration ci   
            CCal.GetPositions(x, ci);   

            // get atom types for configuration ci
            CCal.GetAtomtypes(CCal.nb.atomtype, ci);           

            // form neighbor list
            CCal.NeighborList(x);

            // spherical harmonic descriptors
            ArraySetValue(d, 0.0, M, backend);  
            ArraySetValue(dd, 0.0, N*M, backend);  
            
            // compute descriptors 
            //CCal.RadialSphericalHarmonicDescriptors(d, dd, x, q, param, 0);
            CCal.PotentialDescriptors(d, dd, x, q, param, CCal.common.nmu);
            
            // apply a weight to the descriptor vectors      
            cpuArrayMultiplyScalar(d, CCal.config.we[ci], M);
            
            // form the regression vector b = b + we[ci]*e[ci]*d            
            cpuArrayAXPBY(b, b, d, 1.0, CCal.config.we[ci]*CCal.config.e[ci], M);
            
            // form the regression matrix A = A + d * d^T
            cpuKron(A, d, d, M, M);                        
            
            // apply a weight to the descriptors' derivatives       
            CCal.config.wf[ci] = CCal.config.wf[ci]*10.0;
            cpuArrayMultiplyScalar(dd, -CCal.config.wf[ci], N*M);            
            
            // get DFT forces 
            CCal.GetForces(f, ci);
                                    
            // form the regression vector b = b + wf[ci]*dd^T*f                                    
            PGEMTV(CCal.common.cublasHandle, N, M, &CCal.config.wf[ci], dd, N, f, inc1, &one, b, inc1, backend);    
            
            // form the regression matrix A = A + dd^T * dd
            PGEMTM(CCal.common.cublasHandle, M, M, N, &one, dd, N, dd, N, &one, A, M, backend);             
        }                
    }

    // solve the linear system A*c = b         
    dstype *work = &CCal.tmp.tmpmem[0];  
    Int *ipiv = &CCal.tmp.intmem[0];
    Inverse(CCal.common.cublasHandle, A, work, ipiv, M, 1, backend);     
    PGEMNV(CCal.common.cublasHandle, M, M, &one, A, M, b, inc1, &zero, c, inc1, backend);    
    
    // write coefficients to a binary file
    writearray2file("coefficients.bin", c, M, backend);        
}

void CRegression::GaussianRegression(CCalculation &CCal)
{
}

void CRegression::NeuralNetRegression(CCalculation &CCal)
{
}
        
#endif        

