/***************************************************************************                               
                    Molecular Dynamics Potentials (MDP)
                           CESMIX-MIT Project  
 
 Contributing authors: Ngoc-Cuong Nguyen (cuongng@mit.edu, exapde@gmail.com)
 ***************************************************************************/

#ifndef __OMPFORCEDECOMPOSITION
#define __OMPFORCEDECOMPOSITION

template <typename T> void ompSingleDecomposition(T *e, T *ei, int *ai, int inum, int dim)
{    
    #pragma omp parallel for
    for (int ii=0; ii<inum; ii++) {  // for each atom pair ij in the simulation box     
        int i = ai[ii];       // atom i        
        int l = dim*ii;                
        e[i] += ei[ii];
    }
}
template void ompSingleDecomposition(double*, double*, int*, int, int);
template void ompSingleDecomposition(float*, float*, int*, int, int);

template <typename T> void ompFullNeighPairDecomposition(T *e, T *eij, int *ai, int ijnum, int dim)
{    
    #pragma omp parallel for
    for (int ii=0; ii<ijnum; ii++) {  // for each atom pair ij in the simulation box     
        int i = ai[ii];       // atom i
        int l = dim*ii;
        // use atomicAdd on gpu
        #pragma omp atomic
        e[i] += 0.5*eij[ii];
    }
}
template void ompFullNeighPairDecomposition(double*, double*, int*, int, int);
template void ompFullNeighPairDecomposition(float*, float*, int*, int, int);

template <typename T> void ompCenterAtomPairDecomposition(T *e, T *eij, int *ilist, int *anumsum, int inum, int dim)
{    
    #pragma omp parallel for
    for (int ii=0; ii<inum; ii++) {  // for each atom i in the simulation box     
        int i = ilist[ii];       // atom i        
        int start = anumsum[ii];   
        int m = anumsum[ii+1]-start;        // number of neighbors around i             
        int idim = dim*i;        
        for (int l=0; l<m ; l++)    // loop over each atom around atom i            
            e[i] += 0.5*eij[start + l];                                
    }
}
template void ompCenterAtomPairDecomposition(double*, double*, int*, int*, int, int);
template void ompCenterAtomPairDecomposition(float*, float*, int*, int*, int, int);

template <typename T> void ompHalfNeighPairDecomposition(T *e, T *eij, int *ai, int *aj, int ijnum, int dim)
{    
    #pragma omp parallel for
    for (int ii=0; ii<ijnum; ii++) {  // for each atom pair ij in the simulation box     
        int i = ai[ii];       // atom i
        int j = aj[ii];       // atom j
        int l = dim*ii;        
        // use atomicAdd on gpu
        #pragma omp atomic
        e[i] += 0.5*eij[ii];
        #pragma omp atomic
        e[j] += 0.5*eij[ii];
    }
}
template void ompHalfNeighPairDecomposition(double*, double*, int*, int*, int, int);
template void ompHalfNeighPairDecomposition(float*, float*, int*, int*, int, int);

template <typename T> void ompNeighborAtomPairDecomposition(T *e, T *eij, int *jlist, int *bnumsum, int *index, int jnum, int dim)
{        
    #pragma omp parallel for
    for (int ii=0; ii<jnum; ii++) {  // for each atom j in the simulation box     
        int j = jlist[ii];       // atom j        
        int start = bnumsum[ii];   
        int m = bnumsum[ii+1]-start;        // number of neighbors around i             
        for (int l=0; l<m ; l++) {   // loop over each atom around atom i
            int n = index[start + l];
            e[j] += 0.5*eij[n];                        
        }        
    }    
}
template void ompNeighborAtomPairDecomposition(double*, double*, int*, int*, int*, int, int);
template void ompNeighborAtomPairDecomposition(float*, float*, int*, int*, int*, int, int);

template <typename T> void ompTripletDecomposition(T *e, T *eijk, int *ai, int *aj, int *ak, int ijknum, int dim)
{    
    #pragma omp parallel for
    for (int ii=0; ii<ijknum; ii++) {  // for each atom pair ij in the simulation box     
        int i = ai[ii];       // atom i
        int j = aj[ii];       // atom j
        int k = ak[ii];       // atom k
        // use atomicAdd on gpu
        T ee = eijk[ii]/3.0;
        #pragma omp atomic
        e[i] += ee;
        #pragma omp atomic
        e[j] += ee;
        #pragma omp atomic
        e[k] += ee;
    }
}
template void ompTripletDecomposition(double*, double*, int*, int*, int*, int, int);
template void ompTripletDecomposition(float*, float*, int*, int*, int*, int, int);

template <typename T> void ompCenterAtomTripletDecomposition(T *e, T *eijk, int *ilist, int *anumsum, int inum, int dim)
{    
    #pragma omp parallel for
    for (int ii=0; ii<inum; ii++) {  // for each atom i in the simulation box     
        int i = ilist[ii];       // atom i        
        int start = anumsum[ii];   
        int m = anumsum[ii+1]-start;        // number of neighbors around i             
        for (int l=0; l<m ; l++)    // loop over each atom around atom i            
            e[i] += eijk[start + l]/3.0;                                
    }
}
template void ompCenterAtomTripletDecomposition(double*, double*, int*, int*, int, int);
template void ompCenterAtomTripletDecomposition(float*, float*,  int*, int*, int, int);

template <typename T> void ompNeighborAtomTripletDecomposition(T *e, T *eij, int *jlist, int *bnumsum, int *index, int jnum, int dim)
{        
    #pragma omp parallel for
    for (int ii=0; ii<jnum; ii++) {  // for each atom j in the simulation box     
        int j = jlist[ii];       // atom j        
        int start = bnumsum[ii];   
        int m = bnumsum[ii+1]-start;        // number of neighbors around i             
        for (int l=0; l<m ; l++) {   // loop over each atom around atom i
            int n = index[start + l];
            e[j] += eij[n]/3.0;                        
        }        
    }    
}
template void ompNeighborAtomTripletDecomposition(double*, double*, int*, int*, int*, int, int);
template void ompNeighborAtomTripletDecomposition(float*, float*, int*, int*, int*, int, int);

template <typename T> void ompQuadrupletDecomposition(T *e, T *eijkl, 
        int *ai, int *aj, int *ak, int *al, int ijklnum, int dim)
{    
    #pragma omp parallel for
    for (int ii=0; ii<ijklnum; ii++) {  // for each atom quadruplet ijkl in the simulation box     
        int i = ai[ii];       // atom i
        int j = aj[ii];       // atom j
        int k = ak[ii];       // atom k
        int l = al[ii];       // atom l
        // use atomicAdd on gpu
        T ee = eijkl[ii]/4.0;
        #pragma omp atomic
        e[i] += ee;
        #pragma omp atomic
        e[j] += ee;
        #pragma omp atomic
        e[k] += ee;
        #pragma omp atomic
        e[l] += ee;
    }
}
template void ompQuadrupletDecomposition(double*, double*, int*, int*, int*, int*, int, int);
template void ompQuadrupletDecomposition(float*, float*, int*, int*, int*, int*, int, int);

template <typename T> void ompCenterAtomQuadrupletDecomposition(T *e, T *eijkl, int *ilist, int *anumsum, int inum, int dim)
{    
    #pragma omp parallel for
    for (int ii=0; ii<inum; ii++) {  // for each atom i in the simulation box     
        int i = ilist[ii];       // atom i        
        int start = anumsum[ii];   
        int m = anumsum[ii+1]-start;        // number of neighbors around i             
        for (int l=0; l<m ; l++)    // loop over each atom around atom i            
            e[i] += eijkl[start + l]/4.0;                                
    }
}
template void ompCenterAtomQuadrupletDecomposition(double*, double*, int*, int*, int, int);
template void ompCenterAtomQuadrupletDecomposition(float*, float*, int*, int*, int, int);

template <typename T> void ompNeighborAtomQuadrupletDecomposition(T *e, T *eij, int *jlist, int *bnumsum, int *index, int jnum, int dim)
{        
    #pragma omp parallel for
    for (int ii=0; ii<jnum; ii++) {  // for each atom j in the simulation box     
        int j = jlist[ii];       // atom j        
        int start = bnumsum[ii];   
        int m = bnumsum[ii+1]-start;        // number of neighbors around i                     
        for (int l=0; l<m ; l++) {   // loop over each atom around atom i
            int n = index[start + l];
            e[j] += eij[n]/4.0;                        
        }        
    }    
}
template void ompNeighborAtomQuadrupletDecomposition(double*, double*, int*, int*, int*, int, int);
template void ompNeighborAtomQuadrupletDecomposition(float*, float*, int*, int*, int*, int, int);

template <typename T> void ompSingleDecomposition(T *e, T *f, T *ei, T *fi, int *ai, int inum, int dim)
{    
    #pragma omp parallel for
    for (int ii=0; ii<inum; ii++) {  // for each atom pair ij in the simulation box     
        int i = ai[ii];       // atom i        
        int l = dim*ii;                
        e[i] += ei[ii];
        i = dim*i;
        for (int d=0; d<dim; d++) 
            f[i+d] +=  fi[l+d];        
    }
}
template void ompSingleDecomposition(double*, double*, double*, double*, int*, int, int);
template void ompSingleDecomposition(float*, float*, float*, float*, int*, int, int);

template <typename T> void ompCenterAtomPairDecomposition(T *e, T *f, T *eij, T *fij, int *ilist, int *anumsum, int inum, int dim)
{    
    #pragma omp parallel for
    for (int ii=0; ii<inum; ii++) {  // for each atom i in the simulation box     
        int i = ilist[ii];       // atom i        
        int start = anumsum[ii];   
        int m = anumsum[ii+1]-start;        // number of neighbors around i             
        int idim = dim*i;        
        for (int l=0; l<m ; l++) {   // loop over each atom around atom i
            int n = start + l;
            int ndim = dim*n;
            e[i] += 0.5*eij[n];                        
            for (int d=0; d<dim; d++) 
                f[idim+d] +=  -fij[ndim+d];        
        }        
    }
}
template void ompCenterAtomPairDecomposition(double*, double*, double*, double*, int*, int*, int, int);
template void ompCenterAtomPairDecomposition(float*, float*, float*, float*, int*, int*, int, int);

template <typename T> void ompNeighborAtomPairDecomposition(T *e, T *f, T *eij, T *fij, int *jlist, int *bnumsum, int *index, int jnum, int dim)
{        
    #pragma omp parallel for
    for (int ii=0; ii<jnum; ii++) {  // for each atom j in the simulation box     
        int j = jlist[ii];       // atom j        
        int start = bnumsum[ii];   
        int m = bnumsum[ii+1]-start;        // number of neighbors around i             
        int jdim = dim*j;        
        for (int l=0; l<m ; l++) {   // loop over each atom around atom i
            int n = index[start + l];
            int ndim = dim*n;
            e[j] += 0.5*eij[n];                        
            for (int d=0; d<dim; d++) 
                f[jdim+d] +=  fij[ndim+d];        
        }        
    }    
}
template void ompNeighborAtomPairDecomposition(double*, double*, double*, double*, int*, int*, int*, int, int);
template void ompNeighborAtomPairDecomposition(float*, float*, float*, float*, int*, int*, int*, int, int);

template <typename T> void ompFullNeighPairDecomposition(T *e, T *f, T *eij, T *fij, int *ai, int ijnum, int dim)
{    
    #pragma omp parallel for
    for (int ii=0; ii<ijnum; ii++) {  // for each atom pair ij in the simulation box     
        int i = ai[ii];       // atom i
        int l = dim*ii;
        // use atomicAdd on gpu
        #pragma omp atomic
        e[i] += 0.5*eij[ii];
        i = dim*i;
        for (int d=0; d<dim; d++) 
            #pragma omp atomic
            f[i+d] +=  -fij[l+d];        
    }
}
template void ompFullNeighPairDecomposition(double*, double*, double*, double*, int*, int, int);
template void ompFullNeighPairDecomposition(float*, float*, float*, float*, int*, int, int);

template <typename T> void ompHalfNeighPairDecomposition(T *e, T *f, T *eij, T *fij, int *ai, int *aj, int ijnum, int dim)
{    
    #pragma omp parallel for
    for (int ii=0; ii<ijnum; ii++) {  // for each atom pair ij in the simulation box     
        int i = ai[ii];       // atom i
        int j = aj[ii];       // atom j
        int l = dim*ii;        
        // use atomicAdd on gpu
        #pragma omp atomic
        e[i] += 0.5*eij[ii];
        #pragma omp atomic
        e[j] += 0.5*eij[ii];
        i = dim*i;
        j = dim*j;
        for (int d=0; d<dim; d++) {
            #pragma omp atomic
            f[i+d] += -fij[l+d]; 
            #pragma omp atomic
            f[j+d] +=  fij[l+d]; 
        }        
    }
}
template void ompHalfNeighPairDecomposition(double*, double*, double*, double*, int*, int*, int, int);
template void ompHalfNeighPairDecomposition(float*, float*, float*, float*, int*, int*, int, int);

template <typename T> void ompTripletDecomposition(T *e, T *f, T *eijk, T *fij, T *fik, int *ai, int *aj, int *ak, int ijknum, int dim)
{    
    #pragma omp parallel for
    for (int ii=0; ii<ijknum; ii++) {  // for each atom pair ij in the simulation box     
        int i = ai[ii];       // atom i
        int j = aj[ii];       // atom j
        int k = ak[ii];       // atom k
        int n = dim*ii;        
        // use atomicAdd on gpu
        T ee = eijk[ii]/3.0;
        #pragma omp atomic
        e[i] += ee;
        #pragma omp atomic
        e[j] += ee;
        #pragma omp atomic
        e[k] += ee;
        i = dim*i;
        j = dim*j;
        k = dim*k;
        for (int d=0; d<dim; d++) {
            #pragma omp atomic
            f[i+d] += -(fij[n+d] + fik[n+d]); 
            #pragma omp atomic
            f[j+d] +=  fij[n+d]; 
            #pragma omp atomic
            f[k+d] +=  fik[n+d]; 
        }        
    }
}
template void ompTripletDecomposition(double*, double*, double*, double*, double*, int*, int*, int*, int, int);
template void ompTripletDecomposition(float*, float*, float*, float*, float*, int*, int*, int*, int, int);

template <typename T> void ompCenterAtomTripletDecomposition(T *e, T *f, T *eijk, T *fij, T *fik, int *ilist, int *anumsum, int inum, int dim)
{    
    #pragma omp parallel for
    for (int ii=0; ii<inum; ii++) {  // for each atom i in the simulation box     
        int i = ilist[ii];       // atom i        
        int start = anumsum[ii];   
        int m = anumsum[ii+1]-start;        // number of neighbors around i             
        int idim = dim*i;        
        for (int l=0; l<m ; l++) {   // loop over each atom around atom i
            int n = start + l;
            int ndim = dim*n;
            e[i] += eijk[n]/3.0;                        
            for (int d=0; d<dim; d++) 
                f[idim+d] +=  -(fij[ndim+d]+fik[ndim+d]);        
        }        
    }
}
template void ompCenterAtomTripletDecomposition(double*, double*, double*, double*, double*, int*, int*, int, int);
template void ompCenterAtomTripletDecomposition(float*, float*, float*, float*, float*, int*, int*, int, int);

template <typename T> void ompNeighborAtomTripletDecomposition(T *e, T *f, T *eij, T *fij, int *jlist, int *bnumsum, int *index, int jnum, int dim)
{        
    #pragma omp parallel for
    for (int ii=0; ii<jnum; ii++) {  // for each atom j in the simulation box     
        int j = jlist[ii];       // atom j        
        int start = bnumsum[ii];   
        int m = bnumsum[ii+1]-start;        // number of neighbors around i             
        int jdim = dim*j;        
        for (int l=0; l<m ; l++) {   // loop over each atom around atom i
            int n = index[start + l];
            int ndim = dim*n;
            e[j] += eij[n]/3.0;                        
            for (int d=0; d<dim; d++) 
                f[jdim+d] +=  fij[ndim+d];        
        }        
    }    
}
template void ompNeighborAtomTripletDecomposition(double*, double*, double*, double*, int*, int*, int*, int, int);
template void ompNeighborAtomTripletDecomposition(float*, float*, float*, float*, int*, int*, int*, int, int);

template <typename T> void ompQuadrupletDecomposition(T *e, T *f, T *eijkl, T *fij, T *fik, T *fil, 
        int *ai, int *aj, int *ak, int *al, int ijklnum, int dim)
{    
    #pragma omp parallel for
    for (int ii=0; ii<ijklnum; ii++) {  // for each atom quadruplet ijkl in the simulation box     
        int i = ai[ii];       // atom i
        int j = aj[ii];       // atom j
        int k = ak[ii];       // atom k
        int l = al[ii];       // atom l
        int n = dim*ii;        
        // use atomicAdd on gpu
        T ee = eijkl[ii]/4.0;
        #pragma omp atomic
        e[i] += ee;
        #pragma omp atomic
        e[j] += ee;
        #pragma omp atomic
        e[k] += ee;
        #pragma omp atomic
        e[l] += ee;
        i = dim*i;
        j = dim*j;
        k = dim*k;
        l = dim*l;
        for (int d=0; d<dim; d++) {
            #pragma omp atomic
            f[i+d] += -(fij[n+d] + fik[n+d] + fil[n+d]); 
            #pragma omp atomic
            f[j+d] +=  fij[n+d]; 
            #pragma omp atomic
            f[k+d] +=  fik[n+d]; 
            #pragma omp atomic
            f[l+d] +=  fil[n+d]; 
        }        
    }
}
template void ompQuadrupletDecomposition(double*, double*, double*, double*, double*, double*, int*, int*, int*, int*, int, int);
template void ompQuadrupletDecomposition(float*, float*, float*, float*, float*, float*, int*, int*, int*, int*, int, int);

template <typename T> void ompCenterAtomQuadrupletDecomposition(T *e, T *f, T *eijkl, T *fij, T *fik, T *fil, 
        int *ilist, int *anumsum, int inum, int dim)
{    
    #pragma omp parallel for
    for (int ii=0; ii<inum; ii++) {  // for each atom i in the simulation box     
        int i = ilist[ii];       // atom i        
        int start = anumsum[ii];   
        int m = anumsum[ii+1]-start;        // number of neighbors around i             
        int idim = dim*i;        
        for (int l=0; l<m ; l++) {   // loop over each atom around atom i
            int n = start + l;
            int ndim = dim*n;
            e[i] += eijkl[n]/4.0;                        
            for (int d=0; d<dim; d++) 
                f[idim+d] +=  -(fij[ndim+d]+fik[ndim+d]+fil[ndim+d]);        
        }        
    }
}
template void ompCenterAtomQuadrupletDecomposition(double*, double*, double*, double*, double*, double*, int*, int*, int, int);
template void ompCenterAtomQuadrupletDecomposition(float*, float*, float*, float*, float*, float*, int*, int*, int, int);

template <typename T> void ompNeighborAtomQuadrupletDecomposition(T *e, T *f, T *eij, T *fij, int *jlist, int *bnumsum, int *index, int jnum, int dim)
{        
    #pragma omp parallel for
    for (int ii=0; ii<jnum; ii++) {  // for each atom j in the simulation box     
        int j = jlist[ii];       // atom j        
        int start = bnumsum[ii];   
        int m = bnumsum[ii+1]-start;        // number of neighbors around i             
        int jdim = dim*j;        
        for (int l=0; l<m ; l++) {   // loop over each atom around atom i
            int n = index[start + l];
            int ndim = dim*n;
            e[j] += eij[n]/4.0;                        
            for (int d=0; d<dim; d++) 
                f[jdim+d] +=  fij[ndim+d];        
        }        
    }    
}
template void ompNeighborAtomQuadrupletDecomposition(double*, double*, double*, double*, int*, int*, int*, int, int);
template void ompNeighborAtomQuadrupletDecomposition(float*, float*, float*, float*, int*, int*, int*, int, int);


template <typename T> void ompSingleDecomposition2D(T *f, T *fi, int *ai, int inum)
{    
    #pragma omp parallel for
    for (int ii=0; ii<inum; ii++) {  // for each atom pair ij in the simulation box     
        int i = 2*ai[ii];       // atom i
        int l = 2*ii;
        // use atomicAdd on gpu
        f[i+0] +=  fi[l+0];
        f[i+1] +=  fi[l+1];
    }
}
template void ompSingleDecomposition2D(double*, double*, int*, int);
template void ompSingleDecomposition2D(float*, float*, int*, int);

template <typename T> void ompSingleDecomposition3D(T *f, T *fi, int *ai, int inum)
{    
    #pragma omp parallel for
    for (int ii=0; ii<inum; ii++) {  // for each atom pair ij in the simulation box     
        int i = 3*ai[ii];       // atom i
        int l = 3*ii;
        // use atomicAdd on gpu
        f[i+0] +=  fi[l+0];
        f[i+1] +=  fi[l+1];
        f[i+2] +=  fi[l+2];
    }
}
template void ompSingleDecomposition3D(double*, double*, int*, int);
template void ompSingleDecomposition3D(float*, float*, int*, int);


template <typename T> void ompHalfForceDecomposition2D(T *f, T *fij, int *ai, int *aj, int ijnum)
{    
    #pragma omp parallel for
    for (int ii=0; ii<ijnum; ii++) {  // for each atom pair ij in the simulation box     
        int i = 2*ai[ii];       // atom i
        int j = 2*aj[ii];       // atom j
        int l = 2*ii;
        // use atomicAdd on gpu
        #pragma omp atomic
        f[i+0] +=  -fij[l+0];
        #pragma omp atomic
        f[i+1] +=  -fij[l+1];
        #pragma omp atomic
        f[j+0] +=  fij[l+0];
        #pragma omp atomic
        f[j+1] +=  fij[l+1];
    }
}
template void ompHalfForceDecomposition2D(double*, double*, int*, int*, int);
template void ompHalfForceDecomposition2D(float*, float*, int*, int*, int);

template <typename T> void ompFullForceDecomposition2D(T *f, T *fij, int *ai, int *aj, int ijnum)
{    
    #pragma omp parallel for
    for (int ii=0; ii<ijnum; ii++) {  // for each atom pair ij in the simulation box     
        int i = 2*ai[ii];       // atom i
        int l = 2*ii;
        // use atomicAdd on gpu
        #pragma omp atomic
        f[i+0] +=  -fij[l+0];
        #pragma omp atomic
        f[i+1] +=  -fij[l+1];
    }
}
template void ompFullForceDecomposition2D(double*, double*, int*, int*, int);
template void ompFullForceDecomposition2D(float*, float*, int*, int*, int);

template <typename T> void ompIAtomDecomposition2D(T *f, T *fij, int *ilist, int *anumsum, int inum)
{    
    #pragma omp parallel for
    for (int ii=0; ii<inum; ii++) {  // for each atom i in the simulation box     
        int i = 2*ilist[ii];       // atom i        
        int start = anumsum[ii];   
        int m = anumsum[ii+1]-start;        // number of neighbors around i             
        T f1 = 0.0;
        T f2 = 0.0;
        for (int l=0; l<m ; l++) {   // loop over each atom around atom i
            int k = 2*(start + l);                     
            f1 +=  -fij[k+0];
            f2 +=  -fij[k+1];
        }
        f[i+0] = f1;
        f[i+1] = f2;
    }
}
template void ompIAtomDecomposition2D(double*, double*, int*, int*, int);
template void ompIAtomDecomposition2D(float*, float*, int*, int*, int);

template <typename T> void ompJAtomDecomposition2D(T *f, T *fij, int *jlist, int *bnumsum, int *index, int jnum)
{        
    #pragma omp parallel for
    for (int ii=0; ii<jnum; ii++) {  // for each atom j in the simulation box     
        int j = 2*jlist[ii];       // atom j
        T f1 = 0.0;
        T f2 = 0.0;                
        // need to determine bnumsum and index 
        int start = bnumsum[ii];   
        int m = bnumsum[ii+1]-start;   // number of neighbors i around j  (j < i)           
        for (int l=0; l<m ; l++) { // loop over each atom around i atom j -> pair ij 
            int k = 2*index[start + l];                     
            f1 +=  fij[k+0];
            f2 +=  fij[k+1];
        }
        f[j+0] = f1;
        f[j+1] = f2;
    }
}
template void ompJAtomDecomposition2D(double*, double*, int*, int*, int*, int);
template void ompJAtomDecomposition2D(float*, float*, int*, int*, int*, int);

template <typename T> void ompForceDecompositionTriplet2D(T *f, T *fij, T *fik, int *ai, int *aj, int *ak, int ijknum)
{    
    #pragma omp parallel for
    for (int ii=0; ii<ijknum; ii++) {  // for each atom pair ij in the simulation box     
        int i = 2*ai[ii];       // atom i
        int j = 2*aj[ii];       // atom j
        int k = 2*ak[ii];       // atom k
        int l = 2*ii;
        // use atomicAdd on gpu
        #pragma omp atomic
        f[i+0] +=  -(fij[l+0]+fik[l+0]);
        #pragma omp atomic
        f[i+1] +=  -(fij[l+1]+fik[l+1]);
        #pragma omp atomic
        f[j+0] +=  fij[l+0];
        #pragma omp atomic
        f[j+1] +=  fij[l+1];
        #pragma omp atomic
        f[k+0] +=  fik[l+0];
        #pragma omp atomic
        f[k+1] +=  fik[l+1];
    }
}
template void ompForceDecompositionTriplet2D(double*, double*, double*, int*, int*, int*, int);
template void ompForceDecompositionTriplet2D(float*, float*, float*, int*, int*, int*, int);

template <typename T> void ompAtomDecompositionTriplet2D(T *f, T *fij, T *fik, int *ilist, int *anumsum, int inum)
{    
    #pragma omp parallel for
    for (int ii=0; ii<inum; ii++) {  // for each atom i in the simulation box     
        int i = 2*ilist[ii];       // atom i        
        int start = anumsum[ii];   
        int m = anumsum[ii+1]-start;        // number of neighbors around i             
        T f1 = 0.0;
        T f2 = 0.0;
        for (int l=0; l<m ; l++) {   // loop over each atom pair jk around atom i
            int s = start + l;
            int n = 2*s;                     
            f1 +=  -(fij[n+0] + fik[n+0]);
            f2 +=  -(fij[n+1] + fik[n+1]);
//             // use atomicAdd on gpu
//             int j = 2*aj[s];
//             int k = 2*ak[s];            
//             f[j+0] +=  fij[n+0];
//             f[j+1] +=  fij[n+1];
//             f[k+0] +=  fik[n+0];
//             f[k+1] +=  fik[n+1];            
        }
        f[i+0] = f1;
        f[i+1] = f2;                        
    }
}
template void ompAtomDecompositionTriplet2D(double*, double*, double*, int*, int*, int);
template void ompAtomDecompositionTriplet2D(float*, float*, float*, int*, int*, int);

template <typename T> void ompForceDecompositionQuadruplet2D(T *f, T *fij, T *fik,  T *fil, int *ai, int *aj, int *ak, int *al, int ijklnum)
{    
    #pragma omp parallel for
    for (int ii=0; ii<ijklnum; ii++) {  // for each atom pair ij in the simulation box     
        int i = 2*ai[ii];       // atom i
        int j = 2*aj[ii];       // atom j
        int k = 2*ak[ii];       // atom k
        int l = 2*al[ii];       // atom l
        int n = 2*ii;
        // use atomicAdd on gpu
        #pragma omp atomic
        f[i+0] +=  -(fij[n+0]+fik[n+0]+fil[n+0]);
        #pragma omp atomic
        f[i+1] +=  -(fij[n+1]+fik[n+1]+fil[n+1]);
        #pragma omp atomic
        f[j+0] +=  fij[n+0];
        #pragma omp atomic
        f[j+1] +=  fij[n+1];
        #pragma omp atomic
        f[k+0] +=  fik[n+0];
        #pragma omp atomic
        f[k+1] +=  fik[n+1];
        #pragma omp atomic
        f[l+0] +=  fil[n+0];
        #pragma omp atomic
        f[l+1] +=  fil[n+1];
    }
}
template void ompForceDecompositionQuadruplet2D(double*, double*, double*, double*, int*, int*,  int*, int*, int);
template void ompForceDecompositionQuadruplet2D(float*, float*, float*, float*, int*, int*, int*, int*, int);

template <typename T> void ompAtomDecompositionQuadruplet2D(T *f, T *fij, T *fik, T *fil, int *ilist, int *anumsum, int inum)
{    
    #pragma omp parallel for
    for (int ii=0; ii<inum; ii++) {  // for each atom i in the simulation box     
        int i = 2*ilist[ii];       // atom i        
        int start = anumsum[ii];   
        int m = anumsum[ii+1]-start;        // number of neighbors around i             
        T f1 = 0.0;
        T f2 = 0.0;
        for (int l=0; l<m ; l++) {   // loop over each atom pair jk around atom i
            int s = start + l;
            int n = 2*s;                     
            f1 +=  -(fij[n+0] + fik[n+0] + fil[n+0]);
            f2 +=  -(fij[n+1] + fik[n+1] + fil[n+1]);
        }
        f[i+0] = f1;
        f[i+1] = f2;                        
    }
}
template void ompAtomDecompositionQuadruplet2D(double*, double*, double*, double*, int*, int*, int);
template void ompAtomDecompositionQuadruplet2D(float*, float*, float*, float*, int*, int*, int);

template <typename T> void ompHalfForceDecomposition3D(T *f, T *fij, int *ai, int *aj, int ijnum)
{    
    #pragma omp parallel for
    for (int ii=0; ii<ijnum; ii++) {  // for each atom pair ij in the simulation box     
        int i = 3*ai[ii];       // atom i
        int j = 3*aj[ii];       // atom j
        int l = 3*ii;
        // use atomicAdd on gpu
        #pragma omp atomic
        f[i+0] +=  -fij[l+0];
        #pragma omp atomic
        f[i+1] +=  -fij[l+1];
        #pragma omp atomic
        f[i+2] +=  -fij[l+2];                 
        #pragma omp atomic
        f[j+0] +=  fij[l+0];
        #pragma omp atomic
        f[j+1] +=  fij[l+1];
        #pragma omp atomic
        f[j+2] +=  fij[l+2];                 
    }
}
template void ompHalfForceDecomposition3D(double*, double*, int*, int*, int);
template void ompHalfForceDecomposition3D(float*, float*, int*, int*, int);

template <typename T> void ompFullForceDecomposition3D(T *f, T *fij, int *ai, int *aj, int ijnum)
{    
    #pragma omp parallel for
    for (int ii=0; ii<ijnum; ii++) {  // for each atom pair ij in the simulation box     
        int i = 3*ai[ii];             // atom i
        int l = 3*ii;
        // use atomicAdd on gpu
        #pragma omp atomic
        f[i+0] +=  -fij[l+0];
        #pragma omp atomic
        f[i+1] +=  -fij[l+1];
        #pragma omp atomic
        f[i+2] +=  -fij[l+2];                         
    }
}
template void ompFullForceDecomposition3D(double*, double*, int*, int*, int);
template void ompFullForceDecomposition3D(float*, float*, int*, int*, int);

template <typename T> void ompIAtomDecomposition3D(T *f, T *fij, int *ilist, int *anumsum, int inum)
{    
    #pragma omp parallel for
    for (int ii=0; ii<inum; ii++) {  // for each atom i in the simulation box     
        int i = 3*ilist[ii];       // atom i        
        T f1 = 0.0;
        T f2 = 0.0;        
        T f3 = 0.0;        
        int start = anumsum[ii];   
        int m = anumsum[ii+1]-start;        // number of neighbors around i             
        for (int l=0; l<m ; l++) {   // loop over each atom around atom i
            int k = 3*(start + l);                     
            f1 +=  -fij[k+0];
            f2 +=  -fij[k+1];
            f3 +=  -fij[k+2];             
        }
        f[i+0] = f1;
        f[i+1] = f2;
        f[i+2] = f3;
    }
}
template void ompIAtomDecomposition3D(double*, double*, int*, int*, int);
template void ompIAtomDecomposition3D(float*, float*, int*, int*, int);

template <typename T> void ompJAtomDecomposition3D(T *f, T *fij, int *jlist, int *bnumsum, int *index, int jnum)
{        
    #pragma omp parallel for
    for (int ii=0; ii<jnum; ii++) {  // for each atom j in the simulation box     
        int j = 3*jlist[ii];       // atom j
        T f1 = 0.0;
        T f2 = 0.0;                
        T f3 = 0.0;                
        // need to determine bnumsum and index 
        int start = bnumsum[ii];   
        int m = bnumsum[ii+1]-start;   // number of neighbors i around j  (j < i)           
        for (int l=0; l<m ; l++) { // loop over each atom around i atom j -> pair ij 
            int k = 3*index[start + l];                     
            f1 +=  fij[k+0];
            f2 +=  fij[k+1];
            f3 +=  fij[k+2];
        }
        f[j+0] = f1;
        f[j+1] = f2;
        f[j+2] = f3;
    }
}
template void ompJAtomDecomposition3D(double*, double*, int*, int*, int*, int);
template void ompJAtomDecomposition3D(float*, float*, int*, int*, int*, int);

template <typename T> void ompForceDecompositionTriplet3D(T *f, T *fij, T *fik, int *ai, int *aj, int *ak, int ijknum)
{    
    #pragma omp parallel for
    for (int ii=0; ii<ijknum; ii++) {  // for each atom pair ij in the simulation box     
        int i = 3*ai[ii];       // atom i
        int j = 3*aj[ii];       // atom j
        int k = 3*ak[ii];       // atom k
        int l = 3*ii;
        // use atomicAdd on gpu
        #pragma omp atomic
        f[i+0] +=  -(fij[l+0]+fik[l+0]);
        #pragma omp atomic
        f[i+1] +=  -(fij[l+1]+fik[l+1]);
        #pragma omp atomic
        f[i+2] +=  -(fij[l+2]+fik[l+2]);
        #pragma omp atomic
        f[j+0] +=  fij[l+0];
        #pragma omp atomic
        f[j+1] +=  fij[l+1];
        #pragma omp atomic
        f[j+2] +=  fij[l+2];
        #pragma omp atomic
        f[k+0] +=  fik[l+0];
        #pragma omp atomic
        f[k+1] +=  fik[l+1];
        #pragma omp atomic
        f[k+2] +=  fik[l+2];
    }
}
template void ompForceDecompositionTriplet3D(double*, double*, double*, int*, int*, int*, int);
template void ompForceDecompositionTriplet3D(float*, float*, float*, int*, int*, int*, int);

template <typename T> void ompAtomDecompositionTriplet3D(T *f, T *fij, T *fik, int *ilist, int *anumsum, int inum)
{    
    #pragma omp parallel for
    for (int ii=0; ii<inum; ii++) {  // for each atom i in the simulation box     
        int i = 3*ilist[ii];       // atom i        
        int start = anumsum[ii];   
        int m = anumsum[ii+1]-start;        // number of neighbors around i             
        T f1 = 0.0;
        T f2 = 0.0;
        T f3 = 0.0;
        for (int l=0; l<m ; l++) {   // loop over each atom pair jk around atom i
            int s = start + l;
            int n = 3*s;                     
            f1 +=  -(fij[n+0] + fik[n+0]);
            f2 +=  -(fij[n+1] + fik[n+1]);
            f3 +=  -(fij[n+2] + fik[n+2]);
        }
        f[i+0] = f1;
        f[i+1] = f2;                        
        f[i+2] = f3;                        
    }
}
template void ompAtomDecompositionTriplet3D(double*, double*, double*, int*, int*, int);
template void ompAtomDecompositionTriplet3D(float*, float*, float*, int*, int*, int);

template <typename T> void ompForceDecompositionQuadruplet3D(T *f, T *fij, T *fik,  T *fil, int *ai, int *aj, int *ak, int *al, int ijklnum)
{    
    #pragma omp parallel for
    for (int ii=0; ii<ijklnum; ii++) {  // for each atom pair ij in the simulation box     
        int i = 3*ai[ii];       // atom i
        int j = 3*aj[ii];       // atom j
        int k = 3*ak[ii];       // atom k
        int l = 3*al[ii];       // atom l
        int n = 3*ii;
        // use atomicAdd on gpu
        #pragma omp atomic
        f[i+0] +=  -(fij[n+0]+fik[n+0]+fil[n+0]);
        #pragma omp atomic
        f[i+1] +=  -(fij[n+1]+fik[n+1]+fil[n+1]);
        #pragma omp atomic
        f[i+2] +=  -(fij[n+2]+fik[n+2]+fil[n+2]);
        #pragma omp atomic
        f[j+0] +=  fij[n+0];
        #pragma omp atomic
        f[j+1] +=  fij[n+1];
        #pragma omp atomic
        f[j+2] +=  fij[n+2];
        #pragma omp atomic
        f[k+0] +=  fik[n+0];
        #pragma omp atomic
        f[k+1] +=  fik[n+1];
        #pragma omp atomic
        f[k+2] +=  fik[n+2];
        #pragma omp atomic
        f[l+0] +=  fil[n+0];
        #pragma omp atomic
        f[l+1] +=  fil[n+1];
        #pragma omp atomic
        f[l+2] +=  fil[n+2];
    }
}
template void ompForceDecompositionQuadruplet3D(double*, double*, double*, double*, int*, int*,  int*, int*, int);
template void ompForceDecompositionQuadruplet3D(float*, float*, float*, float*, int*, int*, int*, int*, int);

template <typename T> void ompAtomDecompositionQuadruplet3D(T *f, T *fij, T *fik, T *fil, int *ilist, int *anumsum, int inum)
{    
    #pragma omp parallel for
    for (int ii=0; ii<inum; ii++) {  // for each atom i in the simulation box     
        int i = 3*ilist[ii];       // atom i        
        int start = anumsum[ii];   
        int m = anumsum[ii+1]-start;        // number of neighbors around i             
        T f1 = 0.0;
        T f2 = 0.0;
        T f3 = 0.0;
        for (int l=0; l<m ; l++) {   // loop over each atom pair jk around atom i
            int s = start + l;
            int n = 3*s;                     
            f1 +=  -(fij[n+0] + fik[n+0] + fil[n+0]);
            f2 +=  -(fij[n+1] + fik[n+1] + fil[n+1]);
            f3 +=  -(fij[n+2] + fik[n+2] + fil[n+2]);
        }
        f[i+0] = f1;
        f[i+1] = f2;                        
        f[i+2] = f3;                        
    }
}
template void ompAtomDecompositionQuadruplet3D(double*, double*, double*, double*, int*, int*, int);
template void ompAtomDecompositionQuadruplet3D(float*, float*, float*, float*, int*, int*, int);

//********************************************************************************************//

// template <typename T> void ompGetNeighPairs(T *xij, T *qi, T *qj, T *x, T *q, int *ai, int *aj,  
//       int *ti, int *tj, int *anum, int *anumsum, int *ilist, int *alist, int *neighlist, int *neighnum, 
//       int *atomtype, int inum, int jnum, int ncq, int dim)
// {        
//     for (int ii=0; ii<inum; ii++) {  // for each atom i in the simulation box     
//         int i = ilist[ii];       // atom i
//         int itype = atomtype[i];
//         int m = anum[ii];        // number of neighbors around i             
//         int start = anumsum[ii];   
//         for (int l=0; l<m ; l++) {   // loop over each atom around atom i
//             int g = neighlist[l +jnum*i];  // ghost index of atom j  
//             int j = alist[g];  // atom j
//             int k = start + l;                                     
//             ai[k]        = i;
//             aj[k]        = j; // should be alist[j];         
//             ti[k]        = itype;       
//             tj[k]        = atomtype[j]; // should be atomtype[alist[j]];         
//             for (int p=0; p<dim; p++) 
//                 xij[k*dim+p]   = x[g*dim+p] - x[i*dim+p];  // xj - xi            
//             for (int p=0; p<ncq; p++) {                
//                 qi[k*ncq+p] = q[i*ncq+p];
//                 qj[k*ncq+p] = q[j*ncq+p];
//             }                
//         }
//     }    
// }
// template void ompGetNeighPairs(double*, double*, double*, double*, double*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int, int, int, int);
// template void ompGetNeighPairs(float*, float*, float*, float*, float*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int, int, int, int);
// 
// template <typename T> void ompGetHalfNeighPairs(T *xij, T *qi, T *qj, T *x, T *q, int *ai, int *aj,  
//       int *ti, int *tj, int *anum, int *anumsum, int *ilist, int *alist, int *neighlist, int *neighnum, 
//       int *atomtype, int inum, int jnum, int ncq, int dim)
// {        
//     for (int ii=0; ii<inum; ii++) {  // for each atom i in the simulation box     
//         int i = ilist[ii];       // atom i
//         int itype = atomtype[i];
//         int m = anum[ii];        // number of neighbors around i             
//         int start = anumsum[ii];   
//         int inc=0;
//         for (int l=0; l<m ; l++) {   // loop over each atom around atom i
//             int g = neighlist[l +jnum*i];  // ghost index of atom j  
//             int j = alist[g];  // atom j
//             if (i < j) {
//                 int k = start + inc;                         
//                 ai[k]        = i;
//                 aj[k]        = j; // should be alist[j];         
//                 ti[k]        = itype;       
//                 tj[k]        = atomtype[j]; // should be atomtype[alist[j]];                                        
//                 for (int p=0; p<dim; p++) 
//                     xij[k*dim+p]   = x[g*dim+p] - x[i*dim+p];  // xj - xi                            
//                 for (int p=0; p<ncq; p++) {                
//                     qi[k*ncq+p] = q[i*ncq+p];
//                     qj[k*ncq+p] = q[j*ncq+p];
//                 }                
//                 inc += 1;
//             }
//         }
//     }    
// }
// template void ompGetHalfNeighPairs(double*, double*, double*, double*, double*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int, int, int, int);
// template void ompGetHalfNeighPairs(float*, float*, float*, float*, float*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int, int, int, int);
// 
// template <typename T> void ompGetNeighPairs(T *xij, T *qi, T *qj, T *x, T *q, int *ai, int *aj,  
//       int *ti, int *tj, int *anum, int *anumsum, int *ilist, int *alist, int *neighlist, int *neighnum, 
//       int *atomtype, int typej, int inum, int jnum, int ncq, int dim)
// {        
// //     // form anum
// //     for (int ii=0; ii<inum; ii++) {
// //         int i = ilist[ii];       // atom i
// //         int m = neighnum[i];     // number of neighbors around i             
// //         anum[ii] = 0;              
// //         for (int l=0; l<m ; l++) {   // loop over each atom around atom i
// //             int g = neighlist[l + jnum*i];  // atom j           
// //             int j = alist[g];  // atom j
// //             if (atomtype[j] == typej) 
// //                 anum[ii] += 1;
// //         }                
// //     }
// //     
// //     // a list contains the starting positions of the first neighbor 
// //     ompCumsum(anumsum, anum, inum+1);                 
//     for (int ii=0; ii<inum; ii++) {  // for each atom i in the simulation box     
//         int i = ilist[ii];       // atom i
//         int itype = atomtype[i];
//         int m = anum[ii];        // number of neighbors around i             
//         int start = anumsum[ii];   
//         int inc = 0;
//         for (int l=0; l<m ; l++) {   // loop over each atom j around atom i
//             int g = neighlist[l +jnum*i];  // ghost index of atom j  
//             int j = alist[g];  // atom j
//             if (atomtype[j] == typej)  {
//                 int k = start + inc;                         
//                 ai[k]        = i;
//                 aj[k]        = j; // should be alist[j];         
//                 ti[k]        = itype;       
//                 tj[k]        = atomtype[j]; // should be atomtype[alist[j]];                                        
//                 for (int p=0; p<dim; p++) 
//                     xij[k*dim+p]   = x[g*dim+p] - x[i*dim+p];  // xj - xi                            
//                 for (int p=0; p<ncq; p++) {                
//                     qi[k*ncq+p] = q[i*ncq+p];
//                     qj[k*ncq+p] = q[j*ncq+p];
//                 }                
//                 inc += 1;
//             }
//         }
//     }    
// }
// template void ompGetNeighPairs(double*, double*, double*, double*, double*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int, int, int, int, int);
// template void ompGetNeighPairs(float*, float*, float*, float*, float*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int, int, int, int, int);
// 
// void ompFindAtom(int *tlist, int* ilist, int *atomtype, int typei, int inum)
// {
//     for (int ii=0; ii<inum; ii++) {  // for each atom i in the simulation box     
//         int i = ilist[ii];           // atom i
//         tlist[ii] = (atomtype[i] == typei) ? 1 : 0;        
//     }    
// }
// 
// void ompFindStart(int *start, int* slist, int inum)
// {
//     for (int ii=0; ii<inum; ii++) {        
//         if ((ii==0) && (slist[ii]==1)) {
//             start[0] = 0;
//             start[1] = inum;
//         }
//         if ((ii==inum-1) && (slist[ii]==0)) {
//             start[0] = inum;
//             start[1] = 0;
//         }    
//         if (slist[ii]-slist[ii-1]==1) {
//             start[0] = ii;
//             start[1] = inum - ii;
//         }    
//     }
// }
//         
// void ompCopyAtIndex(int* tlist, int* ilist, int* index, int *start)
// {
//     for (int ii=0; ii<start[1]; ii++) {                  
//         tlist[ii] = ilist[index[start[0]+ii]];
//     }    
// }
// 
// void ompGetIlist(int *tlist, int *start, int *ilist, int *slist, int *index, int *atomtype, int typei, int inum)
// {
//     ompFindAtom(tlist, ilist, atomtype, typei, inum);
//  
//     // sort list    
//     ompMergeSort(slist, index, tlist, inum);
//         
//     ompFindStart(start, slist, inum);        
//         
//     ompCopyAtIndex(tlist, ilist, index, start);
// }
// 
// template <typename T> void ompGetNeighPairs(T *xij, T *qi, T *qj, T *x, T *q, int *ai, int *aj, 
//       int *ti, int *tj, int *anum, int *anumsum, int *ilist, int *tlist, int *alist, int *neighlist, int *neighnum, 
//       int *atomtype, int *slist, int *start, int *index, int typei, int inum, int jnum, int ncq, int dim)
// {    
//     // form tlist from ilist
//     ompGetIlist(tlist, start, ilist, slist, index, atomtype, typei, inum);
//             
//     ompGetNeighPairs(xij, qi, qj, x, q, ai, aj, ti, tj, anum, anumsum, tlist, alist,
//             neighlist, neighnum, atomtype, start[1], jnum, ncq, dim);             
// }
// template void ompGetNeighPairs(double*, double*, double*, double*, double*, int*, int*, int*, int*, int*, 
//         int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int, int, int, int, int);
// template void ompGetNeighPairs(float*, float*, float*, float*, float*, int*, int*, int*, int*, int*, 
//         int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int, int, int, int, int);
// 
// template <typename T> void ompGetNeighPairs(T *xij, T *qi, T *qj, T *x, T *q, int *ai, int *aj, 
//       int *ti, int *tj, int *anum, int *anumsum, int *ilist, int *tlist, int *alist, int *neighlist, int *neighnum, 
//       int *atomtype, int *slist, int *start, int *index, int typei, int typej, int inum, int jnum, int ncq, int dim)
// {    
//     // form tlist from ilist
//     ompGetIlist(tlist, start, ilist, slist, index, atomtype, typei, inum);
//             
//     ompGetNeighPairs(xij, qi, qj, x, q, ai, aj, ti, tj, tlist, alist, anum, anumsum,
//             neighlist, neighnum, atomtype, typej, start[1], jnum, ncq, dim);             
// }
// template void ompGetNeighPairs(double*, double*, double*, double*, double*, int*, int*, int*, int*, int*, 
//         int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int, int, int, int, int, int);
// template void ompGetNeighPairs(float*, float*, float*, float*, float*, int*, int*, int*, int*, int*, 
//         int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int, int, int, int, int, int);
// 
// template <typename T> void ompGetNeighTriplets(T *xij, T *xik, T *qi, T *qj, T *qk, T *x, T *q, int *ai, int *aj,  
//       int *ak, int *ti, int *tj, int *tk, int *anum, int *anumsum, int *ilist, int *alist, int *neighlist, int *neighnum, 
//       int *atomtype, int inum, int jnum, int ncq, int dim)
// {        
//     // form anum
//     for (int ii=0; ii<inum; ii++) {
//         int i = ilist[ii];       // atom i
//         int m = neighnum[i];     // number of neighbors around i             
//         anum[ii] = 0;              
//         for (int lj=0; lj<m ; lj++) {   // loop over each atom around atom i
//             int gj = neighlist[lj + jnum*i];  // atom j           
//             int j = alist[gj];  // atom j
//             for (int lk=0; lk<m; lk++) {
//                 int gk = neighlist[lk + jnum*i];  // atom k
//                 int k = alist[gk];  // atom k
//                 if (j < k)
//                     anum[ii] += 1;
//             }            
//         }                
//     }
//     
//     // a list contains the starting positions of the first neighbor 
//     ompCumsum(anumsum, anum, inum+1);             
//     
//     for (int ii=0; ii<inum; ii++) {  // for each atom i in the simulation box     
//         int i = ilist[ii];       // atom i
//         int itype = atomtype[i];
//         int m = neighnum[i];     // number of neighbors around i         
//         int start = anumsum[ii];   
//         int inc = 0;
//         for (int lj=0; lj<m ; lj++) {   // loop over each atom around atom i
//             int gj = neighlist[lj +jnum*i];  // ghost index of atom j  
//             int j = alist[gj];  // atom j
//             for (int lk=0; lk<m; lk++) {
//                 int gk = neighlist[lk + jnum*i];  // atom k
//                 int k = alist[gk];  // atom k
//                 if (j < k) {
//                     int n = start + inc;                         
//                     for (int p=0; p<dim; p++) {
//                         xij[n*dim+p]   = x[gj*dim+p] - x[i*dim+p];  // xj - xi            
//                         xik[n*dim+p]   = x[gk*dim+p] - x[i*dim+p];  // xj - xi            
//                     }
//                     ai[n]        = i;
//                     aj[n]        = j; // should be alist[j];         
//                     ak[n]        = k; // should be alist[j];         
//                     ti[n]        = itype;       
//                     tj[n]        = atomtype[j]; // should be atomtype[alist[j]];                                        
//                     tk[n]        = atomtype[k]; // should be atomtype[alist[j]];                                        
//                     for (int p=0; p<ncq; p++) {                
//                         qi[n*ncq+p] = q[i*ncq+p];
//                         qj[n*ncq+p] = q[j*ncq+p];
//                         qk[n*ncq+p] = q[k*ncq+p];
//                     }                
//                     inc += 1;
//                 }
//             }                                    
//         }
//     }    
// }
// template void ompGetNeighTriplets(double*, double*, double*, double*, double*, double*, double*, 
//         int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int, int, int, int);
// template void ompGetNeighTriplets(float*, float*, float*, float*, float*, float*, float*, 
//         int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int, int, int, int);
// 
// template <typename T> void ompGetNeighFullTriplets(T *xij, T *xik, T *qi, T *qj, T *qk, T *x, T *q, int *ai, int *aj,  
//       int *ak, int *ti, int *tj, int *tk, int *anum, int *anumsum, int *ilist, 
//       int *alist, int *neighlist, int *neighnum, int *atomtype, int inum, int jnum, int ncq, int dim)
// {        
//     // form anum
//     for (int ii=0; ii<inum; ii++) {
//         int i = ilist[ii];       // atom i
//         int m = neighnum[i];     // number of neighbors around i             
//         anum[ii] = m*m;       
//     }
//     
//     // a list contains the starting positions of the first neighbor 
//     ompCumsum(anumsum, anum, inum+1);             
//     
//     for (int ii=0; ii<inum; ii++) {  // for each atom i in the simulation box     
//         int i = ilist[ii];       // atom i
//         int itype = atomtype[i];
//         int m = neighnum[i];     // number of neighbors around i             
//         int start = anumsum[ii];   
//         int inc = 0;
//         for (int lj=0; lj<m ; lj++) {   // loop over each atom around atom i
//             int gj = neighlist[lj +jnum*i];  // ghost index of atom j  
//             int j = alist[gj];  // atom j
//             for (int lk=0; lk<m; lk++) {
//                 int gk = neighlist[lk + jnum*i];  // atom k
//                 int k = alist[gk];  // atom k                
//                 int n = start + inc;                         
//                 for (int p=0; p<dim; p++) {
//                     xij[n*dim+p]   = x[gj*dim+p] - x[i*dim+p];  // xj - xi            
//                     xik[n*dim+p]   = x[gk*dim+p] - x[i*dim+p];  // xj - xi            
//                 }
//                 ai[n]        = i;
//                 aj[n]        = j; // should be alist[j];         
//                 ak[n]        = k; // should be alist[j];         
//                 ti[n]        = itype;       
//                 tj[n]        = atomtype[j]; // should be atomtype[alist[j]];                                        
//                 tk[n]        = atomtype[k]; // should be atomtype[alist[j]];                                        
//                 for (int p=0; p<ncq; p++) {                
//                     qi[n*ncq+p] = q[i*ncq+p];
//                     qj[n*ncq+p] = q[j*ncq+p];
//                     qk[n*ncq+p] = q[k*ncq+p];
//                 }                
//                 inc += 1;                
//             }                                    
//         }
//     }    
// }
// template void ompGetNeighFullTriplets(double*, double*, double*, double*, double*, double*, double*, 
//         int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int, int, int, int);
// template void ompGetNeighFullTriplets(float*, float*, float*, float*, float*, float*, float*, 
//         int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int, int, int, int);
// 
// template <typename T> void ompGetNeighTriplets(T *xij, T *xik, T *qi, T *qj, T *qk, T *x, T *q, int *ai, int *aj,  
//       int *ak, int *ti, int *tj, int *tk, int *anum, int *anumsum, int *ilist, int *alist, int *neighlist, int *neighnum, 
//       int *atomtype, int typej, int typek, int inum, int jnum, int ncq, int dim)
// {        
//     // form anum
//     for (int ii=0; ii<inum; ii++) {
//         int i = ilist[ii];       // atom i
//         int m = neighnum[i];     // number of neighbors around i             
//         anum[ii] = 0;              
//         for (int lj=0; lj<m ; lj++) {   // loop over each atom around atom i
//             int gj = neighlist[lj + jnum*i];  // atom j           
//             int j = alist[gj];  // atom j
//             for (int lk=0; lk<m; lk++) {
//                 int gk = neighlist[lk + jnum*i];  // atom k
//                 int k = alist[gk];  // atom k
//                 if ((atomtype[j] == typej) && (atomtype[k] == typek) && (j < k))
//                     anum[ii] += 1;
//             }            
//         }                
//     }
//     
//     // a list contains the starting positions of the first neighbor 
//     ompCumsum(anumsum, anum, inum+1);             
//     
//     for (int ii=0; ii<inum; ii++) {  // for each atom i in the simulation box     
//         int i = ilist[ii];       // atom i
//         int itype = atomtype[i];
//         int m = anum[ii];        // number of neighbors around i             
//         int start = anumsum[ii];   
//         int inc = 0;
//         for (int lj=0; lj<m ; lj++) {   // loop over each atom around atom i
//             int gj = neighlist[lj +jnum*i];  // ghost index of atom j  
//             int j = alist[gj];  // atom j
//             for (int lk=0; lk<m; lk++) {
//                 int gk = neighlist[lk + jnum*i];  // atom k
//                 int k = alist[gk];  // atom k
//                 if ((atomtype[j] == typej) && (atomtype[k] == typek) && (j < k)) {
//                     int n = start + inc;               
//                     for (int p=0; p<dim; p++) {
//                         xij[n*dim+p]   = x[gj*dim+p] - x[i*dim+p];  // xj - xi            
//                         xik[n*dim+p]   = x[gk*dim+p] - x[i*dim+p];  // xj - xi            
//                     }
//                     ai[n]        = i;
//                     aj[n]        = j; // should be alist[j];         
//                     ak[n]        = k; // should be alist[j];         
//                     ti[n]        = itype;       
//                     tj[n]        = atomtype[j]; // should be atomtype[alist[j]];                                        
//                     tk[n]        = atomtype[k]; // should be atomtype[alist[j]];                                        
//                     for (int p=0; p<ncq; p++) {                
//                         qi[n*ncq+p] = q[i*ncq+p];
//                         qj[n*ncq+p] = q[j*ncq+p];
//                         qk[n*ncq+p] = q[k*ncq+p];
//                     }                
//                     inc += 1;
//                 }
//             }                                    
//         }
//     }    
// }
// template void ompGetNeighTriplets(double*, double*, double*, double*, double*, double*, double*, 
//         int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int, int, int, int, int, int);
// template void ompGetNeighTriplets(float*, float*, float*, float*, float*, float*, float*, 
//         int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int, int, int, int, int, int);
// 
// template <typename T> void ompGetNeighTriplets(T *xij, T *xik, T *qi, T *qj, T *qk, T *x, T *q, int *ai, int *aj,  
//       int *ak, int *ti, int *tj, int *tk, int *anum, int *anumsum, int *ilist,  int *tlist, int *alist, int *neighlist, int *neighnum, 
//       int *atomtype, int *slist, int *start, int *index, int typei, int typej, int typek, int inum, int jnum, int ncq, int dim)
// {    
//     // form tlist from ilist
//     ompGetIlist(tlist, start, ilist, slist, index, atomtype, typei, inum);
//             
//     ompGetNeighTriplets(xij, xik, qi, qj,qk,  x, q, ai, aj, ak, ti, tj, tk, anum, anumsum, tlist, alist,
//             neighlist, neighnum, atomtype, typej, typek, start[1], jnum, ncq, dim);             
// }
// template void ompGetNeighTriplets(double*, double*, double*, double*, double*, double*, double*, 
//         int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, 
//         int, int, int, int, int, int, int);
// template void ompGetNeighTriplets(float*, float*, float*, float*, float*, float*, float*, 
//         int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, 
//         int, int, int, int, int, int, int);
// 
// template <typename T> void ompGetNeighFullTriplets(T *xij, T *xik, T *qi, T *qj, T *qk, T *x, T *q, int *ai, int *aj,  
//       int *ak, int *ti, int *tj, int *tk, int *anum, int *anumsum, int *ilist, int *alist, int *neighlist, int *neighnum, 
//       int *atomtype, int typej, int typek, int inum, int jnum, int ncq, int dim)
// {        
//     // form anum
//     for (int ii=0; ii<inum; ii++) {
//         int i = ilist[ii];       // atom i
//         int m = neighnum[i];     // number of neighbors around i             
//         anum[ii] = 0;              
//         for (int lj=0; lj<m ; lj++) {   // loop over each atom j around atom i
//             int gj = neighlist[lj + jnum*i];  // atom j           
//             int j = alist[gj];  // atom j
//             for (int lk=0; lk<m; lk++) {
//                 int gk = neighlist[lk + jnum*i];  // atom k
//                 int k = alist[gk];  // atom k
//                 if ((atomtype[j] == typej) && (atomtype[k] == typek))
//                     anum[ii] += 1;
//             }            
//         }                
//     }
//     
//     // a list contains the starting positions of the first neighbor 
//     ompCumsum(anumsum, anum, inum+1);             
//     
//     for (int ii=0; ii<inum; ii++) {  // for each atom i in the simulation box     
//         int i = ilist[ii];       // atom i
//         int itype = atomtype[i];
//         int m = anum[ii];        // number of neighbors around i             
//         int start = anumsum[ii];   
//         int inc = 0;
//         for (int lj=0; lj<m ; lj++) {   // loop over each atom around atom i
//             int gj = neighlist[lj +jnum*i];  // ghost index of atom j  
//             int j = alist[gj];  // atom j
//             for (int lk=0; lk<m; lk++) {
//                 int gk = neighlist[lk + jnum*i];  // atom k
//                 int k = alist[gk];  // atom k
//                 if ((atomtype[j] == typej) && (atomtype[k] == typek)) {
//                     int n = start + inc;               
//                     for (int p=0; p<dim; p++) {
//                         xij[n*dim+p]   = x[gj*dim+p] - x[i*dim+p];  // xj - xi            
//                         xik[n*dim+p]   = x[gk*dim+p] - x[i*dim+p];  // xj - xi            
//                     }
//                     ai[n]        = i;
//                     aj[n]        = j; // should be alist[j];         
//                     ak[n]        = k; // should be alist[j];         
//                     ti[n]        = itype;       
//                     tj[n]        = atomtype[j]; // should be atomtype[alist[j]];                                        
//                     tk[n]        = atomtype[k]; // should be atomtype[alist[j]];                                        
//                     for (int p=0; p<ncq; p++) {                
//                         qi[n*ncq+p] = q[i*ncq+p];
//                         qj[n*ncq+p] = q[j*ncq+p];
//                         qk[n*ncq+p] = q[k*ncq+p];
//                     }                
//                     inc += 1;
//                 }
//             }                                    
//         }
//     }    
// }
// template void ompGetNeighFullTriplets(double*, double*, double*, double*, double*, double*, double*, 
//         int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int, int, int, int, int, int);
// template void ompGetNeighFullTriplets(float*, float*, float*, float*, float*, float*, float*, 
//         int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int, int, int, int, int, int);
// 
// template <typename T> void ompGetNeighFullTriplets(T *xij, T *xik, T *qi, T *qj, T *qk, T *x, T *q, int *ai, int *aj,  
//       int *ak, int *ti, int *tj, int *tk, int *anum, int *anumsum, int *ilist,  int *tlist, int *alist, int *neighlist, int *neighnum, 
//       int *atomtype, int *slist, int *start, int *index, int typei, int typej, int typek, int inum, int jnum, int ncq, int dim)
// {    
//     // form tlist from ilist
//     ompGetIlist(tlist, start, ilist, slist, index, atomtype, typei, inum);
//             
//     ompGetNeighFullTriplets(xij, xik, qi, qj,qk,  x, q, ai, aj, ak, ti, tj, tk, anum, anumsum, tlist, alist,
//             neighlist, neighnum, atomtype, typej, typek, start[1], jnum, ncq, dim);             
// }
// template void ompGetNeighFullTriplets(double*, double*, double*, double*, double*, double*, double*, 
//         int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, 
//         int, int, int, int, int, int, int);
// template void ompGetNeighFullTriplets(float*, float*, float*, float*, float*, float*, float*, 
//         int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, 
//         int, int, int, int, int, int, int);
// 
// 
// template <typename T> void ompGetNeighFullTriplets(T *xij, T *xik, T *qi, T *qj, T *qk, T *x, T *q, int *ai, int *aj,  
//       int *ak, int *ti, int *tj, int *tk, int *anum, int *bnum, int *cnum, int *anumsum, int *ilist, int *alist, int *neighlist, int *neighnum, 
//       int *atomtype, int typej, int typek, int inum, int jnum, int ncq, int dim)
// {        
//     // form anum
//     for (int ii=0; ii<inum; ii++) {
//         int i = ilist[ii];       // atom i
//         int m = neighnum[i];     // number of neighbors around i             
//         anum[ii] = 0;                     
//         for (int lj=0; lj<m ; lj++) {   // loop over each atom j around atom i
//             int gj = neighlist[lj + jnum*i];  // atom j           
//             int j = alist[gj];  // atom j
//             for (int lk=0; lk<m; lk++) {
//                 int gk = neighlist[lk + jnum*i];  // atom k
//                 int k = alist[gk];  // atom k
//                 if ((atomtype[j] == typej) && (atomtype[k] == typek))
//                     anum[ii] += 1;
//             }            
//         }                
//         bnum[ii] = 0;             
//         cnum[ii] = 0;             
//         for (int lj=0; lj<m ; lj++) {   // loop over each atom j around atom i
//             int gj = neighlist[lj + jnum*i];  // atom j           
//             int j = alist[gj];  // atom j
//             if (atomtype[j] == typej)
//                 bnum[ii] += 1;
//             if (atomtype[j] == typek)
//                 cnum[ii] += 1;
//         }
//     }
//     
//     // a list contains the starting positions of the first neighbor 
//     ompCumsum(anumsum, anum, inum+1);             
//     
//     for (int ii=0; ii<inum; ii++) {  // for each atom i in the simulation box     
//         int i = ilist[ii];       // atom i
//         int itype = atomtype[i];
//         int m = anum[ii];        // number of neighbors around i             
//         int start = anumsum[ii];   
//         int inc = 0;
//         for (int lj=0; lj<m ; lj++) {   // loop over each atom around atom i
//             int gj = neighlist[lj +jnum*i];  // ghost index of atom j  
//             int j = alist[gj];  // atom j
//             for (int lk=0; lk<m; lk++) {
//                 int gk = neighlist[lk + jnum*i];  // atom k
//                 int k = alist[gk];  // atom k
//                 if ((atomtype[j] == typej) && (atomtype[k] == typek)) {
//                     int n = start + inc;               
//                     for (int p=0; p<dim; p++) {
//                         xij[n*dim+p]   = x[gj*dim+p] - x[i*dim+p];  // xj - xi            
//                         xik[n*dim+p]   = x[gk*dim+p] - x[i*dim+p];  // xj - xi            
//                     }
//                     ai[n]        = i;
//                     aj[n]        = j; // should be alist[j];         
//                     ak[n]        = k; // should be alist[j];         
//                     ti[n]        = itype;       
//                     tj[n]        = atomtype[j]; // should be atomtype[alist[j]];                                        
//                     tk[n]        = atomtype[k]; // should be atomtype[alist[j]];                                        
//                     for (int p=0; p<ncq; p++) {                
//                         qi[n*ncq+p] = q[i*ncq+p];
//                         qj[n*ncq+p] = q[j*ncq+p];
//                         qk[n*ncq+p] = q[k*ncq+p];
//                     }                
//                     inc += 1;
//                 }
//             }                                    
//         }
//     }    
// }
// template void ompGetNeighFullTriplets(double*, double*, double*, double*, double*, double*, double*, 
//         int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int, int, int, int, int, int);
// template void ompGetNeighFullTriplets(float*, float*, float*, float*, float*, float*, float*, 
//         int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int, int, int, int, int, int);
// 
// // template <typename T> void ompGetNeighPairs(T *xij, T *qi, T *qj, T *x, T *q, int *ai, int *aj,  
// //       int *ti, int *tj, int *alist, int *neighlist, int *neighnumsum, int *atomtype, 
// //       int istart, int iend, int jnum, int ncq, int dim)
// // {    
// //     for (int i=istart; i<iend; i++) {  // for each atom i in the simulation box     
// //         int itype = atomtype[i];        
// //         int m = neighnumsum[i+1] - neighnumsum[i];        // number of neighbors around i             
// //         int start = neighnumsum[i] - neighnumsum[istart];   
// //         for (int l=0; l<m ; l++) {   // loop over each atom around atom i
// //             int g = neighlist[l +jnum*i];  // ghost index of atom j  
// //             int j = alist[g];  // atom j
// //             int k = start + l;                                     
// //             ai[k]        = i;
// //             aj[k]        = j; // should be alist[j];         
// //             ti[k]        = itype;       
// //             tj[k]        = atomtype[j]; // should be atomtype[alist[j]];         
// //             for (int p=0; p<dim; p++) 
// //                 xij[k*dim+p]   = x[g*dim+p] - x[i*dim+p];  // xj - xi            
// //             for (int p=0; p<ncq; p++) {                
// //                 qi[k*ncq+p] = q[i*ncq+p];
// //                 qj[k*ncq+p] = q[j*ncq+p];
// //             }                
// //         }
// //     }    
// // }
// // template void ompGetNeighPairs(double*, double*, double*, double*, double*, int*, int*, int*, int*, int*, int*, int*, int*, int, int, int, int, int);
// // template void ompGetNeighPairs(float*, float*, float*, float*, float*, int*, int*, int*, int*, int*, int*, int*, int*, int, int, int, int, int);
// // 
// // template <typename T> void ompGetNeighPairs(T *xij, T *qi, T *qj, T *x, T *q, int *ai, int *aj,  
// //       int *ti, int *tj, int *anum, int *anumsum, int *alist, int *neighlist, int *neighnumsum, 
// //       int *atomtype, int typej, int istart, int iend, int jnum, int ncq, int dim)
// // {        
// //     // form anum
// //     for (int i=istart; i<iend; i++) {       
// //         int m = neighnumsum[i+1] - neighnumsum[i];        // number of neighbors around i     
// //         anum[i-istart] = 0;              
// //         for (int l=0; l<m ; l++) {   // loop over each atom around atom i
// //             int g = neighlist[l + jnum*i];  // atom j           
// //             int j = alist[g];  // atom j
// //             if (atomtype[j] == typej) 
// //                 anum[i-istart] += 1;
// //         }                
// //     }
// //     
// //     // a list contains the starting positions of the first neighbor 
// //     int inum = iend - istart; // number of atoms i
// //     ompCumsum(anumsum, anum, inum+1);             
// //     
// //     for (int i=istart; i<iend; i++) {  // for each atom i in the simulation box     
// //         int itype = atomtype[i];
// //         int m = anum[i-istart];        // number of neighbors around i             
// //         int start = anumsum[i-istart];   
// //         int inc = 0;
// //         for (int l=0; l<m ; l++) {   // loop over each atom around atom i
// //             int g = neighlist[l +jnum*i];  // ghost index of atom j  
// //             int j = alist[g];  // atom j
// //             if (atomtype[j] == typej)  {
// //                 int k = start + inc;                         
// //                 ai[k]        = i;
// //                 aj[k]        = j; // should be alist[j];         
// //                 ti[k]        = itype;       
// //                 tj[k]        = atomtype[j]; // should be atomtype[alist[j]];                                        
// //                 for (int p=0; p<dim; p++) 
// //                     xij[k*dim+p]   = x[g*dim+p] - x[i*dim+p];  // xj - xi                            
// //                 for (int p=0; p<ncq; p++) {                
// //                     qi[k*ncq+p] = q[i*ncq+p];
// //                     qj[k*ncq+p] = q[j*ncq+p];
// //                 }                
// //                 inc += 1;
// //             }
// //         }
// //     }    
// // }
// // template void ompGetNeighPairs(double*, double*, double*, double*, double*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int, int, int, int, int, int);
// // template void ompGetNeighPairs(float*, float*, float*, float*, float*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int, int, int, int, int, int);
// // 
// // template <typename T> void ompGetNeighTriplets(T *xij, T *xik, T *qi, T *qj, T *qk, T *x, T *q, int *ai, int *aj,  
// //       int *ak, int *ti, int *tj, int *tk, int *anum, int *anumsum, int *alist, int *neighlist, int *neighnumsum, 
// //       int *atomtype, int istart, int iend, int jnum, int ncq, int dim)
// // {        
// //     // form anum
// //     for (int i=istart; i<iend; i++) {  // for each atom i in the simulation box     
// //         int m = neighnumsum[i+1] - neighnumsum[i];        // number of neighbors around i             
// //         anum[i] = 0;              
// //         for (int lj=0; lj<m ; lj++) {   // loop over each atom around atom i
// //             int gj = neighlist[lj + jnum*i];  // atom j           
// //             int j = alist[gj];  // atom j
// //             for (int lk=0; lk<m; lk++) {
// //                 int gk = neighlist[lk + jnum*i];  // atom k
// //                 int k = alist[gk];  // atom k
// //                 if (j < k)
// //                     anum[i] += 1;
// //             }            
// //         }                
// //     }
// //     
// //     // a list contains the starting positions of the first neighbor 
// //     int inum = iend - istart;
// //     ompCumsum(anumsum, anum, inum+1);             
// //     
// //     for (int i=istart; i<iend; i++) {  // for each atom i in the simulation box           
// //         int itype = atomtype[i];
// //         int m = anum[i-istart];        // number of neighbors around i             
// //         int start = anumsum[i-istart];   
// //         int inc = 0;
// //         for (int lj=0; lj<m ; lj++) {   // loop over each atom around atom i
// //             int gj = neighlist[lj +jnum*i];  // ghost index of atom j  
// //             int j = alist[gj];  // atom j
// //             for (int lk=0; lk<m; lk++) {
// //                 int gk = neighlist[lk + jnum*i];  // atom k
// //                 int k = alist[gk];  // atom k
// //                 if (j < k) {
// //                     int n = start + inc;                         
// //                     for (int p=0; p<dim; p++) {
// //                         xij[n*dim+p]   = x[gj*dim+p] - x[i*dim+p];  // xj - xi            
// //                         xik[n*dim+p]   = x[gk*dim+p] - x[i*dim+p];  // xj - xi            
// //                     }
// //                     ai[n]        = i;
// //                     aj[n]        = j; // should be alist[j];         
// //                     ak[n]        = k; // should be alist[j];         
// //                     ti[n]        = itype;       
// //                     tj[n]        = atomtype[j]; // should be atomtype[alist[j]];                                        
// //                     tk[n]        = atomtype[k]; // should be atomtype[alist[j]];                                        
// //                     for (int p=0; p<ncq; p++) {                
// //                         qi[n*ncq+p] = q[i*ncq+p];
// //                         qj[n*ncq+p] = q[j*ncq+p];
// //                         qk[n*ncq+p] = q[k*ncq+p];
// //                     }                
// //                     inc += 1;
// //                 }
// //             }                                    
// //         }
// //     }    
// // }
// // template void ompGetNeighTriplets(double*, double*, double*, double*, double*, double*, double*, 
// //         int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int, int, int, int, int);
// // template void ompGetNeighTriplets(float*, float*, float*, float*, float*, float*, float*, 
// //         int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int, int, int, int, int);
// 
// //********************************************************************************************//
// 
// // template <typename T> void ompFullAtomDecomposition2D(T *fi, T *fij, int *neighnumsum, int istart, int iend)
// // {        
// //     for (int i=istart; i<iend; i++) {  // for each atom i in the simulation box     
// //         int m = neighnumsum[i+1] - neighnumsum[i];        // number of neighbors around i             
// //         int start = neighnumsum[i] - neighnumsum[istart];       
// //         T f1 = 0.0;
// //         T f2 = 0.0;
// //         for (int l=0; l<m ; l++) {   // loop over each atom j around atom i
// //             int k = start + l;                     
// //             f1 +=  -fij[2*k+0];
// //             f2 +=  -fij[2*k+1];
// //         }
// //         fi[2*i+0] = f1;
// //         fi[2*i+1] = f2;
// //     }
// // }
// // template void ompFullAtomDecomposition2D(double*, double*, int*, int, int);
// // template void ompFullAtomDecomposition2D(float*, float*, int*, int, int);
// // 
// // template <typename T> void ompHalfAtomDecomposition2D(T *fi, T *fij, int *neighnumsum, 
// //         int *bnumsum, int *index, int istart, int iend)
// // {    
// //     for (int i=istart; i<iend; i++) {  // for each atom i in the simulation box     
// //         int m = neighnumsum[i+1] - neighnumsum[i];        // number of neighbors around i             
// //         int start = neighnumsum[i] - neighnumsum[istart];               
// //         T f1 = 0.0;
// //         T f2 = 0.0;        
// //         for (int l=0; l<m ; l++) {   // loop over each atom around j atom i -> pair ij 
// //             int k = start + l;                     
// //             f1 +=  -fij[2*k+0];
// //             f2 +=  -fij[2*k+1];
// //         }                
// //         
// //         int ii = i - istart;
// //         start = bnumsum[ii];   
// //         m = bnumsum[ii+1]-start;   // number of neighbors j around i  (j < i)           
// //         for (int l=0; l<m ; l++) { // loop over each atom around j atom i -> pair ji 
// //             int k = index[start + l];                     
// //             f1 +=  fij[2*k+0];
// //             f2 +=  fij[2*k+1];
// //         }
// //         fi[2*i+0] = f1;
// //         fi[2*i+1] = f2;
// //     }
// // }
// // template void ompHalfAtomDecomposition2D(double*, double*, int*, int*, int*, int, int);
// // template void ompHalfAtomDecomposition2D(float*, float*, int*, int*, int*, int, int);
// 
// // template <typename T> void ompHalfAtomDecomposition2D(T *f, T *fij, int *ilist, int *anumsum, 
// //         int *jlist, int *bnumsum, int *index, int inum, int jnum)
// // {    
// //     ompIAtomDecomposition2D(f, fij, ilist, anumsum, inum);
// //     ompJAtomDecomposition2D(f, fij, jlist, bnumsum, index, jnum);   
// // }
// // template void ompHalfAtomDecomposition2D(double*, double*, int*, int*, int*, int*, int*, int, int);
// // template void ompHalfAtomDecomposition2D(float*, float*, int*, int*, int*, int*, int*, int, int);
// 
// // 
// // template <typename T> void ompFullAtomDecomposition3D(T *fi, T *fij, int *neighnumsum, int istart, int iend)
// // {        
// //     for (int i=istart; i<iend; i++) {  // for each atom i in the simulation box     
// //         int m = neighnumsum[i+1] - neighnumsum[i];        // number of neighbors around i             
// //         int start = neighnumsum[i] - neighnumsum[istart];       
// //         T f1 = 0.0;
// //         T f2 = 0.0;
// //         T f3 = 0.0;        
// //         for (int l=0; l<m ; l++) {   // loop over each atom around atom i
// //             int k = start + l;                     
// //             f1 +=  -fij[3*k+0];
// //             f2 +=  -fij[3*k+1];
// //             f3 +=  -fij[3*k+2];             
// //         }
// //         fi[3*i+0] = f1;
// //         fi[3*i+1] = f2;
// //         fi[3*i+2] = f3;        
// //     }
// // }
// // template void ompFullAtomDecomposition3D(double*, double*, int*, int, int);
// // template void ompFullAtomDecomposition3D(float*, float*, int*, int, int);
// // 
// // template <typename T> void ompHalfAtomDecomposition3D(T *fi, T *fij, int *neighnumsum, 
// //         int *bnumsum, int *index, int istart, int iend)
// // {    
// //     for (int i=istart; i<iend; i++) {  // for each atom i in the simulation box     
// //         int m = neighnumsum[i+1] - neighnumsum[i];        // number of neighbors around i             
// //         int start = neighnumsum[i] - neighnumsum[istart];               
// //         T f1 = 0.0;
// //         T f2 = 0.0;        
// //         T f3 = 0.0;        
// //         for (int l=0; l<m ; l++) {   // loop over each atom around j atom i -> pair ij 
// //             int k = start + l;             
// //             f1 +=  -fij[3*k+0];
// //             f2 +=  -fij[3*k+1];
// //             f3 +=  -fij[3*k+2];                         
// //         }                
// //         
// //         // need to determine bnumsum and index 
// //         int ii = i - istart;
// //         start = bnumsum[ii];   
// //         m = bnumsum[ii+1]-start;   // number of neighbors j around i  (j < i)           
// //         for (int l=0; l<m ; l++) { // loop over each atom around j atom i -> pair ji 
// //             int k = index[start + l];     
// //             f1 +=  fij[3*k+0];
// //             f2 +=  fij[3*k+1];
// //             f3 +=  fij[3*k+2];             
// //         }
// //         fi[3*i+0] = f1;
// //         fi[3*i+1] = f2;
// //         fi[3*i+2] = f3;        
// //     }
// // }
// // template void ompHalfAtomDecomposition3D(double*, double*, int*, int*, int*, int, int);
// // template void ompHalfAtomDecomposition3D(float*, float*, int*, int*, int*, int, int);
// 
// // template <typename T> void ompHalfAtomDecomposition3D(T *fi, T *fij, int *ilist, int *anumsum, 
// //         int *bnumsum, int *index, int inum)
// // {    
// //     for (int ii=0; ii<inum; ii++) {  // for each atom i in the simulation box     
// //         int i = ilist[ii];       // atom i
// //         T f1 = 0.0;
// //         T f2 = 0.0;        
// //         T f3 = 0.0;        
// //         int start = anumsum[ii];   
// //         int m = anumsum[ii+1]-start; // number of neighbors j around i  (j > i)             
// //         for (int l=0; l<m ; l++) {   // loop over each atom around j atom i -> pair ij 
// //             int k = start + l;                     
// //             f1 +=  -fij[3*k+0];
// //             f2 +=  -fij[3*k+1];
// //             f3 +=  -fij[3*k+2];             
// //         }                
// //         // need to determine bnumsum and index
// //         start = bnumsum[ii];   
// //         m = bnumsum[ii+1]-start;   // number of neighbors j around i  (j < i)           
// //         for (int l=0; l<m ; l++) { // loop over each atom around j atom i -> pair ji 
// //             int k = index[start + l];                     
// //             f1 +=  fij[3*k+0];
// //             f2 +=  fij[3*k+1];
// //             f3 +=  fij[3*k+2];             
// //         }
// //         fi[3*i+0] = f1;
// //         fi[3*i+1] = f2;
// //         fi[3*i+2] = f3;
// //     }
// // }
// // template void ompHalfAtomDecomposition3D(double*, double*, int*, int*, int*, int*, int);
// // template void ompHalfAtomDecomposition3D(float*, float*, int*, int*, int*, int*, int);


#endif


