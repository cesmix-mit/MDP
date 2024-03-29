/***************************************************************************                               
                    Molecular Dynamics Potentials (MDP)
                           CESMIX-MIT Project  
 
 Contributing authors: Ngoc-Cuong Nguyen (cuongng@mit.edu, exapde@gmail.com)
 ***************************************************************************/

#ifndef __CPUENERGYFORCETALLY
#define __CPUENERGYFORCETALLY

template <typename T> void cpuSingleDecomposition(T *e, T *ei, int *ai, int inum, int dim)
{    
    for (int ii=0; ii<inum; ii++) {  // for each atom pair ij in the simulation box     
        int i = ai[ii];       // atom i        
        e[i] += ei[ii];
    }
}
template void cpuSingleDecomposition(double*, double*, int*, int, int);
template void cpuSingleDecomposition(float*, float*, int*, int, int);

template <typename T> void cpuFullNeighPairDecomposition(T *e, T *eij, int *ai, int ijnum, int dim)
{    
    for (int ii=0; ii<ijnum; ii++) {  // for each atom pair ij in the simulation box     
        int i = ai[ii];       // atom i
        // use atomicAdd on gpu
        e[i] += 0.5*eij[ii];
    }
}
template void cpuFullNeighPairDecomposition(double*, double*, int*, int, int);
template void cpuFullNeighPairDecomposition(float*, float*, int*, int, int);

template <typename T> void cpuHalfNeighPairDecomposition(T *e, T *eij, int *ai, int *aj, int ijnum, int dim)
{    
    for (int ii=0; ii<ijnum; ii++) {  // for each atom pair ij in the simulation box     
        int i = ai[ii];       // atom i
        int j = aj[ii];       // atom j
        // use atomicAdd on gpu
        e[i] += 0.5*eij[ii];
        e[j] += 0.5*eij[ii];
    }
}
template void cpuHalfNeighPairDecomposition(double*, double*, int*, int*, int, int);
template void cpuHalfNeighPairDecomposition(float*, float*, int*, int*, int, int);

template <typename T> void cpuCenterAtomPairDecomposition(T *e, T *eij, int *ilist, int *anumsum, int inum, int dim)
{    
    for (int ii=0; ii<inum; ii++) {  // for each atom i in the simulation box     
        int i = ilist[ii];       // atom i        
        int start = anumsum[ii];   
        int m = anumsum[ii+1]-start;        // number of neighbors around i             
        for (int l=0; l<m ; l++)    // loop over each atom around atom i            
            e[i] += 0.5*eij[start + l];                                
    }
}
template void cpuCenterAtomPairDecomposition(double*, double*, int*, int*, int, int);
template void cpuCenterAtomPairDecomposition(float*, float*, int*, int*, int, int);

template <typename T> void cpuNeighborAtomPairDecomposition(T *e, T *eij, int *jlist, int *bnumsum, int *index, int jnum, int dim)
{        
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
template void cpuNeighborAtomPairDecomposition(double*, double*, int*, int*, int*, int, int);
template void cpuNeighborAtomPairDecomposition(float*, float*, int*, int*, int*, int, int);

template <typename T> void cpuTripletDecomposition(T *e, T *eijk, int *ai, int *aj, int *ak, int ijknum, int dim)
{    
    for (int ii=0; ii<ijknum; ii++) {  // for each atom pair ij in the simulation box     
        int i = ai[ii];       // atom i
        int j = aj[ii];       // atom j
        int k = ak[ii];       // atom k
        // use atomicAdd on gpu
        T ee = eijk[ii]/3.0;
        e[i] += ee;
        e[j] += ee;
        e[k] += ee;
    }
}
template void cpuTripletDecomposition(double*, double*, int*, int*, int*, int, int);
template void cpuTripletDecomposition(float*, float*, int*, int*, int*, int, int);

template <typename T> void cpuCenterAtomTripletDecomposition(T *e, T *eijk, int *ilist, int *anumsum, int inum, int dim)
{    
    for (int ii=0; ii<inum; ii++) {  // for each atom i in the simulation box     
        int i = ilist[ii];       // atom i        
        int start = anumsum[ii];   
        int m = anumsum[ii+1]-start;        // number of neighbors around i             
        for (int l=0; l<m ; l++)    // loop over each atom around atom i            
            e[i] += eijk[start + l]/3.0;                                
    }
}
template void cpuCenterAtomTripletDecomposition(double*, double*, int*, int*, int, int);
template void cpuCenterAtomTripletDecomposition(float*, float*,  int*, int*, int, int);

template <typename T> void cpuNeighborAtomTripletDecomposition(T *e, T *eij, int *jlist, int *bnumsum, int *index, int jnum, int dim)
{        
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
template void cpuNeighborAtomTripletDecomposition(double*, double*, int*, int*, int*, int, int);
template void cpuNeighborAtomTripletDecomposition(float*, float*, int*, int*, int*, int, int);

template <typename T> void cpuQuadrupletDecomposition(T *e, T *eijkl, 
        int *ai, int *aj, int *ak, int *al, int ijklnum, int dim)
{    
    for (int ii=0; ii<ijklnum; ii++) {  // for each atom quadruplet ijkl in the simulation box     
        int i = ai[ii];       // atom i
        int j = aj[ii];       // atom j
        int k = ak[ii];       // atom k
        int l = al[ii];       // atom l
        // use atomicAdd on gpu
        T ee = eijkl[ii]/4.0;
        e[i] += ee;
        e[j] += ee;
        e[k] += ee;
        e[l] += ee;
    }
}
template void cpuQuadrupletDecomposition(double*, double*, int*, int*, int*, int*, int, int);
template void cpuQuadrupletDecomposition(float*, float*, int*, int*, int*, int*, int, int);

template <typename T> void cpuCenterAtomQuadrupletDecomposition(T *e, T *eijkl, int *ilist, int *anumsum, int inum, int dim)
{    
    for (int ii=0; ii<inum; ii++) {  // for each atom i in the simulation box     
        int i = ilist[ii];       // atom i        
        int start = anumsum[ii];   
        int m = anumsum[ii+1]-start;        // number of neighbors around i             
        for (int l=0; l<m ; l++)    // loop over each atom around atom i            
            e[i] += eijkl[start + l]/4.0;                                
    }
}
template void cpuCenterAtomQuadrupletDecomposition(double*, double*, int*, int*, int, int);
template void cpuCenterAtomQuadrupletDecomposition(float*, float*, int*, int*, int, int);

template <typename T> void cpuNeighborAtomQuadrupletDecomposition(T *e, T *eij, int *jlist, int *bnumsum, int *index, int jnum, int dim)
{        
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
template void cpuNeighborAtomQuadrupletDecomposition(double*, double*, int*, int*, int*, int, int);
template void cpuNeighborAtomQuadrupletDecomposition(float*, float*, int*, int*, int*, int, int);

template <typename T> void cpuSingleDecomposition(T *e, T *f, T *ei, T *fi, int *ai, int inum, int dim)
{    
    for (int ii=0; ii<inum; ii++) {  // for each atom pair ij in the simulation box     
        int i = ai[ii];       // atom i        
        int l = dim*ii;                
        e[i] += ei[ii];
        i = dim*i;
        for (int d=0; d<dim; d++) 
            f[i+d] +=  fi[l+d];        
    }
}
template void cpuSingleDecomposition(double*, double*, double*, double*, int*, int, int);
template void cpuSingleDecomposition(float*, float*, float*, float*, int*, int, int);

template <typename T> void cpuCenterAtomPairDecomposition(T *e, T *f, T *eij, T *fij, int *ilist, int *anumsum, int inum, int dim)
{    
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
template void cpuCenterAtomPairDecomposition(double*, double*, double*, double*, int*, int*, int, int);
template void cpuCenterAtomPairDecomposition(float*, float*, float*, float*, int*, int*, int, int);

template <typename T> void cpuNeighborAtomPairDecomposition(T *e, T *f, T *eij, T *fij, int *jlist, int *bnumsum, int *index, int jnum, int dim)
{        
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
template void cpuNeighborAtomPairDecomposition(double*, double*, double*, double*, int*, int*, int*, int, int);
template void cpuNeighborAtomPairDecomposition(float*, float*, float*, float*, int*, int*, int*, int, int);

template <typename T> void cpuFullNeighPairDecomposition(T *e, T *f, T *eij, T *fij, int *ai, int ijnum, int dim)
{    
    for (int ii=0; ii<ijnum; ii++) {  // for each atom pair ij in the simulation box     
        int i = ai[ii];       // atom i
        int l = dim*ii;
        // use atomicAdd on gpu
        e[i] += 0.5*eij[ii];
        i = dim*i;
        for (int d=0; d<dim; d++) 
            f[i+d] +=  -fij[l+d];        
    }
}
template void cpuFullNeighPairDecomposition(double*, double*, double*, double*, int*, int, int);
template void cpuFullNeighPairDecomposition(float*, float*, float*, float*, int*, int, int);

template <typename T> void cpuHalfNeighPairDecomposition(T *e, T *f, T *eij, T *fij, int *ai, int *aj, int ijnum, int dim)
{    
    for (int ii=0; ii<ijnum; ii++) {  // for each atom pair ij in the simulation box     
        int i = ai[ii];       // atom i
        int j = aj[ii];       // atom j
        int l = dim*ii;        
        // use atomicAdd on gpu
        e[i] += 0.5*eij[ii];
        e[j] += 0.5*eij[ii];
        i = dim*i;
        j = dim*j;
        for (int d=0; d<dim; d++) {
            f[i+d] += -fij[l+d]; 
            f[j+d] +=  fij[l+d]; 
        }        
    }
}
template void cpuHalfNeighPairDecomposition(double*, double*, double*, double*, int*, int*, int, int);
template void cpuHalfNeighPairDecomposition(float*, float*, float*, float*, int*, int*, int, int);

template <typename T> void cpuTripletDecomposition(T *e, T *f, T *eijk, T *fij, T *fik, int *ai, int *aj, int *ak, int ijknum, int dim)
{    
    for (int ii=0; ii<ijknum; ii++) {  // for each atom pair ij in the simulation box     
        int i = ai[ii];       // atom i
        int j = aj[ii];       // atom j
        int k = ak[ii];       // atom k
        int n = dim*ii;        
        // use atomicAdd on gpu
        T ee = eijk[ii]/3.0;
        e[i] += ee;
        e[j] += ee;
        e[k] += ee;
        i = dim*i;
        j = dim*j;
        k = dim*k;
        for (int d=0; d<dim; d++) {
            f[i+d] += -(fij[n+d] + fik[n+d]); 
            f[j+d] +=  fij[n+d]; 
            f[k+d] +=  fik[n+d]; 
        }        
    }
}
template void cpuTripletDecomposition(double*, double*, double*, double*, double*, int*, int*, int*, int, int);
template void cpuTripletDecomposition(float*, float*, float*, float*, float*, int*, int*, int*, int, int);

template <typename T> void cpuCenterAtomTripletDecomposition(T *e, T *f, T *eijk, T *fij, T *fik, int *ilist, int *anumsum, int inum, int dim)
{    
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
template void cpuCenterAtomTripletDecomposition(double*, double*, double*, double*, double*, int*, int*, int, int);
template void cpuCenterAtomTripletDecomposition(float*, float*, float*, float*, float*, int*, int*, int, int);

template <typename T> void cpuNeighborAtomTripletDecomposition(T *e, T *f, T *eij, T *fij, int *jlist, int *bnumsum, int *index, int jnum, int dim)
{        
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
template void cpuNeighborAtomTripletDecomposition(double*, double*, double*, double*, int*, int*, int*, int, int);
template void cpuNeighborAtomTripletDecomposition(float*, float*, float*, float*, int*, int*, int*, int, int);

template <typename T> void cpuQuadrupletDecomposition(T *e, T *f, T *eijkl, T *fij, T *fik, T *fil, 
        int *ai, int *aj, int *ak, int *al, int ijklnum, int dim)
{    
    for (int ii=0; ii<ijklnum; ii++) {  // for each atom quadruplet ijkl in the simulation box     
        int i = ai[ii];       // atom i
        int j = aj[ii];       // atom j
        int k = ak[ii];       // atom k
        int l = al[ii];       // atom l
        int n = dim*ii;        
        // use atomicAdd on gpu
        T ee = eijkl[ii]/4.0;
        e[i] += ee;
        e[j] += ee;
        e[k] += ee;
        e[l] += ee;
        i = dim*i;
        j = dim*j;
        k = dim*k;
        l = dim*l;
        for (int d=0; d<dim; d++) {
            f[i+d] += -(fij[n+d] + fik[n+d] + fil[n+d]); 
            f[j+d] +=  fij[n+d]; 
            f[k+d] +=  fik[n+d]; 
            f[l+d] +=  fil[n+d]; 
        }        
    }
}
template void cpuQuadrupletDecomposition(double*, double*, double*, double*, double*, double*, int*, int*, int*, int*, int, int);
template void cpuQuadrupletDecomposition(float*, float*, float*, float*, float*, float*, int*, int*, int*, int*, int, int);

template <typename T> void cpuCenterAtomQuadrupletDecomposition(T *e, T *f, T *eijkl, T *fij, T *fik, T *fil, 
        int *ilist, int *anumsum, int inum, int dim)
{    
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
template void cpuCenterAtomQuadrupletDecomposition(double*, double*, double*, double*, double*, double*, int*, int*, int, int);
template void cpuCenterAtomQuadrupletDecomposition(float*, float*, float*, float*, float*, float*, int*, int*, int, int);

template <typename T> void cpuNeighborAtomQuadrupletDecomposition(T *e, T *f, T *eij, T *fij, int *jlist, int *bnumsum, int *index, int jnum, int dim)
{        
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
template void cpuNeighborAtomQuadrupletDecomposition(double*, double*, double*, double*, int*, int*, int*, int, int);
template void cpuNeighborAtomQuadrupletDecomposition(float*, float*, float*, float*, int*, int*, int*, int, int);

template <typename T> void cpuHalfForceDecomposition(T *f, T *fij, int *ai, int *aj, int dim, int ijnum)
{    
    for (int ii=0; ii<ijnum; ii++) {  // for each atom pair ij in the simulation box     
        int i = dim*ai[ii];       // atom i
        int j = dim*aj[ii];       // atom j
        int l = dim*ii;
        for (int d = 0; d<dim; d++) {
            f[i+d] -= fij[l+d];
            f[j+d] += fij[l+d];
        }
    }
}
template void cpuHalfForceDecomposition(double*, double*, int*, int*, int, int);
template void cpuHalfForceDecomposition(float*, float*, int*, int*, int, int);

template <typename T> void cpuIAtomDecomposition(T *f, T *fij, int *ilist, int *anumsum, int dim, int inum)
{    
    for (int ii=0; ii<inum; ii++) {  // for each atom i in the simulation box     
        int i = dim*ilist[ii];       // atom i        
        int start = anumsum[ii];   
        int m = anumsum[ii+1]-start;        // number of neighbors around i                     
        for (int l=0; l<m ; l++) {   // loop over each atom around atom i
            int k = dim*(start + l); 
            for (int d=0; d<dim; d++)
                f[i+d] -= fij[k+d];
        }
    }
}
template void cpuIAtomDecomposition(double*, double*, int*, int*, int, int);
template void cpuIAtomDecomposition(float*, float*, int*, int*, int, int);

template <typename T> void cpuJAtomDecomposition(T *f, T *fij, int *jlist, int *bnumsum, int *index, int dim, int jnum)
{        
    for (int ii=0; ii<jnum; ii++) {  // for each atom j in the simulation box     
        int j = dim*jlist[ii];       // atom j
        // need to determine bnumsum and index 
        int start = bnumsum[ii];   
        int m = bnumsum[ii+1]-start;   // number of neighbors i around j  (j < i)           
        for (int l=0; l<m ; l++) { // loop over each atom around i atom j -> pair ij 
            int k = dim*index[start + l];        
            for (int d=0; d<dim; d++)
                f[j+d] += fij[k+d];            
        }
    }
}
template void cpuJAtomDecomposition(double*, double*, int*, int*, int*, int, int);
template void cpuJAtomDecomposition(float*, float*, int*, int*, int*, int, int);


// template <typename T> void cpuSingleDecomposition2D(T *f, T *fi, int *ai, int inum)
// {    
//     for (int ii=0; ii<inum; ii++) {  // for each atom pair ij in the simulation box     
//         int i = 2*ai[ii];       // atom i
//         int l = 2*ii;
//         // use atomicAdd on gpu
//         f[i+0] +=  fi[l+0];
//         f[i+1] +=  fi[l+1];
//     }
// }
// template void cpuSingleDecomposition2D(double*, double*, int*, int);
// template void cpuSingleDecomposition2D(float*, float*, int*, int);
// 
// template <typename T> void cpuSingleDecomposition3D(T *f, T *fi, int *ai, int inum)
// {    
//     for (int ii=0; ii<inum; ii++) {  // for each atom pair ij in the simulation box     
//         int i = 3*ai[ii];       // atom i
//         int l = 3*ii;
//         // use atomicAdd on gpu
//         f[i+0] +=  fi[l+0];
//         f[i+1] +=  fi[l+1];
//         f[i+2] +=  fi[l+2];
//     }
// }
// template void cpuSingleDecomposition3D(double*, double*, int*, int);
// template void cpuSingleDecomposition3D(float*, float*, int*, int);
// 
// 
// template <typename T> void cpuHalfForceDecomposition2D(T *f, T *fij, int *ai, int *aj, int ijnum)
// {    
//     for (int ii=0; ii<ijnum; ii++) {  // for each atom pair ij in the simulation box     
//         int i = 2*ai[ii];       // atom i
//         int j = 2*aj[ii];       // atom j
//         int l = 2*ii;
//         // use atomicAdd on gpu
//         f[i+0] +=  -fij[l+0];
//         f[i+1] +=  -fij[l+1];
//         f[j+0] +=  fij[l+0];
//         f[j+1] +=  fij[l+1];
//     }
// }
// template void cpuHalfForceDecomposition2D(double*, double*, int*, int*, int);
// template void cpuHalfForceDecomposition2D(float*, float*, int*, int*, int);
// 
// template <typename T> void cpuFullForceDecomposition2D(T *f, T *fij, int *ai, int *aj, int ijnum)
// {    
//     for (int ii=0; ii<ijnum; ii++) {  // for each atom pair ij in the simulation box     
//         int i = 2*ai[ii];       // atom i
//         int l = 2*ii;
//         // use atomicAdd on gpu
//         f[i+0] +=  -fij[l+0];
//         f[i+1] +=  -fij[l+1];
//     }
// }
// template void cpuFullForceDecomposition2D(double*, double*, int*, int*, int);
// template void cpuFullForceDecomposition2D(float*, float*, int*, int*, int);
// 
// template <typename T> void cpuIAtomDecomposition2D(T *f, T *fij, int *ilist, int *anumsum, int inum)
// {    
//     for (int ii=0; ii<inum; ii++) {  // for each atom i in the simulation box     
//         int i = 2*ilist[ii];       // atom i        
//         int start = anumsum[ii];   
//         int m = anumsum[ii+1]-start;        // number of neighbors around i             
//         T f1 = 0.0;
//         T f2 = 0.0;
//         for (int l=0; l<m ; l++) {   // loop over each atom around atom i
//             int k = 2*(start + l);                     
//             f1 +=  -fij[k+0];
//             f2 +=  -fij[k+1];
//         }
//         f[i+0] = f1;
//         f[i+1] = f2;
//     }
// }
// template void cpuIAtomDecomposition2D(double*, double*, int*, int*, int);
// template void cpuIAtomDecomposition2D(float*, float*, int*, int*, int);
// 
// template <typename T> void cpuJAtomDecomposition2D(T *f, T *fij, int *jlist, int *bnumsum, int *index, int jnum)
// {        
//     for (int ii=0; ii<jnum; ii++) {  // for each atom j in the simulation box     
//         int j = 2*jlist[ii];       // atom j
//         T f1 = 0.0;
//         T f2 = 0.0;                
//         // need to determine bnumsum and index 
//         int start = bnumsum[ii];   
//         int m = bnumsum[ii+1]-start;   // number of neighbors i around j  (j < i)           
//         for (int l=0; l<m ; l++) { // loop over each atom around i atom j -> pair ij 
//             int k = 2*index[start + l];                     
//             f1 +=  fij[k+0];
//             f2 +=  fij[k+1];
//         }
//         f[j+0] = f1;
//         f[j+1] = f2;
//     }
// }
// template void cpuJAtomDecomposition2D(double*, double*, int*, int*, int*, int);
// template void cpuJAtomDecomposition2D(float*, float*, int*, int*, int*, int);

// template <typename T> void cpuForceDecompositionTriplet2D(T *f, T *fij, T *fik, int *ai, int *aj, int *ak, int ijknum)
// {    
//     for (int ii=0; ii<ijknum; ii++) {  // for each atom pair ij in the simulation box     
//         int i = 2*ai[ii];       // atom i
//         int j = 2*aj[ii];       // atom j
//         int k = 2*ak[ii];       // atom k
//         int l = 2*ii;
//         // use atomicAdd on gpu
//         f[i+0] +=  -(fij[l+0]+fik[l+0]);
//         f[i+1] +=  -(fij[l+1]+fik[l+1]);
//         f[j+0] +=  fij[l+0];
//         f[j+1] +=  fij[l+1];
//         f[k+0] +=  fik[l+0];
//         f[k+1] +=  fik[l+1];
//     }
// }
// template void cpuForceDecompositionTriplet2D(double*, double*, double*, int*, int*, int*, int);
// template void cpuForceDecompositionTriplet2D(float*, float*, float*, int*, int*, int*, int);
// 
// template <typename T> void cpuAtomDecompositionTriplet2D(T *f, T *fij, T *fik, int *ilist, int *anumsum, int inum)
// {    
//     for (int ii=0; ii<inum; ii++) {  // for each atom i in the simulation box     
//         int i = 2*ilist[ii];       // atom i        
//         int start = anumsum[ii];   
//         int m = anumsum[ii+1]-start;        // number of neighbors around i             
//         T f1 = 0.0;
//         T f2 = 0.0;
//         for (int l=0; l<m ; l++) {   // loop over each atom pair jk around atom i
//             int s = start + l;
//             int n = 2*s;                     
//             f1 +=  -(fij[n+0] + fik[n+0]);
//             f2 +=  -(fij[n+1] + fik[n+1]);
// //             // use atomicAdd on gpu
// //             int j = 2*aj[s];
// //             int k = 2*ak[s];            
// //             f[j+0] +=  fij[n+0];
// //             f[j+1] +=  fij[n+1];
// //             f[k+0] +=  fik[n+0];
// //             f[k+1] +=  fik[n+1];            
//         }
//         f[i+0] = f1;
//         f[i+1] = f2;                        
//     }
// }
// template void cpuAtomDecompositionTriplet2D(double*, double*, double*, int*, int*, int);
// template void cpuAtomDecompositionTriplet2D(float*, float*, float*, int*, int*, int);
// 
// template <typename T> void cpuForceDecompositionQuadruplet2D(T *f, T *fij, T *fik,  T *fil, int *ai, int *aj, int *ak, int *al, int ijklnum)
// {    
//     for (int ii=0; ii<ijklnum; ii++) {  // for each atom pair ij in the simulation box     
//         int i = 2*ai[ii];       // atom i
//         int j = 2*aj[ii];       // atom j
//         int k = 2*ak[ii];       // atom k
//         int l = 2*al[ii];       // atom l
//         int n = 2*ii;
//         // use atomicAdd on gpu
//         f[i+0] +=  -(fij[n+0]+fik[n+0]+fil[n+0]);
//         f[i+1] +=  -(fij[n+1]+fik[n+1]+fil[n+1]);
//         f[j+0] +=  fij[n+0];
//         f[j+1] +=  fij[n+1];
//         f[k+0] +=  fik[n+0];
//         f[k+1] +=  fik[n+1];
//         f[l+0] +=  fil[n+0];
//         f[l+1] +=  fil[n+1];
//     }
// }
// template void cpuForceDecompositionQuadruplet2D(double*, double*, double*, double*, int*, int*,  int*, int*, int);
// template void cpuForceDecompositionQuadruplet2D(float*, float*, float*, float*, int*, int*, int*, int*, int);
// 
// template <typename T> void cpuAtomDecompositionQuadruplet2D(T *f, T *fij, T *fik, T *fil, int *ilist, int *anumsum, int inum)
// {    
//     for (int ii=0; ii<inum; ii++) {  // for each atom i in the simulation box     
//         int i = 2*ilist[ii];       // atom i        
//         int start = anumsum[ii];   
//         int m = anumsum[ii+1]-start;        // number of neighbors around i             
//         T f1 = 0.0;
//         T f2 = 0.0;
//         for (int l=0; l<m ; l++) {   // loop over each atom pair jk around atom i
//             int s = start + l;
//             int n = 2*s;                     
//             f1 +=  -(fij[n+0] + fik[n+0] + fil[n+0]);
//             f2 +=  -(fij[n+1] + fik[n+1] + fil[n+1]);
//         }
//         f[i+0] = f1;
//         f[i+1] = f2;                        
//     }
// }
// template void cpuAtomDecompositionQuadruplet2D(double*, double*, double*, double*, int*, int*, int);
// template void cpuAtomDecompositionQuadruplet2D(float*, float*, float*, float*, int*, int*, int);
// 
// template <typename T> void cpuHalfForceDecomposition3D(T *f, T *fij, int *ai, int *aj, int ijnum)
// {    
//     for (int ii=0; ii<ijnum; ii++) {  // for each atom pair ij in the simulation box     
//         int i = 3*ai[ii];       // atom i
//         int j = 3*aj[ii];       // atom j
//         int l = 3*ii;
//         // use atomicAdd on gpu
//         f[i+0] +=  -fij[l+0];
//         f[i+1] +=  -fij[l+1];
//         f[i+2] +=  -fij[l+2];                 
//         f[j+0] +=  fij[l+0];
//         f[j+1] +=  fij[l+1];
//         f[j+2] +=  fij[l+2];                 
//     }
// }
// template void cpuHalfForceDecomposition3D(double*, double*, int*, int*, int);
// template void cpuHalfForceDecomposition3D(float*, float*, int*, int*, int);

// template <typename T> void cpuFullForceDecomposition3D(T *f, T *fij, int *ai, int *aj, int ijnum)
// {    
//     for (int ii=0; ii<ijnum; ii++) {  // for each atom pair ij in the simulation box     
//         int i = 3*ai[ii];             // atom i
//         int l = 3*ii;
//         // use atomicAdd on gpu
//         f[i+0] +=  -fij[l+0];
//         f[i+1] +=  -fij[l+1];
//         f[i+2] +=  -fij[l+2];                         
//     }
// }
// template void cpuFullForceDecomposition3D(double*, double*, int*, int*, int);
// template void cpuFullForceDecomposition3D(float*, float*, int*, int*, int);
// 
// template <typename T> void cpuIAtomDecomposition3D(T *f, T *fij, int *ilist, int *anumsum, int inum)
// {    
//     for (int ii=0; ii<inum; ii++) {  // for each atom i in the simulation box     
//         int i = 3*ilist[ii];       // atom i        
//         T f1 = 0.0;
//         T f2 = 0.0;        
//         T f3 = 0.0;        
//         int start = anumsum[ii];   
//         int m = anumsum[ii+1]-start;        // number of neighbors around i             
//         for (int l=0; l<m ; l++) {   // loop over each atom around atom i
//             int k = 3*(start + l);                     
//             f1 +=  -fij[k+0];
//             f2 +=  -fij[k+1];
//             f3 +=  -fij[k+2];             
//         }
//         f[i+0] = f1;
//         f[i+1] = f2;
//         f[i+2] = f3;
//     }
// }
// template void cpuIAtomDecomposition3D(double*, double*, int*, int*, int);
// template void cpuIAtomDecomposition3D(float*, float*, int*, int*, int);
// 
// template <typename T> void cpuJAtomDecomposition3D(T *f, T *fij, int *jlist, int *bnumsum, int *index, int jnum)
// {        
//     for (int ii=0; ii<jnum; ii++) {  // for each atom j in the simulation box     
//         int j = 3*jlist[ii];       // atom j
//         T f1 = 0.0;
//         T f2 = 0.0;                
//         T f3 = 0.0;                
//         // need to determine bnumsum and index 
//         int start = bnumsum[ii];   
//         int m = bnumsum[ii+1]-start;   // number of neighbors i around j  (j < i)           
//         for (int l=0; l<m ; l++) { // loop over each atom around i atom j -> pair ij 
//             int k = 3*index[start + l];                     
//             f1 +=  fij[k+0];
//             f2 +=  fij[k+1];
//             f3 +=  fij[k+2];
//         }
//         f[j+0] = f1;
//         f[j+1] = f2;
//         f[j+2] = f3;
//     }
// }
// template void cpuJAtomDecomposition3D(double*, double*, int*, int*, int*, int);
// template void cpuJAtomDecomposition3D(float*, float*, int*, int*, int*, int);

// template <typename T> void cpuForceDecompositionTriplet3D(T *f, T *fij, T *fik, int *ai, int *aj, int *ak, int ijknum)
// {    
//     for (int ii=0; ii<ijknum; ii++) {  // for each atom pair ij in the simulation box     
//         int i = 3*ai[ii];       // atom i
//         int j = 3*aj[ii];       // atom j
//         int k = 3*ak[ii];       // atom k
//         int l = 3*ii;
//         // use atomicAdd on gpu
//         f[i+0] +=  -(fij[l+0]+fik[l+0]);
//         f[i+1] +=  -(fij[l+1]+fik[l+1]);
//         f[i+2] +=  -(fij[l+2]+fik[l+2]);
//         f[j+0] +=  fij[l+0];
//         f[j+1] +=  fij[l+1];
//         f[j+2] +=  fij[l+2];
//         f[k+0] +=  fik[l+0];
//         f[k+1] +=  fik[l+1];
//         f[k+2] +=  fik[l+2];
//     }
// }
// template void cpuForceDecompositionTriplet3D(double*, double*, double*, int*, int*, int*, int);
// template void cpuForceDecompositionTriplet3D(float*, float*, float*, int*, int*, int*, int);
// 
// template <typename T> void cpuAtomDecompositionTriplet3D(T *f, T *fij, T *fik, int *ilist, int *anumsum, int inum)
// {    
//     for (int ii=0; ii<inum; ii++) {  // for each atom i in the simulation box     
//         int i = 3*ilist[ii];       // atom i        
//         int start = anumsum[ii];   
//         int m = anumsum[ii+1]-start;        // number of neighbors around i             
//         T f1 = 0.0;
//         T f2 = 0.0;
//         T f3 = 0.0;
//         for (int l=0; l<m ; l++) {   // loop over each atom pair jk around atom i
//             int s = start + l;
//             int n = 3*s;                     
//             f1 +=  -(fij[n+0] + fik[n+0]);
//             f2 +=  -(fij[n+1] + fik[n+1]);
//             f3 +=  -(fij[n+2] + fik[n+2]);
//         }
//         f[i+0] = f1;
//         f[i+1] = f2;                        
//         f[i+2] = f3;                        
//     }
// }
// template void cpuAtomDecompositionTriplet3D(double*, double*, double*, int*, int*, int);
// template void cpuAtomDecompositionTriplet3D(float*, float*, float*, int*, int*, int);
// 
// template <typename T> void cpuForceDecompositionQuadruplet3D(T *f, T *fij, T *fik,  T *fil, int *ai, int *aj, int *ak, int *al, int ijklnum)
// {    
//     for (int ii=0; ii<ijklnum; ii++) {  // for each atom pair ij in the simulation box     
//         int i = 3*ai[ii];       // atom i
//         int j = 3*aj[ii];       // atom j
//         int k = 3*ak[ii];       // atom k
//         int l = 3*al[ii];       // atom l
//         int n = 3*ii;
//         // use atomicAdd on gpu
//         f[i+0] +=  -(fij[n+0]+fik[n+0]+fil[n+0]);
//         f[i+1] +=  -(fij[n+1]+fik[n+1]+fil[n+1]);
//         f[i+2] +=  -(fij[n+2]+fik[n+2]+fil[n+2]);
//         f[j+0] +=  fij[n+0];
//         f[j+1] +=  fij[n+1];
//         f[j+2] +=  fij[n+2];
//         f[k+0] +=  fik[n+0];
//         f[k+1] +=  fik[n+1];
//         f[k+2] +=  fik[n+2];
//         f[l+0] +=  fil[n+0];
//         f[l+1] +=  fil[n+1];
//         f[l+2] +=  fil[n+2];
//     }
// }
// template void cpuForceDecompositionQuadruplet3D(double*, double*, double*, double*, int*, int*,  int*, int*, int);
// template void cpuForceDecompositionQuadruplet3D(float*, float*, float*, float*, int*, int*, int*, int*, int);
// 
// template <typename T> void cpuAtomDecompositionQuadruplet3D(T *f, T *fij, T *fik, T *fil, int *ilist, int *anumsum, int inum)
// {    
//     for (int ii=0; ii<inum; ii++) {  // for each atom i in the simulation box     
//         int i = 3*ilist[ii];       // atom i        
//         int start = anumsum[ii];   
//         int m = anumsum[ii+1]-start;        // number of neighbors around i             
//         T f1 = 0.0;
//         T f2 = 0.0;
//         T f3 = 0.0;
//         for (int l=0; l<m ; l++) {   // loop over each atom pair jk around atom i
//             int s = start + l;
//             int n = 3*s;                     
//             f1 +=  -(fij[n+0] + fik[n+0] + fil[n+0]);
//             f2 +=  -(fij[n+1] + fik[n+1] + fil[n+1]);
//             f3 +=  -(fij[n+2] + fik[n+2] + fil[n+2]);
//         }
//         f[i+0] = f1;
//         f[i+1] = f2;                        
//         f[i+2] = f3;                        
//     }
// }
// template void cpuAtomDecompositionQuadruplet3D(double*, double*, double*, double*, int*, int*, int);
// template void cpuAtomDecompositionQuadruplet3D(float*, float*, float*, float*, int*, int*, int);
// 
// //********************************************************************************************//
// 
// template <typename T> void cpuGetNeighPairs(T *xij, T *qi, T *qj, T *x, T *q, int *ai, int *aj,  
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
// template void cpuGetNeighPairs(double*, double*, double*, double*, double*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int, int, int, int);
// template void cpuGetNeighPairs(float*, float*, float*, float*, float*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int, int, int, int);
// 
// template <typename T> void cpuGetHalfNeighPairs(T *xij, T *qi, T *qj, T *x, T *q, int *ai, int *aj,  
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
// template void cpuGetHalfNeighPairs(double*, double*, double*, double*, double*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int, int, int, int);
// template void cpuGetHalfNeighPairs(float*, float*, float*, float*, float*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int, int, int, int);
// 
// template <typename T> void cpuGetNeighPairs(T *xij, T *qi, T *qj, T *x, T *q, int *ai, int *aj,  
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
// //     cpuCumsum(anumsum, anum, inum+1);                 
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
// template void cpuGetNeighPairs(double*, double*, double*, double*, double*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int, int, int, int, int);
// template void cpuGetNeighPairs(float*, float*, float*, float*, float*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int, int, int, int, int);
// 
// void cpuFindAtom(int *tlist, int* ilist, int *atomtype, int typei, int inum)
// {
//     for (int ii=0; ii<inum; ii++) {  // for each atom i in the simulation box     
//         int i = ilist[ii];           // atom i
//         tlist[ii] = (atomtype[i] == typei) ? 1 : 0;        
//     }    
// }
// 
// void cpuFindStart(int *start, int* slist, int inum)
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
// void cpuCopyAtIndex(int* tlist, int* ilist, int* index, int *start)
// {
//     for (int ii=0; ii<start[1]; ii++) {                  
//         tlist[ii] = ilist[index[start[0]+ii]];
//     }    
// }
// 
// void cpuGetIlist(int *tlist, int *start, int *ilist, int *slist, int *index, int *atomtype, int typei, int inum)
// {
//     cpuFindAtom(tlist, ilist, atomtype, typei, inum);
//  
//     // sort list    
//     cpuMergeSort(slist, index, tlist, inum);
//         
//     cpuFindStart(start, slist, inum);        
//         
//     cpuCopyAtIndex(tlist, ilist, index, start);
// }
// 
// template <typename T> void cpuGetNeighPairs(T *xij, T *qi, T *qj, T *x, T *q, int *ai, int *aj, 
//       int *ti, int *tj, int *anum, int *anumsum, int *ilist, int *tlist, int *alist, int *neighlist, int *neighnum, 
//       int *atomtype, int *slist, int *start, int *index, int typei, int inum, int jnum, int ncq, int dim)
// {    
//     // form tlist from ilist
//     cpuGetIlist(tlist, start, ilist, slist, index, atomtype, typei, inum);
//             
//     cpuGetNeighPairs(xij, qi, qj, x, q, ai, aj, ti, tj, anum, anumsum, tlist, alist,
//             neighlist, neighnum, atomtype, start[1], jnum, ncq, dim);             
// }
// template void cpuGetNeighPairs(double*, double*, double*, double*, double*, int*, int*, int*, int*, int*, 
//         int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int, int, int, int, int);
// template void cpuGetNeighPairs(float*, float*, float*, float*, float*, int*, int*, int*, int*, int*, 
//         int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int, int, int, int, int);
// 
// template <typename T> void cpuGetNeighPairs(T *xij, T *qi, T *qj, T *x, T *q, int *ai, int *aj, 
//       int *ti, int *tj, int *anum, int *anumsum, int *ilist, int *tlist, int *alist, int *neighlist, int *neighnum, 
//       int *atomtype, int *slist, int *start, int *index, int typei, int typej, int inum, int jnum, int ncq, int dim)
// {    
//     // form tlist from ilist
//     cpuGetIlist(tlist, start, ilist, slist, index, atomtype, typei, inum);
//             
//     cpuGetNeighPairs(xij, qi, qj, x, q, ai, aj, ti, tj, tlist, alist, anum, anumsum,
//             neighlist, neighnum, atomtype, typej, start[1], jnum, ncq, dim);             
// }
// template void cpuGetNeighPairs(double*, double*, double*, double*, double*, int*, int*, int*, int*, int*, 
//         int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int, int, int, int, int, int);
// template void cpuGetNeighPairs(float*, float*, float*, float*, float*, int*, int*, int*, int*, int*, 
//         int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int, int, int, int, int, int);
// 
// template <typename T> void cpuGetNeighTriplets(T *xij, T *xik, T *qi, T *qj, T *qk, T *x, T *q, int *ai, int *aj,  
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
//     cpuCumsum(anumsum, anum, inum+1);             
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
// template void cpuGetNeighTriplets(double*, double*, double*, double*, double*, double*, double*, 
//         int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int, int, int, int);
// template void cpuGetNeighTriplets(float*, float*, float*, float*, float*, float*, float*, 
//         int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int, int, int, int);
// 
// template <typename T> void cpuGetNeighFullTriplets(T *xij, T *xik, T *qi, T *qj, T *qk, T *x, T *q, int *ai, int *aj,  
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
//     cpuCumsum(anumsum, anum, inum+1);             
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
// template void cpuGetNeighFullTriplets(double*, double*, double*, double*, double*, double*, double*, 
//         int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int, int, int, int);
// template void cpuGetNeighFullTriplets(float*, float*, float*, float*, float*, float*, float*, 
//         int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int, int, int, int);
// 
// template <typename T> void cpuGetNeighTriplets(T *xij, T *xik, T *qi, T *qj, T *qk, T *x, T *q, int *ai, int *aj,  
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
//     cpuCumsum(anumsum, anum, inum+1);             
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
// template void cpuGetNeighTriplets(double*, double*, double*, double*, double*, double*, double*, 
//         int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int, int, int, int, int, int);
// template void cpuGetNeighTriplets(float*, float*, float*, float*, float*, float*, float*, 
//         int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int, int, int, int, int, int);
// 
// template <typename T> void cpuGetNeighTriplets(T *xij, T *xik, T *qi, T *qj, T *qk, T *x, T *q, int *ai, int *aj,  
//       int *ak, int *ti, int *tj, int *tk, int *anum, int *anumsum, int *ilist,  int *tlist, int *alist, int *neighlist, int *neighnum, 
//       int *atomtype, int *slist, int *start, int *index, int typei, int typej, int typek, int inum, int jnum, int ncq, int dim)
// {    
//     // form tlist from ilist
//     cpuGetIlist(tlist, start, ilist, slist, index, atomtype, typei, inum);
//             
//     cpuGetNeighTriplets(xij, xik, qi, qj,qk,  x, q, ai, aj, ak, ti, tj, tk, anum, anumsum, tlist, alist,
//             neighlist, neighnum, atomtype, typej, typek, start[1], jnum, ncq, dim);             
// }
// template void cpuGetNeighTriplets(double*, double*, double*, double*, double*, double*, double*, 
//         int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, 
//         int, int, int, int, int, int, int);
// template void cpuGetNeighTriplets(float*, float*, float*, float*, float*, float*, float*, 
//         int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, 
//         int, int, int, int, int, int, int);
// 
// template <typename T> void cpuGetNeighFullTriplets(T *xij, T *xik, T *qi, T *qj, T *qk, T *x, T *q, int *ai, int *aj,  
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
//     cpuCumsum(anumsum, anum, inum+1);             
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
// template void cpuGetNeighFullTriplets(double*, double*, double*, double*, double*, double*, double*, 
//         int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int, int, int, int, int, int);
// template void cpuGetNeighFullTriplets(float*, float*, float*, float*, float*, float*, float*, 
//         int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int, int, int, int, int, int);
// 
// template <typename T> void cpuGetNeighFullTriplets(T *xij, T *xik, T *qi, T *qj, T *qk, T *x, T *q, int *ai, int *aj,  
//       int *ak, int *ti, int *tj, int *tk, int *anum, int *anumsum, int *ilist,  int *tlist, int *alist, int *neighlist, int *neighnum, 
//       int *atomtype, int *slist, int *start, int *index, int typei, int typej, int typek, int inum, int jnum, int ncq, int dim)
// {    
//     // form tlist from ilist
//     cpuGetIlist(tlist, start, ilist, slist, index, atomtype, typei, inum);
//             
//     cpuGetNeighFullTriplets(xij, xik, qi, qj,qk,  x, q, ai, aj, ak, ti, tj, tk, anum, anumsum, tlist, alist,
//             neighlist, neighnum, atomtype, typej, typek, start[1], jnum, ncq, dim);             
// }
// template void cpuGetNeighFullTriplets(double*, double*, double*, double*, double*, double*, double*, 
//         int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, 
//         int, int, int, int, int, int, int);
// template void cpuGetNeighFullTriplets(float*, float*, float*, float*, float*, float*, float*, 
//         int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, 
//         int, int, int, int, int, int, int);
// 
// 
// template <typename T> void cpuGetNeighFullTriplets(T *xij, T *xik, T *qi, T *qj, T *qk, T *x, T *q, int *ai, int *aj,  
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
//     cpuCumsum(anumsum, anum, inum+1);             
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
// template void cpuGetNeighFullTriplets(double*, double*, double*, double*, double*, double*, double*, 
//         int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int, int, int, int, int, int);
// template void cpuGetNeighFullTriplets(float*, float*, float*, float*, float*, float*, float*, 
//         int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int, int, int, int, int, int);
// 
// template <typename T> void cpuGetNeighPairs(T *xij, T *qi, T *qj, T *x, T *q, int *ai, int *aj,  
//       int *ti, int *tj, int *alist, int *neighlist, int *neighnumsum, int *atomtype, 
//       int istart, int iend, int jnum, int ncq, int dim)
// {    
//     for (int i=istart; i<iend; i++) {  // for each atom i in the simulation box     
//         int itype = atomtype[i];        
//         int m = neighnumsum[i+1] - neighnumsum[i];        // number of neighbors around i             
//         int start = neighnumsum[i] - neighnumsum[istart];   
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
// template void cpuGetNeighPairs(double*, double*, double*, double*, double*, int*, int*, int*, int*, int*, int*, int*, int*, int, int, int, int, int);
// template void cpuGetNeighPairs(float*, float*, float*, float*, float*, int*, int*, int*, int*, int*, int*, int*, int*, int, int, int, int, int);
// 
// template <typename T> void cpuGetNeighPairs(T *xij, T *qi, T *qj, T *x, T *q, int *ai, int *aj,  
//       int *ti, int *tj, int *anum, int *anumsum, int *alist, int *neighlist, int *neighnumsum, 
//       int *atomtype, int typej, int istart, int iend, int jnum, int ncq, int dim)
// {        
//     // form anum
//     for (int i=istart; i<iend; i++) {       
//         int m = neighnumsum[i+1] - neighnumsum[i];        // number of neighbors around i     
//         anum[i-istart] = 0;              
//         for (int l=0; l<m ; l++) {   // loop over each atom around atom i
//             int g = neighlist[l + jnum*i];  // atom j           
//             int j = alist[g];  // atom j
//             if (atomtype[j] == typej) 
//                 anum[i-istart] += 1;
//         }                
//     }
//     
//     // a list contains the starting positions of the first neighbor 
//     int inum = iend - istart; // number of atoms i
//     cpuCumsum(anumsum, anum, inum+1);             
//     
//     for (int i=istart; i<iend; i++) {  // for each atom i in the simulation box     
//         int itype = atomtype[i];
//         int m = anum[i-istart];        // number of neighbors around i             
//         int start = anumsum[i-istart];   
//         int inc = 0;
//         for (int l=0; l<m ; l++) {   // loop over each atom around atom i
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
// template void cpuGetNeighPairs(double*, double*, double*, double*, double*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int, int, int, int, int, int);
// template void cpuGetNeighPairs(float*, float*, float*, float*, float*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int, int, int, int, int, int);
// 
// template <typename T> void cpuGetNeighTriplets(T *xij, T *xik, T *qi, T *qj, T *qk, T *x, T *q, int *ai, int *aj,  
//       int *ak, int *ti, int *tj, int *tk, int *anum, int *anumsum, int *alist, int *neighlist, int *neighnumsum, 
//       int *atomtype, int istart, int iend, int jnum, int ncq, int dim)
// {        
//     // form anum
//     for (int i=istart; i<iend; i++) {  // for each atom i in the simulation box     
//         int m = neighnumsum[i+1] - neighnumsum[i];        // number of neighbors around i             
//         anum[i] = 0;              
//         for (int lj=0; lj<m ; lj++) {   // loop over each atom around atom i
//             int gj = neighlist[lj + jnum*i];  // atom j           
//             int j = alist[gj];  // atom j
//             for (int lk=0; lk<m; lk++) {
//                 int gk = neighlist[lk + jnum*i];  // atom k
//                 int k = alist[gk];  // atom k
//                 if (j < k)
//                     anum[i] += 1;
//             }            
//         }                
//     }
//     
//     // a list contains the starting positions of the first neighbor 
//     int inum = iend - istart;
//     cpuCumsum(anumsum, anum, inum+1);             
//     
//     for (int i=istart; i<iend; i++) {  // for each atom i in the simulation box           
//         int itype = atomtype[i];
//         int m = anum[i-istart];        // number of neighbors around i             
//         int start = anumsum[i-istart];   
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
// template void cpuGetNeighTriplets(double*, double*, double*, double*, double*, double*, double*, 
//         int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int, int, int, int, int);
// template void cpuGetNeighTriplets(float*, float*, float*, float*, float*, float*, float*, 
//         int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int, int, int, int, int);

//********************************************************************************************//

// template <typename T> void cpuFullAtomDecomposition2D(T *fi, T *fij, int *neighnumsum, int istart, int iend)
// {        
//     for (int i=istart; i<iend; i++) {  // for each atom i in the simulation box     
//         int m = neighnumsum[i+1] - neighnumsum[i];        // number of neighbors around i             
//         int start = neighnumsum[i] - neighnumsum[istart];       
//         T f1 = 0.0;
//         T f2 = 0.0;
//         for (int l=0; l<m ; l++) {   // loop over each atom j around atom i
//             int k = start + l;                     
//             f1 +=  -fij[2*k+0];
//             f2 +=  -fij[2*k+1];
//         }
//         fi[2*i+0] = f1;
//         fi[2*i+1] = f2;
//     }
// }
// template void cpuFullAtomDecomposition2D(double*, double*, int*, int, int);
// template void cpuFullAtomDecomposition2D(float*, float*, int*, int, int);
// 
// template <typename T> void cpuHalfAtomDecomposition2D(T *fi, T *fij, int *neighnumsum, 
//         int *bnumsum, int *index, int istart, int iend)
// {    
//     for (int i=istart; i<iend; i++) {  // for each atom i in the simulation box     
//         int m = neighnumsum[i+1] - neighnumsum[i];        // number of neighbors around i             
//         int start = neighnumsum[i] - neighnumsum[istart];               
//         T f1 = 0.0;
//         T f2 = 0.0;        
//         for (int l=0; l<m ; l++) {   // loop over each atom around j atom i -> pair ij 
//             int k = start + l;                     
//             f1 +=  -fij[2*k+0];
//             f2 +=  -fij[2*k+1];
//         }                
//         
//         int ii = i - istart;
//         start = bnumsum[ii];   
//         m = bnumsum[ii+1]-start;   // number of neighbors j around i  (j < i)           
//         for (int l=0; l<m ; l++) { // loop over each atom around j atom i -> pair ji 
//             int k = index[start + l];                     
//             f1 +=  fij[2*k+0];
//             f2 +=  fij[2*k+1];
//         }
//         fi[2*i+0] = f1;
//         fi[2*i+1] = f2;
//     }
// }
// template void cpuHalfAtomDecomposition2D(double*, double*, int*, int*, int*, int, int);
// template void cpuHalfAtomDecomposition2D(float*, float*, int*, int*, int*, int, int);

// template <typename T> void cpuHalfAtomDecomposition2D(T *f, T *fij, int *ilist, int *anumsum, 
//         int *jlist, int *bnumsum, int *index, int inum, int jnum)
// {    
//     cpuIAtomDecomposition2D(f, fij, ilist, anumsum, inum);
//     cpuJAtomDecomposition2D(f, fij, jlist, bnumsum, index, jnum);   
// }
// template void cpuHalfAtomDecomposition2D(double*, double*, int*, int*, int*, int*, int*, int, int);
// template void cpuHalfAtomDecomposition2D(float*, float*, int*, int*, int*, int*, int*, int, int);

// 
// template <typename T> void cpuFullAtomDecomposition3D(T *fi, T *fij, int *neighnumsum, int istart, int iend)
// {        
//     for (int i=istart; i<iend; i++) {  // for each atom i in the simulation box     
//         int m = neighnumsum[i+1] - neighnumsum[i];        // number of neighbors around i             
//         int start = neighnumsum[i] - neighnumsum[istart];       
//         T f1 = 0.0;
//         T f2 = 0.0;
//         T f3 = 0.0;        
//         for (int l=0; l<m ; l++) {   // loop over each atom around atom i
//             int k = start + l;                     
//             f1 +=  -fij[3*k+0];
//             f2 +=  -fij[3*k+1];
//             f3 +=  -fij[3*k+2];             
//         }
//         fi[3*i+0] = f1;
//         fi[3*i+1] = f2;
//         fi[3*i+2] = f3;        
//     }
// }
// template void cpuFullAtomDecomposition3D(double*, double*, int*, int, int);
// template void cpuFullAtomDecomposition3D(float*, float*, int*, int, int);
// 
// template <typename T> void cpuHalfAtomDecomposition3D(T *fi, T *fij, int *neighnumsum, 
//         int *bnumsum, int *index, int istart, int iend)
// {    
//     for (int i=istart; i<iend; i++) {  // for each atom i in the simulation box     
//         int m = neighnumsum[i+1] - neighnumsum[i];        // number of neighbors around i             
//         int start = neighnumsum[i] - neighnumsum[istart];               
//         T f1 = 0.0;
//         T f2 = 0.0;        
//         T f3 = 0.0;        
//         for (int l=0; l<m ; l++) {   // loop over each atom around j atom i -> pair ij 
//             int k = start + l;             
//             f1 +=  -fij[3*k+0];
//             f2 +=  -fij[3*k+1];
//             f3 +=  -fij[3*k+2];                         
//         }                
//         
//         // need to determine bnumsum and index 
//         int ii = i - istart;
//         start = bnumsum[ii];   
//         m = bnumsum[ii+1]-start;   // number of neighbors j around i  (j < i)           
//         for (int l=0; l<m ; l++) { // loop over each atom around j atom i -> pair ji 
//             int k = index[start + l];     
//             f1 +=  fij[3*k+0];
//             f2 +=  fij[3*k+1];
//             f3 +=  fij[3*k+2];             
//         }
//         fi[3*i+0] = f1;
//         fi[3*i+1] = f2;
//         fi[3*i+2] = f3;        
//     }
// }
// template void cpuHalfAtomDecomposition3D(double*, double*, int*, int*, int*, int, int);
// template void cpuHalfAtomDecomposition3D(float*, float*, int*, int*, int*, int, int);

// template <typename T> void cpuHalfAtomDecomposition3D(T *fi, T *fij, int *ilist, int *anumsum, 
//         int *bnumsum, int *index, int inum)
// {    
//     for (int ii=0; ii<inum; ii++) {  // for each atom i in the simulation box     
//         int i = ilist[ii];       // atom i
//         T f1 = 0.0;
//         T f2 = 0.0;        
//         T f3 = 0.0;        
//         int start = anumsum[ii];   
//         int m = anumsum[ii+1]-start; // number of neighbors j around i  (j > i)             
//         for (int l=0; l<m ; l++) {   // loop over each atom around j atom i -> pair ij 
//             int k = start + l;                     
//             f1 +=  -fij[3*k+0];
//             f2 +=  -fij[3*k+1];
//             f3 +=  -fij[3*k+2];             
//         }                
//         // need to determine bnumsum and index
//         start = bnumsum[ii];   
//         m = bnumsum[ii+1]-start;   // number of neighbors j around i  (j < i)           
//         for (int l=0; l<m ; l++) { // loop over each atom around j atom i -> pair ji 
//             int k = index[start + l];                     
//             f1 +=  fij[3*k+0];
//             f2 +=  fij[3*k+1];
//             f3 +=  fij[3*k+2];             
//         }
//         fi[3*i+0] = f1;
//         fi[3*i+1] = f2;
//         fi[3*i+2] = f3;
//     }
// }
// template void cpuHalfAtomDecomposition3D(double*, double*, int*, int*, int*, int*, int);
// template void cpuHalfAtomDecomposition3D(float*, float*, int*, int*, int*, int*, int);


#endif


