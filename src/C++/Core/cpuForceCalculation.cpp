// assign threads to same-colored cells
// each thread gets its cell number
// each thread loops over atoms in its cell
//     compute xij, ti, tj
//     compute forces
//     compute basis function derivatives

// assign threads to each atom
// each thread loops over atoms 
//     compute basis function derivatives
//     assemble forces using atomicAdd to avoid 

// analyze aj

#ifndef __CPUFORCECALCULATION
#define __CPUFORCECALCULATION


template <typename T> void cpuGetNeighPairs(T *xij, T *qi, T *qj, T *x, T *q, int *ai, int *aj,  
      int *ti, int *tj, int *alist, int *neighlist, int *neighnumsum, int *atomtype, 
      int istart, int iend, int jnum, int ncq, int dim)
{    
    for (int i=istart; i<iend; i++) {  // for each atom i in the simulation box     
        int itype = atomtype[i];        
        int m = neighnumsum[i+1] - neighnumsum[i];        // number of neighbors around i             
        int start = neighnumsum[i] - neighnumsum[istart];   
        for (int l=0; l<m ; l++) {   // loop over each atom around atom i
            int g = neighlist[l +jnum*i];  // ghost index of atom j  
            int j = alist[g];  // atom j
            int k = start + l;                                     
            ai[k]        = i;
            aj[k]        = j; // should be alist[j];         
            ti[k]        = itype;       
            tj[k]        = atomtype[j]; // should be atomtype[alist[j]];         
            for (int p=0; p<dim; p++) 
                xij[k*dim+p]   = x[g*dim+p] - x[i*dim+p];  // xj - xi            
            for (int p=0; p<ncq; p++) {                
                qi[k*ncq+p] = q[i*ncq+p];
                qj[k*ncq+p] = q[j*ncq+p];
            }                
        }
    }    
}
template void cpuGetNeighPairs(double*, double*, double*, double*, double*, int*, int*, int*, int*, int*, int*, int*, int*, int, int, int, int, int);
template void cpuGetNeighPairs(float*, float*, float*, float*, float*, int*, int*, int*, int*, int*, int*, int*, int*, int, int, int, int, int);

template <typename T> void cpuGetNeighPairs(T *xij, T *qi, T *qj, T *x, T *q, int *ai, int *aj,  
      int *ti, int *tj, int *anum, int *anumsum, int *alist, int *neighlist, int *neighnumsum, 
      int *atomtype, int typej, int istart, int iend, int jnum, int ncq, int dim)
{        
    // form anum
    for (int i=istart; i<iend; i++) {       
        int m = neighnumsum[i+1] - neighnumsum[i];        // number of neighbors around i     
        anum[i-istart] = 0;              
        for (int l=0; l<m ; l++) {   // loop over each atom around atom i
            int g = neighlist[l + jnum*i];  // atom j           
            int j = alist[g];  // atom j
            if (atomtype[j] == typej) 
                anum[i-istart] += 1;
        }                
    }
    
    // a list contains the starting positions of the first neighbor 
    int inum = iend - istart; // number of atoms i
    cpuCumsum(anumsum, anum, inum+1);             
    
    for (int i=istart; i<iend; i++) {  // for each atom i in the simulation box     
        int itype = atomtype[i];
        int m = anum[i-istart];        // number of neighbors around i             
        int start = anumsum[i-istart];   
        int inc = 0;
        for (int l=0; l<m ; l++) {   // loop over each atom around atom i
            int g = neighlist[l +jnum*i];  // ghost index of atom j  
            int j = alist[g];  // atom j
            if (atomtype[j] == typej)  {
                int k = start + inc;                         
                ai[k]        = i;
                aj[k]        = j; // should be alist[j];         
                ti[k]        = itype;       
                tj[k]        = atomtype[j]; // should be atomtype[alist[j]];                                        
                for (int p=0; p<dim; p++) 
                    xij[k*dim+p]   = x[g*dim+p] - x[i*dim+p];  // xj - xi                            
                for (int p=0; p<ncq; p++) {                
                    qi[k*ncq+p] = q[i*ncq+p];
                    qj[k*ncq+p] = q[j*ncq+p];
                }                
                inc += 1;
            }
        }
    }    
}
template void cpuGetNeighPairs(double*, double*, double*, double*, double*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int, int, int, int, int, int);
template void cpuGetNeighPairs(float*, float*, float*, float*, float*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int, int, int, int, int, int);

template <typename T> void cpuGetNeighTriplets(T *xij, T *xik, T *qi, T *qj, T *qk, T *x, T *q, int *ai, int *aj,  
      int *ak, int *ti, int *tj, int *tk, int *anum, int *anumsum, int *alist, int *neighlist, int *neighnumsum, 
      int *atomtype, int istart, int iend, int jnum, int ncq, int dim)
{        
    // form anum
    for (int i=istart; i<iend; i++) {  // for each atom i in the simulation box     
        int m = neighnumsum[i+1] - neighnumsum[i];        // number of neighbors around i             
        anum[i] = 0;              
        for (int lj=0; lj<m ; lj++) {   // loop over each atom around atom i
            int gj = neighlist[lj + jnum*i];  // atom j           
            int j = alist[gj];  // atom j
            for (int lk=0; lk<m; lk++) {
                int gk = neighlist[lk + jnum*i];  // atom k
                int k = alist[gk];  // atom k
                if (j < k)
                    anum[i] += 1;
            }            
        }                
    }
    
    // a list contains the starting positions of the first neighbor 
    int inum = iend - istart;
    cpuCumsum(anumsum, anum, inum+1);             
    
    for (int i=istart; i<iend; i++) {  // for each atom i in the simulation box           
        int itype = atomtype[i];
        int m = anum[i-istart];        // number of neighbors around i             
        int start = anumsum[i-istart];   
        int inc = 0;
        for (int lj=0; lj<m ; lj++) {   // loop over each atom around atom i
            int gj = neighlist[lj +jnum*i];  // ghost index of atom j  
            int j = alist[gj];  // atom j
            for (int lk=0; lk<m; lk++) {
                int gk = neighlist[lk + jnum*i];  // atom k
                int k = alist[gk];  // atom k
                if (j < k) {
                    int n = start + inc;                         
                    for (int p=0; p<dim; p++) {
                        xij[n*dim+p]   = x[gj*dim+p] - x[i*dim+p];  // xj - xi            
                        xik[n*dim+p]   = x[gk*dim+p] - x[i*dim+p];  // xj - xi            
                    }
                    ai[n]        = i;
                    aj[n]        = j; // should be alist[j];         
                    ak[n]        = k; // should be alist[j];         
                    ti[n]        = itype;       
                    tj[n]        = atomtype[j]; // should be atomtype[alist[j]];                                        
                    tk[n]        = atomtype[k]; // should be atomtype[alist[j]];                                        
                    for (int p=0; p<ncq; p++) {                
                        qi[n*ncq+p] = q[i*ncq+p];
                        qj[n*ncq+p] = q[j*ncq+p];
                        qk[n*ncq+p] = q[k*ncq+p];
                    }                
                    inc += 1;
                }
            }                                    
        }
    }    
}
template void cpuGetNeighTriplets(double*, double*, double*, double*, double*, double*, double*, 
        int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int, int, int, int, int);
template void cpuGetNeighTriplets(float*, float*, float*, float*, float*, float*, float*, 
        int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int, int, int, int, int);

//********************************************************************************************//

template <typename T> void cpuHalfForceDecomposition2D(T *fi, T *fij, int *ai, int *aj, int ijnum)
{    
    for (int ii=0; ii<ijnum; ii++) {  // for each atom pair ij in the simulation box     
        int i = ai[ii];       // atom i
        int j = aj[ii];       // atom j
        // use atomicAdd on gpu
        fi[2*i+0] +=  -fij[2*ii+0];
        fi[2*i+1] +=  -fij[2*ii+1];
        fi[2*j+0] +=  fij[2*ii+0];
        fi[2*j+1] +=  fij[2*ii+1];
    }
}
template void cpuHalfForceDecomposition2D(double*, double*, int*, int*, int);
template void cpuHalfForceDecomposition2D(float*, float*, int*, int*, int);

template <typename T> void cpuFullForceDecomposition2D(T *fi, T *fij, int *ai, int *aj, int ijnum)
{    
    for (int ii=0; ii<ijnum; ii++) {  // for each atom pair ij in the simulation box     
        int i = ai[ii];       // atom i
        // use atomicAdd on gpu
        fi[2*i+0] +=  -fij[2*ii+0];
        fi[2*i+1] +=  -fij[2*ii+1];
    }
}
template void cpuFullForceDecomposition2D(double*, double*, int*, int*, int);
template void cpuFullForceDecomposition2D(float*, float*, int*, int*, int);

template <typename T> void cpuFullAtomDecomposition2D(T *fi, T *fij, int *neighnumsum, int istart, int iend)
{        
    for (int i=istart; i<iend; i++) {  // for each atom i in the simulation box     
        int m = neighnumsum[i+1] - neighnumsum[i];        // number of neighbors around i             
        int start = neighnumsum[i] - neighnumsum[istart];       
        T f1 = 0.0;
        T f2 = 0.0;
        for (int l=0; l<m ; l++) {   // loop over each atom around atom i
            int k = start + l;                     
            f1 +=  -fij[2*k+0];
            f2 +=  -fij[2*k+1];
        }
        fi[2*i+0] = f1;
        fi[2*i+1] = f2;
    }
}
template void cpuFullAtomDecomposition2D(double*, double*, int*, int, int);
template void cpuFullAtomDecomposition2D(float*, float*, int*, int, int);

template <typename T> void cpuHalfAtomDecomposition2D(T *fi, T *fij, int *neighnumsum, 
        int *bnumsum, int *index, int istart, int iend)
{    
    for (int i=istart; i<iend; i++) {  // for each atom i in the simulation box     
        int m = neighnumsum[i+1] - neighnumsum[i];        // number of neighbors around i             
        int start = neighnumsum[i] - neighnumsum[istart];               
        T f1 = 0.0;
        T f2 = 0.0;        
        for (int l=0; l<m ; l++) {   // loop over each atom around j atom i -> pair ij 
            int k = start + l;                     
            f1 +=  -fij[2*k+0];
            f2 +=  -fij[2*k+1];
        }                
        
        int ii = i - istart;
        start = bnumsum[ii];   
        m = bnumsum[ii+1]-start;   // number of neighbors j around i  (j < i)           
        for (int l=0; l<m ; l++) { // loop over each atom around j atom i -> pair ji 
            int k = index[start + l];                     
            f1 +=  fij[2*k+0];
            f2 +=  fij[2*k+1];
        }
        fi[2*i+0] = f1;
        fi[2*i+1] = f2;
    }
}
template void cpuHalfAtomDecomposition2D(double*, double*, int*, int*, int*, int, int);
template void cpuHalfAtomDecomposition2D(float*, float*, int*, int*, int*, int, int);

template <typename T> void cpuFullAtomDecomposition2D(T *fi, T *fij, int *ilist, int *anumsum, int inum)
{    
    for (int ii=0; ii<inum; ii++) {  // for each atom i in the simulation box     
        int i = ilist[ii];       // atom i        
        int start = anumsum[ii];   
        int m = anumsum[ii+1]-start;        // number of neighbors around i             
        T f1 = 0.0;
        T f2 = 0.0;
        for (int l=0; l<m ; l++) {   // loop over each atom around atom i
            int k = start + l;                     
            f1 +=  -fij[2*k+0];
            f2 +=  -fij[2*k+1];
        }
        fi[2*i+0] = f1;
        fi[2*i+1] = f2;
    }
}
template void cpuFullAtomDecomposition2D(double*, double*, int*, int*, int);
template void cpuFullAtomDecomposition2D(float*, float*, int*, int*, int);

template <typename T> void cpuHalfAtomDecomposition2D(T *fi, T *fij, int *ilist, int *anumsum, 
        int *bnumsum, int *index, int inum)
{    
    for (int ii=0; ii<inum; ii++) {  // for each atom i in the simulation box     
        int i = ilist[ii];       // atom i
        T f1 = 0.0;
        T f2 = 0.0;        
        int start = anumsum[ii];   
        int m = anumsum[ii+1]-start; // number of neighbors j around i  (j > i)             
        for (int l=0; l<m ; l++) {   // loop over each atom around j atom i -> pair ij 
            int k = start + l;                     
            f1 +=  -fij[2*k+0];
            f2 +=  -fij[2*k+1];
        }                
        start = bnumsum[ii];   
        m = bnumsum[ii+1]-start;   // number of neighbors j around i  (j < i)           
        for (int l=0; l<m ; l++) { // loop over each atom around j atom i -> pair ji 
            int k = index[start + l];                     
            f1 +=  fij[2*k+0];
            f2 +=  fij[2*k+1];
        }
        fi[2*i+0] = f1;
        fi[2*i+1] = f2;
    }
}
template void cpuHalfAtomDecomposition2D(double*, double*, int*, int*, int*, int*, int);
template void cpuHalfAtomDecomposition2D(float*, float*, int*, int*, int*, int*, int);

// force calculation for half neighbor list
template <typename T> void cpuHalfForceDecomposition3D(T *fi, T *fij, int *ai, int *aj, int ijnum)
{    
    for (int ii=0; ii<ijnum; ii++) {  // for each atom pair ij in the simulation box     
        int i = ai[ii];       // atom i
        int j = aj[ii];       // atom j
        // use atomicAdd on gpu
        fi[3*i+0] +=  -fij[3*ii+0];
        fi[3*i+1] +=  -fij[3*ii+1];
        fi[3*i+2] +=  -fij[3*ii+2];                 
        fi[3*j+0] +=  fij[3*ii+0];
        fi[3*j+1] +=  fij[3*ii+1];
        fi[3*j+2] +=  fij[3*ii+2];                 
    }
}
template void cpuHalfForceDecomposition3D(double*, double*, int*, int*, int);
template void cpuHalfForceDecomposition3D(float*, float*, int*, int*, int);

// force calculation for full neighbor list
template <typename T> void cpuFullForceDecomposition3D(T *fi, T *fij, int *ai, int *aj, int ijnum)
{    
    for (int ii=0; ii<ijnum; ii++) {  // for each atom pair ij in the simulation box     
        int i = ai[ii];       // atom i
        // use atomicAdd on gpu
        fi[3*i+0] +=  -fij[3*ii+0];
        fi[3*i+1] +=  -fij[3*ii+1];
        fi[3*i+2] +=  -fij[3*ii+2];                         
    }
}
template void cpuFullForceDecomposition3D(double*, double*, int*, int*, int);
template void cpuFullForceDecomposition3D(float*, float*, int*, int*, int);

template <typename T> void cpuFullAtomDecomposition3D(T *fi, T *fij, int *neighnumsum, int istart, int iend)
{        
    for (int i=istart; i<iend; i++) {  // for each atom i in the simulation box     
        int m = neighnumsum[i+1] - neighnumsum[i];        // number of neighbors around i             
        int start = neighnumsum[i] - neighnumsum[istart];       
        T f1 = 0.0;
        T f2 = 0.0;
        T f3 = 0.0;        
        for (int l=0; l<m ; l++) {   // loop over each atom around atom i
            int k = start + l;                     
            f1 +=  -fij[3*k+0];
            f2 +=  -fij[3*k+1];
            f3 +=  -fij[3*k+2];             
        }
        fi[3*i+0] = f1;
        fi[3*i+1] = f2;
        fi[3*i+2] = f3;        
    }
}
template void cpuFullAtomDecomposition3D(double*, double*, int*, int, int);
template void cpuFullAtomDecomposition3D(float*, float*, int*, int, int);

// force calculation for half neighbor list
template <typename T> void cpuHalfAtomDecomposition3D(T *fi, T *fij, int *neighnumsum, 
        int *bnumsum, int *index, int istart, int iend)
{    
    for (int i=istart; i<iend; i++) {  // for each atom i in the simulation box     
        int m = neighnumsum[i+1] - neighnumsum[i];        // number of neighbors around i             
        int start = neighnumsum[i] - neighnumsum[istart];               
        T f1 = 0.0;
        T f2 = 0.0;        
        T f3 = 0.0;        
        // i is self and j is neighbor
        for (int l=0; l<m ; l++) {   // loop over each atom around j atom i -> pair ij 
            int k = start + l;             
            f1 +=  -fij[3*k+0];
            f2 +=  -fij[3*k+1];
            f3 +=  -fij[3*k+2];                         
        }                
        
        // need to determine bnumsum and index 
        int ii = i - istart;
        start = bnumsum[ii];   
        m = bnumsum[ii+1]-start;   // number of neighbors j around i  (j < i)           
        // j is self and i is neighbor
        for (int l=0; l<m ; l++) { // loop over each atom around j atom i -> pair ji 
            int k = index[start + l];     
            f1 +=  fij[3*k+0];
            f2 +=  fij[3*k+1];
            f3 +=  fij[3*k+2];             
        }
        fi[3*i+0] = f1;
        fi[3*i+1] = f2;
        fi[3*i+2] = f3;        
    }
}
template void cpuHalfAtomDecomposition3D(double*, double*, int*, int*, int*, int, int);
template void cpuHalfAtomDecomposition3D(float*, float*, int*, int*, int*, int, int);

template <typename T> void cpuFullAtomDecomposition3D(T *fi, T *fij, int *ilist, int *anumsum, int inum)
{    
    for (int ii=0; ii<inum; ii++) {  // for each atom i in the simulation box     
        int i = ilist[ii];       // atom i        
        T f1 = 0.0;
        T f2 = 0.0;        
        T f3 = 0.0;        
        int start = anumsum[ii];   
        int m = anumsum[ii+1]-start;        // number of neighbors around i             
        for (int l=0; l<m ; l++) {   // loop over each atom around atom i
            int k = start + l;                     
            f1 +=  -fij[3*k+0];
            f2 +=  -fij[3*k+1];
            f3 +=  -fij[3*k+2];             
        }
        fi[3*i+0] = f1;
        fi[3*i+1] = f2;
        fi[3*i+2] = f3;
    }
}
template void cpuFullAtomDecomposition3D(double*, double*, int*, int*, int);
template void cpuFullAtomDecomposition3D(float*, float*, int*, int*, int);

template <typename T> void cpuHalfAtomDecomposition3D(T *fi, T *fij, int *ilist, int *anumsum, 
        int *bnumsum, int *index, int inum)
{    
    for (int ii=0; ii<inum; ii++) {  // for each atom i in the simulation box     
        int i = ilist[ii];       // atom i
        T f1 = 0.0;
        T f2 = 0.0;        
        T f3 = 0.0;        
        int start = anumsum[ii];   
        int m = anumsum[ii+1]-start; // number of neighbors j around i  (j > i)             
        // i is self and j is neighbor
        for (int l=0; l<m ; l++) {   // loop over each atom around j atom i -> pair ij 
            int k = start + l;                     
            f1 +=  -fij[3*k+0];
            f2 +=  -fij[3*k+1];
            f3 +=  -fij[3*k+2];             
        }                
        
        // need to determine bnumsum and index
        start = bnumsum[ii];   
        m = bnumsum[ii+1]-start;   // number of neighbors j around i  (j < i)           
        // j is self and i is neighbor
        for (int l=0; l<m ; l++) { // loop over each atom around j atom i -> pair ji 
            int k = index[start + l];                     
            f1 +=  fij[3*k+0];
            f2 +=  fij[3*k+1];
            f3 +=  fij[3*k+2];             
        }
        fi[3*i+0] = f1;
        fi[3*i+1] = f2;
        fi[3*i+2] = f3;
    }
}
template void cpuHalfAtomDecomposition3D(double*, double*, int*, int*, int*, int*, int);
template void cpuHalfAtomDecomposition3D(float*, float*, int*, int*, int*, int*, int);

//********************************************************************************************//

template <typename T> void cpuGetNeighPairs(T *xij, T *qi, T *qj, T *x, T *q, int *ai, int *aj,  
      int *ti, int *tj, int *anum, int *anumsum, int *ilist, int *alist, int *neighlist, int *neighnum, 
      int *atomtype, int inum, int jnum, int ncq, int dim)
{    
    // form anum
    for (int ii=0; ii<inum; ii++) {
        int i = ilist[ii];       // atom i
        int m = neighnum[i];     // number of neighbors around i             
        anum[ii] = m;              
    }
    
    // a list contains the starting positions of the first neighbor 
    cpuCumsum(anumsum, anum, inum+1);             
    
    for (int ii=0; ii<inum; ii++) {  // for each atom i in the simulation box     
        int i = ilist[ii];       // atom i
        int itype = atomtype[i];
        int m = anum[ii];        // number of neighbors around i             
        int start = anumsum[ii];   
        for (int l=0; l<m ; l++) {   // loop over each atom around atom i
            int g = neighlist[l +jnum*i];  // ghost index of atom j  
            int j = alist[g];  // atom j
            int k = start + l;                                     
            ai[k]        = i;
            aj[k]        = j; // should be alist[j];         
            ti[k]        = itype;       
            tj[k]        = atomtype[j]; // should be atomtype[alist[j]];         
            for (int p=0; p<dim; p++) 
                xij[k*dim+p]   = x[g*dim+p] - x[i*dim+p];  // xj - xi            
            for (int p=0; p<ncq; p++) {                
                qi[k*ncq+p] = q[i*ncq+p];
                qj[k*ncq+p] = q[j*ncq+p];
            }                
        }
    }    
}
template void cpuGetNeighPairs(double*, double*, double*, double*, double*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int, int, int, int);
template void cpuGetNeighPairs(float*, float*, float*, float*, float*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int, int, int, int);

template <typename T> void cpuGetHalfNeighPairs(T *xij, T *qi, T *qj, T *x, T *q, int *ai, int *aj,  
      int *ti, int *tj, int *anum, int *anumsum, int *ilist, int *alist, int *neighlist, int *neighnum, 
      int *atomtype, int inum, int jnum, int ncq, int dim)
{    
    // form anum
    for (int ii=0; ii<inum; ii++) {
        int i = ilist[ii];       // atom i
        int m = neighnum[i];     // number of neighbors around i             
        anum[ii] = 0;              
        for (int l=0; l<m ; l++) {   // loop over each atom around atom i
            int g = neighlist[l + jnum*i];  // atom j           
            int j = alist[g];  // atom j
            if (i < j) 
                anum[ii] += 1;
        }                
    }        
    
    // a list contains the starting positions of the first neighbor 
    cpuCumsum(anumsum, anum, inum+1);             
    
    for (int ii=0; ii<inum; ii++) {  // for each atom i in the simulation box     
        int i = ilist[ii];       // atom i
        int itype = atomtype[i];
        int m = anum[ii];        // number of neighbors around i             
        int start = anumsum[ii];   
        int inc=0;
        for (int l=0; l<m ; l++) {   // loop over each atom around atom i
            int g = neighlist[l +jnum*i];  // ghost index of atom j  
            int j = alist[g];  // atom j
            if (i < j) {
                int k = start + inc;                         
                ai[k]        = i;
                aj[k]        = j; // should be alist[j];         
                ti[k]        = itype;       
                tj[k]        = atomtype[j]; // should be atomtype[alist[j]];                                        
                for (int p=0; p<dim; p++) 
                    xij[k*dim+p]   = x[g*dim+p] - x[i*dim+p];  // xj - xi                            
                for (int p=0; p<ncq; p++) {                
                    qi[k*ncq+p] = q[i*ncq+p];
                    qj[k*ncq+p] = q[j*ncq+p];
                }                
                inc += 1;
            }
        }
    }    
}
template void cpuGetHalfNeighPairs(double*, double*, double*, double*, double*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int, int, int, int);
template void cpuGetHalfNeighPairs(float*, float*, float*, float*, float*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int, int, int, int);

template <typename T> void cpuGetNeighPairs(T *xij, T *qi, T *qj, T *x, T *q, int *ai, int *aj,  
      int *ti, int *tj, int *anum, int *anumsum, int *ilist, int *alist, int *neighlist, int *neighnum, 
      int *atomtype, int typej, int inum, int jnum, int ncq, int dim)
{        
    // form anum
    for (int ii=0; ii<inum; ii++) {
        int i = ilist[ii];       // atom i
        int m = neighnum[i];     // number of neighbors around i             
        anum[ii] = 0;              
        for (int l=0; l<m ; l++) {   // loop over each atom around atom i
            int g = neighlist[l + jnum*i];  // atom j           
            int j = alist[g];  // atom j
            if (atomtype[j] == typej) 
                anum[ii] += 1;
        }                
    }
    
    // a list contains the starting positions of the first neighbor 
    cpuCumsum(anumsum, anum, inum+1);             
    
    for (int ii=0; ii<inum; ii++) {  // for each atom i in the simulation box     
        int i = ilist[ii];       // atom i
        int itype = atomtype[i];
        int m = anum[ii];        // number of neighbors around i             
        int start = anumsum[ii];   
        int inc = 0;
        for (int l=0; l<m ; l++) {   // loop over each atom around atom i
            int g = neighlist[l +jnum*i];  // ghost index of atom j  
            int j = alist[g];  // atom j
            if (atomtype[j] == typej)  {
                int k = start + inc;                         
                ai[k]        = i;
                aj[k]        = j; // should be alist[j];         
                ti[k]        = itype;       
                tj[k]        = atomtype[j]; // should be atomtype[alist[j]];                                        
                for (int p=0; p<dim; p++) 
                    xij[k*dim+p]   = x[g*dim+p] - x[i*dim+p];  // xj - xi                            
                for (int p=0; p<ncq; p++) {                
                    qi[k*ncq+p] = q[i*ncq+p];
                    qj[k*ncq+p] = q[j*ncq+p];
                }                
                inc += 1;
            }
        }
    }    
}
template void cpuGetNeighPairs(double*, double*, double*, double*, double*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int, int, int, int, int);
template void cpuGetNeighPairs(float*, float*, float*, float*, float*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int, int, int, int, int);

void cpuFindAtom(int *tlist, int* ilist, int *atomtype, int typei, int inum)
{
    for (int ii=0; ii<inum; ii++) {  // for each atom i in the simulation box     
        int i = ilist[ii];           // atom i
        tlist[ii] = (atomtype[i] == typei) ? 1 : 0;        
    }    
}

void cpuFindStart(int *start, int* slist, int inum)
{
    for (int ii=0; ii<inum; ii++) {        
        if ((ii==0) && (slist[ii]==1)) {
            start[0] = 0;
            start[1] = inum;
        }
        if ((ii==inum-1) && (slist[ii]==0)) {
            start[0] = inum;
            start[1] = 0;
        }    
        if (slist[ii]-slist[ii-1]==1) {
            start[0] = ii;
            start[1] = inum - ii;
        }    
    }
}
        
void cpuCopyAtIndex(int* tlist, int* ilist, int* index, int *start)
{
    for (int ii=0; ii<start[1]; ii++) {                  
        tlist[ii] = ilist[index[start[0]+ii]];
    }    
}

void cpuGetIlist(int *tlist, int *start, int *ilist, int *slist, int *index, int *atomtype, int typei, int inum)
{
    cpuFindAtom(tlist, ilist, atomtype, typei, inum);
 
    // sort list    
    cpuMergeSort(slist, index, tlist, inum);
        
    cpuFindStart(start, slist, inum);        
        
    cpuCopyAtIndex(tlist, ilist, index, start);
}

template <typename T> void cpuGetNeighPairs(T *xij, T *qi, T *qj, T *x, T *q, int *ai, int *aj, 
      int *ti, int *tj, int *anum, int *anumsum, int *ilist, int *tlist, int *alist, int *neighlist, int *neighnum, 
      int *atomtype, int *slist, int *start, int *index, int typei, int inum, int jnum, int ncq, int dim)
{    
    // form tlist from ilist
    cpuGetIlist(tlist, start, ilist, slist, index, atomtype, typei, inum);
            
    cpuGetNeighPairs(xij, qi, qj, x, q, ai, aj, ti, tj, anum, anumsum, tlist, alist,
            neighlist, neighnum, atomtype, start[1], jnum, ncq, dim);             
}
template void cpuGetNeighPairs(double*, double*, double*, double*, double*, int*, int*, int*, int*, int*, 
        int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int, int, int, int, int);
template void cpuGetNeighPairs(float*, float*, float*, float*, float*, int*, int*, int*, int*, int*, 
        int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int, int, int, int, int);

template <typename T> void cpuGetNeighPairs(T *xij, T *qi, T *qj, T *x, T *q, int *ai, int *aj, 
      int *ti, int *tj, int *anum, int *anumsum, int *ilist, int *tlist, int *alist, int *neighlist, int *neighnum, 
      int *atomtype, int *slist, int *start, int *index, int typei, int typej, int inum, int jnum, int ncq, int dim)
{    
    // form tlist from ilist
    cpuGetIlist(tlist, start, ilist, slist, index, atomtype, typei, inum);
            
    cpuGetNeighPairs(xij, qi, qj, x, q, ai, aj, ti, tj, tlist, alist, anum, anumsum,
            neighlist, neighnum, atomtype, typej, start[1], jnum, ncq, dim);             
}
template void cpuGetNeighPairs(double*, double*, double*, double*, double*, int*, int*, int*, int*, int*, 
        int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int, int, int, int, int, int);
template void cpuGetNeighPairs(float*, float*, float*, float*, float*, int*, int*, int*, int*, int*, 
        int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int, int, int, int, int, int);

template <typename T> void cpuGetNeighTriplets(T *xij, T *xik, T *qi, T *qj, T *qk, T *x, T *q, int *ai, int *aj,  
      int *ak, int *ti, int *tj, int *tk, int *anum, int *anumsum, int *ilist, int *alist, int *neighlist, int *neighnum, 
      int *atomtype, int inum, int jnum, int ncq, int dim)
{        
    // form anum
    for (int ii=0; ii<inum; ii++) {
        int i = ilist[ii];       // atom i
        int m = neighnum[i];     // number of neighbors around i             
        anum[ii] = 0;              
        for (int lj=0; lj<m ; lj++) {   // loop over each atom around atom i
            int gj = neighlist[lj + jnum*i];  // atom j           
            int j = alist[gj];  // atom j
            for (int lk=0; lk<m; lk++) {
                int gk = neighlist[lk + jnum*i];  // atom k
                int k = alist[gk];  // atom k
                if (j < k)
                    anum[ii] += 1;
            }            
        }                
    }
    
    // a list contains the starting positions of the first neighbor 
    cpuCumsum(anumsum, anum, inum+1);             
    
    for (int ii=0; ii<inum; ii++) {  // for each atom i in the simulation box     
        int i = ilist[ii];       // atom i
        int itype = atomtype[i];
        int m = anum[ii];        // number of neighbors around i             
        int start = anumsum[ii];   
        int inc = 0;
        for (int lj=0; lj<m ; lj++) {   // loop over each atom around atom i
            int gj = neighlist[lj +jnum*i];  // ghost index of atom j  
            int j = alist[gj];  // atom j
            for (int lk=0; lk<m; lk++) {
                int gk = neighlist[lk + jnum*i];  // atom k
                int k = alist[gk];  // atom k
                if (j < k) {
                    int n = start + inc;                         
                    for (int p=0; p<dim; p++) {
                        xij[n*dim+p]   = x[gj*dim+p] - x[i*dim+p];  // xj - xi            
                        xik[n*dim+p]   = x[gk*dim+p] - x[i*dim+p];  // xj - xi            
                    }
                    ai[n]        = i;
                    aj[n]        = j; // should be alist[j];         
                    ak[n]        = k; // should be alist[j];         
                    ti[n]        = itype;       
                    tj[n]        = atomtype[j]; // should be atomtype[alist[j]];                                        
                    tk[n]        = atomtype[k]; // should be atomtype[alist[j]];                                        
                    for (int p=0; p<ncq; p++) {                
                        qi[n*ncq+p] = q[i*ncq+p];
                        qj[n*ncq+p] = q[j*ncq+p];
                        qk[n*ncq+p] = q[k*ncq+p];
                    }                
                    inc += 1;
                }
            }                                    
        }
    }    
}
template void cpuGetNeighTriplets(double*, double*, double*, double*, double*, double*, double*, 
        int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int, int, int, int);
template void cpuGetNeighTriplets(float*, float*, float*, float*, float*, float*, float*, 
        int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int, int, int, int);

template <typename T> void cpuGetNeighTriplets(T *xij, T *xik, T *qi, T *qj, T *qk, T *x, T *q, int *ai, int *aj,  
      int *ak, int *ti, int *tj, int *tk, int *anum, int *anumsum, int *ilist, int *alist, int *neighlist, int *neighnum, 
      int *atomtype, int typej, int typek, int inum, int jnum, int ncq, int dim)
{        
    // form anum
    for (int ii=0; ii<inum; ii++) {
        int i = ilist[ii];       // atom i
        int m = neighnum[i];     // number of neighbors around i             
        anum[ii] = 0;              
        for (int lj=0; lj<m ; lj++) {   // loop over each atom around atom i
            int gj = neighlist[lj + jnum*i];  // atom j           
            int j = alist[gj];  // atom j
            for (int lk=0; lk<m; lk++) {
                int gk = neighlist[lk + jnum*i];  // atom k
                int k = alist[gk];  // atom k
                if ((atomtype[j] == typej) && (atomtype[k] == typek) && (j < k))
                    anum[ii] += 1;
            }            
        }                
    }
    
    // a list contains the starting positions of the first neighbor 
    cpuCumsum(anumsum, anum, inum+1);             
    
    for (int ii=0; ii<inum; ii++) {  // for each atom i in the simulation box     
        int i = ilist[ii];       // atom i
        int itype = atomtype[i];
        int m = anum[ii];        // number of neighbors around i             
        int start = anumsum[ii];   
        int inc = 0;
        for (int lj=0; lj<m ; lj++) {   // loop over each atom around atom i
            int gj = neighlist[lj +jnum*i];  // ghost index of atom j  
            int j = alist[gj];  // atom j
            for (int lk=0; lk<m; lk++) {
                int gk = neighlist[lk + jnum*i];  // atom k
                int k = alist[gk];  // atom k
                if ((atomtype[j] == typej) && (atomtype[k] == typek) && (j < k)) {
                    int n = start + inc;               
                    for (int p=0; p<dim; p++) {
                        xij[n*dim+p]   = x[gj*dim+p] - x[i*dim+p];  // xj - xi            
                        xik[n*dim+p]   = x[gk*dim+p] - x[i*dim+p];  // xj - xi            
                    }
                    ai[n]        = i;
                    aj[n]        = j; // should be alist[j];         
                    ak[n]        = k; // should be alist[j];         
                    ti[n]        = itype;       
                    tj[n]        = atomtype[j]; // should be atomtype[alist[j]];                                        
                    tk[n]        = atomtype[k]; // should be atomtype[alist[j]];                                        
                    for (int p=0; p<ncq; p++) {                
                        qi[n*ncq+p] = q[i*ncq+p];
                        qj[n*ncq+p] = q[j*ncq+p];
                        qk[n*ncq+p] = q[k*ncq+p];
                    }                
                    inc += 1;
                }
            }                                    
        }
    }    
}
template void cpuGetNeighTriplets(double*, double*, double*, double*, double*, double*, double*, 
        int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int, int, int, int, int, int);
template void cpuGetNeighTriplets(float*, float*, float*, float*, float*, float*, float*, 
        int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int, int, int, int, int, int);

template <typename T> void cpuGetNeighTriplets(T *xij, T *xik, T *qi, T *qj, T *qk, T *x, T *q, int *ai, int *aj,  
      int *ak, int *ti, int *tj, int *tk, int *anum, int *anumsum, int *ilist,  int *tlist, int *alist, int *neighlist, int *neighnum, 
      int *atomtype, int *slist, int *start, int *index, int typei, int typej, int typek, int inum, int jnum, int ncq, int dim)
{    
    // form tlist from ilist
    cpuGetIlist(tlist, start, ilist, slist, index, atomtype, typei, inum);
            
    cpuGetNeighTriplets(xij, xik, qi, qj,qk,  x, q, ai, aj, ak, ti, tj, tk, anum, anumsum, tlist, alist,
            neighlist, neighnum, atomtype, typej, typek, start[1], jnum, ncq, dim);             
}
template void cpuGetNeighTriplets(double*, double*, double*, double*, double*, double*, double*, 
        int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, 
        int, int, int, int, int, int, int);
template void cpuGetNeighTriplets(float*, float*, float*, float*, float*, float*, float*, 
        int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, 
        int, int, int, int, int, int, int);




#endif

