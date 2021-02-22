#ifndef __CPUNEIGHCOUNT
#define __CPUNEIGHCOUNT

int cpuFindAtomType(int *tlist, int* ilist, int *atomtype, int *p, int *q, int typei, int inum)
{
    // form binary array p
    for (int ii=0; ii<inum; ii++) {  // for each atom i in the simulation box     
        int i = ilist[ii];           // atom i
        p[ii] = (atomtype[i] == typei) ? 1 : 0;        
    }    
    
    // parallel prefix sum on p
    cpuCumsum(q, p, inum+1);
    
    // form tlist by ignoring indices at which p are zero
    for (int ii=0; ii<inum; ii++) 
        if (p[ii] == 1)                     
            tlist[q[ii+1]-1] = ilist[ii];        
    
    return q[inum];
    
}

template <typename T> void cpuNeighSingles(T *xi, T *qi, T *x, T *q, int *ai, 
      int *ti, int *ilist, int *atomtype, int inum, int ncq, int dim)
{        
    for (int ii=0; ii<inum; ii++) {  // for each atom i in the simulation box     
        int i = ilist[ii];       // atom i
        int itype = atomtype[i];
        ai[ii]        = i;
        ti[ii]        = itype;       
        for (int d=0; d<dim; d++) 
            xi[ii*dim+d] = x[i*dim+d];           
        for (int d=0; d<ncq; d++) 
            qi[ii*ncq+d] = q[i*ncq+d];        
    }    
}
template void cpuNeighSingles(double*, double*, double*, double*, int*, int*, int*, int*,  int, int, int);
template void cpuNeighSingles(float*, float*, float*, float*, int*, int*, int*, int*, int, int, int);

template <typename T> void cpuFullNeighPairList(int *pairnum, int *pairlist, T *x, T rcutsq, int *ilist, int *neighlist, 
        int *neighnum, int inum, int jnum, int dim)
{    
    for (int ii=0; ii<inum; ii++) {
        int i = ilist[ii];       // atom i
        int m = neighnum[i];     // number of neighbors around i                       
        int count = 0;              
        for (int l=0; l<m ; l++) {   // loop over each atom around atom i {
            int g = neighlist[l + jnum*i];
            T dij = 0.0;
            for (int d=0; d<dim; d++) {
                T xij = (x[g*dim+d] - x[i*dim+d]);  // xj - xi                        
                dij += dij + xij*xij;
            }
            if (dij<=rcutsq) {
                pairlist[l + jnum*ii] = g;  // atom j     
                count += 1;
            }
        }
        pairnum[ii] = count;
    }    
}
template void cpuFullNeighPairList(int *, int *, double*, double, int*, int*, int*, int, int, int);
template void cpuFullNeighPairList(int *, int *, float*, float, int*, int*, int*, int, int, int);

template <typename T> void cpuFullNeighPairList(int *pairnum, int *pairlist, T *x, T rcutsq,  int *atomtype, int *ilist, int *alist, 
        int *neighlist, int *neighnum, int inum, int jnum, int typej, int dim)
{        
    for (int ii=0; ii<inum; ii++) {
        int i = ilist[ii];       // atom i
        int m = neighnum[i];     // number of neighbors around i             
        int count = 0;              
        for (int l=0; l<m ; l++) {   // loop over each atom around atom i
            int g = neighlist[l + jnum*i];  // atom j           
            int j = alist[g];  // atom j
            if (atomtype[j] == typej) {
                T dij = 0.0;
                for (int d=0; d<dim; d++) {
                    T xij = (x[g*dim+d] - x[i*dim+d]);  // xj - xi                        
                    dij += dij + xij*xij;
                }
                if (dij<=rcutsq) {
                    pairlist[count + jnum*ii] = g;
                    count += 1;
                }
            }
        }   
        pairnum[ii] = count;
    }                        
}
template void cpuFullNeighPairList(int *, int *, double*, double, int*, int*, int*, int*, int*, int, int, int, int);
template void cpuFullNeighPairList(int *, int *, float*, float, int*, int*, int*, int*, int*, int, int, int, int);

template <typename T> void cpuHalfNeighPairList(int *pairnum, int *pairlist, T *x, T rcutsq, int *ilist, int *alist, 
        int *neighlist, int *neighnum, int inum, int jnum, int dim)
{        
    for (int ii=0; ii<inum; ii++) {
        int i = ilist[ii];       // atom i
        int m = neighnum[i];     // number of neighbors around i             
        int count = 0;              
        for (int l=0; l<m ; l++) {   // loop over each atom around atom i
            int g = neighlist[l + jnum*i];  // atom j           
            int j = alist[g];  // atom j
            if (i < j) {
                T dij = 0.0;
                for (int d=0; d<dim; d++) {
                    T xij = (x[g*dim+d] - x[i*dim+d]);  // xj - xi                        
                    dij += dij + xij*xij;
                }
                if (dij<=rcutsq) {
                    pairlist[count + jnum*ii] = g;
                    count += 1;                
                }
            }
        }   
        pairnum[ii] = count;
    }                        
}
template void cpuHalfNeighPairList(int *, int *, double*, double, int*, int*, int*, int*, int, int, int);
template void cpuHalfNeighPairList(int *, int *, float*, float, int*, int*, int*, int*, int, int, int);

template <typename T> void cpuHalfNeighPairList(int *pairnum, int *pairlist, T *x, T rcutsq, int *atomtype, int *ilist, 
        int *alist, int *neighlist, int *neighnum, int inum, int jnum, int typej, int dim)
{        
    for (int ii=0; ii<inum; ii++) {
        int i = ilist[ii];       // atom i
        int m = neighnum[i];     // number of neighbors around i             
        int count = 0;              
        for (int l=0; l<m ; l++) {   // loop over each atom around atom i
            int g = neighlist[l + jnum*i];  // atom j           
            int j = alist[g];  // atom j
            if ((atomtype[j] == typej) && (i < j)) {
                T dij = 0.0;
                for (int d=0; d<dim; d++) {
                    T xij = (x[g*dim+d] - x[i*dim+d]);  // xj - xi                        
                    dij += dij + xij*xij;
                }
                if (dij<=rcutsq) {
                    pairlist[count + jnum*ii] = g;
                    count += 1;
                }
            }
        }   
        pairnum[ii] = count;
    }                        
}
template void cpuHalfNeighPairList(int *, int *, double*, double, int*, int*, int*, int*, int*, int, int, int, int);
template void cpuHalfNeighPairList(int *, int *, float*, float, int*, int*, int*, int*, int*, int, int, int, int);

template <typename T> void cpuNeighPairs(T *xij, T *qi, T *qj, T *x, T *q, int *ai, int *aj,  
      int *ti, int *tj, int *pairnum, int *pairlist, int *pairnumsum, int *ilist, int *alist,  
      int *atomtype, int inum, int jnum, int ncq, int dim)
{        
    for (int ii=0; ii<inum; ii++) {  // for each atom i in the simulation box     
        int i = ilist[ii];       // atom i
        int itype = atomtype[i];
        int m = pairnum[ii];        // number of neighbors around i             
        int start = pairnumsum[ii];   
        for (int l=0; l<m ; l++) {   // loop over each atom around atom i
            int g = pairlist[l +jnum*ii];  // ghost index of atom j  
            int j = alist[g];  // atom j
            int k = start + l;                                     
            ai[k]        = i;
            aj[k]        = j;          
            ti[k]        = itype;       
            tj[k]        = atomtype[j];        
            for (int d=0; d<dim; d++) 
                xij[k*dim+d]   = x[g*dim+d] - x[i*dim+d];  // xj - xi            
            for (int d=0; d<ncq; d++) {                
                qi[k*ncq+d] = q[i*ncq+d];
                qj[k*ncq+d] = q[j*ncq+d];
            }                
        }
    }    
}
template void cpuNeighPairs(double*, double*, double*, double*, double*, int*, int*, int*, int*, 
        int*, int*, int*, int*, int*, int*, int, int, int, int);
template void cpuNeighPairs(float*, float*, float*, float*, float*, int*, int*, int*, int*, 
        int*, int*, int*, int*, int*, int*, int, int, int, int);

// void cpuFullNeighTripletList(int *tripletnum, int *tripletlist, int *atomtype, int *ilist, int *alist, 
//         int *neighlist, int *neighnum, int inum, int jnum, int jknum, int typej, int typek)
// {        
//     for (int ii=0; ii<inum; ii++) {
//         int i = ilist[ii];       // atom i
//         int m = neighnum[i];     // number of neighbors around i             
//         int count = 0;                  
//         for (int lj=0; lj<m ; lj++) {   // loop over each atom j around atom i
//             int gj = neighlist[lj + jnum*i];  // atom j           
//             int j = alist[gj];  // atom j
//             if ((atomtype[j] == typej) || (atomtype[j] == typek)) {           
//                 for (int lk=0; (lk<m) && (lk != lj); lk++) { // loop over each atom k around atom i (k != j)
//                     int gk = neighlist[lk + jnum*i];  // atom k
//                     int k = alist[gk];  // atom k                    
//                     if ((atomtype[k] == typej) || (atomtype[k] == typek)) {                    
//                         tripletlist[0 + 2*count + jknum*ii] = gj;
//                         tripletlist[1 + 2*count + jknum*ii] = gk;
//                         count += 1;
//                     }
//                 }            
//             }
//         }                 
//         tripletnum[ii] = count;                                
//     }
// }

template <typename T> void cpuNeighTripletList(int *tripletnum, int *tripletlist, T *x, T rcutsq, int *atomtype, int *ilist,
        int *alist, int *neighlist, int *neighnum, int inum, int jnum, int jknum, int typej, int typek, int dim)
{        
    for (int ii=0; ii<inum; ii++) {
        int i = ilist[ii];       // atom i
        int m = neighnum[i];     // number of neighbors around i             
        int count = 0;                  
        for (int lj=0; lj<m ; lj++) {   // loop over each atom j around atom i
            int gj = neighlist[lj + jnum*i];  // atom j           
            int j = alist[gj];  // atom j
            if ((atomtype[j] == typej) || (atomtype[j] == typek)) {
                T dij = 0.0;
                for (int d=0; d<dim; d++) {
                    T xij = (x[gj*dim+d] - x[i*dim+d]);  // xj - xi                        
                    dij += dij + xij*xij;
                }
                if (dij <= rcutsq) {
                    for (int lk=lj+1; lk<m; lk++) { // loop over each atom k around atom i (k > j)
                        int gk = neighlist[lk + jnum*i];  // atom k
                        int k = alist[gk];  // atom k                    
                        if ((atomtype[k] == typej) || (atomtype[k] == typek)) {        
                            T dik = 0.0;
                            for (int d=0; d<dim; d++) {
                                T xik = (x[gk*dim+d] - x[i*dim+d]);  // xj - xi                        
                                dik += dij + xik*xik;
                            }
                            if (dij <= rcutsq) {
                                tripletlist[0 + 2*count + jknum*ii] = gj;
                                tripletlist[1 + 2*count + jknum*ii] = gk;
                                count += 1;
                            }
                        }
                    }            
                }
            }
        }                 
        tripletnum[ii] = count;                                
    }
}
template void cpuNeighTripletList(int *, int *, double*, double, int*, int*, int*, int*, int*, int, int, int, int, int, int);
template void cpuNeighTripletList(int *, int *, float*, float, int*, int*, int*, int*, int*, int, int, int, int, int, int);

template <typename T> void cpuNeighTriplets(T *xij, T *xik, T *qi, T *qj, T *qk, T *x, T *q, int *ai, int *aj, int *ak,  
      int *ti, int *tj, int *tk, int *tripletnum, int *tripletlist, int *tripletnumsum, int *ilist, int *alist,  
      int *atomtype, int inum, int jnum, int jknum, int ncq, int dim)
{        
    for (int ii=0; ii<inum; ii++) {  // for each atom i in the simulation box     
        int i = ilist[ii];       // atom i
        int itype = atomtype[i];
        int m = tripletnum[ii];        // number of neighbors around i             
        int start = tripletnumsum[ii];   
        for (int l=0; l<m ; l++) {   // loop over each atom pair (j,k) around atom i
            int gj = tripletlist[0 + 2*l + jknum*ii];  // ghost index of atom j  
            int gk = tripletlist[1 + 2*l + jknum*ii];  // ghost index of atom k  
            int j = alist[gj];  // atom j
            int k = alist[gk];  // atom k
            int n = start + l;                                     
            ai[n]        = i;
            aj[n]        = j;    
            ak[n]        = k;    
            ti[n]        = itype;       
            tj[n]        = atomtype[j];     
            tk[n]        = atomtype[k];     
            for (int d=0; d<dim; d++) {
                xij[n*dim+d]   = x[gj*dim+d] - x[i*dim+d];  // xj - xi  
                xik[n*dim+d]   = x[gk*dim+d] - x[i*dim+d];  // xk - xi  
            }
            for (int d=0; d<ncq; d++) {                
                qi[n*ncq+d] = q[i*ncq+d];
                qj[n*ncq+d] = q[j*ncq+d];
                qk[n*ncq+d] = q[k*ncq+d];
            }                
        }
    }    
}
template void cpuNeighTriplets(double*, double*, double*, double*, double*, double*, double*, int*, int*, int*, int*, 
        int*, int*, int*, int*, int*, int*, int*, int*, int, int, int, int, int);
template void cpuNeighTriplets(float*, float*, float*, float*, float*, float*, float*, int*, int*, int*, int*, 
        int*, int*, int*, int*, int*, int*, int*, int*, int, int, int, int, int);

      
void cpuNeighCountFullTriplet(int *anum, int *ilist, int *neighnum, int inum)
{    
    for (int ii=0; ii<inum; ii++) {
        int i = ilist[ii];       // atom i
        int m = neighnum[i];     // number of neighbors around i             
        anum[ii] = m*m;              
    }    
}

void cpuNeighCountFullTriplet(int *anum, int *bnum, int *atomtype,  int *ilist, int *alist, int *neighlist, int *neighnum, 
        int typej, int inum, int jnum)
{        
    for (int ii=0; ii<inum; ii++) {
        int i = ilist[ii];       // atom i
        int m = neighnum[i];     // number of neighbors around i             
        int count = 0;    
        int countj = 0;    
        for (int lj=0; lj<m ; lj++) {   // loop over each atom j around atom i
            int gj = neighlist[lj + jnum*i];  // atom j           
            int j = alist[gj];  // atom j
            if (atomtype[j] == typej) {                
                count += m;  
                countj += 1;
            }
        }                
        anum[ii] = count;
        bnum[ii] = countj;                
    }
}

void cpuNeighCountFullPairTriplet(int *anum, int *bnum, int *cnum, int *atomtype, int *ilist, int *alist, int *neighlist, int *neighnum, 
        int typej, int typek, int inum, int jnum)
{        
    for (int ii=0; ii<inum; ii++) {
        int i = ilist[ii];       // atom i
        int m = neighnum[i];     // number of neighbors around i             
        int count = 0;                  
        for (int lj=0; lj<m ; lj++) {   // loop over each atom j around atom i
            int gj = neighlist[lj + jnum*i];  // atom j           
            int j = alist[gj];  // atom j
            bool bj =  (atomtype[j] == typej);
            for (int lk=0; lk<m; lk++) {
                int gk = neighlist[lk + jnum*i];  // atom k
                int k = alist[gk];  // atom k
                bool bk =  (atomtype[k] == typek);
                if (bj && bk)
                    count += 1;
            }            
        }                
        anum[ii] = count;        
        bnum[ii] = 0;             
        cnum[ii] = 0;             
        for (int lj=0; lj<m ; lj++) {   // loop over each atom j around atom i
            int gj = neighlist[lj + jnum*i];  // atom j           
            int j = alist[gj];  // atom j
            if (atomtype[j] == typej)
                bnum[ii] += 1;
            if (atomtype[j] == typek)
                cnum[ii] += 1;
        }
    }
}

void cpuNeighCountFullTriplet(int *anum, int *bnum, int *cnum, int *atomtype, int *ilist, int *alist, int *neighlist, int *neighnum, 
        int typej, int typek, int inum, int jnum)
{        
    for (int ii=0; ii<inum; ii++) {
        int i = ilist[ii];       // atom i
        int m = neighnum[i];     // number of neighbors around i             
        int count = 0;                  
        for (int lj=0; lj<m ; lj++) {   // loop over each atom j around atom i
            int gj = neighlist[lj + jnum*i];  // atom j           
            int j = alist[gj];  // atom j
            bool bj =  ((atomtype[j] == typej) && (atomtype[j] == typek));
            for (int lk=0; lk<m; lk++) {
                int gk = neighlist[lk + jnum*i];  // atom k
                int k = alist[gk];  // atom k
                bool bk =  ((atomtype[k] == typej) && (atomtype[k] == typek));
                if (bj && bk)
                    count += 1;
            }            
        }                
        anum[ii] = count;        
        bnum[ii] = 0;             
        cnum[ii] = 0;             
        for (int lj=0; lj<m ; lj++) {   // loop over each atom j around atom i
            int gj = neighlist[lj + jnum*i];  // atom j           
            int j = alist[gj];  // atom j
            if (atomtype[j] == typej)
                bnum[ii] += 1;
            if (atomtype[j] == typek)
                cnum[ii] += 1;
        }
    }
}


void cpuNeighCountHalfTriplet(int *anum, int *ilist, int *alist, int *neighlist, int *neighnum, int inum, int jnum)
{        
    for (int ii=0; ii<inum; ii++) {
        int i = ilist[ii];       // atom i
        int m = neighnum[i];     // number of neighbors around i             
        int count = 0;              
        for (int lj=0; lj<m ; lj++) {   // loop over each atom around atom i
            int gj = neighlist[lj + jnum*i];  // atom j           
            int j = alist[gj];  // atom j
            for (int lk=0; lk<m; lk++) {
                int gk = neighlist[lk + jnum*i];  // atom k
                int k = alist[gk];  // atom k
                if (j < k)
                    count += 1;
            }            
        }                
        anum[ii] = count;
    }    
}

void cpuNeighCountHalfTriplet(int *anum, int *bnum, int *atomtype, int *ilist, int *alist, int *neighlist, int *neighnum, 
        int typej, int inum, int jnum)
{        
    for (int ii=0; ii<inum; ii++) {
        int i = ilist[ii];       // atom i
        int m = neighnum[i];     // number of neighbors around i             
        int count = 0;                  
        for (int lj=0; lj<m ; lj++) {   // loop over each atom j around atom i
            int gj = neighlist[lj + jnum*i];  // atom j           
            int j = alist[gj];  // atom j
            for (int lk=0; lk<m; lk++) {
                int gk = neighlist[lk + jnum*i];  // atom k
                int k = alist[gk];  // atom k
                if ((atomtype[j] == typej) && (j < k))
                    count += 1;
            }            
        }                
        anum[ii] = count;
        bnum[ii] = 0;                       
        for (int lj=0; lj<m ; lj++) {   // loop over each atom j around atom i
            int gj = neighlist[lj + jnum*i];  // atom j           
            int j = alist[gj];  // atom j
            if (atomtype[j] == typej)
                bnum[ii] += 1;
        }
    }
}

void cpuNeighCountHalfPairTriplet(int *anum, int *bnum, int *cnum, int *atomtype, int *ilist, int *alist, int *neighlist, int *neighnum, 
        int typej, int typek, int inum, int jnum)
{        
    for (int ii=0; ii<inum; ii++) {
        int i = ilist[ii];       // atom i
        int m = neighnum[i];     // number of neighbors around i             
        int count = 0;                  
        for (int lj=0; lj<m ; lj++) {   // loop over each atom j around atom i
            int gj = neighlist[lj + jnum*i];  // atom j           
            int j = alist[gj];  // atom j
            bool bj =  (atomtype[j] == typej);
            for (int lk=0; lk<m; lk++) {
                int gk = neighlist[lk + jnum*i];  // atom k
                int k = alist[gk];  // atom k
                bool bk = (atomtype[k] == typek);
                if ((bj) && (bk) && (j < k))
                    count += 1;
            }            
        }                
        anum[ii] = count;
        bnum[ii] = 0;             
        cnum[ii] = 0;             
        for (int lj=0; lj<m ; lj++) {   // loop over each atom j around atom i
            int gj = neighlist[lj + jnum*i];  // atom j           
            int j = alist[gj];  // atom j
            if (atomtype[j] == typej)
                bnum[ii] += 1;
            if (atomtype[j] == typek)
                cnum[ii] += 1;
        }
    }
}

void cpuNeighCountHalfTriplet(int *anum, int *bnum, int *cnum, int *atomtype, int *ilist, int *alist, int *neighlist, int *neighnum, 
        int typej, int typek, int inum, int jnum)
{        
    for (int ii=0; ii<inum; ii++) {
        int i = ilist[ii];       // atom i
        int m = neighnum[i];     // number of neighbors around i             
        int count = 0;                  
        for (int lj=0; lj<m ; lj++) {   // loop over each atom j around atom i
            int gj = neighlist[lj + jnum*i];  // atom j           
            int j = alist[gj];  // atom j
            bool bj =  ((atomtype[j] == typej) && (atomtype[j] == typek));
            for (int lk=0; lk<m; lk++) {
                int gk = neighlist[lk + jnum*i];  // atom k
                int k = alist[gk];  // atom k
                bool bk =  ((atomtype[k] == typej) && (atomtype[k] == typek));
                if ((bj) && (bk) && (j < k))
                    count += 1;
            }            
        }                
        anum[ii] = count;
        bnum[ii] = 0;             
        cnum[ii] = 0;             
        for (int lj=0; lj<m ; lj++) {   // loop over each atom j around atom i
            int gj = neighlist[lj + jnum*i];  // atom j           
            int j = alist[gj];  // atom j
            if (atomtype[j] == typej)
                bnum[ii] += 1;
            if (atomtype[j] == typek)
                cnum[ii] += 1;
        }
    }
}

#endif


