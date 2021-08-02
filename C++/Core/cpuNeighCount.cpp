/***************************************************************************                               
                    Molecular Dynamics Potentials (MDP)
                           CESMIX-MIT Project  
 
 Contributing authors: Ngoc-Cuong Nguyen (cuongng@mit.edu, exapde@gmail.com)
 ***************************************************************************/

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

template <typename T> void cpuFullNeighPairList(int *pairnum, int *pairlist, T *x, T *rcutsq, int *ilist, int *neighlist, 
        int *neighnum, int inum, int jnum, int dim)
{    
    for (int ii=0; ii<inum; ii++) {
        int i = ilist[ii];       // atom i
        int m = neighnum[i];     // number of neighbors around i                       
        int count = 0;              
        for (int l=0; l<m ; l++) {   // loop over each atom around atom i {
            int g = neighlist[l + jnum*i];            
            // distance between atom i and atom j                                    
            T xij0 = x[g*dim] - x[i*dim];  // xj - xi
            T xij1 = x[g*dim+1] - x[i*dim+1]; // xj - xi               
            T dij = xij0*xij0 + xij1*xij1;
            if (dim==3)                 
                dij += (x[g*dim+2] - x[i*dim+2])*(x[g*dim+2] - x[i*dim+2]);              
            if (dij <= rcutsq[0]) {
                pairlist[count + jnum*ii] = g;  // atom j     
                count += 1;
            }
        }        
        pairnum[ii] = count;       
    }    
}
template void cpuFullNeighPairList(int *, int *, double*, double*, int*, int*, int*, int, int, int);
template void cpuFullNeighPairList(int *, int *, float*, float*, int*, int*, int*, int, int, int);

template <typename T> void cpuFullNeighPairList(int *pairnum, int *pairlist, T *x, T *rcutsq,  int *atomtype, int *ilist, int *alist, 
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
                // distance between atom i and atom j                                    
                T xij0 = x[g*dim] - x[i*dim];  // xj - xi
                T xij1 = x[g*dim+1] - x[i*dim+1]; // xj - xi               
                T dij = xij0*xij0 + xij1*xij1;
                if (dim==3)                 
                    dij += (x[g*dim+2] - x[i*dim+2])*(x[g*dim+2] - x[i*dim+2]);              
                if (dij<=rcutsq[0]) {
                    pairlist[count + jnum*ii] = g;
                    count += 1;
                }
            }
        }   
        pairnum[ii] = count;
    }                        
}
template void cpuFullNeighPairList(int *, int *, double*, double*, int*, int*, int*, int*, int*, int, int, int, int);
template void cpuFullNeighPairList(int *, int *, float*, float*, int*, int*, int*, int*, int*, int, int, int, int);

template <typename T> void cpuFullNeighPairList(int *pairnum, int *pairlist, T *x, T *rcutsq,  int *atomtype, int *ilist, int *alist, 
        int *neighlist, int *neighnum, int inum, int jnum, int typej, int typek, int dim)
{        
    for (int ii=0; ii<inum; ii++) {
        int i = ilist[ii];       // atom i
        int m = neighnum[i];     // number of neighbors around i             
        int count = 0;              
        for (int l=0; l<m ; l++) {   // loop over each atom around atom i
            int g = neighlist[l + jnum*i];  // atom j           
            int j = alist[g];  // atom j
            if ((atomtype[j] == typej) || (atomtype[j] == typek)) {
                // distance between atom i and atom j                                    
                T xij0 = x[g*dim] - x[i*dim];  // xj - xi
                T xij1 = x[g*dim+1] - x[i*dim+1]; // xj - xi               
                T dij = xij0*xij0 + xij1*xij1;
                if (dim==3)                 
                    dij += (x[g*dim+2] - x[i*dim+2])*(x[g*dim+2] - x[i*dim+2]);              
                if (dij<=rcutsq[0]) {
                    pairlist[count + jnum*ii] = g;
                    count += 1;
                }
            }
        }   
        pairnum[ii] = count;
    }                        
}
template void cpuFullNeighPairList(int *, int *, double*, double*, int*, int*, int*, int*, int*, int, int, int, int, int);
template void cpuFullNeighPairList(int *, int *, float*, float*, int*, int*, int*, int*, int*, int, int, int, int, int);

template <typename T> void cpuFullNeighPairList(int *pairnum, int *pairlist, T *x, T *rcutsq,  int *atomtype, int *ilist, int *alist, 
        int *neighlist, int *neighnum, int inum, int jnum, int typej, int typek, int typel, int dim)
{        
    for (int ii=0; ii<inum; ii++) {
        int i = ilist[ii];       // atom i
        int m = neighnum[i];     // number of neighbors around i             
        int count = 0;              
        for (int l=0; l<m ; l++) {   // loop over each atom around atom i
            int g = neighlist[l + jnum*i];  // atom j           
            int j = alist[g];  // atom j
            if ((atomtype[j] == typej) || (atomtype[j] == typek) || (atomtype[j] == typel)) {
                // distance between atom i and atom j                                    
                T xij0 = x[g*dim] - x[i*dim];  // xj - xi
                T xij1 = x[g*dim+1] - x[i*dim+1]; // xj - xi               
                T dij = xij0*xij0 + xij1*xij1;
                if (dim==3)                 
                    dij += (x[g*dim+2] - x[i*dim+2])*(x[g*dim+2] - x[i*dim+2]);              
                if (dij<=rcutsq[0]) {
                    pairlist[count + jnum*ii] = g;
                    count += 1;
                }
            }
        }   
        pairnum[ii] = count;
    }                        
}
template void cpuFullNeighPairList(int *, int *, double*, double*, int*, int*, int*, int*, int*, int, int, int, int, int, int);
template void cpuFullNeighPairList(int *, int *, float*, float*, int*, int*, int*, int*, int*, int, int, int, int, int, int);

template <typename T> void cpuHalfNeighPairList(int *pairnum, int *pairlist, T *x, T *rcutsq, int *ilist, int *alist, 
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
                // distance between atom i and atom j                                    
                T xij0 = x[g*dim] - x[i*dim];  // xj - xi
                T xij1 = x[g*dim+1] - x[i*dim+1]; // xj - xi               
                T dij = xij0*xij0 + xij1*xij1;
                if (dim==3)                 
                    dij += (x[g*dim+2] - x[i*dim+2])*(x[g*dim+2] - x[i*dim+2]);              
                if (dij<=rcutsq[0]) {
                    pairlist[count + jnum*ii] = g;
                    count += 1;                
                }
            }
        }   
        pairnum[ii] = count;
    }                        
}
template void cpuHalfNeighPairList(int *, int *, double*, double*, int*, int*, int*, int*, int, int, int);
template void cpuHalfNeighPairList(int *, int *, float*, float*, int*, int*, int*, int*, int, int, int);

template <typename T> void cpuHalfNeighPairList(int *pairnum, int *pairlist, T *x, T *rcutsq, int *atomtype, int *ilist, 
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
                // distance between atom i and atom j                                    
                T xij0 = x[g*dim] - x[i*dim];  // xj - xi
                T xij1 = x[g*dim+1] - x[i*dim+1]; // xj - xi               
                T dij = xij0*xij0 + xij1*xij1;
                if (dim==3)                 
                    dij += (x[g*dim+2] - x[i*dim+2])*(x[g*dim+2] - x[i*dim+2]);              
                if (dij<=rcutsq[0]) {
                    pairlist[count + jnum*ii] = g;
                    count += 1;
                }
            }
        }   
        pairnum[ii] = count;
    }                        
}
template void cpuHalfNeighPairList(int *, int *, double*, double*, int*, int*, int*, int*, int*, int, int, int, int);
template void cpuHalfNeighPairList(int *, int *, float*, float*, int*, int*, int*, int*, int*, int, int, int, int);

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

// template <typename T> void cpuNeighTripletList(int *tripletnum, int *tripletlist, T *x, T rcutsq, int *atomtype, int *ilist,
//         int *alist, int *neighlist, int *neighnum, int inum, int jnum, int jknum, int typej, int typek, int dim)
// {        
//     for (int ii=0; ii<inum; ii++) {
//         int i = ilist[ii];       // atom i
//         int m = neighnum[i];     // number of neighbors around i             
//         int count = 0;                  
//         for (int lj=0; lj<m ; lj++) {   // loop over each atom j around atom i
//             int gj = neighlist[lj + jnum*i];  // atom j           
//             int j = alist[gj];  // atom j
//             if ((atomtype[j] == typej) || (atomtype[j] == typek)) {
//                 T dij = 0.0;
//                 for (int d=0; d<dim; d++) {
//                     T xij = (x[gj*dim+d] - x[i*dim+d]);  // xj - xi                        
//                     dij += dij + xij*xij;
//                 }
//                 if (dij <= rcutsq) {
//                     for (int lk=lj+1; lk<m; lk++) { // loop over each atom k around atom i (k > j)
//                         int gk = neighlist[lk + jnum*i];  // atom k
//                         int k = alist[gk];  // atom k                    
//                         if ((atomtype[k] == typej) || (atomtype[k] == typek)) {        
//                             T dik = 0.0;
//                             for (int d=0; d<dim; d++) {
//                                 T xik = (x[gk*dim+d] - x[i*dim+d]);  // xj - xi                        
//                                 dik += dik + xik*xik;
//                             }
//                             if (dik <= rcutsq) {
//                                 tripletlist[0 + 2*count + jknum*ii] = gj;
//                                 tripletlist[1 + 2*count + jknum*ii] = gk;
//                                 count += 1;
//                             }
//                         }
//                     }            
//                 }
//             }
//         }                 
//         tripletnum[ii] = count;                                
//     }
// }
// template void cpuNeighTripletList(int *, int *, double*, double, int*, int*, int*, int*, int*, int, int, int, int, int, int);
// template void cpuNeighTripletList(int *, int *, float*, float, int*, int*, int*, int*, int*, int, int, int, int, int, int);

// template <typename T> void cpuNeighTriplets(T *xij, T *xik, T *qi, T *qj, T *qk, T *x, T *q, int *ai, int *aj, int *ak,  
//       int *ti, int *tj, int *tk, int *tripletnum, int *tripletlist, int *tripletnumsum, int *ilist, int *alist,  
//       int *atomtype, int inum, int jnum, int jknum, int ncq, int dim)
// {        
//     for (int ii=0; ii<inum; ii++) {  // for each atom i in the simulation box     
//         int i = ilist[ii];       // atom i
//         int itype = atomtype[i];
//         int m = tripletnum[ii];        // number of neighbors around i             
//         int start = tripletnumsum[ii];   
//         for (int l=0; l<m ; l++) {   // loop over each atom pair (j,k) around atom i
//             int gj = tripletlist[0 + 2*l + jknum*ii];  // ghost index of atom j  
//             int gk = tripletlist[1 + 2*l + jknum*ii];  // ghost index of atom k  
//             int j = alist[gj];  // atom j
//             int k = alist[gk];  // atom k
//             int n = start + l;                                     
//             ai[n]        = i;
//             aj[n]        = j;    
//             ak[n]        = k;    
//             ti[n]        = itype;       
//             tj[n]        = atomtype[j];     
//             tk[n]        = atomtype[k];     
//             for (int d=0; d<dim; d++) {
//                 xij[n*dim+d]   = x[gj*dim+d] - x[i*dim+d];  // xj - xi  
//                 xik[n*dim+d]   = x[gk*dim+d] - x[i*dim+d];  // xk - xi  
//             }
//             for (int d=0; d<ncq; d++) {                
//                 qi[n*ncq+d] = q[i*ncq+d];
//                 qj[n*ncq+d] = q[j*ncq+d];
//                 qk[n*ncq+d] = q[k*ncq+d];
//             }                
//         }
//     }    
// }
// template void cpuNeighTriplets(double*, double*, double*, double*, double*, double*, double*, int*, int*, int*, int*, 
//         int*, int*, int*, int*, int*, int*, int*, int*, int, int, int, int, int);
// template void cpuNeighTriplets(float*, float*, float*, float*, float*, float*, float*, int*, int*, int*, int*, 
//         int*, int*, int*, int*, int*, int*, int*, int*, int, int, int, int, int);

void cpuNeighTripletList(int *tripletlist, int *tripletnumsum, int *pairnum, 
        int *pairlist, int *ilist, int *alist, int inum, int jnum)
{        
    for (int ii=0; ii<inum; ii++) {
        int i = ilist[ii];       // atom i
        int p = pairnum[ii];      // number of pairs (i,j) around i         
        //int s = pairnumsum[ii];
        //tripletnum[ii] = (p-1)*p/2;       
        int q = tripletnumsum[ii];
        int count = 0;
        for (int lj=0; lj<p ; lj++) {   // loop over each atom j around atom i
            int gj = pairlist[lj + jnum*ii];  // atom j           
            int j = alist[gj];  // atom j    
            for (int lk=lj+1; lk<p; lk++) { // loop over each atom k around atom i (k > j)
                int gk = pairlist[lk + jnum*ii];  // atom k
                int k = alist[gk];  // atom k                    
                tripletlist[0 + 2*(count + q)] = gj;
                tripletlist[1 + 2*(count + q)] = gk;                
                count += 1;
            }
        }                        
    }
}

template <typename T> void cpuNeighTriplets(T *xij, T *xik, T *qi, T *qj, T *qk, T *x, T *q, int *ai, int *aj, int *ak,  
      int *ti, int *tj, int *tk, int *tripletnum, int *tripletlist, int *tripletnumsum, int *ilist, int *alist,  
      int *atomtype, int inum, int ncq, int dim)
{        
    for (int ii=0; ii<inum; ii++) {  // for each atom i in the simulation box     
        int i = ilist[ii];       // atom i
        int itype = atomtype[i];
        int m = tripletnum[ii];        // number of neighbors around i             
        int start = tripletnumsum[ii];   
        for (int l=0; l<m ; l++) {   // loop over each atom pair (j,k) around atom i
            int gj = tripletlist[0 + 2*(l + start)];  // ghost index of atom j  
            int gk = tripletlist[1 + 2*(l + start)];  // ghost index of atom k  
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
        int*, int*, int*, int*, int*, int*, int*, int*, int, int, int);
template void cpuNeighTriplets(float*, float*, float*, float*, float*, float*, float*, int*, int*, int*, int*, 
        int*, int*, int*, int*, int*, int*, int*, int*, int, int, int);

template <typename T> void cpuNeighTripletList(int *tripletnum, int *tripletlist, T *x, T *rcutsq, int *pairnum, int *pairnumsum, 
        int *pairlist, int *atomtype, int *ilist, int *alist, int *neighlist, int *neighnum, int inum, 
        int jnum, int typek, int dim)
{        
    for (int ii=0; ii<inum; ii++) {
        int i = ilist[ii];       // atom i
        int m = neighnum[i];     // number of neighbors around i                     
        int p = pairnum[ii];      // number of pairs (i,j) for atom i
        int s = pairnumsum[ii];
        for (int lj=0; lj<p ; lj++) {   // loop over each atom j around atom i
            int gj = pairlist[lj + jnum*ii];  // atom j           
            int j = alist[gj];  // atom j
            int count = 0;
            for (int lk=0; lk<m; lk++) { // loop over each atom k around atom i 
                int gk = neighlist[lk + jnum*i];  // atom k
                int k = alist[gk];  // atom k                    
                if ((atomtype[k] == typek) && (j != k)) {        
//                     T dik = 0.0;
//                     T djk = 0.0;
//                     for (int d=0; d<dim; d++) {
//                         T xik = (x[gk*dim+d] - x[i*dim+d]);  // xk - xi                        
//                         dik += dik + xik*xik;
//                         xik = (x[gk*dim+d] - x[gj*dim+d]);  // xk - xj                        
//                         djk += djk + xik*xik;
//                     }                                                        
                    T xij0 = x[gk*dim] - x[i*dim];  // xk - xi
                    T xij1 = x[gk*dim+1] - x[i*dim+1]; // xk - xi               
                    T dik = xij0*xij0 + xij1*xij1;
                    xij0 = x[gk*dim] - x[gj*dim];  // xk - xj
                    xij1 = x[gk*dim+1] - x[gj*dim+1]; // xk - xj               
                    T djk = xij0*xij0 + xij1*xij1;
                    if (dim==3) {                 
                        dik += (x[gk*dim+2] - x[i*dim+2])*(x[gk*dim+2] - x[i*dim+2]);              
                        djk += (x[gk*dim+2] - x[gj*dim+2])*(x[gk*dim+2] - x[gj*dim+2]);              
                    }
                    if ((dik <= rcutsq[0]) || (djk <= rcutsq[0]))  {
                        int nn = 3*(count + (lj + s)*jnum); // count < jnum
                        tripletlist[0 + nn] = i;
                        tripletlist[1 + nn] = gj;
                        tripletlist[2 + nn] = gk;
                        count += 1;
                    }
                }
            }
            tripletnum[lj+s] = count;                
        }                                         
    }
}
template void cpuNeighTripletList(int *, int *, double*, double*, int*, int*, int*, int*, int*, int*, int*, int*, int, int, int, int);
template void cpuNeighTripletList(int *, int *, float*, float*, int*, int*, int*, int*, int*, int*, int*, int*, int, int, int, int);

template <typename T> void cpuNeighTriplets(T *xij, T *xik, T *qi, T *qj, T *qk, T *x, T *q, int *ai, int *aj, int *ak,  
      int *ti, int *tj, int *tk, int *tripletnum, int *tripletlist, int *tripletnumsum, int *alist,  
      int *atomtype, int jnum, int ijnum, int ncq, int dim)
{        
    for (int ii=0; ii<ijnum; ii++) {  // loop over each pair (i,j)
        int m = tripletnum[ii];
        int start = tripletnumsum[ii];
        for (int lk=0; lk<m; lk++) { // loop over each atom k around pair (i, j) 
            int nn = 3*(lk + ii*jnum); 
            int i  = tripletlist[0 + nn];
            int gj = tripletlist[1 + nn];
            int gk = tripletlist[2 + nn];            
            int j = alist[gj];  // atom j
            int k = alist[gk];  // atom k
            int n = start + lk;                                     
            ai[n]        = i;
            aj[n]        = j;    
            ak[n]        = k;    
            ti[n]        = atomtype[i];       
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
        int*, int*, int*, int*, int*, int*, int*, int, int, int, int);
template void cpuNeighTriplets(float*, float*, float*, float*, float*, float*, float*, int*, int*, int*, int*, 
        int*, int*, int*, int*, int*, int*, int*, int, int, int, int);


// template <typename T> void cpuNeighQuadrupletList(int *quadrupletnum, int *quadrupletlist, T *x, T rcutsq, int *atomtype, int *ilist,
//         int *alist, int *neighlist, int *neighnum, int inum, int jnum, int jklnum, int typej, int typek, int typel, int dim)
// {        
//     for (int ii=0; ii<inum; ii++) {
//         int i = ilist[ii];       // atom i
//         int m = neighnum[i];     // number of neighbors around i             
//         int count = 0;                  
//         for (int lj=0; lj<m ; lj++) {   // loop over each atom j around atom i
//             int gj = neighlist[lj + jnum*i];  // atom j           
//             int j = alist[gj];  // atom j
//             if ((atomtype[j] == typej) || (atomtype[j] == typek) || (atomtype[j] == typel)) {
//                 T dij = 0.0;
//                 for (int d=0; d<dim; d++) {
//                     T xij = (x[gj*dim+d] - x[i*dim+d]);  // xj - xi                        
//                     dij += dij + xij*xij;
//                 }
//                 if (dij <= rcutsq) {
//                     for (int lk=lj+1; lk<m; lk++) { // loop over each atom k around atom i (k > j)
//                         int gk = neighlist[lk + jnum*i];  // atom k
//                         int k = alist[gk];  // atom k                    
//                         if ((atomtype[k] == typej) || (atomtype[k] == typek) || (atomtype[k] == typel)) {        
//                             T dik = 0.0;
//                             for (int d=0; d<dim; d++) {
//                                 T xik = (x[gk*dim+d] - x[i*dim+d]);  // xk - xi                        
//                                 dik += dij + xik*xik;
//                             }
//                             if (dik <= rcutsq) {
//                                 for (int ll=lk+1; ll<m; ll++) { // loop over each atom k around atom i (l > k)
//                                     int gl = neighlist[ll + jnum*i];  // atom k
//                                     int l = alist[gl];  // atom k                    
//                                     if ((atomtype[l] == typej) || (atomtype[l] == typek) || (atomtype[l] == typel)) {        
//                                         T dil = 0.0;
//                                         for (int d=0; d<dim; d++) {
//                                             T xil = (x[gl*dim+d] - x[i*dim+d]);  // xl - xi                        
//                                             dil += dij + xil*xil;
//                                         }
//                                         if (dil <= rcutsq) {                                
//                                             quadrupletlist[0 + 3*count + jklnum*ii] = gj;
//                                             quadrupletlist[1 + 3*count + jklnum*ii] = gk;
//                                             quadrupletlist[2 + 3*count + jklnum*ii] = gl;
//                                             count += 1;
//                                         }
//                                     }
//                                 }
//                             }
//                         }
//                     }            
//                 }
//             }
//         }                 
//         quadrupletnum[ii] = count;                                
//     }
// }
// template void cpuNeighQuadrupletList(int *, int *, double*, double, int*, int*, int*, int*, int*, int, int, int, int, int, int, int);
// template void cpuNeighQuadrupletList(int *, int *, float*, float, int*, int*, int*, int*, int*, int, int, int, int, int, int, int);
// 
// template <typename T> void cpuNeighQuadruplets(T *xij, T *xik, T *xil, T *qi, T *qj, T *qk, T *ql, T *x, T *q, 
//       int *ai, int *aj, int *ak, int *al, int *ti, int *tj, int *tk, int *tl, int *quadrupletnum, int *quadrupletlist, 
//       int *quadrupletnumsum, int *ilist, int *alist,  int *atomtype, int inum, int jnum, int jklnum, int ncq, int dim)
// {        
//     for (int ii=0; ii<inum; ii++) {  // for each atom i in the simulation box     
//         int i = ilist[ii];       // atom i
//         int itype = atomtype[i];
//         int m = quadrupletnum[ii];        // number of neighbors around i             
//         int start = quadrupletnumsum[ii];   
//         for (int p=0; p<m ; p++) {   // loop over each atom pair (j,k) around atom i
//             int gj = quadrupletlist[0 + 3*p + jklnum*ii];  // ghost index of atom j  
//             int gk = quadrupletlist[1 + 3*p + jklnum*ii];  // ghost index of atom k  
//             int gl = quadrupletlist[2 + 3*p + jklnum*ii];  // ghost index of atom l  
//             int j = alist[gj];  // atom j
//             int k = alist[gk];  // atom k
//             int l = alist[gl];  // atom l
//             int n = start + p;                                     
//             ai[n]        = i;
//             aj[n]        = j;    
//             ak[n]        = k;    
//             al[n]        = l;    
//             ti[n]        = itype;       
//             tj[n]        = atomtype[j];     
//             tk[n]        = atomtype[k];     
//             tl[n]        = atomtype[l];     
//             for (int d=0; d<dim; d++) {
//                 xij[n*dim+d]   = x[gj*dim+d] - x[i*dim+d];  // xj - xi  
//                 xik[n*dim+d]   = x[gk*dim+d] - x[i*dim+d];  // xk - xi  
//                 xil[n*dim+d]   = x[gl*dim+d] - x[i*dim+d];  // xk - xi  
//             }
//             for (int d=0; d<ncq; d++) {                
//                 qi[n*ncq+d] = q[i*ncq+d];
//                 qj[n*ncq+d] = q[j*ncq+d];
//                 qk[n*ncq+d] = q[k*ncq+d];
//                 ql[n*ncq+d] = q[l*ncq+d];
//             }                
//         }
//     }    
// }
// template void cpuNeighQuadruplets(double*, double*, double*, double*, double*, double*, double*, double*, double*, int*, int*, int*, int*, 
//         int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int, int, int, int, int);
// template void cpuNeighQuadruplets(float*, float*, float*, float*, float*, float*, float*, float*, float*, int*, int*, int*, int*, 
//         int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int, int, int, int, int);


void cpuNeighQuadrupletList(int *quadrupletlist, int *quadrupletnumsum, int *pairnum, 
        int *pairlist, int *ilist, int *alist, int inum, int jnum)
{        
    for (int ii=0; ii<inum; ii++) {
        int i = ilist[ii];       // atom i
        int p = pairnum[ii];      // number of pairs (i,j) around i         
        //int s = pairnumsum[ii];
        //quadrupletnum[ii] = (p-2)*(p-1)*p/6;       
        int q = quadrupletnumsum[ii];
        int count = 0;
        for (int lj=0; lj<p ; lj++) {   // loop over each atom j around atom i
            int gj = pairlist[lj + jnum*ii];  // atom j           
            int j = alist[gj];  // atom j    
            for (int lk=lj+1; lk<p; lk++) { // loop over each atom k around atom i (k > j)
                int gk = pairlist[lk + jnum*ii];  // atom k
                int k = alist[gk];  // atom k
                for (int ll=lk+1; ll<p; ll++) { // loop over each atom l around atom i (l > k)
                    int gl = pairlist[ll + jnum*ii];  // atom l
                    int l = alist[gl];  // atom k                
                    quadrupletlist[0 + 3*(count + q)] = gj;
                    quadrupletlist[1 + 3*(count + q)] = gk;                
                    quadrupletlist[2 + 3*(count + q)] = gl;                
                    count += 1;
                }
            }
        }                        
    }
}

template <typename T> void cpuNeighQuadruplets(T *xij, T *xik, T *xil, T *qi, T *qj, T *qk, T *ql, T *x, T *q, 
      int *ai, int *aj, int *ak, int *al, int *ti, int *tj, int *tk, int *tl, int *quadrupletnum, int *quadrupletlist, 
      int *quadrupletnumsum, int *ilist, int *alist,  int *atomtype, int inum, int ncq, int dim)
{        
    for (int ii=0; ii<inum; ii++) {  // for each atom i in the simulation box     
        int i = ilist[ii];       // atom i
        int itype = atomtype[i];
        int m = quadrupletnum[ii];        // number of neighbors around i             
        int start = quadrupletnumsum[ii];   
        for (int p=0; p<m ; p++) {   // loop over each atom pair (j,k) around atom i
            int gj = quadrupletlist[0 + 3*(p + start)];  // ghost index of atom j  
            int gk = quadrupletlist[1 + 3*(p + start)];  // ghost index of atom k  
            int gl = quadrupletlist[2 + 3*(p + start)];  // ghost index of atom l  
            int j = alist[gj];  // atom j
            int k = alist[gk];  // atom k
            int l = alist[gl];  // atom l
            int n = start + p;                                     
            ai[n]        = i;
            aj[n]        = j;    
            ak[n]        = k;    
            al[n]        = l;    
            ti[n]        = itype;       
            tj[n]        = atomtype[j];     
            tk[n]        = atomtype[k];     
            tl[n]        = atomtype[l];     
            for (int d=0; d<dim; d++) {
                xij[n*dim+d]   = x[gj*dim+d] - x[i*dim+d];  // xj - xi  
                xik[n*dim+d]   = x[gk*dim+d] - x[i*dim+d];  // xk - xi  
                xil[n*dim+d]   = x[gl*dim+d] - x[i*dim+d];  // xk - xi  
            }
            for (int d=0; d<ncq; d++) {                
                qi[n*ncq+d] = q[i*ncq+d];
                qj[n*ncq+d] = q[j*ncq+d];
                qk[n*ncq+d] = q[k*ncq+d];
                ql[n*ncq+d] = q[l*ncq+d];
            }                
        }
    }    
}
template void cpuNeighQuadruplets(double*, double*, double*, double*, double*, double*, double*, double*, double*, int*, int*, int*, int*, 
        int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int, int, int);
template void cpuNeighQuadruplets(float*, float*, float*, float*, float*, float*, float*, float*, float*, int*, int*, int*, int*, 
        int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int, int, int);


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


