/***************************************************************************                               
                    Molecular Dynamics Potentials (MDP)
                           CESMIX-MIT Project  
 
 Contributing authors: Ngoc-Cuong Nguyen (cuongng@mit.edu, exapde@gmail.com)
 ***************************************************************************/

#ifndef MDP_OMPNEIGHBORLIST2D
#define MDP_OMPNEIGHBORLIST2D

#include <string>
#include <fstream>
#include <sstream>
#include <iostream>
#include <math.h>
#include <string.h>
using namespace std;

template <typename T> void ompBoundingBox2D(T *vc, T*wc, T *v, T *w, T *a, T *b, T *r, int *pbc)
{
    T a1 = a[0]; 
    T a2 = a[1];
    T b1 = b[0]; 
    T b2 = b[1];
    T norma = sqrt(ompArraySquareSum(a, 2));
    T normb = sqrt(ompArraySquareSum(b, 2));
    
    // vertices of the parallelogram defined by a and b
    v[0] = 0.0;
    v[1] = 0.0;
    v[2] = a1;
    v[3] = a2;
    v[4] = a1 + b1;
    v[5] = a2 + b2;
    v[6] = b1;
    v[7] = b2;
    
    // Find da such that the distance from the origin (0,0)
    // to the line a2*x - a1*y = da is equal to r
    T da = r[1]*norma;
    // Find db such that the distance from the origin (0,0)
    // to the line b2*x - b1*y = -db is equal to r
    T db = r[0]*normb;
    
    // intersection of a2*x - a1*y = da and b2*x - b1*y = -db
    T A[4] = {a2, b2, -a1, -b1};
    T invA[4];    
    ompSmallMatrixInverse(invA, A, 2);    
    w[0] = invA[0]*da - invA[2]*db;
    w[1] = invA[1]*da - invA[3]*db;
        
    // find e = (e1,e2) such that e1*a + e2*b = w1
    A[0] = a1;
    A[1] = a2;
    A[2] = b1;
    A[3] = b2;
    ompSmallMatrixInverse(invA, A, 2);    
    T e[2];    
    e[0] = invA[0]*w[0] + invA[2]*w[1];
    e[1] = invA[1]*w[0] + invA[3]*w[1];    
    
    // distance between w1 and e(1)*a
    T sb = sqrt((w[0]-e[0]*a[0])*(w[0]-e[0]*a[0]) + (w[1]-e[0]*a[1])*(w[1]-e[0]*a[1]));
    // distance between w1 and e(2)*a
    T sa = sqrt((w[0]-e[1]*b[0])*(w[0]-e[1]*b[0]) + (w[1]-e[1]*b[1])*(w[1]-e[1]*b[1]));
   
    // length of the bounding parallelogram along the 1st axis
    T l1 = norma + 2*sa*(pbc[0]);
    // length of the bounding parallelogram along the 2nd axis
    T l2 = normb + 2*sb*(pbc[1]);
    
    // the 1st vertex of  the bounding parallelepiped
    w[0] = pbc[0]*e[0]*a1 + pbc[1]*e[1]*b1;
    w[1] = pbc[0]*e[0]*a2 + pbc[1]*e[1]*b2;
    
    // the 2nd vertex of  the bounding parallelogram
    w[2] = w[0] + l1*a[0]/norma;
    w[3] = w[1] + l1*a[1]/norma;
    
    // the 3rd vertex of  the bounding parallelogram
    w[4] = v[4] - w[0];
    w[5] = v[5] - w[1];
    
    // the 4th vertex of  the bounding parallelogram
    w[6] = w[0] + l2*b[0]/normb;
    w[7] = w[1] + l2*b[1]/normb;    
    
    // bounding box in the reference domain
    for (int i=0; i<4; i++) {
        vc[2*i+0] = invA[0]*v[2*i+0] + invA[2]*v[2*i+1];
        vc[2*i+1] = invA[1]*v[2*i+0] + invA[3]*v[2*i+1];
        wc[2*i+0] = invA[0]*w[2*i+0] + invA[2]*w[2*i+1];
        wc[2*i+1] = invA[1]*w[2*i+0] + invA[3]*w[2*i+1]; 
    }        
}
template void ompBoundingBox2D(double*, double*, double*, double*, double*, double*, double*, int*);
template void ompBoundingBox2D(float*, float*, float*, float*, float*, float*, float*, int*);

template <typename T> int ompPeriodicImages2D(T *pimages, T *a, T *b, int *pbc)
{
    T px[2][9];
    
    px[0][0] = 0.0;             px[1][0] = 0.0;     
    px[0][1] = a[0];            px[1][1] = a[1];    
    px[0][2] =-a[0];            px[1][2] =-a[1];    
    px[0][3] = b[0];            px[1][3] = b[1];      
    px[0][4] =-b[0];            px[1][4] =-b[1];      
    px[0][5] = a[0]+b[0];       px[1][5] = a[1]+b[1]; 
    px[0][6] = a[0]-b[0];       px[1][6] = a[1]-b[1]; 
    px[0][7] =-a[0]+b[0];       px[1][7] =-a[1]+b[1]; 
    px[0][8] =-a[0]-b[0];       px[1][8] =-a[1]-b[1]; 
    
    int nimages;
    int indx[9];
    if ((pbc[0]==0) && (pbc[1]==0)) {
         nimages = 1;
         int ind[] = {0};
         for (int i=0; i<nimages; i++)
            indx[i] = ind[i];
    }
    else if ((pbc[0]==1) && (pbc[1]==0)) {
        nimages = 3;
        int ind[] = {0,1,2};      
        for (int i=0; i<nimages; i++)
            indx[i] = ind[i];
    }
    else if ((pbc[0]==0) && (pbc[1]==1)) {
        nimages = 3;
        int ind[] = {0,3,4};    
        for (int i=0; i<nimages; i++)
            indx[i] = ind[i];
    }
    else if ((pbc[0]==1) && (pbc[1]==1)) {
        nimages = 9;
        int ind[] = {0,1,2,3,4,5,6,7,8};
        for (int i=0; i<nimages; i++)
            indx[i] = ind[i];
    }
    else {
    }
    
    for (int i=0; i<nimages; i++) {
        pimages[2*i+0] = px[0][indx[i]];
        pimages[2*i+1] = px[1][indx[i]];
    }    
    
    return nimages;
}
template int ompPeriodicImages2D(double*, double*, double*, int*);
template int ompPeriodicImages2D(float*, float*, float*, int*);
    
template <typename T> void ompMakeReferenceGrid(T *eta, T smin, T smax, T ds, int nc)
{
    eta[0] = smin;
    eta[nc] = smax;
    // ds = 1.0/(nc-2)
    #pragma omp parallel for
    for (int i=1; i<nc; i++)                   
        eta[i] = (i-1)*ds; // (nc-2)*1.0/(nc-2) = 1.0
}
template void ompMakeReferenceGrid(double*, double, double, double, int);
template void ompMakeReferenceGrid(float*, float, float, float, int);

void ompGridColoring2D(int *cellcolor, int *nc, int *bsize, int dim)
{
    int nc1 = nc[0];
    int nc2 = nc[1];
    int bs1 = bsize[0];
    int bs2 = bsize[1];
    int nb1 = ceil(((double) nc1)/((double) bs1)); // number of blocks
    int nb2 = ceil(((double) nc2)/((double) bs2)); // number of blocks

    int bid[bs1*bs2];
    for (int i=0; i<bs1*bs2; i++)            
         bid[i] = i;

    for (int i=0; i<nb1; i++) {
        int i1, i2;
        if (i==1) 
            i1 = (i==1) ? 0 : (i1 + bs1 - 1);
        i2 = i1 + bs1 - 1;    
        if (i2 >= nc1)
            i2 = nc1 - 1;
        for (int j=0; j<nb2; j++) {
            int j1, j2;
            if (j == 1) 
                j1 = (j==1) ? 0 : (j1 + bs2 - 1);
            j2 = j1 + bs2 - 1;    
            if (j2 >= nc2)
                j2 = nc2 - 1;              
            for (int m2=j1; m2<=j2; m2++)
                for (int m1=i1; m1<=i2; m1++)                    
                    cellcolor[m1 + nc[0]*m2] = bid[m1-i1 + (m2-j1)*bs1];
        }
    }                        
}

/************************************************************************************/

template <typename T> void ompGhostAtoms2D(int *glistnum, int *inside, T *x, T *pimages, T *wc, T *s2rmap, int n, int pnum, int dim)
{    
    #pragma omp parallel for
    for (int i=0; i<n; i++) { // loop over each atom i inside the simulation box
        glistnum[i] = 0; // set the number of ghost atoms for atom i to 0
        for (int j=1; j<pnum; j++) { // loop over each periodic image of atom i
            T xj0 = x[i*dim+0] + pimages[j*dim+0];  // periodic image of x      
            T xj1 = x[i*dim+1] + pimages[j*dim+1];        
            T xc0 = s2rmap[0]*xj0 + s2rmap[2]*xj1;        // map it to the unit square      
            T xc1 = s2rmap[1]*xj0 + s2rmap[3]*xj1;                        
            /// check if the mapped point is inside the bounding box
            if ((wc[0] <= xc0) && (xc0 <= wc[4]) &&  (wc[1] <= xc1) && (xc1 <= wc[5])) {
                glistnum[i] += 1;            // increase the number of ghost atoms by 1 
                inside[j-1 + i*(pnum-1)] = 1;
            }
            else
                inside[j-1 + i*(pnum-1)] = 0;
        }
    }              
}
template void ompGhostAtoms2D(int*, int*, double*, double*, double*, double*, int, int, int);
template void ompGhostAtoms2D(int*, int*, float*, float*, float*, float*, int, int, int);

template <typename T> void ompAtomList2D(int *alist, int *inside, int *glistnumsum, int *glistnum, 
        T *x, T *pimages, T *wc, T *s2rmap, int inum, int pnum, int dim)
{           
    // a list of ghost atoms
    ompGhostAtoms2D(glistnum, inside, x, pimages, wc, s2rmap, inum, pnum, dim);
    
    // a list contains the starting position of the ghost atom of every atom i
    ompCumsum(glistnumsum, glistnum, inum+1); 
    
    #pragma omp parallel for
    for (int i=0; i<inum; i++) { // loop over each atom i inside the simulation box            
        alist[i] = i;         // add atom i to the list        
        int q = inum + glistnumsum[i]; // offset the starting position by n
        int k = 0;            
        for (int j=1; j<pnum; j++) { // loop over each periodic image of atom i            
            if (inside[j-1 + i*(pnum-1)]) {                    
                x[dim*(q+k)+0] = x[i*dim+0] + pimages[j*dim+0]; // add the periodic image as a ghost atom
                x[dim*(q+k)+1] = x[i*dim+1] + pimages[j*dim+1]; //
                // q + k = inum + glistnumsum[i] + k
                alist[q+k] = i;  // add atom i to the list                
                k += 1;          // increase the number of ghost atoms for atom i by 1    
            }
        }
    }    
}
template void ompAtomList2D(int*, int*, int*, int*, double*, double*, double*, double*, int, int, int);
template void ompAtomList2D(int*, int*, int*, int*, float*, float*, float*, float*, int, int, int);

template <typename T> void ompCellList2D(int *clist, T *x, T *eta1, T *eta2, T *eta3, T *s2rmap, 
        int *nc, int inum, int natom, int dim)
{    
    #pragma omp parallel for
    for (int i=0; i<natom; i++) { // loop over each atom i
        int j1, j2, m;
        int k = dim*i;
        m = (i < inum) ? 1 : 0;
        // position of atom i in the unit square/cube
        T xt0 = s2rmap[0]*x[k] + s2rmap[2]*x[k+1];        
        T xt1 = s2rmap[1]*x[k] + s2rmap[3]*x[k+1];     
        // identify a cell containing atom i 
        for (j1=m; j1<nc[0]-m; j1++)
            if ((eta1[j1] <= xt0) && (xt0<= eta1[j1+1]))
                break;
        for (j2=m; j2<nc[1]-m; j2++) 
            if ((eta2[j2] <= xt1) && (xt1<= eta2[j2+1]))
                break;
        clist[i] = j1 + nc[0]*j2; // link that cell to atom i
    }                        
}
template void ompCellList2D(int*, double*, double*, double*, double*, double*, int*, int, int, int);
template void ompCellList2D(int*, float*, float*, float*, float*, float*, int*, int, int, int);

void getcellcounts(int *cellcounts, int* cell, int natom, int ncell)
{
    int i, j;
    #pragma omp parallel for
    for (i=0; i<natom; i++) {                    
        int cid = cell[i]; // cell id of atom i

        // fill with zeros   
        if (i==0) {
            cellcounts[0] = 0;
            for (j=0; j<cid; j++)             
                cellcounts[1+j] = 0;            
        }
        
        // fill with natom
        if (i==(natom-1))
            for (j=cid; j<ncell; j++)
                cellcounts[1+j] = natom;            

        if ((i>0) && (i<=natom-1)) {
            int lid = cell[i-1]; // cell id of atom i-1
            if (lid != cid)
                for (j = lid; j<cid; j++) 
                    cellcounts[1+j] = i;
        }
    }
}
        
void ompCell2AtomList(int *c2alist, int *c2anumsum, int *c2anum, int *clist, int natom, int ncell)
{       
//     for (int i = 0; i<natom; i++)
//         c2alist[i] = i;
    
    // sort atom-to-cell list    
    ompMergeSort(c2anum, c2alist, clist, natom); // c2alist store indices
    
    // count number of atoms for every cell
    getcellcounts(c2anumsum, c2anum, natom, ncell);     
}

template <typename T> void ompFullNeighborList2D(int *neighborlist, int *neighbornum, T *x, T *rcutsq, int *alist, 
        int *clist, int *c2alist, int *c2anumsum, int *nc, int inum, int jnum, int dim)
{                
    #pragma omp parallel for
    for (int ii=0; ii<inum; ii++) {  // for each atom i in the simulation box    
        int i = alist[ii];        
        T xi0 = x[i*dim];        // position of atom i
        T xi1 = x[i*dim+1];      // position of atom i            
        int j = clist[i];         // cell j contains atom i           
        int j1 = j%nc[0];
        int j2 = (j-j1)/nc[0];
        int inc = 0;                // increment
        for (int i1=-1; i1<=1; i1++) {
            int k1 = j1 + i1;
            for (int i2=-1; i2<=1; i2++) {
                int k2 = j2 + i2;
                int k = k1 + nc[0]*k2; // cell k                
                int s = c2anumsum[k];    // starting position of the first atom in cell k
                int m = c2anumsum[k+1]-s;  // number of atoms in cell k
                for (int l=0; l<m ; l++) { // for each atom j in cell k
                    j = c2alist[s+l];  // atom j
                    if (j != i) {                        
                        T xij0 = x[j*dim] - xi0;  // xj - xi
                        T xij1 = x[j*dim+1] - xi1; // xj - xi                        
                        T rijsq = xij0*xij0 + xij1*xij1;  // distance between atom i and atom j                       
                        if (rijsq <= rcutsq[0]) {                                                
                            neighborlist[inc + jnum*i] = j; // add atom j into the list
                            inc += 1;
                        }
                    }
                    if (inc==jnum)
                        break;                    
                }
                if (inc==jnum)
                    break;                
            }
            if (inc==jnum)
                break;            
        }    
        neighbornum[ii] = inc;
    }        
}
template void ompFullNeighborList2D(int*, int*, double*, double*, int*, int*, int*, int*, int*, int, int, int);
template void ompFullNeighborList2D(int*, int*, float*, float*, int*, int*, int*, int*, int*, int, int, int);

template <typename T> void ompFullNeighborList2D(int *neighborlist, int *neighbornum, T *x, T *rcutsq, 
        int anum, int inum, int jnum, int dim)
{                
    #pragma omp parallel for
    for (int i=0; i<inum; i++) {  // for each atom i in the simulation box    
        T xi0 = x[i*dim];        // position of atom i
        T xi1 = x[i*dim+1];      // position of atom i     
        int inc = 0;                // increment
        for (int j=0; j<anum; j++) { // loop over atom j (including both local and ghost atoms)            
            if (j != i) {
                T xij0 = x[j*dim] - xi0;  // xj - xi
                T xij1 = x[j*dim+1] - xi1; // xj - xi                
                T rijsq = xij0*xij0 + xij1*xij1; // distance between atom i and atom j                        
                if (rijsq <= rcutsq[0]) {       
                    neighborlist[inc + jnum*i] = j; // add atom j into the list
                    inc += 1;
                }            
            }
            if (inc==jnum)
                break;            
        }
        neighbornum[i] = inc;        
    }        
}
template void ompFullNeighborList2D(int*, int*, double*, double*, int, int, int, int);
template void ompFullNeighborList2D(int*, int*, float*, float*, int, int, int, int);

#endif


