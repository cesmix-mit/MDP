/***************************************************************************                               
                    Molecular Dynamics Potentials (MDP)
                           CESMIX-MIT Project  
 
 Contributing authors: Ngoc-Cuong Nguyen (cuongng@mit.edu, exapde@gmail.com)
 ***************************************************************************/

#ifndef __CPUNEIGHBORLIST2D
#define __CPUNEIGHBORLIST2D

#include <string>
#include <fstream>
#include <sstream>
#include <iostream>
#include <math.h>
#include <string.h>
using namespace std;

template <typename T> void cpuBoundingBox2D(T *vc, T*wc, T *v, T *w, T *a, T *b, T *r, int *pbc)
{
    T a1 = a[0]; 
    T a2 = a[1];
    T b1 = b[0]; 
    T b2 = b[1];
    T norma = sqrt(cpuArraySquareSum(a, 2));
    T normb = sqrt(cpuArraySquareSum(b, 2));
    
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
    cpuSmallMatrixInverse(invA, A, 2);    
    w[0] = invA[0]*da - invA[2]*db;
    w[1] = invA[1]*da - invA[3]*db;
        
    // find e = (e1,e2) such that e1*a + e2*b = w1
    A[0] = a1;
    A[1] = a2;
    A[2] = b1;
    A[3] = b2;
    cpuSmallMatrixInverse(invA, A, 2);    
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
template void cpuBoundingBox2D(double*, double*, double*, double*, double*, double*, double*, int*);
template void cpuBoundingBox2D(float*, float*, float*, float*, float*, float*, float*, int*);

template <typename T> int cpuPeriodicImages2D(T *pimages, T *a, T *b, int *pbc)
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
template int cpuPeriodicImages2D(double*, double*, double*, int*);
template int cpuPeriodicImages2D(float*, float*, float*, int*);
    
template <typename T> void cpuMakeReferenceGrid(T *eta, T smin, T smax, T ds, int nc)
{
    eta[0] = smin;
    eta[nc] = smax;
    // ds = 1.0/(nc-2)
    for (int i=1; i<nc; i++)                   
        eta[i] = (i-1)*ds; // (nc-2)*1.0/(nc-2) = 1.0
}
template void cpuMakeReferenceGrid(double*, double, double, double, int);
template void cpuMakeReferenceGrid(float*, float, float, float, int);

void cpuGridColoring2D(int *cellcolor, int *nc, int *bsize, int dim)
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

template <typename T> void cpuGhostAtoms2D(int *glistnum, int *inside, T *x, T *pimages, T *wc, T *s2rmap, int n, int pnum, int dim)
{    
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
template void cpuGhostAtoms2D(int*, int*, double*, double*, double*, double*, int, int, int);
template void cpuGhostAtoms2D(int*, int*, float*, float*, float*, float*, int, int, int);

template <typename T> void cpuAtomList2D(int *alist, int *inside, int *glistnumsum, int *glistnum, 
        T *x, T *pimages, T *wc, T *s2rmap, int inum, int pnum, int dim)
{           
    // a list of ghost atoms
    cpuGhostAtoms2D(glistnum, inside, x, pimages, wc, s2rmap, inum, pnum, dim);
    
    // a list contains the starting position of the ghost atom of every atom i
    cpuCumsum(glistnumsum, glistnum, inum+1); 
    
    for (int i=0; i<inum; i++) { // loop over each atom i inside the simulation box            
        alist[i] = i;         // add atom i to the list        
        int q = inum + glistnumsum[i]; // offset the starting position 
        int k = 0;            
        for (int j=1; j<pnum; j++) { // loop over each periodic image of atom i            
            if (inside[j-1 + i*(pnum-1)]) { // if the periodic image is in the extended domain                    
                x[dim*(q+k)+0] = x[i*dim+0] + pimages[j*dim+0]; // add the periodic image as a ghost atom
                x[dim*(q+k)+1] = x[i*dim+1] + pimages[j*dim+1]; //
                // q + k = inum + glistnumsum[i] + k
                alist[q+k] = i;  // add atom i to the list                
                k += 1;          // increase the number of ghost atoms for atom i by 1    
            }
        }
    }    
}
template void cpuAtomList2D(int*, int*, int*, int*, double*, double*, double*, double*, int, int, int);
template void cpuAtomList2D(int*, int*, int*, int*, float*, float*, float*, float*, int, int, int);

template <typename T> void cpuCellList2D(int *clist, T *x, T *eta1, T *eta2, T *eta3, T *s2rmap, 
        int *nc, int inum, int natom, int dim)
{    
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
template void cpuCellList2D(int*, double*, double*, double*, double*, double*, int*, int, int, int);
template void cpuCellList2D(int*, float*, float*, float*, float*, float*, int*, int, int, int);

void getcellcounts(int *cellcounts, int* cell, int natom, int ncell)
{
    int i, j;
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
        
void cpuCell2AtomList(int *c2alist, int *c2anumsum, int *c2anum, int *clist, int natom, int ncell)
{       
//     for (int i = 0; i<natom; i++)
//         c2alist[i] = i;
    
    // sort atom-to-cell list    
    cpuMergeSort(c2anum, c2alist, clist, natom); // c2alist store indices
    
    // count number of atoms for every cell
    getcellcounts(c2anumsum, c2anum, natom, ncell);     
}

template <typename T> void cpuFullNeighborList2D(int *neighborlist, int *neighbornum, T *x, T *rcutsq, int *alist, 
        int *clist, int *c2alist, int *c2anumsum, int *nc, int inum, int jnum, int dim)
{                
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
template void cpuFullNeighborList2D(int*, int*, double*, double*, int*, int*, int*, int*, int*, int, int, int);
template void cpuFullNeighborList2D(int*, int*, float*, float*, int*, int*, int*, int*, int*, int, int, int);

template <typename T> void cpuFullNeighborList2D(int *neighborlist, int *neighbornum, T *x, T *rcutsq, 
        int anum, int inum, int jnum, int dim)
{                
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
template void cpuFullNeighborList2D(int*, int*, double*, double*, int, int, int, int);
template void cpuFullNeighborList2D(int*, int*, float*, float*, int, int, int, int);


// template <typename T> void cpuFullNeighborList2D(int *neighborlist, int *neighbornum, T *x, T *ellipsoid, int *alist, 
//         int *clist, int *c2alist, int *c2anumsum, int *nc, int inum, int jnum, int dim)
// {            
//     T onet = (T) 1.0;
//     T A00 = ellipsoid[0]; // ellipsoid for pair (i,j)   
//     T A10 = ellipsoid[1]; // ellipsoid for pair (i,j)
//     T A01 = ellipsoid[2]; // ellipsoid for pair (i,j)  
//     T A11 = ellipsoid[3]; // ellipsoid for pair (i,j)                                                          
//     for (int ii=0; ii<inum; ii++) {  // for each atom i in the simulation box    
//         int i = alist[ii];
//         neighbornum[i] = 0;
//         T xi0 = x[i*dim];        // position of atom i
//         T xi1 = x[i*dim+1];      // position of atom i            
//         int j = clist[i];         // cell j contains atom i           
//         int j1 = j%nc[0];
//         int j2 = (j-j1)/nc[0];
//         int inc = 0;                // increment
//         for (int i1=-1; i1<=1; i1++) {
//             int k1 = j1 + i1;
//             for (int i2=-1; i2<=1; i2++) {
//                 int k2 = j2 + i2;
//                 int k = k1 + nc[0]*k2; // cell k
//                 //int m = c2anum[k];     // number of atoms in cell k
//                 int s = c2anumsum[k];    // starting position of the first atom in cell k
//                 int m = c2anumsum[k+1]-s;  // number of atoms in cell k
//                 for (int l=0; l<m ; l++) {
//                     j = c2alist[s+l];  // atom j
//                     T xij0 = x[j*dim] - xi0;  // xj - xi
//                     T xij1 = x[j*dim+1] - xi1; // xj - xi
//                     // distance between atom i and atom j 
//                     T rij = xij0*(A00*xij0 + A01*xij1) + xij1*(A10*xij0 + A11*xij1);                        
//                     if (rij <= onet) {                                                
//                         neighborlist[inc + jnum*i] = j; // add atom j into the list
//                         inc += 1;
//                     }
//                 }
//             }
//         }    
//         neighbornum[ii] = inc;
//     }        
// }
// template void cpuFullNeighborList2D(int*, int*, double*, double*, int*, int*, int*, int*, int*, int, int, int);
// template void cpuFullNeighborList2D(int*, int*, float*, float*, int*, int*, int*, int*, int*, int, int, int);

// template <typename T> void cpuFullNeighborList2D(int *neighborlist, int *neighbornum, T *x, T *rcutsq, int *alist, 
//         int *clist, int *c2alist, int *c2anumsum, int *atomtype, int *nc, int ntype, int inum, int jnum, int dim)
// {            
//     for (int ii=0; ii<inum; ii++) {  // for each atom i in the simulation box    
//         int i = alist[ii];
//         int ti = atomtype[i];
//         neighbornum[i] = 0;
//         T xi0 = x[i*dim];        // position of atom i
//         T xi1 = x[i*dim+1];      // position of atom i            
//         int j = clist[i];         // cell j contains atom i           
//         int j1 = j%nc[0];
//         int j2 = (j-j1)/nc[0];
//         int inc = 0;                // increment
//         for (int i1=-1; i1<=1; i1++) {
//             int k1 = j1 + i1;
//             for (int i2=-1; i2<=1; i2++) {
//                 int k2 = j2 + i2;
//                 int k = k1 + nc[0]*k2; // cell k
//                 //int m = c2anum[k];     // number of atoms in cell k
//                 int s = c2anumsum[k];    // starting position of the first atom in cell k
//                 int m = c2anumsum[k+1]-s;  // number of atoms in cell k
//                 for (int l=0; l<m ; l++) {
//                     j = c2alist[s+l];  // atom j
//                     int tj = atomtype[alist[j]];
//                     T xij0 = x[j*dim] - xi0;  // xj - xi
//                     T xij1 = x[j*dim+1] - xi1; // xj - xi
//                     // distance between atom i and atom j 
//                     T rijsq = xij0*xij0 + xij1*xij1;                        
//                     if (rijsq <= rcutsq[ti + tj*ntype]) {                                                
//                         neighborlist[inc + jnum*i] = j; // add atom j into the list
//                         inc += 1;
//                     }
//                 }
//             }
//         }    
//         neighbornum[ii] = inc;
//     }        
// }
// template void cpuFullNeighborList2D(int*, int*, double*, double*, int*, int*, int*, int*, int*, int*, int, int, int, int);
// template void cpuFullNeighborList2D(int*, int*, float*, float*, int*, int*, int*, int*, int*, int*, int, int, int, int);
// 
// template <typename T> void cpuHalfNeighborList2D(int *neighborlist, int *neighbornum, T *x, T *ellipsoid, int *alist, 
//         int *clist, int *c2alist, int *c2anumsum, int *nc, int inum, int jnum, int dim)
// {            
//     T onet = (T) 1.0;
//     T A00 = ellipsoid[0]; // ellipsoid for pair (i,j)   
//     T A10 = ellipsoid[1]; // ellipsoid for pair (i,j)
//     T A01 = ellipsoid[2]; // ellipsoid for pair (i,j)  
//     T A11 = ellipsoid[3]; // ellipsoid for pair (i,j)                                                          
//     for (int ii=0; ii<inum; ii++) {  // for each atom i in the simulation box    
//         int i = alist[ii];
//         neighbornum[i] = 0;
//         T xi0 = x[i*dim];        // position of atom i
//         T xi1 = x[i*dim+1];      // position of atom i            
//         int j = clist[i];         // cell j contains atom i           
//         int j1 = j%nc[0];
//         int j2 = (j-j1)/nc[0];
//         int inc = 0;                // increment
//         for (int i1=-1; i1<=1; i1++) {
//             int k1 = j1 + i1;
//             for (int i2=-1; i2<=1; i2++) {
//                 int k2 = j2 + i2;
//                 int k = k1 + nc[0]*k2; // cell k
//                 //int m = c2anum[k];     // number of atoms in cell k
//                 int s = c2anumsum[k];    // starting position of the first atom in cell k
//                 int m = c2anumsum[k+1]-s;  // number of atoms in cell k
//                 for (int l=0; l<m ; l++) {
//                     j = c2alist[s+l];  // atom j
//                     if (i < alist[j]) { // atom i is connected to atom j only when i < j  
//                         T xij0 = x[j*dim] - xi0;  // xj - xi
//                         T xij1 = x[j*dim+1] - xi1; // xj - xi
//                         // distance between atom i and atom j 
//                         T rij = xij0*(A00*xij0 + A01*xij1) + xij1*(A10*xij0 + A11*xij1);                        
//                         if (rij <= onet) {                                                
//                             neighborlist[inc + jnum*i] = j; // add atom j into the list
//                             inc += 1;
//                         }
//                     }
//                 }
//             }
//         }    
//         neighbornum[ii] = inc;
//     }        
// }
// template void cpuHalfNeighborList2D(int*, int*, double*, double*, int*, int*, int*, int*, int*, int, int, int);
// template void cpuHalfNeighborList2D(int*, int*, float*, float*, int*, int*, int*, int*, int*, int, int, int);
// 
// template <typename T> void cpuHalfNeighborList2D(int *neighborlist, int *neighbornum, T *x, T *rcutsq, int *alist, 
//         int *clist, int *c2alist, int *c2anumsum, int *atomtype, int *nc, int ntype, int inum, int jnum, int dim)        
// {            
//     for (int ii=0; ii<inum; ii++) {  // for each atom i in the simulation box    
//         int i = alist[ii];
//         int ti = atomtype[i];
//         neighbornum[i] = 0;
//         T xi0 = x[i*dim];        // position of atom i
//         T xi1 = x[i*dim+1];      // position of atom i            
//         int j = clist[i];         // cell j contains atom i           
//         int j1 = j%nc[0];
//         int j2 = (j-j1)/nc[0];
//         int inc = 0;                // increment
//         for (int i1=-1; i1<=1; i1++) {
//             int k1 = j1 + i1;
//             for (int i2=-1; i2<=1; i2++) {
//                 int k2 = j2 + i2;
//                 int k = k1 + nc[0]*k2; // cell k
//                 //int m = c2anum[k];     // number of atoms in cell k
//                 int s = c2anumsum[k];    // starting position of the first atom in cell k
//                 int m = c2anumsum[k+1]-s;  // number of atoms in cell k
//                 for (int l=0; l<m ; l++) {
//                     j = c2alist[s+l];  // atom j                    
//                     if (i < alist[j]) { // atom i is connected to atom j only when i < j  
//                         T xij0 = x[j*dim] - xi0;  // xj - xi
//                         T xij1 = x[j*dim+1] - xi1; // xj - xi
//                         // distance between atom i and atom j 
//                         T rijsq = xij0*xij0 + xij1*xij1;                              
//                         if (rijsq <= rcutsq[ti + atomtype[alist[j]]*ntype]) {                                     
//                             neighborlist[inc + jnum*i] = j; // add atom j into the list
//                             inc += 1;
//                         }
//                     }
//                 }
//             }
//         }    
//         neighbornum[ii] = inc;
//     }        
// }
// template void cpuHalfNeighborList2D(int*, int*, double*, double*, int*, int*, int*, int*, int*, int*, int, int, int, int);
// template void cpuHalfNeighborList2D(int*, int*, float*, float*, int*, int*, int*, int*, int*, int*, int, int, int, int);

//     cpuAtomList2D(alist, inside, glistnumsum, glistnum,  x, pimages, wc, s2rmap, inum, pnum, dim);
//     cpuCellList2D(clist, x, eta1, eta2, eta3, s2rmap, nc, natom, dim);
//     cpuCell2AtomList(c2alist, c2anumsum, c2anum, clist, natom, ncell);
//     cpuFullNeighborList2D(ineighborlist, neighbornum, x, ellipsoid, alist, clist, c2alist, c2anumsum, 
//             nc, selfflag, inum, jnum, dim);
// 

// template <typename T> void cpuGetNeighPairs2D(T *xij, T *qi, T *qj, T *x, T *q, int *ai, int *aj,  
//       int *ti, int *tj, int *anum, int *anumsum, int *ilist, int *alist, int *neighlist, int *neighnum, 
//       int *atomtype, int inum, int jnum, int ncq, int dim)
// {    
//     // form anum
//     for (int ii=0; ii<inum; ii++) {
//         int i = ilist[ii];       // atom i
//         int m = neighnum[i];     // number of neighbors around i             
//         anum[ii] = m;              
//     }
//     
//     // a list contains the starting positions of the first neighbor 
//     cpuCumsum(anumsum, anum, inum+1);             
//     
//     for (int ii=0; ii<inum; ii++) {  // for each atom i in the simulation box     
//         int i = ilist[ii];       // atom i
//         T xi0 = x[i*dim];        // position of atom i
//         T xi1 = x[i*dim+1];      // position of atom i                    
//         int itype = atomtype[i];
//         int m = anum[ii];        // number of neighbors around i             
//         int start = anumsum[ii];   
//         for (int l=0; l<m ; l++) {   // loop over each atom around atom i
//             int g = neighlist[l +jnum*i];  // ghost index of atom j  
//             int j = alist[g];  // atom j
//             int k = start + l;                         
//             xij[k*dim]   = x[g*dim] - xi0;  // xj - xi
//             xij[k*dim+1] = x[g*dim+1] - xi1; // xj - xi     
//             ai[k]        = i;
//             aj[k]        = j; // should be alist[j];         
//             ti[k]        = itype;       
//             tj[k]        = atomtype[j]; // should be atomtype[alist[j]];                                        
//             for (int p=0; p<ncq; p++) {                
//                 qi[k*ncq+p] = q[i*ncq+p];
//                 qj[k*ncq+p] = q[j*ncq+p];
//             }                
//         }
//     }    
// }
// template void cpuGetNeighPairs2D(double*, double*, double*, double*, double*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int, int, int, int);
// template void cpuGetNeighPairs2D(float*, float*, float*, float*, float*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int, int, int, int);
// 
// // template <typename T> void cpuGetHalfNeighPairs2D(T *xij, T *qi, T *qj, T *x, T *q, int *ai, int *aj,  
// //       int *ti, int *tj, int *anum, int *anumsum, int *ilist, int *alist, int *neighlist, int *neighnum, 
// //       int *atomtype, int ncq, int inum, int jnum, int dim)
// // {    
// //     // form anum
// //     for (int ii=0; ii<inum; ii++) {
// //         int i = ilist[ii];       // atom i
// //         int m = neighnum[i];     // number of neighbors around i             
// //         anum[ii] = 0;              
// //         for (int l=0; l<m ; l++) {   // loop over each atom around atom i
// //             int g = neighlist[l + jnum*i];  // atom j           
// //             int j = alist[g];  // atom j
// //             if (i < j) 
// //                 anum[ii] += 1;
// //         }                
// //     }        
// //     
// //     // a list contains the starting positions of the first neighbor 
// //     cpuCumsum(anumsum, anum, inum+1);             
// //     
// //     for (int ii=0; ii<inum; ii++) {  // for each atom i in the simulation box     
// //         int i = ilist[ii];       // atom i
// //         T xi0 = x[i*dim];        // position of atom i
// //         T xi1 = x[i*dim+1];      // position of atom i                    
// //         int itype = atomtype[i];
// //         int m = anum[ii];        // number of neighbors around i             
// //         int start = anumsum[ii];   
// //         int inc=0;
// //         for (int l=0; l<m ; l++) {   // loop over each atom around atom i
// //             int g = neighlist[l +jnum*i];  // ghost index of atom j  
// //             int j = alist[g];  // atom j
// //             if (i < j) {
// //                 int k = start + inc;                         
// //                 xij[k*dim]   = x[g*dim] - xi0;  // xj - xi
// //                 xij[k*dim+1] = x[g*dim+1] - xi1; // xj - xi     
// //                 ai[k]        = i;
// //                 aj[k]        = j; // should be alist[j];         
// //                 ti[k]        = itype;       
// //                 tj[k]        = atomtype[j]; // should be atomtype[alist[j]];                                        
// //                 for (int p=0; p<ncq; p++) {                
// //                     qi[k*ncq+p] = q[i*ncq+p];
// //                     qj[k*ncq+p] = q[j*ncq+p];
// //                 }                
// //                 inc += 1;
// //             }
// //         }
// //     }    
// // }
// // template void cpuGetHalfNeighPairs2D(double*, double*, double*, double*, double*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int, int, int, int);
// // template void cpuGetHalfNeighPairs2D(float*, float*, float*, float*, float*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int, int, int, int);
// 
// template <typename T> void cpuGetNeighPairs2D(T *xij, T *qi, T *qj, T *x, T *q, int *ai, int *aj,  
//       int *ti, int *tj, int *anum, int *anumsum, int *ilist, int *alist, int *neighlist, int *neighnum, 
//       int *atomtype, int typej, int inum, int jnum, int ncq, int dim)
// {        
//     // form anum
//     for (int ii=0; ii<inum; ii++) {
//         int i = ilist[ii];       // atom i
//         int m = neighnum[i];     // number of neighbors around i             
//         anum[ii] = 0;              
//         for (int l=0; l<m ; l++) {   // loop over each atom around atom i
//             int g = neighlist[l + jnum*i];  // atom j           
//             int j = alist[g];  // atom j
//             if (atomtype[j] == typej) 
//                 anum[ii] += 1;
//         }                
//     }
//     
//     // a list contains the starting positions of the first neighbor 
//     cpuCumsum(anumsum, anum, inum+1);             
//     
//     for (int ii=0; ii<inum; ii++) {  // for each atom i in the simulation box     
//         int i = ilist[ii];       // atom i
//         T xi0 = x[i*dim];        // position of atom i
//         T xi1 = x[i*dim+1];      // position of atom i                    
//         int itype = atomtype[i];
//         int m = anum[ii];        // number of neighbors around i             
//         int start = anumsum[ii];   
//         int inc = 0;
//         for (int l=0; l<m ; l++) {   // loop over each atom around atom i
//             int g = neighlist[l +jnum*i];  // ghost index of atom j  
//             int j = alist[g];  // atom j
//             if (atomtype[j] == typej)  {
//                 int k = start + inc;                         
//                 xij[k*dim]   = x[g*dim] - xi0;  // xj - xi
//                 xij[k*dim+1] = x[g*dim+1] - xi1; // xj - xi     
//                 ai[k]        = i;
//                 aj[k]        = j; // should be alist[j];         
//                 ti[k]        = itype;       
//                 tj[k]        = atomtype[j]; // should be atomtype[alist[j]];                                        
//                 for (int p=0; p<ncq; p++) {                
//                     qi[k*ncq+p] = q[i*ncq+p];
//                     qj[k*ncq+p] = q[j*ncq+p];
//                 }                
//                 inc += 1;
//             }
//         }
//     }    
// }
// template void cpuGetNeighPairs2D(double*, double*, double*, double*, double*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int, int, int, int, int);
// template void cpuGetNeighPairs2D(float*, float*, float*, float*, float*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int, int, int, int, int);
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
// template <typename T> void cpuGetNeighPairs2D(T *xij, T *qi, T *qj, T *x, T *q, int *ai, int *aj, 
//       int *ti, int *tj, int *anum, int *anumsum, int *ilist, int *tlist, int *alist, int *neighlist, int *neighnum, 
//       int *atomtype, int *slist, int *start, int *index, int typei, int inum, int jnum, int ncq, int dim)
// {    
//     // form tlist from ilist
//     cpuGetIlist(tlist, start, ilist, slist, index, atomtype, typei, inum);
//             
//     cpuGetNeighPairs2D(xij, qi, qj, x, q, ai, aj, ti, tj, anum, anumsum, tlist, alist,
//             neighlist, neighnum, atomtype, start[1], jnum, ncq, dim);             
// }
// template void cpuGetNeighPairs2D(double*, double*, double*, double*, double*, int*, int*, int*, int*, int*, 
//         int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int, int, int, int, int);
// template void cpuGetNeighPairs2D(float*, float*, float*, float*, float*, int*, int*, int*, int*, int*, 
//         int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int, int, int, int, int);
// 
// template <typename T> void cpuGetNeighPairs2D(T *xij, T *qi, T *qj, T *x, T *q, int *ai, int *aj, 
//       int *ti, int *tj, int *anum, int *anumsum, int *ilist, int *tlist, int *alist, int *neighlist, int *neighnum, 
//       int *atomtype, int *slist, int *start, int *index, int typei, int typej, int inum, int jnum, int ncq, int dim)
// {    
//     // form tlist from ilist
//     cpuGetIlist(tlist, start, ilist, slist, index, atomtype, typei, inum);
//             
//     cpuGetNeighPairs2D(xij, qi, qj, x, q, ai, aj, ti, tj, tlist, alist, anum, anumsum,
//             neighlist, neighnum, atomtype, typej, start[1], jnum, ncq, dim);             
// }
// template void cpuGetNeighPairs2D(double*, double*, double*, double*, double*, int*, int*, int*, int*, int*, 
//         int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int, int, int, int, int, int);
// template void cpuGetNeighPairs2D(float*, float*, float*, float*, float*, int*, int*, int*, int*, int*, 
//         int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int, int, int, int, int, int);
// 
// template <typename T> void cpuGetNeighTriplets2D(T *xij, T *xik, T *qi, T *qj, T *qk, T *x, T *q, int *ai, int *aj,  
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
//         T xi0 = x[i*dim];        // position of atom i
//         T xi1 = x[i*dim+1];      // position of atom i                    
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
//                 if (j < k) {
//                     int n = start + inc;                         
//                     xij[n*dim]   = x[gj*dim] - xi0;  // xj - xi
//                     xij[n*dim+1] = x[gj*dim+1] - xi1; // xj - xi     
//                     xik[n*dim]   = x[gk*dim] - xi0;  // xj - xi
//                     xik[n*dim+1] = x[gk*dim+1] - xi1; // xj - xi     
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
// template void cpuGetNeighTriplets2D(double*, double*, double*, double*, double*, double*, double*, 
//         int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int, int, int, int);
// template void cpuGetNeighTriplets2D(float*, float*, float*, float*, float*, float*, float*, 
//         int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int, int, int, int);
// 
// template <typename T> void cpuGetNeighTriplets2D(T *xij, T *xik, T *qi, T *qj, T *qk, T *x, T *q, int *ai, int *aj,  
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
//         T xi0 = x[i*dim];        // position of atom i
//         T xi1 = x[i*dim+1];      // position of atom i                    
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
//                     xij[n*dim]   = x[gj*dim] - xi0;  // xj - xi
//                     xij[n*dim+1] = x[gj*dim+1] - xi1; // xj - xi     
//                     xik[n*dim]   = x[gk*dim] - xi0;  // xj - xi
//                     xik[n*dim+1] = x[gk*dim+1] - xi1; // xj - xi     
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
// template void cpuGetNeighTriplets2D(double*, double*, double*, double*, double*, double*, double*, 
//         int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int, int, int, int, int, int);
// template void cpuGetNeighTriplets2D(float*, float*, float*, float*, float*, float*, float*, 
//         int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int, int, int, int, int, int);
// 
// template <typename T> void cpuGetNeighTriplets2D(T *xij, T *xik, T *qi, T *qj, T *qk, T *x, T *q, int *ai, int *aj,  
//       int *ak, int *ti, int *tj, int *tk, int *anum, int *anumsum, int *ilist,  int *tlist, int *alist, int *neighlist, int *neighnum, 
//       int *atomtype, int *slist, int *start, int *index, int typei, int typej, int typek, int inum, int jnum, int ncq, int dim)
// {    
//     // form tlist from ilist
//     cpuGetIlist(tlist, start, ilist, slist, index, atomtype, typei, inum);
//             
//     cpuGetNeighTriplets2D(xij, xik, qi, qj,qk,  x, q, ai, aj, ak, ti, tj, tk, anum, anumsum, tlist, alist,
//             neighlist, neighnum, atomtype, typej, typek, start[1], jnum, ncq, dim);             
// }
// template void cpuGetNeighTriplets2D(double*, double*, double*, double*, double*, double*, double*, 
//         int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, 
//         int, int, int, int, int, int, int);
// template void cpuGetNeighTriplets2D(float*, float*, float*, float*, float*, float*, float*, 
//         int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, 
//         int, int, int, int, int, int, int);
// 
// 
// // template <typename T> void cpuFullNeighList2D(int *neighlist, int *neighnum, T *x, T* ellipsoid, 
// //         int *alist, int *verletlist, int *verletnum, int *glistnumsum, int inum, int jnum, int dim)
// // {
// //     T A00 = ellipsoid[0]; // ellipsoid for pair (i,j)   
// //     T A10 = ellipsoid[1]; // ellipsoid for pair (i,j)
// //     T A01 = ellipsoid[2]; // ellipsoid for pair (i,j)  
// //     T A11 = ellipsoid[3]; // ellipsoid for pair (i,j)                                                                              
// //     for (int ii=0; ii<inum; ii++) {  // for each atom i in the simulation box     
// //         int i = alist[ii];  // atom i
// //         T xi0 = x[i*dim];        // position of atom i
// //         T xi1 = x[i*dim+1];      // position of atom i              
// //         //int start = verletnumsum[i];     // starting 
// //         int m = verletnum[i]; // number of verlet atoms around i 
// //         int ninc = 0;                // increment 
// //         for (int l=0; l<m ; l++) {   // loop over each verlet atom around atom i
// //             int j = verletlist[l + jnum*i];  // verlet atom j  
// //             T xij0 = x[j*dim] - xi0;  // xj - xi
// //             T xij1 = x[j*dim+1] - xi1; // xj - xi
// //             // distance between atom i and atom j 
// //             T rij = xij0*(A00*xij0 + A01*xij1) + xij1*(A10*xij0 + A11*xij1);                        
// //             if (rij <= 1.0) {// if atom j is inside the ellipsoid                     
// //                 neighlist[ninc + jnum*i] = j;                   
// //                 ninc += 1;  // increase the number of neighbors by 1                                           
// //             }                            
// //             int g = glistnumsum[j+1] - glistnumsum[j]; // number of ghost atoms of atom j
// //             for (int k=0; k<g; k++) {
// //                 int jg = inum + glistnumsum[j] + k; // ghost atom jg of atom j
// //                 xij0 = x[jg*dim] - xi0;  // xj - xi
// //                 xij1 = x[jg*dim+1] - xi1; // xj - xi
// //                 // distance between atom i and atom jj 
// //                 rij = xij0*(A00*xij0 + A01*xij1) + xij1*(A10*xij0 + A11*xij1);                        
// //                 if (rij <= 1.0) {// if ghost atom jg is inside the ellipsoid                     
// //                     neighlist[ninc + jnum*i] = jg;                   
// //                     ninc += 1;  // increase the number of neighbors by 1                                           
// //                 }                
// //             }
// //         }
// //         neighnum[i] = ninc; // number of neighbors of atom i
// //     }
// // }
// // template void cpuFullNeighList2D(int*, int*, double*, double*, int*, int*, int*, int*, int, int, int);
// // template void cpuFullNeighList2D(int*, int*, float*, float*, int*, int*, int*, int*, int, int, int);
// // 
// // template <typename T> void cpuHalfNeighList2D(int *neighlist, int *neighnum, T *x, T* ellipsoid, 
// //         int *alist, int *verletlist, int *verletnum, int *glistnumsum, int inum, int jnum, int dim)
// // {
// //     T A00 = ellipsoid[0]; // ellipsoid for pair (i,j)   
// //     T A10 = ellipsoid[1]; // ellipsoid for pair (i,j)
// //     T A01 = ellipsoid[2]; // ellipsoid for pair (i,j)  
// //     T A11 = ellipsoid[3]; // ellipsoid for pair (i,j)                                                                              
// //     for (int ii=0; ii<inum; ii++) {  // for each atom i in the simulation box     
// //         int i = alist[ii];  // atom i
// //         T xi0 = x[i*dim];        // position of atom i
// //         T xi1 = x[i*dim+1];      // position of atom i              
// //         //int start = verletnumsum[i];     // starting 
// //         int m = verletnum[i]; // number of verlet atoms around i 
// //         int ninc = 0;                // increment 
// //         for (int l=0; l<m ; l++) {   // loop over each verlet atom around atom i
// //             int j = verletlist[l + jnum*i];  // verlet atom j  
// //             if (i <= j) { // atom i is connected to atom j only when i <= j  
// //                 T xij0 = x[j*dim] - xi0;  // xj - xi
// //                 T xij1 = x[j*dim+1] - xi1; // xj - xi
// //                 // distance between atom i and atom j 
// //                 T rij = xij0*(A00*xij0 + A01*xij1) + xij1*(A10*xij0 + A11*xij1);                        
// //                 if (rij <= 1.0) {// if atom j is inside the ellipsoid                     
// //                     neighlist[ninc + jnum*i] = j;                   
// //                     ninc += 1;  // increase the number of neighbors by 1                                           
// //                 }                            
// //                 int g = glistnumsum[j+1] - glistnumsum[j]; // number of ghost atoms of atom j
// //                 for (int k=0; k<g; k++) {
// //                     int gj = inum + glistnumsum[j] + k; // ghost atom jg of atom j
// //                     xij0 = x[gj*dim] - xi0;  // xj - xi
// //                     xij1 = x[gj*dim+1] - xi1; // xj - xi
// //                     // distance between atom i and atom jj 
// //                     rij = xij0*(A00*xij0 + A01*xij1) + xij1*(A10*xij0 + A11*xij1);                        
// //                     if (rij <= 1.0) {// if atom jj is inside the ellipsoid                     
// //                         neighlist[ninc + jnum*i] = gj;                   
// //                         ninc += 1;  // increase the number of neighbors by 1                                           
// //                     }                
// //                 }
// //             }
// //         }
// //         neighnum[i] = ninc; // number of neighbors of atom i
// //     }
// // }
// // template void cpuHalfNeighList2D(int*, int*, double*, double*, int*, int*, int*, int*, int, int, int);
// // template void cpuHalfNeighList2D(int*, int*, float*, float*, int*, int*, int*, int*, int, int, int);
// // 
// 
// // template <typename T> void cpuVerletAtoms2D(int *verletnum, T *x, T *ellipsoid, int *alist, 
// //         int *clist, int *c2alist, int *c2anumsum, int *nc, int inum, int dim)
// // {        
// //     T onet = (T) 1.0;
// //     T A00 = ellipsoid[0]; // ellipsoid for pair (i,j)   
// //     T A10 = ellipsoid[1]; // ellipsoid for pair (i,j)
// //     T A01 = ellipsoid[2]; // ellipsoid for pair (i,j)  
// //     T A11 = ellipsoid[3]; // ellipsoid for pair (i,j)                                                          
// //     for (int ii=0; ii<inum; ii++) {  // for each atom i in the simulation box    
// //         int i = alist[ii];
// //         verletnum[i] = 0;
// //         T xi0 = x[i*dim];        // position of atom i
// //         T xi1 = x[i*dim+1];      // position of atom i            
// //         int j = clist[i];         // cell j contains atom i           
// //         int j1 = j%nc[0];
// //         int j2 = (j-j1)/nc[0];
// //         for (int i1=-1; i1<=1; i1++) {
// //             int k1 = j1 + i1;
// //             for (int i2=-1; i2<=1; i2++) {
// //                 int k2 = j2 + i2;
// //                 int k = k1 + nc[0]*k2; // cell k
// //                 //int m = c2anum[k];     // number of atoms in cell k
// //                 int s = c2anumsum[k];    // starting position of the first atom in cell k
// //                 int m = c2anumsum[k+1]-s;  // number of atoms in cell k
// //                 for (int l=0; l<m ; l++) {
// //                     j = c2alist[s+l];  // atom j
// //                     T xij0 = x[j*dim] - xi0;  // xj - xi
// //                     T xij1 = x[j*dim+1] - xi1; // xj - xi
// //                     // distance between atom i and atom j 
// //                     T rij = xij0*(A00*xij0 + A01*xij1) + xij1*(A10*xij0 + A11*xij1);                        
// //                     if (rij <= onet)
// //                         verletnum[i] += 1;
// //                 }
// //             }
// //         }                
// //     }        
// // }
// // template void cpuVerletAtoms2D(int*, double*, double*, int*, int*, int*, int*, int*, int, int);
// // template void cpuVerletAtoms2D(int*, float*, float*, int*, int*, int*, int*, int*, int, int);
// // 
// // template <typename T> void cpuCreateVerletList2D(int *verletlist, T *x, T *ellipsoid, int *verletnum, 
// //         int *verletnumsum, int *alist, int *clist, int *c2alist, int *c2anumsum, int *nc, 
// //         int selfflag, int inum, int dim)
// // {    
// //     // a list of verlet atoms
// //     cpuVerletAtoms2D(verletnum, x, ellipsoid, alist, clist, c2alist, c2anumsum, nc, inum, dim);
// // 
// //     // a list contains the starting positions of the first j atom 
// //     cpuCumsum(verletnumsum, verletnum, inum+1); 
// //     
// //     T onet = (T) 1.0;
// //     T A00 = ellipsoid[0]; // ellipsoid for pair (i,j)   
// //     T A10 = ellipsoid[1]; // ellipsoid for pair (i,j)
// //     T A01 = ellipsoid[2]; // ellipsoid for pair (i,j)  
// //     T A11 = ellipsoid[3]; // ellipsoid for pair (i,j)                                                                  
// //     for (int ii=0; ii<inum; ii++) {  // for each atom i in the simulation box    
// //         int i = alist[ii];
// //         T xi0 = x[i*dim];        // position of atom i
// //         T xi1 = x[i*dim+1];      // position of atom i
// //         int nstart = verletnumsum[i];     // starting 
// //         int ninc = 0;                // increment
// //         int j = clist[i];         // cell j contains atom i           
// //         int j1 = j%nc[0];
// //         int j2 = (j-j1)/nc[0];
// //         for (int i1=-1; i1<=1; i1++) {
// //             int k1 = j1 + i1;
// //             for (int i2=-1; i2<=1; i2++) {
// //                 int k2 = j2 + i2;
// //                 int k = k1 + nc[0]*k2; // cell k
// //                 //int m = c2anum[k];     // number of atoms in cell k
// //                 int s = c2anumsum[k];  // starting position of the first atom in cell k
// //                 int m = c2anumsum[k+1]-s;  // number of atoms in cell k
// //                 for (int l=0; l<m ; l++) { //loop over each atom j in cell k
// //                     j = c2alist[s+l];  // atom j
// //                     T xij0 = x[j*dim] - xi0;  // xj - xi
// //                     T xij1 = x[j*dim+1] - xi1; // xj - xi
// //                     // distance between atom i and atom j 
// //                     T rij = xij0*(A00*xij0 + A01*xij1) + xij1*(A10*xij0 + A11*xij1);                        
// //                     if (rij <= onet) {
// //                         //verletlist[nstart + ninc] = j; // add atom j into the list
// //                         verletlist[nstart + ninc] = (selfflag==1) ? alist[j] : j; // add atom j into the list
// //                         ninc += 1;
// //                     }
// //                 }
// //             }
// //         }               
// //     }        
// // }
// // template void cpuCreateVerletList2D(int*, double*, double*, int*, int*, int*, int*, int*, int*, int*, int, int, int);
// // template void cpuCreateVerletList2D(int*, float*, float*, int*, int*, int*, int*, int*, int*, int*, int, int, int);
// // 
// // template <typename T> void cpuFullNeighNum2D(int *neighnum, int *inside, T *x, T* ellipsoid, 
// //         int *alist, int *verletlist, int *verletnumsum, int inum, int dim)
// // {
// //     T A00 = ellipsoid[0]; // ellipsoid for pair (i,j)   
// //     T A10 = ellipsoid[1]; // ellipsoid for pair (i,j)
// //     T A01 = ellipsoid[2]; // ellipsoid for pair (i,j)  
// //     T A11 = ellipsoid[3]; // ellipsoid for pair (i,j)                                                                              
// //     for (int ii=0; ii<inum; ii++) {  // for each atom i in the simulation box     
// //         int i = alist[ii];  // atom i
// //         T xi0 = x[i*dim];        // position of atom i
// //         T xi1 = x[i*dim+1];      // position of atom i        
// //         int start = verletnumsum[i];     // starting 
// //         int m = verletnumsum[i+1]-start; // number of verlet atoms around i 
// //         int ninc = 0;                // increment 
// //         for (int l=0; l<m ; l++) {   // loop over each verlet atom around atom i
// //             int j = verletlist[start+l];  // verlet atom j  
// //             //inum + glistnumsum[j] + k
// //             T xij0 = x[j*dim] - xi0;  // xj - xi
// //             T xij1 = x[j*dim+1] - xi1; // xj - xi
// //             // distance between atom i and atom j 
// //             T rij = xij0*(A00*xij0 + A01*xij1) + xij1*(A10*xij0 + A11*xij1);                        
// //             if (rij <= 1.0) {// if atom j is inside the ellipsoid 
// //                 ninc += 1;  // increase the number of neighbors by 1               
// //                 inside[start+l] = 1;
// //             }
// //             else
// //                 inside[start+l] = 0;
// //         }
// //         neighnum[i] = ninc; // number of neighbors of atom i
// //     }
// // }
// // template void cpuFullNeighNum2D(int*, int*, double*, double*, int*, int*, int*, int, int);
// // template void cpuFullNeighNum2D(int*, int*, float*, float*, int*, int*, int*, int, int);
// // 
// // template <typename T> void cpuFullNeighList2D(int *neighlist, T *x, T* ellipsoid, int *inside, int *neighnum, 
// //         int *neighnumsum, int *alist, int *verletlist, int *verletnumsum, int inum, int dim)
// // {
// //     // get neighbor counts
// //     cpuFullNeighNum2D(neighnum, inside, x, ellipsoid, alist, verletlist, verletnumsum, inum, dim);
// //         
// //     // a list contains the starting positions of the first neighbor 
// //     cpuCumsum(neighnumsum, neighnum, inum+1); 
// //             
// //     for (int ii=0; ii<inum; ii++) {  // for each atom i in the simulation box     
// //         int i = alist[ii];  // atom i
// //         int jstart = verletnumsum[i];  // verlet starting position
// //         int nstart = neighnumsum[i];   // neighbor starting position  
// //         int m = verletnumsum[i+1]-jstart; // number of atoms around i         
// //         int ninc = 0;                // increment 
// //         for (int l=0; l<m ; l++) {   // loop over each atom around atom i            
// //             if (inside[jstart+l]) { // if atom j is inside the ellipsoid 
// //                 ninc += 1;  // increase the number of neighbors by 1               
// //                 neighlist[nstart+ninc] = verletlist[jstart+l];            
// //             }            
// //         }
// //     }    
// // }
// // template void cpuFullNeighList2D(int*, double*, double*, int*, int*, int*, int*, int*, int*, int, int);
// // template void cpuFullNeighList2D(int*, float*, float*, int*, int*, int*, int*, int*, int*, int, int);
// // 
// // template <typename T> void cpuHalfNeighNum2D(int *neighnum, int *inside, T *x, T* ellipsoid, 
// //         int *alist, int *verletlist, int *verletnumsum, int inum, int dim)
// // {    
// //     T A00 = ellipsoid[0]; // ellipsoid for pair (i,j)   
// //     T A10 = ellipsoid[1]; // ellipsoid for pair (i,j)
// //     T A01 = ellipsoid[2]; // ellipsoid for pair (i,j)  
// //     T A11 = ellipsoid[3]; // ellipsoid for pair (i,j)                                                                                                  
// //     for (int ii=0; ii<inum; ii++) {  // for each atom i in the simulation box     
// //         int i = alist[ii];  // atom i
// //         T xi0 = x[i*dim];        // position of atom i
// //         T xi1 = x[i*dim+1];      // position of atom i        
// //         int start = verletnumsum[i];     // starting 
// //         int m = verletnumsum[i+1]-start; // number of atoms around i 
// //         int ninc = 0;                // increment 
// //         for (int l=0; l<m ; l++) {   // loop over each atom around atom i
// //             int j = verletlist[start+l];  // atom j        
// //             if (j >= i) {
// //                 T xij0 = x[j*dim] - xi0;  // xj - xi
// //                 T xij1 = x[j*dim+1] - xi1; // xj - xi
// //                 // distance between atom i and atom j 
// //                 T rij = xij0*(A00*xij0 + A01*xij1) + xij1*(A10*xij0 + A11*xij1);                        
// //                 if (rij <= 1.0) {// if atom j is inside the ellipsoid 
// //                     ninc += 1;  // increase the number of neighbors by 1            
// //                     inside[start+l] = 1;
// //                 }
// //                 else
// //                     inside[start+l] = 0;
// //             }
// //             else 
// //                 inside[start+l] = 0;
// //         }
// //         neighnum[i] = ninc; // number of neighbors of atom i
// //     }
// // }
// // template void cpuHalfNeighNum2D(int*, int*, double*, double*, int*, int*, int*, int, int);
// // template void cpuHalfNeighNum2D(int*, int*, float*, float*, int*, int*, int*, int, int);
// // 
// // template <typename T> void cpuHalfNeighList2D(int *neighlist, T *x, T* ellipsoid, int *inside, int *neighnum, 
// //         int *neighnumsum, int *alist, int *verletlist, int *verletnumsum, int inum, int dim)
// // { 
// //     // get neighbor counts
// //     cpuHalfNeighNum2D(neighnum, inside, x, ellipsoid, alist, verletlist, verletnumsum, inum, dim);
// //     
// //     // a list contains the starting positions of the first neighbor 
// //     cpuCumsum(neighnumsum, neighnum, inum+1);      
// //     
// //     for (int ii=0; ii<inum; ii++) {  // for each atom i in the simulation box     
// //         int i = alist[ii];  // atom i
// //         int jstart = verletnumsum[i];  // verlet starting position
// //         int nstart = neighnumsum[i];   // neighbor starting position  
// //         int m = verletnumsum[i+1]-jstart; // number of atoms around i         
// //         int ninc = 0;                // increment 
// //         for (int l=0; l<m ; l++) {   // loop over each atom around atom i            
// //             if (inside[jstart+l]) { // if atom j is inside the ellipsoid 
// //                 ninc += 1;  // increase the number of neighbors by 1               
// //                 neighlist[nstart+ninc] = verletlist[jstart+l];            
// //             }            
// //         }
// //     }
// // }
// // template void cpuHalfNeighList2D(int*, double*, double*, int*, int*, int*, int*, int*, int*, int, int);
// // template void cpuHalfNeighList2D(int*, float*, float*, int*, int*, int*, int*, int*, int*, int, int);
// // 
// 
// //    
// // template <typename T> int cpuGetNeighPairs(T *xij, T *x, int *anum, int *anumsum, int *ai, int *aj, int *ti, 
// //       int *tj, int *ilist, int *tlist, int *neighlist, int *neighnum, int *neighnumsum, int *atomtype, 
// //        int typei, int inum, int dim)
// // {    
// //     // form ilist from olist
// //     int tnum = 0;
// //     for (int ii=0; ii<inum; ii++) {  // for each atom i in the simulation box     
// //         int i = ilist[ii];           // atom i
// //         if (atomtype[i] == typei)
// //             tlist[tnum++] = i;
// //     }
// //     
// //     cpuGetNeighPairs(xij, x, anum, anumsum, ai, aj, ti, tj, tlist, 
// //             neighlist, neighnum, neighnumsum, atomtype, tnum, dim);     
// //     
// //     return tnum;
// // }
// // template int cpuGetNeighPairs2D(double*, double*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int, int, int);
// // template int cpuGetNeighPairs2D(float*, float*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int, int, int);
// // 
// 
// // template <typename T> void cpuGetNeighPairs(T *xij, T *x, int *anum, int *anumsum, int *ai, int *aj, int *ti, 
// //       int *tj, int *ilist, int *neighlist, int *neighnum, int *neighnumsum, int *atomtype, int typej, int inum, int dim)
// // {    
// //     // form anum
// //     for (int ii=0; ii<inum; ii++) {
// //         int i = ilist[ii];       // atom i
// //         int istart = neighnumsum[i];  
// //         int m = neighnum[i];     // number of neighbors around i             
// //         anum[ii] = 0;              
// //         for (int l=0; l<m ; l++) {   // loop over each atom around atom i
// //             int j = neighlist[istart+l];  // atom j              
// //             if (atomtype[j] == typej) 
// //                 anum[ii] += 1;
// //         }                
// //     }
// //     
// //     // a list contains the starting positions of the first neighbor 
// //     cpuCumsum(anumsum, anum, inum+1);             
// //     
// //     if (dim==2) {    
// //         for (int ii=0; ii<inum; ii++) {  // for each atom i in the simulation box     
// //             int i = ilist[ii];       // atom i
// //             T xi0 = x[i*dim];        // position of atom i
// //             T xi1 = x[i*dim+1];      // position of atom i
// //             int itype = atomtype[i];
// //             int istart = neighnumsum[i];   
// //             int m = neighnum[i];     // number of neighbors around i                     
// //             int nstart = anumsum[ii];      
// //             for (int l=0; l<m ; l++) {   // loop over each atom around atom i
// //                 int j = neighlist[istart+l];  // atom j  
// //                 int jtype = atomtype[j];
// //                 if (jtype == typej) {                    
// //                     int k = nstart + anum[ii];                   
// //                     xij[k*dim]   = x[j*dim] - xi0;  // xj - xi
// //                     xij[k*dim+1] = x[j*dim+1] - xi1; // xj - xi     
// //                     ai[k]        = i;
// //                     aj[k]        = j;          
// //                     ti[k]        = itype;       
// //                     tj[k]        = jtype;    
// //                 }
// //             }
// //         }
// //     }
// //     else if (dim==3) {
// //         for (int ii=0; ii<inum; ii++) {  // for each atom i in the simulation box         
// //             int i = ilist[ii];  // atom i
// //             T xi0 = x[i*dim];        // position of atom i
// //             T xi1 = x[i*dim+1];      // position of atom i
// //             T xi2 = x[i*dim+2];      // position of atom i
// //             int itype = atomtype[i];
// //             int istart = neighnumsum[i];   
// //             int m = neighnum[i];     // number of neighbors around i             
// //             int nstart = anumsum[ii];   
// //             for (int l=0; l<m ; l++) {   // loop over each atom around atom i
// //                 int j = neighlist[istart+l];  // atom j    
// //                 int jtype = atomtype[j];
// //                 if (jtype == typej) {                    
// //                     int k = nstart + anum[ii];                 
// //                     xij[k*dim]   = x[j*dim] - xi0;  // xj - xi
// //                     xij[k*dim+1] = x[j*dim+1] - xi1; // xj - xi                                                
// //                     xij[k*dim+2] = x[j*dim+2] - xi2; // xj - xi      
// //                     ai[k]        = i;
// //                     aj[k]        = j;                
// //                     ti[k]        = itype;       
// //                     tj[k]        = jtype;    
// //                 }
// //             }            
// //         }
// //     }    
// // }
// // template void cpuGetNeighPairs(double*, double*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int, int, int);
// // template void cpuGetNeighPairs(float*, float*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int, int, int);   
// // 
// // template <typename T> int cpuGetNeighPairs(T *xij, T *x, int *anum, int *anumsum, int *ai, int *aj, int *ti, 
// //       int *tj, int *ilist, int *tlist, int *neighlist, int *neighnum, int *neighnumsum, int *atomtype, 
// //        int typei, int typej, int inum, int dim)
// // {    
// //     // form ilist from olist
// //     int tnum = 0;
// //     for (int ii=0; ii<inum; ii++) {  // for each atom i in the simulation box     
// //         int i = ilist[ii];           // atom i
// //         if (atomtype[i] == typei)
// //             tlist[tnum++] = i;
// //     }
// //     
// //     cpuGetNeighPairs(xij, x, anum, anumsum, ai, aj, ti, tj, tlist, 
// //             neighlist, neighnum, neighnumsum, atomtype, typej, tnum, dim);
// //     
// //     return tnum;
// // }
// // template int cpuGetNeighPairs(double*, double*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int, int, int, int);
// // template int cpuGetNeighPairs(float*, float*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int, int, int, int);
// // 
// // 
// // // template <typename T> void cpuGetNeighPairs(T *xij, T *x, int *ti, int *tj, int *ilist, int *neighlist,  
// // //         int *neighnum, int *neighnumsum, int *atomtype, int ntype, int inum, int dim)
// // // {    
// // //     if (dim==2) {    
// // //         for (int ii=0; ii<inum; ii++) {  // for each atom i in the simulation box     
// // //             int i = ilist[ii];  // atom i
// // //             int itype = atomtype[i];
// // //             T xi0 = x[i*dim];        // position of atom i
// // //             T xi1 = x[i*dim+1];      // position of atom i
// // //             int m = neighnum[i];     // number of neighbors around i 
// // //             int nstart = neighnumsum[i];   
// // //             for (int l=0; l<m ; l++) {   // loop over each atom around atom i
// // //                 int k = nstart+l;
// // //                 int j = neighlist[k];  // atom j     
// // //                 xij[k*dim]   = x[j*dim] - xi0;  // xj - xi
// // //                 xij[k*dim+1] = x[j*dim+1] - xi1; // xj - xi      
// // //                 ti[k]       = itype;
// // //                 tj[k]       = atomtype[j];
// // //             }
// // //         }
// // //     }
// // //     else if (dim==3) {
// // //         for (int ii=0; ii<inum; ii++) {  // for each atom i in the simulation box         
// // //             int i = ilist[ii];  // atom i
// // //             int itype = atomtype[i];
// // //             T xi0 = x[i*dim];        // position of atom i
// // //             T xi1 = x[i*dim+1];      // position of atom i
// // //             T xi2 = x[i*dim+2];      // position of atom i
// // //             int m = neighnum[i];     // number of neighbors around i 
// // //             int nstart = neighnumsum[i];               
// // //             for (int l=0; l<m ; l++) {   // loop over each atom around atom i
// // //                 int k = nstart+l;
// // //                 int j = neighlist[k];  // atom j     
// // //                 xij[k*dim]   = x[j*dim] - xi0;  // xj - xi
// // //                 xij[k*dim+1] = x[j*dim+1] - xi1; // xj - xi                                                
// // //                 xij[k*dim+2] = x[j*dim+2] - xi2; // xj - xi       
// // //                 ti[k]       = itype;
// // //                 tj[k]       = atomtype[j];
// // //             }            
// // //         }
// // //     }
// // // }
// // // template void cpuGetNeighPairs(double*, double*, int*, int*, int*, int*, int*, int*, int*, int, int, int);
// // // template void cpuGetNeighPairs(float*, float*, int*, int*, int*, int*, int*, int*, int*, int, int, int);
// // 
// // // template <typename T> void cpuGetNeighPairs(T *xi, T *xj, T *x, int *ti, int *tj, int *ilist, int *neighlist,  
// // //         int *neighnum, int *neighnumsum, int *atomtype, int ntype, int inum, int dim)
// // // {    
// // //     if (dim==2) {    
// // //         for (int ii=0; ii<inum; ii++) {  // for each atom i in the simulation box     
// // //             int i = ilist[ii];  // atom i
// // //             int itype = atomtype[i];
// // //             T xi0 = x[i*dim];        // position of atom i
// // //             T xi1 = x[i*dim+1];      // position of atom i            
// // //             int m = neighnum[i];     // number of neighbors around i 
// // //             int nstart = neighnumsum[i];   
// // //             for (int l=0; l<m ; l++) {   // loop over each atom around atom i
// // //                 int k = nstart+l;
// // //                 int j = neighlist[k];  // atom j     
// // //                 xi[k*dim]   = xi0;        // position of atom i
// // //                 xi[k*dim+1] = xi1;       // position of atom i                
// // //                 xj[k*dim]   = x[j*dim];  // xj 
// // //                 xj[k*dim+1] = x[j*dim+1]; // xj    
// // //                 ti[k]       = itype;
// // //                 tj[k]       = atomtype[j];
// // //             }
// // //         }
// // //     }
// // //     else if (dim==3) {
// // //         for (int ii=0; ii<inum; ii++) {  // for each atom i in the simulation box         
// // //             int i = ilist[ii];  // atom i
// // //             int itype = atomtype[i];
// // //             T xi0 = x[i*dim];        // position of atom i
// // //             T xi1 = x[i*dim+1];      // position of atom i
// // //             T xi2 = x[i*dim+2];      // position of atom i
// // //             int m = neighnum[i];     // number of neighbors around i 
// // //             int nstart = neighnumsum[i];               
// // //             for (int l=0; l<m ; l++) {   // loop over each atom around atom i
// // //                 int k = nstart+l;
// // //                 int j = neighlist[k];  // atom j     
// // //                 xi[k*dim]   = xi0;        // position of atom i
// // //                 xi[k*dim+1] = xi1;       // position of atom i                
// // //                 xi[k*dim+2] = xi2;       // position of atom i                
// // //                 xj[k*dim]   = x[j*dim];  // xj - xi
// // //                 xj[k*dim+1] = x[j*dim+1]; // xj - xi                                                
// // //                 xj[k*dim+2] = x[j*dim+2]; // xj - xi       
// // //                 ti[k]       = itype;
// // //                 tj[k]       = atomtype[j];
// // //             }            
// // //         }
// // //     }
// // // }
// // // template void cpuGetNeighPairs(double*, double*, double*, int*, int*, int*, int*, int*, int*, int*, int, int, int);
// // // template void cpuGetNeighPairs(float*, float*, float*, int*, int*, int*, int*, int*, int*, int*, int, int, int);
// // 
// // 
// // // ***************************************************************/
// // template <typename T> void cpuVerletAtoms(int *verletnum, T *x, T *ellipsoid, int *atomtype, int *ilist, 
// //         int *clist, int *c2alist, int *c2anum, int *c2anumsum, int *nc, int ntype, int inum, int dim)
// // {
// //     int dimsq = dim*dim;
// //     T onet = (T) 1.0;
// //     if (dim==2) {
// //         for (int ii=0; ii<inum; ii++) {  // for each atom i in the simulation box    
// //             int i = ilist[ii];
// //             verletnum[i] = 0;
// //             T xi0 = x[i*dim];        // position of atom i
// //             T xi1 = x[i*dim+1];      // position of atom i
// //             int ti = atomtype[i];      // element of atom i 
// //             int j = clist[i];         // cell j contains atom i           
// //             int j1 = j%nc[0];
// //             int j2 = (j-j1)/nc[0];
// //             for (int i1=-1; i1<=1; i1++) {
// //                 int k1 = j1 + i1;
// //                 for (int i2=-1; i2<=1; i2++) {
// //                     int k2 = j2 + i2;
// //                     int k = k1 + nc[0]*k2; // cell k
// //                     int m = c2anum[k];     // number of atoms in cell k
// //                     int s = c2anumsum[k];  // starting position of the first atom in cell k
// //                     for (int l=0; l<m ; l++) {
// //                         j = c2alist[s+l];  // atom j
// //                         int tj = atomtype[j];      // element of atom j
// //                         int tij = dimsq*(ti + tj*ntype);
// //                         T A00 = ellipsoid[tij]; // ellipsoid for pair (i,j)   
// //                         T A10 = ellipsoid[tij+1]; // ellipsoid for pair (i,j)
// //                         T A01 = ellipsoid[tij+2]; // ellipsoid for pair (i,j)  
// //                         T A11 = ellipsoid[tij+3]; // ellipsoid for pair (i,j)                          
// //                         T xij0 = x[j*dim] - xi0;  // xj - xi
// //                         T xij1 = x[j*dim+1] - xi1; // xj - xi
// //                         // distance between atom i and atom j 
// //                         T rij = xij0*(A00*xij0 + A01*xij1) + xij1*(A10*xij0 + A11*xij1);                        
// //                         if (rij <= onet)
// //                             verletnum[i] += 1;
// //                     }
// //                 }
// //             }                
// //         }        
// //     }
// //     else if (dim==3) {
// //         for (int ii=0; ii<inum; ii++) {  // for each atom i in the simulation box    
// //             int i = ilist[ii];
// //             verletnum[i] = 0;
// //             T xi0 = x[i*dim];        // position of atom i
// //             T xi1 = x[i*dim+1];      // position of atom i
// //             T xi2 = x[i*dim+2];      // position of atom i
// //             int ti = atomtype[i];      // element of atom i 
// //             int j = clist[i];         // cell j contains atom i       
// //             int n = j%(nc[0]*nc[1]);
// //             int j3 = (j-n)/(nc[0]*nc[1]);            
// //             int j1 = n%nc[0];
// //             int j2 = (j-j1)/nc[0];
// //             for (int i1=-1; i1<=1; i1++) {
// //                 int k1 = j1 + i1;
// //                 for (int i2=-1; i2<=1; i2++) {
// //                     int k2 = j2 + i2;
// //                     for (int i3=-1; i3<=1; i3++) {
// //                         int k3 = j3 + i3;                    
// //                         int k = k1 + nc[0]*k2 + nc[0]*nc[1]*k3; // cell k
// //                         int m = c2anum[k];     // number of atoms in cell k
// //                         int s = c2anumsum[k];  // starting position of the first atom in cell k
// //                         for (int l=0; l<m ; l++) {
// //                             j = c2alist[s+l];  // atom j
// //                             int tj = atomtype[j];      // element of atom j
// //                             int tij = dimsq*(ti + tj*ntype);
// //                             T A00 = ellipsoid[tij]; // ellipsoid for element t   
// //                             T A10 = ellipsoid[tij+1]; // ellipsoid for element t   
// //                             T A20 = ellipsoid[tij+2]; // ellipsoid for element t   
// //                             T A01 = ellipsoid[tij+3]; // ellipsoid for element t   
// //                             T A11 = ellipsoid[tij+4]; // ellipsoid for element t   
// //                             T A21 = ellipsoid[tij+5]; // ellipsoid for element t   
// //                             T A02 = ellipsoid[tij+6]; // ellipsoid for element t   
// //                             T A12 = ellipsoid[tij+7]; // ellipsoid for element t   
// //                             T A22 = ellipsoid[tij+8]; // ellipsoid for element t                               
// //                             T xij0 = x[j*dim] - xi0;  // xj - xi
// //                             T xij1 = x[j*dim+1] - xi1; // xj - xi
// //                             T xij2 = x[j*dim+2] - xi2; // xj - xi
// //                             // distance between atom i and atom j 
// //                             T rij = xij0*(A00*xij0 + A01*xij1 + A02*xij2) + xij1*(A10*xij0 + A11*xij1 + A12*xij2) + xij2*(A20*xij0 + A21*xij1 + A22*xij2);
// //                             if (rij <= onet)
// //                                 verletnum[i] += 1;
// //                         }
// //                     }
// //                 }
// //             }                
// //         }                
// //     }        
// // }
// // template void cpuVerletAtoms(int*, double*, double*, int*, int*, int*, int*, int*, int*, int*, int, int, int);
// // template void cpuVerletAtoms(int*, float*, float*, int*, int*, int*, int*, int*, int*, int*, int, int, int);
// // 
// // template <typename T> void cpuCreateVerletList(int *verletlist, T *x, T *ellipsoid, int *verletnum,  int *verletnumsum, 
// //      int *atomtype, int *ilist, int *clist, int *c2alist, int *c2anum, int *c2anumsum, int *nc, int ntype, int inum, int dim)
// // {
// //     // a list contains the starting positions of the first j atom 
// //     cpuCumsum(verletnumsum, verletnum, inum+1); 
// //     
// //     int dimsq = dim*dim;
// //     T onet = (T) 1.0;
// //     if (dim==2) {
// //         for (int ii=0; ii<inum; ii++) {  // for each atom i in the simulation box    
// //             int i = ilist[ii];
// //             T xi0 = x[i*dim];        // position of atom i
// //             T xi1 = x[i*dim+1];      // position of atom i
// //             int ti = atomtype[i];      // element of atom i             
// //             int nstart = verletnumsum[i];     // starting 
// //             int ninc = 0;                // increment
// //             int j = clist[i];         // cell j contains atom i           
// //             int j1 = j%nc[0];
// //             int j2 = (j-j1)/nc[0];
// //             for (int i1=-1; i1<=1; i1++) {
// //                 int k1 = j1 + i1;
// //                 for (int i2=-1; i2<=1; i2++) {
// //                     int k2 = j2 + i2;
// //                     int k = k1 + nc[0]*k2; // cell k
// //                     int m = c2anum[k];     // number of atoms in cell k
// //                     int s = c2anumsum[k];  // starting position of the first atom in cell k
// //                     for (int l=0; l<m ; l++) { //loop over each atom j in cell k
// //                         j = c2alist[s+l];  // atom j
// //                         int tj = atomtype[j];      // element of atom j
// //                         int tij = dimsq*(ti + tj*ntype);
// //                         T A00 = ellipsoid[tij]; // ellipsoid for pair (i,j)   
// //                         T A10 = ellipsoid[tij+1]; // ellipsoid for pair (i,j)
// //                         T A01 = ellipsoid[tij+2]; // ellipsoid for pair (i,j)  
// //                         T A11 = ellipsoid[tij+3]; // ellipsoid for pair (i,j)                                                  
// //                         T xij0 = x[j*dim] - xi0;  // xj - xi
// //                         T xij1 = x[j*dim+1] - xi1; // xj - xi
// //                         // distance between atom i and atom j 
// //                         T rij = xij0*(A00*xij0 + A01*xij1) + xij1*(A10*xij0 + A11*xij1);                        
// //                         if (rij <= onet) {
// //                             verletlist[nstart + ninc] = j; // add atom j into the list
// //                             ninc += 1;
// //                         }
// //                     }
// //                 }
// //             }               
// //             verletnum[i] = ninc;
// //         }        
// //     }
// //     else if (dim==3) {
// //         for (int ii=0; ii<inum; ii++) {  // for each atom i in the simulation box    
// //             int i = ilist[ii];
// //             T xi0 = x[i*dim];        // position of atom i
// //             T xi1 = x[i*dim+1];      // position of atom i
// //             T xi2 = x[i*dim+2];      // position of atom i
// //             int ti = atomtype[i];      // element of atom i  
// //             int nstart = verletnumsum[i];     // starting 
// //             int ninc = 0;                // increment            
// //             int j = clist[i];         // cell j contains atom i       
// //             int n = j%(nc[0]*nc[1]);
// //             int j3 = (j-n)/(nc[0]*nc[1]);            
// //             int j1 = n%nc[0];
// //             int j2 = (j-j1)/nc[0];
// //             for (int i1=-1; i1<=1; i1++) {
// //                 int k1 = j1 + i1;
// //                 for (int i2=-1; i2<=1; i2++) {
// //                     int k2 = j2 + i2;
// //                     for (int i3=-1; i3<=1; i3++) {
// //                         int k3 = j3 + i3;                    
// //                         int k = k1 + nc[0]*k2 + nc[0]*nc[1]*k3; // cell k
// //                         int m = c2anum[k];     // number of atoms in cell k
// //                         int s = c2anumsum[k];  // starting position of the first atom in cell k
// //                         for (int l=0; l<m ; l++) {
// //                             j = c2alist[s+l];  // atom j
// //                             int tj = atomtype[j];      // element of atom j
// //                             int tij = dimsq*(ti + tj*ntype);
// //                             T A00 = ellipsoid[tij]; // ellipsoid for element t   
// //                             T A10 = ellipsoid[tij+1]; // ellipsoid for element t   
// //                             T A20 = ellipsoid[tij+2]; // ellipsoid for element t   
// //                             T A01 = ellipsoid[tij+3]; // ellipsoid for element t   
// //                             T A11 = ellipsoid[tij+4]; // ellipsoid for element t   
// //                             T A21 = ellipsoid[tij+5]; // ellipsoid for element t   
// //                             T A02 = ellipsoid[tij+6]; // ellipsoid for element t   
// //                             T A12 = ellipsoid[tij+7]; // ellipsoid for element t   
// //                             T A22 = ellipsoid[tij+8]; // ellipsoid for element t                                                           
// //                             T xij0 = x[j*dim] - xi0;  // xj - xi
// //                             T xij1 = x[j*dim+1] - xi1; // xj - xi
// //                             T xij2 = x[j*dim+2] - xi2; // xj - xi
// //                             // distance between atom i and atom j 
// //                             T rij = xij0*(A00*xij0 + A01*xij1 + A02*xij2) + xij1*(A10*xij0 + A11*xij1 + A12*xij2) + xij2*(A20*xij0 + A21*xij1 + A22*xij2);
// //                             if (rij <= onet) {
// //                                 verletlist[nstart + ninc] = j; // add atom j into the list
// //                                 ninc += 1;
// //                             }                            
// //                         }
// //                     }
// //                 }
// //             }                
// //             verletnum[i] = ninc;
// //         }        
// //     }        
// // }
// // template void cpuCreateVerletList(int*, double*, double*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int, int, int);
// // template void cpuCreateVerletList(int*, float*, float*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int, int, int);
// // 
// // template <typename T> void cpuFullNeighNum(int *neighnum, T *x, T* ellipsoid, int *atomtype, 
// //         int *ilist, int *verletlist, int *verletnum, int *verletnumsum, int ntype, int inum, int dim)
// // {
// //     int dimsq = dim*dim;
// //     
// //     if (dim==2) {    
// //         for (int ii=0; ii<inum; ii++) {  // for each atom i in the simulation box     
// //             int i = ilist[ii];  // atom i
// //             T xi0 = x[i*dim];        // position of atom i
// //             T xi1 = x[i*dim+1];      // position of atom i
// //             int ti = atomtype[i];      // element of atom i           
// //             int m = verletnum[i]; // number of atoms around i 
// //             int start = verletnumsum[i];     // starting 
// //             int ninc = 0;                // increment 
// //             for (int l=0; l<m ; l++) {   // loop over each atom around atom i
// //                 int j = verletlist[start+l];  // atom j              
// //                 int tj = atomtype[j];      // element of atom j
// //                 int tij = dimsq*(ti + tj*ntype);
// //                 T A00 = ellipsoid[tij]; // ellipsoid for pair (i,j)   
// //                 T A10 = ellipsoid[tij+1]; // ellipsoid for pair (i,j)
// //                 T A01 = ellipsoid[tij+2]; // ellipsoid for pair (i,j)  
// //                 T A11 = ellipsoid[tij+3]; // ellipsoid for pair (i,j)                                                                  
// //                 T xij0 = x[j*dim] - xi0;  // xj - xi
// //                 T xij1 = x[j*dim+1] - xi1; // xj - xi
// //                 // distance between atom i and atom j 
// //                 T rij = xij0*(A00*xij0 + A01*xij1) + xij1*(A10*xij0 + A11*xij1);                        
// //                 if (rij <= 1.0) // if atom j is inside the ellipsoid 
// //                     ninc += 1;  // increase the number of neighbors by 1               
// //             }
// //             neighnum[i] = ninc; // number of neighbors of atom i
// //         }
// //     }
// //     else if (dim==3) {
// //         for (int ii=0; ii<inum; ii++) {  // for each atom i in the simulation box         
// //             int i = ilist[ii];  // atom i
// //             T xi0 = x[i*dim];        // position of atom i
// //             T xi1 = x[i*dim+1];      // position of atom i
// //             T xi2 = x[i*dim+2];      // position of atom i
// //             int ti = atomtype[i];      // element of atom i           
// //             int m = verletnum[i]; // number of atoms around i 
// //             int start = verletnumsum[i];     // starting 
// //             int ninc = 0;                // increment 
// //             for (int l=0; l<m ; l++) {   // loop over each atom around atom i
// //                 int j = verletlist[start+l];  // atom j                
// //                 int tj = atomtype[j];      // element of atom j
// //                 int tij = dimsq*(ti + tj*ntype);
// //                 T A00 = ellipsoid[tij]; // ellipsoid for element t   
// //                 T A10 = ellipsoid[tij+1]; // ellipsoid for element t   
// //                 T A20 = ellipsoid[tij+2]; // ellipsoid for element t   
// //                 T A01 = ellipsoid[tij+3]; // ellipsoid for element t   
// //                 T A11 = ellipsoid[tij+4]; // ellipsoid for element t   
// //                 T A21 = ellipsoid[tij+5]; // ellipsoid for element t   
// //                 T A02 = ellipsoid[tij+6]; // ellipsoid for element t   
// //                 T A12 = ellipsoid[tij+7]; // ellipsoid for element t   
// //                 T A22 = ellipsoid[tij+8]; // ellipsoid for element t                                                                           
// //                 T xij0 = x[j*dim] - xi0;  // xj - xi
// //                 T xij1 = x[j*dim+1] - xi1; // xj - xi
// //                 T xij2 = x[j*dim+2] - xi2; // xj - xi
// //                 // distance between atom i and atom j 
// //                 T rij = xij0*(A00*xij0 + A01*xij1 + A02*xij2) + xij1*(A10*xij0 + A11*xij1 + A12*xij2) + xij2*(A20*xij0 + A21*xij1 + A22*xij2);
// //                 if (rij <= 1.0) // if atom j is inside the ellipsoid 
// //                     ninc += 1;  // increase the number of neighbors by 1                               
// //             }            
// //             neighnum[i] = ninc; // number of neighbors of atom i
// //         }
// //     }
// // }
// // template void cpuFullNeighNum(int*, double*, double*, int*, int*, int*, int*, int*, int, int, int);
// // template void cpuFullNeighNum(int*, float*, float*, int*, int*, int*, int*, int*, int, int, int);
// // 
// // template <typename T> void cpuFullNeighList(int *neighlist, T *x, T* ellipsoid, int *neighnum, int *neighnumsum,
// //         int *atomtype, int *ilist, int *verletlist, int *verletnum, int *verletnumsum, int ntype, int inum, int dim)
// // {
// //     // a list contains the starting positions of the first neighbor 
// //     cpuCumsum(neighnumsum, neighnum, inum+1); 
// //     
// //     int dimsq = dim*dim;
// //     
// //     if (dim==2) {    
// //         for (int ii=0; ii<inum; ii++) {  // for each atom i in the simulation box     
// //             int i = ilist[ii];  // atom i
// //             T xi0 = x[i*dim];        // position of atom i
// //             T xi1 = x[i*dim+1];      // position of atom i
// //             int ti = atomtype[i];      // element of atom i
// //             int m = verletnum[i]; // number of atoms around i 
// //             int jstart = verletnumsum[i];     // starting 
// //             int nstart = neighnumsum[i];   
// //             int ninc = 0;                // increment 
// //             for (int l=0; l<m ; l++) {   // loop over each atom around atom i
// //                 int j = verletlist[jstart+l];  // atom j                
// //                 int tj = atomtype[j];      // element of atom j
// //                 int tij = dimsq*(ti + tj*ntype);
// //                 T A00 = ellipsoid[tij]; // ellipsoid for pair (i,j)   
// //                 T A10 = ellipsoid[tij+1]; // ellipsoid for pair (i,j)
// //                 T A01 = ellipsoid[tij+2]; // ellipsoid for pair (i,j)  
// //                 T A11 = ellipsoid[tij+3]; // ellipsoid for pair (i,j)                                                                                  
// //                 T xij0 = x[j*dim] - xi0;  // xj - xi
// //                 T xij1 = x[j*dim+1] - xi1; // xj - xi
// //                 // distance between atom i and atom j 
// //                 T rij = xij0*(A00*xij0 + A01*xij1) + xij1*(A10*xij0 + A11*xij1);                        
// //                 if (rij <= 1.0) { // if atom j is inside the ellipsoid 
// //                     neighlist[nstart+ninc] = j;
// //                     ninc += 1;  // increase the number of neighbors by 1               
// //                 }
// //             }
// //             neighnum[i] = ninc; // number of neighbors of atom i
// //         }
// //     }
// //     else if (dim==3) {
// //         for (int ii=0; ii<inum; ii++) {  // for each atom i in the simulation box         
// //             int i = ilist[ii];  // atom i
// //             T xi0 = x[i*dim];        // position of atom i
// //             T xi1 = x[i*dim+1];      // position of atom i
// //             T xi2 = x[i*dim+2];      // position of atom i
// //             int ti = atomtype[i];      // element of atom i
// //             int m = verletnum[i]; // number of atoms around i 
// //             int jstart = verletnumsum[i];     // starting 
// //             int nstart = neighnumsum[i];   
// //             int ninc = 0;                // increment 
// //             for (int l=0; l<m ; l++) {   // loop over each atom around atom i
// //                 int j = verletlist[jstart+l];  // atom j              
// //                 int tj = atomtype[j];      // element of atom j
// //                 int tij = dimsq*(ti + tj*ntype);
// //                 T A00 = ellipsoid[tij]; // ellipsoid for element t   
// //                 T A10 = ellipsoid[tij+1]; // ellipsoid for element t   
// //                 T A20 = ellipsoid[tij+2]; // ellipsoid for element t   
// //                 T A01 = ellipsoid[tij+3]; // ellipsoid for element t   
// //                 T A11 = ellipsoid[tij+4]; // ellipsoid for element t   
// //                 T A21 = ellipsoid[tij+5]; // ellipsoid for element t   
// //                 T A02 = ellipsoid[tij+6]; // ellipsoid for element t   
// //                 T A12 = ellipsoid[tij+7]; // ellipsoid for element t   
// //                 T A22 = ellipsoid[tij+8]; // ellipsoid for element t                                                                                           
// //                 T xij0 = x[j*dim] - xi0;  // xj - xi
// //                 T xij1 = x[j*dim+1] - xi1; // xj - xi
// //                 T xij2 = x[j*dim+2] - xi2; // xj - xi
// //                 // distance between atom i and atom j 
// //                 T rij = xij0*(A00*xij0 + A01*xij1 + A02*xij2) + xij1*(A10*xij0 + A11*xij1 + A12*xij2) + xij2*(A20*xij0 + A21*xij1 + A22*xij2);
// //                 if (rij <= 1.0) { // if atom j is inside the ellipsoid 
// //                     neighlist[nstart+ninc] = j;
// //                     ninc += 1;  // increase the number of neighbors by 1               
// //                 }
// //             }            
// //             neighnum[i] = ninc; // number of neighbors of atom i
// //         }
// //     }
// // }
// // template void cpuFullNeighList(int*, double*, double*, int*, int*, int*, int*, int*, int*, int*, int, int, int);
// // template void cpuFullNeighList(int*, float*, float*, int*, int*, int*, int*, int*, int*, int*, int, int, int);
// // 
// // template <typename T> void cpuHalfNeighNum(int *neighnum, T *x, T* ellipsoid, int *atomtype, 
// //         int *ilist, int *verletlist, int *verletnum, int *verletnumsum, int ntype, int inum, int dim)
// // {
// //     int dimsq = dim*dim;
// //     
// //     if (dim==2) {    
// //         for (int ii=0; ii<inum; ii++) {  // for each atom i in the simulation box     
// //             int i = ilist[ii];  // atom i
// //             T xi0 = x[i*dim];        // position of atom i
// //             T xi1 = x[i*dim+1];      // position of atom i
// //             int ti = atomtype[i];      // element of atom i           
// //             int m = verletnum[i]; // number of atoms around i 
// //             int start = verletnumsum[i];     // starting 
// //             int ninc = 0;                // increment 
// //             for (int l=0; l<m ; l++) {   // loop over each atom around atom i
// //                 int j = verletlist[start+l];  // atom j        
// //                 if (j >= i) {
// //                     int tj = atomtype[j];      // element of atom j
// //                     int tij = dimsq*(ti + tj*ntype);
// //                     T A00 = ellipsoid[tij]; // ellipsoid for pair (i,j)   
// //                     T A10 = ellipsoid[tij+1]; // ellipsoid for pair (i,j)
// //                     T A01 = ellipsoid[tij+2]; // ellipsoid for pair (i,j)  
// //                     T A11 = ellipsoid[tij+3]; // ellipsoid for pair (i,j)                                                                  
// //                     T xij0 = x[j*dim] - xi0;  // xj - xi
// //                     T xij1 = x[j*dim+1] - xi1; // xj - xi
// //                     // distance between atom i and atom j 
// //                     T rij = xij0*(A00*xij0 + A01*xij1) + xij1*(A10*xij0 + A11*xij1);                        
// //                     if (rij <= 1.0) // if atom j is inside the ellipsoid 
// //                         ninc += 1;  // increase the number of neighbors by 1            
// //                 }
// //             }
// //             neighnum[i] = ninc; // number of neighbors of atom i
// //         }
// //     }
// //     else if (dim==3) {
// //         for (int ii=0; ii<inum; ii++) {  // for each atom i in the simulation box         
// //             int i = ilist[ii];  // atom i
// //             T xi0 = x[i*dim];        // position of atom i
// //             T xi1 = x[i*dim+1];      // position of atom i
// //             T xi2 = x[i*dim+2];      // position of atom i
// //             int ti = atomtype[i];      // element of atom i           
// //             int m = verletnum[i]; // number of atoms around i 
// //             int start = verletnumsum[i];     // starting 
// //             int ninc = 0;                // increment 
// //             for (int l=0; l<m ; l++) {   // loop over each atom around atom i
// //                 int j = verletlist[start+l];  // atom j                
// //                 if (j >= i) {
// //                     int tj = atomtype[j];      // element of atom j
// //                     int tij = dimsq*(ti + tj*ntype);
// //                     T A00 = ellipsoid[tij]; // ellipsoid for element t   
// //                     T A10 = ellipsoid[tij+1]; // ellipsoid for element t   
// //                     T A20 = ellipsoid[tij+2]; // ellipsoid for element t   
// //                     T A01 = ellipsoid[tij+3]; // ellipsoid for element t   
// //                     T A11 = ellipsoid[tij+4]; // ellipsoid for element t   
// //                     T A21 = ellipsoid[tij+5]; // ellipsoid for element t   
// //                     T A02 = ellipsoid[tij+6]; // ellipsoid for element t   
// //                     T A12 = ellipsoid[tij+7]; // ellipsoid for element t   
// //                     T A22 = ellipsoid[tij+8]; // ellipsoid for element t                                                                           
// //                     T xij0 = x[j*dim] - xi0;  // xj - xi
// //                     T xij1 = x[j*dim+1] - xi1; // xj - xi
// //                     T xij2 = x[j*dim+2] - xi2; // xj - xi
// //                     // distance between atom i and atom j 
// //                     T rij = xij0*(A00*xij0 + A01*xij1 + A02*xij2) + xij1*(A10*xij0 + A11*xij1 + A12*xij2) + xij2*(A20*xij0 + A21*xij1 + A22*xij2);
// //                     if (rij <= 1.0) // if atom j is inside the ellipsoid 
// //                         ninc += 1;  // increase the number of neighbors by 1        
// //                 }
// //             }            
// //             neighnum[i] = ninc; // number of neighbors of atom i
// //         }
// //     }
// // }
// // template void cpuHalfNeighNum(int*, double*, double*, int*, int*, int*, int*, int*, int, int, int);
// // template void cpuHalfNeighNum(int*, float*, float*, int*, int*, int*, int*, int*, int, int, int);
// // 
// // template <typename T> void cpuHalfNeighList(int *neighlist, T *x, T* ellipsoid, int *neighnum, int *neighnumsum,
// //         int *atomtype, int *ilist, int *verletlist, int *verletnum, int *verletnumsum, int ntype, int inum, int dim)
// // {
// //     // a list contains the starting positions of the first neighbor 
// //     cpuCumsum(neighnumsum, neighnum, inum+1); 
// //     
// //     int dimsq = dim*dim;
// //     
// //     if (dim==2) {    
// //         for (int ii=0; ii<inum; ii++) {  // for each atom i in the simulation box     
// //             int i = ilist[ii];  // atom i
// //             T xi0 = x[i*dim];        // position of atom i
// //             T xi1 = x[i*dim+1];      // position of atom i
// //             int ti = atomtype[i];      // element of atom i
// //             int m = verletnum[i]; // number of atoms around i 
// //             int jstart = verletnumsum[i];     // starting 
// //             int nstart = neighnumsum[i];   
// //             int ninc = 0;                // increment 
// //             for (int l=0; l<m ; l++) {   // loop over each atom around atom i
// //                 int j = verletlist[jstart+l];  // atom j     
// //                 if (j >= i) {
// //                     int tj = atomtype[j];      // element of atom j
// //                     int tij = dimsq*(ti + tj*ntype);
// //                     T A00 = ellipsoid[tij]; // ellipsoid for pair (i,j)   
// //                     T A10 = ellipsoid[tij+1]; // ellipsoid for pair (i,j)
// //                     T A01 = ellipsoid[tij+2]; // ellipsoid for pair (i,j)  
// //                     T A11 = ellipsoid[tij+3]; // ellipsoid for pair (i,j)                                                                                  
// //                     T xij0 = x[j*dim] - xi0;  // xj - xi
// //                     T xij1 = x[j*dim+1] - xi1; // xj - xi
// //                     // distance between atom i and atom j 
// //                     T rij = xij0*(A00*xij0 + A01*xij1) + xij1*(A10*xij0 + A11*xij1);                        
// //                     if (rij <= 1.0) { // if atom j is inside the ellipsoid 
// //                         neighlist[nstart+ninc] = j;
// //                         ninc += 1;  // increase the number of neighbors by 1               
// //                     }
// //                 }
// //             }
// //             neighnum[i] = ninc; // number of neighbors of atom i
// //         }
// //     }
// //     else if (dim==3) {
// //         for (int ii=0; ii<inum; ii++) {  // for each atom i in the simulation box         
// //             int i = ilist[ii];  // atom i
// //             T xi0 = x[i*dim];        // position of atom i
// //             T xi1 = x[i*dim+1];      // position of atom i
// //             T xi2 = x[i*dim+2];      // position of atom i
// //             int ti = atomtype[i];      // element of atom i
// //             int m = verletnum[i]; // number of atoms around i 
// //             int jstart = verletnumsum[i];     // starting 
// //             int nstart = neighnumsum[i];   
// //             int ninc = 0;                // increment 
// //             for (int l=0; l<m ; l++) {   // loop over each atom around atom i
// //                 int j = verletlist[jstart+l];  // atom j           
// //                 if (j >= i) {
// //                     int tj = atomtype[j];      // element of atom j
// //                     int tij = dimsq*(ti + tj*ntype);
// //                     T A00 = ellipsoid[tij]; // ellipsoid for element t   
// //                     T A10 = ellipsoid[tij+1]; // ellipsoid for element t   
// //                     T A20 = ellipsoid[tij+2]; // ellipsoid for element t   
// //                     T A01 = ellipsoid[tij+3]; // ellipsoid for element t   
// //                     T A11 = ellipsoid[tij+4]; // ellipsoid for element t   
// //                     T A21 = ellipsoid[tij+5]; // ellipsoid for element t   
// //                     T A02 = ellipsoid[tij+6]; // ellipsoid for element t   
// //                     T A12 = ellipsoid[tij+7]; // ellipsoid for element t   
// //                     T A22 = ellipsoid[tij+8]; // ellipsoid for element t                                                                                           
// //                     T xij0 = x[j*dim] - xi0;  // xj - xi
// //                     T xij1 = x[j*dim+1] - xi1; // xj - xi
// //                     T xij2 = x[j*dim+2] - xi2; // xj - xi
// //                     // distance between atom i and atom j 
// //                     T rij = xij0*(A00*xij0 + A01*xij1 + A02*xij2) + xij1*(A10*xij0 + A11*xij1 + A12*xij2) + xij2*(A20*xij0 + A21*xij1 + A22*xij2);
// //                     if (rij <= 1.0) { // if atom j is inside the ellipsoid 
// //                         neighlist[nstart+ninc] = j;
// //                         ninc += 1;  // increase the number of neighbors by 1               
// //                     }
// //                 }
// //             }            
// //             neighnum[i] = ninc; // number of neighbors of atom i
// //         }
// //     }
// // }
// // template void cpuHalfNeighList(int*, double*, double*, int*, int*, int*, int*, int*, int*, int*, int, int, int);
// // template void cpuHalfNeighList(int*, float*, float*, int*, int*, int*, int*, int*, int*, int*, int, int, int);
// // 
// 
// // template <typename T> void cpuVerletAtoms(int *verletnum, T *x, T *ellipsoid, int *atomtype, 
// //         int *ilist, int *clist, int *c2alist, int *c2anum, int *c2anumsum, int *nc, int inum, int dim)
// // {
// //     int dimsq = dim*dim;
// //     T onet = (T) 1.0;
// //     if (dim==2) {
// //         for (int ii=0; ii<inum; ii++) {  // for each atom i in the simulation box    
// //             int i = ilist[ii];
// //             verletnum[i] = 0;
// //             T xi0 = x[i*dim];        // position of atom i
// //             T xi1 = x[i*dim+1];      // position of atom i
// //             int t = atomtype[i];      // element of atom i 
// //             T A00 = ellipsoid[dimsq*t]; // ellipsoid for element t   
// //             T A10 = ellipsoid[dimsq*t+1]; // ellipsoid for element t   
// //             T A01 = ellipsoid[dimsq*t+2]; // ellipsoid for element t   
// //             T A11 = ellipsoid[dimsq*t+3]; // ellipsoid for element t   
// //             int j = clist[i];         // cell j contains atom i           
// //             int j1 = j%nc[0];
// //             int j2 = (j-j1)/nc[0];
// //             for (int i1=-1; i1<=1; i1++) {
// //                 int k1 = j1 + i1;
// //                 for (int i2=-1; i2<=1; i2++) {
// //                     int k2 = j2 + i2;
// //                     int k = k1 + nc[0]*k2; // cell k
// //                     int m = c2anum[k];     // number of atoms in cell k
// //                     int s = c2anumsum[k];  // starting position of the first atom in cell k
// //                     for (int l=0; l<m ; l++) {
// //                         j = c2alist[s+l];  // atom j
// //                         T xij0 = x[j*dim] - xi0;  // xj - xi
// //                         T xij1 = x[j*dim+1] - xi1; // xj - xi
// //                         // distance between atom i and atom j 
// //                         T rij = xij0*(A00*xij0 + A01*xij1) + xij1*(A10*xij0 + A11*xij1);                        
// //                         if (rij <= onet)
// //                             verletnum[i] += 1;
// //                     }
// //                 }
// //             }                
// //         }        
// //     }
// //     else if (dim==3) {
// //         for (int ii=0; ii<inum; ii++) {  // for each atom i in the simulation box    
// //             int i = ilist[ii];
// //             T xi0 = x[i*dim];        // position of atom i
// //             T xi1 = x[i*dim+1];      // position of atom i
// //             T xi2 = x[i*dim+2];      // position of atom i
// //             int t = atomtype[i];      // element of atom i 
// //             T A00 = ellipsoid[dimsq*t]; // ellipsoid for element t   
// //             T A10 = ellipsoid[dimsq*t+1]; // ellipsoid for element t   
// //             T A20 = ellipsoid[dimsq*t+2]; // ellipsoid for element t   
// //             T A01 = ellipsoid[dimsq*t+3]; // ellipsoid for element t   
// //             T A11 = ellipsoid[dimsq*t+4]; // ellipsoid for element t   
// //             T A21 = ellipsoid[dimsq*t+5]; // ellipsoid for element t   
// //             T A02 = ellipsoid[dimsq*t+6]; // ellipsoid for element t   
// //             T A12 = ellipsoid[dimsq*t+7]; // ellipsoid for element t   
// //             T A22 = ellipsoid[dimsq*t+8]; // ellipsoid for element t   
// //             int j = clist[i];         // cell j contains atom i       
// //             int n = j%(nc[0]*nc[1]);
// //             int j3 = (j-n)/(nc[0]*nc[1]);            
// //             int j1 = n%nc[0];
// //             int j2 = (j-j1)/nc[0];
// //             for (int i1=-1; i1<=1; i1++) {
// //                 int k1 = j1 + i1;
// //                 for (int i2=-1; i2<=1; i2++) {
// //                     int k2 = j2 + i2;
// //                     for (int i3=-1; i3<=1; i3++) {
// //                         int k3 = j3 + i3;                    
// //                         int k = k1 + nc[0]*k2 + nc[0]*nc[1]*k3; // cell k
// //                         int m = c2anum[k];     // number of atoms in cell k
// //                         int s = c2anumsum[k];  // starting position of the first atom in cell k
// //                         for (int l=0; l<m ; l++) {
// //                             j = c2alist[s+l];  // atom j
// //                             T xij0 = x[j*dim] - xi0;  // xj - xi
// //                             T xij1 = x[j*dim+1] - xi1; // xj - xi
// //                             T xij2 = x[j*dim+2] - xi2; // xj - xi
// //                             // distance between atom i and atom j 
// //                             T rij = xij0*(A00*xij0 + A01*xij1 + A02*xij2) + xij1*(A10*xij0 + A11*xij1 + A12*xij2) + xij2*(A20*xij0 + A21*xij1 + A22*xij2);
// //                             if (rij <= onet)
// //                                 verletnum[i] += 1;
// //                         }
// //                     }
// //                 }
// //             }                
// //         }                
// //     }        
// // }
// // template void cpuVerletAtoms(int*, double*, double*, int*, int*, int*, int*, int*, int*, int*, int, int);
// // template void cpuVerletAtoms(int*, float*, float*, int*, int*, int*, int*, int*, int*, int*, int, int);
// // 
// // template <typename T> void cpuCreateVerletList(int *verletlist, T *x, T *ellipsoid, int *verletnum,  int *verletnumsum, 
// //      int *atomtype, int *ilist, int *clist, int *c2alist, int *c2anum, int *c2anumsum, int *nc, int inum, int dim)
// // {
// //     // a list contains the starting positions of the first j atom 
// //     cpuCumsum(verletnumsum, verletnum, inum+1); 
// //     
// //     int dimsq = dim*dim;
// //     T onet = (T) 1.0;
// //     if (dim==2) {
// //         for (int ii=0; ii<inum; ii++) {  // for each atom i in the simulation box    
// //             int i = ilist[ii];
// //             T xi0 = x[i*dim];        // position of atom i
// //             T xi1 = x[i*dim+1];      // position of atom i
// //             int t = atomtype[i];      // element of atom i 
// //             T A00 = ellipsoid[dimsq*t]; // ellipsoid for element t   
// //             T A10 = ellipsoid[dimsq*t+1]; // ellipsoid for element t   
// //             T A01 = ellipsoid[dimsq*t+2]; // ellipsoid for element t   
// //             T A11 = ellipsoid[dimsq*t+3]; // ellipsoid for element t   
// //             int nstart = verletnumsum[i];     // starting 
// //             int ninc = 0;                // increment
// //             int j = clist[i];         // cell j contains atom i           
// //             int j1 = j%nc[0];
// //             int j2 = (j-j1)/nc[0];
// //             for (int i1=-1; i1<=1; i1++) {
// //                 int k1 = j1 + i1;
// //                 for (int i2=-1; i2<=1; i2++) {
// //                     int k2 = j2 + i2;
// //                     int k = k1 + nc[0]*k2; // cell k
// //                     int m = c2anum[k];     // number of atoms in cell k
// //                     int s = c2anumsum[k];  // starting position of the first atom in cell k
// //                     for (int l=0; l<m ; l++) { //loop over each atom j in cell k
// //                         j = c2alist[s+l];  // atom j
// //                         T xij0 = x[j*dim] - xi0;  // xj - xi
// //                         T xij1 = x[j*dim+1] - xi1; // xj - xi
// //                         // distance between atom i and atom j 
// //                         T rij = xij0*(A00*xij0 + A01*xij1) + xij1*(A10*xij0 + A11*xij1);                        
// //                         if (rij <= onet) {
// //                             verletlist[nstart + ninc] = j; // add atom j into the list
// //                             ninc += 1;
// //                         }
// //                     }
// //                 }
// //             }               
// //             verletnum[i] = ninc;
// //         }        
// //     }
// //     else if (dim==3) {
// //         for (int ii=0; ii<inum; ii++) {  // for each atom i in the simulation box    
// //             int i = ilist[ii];
// //             T xi0 = x[i*dim];        // position of atom i
// //             T xi1 = x[i*dim+1];      // position of atom i
// //             T xi2 = x[i*dim+2];      // position of atom i
// //             int t = atomtype[i];      // element of atom i 
// //             T A00 = ellipsoid[dimsq*t]; // ellipsoid for element t   
// //             T A10 = ellipsoid[dimsq*t+1]; // ellipsoid for element t   
// //             T A20 = ellipsoid[dimsq*t+2]; // ellipsoid for element t   
// //             T A01 = ellipsoid[dimsq*t+3]; // ellipsoid for element t   
// //             T A11 = ellipsoid[dimsq*t+4]; // ellipsoid for element t   
// //             T A21 = ellipsoid[dimsq*t+5]; // ellipsoid for element t   
// //             T A02 = ellipsoid[dimsq*t+6]; // ellipsoid for element t   
// //             T A12 = ellipsoid[dimsq*t+7]; // ellipsoid for element t   
// //             T A22 = ellipsoid[dimsq*t+8]; // ellipsoid for element t   
// //             int nstart = verletnumsum[i];     // starting 
// //             int ninc = 0;                // increment            
// //             int j = clist[i];         // cell j contains atom i       
// //             int n = j%(nc[0]*nc[1]);
// //             int j3 = (j-n)/(nc[0]*nc[1]);            
// //             int j1 = n%nc[0];
// //             int j2 = (j-j1)/nc[0];
// //             for (int i1=-1; i1<=1; i1++) {
// //                 int k1 = j1 + i1;
// //                 for (int i2=-1; i2<=1; i2++) {
// //                     int k2 = j2 + i2;
// //                     for (int i3=-1; i3<=1; i3++) {
// //                         int k3 = j3 + i3;                    
// //                         int k = k1 + nc[0]*k2 + nc[0]*nc[1]*k3; // cell k
// //                         int m = c2anum[k];     // number of atoms in cell k
// //                         int s = c2anumsum[k];  // starting position of the first atom in cell k
// //                         for (int l=0; l<m ; l++) {
// //                             j = c2alist[s+l];  // atom j
// //                             T xij0 = x[j*dim] - xi0;  // xj - xi
// //                             T xij1 = x[j*dim+1] - xi1; // xj - xi
// //                             T xij2 = x[j*dim+2] - xi2; // xj - xi
// //                             // distance between atom i and atom j 
// //                             T rij = xij0*(A00*xij0 + A01*xij1 + A02*xij2) + xij1*(A10*xij0 + A11*xij1 + A12*xij2) + xij2*(A20*xij0 + A21*xij1 + A22*xij2);
// //                             if (rij <= onet) {
// //                                 verletlist[nstart + ninc] = j; // add atom j into the list
// //                                 ninc += 1;
// //                             }                            
// //                         }
// //                     }
// //                 }
// //             }                
// //             verletnum[i] = ninc;
// //         }        
// //     }        
// // }
// // template void cpuCreateVerletList(int*, double*, double*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int, int);
// // template void cpuCreateVerletList(int*, float*, float*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int, int);
// // 
// // template <typename T> void cpuFullNeighNum(int *neighnum, T *x, T* ellipsoid, int *atomtype, 
// //         int *ilist, int *verletlist, int *verletnum, int *verletnumsum, int inum, int dim)
// // {
// //     int dimsq = dim*dim;
// //     
// //     if (dim==2) {    
// //         for (int ii=0; ii<inum; ii++) {  // for each atom i in the simulation box     
// //             int i = ilist[ii];  // atom i
// //             T xi0 = x[i*dim];        // position of atom i
// //             T xi1 = x[i*dim+1];      // position of atom i
// //             int t = atomtype[i];      // element of atom i 
// //             T A00 = ellipsoid[dimsq*t]; // ellipsoid for element t   
// //             T A10 = ellipsoid[dimsq*t+1]; // ellipsoid for element t   
// //             T A01 = ellipsoid[dimsq*t+2]; // ellipsoid for element t   
// //             T A11 = ellipsoid[dimsq*t+3]; // ellipsoid for element t      
// //             int m = verletnum[i]; // number of atoms around i 
// //             int start = verletnumsum[i];     // starting 
// //             int ninc = 0;                // increment 
// //             for (int l=0; l<m ; l++) {   // loop over each atom around atom i
// //                 int j = verletlist[start+l];  // atom j                
// //                 T xij0 = x[j*dim] - xi0;  // xj - xi
// //                 T xij1 = x[j*dim+1] - xi1; // xj - xi
// //                 // distance between atom i and atom j 
// //                 T rij = xij0*(A00*xij0 + A01*xij1) + xij1*(A10*xij0 + A11*xij1);                        
// //                 if (rij <= 1.0) // if atom j is inside the ellipsoid 
// //                     ninc += 1;  // increase the number of neighbors by 1               
// //             }
// //             neighnum[i] = ninc; // number of neighbors of atom i
// //         }
// //     }
// //     else if (dim==3) {
// //         for (int ii=0; ii<inum; ii++) {  // for each atom i in the simulation box         
// //             int i = ilist[ii];  // atom i
// //             T xi0 = x[i*dim];        // position of atom i
// //             T xi1 = x[i*dim+1];      // position of atom i
// //             T xi2 = x[i*dim+2];      // position of atom i
// //             int t = atomtype[i];      // element of atom i 
// //             T A00 = ellipsoid[dimsq*t]; // ellipsoid for element t   
// //             T A10 = ellipsoid[dimsq*t+1]; // ellipsoid for element t   
// //             T A20 = ellipsoid[dimsq*t+2]; // ellipsoid for element t   
// //             T A01 = ellipsoid[dimsq*t+3]; // ellipsoid for element t   
// //             T A11 = ellipsoid[dimsq*t+4]; // ellipsoid for element t   
// //             T A21 = ellipsoid[dimsq*t+5]; // ellipsoid for element t   
// //             T A02 = ellipsoid[dimsq*t+6]; // ellipsoid for element t   
// //             T A12 = ellipsoid[dimsq*t+7]; // ellipsoid for element t   
// //             T A22 = ellipsoid[dimsq*t+8]; // ellipsoid for element t   
// //             int m = verletnum[i]; // number of atoms around i 
// //             int start = verletnumsum[i];     // starting 
// //             int ninc = 0;                // increment 
// //             for (int l=0; l<m ; l++) {   // loop over each atom around atom i
// //                 int j = verletlist[start+l];  // atom j                
// //                 T xij0 = x[j*dim] - xi0;  // xj - xi
// //                 T xij1 = x[j*dim+1] - xi1; // xj - xi
// //                 T xij2 = x[j*dim+2] - xi2; // xj - xi
// //                 // distance between atom i and atom j 
// //                 T rij = xij0*(A00*xij0 + A01*xij1 + A02*xij2) + xij1*(A10*xij0 + A11*xij1 + A12*xij2) + xij2*(A20*xij0 + A21*xij1 + A22*xij2);
// //                 if (rij <= 1.0) // if atom j is inside the ellipsoid 
// //                     ninc += 1;  // increase the number of neighbors by 1                               
// //             }            
// //             neighnum[i] = ninc; // number of neighbors of atom i
// //         }
// //     }
// // }
// // template void cpuFullNeighNum(int*, double*, double*, int*, int*, int*, int*, int*, int, int);
// // template void cpuFullNeighNum(int*, float*, float*, int*, int*, int*, int*, int*, int, int);
// // 
// // template <typename T> void cpuFullNeighList(int *neighlist, T *x, T* ellipsoid, int *neighnum, int *neighnumsum,
// //         int *atomtype, int *ilist, int *verletlist, int *verletnum, int *verletnumsum, int inum, int dim)
// // {
// //     // a list contains the starting positions of the first neighbor 
// //     cpuCumsum(neighnumsum, neighnum, inum+1); 
// //     
// //     int dimsq = dim*dim;
// //     
// //     if (dim==2) {    
// //         for (int ii=0; ii<inum; ii++) {  // for each atom i in the simulation box     
// //             int i = ilist[ii];  // atom i
// //             T xi0 = x[i*dim];        // position of atom i
// //             T xi1 = x[i*dim+1];      // position of atom i
// //             int t = atomtype[i];      // element of atom i 
// //             T A00 = ellipsoid[dimsq*t]; // ellipsoid for element t   
// //             T A10 = ellipsoid[dimsq*t+1]; // ellipsoid for element t   
// //             T A01 = ellipsoid[dimsq*t+2]; // ellipsoid for element t   
// //             T A11 = ellipsoid[dimsq*t+3]; // ellipsoid for element t      
// //             int m = verletnum[i]; // number of atoms around i 
// //             int jstart = verletnumsum[i];     // starting 
// //             int nstart = neighnumsum[i];   
// //             int ninc = 0;                // increment 
// //             for (int l=0; l<m ; l++) {   // loop over each atom around atom i
// //                 int j = verletlist[jstart+l];  // atom j                
// //                 T xij0 = x[j*dim] - xi0;  // xj - xi
// //                 T xij1 = x[j*dim+1] - xi1; // xj - xi
// //                 // distance between atom i and atom j 
// //                 T rij = xij0*(A00*xij0 + A01*xij1) + xij1*(A10*xij0 + A11*xij1);                        
// //                 if (rij <= 1.0) { // if atom j is inside the ellipsoid 
// //                     neighlist[nstart+ninc] = j;
// //                     ninc += 1;  // increase the number of neighbors by 1               
// //                 }
// //             }
// //             neighnum[i] = ninc; // number of neighbors of atom i
// //         }
// //     }
// //     else if (dim==3) {
// //         for (int ii=0; ii<inum; ii++) {  // for each atom i in the simulation box         
// //             int i = ilist[ii];  // atom i
// //             T xi0 = x[i*dim];        // position of atom i
// //             T xi1 = x[i*dim+1];      // position of atom i
// //             T xi2 = x[i*dim+2];      // position of atom i
// //             int t = atomtype[i];      // element of atom i 
// //             T A00 = ellipsoid[dimsq*t]; // ellipsoid for element t   
// //             T A10 = ellipsoid[dimsq*t+1]; // ellipsoid for element t   
// //             T A20 = ellipsoid[dimsq*t+2]; // ellipsoid for element t   
// //             T A01 = ellipsoid[dimsq*t+3]; // ellipsoid for element t   
// //             T A11 = ellipsoid[dimsq*t+4]; // ellipsoid for element t   
// //             T A21 = ellipsoid[dimsq*t+5]; // ellipsoid for element t   
// //             T A02 = ellipsoid[dimsq*t+6]; // ellipsoid for element t   
// //             T A12 = ellipsoid[dimsq*t+7]; // ellipsoid for element t   
// //             T A22 = ellipsoid[dimsq*t+8]; // ellipsoid for element t   
// //             int m = verletnum[i]; // number of atoms around i 
// //             int jstart = verletnumsum[i];     // starting 
// //             int nstart = neighnumsum[i];   
// //             int ninc = 0;                // increment 
// //             for (int l=0; l<m ; l++) {   // loop over each atom around atom i
// //                 int j = verletlist[jstart+l];  // atom j                
// //                 T xij0 = x[j*dim] - xi0;  // xj - xi
// //                 T xij1 = x[j*dim+1] - xi1; // xj - xi
// //                 T xij2 = x[j*dim+2] - xi2; // xj - xi
// //                 // distance between atom i and atom j 
// //                 T rij = xij0*(A00*xij0 + A01*xij1 + A02*xij2) + xij1*(A10*xij0 + A11*xij1 + A12*xij2) + xij2*(A20*xij0 + A21*xij1 + A22*xij2);
// //                 if (rij <= 1.0) { // if atom j is inside the ellipsoid 
// //                     neighlist[nstart+ninc] = j;
// //                     ninc += 1;  // increase the number of neighbors by 1               
// //                 }
// //             }            
// //             neighnum[i] = ninc; // number of neighbors of atom i
// //         }
// //     }
// // }
// // template void cpuFullNeighList(int*, double*, double*, int*, int*, int*, int*, int*, int*, int*, int, int);
// // template void cpuFullNeighList(int*, float*, float*, int*, int*, int*, int*, int*, int*, int*, int, int);
// 
// 

#endif


