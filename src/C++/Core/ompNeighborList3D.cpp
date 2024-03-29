/***************************************************************************                               
                    Molecular Dynamics Potentials (MDP)
                           CESMIX-MIT Project  
 
 Contributing authors: Ngoc-Cuong Nguyen (cuongng@mit.edu, exapde@gmail.com)
 ***************************************************************************/

#ifndef MDP_OMPNEIGHBORLIST3D
#define MDP_OMPNEIGHBORLIST3D

template <typename T> void ompBoundingBox3D(T *vc, T *wc, T *v, T *w, T *a, T *b, T *c, T *r, int *pbc)
{
    T a1 = a[0]; 
    T a2 = a[1];
    T a3 = a[2];
    T b1 = b[0]; 
    T b2 = b[1];
    T b3 = b[2];
    T c1 = c[0]; 
    T c2 = c[1];
    T c3 = c[2];
    
    T norma = sqrt(ompArraySquareSum(a, 3));
    T normb = sqrt(ompArraySquareSum(b, 3));
    T normc = sqrt(ompArraySquareSum(c, 3));
    
    // 8 vertices of the parallelogram defined by a, b, c
    v[0] = 0.0;
    v[1] = 0.0;
    v[2] = 0.0;
    v[3] = a1;
    v[4] = a2;
    v[5] = a3;    
    v[6] = a1 + b1;
    v[7] = a2 + b2;
    v[8] = a3 + b3;    
    v[9] = b1;
    v[10] = b2;
    v[11] = b3;
    v[12] = c1;
    v[13] = c2;
    v[14] = c3;
    v[15] = a1 + c1;
    v[16] = a2 + c2;
    v[17] = a3 + c3;    
    v[18] = a1 + b1 + c1;
    v[19] = a2 + b2 + c2;
    v[20] = a3 + b3 + c3;    
    v[21] = b1 + c1;
    v[22] = b2 + c2;
    v[23] = b3 + c3;
    
    // the 1st plane defined by a and b
    T p1[3] = {a2*b3-a3*b2, a3*b1-a1*b3, a1*b2-a2*b1};
    T normp1 = sqrt(ompArraySquareSum(p1, 3));    
    // Since the distance between (0,0,0) and the 1st plane p1(1)*x + p1(2)*y + p1(3)*z = d1
    // is equal to r, we have
    T d1 = -r[2]*normp1;
    
    // the 2nd plane defined by b and c
    T p2[3] = {b2*c3-b3*c2, b3*c1-b1*c3, b1*c2-b2*c1};
    T normp2 = sqrt(ompArraySquareSum(p2, 3));
    T d2 = -r[0]*normp2;
    
    // the 3rd plane defined by c and a
    T p3[3] = {c2*a3-c3*a2, c3*a1-c1*a3, c1*a2-c2*a1};
    T normp3 = sqrt(ompArraySquareSum(p3, 3));
    T d3 = -r[1]*normp3;
        
    // intersection of the above three planes
    //w1 = [p1; p2; p3]\[d1; d2; d3];
    T A[9] = {p1[0], p1[1], p1[2], p2[0], p2[1], p2[2], p3[0], p3[1], p3[2]};
    T invA[9];    
    ompSmallMatrixInverse(invA, A, 3);        
    w[0] = invA[0]*d1 + invA[3]*d2 + invA[6]*d3;
    w[1] = invA[1]*d1 + invA[4]*d2 + invA[7]*d3;
    w[2] = invA[2]*d1 + invA[5]*d2 + invA[8]*d3;
    
    // find e = (e1,e2,e3) such that e1*a + e2*b + e3*c = w1
    //e = [a(:) b(:) c(:)]\w1;    
    ompSmallMatrixInverse(invA, a, b, c);
    T e1 = invA[0]*w[0] + invA[3]*w[1] + invA[6]*w[2];
    T e2 = invA[1]*w[0] + invA[4]*w[1] + invA[7]*w[2];
    T e3 = invA[2]*w[0] + invA[5]*w[1] + invA[8]*w[2];
    
    T onet = (T) 1.0;
    T sa[3], sb[3], sc[3];
    // distance between w1 and the point e2*b + e3*c    
    ompArrayAdd3Vectors(sa, w, b, c, onet, -e2, -e3, 3);
    
    // distance between w1 and the point e1*a + e3*c
    ompArrayAdd3Vectors(sb, w, a, c, onet, -e1, -e3, 3);

    // distance between w1 and the point e1*a + e2*b
    ompArrayAdd3Vectors(sc, w, a, b, onet, -e1, -e2, 3);
    
    // length of the bounding parallelepiped along the 1st axis
    T la = norma + (2.0*pbc[0])*sqrt(ompArraySquareSum(sa, 3));
    // length of the bounding parallelepiped along the 2nd axis
    T lb = normb + (2.0*pbc[1])*sqrt(ompArraySquareSum(sb, 3));
    // length of the bounding parallelepiped along the 3rd axis
    T lc = normc + (2.0*pbc[2])*sqrt(ompArraySquareSum(sc, 3));
    
    // the 1st vertex of the bounding parallelepiped   
    ompArrayAdd3Vectors(&w[0], a, b, c, pbc[0]*e1, pbc[1]*e2, pbc[2]*e3, 3);    
//     w[0] = pbc[0]*e[0]*a1 + pbc[1]*e[1]*b1 + pbc[c]*e[2]*c1;
//     w[1] = pbc[0]*e[0]*a2 + pbc[1]*e[1]*b2 + pbc[c]*e[2]*c2;
//     w[2] = pbc[0]*e[0]*a3 + pbc[1]*e[1]*b3 + pbc[c]*e[2]*c3;    
    // the 2nd vertex of the bounding parallelepiped: w2 = w1 + la*a/norma;
    ompArrayAXPBY(&w[3], &w[0], a, onet, la/norma, 3);
    // the 4th vertex of the bounding parallelepiped: w4 = w1 + lb*b/normb;
    ompArrayAXPBY(&w[9], &w[0], b, onet, lb/normb, 3);
    // the 3rd vertex of the bounding parallelepiped: w3 = w2 + w4 - w1;
    ompArrayAdd3Vectors(&w[6], &w[0], &w[3], &w[9], -onet, onet, onet, 3);    
    // the 5th vertex of the bounding parallelepiped: w5 = w1 + lc*c/normc;
    ompArrayAXPBY(&w[12], &w[0], c, onet, lc/normc, 3);
    // the 6th vertex of the bounding parallelepiped: w6 = w5 + la*a/norma;
    ompArrayAXPBY(&w[15], &w[12], a, onet, la/norma, 3);
    // the 8th vertex of the bounding parallelepiped: w8 = w5 + lb*b/normb;
    ompArrayAXPBY(&w[21], &w[12], b, onet, lb/normb, 3);
    // the 3rd vertex of the bounding parallelepiped: w7 = w6 + w8 - w5;
    ompArrayAdd3Vectors(&w[18], &w[12], &w[15], &w[21], -onet, onet, onet, 3);    
    
    for (int i=0; i<8; i++) {
        vc[3*i+0] = invA[0]*v[3*i+0] + invA[3]*v[3*i+1] + invA[6]*v[3*i+2];
        vc[3*i+1] = invA[1]*v[3*i+0] + invA[4]*v[3*i+1] + invA[7]*v[3*i+2];
        vc[3*i+2] = invA[2]*v[3*i+0] + invA[5]*v[3*i+1] + invA[8]*v[3*i+2];
        wc[3*i+0] = invA[0]*w[3*i+0] + invA[3]*w[3*i+1] + invA[6]*w[3*i+2];
        wc[3*i+1] = invA[1]*w[3*i+0] + invA[4]*w[3*i+1] + invA[7]*w[3*i+2];
        wc[3*i+2] = invA[2]*w[3*i+0] + invA[5]*w[3*i+1] + invA[8]*w[3*i+2];        
    }    
}
template void ompBoundingBox3D(double*, double*, double*, double*, double*, double*, double*, double*, int*);
template void ompBoundingBox3D(float*, float*, float*, float*, float*, float*, float*, float*, int*);
    
template <typename T> int ompPeriodicImages3D(T *pimages, T *a, T *b, T *c, int *pbc)
{
    T px[3][27];
    
    px[0][0] = 0.0;             px[1][0] = 0.0;             px[2][0] = 0.0;
    px[0][1] = a[0];            px[1][1] = a[1];            px[2][1] = a[2];
    px[0][2] =-a[0];            px[1][2] =-a[1];            px[2][2] =-a[2];
    px[0][3] = b[0];            px[1][3] = b[1];            px[2][3] = b[2];
    px[0][4] =-b[0];            px[1][4] =-b[1];            px[2][4] =-b[2];
    px[0][5] = a[0]+b[0];       px[1][5] = a[1]+b[1];       px[2][5] = a[2]+b[2];
    px[0][6] = a[0]-b[0];       px[1][6] = a[1]-b[1];       px[2][6] = a[2]-b[2];
    px[0][7] =-a[0]+b[0];       px[1][7] =-a[1]+b[1];       px[2][7] =-a[2]+b[2];
    px[0][8] =-a[0]-b[0];       px[1][8] =-a[1]-b[1];       px[2][8] =-a[2]-b[2];
    
    px[0][9] =  c[0]+0.0;       px[1][9] =  c[1]+0.0;       px[2][9] =  c[2]+0.0;
    px[0][10] = c[0]+a[0];      px[1][10] = c[1]+a[1];      px[2][10] = c[2]+a[2];
    px[0][11] = c[0]-a[0];      px[1][11] = c[1]-a[1];      px[2][11] = c[2]-a[2];
    px[0][12] = c[0]+b[0];      px[1][12] = c[1]+b[1];      px[2][12] = c[2]+b[2];
    px[0][13] = c[0]-b[0];      px[1][13] = c[1]-b[1];      px[2][13] = c[2]-b[2];
    px[0][14] = c[0]+a[0]+b[0]; px[1][14] = c[1]+a[1]+b[1]; px[2][14] = c[2]+a[2]+b[2];
    px[0][15] = c[0]+a[0]-b[0]; px[1][15] = c[1]+a[1]-b[1]; px[2][15] = c[2]+a[2]-b[2];
    px[0][16] = c[0]-a[0]+b[0]; px[1][16] = c[1]-a[1]+b[1]; px[2][16] = c[2]-a[2]+b[2];
    px[0][17] = c[0]-a[0]-b[0]; px[1][17] = c[1]-a[1]-b[1]; px[2][17] = c[2]-a[2]-b[2];

    px[0][18] =-c[0]+0.0;       px[1][18] =-c[1]+0.0;       px[2][18] =-c[2]+0.0;
    px[0][19] =-c[0]+a[0];      px[1][19] =-c[1]+a[1];      px[2][19] =-c[2]+a[2];
    px[0][20] =-c[0]-a[0];      px[1][20] =-c[1]-a[1];      px[2][20] =-c[2]-a[2];
    px[0][21] =-c[0]+b[0];      px[1][21] =-c[1]+b[1];      px[2][21] =-c[2]+b[2];
    px[0][22] =-c[0]-b[0];      px[1][22] =-c[1]-b[1];      px[2][22] =-c[2]-b[2];
    px[0][23] =-c[0]+a[0]+b[0]; px[1][23] =-c[1]+a[1]+b[1]; px[2][23] =-c[2]+a[2]+b[2];
    px[0][24] =-c[0]+a[0]-b[0]; px[1][24] =-c[1]+a[1]-b[1]; px[2][24] =-c[2]+a[2]-b[2];
    px[0][25] =-c[0]-a[0]+b[0]; px[1][25] =-c[1]-a[1]+b[1]; px[2][25] =-c[2]-a[2]+b[2];
    px[0][26] =-c[0]-a[0]-b[0]; px[1][26] =-c[1]-a[1]-b[1]; px[2][26] =-c[2]-a[2]-b[2];
    
    int nimages;
    int indx[27];
    if ((pbc[0]==0) && (pbc[1]==0) && (pbc[2]==0)) {
         nimages = 1;
         int ind[] = {0};
         for (int i=0; i<nimages; i++)
            indx[i] = ind[i];
    }
    else if ((pbc[0]==1) && (pbc[1]==0) && (pbc[2]==0)) {
        nimages = 3;
        int ind[] = {0,1,2};      
        for (int i=0; i<nimages; i++)
            indx[i] = ind[i];
    }
    else if ((pbc[0]==0) && (pbc[1]==1) && (pbc[2]==0)) {
        nimages = 3;
        int ind[] = {0,3,4};    
        for (int i=0; i<nimages; i++)
            indx[i] = ind[i];
    }
    else if ((pbc[0]==0) && (pbc[1]==0) && (pbc[2]==1)) {
        nimages = 3;
        int ind[] = {0,9,18};
        for (int i=0; i<nimages; i++)
            indx[i] = ind[i];
    }
    else if ((pbc[0]==1) && (pbc[1]==1) && (pbc[2]==0)) {
        nimages = 9;
        int ind[] = {0,1,2,3,4,5,6,7,8};
        for (int i=0; i<nimages; i++)
            indx[i] = ind[i];
    }
    else if ((pbc[0]==1) && (pbc[1]==0) && (pbc[2]==1)) {
        nimages = 9;
        int ind[] = {0,1,2,9,10,11,18,19,20};
        for (int i=0; i<nimages; i++)
            indx[i] = ind[i];
    }
    else if ((pbc[0]==0) && (pbc[1]==1) && (pbc[2]==1)) {
        nimages = 9;
        int ind[] = {0,3,4,9,12,13,18,21,22};
        for (int i=0; i<nimages; i++)
            indx[i] = ind[i];
    }
    else if ((pbc[0]==1) && (pbc[1]==1) && (pbc[2]==1)) {
        nimages = 27;        
        for (int i=0; i<nimages; i++)
            indx[i] = i;
    }
    else {
    }
    
    for (int i=0; i<nimages; i++) {
        pimages[3*i+0] = px[0][indx[i]];
        pimages[3*i+1] = px[1][indx[i]];
        pimages[3*i+2] = px[2][indx[i]];
    }    
    
    return nimages;
}
template int ompPeriodicImages3D(double*, double*, double*, double*, int*);
template int ompPeriodicImages3D(float*, float*, float*, float*, int*);

/************************************************************************************/

void ompGridColoring3D(int *cellcolor, int *nc, int *bsize, int dim)
{
    int nc1 = nc[0];
    int nc2 = nc[1];
    int nc3 = nc[2];
    int bs1 = bsize[0];
    int bs2 = bsize[1];
    int bs3 = bsize[2];        
    int nb1 = ceil(((double) nc1)/((double) bs1)); // number of blocks
    int nb2 = ceil(((double) nc2)/((double) bs2)); // number of blocks
    int nb3 = ceil(((double) nc3)/((double) bs3)); // number of blocks

    int bid[bs1*bs2*bs3];
    for (int i=0; i<bs1*bs2*bs3; i++)            
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
            if (j==1) 
                j1 = (j==1) ? 0 : (j1 + bs2 - 1);
            j2 = j1 + bs2 - 1;    
            if (j2 >= nc2)
                j2 = nc2 - 1;                              
            for (int k=0; k<nb3; k++) {
                int k1, k2;
                if (k==1) 
                    k1 = (k==1) ? 0 : (k1 + bs3 - 1);
                k2 = k1 + bs3 - 1;    
                if (k2 >= nc3)
                    k2 = nc3 - 1;          

                for (int m3=k1; m3<=k2; m3++)
                    for (int m2=j1; m2<=j2; m2++)
                        for (int m1=i1; m1<=i2; m1++)                    
                            cellcolor[m1 + nc[0]*m2 + nc[0]*nc[1]*m3] = bid[m1-i1 + (m2-j1)*bs1 + (m3-k1)*bs1*bs2];                
            }
        }
    }                                    
}

template <typename T> void ompGhostAtoms3D(int *glistnum, int *inside, T *x, T *pimages, T *wc, T *s2rmap, int n, int pnum, int dim)
{        
    #pragma omp parallel for
    for (int i=0; i<n; i++) { // loop over each atom i inside the simulation box
        glistnum[i] = 0; // set the number of ghost atoms for atom i to 0
        for (int j=1; j<pnum; j++) { // loop over each periodic image of atom i
            T xj0 = x[i*dim+0] + pimages[j*dim+0];    // periodic image of x          
            T xj1 = x[i*dim+1] + pimages[j*dim+1];        
            T xj2 = x[i*dim+2] + pimages[j*dim+2];        
            T xc0 = s2rmap[0]*xj0 + s2rmap[3]*xj1 + s2rmap[6]*xj2;  // map it to the unit square            
            T xc1 = s2rmap[1]*xj0 + s2rmap[4]*xj1 + s2rmap[7]*xj2;        
            T xc2 = s2rmap[2]*xj0 + s2rmap[5]*xj1 + s2rmap[8]*xj2;   
                        /// check if the mapped point is inside the bounding box
            if ((wc[0] <= xc0) && (xc0 <= wc[18]) && (wc[1] <= xc1) && (xc1 <= wc[19]) && (wc[2] <= xc2) && (xc2 <= wc[20])) { 
                glistnum[i] += 1;            // increase the number of ghost atoms by 1 
                inside[j-1 + i*(pnum-1)] = 1;
            }
            else
                inside[j-1 + i*(pnum-1)] = 0;
        }
    }              
}
template void ompGhostAtoms3D(int*, int*, double*, double*, double*, double*, int, int, int);
template void ompGhostAtoms3D(int*, int*, float*, float*, float*, float*, int, int, int);

template <typename T> void ompAtomList3D(int *alist, int *inside, int *glistnumsum, int *glistnum, 
        T *x, T *pimages, T *wc, T *s2rmap, int inum, int pnum, int dim)
{           
    // a list of ghost atoms
    ompGhostAtoms3D(glistnum, inside, x, pimages, wc, s2rmap, inum, pnum, dim);
    
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
                x[dim*(q+k)+2] = x[i*dim+2] + pimages[j*dim+2]; //
                // q + k = inum + glistnumsum[i] + k
                alist[q+k] = i;  // add atom i to the list                
                k += 1;          // increase the number of ghost atoms for atom i by 1    
            }
        }
    }    
}
template void ompAtomList3D(int*, int*, int*, int*, double*, double*, double*, double*, int, int, int);
template void ompAtomList3D(int*, int*, int*, int*, float*, float*, float*, float*, int, int, int);

template <typename T> void ompCellList3D(int *clist, T *x, T *eta1, T *eta2, T *eta3, T *s2rmap, 
        int *nc, int inum, int natom, int dim)
{    
    #pragma omp parallel for
    for (int i=0; i<natom; i++) { // loop over each atom i
        int j1, j2, j3, m;
        int k = dim*i;
        m = (i < inum) ? 1 : 0;
        // position of atom i in the unit square/cube
        T xt0 = s2rmap[0]*x[k+0] + s2rmap[3]*x[k+1] + s2rmap[6]*x[k+2];        
        T xt1 = s2rmap[1]*x[k+0] + s2rmap[4]*x[k+1] + s2rmap[7]*x[k+2];                                    
        T xt2 = s2rmap[2]*x[k+0] + s2rmap[5]*x[k+1] + s2rmap[8]*x[k+2];                
        // identify a cell containing atom i 
        for (j1=m; j1<nc[0]-m; j1++)
            if ((eta1[j1] <= xt0) && (xt0<= eta1[j1+1]))
                break;
        for (j2=m; j2<nc[1]-m; j2++) 
            if ((eta2[j2] <= xt1) && (xt1<= eta2[j2+1]))
                break;
        for (j3=m; j3<nc[2]-m; j3++) 
            if ((eta3[j3] <= xt2) && (xt2<= eta3[j3+1]))
                break;
        clist[i] = j1 + nc[0]*j2 + nc[0]*nc[1]*j3; // link that cell to atom i                
    }                        
}
template void ompCellList3D(int*, double*, double*, double*, double*, double*, int*, int, int, int);
template void ompCellList3D(int*, float*, float*, float*, float*, float*, int*, int, int, int);

template <typename T> void ompFullNeighborList3D(int *neighborlist, int *neighbornum, T *x, T *rcutsq, int *alist, 
        int *clist, int *c2alist, int *c2anumsum, int *nc, int inum, int jnum, int dim)
{            
    #pragma omp parallel for
    for (int ii=0; ii<inum; ii++) {  // for each atom i in the simulation box    
        int i = alist[ii];        
        T xi0 = x[i*dim];        // position of atom i
        T xi1 = x[i*dim+1];      // position of atom i        
        T xi2 = x[i*dim+2];      // position of atom i
        int ja = clist[i];         // cell j contains atom i           
        int jb = ja%(nc[0]*nc[1]);
        int j3 = (ja-jb)/(nc[0]*nc[1]);
        int j1 = jb%nc[0];
        int j2 = (jb-j1)/nc[0];
        
        int inc = 0;                // increment
        for (int i1=-1; i1<=1; i1++) {
            int k1 = j1 + i1;
            for (int i2=-1; i2<=1; i2++) {
                int k2 = j2 + i2;
                for (int i3=-1; i3<=1; i3++) {
                    int k3 = j3 + i3;                
                    int k = k1 + nc[0]*k2 + nc[0]*nc[1]*k3; // cell k
                    //int m = c2anum[k];     // number of atoms in cell k
                    int s = c2anumsum[k];    // starting position of the first atom in cell k
                    int m = c2anumsum[k+1]-s;  // number of atoms in cell k
                    for (int l=0; l<m ; l++) {
                        int j = c2alist[s+l];  // atom j
                        if (j != i) {
                            T xij0 = x[j*dim] - xi0;  // xj - xi
                            T xij1 = x[j*dim+1] - xi1; // xj - xi
                            T xij2 = x[j*dim+2] - xi2; // xj - xi
                            // distance between atom i and atom j 
                            T rijsq = xij0*xij0 + xij1*xij1 + xij2*xij2;          
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
            if (inc==jnum)
                break;                                
        }    
        neighbornum[ii] = inc;
    }        
}
template void ompFullNeighborList3D(int*, int*, double*, double*, int*, int*, int*, int*, int*, int, int, int);
template void ompFullNeighborList3D(int*, int*, float*, float*, int*, int*, int*, int*, int*, int, int, int);

template <typename T> void ompFullNeighborList3D(int *neighborlist, int *neighbornum, T *x, T *rcutsq, 
        int anum, int inum, int jnum, int dim)
{                
    #pragma omp parallel for
    for (int i=0; i<inum; i++) {  // for each atom i in the simulation box    
        T xi0 = x[i*dim];        // position of atom i
        T xi1 = x[i*dim+1];      // position of atom i     
        T xi2 = x[i*dim+2];      // position of atom i
        int inc = 0;                // increment
        for (int j=0; j<anum; j++) { // loop over atom j (including both local and ghost atoms)            
            if (j != i) {
                T xij0 = x[j*dim] - xi0;  // xj - xi
                T xij1 = x[j*dim+1] - xi1; // xj - xi               
                T xij2 = x[j*dim+2] - xi2; // xj - xi
                T rijsq = xij0*xij0 + xij1*xij1 + xij2*xij2;  // distance between atom i and atom j                        
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
template void ompFullNeighborList3D(int*, int*, double*, double*, int, int, int, int);
template void ompFullNeighborList3D(int*, int*, float*, float*, int, int, int, int);


#endif


