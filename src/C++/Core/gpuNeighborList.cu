#ifndef __GPUNEIGHBORLIST
#define __GPUNEIGHBORLIST

//template <typename T> void cpuGhostAtoms2D(int *glistnum, int *inside, T *x, T *pimages, T *wc, T *s2rmap, int n, int pnum, int dim)
template <typename T>
__global__ void gpuKernelGhostAtoms2D(int *glistnum, T *x, T *pimages, T *wc, T *B2C, int n, int m, int dim)
{
    int i = threadIdx.x + blockIdx.x * blockDim.x;
    while (i < n) {
        glistnum[i] = 0; // set the number of ghost atoms for atom i to 0
        for (int j=1; j<m; j++) { // loop over each periodic image of atom i
            T xj0 = x[i*dim+0] + pimages[j*dim+0];  // periodic image of x      
            T xj1 = x[i*dim+1] + pimages[j*dim+1];        
            T xc0 = B2C[0]*xj0 + B2C[2]*xj1;        // map it to the unit square      
            T xc1 = B2C[1]*xj0 + B2C[3]*xj1;                        
            /// check if the mapped point is inside the bounding box
            if ((wc[0] <= xc0) && (xc0 <= wc[4]) &&  (wc[1] <= xc1) && (xc1 <= wc[5])) 
                glistnum[i] += 1;            // increase the number of ghost atoms by 1                    
        }        
        i += blockDim.x * gridDim.x;
    }
}

template <typename T>
__global__ void gpuKernelGhostAtoms3D(int *glistnum, T *x, T *pimages, T *wc, T *B2C, int n, int m, int dim)
{
    int i = threadIdx.x + blockIdx.x * blockDim.x;
    while (i < n) {
        glistnum[i] = 0; // set the number of ghost atoms for atom i to 0
        for (int j=1; j<m; j++) { // loop over each periodic image of atom i
            T xj0 = x[i*dim+0] + pimages[j*dim+0];    // periodic image of x          
            T xj1 = x[i*dim+1] + pimages[j*dim+1];        
            T xj2 = x[i*dim+2] + pimages[j*dim+2];        
            T xc0 = B2C[0]*xj0 + B2C[3]*xj1 + B2C[6]*xj2;  // map it to the unit square            
            T xc1 = B2C[1]*xj0 + B2C[4]*xj1 + B2C[7]*xj2;        
            T xc2 = B2C[2]*xj0 + B2C[5]*xj1 + B2C[8]*xj2;   
            /// check if the mapped point is inside the bounding box
            if ((wc[0] <= xc0) && (xc0 <= wc[18]) && (wc[1] <= xc1) && (xc1 <= wc[19]) && (wc[2] <= xc2) && (xc2 <= wc[20])) 
                glistnum[i] += 1;            // increase the number of ghost atoms by 1                      
        }        
        i += blockDim.x * gridDim.x;
    }
}

template <typename T> void gpuGhostAtoms(int *glistnum, T *x, T *pimages, T *wc, T *B2C, int n, int m, int dim)
{        
    int blockDim = 256;
    int gridDim = (n + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;
    if (dim==2)
        gpuKernelGhostAtoms2D<<<gridDim, blockDim>>>(glistnum, x, pimages, wc, B2C, n, m, dim);
    else
        gpuKernelGhostAtoms3D<<<gridDim, blockDim>>>(glistnum, x, pimages, wc, B2C, n, m, dim);
}
template void gpuGhostAtoms(int*, double*, double*, double*, double*, int, int, int);
template void gpuGhostAtoms(int*, float*, float*, float*, float*, int, int, int);


template <typename T>
__global__ void gpuKernelCreateAtomList2D(int *ilist, int *glistnumsum, int *atomtype, T *x, T *pimages, T *wc, T *B2C, int n, int m, int dim)
{
    int i = threadIdx.x + blockIdx.x * blockDim.x;
    while (i < n) {
        ilist[i] = i;         // add atom i to the list
        int q = n + glistnumsum[i]; // offset the starting position by n
        int k = 0;            
        for (int j=1; j<m; j++) { // loop over each periodic image of atom i
            T xj0 = x[i*dim+0] + pimages[j*dim+0];  // periodic image of x      
            T xj1 = x[i*dim+1] + pimages[j*dim+1];        
            T xc0 = B2C[0]*xj0 + B2C[2]*xj1;        // map it to the unit square      
            T xc1 = B2C[1]*xj0 + B2C[3]*xj1;        
            /// check if the mapped point is inside the bounding box
            if ((wc[0] <= xc0) && (xc0 <= wc[4]) &&  (wc[1] <= xc1) && (xc1 <= wc[5])) {                    
                x[dim*(q+k)+0] = xj0; // add the periodic image as a ghost atom
                x[dim*(q+k)+1] = xj1; //
                atomtype[q+k] = atomtype[i];
                ilist[q+k] = i;       // add atom i to the list
                k += 1;            // increase the number of ghost atoms for atom i by 1    
            }
        }        
        i += blockDim.x * gridDim.x;
    }
}

template <typename T>
__global__ void gpuKernelCreateAtomList3D(int *ilist, int *glistnumsum, int *atomtype, T *x, T *pimages, T *wc, T *B2C, int n, int m, int dim)
{
    int i = threadIdx.x + blockIdx.x * blockDim.x;
    while (i < n) {
        ilist[i] = i; // add atom i to the list
        int q = n + glistnumsum[i]; // offset the starting position by n
        int k = 0;            
        for (int j=1; j<m; j++) { // loop over each periodic image of atom i
            T xj0 = x[i*dim+0] + pimages[j*dim+0];    // periodic image of x          
            T xj1 = x[i*dim+1] + pimages[j*dim+1];        
            T xj2 = x[i*dim+2] + pimages[j*dim+2];        
            T xc0 = B2C[0]*xj0 + B2C[3]*xj1 + B2C[6]*xj2;  // map it to the unit square            
            T xc1 = B2C[1]*xj0 + B2C[4]*xj1 + B2C[7]*xj2;        
            T xc2 = B2C[2]*xj0 + B2C[5]*xj1 + B2C[8]*xj2;   
            /// check if the mapped point is inside the bounding box
            if ((wc[0] <= xc0) && (xc0 <= wc[18]) && (wc[1] <= xc1) && (xc1 <= wc[19]) && (wc[2] <= xc2) && (xc2 <= wc[20])) {
                x[dim*(q+k)+0] = xj0; // add the periodic image as a ghost atom
                x[dim*(q+k)+1] = xj1;  
                x[dim*(q+k)+2] = xj2;
                atomtype[q+k] = atomtype[i];
                ilist[q+k] = i;      // add atom i to the list
                k += 1;           // increase the number of ghost atoms by 1              
            }
        }       
        i += blockDim.x * gridDim.x;
    }
}

template <typename T> void gpuCreateAtomList(int *ilist, int *glistnumsum, int *glistnum, int *atomtype,
        int *d_sums, int *d_incr, T *x, T *pimages, T *wc, T *B2C, int n, int m, int dim)
{                        
    // a list contains the starting position of the ghost atom of every atom i
    gpuCumsum(glistnumsum, glistnum, d_sums, d_incr, n+1); 

    int blockDim = 256;
    int gridDim = (n + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;
    if (dim==2)
        gpuKernelCreateAtomList2D<<<gridDim, blockDim>>>(ilist, glistnumsum, atomtype, x, pimages, wc, B2C, n, m, dim);
    else
        gpuKernelCreateAtomList3D<<<gridDim, blockDim>>>(ilist, glistnumsum, atomtype, x, pimages, wc, B2C, n, m, dim);
}
template void gpuCreateAtomList(int*, int*, int*, int*, int*, int*, double*, double*, double*, double*, int, int, int);
template void gpuCreateAtomList(int*, int*, int*, int*, int*, int*, float*, float*, float*, float*, int, int, int);


template <typename T>
__global__ void gpuKernelCellList2D(int *clist, int *c2anum, T *xi, T *eta1, T *eta2, T *eta3, T *B2C, int *nc, int n, int dim)
{
    int i = threadIdx.x + blockIdx.x * blockDim.x;
    while (i < n) {
        int j1, j2;
        int k = dim*i;
        // position of atom i 
        T xt0 = B2C[0]*xi[k] + B2C[2]*xi[k+1];        
        T xt1 = B2C[1]*xi[k] + B2C[3]*xi[k+1];     
        // identify a cell containing atom i 
        for (j1=0; j1<nc[0]; j1++)
            if ((eta1[j1] <= xt0) && (xt0<= eta1[j1+1]))
                break;
        for (j2=0; j2<nc[1]; j2++) 
            if ((eta2[j2] <= xt1) && (xt1<= eta2[j2+1]))
                break;
        int c = j1 + nc[0]*j2;
        clist[i] = c; // link that cell to atom i
        // use atomicAdd on GPU to avoid race condition                    
        atomicAdd(&c2anum[c], 1); // increase the number of atoms in that cell by 1         
        i += blockDim.x * gridDim.x;
    }
}

template <typename T>
__global__ void gpuKernelCellList3D(int *clist, int *c2anum, T *xi, T *eta1, T *eta2, T *eta3, T *B2C, int *nc, int n, int dim)
{
    int i = threadIdx.x + blockIdx.x * blockDim.x;
    while (i < n) {
        int j1, j2, j3;
        int k = dim*i;
        // position of atom i in the unit square/cube
        T xt0 = B2C[0]*xi[k+0] + B2C[3]*xi[k+1] + B2C[6]*xi[k+2];        
        T xt1 = B2C[1]*xi[k+0] + B2C[4]*xi[k+1] + B2C[7]*xi[k+2];                                    
        T xt2 = B2C[2]*xi[k+0] + B2C[5]*xi[k+1] + B2C[8]*xi[k+2];        
        // identify a cell containing atom i 
        for (j1=0; j1<nc[0]; j1++)
            if ((eta1[j1] <= xt0) && (xt0<= eta1[j1+1]))
                break;
        for (j2=0; j2<nc[1]; j2++) 
            if ((eta2[j2] <= xt1) && (xt1<= eta2[j2+1]))
                break;
        for (j3=0; j3<nc[2]; j3++) 
            if ((eta3[j3] <= xt2) && (xt2<= eta3[j3+1]))
                break;
        int c = j1 + nc[0]*j2 + nc[0]*nc[1]*j3; // cell c
        clist[i] = c; // link that cell to atom i
        // use atomicAdd on GPU to avoid race condition                    
        atomicAdd(&c2anum[c], 1); // increase the number of atoms in that cell by 1                                 
        i += blockDim.x * gridDim.x;
    }
}


template <typename T> void gpuCellList(int *clist, int *c2anum, T *xi, T *eta1, T *eta2, T *eta3, T *B2C, int *nc, int inum, int gnum, int dim)
{                        
    int n = inum + gnum;

    int blockDim = 256;
    int gridDim = (n + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;
    if (dim==2) {
        int k = nc[0]*nc[1];   // the total number of cells        
        gpuArraySetValue(c2anum, 0, k);        
        gpuKernelCellList2D<<<gridDim, blockDim>>>(clist, c2anum, xi, eta1, eta2, eta3, B2C, nc, n, dim);
    }
    else {
        int k = nc[0]*nc[1]*nc[2];   // the total number of cells        
        gpuArraySetValue(c2anum, 0, k);        
        gpuKernelCellList3D<<<gridDim, blockDim>>>(clist, c2anum, xi, eta1, eta2, eta3, B2C, nc, n, dim);
    }
}
template void gpuCellList(int*, int*, double*, double*, double*, double*, double*, int*, int, int, int);
template void gpuCellList(int*, int*, float*, float*, float*, float*, float*, int*, int, int, int);

__global__ void gpuKernelCell2AtomList(int *c2alist, int *c2anumsum, int *c2anum, int *clist, int n)
{
    int i = threadIdx.x + blockIdx.x * blockDim.x;
    while (i < n) {
        int k = clist[i];         // cell k contains atom i           
        int i1 = c2anumsum[k];    // starting position for cell k
        int i2 = c2anum[k];       // position of the atom i in cell k
        c2alist[i1+i2] = i;       // add atom i to the list of cell k
        // (using atomicAdd on GPU to avoid race condition)                
        atomicAdd(&c2anum[k], 1);           // increase the number of atoms in cell k by 1 
        i += blockDim.x * gridDim.x;
    }
}

void gpuCell2AtomList(int *c2alist, int *c2anumsum, int *c2anum, int *clist, int *d_sums, int *d_incr, int *nc, int inum, int gnum, int dim)
{   
    // number of cells
    int ncell = (dim==2) ? nc[0]*nc[1] : nc[0]*nc[1]*nc[2];
    
    // number of atoms
    int natom = inum+gnum;        

    // a list contains the starting positions of the first atom of a cell
    gpuCumsum(c2anumsum, c2anum, d_sums, d_incr, ncell+1); 
    
    // initialize the list containting the number of atoms in each cell
    gpuArraySetValue(c2anum, 0, ncell);        

    int blockDim = 256;
    int gridDim = (natom + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;
    
    gpuKernelCell2AtomList<<<gridDim, blockDim>>>(c2alist, c2anumsum, c2anum, clist, natom);
}

template <typename T>
__global__ void gpuKernelVerletAtoms2D(int *verletnum, T *x, T *ellipsoid, int *atomtype, int *ilist, 
        int *clist, int *c2alist, int *c2anum, int *c2anumsum, int *nc, int ntype, int inum, int dim, int dimsq)
{
    int ii = threadIdx.x + blockIdx.x * blockDim.x;
    while (ii < inum) {    
        int i = ilist[ii];
        verletnum[i] = 0;
        T xi0 = x[i*dim];        // position of atom i
        T xi1 = x[i*dim+1];      // position of atom i
        int ti = atomtype[i];      // element of atom i 
        int j = clist[i];         // cell j contains atom i           
        int j1 = j%nc[0];
        int j2 = (j-j1)/nc[0];
        for (int i1=-1; i1<=1; i1++) {
            int k1 = j1 + i1;
            for (int i2=-1; i2<=1; i2++) {
                int k2 = j2 + i2;
                int k = k1 + nc[0]*k2; // cell k
                int m = c2anum[k];     // number of atoms in cell k
                int s = c2anumsum[k];  // starting position of the first atom in cell k
                for (int l=0; l<m ; l++) {
                    j = c2alist[s+l];  // atom j
                    int tj = atomtype[j];      // element of atom j
                    int tij = dimsq*(ti + tj*ntype);
                    T A00 = ellipsoid[tij]; // ellipsoid for pair (i,j)   
                    T A10 = ellipsoid[tij+1]; // ellipsoid for pair (i,j)
                    T A01 = ellipsoid[tij+2]; // ellipsoid for pair (i,j)  
                    T A11 = ellipsoid[tij+3]; // ellipsoid for pair (i,j)                                                      
                    T xij0 = x[j*dim] - xi0;  // xj - xi
                    T xij1 = x[j*dim+1] - xi1; // xj - xi
                    // distance between atom i and atom j 
                    T rij = xij0*(A00*xij0 + A01*xij1) + xij1*(A10*xij0 + A11*xij1);                        
                    if (rij <= 1.0)
                        verletnum[i] += 1;
                }
            }
        }                
        ii += blockDim.x * gridDim.x;
    }    
}

template <typename T>
__global__ void gpuKernelVerletAtoms3D(int *verletnum, T *x, T *ellipsoid, int *atomtype, int *ilist, 
        int *clist, int *c2alist, int *c2anum, int *c2anumsum, int *nc, int ntype, int inum, int dim, int dimsq)
{
    int ii = threadIdx.x + blockIdx.x * blockDim.x;
    while (ii < inum) {    
        int i = ilist[ii];
        verletnum[i] = 0;
        T xi0 = x[i*dim];        // position of atom i
        T xi1 = x[i*dim+1];      // position of atom i
        T xi2 = x[i*dim+2];      // position of atom i
        int ti = atomtype[i];      // element of atom i 
        int j = clist[i];         // cell j contains atom i       
        int n = j%(nc[0]*nc[1]);
        int j3 = (j-n)/(nc[0]*nc[1]);            
        int j1 = n%nc[0];
        int j2 = (j-j1)/nc[0];
        for (int i1=-1; i1<=1; i1++) {
            int k1 = j1 + i1;
            for (int i2=-1; i2<=1; i2++) {
                int k2 = j2 + i2;
                for (int i3=-1; i3<=1; i3++) {
                    int k3 = j3 + i3;                    
                    int k = k1 + nc[0]*k2 + nc[0]*nc[1]*k3; // cell k
                    int m = c2anum[k];     // number of atoms in cell k
                    int s = c2anumsum[k];  // starting position of the first atom in cell k
                    for (int l=0; l<m ; l++) {
                        j = c2alist[s+l];  // atom j
                        int tj = atomtype[j];      // element of atom j
                        int tij = dimsq*(ti + tj*ntype);
                        T A00 = ellipsoid[tij]; // ellipsoid for element t   
                        T A10 = ellipsoid[tij+1]; // ellipsoid for element t   
                        T A20 = ellipsoid[tij+2]; // ellipsoid for element t   
                        T A01 = ellipsoid[tij+3]; // ellipsoid for element t   
                        T A11 = ellipsoid[tij+4]; // ellipsoid for element t   
                        T A21 = ellipsoid[tij+5]; // ellipsoid for element t   
                        T A02 = ellipsoid[tij+6]; // ellipsoid for element t   
                        T A12 = ellipsoid[tij+7]; // ellipsoid for element t   
                        T A22 = ellipsoid[tij+8]; // ellipsoid for element t                                                               
                        T xij0 = x[j*dim] - xi0;  // xj - xi
                        T xij1 = x[j*dim+1] - xi1; // xj - xi
                        T xij2 = x[j*dim+2] - xi2; // xj - xi
                        // distance between atom i and atom j 
                        T rij = xij0*(A00*xij0 + A01*xij1 + A02*xij2) + xij1*(A10*xij0 + A11*xij1 + A12*xij2) + xij2*(A20*xij0 + A21*xij1 + A22*xij2);
                        if (rij <= 1.0)
                            verletnum[i] += 1;
                    }
                }
            }
        }                        
        ii += blockDim.x * gridDim.x;
    }    
}

template <typename T> void gpuVerletAtoms(int *verletnum, T *x, T *ellipsoid, int *atomtype, int *ilist, 
        int *clist, int *c2alist, int *c2anum, int *c2anumsum, int *nc, int ntype, int inum, int dim)
{                        
    int n = inum;
    int dimsq = dim*dim;
    int blockDim = 256;
    int gridDim = (n + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;
    if (dim==2) {
        gpuKernelVerletAtoms2D<<<gridDim, blockDim>>>(verletnum, x, ellipsoid, atomtype, ilist, clist, c2alist, 
                c2anum, c2anumsum, nc, ntype, n, dim, dimsq);
    }
    else {
        gpuKernelVerletAtoms3D<<<gridDim, blockDim>>>(verletnum, x, ellipsoid, atomtype, ilist, clist, c2alist, 
                c2anum, c2anumsum, nc, ntype, n, dim, dimsq);
    }
}
template void gpuVerletAtoms(int*, double*, double*, int*, int*, int*, int*, int*, int*, int*, int, int, int);
template void gpuVerletAtoms(int*, float*, float*, int*, int*, int*, int*, int*, int*, int*, int, int, int);

template <typename T>
__global__ void gpuKernelCreateVerletList2D(int *verletlist, T *x, T *ellipsoid, int *verletnum, int *verletnumsum, int *atomtype, 
      int *ilist, int *clist, int *c2alist, int *c2anum, int *c2anumsum, int *nc, int ntype, int inum, int dim, int dimsq)
{
    int ii = threadIdx.x + blockIdx.x * blockDim.x;
    while (ii < inum) {   
        int i = ilist[ii];
        T xi0 = x[i*dim];        // position of atom i
        T xi1 = x[i*dim+1];      // position of atom i
        int ti = atomtype[i];      // element of atom i 
        int nstart = verletnumsum[i];     // starting 
        int ninc = 0;                // increment
        int j = clist[i];         // cell j contains atom i           
        int j1 = j%nc[0];
        int j2 = (j-j1)/nc[0];
        for (int i1=-1; i1<=1; i1++) {
            int k1 = j1 + i1;
            for (int i2=-1; i2<=1; i2++) {
                int k2 = j2 + i2;
                int k = k1 + nc[0]*k2; // cell k
                int m = c2anum[k];     // number of atoms in cell k
                int s = c2anumsum[k];  // starting position of the first atom in cell k
                for (int l=0; l<m ; l++) { //loop over each atom j in cell k
                    j = c2alist[s+l];  // atom j
                    int tj = atomtype[j];      // element of atom j
                    int tij = dimsq*(ti + tj*ntype);
                    T A00 = ellipsoid[tij]; // ellipsoid for pair (i,j)   
                    T A10 = ellipsoid[tij+1]; // ellipsoid for pair (i,j)
                    T A01 = ellipsoid[tij+2]; // ellipsoid for pair (i,j)  
                    T A11 = ellipsoid[tij+3]; // ellipsoid for pair (i,j)                                                      
                    T xij0 = x[j*dim] - xi0;  // xj - xi
                    T xij1 = x[j*dim+1] - xi1; // xj - xi
                    // distance between atom i and atom j 
                    T rij = xij0*(A00*xij0 + A01*xij1) + xij1*(A10*xij0 + A11*xij1);                        
                    if (rij <= 1.0) {
                        verletlist[nstart + ninc] = j; // add atom j into the list
                        ninc += 1;
                    }
                }
            }
        }               
        verletnum[i] = ninc;
        ii += blockDim.x * gridDim.x;
    }
}

template <typename T>
__global__ void gpuKernelCreateVerletList3D(int *verletlist, T *x, T *ellipsoid, int *verletnum, int *verletnumsum, int *atomtype, 
      int *ilist, int *clist, int *c2alist, int *c2anum, int *c2anumsum, int *nc, int ntype, int inum, int dim, int dimsq)
{
    int ii = threadIdx.x + blockIdx.x * blockDim.x;
    while (ii < inum) {    
        int i = ilist[ii];
        T xi0 = x[i*dim];        // position of atom i
        T xi1 = x[i*dim+1];      // position of atom i
        T xi2 = x[i*dim+2];      // position of atom i
        int ti = atomtype[i];      // element of atom i         
        int nstart = verletnumsum[i];     // starting 
        int ninc = 0;                // increment            
        int j = clist[i];         // cell j contains atom i       
        int n = j%(nc[0]*nc[1]);
        int j3 = (j-n)/(nc[0]*nc[1]);            
        int j1 = n%nc[0];
        int j2 = (j-j1)/nc[0];
        for (int i1=-1; i1<=1; i1++) {
            int k1 = j1 + i1;
            for (int i2=-1; i2<=1; i2++) {
                int k2 = j2 + i2;
                for (int i3=-1; i3<=1; i3++) {
                    int k3 = j3 + i3;                    
                    int k = k1 + nc[0]*k2 + nc[0]*nc[1]*k3; // cell k
                    int m = c2anum[k];     // number of atoms in cell k
                    int s = c2anumsum[k];  // starting position of the first atom in cell k
                    for (int l=0; l<m ; l++) {
                        j = c2alist[s+l];  // atom j
                        int tj = atomtype[j];      // element of atom j
                        int tij = dimsq*(ti + tj*ntype);
                        T A00 = ellipsoid[tij]; // ellipsoid for element t   
                        T A10 = ellipsoid[tij+1]; // ellipsoid for element t   
                        T A20 = ellipsoid[tij+2]; // ellipsoid for element t   
                        T A01 = ellipsoid[tij+3]; // ellipsoid for element t   
                        T A11 = ellipsoid[tij+4]; // ellipsoid for element t   
                        T A21 = ellipsoid[tij+5]; // ellipsoid for element t   
                        T A02 = ellipsoid[tij+6]; // ellipsoid for element t   
                        T A12 = ellipsoid[tij+7]; // ellipsoid for element t   
                        T A22 = ellipsoid[tij+8]; // ellipsoid for element t           
                        T xij0 = x[j*dim] - xi0;  // xj - xi
                        T xij1 = x[j*dim+1] - xi1; // xj - xi
                        T xij2 = x[j*dim+2] - xi2; // xj - xi
                        // distance between atom i and atom j 
                        T rij = xij0*(A00*xij0 + A01*xij1 + A02*xij2) + xij1*(A10*xij0 + A11*xij1 + A12*xij2) + xij2*(A20*xij0 + A21*xij1 + A22*xij2);
                        if (rij <= 1.0) {
                            verletlist[nstart + ninc] = j; // add atom j into the list
                            ninc += 1;
                        }                            
                    }
                }
            }
        }                
        verletnum[i] = ninc;
        ii += blockDim.x * gridDim.x;
    }
}

template <typename T> void gpuCreateVerletList(int *verletlist, T *x, T *ellipsoid, int *verletnum,  int *verletnumsum, int *tm1, int *tm2,
     int *atomtype, int *ilist, int *clist, int *c2alist, int *c2anum, int *c2anumsum, int *nc, int ntype, int inum, int dim)
{                        
    int n = inum;
    int blockDim = 256;
    int gridDim = (n + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;

    int dimsq = dim*dim;

    // a list contains the starting positions of the first j atom 
    gpuCumsum(verletnumsum, verletnum, tm1, tm2, inum+1); 

    if (dim==2) {
        gpuKernelCreateVerletList2D<<<gridDim, blockDim>>>(verletlist, x, ellipsoid, verletnum, verletnumsum, atomtype, 
                ilist, clist, c2alist, c2anum, c2anumsum, nc, ntype, n, dim, dimsq);
    }
    else {
        gpuKernelCreateVerletList3D<<<gridDim, blockDim>>>(verletlist, x, ellipsoid, verletnum, verletnumsum, atomtype, 
                ilist, clist, c2alist, c2anum, c2anumsum, nc, ntype, n, dim, dimsq);
    }
}
template void gpuCreateVerletList(int*, double*, double*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int, int, int);
template void gpuCreateVerletList(int*, float*, float*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int, int, int);

template <typename T> __global__ void gpuKernelFullNeighNum2D(int *neighnum, T *x, T* ellipsoid, int *atomtype, 
        int *ilist, int *verletlist, int *verletnum, int *verletnumsum, int ntype, int inum, int dim, int dimsq)
{
    int ii = threadIdx.x + blockIdx.x * blockDim.x;
    while (ii < inum) {            
        int i = ilist[ii];  // atom i
        T xi0 = x[i*dim];        // position of atom i
        T xi1 = x[i*dim+1];      // position of atom i
        int ti = atomtype[i];      // element of atom i 
        int m = verletnum[i]; // number of atoms around i 
        int start = verletnumsum[i];     // starting 
        int ninc = 0;                // increment 
        for (int l=0; l<m ; l++) {   // loop over each atom around atom i
            int j = verletlist[start+l];  // atom j    
            int tj = atomtype[j];      // element of atom j
            int tij = dimsq*(ti + tj*ntype);
            T A00 = ellipsoid[tij]; // ellipsoid for pair (i,j)   
            T A10 = ellipsoid[tij+1]; // ellipsoid for pair (i,j)
            T A01 = ellipsoid[tij+2]; // ellipsoid for pair (i,j)  
            T A11 = ellipsoid[tij+3]; // ellipsoid for pair (i,j)                                  
            T xij0 = x[j*dim] - xi0;  // xj - xi
            T xij1 = x[j*dim+1] - xi1; // xj - xi
            // distance between atom i and atom j 
            T rij = xij0*(A00*xij0 + A01*xij1) + xij1*(A10*xij0 + A11*xij1);                        
            if (rij <= 1.0) // if atom j is inside the ellipsoid 
                ninc += 1;  // increase the number of neighbors by 1               
        }
        neighnum[i] = ninc; // number of neighbors of atom i

        ii += blockDim.x * gridDim.x;                
    }
}

template <typename T> __global__ void gpuKernelFullNeighNum3D(int *neighnum, T *x, T* ellipsoid, int *atomtype, 
        int *ilist, int *verletlist, int *verletnum, int *verletnumsum, int ntype, int inum, int dim, int dimsq)
{
    int ii = threadIdx.x + blockIdx.x * blockDim.x;
    while (ii < inum) {            
        int i = ilist[ii];  // atom i
        T xi0 = x[i*dim];        // position of atom i
        T xi1 = x[i*dim+1];      // position of atom i
        T xi2 = x[i*dim+2];      // position of atom i
        int ti = atomtype[i];      // element of atom i           
        int m = verletnum[i]; // number of atoms around i 
        int start = verletnumsum[i];     // starting 
        int ninc = 0;                // increment 
        for (int l=0; l<m ; l++) {   // loop over each atom around atom i
            int j = verletlist[start+l];  // atom j  
            int tj = atomtype[j];      // element of atom j
            int tij = dimsq*(ti + tj*ntype);
            T A00 = ellipsoid[tij]; // ellipsoid for element t   
            T A10 = ellipsoid[tij+1]; // ellipsoid for element t   
            T A20 = ellipsoid[tij+2]; // ellipsoid for element t   
            T A01 = ellipsoid[tij+3]; // ellipsoid for element t   
            T A11 = ellipsoid[tij+4]; // ellipsoid for element t   
            T A21 = ellipsoid[tij+5]; // ellipsoid for element t   
            T A02 = ellipsoid[tij+6]; // ellipsoid for element t   
            T A12 = ellipsoid[tij+7]; // ellipsoid for element t   
            T A22 = ellipsoid[tij+8]; // ellipsoid for element t                                                                           
            T xij0 = x[j*dim] - xi0;  // xj - xi
            T xij1 = x[j*dim+1] - xi1; // xj - xi
            T xij2 = x[j*dim+2] - xi2; // xj - xi
            // distance between atom i and atom j 
            T rij = xij0*(A00*xij0 + A01*xij1 + A02*xij2) + xij1*(A10*xij0 + A11*xij1 + A12*xij2) + xij2*(A20*xij0 + A21*xij1 + A22*xij2);
            if (rij <= 1.0) // if atom j is inside the ellipsoid 
                ninc += 1;  // increase the number of neighbors by 1                               
        }            
        neighnum[i] = ninc; // number of neighbors of atom i        

        ii += blockDim.x * gridDim.x;                
    }
}

template <typename T> void gpuFullNeighNum(int *neighnum, T *x, T* ellipsoid, int *atomtype, 
        int *ilist, int *verletlist, int *verletnum, int *verletnumsum, int ntype, int inum, int dim)
{                        
    int n = inum;
    int blockDim = 256;
    int gridDim = (n + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;

    int dimsq = dim*dim;
    if (dim==2) {
        gpuKernelFullNeighNum2D<<<gridDim, blockDim>>>(neighnum, x, ellipsoid, atomtype, 
                ilist, verletlist, verletnum, verletnumsum, ntype, n, dim, dimsq);
    }
    else {
        gpuKernelFullNeighNum3D<<<gridDim, blockDim>>>(neighnum, x, ellipsoid, atomtype, 
                ilist, verletlist, verletnum, verletnumsum, ntype, n, dim, dimsq);
    }
}
template void gpuFullNeighNum(int*, double*, double*, int*, int*, int*, int*, int*, int, int, int);
template void gpuFullNeighNum(int*, float*, float*, int*, int*, int*, int*, int*, int, int, int);

template <typename T> __global__ void gpuKernelFullNeighList2D(int *neighlist, T *x, T* ellipsoid, int *neighnum, 
   int *neighnumsum, int *atomtype, int *ilist, int *verletlist, int *verletnum, int *verletnumsum, int ntype, int inum, int dim, int dimsq)
{
    int ii = threadIdx.x + blockIdx.x * blockDim.x;
    while (ii < inum) {            
        int i = ilist[ii];  // atom i
        T xi0 = x[i*dim];        // position of atom i
        T xi1 = x[i*dim+1];      // position of atom i
        int ti = atomtype[i];      // element of atom i         
        int m = verletnum[i]; // number of atoms around i 
        int jstart = verletnumsum[i];     // starting 
        int nstart = neighnumsum[i];   
        int ninc = 0;                // increment 
        for (int l=0; l<m ; l++) {   // loop over each atom around atom i
            int j = verletlist[jstart+l];  // atom j              
            int tj = atomtype[j];      // element of atom j
            int tij = dimsq*(ti + tj*ntype);
            T A00 = ellipsoid[tij]; // ellipsoid for pair (i,j)   
            T A10 = ellipsoid[tij+1]; // ellipsoid for pair (i,j)
            T A01 = ellipsoid[tij+2]; // ellipsoid for pair (i,j)  
            T A11 = ellipsoid[tij+3]; // ellipsoid for pair (i,j)                                                                                      
            T xij0 = x[j*dim] - xi0;  // xj - xi
            T xij1 = x[j*dim+1] - xi1; // xj - xi
            // distance between atom i and atom j 
            T rij = xij0*(A00*xij0 + A01*xij1) + xij1*(A10*xij0 + A11*xij1);                        
            if (rij <= 1.0) { // if atom j is inside the ellipsoid 
                neighlist[nstart+ninc] = j;
                ninc += 1;  // increase the number of neighbors by 1               
            }
        }
        neighnum[i] = ninc; // number of neighbors of atom i

        ii += blockDim.x * gridDim.x;                
    }
}

template <typename T> __global__ void gpuKernelFullNeighList3D(int *neighlist, T *x, T* ellipsoid, int *neighnum, 
   int *neighnumsum, int *atomtype, int *ilist, int *verletlist, int *verletnum, int *verletnumsum, int ntype, int inum, int dim, int dimsq)
{
    int ii = threadIdx.x + blockIdx.x * blockDim.x;
    while (ii < inum) {            
        int i = ilist[ii];  // atom i
        T xi0 = x[i*dim];        // position of atom i
        T xi1 = x[i*dim+1];      // position of atom i
        T xi2 = x[i*dim+2];      // position of atom i
        int ti = atomtype[i];      // element of atom i 
        int m = verletnum[i]; // number of atoms around i 
        int jstart = verletnumsum[i];     // starting 
        int nstart = neighnumsum[i];   
        int ninc = 0;                // increment 
        for (int l=0; l<m ; l++) {   // loop over each atom around atom i
            int j = verletlist[jstart+l];  // atom j          
            int tj = atomtype[j];      // element of atom j
            int tij = dimsq*(ti + tj*ntype);
            T A00 = ellipsoid[tij]; // ellipsoid for element t   
            T A10 = ellipsoid[tij+1]; // ellipsoid for element t   
            T A20 = ellipsoid[tij+2]; // ellipsoid for element t   
            T A01 = ellipsoid[tij+3]; // ellipsoid for element t   
            T A11 = ellipsoid[tij+4]; // ellipsoid for element t   
            T A21 = ellipsoid[tij+5]; // ellipsoid for element t   
            T A02 = ellipsoid[tij+6]; // ellipsoid for element t   
            T A12 = ellipsoid[tij+7]; // ellipsoid for element t   
            T A22 = ellipsoid[tij+8]; // ellipsoid for element t                                                                                               
            T xij0 = x[j*dim] - xi0;  // xj - xi
            T xij1 = x[j*dim+1] - xi1; // xj - xi
            T xij2 = x[j*dim+2] - xi2; // xj - xi
            // distance between atom i and atom j 
            T rij = xij0*(A00*xij0 + A01*xij1 + A02*xij2) + xij1*(A10*xij0 + A11*xij1 + A12*xij2) + xij2*(A20*xij0 + A21*xij1 + A22*xij2);
            if (rij <= 1.0) { // if atom j is inside the ellipsoid 
                neighlist[nstart+ninc] = j;
                ninc += 1;  // increase the number of neighbors by 1               
            }
        }            
        neighnum[i] = ninc; // number of neighbors of atom i

        ii += blockDim.x * gridDim.x;                
    }
}

template <typename T> void gpuFullNeighList(int *neighlist, T *x, T* ellipsoid, int *neighnum, int *neighnumsum, 
      int *atomtype, int *ilist, int *verletlist, int *verletnum, int *verletnumsum, int *tm1, int *tm2, int ntype, int inum, int dim)
{                        
    int n = inum;
    int blockDim = 256;
    int gridDim = (n + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;

    // a list contains the starting positions of the first neighbor 
    gpuCumsum(neighnumsum, neighnum, tm1, tm2, inum+1); 

    int dimsq = dim*dim;
    if (dim==2) {
        gpuKernelFullNeighList2D<<<gridDim, blockDim>>>(neighlist, x, ellipsoid, neighnum, neighnumsum, atomtype, 
                ilist, verletlist, verletnum, verletnumsum, ntype, n, dim, dimsq);
    }
    else {
        gpuKernelFullNeighList3D<<<gridDim, blockDim>>>(neighlist, x, ellipsoid, neighnum, neighnumsum, atomtype, 
                ilist, verletlist, verletnum, verletnumsum, ntype, n, dim, dimsq);
    }
}
template void gpuFullNeighList(int*, double*, double*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int, int, int);
template void gpuFullNeighList(int*, float*, float*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int, int, int);


template <typename T> __global__ void gpuKernelHalfNeighNum2D(int *neighnum, T *x, T* ellipsoid, int *atomtype, 
        int *ilist, int *verletlist, int *verletnum, int *verletnumsum, int ntype, int inum, int dim, int dimsq)
{
    int ii = threadIdx.x + blockIdx.x * blockDim.x;
    while (ii < inum) {            
        int i = ilist[ii];  // atom i
        T xi0 = x[i*dim];        // position of atom i
        T xi1 = x[i*dim+1];      // position of atom i
        int ti = atomtype[i];      // element of atom i 
        int m = verletnum[i]; // number of atoms around i 
        int start = verletnumsum[i];     // starting 
        int ninc = 0;                // increment 
        for (int l=0; l<m ; l++) {   // loop over each atom around atom i
            int j = verletlist[start+l];  // atom j    
            if (j >= i) {        
                int tj = atomtype[j];      // element of atom j
                int tij = dimsq*(ti + tj*ntype);
                T A00 = ellipsoid[tij]; // ellipsoid for pair (i,j)   
                T A10 = ellipsoid[tij+1]; // ellipsoid for pair (i,j)
                T A01 = ellipsoid[tij+2]; // ellipsoid for pair (i,j)  
                T A11 = ellipsoid[tij+3]; // ellipsoid for pair (i,j)                                  
                T xij0 = x[j*dim] - xi0;  // xj - xi
                T xij1 = x[j*dim+1] - xi1; // xj - xi
                // distance between atom i and atom j 
                T rij = xij0*(A00*xij0 + A01*xij1) + xij1*(A10*xij0 + A11*xij1);                        
                if (rij <= 1.0) // if atom j is inside the ellipsoid 
                    ninc += 1;  // increase the number of neighbors by 1    
            }
        }
        neighnum[i] = ninc; // number of neighbors of atom i

        ii += blockDim.x * gridDim.x;                
    }
}

template <typename T> __global__ void gpuKernelHalfNeighNum3D(int *neighnum, T *x, T* ellipsoid, int *atomtype, 
        int *ilist, int *verletlist, int *verletnum, int *verletnumsum, int ntype, int inum, int dim, int dimsq)
{
    int ii = threadIdx.x + blockIdx.x * blockDim.x;
    while (ii < inum) {            
        int i = ilist[ii];  // atom i
        T xi0 = x[i*dim];        // position of atom i
        T xi1 = x[i*dim+1];      // position of atom i
        T xi2 = x[i*dim+2];      // position of atom i
        int ti = atomtype[i];      // element of atom i           
        int m = verletnum[i]; // number of atoms around i 
        int start = verletnumsum[i];     // starting 
        int ninc = 0;                // increment 
        for (int l=0; l<m ; l++) {   // loop over each atom around atom i
            int j = verletlist[start+l];  // atom j  
            if (j >= i) {        
                int tj = atomtype[j];      // element of atom j
                int tij = dimsq*(ti + tj*ntype);
                T A00 = ellipsoid[tij]; // ellipsoid for element t   
                T A10 = ellipsoid[tij+1]; // ellipsoid for element t   
                T A20 = ellipsoid[tij+2]; // ellipsoid for element t   
                T A01 = ellipsoid[tij+3]; // ellipsoid for element t   
                T A11 = ellipsoid[tij+4]; // ellipsoid for element t   
                T A21 = ellipsoid[tij+5]; // ellipsoid for element t   
                T A02 = ellipsoid[tij+6]; // ellipsoid for element t   
                T A12 = ellipsoid[tij+7]; // ellipsoid for element t   
                T A22 = ellipsoid[tij+8]; // ellipsoid for element t                                                                           
                T xij0 = x[j*dim] - xi0;  // xj - xi
                T xij1 = x[j*dim+1] - xi1; // xj - xi
                T xij2 = x[j*dim+2] - xi2; // xj - xi
                // distance between atom i and atom j 
                T rij = xij0*(A00*xij0 + A01*xij1 + A02*xij2) + xij1*(A10*xij0 + A11*xij1 + A12*xij2) + xij2*(A20*xij0 + A21*xij1 + A22*xij2);
                if (rij <= 1.0) // if atom j is inside the ellipsoid 
                    ninc += 1;  // increase the number of neighbors by 1      
            }
        }            
        neighnum[i] = ninc; // number of neighbors of atom i        

        ii += blockDim.x * gridDim.x;                
    }
}

template <typename T> void gpuHalfNeighNum(int *neighnum, T *x, T* ellipsoid, int *atomtype, 
        int *ilist, int *verletlist, int *verletnum, int *verletnumsum, int ntype, int inum, int dim)
{                        
    int n = inum;
    int blockDim = 256;
    int gridDim = (n + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;

    int dimsq = dim*dim;
    if (dim==2) {
        gpuKernelHalfNeighNum2D<<<gridDim, blockDim>>>(neighnum, x, ellipsoid, atomtype, 
                ilist, verletlist, verletnum, verletnumsum, ntype, n, dim, dimsq);
    }
    else {
        gpuKernelHalfNeighNum3D<<<gridDim, blockDim>>>(neighnum, x, ellipsoid, atomtype, 
                ilist, verletlist, verletnum, verletnumsum, ntype, n, dim, dimsq);
    }
}
template void gpuHalfNeighNum(int*, double*, double*, int*, int*, int*, int*, int*, int, int, int);
template void gpuHalfNeighNum(int*, float*, float*, int*, int*, int*, int*, int*, int, int, int);

template <typename T> __global__ void gpuKernelHalfNeighList2D(int *neighlist, T *x, T* ellipsoid, int *neighnum, 
   int *neighnumsum, int *atomtype, int *ilist, int *verletlist, int *verletnum, int *verletnumsum, int ntype, int inum, int dim, int dimsq)
{
    int ii = threadIdx.x + blockIdx.x * blockDim.x;
    while (ii < inum) {            
        int i = ilist[ii];  // atom i
        T xi0 = x[i*dim];        // position of atom i
        T xi1 = x[i*dim+1];      // position of atom i
        int ti = atomtype[i];      // element of atom i         
        int m = verletnum[i]; // number of atoms around i 
        int jstart = verletnumsum[i];     // starting 
        int nstart = neighnumsum[i];   
        int ninc = 0;                // increment 
        for (int l=0; l<m ; l++) {   // loop over each atom around atom i
            int j = verletlist[jstart+l];  // atom j     
            if (j >= i) {        
                int tj = atomtype[j];      // element of atom j
                int tij = dimsq*(ti + tj*ntype);
                T A00 = ellipsoid[tij]; // ellipsoid for pair (i,j)   
                T A10 = ellipsoid[tij+1]; // ellipsoid for pair (i,j)
                T A01 = ellipsoid[tij+2]; // ellipsoid for pair (i,j)  
                T A11 = ellipsoid[tij+3]; // ellipsoid for pair (i,j)                                                                                      
                T xij0 = x[j*dim] - xi0;  // xj - xi
                T xij1 = x[j*dim+1] - xi1; // xj - xi
                // distance between atom i and atom j 
                T rij = xij0*(A00*xij0 + A01*xij1) + xij1*(A10*xij0 + A11*xij1);                        
                if (rij <= 1.0) { // if atom j is inside the ellipsoid 
                    neighlist[nstart+ninc] = j;
                    ninc += 1;  // increase the number of neighbors by 1               
                }
            }
        }
        neighnum[i] = ninc; // number of neighbors of atom i

        ii += blockDim.x * gridDim.x;                
    }
}

template <typename T> __global__ void gpuKernelHalfNeighList3D(int *neighlist, T *x, T* ellipsoid, int *neighnum, 
   int *neighnumsum, int *atomtype, int *ilist, int *verletlist, int *verletnum, int *verletnumsum, int ntype, int inum, int dim, int dimsq)
{
    int ii = threadIdx.x + blockIdx.x * blockDim.x;
    while (ii < inum) {            
        int i = ilist[ii];  // atom i
        T xi0 = x[i*dim];        // position of atom i
        T xi1 = x[i*dim+1];      // position of atom i
        T xi2 = x[i*dim+2];      // position of atom i
        int ti = atomtype[i];      // element of atom i 
        int m = verletnum[i]; // number of atoms around i 
        int jstart = verletnumsum[i];     // starting 
        int nstart = neighnumsum[i];   
        int ninc = 0;                // increment 
        for (int l=0; l<m ; l++) {   // loop over each atom around atom i
            int j = verletlist[jstart+l];  // atom j          
            if (j >= i) {       
                int tj = atomtype[j];      // element of atom j
                int tij = dimsq*(ti + tj*ntype);
                T A00 = ellipsoid[tij]; // ellipsoid for element t   
                T A10 = ellipsoid[tij+1]; // ellipsoid for element t   
                T A20 = ellipsoid[tij+2]; // ellipsoid for element t   
                T A01 = ellipsoid[tij+3]; // ellipsoid for element t   
                T A11 = ellipsoid[tij+4]; // ellipsoid for element t   
                T A21 = ellipsoid[tij+5]; // ellipsoid for element t   
                T A02 = ellipsoid[tij+6]; // ellipsoid for element t   
                T A12 = ellipsoid[tij+7]; // ellipsoid for element t   
                T A22 = ellipsoid[tij+8]; // ellipsoid for element t                                                                                               
                T xij0 = x[j*dim] - xi0;  // xj - xi
                T xij1 = x[j*dim+1] - xi1; // xj - xi
                T xij2 = x[j*dim+2] - xi2; // xj - xi
                // distance between atom i and atom j 
                T rij = xij0*(A00*xij0 + A01*xij1 + A02*xij2) + xij1*(A10*xij0 + A11*xij1 + A12*xij2) + xij2*(A20*xij0 + A21*xij1 + A22*xij2);
                if (rij <= 1.0) { // if atom j is inside the ellipsoid 
                    neighlist[nstart+ninc] = j;
                    ninc += 1;  // increase the number of neighbors by 1               
                }
            }
        }            
        neighnum[i] = ninc; // number of neighbors of atom i

        ii += blockDim.x * gridDim.x;                
    }
}

template <typename T> void gpuHalfNeighList(int *neighlist, T *x, T* ellipsoid, int *neighnum, int *neighnumsum, 
      int *atomtype, int *ilist, int *verletlist, int *verletnum, int *verletnumsum, int *tm1, int *tm2, int ntype, int inum, int dim)
{                        
    int n = inum;
    int blockDim = 256;
    int gridDim = (n + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;

    // a list contains the starting positions of the first neighbor 
    gpuCumsum(neighnumsum, neighnum, tm1, tm2, inum+1); 

    int dimsq = dim*dim;
    if (dim==2) {
        gpuKernelHalfNeighList2D<<<gridDim, blockDim>>>(neighlist, x, ellipsoid, neighnum, neighnumsum, atomtype, 
                ilist, verletlist, verletnum, verletnumsum, ntype, n, dim, dimsq);
    }
    else {
        gpuKernelHalfNeighList3D<<<gridDim, blockDim>>>(neighlist, x, ellipsoid, neighnum, neighnumsum, atomtype, 
                ilist, verletlist, verletnum, verletnumsum, ntype, n, dim, dimsq);
    }
}
template void gpuHalfNeighList(int*, double*, double*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int, int, int);
template void gpuHalfNeighList(int*, float*, float*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int, int, int);


template <typename T> __global__ void gpuKernelGetNeighPairs2D(T *xij, T *x, int *ti, int *tj, int *ilist, 
        int *neighlist,  int *neighnum, int *neighnumsum, int *atomtype, int ntype, int inum, int dim)
{
    int ii = threadIdx.x + blockIdx.x * blockDim.x;
    while (ii < inum) {            
        int i = ilist[ii];  // atom i
        int itype = atomtype[i];
        T xi0 = x[i*dim];        // position of atom i
        T xi1 = x[i*dim+1];      // position of atom i
        int m = neighnum[i];     // number of neighbors around i 
        int nstart = neighnumsum[i];   
        for (int l=0; l<m ; l++) {   // loop over each atom around atom i
            int k = nstart+l;
            int j = neighlist[k];  // atom j     
            xij[k*dim]   = x[j*dim] - xi0;  // xj - xi
            xij[k*dim+1] = x[j*dim+1] - xi1; // xj - xi      
            ti[k]       = itype;
            tj[k]       = atomtype[j];
        }
        ii += blockDim.x * gridDim.x;      
    }
}

template <typename T> __global__ void gpuKernelGetNeighPairs3D(T *xij, T *x, int *ti, int *tj, int *ilist, 
        int *neighlist,  int *neighnum, int *neighnumsum, int *atomtype, int ntype, int inum, int dim)
{
    int ii = threadIdx.x + blockIdx.x * blockDim.x;
    while (ii < inum) {            
        int i = ilist[ii];  // atom i
        int itype = atomtype[i];
        T xi0 = x[i*dim];        // position of atom i
        T xi1 = x[i*dim+1];      // position of atom i
        T xi2 = x[i*dim+2];      // position of atom i
        int m = neighnum[i];     // number of neighbors around i 
        int nstart = neighnumsum[i];               
        for (int l=0; l<m ; l++) {   // loop over each atom around atom i
            int k = nstart+l;
            int j = neighlist[k];  // atom j     
            xij[k*dim]   = x[j*dim] - xi0;  // xj - xi
            xij[k*dim+1] = x[j*dim+1] - xi1; // xj - xi                                                
            xij[k*dim+2] = x[j*dim+2] - xi2; // xj - xi       
            ti[k]       = itype;
            tj[k]       = atomtype[j];
        }
        ii += blockDim.x * gridDim.x;      
    }
}

template <typename T> void gpuGetNeighPairs(T *xij, T *x, int *ti, int *tj, int *ilist, int *neighlist,  
        int *neighnum, int *neighnumsum, int *atomtype, int ntype, int inum, int dim)
{                        
    int n = inum;
    int blockDim = 256;
    int gridDim = (n + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;

    if (dim==2) {
        gpuKernelGetNeighPairs2D<<<gridDim, blockDim>>>(xij, x, ti, tj, ilist, neighlist, neighnum, 
                neighnumsum, atomtype, ntype, n, dim);
    }
    else {
        gpuKernelGetNeighPairs3D<<<gridDim, blockDim>>>(xij, x, ti, tj, ilist, neighlist, neighnum, 
                neighnumsum, atomtype, ntype, n, dim);
    }
}
template void gpuGetNeighPairs(double*, double*, int*, int*, int*, int*, int*, int*, int*, int, int, int);
template void gpuGetNeighPairs(float*, float*, int*, int*, int*, int*, int*, int*, int*, int, int, int);


template <typename T> __global__ void gpuKernelGetNeighPairs2D(T *xi, T *xj, T *x, int *ti, int *tj, int *ilist, 
        int *neighlist,  int *neighnum, int *neighnumsum, int *atomtype, int ntype, int inum, int dim)
{
    int ii = threadIdx.x + blockIdx.x * blockDim.x;
    while (ii < inum) {            
        int i = ilist[ii];  // atom i
        int itype = atomtype[i];
        T xi0 = x[i*dim];        // position of atom i
        T xi1 = x[i*dim+1];      // position of atom i
        int m = neighnum[i];     // number of neighbors around i 
        int nstart = neighnumsum[i];   
        for (int l=0; l<m ; l++) {   // loop over each atom around atom i
            int k = nstart+l;
            int j = neighlist[k];  // atom j     
            xi[k*dim]   = xi0;        // position of atom i
            xi[k*dim+1] = xi1;       // position of atom i                        
            xj[k*dim]   = x[j*dim];  // xj - xi
            xj[k*dim+1] = x[j*dim+1]; // xj - xi      
            ti[k]       = itype;
            tj[k]       = atomtype[j];
        }
        ii += blockDim.x * gridDim.x;      
    }
}

template <typename T> __global__ void gpuKernelGetNeighPairs3D(T *xi, T *xj, T *x, int *ti, int *tj, int *ilist, 
        int *neighlist,  int *neighnum, int *neighnumsum, int *atomtype, int ntype, int inum, int dim)
{
    int ii = threadIdx.x + blockIdx.x * blockDim.x;
    while (ii < inum) {            
        int i = ilist[ii];  // atom i
        int itype = atomtype[i];
        T xi0 = x[i*dim];        // position of atom i
        T xi1 = x[i*dim+1];      // position of atom i
        T xi2 = x[i*dim+2];      // position of atom i
        int m = neighnum[i];     // number of neighbors around i 
        int nstart = neighnumsum[i];               
        for (int l=0; l<m ; l++) {   // loop over each atom around atom i
            int k = nstart+l;
            int j = neighlist[k];  // atom j     
            xi[k*dim]   = xi0;        // position of atom i
            xi[k*dim+1] = xi1;       // position of atom i                
            xi[k*dim+2] = xi2;       // position of atom i                        
            xj[k*dim]   = x[j*dim];  // xj - xi
            xj[k*dim+1] = x[j*dim+1]; // xj - xi                                                
            xj[k*dim+2] = x[j*dim+2]; // xj - xi       
            ti[k]       = itype;
            tj[k]       = atomtype[j];
        }
        ii += blockDim.x * gridDim.x;      
    }
}

template <typename T> void gpuGetNeighPairs(T *xi, T *xj, T *x, int *ti, int *tj, int *ilist, int *neighlist,  
        int *neighnum, int *neighnumsum, int *atomtype, int ntype, int inum, int dim)
{                        
    int n = inum;
    int blockDim = 256;
    int gridDim = (n + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;

    if (dim==2) {
        gpuKernelGetNeighPairs2D<<<gridDim, blockDim>>>(xi, xj, x, ti, tj, ilist, neighlist, neighnum, 
                neighnumsum, atomtype, ntype, n, dim);
    }
    else {
        gpuKernelGetNeighPairs3D<<<gridDim, blockDim>>>(xi, xj, x, ti, tj, ilist, neighlist, neighnum, 
                neighnumsum, atomtype, ntype, n, dim);
    }
}
template void gpuGetNeighPairs(double*, double*, double*, int*, int*, int*, int*, int*, int*, int*, int, int, int);
template void gpuGetNeighPairs(float*, float*, float*, int*, int*, int*, int*, int*, int*, int*, int, int, int);


// template <typename T>
// __global__ void gpuKernelVerletAtoms2D(int *verletnum, T *x, T *ellipsoid, int *atomtype, 
//         int *clist, int *c2alist, int *c2anum, int *c2anumsum, int *nc, int inum, int dim)
// {
//     int i = threadIdx.x + blockIdx.x * blockDim.x;
//     while (i < inum) {    
//         verletnum[i] = 0;
//         T xi0 = x[i*dim];        // position of atom i
//         T xi1 = x[i*dim+1];      // position of atom i
//         int t = atomtype[i];      // element of atom i 
//         T A00 = ellipsoid[dim*dim*t]; // ellipsoid for element t   
//         T A10 = ellipsoid[dim*dim*t+1]; // ellipsoid for element t   
//         T A01 = ellipsoid[dim*dim*t+2]; // ellipsoid for element t   
//         T A11 = ellipsoid[dim*dim*t+3]; // ellipsoid for element t   
//         int j = clist[i];         // cell j contains atom i           
//         int j1 = j%nc[0];
//         int j2 = (j-j1)/nc[0];
//         for (int i1=-1; i1<=1; i1++) {
//             int k1 = j1 + i1;
//             for (int i2=-1; i2<=1; i2++) {
//                 int k2 = j2 + i2;
//                 int k = k1 + nc[0]*k2; // cell k
//                 int m = c2anum[k];     // number of atoms in cell k
//                 int s = c2anumsum[k];  // starting position of the first atom in cell k
//                 for (int l=0; l<m ; l++) {
//                     j = c2alist[s+l];  // atom j
//                     T xij0 = x[j*dim] - xi0;  // xj - xi
//                     T xij1 = x[j*dim+1] - xi1; // xj - xi
//                     // distance between atom i and atom j 
//                     T rij = xij0*(A00*xij0 + A01*xij1) + xij1*(A10*xij0 + A11*xij1);                        
//                     if (rij <= 1.0)
//                         verletnum[i] += 1;
//                 }
//             }
//         }                
//         i += blockDim.x * gridDim.x;
//     }    
// }
// 
// template <typename T>
// __global__ void gpuKernelVerletAtoms3D(int *verletnum, T *x, T *ellipsoid, int *atomtype, 
//         int *clist, int *c2alist, int *c2anum, int *c2anumsum, int *nc, int inum, int dim)
// {
//     int i = threadIdx.x + blockIdx.x * blockDim.x;
//     while (i < inum) {    
//         verletnum[i] = 0;
//         T xi0 = x[i*dim];        // position of atom i
//         T xi1 = x[i*dim+1];      // position of atom i
//         T xi2 = x[i*dim+2];      // position of atom i
//         int t = atomtype[i];      // element of atom i 
//         T A00 = ellipsoid[dim*dim*t]; // ellipsoid for element t   
//         T A10 = ellipsoid[dim*dim*t+1]; // ellipsoid for element t   
//         T A20 = ellipsoid[dim*dim*t+2]; // ellipsoid for element t   
//         T A01 = ellipsoid[dim*dim*t+3]; // ellipsoid for element t   
//         T A11 = ellipsoid[dim*dim*t+4]; // ellipsoid for element t   
//         T A21 = ellipsoid[dim*dim*t+5]; // ellipsoid for element t   
//         T A02 = ellipsoid[dim*dim*t+6]; // ellipsoid for element t   
//         T A12 = ellipsoid[dim*dim*t+7]; // ellipsoid for element t   
//         T A22 = ellipsoid[dim*dim*t+8]; // ellipsoid for element t   
//         int j = clist[i];         // cell j contains atom i       
//         int n = j%(nc[0]*nc[1]);
//         int j3 = (j-n)/(nc[0]*nc[1]);            
//         int j1 = n%nc[0];
//         int j2 = (j-j1)/nc[0];
//         for (int i1=-1; i1<=1; i1++) {
//             int k1 = j1 + i1;
//             for (int i2=-1; i2<=1; i2++) {
//                 int k2 = j2 + i2;
//                 for (int i3=-1; i3<=1; i3++) {
//                     int k3 = j3 + i3;                    
//                     int k = k1 + nc[0]*k2 + nc[0]*nc[1]*k3; // cell k
//                     int m = c2anum[k];     // number of atoms in cell k
//                     int s = c2anumsum[k];  // starting position of the first atom in cell k
//                     for (int l=0; l<m ; l++) {
//                         j = c2alist[s+l];  // atom j
//                         T xij0 = x[j*dim] - xi0;  // xj - xi
//                         T xij1 = x[j*dim+1] - xi1; // xj - xi
//                         T xij2 = x[j*dim+2] - xi2; // xj - xi
//                         // distance between atom i and atom j 
//                         T rij = xij0*(A00*xij0 + A01*xij1 + A02*xij2) + xij1*(A10*xij0 + A11*xij1 + A12*xij2) + xij2*(A20*xij0 + A21*xij1 + A22*xij2);
//                         if (rij <= 1.0)
//                             verletnum[i] += 1;
//                     }
//                 }
//             }
//         }                        
//         i += blockDim.x * gridDim.x;
//     }    
// }
// 
// template <typename T> void gpuVerletAtoms(int *verletnum, T *x, T *ellipsoid, int *atomtype, 
//         int *clist, int *c2alist, int *c2anum, int *c2anumsum, int *nc, int inum, int dim)
// {                        
//     int n = inum;
// 
//     int blockDim = 256;
//     int gridDim = (n + blockDim - 1) / blockDim;
//     gridDim = (gridDim>1024)? 1024 : gridDim;
//     if (dim==2) {
//         gpuKernelVerletAtoms2D<<<gridDim, blockDim>>>(verletnum, x, ellipsoid, atomtype, clist, c2alist, 
//                 c2anum, c2anumsum, nc, n, dim);
//     }
//     else {
//         gpuKernelVerletAtoms3D<<<gridDim, blockDim>>>(verletnum, x, ellipsoid, atomtype, clist, c2alist, 
//                 c2anum, c2anumsum, nc, n, dim);
//     }
// }
// template void gpuVerletAtoms(int*, double*, double*, int*, int*, int*, int*, int*, int*, int, int);
// template void gpuVerletAtoms(int*, float*, float*, int*, int*, int*, int*, int*, int*, int, int);
// 
// template <typename T>
// __global__ void gpuKernelCreateVerletList2D(int *verletlist, T *x, T *ellipsoid, int *verletnum, int *verletnumsum, int *atomtype, 
//       int *ilist, int *clist, int *c2alist, int *c2anum, int *c2anumsum, int *nc, int inum, int dim)
// {
//     int i = threadIdx.x + blockIdx.x * blockDim.x;
//     while (i < inum) {    
//         T xi0 = x[i*dim];        // position of atom i
//         T xi1 = x[i*dim+1];      // position of atom i
//         int t = atomtype[i];      // element of atom i 
//         T A00 = ellipsoid[dim*dim*t]; // ellipsoid for element t   
//         T A10 = ellipsoid[dim*dim*t+1]; // ellipsoid for element t   
//         T A01 = ellipsoid[dim*dim*t+2]; // ellipsoid for element t   
//         T A11 = ellipsoid[dim*dim*t+3]; // ellipsoid for element t   
//         int nstart = verletnumsum[i];     // starting 
//         int ninc = 0;                // increment
//         int j = clist[i];         // cell j contains atom i           
//         int j1 = j%nc[0];
//         int j2 = (j-j1)/nc[0];
//         for (int i1=-1; i1<=1; i1++) {
//             int k1 = j1 + i1;
//             for (int i2=-1; i2<=1; i2++) {
//                 int k2 = j2 + i2;
//                 int k = k1 + nc[0]*k2; // cell k
//                 int m = c2anum[k];     // number of atoms in cell k
//                 int s = c2anumsum[k];  // starting position of the first atom in cell k
//                 for (int l=0; l<m ; l++) { //loop over each atom j in cell k
//                     j = c2alist[s+l];  // atom j
//                     T xij0 = x[j*dim] - xi0;  // xj - xi
//                     T xij1 = x[j*dim+1] - xi1; // xj - xi
//                     // distance between atom i and atom j 
//                     T rij = xij0*(A00*xij0 + A01*xij1) + xij1*(A10*xij0 + A11*xij1);                        
//                     if (rij <= 1.0) {
//                         verletlist[nstart + ninc] = j; // add atom j into the list
//                         ninc += 1;
//                     }
//                 }
//             }
//         }               
//         verletnum[i] = ninc;
//         i += blockDim.x * gridDim.x;
//     }
// }
// 
// template <typename T>
// __global__ void gpuKernelCreateVerletList3D(int *verletlist, T *x, T *ellipsoid, int *verletnum, int *verletnumsum, int *atomtype, 
//       int *ilist, int *clist, int *c2alist, int *c2anum, int *c2anumsum, int *nc, int inum, int dim)
// {
//     int i = threadIdx.x + blockIdx.x * blockDim.x;
//     while (i < inum) {    
//         T xi0 = x[i*dim];        // position of atom i
//         T xi1 = x[i*dim+1];      // position of atom i
//         T xi2 = x[i*dim+2];      // position of atom i
//         int t = atomtype[i];      // element of atom i 
//         T A00 = ellipsoid[dim*dim*t]; // ellipsoid for element t   
//         T A10 = ellipsoid[dim*dim*t+1]; // ellipsoid for element t   
//         T A20 = ellipsoid[dim*dim*t+2]; // ellipsoid for element t   
//         T A01 = ellipsoid[dim*dim*t+3]; // ellipsoid for element t   
//         T A11 = ellipsoid[dim*dim*t+4]; // ellipsoid for element t   
//         T A21 = ellipsoid[dim*dim*t+5]; // ellipsoid for element t   
//         T A02 = ellipsoid[dim*dim*t+6]; // ellipsoid for element t   
//         T A12 = ellipsoid[dim*dim*t+7]; // ellipsoid for element t   
//         T A22 = ellipsoid[dim*dim*t+8]; // ellipsoid for element t   
//         int nstart = verletnumsum[i];     // starting 
//         int ninc = 0;                // increment            
//         int j = clist[i];         // cell j contains atom i       
//         int n = j%(nc[0]*nc[1]);
//         int j3 = (j-n)/(nc[0]*nc[1]);            
//         int j1 = n%nc[0];
//         int j2 = (j-j1)/nc[0];
//         for (int i1=-1; i1<=1; i1++) {
//             int k1 = j1 + i1;
//             for (int i2=-1; i2<=1; i2++) {
//                 int k2 = j2 + i2;
//                 for (int i3=-1; i3<=1; i3++) {
//                     int k3 = j3 + i3;                    
//                     int k = k1 + nc[0]*k2 + nc[0]*nc[1]*k3; // cell k
//                     int m = c2anum[k];     // number of atoms in cell k
//                     int s = c2anumsum[k];  // starting position of the first atom in cell k
//                     for (int l=0; l<m ; l++) {
//                         j = c2alist[s+l];  // atom j
//                         T xij0 = x[j*dim] - xi0;  // xj - xi
//                         T xij1 = x[j*dim+1] - xi1; // xj - xi
//                         T xij2 = x[j*dim+2] - xi2; // xj - xi
//                         // distance between atom i and atom j 
//                         T rij = xij0*(A00*xij0 + A01*xij1 + A02*xij2) + xij1*(A10*xij0 + A11*xij1 + A12*xij2) + xij2*(A20*xij0 + A21*xij1 + A22*xij2);
//                         if (rij <= 1.0) {
//                             verletlist[nstart + ninc] = j; // add atom j into the list
//                             ninc += 1;
//                         }                            
//                     }
//                 }
//             }
//         }                
//         verletnum[i] = ninc;
//         i += blockDim.x * gridDim.x;
//     }
// }
// 
// 
// template <typename T> void gpuCreateVerletList(int *verletlist, T *x, T *ellipsoid, int *verletnum,  int *verletnumsum, int *tm1, int *tm2,
//      int *atomtype, int *ilist, int *clist, int *c2alist, int *c2anum, int *c2anumsum, int *nc, int inum, int dim)
// {                        
//     int n = inum;
//     int blockDim = 256;
//     int gridDim = (n + blockDim - 1) / blockDim;
//     gridDim = (gridDim>1024)? 1024 : gridDim;
// 
//     // a list contains the starting positions of the first j atom 
//     gpuCumsum(verletnumsum, verletnum, tm1, tm2, inum+1); 
// 
//     if (dim==2) {
//         gpuKernelCreateVerletList2D<<<gridDim, blockDim>>>(verletlist, x, ellipsoid, verletnum, verletnumsum, atomtype, 
//                 ilist, clist, c2alist, c2anum, c2anumsum, nc, n, dim);
//     }
//     else {
//         gpuKernelCreateVerletList3D<<<gridDim, blockDim>>>(verletlist, x, ellipsoid, verletnum, verletnumsum, atomtype, 
//                 ilist, clist, c2alist, c2anum, c2anumsum, nc, n, dim);
//     }
// }
// template void gpuCreateVerletList(int*, double*, double*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int, int);
// template void gpuCreateVerletList(int*, float*, float*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int, int);
//       
// template <typename T> __global__ void gpuKernelFullNeighNum2D(int *neighnum, T *x, T* ellipsoid, int *atomtype, 
//         int *ilist, int *verletlist, int *verletnum, int *verletnumsum, int inum, int dim, int dimsq)
// {
//     int ii = threadIdx.x + blockIdx.x * blockDim.x;
//     while (ii < inum) {            
//         int i = ilist[ii];  // atom i
//         T xi0 = x[i*dim];        // position of atom i
//         T xi1 = x[i*dim+1];      // position of atom i
//         int t = atomtype[i];      // element of atom i 
//         T A00 = ellipsoid[dimsq*t]; // ellipsoid for element t   
//         T A10 = ellipsoid[dimsq*t+1]; // ellipsoid for element t   
//         T A01 = ellipsoid[dimsq*t+2]; // ellipsoid for element t   
//         T A11 = ellipsoid[dimsq*t+3]; // ellipsoid for element t      
//         int m = verletnum[i]; // number of atoms around i 
//         int start = verletnumsum[i];     // starting 
//         int ninc = 0;                // increment 
//         for (int l=0; l<m ; l++) {   // loop over each atom around atom i
//             int j = verletlist[start+l];  // atom j                
//             T xij0 = x[j*dim] - xi0;  // xj - xi
//             T xij1 = x[j*dim+1] - xi1; // xj - xi
//             // distance between atom i and atom j 
//             T rij = xij0*(A00*xij0 + A01*xij1) + xij1*(A10*xij0 + A11*xij1);                        
//             if (rij <= 1.0) // if atom j is inside the ellipsoid 
//                 ninc += 1;  // increase the number of neighbors by 1               
//         }
//         neighnum[i] = ninc; // number of neighbors of atom i
// 
//         ii += blockDim.x * gridDim.x;                
//     }
// }
// 
// template <typename T> __global__ void gpuKernelFullNeighNum3D(int *neighnum, T *x, T* ellipsoid, int *atomtype, 
//         int *ilist, int *verletlist, int *verletnum, int *verletnumsum, int inum, int dim, int dimsq)
// {
//     int ii = threadIdx.x + blockIdx.x * blockDim.x;
//     while (ii < inum) {            
//         int i = ilist[ii];  // atom i
//         T xi0 = x[i*dim];        // position of atom i
//         T xi1 = x[i*dim+1];      // position of atom i
//         T xi2 = x[i*dim+2];      // position of atom i
//         int t = atomtype[i];      // element of atom i 
//         T A00 = ellipsoid[dimsq*t]; // ellipsoid for element t   
//         T A10 = ellipsoid[dimsq*t+1]; // ellipsoid for element t   
//         T A20 = ellipsoid[dimsq*t+2]; // ellipsoid for element t   
//         T A01 = ellipsoid[dimsq*t+3]; // ellipsoid for element t   
//         T A11 = ellipsoid[dimsq*t+4]; // ellipsoid for element t   
//         T A21 = ellipsoid[dimsq*t+5]; // ellipsoid for element t   
//         T A02 = ellipsoid[dimsq*t+6]; // ellipsoid for element t   
//         T A12 = ellipsoid[dimsq*t+7]; // ellipsoid for element t   
//         T A22 = ellipsoid[dimsq*t+8]; // ellipsoid for element t   
//         int m = verletnum[i]; // number of atoms around i 
//         int start = verletnumsum[i];     // starting 
//         int ninc = 0;                // increment 
//         for (int l=0; l<m ; l++) {   // loop over each atom around atom i
//             int j = verletlist[start+l];  // atom j                
//             T xij0 = x[j*dim] - xi0;  // xj - xi
//             T xij1 = x[j*dim+1] - xi1; // xj - xi
//             T xij2 = x[j*dim+2] - xi2; // xj - xi
//             // distance between atom i and atom j 
//             T rij = xij0*(A00*xij0 + A01*xij1 + A02*xij2) + xij1*(A10*xij0 + A11*xij1 + A12*xij2) + xij2*(A20*xij0 + A21*xij1 + A22*xij2);
//             if (rij <= 1.0) // if atom j is inside the ellipsoid 
//                 ninc += 1;  // increase the number of neighbors by 1                               
//         }            
//         neighnum[i] = ninc; // number of neighbors of atom i        
// 
//         ii += blockDim.x * gridDim.x;                
//     }
// }
// 
// template <typename T> void gpuFullNeighNum(int *neighnum, T *x, T* ellipsoid, int *atomtype, 
//         int *ilist, int *verletlist, int *verletnum, int *verletnumsum, int inum, int dim)
// {                        
//     int n = inum;
//     int blockDim = 256;
//     int gridDim = (n + blockDim - 1) / blockDim;
//     gridDim = (gridDim>1024)? 1024 : gridDim;
// 
//     int dimsq = dim*dim;
//     if (dim==2) {
//         gpuKernelFullNeighNum2D<<<gridDim, blockDim>>>(neighnum, x, ellipsoid, atomtype, 
//                 ilist, verletlist, verletnum, verletnumsum, n, dim, dimsq);
//     }
//     else {
//         gpuKernelFullNeighNum3D<<<gridDim, blockDim>>>(neighnum, x, ellipsoid, atomtype, 
//                 ilist, verletlist, verletnum, verletnumsum, n, dim, dimsq);
//     }
// }
// template void gpuFullNeighNum(int*, double*, double*, int*, int*, int*, int*, int*, int, int);
// template void gpuFullNeighNum(int*, float*, float*, int*, int*, int*, int*, int*, int, int);
// 
// template <typename T> __global__ void gpuKernelFullNeighList2D(int *neighlist, T *x, T* ellipsoid, int *neighnum, 
//    int *neighnumsum, int *atomtype, int *ilist, int *verletlist, int *verletnum, int *verletnumsum, int inum, int dim, int dimsq)
// {
//     int ii = threadIdx.x + blockIdx.x * blockDim.x;
//     while (ii < inum) {            
//         int i = ilist[ii];  // atom i
//         T xi0 = x[i*dim];        // position of atom i
//         T xi1 = x[i*dim+1];      // position of atom i
//         int t = atomtype[i];      // element of atom i 
//         T A00 = ellipsoid[dimsq*t]; // ellipsoid for element t   
//         T A10 = ellipsoid[dimsq*t+1]; // ellipsoid for element t   
//         T A01 = ellipsoid[dimsq*t+2]; // ellipsoid for element t   
//         T A11 = ellipsoid[dimsq*t+3]; // ellipsoid for element t      
//         int m = verletnum[i]; // number of atoms around i 
//         int jstart = verletnumsum[i];     // starting 
//         int nstart = neighnumsum[i];   
//         int ninc = 0;                // increment 
//         for (int l=0; l<m ; l++) {   // loop over each atom around atom i
//             int j = verletlist[jstart+l];  // atom j                
//             T xij0 = x[j*dim] - xi0;  // xj - xi
//             T xij1 = x[j*dim+1] - xi1; // xj - xi
//             // distance between atom i and atom j 
//             T rij = xij0*(A00*xij0 + A01*xij1) + xij1*(A10*xij0 + A11*xij1);                        
//             if (rij <= 1.0) { // if atom j is inside the ellipsoid 
//                 neighlist[nstart+ninc] = j;
//                 ninc += 1;  // increase the number of neighbors by 1               
//             }
//         }
//         neighnum[i] = ninc; // number of neighbors of atom i
// 
//         ii += blockDim.x * gridDim.x;                
//     }
// }
// 
// template <typename T> __global__ void gpuKernelFullNeighList3D(int *neighlist, T *x, T* ellipsoid, int *neighnum, 
//    int *neighnumsum, int *atomtype, int *ilist, int *verletlist, int *verletnum, int *verletnumsum, int inum, int dim, int dimsq)
// {
//     int ii = threadIdx.x + blockIdx.x * blockDim.x;
//     while (ii < inum) {            
//         int i = ilist[ii];  // atom i
//         T xi0 = x[i*dim];        // position of atom i
//         T xi1 = x[i*dim+1];      // position of atom i
//         T xi2 = x[i*dim+2];      // position of atom i
//         int t = atomtype[i];      // element of atom i 
//         T A00 = ellipsoid[dimsq*t]; // ellipsoid for element t   
//         T A10 = ellipsoid[dimsq*t+1]; // ellipsoid for element t   
//         T A20 = ellipsoid[dimsq*t+2]; // ellipsoid for element t   
//         T A01 = ellipsoid[dimsq*t+3]; // ellipsoid for element t   
//         T A11 = ellipsoid[dimsq*t+4]; // ellipsoid for element t   
//         T A21 = ellipsoid[dimsq*t+5]; // ellipsoid for element t   
//         T A02 = ellipsoid[dimsq*t+6]; // ellipsoid for element t   
//         T A12 = ellipsoid[dimsq*t+7]; // ellipsoid for element t   
//         T A22 = ellipsoid[dimsq*t+8]; // ellipsoid for element t   
//         int m = verletnum[i]; // number of atoms around i 
//         int jstart = verletnumsum[i];     // starting 
//         int nstart = neighnumsum[i];   
//         int ninc = 0;                // increment 
//         for (int l=0; l<m ; l++) {   // loop over each atom around atom i
//             int j = verletlist[jstart+l];  // atom j                
//             T xij0 = x[j*dim] - xi0;  // xj - xi
//             T xij1 = x[j*dim+1] - xi1; // xj - xi
//             T xij2 = x[j*dim+2] - xi2; // xj - xi
//             // distance between atom i and atom j 
//             T rij = xij0*(A00*xij0 + A01*xij1 + A02*xij2) + xij1*(A10*xij0 + A11*xij1 + A12*xij2) + xij2*(A20*xij0 + A21*xij1 + A22*xij2);
//             if (rij <= 1.0) { // if atom j is inside the ellipsoid 
//                 neighlist[nstart+ninc] = j;
//                 ninc += 1;  // increase the number of neighbors by 1               
//             }
//         }            
//         neighnum[i] = ninc; // number of neighbors of atom i
// 
//         ii += blockDim.x * gridDim.x;                
//     }
// }
// 
// template <typename T> void gpuFullNeighList(int *neighlist, T *x, T* ellipsoid, int *neighnum, int *neighnumsum, 
//       int *atomtype, int *ilist, int *verletlist, int *verletnum, int *verletnumsum, int *tm1, int *tm2, int inum, int dim)
// {                        
//     int n = inum;
//     int blockDim = 256;
//     int gridDim = (n + blockDim - 1) / blockDim;
//     gridDim = (gridDim>1024)? 1024 : gridDim;
// 
//     // a list contains the starting positions of the first neighbor 
//     gpuCumsum(neighnumsum, neighnum, tm1, tm2, inum+1); 
// 
//     int dimsq = dim*dim;
//     if (dim==2) {
//         gpuKernelFullNeighList2D<<<gridDim, blockDim>>>(neighlist, x, ellipsoid, neighnum, neighnumsum, atomtype, 
//                 ilist, verletlist, verletnum, verletnumsum, n, dim, dimsq);
//     }
//     else {
//         gpuKernelFullNeighList3D<<<gridDim, blockDim>>>(neighlist, x, ellipsoid, neighnum, neighnumsum, atomtype, 
//                 ilist, verletlist, verletnum, verletnumsum, n, dim, dimsq);
//     }
// }
// template void gpuFullNeighList(int*, double*, double*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int, int);
// template void gpuFullNeighList(int*, float*, float*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int, int);


#endif


