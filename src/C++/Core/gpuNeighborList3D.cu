#ifndef __GPUNEIGHBORLIST3D
#define __GPUNEIGHBORLIST3D

template <typename T>
__global__ void gpuKernelGhostAtoms3D(int *glistnum, int *inside, T *x, T *pimages, T *wc, T *s2rmap, int n, int m, int dim)
{
    int i = threadIdx.x + blockIdx.x * blockDim.x;
    while (i < n) {
        glistnum[i] = 0; // set the number of ghost atoms for atom i to 0
        for (int j=1; j<m; j++) { // loop over each periodic image of atom i
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
        i += blockDim.x * gridDim.x;
    }
}

template <typename T>
__global__ void gpuKernelAtomList3D(int *alist, int *inside, int *glistnumsum, 
        T *x, T *pimages, int n, int m, int dim)
{
    int i = threadIdx.x + blockIdx.x * blockDim.x;
    while (i < n) {
        alist[i] = i;         // add atom i to the list
        int q = n + glistnumsum[i]; // offset the starting position by n
        int k = 0;            
        for (int j=1; j<m; j++) { // loop over each periodic image of atom i
            if (inside[j-1 + i*(pnum-1)]) {      
                x[dim*(q+k)+0] = x[i*dim+0] + pimages[j*dim+0]; // add the periodic image as a ghost atom
                x[dim*(q+k)+1] = x[i*dim+1] + pimages[j*dim+1]; //
                x[dim*(q+k)+2] = x[i*dim+2] + pimages[j*dim+2]; //                
                alist[q+k] = i;       // add atom i to the list
                k += 1;            // increase the number of ghost atoms for atom i by 1    
            }
        }        
        i += blockDim.x * gridDim.x;
    }
}

template <typename T> void gpuAtomList3D(int *alist,  int *inside, int *glistnumsum, int *glistnum, 
        int *d_sums, int *d_incr, T *x, T *pimages, T *wc, T *s2rmap, int inum, int pnum, int dim)
{                        
    int blockDim = 256;
    int gridDim = (inum + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;
    
    // a list of ghost atoms
    gpuGhostAtoms3D<<<gridDim, blockDim>>>(glistnum, inside, x, pimages, wc, s2rmap, inum, pnum, dim);
    
    // a list contains the starting position of the ghost atom of every atom i
    gpuCumsum(glistnumsum, glistnum, d_sums, d_incr, n+1); 

    gpuKernelAtomList3D<<<gridDim, blockDim>>>(alist, inside, glistnumsum, atomtype, x, pimages, 
            wc, s2rmap, inum, pnum, dim);
}
template void gpuAtomList3D(int*, int*, int*, int*, int*, int*, double*, double*, double*, double*, int, int, int);
template void gpuAtomList3D(int*, int*, int*, int*, int*, int*, float*, float*, float*, float*, int, int, int);

template <typename T>
__global__ void gpuKernelCellList3D(int *clist, T *xi, T *eta1, T *eta2, T *eta3, T *s2rmap, 
        int *nc, int inum, int natom, int dim)
{
    int i = threadIdx.x + blockIdx.x * blockDim.x;
    while (i < natom) {
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
        i += blockDim.x * gridDim.x;
    }
}

template <typename T> void gpuCellList3D(int *clist, T *xi, T *eta1, T *eta2, T *eta3, T *s2rmap, 
        int *nc, int inum, int natom, int dim)
{                        
    int blockDim = 256;
    int gridDim = (natom + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;
      
    gpuKernelCellList3D<<<gridDim, blockDim>>>(clist, xi, eta1, eta2, eta3, s2rmap, 
        nc, inum, natom, dim);
}
template void gpuCellList3D(int*, double*, double*, double*, double*, double*, int*, int, int, int);
template void gpuCellList3D(int*, float*, float*, float*, float*, float*, int*, int, int, int);


template <typename T> 
__global__ void gpuKernelFullNeighborList3D(int *neighborlist, int *neighbornum, T *x, T *rcutsq, int *alist, 
        int *clist, int *c2alist, int *c2anumsum, int *nc, int inum, int jnum, int dim)
{            
    int ii = threadIdx.x + blockIdx.x * blockDim.x;
    while (ii < inum) {          // for each atom i in the simulation box    
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
        ii += blockDim.x * gridDim.x;
    }        
}

template <typename T> void gpuFullNeighborList3D(int *neighborlist, int *neighbornum, T *x, T *rcutsq,
        int *alist, int *clist, int *c2alist, int *c2anumsum, int *nc, int inum, int jnum, int dim)
{                        
    int blockDim = 256;
    int gridDim = (inum + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;
    
    gpuKernelFullNeighborList3D<<<gridDim, blockDim>>>(neighborlist, neighbornum, 
            x, rcutsq, alist, clist, c2alist, c2anumsum, nc, inum, jnum, dim);
}
template void gpuFullNeighborList3D(int*, int*, double*, double*, int*, int*, int*, int*, int*, int, int, int);
template void gpuFullNeighborList3D(int*, int*, float*, float*, int*, int*, int*, int*, int*, int, int, int);


template <typename T> 
__global__ void gpuKernelFullNeighborList3D(int *neighborlist, int *neighbornum, T *x, T *rcutsq, 
        int anum, int inum, int jnum, int dim)
{                
    int i = threadIdx.x + blockIdx.x * blockDim.x;
    while (i < inum) {          // for each atom ii in the simulation box    
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
        i += blockDim.x * gridDim.x;
    }        
}

template <typename T> void gpuFullNeighborList3D(int *neighborlist, int *neighbornum, T *x, T *rcutsq,
        int anum, int inum, int jnum, int dim)
{                        
    int blockDim = 256;
    int gridDim = (inum + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;
    
    gpuKernelFullNeighborList3D<<<gridDim, blockDim>>>(neighborlist, neighbornum, 
            x, rcutsq, anum, inum, jnum, dim);
}
template void gpuFullNeighborList3D(int*, int*, double*, double*, int, int, int, int);
template void gpuFullNeighborList3D(int*, int*, float*, float*, int, int, int, int);



#endif


