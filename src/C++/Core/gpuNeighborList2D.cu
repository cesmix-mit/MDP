#ifndef __GPUNEIGHBORLIST2D
#define __GPUNEIGHBORLIST2D

#include "cub112/device/device_radix_sort.cuh"

template <typename T>
__global__ void gpuKernelGhostAtoms2D(int *glistnum, int *inside, T *x, T *pimages, T *wc, 
    T *s2rmap, int inum, int pnum, int dim)
{
    int i = threadIdx.x + blockIdx.x * blockDim.x;
    while (i < inum) {
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
        i += blockDim.x * gridDim.x;
    }
}


template <typename T>
__global__ void gpuKernelAtomList2D(int *alist, int *inside, int *glistnumsum,
        T *x, T *pimages, int inum, int pnum, int dim)
{
    int i = threadIdx.x + blockIdx.x * blockDim.x;
    while (i < inum) {
        alist[i] = i;         // add atom i to the list
        int q = inum + glistnumsum[i]; // offset the starting position by n
        int k = 0;            
        for (int j=1; j<pnum; j++) { // loop over each periodic image of atom i            
            /// check if the mapped point is inside the bounding box
            if (inside[j-1 + i*(pnum-1)]) {  
                x[dim*(q+k)+0] = x[i*dim+0] + pimages[j*dim+0]; // add the periodic image as a ghost atom
                x[dim*(q+k)+1] = x[i*dim+1] + pimages[j*dim+1]; //
                alist[q+k] = i;       // add atom i to the list
                k += 1;            // increase the number of ghost atoms for atom i by 1    
            }
        }        
        i += blockDim.x * gridDim.x;
    }
}

template <typename T> void gpuAtomList2D(int *alist,  int *inside, int *glistnumsum, int *glistnum, 
        int *d_sums, int *d_incr, T *x, T *pimages, T *wc, T *s2rmap, int inum, int pnum, int dim)
{                        
    int blockDim = 256;
    int gridDim = (inum + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;
    
    // a list of ghost atoms
    gpuKernelGhostAtoms2D<<<gridDim, blockDim>>>(glistnum, inside, x, pimages, wc, s2rmap, inum, pnum, dim);
    
    // a list contains the starting position of the ghost atom of every atom i
    gpuCumsum(glistnumsum, glistnum, d_sums, d_incr, inum+1); 

    gpuKernelAtomList2D<<<gridDim, blockDim>>>(alist, inside, glistnumsum, x, pimages, inum, pnum, dim);
}
template void gpuAtomList2D(int*, int*, int*, int*, int*, int*, double*, double*, double*, double*, int, int, int);
template void gpuAtomList2D(int*, int*, int*, int*, int*, int*, float*, float*, float*, float*, int, int, int);

template <typename T>
__global__ void gpuKernelCellList2D(int *clist, T *x, T *eta1, T *eta2, T *eta3, T *s2rmap, 
        int *nc, int inum, int natom, int dim)
{
    int i = threadIdx.x + blockIdx.x * blockDim.x;
    while (i < natom) {
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
        i += blockDim.x * gridDim.x;
    }
}

template <typename T> void gpuCellList2D(int *clist, T *xi, T *eta1, T *eta2, T *eta3, T *s2rmap, 
        int *nc, int inum, int natom, int dim)
{                        
    int blockDim = 256;
    int gridDim = (natom + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;
      
    gpuKernelCellList2D<<<gridDim, blockDim>>>(clist, xi, eta1, eta2, eta3, s2rmap, 
        nc, inum, natom, dim);
}
template void gpuCellList2D(int*, double*, double*, double*, double*, double*, int*, int, int, int);
template void gpuCellList2D(int*, float*, float*, float*, float*, float*, int*, int, int, int);

__global__ void gpuKernelCellCounts(int *cellcounts, int* cell, int natom, int ncell)
{
    int i = threadIdx.x + blockIdx.x * blockDim.x;
    while (i < natom) {              
        int cid = cell[i]; // cell id of atom i
        int j;

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

        i += blockDim.x * gridDim.x;
    }
}
        
void gpuCellCounts(int *cellcounts, int* cell, int natom, int ncell)
{        
    int blockDim = 256;
    int gridDim = (natom + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;
    gpuKernelCellCounts<<<gridDim, blockDim>>>(cellcounts, cell, natom, ncell);
}

void gpuCell2AtomList(int *c2alist, int *c2anumsum, int *c2anum, int *clist, int *d_temp, int natom, int ncell)
{           
    // sort atom-to-cell list    
    size_t  temp_storage_bytes  = 0;
    void  *d_temp_storage = NULL;
    cub::DeviceRadixSort::SortPairs(d_temp_storage, temp_storage_bytes, clist, c2anum, c2anumsum, c2alist, natom);
            
    //cudaMalloc( (void**)&d_temp_storage, temp_storage_bytes);
    gpuIndexInit(c2anumsum, natom);
    cub::DeviceRadixSort::SortPairs((void*) d_temp, temp_storage_bytes, clist, c2anum, c2anumsum, c2alist, natom);

    // count number of atoms for every cell
    gpuCellCounts(c2anumsum, c2anum, natom, ncell);     
}

template <typename T> 
__global__ void gpuKernelFullNeighborList2D(int *neighborlist, int *neighbornum, T *x, T *rcutsq, int *alist, 
        int *clist, int *c2alist, int *c2anumsum, int *nc, int inum, int jnum, int dim)
{                
    int ii = threadIdx.x + blockIdx.x * blockDim.x;
    while (ii < inum) {           // for each atom i in the simulation box           
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
        ii += blockDim.x * gridDim.x;
    }        
}

template <typename T> void gpuFullNeighborList2D(int *neighborlist, int *neighbornum, T *x, T *rcutsq,
        int *alist, int *clist, int *c2alist, int *c2anumsum, int *nc, int inum, int jnum, int dim)
{                        
    int blockDim = 256;
    int gridDim = (inum + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;
    
    gpuKernelFullNeighborList2D<<<gridDim, blockDim>>>(neighborlist, neighbornum, 
            x, rcutsq, alist, clist, c2alist, c2anumsum, nc, inum, jnum, dim);
}
template void gpuFullNeighborList2D(int*, int*, double*, double*, int*, int*, int*, int*, int*, int, int, int);
template void gpuFullNeighborList2D(int*, int*, float*, float*, int*, int*, int*, int*, int*, int, int, int);

template <typename T> 
__global__ void gpuKernelFullNeighborList2D(int *neighborlist, int *neighbornum, T *x, T *rcutsq, 
        int anum, int inum, int jnum, int dim)
{                
    int i = threadIdx.x + blockIdx.x * blockDim.x;
    while (i < inum) {          // for each atom ii in the simulation box    
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
        i += blockDim.x * gridDim.x;
    }        
}

template <typename T> void gpuFullNeighborList2D(int *neighborlist, int *neighbornum, T *x, T *rcutsq,
        int anum, int inum, int jnum, int dim)
{                        
    int blockDim = 256;
    int gridDim = (inum + blockDim - 1) / blockDim;
    gridDim = (gridDim>1024)? 1024 : gridDim;
    
    gpuKernelFullNeighborList2D<<<gridDim, blockDim>>>(neighborlist, neighbornum, 
            x, rcutsq, anum, inum, jnum, dim);
}
template void gpuFullNeighborList2D(int*, int*, double*, double*, int, int, int, int);
template void gpuFullNeighborList2D(int*, int*, float*, float*, int, int, int, int);


#endif


