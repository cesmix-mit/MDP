/***************************************************************************                               
                    Molecular Dynamics Potentials (MDP)
                           CESMIX-MIT Project  
 
 Contributing authors: Ngoc-Cuong Nguyen (cuongng@mit.edu, exapde@gmail.com)
 ***************************************************************************/

#ifndef __CPUSORT
#define __CPUSORT

// #include <fstream>
// #include <sstream>
// #include <iostream>
// 
// void print2d(int* a, int m, int n)
// {
//     for (int i=0; i<m; i++) {
//         for (int j=0; j<n; j++)
//             std::cout << a[j*m+i] << "   ";
//         std::cout << std::endl;
//     }
//     std::cout << std::endl;
// }


void merge(int *output, int *index, int *input, int lo, int mid, int hi) 
{
    int i = lo, j = mid + 1;

    for (int k = lo; k <= hi; k++) {
        output[k] = index[k];
    }

    for (int k = lo; k <= hi; k++) {
        if (i > mid) {
            index[k] = output[j++];
        } else if (j > hi) {
            index[k] = output[i++];
        } else if (input[output[i]] <= input[output[j]]) {
            index[k] = output[i++];
        } else {
            index[k] = output[j++];
        }
    }
}

void mergeSort(int *output, int *index, int *input, int lo, int hi) 
{
    if (hi <= lo)
        return;
    int mid = (hi + lo) / 2;
    mergeSort(output, index, input, lo, mid);
    mergeSort(output, index, input, mid + 1, hi);
    merge(output, index, input, lo, mid, hi);
}

void cpuMergeSort(int *output, int *index, int *input, int length) 
{
    for (int i = 0; i<length; i++)
        index[i] = i;
    
    mergeSort(output, index, input, 0, length-1);
    for (int i=0; i<length; i++)
        output[i] = input[index[i]];
}

void cpuMergeSortChunk(int *output, int *index, int *input, int n, int chunk){
    int chunk_id;
    for(chunk_id=0; chunk_id*chunk<=n; chunk_id++){
        int start = chunk_id * chunk, end, mid;
        if(start >= n) return;
        mid = min(start + chunk/2, n);
        end = min(start + chunk, n);
        merge(output, index, input, start, mid, end);
    }
}

int cpuUniqueElements(int *b, int *c, int *e, int *t, int *p, int *q, int n)
{   
  
// Start with array in memory:
// X = 0, 0, 5, 10, 11, 12, 0, 0, 0, 20, 30, 40
// Compute an integer array of same length as X with 0 if the corresponding entry in X is 0, 1 otherwise:
// P = 0, 0, 1, 1, 1, 1, 0, 0, 0, 1, 1, 1
// Perform a parallel prefix sum on P (this is very similar to how a reduction is done):
// P = 0, 0, 1, 2, 3, 4, 4, 4, 4, 5, 6, 7
// Read in elements of X and P together. If X non-zero, write it to the output 
//    array at the offset indicated by the corresponding element in P.        
        
    // form array t and array p 
    for (int i=0; i<n; i++) {                
        if (i < n-1) {
            if (e[i+1] - e[i] > 0) {
                t[i] = i+1;
                p[i] = 1;
            }
            else {
                t[i] = 0;
                p[i] = 0;
            }                            
        }        
        
        if (i==n-1) {
            t[i] = n;
            p[i] = 1;
        }        
    }
            
    // parallel prefix sum on p
    cpuCumsum(q, p, n+1);
        
    // remove zeros for the array t to obtain c
    for (int i=0; i<n; i++) {
        if (t[i] > 0) {
            int j = q[i+1];
            int k = t[i];
            c[j] = k;
            b[j-1] = e[k-1];
        }            
    }
        
    return q[n];    
}

int cpuUniqueSort(int *b, int *c, int *d, int *e, int *a, int *p, int *t, int *q, int n)
{       
    // a  (input): an array of n integer elements
    // b (output): unique elements of the original array a
    // c (output): the number of counts for each element of the output array b
    // d (output): sorted indices of the original array a
    // e (output): sorted elements of the original array a
    // m (output): the length of the output array b
    // Example: a = [9  0  6  2  2  3  1  5  9  3  5  6  7  4  8  9  7  2  7  1]
    //          e = [0  1  1  2  2  2  3  3  4  5  5  6  6  7  7  7  8  9  9  9]
    //          d = [1  6  19  3  4  17  5  9  13  7  10  2  11  12  16  18  14  0  8  15]
    //          b = [0  1  2  3  4  5  6  7  8  9]
    //          c = [1  2  3  2  1  2  2  3  1  3] -> c = [0  1  3  6  8  9  11  13  16  17  20]        
    //          m = 10        
    
    // sort array a
    cpuMergeSort(e, d, a, n);
            
    // make a new array of unique elements and their counts
    return cpuUniqueElements(b, c, e, p, t, q, n);                          
}

#endif


