/***************************************************************************                               
                    Molecular Dynamics Potentials (MDP)
                           CESMIX-MIT Project  
 
 Contributing authors: Ngoc-Cuong Nguyen (cuongng@mit.edu, exapde@gmail.com)
 ***************************************************************************/

#ifndef __HIPSORT
#define __HIPSORT

#include "hipcub/device/device_radix_sort.hpp"
#include "hipcub/device/device_run_length_encode.hpp"
#include "hipcub/device/device_scan.hpp"
#include "hipcub/device/device_reduce.hpp"

template <typename T> T hipArraySum(T *a, T *b, int n)
{       
    size_t  temp_storage_bytes  = 0;
    void  *d_temp_storage = NULL;
    hipcub::DeviceReduce::Sum(d_temp_storage, temp_storage_bytes, a, b, n);                
    hipcub::DeviceReduce::Sum((void*) &b[1], temp_storage_bytes, a, b, n);

    return hipArrayGetValueAtIndex(b, 0);    
}
template int hipArraySum(int*, int*,  int);
template float hipArraySum(float*, float*,  int);
template double hipArraySum(double*, double*, int);

//     static hipError_t Sum(
//         void                        *d_temp_storage,                    ///< [in] %Device-accessible allocation of temporary storage.  When NULL, the required allocation size is written to \p temp_storage_bytes and no work is done.
//         size_t                      &temp_storage_bytes,                ///< [in,out] Reference to size in bytes of \p d_temp_storage allocation
//         InputIteratorT              d_in,                               ///< [in] Pointer to the input sequence of data items
//         OutputIteratorT             d_out,                              ///< [out] Pointer to the output aggregate
//         int                         num_items,                          ///< [in] Total number of input items (i.e., length of \p d_in)
//         hipStream_t                stream              = 0,            ///< [in] <b>[optional]</b> CUDA stream to launch kernels within.  Default is stream<sub>0</sub>.
//         bool                        debug_synchronous   = false)        ///< [in] <b>[optional]</b> Whether or not to synchronize the stream after every kernel launch to check for errors.  Also causes launch configurations to be printed to the console.  Default is \p false.


int hipUniqueSort(int *b, int *c, int *d, int *e, int *a, int *p, int *t, int *q, int n)
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
    
    // p = [0,1,2,...,n]        
    hipIndexInit(p, n);
            
    // sort array a
    // e is the sorted array
    // d is the sorted index        
    size_t  temp_storage_bytes  = 0;
    void  *d_temp_storage = NULL;
    hipcub::DeviceRadixSort::SortPairs(d_temp_storage, temp_storage_bytes, a, e, p, d, n);                
    hipcub::DeviceRadixSort::SortPairs((void*) q, temp_storage_bytes, a, e, p, d, n);

    // b contains unique elements of e
    // p contains the counts for unique elements                
    // t[0] is the number of unique elements            
    temp_storage_bytes  = 0;
    hipcub::DeviceRunLengthEncode::Encode(d_temp_storage, temp_storage_bytes, e, b, p, t, n);                
    hipcub::DeviceRunLengthEncode::Encode((void*) q, temp_storage_bytes, e, b, p, t, n);

    // number of unique elements
    int m = hipArrayGetValueAtIndex(t, 0);
    hipArraySetValueAtIndex(c, (int) 0, (int) 0);

    // c is inclusive scan of p
    temp_storage_bytes  = 0;
    hipcub::DeviceScan::InclusiveSum(d_temp_storage, temp_storage_bytes, p, &c[1], m);                
    hipcub::DeviceScan::InclusiveSum((void*) q, temp_storage_bytes, p, &c[1], m);
    
    return m;                          
}

#endif


