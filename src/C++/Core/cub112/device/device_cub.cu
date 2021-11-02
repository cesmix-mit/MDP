
/******************************************************************************
 * Copyright (c) 2011, Duane Merrill.  All rights reserved.
 * Copyright (c) 2011-2016, NVIDIA CORPORATION.  All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *     * Redistributions of source code must retain the above copyright
 *       notice, this list of conditions and the following disclaimer.
 *     * Redistributions in binary form must reproduce the above copyright
 *       notice, this list of conditions and the following disclaimer in the
 *       documentation and/or other materials provided with the distribution.
 *     * Neither the name of the NVIDIA CORPORATION nor the
 *       names of its contributors may be used to endorse or promote products
 *       derived from this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL NVIDIA CORPORATION BE LIABLE FOR ANY
 * DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 ******************************************************************************/

//!nvcc -std=c++11 -D_FORCE_INLINES -O3 -c --compiler-options '-fPIC' -w device_cub.cu -o cub.o

/**
 * \file
 * cub::DeviceRadixSort provides device-wide, parallel operations for computing a radix sort across a sequence of data items residing within device-accessible memory.
 */

#pragma once

#include <stdio.h>
#include <iterator>

#include "dispatch/dispatch_radix_sort.cuh"
#include "dispatch/dispatch_scan.cuh"
#include "dispatch/dispatch_reduce.cuh"
#include "dispatch/dispatch_rle.cuh"
#include "dispatch/dispatch_reduce_by_key.cuh"        
#include "dispatch/dispatch_select_if.cuh"        
#include "dispatch/dispatch_spmv_orig.cuh"

template <
    typename            KeyT,
    typename            ValueT>
    static void         RadixSortPairs(
    void                *d_temp_storage,                        ///< [in] %Device-accessible allocation of temporary storage.  When NULL, the required allocation size is written to \p temp_storage_bytes and no work is done.
    size_t              &temp_storage_bytes,                    ///< [in,out] Reference to size in bytes of \p d_temp_storage allocation
    KeyT          *d_keys_in,                             ///< [in] Pointer to the input data of key data to sort
    KeyT                *d_keys_out,                            ///< [out] Pointer to the sorted output sequence of key data
    ValueT        *d_values_in,                           ///< [in] Pointer to the corresponding input sequence of associated value items
    ValueT              *d_values_out,                          ///< [out] Pointer to the correspondingly-reordered output sequence of associated value items
    int                 num_items,                              ///< [in] Number of items to sort
    int                 begin_bit           = 0,                ///< [in] <b>[optional]</b> The least-significant bit index (inclusive)  needed for key comparison
    int                 end_bit             = sizeof(KeyT) * 8, ///< [in] <b>[optional]</b> The most-significant bit index (exclusive) needed for key comparison (e.g., sizeof(unsigned int) * 8)
    cudaStream_t        stream              = 0,                ///< [in] <b>[optional]</b> CUDA stream to launch kernels within.  Default is stream<sub>0</sub>.
    bool                debug_synchronous   = false)            ///< [in] <b>[optional]</b> Whether or not to synchronize the stream after every kernel launch to check for errors.  Also causes launch configurations to be printed to the console.  Default is \p false.
{
    // Signed integer type for global offsets
    typedef int OffsetT;

    cub::DoubleBuffer<KeyT>       d_keys(const_cast<KeyT*>(d_keys_in), d_keys_out);
    cub::DoubleBuffer<ValueT>     d_values(const_cast<ValueT*>(d_values_in), d_values_out);

    cub::DispatchRadixSort<false, KeyT, ValueT, OffsetT>::Dispatch(
        d_temp_storage,
        temp_storage_bytes,
        d_keys,
        d_values,
        num_items,
        begin_bit,
        end_bit,
        false,
        stream,
        debug_synchronous);
};
//template void RadixSortPairs(void*, size_t&, float*, float*, int*, int*, int, int, int, cudaStream_t, bool);

template <
    typename            KeyT,
    typename            ValueT>
static void RadixSortPairs(
    void                    *d_temp_storage,                        ///< [in] %Device-accessible allocation of temporary storage.  When NULL, the required allocation size is written to \p temp_storage_bytes and no work is done.
    size_t                  &temp_storage_bytes,                    ///< [in,out] Reference to size in bytes of \p d_temp_storage allocation
    cub::DoubleBuffer<KeyT>      &d_keys,                                ///< [in,out] Reference to the double-buffer of keys whose "current" device-accessible buffer contains the unsorted input keys and, upon return, is updated to point to the sorted output keys
    cub::DoubleBuffer<ValueT>    &d_values,                              ///< [in,out] Double-buffer of values whose "current" device-accessible buffer contains the unsorted input values and, upon return, is updated to point to the sorted output values
    int                     num_items,                              ///< [in] Number of items to sort
    int                     begin_bit           = 0,                ///< [in] <b>[optional]</b> The least-significant bit index (inclusive)  needed for key comparison
    int                     end_bit             = sizeof(KeyT) * 8, ///< [in] <b>[optional]</b> The most-significant bit index (exclusive) needed for key comparison (e.g., sizeof(unsigned int) * 8)
    cudaStream_t            stream              = 0,                ///< [in] <b>[optional]</b> CUDA stream to launch kernels within.  Default is stream<sub>0</sub>.
    bool                    debug_synchronous   = false)            ///< [in] <b>[optional]</b> Whether or not to synchronize the stream after every kernel launch to check for errors.  Also causes launch configurations to be printed to the console.  Default is \p false.
{
    // Signed integer type for global offsets
    typedef int OffsetT;

    cub::DispatchRadixSort<false, KeyT, ValueT, OffsetT>::Dispatch(
        d_temp_storage,
        temp_storage_bytes,
        d_keys,
        d_values,
        num_items,
        begin_bit,
        end_bit,
        true,
        stream,
        debug_synchronous);
}
        
template <
    typename            KeyT,
    typename            ValueT>
    static void         RadixSortPairsDescending(
    void                *d_temp_storage,                        ///< [in] %Device-accessible allocation of temporary storage.  When NULL, the required allocation size is written to \p temp_storage_bytes and no work is done.
    size_t              &temp_storage_bytes,                    ///< [in,out] Reference to size in bytes of \p d_temp_storage allocation
    const KeyT          *d_keys_in,                             ///< [in] Pointer to the input data of key data to sort
    KeyT                *d_keys_out,                            ///< [out] Pointer to the sorted output sequence of key data
    const ValueT        *d_values_in,                           ///< [in] Pointer to the corresponding input sequence of associated value items
    ValueT              *d_values_out,                          ///< [out] Pointer to the correspondingly-reordered output sequence of associated value items
    int                 num_items,                              ///< [in] Number of items to sort
    int                 begin_bit           = 0,                ///< [in] <b>[optional]</b> The least-significant bit index (inclusive)  needed for key comparison
    int                 end_bit             = sizeof(KeyT) * 8, ///< [in] <b>[optional]</b> The most-significant bit index (exclusive) needed for key comparison (e.g., sizeof(unsigned int) * 8)
    cudaStream_t        stream              = 0,                ///< [in] <b>[optional]</b> CUDA stream to launch kernels within.  Default is stream<sub>0</sub>.
    bool                debug_synchronous   = false)            ///< [in] <b>[optional]</b> Whether or not to synchronize the stream after every kernel launch to check for errors.  Also causes launch configurations to be printed to the console.  Default is \p false.
{
    // Signed integer type for global offsets
    typedef int OffsetT;

    cub::DoubleBuffer<KeyT>       d_keys(const_cast<KeyT*>(d_keys_in), d_keys_out);
    cub::DoubleBuffer<ValueT>     d_values(const_cast<ValueT*>(d_values_in), d_values_out);

    cub::DispatchRadixSort<true, KeyT, ValueT, OffsetT>::Dispatch(
        d_temp_storage,
        temp_storage_bytes,
        d_keys,
        d_values,
        num_items,
        begin_bit,
        end_bit,
        false,
        stream,
        debug_synchronous);
}
//template void RadixSortPairsDescending(void*, size_t&, float*, float*, int*, int*, int, int, int, cudaStream_t, bool);        

template <
    typename            KeyT,
    typename            ValueT>
static void RadixSortPairsDescending(
    void                    *d_temp_storage,                        ///< [in] %Device-accessible allocation of temporary storage.  When NULL, the required allocation size is written to \p temp_storage_bytes and no work is done.
    size_t                  &temp_storage_bytes,                    ///< [in,out] Reference to size in bytes of \p d_temp_storage allocation
    cub::DoubleBuffer<KeyT>      &d_keys,                                ///< [in,out] Reference to the double-buffer of keys whose "current" device-accessible buffer contains the unsorted input keys and, upon return, is updated to point to the sorted output keys
    cub::DoubleBuffer<ValueT>    &d_values,                              ///< [in,out] Double-buffer of values whose "current" device-accessible buffer contains the unsorted input values and, upon return, is updated to point to the sorted output values
    int                     num_items,                              ///< [in] Number of items to sort
    int                     begin_bit           = 0,                ///< [in] <b>[optional]</b> The least-significant bit index (inclusive)  needed for key comparison
    int                     end_bit             = sizeof(KeyT) * 8, ///< [in] <b>[optional]</b> The most-significant bit index (exclusive) needed for key comparison (e.g., sizeof(unsigned int) * 8)
    cudaStream_t            stream              = 0,                ///< [in] <b>[optional]</b> CUDA stream to launch kernels within.  Default is stream<sub>0</sub>.
    bool                    debug_synchronous   = false)            ///< [in] <b>[optional]</b> Whether or not to synchronize the stream after every kernel launch to check for errors.  Also causes launch configurations to be printed to the console.  Default is \p false.
{
    // Signed integer type for global offsets
    typedef int OffsetT;

    cub::DispatchRadixSort<true, KeyT, ValueT, OffsetT>::Dispatch(
        d_temp_storage,
        temp_storage_bytes,
        d_keys,
        d_values,
        num_items,
        begin_bit,
        end_bit,
        true,
        stream,
        debug_synchronous);
}
        
template <typename KeyT>
    static void         RadixSortKeys(
    void                *d_temp_storage,                        ///< [in] %Device-accessible allocation of temporary storage.  When NULL, the required allocation size is written to \p temp_storage_bytes and no work is done.
    size_t              &temp_storage_bytes,                    ///< [in,out] Reference to size in bytes of \p d_temp_storage allocation
    const KeyT          *d_keys_in,                             ///< [in] Pointer to the input data of key data to sort
    KeyT                *d_keys_out,                            ///< [out] Pointer to the sorted output sequence of key data
    int                 num_items,                              ///< [in] Number of items to sort
    int                 begin_bit           = 0,                ///< [in] <b>[optional]</b> The least-significant bit index (inclusive)  needed for key comparison
    int                 end_bit             = sizeof(KeyT) * 8, ///< [in] <b>[optional]</b> The most-significant bit index (exclusive) needed for key comparison (e.g., sizeof(unsigned int) * 8)
    cudaStream_t        stream              = 0,                ///< [in] <b>[optional]</b> CUDA stream to launch kernels within.  Default is stream<sub>0</sub>.
    bool                debug_synchronous   = false)            ///< [in] <b>[optional]</b> Whether or not to synchronize the stream after every kernel launch to check for errors.  Also causes launch configurations to be printed to the console.  Default is \p false.
{
    // Signed integer type for global offsets
    typedef int OffsetT;

    // Null value type
    cub::DoubleBuffer<KeyT>      d_keys(const_cast<KeyT*>(d_keys_in), d_keys_out);
    cub::DoubleBuffer<cub::NullType>  d_values;

    cub::DispatchRadixSort<false, KeyT, cub::NullType, OffsetT>::Dispatch(
        d_temp_storage,
        temp_storage_bytes,
        d_keys,
        d_values,
        num_items,
        begin_bit,
        end_bit,
        false,
        stream,
        debug_synchronous);
}

template <typename KeyT>
static void RadixSortKeys(
    void                *d_temp_storage,                        ///< [in] %Device-accessible allocation of temporary storage.  When NULL, the required allocation size is written to \p temp_storage_bytes and no work is done.
    size_t              &temp_storage_bytes,                    ///< [in,out] Reference to size in bytes of \p d_temp_storage allocation
    cub::DoubleBuffer<KeyT>  &d_keys,                                ///< [in,out] Reference to the double-buffer of keys whose "current" device-accessible buffer contains the unsorted input keys and, upon return, is updated to point to the sorted output keys
    int                 num_items,                              ///< [in] Number of items to sort
    int                 begin_bit           = 0,                ///< [in] <b>[optional]</b> The least-significant bit index (inclusive)  needed for key comparison
    int                 end_bit             = sizeof(KeyT) * 8, ///< [in] <b>[optional]</b> The most-significant bit index (exclusive) needed for key comparison (e.g., sizeof(unsigned int) * 8)
    cudaStream_t        stream              = 0,                ///< [in] <b>[optional]</b> CUDA stream to launch kernels within.  Default is stream<sub>0</sub>.
    bool                debug_synchronous   = false)            ///< [in] <b>[optional]</b> Whether or not to synchronize the stream after every kernel launch to check for errors.  Also causes launch configurations to be printed to the console.  Default is \p false.
{
    // Signed integer type for global offsets
    typedef int OffsetT;

    // Null value type
    cub::DoubleBuffer<cub::NullType> d_values;

    cub::DispatchRadixSort<false, KeyT, cub::NullType, OffsetT>::Dispatch(
        d_temp_storage,
        temp_storage_bytes,
        d_keys,
        d_values,
        num_items,
        begin_bit,
        end_bit,
        true,
        stream,
        debug_synchronous);
}

template <typename KeyT>
    static void RadixSortKeysDescending(
    void                *d_temp_storage,                        ///< [in] %Device-accessible allocation of temporary storage.  When NULL, the required allocation size is written to \p temp_storage_bytes and no work is done.
    size_t              &temp_storage_bytes,                    ///< [in,out] Reference to size in bytes of \p d_temp_storage allocation
    const KeyT          *d_keys_in,                             ///< [in] Pointer to the input data of key data to sort
    KeyT                *d_keys_out,                            ///< [out] Pointer to the sorted output sequence of key data
    int                 num_items,                              ///< [in] Number of items to sort
    int                 begin_bit           = 0,                ///< [in] <b>[optional]</b> The least-significant bit index (inclusive)  needed for key comparison
    int                 end_bit             = sizeof(KeyT) * 8, ///< [in] <b>[optional]</b> The most-significant bit index (exclusive) needed for key comparison (e.g., sizeof(unsigned int) * 8)
    cudaStream_t        stream              = 0,                ///< [in] <b>[optional]</b> CUDA stream to launch kernels within.  Default is stream<sub>0</sub>.
    bool                debug_synchronous   = false)            ///< [in] <b>[optional]</b> Whether or not to synchronize the stream after every kernel launch to check for errors.  Also causes launch configurations to be printed to the console.  Default is \p false.
{
    // Signed integer type for global offsets
    typedef int OffsetT;

    cub::DoubleBuffer<KeyT>      d_keys(const_cast<KeyT*>(d_keys_in), d_keys_out);
    cub::DoubleBuffer<cub::NullType>  d_values;

    cub::DispatchRadixSort<true, KeyT, cub::NullType, OffsetT>::Dispatch(
        d_temp_storage,
        temp_storage_bytes,
        d_keys,
        d_values,
        num_items,
        begin_bit,
        end_bit,
        false,
        stream,
        debug_synchronous);
}

template <typename KeyT>
static void RadixSortKeysDescending(
    void                *d_temp_storage,                        ///< [in] %Device-accessible allocation of temporary storage.  When NULL, the required allocation size is written to \p temp_storage_bytes and no work is done.
    size_t              &temp_storage_bytes,                    ///< [in,out] Reference to size in bytes of \p d_temp_storage allocation
    cub::DoubleBuffer<KeyT>  &d_keys,                                ///< [in,out] Reference to the double-buffer of keys whose "current" device-accessible buffer contains the unsorted input keys and, upon return, is updated to point to the sorted output keys
    int                 num_items,                              ///< [in] Number of items to sort
    int                 begin_bit           = 0,                ///< [in] <b>[optional]</b> The least-significant bit index (inclusive)  needed for key comparison
    int                 end_bit             = sizeof(KeyT) * 8, ///< [in] <b>[optional]</b> The most-significant bit index (exclusive) needed for key comparison (e.g., sizeof(unsigned int) * 8)
    cudaStream_t        stream              = 0,                ///< [in] <b>[optional]</b> CUDA stream to launch kernels within.  Default is stream<sub>0</sub>.
    bool                debug_synchronous   = false)            ///< [in] <b>[optional]</b> Whether or not to synchronize the stream after every kernel launch to check for errors.  Also causes launch configurations to be printed to the console.  Default is \p false.
{
    // Signed integer type for global offsets
    typedef int OffsetT;

    // Null value type
    cub::DoubleBuffer<cub::NullType> d_values;

    cub::DispatchRadixSort<true, KeyT, cub::NullType, OffsetT>::Dispatch(
        d_temp_storage,
        temp_storage_bytes,
        d_keys,
        d_values,
        num_items,
        begin_bit,
        end_bit,
        true,
        stream,
        debug_synchronous);
}

template <
    typename        InputIteratorT,
    typename        OutputIteratorT>
    static void     ExclusiveSum(
    void            *d_temp_storage,                    ///< [in] %Device-accessible allocation of temporary storage.  When NULL, the required allocation size is written to \p temp_storage_bytes and no work is done.
    size_t          &temp_storage_bytes,                ///< [in,out] Reference to size in bytes of \p d_temp_storage allocation
    InputIteratorT  d_in,                               ///< [in] Pointer to the input sequence of data items
    OutputIteratorT d_out,                              ///< [out] Pointer to the output sequence of data items
    int             num_items,                          ///< [in] Total number of input items (i.e., the length of \p d_in)
    cudaStream_t    stream              = 0,            ///< [in] <b>[optional]</b> CUDA stream to launch kernels within.  Default is stream<sub>0</sub>.
    bool            debug_synchronous   = false)        ///< [in] <b>[optional]</b> Whether or not to synchronize the stream after every kernel launch to check for errors.  May cause significant slowdown.  Default is \p false.
{
    // Signed integer type for global offsets
    typedef int OffsetT;

    // The output value type
    typedef typename cub::If<(cub::Equals<typename std::iterator_traits<OutputIteratorT>::value_type, void>::VALUE),  // OutputT =  (if output iterator's value type is void) ?
        typename std::iterator_traits<InputIteratorT>::value_type,                                          // ... then the input iterator's value type,
        typename std::iterator_traits<OutputIteratorT>::value_type>::Type OutputT;                          // ... else the output iterator's value type

    // Initial value
    OutputT init_value = 0;

    cub::DispatchScan<InputIteratorT, OutputIteratorT, cub::Sum, OutputT, OffsetT>::Dispatch(
        d_temp_storage,
        temp_storage_bytes,
        d_in,
        d_out,
        cub::Sum(),
        init_value,
        num_items,
        stream,
        debug_synchronous);
}

template <
    typename        InputIteratorT,
    typename        OutputIteratorT,
    typename        ScanOpT,
    typename        InitValueT>
    static void     ExclusiveScan(
    void            *d_temp_storage,                    ///< [in] %Device-accessible allocation of temporary storage.  When NULL, the required allocation size is written to \p temp_storage_bytes and no work is done.
    size_t          &temp_storage_bytes,                ///< [in,out] Reference to size in bytes of \p d_temp_storage allocation
    InputIteratorT  d_in,                               ///< [in] Pointer to the input sequence of data items
    OutputIteratorT d_out,                              ///< [out] Pointer to the output sequence of data items
    ScanOpT         scan_op,                            ///< [in] Binary scan functor
    InitValueT      init_value,                         ///< [in] Initial value to seed the exclusive scan (and is assigned to *d_out)
    int             num_items,                          ///< [in] Total number of input items (i.e., the length of \p d_in)
    cudaStream_t    stream              = 0,            ///< [in] <b>[optional]</b> CUDA stream to launch kernels within.  Default is stream<sub>0</sub>.
    bool            debug_synchronous   = false)        ///< [in] <b>[optional]</b> Whether or not to synchronize the stream after every kernel launch to check for errors.  May cause significant slowdown.  Default is \p false.
{
    // Signed integer type for global offsets
    typedef int OffsetT;

    cub::DispatchScan<InputIteratorT, OutputIteratorT, ScanOpT, InitValueT, OffsetT>::Dispatch(
        d_temp_storage,
        temp_storage_bytes,
        d_in,
        d_out,
        scan_op,
        init_value,
        num_items,
        stream,
        debug_synchronous);
}
        
template <typename            InputIteratorT, typename            OutputIteratorT>
static void         InclusiveSum(
    void*               d_temp_storage,                 ///< [in] %Device-accessible allocation of temporary storage.  When NULL, the required allocation size is written to \p temp_storage_bytes and no work is done.
    size_t&             temp_storage_bytes,             ///< [in,out] Reference to size in bytes of \p d_temp_storage allocation
    InputIteratorT      d_in,                           ///< [in] Pointer to the input sequence of data items
    OutputIteratorT     d_out,                          ///< [out] Pointer to the output sequence of data items
    int                 num_items,                      ///< [in] Total number of input items (i.e., the length of \p d_in)
    cudaStream_t        stream             = 0,         ///< [in] <b>[optional]</b> CUDA stream to launch kernels within.  Default is stream<sub>0</sub>.
    bool                debug_synchronous  = false)     ///< [in] <b>[optional]</b> Whether or not to synchronize the stream after every kernel launch to check for errors.  May cause significant slowdown.  Default is \p false.
{
    // Signed integer type for global offsets
    typedef int OffsetT;

    cub::DispatchScan<InputIteratorT, OutputIteratorT, cub::Sum, cub::NullType, OffsetT>::Dispatch(
        d_temp_storage,
        temp_storage_bytes,
        d_in,
        d_out,
        cub::Sum(),
        cub::NullType(),
        num_items,
        stream,
        debug_synchronous);
}

template <
    typename        InputIteratorT,
    typename        OutputIteratorT,
    typename        ScanOpT>
static void InclusiveScan(
    void            *d_temp_storage,                    ///< [in] %Device-accessible allocation of temporary storage.  When NULL, the required allocation size is written to \p temp_storage_bytes and no work is done.
    size_t          &temp_storage_bytes,                ///< [in,out] Reference to size in bytes of \p d_temp_storage allocation
    InputIteratorT  d_in,                               ///< [in] Pointer to the input sequence of data items
    OutputIteratorT d_out,                              ///< [out] Pointer to the output sequence of data items
    ScanOpT         scan_op,                            ///< [in] Binary scan functor
    int             num_items,                          ///< [in] Total number of input items (i.e., the length of \p d_in)
    cudaStream_t    stream             = 0,             ///< [in] <b>[optional]</b> CUDA stream to launch kernels within.  Default is stream<sub>0</sub>.
    bool            debug_synchronous  = false)         ///< [in] <b>[optional]</b> Whether or not to synchronize the stream after every kernel launch to check for errors.  May cause significant slowdown.  Default is \p false.
{
    // Signed integer type for global offsets
    typedef int OffsetT;

    cub::DispatchScan<InputIteratorT, OutputIteratorT, ScanOpT, cub::NullType, OffsetT>::Dispatch(
        d_temp_storage,
        temp_storage_bytes,
        d_in,
        d_out,
        scan_op,
        cub::NullType(),
        num_items,
        stream,
        debug_synchronous);
}

template <
    typename                    InputIteratorT,
    typename                    OutputIteratorT,
    typename                    ReductionOpT,
    typename                    T>
static void Reduce(
    void                        *d_temp_storage,                    ///< [in] %Device-accessible allocation of temporary storage.  When NULL, the required allocation size is written to \p temp_storage_bytes and no work is done.
    size_t                      &temp_storage_bytes,                ///< [in,out] Reference to size in bytes of \p d_temp_storage allocation
    InputIteratorT              d_in,                               ///< [in] Pointer to the input sequence of data items
    OutputIteratorT             d_out,                              ///< [out] Pointer to the output aggregate
    int                         num_items,                          ///< [in] Total number of input items (i.e., length of \p d_in)
    ReductionOpT                reduction_op,                       ///< [in] Binary reduction functor
    T                           init,                               ///< [in] Initial value of the reduction
    cudaStream_t                stream              = 0,            ///< [in] <b>[optional]</b> CUDA stream to launch kernels within.  Default is stream<sub>0</sub>.
    bool                        debug_synchronous   = false)        ///< [in] <b>[optional]</b> Whether or not to synchronize the stream after every kernel launch to check for errors.  Also causes launch configurations to be printed to the console.  Default is \p false.
{
    // Signed integer type for global offsets
    typedef int OffsetT;

    cub::DispatchReduce<InputIteratorT, OutputIteratorT, OffsetT, ReductionOpT>::Dispatch(
        d_temp_storage,
        temp_storage_bytes,
        d_in,
        d_out,
        num_items,
        reduction_op,
        init,
        stream,
        debug_synchronous);
}

template <
    typename                    InputIteratorT,
    typename                    OutputIteratorT>
static void Sum(
    void                        *d_temp_storage,                    ///< [in] %Device-accessible allocation of temporary storage.  When NULL, the required allocation size is written to \p temp_storage_bytes and no work is done.
    size_t                      &temp_storage_bytes,                ///< [in,out] Reference to size in bytes of \p d_temp_storage allocation
    InputIteratorT              d_in,                               ///< [in] Pointer to the input sequence of data items
    OutputIteratorT             d_out,                              ///< [out] Pointer to the output aggregate
    int                         num_items,                          ///< [in] Total number of input items (i.e., length of \p d_in)
    cudaStream_t                stream              = 0,            ///< [in] <b>[optional]</b> CUDA stream to launch kernels within.  Default is stream<sub>0</sub>.
    bool                        debug_synchronous   = false)        ///< [in] <b>[optional]</b> Whether or not to synchronize the stream after every kernel launch to check for errors.  Also causes launch configurations to be printed to the console.  Default is \p false.
{
    // Signed integer type for global offsets
    typedef int OffsetT;

    // The output value type
    typedef typename cub::If<(cub::Equals<typename std::iterator_traits<OutputIteratorT>::value_type, void>::VALUE),  // OutputT =  (if output iterator's value type is void) ?
        typename std::iterator_traits<InputIteratorT>::value_type,                                          // ... then the input iterator's value type,
        typename std::iterator_traits<OutputIteratorT>::value_type>::Type OutputT;                          // ... else the output iterator's value type

    cub::DispatchReduce<InputIteratorT, OutputIteratorT, OffsetT, cub::Sum>::Dispatch(
        d_temp_storage,
        temp_storage_bytes,
        d_in,
        d_out,
        num_items,
        cub::Sum(),
        OutputT(),            // zero-initialize
        stream,
        debug_synchronous);
}

template <
    typename                    InputIteratorT,
    typename                    OutputIteratorT>
static void Min(
    void                        *d_temp_storage,                    ///< [in] %Device-accessible allocation of temporary storage.  When NULL, the required allocation size is written to \p temp_storage_bytes and no work is done.
    size_t                      &temp_storage_bytes,                ///< [in,out] Reference to size in bytes of \p d_temp_storage allocation
    InputIteratorT              d_in,                               ///< [in] Pointer to the input sequence of data items
    OutputIteratorT             d_out,                              ///< [out] Pointer to the output aggregate
    int                         num_items,                          ///< [in] Total number of input items (i.e., length of \p d_in)
    cudaStream_t                stream              = 0,            ///< [in] <b>[optional]</b> CUDA stream to launch kernels within.  Default is stream<sub>0</sub>.
    bool                        debug_synchronous   = false)        ///< [in] <b>[optional]</b> Whether or not to synchronize the stream after every kernel launch to check for errors.  Also causes launch configurations to be printed to the console.  Default is \p false.
{
    // Signed integer type for global offsets
    typedef int OffsetT;

    // The input value type
    typedef typename std::iterator_traits<InputIteratorT>::value_type InputT;

    cub::DispatchReduce<InputIteratorT, OutputIteratorT, OffsetT, cub::Min>::Dispatch(
        d_temp_storage,
        temp_storage_bytes,
        d_in,
        d_out,
        num_items,
        cub::Min(),
        cub::Traits<InputT>::Max(), // replace with std::numeric_limits<T>::max() when C++11 support is more prevalent
        stream,
        debug_synchronous);
}

template <
    typename                    InputIteratorT,
    typename                    OutputIteratorT>
static void ArgMin(
    void                        *d_temp_storage,                    ///< [in] %Device-accessible allocation of temporary storage.  When NULL, the required allocation size is written to \p temp_storage_bytes and no work is done.
    size_t                      &temp_storage_bytes,                ///< [in,out] Reference to size in bytes of \p d_temp_storage allocation
    InputIteratorT              d_in,                               ///< [in] Pointer to the input sequence of data items
    OutputIteratorT             d_out,                              ///< [out] Pointer to the output aggregate
    int                         num_items,                          ///< [in] Total number of input items (i.e., length of \p d_in)
    cudaStream_t                stream              = 0,            ///< [in] <b>[optional]</b> CUDA stream to launch kernels within.  Default is stream<sub>0</sub>.
    bool                        debug_synchronous   = false)        ///< [in] <b>[optional]</b> Whether or not to synchronize the stream after every kernel launch to check for errors.  Also causes launch configurations to be printed to the console.  Default is \p false.
{
    // Signed integer type for global offsets
    typedef int OffsetT;

    // The input type
    typedef typename std::iterator_traits<InputIteratorT>::value_type InputValueT;

    // The output tuple type
    typedef typename cub::If<(cub::Equals<typename std::iterator_traits<OutputIteratorT>::value_type, void>::VALUE),  // OutputT =  (if output iterator's value type is void) ?
        cub::KeyValuePair<OffsetT, InputValueT>,                                                                 // ... then the key value pair OffsetT + InputValueT
        typename std::iterator_traits<OutputIteratorT>::value_type>::Type OutputTupleT;                     // ... else the output iterator's value type

    // The output value type
    typedef typename OutputTupleT::Value OutputValueT;

    // Wrapped input iterator to produce index-value <OffsetT, InputT> tuples
    typedef cub::ArgIndexInputIterator<InputIteratorT, OffsetT, OutputValueT> ArgIndexInputIteratorT;
    ArgIndexInputIteratorT d_indexed_in(d_in);

    // Initial value
    OutputTupleT initial_value(1, cub::Traits<InputValueT>::Max());   // replace with std::numeric_limits<T>::max() when C++11 support is more prevalent

    cub::DispatchReduce<ArgIndexInputIteratorT, OutputIteratorT, OffsetT, cub::ArgMin>::Dispatch(
        d_temp_storage,
        temp_storage_bytes,
        d_indexed_in,
        d_out,
        num_items,
        cub::ArgMin(),
        initial_value,
        stream,
        debug_synchronous);
}

template <
    typename                    InputIteratorT,
    typename                    OutputIteratorT>
static void Max(
    void                        *d_temp_storage,                    ///< [in] %Device-accessible allocation of temporary storage.  When NULL, the required allocation size is written to \p temp_storage_bytes and no work is done.
    size_t                      &temp_storage_bytes,                ///< [in,out] Reference to size in bytes of \p d_temp_storage allocation
    InputIteratorT              d_in,                               ///< [in] Pointer to the input sequence of data items
    OutputIteratorT             d_out,                              ///< [out] Pointer to the output aggregate
    int                         num_items,                          ///< [in] Total number of input items (i.e., length of \p d_in)
    cudaStream_t                stream              = 0,            ///< [in] <b>[optional]</b> CUDA stream to launch kernels within.  Default is stream<sub>0</sub>.
    bool                        debug_synchronous   = false)        ///< [in] <b>[optional]</b> Whether or not to synchronize the stream after every kernel launch to check for errors.  Also causes launch configurations to be printed to the console.  Default is \p false.
{
    // Signed integer type for global offsets
    typedef int OffsetT;

    // The input value type
    typedef typename std::iterator_traits<InputIteratorT>::value_type InputT;

    cub::DispatchReduce<InputIteratorT, OutputIteratorT, OffsetT, cub::Max>::Dispatch(
        d_temp_storage,
        temp_storage_bytes,
        d_in,
        d_out,
        num_items,
        cub::Max(),
        cub::Traits<InputT>::Lowest(),    // replace with std::numeric_limits<T>::lowest() when C++11 support is more prevalent
        stream,
        debug_synchronous);
}

template <
    typename                    InputIteratorT,
    typename                    OutputIteratorT>
static void ArgMax(
    void                        *d_temp_storage,                    ///< [in] %Device-accessible allocation of temporary storage.  When NULL, the required allocation size is written to \p temp_storage_bytes and no work is done.
    size_t                      &temp_storage_bytes,                ///< [in,out] Reference to size in bytes of \p d_temp_storage allocation
    InputIteratorT              d_in,                               ///< [in] Pointer to the input sequence of data items
    OutputIteratorT             d_out,                              ///< [out] Pointer to the output aggregate
    int                         num_items,                          ///< [in] Total number of input items (i.e., length of \p d_in)
    cudaStream_t                stream              = 0,            ///< [in] <b>[optional]</b> CUDA stream to launch kernels within.  Default is stream<sub>0</sub>.
    bool                        debug_synchronous   = false)        ///< [in] <b>[optional]</b> Whether or not to synchronize the stream after every kernel launch to check for errors.  Also causes launch configurations to be printed to the console.  Default is \p false.
{
    // Signed integer type for global offsets
    typedef int OffsetT;

    // The input type
    typedef typename std::iterator_traits<InputIteratorT>::value_type InputValueT;

    // The output tuple type
    typedef typename cub::If<(cub::Equals<typename std::iterator_traits<OutputIteratorT>::value_type, void>::VALUE),  // OutputT =  (if output iterator's value type is void) ?
        cub::KeyValuePair<OffsetT, InputValueT>,                                                                 // ... then the key value pair OffsetT + InputValueT
        typename std::iterator_traits<OutputIteratorT>::value_type>::Type OutputTupleT;                     // ... else the output iterator's value type

    // The output value type
    typedef typename OutputTupleT::Value OutputValueT;

    // Wrapped input iterator to produce index-value <OffsetT, InputT> tuples
    typedef cub::ArgIndexInputIterator<InputIteratorT, OffsetT, OutputValueT> ArgIndexInputIteratorT;
    ArgIndexInputIteratorT d_indexed_in(d_in);

    // Initial value
    OutputTupleT initial_value(1, cub::Traits<InputValueT>::Lowest());     // replace with std::numeric_limits<T>::lowest() when C++11 support is more prevalent

    cub::DispatchReduce<ArgIndexInputIteratorT, OutputIteratorT, OffsetT, cub::ArgMax>::Dispatch(
        d_temp_storage,
        temp_storage_bytes,
        d_indexed_in,
        d_out,
        num_items,
        cub::ArgMax(),
        initial_value,
        stream,
        debug_synchronous);
}

template <
    typename                    KeysInputIteratorT,
    typename                    UniqueOutputIteratorT,
    typename                    ValuesInputIteratorT,
    typename                    AggregatesOutputIteratorT,
    typename                    NumRunsOutputIteratorT,
    typename                    ReductionOpT>
__forceinline__
static void ReduceByKey(
    void                        *d_temp_storage,                ///< [in] %Device-accessible allocation of temporary storage.  When NULL, the required allocation size is written to \p temp_storage_bytes and no work is done.
    size_t                      &temp_storage_bytes,            ///< [in,out] Reference to size in bytes of \p d_temp_storage allocation
    KeysInputIteratorT          d_keys_in,                      ///< [in] Pointer to the input sequence of keys
    UniqueOutputIteratorT       d_unique_out,                   ///< [out] Pointer to the output sequence of unique keys (one key per run)
    ValuesInputIteratorT        d_values_in,                    ///< [in] Pointer to the input sequence of corresponding values
    AggregatesOutputIteratorT   d_aggregates_out,               ///< [out] Pointer to the output sequence of value aggregates (one aggregate per run)
    NumRunsOutputIteratorT      d_num_runs_out,                 ///< [out] Pointer to total number of runs encountered (i.e., the length of d_unique_out)
    ReductionOpT                reduction_op,                   ///< [in] Binary reduction functor
    int                         num_items,                      ///< [in] Total number of associated key+value pairs (i.e., the length of \p d_in_keys and \p d_in_values)
    cudaStream_t                stream             = 0,         ///< [in] <b>[optional]</b> CUDA stream to launch kernels within.  Default is stream<sub>0</sub>.
    bool                        debug_synchronous  = false)     ///< [in] <b>[optional]</b> Whether or not to synchronize the stream after every kernel launch to check for errors.  May cause significant slowdown.  Default is \p false.
{
    // Signed integer type for global offsets
    typedef int OffsetT;

    // FlagT iterator type (not used)

    // Selection op (not used)

    // Default == operator
    typedef cub::Equality EqualityOp;

    cub::DispatchReduceByKey<KeysInputIteratorT, UniqueOutputIteratorT, ValuesInputIteratorT, AggregatesOutputIteratorT, NumRunsOutputIteratorT, EqualityOp, ReductionOpT, OffsetT>::Dispatch(
        d_temp_storage,
        temp_storage_bytes,
        d_keys_in,
        d_unique_out,
        d_values_in,
        d_aggregates_out,
        d_num_runs_out,
        EqualityOp(),
        reduction_op,
        num_items,
        stream,
        debug_synchronous);
}

template <
    typename                    InputIteratorT,
    typename                    FlagIterator,
    typename                    OutputIteratorT,
    typename                    NumSelectedIteratorT>
 __forceinline__
static void PartitionFlag(
    void*               d_temp_storage,                ///< [in] %Device-accessible allocation of temporary storage.  When NULL, the required allocation size is written to \p temp_storage_bytes and no work is done.
    size_t                      &temp_storage_bytes,            ///< [in,out] Reference to size in bytes of \p d_temp_storage allocation
    InputIteratorT              d_in,                           ///< [in] Pointer to the input sequence of data items
    FlagIterator                d_flags,                        ///< [in] Pointer to the input sequence of selection flags
    OutputIteratorT             d_out,                          ///< [out] Pointer to the output sequence of partitioned data items
    NumSelectedIteratorT        d_num_selected_out,             ///< [out] Pointer to the output total number of items selected (i.e., the offset of the unselected partition)
    int                         num_items,                      ///< [in] Total number of items to select from
    cudaStream_t                stream             = 0,         ///< [in] <b>[optional]</b> CUDA stream to launch kernels within.  Default is stream<sub>0</sub>.
    bool                        debug_synchronous  = false)     ///< [in] <b>[optional]</b> Whether or not to synchronize the stream after every kernel launch to check for errors.  May cause significant slowdown.  Default is \p false.
{
    typedef int                     OffsetT;         // Signed integer type for global offsets
    typedef cub::NullType                SelectOp;       // Selection op (not used)
    typedef cub::NullType                EqualityOp;     // Equality operator (not used)

    cub::DispatchSelectIf<InputIteratorT, FlagIterator, OutputIteratorT, NumSelectedIteratorT, SelectOp, EqualityOp, OffsetT, true>::Dispatch(
        d_temp_storage,
        temp_storage_bytes,
        d_in,
        d_flags,
        d_out,
        d_num_selected_out,
        SelectOp(),
        EqualityOp(),
        num_items,
        stream,
        debug_synchronous);
}

template <
    typename                    InputIteratorT,
    typename                    OutputIteratorT,
    typename                    NumSelectedIteratorT,
    typename                    SelectOp>
__forceinline__
static void PartitionIf(
    void*               d_temp_storage,                ///< [in] %Device-accessible allocation of temporary storage.  When NULL, the required allocation size is written to \p temp_storage_bytes and no work is done.
    size_t                      &temp_storage_bytes,            ///< [in,out] Reference to size in bytes of \p d_temp_storage allocation
    InputIteratorT              d_in,                           ///< [in] Pointer to the input sequence of data items
    OutputIteratorT             d_out,                          ///< [out] Pointer to the output sequence of partitioned data items
    NumSelectedIteratorT        d_num_selected_out,             ///< [out] Pointer to the output total number of items selected (i.e., the offset of the unselected partition)
    int                         num_items,                      ///< [in] Total number of items to select from
    SelectOp                    select_op,                      ///< [in] Unary selection operator
    cudaStream_t                stream             = 0,         ///< [in] <b>[optional]</b> CUDA stream to launch kernels within.  Default is stream<sub>0</sub>.
    bool                        debug_synchronous  = false)     ///< [in] <b>[optional]</b> Whether or not to synchronize the stream after every kernel launch to check for errors.  May cause significant slowdown.  Default is \p false.
{
    typedef int                     OffsetT;         // Signed integer type for global offsets
    typedef cub::NullType*               FlagIterator;   // FlagT iterator type (not used)
    typedef cub::NullType                EqualityOp;     // Equality operator (not used)

    cub::DispatchSelectIf<InputIteratorT, FlagIterator, OutputIteratorT, NumSelectedIteratorT, SelectOp, EqualityOp, OffsetT, true>::Dispatch(
        d_temp_storage,
        temp_storage_bytes,
        d_in,
        NULL,
        d_out,
        d_num_selected_out,
        select_op,
        EqualityOp(),
        num_items,
        stream,
        debug_synchronous);
}

template <
    typename                    InputIteratorT,
    typename                    FlagIterator,
    typename                    OutputIteratorT,
    typename                    NumSelectedIteratorT>
__forceinline__
static void SelectFlag(
    void*               d_temp_storage,                ///< [in] %Device-accessible allocation of temporary storage.  When NULL, the required allocation size is written to \p temp_storage_bytes and no work is done.
    size_t                      &temp_storage_bytes,            ///< [in,out] Reference to size in bytes of \p d_temp_storage allocation
    InputIteratorT              d_in,                           ///< [in] Pointer to the input sequence of data items
    FlagIterator                d_flags,                        ///< [in] Pointer to the input sequence of selection flags
    OutputIteratorT             d_out,                          ///< [out] Pointer to the output sequence of selected data items
    NumSelectedIteratorT         d_num_selected_out,                 ///< [out] Pointer to the output total number of items selected (i.e., length of \p d_out)
    int                         num_items,                      ///< [in] Total number of input items (i.e., length of \p d_in)
    cudaStream_t                stream             = 0,         ///< [in] <b>[optional]</b> CUDA stream to launch kernels within.  Default is stream<sub>0</sub>.
    bool                        debug_synchronous  = false)     ///< [in] <b>[optional]</b> Whether or not to synchronize the stream after every kernel launch to check for errors.  May cause significant slowdown.  Default is \p false.
{
    typedef int                     OffsetT;         // Signed integer type for global offsets
    typedef cub::NullType                SelectOp;       // Selection op (not used)
    typedef cub::NullType                EqualityOp;     // Equality operator (not used)

    cub::DispatchSelectIf<InputIteratorT, FlagIterator, OutputIteratorT, NumSelectedIteratorT, SelectOp, EqualityOp, OffsetT, false>::Dispatch(
        d_temp_storage,
        temp_storage_bytes,
        d_in,
        d_flags,
        d_out,
        d_num_selected_out,
        SelectOp(),
        EqualityOp(),
        num_items,
        stream,
        debug_synchronous);
}

template <
    typename                    InputIteratorT,
    typename                    OutputIteratorT,
    typename                    NumSelectedIteratorT,
    typename                    SelectOp>
__forceinline__
static void SelectIf(
    void*               d_temp_storage,                ///< [in] %Device-accessible allocation of temporary storage.  When NULL, the required allocation size is written to \p temp_storage_bytes and no work is done.
    size_t                      &temp_storage_bytes,            ///< [in,out] Reference to size in bytes of \p d_temp_storage allocation
    InputIteratorT              d_in,                           ///< [in] Pointer to the input sequence of data items
    OutputIteratorT             d_out,                          ///< [out] Pointer to the output sequence of selected data items
    NumSelectedIteratorT         d_num_selected_out,                 ///< [out] Pointer to the output total number of items selected (i.e., length of \p d_out)
    int                         num_items,                      ///< [in] Total number of input items (i.e., length of \p d_in)
    SelectOp                    select_op,                      ///< [in] Unary selection operator
    cudaStream_t                stream             = 0,         ///< [in] <b>[optional]</b> CUDA stream to launch kernels within.  Default is stream<sub>0</sub>.
    bool                        debug_synchronous  = false)     ///< [in] <b>[optional]</b> Whether or not to synchronize the stream after every kernel launch to check for errors.  May cause significant slowdown.  Default is \p false.
{
    typedef int                     OffsetT;         // Signed integer type for global offsets
    typedef cub::NullType*               FlagIterator;   // FlagT iterator type (not used)
    typedef cub::NullType                EqualityOp;     // Equality operator (not used)

    cub::DispatchSelectIf<InputIteratorT, FlagIterator, OutputIteratorT, NumSelectedIteratorT, SelectOp, EqualityOp, OffsetT, false>::Dispatch(
        d_temp_storage,
        temp_storage_bytes,
        d_in,
        NULL,
        d_out,
        d_num_selected_out,
        select_op,
        EqualityOp(),
        num_items,
        stream,
        debug_synchronous);
}

template <
    typename                    InputIteratorT,
    typename                    OutputIteratorT,
    typename                    NumSelectedIteratorT>
__forceinline__
static void Unique(
    void*               d_temp_storage,                ///< [in] %Device-accessible allocation of temporary storage.  When NULL, the required allocation size is written to \p temp_storage_bytes and no work is done.
    size_t                      &temp_storage_bytes,            ///< [in,out] Reference to size in bytes of \p d_temp_storage allocation
    InputIteratorT              d_in,                           ///< [in] Pointer to the input sequence of data items
    OutputIteratorT             d_out,                          ///< [out] Pointer to the output sequence of selected data items
    NumSelectedIteratorT         d_num_selected_out,             ///< [out] Pointer to the output total number of items selected (i.e., length of \p d_out)
    int                         num_items,                      ///< [in] Total number of input items (i.e., length of \p d_in)
    cudaStream_t                stream             = 0,         ///< [in] <b>[optional]</b> CUDA stream to launch kernels within.  Default is stream<sub>0</sub>.
    bool                        debug_synchronous  = false)     ///< [in] <b>[optional]</b> Whether or not to synchronize the stream after every kernel launch to check for errors.  May cause significant slowdown.  Default is \p false.
{
    typedef int                     OffsetT;         // Signed integer type for global offsets
    typedef cub::NullType*               FlagIterator;   // FlagT iterator type (not used)
    typedef cub::NullType                SelectOp;       // Selection op (not used)
    typedef cub::Equality                EqualityOp;     // Default == operator

    cub::DispatchSelectIf<InputIteratorT, FlagIterator, OutputIteratorT, NumSelectedIteratorT, SelectOp, EqualityOp, OffsetT, false>::Dispatch(
        d_temp_storage,
        temp_storage_bytes,
        d_in,
        NULL,
        d_out,
        d_num_selected_out,
        SelectOp(),
        EqualityOp(),
        num_items,
        stream,
        debug_synchronous);
}

template <
    typename                    InputIteratorT,
    typename                    UniqueOutputIteratorT,
    typename                    LengthsOutputIteratorT,
    typename                    NumRunsOutputIteratorT>
__forceinline__
static void Encode(
    void*                       d_temp_storage,                ///< [in] %Device-accessible allocation of temporary storage.  When NULL, the required allocation size is written to \p temp_storage_bytes and no work is done.
    size_t                      &temp_storage_bytes,            ///< [in,out] Reference to size in bytes of \p d_temp_storage allocation
    InputIteratorT              d_in,                           ///< [in] Pointer to the input sequence of keys
    UniqueOutputIteratorT       d_unique_out,                   ///< [out] Pointer to the output sequence of unique keys (one key per run)
    LengthsOutputIteratorT      d_counts_out,                   ///< [out] Pointer to the output sequence of run-lengths (one count per run)
    NumRunsOutputIteratorT      d_num_runs_out,                     ///< [out] Pointer to total number of runs
    int                         num_items,                      ///< [in] Total number of associated key+value pairs (i.e., the length of \p d_in_keys and \p d_in_values)
    cudaStream_t                stream             = 0,         ///< [in] <b>[optional]</b> CUDA stream to launch kernels within.  Default is stream<sub>0</sub>.
    bool                        debug_synchronous  = false)     ///< [in] <b>[optional]</b> Whether or not to synchronize the stream after every kernel launch to check for errors.  May cause significant slowdown.  Default is \p false.
{
    typedef int         OffsetT;                    // Signed integer type for global offsets
    typedef cub::NullType*   FlagIterator;               // FlagT iterator type (not used)
    typedef cub::NullType    SelectOp;                   // Selection op (not used)
    typedef cub::Equality    EqualityOp;                 // Default == operator
    typedef cub::Sum    ReductionOp;                // Value reduction operator

    // The lengths output value type
    typedef typename cub::If<(cub::Equals<typename std::iterator_traits<LengthsOutputIteratorT>::value_type, void>::VALUE),   // LengthT =  (if output iterator's value type is void) ?
        OffsetT,                                                                                                    // ... then the OffsetT type,
        typename std::iterator_traits<LengthsOutputIteratorT>::value_type>::Type LengthT;                           // ... else the output iterator's value type

    // Generator type for providing 1s values for run-length reduction
    typedef cub::ConstantInputIterator<LengthT, OffsetT> LengthsInputIteratorT;

    cub::DispatchReduceByKey<InputIteratorT, UniqueOutputIteratorT, LengthsInputIteratorT, LengthsOutputIteratorT, NumRunsOutputIteratorT, EqualityOp, ReductionOp, OffsetT>::Dispatch(
        d_temp_storage,
        temp_storage_bytes,
        d_in,
        d_unique_out,
        LengthsInputIteratorT((LengthT) 1),
        d_counts_out,
        d_num_runs_out,
        EqualityOp(),
        ReductionOp(),
        num_items,
        stream,
        debug_synchronous);
}

template <
    typename                InputIteratorT,
    typename                OffsetsOutputIteratorT,
    typename                LengthsOutputIteratorT,
    typename                NumRunsOutputIteratorT>
__forceinline__
static void NonTrivialRuns(
    void*               d_temp_storage,                ///< [in] %Device-accessible allocation of temporary storage.  When NULL, the required allocation size is written to \p temp_storage_bytes and no work is done.
    size_t                  &temp_storage_bytes,            ///< [in,out] Reference to size in bytes of \p d_temp_storage allocation
    InputIteratorT          d_in,                           ///< [in] Pointer to input sequence of data items
    OffsetsOutputIteratorT  d_offsets_out,                  ///< [out] Pointer to output sequence of run-offsets (one offset per non-trivial run)
    LengthsOutputIteratorT  d_lengths_out,                  ///< [out] Pointer to output sequence of run-lengths (one count per non-trivial run)
    NumRunsOutputIteratorT  d_num_runs_out,                 ///< [out] Pointer to total number of runs (i.e., length of \p d_offsets_out)
    int                     num_items,                      ///< [in] Total number of associated key+value pairs (i.e., the length of \p d_in_keys and \p d_in_values)
    cudaStream_t            stream             = 0,         ///< [in] <b>[optional]</b> CUDA stream to launch kernels within.  Default is stream<sub>0</sub>.
    bool                    debug_synchronous  = false)     ///< [in] <b>[optional]</b> Whether or not to synchronize the stream after every kernel launch to check for errors.  May cause significant slowdown.  Default is \p false.
{
    typedef int         OffsetT;                    // Signed integer type for global offsets
    typedef cub::Equality    EqualityOp;                 // Default == operator

    cub::DeviceRleDispatch<InputIteratorT, OffsetsOutputIteratorT, LengthsOutputIteratorT, NumRunsOutputIteratorT, EqualityOp, OffsetT>::Dispatch(
        d_temp_storage,
        temp_storage_bytes,
        d_in,
        d_offsets_out,
        d_lengths_out,
        d_num_runs_out,
        EqualityOp(),
        num_items,
        stream,
        debug_synchronous);
}

template <
    typename            KeyT,
    typename            ValueT,
    typename            OffsetIteratorT>
static void SegmentedRadixSortPairs(
    void                *d_temp_storage,                        ///< [in] %Device-accessible allocation of temporary storage.  When NULL, the required allocation size is written to \p temp_storage_bytes and no work is done.
    size_t              &temp_storage_bytes,                    ///< [in,out] Reference to size in bytes of \p d_temp_storage allocation
    const KeyT          *d_keys_in,                             ///< [in] %Device-accessible pointer to the input data of key data to sort
    KeyT                *d_keys_out,                            ///< [out] %Device-accessible pointer to the sorted output sequence of key data
    const ValueT        *d_values_in,                           ///< [in] %Device-accessible pointer to the corresponding input sequence of associated value items
    ValueT              *d_values_out,                          ///< [out] %Device-accessible pointer to the correspondingly-reordered output sequence of associated value items
    int                 num_items,                              ///< [in] The total number of items to sort (across all segments)
    int                 num_segments,                           ///< [in] The number of segments that comprise the sorting data
    OffsetIteratorT     d_begin_offsets,                        ///< [in] Pointer to the sequence of beginning offsets of length \p num_segments, such that <tt>d_begin_offsets[i]</tt> is the first element of the <em>i</em><sup>th</sup> data segment in <tt>d_keys_*</tt> and <tt>d_values_*</tt>
    OffsetIteratorT     d_end_offsets,                          ///< [in] Pointer to the sequence of ending offsets of length \p num_segments, such that <tt>d_end_offsets[i]-1</tt> is the last element of the <em>i</em><sup>th</sup> data segment in <tt>d_keys_*</tt> and <tt>d_values_*</tt>.  If <tt>d_end_offsets[i]-1</tt> <= <tt>d_begin_offsets[i]</tt>, the <em>i</em><sup>th</sup> is considered empty.
    int                 begin_bit           = 0,                ///< [in] <b>[optional]</b> The least-significant bit index (inclusive)  needed for key comparison
    int                 end_bit             = sizeof(KeyT) * 8, ///< [in] <b>[optional]</b> The most-significant bit index (exclusive) needed for key comparison (e.g., sizeof(unsigned int) * 8)
    cudaStream_t        stream              = 0,                ///< [in] <b>[optional]</b> CUDA stream to launch kernels within.  Default is stream<sub>0</sub>.
    bool                debug_synchronous   = false)            ///< [in] <b>[optional]</b> Whether or not to synchronize the stream after every kernel launch to check for errors.  Also causes launch configurations to be printed to the console.  Default is \p false.
{
    // Signed integer type for global offsets
    typedef int OffsetT;

    cub::DoubleBuffer<KeyT>       d_keys(const_cast<KeyT*>(d_keys_in), d_keys_out);
    cub::DoubleBuffer<ValueT>     d_values(const_cast<ValueT*>(d_values_in), d_values_out);

    cub::DispatchSegmentedRadixSort<false, KeyT, ValueT, OffsetIteratorT, OffsetT>::Dispatch(
        d_temp_storage,
        temp_storage_bytes,
        d_keys,
        d_values,
        num_items,
        num_segments,
        d_begin_offsets,
        d_end_offsets,
        begin_bit,
        end_bit,
        false,
        stream,
        debug_synchronous);
}

template <
    typename                KeyT,
    typename                ValueT,
    typename                OffsetIteratorT>    
static void SegmentedRadixSortPairs(
    void                    *d_temp_storage,                        ///< [in] %Device-accessible allocation of temporary storage.  When NULL, the required allocation size is written to \p temp_storage_bytes and no work is done.
    size_t                  &temp_storage_bytes,                    ///< [in,out] Reference to size in bytes of \p d_temp_storage allocation
    cub::DoubleBuffer<KeyT>      &d_keys,                                ///< [in,out] Reference to the double-buffer of keys whose "current" device-accessible buffer contains the unsorted input keys and, upon return, is updated to point to the sorted output keys
    cub::DoubleBuffer<ValueT>    &d_values,                              ///< [in,out] Double-buffer of values whose "current" device-accessible buffer contains the unsorted input values and, upon return, is updated to point to the sorted output values
    int                     num_items,                              ///< [in] The total number of items to sort (across all segments)
    int                     num_segments,                           ///< [in] The number of segments that comprise the sorting data
    OffsetIteratorT         d_begin_offsets,                        ///< [in] Pointer to the sequence of beginning offsets of length \p num_segments, such that <tt>d_begin_offsets[i]</tt> is the first element of the <em>i</em><sup>th</sup> data segment in <tt>d_keys_*</tt> and <tt>d_values_*</tt>
    OffsetIteratorT         d_end_offsets,                          ///< [in] Pointer to the sequence of ending offsets of length \p num_segments, such that <tt>d_end_offsets[i]-1</tt> is the last element of the <em>i</em><sup>th</sup> data segment in <tt>d_keys_*</tt> and <tt>d_values_*</tt>.  If <tt>d_end_offsets[i]-1</tt> <= <tt>d_begin_offsets[i]</tt>, the <em>i</em><sup>th</sup> is considered empty.
    int                     begin_bit           = 0,                ///< [in] <b>[optional]</b> The least-significant bit index (inclusive)  needed for key comparison
    int                     end_bit             = sizeof(KeyT) * 8, ///< [in] <b>[optional]</b> The most-significant bit index (exclusive) needed for key comparison (e.g., sizeof(unsigned int) * 8)
    cudaStream_t            stream              = 0,                ///< [in] <b>[optional]</b> CUDA stream to launch kernels within.  Default is stream<sub>0</sub>.
    bool                    debug_synchronous   = false)            ///< [in] <b>[optional]</b> Whether or not to synchronize the stream after every kernel launch to check for errors.  Also causes launch configurations to be printed to the console.  Default is \p false.
{
    // Signed integer type for global offsets
    typedef int OffsetT;

    cub::DispatchSegmentedRadixSort<false, KeyT, ValueT, OffsetIteratorT, OffsetT>::Dispatch(
        d_temp_storage,
        temp_storage_bytes,
        d_keys,
        d_values,
        num_items,
        num_segments,
        d_begin_offsets,
        d_end_offsets,
        begin_bit,
        end_bit,
        true,
        stream,
        debug_synchronous);
}

template <
    typename            KeyT,
    typename            ValueT,
    typename            OffsetIteratorT>
static void SegmentedRadixSortPairsDescending(
    void                *d_temp_storage,                        ///< [in] %Device-accessible allocation of temporary storage.  When NULL, the required allocation size is written to \p temp_storage_bytes and no work is done.
    size_t              &temp_storage_bytes,                    ///< [in,out] Reference to size in bytes of \p d_temp_storage allocation
    const KeyT          *d_keys_in,                             ///< [in] %Device-accessible pointer to the input data of key data to sort
    KeyT                *d_keys_out,                            ///< [out] %Device-accessible pointer to the sorted output sequence of key data
    const ValueT        *d_values_in,                           ///< [in] %Device-accessible pointer to the corresponding input sequence of associated value items
    ValueT              *d_values_out,                          ///< [out] %Device-accessible pointer to the correspondingly-reordered output sequence of associated value items
    int                 num_items,                              ///< [in] The total number of items to sort (across all segments)
    int                 num_segments,                           ///< [in] The number of segments that comprise the sorting data
    OffsetIteratorT     d_begin_offsets,                        ///< [in] Pointer to the sequence of beginning offsets of length \p num_segments, such that <tt>d_begin_offsets[i]</tt> is the first element of the <em>i</em><sup>th</sup> data segment in <tt>d_keys_*</tt> and <tt>d_values_*</tt>
    OffsetIteratorT     d_end_offsets,                          ///< [in] Pointer to the sequence of ending offsets of length \p num_segments, such that <tt>d_end_offsets[i]-1</tt> is the last element of the <em>i</em><sup>th</sup> data segment in <tt>d_keys_*</tt> and <tt>d_values_*</tt>.  If <tt>d_end_offsets[i]-1</tt> <= <tt>d_begin_offsets[i]</tt>, the <em>i</em><sup>th</sup> is considered empty.
    int                 begin_bit           = 0,                ///< [in] <b>[optional]</b> The least-significant bit index (inclusive)  needed for key comparison
    int                 end_bit             = sizeof(KeyT) * 8, ///< [in] <b>[optional]</b> The most-significant bit index (exclusive) needed for key comparison (e.g., sizeof(unsigned int) * 8)
    cudaStream_t        stream              = 0,                ///< [in] <b>[optional]</b> CUDA stream to launch kernels within.  Default is stream<sub>0</sub>.
    bool                debug_synchronous   = false)            ///< [in] <b>[optional]</b> Whether or not to synchronize the stream after every kernel launch to check for errors.  Also causes launch configurations to be printed to the console.  Default is \p false.
{
    // Signed integer type for global offsets
    typedef int OffsetT;

    cub::DoubleBuffer<KeyT>       d_keys(const_cast<KeyT*>(d_keys_in), d_keys_out);
    cub::DoubleBuffer<ValueT>     d_values(const_cast<ValueT*>(d_values_in), d_values_out);

    cub::DispatchSegmentedRadixSort<true, KeyT, ValueT, OffsetIteratorT, OffsetT>::Dispatch(
        d_temp_storage,
        temp_storage_bytes,
        d_keys,
        d_values,
        num_items,
        num_segments,
        d_begin_offsets,
        d_end_offsets,
        begin_bit,
        end_bit,
        false,
        stream,
        debug_synchronous);
}

template <
    typename                KeyT,
    typename                ValueT,
    typename                OffsetIteratorT>
static void SegmentedRadixSortPairsDescending(
    void                    *d_temp_storage,                        ///< [in] %Device-accessible allocation of temporary storage.  When NULL, the required allocation size is written to \p temp_storage_bytes and no work is done.
    size_t                  &temp_storage_bytes,                    ///< [in,out] Reference to size in bytes of \p d_temp_storage allocation
    cub::DoubleBuffer<KeyT>      &d_keys,                                ///< [in,out] Reference to the double-buffer of keys whose "current" device-accessible buffer contains the unsorted input keys and, upon return, is updated to point to the sorted output keys
    cub::DoubleBuffer<ValueT>    &d_values,                              ///< [in,out] Double-buffer of values whose "current" device-accessible buffer contains the unsorted input values and, upon return, is updated to point to the sorted output values
    int                     num_items,                              ///< [in] The total number of items to sort (across all segments)
    int                     num_segments,                           ///< [in] The number of segments that comprise the sorting data
    OffsetIteratorT         d_begin_offsets,                        ///< [in] Pointer to the sequence of beginning offsets of length \p num_segments, such that <tt>d_begin_offsets[i]</tt> is the first element of the <em>i</em><sup>th</sup> data segment in <tt>d_keys_*</tt> and <tt>d_values_*</tt>
    OffsetIteratorT         d_end_offsets,                          ///< [in] Pointer to the sequence of ending offsets of length \p num_segments, such that <tt>d_end_offsets[i]-1</tt> is the last element of the <em>i</em><sup>th</sup> data segment in <tt>d_keys_*</tt> and <tt>d_values_*</tt>.  If <tt>d_end_offsets[i]-1</tt> <= <tt>d_begin_offsets[i]</tt>, the <em>i</em><sup>th</sup> is considered empty.
    int                     begin_bit           = 0,                ///< [in] <b>[optional]</b> The least-significant bit index (inclusive)  needed for key comparison
    int                     end_bit             = sizeof(KeyT) * 8, ///< [in] <b>[optional]</b> The most-significant bit index (exclusive) needed for key comparison (e.g., sizeof(unsigned int) * 8)
    cudaStream_t            stream              = 0,                ///< [in] <b>[optional]</b> CUDA stream to launch kernels within.  Default is stream<sub>0</sub>.
    bool                    debug_synchronous   = false)            ///< [in] <b>[optional]</b> Whether or not to synchronize the stream after every kernel launch to check for errors.  Also causes launch configurations to be printed to the console.  Default is \p false.
{
    // Signed integer type for global offsets
    typedef int OffsetT;

    cub::DispatchSegmentedRadixSort<true, KeyT, ValueT, OffsetIteratorT, OffsetT>::Dispatch(
        d_temp_storage,
        temp_storage_bytes,
        d_keys,
        d_values,
        num_items,
        num_segments,
        d_begin_offsets,
        d_end_offsets,
        begin_bit,
        end_bit,
        true,
        stream,
        debug_synchronous);
}

template <
    typename            KeyT,
    typename            OffsetIteratorT>
static void SegmentedRadixSortKeys(
    void                *d_temp_storage,                        ///< [in] %Device-accessible allocation of temporary storage.  When NULL, the required allocation size is written to \p temp_storage_bytes and no work is done.
    size_t              &temp_storage_bytes,                    ///< [in,out] Reference to size in bytes of \p d_temp_storage allocation
    const KeyT          *d_keys_in,                             ///< [in] %Device-accessible pointer to the input data of key data to sort
    KeyT                *d_keys_out,                            ///< [out] %Device-accessible pointer to the sorted output sequence of key data
    int                 num_items,                              ///< [in] The total number of items to sort (across all segments)
    int                 num_segments,                           ///< [in] The number of segments that comprise the sorting data
    OffsetIteratorT     d_begin_offsets,                        ///< [in] Pointer to the sequence of beginning offsets of length \p num_segments, such that <tt>d_begin_offsets[i]</tt> is the first element of the <em>i</em><sup>th</sup> data segment in <tt>d_keys_*</tt> and <tt>d_values_*</tt>
    OffsetIteratorT     d_end_offsets,                          ///< [in] Pointer to the sequence of ending offsets of length \p num_segments, such that <tt>d_end_offsets[i]-1</tt> is the last element of the <em>i</em><sup>th</sup> data segment in <tt>d_keys_*</tt> and <tt>d_values_*</tt>.  If <tt>d_end_offsets[i]-1</tt> <= <tt>d_begin_offsets[i]</tt>, the <em>i</em><sup>th</sup> is considered empty.
    int                 begin_bit           = 0,                ///< [in] <b>[optional]</b> The least-significant bit index (inclusive)  needed for key comparison
    int                 end_bit             = sizeof(KeyT) * 8, ///< [in] <b>[optional]</b> The most-significant bit index (exclusive) needed for key comparison (e.g., sizeof(unsigned int) * 8)
    cudaStream_t        stream              = 0,                ///< [in] <b>[optional]</b> CUDA stream to launch kernels within.  Default is stream<sub>0</sub>.
    bool                debug_synchronous   = false)            ///< [in] <b>[optional]</b> Whether or not to synchronize the stream after every kernel launch to check for errors.  Also causes launch configurations to be printed to the console.  Default is \p false.
{
    // Signed integer type for global offsets
    typedef int OffsetT;

    // Null value type
    cub::DoubleBuffer<KeyT>      d_keys(const_cast<KeyT*>(d_keys_in), d_keys_out);
    cub::DoubleBuffer<cub::NullType>  d_values;

    cub::DispatchSegmentedRadixSort<false, KeyT, cub::NullType, OffsetIteratorT, OffsetT>::Dispatch(
        d_temp_storage,
        temp_storage_bytes,
        d_keys,
        d_values,
        num_items,
        num_segments,
        d_begin_offsets,
        d_end_offsets,
        begin_bit,
        end_bit,
        false,
        stream,
        debug_synchronous);
}

template <
    typename            KeyT,
    typename            OffsetIteratorT>    
static void SegmentedRadixSortKeys(
    void                *d_temp_storage,                        ///< [in] %Device-accessible allocation of temporary storage.  When NULL, the required allocation size is written to \p temp_storage_bytes and no work is done.
    size_t              &temp_storage_bytes,                    ///< [in,out] Reference to size in bytes of \p d_temp_storage allocation
    cub::DoubleBuffer<KeyT>  &d_keys,                                ///< [in,out] Reference to the double-buffer of keys whose "current" device-accessible buffer contains the unsorted input keys and, upon return, is updated to point to the sorted output keys
    int                 num_items,                              ///< [in] The total number of items to sort (across all segments)
    int                 num_segments,                           ///< [in] The number of segments that comprise the sorting data
    OffsetIteratorT     d_begin_offsets,                        ///< [in] Pointer to the sequence of beginning offsets of length \p num_segments, such that <tt>d_begin_offsets[i]</tt> is the first element of the <em>i</em><sup>th</sup> data segment in <tt>d_keys_*</tt> and <tt>d_values_*</tt>
    OffsetIteratorT     d_end_offsets,                          ///< [in] Pointer to the sequence of ending offsets of length \p num_segments, such that <tt>d_end_offsets[i]-1</tt> is the last element of the <em>i</em><sup>th</sup> data segment in <tt>d_keys_*</tt> and <tt>d_values_*</tt>.  If <tt>d_end_offsets[i]-1</tt> <= <tt>d_begin_offsets[i]</tt>, the <em>i</em><sup>th</sup> is considered empty.
    int                 begin_bit           = 0,                ///< [in] <b>[optional]</b> The least-significant bit index (inclusive)  needed for key comparison
    int                 end_bit             = sizeof(KeyT) * 8, ///< [in] <b>[optional]</b> The most-significant bit index (exclusive) needed for key comparison (e.g., sizeof(unsigned int) * 8)
    cudaStream_t        stream              = 0,                ///< [in] <b>[optional]</b> CUDA stream to launch kernels within.  Default is stream<sub>0</sub>.
    bool                debug_synchronous   = false)            ///< [in] <b>[optional]</b> Whether or not to synchronize the stream after every kernel launch to check for errors.  Also causes launch configurations to be printed to the console.  Default is \p false.
{
    // Signed integer type for global offsets
    typedef int OffsetT;

    // Null value type
    cub::DoubleBuffer<cub::NullType> d_values;

    cub::DispatchSegmentedRadixSort<false, KeyT, cub::NullType, OffsetIteratorT, OffsetT>::Dispatch(
        d_temp_storage,
        temp_storage_bytes,
        d_keys,
        d_values,
        num_items,
        num_segments,
        d_begin_offsets,
        d_end_offsets,
        begin_bit,
        end_bit,
        true,
        stream,
        debug_synchronous);
}

template <
    typename            KeyT,
    typename            OffsetIteratorT>
static void SegmentedRadixSortKeysDescending(
    void                *d_temp_storage,                        ///< [in] %Device-accessible allocation of temporary storage.  When NULL, the required allocation size is written to \p temp_storage_bytes and no work is done.
    size_t              &temp_storage_bytes,                    ///< [in,out] Reference to size in bytes of \p d_temp_storage allocation
    const KeyT          *d_keys_in,                             ///< [in] %Device-accessible pointer to the input data of key data to sort
    KeyT                *d_keys_out,                            ///< [out] %Device-accessible pointer to the sorted output sequence of key data
    int                 num_items,                              ///< [in] The total number of items to sort (across all segments)
    int                 num_segments,                           ///< [in] The number of segments that comprise the sorting data
    OffsetIteratorT     d_begin_offsets,                        ///< [in] Pointer to the sequence of beginning offsets of length \p num_segments, such that <tt>d_begin_offsets[i]</tt> is the first element of the <em>i</em><sup>th</sup> data segment in <tt>d_keys_*</tt> and <tt>d_values_*</tt>
    OffsetIteratorT     d_end_offsets,                          ///< [in] Pointer to the sequence of ending offsets of length \p num_segments, such that <tt>d_end_offsets[i]-1</tt> is the last element of the <em>i</em><sup>th</sup> data segment in <tt>d_keys_*</tt> and <tt>d_values_*</tt>.  If <tt>d_end_offsets[i]-1</tt> <= <tt>d_begin_offsets[i]</tt>, the <em>i</em><sup>th</sup> is considered empty.
    int                 begin_bit           = 0,                ///< [in] <b>[optional]</b> The least-significant bit index (inclusive)  needed for key comparison
    int                 end_bit             = sizeof(KeyT) * 8, ///< [in] <b>[optional]</b> The most-significant bit index (exclusive) needed for key comparison (e.g., sizeof(unsigned int) * 8)
    cudaStream_t        stream              = 0,                ///< [in] <b>[optional]</b> CUDA stream to launch kernels within.  Default is stream<sub>0</sub>.
    bool                debug_synchronous   = false)            ///< [in] <b>[optional]</b> Whether or not to synchronize the stream after every kernel launch to check for errors.  Also causes launch configurations to be printed to the console.  Default is \p false.
{
    // Signed integer type for global offsets
    typedef int OffsetT;

    cub::DoubleBuffer<KeyT>      d_keys(const_cast<KeyT*>(d_keys_in), d_keys_out);
    cub::DoubleBuffer<cub::NullType>  d_values;

    cub::DispatchSegmentedRadixSort<true, KeyT, cub::NullType, OffsetIteratorT, OffsetT>::Dispatch(
        d_temp_storage,
        temp_storage_bytes,
        d_keys,
        d_values,
        num_items,
        num_segments,
        d_begin_offsets,
        d_end_offsets,
        begin_bit,
        end_bit,
        false,
        stream,
        debug_synchronous);
}

template <
    typename            KeyT,
    typename            OffsetIteratorT>
static void SegmentedRadixSortKeysDescending(
    void                *d_temp_storage,                        ///< [in] %Device-accessible allocation of temporary storage.  When NULL, the required allocation size is written to \p temp_storage_bytes and no work is done.
    size_t              &temp_storage_bytes,                    ///< [in,out] Reference to size in bytes of \p d_temp_storage allocation
    cub::DoubleBuffer<KeyT>  &d_keys,                                ///< [in,out] Reference to the double-buffer of keys whose "current" device-accessible buffer contains the unsorted input keys and, upon return, is updated to point to the sorted output keys
    int                 num_items,                              ///< [in] The total number of items to sort (across all segments)
    int                 num_segments,                           ///< [in] The number of segments that comprise the sorting data
    OffsetIteratorT     d_begin_offsets,                        ///< [in] Pointer to the sequence of beginning offsets of length \p num_segments, such that <tt>d_begin_offsets[i]</tt> is the first element of the <em>i</em><sup>th</sup> data segment in <tt>d_keys_*</tt> and <tt>d_values_*</tt>
    OffsetIteratorT     d_end_offsets,                          ///< [in] Pointer to the sequence of ending offsets of length \p num_segments, such that <tt>d_end_offsets[i]-1</tt> is the last element of the <em>i</em><sup>th</sup> data segment in <tt>d_keys_*</tt> and <tt>d_values_*</tt>.  If <tt>d_end_offsets[i]-1</tt> <= <tt>d_begin_offsets[i]</tt>, the <em>i</em><sup>th</sup> is considered empty.
    int                 begin_bit           = 0,                ///< [in] <b>[optional]</b> The least-significant bit index (inclusive)  needed for key comparison
    int                 end_bit             = sizeof(KeyT) * 8, ///< [in] <b>[optional]</b> The most-significant bit index (exclusive) needed for key comparison (e.g., sizeof(unsigned int) * 8)
    cudaStream_t        stream              = 0,                ///< [in] <b>[optional]</b> CUDA stream to launch kernels within.  Default is stream<sub>0</sub>.
    bool                debug_synchronous   = false)            ///< [in] <b>[optional]</b> Whether or not to synchronize the stream after every kernel launch to check for errors.  Also causes launch configurations to be printed to the console.  Default is \p false.
{
    // Signed integer type for global offsets
    typedef int OffsetT;

    // Null value type
    cub::DoubleBuffer<cub::NullType> d_values;

    cub::DispatchSegmentedRadixSort<true, KeyT, cub::NullType, OffsetIteratorT, OffsetT>::Dispatch(
        d_temp_storage,
        temp_storage_bytes,
        d_keys,
        d_values,
        num_items,
        num_segments,
        d_begin_offsets,
        d_end_offsets,
        begin_bit,
        end_bit,
        true,
        stream,
        debug_synchronous);
}

template <
    typename            InputIteratorT,
    typename            OutputIteratorT,
    typename            OffsetIteratorT,
    typename            ReductionOp,
    typename            T>
static void SegmentedReduce(
    void                *d_temp_storage,                    ///< [in] %Device-accessible allocation of temporary storage.  When NULL, the required allocation size is written to \p temp_storage_bytes and no work is done.
    size_t              &temp_storage_bytes,                ///< [in,out] Reference to size in bytes of \p d_temp_storage allocation
    InputIteratorT      d_in,                               ///< [in] Pointer to the input sequence of data items
    OutputIteratorT     d_out,                              ///< [out] Pointer to the output aggregate
    int                 num_segments,                       ///< [in] The number of segments that comprise the sorting data
    OffsetIteratorT     d_begin_offsets,                    ///< [in] Pointer to the sequence of beginning offsets of length \p num_segments, such that <tt>d_begin_offsets[i]</tt> is the first element of the <em>i</em><sup>th</sup> data segment in <tt>d_keys_*</tt> and <tt>d_values_*</tt>
    OffsetIteratorT     d_end_offsets,                      ///< [in] Pointer to the sequence of ending offsets of length \p num_segments, such that <tt>d_end_offsets[i]-1</tt> is the last element of the <em>i</em><sup>th</sup> data segment in <tt>d_keys_*</tt> and <tt>d_values_*</tt>.  If <tt>d_end_offsets[i]-1</tt> <= <tt>d_begin_offsets[i]</tt>, the <em>i</em><sup>th</sup> is considered empty.
    ReductionOp         reduction_op,                       ///< [in] Binary reduction functor 
    T                   initial_value,                      ///< [in] Initial value of the reduction for each segment
    cudaStream_t        stream              = 0,            ///< [in] <b>[optional]</b> CUDA stream to launch kernels within.  Default is stream<sub>0</sub>.
    bool                debug_synchronous   = false)        ///< [in] <b>[optional]</b> Whether or not to synchronize the stream after every kernel launch to check for errors.  Also causes launch configurations to be printed to the console.  Default is \p false.
{
    // Signed integer type for global offsets
    typedef int OffsetT;

    cub::DispatchSegmentedReduce<InputIteratorT, OutputIteratorT, OffsetIteratorT, OffsetT, ReductionOp>::Dispatch(
        d_temp_storage,
        temp_storage_bytes,
        d_in,
        d_out,
        num_segments,
        d_begin_offsets,
        d_end_offsets,
        reduction_op,
        initial_value,
        stream,
        debug_synchronous);
}

template <
    typename            InputIteratorT,
    typename            OutputIteratorT,
    typename            OffsetIteratorT>
static void SegmentedSum(
    void                *d_temp_storage,                    ///< [in] %Device-accessible allocation of temporary storage.  When NULL, the required allocation size is written to \p temp_storage_bytes and no work is done.
    size_t              &temp_storage_bytes,                ///< [in,out] Reference to size in bytes of \p d_temp_storage allocation
    InputIteratorT      d_in,                               ///< [in] Pointer to the input sequence of data items
    OutputIteratorT     d_out,                              ///< [out] Pointer to the output aggregate
    int                 num_segments,                       ///< [in] The number of segments that comprise the sorting data
    OffsetIteratorT     d_begin_offsets,                    ///< [in] Pointer to the sequence of beginning offsets of length \p num_segments, such that <tt>d_begin_offsets[i]</tt> is the first element of the <em>i</em><sup>th</sup> data segment in <tt>d_keys_*</tt> and <tt>d_values_*</tt>
    OffsetIteratorT     d_end_offsets,                      ///< [in] Pointer to the sequence of ending offsets of length \p num_segments, such that <tt>d_end_offsets[i]-1</tt> is the last element of the <em>i</em><sup>th</sup> data segment in <tt>d_keys_*</tt> and <tt>d_values_*</tt>.  If <tt>d_end_offsets[i]-1</tt> <= <tt>d_begin_offsets[i]</tt>, the <em>i</em><sup>th</sup> is considered empty.
    cudaStream_t        stream              = 0,            ///< [in] <b>[optional]</b> CUDA stream to launch kernels within.  Default is stream<sub>0</sub>.
    bool                debug_synchronous   = false)        ///< [in] <b>[optional]</b> Whether or not to synchronize the stream after every kernel launch to check for errors.  Also causes launch configurations to be printed to the console.  Default is \p false.
{
    // Signed integer type for global offsets
    typedef int OffsetT;

    // The output value type
    typedef typename cub::If<(cub::Equals<typename std::iterator_traits<OutputIteratorT>::value_type, void>::VALUE),  // OutputT =  (if output iterator's value type is void) ?
        typename std::iterator_traits<InputIteratorT>::value_type,                                          // ... then the input iterator's value type,
        typename std::iterator_traits<OutputIteratorT>::value_type>::Type OutputT;                          // ... else the output iterator's value type

    cub::DispatchSegmentedReduce<InputIteratorT,  OutputIteratorT, OffsetIteratorT, OffsetT, cub::Sum>::Dispatch(
        d_temp_storage,
        temp_storage_bytes,
        d_in,
        d_out,
        num_segments,
        d_begin_offsets,
        d_end_offsets,
        cub::Sum(),
        OutputT(),            // zero-initialize
        stream,
        debug_synchronous);
}

template <
    typename            InputIteratorT,
    typename            OutputIteratorT,
    typename            OffsetIteratorT>
static void SegmentedMin(
    void                *d_temp_storage,                    ///< [in] %Device-accessible allocation of temporary storage.  When NULL, the required allocation size is written to \p temp_storage_bytes and no work is done.
    size_t              &temp_storage_bytes,                ///< [in,out] Reference to size in bytes of \p d_temp_storage allocation
    InputIteratorT      d_in,                               ///< [in] Pointer to the input sequence of data items
    OutputIteratorT     d_out,                              ///< [out] Pointer to the output aggregate
    int                 num_segments,                       ///< [in] The number of segments that comprise the sorting data
    OffsetIteratorT     d_begin_offsets,                    ///< [in] Pointer to the sequence of beginning offsets of length \p num_segments, such that <tt>d_begin_offsets[i]</tt> is the first element of the <em>i</em><sup>th</sup> data segment in <tt>d_keys_*</tt> and <tt>d_values_*</tt>
    OffsetIteratorT     d_end_offsets,                      ///< [in] Pointer to the sequence of ending offsets of length \p num_segments, such that <tt>d_end_offsets[i]-1</tt> is the last element of the <em>i</em><sup>th</sup> data segment in <tt>d_keys_*</tt> and <tt>d_values_*</tt>.  If <tt>d_end_offsets[i]-1</tt> <= <tt>d_begin_offsets[i]</tt>, the <em>i</em><sup>th</sup> is considered empty.
    cudaStream_t        stream              = 0,            ///< [in] <b>[optional]</b> CUDA stream to launch kernels within.  Default is stream<sub>0</sub>.
    bool                debug_synchronous   = false)        ///< [in] <b>[optional]</b> Whether or not to synchronize the stream after every kernel launch to check for errors.  Also causes launch configurations to be printed to the console.  Default is \p false.
{
    // Signed integer type for global offsets
    typedef int OffsetT;

    // The input value type
    typedef typename std::iterator_traits<InputIteratorT>::value_type InputT;

    cub::DispatchSegmentedReduce<InputIteratorT,  OutputIteratorT, OffsetIteratorT, OffsetT, cub::Min>::Dispatch(
        d_temp_storage,
        temp_storage_bytes,
        d_in,
        d_out,
        num_segments,
        d_begin_offsets,
        d_end_offsets,
        cub::Min(),
        cub::Traits<InputT>::Max(),    // replace with std::numeric_limits<T>::max() when C++11 support is more prevalent
        stream,
        debug_synchronous);
}

template <
    typename            InputIteratorT,
    typename            OutputIteratorT,
    typename            OffsetIteratorT>
static void SegmentedArgMin(
    void                *d_temp_storage,                    ///< [in] %Device-accessible allocation of temporary storage.  When NULL, the required allocation size is written to \p temp_storage_bytes and no work is done.
    size_t              &temp_storage_bytes,                ///< [in,out] Reference to size in bytes of \p d_temp_storage allocation
    InputIteratorT      d_in,                               ///< [in] Pointer to the input sequence of data items
    OutputIteratorT     d_out,                              ///< [out] Pointer to the output aggregate
    int                 num_segments,                       ///< [in] The number of segments that comprise the sorting data
    OffsetIteratorT     d_begin_offsets,                    ///< [in] Pointer to the sequence of beginning offsets of length \p num_segments, such that <tt>d_begin_offsets[i]</tt> is the first element of the <em>i</em><sup>th</sup> data segment in <tt>d_keys_*</tt> and <tt>d_values_*</tt>
    OffsetIteratorT     d_end_offsets,                      ///< [in] Pointer to the sequence of ending offsets of length \p num_segments, such that <tt>d_end_offsets[i]-1</tt> is the last element of the <em>i</em><sup>th</sup> data segment in <tt>d_keys_*</tt> and <tt>d_values_*</tt>.  If <tt>d_end_offsets[i]-1</tt> <= <tt>d_begin_offsets[i]</tt>, the <em>i</em><sup>th</sup> is considered empty.
    cudaStream_t        stream              = 0,            ///< [in] <b>[optional]</b> CUDA stream to launch kernels within.  Default is stream<sub>0</sub>.
    bool                debug_synchronous   = false)        ///< [in] <b>[optional]</b> Whether or not to synchronize the stream after every kernel launch to check for errors.  Also causes launch configurations to be printed to the console.  Default is \p false.
{
    // Signed integer type for global offsets
    typedef int OffsetT;

    // The input type
    typedef typename std::iterator_traits<InputIteratorT>::value_type InputValueT;

    // The output tuple type
    typedef typename cub::If<(cub::Equals<typename std::iterator_traits<OutputIteratorT>::value_type, void>::VALUE),  // OutputT =  (if output iterator's value type is void) ?
        cub::KeyValuePair<OffsetT, InputValueT>,                                                                 // ... then the key value pair OffsetT + InputValueT
        typename std::iterator_traits<OutputIteratorT>::value_type>::Type OutputTupleT;                     // ... else the output iterator's value type

    // The output value type
    typedef typename OutputTupleT::Value OutputValueT;

    // Wrapped input iterator to produce index-value <OffsetT, InputT> tuples
    typedef cub::ArgIndexInputIterator<InputIteratorT, OffsetT, OutputValueT> ArgIndexInputIteratorT;
    ArgIndexInputIteratorT d_indexed_in(d_in);

    // Initial value
    OutputTupleT initial_value(1, cub::Traits<InputValueT>::Max());   // replace with std::numeric_limits<T>::max() when C++11 support is more prevalent

    cub::DispatchSegmentedReduce<ArgIndexInputIteratorT,  OutputIteratorT, OffsetIteratorT, OffsetT, cub::ArgMin>::Dispatch(
        d_temp_storage,
        temp_storage_bytes,
        d_indexed_in,
        d_out,
        num_segments,
        d_begin_offsets,
        d_end_offsets,
        cub::ArgMin(),
        initial_value,
        stream,
        debug_synchronous);
}

template <
    typename            InputIteratorT,
    typename            OutputIteratorT,
    typename            OffsetIteratorT>
static void SegmentedMax(
    void                *d_temp_storage,                    ///< [in] %Device-accessible allocation of temporary storage.  When NULL, the required allocation size is written to \p temp_storage_bytes and no work is done.
    size_t              &temp_storage_bytes,                ///< [in,out] Reference to size in bytes of \p d_temp_storage allocation
    InputIteratorT      d_in,                               ///< [in] Pointer to the input sequence of data items
    OutputIteratorT     d_out,                              ///< [out] Pointer to the output aggregate
    int                 num_segments,                       ///< [in] The number of segments that comprise the sorting data
    OffsetIteratorT     d_begin_offsets,                    ///< [in] Pointer to the sequence of beginning offsets of length \p num_segments, such that <tt>d_begin_offsets[i]</tt> is the first element of the <em>i</em><sup>th</sup> data segment in <tt>d_keys_*</tt> and <tt>d_values_*</tt>
    OffsetIteratorT     d_end_offsets,                      ///< [in] Pointer to the sequence of ending offsets of length \p num_segments, such that <tt>d_end_offsets[i]-1</tt> is the last element of the <em>i</em><sup>th</sup> data segment in <tt>d_keys_*</tt> and <tt>d_values_*</tt>.  If <tt>d_end_offsets[i]-1</tt> <= <tt>d_begin_offsets[i]</tt>, the <em>i</em><sup>th</sup> is considered empty.
    cudaStream_t        stream              = 0,            ///< [in] <b>[optional]</b> CUDA stream to launch kernels within.  Default is stream<sub>0</sub>.
    bool                debug_synchronous   = false)        ///< [in] <b>[optional]</b> Whether or not to synchronize the stream after every kernel launch to check for errors.  Also causes launch configurations to be printed to the console.  Default is \p false.
{
    // Signed integer type for global offsets
    typedef int OffsetT;

    // The input value type
    typedef typename std::iterator_traits<InputIteratorT>::value_type InputT;

    cub::DispatchSegmentedReduce<InputIteratorT,  OutputIteratorT, OffsetIteratorT, OffsetT, cub::Max>::Dispatch(
        d_temp_storage,
        temp_storage_bytes,
        d_in,
        d_out,
        num_segments,
        d_begin_offsets,
        d_end_offsets,
        cub::Max(),
        cub::Traits<InputT>::Lowest(),    // replace with std::numeric_limits<T>::lowest() when C++11 support is more prevalent
        stream,
        debug_synchronous);
}

template <
    typename            InputIteratorT,
    typename            OutputIteratorT,
    typename            OffsetIteratorT>
static void SegmentedArgMax(
    void                *d_temp_storage,                    ///< [in] %Device-accessible allocation of temporary storage.  When NULL, the required allocation size is written to \p temp_storage_bytes and no work is done.
    size_t              &temp_storage_bytes,                ///< [in,out] Reference to size in bytes of \p d_temp_storage allocation
    InputIteratorT      d_in,                               ///< [in] Pointer to the input sequence of data items
    OutputIteratorT     d_out,                              ///< [out] Pointer to the output aggregate
    int                 num_segments,                       ///< [in] The number of segments that comprise the sorting data
    OffsetIteratorT     d_begin_offsets,                    ///< [in] Pointer to the sequence of beginning offsets of length \p num_segments, such that <tt>d_begin_offsets[i]</tt> is the first element of the <em>i</em><sup>th</sup> data segment in <tt>d_keys_*</tt> and <tt>d_values_*</tt>
    OffsetIteratorT     d_end_offsets,                      ///< [in] Pointer to the sequence of ending offsets of length \p num_segments, such that <tt>d_end_offsets[i]-1</tt> is the last element of the <em>i</em><sup>th</sup> data segment in <tt>d_keys_*</tt> and <tt>d_values_*</tt>.  If <tt>d_end_offsets[i]-1</tt> <= <tt>d_begin_offsets[i]</tt>, the <em>i</em><sup>th</sup> is considered empty.
    cudaStream_t        stream              = 0,            ///< [in] <b>[optional]</b> CUDA stream to launch kernels within.  Default is stream<sub>0</sub>.
    bool                debug_synchronous   = false)        ///< [in] <b>[optional]</b> Whether or not to synchronize the stream after every kernel launch to check for errors.  Also causes launch configurations to be printed to the console.  Default is \p false.
{
    // Signed integer type for global offsets
    typedef int OffsetT;

    // The input type
    typedef typename std::iterator_traits<InputIteratorT>::value_type InputValueT;

    // The output tuple type
    typedef typename cub::If<(cub::Equals<typename std::iterator_traits<OutputIteratorT>::value_type, void>::VALUE),  // OutputT =  (if output iterator's value type is void) ?
        cub::KeyValuePair<OffsetT, InputValueT>,                                                                 // ... then the key value pair OffsetT + InputValueT
        typename std::iterator_traits<OutputIteratorT>::value_type>::Type OutputTupleT;                     // ... else the output iterator's value type

    // The output value type
    typedef typename OutputTupleT::Value OutputValueT;

    // Wrapped input iterator to produce index-value <OffsetT, InputT> tuples
    typedef cub::ArgIndexInputIterator<InputIteratorT, OffsetT, OutputValueT> ArgIndexInputIteratorT;
    ArgIndexInputIteratorT d_indexed_in(d_in);

    // Initial value
    OutputTupleT initial_value(1, cub::Traits<InputValueT>::Lowest());     // replace with std::numeric_limits<T>::lowest() when C++11 support is more prevalent

    cub::DispatchSegmentedReduce<ArgIndexInputIteratorT, OutputIteratorT, OffsetIteratorT, OffsetT, cub::ArgMax>::Dispatch(
        d_temp_storage,
        temp_storage_bytes,
        d_indexed_in,
        d_out,
        num_segments,
        d_begin_offsets,
        d_end_offsets,
        cub::ArgMax(),
        initial_value,
        stream,
        debug_synchronous);
}

template <
    typename            ValueT>
static void CsrMV(
    void*               d_temp_storage,                     ///< [in] %Device-accessible allocation of temporary storage.  When NULL, the required allocation size is written to \p temp_storage_bytes and no work is done.
    size_t&             temp_storage_bytes,                 ///< [in,out] Reference to size in bytes of \p d_temp_storage allocation
    ValueT*             d_values,                           ///< [in] Pointer to the array of \p num_nonzeros values of the corresponding nonzero elements of matrix <b>A</b>.
    int*                d_row_offsets,                      ///< [in] Pointer to the array of \p m + 1 offsets demarcating the start of every row in \p d_column_indices and \p d_values (with the final entry being equal to \p num_nonzeros)
    int*                d_column_indices,                   ///< [in] Pointer to the array of \p num_nonzeros column-indices of the corresponding nonzero elements of matrix <b>A</b>.  (Indices are zero-valued.)
    ValueT*             d_vector_x,                         ///< [in] Pointer to the array of \p num_cols values corresponding to the dense input vector <em>x</em>
    ValueT*             d_vector_y,                         ///< [out] Pointer to the array of \p num_rows values corresponding to the dense output vector <em>y</em>
    int                 num_rows,                           ///< [in] number of rows of matrix <b>A</b>.
    int                 num_cols,                           ///< [in] number of columns of matrix <b>A</b>.
    int                 num_nonzeros,                       ///< [in] number of nonzero elements of matrix <b>A</b>.
    cudaStream_t        stream                  = 0,        ///< [in] <b>[optional]</b> CUDA stream to launch kernels within.  Default is stream<sub>0</sub>.
    bool                debug_synchronous       = false)    ///< [in] <b>[optional]</b> Whether or not to synchronize the stream after every kernel launch to check for errors.  May cause significant slowdown.  Default is \p false.
{
    cub::SpmvParams<ValueT, int> spmv_params;
    spmv_params.d_values             = d_values;
    spmv_params.d_row_end_offsets    = d_row_offsets + 1;
    spmv_params.d_column_indices     = d_column_indices;
    spmv_params.d_vector_x           = d_vector_x;
    spmv_params.d_vector_y           = d_vector_y;
    spmv_params.num_rows             = num_rows;
    spmv_params.num_cols             = num_cols;
    spmv_params.num_nonzeros         = num_nonzeros;
    spmv_params.alpha                = 1.0;
    spmv_params.beta                 = 0.0;

    cub::DispatchSpmv<ValueT, int>::Dispatch(
        d_temp_storage,
        temp_storage_bytes,
        spmv_params,
        stream,
        debug_synchronous);
}

