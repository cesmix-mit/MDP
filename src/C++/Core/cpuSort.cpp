#ifndef __CPUSORT
#define __CPUSORT
        
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


// void merge(int *sorted, int *index, int *list, int start, int mid, int end)
// {
//     int ti=start, i=start, j=mid;
//     
//     while (i<mid || j<end)
//     {
//         if (j==end) index[ti] = outpt[i++];
//         else if (i==mid) index[ti] = output[j++];
//         else if (list[i]<list[j]) index[ti] =output[i++];
//         else index[ti] = output[j++];
//         ti++;
//     }
// 
//     for (ti=start; ti<end; ti++)
//         output[ti] = index[ti];
// }
// 
// void mergesort_recur(int *list, int *sorted, int start, int end)
// {
//     if (end-start<2)
//         return;
// 
//     mergesort_recur(list, sorted, start, start + (end-start)/2);
//     mergesort_recur(list, sorted, start + (end-start)/2, end);
//     merge(list, sorted, start, start + (end-start)/2, end);
// }
// 
// int mergesort_cpu(int *list, int *sorted, int n)
// {
//     mergesort_recur(list, sorted, 0, n);
//     return 1;
// }

// void merge(int *list, int *sorted, int start, int mid, int end)
// {
//     int ti=start, i=start, j=mid;
//     while (i<mid || j<end)
//     {
//         if (j==end) sorted[ti] = list[i++];
//         else if (i==mid) sorted[ti] = list[j++];
//         else if (list[i]<list[j]) sorted[ti] = list[i++];
//         else sorted[ti] = list[j++];
//         ti++;
//     }
// 
//     for (ti=start; ti<end; ti++)
//         list[ti] = sorted[ti];
// }
// 
// void mergesort_recur(int *list, int *sorted, int start, int end)
// {
//     if (end-start<2)
//         return;
// 
//     mergesort_recur(list, sorted, start, start + (end-start)/2);
//     mergesort_recur(list, sorted, start + (end-start)/2, end);
//     merge(list, sorted, start, start + (end-start)/2, end);
// }
// 
// int cpuMergeSort(int *sorted, int *list, int n)
// {
//     mergesort_recur(list, sorted, 0, n);
//     return 1;
// }
// 
// 
// // void merge(int *list, int *sorted, int *key, int start, int mid, int end)
// // {
// // //     int i = start, j = mid + 1;
// // // 
// // //     for (int k = start; k <= end; k++) {
// // //         list[k] = key[k];
// // //     }
// // // 
// // //     for (int k = start; k <= end; k++) {
// // //         if (i > mid) {
// // //             key[k] = list[j++];
// // //         } else if (j > end) {
// // //             key[k] = list[i++];
// // //         } else if (sorted[key[i]] <= sorted[key[j]]) {
// // //             key[k] = list[i++];
// // //         } else {
// // //             key[k] = list[j++];
// // //         }
// // //     }
// // //     
// // //     int ti=start, i=start, j=mid;
// // //     while (i<mid || j<end)
// // //     {
// // //         if (j==end) sorted[ti] = list[i++];
// // //         else if (i==mid) sorted[ti] = list[j++];
// // //         else if (list[i]<list[j]) sorted[ti] = list[i++];
// // //         else sorted[ti] = list[j++];
// // //         ti++;
// // //     }
// // //     for (ti=start; ti<end; ti++)
// // //         list[ti] = sorted[ti];
// //     
// // //     int ti=start, i=start, j=mid;
// // //     while (i<mid || j<end)
// // //     {
// // //         if (j==end) {            
// // //             sorted[ti] = list[i++];
// // //             key[i] = ti;
// // //             //key[ti] = i;            
// // //         }
// // //         else if (i==mid) {
// // //             sorted[ti] = list[j++];
// // //             key[j] = ti;
// // //             //key[ti] = j;
// // //         }
// // //         else if (list[i]<list[j]) {
// // //             sorted[ti] = list[i++];
// // //             key[i] = ti;
// // //             //key[ti] = i;
// // //         }
// // //         else {
// // //             sorted[ti] = list[j++];
// // //             key[j] = ti;
// // //             //key[ti] = j;
// // //         }
// // //         ti++;
// // //     }
// // // 
// // //     for (ti=start; ti<end; ti++) {
// // //         list[ti] = sorted[ti];
// // //         //key[ti] = end-ti;
// // //     }
// // }
// // 
// // void mergesort_recur(int *list, int *sorted, int *key, int start, int end)
// // {
// //     if (end-start<2)
// //         return;
// // 
// //     mergesort_recur(list, sorted, key, start, start + (end-start)/2);
// //     mergesort_recur(list, sorted, key, start + (end-start)/2, end);
// //     merge(list, sorted, key, start, start + (end-start)/2, end);
// // }
// // 
// // int cpuMergeSort(int *sorted, int *key, int *list, int n)
// // {
// //     mergesort_recur(list, sorted, key, 0, n);
// //     return 1;
// // }

#endif


