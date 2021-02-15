#include <time.h>
#include <string>
#include <fstream>
#include <sstream>
#include <iostream>
#include <math.h>
#include <string.h>

void printiarray(int* a, int m)
{    
    for (int i=0; i<m; i++)
        printf("%i  ", a[i]);
    printf("\n");
}

#include "cpuSort.cpp"
        
int main(int argc, char const *argv[]) {

    struct timespec start, stop;
    //int n_list[] = {10, 100, 1000, 10000, 100000, 1000000, 10000000, 100000000};
    int n_list[] = {10, 20, 100, 1000, 10000, 100000, 1000000, 10000000};
    int i, j;
    for(j=0; j<8; j++){
        printf("############ LENGTH OF LIST: %d ############\n", n_list[j]);
        int *sorted = (int *) malloc(n_list[j]*sizeof(int));
        int *list = (int *) malloc(n_list[j]*sizeof(int));
        int *sorted_s = (int *) malloc(n_list[j]*sizeof(int));
        int *list_s = (int *) malloc(n_list[j]*sizeof(int));
        int *key = (int *) malloc(n_list[j]*sizeof(int));
        for(i=0; i<n_list[j]; i++){
            list[i] = rand()%10000;
            list_s[i] = list[i];
            key[i] = i;
        }
                        
//         if (n_list[j]<=20) 
//             printiarray(list_s, n_list[j]);
                
        clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &start);
        //cpuMergeSort(sorted_s, list_s, n_list[j]);
        cpuMergeSort(sorted_s, key, list_s, n_list[j]);
        //cpuMergeSort(list_s, sorted_s, key, n_list[j]);
        clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &stop);
        double result = (stop.tv_sec - start.tv_sec) * 1e3 + (stop.tv_nsec - start.tv_nsec) / 1e6;
        printf("TIME TAKEN(Sequential CPU): %fms\n", result);
        
        if (n_list[j]<=20)  {
            printiarray(list_s, n_list[j]);
            printiarray(sorted_s, n_list[j]);
            printiarray(key, n_list[j]);
        }
        
        for(i=1; i<n_list[j]; i++){
            if(sorted_s[i-1]>sorted_s[i]){
                printf("WRONG ANSWER _1\n");
                return -1;
            }
        }
        printf("CORRECT ANSWER\n");


//         clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &start);
//         mergesort(list, sorted, n_list[j]);
//         clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &stop);
//         double result = (stop.tv_sec - start.tv_sec) * 1e3 + (stop.tv_nsec - start.tv_nsec) / 1e6;
//         printf("TIME TAKEN(Parallel GPU): %fms\n", result);
// 
// 
//         clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &start);
//         mergesort_cpu(list_s, sorted_s, n_list[j]);
//         clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &stop);
//         result = (stop.tv_sec - start.tv_sec) * 1e3 + (stop.tv_nsec - start.tv_nsec) / 1e6;
//         printf("TIME TAKEN(Sequential CPU): %fms\n", result);
// 
//         for(i=1; i<n_list[j]; i++){
//             if(sorted[i-1]>sorted[i]){
//                 printf("WRONG ANSWER _1\n");
//                 return -1;
//             }
//         }
//         for(i=0; i<n_list[j]; i++){
//             if(sorted_s[i]!=sorted[i]){
//                 printf("WRONG ANSWER _2\n");
//                 printf("P:%d, S:%d, Index:%d\n", sorted[i], sorted_s[i], i);
//                 return -1;
//             }
//         }
//         printf("CORRECT ANSWER\n");

        free(sorted);
        free(list);
        free(sorted_s);
        free(list_s);
        free(key);
        printf("##################################################\n");
    }
    return 0;
}