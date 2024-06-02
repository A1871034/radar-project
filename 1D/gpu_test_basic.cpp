#include <stdio.h>
#include <math.h>
#include <iostream>

#include "omp.h"

using namespace std;

int main()
{
    int time = 100;
    int size = 200;
    int scalar = 0;

    double * ss = (double *) calloc(size, sizeof(double));

    /* do time stepping */
    int * gpu_step = (int *) omp_target_alloc(sizeof(int), 1);
    #pragma omp target data map(tofrom: ss[0:size], scalar)
    {
        
        for (int i = 0; i < time; i++) {
            #pragma omp target teams distribute parallel for simd
            for (int j = 0; j < size-1; j++) {
                ss[j]++;
            }
            #pragma omp target update to(time)
            #pragma omp single
            {
                ss[size-1] = time;
            }
        }
    }
    omp_target_free(gpu_step, 1);
    printf("scalar: %d\n", scalar);

    for (int i = 0; i < size; i++) {
        if (ss[i] != time) {
            printf("%d:%f\n", i, ss[i]);
        }
    }

    free(ss);

    return 0;
}