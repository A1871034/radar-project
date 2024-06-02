#include <stdio.h>
#include <math.h>
#include <iostream>
#include <fstream>
#include <sstream>

#include "omp.h"

#define SIZE 200000

using namespace std;

int main()
{
    double * ez = (double*) calloc(SIZE, sizeof(double));
    double * hy = (double*) calloc(SIZE, sizeof(double));
    double imp0 = 377.0;
    int maxTime = 200;
    int qTime = 0;

    double * ss = (double *) malloc(maxTime*sizeof(double));
    

    /* do time stepping */
    #pragma omp target data map(to: ez[0:SIZE], hy[0:SIZE]) map(tofrom: ss[0:maxTime])
    {
        for (qTime = 0; qTime < maxTime; qTime++) {            
            /* update magnetic field */
            #pragma omp target teams distribute parallel for simd
            for (int mm = 0; mm < SIZE - 1; mm++) {
                hy[mm] = hy[mm] + (ez[mm + 1] - ez[mm]) / imp0;
            }

            #pragma omp target nowait
            {
                /* simple ABC for ez[0] */
                ez[0] = ez[1];
            }

            /* update electric field */
            #pragma omp target teams distribute parallel for simd
            for (int mm = 1; mm < SIZE; mm++) {
                ez[mm] = ez[mm] + (hy[mm] - hy[mm - 1]) * imp0;
            }

            #pragma omp target map(to: qTime)
            {
                /* hardwire a source node */
                ez[67] += exp(-(qTime - 30.) * (qTime - 30.) / 100.);
                ez[133] += exp(-(qTime - 30.) * (qTime - 30.) / 100.);
                ss[qTime] = ez[67];
    
                /* simple ABC for hy[Size-1] */
                hy[SIZE-1] = hy[SIZE-2];
            }
        } /* end of time-stepping */
    }

    printf("step: %d\n", qTime);

    for (int i = 0; i < maxTime; i++) {
        if (ss[i] != 0) {
            printf("%d:%f\n", i, ss[i]);
        }
    }

    free(ez);
    free(hy);
    free(ss);

    return 0;
}