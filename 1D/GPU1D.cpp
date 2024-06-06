#include <stdio.h>
#include <math.h>
#include <iostream>
#include <fstream>
#include <sstream>

#include "omp.h"

#define SIZE 20000000

using namespace std;

int main()
{
    float * ez = (float*) calloc(SIZE, sizeof(float));
    float * hy = (float*) calloc(SIZE, sizeof(float));
    float imp0 = 377.0;
    int qTime, maxTime = 10000, mm, sample_period = 10;

    char * ss = (char *) malloc((maxTime/sample_period + 1)*sizeof(char));

    /* do time stepping */
    #pragma omp target teams map(ez[0:SIZE], hy[0:SIZE], imp0, qTime, maxTime, mm)
    {
        for (qTime = 0; qTime < maxTime; qTime++) {
            /* simple ABC for hy[Size-1] */
            hy[SIZE-1] = hy[SIZE-2];

            /* update magnetic field */
            #pragma omp distribute parallel for
            for (mm = 0; mm < SIZE - 1; mm++)
            hy[mm] = hy[mm] + (ez[mm + 1] - ez[mm]) / imp0;

            /* simple ABC for ez[0] */
            ez[0] = ez[1];

            /* update electric field */
            #pragma omp distribute parallel for 
            for (mm = 1; mm < SIZE; mm++)
            ez[mm] = ez[mm] + (hy[mm] - hy[mm - 1]) * imp0;

            /* hardwire a source node */
            ez[67] += exp(-(qTime - 30.) * (qTime - 30.) / 100.);
            ez[133] += exp(-(qTime - 30.) * (qTime - 30.) / 100.);

            hy[49] -= exp(-(qTime - 200.) * (qTime - 200.) / 100.) / imp0;
            ez[50] += exp(-(qTime - 200.) * (qTime - 200.) / 100.);
            ez[150] -= exp(-(qTime - 200.) * (qTime - 200.) / 100.);
            hy[150] -= exp(-(qTime - 200.) * (qTime - 200.) / 100.) / imp0;

            ez[50] += exp(-(qTime - 400.) * (qTime - 400.) / 100.);
            ez[100] -= exp(-(qTime - 400.) * (qTime - 400.) / 100.);
            ez[150] += exp(-(qTime - 400.) * (qTime - 400.) / 100.);
            
        } /* end of time-stepping */
    }

    free(ez);
    free(hy);
    free(ss);

    return 0;
}