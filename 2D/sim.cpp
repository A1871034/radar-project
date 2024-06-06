#include "sim.h"

#include <cstdlib>
#include <omp.h>
#include <math.h>

// Update Equations and ABCs from https://eecs.wsu.edu/~schneidj/ufdtd/ufdtd.pdf


#include <string>
#include <fstream>
#include <iostream>
#include <stdio.h>
#include <cstdint>

sim::sim(unsigned int sizeX, unsigned int sizeY, unsigned int PPW) : SIZE_X(sizeX), SIZE_Y(sizeY), PPW(PPW) {
    hx = (double *) calloc(SIZE_X * (SIZE_Y - 1), sizeof(double)); 
    hy = (double *) calloc((SIZE_X - 1) * SIZE_Y, sizeof(double));
    ez = (double *) calloc(SIZE_X * SIZE_Y, sizeof(double)); 
    abcInit();
}

sim::~sim() {
    #if DEBUG_SAMPLE
    if (FILE_HANDLE != nullptr) {
        fclose(FILE_HANDLE);
    }
    #endif
    free(hx);
    free(hy);
    free(ez);
    free(ezTop);
    free(ezRight);
    free(ezBot);
    free(ezLeft);
}

void sim::abcInit() {
    ezTop = (double *) calloc(SIZE_X*6, sizeof(double));
    ezRight = (double *) calloc(SIZE_Y*6, sizeof(double));
    ezBot = (double *) calloc(SIZE_X*6, sizeof(double));
    ezLeft = (double *) calloc(SIZE_Y*6, sizeof(double));

    double temp0 = sqrt(CONST_E_DUE_TO_H * CONST_H_DUE_TO_E); // in case of our solely free space this just becomes CDTDS
    double temp1 = 1.0 / temp0 + 2.0 + temp0;
    ABC_COEF0 = -(1.0 / temp0 - 2.0 + temp0) / temp1;
    ABC_COEF1 = -2.0 * (temp0 - 1.0 / temp0) / temp1;
    ABC_COEF2 = 4.0 * (temp0 + 1.0 / temp0) / temp1;
}

void sim::pecInit(int *data) {
    PEC_HEIGHTS = data;
}

const double* sim::get_ez() {
    return ez;
}

void sim::run(unsigned int timesteps) {
    unsigned int SIZE_X = SIZE_X;
    unsigned int SIZE_Y = SIZE_Y;
    double * hx = hx;
    double * hy = hy;
    double * ez = ez;
    int * PEC_HEIGHTS = PEC_HEIGHTS;

    long s = 0;
    #if DEBUG_SAMPLE
        unsigned int samples = (timesteps/SAMPLE_PERIOD + 1);
        float * sample_data = (float *) calloc(SIZE_X * SIZE_Y * samples, sizeof(float));
        #pragma omp target data map(to: hx[0:SIZE_X * (SIZE_Y - 1)], \
                hy[0:(SIZE_X - 1) * SIZE_Y], \
                ez[0:SIZE_X * SIZE_Y], \
                ezTop[0:SIZE_X*6], \
                ezRight[0:SIZE_Y*6], \
                ezBot[0:SIZE_X*6], \
                ezLeft[0:SIZE_Y*6], \
                PEC_HEIGHTS[0:SIZE_X]) \
            map(tofrom: sample_data[0:SIZE_X * SIZE_Y * samples])
    #else
    #pragma omp target data map(to: hx[0:(SIZE_X*(SIZE_Y-1))], \
            hy[0:((SIZE_X-1)*SIZE_Y)], \
            ez[0:SIZE_X*SIZE_Y], \
            ezTop[0:SIZE_X*6], \
            ezRight[0:SIZE_Y*6], \
            ezBot[0:SIZE_X*6], \
            ezLeft[0:SIZE_Y*6], \
            PEC_HEIGHTS[0:SIZE_X], \
            SIZE_X, SIZE_Y)
    #endif
    {
        #pragma omp target
        { DEVICE_N = omp_get_device_num(); }
        for (unsigned int step = 0; step < timesteps; step++) {
            // updateH();
            #pragma omp target teams distribute parallel for simd simdlen(32) collapse(2)
            for (unsigned int m = 0; m < SIZE_X; m++) {
                for (unsigned int n = 0; n < SIZE_Y - 1; n++) {                    
                    //Hx(mm, nn) = Chxh(mm, nn) * Hx(mm, nn) - Chxe(mm, nn) * (Ez(mm, nn + 1) - Ez(mm, nn));
                    I_hx(m, n) = I_hx(m, n) - (CDTDS/IMP0)*(I_ez(m, n + 1) - I_ez(m, n));
                }
            }

            #pragma omp target teams distribute parallel for simd simdlen(32) collapse(2)
            for (unsigned int m = 0; m < SIZE_X - 1; m++) {
                for (unsigned int n = 0; n < SIZE_Y; n++) {
                    //Hy(mm, nn) = Chyh(mm, nn) * Hy(mm, nn) + Chye(mm, nn) * (Ez(mm + 1, nn) - Ez(mm, nn));
                    I_hy(m, n) = I_hy(m, n) + (CDTDS/IMP0)*(I_ez(m + 1, n) - I_ez(m, n));
                }
            }   

            // updateE();
            // For reasons unknown
            // Doing (SIZE_X - 1) AND (SIZE_Y - 1) in both iterations causes a memory access error
            // Individually is fine
            // This should just mean some extra work is done and overwritten by ABCs
            #pragma omp target teams distribute parallel for simd simdlen(32) collapse(2)
            for (unsigned int m = 1; m < SIZE_X - 1; m++) {
                for (unsigned int n = 1; n < SIZE_Y; n++) {
                    // Ez(mm, nn) = Ceze(mm, nn) * Ez(mm, nn) +
                    // Cezh(mm, nn) * ((Hy(mm, nn) - Hy(mm - 1, nn)) -
                    // (Hx(mm, nn) - Hx(mm, nn - 1)));
                    I_ez(m, n) = I_ez(m, n) + (CDTDS*IMP0) * ((I_hy(m, n) - I_hy(m - 1, n)) - (I_hx(m, n) - I_hx(m, n - 1)));
                }
            }

            // pec();
            //#pragma omp target teams distribute simd
            //for (unsigned int i = 1; i < SIZE_X - 1; i++) {
            //    I_ez(PEC_HEIGHTS[i], i) = 0.0;
            //}
            // causing segfaults on I_ez accessing even just ez[0] = 0.0;
            #pragma omp target update from(ez[0:SIZE_X*SIZE_Y])
            //// abc();
            ///* ABC at left/right side of grid */
            for (unsigned int n = 0; n < SIZE_Y; n++) {
                I_ez(0, n) = ABC_COEF0 * (I_ez(2, n) + I_ezLeft(0, 1, n))
                    + ABC_COEF1 * (I_ezLeft(0, 0, n) + I_ezLeft(2, 0, n)
                    - I_ez(1, n) - I_ezLeft(1, 1, n))
                    + ABC_COEF2 * I_ezLeft(1, 0, n) - I_ezLeft(2, 1, n);
                /* memorize old fields */
                for (unsigned int m = 0; m < 3; m++) {
                    I_ezLeft(m, 1, n) = I_ezLeft(m, 0, n);
                    I_ezLeft(m, 0, n) = I_ez(m, n);
                }
                I_ez(SIZE_X - 1, n) = ABC_COEF0 * (I_ez(SIZE_X - 3, n) + I_ezRight(0, 1, n))
                + ABC_COEF1 * (I_ezRight(0, 0, n) + I_ezRight(2, 0, n)
                - I_ez(SIZE_X - 2, n) - I_ezRight(1, 1, n))
                + ABC_COEF2 * I_ezRight(1, 0, n) - I_ezRight(2, 1, n);
                /* memorize old fields */
                for (unsigned int m = 0; m < 3; m++) {
                    I_ezRight(m, 1, n) = I_ezRight(m, 0, n);
                    I_ezRight(m, 0, n) = I_ez(SIZE_X - 1 - m, n);
                }
            }
            ///* ABC at bottom/top of grid */
            //#pragma omp target teams distribute parallel for nowait
            for (unsigned int m = 0; m < SIZE_X; m++) {
                I_ez(m, 0) = ABC_COEF0 * (I_ez(m, 2) + I_ezBot(0, 1, m))
                + ABC_COEF1 * (I_ezBot(0, 0, m) + I_ezBot(2, 0, m)
                - I_ez(m, 1) - I_ezBot(1, 1, m))
                + ABC_COEF2 * I_ezBot(1, 0, m) - I_ezBot(2, 1, m);
                /* memorize old fields */
                for (unsigned int n = 0; n < 3; n++) {
                    I_ezBot(n, 1, m) = I_ezBot(n, 0, m);
                    I_ezBot(n, 0, m) = I_ez(m, n);
                }
                I_ez(m, SIZE_Y - 1) = ABC_COEF0 * (I_ez(m, SIZE_Y - 3) + I_ezTop(0, 1, m))
                + ABC_COEF1 * (I_ezTop(0, 0, m) + I_ezTop(2, 0, m)
                - I_ez(m, SIZE_Y - 2) - I_ezTop(1, 1, m))
                + ABC_COEF2 * I_ezTop(1, 0, m) - I_ezTop(2, 1, m);
                /* memorize old fields */
                for (unsigned int n = 0; n < 3; n++) {
                        I_ezTop(n, 1, m) = I_ezTop(n, 0, m);
                        I_ezTop(n, 0, m) = I_ez(m, SIZE_Y - 1 - n);
                }
            }

            //// Hardwired additive sourcenode in center
            // Also somehow causing segfault
            //double arg = M_PI * ((CDTDS * step - 0.0) / PPW - 1.0);
            //arg = arg * arg;
            //I_ez(100, 150) += (1.0 - 2.0 * arg) * exp(-arg);
            #pragma omp target update to(ez[0:SIZE_X*SIZE_Y])

            // Debug sample
            #if DEBUG_SAMPLE
            if (step % SAMPLE_PERIOD == 0) {
                unsigned int step_offset = SIZE_X * SIZE_Y * (unsigned int) (step/SAMPLE_PERIOD);
                #pragma omp target teams distribute parallel for simd simdlen(32) collapse(2) map(to: step_offset)
                for (int i = 0; i < SIZE_X; i++) {
                    for (int j = 0; j < SIZE_Y; j++) {
                        sample_data[step_offset + ((i) * SIZE_Y + j)] = (float) I_ez(i, j);
                    }
                }
            }
            #endif
        }
        
    }
    printf("%ld\n", s);

    #if DEBUG_SAMPLE
    FILE_HANDLE = fopen("sim.debug", "w");
    fprintf(FILE_HANDLE, "%u,%u,%lu\n", SIZE_X, SIZE_Y, sizeof(float));
    fclose(FILE_HANDLE);
    printf("Writing %u floats\n", SIZE_X * SIZE_Y * samples);
        FILE_HANDLE = fopen("sim.debug", "ab");
        fwrite(sample_data, sizeof(float), SIZE_X * SIZE_Y * samples, FILE_HANDLE);
    fclose(FILE_HANDLE);
        free(sample_data);
    FILE_HANDLE = nullptr;
    #endif
}

int sim::getNumThreads() {
    return N_THREADS;
}

int sim::getDeviceNum() {
    return DEVICE_N;
}

unsigned int sim::getSizeX() {
    return SIZE_X;
}

unsigned int sim::getSizeY() {
    return SIZE_Y;
}

#if DEBUG_SAMPLE
void sim::setSamplePeriod(unsigned int period) {
    SAMPLE_PERIOD = period;
}
#endif