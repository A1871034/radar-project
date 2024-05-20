#include "sim.h"

#include <cstdlib>
#include <omp.h>

// Update Equations and ABCs from https://eecs.wsu.edu/~schneidj/ufdtd/ufdtd.pdf


sim::sim(unsigned int sizeX, unsigned int sizeY) : SIZE_X(sizeX), SIZE_Y(sizeY) {
    hx = (double *) calloc(SIZE_X * (SIZE_Y - 1), sizeof(double)); 
    hy = (double *) calloc((SIZE_X - 1) * SIZE_Y, sizeof(double));
    ez = (double *) calloc(SIZE_X * SIZE_Y, sizeof(double)); 
    abcInit();
}

sim::~sim() {
    free(hx);
    free(hy);
    free(ez);
    free(ezTop);
    free(ezRight);
    free(ezBot);
    free(ezLeft);
}

void sim::setSampler(sampler * _s) {
    s = _s;
}

void sim::setSource(sources *_so) {
    source = _so;
}

void sim::abcInit() {
    ezTop = (double *) calloc(SIZE_X*6, sizeof(double));
    ezRight = (double *) calloc(SIZE_Y*6, sizeof(double));
    ezBot = (double *) calloc(SIZE_X*6, sizeof(double));
    ezLeft = (double *) calloc(SIZE_Y*6, sizeof(double));

    double temp0 = sqrt(CONST_E_DUE_TO_H * CONST_H_DUE_TO_E); // in case of our solely free space just becomes CDTDS
    double temp1 = 1.0 / temp0 + 2.0 + temp0;
    ABC_COEF0 = -(1.0 / temp0 - 2.0 + temp0) / temp1;
    ABC_COEF1 = -2.0 * (temp0 - 1.0 / temp0) / temp1;
    ABC_COEF2 = 4.0 * (temp0 + 1.0 / temp0) / temp1;
}

inline void sim::abc() {
    int m, n;
    /* ABC at left side of grid */
    #pragma omp teams distribute parallel for
    for (n = 0; n < SIZE_Y; n++) {
        I_ez(0, n) = ABC_COEF0 * (I_ez(2, n) + I_ezLeft(0, 1, n))
            + ABC_COEF1 * (I_ezLeft(0, 0, n) + I_ezLeft(2, 0, n)
            - I_ez(1, n) - I_ezLeft(1, 1, n))
            + ABC_COEF2 * I_ezLeft(1, 0, n) - I_ezLeft(2, 1, n);
        /* memorize old fields */
        for (m = 0; m < 3; m++) {
            I_ezLeft(m, 1, n) = I_ezLeft(m, 0, n);
            I_ezLeft(m, 0, n) = I_ez(m, n);
        }
    }
    /* ABC at right side of grid */
    #pragma omp teams distribute parallel for
    for (n = 0; n < SIZE_Y; n++) {
        I_ez(SIZE_X - 1, n) = ABC_COEF0 * (I_ez(SIZE_X - 3, n) + I_ezRight(0, 1, n))
        + ABC_COEF1 * (I_ezRight(0, 0, n) + I_ezRight(2, 0, n)
        - I_ez(SIZE_X - 2, n) - I_ezRight(1, 1, n))
        + ABC_COEF2 * I_ezRight(1, 0, n) - I_ezRight(2, 1, n);
        /* memorize old fields */
        for (m = 0; m < 3; m++) {
            I_ezRight(m, 1, n) = I_ezRight(m, 0, n);
            I_ezRight(m, 0, n) = I_ez(SIZE_X - 1 - m, n);
        }
    }
    /* ABC at bottom of grid */
    #pragma omp teams distribute parallel for
    for (m = 0; m < SIZE_X; m++) {
        I_ez(m, 0) = ABC_COEF0 * (I_ez(m, 2) + I_ezBot(0, 1, m))
        + ABC_COEF1 * (I_ezBot(0, 0, m) + I_ezBot(2, 0, m)
        - I_ez(m, 1) - I_ezBot(1, 1, m))
        + ABC_COEF2 * I_ezBot(1, 0, m) - I_ezBot(2, 1, m);
        /* memorize old fields */
        for (n = 0; n < 3; n++) {
            I_ezBot(n, 1, m) = I_ezBot(n, 0, m);
            I_ezBot(n, 0, m) = I_ez(m, n);
        }
    }
    /* ABC at top of grid */
    #pragma omp teams distribute parallel for
    for (m = 0; m < SIZE_X; m++) {
        I_ez(m, SIZE_Y - 1) = ABC_COEF0 * (I_ez(m, SIZE_Y - 3) + I_ezTop(0, 1, m))
        + ABC_COEF1 * (I_ezTop(0, 0, m) + I_ezTop(2, 0, m)
        - I_ez(m, SIZE_Y - 2) - I_ezTop(1, 1, m))
        + ABC_COEF2 * I_ezTop(1, 0, m) - I_ezTop(2, 1, m);
        /* memorize old fields */
        for (n = 0; n < 3; n++) {
            I_ezTop(n, 1, m) = I_ezTop(n, 0, m);
            I_ezTop(n, 0, m) = I_ez(m, SIZE_Y - 1 - n);
        }
    }
}

void sim::pecInit(int *data) {
    PEC_HEIGHTS = data;
}

inline void sim::pec() {
    #pragma omp teams distribute parallel for
    for (int i = 0; i < SIZE_X; i++) {
        I_ez(i, PEC_HEIGHTS[i]) = 0.0;
    }
}

inline void sim::updateE() {
    #pragma omp teams distribute parallel for collapse(2)
    for (int m = 1; m < SIZE_X - 1; m++) {
        for (int n = 1; n < SIZE_Y - 1; n++) {
            // Ez(mm, nn) = Ceze(mm, nn) * Ez(mm, nn) +
            // Cezh(mm, nn) * ((Hy(mm, nn) - Hy(mm - 1, nn)) -
            // (Hx(mm, nn) - Hx(mm, nn - 1)));
            I_ez(m, n) = CONST_SAME_FIELD * I_ez(m, n) + 
                CONST_E_DUE_TO_H * ((I_hy(m, n) - I_hy(m - 1, n)) -
                (I_hx(m, n) - I_hx(m, n - 1)));
        }
    }
}

inline void sim::updateH() { 
    int m, n;

    #pragma omp teams distribute parallel for collapse(2)
    for (m = 0; m < SIZE_X; m++) {
        for (n = 0; n < SIZE_Y - 1; n++) {
            //Hx(mm, nn) = Chxh(mm, nn) * Hx(mm, nn) - Chxe(mm, nn) * (Ez(mm, nn + 1) - Ez(mm, nn));
            I_hx(m, n) = CONST_SAME_FIELD*I_hx(m, n) - CONST_H_DUE_TO_E*(I_ez(m, n + 1) - I_ez(m, n));
        }
    }

    #pragma omp teams distribute parallel for collapse(2)
    for (m = 0; m < SIZE_X - 1; m++) {
        for (n = 0; n < SIZE_Y; n++) {
            //Hy(mm, nn) = Chyh(mm, nn) * Hy(mm, nn) + Chye(mm, nn) * (Ez(mm + 1, nn) - Ez(mm, nn));
            I_hy(m, n) = CONST_SAME_FIELD*I_hy(m, n) + CONST_H_DUE_TO_E*(I_ez(m + 1, n) - I_ez(m, n));
        }
    }
}

const double* sim::get_ez() {
    return ez;
}

void sim::run(unsigned int timesteps) {
    #pragma omp target map(to: hx[0:SIZE_X * (SIZE_Y - 1)], \
                                          hy[0:(SIZE_X - 1) * SIZE_Y], \
                                          ez[0:SIZE_X * SIZE_Y], \
                                          ezTop[0:SIZE_X*6], \
                                          ezRight[0:SIZE_Y*6], \
                                          ezBot[0:SIZE_X*6], \
                                          ezLeft[0:SIZE_Y*6])
    {
        N_THREADS = omp_get_num_threads();
        DEVICE_N = omp_get_device_num();
        printf("Devices: %d | Default Device: %d\nDevice #: %d\n# of Threads: %d\nDynamic Thread Count: %d\n", omp_get_num_devices(), omp_get_default_device(), DEVICE_N, N_THREADS, omp_get_dynamic());
        #pragma omp teams num_teams(12)
        {
            for (int step = 0; step < timesteps; step++) {
                // updateH();
                #pragma omp distribute parallel for collapse(2)
                for (int m = 0; m < SIZE_X; m++) {
                    for (int n = 0; n < SIZE_Y - 1; n++) {
                        //Hx(mm, nn) = Chxh(mm, nn) * Hx(mm, nn) - Chxe(mm, nn) * (Ez(mm, nn + 1) - Ez(mm, nn));
                        I_hx(m, n) = CONST_SAME_FIELD*I_hx(m, n) - CONST_H_DUE_TO_E*(I_ez(m, n + 1) - I_ez(m, n));
                    }
                }
                #pragma omp distribute parallel for collapse(2)
                for (int m = 0; m < SIZE_X - 1; m++) {
                    for (int n = 0; n < SIZE_Y; n++) {
                        //Hy(mm, nn) = Chyh(mm, nn) * Hy(mm, nn) + Chye(mm, nn) * (Ez(mm + 1, nn) - Ez(mm, nn));
                        I_hy(m, n) = CONST_SAME_FIELD*I_hy(m, n) + CONST_H_DUE_TO_E*(I_ez(m + 1, n) - I_ez(m, n));
                    }
                }   
                #pragma barrier

                // updateE();
                #pragma omp distribute parallel for collapse(2)
                for (int m = 1; m < SIZE_X - 1; m++) {
                    for (int n = 1; n < SIZE_Y - 1; n++) {
                        // Ez(mm, nn) = Ceze(mm, nn) * Ez(mm, nn) +
                        // Cezh(mm, nn) * ((Hy(mm, nn) - Hy(mm - 1, nn)) -
                        // (Hx(mm, nn) - Hx(mm, nn - 1)));
                        I_ez(m, n) = CONST_SAME_FIELD * I_ez(m, n) + 
                            CONST_E_DUE_TO_H * ((I_hy(m, n) - I_hy(m - 1, n)) -
                            (I_hx(m, n) - I_hx(m, n - 1)));
                    }
                }
                #pragma barrier

                //pec();
                #pragma omp distribute parallel for
                for (int i = 1; i < SIZE_X - 1; i++) {
                    I_ez(i, PEC_HEIGHTS[i]) = 0.0;
                }                

                //abc();
                /* ABC at left side of grid */
                #pragma omp distribute parallel for
                for (int n = 0; n < SIZE_Y; n++) {
                    I_ez(0, n) = ABC_COEF0 * (I_ez(2, n) + I_ezLeft(0, 1, n))
                        + ABC_COEF1 * (I_ezLeft(0, 0, n) + I_ezLeft(2, 0, n)
                        - I_ez(1, n) - I_ezLeft(1, 1, n))
                        + ABC_COEF2 * I_ezLeft(1, 0, n) - I_ezLeft(2, 1, n);
                    /* memorize old fields */
                    for (int m = 0; m < 3; m++) {
                        I_ezLeft(m, 1, n) = I_ezLeft(m, 0, n);
                        I_ezLeft(m, 0, n) = I_ez(m, n);
                    }
                }
                /* ABC at right side of grid */
                #pragma omp distribute parallel for
                for (int n = 0; n < SIZE_Y; n++) {
                    I_ez(SIZE_X - 1, n) = ABC_COEF0 * (I_ez(SIZE_X - 3, n) + I_ezRight(0, 1, n))
                    + ABC_COEF1 * (I_ezRight(0, 0, n) + I_ezRight(2, 0, n)
                    - I_ez(SIZE_X - 2, n) - I_ezRight(1, 1, n))
                    + ABC_COEF2 * I_ezRight(1, 0, n) - I_ezRight(2, 1, n);
                    /* memorize old fields */
                    for (int m = 0; m < 3; m++) {
                        I_ezRight(m, 1, n) = I_ezRight(m, 0, n);
                        I_ezRight(m, 0, n) = I_ez(SIZE_X - 1 - m, n);
                    }
                }
                /* ABC at bottom of grid */
                #pragma omp distribute parallel for
                for (int m = 0; m < SIZE_X; m++) {
                    I_ez(m, 0) = ABC_COEF0 * (I_ez(m, 2) + I_ezBot(0, 1, m))
                    + ABC_COEF1 * (I_ezBot(0, 0, m) + I_ezBot(2, 0, m)
                    - I_ez(m, 1) - I_ezBot(1, 1, m))
                    + ABC_COEF2 * I_ezBot(1, 0, m) - I_ezBot(2, 1, m);
                    /* memorize old fields */
                    for (int n = 0; n < 3; n++) {
                        I_ezBot(n, 1, m) = I_ezBot(n, 0, m);
                        I_ezBot(n, 0, m) = I_ez(m, n);
                    }
                }
                /* ABC at top of grid */
                #pragma omp distribute parallel for
                for (int m = 0; m < SIZE_X; m++) {
                I_ez(m, SIZE_Y - 1) = ABC_COEF0 * (I_ez(m, SIZE_Y - 3) + I_ezTop(0, 1, m))
                + ABC_COEF1 * (I_ezTop(0, 0, m) + I_ezTop(2, 0, m)
                - I_ez(m, SIZE_Y - 2) - I_ezTop(1, 1, m))
                + ABC_COEF2 * I_ezTop(1, 0, m) - I_ezTop(2, 1, m);
                /* memorize old fields */
                for (int n = 0; n < 3; n++) {
                        I_ezTop(n, 1, m) = I_ezTop(n, 0, m);
                        I_ezTop(n, 0, m) = I_ez(m, SIZE_Y - 1 - n);
                    }
                }
                #pragma barrier

                // Hardwired additive sourcenode in center
                I_ez(100, 150) += source->source(step, 0.0);

                //s->run(ez);
            }
        }
    }
}

void sim::reset() {
    // Reset hx
    for (int i = 0; i < SIZE_X * (SIZE_Y - 1); i++) {
        hx[i] = 0;
    }
    // Reset hy
    for (int i = 0; i < (SIZE_X - 1) * SIZE_Y; i++) {
        hy[i] = 0;
    }
    // Reset ez
    for (int i = 0; i < SIZE_X * SIZE_Y; i++) {
        ez[i] = 0;
    }

    // Reset horizontal ABCs' buffers
    for (int i = 0; i < SIZE_X*6; i++) {
        ezTop[i] == 0;
        ezBot[i] == 0;
    }
    // Reset vertical ABCs' buffers
    for (int i = 0; i < SIZE_Y*6; i++) {
        ezLeft[i] == 0;
        ezRight[i] == 0;
    }
}

int sim::getNumThreads() {
    return N_THREADS;
}

void sim::setDesiredThreads(int num = 0) {
    if (num == 0) {
        omp_set_dynamic(0);
        //omp_set_dynamic(1); // Allow thread count to change
        int num = omp_get_num_procs();
    } else {
        omp_set_dynamic(0); // Ensure thread count is not changing throughout
    }
    printf("Requesting %d threads\n", num);
    omp_set_num_threads(num);
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
