#include "sim.h"

#include <cstdlib>
#include <omp.h>
#include <math.h>

// Update Equations and ABCs from https://eecs.wsu.edu/~schneidj/ufdtd/ufdtd.pdf


sim::sim(unsigned int sizeX, unsigned int sizeY, unsigned int PPW) : SIZE_X(sizeX), SIZE_Y(sizeY), PPW(PPW) {
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

void sim::pecInit(int *data) {
    PEC_HEIGHTS = data;
}

const double* sim::get_ez() {
    return ez;
}

void sim::run(unsigned int timesteps) {
    for (unsigned int step = 0; step < timesteps; step++) {
        // updateH();
        #pragma omp parallel for collapse(2)
        for (unsigned int m = 0; m < SIZE_X; m++) {
            for (unsigned int n = 0; n < SIZE_Y - 1; n++) {        
                //Hx(mm, nn) = Chxh(mm, nn) * Hx(mm, nn) - Chxe(mm, nn) * (Ez(mm, nn + 1) - Ez(mm, nn));
                I_hx(m, n) = CONST_SAME_FIELD*I_hx(m, n) - CONST_H_DUE_TO_E*(I_ez(m, n + 1) - I_ez(m, n));
            }
        }

        #pragma omp parallel for collapse(2)
        for (unsigned int m = 0; m < SIZE_X - 1; m++) {
            for (unsigned int n = 0; n < SIZE_Y; n++) {
                //Hy(mm, nn) = Chyh(mm, nn) * Hy(mm, nn) + Chye(mm, nn) * (Ez(mm + 1, nn) - Ez(mm, nn));
                I_hy(m, n) = CONST_SAME_FIELD*I_hy(m, n) + CONST_H_DUE_TO_E*(I_ez(m + 1, n) - I_ez(m, n));
            }
        }   

        // updateE();
        #pragma omp parallel for collapse(2)
        for (unsigned int m = 1; m < SIZE_X - 1; m++) {
            for (unsigned int n = 1; n < SIZE_Y - 1; n++) {
                // Ez(mm, nn) = Ceze(mm, nn) * Ez(mm, nn) +
                // Cezh(mm, nn) * ((Hy(mm, nn) - Hy(mm - 1, nn)) -
                // (Hx(mm, nn) - Hx(mm, nn - 1)));
                I_ez(m, n) = CONST_SAME_FIELD * I_ez(m, n) + 
                    CONST_E_DUE_TO_H * ((I_hy(m, n) - I_hy(m - 1, n)) -
                    (I_hx(m, n) - I_hx(m, n - 1)));
            }
        }

        // pec();
        #pragma omp parallel for
        for (unsigned int i = 1; i < SIZE_X - 1; i++) {
            I_ez(i, PEC_HEIGHTS[i]) = 0.0;
        }                

        // abc();
        /* ABC at left/right side of grid */
        #pragma omp parallel for
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
        /* ABC at bottom/top of grid */
        #pragma omp parallel for
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

        // Hardwired additive sourcenode
        double arg = M_PI * ((CDTDS * step - 0.0) / PPW - 1.0);
        arg = arg * arg;
        I_ez(100, 150) += (1.0 - 2.0 * arg) * exp(-arg); 

        // Sample
        //s->run(ez);
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

void sim::setDesiredThreads(int num = 0) {
    omp_set_dynamic(0); // Ensure thread count is not changing throughout
    printf("Requesting %d threads\n", num);
    omp_set_num_threads(num);
    THREADS = num;
    // Verifying here would make sense but don't want to work too much on CPU implementation
}

unsigned int sim::getSizeX() {
    return SIZE_X;
}

unsigned int sim::getSizeY() {
    return SIZE_Y;
}
