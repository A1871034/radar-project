#pragma once

#include "sampler.h"
#include "constants.h"

#include <math.h>
#include <string>
#include <fstream>
#include <iostream>
#include <stdio.h>

#define I_hx(m, n) hx[(m) * (SIZE_Y - 1) + n]
#define I_hy(m, n) hy[(m) * SIZE_Y + n]
#define I_ez(m, n) ez[(m) * SIZE_Y + n]

#define I_ezLeft(M, Q, N) ezLeft[(N) * 6 + (Q) * 3 + (M)]
#define I_ezRight(M, Q, N) ezRight[(N) * 6 + (Q) * 3 + (M)]
#define I_ezTop(N, Q, M) ezTop[(M) * 6 + (Q) * 3 + (N)]
#define I_ezBot(N, Q, M) ezBot[(M) * 6 + (Q) * 3 + (N)]

#define DEBUG_SAMPLE true

class sim {
    friend class sampler;

    double * hx;
    double * hy;
    double * ez;

    unsigned int SIZE_X;
    unsigned int SIZE_Y;
    unsigned int PPW;
    #if DEBUG_SAMPLE
    FILE * FILE_HANDLE;
    unsigned int SAMPLE_PERIOD = 199;
    #endif

    // Update E_CONST, H_CONST with a grid and stored values to have different materials at locations
    double CONST_SAME_FIELD = 1.0;
    double CONST_E_DUE_TO_H = CDTDS*IMP0;
    double CONST_H_DUE_TO_E = CDTDS/IMP0;

    // 2nd Order ABC Variables
    double ABC_COEF0;
    double ABC_COEF1;
    double ABC_COEF2;
    double * ezTop;
    double * ezRight;
    double * ezBot;
    double * ezLeft;

    // Terrain reflection
    int * PEC_HEIGHTS = nullptr;

    int N_THREADS = -1;
    int DEVICE_N = -1;

    sampler * s;

    public:
    sim(unsigned int sizeX, unsigned int sizeY, unsigned int PPW);
    ~sim();
    void setSampler(sampler * _s);
    void abcInit();
    void pecInit(int * data);
    const double *get_ez();
    void run(unsigned int timesteps);
    int getNumThreads();
    int getDeviceNum();
    unsigned int getSizeX();
    unsigned int getSizeY();
    #if DEBUG_SAMPLE
    void initSampler();
    void setSamplePeriod(unsigned int period);
    #endif
};