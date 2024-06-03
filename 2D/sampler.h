#pragma once

#include <string>
#include <fstream>
#include <iostream>
#include <stdio.h>

#define DEBUG_SAMPLE false

class sampler {
    unsigned int sampleEvery = 2;
    unsigned int step = 0;
    unsigned int SIZE_X;
    unsigned int SIZE_Y;
    unsigned int SAMPLE_X;
    unsigned int SAMPLE_Y;
    float * sample_buffer = nullptr;
    std::string fileName;
    FILE * fileHandle;
    bool fileHeaderExists = false;

    public:
        sampler(unsigned int sampleEvery, unsigned int sizeY, unsigned int sizeX, std::string fileName=(DEBUG_SAMPLE ? "sim.debug" : "sim.data"), \
                unsigned int sample_x = 0, unsigned int sample_y = 0);
        sampler() {}
        ~sampler();
        void setSampleEvery(unsigned int _sampleEvery);
        void prepareWriting();
        void run(const double *ez);
        void setSampleSize(unsigned int sample_x, unsigned int sample_y);
};