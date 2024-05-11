#pragma once

#include <string>
#include <fstream>
#include <iostream>
#include <stdio.h>

class sampler {
    unsigned int sampleEvery = 2;
    unsigned int step = 0;
    unsigned int SIZE_X;
    unsigned int SIZE_Y;
    std::string fileName;
    FILE * fileHandle;
    bool fileHeaderExists = false;

    public:
        sampler(unsigned int sampleEvery, unsigned int sizeY, unsigned int sizeX, std::string fileName);
        sampler() {}
        ~sampler();
        void setSampleEvery(unsigned int _sampleEvery);
        void prepareWriting();
        void run(const double *ez);
};