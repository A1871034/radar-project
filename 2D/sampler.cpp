#include "sampler.h"

#include "sim.h"

sampler::sampler(unsigned int sampleEvery, unsigned int sizeX, unsigned int sizeY, std::string fileName)
    : sampleEvery(sampleEvery), SIZE_X(sizeX), SIZE_Y(sizeY), fileName(fileName)
{
    // Gate no writing
    if (sampleEvery == 0) {
        return;
    }

    prepareWriting();
}

sampler::~sampler() {
    if (sampleEvery != 0) {
        fclose(fileHandle);
    }
}

void sampler::setSampleEvery(unsigned int _sampleEvery)
{
    if (sampleEvery == 0 && _sampleEvery != 0) {
        prepareWriting();
    } else if (sampleEvery != 0 && _sampleEvery == 0) {
        fclose(fileHandle);
    }
    sampleEvery = _sampleEvery;
}

void sampler::prepareWriting() {
    #if DEBUG_SAMPLE
    // Clear previous data and write header
    if (!fileHeaderExists) {
        fileHandle = fopen(fileName.c_str(), "w");
        fprintf(fileHandle, "%u,%u,%lu\n", SIZE_X, SIZE_Y, sizeof(double));
        fileHeaderExists = true;
        fclose(fileHandle);
    }
    #endif
    // Prepare for byte rights by sampler::run
    fileHandle = fopen(fileName.c_str(), "ab");
}

void sampler::run(const double* ez)
{
    if ((sampleEvery != 0) && (step++ % sampleEvery == 0)) {
        printf("step: %d\n", step - 1);
        #if DEBUG_SAMPLE
        fwrite(ez, sizeof(double), SIZE_X * SIZE_Y, fileHandle);
        #else
        fwrite(&I_ez(100, 150), sizeof(double), 1, fileHandle);
        #endif
    }

}