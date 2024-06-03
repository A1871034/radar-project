#include "sampler.h"

#include "sim.h"

sampler::sampler(unsigned int sampleEvery, unsigned int sizeX, unsigned int sizeY, std::string fileName, unsigned int sample_x, unsigned int sample_y)
    : sampleEvery(sampleEvery), SIZE_X(sizeX), SIZE_Y(sizeY), fileName(fileName)
{
    // TODO: Fix backwards x and y throughout sim
    if (sample_x == 0 || sample_y == 0) {
        SAMPLE_X = SIZE_X;
        SAMPLE_Y = SIZE_Y;
    } else {
        SAMPLE_X = sample_x;
        SAMPLE_Y = sample_y;
    }
    sample_buffer = (float *) malloc(SAMPLE_X*SAMPLE_Y*sizeof(float));
    
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
        fprintf(fileHandle, "%u,%u,%lu\n", SAMPLE_X, SAMPLE_Y, sizeof(float));
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
        #pragma omp parallel for collapse(2)
        for (int i = 0; i < SAMPLE_X; i++) {
            for (int j = 0; j < SAMPLE_Y; j++) {
                sample_buffer[i * SAMPLE_Y + j] = I_ez(i, j);
            }
        }
        #if DEBUG_SAMPLE
        fwrite(sample_buffer, sizeof(float), SAMPLE_Y * SAMPLE_X, fileHandle);
        #else
        fwrite(&I_ez(100, 150), sizeof(double), 1, fileHandle);
        #endif
    }

}


// Don't change sample size on the fly
void sampler::setSampleSize(unsigned int sample_x, unsigned int sample_y) {
    SAMPLE_X = sample_x;
    SAMPLE_Y = sample_y;


}
