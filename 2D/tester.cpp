#include "tester.h"

#include <chrono>
#include <stdio.h>

using namespace std::chrono;

void Tester::initFile(std::string fileName) {
    fileHandle = fopen(fileName.c_str(), "w");
    fprintf(fileHandle, "Steps, Width, Height, Time (ns), Time Per Step (ns), Time Per Cell (ns)\n");
    OUTPUT_TO_FILE = true;
}

/* Nanoseconds */
uint64_t Tester::timeRun() {
    auto start = steady_clock::now();
    s->run(STEPS);
    return (steady_clock::now() - start).count();
}   

uint64_t Tester::timePerStep(uint64_t runTime) {
    return runTime/STEPS;
}

float Tester::timePerCell(uint64_t runTime) {
    return float(runTime)/float(STEPS*(s->getSizeX())*(s->getSizeY()));
}

void Tester::test() {
    s->reset();
    uint64_t duration = timeRun();

    // Device cannot be GPU here
    unsigned int width = s->getSizeX();
    unsigned int height = s->getSizeY();
    uint64_t tps = timePerStep(duration);
    float tpc = timePerCell(duration);

    if (VERBOSE) {
        printf("ns: %lu\n", duration);
        printf("Steps: %d | Dim: %ux%u\nTime (s): %.2f\nTime Per Step (ms): %.2f\nTime Per Cell (ns): %.4f\n--------------\n",
            STEPS, width, height, float(duration)/float(pow(10, 9)), float(tps)/float(pow(10,6)), tpc);
    }

    if (OUTPUT_TO_FILE) {
        fprintf(fileHandle, "%d, %u, %u, %lu, %lu, %f\n", STEPS, width, height, duration, tps, tpc);
    }
}

void Tester::setSim(sim *_s) {
    s = _s;
}
