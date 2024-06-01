#include "tester.h"

#include <omp.h>
#include <chrono>
#include <stdio.h>

using namespace std::chrono;

void Tester::initFile(std::string fileName) {
    fileHandle = fopen(fileName.c_str(), "w");
    fprintf(fileHandle, "Device, Threads, Steps, Width, Height, Time (ns), Time Per Step (ns), Time Per Cell (ns)\n");
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

double Tester::timePerCell(uint64_t runTime) {
    return double(runTime)/double(STEPS*(s->getSizeX())*(s->getSizeY()));
}

void Tester::test(int numThreads, int deviceNum = -1) {
    if (deviceNum != -1) {
        omp_set_default_device(deviceNum);
    }

    if (numThreads == 0) {
        omp_set_dynamic(1); // Allow thread count to change
        numThreads = omp_get_num_procs();
    }
    s->setDesiredThreads(numThreads);

    // Kind of verify we can get that thread count
    int max_thread_count_got = omp_get_max_threads(); 
    if (max_thread_count_got != numThreads) {
        printf("Undesired thread-count: Requested %d, got %d\nSkipping current test\n", numThreads, max_thread_count_got);
        return;
    }

    s->reset();
    uint64_t duration = timeRun();

    // Device cannot be GPU here
    unsigned int width = s->getSizeX();
    unsigned int height = s->getSizeY();
    uint64_t tps = timePerStep(duration);
    double tpc = timePerCell(duration);

    if (VERBOSE) {
        printf("Device: CPU | Threads: %d\nSteps: %d | Dim: %ux%u\nTime (s): %.2f\nTime Per Step (ms): %.2f\nTime Per Cell (ns): %.4f\n--------------\n",
            numThreads, STEPS, width, height, double(duration)/double(pow(10, 9)), double(tps)/double(pow(10,6)), tpc);
    }

    if (OUTPUT_TO_FILE) {
        fprintf(fileHandle, "CPU, %d, %d, %u, %u, %lu, %lu, %f\n", numThreads, STEPS, width, height, duration, tps, tpc);
    }
}

void Tester::setSim(sim *_s) {
    s = _s;
}
