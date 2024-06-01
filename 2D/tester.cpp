#include "tester.h"

#include <omp.h>
#include <chrono>
#include <stdio.h>

using namespace std::chrono;

void Tester::initFile(std::string fileName) {
    fileHandle = fopen(fileName.c_str(), "w");
    fprintf(fileHandle, "Device, Threads, Teams, Steps, Width, Height, Time (ns), Time Per Step (ns), Time Per Cell (ns)\n");
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
    printf("Default device: %d\n", omp_get_default_device());
    printf("Devices: %d\n", omp_get_num_devices());
    if (deviceNum != -1) {
        // DO something idk yet, will do once actually on gpu
        omp_set_default_device(deviceNum);
    }

    if (numThreads == 0) {
        omp_set_dynamic(1); // Allow thread count to change
        //numThreads = omp_get_num_procs();
    } else {
        s->setDesiredThreads(numThreads);
    }

    uint64_t duration = timeRun();

    if (s->getNumThreads() != numThreads) {
        // Doing this after is a bit silly but whatever for now
        printf("Undesired thread-count: Requested %d, got %d\n", numThreads, s->getNumThreads());
    }

    int device = s->getDeviceNum();
    int threads = s->getNumThreads();
    int teams = s->getNumTeams();
    unsigned int width = s->getSizeX();
    unsigned int height = s->getSizeY();
    uint64_t tps = timePerStep(duration);
    double tpc = timePerCell(duration);

    if (VERBOSE) {
        printf("Device: %d | Threads: %d | Teams: %d\nSteps: %d | Dim: %ux%u\nTime (s): %.2f\nTime Per Step (ms): %.2f\nTime Per Cell (ns): %.4f\n--------------\n",
            device, threads, teams, STEPS, width, height, double(duration)/double(pow(10, 9)), double(tps)/double(pow(10,6)), tpc);
    }

    if (OUTPUT_TO_FILE) {
        fprintf(fileHandle, "%d, %d, %d, %d, %u, %u, %lu, %lu, %f\n", device, threads, teams, STEPS, width, height, duration, tps, tpc);
    }
}
