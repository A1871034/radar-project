#pragma once

#include "sim.h"

/* 
Can output to CSV like
Device, Threads, Steps, Width, Height, Time, Time Per Step, Time Per Cell 
*/
class Tester {
    sim* s;
    int STEPS;
    bool VERBOSE = false;
    bool OUTPUT_TO_FILE = false;
    FILE * fileHandle;
    
    uint64_t timeRun(); // Nanoseconds
    uint64_t timePerStep(uint64_t runTime);
    double timePerCell(uint64_t runTime);
    public:
    Tester(sim* s, int steps, bool verbose = false) : s(s), STEPS(steps), VERBOSE(verbose) {}
    ~Tester() {}
    void initFile(std::string fileName);
    void test(int numThreads, int deviceNum);
};