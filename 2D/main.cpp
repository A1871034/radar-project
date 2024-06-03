#include "sim.h"
#include "sampler.h"
#include "height_data.h"
#include "tester.h"

#include <omp.h>

int main() {
    unsigned int PPW = 20;
    int gpu_device = 0;

    // Setup Height
    height_data hd ("./heights.data");
    int width = hd.get_x();
    int height = hd.get_min_y() + PPW*10;
    height = 1000;
    int steps = 100;//width;

    // Setup Sampler
    sampler sa = sampler(0, width, height, "sim.data");

    // OpenMP Device Info
    printf("Default device: %d\n", omp_get_default_device());
    printf("Devices: %d\n", omp_get_num_devices());

    // Setup tester
    Tester tester(nullptr, steps, true);
    tester.initFile("time_test.csv");

    #if DEBUG_SAMPLE
        sim si(width, height, PPW);
        si.pecInit(hd.get_data());
        si.initSampler();
        si.setSamplePeriod(20);
        tester.setSim(&si);
        // Test Sim
        tester.test(gpu_device);
        return 0;
    #endif

    // Initialise Testing Parameters
    int heights[] = {1000, 2500, 5000, 7500, 10000, 12500, 15000};

    // Run test    
    for (int height: heights) {
        // Setup simulation
        sim si(width, height, PPW);
        si.pecInit(hd.get_data());
        tester.setSim(&si);
        // Test Sim
        tester.test(gpu_device);
    }

    return 0;
}