#include "sim.h"
#include "sampler.h"
#include "height_data.h"
#include "tester.h"

#include <omp.h>

int main() {
    unsigned int PPW = 20;

    // Setup Height
    height_data hd ("./heights.data");
    int width = hd.get_x();
    int height = hd.get_min_y() + PPW*10; // 10 Wavelengths over top of highest point
    height = 500;
    int steps = 1050;
    
    // Setup Sampler
    sampler sa = sampler(8, 1000, height, std::string((DEBUG_SAMPLE ? "sim.debug" : "sim.data")), 1001, 501);

    // OpenMP Device Info
    printf("Default device: %d\n", omp_get_default_device());
    printf("Devices: %d\n", omp_get_num_devices());
    
    // Setup Tester
    Tester tester(nullptr, steps, true);
    tester.initFile("time_test.csv");

    #if DEBUG_SAMPLE
    sim si(width, height, PPW);
    si.pecInit(hd.get_data());
    si.setSampler(&sa);
    tester.setSim(&si);
    // Test Sim
    tester.test(12,1);
    return 0;
    #endif

    // Initialise Testing Parameters
    int heights[] = {1000, 2500, 5000, 7500, 10000, 12500, 15000};
    int threads[] = {1,2,4,8,12};

    for (int height: heights) {
        // Setup simulation
        sim si(width, height, PPW);
        si.pecInit(hd.get_data());
        si.setSampler(&sa);
        tester.setSim(&si);
        for (int thread_count: threads) {
            // Test Sim
            tester.test(thread_count,1);
        }
    }

    return 0;
}