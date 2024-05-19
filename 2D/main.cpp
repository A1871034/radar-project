#include "sim.h"
#include "sampler.h"
#include "sources.h"
#include "height_data.h"
#include "tester.h"

#include <omp.h>


int main() {

    height_data hd ("./heights.data");
    int width = hd.get_x();
    int height = hd.get_min_y() + 200;
    int steps = 1000;

    sampler sa = sampler(0, width, height, "sim.data");
    sources so = sources(20);

    sim si(width, height);
    si.pecInit(hd.get_data());
    si.setSampler(&sa);
    si.setSource(&so);
    
    Tester tester(&si, steps, true);
    tester.initFile("time_test.csv");
    tester.test(12, 0);
    tester.test(8, 0);
    tester.test(6, 0);
    tester.test(4, 0);
    tester.test(2, 0);
    tester.test(1, 0);

    //printf("Default %d\nNum Devices %d\nInitial/Host %d",omp_get_default_device(), omp_get_num_devices(), omp_get_initial_device());
    return 0;
}