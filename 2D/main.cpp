#include "sim.h"
#include "sampler.h"
#include "sources.h"
#include "height_data.h"

int main() {

    height_data hd ("./heights.data");
    int width = hd.get_x();
    int height = hd.get_min_y() + 200;

    sampler sa = sampler(1, width, height, "sim.data");
    sources so = sources(20);

    sim si(width, height);
    si.pecInit(hd.get_data());
    si.setSampler(&sa);
    si.setSource(&so);
    
    si.run(20);
    return 0;
}