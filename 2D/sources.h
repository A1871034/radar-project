#pragma once

#include "constants.h"

class sources {
    double ppw = 0;
    public:
    sources();
    sources(double ppw) : ppw(ppw) {}
    double source(unsigned int step, double location);
};