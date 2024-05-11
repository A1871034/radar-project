#include "sources.h"

#include <cstdio>
#include <math.h>
#include <stdexcept>

/* initialize source-function variables */
sources::sources() {
    printf("Enter the points per wavelength for Ricker source: ");
    scanf(" %lf", &ppw);
    if (ppw < 0) {
        std::runtime_error("Points per wavelength must be positive.");
    }
}

/* calculate source function at given time and location */
double sources::source(unsigned int step, double location) {
    double arg = M_PI * ((CDTDS * step - location) / ppw - 1.0);
    arg = arg * arg;
    return (1.0 - 2.0 * arg) * exp(-arg);
}
