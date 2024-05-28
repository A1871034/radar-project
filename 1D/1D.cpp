#include <stdio.h>
#include <math.h>
#include <iostream>
#include <fstream>
#include <sstream>

#define SIZE 200

using namespace std;

int main()
{
    double ez[SIZE] = {0.}, hy[SIZE] = {0.}, imp0 = 377.0;
    int qTime, maxTime = 580, mm;

    stringstream ss;

    ofstream myfile;
    myfile.open("test.csv");
    //myfile << "T,Ez[25]\n";

    /* do time stepping */
    for (qTime = 0; qTime < maxTime; qTime++) {

        /* simple ABC for hy[Size-1] */
        hy[SIZE-1] = hy[SIZE-2];

        /* update magnetic field */
        for (mm = 0; mm < SIZE - 1; mm++)
        hy[mm] = hy[mm] + (ez[mm + 1] - ez[mm]) / imp0;

        /* simple ABC for ez[0] */
        ez[0] = ez[1];

        /* update electric field */
        for (mm = 1; mm < SIZE; mm++)
        ez[mm] = ez[mm] + (hy[mm] - hy[mm - 1]) * imp0;

        /* hardwire a source node */
        ez[67] += exp(-(qTime - 30.) * (qTime - 30.) / 100.);
        ez[133] += exp(-(qTime - 30.) * (qTime - 30.) / 100.);

        hy[49] -= exp(-(qTime - 200.) * (qTime - 200.) / 100.) / imp0;
        ez[50] += exp(-(qTime - 200.) * (qTime - 200.) / 100.);
        ez[150] -= exp(-(qTime - 200.) * (qTime - 200.) / 100.);
        hy[150] -= exp(-(qTime - 200.) * (qTime - 200.) / 100.) / imp0;

        ez[50] += exp(-(qTime - 400.) * (qTime - 400.) / 100.);
        ez[100] -= exp(-(qTime - 400.) * (qTime - 400.) / 100.);
        ez[150] += exp(-(qTime - 400.) * (qTime - 400.) / 100.);

        
        //printf("%g\n", ez[50]);
        if (qTime%2 == 0) {
            for (int i = 0;i < SIZE-1; i++) {
                ss << std::to_string(ez[i]) + ",";
            }
            ss << std::to_string(ez[-1]) << "\n";
        }
        
    } /* end of time-stepping */

    myfile << ss.str();

    return 0;
}