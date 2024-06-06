#pragma once

#include <iostream>
#include <fstream>

class height_data {
    std::ifstream file;
    std::streampos size;
    char * raw_data;
    int * data;
    public:
    height_data(const char * filepath);
    ~height_data();
    void parse();
    int * get_data();
    int get_x();
    int get_min_y();
    void set_max_y(int max_height);
};