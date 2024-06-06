#include "height_data.h"

height_data::height_data(const char *filepath) {
    file.open(filepath, std::ios::in|std::ios::binary|std::ios::ate);
    if (!file.is_open()) {
        std::cout << "File Not Open" << std::endl;
        exit(1);
    }
    size = file.tellg();
    raw_data = new char [size];
    file.seekg(0, std::ios::beg);
    file.read(raw_data, size);
    file.close();
    parse();
}

height_data::~height_data() {
    delete data;
}

void height_data::parse() {
    data = (int *) raw_data;
    size = size/4;
}

int * height_data::get_data() {
    return data;
}

int height_data::get_x() {
    return int(size);
}

int height_data::get_min_y() {
    int ma = 1;
    for (int i = 0; i < size; i++) {
        ma = data[i] > ma ? data[i] : ma;
    }
    return ma;
}

void height_data::set_max_y(int max_height) {
    for (int i = 0; i < size; i++) {
        data[i] = (data[i] > max_height) ? max_height : data[i];
    }
}