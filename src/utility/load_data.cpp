#include "load_data.h"
#include <iostream>
#include <iomanip>
#include <fstream>


bool LoadBinary(const std::string& filename, c2vector& data, int index) {
    std::ifstream file(filename, std::ios::in | std::ios::binary);
    if (!file.is_open())
        return false;
            
    size_t size;
    file.seekg(0, std::ios::end);
    size = file.tellg();
    file.seekg(0, std::ios::beg);

    if (index == 0) {
        data.resize(size / sizeof(complex));
        int nf = int(data.size());
        for (int iff = 0; iff < nf; iff++) {
            file.read((char*)&data[iff].x, sizeof(complex));
            data[iff].y = 0;
        }
    } else if (index == 1) {
        data.resize(size / sizeof(complex));
        int nf = int(data.size());
        for (int iff = 0; iff < nf; iff++) {
            data[iff].x = 0;
            file.read((char*)&data[iff].y, sizeof(complex));
        }
    } else {
        data.resize(size / sizeof(cvec2));
        file.read((char*)data.data(), data.size() * sizeof(cvec2));
    }
    return true;
}

bool LoadBinary(const std::string& filename, dvector& data) {
    int nf = int(data.size());
    std::ifstream file(filename, std::ios::in | std::ios::binary);
    if (!file.is_open())
        return false;
    
    size_t size;
    file.seekg(0, std::ios::end);
    size = file.tellg();
    file.seekg(0, std::ios::beg);

    data.resize(size / sizeof(double));
    file.read((char*)data.data(), data.size() * sizeof(double));
    
    return true;
}