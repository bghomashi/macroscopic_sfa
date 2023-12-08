#include "save_data.h"
#include <iostream>
#include <iomanip>
#include <fstream>

bool Store(const std::string& filename, const dvector& domain, const c2vector& data) {
    int nf = int(domain.size());
    std::ofstream file(filename);
    if (!file.is_open())
        return false;
    file << std::setprecision(8) << std::scientific;

    for (int iff = 0; iff < nf; iff++) {
        file << domain[iff] << "\t"
             << std::real(data[iff].x) << "\t"
             << std::imag(data[iff].x) << "\t"
             << std::real(data[iff].y) << "\t"
             << std::imag(data[iff].y) << "\n";

    }
    return true;
}
bool Store(const std::string& filename, const dvector& domain, const d2vector& data) {
    int nf = int(domain.size());
    std::ofstream file(filename);
    if (!file.is_open())
        return false;
    file << std::setprecision(8) << std::scientific;

    for (int iff = 0; iff < nf; iff++) {
        file << domain[iff] << "\t"
             << data[iff].x << "\t"
             << data[iff].y << "\n";
    }
    return true;
}
bool Store(const std::string& filename, const dvector& domain, const cvector& data) {
    int nf = int(domain.size());
    std::ofstream file(filename);
    if (!file.is_open())
        return false;
    file << std::setprecision(8) << std::scientific;

    for (int iff = 0; iff < nf; iff++) {
        file << domain[iff] << "\t"
             << std::real(data[iff]) << "\t"
             << std::imag(data[iff]) << "\n";
    }
    return true;
}
bool Store(const std::string& filename, const dvector& domain, const dvector& data) {
    int nf = int(domain.size());
    std::ofstream file(filename);
    if (!file.is_open())
        return false;
    file << std::setprecision(8) << std::scientific;

    for (int iff = 0; iff < nf; iff++) {
        file << domain[iff] << "\t"
             << data[iff] << "\n";
    }
    return true;
}