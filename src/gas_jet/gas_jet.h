#pragma once

#include <random>
#include "maths/maths.h"
#include "maths/vec3.h"
#include "laser/laser.h"

class CylindricalGasJet {
    double _density, _sigma;
    double _radius, _length;
    size_t _n;

    std::random_device rd;
    std::mt19937 gen;
    // std::uniform_real_distribution<double> dis_xz;
    std::normal_distribution<double> dis_xz;
    std::uniform_real_distribution<double> dis_y;
    std::uniform_real_distribution<double> dis_t;
public:
    struct Cell {
        dvec3 pos;
        double density;

        std::vector<double> intensity;
        std::vector<double> phase;
    };

    std::vector<Cell> cells;

    CylindricalGasJet();
    CylindricalGasJet(double density_cm3, double sigma_um, double radius_um, double length_um, size_t num);

    void SampleCylinder(const std::vector<Laser>& lasers);
    CylindricalGasJet& operator=(const CylindricalGasJet& o);
};