#include "gas_jet.h"
#include "maths/maths.h"


CylindricalGasJet::CylindricalGasJet() : gen(rd()), dis_xz(0., 1.), dis_y(-_length / 2., _length / 2.), dis_t(0., 2. * Pi), _density(0), _sigma(0) {}
CylindricalGasJet::CylindricalGasJet(double density_cm3, double sigma_um, double radius_um, double length_um, size_t n) :
    _n(n), _length(length_um * umToAU), _radius(radius_um * umToAU), _sigma(sigma_um * umToAU), _density(density_cm3 / cmToAU / cmToAU / cmToAU), 
    gen(rd()), dis_xz(0., _sigma), dis_y(-_length / 2., _length / 2.), dis_t(0., 2. * Pi) {
}



void CylindricalGasJet::SampleCylinder(const std::vector<Laser>& lasers) {
    if (_n == 1) {
        Cell cell;
        cell.pos.x = 0;
        cell.pos.y = 0;
        cell.pos.z = 0;
        cell.density = 1;
        double rho = sqrt(cell.pos.x * cell.pos.x + cell.pos.y * cell.pos.y);

        cell.intensity.resize(lasers.size());
        cell.phase.resize(lasers.size());
        for (int i = 0; i < lasers.size(); i++) {
            cell.intensity[i] = lasers[i].IntensityAt(rho, cell.pos.z);
            cell.phase[i] = lasers[i].Phase(rho, cell.pos.z);
        }
        cells.push_back(cell);
    }
    else {
        for (int i = 0; i < _n;) {
            Cell cell;
            cell.pos.z = dis_xz(gen);
            cell.pos.x = dis_xz(gen);
            cell.pos.y = dis_y(gen);
            cell.density = 1; //_density*exp(-0.5*R*R/_sigma/_sigma);
            double rho = sqrt(cell.pos.x * cell.pos.x + cell.pos.y * cell.pos.y);

            cell.intensity.resize(lasers.size());
            cell.phase.resize(lasers.size());
            for (int i = 0; i < lasers.size(); i++) {
                cell.intensity[i] = lasers[i].IntensityAt(rho, cell.pos.z);
                cell.phase[i] = lasers[i].Phase(rho, cell.pos.z);
            }

            for (int j = 0; j < lasers.size(); j++) {
                if (cell.intensity[j] / lasers[j].PeakI0() > 0.032) {
                    cells.push_back(cell);
                    i++;
                    break;
                }
            }
        }
    }
}

CylindricalGasJet& CylindricalGasJet::operator=(const CylindricalGasJet& o) {
    _density = o._density;
    _sigma = o._sigma;
    _radius = o._radius;
    _length = o._length;
    _n = o._n;
    dis_xz = o.dis_xz;
    dis_y = o.dis_y;
    dis_t = o.dis_t;
    cells = o.cells;

    return *this;
}