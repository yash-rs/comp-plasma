#include "fd/poisson1d/scheme.hpp"
#include "fd/grid1d.hpp"


namespace fd::poisson1d{
    SecondOrderScheme::SecondOrderScheme() {}

    // center coefficient
    double SecondOrderScheme::centerCoeff(double h) const{
        return 2.0 / (h*h);
    }

    double SecondOrderScheme::neighborCoeff(double h) const{
        return -1.0 / (h*h);
    }

    void SecondOrderScheme::addInteriorRows(std::vector<Eigen::Triplet<double>>& triplets, const Grid1D& grid, int firstRow, int lastRow) const{
        
        const double h = grid.h();
        for (int row=firstRow; row<=lastRow; row++){
            triplets.emplace_back(row, row, centerCoeff(h));
            triplets.emplace_back(row, row-1, neighborCoeff(h));
            triplets.emplace_back(row, row+1, neighborCoeff(h));
        }
    }
}
