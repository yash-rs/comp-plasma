#include "fd/poisson1d/scheme.hpp"
#include "fd/grid1d.hpp"
#include "fd/stencil.hpp"
#include <stdexcept>
#include <utility>
#include <vector>


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

    Eigen::VectorXd GenScheme::computeWeights(const std::vector<int>& offsets){
        Eigen::VectorXd o(offsets.size());
        for (int i=0;i<static_cast<int>(offsets.size());i++){
            o(i) = offsets[i];
        }
        return stencilWeights(o, 2);
    }

    GenScheme::GenScheme(const std::vector<int>& intOffsets,
                         const std::vector<int>& leftOffsets,
                         const std::vector<int>& rightOffsets)
            : intOffsets_(intOffsets),
              leftOffsets_(leftOffsets),
              rightOffsets_(rightOffsets),
              intWeights_(computeWeights(intOffsets)),
              leftWeights_(computeWeights(leftOffsets)),
              rightWeights_(computeWeights(rightOffsets)) {}

    std::pair<std::vector<int>, Eigen::VectorXd> GenScheme::nodeWeights(int i, const Grid1D& grid) const{
        int N = grid.N();
        double h = grid.h();
        double scale = -1 / (h*h);
        if (i<1 || i>N-1){
            throw std::out_of_range("only nodes 1,...,N-1 are handled by the scheme");
        }
        if (i == 1){
            return {leftOffsets_, scale*leftWeights_};
        }
        else if (i == N-1) {
            return {rightOffsets_, scale*rightWeights_};
        }
        else {
            return {intOffsets_, scale*intWeights_};
        }

    }
    std::pair<std::vector<int>, Eigen::VectorXd> GenScheme::intOffsetsAndWeights() const{
        return {intOffsets_, intWeights_};
    }
}
