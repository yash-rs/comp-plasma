#include "fd/grid1d.hpp"
#include "fd/error_norms.hpp"
#include <cmath>

namespace fd {
    double l1Norm(const Eigen::VectorXd& numeric, const Eigen::VectorXd& exact, const Grid1D& grid){
        double h = grid.h();
        double norm = 0;
        for (int i=0; i<grid.N()+1; i++){
            norm += std::abs(numeric(i) - exact(i));
        }
        norm *= h;
        return norm;
    }

    double l2Norm(const Eigen::VectorXd& numeric, const Eigen::VectorXd& exact, const Grid1D& grid){
        double h = grid.h();
        double norm = 0;
        for (int i=0; i<grid.N()+1; i++){
            norm += (numeric(i) - exact(i))*(numeric(i) - exact(i));
        }
        norm *= h;
        norm = std::sqrt(norm);
        return norm;
    }

    double linfNorm(const Eigen::VectorXd& numeric, const Eigen::VectorXd& exact, const Grid1D& grid){
        double norm = 0;
        double norm_ = 0;
        for (int i=0; i<grid.N()+1; i++){
            norm_ = std::abs(numeric(i) - exact(i));
            if (norm_ > norm) norm = norm_;

        }
        return norm;
    }

}
