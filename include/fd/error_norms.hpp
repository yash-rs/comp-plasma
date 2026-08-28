#pragma once

#include "fd/grid1d.hpp"
#include <Eigen/Dense>

namespace fd {
    double l1Norm(const Eigen::VectorXd&, const Eigen::VectorXd&, const Grid1D&);
    double l2Norm(const Eigen::VectorXd&, const Eigen::VectorXd&, const Grid1D&);
    double linfNorm(const Eigen::VectorXd&, const Eigen::VectorXd&, const Grid1D&);

}