#pragma once

#include <Eigen/Dense>

namespace fd{
    Eigen::VectorXd stencilWeights(const Eigen::VectorXd& offsets, int order);
}