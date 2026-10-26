#pragma once

#include <Eigen/Dense>

namespace core{
    double convergenceRate(const Eigen::VectorXd& h, const Eigen::VectorXd& error);
}