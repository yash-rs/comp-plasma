#pragma once

#include <Eigen/Sparse>
#include <Eigen/Dense>

namespace fd{
    Eigen::VectorXd solveSPD(const Eigen::SparseMatrix<double>&, const Eigen::VectorXd&);
}