#pragma once

#include <Eigen/Sparse>
#include <Eigen/Dense>

namespace fd{
    Eigen::VectorXd solveSPD(const Eigen::SparseMatrix<double>&, const Eigen::VectorXd&);
    Eigen::VectorXd solveSparseLU(const Eigen::SparseMatrix<double>&, const Eigen::VectorXd&);
}