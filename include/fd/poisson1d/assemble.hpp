#pragma once

#include "fd/grid1d.hpp"
#include "core/boundary_condition.hpp"
#include <Eigen/Sparse>
#include <functional>

namespace fd::poisson1d{
    // assembling the A matrix
    Eigen::SparseMatrix<double> assembleA(const Grid1D&);
    //assemble RHS
    Eigen::VectorXd assembleRHS(const Grid1D&, const core::DirichletBC&, std::function<double(double)>);

}