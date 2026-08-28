#pragma once

#include "fd/grid1d.hpp"
#include "core/boundary_condition.hpp"
#include <Eigen/Sparse>
#include <Eigen/Dense>

namespace fd {
    Eigen::VectorXd solveDirichletFD(const Eigen::SparseMatrix<double>&, const Eigen::VectorXd&, const Grid1D&, const core::DirichletBC&);
}