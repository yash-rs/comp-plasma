#include "fd/sparse_solve.hpp"
#include <Eigen/Sparse>

namespace fd {
    Eigen::VectorXd solveSPD(const Eigen::SparseMatrix<double>& A, const Eigen::VectorXd& rhs){
        Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver;
        solver.compute(A);
        return solver.solve(rhs);
    }

    Eigen::VectorXd solveSparseLU(const Eigen::SparseMatrix<double>& A, const Eigen::VectorXd& rhs){
        Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
        solver.compute(A);
        return solver.solve(rhs);
    }
}