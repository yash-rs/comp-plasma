#pragma once
#include "fd/grid1d.hpp"
#include <Eigen/Sparse>
#include <Eigen/Dense>
#include <functional>

namespace fd::poisson1d{
    class BoundaryCondition{
        public:
            virtual ~BoundaryCondition() = default;
            virtual Eigen::SparseMatrix<double> assembleA(const Grid1D&) const = 0;
            virtual Eigen::VectorXd assembleRHS(const Grid1D&, std::function<double(double)>) const = 0;
            virtual Eigen::VectorXd reconstructSolution(const Grid1D&, const Eigen::VectorXd&) const = 0;

    };
    //DirichletBC
    class DirichletBC : public BoundaryCondition{
        double alpha_;
        double beta_;
        public:
            DirichletBC(double alpha, double beta);
            Eigen::SparseMatrix<double> assembleA(const Grid1D&) const override;
            Eigen::VectorXd assembleRHS(const Grid1D&, std::function<double(double)>) const override;
            Eigen::VectorXd reconstructSolution(const Grid1D&, const Eigen::VectorXd&) const override;
    };
}