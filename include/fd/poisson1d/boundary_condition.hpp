#pragma once
#include "fd/grid1d.hpp"
#include "fd/poisson1d/scheme.hpp"
#include <Eigen/Sparse>
#include <Eigen/Dense>
#include <functional>

namespace fd::poisson1d{
    class BoundaryCondition{
        public:
            virtual ~BoundaryCondition() = default;
            virtual Eigen::SparseMatrix<double> assembleA(const Grid1D&, const Scheme&) const = 0;
            virtual Eigen::VectorXd assembleRHS(const Grid1D&, const Scheme&, std::function<double(double)>) const = 0;
            virtual Eigen::VectorXd reconstructSolution(const Grid1D&, const Eigen::VectorXd&) const = 0;

    };
    //DirichletBC
    class DirichletBC : public BoundaryCondition{
        double alpha_;
        double beta_;
        public:
            DirichletBC(double alpha, double beta);
            Eigen::SparseMatrix<double> assembleA(const Grid1D&, const Scheme&) const override;
            Eigen::VectorXd assembleRHS(const Grid1D&, const Scheme&, std::function<double(double)>) const override;
            Eigen::VectorXd reconstructSolution(const Grid1D&, const Eigen::VectorXd&) const override;
    };

    class NeumannBC : public BoundaryCondition{
        double alpha_; //left boundary dirichlet
        double beta_; //right boundary neumann
        public:
            NeumannBC(double alpha, double beta);
            // scheme related functions
            Eigen::SparseMatrix<double> assembleA(const Grid1D&, const Scheme&) const override;
            Eigen::VectorXd assembleRHS(const Grid1D&, const Scheme&, std::function<double(double)>) const override;
            Eigen::VectorXd reconstructSolution(const Grid1D&, const Eigen::VectorXd&) const override;
    };
}