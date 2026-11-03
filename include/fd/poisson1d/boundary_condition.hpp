#pragma once
#include "fd/grid1d.hpp"
#include "fd/poisson1d/scheme.hpp"
#include <Eigen/Sparse>
#include <Eigen/Dense>
#include <functional>
#include <vector>

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

    // Gen BC, to be refactored later
    class MixedDirichletNeumannBC{
        double alpha_;
        double beta_;
        std::vector<int> rbo_;
        public:
            MixedDirichletNeumannBC(double alpha, double beta, const std::vector<int>& rightBoundaryOffsets);
            Eigen::SparseMatrix<double> assembleA(const Grid1D&, const GenScheme&) const;
            Eigen::VectorXd assembleRHS(const Grid1D&, const GenScheme&, std::function<double(double)>) const;
            Eigen::VectorXd reconstructSolution(const Grid1D&, const Eigen::VectorXd&) const;

    };

    class PeriodicBC{
        Eigen::VectorXcd computeEigs(const Grid1D&, const GenScheme&) const; // private function to calculate eigenvalues of Ah as they might differ for diffferent schemes
        public:
            PeriodicBC();
            Eigen::VectorXd assembleRHS(const Grid1D&, std::function<double(double)>) const;
            Eigen::VectorXd solveWithDFT(const Grid1D&, const GenScheme&, const Eigen::VectorXd&) const;
            Eigen::VectorXd reconstructSolution(const Grid1D&, const Eigen::VectorXd&) const;
    };
}