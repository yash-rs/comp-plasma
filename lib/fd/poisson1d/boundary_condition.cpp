#include "fd/poisson1d/boundary_condition.hpp"
#include "fd/poisson1d/scheme.hpp"
#include "fd/grid1d.hpp"
#include <vector>

namespace fd::poisson1d{
    DirichletBC::DirichletBC(double alpha, double beta) : alpha_(alpha), beta_(beta) {}

    Eigen::SparseMatrix<double> DirichletBC::assembleA(const Grid1D& grid, const Scheme& scheme) const{
        const int n = grid.N() - 1; //interior points for N-1 x N-1 matrix
        const double h = grid.h(); //mesh size
        Eigen::SparseMatrix<double> A(n, n);
        std::vector<Eigen::Triplet<double>> triplets;
        scheme.addInteriorRows(triplets, grid, 1, n-2);
        //row 0
        triplets.emplace_back(0, 0, scheme.centerCoeff(h));
        triplets.emplace_back(0, 1, scheme.neighborCoeff(h));
        //row n-1
        triplets.emplace_back(n-1, n-1, scheme.centerCoeff(h));
        triplets.emplace_back(n-1, n-2, scheme.neighborCoeff(h));
        A.setFromTriplets(triplets.begin(), triplets.end());
        A.makeCompressed();

        return A;
    }

    Eigen::VectorXd DirichletBC::assembleRHS(const Grid1D& grid, const Scheme& scheme, std::function<double(double)> rho) const{
        const int n = grid.N() - 1;
        const double h = grid.h();
        const Eigen::VectorXd nodes = grid.nodes();
        Eigen::VectorXd rhs(n);
        for (int i = 0; i<n; i++){
            rhs(i) = rho(nodes(i+1));
            if(i == 0) rhs(i) += alpha_*-scheme.neighborCoeff(h);
            if(i == n-1) rhs(i) += beta_*-scheme.neighborCoeff(h);
        }

        return rhs;
    }

    Eigen::VectorXd DirichletBC::reconstructSolution(const Grid1D& grid, const Eigen::VectorXd& reduced_sol) const{
        const int N = grid.N();
        Eigen::VectorXd sol(N+1);
        sol(0) = alpha_;
        sol(N) = beta_;
        for (int i=1;i<N;i++){
            sol(i) = reduced_sol(i-1);
        }

        return sol;
    }

    // Neumann BC
    NeumannBC::NeumannBC(double alpha, double beta) : alpha_(alpha), beta_(beta){}
    Eigen::SparseMatrix<double> NeumannBC::assembleA(const Grid1D& grid, const Scheme& scheme) const {
        const int n = grid.N();
        const double h = grid.h();
        Eigen::SparseMatrix<double> A(n, n);
        std::vector<Eigen::Triplet<double>> triplets;
        scheme.addInteriorRows(triplets, grid, 1, n-2);
        //row 0
        triplets.emplace_back(0, 0, scheme.centerCoeff(h));
        triplets.emplace_back(0, 1, scheme.neighborCoeff(h));
        //row n-1
        triplets.emplace_back(n-1, n-1, scheme.centerCoeff(h));
        triplets.emplace_back(n-1, n-2, 2*scheme.neighborCoeff(h));
        A.setFromTriplets(triplets.begin(), triplets.end());
        A.makeCompressed();

        return A;
    }

    Eigen::VectorXd NeumannBC::assembleRHS(const Grid1D& grid, const Scheme& scheme, std::function<double(double)> rho) const{
        const int n = grid.N();
        const double h = grid.h();
        const Eigen::VectorXd nodes = grid.nodes();
        Eigen::VectorXd rhs(n);
        for (int i = 0; i<n; i++){
            rhs(i) = rho(nodes(i+1));
            if(i == 0) rhs(i) += alpha_*-scheme.neighborCoeff(h);
            if(i == n-1) rhs(i) += 2*beta_/h;
        }

        return rhs;
    }

    Eigen::VectorXd NeumannBC::reconstructSolution(const Grid1D& grid, const Eigen::VectorXd& reduced_sol) const{
        const int N = grid.N();
        Eigen::VectorXd sol(N+1);
        sol(0) = alpha_;
        for (int i=1;i<N+1;i++){
            sol(i) = reduced_sol(i-1);
        }

        return sol;
    }
}