#include "fd/poisson1d/assemble.hpp"
#include "fd/grid1d.hpp"
#include "core/boundary_condition.hpp"
#include <vector>
#include <Eigen/Sparse>
#include <Eigen/Dense>

namespace fd::poisson1d{
    Eigen::SparseMatrix<double> assembleA(const Grid1D& grid){
        const int n = grid.N() - 1; //interior points for N-1 x N-1 matrix
        const double h = grid.h(); //mesh size
        Eigen::SparseMatrix<double> A(n, n);
        std::vector<Eigen::Triplet<double>> triplets;
        for (int i=0; i<n; i++){
            triplets.emplace_back(i, i, 2.0/(h*h));
            if (i > 0){
                triplets.emplace_back(i, i-1, -1.0/(h*h));
            }
            if (i < n-1){
                triplets.emplace_back(i, i+1, -1.0/(h*h));
            }
        }
        A.setFromTriplets(triplets.begin(), triplets.end());
        A.makeCompressed();

        return A;
    }

    Eigen::VectorXd assembleRHS(const Grid1D& grid, const core::DirichletBC& bc, std::function<double(double)> rho){
        const int n = grid.N() - 1;
        const double h = grid.h();
        const double alpha = bc.alpha; // left dirichlet BC
        const double beta = bc.beta; // right dirichlet BC
        const Eigen::VectorXd nodes = grid.nodes();
        Eigen::VectorXd rhs(n);
        for (int i = 0; i<n; i++){
            rhs(i) = rho(nodes(i+1));
            if(i == 0) rhs(i) += alpha/(h*h);
            if(i == n-1) rhs(i) += beta/(h*h);
        }

        return rhs;
    }
}