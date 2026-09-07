#include "fd/poisson1d/boundary_condition.hpp"
#include "fd/grid1d.hpp"
#include <vector>

namespace fd::poisson1d{
    DirichletBC::DirichletBC(double alpha, double beta) : alpha_(alpha), beta_(beta) {}

    Eigen::SparseMatrix<double> DirichletBC::assembleA(const Grid1D& grid) const{
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

    Eigen::VectorXd DirichletBC::assembleRHS(const Grid1D& grid, std::function<double(double)> rho) const{
        const int n = grid.N() - 1;
        const double h = grid.h();
        const Eigen::VectorXd nodes = grid.nodes();
        Eigen::VectorXd rhs(n);
        for (int i = 0; i<n; i++){
            rhs(i) = rho(nodes(i+1));
            if(i == 0) rhs(i) += alpha_/(h*h);
            if(i == n-1) rhs(i) += beta_/(h*h);
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
}