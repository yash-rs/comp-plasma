#include "fd/grid1d.hpp"
#include "core/boundary_condition.hpp"
#include <Eigen/Sparse>
#include <Eigen/Dense>


namespace fd {
    Eigen::VectorXd solveDirichletFD(const Eigen::SparseMatrix<double>& A, const Eigen::VectorXd& rhs, const Grid1D& grid, const core::DirichletBC& bc){
        const int N = grid.N(); //full grid needed
        const double alpha = bc.alpha;
        const double beta = bc.beta;
        Eigen::VectorXd sol(N+1);
        //solve the sparse system
        Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver;
        solver.compute(A);
        Eigen::VectorXd sol_ = solver.solve(rhs); //interior points
        sol(0) = alpha;
        sol(N) = beta;
        for (int i=1; i<N; i++){
            sol(i) = sol_(i-1);
        }

        return sol;

    }

}