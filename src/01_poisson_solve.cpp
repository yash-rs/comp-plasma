#include "fd/grid1d.hpp"
#include "fd/poisson1d/assemble.hpp"
#include "core/boundary_condition.hpp"
#include "core/manufactured_solutions.hpp"
#include "fd/dirichlet_solve.hpp"
#include <cmath>
#include <functional>
#include <iostream>
#include <Eigen/Sparse>
#include <Eigen/Dense>
#include <fstream>

using namespace core;
using namespace fd;
using namespace fd::poisson1d;

int main(){
    /* Setup 1d Poisson for a test case and
    compare against manufactured solutions*/
    std::function<double(double)> rho_fn = [](double x){return 2*std::sin(x)+x*std::cos(x);};
    std::function<double(double)> phi_fn = [](double x){return x*std::cos(x);};
    ManufacturedSolution ms{rho_fn, phi_fn};
    // steup domain and mesh
    const double a = 0.;
    const double b = 2*std::acos(-1.0); //2pi
    // bc values
    const double alpha = phi_fn(a);
    const double beta = phi_fn(b);
    DirichletBC bc{alpha, beta};
    const int N = 500;
    const Grid1D grid(a, b, N);
    Eigen::SparseMatrix<double> A = assembleA(grid);
    Eigen::VectorXd rhs = assembleRHS(grid, bc, ms.rho);
    Eigen::VectorXd sol = solveDirichletFD(A, rhs, grid, bc);
    // exact solution
    Eigen::VectorXd sol_exact(grid.N()+1);
    Eigen::VectorXd nodes = grid.nodes();
    for (int i = 0; i<N+1; i++){
        sol_exact(i) = ms.phi_exact(nodes(i));
    }
    // write solution
    std::ofstream file("data/poisson_solution.csv");
    file << "x,phi_numeric,phi_exact\n";
    for (int i=0; i<N+1; i++){
        file << nodes(i)<<","<<sol(i)<<","<< sol_exact(i)<<"\n";
    }

    

}