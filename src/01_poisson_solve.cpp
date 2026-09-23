#include "fd/grid1d.hpp"
#include "fd/sparse_solve.hpp"
#include "fd/poisson1d/boundary_condition.hpp"
#include "core/manufactured_solutions.hpp"
#include "fd/poisson1d/scheme.hpp"
#include <cmath>
#include <functional>
#include <iostream>
#include <Eigen/Sparse>
#include <Eigen/Dense>
#include <fstream>

using namespace core;
using namespace fd;

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
    poisson1d::DirichletBC bc(alpha, beta);
    const int N = 500;
    const Grid1D grid(a, b, N);
    //scheme
    const poisson1d::SecondOrderScheme scheme;
    Eigen::SparseMatrix<double> A = bc.assembleA(grid, scheme);
    Eigen::VectorXd rhs = bc.assembleRHS(grid, scheme, rho_fn);
    Eigen::VectorXd reduced_sol = solveSPD(A, rhs);
    Eigen::VectorXd sol = bc.reconstructSolution(grid, reduced_sol);
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